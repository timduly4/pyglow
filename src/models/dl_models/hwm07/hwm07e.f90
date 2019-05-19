!
!  Horizontal Wind Model 07 (HWM07)
!  Version HWM071308E, 13 July 2008
!  See readme.txt file for detailed release notes.
!
!  AUTHOR
!    Douglas P. Drob
!    Space Science Division
!    Naval Research Laboratory
!    4555 Overlook Ave.
!    Washington, DC 20375
!
!  Point of Contact
!    msishwmhelp@nrl.navy.mil
!
!  DATE
!    19 August 2008
!
!  REFERENCE
!    Drob, D. P, J. T. Emmert, G. Crowley, J. M. Picone, G. G. Shepherd, 
!      W. Skinner, Paul Hayes, R. J. Niciejewski, M. Larsen, C.Y. She, 
!      J. W. Meriwether, G. Hernandez, M. J. Jarvis, D. P. Sipler, C. A. Tepley,
!      M. S. O�Brien, J. R. Bowman, Q. Wu, Y. Murayama, S. Kawamura, I.M. Reid,
!      and R.A. Vincent (2008), An Empirical Model of the Earth�s Horizontal 
!      Wind Fields: HWM07, J. Geophy. Res., doi:10.1029/2008JA013668.
!
!==================================================================================
! Input arguments:
!        iyd - year and day as yyddd
!        sec - ut(sec)
!        alt - altitude(km)
!        glat - geodetic latitude(deg)
!        glon - geodetic longitude(deg)
!        stl - not used
!        f107a - not used
!        f107 - not used
!        ap - two element array with
!             ap(1) = not used
!             ap(2) = current 3hr ap index
!
! Output argument:
!        w(1) = meridional wind (m/sec + northward)
!        w(2) = zonal wind (m/sec + eastward)
!
!================================================================================


subroutine hwm07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)

    implicit none
    integer(4),intent(in)   :: iyd
    real(4),intent(in)      :: sec,alt,glat,glon,stl,f107a,f107
    real(4),intent(in)      :: ap(2)
    real(4),intent(out)     :: w(2)

    real(4)                 :: qw(2),dw(2)

    call hwmqt(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,qw)
    
    if (ap(2) .ge. 0.0) then
      call dwm07b_hwm_interface(iyd,sec,alt,glat,glon,ap,dw)
      w = qw + dw
    else
      w = qw
    endif
    
    return
    
end subroutine HWM07 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                           Global Data Module NEWmodel
!
! Description: This is a common data module for model definition.  These 
!  parameters set by the first calling the subroutine loadmodel().
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module NEWmodel

    implicit none

    integer                 :: nbf              ! Count of basis terms per model level   
    integer                 :: maxs             ! s seasonal
    integer                 :: maxm             ! m stationary
    integer                 :: maxl             ! l migrating
    integer                 :: maxn             ! n latitude
    
    integer                 :: p                ! B-splines order, p=4 cubic, p=3 quadratic 
    integer                 :: nlev             ! e.g. Number of B-spline nodes
    integer                 :: nnode            ! nlev + p

    real(8)                 :: alttns           ! Transition 1
    real(8)                 :: altsym           ! Transition 2
    real(8)                 :: altiso           ! Constant Limit

    integer,allocatable     :: nb(:)            ! total number of basis functions @ level
    integer,allocatable     :: order(:,:)       ! spectral content @ level
    real(8),allocatable     :: vnode(:)         ! Vertical Altitude Nodes
    real(8),allocatable     :: mparm(:,:)       ! Model Parameters
    
    ! Global store for quasi-static model space parameters
    ! These will change internally depending on the input parameters
    
    real(8),allocatable     :: gfs(:,:),gfm(:,:),gfl(:,:)
    real(8),allocatable     :: gvbar(:,:),gwbar(:,:)
    real(8),allocatable     :: gbz(:,:),gbm(:,:)
    
    real(8),allocatable     :: gzwght(:)
    integer                 :: glev
     
    ! Miscellaneous flags and indicies
    
    integer                 :: maxo
    integer                 :: cseason = 0
    integer                 :: cwave = 0
    integer                 :: ctide = 0
    
    logical                 :: content(5) = .true.          ! Season/Waves/Tides
    logical                 :: component(0:1) = .true.      ! Compute zonal/meridional
    
    ! Initialization flags and information
    
    logical                 :: modelinit = .true.
    logical                 :: reset = .true.
    character(128)          :: defaultdata = 'hwm071308e.dat'
    
end module NEWmodel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                       HWM-07 Legacy Wrapper
!
!  Description: Emulate HWM-93 subroutine calling convension
!
!  Programming Notes:
!
!
!  Required Subroutines:
!
!    loadmodel()
!    HWMupdate()
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine HWMQT(IYD,SEC,ALT,GLAT,GLON,STL,F107A,F107,AP,W)
    
    use NEWmodel
    implicit none

    INTEGER,intent(in)      :: IYD
    REAL(4),intent(in)      :: SEC,ALT,GLAT,GLON,STL,F107A,F107
    REAL(4),intent(in)      :: AP(2)
    REAL(4),intent(out)     :: W(2)

    real(8)                 :: last(5)
    real(8)                 :: input(5)
    real(8)                 :: u,v
    
    input(1) = dble(mod(IYD,1000))
    input(2) = dble(sec)
    input(3) = dble(glon)
    input(4) = dble(glat)
    input(5) = dble(alt)

    if (modelinit) then
        call loadmodel(defaultdata)
        call HWMupdate(input,last,gfs,gfl,gfm,gvbar,gwbar,gbz,gbm,gzwght,glev,u,v)
    endif
    
    if (reset) then
        last = 1.0d-32
        reset = .false.
    endif
   
    call HWMupdate(input,last,gfs,gfl,gfm,gvbar,gwbar,gbz,gbm,gzwght,glev,u,v)
              
    w(1) = -v
    w(2) = u

    return

end subroutine HWMQT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   This subroutine is used to calculate the Vector Spherical
!   Harmonic basis functions for a given observation.
!
!  Programming Notes:
! 
!   This subroutine is only OPENMP/THREAD SAFE when no calls to
!   loadmodel() are made.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine HWMupdate(input,last,fs,fl,fm,vbar,wbar,ebz,ebm,zwght,lev,u,v)

    use NewModel
    implicit none
    
    real(8),intent(in)              :: input(5) ! jday,utsec,glon,glat,alt
    real(8),intent(inout)           :: last(5)
    real(8),intent(inout)           :: fs(0:maxs,2)
    real(8),intent(inout)           :: fm(0:maxm,2)
    real(8),intent(inout)           :: fl(0:maxl,2)
    
    real(8),intent(inout)           :: vbar(0:maxn,0:maxo)
    real(8),intent(inout)           :: wbar(0:maxn,0:maxo)
    real(8),intent(inout),target    :: ebz(nbf,0:p)
    real(8),intent(inout),target    :: ebm(nbf,0:p)

    real(8),intent(inout)           :: zwght(0:p)
    integer,intent(inout)           :: lev
   
    real(8),intent(out)             :: u,v
 
    ! Local variables

    real(8),pointer           :: bz(:)
    real(8),pointer           :: bm(:)
    
    real(8)                   :: cs,ss,cm,sm,cl,sl
    real(8)                   :: cmcs,smcs,cmss,smss
    real(8)                   :: clcs,slcs,clss,slss	
    real(8)                   :: AA,BB,CC,DD
    real(8)                   :: vb,wb
    
    integer                   :: b,c,d,m,n,s,l
    
    integer                   :: amaxs,amaxn
    integer                   :: pmaxm,pmaxs,pmaxn
    integer                   :: tmaxl,tmaxs,tmaxn

    logical                   :: refresh(5)
       
    real(8),parameter         :: twoPi = 2.0d0*3.1415926535897932384626433832795d0
    real(8),parameter         :: deg2rad = twoPi/360.0d0
    
! ====================================================================
! Update VSH model terms based on any change in the input parameters
! ====================================================================
    
    refresh = .false.

    ! Seasonal variations

    if (input(1) .ne. last(1)) then
        AA = input(1)*twoPi/365.25d0
        do s = 0,MAXS
            BB = dble(s)*AA
            fs(s,1) = dcos(BB)
            fs(s,2) = dsin(BB)
        enddo
        refresh(1:5) = .true.
        last(1) = input(1)
    endif

    ! Hourly time changes, tidal variations

    if (input(2) .ne. last(2) .or. input(3) .ne. last(3)) then

        AA = mod(input(2)/3600.d0 + input(3)/15.d0 + 48.d0,24.d0)
        BB = AA*twoPi/24.d0
        do l = 0,MAXL
            CC = dble(l)*BB
            fl(l,1) = dcos(CC)
            fl(l,2) = dsin(CC)
        enddo

        refresh(3) = .true.
        last(2) = input(2)
    endif

    ! Longitudinal variations, planetary waves
    
    if (input(3) .ne. last(3)) then
        AA = input(3)*deg2rad
        do m = 0,MAXM
            BB = dble(m)*AA
            fm(m,1) = dcos(BB)
            fm(m,2) = dsin(BB)
        enddo
        refresh(2) = .true.
        last(3) = input(3)
    endif

    ! Latitude

    if (input(4) .ne. last(4)) then
        AA = (90.0d0 - input(4))*deg2rad        ! theta = colatitude in radians
        call vshbasis(maxn,maxo,AA,vbar,wbar)
        refresh(1) = .true.
        refresh(2) = .true.
        refresh(3) = .true.
        refresh(4) = .true.
        last(4) = input(4)
    endif

    ! Altitude
    
    if (input(5) .ne. last(5)) then
        call vertwght(input(5),zwght,lev)
        last(5) = input(5)
    endif

    ! ====================================================================
    ! Linearize the VSH functions
    ! ====================================================================
       
    u = 0.0d0
    v = 0.0d0
       
    refresh = .true.
       
    ebz = 0.0
    ebm = 0.0
    
    do b = 0,p

        if (zwght(b) .eq. 0.d0) cycle
        
        d = b + lev
        
        bz => ebz(:,b)
        bm => ebm(:,b)

        amaxs = order(1,d)
        amaxn = order(2,d)
        pmaxm = order(3,d)
        pmaxs = order(4,d)
        pmaxn = order(5,d)
        tmaxl = order(6,d)
        tmaxs = order(7,d)
        tmaxn = order(8,d)

        c = 1
 
        ! ------------- Seasonal - Zonal average (m = 0) ----------------

        if (refresh(1) .and. content(1)) then

            c = 1                           ! global constants
            do n = 1,AMAXN
                bz(c)   = -0.5d0*vbar(n,0)  ! Cr
                bz(c+1) =  0.0d0            ! Br
                bm(c)   =  0.0d0            ! Cr
                bm(c+1) =  0.5d0*vbar(n,0)  ! Br
                c = c + 2
            enddo
                
            do s = 1,AMAXS                   ! seasonal variations
                cs = fs(s,1)
                ss = fs(s,2)
                do n = s,AMAXN
                    vb = vbar(n,s)
                    wb = wbar(n,s)
                    AA =  vb*cs
                    BB =  vb*ss
                    CC = -wb*ss
                    DD = -wb*cs
                      bz(c) = -AA   ! Cr
                    bz(c+1) =  BB   ! Ci
                    bz(c+2) =  CC   ! Br
                    bz(c+3) =  DD   ! Bi
                      bm(c) =  CC   ! Cr
                    bm(c+1) =  DD   ! Ci
                    bm(c+2) =  AA   ! Br
                    bm(c+3) = -BB   ! Bi
                    c = c + 4
                enddo
            enddo
            cseason = c
        else
            c = cseason
        endif
            
        ! ---------------- Stationary planetary waves --------------------

        if (refresh(2) .and. content(2)) then

            do m = 1,pmaxm

               cm = fm(m,1)
               sm = fm(m,2)

               do n = m,pmaxn           ! s = 0
               
                    vb = vbar(n,m)
                    wb = wbar(n,m)
               
                    bz(c) =   -vb*cm    ! Cr * (cm) * -vb
                    bz(c+1) =  vb*sm    ! Ci * (sm) *  vb
                    bz(c+2) = -wb*sm    ! Br * (sm) * -wb
                    bz(c+3) = -wb*cm    ! Bi * (cm) * -wb
                    
                    bm(c) =   -wb*sm    ! Cr * (sm) * -wb
                    bm(c+1) = -wb*cm    ! Ci * (sm) * -wb
                    bm(c+2) =  vb*cm    ! Br * (cm) *  vb
                    bm(c+3) = -vb*sm    ! Bi * (sm) * -vb
                    
                    c = c + 4
               
               enddo
               
               do s = 1,pmaxs
               
                  cs = fs(s,1)
                  ss = fs(s,2)
               
                  do n = m,pmaxn
                     vb = vbar(n,m)
                     wb = wbar(n,m)
                     
                     bz(c) =   -vb*cm*cs    ! Crc * (cmcs) * -vb
                     bz(c+1) =  vb*sm*cs    ! Cic * (smcs) *  vb
                     bz(c+2) = -wb*sm*cs    ! Brc * (smcs) * -wb
                     bz(c+3) = -wb*cm*cs    ! Bic * (cmcs) * -wb
                     bz(c+4) = -vb*cm*ss    ! Crs * (cmss) * -vb
                     bz(c+5) =  vb*sm*ss    ! Cis * (smss) *  vb
                     bz(c+6) = -wb*sm*ss    ! Brs * (smss) * -wb
                     bz(c+7) = -wb*cm*ss    ! Bis * (cmss) * -wb
                     
                     bm(c) =   -wb*sm*cs    ! Crc * (smcs) * -wb
                     bm(c+1) = -wb*cm*cs    ! Cic * (smcs) * -wb
                     bm(c+2) =  vb*cm*cs    ! Brc * (cmcs) *  vb
                     bm(c+3) = -vb*sm*cs    ! Bic * (smcs) * -vb
                     bm(c+4) = -wb*sm*ss    ! Crs * (smss) * -wb
                     bm(c+5) = -wb*cm*ss    ! Cis * (smss) * -wb
                     bm(c+6) =  vb*cm*ss    ! Brs * (cmss) *  vb
                     bm(c+7) = -vb*sm*ss    ! Bis * (smss) * -vb
                     
                     c = c + 8
                  
                  enddo
            
               enddo
               cwave = c      
            enddo        
        else
            c = cwave  
        endif

        ! ---------------- Migrating Solar Tides ---------------------

        if (refresh(3) .and. content(3)) then
             do l = 1,tmaxl
            
               cl = fl(l,1)
               sl = fl(l,2)

               s = 0
               do n = l,tmaxn
               
                    vb = vbar(n,l)
                    wb = wbar(n,l)                                          
               
                    bz(c) =   -vb*cl    ! Cr * (cl) * -vb
                    bz(c+1) =  vb*sl    ! Ci * (sl) *  vb
                    bz(c+2) = -wb*sl    ! Br * (sl) * -wb
                    bz(c+3) = -wb*cl    ! Bi * (cl) * -wb
                    
                    bm(c) =   -wb*sl    ! Cr * (sl) * -wb
                    bm(c+1) = -wb*cl    ! Ci * (sl) * -wb
                    bm(c+2) =  vb*cl    ! Br * (cl) *  vb
                    bm(c+3) = -vb*sl    ! Bi * (sl) * -vb
                    
                    c = c + 4
               
               enddo
               
               do s = 1,tmaxs
                  
                  cs = fs(s,1)
                  ss = fs(s,2)
                  
                  do n = l,tmaxn
                  
                     vb = vbar(n,l)
                     wb = wbar(n,l)
                                   
                     bz(c) =   -vb*cl*cs	! Crc * (clcs) * -vb
                     bz(c+1) =  vb*sl*cs    ! Cic * (slcs) *  vb
                     bz(c+2) = -wb*sl*cs    ! Brc * (slcs) * -wb
                     bz(c+3) = -wb*cl*cs    ! Bic * (clcs) * -wb
                     
                     bz(c+4) = -vb*cl*ss    ! Crs * (clss) * -vb
                     bz(c+5) =  vb*sl*ss    ! Cis * (slss) *  vb
                     bz(c+6) = -wb*sl*ss    ! Brs * (slss) * -wb
                     bz(c+7) = -wb*cl*ss    ! Bis * (clss) * -wb
                                         
                     bm(c) =   -wb*sl*cs    ! Crc * (slcs) * -wb
                     bm(c+1) = -wb*cl*cs    ! Cic * (slcs) * -wb
                     bm(c+2) =  vb*cl*cs    ! Brc * (clcs) *  vb
                     bm(c+3) = -vb*sl*cs    ! Bic * (slcs) * -vb

                     bm(c+4) = -wb*sl*ss    ! Crs * (slss) * -wb
                     bm(c+5) = -wb*cl*ss    ! Cis * (slss) * -wb
                     bm(c+6) =  vb*cl*ss    ! Brs * (clss) *  vb
                     bm(c+7) = -vb*sl*ss    ! Bis * (slss) * -vb
                     
                     c = c + 8
                  enddo
            
               enddo
               ctide = c
            enddo
        else
            c = ctide  
        endif
            
        ! ---------------- Non-Migrating Solar Tides ------------------
        
        ! TBD
            
        c = c - 1
        
        ! ====================================================================
        ! Calculate the wind components 
        ! ====================================================================

        if (component(0)) u = u + zwght(b)*dot_product(bz(1:c),mparm(1:c,d))
        if (component(1)) v = v + zwght(b)*dot_product(bm(1:c),mparm(1:c,d))

    enddo

    return 

end subroutine HWMupdate

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine initalizes the NEWmodel
!
!
! Required Subroutines:
!
!   vshengineinit()
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine loadmodel(datafile)

    use NEWmodel
    implicit none

    character(128),intent(in)   :: datafile
    
    integer                     :: i,j
    integer                     :: ncomp

    if (allocated(vnode)) then
        deallocate(order,nb,vnode,mparm)
        deallocate(gfs,gfm,gfl,gvbar,gwbar,gzwght,gbz,gbm)
    endif

    open(unit=23,file=trim(datafile),form='unformatted')
    read(23) nbf,maxs,maxm,maxl,maxn,ncomp
    read(23) nlev,p
    nnode = nlev + p
    allocate(nb(0:nnode))
    allocate(order(ncomp,0:nnode))
    allocate(vnode(0:nnode))
    read(23) vnode
    vnode(3) = 0.0
    allocate(mparm(nbf,0:nlev))
    mparm = 0.0d0
    do i = 0,nlev-p+1-2
        read(23) order(1:ncomp,i)
        read(23) nb(i)
        read(23) mparm(1:nbf,i)
    enddo
    close(23)

    ! Set transition levels
    
    alttns = vnode(nlev-2)
    altsym = vnode(nlev-1)
    altiso = vnode(nlev)

    ! Initialize the vector spherical harmonics engine

    call vshengineinit()

    ! Allocate the global store of quasi-static parameters
    
    maxo = max(maxs,maxm,maxl)

    allocate(gfs(0:maxs,2),gfm(0:maxm,2),gfl(0:maxl,2))
    allocate(gvbar(0:maxn,0:maxo),gwbar(0:maxn,0:maxo))
    allocate(gbz(nbf,0:p),gbm(nbf,0:p))
    allocate(gzwght(0:p))
    
    gvbar = 0.0d0
    gwbar = 0.0d0
    gbz = 0.0d0
    gbm = 0.0d0
  
    ! Signal that the model has been initalized
    
    modelinit = .false.
    
    ! Signal a reset of the input variable comparison flags
    
    reset = .true.
    
    return

end subroutine loadmodel
 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine updates the vertical weights
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine vertwght(alt,wght,iz)

    use NEWmodel
    implicit none
    
    real(8),intent(in)  :: alt
    real(8),intent(out) :: wght(4)
    integer,intent(out) :: iz
    
    real(8)             :: we(0:4)

    real(8)             :: e1(0:4) = &
        (/1.d0, 0.428251121076233d0,0.192825112107623d0,0.484304932735426d0,0.0d0/)
    real(8)             :: e2(0:4) = &
        (/0.d0, 0.571748878923767d0,0.807174887892377d0,-0.484304932735426d0,1.0d0/)

    real(8),parameter   :: H = 60.0d0
        
    iz = findspan(nnode-p-1,p,alt,vnode) - p

    iz = min(iz,26)
    
    wght(1) = bspline(p,nnode,vnode,iz,alt)
    wght(2) = bspline(p,nnode,vnode,iz+1,alt)
    if (iz .le. 25) then
        wght(3) = bspline(p,nnode,vnode,iz+2,alt)
        wght(4) = bspline(p,nnode,vnode,iz+3,alt)
        return
    endif
        if (alt .gt. 250.0d0) then
            we(0) = 0.0d0
            we(1) = 0.0d0
            we(2) = 0.0d0
            we(3) = exp(-(alt - 250.0d0)/H) 
            we(4) = 1.0d0
        else
            we(0) = bspline(p,nnode,vnode,iz+2,alt)
            we(1) = bspline(p,nnode,vnode,iz+3,alt)
            we(2) = bspline(p,nnode,vnode,iz+4,alt)
            we(3) = 0.0d0
            we(4) = 0.0d0
        endif
        wght(3) = dot_product(we,e1)
        wght(4) = dot_product(we,e2)
        
    return
    
contains

    function bspline(p,m,V,i,u)

        implicit none
        
        real(8)     :: bspline
        integer     :: p,m
        real(8)     :: V(0:m)
        integer     :: i
        real(8)     :: u
        
        real(8)     :: N(0:p+1)
        real(8)     :: Vleft,Vright
        real(8)     :: saved,temp
        integer     :: j,k
                
        if ((i .eq. 0) .and. (u .eq. V(0))) then
            bspline = 1.d0
            return
        endif
        
        if ((i .eq. (m-p-1)) .and. (u .eq. V(m))) then
            bspline = 1.d0
            return
        endif

        if (u .lt. V(i) .or. u .ge. V(i+p+1)) then
            bspline = 0.d0
            return
        endif
        
        N = 0.0d0
        do j = 0,p
            if (u .ge. V(i+j) .and. u .lt. V(i+j+1)) then
                N(j) = 1.0d0
            else
                N(j) = 0.0d0
            endif
        enddo
        
        do k = 1,p
            if (N(0) .eq. 0.d0) then
                saved = 0.d0
            else
                saved = ((u - V(i))*N(0))/(V(i+k) - V(i))
            endif
            do j = 0,p-k
                Vleft = V(i+j+1)
                Vright = V(i+j+k+1)
                if (N(j+1) .eq. 0.d0) then
                    N(j) = saved
                    saved = 0.d0
                else
                    temp = N(j+1)/(Vright - Vleft)
                    N(j) = saved + (Vright - u)*temp
                    saved = (u - Vleft)*temp
                endif
            enddo
        enddo
        
        bspline = N(0)

        return

    end function bspline

    ! =====================================================
    ! Function to locate the knot span
    ! =====================================================

    integer function findspan(n,p,u,V)

        implicit none
        
        integer,intent(in)      :: n,p
        real(8),intent(in)      :: u
        real(8),intent(in)      :: V(0:n+1)
        integer                 :: low,mid,high
        
        if (u .ge. V(n+1)) then
            findspan = n
            return
        endif
        
        low = p
        high = n+1
        mid = (low + high)/2

        do while (u .lt. V(mid) .or. u .ge. V(mid + 1))
            if (u .lt. V(mid)) then
                high = mid
            else
                low = mid
            endif
            mid = (low + high)/2
        end do
    
        findspan = mid
        return

    end function findspan
    
end subroutine vertwght

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Description:  Computational Engine for calculating Normalized Associated Vector
!                 Spherical Harmonic basis functions Vbar and Wbar.
!
!   Author: Douglas P. Drob
!           Space Science Division
!           Naval Research Laboratory
!           4555 Overlook Ave
!           Washington, DC
!
!   Date: January 2007.
!
!   Notes:
!
!   The routines for the calculation Pbar have been adapted from ALFPACK developed 
!   by Paul Swarztrauber at The National Center for Atmospheric Research, Boulder, 
!   Colorado (80307) U.S.A. which is sponsored by The National Science Foundation.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module vshengine
    implicit none
    integer,parameter          :: mmax = 11
    integer,parameter          :: nmax = 11
    real(8)                    :: p0(0:nmax,0:mmax),p1(0:nmax,0:mmax)
    real(8)                    :: sf(0:nmax,0:mmax)
    real(8)                    :: a(0:nmax,0:mmax),b(0:nmax,0:mmax)
    real(8)                    :: c(0:nmax,0:mmax),d(0:nmax,0:mmax)    
end module vshengine

! ==========================================================================
! Description: Subroutine to calculate the VSH coeffiecents @ theta
!
! Notes: Before calling vshbasis() the static coeffiecents must calculated
!        by calling vshengineinit().
!
! ==========================================================================

subroutine vshbasis(maxn,maxo,theta,vbar,wbar)

    use vshengine
    implicit none

    integer,intent(in)  :: maxn,maxo
    real(8),intent(in)  :: theta
    real(8),intent(out) :: vbar(0:maxn,0:maxo)
    real(8),intent(out) :: wbar(0:maxn,0:maxo)

    real(8)             :: pbar(0:nmax,0:maxo+1)
    real(8)             :: pb(0:nmax),td(0:nmax)
    real(8)             :: p0i,p1i,cost,r
    integer             :: i,m,mm1,mp1,n,nm1,nm2,nmm

    ! Normalized Associated Legendre Polynomials

    pbar = 0.0d0

    nm1 = nmax - 1
    nm2 = nmax - 2
    do m = 0,maxo+1
        nmm = nmax - m
        p0i = lfpt(nm1,m,theta,p0(:,m))
        p1i = lfpt(nm2,m,theta,p1(:,m))
        pbar(nm1,m) = p0i
        if (nmm .le. 0) cycle
        pbar(nm2,m) = p1i
        if (nmm .eq. 1) cycle
        cost = dcos(theta)
        do n = 0,nmm-1
            pb(n) = -cost
            td(n) = sf(n,m)
        enddo
        if (abs(p0i) .ge. abs(p1i)) then
            pb(0) = p0i
            r = -td(0)*pb(0)
            call tridiag(nmm-1,r,td(0),pb(1),td(1))
        else
            pb(0) = p0i
            pb(1) = p1i
            r = -td(1)*pb(1)
            call tridiag(nmm-2,r,td(1),pb(2),td(2))
        endif
        do n = m,nmax-1
            i = nmax-n-1
            pbar(n,m) = pb(i)
        enddo
    enddo

    ! Vector Spherical Harmonic Basis Functions from Pbar

    do n = 0,maxn
        vbar(n,0) = -pbar(n,1)
    enddo

    do m = 1,maxo
        mm1 = m - 1
        mp1 = m + 1
         do n = m,maxn
            nm1 = n-1
            vbar(n,m) = a(n,m)*pbar(n,mp1) + b(n,m)*pbar(n,mm1)
            wbar(n,m) = c(n,m)*pbar(nm1,mm1) + d(n,m)*pbar(nm1,mp1)
      enddo
    enddo

    return

contains

    ! --------------------------------------------------------------
    ! Function for the normalized associated legendre polynomial 
    ! along the diagonal from a recursion relation.
    ! --------------------------------------------------------------

    real(8) function lfpt(n,m,theta,cp)
        implicit none
        integer,intent(in)      :: n,m
        real(8),intent(in)      :: theta
        real(8),intent(in)      :: cp(n)
        real(8)                 :: cdt,sdt
        real(8)                 :: ct,st,cth
        integer                 :: kdo,k
        lfpt = 0.0d0
        if (m .gt. n) return
        if (n .le. 0 .and. m .le. 0) then
            lfpt = dsqrt(0.5d0)
            return
        endif
        cdt = dcos(theta+theta)
        sdt = dsin(theta+theta)
        if (mod(n,2) .le. 0) then
            ct = 1.0d0
            st = 0.0d0
            if (mod(m,2) .le. 0) then
                kdo = n/2+1
                lfpt = 0.5d0*cp(1)
                do k = 2,kdo
                    cth = cdt*ct - sdt*st
                    st = sdt*ct + cdt*st
                    ct = cth
                    lfpt = lfpt + cp(k)*ct
                enddo
            else
                kdo = n/2
                do k = 1,kdo
                    cth = cdt*ct - sdt*st
                    st = sdt*ct + cdt*st
                    ct = cth
                    lfpt = lfpt + cp(k)*st
                enddo
            endif
        else
            kdo = (n+1)/2
            ct = dcos(theta)
            st = -dsin(theta)
            if (mod(m,2) .le. 0) then
                do k = 1,kdo
                    cth = cdt*ct - sdt*st
                    st = sdt*ct + cdt*st
                    ct = cth
                    lfpt = lfpt + cp(k)*ct
                enddo
            else
                do k = 1,kdo
                    cth = cdt*ct - sdt*st
                    st = sdt*ct + cdt*st
                    ct = cth
                    lfpt = lfpt + cp(k)*st
                enddo
            endif
        endif
        return
    end function lfpt

    ! -----------------------------
    !  Tri-diagonal matrix solver 
    ! -----------------------------

    subroutine tridiag(n,r,a,b,c)
        implicit none
        integer,intent(in)      :: n
        real(8),intent(in)      :: r,a(n)
        real(8),intent(inout)   :: b(n),c(n)
        real(8)                 :: bih,bih1,b1
        real(8)                 :: qih,q2,ratio,rih
        integer                 :: i,j
        select case (n)
            case(:0)
                return
            case(1)
                b(1) = r/b(1)
                return
            case(2)
                qih = a(2)
                bih = b(2)
            case(3:)
                qih = a(n)
                bih = b(n)
                do j = 3,n
                    i = n - j + 2
                    if (abs(bih) .ge. abs(c(i))) then
                        ratio = c(i)/bih
                        c(i) = 0.0d0
                        b(i+1) = qih/bih
                        bih = b(i) - ratio*qih
                        qih = a(i)
                    else
                        b(i+1) = b(i)/c(i)
                        c(i) = a(i)/c(i)
                        bih1 = qih - bih*b(i+1)
                        qih = -bih*c(i)
                        bih = bih1
                    endif
                enddo
        end select
        if (abs(bih) .ge. abs(c(1))) then
            q2 = qih/bih
            bih = b(1) - c(1)/bih*qih
            b(1) = r/bih
            b(2) = -q2*b(1)
        else
            ratio = bih/c(1)
            bih = qih - ratio*b(1)
            rih = -ratio*r
            b1 = rih/bih
            b(2) = (r - b(1)*b1)/c(1)
            b(1) = b1
        endif
        if (n-3 .ge. 0) then
            do i = 3,n
                b(i) = -b(i)*b(i-1) - c(i-1)*b(i-2)
            enddo
        endif
        return
    end subroutine tridiag

end subroutine vshbasis
    
! ========================================================================
! This subroutine intializes the scaling and recursion coeffiecents
! for normalized associated vector spherical harmonic basis funcions.
! ========================================================================

subroutine vshengineinit()

    use vshengine
    implicit none
    real(8)             :: fnmm,fnpn,fnpm
    real(8)             :: nfac1,nfac2
    integer             :: m,n,npm,nmm,nm1,nm2,i
    
    ! Factors for the normalized associated legendre polynomial

    sf = 0.0d0
    nm1 = nmax - 1
    nm2 = nmax - 2
    do m = 0,mmax
        call alfk(nm1,m,p0(:,m))
        call alfk(nm2,m,p1(:,m))
        fnmm = dble(nmax - m)
        fnpn = dble(nmax + nm1)
        fnpm = dble(nmax + m)
        i = 0
        do n = m,nm1
            fnmm = fnmm - 1.d0
            fnpn = fnpn - 2.d0
            fnpm = fnpm - 1.d0
            sf(i,m) = dsqrt(dble(fnmm*fnpm/(fnpn*(fnpn+2.d0))))
            i = i + 1
        end do
    enddo
        
    ! Theta indepedent scaling factors for calculation of Vbar and Wbar 
    ! from Pbar.
    
    do n = 1,nmax
        nfac1 = 0.5d0/dsqrt(dble(n*(n+1)))
        nfac2 = nfac1*dsqrt(dble(2*n+1)/dble(2*n-1))
        do m = 1,n
            npm = n + m 
            nmm = n - m
            a(n,m) = -nfac1*dsqrt(dble(nmm*(npm+1)))
            b(n,m) =  nfac1*dsqrt(dble(npm*(nmm+1)))
            c(n,m) =  nfac2*dsqrt(dble(npm*(npm-1)))
            d(n,m) =  nfac2*dsqrt(dble(nmm*(nmm-1)))
        enddo
    enddo

    return

contains

    ! -------------------------------------------------------------
    ! For a given M and N, the unnormalized Legendre polynomial 
    ! are calculated for all cases where (m .ge. 0).
    ! ------------------------------------------------------------

    subroutine alfk(n,m,cp)
        implicit none
        integer,intent(in)      :: n,m
        real(8),intent(out)     :: cp(n)
        real(8)                 :: fnum,fden,fnmh,fnnp1,fnmsq,fk
        real(8)                 :: pm1,t1,t2,cp2,a1,b1,c1
        integer                 :: nmms2,l,i
        if (m .gt. n) then
            cp(1) = 0.0d0
            return
        endif
        if (n .le. 0) then
            cp(1) = dsqrt(2.d0)
            return
        endif
        if (n .eq. 1) then
            if (m .eq. 0) then
                cp(1) = dsqrt(1.5d0)
            else
                cp(1) = dsqrt(0.75d0)
            endif 
            return
        endif
        if (mod(n+m,2) .eq. 0) then
            nmms2 = (n - m)/2
            fnum = dble(n + m + 1)
            fnmh = dble(n - m + 1)
            pm1 = 1.0d0
        else
            nmms2 = (n - m - 1)/2
            fnum = dble(n + m + 2)
            fnmh = dble(n - m + 2)
            pm1 = -1.0d0
        endif 
        t1 = 1.0d0
        t2 = 1.0d0
        if (nmms2 .ge. 1) then
            fden = 2.0d0
            do i = 1,nmms2
                t1 = fnum*t1/fden
                fnum = fnum + 2.0d0
                fden = fden + 2.0d0
            enddo
        endif
        if (m .ne. 0) then
            do i = 1,m
                t2 = fnmh*t2/(fnmh + pm1)
                fnmh = fnmh + 2.0d0
            enddo
        endif
        if (mod(m/2,2) .ne. 0) t1 = -t1
        cp2 = t1 * dsqrt( (dble(n) + 0.5d0) *t2 ) / (2.0d0**(n-1))
        fnnp1 = dble(n*(n+1))
        fnmsq = fnnp1 - 2.0d0*dble(m*m)
        l = (n+1)/2
        if (mod(n,2) .eq. 0 .and. mod(m,2) .eq. 0) l = l + 1
        cp(l) = cp2
        if (l .le. 1) return
        fk = dble(n)
        a1 = (fk - 2.0d0)*(fk - 1.0d0) - fnnp1
        b1 = 2.0d0*(fk*fk - fnmsq)
        cp(l-1) = b1*cp(l)/a1
        l = l - 1
        do while (l .gt. 1)
            fk = fk - 2.0d0
            a1 = (fk - 2.0d0)*(fk - 1.0d0) - fnnp1
            b1 = -2.0d0*(fk*fk - fnmsq)
            c1 = (fk + 1.0d0)*(fk + 2.0d0) - fnnp1
            cp(l-1) = -(b1*cp(l) + c1*cp(l+1))/a1
            l = l - 1
        enddo
        return
     end subroutine alfk
    
end subroutine vshengineinit
