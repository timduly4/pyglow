!
!  Code used in HWM07 for computing quasi-dipole magnetic 
!  coordinates. Modified from apex.zip code package obtained 
!  from the NSF-CEDAR database. 
!  See readme.txt file for detailed release notes.
!
!  AUTHOR
!    A. D. Richmond (NCAR-HAO)
!    Adapted for HWM07 by Douglas P. Drob
!
! Point of Contact
!    msishwmhelp@nrl.navy.mil
!
!  DATE
!    19 August 2008
!
!  REFERENCE
!    Richmond, A. D. (1995), Ionospheric electrodynamics using 
!    Magnetic Apex Coordinates, J. Geomag. & Geoelec., 47, 
!    191пїЅ212.
!
! ===============================================================

Module apexcord

    implicit none
    
    real(4),parameter        :: xmiss=-32767.0
    real(4),parameter        :: glatlim=89.9
    real(4),parameter        :: precise=7.6e-11
    real(4),parameter        :: datdmx=1.
    real(4),parameter        :: datimx=2.5
    integer(4),parameter     :: irlf=4

    integer(4)               :: kgma = 0
    real(4)                  :: glalmn,glalmx
    integer(4)               :: nla,nlo,nal
    integer(4)               :: lbx,lby,lbz,lbv
    integer(4)               :: lla,llo,lal

    integer(4)               :: lwk
    real(4),allocatable      :: wk(:)

    real(4)                  :: colat
    real(4)                  :: elon
    real(4)                  :: vp
    real(4)                  :: ctp
    real(4)                  :: stp
      
    real(4)                  :: rtod
    real(4)                  :: dtor
    real(4)                  :: re
    real(4)                  :: req
    real(4)                  :: msgu
    real(4)                  :: pola

    integer(4)               :: io = 1
    integer(4)               :: jo = 1
    integer(4)               :: ko = 1

    logical                  :: loaddata = .true.

end module apexcord
 
! ===============================================================
! Initalize the apex coordinate module variables
! ===============================================================

subroutine apxrda()

    use apexcord
    implicit none
    
    integer(4) :: i

    open(unit=77,file='apexgrid.dat',form='unformatted')
    
    read(77) kgma,glalmn,glalmx,nla,nlo,nal, &
        lbx,lby,lbz,lbv,lla,llo,lal
        
    read(77) colat,elon,vp,ctp,stp
    read(77) rtod,dtor,re,req,msgu,pola

    lwk = nla*nlo*nal*5 + nla + nlo + nal
    
    allocate(wk(lwk))
    
    read(77) (wk(i),i=1,lwk)
    close(77)
    
    loaddata = .false.
    
    return

end subroutine apxrda

! =======================================================================
! Convert from (glat,glon) to apex coordinates
! =======================================================================

subroutine apex(glat,glon,alt,hr,alon,xlatqd,f1,f2,ist)
     
    use apexcord
    implicit none

    real(4),intent(in)  :: glat,glon,alt,hr

    real(4)        :: b(3),bhat(3)
    real(4)        :: bmag,si
    real(4)        :: alon
    real(4)        :: xlatm,vmp,w,d
    real(4)        :: be3,sim
    real(4)        :: xlatqd
    real(4)        :: d1(3),d2(3),d3(3)
    real(4)        :: e1(3),e2(3),e3(3)
    real(4)        :: f1(2),f2(2)
    
    integer(4)     :: ist
  
    real(4)        :: dfxdth,dfydth,dfzdth,dfvdth,dfxdln
    real(4)        :: dfydln,dfzdln,dfvdln,dfxdh,dfydh,dfzdh,dfvdh
    
    real(4)        :: fx,fy,fz,fv
    
    real(4)        :: cth,sth
    real(4)        :: gradx(3),grady(3),gradz(3),gradv(3)
    real(4)        :: glatx
    
    real(4)        :: fxdum,fydum,fzdum,fvdum
    real(4)        :: dmxdth,dmydth,dmzdth,dmvdth
    real(4)        :: dmxdh,dmydh,dmzdh,dmvdh
    
    real(4)        :: clm,r3_2
    real(4)        :: grclm(3),clmgrp(3),rgrlp(3)
        
    real(4)        :: f
    
! Initalize the module if needed
    
    if (loaddata) then
        call apxrda()
        loaddata = .false.
    endif
                   
! Interpolate the 3d data cube
                   
    call intrp(glat,glon,alt,wk(lbx),wk(lby),wk(lbz),wk(lbv), &
               nla,nlo,nal,wk(lla),wk(llo),wk(lal),fx,fy,fz,fv, &
               dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln, &
               dfvdln,dfxdh,dfydh,dfzdh,dfvdh,ist)

    call adpl(glat,glon,cth,sth,fx,fy,fz,fv, &
         dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
     
    call gradxyzv(alt,cth,sth,  &
        dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
        dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)

    ! If the point is very close to either the North or South
    ! geographic pole, recompute the east-west gradients after
    ! stepping a small distance from the pole.
    
    if (glat .gt. glalmx .or. glat .lt. glalmn) then
            
        glatx = glalmx
          
        if (glat .lt. 0.0) glatx = glalmn
          
        call intrp(glatx,glon,alt,wk(lbx),wk(lby),wk(lbz),wk(lbv), &
              nla,nlo,nal,wk(lla),wk(llo),wk(lal),fxdum,fydum,fzdum,fvdum, &
              dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln, &
              dfvdln,dmxdh,dmydh,dmzdh,dmvdh,ist)
          
        call adpl(glatx,glon,cth,sth,fxdum,fydum,fzdum,fvdum,&
            dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,dfvdln)
         
        call grapxyzv(alt,cth,sth,dfxdln, &
            dfydln,dfzdln,dfvdln,gradx,grady,gradz,gradv)
      
    endif

    call gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
        xlatm,alon,vmp,grclm,clmgrp,xlatqd,rgrlp,b,clm,r3_2)
     
    call basvec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
        bmag,sim,si,f,d,w,bhat,d1,d2,d3,e1,e2,e3,f1,f2)

    be3 = bmag/d

    ist = 0

    return

end subroutine apex

! =========================================================================
! Interpolation subroutine
! =========================================================================

subroutine intrp(glat,glon,alt,x,y,z,v,nlat,nlon,nalt, &
             gplat,gplon,gpalt,fx,fy,fz,fv, &
             dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln, &
             dfvdln,dfxdh,dfydh,dfzdh,dfvdh,ist)
     
    use apexcord
    implicit none

    real(4)        :: glat,glon,alt
    integer(4)     :: nlat,nlon,nalt
    real(4)        :: x(nlat,nlon,nalt)
    real(4)        :: y(nlat,nlon,nalt)
    real(4)        :: z(nlat,nlon,nalt)
    real(4)        :: v(nlat,nlon,nalt)
    real(4)        :: gplat(nlat),gplon(nlon),gpalt(nalt)
    real(4)        :: fx,fy,fz,fv
    real(4)        :: dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln, &
                   dfvdln,dfxdh,dfydh,dfzdh,dfvdh
    integer(4)     :: ist

    
    real(4)        :: dlat,dlon
    real(4)        :: xi,yj
    real(4)        :: hti,diht
    real(4)        :: zk

    real(4)        :: dfxdn,dfxde,dfxdd
    real(4)        :: dfydn,dfyde,dfydd
    real(4)        :: dfzdn,dfzde,dfzdd
    real(4)        :: dfvdn,dfvde,dfvdd
    real(4)        :: dmf,dmdfdn,dmdfde,dmdfdd
    
    real(4)        :: fac
    real(4)        :: omfac

    integer(4)     :: i,j,k

    ist = 0

! input bounds checks

    if (glat .lt. gplat(1) .or. glat .gt. gplat(nlat)) stop 'glat out of bounds'
    if (alt  .lt. gpalt(1) .or. alt  .gt. gpalt(nalt)) stop 'alt out of bounds'
    if (glon .lt. gplon(1))    glon = glon + 360.0
    if (glon .gt. gplon(nlon)) glon = glon - 360.0
    if (glon .lt. gplon(1) .or. glon .gt. gplon(nlon)) stop 'glon out of bounds'

! Locate the closest indicies of the points in the array and
! and calculate the fractional offset for linear interpolation

    ! Find I
    
    i = io
    if (gplat(i+1) .lt. glat) then
        do while (gplat(i+1) .lt. glat)
            i = i + 1
        enddo
    endif
    if (gplat(i) .gt. glat) then
        do while (gplat(i) .gt. glat)
            i = i - 1
        enddo
    endif
    io = i
    dlat = gplat(i+1) - gplat(i)
    xi = (glat - gplat(i)) / dlat
        
    ! Find J

    j = jo
    if (gplon(j+1) .lt. glon) then
        do while (gplon(j+1) .lt. glon)
            j = j + 1
        enddo
    endif
    if (gplon(j) .gt. glon) then
        do while (gplon(j) .gt. glon)
            j = j - 1
        enddo
    endif
    jo = j
    dlon = gplon(j+1) - gplon(j)
    yj = (glon - gplon(j)) / dlon

    ! Find K

    k = ko
    if (gpalt(k+1) .lt. alt) then
        do while (gpalt(k+1) .lt. alt)
            k = k + 1
        enddo
    endif
    if (gpalt(k) .gt. alt) then
        do while (gpalt(k) .gt. alt)
            k = k - 1
        enddo
    endif
    ko = k
   
    hti  = re/(re+alt)
    diht = re/(re+gpalt(k+1)) - re/(re+gpalt(k))
    zk  = (hti - re/(re+gpalt(k))) / diht
  
!  For intrp:

    call trilin(x(i,j,k),nlat,nlon,xi,yj,zk,fx,dfxdn,dfxde,dfxdd)
    dfxdth = -dfxdn*rtod/dlat
    dfxdln =  dfxde*rtod/dlon
    dfxdh  = -hti*hti*dfxdd/(re*diht)

    call trilin(y(i,j,k),nlat,nlon,xi,yj,zk,fy,dfydn,dfyde,dfydd)
    dfydth = -dfydn*rtod/dlat
    dfydln =  dfyde*rtod/dlon
    dfydh  = -hti*hti*dfydd/(re*diht)

    call trilin(z(i,j,k),nlat,nlon,xi,yj,zk,fz,dfzdn,dfzde,dfzdd)
    dfzdth = -dfzdn*rtod/dlat
    dfzdln =  dfzde*rtod/dlon
    dfzdh  = -hti*hti*dfzdd/(re*diht)

    call trilin(v(i,j,k),nlat,nlon,xi,yj,zk,fv,dfvdn,dfvde,dfvdd)
    dfvdth = -dfvdn*rtod/dlat
    dfvdln =  dfvde*rtod/dlon
    dfvdh  = -hti*hti*dfvdd/(re*diht)

! Improve calculation of longitudinal derivatives near poles

    if (glat .lt. dlat - 90.0) then
        fac = 0.5*xi
        omfac = 1.0 - fac
        xi = xi - 1.0
        i = i + 1
        call trilin(x(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
        dfxdln = dfxdln*omfac + fac*dmdfde*rtod/dlon
        call trilin(y(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
        dfydln = dfydln*omfac + fac*dmdfde*rtod/dlon
        call trilin(v(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
        dfvdln = dfvdln*omfac + fac*dmdfde*rtod/dlon
        return
    endif

    if (glat .gt. 90.0-dlat) then
        fac = 0.5*(1.0- xi)
        omfac = 1.0 - fac
        xi = xi + 1.0
        i = i - 1
        call trilin(x(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
        dfxdln = dfxdln*omfac + fac*dmdfde*rtod/dlon
        call trilin(y(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
        dfydln = dfydln*omfac + fac*dmdfde*rtod/dlon
        call trilin(v(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
        dfvdln = dfvdln*omfac + fac*dmdfde*rtod/dlon
    endif
      
    return

end subroutine intrp

! ============================================================================
!  Trilinear interpolation of u and its derivatives
! 940803 A. D. Richmond
! Inputs:
!   u(1,1,1) = address of lower corner of interpolation box 
!   nlat = first dimension of u from calling routine
!   nlon = second dimension of u from calling routine
!   xi = fractional distance across box in x direction 
!   yj = fractional distance across box in y direction 
!   zk = fractional distance across box in z direction 
! Outputs:
!   fu = interpolated value of u
!   dfudx = interpolated derivative of u with respect to i (x direction)
!   dfudy = interpolated derivative of u with respect to j (y direction)
!   dfudz = interpolated derivative of u with respect to k (z direction)
! ============================================================================

subroutine trilin(u,nlat,nlon,xi,yj,zk,fu,dfudx,dfudy,dfudz)

    implicit none

    integer(4) :: nlat,nlon
    real(4)    :: u(nlat,nlon,2)

    real(4)    :: xi,yj,zk
    real(4)    :: fu
    real(4)    :: dfudx,dfudy,dfudz
    real(4)    :: omxi,omyj,omzk

    omxi = 1.0 - xi
    omyj = 1.0 - yj
    omzk = 1.0 - zk

    fu = u(1,1,1)*omxi*omyj*omzk &
        + u(2,1,1)*xi*omyj*omzk &
        + u(1,2,1)*omxi*yj*omzk &
        + u(1,1,2)*omxi*omyj*zk &
        + u(2,2,1)*xi*yj*omzk &
        + u(2,1,2)*xi*omyj*zk &
        + u(1,2,2)*omxi*yj*zk &
        + u(2,2,2)*xi*yj*zk

    dfudx = (u(2,1,1)-u(1,1,1))*omyj*omzk + (u(2,2,1)-u(1,2,1))*yj*omzk &
        + (u(2,1,2)-u(1,1,2))*omyj*zk &
        + (u(2,2,2)-u(1,2,2))*yj*zk
        
    dfudy = (u(1,2,1)-u(1,1,1))*omxi*omzk &
        + (u(2,2,1)-u(2,1,1))*xi*omzk &
        + (u(1,2,2)-u(1,1,2))*omxi*zk &
        + (u(2,2,2)-u(2,1,2))*xi*zk
        
    dfudz = (u(1,1,2)-u(1,1,1))*omxi*omyj &
        + (u(2,1,2)-u(2,1,1))*xi*omyj &
        + (u(1,2,2)-u(1,2,1))*omxi*yj &
        + (u(2,2,2)-u(2,2,1))*xi*yj
        
    return
      
end subroutine trilin

! =============================================================================
!C  v is used for vr2n
!C  Add-back of pseudodipole component to x,y,z,v and their derivatives.
!C  940715 A. D. Richmond
!C Inputs:
!C   glat = latitude, degrees
!C   glon = longitude, degrees
!C   fx = interpolated value of x
!C   fy = interpolated value of y
!C   fz = interpolated value of z
!C   fv = interpolated value of v
!C   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
!C	colatitude, in radians-1
!C   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
!C	longitude, in radians-1
!C Output:
!C   cth,sth = cos(colatitude), sin(colatitude)
!C   fx = interpolated value of x
!C   fy = interpolated value of y
!C   fz = interpolated value of z
!C   fv = interpolated value of v
!C   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
!C	colatitude, in radians-1
!C   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
!C	longitude, in radians-1
! ==============================================================================

subroutine adpl(glat,glon,cth,sth,fx,fy,fz,fv, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)

    use apexcord
    implicit none
    
    real(4)    :: glat,glon
    real(4)    :: cth,sth
    real(4)    :: fx,fy,fz,fv
    real(4)    :: dfxdth,dfydth,dfzdth,dfvdth
    real(4)    :: dfxdln,dfydln,dfzdln,dfvdln
    
    real(4)    :: cph,sph
    real(4)    :: ctm

    cph = cos((glon-elon)*dtor)
    sph = sin((glon-elon)*dtor)
    cth = sin(glat*dtor)
    sth = cos(glat*dtor)
    ctm = ctp*cth + stp*sth*cph

    fx = fx + sth*ctp*cph - cth*stp
    fy = fy + sth*sph
    fz = fz + ctm
    fv = fv - ctm

    dfxdth = dfxdth + ctp*cth*cph + stp*sth
    dfydth = dfydth + cth*sph
    dfzdth = dfzdth - ctp*sth + stp*cth*cph
    dfvdth = dfvdth + ctp*sth - stp*cth*cph
    dfxdln = dfxdln - ctp*sth*sph
    dfydln = dfydln + sth*cph
    dfzdln = dfzdln - stp*sth*sph
    dfvdln = dfvdln + stp*sth*sph

    return
end subroutine adpl

!*******************************************************************************
!          Calculates east, north, and up components of gradients of x,y,z,v in
!          geodetic coordinates.  All gradients are in inverse km.  Assumes
!          flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
!         940803 A. D. Richmond
!
!
!          40680925. = req**2 (rounded off)
!          272340.   = req**2 * e2, where e2 = (2. - 1./298.25)/298.25
!                      is square of eccentricity of ellipsoid.
!
! ==============================================================================

subroutine gradxyzv(alt,cth,sth, &
            dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
            dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)
            
            
    implicit none  
    
    real(4)    :: alt,cth,sth
    real(4)    :: dfxdth,dfydth,dfzdth,dfvdth
    real(4)    :: dfxdln,dfydln,dfzdln,dfvdln
    real(4)    :: dfxdh,dfydh,dfzdh,dfvdh
    real(4)    :: gradx(3),grady(3),gradz(3),gradv(3)

    real(4)    :: d2,d,rho
    real(4)    :: dddthod,drhodth,dzetdth,ddisdth

    d2 = 40680925.e0 - 272340.e0*cth*cth
    d = sqrt(d2)
    rho = sth*(alt + 40680925.e0/d)
    dddthod = 272340.e0*cth*sth/d2
    drhodth = alt*cth + (40680925.e0/d)*(cth-sth*dddthod)
    dzetdth =-alt*sth - (40408585.e0/d)*(sth+cth*dddthod)
    ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)
    
    gradx(1) = dfxdln/rho
    grady(1) = dfydln/rho
    gradz(1) = dfzdln/rho
    gradv(1) = dfvdln/rho

    gradx(2) = -dfxdth/ddisdth
    grady(2) = -dfydth/ddisdth
    gradz(2) = -dfzdth/ddisdth
    gradv(2) = -dfvdth/ddisdth
    
    gradx(3) = dfxdh
    grady(3) = dfydh
    gradz(3) = dfzdh
    gradv(3) = dfvdh

    return
    
end subroutine gradxyzv

!*******************************************************************************
!          Calculates east, north, and up components of gradients of x,y,z,v in
!          geodetic coordinates.  All gradients are in inverse km.  Assumes
!          flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
!         940803 A. D. Richmond
!
!          40680925. = req**2 (rounded off)
!          272340.   = req**2 * e2, where e2 = (2. - 1./298.25)/298.25
!                      is square of eccentricity of ellipsoid.
!
! ==============================================================================

subroutine grapxyzv(alt,cth,sth, &
        dfxdln,dfydln,dfzdln,dfvdln,gradx,grady,gradz,gradv)

    implicit none

    real(4)    :: alt,cth,sth
    real(4)    :: dfxdln,dfydln,dfzdln,dfvdln
    real(4)    :: gradx(3),grady(3),gradz(3),gradv(3)

    real(4)    :: d2,d,rho
    real(4)    :: dddthod,drhodth,dzetdth,ddisdth


    d2 = 40680925.e0 - 272340.e0*cth*cth  
    d = sqrt(d2)
    rho = sth*(alt + 40680925.e0/d)
    dddthod = 272340.e0*cth*sth/d2
    drhodth = alt*cth + (40680925.e0/d)*(cth-sth*dddthod)
    dzetdth =-alt*sth - (40408585.e0/d)*(sth+cth*dddthod)
    ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)
    
    gradx(1) = dfxdln/rho
    grady(1) = dfydln/rho
    gradz(1) = dfzdln/rho
    gradv(1) = dfvdln/rho

    return

end subroutine grapxyzv
      
!*******************************************************************************
!          Uses gradients of x,y,z,v to compute geomagnetic field and
!          gradients of apex latitude, longitude.
!          940819 A. D. Richmond
!          INPUT:
!            HR     = reference altitude
!            ALT    = altitude
!            FX,FY,FZ,FV = interpolated values of x,y,z,v, plus pseudodipole
!                     component
!            GRADX,GRADY,GRADZ,GRADV = interpolated gradients of x,y,z,v,
!                     including pseudodipole components (east,north,up)
!          OUTPUT:
!            XLATM  = modified apex latitude (lambda_m), degrees
!            XLONM  = apex longitude (phi_a), degrees
!            VMP    = magnetic potential, in T.m.
!            GRCLM  = grad(cos(lambda_m)), in km-1
!            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
!            QDLAT  = quasi-dipole latitude, degrees
!            RGRLP  = (re + alt)*grad(lambda')
!            B      = magnetic field, in nT
!            CLM    = cos(lambda_m)
!            R3_2   = ((re + alt)/(re + hr))**(3/2)
! ===========================================================================
      
subroutine gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
            xlatm,xlonm,vmp,grclm,clmgrp,qdlat,rgrlp,b,clm,r3_2)

    use apexcord
    implicit none

    real(4)    :: hr,alt
    real(4)    :: fx,fy,fz,fv
    real(4)    :: gradx(3),grady(3),gradz(3),gradv(3)
    real(4)    :: xlatm,xlonm,vmp
    real(4)    :: grclm(3),clmgrp(3),rgrlp(3)
    real(4)    :: b(3)
    real(4)    :: r3_2

    real(4)    :: rr,r,rn
    real(4)    :: sqrror
    
    real(4)    :: cpm,spm
    real(4)    :: bo
    real(4)    :: rn2
    real(4)    :: x2py2
    real(4)    :: xnorm
    real(4)    :: xlp,slp,clp
    real(4)    :: qdlat,clm
    real(4)    :: grclp

    integer(4) :: i

    rr = re + hr
    r  = re + alt
    rn = r/re
    sqrror = sqrt(rr/r)
    r3_2 = 1.0/sqrror/sqrror/sqrror
    xlonm = atan2(fy,fx)
    cpm = cos(xlonm)
    spm = sin(xlonm)
    xlonm = rtod*xlonm
    bo = vp*1.e6

! 1.E6 converts T to nT and km-1 to m-1.

    rn2 = rn*rn
    vmp = vp*fv/rn2
    b(1) = -bo*gradv(1)/rn2
    b(2) = -bo*gradv(2)/rn2
    b(3) = -bo*(gradv(3) - 2.0*fv/r)/rn2

    x2py2 = fx*fx + fy*fy
    xnorm = sqrt(x2py2 + fz*fz)
    xlp = atan2(fz,sqrt(x2py2))
    slp = sin(xlp)
    clp = cos(xlp)
    qdlat = xlp*rtod
    clm = sqrror*clp

    if (clm .gt. 1.0) stop &
        'gradlpv, point lies below field line that peaks at reference height.'
                
    xlatm = rtod*acos(clm)
    
!  If southern magnetic hemisphere, reverse sign of xlatm

    if (slp .lt. 0.0) xlatm = - xlatm

    do i = 1,3
        grclp = cpm*gradx(i) + spm*grady(i)
        rgrlp(i) = r*(clp*gradz(i) - slp*grclp)
        grclm(i) = sqrror*grclp
        clmgrp(i) = sqrror*(cpm*grady(i)-spm*gradx(i))
    enddo
      
    grclm(3) = grclm(3) - sqrror*clp/(2.0*r)

    return

end subroutine gradlpv

!*******************************************************************************
!          Computes base vectors and other parameters for apex coordinates.
!          Vector components:  east, north, up
!          940801 A. D. Richmond
!          Reference:
!            Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
!            Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
!          INPUTS:
!            HR     = reference altitude
!            XLATM  = modified apex latitude, degrees
!            GRCLM  = grad(cos(lambda_m)), in km-1
!            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
!            RGRLP  = (re + altitude)*grad(lambda')
!            B      = magnetic field, in nT
!            CLM    = cos(lambda_m)
!            R3_2   = ((re + altitude)/(re + hr))**(3/2)
!          RETURNS:
!            BMAG    = magnitude of magnetic field, in nT
!            SIM     = sin(I_m) of article
!            SI      = sin(I)
!            F       = F of article
!            D       = D of article
!            W       = W of article
!            BHAT    = unit vector along geomagnetic field direction
!            D1...F2 = base vectors of article
!
!=============================================================================

subroutine basvec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
            bmag,sim,si,f,d,w,bhat,d1,d2,d3,e1,e2,e3,f1,f2)

    use apexcord
    implicit none

    real(4)    :: hr,xlatm
    real(4)    :: grclm(3),clmgrp(3),rgrlp(3),b(3)
    real(4)    :: clm,r3_2
    real(4)    :: bmag
    real(4)    :: sim
    real(4)    :: si
    real(4)    :: f
    real(4)    :: d
    real(4)    :: w
    real(4)    :: bhat(3)
    real(4)    :: d1(3),d2(3),d3(3)
    real(4)    :: e1(3),e2(3),e3(3)
    real(4)    :: f1(2),f2(2)
    
    real(4)    :: rr
    real(4)    :: simoslm
    
    integer(4) :: i
    real(4)    :: d1db,d2db
    
    rr = re + hr
    simoslm = 2.0/sqrt(4.0 - 3.0*clm*clm)
    sim = simoslm*sin(xlatm*dtor)
    bmag = sqrt(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))
    d1db = 0.0
    d2db = 0.0
    do i = 1,3
        bhat(i) = b(i)/bmag
        d1(i) = rr*clmgrp(i)
        d1db = d1db + d1(i)*bhat(i)
        d2(i) = rr*simoslm*grclm(i)
        d2db = d2db + d2(i)*bhat(i)
    enddo

!   Ensure that d1,d2 are exactly perpendicular to B:

    do i = 1,3
        d1(i) = d1(i) - d1db*bhat(i)
        d2(i) = d2(i) - d2db*bhat(i)
    enddo  

    e3(1) = d1(2)*d2(3) - d1(3)*d2(2)
    e3(2) = d1(3)*d2(1) - d1(1)*d2(3)
    e3(3) = d1(1)*d2(2) - d1(2)*d2(1)
    d = bhat(1)*e3(1) + bhat(2)*e3(2) + bhat(3)*e3(3)

! Following step may be unnecessary, but it ensures that e3 lies along bhat.
    
    do i = 1,3
        d3(i) = bhat(i)/d
        e3(i) = bhat(i)*d
    enddo
   
    e1(1) = d2(2)*d3(3) - d2(3)*d3(2)
    e1(2) = d2(3)*d3(1) - d2(1)*d3(3)
    e1(3) = d2(1)*d3(2) - d2(2)*d3(1)
    e2(1) = d3(2)*d1(3) - d3(3)*d1(2)
    e2(2) = d3(3)*d1(1) - d3(1)*d1(3)
    e2(3) = d3(1)*d1(2) - d3(2)*d1(1)

    w = rr*rr*clm*abs(sim)/(bmag*d)

    si = -bhat(3)

    f1(1) =  rgrlp(2) 
    f1(2) = -rgrlp(1)
    f2(1) = -d1(2)*r3_2
    f2(2) =  d1(1)*r3_2
    f = f1(1)*f2(2) - f1(2)*f2(1)

    return

end subroutine basvec
