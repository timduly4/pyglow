!
!  Disturbance wind part of Horizontal Wind Model HWM07
!  Version DWM07B104i
!  See readme.txt file for detailed release notes.
!
!  AUTHOR
!    John Emmert
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
!    Emmert, J. T., D. P. Drob, G. G. Shepherd, G. Hernandez, M. J. Jarvis, J. W. 
!      Meriwether, R. J. Niciejewski, D. P. Sipler, and C. A. Tepley (2008),
!      DWM07 global empirical model of upper thermospheric storm-induced 
!      disturbance winds, J. Geophys Res., 113, doi:10.1029/2008JA013541.
!


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                           Module DWM
!
! Description: This is a common data module for model definition.  These
!  parameters are set by the first calling the subroutine loaddwm().
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module dwm_module

    implicit none

    integer(4)             :: nterm          ! Number of terms in the model
    integer(4)             :: lmax           ! Max latitudinal degree
    integer(4)             :: mmax           ! Max order of MLT var.

    integer(4),allocatable :: termarr(:,:)   ! 3 x nterm index of coupled terms
    real(4),allocatable    :: coeff(:)       ! Model coefficients
    real(4),allocatable    :: vsh_terms(:,:) ! VSH basis values
    integer(4)             :: nvshfn         ! # of VSH basis functions
    real(4)                :: twidth         ! Transition width of high-lat mask
    real(4),allocatable    :: termval(:,:)   ! Term values to which coefficients are applied

    ! Initialization flags and information

    logical                 :: modelinit = .true.
    character(128)          :: defaultdata = 'dwm07b_104i.dat'

end module dwm_module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                       DWM-HWM Interface
!
!  Author: John Emmert
!  Date: 7/10/2007
!  Description: Using HWM inputs, computes Quasi-dipole latitude and local time,
!               and Kp. Retrieves DWM results for these conditions, converts
!               to geographic directions, and applies artificial height profile.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dwm07b_hwm_interface(IYD,SEC,ALT,GLAT,GLON,AP,DW)

    implicit none

    INTEGER,intent(in)      :: IYD
    REAL(4),intent(in)      :: SEC,ALT,GLAT,GLON
    REAL(4),intent(in)      :: AP(2)
    REAL(4),intent(out)     :: DW(2)

    real(4)                 :: ut, mlat, mlt, kp, mmpwind, mzpwind
    real(4)                 :: f1e, f1n, f2e, f2n
    real(4)                 :: dummy
    real(4)                 :: day, mlon, asun_glat, asun_glon, asun_mlat, asun_mlon
    real(4), external       :: ap_to_kp, dwm_altwgt

    real(8), parameter      :: pi=3.141592653590, dtor=pi/180D0, sin_eps=0.39781868

    !CONVERT AP TO KP
    kp = ap_to_kp(ap(2))

    !CONVERT GEO LAT/LON TO QD LAT/LON

    call gd2qd(glat,glon,mlat,mlon,f1e,f1n,f2e,f2n)

    !COMPUTE QD MAGNETIC LOCAL TIME (LOW-PRECISION)
    day = real(mod(iyd,1000))
    ut = sec / 3600.0
    asun_glat = -asin(sin((day-80.0)*dtor) * sin_eps) / dtor
    asun_glon = -ut * 15.0
    call gd2qd(asun_glat, asun_glon, asun_mlat, asun_mlon, &
               dummy,dummy,dummy,dummy)
    mlt = (mlon - asun_mlon) / 15.0

    !RETRIEVE DWM WINDS
    call dwm07b(mlt, mlat, kp, mmpwind, mzpwind)

    !CONVERT TO GEOGRAPHIC COORDINATES
    dw(1) = f2n*mmpwind + f1n*mzpwind
    dw(2) = f2e*mmpwind + f1e*mzpwind

    !APPLY HEIGHT PROFILE
    dw = dw * dwm_altwgt(alt)

    return

end subroutine dwm07b_hwm_interface



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   This subroutine is used to evaluate DWM
!
!  Programming Notes:
!
!   This subroutine is only OPENMP/THREAD SAFE when no calls to
!   loaddwm() are made.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dwm07b(mlt, mlat, kp, mmpwind, mzpwind)

    use dwm_module
    implicit none

    real(4),intent(in)        :: mlt            !Magnetic local time (hours)
    real(4),intent(in)        :: mlat           !Magnetic latitude (degrees)
    real(4),intent(in)        :: kp             !3-hour Kp

    real(4),intent(out)       :: mmpwind        !Mer. disturbance wind (+north, QD coordinates)
    real(4),intent(out)       :: mzpwind        !Zon. disturbance wind (+east, QD coordinates)

    ! Local variables
    real(4)                   :: mltdeg
    real(4)                   :: kp_terms(0:2)
    real(4)                   :: latwgt_terms
    real(4)                   :: termval_temp(0:1)
    integer(4)                :: iterm

    external                  :: loaddwm, vsh_basis
    real(4), external         :: dwm_latwgt2

    
    !LOAD MODEL PARAMETERS IF NECESSARY
    if (modelinit) call loaddwm(defaultdata)

    !COMPUTE VSH TERMS
    mltdeg = 15.0*mlt
    call vsh_basis(mlat, mltdeg)

    !COMPUTE KP TERMS
    call dwm_kpspl3_calc(kp, kp_terms)

    !COMPUTE LATITUDINAL WEIGHTING TERMS
    latwgt_terms = dwm_latwgt2(mlat, mlt, kp, twidth)

    !GENERATE COUPLED TERMS
    do iterm = 0, nterm-1
      termval_temp = 1
      if (termarr(0,iterm) .ne. 999) termval_temp = termval_temp * vsh_terms(0:1,termarr(0,iterm))
      if (termarr(1,iterm) .ne. 999) termval_temp = termval_temp * kp_terms(termarr(1,iterm))
      if (termarr(2,iterm) .ne. 999) termval_temp = termval_temp * latwgt_terms
      termval(0:1,iterm) = termval_temp
    end do

    !APPLY COEFFICIENTS
    mmpwind = dot_product(coeff, termval(0,0:nterm-1))
    mzpwind = dot_product(coeff, termval(1,0:nterm-1))

    return

end subroutine dwm07b



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine loads the disturbance wind model parameters
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine loaddwm(datafile)

    use dwm_module, maxmd=>mmax, maxld=>lmax, modelinitd=>modelinit

    implicit none

    character(128),intent(in)   :: datafile

    external vsh_basis_init

    open(unit=23,file=trim(datafile),form='unformatted')

    if (allocated(termarr)) deallocate(termarr,coeff)
    read(23) nterm, maxmd, maxld
    allocate(termarr(0:2, 0:nterm-1))
    read(23) termarr
    allocate(coeff(0:nterm-1))
    allocate(termval(0:1, 0:nterm-1))
    read(23) coeff
    read(23) twidth

    close(23)

    call vsh_basis_init

    modelinitd = .false.

    return

end subroutine loaddwm



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                           Module vsh_basis_module
!
! Description: This is a common data module for model definition.  These
!  parameters set by the first calling the subroutine loaddwm().
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module vsh_basis_module

    implicit none

    type :: fn_map                              ! Structure of VSH basis function map
      integer(4) :: realimag
      integer(4) :: irrotational
      integer(4) :: M
      integer(4) :: N
    end type
    type (fn_map), allocatable :: fn_map1(:)    ! Working array of VSH map structure

    real(8), allocatable :: A(:,:)  !Derivative of the Associated Legendre Functions
    real(8), allocatable :: B(:,:)  !m/sin(theta) times the Associated Legendre Functions
    real(8), allocatable :: anm(:,:)
    real(8), allocatable :: bnm(:,:)
    real(8), allocatable :: cm(:)
    real(8), allocatable :: cn(:)
    real(8), allocatable :: e0n(:)
    real(8), allocatable :: fnm(:,:)
    real(8), allocatable :: m_arr(:)
    real(8), allocatable :: n_arr(:)
    real(8), allocatable :: cosmz(:)
    real(8), allocatable :: sinmz(:)

end module vsh_basis_module



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine initializes the VSH basis recursion coefficients
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine vsh_basis_init

    use dwm_module, only:m_max=>mmax, n_max=>lmax, nvshfn, vsh_terms
    use vsh_basis_module

    implicit none

    integer(4)       :: i, j, n, m, ifn
    type (fn_map)    :: fn_map0(0:2*2*(m_max+1)*(n_max+1))


    !CREATE VSH FUNCTION MAP
    ifn = -1
    do n = 0, n_max
      do m = 0, m_max
        if ((m .eq. 0) .and. (n .eq. 0)) cycle
        if (m .gt. n) cycle
        do j = 0, 1
          do i = 0, 1
            if ((m .eq. 0) .and. (i .eq. 1)) cycle
            ifn = ifn + 1
            fn_map0(ifn)%realimag     = i
            fn_map0(ifn)%irrotational = j
            fn_map0(ifn)%M            = m
            fn_map0(ifn)%N            = n
          enddo
        enddo
      enddo
    enddo
    nvshfn = ifn + 1
    allocate(fn_map1(0:nvshfn-1))
    do ifn = 0, nvshfn-1
      fn_map1(ifn)%realimag     = fn_map0(ifn)%realimag
      fn_map1(ifn)%irrotational = fn_map0(ifn)%irrotational
      fn_map1(ifn)%M            = fn_map0(ifn)%M
      fn_map1(ifn)%N            = fn_map0(ifn)%N
    end do

    !CREATE ARRAY THAT WILL CONTAIN VSH BASIS VALUES
    allocate(vsh_terms(0:1, 0:nvshfn-1))

    !CREATE RECURSION COEFFICIENT ARRAYS
    allocate( A(0:m_max, 0:n_max) )
    allocate( B(0:m_max, 0:n_max) )
    allocate( anm(0:m_max, 0:n_max) )
    allocate( bnm(0:m_max, 0:n_max) )
    allocate( fnm(0:m_max, 0:n_max) )
    allocate( cm(0:m_max) )
    allocate( cn(0:n_max) )
    allocate( e0n(0:n_max) )
    allocate( m_arr(0:m_max) )
    allocate( n_arr(0:n_max) )
    allocate( cosmz(0:m_max) )
    allocate( sinmz(0:m_max) )
    A = 0
    B = 0

    !COMPUTE RECURSION COEFFICIENTS
    do m = 0, m_max
      cm(m) =  dsqrt(1 + 0.5/dble(max(m,1)))
      m_arr(m) = dble(m)
    end do
    do n = 0, n_max
      n_arr(n) = dble(n)
      cn(n) = 1 / dsqrt(dble(max(n,1)) * dble(n+1))
      e0n(n) = dsqrt(dble(n*(n+1)) / 2.0)
    end do
    anm = 0
    bnm = 0
    fnm = 0
    do m = 1, m_max
      if (m .eq. n_max) cycle
      do n = m+1, n_max
        anm(m,n) = dsqrt( dble((2*n-1)*(2*n+1)) / dble((n-m)*(n+m)) )
        bnm(m,n) = dsqrt( dble((2*n+1)*(n+m-1)*(n-m-1)) / dble((n-m)*(n+m)*(2*n-3)) )
        fnm(m,n) = dsqrt( dble((n-m)*(n+m)*(2*n+1)) / dble(2*n-1) )
      end do
    enddo

    return

end subroutine vsh_basis_init



!***************************************************************************************************
!
!JOHN EMMERT  3/26/07
!EVALUATES A SERIES OF 2-COMPONENT (HORIZONTAL) VECTOR SPHERICAL HARMONIC FUNCTIONS, UP TO SPECIFIED
!DEGREE AND ORDER. THE CONVENTION OF KILLEEN ET AL. [ADV. SPACE RES., 1987, P. 207]
!IS USED, EXCEPT THAT NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS ARE USED.
!
!7/7/07 TRANSLATED TO F90 FROM IDL VERSION. TRANSLATION IS LIMITED, AND ASSUMES THAT THE INPUT
!ARE IN DEGREES, THE THETA ARGUMENT IS A LATITUDE, AND ONLY VALID (NON-ZERO) FUNCTIONS ARE TO BE 
!RETURNED. MUST BE INITIALIZED WITH vsh_basis_init
!
!CALLING SEQUENCE
!        Result = vsh_basis(theta,phi)
!
!ARGUMENTS
!        theta        The latitudinal coordinate, specified in degrees.
!        phi            The azimuthal coordinate, specified in degrees.
!
!RESULT
!       A 2 x nfn array, where nfn is the number of valid functions in the basis.
!        The two elements of the first dimension correspond to the meridional and zonal components
!       (in that order) of the horizontal spherical harmonic vectors.
!
!ROUTINES USED
!        None.


subroutine vsh_basis(theta, phi)

    use dwm_module, only:m_max=>mmax, n_max=>lmax, nvshfn, vsh_terms
    use vsh_basis_module

    implicit none

    real(4)            :: theta, phi
    real(8)            :: x, y, z, mz, norm_m
    integer(4)         :: n, m, ifn
    integer(4)         :: fn_id

    real(8), parameter :: pi=3.141592653590, dtor=pi/180D0

    !PROCESS INPUT ARGUMENTS

    x = cos((90. - dble(theta)) * dtor)
    y = dsqrt(1-x*x)                        !y = sin(theta)
    z = dble(phi) * dtor

    !CALCULATE Pnm / sin(theta)
    if (m_max .ge. 1) B(1,1) = dsqrt(3D0)
    do m = 2, m_max
      B(m,m) = y * cm(m) * B(m-1,m-1)
    enddo
    do m = 1, m_max
      do n = m+1, n_max
        B(m,n) = anm(m,n) * x * B(m,n-1) - bnm(m,n) * B(m,n-2)
      end do
    end do

    !CALCULATE d(Pnm) / d(theta)
    do m = 1, m_max
      do n = m, n_max
        A(m,n) = n_arr(n) * x * B(m,n) - fnm(m,n) * B(m,n-1)
      end do
    end do
    do n = 1, n_max
      A(0,n) = -e0n(n) * y * B(1,n)
    end do

    !CALCULATE m(Pnm) / sin(theta) AND APPLY SECTORAL NORMALIZATION
    do m = 0, m_max
      if (m .eq. 0) then
        norm_m = 1D0/sqrt(2D0)  !Holmes and Featherstone norm factor adjusted to match Swartztrauber
      else
        norm_m = 0.5D0          !Holmes and Featherstone norm factor adjusted to match Swartztrauber
      end if
      do n = m, n_max
        B(m,n) = B(m,n) * m_arr(m) * cn(n) * norm_m
        A(m,n) = A(m,n) *            cn(n) * norm_m
      end do
    end do

    !CALCULATE VECTOR SPHERICAL HARMONIC FUNCTIONS

    do m = 0, m_max
      mz = dble(m)*z
      cosmz(m) = cos(mz)
      sinmz(m) = sin(mz)
    end do
    do ifn = 0, nvshfn-1
      m = fn_map1(ifn)%M
      n = fn_map1(ifn)%N
      fn_id = fn_map1(ifn)%realimag + 2*fn_map1(ifn)%irrotational
      select case (fn_id)
      case (0)
        vsh_terms(0,ifn) = real( -A(m,n) * cosmz(m) )  !Multiplies real  pt of irrot. coeff. (b)
        vsh_terms(1,ifn) = real( -B(m,n) * sinmz(m) )  !Multiplies real  pt of irrot. coeff. (b)
      case (1)
        vsh_terms(0,ifn) = real(  A(m,n) * sinmz(m) )  !Multiplies imag. pt of irrot. coeff. (b)
        vsh_terms(1,ifn) = real( -B(m,n) * cosmz(m) )  !Multiplies imag. pt of irrot. coeff. (b)
      case (2)
        vsh_terms(0,ifn) = real(  B(m,n) * sinmz(m) )  !Multiplies real  pt of solen. coeff. (c)
        vsh_terms(1,ifn) = real( -A(m,n) * cosmz(m) )  !Multiplies real  pt of solen. coeff. (c)
      case (3)
        vsh_terms(0,ifn) = real(  B(m,n) * cosmz(m) )  !Multiplies imag. pt of solen. coeff. (c)
        vsh_terms(1,ifn) = real(  A(m,n) * sinmz(m) )  !Multiplies imag. pt of solen. coeff. (c)
      end select
    end do

    return

end subroutine vsh_basis



!***************************************************************************************************
!
!JOHN EMMERT 5/3/07
!EVALUATES A BASIS OF QUADRATIC KP SPLINES WITH NODES AT [-10,-8,0,2,5,8,18,20]
!AND SUMS THE FIRST TWO AND LAST TWO BASIS FUNCTIONS. THIS IS EQUIVALENT TO IMPOSING
!A CONSTRAINT OF ZERO SLOPE AT Kp=0 and Kp=8. THE FUNCTION RETURNS THREE TERMS. INPUT
!Kp VALUES GREATER THAN 8 ARE TRUNCATED TO 8.
!
!TRANSLATED TO F90 FROM IDL VERSION, 7/7/07
!
!CALLING SEQUENCE
!   Result = CALL DWM_KPSPL3_CALC( kp, dwm_kpspl3 )
!
!ARGUMENTS
!   kp          Kp index (0-8)
!   dwm_kpspl3    A 3-element real array containing the resulting basis functions.
!
!ROUTINES USED
!    bspline_calc


subroutine dwm_kpspl3_calc(kp0, dwm_kpspl3)

    implicit none

    real(4), intent(in)       :: kp0
    real(4), intent(out)      :: dwm_kpspl3(0:2)

    real(4)                   :: kp
    real(4)                   :: kpspl(0:4)

    external :: bspline_calc

    kp = max(kp0, 0.)
    kp = min(kp,  8.)
    call bspline_calc(8,  kp, (/-10., -8., 0., 2., 5., 8., 18., 20./), kpspl, 2, 0)
    dwm_kpspl3(0) = kpspl(0) + kpspl(1)
    dwm_kpspl3(1) = kpspl(2)
    dwm_kpspl3(2) = kpspl(3) + kpspl(4)

    return

end subroutine dwm_kpspl3_calc

!***************************************************************************************************

!JOHN EMMERT 5/4/07
!COMPUTES A LATITUDE DEPENDENT WEIGHTING FUNCTION THAT GOES TO
!ZERO AT LOW LATITUDES AND ONE AT HIGH LATITUDES. THE TRANSITION
!IS AN EXPONENTIAL S-CURVE, WITH THE TRANSITION LATITUDE DETERMINED
!BY THE MLT/Kp MODEL GENERATED BY dwm_maxgrad1_02.pro, AND A
!TRANSITION WITH SPECIFIED BY THE CALLING ROUTINE.
!
!TRANSLATED TO F90 FROM IDL VERSION, 7/7/07
!
!
!CALLING SEQUENCE
!   Result = DWM_LATWGT2(mlat, mlt, kp, twidth)
!
!ARGUMENTS
!     mlt:       An n-element array of magnetic local times (hours)
!     mlat:      Magnetic latitude (scalar value)
!     kp:        3-hour Kp index
!     twidth:    Latitude transition width
!
!ROUTINES USED
!    None


function dwm_latwgt2(mlat, mlt, kp0, twidth)

    implicit none

    real(4)                   :: dwm_latwgt2
    real(4)                   :: mlat, mlt, kp0, kp, twidth
    real(4)                   :: mltrad, sinmlt, cosmlt, tlat

    real(4), parameter :: coeff(0:5) = (/ 65.7633,  -4.60256,  -3.53915,  &
                                         -1.99971,  -0.752193,  0.972388 /)
    real(4), parameter :: pi=3.141592653590, dtor=pi/180D0


    mltrad = mlt * 15.0 * dtor
    sinmlt = sin(mltrad)
    cosmlt = cos(mltrad)
    kp = max(kp0, 0.)
    kp = min(kp,  8.)
    tlat = coeff(0) + coeff(1)*cosmlt + coeff(2)*sinmlt +   &
           kp*(coeff(3) + coeff(4)*cosmlt + coeff(5)*sinmlt)
    dwm_latwgt2 = 1.0 / ( 1 + exp(-(abs(mlat)-tlat)/twidth) )

    return

end function dwm_latwgt2



!******************************************************************************
!******************************************************************************
!
!BSPLINE
!JOHN EMMERT 3/31/05
!TRANSLATED TO FORTRAN-90 10/4/06. FORTRAN VERSION ONLY ALLOWS SCALAR X ARGUMENT
!CALCULATES A BASIS OF B-SPLINES OF SPECIFIED ORDER
!USES DE BOOR'S RECURSION RELATION
!DeBoor, C. A. (1978), A practical guide to splines, Appl. Math. Sci., 27,
!  Springer, New York.
!
!CALLING SEQUENCE:   CALL BSPLINE_CALC(m, x, nodes, bspline, order, periodic)
!
!ARGUMENTS
!     m:         The number of elements in the 'nodes' argument array (integer).
!     x:         The dependent variable location at which the spline basis is
!                to be evaluated (real).
!     nodes:     Vector of ordered node positions (real array).
!
!                For a non-periodic basis:
!                  Add k-1 (=the value of the order keyword, below) nodes to
!                  either end of the range (domain) of the data. The extra
!                  nodes on the left are the starting points of splines that
!                  end inside the domain.  The extra nodes on the right are the
!                  end points of splines that start inside the domain. The
!                  splines actually used in the basis are those that begin at
!                  the first m-k nodes. The sum of splines is then one (i.e.,
!                  normalized) over the domain.
!
!                For a periodic basis:
!                  The last node is identified with the first, so that there are
!                  m-1 splines in the basis. For example, for a periodic basis
!                  of 6 evenly spaced splines over the interval [0,24], with the
!                  first node at 2, use the following array of 7 nodes:
!                  [2,6,10,14,18,22,26]
!
!     order:     Set this integer argument to the desired spline order (=k-1
!                in De Boor's notation).
!                  Simple bins:           order=0
!                  Linear interpolation:  order=1
!                  Quadratic splines:     order=2
!                  Cubic splines:         order=3 (default)
!                  etc.
!                Note that the number of nodes must be greater than k.
!
!     periodic:  Set this integer argument to 1 if a basis of periodic
!                splines is desired; 0 otherwise.

!RESULT
!     The function returns bspline, a m1-element array, where m1 = m-k in the
!     case of a non-periodic basis and m-1 in the case of a periodic basis. The
!     elements of the array represent each of the basis splines.

!EXAMPLE
!     x = 18.0
!     nodes = (/1.0, 3.5, 6.0, 18.0, 20.5, 23.0, 25.0/)
!     m = size(nodes)
!     order = 3
!     periodic = 1
!     call bspline_calc(m, x, nodes, bspline, order, periodic)
!     print *, bspline
!       0.025355    0.390467    0.584179    0.000000    0.000000    0.000000
!
!ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
!     pershift
!


subroutine bspline_calc(nnode0, x0, node0, bspline, order, periodic)

  implicit none

  integer(4), intent(in) :: nnode0, order, periodic
  real(4), intent(in)    :: x0, node0(0:nnode0-1)
  integer(4)             :: i, j, k, nnode, nspl
  real(4)                :: x, perint(0:1), perspan, dx1, dx2
  real(4), allocatable   :: node(:), node1(:)
  real(4), allocatable   :: bspline0(:)
  real(4), intent(out)   :: bspline(0:nnode0-1-periodic-(1-periodic)*(order+1))
  real(4), external      :: pershift


  !PREPARE WORKING ARRAYS
  k = order + 1
  nnode=nnode0 + periodic*order
  nspl = nnode - k
  x = x0
  allocate(node(0:nnode-1), node1(0:nnode-1))
  allocate(bspline0(0:nnode-2))
  if (periodic .eq. 1) then
    perint = (/ node0(0), node0(nnode0-1) /)
    perspan = perint(1) - perint(0)
    x = pershift(x, perint)
    node = (/ node0, node0(1:order) + perspan /)
    do i=0, nnode0+order-1
      node1(i) = pershift(node(i), perint)
    end do
  else
    node = node0
    node1 = node
  end if

  !COMPUTE SPLINES
  do i = 0, nnode-2
    bspline0(i) = 0.0
    if (node1(i+1) .gt. node1(i)) then
      if ((x .ge. node1(i)) .and. (x .lt. node1(i+1))) bspline0(i) = 1.0
    else
      if ((x .ge. node1(i)) .or.  (x .lt. node1(i+1))) bspline0(i) = 1.0
    end if
  end do
  do j = 2, k
    do i = 0, nnode-j-1
      dx1 = x - node1(i)
      dx2 = node1(i+j) - x
      if (periodic .eq. 1) then
        if (dx1 .lt. 0) dx1 = dx1 + perspan
        if (dx2 .lt. 0) dx2 = dx2 + perspan
      end if
      bspline0(i) =    bspline0(i)   * dx1 / (node(i+j-1) - node(i)) &
                     + bspline0(i+1) * dx2 / (node(i+j)   - node(i+1))
    end do
  end do
  bspline = bspline0(0:nspl-1)
  deallocate(node, node1, bspline0)

  return

end subroutine bspline_calc



!******************************************************************************
!******************************************************************************
!
!PERSHIFT
!JOHN EMMERT   9/12/03
!TRANSLATED TO FORTRAN-90 10/4/06. FORTRAN VERSION ONLY ALLOWS SCALAR INPUTS
!SHIFTS INPUT VALUES INTO A SPECIFIED PERIODIC INTERVAL
!
!CALLING SEQUENCE:   Result = PERSHIFT(x, range)
!
!ARGUMENTS
!      x:        The value to be shifted
!      perint:   2-element vector containing the start and end values
!                of the desired periodic interval.  The periodicity is
!                determined by the span of the range.
!
!ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
!      None

function pershift(x, perint)

  real(4), parameter :: tol=1e-4
  real(4)            :: x, perint(0:1)
  real(4)            :: a, span, offset, offset1, pershift

  pershift = x
  a = perint(0)
  span = perint(1) - perint(0)
  if (span .ne. 0) then
    offset = x-a
    offset1 = mod(offset,span)
    if (abs(offset1) .lt. tol) offset1 = 0
  endif
  pershift = a + offset1
  if ((offset .lt. 0) .and. (offset1 .ne. 0)) pershift = pershift + span

  return

end function pershift



!******************************************************************************
!******************************************************************************
!
!AP_TO_KP
!JOHN EMMERT   7/10/07
!CONVERTS AP VALUES TO KP VALUES, VIA LINEAR INTERPOLATION ON THE LOOKUP
!TABLE
!
!CALLING SEQUENCE:   Result = AP_TO_KP(ap)
!
!ARGUMENTS
!      ap:       The ap index. Values <0 or >400 are truncated
!                to 0 and 400, respectively
!
!ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
!      None


function ap_to_kp(ap0)

  real(4), parameter :: apgrid(0:27) = (/0.,2.,3.,4.,5.,6.,7.,9.,12.,15.,18., &
                                         22.,27.,32.,39.,48.,56.,67.,80.,94., &
                                       111.,132.,154.,179.,207.,236.,300.,400./)
  real(4), parameter :: kpgrid(0:27) = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11., &
                                         12.,13.,14.,15.,16.,17.,18.,19.,20.,21., &
                                         22.,23.,24.,25.,26.,27./) / 3.0
  real(4)            :: ap0, ap, ap_to_kp
  integer(4)         :: i


  ap = ap0
  if (ap .lt. 0) ap=0
  if (ap .gt. 400) ap=400

  i = 1
  do while (ap .gt. apgrid(i))
    i = i + 1
  end do
  if (ap .eq. apgrid(i)) then
    ap_to_kp = kpgrid(i)
  else
    ap_to_kp = kpgrid(i-1) + (ap - apgrid(i-1)) / (3.0 * (apgrid(i) - apgrid(i-1)))
  end if

  return

end function ap_to_kp

!******************************************************************************
!******************************************************************************
!
!DWM_ALTWGT
!JOHN EMMERT   7/10/07
!COMPUTES AN EXPONENTIAL STEP FUNCTION IN HEIGHT TO BE APPLIED TO DWM WINDS (WHICH
!HAVE NO HEIGHT DEPENDENCE). THE FUNCTION CUTS OFF DISTURBANCE WINDS BELOW 125 KM.
!
!CALLING SEQUENCE:   Result = DWM_ALTWGT(alt)
!
!ARGUMENTS
!      alt:      Height in km
!
!ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
!      None

function dwm_altwgt(alt)

  real(4), parameter :: talt=125.0, twidth=5.0
  real(4)            :: alt, dwm_altwgt

  dwm_altwgt = 1.0 / (1 + exp(-(alt - talt)/twidth))

  return

end function dwm_altwgt

!*********************************************************************

subroutine gd2qd(glat,glon,qdlat,qdlon,f1e,f1n,f2e,f2n)

    implicit none

    real(4), intent(in)      :: glat,glon
    real(4), intent(out)     :: qdlat, qdlon
    real(4), intent(out)     :: f1e,f1n,f2e,f2n
    real(4)                  :: f1(2),f2(2)
    integer(4)               :: ist
    real(4), parameter       :: alt = 250.0
    real(4)                  :: hr
    
    call apex(glat,glon,alt,hr,qdlon,qdlat,f1,f2,ist)
                    
    if (ist .gt. 0) stop

    f1e = f1(1)
    f1n = f1(2)
    f2e = f2(1)
    f2n = f2(2)

    return

end subroutine gd2qd
