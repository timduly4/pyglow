!
!  Test driver for HWM07 subroutines
!  The output of the program is given at the end of the file
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
!*******************************************************************

program checkhwm07

  implicit none
  INTEGER      :: IYD
  REAL(4)      :: SEC,ALT,GLAT,GLON,STL,F107A,F107,AP(2)
  REAL(4)      :: W(2), QW(2), DW(2)
  real(4)      :: mlt, mlat, kp, mmpwind, mzpwind
  real(4)      :: ut, apqt(2)
  integer      :: day
  real(4)      :: pershift
  integer      :: ialt,istl,ilat,ilon,iaptemp
  integer      :: imlat,imlt,ikp

! HEIGHT PROFILE
  day = 150
  iyd = 95000 + day
  ut = 12.0
  sec = ut * 3600.0
  glat = -45.0
  glon = -85.0
  stl = pershift(ut + glon/15.0, (/0.0, 24.0/) )
  ap(2) = 80.0
  apqt(2) = -1.0
  print *, 'HEIGHT PROFILE'
  print '(a5,i3, a5,f4.1, a7,f5.1, a7,f6.1, a6,f4.1, a5,f5.1)', &
            'DAY=',day, ', UT=',ut, ', GLAT=',glat,  &
            ', GLON=',glon, ', STL=',stl, ', ap=',ap(2)
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'ALT', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do ialt = 0, 400, 25
    alt = float(ialt)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f6.0,3(f12.3,f10.3))', alt, qw, dw, w
  end do
  print *
  print *

! LATITUDE PROFILE
  day = 305
  iyd = 95000 + day
  ut = 18.0
  sec = ut * 3600.0
  alt = 250.0
  glon = 30.0
  stl = pershift(ut + glon/15.0, (/0.0, 24.0/) )
  ap(2) = 48.0
  apqt(2) = -1.0
  print *, 'LATITUDE PROFILE'
  print '(a5,i3, a5,f4.1, a6,f5.1, a7,f6.1, a6,f4.1, a5,f5.1)', &
            'DAY=',day, ', UT=',ut, ', ALT=',alt,  &
            ', GLON=',glon, ', STL=',stl, ', ap=',ap(2)
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'GLAT', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do ilat = -90, 90, 10
    glat = float(ilat)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f6.1,3(f12.3,f10.3))', glat, qw, dw, w
  end do
  print *
  print *

! LOCAL TIME PROFILE
  day = 75
  iyd = 95000 + day
  alt = 125.0
  glat = 45.0
  glon = -70.0
  ap(2) = 30.0
  apqt(2) = -1.0
  print *, 'LOCAL TIME PROFILE'
  print '(a5,i3, a6,f5.1, a7,f5.1, a7,f6.1, a5,f5.1)', &
            'DAY=',day, ', ALT=',alt, ', GLAT=',glat,  &
            ', GLON=',glon, ', ap=',ap(2)
  print '(5x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a5,3(a12,a10))', 'STL', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do istl = 0,16
    stl = 1.5 * float(istl)
    sec = (stl - glon/15.0) * 3600.0
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f5.1,3(f12.3,f10.3))', stl, qw, dw, w
  end do
  print *
  print *

! LONGITUDE PROFILE
  day = 330
  iyd = 95000 + day
  ut = 6.0
  sec = ut * 3600.0
  alt = 40.0
  glat = -5.0
  ap(2) = 4.0
  apqt(2) = -1.0
  print *, 'LONGITUDE PROFILE'
  print '(a5,i3, a5,f4.1, a6,f5.1, a7,f5.1, a7,f6.1, a5,f5.1)', &
            'DAY=',day, ', UT=',ut, ', ALT=',alt, ', GLAT=',glat,  &
            ', GLON=',glon, ', ap=',ap(2)
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'GLON', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do ilon = -180, 180, 20
    glon = float(ilon)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f6.0,3(f12.3,f10.3))', glon, qw, dw, w
  end do
  print *
  print *

! DAY OF YEAR PROFILE
  ut = 21.0
  sec = ut * 3600.0
  alt = 200.0
  glat = -65.0
  glon = -135.0
  stl = pershift(ut + glon/15.0, (/0.0, 24.0/) )
  ap(2) = 15.0
  apqt(2) = -1.0
  print *, 'DAY OF YEAR PROFILE'
  print '(a4,f4.1, a6,f5.1, a7,f5.1, a7,f6.1, a6,f4.1, a5,f5.1)', &
            'UT=',ut, ', ALT=',alt, ', GLAT=',glat,  &
            ', GLON=',glon, ', STL=',stl, ', ap=',ap(2)
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'DAY', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do day = 0, 360, 20
    iyd = 95000 + day
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(i6,3(f12.3,f10.3))', day, qw, dw, w
  end do
  print *
  print *

! AP PROFILE
  day = 280
  iyd = 95000 + day
  ut = 21.0
  sec = ut * 3600.0
  alt = 350.0
  glat = 38.0
  glon = 125.0
  stl = pershift(ut + glon/15.0, (/0.0, 24.0/) )
  ap(2) = 48.0
  apqt(2) = -1.0
  print *, 'MAGNETIC ACTIVITY PROFILE'
  print '(a5,i3, a5,f4.1, a6,f5.1, a7,f5.1, a7,f6.1, a6,f4.1)', &
            'DAY=',day, ', UT=',ut, ', ALT=',alt,  &
            ', GLAT=',glat, ', GLON=',glon, ', STL=',stl
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'ap', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do iaptemp = 0, 260, 20
    ap(2) = float(iaptemp)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f6.1,3(f12.3,f10.3))', ap(2), qw, dw, w
  end do
  print *
  print *

! DWM: MLAT PROFILE
  kp = 6.0
  mlt = 3.0
  print *, 'DWM: MAGNETIC LATITUDE PROFILE'
  print '(a5,f4.1, a5,f3.1)', 'MLT=',mlt, ', Kp=',kp
  print '(a6,a12,a10)', 'MLAT', 'MAG MER', 'MAG ZON'
  do imlat = -90, 90, 10
    mlat = float(imlat)
    call dwm07b(mlt, mlat, kp, mmpwind, mzpwind)
    print '(f6.1,f12.3,f10.3)', mlat, mmpwind, mzpwind
  end do
  print *
  print *

! DWM: MLT PROFILE
  kp = 6.0
  mlat = 45.0
  print *, 'DWM: MAGNETIC LOCAL TIME PROFILE'
  print '(a6,f5.1, a5,f3.1)', 'MLAT=',mlat, ', Kp=',kp
  print '(a6,a12,a10)', 'MLT', 'MAG MER', 'MAG ZON'
  do imlt = 0, 16
    mlt = float(imlt)*1.5
    call dwm07b(mlt, mlat, kp, mmpwind, mzpwind)
    print '(f6.1,f12.3,f10.3)', mlt, mmpwind, mzpwind
  end do
  print *
  print *

! DWM: KP PROFILE
  mlat = -50.0
  mlt = 3.0
  print *, 'DWM: Kp PROFILE'
  print '(a6,f5.1, a6,f4.1)', 'MLAT=',mlat, ', MLT=',mlt
  print '(a6,a12,a10)', 'Kp', 'MAG MER', 'MAG ZON'
  do ikp = 0, 18
    kp = float(ikp)*0.5
    call dwm07b(mlt, mlat, kp, mmpwind, mzpwind)
    print '(f6.1,f12.3,f10.3)', kp, mmpwind, mzpwind
  end do
  print *
  print *

end program checkhwm07


!******************************************************************************
!  TEST OUTPUT
!******************************************************************************


! HEIGHT PROFILE
! DAY=150, UT=12.0, GLAT=-45.0, GLON= -85.0, STL= 6.3, ap= 80.0
!                       QUIET             DISTURBED                 TOTAL
!   ALT         MER       ZON         MER       ZON         MER       ZON
!    0.       0.484     7.859       0.000     0.000       0.484     7.859
!   25.       2.502    24.992       0.000     0.000       2.502    24.992
!   50.      -7.478    98.432       0.000     0.000      -7.478    98.432
!   75.      15.985    49.706       0.002    -0.001      15.987    49.705
!  100.     -38.402    45.775       0.304    -0.126     -38.098    45.648
!  125.      36.199   -27.612      22.679    -9.436      58.878   -37.048
!  150.      -3.955   -32.830      45.054   -18.745      41.099   -51.575
!  175.     -23.771   -38.211      45.355   -18.871      21.584   -57.081
!  200.     -11.783   -36.996      45.357   -18.872      33.575   -55.868
!  225.      -0.942   -38.801      45.358   -18.872      44.415   -57.673
!  250.       4.802   -42.668      45.358   -18.872      50.160   -61.540
!  275.       8.318   -45.691      45.358   -18.872      53.675   -64.563
!  300.      10.635   -47.684      45.358   -18.872      55.993   -66.556
!  325.      12.163   -48.998      45.358   -18.872      57.520   -67.870
!  350.      13.170   -49.864      45.358   -18.872      58.528   -68.736
!  375.      13.834   -50.435      45.358   -18.872      59.192   -69.307
!  400.      14.272   -50.812      45.358   -18.872      59.629   -69.684
! 
! 
! LATITUDE PROFILE
! DAY=305, UT=18.0, ALT=250.0, GLON=  30.0, STL=20.0, ap= 48.0
!                       QUIET             DISTURBED                 TOTAL
!  GLAT         MER       ZON         MER       ZON         MER       ZON
! -90.0     -46.412   110.719     -86.605    45.868    -133.095   156.595
! -80.0     -26.264    68.190    -134.180   -11.922    -160.444    56.268
! -70.0       9.717    28.536    -137.046  -109.312    -127.328   -80.777
! -60.0      32.894    11.525     -58.228  -156.647     -25.334  -145.122
! -50.0      32.972    20.006     -15.470  -112.360      17.502   -92.354
! -40.0      25.222    44.460     -10.996   -52.103      14.226    -7.643
! -30.0      26.510    67.704       7.172   -36.536      33.682    31.168
! -20.0      34.481    75.179      14.921   -34.819      49.402    40.360
! -10.0      35.540    68.603      11.329   -27.280      46.869    41.322
!   0.0      26.231    64.273       6.125   -17.274      32.356    46.999
!  10.0      18.178    72.452       1.009   -13.377      19.187    59.075
!  20.0      22.910    84.147      -4.819   -20.219      18.091    63.928
!  30.0      37.805    84.437     -12.107   -34.870      25.698    49.567
!  40.0      47.879    72.691     -15.141   -42.692      32.737    29.998
!  50.0      38.036    60.953      -6.580   -54.902      31.456     6.051
!  60.0       4.256    57.496      -7.254  -201.700      -2.998  -144.204
!  70.0     -42.170    63.602     -36.408   -52.818     -78.579    10.783
!  80.0     -79.970    83.380    -103.037    81.185    -183.006   164.565
!  90.0     -92.092   120.221     -75.803    89.986    -167.896   210.207
! 
! 
! LOCAL TIME PROFILE
! DAY= 75, ALT=125.0, GLAT= 45.0, GLON= -70.0, ap= 30.0
!                      QUIET             DISTURBED                 TOTAL
!  STL         MER       ZON         MER       ZON         MER       ZON
!  0.0     -25.523   -25.139      -8.371   -20.333     -33.894   -45.472
!  1.5     -29.545   -53.409      -7.050     2.629     -36.654   -50.787
!  3.0     -21.723   -65.269      -5.992    20.035     -27.715   -45.234
!  4.5      -8.261   -47.142      -8.483    21.325     -16.744   -25.817
!  6.0      -0.534   -15.152     -13.947     8.969     -14.481    -6.182
!  7.5       1.132     2.634     -16.888    -1.507     -15.682     1.143
!  9.0       8.180    -2.976     -16.086    -4.841      -7.906    -7.818
! 10.5      27.106   -14.648     -13.948    -4.077      13.158   -18.725
! 12.0      46.642   -11.904     -11.530    -1.994      35.112   -13.898
! 13.5      46.855     3.871      -8.293     0.136      38.562     4.007
! 15.0      21.394    12.923      -4.773    -0.880      16.621    12.043
! 16.5     -13.202     3.319      -3.926   -13.039     -17.128    -9.720
! 18.0     -33.181   -13.810      -7.930   -40.623     -41.111   -54.433
! 19.5     -31.037   -19.100     -12.826   -64.625     -43.863   -83.725
! 21.0     -20.447   -10.681     -12.427   -60.116     -32.874   -70.797
! 22.5     -18.238    -7.532      -9.901   -41.765     -28.138   -49.298
! 24.0     -25.523   -25.139      -8.371   -20.333     -33.894   -45.472
! 
! 
! LONGITUDE PROFILE
! DAY=330, UT= 6.0, ALT= 40.0, GLAT= -5.0, GLON= -70.0, ap=  4.0
!                       QUIET             DISTURBED                 TOTAL
!  GLON         MER       ZON         MER       ZON         MER       ZON
! -180.       0.125   -15.666       0.000     0.000       0.125   -15.666
! -160.       0.669   -18.202       0.000     0.000       0.669   -18.202
! -140.       0.871   -20.365       0.000     0.000       0.871   -20.365
! -120.       0.679   -21.247       0.000     0.000       0.679   -21.247
! -100.       0.227   -20.376       0.000     0.000       0.227   -20.376
!  -80.      -0.229   -17.950       0.000     0.000      -0.229   -17.950
!  -60.      -0.440   -14.761       0.000     0.000      -0.440   -14.761
!  -40.      -0.284   -11.873       0.000     0.000      -0.284   -11.873
!  -20.       0.175   -10.172       0.000     0.000       0.175   -10.172
!    0.       0.716   -10.009       0.000     0.000       0.716   -10.009
!   20.       1.063   -11.088       0.000     0.000       1.063   -11.088
!   40.       1.020   -12.651       0.000     0.000       1.020   -12.651
!   60.       0.564   -13.864       0.000     0.000       0.564   -13.864
!   80.      -0.137   -14.216       0.000     0.000      -0.137   -14.216
!  100.      -0.798   -13.754       0.000     0.000      -0.798   -13.754
!  120.      -1.146   -13.036       0.000     0.000      -1.146   -13.036
!  140.      -1.042   -12.827       0.000     0.000      -1.042   -12.827
!  160.      -0.543   -13.690       0.000     0.000      -0.543   -13.690
!  180.       0.125   -15.666       0.000     0.000       0.125   -15.666
! 
! 
! DAY OF YEAR PROFILE
! UT=21.0, ALT=200.0, GLAT=-65.0, GLON=-135.0, STL=12.0, ap= 15.0
!                       QUIET             DISTURBED                 TOTAL
!   DAY         MER       ZON         MER       ZON         MER       ZON
!     0      -1.953   -50.220       4.688     9.403       2.735   -40.817
!    20      -3.100   -45.982       4.689     9.400       1.589   -36.582
!    40      -8.144   -44.454       4.691     9.392      -3.453   -35.062
!    60     -16.767   -43.867       4.692     9.386     -12.075   -34.481
!    80     -28.168   -41.056       4.689     9.401     -23.479   -31.655
!   100     -41.057   -33.308       4.673     9.468     -36.384   -23.840
!   120     -53.790   -20.060       4.645     9.595     -49.145   -10.464
!   140     -64.620    -3.635       4.617     9.735     -60.003     6.100
!   160     -72.017    11.376       4.601     9.823     -67.416    21.199
!   180     -74.968    19.873       4.601     9.823     -70.367    29.696
!   200     -73.161    18.359       4.617     9.735     -68.544    28.094
!   220     -67.005     6.490       4.645     9.595     -62.360    16.086
!   240     -57.496   -12.646       4.673     9.468     -52.822    -3.178
!   260     -45.970   -33.666       4.689     9.401     -41.281   -24.265
!   280     -33.850   -50.949       4.692     9.386     -29.158   -41.563
!   300     -22.446   -60.734       4.691     9.392     -17.755   -51.342
!   320     -12.849   -62.349       4.689     9.400      -8.160   -52.949
!   340      -5.917   -58.075       4.688     9.403      -1.228   -48.672
!   360      -2.295   -51.752       4.688     9.403       2.393   -42.349
! 
! 
! MAGNETIC ACTIVITY PROFILE
! DAY=280, UT=21.0, ALT=350.0, GLAT= 38.0, GLON= 125.0, STL= 5.3
!                       QUIET             DISTURBED                 TOTAL
!    ap         MER       ZON         MER       ZON         MER       ZON
!   0.0      21.580   -43.920      -1.446    -5.351      20.134   -49.272
!  20.0      21.580   -43.920      -9.171    -2.526      12.409   -46.447
!  40.0      21.580   -43.920     -20.753   -16.908       0.827   -60.829
!  60.0      21.580   -43.920     -29.600   -30.261      -8.020   -74.181
!  80.0      21.580   -43.920     -34.510   -38.143     -12.931   -82.063
! 100.0      21.580   -43.920     -37.454   -43.207     -15.874   -87.128
! 120.0      21.580   -43.920     -39.084   -46.293     -17.505   -90.213
! 140.0      21.580   -43.920     -40.001   -48.296     -18.422   -92.216
! 160.0      21.580   -43.920     -40.448   -49.598     -18.868   -93.518
! 180.0      21.580   -43.920     -40.516   -50.297     -18.936   -94.217
! 200.0      21.580   -43.920     -40.312   -50.531     -18.732   -94.452
! 220.0      21.580   -43.920     -40.183   -50.525     -18.603   -94.445
! 240.0      21.580   -43.920     -40.183   -50.525     -18.603   -94.445
! 260.0      21.580   -43.920     -40.183   -50.525     -18.603   -94.445
! 
! 
! DWM: MAGNETIC LATITUDE PROFILE
! MLT= 3.0, Kp=6.0
!  MLAT     MAG MER   MAG ZON
! -90.0     157.613  -158.954
! -80.0     158.138  -141.865
! -70.0      56.293   -41.271
! -60.0      60.878    69.305
! -50.0     121.506    -5.004
! -40.0      40.882   -15.775
! -30.0      21.850   -25.517
! -20.0      11.640   -34.856
! -10.0       5.533   -40.228
!   0.0       5.327   -45.354
!  10.0       8.204   -50.906
!  20.0       9.422   -55.711
!  30.0       3.494   -56.017
!  40.0      -7.319   -49.619
!  50.0      -7.873    18.973
!  60.0     -33.112    80.795
!  70.0     -20.228    -8.189
!  80.0     -55.275   -64.783
!  90.0    -137.149  -170.526
! 
! 
! DWM: MAGNETIC LOCAL TIME PROFILE
! MLAT= 45.0, Kp=6.0
!   MLT     MAG MER   MAG ZON
!   0.0     -22.332  -129.270
!   1.5     -28.861   -78.026
!   3.0      -6.246   -33.077
!   4.5      -5.592   -33.656
!   6.0     -32.813   -39.796
!   7.5     -49.672   -34.709
!   9.0     -48.676   -26.754
!  10.5     -41.507   -20.200
!  12.0     -35.088   -15.157
!  13.5     -28.140   -12.072
!  15.0     -18.344   -13.870
!  16.5      -9.103   -27.605
!  18.0     -12.015   -57.829
!  19.5     -27.512   -86.677
!  21.0     -24.814  -103.292
!  22.5     -11.058  -127.315
!  24.0     -22.332  -129.270
! 
! 
! DWM: Kp PROFILE
! MLAT=-50.0, MLT= 3.0
!    Kp     MAG MER   MAG ZON
!   0.0      -8.066     5.787
!   0.5      -6.955     5.282
!   1.0      -3.365     3.849
!   1.5       3.202     1.647
!   2.0      13.413    -1.112
!   2.5      26.175    -3.461
!   3.0      39.805    -4.712
!   3.5      54.121    -5.081
!   4.0      68.787    -4.864
!   4.5      83.302    -4.428
!   5.0      97.009    -4.194
!   5.5     109.699    -4.429
!   6.0     121.506    -5.004
!   6.5     132.195    -5.750
!   7.0     141.548    -6.483
!   7.5     149.386    -7.019
!   8.0     155.577    -7.193
!   8.5     155.577    -7.193
!   9.0     155.577    -7.193
!
