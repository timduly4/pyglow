c irisub.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------        
C Includes subroutines IRI_SUB and IRI_WEB to compute IRI parameters 
C for specified location, date, time, and altitude range and subroutine 
C IRI_WEB to computes IRI parameters for specified location, date, time 
C and variable range; variable can be altitude, latitude, longitude, 
C year, month, day of month, day of year, or hour (UT or LT). 
C IRI_WEB requires IRI_SUB. Both subroutines require linking with the 
c following library files IRIFUN.FOR, IRITEC.FOR, IRIDREG.FOR, 
c CIRA.FOR, IGRF.FOR
c-----------------------------------------------------------------------        
c Programs using IRISUB need to include (see IRITEST):
c 
c       COMMON/const2/icalls,nmono,iyearo,idaynro,rzino,igino,ut0
c
c       icalls=0
c       nmono=-1
c       iyearo=-1
c       idaynro=-1
c       rzino=-1
c       igino=-1
c       ut0=-1
c        
c       do 6249 i=1,50
c 6249  oar(i,1)=-1.0
c-----------------------------------------------------------------------        
c Required i/o units:  
c  KONSOL= 6 IRISUB: Program messages (used when jf(12)=.true. -> konsol)
c  IUCCIR=10 IRISUB: CCIR and URSI coefficients (CCIR%%.ASC, %%=month+10)
c  KONSOL=11 IRISUB: Program messages (used when jf(12)=.false. -> MESSAGES.TXT)
c     KONSOL=6/11 is also used in IRIFUN and IGRF. COMMON/iounit/konsol 
c     is used to pass the value of KONSOL. If KONSOL=1 messages are turned off.
c  UNIT=12 IRIFUN/TCON:  Solar/ionospheric indices IG12, R12 (IG_RZ.DAT) 
c  UNIT=13 IRIFUN/APF..: Magnetic indices and F10.7 (APF107.DAT 
c  UNIT=14 IGRF/GETSHC:  IGRF coeff. (DGRF%%%%.DAT or IGRF%%%%.DAT, %%%%=year)
c-----------------------------------------------------------------------        
C CHANGES FROM  IRIS11.FOR  TO   IRIS12.FOR:
C    - CIRA-1986 INSTEAD OF CIRA-1972 FOR NEUTRAL TEMPERATURE
C    - 10/30/91 VNER FOR NIGHTTIME LAY-VERSION:  ABS(..)
C    - 10/30/91 XNE(..) IN CASE OF LAY-VERSION
C    - 10/30/91 CHANGE SSIN=F/T TO IIQU=0,1,2
C    - 10/30/91 Te > Ti > Tn ENFORCED IN FINAL PROFILE
C    - 10/30/91 SUB ALL NAMES WITH 6 OR MORE CHARACTERS
C    - 10/31/91 CORRECTED HF1 IN HST SEARCH:  NE(HF1)>NME
C    - 11/14/91 C1=0 IF NO F1-REGION
C    - 11/14/91 CORRECTED HHMIN AND HZ FOR LIN. APP.
C    -  1/28/92 RZ12=0 included
C    -  1/29/92 NEQV instead of NE between URSIF2 and URSIFO
C    -  5/ 1/92 CCIR and URSI input as in IRID12
C    -  9/ 2/94 Decimal month (ZMONTH) for IONCOM
C    -  9/ 2/94 Replace B0POL with B0_TAB; better annually
C    -  1/ 4/95 DY for h>hmF2
C    -  2/ 2/95 IG for foF2, topside; RZ for hmF2, B0_TAB, foF1, NmD
C    -  2/ 2/95 winter no longer exclusive for F1 occurrrence
C    -  2/ 2/95 RZ and IG incl as DATA statement; smooth annual var.
C CHANGES FROM  IRIS12.FOR  TO   IRIS13.FOR:
C    - 10/26/95 incl year as input and corrected MODA; nrm for zmonth
C    - 10/26/95 use TCON and month-month interpolation in foF2, hmF2
C    - 10/26/95 TCON only if date changes
C    - 11/25/95 take out logicals TOPSI, BOTTO, and BELOWE
C    - 12/ 1/95 UT_LT for (date-)correct UT<->LT conversion
C    - 12/22/95 Change ZETA cov term to cov < 180; use cov inst covsat
C    -  2/23/96 take covmax(R<150) for topside; lyear,.. for lt
C    -  3/26/96 topside: 94.5/BETA inst 94.45/..; cov -> covsat(<=188)
C    -  5/01/96 No longer DY for h>hmF2 (because of discontinuity)
C    - 12/01/96 IRIV13: HOUR for IVAR=1 (height)
C    -  4/25/97 D-region: XKK le 10 with D1 calc accordingly.
C    -  1/12/97 DS model for lower ion compoistion DY model
C    -  5/19/98 seamon=zmonth if lati>0; zmonth= ...(1.0*iday)/..
C    -  5/19/98 DY ion composition model below 300 km now DS model
C    -  5/19/98 DS model includes N+, Cl down to 75 km HNIA changed
C    -  5/28/98 User input for Rz12, foF1/NmF1, hmF1, foE/NmE, hmE
C    -  9/ 2/98 1 instead of 0 in MODA after UT_LT call
C    -  4/30/99 constants moved from DATA statement into program
C    -  4/30/99 changed konsol-unit to 13 (12 is for IG_RZ).
C    -  5/29/99 the limit for IG comp. from Rz12-input is 174 not 274
C    - 11/08/99 jf(18)=t simple UT to LT conversion, otherwise UT_LT
C    - 11/09/99 added COMMON/const1/humr,dumr also for CIRA86
C CHANGES FROM  IRIS13.FOR  TO   IRISUB.FOR:
c-----------------------------------------------------------------------        
C-Version-MM/DD/YY-Description (person reporting correction)
C 2000.01 05/09/00 B0_98 replaces B0_TAB and B1: 1.9/day to 2.6/night
C 2000.02 06/11/00 including new F1 and indermediate region
C 2000.03 10/15/00 include Scherliess-Fejer drift model
C 2000.04 10/29/00 include special option for D region models
C 2000.05 12/07/00 change name IRIS13 to IRISUB
C 2000.06 12/14/00 jf(30),outf(20,100),oarr(50)
C 2000.07 03/17/01 include Truhlik-Triskova Te model and IGRF
C 2000.08 05/07/01 include Fuller-Rowell-Condrescu storm model 
C 2000.09 07/09/01 LATI instead of LAT1 in F00 call -------- M. Torkar
C 2000.10 07/09/01 sdte instead of dte in ELTEIK call --- P. Wilkinson
C 2000.11 09/18/01 correct computation of foF2 for Rz12 user input
C 2000.12 09/19/01 Call APF only if different date and time -- P. Webb
c 2000.13 10/28/02 replace TAB/6 blanks, enforce 72/line -- D. Simpson
C 2000.14 11/08/02 change unit for message file to 11 (13 is Kp)
C 2000.15 01/27/03 change F1_prob output; Te-IK for fix h and ELTE(h)
C 2000.16 02/04/03 along<0 -> along=along+360; F1 occ for hmf1&foF1
C 2000.17 02/05/03 zyear =12.97 (Dec 31); idayy=#days per year
C 2000.18 02/06/03 jf(27) for IG12 user input; all F1 prob in oar
C 2000.19 07/14/04 covsat<188 instead of covsat=<f(IG)<188
C 2000.19 02/09/05 declare INVDIP as real ------------------ F. Morgan
C 2000.20 11/09/05 replace B0B1 with BCOEF --------------- T. Gulyaeva
C 2005.01 11/09/05 new topside ion composition; F107D from file 
C 2005.02 11/14/05 jf(18)=T: dip,mlat IGRF10 (igrf_dip igrf.for); F:POGO-75;
C 2005.03 11/15/05 sunrise/sunset/night for D,E,F1,F2; UT_LT removed
C 2005.04 05/06/06 FIRI D-region option not tied to peak
C 2005.04 05/06/06 Spread-F included, NeQuick included
C 2005.05 01/15/07 NeQuick uses CCIR-M3000F2 even if user-hmF2  
C 2007.00 05/18/07 Release of IRI-2007
C 2007.01 01/23/08 ryear = .. (daynr-1.0)/idayy ---------- R. Scharroo
C 2007.02 10/31/08 outf(100) -> outf(500), numhei=numstp=500 
C 2007.03 02/12/09 Jf(24)=.false.-> outf(1,60-140km)=FIRI- M. Friedrich
C 2007.04 03/14/09 SOCO(70->80;500->300km) --------------- R. Davidson
C 2007.05 03/26/09 call for APF_ONLY includes F107M
C 2007.09 08/17/09 STROM off if input; fof2in, fof1in,foein corr
C 2007.10 02/03/10 F10.7D = F10.7M = COV if EOF
C 2007.11 04/19/10 Corrections in irifun.for, cira.for 
C 2007.12 11/23/10 FNIGHT computed twice at 8334 --------- C. Vasly 
C
C 2012.00 10/05/11 IRI-2012: bottomside B0 B1 model (SHAMDB0D, SHAB1D),
C 2012.00 10/05/11  bottomside Ni model (iriflip.for), auroral foE
C 2012.00 10/05/11  storm model (storme_ap), Te with PF10.7 (elteik),
C 2012.00 10/05/11  oval kp model (auroral_boundary),IGRF-11(igrf.for), 
C 2012.00 10/05/11  NRLMSIS00 (cira.for), CGM coordinates, F10.7 daily
C 2012.00 10/05/11  81-day 365-day indices (apf107.dat), ap->kp (ckp),
C 2012.00 10/05/11  array size change jf(50) outf(20,1000), oarr(100).
C 2012.01 11/01/11 delete TEDER from EXTERNAL; GTD7 call 0 to 0.0
C 2012.01 12/12/11 put FMODIP in EXTERNAL; cgn_lon -> cgm_lon
C 2012.01 01/24/12 Change FLAT to LATI in SHAB1D call [D. Altadill]
C 2012.01 08/09/12 add jf(36)=t/f foF2 for hmF2 wout/with storm
C 2012.01 08/09/12 replace foF2_storm with foF2 for topside (NeQ, corr)
C 2012.01 08/09/12 call stormE_ap only if ap available
C 2012.01 08/09/12 If ap not available then auroral boundary for Kp=3
C 2012.02 12/17/12 Add magnetic declination as oarr(84) output
C 2012.03 02/13/13 Move B1 before B0 for Gulyaeva-1987
C 2012.03 02/20/13 Use foot-point for CGM to be closer to AACGM
C 2012.03 02/20/13 DAT(11,*) is UT time of MLT=0
C
C*****************************************************************
C********* INTERNATIONAL REFERENCE IONOSPHERE (IRI). *************
C*****************************************************************
C**************** ALL-IN-ONE SUBROUTINE  *************************
C*****************************************************************
C
C
       SUBROUTINE IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
C-----------------------------------------------------------------
C
C INPUT:  JF(1:50)      true/false switches for several options
C         JMAG          =0 geographic   = 1 geomagnetic coordinates
C         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C         IYYYY         Year as YYYY, e.g. 1985
C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL 
C                          HOURS
C         HEIBEG,       HEIGHT RANGE IN KM; maximal 100 heights, i.e.
C          HEIEND,HEISTP        int((heiend-heibeg)/heistp)+1.le.100
C
C    JF switches to turn off/on (.true./.false.) several options
C
C    i       .true.                  .false.          standard version
C    -----------------------------------------------------------------
C    1    Ne computed            Ne not computed                     t
C    2    Te, Ti computed        Te, Ti not computed                 t
C    3    Ne & Ni computed       Ni not computed                     t
C    4    B0,B1 - Bil-2000       B0,B1 - other models jf(31)     false
C    5    foF2 - CCIR            foF2 - URSI                     false
C    6    Ni - DS-1995 & DY-1985 Ni - RBV-2010 & TTS-2003        false
C    7    Ne - Tops: f10.7<188   f10.7 unlimited                     t            
C    8    foF2 from model        foF2 or NmF2 - user input           t
C    9    hmF2 from model        hmF2 or M3000F2 - user input        t
C   10    Te - Standard          Te - Using Te/Ne correlation        t
C   11    Ne - Standard Profile  Ne - Lay-function formalism         t
C   12    Messages to unit 6     to messages.text on unit 11         t
C   13    foF1 from model        foF1 or NmF1 - user input           t
C   14    hmF1 from model        hmF1 - user input (only Lay version)t
C   15    foE  from model        foE or NmE - user input             t
C   16    hmE  from model        hmE - user input                    t
C   17    Rz12 from file         Rz12 - user input                   t
C   18    IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
C   19    F1 probability model   critical solar zenith angle (old)   t
C   20    standard F1            standard F1 plus L condition        t
C   21    ion drift computed     ion drift not computed          false
C   22    ion densities in %     ion densities in m-3                t
C   23    Te_tops (Bil-1985)     Te_topside (TBT-2012)           false
C   24    D-region: IRI-1990     FT-2001 and all opts in OUTF(14,*)  t
C   25    F107D from APF107.DAT  F107D user input (oarr(41))         t
C   26    foF2 storm model       no storm updating                   t
C   27    IG12 from file         IG12 - user                         t
C   28    spread-F probability 	 not computed                    false
C   29    IRI01-topside          new options as def. by JF(30)   false
C   30    IRI01-topside corr.    NeQuick topside model   	     false 
C (29,30) = (t,t) IRIold, (f,t) IRIcor, (f,f) NeQuick, (t,f) Gulyaeva
C   31    B0,B1 ABT-2009	     B0 Gulyaeva-1987 h0.5               t   
C (4,31) = (t,t) Bil-00, (f,t) ABT-09, (f,f) Gul-87, (t,f) not used
C   32    F10.7_81 from file     F10.7_81 - user input (oarr(46))    t
C   33    Auroral boundary model on/off  true/false	             false
C   34    Messages on            Messages off                        t
C   35    foE storm model        no foE storm updating           false
C   36    hmF2 w/out foF2_storm  with foF2-storm                     t
C   37    topside w/out foF2-storm  with foF2-storm                  t
C   38    turn WRITEs off in IRIFLIP   turn WRITEs on                t
C   ..    ....
C   50    ....
C   ------------------------------------------------------------------
C
C  Depending on the jf() settings additional INPUT parameters may 
c  be required:
C
C       Setting              INPUT parameter
C    -----------------------------------------------------------------
C    jf(8)  =.false.     OARR(1)=user input for foF2/MHz or NmF2/m-3
C    jf(9)  =.false.     OARR(2)=user input for hmF2/km or M(3000)F2
C    jf(10 )=.false.     OARR(15),OARR(16)=user input for Ne(300km),
C       Ne(400km)/m-3. Use OARR()=-1 if one of these values is not 
C       available. If jf(23)=.false. then Ne(300km), Ne(550km)/m-3.
C    jf(13) =.false.     OARR(3)=user input for foF1/MHz or NmF1/m-3 
C    jf(14) =.false.     OARR(4)=user input for hmF1/km
C    jf(15) =.false.     OARR(5)=user input for foE/MHz or NmE/m-3 
C    jf(16) =.false.     OARR(6)=user input for hmE/km
C    jf(17) =.flase.     OARR(33)=user input for Rz12
C    jf(21) =.true.      OARR(41)=user input for daily F10.7 index
C    jf(23) =.false.     OARR(41)=user input for daily F10.7 index
C    jf(24) =.false.     OARR(41)=user input for daily F10.7 index
C          optional for jf(21:24); default is F10.7D=COV
C    jf(25) =.false.     OARR(41)=user input for daily F10.7 index
C          if oarr(41).le.0 then 12-month running mean is 
C          taken from internal file]
C    jf(27) =.flase.     OARR(39)=user input for IG12
C    jf(28) =.flase.     OARR(41)=user input for daily F10.7 index
C
C
C  OUTPUT:  OUTF(1:20,1:1000)
C               OUTF(1,*)  ELECTRON DENSITY/M-3
C               OUTF(2,*)  NEUTRAL TEMPERATURE/K
C               OUTF(3,*)  ION TEMPERATURE/K
C               OUTF(4,*)  ELECTRON TEMPERATURE/K
C               OUTF(5,*)  O+ ION DENSITY/% or /M-3 if jf(22)=f 
C               OUTF(6,*)  H+ ION DENSITY/% or /M-3 if jf(22)=f
C               OUTF(7,*)  HE+ ION DENSITY/% or /M-3 if jf(22)=f
C               OUTF(8,*)  O2+ ION DENSITY/% or /M-3 if jf(22)=f
C               OUTF(9,*)  NO+ ION DENSITY/% or /M-3 if jf(22)=f
C                 AND, IF JF(6)=.FALSE.:
C               OUTF(10,*)  CLUSTER IONS DEN/% or /M-3 if jf(22)=f
C               OUTF(11,*)  N+ ION DENSITY/% or /M-3 if jf(22)=f
C               OUTF(12,*)  
C               OUTF(13,*)  
C  if(jf(24)    OUTF(14,1:11) standard IRI-Ne for 60,65,..,110km 
C     =.false.)        12:22) Friedrich (FIRI) model at these heights 
C                      23:33) standard Danilov (SW=0, WA=0) 
C                      34:44) for minor Stratospheric Warming (SW=0.5) 
C                      45:55) for major Stratospheric Warming (SW=1) 
C                      56:66) weak Winter Anomaly (WA=0.5) conditions
C                      67:77) strong Winter Anomaly (WA=1) conditions
C               OUTF(15-20,*)  free
c
C            OARR(1:100)   ADDITIONAL OUTPUT PARAMETERS         
C
C      #OARR(1) = NMF2/M-3           #OARR(2) = HMF2/KM
C      #OARR(3) = NMF1/M-3           #OARR(4) = HMF1/KM
C      #OARR(5) = NME/M-3            #OARR(6) = HME/KM
C       OARR(7) = NMD/M-3             OARR(8) = HMD/KM
C       OARR(9) = HHALF/KM            OARR(10) = B0/KM
C       OARR(11) =VALLEY-BASE/M-3     OARR(12) = VALLEY-TOP/KM
C       OARR(13) = TE-PEAK/K          OARR(14) = TE-PEAK HEIGHT/KM
C      #OARR(15) = TE-MOD(300KM)     #OARR(16) = TE-MOD(400KM)/K
C       OARR(17) = TE-MOD(600KM)      OARR(18) = TE-MOD(1400KM)/K
C       OARR(19) = TE-MOD(3000KM)     OARR(20) = TE(120KM)=TN=TI/K
C       OARR(21) = TI-MOD(430KM)      OARR(22) = X/KM, WHERE TE=TI
C       OARR(23) = SOL ZENITH ANG/DEG OARR(24) = SUN DECLINATION/DEG
C       OARR(25) = DIP/deg            OARR(26) = DIP LATITUDE/deg
C       OARR(27) = MODIFIED DIP LAT.  OARR(28) = Geographic latitude
C       OARR(29) = sunrise/dec. hours OARR(30) = sunset/dec. hours
C       OARR(31) = ISEASON (1=spring) OARR(32) = Geographic longitude
C      #OARR(33) = Rz12               OARR(34) = Covington Index
C       OARR(35) = B1                 OARR(36) = M(3000)F2
C      $OARR(37) = TEC/m-2           $OARR(38) = TEC_top/TEC*100.
C      #OARR(39) = gind (IG12)        OARR(40) = F1 probability 
C      #OARR(41) = F10.7 daily        OARR(42) = c1 (F1 shape)
C       OARR(43) = daynr              OARR(44) = equatorial vertical 
C       OARR(45) = foF2_storm/foF2_quiet         ion drift in m/s
C      #OARR(46) = F10.7_81           OARR(47) = foE_storm/foE_quiet 
C       OARR(48) = spread-F probability          
C       OARR(49) = Geomag. latitude   OARR(50) = Geomag. longitude  
C       OARR(51) = ap at current time OARR(52) = daily ap
C       OARR(53) = invdip/degree      OARR(54) = MLT-Te
C       OARR(55) = CGM-latitude       OARR(56) = CGM-longitude
C       OARR(57) = CGM-MLT            OARR(58) = CGM lat eq. aurl bodry
C       OARR(59) = CGM-lati(MLT=0)    OARR(60) = CGM-lati for MLT=1
C       OARR(61) = CGM-lati(MLT=2)    OARR(62) = CGM-lati for MLT=3
C       OARR(63) = CGM-lati(MLT=4)    OARR(64) = CGM-lati for MLT=5
C       OARR(65) = CGM-lati(MLT=6)    OARR(66) = CGM-lati for MLT=7
C       OARR(67) = CGM-lati(MLT=8)    OARR(68) = CGM-lati for MLT=9
C       OARR(69) = CGM-lati(MLT=10)   OARR(70) = CGM-lati for MLT=11
C       OARR(71) = CGM-lati(MLT=12)   OARR(72) = CGM-lati for MLT=13
C       OARR(73) = CGM-lati(MLT=14)   OARR(74) = CGM-lati for MLT=15
C       OARR(75) = CGM-lati(MLT=16)   OARR(76) = CGM-lati for MLT=17
C       OARR(77) = CGM-lati(MLT=18)   OARR(78) = CGM-lati for MLT=19
C       OARR(79) = CGM-lati(MLT=20)   OARR(80) = CGM-lati for MLT=21
C       OARR(81) = CGM-lati(MLT=22)   OARR(82) = CGM-lati for MLT=23
C       OARR(83) = Kp at current time OARR(84) = magnetic declination 
C                # INPUT as well as OUTPUT parameter
C                $ special for IRIWeb (only place-holders)
c-----------------------------------------------------------------------        
C*****************************************************************
C*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
C***     ELECTRON DENSITY         60/80 KM       1500 KM       ***
C***     TEMPERATURES               60 KM        2500/3000 KM  ***
C***     ION DENSITIES             100 KM        1500 KM       ***
C*****************************************************************
C*****************************************************************
C*********            INTERNALLY                    **************
C*********       ALL ANGLES ARE IN DEGREE           **************
C*********       ALL DENSITIES ARE IN M-3           **************
C*********       ALL ALTITUDES ARE IN KM            **************
C*********     ALL TEMPERATURES ARE IN KELVIN       **************
C*********     ALL TIMES ARE IN DECIMAL HOURS       **************
C*****************************************************************
C*****************************************************************
C*****************************************************************
      INTEGER    DAYNR,DDO,DO2,SEASON,SEADAY
      REAL       LATI,LONGI,MO2,MO,MODIP,NMF2,MAGBR,INVDIP,IAPO,  
     &           NMF1,NME,NMD,MM,MLAT,MLONG,NMF2S,NMES
      CHARACTER  FILNAM*12
c-web-for webversion
c      CHARACTER FILNAM*53

      DIMENSION  ARIG(3),RZAR(3),F(3),E(4),XDELS(4),DNDS(4),
     &  FF0(988),XM0(441),F2(13,76,2),FM3(9,49,2),ddens(5,11),
     &  elg(7),FF0N(988),XM0N(441),F2N(13,76,2),FM3N(9,49,2),
     &  INDAP(13),AMP(4),HXL(4),SCL(4),XSM(4),MM(5),DTI(4),AHH(7),
     &  STTE(6),DTE(5),ATE(7),TEA(6),XNAR(2),param(2),OARR(100),
     &  OUTF(20,1000),DDO(4),DO2(2),DION(7),
     &  osfbr(25),D_MSIS(9),T_MSIS(2),IAPO(7),SWMI(25),ab_mlat(48),
     &  DAT(11,4), PLA(4), PLO(4)

      LOGICAL  EXT,SCHALT,TECON(2),sam_mon,sam_yea,sam_ut,sam_date,
     &  F1REG,FOF2IN,HMF2IN,URSIF2,LAYVER,DY,DREG,rzino,FOF1IN,
     &  HMF1IN,FOEIN,HMEIN,RZIN,sam_doy,F1_OCPRO,F1_L_COND,NODEN,
     &  NOTEM,NOION,TENEOP,OLD79,JF(50),URSIFO,igin,igino,
     &  dnight,enight,fnight,TOPO,TOPC,fstorm_on,estorm_on

      COMMON /CONST/UMR  /const1/humr,dumr   /ARGEXP/ARGMAX
     &	 /const2/icalls,nmono,iyearo,idaynro,rzino,igino,ut0
     &   /BLOCK1/HMF2,NMF2,HMF1,F1REG  /BLOCK2/B0,B1,C1  
     &   /BLOCK3/HZ,T,HST              /BLOCK4/HME,NME,HEF 
     &   /BLOCK5/ENIGHT,E              /BLOCK6/HMD,NMD,HDX
     &   /BLOCK7/D1,XKK,FP30,FP3U,FP1,FP2
     &   /BLOCK8/HS,TNHS,XSM,MM,DTI,MXSM       
c     &   /BLOTN/XSM1,TEXOS,TLBDH,SIGMA /BLOTE/AHH,ATE1,STTE,DTE
     &   /BLOTE/AHH,ATE1,STTE,DTE
     &   /BLO10/BETA,ETA,DELTA,ZETA    /findRLAT/FLON,RYEAR   
     &   /BLO11/B2TOP,TC3,itopn,alg10,hcor1       
     &   /iounit/konsol/QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU

      DATA SWMI/25*1./

      EXTERNAL          XE1,XE2,XE3_1,XE4_1,XE5,XE6,FMODIP

        save
        
        nummax=1000
        DO 7397 KI=1,20
        do 7397 kk=1,nummax
7397    OUTF(KI,kk)=-1.
C
C oarr(1:6,15,16,33,39:41) is used for inputs
C 
        do 8398 kind=7,14,1
8398    oarr(kind)=-1.
        do 8378 kind=17,32,1
8378    oarr(kind)=-1.
        do 8478 kind=34,38,1
8478    oarr(kind)=-1.
        oarr(40)=-1.
        do 8428 kind=42,100,1
8428    if(kind.ne.46) oarr(kind)=-1.

C
C PROGRAM CONSTANTS
C
        if(icalls.lt.1) then
        	ARGMAX=88.0
        	pi=ATAN(1.0)*4.
        	UMR=pi/180.
        	humr=pi/12.
        	dumr=pi/182.5
        	ALOG2=ALOG(2.)
        	ALG10=ALOG(10.)
        	ALG100=ALOG(100.)
        	endif
 
        numhei=int(abs(heiend-heibeg)/abs(heistp))+1
        if(numhei.gt.nummax) numhei=nummax
C
C NEW-GUL------------------------------
         Y05=.6931473
         QF=1.
	     h05top=0.
C NEW-GUL------------------------------

C
C Code inserted to aleviate block data problem for PC version.
C Thus avoiding DATA statement with parameters from COMMON block.
C
        XDELS(1)=5.
        XDELS(2)=5.
        XDELS(3)=5.
        XDELS(4)=10.
        DNDS(1)=.016
        DNDS(2)=.01
        DNDS(3)=.016
        DNDS(4)=.016
        DDO(1)=9
        DDO(2)=5
        DDO(3)=5
        DDO(4)=25
        DO2(1)=5
        DO2(2)=5
        XNAR(1)=0.0
        XNAR(2)=0.0
        DTE(1)=5.
        DTE(2)=5.
        DTE(3)=10.
        DTE(4)=20.
        DTE(5)=20.
        DTI(1)=10.
        DTI(2)=10.
        DTI(3)=20.
        DTI(4)=20.
C
C FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
C AGNR=OUTPUT (OUTPUT IS DISPLAYED OR STORED IN FILE OUTPUT.IRI)...
C IUCCIR=UNIT NUMBER FOR CCIR COEFFICIENTS ........................
c
        IUCCIR=10
c-web- special for web version
c-web- messages should be turned off with jf(34)=.false. and 
c-web- konsol=1 used in IRIFUN

        KONSOL=6
        if(.not.jf(12)) then
                konsol=11
                open(11,file='messages.txt')
                endif
        if(.not.jf(34)) KONSOL=1
c
c selection of density, temperature and ion composition options ......
c

      NODEN=(.not.jf(1))
      NOTEM=(.not.jf(2))
      NOION=(.not.jf(3))
      if(.not.NOION) NODEN=.false.
      DY=(.not.jf(6))
      LAYVER=(.not.jf(11))
      OLD79=(.not.jf(7))
      F1_OCPRO=(jf(19))
      F1_L_COND=(.not.jf(20))
      DREG=jf(24)
      TOPO=jf(29)
      TOPC=jf(30)
c
c rz12, IG12, F10.7D, PF10.7 input option ............................
c
      RZIN=(.not.jf(17))
      IF(RZIN) THEN
          ARZIN=OARR(33)
      else
          oarr(33)=-1.
      ENDIF
      
      IGIN=(.not.jf(27))
      IF(IGIN) THEN
          AIGIN=OARR(39)
      else
          oarr(39)=-1.
      ENDIF

      IF(.not.jf(25)) THEN
          f107din=OARR(41)
      else
          oarr(41)=-1.
      ENDIF

      IF(.not.jf(32)) THEN
          f10781in=OARR(46)
      else
          oarr(46)=-1.
      ENDIF
c
c Topside density ....................................................
c
        if(TOPO) then
             if (TOPC) then 
                 itopn=0
             else
                 itopn=3
             endif
        else 
             if (TOPC) then
                 itopn=1
             else
                 itopn=2
             endif
        endif
c
c F2 peak density ....................................................
c
      FOF2IN=(.not.jf(8))
       IF(FOF2IN) THEN
          OARR1=OARR(1)
       	  AFOF2=OARR1
       	  ANMF2=OARR1
          IF(OARR1.LT.100.) ANMF2=1.24E10*AFOF2*AFOF2
          IF(OARR1.GE.100.) AFOF2=SQRT(ANMF2/1.24E10)
        else
          oarr(1)=-1.
          ENDIF
      URSIF2=(.not.jf(5))
c
c F2 peak altitude ..................................................
c
      HMF2IN=(.not.jf(9))
       IF(HMF2IN) then
                AHMF2=OARR(2)
        else
                oarr(2)=-1.
        endif
c
c F1 peak density ...................................................
c
      FOF1IN=(.not.jf(13))
       IF(FOF1IN) THEN
          OARR3=OARR(3)
       	  AFOF1=OARR3
       	  ANMF1=OARR3
          IF(OARR3.LT.100.) ANMF1=1.24E10*AFOF1*AFOF1
          IF(OARR3.GE.100.) AFOF1=SQRT(ANMF1/1.24E10)
        else
          oarr(3)=-1.
          ENDIF
c
c F1 peak altitude ..................................................
c
      HMF1IN=(.not.jf(14))
       IF(HMF1IN) then
                AHMF1=OARR(4)
                if(.not.layver.and.jf(34)) write(konsol,1939)
1939  format(' *Ne* User input of hmF1 is only possible for the LAY-',
     &          'version')
        else
                oarr(4)=-1.
        endif
c
c E peak density ....................................................
c
      FOEIN=(.not.jf(15))
       IF(FOEIN) THEN
          OARR5=OARR(5)
       	  AFOE=OARR5
       	  ANME=OARR5
          IF(OARR5.LT.100.) ANME=1.24E10*AFOE*AFOE
          IF(OARR5.GE.100.) AFOE=SQRT(ANME/1.24E10)
        else
          oarr(5)=-1.
        ENDIF
c
c E peak altitude ..................................................
c
      HMEIN=(.not.jf(16))
       IF(HMEIN) then
                AHME=OARR(6)
        else
                oarr(6)=-1.
        endif
c
C TE-NE MODEL OPTION ..............................................
C
      TENEOP=(.not.jf(10))
        IF(TENEOP) THEN
           DO 8154 JXNAR=1,2
              XNAR(JXNAR)=OARR(JXNAR+14)
              TECON(JXNAR)=.FALSE. 
8154          IF(XNAR(JXNAR).GT.0.) TECON(JXNAR)=.TRUE. 
        else
           oarr(15)=-1.
           oarr(16)=-1.
           ENDIF
c
c lists the selected options before starting the table
c

      if(icalls.gt.1.or.(.not.jf(34))) goto 8201
          write(konsol,2911) 
        if(NODEN) goto 2889
          if(LAYVER) write(konsol,9012) 
          if(OLD79) write(konsol,9014) 
          if (itopn.eq.0) write(konsol,9207)
          if (itopn.eq.1) write(konsol,9204)
          if (itopn.eq.2) write(konsol,9205)
          if (itopn.eq.3) write(konsol,9206)
          if(FOF2IN) then
                write(konsol,9015) 
                goto 2889
                endif
        if(URSIF2) then
                write(konsol,9016) 
        else
                write(konsol,9017) 
        endif
        if(HMF2IN) write(konsol,9018) 
        if(fof1in) write(konsol,9019) 
        if(HMF1IN.and.LAYVER) write(konsol,9021) 
        if(jf(4)) then
             write(konsol,9214)
        else 
             if (jf(31)) then
                 write(konsol,9216)
             else
                 write(konsol,9215)
             endif
        endif
        if(foein) write(konsol,9022) 
        if(HMEIN) write(konsol,9023) 
        if(F1_OCPRO) write(konsol,9024) 
        if(F1_L_COND) write(konsol,9025) 
        if(DREG) then 
            write(konsol,9026) 
        else
            write(konsol,9027) 
        endif
        if(jf(26)) then
            if(fof2in) then 
                  write(konsol,9028) 
                  jf(26)=.false.
            else
                  write(konsol,9029) 
            endif
            endif

2889    continue

        if((.not.NOION).and.(DY)) write(konsol,9031) 
        if((.not.NOION).and.(.not.DY)) write(konsol,9039) 

        if(NOTEM) goto 8201
          if(TENEOP) write(konsol,9032) 
          if(jf(23)) then 
            write(konsol,9033) 
          else
            write(konsol,9034) 
          endif

        if(jf(33)) then
            write(konsol,4031)
        else
            write(konsol,4032)
        endif

        if(jf(35)) then
            write(konsol,9128)
        else
            write(konsol,9129)
        endif

2911    format('*** IRI parameters are being calculated ***')
9012    format('Ne, E-F: The LAY-Version is prelimenary.',
     &          ' Erroneous profile features can occur.')
9014    format('Ne: No upper limit for F10.7 in',
     &          ' topside formula.')
9204    format('Ne: Corrected Topside Formula')
9205    format('Ne: NeQuick for Topside')
9206    format('Ne: Gul-h0.5 for Topside')
9207    format('Ne: IRI-2001 for Topside')
9214    format('Ne: B0,B1 Bil-2000')
9215    format('Ne: B0 Gul-1987')
9216    format('Ne: B0,B1-ABT-2009')
9015    format('Ne, foF2/NmF2: provided by user.')
9016    format('Ne, foF2: URSI model is used.')
9017    format('Ne, foF2: CCIR model is used.')
9018    format('Ne, hmF2/M3000F2: provided by user.')
9019    format('Ne, foF1/NmF1: provided by user.')
9021    format('Ne, hmF1: provided by user.')
9022    format('Ne, foE/NmE: provided by user.')
9023    format('Ne, hmE: provided by user.')
9024    format('Ne, foF1: probability function used.')
9025    format('Ne, foF1: L condition cases included.')
9026    format('Ne, D: IRI1990')
9027    format('Ne, D: FT2001; IRI-90, FT-01, DRS-95)')
9028    format('Ne, foF2: Storm model turned off if foF2 or',
     &          ' NmF2 user input')
9029    format('Ne, foF2: storm model included')
9128    format('Ne, foE: storm model on')
9129    format('Ne, foE: storm model off')
9039    format('Ion Com.: DS-95 & DY-85')
9031    format('Ion Com.: RBV-10 & TTS-03')
9032    format('Te: Temperature-density correlation is used.')
9033    format('Te: Aeros/AE/ISIS model')
9034    format('Te: TBT-2012 model')
4031    format('Auroral boundary model on')
4032    format('Auroral boundary model off')

8201    continue

C
C CALCULATION OF DAY OF YEAR OR MONTH/DAY AND DECIMAL YEAR 
c NRDAYM is the number of days in the current month 
c IDAYY is the number of days in the current year
c
c  leap year rule: years evenly divisible by 4 are leap years, except
c  years also evenly divisible by 100 are not leap years, except years 
c  also evenly divisible by 400 are leap years. The year 2000 is a 100 
c  and 400 year exception and therefore it is a normal leap year. 
c  The next 100 year exception will be in the year 2100!
c

        iyear=iyyyy
        if(iyear.lt.100) iyear=iyear+1900
        if(iyear.lt.30) iyear=iyear+2000
        idayy=365
        if(iyear/4*4.eq.iyear) idayy=366    ! leap year

        if(MMDD.lt.0) then
                DAYNR=-MMDD
                call MODA(1,iyear,MONTH,IDAY,DAYNR,nrdaym)
        else
                MONTH=MMDD/100
                IDAY=MMDD-MONTH*100
                call MODA(0,iyear,MONTH,IDAY,DAYNR,nrdaym)
        endif

        ryear = iyear + (daynr-1.0)/idayy
        iyd = iyear *1000 + daynr

C
C calculate center height for CGM computation
C

        height_center=(HEIBEG+HEIEND)/2.
        

C
C CALCULATION OF GEODETIC/GEOMAGNETIC COORDINATES (LATI, LONGI AND 
C MLAT, MLONG), MAGNETIC INCLINATION (DIP), DIP LATITUDE (MAGBR) 
C AND MODIFIED DIP (MODIP), ALL IN DEGREES
C

        if(along.lt.0.) along = along + 360. ! -180/180 to 0-360
        
        IF(JMAG.GT.0) THEN
           MLAT=ALATI
           MLONG=ALONG
        ELSE
           LATI=ALATI
           LONGI=ALONG
        ENDIF
        CALL GEODIP(IYEAR,LATI,LONGI,MLAT,MLONG,JMAG)
        call igrf_dip(lati,longi,ryear,300.0,dec,dip,magbr,modip)
        if(.not.jf(18)) then
        CALL FIELDG(LATI,LONGI,300.0,XMA,YMA,ZMA,BET,DIP,DEC,MODIP)
        	MAGBR=ATAN(0.5*TAN(DIP*UMR))/UMR
        	endif

        ABSLAT=ABS(LATI)
        ABSMLT=ABS(MLAT)
        ABSMDP=ABS(MODIP)
        ABSMBR=ABS(MAGBR)

c
C CALCULATION OF UT/LT and XMLT  ...............
c
        IF(DHOUR.le.24.0) then
        	HOUR=DHOUR					! dhour =< 24 is LT
            hourut=hour-longi/15.
            if(hourut.lt.0) hourut=hourut+24.
		else
        	hourut=DHOUR-25.				 ! dhour>24 is UT+25
        	hour=hourut+longi/15.	 		 ! hour becomes LT
        	if(hour.gt.24.) hour=hour-24.
        endif
        
        CALL CLCMLT(IYEAR,DAYNR,HOURUT,LATI,LONGI,XMLT)

c
c SEASON assumes equal length seasons (92 days) with spring 
c (SEASON=1) starting at day-of-year=45; for lati < 0 adjustment 
c for southern hemisphere is made. Some models require the
c seasonal month (ISEAMON) or the seasonal day-of year (SEADAY)
c ZMONTH is decimal month (Jan 1 = 1.0 and Dec 31 = 12.97)
c SDAY is the day number reduced to a 360 day year (TOPH05)
c NRDAYM is the number of days in the current month 
c IDAYY is the number of days in the current year
c 
      
      SEASON=INT((DAYNR+45.0)/92.0)
      IF(SEASON.LT.1) SEASON=4
      NSEASN=SEASON				! Northern hemisphere season
      zmonth = month + (iday-1)*1./nrdaym
C NEW-GUL------------------------------
      sday=daynr/idayy*360.			 
C NEW-GUL------------------------------
      seaday=daynr
      iseamon=month
      IF(LATI.GE.0.0) GOTO 5592
        	SEASON=SEASON-2			
        	IF(SEASON.LT.1) SEASON=SEASON+4
        	iseamon=month+6
        	if(iseamon.gt.12) iseamon=iseamon-12
        	seaday=daynr+idayy/2.
        	if(seaday.gt.idayy) seaday=seaday-idayy
C NEW-GUL------------------------------
            sday=sday+180.						
            if (sday.gt.360.) sday=sday-360.	
C NEW-GUL------------------------------

C
C 12-month running mean sunspot number (rssn) and Ionospheric Global 
C index (gind), daily F10.7 cm solar radio flux (f107d) and monthly 
C F10.7 (cov) index   
C

5592    continue
        sam_mon=(month.eq.montho)
        sam_yea=(iyear.eq.iyearo)
        sam_doy=(daynr.eq.idaynro)
        sam_date=(sam_yea.and.sam_doy)
        sam_ut=(hourut.eq.ut0)
        
        if(sam_date.and..not.rzino.and..not.rzin.
     &                   and..not.igin.and..not.igino) goto 2910
     
        call tcon(iyear,month,iday,daynr,rzar,arig,ttt,nmonth)
        if(nmonth.lt.0) goto 3330		! jump to end of program

        if(RZIN) then
        	rrr = arzin
        	rzar(1) = rrr
        	rzar(2) = rrr
        	rzar(3) = rrr
c       	zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
c       	if(zi.gt.174.0) zi=174.0
c       	arig(1) = zi
c       	arig(2) = zi
c       	arig(3) = zi
        	endif
        if(IGIN) then
        	zi = aigin
        	arig(1) = zi
        	arig(2) = zi
        	arig(3) = zi
        	endif
        rssn=rzar(3)
        gind=arig(3)
        COV=63.75+RSSN*(0.728+RSSN*0.00089)
c        rlimit=gind
c        COVSAT=63.75+rlimit*(0.728+rlimit*0.00089)
        COVSAT=cov
        if(covsat.gt.188.) covsat=188
        
        f107d=cov
        f107pd=cov
        f10781=cov
        f107365=cov
        call APF_ONLY(iyear,month,iday,F107_daily,F107PD,F107_81,
     &        F107_365,IAP_daily)
        if(F107_daily.gt.-11.1) then
            f107d=f107_daily
            f107y=f107PD
            f10781=f107_81
            f107365=f107_365
			endif
        if(.not.jf(25)) then
        		f107d=f107din 		! F10.7D user input
        		f107y=f107din 		! F10.7PD user input
				endif
        if(.not.jf(32)) f10781=f10781in 	! F10.7_81 user input

		pf107 = (f107d+f10781)/2.
        				

C
C CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG), SUN DECLINATION ANGLE 
C (SUNDEC),SOLAR ZENITH ANGLE AT NOON (XHINON) AND TIME OF LOCAL 
C SUNRISE/SUNSET (SAX, SUX; dec. hours) AT 70 KM (D-REGION), 110 KM
C (E-REGION), 200 KM (F1-REGION), AND 500 KM (TE, TI).
C

2910    continue
        CALL SOCO(daynr,HOUR,LATI,LONGI,80.,SUNDEC,XHI1,SAX80,SUX80)
        CALL SOCO(daynr,HOUR,LATI,LONGI,110.,SUD1,XHI2,SAX110,SUX110)
        CALL SOCO(daynr,HOUR,LATI,LONGI,200.,SUD1,XHI,SAX200,SUX200)
        CALL SOCO(daynr,HOUR,LATI,LONGI,300.,SUD1,XHI3,SAX300,SUX300)
        CALL SOCO(daynr,12.0,LATI,LONGI,110.,SUNDE1,XHINON,SAX1,SUX1)
        DNIGHT=.FALSE.
        if(abs(sax80).gt.25.0) then
                if(sax80.lt.0.0) DNIGHT=.TRUE.
                goto 9334
                endif
        if(SAX80.le.SUX80) goto 1386
        if((hour.gt.sux80).and.(hour.lt.sax80)) dnight=.true.
        goto 9334
1386    IF((HOUR.GT.SUX80).OR.(HOUR.LT.SAX80)) DNIGHT=.TRUE.

9334    ENIGHT=.FALSE.
        if(abs(sax110).gt.25.0) then
                if(sax110.lt.0.0) ENIGHT=.TRUE.
                goto 8334
                endif
        if(SAX110.le.SUX110) goto 9386
        if((hour.gt.sux110).and.(hour.lt.sax110)) enight=.true.
        goto 8334
9386    IF((HOUR.GT.SUX110).OR.(HOUR.LT.SAX110)) ENIGHT=.TRUE.

8334    FNIGHT=.FALSE.
        if(abs(sax200).gt.25.0) then
                if(sax200.lt.0.0) FNIGHT=.TRUE.
                goto 1334
                endif
        if(SAX200.le.SUX200) goto 7386
        if((hour.gt.sux200).and.(hour.lt.sax200)) fnight=.true.
        goto 1334
7386    IF((HOUR.GT.SUX200).OR.(HOUR.LT.SAX200)) FNIGHT=.TRUE.

C
C CALCULATION OF ELECTRON DENSITY PARAMETERS................
C lower height boundary (HNEA), upper boundary (HNEE)
C
      
1334  continue
      HNEA=65.
      IF(DNIGHT) HNEA=80.
      HNEE=2000.
      IF(NODEN) GOTO 4933

      DELA=4.32
      IF(ABSMDP.GE.18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)
      DELL=1+EXP(-(ABSLAT-20.)/10.)
c
c E peak critical frequency (foE), density (NmE), and height (hmE)
c
        IF(FOEIN) THEN
          FOE=AFOE
          NME=ANME
        ELSE
          FOE=FOEEDI(COV,XHI,XHINON,ABSLAT)
          NME=1.24E10*FOE*FOE
        ENDIF
        IF(HMEIN) THEN
          HME=AHME
        ELSE
          HME=110.0
        ENDIF
c
c F2 peak critical frequency foF2, density NmF2, and height hmF2
c
C READ CCIR AND URSI COEFFICIENT SET FOR CHOSEN MONTH 
C
      IF((FOF2IN).AND.(HMF2IN).and.(itopn.ne.2)) GOTO 501
      IF(URSIF2.NEQV.URSIFO) GOTO 7797
      if(.not.rzin.and..not.rzino.and..not.igin.and..not.igino) then
          IF(sam_mon.AND.(nmonth.EQ.nmono).and.sam_yea) GOTO 4292
          IF(sam_mon) GOTO 4293
          endif
c
c the program expects the coefficients files in ASCII format; if you
C want to use the binary version of the coefficients, please use the
C the statements that are commented-out below and comment-out the
C ASCII-related statements.
c
7797    URSIFO=URSIF2
        WRITE(FILNAM,104) MONTH+10
104         FORMAT('ccir',I2,'.asc')
c-binary- if binary files than use:
c-binary-104   FORMAT('ccir',I2,'.bin')
c-web- special for web-version:
c104   FORMAT('/usr/local/etc/httpd/cgi-bin/models/IRI/ccir',I2,'.asc')
c-web- special for web-version:
c
        OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=8448,
     &          FORM='FORMATTED')
c-binary- if binary files than use:
c-binary-     &          FORM='UNFORMATTED')
c
        READ(IUCCIR,4689) F2,FM3
4689    FORMAT(1X,4E15.8)
c-binary- if binary files than use:
c-binary-        READ(IUCCIR) F2,FM3

c
        CLOSE(IUCCIR)
C
C then URSI if chosen ....................................
C
        if(URSIF2) then
          WRITE(FILNAM,1144) MONTH+10
1144          FORMAT('ursi',I2,'.asc')
c-web- special for web-version:
c1144  FORMAT('/usr/local/etc/httpd/cgi-bin/models/IRI/ursi',I2,'.asc')
c-binary- if binary files than use:
c-binary-1144          FORMAT('ursi',I2,'.bin')

          OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=8448,
     &         FORM='FORMATTED')
c-binary- if binary files than use:
c-binary-     &         FORM='UNFORMATTED')

          READ(IUCCIR,4689) F2
c-binary- if binary files than use:
c-binary-          READ(IUCCIR) F2

          CLOSE(IUCCIR)
        endif

C
C READ CCIR AND URSI COEFFICIENT SET FOR NMONTH, i.e. previous 
c month if day is less than 15 and following month otherwise 
C

4293    continue

c
c first CCIR ..............................................
c

        WRITE(FILNAM,104) NMONTH+10
        OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=8448,
     &          FORM='FORMATTED')
c-binary- if binary files than use:
c-binary-     &          FORM='unFORMATTED')

        READ(IUCCIR,4689) F2N,FM3N
c-binary- if binary files than use:
c-binary-        READ(IUCCIR) F2N,FM3N

        CLOSE(IUCCIR)

C
C then URSI if chosen .....................................
C
        if(URSIF2) then
          WRITE(FILNAM,1144) NMONTH+10
          OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=8448,
     &         FORM='FORMATTED')
c-binary- if binary files than use:
c-binary-     &         FORM='unFORMATTED')

          READ(IUCCIR,4689) F2N
c-binary- if binary files than use:
c-binary-          READ(IUCCIR) F2N

          CLOSE(IUCCIR)
          endif

        GOTO 4291
        
8448    WRITE(konsol,8449) FILNAM
8449    FORMAT(1X////,
     &    ' The file ',A30,'is not in your directory.')
        GOTO 3330
C
C LINEAR INTERPOLATION IN SOLAR ACTIVITY. IG12 used for foF2
C

4291    continue
        RR2=ARIG(1)/100.
        RR2N=ARIG(2)/100.
        RR1=1.-RR2
        RR1N=1.-RR2N
        DO 20 I=1,76
        DO 20 J=1,13
              K=J+13*(I-1)
              FF0N(K)=F2N(J,I,1)*RR1N+F2N(J,I,2)*RR2N
20            FF0(K)=F2(J,I,1)*RR1+F2(J,I,2)*RR2

        RR2=RZAR(1)/100.
        RR2N=RZAR(2)/100.
        RR1=1.-RR2
        RR1N=1.-RR2N
        DO 30 I=1,49
        DO 30 J=1,9
              K=J+9*(I-1)
              XM0N(K)=FM3N(J,I,1)*RR1N+FM3N(J,I,2)*RR2N
30            XM0(K)=FM3(J,I,1)*RR1+FM3(J,I,2)*RR2

4292    zfof2  =  FOUT(MODIP,LATI,LONGI,HOURUT,FF0)
        fof2n  =  FOUT(MODIP,LATI,LONGI,HOURUT,FF0N)
        zm3000 = XMOUT(MODIP,LATI,LONGI,HOURUT,XM0)
        xm300n = XMOUT(MODIP,LATI,LONGI,HOURUT,XM0N)
        midm=15
        if(month.eq.2) midm=14
        if (iday.lt.midm) then
                yfof2 = fof2n + ttt * (zfof2-fof2n)
                xm3000= xm300n+ ttt * (zm3000-xm300n)
        else
                yfof2 = zfof2 + ttt * (fof2n-zfof2)
                xm3000= zm3000+ ttt * (xm300n-zm3000)
        endif
			
501     IF(FOF2IN) THEN
          FOF2=AFOF2
          NMF2=ANMF2
        ELSE
          FOF2=YFOF2
          NMF2=1.24E10*FOF2*FOF2
        ENDIF
c
c stormtime updating for foF2 (foF2s, NmF2s) and foE (foEs,
c NmEs) and auroral boundary computation.
c
        foF2s=foF2
        foEs=foE
        NmF2s=NmF2
        NmEs=NmE
        stormcorr=-1.
        estormcor=-1.
        fstorm_on=jf(26).and.jf(8)
        estorm_on=jf(35).and.jf(15)
        if((fstorm_on.or.jf(33).or.estorm_on).and.
     &      (.not.sam_date.or..not.sam_ut)) 
     &      call apf(iyear,month,iday,hourut,indap)
c
c stormtime updating for foF2 (foF2s, NmF2s) 
c
        if(fstorm_on.and.(indap(1).gt.-1)) then    
            icoord=1
            kut=int(hourut)
            call STORM(indap,lati,longi,icoord,cglat,kut,
     &              daynr,stormcorr)
            fof2s=fof2*stormcorr
            NMF2S=1.24E10*FOF2S*FOF2S
            endif
c
c stormtime updating for foE (foEs, NmEs)
c
        if(estorm_on.and.(indap(1).gt.-1)) then    
            estormcor=STORME_AP(DAYNR,MLAT,INDAP(13)*1.0)
            if(estromcor.gt.-2.0) foes=foe*estormcor
            NMES=1.24E10*FOES*FOES
            endif
c
c calculation of equatorward auroral boundary
c
        if(jf(33)) then
            if(indap(1).gt.-1) then
 	            xkp=ckp(indap(13))
            else     
                xkp=3.0
	        endif   
c Corrected magnetic latitude CGM of equatorward boundary, 
c ab_mlat(48), for MLT=0.0,0.5,1.0 ... 23.5 h and kp=xkp
            call auroral_boundary(xkp,-1.0,cgmlat,ab_mlat)
	        DAT(1,1)=lati
	        DAT(2,1)=longi
            call GEOCGM01(1,IYEAR,height_center,DAT,PLA,PLO)
            cgm_lat=DAT(3,1)
            cgm_lon=DAT(4,1)
            cgm_mlt00_ut=DAT(11,1)
            cgm_mlt=hourut-cgm_mlt00_ut
            if(cgm_mlt.lt.0.) cgm_mlt=24.+hourut-cgm_mlt00_ut
C
C CGM latitude of boundary (cgmlat) for present MLT value
C 2012.02 12/17/12 Add magnetic declination as oarr(84) output
            zmlt=xmlt
C            zmlt=cgm_mlt
            cgmlat=100.0
            if(zmlt.ge.0.0.and.zmlt.le.24.0) 
     &          call auroral_boundary(xkp,zmlt,cgmlat,ab_mlat)
            endif
c
c computing hmF2. set jf(36)=false to use foF2_storm in hmF2 formula
c  
        IF(HMF2IN) THEN
          IF(AHMF2.LT.50.0) THEN
          	XM3000=AHMF2
          	HMF2=HMF2ED(MAGBR,RSSN,FOF2/FOE,XM3000)
          ELSE
          	HMF2=AHMF2
c          	XM3000=XM3000HM(MAGBR,RSSN,FOF2/FOE,HMF2)
          ENDIF
        ELSE
          ratf=fof2s/foe
          if(jf(36)) ratf=fof2/foe
          HMF2=HMF2ED(MAGBR,RSSN,RATF,XM3000)
        ENDIF


        nmono=nmonth
        MONTHO=MONTH
        iyearo=iyear
        idaynro=daynr
        rzino=rzin
        igino=igin
        ut0=hourut

c
c topside profile parameters .............................
c
        COS2=COS(MLAT*UMR)
        COS2=COS2*COS2
        FLU=(COVSAT-40.0)/30.0
c
c option to use unlimited F10.7M for the topside
c previously: IF(OLD79) ETA1=-0.0070305*COS2
c
        IF(OLD79) FLU=(COV-40.0)/30.0
        FO1 = FOF2S
        IF(JF(37)) FO1 = FOF2
        EX=EXP(-MLAT/15.)
        EX1=EX+1
        EPIN=4.*EX/(EX1*EX1)
        ETA1=-0.02*EPIN
        ETA = 0.058798 + ETA1 - 
     &    FLU * (0.014065  - 0.0069724 * COS2) + 
     &    FO1* (0.0024287 + 0.0042810 * COS2  - 0.0001528 * FO1)
        ZETA = 0.078922 - 0.0046702 * COS2 -  
     &    FLU * (0.019132  - 0.0076545 * COS2) +
     &    FO1* (0.0032513 + 0.0060290 * COS2  - 0.00020872 * FO1)
        BETA=-128.03 + 20.253 * COS2 -
     &    FLU * (8.0755  + 0.65896 * COS2) +
     &    FO1* (0.44041 + 0.71458 * COS2 - 0.042966 * FO1)
        Z=EXP(94.5/BETA)
        Z1=Z+1
        Z2=Z/(BETA*Z1*Z1)
        DELTA=(ETA/Z1-ZETA/2.0)/(ETA*Z2+ZETA/400.0)
c
c Correction term for topside (Bilitza)
c
      if(itopn.eq.1) then
          zmp1 = exp(modip / 10.)
          zmp11 = 1. + zmp1
          zmp111 = zmp1 / (zmp11 * zmp11)
          zmp2 = exp(modip / 19.)
          zmp22 = 1. + zmp2
          zmp222 = zmp2 / (zmp22 * zmp22) 
          r2n = -0.84 - 1.6 * zmp111  
          r2d = -0.84 - 0.64 * zmp111 
          x1n = 230. - 700. * zmp222  
          x1d = 550. - 1900. * zmp222
          r2 = HPOL(HOUR,r2d,r2n,SAX300,SUX300,1.,1.)
          x1 = HPOL(HOUR,x1d,x1n,SAX300,SUX300,1.,1.)
          hcor1 = hmF2 + x1
          x12 = 1500. - x1
	      tc3 = r2 / x12
          endif
c NEW-GUL--------------------------------
c
c Correction term for topside (Gulyaeva)
c
      if(itopn.eq.3) then
          hei05=0.
          CALL TOPH05(COV,MLAT,HOUR,HMF2,HEI05,SDAY)
          h05top=hei05
          xnetop=XE_1(H05TOP)
          endif
c NEW-GUL--------------------------------

c
c NeQuick topside parameters (use CCIR-M3000F2 even if user-hmF2)
c
      if (itopn.eq.2) then
         fo2=foF2s
         if(jf(37)) fo2=foF2
         dNdHmx=-3.467+1.714*log(fo2)+2.02*log(xM3000)
         dNdHmx=exp(dNdHmx)*0.01
         B2bot=0.04774*fo2*fo2/dNdHmx
         b2k = 3.22-0.0538*fo2-0.00664*hmF2+0.113*hmF2/B2bot
     &   	+0.00257*rssn
         ee=exp(2.0*(b2k-1.0))
         b2k=(b2k*ee+1.0)/(ee+1.0)
         B2TOP=b2k*B2bot
      endif
c
c Bottomside thickness parameter B0 and shape parameters B1
c                              
        if(jf(4)) then
          B0=B0_98(HOUR,SAX200,SUX200,NSEASN,RSSN,LONGI,MODIP)
	      B1=HPOL(HOUR,1.9,2.6,SAX200,SUX200,1.,1.)
        else if (jf(31)) then
          FLON=LONGI+15.*hourut
          if(FLON.gt.360.) FLON=FLON-360.
          X11=-90
          X22=90
          FX11=fmodip(x11)
          FX22=fmodip(x22)
          CALL REGFA1(X11,X22,FX11,FX22,0.001,MODIP,FMODIP,SCHALT,XRLAT)
          IF(SCHALT) THEN
             XRLAT=LATI
             if(jf(34)) WRITE(KONSOL,656) LATI
656     FORMAT(1X,'*NE* ABT-B0 computed with RLAT=LATI=',F6.2)
             endif
          RLAT=XRLAT
 	      CALL SHAMDB0D (RLAT,FLON,ZMONTH,RSSN,B0)
          CALL SHAB1D (LATI,FLON,ZMONTH,RSSN,B1)
        else
          CALL ROGUL(SEADAY,XHI,SEAX,GRAT)
          IF (FNIGHT) GRAT = 0.91D0 - HMF2/4000.D0
	      B1=HPOL(HOUR,1.9,2.6,SAX200,SUX200,1.,1.)
          BCOEF = B1*(B1*(0.0046D0*B1-0.0548D0)+0.2546D0) + 0.3606D0
          B0CNEW = HMF2*(1.D0-GRAT)
          B0 = B0CNEW/BCOEF
        endif
c
c F1 layer height hmF1, critical frequency foF1, peak density NmF1
c
        IF(FOF1IN) THEN
            FOF1=AFOF1
            NMF1=ANMF1
        ELSE
            FOF1=FOF1ED(ABSMBR,RSSN,XHI)
            NMF1=1.24E10*FOF1*FOF1
        ENDIF
c
c F1 layer thickness parameter c1
c
        c1 = f1_c1(modip,hour,sux200,sax200)
c
c F1 occurrence probability with Scotto et al. 1997 or Ducharme et al. 
c if jf(19)=f1_ocpro=.true. or .false.
c If .not.jf(20)=f1_l_cond=.true. than Scotto model with L-condition
c
        if(f1_ocpro) then
        	call f1_prob(xhi,mlat,rssn,f1pbw,f1pbl)
            f1pb = f1pbw
            if(f1_l_cond) f1pb = f1pbl
            f1reg=.false.
            if((fof1in).or.(f1pb.ge.0.5)) f1reg=.true.
        else			
        	f1pb = 0.0 
        	if(fof1in.or.((.not.fnight).and.(fof1.gt.0.0))) f1pb=1. 
            f1reg=.false.
            if(f1pb.gt.0.0) f1reg=.true.
        endif
            
c
c E-valley: depth, width, height of deepest point (HDEEP),
c height of valley top (HEF)
c
      XDEL=XDELS(SEASON)/DELA
      DNDHBR=DNDS(SEASON)/DELA
      HDEEP=HPOL(HOUR,10.5/DELA,28.,SAX110,SUX110,1.,1.)
      WIDTH=HPOL(HOUR,17.8/DELA,45.+22./DELA,SAX110,SUX110,1.,1.)
      DEPTH=HPOL(HOUR,XDEL,81.,SAX110,SUX110,1.,1.)
      DLNDH=HPOL(HOUR,DNDHBR,.06,SAX110,SUX110,1.,1.)
      IF(DEPTH.LT.1.0) GOTO 600
        IF(ENIGHT) DEPTH=-DEPTH
        CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
        IF(.NOT.EXT) GOTO 667
        if(jf(34)) WRITE(KONSOL,650)
650     FORMAT(1X,'*NE* E-REGION VALLEY CAN NOT BE MODELLED')
600   WIDTH=.0
667   HEF=HME+WIDTH
      hefold=hef
      VNER = (1. - ABS(DEPTH) / 100.) * NMES

c
c Parameters below E  .............................
c

2727  continue
      hmex=hme-9.
      NMD=XMDED(XHI,RSSN,4.0E8)
      HMD=HPOL(HOUR,81.0,88.0,SAX80,SUX80,1.,1.)
      F(1)=HPOL(HOUR,0.02+0.03/DELA,0.05,SAX80,SUX80,1.,1.)
      F(2)=HPOL(HOUR,4.6,4.5,SAX80,SUX80,1.,1.)
      F(3)=HPOL(HOUR,-11.5,-4.0,SAX80,SUX80,1.,1.)
      FP1=F(1)
      FP2=-FP1*FP1/2.0
      FP30=(-F(2)*FP2-FP1+1.0/F(2))/(F(2)*F(2))
      FP3U=(-F(3)*FP2-FP1-1.0/F(3))/(F(3)*F(3))
      HDX=HMD+F(2)
c
c indermediate region between D and E region; parameters xkk
c and d1 are found such that the function reaches hdx/xdx/dxdh
c
       X=HDX-HMD
       XDX=NMD*EXP(X*(FP1+X*(FP2+X*FP30)))
       DXDX=XDX*(FP1+X*(2.0*FP2+X*3.0*FP30))
       X=HME-HDX
       XKK=-DXDX*X/(XDX*ALOG(XDX/NMES))
c
c if exponent xkk is larger than xkkmax, then xkk will be set to 
c xkkmax and d1 will be determined such that the point hdx/xdx is 
c reached; derivative is no longer continuous.
c
        xkkmax=5.
        if(xkk.gt.xkkmax) then
                xkk=xkkmax
                d1=-alog(xdx/nmes)/(x**xkk)
        else
                D1=DXDX/(XDX*XKK*X**(XKK-1.0))
        endif
c
c compute Danilov et al. (1995) D-region model values
c
      if(.not.dreg) then
          vKp=1.
          f5sw=0.
          f6wa=0.
          call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
          do ii=1,11
            ddens(1,ii)=-1.
            if(ii.lt.8) ddens(1,ii)=10**(elg(ii)+6)
            enddo
          f5sw=0.5
          f6wa=0.
          call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
          do ii=1,11
            ddens(2,ii)=-1.
            if(ii.lt.8) ddens(2,ii)=10**(elg(ii)+6)
            enddo
          f5sw=1.
          f6wa=0.
          call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
          do ii=1,11
            ddens(3,ii)=-1.
            if(ii.lt.8) ddens(3,ii)=10**(elg(ii)+6)
            enddo
          f5sw=0.
          f6wa=0.5
          call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
          do ii=1,11
            ddens(4,ii)=-1.
            if(ii.lt.8) ddens(4,ii)=10**(elg(ii)+6)
            enddo          
          f5sw=0.
          f6wa=1.
          call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
          do ii=1,11
            ddens(5,ii)=-1.
            if(ii.lt.8) ddens(5,ii)=10**(elg(ii)+6)
            enddo
          endif
C
C SEARCH FOR HMF1 ..................................................
C

       if(LAYVER) goto 6153
       hmf1=0
       IF(.not.F1REG) GOTO 380

c omit F1 feature if nmf1*0.9 is smaller than nme
       bnmf1=0.9*nmf1
       if(nmes.ge.bnmf1) goto 9427

9245   XE2H=XE2(HEF)
       if(xe2h.gt.bnmf1) then
            hef=hef-1
            if(hef.le.hme) goto 9427
            goto 9245
            endif      
        CALL REGFA1(HEF,HMF2,XE2H,NMF2S,0.001,NMF1,XE2,SCHALT,HMF1)
        IF(.not.SCHALT) GOTO 3801

c
c omit F1 feature ....................................................
c

9427    if(jf(34)) WRITE(KONSOL,11) 
11      FORMAT(1X,'*NE* HMF1 IS NOT EVALUATED BY THE FUNCTION XE2'/
     &        1X,'CORR.: NO F1 REGION, B1=3, C1=0.0')
        HMF1=0.
        F1REG=.FALSE.
c        NMF1=0.
c        C1=0.0
c
c Determine E-valley parameters if HEF was changed
c

3801     continue
         if(hef.ne.hefold) then
            width=hef-hme
            IF(ENIGHT) DEPTH=-DEPTH
            CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
            IF(.NOT.EXT) GOTO 380
            if(jf(34)) WRITE(KONSOL,650)
            WIDTH=.0
            hef=hme
            hefold=hef
            goto 9245
            endif

C
C SEARCH FOR HST [NE3(HST)=NMEs] ......................................
C

380     continue

        IF(F1REG) then
            hf1=hmf1
            xf1=nmf1
        else
            hf1=(hmf2+hef)/2.
            xf1=xe2(hf1)
        endif

        hf2=hef
        xf2=xe3_1(hf2)
        if(xf2.gt.nmes) goto 3885
            
        CALL REGFA1(hf1,HF2,XF1,XF2,0.001,NMES,XE3_1,SCHALT,HST)
        if(schalt) goto 3885

        HZ=(HST+HF1)/2.0
        D=HZ-HST
        T=D*D/(HZ-HEF-D)
        GOTO 4933

c
c assume linear interpolation between HZ and HEF ..................
c

3885    if(jf(34)) WRITE(KONSOL,100)
100     FORMAT(1X,'*NE* HST IS NOT EVALUATED BY THE FUNCTION XE3')
        HZ=(HEF+HF1)/2.
        xnehz=xe3_1(hz)
        if(jf(34)) WRITE(KONSOL,901) HZ,HEF
901     FORMAT(6X,'CORR.: LIN. APP. BETWEEN HZ=',F5.1,
     &          ' AND HEF=',F5.1)
        T=(XNEHZ-NMES)/(HZ-HEF)
        HST=-333.
        GOTO 4933

C
C LAY-functions for middle ionosphere
C

6153    IF(HMF1IN) THEN
          HMF1M=AHMF1
        ELSE
          HMF1M=165.+0.6428*XHI
        ENDIF
        HHALF = GRAT * HMF2
        HV1R = HME + WIDTH
        HV2R = HME + HDEEP
        HHMF2 = HMF2
        CALL INILAY(FNIGHT,F1REG,NMF2S,NMF1,NMES,VNER,HHMF2,HMF1M,HME,
     &                  HV1R,HV2R,HHALF,HXL,SCL,AMP,IIQU)
        IF((IIQU.EQ.1).and.(jf(34))) WRITE(KONSOL,7733)
7733   FORMAT('*NE* LAY amplitudes found with 2nd choice of HXL(1).')
        IF((IIQU.EQ.2).and.(jf(34))) WRITE(KONSOL,7722)
7722   FORMAT('*NE* LAY amplitudes could not be found.')

C
C---------- CALCULATION OF NEUTRAL TEMPERATURE PARAMETER-------
C

4933  HTA=60.0
      HEQUI=120.0
      IF(NOTEM) GOTO 240
      SEC=hourut*3600.
      CALL APFMSIS(IYEAR,MONTH,IDAY,HOUR,IAPO)
c           CALL TRETRV(SWMI)
           SWMI(9)=-1.0
           CALL TSELEC(SWMI)
      if(iapo(2).lt.0.0) then
           SWMI(9)=0.
           CALL TSELEC(SWMI)
           endif      
      CALL GTD7(IYD,SEC,HEQUI,LATI,LONGI,HOUR,F10781,F107Y,IAPO,0,
     &        D_MSIS,T_MSIS)
      TN120=T_MSIS(2)
      IF(HOUR.NE.0.0) THEN
         if(jf(18)) then
             secni=(24.-longi/15)*3600.
         else
             iyz=iyear
             idz=daynr
             call ut_lt(1,utni,0.0,longi,iyz,idz)
             secni=utni*3600.
         endif
         CALL GTD7(IYD,SECNI,HEQUI,LATI,LONGI,0.0,F10781,F107Y,IAPO,0,
     &        D_MSIS,T_MSIS)
         TN1NI=T_MSIS(2)         
      ELSE
         TN1NI=T_MSIS(2)
      ENDIF

C
C--------- CALCULATION OF ELECTRON TEMPERATURE PARAMETER--------
C

881   CONTINUE

c Te(120km) = Tn(120km)

            AHH(1)=120.
            ATE(1)=TN120

C Te-MAXIMUM based on JICAMARCA and ARECIBO data 

      HMAXD=60.*EXP(-(MLAT/22.41)**2)+210.
      HMAXN=150.
      AHH(2)=HPOL(HOUR,HMAXD,HMAXN,SAX200,SUX200,1.,1.)
      TMAXD=800.*EXP(-(MLAT/33.)**2)+1500.
      CALL GTD7(IYD,SEC,HMAXN,LATI,LONGI,HOUR,F10781,F107Y,IAPO,0,
     &        D_MSIS,T_MSIS)
      TMAXN=T_MSIS(2)
      ATE(2)=HPOL(HOUR,TMAXD,TMAXN,SAX200,SUX200,1.,1.)

c Te(300km), Te(400km) from AE-C, Te(1400km), Te(3000km) from 
c ISIS, Brace and Theis

              DIPLAT=MAGBR
              CALL TEBA(DIPLAT,HOUR,NSEASN,TEA)

              icd=0              
          if(jf(23)) then

c Te at fixed heights taken from Brace and Theis

              AHH(3)=300.
              AHH(4)=400.
              AHH(5)=600.
              AHH(6)=1400.
              AHH(7)=3000.
              hte=3000
              ATE(3)=TEA(1)
              ATE(4)=TEA(2)
              ATE(6)=TEA(3)
              ATE(7)=TEA(4)

c Te(600km) from AEROS, Spenner and Plugge (1979)

              ETT=EXP(-MLAT/11.35)
              TET=2900.-5600.*ETT/((ETT+1)**2.)
              TEN=839.+1161./(1.+EXP(-(ABSMLT-45.)/5.))
              ATE(5)=HPOL(HOUR,TET,TEN,SAX300,SUX300,1.5,1.5)
          else

c New model with solar activity effects included (Truhlik et al., 2011)
c Te at fixed heights 350, 550, 650, 1400, and 2000 km

              AHH(3)=350.
              AHH(4)=550.
              AHH(5)=850.
              AHH(6)=1400.
              AHH(7)=2000.
              hte=2500
              dimo=0.311653
c              icd=1    ! compute INVDIP
c              isa=0    ! solar activity correction off
c              ise=0    ! season correction off
              do 3491 ijk=3,7
                 call igrf_sub(lati,longi,ryear,ahh(ijk),
     &                xl,icode,dipl,babs)
                 if(xl.gt.10.) xl=10.
                 call elteik(1,1,invdip,xl,dimo,babs,dipl,
     &              xmlt,ahh(ijk),daynr,pf107,teh2,sdte)
3491             ate(ijk)=teh2
          endif

c Option to use Te = f(Ne) relation at ahh(3), ahh(4)

          IF(TENEOP) THEN
              DO 3395 I=1,2
3395              IF(TECON(I)) ATE(I+2)=TEDE(AHH(I+2),XNAR(I),-COV)
              endif

c Te corrected and Te > Tn enforced

      CALL GTD7(IYD,SEC,AHH(2),LATI,LONGI,HOUR,F10781,F107Y,IAPO,0,
     &        D_MSIS,T_MSIS)
      TNAHH2=T_MSIS(2)
      IF(ATE(2).LT.TNAHH2) ATE(2)=TNAHH2
      STTE1=(ATE(2)-ATE(1))/(AHH(2)-AHH(1))
      DO 1901 I=2,6
         CALL GTD7(IYD,SEC,AHH(I+1),LATI,LONGI,HOUR,F10781,F107Y,
     &        IAPO,0,D_MSIS,T_MSIS)
         TNAHHI=T_MSIS(2)
         IF(ATE(I+1).LT.TNAHHI) ATE(I+1)=TNAHHI
         STTE2=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
         ATE(I)=ATE(I)-(STTE2-STTE1)*DTE(I-1)*ALOG2
1901  STTE1=STTE2

c Te gradients STTE are computed for each segment

      DO 1902 I=1,6
1902     STTE(I)=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
      ATE1=ATE(1)
887   CONTINUE

C
C------------ CALCULATION OF ION TEMPERATURE PARAMETERS--------
C

c Ti(430km) during daytime from AEROS data

      XSM1=430.0
      XSM(1)=XSM1
      Z1=EXP(-0.09*MLAT)
      Z2=Z1+1.
      TID1 = 1240.0 - 1400.0 * Z1 / ( Z2 * Z2 )
      MM(2)=HPOL(HOUR,3.0,0.0,SAX300,SUX300,1.,1.)

c Ti(430km) duirng nighttime from AEROS data

      Z1=ABSMLT
      Z2=Z1*(0.47+Z1*0.024)*UMR
      Z3=COS(Z2)
      TIN1=1200.0-300.0*SIGN(1.0,Z3)*SQRT(ABS(Z3))

c Ti(430km) for specified time using HPOL

      TI1=TIN1  
      IF(TID1.GT.TIN1) TI1=HPOL(HOUR,TID1,TIN1,SAX300,SUX300,1.,1.)

c Tn < Ti < Te enforced

      TEN1=ELTE(XSM1)
      CALL GTD7(IYD,SECNI,XSM1,LATI,LONGI,0.0,F10781,F107Y,
     &        IAPO,0,D_MSIS,T_MSIS)
      TNN1=T_MSIS(2)
      IF(TEN1.LT.TNN1) TEN1=TNN1
      IF(TI1.GT.TEN1) TI1=TEN1
      IF(TI1.LT.TNN1) TI1=TNN1

c Tangent on Tn profile determines HS

      HS=200.
      CALL GTD7(IYD,SEC,HS,LATI,LONGI,HOUR,F10781,F107Y,
     &        IAPO,0,D_MSIS,T_MSIS)
      TNHS=T_MSIS(2)
      MM(1)=(TI1-TNHS)/(XSM1-HS)
      MXSM=2

c XTETI is altitude where Te=Ti

2391    XTTS=500.
        X=500.
2390    X=X+XTTS
        IF(X.GE.AHH(7)) GOTO 240
        TEX=ELTE(X)
        TIX=TI(X)
        IF(TIX.LT.TEX) GOTO 2390
        X=X-XTTS
        XTTS=XTTS/10.
        IF(XTTS.GT.0.1) GOTO 2390
        XTETI=X+XTTS*5.

c Ti=Te above XTETI 

        MXSM=3
        MM(3)=STTE(6)
        XSM(2)=XTETI
        IF(XTETI.GT.AHH(6)) GOTO 240
        MXSM=4
        MM(3)=STTE(5)
        MM(4)=STTE(6)
        XSM(3)=AHH(6)
        IF(XTETI.GT.AHH(5)) GOTO 240
        MXSM=5
        DTI(1)=5.
        DTI(2)=5.
        MM(3)=STTE(4)
        MM(4)=STTE(5)
        MM(5)=STTE(6)
        XSM(3)=AHH(5)
        XSM(4)=AHH(6)

C
C CALCULATION OF ION DENSITY PARAMETER..................
C

240   IF(NOION) GOTO 141
      HNIA=75.
      if(DY) HNIA=80.
      HNIE=2000.
C
C CALCULATION FOR THE REQUIRED HEIGHT RANGE.......................
C In the absence of an F1 layer hmf1=hz since hmf1 is used in XE
C

141     xhmf1=hmf1
        IF(hmf1.le.0.0) HMF1=HZ

        height=heibeg
        kk=1

300   IF(NODEN) GOTO 330
      IF((HEIGHT.GT.HNEE).OR.(HEIGHT.LT.HNEA)) GOTO 330
      IF(LAYVER) THEN
          ELEDE=-9.
          IF(IIQU.LT.2) ELEDE=
     &            XEN(HEIGHT,HMF2,NMF2S,HME,4,HXL,SCL,AMP)
          outf(1,kk)=elede
          goto 330
          endif
      ELEDE=XE_1(HEIGHT)

c
c electron density in m-3 in outf(1,*)
c

      OUTF(1,kk)=ELEDE

c
c plasma temperatures
c

330   IF(NOTEM) GOTO 7108
      IF((HEIGHT.GT.HTE).OR.(HEIGHT.LT.HTA)) GOTO 7108
      CALL GTD7(IYD,SEC,HEIGHT,LATI,LONGI,HOUR,F10781,F107Y,
     &        IAPO,0,D_MSIS,T_MSIS)
      TNH=T_MSIS(2)
      TIH=TNH
      if(HEIGHT.GT.HS) then
      	TIH=TI(HEIGHT)
      	if(TIH.lt.TNH) TIH=TNH
      	endif
      TEH=TNH
      if(HEIGHT.GT.HEQUI) then 
         TEH=ELTE(HEIGHT)
      	 if(TEH.lt.TIH) TEH=TIH
      	 endif

        OUTF(2,kk)=TNH
        OUTF(3,kk)=TIH
        OUTF(4,kk)=TEH

c
c ion composition
c
7108  IF(NOION) GOTO 7118
      IF((HEIGHT.GT.HNIE).OR.(HEIGHT.LT.HNIA)) GOTO 7118

            ROX=-1.
            RHX=-1.
            RNX=-1.
            RHEX=-1.
            RNOX=-1.
            RO2X=-1.
            RCLUST=-1.
      if(DY) then
        if (height.gt.300.) then
c Triskova-Truhlik-Smilauer-2003 model
       		call igrf_sub(lati,longi,ryear,height,
     &  	    xl,icode,dipl,babs)
        	if(xl.gt.10.) xl=10.
        	dimo=0.311653
			call CALION(1,xinvdip,xl,dimo,babs,dipl,xmlt,
     &  	   height,daynr,f107d,xic_O,xic_H,xic_He,xic_N)
       		rox=xic_O*100.
        	rhx=xic_H*100.
        	rnx=xic_N*100.
        	rhex=xic_He*100.
        	rnox=0.
        	ro2x=0.
        else
c Richards-Bilitza-Voglozin-2010 IDC model
            CALL GTD7(IYD,SEC,height,lati,longi,HOUR,f10781,f107y,
     &        IAPO,48,D_MSIS,T_MSIS)
			XN4S = 0.5 * D_MSIS(8)
			EDENS=ELEDE/1.e6
			jprint=1
			if(jf(38)) jprint=0
            CALL CHEMION(jprint,height,F107Y,F10781,TEH,TIH,TNH,
     &       	D_MSIS(2),D_MSIS(4),D_MSIS(3),D_MSIS(1),-1.0,XN4S,
     &       	EDENS,-1.0,xhi,ro,ro2,rno,rn2,rn,Den_NO,Den_N2D,INEWT)                              

c			if(INEWT.gt.0) then
				sumion = edens/100.
        		rox=ro/sumion
        	    rhx=0.
        	    rhex=0.
        		rnx=rn/sumion
        		rnox=rno/sumion
        		ro2x=ro2/sumion
c        		endif
        endif
      else
c Danilov-Smirnova-1995 model and Danilov-Yaichnikov-1985 model (upper)
            call iondani(iday,iseamon,height,xhi,
     &            lati,f107365,dion)
            ROX=DION(1)
            RHX=DION(2)
            RNX=DION(3)
            RHEX=DION(4)
            RNOX=DION(5)
            RO2X=DION(6)
            RCLUST=DION(7)
      endif

c
c ion densities are given in percent of total electron density;
c

      if(jf(22)) then 
            xnorm=1
      else
            xnorm=elede/100.
      endif
      OUTF(5,kk)=ROX*xnorm
      OUTF(6,kk)=RHX*xnorm
      OUTF(7,kk)=RHEX*xnorm
      OUTF(8,kk)=RO2X*xnorm
      OUTF(9,kk)=RNOX*xnorm
      OUTF(10,kk)=RCLUST*xnorm
      OUTF(11,kk)=RNX*xnorm

c
c D region special: Friedrich&Torkar model in outf(13,*)
c

7118    if(.not.dreg.and.height.le.140.) then
            outf(1,kk)=-1.
            call F00(HEIGHT,LATI,DAYNR,XHI,F107D,EDENS,IERROR)
            if(ierror.eq.0.or.ierror.eq.2) outf(1,kk)=edens
            endif


        height=height+heistp
        kk=kk+1
        if(kk.le.numhei) goto 300

C
C END OF PARAMETER COMPUTATION LOOP 
C

c
c D region special: densities for 11 heights (60,65,70,..,110km)
c outf(14,1:11)=IRI-07, outf(14,12:22)=FIRI, 
c outf(14,23:33)= Danilov et al.(1995) with SW=0,WA=0 
c outf(14,34:44)= with SW=0.5,WA=0, 
c outf(14,45:55)= with SW=1,WA=0,  
c outf(14,56:66)= with SW=0,WA=0.5, 
c outf(14,67:77)= with SW=0,WA=1,  
c

      if(.not.dreg) then
            do ii=1,11
                  Htemp=55+ii*5  
                  outf(14,ii)=-1.     
                  if(Htemp.ge.65.) outf(14,ii)=XE6(Htemp)     
                  outf(14,11+ii)=-1.
                  call F00(Htemp,LATI,DAYNR,XHI,F107D,EDENS,IERROR)
                  if(ierror.eq.0.or.ierror.eq.2) outf(14,11+ii)=edens
                  outf(14,22+ii)=ddens(1,ii)      
                  outf(14,33+ii)=ddens(2,ii)      
                  outf(14,44+ii)=ddens(3,ii)      
                  outf(14,55+ii)=ddens(4,ii)      
                  outf(14,66+ii)=ddens(5,ii)      
                  enddo
            endif

c
c equatorial vertical ion drift
c

      drift=-1.
      if(jf(21).and.abs(magbr).lt.25.0) then
            param(1)=daynr
            param(2)=f107d
            call vdrift(hour,longi,param,drift)
            endif
c
c spread-F occurrence probability
c
      spreadf=-1.
      if(.not.jf(28)) goto 1937
      if(hour.gt.7.25.and.hour.lt.17.75) goto 1937
      if(abs(lati).gt.25.0) goto 1937
			spfhour=hour
			daynr1=daynr
			if(hour.lt.12.0) then
				spfhour=hour+24.0
				daynr1=daynr-1
				if(daynr1.lt.1) daynr1=idayy
            	endif
            call spreadf_brazil(daynr,idayy,f107d,lati,osfbr)
			ispf=int((spfhour-17.75)/0.5)+1
 			if(ispf.gt.0.and.ispf.lt.26) spreadf=osfbr(ispf)
1937   continue
C
C ADDITIONAL PARAMETER FIELD OARR
C

        IF(NODEN) GOTO 6192
      OARR(1)=NMF2S
      OARR(2)=HMF2
      if(f1reg) OARR(3)=NMF1
      if(f1reg) OARR(4)=XHMF1
      OARR(5)=NMES
      OARR(6)=HME
      OARR(7)=NMD
      OARR(8)=HMD
      OARR(9)=HHALF
      OARR(10)=B0
      OARR(11)=VNER
      OARR(12)=HEF
6192    IF(NOTEM) GOTO 6092
      OARR(13)=ATE(2)
      OARR(14)=AHH(2)
      OARR(15)=ATE(3)
      OARR(16)=ATE(4)
      OARR(17)=ATE(5)
      OARR(18)=ATE(6)
      OARR(19)=ATE(7)
      OARR(20)=ATE(1)
      OARR(21)=TI1
      OARR(22)=XTETI
6092  OARR(23)=XHI
      OARR(24)=SUNDEC
      OARR(25)=DIP
      OARR(26)=MAGBR
      OARR(27)=MODIP
      OARR(28)=LATI		
      OARR(29)=SAX200
      OARR(30)=SUX200
      OARR(31)=SEASON
      OARR(32)=LONGI		
      OARR(33)=rssn
      OARR(34)=COV
      OARR(35)=B1
      OARR(36)=xm3000
C OARR(37) used for TEC and 38 for TEC-top
      OARR(39)=gind
      OARR(40)=f1pb
      OARR(41)=f107d
      OARR(42)=c1
      OARR(43)=daynr
      OARR(44)=drift
      OARR(45)=stormcorr
      OARR(46)=f10781
      OARR(47)=estormcor
      OARR(48)=spreadf
      OARR(49)=MLAT
      OARR(50)=MLONG
      OARR(51)=indap(13)   ! ap for current UT
      OARR(52)=IAP_daily   ! daily ap
      OARR(53)=invdip
      OARR(54)=XMLT
      OARR(55)=cgm_lat
      OARR(56)=cgm_lon
      OARR(57)=cgm_mlt
      OARR(58)=cgmlat   ! CGM latitude of equatorward boundary
c include only every second auroral boundary point (MLT=0,1,2..23)
      jjj=58
      do iii=1,47,2
         jjj=jjj+1 
         oarr(jjj)=ab_mlat(iii)
         enddo
      OARR(83)=xkp
      OARR(84)=dec

3330  CONTINUE

       icalls=icalls+1

      RETURN
      END
c
c
        subroutine iri_web(jmag,jf,alati,along,iyyyy,mmdd,iut,dhour,
     &          height,h_tec_max,ivar,vbeg,vend,vstp,a,b)
c-----------------------------------------------------------------------        
c changes:
c       11/16/99 jf(30) instead of jf(17)
c       10/31/08 outf, a, b (100 -> 500)
c
c-----------------------------------------------------------------------        
c input:   jmag,alati,along,iyyyy,mmdd,dhour  see IRI_SUB
c          height  height in km
c          h_tec_max  =0 no TEC otherwise upper boundary for integral
c          iut     =1 for UT       =0 for LT
c          ivar    =1      altitude
c                  =2,3    latitude,longitude
c                  =4,5,6  year,month,day
c                  =7      day of year
c                  =8      hour (UT or LT)
c          vbeg,vend,vstp  variable range (begin,end,step)
c output:  a       similar to outf in IRI_SUB
c          b       similar to oarr in IRI_SUB
c
c          numstp  number of steps; maximal 1000
c-----------------------------------------------------------------------        
        dimension   outf(20,1000),oar(100),oarr(100),a(20,1000)
        dimension   xvar(8),b(100,1000)
        logical     jf(50)

		nummax=1000
        numstp=int((vend-vbeg)/vstp)+1
        if(numstp.gt.nummax) numstp=nummax

        do 6249 i=1,100
6249          oar(i)=b(i,1) 

        if(ivar.eq.1) then
            do 1249 i=1,100
1249            oarr(i)=oar(i) 
            xhour=dhour+iut*25.
            call IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,XHOUR,
     &                  VBEG,VEND,VSTP,a,OARR)
            if(h_tec_max.gt.50.) then 
                call iri_tec (50.,h_tec_max,2,tec,tect,tecb)
                oarr(37)=tec
                oarr(38)=tect
                endif
            do 1111 i=1,100
1111            b(i,1)=oarr(i)
            return
            endif

        if(height.le.0.0) height=100
        xvar(2)=alati
        xvar(3)=along
        xvar(4)=iyyyy
        xvar(5)=mmdd/100
        xvar(6)=mmdd-xvar(5)*100
        xvar(7)=abs(mmdd*1.)
        xvar(8)=dhour

        xvar(ivar)=vbeg

        alati=xvar(2)
        along=xvar(3)
        iyyyy=int(xvar(4))
        if(ivar.eq.7) then
                mmdd=-int(vbeg)
        else
                mmdd=int(xvar(5)*100+xvar(6))
        endif
        dhour=xvar(8)+iut*25.

        do 1 i=1,numstp
                do 1349 iii=1,100
1349                    oarr(iii)=b(iii,i) 
                call IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &                  height,height,1.,OUTF,OARR)
                if(h_tec_max.gt.50.) then
                        call iri_tec (50.,h_tec_max,2,tec,tect,tecb)
                        oarr(37)=tec
                        oarr(38)=tect
                        endif
                do 2 ii=1,20
2                       a(ii,i)=outf(ii,1)
                do 2222 ii=1,100
2222                    b(ii,i)=oarr(ii)
                xvar(ivar)=xvar(ivar)+vstp

                alati=xvar(2)
                along=xvar(3)
                iyyyy=int(xvar(4))
                if(ivar.eq.7) then
                        mmdd=-xvar(7)
                else
                        mmdd=int(xvar(5)*100+xvar(6))
                endif
                dhour=xvar(8)+iut*25.
1       continue

        return
        end
