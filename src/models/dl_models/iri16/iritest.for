c iritest.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
c
c test program for the iri_web subroutine
c
c-version-mm/dd/yy ----------corrections--------------------------
c 2000.01 05/07/01 initial version
c 2000.02 07/11/01 line 210: do i=1,100 instead of i=2,100 (K. Tokar)
c 2000.03 28/12/01 output oar(39) for IG12 (R. Conkright, NGDC, NOAA)
c 2000.04 28/10/02 replace TAB/6 blanks, enforce 72/line (D. Simpson)
c 2000.05 02/06/03 Ne(Te) only 300,400; foF1 and hmF1 output corr.
c 2000.06 01/19/05 (.not.jf(20)) instead of (..jf(2)) (G. Schiralli)
c 2005.01 05/06/06 included spread-F (jf(28)) and topside (jf(29)) options
C 2007.00 05/18/07 Release of IRI-2007
c 2007.02 10/31/08 outf(100) -> outf(500), numhei=numstp=500
c 2007.03 02/12/09 added new D-region option (h=-3)
c 2007.11 04/19/10 correct TEC for normal output  [Shunrong Zhang] 
c
C 2012.00 10/05/11 IRI-2012: bottomside B0 B1 model (SHAMDB0D, SHAB1D),
C 2012.00 10/05/11    bottomside Ni model (iriflip.for), auroral foE
C 2012.00 10/05/11    storm model (storme_ap), Te with PF10.7 (elteik),
C 2012.00 10/05/11    oval kp model (auroral_boundary), IGRF-11(igrf.for), 
C 2012.00 10/05/11    NRLMSIS00 (cira.for), CGM coordinates, F10.7 daily
C 2012.00 10/05/11    81-day 365-day indices (apf107.dat), ap->kp (ckp),
C 2012.00 10/05/11    array size change jf(50) outf(20,1000), oarr(100).
C 2012.00 03/21/12    PIKTAB=4 output for D-region
C 2012.01 09/16/12 Corrected UT output (hour-25)
C 2012.03 09/18/14 input JF(18): FIELDG not UT_LT          (A. Mazzella)
C 2012.03 09/18/14 rzino, igino logical not real; and more (A. Mazzella)
C 2012.03 09/18/14 jf defaults and special for piktab=3    (A. Mazzella)
C 2012.03 09/23/14 jf(22): output option Ni in [m-3]/1.e9
C 2012.03 09/23/14 jf(36:38) include in input choices
C 2012.03 09/23/14 added output options for parameters oar(59:86) 
C 2012.05 04/27/15 delete double line CHARACTER*9  pname(6)
C 2015.01 06/30/15 scid=1.0E-9                              (A. Charisi)
C 2015.01 06/30/15 jino-outf(9)&jio2-outf(8) write.jio2,jino(A. Charisi)
C 2015.01 07/12/15 output of ion density/composition now with 3 digits
C 2015.01 07/12/15 call read_ig_rz readapf107
C 2015.02 08/13/15 delete COMMON/CONST2
C 2015.02 08/13/15 ursifo=jf(5) just before IRI_SUB call
C 2016.01 06/01/16 User specified B0 when jf(43)=false
C 2016.02 08/22/16 Allow user input options for F10.7D and F10.7_81
C
      INTEGER           pad1(6),jdprof(77),piktab
      DIMENSION         outf(20,1000),oar(100,1000),jfi(6)
      LOGICAL		    jf(50),rzino,igino
      CHARACTER*2       timev(2)
      CHARACTER*3       uni(86),sopt,seopt
      CHARACTER*4       IMZ(8),MAP,xtex,coorv(2)
      CHARACTER*5       ITEXT(8),tsopt
      CHARACTER*6       pna(86)
      CHARACTER*7       popt
      CHARACTER*8       bopt,topt
      CHARACTER*9       pname(7)
      CHARACTER*10      dopt,hopt
      CHARACTER*11      iopt,rzopt,igopt,f8opt,fdopt
      CHARACTER*16      f1opt

      DATA  IMZ  /' km ','GEOD','GEOD','yyyy',' mm ',' dd ','YEAR',
     &      'L.T.'/, ITEXT/'  H  ',' LATI',
     &      ' LONG',' YEAR','MONTH',' DAY ','DAYOF',' HOUR'/

      data pna/'NmF2','hmF2','NmF1','hmF1','NmE','hmE','NmD','hmD',
     &  'h05','B0','NVmin','hVtop','Tpeak','hTpek','T300','T400','T600',
     &  'T1400','T3000','T120','Ti450','hTeTi','sza','sundec','dip',
     &  'diplat','modip','Lati','Srise','Sset','season','Longi',
     &  'Rz12','cov','B1','M3000','TEC','TECtop','IG12','F1prb','F107d',
     &  'C1','daynr','vdrft','foF2r','F10781','foEr','sprd_F','MLAT',
     &  'MLON','Ap_t','Ap_d','invdip','MLTinv','CGMlat','CGMlon',
     &  'CGMmlt','CGM_AB','CGMm0','CGMm1','CGMm2','CGMm3','CGMm4',
     &  'CGMm5','CGMm6','CGMm7','CGMm8','CGMm9','CGMm10','CGMm11',
     &  'CGMm12','CGMm13','CGMm14','CGMm15','CGMm16','CGMm17','CGMm18',
     &  'CGMm19','CGMm20','CGMm21','CGMm22','CGMm23','kp_t','dec','L',
     &  'DIMO'/
      data uni/'m-3','km','m-3','km','m-3','km','m-3','km','km','km',
     &   'm-3','km','K','km',7*'K','km',6*'deg',2*'h',' ','deg',4*' ',
     &   'm-2','%',5*' ','m/s',4*' ',2*'deg',2*' ','deg','h',2*'deg',
     &   'h',25*'deg',' ','deg',' ','Gau'/,
     &   timev/'LT','UT'/,coorv/'geog','geom'/

      data jfi/8,9,13,14,15,16/

c	  COMMON/const2/icalls,montho,nmono,iyearo,idaynro,ursifo,rzino,
c     &	            igino,ut0
c
c		icalls=0
c		montho=-1
c		nmono=-1
c		iyearo=-1
c		idaynro=-1
c		rzino=.true.
c		igino=.true.
c		ut0=-1

		call read_ig_rz
        call readapf107
        
        nummax=1000
        
        do 6249 i=1,100
6249    oar(i,1)=-1.0

c user input of IRI input parameters
c
1       print *,'jmag(=0/1,geog/geom),lati/deg,long/deg'
        read(5,*) jm,xlat,xlon
        print *,'year(yyyy),mmdd(or -ddd),iut(=0/1,LT/UT),hour'
        read(5,*) iy,imd,iut,hour
        print *,'height/km'
        read(5,*) hx

        print *,'output-option'
        print *,'(enter 0 for standard table of IRI parameters)'
        print *,'(enter 1 for list of peak heights and densities)'
        print *,'(enter 2 for plasma frequencies, B0, M3000, ',
     &                        'valley, width and depth,)'
        print *,'(enter 3 for 6 parameters of your choice)'
        print *,'(enter 4 for D-region models at 60,65,..,110 km)'
        print *,'(enter 5 special test output)'
        read(5,*) piktab

        print *,'upper height [km] for TEC integration (0 for no TEC)'
        read(5,*) htec_max

        print *,'variable? (1/2/../8 for height/lat/long/year/month/',
     &                        'day/day of year/hour)'
        read(5,*) ivar
        print *,'begin, end, and stepsize for the selected variable'
        read(5,*) vbeg,vend,vstp

        print *,'Options: t(rue) or f(alse)'
        print *,'Enter 0 to use standard or 1 to enter your own'
        read(5,*) jchoice
          do i=1,50 
                jf(i)=.true.
                enddo
          if(piktab.eq.4) jf(24)=.false.
        if(jchoice.eq.0) then
c defaults for jf(1:50)
c          jf(1)=.false.      ! f=no electron densities (t) 
c          jf(2)=.false.      ! f=no temperatures (t)
c          jf(3)=.false.      ! f=no ion composition (t)
          jf(4)=.false.      ! t=B0table f=other models (f)
          jf(5)=.false.      ! t=CCIR  f=URSI foF2 model (f)
          jf(6)=.false.      ! t=DS95+DY85   f=RBV10+TBT15 (f)
c          jf(7)=.false.      ! t=tops f10.7<188 f=unlimited (t)
          jf(21)=.false.     ! f=ion drift not computed (f)
c          jf(22)=.false.     ! ion densities in m-3 (t)
          jf(23)=.false.     ! t=AEROS/ISIS f=TTS Te with PF10.7 (f)
c          jf(24)=.false.     ! t=D-reg-IRI-1990 f=FT-2001 (t)
c          jf(26)=.false.	  ! f=STORM model turned off (t)
          jf(28)=.false.	  ! f=spread-F not computed (f)
          jf(29)=.false.     ! t=old  f=New Topside options (f)
          jf(30)=.false.     ! t=corr f=NeQuick topside (f)
c          jf(31)=.false.     ! t=B0ABT f=Gulyaeva (t)
c          jf(33)=.false. 	  ! f=auroral boundary off (f)
          jf(33)=.true. 	  ! f=auroral boundary off (f)
c          jf(34)=.false. 	  ! t=messages on f= off (t)
          jf(35)=.false. 	  ! f=auroral E-storm model off (f)
c          jf(36)=.false. 	  ! t=hmF2 w/out foF2_storm f=with (t)
c          jf(37)=.false. 	  ! t=topside w/out foF2_storm f=with (t)
c          jf(38)=.false. 	  ! t=WRITEs off in IRIFLIP f=on (t)
          jf(39)=.false. 	  ! t=M3000F2 model f=new hmF2 models (f)
c          jf(40)=.false. 	  ! t=AMTB-model, f=Shubin-COSMIC model (t) 
c          jf(41)=.false. 	  ! t:COV=F10.7_386 f:COV=f(IG12) (t) 
c          jf(42)=.false. 	  ! t/f=Te w/o PF10.7 dependance (t)
c          jf(43)=.false. 	  ! t= B0 model f= B0 user input (t)
c          jf(44)=.false. 	  ! t= B1 model f= B1 user input (t)
        else
          print *,'Compute Ne, T, Ni? (enter: t,t,t  if you want all)'
          read(5,*) jf(1),jf(2),jf(3)
        if(jf(1)) then 
              print *,'LAY version: t=standard ver., f=LAY version.',
     &              ' {standard:t}'
              read(5,*) jf(11)
              print *,'Ne Topside: t=IRI-2001/h0.5, f=new options {f}'
              read(5,*) jf(29)
              print *,'Ne Topside: t=IRI01_corrt, f=NeQuick/h0.5 {f}'
              read(5,*) jf(30)
              print *,'Ne Topside: t=F10.7<188, f=unlimited {t}'
              read(5,*) jf(7)
              print *,'Ne topside: t=w/o foF2 storm model, f=with {t}'
              read(5,*) jf(37)
              print *,'foF2 model: t=CCIR, f=URSI-88 {f}'
              read(5,*) jf(5)
              print *,'foF2: t=with storm model, f=without {t}'
              read(5,*) jf(26)
              print *,'hmF2: t=f(M3000F2), f=new models {f}'
              read(5,*) jf(39)              
              print *,'hmF2: t=AMTB-model, f=Shubin-COSMIC model {t}'
              read(5,*) jf(40)              
              if(jf(39)) then
                 print *,'hmF2: t=w/o foF2 storm model, f=with {t}'
                 read(5,*) jf(36)
                 endif              
              print *,'F2 peak density or foF2: t=model, ',
     &              'f=user input {t}'
              read(5,*) jf(8)
              print *,'F2 peak height or M3000F2: t=model, ',
     &              'f=user input {t}'
              read(5,*) jf(9)
              print *,'Auroral boundary model: t=on, f=off {f}'
              read(5,*) jf(33)
              print *,'Bottomside thickness B0: t=Bil-2000, ',
     &            'f=other options {f}.'
              read(5,*) jf(4)
              print *,'Bottomside thickness B0: t=ABT-2009, ',
     &            'f= Gul-1987 {t}.'     
              read(5,*) jf(31)
              print *,'Bottomside thickness B0: t=model, ',
     &            'f=user input {t}.'     
              read(5,*) jf(43)
              print *,'F1 peak density or foF1: t=model, ',
     &            'f=user input {t}'
              read(5,*) jf(13)
            if(.not.jf(11)) then
              print *,'F1 peak height: t=model, f=user input {t}'
              read(5,*) jf(14)
            endif
              print *,'F1: t=with probability model, f=without   {t}'
              read(5,*) jf(19)
              print *,'F1: t=standard probability, f=with L ',
     &            'condition {t}'
              read(5,*) jf(20)
              print *,'E peak density or foE: t=model, f=user input {t}'
              read(5,*) jf(15)
              print *,'E peak height: t=model, f=user input {t}'
              read(5,*) jf(16)
              print *,'E peak auroral storm model: t=on, f=off {f}'
              read(5,*) jf(35)
              print *,'D: t=IRI-1990, f= FT-2001 {t}'
              read(5,*) jf(24)
          endif
        if(jf(2)) then
              print *,'Te(Ne) model: t=not used, f=correlation is',
     &          ' used. {t}'
              read(5,*) jf(10)
              print *,'Te: t=Bil-1985, f=TBT-2012 {f}'
              read(5,*) jf(23)
              print *,'Te: t=TBT-2012 with PF107 dep., f=w/o {t}'
              read(5,*) jf(42)
          endif
        if(jf(3)) then
              print *,'Ion comp. model: t=DS95/DY85, f=RBV10/TBT15 {f}' 
              read(5,*) jf(6)
              if(.not.jf(6)) then
              	print *,'IRIFLIP: t=messages off, f=on {t}' 
              	read(5,*) jf(38)
              	endif              
              print *,'Ni: t=ion composition in %, f=ion densities',
     &             'in cm-3 {t}'
              read(5,*) jf(22)
              endif
           print *,'Equat. Vert. Ion Drift: t=computed, ',
     &            'f=not computed {t}'
              read(5,*) jf(21)
           print *,'Spread-F probability: t=computed, ',
     &            'f=not computed {t}'
              read(5,*) jf(28)
           print *,'COV: t: COV=F10.7_365, f: COV=func(IG12).  {t}'
              read(5,*) jf(41)
           print *,'Sunspot index: t=from file, f=user input.  {t}'
              read(5,*) jf(17)
           print *,'Ionospheric index: t=from file, f=user input. {t}'
              read(5,*) jf(27)
           print *,'F10.7D Index: t=from file, f=user input {t}'
              read(5,*) jf(25)
           print *,'F10.7_81 Index: t=from file, f=user input {t}'
              read(5,*) jf(32)
           print *,'dip, magbr, modip: t=IGRF, f=old FIELDG using ',
     &             'POGO68/10 for 1973 {t}'
              read(5,*) jf(18)
           print *,'Messages on (t) off (f) {t}'
              read(5,*) jf(34)
           print *,'Message output unit: t=(UNIT=6), f=(UNIT=11). {t}'
              read(5,*) jf(12)
       endif

c       if(piktab.gt.3) jf(24)=.false.
c option to enter six additional parameters 
c

      if(PIKTAB.eq.3) then
        print *,'6 Parameters of your choice (number:1-86)'
        print *,(pna(j),j=1,10)
        print *,(pna(j),j=11,20)
        print *,(pna(j),j=21,30)
        print *,(pna(j),j=31,40)
        print *,(pna(j),j=41,50)
        print *,(pna(j),j=51,60)
        print *,(pna(j),j=61,70)
        print *,(pna(j),j=71,80)
        print *,(pna(j),j=81,86)
        print *,'or 0,0,0,0,0,0 for default:'
        print *,'      spread-F probability [48]'
        print *,'      equatorial vertical ion drift [44]'
        print *,'      foF2_storm/foF2_quiet [45]'
        print *,'      foE_storm/foE_quiet [47]' 
        print *,'      eqward auroral boundy CGM-Lat [58]'
c        print *,'      solar zenith angle [23]'
c        print *,'      modified dip latitude [27]' 
        print *,'      Ap for current UT [51]' 
        read(5,*) (pad1(j),j=1,6)
c change defaults for computation of specific parameters
        	jf(21)=.true.  ! spread-F prob. computed
        	jf(28)=.true.  ! vertical ion drift computed
        	jf(33)=.true.  ! auroral boundary computed
        	jf(35)=.true.  ! foE_storm computed
        if(pad1(1).eq.0) then
        	pad1(1)=48     ! spread-F probability
            pad1(2)=44     ! equatorial vertical ion drift
            pad1(3)=45     ! fof2_storm/foF2_quiet
            pad1(4)=47     ! foE_storm/foE_quiet
            pad1(5)=58     ! CGM_lat auroral boundary
            pad1(6)=51     ! ap for current UT
            endif
      endif
       
c option to enter measured values for NmF2, hmF2, NmF1, hmF1, NmE, hmE,
c B0, N(300), N(400), N(600) if available; 
c
          print *,' '
          print *,' '
          print *,' '
          numstp=int((vend-vbeg)/vstp)+1            
              if(ivar.eq.1) numstp=1
       if(jf(1)) then
         if(.not.jf(8).or..not.jf(9).or..not.jf(13).or..not.jf(14).or.
     &      .not.jf(15).or..not.jf(16).or..not.jf(43)) then
            var=vbeg
            i=1
2234        if(.not.jf(8)) then
              jf(26)=.false.    ! storm model off, if user input
              print *,'foF2/Mhz or NmF2/m-3 for ',itext(ivar),
     &             '=',var
              read(5,*) oar(1,i)
              pname(1)='foF2/MHz'
              if(oar(1,i).gt.30.) pname(1)='NmF2/m-3'
              endif
            if(.not.jf(9)) then
              print *,'hmF2/km or M3000F2 for ',itext(ivar),
     &              '=',var
              read(5,*) oar(2,i)
              pname(2)='M(3000)F2'
              if(oar(2,i).gt.50.) pname(2)='hmF2/km'
              endif
            if(.not.jf(13)) then
              print *,'foF1/MHz or NmF1/m-3 for ',itext(ivar),
     &               '=',var
              read(5,*) oar(3,i)
              pname(3)='foF1/MHz'
              if(oar(3,i).gt.30.) pname(3)='NmF1/m-3'
              endif
            if(.not.jf(14)) then
              print *,'hmF1/km for ',itext(ivar),'=',var
              read(5,*) oar(4,i)
              pname(4)='hmF1/km'
              endif
            if(.not.jf(15)) then
              print *,'foE/MHz or NmE/m-3 for ',itext(ivar),
     &                '=',var
              read(5,*) oar(5,i)
              pname(5)='foE/MHz'
              if(oar(5,i).gt.30.) pname(5)='NmE/m-3'
              endif
            if(.not.jf(16)) then
              print *,'hmE/km for ',itext(ivar),'=',var
              read(5,*) oar(6,i)
              pname(6)='hmE/km'
              endif
            if(.not.jf(43)) then
              print *,'B0/km for ',itext(ivar),'=',var
              read(5,*) oar(10,i)
              pname(7)='B0/km '
              endif
            if(.not.jf(44)) then
              print *,'B1 for ',itext(ivar),'=',var
              read(5,*) oar(87,i)
              pname(7)='B1    '
              endif
            i=i+1
            var=var+vstp
            if(ivar.gt.1.and.var.le.vend) goto 2234
         endif
       endif

c option to enter Ne for Te-Ne relationship
c
        if(jf(2).and..not.jf(10)) then
          var=vbeg
          do 1235 i=1,numstp 
            print *,'Ne(300km),Ne(400km)/m-3',
     &         ' for ',itext(ivar),'=',var,' [-1 if not]'
            read(5,*) oar(15,i),oar(16,i)
1235        var=var+vstp
          endif

c option to enter F107D and/or PF107 
c
            if(.not.jf(25)) then
                        print *,'User input for F107D:'
                        read(5,*) f107d
                        do i=1,100
                                    oar(41,i)=f107d
                                    enddo
                        endif

            if(.not.jf(32)) then
                        print *,'User input for PF107:'
                        read(5,*) pf107d
                        do i=1,100
                                    oar(46,i)=pf107d
                                    enddo
                        endif

c option to enter Rz12 and/or IG12
c
            if(.not.jf(17)) then
                        print *,'User input for Rz12'
                        read(5,*) oar(33,1)
                        do i=2,100
                                    oar(33,i)=oar(33,1)
                                    enddo
                        endif

            if(.not.jf(27)) then
                        print *,'User input for IG12'
                        read(5,*) oar(39,1)
                        do i=2,100
                                    oar(39,i)=oar(39,1)
                                    enddo
                        endif

            if(.not.jf(25)) then
                        print *,'User input for F10.7D'
                        read(5,*) oar(41,1)
                        do i=2,100
                                    oar(41,i)=oar(41,1)
                                    enddo
                        endif

            if(.not.jf(32)) then
                        print *,'User input for F10.7_81d'
                        read(5,*) oar(46,1)
                        do i=2,100
                                    oar(46,i)=oar(46,1)
                                    enddo
                        endif

c end of user input
c

        num1=(vend-vbeg)/vstp+1
        numstp=iabs(num1)
        if(numstp.GT.nummax) numstp=nummax

        if(jf(29)) then
             popt='IRI2001'
        else
             if(jf(30)) then
                   popt='IRIcorr'
             else
                   popt='NeQuick'
             endif
        endif
        map='URSI'
        if(jf(5)) map='CCIR'

        if(jf(39)) then
             hopt='CCIR-M3000'
        else
             if(jf(40)) then
                   hopt='AMTB-2013'
             else
                   hopt='Shubin2015'
             endif
        endif

        if(jf(4)) then
             bopt='BIl-2000'
        else
             if(jf(31)) then
                   bopt='ABT-2009'
             else
                   bopt='Gul-1987'
             endif
        endif

        iopt='RBV10+TBT15'
        if(jf(6)) iopt='DS95 + DY85'

        dopt='FT01+DRS95'
        if(jf(24)) dopt='IRI-1990'

        sopt='off'
        if(jf(26)) sopt='on '
        seopt='off'
        if(jf(35)) seopt='on '

        topt='TBT-2012'
        if(jf(23)) topt='BIl-1985'
        tsopt=' with'
        if(jf(23)) tsopt='  w/o'

        if(jf(19)) then
              f1opt='Scotto-97 no L'
              if(.not.jf(20)) f1opt='Scotto-97 with L'
        else
              f1opt='IRI-95'
        endif

        rzopt=' user input'
        if(jf(17)) rzopt=' '
        igopt=' user input'
        if(jf(27)) igopt=' '
        fdopt=' user input'
        if(jf(25)) fdopt=' '
        f8opt=' user input'
        if(jf(32)) f8opt=' '
        
        hxx=hx
        jmag=jm
        mmdd=imd

c calling IRI subroutine
c 
        phour=hour
        call iri_web(jmag,jf,xlat,xlon,iy,mmdd,iut,hour,
     &          hxx,htec_max,ivar,vbeg,vend,vstp,outf,oar)

c preparation of results page
c
        write(7,3991) iy,mmdd,phour,timev(iut+1),
     &        coorv(jmag+1),xlat,xlon,hxx
        if(jf(1)) then
       		write(7,3314) popt
            if(jf(8)) then
            	write(7,301) map
                write(7,3291) sopt
                endif
            if(jf(9)) write(7,303) hopt
            write(7,309) bopt
            write(7,3295) f1opt
            write(7,3299) seopt
            write(7,3081) dopt
            numi=numstp
            if(ivar.eq.1) numi=1
            do j=1,6
                ij=jfi(j)
                if(.not.jf(ij)) then
                    write(7,302) pname(j)
                    write(7,402) (oar(j,i),i=1,numi)
                    endif
                enddo
                if(.not.jf(43)) then
                    write(7,302) pname(7)
                    write(7,402) (oar(10,i),i=1,numi)
                    endif                
            endif 

        if(jf(2)) write(7,3292) topt,tsopt
        if(jf(3)) write(7,329) iopt

        if(ivar.eq.1) then
                if(oar(3,1).lt.1.) oar(4,1)=0.
                yp2=0
                if(oar(3,1).gt.0.0) yp2=oar(3,1)/1.e6
                write(7,213) oar(1,1)/1.E6,yp2,oar(5,1)/1.E6
                write(7,214) oar(2,1),oar(4,1),oar(6,1)
        else
                write(7,307)
        endif

        write(7,211) oar(23,1),oar(25,1),oar(27,1)

        write(7,223) oar(33,1),rzopt
        write(7,2231) oar(39,1),igopt
        write(7,2238) oar(41,1),fdopt
        write(7,2237) oar(46,1),f8opt

        if(htec_max.gt.50.0) write(7,3914) htec_max

3991    format(///'yyyy/mmdd(or -ddd)/hh.h):',I4,'/',I4,'/',F4.1,
     &      A2,2X,A4,' Lat/Long/Alt=',F5.1,'/',F6.1,'/',F6.1/)
3914    format(/'TEC [1.E16 m-2] is obtained by numerical integ',
     &     'ration in 1km steps'/'  from 50 to ',f6.1,' km. ',
     &     't is the percentage of TEC above the F peak.') 
3916    format(/'M3000F2: Propagation factor related to hmF2'/
     &     'B0: bottomside thickness parameter.') 
301     format(A4,' maps are used for the F2 peak density (NmF2)')
302     format(A9,' provided by user:')
402     format(7(1PE10.3))
303     format(A10,'model is used for F2 peak height (hmF2)')
304     format(A25,'=',F5.1,') provided by user')
307     format(1x/'Solar and magnetic parameter for the 1st profile',
     &          ' point:')
309     format(A8,' option is used for the bottomside thickness ',
     &          'parameter B0')
3081    format(A10,' option is used for D-region')
329     format(A11,' option is used for ion composition')
3291    format('foF2 STORM model is turned ',A3)
3299    format('foE auroral storm model is turned ',A3)
3292    format(A8,A5,' solar dependence is used for the electron',
     &          ' temperature')
3293    format(A8,' option is used for the D-region Ne')
3295    format(A16,' option is used for the F1 occurrence probability')
211     format('Solar Zenith Angle/degree',28X,F6.1/
     &          'Dip (Magnetic Inclination)/degree',21X,
     &          F6.2/'Modip (Modified Dip)/degree',27X,F6.2)
223     format('Solar Sunspot Number (12-months running mean) Rz12',
     &          4X,F5.1,A11)
2231    format('Ionospheric-Effective Solar Index IG12',16X,
     &          F5.1,A11)
2238    format('Solar radio flux F10.7 (daily)',24X,
     &          F5.1,A11)
2237    format('Solar radio flux F10.7 (81-day average)',15X,
     &          F5.1,A11)
213     format(/'Peak Densities/cm-3: NmF2=',F9.1,'   NmF1=',F9.1,
     &          '   NmE=',F9.1)
214     format('Peak Heights/km:     hmF2=',F9.2,'   hmF1=',F9.2,
     &          '   hmE=',F9.2/)
3314    format(A7,' is used for topside Ne profile')
c
c table head .......................................................
c

        agnr=7          !output unit number
        xtex=imz(ivar)
        if(jmag.gt.0.and.(ivar.eq.2.or.ivar.eq.3)) xtex='GEOM'
        if(iut.gt.0.and.ivar.eq.8) xtex='U.T.'

        IF(PIKTAB.EQ.4) WRITE(7,8199) 
        IF(PIKTAB.EQ.3) WRITE(7,8191) ITEXT(IVAR),
     &    (pna(pad1(j)),j=1,6),xtex,(uni(pad1(j)),j=1,6)
        IF(PIKTAB.EQ.2) WRITE(7,8194) ITEXT(IVAR),xtex
        IF(PIKTAB.EQ.1) WRITE(7,8192) ITEXT(IVAR),xtex

        IF(PIKTAB.EQ.0) THEN
        	if(jf(22)) then
        		WRITE(7,8193) ITEXT(IVAR),xtex
        	else
        		WRITE(7,9193) ITEXT(IVAR),xtex        	
        	endif
			ENDIF
8191  FORMAT(/'-'/2X,A5,6A10/3X,A4,6A10)
8192  FORMAT(/'-'/2X,A5,6X,'PEAK ALTITUDES IN KM',8X,'PEAK DEN',
     &  'SITIES IN cm-3  TEC top/%'/3X,A4,'    hmF2  hmF1   hmE   ',
     &  'hmD      NmF2   NmF1    NmE    NmD  1E16m-2')
8194  FORMAT(/'-'/2X,A5,3X,'M3000   B0',3X,'B1   E-VALLEY',7X,
     &  'PLASMA ','FREQUENCIES / MHz'/3X,A4,11X,'km       W/km ',
     &  ' Depth',5X,'foF2   foF1   foE   foD')
8193  FORMAT(/'-'/1X,A5,' ELECTRON DENSITY   TEMPERATURES ',
     &  5X,'ION PERCENTAGES[%]*10',4X,'1E16m-2'/2X,A4,' Ne/cm-3 Ne/',
     &  'NmF2 Tn/K  Ti/K  Te/K  O+  N+  H+ He+ O2+ NO+ Clust TEC t/%')
9193  FORMAT(/'-'/1X,A5,' ELECTRON DENSITY   TEMPERATURES ',
     &  4X,'ION DENSITIES[cm-3]/100',4x,'1E16m-2'/2X,A4,' Ne/cm-3 Ne/',
     &  'NmF2 Tn/K  Ti/K  Te/K  O+  N+  H+ He+ O2+ NO+ Clust TEC t/%')
8199  FORMAT(/'-'/1X,'h',8X,' D-REGION ELECTRON DENSITY IN CM-3'/
     &  1X,'km',18X,'DRS-95: Stratos Warming/Winter Anomaly'/5X,
     &  'IRI-07',4x,'FIRI  SW/WA=0/0  0.5/0   1/0    0/0.5    0/1')

c
c output: D-region PIKTAB=4
c
c D-REGION ELECTRON DENSITY IN CM-3: 
c    IRI-07    FIRI  Danilov:SW/WA=0/0  0.5/0   1/0    0/0.5    0/1 
c                    DRS-95: Stratos Warming/Winter Anomaly
c

		if(piktab.eq.4) then
            do 2591 lix=1,77 
            	jdprof(lix)=-1
            	dichte=outf(14,lix)
2591            if(dichte.gt.0.) jdprof(lix)=int(dichte/1.e6+0.5)
			do 2592 lix=1,11
				ihtemp=55+lix*5
            	WRITE(7,3810) ihtemp,jdprof(lix),jdprof(lix+11),
     &  			jdprof(lix+22),jdprof(lix+33),jdprof(lix+44),
     &  			jdprof(lix+55),jdprof(lix+66)
2592		    continue			
3810    FORMAT(I3,7I8)
			goto 2357
		 	endif
		
        xcor=vbeg

        do 1234 li=1,numstp

c
c output: peak densities and altitudes PIKTAB=1
c

      IF(PIKTAB.eq.1) THEN
        if(oar(3,li).lt.1.) oar(4,li)=0.
        iyp1=int(oar(1,li)/1.e6+.5)
        iyp2=0
        if(oar(3,li).gt.0.0) iyp2=int(oar(3,li)/1.e6+.5)
        iyp3=int(oar(5,li)/1.e6+.5)
        iyp4=int(oar(7,li)/1.e6+.5)
            tec=oar(37,li)
        if(tec.gt.0.0) then
            tec=tec/1.e16
            itopp=int(oar(38,li)+.5)
        else
            tec=-1.0
            itopp=-1
        endif
        WRITE(7,3910) XCOR,oar(2,li),oar(4,li),oar(6,li),oar(8,li),
     &    iyp1,iyp2,iyp3,iyp4,tec,itopp
3910    FORMAT(F7.1,2X,4F6.1,1X,I9,3I7,1X,F6.2,I4)
        GOTO 1234
      ENDIF
c
c output: plasma frequencies and profile parameters  PIKTAB=2
c

      IF(PIKTAB.eq.2) THEN
        if(oar(3,li).lt.1.) oar(4,li)=0.
        fyp1=SQRT(oar(1,li)/1.24E10)
        fyp2=0
        if(oar(3,li).gt.0.0) fyp2=SQRT(oar(3,li)/1.24E10)
        fyp3=SQRT(oar(5,li)/1.24E10)
        fyp4=SQRT(oar(7,li)/1.24E10)
            tec=oar(37,li)
        if(tec.gt.0.0) then
            tec=tec/1.e16
            itopp=int(oar(38,li)+.5)
        else
            tec=-1.0
            itopp=-1
        endif
        wvalley=oar(12,li)-oar(6,li)
        dvalley=0.0
        if(oar(5,li).gt.0.0) dvalley=oar(11,li)/oar(5,li)
        WRITE(7,3950) XCOR,oar(36,li),oar(10,li),oar(35,li),
     &    wvalley,dvalley,fyp1,fyp2,fyp3,fyp4
3950    FORMAT(F7.1,2X,F6.4,F6.1,F4.1,F6.1,F8.4,1X,4F7.3)
        GOTO 1234
      ENDIF
c
c output: 6 parameters of your choice    PIKTAB=3
c

      IF(PIKTAB.eq.3) THEN
c        if(pad1.eq.45.and.oar(pad1,li).le.0.0) oar(pad1,li)=-1.
c        if(pad2.eq.45.and.oar(pad2,li).le.0.0) oar(pad2,li)=-1.
c        if(pad3.eq.45.and.oar(pad3,li).le.0.0) oar(pad3,li)=-1.
        WRITE(7,3919) XCOR,oar(pad1(1),li),oar(pad1(2),li),
     &        oar(pad1(3),li),oar(pad1(4),li),oar(pad1(5),li),
     &        oar(pad1(6),li)
3919    FORMAT(F7.1,6(1X,1PE9.2))
        GOTO 1234
      ENDIF
c
c output: special for test purposes    PIKTAB=5
c

      IF(PIKTAB.eq.5) THEN
c ----------- B0, B1 ----------------
c        WRITE(8,4919) XCOR,oar(10,li),oar(35,li)
c4919    FORMAT(F7.1,2X,F6.2,2X,F5.3)
c ----------- Te ----------------
        WRITE(8,4919) XCOR,outf(2,li),outf(3,li),outf(4,li)
4919    FORMAT(F7.1,2X,F6.1,2X,F6.1,2X,F6.1)
c ----------- Ne, TEC ----------------
c        WRITE(8,4919) XCOR,outf(1,li),oar(37,li)
c4919    FORMAT(F7.1,2X,E12.5,2X,E12.5)
c ----------- Ni ----------------
c        type*,XCOR,outf(1,li),outf(5,li),outf(6,li),
c     &   outf(7,li),outf(8,li),outf(9,li),outf(10,li),outf(11,li)
c        WRITE(8,4919) XCOR,outf(1,li),outf(5,li),outf(6,li),
c     &   outf(7,li),outf(8,li),outf(9,li),outf(10,li),outf(11,li)
c4919    FORMAT(F7.1,2X,E12.5,7F10.4)
c ----------- hmF2 ----------------
c        type*,XCOR,oar(26,li),oar(2,li),oar(36,li),
c     &   sqrt(oar(1,li)/oar(5,li)),oar(33,li),oar(39,li)
c        WRITE(8,4919) XCOR,outf(1,li),outf(5,li),outf(6,li),
c     &   outf(7,li),outf(8,li),outf(9,li),outf(10,li),outf(11,li)
c4919    FORMAT(F7.1,2X,E12.5,7F10.4)
c ----------- auroral boundary ----------------
c        print *,XCOR,oar(55,li),oar(56,li),oar(54,li),
c     &     oar(57,li),oar(54,li)-oar(57,li),oar(58,li)
cc        WRITE(8,4919) XCOR,outf(1,li),outf(5,li),outf(6,li),
cc     &   outf(7,li),outf(8,li),outf(9,li),outf(10,li),outf(11,li)
cc4919    FORMAT(F7.1,2X,E12.5,7F10.4)

        GOTO 1234
      ENDIF
c
c output: standard
c
        if(ivar.eq.1) then
                oar(1,li)=oar(1,1)
                oar(37,li)=oar(37,1)
                oar(38,li)=oar(38,1)
                endif
        jne=int(outf(1,li)/1.e6+.5)
        xner=outf(1,li)/oar(1,li)
        jtn=int(outf(2,li)+.5)
        jti=int(outf(3,li)+.5)
        jte=int(outf(4,li)+.5)
        scid=1.0E-8
        if(jf(22)) scid=10.
        jio=INT(OUTF(5,li)*scid+.5)
        jih=INT(OUTF(6,li)*scid+.5)
        jihe=INT(OUTF(7,li)*scid+.5)
        jio2=INT(OUTF(8,li)*scid+.5)
        jino=INT(OUTF(9,li)*scid+.5)
        jicl=INT(OUTF(10,li)*scid+.5)
        jin=INT(OUTF(11,li)*scid+.5)
c        print *,'O+,N+,O2+,NO+=',outf(5,li),
c     &    outf(11,li),outf(9,li),outf(8,li)
        if(outf(1,li).lt.0) jne=-1
        if(outf(1,li).lt.0) xner=-1.
        if(outf(2,li).lt.0) jtn=-1
        if(outf(3,li).lt.0) jti=-1
        if(outf(4,li).lt.0) jte=-1
        if(outf(5,li).lt.0) jio=-1
        if(outf(6,li).lt.0) jih=-1
        if(outf(7,li).lt.0) jihe=-1
        if(outf(8,li).lt.0) jio2=-1
        if(outf(9,li).lt.0) jino=-1
        if(outf(10,li).lt.0) jicl=-1
        if(outf(11,li).lt.0) jin=-1
            tec=oar(37,li)
        if(tec.gt.0.0) then
            tec=tec/1.e16
            itopp=int(oar(38,li)+.5)
        else
            tec=-1.0
            itopp=-1
        endif
c        print *, XCOR,jne,xner,jtn,jti,jte,jio,jin,
c     &        jih,jihe,jino,jio2,jicl,tec,itopp
        WRITE(7,7117) XCOR,jne,xner,jtn,jti,jte,jio,jin,
     &        jih,jihe,jio2,jino,jicl,tec,itopp
7117    FORMAT(F6.1,I8,1x,F6.3,3I6,7I4,f6.1,i4)

1234    xcor=xcor+vstp

2357    print *,'Enter 0 to exit or 1 to generate another profile?' 
        read(5,*) icontinue
        if (icontinue.gt.0) goto 1
		    
            stop
            end
