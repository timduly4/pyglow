c irifun.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
C Functions and subroutines for the International Reference Ionosphere 
C (IRI) model. These functions and subroutines are called by the main
C IRI subroutine IRI_SUB in IRISUB.FOR.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Required i/o units:  
c  KONSOL= 6 Program messages (used when jf(12)=.true. -> konsol)
c  KONSOL=11 Program messages (used when jf(12)=.false. -> MESSAGES.TXT)
c
c     COMMON/iounit/konsol,mess is used to pass the value of KONSOL from 
c     IRISUB to IRIFUN and IGRF. If mess=false then messages are turned off.
c     
c  UNIT=12 TCON: Solar/ionospheric indices IG12, R12 (IG_RZ.DAT) 
c  UNIT=13 APF,APFMSIS,APF_ONLY: Magnetic indices and F10.7 (APF107.DAT) 
c
c I/o Units used in other programs:
c  IUCCIR=10 in IRISUB for CCIR and URSI coefficients (CCIR%%.ASC, %%=month+10)
c  UNIT=14 in IGRF/GETSHC for IGRF coeff. (DGRF%%%%.DAT, %%%%=year)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c changes from IRIFU9 to IRIF10:
c       SOCO for solar zenith angle 
c       ACOS and ASIN argument forced to be within -1 / +1
c       EPSTEIN functions corrected for large arguments
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c changes from IRIF10 to IRIF11: 
c       LAY subroutines introduced
c       TEBA corrected for 1400 km
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c changes from IRIF11 to IRIF12:
C       Neutral temperature subroutines now in CIRA86.FOR 
C       TEDER changed
C       All names with 6 or more characters replaced 
C       10/29/91 XEN: 10^ in loop, instead of at the end
C       1/21/93 B0_TAB instead of B0POL
C       9/22/94 Alleviate underflow condition in IONCOM exp()
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c changes from IRIF12 to IRIF13:
C        9/18/95 MODA: add leap year and number of days in month
C        9/29/95 replace F2out with FOUT and XMOUT.
C       10/ 5/95 add TN and DTNDH; earlier in CIRA86.FOR
C       10/ 6/95 add TCON for reading indices
C       10/20/95 MODA: IN=1 MONTH=IMO
C       10/20/95 TCON: now includes RZ interpolation
C       11/05/95 IONCOM->IONCO1, added IONCOM_new, IONCO2
C       11/05/95 LSTID added for strom-time updating
C       11/06/95 ROGUL: transition 20. instead of 15.
C       12/01/95 add UT_LT for (date-)correct UT<->LT conversion
C       01/16/96 TCON: add IMST to SAVE statement
C       02/02/96 ROGUL: 15. reinstated
C       02/07/96 UT_LT: ddd, dddend integer, no leap year 2000
C       03/15/96 ZERO: finding delta for topside
C       03/18/96 UT_LT: mode=1, change of year
C       12/09/96 since 2000 is leap, delete y/100*100 condition
C       04/25/97 XMDED: minimal value also daytime
C       05/18/98 TCON: changes to IG_RZ (update date); -R = Cov
C       05/19/98 Replaced IONCO2&APROK; HEI,XHI in IONCOM_NEW
C       10/01/98 added INITIALIZE
C       04/30/99 MODA: reset bb(2)=28
C       11/08/99 avoid negative x value in function XE2. Set x=0.
C       11/09/99 added COMMON/const1/humr,dumr also for CIRA86
C       11/30/99 EXIT in APROK replaced with GOTO (N. Smirnova)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c changes from IRIF13 to IRIFUN:
C-Version-mm/dd/yy-description [person reporting correction]
C 2000.01 05/09/00 Include B0_98 subroutine to replace B0_NEW and B0POL
C 2000.02 05/18/00 Include Elteik and spharm_ik for Te
C 2000.03 06/09/00 Include xe3_1, xe4_1, xe_1
C 2000.04 06/11/00 Include f1_c1, f1_prob, modified fof1ed
C 2000.05 10/30/00 Include vdrift
C 2000.06 04/15/01 Include IGRF_SUB subroutine for IK Te model
C 2000.07 05/07/01 Include storm subroutine STORM and Ap access s/w
C 2000.08 09/07/01 APF: if(j1.eq.j2) -> if(IY.eq.j2) [P. Wilkinson]
C 2000.09 09/07/01 CONVER: LO2 = MOD(LO1,20)+1 [P. Webb,D. Pesnell]
C 2000.10 02/20/02 CONVER/DATA: 105.78 -> 015.78 [A. Shovkoplyas] 
C 2000.11 10/28/02 replace TAB/6 blanks, enforce 72/line [D. Simpson]
C 2000.12 11/08/02 removing unused variables (corr); apf0 removed
C 2000.13 11/26/02 apf() using keyed access to ap.dat file; apf->apf1
C 2000.14 11/27/02 changed F1_PROB; always 6 preceeding spaces 
C 2005.01 03/09/05 CALION,INVDPC,CALNE for new Ne, Ni models
C 2005.01 11/14/05 APF_ONLY for F107D;  
C 2005.01 11/14/05 spreadf_brazil; added constraint 0<=P<=1 
C 2005.02 05/11/06 NeQuick: XE1,TOPQ, M3000HM; stormvd,
C 2005.02 03/27/07 STORM: hourly interpolation of Ap  [A. Oinats]
C 2007.00 05/18/07 Release of IRI-2007
C 2007.01 09/19/07 vdrift et al.: without *8 (no change in results)
C 2007.04 02/07/09 IONLOW: N+ correction [V. Truhlik]
C 2007.05 03/30/09 NMDED: avoid exp underflow [K. Choi] 
C 2007.05 03/30/09 spreadf_brazil: bspl2f et al b(20->30) [Tab Ji]
C 2007.05 03/30/09 APF_ONLY: Compute monthly F10.7
C 2007.06 05/26/09 APF_ONLY: replace i with 1 and IMN with ID [R.Conde]
C 2007.07 07/10/09 CONVER/DATA: 015.78 -> 005.78 [E. Araujo] 
C 2007.08 07/23/09 STORM/CONVER: long. discont. [R. Conde, E. Araujo] 
C 2007.08 07/23/09 APF,APF_ONLY: use YearBegin from ap.dat [R. Conde] 
C 2007.10 02/03/10 APF: eof error message; clean-up APF and APF_only
C 2007.11 04/19/10 ELTEIK:     IF (ALT .GE. 900) THEN      [A. Senior]
C 2007.11 04/19/10 INILAY: HFFF,XFFF when NIGHT=F1REG=f    [A. Senior]
C
C 2012.00 10/05/11 IRI-2012: bottomside B0 B1 model (SHAMDB0D, SHAB1D),
C 2012.00 10/05/11   bottomside Ni model (iriflip.for), auroral foE
C 2012.00 10/05/11   storm model (storme_ap), Te with PF10.7 (elteik),
C 2012.00 10/05/11   oval kp model (auroral_boundary),IGRF-11(igrf.for), 
C 2012.00 10/05/11   NRLMSIS00 (cira.for), CGM coordinates, F10.7 daily
C 2012.00 10/05/11   81-day 365-day indices (apf107.dat), ap->kp (ckp),
C 2012.00 10/05/11   array size change jf(50) outf(20,1000), oarr(100).
C 2012.00 10/05/11   No longer needed: TEDER,DTNDH
C 2012.01 11/08/11 SHAMDB0D+: Initialize COMMON variables outside DATA
C 2012.01 11/08/11 SCHNEVPD,LEGFUN: COSD(x) -> COS(x*UMR)
C 2012.01 11/08/11 auroral_boundary: i1.gt/ge.48
C 2012.01 11/08/11 CALION: F107D upper/lower boundary 220/65
C 2012.01 01/06/12 delete PAUSE statement in SPLINT
C 2012.01 01/24/12 STORME_AP: change 365 to 366 leap year    [F. Simoes]
C 2012.01 06/30/12 HMF2ED: hmF2 too low foF2/foE ge 1.7 [I.Zakharenkova]
C 2012.01 09/14/12 SHAMDB0D,SHAB1D: initializing CONS2      [P. Coisson]
C 2012.02 12/12/12 STORME_AP: add KONSOL and ERROR ouput STORME_AP=-5.
C 2014.00 01/22/14 TPCORR: INVDIP --> INVDP                [J.K. Knight]
C 2014.02 07/24/14 COMMON/iounit/ added 'mess'
C 2014.03 08/25/14 ELTEIK,INVDPC,SOCO: ACOS: if(abs(x).gt.1) x=sign(1,x)
C 2014.03 10/02/14 IONLOW,IONHIGH: C1(82) added for SPHARM_IK [J.C. Xue]
C 2014.05 12/22/14 APFMSIS: changed text if out of range of APF107.DAT
C 2015.01 04/27/15 TCON: ionoindx(806),indrz(806)
C 2015.02 07/12/15 add read_ig_rz, readapf107 change TCON,APF* 
C 2015.02 07/20/15 APFMSIS: remove duplicate iiap(8)         [E. Blanch]
C 2015.03 08/13/15 DATA in subr: DATA values only used first CALL 
C 2015.03 08/13/15 COMMON/CONST/UMR,PI 
C 2015.03 09/30/15 hmF2 new: SHAMDHMF2, SDMF2 and associated subroutines 
C 2015.03 09/30/15 ELTEIK,CALION,IONLOW,IONHIGH w/o invdip calc; INVDPC
C 2015.03 10/12/15 READAPF107: F365 -> F107_365             [M. Hausman]
C 2015.03 10/14/15 LEGFUN: replace print * with write(konsol,..)
C 2015.03 10/14/15 SHAMDB0D,SHAB1D,SCHNEVPD: COMMON/ATB/
C 2015.03 10/14/15 CLCMLT,DPMTRX --> IGRF.FOR
C 2015.04 02/01/16 TAL: if(SHBR.le.0.0) -> RETURN
C 2015.04 02/01/16 IONLOW,IONHIGH: ALT -> ALTI extrapolation  [S.R. Zhang]
C 2016.01 02/24/16 SCHNEVPH: COSD(X) -> COS(X*UMR); CONST/UMR   [W. Toler]
C 2016.01 02/24/16 SCHNEVP,LEGFUN: COMMON/COST/dfarg,PI         [W. Toler]
C 2016.02 03/23/16 CALION, IONLOW, IONHIGH, INVDPC  revised   [V. Truhlik]
C 2016.03 07/19/16 XE3_1: change D1F1 to C1                    [I. Galkin]
C 2016.04 09/08/16 CALION: Version 2.5 C/NOFS correction      [V. Truhlik]
C 2016.04 09/08/16 NEW: model_hmF2                             [V. Shubin]
C 2016.05 10/19/16 read_ig_rz: *0.7 for r12_new starting 01/2014
C 2017.01 02/23/17 SHAB1D: new SCHNEVPDB1 and COMMON/ATB1/ 
C 2018.01 03/22/18 INVDPC= ... ALFA*SIGN(1.0,DIPL)*INVL       [V. Truhlik]
C 2018.01 03/22/18 INVDPC_OLD for ELTEIK                      [V. Truhlik]
C                  
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c IRI functions and subroutines:
C Ne:       XE1,TOPQ,ZERO,DXE1N,XE2,XE3_1,XE4_1,XE5,XE6,XE_1
C Te, Ti:   ELTEIK{INTERP,KODERR,KOEFD,KOF107,LOCATE,SPHARM_IK, 
C		    SPLINE,SPLINT,SWAPEL,TEDIFI,TPCAS,TPCORR},TEBA,SPHARM,
C           ELTE,TEDE,TI,TN
C Ni:       RPID,RDHHE,RDNO,KOEFP1,KOEFP2,KOEFP3,SUFE,IONDANI,IONCO1, 
C           IONCO2,APROK,CALION,IONLOW,IONHIGH,INVDPC
C PEAKS:    FOUT,XMOUT,HMF2ED,XM3000HM,SHAMDHMF2,SCHNEVPDH,model_hmF2,
C     		SDMF2,hmF2_med_SD,read_data_SD,fun_hmF2_SD,fun_Gk,Legendre, 
C     		fun_hmF2UT,Koeff_UT,fun_Akp_UT,fun_Fk_UT,fun_Gk_UT
C     		FOF1ED,f1_c1,f1_prob,FOEEDI,XMDED,GAMMA1
C PROFILE:  TOPH05,CHEBISH,SHAMDB0D,SHAB1D,SCHNEVPD,TBFIT,LEGFUN,  
C           B0_98,TAL,VALGUL,DREGION
C MAG. FIELD: FIELDG, CONVER(Geom. Corrected Latitude)
C TIME:     SOCO,HPOL,MODA,UT_LT,SUN
C EPSTEIN:  RLAY,D1LAY,D2LAY,EPTR,EPST,EPSTEP,EPLA
C LAY:      XE2TO5,XEN,ROGUL,LNGLSN,LSKNM,INILAY
C INDICES:  TCON,APF,APFMSIS,APF_ONLY
C STORM:   	CONVER, STORM, STORME_AP
C Vi:       vdrift,bspl4_time,bspl4_long,g,stormvd,bspl4_ptime
C Spread-F:	spreadf_brazil,bspl2f,bspl2l,bspl2s,bspl4t
C Auroral: 	auroral_boundary, ckp
C Misc:     REGFA1
c-----------------------------------------------------------------------
C
C        
C*************************************************************   
C*************** ELECTRON DENSITY ****************************   
C*************************************************************   
C

        FUNCTION XE1(H)    
c----------------------------------------------------------------
C DETERMINING ELECTRON DENSITY(M-3) IN THE TOPSIDE IONOSPHERE   
C (H=HMF2....2000 KM) BY HARMONIZED BENT-MODEL ADMITTING 
C VARIABILITY OF THE GLOBAL PARAMETERS BETA,ETA,DELTA,ZETA WITH        
C GEOM. LATITUDE, SMOOTHED SOLAR FLUX AND CRITICAL FREQUENCY.     
C BETA,ETA,DELTA,ZETA are computed in IRISUB program and 
C communicated via COMMON /BLO10. This is the IRI-2001 approach
C [REF.:K.RAWER,S.RAMAKRISHNAN,1978] 
C New options include:
C (1) IRI-corrected: TC3,alg10,hcor1 in COMMON /BLO11. 
C   TC3     correction term divided by (1500-(hcor1-hmF2))
C   alg10   = alog(10.)
C	hcor1	lower height boundary for correction
C (2) NeQuick:  B2TOP  in COMMON /BLO11.
C	B2TOP   is the topside scale height that depends on foF2 and 
C           hmF2. 
C Switch for choosing the desired option is itopn in COMMON /BLO11
C   itopn   =0 IRI-2001, =1 IRI-2001-corrected, =2 NeQuick
C           =3 Gulyaeva-0.5 is not yet implemented. 
c----------------------------------------------------------------
        COMMON  /BLOCK1/HMF2,XNMF2,HMF1,F1REG
     &          /BLO10/BETA,ETA,DELTA,ZETA
     &          /BLO11/B2TOP,TC3,itopn,alg10,hcor1
     &          /QTOP/Y05,H05TOP,QF,XNETOP,xm3000,hhalf,tau
     &          /ARGEXP/ARGMAX

        logical 	f1reg              

        IF(itopn.eq.2) THEN
          XE1=TOPQ(H,XNMF2,HMF2,B2TOP)
          RETURN
        ENDIF
      
        DXDH = (1000.-HMF2)/700.
        x0 = 300. - delta
        xmx0 = (H-HMF2)/DXDH
        x = xmx0 + x0
        eptr1 = eptr(x,beta,394.5) - eptr(x0,beta,394.5)
        eptr2 = eptr(x,100.,300.0) - eptr(x0,100.,300.0) 
        y = BETA * ETA * eptr1 + ZETA * (100. * eptr2 - xmx0)
        Y = y * dxdh
        if(abs(Y).gt.argmax) Y = sign(argmax,Y)

        IF(itopn.eq.3) then
	         IF((QF.EQ.1.).AND.(ABS(H-H05TOP).LT.1.)) QF=Y05/Y
             XE1 = XNMF2 * EXP(-Y*QF)                             
             RETURN          
             endif
        TCOR = 0.
        IF(itopn.eq.1.and.h.gt.hcor1) then
             xred = h - hcor1
             rco = tc3 * xred
             TCOR = rco * alg10
             endif
        XE1 = XNMF2 * EXP(-Y+TCOR)                             
        RETURN          
        END             
C
C

        REAL FUNCTION TOPQ(h,No,hmax,Ho)
c----------------------------------------------------------------
c  NeQuick formula
c----------------------------------------------------------------
        REAL No
        PARAMETER (g=0.125,rfac=100.0)
          dh=h-hmax
          g1=g*dh
          z=dh/(Ho*(1.0+rfac*g1/(rfac*Ho+g1)))
          if(z.gt.40) then
            topq=0.0
            return
          endif
          ee=exp(z)
          if (ee.gt.1.0e7) then
            ep=4.0/ee
          else
            ep=4.0*ee/(1.0+ee)**2
          endif
          TOPQ=No*ep
        RETURN
        END

C 
C
        REAL FUNCTION ZERO(DELTA)
C FOR A PEAK AT X0 THE FUNCTION ZERO HAS TO BE EQUAL TO 0.
        COMMON  /BLO10/         BETA,ETA,DEL,ZETA
     &          /ARGEXP/        ARGMAX

        arg1=delta/100.
        if (abs(arg1).lt.argmax) then
                z1=1./(1.+exp(arg1))
        else if (arg1.lt.0) then
                z1=1.
        else
                z1=0.
        endif
        arg1=(delta+94.5)/beta
        if (abs(arg1).lt.argmax) then
                z2=1./(1.+exp(arg1))
        else if (arg1.lt.0) then
                z2=1.
        else
                z2=0.
        endif
        zero=zeta*(1.-z1) - eta*z2
        return
        end
C
C
        FUNCTION DXE1N(H)                            
C LOGARITHMIC DERIVATIVE OF FUNCTION XE1 (KM-1).   
        COMMON    /BLOCK1/HMF2,XNMF2,HMF1,F1REG
     &            /BLO10/BETA,ETA,DELTA,ZETA                    
	    logical f1reg

        x0 = 300. - delta
        X=(H-HMF2)/(1000.0-HMF2)*700.0 + x0
        epst2 = epst(x,100.0,300.0)
        epst1 = epst(x,beta ,394.5)
        DXE1N = - ETA * epst1 + ZETA * (1. - epst2)             
        RETURN          
        END             
C
C
        REAL FUNCTION XE2(H)                         
C ELECTRON DENSITY FOR THE BOTTOMSIDE F-REGION (HMF1...HMF2).                   
        COMMON    /BLOCK1/HMF2,XNMF2,HMF1,F1REG
     &          /BLOCK2/B0,B1,C1        /ARGEXP/ARGMAX
	    logical	f1reg

        X=(HMF2-H)/B0
        if(x.le.0.0) x=0.0
        z=x**b1
        if(z.gt.argmax) z=argmax
        XE2=XNMF2*EXP(-z)/COSH(X)                 
        RETURN          
        END             
C
C
        REAL FUNCTION XE3_1(H)
C ELECTRON DENSITY FOR THE F1-LAYER (HZ.....HMF1)
C USING THE NEW DEFINED F1-LAYER FUNCTION (Reinisch and Huang, Advances 
C in Space Research, Volume 25, Number 1, 81-88, 2000)
        COMMON	/BLOCK1/	HMF2,XNMF2,HMF1,F1REG
     &		/BLOCK2/	B0,B1,C1
	    logical	f1reg
C
	    h1bar=h
        if (f1reg) H1BAR=HMF1*(1.0-((HMF1-H)/HMF1)**(1.0+C1))
        XE3_1=XE2(H1BAR)
        RETURN
        END
C
C
        REAL FUNCTION XE4_1(H)
C ELECTRON DENSITY FOR THE INTERMEDIATE REGION (HEF...HZ)
C USING THE NEW DEFINED FUNCTION
        COMMON	/BLOCK3/	HZ,T,HST
     &		/BLOCK4/	HME,XNME,HEF
C
	    if(hst.lt.0.0) then
		xe4_1=xnme+t*(h-hef)
		return
		endif
        IF(HST.EQ.HEF) THEN
           H1BAR=H
        ELSE
           H1BAR=HZ+0.5*T-SIGN(1.0,T)*SQRT(T*(0.25*T+HZ-H))
        ENDIF
        XE4_1=XE3_1(H1BAR)
        RETURN
        END
C
C
        REAL FUNCTION XE5(H)                         
C ELECTRON DENSITY FOR THE E AND VALLEY REGION (HME..HEF).   
        LOGICAL NIGHT   
        COMMON    /BLOCK4/        HME,XNME,HEF
     &          /BLOCK5/        NIGHT,E(4)                    
        T3=H-HME        
        T1=T3*T3*(E(1)+T3*(E(2)+T3*(E(3)+T3*E(4))))  
        IF(NIGHT) GOTO 100                           
        XE5=XNME*(1+T1)  
        RETURN          
100     XE5=XNME*EXP(T1)                              
        RETURN          
        END             
C
C
        REAL FUNCTION XE6(H)                         
C ELECTRON DENSITY FOR THE D REGION (HA...HME).    
        COMMON    /BLOCK4/        HME,XNME,HEF
     &          /BLOCK6/        HMD,XNMD,HDX
     &        /BLOCK7/        D1,XKK,FP30,FP3U,FP1,FP2    
        IF(H.GT.HDX) GOTO 100                        
        Z=H-HMD         
        FP3=FP3U        
        IF(Z.GT.0.0) FP3=FP30                        
        XE6=XNMD*EXP(Z*(FP1+Z*(FP2+Z*FP3)))           
        RETURN          
100     Z=HME-H         
        XE6=XNME*EXP(-D1*Z**XKK)
        RETURN          
        END             
C
C
        REAL FUNCTION XE_1(H)                          
C ELECTRON DENSITY BEETWEEN HA(KM) AND 1000 KM     
C SUMMARIZING PROCEDURES  NE1....6;                
        COMMON    /BLOCK1/HMF2,XNMF2,XHMF1,F1REG         
     &          /BLOCK3/HZ,T,HST
     &        /BLOCK4/HME,XNME,HEF
	    logical 	f1reg
	    if(f1reg) then
		   hmf1=xhmf1
	    else
		   hmf1=hmf2
	    endif
        IF(H.LT.HMF2) GOTO 100                       
        XE_1=XE1(H)     
        RETURN          

100     IF(H.LT.HMF1) GOTO 300                       
        XE_1=XE2(H)       
        RETURN          

300     IF(H.LT.HZ) GOTO 400                         
        XE_1=XE3_1(H)       
        RETURN          

400     IF(H.LT.HEF) GOTO 500                        
        XE_1=XE4_1(H)       
        RETURN          

500     IF(H.LT.HME) GOTO 600                        
        XE_1=XE5(H)       
        RETURN          

600     XE_1=XE6(H)       
        RETURN          
        END             
C
C                     
C**********************************************************                     
C***************** ELECTRON TEMPERATURE ********************                    
C**********************************************************                     
C
      SUBROUTINE ELTEIK(PF107Y,INVDIP,MLT,ALT,DDD,PF107,TE,SIGTE)
c      SUBROUTINE ELTEIK(CRD,PF107Y,INVDIP,FL,DIMO,B0,
c     &                   DIPL,MLT,ALT,DDD,PF107,TE,SIGTE)
C----------------------------------------------------------------------
C Empirical model of electron temperature (Te) in the outer ionosphere
C with inclusion of solar activity.
C Based on spherical harmonics approximation of measured
C Te (all available satellites) at altitudes centred on 350km, 550km,
C 850km, 1400km, and 2000km. For intermediate altitudes a linear
C interpolation is used. Recommended altitude range: 300-2500 km!!!
C Linear extrapolation is used for altitude ranges <300;350)km
C and (2000;2500> km. For days between seasons centred at
C (21.3. = 79; 21.6. = 171; 23.9. 265; 21.12. = 354) Te is
C interpolated by a harmonic function.
Cc Inputs: CRD - 0 .. INVDIP
Cc               1 .. FL, DIMO, B0, DIPL (used for calc. INVDIP inside)
C         PF107Y - 0 .. PF107 correction NOT included
C                  1 .. PF107 correction included
C         INVDIP - "mix" coordinate of the dip latitude and of
C                the invariant latitude (can be computed with INVDPC)
C                positive northward, in deg, range <-90.0;90.0>
Cc         FL, DIMO, BO - McIlwain L parameter, dipole moment in
Cc                        Gauss, magnetic field strength in Gauss -
Cc                        parameters needed for invariant latitude
Cc                        calculation
Cc         DIPL - dip latitude
Cc                positive northward, in deg, range <-90.0;90.0>
C         MLT - magnetic local time (central dipole)
C               in hours, range <0;24)
C         ALT - altitude above the Earth's surface;
C               in km, range <500;3000>
C         DDD - day of year; range <0;365>
C         PF107 - Phil Richard's solar radio flux;
C Output: TE - electron temperature in K
C         SIGTE - standard deviation (or model error) of TE in K
C Versions: 1.00 (IDL) the first version Te=Te(invl,mlt,alt,season)
C           1.50 (IDL) corrected IK19 Te at 900km for possible Ne>2E11 m-3
C           2.00 (IDL) F107 included as a linear perturbation on global Te pattern
C                      Te=Te(invlat,mlt,alt,season,F107)
C           3.00 (IDL) invdipl introduced
C           2000 (IDL,FORTRAN) correction for seasons included
C           2010 (IDL,FORTRAN) completely new version 
C Authors of the model (v 2011)
C                V. Truhlik, D. Bilitza, and L. Triskova
C Author of the code:
C         Vladimir Truhlik
C         Institute of Atm. Phys.
C         Bocni II.
C         141 31 Praha 4, Sporilov
C         Czech Republic
C         e-mail: vtr@ufa.cas.cz
C----------------------------------------------------------------------
c      REAL INVDIP,FL,DIMO,B0,DIPL,MLT,ALT,PF107,TE,SIGTE
c      INTEGER CRD,PF107Y,DDD,SEZDAY,XDAY
      REAL INVDIP,MLT,ALT,PF107,TE,SIGTE
      INTEGER PF107Y,DDD,SEZDAY,XDAY
      INTEGER MIRREQ(81)
      REAL D(5,3,81),DERRTE(5,3,81),DPF107(5,3,81)
      DOUBLE PRECISION B(8),A
      REAL DPI,DTOR,ASA,INVL,RINVL,INVDP,RDIPL,ALFA,BETA
      REAL RMLT,RCOLAT
      REAL C(82)
      INTEGER SEZA,SEZB,DDDA,DDDB,DDDD
      REAL T350,T350A,T350B,T550,T550A,T550B,T850,T850A,T850B,
     &     T1400,T1400A,T1400B,T2000,T2000A,T2000B 
      REAL P350A,P350B,P550A,P550B,P850A,P850B,
     &     P1400A,P1400B,P2000A,P2000B 
      REAL E350,E350A,E350B,E550,E550A,E550B,E850,E850A,E850B,
     &     E1400,E1400A,E1400B,E2000,E2000A,E2000B 
      REAL TP350A,TP350B,TP550A,TP550B,TP850A,TP850B,
     &     TP140A,TP140B,TP200A,TP200B
      INTEGER FUN
      INTEGER I
      DATA B/1.259921D0  ,-0.1984259D0 ,-0.04686632D0,-0.01314096D0,
     &      -0.00308824D0, 0.00082777D0,-0.00105877D0, 0.00183142D0/
C////////////////////////////////coefficients - main model part//////////////////////
      DATA (MIRREQ(J),J=1,81)/
     &  1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,
     &  1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,
     & -1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1, 1,-1,
     &  1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1, 1,-1, 1,-1, 1, 1/
      CALL KOEFD(MIRREQ,D)
      CALL KODERR(MIRREQ,DERRTE)
      CALL KOF107(MIRREQ,DPF107)
C//////////////////////thresholds of solar activity/////////////////////////////////////
      IF (PF107 .GT. 250) PF107=250
      IF (PF107 .LT. 80) PF107=80
C////////////////////////////////////////////////////////////////////////////////////
      DPI=3.1415926535897
      DTOR=DPI/180.0
c      IF (CRD .EQ. 1) THEN
C      calculation of INVDIP from FL, DIMO, BO, and DIPL
C      invariant latitude calculated by highly
C      accurate polynomial expansion
c       A=(DIMO/B0)**(1.0D0/3.0D0)/FL
c       ASA=A*(B(1)+B(2)*A+B(3)*A**2+B(4)*A**3+B(5)*A**4+
c     &        B(6)*A**5+B(7)*A**6+B(8)*A**7)
c       IF (ASA .GT. 1.0) ASA=1.0
c       IF (ASA .LT. 0.0) ASA=0.0
C      invariant latitude (absolute value)
c       RINVL=ACOS(SQRT(ASA))
c       INVL=RINVL/DTOR
c       RDIPL=DIPL*DTOR
c       ALFA=SIN(ABS(RDIPL))**3
c       BETA=COS(RINVL)**3
c       INVDP=(ALFA*SIGN(1.0,DIPL)*INVL+BETA*DIPL)/(ALFA+BETA)
c      ELSE IF	(CRD .EQ. 0) THEN
c       INVDP=INVDIP
c      ELSE
c       RETURN
c      END IF

      INVDP=INVDIP
      RMLT=MLT*DTOR*15.0
      RCOLAT=(90.0-INVDP)*DTOR
      CALL SPHARM_IK(C,8,8,RCOLAT,RMLT)
C     21.3. - 20.6.
      IF ((DDD .GE. 79) .AND. (DDD .LT. 171)) THEN
       SEZA=1
       SEZB=2
       DDDA=79
       DDDB=171
       DDDD=DDD
       FUN=0
      END IF
C     21.6. - 22.9.
      IF ((DDD .GE. 171) .AND. (DDD .LT. 265)) THEN
       SEZA=2
       SEZB=1
       DDDA=171
       DDDB=265
       DDDD=DDD
       FUN=1
      END IF
C     23.9. - 20.12.
      IF ((DDD .GE. 265) .AND. (DDD .LT. 354)) THEN
       SEZA=1
       SEZB=3
       DDDA=265
       DDDB=354
       DDDD=DDD
       FUN=0
      END IF
C     21.12. - 20.3.
      IF ((DDD .GE. 354) .OR. (DDD .LT. 79)) THEN
       SEZA=3
       SEZB=1
       DDDA=354
       DDDB=365+79
       DDDD=DDD
        IF (DDD .GE. 354) THEN
         DDDD=DDD
        ELSE
         DDDD=DDD+365
        END IF
       FUN=1 
      END IF
C     model Te
      T350A=0.0
      T350B=0.0
      T550A=0.0
      T550B=0.0
      T850A=0.0
      T850B=0.0
      T1400A=0.0
      T1400B=0.0
      T2000A=0.0
      T2000B=0.0                
      DO 30 I=1,81
       T350A=T350A+C(I)*D(1,SEZA,I)
       T350B=T350B+C(I)*D(1,SEZB,I)
       T550A=T550A+C(I)*D(2,SEZA,I)
       T550B=T550B+C(I)*D(2,SEZB,I)
       T850A=T850A+C(I)*D(3,SEZA,I)
       T850B=T850B+C(I)*D(3,SEZB,I)
       T1400A=T1400A+C(I)*D(4,SEZA,I)
       T1400B=T1400B+C(I)*D(4,SEZB,I)
       T2000A=T2000A+C(I)*D(5,SEZA,I)
30     T2000B=T2000B+C(I)*D(5,SEZB,I)
      T350A=10**T350A
      T350B=10**T350B
      T550A=10**T550A
      T550B=10**T550B
      T850A=10**T850A
      T850B=10**T850B
      T1400A=10**T1400A
      T1400B=10**T1400B
      T2000A=10**T2000A
      T2000B=10**T2000B
C     model PF107
      P350A=0.0
      P350B=0.0
      P550A=0.0
      P550B=0.0
      P850A=0.0
      P850B=0.0
      P1400A=0.0
      P1400B=0.0
      P2000A=0.0
      P2000B=0.0                
      DO 40 I=1,81
       P350A=P350A+C(I)*DPF107(1,SEZA,I)
       P350B=P350B+C(I)*DPF107(1,SEZB,I)
       P550A=P550A+C(I)*DPF107(2,SEZA,I)
       P550B=P550B+C(I)*DPF107(2,SEZB,I)
       P850A=P850A+C(I)*DPF107(3,SEZA,I)
       P850B=P850B+C(I)*DPF107(3,SEZB,I)
       P1400A=P1400A+C(I)*DPF107(4,SEZA,I)
       P1400B=P1400B+C(I)*DPF107(4,SEZB,I)
       P2000A=P2000A+C(I)*DPF107(5,SEZA,I)
40     P2000B=P2000B+C(I)*DPF107(5,SEZB,I)
      P350A=10**P350A
      P350B=10**P350B
      P550A=10**P550A
      P550B=10**P550B
      P850A=10**P850A
      P850B=10**P850B
      P1400A=10**P1400A
      P1400B=10**P1400B
      P2000A=10**P2000A
      P2000B=10**P2000B
C     model errTe
      E350A=0.0
      E350B=0.0
      E550A=0.0
      E550B=0.0
      E850A=0.0
      E850B=0.0
      E1400A=0.0
      E1400B=0.0
      E2000A=0.0
      E2000B=0.0                
      DO 50 I=1,81
       E350A=E350A+C(I)*DERRTE(1,SEZA,I)
       E350B=E350B+C(I)*DERRTE(1,SEZB,I)
       E550A=E550A+C(I)*DERRTE(2,SEZA,I)
       E550B=E550B+C(I)*DERRTE(2,SEZB,I)
       E850A=E850A+C(I)*DERRTE(3,SEZA,I)
       E850B=E850B+C(I)*DERRTE(3,SEZB,I)
       E1400A=E1400A+C(I)*DERRTE(4,SEZA,I)
       E1400B=E1400B+C(I)*DERRTE(4,SEZB,I)
       E2000A=E2000A+C(I)*DERRTE(5,SEZA,I)
50     E2000B=E2000B+C(I)*DERRTE(5,SEZB,I)
      E350A=10**E350A
      E350B=10**E350B
      E550A=10**E550A
      E550B=10**E550B
      E850A=10**E850A
      E850B=10**E850B
      E1400A=10**E1400A
      E1400B=10**E1400B
      E2000A=10**E2000A
      E2000B=10**E2000B 
C
      IF (PF107Y .EQ. 1) THEN
       CALL TPCORR(INVDP,MLT,DDD,PF107,
     &             P350A,P350B,P550A,P550B,P850A,P850B,
     &             P1400A,P1400B,P2000A,P2000B, 
     &             TP350A,TP350B,TP550A,TP550B,TP850A,TP850B,
     &             TP140A,TP140B,TP200A,TP200B) 
       T350A=T350A+TP350A
       T350B=T350B+TP350B
       T550A=T550A+TP550A
       T550B=T550B+TP550B
       T850A=T850A+TP850A
       T850B=T850B+TP850B
       T1400A=T1400A+TP140A
       T1400B=T1400B+TP140B
       T2000A=T2000A+TP200A
       T2000B=T2000B+TP200B
      END IF 
C     Te
      IF (FUN .EQ. 0) THEN
       SEZDAY=(DDDB-DDDA)
       XDAY=DDDD-DDDA
       T350=(T350B-T350A)*SIN(DPI/2.0*XDAY/SEZDAY)+T350A
       T550=(T550B-T550A)*SIN(DPI/2.0*XDAY/SEZDAY)+T550A
       T850=(T850B-T850A)*SIN(DPI/2.0*XDAY/SEZDAY)+T850A
       T1400=(T1400B-T1400A)*SIN(DPI/2.0*XDAY/SEZDAY)+T1400A
       T2000=(T2000B-T2000A)*SIN(DPI/2.0*XDAY/SEZDAY)+T2000A
      ELSE
       SEZDAY=(DDDB-DDDA)
       XDAY=DDDD-DDDA
       T350=(T350A-T350B)*COS(DPI/2.0*XDAY/SEZDAY)+T350B
       T550=(T550A-T550B)*COS(DPI/2.0*XDAY/SEZDAY)+T550B
       T850=(T850A-T850B)*COS(DPI/2.0*XDAY/SEZDAY)+T850B
       T1400=(T1400A-T1400B)*COS(DPI/2.0*XDAY/SEZDAY)+T1400B
       T2000=(T2000A-T2000B)*COS(DPI/2.0*XDAY/SEZDAY)+T2000B
      END IF
C     error Te
      IF (FUN .EQ. 0) THEN
       SEZDAY=(DDDB-DDDA)
       XDAY=DDDD-DDDA
       E350=(E350B-E350A)*SIN(DPI/2.0*XDAY/SEZDAY)+E350A
       E550=(E550B-E550A)*SIN(DPI/2.0*XDAY/SEZDAY)+E550A
       E850=(E850B-E850A)*SIN(DPI/2.0*XDAY/SEZDAY)+E850A
       E1400=(E1400B-E1400A)*SIN(DPI/2.0*XDAY/SEZDAY)+E1400A
       E2000=(E2000B-E2000A)*SIN(DPI/2.0*XDAY/SEZDAY)+E2000A
      ELSE
       SEZDAY=(DDDB-DDDA)
       XDAY=DDDD-DDDA
       E350=(E350A-E350B)*COS(DPI/2.0*XDAY/SEZDAY)+E350B
       E550=(E550A-E550B)*COS(DPI/2.0*XDAY/SEZDAY)+E550B
       E850=(E850A-E850B)*COS(DPI/2.0*XDAY/SEZDAY)+E850B
       E1400=(E1400A-E1400B)*COS(DPI/2.0*XDAY/SEZDAY)+E1400B
       E2000=(E2000A-E2000B)*COS(DPI/2.0*XDAY/SEZDAY)+E2000B
      END IF      
C ////////////////////////////////////////////////////////
C     Te linear interpolation for altitude
      IF (ALT .LT. 550) THEN
       TE=(T550-T350)/200.0*(ALT-350)+T350
       SIGTE=(E550-E350)/200.0*(ALT-350)+E350    
      END IF
      IF ((ALT .GE. 550) .AND. (ALT .LT. 850)) THEN       
       TE=(T850-T550)/300.0*(ALT-550)+T550
       SIGTE=(E850-E550)/300.0*(ALT-550)+E550    
      END IF
      IF ((ALT .GE. 850) .AND. (ALT .LT. 1400)) THEN       
       TE=(T1400-T850)/550.0*(ALT-850)+T850
       SIGTE=(E1400-E850)/550.0*(ALT-850)+E850    
      END IF
      IF (ALT .GE. 1400) THEN       
       TE=(T2000-T1400)/600.0*(ALT-1400)+T1400
       SIGTE=(E2000-E1400)/600.0*(ALT-1400)+E1400    
      END IF 
      
      INVDIP=INVDP
           
      RETURN
      END
C
C
       REAL FUNCTION INTERP(N,L,V,X,XOUT)
C---------------------------------------------------------------------------------
        INTEGER N,L,S,S0,I
        REAL V(N),X(N),XOUT,Y2(N),YOUT,X0(4),V0(4)
        REAL XA,XB,XC,VA,VB,VC
        CALL locate(X,N,XOUT,S)
        IF (L .EQ. 0) THEN   
C       Spline interpolation (L=0)       
         IF (S .LT. 2) S=2
         IF (S .GT. (N-2)) S=N-2    
         S0=S-1 
         DO 10 I=1,4
          X0(I)=X(S0+I-1)
10        V0(I)=V(S0+I-1)
         CALL SPLINE(X0,V0,4,1e30,1e30,Y2)              
         CALL SPLINT(X0,V0,Y2,4,XOUT,YOUT)
        END IF 
        IF (L .EQ. 1) THEN
C       Linear interpolation (L=1)       
         IF ((S .GE.1) .AND. (S .LT. N)) THEN
          YOUT=(V(S+1)-V(S))/(X(S+1)-X(S))*(XOUT-X(S))+V(S)
         END IF
         IF (S .EQ. 0) THEN 
          YOUT=(V(2)-V(1))/(X(2)-X(1))*(XOUT-X(1))+V(1)
         END IF
         IF (S .EQ. N) THEN 
          YOUT=(V(N)-V(N-1))/(X(N)-X(N-1))*(XOUT-X(N))+V(N)
         END IF  
        END IF        
         IF (L .EQ. 2) THEN
C       Quadratic interpolation (L=2)    
         IF (S .LT. 2) S=2
         IF (S .GT. (N-1)) S=N-1    
          XA=X(S-1)
          XB=X(S)
          XC=X(S+1)
          VA=V(S-1)
          VB=V(S)
          VC=V(S+1)
          YOUT=VA*(XOUT-XB)*(XOUT-XC)/((XA-XB)*(XA-XC))+
     &         VB*(XOUT-XA)*(XOUT-XC)/((XB-XA)*(XB-XC))+
     &         VC*(XOUT-XA)*(XOUT-XB)/((XC-XA)*(XC-XB))        
         END IF         
        INTERP=YOUT
        RETURN
       END
C
C
      SUBROUTINE KODERR(MIRREQ,DOUT)
C------------------------------------------------------------------------------------
C coefficients - error model part
C------------------------------------------------------------------------------------
      REAL DOUT(5,3,81)
      INTEGER MIRREQ(81),I,J,K
      REAL DERRTE(5,3,81)
C     350km equinox
      DATA (DERRTE(1,1,J),J=1,81)/ 2.1178E+00, 1.3114E-07, 2.6148E-01,
     &                            -2.8483E-07,-8.3770E-02,-2.3135E-08,
     &                             5.7835E-02,-9.9532E-08,-2.0735E-02,
     &                            -1.9633E-01, 4.4717E-08, 5.0244E-02,
     &                            -9.1274E-09, 4.8611E-02, 2.5790E-08,
     &                             1.2098E-02,-3.3205E-08, 6.4414E-02,
     &                             2.4963E-08, 2.6340E-03,-3.3202E-08,
     &                             1.0058E-02, 1.0816E-09,-1.3928E-02,
     &                             1.4163E-08,-1.5039E-01, 1.0533E-08,
     &                            -4.5060E-03, 8.3297E-09, 2.4488E-03,
     &                            -3.7760E-09,-1.8094E-03, 9.6798E-02,
     &                            -1.8772E-09, 1.0867E-02,-1.9474E-09,
     &                            -6.8032E-03, 3.3603E-09, 7.8251E-03,
     &                             1.6945E-01,-1.1341E-08,-2.0581E-02,
     &                            -1.9111E-09,-5.8315E-03,-5.3837E-11,
     &                            -3.5216E-02,-3.3668E-09, 1.1079E-03,
     &                             8.5225E-09, 4.8433E-03,-2.6840E-09,
     &                             9.9141E-02,-1.0567E-08, 1.8306E-03,
     &                             1.9576E-10, 1.5468E-03, 1.6529E-01,
     &                            -8.3743E-09,-1.0179E-02,-1.0826E-09,
     &                            -1.4207E-03,-4.8573E-02, 1.1394E-08,
     &                             1.4470E-02, 2.6626E-10, 9.4379E-02,
     &                            -1.8262E-08,-2.0741E-03, 1.3797E-09,
     &                             5.9985E-03, 1.8945E-10, 3.4170E-04,
     &                             6.8581E-02, 5.7838E-09, 6.9529E-03,
     &                             1.5896E-02,-5.0075E-09, 2.8420E-02,
     &                            -4.8417E-09, 2.9191E-02, 4.8866E-02/
C     350km June solstice
      DATA (DERRTE(1,2,J),J=1,81)/ 2.1020E+00, 6.4074E-02, 8.8979E-02,
     &                            -1.6295E-01,-6.8680E-02,-4.9903E-02,
     &                             2.7634E-02, 5.3483E-02,-3.6534E-02,
     &                            -2.3004E-01, 5.3630E-02, 9.0376E-02,
     &                             2.6442E-02, 3.7594E-02, 1.7345E-02,
     &                            -4.3779E-02,-4.7004E-03, 9.0849E-02,
     &                             2.6102E-02,-3.2626E-02,-1.1288E-03,
     &                             2.4115E-02,-1.8059E-03,-3.9570E-03,
     &                             5.6676E-03,-6.3478E-02, 1.5132E-02,
     &                             3.1012E-02,-1.7904E-04,-1.7065E-04,
     &                            -6.6273E-03, 5.2600E-03,-1.5202E-01,
     &                             2.2061E-02, 1.0633E-02,-2.8630E-03,
     &                            -9.6908E-03,-1.3059E-03, 2.4577E-03,
     &                             9.1133E-02, 1.0057E-02, 7.1922E-04,
     &                             2.9878E-03, 7.7424E-04,-4.4353E-03,
     &                            -7.6271E-03, 1.4035E-02,-1.3870E-02,
     &                             1.6514E-03, 1.3885E-03, 7.7953E-04,
     &                             1.1731E-01, 2.3511E-02,-8.2554E-04,
     &                            -2.4741E-03, 1.5413E-04, 3.3407E-02,
     &                            -5.2112E-04,-6.3114E-03, 5.3990E-03,
     &                            -9.1618E-04, 1.1605E-02, 1.2671E-03,
     &                            -2.1350E-03, 2.6276E-03, 8.2348E-02,
     &                            -2.2401E-06,-1.4076E-03, 1.6459E-03,
     &                            -4.3340E-02, 1.4798E-02, 4.6935E-03,
     &                             3.0577E-02, 6.6725E-03, 9.0547E-03,
     &                            -1.5337E-02,-4.3676E-03,-6.5343E-02,
     &                             1.1111E-02, 2.3191E-03,-3.1603E-02/
C     550km equinox
      DATA (DERRTE(2,1,J),J=1,81)/ 2.0812E+00, 2.4308E-07, 5.6584E-01,
     &                             1.9086E-08,-7.6095E-02,-5.1032E-07,
     &                            -1.7049E-01,-1.0472E-08, 6.7164E-02,
     &                            -2.2595E-01, 4.2122E-09, 6.1906E-03,
     &                             3.9518E-08, 4.5763E-02,-2.2074E-08,
     &                            -1.1397E-02,-1.0607E-08, 1.6385E-01,
     &                             4.4170E-08,-2.2410E-02,-1.5556E-08,
     &                             1.2078E-02,-2.9345E-08,-2.3789E-03,
     &                             5.0526E-08,-2.4343E-01, 1.6201E-08,
     &                             1.9976E-02,-2.2397E-11, 4.6519E-03,
     &                            -1.1324E-09,-2.2603E-03,-1.2091E-01,
     &                             6.5523E-09, 2.3472E-02, 1.2557E-10,
     &                             3.5617E-03, 2.1655E-09,-1.3878E-03,
     &                            -2.4845E-02,-3.1878E-09,-1.3496E-02,
     &                             3.2648E-09,-8.7088E-04,-3.2223E-09,
     &                            -8.9965E-02, 1.6719E-09, 5.9083E-03,
     &                             7.0153E-10, 6.3045E-06,-1.7459E-09,
     &                             1.1872E-01, 3.4154E-09, 8.3146E-03,
     &                             7.3305E-11,-1.5977E-04, 2.3161E-02,
     &                             1.6234E-09,-9.7428E-03,-3.1068E-10,
     &                             1.1073E-03, 6.8273E-02, 7.3145E-09,
     &                             4.3073E-03, 1.3123E-10, 1.7105E-02,
     &                             1.9378E-09,-1.6849E-03, 3.9189E-10,
     &                             1.4436E-02,-7.5937E-10, 5.5921E-04,
     &                             3.0434E-02,-5.3588E-10, 1.7285E-03,
     &                            -4.6172E-02, 1.0424E-09,-3.3645E-02,
     &                             4.5595E-10,-1.0524E-02,-6.6919E-02/
C     550km June solstice
      DATA (DERRTE(2,2,J),J=1,81)/ 2.0903E+00, 4.9465E-02, 4.1834E-01,
     &                            -3.0721E-02,-1.0157E-01,-1.5666E-01,
     &                            -7.6741E-02, 1.2230E-01, 3.1210E-02,
     &                            -1.5269E-01, 4.5530E-02,-2.1376E-03,
     &                            -1.8327E-02, 2.3483E-02, 3.1872E-02,
     &                            -2.3882E-02,-1.4407E-02, 4.4489E-02,
     &                             1.4618E-02,-2.3584E-04,-6.5495E-03,
     &                             1.6602E-02,-7.5274E-03,-1.3698E-02,
     &                            -8.2527E-03,-2.1590E-01, 2.5813E-02,
     &                             5.0590E-02, 2.1022E-03, 1.0141E-02,
     &                            -5.1961E-03, 5.4961E-03, 9.7779E-03,
     &                             1.5319E-02,-1.0422E-02,-9.1045E-03,
     &                            -8.1248E-03,-5.5846E-03, 1.2857E-03,
     &                            -3.1836E-02,-4.4886E-03,-1.8187E-02,
     &                            -1.0253E-03,-4.4267E-04,-1.6932E-03,
     &                            -8.0883E-02, 2.2869E-02, 3.6385E-03,
     &                             6.0430E-03, 1.8587E-03,-1.5704E-03,
     &                             5.2741E-02, 2.3327E-03,-8.1670E-03,
     &                             2.3873E-03,-4.6541E-04,-5.5659E-02,
     &                             1.0942E-02,-1.0447E-02,-4.5905E-04,
     &                             2.9831E-04, 9.2199E-02,-8.8334E-03,
     &                             4.9657E-03,-1.2876E-03, 1.5785E-02,
     &                            -1.1838E-02,-2.4461E-04, 8.0350E-04,
     &                             8.4769E-03, 7.2212E-03, 1.2089E-03,
     &                             4.0519E-02,-8.4599E-03, 1.1588E-03,
     &                            -3.3485E-02, 5.8463E-03,-2.8652E-03,
     &                             8.8555E-03,-2.6834E-02,-5.0552E-02/
C     850km equinox
      DATA (DERRTE(3,1,J),J=1,81)/ 2.1741E+00,-3.0713E-07, 7.0198E-02,
     &                             4.9868E-07,-1.4427E-01, 2.8235E-08,
     &                             2.9561E-02,-8.1150E-07, 1.8500E-03,
     &                             8.5103E-02, 3.3995E-09, 3.7822E-02,
     &                             4.4995E-08, 3.8160E-02, 6.1687E-09,
     &                            -7.6820E-03,-3.5958E-08, 5.4672E-02,
     &                             1.5964E-10, 1.6325E-02,-2.9809E-08,
     &                            -2.8099E-03, 2.0887E-08,-1.6249E-03,
     &                             6.0045E-09,-4.3312E-02, 1.2503E-08,
     &                             3.9802E-02, 1.6466E-08, 1.4980E-02,
     &                            -7.7490E-09, 4.3219E-03, 1.0214E-01,
     &                            -1.6204E-08, 8.0035E-03, 3.9310E-09,
     &                            -2.1047E-04,-2.0211E-10, 2.3744E-03,
     &                            -1.0833E-01,-2.6971E-09,-1.8076E-02,
     &                             3.9598E-09, 8.3868E-05,-1.0948E-10,
     &                             7.6547E-02, 1.1306E-09,-1.5295E-02,
     &                            -8.6391E-10,-9.3252E-04,-9.1719E-12,
     &                             6.3563E-02,-1.2694E-08,-6.3280E-03,
     &                             5.2632E-10,-7.9613E-04,-1.1528E-02,
     &                             1.0676E-08, 1.0950E-02,-5.2360E-09,
     &                            -1.0426E-03,-8.4653E-03, 2.7643E-09,
     &                            -1.0339E-03, 1.0439E-09, 8.2576E-02,
     &                            -1.0238E-08,-5.2244E-03,-2.8944E-10,
     &                             1.7164E-02,-3.7405E-09, 1.3255E-03,
     &                            -1.6994E-02, 4.5616E-09, 1.3458E-03,
     &                             2.5439E-02,-5.2198E-09, 3.0276E-02,
     &                            -3.3313E-09,-5.7659E-02, 2.8990E-02/
C     850km June solstice
      DATA (DERRTE(3,2,J),J=1,81)/ 2.1685E+00,-9.1071E-02,-1.2046E-01,
     &                            -2.0084E-01,-7.3145E-02,-1.2735E-01,
     &                            -5.5899E-02, 5.2264E-02, 2.3232E-02,
     &                             2.8355E-02,-3.9715E-02, 6.7638E-02,
     &                             2.3396E-02,-8.2230E-03, 1.1384E-02,
     &                            -1.6781E-02,-1.7576E-03,-8.0253E-02,
     &                            -3.9270E-02,-1.0264E-02,-1.7191E-03,
     &                            -4.4668E-03,-1.2252E-04,-7.9374E-03,
     &                             2.3770E-03, 5.3901E-02,-1.9019E-03,
     &                             3.2748E-02, 8.1218E-03, 3.1385E-03,
     &                             6.5387E-03, 1.9147E-03, 9.6127E-02,
     &                            -2.7965E-02,-7.7649E-03, 3.7879E-03,
     &                            -4.0496E-03, 6.9544E-03,-5.2783E-04,
     &                            -1.5082E-01, 1.2699E-02,-1.2565E-03,
     &                             6.7728E-03, 3.6761E-03,-1.7899E-03,
     &                            -6.4801E-02,-1.7940E-02,-1.3257E-02,
     &                             5.2332E-03,-6.6084E-03,-2.9180E-03,
     &                             3.0122E-02, 2.1713E-02,-8.3407E-03,
     &                             1.1084E-03, 9.2082E-04,-8.5997E-02,
     &                            -1.4387E-02, 9.3689E-03, 6.2303E-04,
     &                             3.5905E-04, 1.7038E-02,-6.3187E-03,
     &                             3.6498E-03, 5.2909E-04,-7.0981E-03,
     &                             5.5184E-03, 6.6041E-04, 3.3276E-04,
     &                             6.9458E-02,-4.9348E-03, 6.0775E-03,
     &                            -5.0934E-02,-4.4936E-03, 1.3723E-03,
     &                            -8.1570E-03, 1.3748E-02, 7.0236E-03,
     &                             7.9914E-04,-4.6833E-02,-1.6385E-02/
C     1400km equinox
      DATA (DERRTE(4,1,J),J=1,81)/ 2.1822E+00, 3.3549E-08, 1.2164E-01,
     &                             2.4764E-07,-8.7502E-02,-4.0691E-07,
     &                             2.3355E-03,-2.8158E-07,-1.3242E-03,
     &                             1.4782E-02, 1.5416E-08, 8.8183E-02,
     &                             2.7046E-08,-2.2398E-02,-4.4919E-09,
     &                            -7.8416E-03,-1.7223E-08,-1.3115E-02,
     &                             2.6125E-08, 5.3183E-02,-2.9019E-08,
     &                            -3.7068E-04, 5.4685E-10, 1.7274E-03,
     &                             1.8592E-08,-6.2685E-02,-5.3055E-09,
     &                             3.1005E-02, 1.7579E-09, 2.4121E-03,
     &                             4.4444E-11, 4.1500E-03,-2.4028E-01,
     &                            -2.3704E-08, 8.7659E-03, 4.2252E-09,
     &                             4.8969E-03, 4.0202E-10, 1.9055E-03,
     &                            -1.6350E-01, 3.2289E-09, 1.4824E-03,
     &                            -3.6081E-09, 1.3517E-03, 2.4239E-09,
     &                            -1.6658E-01, 1.0138E-09, 1.3305E-02,
     &                             4.1669E-09, 1.8869E-03,-2.0831E-09,
     &                             2.0076E-02,-1.0135E-08, 2.2258E-03,
     &                            -1.2563E-09,-1.3067E-03,-1.8645E-01,
     &                             2.1705E-09,-7.8197E-03,-9.9205E-10,
     &                             9.4107E-04, 1.7436E-01,-2.0898E-09,
     &                             7.8659E-03, 2.7303E-10, 7.1875E-02,
     &                            -2.0799E-09,-3.6387E-03, 2.5309E-09,
     &                            -9.2385E-02, 1.1905E-09, 2.6655E-03,
     &                            -2.4122E-02, 6.0810E-10, 5.8355E-03,
     &                            -5.1302E-02,-7.9405E-09,-7.6189E-02,
     &                            -6.7824E-10,-3.1716E-02,-4.1055E-02/
C     1400km June solstice
      DATA (DERRTE(4,2,J),J=1,81)/ 2.1442E+00,-1.3864E-01, 1.7112E-01,
     &                            -2.1389E-01,-1.8186E-02,-1.0237E-01,
     &                             5.9302E-02, 6.0929E-02, 2.4732E-02,
     &                             7.8690E-02,-2.1925E-02, 1.2199E-01,
     &                            -2.8987E-02,-5.3213E-02,-4.8915E-03,
     &                             4.6986E-03,-1.1799E-02,-2.4101E-02,
     &                            -2.2908E-02, 2.3052E-02,-8.2149E-03,
     &                             3.1831E-03,-7.8876E-04, 9.0815E-03,
     &                             2.5709E-03,-1.9943E-01, 1.6644E-02,
     &                             3.6662E-02, 8.7738E-03,-2.2514E-03,
     &                            -3.0569E-04, 3.6344E-04, 5.6322E-02,
     &                             1.1496E-02, 1.8602E-02,-1.9180E-03,
     &                             6.0620E-03,-6.3761E-04, 3.2952E-03,
     &                            -1.4087E-01, 8.0978E-03, 4.2363E-03,
     &                             5.9284E-03, 9.9418E-04, 1.9006E-03,
     &                            -9.8914E-02, 8.8500E-03, 2.2600E-03,
     &                             4.5221E-05,-1.4855E-04, 6.4891E-04,
     &                             3.2002E-02,-1.8770E-02, 2.5766E-04,
     &                             6.8798E-04, 7.8152E-04,-5.2227E-02,
     &                             1.3226E-03,-1.6813E-03, 4.7524E-05,
     &                            -3.5793E-04, 6.7469E-02,-1.0750E-02,
     &                            -2.0372E-03, 2.9722E-04, 4.0765E-02,
     &                            -1.4128E-02,-6.7906E-03, 6.6127E-04,
     &                             2.0483E-02, 2.7251E-03,-3.2169E-04,
     &                             7.7473E-02,-6.8423E-03,-1.5793E-03,
     &                            -2.6085E-02, 8.8083E-03, 4.2452E-02,
     &                             8.2699E-04, 2.1088E-02, 1.7894E-02/
C     2000km equinox
      DATA (DERRTE(5,1,J),J=1,81)/ 2.4436E+00,-4.3843E-07, 1.5425E-01,
     &                            -5.5662E-08, 4.0385E-03, 8.1161E-07,
     &                             5.3553E-02,-6.9405E-07,-1.8544E-02,
     &                             1.7362E-01, 4.5975E-08,-6.9333E-03,
     &                             6.9511E-10,-1.9003E-02, 2.0966E-08,
     &                            -2.9413E-04,-4.3705E-08, 1.4120E-02,
     &                            -4.4319E-08, 2.6738E-02,-1.3382E-09,
     &                            -7.7713E-03,-2.1625E-08, 6.7646E-03,
     &                             2.9648E-08,-4.9498E-03,-6.9172E-09,
     &                             1.8036E-02, 8.0627E-09, 1.6950E-03,
     &                            -8.2902E-11, 2.3367E-04, 1.2521E-02,
     &                            -2.9426E-08,-1.0396E-02, 8.0042E-10,
     &                             3.0297E-03, 2.7689E-10,-9.5301E-04,
     &                            -2.3612E-02,-7.9111E-10,-7.6654E-03,
     &                             2.9311E-10,-2.8914E-04, 8.3715E-11,
     &                            -7.8202E-02, 6.6501E-09, 3.8289E-03,
     &                             4.0403E-09, 1.6187E-03,-5.7231E-10,
     &                             6.4240E-02,-3.3474E-09,-1.3766E-03,
     &                            -5.7988E-10, 6.2600E-04, 3.8725E-03,
     &                             1.5940E-08,-5.3759E-03,-2.1427E-09,
     &                            -8.9545E-04, 3.5896E-02,-4.3521E-09,
     &                             2.9789E-03, 1.0916E-09, 7.4931E-02,
     &                            -4.2211E-09,-8.2421E-03, 1.4316E-09,
     &                            -2.6314E-02,-4.4286E-10,-3.5706E-03,
     &                             6.3435E-02,-2.5381E-09,-2.2050E-03,
     &                            -6.6893E-02, 1.6210E-09, 4.0881E-02,
     &                            -5.1892E-09, 9.0158E-03, 2.3495E-02/
C     2000km June solstice
      DATA (DERRTE(5,2,J),J=1,81)/ 2.4223E+00,-1.6581E-01, 2.0883E-01,
     &                            -8.0645E-02,-3.9722E-02, 1.5954E-02,
     &                             3.0088E-02,-9.6053E-03, 7.4229E-03,
     &                             1.4053E-01, 7.9878E-03, 4.2976E-02,
     &                            -2.6882E-02,-1.1069E-02,-5.6813E-03,
     &                             3.1338E-03,-7.3863E-03,-9.6881E-03,
     &                            -4.0146E-02, 1.3978E-02,-1.6224E-02,
     &                             9.3162E-03,-1.0888E-02,-1.4472E-02,
     &                            -5.4611E-03,-5.0694E-02, 5.5719E-02,
     &                             6.4098E-03, 2.0679E-03,-3.4277E-03,
     &                             2.9826E-03, 2.2637E-03, 4.8907E-02,
     &                             8.6854E-03, 6.6625E-03, 3.3814E-03,
     &                             7.7458E-04, 2.9699E-03,-6.9636E-04,
     &                            -9.5779E-03, 1.6612E-03, 1.3018E-02,
     &                            -5.0237E-03,-1.6025E-03, 1.3429E-03,
     &                             2.9054E-02, 7.2855E-03,-1.3993E-02,
     &                            -5.2888E-03,-2.7511E-05, 2.6011E-04,
     &                             1.5534E-02,-1.6128E-02, 6.3239E-03,
     &                            -3.7710E-03, 7.4519E-04,-5.6335E-02,
     &                             1.4123E-02, 7.2038E-03,-8.2835E-04,
     &                            -3.9963E-04, 2.0466E-02, 1.4346E-02,
     &                            -3.4267E-03, 5.6522E-04,-6.6495E-02,
     &                            -6.4634E-04, 4.2662E-04, 1.3420E-03,
     &                             2.0465E-02, 6.0863E-04, 4.9604E-04,
     &                            -2.2755E-02,-7.0387E-03, 1.3109E-03,
     &                             3.6849E-02, 2.2601E-03, 2.2893E-02,
     &                            -1.1385E-02, 4.4417E-02,-6.5754E-03/
      DO 10 I=1,81
       DERRTE(1,3,I)=DERRTE(1,2,I)*MIRREQ(I)
       DERRTE(2,3,I)=DERRTE(2,2,I)*MIRREQ(I)
       DERRTE(3,3,I)=DERRTE(3,2,I)*MIRREQ(I)
       DERRTE(4,3,I)=DERRTE(4,2,I)*MIRREQ(I)
10     DERRTE(5,3,I)=DERRTE(5,2,I)*MIRREQ(I)
      DO 40 K=1,81
       DO 30 J=1,3
        DO 20 I=1,5 
         DOUT(I,J,K)=DERRTE(I,J,K)
20      CONTINUE
30     CONTINUE
40    CONTINUE
C////////////////////////////////////////////////////////////////////////////////////
      RETURN
      END
C
C
      SUBROUTINE KOEFD(MIRREQ,DOUT)
C------------------------------------------------------------------------------------
C    coefficients - main model part
C------------------------------------------------------------------------------------
      REAL DOUT(5,3,81)
      INTEGER MIRREQ(81),I,J,K
      REAL D(5,3,81)
C     350km equinox
      DATA (D(1,1,J),J=1,81)/ 3.1753E+00, 5.7689E-07, 1.9586E-01,
     &                       -6.8275E-07, 1.8921E-02, 3.6213E-07,
     &                       -1.8628E-02, 5.0105E-08,-8.5600E-03,
     &                       -1.6010E-01,-9.6605E-08,-1.9714E-02,
     &                        2.1637E-07, 4.8494E-03,-3.8884E-07,
     &                       -3.4967E-03, 2.2421E-07, 2.1263E-02,
     &                        9.3647E-08,-7.9466E-03, 1.0780E-07,
     &                        9.0129E-03,-1.2648E-07, 9.2017E-05,
     &                        5.0372E-08,-7.7926E-02,-5.7641E-08,
     &                       -9.1840E-03,-2.8697E-08, 8.4038E-04,
     &                       -3.3483E-09, 9.8633E-04,-3.6580E-02,
     &                        2.6765E-07, 9.6558E-03, 3.1708E-09,
     &                       -3.7801E-04,-1.4056E-08,-7.5288E-04,
     &                        5.8945E-02,-6.4803E-08, 6.6372E-03,
     &                       -5.1910E-08, 6.6456E-05, 5.3576E-09,
     &                       -5.9863E-02,-2.4909E-07,-1.1843E-03,
     &                        4.9602E-08, 2.4969E-04,-7.2411E-09,
     &                        6.1304E-02,-9.0388E-08, 1.4235E-03,
     &                        9.6584E-09, 4.8905E-04, 3.4504E-02,
     &                        3.9764E-08, 2.2885E-04,-2.0067E-08,
     &                       -2.3446E-04,-5.3749E-02, 5.5106E-08,
     &                       -2.5829E-05,-3.6483E-10, 3.9294E-02,
     &                       -1.1888E-07,-1.1182E-04, 8.2427E-09,
     &                       -3.4563E-02,-4.2396E-08,-6.5573E-04,
     &                       -2.7676E-02,-1.3017E-08, 5.8774E-04,
     &                        1.6090E-02, 4.6780E-08,-2.0595E-02,
     &                        2.5706E-08, 1.2541E-02, 1.7794E-02/
C     350km June solstice
      DATA (D(1,2,J),J=1,81)/ 3.1629E+00, 6.1302E-02, 2.0063E-01,
     &                        4.5285E-02,-4.1826E-02, 1.9260E-03,
     &                        2.4969E-03, 8.2077E-03,-2.8539E-02,
     &                       -1.6297E-01, 1.3558E-02,-2.4389E-02,
     &                       -3.1286E-03, 2.2763E-02, 2.1271E-03,
     &                       -6.5418E-03, 9.5738E-04, 1.6215E-02,
     &                        1.8507E-03,-1.0563E-03, 1.1742E-03,
     &                       -1.0133E-03,-1.9169E-03, 3.6940E-03,
     &                       -1.5715E-04,-5.6835E-02,-4.2398E-03,
     &                        6.2377E-03,-1.5667E-03, 4.1399E-04,
     &                       -2.1687E-03, 2.0367E-04,-1.4575E-02,
     &                        8.7408E-03,-1.5991E-03, 2.2254E-03,
     &                        1.3953E-03,-6.9929E-04,-1.0188E-03,
     &                        2.6039E-02,-7.4757E-04, 2.1545E-03,
     &                        1.0206E-03, 9.9516E-05,-5.1330E-04,
     &                       -4.3098E-03, 6.8083E-04, 1.7507E-03,
     &                        1.8198E-03,-5.2904E-04,-5.0768E-04,
     &                        3.1658E-02, 4.7564E-03,-9.2531E-04,
     &                        1.7889E-03,-3.7194E-05, 1.6896E-02,
     &                       -5.0829E-04,-8.1960E-04, 4.8322E-04,
     &                       -2.6520E-05,-1.0728E-02, 4.1619E-03,
     &                        3.8523E-04, 4.9239E-05, 1.1805E-02,
     &                       -2.9327E-03, 2.9543E-04, 2.9735E-04,
     &                       -9.7455E-03,-1.6474E-03, 1.0511E-03,
     &                       -1.7878E-02, 5.3427E-03, 1.7965E-05,
     &                        3.5806E-03, 6.6563E-04,-1.2012E-02,
     &                       -5.7183E-03, 2.4944E-03, 4.3333E-03/
C     550km equinox
      DATA (D(2,1,J),J=1,81)/ 3.2751E+00, 1.3545E-06, 1.7150E-01,
     &                       -7.6203E-07, 4.0429E-02, 5.0117E-07,
     &                       -2.2589E-02, 3.0332E-08,-2.0861E-02,
     &                       -1.4168E-01, 2.7793E-07,-1.8331E-02,
     &                        2.2669E-07, 4.6414E-03, 6.9158E-08,
     &                        1.3312E-03,-1.0690E-07, 1.7280E-02,
     &                        3.6689E-07,-1.0921E-03,-4.8798E-07,
     &                        2.0699E-03, 1.5270E-07, 4.0919E-03,
     &                       -4.0587E-08,-8.4276E-02,-1.8326E-07,
     &                       -7.7966E-03, 2.9950E-07, 1.3969E-03,
     &                        2.2649E-08, 1.0082E-03,-3.1392E-02,
     &                       -5.0450E-08, 1.7356E-03,-1.4810E-08,
     &                        2.4669E-03, 2.0049E-08,-1.6503E-04,
     &                        5.6435E-02,-2.2274E-07, 4.4061E-03,
     &                       -1.8825E-08, 3.1341E-04, 1.3638E-08,
     &                       -4.7964E-02,-3.1192E-07,-3.8772E-04,
     &                        2.2998E-08,-5.8125E-05, 1.5519E-09,
     &                        4.9958E-02, 7.1110E-08, 1.6467E-03,
     &                       -3.5644E-08, 1.8249E-04, 2.7835E-02,
     &                       -1.2413E-07,-1.5793E-03, 2.2114E-08,
     &                       -3.9754E-04,-2.2227E-02, 2.2094E-07,
     &                        2.8621E-04, 5.3520E-09, 1.9280E-02,
     &                        5.1741E-08, 3.2861E-04,-1.7477E-08,
     &                       -1.9664E-02,-6.0668E-08,-3.8234E-05,
     &                       -1.8001E-02, 1.3656E-07, 4.1525E-04,
     &                        7.9458E-03,-1.0954E-07,-1.2869E-02,
     &                       -6.1937E-08, 9.2866E-03, 5.3834E-03/
C     550km June solstice
      DATA (D(2,2,J),J=1,81)/ 3.2706E+00, 3.5202E-02, 1.6153E-01,
     &                        7.5517E-02,-1.2975E-02,-1.0202E-02,
     &                       -1.4402E-02, 8.0992E-03,-2.4350E-02,
     &                       -1.5168E-01, 2.1369E-02,-1.2387E-02,
     &                       -6.1981E-03, 9.1958E-03, 2.4677E-03,
     &                       -4.9595E-04,-9.1671E-04, 1.1726E-02,
     &                        1.3928E-04,-6.8449E-03, 3.9822E-03,
     &                        6.7074E-03,-2.3986E-03, 9.3211E-04,
     &                       -1.2521E-03,-7.7734E-02, 5.3710E-03,
     &                        3.3728E-04,-6.9879E-04, 1.9755E-03,
     &                       -2.3413E-04, 7.7060E-04,-2.7894E-02,
     &                        1.2731E-02, 1.9542E-03, 2.6609E-03,
     &                       -3.5695E-04,-1.1310E-03,-6.4754E-04,
     &                        4.2322E-02,-6.9529E-03,-7.8122E-04,
     &                        1.9300E-05,-6.5241E-04,-2.0369E-05,
     &                       -3.0919E-02, 1.4165E-03, 1.9683E-03,
     &                        3.6586E-05, 2.3327E-04, 7.7364E-05,
     &                        3.6784E-02, 2.2045E-03, 2.1286E-04,
     &                        6.5019E-04, 1.3530E-04, 1.2010E-02,
     &                       -1.1533E-03,-1.0365E-03, 5.2595E-04,
     &                       -1.9816E-04,-6.4129E-03,-2.9187E-04,
     &                        6.4120E-04,-4.3814E-04, 1.9807E-02,
     &                       -1.7670E-03,-5.2232E-04,-2.5243E-05,
     &                       -1.4894E-02, 7.9595E-04, 1.8406E-04,
     &                        9.1260E-04, 5.9450E-04, 3.9428E-04,
     &                        1.7175E-03,-7.0834E-04,-7.3342E-03,
     &                        9.3208E-04, 8.0461E-03, 1.7211E-03/
C     850km equinox
      DATA (D(3,1,J),J=1,81)/ 3.4090E+00,-8.7890E-07, 1.4192E-01,
     &                        2.2022E-06,-2.3026E-02, 2.8080E-07,
     &                       -2.0490E-02,-2.4588E-06, 8.0287E-03,
     &                       -1.7872E-01,-5.9833E-07,-3.1967E-03,
     &                        5.7087E-07, 5.6724E-03,-3.6116E-08,
     &                       -3.9539E-03,-1.1406E-07, 2.5913E-02,
     &                       -4.3955E-07,-5.9265E-03, 3.8624E-07,
     &                        4.3640E-03,-2.1767E-07,-1.6975E-03,
     &                        1.4715E-08,-6.7945E-02, 6.7140E-07,
     &                       -5.8461E-03, 1.1914E-07, 2.4496E-04,
     &                       -4.9921E-08,-3.8830E-04,-1.2488E-02,
     &                       -3.9557E-07, 9.4377E-04, 1.5100E-07,
     &                        4.1863E-04,-2.6552E-08,-9.4832E-04,
     &                        3.7187E-02, 7.8305E-07,-3.6193E-03,
     &                        1.9130E-07,-5.5284E-04,-7.0234E-08,
     &                       -1.1616E-02, 7.5494E-08, 1.4075E-03,
     &                        7.9285E-08, 3.1501E-04,-3.4879E-08,
     &                        3.8138E-02, 4.5561E-07,-5.2993E-06,
     &                        7.2506E-08,-4.6427E-04, 1.1952E-02,
     &                        1.9508E-07,-7.3901E-04, 1.6083E-08,
     &                        4.6821E-04, 2.2138E-03, 1.9836E-08,
     &                        9.1512E-05, 2.5883E-10, 1.6160E-02,
     &                        9.8914E-08, 2.2891E-04, 3.5686E-08,
     &                       -2.3945E-02,-9.7206E-08, 4.1304E-04,
     &                       -1.1323E-02,-1.0595E-07, 2.4607E-06,
     &                       -1.1919E-02,-1.2954E-07,-1.5160E-02,
     &                       -9.8974E-08, 9.1103E-03,-9.8102E-04/
C     850km June solstice
      DATA (D(3,2,J),J=1,81)/ 3.4298E+00, 2.6491E-02, 1.3842E-01,
     &                        4.0680E-02,-5.6606E-02,-5.7153E-03,
     &                       -1.4690E-02, 7.5677E-03, 3.2216E-04,
     &                       -1.8092E-01, 1.2032E-02, 8.3830E-03,
     &                       -1.6994E-03, 7.3073E-03, 4.6493E-03,
     &                       -5.9776E-03,-1.0911E-03, 1.8614E-02,
     &                       -6.0207E-03,-4.0287E-03,-2.4354E-03,
     &                        3.2668E-03, 2.0486E-03,-2.4568E-03,
     &                        1.7932E-04,-7.2842E-02, 1.2262E-03,
     &                       -2.4092E-03, 2.1826E-03, 1.1999E-03,
     &                       -5.7385E-04,-7.9356E-04,-4.8885E-03,
     &                        6.0397E-03, 1.5036E-03,-6.6850E-04,
     &                       -1.5557E-03, 7.0217E-04,-6.6930E-04,
     &                        2.7146E-02,-1.7554E-04,-5.7637E-03,
     &                       -4.8048E-05, 1.7230E-04, 2.7707E-04,
     &                       -1.8722E-02,-2.5534E-03, 2.7433E-03,
     &                       -9.4675E-04, 3.5885E-04, 1.6066E-04,
     &                        4.7142E-02,-3.8642E-03,-1.8094E-03,
     &                        3.7500E-04,-8.2361E-05, 3.4841E-04,
     &                       -2.1525E-03,-1.6265E-03,-5.6589E-04,
     &                        3.9127E-04, 2.3966E-02, 5.8480E-04,
     &                        1.1863E-05,-4.6920E-04,-3.2951E-03,
     &                       -1.3868E-03,-3.2869E-04,-4.1003E-05,
     &                       -8.5024E-03, 6.8306E-04, 1.2730E-03,
     &                        3.5803E-03,-6.7022E-05,-9.2460E-04,
     &                       -5.0650E-03, 2.3747E-03,-1.1690E-02,
     &                        2.7210E-04,-1.2308E-02,-2.4554E-03/
C     1400km equinox
      DATA (D(4,1,J),J=1,81)/ 3.4323E+00, 5.9891E-08, 1.4440E-01,
     &                       -2.5296E-06,-3.0483E-02, 8.7633E-07,
     &                       -2.0333E-02, 1.4639E-06, 8.3905E-03,
     &                       -2.4661E-01,-3.3357E-07, 2.3371E-02,
     &                       -2.6009E-07, 7.4049E-03,-4.2111E-07,
     &                       -4.0153E-03, 3.2564E-07,-3.3975E-02,
     &                       -3.1496E-07, 5.0539E-04,-2.1844E-07,
     &                        6.6046E-03, 3.7777E-07, 1.0884E-03,
     &                       -1.7690E-07,-9.7080E-02,-1.0183E-06,
     &                        3.2316E-04, 2.3879E-07, 1.0714E-03,
     &                       -1.5921E-07,-2.2353E-04,-7.8730E-02,
     &                       -4.9271E-07,-3.3524E-03,-9.7695E-08,
     &                        4.0369E-06, 9.1746E-08, 4.3697E-04,
     &                        5.4799E-02,-1.5613E-06, 2.0376E-05,
     &                        1.5198E-07, 5.2388E-05,-4.3276E-08,
     &                       -1.9427E-02,-4.8251E-07,-3.6941E-03,
     &                        1.0404E-08,-3.8792E-04, 1.9478E-08,
     &                        2.7545E-02,-7.9193E-07, 1.9158E-05,
     &                       -3.8582E-08, 4.3319E-06, 2.1105E-02,
     &                        5.6969E-08,-1.4642E-03, 1.7346E-08,
     &                       -2.0981E-04,-1.6940E-02,-3.2828E-07,
     &                        9.3221E-04,-5.8952E-08, 2.2700E-02,
     &                        8.6772E-07,-3.0788E-04,-3.1938E-08,
     &                       -3.5342E-02,-4.1651E-08,-4.6744E-04,
     &                        3.2506E-03, 9.8322E-07, 2.4373E-04,
     &                        1.0563E-02,-4.2559E-08,-1.4966E-02,
     &                        6.3242E-07, 3.0936E-02, 5.0785E-04/
C     1400km June solstice
      DATA (D(4,2,J),J=1,81)/ 3.4432E+00, 1.4031E-02, 1.3160E-01,
     &                        2.5090E-02,-4.2982E-02,-4.1350E-03,
     &                       -2.3576E-02, 4.5813E-04, 9.9480E-03,
     &                       -2.4132E-01, 3.9410E-03, 2.5328E-02,
     &                        2.5753E-03, 4.7467E-03, 2.4033E-03,
     &                       -4.5280E-03,-1.6768E-03,-2.4686E-02,
     &                        1.3680E-03, 4.5365E-03,-2.2588E-04,
     &                        4.7237E-03,-3.4808E-03, 1.9872E-03,
     &                        2.0111E-05,-1.0665E-01, 1.2940E-03,
     &                        5.1466E-04, 1.1414E-03,-2.6179E-04,
     &                        1.6622E-03,-2.9043E-04,-4.1638E-02,
     &                        2.5988E-03, 2.4967E-03, 1.1947E-03,
     &                        9.8563E-04,-1.1244E-03, 4.6249E-04,
     &                        2.2485E-02,-1.8404E-03,-2.7030E-03,
     &                       -5.6671E-04, 6.2775E-05, 2.3428E-04,
     &                       -1.9091E-02, 3.9310E-03,-1.3812E-03,
     &                        6.2219E-04, 4.0298E-04,-2.5167E-04,
     &                        3.6186E-02,-1.6894E-03,-9.2292E-04,
     &                       -3.5095E-04,-1.6027E-05, 5.0279E-03,
     &                        1.6496E-03,-1.2996E-03, 1.0942E-04,
     &                       -8.8813E-05,-9.8253E-04, 6.0205E-05,
     &                        3.1574E-04,-2.5650E-04, 1.6418E-02,
     &                       -7.4029E-04, 8.5019E-05, 7.1348E-05,
     &                       -8.8478E-03, 7.9330E-04, 4.9878E-04,
     &                        6.8566E-03,-1.3812E-03, 3.7630E-04,
     &                       -4.2406E-04,-4.0773E-05,-7.0264E-03,
     &                       -4.8997E-04, 3.2745E-03,-4.3985E-03/
C     2000km equinox
      DATA (D(5,1,J),J=1,81)/ 3.5267E+00, 3.1496E-06, 1.0072E-01,
     &                       -4.0719E-06,-2.9570E-02, 3.1967E-06,
     &                       -1.6759E-02,-2.2291E-06, 1.4529E-02,
     &                       -2.0303E-01, 2.6118E-06, 2.3639E-02,
     &                       -1.5963E-06, 2.8544E-03, 6.3910E-07,
     &                       -6.4201E-03,-1.4241E-07,-2.2074E-02,
     &                        1.3024E-07,-6.5225E-04,-7.5969E-09,
     &                        6.8229E-03,-1.1600E-08, 6.5430E-04,
     &                       -3.1299E-09,-8.9155E-02, 1.0420E-06,
     &                        2.2701E-04,-2.8037E-07, 7.6229E-05,
     &                        9.0789E-08,-8.3574E-04,-4.3006E-02,
     &                        9.5987E-08,-2.5486E-03, 1.2764E-08,
     &                        1.6634E-03,-9.8076E-09,-2.0750E-05,
     &                        4.1638E-02, 1.2831E-07, 3.2205E-03,
     &                       -7.2843E-08, 4.7908E-04, 1.0592E-09,
     &                       -8.8335E-03,-8.6869E-08,-2.5895E-03,
     &                        8.5979E-09, 8.1594E-04, 4.8263E-09,
     &                        4.4205E-02,-4.5836E-08, 2.4211E-03,
     &                       -1.7175E-08, 1.3454E-04, 2.9288E-02,
     &                       -6.4336E-08,-4.7940E-04,-1.7677E-08,
     &                        2.7104E-04, 3.8585E-03,-1.9495E-08,
     &                        9.5007E-04, 1.4881E-10, 3.5952E-02,
     &                        1.2439E-07,-2.6259E-04,-3.3433E-08,
     &                       -4.5002E-03,-9.7622E-08, 1.0627E-04,
     &                        1.7064E-02, 2.6949E-07, 3.5386E-05,
     &                        1.3508E-03,-9.3600E-08, 2.2163E-03,
     &                        1.8302E-07, 1.4533E-02, 6.8208E-03/
C     2000km June solstice
      DATA (D(5,2,J),J=1,81)/ 3.5455E+00, 2.3030E-02, 9.1028E-02,
     &                        5.6686E-03,-4.3555E-02, 7.5223E-03,
     &                       -5.3314E-03,-2.3664E-03, 9.2735E-03,
     &                       -1.7562E-01, 6.8722E-03, 2.6142E-02,
     &                       -7.5933E-03,-2.6850E-03, 1.5072E-03,
     &                       -3.8334E-03,-1.8455E-03,-8.3019E-03,
     &                       -3.1750E-03,-3.4582E-03,-3.0116E-03,
     &                        6.4951E-03,-1.9833E-04, 1.9620E-03,
     &                       -1.6688E-03,-9.1850E-02, 2.8064E-03,
     &                        1.1347E-03, 5.5757E-04,-1.7941E-03,
     &                        4.4038E-04,-7.6030E-04,-1.7519E-02,
     &                        5.5593E-03,-4.0485E-04, 1.0403E-04,
     &                        9.7795E-05,-4.4813E-04, 1.0662E-04,
     &                       -4.4894E-03,-4.3003E-03,-1.3180E-03,
     &                        1.2638E-03,-2.0937E-04,-1.0982E-04,
     &                       -2.8474E-02, 5.8161E-03,-5.9371E-04,
     &                        1.4721E-03,-1.8707E-04, 6.4902E-05,
     &                        6.2783E-03,-1.3854E-03, 2.9896E-04,
     &                       -7.7301E-04,-1.6992E-05,-2.8036E-02,
     &                        1.9110E-03,-1.7128E-04, 6.0156E-04,
     &                       -1.7487E-04, 2.9229E-03,-1.4161E-03,
     &                        5.1841E-04,-8.2856E-05,-1.0698E-02,
     &                       -4.2081E-03,-4.5016E-04, 5.1778E-07,
     &                        1.4779E-02, 5.7698E-04,-3.6797E-04,
     &                       -4.7711E-03,-1.6291E-03,-4.5695E-04,
     &                        1.1890E-02,-1.6669E-04,-5.5450E-03,
     &                       -1.0370E-03,-4.2745E-03, 1.8717E-03/
      DO 10 I=1,81
       D(1,3,I)=D(1,2,I)*MIRREQ(I)
       D(2,3,I)=D(2,2,I)*MIRREQ(I)
       D(3,3,I)=D(3,2,I)*MIRREQ(I)
       D(4,3,I)=D(4,2,I)*MIRREQ(I)
10     D(5,3,I)=D(5,2,I)*MIRREQ(I)
      DO 40 K=1,81
       DO 30 J=1,3
        DO 20 I=1,5 
         DOUT(I,J,K)=D(I,J,K)
20      CONTINUE
30     CONTINUE
40    CONTINUE
C////////////////////////////////////////////////////////////////////////////////////
      RETURN
      END
C
C
      SUBROUTINE KOF107(MIRREQ,DOUT)
C------------------------------------------------------------------------------------
C   coefficients - F107 model part
C------------------------------------------------------------------------------------
      REAL DOUT(5,3,81)
      INTEGER MIRREQ(81),I,J,K
      REAL DPF107(5,3,81)
C     350km equinox
      DATA (DPF107(1,1,J),J=1,81)/ 2.1179E+00, 2.6093E-07, 4.4282E-02,
     &                            -1.3876E-06, 2.7803E-02, 1.3986E-06,
     &                             2.9846E-02, 5.5440E-07, 1.3094E-02,
     &                            -2.5151E-02,-1.0646E-07,-3.7429E-03,
     &                             1.6883E-08,-8.7810E-05, 2.4570E-08,
     &                             5.3157E-04, 2.4498E-09,-7.4921E-03,
     &                             3.4559E-07, 2.7384E-03,-2.9301E-07,
     &                            -6.4160E-03,-2.9594E-08,-1.2326E-03,
     &                             2.0038E-07, 1.4916E-01, 8.9752E-08,
     &                             1.7673E-02,-3.3613E-08, 3.8401E-03,
     &                             5.1558E-08, 4.3661E-05,-1.3354E-01,
     &                            -4.2996E-08,-1.6325E-02,-5.9756E-08,
     &                            -7.8924E-04, 9.0574E-09, 9.4891E-05,
     &                            -4.6522E-02, 1.6467E-07, 3.0933E-03,
     &                            -2.4212E-08,-8.5264E-04,-2.8060E-09,
     &                             2.0465E-02, 3.9607E-08, 2.4059E-03,
     &                             2.3086E-08,-9.3343E-04, 2.1124E-08,
     &                             3.5658E-02,-4.9285E-08,-9.9046E-04,
     &                            -2.8411E-08, 4.9889E-04,-2.0836E-02,
     &                             2.7124E-07, 3.9610E-03,-5.1577E-08,
     &                             1.2194E-03, 2.0681E-02, 2.2079E-08,
     &                            -1.5148E-03,-1.4727E-08, 2.9282E-02,
     &                            -7.1036E-08,-1.0708E-03,-7.0356E-09,
     &                             3.5476E-02, 1.4627E-08, 1.6406E-03,
     &                             2.0173E-02, 1.1973E-07, 5.7493E-06,
     &                             9.8510E-04,-5.6858E-08,-6.9055E-03,
     &                             8.8015E-08,-5.2797E-03, 4.9785E-02/
C     350km June solstice
      DATA (DPF107(1,2,J),J=1,81)/ 2.1118E+00, 1.2667E-02, 7.8532E-02,
     &                             2.8013E-03, 4.6247E-03,-1.8002E-02,
     &                             1.7402E-02, 2.1645E-02, 9.7880E-04,
     &                             1.6106E-03, 1.5252E-02, 4.2830E-04,
     &                             2.2738E-03, 5.1962E-03,-1.2009E-03,
     &                            -2.4580E-03, 1.1339E-03, 2.4005E-03,
     &                            -9.1238E-03,-6.2242E-03, 7.1732E-04,
     &                             4.7369E-03,-2.0235E-03,-2.8728E-03,
     &                             1.5519E-03,-1.5735E-01, 5.1640E-03,
     &                            -1.0377E-02, 1.5376E-03,-1.2147E-04,
     &                             1.1153E-03,-4.0444E-04, 6.2152E-02,
     &                             4.3814E-03, 8.5715E-03,-1.1794E-03,
     &                            -1.4507E-03,-2.9136E-04, 4.8678E-05,
     &                            -2.5783E-02,-6.4352E-03, 2.8888E-03,
     &                            -2.3923E-04,-7.8468E-04, 1.0567E-04,
     &                             1.8717E-02,-3.7789E-03, 1.9293E-03,
     &                             2.7927E-03, 6.1046E-04,-6.4161E-04,
     &                             2.9512E-02, 1.8061E-03,-6.0273E-04,
     &                            -1.4189E-03,-1.2866E-04, 3.4173E-02,
     &                             7.0098E-03, 2.7884E-03, 1.3619E-03,
     &                             1.0141E-03, 4.4472E-02,-3.6668E-03,
     &                            -6.6634E-04,-9.1539E-04,-4.3310E-02,
     &                             4.6828E-03,-1.0213E-03, 9.6850E-04,
     &                            -8.9240E-04, 3.0694E-03, 1.1888E-03,
     &                            -3.2205E-02, 2.4068E-04, 1.5491E-03,
     &                            -1.3486E-03, 1.5805E-03, 1.2679E-02,
     &                             1.5963E-03,-9.2296E-03, 3.3470E-02/
C     550km equinox
      DATA (DPF107(2,1,J),J=1,81)/ 2.2596E+00,-1.0250E-06,-6.0469E-02,
     &                            -5.1717E-08, 2.4081E-02, 1.3561E-06,
     &                             2.2501E-02, 1.8975E-07,-1.9328E-03,
     &                             1.8631E-03,-6.9838E-08, 8.5643E-03,
     &                             2.2625E-07,-9.6354E-04, 1.4344E-07,
     &                            -9.9794E-04,-1.4888E-07,-2.1439E-02,
     &                             1.6757E-07, 2.2087E-02, 1.9753E-07,
     &                             3.3401E-03, 1.7510E-08,-4.4513E-03,
     &                            -1.9434E-07, 2.1312E-02, 4.1158E-07,
     &                             4.1293E-03, 1.2532E-07, 1.6378E-03,
     &                            -6.0591E-08, 3.8969E-04,-8.2128E-04,
     &                            -3.7332E-07,-9.0674E-03,-1.6562E-07,
     &                            -6.4580E-04,-2.4214E-08, 1.7187E-04,
     &                            -7.8323E-03,-3.4188E-09,-3.3704E-03,
     &                            -3.7064E-08,-8.1990E-04,-3.4103E-08,
     &                             2.8832E-04,-1.3091E-07, 5.9580E-04,
     &                             3.5007E-09,-2.4893E-05, 3.7407E-08,
     &                             9.8550E-03, 5.2340E-08, 2.4715E-03,
     &                            -2.2799E-09, 6.0276E-04,-2.3370E-03,
     &                             1.5897E-07,-7.8384E-04, 6.2186E-08,
     &                             1.5512E-04, 1.0034E-03, 6.2956E-08,
     &                             1.2935E-03, 1.1142E-08,-5.6848E-03,
     &                             5.7028E-08,-1.1159E-03, 5.2481E-09,
     &                            -6.4583E-04,-4.9024E-08,-8.7256E-05,
     &                            -3.6805E-03,-4.2162E-11, 5.9164E-04,
     &                            -1.1688E-03,-6.3884E-09,-3.4464E-03,
     &                            -2.9236E-08,-3.2531E-03, 7.3254E-04/
C     550km June solstice
      DATA (DPF107(2,2,J),J=1,81)/ 2.2443E+00,-7.5467E-03, 4.6154E-03,
     &                            -3.2706E-02,-9.4369E-03,-2.9842E-03,
     &                             9.5202E-03, 2.6462E-02, 8.1390E-03,
     &                             2.2865E-03,-8.1542E-03, 1.1328E-03,
     &                            -4.8105E-03, 1.5307E-03, 1.8632E-03,
     &                             1.1856E-03, 9.7583E-04,-1.4236E-03,
     &                             1.0665E-02, 1.1196E-02, 1.0141E-02,
     &                             2.7029E-03, 7.6283E-04,-2.2904E-04,
     &                            -7.3431E-04,-6.9053E-03, 2.8684E-04,
     &                             8.7857E-04, 8.9792E-04, 1.8426E-03,
     &                             1.0182E-05, 9.9144E-04,-6.8376E-03,
     &                             4.1889E-03,-8.3520E-04, 4.7692E-03,
     &                            -2.1779E-04,-3.2160E-04,-2.3404E-04,
     &                            -3.9832E-03,-8.0937E-04,-2.5324E-03,
     &                            -1.1029E-04,-3.9410E-04,-2.9127E-04,
     &                             1.4673E-02, 1.4644E-03, 4.4913E-03,
     &                             1.2082E-03, 1.1461E-03, 9.7184E-04,
     &                             8.7871E-03, 2.9808E-05, 2.7913E-03,
     &                             7.6192E-04, 9.4134E-04, 5.0748E-04,
     &                             3.2757E-03, 1.3459E-03, 7.4301E-04,
     &                             5.2271E-05,-1.1622E-03, 2.0424E-03,
     &                             7.4028E-04, 3.1044E-04,-7.8888E-04,
     &                             7.3294E-04, 4.5496E-04, 3.7781E-04,
     &                            -6.7759E-03,-1.3669E-03, 2.5073E-04,
     &                             7.8526E-03, 2.1084E-04,-6.8010E-06,
     &                            -1.2768E-04, 9.3358E-05, 5.0014E-03,
     &                             5.5091E-04, 6.4201E-04, 2.6821E-03/
C     850km equinox
      DATA (DPF107(3,1,J),J=1,81)/ 2.2430E+00, 6.3529E-07,-2.1093E-02,
     &                             2.3192E-06, 3.9848E-03,-4.0399E-06,
     &                            -2.0864E-03, 1.4147E-06,-2.7632E-03,
     &                             3.8821E-03, 4.3770E-07,-6.7410E-04,
     &                             3.7620E-07, 2.5406E-04,-5.5479E-07,
     &                             4.2678E-04, 3.7791E-07,-2.1396E-03,
     &                             4.3688E-07,-1.3975E-04, 5.5689E-08,
     &                            -7.7591E-04,-2.3032E-07,-7.6251E-06,
     &                             1.1410E-07, 3.5107E-02, 3.3463E-07,
     &                             5.8991E-03,-1.6889E-07, 6.4584E-04,
     &                             1.6902E-07, 1.7415E-04, 7.1402E-03,
     &                             4.6111E-07,-6.6733E-04, 1.5319E-07,
     &                            -2.8179E-04,-4.4107E-08, 1.8358E-04,
     &                             3.9796E-03,-4.7690E-08,-2.6286E-04,
     &                            -1.6931E-07, 5.3136E-05, 4.6929E-08,
     &                             4.6506E-06, 1.0615E-07,-2.2795E-04,
     &                            -1.9598E-08,-3.4348E-05,-1.8839E-08,
     &                            -6.1345E-03,-1.5260E-08, 5.7992E-04,
     &                            -7.6335E-08, 2.2449E-04, 4.1926E-03,
     &                             4.8940E-08, 1.2125E-03,-7.3629E-08,
     &                             3.8057E-04, 1.9487E-03,-1.6906E-09,
     &                             2.1031E-04,-8.3425E-09,-2.3397E-03,
     &                             3.7214E-08,-2.0578E-04,-7.7983E-08,
     &                             1.9869E-02,-1.5178E-07, 7.1451E-05,
     &                            -1.6288E-02,-2.8123E-08,-6.2559E-04,
     &                             1.6419E-04,-8.6871E-08,-3.8220E-03,
     &                            -3.4288E-08,-1.6062E-02,-2.3689E-03/
C     850km June solstice
      DATA (DPF107(3,2,J),J=1,81)/ 2.2234E+00, 8.7471E-03,-3.6271E-02,
     &                            -1.8205E-02, 1.4310E-02, 1.4137E-02,
     &                             1.1152E-03,-5.6408E-03,-2.3086E-04,
     &                             9.8994E-03, 3.6786E-04,-7.5900E-03,
     &                             3.6605E-03, 4.2347E-03, 9.9788E-04,
     &                            -1.0637E-03, 6.0786E-05, 3.3041E-03,
     &                             3.2600E-03,-3.9354E-03,-1.1438E-04,
     &                             2.8440E-03, 6.3044E-04,-1.1396E-03,
     &                             1.1203E-04, 3.9674E-02, 1.7704E-03,
     &                            -1.3756E-03,-2.0221E-04,-5.5309E-04,
     &                             2.5297E-04, 2.3301E-04, 2.1051E-02,
     &                             8.2900E-03,-9.2282E-04,-2.2149E-03,
     &                            -1.2940E-03, 5.0482E-04, 3.2355E-04,
     &                             3.3826E-03,-6.9424E-03, 4.4450E-04,
     &                             3.5048E-04,-2.8501E-04, 8.5714E-05,
     &                             5.0172E-05,-6.0011E-03,-1.9778E-03,
     &                             7.9132E-04, 4.8467E-04,-8.8136E-06,
     &                             5.5216E-03,-3.1276E-03,-1.6523E-03,
     &                             2.7189E-04, 1.2939E-04, 1.6626E-02,
     &                            -3.9681E-04, 3.1331E-03, 1.1602E-05,
     &                            -2.8408E-04, 1.3352E-02, 2.8371E-03,
     &                             4.9534E-04,-1.7470E-04,-6.1464E-04,
     &                            -7.3693E-03,-2.1311E-04, 8.8980E-05,
     &                             2.0388E-02, 2.5072E-04,-1.3785E-03,
     &                            -2.1501E-02,-2.0520E-03, 3.9847E-04,
     &                            -9.8123E-04, 3.0419E-03, 6.5642E-03,
     &                            -1.3463E-03,-4.1453E-03,-3.4043E-03/
C     1400km equinox
      DATA (DPF107(4,1,J),J=1,81)/ 2.0721E+00,-1.2810E-06,-2.2446E-02,
     &                             4.0113E-06, 1.5405E-02,-4.0985E-06,
     &                            -4.4873E-03, 1.6972E-06, 1.0420E-03,
     &                            -1.0762E-03, 4.3349E-08, 1.2188E-03,
     &                             2.3954E-08,-2.0416E-03,-2.3801E-08,
     &                             2.9302E-04,-2.3188E-08,-7.1780E-04,
     &                            -1.8559E-08,-2.6530E-03,-9.1704E-08,
     &                             1.9664E-03, 7.4905E-08, 9.7068E-06,
     &                             3.6472E-08, 9.6520E-03,-1.1466E-06,
     &                            -5.4198E-03, 4.4641E-07,-1.3638E-04,
     &                            -8.7144E-08, 4.4808E-05,-3.3499E-02,
     &                             1.2687E-06, 2.7432E-03,-6.2446E-07,
     &                            -4.9937E-04, 1.3601E-07,-1.1713E-04,
     &                             7.6744E-03,-2.2096E-08,-5.5947E-05,
     &                            -3.6655E-08,-2.3160E-04, 1.5485E-09,
     &                             1.9511E-03, 3.2565E-08,-5.0748E-04,
     &                            -2.2702E-08, 2.1519E-04, 4.2344E-09,
     &                            -2.4070E-02, 4.2536E-08,-1.0566E-03,
     &                            -4.8657E-08,-2.4020E-05,-5.3805E-02,
     &                             1.2672E-06, 1.9215E-03,-1.5258E-07,
     &                             1.4154E-04,-3.3868E-03, 4.4217E-08,
     &                            -6.1293E-04,-1.1966E-08,-7.1448E-03,
     &                            -1.3842E-08,-4.8955E-04, 1.9888E-08,
     &                            -5.3833E-02, 5.0588E-07, 1.7924E-04,
     &                            -1.0013E-02, 4.5334E-07, 4.6489E-04,
     &                            -1.7219E-05, 2.8226E-08, 1.2155E-03,
     &                            -5.5547E-08,-2.5776E-02, 3.7582E-02/
C     1400km June solstice
      DATA (DPF107(4,2,J),J=1,81)/ 2.0715E+00, 1.1060E-02,-7.2558E-03,
     &                             1.6759E-02, 4.4526E-03,-3.1949E-03,
     &                            -1.1849E-03, 1.9253E-03, 1.6289E-03,
     &                             2.1829E-04,-3.3721E-03, 2.6217E-03,
     &                             7.8261E-04, 5.2241E-04,-1.0364E-03,
     &                             1.4120E-05,-1.7664E-04, 2.6144E-03,
     &                             4.9622E-04,-1.1139E-03, 2.5938E-03,
     &                            -1.3858E-03,-4.7107E-04,-2.2748E-04,
     &                            -6.5263E-04, 1.5465E-04,-1.8961E-04,
     &                             1.9626E-03,-1.3251E-03,-3.5464E-04,
     &                            -3.6056E-04, 7.2928E-04,-3.4355E-02,
     &                            -5.4377E-04,-3.7358E-04,-5.4676E-04,
     &                            -6.0062E-04,-8.5823E-05,-4.2649E-04,
     &                            -2.5657E-03, 2.0810E-03,-5.0775E-04,
     &                            -6.9662E-05,-1.1269E-04,-2.9248E-05,
     &                            -4.2538E-04, 3.1001E-03, 9.3147E-05,
     &                            -2.0474E-04, 9.2787E-05, 1.5224E-04,
     &                             1.6699E-02,-1.0915E-03, 3.9374E-04,
     &                             6.2488E-05, 5.3697E-05, 4.2104E-03,
     &                             2.6857E-04,-5.6046E-04,-1.7589E-04,
     &                            -1.6645E-04,-1.0895E-03, 3.1274E-04,
     &                            -7.0907E-05,-1.8092E-04,-3.0431E-03,
     &                            -4.6234E-04,-2.6949E-05,-3.1150E-06,
     &                             4.6148E-03,-3.4373E-04, 9.9880E-04,
     &                            -2.2815E-03, 3.8619E-04,-8.9996E-05,
     &                            -2.2212E-05, 2.8260E-04, 3.3863E-03,
     &                            -6.4781E-04, 1.0682E-02,-5.7521E-03/
C     2000km equinox
      DATA (DPF107(5,1,J),J=1,81)/ 2.2225E+00,-1.8465E-07, 3.9702E-02,
     &                             7.4498E-08, 6.7294E-03,-1.5623E-07,
     &                             8.6264E-04, 4.9782E-07,-1.0457E-02,
     &                            -1.7603E-02,-4.2609E-07, 1.1957E-02,
     &                             2.8967E-07, 4.0676E-03,-4.8392E-08,
     &                            -2.0672E-03,-2.7825E-08, 2.3117E-02,
     &                            -8.1847E-07, 6.3865E-04, 1.2350E-07,
     &                            -3.1363E-03, 2.3973E-08, 7.4011E-04,
     &                             2.2201E-08, 1.5339E-02,-9.2354E-07,
     &                            -5.0780E-04, 2.4487E-07,-1.4002E-03,
     &                            -9.5180E-08,-2.2099E-04,-2.6709E-02,
     &                             2.0174E-07, 5.2254E-03,-6.4659E-08,
     &                            -4.5011E-04,-1.3226E-08,-5.0770E-04,
     &                             2.3246E-02, 4.0279E-07, 4.2692E-03,
     &                            -8.2578E-08, 2.8635E-04, 1.9589E-08,
     &                             6.1547E-04,-1.6112E-07,-3.4024E-03,
     &                            -3.3185E-08,-8.0475E-04, 6.1774E-09,
     &                             3.9677E-02, 1.8432E-07, 2.2421E-03,
     &                            -3.3228E-09,-2.2705E-05,-2.1308E-02,
     &                            -1.1246E-07, 2.6961E-03, 2.2048E-08,
     &                             1.6385E-05, 1.0901E-02, 4.6221E-07,
     &                            -3.8862E-06,-4.9407E-08, 1.3349E-02,
     &                             9.8339E-08,-5.5719E-04,-1.8928E-08,
     &                             3.9061E-02, 4.3038E-07, 1.0085E-03,
     &                             1.1895E-02, 4.6967E-07, 4.3363E-04,
     &                            -5.3628E-03, 4.2370E-07, 1.2707E-03,
     &                             1.0193E-07, 1.2395E-02, 1.0030E-02/
C     2000km June solstice
      DATA (DPF107(5,2,J),J=1,81)/ 2.2631E+00, 2.1866E-02, 4.6921E-03,
     &                            -3.3329E-02, 1.0108E-02,-1.2195E-02,
     &                             4.8366E-03, 2.1350E-02,-8.6916E-03,
     &                             3.4443E-03,-7.1057E-03,-1.4225E-03,
     &                            -9.9691E-03,-6.0693E-04,-2.9332E-03,
     &                            -1.0505E-03,-1.9665E-03, 2.1469E-02,
     &                            -1.9500E-02,-1.0986E-02, 5.2916E-03,
     &                             2.7388E-03,-4.6114E-05,-2.6742E-03,
     &                             2.0601E-06, 7.5298E-03,-1.1021E-02,
     &                            -3.8100E-03,-1.2864E-03,-1.0705E-03,
     &                            -4.8383E-04,-1.2329E-03, 3.4064E-02,
     &                            -1.6264E-02, 2.6406E-04, 1.1486E-03,
     &                             2.0181E-03, 4.7744E-04,-2.2118E-04,
     &                             8.5559E-03, 4.6544E-03,-3.0707E-03,
     &                             1.0697E-03, 7.4102E-05,-6.8598E-04,
     &                            -3.7443E-02,-4.1490E-03, 4.8394E-03,
     &                             2.8449E-04,-8.1610E-04, 4.1996E-04,
     &                            -2.3120E-02,-2.2430E-03,-1.6598E-03,
     &                             5.4323E-04,-3.9800E-04, 8.7947E-03,
     &                             5.9441E-03, 1.2565E-04, 8.0901E-04,
     &                            -7.3510E-04,-2.1540E-02,-2.6029E-04,
     &                             1.1616E-03,-1.5312E-04, 6.3625E-03,
     &                            -4.0037E-03,-1.5571E-03,-5.7964E-04,
     &                             8.1764E-03,-2.3762E-05,-1.2959E-04,
     &                             6.2361E-03,-4.9235E-03, 6.0991E-04,
     &                             1.6101E-03,-3.9088E-03,-1.3380E-02,
     &                            -4.2837E-04,-1.1667E-02, 3.9335E-03/
      DO 10 I=1,81
       DPF107(1,3,I)=DPF107(1,2,I)*MIRREQ(I)
       DPF107(2,3,I)=DPF107(2,2,I)*MIRREQ(I)
       DPF107(3,3,I)=DPF107(3,2,I)*MIRREQ(I)
       DPF107(4,3,I)=DPF107(4,2,I)*MIRREQ(I)
10     DPF107(5,3,I)=DPF107(5,2,I)*MIRREQ(I)
      DO 40 K=1,81
       DO 30 J=1,3
        DO 20 I=1,5 
         DOUT(I,J,K)=DPF107(I,J,K)
20      CONTINUE
30     CONTINUE
40    CONTINUE
C////////////////////////////////////////////////////////////////////////////////////
      RETURN
      END
C     
C
      SUBROUTINE locate(xx,n,x,j)
C------------------------------------------------------------------------------------
      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
C
C
        SUBROUTINE SPHARM_IK(C,L,M,COLAT,AZ)
C------------------------------------------------------------------------------------
C CALCULATES THE COEFFICIENTS OF THE SPHERICAL HARMONIC
C FROM IRI 95 MODEL
C NOTE: COEFFICIENTS CORRESPONDING TO COS, SIN SWAPPED!!!
C------------------------------------------------------------------------------------
      DIMENSION C(82)
      C(1)=1.
      K=2
      X=COS(COLAT)
      C(K)=X
      K=K+1
      DO 10 I=2,L
      C(K)=((2*I-1)*X*C(K-1)-(I-1)*C(K-2))/I
10    K=K+1
      Y=SIN(COLAT)
      DO 20 MT=1,M
      CAZ=COS(MT*AZ)
      SAZ=SIN(MT*AZ)
      C(K)=Y**MT
      K=K+1
      IF(MT.EQ.L) GOTO 16
      C(K)=C(K-1)*X*(2*MT+1)
      K=K+1
      IF((MT+1).EQ.L) GOTO 16
      DO 15 I=2+MT,L
      C(K)=((2*I-1)*X*C(K-1)-(I+MT-1)*C(K-2))/(I-MT)
15    K=K+1
16    N=L-MT+1
      DO 18 I=1,N
      C(K)=C(K-N)*SAZ
      C(K-N)=C(K-N)*CAZ
18    K=K+1
20    CONTINUE
      RETURN
      END
C     
C
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
C------------------------------------------------------------------------------------
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C
C
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
C------------------------------------------------------------------------------------
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END
C
C
       SUBROUTINE SWAPEL(N,A)
C------------------------------------------------------------------------------------
C      swaps elements of array
C------------------------------------------------------------------------------------
       INTEGER N,I  
       REAL A(N),AT(N)
       DO 10 I=1,N
10      AT(I)=A(I) 
       DO 20 I=0,N-1
20      A(I+1)=AT(N-I)                   
       RETURN
       END
C
C
       SUBROUTINE TEDIFI(F107IN,TEXN,TEDN,F107DF,TEDIF)
C------------------------------------------------------------------------------------
C      interpolation for solar activity        
C------------------------------------------------------------------------------------
       REAL F107IN,TEXN,TEDN,F107DF(3),TEDIF
       REAL T0DNXN(3),T0DN(2),TDNXN(2)
       REAL INTERP
       IF ((F107IN .GE. F107DF(1)) .AND. (F107IN .LE. F107DF(3))) THEN
        T0DNXN(1)=0.
        T0DNXN(2)=TEDN
        T0DNXN(3)=TEXN
        TEDIF=INTERP(3,2,T0DNXN,F107DF,F107IN)
       END IF    
       IF (F107IN .LT. F107DF(1)) THEN
        T0DN(1)=0.
        T0DN(2)=TEDN
        TEDIF=INTERP(2,1,T0DN,F107DF(1),F107IN)     
       END IF 
       IF (F107IN .GT. F107DF(3)) THEN
        TDNXN(1)=TEDN
        TDNXN(2)=TEXN
        TEDIF=INTERP(2,1,TDNXN,F107DF(2),F107IN)     
       END IF        
       RETURN
       END
C
C
       SUBROUTINE TPCAS(MLTRAD,PF107,PF107M,
     &                  XNDI,DNDI,PD,XNNI,DNNI,PN,TPASEA)
C------------------------------------------------------------------------------------
C      correction for season at fixed altitude
C------------------------------------------------------------------------------------
       REAL MLTRAD,PF107,PF107M,XNDI,DNDI,PD(3),XNNI,DNNI,PN(3),TPASEA
       REAL TA,TM,TDC,TNC
       CALL TEDIFI(PF107,XNDI,DNDI,PD,TA)
       CALL TEDIFI(PF107M,XNDI,DNDI,PD,TM)  
       TDC=TA-TM
        TDC=AMAX1(TDC,-1250.)
        TDC=AMIN1(TDC,1250.)
       CALL TEDIFI(PF107,XNNI,DNNI,PN,TA)
       CALL TEDIFI(PF107M,XNNI,DNNI,PN,TM)         
       TNC=TA-TM
        TNC=AMAX1(TNC,-1250.)
        TNC=AMIN1(TNC,1250.)       
C      harmonic interpolation for local time       
       TPASEA=(1.-COS(MLTRAD))/2.*(TDC-TNC)+TNC
       RETURN
       END

C
C
      SUBROUTINE   TPCORR(INVDIP,MLT,DDD,PF107,
     &             P350A,P350B,P550A,P550B,P850A,P850B,
     &             P1400A,P1400B,P2000A,P2000B, 
     &             TP350A,TP350B,TP550A,TP550B,TP850A,TP850B,
     &             TP140A,TP140B,TP200A,TP200B)
C------------------------------------------------------------------------------------
       REAL INVDIP,MLT,PF107
       INTEGER DDD
       REAL        P350A,P350B,P550A,P550B,P850A,P850B,
     &             P1400A,P1400B,P2000A,P2000B, 
     &             TP350A,TP350B,TP550A,TP550B,TP850A,TP850B,
     &             TP140A,TP140B,TP200A,TP200B
       REAL INTERP
       REAL MLTRAD
C      Constants    
       REAL INVDPQ(13)
C      PF107 Day Equinox       
       REAL P2DE(13,3),P1DE(13,3),P8DE(13,3),P5DE(13,3),P3DE(13,3)
C      Te max-min dif Day Equinox       
       REAL CXN2DE(13),CXN1DE(13),CXN8DE(13),CXN5DE(13),CXN3DE(13)
C      Te med-min dif Day Equinox       
       REAL CDN2DE(13),CDN1DE(13),CDN8DE(13),CDN5DE(13),CDN3DE(13)
C
C      PF107 Night Equinox       
       REAL P2NE(13,3),P1NE(13,3),P8NE(13,3),P5NE(13,3),P3NE(13,3)
C      Te max-min dif Night Equinox       
       REAL CXN2NE(13),CXN1NE(13),CXN8NE(13),CXN5NE(13),CXN3NE(13)
C      Te med-min dif Night Equinox       
       REAL CDN2NE(13),CDN1NE(13),CDN8NE(13),CDN5NE(13),CDN3NE(13)
C       
C
C      PF107 Day Solstice       
       REAL P2DS(13,3),P1DS(13,3),P8DS(13,3),P5DS(13,3),P3DS(13,3)
C      Te max-min dif Day Solstice       
       REAL CXN2DS(13),CXN1DS(13),CXN8DS(13),CXN5DS(13),CXN3DS(13)
C      Te med-min dif Day Solstice       
       REAL CDN2DS(13),CDN1DS(13),CDN8DS(13),CDN5DS(13),CDN3DS(13)
C
C      PF107 Night Solstice       
       REAL P2NS(13,3),P1NS(13,3),P8NS(13,3),P5NS(13,3),P3NS(13,3)
C      Te max-min dif Night Solstice       
       REAL CXN2NS(13),CXN1NS(13),CXN8NS(13),CXN5NS(13),CXN3NS(13)
C      Te med-min dif Night Solstice       
       REAL CDN2NS(13),CDN1NS(13),CDN8NS(13),CDN5NS(13),CDN3NS(13)
C      working variables
       REAL TXN2DE(13),TXN1DE(13),TXN8DE(13),TXN5DE(13),TXN3DE(13)
       REAL TDN2DE(13),TDN1DE(13),TDN8DE(13),TDN5DE(13),TDN3DE(13)
       REAL TXN2NE(13),TXN1NE(13),TXN8NE(13),TXN5NE(13),TXN3NE(13)
       REAL TDN2NE(13),TDN1NE(13),TDN8NE(13),TDN5NE(13),TDN3NE(13)
       REAL TXN2DS(13),TXN1DS(13),TXN8DS(13),TXN5DS(13),TXN3DS(13)
       REAL TDN2DS(13),TDN1DS(13),TDN8DS(13),TDN5DS(13),TDN3DS(13)
       REAL TXN2NS(13),TXN1NS(13),TXN8NS(13),TXN5NS(13),TXN3NS(13)
       REAL TDN2NS(13),TDN1NS(13),TDN8NS(13),TDN5NS(13),TDN3NS(13)      
C
C      Interpolated PF107
       REAL P2DEI(3),P1DEI(3),P8DEI(3),P5DEI(3),P3DEI(3)
       REAL P2NEI(3),P1NEI(3),P8NEI(3),P5NEI(3),P3NEI(3)
       REAL P2DSI(3),P1DSI(3),P8DSI(3),P5DSI(3),P3DSI(3)
       REAL P2NSI(3),P1NSI(3),P8NSI(3),P5NSI(3),P3NSI(3) 
C      Additional local and temporary variables
       REAL XN2DEI,XN1DEI,XN8DEI,XN5DEI,XN3DEI
       REAL XN2NEI,XN1NEI,XN8NEI,XN5NEI,XN3NEI
       REAL DN2DEI,DN1DEI,DN8DEI,DN5DEI,DN3DEI
       REAL DN2NEI,DN1NEI,DN8NEI,DN5NEI,DN3NEI
       REAL XN2DSI,XN1DSI,XN8DSI,XN5DSI,XN3DSI
       REAL XN2NSI,XN1NSI,XN8NSI,XN5NSI,XN3NSI
       REAL DN2DSI,DN1DSI,DN8DSI,DN5DSI,DN3DSI
       REAL DN2NSI,DN1NSI,DN8NSI,DN5NSI,DN3NSI
       REAL MLTTMP
       INTEGER I      
C
       DATA (INVDPQ(I),I=1,13) /-90.0,-75.0,-60.0,-45.0,-30.0,-15.0,
     &                      0.0, 15.0, 30.0, 45.0, 60.0, 75.0,90.0/
C      Equinox
       DATA ((P2DE(I,J),I=1,13),J=1,3) /
     &  125., 120., 115., 113., 118., 119., 118., 119., 118., 113.,
     &  115., 120., 125.,
     &  150., 151., 152., 158., 154., 156., 147., 156., 154., 158.,
     &  152., 151., 150.,
     &  181., 197., 213., 214., 205., 199., 197., 199., 205., 214.,
     &  213., 197., 181./
       DATA (CXN2DE(I),I=1,13) /  587.,  466.,  345.,  264.,  355.,
     &   335.,  229.,  335.,  355.,  264.,  345.,  466.,  587./
       DATA (CDN2DE(I),I=1,13) /  244.,   93.,  -59., -103.,   57.,
     &   170.,   16.,  170.,   57., -103.,  -59.,   93.,  244./
       DATA ((P1DE(I,J),I=1,13),J=1,3) /
     &  111., 110., 109., 109., 109., 110., 109., 110., 109., 109.,
     &  109., 110., 111.,
     &  135., 135., 134., 128., 131., 145., 148., 145., 131., 128.,
     &  134., 135., 135.,
     &  206., 201., 196., 200., 198., 197., 211., 197., 198., 200.,
     &  196., 201., 206./
       DATA (CXN1DE(I),I=1,13) /  929.,  590.,  250.,  375.,  126.,
     &   471.,  427.,  471.,  126.,  375.,  250.,  590.,  929./
       DATA (CDN1DE(I),I=1,13) /  -36.,   19.,   75.,   77.,   51.,
     &   100.,  105.,  100.,   51.,   77.,   75.,   19.,  -36./
       DATA ((P8DE(I,J),I=1,13),J=1,3) /
     &  125., 124., 124., 119., 102.,  89.,  98.,  89., 102., 119.,
     &  124., 124., 125.,
     &  147., 149., 151., 158., 164., 168., 166., 168., 164., 158.,
     &  151., 149., 147.,
     &  195., 194., 193., 195., 196., 191., 190., 191., 196., 195.,
     &  193., 194., 195./
       DATA (CXN8DE(I),I=1,13) /  159.,   75.,   -9., -151.,  511.,
     &   694.,  721.,  694.,  511., -151.,   -9.,   75.,  159./
       DATA (CDN8DE(I),I=1,13) /   53.,   47.,   42.,  -45.,  594.,
     &  1002.,  733., 1002.,  594.,  -45.,   42.,   47.,   53./
       DATA ((P5DE(I,J),I=1,13),J=1,3) /
     &  109., 102.,  95.,  93.,  95.,  92.,  86.,  92.,  95.,  93.,
     &   95., 102., 109.,
     &  161., 165., 169., 167., 160., 162., 164., 162., 160., 167.,
     &  169., 165., 161.,
     &  197., 200., 204., 203., 209., 214., 218., 214., 209., 203.,
     &  204., 200., 197./
       DATA (CXN5DE(I),I=1,13) /  703.,  472.,  242., -524., -241.,
     &   364.,  532.,  364., -241., -524.,  242.,  472.,  703./
       DATA (CDN5DE(I),I=1,13) /  533.,  416.,  300., -591., -398.,
     &   115.,  201.,  115., -398., -591.,  300.,  416.,  533./
       DATA ((P3DE(I,J),I=1,13),J=1,3) /
     &   70.,  78.,  86.,  89.,  90.,  92.,  89.,  92.,  90.,  89.,
     &   86.,  78.,  70.,
     &  159., 155., 152., 153., 152., 147., 143., 147., 152., 153.,
     &  152., 155., 159.,
     &  212., 212., 212., 208., 208., 207., 205., 207., 208., 208.,
     &  212., 212., 212./
       DATA (CXN3DE(I),I=1,13) /  277.,   77., -123., -364., -136.,
     &   317.,  411.,  317., -136., -364., -123.,   77.,  277./
       DATA (CDN3DE(I),I=1,13) /  -24.,  -15.,   -6.,  -74., -217.,
     &   188.,  269.,  188., -217.,  -74.,   -6.,  -15.,  -24./
       DATA ((P2NE(I,J),I=1,13),J=1,3) /
     &  113., 111., 109., 107., 103., 106., 102., 106., 103., 107.,
     &  109., 111., 113.,
     &  155., 157., 160., 158., 157., 146., 137., 146., 157., 158.,
     &  160., 157., 155.,
     &  188., 195., 202., 208., 209., 207., 201., 207., 209., 208.,
     &  202., 195., 188./
       DATA (CXN2NE(I),I=1,13) /   76.,  134.,  192.,  821.,  801.,
     &   583.,  551.,  583.,  801.,  821.,  192.,  134.,   76./
C       DATA (CDN2NE(I),I=1,13) /  -19.,  -83., -148.,  399.,  289.,
C     &   340.,  542.,  340.,  289.,  399., -148.,  -83.,  -19./  
C        equator corrected  
        DATA (CDN2NE(I),I=1,13) /  -19.,  -83., -148.,  399.,  289.,
     &   340.,  340.,  340.,  289.,  399., -148.,  -83.,  -19./
       DATA ((P1NE(I,J),I=1,13),J=1,3) /
     &  112., 111., 109., 108., 109., 110., 109., 110., 109., 108.,
     &  109., 111., 112.,
     &  132., 132., 133., 131., 134., 143., 143., 143., 134., 131.,
     &  133., 132., 132.,
     &  220., 211., 201., 195., 196., 204., 213., 204., 196., 195.,
     &  201., 211., 220./
       DATA (CXN1NE(I),I=1,13) /  806.,  689.,  572.,  607.,  540.,
     &   535.,  443.,  535.,  540.,  607.,  572.,  689.,  806./
       DATA (CDN1NE(I),I=1,13) /  124.,   23.,  -77.,   50.,  128.,
     &   120.,   79.,  120.,  128.,   50.,  -77.,   23.,  124./
       DATA ((P8NE(I,J),I=1,13),J=1,3) /
     &  124., 124., 124., 123., 119., 108., 110., 108., 119., 123.,
     &  124., 124., 124.,
     &  149., 149., 149., 149., 155., 164., 170., 164., 155., 149.,
     &  149., 149., 149.,
     &  195., 195., 194., 193., 194., 195., 194., 195., 194., 193.,
     &  194., 195., 195./
       DATA (CXN8NE(I),I=1,13) /   57.,  -17.,  -90.,  -19.,  260.,
     &   431.,  417.,  431.,  260.,  -19.,  -90.,  -17.,   57./
       DATA (CDN8NE(I),I=1,13) /   18.,   -7.,  -32.,   31.,  139.,
     &   320.,  351.,  320.,  139.,   31.,  -32.,   -7.,   18./
       DATA ((P5NE(I,J),I=1,13),J=1,3) /
     &  112., 105.,  98., 100.,  98.,  91.,  90.,  91.,  98., 100.,
     &   98., 105., 112.,
     &  162., 161., 160., 163., 166., 165., 165., 165., 166., 163.,
     &  160., 161., 162.,
     &  208., 205., 203., 201., 202., 207., 208., 207., 202., 201.,
     &  203., 205., 208./
       DATA (CXN5NE(I),I=1,13) /  575.,  437.,  299.,  198.,  293.,
     &   384.,  403.,  384.,  293.,  198.,  299.,  437.,  575./
       DATA (CDN5NE(I),I=1,13) /  353.,  256.,  159.,   58.,  203.,
     &   253.,  230.,  253.,  203.,   58.,  159.,  256.,  353./
       DATA ((P3NE(I,J),I=1,13),J=1,3) /
     &   72.,  78.,  85.,  90.,  94.,  89.,  88.,  89.,  94.,  90.,
     &   85.,  78.,  72.,
     &  163., 159., 154., 160., 164., 161., 157., 161., 164., 160.,
     &  154., 159., 163.,
     &  212., 211., 210., 205., 208., 224., 231., 224., 208., 205.,
     &  210., 211., 212./
       DATA (CXN3NE(I),I=1,13) / 1457.,  763.,   68.,  153.,  187.,
     &   390.,  364.,  390.,  187.,  153.,   68.,  763., 1457./
       DATA (CDN3NE(I),I=1,13) /  722.,  329.,  -63.,   40.,   31.,
     &   102.,  111.,  102.,   31.,   40.,  -63.,  329.,  722./
C      Solstice
       DATA ((P2DS(I,J),I=1,13),J=1,3) /
     &   90.,  90.,  91.,  91.,  98., 105., 112., 119., 125., 124.,
     &  120., 126., 131.,
     &  155., 156., 156., 161., 163., 158., 152., 155., 154., 147.,
     &  152., 148., 145.,
     &  205., 202., 199., 191., 195., 206., 212., 210., 206., 205.,
     &  203., 201., 199./
       DATA (CXN2DS(I),I=1,13) / -714., -294.,  126.,  546.,  516.,
     &   486.,  456.,  426.,  396.,  225.,  366.,  710., 1055./
       DATA (CDN2DS(I),I=1,13) / -732., -357.,   17.,  392.,  327.,
     &   262.,  197.,  132.,   67.,  -80.,   62.,  330.,  598./
       DATA ((P1DS(I,J),I=1,13),J=1,3) /
     &  105., 107., 108., 112., 111., 111., 112., 113., 113., 113.,
     &  114., 116., 118.,
     &  125., 127., 130., 132., 135., 134., 135., 139., 134., 133.,
     &  135., 136., 137.,
     &  189., 193., 198., 211., 218., 202., 199., 200., 206., 213.,
     &  217., 224., 231./
       DATA (CXN1DS(I),I=1,13) /  706.,  480.,  255.,   56.,  310.,
     &   111.,  164.,  249.,  324.,  629.,  837.,  938., 1038./
       DATA (CDN1DS(I),I=1,13) /   25.,   51.,   77.,   26.,  107.,
     &   131.,  118.,  141.,  130.,   68.,   45.,  123.,  202./
       DATA ((P8DS(I,J),I=1,13),J=1,3) /
     &  125., 122., 120., 121., 119., 116., 113., 110., 107., 114.,
     &  125., 127., 129.,
     &  148., 147., 146., 149., 159., 163., 162., 160., 171., 160.,
     &  154., 152., 150.,
     &  199., 199., 198., 200., 198., 200., 214., 221., 215., 192.,
     &  193., 196., 200./
       DATA (CXN8DS(I),I=1,13) /  201.,   47., -106., -130., -106.,
     &    26.,  158.,  290.,  422.,  214.,  185.,  224.,  263./
       DATA (CDN8DS(I),I=1,13) /  104.,   63.,   22.,   -3.,   -1.,
     &    60.,  121.,  182.,  242.,  170.,   51.,   67.,   84./
       DATA ((P5DS(I,J),I=1,13),J=1,3) /
     &   91.,  88.,  86.,  95.,  94.,  88.,  85.,  90.,  95.,  96.,
     &   96.,  97.,  98.,
     &  165., 167., 168., 150., 159., 161., 161., 159., 158., 158.,
     &  145., 149., 154.,
     &  189., 191., 194., 201., 204., 200., 196., 193., 191., 196.,
     &  204., 201., 199./
       DATA (CXN5DS(I),I=1,13) /-1839.,-1044., -248., -253., -416.,
     &    17.,  156.,  118., -110.,  372.,  934.,  895.,  856./
       DATA (CDN5DS(I),I=1,13) /-1603.,-1092., -581., -656., -433.,
     &  -166.,   51.,  -29., -325.,  117.,  -36.,  321.,  678./
       DATA ((P3DS(I,J),I=1,13),J=1,3) /
     &   97.,  94.,  92.,  94.,  90.,  89.,  90.,  92.,  93.,  91.,
     &   95.,  85.,  75.,
     &  155., 148., 141., 142., 142., 140., 139., 142., 145., 141.,
     &  144., 155., 167.,
     &  197., 200., 204., 207., 200., 203., 206., 209., 212., 219.,
     &  219., 204., 190./
       DATA (CXN3DS(I),I=1,13) / -894., -797., -700., -604.,  139.,
     &   246.,  354.,  462.,  569.,  142.,   81.,  293.,  505./
       DATA (CDN3DS(I),I=1,13) / -740., -577., -414., -140.,   55.,
     &   -27.,  110.,   74.,  -10., -122.,   52.,  209.,  366./
       DATA ((P2NS(I,J),I=1,13),J=1,3) /
     &  105., 107., 110., 112., 115., 110., 104.,  99.,  95., 103.,
     &  117., 127., 136.,
     &  149., 146., 144., 147., 143., 139., 138., 143., 152., 154.,
     &  154., 148., 142.,
     &  198., 200., 201., 206., 211., 210., 209., 200., 196., 197.,
     &  205., 205., 206./
       DATA (CXN2NS(I),I=1,13) / -817., -511., -205.,  100.,  406.,
     &   312.,  217.,  123.,  640., 1037.,  580.,  677.,  773./
       DATA (CDN2NS(I),I=1,13) /  755.,  684.,  612.,  540.,  469.,
     &   191.,  -87., -365., -171., -224., -112.,  146.,  403./
       DATA ((P1NS(I,J),I=1,13),J=1,3) /
     &  104., 108., 112., 112., 114., 113., 112., 112., 112., 114.,
     &  114., 116., 119.,
     &  137., 135., 133., 131., 136., 137., 133., 134., 131., 133.,
     &  135., 135., 134.,
     &  240., 224., 208., 207., 204., 201., 199., 208., 208., 211.,
     &  196., 198., 200./
       DATA (CXN1NS(I),I=1,13) / 1440.,  430., -580.,  201.,  354.,
     &   470.,  281.,  329.,  511.,  546.,  741.,  911., 1082./
       DATA (CDN1NS(I),I=1,13) / -236., -256., -276.,  182.,  155.,
     &   100.,  128.,  148.,  100.,  -61., -102.,   33.,  169./
       DATA ((P8NS(I,J),I=1,13),J=1,3) /
     &  128., 127., 127., 127., 126., 122., 102., 115., 116., 118.,
     &  121., 121., 121.,
     &  153., 153., 153., 153., 154., 158., 163., 156., 150., 146.,
     &  147., 147., 147.,
     &  198., 197., 197., 195., 193., 194., 197., 198., 198., 198.,
     &  198., 198., 199./
       DATA (CXN8NS(I),I=1,13) /  199.,   -8., -215.,   92.,  275.,
     &   284.,  304.,  322.,  284.,  161.,  121.,  200.,  280./
       DATA (CDN8NS(I),I=1,13) /  123.,    1., -121.,   -6.,  151.,
     &   146.,  -19.,  169.,  149.,   75.,   22.,   78.,  134./
       DATA ((P5NS(I,J),I=1,13),J=1,3) /
     &   85.,  89.,  93., 102.,  95.,  85.,  83.,  86.,  97., 101.,
     &   93.,  92.,  90.,
     &  166., 163., 160., 167., 168., 167., 165., 167., 169., 164.,
     &  157., 162., 166.,
     &  217., 207., 196., 183., 193., 194., 197., 205., 209., 205.,
     &  205., 201., 198./
       DATA (CXN5NS(I),I=1,13) /  274.,  246.,  218.,   90.,  248.,
     &   337.,  370.,  433.,  361.,  461.,  541.,  615.,  689./
       DATA (CDN5NS(I),I=1,13) / -153.,  -61.,   30.,   48.,  169.,
     &   241.,  279.,  260.,  187.,  380.,  316.,  436.,  557./
       DATA ((P3NS(I,J),I=1,13),J=1,3) /
     &   74.,  80.,  87.,  89.,  95.,  94.,  92.,  91.,  90.,  93.,
     &   93.,  87.,  81.,
     &  151., 146., 141., 142., 151., 158., 159., 159., 160., 155.,
     &  149., 153., 157.,
     &  221., 205., 189., 184., 180., 175., 176., 179., 191., 199.,
     &  210., 204., 198./
       DATA (CXN3NS(I),I=1,13) /   68.,   34.,    1.,  154.,  151.,
     &   178.,  224.,  198.,  220.,  380.,  559.,  228., -104./
       DATA (CDN3NS(I),I=1,13) /   97.,    1.,  -95.,   47.,   37.,
     &   112.,   96.,  106.,  129.,  252.,  470.,  220.,  -30./
C
       DO 5 I=1,13
        TXN2DE(I)=CXN2DE(I)
        TDN2DE(I)=CDN2DE(I)
        TXN1DE(I)=CXN1DE(I)
        TDN1DE(I)=CDN1DE(I)
        TXN8DE(I)=CXN8DE(I)
        TDN8DE(I)=CDN8DE(I)
        TXN5DE(I)=CXN5DE(I)
        TDN5DE(I)=CDN5DE(I)
        TXN3DE(I)=CXN3DE(I)
        TDN3DE(I)=CDN3DE(I)
        TXN2NE(I)=CXN2NE(I)
        TDN2NE(I)=CDN2NE(I)
        TXN1NE(I)=CXN1NE(I)
        TDN1NE(I)=CDN1NE(I)
        TXN8NE(I)=CXN8NE(I)
        TDN8NE(I)=CDN8NE(I)
        TXN5NE(I)=CXN5NE(I)
        TDN5NE(I)=CDN5NE(I)
        TXN3NE(I)=CXN3NE(I)
        TDN3NE(I)=CDN3NE(I)       
        TXN2DS(I)=CXN2DS(I)
        TDN2DS(I)=CDN2DS(I)
        TXN1DS(I)=CXN1DS(I)
        TDN1DS(I)=CDN1DS(I)
        TXN8DS(I)=CXN8DS(I)
        TDN8DS(I)=CDN8DS(I)
        TXN5DS(I)=CXN5DS(I)
        TDN5DS(I)=CDN5DS(I)
        TXN3DS(I)=CXN3DS(I)
        TDN3DS(I)=CDN3DS(I)
        TXN2NS(I)=CXN2NS(I)
        TDN2NS(I)=CDN2NS(I)
        TXN1NS(I)=CXN1NS(I)
        TDN1NS(I)=CDN1NS(I)
        TXN8NS(I)=CXN8NS(I)
        TDN8NS(I)=CDN8NS(I)
        TXN5NS(I)=CXN5NS(I)
        TDN5NS(I)=CDN5NS(I)
        TXN3NS(I)=CXN3NS(I)
5       TDN3NS(I)=CDN3NS(I)       
C
       IF (((DDD .GE. 265) .AND. (DDD .LT. 354)) .OR. 
     &    ((DDD .GE. 354) .OR. (DDD .LT. 79))) THEN
         CALL SWAPEL(13,TXN2DS)
         CALL SWAPEL(13,TDN2DS)
         CALL SWAPEL(13,TXN1DS)
         CALL SWAPEL(13,TDN1DS)
         CALL SWAPEL(13,TXN8DS)
         CALL SWAPEL(13,TDN8DS)
         CALL SWAPEL(13,TXN5DS)
         CALL SWAPEL(13,TDN5DS)
         CALL SWAPEL(13,TXN3DS)
         CALL SWAPEL(13,TDN3DS)
         CALL SWAPEL(13,TXN2NS)
         CALL SWAPEL(13,TDN2NS)
         CALL SWAPEL(13,TXN1NS)
         CALL SWAPEL(13,TDN1NS)
         CALL SWAPEL(13,TXN8NS)
         CALL SWAPEL(13,TDN8NS)
         CALL SWAPEL(13,TXN5NS)
         CALL SWAPEL(13,TDN5NS)
         CALL SWAPEL(13,TXN3NS)
         CALL SWAPEL(13,TDN3NS)
       END IF
       
C      interpolated Te values for invdip
C      Te max-min day equinox
       XN2DEI=INTERP(13,0,TXN2DE,INVDPQ,INVDIP)
       XN1DEI=INTERP(13,0,TXN1DE,INVDPQ,INVDIP)
       XN8DEI=INTERP(13,0,TXN8DE,INVDPQ,INVDIP)
       XN5DEI=INTERP(13,0,TXN5DE,INVDPQ,INVDIP)
       XN3DEI=INTERP(13,0,TXN3DE,INVDPQ,INVDIP)
C      Te max-min night equinox       
       XN2NEI=INTERP(13,0,TXN2NE,INVDPQ,INVDIP)
       XN1NEI=INTERP(13,0,TXN1NE,INVDPQ,INVDIP)
       XN8NEI=INTERP(13,0,TXN8NE,INVDPQ,INVDIP)
       XN5NEI=INTERP(13,0,TXN5NE,INVDPQ,INVDIP)
       XN3NEI=INTERP(13,0,TXN3NE,INVDPQ,INVDIP)
C      Te med-min day equinox       
       DN2DEI=INTERP(13,0,TDN2DE,INVDPQ,INVDIP)
       DN1DEI=INTERP(13,0,TDN1DE,INVDPQ,INVDIP)
       DN8DEI=INTERP(13,0,TDN8DE,INVDPQ,INVDIP)
       DN5DEI=INTERP(13,0,TDN5DE,INVDPQ,INVDIP)
       DN3DEI=INTERP(13,0,TDN3DE,INVDPQ,INVDIP)
C      Te med-min night equinox       
       DN2NEI=INTERP(13,0,TDN2NE,INVDPQ,INVDIP)
       DN1NEI=INTERP(13,0,TDN1NE,INVDPQ,INVDIP)
       DN8NEI=INTERP(13,0,TDN8NE,INVDPQ,INVDIP)
       DN5NEI=INTERP(13,0,TDN5NE,INVDPQ,INVDIP)
       DN3NEI=INTERP(13,0,TDN3NE,INVDPQ,INVDIP)
C       
C      Te max-min day solstice       
       XN2DSI=INTERP(13,0,TXN2DS,INVDPQ,INVDIP)
       XN1DSI=INTERP(13,0,TXN1DS,INVDPQ,INVDIP)
       XN8DSI=INTERP(13,0,TXN8DS,INVDPQ,INVDIP)
       XN5DSI=INTERP(13,0,TXN5DS,INVDPQ,INVDIP)
       XN3DSI=INTERP(13,0,TXN3DS,INVDPQ,INVDIP)
C      Te max-min night solstice            
       XN2NSI=INTERP(13,0,TXN2NS,INVDPQ,INVDIP)
       XN1NSI=INTERP(13,0,TXN1NS,INVDPQ,INVDIP)
       XN8NSI=INTERP(13,0,TXN8NS,INVDPQ,INVDIP)
       XN5NSI=INTERP(13,0,TXN5NS,INVDPQ,INVDIP)
       XN3NSI=INTERP(13,0,TXN3NS,INVDPQ,INVDIP)
C      Te med-min day solstice      
       DN2DSI=INTERP(13,0,TDN2DS,INVDPQ,INVDIP)
       DN1DSI=INTERP(13,0,TDN1DS,INVDPQ,INVDIP)
       DN8DSI=INTERP(13,0,TDN8DS,INVDPQ,INVDIP)
       DN5DSI=INTERP(13,0,TDN5DS,INVDPQ,INVDIP)
       DN3DSI=INTERP(13,0,TDN3DS,INVDPQ,INVDIP)
C      Te med-min night solstice           
       DN2NSI=INTERP(13,0,TDN2NS,INVDPQ,INVDIP)
       DN1NSI=INTERP(13,0,TDN1NS,INVDPQ,INVDIP)
       DN8NSI=INTERP(13,0,TDN8NS,INVDPQ,INVDIP)
       DN5NSI=INTERP(13,0,TDN5NS,INVDPQ,INVDIP)
       DN3NSI=INTERP(13,0,TDN3NS,INVDPQ,INVDIP)
       DO 10 I=1,3 
        P2DEI(I)=INTERP(13,1,P2DE(1,I),INVDPQ,INVDIP)
        P2NEI(I)=INTERP(13,1,P2NE(1,I),INVDPQ,INVDIP)
        P1DEI(I)=INTERP(13,1,P1DE(1,I),INVDPQ,INVDIP)
        P1NEI(I)=INTERP(13,1,P1NE(1,I),INVDPQ,INVDIP)
        P8DEI(I)=INTERP(13,1,P8DE(1,I),INVDPQ,INVDIP)
        P8NEI(I)=INTERP(13,1,P8NE(1,I),INVDPQ,INVDIP)
        P5DEI(I)=INTERP(13,1,P5DE(1,I),INVDPQ,INVDIP)
        P5NEI(I)=INTERP(13,1,P5NE(1,I),INVDPQ,INVDIP)
        P3DEI(I)=INTERP(13,1,P3DE(1,I),INVDPQ,INVDIP)
        P3NEI(I)=INTERP(13,1,P3NE(1,I),INVDPQ,INVDIP)
                                            
        P2DSI(I)=INTERP(13,1,P2DS(1,I),INVDPQ,INVDIP)
        P2NSI(I)=INTERP(13,1,P2NS(1,I),INVDPQ,INVDIP)
        P1DSI(I)=INTERP(13,1,P1DS(1,I),INVDPQ,INVDIP)
        P1NSI(I)=INTERP(13,1,P1NS(1,I),INVDPQ,INVDIP)
        P8DSI(I)=INTERP(13,1,P8DS(1,I),INVDPQ,INVDIP)
        P8NSI(I)=INTERP(13,1,P8NS(1,I),INVDPQ,INVDIP)
        P5DSI(I)=INTERP(13,1,P5DS(1,I),INVDPQ,INVDIP)
        P5NSI(I)=INTERP(13,1,P5NS(1,I),INVDPQ,INVDIP)
        P3DSI(I)=INTERP(13,1,P3DS(1,I),INVDPQ,INVDIP)
10      P3NSI(I)=INTERP(13,1,P3NS(1,I),INVDPQ,INVDIP)
       MLTTMP=MLT-1
       IF (MLTTMP .LT. 0) MLTTMP=MLTTMP+24.0
       MLTRAD=MLTTMP/24.0*2*3.1415927
      IF (((DDD .GE. 79) .AND. (DDD .LT. 171)) .OR. 
     &    ((DDD .GE. 265) .AND. (DDD .LT. 354))) THEN
C      Equinox     
       CALL TPCAS(MLTRAD,PF107,P2000A,
     &             XN2DEI,DN2DEI,P2DEI,XN2NEI,DN2NEI,P2NEI,TP200A)
       CALL TPCAS(MLTRAD,PF107,P1400A,
     &             XN1DEI,DN1DEI,P1DEI,XN1NEI,DN1NEI,P1NEI,TP140A)
       CALL TPCAS(MLTRAD,PF107,P850A,
     &             XN8DEI,DN8DEI,P8DEI,XN8NEI,DN8NEI,P8NEI,TP850A)
       CALL TPCAS(MLTRAD,PF107,P550A,
     &             XN5DEI,DN5DEI,P5DEI,XN5NEI,DN5NEI,P5NEI,TP550A)
       CALL TPCAS(MLTRAD,PF107,P350A,
     &             XN3DEI,DN3DEI,P3DEI,XN3NEI,DN3NEI,P3NEI,TP350A)
C       Solstice     
       CALL TPCAS(MLTRAD,PF107,P2000B,
     &             XN2DSI,DN2DSI,P2DSI,XN2NSI,DN2NSI,P2NSI,TP200B)
       CALL TPCAS(MLTRAD,PF107,P1400B,
     &             XN1DSI,DN1DSI,P1DSI,XN1NSI,DN1NSI,P1NSI,TP140B)
       CALL TPCAS(MLTRAD,PF107,P850B,
     &             XN8DSI,DN8DSI,P8DSI,XN8NSI,DN8NSI,P8NSI,TP850B)
       CALL TPCAS(MLTRAD,PF107,P550B,
     &             XN5DSI,DN5DSI,P5DSI,XN5NSI,DN5NSI,P5NSI,TP550B)
       CALL TPCAS(MLTRAD,PF107,P350B,
     &             XN3DSI,DN3DSI,P3DSI,XN3NSI,DN3NSI,P3NSI,TP350B)      
      END IF
      IF (((DDD .GE. 171) .AND. (DDD .LT. 265)) .OR. 
     &    ((DDD .GE. 354) .OR. (DDD .LT. 79)))  THEN
C       Solstice     
       CALL TPCAS(MLTRAD,PF107,P2000A,
     &             XN2DSI,DN2DSI,P2DSI,XN2NSI,DN2NSI,P2NSI,TP200A)
       CALL TPCAS(MLTRAD,PF107,P1400A,
     &             XN1DSI,DN1DSI,P1DSI,XN1NSI,DN1NSI,P1NSI,TP140A)
       CALL TPCAS(MLTRAD,PF107,P850A,
     &             XN8DSI,DN8DSI,P8DSI,XN8NSI,DN8NSI,P8NSI,TP850A)
       CALL TPCAS(MLTRAD,PF107,P550A,
     &             XN5DSI,DN5DSI,P5DSI,XN5NSI,DN5NSI,P5NSI,TP550A)
       CALL TPCAS(MLTRAD,PF107,P350A,
     &             XN3DSI,DN3DSI,P3DSI,XN3NSI,DN3NSI,P3NSI,TP350A)  
C      Equinox     
       CALL TPCAS(MLTRAD,PF107,P2000B,
     &             XN2DEI,DN2DEI,P2DEI,XN2NEI,DN2NEI,P2NEI,TP200B)
       CALL TPCAS(MLTRAD,PF107,P1400B,
     &             XN1DEI,DN1DEI,P1DEI,XN1NEI,DN1NEI,P1NEI,TP140B)
       CALL TPCAS(MLTRAD,PF107,P850B,
     &             XN8DEI,DN8DEI,P8DEI,XN8NEI,DN8NEI,P8NEI,TP850B)
       CALL TPCAS(MLTRAD,PF107,P550B,
     &             XN5DEI,DN5DEI,P5DEI,XN5NEI,DN5NEI,P5NEI,TP550B)
       CALL TPCAS(MLTRAD,PF107,P350B,
     &             XN3DEI,DN3DEI,P3DEI,XN3NEI,DN3NEI,P3NEI,TP350B)     
      END IF 
      RETURN
      END
c
c
      SUBROUTINE TEBA(DIPL,SLT,NS,TE) 
C CALCULATES ELECTRON TEMPERATURES TE(1) TO TE(4) AT ALTITUDES
C 300, 400, 1400 AND 3000 KM FOR DIP-LATITUDE DIPL/DEG AND 
C LOCAL SOLAR TIME SLT/H USING THE BRACE-THEIS-MODELS (J. ATMOS.
C TERR. PHYS. 43, 1317, 1981); NS IS SEASON IN NORTHERN
C HEMISOHERE: IS=1 SPRING, IS=2 SUMMER ....
C ALSO CALCULATED ARE THE TEMPERATURES AT 400 KM ALTITUDE FOR
C MIDNIGHT (TE(5)) AND NOON (TE(6)).   
      DIMENSION C(4,2,81),A(82),TE(6)
      COMMON    /CONST/UMR,PI      /const1/humr,dumr
      DATA (C(1,1,J),J=1,81)/                      
     &.3100E1,-.3215E-2,.2440E+0,-.4613E-3,-.1711E-1,.2605E-1,                  
     &-.9546E-1,.1794E-1,.1270E-1,.2791E-1,.1536E-1,-.6629E-2,                  
     &-.3616E-2,.1229E-1,.4147E-3,.1447E-2,-.4453E-3,-.1853,                    
     &-.1245E-1,-.3675E-1,.4965E-2,.5460E-2,.8117E-2,-.1002E-1,                 
     &.5466E-3,-.3087E-1,-.3435E-2,-.1107E-3,.2199E-2,.4115E-3,                 
     &.6061E-3,.2916E-3,-.6584E-1,.4729E-2,-.1523E-2,.6689E-3,                  
     &.1031E-2,.5398E-3,-.1924E-2,-.4565E-1,.7244E-2,-.8543E-4,                 
     &.1052E-2,-.6696E-3,-.7492E-3,.4405E-1,.3047E-2,.2858E-2,                  
     &-.1465E-3,.1195E-2,-.1024E-3,.4582E-1,.8749E-3,.3011E-3,                  
     &.4473E-3,-.2782E-3,.4911E-1,-.1016E-1,.27E-2,-.9304E-3,                   
     &-.1202E-2,.2210E-1,.2566E-2,-.122E-3,.3987E-3,-.5744E-1,                  
     &.4408E-2,-.3497E-2,.83E-3,-.3536E-1,-.8813E-2,.2423E-2,                   
     &-.2994E-1,-.1929E-2,-.5268E-3,-.2228E-1,.3385E-2,                         
     &.413E-1,.4876E-2,.2692E-1,.1684E-2/          
      DATA (C(1,2,J),J=1,81)/.313654E1,.6796E-2,.181413,.8564E-1,               
     &-.32856E-1,-.3508E-2,-.1438E-1,-.2454E-1,.2745E-2,.5284E-1,               
     &.1136E-1,-.1956E-1,-.5805E-2,.2801E-2,-.1211E-2,.4127E-2,                 
     &.2909E-2,-.25751,-.37915E-2,-.136E-1,-.13225E-1,.1202E-1,                 
     &.1256E-1,-.12165E-1,.1326E-1,-.7123E-1,.5793E-3,.1537E-2,                 
     &.6914E-2,-.4173E-2,.1052E-3,-.5765E-3,-.4041E-1,-.1752E-2,                
     &-.542E-2,-.684E-2,.8921E-3,-.2228E-2,.1428E-2,.6635E-2,-.48045E-2,        
     &-.1659E-2,-.9341E-3,.223E-3,-.9995E-3,.4285E-1,-.5211E-3,                 
     &-.3293E-2,.179E-2,.6435E-3,-.1891E-3,.3844E-1,.359E-2,-.8139E-3,          
     &-.1996E-2,.2398E-3,.2938E-1,.761E-2,.347655E-2,.1707E-2,.2769E-3,         
     &-.157E-1,.983E-3,-.6532E-3,.929E-4,-.2506E-1,.4681E-2,.1461E-2,           
     &-.3757E-5,-.9728E-2,.2315E-2,.6377E-3,-.1705E-1,.2767E-2,                 
     &-.6992E-3,-.115E-1,-.1644E-2,.3355E-2,-.4326E-2,.2035E-1,.2985E-1/        
      DATA (C(2,1,J),J=1,81)/.3136E1,.6498E-2,.2289,.1859E-1,-.3328E-1,         
     &-.4889E-2,-.3054E-1,-.1773E-1,-.1728E-1,.6555E-1,.1775E-1,                
     &-.2488E-1,-.9498E-2,.1493E-1,.281E-2,.2406E-2,.5436E-2,-.2115,            
     &.7007E-2,-.5129E-1,-.7327E-2,.2402E-1,.4772E-2,-.7374E-2,                 
     &-.3835E-3,-.5013E-1,.2866E-2,.2216E-2,.2412E-3,.2094E-2,.122E-2           
     &,-.1703E-3,-.1082,-.4992E-2,-.4065E-2,.3615E-2,-.2738E-2,                 
     &-.7177E-3,.2173E-3,-.4373E-1,-.375E-2,.5507E-2,-.1567E-2,                 
     &-.1458E-2,-.7397E-3,.7903E-1,.4131E-2,.3714E-2,.1073E-2,                  
     &-.8991E-3,.2976E-3,.2623E-1,.2344E-2,.5608E-3,.4124E-3,.1509E-3,          
     &.5103E-1,.345E-2,.1283E-2,.7238E-3,-.3464E-4,.1663E-1,-.1644E-2,          
     &-.71E-3,.5281E-3,-.2729E-1,.3556E-2,-.3391E-2,-.1787E-3,.2154E-2,         
     &.6476E-2,-.8282E-3,-.2361E-1,.9557E-3,.3205E-3,-.2301E-1,                 
     &-.854E-3,-.1126E-1,-.2323E-2,-.8582E-2,.2683E-1/                          
      DATA (C(2,2,J),J=1,81)/.3144E1,.8571E-2,.2539,.6937E-1,-.1667E-1,         
     &.2249E-1,-.4162E-1,.1201E-1,.2435E-1,.5232E-1,.2521E-1,-.199E-1,          
     &-.7671E-2,.1264E-1,-.1551E-2,-.1928E-2,.3652E-2,-.2019,.5697E-2,          
     &-.3159E-1,-.1451E-1,.2868E-1,.1377E-1,-.4383E-2,.1172E-1,                 
     &-.5683E-1,.3593E-2,.3571E-2,.3282E-2,.1732E-2,-.4921E-3,-.1165E-2         
     &,-.1066,-.1892E-1,.357E-2,-.8631E-3,-.1876E-2,-.8414E-4,.2356E-2,         
     &-.4259E-1,-.322E-2,.4641E-2,.6223E-3,-.168E-2,-.1243E-3,.7393E-1,         
     &-.3143E-2,-.2362E-2,.1235E-2,-.1551E-2,.2099E-3,.2299E-1,.5301E-2         
     &,-.4306E-2,-.1303E-2,.7687E-5,.5305E-1,.6642E-2,-.1686E-2,                
     &.1048E-2,.5958E-3,.4341E-1,-.8819E-4,-.333E-3,-.2158E-3,-.4106E-1         
     &,.4191E-2,.2045E-2,-.1437E-3,-.1803E-1,-.8072E-3,-.424E-3,                
     &-.26E-1,-.2329E-2,.5949E-3,-.1371E-1,-.2188E-2,.1788E-1,                  
     &.6405E-3,.5977E-2,.1333E-1/                  
      DATA (C(3,1,J),J=1,81)/.3372E1,.1006E-1,.1436,.2023E-2,-.5166E-1,         
     &.9606E-2,-.5596E-1,.4914E-3,-.3124E-2,-.4713E-1,-.7371E-2,                
     &-.4823E-2,-.2213E-2,.6569E-2,-.1962E-3,.3309E-3,-.3908E-3,                
     &-.2836,.7829E-2,.1175E-1,.9919E-3,.6589E-2,.2045E-2,-.7346E-2             
     &,-.89E-3,-.347E-1,-.4977E-2,.147E-2,-.2823E-5,.6465E-3,                   
     &-.1448E-3,.1401E-2,-.8988E-1,-.3293E-4,-.1848E-2,.4439E-3,                
     &-.1263E-2,.317E-3,-.6227E-3,.1721E-1,-.199E-2,-.4627E-3,                  
     &.2897E-5,-.5454E-3,.3385E-3,.8432E-1,-.1951E-2,.1487E-2,                  
     &.1042E-2,-.4788E-3,-.1276E-3,.2373E-1,.2409E-2,.5263E-3,                  
     &.1301E-2,-.4177E-3,.3974E-1,.1418E-3,-.1048E-2,-.2982E-3,                 
     &-.3396E-4,.131E-1,.1413E-2,-.1373E-3,.2638E-3,-.4171E-1,                  
     &-.5932E-3,-.7523E-3,-.6883E-3,-.2355E-1,.5695E-3,-.2219E-4,               
     &-.2301E-1,-.9962E-4,-.6761E-3,.204E-2,-.5479E-3,.2591E-1,                 
     &-.2425E-2,.1583E-1,.9577E-2/                 
      DATA (C(3,2,J),J=1,81)/.3367E1,.1038E-1,.1407,.3622E-1,-.3144E-1,         
     &.112E-1,-.5674E-1,.3219E-1,.1288E-2,-.5799E-1,-.4609E-2,                  
     &.3252E-2,-.2859E-3,.1226E-1,-.4539E-2,.1310E-2,-.5603E-3,                 
     &-.311,-.1268E-2,.1539E-1,.3146E-2,.7787E-2,-.143E-2,-.482E-2              
     &,.2924E-2,-.9981E-1,-.7838E-2,-.1663E-3,.4769E-3,.4148E-2,                
     &-.1008E-2,-.979E-3,-.9049E-1,-.2994E-2,-.6748E-2,-.9889E-3,               
     &.1488E-2,-.1154E-2,-.8412E-4,-.1302E-1,-.4859E-2,-.7172E-3,               
     &-.9401E-3,.9101E-3,-.1735E-3,.7055E-1,.6398E-2,-.3103E-2,                 
     &-.938E-3,-.4E-3,-.1165E-2,.2713E-1,-.1654E-2,.2781E-2,                    
     &-.5215E-5,.2258E-3,.5022E-1,.95E-2,.4147E-3,.3499E-3,                     
     &-.6097E-3,.4118E-1,.6556E-2,.3793E-2,-.1226E-3,-.2517E-1,                 
     &.1491E-3,.1075E-2,.4531E-3,-.9012E-2,.3343E-2,.3431E-2,                   
     &-.2519E-1,.3793E-4,.5973E-3,-.1423E-1,-.132E-2,-.6048E-2,                 
     &-.5005E-2,-.115E-1,.2574E-1/                 
      DATA (C(4,1,J),J=1,81)/.3574E1,.0,.7537E-1,.0,-.8459E-1,                  
     &0.,-.294E-1,0.,.4547E-1,-.5321E-1,0.,.4328E-2,0.,.6022E-2,                
     &.0,-.9168E-3,.0,-.1768,.0,.294E-1,.0,.5902E-3,.0,-.9047E-2,               
     &.0,-.6555E-1,.0,-.1033E-2,.0,.1674E-2,.0,.2802E-3,-.6786E-1               
     &,.0,.4193E-2,.0,-.6448E-3,.0,.9277E-3,-.1634E-1,.0,-.2531E-2              
     &,.0,.193E-4,.0,.528E-1,.0,.2438E-2,.0,-.5292E-3,.0,.1555E-1               
     &,.0,-.3259E-2,.0,-.5998E-3,.3168E-1,.0,.2382E-2,.0,-.4078E-3              
     &,.2312E-1,.0,.1481E-3,.0,-.1885E-1,.0,.1144E-2,.0,-.9952E-2               
     &,.0,-.551E-3,-.202E-1,.0,-.7283E-4,-.1272E-1,.0,.2224E-2,                 
     &.0,-.251E-2,.2434E-1/                        
      DATA (C(4,2,J),J=1,81)/.3574E1,-.5639E-2,.7094E-1,                        
     &-.3347E-1,-.861E-1,-.2877E-1,-.3154E-1,-.2847E-2,.1235E-1,                
     &-.5966E-1,-.3236E-2,.3795E-3,-.8634E-3,.3377E-2,-.1071E-3,                
     &-.2151E-2,-.4057E-3,-.1783,.126E-1,.2835E-1,-.242E-2,                     
     &.3002E-2,-.4684E-2,-.6756E-2,-.7493E-3,-.6147E-1,-.5636E-2                
     &,-.1234E-2,-.1613E-2,-.6353E-4,-.2503E-3,-.1729E-3,-.7148E-1              
     &,.5326E-2,.4006E-2,.6484E-3,-.1046E-3,-.6034E-3,-.9435E-3,                
     &-.2385E-2,.6853E-2,.151E-2,.1319E-2,.9049E-4,-.1999E-3,                   
     &.3976E-1,.2802E-2,-.103E-2,.5599E-3,-.4791E-3,-.846E-4,                   
     &.2683E-1,.427E-2,.5911E-3,.2987E-3,-.208E-3,.1396E-1,                     
     &-.1922E-2,-.1063E-2,.3803E-3,.1343E-3,.1771E-1,-.1038E-2,                 
     &-.4645E-3,-.2481E-3,-.2251E-1,-.29E-2,-.3977E-3,-.516E-3,                 
     &-.8079E-2,-.1528E-2,.306E-3,-.1582E-1,-.8536E-3,.1565E-3,                 
     &-.1252E-1,.2319E-3,.4311E-2,.1024E-2,.1296E-5,.179E-1/ 
                        
        IF(NS.LT.3) THEN
           IS=NS
        ELSE IF(NS.GT.3) THEN
           IS=2
           DIPL=-DIPL
        ELSE
           IS=1
        ENDIF
      COLAT=UMR*(90.-DIPL)                    
      AZ=humr*SLT    
      CALL SPHARM(A,8,8,COLAT,AZ)
        IF(IS.EQ.2) THEN
           KEND=3
        ELSE
           KEND=4
        ENDIF                  
      DO 2 K=1,KEND      
      STE=0.          
      DO 1 I=1,81     
1       STE=STE+A(I)*C(K,IS,I)                       
2     TE(K)=10.**STE
        IF(IS.EQ.2) THEN
           DIPL=-DIPL
           COLAT=UMR*(90.-DIPL)                    
           CALL SPHARM(A,8,8,COLAT,AZ)
           STE=0.          
           DO 11 I=1,81     
11            STE=STE+A(I)*C(4,2,I)                       
           TE(4)=10.**STE
        ENDIF

C---------- TEMPERATURE AT 400KM AT MIDNIGHT AND NOON
      DO 4 J=1,2      
        STE=0.          
        AZ=humr*(J-1)*12.                           
        CALL SPHARM(A,8,8,COLAT,AZ)                  
        DO 3 I=1,81     
3         STE=STE+A(I)*C(2,IS,I)                       
4       TE(J+4)=10.**STE                             
      RETURN          
      END             
C
C
      SUBROUTINE SPHARM(C,L,M,COLAT,AZ)            
C CALCULATES THE COEFFICIENTS OF THE SPHERICAL HARMONIC                         
C EXPANSION THAT WAS USED FOR THE BRACE-THEIS-MODELS.                           
      DIMENSION C(82)                              
      C(1)=1.         
      K=2             
      X=COS(COLAT)    
      C(K)=X          
      K=K+1           
      DO 10 I=2,L     
      C(K)=((2*I-1)*X*C(K-1)-(I-1)*C(K-2))/I       
10    K=K+1           
      Y=SIN(COLAT)    
      DO 20 MT=1,M    
      CAZ=COS(MT*AZ)  
      SAZ=SIN(MT*AZ)  
      C(K)=Y**MT      
      K=K+1           
      IF(MT.EQ.L) GOTO 16                          
      C(K)=C(K-1)*X*(2*MT+1)                       
      K=K+1           
      IF((MT+1).EQ.L) GOTO 16                      
      DO 15 I=2+MT,L  
      C(K)=((2*I-1)*X*C(K-1)-(I+MT-1)*C(K-2))/(I-MT)                            
15    K=K+1           
16    N=L-MT+1        
      DO 18 I=1,N     
      C(K)=C(K-N)*CAZ                              
      C(K-N)=C(K-N)*SAZ                            
18    K=K+1           
20    CONTINUE        
      RETURN          
      END             
C
C
      REAL FUNCTION ELTE(H)
c----------------------------------------------------------------
C ELECTRON TEMPERATURE PROFILE BASED ON THE TEMPERATURES AT 7 FIXED
C HEIGHTS (AH(7)) AND THE TEMPERATURE GRADIENTS BETWEEN THESE THESE 
C HEIGHTS (ST(6)) GIVEN IN THE COMMON BLOCK. ATE1 IS THE TEMPERATURE
C AT THE STARTING HEIGHT 120 KM. D(5) DEFINE THE TRANSITION SPAN FROM
C ONE CONSTANT GRADIENT REGION TO THE NEXT.
c----------------------------------------------------------------
      COMMON /BLOTE/AH(7),ATE1,ST(6),D(5)
C
      SUM=ATE1+ST(1)*(H-AH(1))                     
      DO 1 I=1,5
        aa = eptr(h    ,d(i),ah(i+1))
        bb = eptr(ah(1),d(i),ah(i+1))
1     SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*D(I)                
      ELTE=SUM        
      RETURN          
      END             
C
C
      FUNCTION TEDE(H,DEN,COV)                     
C ELECTRON TEMEPERATURE MODEL AFTER BRACE,THEIS .  
C FOR NEG. COV THE MEAN COV-INDEX (3 SOLAR ROT.) IS EXPECTED.                   
C DEN IS THE ELECTRON DENSITY IN M-3.              
      Y=1051.+(17.01*H-2746.)*                     
     &EXP(-5.122E-4*H+(6.094E-12-3.353E-14*H)*DEN) 
      ACOV=ABS(COV)   
      YC=1.+(.117+2.02E-3*ACOV)/(1.+EXP(-(ACOV-102.5)/5.))                      
      IF(COV.LT.0.)   
     &YC=1.+(.123+1.69E-3*ACOV)/(1.+EXP(-(ACOV-115.)/10.))                      
      TEDE=Y*YC       
      RETURN          
      END             
C
C                     
C*************************************************************                  
C**************** ION TEMPERATURE ****************************
C*************************************************************                  
C
C
      REAL FUNCTION TI(H)
c----------------------------------------------------------------
C ION TEMPERATURE FOR HEIGHTS NOT GREATER 1000 KM AND NOT LESS HS               
C EXPLANATION SEE FUNCTION RPID.                   
c----------------------------------------------------------------
      REAL              MM
      COMMON  /BLOCK8/  HS,TNHS,XSM(4),MM(5),G(4),M

      SUM=MM(1)*(H-HS)+TNHS                        
      DO 100 I=1,M-1  
        aa = eptr(h ,g(i),xsm(i))
        bb = eptr(hs,g(i),xsm(i))
100     SUM=SUM+(MM(I+1)-MM(I))*(AA-BB)*G(I)                
      TI=SUM          
      RETURN          
      END             
C
C                     
C*************************************************************                  
C************* ION RELATIVE PRECENTAGE DENSITY *****************                
C*************************************************************                  
C
C
      REAL FUNCTION RPID (H, H0, N0, M, ST, ID, XS)
c------------------------------------------------------------------
C D.BILITZA,1977,THIS ANALYTIC FUNCTION IS USED TO REPRESENT THE                
C RELATIVE PRECENTAGE DENSITY OF ATOMAR AND MOLECULAR OXYGEN IONS.              
C THE M+1 HEIGHT GRADIENTS ST(M+1) ARE CONNECTED WITH EPSTEIN-                  
C STEP-FUNCTIONS AT THE STEP HEIGHTS XS(M) WITH TRANSITION                      
C THICKNESSES ID(M). RPID(H0,H0,N0,....)=N0.       
C ARGMAX is the highest allowed argument for EXP in your system.
c------------------------------------------------------------------
      REAL              N0         
      DIMENSION         ID(4), ST(5), XS(4)                
      COMMON  /ARGEXP/  ARGMAX

      SUM=(H-H0)*ST(1)                             
      DO 100  I=1,M   
              XI=ID(I)
                aa = eptr(h ,xi,xs(i))
                bb = eptr(h0,xi,xs(i))
100           SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*XI 
      IF(ABS(SUM).LT.ARGMAX) then
        SM=EXP(SUM)
      else IF(SUM.Gt.0.0) then
        SM=EXP(ARGMAX)
      else
        SM=0.0
      endif
      RPID= n0 * SM        
      RETURN          
      END             
C
c
      SUBROUTINE RDHHE (H,HB,RDOH,RDO2H,RNO,PEHE,RDH,RDHE)                      
C BILITZA,FEB.82,H+ AND HE+ RELATIVE PERECENTAGE DENSITY BELOW                  
C 1000 KM. THE O+ AND O2+ REL. PER. DENSITIES SHOULD BE GIVEN                   
C (RDOH,RDO2H). HB IS THE ALTITUDE OF MAXIMAL O+ DENSITY. PEHE                  
C IS THE PRECENTAGE OF HE+ IONS COMPARED TO ALL LIGHT IONS.                     
C RNO IS THE RATIO OF NO+ TO O2+DENSITY AT H=HB.   
      RDHE=0.0        
      RDH=0.0         
      IF(H.LE.HB) GOTO 100                         
      REST=100.0-RDOH-RDO2H-RNO*RDO2H              
      RDH=REST*(1.-PEHE/100.)                      
      RDHE=REST*PEHE/100.                          
100   RETURN          
      END             
C
C
      REAL FUNCTION RDNO(H,HB,RDO2H,RDOH,RNO)      
C D.BILITZA, 1978. NO+ RELATIVE PERCENTAGE DENSITY ABOVE 100KM.                 
C FOR MORE INFORMATION SEE SUBROUTINE RDHHE.       
      IF (H.GT.HB) GOTO 200                        
      RDNO=100.0-RDO2H-RDOH                        
      RETURN          
200   RDNO=RNO*RDO2H  
      RETURN          
      END
C
C
      SUBROUTINE  KOEFP1(PG1O)                     
C THIEMANN,1979,COEFFICIENTS PG1O FOR CALCULATING  O+ PROFILES                  
C BELOW THE F2-MAXIMUM. CHOSEN TO APPROACH DANILOV-                             
C SEMENOV'S COMPILATION.                           
      DIMENSION PG1O(80), FELD (80)  
      DATA FELD/-11.0,-11.0,4.0,-11.0,0.08018,0.13027,0.04216,     
     &0.25,-0.00686,0.00999,5.113,0.1 ,170.0,180.0,0.1175,0.15,    
     &-11.0,1.0 ,2.0,-11.0,0.069,0.161,0.254,0.18,0.0161,0.0216,                             
     &0.03014,0.1,152.0,167.0,0.04916,0.17,-11.0,2.0,2.0,-11.0,                          
     &0.072,0.092,0.014,0.21,0.01389,0.03863,0.05762,0.12,165.0,                           
     &168.0,0.008,0.258,-11.0,1.0,3.0,-11.0,0.091,0.088,0.008,       
     &0.34,0.0067,0.0195,0.04,0.1,158.0,172.0,0.01,0.24,-11.0,                            
     &2.0,3.0,-11.0,0.083,0.102,0.045,0.03,0.00127,0.01,0.05,                            
     &0.09,167.0,185.0,0.015,0.18/     
      DO 10 I=1,80    
10    	PG1O(I)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE KOEFP2(PG2O)                      
C THIEMANN,1979,COEFFICIENTS FOR CALCULATION OF O+ PROFILES                     
C ABOVE THE F2-MAXIMUM (DUMBS,SPENNER:AEROS-COMPILATION)                        
      DIMENSION PG2O(32), FELD(32)   
      DATA FELD/1.0,-11.0,-11.0,1.0,695.0,-.000781,-.00264,                             
     &2177.0,1.0,-11.0,-11.0,2.0,570.0,-.002,-.0052,1040.0,    
     &2.0,-11.0,-11.0,1.0,695.0,-.000786,-.00165,3367.0,2.0, 
     &-11.0,-11.0,2.0,575.0,-.00126,-.00524,1380.0/                
      DO 10 I=1,32    
10    	PG2O(I)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE  KOEFP3(PG3O)                     
C THIEMANN,1979,COEFFICIENTS FOR CALCULATING O2+ PROFILES.                      
C CHOSEN AS TO APPROACH DANILOV-SEMENOV'S COMPILATION.                          
      DIMENSION PG3O(80), FELD(80)   
      DATA FELD/-11.0,1.0,2.0,-11.0,160.0,31.0,130.0,-10.0,                           
     &198.0,0.0,0.05922,-0.07983,-0.00397,0.00085,-0.00313,            
     &0.0,-11.0,2.0,2.0,-11.0,140.0,30.0,130.0,-10.0,190.0,                      
     &0.0,0.05107,-0.07964,0.00097,-0.01118,-0.02614,-0.09537,       
     &-11.0,1.0,3.0,-11.0,140.0,37.0,125.0,0.0,182.0,0.0,                           
     &0.0307,-0.04968,-0.00248,-0.02451,-0.00313,0.0,-11.0,   
     &2.0,3.0,-11.0,140.0,37.0,125.0,0.0,170.0,0.0,0.02806,              
     &-0.04716,0.00066,-0.02763,-0.02247,-0.01919,-11.0,-11.0,                      
     &4.0,-11.0,140.0,45.0,136.0,-9.0,181.0,-26.0,0.02994,                
     &-0.04879,-0.01396,0.00089,-0.09929,0.05589/           
      DO 10 I=1,80    
10    	PG3O(I)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE SUFE (FIELD,RFE,M,FE)             
C SELECTS THE REQUIRED ION DENSITY PARAMETER SET.
C THE INPUT FIELD INCLUDES DIFFERENT SETS OF DIMENSION M EACH                
C CARACTERISED BY 4 HEADER NUMBERS. RFE(4) SHOULD CONTAIN THE                   
C CHOSEN HEADER NUMBERS.FE(M) IS THE CORRESPONDING SET.                         
      DIMENSION RFE(4),FE(12),FIELD(80),EFE(4)     
      K=0             
100   DO 101 I=1,4    
      K=K+1           
101   EFE(I)=FIELD(K)                              
      DO 111 I=1,M    
      K=K+1           
111   FE(I)=FIELD(K)  
      DO 120 I=1,4    
      IF((EFE(I).GT.-10.0).AND.(RFE(I).NE.EFE(I))) GOTO 100                     
120   CONTINUE        
      RETURN          
      END             
C
C
        subroutine iondani(id,ismo,hx,zd,fd,fs,dion)
c-------------------------------------------------------
c       id      day of month
c       ismo    seasonal month (Northern Hemisphere January 
c                   is ismo=1 and so is Southern H. July)
c       hx      altitude in km
c       zd      solar zenith angle in degrees
c       fd      latitude in degrees
c       fs      10.7cm solar radio flux (12-month running mean)
c       dion(1)   O+  relative density in percent
c       dion(2)   H+  relative density in percent
c       dion(3)   N+  relative density in percent
c       dion(4)   He+ relative density in percent
c       dion(5)   NO+ relative density in percent
c       dion(6)   O2+ relative density in percent
c       dion(7)   Cluster+ relative density in percent
c
c Uses ionco2 (DS-95) for the molecular ions and ionco1 (DY-85)
c for the atomic ions.
c-------------------------------------------------------
        dimension       dion(7)
        common  /const/umr,pi

        do 1122 i=1,7
1122    dion(i)=0.

        h = hx
        xhi = zd
        xlati = fd
        f107 = fs
        deci_month = ismo + id/29.0
        if (h.gt.300.) then
        	call ionco1(h,xhi,xlati,f107,deci_month,dion)
			dion(5)=0.0
			dion(6)=0.0
			dion(7)=0.0
        else
        	call ionco2(h,xhi,ismo,f107,rno,ro2,rcl,ro)
        	dion(5)=rno
        	dion(6)=ro2
       		dion(7)=rcl
        	dion(1)=ro
        endif

        return
        end
c
c
        subroutine ionco1(h,zd,fd,fs,t,cn)
c---------------------------------------------------------------
c ion composition model
c   A.D. Danilov and A.P. Yaichnikov, A New Model of the Ion
c   Composition at 75 to 1000 km for IRI, Adv. Space Res. 5, #7,
c   75-79, 107-108, 1985
c
c       h       altitude in km
c       zd      solar zenith angle in degrees
c       fd      latitude in degrees (same result for fd and -fd)
c       fs      10.7cm solar radio flux
c       t       seasonal decimal month (Northern Hemisphere January 
c                   15 is t=1.5 and so is Southern Hemisphere July 15)
c       cn(1)   O+  relative density in percent
c       cn(2)   H+  relative density in percent
c       cn(3)   N+  relative density in percent
c       cn(4)   He+ relative density in percent
c Please note: molecular ions are now computed in IONCO2
c       [cn(5)   NO+ relative density in percent
c       [cn(6)   O2+ relative density in percent
c       [cn(7)   cluster ions  relative density in percent
c---------------------------------------------------------------
c
c        dimension       cn(7),cm(7),hm(7),alh(7),all(7),beth(7),
c     &                  betl(7),p(5,6,7),var(6),po(5,6),ph(5,6),
c     &                  pn(5,6),phe(5,6),pno(5,6),po2(5,6),pcl(5,6)
        dimension       cn(4),cm(4),hm(4),alh(4),all(4),beth(4),
     &                  betl(4),p(5,6,4),var(6),po(5,6),ph(5,6),
     &                  pn(5,6),phe(5,6)

        common  /argexp/argmax
        common  /CONST/UMR,PI
        data po/4*0.,98.5,4*0.,320.,4*0.,-2.59E-4,2.79E-4,-3.33E-3,
     &          -3.52E-3,-5.16E-3,-2.47E-2,4*0.,-2.5E-6,1.04E-3,
     &          -1.79E-4,-4.29E-5,1.01E-5,-1.27E-3/
        data ph/-4.97E-7,-1.21E-1,-1.31E-1,0.,98.1,355.,-191.,
     &          -127.,0.,2040.,4*0.,-4.79E-6,-2.E-4,5.67E-4,
     &          2.6E-4,0.,-5.08E-3,10*0./
        data pn/7.6E-1,-5.62,-4.99,0.,5.79,83.,-369.,-324.,0.,593.,
     &          4*0.,-6.3E-5,-6.74E-3,-7.93E-3,-4.65E-3,0.,-3.26E-3,
     &          4*0.,-1.17E-5,4.88E-3,-1.31E-3,-7.03E-4,0.,-2.38E-3/
        data phe/-8.95E-1,6.1,5.39,0.,8.01,4*0.,1200.,4*0.,-1.04E-5,
     &          1.9E-3,9.53E-4,1.06E-3,0.,-3.44E-3,10*0./ 
c       data pno/-22.4,17.7,-13.4,-4.88,62.3,32.7,0.,19.8,2.07,115.,
c    &          5*0.,3.94E-3,0.,2.48E-3,2.15E-4,6.67E-3,5*0.,
c    &          -8.4E-3,0.,-3.64E-3,2.E-3,-2.59E-2/
c       data po2/8.,-12.2,9.9,5.8,53.4,-25.2,0.,-28.5,-6.72,120.,
c    &          5*0.,-1.4E-2,0.,-9.3E-3,3.3E-3,2.8E-2,5*0.,4.25E-3,
c    &          0.,-6.04E-3,3.85E-3,-3.64E-2/
c       data pcl/4*0.,100.,4*0.,75.,10*0.,4*0.,-9.04E-3,-7.28E-3,
c    &          2*0.,3.46E-3,-2.11E-2/

        z=zd*umr
        f=fd*umr

        DO 8 I=1,5
        DO 8 J=1,6
                p(i,j,1)=po(i,j)
                p(i,j,2)=ph(i,j)
                p(i,j,3)=pn(i,j)
                p(i,j,4)=phe(i,j)
c               p(i,j,5)=pno(i,j)
c               p(i,j,6)=po2(i,j)
c               p(i,j,7)=pcl(i,j)
8       continue

        s=0.
c       do 5 i=1,7
        do 5 i=1,4
          do 7 j=1,6
                var(j) = p(1,j,i)*cos(z) + p(2,j,i)*cos(f) +
     &                   p(3,j,i)*cos(0.013*(300.-fs)) +
     &                   p(4,j,i)*cos(0.52*(t-6.)) + p(5,j,i)
7         continue
          cm(i)  = var(1)
          hm(i)  = var(2)
          all(i) = var(3)
          betl(i)= var(4)
          alh(i) = var(5)
          beth(i)= var(6)
          hx=h-hm(i)
          if(hx) 1,2,3
1               arg = hx * (hx * all(i) + betl(i)) 
                cn(i) = 0.
                if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
                goto 4
2               cn(i) = cm(i)
                goto 4
3               arg = hx * (hx * alh(i) + beth(i)) 
                cn(i) = 0.
                if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
4         continue
          if(cn(i).LT.0.005*cm(i)) cn(i)=0.
          if(cn(i).GT.cm(i)) cn(i)=cm(i)
          s=s+cn(i)
5       continue
c       do 6 i=1,7
        do 6 i=1,4
6               cn(i)=cn(i)/s*100.
        return
        end
C
C
      Subroutine ionco2(hei,xhi,it,F,R1,R2,R3,R4)
*----------------------------------------------------------------
*     INPUT:
*      hei  -  altitude in km
*      xhi  -  solar zenith angle in degree
*      it   -  seasonal month (Northern Hemisphere January 
*              is ismo=1 and so is Southern Hemisohere July)
*      F    -  10.7cm solar radio flux (12-month running mean)
*     OUTPUT:
*     R1 -  NO+ concentration (in percent)
*     R2 -  O2+ concentration (in percent) 
*     R3 -  Cb+ concentration (in percent) 
*     R4 -  O+  concentration (in percent) 
*
*  A.D. Danilov and N.V. Smirnova, Improving the 75 to 300 km ion 
*  composition model of the IRI, Adv. Space Res. 15, #2, 171-177, 1995.
*
*-----------------------------------------------------------------
      dimension j1ms70(7),j2ms70(7),h1s70(13,7),h2s70(13,7),
     *       R1ms70(13,7),R2ms70(13,7),rk1ms70(13,7),rk2ms70(13,7),
     *       j1ms140(7),j2ms140(7),h1s140(13,7),h2s140(13,7), 
     *       R1ms140(13,7),R2ms140(13,7),rk1ms140(13,7),rk2ms140(13,7),
     *       j1mw70(7),j2mw70(7),h1w70(13,7),h2w70(13,7),
     *       R1mw70(13,7),R2mw70(13,7),rk1mw70(13,7),rk2mw70(13,7),
     *       j1mw140(7),j2mw140(7),h1w140(13,7),h2w140(13,7), 
     *       R1mw140(13,7),R2mw140(13,7),rk1mw140(13,7),rk2mw140(13,7),
     *       j1mr70(7),j2mr70(7),h1r70(13,7),h2r70(13,7),
     *       R1mr70(13,7),R2mr70(13,7),rk1mr70(13,7),rk2mr70(13,7),
     *       j1mr140(7),j2mr140(7),h1r140(13,7),h2r140(13,7), 
     *       R1mr140(13,7),R2mr140(13,7),rk1mr140(13,7),rk2mr140(13,7)
      data j1ms70/11,11,10,10,11,9,11/ 
      data j2ms70/13,11,10,11,11,9,11/
      data h1s70/75,85,90,95,100,120,130,200,220,250,270,0,0,
     *        75,85,90,95,100,120,130,200,220,250,270,0,0, 
     *        75,85,90,95,100,115,200,220,250,270,0,0,0,
     *        75,80,95,100,120,140,200,220,250,270,0,0,0,
     *        75,80,95,100,120,150,170,200,220,250,270,0,0,
     *        75,80,95,100,140,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,110,145,200,220,250,270,0,0/
      data h2s70/75,80,90,95,100,120,130,140,150,200,220,250,270,
     *        75,80,90,95,100,120,130,200,220,250,270,0,0, 
     *        75,80,90,95,100,115,200,220,250,270,0,0,0,
     *        75,80,95,100,120,140,150,200,220,250,270,0,0,
     *        75,80,95,100,120,150,170,200,220,250,270,0,0,
     *        75,80,95,100,140,200,220,250,270,0,0,0,0,
     *        75,80,90,95,100,110,145,200,220,250,270,0,0/
      data R1ms70/6,30,60,63,59,59,66,52,20,4,2,0,0,
     *         6,30,60,63,69,62,66,52,20,4,2,0,0, 
     *         6,30,60,63,80,68,53,20,4,2,0,0,0,
     *         4,10,60,85,65,65,52,25,12,4,0,0,0, 
     *         4,10,60,89,72,60,60,52,30,20,10,0,0, 
     *         4,10,60,92,68,54,40,25,13,0,0,0,0, 
     *         1,8,20,60,95,93,69,65,45,30,20,0,0/ 
      data R2ms70/4,10,30,32,41,41,32,29,34,28,15,3,1,
     *         4,10,30,32,31,38,32,28,15,3,1,0,0,
     *         4,10,30,32,20,32,28,15,3,1,0,0,0,
     *         2,6,30,15,35,30,34,26,19,8,3,0,0,
     *         2,6,30,11,28,38,29,29,25,12,5,0,0,
     *         2,6,30,8,32,30,20,14,8,0,0,0,0,
     *         1,2,10,20,5,7,31,23,18,15,10,0,0/ 
      data rk1ms70/2.4,6.,.6,-.8,0,.7,-.2,-1.6,-.533,-.1,-.067,0,0,
     *         2.4,6.,.6,1.2,-.35,.4,-.2,-1.6,-.533,-.1,-.067,0,0, 
     *         2.4,6.,.6,3.4,-.8,-.176,-1.65,-.533,-.1,-.067,0,0,0,
     *         1.2,3.333,5.,-1.,0,-.216,-1.35,-.433,-.4,-.1,0,0,0,
     *         1.2,3.333,5.8,-.85,-.4,0,-.267,-1.1,-.333,-.4,-.2,0,0, 
     *         1.2,3.333,6.4,-.6,-.233,-.7,-.5,-.6,-.267,0,0,0,0, 
     *         1.4,2.4,4.,7.,-.2,-.686,-.072,-1.,-.5,-.5,-.5,0,0/
      data rk2ms70/1.2,2.,.4,1.8,0,-.9,-.3,.5,-.12,-.65,-.4,-.1,-.033,
     *         1.2,2.,.4,-.2,.35,-.6,-.057,-.65,-.4,-.1,-.033,0,0,
     *         1.2,2.,.4,-2.4,.8,-.047,-.65,-.4,-.1,-.033,0,0,0,
     *         .8,1.6,-3.,1.,-.25,.4,-.16,-.35,-.367,-.25,-.1,0,0,
     *         .8,1.6,-3.8,.85,.333,-.45,0,-.2,-.433,-.35,-.1,0,0,
     *         .8,1.6,-4.4,.6,-.033,-.5,-.2,-.3,-.2,0,0,0,0,
     *         .2,.8,2.,-3.,.2,.686,-.145,-.25,-.1,-.25,-.2,0,0/
      data j1ms140/11,11,10,10,9,9,12/ 
      data j2ms140/11,11,10,9,10,10,12/
      data h1s140/75,85,90,95,100,120,130,140,200,220,250,0,0,
     *        75,85,90,95,100,120,130,140,200,220,250,0,0,
     *        75,85,90,95,100,120,140,200,220,250,0,0,0,
     *        75,80,95,100,120,140,200,220,250,270,0,0,0,
     *        75,80,95,100,120,200,220,250,270,0,0,0,0,
     *        75,80,95,100,130,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,110,140,180,200,220,250,270,0/
      data h2s140/75,80,90,95,100,120,130,155,200,220,250,0,0,
     *        75,80,90,95,100,120,130,160,200,220,250,0,0,
     *        75,80,90,95,100,120,165,200,220,250,0,0,0,
     *        75,80,95,100,120,180,200,250,270,0,0,0,0,
     *        75,80,95,100,120,160,200,220,250,270,0,0,0,
     *        75,80,95,100,130,160,200,220,250,270,0,0,0,
     *        75,80,90,95,100,110,140,180,200,220,250,270,0/
      data R1ms140/6,30,60,63,59,59,66,66,38,14,1,0,0,
     *         6,30,60,63,69,62,66,66,38,14,1,0,0,
     *         6,30,60,63,80,65,65,38,14,1,0,0,0,
     *         4,10,60,85,66,66,38,22,9,1,0,0,0,
     *         4,10,60,89,71,42,26,17,10,0,0,0,0,
     *         4,10,60,93,71,48,35,22,10,0,0,0,0,
     *         1,8,20,60,95,93,72,60,58,40,26,13,0/ 
      data R2ms140/4,10,30,32,41,41,30,30,10,6,1,0,0,
     *         4,10,30,32,31,38,31,29,9,6,1,0,0,
     *         4,10,30,32,20,35,26,9,6,1,0,0,0,
     *         2,6,30,15,34,24,10,5,1,0,0,0,0,
     *         2,6,30,11,28,37,21,14,8,5,0,0,0,
     *         2,6,30,7,29,36,29,20,13,5,0,0,0,
     *         1,2,10,20,5,7,28,32,28,20,14,7,0/ 
      data rk1ms140/2.4,6.,.6,-.8,0,.7,0,-.467,-1.2,-.433,0,0,0,
     *         2.4,6.,.6,1.2,-.35,.4,0,-.467,-1.2,-.433,0,0,0,    
     *         2.4,6.,.6,3.4,-.75,0,-.45,-1.2,-.433,0,0,0,0,
     *         1.2,3.333,5.,-.95,0,-.467,-.8,-.433,-.4,0,0,0,0,
     *         1.2,3.333,5.8,-.9,-.363,-.8,-.3,-.35,-.3,0,0,0,0,
     *         1.2,3.333,6.6,-.733,-.329,-.65,-.433,-.6,-.267,0,0,0,0,
     *         1.4,2.4,4.,7.,-.2,-.7,-.3,-.1,-.9,-.467,-.65,-.333,0/
      data rk2ms140/1.2,2.,.4,1.8,0,-1.1,0,-.444,-.2,-.166,0,0,0,
     *         1.2,2.,.4,-.2,.35,-.7,-.067,-.5,-.15,-.166,0,0,0,
     *         1.2,2.,.4,-2.4,.75,-.2,-.486,-.15,-.166,0,0,0,0,
     *         .8,1.6,-3.,.95,-.167,-.7,-.1,-.2,0,0,0,0,0,
     *         .8,1.6,-3.8,.85,.225,-.4,-.35,-.2,-.15,-.133,0,0,0,
     *         .8,1.6,-4.6,.733,.233,-.175,-.45,-.233,-.4,-.1,0,0,0, 
     *         .2,.8,2.,-3.,.2,.7,.1,-.2,-.4,-.2,-.35,-.167,0/
      data j1mr70/12,12,12,9,10,11,13/ 
      data j2mr70/9,9,10,13,12,11,11/
      data h1r70/75,80,90,95,100,120,140,180,200,220,250,270,0,
     *        75,80,90,95,100,120,145,180,200,220,250,270,0, 
     *        75,80,90,95,100,120,145,180,200,220,250,270,0,  
     *        75,95,100,110,140,180,200,250,270,0,0,0,0,
     *        75,95,125,150,185,195,200,220,250,270,0,0,0,
     *        75,95,100,150,160,170,190,200,220,250,270,0,0,
     *        75,80,85,95,100,140,160,170,190,200,220,250,270/
      data h2r70/75,95,100,120,180,200,220,250,270,0,0,0,0,
     *        75,95,100,120,180,200,220,250,270,0,0,0,0, 
     *        75,95,100,120,130,190,200,220,250,270,0,0,0, 
     *        75,80,85,95,100,110,130,180,190,200,220,250,270,
     *        75,80,85,95,100,125,150,190,200,220,250,270,0,
     *        75,80,85,95,100,150,190,200,220,250,270,0,0, 
     *        75,85,95,100,140,180,190,200,220,250,270,0,0/
      data R1mr70/13,17,57,57,30,53,58,38,33,14,6,2,0,
     *         13,17,57,57,37,56,56,38,33,14,6,2,0, 
     *         13,17,57,57,47,58,55,37,33,14,6,2,0, 
     *         5,65,54,58,58,38,33,9,1,0,0,0,0, 
     *         5,65,65,54,40,40,45,26,17,10,0,0,0,    
     *         5,65,76,56,57,48,44,51,35,22,10,0,0, 
     *         3,11,35,75,90,65,63,54,54,50,40,26,13/ 
      data R2mr70/7,43,70,47,15,17,10,4,0,0,0,0,0,
     *         7,43,63,44,17,17,10,4,0,0,0,0,0, 
     *         7,43,53,42,42,13,17,10,4,0,0,0,0,
     *         3,5,26,34,46,42,41,23,16,16,10,1,0,
     *         3,5,26,34,35,35,42,25,22,14,8,5,0, 
     *         3,5,26,34,24,41,31,26,20,13,5,0,0,
     *         3,15,15,10,35,35,30,34,20,14,7,0,0/ 
      data rk1mr70/.8,4.,0,-5.4,1.15,.25,-.5,-.25,-.95,-.267,-.2,
     *             -.067,0,
     *         .8,4.,0,-4.,.95,0,-.514,-.25,-.95,-.267,-.2,-.067,0, 
     *         .8,4.,0,-2.,.55,-.12,-.514,-.2,-.95,-.267,-.2,-.067,0, 
     *         3.,-2.2,.4,0,-.5,-.25,-.48,-.4,-.033,0,0,0,0,   
     *         3.,0,-.44,-.466,0,1.0,-.95,-.3,-.35,-.3,0,0,0, 
     *         3.,2.2,-.4,0.1,-.9,-.2,.7,-.8,-.433,-.6,-.267,0,0, 
     *         1.6,4.8,4.,3.,-.625,-.1,-.9,0,-.4,-.5,-.467,-.65,-.3/
      data rk2mr70/1.8,5.4,-1.15,-.533,.1,-.35,-.2,-.2,0,0,0,0,0,
     *         1.8,4.,-.95,-.45,0,-.35,-.2,-.2,0,0,0,0,0,
     *         1.8,2.,-.55,0,-.483,.4,-.35,-.2,-.2,0,0,0,0,    
     *         .4,4.2,.8,2.4,-.4,-.05,-.36,-.7,0,-.3,-.3,-.05,0,
     *         .4,4.2,.8,.2,0,.28,-.425,-.3,-.4,-.2,-.15,-.133,0,  
     *         .4,4.2,.8,-2.,.34,-.25,-.5,-.3,-.233,-.4,-.1,0,0,  
     *         1.2,0,-1.,.625,0,-.5,.4,-.7,-.2,-.35,-.167,0,0/
      data j1mr140/12,12,11,12,9,9,13/ 
      data j2mr140/10,9,10,12,13,13,12/
      data h1r140/75,80,90,95,100,115,130,145,200,220,250,270,0,
     *        75,80,90,95,100,110,120,145,200,220,250,270,0, 
     *        75,80,90,95,100,115,150,200,220,250,270,0,0,
     *        75,95,100,120,130,140,150,190,200,220,250,270,0,
     *        75,95,120,150,190,200,220,250,270,0,0,0,0,
     *        75,95,100,145,190,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,120,160,170,190,200,220,250,270/
      data h2r140/75,95,100,115,130,175,200,220,250,270,0,0,0,
     *        75,95,100,110,175,200,220,250,270,0,0,0,0, 
     *        75,95,100,115,130,180,200,220,250,270,0,0,0, 
     *        75,80,85,95,100,120,130,190,200,220,250,270,0,
     *        75,80,85,95,100,120,140,160,190,200,220,250,270,
     *        75,80,85,95,100,145,165,180,190,200,220,250,270,
     *        75,85,95,100,120,145,170,190,200,220,250,270,0/
      data R1mr140/13,17,57,57,28,51,56,56,12,8,1,0,0,
     *         13,17,57,57,36,46,55,56,10,8,1,0,0, 
     *         13,17,57,57,46,56,55,12,8,1,0,0,0,
     *         5,65,54,59,56,56,53,23,16,13,3,1,0,
     *         5,65,65,54,29,16,16,10,2,0,0,0,0,
     *         5,65,76,58,36,25,20,12,7,0,0,0,0,
     *         3,11,35,75,91,76,58,49,45,32,28,20,12/ 
      data R2mr140/7,43,72,49,44,14,7,4,1,0,0,0,0,
     *         7,43,64,51,14,7,4,1,0,0,0,0,0, 
     *         7,43,54,44,44,13,7,4,1,0,0,0,0,
     *         3,5,26,34,46,41,44,9,11,7,2,1,0,
     *         3,5,26,34,35,35,40,40,16,14,9,5,2, 
     *         3,5,26,34,24,40,40,32,19,20,10,7,3,
     *         3,15,15,9,24,35,40,28,28,20,10,8,0/ 
      data rk1mr140/.8,4.,0,-5.8,1.533,.333,0,-.8,-.2,-.233,-.05,0,0,
     *         .8,4.,0,-4.2,1.3,.6,.04,-.836,-.1,-.233,-.05,0,0,
     *         .8,4.,0,-2.2,.667,-.029,-.86,-.2,-.233,-.05,0,0,0,         
     *         3.,-2.2,.25,-.3,0,-.3,-.75,-.7,-.15,-.333,-.1,-.033,0,
     *         3.,0,-.367,-.625,-1.3,0,-.2,-.4,-.067,0,0,0,0, 
     *         3.,2.2,-.4,-.489,-1.1,-.25,-.267,-.25,-.2,0,0,0,0,
     *         1.6,4.8,4.,3.2,-.75,-.45,-.9,-.2,-1.3,-.2,-.267,-.4,-.3/
      data rk2mr140/1.8,5.8,-1.533,-.333,-.667,-.28,-.15,-.1,-.05,
     *              0,0,0,0,
     *         1.8,4.2,-1.3,-.569,-.28,-.15,-.1,-.05,0,0,0,0,0,
     *         1.8,2.2,-.667,0,-.62,-.3,-.15,-.1,-.05,0,0,0,0,    
     *         .4,4.2,.8,2.4,-.25,.3,-.583,.2,-.2,-.167,-.05,-.033,0,
     *         .4,4.2,.8,.02,0,.25,0,-.6,-.2,-.25,-.133,-.15,-.067,  
     *         .4,4.2,.8,-2.,.356,0,-.533,-1.3,.1,-.5,-.1,-.2,-.1,  
     *         1.2,0,-1.2,.75,.44,.2,-.6,0,-.4,-.333,-.1,-.2,0/
      data j1mw70/13,13,13,13,9,8,9/ 
      data j2mw70/10,10,11,11,9,8,11/
      data h1w70/75,80,85,95,100,110,125,145,180,200,220,250,270,
     *        75,80,85,95,100,110,120,150,180,200,220,250,270,
     *        75,80,85,95,100,110,120,155,180,200,220,250,270,
     *        75,80,90,100,110,120,140,160,190,200,220,250,270, 
     *        75,80,90,110,150,200,220,250,270,0,0,0,0,
     *        75,80,90,100,150,200,250,270,0,0,0,0,0,
     *        75,80,90,100,120,130,140,200,270,0,0,0,0/
      data h2w70/75,90,95,100,110,125,190,200,250,270,0,0,0,
     *        75,90,95,100,110,125,190,200,250,270,0,0,0, 
     *        75,90,95,100,110,120,145,190,200,250,270,0,0,
     *        75,80,95,100,110,120,150,200,220,250,270,0,0,
     *        75,80,90,95,110,145,200,250,270,0,0,0,0,
     *        75,80,90,100,140,150,200,250,0,0,0,0,0,
     *        75,80,85,90,100,120,130,140,160,200,270,0,0/ 
      data R1mw70/28,35,65,65,28,44,46,50,25,25,10,5,0,
     *         28,35,65,65,36,49,47,47,25,25,10,5,0,
     *         28,35,65,65,48,54,51,43,25,25,10,5,0,
     *         16,24,66,54,58,50,50,38,25,25,10,5,0, 
     *         16,24,66,66,46,30,20,6,3,0,0,0,0,  
     *         16,24,66,76,49,32,12,7,0,0,0,0,0,      
     *         6,19,67,91,64,68,60,40,12,0,0,0,0/ 
      data R2mw70/5,35,35,72,56,54,12,12,2,0,0,0,0,
     *         5,35,35,64,51,53,12,12,2,0,0,0,0, 
     *         5,35,35,52,46,49,41,12,12,2,0,0,0,
     *         4,10,40,46,42,50,41,12,7,2,0,0,0,
     *         4,10,30,34,34,51,14,4,2,0,0,0,0,
     *         4,10,30,24,45,48,20,5,0,0,0,0,0,
     *         2,6,17,23,9,36,32,40,40,20,6,0,0/ 
      data rk1mw70/1.4,6.,0,-7.4,1.6,.133,.2,-.714,0,-.75,-.167,-.25,0,
     *         1.4,6.,0,-5.8,1.3,-.2,0,-.733,0,-.75,-.167,-.25,0,
     *         1.4,6.,0,-3.4,.6,-.3,-.229,-.72,0,-.75,-.167,-.25,0,
     *         1.6,4.2,-1.2,.4,-.8,0,-.6,-.433,0,-.75,-.167,-.25,0,   
     *         1.6,4.2,0,-.5,-.32,-.5,-.467,-.15,-.1,0,0,0,0,
     *         1.6,4.2,1.,-.54,-.34,-.4,-.25,-.2,0,0,0,0,0,      
     *         2.6,4.8,2.4,-1.35,.4,-.8,-.333,-.4,-.3,0,0,0,0/
      data rk2mw70/2.,0,7.4,-1.6,-.133,-.646,0,-.2,-.1,0,0,0,0,
     *         2.,0,5.8,-1.3,.133,-.631,0,-.2,-.1,0,0,0,0,    
     *         2.,0,3.4,-.6,.3,-.32,-.644,0,-.2,-.1,0,0,0,
     *         1.2,2.,1.2,-.4,.8,-.3,-.58,-.25,-.167,-.1,0,0,0,
     *         1.2,2.,.8,0,.486,-.673,-.2,-.1,-.066,0,0,0,0,  
     *         1.2,2.,-.6,.525,.3,-.56,-.3,-.1,0,0,0,0,0,  
     *         .8,2.2,1.2,-1.4,1.35,-.4,.8,0,-.5,-.2,-.167,0,0/
      data j1mw140/12,11,11,11,11,10,12/ 
      data j2mw140/10,11,11,11,11,10,12/
      data h1w140/75,80,85,95,100,110,125,145,190,200,220,250,0,
     *        75,80,85,95,100,110,120,150,190,220,250,0,0,
     *        75,80,85,95,100,110,120,155,190,220,250,0,0,
     *        75,80,90,100,110,120,140,160,190,220,250,0,0,
     *        75,80,90,110,150,160,190,200,220,250,270,0,0,
     *        75,80,90,100,150,160,190,200,250,270,0,0,0,
     *        75,80,90,100,120,130,140,160,190,200,250,270,0/
      data h2w140/75,90,95,100,110,125,190,200,220,250,0,0,0,
     *        75,90,95,100,110,120,125,190,200,220,250,0,0,
     *        75,90,95,100,110,120,145,190,200,220,250,0,0,  
     *        75,80,95,100,110,120,150,190,200,220,250,0,0,
     *        75,80,90,95,110,145,190,200,220,250,270,0,0,
     *        75,80,90,100,140,150,200,220,250,270,0,0,0,
     *        75,80,85,90,100,120,130,140,160,180,200,220,0/ 
      data R1mw140/28,35,65,65,28,44,46,50,9,6,2,0,0,
     *         28,35,65,65,36,49,47,47,8,2,0,0,0,
     *         28,35,65,65,48,54,51,43,8,2,0,0,0,
     *         16,24,66,54,58,50,50,42,8,2,0,0,0, 
     *         16,24,66,66,46,49,9,10,7,2,0,0,0,
     *         16,24,66,76,49,54,10,14,4,1,0,0,0,
     *         6,19,67,91,64,68,60,58,11,20,5,2,0/ 
      data R2mw140/5,35,35,72,56,54,5,5,1,0,0,0,0,
     *         5,35,35,64,51,53,53,5,5,1,0,0,0,
     *         5,35,35,52,46,49,41,5,5,1,0,0,0,    
     *         4,10,40,46,42,50,41,5,5,1,0,0,0,   
     *         4,10,30,34,34,51,10,5,3,1,0,0,0,  
     *         4,10,30,24,45,48,4,2,1,0,0,0,0,
     *         2,6,17,23,9,36,32,40,39,29,1,0,0/ 
      data rk1mw140/1.4,6.,0,-7.4,1.6,.133,.2,-.911,-.3,-.2,-.066,0,0,
     *         1.4,6.,0,-5.8,1.3,-.2,0,-.975,-.2,-.066,0,0,0,
     *         1.4,6.,0,-3.4,.6,-.3,-.229,-1.,-.2,-.066,0,0,0,
     *         1.6,4.2,-1.2,.4,-.8,0,-.4,-1.133,-.2,-.066,0,0,0,
     *         1.6,4.2,0,-.5,.3,-1.133,.1,-.15,-.166,-.1,0,0,0,
     *         1.6,4.2,1.,-.54,.5,-1.466,.4,-.2,-.15,-.0333,0,0,0,
     *         2.6,4.8,2.4,-1.35,.4,-.8,-.1,-1.566,.9,-.3,-.15,-.05,0/
      data rk2mw140/2.,0,7.4,-1.6,-.133,-.754,0,-.2,-.033,0,0,0,0,
     *         2.,0,5.8,-1.3,.2,0,-.738,0,-.2,-.033,0,0,0,
     *         2.,0,3.4,-.6,.3,-.32,-.8,0,-.2,-.033,0,0,0,
     *         1.2,2.,1.2,-.4,.8,-.3,-.9,0,-.2,-.033,0,0,0,
     *         1.2,2.,.8,0,.486,-.911,-.5,-.1,-.066,-.05,0,0,0,
     *         1.2,2.,-.6,.525,.3,-.88,-.1,-.033,-.05,0,0,0,0,
     *         .8,2.2,1.2,-1.4,1.35,-.4,.8,-.05,-.5,-1.4,-.05,0,0/

        h = hei
        z = xhi

         if(z.lt.20)z=20
         if(z.gt.90)z=90
        if((it.eq.1).or.(it.eq.2).or.(it.eq.11).or.(it.eq.12))then
         if(f.lt.140)then
           Call aprok(j1mw70,j2mw70,h1w70,h2w70,R1mw70,R2mw70,        
     *                rk1mw70,rk2mw70,h,z,R1,R2) 
           R170=R1
           R270=R2
           endif
         if(f.gt.70)then
           Call aprok(j1mw140,j2mw140,h1w140,h2w140,R1mw140,R2mw140,        
     *                rk1mw140,rk2mw140,h,z,R1,R2) 
           R1140=R1
           R2140=R2
           endif
         if((f.gt.70).and.(f.lt.140))then
           R1=R170+(R1140-R170)*(f-70)/70
           R2=R270+(R2140-R270)*(f-70)/70
           endif
         endif
        if((it.eq.5).or.(it.eq.6).or.(it.eq.7).or.(it.eq.8))then
         if(f.lt.140)then
           Call aprok(j1ms70,j2ms70,h1s70,h2s70,R1ms70,R2ms70,        
     *                rk1ms70,rk2ms70,h,z,R1,R2) 
           R170=R1
           R270=R2
           endif
         if(f.gt.70)then
           Call aprok(j1ms140,j2ms140,h1s140,h2s140,R1ms140,R2ms140,        
     *                rk1ms140,rk2ms140,h,z,R1,R2) 
           R1140=R1
           R2140=R2
           endif
         if((f.gt.70).and.(f.lt.140))then
           R1=R170+(R1140-R170)*(f-70)/70
           R2=R270+(R2140-R270)*(f-70)/70
           endif
         endif
        if((it.eq.3).or.(it.eq.4).or.(it.eq.9).or.(it.eq.10))then
         if(f.lt.140)then
           Call aprok(j1mr70,j2mr70,h1r70,h2r70,R1mr70,R2mr70,        
     *                rk1mr70,rk2mr70,h,z,R1,R2) 
           R170=R1
           R270=R2
           endif
         if(f.gt.70)then
           Call aprok(j1mr140,j2mr140,h1r140,h2r140,R1mr140,R2mr140,
     *                rk1mr140,rk2mr140,h,z,R1,R2) 
           R1140=R1
           R2140=R2
           endif
         if((f.gt.70).and.(f.lt.140))then
           R1=R170+(R1140-R170)*(f-70)/70
           R2=R270+(R2140-R270)*(f-70)/70
           endif
         endif
        R3=0
        R4=0
        if (h.lt.100) R3=100-(R1+R2)
        if (h.ge.100) R4=100-(R1+R2)
         if(R3.lt.0) R3=0
         if(R4.lt.0) R4=0
        R1=ANINT(R1)
        R2=ANINT(R2)
        R3=ANINT(R3)
        R4=ANINT(R4)
 300   continue
        end
c
c
      Subroutine aprok(j1m,j2m,h1,h2,R1m,R2m,rk1m,rk2m,hei,xhi,R1,R2)
c----------------------------------------------------------------- 
      dimension   zm(7),j1m(7),j2m(7),h1(13,7),h2(13,7),R1m(13,7),
     *            R2m(13,7),rk1m(13,7),rk2m(13,7)
      data        zm/20,40,60,70,80,85,90/
      
        h=hei
        z=xhi

         j1=1
         j2=1
         i1=1
       do 1 i=1,7
         i1=i
        if(z.eq.zm(i)) j1=0
        if(z.le.zm(i)) goto 11
 1     continue
 11    continue
          i2=1
         do 2 i=2,j1m(i1)
          i2=i-1
          if(h.lt.h1(i,i1)) goto 22
          i2=j1m(i1)
 2       continue
 22      continue
          i3=1
         do 3 i=2,j2m(i1)
          i3=i-1
          if(h.lt.h2(i,i1)) goto 33
          i3=j2m(i1)
 3       continue
 33      continue
        R01=R1m(i2,i1)
        R02=R2m(i3,i1)
        rk1=rk1m(i2,i1)
        rk2=rk2m(i3,i1)
        h01=h1(i2,i1)
        h02=h2(i3,i1)
        R1=R01+rk1*(h-h01)
        R2=R02+rk2*(h-h02)
        if(j1.eq.1)then
          j1=0
          j2=0
          i1=i1-1
          R11=R1
          R12=R2
          goto 11
        endif
        if(j2.eq.0)then
          rk=(z-zm(i1))/(zm(i1+1)-zm(i1)) 
          R1=R1+(R11-R1)*rk
          R2=R2+(R12-R2)*rk
        endif
       end
C
C
	SUBROUTINE CALION(INVDIP,MLT,ALT,DDD,PF107,NO,NH,NHE,NN)
C----------------------------------------------------------------------
C Version 2.5 (released 6.4.2016) (Includes a correction for PF107<87.5 
C based on C/NOFS CINDI data)
C CALION calculates ion density of O+, H+, He+ and N+  in the outer
C ionosphere with regard to solar activity (F107 index).
C CALION uses subroutines IONLOW and IONHIGH.
C Linearly interpolates the logarithmic ion densities in solar activity.
C Inputs: INVDIP - "mix" coordinate of the dip latitude and of
C                    the invariant latitude;
C                    positive northward, in deg, range <-90.0;90.0>
C         MLT - magnetic local time (central dipole)
C               in hours, range <0;24)
C         ALT - altitude above the Earth's surface;
C               in km, range <350;2000>
C         DDD - day of year; range <0;365>
C         PF107 - PF107 index
C Output: NO,NH,NHE,NN - relative density of O+, H+, He+, and N+ in %
C Versions:  2.5 FORTRAN
C
C REFERENCES:
C   Triskova, L., Truhlik, V., Smilauer, J. An empirical model of ion
C      composition in the outer ionosphere. Adv. Space Res. 31 (3), 
C      653–663, 2003.
C   Truhlik, V., Triskova, L., Smilauer, J. New advances in empirical
C      modeling of ion composition in the outer ionosphere. Adv. Space
C      Res. 33, 844–849, 2004.
C   Truhlik V., D. Bilitza, L. Triskova, Towards better description of 
C     solar activity variation in the International Reference Ionosphere 
C     topside ion composition model, Advances in Space Research, 
C     Volume 55, Issue 8, 15 April 2015, Pages 2099-2105, ISSN 0273-1177, 
C     http://dx.doi.org/10.1016/j.asr.2014.07.033.
C
C Author of the code:
C         Vladimir Truhlik
C         Institute of Atm. Phys.
C         Bocni II.
C         141 31 Praha 4, Sporilov
C         Czech Republic
C         e-mail: vtr@ufa.cas.cz
C         tel/fax: +420 267103058, +420 728073539 / +420 272 762528
C----------------------------------------------------------------------
      REAL INVDIP,MLT,ALT,PF107
      INTEGER DDD,ION
      REAL NO,NH,NHE,NN,NOH,NHH,NHEH,NNH,NOL,NHL,NHEL,NNL,NTOT
      REAL NOCORR,NHCORR
      DIMENSION  DOL(4,3,49),DHL(4,3,49),DHEL(4,3,49),DNL(4,3,49)
      DIMENSION  DOH(4,3,49),DHH(4,3,49),DHEH(4,3,49),DNH(4,3,49)
      DIMENSION  CORRO(3),CORRH(3)
C/////////////////////coefficients high solar activity////////////////////////
C//////////////////////////////////O+/////////////////////////////////////////
C     550km equinox
      DATA (DOH(1,1,J),J=1,49)/ 1.1638E+001, 2.1979E-006,-2.3248E-001,
     &                       -3.3214E-007, 1.1898E-001, 3.0748E-006,
     &                        2.6115E-003,-2.4560E-001, 1.2462E-007,
     &                        2.0023E-002, 1.6400E-007, 6.2374E-003,
     &                       -6.4331E-008,-2.6353E-001, 3.3890E-008,
     &                        3.6762E-002, 1.8064E-007,-2.2751E-002,
     &                        2.8274E-008, 1.1224E-001,-1.9642E-007,
     &                       -8.7787E-003,-6.0584E-009,-4.6188E-003,
     &                        1.0754E-002, 1.0973E-007, 8.7284E-003,
     &                        5.2698E-008, 5.1970E-003,-6.6284E-002,
     &                       -5.0381E-007, 2.5389E-003, 1.7134E-007,
     &                        2.2995E-002, 4.7651E-007, 4.8065E-004,
     &                       -1.5021E-008, 2.1538E-003,-1.8609E-007,
     &                        5.0122E-003, 3.2356E-003, 2.1760E-008,
     &                        9.9773E-004,-1.6112E-003,-1.7006E-007,
     &                        4.5044E-003, 8.6730E-008, 1.4997E-003,
     &                       -4.2748E-003/
C     550km June solstice
      DATA (DOH(1,2,J),J=1,49)/ 1.1466E+001, 4.0119E-001,-4.1736E-001,
     &                       -1.9200E-001, 3.7017E-002,-2.6320E-002,
     &                       -7.3409E-002,-1.6393E-001,-7.5744E-003,
     &                        4.7543E-002,-5.9249E-003,-9.1045E-003,
     &                        7.5707E-003,-6.6604E-002, 3.3164E-002,
     &                        8.8713E-003, 6.9505E-003,-3.7147E-003,
     &                       -9.4076E-003,-4.7440E-002, 5.1116E-004,
     &                       -2.0023E-002, 1.1603E-002, 5.8822E-003,
     &                       -1.8973E-002,-3.2266E-004,-1.0472E-002,
     &                        8.7367E-003, 8.7078E-004, 3.0624E-002,
     &                       -1.4384E-002, 1.1849E-003, 1.6199E-003,
     &                        5.2375E-002,-5.1855E-003, 5.0677E-003,
     &                        5.8491E-004,-2.6023E-002, 2.1862E-003,
     &                        1.7057E-004, 1.4512E-003, 9.2850E-003,
     &                       -1.8361E-003,-1.5087E-002, 1.6268E-003,
     &                       -5.1916E-003, 3.1410E-003,-7.5738E-003,
     &                        2.0747E-003/
C     900km equinox
      DATA (DOH(2,1,J),J=1,49)/ 1.0822E+001, 1.0243E-006,-1.1612E-002,
     &                       -3.9090E-006, 1.7949E-001, 5.2434E-006,
     &                       -5.5277E-002,-4.2250E-001, 2.2743E-007,
     &                       -1.7850E-002,-2.1446E-007, 2.0109E-002,
     &                        2.5441E-008,-2.8437E-001, 2.8675E-007,
     &                        4.6531E-002,-1.7470E-007,-1.2887E-002,
     &                       -7.0467E-008,-1.1997E-001, 2.0672E-007,
     &                       -2.1904E-002, 3.2349E-008, 1.4037E-002,
     &                        5.9543E-002, 4.1880E-007,-3.2790E-004,
     &                       -9.1809E-008, 8.8505E-003,-1.0878E-002,
     &                       -4.1668E-007,-1.6918E-002, 4.1870E-008,
     &                       -3.1914E-002, 2.2837E-007,-1.0269E-002,
     &                       -3.4403E-008, 6.6131E-004,-1.4563E-007,
     &                       -6.0043E-003,-1.7662E-002,-1.6897E-009,
     &                       -7.6391E-003,-3.0607E-003,-1.5907E-007,
     &                       -8.2787E-003, 2.1356E-009,-7.6774E-003,
     &                        1.2753E-002/
C     900km June solstice
      DATA (DOH(2,2,J),J=1,49)/ 1.0620E+001, 4.1405E-001, 4.7224E-002,
     &                       -1.7052E-001, 2.0812E-002,-1.6568E-002,
     &                       -1.7739E-002,-3.6758E-001, 1.2945E-001,
     &                       -5.7690E-002, 1.6517E-002, 4.3426E-004,
     &                       -8.4872E-003,-1.7111E-001,-4.3498E-002,
     &                        9.6720E-002, 6.6361E-003,-1.7696E-002,
     &                       -7.3752E-003,-5.9474E-002, 7.0600E-002,
     &                       -2.1261E-002,-4.5486E-003,-2.3596E-003,
     &                       -6.0694E-002,-3.4020E-002, 2.5395E-002,
     &                       -1.0267E-002, 5.1524E-003, 2.4349E-002,
     &                        1.2764E-002, 4.7153E-003, 3.0511E-003,
     &                        7.1283E-002, 7.9275E-003, 9.1649E-003,
     &                       -1.2318E-003,-4.2938E-002,-3.0165E-003,
     &                        2.9838E-003,-1.6615E-002, 1.2946E-003,
     &                        2.7673E-003,-3.7213E-002, 6.0640E-004,
     &                       -2.7206E-002, 1.4542E-003,-2.8307E-002,
     &                       -7.3467E-003/
C     1500km equinox
      DATA (DOH(3,1,J),J=1,49)/ 9.7937E+000,-1.9164E-006, 3.5453E-001,
     &                        1.6223E-006, 2.0156E-001,-1.0210E-006,
     &                       -1.1503E-001,-8.2204E-001, 5.7192E-008,
     &                       -4.8372E-002, 9.9279E-008, 6.0963E-002,
     &                        4.4030E-008,-4.0169E-001,-1.2124E-006,
     &                        1.1542E-001, 7.4444E-007,-1.6794E-002,
     &                       -7.7795E-008,-2.3478E-001, 1.8889E-006,
     &                       -5.7216E-002,-8.5613E-008, 2.2638E-003,
     &                       -1.8245E-001,-6.7741E-007, 5.6764E-004,
     &                        1.0843E-007, 8.6873E-003,-3.2154E-002,
     &                        1.2769E-006,-2.3455E-002, 2.3495E-008,
     &                        4.4324E-002, 5.8249E-007,-9.7873E-003,
     &                       -4.9479E-008,-4.3582E-002, 2.4400E-007,
     &                       -4.7593E-003, 1.1520E-002, 6.1115E-007,
     &                       -6.1079E-003, 4.6842E-003,-7.7660E-008,
     &                       -1.4542E-002, 1.3010E-007, 1.2669E-004,
     &                        2.6225E-002/
C     1500km June solstice
      DATA (DOH(3,2,J),J=1,49)/ 9.7559E+000, 7.3006E-001, 2.3915E-002,
     &                       -3.5479E-001, 2.3232E-001, 4.3105E-002,
     &                       -8.3276E-002,-6.1803E-001, 2.4574E-001,
     &                        1.0126E-002,-1.3819E-002, 5.2543E-002,
     &                       -2.4424E-002,-1.6943E-001, 4.6190E-002,
     &                        5.4503E-002,-4.7693E-002, 2.8382E-002,
     &                       -4.0051E-002, 3.6658E-002,-7.1019E-002,
     &                        3.3261E-002, 1.4489E-003, 2.5216E-003,
     &                       -2.6986E-003, 2.8257E-003, 1.1303E-002,
     &                        1.3446E-002, 6.4292E-004, 7.0927E-002,
     &                       -4.8257E-002, 1.1750E-002, 2.9293E-004,
     &                        3.0686E-002,-2.9274E-002,-1.4456E-004,
     &                        2.2379E-003,-1.6401E-002,-1.1951E-002,
     &                       -8.3522E-005,-2.2854E-003,-5.1928E-003,
     &                       -1.0387E-003,-2.4951E-002, 2.3471E-003,
     &                       -1.3592E-002,-2.2130E-003,-7.8338E-003,
     &                        7.0712E-004/
C     2250km equinox
      DATA (DOH(4,1,J),J=1,49)/ 9.1209E+000,-3.5060E-008, 5.6220E-001,
     &                        7.2825E-007, 1.8242E-001, 1.6496E-006,
     &                       -1.3481E-001,-6.2102E-001, 1.3313E-006,
     &                       -7.2196E-002, 2.8528E-007, 3.9029E-002,
     &                        5.3745E-008,-6.0644E-001, 3.8846E-007,
     &                        1.0864E-001,-2.9182E-007, 1.9555E-002,
     &                        1.1610E-008,-1.1778E-001, 6.6455E-009,
     &                       -4.3774E-002, 2.4623E-007,-7.3456E-003,
     &                       -1.0794E-001, 2.8097E-007,-6.1924E-003,
     &                        5.5108E-008, 6.5963E-003, 9.7996E-003,
     &                       -1.4108E-007, 1.3135E-002, 3.9136E-009,
     &                        2.0612E-001, 1.5759E-007, 2.9679E-003,
     &                        1.9150E-007, 5.5943E-002,-9.0442E-009,
     &                        8.7562E-003, 7.1924E-002,-8.2087E-008,
     &                        3.7554E-003, 3.5717E-003, 7.9697E-008,
     &                        7.6157E-003,-8.0663E-008,-3.9358E-003,
     &                        1.0457E-002/
C     2250km June solstice
      DATA (DOH(4,2,J),J=1,49)/ 8.9497E+000, 6.7881E-001, 3.9055E-001,
     &                       -3.2002E-001, 1.3908E-001, 2.8471E-002,
     &                       -1.7994E-001,-7.3377E-001, 1.1401E-001,
     &                       -1.6313E-002, 5.8219E-002,-2.4831E-003,
     &                       -1.8989E-002,-6.0576E-001, 1.2761E-001,
     &                        7.9504E-002,-5.2144E-002, 1.1163E-002,
     &                        6.6893E-004,-3.4924E-002, 1.6598E-002,
     &                       -7.1152E-003, 4.5419E-003,-2.5894E-003,
     &                       -1.7938E-001,-9.1538E-003,-5.5001E-003,
     &                        8.3203E-003, 1.2648E-003, 1.8489E-001,
     &                        2.9295E-003, 7.8076E-003,-6.2127E-003,
     &                       -6.9425E-002, 6.4065E-003,-3.1037E-003,
     &                        2.8284E-003,-1.5089E-002, 1.3666E-003,
     &                        2.2100E-003, 3.6826E-003, 9.4428E-003,
     &                       -2.4497E-003,-1.0100E-002,-4.1794E-003,
     &                       -2.4757E-002, 8.2641E-005,-9.0907E-003,
     &                       -5.7800E-003/
C//////////////////////////////////////////////////////////////////////
C//////////////////////////////////H+//////////////////////////////////
C     550km equinox
      DATA (DHH(1,1,J),J=1,49)/ 8.4599E+000,-9.9884E-007,-3.5443E-001,
     &                       -6.5810E-007,-2.9651E-002, 5.2904E-007,
     &                        9.3573E-003, 4.5701E-001, 1.9996E-007,
     &                       -2.4248E-002, 7.7567E-008, 1.1724E-002,
     &                       -2.3639E-008, 4.2217E-003,-1.1230E-007,
     &                        2.3361E-002, 1.2442E-007,-7.5842E-003,
     &                        2.9140E-008, 5.0338E-002,-5.2187E-008,
     &                        7.0107E-003,-6.2618E-008,-1.1794E-003,
     &                       -1.2422E-001,-7.0145E-008,-6.0802E-003,
     &                       -1.7967E-008,-3.5129E-003,-5.8675E-002,
     &                        1.8341E-008, 1.2289E-003,-1.1887E-008,
     &                        1.1607E-001, 2.2739E-008, 1.4341E-002,
     &                       -8.2486E-009,-3.2079E-002,-4.4285E-008,
     &                        2.1810E-003, 4.0119E-002, 3.4542E-009,
     &                        1.4022E-003,-9.4002E-003, 2.0700E-008,
     &                        6.9210E-003, 1.3041E-008,-3.3684E-004,
     &                        8.8583E-003/
C     550km June solstice
      DATA (DHH(1,2,J),J=1,49)/ 8.7052E+000,-1.0386E-002,-3.9422E-001,
     &                        2.2079E-002, 1.0061E-001,-1.0053E-001,
     &                       -3.8012E-002, 4.0303E-001, 1.1450E-001,
     &                       -3.6893E-002, 7.1546E-003, 4.8412E-003,
     &                       -3.1532E-003, 5.3452E-002,-7.1005E-002,
     &                        1.1385E-002, 2.0830E-004,-3.4358E-003,
     &                       -9.9962E-004, 1.0466E-001, 2.6587E-002,
     &                       -5.0469E-004, 1.0174E-002,-3.1445E-004,
     &                       -1.4247E-001, 1.0854E-002,-8.3255E-003,
     &                       -1.0450E-002, 3.3930E-003,-4.1298E-002,
     &                        8.3111E-003, 2.4336E-003,-2.3697E-003,
     &                       -2.1893E-002, 1.5365E-002, 3.0954E-003,
     &                       -5.5201E-003,-1.1731E-002, 1.1634E-002,
     &                       -1.4933E-003, 5.2614E-004, 7.8891E-003,
     &                       -1.0120E-003,-8.6891E-003, 8.2004E-003,
     &                        6.2137E-003, 4.2136E-003,-1.3120E-002,
     &                        9.2474E-003/
C     900km equinox
      DATA (DHH(2,1,J),J=1,49)/ 9.3870E+000,-1.3719E-007,-8.3508E-001,
     &                       -5.9282E-007, 6.1012E-002,-5.2218E-007,
     &                        1.9141E-001, 4.1174E-001, 3.8016E-007,
     &                       -2.5064E-002,-9.7075E-008,-3.3704E-002,
     &                       -3.0830E-008, 6.3980E-002, 5.1005E-007,
     &                       -2.2929E-002, 3.6584E-007,-1.2908E-003,
     &                       -4.1545E-008, 2.3400E-002,-3.5827E-008,
     &                        1.6850E-002, 1.4133E-007, 4.5264E-003,
     &                       -1.6435E-001, 1.7264E-007, 4.2052E-003,
     &                        1.6665E-007, 5.7113E-003, 2.6847E-002,
     &                        2.0150E-007, 7.1160E-004, 2.6940E-008,
     &                        5.9238E-002, 3.4505E-008, 2.3152E-003,
     &                       -1.4563E-010,-3.0566E-002, 9.2142E-009,
     &                        3.3032E-004,-2.0039E-002, 1.0352E-008,
     &                        1.9382E-003,-5.5282E-004,-6.9643E-008,
     &                       -4.2814E-003, 8.0598E-008, 1.1031E-002,
     &                       -2.6108E-002/
C     900km June solstice
      DATA (DHH(2,2,J),J=1,49)/ 9.4299E+000,-5.1759E-002,-8.1008E-001,
     &                        1.3248E-001, 1.4374E-001, 2.8595E-002,
     &                       -3.6756E-003, 2.4542E-001, 8.4647E-002,
     &                        2.6415E-002, 3.3343E-002,-1.5748E-002,
     &                       -7.9269E-003, 8.0631E-002, 1.0277E-002,
     &                       -1.8242E-002,-2.3267E-002,-1.0559E-003,
     &                        8.2585E-003, 2.1509E-002,-8.1102E-003,
     &                        3.9426E-003,-1.0183E-003, 6.8868E-003,
     &                       -5.1945E-002, 1.8767E-002, 7.0042E-003,
     &                        3.0618E-003,-1.3121E-003,-6.8901E-002,
     &                        2.3944E-002, 4.5513E-003,-5.0356E-003,
     &                       -1.6834E-004, 3.4679E-003, 5.1671E-003,
     &                       -5.0158E-004, 8.3745E-003, 1.1080E-002,
     &                        2.6952E-003, 1.5059E-002, 4.1633E-003,
     &                        1.9192E-003,-2.1979E-003, 1.0921E-002,
     &                       -2.1533E-003, 6.5616E-003, 3.8826E-004,
     &                        1.1183E-004/
C     1500km equinox
      DATA (DHH(3,1,J),J=1,49)/ 9.5205E+000,-2.0222E-006,-1.1864E+000,
     &                        7.5985E-007,-2.6909E-002,-4.4005E-007,
     &                        3.0973E-001, 1.3675E-001,-4.8308E-007,
     &                       -2.2609E-002,-3.7822E-007,-3.5004E-002,
     &                        4.5696E-008, 6.1285E-002, 2.5084E-007,
     &                       -1.2329E-002,-4.2144E-008,-8.3540E-003,
     &                       -2.0386E-007, 4.9594E-002,-2.5843E-007,
     &                       -1.3798E-002,-1.6222E-007,-1.1448E-002,
     &                       -4.3079E-002,-2.8827E-007, 7.7978E-003,
     &                       -2.7406E-008, 1.7355E-003,-3.0989E-002,
     &                       -1.1889E-007,-1.2109E-002,-1.3616E-008,
     &                        1.4991E-002, 3.3699E-008,-7.4128E-005,
     &                       -5.9458E-009, 7.1308E-003, 6.3098E-008,
     &                        1.0935E-003,-2.3455E-002,-1.2795E-007,
     &                        1.5102E-003,-1.3092E-003, 3.8446E-008,
     &                       -3.8817E-004, 5.9744E-008,-3.7063E-003,
     &                       -8.6531E-003/
C     1500km June solstice
      DATA (DHH(3,2,J),J=1,49)/ 9.4197E+000,-1.4545E-002,-1.3593E+000,
     &                        3.3999E-001,-1.9479E-001, 4.7850E-002,
     &                        3.7597E-001, 1.5466E-001,-5.5988E-002,
     &                        5.9355E-002,-4.2765E-002,-2.5813E-002,
     &                        3.1549E-002, 3.0235E-002,-3.1327E-002,
     &                        7.0621E-002,-3.5920E-002, 8.6809E-003,
     &                        2.1371E-002,-1.5567E-002,-2.0810E-002,
     &                        3.6965E-002,-1.1192E-002, 4.8831E-003,
     &                       -1.1065E-003, 1.4821E-002, 1.0736E-003,
     &                        3.4666E-003,-2.4916E-003,-1.3773E-002,
     &                        7.9809E-003,-3.6052E-003, 4.3504E-003,
     &                       -3.6110E-002, 8.2627E-003,-4.7435E-003,
     &                        7.3676E-003, 2.1851E-002, 1.0228E-002,
     &                       -2.6510E-003,-2.3925E-003, 1.1321E-005,
     &                        2.1796E-003, 9.6299E-003, 5.5079E-003,
     &                        5.9111E-003, 4.1776E-003, 2.1261E-003,
     &                        2.4663E-003/
C     2250km equinox
      DATA (DHH(4,1,J),J=1,49)/ 9.3714E+000,-3.4200E-006,-1.6846E+000,
     &                       -8.1017E-008,-1.7081E-002, 3.8242E-006,
     &                        4.5205E-001, 4.3324E-002,-1.2896E-007,
     &                        1.0528E-003,-1.8358E-007,-4.1534E-002,
     &                       -1.3077E-007, 1.5865E-001,-5.4277E-007,
     &                        3.9447E-002,-3.1248E-007,-6.7779E-003,
     &                        2.9556E-007,-7.9049E-003, 1.2398E-007,
     &                       -5.0936E-003, 6.8695E-008,-4.9243E-003,
     &                       -3.9880E-002,-1.0480E-007,-2.4017E-002,
     &                       -5.0755E-008,-3.1475E-003, 1.3985E-002,
     &                       -1.2591E-008,-1.1025E-003,-9.3741E-009,
     &                       -2.9468E-002,-1.2756E-008,-7.7319E-003,
     &                       -1.9948E-008, 1.7882E-003,-3.5214E-009,
     &                        9.4125E-004,-4.4217E-003,-6.8571E-009,
     &                        1.5396E-003,-3.7219E-004, 2.5287E-009,
     &                        6.3584E-004, 8.8245E-009,-1.0107E-003,
     &                       -5.4499E-003/
C     2250km June solstice
      DATA (DHH(4,2,J),J=1,49)/ 9.2865E+000,-3.1689E-002,-1.6086E+000,
     &                        1.1140E-001,-1.7424E-001, 1.2264E-001,
     &                        3.9063E-001,-4.0969E-002, 6.4995E-002,
     &                       -3.0890E-002, 2.6063E-002,-4.6020E-002,
     &                        1.0435E-002, 4.0526E-002, 2.4349E-003,
     &                       -7.3789E-003, 3.2772E-002,-5.0331E-003,
     &                       -1.7989E-002, 1.3328E-002,-1.0732E-002,
     &                       -1.5973E-002, 5.7027E-003,-5.3769E-003,
     &                       -6.0208E-002, 3.3963E-002,-1.2802E-002,
     &                        6.6097E-003,-2.5335E-003,-4.0039E-002,
     &                       -5.7030E-003,-8.9255E-003,-1.2960E-003,
     &                       -9.5029E-003, 2.6328E-002,-5.4759E-003,
     &                        1.8651E-003,-8.9029E-003, 7.3539E-003,
     &                       -1.2681E-003, 2.2857E-003,-4.7377E-003,
     &                       -8.9554E-004, 1.3359E-003, 1.3630E-003,
     &                        3.7962E-003, 2.9694E-003,-4.0318E-003,
     &                        4.5768E-003/
C//////////////////////////////////////////////////////////////////////
C//////////////////////////////////He+/////////////////////////////////
C     550km equinox
      DATA (DHEH(1,1,J),J=1,49)/ 8.5562E+000,-6.7753E-008,-5.9333E-001,
     &                        1.0384E-006,-2.8833E-001,-8.8932E-007,
     &                        1.3836E-001, 2.0484E-001, 1.3434E-008,
     &                       -5.7765E-002,-4.4369E-008, 2.7521E-002,
     &                       -8.8706E-008, 3.9828E-001, 9.9476E-008,
     &                       -1.9442E-002, 3.6408E-007,-2.6872E-002,
     &                       -3.1503E-007,-1.9197E-001,-4.7243E-008,
     &                        8.8640E-004,-1.3824E-008,-1.3470E-002,
     &                       -2.7662E-001, 1.2720E-007,-9.6029E-003,
     &                        4.5527E-008, 1.0784E-003,-9.8919E-002,
     &                       -3.6199E-009,-7.8780E-003, 1.1979E-008,
     &                       -6.3545E-002, 5.1229E-009, 2.1340E-002,
     &                        6.9696E-008,-2.0368E-003, 1.9672E-008,
     &                        1.3091E-003, 1.4277E-002, 8.1100E-008,
     &                        1.5916E-003,-6.1413E-003, 4.4252E-009,
     &                        5.2119E-003, 4.2404E-008,-1.5831E-002,
     &                        8.7771E-003/
C     550km June solstice
      DATA (DHEH(1,2,J),J=1,49)/ 8.4548E+000,-4.6108E-001,-2.8570E-001,
     &                        3.1867E-001,-1.3751E-003,-2.2585E-001,
     &                        2.7861E-002, 5.3882E-002, 2.4531E-001,
     &                       -1.6614E-002,-2.1614E-002, 5.2188E-004,
     &                       -7.1877E-004, 1.0749E-001, 7.9494E-003,
     &                       -5.7983E-002,-3.5115E-002, 2.1946E-003,
     &                        4.4208E-003,-1.2911E-001, 6.0450E-002,
     &                        1.7482E-002, 5.5533E-003,-4.0699E-003,
     &                       -2.3475E-001, 6.2598E-002,-2.5732E-002,
     &                       -6.0082E-003, 4.0616E-004,-3.5725E-002,
     &                       -7.0466E-003, 3.2939E-003,-4.0902E-003,
     &                       -1.5459E-002, 2.4557E-003,-1.6614E-003,
     &                       -8.6784E-004, 2.9412E-002,-7.3953E-004,
     &                       -3.2179E-003,-4.6908E-003,-1.4583E-003,
     &                        9.0449E-004,-7.5250E-003,-1.0508E-003,
     &                        1.8460E-003,-1.3508E-003,-6.8072E-003,
     &                       -1.3448E-003/
C     900km equinox
      DATA (DHEH(2,1,J),J=1,49)/ 9.2145E+000, 1.6979E-007,-9.1857E-001,
     &                        2.0716E-006,-2.3914E-001,-3.5727E-006,
     &                        3.1041E-001,-4.7453E-002, 1.0117E-006,
     &                        6.5570E-002, 3.0985E-007,-6.6558E-002,
     &                       -4.5164E-007, 7.7242E-002, 2.8717E-007,
     &                       -8.9984E-002, 5.2001E-008,-4.5497E-002,
     &                       -1.4626E-007,-1.0729E-001,-3.8360E-008,
     &                        1.6964E-002, 8.6957E-008, 1.8799E-004,
     &                       -1.8424E-002,-1.4543E-007, 3.3958E-002,
     &                        1.0099E-007,-2.6298E-003,-1.6571E-001,
     &                        6.3739E-008,-1.7494E-003, 1.6987E-007,
     &                        5.1204E-002, 2.6500E-007, 5.0013E-003,
     &                        3.6144E-008, 8.7354E-003,-3.7456E-008,
     &                        2.8601E-003,-1.7103E-003,-7.1480E-009,
     &                       -2.5359E-003,-3.0116E-003, 5.7501E-008,
     &                        5.6945E-003, 6.9904E-008,-1.9349E-002,
     &                       -1.2641E-002/
C     900km June solstice
      DATA (DHEH(2,2,J),J=1,49)/ 8.7717E+000,-4.8916E-001,-5.3663E-001,
     &                        2.0667E-001,-1.1371E-001, 2.6697E-003,
     &                        1.4339E-002,-2.7105E-001, 3.6089E-002,
     &                        1.6995E-002, 3.0412E-002,-5.3812E-002,
     &                       -5.6648E-003,-4.5882E-002,-1.2868E-002,
     &                       -1.8613E-002,-9.6389E-003,-5.0658E-003,
     &                       -5.5670E-003, 6.7055E-003,-1.0661E-002,
     &                       -8.0307E-003,-4.6221E-003, 8.3254E-003,
     &                       -1.2762E-001, 2.5276E-002,-2.4689E-002,
     &                        1.2249E-002,-1.0867E-002, 2.5041E-002,
     &                        1.8546E-002, 1.0281E-002,-2.5269E-003,
     &                       -5.5967E-002, 1.4229E-002, 1.0142E-002,
     &                        2.1503E-003, 3.5215E-002, 5.8207E-003,
     &                       -4.2017E-004, 2.7797E-003,-1.4003E-003,
     &                        1.5366E-003, 1.4103E-002, 6.1949E-003,
     &                        1.0431E-002, 1.7503E-003, 1.6873E-002,
     &                        1.0208E-002/
C     1500km equinox
      DATA (DHEH(3,1,J),J=1,49)/ 9.1382E+000,-1.6957E-006,-1.0442E+000,
     &                        1.2347E-006,-3.4358E-001, 6.8230E-007,
     &                        4.4436E-001, 1.6324E-002,-4.3446E-007,
     &                       -2.3140E-002,-1.8743E-007,-3.8261E-002,
     &                        1.6702E-007,-3.1722E-002,-5.1618E-007,
     &                        4.6237E-002,-1.1968E-007,-1.7095E-002,
     &                        2.3684E-007,-1.0912E-001, 8.5906E-008,
     &                       -7.5780E-003,-5.7318E-008, 3.4197E-003,
     &                       -1.6487E-001,-5.7436E-007, 3.4923E-002,
     &                       -4.8531E-010,-1.2999E-002,-5.7078E-002,
     &                        3.2258E-008,-4.1749E-003, 4.3545E-008,
     &                       -4.3926E-002, 8.9899E-009, 9.0836E-003,
     &                       -1.0023E-008, 2.6406E-002,-4.3460E-008,
     &                        9.7611E-004,-4.5719E-002, 9.5360E-008,
     &                        3.2015E-003,-3.0497E-003,-4.9917E-008,
     &                       -5.2358E-003, 1.5254E-008,-2.0919E-002,
     &                       -4.6316E-003/
C     1500km June solstice
      DATA (DHEH(3,2,J),J=1,49)/ 8.6956E+000,-3.7513E-001,-7.6166E-001,
     &                        5.8258E-002,-4.4109E-001, 1.1952E-001,
     &                        4.4731E-001, 6.1168E-002,-5.0179E-002,
     &                        5.2661E-003,-1.0086E-002,-1.1434E-002,
     &                        2.4774E-002, 3.4482E-002,-1.3156E-001,
     &                        3.0048E-002,-6.1629E-002,-8.0153E-003,
     &                       -3.4766E-003, 1.7056E-002,-2.5625E-002,
     &                        1.6715E-002, 6.7878E-003,-2.1396E-003,
     &                       -4.9347E-002, 1.5841E-002, 1.7152E-003,
     &                        7.2410E-004, 1.6179E-003,-3.0269E-002,
     &                        1.6246E-002, 1.4672E-003, 2.1578E-003,
     &                       -2.5516E-002, 1.2958E-002,-6.5452E-003,
     &                        3.3270E-003,-1.6907E-003, 7.4365E-003,
     &                       -1.1628E-003, 1.0449E-002,-1.2092E-003,
     &                        2.5582E-003, 7.1529E-003, 1.1959E-003,
     &                        9.8574E-003, 1.6045E-003, 3.2251E-003,
     &                        7.2336E-003/
C     2250km equinox
      DATA (DHEH(4,1,J),J=1,49)/ 8.8817E+000,-2.4391E-006,-1.3655E+000,
     &                        8.4386E-007,-2.2789E-001, 8.7854E-007,
     &                        4.3537E-001, 1.8615E-002,-7.1575E-007,
     &                        1.3645E-002,-4.9988E-007,-3.2013E-002,
     &                        3.6249E-007,-1.4024E-002,-4.6747E-007,
     &                        4.3299E-002,-2.0939E-007,-1.9031E-002,
     &                        5.3640E-007,-5.7979E-002,-1.1916E-007,
     &                       -9.7389E-003,-7.0552E-008,-3.8403E-003,
     &                       -1.8747E-002,-4.5508E-007, 5.8662E-004,
     &                       -1.4839E-007,-2.3248E-003, 3.5164E-002,
     &                       -3.8530E-008, 1.4493E-003, 9.0911E-009,
     &                        6.5064E-002,-7.7367E-008,-7.6839E-003,
     &                       -4.4329E-008, 2.6307E-002, 6.0727E-008,
     &                        2.6258E-003,-8.0157E-003,-8.9614E-009,
     &                        5.7615E-004,-1.3741E-003, 1.9524E-008,
     &                        3.2822E-003,-5.1476E-009,-1.3574E-003,
     &                        1.7923E-003/
C     2250km June solstice
      DATA (DHEH(4,2,J),J=1,49)/ 8.5683E+000,-1.4807E-001,-8.1147E-001,
     &                       -1.8822E-002,-4.5791E-001, 6.8315E-002,
     &                        4.1478E-001,-3.7850E-002, 1.2571E-001,
     &                        1.2702E-002, 7.4260E-003,-3.1216E-002,
     &                        1.3414E-002,-9.9827E-002, 4.7220E-002,
     &                       -2.3880E-003,-6.6048E-003,-6.6849E-003,
     &                       -1.4619E-002,-1.1156E-002, 3.6579E-003,
     &                        4.1293E-003, 1.4703E-002, 3.9578E-003,
     &                       -5.7658E-002, 3.0119E-003,-1.0730E-002,
     &                        3.4455E-003,-5.9462E-003,-1.8327E-002,
     &                       -1.0196E-002,-5.1724E-003, 2.0240E-003,
     &                       -2.8124E-002, 1.0821E-002,-4.7065E-003,
     &                        3.4806E-003,-1.4711E-002, 4.2184E-003,
     &                       -2.0076E-003, 9.1568E-003,-2.7745E-003,
     &                        9.3070E-004,-1.1010E-003,-2.9865E-004,
     &                       -1.0161E-003, 1.0107E-003, 7.9248E-004,
     &                        4.3800E-003/      
C//////////////////////////////////////////////////////////////////////
C///////////////////////////////////N+/////////////////////////////////
C     550km equinox
      DATA (DNH(1,1,J),J=1,49)/ 1.0019E+001, 2.5723E-006,-8.4960E-004,
     &                        8.1102E-007, 9.7779E-002, 1.0615E-006,
     &                        7.6222E-002,-1.4428E-001, 3.2342E-007,
     &                       -1.4089E-002,-3.4285E-007,-1.3037E-002,
     &                        7.6984E-008,-2.0229E-001, 8.1544E-007,
     &                       -4.2526E-002, 1.9650E-007, 8.5229E-003,
     &                        2.0090E-007, 4.4525E-002,-7.3598E-008,
     &                        1.1851E-002,-7.3302E-008, 1.0225E-003,
     &                        2.3862E-002, 2.6317E-008, 1.7643E-002,
     &                       -6.8800E-008,-1.8743E-003,-2.5839E-002,
     &                        2.0570E-008, 3.7364E-003,-7.4188E-008,
     &                       -2.1742E-002,-9.5189E-008, 5.7009E-003,
     &                       -2.2182E-008, 2.5205E-003, 6.1330E-008,
     &                        1.2046E-003,-1.4426E-002, 2.2803E-008,
     &                        1.8722E-003,-2.2612E-003,-8.6377E-008,
     &                       -2.3099E-003, 6.1868E-008,-2.2280E-003,
     &                       -7.2127E-003/
C     550km June solstice
      DATA (DNH(1,2,J),J=1,49)/ 9.8950E+000, 5.3898E-001,-1.8438E-001,
     &                       -3.3234E-002, 6.1092E-002,-4.4223E-002,
     &                        4.9314E-002,-1.6885E-001, 2.3645E-002,
     &                        4.9904E-002,-1.3272E-002,-1.0006E-002,
     &                        9.8171E-003,-1.2770E-001, 5.4348E-002,
     &                       -9.4108E-003,-1.6097E-002, 2.5746E-003,
     &                       -5.6109E-003,-9.5980E-002, 2.4506E-002,
     &                       -1.9320E-002, 4.4041E-003, 5.8621E-004,
     &                        1.0298E-002, 6.6383E-003,-3.6596E-003,
     &                        3.1693E-003,-1.2114E-003, 2.7192E-002,
     &                       -1.2499E-002, 2.2087E-003, 7.1793E-005,
     &                        2.9456E-002,-1.8417E-002, 6.6479E-003,
     &                       -1.0839E-003,-3.2662E-002, 2.4478E-003,
     &                       -1.3406E-003,-2.3848E-002, 3.9534E-003,
     &                       -1.9607E-003,-2.4437E-002,-2.5257E-004,
     &                       -1.2976E-002, 5.7891E-004,-1.3725E-002,
     &                        3.9792E-004/
C     900km equinox
      DATA (DNH(2,1,J),J=1,49)/ 9.4694E+000,-1.6544E-006, 1.2337E-001,
     &                        7.6922E-006, 4.0844E-002,-6.2726E-006,
     &                       -1.2085E-002,-3.1155E-001,-2.2502E-006,
     &                        5.2141E-003, 1.5292E-006, 4.6169E-003,
     &                       -3.4430E-007,-2.0980E-001,-6.0445E-007,
     &                        4.9811E-002, 5.3263E-007, 6.5131E-004,
     &                       -1.6033E-007,-1.0316E-001,-1.1576E-006,
     &                       -1.6319E-002, 1.8240E-007, 1.0867E-002,
     &                        7.1700E-002,-1.0400E-006, 1.5320E-003,
     &                        7.3160E-007, 5.7948E-003,-7.4410E-002,
     &                       -5.4983E-007,-1.1335E-002,-9.1272E-009,
     &                       -4.2801E-002,-6.1656E-007,-1.0028E-002,
     &                        1.2323E-007,-2.0959E-003,-2.8620E-007,
     &                       -4.1891E-003,-1.0114E-002,-2.4946E-007,
     &                       -5.4484E-003,-1.0014E-004,-6.1357E-008,
     &                        2.5335E-003,-1.4098E-007,-2.7616E-002,
     &                        1.7588E-002/
C     900km June solstice
      DATA (DNH(2,2,J),J=1,49)/ 9.2850E+000, 5.5985E-001, 2.2704E-001,
     &                       -7.0123E-002,-1.1388E-001,-9.4701E-002,
     &                        2.8391E-002,-1.9589E-001, 6.8225E-002,
     &                       -5.3756E-002, 1.4622E-002, 5.7994E-004,
     &                       -2.9261E-003,-2.4712E-001,-1.7103E-002,
     &                        6.1000E-002, 1.5448E-002,-1.2181E-002,
     &                       -8.6816E-003,-2.4797E-002, 6.9824E-002,
     &                       -4.3173E-004,-2.9083E-003,-1.6011E-003,
     &                       -2.1840E-002,-3.5010E-002, 1.4547E-002,
     &                       -9.2115E-003, 4.7761E-003, 1.2071E-002,
     &                        7.7499E-003, 8.1855E-003, 2.2545E-003,
     &                        2.4944E-002,-6.0505E-003, 7.8059E-003,
     &                        2.3805E-003,-5.1867E-002,-8.0380E-003,
     &                       -9.7830E-004,-2.9220E-002,-1.9427E-003,
     &                        3.0433E-003,-3.3815E-002,-2.2415E-003,
     &                       -2.9646E-002,-3.4173E-003,-2.9588E-002,
     &                       -9.4138E-003/
C     1500km equinox
      DATA (DNH(3,1,J),J=1,49)/ 8.7155E+000, 1.4741E-006, 3.3257E-001,
     &                       -3.9187E-007, 9.1878E-003,-1.6825E-006,
     &                       -4.4178E-002,-4.9152E-001, 1.8469E-006,
     &                       -1.0781E-002, 4.6555E-007, 3.2247E-002,
     &                       -4.2322E-007,-2.7076E-001, 5.7315E-007,
     &                        1.0441E-001, 1.5751E-007,-1.9427E-002,
     &                       -2.5181E-007,-1.5498E-001, 4.8205E-007,
     &                       -2.9731E-002, 8.4155E-008, 2.7337E-003,
     &                       -1.5423E-001, 4.0527E-007, 1.0853E-002,
     &                        1.7637E-007, 2.9832E-004,-3.2027E-002,
     &                        4.2496E-008,-1.0391E-002, 1.6182E-008,
     &                        2.0184E-002, 1.1780E-007,-3.6163E-003,
     &                        7.9778E-008,-3.7631E-002, 3.5515E-008,
     &                       -3.8295E-003, 2.5399E-003,-4.9073E-009,
     &                       -2.8015E-003, 2.1214E-004, 2.7217E-008,
     &                       -1.1954E-002,-9.2222E-008,-9.6978E-003,
     &                        6.5686E-003/
C     1500km June solstice
      DATA (DNH(3,2,J),J=1,49)/ 8.6689E+000, 7.7532E-001, 1.6981E-001,
     &                       -1.6926E-001,-1.6344E-002,-4.9955E-002,
     &                       -3.0453E-002,-2.4772E-001, 1.2450E-001,
     &                        1.3351E-002, 8.7996E-004, 1.8177E-002,
     &                       -2.5308E-002,-1.0649E-001, 7.8115E-002,
     &                        7.8976E-003,-3.8263E-002, 2.3854E-002,
     &                       -3.1419E-002,-1.3547E-003,-7.1151E-002,
     &                        2.6453E-002,-6.8153E-004, 1.0626E-003,
     &                        2.1740E-002, 3.1895E-002, 1.2209E-003,
     &                        6.4036E-003, 2.8330E-003, 5.4547E-003,
     &                       -2.6100E-002, 6.7990E-003,-1.1911E-003,
     &                       -7.8806E-003,-1.5986E-002,-6.6544E-003,
     &                        2.2161E-003,-1.6487E-002,-8.9076E-003,
     &                        4.5130E-004,-1.1831E-002,-1.1070E-002,
     &                       -2.2416E-003,-1.4739E-002,-3.8693E-003,
     &                       -1.1053E-002,-1.5497E-003,-1.3144E-002,
     &                       -6.7004E-003/
C     2250km equinox
      DATA (DNH(4,1,J),J=1,49)/ 8.2446E+000,-3.9022E-006, 4.2086E-001,
     &                       -8.4724E-008,-4.5965E-002, 5.3823E-007,
     &                       -1.3220E-002,-2.2697E-001,-5.8269E-007,
     &                       -1.4812E-002,-2.9880E-007,-3.9232E-003,
     &                        1.7290E-008,-4.6253E-001,-9.6196E-008,
     &                        6.2313E-002, 4.2635E-007, 1.7812E-002,
     &                        1.2969E-007,-6.4979E-002, 1.4935E-007,
     &                       -3.5472E-002, 6.1687E-008,-7.7464E-003,
     &                       -4.4305E-002,-3.3255E-007,-2.1269E-003,
     &                       -2.9698E-007, 4.2903E-003, 1.0216E-002,
     &                        2.0950E-008, 3.3693E-003, 1.5003E-008,
     &                        1.4030E-001, 9.8769E-008,-4.3899E-003,
     &                       -6.7915E-008, 6.9364E-002, 4.3913E-008,
     &                        3.6741E-003, 2.6018E-002, 8.9366E-008,
     &                        3.9212E-004, 8.6551E-003, 6.6664E-008,
     &                        9.1264E-003, 1.4846E-008,-1.1903E-003,
     &                        6.8326E-005/
C     2250km June solstice
      DATA (DNH(4,2,J),J=1,49)/ 8.1232E+000, 6.3400E-001, 3.4451E-001,
     &                       -1.7706E-001,-9.6395E-002,-3.3631E-002,
     &                       -3.8348E-002,-3.7099E-001, 4.3368E-002,
     &                       -2.2563E-002, 4.1095E-002,-2.2098E-002,
     &                        3.7783E-003,-3.8367E-001, 1.1212E-001,
     &                        3.9900E-002,-3.4603E-002, 9.0102E-004,
     &                        1.7867E-004,-4.9073E-002, 5.7379E-003,
     &                       -1.6769E-002, 9.4616E-003,-4.8618E-003,
     &                       -2.0081E-001, 7.2646E-003,-9.1609E-003,
     &                        8.1090E-003,-5.8271E-004, 1.3650E-001,
     &                       -2.4914E-003,-1.2280E-003,-4.3112E-004,
     &                       -9.0520E-002, 2.2466E-002,-7.9454E-003,
     &                        2.7647E-003, 8.3267E-003,-6.3168E-003,
     &                        2.2216E-003,-8.6691E-003, 2.5916E-003,
     &                       -1.9514E-003,-4.9024E-003,-3.6766E-003,
     &                       -2.3492E-002,-7.9399E-004, 3.7814E-003,
     &                       -1.6833E-002/
C//////////////////////////////////////////////////////////////////////
C/////////////////coefficients low solar activity///////////////////// 
C//////////////////////////////////O+//////////////////////////////////
C     400km equinox
      DATA (DOL(1,1,J),J=1,49)/ 1.1028E+001,-5.9947E-007,-3.3742E-001,
     &                        7.8683E-007, 1.3937E-001, 8.3175E-008,
     &                       -4.8640E-002,-2.4412E-001, 3.1765E-007,
     &                       -3.0586E-002,-3.2124E-008,-2.0091E-002,
     &                        9.5865E-008,-2.4533E-001,-3.1051E-007,
     &                       -1.3402E-002, 1.5715E-007,-3.5108E-002,
     &                       -7.0552E-009, 7.2451E-002,-5.4795E-008,
     &                        4.6388E-003, 2.6538E-008,-1.7863E-003,
     &                       -3.2795E-002, 2.1667E-007,-1.6886E-003,
     &                       -1.1351E-007, 4.8819E-004,-1.8702E-002,
     &                        1.5855E-008, 2.4774E-003,-6.0513E-011,
     &                        2.6071E-002,-1.5258E-008, 5.4362E-003,
     &                        2.6064E-008,-5.8978E-003,-5.6781E-009,
     &                        1.7888E-003, 7.9854E-003, 4.3021E-008,
     &                        7.5171E-004,-3.4351E-003, 2.3671E-008,
     &                        1.8393E-003, 4.3624E-011,-3.1620E-003,
     &                        2.2396E-003/
C     400km June solstice
      DATA (DOL(1,2,J),J=1,49)/ 1.0935E+001, 4.3766E-001,-3.1289E-001,
     &                       -1.3726E-001, 1.3358E-001, 4.3905E-002,
     &                        5.0160E-002,-2.5849E-001, 4.5659E-002,
     &                        5.4857E-003, 2.2218E-002,-2.8317E-002,
     &                        8.2630E-003,-2.8326E-001,-9.6791E-003,
     &                        4.1894E-002, 7.7419E-003,-1.9546E-002,
     &                       -1.0279E-002, 3.9292E-002,-1.5802E-002,
     &                        9.4275E-003,-8.1618E-004, 1.8252E-003,
     &                       -4.4936E-002,-9.5925E-003, 1.0290E-002,
     &                       -4.3728E-003, 2.6103E-003, 1.7584E-002,
     &                       -1.5570E-004,-5.0717E-003, 2.2696E-004,
     &                        7.7040E-003, 3.7456E-003,-1.7878E-003,
     &                        8.0229E-004,-1.4827E-002, 6.1296E-003,
     &                       -1.2200E-003, 1.6822E-003, 6.2744E-003,
     &                       -2.1599E-003,-2.2537E-002, 3.3963E-003,
     &                       -8.4230E-003, 2.8959E-003,-1.7315E-002,
     &                       -4.7823E-003/
C     550km equinox
      DATA (DOL(2,1,J),J=1,49)/ 1.0441E+001,-1.8583E-007,-1.5392E-001,
     &                        1.3279E-006,-2.8910E-002, 4.6595E-008,
     &                       -4.0713E-002,-4.8280E-001,-4.4062E-007,
     &                        8.8233E-003, 1.5207E-007,-2.6193E-002,
     &                       -5.6026E-008,-2.5563E-001, 3.7636E-007,
     &                       -5.3640E-002,-2.5118E-007,-3.9499E-002,
     &                        5.6905E-008, 6.3005E-002,-1.0823E-007,
     &                       -4.6264E-003, 1.3677E-008, 5.7388E-003,
     &                       -2.6456E-002, 2.1980E-007,-6.3113E-003,
     &                       -9.6352E-008, 6.3825E-003, 2.1161E-002,
     &                       -1.1267E-007,-3.3647E-003, 1.6196E-008,
     &                        2.6800E-002, 1.7619E-007, 9.5012E-003,
     &                       -4.3477E-009, 1.9645E-002, 6.8076E-009,
     &                        8.3642E-004, 2.0048E-002,-2.6173E-008,
     &                        2.8416E-003, 1.2859E-003, 2.4308E-008,
     &                       -3.1225E-003,-7.8282E-009,-4.2809E-003,
     &                        6.4016E-003/
C     550km June solstice
      DATA (DOL(2,2,J),J=1,49)/ 1.0428E+001, 5.0284E-001,-1.0415E-001,
     &                       -2.3068E-001, 9.9921E-002, 2.4525E-002,
     &                        1.0017E-001,-4.3949E-001, 7.7712E-002,
     &                        2.7670E-002, 9.8040E-003,-4.0890E-002,
     &                        1.1955E-002,-2.3115E-001,-1.7800E-002,
     &                        2.9761E-002, 3.2814E-003,-1.0608E-002,
     &                       -1.1838E-002, 4.6797E-002,-2.5308E-002,
     &                       -2.1772E-003,-6.2791E-004, 3.9284E-004,
     &                       -9.5448E-002,-7.5227E-003, 4.7040E-003,
     &                        5.5001E-003, 1.1420E-003, 5.9417E-002,
     &                       -1.0356E-002,-1.1835E-002, 4.9372E-003,
     &                       -1.6572E-002, 7.7212E-003, 2.5303E-003,
     &                       -2.1244E-003, 1.8799E-003, 7.5457E-003,
     &                       -3.9116E-003, 2.9587E-003, 1.0549E-003,
     &                       -1.6082E-004,-2.2565E-002,-3.8407E-004,
     &                       -6.4016E-003, 1.9367E-003,-2.4715E-002,
     &                       -4.3424E-003/
C     750km equinox
      DATA (DOL(3,1,J),J=1,49)/ 9.7046E+000,-1.9504E-006, 1.6453E-001,
     &                        4.1395E-006,-2.8124E-002,-5.0045E-006,
     &                        1.6867E-001,-9.4930E-001,-1.9466E-006,
     &                        7.2760E-004, 4.5795E-007,-4.0390E-002,
     &                        3.8333E-009,-4.1313E-001, 7.3783E-007,
     &                       -9.0851E-002, 8.4989E-007,-9.4220E-003,
     &                       -8.8451E-007, 8.8378E-002,-2.2731E-007,
     &                       -1.3118E-003,-2.4401E-007, 6.7618E-003,
     &                       -7.5769E-002, 3.1554E-007,-9.4427E-003,
     &                       -3.2472E-008,-8.2261E-004, 1.2120E-001,
     &                        9.9090E-008,-9.4052E-003, 7.1322E-008,
     &                        9.7052E-002, 1.2635E-007, 1.2328E-002,
     &                       -6.2071E-008, 2.3364E-002, 2.8730E-009,
     &                       -4.0942E-003, 4.4158E-002,-8.1883E-008,
     &                        9.6001E-004, 3.4871E-003,-8.4680E-008,
     &                       -2.9931E-003,-8.0600E-008, 3.5474E-004,
     &                        5.9745E-003/
C     750km June solstice
      DATA (DOL(3,2,J),J=1,49)/ 9.7665E+000, 7.0796E-001, 9.7141E-002,
     &                       -3.4388E-001, 1.5406E-001,-9.9645E-002,
     &                        8.6633E-002,-6.6761E-001, 1.9789E-001,
     &                        1.2181E-002, 1.9043E-002,-1.4917E-002,
     &                       -1.3319E-002,-2.4644E-001, 1.9972E-002,
     &                        7.0064E-003, 2.0894E-002,-4.9076E-003,
     &                       -1.8647E-002, 2.7701E-002,-7.5134E-003,
     &                       -1.0716E-003,-3.5230E-003, 3.1439E-003,
     &                       -1.9860E-001, 1.9707E-002, 3.3487E-003,
     &                        2.0649E-002, 2.9560E-004, 9.4027E-002,
     &                       -3.3459E-002, 3.2062E-003, 2.2061E-003,
     &                       -4.4366E-002, 8.8768E-003, 3.6112E-003,
     &                        6.0056E-004,-9.2459E-003,-4.6701E-003,
     &                        1.9506E-003,-1.4880E-003,-9.6133E-003,
     &                        2.2348E-003,-3.2551E-002,-5.1712E-003,
     &                       -2.8312E-002,-1.5611E-004,-3.9494E-002,
     &                       -1.3699E-002/
C     1000km equinox
      DATA (DOL(4,1,J),J=1,49)/ 9.1361E+000,-9.7430E-007, 1.6267E-001,
     &                        4.3277E-006, 1.1540E-001,-3.0602E-006,
     &                        5.8135E-002,-1.1270E+000,-6.2546E-007,
     &                       -4.8728E-002, 1.7487E-006,-2.1179E-002,
     &                       -4.7905E-007,-4.0648E-001, 1.3901E-006,
     &                       -1.2238E-001,-4.2540E-007,-3.9834E-002,
     &                        1.7409E-007,-3.8208E-002,-1.0080E-006,
     &                        3.4987E-003, 5.4186E-007,-8.0568E-004,
     &                       -5.4931E-002, 8.2682E-007, 2.1062E-002,
     &                        2.2819E-007,-8.6679E-003, 7.2199E-002,
     &                        5.2536E-008,-5.0089E-003,-2.0015E-008,
     &                        8.7559E-002, 3.1601E-007, 8.5504E-003,
     &                        1.1819E-007,-3.0228E-003, 8.0826E-008,
     &                       -2.8599E-003, 2.6641E-002, 2.1867E-007,
     &                       -2.9991E-003,-1.8339E-003, 4.6824E-008,
     &                       -1.7178E-003, 1.4137E-007, 2.2711E-003,
     &                       -3.6990E-004/
C     1000km June solstice
      DATA (DOL(4,2,J),J=1,49)/ 9.2203E+000, 6.2242E-001, 3.7901E-001,
     &                       -3.7432E-001, 2.3634E-001,-2.0471E-001,
     &                        3.3947E-002,-7.8196E-001, 7.2937E-002,
     &                       -4.0200E-003, 1.4082E-002, 2.6208E-002,
     &                       -2.7142E-002,-2.9473E-001,-2.0525E-002,
     &                        9.2248E-002, 7.9193E-003, 4.4462E-002,
     &                       -4.4945E-002,-2.9301E-002,-7.1699E-002,
     &                       -5.8144E-003,-4.6876E-003,-5.4680E-003,
     &                       -2.1563E-001, 2.8593E-002, 9.6277E-003,
     &                        1.9712E-002, 1.3128E-002, 7.2789E-002,
     &                       -1.0102E-002, 1.0531E-003,-5.9706E-003,
     &                        1.7574E-002,-3.6999E-003,-2.9097E-004,
     &                        5.3090E-003,-1.2282E-002,-4.8710E-004,
     &                        6.4886E-004, 8.4747E-003, 8.8815E-005,
     &                       -3.6183E-003,-1.4276E-002,-5.3265E-003,
     &                       -2.2889E-002, 3.0205E-003,-1.1674E-002,
     &                       -4.9931E-003/
C//////////////////////////////////////////////////////////////////////
C//////////////////////////////////H+//////////////////////////////////
C     400km equinox
      DATA (DHL(1,1,J),J=1,49)/ 8.6764E+000,-5.8497E-007,-4.5196E-001,
     &                        7.6173E-007, 8.3042E-002, 7.3664E-008,
     &                        1.8834E-003, 9.2107E-003,-1.1006E-007,
     &                       -2.3918E-003, 9.0788E-008,-1.0721E-002,
     &                       -1.2242E-008,-3.9583E-002,-1.9542E-007,
     &                        3.4700E-002, 2.6014E-007,-2.6578E-002,
     &                       -6.0613E-008, 1.1531E-001, 3.2992E-008,
     &                       -9.5224E-003, 1.9769E-008,-7.9370E-004,
     &                       -6.5651E-002,-1.3860E-008, 7.9121E-003,
     &                        2.2614E-008,-2.7439E-003,-1.8201E-002,
     &                        5.4429E-008, 1.4850E-003, 1.3296E-008,
     &                        2.5910E-002,-2.0974E-008, 3.0923E-003,
     &                        9.1638E-009,-1.6210E-002,-1.6375E-008,
     &                        7.0846E-004, 2.0751E-002, 5.2109E-009,
     &                       -7.3979E-004,-4.2000E-003,-1.3623E-008,
     &                        3.9613E-003, 2.3420E-008,-3.6036E-003,
     &                        3.7577E-003/
C     400km June solstice
      DATA (DHL(1,2,J),J=1,49)/ 8.7464E+000, 1.3755E-001,-4.9056E-001,
     &                       -4.8260E-002, 1.6354E-001, 6.8999E-002,
     &                        8.5903E-002, 1.1020E-001,-1.7707E-002,
     &                        1.2783E-002, 9.1457E-003,-1.7928E-002,
     &                        1.1734E-002,-3.8862E-002, 2.3903E-002,
     &                        7.1281E-003, 1.1954E-002,-2.2804E-003,
     &                       -8.6807E-003, 6.1365E-002, 2.8978E-003,
     &                       -9.6905E-004,-4.6418E-004, 3.6991E-003,
     &                       -6.0207E-002, 1.2697E-002,-4.1056E-003,
     &                       -4.5619E-003, 1.2074E-003,-1.6589E-002,
     &                        1.6563E-002,-2.1209E-004, 2.4490E-003,
     &                        2.9887E-002, 3.8303E-003, 1.9922E-003,
     &                        2.0211E-004,-2.0114E-002, 5.2825E-003,
     &                        2.5115E-004, 2.2024E-003, 8.1118E-003,
     &                       -9.7042E-004,-7.7664E-003, 6.0330E-003,
     &                        6.8127E-004, 3.9275E-003,-3.8098E-003,
     &                        5.7959E-003/
C     550km equinox
      DATA (DHL(2,1,J),J=1,49)/ 9.3142E+000,-4.9414E-007,-6.4289E-001,
     &                       -2.0151E-006,-7.8257E-002,-1.3701E-006,
     &                       -1.2501E-001, 6.5449E-002,-4.8812E-007,
     &                        1.4953E-002, 3.4548E-007,-3.3211E-002,
     &                       -5.4121E-007, 1.7444E-003,-3.9167E-007,
     &                        1.3480E-002, 5.4127E-007,-6.9938E-002,
     &                       -1.6812E-007, 1.0559E-001,-2.5102E-008,
     &                       -1.0793E-003,-1.5060E-007, 5.9427E-003,
     &                       -7.9428E-004, 2.9945E-007, 1.1037E-002,
     &                        1.2265E-008, 1.8328E-003,-3.5544E-002,
     &                       -1.6020E-007,-6.0623E-004, 1.8470E-008,
     &                        8.4732E-003,-5.0107E-008, 5.1421E-003,
     &                        1.0662E-008,-2.9206E-002, 4.8502E-008,
     &                       -2.6379E-003, 1.4690E-002, 5.0003E-008,
     &                        1.8659E-004,-3.9958E-003, 2.9606E-008,
     &                        6.9968E-003, 1.1040E-008,-5.9425E-004,
     &                        1.2121E-003/
C     550km June solstice
      DATA (DHL(2,2,J),J=1,49)/ 9.4313E+000,-7.1878E-002,-6.9971E-001,
     &                        1.6547E-001, 1.0161E-001, 1.5780E-001,
     &                        2.9380E-001, 1.4776E-001,-6.2061E-002,
     &                        1.8808E-002,-1.2118E-002,-2.0569E-002,
     &                        3.2844E-002,-3.5689E-003, 1.0216E-002,
     &                       -1.7151E-002,-1.6616E-002,-1.2647E-002,
     &                        1.0205E-002, 6.1252E-002, 1.0258E-002,
     &                       -1.0887E-002, 3.1761E-003,-2.8283E-003,
     &                       -4.9813E-002, 3.4723E-002,-3.1727E-002,
     &                        1.5177E-002,-6.5093E-003,-2.4869E-002,
     &                        3.0012E-003,-1.5630E-003, 6.6333E-004,
     &                        7.3866E-003, 1.8974E-002,-1.6075E-003,
     &                        3.4999E-003,-2.2340E-003, 5.1294E-004,
     &                        2.3417E-003, 1.1325E-002, 4.8485E-003,
     &                       -1.3997E-003, 7.1165E-003, 5.1959E-003,
     &                        7.3597E-003, 4.3606E-003, 1.0201E-002,
     &                        6.5080E-003/
C     750km equinox
      DATA (DHL(3,1,J),J=1,49)/ 9.5154E+000, 4.4387E-007,-8.3777E-001,
     &                       -1.4337E-007,-3.6703E-002, 1.1904E-006,
     &                        5.7850E-002, 8.9301E-002,-1.0837E-007,
     &                        6.9043E-002,-2.5006E-007,-3.6223E-002,
     &                        3.6250E-007,-3.6352E-002,-1.4407E-007,
     &                       -2.8911E-002,-7.3010E-008,-4.7333E-002,
     &                        9.3189E-008, 6.4058E-002,-2.5044E-008,
     &                        1.5224E-002, 5.7261E-008, 1.0738E-002,
     &                        1.4165E-002,-2.4750E-007, 9.8615E-003,
     &                       -2.6232E-007, 1.7438E-003,-4.4838E-003,
     &                        1.1605E-007, 1.5004E-003, 6.3422E-008,
     &                        2.0907E-002,-6.8549E-008, 3.1254E-003,
     &                       -6.0456E-008,-2.6752E-002, 3.1056E-009,
     &                       -3.9871E-003, 2.9555E-003,-2.4656E-009,
     &                        1.4655E-003,-5.2332E-004,-1.1751E-007,
     &                        1.6169E-003, 3.9033E-008,-1.3868E-003,
     &                        7.0440E-004/
C     750km June solstice
      DATA (DHL(3,2,J),J=1,49)/ 9.5243E+000,-5.9807E-002,-8.7275E-001,
     &                        2.2117E-001, 1.5301E-001, 1.0504E-001,
     &                        3.3841E-001, 1.0576E-001,-7.2307E-002,
     &                        9.6274E-002,-2.9255E-002,-1.4589E-002,
     &                        2.0619E-002, 3.9283E-002,-5.8550E-002,
     &                       -9.2836E-003,-4.7394E-002,-1.7612E-002,
     &                        4.0558E-002, 7.2842E-002,-1.3231E-002,
     &                        7.7768E-003,-9.7764E-003, 2.0915E-003,
     &                       -1.8160E-002,-1.9791E-002,-8.6255E-003,
     &                        1.2757E-003, 1.8743E-003,-6.6028E-003,
     &                       -5.7812E-003,-1.1538E-002, 4.6084E-003,
     &                        2.8366E-003,-2.2018E-003, 4.9727E-003,
     &                       -5.6196E-003,-1.0703E-003, 6.8549E-003,
     &                       -5.9359E-004,-1.0457E-003, 2.7387E-003,
     &                       -2.9393E-003, 6.8835E-003, 4.3346E-003,
     &                        8.0149E-003, 1.5508E-003, 8.8104E-003,
     &                        4.8858E-003/
C     1000km equinox
      DATA (DHL(4,1,J),J=1,49)/ 9.5646E+000, 7.0950E-007,-8.9711E-001,
     &                        1.9288E-006,-1.6734E-001, 4.3813E-006,
     &                        1.1898E-001,-2.0680E-003, 7.7446E-007,
     &                        7.0217E-002, 3.5016E-007, 1.7233E-002,
     &                        5.7925E-007,-1.1643E-001, 2.5154E-007,
     &                       -6.3259E-002, 5.6847E-007,-1.1006E-002,
     &                        5.3229E-008, 7.9520E-002, 3.0792E-007,
     &                        2.3250E-002, 1.2250E-007, 5.4205E-003,
     &                        6.0078E-003, 5.5043E-007,-1.1446E-002,
     &                        2.6075E-007,-2.6486E-003,-1.7948E-002,
     &                        5.7737E-008,-7.7932E-003, 7.3237E-009,
     &                        1.3457E-002, 1.4871E-007,-1.8273E-003,
     &                        9.6374E-008,-1.4154E-002, 3.4815E-008,
     &                       -5.0702E-003, 1.3526E-002, 3.8902E-008,
     &                       -1.5015E-003, 2.5405E-003, 8.9053E-008,
     &                        5.4088E-003,-2.6652E-008, 3.6377E-003,
     &                        3.5619E-004/
C     1000km June solstice
      DATA (DHL(4,2,J),J=1,49)/ 9.4220E+000, 1.2025E-001,-7.9687E-001,
     &                        1.4928E-001, 2.1997E-001,-2.5538E-001,
     &                        3.6075E-001, 1.6937E-001,-1.3366E-002,
     &                        1.9409E-001,-4.7118E-002,-1.3312E-002,
     &                        1.9160E-002,-4.4608E-003,-1.3351E-001,
     &                        2.7297E-002,-4.4798E-002,-1.2571E-002,
     &                        3.3471E-002, 6.4905E-002,-2.7860E-002,
     &                        4.0406E-002,-1.1265E-002, 1.3243E-002,
     &                        3.2774E-004,-4.9982E-002, 8.1220E-003,
     &                        3.5041E-003, 7.9319E-004, 9.8909E-003,
     &                       -5.4041E-003, 7.0907E-004, 4.7071E-004,
     &                        4.4426E-002,-5.3094E-003, 1.0420E-003,
     &                       -1.3821E-003,-2.5690E-002, 2.5739E-003,
     &                       -1.2000E-003, 3.9956E-003,-3.3401E-003,
     &                        3.1321E-003,-1.0688E-002, 1.5289E-003,
     &                       -4.8343E-003, 6.6546E-004,-9.5524E-003,
     &                       -2.9806E-003/
C//////////////////////////////////////////////////////////////////////
C//////////////////////////////////He+/////////////////////////////////
C     400km equinox
      DATA (DHEL(1,1,J),J=1,49)/ 8.0470E+000,-1.2106E-006,-1.8765E-001,
     &                       -4.5251E-007,-1.8951E-001,-1.9758E-007,
     &                        1.3990E-001,-7.3152E-002,-7.7520E-007,
     &                       -1.8854E-002, 3.0364E-007,-3.3138E-002,
     &                       -3.4716E-007, 4.1146E-001,-2.3661E-007,
     &                        3.9635E-002,-3.3287E-009, 2.9161E-003,
     &                       -6.2000E-008, 3.1584E-002,-1.7889E-007,
     &                       -1.1569E-002, 2.1482E-008, 3.2708E-003,
     &                        3.1773E-002, 2.0077E-008,-1.2848E-002,
     &                       -8.0711E-008,-1.3686E-002,-1.8028E-002,
     &                       -8.6107E-008,-8.6541E-003, 1.6977E-008,
     &                        2.1756E-002,-1.0723E-007, 5.4701E-003,
     &                       -9.9269E-009, 1.9697E-002, 3.8353E-008,
     &                       -3.1468E-004,-8.2770E-003,-3.4043E-008,
     &                       -2.7124E-003, 2.9360E-003, 8.8432E-009,
     &                        1.9854E-003,-3.0817E-009, 1.6478E-003,
     &                        2.2964E-003/
C     400km June solstice
      DATA (DHEL(1,2,J),J=1,49)/ 7.7686E+000,-7.0994E-001,-2.9210E-002,
     &                        4.2677E-003,-2.8787E-003, 2.4734E-002,
     &                       -4.2253E-002, 1.6763E-003,-2.2820E-002,
     &                        1.9391E-002, 1.7941E-002, 2.3819E-002,
     &                       -3.7628E-003,-9.3643E-002, 1.0094E-001,
     &                       -3.9462E-002,-4.3404E-002,-1.5448E-003,
     &                       -1.2831E-002,-7.4652E-002, 4.3650E-003,
     &                        1.4048E-002, 4.5125E-003,-7.4554E-004,
     &                       -1.7399E-001, 1.3850E-002, 3.1521E-003,
     &                       -1.1665E-002,-6.9860E-004, 2.9414E-003,
     &                        9.6906E-003,-6.5594E-003, 7.0166E-004,
     &                        3.1185E-002,-1.2213E-002, 1.0771E-004,
     &                        2.6551E-003, 3.5949E-002,-3.5913E-004,
     &                       -9.4525E-004, 1.8059E-002,-3.3219E-003,
     &                       -5.0513E-004, 2.7967E-002,-5.3283E-005,
     &                        1.7352E-002,-9.0916E-004, 2.8545E-002,
     &                        1.3304E-002/
C     550km equinox
      DATA (DHEL(2,1,J),J=1,49)/ 8.5471E+000, 1.0036E-006,-2.1407E-001,
     &                        1.1825E-006,-3.0424E-001, 1.5359E-006,
     &                        1.3422E-001,-1.6014E-002, 1.8552E-006,
     &                        6.7452E-002,-1.4486E-006,-3.6290E-002,
     &                        1.3742E-006, 1.4981E-001, 4.8550E-008,
     &                       -2.4666E-003, 1.4397E-007,-2.4247E-002,
     &                        1.2519E-008,-3.3219E-002,-1.2566E-007,
     &                        3.8258E-003, 1.2039E-007, 6.3462E-003,
     &                        3.6086E-002, 3.5782E-007, 1.8118E-002,
     &                       -3.0822E-008,-1.0032E-002,-3.1399E-002,
     &                        3.2563E-007, 7.5683E-004,-8.9661E-008,
     &                        1.0092E-002,-1.6976E-007, 7.8276E-003,
     &                        5.0854E-008, 5.5042E-003,-1.3474E-007,
     &                       -2.6094E-003, 2.5317E-004, 1.2753E-007,
     &                        2.9530E-003, 5.2850E-003,-3.6304E-008,
     &                       -1.7770E-003,-7.7952E-008, 7.2939E-004,
     &                       -5.0116E-003/
C     550km June solstice
      DATA (DHEL(2,2,J),J=1,49)/ 8.1116E+000,-5.9186E-001,-1.3615E-001,
     &                       -1.8696E-001, 9.5077E-002, 9.0620E-002,
     &                        2.3509E-001,-1.4233E-001,-2.2280E-002,
     &                        7.2726E-002,-1.7002E-002, 1.3668E-002,
     &                       -1.6439E-002,-1.8320E-001,-2.9248E-002,
     &                       -8.7976E-002,-5.1588E-002, 2.7987E-002,
     &                       -2.0067E-002,-6.6101E-002,-1.2527E-002,
     &                        1.6175E-002,-5.1555E-003, 7.3867E-004,
     &                       -3.2870E-001,-1.9420E-002,-1.6805E-002,
     &                        2.5031E-003, 8.4785E-003, 1.3727E-002,
     &                       -3.3331E-002, 4.1629E-003, 3.6778E-004,
     &                        1.4714E-002,-4.3042E-003,-1.2047E-003,
     &                        2.5714E-004, 2.6942E-002,-1.2786E-002,
     &                        5.0574E-004, 2.9082E-002,-9.3908E-003,
     &                        5.0755E-003, 1.5606E-002,-5.0868E-003,
     &                        1.8251E-002,-2.6124E-003, 1.5416E-002,
     &                        1.8188E-002/
C     750km equinox
      DATA (DHEL(3,1,J),J=1,49)/ 8.5395E+000,-1.2355E-006,-3.5963E-001,
     &                       -8.3898E-007,-3.0322E-001, 1.7365E-006,
     &                        1.7334E-001, 7.2227E-002,-1.0290E-006,
     &                        1.0917E-001, 3.3145E-008,-2.0338E-002,
     &                       -2.8147E-007, 1.9255E-002,-3.5607E-008,
     &                       -4.7856E-002,-2.9990E-007, 1.7713E-002,
     &                        8.7537E-008,-3.0223E-002,-2.1642E-007,
     &                        4.9244E-003,-1.5008E-007, 9.0134E-003,
     &                       -5.4076E-003,-1.7635E-007, 2.6828E-003,
     &                        1.8278E-009,-4.6448E-003,-2.0650E-004,
     &                       -1.1633E-007,-2.3058E-004,-7.8885E-008,
     &                        5.5645E-003, 2.2420E-008, 3.6199E-003,
     &                        3.1934E-008,-1.4958E-002,-4.0761E-008,
     &                       -2.4953E-003, 1.9343E-002,-2.4736E-008,
     &                        4.8865E-003, 2.6303E-004, 2.6057E-008,
     &                       -7.6145E-003,-3.5056E-008, 4.4221E-003,
     &                       -1.1743E-002/
C     750km June solstice
      DATA (DHEL(3,2,J),J=1,49)/ 8.0106E+000,-4.4548E-001,-1.8174E-001,
     &                       -1.4286E-001, 2.3503E-002, 1.5303E-001,
     &                        3.0849E-001,-1.3076E-001, 6.3870E-002,
     &                        8.5933E-002,-2.7825E-002,-1.9401E-002,
     &                       -1.1111E-002,-1.5264E-001,-1.1578E-001,
     &                       -1.5510E-001,-3.5013E-002, 2.8185E-002,
     &                       -2.6338E-003, 4.7758E-003, 1.1463E-002,
     &                        1.0401E-002,-6.6656E-003,-1.2770E-003,
     &                       -2.7712E-001,-2.2595E-002,-1.7372E-002,
     &                        3.6625E-003, 1.2089E-003, 3.6612E-002,
     &                       -3.5766E-002,-4.7622E-003, 8.0404E-004,
     &                        2.0356E-002,-1.2600E-002,-3.8058E-003,
     &                       -3.0403E-003, 2.1958E-002,-9.7326E-003,
     &                       -4.7118E-003,-1.8127E-002,-6.5368E-003,
     &                        1.8119E-003, 5.6672E-004,-7.0145E-003,
     &                       -2.2251E-005, 3.0327E-004, 5.2279E-003,
     &                        7.1842E-003/
C     1000km equinox
      DATA (DHEL(4,1,J),J=1,49)/ 8.5126E+000, 6.6638E-007,-4.3467E-001,
     &                        2.6349E-006,-2.1941E-001, 1.2981E-006,
     &                        7.8944E-002,-2.4220E-002,-5.7262E-007,
     &                        9.9112E-002,-1.2219E-008, 3.8104E-002,
     &                        6.1122E-007,-1.3751E-001, 2.0071E-007,
     &                       -5.3385E-002, 6.6087E-007, 3.8775E-002,
     &                        3.6086E-007,-6.2806E-002,-4.5958E-008,
     &                       -5.2240E-003, 2.7794E-008,-1.0112E-002,
     &                       -2.0386E-002,-4.1445E-007,-4.7492E-002,
     &                       -3.1373E-008,-1.0906E-002, 3.1090E-002,
     &                       -7.2809E-009,-7.1495E-003, 3.8796E-008,
     &                        5.3399E-002, 1.9761E-007,-2.0636E-003,
     &                       -9.4328E-008,-5.1179E-004,-5.4058E-008,
     &                       -1.2413E-002, 2.4460E-002,-8.3003E-008,
     &                       -8.0597E-003, 1.0507E-002,-3.2739E-008,
     &                        4.9636E-003,-1.8945E-008, 1.3645E-002,
     &                       -7.1743E-003/
C     1000km June solstice
      DATA (DHEL(4,2,J),J=1,49)/ 7.8861E+000,-2.5140E-001,-2.4407E-002,
     &                       -2.0198E-001,-2.5291E-002, 1.1063E-001,
     &                        1.3606E-001,-5.0165E-002, 1.4380E-001,
     &                        1.2596E-001,-2.5094E-002,-4.4737E-003,
     &                       -1.0152E-002,-2.4527E-001,-2.1181E-001,
     &                       -1.1042E-001,-9.5258E-003, 1.1105E-002,
     &                       -1.8623E-002, 1.1621E-001, 1.3438E-002,
     &                        2.1396E-002, 6.3530E-003, 1.0809E-002,
     &                       -1.5715E-001,-4.7383E-002,-1.0913E-003,
     &                       -1.4630E-002,-3.6688E-003, 3.7541E-002,
     &                       -1.2694E-002, 8.7360E-003,-1.1349E-003,
     &                        4.0484E-002,-2.6862E-002,-3.4730E-003,
     &                       -7.2372E-003,-3.3460E-003,-8.1202E-003,
     &                        5.1497E-003,-1.3913E-002,-3.2409E-003,
     &                        2.1047E-003,-8.3477E-003,-8.1966E-003,
     &                        4.2594E-003,-2.9227E-003,-9.2693E-003,
     &                        2.2682E-003/
C/////////////////////////////////////////////////////////////////////
C//////////////////////////////////N+/////////////////////////////////
C     400km equinox
      DATA (DNL(1,1,J),J=1,49)/ 9.2765E+000,-1.8499E-007,-3.0195E-002,
     &                       -2.3630E-007, 1.2337E-002,-5.3134E-007,
     &                        1.0174E-002,-3.0143E-001,-4.4969E-007,
     &                        4.0224E-002, 2.5418E-007,-2.2113E-002,
     &                       -9.3834E-008,-2.6967E-001,-2.2087E-007,
     &                       -3.0670E-002, 1.8923E-007,-1.7787E-002,
     &                       -7.4715E-008, 5.8658E-002, 1.0175E-007,
     &                       -4.9642E-003, 8.0931E-009,-5.8069E-004,
     &                       -9.6466E-003,-1.9936E-007, 8.5438E-005,
     &                        9.6678E-008,-1.4405E-003,-8.1367E-003,
     &                        1.2139E-007,-6.4643E-004,-2.8006E-009,
     &                        2.1228E-002,-8.6142E-008, 6.2702E-003,
     &                        3.6155E-008, 8.2411E-003, 4.7083E-008,
     &                        3.8409E-004, 1.2262E-002,-1.7636E-008,
     &                        1.5587E-003, 1.9816E-003,-3.3066E-008,
     &                       -2.7297E-003,-1.8297E-008, 1.2335E-004,
     &                        1.1364E-003/
C     400km June solstice
      DATA (DNL(1,2,J),J=1,49)/ 9.2392E+000, 4.0188E-001, 6.3930E-002,
     &                       -8.5958E-002,-4.7951E-003,-2.4724E-002,
     &                       -3.2147E-002,-3.5425E-001, 9.1141E-003,
     &                        5.4439E-002, 2.1320E-002,-2.7290E-002,
     &                       -3.2257E-003,-3.2031E-001,-4.7704E-003,
     &                        4.1115E-002,-3.6383E-003,-1.6398E-002,
     &                       -2.4104E-003, 8.2126E-004,-1.1835E-002,
     &                        5.9562E-003, 1.7010E-003, 5.0395E-004,
     &                       -9.3262E-002,-1.1371E-002, 1.5332E-002,
     &                       -2.0002E-003,-2.8968E-004, 3.9653E-002,
     &                       -5.2539E-003,-5.2158E-003,-3.2324E-004,
     &                        4.3222E-003,-5.3072E-003, 6.3232E-004,
     &                        6.9118E-004,-6.0341E-003,-6.2083E-004,
     &                       -1.3747E-003,-4.2513E-005,-4.3487E-003,
     &                       -2.0489E-004,-1.6498E-002,-3.4412E-003,
     &                       -1.2307E-002,-3.0081E-003,-1.4694E-002,
     &                       -9.0169E-003/
C     550km equinox
      DATA (DNL(2,1,J),J=1,49)/ 9.0063E+000, 8.5692E-007, 1.3438E-001,
     &                        7.5368E-007,-1.0519E-001,-5.5276E-006,
     &                        9.9278E-002,-3.0846E-001, 4.0235E-007,
     &                        6.3959E-002, 9.1603E-007,-2.7028E-002,
     &                       -6.3298E-007,-3.8941E-001, 6.3463E-007,
     &                       -8.3377E-002, 5.6084E-007,-2.0251E-002,
     &                       -5.4261E-007, 8.9982E-002,-9.5019E-008,
     &                       -5.2807E-003, 4.6354E-008, 2.9355E-003,
     &                       -4.6802E-003, 9.8816E-008,-2.7926E-003,
     &                        2.2151E-007, 1.3646E-003, 2.7718E-002,
     &                       -3.4461E-009,-5.5476E-003,-2.4872E-008,
     &                        3.4413E-002,-8.8181E-008, 9.3552E-003,
     &                        1.3853E-008, 2.2646E-002, 2.8251E-008,
     &                       -1.5192E-003, 1.3859E-002, 1.6226E-008,
     &                        3.6852E-003, 2.4081E-003, 3.9640E-008,
     &                       -9.0767E-003, 7.0011E-008,-1.5147E-003,
     &                        7.8946E-004/
C     550km June solstice
      DATA (DNL(2,2,J),J=1,49)/ 8.9712E+000, 6.2259E-001, 2.5911E-001,
     &                       -2.9261E-001, 1.1094E-001,-2.2072E-001,
     &                       -6.3167E-002,-3.7597E-001, 1.1875E-001,
     &                        6.3175E-004, 9.1767E-003,-1.2431E-002,
     &                       -1.0901E-002,-2.7977E-001,-1.3503E-002,
     &                       -1.1409E-002, 2.6856E-002, 9.0748E-003,
     &                       -3.0204E-002, 5.3404E-002,-2.3036E-002,
     &                        3.1147E-004, 5.5788E-003,-5.4323E-003,
     &                       -9.7330E-002,-2.3986E-002, 4.3029E-003,
     &                        2.0047E-002,-3.7340E-003, 2.3783E-002,
     &                       -1.5343E-002,-2.8604E-003, 2.7675E-003,
     &                       -1.5551E-002,-7.1283E-003, 5.3254E-004,
     &                        1.5451E-003,-1.5329E-002, 8.1201E-004,
     &                       -2.0467E-003,-1.3033E-002,-5.4349E-003,
     &                        1.3550E-003,-2.7607E-002,-5.4011E-003,
     &                       -1.6276E-002,-1.6858E-003,-3.5652E-002,
     &                       -1.3964E-002/
C     750km equinox
      DATA (DNL(3,1,J),J=1,49)/ 8.5862E+000,-3.2389E-007, 2.5327E-001,
     &                       -5.5733E-006,-1.0952E-001, 1.3457E-006,
     &                        1.5412E-001,-5.1416E-001, 2.0209E-006,
     &                       -1.1427E-002,-9.9065E-007,-3.0871E-002,
     &                       -1.4579E-007,-3.5447E-001,-1.8618E-006,
     &                       -1.2044E-001, 2.7905E-007,-1.8735E-002,
     &                       -6.8531E-007, 8.9785E-002, 1.2276E-006,
     &                       -7.8418E-003, 1.3500E-008, 4.7094E-003,
     &                        9.3388E-002, 1.3278E-008,-8.3880E-003,
     &                        4.8633E-007,-5.9264E-003, 4.2086E-002,
     &                        9.9745E-008,-1.2530E-002,-1.5369E-007,
     &                        1.0372E-001,-1.7825E-008, 1.1288E-002,
     &                        2.3850E-007, 8.1607E-003,-1.1460E-009,
     &                       -2.8998E-003, 1.5298E-002,-1.7277E-008,
     &                        1.5876E-003, 4.5693E-003, 8.3482E-008,
     &                       -5.8085E-003, 1.1993E-007,-4.8027E-004,
     &                       -2.8110E-005/
C     750km June solstice
      DATA (DNL(3,2,J),J=1,49)/ 8.6860E+000, 6.2474E-001, 3.9864E-001,
     &                       -3.5056E-001, 6.7989E-002,-2.0349E-001,
     &                       -4.3468E-002,-2.8985E-001, 6.9687E-002,
     &                        3.8408E-002,-3.9263E-002, 2.4613E-002,
     &                       -1.6681E-002,-2.0418E-001,-4.9327E-002,
     &                       -1.7733E-002, 6.3639E-002,-1.9788E-003,
     &                       -4.0926E-002, 4.2652E-002,-2.0606E-002,
     &                        1.2350E-002,-2.4292E-002, 1.6044E-002,
     &                       -8.9173E-002,-2.0330E-002, 2.0639E-002,
     &                        6.8728E-003, 2.5223E-003, 3.1266E-002,
     &                       -1.4652E-002, 5.5005E-003,-5.0159E-003,
     &                       -2.7642E-004, 2.0299E-003, 3.6913E-003,
     &                       -1.7547E-003,-1.6781E-002,-1.4129E-002,
     &                        7.0872E-003,-1.0006E-002,-2.8588E-003,
     &                       -6.3962E-004,-3.6016E-002,-7.0462E-003,
     &                       -2.0823E-002,-1.2145E-003,-3.6628E-002,
     &                       -1.7476E-002/
C     1000km equinox
      DATA (DNL(4,1,J),J=1,49)/ 8.2783E+000,-1.3145E-006, 1.6634E-001,
     &                       -4.4507E-007,-9.3753E-003, 2.7889E-006,
     &                       -1.1895E-002,-4.9691E-001,-8.0933E-007,
     &                       -5.0906E-002,-9.0499E-007,-6.5566E-003,
     &                        5.0452E-007,-3.2368E-001, 4.0711E-007,
     &                       -1.1566E-001, 2.2821E-008,-4.6834E-002,
     &                       -1.9870E-008, 8.9269E-002, 5.2219E-007,
     &                        9.2490E-003,-3.0279E-007,-2.4639E-003,
     &                        1.2319E-001, 3.4657E-007, 9.8288E-003,
     &                       -2.0789E-008,-3.0805E-003, 1.2434E-002,
     &                        2.7959E-007,-5.3884E-003,-4.8629E-008,
     &                        8.6375E-002,-1.4080E-007, 1.0559E-002,
     &                        5.5476E-008,-5.8048E-004, 3.2727E-008,
     &                       -9.9461E-004,-5.6637E-003, 1.7444E-007,
     &                       -3.0875E-003, 4.4380E-003,-2.7864E-008,
     &                       -3.0832E-003, 2.9972E-008,-1.9701E-003,
     &                       -7.5436E-003/
C     1000km June solstice
      DATA (DNL(4,2,J),J=1,49)/ 8.3694E+000, 5.9656E-001, 4.4055E-001,
     &                       -2.3417E-001,-1.1354E-001, 1.7281E-002,
     &                       -5.3018E-002,-3.2389E-001, 4.6938E-003,
     &                        8.4262E-003,-8.3879E-003, 7.7294E-003,
     &                        4.6891E-003,-2.2034E-001,-3.4145E-002,
     &                        6.2052E-002, 4.8528E-002,-5.2128E-003,
     &                       -2.4769E-002, 1.4277E-002,-3.0019E-002,
     &                       -1.5376E-002,-6.8305E-003,-1.7002E-003,
     &                       -4.0215E-002, 1.6609E-002, 2.3233E-002,
     &                        5.4203E-003, 1.4267E-003, 5.5382E-002,
     &                       -1.7332E-002,-2.7946E-003,-1.3890E-003,
     &                        3.7206E-002,-4.1274E-003, 4.2794E-003,
     &                        2.7750E-004,-7.0599E-003,-1.0803E-004,
     &                        1.3258E-003,-1.2112E-002, 7.3725E-004,
     &                        6.3315E-004,-1.2200E-002,-5.4789E-003,
     &                       -1.6331E-002,-1.5537E-003,-5.9603E-003,
     &                       -3.9097E-003/
      DATA (CORRH(J),J=1,3)/ 0.762,0.836,1.033/
      DATA (CORRO(J),J=1,3)/ 1.872,1.640,1.234/
C//////////////////////////////////////////////////////////////////////
C/////////////////////////solar minimum////////////////////////////////
      CALL IONLOW(INVDIP,MLT,ALT,DDD,DOL,0,NOL)
      CALL IONLOW(INVDIP,MLT,ALT,DDD,DHL,1,NHL)
      CALL IONLOW(INVDIP,MLT,ALT,DDD,DHEL,2,NHEL)
      CALL IONLOW(INVDIP,MLT,ALT,DDD,DNL,3,NNL)
C///////////////////////////solar maximum//////////////////////////////
      CALL IONHIGH(INVDIP,MLT,ALT,DDD,DOH,0,NOH)
      CALL IONHIGH(INVDIP,MLT,ALT,DDD,DHH,1,NHH)
      CALL IONHIGH(INVDIP,MLT,ALT,DDD,DHEH,2,NHEH)
      CALL IONHIGH(INVDIP,MLT,ALT,DDD,DNH,3,NNH)
C     interpolation (in logarithm)
      IF (PF107 .GT. 260) PF107=260
	IF (PF107 .LT. 65) PF107=65
      NO=(ALOG10(NOH)-ALOG10(NOL))/(210.0-85.0)*(PF107-85.0)+ALOG10(NOL)
      NH=(ALOG10(NHH)-ALOG10(NHL))/(210.0-85.0)*(PF107-85.0)+ALOG10(NHL)
      NHE=(ALOG10(NHEH)-ALOG10(NHEL))/(210.0-85.0)*(PF107-85.0)+
     &     ALOG10(NHEL)
      NN=(ALOG10(NNH)-ALOG10(NNL))/(210.0-85.0)*(PF107-85.0)+ALOG10(NNL)
C     percentages
      NO=10**NO
      NH=10**NH
      NHE=10**NHE
      NN=10**NN
C     normalization
      NTOT=NO+NH+NHE+NN
      NO=NO/NTOT
      NH=NH/NTOT
      NHE=NHE/NTOT
      NN=NN/NTOT
C     C/NOFS correction
      IF (PF107 .LE. 67.5) THEN
       NOCORR=CORRO(1)
       NHCORR=CORRH(1)
      ENDIF
      IF ((PF107 .GT. 67.5) .AND. (PF107 .LT. 87.5)) THEN
       IL=IFIX((PF107-57.5)/10.0)
       NOCORR=(CORRO(IL+1)-CORRO(IL))/10.0*(PF107-57.5-10*IL)+CORRO(IL)
       NHCORR=(CORRH(IL+1)-CORRH(IL))/10.0*(PF107-57.5-10*IL)+CORRH(IL)
      ENDIF 
      IF (PF107 .GE. 87.5) GOTO 10
      NO=NO/NOCORR
      NH=NH/NHCORR
      NHE=NHE/NHCORR
      NN=NN/NOCORR
C     normalization
      NTOT=NO+NH+NHE+NN
      NO=NO/NTOT
      NH=NH/NTOT
      NHE=NHE/NTOT
      NN=NN/NTOT
10    CONTINUE      
      RETURN
      END
C
C
      SUBROUTINE IONLOW(INVDIP,MLT,ALT,DDD,D,ION,NION)
C---------------------------------------------------------------------------
C IONLOW calculates absolute density of O+, H+, He+ or N+  in the outer
C ionosphere for a low solar activity (F107 < 100).
C Based on spherical harmonics approximation of relative ion density
C (by AE-C, and AE-E) at altitudes centred on 400km, 550km, 750km, and 1000km.
C For intermediate altitudes an interpolation is used. 
C Recommended altitude range: 350-2500 km!!!
C For days between seasons centred at (21.3. = 79; 21.6. = 171;
C 23.9. 265; 21.12. = 354) relative ion density is linearly interpolated.
Cc Inputs: INVDIP - "mix" coordinate of the dip latitude and of
C                    the invariant latitude;
C                    positive northward, in deg, range <-90.0;90.0>
C         MLT - magnetic local time (central dipole)
C               in hours, range <0;24)
C         ALT - altitude above the Earth's surface;
C               in km, range <350;2000>
C         DDD - day of year; range <0;365>
C         D - coefficints of spherical harmonics for a given ion
C         ION - ion species (0...O+, 1...H+, 2...He+, 3...N+)
C Output: NION - absolute density for given ion 
C---------------------------------------------------------------------------
      REAL INVDIP,MLT,ALT,NION
      INTEGER DDD,ION
      DIMENSION  D(4,3,49),MIRREQ(49)
      REAL INVDP,DTOR
      REAL RMLT,RCOLAT
      REAL C(49),C1(82)
      INTEGER SEZA,SEZB,SEZAI,SEZBI,DDDA,DDDB,DDDD
      REAL N0A400,N0B400,N400A,N400B,N400
      REAL N0A550,N0B550,N550A,N550B,N550
      REAL N0A750,N0B750,N750A,N750B,N750
      REAL N0A100,N0B100,N100A,N100B,N1000
	REAL ANO(4),AH(4),DNO(2),ST(3)
	COMMON/CONST/DTOR,PI
	DATA (MIRREQ(J),J=1,49)/
     &            1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,
     &            1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,
     &           -1, 1,-1, 1,-1, 1, 1,-1, 1, 1,-1, 1,-1, 1, 1/
C////////////////////////////////////////////////////////////////////////////////////
C     coefficients for mirroring
      DO 10 I=1,49
       D(1,3,I)=D(1,2,I)*MIRREQ(I)
       D(2,3,I)=D(2,2,I)*MIRREQ(I)
       D(3,3,I)=D(3,2,I)*MIRREQ(I)
10     D(4,3,I)=D(4,2,I)*MIRREQ(I)
      INVDP=INVDIP

      RMLT=MLT*DTOR*15.0
      RCOLAT=(90.0-INVDP)*DTOR
      CALL SPHARM_IK(C1,6,6,RCOLAT,RMLT)
      do i=1,49
      	c(i)=c1(i)
      	enddo
C     21.3. - 20.6.
      IF ((DDD .GE. 79) .AND. (DDD .LT. 171)) THEN
       SEZA=1
       SEZB=2
       DDDA=79
       DDDB=171
       DDDD=DDD
      END IF
C     21.6. - 22.9.
      IF ((DDD .GE. 171) .AND. (DDD .LT. 265)) THEN
       SEZA=2
       SEZB=4
       DDDA=171
       DDDB=265
       DDDD=DDD
      END IF
C     23.9. - 20.12.
      IF ((DDD .GE. 265) .AND. (DDD .LT. 354)) THEN
       SEZA=4
       SEZB=3
       DDDA=265
       DDDB=354
       DDDD=DDD
      END IF
C     21.12. - 20.3.
      IF ((DDD .GE. 354) .OR. (DDD .LT. 79)) THEN
       SEZA=3
       SEZB=1
       DDDA=354
       DDDB=365+79
       DDDD=DDD
        IF (DDD .GE. 354) THEN
         DDDD=DDD
        ELSE
         DDDD=DDD+365
        END IF
      END IF
       SEZAI=MOD(SEZA-1,3)+1
       SEZBI=MOD(SEZB-1,3)+1
C     400km level
        N0A400=0.0
        N0B400=0.0
        DO 30 I=1,49
         N0A400=N0A400+C(I)*D(1,SEZAI,I)
30       N0B400=N0B400+C(I)*D(1,SEZBI,I)
        N400A=N0A400
        N400B=N0B400
        N400=(N400B-N400A)/(DDDB-DDDA)*(DDDD-DDDA)+N400A
C     550km level
      N0A550=0.0
	N0B550=0.0
	DO 70 I=1,49
         N0A550=N0A550+C(I)*D(2,SEZAI,I)
70       N0B550=N0B550+C(I)*D(2,SEZBI,I)
	N550A=N0A550
	N550B=N0B550
	N550=(N550B-N550A)/(DDDB-DDDA)*(DDDD-DDDA)+N550A
C     750km level
      N0A750=0.0
      N0B750=0.0
	DO 110 I=1,49
         N0A750=N0A750+C(I)*D(3,SEZAI,I)
110      N0B750=N0B750+C(I)*D(3,SEZBI,I)
	N750A=N0A750
	N750B=N0B750
	N750=(N750B-N750A)/(DDDB-DDDA)*(DDDD-DDDA)+N750A
C     1000km level
      N0A100=0.0
	N0B100=0.0
       DO 150 I=1,49
        N0A100=N0A100+C(I)*D(4,SEZAI,I)
150     N0B100=N0B100+C(I)*D(4,SEZBI,I)
       N100A=N0A100
       N100B=N0B100
       N1000=(N100B-N100A)/(DDDB-DDDA)*(DDDD-DDDA)+N100A

C     n(O+) AND n(N+) must not increase above 750km
      IF (((ION .EQ. 0) .OR. (ION .EQ. 3)) .AND. (N1000 .GT. N750))
     &      N1000=N750
C     n(H+) and n(He+) must not decrease above 750km
      IF (((ION .EQ. 1) .OR. (ION .EQ. 2)) .AND. (N1000 .LT. N750)) 
     &      N1000=N750

      IF (ALT .GE. 960) SUM=(N1000-N750)/220.0*(ALT-740.0)+N750     
      IF (ALT .GE. 960) GOTO 240
                
      ANO(1)=N400
	 ANO(2)=N550
	  ANO(3)=N750
	   ANO(4)=N1000

C	AH(1)=400.
C       AH(2)=550.
C        AH(3)=750.
C         AH(4)=1000.

C     exact values
	AH(1)=390.
       AH(2)=550.
        AH(3)=740.
         AH(4)=960.
      DNO(1)=20.
       DNO(2)=20.

      ST1=(ANO(2)-ANO(1))/(AH(2)-AH(1))
      DO 200 I=2,3
       ST2=(ANO(I+1)-ANO(I))/(AH(I+1)-AH(I))
       ANO(I)=ANO(I)-(ST2-ST1)*DNO(I-1)*ALOG(2.)
200   ST1=ST2

      DO 220 I=1,3
220   ST(I)=(ANO(I+1)-ANO(I))/(AH(I+1)-AH(I))

c      ARGMAX=88.0
      SUM=ANO(1)+ST(1)*(ALT-AH(1))                     
     
      DO 230 I=1,2
	aa = eptr(alt  ,dno(i),ah(i+1))
	bb = eptr(ah(1),dno(i),ah(i+1))
230   SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*DNO(I)
                
240   NION=10**SUM       
      RETURN
      END     
C
C
      SUBROUTINE IONHIGH(INVDIP,MLT,ALT,DDD,D,ION,NION)
C---------------------------------------------------------------------------
C IONHIGH calculates absolute density of O+, H+, He+ or N+  in the outer
C ionosphere for high solar activity conditions (PF10.7 = 210).
C Based on spherical harmonics approximation of relative ion density
C (by IK24) at altitudes centred on 550km, 900km, 1500km, and 2250km.
C For intermediate altitudes a interpolation is used. 
C Recommended altitude range: 500-2500 km!!!
C For days between seasons centred at (21.3. = 79; 21.6. = 171;
C 23.9. 265; 21.12. = 354) relative ion density is linearly interpolated.
Cc Inputs: INVDIP - "mix" coordinate of the dip latitude and of
C                    the invariant latitude;
C                    positive northward, in deg, range <-90.0;90.0>
C         MLT - magnetic local time (central dipole)
C               in hours, range <0;24)
C         ALT - altitude above the Earth's surface;
C               in km, range <500;3000>
C         DDD - day of year; range <0;365>
C         D - coefficints of spherical harmonics for a given ion
C         ION - ion species (0...O+, 1...H+, 2...He+, 3...N+)
C Output: NION - absolute density for given ion 
C---------------------------------------------------------------------------
      REAL INVDIP,MLT,ALT,NION
	  INTEGER DDD,ION
      DIMENSION  D(4,3,49),MIRREQ(49)
      REAL INVDP,DTOR
      REAL RMLT,RCOLAT
      REAL C(49),C1(82)
      INTEGER SEZA,SEZB,SEZAI,SEZBI,DDDA,DDDB,DDDD
      REAL N0A550,N0B550,N550A,N550B,N550
      REAL N0A900,N0B900,N900A,N900B,N900
      REAL N0A150,N0B150,N150A,N150B,N1500
      REAL N0A250,N0B250,N250A,N250B,N2500
	REAL ANO(4),AH(4),DNO(2),ST(3)
	COMMON/CONST/DTOR,PI
	DATA (MIRREQ(J),J=1,49)/
     &            1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,
     &            1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,
     &           -1, 1,-1, 1,-1, 1, 1,-1, 1, 1,-1, 1,-1, 1, 1/
C////////////////////////////////////////////////////////////////////////////////////
C     coefficients for mirroring
      DO 10 I=1,49
       D(1,3,I)=D(1,2,I)*MIRREQ(I)
       D(2,3,I)=D(2,2,I)*MIRREQ(I)
       D(3,3,I)=D(3,2,I)*MIRREQ(I)
10     D(4,3,I)=D(4,2,I)*MIRREQ(I)
      INVDP=INVDIP
      
      RMLT=MLT*DTOR*15.0
      RCOLAT=(90.0-INVDP)*DTOR
      CALL SPHARM_IK(C1,6,6,RCOLAT,RMLT)
      do i=1,49
      	c(i)=c1(i)
      	enddo
C     21.3. - 20.6.
      IF ((DDD .GE. 79) .AND. (DDD .LT. 171)) THEN
       SEZA=1
       SEZB=2
       DDDA=79
       DDDB=171
       DDDD=DDD
      END IF
C     21.6. - 22.9.
      IF ((DDD .GE. 171) .AND. (DDD .LT. 265)) THEN
       SEZA=2
       SEZB=4
       DDDA=171
       DDDB=265
       DDDD=DDD
      END IF
C     23.9. - 20.12.
      IF ((DDD .GE. 265) .AND. (DDD .LT. 354)) THEN
       SEZA=4
       SEZB=3
       DDDA=265
       DDDB=354
       DDDD=DDD
      END IF
C     21.12. - 20.3.
      IF ((DDD .GE. 354) .OR. (DDD .LT. 79)) THEN
       SEZA=3
       SEZB=1
       DDDA=354
       DDDB=365+79
       DDDD=DDD
        IF (DDD .GE. 354) THEN
         DDDD=DDD
        ELSE
         DDDD=DDD+365
        END IF
      END IF
       SEZAI=MOD(SEZA-1,3)+1
       SEZBI=MOD(SEZB-1,3)+1
C     550km level
        N0A550=0.0
        N0B550=0.0
        DO 30 I=1,49
         N0A550=N0A550+C(I)*D(1,SEZAI,I)
30       N0B550=N0B550+C(I)*D(1,SEZBI,I)
        N550A=N0A550
        N550B=N0B550
        N550=(N550B-N550A)/(DDDB-DDDA)*(DDDD-DDDA)+N550A
C     900km level
      N0A900=0.0
	N0B900=0.0
	DO 70 I=1,49
         N0A900=N0A900+C(I)*D(2,SEZAI,I)
70       N0B900=N0B900+C(I)*D(2,SEZBI,I)
	N900A=N0A900
	N900B=N0B900
	N900=(N900B-N900A)/(DDDB-DDDA)*(DDDD-DDDA)+N900A
C     1500km level
      N0A150=0.0
      N0B150=0.0
	DO 110 I=1,49
         N0A150=N0A150+C(I)*D(3,SEZAI,I)
110      N0B150=N0B150+C(I)*D(3,SEZBI,I)
	N150A=N0A150
	N150B=N0B150
	N1500=(N150B-N150A)/(DDDB-DDDA)*(DDDD-DDDA)+N150A
C     2500km level
      N0A250=0.0
	  N0B250=0.0
       DO 150 I=1,49
        N0A250=N0A250+C(I)*D(4,SEZAI,I)
150     N0B250=N0B250+C(I)*D(4,SEZBI,I)
       N250A=N0A250
       N250B=N0B250
       N2500=(N250B-N250A)/(DDDB-DDDA)*(DDDD-DDDA)+N250A

C     O+ and N+ may not increase above 1500km 
      IF (((ION .EQ. 0) .OR. (ION .EQ. 3)) .AND. (N2500 .GT. N1500)) 
     & N2500=N1500
C     H+ and He+ may not decrease above 1500km 
      IF (((ION .EQ. 1) .OR. (ION .EQ. 2)) .AND. (N2500 .LT. N1500)) 
     & N2500=N1500
              
      IF (ALT .GE. 2250.0) SUM=(N2500-N1500)/750.0*(ALT-2250.0)+N2500
      IF (ALT .GE. 2250.0) GOTO 240
      
      ANO(1)=N550
	  ANO(2)=N900
	  ANO(3)=N1500
	  ANO(4)=N2500

	  AH(1)=550.
      AH(2)=900.
      AH(3)=1500.
      AH(4)=2250.
      DNO(1)=20.
      DNO(2)=20.

      ST1=(ANO(2)-ANO(1))/(AH(2)-AH(1))
      DO 200 I=2,3
       ST2=(ANO(I+1)-ANO(I))/(AH(I+1)-AH(I))
       ANO(I)=ANO(I)-(ST2-ST1)*DNO(I-1)*ALOG(2.)
200   ST1=ST2

      DO 220 I=1,3
220   ST(I)=(ANO(I+1)-ANO(I))/(AH(I+1)-AH(I))

c      ARGMAX=88.0
      SUM=ANO(1)+ST(1)*(ALT-AH(1))                     
     
      DO 230 I=1,2
		aa = eptr(alt  ,dno(i),ah(i+1))
		bb = eptr(ah(1),dno(i),ah(i+1))
230   	SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*DNO(I)
                
240   NION=10**SUM       
      RETURN
      END
C
C
      REAL FUNCTION INVDPC(FL,DIMO,B0,DIPL)
C---------------------------------------------------------------------------
C INPUT:	FL		McIlwain L parameter, 
C			DIMO	dipole moment in Gauss
C           B0      magnetic field strength in Gauss 
C           parameters FL, DIMO and B0 are needed for invariant latitude
C			computation; uses a highly accurate polynomial expansion
C			DIPL	dip latitude in degree
C                	positive northward, in deg, range <-90.0;90.0>
C RESULT:   invdip 	a "mix" coordinate of the dip latitude (DIPL) and 
C 					of the invariant latitude; INVDIP is positive northward 
C					in degree and ranges from -90.0 to 90.0.
C---------------------------------------------------------------------------
      REAL FL,DIMO,B0,DIPL,DTOR,ASA,INVL,RINVL,RDIPL,ALFA,BETA
	  DOUBLE PRECISION B(8),A
	  COMMON/CONST/DTOR,PI
      DATA B/1.259921D0  ,-0.1984259D0 ,-0.04686632D0,-0.01314096D0,
     &      -0.00308824D0, 0.00082777D0,-0.00105877D0, 0.00183142D0/
       A=(DIMO/B0)**(1.0D0/3.0D0)/FL
       ASA=A*(B(1)+B(2)*A+B(3)*A**2+B(4)*A**3+B(5)*A**4+
     &        B(6)*A**5+B(7)*A**6+B(8)*A**7)
       IF (ASA .GT. 1.0) ASA=1.0
       IF (ASA .LT. 0.0) ASA=0.0
C      invariant latitude (absolute value)
       RINVL=ACOS(SQRT(ASA))
       INVL=RINVL/DTOR
       ALFA=2.0-(1.0/(exp((DIPL-25.0)/2.0)+1.0)+
     &           1.0/(exp((-DIPL-25.0)/2.0)+1.0))
       BETA=1.0/(exp((INVL-25.0)/2.0)+1.0)+
     &      1.0/(exp((-INVL-25.0)/2.0)+1.0)-1.0
c       INVDPC=(ALFA*INVL+BETA*DIPL)/(ALFA+BETA)
       INVDPC=(ALFA*SIGN(1.0,DIPL)*INVL+BETA*DIPL)/(ALFA+BETA)
       RETURN
       END
C
C       
      REAL FUNCTION INVDPC_OLD(FL,DIMO,B0,DIPL)
C---------------------------------------------------------------------------
C     calculation of INVDIP from FL, DIMO, BO, and DIPL
C     invariant latitude calculated by highly
C     accurate polynomial expansion
C---------------------------------------------------------------------------
      REAL FL,DIMO,B0,DIPL
      DOUBLE PRECISION B(8),A
      REAL DTOR,ASA,INVL,RINVL,RDIPL,ALFA,BETA
	  COMMON/CONST/DTOR,PI
      DATA B/1.259921D0  ,-0.1984259D0 ,-0.04686632D0,-0.01314096D0,
     &      -0.00308824D0, 0.00082777D0,-0.00105877D0, 0.00183142D0/
      A=(DIMO/B0)**(1.0D0/3.0D0)/FL
      ASA=A*(B(1)+B(2)*A+B(3)*A**2+B(4)*A**3+B(5)*A**4+
     &        B(6)*A**5+B(7)*A**6+B(8)*A**7)
      IF (ASA .GT. 1.0) ASA=1.0
      IF (ASA .LT. 0.0) ASA=0.0
c      invariant latitude (absolute value)
      RINVL=ACOS(SQRT(ASA))
      INVL=RINVL/DTOR
      RDIPL=DIPL*DTOR
      ALFA=SIN(ABS(RDIPL))**3
      BETA=COS(RINVL)**3
      INVDPC_OLD=(ALFA*SIGN(1.0,DIPL)*INVL+BETA*DIPL)/(ALFA+BETA)
      RETURN
      END

C
C                     
C*************************************************************                  
C************* PEAK VALUES ELECTRON DENSITY ******************                  
C*************************************************************                  
C
C
      real function FOUT(XMODIP,XLATI,XLONGI,UT,FF0)
c--------------------------------------------------------------
C CALCULATES CRITICAL FREQUENCY FOF2/MHZ USING SUBROUTINE GAMMA1.      
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME 
C (DEC. HOURS), FF0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
c--------------------------------------------------------------
      DIMENSION FF0(988)
      INTEGER QF(9)
      DATA QF/11,11,8,4,1,0,0,0,0/
      FOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,6,QF,9,76,13,988,FF0)
      RETURN
      END
C
C
      real function XMOUT(XMODIP,XLATI,XLONGI,UT,XM0)
c--------------------------------------------------------------
C CALCULATES PROPAGATION FACTOR M3000 USING THE SUBROUTINE GAMMA1.
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME 
C (DEC. HOURS), XM0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
c--------------------------------------------------------------
      DIMENSION XM0(441)
      INTEGER QM(7)
      DATA QM/6,7,5,2,1,0,0/
      XMOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,4,QM,7,49,9,441,XM0)
      RETURN
      END
C
C
      REAL FUNCTION HMF2ED(XMAGBR,R,X,XM3)         
c--------------------------------------------------------------
C CALCULATES THE PEAK HEIGHT HMF2/KM FOR THE MAGNETIC                           
C LATITUDE XMAGBR/DEGREE AND THE SMOOTHED ZUERICH SUNSPOT                         
C NUMBER R USING CCIR-M3000 XM3 AND THE RATIO X=FOF2/FOE.
C FOLLOWING CCIR RECOMMENDATION X IS LIMITED TO VALUE
C GREATER OR EQUAL TO 1.7 .                       
C [REF. D.BILITZA ET AL., TELECOMM.J., 46, 549-553, 1979]                       
C D.BILITZA,1980.     
c--------------------------------------------------------------
      F1=0.00232*R+0.222                         
      F2=1.2-0.0116*EXP(0.0239*R)            
      F3=0.096*(R-25.0)/150.0                      
      F4=1.0-R/150.0*EXP(-XMAGBR*XMAGBR/1600.0)
      if(x.lt.1.7) x=1.7
      DELM=F1*F4/(X-F2)+F3                
      HMF2ED=1490.0/(XM3+DELM)-176.0 
      RETURN          
      END             
C
C
      REAL FUNCTION XM3000HM(XMAGBR,R,X,HMF2)         
c--------------------------------------------------------------
C CALCULATES THE PROPAGATION FACTOR M3000 FOR THE MAGNETIC LATITUDE
C XMAGBR/DEG. AND THE SMOOTHED ZUERICH SUNSPOT NUMBER R USING THE                        
C PEAK HEIGHT HMF2/KM AND THE RATIO X=FOF2/FOE. Reverse of HMF2ED.                      
C [REF. D.BILITZA ET AL., TELECOMM.J., 46, 549-553, 1979]                       
C D.BILITZA,1980. ----- no longer used    
c--------------------------------------------------------------
      F1=0.00232*R+0.222                         
      F2=1.2-0.0116*EXP(0.0239*R)            
      F3=0.096*(R-25.0)/150.0                      
      F4=1.0-R/150.0*EXP(-XMAGBR*XMAGBR/1600.0)
      if(x.lt.1.7) x=1.7
      DELM=F1*F4/(X-F2)+F3                
	  XM3000HM=1490.0/(HMF2+176.0)-DELM
      RETURN          
      END             
C
C
	SUBROUTINE  SHAMDHMF2 (RLAT,FLON,T,RZ,HMF2)
C-------------------------------------------------------------------
C	COMPUTES THE HOURLY VALUES OF hmF2 FROM A SET OF SH COEFFICIENTS
C	IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE OF THE 
C	EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUNSPOT NUMBER.
C   PARAMETERS AND COEFFICIENTS ARE GIVEN IN DATA STATEMENTS.
C
C INPUT:	RLAT    The geographic latitude on the FLON meridian  
C                   where the modified dip latitude is equal to the 
C                   modip value for the considered point in space.
C           FLON    =LONGITUDE+15.*UT(hours)
C           T       Month as a REAL number (1.0 to 12.0)
C           RZ      12-month running mean
C OUTPUT:   HMF2    F2 peak altitude in km
C
C Altadill, D., S. Magdaleno, J.M. Torta, and E. Blanch 
C     Adv. Space Res. 52, 1756-1769, 2013.
C-------------------------------------------------------------------
      PARAMETER   (IBO=0,JBO=1,KDIM=8,LDIM=4,L=-1)
      DIMENSION   FN(0:KDIM,0:KDIM), CONST(0:KDIM,0:KDIM)
      DIMENSION   GNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			  HNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
	  DIMENSION     GANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			  HANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            GBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			  HBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      CHARACTER*1 IE       
      COMMON/AMTB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,
     *              LEXT,KMAX,FN
c     ,CONST

      DATA THETA,RE,TZERO,IFIT,ICEN,IREF,IB,KINT,LINT,KEXT,LEXT
     *     /180.,6371.2,1.0,  -1,  0,   0,   2, 8,   4,   0,   -1/

  	  DATA ((CONST(N,M), M=0,N), N=0,KDIM)
     *	 /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,
     *	  3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,
     *	  4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,
     *      2.32681,0.671693,1.,5.29150,9.72111,11.4564,9.49918,
     *      5.69951,2.42182,0.647260,1.,6.,12.5499,16.9926,16.4531,
     *      11.8645,6.40755,2.50683,0.626707/

	  DATA IE /"I"/

      DATA (((GANM(N,M,J),GBNM(N,M,J),HANM(N,M,J),HBNM(N,M,J),
     *    	J=0,LDIM), M=0,N), N=0,KDIM)
     * / 340.552,-0.376,   0.000, 0.000,
     *   -70.823, 1.150,   0.000, 0.000, -49.250, 0.591,   0.000, 0.000,
     *    24.037,-0.401,   0.000, 0.000,  -8.780, 0.282,   0.000, 0.000,
     *   -38.548, 0.177,   0.000, 0.000,
     *    52.811,-0.282,   0.000, 0.000,  35.165,-0.265,   0.000, 0.000,
     *   -21.227, 0.273,   0.000, 0.000,  10.183,-0.064,   0.000, 0.000,
     *    67.946, 1.288,  67.397, 0.850,
     *   -47.987,-1.643, -93.888,-1.257, -40.298,-1.075, -56.623,-0.768,
     *    19.471, 0.532,  37.264, 0.400,  -0.830,-0.267, -17.603,-0.243,
     *  -139.733,-1.204,   0.000, 0.000,
     *   140.567, 1.145,   0.000, 0.000,  93.356, 0.980,   0.000, 0.000,
     *   -48.353,-0.427,   0.000, 0.000,  21.370,-0.016,   0.000, 0.000,
     *   -64.759, 0.228, -83.194, 0.121,
     *    83.857,-0.228, 107.345, 0.013,  63.040,-0.268,  80.004,-0.254,
     *   -48.692, 0.109, -52.767, 0.079,  15.922, 0.004,  14.904, 0.172,
     *  -134.112, 3.405, 152.622,-2.016,
     *   192.632,-4.277,-188.201, 2.383, 114.890,-2.568,-121.624, 1.712,
     *   -69.827, 1.414,  76.826,-0.902,  35.669,-0.723, -33.598, 0.290,
     *  -241.557,-0.964,   0.000, 0.000,
     *   303.719, 1.393,   0.000, 0.000, 198.143, 0.731,   0.000, 0.000,
     *  -136.985,-0.557,   0.000, 0.000,  56.365, 0.420,   0.000, 0.000,
     *  -173.032, 0.275, -59.684,-0.368,
     *   251.737,-0.285,  79.946, 0.667, 184.674,-0.313,  30.894, 0.828,
     *  -100.746, 0.190, -14.771,-0.389,  30.175, 0.002,  21.618,-0.069,
     *     6.871, 0.189,  70.508,-0.893,
     *   -12.804, 0.018, -91.499, 1.128,  -5.759,-0.284, -60.616, 0.712,
     *    21.336, 0.036,  37.937,-0.537,  -9.778, 0.175, -15.741, 0.211,
     *    28.298, 1.499,  40.261,-2.376,
     *   -54.316,-1.889, -45.615, 2.789, -34.102,-1.211, -36.189, 2.006,
     *    20.337, 0.617,  15.090,-1.031, -10.873,-0.285,  -4.738, 0.324,
     *   -90.022, 1.611,   0.000, 0.000,
     *   139.241,-1.741,   0.000, 0.000,  69.750,-0.769,   0.000, 0.000,
     *   -28.658, 0.280,   0.000, 0.000,  22.881,-0.282,   0.000, 0.000,
     *   -44.702,-0.560,  87.636,-0.530,
     *    53.150, 0.510,-119.895, 0.539,  31.695, 0.735, -87.077, 0.547,
     *     6.348,-0.394,  54.557,-0.382,   0.875,-0.114, -15.031,-0.032,
     *  -159.602, 0.277,-144.746, 0.246,
     *   188.366,-0.118, 186.080,-0.124, 123.670,-0.172, 107.109,-0.067,
     *   -58.392, 0.093, -64.485,-0.014,  30.929,-0.054,  37.817, 0.010,
     *    19.949,-0.483,  29.250, 0.916,
     *   -30.468, 0.725, -36.526,-1.320, -13.301, 0.259, -31.395,-0.638,
     *    12.609,-0.244,  27.485, 0.374, -10.153, 0.253,  -5.509,-0.324,
     *    42.493,-0.379, 138.667,-1.337,
     *   -50.716, 0.303,-181.247, 1.744, -47.625, 0.434,-110.876, 1.038,
     *    27.315,-0.192,  57.106,-0.545,   0.778,-0.106, -30.096, 0.284,
     *   125.336, 0.325,   0.000, 0.000,
     *  -164.302,-0.458,   0.000, 0.000,-119.173, 0.120,   0.000, 0.000,
     *    78.273,-0.125,   0.000, 0.000, -26.471,-0.234,   0.000, 0.000,
     *  -164.084, 1.528,  68.950,-2.145,
     *   187.815,-1.754, -97.028, 2.521, 123.385,-1.348, -54.733, 1.333,
     *   -58.316, 0.589,  33.701,-0.805,  27.475,-0.118, -22.520, 0.539,
     *    30.638,-1.477, -71.140, 1.204,
     *   -43.070, 1.448,  86.162,-1.443, -31.966, 1.202,  68.472,-1.014,
     *    23.982,-0.552, -29.584, 0.633,  -3.900, 0.025,   4.363,-0.167,
     *    81.222, 0.150,  66.889, 0.006,
     *   -98.117,-0.095, -82.546, 0.149, -65.414,-0.057, -59.117,-0.113,
     *    33.457, 0.113,  34.682,-0.021, -13.056,-0.019, -12.982, 0.178,
     *    39.267,-0.337, -33.557, 0.109,
     *   -42.755, 0.319,  34.263,-0.255, -26.894, 0.297,  30.837,-0.086,
     *     6.342,-0.172, -17.113, 0.082,  -4.386,-0.001,   2.650,-0.097,
     *  20*0.000,
     *   -56.482,-0.951,   0.000, 0.000,
     *    49.870, 1.103,   0.000, 0.000,  47.355, 0.689,   0.000, 0.000,
     *   -21.830,-0.374,   0.000, 0.000,   0.708, 0.176,   0.000, 0.000,
     *     8.535, 0.422, -56.552, 0.790,
     *   -27.737,-0.208,  76.156,-0.834,  -9.129,-0.660,  44.140,-0.477,
     *    -5.598, 0.234, -24.694, 0.257,  -2.797, 0.226,   9.278,-0.054,
     *   -97.867, 1.578,  44.857,-0.495,
     *   116.177,-1.680, -49.651, 0.408,  93.725,-1.524, -42.861, 0.447,
     *   -55.150, 0.740,  28.612,-0.210,   7.615, 0.039,  -7.158,-0.040,
     *    84.841,-1.243, -74.812, 0.103,
     *  -100.335, 1.300,  87.349,-0.037, -78.486, 1.059,  70.399,-0.062,
     *    50.717,-0.613, -36.938, 0.044, -12.384, 0.093,   7.046, 0.005,
     *    50.573,-0.185,  47.967,-0.154,
     *   -59.881, 0.390, -55.625, 0.206, -35.187, 0.074, -39.209, 0.028,
     *    12.490,-0.039,  18.874,-0.072, -11.613, 0.191,  -8.454, 0.127,
     *  40*0.000,
     *  -243.265, 1.567,   0.000, 0.000,
     *   306.183,-1.910,   0.000, 0.000, 183.262,-0.922,   0.000, 0.000,
     *   -97.132, 0.553,   0.000, 0.000,  53.004,-0.403,   0.000, 0.000,
     *  -212.916, 1.108,   8.322,-0.089,
     *   257.690,-1.278, -15.853, 0.247, 179.077,-0.962, -15.457, 0.243,
     *   -99.826, 0.628,   1.487,-0.098,  34.120,-0.178,   1.112,-0.006,
     *   -61.739, 2.205,  62.259,-1.018,
     *    69.646,-2.288, -77.640, 1.293,  46.838,-1.850, -65.274, 0.858,
     *   -20.796, 0.919,  33.662,-0.542,   7.719,-0.165,  -5.840, 0.227,
     *    57.782, 0.176,  30.904,-0.524,
     *   -64.106,-0.272, -32.527, 0.468, -42.125,-0.123, -22.160, 0.471,
     *    21.006, 0.034,  12.768,-0.214,  -9.326,-0.077,  -3.047,-0.047,
     *   -15.644,-0.224, -96.643, 1.478,
     *    30.526, 0.117, 116.896,-1.707,   8.932, 0.254,  72.474,-1.045,
     *    -8.442,-0.098, -38.557, 0.559,  12.040,-0.100,  20.158,-0.290,
     *  60*0.000,
     *   -93.235, 1.670,   0.000, 0.000,
     *   115.863,-1.904,   0.000, 0.000,  42.475,-0.878,   0.000, 0.000,
     *   -16.805, 0.317,   0.000, 0.000,  26.766,-0.353,   0.000, 0.000,
     *   -10.938,-0.552, 229.921,-3.050,
     *     8.234, 0.386,-289.115, 3.548,   3.405, 0.746,-184.290, 2.115,
     *    12.877,-0.410, 116.519,-1.229,  -3.016,-0.147, -49.942, 0.564,
     *   -40.203,-0.273,   4.801, 0.004,
     *    46.012, 0.180, -11.896, 0.129,  38.085, 0.142,   2.665,-0.054,
     *    -6.701,-0.111,  -8.326, 0.065,   2.219,-0.065,  -3.252, 0.071,
     *    -5.474, 0.833,  61.060,-0.178,
     *     0.739,-0.777, -79.078, 0.293,  17.611,-0.764, -55.985, 0.102,
     *    -7.107, 0.359,  40.779,-0.184,  -8.184, 0.030, -14.352, 0.127,
     *   -51.003, 0.485,  36.183, 0.144,
     *    65.647,-0.752, -42.617,-0.228,  33.590,-0.298, -27.554,-0.093,
     *   -15.194, 0.193,  20.247, 0.033,  14.304,-0.227,  -6.789,-0.088,
     *  80*0.000/                                

      KMAX = MAX(KINT,KEXT)
      IF (KMAX .GT. KDIM)  GO TO 9999
      KT = MAX(LINT,LEXT)
      IF (KT .GT. LDIM)  GO TO 9999

	  DO 500 N=0,KMAX
	  DO 500 M=0,N

	  DO J=0,KT
	    GNM(N,M,J)=GANM(N,M,J)+GBNM(N,M,J)*rz
	    HNM(N,M,J)=HANM(N,M,J)+HBNM(N,M,J)*rz
	    ENDDO
 
      IF (IE .EQ. 'I')  THEN
         IF (N .GT. KINT)  GO TO 500
         LJ = LINT
      ELSE
         IF (N .GT. KEXT)  GO TO 500
         LJ = LEXT
      END IF

      FN(N,M) = FLOAT(N)

      IF (M .GT. 0)  GO TO 300
      DO 301 J=1-IBO-JBO,KT
      IF (IE .EQ. 'I')  THEN
         BINT(N,M,J)   = GNM(N,M,J)
      ELSE
         BEXT(N,M,J)   = GNM(N,M,J)
         END IF
  301 CONTINUE
      GO TO 500
  300 continue
      DO 302 J=1-IBO-JBO,LJ
      IF (IE .EQ. 'I')  THEN
         BINT(N,M,J)   = GNM(N,M,J)
         BINT(M-1,N,J) = HNM(N,M,J)
      ELSE
         BEXT(N,M,J)   = GNM(N,M,J)
         BEXT(M-1,N,J) = HNM(N,M,J)
         END IF
  302 CONTINUE
C
  500 CONTINUE
C     **********************************************************
C     SYNTHESIZES THE VALUE OF HMF2 FROM THE MODEL	
C     **********************************************************
      CALL SCHNEVPDH (RZ,RLAT,FLON,dum,T,L,dum,dum,HMF2)
	  RETURN
9999  STOP
      END
C
C
      SUBROUTINE SCHNEVPDH (RZ,FLAT,FLON,R,T,L,BN,BE,BV)
C------------------------------------------------------------------------
C WHEN L IS POSITIVE:
C COMPUTES SPHERICAL CAP HARMONIC (GEOCENTRIC) FIELD COMPONENTS
C HORIZONTAL NORTH BN,HORIZONTAL EAST BE,AND VERTICAL DOWNWARD BV.
C WHEN L IS NEGATIVE:
C COMPUTES GENERAL FUNCTION BV, ITS HORIZONTAL NORTH DERIVATIVE BN,
C AND ITS HORIZONTAL EAST DERIVATIVE BE, ON SPHERICAL CAP SURFACE.
C NOTE THAT THESE ARE METRICAL DERIVATIVES, AND BE IS THE
C LONGITUDINAL DERIVATIVE DIVIDED BY SIN(COLATITUDE).
C	  INPUT:
C FLAT,FLON,R ARE GEOCENTRIC SPHERICAL CAP LATITUDE,LONGITUDE,RADIAL
C DISTANCE; T IS TIME.
C L =  0  ON FIRST CALL:  RETURNS SPHERICAL CAP POLE POSITION FLATO,FLONO
C         AND HALF-ANGLE THETA AS BN,BE, AND BV AFTER INITIALIZATION.
C         ON SUBSEQUENT CALLS:  ACTS AS L=1.
C      1  COMPUTES POTENTIAL FIELD COMPONENTS FROM INTERNAL COEFFICIENTS.
C      2  COMPUTES POTENTIAL FIELD COMPONENTS FROM EXTERNAL COEFFICIENTS.
C      3  COMPUTES FIELD FROM BOTH INTERNAL AND EXTERNAL COEFFICIENTS.
C     -1  COMPUTES GENERAL FUNCTION BV AND DERIVATIVES BN WITH RESPECT TO
C         LATITUDE AND BE WITH RESPECT TO LONGITUDE DIVIDED BY COS(LAT)
C         (R IS DUMMY VARIABLE IN THIS CASE).
C NOTE:   SUBROUTINE IS INITIALIZED DURING FIRST CALL REGARDLESS OF L.
C
C SUBROUTINE USED:  LEGFUN
C
C PARAMS & COEFFS TRANSFERRED FROM MAIN PROGRAM IN COMMON/AMTB/
C
C ADAPTED FROM SUBROUTINE SCHNEV OF G.V. HAINES (COMPUTERS & GEOSCIENCES, 
C 14, 413-447, 1988)
C------------------------------------------------------------------------

      PARAMETER   (IBO=0,JBO=1,KDIM=8,LDIM=4)                           
      DIMENSION   FN(0:KDIM,0:KDIM), CONSTP(0:KDIM,0:KDIM)
      DIMENSION   CML(KDIM), SML(KDIM)
      DIMENSION   DELT(0:LDIM)
      DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      COMMON	/CONST/UMR,PI
     * 			/AMTB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,
     *              LEXT,KMAX,FN
C     ,CONST
      CHARACTER*1 IE,RESP
  	  DATA ((CONSTP(N,M), M=0,N), N=0,KDIM)
     *	 /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,
     *	  3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,
     *	  4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,
     *      2.32681,0.671693,1.,5.29150,9.72111,11.4564,9.49918,
     *      5.69951,2.42182,0.647260,1.,6.,12.5499,16.9926,16.4531,
     *      11.8645,6.40755,2.50683,0.626707/

C     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS BASIS FUNCTIONS
C              1          LEGENDRE POLYNOMIALS
C              2          FOURIER SERIES
C              3          COSINE SERIES
C              4          SINE SERIES
C     NOTE:    TZERO AND THINT MAY DEPEND ON IBF.
      IBF   =  2                                                        
      T1=1.
      T2=12.
      CALL TBFIT (T1,T2,IBF,THINT,TZERO)                                

C      IF (L .NE. 0)  GO TO 100
C      BN = FLATO
C      BE = FLONO
C      BV = THETA
C      RETURN

  100 IF (L .GE. 0)  THEN
         IF (IFIT .LT. 0)  GO TO 999
         AOR = RE/R
         AR = AOR**2
         IF (L .GT. 1)  GO TO 107
      ELSE
         IF (IFIT .GE. 0)  GO TO 999
         AR = -1.
         END IF
      KT = LINT
      GO TO 109
  107 IF (KEXT .GT. 0)  AOR3 = AOR*AR
      IF (L .GT. 2)  GO TO 108
      KT = LEXT
      GO TO 109
  108 KT = MAX (LINT,LEXT)
  109 DELT(0) = 1.
      IF (KT .LE. 0)  GO TO 103
      DEL = (T - TZERO)/THINT
      DO 102 I=1,KT
      IF (I .EQ. 1)  THEN
         IF (IBF .LE. 1) THEN
             DELT(I) = DEL
             ELSE IF (IBF .EQ. 2)  THEN
             ST = SIN(DEL)
             DELT(I) = ST
             ELSE IF (IBF .EQ. 3)  THEN
             DELT(I) = COS(DEL)
             ELSE
             DELT(I) = SIN(DEL)
             ENDIF
         GO TO 102
         ENDIF
      IF (IBF .EQ. 0)  THEN
          DELT(I) = DELT(I-1)*DEL
      ELSE IF (IBF .EQ. 1)  THEN
          RECIP = 1./FLOAT(I)
          DELT(I) = (2.-RECIP)*DELT(I-1)*DEL - (1.-RECIP)*DELT(I-2)
      ELSE IF (IBF .EQ. 2)  THEN
           IF ((I/2)*2 .EQ. I)  THEN
             IF (I .EQ. 2)  THEN
                CT = COS(DEL)
                DELT(I) = CT
                ELSE
                DELT(I) = DELT(I-2)*CT - DELT(I-3)*ST
                ENDIF
             ELSE
             DELT(I) = DELT(I-2)*CT + DELT(I-1)*ST
             ENDIF
      ELSE IF (IBF .EQ. 3)  THEN
          DELT(I) = COS(I*DEL)
      ELSE IF (IBF .EQ. 4)  THEN
          DELT(I) = SIN(I*DEL)
      ELSE
          GO TO 999
      ENDIF
  102 CONTINUE
      incept = 0                                                            set up
      if ((ibf.eq.2 .or. ibf.eq.3) .and. incept .eq. 1)  then
c     change to intercept form of fourier series.
          do i=2,lint,4-ibf
          delt(i) = 1. - delt(i)
          enddo
          endif
  103 X = 0.
      Y = 0.
      Z = 0.
      IF (L .EQ. 2)  GO TO 106
      IF (KINT .LT. 0)  GO TO 106
      GTI = 0.
      DO 105 I=1-IBO-JBO,LINT
  105 GTI = GTI + BINT(0,0,I)*DELT(I)
      Z = -AR*GTI
      N =  0
  106 COLAT = 90. - FLAT
      DO 150 N=1,KMAX
      IF (N .GT. 1)  GO TO 115
      CL = COS(FLON*UMR)
      SL = SIN(FLON*UMR)
      CML(1) = CL
      SML(1) = SL
      GO TO 120
  115 SML(N) = SL*CML(N-1) + CL*SML(N-1)
      CML(N) = CL*CML(N-1) - SL*SML(N-1)
  120 CONTINUE
      DO 150 M=0,N
      IF (IB .EQ. 2)  GO TO 121
      NMM = N - M
      IF ((NMM/2)*2 .NE. NMM)  GO TO 150
  121 FFN = FN(N,M)
      CALL LEGFUN (M,FFN,CONSTP(N,M),COLAT,P,DP,PMS,0)
      IF (L .GE. 0)  THEN
         AR = AOR**(FFN+2.)
      ELSE
         AR = 1.
         FFN = -2.
         DP = -DP
         PMS = -PMS
         END IF
      IF (M .NE. 0)  GO TO 130
      BT1 = 0.
      BT3 = 0.
      BT  = 0.
      IF (L .EQ. 2)  GO TO 123
      IF (N .GT. KINT)  GO TO 123
      GTI = 0.
      DO 122 I=1-IBO-JBO,LINT
  122 GTI  = GTI  + BINT(N,M,I)*DELT(I)
      BT1  = AR*GTI
      BT3  = BT1
  123 IF (L .LE. 1)  GO TO 125
      IF (N .GT. KEXT)  GO TO 125
      GTE = 0.
      DO 124 I=1-IBO-JBO,LEXT
  124 GTE = GTE + BEXT(N,M,I)*DELT(I)
      BT  = AOR3/AR*GTE
      BT1 = BT1 + BT
  125 X = X + BT1*DP
      Z = Z - (FFN*(BT3-BT)+BT3)*P
      GO TO 150
  130 BT1 = 0.
      BT2 = 0.
      BT3 = 0.
      BT  = 0.
      IF (L .EQ. 2)  GO TO 133
      IF (N .GT. KINT)  GO TO 133
      GTI = 0.
      HTI = 0.
      DO 132 I=1-IBO-JBO,LINT
      GTI = GTI + BINT(N,M,I)*DELT(I)
  132 HTI = HTI + BINT(M-1,N,I)*DELT(I)
      BT1 = AR*(GTI*CML(M) + HTI*SML(M))
      BT2 = AR*(GTI*SML(M) - HTI*CML(M))
      BT3 = BT1
  133 IF (L .LE. 1)  GO TO 135
      IF (N .GT. KEXT)  GO TO 135
      GTE = 0.
      HTE = 0.
      DO 134 I=1-IBO-JBO,LEXT
      GTE = GTE + BEXT(N,M,I)*DELT(I)
  134 HTE = HTE + BEXT(M-1,N,I)*DELT(I)
      RA = AOR3/AR
      BT = RA*(GTE*CML(M) + HTE*SML(M))
      BT1 = BT1 + BT
      BT2 = BT2 + RA*(GTE*SML(M) - HTE*CML(M))
  135 X = X + BT1*DP
      Y = Y + BT2*PMS
      Z = Z - (FFN*(BT3-BT)+BT3)*P
  150 CONTINUE
      BN = X
      BE = Y
      BV = Z
      RETURN
  999 STOP
      END
C
C
      subroutine model_hmF2(day,month,UT,xmodip,long,F107_81,hmF2)
c--------------------------------------------------------------------------
c  Input:
c    day, month (integer)
c    UT         (real) - universal time
c    xmodip     (real) - modified dip latitude in degrees
c    long       (real) - geodatic longitude    in degrees
c    F107_81    (real) - F10.7 average over 3 solar rotations
c                        (81 days, centered on the current day)
c
c Output:
c    hmF2       (real) - F2-layer peak height in km
c--------------------------------------------------------------------------
	implicit none
c     .. scalar arguments ..
      integer day, month
      real UT, xmodip, long, F107_81, hmF2
c     .. local scalars ..
      integer monthr, montha
	real hmF2_0, hmF2_m, hmF2_p, hmF2med
c     ..
	call SDMF2(UT,month,F107_81,xmodip,long,hmF2_0)

      if (day.le.15) then
         if (day.eq.15) then
             hmF2med = hmF2_0
         else
             monthr = month-1
             if (monthr.eq.0) monthr = 12
	       call SDMF2(UT,monthr,F107_81,xmodip,long,hmF2_m)
             hmF2med = (hmF2_0-(day-15)*(hmF2_m-hmF2_0)/30.)
         end if
      else
         montha = mod(month,12) + 1
	   call SDMF2(UT,montha,F107_81,xmodip,long,hmF2_p)
         hmF2med  = (hmF2_0+(day-15)*(hmF2_p-hmF2_0)/30.)
      end if

	hmF2 = hmF2med
      return
      end
C
C
      subroutine SDMF2(UT,monthut,F107A,xmodip,long,hmF2)
c--------------------------------------------------------------------------
c    Global median model of the F2-layer peak height
c
c  Requires the following subroutines and functions:
c     SDMF2, hmF2_med_SD, read_data_SD, fun_hmF2_SD,
c     fun_Gk, Legendre, fun_hmF2UT, Koeff_UT, fun_Akp_UT,
c     fun_Fk_UT, fun_Gk_UT  
c
c Author of the code:
c         Valentin Shubin
c         Pushkov Institute of Terrestrial Magnetism,
c         Ionosphere and Radio wave propagation (IZMIRAN)
c         Moscow, Troitsk, 142190, Russia
c         e-mail: shubin@izmiran.ru
c         
c     [Ref. V.N. Shubin. Global median model of the F2-layer
c     peak height based on ionospheric radio-occultation and
c     ground-based Digisonde observations. Advances in Space Research (2015)
c     http://dx.doi.org/10.1016/j.asr.2015.05.029]
c
c  Input:
c    UT      - universal time (real)
c    monthut - month (integer)
c    xmodip  - modified dip latitude in degrees (real)
c    long    - geodatic longitude    in degrees (real)
c    F107A   - F10.7 index averaged over the 3 Sun rotations 
c              in units of 10^-22 W/(m^2 Hz) (real)
c
c  Output:
c    hmF2  - F2-layer peak height in km (real)
c--------------------------------------------------------------------------
      implicit none
c     .. scalar arguments ..
      integer monthut
      real F107A
      real UT, hmF2
      real xmodip, long
c     .. local scalars ..
      integer i
      double precision T
c     .. local arrays ..
	  double precision xUT(0:23)
c	.. array in common ..
	  double precision hmF2_UT(0:23)
	  common/hmF2UT/hmF2_UT
c     .. function references .
      real hmF2_med_SD, fun_hmF2UT
c
      hmF2_UT = 0.0
	  do i=0,23
         hmF2_UT(i) = hmF2_med_SD(i,monthut,F107A,xmodip,long)
	     xUT(i) = dble(i)
         end do
c 
      T = dble(UT)
      hmF2 = fun_hmF2UT(T) 
c
      return
      end
c
c
      real function hmF2_med_SD(iUT,monthut,F107A,xmodip,long)
c---------------------------------------------------------------------
c    Input: 
c      iUT     - universal time (real)
c      monthut - month (integer)
c      F107A   - F10.7 index averaged over the 3 Sun rotations 
c                in units of 10^-22 W/(m^2 Hz) (real)
c      xmodip  - modified dip latitude in degrees (real)
c      long    - geodatic longitude    in degrees (real)
c
c    function to interpolate hmF2 between the two levels of solar activity
c    used the following auxiliary subroutines and functions:
c    read_data_SD, fun_hmF2_SD
c---------------------------------------------------------------------
      implicit none
c	..   scalar arguments ..
      integer monthut, iUT
      real F107A
      real xmodip, long
c	real xmodip, long, umr, pi
c	..   local scalars ..
      real cov, cov1, cov2
      real a, b, hmF2_1, hmF2_2
      double precision teta
c	..   local arrays ..
      double precision coeff_month(0:148,0:47)
      double precision Kf(0:148)
      real ft1(12), ft2(12)
c
c    Arrays ft1 (12) and ft2 (12) are the median values of F10.7A,
c    which were used as the margins for
c    the low and high solar activity levels      
c
c               Jan   Feb   Mar   Apr   May  Jun
	data ft1/ 73.6, 72.3, 71.8, 70.9, 73.6,73.0,
     *          71.1, 69.4, 69.1, 70.9, 72.3,74.6/
	data ft2/144.2,142.9,167.2,125.3,124.4,127.9,
     *         142.0,165.9,132.6,142.0,145.6,143.0/
c               Jul   Aug   Sep   Oct   Nov  Dec
c	.. local in common ..
	double precision umr
	common/constt/umr
c	common/const/umr,pi
c     .. function references ..
	real fun_hmF2_SD
c     .. subroutine references ..
c     read_data_SD
c
	umr=atan(1.0)*4./180
      teta = 90.0-xmodip
c
      call read_data_SD(monthut,coeff_month)
      Kf = coeff_month(0:148,iUT)
	hmF2_1 = fun_hmF2_SD(teta,long,Kf)
      Kf = coeff_month(0:148,iUT+24)
	hmF2_2 = fun_hmF2_SD(teta,long,Kf)
c
      cov = F107A
	cov1 = ft1(monthut) 
	cov2 = ft2(monthut)
c 
 	a = (hmF2_2 - hmF2_1)/log(cov2/cov1)
	b =  hmF2_2 - a*log(cov2)
	hmF2_med_SD = a*log(cov) + b
c
      return
      end
c
c
      subroutine read_data_SD(month,coeff_month)
c------------------------------------------------------------------
c    subroutine to read arrays mcsat11.datÖ mcsat22.dat
c    with coefficients of hmF2 spatial decomposition
c    for for 12 month, 24 UT hour and two solar activity levels
c------------------------------------------------------------------
	implicit none
c     .. scalar arguments ..
	integer month
c     .. array arguments ..
	double precision coeff_month(0:148,0:47)
c     .. local scalars ..
	integer coeff_month_read(1:12)
	character(256) filedata
	integer i, j
c     .. local arrays ..
	double precision coeff_month_all(0:148,0:47,1:12)
	save coeff_month_all
	data coeff_month_read /12*0/
c
      if (coeff_month_read(month) .eq. 0) then
        write(filedata, 10) month+10
        open(10, File=filedata, status='old')
	  do j=0,47
	    read(10,20) (coeff_month_all(i,j,month),i=0,148)
        end do
	  close(10)
	  coeff_month_read(month) = 1
	end if
c
      coeff_month = coeff_month_all(0:148,0:47,month)	
c
c    Previous operator can be replaced by
c
c      do i=0,148
c	   do j=0,47 
c	      coeff_month(i,j) = coeff_month_all(i,j,month)
c        end do
c	end do	   
c
	  return
	
c 10   format('data\mcsat',i2,'.dat')
 10   format('mcsat',i2,'.dat')
 20   format(6(d12.5))
      end
c
c      
      real function fun_hmF2_SD(teta,long,Kf)
c---------------------------------------------------------------------
c    Input: 
c      teta - (90-modip) in degrees
c      long - geodatic longitude in degrees
c      Kf   - coefficients of hmF2 spatial decomposition 
c
c    function to calculate spherical harmonics decomposition
c    for the spatial dependencies of hmF2
c    used the following auxiliary subroutines and functions:
c    fun_Gk, Legendre
c---------------------------------------------------------------------
      implicit none
c	.. scalar arguments ..
      double precision teta
      real long
c     .. array arguments ..
      double precision Kf(0:148)
c	.. local scalars ..
      integer k
      double precision hmF2
c	 .. local arrays ..
      double precision Gk(0:148)
c     .. subroutine references ..
c     fun_Gk
c
	call fun_Gk(teta,long,Gk)
	hmF2 = 0.d0
	do k=0,148
	   hmF2 = hmF2 + Kf(k)*Gk(k) 
	end do
      fun_hmF2_SD = hmF2
c
      return
      end
c
c
      subroutine fun_Gk(teta,long,Gk)
c---------------------------------------------------------------------
      implicit none
c	.. scalar arguments ..
      double precision teta
      real long
c      real long, umr, pi
c     .. array arguments ..
      double precision Gk(0:148)
c	.. local scalars ..
      integer mm, nn, m, n, k
c	.. local arrays ..
	double precision Pl_mn(0:8,0:12)
c	.. local in common ..
	double precision umr
	common/constt/umr
c	common/const/umr,pi
c     .. subroutine references ..
c     Legendre
c
      Pl_mn = 0.d0
	mm = 8
	nn = 12
      call Legendre(mm,nn,Pl_mn,teta)
      Gk = 0.d0
      k = 0
      do m=0,mm
         if (m==0) then
	      do n=0,nn
              Gk(k) = Pl_mn(m,n)
			k = k + 1
            end do
         else 
	      do n=m,nn
              Gk(k)   = Pl_mn(m,n)*cos(m*long*umr)
		    Gk(k+1) = Pl_mn(m,n)*sin(m*long*umr)
			k = k + 2
           end do
         end if
      end do
c
      return
      end
c
c
      subroutine Legendre(mm,nn,p,teta)
c---------------------------------------------------------------------
c     Input:
c       mm       - harmonics for longitude
c       nn       - harmonics for the modified dip latitude (modip) 
c       teta     - (90-modip) in degrees
c    Output:
c       P(mm,nn) - associated Legendre function    
c
c    subroutine to calculate associated Legendre function P(mm,nn)
c    with Schmidt normalization
c---------------------------------------------------------------------
	  implicit none
c     .. scalar arguments ..
      integer mm, nn
c      real umr, pi
      double precision teta
c     .. array arguments ..
      double precision p(0:mm,0:nn)
c     .. local scalars ..
      integer j,l,m,n
      double precision z, s
c	.. local in common ..
	  double precision umr
	  common/constt/umr
c	  common/const/umr,pi
c
      p = 0.0
	  z=cos(umr*teta)
      p(0,0)=1.
      p(0,1)=z
      if (mm.ne.0) p(1,1)=sin(umr*teta)
c
      do j=2,mm
         p(j,j)=(2*j-1)*p(j-1,j-1)*p(1,1)
      end do
c
      do m=0,mm
	   do n=1,nn
		  if (m.gt.n) then
		     p(m,n) = 0.0
			 cycle
            end if
            if ((n+1).gt.nn) exit
	      if (n+1.eq.m) cycle
		  if (m.gt.(n-1)) then
		     p(m,n+1)= (2*n+1)*z*p(m,n)/(n+1-m)
            else
		     p(m,n+1)=((2*n+1)*z*p(m,n)-(n+m)*p(m,n-1))/(n+1-m)
            end if
         end do
      end do
c
      do n=1,nn
         do m=1,mm
		  if (m.gt.n) then
		     p(m,n) = 0.0
			 exit
            end if
            s=1
	      do l=n-m+1,n+m
   	         s=s*l
            end do
	      p(m,n)=p(m,n)*sqrt(2./s)
        end do
      end do
c
      return
      end
c
c
      real function fun_hmF2UT(T)
c---------------------------------------------------------------------
c    Input: T - universal time        
c 
c    function to calculate Fourier  decomposition
c    for the temporal variations of hmF2
c    used the following auxiliary subroutines:
c    Koeff_UT, fun_Akp_UT, fun_Fk_UT, fun_Gk_UT
c---------------------------------------------------------------------
      implicit none
c	.. scalar arguments ..
      double precision T
c	.. local scalars ..
      integer k, mm, mk
c      real dtr, dumr
      double precision hmF2
c   .. local arrays ..	
	  double precision Gk_UT(0:6), Kf_UT(0:6)
c	.. local in common ..	
	  double precision dtr
	  common/radUT/dtr
c	  common/const1/dtr,dumr
c   .. subroutine references ..
c     Koeff_UT, fun_Gk_UT
c
      dtr=atan(1.0)*4.0/12.0
      mm = 3
	  mk = 2*mm
c
	  call Koeff_UT(mm,mk,Kf_UT)
      call fun_Gk_UT(mm,mk,t,Gk_UT)
	  hmF2 = 0.d0
	  do k=0,mk
	    hmF2 = hmF2 + Kf_UT(k)*Gk_UT(k) 
	  end do
      fun_hmF2UT = hmF2
c
      return
      end
c
c
      subroutine Koeff_UT(mm,mk,Kf_UT)
c---------------------------------------------------------------------
      implicit none
c	..  scalar arguments ..
      integer mm, mk
c   .. array arguments ..
	  double precision Kf_UT(0:mk)
c	.. local scalars ..
      integer k, m
      double precision sum_D
c   .. local arrays ..
	  double precision Akp_UT(0:mk,0:mk)
	  double precision Dk_UT(0:mk)
c   .. subroutine references ..
c   fun_Akp_UT 
c
      call fun_Akp_UT(mm,mk,Akp_UT,Dk_UT)
      Kf_UT = 0.d0
      do k=mk,0,-1
         sum_D = 0.d0
         do m=k+1,mk
            sum_D = sum_D + Akp_UT(m,k)*Kf_UT(m)
         end do
         Kf_UT(k) = sum_D + Dk_UT(k)
      end do   

      return
	  end
c
c
      subroutine fun_Akp_UT(mm,mk,Akp_UT,Dk_UT)
c---------------------------------------------------------------------
      implicit none
c	.. scalar arguments ..
      integer mm, mk
c   .. array arguments ..
      double precision Akp_UT(0:mk,0:mk), Dk_UT(0:mk)
c	.. local scalars ..
	  integer i, k, p
	  double precision t
	  double precision sum_An, sum_Dn
	  double precision sum_Ad, sum_Dd
c   .. local arrays ..
	  double precision Gk_UT(0:mk), Fk_UT(0:mk)
c	.. array in common ..	
	  double precision hmF2_UT(0:23)
	  common/hmF2UT/hmF2_UT
c   .. subroutine references ..
c	fun_Gk_UT, fun_Fk_UT
c
      Gk_UT = 0.d0
      Gk_UT(0) = 1.0
	  Fk_UT = 0.d0
      Fk_UT(0) = 1.d0
	  Akp_UT = 0.d0
	  Dk_UT = 0.d0
	  do p=0,mk
         sum_Dn=0.d0
         sum_Dd=0.d0
         do k=p+1,mk
            sum_An=0.d0
	        sum_Ad=0.d0
	        do i=0,23
               t = dble(i)
               call fun_Gk_UT(mm,mk,t,Gk_UT)
	           call fun_Fk_UT(mk,Gk_UT,Akp_UT,Fk_UT)
		       sum_An = sum_An + Gk_UT(k)*Fk_UT(p)
		       sum_Ad = sum_Ad + Fk_UT(p)*Fk_UT(p)
			   if (p.eq.(k-1)) then
			     sum_Dn = sum_Dn + hmF2_UT(i)*Fk_UT(p)
                 sum_Dd = sum_Dd + Fk_UT(p)*Fk_UT(p)
                 end if
            end do
            Akp_UT(k,p) = - sum_An/sum_Ad
	     end do
         if (p.lt.mk) then
	        Dk_UT(p) = sum_Dn/sum_Dd
            end if
	  end do
c
      p=mk
      sum_Dn=0.d0
      sum_Dd=0.d0
      do i=0,23
         t = dble(i)
         call fun_Gk_UT(mm,mk,t,Gk_UT)
         call fun_Fk_UT(mk,Gk_UT,Akp_UT,Fk_UT)
         sum_Dn = sum_Dn + hmF2_UT(i)*Fk_UT(p)
         sum_Dd = sum_Dd + Fk_UT(p)*Fk_UT(p)
      end do
      Dk_UT(p) = sum_Dn/sum_Dd
c
      return
      end
c
c
      subroutine fun_Fk_UT(mk,Gk_UT,Akp_UT,Fk_UT)
c---------------------------------------------------------------------
      implicit none
c	.. scalar arguments ..
	  integer mk
c   .. array arguments ..
	  double precision Gk_UT(0:mk)
	  double precision Akp_UT(0:mk,0:mk)
	  double precision Fk_UT(0:mk)
c	.. local scalars ..
      integer k, p
      double precision sum_G
c
      Fk_UT = 0.d0
      do k=0,mk
         sum_G = 0.d0
         do p=0,k
            if (k.eq.p) cycle
            sum_G = sum_G + Akp_UT(k,p)*Fk_UT(p)
         end do
         Fk_UT(k) = sum_G + Gk_UT(k)
      end do
c
      return
      end
c
c
	subroutine fun_Gk_UT(mm,mk,t,Gk_UT)
c---------------------------------------------------------------------
      implicit none
c	.. scalar arguments ..
      integer mm, mk
c      real dtr, dumr
      double precision t
c   .. array arguments ..
	  double precision Gk_UT(0:mk)
c	.. local scalars ..
      integer m, k
c	.. local in common ..
	  double precision dtr
	  common/radUT/dtr
c	  common/const1/dtr, dumr
c
	  Gk_UT = 0.d0
        k = 0
        do m=0,mm
          if (m.eq.0) then
            Gk_UT(k) = 1
            k = k + 1
          else 
            Gk_UT(k)   = cos(m*t*dtr)
            Gk_UT(k+1) = sin(m*t*dtr)
		    k = k + 2
          end if
       end do
c
      return
      end
C
C
      REAL FUNCTION FOF1ED(YLATI,R,CHI)
c--------------------------------------------------------------
C CALCULATES THE F1 PEAK PLASMA FREQUENCY (FOF1/MHZ)
C INPUT:   
C		YLATI	ABSOLUT VALUE OF DIP-LATITUDE IN DEGREE
c       R		12-MONTH RUNNING MEAN OF SUNSPOT NUMBER 
c       CHI		SOLAR ZENITH ANGLE IN DEGREE
C REFERENCE: 
c       E.D.DUCHARME ET AL., RADIO SCIENCE 6, 369-378, 1971
C                                      AND 8, 837-839, 1973
c       HOWEVER WITH MAGNETIC DIP LATITUDE INSTEAD OF GEOMAGNETIC
c       DIPOLE LATITUDE, EYFRIG, 1979                    
C--------------------------------------------- D. BILITZA, 1988.   
        COMMON/CONST/UMR,PI
	    fof1ed=0.0
	    if (chi.gt.90.0) return

        DLA =  YLATI
        F0 = 4.35 + DLA * ( 0.0058 - 1.2E-4 * DLA ) 
        F100 = 5.348 + DLA * ( 0.011 - 2.3E-4 * DLA )
        FS = F0 + ( F100 - F0 ) * R / 100.0
        XMUE = 0.093 + DLA * ( 0.0046 - 5.4E-5 * DLA ) + 3.0E-4 * R
        FOF1 = FS * COS( CHI * UMR ) ** XMUE
                CHI0 = 49.84733 + 0.349504 * DLA
                CHI100 = 38.96113 + 0.509932 * DLA
                CHIM = ( CHI0 + ( CHI100 - CHI0 ) * R / 100. )
                IF(CHI.GT.CHIM) FOF1=-FOF1 
        FOF1ED = FOF1     
        RETURN
        END             
C
C
	real function f1_c1(xmodip,hour,saxnon,suxnon)
c F1 layer shape parameter C1 after Reinisch and Huang, Advances in
c Space Research, Volume 25, Number 1, 81-88, 2000.

        common	/const/umr,pi
        pi = umr * 180.
	
        ABSMDP=ABS(XMODIP)
      	DELA=4.32
      	IF(ABSMDP.GE.18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)

      	C1OLD = 0.09 + 0.11/DELA
        if(suxnon.eq.saxnon) then
            c1 = 2.5 * c1old
        else
            c1 = 2.5*c1old*cos((HOUR-12.)/(suxnon-saxnon)*pi)
        endif
      	if(c1.lt.0.0) c1=0.0
	    f1_c1=c1
	    return
	    end
c
c
        subroutine f1_prob (sza,glat,rz12,f1prob,f1probl)
c--------------------------------------------------------------------------
c Occurrence probability of F1 layer after Scotto et al., Advances in
c Space Research, Volume 20, Number 9, 1773-1775, 1997.
c
c Input: 	sza		solar zenith angle in degrees 
c 			glat	geomagnetic latitude in degrees
C			rz12	12-month running mean of sunspot number
c Output: 	f1prob	F1 occurrence probability without L-condition cases 
c 			f1probl	F1 occurrence probability with L-condition cases
c--------------------------------------------------------------------------
c
        common /const/umr,pi

	    xarg = 0.5 + 0.5 * cos(sza*umr)
		a = 2.98 + 0.0854 * rz12
		b = 0.0107 - 0.0022 * rz12
		c = -0.000256 + 0.0000147 * rz12
		gamma = a + ( b + c * glat) * glat
	    f1pr = xarg ** gamma
        if(f1pr.lt.1.e-3) f1pr=0.0
        f1prob=f1pr
	    f1prl = xarg ** 2.36
        if(f1prl.lt.1.e-3) f1prl=0.0
        f1probl=f1prl
	    return
	    end
C
C
        REAL FUNCTION FOEEDI(COV,XHI,XHIM,XLATI)
C-------------------------------------------------------
C CALCULATES FOE/MHZ BY THE EDINBURGH-METHOD.      
C INPUT: 
C 	COV		MONTHLY MEAN 10.7CM SOLAR RADIO FLUX measured at 
C           ground level  
C   XHI		SOLAR ZENITH ANGLE IN DEGREE 
C   XHIM	SOLAR ZENITH ANGLE AT NOON IN DEGREE
C   XLATI 	ABSOLUTE VALUE OF GEOGRAPHIC LATITUDE IN DEGREE, 
C REFERENCE: 
C       KOURIS-MUGGELETON, CCIR DOC. 6/3/07, 1973
C       TROST, J. GEOPHYS. RES. 84, 2736, 1979 (was used
C               to improve the nighttime varition)
C       RAWER AND BILITZA, Adv. Space Res. 10(8), 5-14, 1990
C D.BILITZA--------------------------------- AUGUST 1986.    
        COMMON/CONST/UMR,PI
C variation with solar activity (factor A) ...............
        A=1.0+0.0094*(COV-66.0)                      
C variation with noon solar zenith angle (B) and with latitude (C)
        SL=COS(XLATI*UMR)
        IF(XLATI.LT.32.0) THEN
                SM=-1.93+1.92*SL                             
                C=23.0+116.0*SL                              
        ELSE
                SM=0.11-0.49*SL                              
                C=92.0+35.0*SL  
        ENDIF
        if(XHIM.ge.90.) XHIM=89.999
        B = COS(XHIM*UMR) ** SM
C variation with solar zenith angle (D) ..........................        
        IF(XLATI.GT.12.0) THEN
                SP=1.2
        ELSE
                SP=1.31         
        ENDIF
C adjusted solar zenith angle during nighttime (XHIC) .............
        XHIC=XHI-3.*ALOG(1.+EXP((XHI-89.98)/3.))   
        D=COS(XHIC*UMR)**SP       
C determine foE**4 ................................................
        R4FOE=A*B*C*D     
C minimum allowable foE (foe_min=sqrt[SMIN])...............................
        SMIN=0.121+0.0015*(COV-60.)
        SMIN=SMIN*SMIN
        IF(R4FOE.LT.SMIN) R4FOE=SMIN                     
        FOEEDI=R4FOE**0.25                           
        RETURN          
        END   
C
C
        REAL FUNCTION XMDED(XHI,R,YW)                
C D. BILITZA, 1978, CALCULATES ELECTRON DENSITY OF D MAXIMUM.                   
C XHI/DEG. IS SOLAR ZENITH ANGLE, R SMOOTHED ZURICH SUNSPOT NUMBER              
C AND YW/M-3 THE ASSUMED CONSTANT NIGHT VALUE.     
C [REF.: D.BILITZA, WORLD DATA CENTER A REPORT UAG-82,7,BOULDER,1981]
C corrected 4/25/97 - D. Bilitza
c
        COMMON/CONST/UMR,PI
c
        if(xhi.ge.90) goto 100
        Y = 6.05E8 + 0.088E8 * R
        yy = cos ( xhi * umr )
        yyy = -0.1 / ( yy**2.7 ) 
        if (yyy.lt.-40.) then 
        	ymd=0.0
        else
        	ymd = y * exp(yyy)
        endif
        if (ymd.lt.yw) ymd = yw
        xmded=ymd
        RETURN          

100     XMDED=YW        
        RETURN          
        END
C
C
        REAL FUNCTION GAMMA1(SMODIP,SLAT,SLONG,HOUR,
     &                          IHARM,NQ,K1,M,MM,M3,SFE)      
C---------------------------------------------------------------
C CALCULATES GAMMA1=FOF2 OR M3000 USING CCIR NUMERICAL MAP                      
C COEFFICIENTS SFE(M3) FOR MODIFIED DIP LATITUDE (SMODIP/DEG)
C GEOGRAPHIC LATITUDE (SLAT/DEG) AND LONGITUDE (SLONG/DEG)  
C AND UNIVERSIAL TIME (HOUR/DECIMAL HOURS). IHARM IS THE MAXIMUM
C NUMBER OF HARMONICS USED FOR DESCRIBING DIURNAL VARIATION.
C NQ(K1) IS AN INTEGER ARRAY GIVING THE HIGHEST DEGREES IN 
C LATITUDE FOR EACH LONGITUDE HARMONIC WHERE K1 GIVES THE NUMBER 
C OF LONGITUDE HARMONICS. M IS THE NUMBER OF COEFFICIENTS FOR 
C DESCRIBING VARIATIONS WITH SMODIP, SLAT, AND SLONG. MM IS THE
C NUMBER OF COEFFICIENTS FOR THE FOURIER TIME SERIES DESCRIBING
C VARIATIONS WITH UT.
C M=1+NQ(1)+2*[NQ(2)+1]+2*[NQ(3)+1]+... , MM=2*IHARM+1, M3=M*MM  
C SHEIKH,4.3.77.      
C---------------------------------------------------------------
      REAL*8 C(12),S(12),COEF(100),SUM             
      DIMENSION NQ(K1),XSINX(13),SFE(M3)           
      COMMON/CONST/UMR,PI
      HOU=(15.0*HOUR-180.0)*UMR                    
      S(1)=SIN(HOU)   
      C(1)=COS(HOU)   

      DO 250 I=2,IHARM                             
        C(I)=C(1)*C(I-1)-S(1)*S(I-1)                 
        S(I)=C(1)*S(I-1)+S(1)*C(I-1)                 
250     CONTINUE        

      DO 300 I=1,M    
        MI=(I-1)*MM     
        COEF(I)=SFE(MI+1)                            
        DO 300 J=1,IHARM                             
          COEF(I)=COEF(I)+SFE(MI+2*J)*S(J)+SFE(MI+2*J+1)*C(J)                       
300       CONTINUE        

      SUM=COEF(1)     
      SS=SIN(SMODIP*UMR)                           
      S3=SS           
      XSINX(1)=1.0    
      INDEX=NQ(1)     

      DO 350 J=1,INDEX                             
        SUM=SUM+COEF(1+J)*SS                         
        XSINX(J+1)=SS   
        SS=SS*S3        
350     CONTINUE        

      XSINX(NQ(1)+2)=SS                            
      NP=NQ(1)+1      
      SS=COS(SLAT*UMR)                             
      S3=SS           

      DO 400 J=2,K1   
        S0=SLONG*(J-1.)*UMR                          
        S1=COS(S0)      
        S2=SIN(S0)      
        INDEX=NQ(J)+1   
        DO 450 L=1,INDEX                             
          NP=NP+1         
          SUM=SUM+COEF(NP)*XSINX(L)*SS*S1              
          NP=NP+1         
          SUM=SUM+COEF(NP)*XSINX(L)*SS*S2              
450       CONTINUE        
        SS=SS*S3        
400     CONTINUE
        
      GAMMA1=SUM      

      RETURN          
      END 
C
C                     
C************************************************************                   
C***************** PROFILE PARAMETERS ***********************                   
C************************************************************                 
C
C

	      SUBROUTINE TOPH05(COVI,AMLAT,TIME,HMAX,HT05,SG)
C---------------------------------------------------------------------------------
C Gulyaeva T.L. (2003) Variations in the half-width of the topside ionosphere 
C    according to the observations by space ionosondes ISIS 1,ISIS 2, and IK19.
C    International J. of Geomagnetism and Aeronomy, 4(3), 201-207.
C Gulyaeva T.L., Titheridge J.E. (2006) Advanced specification of electron density 
C    and temperature in the IRI ionosphere-plasmasphere model. 
C    Adv. Space Res. 38(11), 2587-2595, doi:10.1016/j.asr.2005.08.045.
C
C  Implementation of empirical RAT=(h05top-hmF2)/hmF2 derived from ISIS and IK19
C  topside electron density profiles to obtain half peak density topside height
C  h05top  from the Chebishev polinomial coefficients given for 
C  (1) 4 levels of solar activity: Rz= 0,  50, 100, 150 replaced by
C      solar radio flux          covi=60, 106, 152, 198
C  (2) 10 selected grids of geomagnetic latitude (N=S):0,10,20,30,40,50,60,70,80,90
C  (3) 5 selected grids of local time: 0, 6, 12, 18, 24.
C  (4) 4 seasonal grids: 1 equinox(SG=90deg), 2 summer (SG=180), 
C                        3 equinox (SG=270), 4 winter(SG=360)
C   SG=season grids=90,180,270,360
C---------------------------------------------------------------------------------
      DIMENSION CVLEV(4)	
	  COMMON   /BLOCK1/HMF2,XNMF2,XHMF1,F1REG         
     *         /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU
	  DATA CVLEV/60.,106.,152.,198./
	  LOGICAL F1REG

	  ABMLAT=ABS(AMLAT)
       IR=IFIX((covi-60.)/46.)+1	
		 M1=IFIX(ABMLAT/10.)+1
			L1=IFIX(TIME/6.)+1
	  M2=M1+1
	       IF(M1.EQ.10) M2=10
        L2=L1+1
	       IF(L1.EQ.5) L2=5
C
C INTERPOLATE RAT FOR GIVEN RZI
C Call Chebishev approximation to interpolate for given ABMLAT, HRLT
C
	  call CHEBISH(CVLEV(IR),TIME,ABMLAT,XX,SG)
        IF (IR.EQ.4) THEN
	          RAT05=XX
	          GOTO 10
	          ENDIF
	  call CHEBISH(CVLEV(IR+1),TIME,ABMLAT,YY,SG)
	  RAT05=XX+(YY-XX)*(COVI-CVLEV(IR))/46.
10    HT05=HMAX*(1.+RAT05)
	  RETURN
	  END
C
C	  
	SUBROUTINE CHEBISH(COVS,HOURLT,ABMLAT,RATCH,SG)
C---------------------------------------------------------------------------------
C CHEBISHEV POLINOMIALS FOR ABMLAT(10),HOURLT(5)
C CR((C0...C5),(LT=0,6,...24),(SG=season grids=90,180,270,360)
C							(COV=60,106,152,198)
C---------------------------------------------------------------------------------
c      REAL UK(0:10),CR(0:5,5,3,4),YI(5),YY(5,3)
      REAL BR(6,5,3,4),YI(5),YY(5,3)
      REAL PL1(5),PL2(5),PL3(5),CL(0:3)
C
      COMMON /CONST/rad,pi  
      DATA PL1/-2.,-1.,0.,1.,2./
  	  DATA PL2/2.,-1.,-2.,-1.,2./
  	  DATA PL3/-1.,2.,0.,-2.,1./
      DATA BR/
C Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=60 (mlat/10=0,1,...,9)
C Equinox	B0MLAT:
     *  -1.5859,  3.5789, -3.7884, 2.7094,-1.2962,.6759
     *, -5.3449, 12.8379,-12.0165, 5.9746,-1.9084,.7669
     *,-12.8000, 35.3084,-38.0043,19.6004,-4.4974,.6975
     *,  5.8282,-13.3538,  9.1674,-0.9593,-0.8909,.6062
     *, -1.5859,  3.5789, -3.7884, 2.7094,-1.2962,.6759
C Summer	B0MLAT    
     *, -7.1103, 21.0389,-24.5539,14.1607,-3.8537,.7266
     *,  5.5333,-10.6242,  4.8751, 1.0587,-1.0821,.7527
     *,-15.4487, 42.9269,-45.0314,21.4718,-4.2116,.6026
     *, -6.6436, 16.4533,-15.5142, 6.8287,-1.2871,.4976
     *, -7.1103, 21.0389,-24.5539,14.1607,-3.8537,.7266
C Winter   B0MLAT                                                                       
     *, 14.9103,-35.2337, 27.3078,-6.5362,-0.6265,.7509
     *,  2.3846, -5.8840,  3.7023, 0.8525,-1.2663,.7086
     *, -9.8846, 26.6649,-27.0173,12.6959,-2.6536,.6295
     *,  1.7692, -2.3578, -0.7945, 2.2477,-0.9691,.5719
     *, 14.9103,-35.2337, 27.3078,-6.5362,-0.6265,.7509
C Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=106 (mlat=0,10,...,90)
C Equinox	B1MLAT
     *, -4.1218, 10.6136,-11.4922, 6.0470,-1.3620,.5563
     *,  0.9077,  2.9562, -8.9880, 6.8680,-1.9621,.7737
     *,-16.2744, 42.8047,-43.7009,20.7965,-4.0697,.6619
     *,-17.3038, 44.3336,-40.9249,15.9042,-2.1554,.4796
     *, -4.1218, 10.6136,-11.4922, 6.0470,-1.3620,.5563
C Summer	B1MLAT  
     *, -4.9692, 16.5753,-21.3543,12.7061,-3.1758,.6446
     *,  1.9000, -2.8167, -0.9962, 3.0687,-1.3454,.6859
     *,  7.6769,-14.8343,  6.7030, 1.5578,-1.0626,.4291
     *,  5.4833,-10.6322,  4.7571, 1.2178,-0.8223,.4615
     *, -4.9692, 16.5753,-21.3543,12.7061,-3.1758,.6446
C Winter	B1MLAT  
     *, -4.7282, 13.4491,-15.6931, 8.8388,-1.9732,.5874
     *,  5.6756,-14.8458, 11.8927,-2.2632,-0.6122,.6948
     *,-14.2872, 40.0829,-41.2716,18.1696,-2.7203,.4916
     *,-13.6128, 33.4657,-29.7231,11.0972,-1.2884,.5034
     *, -4.7282, 13.4491,-15.6931, 8.8388,-1.9732,.5874
C Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=152 (mlat=0,10,...,90)
C Equinox	B2MLAT
     *, -3.3282, 10.4296,-12.4722, 6.7623,-1.5172,.4931
     *, -8.9744, 20.1311,-17.4552, 7.6518,-1.7371,.6702
     *, 12.0462,-27.8932, 20.6241,-4.5781, 0.0814,.3501
     *,-17.0551, 42.3258,-37.1874,13.3608,-1.4804,.4216
     *, -3.3282, 10.4296,-12.4722, 6.7623,-1.5172,.4931
C Summer	B2MLAT  
     *,  7.3077,-17.1579, 11.6872,-0.7405,-1.0298,.5754
     *, 19.2641,-45.1886, 34.3297,-8.1879,-0.1875,.6562
     *,  6.0987,-11.0903,  4.3569, 1.4001,-0.7309,.3885
     *,  5.9295,-13.9205, 10.2347,-2.2818, 0.0853,.3915
     *,  7.3077,-17.1579, 11.6872,-0.7405,-1.0298,.5754
C Winter	B2MLAT  
     *, -1.6821,  8.6010,-13.6570, 8.6307,-1.9846,.5635
     *,  5.4679,-12.3750,  7.5620, 0.5394,-1.4415,.6659
     *, -8.0821, 21.9288,-21.8597, 9.3455,-1.4644,.3599
     *, -8.3000, 19.3076,-16.3295, 6.1619,-0.9144,.3846
     *, -1.6821,  8.6010,-13.6570, 8.6307,-1.9846,.5635
C Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=198 (mlat=0,10,...,90)
C Equinox	B3MLAT
     *,-16.4051, 28.2880,-16.0982, 4.6328,-1.0405,.5486
     *, 13.0846,-34.8291, 30.0074,-8.6402, 0.1529,.6165
     *, 19.7474,-42.7116, 28.9430,-6.0487, 0.1492,.3748
     *, 16.2795,-36.6982, 26.5094,-6.3492, 0.2926,.3946
     *,-16.4051, 28.2880,-16.0982, 4.6328,-1.0405,.5486
C Summer	B3MLAT
     *,  4.6410,-13.7931, 11.6548,-1.9248,-0.7246,.5264
     *, -2.4090,  3.1805, -2.8423, 2.8861,-0.9937,.5234
     *,  6.3410,-13.9643,  8.2461,-0.0186,-0.7009,.3582
     *,  9.0987,-20.8618, 14.7262,-2.8798,-0.0512,.3662
     *,  4.6410,-13.7931, 11.6548,-1.9248,-0.7246,.5264
C Winter	B3MLAT
     *, -4.6526, 12.1878,-14.4047, 8.5226,-2.0493,.5903
     *,  3.9821, -6.9477,  0.8382, 3.4898,-1.5694,.6283
     *, -7.0474, 17.3974,-17.3465, 8.3671,-1.5708,.3759
     *,  4.2782, -9.9880,  5.9834, 0.0975,-0.4900,.3842
     *, -4.6526, 12.1878,-14.4047, 8.5226,-2.0493,.5903/

C	DATA UL/-2.,-1.,0.,1.,2./

	do k=0,3
	cl(k)=0.
	enddo
C
	IR=IFIX((covs-60.)/46.)+1	
C Given geomagnetic latitude parameter:
	xi=abmlat/100.
	DO LS=1,3
          DO LL=1,5
      B1=BR(6,LL,LS,IR)      
      B2=BR(5,LL,LS,IR)        
      B3=BR(4,LL,LS,IR)       
      B4=BR(3,LL,LS,IR)        
	B5=BR(2,LL,LS,IR)        
	B6=BR(1,LL,LS,IR)        
       HLT=(LL-1)*6.0

      YY(LL,LS)=B1+xi*(B2+xi*(B3+xi*(B4+xi*(B5+xi*B6))))
	ENDDO
	ENDDO            ! end of season/day cycle
C Apply seasonal interpolation
	 do i=1,5
	p0=(2.*YY(i,1)+YY(i,2)+YY(i,3))/4.
	p1=(YY(i,3)-YY(i,2))/2.
	p2=(YY(i,2)+YY(i,3)-2.*YY(i,1))/4.
	YI(i)=p0+p1*cos(sg*rad)+p2*cos(2.*sg*rad)
	 enddo
      DO K=1,5
      CL(0)=CL(0)+YI(K)
      CL(1)=CL(1)+YI(K)*PL1(K)
      CL(2)=CL(2)+YI(K)*PL2(K)
      CL(3)=CL(3)+YI(K)*PL3(K)
	ENDDO
      CL(0)=CL(0)/5.
      CL(1)=CL(1)/10.
      CL(2)=CL(2)/14.
      CL(3)=CL(3)/12.
      ULL=(HOURLT-12.)/6.
      ZA=CL(0)-2.*CL(2)
      RATCH=ZA+ULL*(CL(1)-3.4*CL(3)+ULL*(CL(2)+ULL*CL(3)))

	RETURN
	END	
C
C
	SUBROUTINE  SHAMDB0D (RLAT,FLON,T,RZ,B)
C-------------------------------------------------------------------
C	COMPUTES THE HOURLY VALUES OF B0 FROM A SET OF SH COEFFICIENTS
C	IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE
C	OF THE EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUSPOT NUMER
C
C INPUT:	RLAT    The geogrphic latitude on the meridian given by 
C					the local time (FLON), where the modified dip
C                   latitude is the same as of the orginal site.
C			FLON	=LONGITUDE+15.*UT(hours)
C			T		Month as a REAL number (1.0 to 12.0)
C			RZ		12-month running mean
C OUTOUT	B		=B0
C
C  Blanch E., D. Arrazola, D. Altadill, D. Buresova, M. Mosert, 
C     Adv. Space Res. 39, 701-710, 2007.
C  Altadill, D., D. Arrazola, E. Blanch, D. Buresova, 
C     Adv. Space Res. 42, 610-616, 2008.
C  Altadill, D., J.M. Torta, and E. Blanch, 
C     Adv. Space Res. 43,1825-1834, 2009.
C-------------------------------------------------------------------
      PARAMETER   (IBO=0,JBO=1,KDIM=6,LDIM=4,L=-1)
      DIMENSION   FN(0:KDIM,0:KDIM), CONST(0:KDIM,0:KDIM)
      DIMENSION   GNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			  HNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
	  DIMENSION   GANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			  HANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            GBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			  HBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      CHARACTER*1 IE       
      COMMON/ATB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,
     *              LEXT,KMAX,FN

      DATA THETA,RE,TZERO,IFIT,ICEN,IREF,IB,KINT,LINT,KEXT,LEXT
     *     /180.,6371.2,1.0,  -1,  0,   0,   2, 6,   4,   0,   -1/

  	  DATA ((CONST(N,M), M=0,N), N=0,KDIM)
     *	 /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,
     *	  3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,
     *	  4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,
     *      2.32681,0.671693/

	  DATA IE /"I"/

      DATA (((GANM(N,M,J),GBNM(N,M,J),HANM(N,M,J),HBNM(N,M,J),
     *    	J=0,LDIM), M=0,N), N=0,KDIM)
     * /176.746, 0.233,   0.000, 0.000,
     * -109.413, 0.072,   0.000, 0.000, -66.000,-0.004,   0.000, 0.000,
     *   36.874,-0.018,   0.000, 0.000, -19.515, 0.040,   0.000, 0.000,
     *   94.998, 0.724,   0.000, 0.000,
     * -116.900,-0.971,   0.000, 0.000, -93.849,-0.590,   0.000, 0.000,
     *   80.579, 0.425,   0.000, 0.000, -19.205,-0.220,   0.000, 0.000,
     *  -94.824, 0.115,   6.055, 0.265,
     *   84.720,-0.161,  -7.101,-0.374,  35.200,-0.138,   1.043,-0.350,
     *  -23.960, 0.109,  -2.478, 0.133,  25.550,-0.049,  -3.143, 0.003,
     *  -29.418,-0.823,   0.000, 0.000,
     *   40.070, 0.800,   0.000, 0.000,  24.019, 0.478,   0.000, 0.000,
     *  -13.175,-0.265,   0.000, 0.000,   8.799, 0.090,   0.000, 0.000,
     *  -31.099,-1.117,  -1.906, 0.498,
     *   43.807, 1.406,  -3.216,-0.520,  37.957, 0.902,  -0.717,-0.570,
     *  -40.232,-0.583,  12.171, 0.244,  11.595, 0.241,  -4.890, 0.054,
     *  -87.665, 1.635, 117.581,-1.765,
     *  123.444,-2.119,-146.917, 2.131,  81.137,-1.173, -99.063, 1.548,
     *  -42.646, 0.681,  61.263,-0.811,  17.550,-0.408, -24.374, 0.260,
     *   54.538,-0.170,   0.000, 0.000,
     *  -71.552, 0.361,   0.000, 0.000, -50.565,-0.077,   0.000, 0.000,
     *   36.653,-0.071,   0.000, 0.000, -10.816, 0.236,   0.000, 0.000,
     *  -31.138, 1.156,  37.307,-1.407,
     *   40.390,-1.390, -34.573, 1.730,  41.597,-0.835, -41.318, 1.550,
     *  -19.779, 0.404,  15.954,-0.696,  -1.706,-0.220,   5.084, 0.040,
     *  -57.671, 0.045,  42.683,-0.800,
     *   71.491, 0.048, -49.441, 0.980,  47.893,-0.037, -36.191, 0.562,
     *  -26.638,-0.029,  20.346,-0.384,   9.998, 0.067,  -6.787, 0.213,
     *   90.187,-1.198, -66.032,-0.056,
     * -119.148, 1.428,  81.202, 0.022, -63.375, 0.754,  53.070, 0.149,
     *   39.249,-0.436, -30.898,-0.052, -27.293, 0.301,  12.838,-0.067,
     * -110.261, 1.509,   0.000, 0.000,
     *  164.956,-1.761,   0.000, 0.000, 103.699,-1.005,   0.000, 0.000,
     *  -55.127, 0.569,   0.000, 0.000,  25.376,-0.315,   0.000, 0.000,
     * -104.655, 1.341, 109.057,-1.367,
     *  139.129,-1.730,-127.325, 1.532,  88.526,-1.068,-106.461, 1.397,
     *  -38.306, 0.508,  56.240,-0.798,  17.239,-0.267,  -7.766, 0.058,
     *   -6.494,-1.253,   5.714, 0.132,
     *    3.547, 1.545,  -5.372,-0.106,  -4.343, 1.103,  -3.393,-0.017,
     *    2.454,-0.626,  -3.297,-0.025,   5.871, 0.160,   2.040,-0.036,
     *   50.814,-0.230, -25.094, 0.817,
     *  -65.502, 0.304,  32.267,-1.075, -44.176, 0.019,  14.606,-0.605,
     *   27.869,-0.009,  -5.147, 0.387, -11.041, 0.131,   5.922,-0.225,
     *   77.825,-0.728, 128.501,-0.810,
     *  -87.685, 0.838,-164.016, 1.103, -74.431, 0.807, -95.539, 0.498,
     *   40.631,-0.454,  49.950,-0.292,  -4.229, 0.000, -29.666, 0.272,
     *  152.380,-1.232,   0.000, 0.000,
     * -192.098, 1.514,   0.000, 0.000,-132.417, 1.370,   0.000, 0.000,
     *   82.894,-0.709,   0.000, 0.000, -28.162, 0.050,   0.000, 0.000,
     *  -12.633, 1.192,  47.246,-1.193,
     *   -5.488,-1.387, -67.206, 1.486,  -9.917,-0.914, -34.438, 0.552,
     *   13.185, 0.477,  21.225,-0.387,   0.586,-0.208, -15.426, 0.419,
     *   -4.478,-0.118,  17.908, 0.175,
     *   -0.417, 0.067, -27.047,-0.241,   7.636, 0.028, -10.075,-0.109,
     *  -10.582, 0.005,  14.496, 0.086,   0.421, 0.001, -12.200,-0.041,
     *   16.086, 0.321,  47.044,-0.126,
     *  -24.823,-0.280, -62.615, 0.210, -12.030,-0.136, -44.003,-0.023,
     *    4.929, 0.137,  28.340,-0.009,  -4.688,-0.057,  -9.315, 0.103,
     *   28.023,-0.031, -21.535, 0.115,
     *  -31.946, 0.011,  24.143,-0.180, -21.019,-0.057,  24.108,-0.116,
     *   13.969, 0.004, -13.823, 0.042,  -6.860, 0.031,   0.546,-0.035,
     * 20*0.000,
     *  -31.994, 0.409,   0.000, 0.000,
     *   18.217,-0.458,   0.000, 0.000,  39.280,-0.754,   0.000, 0.000,
     *  -20.453, 0.324,   0.000, 0.000,  -8.111, 0.139,   0.000, 0.000,
     *   28.765,-0.477, -28.368, 0.516,
     *  -50.604, 0.751,  25.725,-0.471, -23.444, 0.283,  29.966,-0.558,
     *    2.759,-0.146, -10.824, 0.341,  -7.419, 0.206,  -3.711, 0.056,
     *   42.429,-0.415,   1.993, 0.117,
     *  -53.162, 0.555,   7.229,-0.246, -19.307, 0.039,  -8.028, 0.028,
     *    9.849,-0.035,   6.834, 0.033, -17.010, 0.272,   4.668,-0.129,
     *    4.546,-0.359, -57.796, 0.359,
     *    0.738, 0.343,  73.027,-0.423,  -7.421, 0.420,  56.067,-0.327,
     *    5.093,-0.279, -37.581, 0.226,   3.636,-0.041,  10.910,-0.059,
     *   88.440,-0.393, -69.598, 0.643,
     * -109.481, 0.532,  82.266,-0.765, -59.229, 0.182,  55.279,-0.580,
     *   28.514,-0.057, -30.282, 0.326, -22.924, 0.164,  11.602,-0.073,
     * 40*0.000/                                

      KMAX = MAX(KINT,KEXT)
      IF (KMAX .GT. KDIM)  GO TO 9999
      KT = MAX(LINT,LEXT)
      IF (KT .GT. LDIM)  GO TO 9999

	DO 500 N=0,KMAX
	DO 500 M=0,N

	DO J=0,KT
	 GNM(N,M,J)=GANM(N,M,J)+GBNM(N,M,J)*rz
	 HNM(N,M,J)=HANM(N,M,J)+HBNM(N,M,J)*rz
	ENDDO
 
      IF (IE .EQ. 'I')  THEN
         IF (N .GT. KINT)  GO TO 500
         LJ = LINT
      ELSE
         IF (N .GT. KEXT)  GO TO 500
         LJ = LEXT
         END IF

      FN(N,M) = FLOAT(N)

      IF (M .GT. 0)  GO TO 300
      DO 301 J=1-IBO-JBO,KT
      IF (IE .EQ. 'I')  THEN
         BINT(N,M,J)   = GNM(N,M,J)
      ELSE
         BEXT(N,M,J)   = GNM(N,M,J)
         END IF
  301 CONTINUE
      GO TO 500
  300 continue
      DO 302 J=1-IBO-JBO,LJ
      IF (IE .EQ. 'I')  THEN
         BINT(N,M,J)   = GNM(N,M,J)
         BINT(M-1,N,J) = HNM(N,M,J)
      ELSE
         BEXT(N,M,J)   = GNM(N,M,J)
         BEXT(M-1,N,J) = HNM(N,M,J)
         END IF
  302 CONTINUE
C
  500 CONTINUE
C     **********************************************************
C     SYNTHESIZES THE VALUE OF B0 FROM THE MODEL	
C     **********************************************************
      CALL SCHNEVPD(RZ,RLAT,FLON,dum,T,L,dum,dum,B)
	  RETURN
9999  STOP
      END
C
C
      SUBROUTINE SCHNEVPD (RZ,FLAT,FLON,R,T,L,BN,BE,BV)
C-------------------------------------------------------------------
C     WHEN L IS POSITIVE:
C     COMPUTES SPHERICAL CAP HARMONIC (GEOCENTRIC) FIELD COMPONENTS
C     HORIZONTAL NORTH BN,HORIZONTAL EAST BE,AND VERTICAL DOWNWARD BV.
C     WHEN L IS NEGATIVE:
C     COMPUTES GENERAL FUNCTION BV, ITS HORIZONTAL NORTH DERIVATIVE BN,
C     AND ITS HORIZONTAL EAST DERIVATIVE BE, ON SPHERICAL CAP SURFACE.
C     NOTE THAT THESE ARE METRICAL DERIVATIVES, AND BE IS THE
C     LONGITUDINAL DERIVATIVE DIVIDED BY SIN(COLATITUDE).

C     FLAT,FLON,R ARE GEOCENTRIC SPHERICAL CAP LATITUDE,LONGITUDE,RADIAL
C     DISTANCE; T IS TIME.

C     L =  0  ON FIRST CALL:  RETURNS SPHERICAL CAP POLE POSITION FLATO,FLONO
C             AND HALF-ANGLE THETA AS BN,BE, AND BV AFTER INITIALIZATION.
C             ON SUBSEQUENT CALLS:  ACTS AS L=1.
C          1  COMPUTES POTENTIAL FIELD COMPONENTS FROM INTERNAL COEFFICIENTS.
C          2  COMPUTES POTENTIAL FIELD COMPONENTS FROM EXTERNAL COEFFICIENTS.
C          3  COMPUTES FIELD FROM BOTH INTERNAL AND EXTERNAL COEFFICIENTS.
C         -1  COMPUTES GENERAL FUNCTION BV AND DERIVATIVES BN WITH RESPECT TO
C             LATITUDE AND BE WITH RESPECT TO LONGITUDE DIVIDED BY COS(LAT)
C             (R IS DUMMY VARIABLE IN THIS CASE).
C     NOTE:   SUBROUTINE IS INITIALIZED DURING FIRST CALL REGARDLESS OF L.

C     SUBPROGRAM USED:  LEGFUN

C	***** PARAMS & COEFFS TRANSFERRED FROM MAIN PROGRAM *****

C	ADAPTED FROM SUBROUTINE SCHNEV OF G.V. HAINES (COMPUTERS & GEOSCIENCES, 
C      14, 413-447, 1988)
C-------------------------------------------------------------------

      PARAMETER   (IBO=0,JBO=1,KDIM=6,LDIM=4)                                     
      DIMENSION   FN(0:KDIM,0:KDIM), CONSTP(0:KDIM,0:KDIM)
      DIMENSION   CML(KDIM), SML(KDIM)
      DIMENSION   DELT(0:LDIM)
      DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      COMMON	/CONST/dfarg,PI 
     * 			/ATB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,
     *              LEXT,KMAX,FN
      CHARACTER*1 IE,RESP
C
      DATA ((CONSTP(N,M), M=0,N), N=0,KDIM)
     *	 /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,
     *	  3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,
     *	  4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,
     *      2.32681,0.671693/

C     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS BASIS FUNCTIONS
C              1          LEGENDRE POLYNOMIALS
C              2          FOURIER SERIES
C              3          COSINE SERIES
C              4          SINE SERIES
C     NOTE:    TZERO AND THINT MAY DEPEND ON IBF.
      IBF   =  2                                                       
      T1=1.
      T2=12.
      CALL TBFIT (T1,T2,IBF,THINT,TZERO)                                

C      IF (L .NE. 0)  GO TO 100
C      BN = FLATO
C      BE = FLONO
C      BV = THETA
C      RETURN

  100 IF (L .GE. 0)  THEN
         IF (IFIT .LT. 0)  GO TO 999
         AOR = RE/R
         AR = AOR**2
         IF (L .GT. 1)  GO TO 107
      ELSE
         IF (IFIT .GE. 0)  GO TO 999
         AR = -1.
         END IF
      KT = LINT
      GO TO 109
  107 IF (KEXT .GT. 0)  AOR3 = AOR*AR
      IF (L .GT. 2)  GO TO 108
      KT = LEXT
      GO TO 109
  108 KT = MAX (LINT,LEXT)
  109 DELT(0) = 1.
      IF (KT .LE. 0)  GO TO 103
      DEL = (T - TZERO)/THINT
      DO 102 I=1,KT
      IF (I .EQ. 1)  THEN
         IF (IBF .LE. 1) THEN
             DELT(I) = DEL
             ELSE IF (IBF .EQ. 2)  THEN
             ST = SIN(DEL)
             DELT(I) = ST
             ELSE IF (IBF .EQ. 3)  THEN
             DELT(I) = COS(DEL)
             ELSE
             DELT(I) = SIN(DEL)
             ENDIF
         GO TO 102
         ENDIF
      IF (IBF .EQ. 0)  THEN
          DELT(I) = DELT(I-1)*DEL
      ELSE IF (IBF .EQ. 1)  THEN
          RECIP = 1./FLOAT(I)
          DELT(I) = (2.-RECIP)*DELT(I-1)*DEL - (1.-RECIP)*DELT(I-2)
      ELSE IF (IBF .EQ. 2)  THEN
           IF ((I/2)*2 .EQ. I)  THEN
             IF (I .EQ. 2)  THEN
                CT = COS(DEL)
                DELT(I) = CT
                ELSE
                DELT(I) = DELT(I-2)*CT - DELT(I-3)*ST
                ENDIF
             ELSE
             DELT(I) = DELT(I-2)*CT + DELT(I-1)*ST
             ENDIF
      ELSE IF (IBF .EQ. 3)  THEN
          DELT(I) = COS(I*DEL)
      ELSE IF (IBF .EQ. 4)  THEN
          DELT(I) = SIN(I*DEL)
      ELSE
          GO TO 999
      ENDIF
  102 CONTINUE
      incept = 0                                                              
      if ((ibf.eq.2 .or. ibf.eq.3) .and. incept .eq. 1)  then
c     change to intercept form of fourier series.
          do i=2,lint,4-ibf
          delt(i) = 1. - delt(i)
          enddo
          endif
  103 X = 0.
      Y = 0.
      Z = 0.
      IF (L .EQ. 2)  GO TO 106
      IF (KINT .LT. 0)  GO TO 106
      GTI = 0.
      DO 105 I=1-IBO-JBO,LINT
  105 GTI = GTI + BINT(0,0,I)*DELT(I)
      Z = -AR*GTI
      N =  0
  106 COLAT = 90. - FLAT
      DO 150 N=1,KMAX
      IF (N .GT. 1)  GO TO 115
      CL = COS(FLON*dfarg)
      SL = SIN(FLON*dfarg)
      CML(1) = CL
      SML(1) = SL
      GO TO 120
  115 SML(N) = SL*CML(N-1) + CL*SML(N-1)
      CML(N) = CL*CML(N-1) - SL*SML(N-1)
  120 CONTINUE
      DO 150 M=0,N
      IF (IB .EQ. 2)  GO TO 121
      NMM = N - M
      IF ((NMM/2)*2 .NE. NMM)  GO TO 150
  121 FFN = FN(N,M)
      CALL LEGFUN (M,FFN,CONSTP(N,M),COLAT,P,DP,PMS,0)
      IF (L .GE. 0)  THEN
         AR = AOR**(FFN+2.)
      ELSE
         AR = 1.
         FFN = -2.
         DP = -DP
         PMS = -PMS
         END IF
      IF (M .NE. 0)  GO TO 130
      BT1 = 0.
      BT3 = 0.
      BT  = 0.
      IF (L .EQ. 2)  GO TO 123
      IF (N .GT. KINT)  GO TO 123
      GTI = 0.
      DO 122 I=1-IBO-JBO,LINT
  122 GTI  = GTI  + BINT(N,M,I)*DELT(I)
      BT1  = AR*GTI
      BT3  = BT1
  123 IF (L .LE. 1)  GO TO 125
      IF (N .GT. KEXT)  GO TO 125
      GTE = 0.
      DO 124 I=1-IBO-JBO,LEXT
  124 GTE = GTE + BEXT(N,M,I)*DELT(I)
      BT  = AOR3/AR*GTE
      BT1 = BT1 + BT
  125 X = X + BT1*DP
      Z = Z - (FFN*(BT3-BT)+BT3)*P
      GO TO 150
  130 BT1 = 0.
      BT2 = 0.
      BT3 = 0.
      BT  = 0.
      IF (L .EQ. 2)  GO TO 133
      IF (N .GT. KINT)  GO TO 133
      GTI = 0.
      HTI = 0.
      DO 132 I=1-IBO-JBO,LINT
      GTI = GTI + BINT(N,M,I)*DELT(I)
  132 HTI = HTI + BINT(M-1,N,I)*DELT(I)
      BT1 = AR*(GTI*CML(M) + HTI*SML(M))
      BT2 = AR*(GTI*SML(M) - HTI*CML(M))
      BT3 = BT1
  133 IF (L .LE. 1)  GO TO 135
      IF (N .GT. KEXT)  GO TO 135
      GTE = 0.
      HTE = 0.
      DO 134 I=1-IBO-JBO,LEXT
      GTE = GTE + BEXT(N,M,I)*DELT(I)
  134 HTE = HTE + BEXT(M-1,N,I)*DELT(I)
      RA = AOR3/AR
      BT = RA*(GTE*CML(M) + HTE*SML(M))
      BT1 = BT1 + BT
      BT2 = BT2 + RA*(GTE*SML(M) - HTE*CML(M))
  135 X = X + BT1*DP
      Y = Y + BT2*PMS
      Z = Z - (FFN*(BT3-BT)+BT3)*P
  150 CONTINUE
      BN = X
      BE = Y
      BV = Z
      RETURN
  999 STOP
      END
C
C
      SUBROUTINE  SHAB1D (FLAT,FLON,T,RZ,B)
C-------------------------------------------------------------------
C	COMPUTES THE HOURLY VALUES OF B1 FROM A SET OF SH COEFFICIENTS
C	IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE
C	OF THE EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUSPOT NUMER
C
C   PARAMETERS ARE THE SAME AS IN SHAMDB0D, EXCEPT:
C		FLAT	Geographic latitude
C		B		=B1
C
C	***** PARAMS & COEFFS IN DATA SENTENCES *****
C-------------------------------------------------------------------
C
      PARAMETER   (IBO=0,JBO=1,KDIM=6,LDIM=4,L=-1)
      DIMENSION   FN(0:KDIM,0:KDIM), CONST(0:KDIM,0:KDIM)
      DIMENSION   GNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			HNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
	  DIMENSION   GANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			HANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            GBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *			HBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      CHARACTER*1 IE       
      COMMON/ATB1/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,
     *              LEXT,KMAX,FN

	  DATA ALT /300./

      DATA THETA,RE,TZERO,IFIT,ICEN,IREF,IB,KINT,LINT,KEXT,LEXT
     *     /180.,6371.2,1.0,  -1,  0,   0,   2, 6,   4,   0,   -1/

  	  DATA ((CONST(N,M), M=0,N), N=0,KDIM)
     *	 /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,
     *	  3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,
     *	  4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,
     *      2.32681,0.671693/

	  DATA IE /"I"/

      DATA (((GANM(N,M,J),GBNM(N,M,J),HANM(N,M,J),HBNM(N,M,J),
     *    	J=0,LDIM), M=0,N), N=0,KDIM)
     *   /1.156, 0.039,  0.000, 0.000,
     *    1.725,-0.053,  0.000, 0.000,  1.097,-0.032,  0.000, 0.000,
     *   -0.579, 0.019,  0.000, 0.000,  0.265,-0.010,  0.000, 0.000,
     *   -2.895, 0.023,  0.000, 0.000,
     *    3.269,-0.025,  0.000, 0.000,  2.278,-0.015,  0.000, 0.000,
     *   -1.789, 0.009,  0.000, 0.000,  0.653,-0.005,  0.000, 0.000,
     *   -3.240, 0.052, -2.645, 0.030,
     *    4.404,-0.062,  3.283,-0.038,  2.827,-0.038,  2.181,-0.019,
     *   -1.496, 0.020, -1.143, 0.011,  0.688,-0.011,  0.512,-0.008,
     *   -0.023,-0.025,  0.000, 0.000,
     *   -0.370, 0.031,  0.000, 0.000, -0.385, 0.017,  0.000, 0.000,
     *    0.150,-0.009,  0.000, 0.000,  0.019, 0.007,  0.000, 0.000,
     *    1.704, 0.006, -1.766, 0.022,
     *   -2.115,-0.007,  2.309,-0.028, -1.610,-0.003,  1.314,-0.015,
     *    0.806, 0.002, -0.873, 0.010, -0.249,-0.002,  0.445,-0.005,
     *    0.256,-0.009, -2.959, 0.023,
     *   -0.662, 0.013,  3.643,-0.030, -0.208, 0.001,  2.208,-0.015,
     *    0.144,-0.001, -1.326, 0.010, -0.205, 0.005,  0.678,-0.007,
     *    1.730,-0.030,  0.000, 0.000,
     *   -2.072, 0.035,  0.000, 0.000, -1.027, 0.016,  0.000, 0.000,
     *    0.946,-0.008,  0.000, 0.000, -0.644, 0.009,  0.000, 0.000,
     *    3.925,-0.060,  1.613,-0.015,
     *   -4.790, 0.070, -2.021, 0.019, -3.293, 0.048, -1.273, 0.006,
     *    1.829,-0.026,  0.797,-0.005, -0.767, 0.011, -0.395, 0.006,
     *   -1.988, 0.027, -0.761, 0.003,
     *    2.258,-0.031,  0.978,-0.004,  1.772,-0.026,  0.459,-0.003,
     *   -0.813, 0.014, -0.257, 0.004,  0.173,-0.003,  0.223,-0.001,
     * 20*0.000,
     *   -1.356, 0.011,  0.000, 0.000,
     *    0.079,-0.003,  0.000, 0.000,  0.415,-0.012,  0.000, 0.000,
     *   -0.249, 0.004,  0.000, 0.000, -0.087, 0.004,  0.000, 0.000,
     *   -1.155,-0.012,  0.261,-0.026,
     *    1.619, 0.011, -0.217, 0.032,  0.989, 0.009, -0.361, 0.019,
     *   -0.745,-0.001,  0.009,-0.009,  0.347, 0.000,  0.178, 0.005,
     *    4.672,-0.032,  0.562, 0.017,
     *   -5.868, 0.040, -0.850,-0.020, -3.798, 0.026, -0.145,-0.020,
     *    2.079,-0.014,  0.160, 0.008, -0.943, 0.007, -0.417, 0.001,
     * 40*0.000,
     *    2.477, 0.000,  0.000, 0.000,
     *   -1.815,-0.011,  0.000, 0.000, -1.571, 0.002,  0.000, 0.000,
     *    0.551, 0.004,  0.000, 0.000, -0.044,-0.008,  0.000, 0.000,
     *    2.160,-0.010, -0.305, 0.012,
     *   -2.618, 0.010,  0.529,-0.016, -1.782, 0.006,  0.634,-0.007,
     *    0.976,-0.003, -0.449, 0.006, -0.331, 0.001, -0.004,-0.004,
     *   -0.394, 0.002,  0.851,-0.020,
     *    0.359,-0.003, -1.051, 0.024,  0.357, 0.002, -0.239, 0.012,
     *    0.005,-0.001,  0.210,-0.008, -0.028,-0.003, -0.322, 0.005,
     * 60*0.000,
     *   -6.760, 0.064,  0.000, 0.000,
     *    7.700,-0.073,  0.000, 0.000,  5.394,-0.054,  0.000, 0.000,
     *   -2.788, 0.026,  0.000, 0.000,  0.923,-0.007,  0.000, 0.000,
     *   -2.328, 0.024, -0.463, 0.020,
     *    2.923,-0.027,  0.490,-0.025,  1.768,-0.019,  0.711,-0.017,
     *   -1.068, 0.009, -0.363, 0.010,  0.596,-0.004, -0.073,-0.004,
     *   -1.911, 0.016, -4.519, 0.041,
     *    2.644,-0.024,  5.569,-0.050,  1.287,-0.009,  3.707,-0.031,
     *   -0.894, 0.007, -2.121, 0.019,  0.669,-0.007,  0.933,-0.010,
     * 80*0.000/ 

      KMAX = MAX(KINT,KEXT)
      IF (KMAX .GT. KDIM)  GO TO 9999
      KT = MAX(LINT,LEXT)
      IF (KT .GT. LDIM)  GO TO 9999

	DO 500 N=0,KMAX
	DO 500 M=0,N

	DO J=0,KT
	 GNM(N,M,J)=GANM(N,M,J)+GBNM(N,M,J)*rz
	 HNM(N,M,J)=HANM(N,M,J)+HBNM(N,M,J)*rz
	ENDDO
 
      IF (IE .EQ. 'I')  THEN
         IF (N .GT. KINT)  GO TO 500
         LJ = LINT
      ELSE
         IF (N .GT. KEXT)  GO TO 500
         LJ = LEXT
         END IF

      FN(N,M) = FLOAT(N)

      IF (M .GT. 0)  GO TO 300
  255 FORMAT (1X,A1,2I3,F9.4,E15.6,F10.3,4F20.3:/(22X,5F20.3))
      DO 301 J=1-IBO-JBO,KT
      IF (IE .EQ. 'I')  THEN
         BINT(N,M,J)   = GNM(N,M,J)
      ELSE
         BEXT(N,M,J)   = GNM(N,M,J)
         END IF
  301 CONTINUE
      GO TO 500
  300 continue
  260 FORMAT (1X,A1,2I3,F9.4,E15.6,10F10.3:/(32X,10F10.3))
      DO 302 J=1-IBO-JBO,LJ
      IF (IE .EQ. 'I')  THEN
         BINT(N,M,J)   = GNM(N,M,J)
         BINT(M-1,N,J) = HNM(N,M,J)
      ELSE
         BEXT(N,M,J)   = GNM(N,M,J)
         BEXT(M-1,N,J) = HNM(N,M,J)
         END IF
  302 CONTINUE

  500 CONTINUE
C     **********************************************************
C     SYNTHESIZES THE VALUE OF B1 FROM THE MODEL	
C     **********************************************************
      CALL SCHNEVPDB1(RZ,FLAT,FLON,dum,T,L,dum,dum,B)
C
      RETURN
9999  STOP
      END
C
C
      SUBROUTINE SCHNEVPDB1 (RZ,FLAT,FLON,R,T,L,BN,BE,BV)
C-------------------------------------------------------------------
C     WHEN L IS POSITIVE:
C     COMPUTES SPHERICAL CAP HARMONIC (GEOCENTRIC) FIELD COMPONENTS
C     HORIZONTAL NORTH BN,HORIZONTAL EAST BE,AND VERTICAL DOWNWARD BV.
C     WHEN L IS NEGATIVE:
C     COMPUTES GENERAL FUNCTION BV, ITS HORIZONTAL NORTH DERIVATIVE BN,
C     AND ITS HORIZONTAL EAST DERIVATIVE BE, ON SPHERICAL CAP SURFACE.
C     NOTE THAT THESE ARE METRICAL DERIVATIVES, AND BE IS THE
C     LONGITUDINAL DERIVATIVE DIVIDED BY SIN(COLATITUDE).

C     FLAT,FLON,R ARE GEOCENTRIC SPHERICAL CAP LATITUDE,LONGITUDE,RADIAL
C     DISTANCE; T IS TIME.

C     L =  0  ON FIRST CALL:  RETURNS SPHERICAL CAP POLE POSITION FLATO,FLONO
C             AND HALF-ANGLE THETA AS BN,BE, AND BV AFTER INITIALIZATION.
C             ON SUBSEQUENT CALLS:  ACTS AS L=1.
C          1  COMPUTES POTENTIAL FIELD COMPONENTS FROM INTERNAL COEFFICIENTS.
C          2  COMPUTES POTENTIAL FIELD COMPONENTS FROM EXTERNAL COEFFICIENTS.
C          3  COMPUTES FIELD FROM BOTH INTERNAL AND EXTERNAL COEFFICIENTS.
C         -1  COMPUTES GENERAL FUNCTION BV AND DERIVATIVES BN WITH RESPECT TO
C             LATITUDE AND BE WITH RESPECT TO LONGITUDE DIVIDED BY COS(LAT)
C             (R IS DUMMY VARIABLE IN THIS CASE).
C     NOTE:   SUBROUTINE IS INITIALIZED DURING FIRST CALL REGARDLESS OF L.

C     SUBPROGRAM USED:  LEGFUN

C	***** PARAMS & COEFFS TRANSFERRED FROM MAIN PROGRAM *****

C	ADAPTED FROM SUBROUTINE SCHNEV OF G.V. HAINES (COMPUTERS & GEOSCIENCES, 
C      14, 413-447, 1988)
C-------------------------------------------------------------------

      PARAMETER   (IBO=0,JBO=1,KDIM=6,LDIM=4)                                     
      DIMENSION   FN(0:KDIM,0:KDIM), CONSTP(0:KDIM,0:KDIM)
      DIMENSION   CML(KDIM), SML(KDIM)
      DIMENSION   DELT(0:LDIM)
      DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),
     *            BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
      COMMON	/CONST/dfarg,PI 
     * 			/ATB1/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,
     *              LEXT,KMAX,FN
      CHARACTER*1 IE,RESP
C
      DATA ((CONSTP(N,M), M=0,N), N=0,KDIM)
     *	 /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,
     *	  3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,
     *	  4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,
     *      2.32681,0.671693/

C     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS BASIS FUNCTIONS
C              1          LEGENDRE POLYNOMIALS
C              2          FOURIER SERIES
C              3          COSINE SERIES
C              4          SINE SERIES
C     NOTE:    TZERO AND THINT MAY DEPEND ON IBF.
      IBF   =  2                                                       
      T1=1.
      T2=12.
      CALL TBFIT (T1,T2,IBF,THINT,TZERO)                                

C      IF (L .NE. 0)  GO TO 100
C      BN = FLATO
C      BE = FLONO
C      BV = THETA
C      RETURN

  100 IF (L .GE. 0)  THEN
         IF (IFIT .LT. 0)  GO TO 999
         AOR = RE/R
         AR = AOR**2
         IF (L .GT. 1)  GO TO 107
      ELSE
         IF (IFIT .GE. 0)  GO TO 999
         AR = -1.
         END IF
      KT = LINT
      GO TO 109
  107 IF (KEXT .GT. 0)  AOR3 = AOR*AR
      IF (L .GT. 2)  GO TO 108
      KT = LEXT
      GO TO 109
  108 KT = MAX (LINT,LEXT)
  109 DELT(0) = 1.
      IF (KT .LE. 0)  GO TO 103
      DEL = (T - TZERO)/THINT
      DO 102 I=1,KT
      IF (I .EQ. 1)  THEN
         IF (IBF .LE. 1) THEN
             DELT(I) = DEL
             ELSE IF (IBF .EQ. 2)  THEN
             ST = SIN(DEL)
             DELT(I) = ST
             ELSE IF (IBF .EQ. 3)  THEN
             DELT(I) = COS(DEL)
             ELSE
             DELT(I) = SIN(DEL)
             ENDIF
         GO TO 102
         ENDIF
      IF (IBF .EQ. 0)  THEN
          DELT(I) = DELT(I-1)*DEL
      ELSE IF (IBF .EQ. 1)  THEN
          RECIP = 1./FLOAT(I)
          DELT(I) = (2.-RECIP)*DELT(I-1)*DEL - (1.-RECIP)*DELT(I-2)
      ELSE IF (IBF .EQ. 2)  THEN
           IF ((I/2)*2 .EQ. I)  THEN
             IF (I .EQ. 2)  THEN
                CT = COS(DEL)
                DELT(I) = CT
                ELSE
                DELT(I) = DELT(I-2)*CT - DELT(I-3)*ST
                ENDIF
             ELSE
             DELT(I) = DELT(I-2)*CT + DELT(I-1)*ST
             ENDIF
      ELSE IF (IBF .EQ. 3)  THEN
          DELT(I) = COS(I*DEL)
      ELSE IF (IBF .EQ. 4)  THEN
          DELT(I) = SIN(I*DEL)
      ELSE
          GO TO 999
      ENDIF
  102 CONTINUE
      incept = 0                                                              
      if ((ibf.eq.2 .or. ibf.eq.3) .and. incept .eq. 1)  then
c     change to intercept form of fourier series.
          do i=2,lint,4-ibf
          delt(i) = 1. - delt(i)
          enddo
          endif
  103 X = 0.
      Y = 0.
      Z = 0.
      IF (L .EQ. 2)  GO TO 106
      IF (KINT .LT. 0)  GO TO 106
      GTI = 0.
      DO 105 I=1-IBO-JBO,LINT
  105 GTI = GTI + BINT(0,0,I)*DELT(I)
      Z = -AR*GTI
      N =  0
  106 COLAT = 90. - FLAT
      DO 150 N=1,KMAX
      IF (N .GT. 1)  GO TO 115
      CL = COS(FLON*dfarg)
      SL = SIN(FLON*dfarg)
      CML(1) = CL
      SML(1) = SL
      GO TO 120
  115 SML(N) = SL*CML(N-1) + CL*SML(N-1)
      CML(N) = CL*CML(N-1) - SL*SML(N-1)
  120 CONTINUE
      DO 150 M=0,N
      IF (IB .EQ. 2)  GO TO 121
      NMM = N - M
      IF ((NMM/2)*2 .NE. NMM)  GO TO 150
  121 FFN = FN(N,M)
      CALL LEGFUN (M,FFN,CONSTP(N,M),COLAT,P,DP,PMS,0)
      IF (L .GE. 0)  THEN
         AR = AOR**(FFN+2.)
      ELSE
         AR = 1.
         FFN = -2.
         DP = -DP
         PMS = -PMS
         END IF
      IF (M .NE. 0)  GO TO 130
      BT1 = 0.
      BT3 = 0.
      BT  = 0.
      IF (L .EQ. 2)  GO TO 123
      IF (N .GT. KINT)  GO TO 123
      GTI = 0.
      DO 122 I=1-IBO-JBO,LINT
  122 GTI  = GTI  + BINT(N,M,I)*DELT(I)
      BT1  = AR*GTI
      BT3  = BT1
  123 IF (L .LE. 1)  GO TO 125
      IF (N .GT. KEXT)  GO TO 125
      GTE = 0.
      DO 124 I=1-IBO-JBO,LEXT
  124 GTE = GTE + BEXT(N,M,I)*DELT(I)
      BT  = AOR3/AR*GTE
      BT1 = BT1 + BT
  125 X = X + BT1*DP
      Z = Z - (FFN*(BT3-BT)+BT3)*P
      GO TO 150
  130 BT1 = 0.
      BT2 = 0.
      BT3 = 0.
      BT  = 0.
      IF (L .EQ. 2)  GO TO 133
      IF (N .GT. KINT)  GO TO 133
      GTI = 0.
      HTI = 0.
      DO 132 I=1-IBO-JBO,LINT
      GTI = GTI + BINT(N,M,I)*DELT(I)
  132 HTI = HTI + BINT(M-1,N,I)*DELT(I)
      BT1 = AR*(GTI*CML(M) + HTI*SML(M))
      BT2 = AR*(GTI*SML(M) - HTI*CML(M))
      BT3 = BT1
  133 IF (L .LE. 1)  GO TO 135
      IF (N .GT. KEXT)  GO TO 135
      GTE = 0.
      HTE = 0.
      DO 134 I=1-IBO-JBO,LEXT
      GTE = GTE + BEXT(N,M,I)*DELT(I)
  134 HTE = HTE + BEXT(M-1,N,I)*DELT(I)
      RA = AOR3/AR
      BT = RA*(GTE*CML(M) + HTE*SML(M))
      BT1 = BT1 + BT
      BT2 = BT2 + RA*(GTE*SML(M) - HTE*CML(M))
  135 X = X + BT1*DP
      Y = Y + BT2*PMS
      Z = Z - (FFN*(BT3-BT)+BT3)*P
  150 CONTINUE
      BN = X
      BE = Y
      BV = Z
      RETURN
  999 STOP
      END
C
C
      SUBROUTINE TBFIT (T1,T2,IBF,THINT,TZERO)
C-------------------------------------------------------------------
C	COURTESY OF G.V. HAINES
C
C     T2    =  BEGINNING OF TIME INTERVAL.
C     T1    =  END OF TIME INTERVAL.
C     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS TEMPORAL BASIS FUNCTIONS
C              1          LEGENDRE POLYNOMIALS
C              2          FOURIER SERIES
C              3          COSINE SERIES
C              4          SINE SERIES
C              5          COSINE + SINE SERIES
C     TZERO =  TIME-TRANSLATION PARAMETER:
C              FOR IBF.LE.1, CHOOSE TZERO = CENTER OF TIME INTERVAL
C              FOR IBF.GE.2, CHOOSE TZERO = BEGINNING OF TIME INTERVAL.
C     THINT =  TIME-SCALING PARAMETER. THINT = HALF OF TIME INTERVAL T2-T1
C              FOR IBF.LE.1; PI*HALF OF TIME INTERVAL FOR IBF.EQ.2;
C              AND PI*TIME INTERVAL FOR IBF.GE.3.
C     NOTE:    CHOOSING TZERO AND THINT IN THIS WAY SCALES TIME
C              TO (-1,1) FOR IBF.LE.1;  TO (0,2PI) FOR IBF.EQ.2;
C              AND TO (0,PI) FOR IBF.GE.3.
C-------------------------------------------------------------------

      IBFI  = IBF
      IF (IBFI .LE. 1)  THEN
          TZERO = (T2+T1)/2.D0
         ELSE
          TZERO =  T1
         ENDIF
      THINT = T2 - T1
      IF (IBFI .LE. 2)  THINT = THINT/2.D0
C      IF (IBFI .EQ. 4)  THEN
C          DELT(0) = 0.D0
C         ELSE
C          DELT(0) = 1.D0
C         ENDIF
      RETURN
      END
C
C
      SUBROUTINE LEGFUN (M,FN,CONST,COLAT,P,DP,PMS,IPRT)
C-------------------------------------------------------------------
C     SERIES FORM FOR ASSOCIATED LEGENDRE FUNCTION P, ITS DERIVATIVE DP,
C     AND THE FUNCTION PMS=P*M/SIN(COLAT), IN POWERS OF (1-COS(COLAT))/2.
C     INTEGRAL ORDER M, REAL DEGREE FN, NORMALIZING CONSTANT CONST.
C     COLATITUDE COLAT IN DEGREES.
C     IPRT = 0     NO PRINT-OUT
C            1     PRINT PARAMETERS AND P SERIES
C            2     PRINT PARAMETERS AND DP SERIES
C            3     PRINT PARAMETERS AND BOTH P AND DP SERIES
C           -1     PRINT PARAMETERS ONLY
C     INPUT M,FN,CONST,COLAT,IPRT.   OUTPUT P,DP,PMS
C	ADAPTED FROM G.V. HAINES (COMPUTERS & GEOSCIENCES, 14, 413-447, 1988).
C-------------------------------------------------------------------
      REAL*8  		FNN,AL,A,B,PNM,DPNM
      DIMENSION  	AM(60), BM(60)
      LOGICAL		mess
	  COMMON/iounit/konsol,mess	/CONST/dfarg,PI

      DATA   JMAX/60/

c      dfarg=(atan(1.0)*4.)/180.
      FNN = FN*(FN+1.)
      IF (COLAT .LT. 60.)  THEN
          X = SIN(dfarg*COLAT/2.)**2
          C = 1. - 2.*X
      ELSE
          C = COS(COLAT*dfarg)
          X = (1. - C)/2.
          END IF
      S = SIN(COLAT*dfarg)
      IF (M .GT.1)  GO TO 20
      IF (M .LT. 0)  STOP
      AL = CONST
      GO TO 50
   20 AL = CONST*S**(M-1)
   50 PNM = AL
      DPNM = 0.
      J = 0
  100 J = J + 1
      JPM = J + M
      B = AL*((JPM-1)-FNN/JPM)
      DPNM = DPNM + B
      A = (B*X)/J
      PNM = PNM + A
      AL = A
C     STORE P OR DP SERIES FOR PRINTOUT.
      IF (IPRT .LE. 0)  GO TO 150
      IF (IPRT .EQ. 2)  GO TO 145
      AM(J) = A
      IF (IPRT .EQ. 1)  GO TO 150
  145 BM(J) = B
C     CHECK FOR TERMINATION OF SERIES.
  150 ABSA = ABS(A)
      ABSB = ABS(B)
      IF (ABSB .GE. 1.E-7)  GO TO 160
      IF (ABSA .LT. 1.E-7)  GO TO 110
  160 IF (ABSB .GE. 1.E+13)  GO TO 105
      IF (ABSA .GE. 1.E+13)  GO TO 105
C     CHANGE CHECK LIMITS ACCORDING TO ACCURACY DESIRED AND ACCORDING
C     TO WORD SIZE OF COMPUTER.
C     FOR 32-BIT WORD, DOUBLE PRECISION, E-8 AND E+8 GIVE 7 DIGITS ACCURACY.
C     FOR 60-BIT WORD, DOUBLE PRECISION, E-15 AND E+15 GIVE 14 DIGITS ACCURACY.
C     FOR 60-BIT WORD, SINGLE PRECISION, E-8 AND E+7 GIVE 7 DIGITS ACCURACY.
C     (DOUBLE OR SINGLE PRECISION REFER TO FNN,AL,A,B,PNM,DPNM)
      IF (J .LT. JMAX)  GO TO 100
C     CONVERGENCE SLOW OR JMAX TOO SMALL
  105 CONTINUE
C     NUMERICAL ERROR UNACCEPTABLY LARGE DUE TO ADDING OF
C     LARGE AND SMALL NUMBERS.
      if (mess) WRITE(konsol,108) M,FN,CONST,J,A,B
  108 FORMAT (//12H ** ERROR **/1X,I5,F10.5,E15.7,I5,2D15.7)
      STOP
C     SERIES TRUNCATED SUCCESSFULLY.
  110 PS = PNM
      DPS = DPNM
      IF (M .NE. 0)  GO TO 115
      PMS = 0.
      P = PS
      DP = DPS*S/2.
      GO TO 120
  115 PMS = PS*M
      P = PS*S
      DP = DPS*S*S/2. + C*PMS
  120 CONTINUE
C     PRINT TERMS OF SERIES
      IF (IPRT .EQ. 0)  RETURN
      if (mess) WRITE(konsol,125) M,FN,CONST,COLAT,P,DP,PMS,J
  125 FORMAT (/1X,I5,F10.5,E20.12,F10.2,3F25.14,I5)
      IF (IPRT .LT. 0)  RETURN
      IF (IPRT .EQ. 2)  GO TO 135
      if (mess) write(konsol,130) (AM(I),I=1,J)
  130 FORMAT (1X,16E8.1)
      IF (IPRT .EQ. 1)  RETURN
  135 CONTINUE
C  135 PRINT *, (BM(I),I=1,J)
      RETURN
      END
C      
C      
        REAL FUNCTION B0_98 ( HOUR, SAX, SUX, NSEASN, R, ZLO, ZMODIP)
C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON is northern season with
C ISEASON=1 northern spring), low and high solar activity Rz12=10,
C 100 (IR=1,2), and modified dip latitudes of 0, 18 and 45
C degress (ILATI=1,2,3). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C
C input:
C       hour    LT in decimal hours
C       SAX     time of sunrise in decimal hours
C       SUX     time of sunset in decimal hours
C       nseasn  season in northern hemisphere (1=spring)
C       R       12-month running mean of sunspot number
C       ZLO     longitude
C       ZMODIP  modified dip latitude
C
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
C Updates (B0_new -> B0_98):
C
C 01/98 corrected to include a smooth transition at the modip equator
C       and no discontinuity at the equatorial change in season.
C 09/98 new B0 values incl values at the magnetic equator
C 10/98 longitude as input to determine if magnetic equator in northern 
C         or southern hemisphere
C-------------------------------------------------------------------
      REAL      NITVAL
      DIMENSION B0F(2,4,2,3),bfr(2,2,3),bfd(2,3),zx(5),g(6),dd(5)
      DATA      B0F/201,68,210,61,192,68,199,67,240,80,245,83,
     &              233,71,230,65,108,65,142,81,110,68,77,75,
     &              124,98,164,100,120,94,96,112,78,81,94,84,
     &              81,81,65,70,102,87,127,91,109,88,81,78/
      DATA      zx/45.,72.,90.,108.,135./,dd/5*3.0/

        num_lat=3

C jseasn is southern hemisphere season
        jseasn=nseasn+2
        if(jseasn.gt.4) jseasn=jseasn-4

        zz = zmodip + 90.
        zz0 = 0.

C Interpolation in Rz12: linear from 10 to 100
        DO 7035 ISL=1,num_lat
          DO 7034 ISD=1,2
            bfr(isd,1,isl) = b0f(isd,nseasn,1,isl) +
     &      (b0f(isd,nseasn,2,isl) - b0f(isd,nseasn,1,isl))/90.*(R-10.)
            bfr(isd,2,isl) = b0f(isd,jseasn,1,isl) +
     &      (b0f(isd,jseasn,2,isl) - b0f(isd,jseasn,1,isl))/90.*(R-10.)
7034      continue
C Interpolation day/night with transitions at SAX (sunrise)
C and SUX (sunset) for northern/southern hemisphere iss=1/2
          do 7033 iss=1,2
                DAYVAL = BFR(1,ISS,ISL)
                NITVAL = BFR(2,ISS,ISL)
                BFD(iss,ISL) = HPOL(HOUR,DAYVAL,NITVAL,SAX,SUX,1.,1.)
7033      continue
7035    continue

C Interpolation with epstein-transitions in modified dip latitude.
C Transitions at +/-18 and +/-45 degrees; constant above +/-45.
C
C g(1:5) are the latitudinal slopes of B0;
C       g(1) is for the region from -90 to -45 degrees
C       g(2) is for the region from -45 to -18 degrees
C       g(3) is for the region from -18 to   0 degrees
C       g(4) is for the region from   0 to  18 degrees
C       g(5) is for the region from  18 to  45 degrees
C       g(6) is for the region from  45 to  90 degrees
C
C B0 =  bfd(2,3) at modip = -45,
C       bfd(2,2) at modip = -18,
C       bfd(2,1) or bfd(1,1) at modip = 0,
C       bfd(1,2) at modip = 20,
C       bfd(1,3) at modip = 45.
C If the Longitude is between 200 and 320 degrees then the modip 
C equator is in the southern hemisphere and bfd(2,1) is used at the 
C equator, otherwise bfd(1,1) is used.
c
        zx1=bfd(2,3)
        zx2=bfd(2,2)
        zx3=bfd(1,1)
        if(zlo.gt.200.0.and.zlo.lt.320) zx3=bfd(2,1)
        zx4=bfd(1,2)
        zx5=bfd(1,3)
        g(1) = 0.
        g(2) = ( zx2 - zx1 ) / 27.
        g(3) = ( zx3 - zx2 ) / 18.
        g(4) = ( zx4 - zx3 ) / 18.
        g(5) = ( zx5 - zx4 ) / 27.
        g(6) = 0.

c        bb0 = bfd(2,3)
c      SUM = bb0
        sum=zx1
      DO 1 I=1,5
        aa = eptr(zz ,dd(i),zx(i))
        bb = eptr(zz0,dd(i),zx(i))
        DSUM = (G(I+1) - G(I)) * (AA-BB) * dd(i)
        SUM = SUM + DSUM
1       continue
      B0_98 = SUM

        RETURN
        END
C
C
      SUBROUTINE TAL(SHABR,SDELTA,SHBR,SDTDH0,AUS6,SPT)                         
C-----------------------------------------------------------
C CALCULATES THE COEFFICIENTS SPT FOR THE POLYNOMIAL
C Y(X)=1+SPT(1)*X**2+SPT(2)*X**3+SPT(3)*X**4+SPT(4)*X**5               
C TO FIT THE VALLEY IN Y, REPRESENTED BY:                
C Y(X=0)=1, THE X VALUE OF THE DEEPEST VALLEY POINT (SHABR),                    
C THE PRECENTAGE DEPTH (SDELTA), THE WIDTH (SHBR) AND THE                       
C DERIVATIVE DY/DX AT THE UPPER VALLEY BOUNDRY (SDTDH0).                        
C IF THERE IS AN UNWANTED ADDITIONAL EXTREMUM IN THE VALLEY                     
C REGION, THEN AUS6=.TRUE., ELSE AUS6=.FALSE..     
C FOR -SDELTA THE COEFF. ARE CALCULATED FOR THE FUNCTION                        
C Y(X)=EXP(SPT(1)*X**2+...+SPT(4)*X**5).           
C-----------------------------------------------------------
      DIMENSION SPT(4)                             
      LOGICAL AUS6    
      AUS6=.FALSE.
      if(SHBR.le.0.0) then
          AUS6=.TRUE.
          RETURN
          ENDIF
      Z1=-SDELTA/(100.0*SHABR*SHABR)               
      IF(SDELTA.GT.0.) GOTO 500                    
      SDELTA=-SDELTA  
      Z1=ALOG(1.-SDELTA/100.)/(SHABR*SHABR)        
500   Z3=SDTDH0/(2.*SHBR)                          
      Z4=SHABR-SHBR   
      SPT(4)=2.0*(Z1*(SHBR-2.0*SHABR)*SHBR+Z3*Z4*SHABR)/                        
     &  (SHABR*SHBR*Z4*Z4*Z4)                        
      SPT(3)=Z1*(2.0*SHBR-3.0*SHABR)/(SHABR*Z4*Z4)-
     &  (2.*SHABR+SHBR)*SPT(4)          
      SPT(2)=-2.0*Z1/SHABR-2.0*SHABR*SPT(3)-3.0*SHABR*SHABR*SPT(4)              
      SPT(1)=Z1-SHABR*(SPT(2)+SHABR*(SPT(3)+SHABR*SPT(4)))                      
      B=4.*SPT(3)/(5.*SPT(4))+SHABR                
      C=-2.*SPT(1)/(5*SPT(4)*SHABR)                
      Z2=B*B/4.-C     
      IF(Z2.LT.0.0) GOTO 300                       
      Z3=SQRT(Z2)     
      Z1=B/2.         
      Z2=-Z1+Z3       
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
      IF (ABS(Z3).GT.1.E-15) GOTO 400              
      Z2=C/Z2         
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
      RETURN          
400   Z2=-Z1-Z3       
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
300   RETURN          
      END             
C
C
        SUBROUTINE VALGUL(XHI,HVB,VWU,VWA,VDP)
C --------------------------------------------------------------------- 
C   CALCULATES E-F VALLEY PARAMETERS; T.L. GULYAEVA, ADVANCES IN
C   SPACE RESEARCH 7, #6, 39-48, 1987.
C
C       INPUT:  XHI     SOLAR ZENITH ANGLE [DEGREE]
C       
C       OUTPUT: VDP     VALLEY DEPTH  (NVB/NME)
C               VWU     VALLEY WIDTH  [KM]
C               VWA     VALLEY WIDTH  (SMALLER, CORRECTED BY RAWER)
C               HVB     HEIGHT OF VALLEY BASE [KM]
C -----------------------------------------------------------------------
C
        COMMON  /CONST/UMR,PI
C
        CS = 0.1 + COS(UMR*XHI)
        ABC = ABS(CS)
        VDP = 0.45 * CS / (0.1 + ABC ) + 0.55
        ARL = ( 0.1 + ABC + CS ) / ( 0.1 + ABC - CS)
        ZZZ = ALOG( ARL )
        VWU = 45. - 10. * ZZZ
        VWA = 45. -  5. * ZZZ
        HVB = 1000. / ( 7.024 + 0.224 * CS + 0.966 * ABC )
        RETURN
        END
C
C
      Subroutine DRegion(z,it,f,vKp,f5SW,f6WA,elg)
c-----------------------------------------------------------------------
c Reference: Danilov, Rodevich, and Smirnova, Adv. Space Res.  
C     15, #2, 165, 1995.
C
C Input:     z    - solar zenith angle in degrees
C            it   - season (month)
C            f    - F10.7 solar radio flux (daily)
C            vKp  - Kp magnetic index (3-hour)
C            f5SW - indicator for Stratospheric Warming (SW) conditions
C                   =0 no SW, =0.5 minor SW, =1 major SW
C            f6WA - indicator for Winter Anomaly (WA) conditions
C                   =0 no WA, =0.5 weak WA, =1 strong WA
C Criteria for SW and WA indicators:
C      SW minor:  Temperature increase at the 30 hPa level by 10 deg.
C      SA major:  The same but by 20 degrees.
C         Temperature data for each year are published  
C         in Beilage zur Berliner Wetterkarte (K. Labitzke et al.).
C      WA weak:   An increase of the absorption in the 2-2.8 MHz  
C                 range at short A3 paths by 15 dB
C      WA strong: The same by 30 dB.
C 
C       Only for month 12 to 2 (winter).
C
C Output:      elg(7)  alog10 of electron density [cm-3] at h=60,65,
C                  70,75,80,85, and 90km
c-----------------------------------------------------------------------
c            
cor   dimension h(7),A0(7),A1(7),A2(7),A3(7),A4(7),A5(7),A6(7),elg(7)
      dimension A0(7),A1(7),A2(7),A3(7),A4(7),A5(7),A6(7),elg(7)
      data A0/1.0,1.2,1.4,1.5,1.6,1.7,3.0/
      data A1/0.6,0.8,1.1,1.2,1.3,1.4,1.0/
      data A2/0.,0.,0.08,0.12,0.05,0.2,0./
      data A3/0.,0.,0.,0.,0.,0.,1./
      data A4/0.,0.,-0.30,0.10,0.20,0.30,0.15/
      data A5/0.,-0.10,-0.20,-0.25,-0.30,-.30,0./
      data A6/0.,0.1,0.3,0.6,1.,1.,0.7/
        pi=3.14159265
         if(z.le.45) then
           f1z=1.
         else
           if(z.lt.90) then
             f1z=1.1892*(cos(z*pi/180))**0.5
           else
             f1z=0.
           endif
         endif
         f4S=1.
       if((it.ge.5).and.(it.le.9))then
         f4S=0.
         f5SW=0
         f6WA=0
       endif
       if((it.eq.3).or.(it.eq.4).or.(it.eq.10).or.(it.eq.11))then
         f4S=0.5
         f5SW=0
         f6WA=0
       endif
         f2Kp=vKp
         if(vKp.gt.2) f2Kp=2.
         f3F=(f-60.)/300.*f1z
         do 1 i=1,7
         elg(i)=A0(i)+A1(i)*f1z+A2(i)*f2Kp+A3(i)*f3F+A4(i)*f4S
     *         +A5(i)*f5SW+A6(i)*f6WA
   1     continue
         end
C
C
C
C                     
C************************************************************                   
C*************** EARTH MAGNETIC FIELD ***********************                   
C**************************************************************                 
C
C
      SUBROUTINE FIELDG(DLAT,DLONG,ALT,X,Y,Z,F,DIP,DEC,SMODIP)                  
C THIS IS A SPECIAL VERSION OF THE POGO 68/10 MAGNETIC FIELD                    
C LEGENDRE MODEL. TRANSFORMATION COEFF. G(144) VALID FOR 1973.                  
C INPUT: DLAT, DLONG=GEOGRAPHIC COORDINATES/DEG.(-90/90,0/360),                 
C        ALT=ALTITUDE/KM.                          
C OUTPUT: F TOTAL FIELD (GAUSS), Z DOWNWARD VERTICAL COMPONENT                  
C        X,Y COMPONENTS IN THE EQUATORIAL PLANE (X TO ZERO LONGITUDE).          
C        DIP INCLINATION ANGLE(DEGREE). SMODIP RAWER'S MODFIED DIP.             
C SHEIK,1977.         
      DIMENSION H(144),XI(3),G(144),FEL1(72),FEL2(72)
      COMMON/CONST/UMR,PI                           
      DATA FEL1/0.0, 0.1506723,0.0101742, -0.0286519, 0.0092606,                
     & -0.0130846, 0.0089594, -0.0136808,-0.0001508, -0.0093977,                
     & 0.0130650, 0.0020520, -0.0121956, -0.0023451, -0.0208555,                
     & 0.0068416,-0.0142659, -0.0093322, -0.0021364, -0.0078910,                
     & 0.0045586,  0.0128904, -0.0002951, -0.0237245,0.0289493,                 
     & 0.0074605, -0.0105741, -0.0005116, -0.0105732, -0.0058542,               
     &0.0033268, 0.0078164,0.0211234, 0.0099309, 0.0362792,                     
     &-0.0201070,-0.0046350,-0.0058722,0.0011147,-0.0013949,                    
     & -0.0108838,  0.0322263, -0.0147390,  0.0031247, 0.0111986,               
     & -0.0109394,0.0058112,  0.2739046, -0.0155682, -0.0253272,                
     &  0.0163782, 0.0205730,  0.0022081, 0.0112749,-0.0098427,                 
     & 0.0072705, 0.0195189, -0.0081132, -0.0071889, -0.0579970,                
     & -0.0856642, 0.1884260,-0.7391512, 0.1210288, -0.0241888,                 
     & -0.0052464, -0.0096312, -0.0044834, 0.0201764,  0.0258343,               
     &0.0083033,  0.0077187/                       
      DATA FEL2/0.0586055,0.0102236,-0.0396107,    
     & -0.0167860, -0.2019911, -0.5810815,0.0379916,  3.7508268,                
     & 1.8133030, -0.0564250, -0.0557352, 0.1335347, -0.0142641,                
     & -0.1024618,0.0970994, -0.0751830,-0.1274948, 0.0402073,                  
     &  0.0386290, 0.1883088,  0.1838960, -0.7848989,0.7591817,                 
     & -0.9302389,-0.8560960, 0.6633250, -4.6363869, -13.2599277,               
     & 0.1002136,  0.0855714,-0.0991981, -0.0765378,-0.0455264,                 
     &  0.1169326, -0.2604067, 0.1800076, -0.2223685, -0.6347679,               
     &0.5334222, -0.3459502,-0.1573697,  0.8589464, 1.7815990,                  
     &-6.3347645, -3.1513653, -9.9927750,13.3327637, -35.4897308,               
     &37.3466339, -0.5257398,  0.0571474, -0.5421217,  0.2404770,               
     & -0.1747774,-0.3433644, 0.4829708,0.3935944, 0.4885033,                   
     &  0.8488121, -0.7640999, -1.8884945, 3.2930784,-7.3497229,                
     & 0.1672821,-0.2306652, 10.5782146, 12.6031065, 8.6579742,                 
     & 215.5209961, -27.1419220,22.3405762,1108.6394043/                        
      K=0             
      DO 10 I=1,72    
      K=K+1           
      G(K)=FEL1(I)    
10    G(72+K)=FEL2(I)                              
      RLAT=DLAT*UMR   
      CT=SIN(RLAT)    
      ST=COS(RLAT)    
      NMAX=11         
      D=SQRT(40680925.0-272336.0*CT*CT)            
      RLONG=DLONG*UMR                              
      CP=COS(RLONG)   
      SP=SIN(RLONG)   
      ZZZ=(ALT+40408589.0/D)*CT/6371.2             
      RHO=(ALT+40680925.0/D)*ST/6371.2             
      XXX=RHO*CP      
      YYY=RHO*SP      
      RQ=1.0/(XXX*XXX+YYY*YYY+ZZZ*ZZZ)             
      XI(1)=XXX*RQ    
      XI(2)=YYY*RQ    
      XI(3)=ZZZ*RQ    
      IHMAX=NMAX*NMAX+1                            
      LAST=IHMAX+NMAX+NMAX                         
      IMAX=NMAX+NMAX-1                             
      DO 100 I=IHMAX,LAST                          
100   H(I)=G(I)       
      DO 200 K=1,3,2  
      I=IMAX          
      IH=IHMAX        
300   IL=IH-I         
      F1=2./(I-K+2.)  
      X1=XI(1)*F1     
      Y1=XI(2)*F1     
      Z1=XI(3)*(F1+F1)                             
      I=I-2           
      IF((I-1).LT.0) GOTO 400                      
      IF((I-1).EQ.0) GOTO 500                      
      DO 600 M=3,I,2  
      H(IL+M+1)=G(IL+M+1)+Z1*H(IH+M+1)+X1*(H(IH+M+3)-H(IH+M-1))-                
     &Y1*(H(IH+M+2)+H(IH+M-2))                     
      H(IL+M)=G(IL+M)+Z1*H(IH+M)+X1*(H(IH+M+2)-H(IH+M-2))+                      
     &Y1*(H(IH+M+3)+H(IH+M-1))                     
600   CONTINUE        
500   H(IL+2)=G(IL+2)+Z1*H(IH+2)+X1*H(IH+4)-Y1*(H(IH+3)+H(IH))                  
      H(IL+1)=G(IL+1)+Z1*H(IH+1)+Y1*H(IH+4)+X1*(H(IH+3)-H(IH))                  
400   H(IL)=G(IL)+Z1*H(IH)+2.0*(X1*H(IH+1)+Y1*H(IH+2))                          
700   IH=IL           
      IF(I.GE.K) GOTO 300                          
200   CONTINUE        
      S=0.5*H(1)+2.0*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))                         
      XT=(RQ+RQ)*SQRT(RQ)                          
      X=XT*(H(3)-S*XXX)                            
      Y=XT*(H(4)-S*YYY)                            
      Z=XT*(H(2)-S*ZZZ)                            
      F=SQRT(X*X+Y*Y+Z*Z)                          
      BRH0=Y*SP+X*CP  
      Y=Y*CP-X*SP     
      X=Z*ST-BRH0*CT  
      Z=-Z*CT-BRH0*ST 
        zdivf=z/f
        IF(ABS(zdivf).GT.1.) zdivf=SIGN(1.,zdivf)
      DIP=ASIN(zdivf)
        ydivs=y/sqrt(x*x+y*y)  
        IF(ABS(ydivs).GT.1.) ydivs=SIGN(1.,ydivs)
      DEC=ASIN(ydivs)
        dipdiv=DIP/SQRT(DIP*DIP+ST)
        IF(ABS(dipdiv).GT.1.) dipdiv=SIGN(1.,dipdiv)
      SMODIP=ASIN(dipdiv)
      DIP=DIP/UMR     
      DEC=DEC/UMR     
      SMODIP=SMODIP/UMR                            
      RETURN          
      END             
C
C
C************************************************************                   
C*********** INTERPOLATION AND REST ***************************                 
C**************************************************************                 
C
C
      SUBROUTINE REGFA1(X11,X22,FX11,FX22,EPS,FW,F,SCHALT,X) 
C REGULA-FALSI-PROCEDURE TO FIND X WITH F(X)-FW=0. X1,X2 ARE THE                
C STARTING VALUES. THE COMUTATION ENDS WHEN THE X-INTERVAL                      
C HAS BECOME LESS THEN EPS . IF SIGN(F(X1)-FW)= SIGN(F(X2)-FW)                  
C THEN SCHALT=.TRUE.  
      LOGICAL L1,LINKS,K,SCHALT                    
      SCHALT=.FALSE.
      EP=EPS  
      X1=X11          
      X2=X22          
      F1=FX11-FW     
      F2=FX22-FW     
      K=.FALSE.       
      NG=2       
      LFD=0  
      IF(F1*F2.LE.0.0) GOTO 200
        X=0.0           
        SCHALT=.TRUE.   
        RETURN
200   X=(X1*F2-X2*F1)/(F2-F1) 
      GOTO 400        
300     L1=LINKS        
        DX=(X2-X1)/NG
        IF(.NOT.LINKS) DX=DX*(NG-1)
        X=X1+DX
400   FX=F(X)-FW
      LFD=LFD+1
      IF(LFD.GT.20) THEN
        EP=EP*10.
        LFD=0
      ENDIF 
      LINKS=(F1*FX.GT.0.0)
      K=.NOT.K        
      IF(LINKS) THEN
        X1=X            
        F1=FX           
      ELSE
        X2=X 
        F2=FX 
      ENDIF   
      IF(ABS(X2-X1).LE.EP) GOTO 800               
      IF(K) GOTO 300  
      IF((LINKS.AND.(.NOT.L1)).OR.(.NOT.LINKS.AND.L1)) NG=2*NG                  
      GOTO 200        
800   RETURN          
      END             
C
C
C******************************************************************
C********** ZENITH ANGLE, DAY OF YEAR, TIME ***********************
C******************************************************************
C
C
        subroutine soco (ld,t,flat,Elon,height,
     &          DECLIN, ZENITH, SUNRSE, SUNSET)
c--------------------------------------------------------------------
c       s/r to calculate the solar declination, zenith angle, and
c       sunrise & sunset times  - based on Newbern Smith's algorithm
c       [leo mcnamara, 1-sep-86, last modified 16-jun-87]
c       {dieter bilitza, 30-oct-89, modified for IRI application}
c
c in:   ld      local day of year
c       t       local hour (decimal)
c       flat    northern latitude in degrees
c       elon    east longitude in degrees
c		height	height in km
c
c out:  declin      declination of the sun in degrees
c       zenith      zenith angle of the sun in degrees
c       sunrse      local time of sunrise in hours 
c       sunset      local time of sunset in hours 
c-------------------------------------------------------------------
c
        common/const/dtr,pi     /const1/humr,dumr
c amplitudes of Fourier coefficients  --  1955 epoch.................
        data    p1,p2,p3,p4,p6 /
     &  0.017203534,0.034407068,0.051610602,0.068814136,0.103221204 /
c
c s/r is formulated in terms of WEST longitude.......................
        wlon = 360. - Elon
c
c time of equinox for 1980...........................................
        td = ld + (t + Wlon/15.) / 24.
        te = td + 0.9369
c
c declination of the sun..............................................
        dcl = 23.256 * sin(p1*(te-82.242)) + 0.381 * sin(p2*(te-44.855))
     &      + 0.167 * sin(p3*(te-23.355)) - 0.013 * sin(p4*(te+11.97))
     &      + 0.011 * sin(p6*(te-10.41)) + 0.339137
        DECLIN = dcl
        dc = dcl * dtr
c
c the equation of time................................................
        tf = te - 0.5
        eqt = -7.38*sin(p1*(tf-4.)) - 9.87*sin(p2*(tf+9.))
     &      + 0.27*sin(p3*(tf-53.)) - 0.2*cos(p4*(tf-17.))
        et = eqt * dtr / 4.
c
        fa = flat * dtr
        phi = humr * ( t - 12.) + et
c
        a = sin(fa) * sin(dc)
        b = cos(fa) * cos(dc)
        	cosx = a + b * cos(phi)
        	if(abs(cosx).gt.1.) cosx=sign(1.,cosx)
        zenith = acos(cosx) / dtr
c
c calculate sunrise and sunset times --  at the ground...........
c see Explanatory Supplement to the Ephemeris (1961) pg 401......
c sunrise at height h metres is at...............................
		h=height*1000.
        chih = 90.83 + 0.0347 * sqrt(h)
c this includes corrections for horizontal refraction and........
c semi-diameter of the solar disk................................
        ch = cos(chih * dtr)
        cosphi = (ch -a ) / b
c if abs(secphi) > 1., sun does not rise/set.....................
c allow for sun never setting - high latitude summer.............
        secphi = 999999.
        if(cosphi.ne.0.) secphi = 1./cosphi
        sunset = 99.
        sunrse = 99.
        if(secphi.gt.-1.0.and.secphi.le.0.) return
c allow for sun never rising - high latitude winter..............
        sunset = -99.
        sunrse = -99.
        if(secphi.gt.0.0.and.secphi.lt.1.) return
c
        	cosx = cosphi
        	if(abs(cosx).gt.1.) cosx=sign(1.,cosx)
        phi = acos(cosx)
        et = et / humr
        phi = phi / humr
        sunrse = 12. - phi - et
        sunset = 12. + phi - et
        if(sunrse.lt.0.) sunrse = sunrse + 24.
        if(sunset.ge.24.) sunset = sunset - 24.
c
        return
        end
c
C
      FUNCTION HPOL(HOUR,TW,XNW,SA,SU,DSA,DSU)            
C-------------------------------------------------------
C PROCEDURE FOR SMOOTH TIME-INTERPOLATION USING EPSTEIN  
C STEP FUNCTION AT SUNRISE (SA) AND SUNSET (SU). THE 
C STEP-WIDTH FOR SUNRISE IS DSA AND FOR SUNSET DSU.
C TW,NW ARE THE DAY AND NIGHT VALUE OF THE PARAMETER TO 
C BE INTERPOLATED. SA AND SU ARE TIME OF SUNRIES AND 
C SUNSET IN DECIMAL HOURS.
C BILITZA----------------------------------------- 1979.
        IF(ABS(SU).GT.25.) THEN
                IF(SU.GT.0.0) THEN
                        HPOL=TW
                ELSE
                        HPOL=XNW
                ENDIF
                RETURN
        ENDIF
      HPOL=XNW+(TW-XNW)*EPST(HOUR,DSA,SA)+
     &  (XNW-TW)*EPST(HOUR,DSU,SU) 
      RETURN          
      END       
C      
C
        SUBROUTINE MODA(IN,IYEAR,MONTH,IDAY,IDOY,NRDAYMO)
C-------------------------------------------------------------------
C CALCULATES DAY OF YEAR (IDOY, ddd) FROM YEAR (IYEAR, yy or yyyy), 
C MONTH (MONTH, mm) AND DAY OF MONTH (IDAY, dd) IF IN=0, OR MONTH 
C AND DAY FROM YEAR AND DAY OF YEAR IF IN=1. NRDAYMO is an output 
C parameter providing the number of days in the specific month.
C-------------------------------------------------------------------
        DIMENSION       MM(12)
        DATA            MM/31,28,31,30,31,30,31,31,30,31,30,31/

        IMO=0
        MOBE=0
c
c  leap year rule: years evenly divisible by 4 are leap years, except
c  years also evenly divisible by 100 are not leap years, except years 
c  also evenly divisible by 400 are leap years. The year 2000 therefore 
C  is a leap year. The 100 and 400 year exception rule
c     if((iyear/4*4.eq.iyear).and.(iyear/100*100.ne.iyear)) mm(2)=29
c  will become important again in the year 2100 which is not a leap 
C  year.
c
        mm(2)=28
        if(iyear/4*4.eq.iyear) mm(2)=29

        IF(IN.GT.0) GOTO 5
                mosum=0
                if(month.gt.1) then
                        do 1234 i=1,month-1 
1234                            mosum=mosum+mm(i)
                        endif
                idoy=mosum+iday
                nrdaymo=mm(month)
                RETURN

5       IMO=IMO+1
                IF(IMO.GT.12) GOTO 55
                MOOLD=MOBE
                nrdaymo=mm(imo)
                MOBE=MOBE+nrdaymo
                IF(MOBE.LT.IDOY) GOTO 5
55              MONTH=IMO
                IDAY=IDOY-MOOLD
        RETURN
        END             
c
c
        subroutine ut_lt(mode,ut,slt,glong,iyyy,ddd)
c -----------------------------------------------------------------
c Converts Universal Time UT (decimal hours) into Solar Local Time
c SLT (decimal hours) for given date (iyyy is year, e.g. 1995; ddd
c is day of year, e.g. 1 for Jan 1) and geodatic longitude in degrees.
C For mode=0 UT->LT and for mode=1 LT->UT
c Please NOTE that iyyy and ddd are input as well as output parameters
c since the determined LT may be for a day before or after the UT day.
c ------------------------------------------------- bilitza nov 95
        integer         ddd,dddend

        xlong=glong
        if(glong.gt.180) xlong=glong-360
        if(mode.ne.0) goto 1
c
c UT ---> LT
c
        SLT=UT+xlong/15.
        if((SLT.ge.0.).and.(SLT.le.24.)) goto 2
        if(SLT.gt.24.) goto 3
                SLT=SLT+24.
                ddd=ddd-1
                if(ddd.lt.1.) then
                        iyyy=iyyy-1
                        ddd=365
c
c leap year if evenly divisible by 4 and not by 100, except if evenly
c divisible by 400. Thus 2000 will be a leap year.
c
                        if(iyyy/4*4.eq.iyyy) ddd=366
                        endif
                goto 2
3               SLT=SLT-24.
                ddd=ddd+1
                dddend=365
                if(iyyy/4*4.eq.iyyy) dddend=366
                if(ddd.gt.dddend) then
                        iyyy=iyyy+1
                        ddd=1
                        endif
                goto 2
c
c LT ---> UT
c
1       UT=SLT-xlong/15.
        if((UT.ge.0.).and.(UT.le.24.)) goto 2
        if(UT.gt.24.) goto 5
                UT=UT+24.
                ddd=ddd-1
                if(ddd.lt.1.) then
                        iyyy=iyyy-1
                        ddd=365
                        if(iyyy/4*4.eq.iyyy) ddd=366
                        endif
                goto 2
5               UT=UT-24.
                ddd=ddd+1
                dddend=365
                if(iyyy/4*4.eq.iyyy) dddend=366
                if(ddd.gt.dddend) then
                        iyyy=iyyy+1
                        ddd=1
                        endif
2       return
        end
C
C

      SUBROUTINE SUN (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C-----------------------------------------------------------------------------
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C  LAST MODIFICATION:  MARCH 31, 2003 (ONLY SOME NOTATION CHANGES)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION DJ,FDAY
      COMMON /CONST/UMR,PI
C
      IF(IYEAR.LT.1901.OR.IYEAR.GT.2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.
      VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
      GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)*UMR
      G=DMOD(358.475845+0.985600267*DJ,360.D0)*UMR
      SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))*UMR
      IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
      IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
      OBLIQ=(23.45229-0.0130125*T)*UMR
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924E-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
C   THE ORBITAL MOTION OF THE EARTH
C
      SIN1=SOB*SIN(SLP)
      COS1=SQRT(1.-SIN1**2)
      SC=SIN1/COS1
      SDEC=ATAN(SC)
      SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COS1)
      RETURN
      END
C      
C
C *********************************************************************
C ************************ EPSTEIN FUNCTIONS **************************
C *********************************************************************
C REF:  H. G. BOOKER, J. ATMOS. TERR. PHYS. 39, 619-623, 1977
C       K. RAWER, ADV. SPACE RES. 4, #1, 11-15, 1984
C *********************************************************************
C
C
        REAL FUNCTION  RLAY ( X, XM, SC, HX )
C -------------------------------------------------------- RAWER  LAYER
        Y1  = EPTR ( X , SC, HX )
        Y1M = EPTR ( XM, SC, HX )
        Y2M = EPST ( XM, SC, HX )
        RLAY = Y1 - Y1M - ( X - XM ) * Y2M / SC
        RETURN
        END
C
C
        REAL FUNCTION D1LAY ( X, XM, SC, HX )
C ------------------------------------------------------------ dLAY/dX
        D1LAY = ( EPST(X,SC,HX) - EPST(XM,SC,HX) ) /  SC
        RETURN
        END
C
C
        REAL FUNCTION D2LAY ( X, XM, SC, HX )
C ---------------------------------------------------------- d2LAY/dX2
        D2LAY = EPLA(X,SC,HX) /  (SC * SC)
        RETURN
        END
C
C
        REAL FUNCTION EPTR ( X, SC, HX )
C --------------------------------------------------------- TRANSITION
        COMMON/ARGEXP/ARGMAX
        D1 = ( X - HX ) / SC
        IF (ABS(D1).LT.ARGMAX) GOTO 1
        IF (D1.GT.0.0) THEN
          EPTR = D1
        ELSE
          EPTR = 0.0
        ENDIF
        RETURN
1       EPTR = ALOG ( 1. + EXP( D1 ))
        RETURN
        END
C
C
        REAL FUNCTION EPST ( X, SC, HX )
C -------------------------------------------------------------- STEP
        COMMON/ARGEXP/ARGMAX
        D1 = ( X - HX ) / SC
        IF (ABS(D1).LT.ARGMAX) GOTO 1
        IF (D1.GT.0.0) THEN
          EPST = 1.
        ELSE
          EPST = 0.
        ENDIF
        RETURN
1       EPST = 1. / ( 1. + EXP( -D1 ))
        RETURN
        END
C
C
        REAL FUNCTION EPSTEP ( Y2, Y1, SC, HX, X)
C---------------------------------------------- STEP FROM Y1 TO Y2      
        EPSTEP = Y1 + ( Y2 - Y1 ) * EPST ( X, SC, HX)
        RETURN
        END
C
C
        REAL FUNCTION EPLA ( X, SC, HX )
C ------------------------------------------------------------ PEAK 
        COMMON/ARGEXP/ARGMAX
        D1 = ( X - HX ) / SC
        IF (ABS(D1).LT.ARGMAX) GOTO 1
                EPLA = 0
                RETURN  
1       D0 = EXP ( D1 )
        D2 = 1. + D0
        EPLA = D0 / ( D2 * D2 )
        RETURN
        END
c
c
        FUNCTION XE2TO5(H,HMF2,NL,HX,SC,AMP)
C----------------------------------------------------------------------
C NORMALIZED ELECTRON DENSITY (N/NMF2) FOR THE MIDDLE IONOSPHERE FROM 
C HME TO HMF2 USING LAY-FUNCTIONS.
C----------------------------------------------------------------------
        DIMENSION       HX(NL),SC(NL),AMP(NL)
        SUM = 1.0
        DO 1 I=1,NL
           YLAY = AMP(I) * RLAY( H, HMF2, SC(I), HX(I) )
           zlay=10.**ylay
1          sum=sum*zlay
        XE2TO5 = sum
        RETURN
        END
C
C
        REAL FUNCTION XEN(H,HMF2,XNMF2,HME,NL,HX,SC,AMP)
C----------------------------------------------------------------------
C ELECTRON DENSITY WITH NEW MIDDLE IONOSPHERE
C----------------------------------------------------------------------
        DIMENSION       HX(NL),SC(NL),AMP(NL)
C
        IF(H.LT.HMF2) GOTO 100
                XEN = XE1(H)
                RETURN
100     IF(H.LT.HME) GOTO 200
                XEN = XNMF2 * XE2TO5(H,HMF2,NL,HX,SC,AMP)
                RETURN
200     XEN = XE6(H)
        RETURN
        END
C
C
        SUBROUTINE ROGUL(IDAY,XHI,SX,GRO)
C --------------------------------------------------------------------- 
C   CALCULATES RATIO H0.5/HMF2 FOR HALF-DENSITY POINT (NE(H0.5)=0.5*
C   NMF2) T. GULYAEVA, ADVANCES IN SPACE RESEARCH 7, #6, 39-48, 1987.
C
C       INPUT:  IDAY    DAY OF YEAR
C               XHI     SOLAR ZENITH ANGLE [DEGREE]
C       
C       OUTPUT: GRO     RATIO OF HALF DENSITY HEIGHT TO F PEAK HEIGHT
C               SX      SMOOTHLY VARYING SEASON PARAMTER (SX=1 FOR 
C                       DAY=1; SX=3 FOR DAY=180; SX=2 FOR EQUINOX)
C ---------------------------------------------------------------------
C
        common  /const1/humr,dumr
        SX = 2. - COS ( IDAY * dumr )
        XS = ( XHI - 20. * SX) / 15.
        GRO = 0.8 - 0.2 / ( 1. + EXP(XS) )
c same as gro=0.6+0.2/(1+exp(-xs))
        RETURN
        END
C
C
        SUBROUTINE LNGLSN ( N, A, B, AUS)
C --------------------------------------------------------------------
C SOLVES QUADRATIC SYSTEM OF LINEAR EQUATIONS:
C
C       INPUT:  N       NUMBER OF EQUATIONS (= NUMBER OF UNKNOWNS)
C               A(N,N)  MATRIX (LEFT SIDE OF SYSTEM OF EQUATIONS)
C               B(N)    VECTOR (RIGHT SIDE OF SYSTEM)
C
C       OUTPUT: AUS     =.TRUE.   NO SOLUTION FOUND
C                       =.FALSE.  SOLUTION IS IN  A(N,J) FOR J=1,N
C --------------------------------------------------------------------
C
        DIMENSION       A(5,5), B(5), AZV(10)
        LOGICAL         AUS
C
        NN = N - 1
        AUS = .FALSE.
        DO 1 K=1,N-1
                IMAX = K
                L    = K
                IZG  = 0
                AMAX  = ABS( A(K,K) )
110             L = L + 1
                IF (L.GT.N) GOTO 111
                HSP = ABS( A(L,K) )
                IF (HSP.LT.1.E-8) IZG = IZG + 1
                IF (HSP.LE.AMAX) GOTO 110
111             IF (ABS(AMAX).GE.1.E-10) GOTO 133
                        AUS = .TRUE.
                        RETURN
133             IF (IMAX.EQ.K) GOTO 112
                DO 2 L=K,N
                        AZV(L+1)  = A(IMAX,L)
                        A(IMAX,L) = A(K,L)
2                       A(K,L)    = AZV(L+1)
                AZV(1)  = B(IMAX)
                B(IMAX) = B(K)
                B(K)    = AZV(1)
112             IF (IZG.EQ.(N-K)) GOTO 1
                AMAX = 1. / A(K,K)
                AZV(1) = B(K) * AMAX
                DO 3 M=K+1,N
3                       AZV(M+1) = A(K,M) * AMAX
                DO 4 L=K+1,N
                        AMAX = A(L,K)
                        IF (ABS(AMAX).LT.1.E-8) GOTO 4
                        A(L,K) = 0.0
                        B(L) = B(L) - AZV(1) * AMAX
                        DO 5 M=K+1,N
5                               A(L,M) = A(L,M) - AMAX * AZV(M+1)
4               CONTINUE
1       CONTINUE
        DO 6 K=N,1,-1
                AMAX = 0.0
                IF (K.LT.N) THEN
                        DO 7 L=K+1,N
7                               AMAX = AMAX + A(K,L) * A(N,L)
                        ENDIF
                IF (ABS(A(K,K)).LT.1.E-6) THEN
                        A(N,K) = 0.0
                ELSE
                        A(N,K) = ( B(K) - AMAX ) / A(K,K)
                ENDIF
6       CONTINUE
        RETURN
        END
C
C
        SUBROUTINE LSKNM ( N, M, M0, M1, HM, SC, HX, W,X,Y,VAR,SING)
C --------------------------------------------------------------------
C   DETERMINES LAY-FUNCTIONS AMPLITUDES FOR A NUMBER OF CONSTRAINTS:
C
C       INPUT:  N       NUMBER OF AMPLITUDES ( LAY-FUNCTIONS)
C               M       NUMBER OF CONSTRAINTS
C               M0      NUMBER OF POINT CONSTRAINTS
C               M1      NUMBER OF FIRST DERIVATIVE CONSTRAINTS
C               HM      F PEAK ALTITUDE  [KM]
C               SC(N)   SCALE PARAMETERS FOR LAY-FUNCTIONS  [KM]
C               HX(N)   HEIGHT PARAMETERS FOR LAY-FUNCTIONS  [KM]
C               W(M)    WEIGHT OF CONSTRAINTS
C               X(M)    ALTITUDES FOR CONSTRAINTS  [KM]
C               Y(M)    LOG(DENSITY/NMF2) FOR CONSTRAINTS
C
C       OUTPUT: VAR(M)  AMPLITUDES
C               SING    =.TRUE.   NO SOLUTION
C ---------------------------------------------------------------------
C
        LOGICAL         SING
        DIMENSION       VAR(N), HX(N), SC(N), W(M), X(M), Y(M),
     &                  BLI(5), ALI(5,5), XLI(5,10)
C
        M01=M0+M1
        SCM=0
        DO 1 J=1,5
                BLI(J) = 0.
                DO 1 I=1,5
1                       ALI(J,I) = 0. 
        DO 2 I=1,N
                DO 3 K=1,M0
3                       XLI(I,K) = RLAY( X(K), HM, SC(I), HX(I) )
                DO 4 K=M0+1,M01
4                       XLI(I,K) = D1LAY( X(K), HM, SC(I), HX(I) )
                DO 5 K=M01+1,M
5                       XLI(I,K) = D2LAY( X(K), HM, SC(I), HX(I) )
2       CONTINUE
                DO 7 J=1,N
                DO 6 K=1,M
                        BLI(J) = BLI(J) + W(K) * Y(K) * XLI(J,K)
                        DO 6 I=1,N
6                               ALI(J,I) = ALI(J,I) + W(K) * XLI(I,K) 
     &                                  * XLI(J,K)
7       CONTINUE
        CALL LNGLSN( N, ALI, BLI, SING )
        IF (.NOT.SING) THEN
                DO 8 I=1,N
8                       VAR(I) = ALI(N,I)
                ENDIF
        RETURN
        END
C
C
        SUBROUTINE INILAY(NIGHT,F1REG,XNMF2,XNMF1,XNME,VNE,HMF2,HMF1, 
     &                          HME,HV1,HV2,HHALF,HXL,SCL,AMP,IQUAL)
C-------------------------------------------------------------------
C CALCULATES AMPLITUDES FOR LAY FUNCTIONS
C D. BILITZA, DECEMBER 1988
C
C INPUT:        NIGHT   LOGICAL VARIABLE FOR DAY/NIGHT DISTINCTION
C               F1REG   LOGICAL VARIABLE FOR F1 OCCURRENCE
C               XNMF2   F2 PEAK ELECTRON DENSITY [M-3]
C               XNMF1   F1 PEAK ELECTRON DENSITY [M-3]
C               XNME    E  PEAK ELECTRON DENSITY [M-3]
C               VNE     ELECTRON DENSITY AT VALLEY BASE [M-3]
C               HMF2    F2 PEAK ALTITUDE [KM]
C               HMF1    F1 PEAK ALTITUDE [KM]
C               HME     E  PEAK ALTITUDE [KM]
C               HV1     ALTITUDE OF VALLEY TOP [KM]
C               HV2     ALTITUDE OF VALLEY BASE [KM]
C               HHALF   ALTITUDE OF HALF-F2-PEAK-DENSITY [KM]
C
C OUTPUT:       HXL(4)  HEIGHT PARAMETERS FOR LAY FUNCTIONS [KM] 
C               SCL(4)  SCALE PARAMETERS FOR LAY FUNCTIONS [KM]
C               AMP(4)  AMPLITUDES FOR LAY FUNCTIONS
C               IQUAL   =0 ok, =1 ok using second choice for HXL(1)
C                       =2 NO SOLUTION
C---------------------------------------------------------------  
        DIMENSION       XX(8),YY(8),WW(8),AMP(4),HXL(4),SCL(4)
        LOGICAL         SSIN,NIGHT,F1REG
c
c constants --------------------------------------------------------
                NUMLAY=4
                NC1 = 2
                ALG102=ALOG10(2.)
c
c constraints: xx == height     yy == log(Ne/NmF2)    ww == weights
c -----------------------------------------------------------------
                ALOGF = ALOG10(XNMF2)
                ALOGEF = ALOG10(XNME) - ALOGF
                XHALF=XNMF2/2.
                XX(1) = HHALF
                XX(2) = HV1
                XX(3) = HV2
                XX(4) = HME
                XX(5) = HME - ( HV2 - HME )
                YY(1) = -ALG102
                YY(2) = ALOGEF
                YY(3) = ALOG10(VNE) - ALOGF
                YY(4) = ALOGEF
                YY(5) = YY(3)
                YY(7) = 0.0
                WW(2) = 1.
                WW(3) = 2.
                WW(4) = 5.
c
c geometric paramters for LAY -------------------------------------
c difference to earlier version:  HXL(3) = HV2 + SCL(3)
c
                SCL0 = 0.7 * ( 0.216 * ( HMF2 - HHALF ) + 56.8 )
                SCL(1) = 0.8 * SCL0
                SCL(2) = 10.
                SCL(3) = 9.
                SCL(4) = 6.
                HXL(3) = HV2
                HFFF=HHALF
                XFFF=XHALF
c
C DAY CONDITION--------------------------------------------------
c earlier tested:       HXL(2) = HMF1 + SCL(2)
c 
            IF(NIGHT) GOTO 7711
                NUMCON = 8
                HXL(1) = 0.9 * HMF2
                  HXL1T  = HHALF
                HXL(2) = HMF1
                HXL(4) = HME - SCL(4)
                XX(6) = HMF1
                XX(7) = HV2
                XX(8) = HME
                YY(8) = 0.0
                WW(5) = 1.
                WW(7) = 50.
                WW(8) = 500.
c without F-region ----------------------------------------------
                IF(F1REG) GOTO 100
                        HXL(2)=(HMF2+HHALF)/2.
                        YY(6) = 0.
                        WW(6) = 0.
                        WW(1) = 1.
                        GOTO 7722
c with F-region --------------------------------------------
100             YY(6) = ALOG10(XNMF1) - ALOGF
                WW(6) = 3.
                IF((XNMF1-XHALF)*(HMF1-HHALF).LT.0.0) THEN
                  WW(1)=0.5
                ELSE
                  ZET = YY(1) - YY(6)
                  WW(1) = EPST( ZET, 0.1, 0.15)
                ENDIF
                IF(HHALF.GT.HMF1) THEN
                  HFFF=HMF1
                  XFFF=XNMF1
                ELSE
                  HFFF=HHALF
                  XFFF=XHALF
                ENDIF
                GOTO 7722
c
C NIGHT CONDITION---------------------------------------------------
c different HXL,SCL values were tested including: 
c       SCL(1) = HMF2 * 0.15 - 27.1     HXL(2) = 200.   
c       HXL(2) = HMF1 + SCL(2)          HXL(3) = 140.
c       SCL(3) = 5.                     HXL(4) = HME + SCL(4)
c       HXL(4) = 105.                   
c
7711            NUMCON = 7
                HXL(1) = HHALF
                  HXL1T  = 0.4 * HMF2 + 30.
                HXL(2) = ( HMF2 + HV1 ) / 2.
                HXL(4) = HME
                XX(6) = HV2
                XX(7) = HME
                YY(6) = 0.0
                WW(1) = 1.
                WW(3) = 3.
                WW(5) = 0.5
                WW(6) = 50.
                WW(7) = 500.
                HFFF=HHALF
                XFFF=XHALF
c
C are valley-top and bottomside point compatible ? -------------
C
7722    IF((HV1-HFFF)*(XNME-XFFF).LT.0.0) WW(2)=0.5
        IF(HV1.LE.HV2+5.0) WW(2)=0.5
c
C DETERMINE AMPLITUDES-----------------------------------------
C
            NC0=NUMCON-NC1
            IQUAL=0
2299        CALL LSKNM(NUMLAY,NUMCON,NC0,NC1,HMF2,SCL,HXL,WW,XX,YY,
     &          AMP,SSIN)
                IF(IQUAL.gt.0) GOTO 1937
            IF((ABS(AMP(1)).GT.10.0).OR.(SSIN)) THEN
                IQUAL=1
                HXL(1)=HXL1T
                GOTO 2299
                ENDIF
1937        IF(SSIN) IQUAL=2
            RETURN
            END
c
c
           subroutine read_ig_rz 
c----------------------------------------------------------------
c Reads the Rz12 and IG12 indices file IG_RZ.DAT from I/O UNIT=12 
c and stores the indices in COMMON:
c		common/igrz/aig,arziyst,iyed   with aig(806),arz(806),
c											start year (iyst)
c                                           end year (iyed)
c 
c The indices file IG_RZ.DAT is structured as follows (values are 
c separated by comma): 
c   day, month, year of the last update of this file,
c   a blank line
c   start month, start year, end month, end year,
c   a blank line
c   the IG index for December of start year minus 1 (this value is 
c		needed for interpolating from 1st to 15th of first year)
c   the 12 IG indices (13-months running mean) for start year, 
c   the 12 IG indices for the second year 
c       .. and so on until the last year,
c   the 12 IG indices for the last year 
c   the IG index for January of end year plus 1 (needed for interpolation)
c   a blank line
c   the Rz index for December of start year minus 1 
c   the 12 Rz indices (13-months running mean) for the start year,
c   the 12 Rz indices for the second year 
c       .. and so on until the last year.
c   the 12 Rz indices for the last year 
c   the Rz index for January of end year plus 1 
c 
c A negative Rz index means that the given index is the 13-months-
C running mean of the solar radio flux (F10.7). The close correlation 
C between (Rz)12 and (F10.7)12 is used to compute the (Rz)12 indices.
c
c An IG index of -111 indicates that no IG values are available for the
c time period. In this case a correlation function between (IG)12 and 
C (Rz)12 is used to obtain (IG)12.
c
c The computation of the 13-month-running mean for month M requires the
c indices for the six months preceeding M and the six months following 
C M (month: M-6, ..., M+6). To calculate the current running mean one 
C therefore requires predictions of the indix for the next six months. 
C Starting from six months before the UPDATE DATE (listed at the top of 
c the file) and onward the indices are therefore based on indices 
c predictions.
c----------------------------------------------------------------

           integer	iyst,iyend,iymst,iupd,iupm,iupy,imst,imend
           real		aig(806),arz(806)
           
           common /igrz/aig,arz,iymst,iymend

           open(unit=12,file='ig_rz.dat',FORM='FORMATTED',status='old')

c-web- special for web version
c            open(unit=12,file=
c     *         '/var/www/omniweb/cgi/vitmo/IRI/ig_rz.dat',
c     *         FORM='FORMATTED',status='old')

c Read the update date, the start date and the end date (mm,yyyy), and
c get number of data points to read.

            read(12,*) iupd,iupm,iupy
            read(12,*) imst,iyst,imend,iyend
            iymst=iyst*100+imst
            iymend=iyend*100+imend
     
c inum_vals= 12-imst+1+(iyend-iyst-1)*12 +imend + 2
c 1st year \ full years       \last y\ before & after

            inum_vals= 3-imst+(iyend-iyst)*12 +imend

c read all the IG12 (ionoindx) and Rz12 (indrz) values
            read(12,*) (aig(i),i=1,inum_vals)
            read(12,*) (arz(i),i=1,inum_vals)

c use scale factor 0.7 for new sunspot number starting from Jan 2014
c and starting with ig_rz file for Oct 2016.
            if(iupy*100+iupm.gt.201609) then
            	inum_chan= 3-imst+(2014-iyst)*12 
            	do 1 jj=inum_chan,inum_vals
                	arz(jj)=arz(jj)*0.7
c                ggg=aig(jj)
c                if(rrr.lt.0.0) then
c                    covr=abs(rrr)
c                    rrr=33.52*sqrt(covr+85.12)-408.99
c                    if(rrr.lt.0.0) rrr=0.0
c                    endif
c                if(ggg.le.-90.) then
c                    zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
c                    if(zi.gt.274.0) zi=274.0
c                    ggg=zi
c                    endif
c                arz(jj)=rrr
c                aig(jj)=ggg    
1               	continue
				endif
            close(unit=12)

            return
            end
c
c
           subroutine tcon(yr,mm,day,idn,rz,ig,rsn,nmonth)
c----------------------------------------------------------------
c input:        yr,mm,day       year(yyyy),month(mm),day(dd)
c               idn             day of year(ddd)
c output:       rz(3)           12-month-smoothed solar sunspot number
c               ig(3)           12-month-smoothed IG index
c               rsn             interpolation parameter
c               nmonth          previous or following month depending
c                               on day
c
c Uses read_ig_rz and common/igrz/ to get indices 
c 
c rz(1) & ig(1) contain the indices for the month mm and rz(2) & ig(2)
c for the previous month (if day less than 15) or for the following
c month (if day greater than 15). These indices are for the mid of the 
c month. The indices for the given day are obtained by linear 
c interpolation and are stored in rz(3) and ig(3).
c----------------------------------------------------------------

           integer	yr,mm,day,iyst,iyend,iymst
           integer	imst,iymend
           real		ionoindx(806),indrz(806)
           real		ig(3),rz(3)
           logical	mess
           
           common 	/iounit/konsol,mess  
           common	/igrz/ionoindx,indrz,iymst,iymend

        iytmp=yr*100+mm
        if (iytmp.lt.iymst.or.iytmp.gt.iymend) then
               if(mess) write(konsol,8000) iytmp,iymst,iymend
8000           format(1x,I10,'** OUT OF RANGE **'/,5x,
     &  'The file IG_RZ.DAT which contains the indices Rz12',
     &  ' and IG12'/5x,'currently only covers the time period',
     &  ' (yymm) : ',I6,'-',I6)
               nmonth=-1
               return
               endif

		iyst=iymst/100
		imst=iymst-iyst*100
c       num=12-imst+1+(yr-iyst-1)*12+mm+1
        num=2-imst+(yr-iyst)*12+mm

        rz(1)=indrz(num)
        ig(1)=ionoindx(num)
        midm=15
        if(mm.eq.2) midm=14
        call MODA(0,yr,mm,midm,idd1,nrdaym)

        if(day.lt.midm) goto 1926
c day is at or after mid of month
                imm2=mm+1
                if(imm2.gt.12) then
                        imm2=1
                        iyy2=yr+1
                        idd2=380            ! =365+15 mid-January
c               if((yr/4*4.eq.yr).and.(yr/100*100.ne.yr)) idd2=381
                        if(yr/4*4.eq.yr) idd2=381
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                rz(2)=indrz(num+1)
                ig(2)=ionoindx(num+1)
                rsn=(idn-idd1)*1./(idd2-idd1)                
                rz(3)=rz(1)+(rz(2)-rz(1))*rsn
                ig(3)=ig(1)+(ig(2)-ig(1))*rsn
                goto 1927
1926            imm2=mm-1
                if(imm2.lt.1) then
                        imm2=12
                        idd2=-16
                        iyy2=yr-1
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                rz(2)=indrz(num-1)
                ig(2)=ionoindx(num-1)
                rsn=(idn-idd2)*1./(idd1-idd2)
                rz(3)=rz(2)+(rz(1)-rz(2))*rsn
                ig(3)=ig(2)+(ig(1)-ig(2))*rsn

1927    nmonth=imm2
            return
            end
C
C
		subroutine readapf107
C-------------------------------------------------------------------------
c Reads APF107.DAT file (on UNIT=13) and stores contents in COMMON block:
C 	COMMON/AAP,AF107,N/ with  AAP(23000,9) and AF107(23000,3)
C		AAP(*,1)	3-hour Ap indices for the UT interval )0-3)
C		AAP(*,2)	3-hour Ap indices for the UT interval )3-6)
C          ....                       ....
C		AAP(*,8)	3-hour Ap indices for the UT interval )21-6)
C		AAP(*,9)	daily Ap
C		AF107(*,1)	F10.7 radio flux for the day
C		AF107(*,2)	81-day average of F10.7 radio flux 
C		AF107(*,3)	365-day average of F10.7
C       N           total number of records
c
c APF107.DAT is structured as follows:
c 		JY(I3),JMN(I3),JD(I3)	year, month, day 
c		IIAP(8)	(8I3)			3-hour Ap indices for the UT intervals 
c								(0-3(,(3-6(,(6-9(, .., (18-21(,(21-24(
c		IAPD (I3)				daily Ap
c		IR (I3)					sunspot number for the day (empty)
c		F107 (F5.1)				F10.7 radio flux for the day
c		F107_81 (F5.1)			81-day average of F10.7 radio flux 
c       F107_365 (F5.1)         365-day average of F10.7 centered on 
c                               the date of interest. At start and end  
c								of index file it takes all available  
c                               indices, e.g. for the first date the 
c                               average is only over 40 F10.7 values  
c                               and over 41 values on the 2nd date.  
c
c If date is outside the range of the Ap indices file then IAP(1)=-5  
C-------------------------------------------------------------------------
C
        INTEGER		aap(23000,9),iiap(8)
        DIMENSION 	af107(23000,3)
        COMMON		/apfa/aap,af107,n

        Open(13,FILE='apf107.dat',FORM='FORMATTED',STATUS='OLD')
c-web-sepcial vfor web version
c      OPEN(13,FILE='/var/www/omniweb/cgi/vitmo/IRI/apf107.dat',
c     *    FORM='FORMATTED',STATUS='OLD')

		i=1
1       READ(13,10,END=21) JY,JMN,JD,iiap,iapda,IR,F107D,F107_81,
     *        F107_365
10      FORMAT(3I3,9I3,I3,3F5.1)

c 		adate(i)=jy*10000+jmn*100+jd 
		do j=1,8 
             aap(i,j)=iiap(j)
             enddo
		aap(i,9)=iapda
c		irza(i)=ir
 		if(F107_81.lt.-4.) F107_81=F107D
 		if(F107_365.lt.-4.) F107_365=F107D
		af107(i,1)=f107d 
        af107(i,2)=f107_81 
        af107(i,3)=f107_365
		i=i+1 
        goto 1

21		n=i-1
 
        CLOSE(13)
		return
		end
C
C

        SUBROUTINE APF(ISDATE,HOUR,IAP)
c-----------------------------------------------------------------------
c Finds 3-hourly Ap indices for IRI-STORM model
c    INPUTS: 	ISDATE		Array-index from APF_ONLY
c			    HOUR		UT in decimal hours
c    OUTPUT:    IAP(1:13)	3-hourly Ap index
c							IAP(13) Ap index for current UT
c							IAP(1) AP index for UT-39 hours.
c
c Gets indices from COMMON/APFA/
c
c If date is outside the range of the Ap indices file than IAP(1)=-5  
c-----------------------------------------------------------------------

        INTEGER		aap(23000,9),iiap(8),iap(13)
        DIMENSION 	af107(23000,3)
        LOGICAL 	mess
        COMMON 		/iounit/konsol,mess	/apfa/aap,af107,nf107
       
        do i=1,8
              iap(i)=-1
              enddo

	   IS = ISDATE

        ihour=int(hour/3.)+1
        if(ihour.gt.8) ihour=8

        if(is*8+ihour.lt.13) goto 21   ! less then 13 indices available	

        j1=13-ihour
        do i=1,ihour
           iapi=aap(is,i)
           if(iapi.lt.-2) goto 21
           iap(j1+i)=iapi
           enddo

        if(ihour.gt.4) then
           do i=1,j1
              iapi=aap(is-1,8-j1+i)
              if(iapi.lt.-2) goto 21
              iap(i)=iapi
              enddo
        else           
           j2=5-ihour
           do i=1,8
              iapi=aap(is-1,i)
              if(iapi.lt.-2) goto 21
              iap(j2+i)=iapi
              enddo
           do i=1,j2
              iapi=aap(is-2,8-j2+i)
              if(iapi.lt.-2) goto 21
              iap(i)=iapi
              enddo
        endif         
        goto 20
        
21      if(mess) write(konsol,100)
100     format(1X,'One of the ap indeces is negative.',
     &     ' STORM model is turned off.')
        IAP(1)=-5
      
20    RETURN
      END
C
C
        SUBROUTINE APFMSIS(ISDATE,HOUR,IAPO)
c-----------------------------------------------------------------------
c Finds 3-hourly Ap indices for NRLMSIS00 model 
c    INPUTS: 	ISDATE		Array-index from APF_ONLY
c				HOUR		UT in decimal hours
c    OUTPUT:   	IAPO(1:7)	3-hourly Ap index
C  
C IAPO(1) DAILY AP
C IAPO(2) 3-HR AP INDEX FOR CURRENT TIME			   	
C IAPO(3) 3-HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME	
C IAPO(4) 3-HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME	
C IAPO(5) 3-HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME  
C IAPO(6) AVERAGE OF EIGHT 3-HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C         TO CURRENT TIME
C IAPO(7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C         TO CURRENT TIME
c
c The 3-hour UT intervals during the day are: (0-3),)3-6),)6-9),)9-12),
c )12-15),)15-18),)18-21),)21-24(.
c 
c If date is outside the range of the Ap indices file then IAPO(2)=-5  
c-----------------------------------------------------------------------
c
		REAL 		IAPO
        INTEGER		aap(23000,9),iiap(8)
        DIMENSION 	af107(23000,3),iap(20),lm(12),iapo(7)
        LOGICAL  	mess

        COMMON 		/iounit/konsol,mess	/apfa/aap,af107,nf107

        IS=ISDATE

        ihour=int(hour/3.)+1
        if(ihour.gt.8) ihour=8

		iapo(1)=aap(is,9)		

C There must be at least 20 indices available
        if((is-1)*8+ihour.lt.20) goto 21      	

C assemble Ap values as needed by MSIS
        j1=ihour+1
        do i=1,ihour
           iap(j1-i)=aap(is,i)
           enddo

        j1=ihour+9
        do i=1,8
            iap(j1-i)=aap(is-1,i)
        	enddo

        j1=ihour+17
        j2=8-(20-ihour-8)+1
        if(j2.lt.1) j2=1
        do i=j2,8
            iap(j1-i)=aap(is-2,i)
            enddo

        if(ihour.lt.4) then
          j1=ihour+25
          j2=8-(20-ihour-16)+1
          do i=j2,8
            iap(j1-i)=aap(is-3,i)
            enddo
          endif

        do 25 i=1,4 
25         iapo(i+1)=iap(i)*1.0

        sum1=0.
        sum2=0.
        do 26 i=1,8
           sum1=sum1+iap(4+i)
26         sum2=sum2+iap(12+i)
c        iapo(6)=int(sum1/8.+.5)
c        iapo(7)=int(sum2/8.+.5)
        iapo(6)=sum1/8.
        iapo(7)=sum2/8.
        goto 20
        
21      if(mess) write(konsol,100)
100     format(1X,'APFMSIS: No Ap dependence because date is not',
     &       ' covered by APF107.DAT indices file')
        IAPO(2)=-5.0
20    	continue

	  RETURN
      END
C
C
        SUBROUTINE APF_ONLY(IYYYY,IMN,ID,F107D,F107PD,F107_81,F107_365,
     *        IAPDA,ISDATE)
c-----------------------------------------------------------------------
c Finds daily F10.7, daily Ap, and 81-day and 365-day F10.7 index: 
c
c    INPUTS:   IYYYY (yyyy)	year 
c              IMN (mm)		month 
c	           ID (dd)		day 
c    OUTPUT:   F107D		F10.7 index for the day (adjusted 
c								to 1AU)
C              F107PD  		F10.7 index for one day prior (used in MSIS)
c              F107_81		F10.7 average over 3 solar rotations
c                               (81 days, centered on the current day) 
c              F107_365 	F10.7 12-month running mean
c              IAPDA		Daily Ap
c              ISDATE		Array-index for the specified date (for
c                                use in APF subroutine.
c 
c Using COMMON/apfa/ for indices
c
c Is used for vdrift and foeedi.
c
c If date is outside the range of indices file than F107D=F107_81=-11.1  
c-----------------------------------------------------------------------

        INTEGER		aap(23000,9),iiap(8),lm(12)
        DIMENSION 	af107(23000,3)
        LOGICAL 	mess

        common 		/iounit/konsol,mess /apfa/aap,af107,nf107

        DATA LM/31,28,31,30,31,30,31,31,30,31,30,31/

        IYBEG=1958
        if(iyyyy.lt.IYBEG) goto 21   ! APF107.DAT starts at Jan 1, 1958

        is=0
        do i=IYBEG,iyyyy-1
            nyd=365
            if(i/4*4.eq.i) nyd=366	! leap year
            IS=IS+nyd
            enddo
        lm(2)=28
        if(iyyyy/4*4.eq.iyyyy) lm(2)=29	  ! leap year
        do i=1,IMN-1
            IS=IS+LM(i)
            ENDDO
        
        IS=IS+ID
        ISDATE = IS
        if(IS.gt.nf107) goto 21

        F107D = AF107(IS,1)
        F107PD = F107D
		if(IS.gt.1) F107PD = AF107(IS-1,1)
        F107_81=AF107(IS,2)
 		if(F107_81.lt.-4.) F107_81=F107D
        F107_365=AF107(IS,3)
  		if(F107_365.lt.-4.) F107_365=F107D
        IAPDA=AAP(is,9)
        goto 20
                 
21      if(mess) write(konsol,100)
100     format(1X,'APF_ONLY: Date is outside range of F10.7D indices',
     &    ' file (F10.7D = F10.7_81 = F10.7RM12).')
        F107D = -11.1
        F107_81 = -11.1
        F107_365 = -11.1
        IAPDA = -11
    
20    RETURN
      END
C      
C
C----------------------STORM MODEL --------------------------------
C
      SUBROUTINE CONVER(rga,rgo,rgma)

C     This subroutine converts a geographic latitude and longitude
C     location to a corrected geomagnetic latitude.
C
C     INPUT: 
C       geographic latitude   -90. to +90.
C       geographic longitude  0. to 360. positive east from Greenwich.
C
C     OUTPUT:
C       corrected geomagnetic latitude	-90. to +90.


      DIMENSION CORMAG(20,91)      
      DATA ((CORMAG(i,j),i=1,20),j=1,31)/
     +163.68,163.68,163.68,163.68,163.68,163.68,
     +163.68,163.68,163.68,163.68,163.68,163.68,163.68,163.68,
     +163.68,163.68,163.68,163.68,163.68,163.68,162.60,163.12,
     +163.64,164.18,164.54,164.90,165.16,165.66,166.00,165.86,
     +165.20,164.38,163.66,162.94,162.42,162.00,161.70,161.70,
     +161.80,162.14,161.20,162.18,163.26,164.44,165.62,166.60,
     +167.42,167.80,167.38,166.82,166.00,164.66,163.26,162.16,
     +161.18,160.40,159.94,159.80,159.98,160.44,159.80,161.14,
     +162.70,164.50,166.26,167.90,169.18,169.72,169.36,168.24,
     +166.70,164.80,162.90,161.18,159.74,158.60,157.94,157.80,
     +157.98,158.72,158.40,160.10,162.02,164.28,166.64,169.00,
     +170.80,171.72,171.06,169.46,167.10,164.64,162.18,160.02,
     +158.20,156.80,156.04,155.80,156.16,157.02,157.00,158.96,
     +161.24,163.86,166.72,169.80,172.42,173.72,172.82,170.34,
     +167.30,164.22,161.34,158.74,156.60,155.00,154.08,153.90,
     +154.36,155.36,155.50,157.72,160.36,163.32,166.60,170.20,
     +173.70,175.64,174.18,170.80,167.10,163.56,160.24,157.36,
     +154.96,153.10,152.08,151.92,152.46,153.76,154.10,156.52,
     +159.36,162.52,166.24,170.30,174.62,177.48,175.04,170.82,
     +166.60,162.70,159.02,155.88,153.22,151.20,150.08,149.92,
     +150.64,152.20,152.80,155.32,158.28,161.70,165.58,170.00,
     +174.84,178.46,175.18,170.38,165.80,161.64,157.80,154.38,
     +151.52,149.30,148.18,148.02,148.92,150.60,151.40,154.08,
     +157.18,160.68,164.78,169.40,174.34,177.44,174.28,169.44,
     +164.70,160.34,156.30,152.78,149.72,147.40,146.18,146.04,
     +147.12,149.04,150.10,152.88,156.00,159.58,163.78,168.50,
     +173.28,175.60,172.86,168.14,163.40,158.98,154.88,151.10,
     +147.98,145.50,144.18,144.14,145.40,147.48,148.80,151.68,
     +154.88,158.48,162.68,167.40,171.76,173.60,171.12,166.68,
     +162.00,157.48,153.28,149.50,146.18,143.50,142.18,142.24,
     +143.68,145.98,147.50,150.54,153.68,157.28,161.42,166.10,
     +170.10,171.48,169.22,164.98,160.40,155.88,151.68,147.80,
     +144.34,141.60,140.18,140.26,141.98,144.62,146.30,149.34,
     +152.48,155.98,160.08,164.60,168.34,169.38,167.20,163.18,
     +158.60,154.18,149.98,146.02,142.54,139.70,138.18,138.46,
     +140.26,143.16,145.10,148.14,151.18,154.60,158.68,163.10,
     +166.48,167.28,165.18,161.32,156.90,152.48,148.28,144.32,
     +140.74,137.80,136.22,136.48,138.64,141.76,143.90,146.98,
     +149.98,153.30,157.24,161.40,164.52,165.16,162.86,159.42,
     +155.00,150.68,146.48,142.52,138.94,135.90,134.22,134.68,
     +137.02,140.40,142.70,145.84,148.76,151.92,155.74,159.70,
     +162.52,162.96,160.98,157.42,153.10,148.84,144.68,140.82,
     +137.20,134.00,132.32,132.80,135.42,139.10,141.60,144.74,
     +147.46,150.52,154.20,158.00,160.46,160.76,158.86,155.36,
     +151.20,146.94,142.88,139.02,135.40,132.10,130.32,131.00,
     +133.80,137.74,140.50,143.58,146.24,149.12,152.60,156.20,
     +158.40,158.66,156.76,153.36,149.30,145.04,141.08,137.30,
     +133.60,130.30,128.42,129.12,132.28,136.44,139.30,142.48,
     +144.94,147.64,150.48,154.30,156.34,156.36,154.56,151.26,
     +147.30,143.14,139.20,135.50,131.90,128.40,126.52,127.32,
     +130.76,135.18,138.20,141.28,143.72,146.24,149.26,152.40,
     +154.24,154.16,152.36,149.16,145.30,141.24,137.30,133.70,
     +130.10,126.60,124.62,125.54,129.16,133.92,137.10,140.18,
     +142.42,144.66,147.62,150.50,152.18,151.96,150.16,147.10,
     +143.30,139.24,135.50,131.90,128.36,124.80,122.72,123.74,
     +127.64,132.62,135.90,139.02,141.12,143.18,145.92,148.60,
     +149.98,149.76,148.04,145.00,141.20,137.30,133.60,130.10,
     +126.60,123.00,120.86,121.96,126.12,131.36,134.80,137.88,
     +139.80,141.68,144.08,146.60,147.88,147.56,145.84,142.90,
     +139.20,135.30,131.70,128.28,124.86,121.30,118.96,120.18,
     +124.70,130.16,133.60,136.72,138.48,140.10,142.38,144.60,
     +145.72,145.34,143.64,140.80,137.10,133.30,129.72,126.48,
     +123.10,119.50,117.16,118.48,123.18,128.86,132.40,135.42,
     +137.08,138.50,140.54,142.60,143.52,143.06,141.44,138.70,
     +135.10,131.30,127.82,124.58,121.40,117.70,115.26,116.70,
     +121.66,127.60,131.20,134.22,135.66,136.82,138.70,140.60,
     +141.36,140.86,139.24,136.50,133.00,129.30,125.92,122.78,
     +119.60,116.00,113.40,114.92,120.16,126.30,130.00,132.92,
     +134.24,135.14,136.80,138.60,139.16,138.64,137.12,134.40,
     +130.90,127.20,123.92,120.96,117.90,114.20,111.56,113.12,
     +118.64,124.90,128.70,131.56,132.74,133.44,134.90,136.50,
     +137.00,136.36,134.82,132.30,128.70,125.16,121.94,119.06,
     +116.10,112.50,109.70,111.42,117.14,123.60,127.30,130.16,
     +131.22,131.66,133.00,134.50,134.80,134.14,132.62,130.14,
     +126.60,123.06,119.94,117.16,114.30,110.70,107.80,109.64,
     +115.62,122.24,125.90,128.76,129.62,129.96,131.06,132.40,
     +132.60,131.86,130.42,128.00,124.50,120.96,117.96,115.26,
     +112.54,108.90,105.94,107.86,114.02,120.84/

      DATA ((CORMAG(i,j),i=1,20),j=32,61)/
     +124.05,126.79,
     +127.55,127.83,128.90,130.21,130.41,129.71,128.33,125.96,
     +122.49,118.96,115.97,113.26,110.52,106.89,104.01,106.00,
     +112.21,119.06,122.19,124.82,125.48,125.69,126.73,128.03,
     +128.22,127.55,126.23,123.92,120.47,116.97,113.97,111.26,
     +108.50,104.89,102.08,104.14,110.41,117.29,120.34,122.85,
     +123.40,123.56,124.57,125.84,126.03,125.40,124.14,121.88,
     +118.46,114.97,111.98,109.26,106.48,102.88,100.15,102.28,
     +108.60,115.51,118.49,120.88,121.33,121.42,122.40,123.65,
     +123.84,123.24,122.04,119.83,116.45,112.97,109.98,107.26,
     +104.46,100.87,098.22,100.42,106.79,113.74,116.63,118.91,
     +119.26,119.29,120.24,121.47,121.65,121.09,119.95,117.79,
     +114.43,110.98,107.99,105.26,102.44,098.87,096.29,098.56,
     +104.98,111.96,114.78,116.94,117.19,117.15,118.07,119.28,
     +119.46,118.93,117.86,115.75,112.42,108.98,106.00,103.26,
     +100.42,096.86,094.36,096.70,103.18,110.19,112.93,114.97,
     +115.12,115.02,115.91,117.09,117.27,116.78,115.76,113.71,
     +110.41,106.98,104.00,101.26,098.40,094.85,092.43,094.84,
     +101.37,108.41,111.07,113.00,113.04,112.88,113.74,114.91,
     +115.08,114.62,113.67,111.67,108.39,104.99,102.01,099.26,
     +096.38,092.85,090.51,092.97,099.56,106.64,109.22,111.03,
     +110.97,110.75,111.58,112.72,112.89,112.47,111.57,109.63,
     +106.38,102.99,100.01,097.26,094.36,090.84,088.58,091.11,
     +097.75,104.86,107.37,109.06,108.90,108.61,109.41,110.53,
     +110.70,110.31,109.48,107.59,104.37,100.99,098.02,095.26,
     +092.34,088.83,086.65,089.25,095.95,103.09,105.51,107.09,
     +106.83,106.48,107.25,108.35,108.51,108.16,107.39,105.55,
     +102.35,099.00,096.03,093.26,090.32,086.83,084.72,087.39,
     +094.14,101.31,103.66,105.12,104.76,104.34,105.08,106.16,
     +106.32,106.00,105.29,103.50,100.34,097.00,094.03,091.26,
     +088.30,084.82,082.79,085.53,092.33,099.54,101.81,103.15,
     +102.68,102.21,102.92,103.97,104.13,103.85,103.20,101.46,
     +098.33,095.00,092.04,089.26,086.28,082.81,080.86,083.67,
     +090.52,097.76,099.95,101.18,100.61,100.07,100.75,101.79,
     +101.94,101.69,101.10,099.42,096.31,093.01,090.04,087.26,
     +084.26,080.81,078.93,081.81,088.72,095.99,098.10,099.21,
     +098.54,097.94,098.59,099.60,099.75,099.54,099.01,097.38,
     +094.30,091.01,088.05,085.26,082.24,078.80,077.00,079.95,
     +086.91,094.21,096.25,097.24,096.47,095.81,096.43,097.41,
     +097.56,097.39,096.92,095.34,092.29,089.01,086.06,083.26,
     +080.22,076.79,075.07,078.09,085.10,092.43,094.39,095.27,
     +094.40,093.67,094.26,095.23,095.37,095.23,094.82,093.30,
     +090.27,087.02,084.06,081.26,078.20,074.79,073.14,076.23,
     +083.30,090.66,092.54,093.30,092.32,091.54,092.10,093.04,
     +093.18,093.08,092.73,091.26,088.26,085.02,082.07,079.26,
     +076.18,072.78,071.21,074.37,081.49,088.88,090.69,091.33,
     +090.25,089.40,089.93,090.85,090.99,090.92,090.63,089.21,
     +086.25,083.02,080.07,077.26,074.16,070.77,069.28,072.51,
     +079.68,087.11,088.83,089.36,088.18,087.27,087.77,088.67,
     +088.80,088.77,088.54,087.17,084.23,081.03,078.08,075.26,
     +072.14,068.77,067.35,070.65,077.87,085.33,086.98,087.39,
     +086.11,085.13,085.60,086.48,086.61,086.61,086.45,085.13,
     +082.22,079.03,076.09,073.26,070.12,066.76,065.42,068.79,
     +076.07,083.56,085.13,085.42,084.04,083.00,083.44,084.29,
     +084.42,084.46,084.35,083.09,080.21,077.03,074.09,071.26,
     +068.10,064.75,063.49,066.93,074.26,081.78,083.27,083.45,
     +081.96,080.86,081.27,082.11,082.23,082.30,082.26,081.05,
     +078.19,075.04,072.10,069.26,066.08,062.75,061.57,065.06,
     +072.45,080.01,081.42,081.48,079.89,078.73,079.11,079.92,
     +080.04,080.15,080.16,079.01,076.18,073.04,070.10,067.26,
     +064.06,060.74,059.64,063.20,070.64,078.23,079.57,079.51,
     +077.82,076.59,076.94,077.73,077.85,077.99,078.07,076.97,
     +074.17,071.04,068.11,065.26,062.04,058.73,057.71,061.34,
     +068.84,076.46,077.71,077.54,075.75,074.46,074.78,075.55,
     +075.66,075.84,075.98,074.93,072.15,069.05,066.12,063.26,
     +060.02,056.73,055.78,059.48,067.03,074.68,075.86,075.57,
     +073.68,072.32,072.61,073.36,073.47,073.68,073.88,072.88,
     +070.14,067.05,064.12,061.26,058.00,054.72,053.85,057.62,
     +065.22,072.91,074.01,073.60,071.60,070.19,070.45,071.17,
     +071.28,071.53,071.79,070.84,068.13,065.05,062.13,059.26,
     +055.98,052.71,051.92,055.76,063.41,071.13,072.15,071.63,
     +069.53,068.05,068.28,068.99,069.09,069.37,069.69,068.80,
     +066.11,063.06,060.13,057.26,053.96,050.71,049.99,053.90,
     +061.61,069.36,070.30,069.66,067.46,065.92,066.12,066.80,
     +066.90,067.22,067.60,066.76,064.10,061.06,058.14,055.26,
     +051.94,048.70,048.06,052.04,059.80,067.58/

      DATA ((CORMAG(i,j),i=1,20),j=62,91)/
     +067.70,067.06,
     +065.08,063.72,063.98,064.60,064.80,065.12,065.60,064.86,
     +062.40,059.26,056.24,053.18,049.84,046.60,046.12,050.12,
     +057.52,064.80,064.90,064.42,062.70,061.62,061.78,062.40,
     +062.60,063.04,063.58,063.00,060.60,057.46,054.42,051.18,
     +047.70,044.60,044.22,048.02,055.06,061.92,062.10,061.72,
     +060.32,059.50,059.68,060.20,060.46,060.94,061.58,061.00,
     +058.70,055.66,052.52,049.18,045.60,042.50,042.22,046.00,
     +052.60,058.98,059.20,059.18,058.12,057.32,057.48,058.00,
     +058.30,058.84,059.48,059.04,056.90,053.86,050.62,047.10,
     +043.50,040.50,040.28,043.98,050.22,056.18,056.40,056.64,
     +055.84,055.20,055.38,055.80,056.16,056.84,057.48,057.04,
     +055.10,052.06,048.70,045.10,041.40,038.40,038.28,041.88,
     +047.94,053.44,053.70,054.14,053.56,053.10,053.24,053.70,
     +054.06,054.74,055.38,055.14,053.20,050.26,046.80,043.10,
     +039.34,036.40,036.38,039.96,045.56,050.84,051.10,051.70,
     +051.36,051.00,051.14,051.50,051.96,052.64,053.38,053.08,
     +051.30,048.36,044.90,041.02,037.24,034.40,034.38,037.86,
     +043.28,048.20,048.50,049.26,049.18,048.90,049.04,049.40,
     +049.86,050.64,051.28,051.08,049.40,046.46,042.98,039.02,
     +035.14,032.40,032.48,035.72,041.00,045.70,046.00,046.96,
     +046.98,046.80,046.94,047.30,047.76,048.54,049.28,049.08,
     +047.40,044.56,041.08,037.02,033.14,030.40,030.58,033.84,
     +038.72,043.20,043.50,044.62,044.80,044.80,044.94,045.20,
     +045.76,046.54,047.18,046.98,045.50,042.66,039.08,035.02,
     +031.14,028.40,028.58,031.82,036.52,040.80,041.20,042.32,
     +042.54,042.70,042.84,043.20,043.66,044.44,045.08,044.98,
     +043.50,040.76,037.08,033.04,029.04,026.40,026.68,029.82,
     +034.34,038.40,038.80,040.12,040.60,040.70,040.84,041.10,
     +041.62,042.34,042.98,042.88,041.50,038.76,035.18,031.04,
     +027.14,024.50,024.78,027.70,032.14,036.06,036.50,037.88,
     +038.50,038.68,038.84,039.10,039.56,040.34,040.88,040.82,
     +039.40,036.76,033.18,029.12,025.14,022.50,022.88,025.90,
     +029.96,033.86,034.30,035.68,036.42,036.68,036.84,037.10,
     +037.56,038.24,038.88,038.72,037.40,034.76,031.18,027.12,
     +023.14,020.60,020.98,023.90,027.88,031.66,032.10,033.58,
     +034.32,034.68,034.84,035.10,035.56,036.24,036.78,036.62,
     +035.30,032.72,029.18,025.14,021.24,018.70,019.08,021.90,
     +025.88,029.42,029.90,031.48,032.32,032.68,032.84,033.10,
     +033.56,034.22,034.68,034.42,033.20,030.72,027.28,023.22,
     +019.34,016.80,017.24,020.00,023.78,027.32,027.70,029.38,
     +030.24,030.68,030.94,031.20,031.66,032.22,032.58,032.32,
     +031.10,028.62,025.28,021.32,017.48,015.00,015.38,018.18,
     +021.80,025.22,025.70,027.28,028.24,028.78,029.04,029.30,
     +029.66,030.22,030.50,030.22,029.00,026.62,023.30,019.42,
     +015.64,013.10,013.54,016.28,019.80,023.12,023.60,025.24,
     +026.24,026.78,027.14,027.40,027.76,028.22,028.40,028.12,
     +026.80,024.52,021.30,017.52,013.78,011.30,011.74,014.48,
     +017.90,021.12,021.60,023.24,024.34,024.88,025.24,025.50,
     +025.86,026.22,026.40,025.98,024.70,022.48,019.40,015.72,
     +012.04,009.50,009.94,012.58,016.02,019.12,019.60,021.24,
     +022.34,022.98,023.34,023.70,024.00,024.30,024.40,023.88,
     +022.60,020.48,017.52,014.00,010.34,007.80,008.18,010.88,
     +014.22,017.18,017.60,019.34,020.44,021.16,021.54,021.90,
     +022.16,022.40,022.32,021.78,020.60,018.48,015.62,012.20,
     +008.68,006.00,006.44,009.18,012.42,015.28,015.80,017.44,
     +018.54,019.26,019.74,020.10,020.30,020.50,020.32,019.72,
     +018.50,016.54,013.84,010.68,007.14,004.40,004.74,007.58,
     +010.74,013.48,014.00,015.54,016.74,017.46,017.94,018.30,
     +018.50,018.58,018.32,017.72,016.50,014.64,012.24,009.18,
     +005.84,002.90,003.30,006.16,009.14,011.84,012.30,013.78,
     +014.94,015.66,016.24,016.50,016.70,016.70,016.42,005.78,
     +014.60,012.90,010.66,007.86,004.88,001.60,001.72,004.96,
     +007.84,010.24,010.70,012.14,013.24,013.96,014.44,014.80,
     +014.90,014.88,014.52,013.92,012.80,011.30,009.28,006.94,
     +004.32,001.80,001.94,004.34,006.78,008.94,009.40,010.58,
     +011.64,012.36,012.74,013.10,013.20,013.08,012.72,012.12,
     +011.10,009.86,008.30,006.50,004.60,003.10,003.16,004.50,
     +006.20,007.90,008.40,009.42,010.14,010.76,011.14,011.40,
     +011.40,011.38,011.02,010.46,009.70,008.72,007.64,006.46,
     +005.42,004.60,004.70,005.34,006.24,007.36,007.90,008.46,
     +008.92,009.28,009.54,009.70,009.70,009.68,009.42,009.06,
     +008.60,008.08,007.56,007.02,006.56,006.30,006.30,006.52,
     +006.96,007.38,008.15,008.15,008.15,008.15,008.15,008.15,
     +008.15,008.15,008.15,008.15,008.15,008.15,008.15,008.15,
     +008.15,008.15,008.15,008.15,008.15,008.15/

C     Data Input      
      rlan = rga
      rlo = rgo      
      
C     From "normal" geographic latitude 
C     to angle from South Pole.       
      rla = rlan + 90

      IF (rlo .EQ. 360) THEN
      	rlo = 0
        END IF

C     PROXIMITY

C     coefficients of the latitudinal points		
      LA1 = (INT(rla/2)+1)
      LA2 = LA1 + 1
      if(la2.gt.91) la2=91

C     coefficients of the longitudinal points		
      LO1 = (INT(rlo/18)+1)
corr      LO2 = LO1 + 1
      LO2 = MOD(LO1,20) + 1 

C     Four points of Geomagnetic Coordinates
      gm1 = CORMAG(LO1,LA1)
      gm2 = CORMAG(LO1,LA2) 
      gm3 = CORMAG(LO2,LA1)
      gm4 = CORMAG(LO2,LA2)

C     latitudinal points		
C      X1 = ABS(rla - (INT(rla)))                        
C      X2 = 2. - X1
	  x = (rla/2.0 - (INT(rla/2.0)))

C     longitudinal points		
C      Y1 = ABS(rlo - (INT(rlo)))                        
C      Y2 = 18. - Y1
      y =(rlo/18.0 - (INT(rlo/18.0))) 

C     X AND Y VALUES
C      x = X1 / (X1 + X2)
C      y = Y1 / (Y1 + Y2)

C     INTERPOLATION
      gmla = gm1 * (1 - x) * (1 - y) + gm2 * (1 - y) * (x) + gm3 * (y)
     1 * (1 - x) + gm4 * (x) * (y)

C     OUTPUT OF THE PROGRAM
C     From corrected geomagnetic latitude from North Pole
C     to "normal"  geomagnetic latitude.       
      rgma = 90. - gmla

      END
c
c
      SUBROUTINE STORM(ap,rga,rgo,coor,rgma,ut,doy,cf)
C----------------------------------------------------------------------
C      Fortran code to obtain the foF2 storm-time correction factor at 
C      a given location and time, using the current and the 12 previous
C      ap values as input.
C
C      ap ---> (13 elements integer array). Array with the preceeding
C              13 value of the 3-hourly ap index. The 13th value
C              in the array will contain the ap at the UT of interest,
C              the 12th value will contain the 1st three hourly interval
C              preceeding the time of interest, and so on to the 1st
C              ap value at the earliest time.
C     coor --> (integer). If coor = 2, rga should contain the 
C                         geomagnetic latitude.
C                         If coor = 1, rga should contain the 
C                         geographic latitude.
C     rga ---> (real, -90 to 90) geographic or geomagnetic latitude.
C     rgo ---> (real, 0 to 360, positive east from Greenwich.)
C                           geographic longitude, only used if coor=1.
C     ut  ---> (integer, hours 00 to 23) Universal Time of interest.
C     doy ---> (integer, 1 to 366)Day of the year.
C     cf  ---> (real) The output; the storm-time correction factor used
C              to scale foF2, foF2 * cf.
C     rgma --> corrected magnetic latitude calculated from rga and rgo
C
C     This model and computer code was developed by E. Araujo-Pradere,
C     T. Fuller-Rowell and M. Condrescu, SEC, NOAA, Boulder, USA
C     Ref: 
C     T. Fuller-Rowell, E. Araujo-Pradere, and M. Condrescu, An 
C       Empirical Ionospheric Storm-Time Ionospheric Correction Model,
C       Adv. Space Res. 8, 8, 15-24, 2000.
C----------------------------------------------------------------------
C     DIMENSIONS AND COEFFICIENTS VALUES

      DIMENSION c4(20)
      DATA c4/0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,
     +0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,
     +0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00/

      DIMENSION c3(20)
      DATA c3/0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,-9.44E-12,
     +0.00E+00,3.04E-12,0.00E+00,9.32E-12,-1.07E-11,0.00E+00,0.00E+00,
     +0.00E+00,1.09E-11,0.00E+00,0.00E+00,0.00E+00,0.00E+00,-1.01E-11/

      DIMENSION c2(20)
      DATA c2/1.16E-08,0.00E+00,0.00E+00,-1.46E-08,0.00E+00,9.86E-08,
     +2.25E-08,-1.67E-08,-1.62E-08,-9.42E-08,1.17E-07,4.32E-08,3.97E-08,
     +3.13E-08,-8.04E-08,3.91E-08,2.58E-08,3.45E-08,4.76E-08,1.13E-07/

      DIMENSION c1(20)
      DATA c1/-9.17E-05,-1.37E-05,0.00E+00,7.14E-05,0.00E+00,-3.21E-04,
     +-1.66E-04,-4.10E-05,1.36E-04,2.29E-04,-3.89E-04,-3.08E-04,
     +-2.81E-04,-1.90E-04,4.76E-05,-2.80E-04,-2.07E-04,-2.91E-04,
     +-3.30E-04,-4.04E-04/

      DIMENSION c0(20)
      DATA c0/1.0136E+00,1.0478E+00,1.00E+00,1.0258E+00,1.00E+00,
     +1.077E+00,1.0543E+00,1.0103E+00,9.9927E-01,9.6876E-01,1.0971E+00,
     +1.0971E+00,1.0777E+00,1.1134E+00,1.0237E+00,1.0703E+00,1.0248E+00,
     +1.0945E+00,1.1622E+00,1.1393E+00/

      DIMENSION fap(36)
      DATA fap/0.,0.,0.037037037,0.074074074,0.111111111,0.148148148,
     10.185185185,0.222222222,0.259259259,0.296296296,0.333333333,
     20.37037037,0.407407407,0.444444444,0.481481481,0.518518519,
     30.555555556,0.592592593,0.62962963,0.666666667,0.703703704,
     40.740740741,0.777777778,0.814814815,0.851851852,0.888888889,
     50.925925926,0.962962963,1.,0.66666667,0.33333334,0.,0.333333,
     60.666666,1.,0.7/

      integer code(8,6)
      data code/3,4,5,4,3,2,1,2,3,2,1,2,3,4,5,4,8,7,6,7,8,9,10,9,
     *13,12,11,12,13,14,15,14,18,17,16,17,18,19,20,19,18,17,16,17,
     *18,19,20,19/

      INTEGER ape(39)
      INTEGER ap(13)
      INTEGER ut,doy,dayno,coor,s1,s2,l1,l2
      REAL rgma, rap, rga, rgo, rs, rl

C      CALLING THE PROGRAM TO CONVERT TO GEOMAGNETIC COORDINATES

       IF (coor .EQ. 1) THEN

           CALL CONVER (rga,rgo,rgma)

       ELSE IF (coor .EQ. 2) THEN
                rgma = rga

       ELSE

          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)'   Wrong Coordinates Selection -------- >>', coor
          WRITE (6,*)' '
          GOTO 100
       ENDIF

C FROM 3-HOURLY TO HOURLY ap (New, interpolates between the three hourly ap values)

       ape(1)=ap(1)
       ape(2)=ap(1)
       ape(38)=ap(13)
       ape(39)=ap(13)

       DO k = 1,13
          i = (k * 3) - 1
          ape(i) = ap(k)
          END DO

       DO k = 1,12
          i = k * 3
          ape(i) = (ap(k)*2 + ap(k+1))/3.0
          END DO

       DO k = 2,13
          i = (k * 3) - 2
          ape(i) = (ap(k-1) + ap(k)*2)/3.0
          END DO

C     FROM 3-HOURLY TO HOURLY ap (old version without interpolation)
c      i = 1
c      DO 10 k = 1,13
c         DO j = 1,3
c            ape(i) = ap(k)
c            i = i + 1
c            END DO
c10    CONTINUE

C     TO OBTAIN THE INTEGRAL OF ap.
C     INTEGRAL OF ap

      if(ut.eq.24) ut=0
      IF (ut .EQ. 0 .OR. ut .EQ. 3 .OR. ut .EQ. 6 .OR. ut .EQ. 9 .OR.
     1ut .EQ. 12 .OR. ut .EQ. 15 .OR. ut .EQ. 18 .OR. ut .EQ. 21) THEN
          k = 1
      ELSE IF (ut .EQ. 1 .OR. ut .EQ. 4 .OR. ut .EQ. 7 .OR. ut .EQ. 10
     1.OR.ut .EQ. 13 .OR. ut .EQ. 16 .OR. ut .EQ. 19 .OR. ut .EQ. 22)
     2THEN
          k = 2
      ELSE IF (ut .EQ. 2 .OR. ut .EQ. 5 .OR. ut .EQ. 8 .OR. ut .EQ. 11
     1.OR. ut .EQ. 14 .OR. ut .EQ. 17 .OR. ut .EQ. 20 .OR. ut .EQ. 23)
     2THEN
          k = 3

      ELSE

          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)'  Wrong Universal Time value -------- >>', ut
          WRITE (6,*)' '
          GOTO 100

      END IF

      rap = 0

      DO j = 1,36
      rap = rap + fap(j) * ape(k+j)
      END DO

      if(rap.le.200.)then
      cf=1.0
      goto 100
      end if

      if(doy.gt.366.or.doy.lt.1)then
          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)'      Wrong Day of Year value --- >>', doy
          WRITE (6,*)' '
          GOTO 100
      end if

      if(rgma.gt.90.0.or.rgma.lt.-90.0)then
          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)'   Wrong GEOMAGNETIC LATITUDE value --- >>', rgma
          WRITE (6,*)' '
          GOTO 100
      end if

c      write(6,*)rgma

      dayno=doy
      if(rgma.lt.0.0)then
      dayno=doy+172
      if(dayno.gt.365)dayno=dayno-365
      end if

      if (dayno.ge.82) rs=(dayno-82.)/45.6+1.
      if (dayno.lt.82) rs=(dayno+283.)/45.6+1.
      s1=rs
      facs=rs-s1
      s2=s1+1
      if(s2.eq.9) s2=1
c      write(6,*)s1,s2,rs

      rgma = abs(rgma)

      rl=(rgma+10.)/20.+1
      if(rl.eq.6.0)rl=5.9
      l1=rl
      facl=rl-l1
      l2=l1+1
c      write(6,*)l1,l2,rl

C     FACTORS CALCULATIONS

      if(rap.lt.300.)then
      rapf=300.
      n1=code(s1,l1)
      cf1=c4(n1)*(rapf**4)+c3(n1) * (rapf**3) + c2(n1) * (rapf**2) +
     1c1(n1) * rapf + c0(n1)
      n2=code(s1,l2)
      cf2=c4(n2)*(rapf**4)+c3(n2) * (rapf**3) + c2(n2) * (rapf**2) +
     1c1(n2) * rapf + c0(n2)
      n3=code(s2,l1)
      cf3=c4(n3)*(rapf**4)+c3(n3) * (rapf**3) + c2(n3) * (rapf**2) +
     1c1(n3) * rapf + c0(n3)
      n4=code(s2,l2)
      cf4=c4(n4)*(rapf**4)+c3(n4) * (rapf**3) + c2(n4) * (rapf**2) +
     1c1(n4) * rapf + c0(n4)

C     INTERPOLATION

      cf300=cf1*(1 - facs) * (1 - facl) + cf2 * (1 - facs) * (facl) +
     *cf3 * (facs) * (1 - facl) + cf4 * (facs) * (facl)

      cf = (cf300-1.0)*rap/100.-2.*cf300+3.
      goto 100
      end if

      n1=code(s1,l1)
c      write(6,*)n1
      cf1 = c4(n1) * (rap**4) + c3(n1) * (rap**3) + c2(n1) * (rap**2) +
     1c1(n1) * rap + c0(n1)
      n2=code(s1,l2)
      cf2 = c4(n2) * (rap**4) + c3(n2) * (rap**3) + c2(n2) * (rap**2) +
     1c1(n2) * rap + c0(n2)
      n3=code(s2,l1)
      cf3 = c4(n3) * (rap**4) + c3(n3) * (rap**3) + c2(n3) * (rap**2) +
     1c1(n3) * rap + c0(n3)
      n4=code(s2,l2)
      cf4 = c4(n4) * (rap**4) + c3(n4) * (rap**3) + c2(n4) * (rap**2) +
     1c1(n4) * rap + c0(n4)

C     INTERPOLATION

      cf = cf1 * (1 - facs) * (1 - facl) + cf2 * (1 - facs) * (facl) +
     *cf3 * (facs) * (1 - facl) + cf4 * (facs) * (facl)

100   CONTINUE

      RETURN

      END
C
C
      FUNCTION STORME_AP(JDOY,XMLAT,AP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C EMPIRICAL STORM-E MODEL: COMPUTES A STORM-TO-QUIET RATIO (SQR) FACTOR
C TO ADJUST THE QUIESCENT E-REGION PEAK ELECTRON DENSITY TO ACCOUNT FOR
C ENHANCEMENTS DUE TO GEOMAGNETIC ACTIVITY. THE SQR FACTORS WERE
C COMPUTED FROM NO+ 4.3 UM VOLUME EMISSION RATES DERIVED FROM 
C TIMED/SABER LIMB RADIANCE MEASUREMENTS. THE SABER-DERIVED SQR FACTORS
C WERE FIT TO POWE-LAW IN THE ap INDEX.   
C
C INPUT PARAMETERS:
C
C  JDOY      --- DAY OF YEAR (1-365) 
C  XMLAT     --- MAGNETIC LATITUDE (DEGREES)
C  AP        --- ap INDEX
C
C OUTPUT PARAMETER
C 
C  STORME_AP --- STORM-TO-QUIET RATIO (SQR) TO ADJUST QUIESCENT E-REGION
C                PEAK ELECTRON DENSITY TO ACCOUNT FOR GEOMAGNETIC
C                ENHANCEMENTS. SQR COMPUTED FROM A POWER-LAW FIT
C                IN AP-INDEX: SQR=C1*AP**C2+C3
C
C REFERENCES:
C 
C  (1) Mertens et al. [submitted to JASR, 2011]
C  (2) Fernandez et al. [JASR, Vol. 46, 2010]
C  (3) Mertens et al. [Proc. of SPIE, Vol. 7475, 2009]
C  (4) Mertens et al. [Proc. of SPIE, Vol. 6745, 2007]
C  (5) Mertens et al. [JASR, Vol. 39, 2007] 
C 
C SOFTWARE WRITTEN BY Christopher J. Mertens
C                     NASA Langley Research Center
C                     Atmospheric Sciences Competency
C                     21 Langley Blvd., Mail Stop 401B
C                     Hampton, VA 23681-2199
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER(NMLG=37,NDBD=5)
      DIMENSION C3(NMLG,NDBD)
      DIMENSION XMLG(NMLG),IDBD(NDBD),C1(NMLG,NDBD),C2(NMLG,NDBD)
      LOGICAL mess

      COMMON /iounit/konsol,mess

      DATA XMLG/-90.0,-85.0,-80.0,-75.0,-70.0,-65.0,-60.0,-55.0,-50.0,
     &          -45.0,-40.0,-35.0,-30.0,-25.0,-20.0,-15.0,-10.0,-5.0,
     &          0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,
     &          55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0/
      DATA IDBD/79,171,264,354,366/
C
C JDOY        | 0-79 ||80-171||172-264||265-354||355-365|
C
      DATA ((C1(IML,IDB),IDB=1,NDBD),IML=1,NMLG)
     &     / 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-90.0 
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-85.0
     &       0.00000, 0.00508, 0.17360, 0.00000, 0.00000, !-80.0
     &       0.00000, 0.31576, 0.31498, 0.00000, 0.00000, !-75.0
     &       0.00000, 0.39217, 0.40121, 0.00000, 0.00000, !-70.0
     &       0.00000, 0.32634, 0.30179, 0.00000, 0.00000, !-65.0
     &       0.11573, 0.06211, 0.16230, 0.20233, 0.11573, !-60.0
     &       0.00526, 0.00013, 0.00204, 0.04965, 0.00526, !-55.0
     &       0.00011, 0.00013, 0.00018, 0.00040, 0.00011, !-50.0
     &       0.00001, 0.00002, 0.00040, 0.00001, 0.00001, !-45.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-40.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-35.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-30.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-25.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-20.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-15.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-10.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! -5.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !  0.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !  5.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 10.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 15.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 20.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 25.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 30.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 35.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 40.0
     &       0.00000,-0.02738, 0.00004, 0.00001, 0.00000, ! 45.0
     &       0.00001,-0.00022, 0.00001, 0.00004, 0.00001, ! 50.0
     &       0.00126, 0.00062, 0.00011, 0.00345, 0.00126, ! 55.0
     &       0.14923, 0.05483, 0.07113, 0.18282, 0.14923, ! 60.0
     &       0.37361, 0.00000, 0.00000, 0.44592, 0.37361, ! 65.0
     &       0.27792, 0.00000, 0.00000, 0.00804, 0.27792, ! 70.0
     &       0.06445, 0.00000, 0.00000, 0.10315, 0.06445, ! 75.0
     &       0.00149, 0.00000, 0.00000, 0.00073, 0.00149, ! 80.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 85.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000/ ! 90.0  
C
C JDOY        | 0-79 ||80-171||172-264||265-354||355-365|
C
      DATA ((C2(IML,IDB),IDB=1,NDBD),IML=1,NMLG)
     &     / 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-90.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-85.0
     &       0.00000, 1.00000, 0.32415, 0.00000, 0.00000, !-80.0
     &       0.00000, 0.36538, 0.35455, 0.00000, 0.00000, !-75.0
     &       0.00000, 0.41287, 0.38062, 0.00000, 0.00000, !-70.0
     &       0.00000, 0.52224, 0.52810, 0.00000, 0.00000, !-65.0
     &       0.73025, 0.90723, 0.68107, 0.64815, 0.73025, !-60.0
     &       1.29410, 2.06038, 1.47332, 0.84843, 1.29410, !-55.0
     &       1.79442, 1.77511, 1.59906, 1.59141, 1.79442, !-50.0
     &       1.84434, 1.70607, 1.03056, 1.92168, 1.84434, !-45.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-40.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-35.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-30.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-25.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-20.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-15.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !-10.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! -5.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !  0.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, !  5.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 10.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 15.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 20.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 25.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 30.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 35.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 40.0
     &       2.04797, 0.27034, 1.00000, 1.90792, 2.04797, ! 45.0
     &       2.24520, 1.00000, 1.88721, 2.01535, 2.24520, ! 50.0
     &       1.56493, 1.52468, 1.93389, 1.38532, 1.56493, ! 55.0
     &       0.71442, 0.87492, 0.78890, 0.66828, 0.71442, ! 60.0
     &       0.53546, 0.00000, 0.00000, 0.42597, 0.53546, ! 65.0
     &       0.48647, 0.00000, 0.00000, 1.00000, 0.48647, ! 70.0
     &       0.67340, 0.00000, 0.00000, 0.36809, 0.67340, ! 75.0
     &       1.44025, 0.00000, 0.00000, 1.13529, 1.44025, ! 80.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, ! 85.0
     &       0.00000, 0.00000, 0.00000, 0.00000, 0.00000/ ! 90.0 
C
C JDOY        | 0-79 ||80-171||172-264||265-354||355-365|
C
      DATA ((C3(IML,IDB),IDB=1,NDBD),IML=1,NMLG)
     &     / 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-90.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-85.0
     &       1.00000, 1.03132, 0.76703, 1.00000, 1.00000, !-80.0
     &       1.00000, 0.57588, 0.56324, 1.00000, 1.00000, !-75.0
     &       1.00000, 0.41370, 0.38549, 1.00000, 1.00000, !-70.0
     &       1.00000, 0.51704, 0.50217, 1.00000, 1.00000, !-65.0
     &       0.55236, 0.80162, 0.60824, 0.46999, 0.55236, !-60.0
     &       0.90923, 0.99688, 0.96752, 0.67312, 0.90923, !-55.0
     &       0.99338, 0.98486, 0.99503, 0.87473, 0.99338, !-50.0
     &       1.00031, 1.00369, 1.00225, 0.91242, 1.00031, !-45.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-40.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-35.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-30.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-25.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-20.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-15.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !-10.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! -5.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !  0.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, !  5.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! 10.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! 15.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! 20.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! 25.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! 30.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! 35.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! 40.0
     &       1.04797, 1.02319, 0.98581, 1.01457, 1.04797, ! 45.0
     &       1.03332, 0.97016, 0.97807, 0.99044, 1.03332, ! 50.0
     &       1.00633, 0.94822, 0.96340, 0.95363, 1.00633, ! 55.0
     &       0.67902, 0.71540, 0.70230, 0.60821, 0.67902, ! 60.0
     &       0.35017, 1.00000, 1.00000, 0.51033, 0.35017, ! 65.0
     &       0.63358, 1.00000, 1.00000, 1.37782, 0.63358, ! 70.0
     &       0.85724, 1.00000, 1.00000, 0.91942, 0.85724, ! 75.0
     &       0.92703, 1.00000, 1.00000, 1.00502, 0.92703, ! 80.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000, ! 85.0
     &       1.00000, 1.00000, 1.00000, 1.00000, 1.00000/ ! 90.0
C
C ... Find Season-Averaged Coefficient Index 
C
      IDXS=0
      IF(JDOY.LE.IDBD(1)) IDXS=1
      DO IS=2,NDBD
        IF((JDOY.GT.IDBD(IS-1)).AND.(JDOY.LE.IDBD(IS))) IDXS=IS
      ENDDO
      IF(IDXS.EQ.0) THEN
        if(mess) WRITE(konsol,*) 'ERROR IN STORME_AP: ',
     &   'PROBLEM FINDING SEASON-AVERAGED COEFFICIENT',
     &   'DAY OF YEAR = ',JDOY
        STORME_AP=-5.0
        GOTO 222 
      ENDIF
C
C ... Find Magnetic Latitude Coefficient Index
C 
      IDXL=0
      DELG=ABS(XMLG(1)-XMLG(2))
      DELD=DELG/2.0
      YMP=XMLG(1)+DELD
      YMM=XMLG(NMLG)-DELD
      IF((XMLAT.GE.XMLG(1)).AND.(XMLAT.LE.YMP)) IDXL=1
      IF((XMLAT.GT.YMM).AND.(XMLAT.LE.XMLG(NMLG))) IDXL=NMLG
      DO IL=2,NMLG-1
        YMP=XMLG(IL)+DELD
        YMM=XMLG(IL)-DELD
        IF((XMLAT.GT.YMM).AND.(XMLAT.LE.YMP)) IDXL=IL
      ENDDO
      IF(IDXL.EQ.0) THEN
        if(mess) WRITE(konsol,*) 'ERROR IN STORME_AP: ', 
     &     'PROBLEM FINDING MAGNETIC LATITUDE COEFFICIENT',
     &     'MAGNETIC LATITUDE(DEGREES) = ',XMLAT
        STORME_AP=-5.0
        GOTO 222
      ENDIF
C
C ... COMPUTE E-REGION ELECTRON DENSITY GEOMAGNETIC STORM ENHANCEMET
C ... FACTOR (i.e., THE STORM-TO-QUIET RATIO (SQR)) 
C
      SQR=C1(IDXL,IDXS)*AP**(C2(IDXL,IDXS))+C3(IDXL,IDXS)
      IF(SQR.LT.1.0) SQR=1.0
      STORME_AP=SQR
   
222   RETURN
      END    
C
C
C****************************************************************************
C
        subroutine vdrift(xt,xl,param,y)
C-------------------------------------------------------------------
C       SUBROUTINE CALCULATES EQUATORIAL VERTICAL DRIFT AS DESCRIBED 
C       IN SCHERLIESS AND FEJER, JGR, 104, 6829-6842, 1999
C
C       INPUT:   XT: SOLAR LOCAL TIME  [h]
C                XL: GEOGRAPHIC LONGITUDE (+ EAST) [degrees]
C               
C           PARAM: 2-DIM ARRAY (DOY,F10.7CM)
C                  DOY     :Day of Year has to run from 1 to 365(366)
C                  F10.7cm : F10.7cm solar flux (daily value)
C             
C       OUTPUT:   Y: EQUATORIAL VERTICAL DRIFT [m/s]
C
C-------------------------------------------------------------------
c        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT REAL (A-H,O-Z)
       
c        real*8 param(2),coeff(624),coeff1(594),coeff2(30),funct(6)
c        real*8 xt,xl,y
c        real*8 bspl4,bspl4_time,bspl4_long

        real param(2),coeff(624),coeff1(594),coeff2(30),funct(6)
        real xt,xl,y
        real bspl4,bspl4_time,bspl4_long

        integer i,j,ind,il,kk
        integer index_t,dim_t,index_l,dim_l,index,dim,nfunc
        
        data index_t/13/,dim_t/78/,index_l/8/,dim_l/48/,index/104/,
     *   dim/624/,nfunc/6/

        data coeff1/
     *  -10.80592, -9.63722,-11.52666, -0.05716, -0.06288,  0.03564,
     *   -5.80962, -7.86988, -8.50888, -0.05194, -0.05798, -0.00138,
     *    2.09876,-19.99896, -5.11393, -0.05370, -0.06585,  0.03171,
     *  -10.22653, -3.62499,-14.85924, -0.04023, -0.01190, -0.09656,
     *   -4.85180,-26.26264, -6.20501, -0.05342, -0.05174,  0.02419,
     *  -13.98936,-18.10416, -9.30503, -0.01969, -0.03132, -0.01984,
     *  -18.36633,-24.44898,-16.69001,  0.02033, -0.03414, -0.02062,
     *  -20.27621,-16.95623,-36.58234,  0.01445, -0.02044, -0.08297,
     *    1.44450,  5.53004,  4.55166, -0.02356, -0.04267,  0.05023,
     *    5.50589,  7.05381,  1.94387, -0.03147, -0.03548,  0.01166,
     *    3.24165, 10.05002,  4.26218, -0.03419, -0.02651,  0.07456,
     *    7.02218,  0.06708,-11.31012, -0.03252, -0.01021, -0.09008,
     *   -3.47588, -2.82534, -4.17668, -0.03719, -0.01519,  0.06507,
     *   -4.02607,-11.19563,-10.52923, -0.00592, -0.01286, -0.00477,
     *  -11.47478, -9.57758,-10.36887,  0.04555, -0.02249,  0.00528,
     *  -14.19283,  7.86422, -8.76821,  0.05758, -0.02398, -0.04075,
     *   14.58890, 36.63322, 27.57497,  0.01358, -0.02316,  0.04723,
     *   12.53122, 29.38367, 21.40356, -0.00071, -0.00553,  0.01484,
     *   18.64421, 26.27327, 18.32704,  0.00578,  0.03349,  0.11249,
     *    4.53014,  6.15099,  7.41935, -0.02860, -0.00395, -0.08394,
     *   14.29422,  9.77569,  2.85689, -0.00107,  0.04263,  0.10739,
     *    7.17246,  4.40242, -1.00794,  0.00089,  0.01436,  0.00626,
     *    7.75487,  5.01928,  4.36908,  0.03952, -0.00614,  0.03039,
     *   10.25556,  8.82631, 24.21745,  0.05492, -0.02968,  0.00177,
     *   21.86648, 24.03218, 39.82008,  0.00490, -0.01281, -0.01715,
     *   19.18547, 23.97403, 34.44242,  0.01978,  0.01564, -0.02434,
     *   26.30614, 14.22662, 31.16844,  0.06495,  0.19590,  0.05631,
     *   21.09354, 25.56253, 29.91629, -0.04397, -0.08079, -0.07903,
     *   28.30202, 16.80567, 38.63945,  0.05864,  0.16407,  0.07622,
     *   22.68528, 25.91119, 40.45979, -0.03185, -0.01039, -0.01206,
     *   31.98703, 24.46271, 38.13028, -0.08738, -0.00280,  0.01322,
     *   46.67387, 16.80171, 22.77190, -0.13643, -0.05277, -0.01982,
     *   13.87476, 20.52521,  5.22899,  0.00485, -0.04357,  0.09970,
     *   21.46928, 13.55871, 10.23772, -0.04457,  0.01307,  0.06589,
     *   16.18181, 16.02960,  9.28661, -0.01225,  0.14623, -0.01570,
     *   18.16289, -1.58230, 14.54986, -0.00375, -0.00087,  0.04991,
     *   10.00292, 11.82653,  0.44417, -0.00768,  0.15940, -0.01775,
     *   12.15362,  5.65843, -1.94855, -0.00689,  0.03851,  0.04851,
     *   -1.25167,  9.05439,  0.74164,  0.01065,  0.03153,  0.02433,
     *  -15.46799, 18.23132, 27.45320,  0.00899, -0.00017,  0.03385,
     *    2.70396, -0.87077,  6.11476, -0.00081,  0.05167, -0.08932,
     *    3.21321, -1.06622,  5.43623,  0.01942,  0.05449, -0.03084,
     *   17.79267, -3.44694,  7.10702,  0.04734, -0.00945,  0.11516,
     *    0.46435,  6.78467,  4.27231, -0.02122,  0.10922, -0.03331,
     *   15.31708,  1.70927,  7.99584,  0.07462,  0.07515,  0.08934,
     *    4.19893,  6.01231,  8.04861,  0.04023,  0.14767, -0.04308,
     *    9.97541,  5.99412,  5.93588,  0.06611,  0.12144, -0.02124,
     *   13.02837, 10.29950, -4.86200,  0.04521,  0.10715, -0.05465,
     *    5.26779,  7.09019,  1.76617,  0.09339,  0.22256,  0.09222,
     *    9.17810,  5.27558,  5.45022,  0.14749,  0.11616,  0.10418,
     *    9.26391,  4.19982, 12.66250,  0.11334,  0.02532,  0.18919,
     *   13.18695,  6.06564, 11.87835,  0.26347,  0.02858,  0.14801,
     *   10.08476,  6.14899, 17.62618,  0.09331,  0.08832,  0.28208,
     *   10.75302,  7.09244, 13.90643,  0.09556,  0.16652,  0.22751,
     *    6.70338, 11.97698, 18.51413,  0.15873,  0.18936,  0.15705,
     *    5.68102, 23.81606, 20.65174,  0.19930,  0.15645,  0.08151,
     *   29.61644,  5.49433, 48.90934,  0.70710,  0.40791,  0.26325,
     *   17.11994, 19.65380, 44.88810,  0.45510,  0.41689,  0.22398,
     *    8.45700, 34.54442, 27.25364,  0.40867,  0.37223,  0.22374,
     *   -2.30305, 32.00660, 47.75799,  0.02178,  0.43626,  0.30187,
     *    8.98134, 33.01820, 33.09674,  0.33703,  0.33242,  0.41156,
     *   14.27619, 20.70858, 50.10005,  0.30115,  0.32570,  0.45061,
     *   14.44685, 16.14272, 45.40065,  0.37552,  0.31419,  0.30129,
     *    6.19718, 18.89559, 28.24927,  0.08864,  0.41627,  0.19993,
     *    7.70847, -2.36281,-21.41381,  0.13766,  0.05113, -0.11631,
     *   -9.07236,  3.76797,-20.49962,  0.03343,  0.08630,  0.00188,
     *   -8.58113,  5.06009, -6.23262,  0.04967,  0.03334,  0.24214,
     *  -27.85742,  8.34615,-27.72532, -0.08935,  0.15905, -0.03655,
     *    2.77234,  0.14626, -4.01786,  0.22338, -0.04478,  0.18650,
     *    5.61364, -3.82235,-16.72282,  0.26456, -0.03119, -0.08376,
     *   13.35847, -6.11518,-16.50327,  0.28957, -0.01345, -0.19223,
     *   -5.37290, -0.09562,-27.27889,  0.00266,  0.22823, -0.35585,
     *  -15.29676,-18.36622,-24.62948, -0.31299, -0.23832, -0.08463,
     *  -23.37099,-13.69954,-26.71177, -0.19654, -0.18522, -0.20679,
     *  -26.33762,-15.96657,-42.51953, -0.13575, -0.00329, -0.28355,
     *  -25.42140,-14.14291,-21.91748, -0.20960, -0.19176, -0.32593,
     *  -23.36042,-23.89895,-46.05270, -0.10336,  0.03030, -0.21839,
     *  -19.46259,-21.27918,-32.38143, -0.17673, -0.15484, -0.11226,
     *  -19.06169,-21.13240,-34.01677, -0.25497, -0.16878, -0.11004,
     *  -18.39463,-16.11516,-19.55804, -0.19834, -0.23271, -0.25699,
     *  -19.93482,-17.56433,-18.58818,  0.06508, -0.18075,  0.02796,
     *  -23.64078,-18.77269,-22.77715, -0.02456, -0.12238,  0.02959,
     *  -12.44508,-21.06941,-19.36011,  0.02746, -0.16329,  0.19792,
     *  -26.34187,-19.78854,-24.06651, -0.07299, -0.03082, -0.03535,
     *  -10.71667,-26.04401,-16.59048,  0.02850, -0.09680,  0.15143,
     *  -18.40481,-23.37770,-16.31450, -0.03989, -0.00729, -0.01688,
     *   -9.68886,-20.59304,-18.46657,  0.01092, -0.07901,  0.03422,
     *   -0.06685,-19.24590,-29.35494,  0.12265, -0.24792,  0.05978,
     *  -15.32341, -9.07320,-13.76101, -0.17018, -0.15122, -0.06144,
     *  -14.68939,-14.82251,-13.65846, -0.11173, -0.14410, -0.07133,
     *  -18.38628,-18.94631,-19.00893, -0.08062, -0.14481, -0.12949,
     *  -16.15328,-17.40999,-14.08705, -0.08485, -0.06896, -0.11583,
     *  -14.50295,-16.91671,-25.25793, -0.06814, -0.13727, -0.12213,
     *  -10.92188,-14.10852,-24.43877, -0.09375, -0.11638, -0.09053,
     *  -11.64716,-14.92020,-19.99063, -0.14792, -0.08681, -0.12085,
     *  -24.09766,-16.14519, -8.05683, -0.24065, -0.05877, -0.23726,
     *  -25.18396,-15.02034,-15.50531, -0.12236, -0.09610, -0.00529,
     *  -15.27905,-19.36708,-12.94046, -0.08571, -0.09560, -0.03544,
     *   -7.48927,-16.00753,-13.02842, -0.07862, -0.10110, -0.05807/
        data coeff2/
     *  -13.06383,-27.98698,-18.80004, -0.05875, -0.03737, -0.11214,
     *  -13.67370,-16.44925,-16.12632, -0.07228, -0.09322, -0.05652,
     *  -22.61245,-21.24717,-18.09933, -0.05197, -0.07477, -0.05235,
     *  -27.09189,-21.85181,-20.34676, -0.05123, -0.05683, -0.07214,
     *  -27.09561,-22.76383,-25.41151, -0.10272, -0.02058, -0.16720/

        do i=1,594 
        	coeff(i)=coeff1(i)
        	enddo
        do i=1,30 
        	coeff(i+594)=coeff2(i)
        	enddo

        call g(param,funct,xl)

        y=0.

        do i=1,index_t
          do il=1,index_l
            kk=index_l*(i-1)+il
            do j=1,nfunc
               ind=nfunc*(kk-1)+j
               bspl4=bspl4_time(i,xt)*bspl4_long(il,xl)
               y=y+bspl4*funct(j)*coeff(ind)
            end do
          end do
        end do

       end
C
C
c        real*8 function bspl4_time(i,x1)
        real function bspl4_time(i,x1)
c       *************************************************
c        implicit REAL*8 (A-H,O-Z)
        implicit REAL (A-H,O-Z)
		 
        integer i,order,j,k
c        real*8 t_t(0:39)
c        real*8 x,b(20,20),x1
        real t_t(0:39)
        real x,b(20,20),x1

        data t_t/
     *          0.00,2.75,4.75,5.50,6.25,
     *          7.25,10.00,14.00,17.25,18.00,
     *          18.75,19.75,21.00,24.00,26.75,
     *          28.75,29.50,30.25,31.25,34.00,
     *          38.00,41.25,42.00,42.75,43.75,
     *          45.00,48.00,50.75,52.75,53.50,
     *          54.25,55.25,58.00,62.00,65.25,
     *          66.00,66.75,67.75,69.00,72.00/

        order=4
        x=x1
        if(i.ge.0) then
          if (x.lt.t_t(i-0)) then
              x=x+24
              end if
          end if
        do j=i,i+order-1
           if(x.ge.t_t(j).and.x.lt.t_t(j+1)) then
               b(j,1)=1
           else
	           b(j,1)=0
           end if
        end do

        do j=2,order
          do k=i,i+order-j
            b(k,j)=(x-t_t(k))/(t_t(k+j-1)-t_t(k))*b(k,j-1)
            b(k,j)=b(k,j)+(t_t(k+j)-x)/(t_t(k+j)-t_t(k+1))*
     &                b(k+1,j-1)
          end do
        end do

        bspl4_time=b(i,order)
        end
C
C
        real function bspl4_long(i,x1)
c        real*8 function bspl4_long(i,x1)
c       *************************************************
c       implicit real*8 (A-H,O-Z) 
       implicit real (A-H,O-Z) 

        integer i,order,j,k
c        real*8 t_l(0:24)
c        real*8 x,b(20,20),x1
        real t_l(0:24)
        real x,b(20,20),x1

        data t_l/
     *          0,10,100,190,200,250,280,310,
     *          360,370,460,550,560,610,640,670,
     *          720,730,820,910,920,970,1000,1030,1080/
       
        order=4
        x=x1
        if(i.ge.0) then
          if (x.lt.t_l(i-0)) then
              x=x+360
              end if
          end if
        do j=i,i+order-1
           if(x.ge.t_l(j).and.x.lt.t_l(j+1)) then
              b(j,1)=1
           else
              b(j,1)=0
           end if
        end do

        do j=2,order
          do k=i,i+order-j
             b(k,j)=(x-t_l(k))/(t_l(k+j-1)-t_l(k))*b(k,j-1)
             b(k,j)=b(k,j)+(t_l(k+j)-x)/(t_l(k+j)-t_l(k+1))*
     &                 b(k+1,j-1)
          end do
        end do

        bspl4_long=b(i,order)
        end
C
C
        subroutine g(param,funct,x)
c       *************************************************
c        implicit real*8 (A-H,O-Z)
        implicit real (A-H,O-Z)

        integer i
c        real*8 param(2),funct(6)
c        real*8 x,a,sigma,gauss,flux,cflux
        real param(2),funct(6)
        real x,a,sigma,gauss,flux,cflux

c       *************************************************
        flux=param(2)
        if(param(2).le.75)  flux=75.
        if(param(2).ge.230) flux=230.
        cflux=flux

        a=0.
        if((param(1).ge.120).and.(param(1).le.240)) a=170.
        if((param(1).ge.120).and.(param(1).le.240)) sigma=60
        if((param(1).le.60).or.(param(1).ge.300)) a=170.
        if((param(1).le.60).or.(param(1).ge.300)) sigma=40

        if((flux.le.95).and.(a.ne.0)) then
           gauss=exp(-0.5*((x-a)**2)/sigma**2)
           cflux=gauss*95.+(1-gauss)*flux
           end if
c       *************************************************

c       *************************************************
        do i=1,6
         funct(i)=0.
        end do
c       *************************************************

c       *************************************************
        if((param(1).ge.135).and.(param(1).le.230)) funct(1)=1
        if((param(1).le.45).or.(param(1).ge.320)) funct(2)=1
        if((param(1).gt.75).and.(param(1).lt.105)) funct(3)=1
        if((param(1).gt.260).and.(param(1).lt.290)) funct(3)=1
c       *************************************************

        if((param(1).ge.45).and.(param(1).le.75)) then  ! W-E
            funct(2)=1.-(param(1)-45.)/30.
            funct(3)=1-funct(2)
            end if
        if((param(1).ge.105).and.(param(1).le.135)) then  ! E-S
            funct(3)=1.-(param(1)-105.)/30.
            funct(1)=1-funct(3)
            end if
        if((param(1).ge.230).and.(param(1).le.260)) then  ! S-E
            funct(1)=1.-(param(1)-230.)/30.
            funct(3)=1-funct(1)
            end if
        if((param(1).ge.290).and.(param(1).le.320)) then  ! E-W
            funct(3)=1.-(param(1)-290.)/30.
            funct(2)=1-funct(3)
            end if

c       *************************************************
        funct(4)=(cflux-140)*funct(1)
        funct(5)=(cflux-140)*funct(2)
        funct(6)=(flux-140)*funct(3)
c       *************************************************

        end
c
c
       SUBROUTINE StormVd(FLAG,iP,AE,SLT,PromptVd,DynamoVd,Vd)
C *******************************************************************
C  Empirical vertical disturbance drifts model
C  After Fejer and Scherliess, JGR, 102, 24047-24056,1997
C*********************************************************************
C  INPUT:
C    AE: AE(in nT) in 1 hour or 15 minute resolution;
C    SLT: Local time(in hrs) for wanted Vd;
C  OUTPUT:
C    PromptVd: Prompt penetration vertical drifts at given conditions;
C    DynamoVd: Disturbane dynamo vertical drifts at given conditions;
C    Vd: PromptVd+DynamoVd;
C*********************************************************************

c       IMPLICIT REAL*8(A-H,O-Z)
       IMPLICIT REAL(A-H,O-Z)
c       REAL*8 AE(1:366*24*4),Coff1(1:5,1:9),Coff15(1:6,1:9)
       REAL AE(1:366*24*4),Coff1(1:5,1:9),Coff15(1:6,1:9)
       INTEGER FLAG 
       DATA Coff1/
     @           0.0124,-0.0168,-0.0152,-0.0174,-0.0704,
     @          -0.0090,-0.0022,-0.0107, 0.0152,-0.0674,
     @           0.0275, 0.0051,-0.0132, 0.0020,-0.0110,
     @          -0.0022, 0.0044, 0.0095, 0.0036,-0.0206,
     @           0.0162, 0.0007, 0.0085,-0.0140, 0.0583,
     @           0.0181, 0.0185,-0.0109,-0.0031,-0.0427,
     @          -0.0057, 0.0002, 0.0086, 0.0149, 0.2637,
     @          -0.0193, 0.0035, 0.0117, 0.0099, 0.3002,
     @          -0.0492,-0.0201, 0.0338, 0.0099, 0.0746/

	 DATA Coff15/
     @	        0.0177, 0.0118,-0.0006,-0.0152,-0.0174,-0.0704,
     @	        0.0051,-0.0074,-0.0096,-0.0107, 0.0152,-0.0674,
     @	        0.0241, 0.0183, 0.0122,-0.0132, 0.0020,-0.0110,
     @	        0.0019,-0.0010, 0.0001, 0.0095, 0.0036,-0.0206,
     @	        0.0170, 0.0183, 0.0042, 0.0085,-0.0140, 0.0583,
     @          0.0086, 0.0189, 0.0200,-0.0109,-0.0031,-0.0427,
     @	       -0.0070,-0.0053,-0.0090, 0.0086, 0.0149, 0.2637,
     @	       -0.0326,-0.0101, 0.0076, 0.0117, 0.0099, 0.3002,
     @	       -0.0470,-0.0455,-0.0274, 0.0338, 0.0099, 0.0746/

CCCCCCCCCCCCCCCCC**Define to variables**CCCCCCCCCCCCCCCCCCCCC
C To 1 h time resolution:
C dAEt_30=AE(t)-AE(t-1 hour);
C dAEt_90=AE(t-1 hour)-AE(t-2 hour);
CC
C To 15 MIN time resolution :
C dAEt_7P5=AE(t)-AE(t-15min);
C dAEt_30=AE(t-15)-AE(t-45min);
C dAEt_75=AE(t-45)-AE(t-105min);
CC
C  Following variables are the same to two resolution: 
C AE1_6=average(AE(1-6hours));
C AE7_12=average(AE(7-12hours));
C AE1_12=average(AE(1-12hours));
C AEd1_6=average(X(AE(1-6hours)-130 nT));
C AEd7_12=average(X(AE(7-12hours)-130 nT));
C AEd1_12=average(X(AE(1-12hours)-130 nT));
C AEd22_28=average(X(AE(22-28hours)-130 nT));
C Here X(a)=a, a>0; =0, a<=0;
C Alfa=0,            AE1_6<200 nT;
C      AE1_6/100-2, 200 nT<AE1_6<200 nT;
C      1,            AE1_6>300 nT;
C Beta=exp(-AE1_12/90),  AE1_12>=70nT;
C      0.46,              AE1_12<70 nT;
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccccccc
C*****************************************************
CC        FLAG>0--> 1 h time resolution
C**************************************************** 

       IF (FLAG.GT.0) THEN
C 
         dAEt_30=AE(iP)-AE(iP-1)
         dAEt_90=AE(iP-1)-AE(iP-2)
C
         AE1_6=0.0D0
         AEd1_6=0.0D0

         DO i=-1,-6,-1
             AE1_6=AE1_6+AE(iP+i)
             AEd1_6S=AE(iP+i)-130.0D0
             IF (AEd1_6S.LE.0.0D0) AEd1_6S=0.0D0
             AEd1_6=AEd1_6+AEd1_6S
         END DO

         AE1_6=AE1_6/6.0D0
         AEd1_6=AEd1_6/6.0D0
C
         AEd7_12=0.0D0

         DO i=-7,-12,-1
            AEd7_12S=AE(iP+i)-130.0D0
            IF (AEd7_12S.LE.0.0D0) AE7_12S=0.0D0
            AEd7_12=AEd7_12+AEd7_12S
         END DO

         AEd7_12=AEd7_12/6.0D0
C
         AE1_12=0.0D0

         DO i=-1,-12,-1
            AE1_12=AE1_12+AE(iP+i)
         END DO

         AE1_12=AE1_12/12.0D0
C
         AEd22_28=0.0D0

         DO i=-22,-28,-1
            AEd22_28S=AE(iP+i)-130.0D0
            IF (AED22_28S.LE.0.0D0) AEd22_28S=0.0D0
            AEd22_28=AEd22_28+AEd22_28S 
         END DO

         AEd22_28=AEd22_28/7.0D0
         AEd22_28P=AEd22_28-200.0D0
         IF (AEd22_28P.LE.0.0D0) AEd22_28P=0.0D0
CC
         IF (AE1_6.GT.300.0D0) THEN  
           Alfa=1.0D0
         ELSE IF (AE1_6.GT.200.0D0) THEN
           ALfa=AE1_6/100.0D0-2.0D0
         ELSE 
           ALfa=0.0D0
         ENDIF
CC
         IF (AE1_12.GE.70.0D0) THEN
            Beta=dexp(-AE1_12/90.0D0)
         ELSE
            Beta=0.46D0
         END IF
         PromptVd=0.0D0
         DO J=1,9
            PromptVd=PromptVd +(Coff1(1,J)*dAEt_30 +Coff1(2,J)*dAEt_90
     #                         )*bspl4_ptime(J,SLT)
         END DO
         DynamoVd=0.0D0
         DO J=1,9
            DynamoVd=DynamoVd+
     #               (Coff1(3,J)*AEd1_6+Coff1(4,J)*Alfa*AEd7_12
     #                +Coff1(5,J)*Beta*AEd22_28P)*bspl4_ptime(J,SLT)
         END DO
         Vd=PromptVd+DynamoVd
         RETURN

C 1 h time resolution end;
C********************************************************************
C                  15 min time resolution
C********************************************************************

       ELSE
         dAEt_7P5=AE(iP)-AE(iP-1)
         dAEt_30=AE(iP-1)-AE(iP-3)
         dAEt_75=AE(iP-3)-AE(iP-7)

         AE1_6=0.0D0
         AEd1_6=0.0D0

         DO i=-4,-24,-1
           AE1_6=AE1_6+AE(iP+i)
           AEd1_6s=AE(iP+i)-130.
           IF (AEd1_6s.LE.0.0) AEd1_6s=0.0
           AEd1_6=AEd1_6+AEd1_6S
         ENDDO

         AE1_6=AE1_6/21.0D0
         AEd1_6=AEd1_6/21.0D0
CC
         AEd7_12=0.0D0 
         DO i=-28,-48,-1
            AEd7_12s=AE(iP+i)-130.0
            IF (AEd7_12s.LE.0) AEd7_12s=0.0
            AEd7_12=AEd7_12+AEd7_12S
         ENDDO
         AEd7_12=AEd7_12/21.0D0
CC
         AE1_12=0.0D0
         DO i=-4,-48,-1
            AE1_12=AE1_12+AE(iP+i)
         END DO
         AE1_12=AE1_12/45.0D0
CC 
	  AEd22_28=0.0D0
          DO i=-88,-112,-1
             AEd22_28s=AE(iP+i)-130.
             IF (AEd22_28s.LE.0) AEd22_28s=0.0
	     AEd22_28=AEd22_28+AEd22_28s
	  ENDDO
          AEd22_28=AEd22_28/25.0D0
	  AEd22_28P=AEd22_28-200.0D0
	  IF (AEd22_28P.LE.0.0D0) AEd22_28P=0.0D0

c         AE1_6=0.0D0
c         AEd1_6=0.0D0
c         AEd7_12=0.0D0 
c         AEd22_28P=0.0D0
c         AE1_12=0.0D0
c         dAEt_7P5=400.D0
c         dAEt_30=0.D0
c         dAEt_75=0.D0
CC
  	  IF (AE1_6.GT.300.0D0) THEN  
             Alfa=1.0D0
          ELSE IF (AE1_6.GT.200.0D0) THEN
             ALfa=AE1_6/100.0D0-2.0D0
          ELSE 
             ALfa=0.0D0
          ENDIF
CC
          IF (AE1_12.GE.70.0D0) THEN
             Beta=dexp(-AE1_12/90.0D0)
          ELSE
             Beta=0.46D0
          END IF
CC
          PromptVd=0.0D0
          DO J=1,9
            PromptVd=PromptVd+(Coff15(1,J)*dAEt_7P5+Coff15(2,J)*dAEt_30
     #                      +Coff15(3,J)*dAEt_75)*bspl4_ptime(J,SLT)
          END DO
          DynamoVd=0.0D0
          DO J=1,9
             DynamoVd=DynamoVd +(Coff15(4,J)*AEd1_6+
     #                           Coff15(5,J)*Alfa*AEd7_12+
     #                           Coff15(6,J)*Beta*AEd22_28P
     #                          )*bspl4_ptime(J,SLT)
          END DO
          Vd=PromptVd+DynamoVd
       ENDIF
       RETURN
       END                     
C
C
       real function bspl4_ptime(i,x1)
c       real*8 function bspl4_ptime(i,x1)
C *************************************************

c       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT REAL (A-H,O-Z)

       integer i,order,j,k
c       real*8 t_t(0:27)
c       real*8 x,b(20,20),x1
       real t_t(0:27)
       real x,b(20,20),x1

       data t_t/0.00,3.00,4.50,6.00,9.00,12.0,15.0,18.0,21.0,
     *          24.0,27.0,28.5,30.0,33.0,36.0,39.0,42.0,45.0,
     *          48.0,51.0,52.5,54.0,57.0,60.0,63.0,66.0,69.0,72.0/
C
       order=4
       x=x1
       if(i.ge.0) then
          if (x.lt.t_t(i-0)) then
	       x=x+24
          end if
       end if
       do j=i,i+order-1
	   if(x.ge.t_t(j).and.x.lt.t_t(j+1)) then
	       b(j,1)=1
	   else
	       b(j,1)=0
	   end if
       end do
c
       do j=2,order
          do k=i,i+order-j
             b(k,j)=(x-t_t(k))/(t_t(k+j-1)-t_t(k))*b(k,j-1)
             b(k,j)=b(k,j)+(t_t(k+j)-x)/(t_t(k+j)-t_t(k+1))*b(k+1,j-1)
          end do
       end do
       bspl4_ptime=b(i,order)
       return
       end

C
C***************************************************************************
C

       subroutine spreadf_brazil(idoy,idiy,f107,geolat,osfbr)
**********************************************************************       
*
*       SUBROUTINE CALCULATES PERCENTAGE OF SPREAD F OCCURRENCE OVER 
*       BRAZILIAN SECTOR AS DESCRIBED IN:
*       ABDU ET AL., Advances in Space Research, 31(3), 
*       703-716, 2003
*
*    INPUT:
*         IDOY: DAY OF YEAR (1 TO 365/366)
*         IDIY: DAYS IN YEAR (365 OR 366)
*         F107: F10.7 cm SOLAR FLUX (DAILY VALUE)
*         GEOLAT: BRAZILIAN GEOGRAPHIC LATITUDE BETWEEN -4 AND -22.5
*
*    OUTPUT:         
*         OSFBR(25): PERCENTAGE OF SPREAD F OCCURRENCE FOR 25 TIME 
*                    STEPS FROM LT=18 TO LT=7 ON THE NEXT DAY IN
*                    STEPS OF 0.5 HOURS.
*
**********************************************************************
*
         dimension param(3),osfbr(25),coef_sfa(684),coef_sfb(684),
     &             sosf(2,32,3,12)
         common/mflux/kf,n     
         data coef_sfa/
     *   0.07,0.13,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.05,0.04,0.03
     *  ,0.06,0.07,0.02,0.03,0.03,0.07,0.06,0.07,0.21,0.28,0.34,0.16
     *  ,0.12,0.00,0.02,0.02,0.04,0.05,0.02,0.11,0.19,0.31,0.31,0.11
     *  ,0.14,0.16,0.03,0.00,0.00,0.02,0.00,0.00,0.05,0.55,0.61,0.28
     *  ,0.26,0.10,0.15,0.23,0.07,0.06,0.03,0.03,0.41,0.88,0.89,0.65
     *  ,0.19,0.18,0.17,0.10,0.14,0.15,0.03,0.14,0.46,0.72,0.71,0.53
     *  ,0.57,0.38,0.30,0.14,0.00,0.04,0.03,0.02,0.21,0.84,0.87,0.72
     *  ,0.79,0.60,0.65,0.70,0.29,0.19,0.19,0.32,0.73,0.96,0.99,0.84
     *  ,0.75,0.78,0.79,0.70,0.63,0.24,0.28,0.53,0.75,0.77,0.75,0.85
     *  ,0.78,0.51,0.59,0.24,0.00,0.07,0.05,0.06,0.33,0.92,0.96,0.89
     *  ,0.90,0.84,0.86,0.81,0.33,0.27,0.23,0.47,0.90,1.00,1.00,0.96
     *  ,0.96,0.89,0.92,0.84,0.80,0.27,0.35,0.61,0.81,0.93,0.86,0.97
     *  ,0.84,0.65,0.75,0.25,0.00,0.04,0.08,0.06,0.53,0.93,0.96,0.94
     *  ,0.95,0.84,0.91,0.71,0.18,0.17,0.21,0.42,0.92,0.99,0.97,0.92
     *  ,0.92,0.93,0.92,0.67,0.58,0.21,0.38,0.55,0.83,0.90,0.89,0.97
     *  ,0.84,0.71,0.91,0.21,0.02,0.07,0.03,0.03,0.60,0.95,0.96,0.92
     *  ,0.97,0.91,0.92,0.67,0.11,0.08,0.09,0.23,0.90,0.99,0.99,0.96
     *  ,0.96,0.93,0.98,0.63,0.25,0.08,0.12,0.41,0.79,0.95,0.98,0.99
     *  ,0.86,0.80,0.94,0.22,0.02,0.04,0.03,0.03,0.63,0.95,0.96,0.94
     *  ,0.98,0.90,0.91,0.59,0.10,0.04,0.07,0.15,0.83,0.97,0.97,0.90
     *  ,0.92,0.93,0.95,0.57,0.12,0.03,0.05,0.23,0.74,0.94,0.94,0.99
     *  ,0.84,0.84,0.90,0.24,0.02,0.07,0.07,0.03,0.60,0.95,0.96,0.97
     *  ,0.93,0.82,0.83,0.51,0.08,0.07,0.09,0.09,0.71,0.95,0.92,0.87
     *  ,0.91,0.91,0.89,0.50,0.14,0.03,0.06,0.14,0.61,0.84,0.89,0.94
     *  ,0.77,0.82,0.84,0.34,0.10,0.11,0.12,0.06,0.43,0.87,0.94,0.97
     *  ,0.91,0.77,0.68,0.42,0.06,0.08,0.10,0.04,0.51,0.78,0.71,0.77
     *  ,0.85,0.88,0.77,0.35,0.16,0.05,0.08,0.15,0.53,0.70,0.60,0.89
     *  ,0.85,0.71,0.72,0.26,0.16,0.17,0.08,0.15,0.38,0.73,0.91,0.91
     *  ,0.89,0.68,0.53,0.26,0.06,0.12,0.08,0.09,0.32,0.63,0.67,0.77
     *  ,0.81,0.79,0.59,0.21,0.14,0.03,0.06,0.09,0.23,0.51,0.34,0.79
     *  ,0.88,0.66,0.59,0.16,0.18,0.15,0.16,0.16,0.33,0.67,0.75,0.88
     *  ,0.80,0.64,0.52,0.16,0.04,0.09,0.04,0.09,0.24,0.47,0.53,0.50
     *  ,0.73,0.69,0.48,0.11,0.14,0.03,0.03,0.03,0.20,0.37,0.28,0.54
     *  ,0.81,0.64,0.49,0.18,0.12,0.17,0.16,0.19,0.31,0.57,0.70,0.83
     *  ,0.76,0.57,0.52,0.13,0.04,0.06,0.05,0.08,0.21,0.49,0.47,0.39
     *  ,0.69,0.66,0.43,0.11,0.10,0.02,0.00,0.03,0.16,0.39,0.24,0.35
     *  ,0.77,0.45,0.39,0.10,0.10,0.13,0.15,0.18,0.29,0.57,0.70,0.69
     *  ,0.71,0.49,0.54,0.20,0.05,0.06,0.05,0.06,0.27,0.42,0.36,0.42
     *  ,0.61,0.59,0.50,0.08,0.06,0.02,0.03,0.02,0.16,0.40,0.17,0.31
     *  ,0.68,0.30,0.28,0.13,0.10,0.16,0.14,0.08,0.19,0.50,0.63,0.62
     *  ,0.63,0.45,0.51,0.13,0.06,0.07,0.04,0.06,0.27,0.42,0.28,0.35
     *  ,0.68,0.53,0.57,0.15,0.05,0.00,0.00,0.05,0.31,0.33,0.18,0.22
     *  ,0.59,0.32,0.21,0.06,0.10,0.16,0.12,0.10,0.19,0.41,0.55,0.54
     *  ,0.69,0.43,0.43,0.15,0.06,0.05,0.05,0.08,0.29,0.39,0.23,0.29
     *  ,0.57,0.51,0.56,0.13,0.06,0.00,0.00,0.05,0.34,0.27,0.19,0.24
     *  ,0.49,0.16,0.13,0.09,0.04,0.11,0.11,0.05,0.17,0.32,0.49,0.49
     *  ,0.60,0.42,0.38,0.11,0.06,0.04,0.07,0.07,0.25,0.36,0.21,0.25
     *  ,0.65,0.48,0.53,0.17,0.05,0.00,0.00,0.11,0.29,0.14,0.20,0.22
     *  ,0.44,0.16,0.18,0.07,0.04,0.04,0.07,0.03,0.12,0.23,0.39,0.43
     *  ,0.57,0.40,0.35,0.14,0.06,0.03,0.04,0.07,0.18,0.27,0.14,0.15
     *  ,0.45,0.50,0.50,0.19,0.06,0.00,0.02,0.05,0.26,0.19,0.15,0.18
     *  ,0.23,0.09,0.12,0.06,0.04,0.02,0.02,0.02,0.10,0.03,0.14,0.26
     *  ,0.39,0.34,0.22,0.07,0.03,0.00,0.04,0.01,0.15,0.01,0.04,0.14
     *  ,0.41,0.39,0.35,0.13,0.02,0.00,0.00,0.06,0.17,0.07,0.06,0.14
     *  ,0.07,0.02,0.03,0.00,0.00,0.00,0.00,0.00,0.00,0.01,0.03,0.08
     *  ,0.19,0.14,0.14,0.00,0.03,0.01,0.02,0.00,0.09,0.00,0.01,0.00
     *  ,0.18,0.09,0.16,0.08,0.01,0.00,0.02,0.02,0.15,0.00,0.03,0.04/
*
        data coef_sfb/
     *   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
     *  ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.03,0.00,0.00,0.00,0.00
     *  ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.02,0.01,0.00,0.00,0.00
     *  ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
     *  ,0.00,0.01,0.00,0.00,0.00,0.00,0.00,0.01,0.01,0.00,0.00,0.00
     *  ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.05,0.03,0.00,0.02,0.00
     *  ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
     *  ,0.00,0.04,0.00,0.01,0.00,0.00,0.00,0.01,0.01,0.05,0.00,0.00
     *  ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.04,0.00,0.03,0.03,0.00
     *  ,0.01,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.01,0.00
     *  ,0.01,0.04,0.04,0.03,0.00,0.01,0.00,0.01,0.00,0.27,0.14,0.06
     *  ,0.05,0.04,0.02,0.00,0.00,0.00,0.00,0.04,0.09,0.48,0.43,0.27
     *  ,0.05,0.04,0.01,0.00,0.00,0.00,0.00,0.00,0.00,0.13,0.16,0.06
     *  ,0.26,0.12,0.29,0.04,0.01,0.02,0.00,0.01,0.08,0.65,0.56,0.45
     *  ,0.43,0.42,0.42,0.09,0.00,0.02,0.00,0.00,0.34,0.67,0.73,0.72
     *  ,0.10,0.05,0.04,0.00,0.01,0.00,0.00,0.00,0.00,0.18,0.39,0.15
     *  ,0.61,0.37,0.51,0.06,0.01,0.02,0.01,0.01,0.18,0.72,0.63,0.80
     *  ,0.77,0.66,0.70,0.19,0.00,0.02,0.02,0.02,0.41,0.68,0.88,0.85
     *  ,0.24,0.11,0.08,0.00,0.01,0.00,0.00,0.00,0.00,0.28,0.51,0.29
     *  ,0.75,0.48,0.57,0.11,0.00,0.02,0.01,0.01,0.19,0.77,0.77,0.88
     *  ,0.89,0.81,0.74,0.21,0.02,0.02,0.02,0.02,0.42,0.71,0.93,0.95
     *  ,0.49,0.30,0.19,0.00,0.00,0.00,0.00,0.01,0.06,0.38,0.64,0.48
     *  ,0.86,0.60,0.62,0.12,0.00,0.02,0.01,0.00,0.18,0.81,0.84,0.94
     *  ,0.88,0.79,0.70,0.26,0.03,0.02,0.02,0.02,0.36,0.61,0.98,0.93
     *  ,0.60,0.46,0.31,0.03,0.00,0.01,0.00,0.00,0.09,0.50,0.71,0.58
     *  ,0.90,0.65,0.66,0.10,0.00,0.02,0.01,0.02,0.15,0.69,0.80,0.86
     *  ,0.84,0.75,0.64,0.09,0.03,0.00,0.00,0.04,0.26,0.54,0.78,0.92
     *  ,0.62,0.59,0.44,0.01,0.00,0.01,0.00,0.00,0.13,0.52,0.77,0.63
     *  ,0.84,0.67,0.63,0.11,0.00,0.00,0.03,0.03,0.18,0.65,0.75,0.84
     *  ,0.81,0.63,0.47,0.06,0.02,0.00,0.00,0.05,0.14,0.49,0.76,0.91
     *  ,0.58,0.63,0.47,0.09,0.00,0.07,0.01,0.04,0.15,0.48,0.68,0.61                     
     *  ,0.79,0.63,0.55,0.12,0.01,0.01,0.02,0.05,0.13,0.57,0.51,0.63
     *  ,0.72,0.54,0.43,0.11,0.02,0.00,0.00,0.09,0.16,0.39,0.59,0.72
     *  ,0.46,0.55,0.39,0.07,0.01,0.03,0.03,0.06,0.15,0.37,0.51,0.50
     *  ,0.61,0.43,0.38,0.11,0.01,0.03,0.02,0.03,0.10,0.38,0.38,0.60
     *  ,0.58,0.42,0.38,0.15,0.02,0.00,0.00,0.11,0.13,0.24,0.41,0.51
     *  ,0.36,0.36,0.21,0.04,0.04,0.03,0.06,0.05,0.06,0.26,0.39,0.43
     *  ,0.43,0.31,0.24,0.09,0.02,0.00,0.02,0.02,0.06,0.24,0.24,0.40
     *  ,0.53,0.19,0.28,0.13,0.02,0.02,0.02,0.09,0.13,0.17,0.24,0.40
     *  ,0.32,0.27,0.17,0.03,0.04,0.02,0.04,0.03,0.06,0.13,0.34,0.36
     *  ,0.42,0.31,0.20,0.09,0.03,0.00,0.02,0.01,0.07,0.19,0.24,0.32
     *  ,0.44,0.10,0.23,0.13,0.03,0.02,0.00,0.09,0.12,0.17,0.21,0.33
     *  ,0.32,0.23,0.16,0.00,0.02,0.04,0.03,0.03,0.06,0.15,0.29,0.34
     *  ,0.36,0.26,0.28,0.07,0.01,0.00,0.01,0.02,0.04,0.19,0.17,0.27
     *  ,0.34,0.14,0.26,0.09,0.03,0.02,0.00,0.06,0.13,0.09,0.16,0.22
     *  ,0.29,0.21,0.15,0.00,0.02,0.02,0.02,0.03,0.11,0.16,0.26,0.28
     *  ,0.29,0.22,0.27,0.05,0.01,0.00,0.01,0.01,0.02,0.14,0.09,0.19
     *  ,0.25,0.19,0.25,0.07,0.02,0.02,0.00,0.00,0.09,0.07,0.12,0.15
     *  ,0.23,0.20,0.16,0.00,0.03,0.04,0.00,0.00,0.08,0.09,0.21,0.18
     *  ,0.22,0.21,0.19,0.02,0.02,0.00,0.01,0.03,0.04,0.08,0.06,0.14
     *  ,0.20,0.12,0.23,0.02,0.00,0.02,0.00,0.00,0.05,0.05,0.09,0.11
     *  ,0.14,0.16,0.13,0.00,0.03,0.04,0.00,0.00,0.05,0.05,0.04,0.09
     *  ,0.09,0.13,0.16,0.03,0.01,0.00,0.01,0.03,0.01,0.03,0.04,0.10
     *  ,0.14,0.09,0.17,0.02,0.02,0.00,0.00,0.02,0.04,0.04,0.03,0.07
     *  ,0.00,0.11,0.09,0.00,0.02,0.00,0.00,0.00,0.01,0.00,0.02,0.02
     *  ,0.02,0.06,0.11,0.00,0.00,0.00,0.00,0.01,0.00,0.00,0.01,0.02
     *  ,0.06,0.09,0.13,0.00,0.02,0.00,0.03,0.02,0.03,0.01,0.02,0.01/
*
        param(1)=idoy
        param(2)=f107
        param(3)=geolat
        n=idiy-365
*
        if(param(1).le.31.)kf=1
         if(param(1).gt.31..and.param(1).le.(59+n))kf=2
	      if(param(1).gt.(59+n).and.param(1).le.(90+n))kf=3
           if(param(1).gt.(90+n).and.param(1).le.(120+n))kf=4
            if(param(1).gt.(120+n).and.param(1).le.(151+n))kf=5
	         if(param(1).gt.(151+n).and.param(1).le.(181+n))kf=6
             if(param(1).gt.(181+n).and.param(1).le.(212+n))kf=7
            if(param(1).gt.(212+n).and.param(1).le.(243+n))kf=8
           if(param(1).gt.(243+n).and.param(1).le.(273+n))kf=9
	      if(param(1).gt.(273+n).and.param(1).le.(304+n))kf=10
         if(param(1).gt.(304+n).and.param(1).le.(334+n))kf=11
        if(param(1).gt.(334+n).and.param(1).le.(365+n))kf=12
*
      do i=1,32
       do j=1,3
        do k=1,12
           sosf(1,i,j,k)=0.
           sosf(2,i,j,k)=0.           
        enddo
       enddo
      enddo
*         
      kc=0
      do i=5,23
       do j=1,3
        do k=1,12
           kc=kc+1
           sosf(1,i,j,k)=coef_sfa(kc)
           sosf(2,i,j,k)=coef_sfb(kc) 
        enddo
       enddo
      enddo
           
      kk=0    
      do it=1600,3200,50      
        slt=it/100.
        osft=0.
        do i=1,23
          il=i+3
	      if(il.gt.23)il=il-23
          do j=1,12
	        jl=j+2
	        if(jl.gt.12)jl=jl-12
            do m=1,3
	          ml=m+1
	          if(ml.gt.3)ml=ml-3
              do  l=1,2
                bspl4=bspl4t(i,slt)*bspl2s(j,param(1))*
     &	         bspl2l(l,param(3))*bspl2f(m,param(2))
                osft=osft+bspl4*sosf(l,il,ml,jl)           
              enddo
            enddo
          enddo
        enddo
        if(slt.gt.17.98.and.slt.lt.30.01)then
          kk=kk+1
          osfbr(kk)=osft 
          endif
      enddo
*
*
      do iii=1,25 
         if(osfbr(iii).gt.1.) osfbr(iii)=1.
         if(osfbr(iii).lt.0.) osfbr(iii)=0.
         enddo
      return
      end
*
**********************************************************************
      function bspl4t(i,t1)
**********************************************************************
      dimension tt(0:78),b(30,30)
*
      data tt/16.00,16.50,17.00,17.50,18.00,18.50,19.00,19.50,20.00,
     &  20.50,21.00,22.00,23.00,24.00,25.00,26.00,27.00,27.50,28.00,
     &  28.50,29.00,29.50,30.00,30.50,31.00,32.00,40.00,40.50,41.00,
     &  41.50,42.00,42.50,43.00,43.50,44.00,44.50,45.00,46.00,47.00,
     &  48.00,49.00,50.00,51.00,51.50,52.00,52.50,53.00,53.50,54.00,
     &  54.50,55.00,56.00,64.00,64.50,65.00,65.50,66.00,66.50,67.00,
     &  67.50,68.00,68.50,69.00,70.00,71.00,72.00,73.00,74.00,75.00,
     &  75.50,76.00,76.50,77.00,77.50,78.00,78.50,79.00,80.00,88.00/
*
      t=t1
      if(i.ge.0.and.t.lt.tt(i)) then
         t=t+24.
      endif
      do j=i,i+4-1
        if(t.ge.tt(j).and.t.lt.tt(j+1)) then
           b(j,1)=1.
        else
           b(j,1)=0.
        endif
      enddo
      do j=2,4
        do k=i,i+4-j
          b(k,j)=(t-tt(k))/(tt(k+j-1)-tt(k))*b(k,j-1)
          b(k,j)=b(k,j)+(tt(k+j)-t)/(tt(k+j)-tt(k+1))*b(k+1,j-1)
        enddo
      enddo
*
      bspl4t=b(i,4)
*
      return
*
      end
*
******************************************************************
      function bspl2s(i,t1)
******************************************************************
      dimension ts(0:36),b(30,30)
*
      data ts/ 15,46,74,105,135,166,196,227,258,288,319,349,
     *        380,411,439,470,500,531,561,592,623,653,684,714,
     *        745,776,804,835,865,896,926,957,988,1018,1049,
     *        1079,1110/
*
      t=t1
      if(i.ge.0.and.t.lt.ts(i)) then
         t=t+365.
      endif
      do j=i,i+2-1
        if(t.ge.ts(j).and.t.lt.ts(j+1)) then
           b(j,1)=1.
        else
           b(j,1)=0.
        endif
      enddo
*
      do j=2,2
        do k=i,i+4-j
          b(k,j)=(t-ts(k))/(ts(k+j-1)-ts(k))*b(k,j-1)
          b(k,j)=b(k,j)+(ts(k+j)-t)/(ts(k+j)-ts(k+1))*b(k+1,j-1)
       enddo
      enddo
*
      bspl2s=b(i,2)
      return
      end
*
*******************************************************************
      function bspl2l(i,t1)
*******************************************************************
      dimension ts(0:6),b(30,30)
*
      data ts/ 94.,112.5,454.,472.5,814.,832.5,1174./
*
      t=t1
      if(i.ge.0.and.t.lt.ts(i)) then
         t=t+360.
      endif
      do j=i,i+2-1
        if(t.ge.ts(j).and.t.lt.ts(j+1)) then
           b(j,1)=1.
        else
           b(j,1)=0.
        endif
      enddo
*
      do j=2,2
        do k=i,i+2-j
          b(k,j)=(t-ts(k))/(ts(k+j-1)-ts(k))*b(k,j-1)
          b(k,j)=b(k,j)+(ts(k+j)-t)/(ts(k+j)-ts(k+1))*b(k+1,j-1)
        enddo
      enddo
*
      bspl2l=b(i,2)
*
      return
*
      end
*
*************************************************************************
      function bspl2f(i,t1)
*************************************************************************
       dimension ts(0:9),b(30,30),
     & ifnodes1(12),ifnodes2(12),ifnodes3(12)
       common/mflux/kf,n
*
      data ifnodes1 / 78, 77, 75, 79, 80, 77, 78, 80, 76, 81, 78, 78/
      data ifnodes2 /144,140,139,142,139,146,142,139,150,151,150,157/
      data ifnodes3 /214,211,201,208,213,220,203,209,213,215,236,221/ 
*
	ts(0)=ifnodes1(kf)
        ts(1)=ifnodes2(kf)
	ts(2)=ifnodes3(kf)
	ts(3)=ts(1)+367
        ts(4)=ts(2)+367
	ts(5)=ts(3)+367
	ts(6)=ts(4)+367
        ts(7)=ts(5)+367
	ts(8)=ts(6)+367
        ts(9)=ts(7)+367
*
      t=t1
      if(i.ge.0.and.t.lt.ts(i)) then
         t=t+367.
      endif
      do j=i,i+2-1
        if(t.ge.ts(j).and.t.lt.ts(j+1)) then
           b(j,1)=1.
        else
           b(j,1)=0.
        endif
      enddo
*
      do j=2,2
        do k=i,i+2-j
          b(k,j)=(t-ts(k))/(ts(k+j-1)-ts(k))*b(k,j-1)
          b(k,j)=b(k,j)+(ts(k+j)-t)/(ts(k+j)-ts(k+1))*b(k+1,j-1)
        enddo
      enddo
*
      bspl2f=b(i,2)
      return
      end
C
C
        function ckp(ap)
C-----------------------------------------------------------------------
C Converts ap index (ap is integer variable varying from 0 to 400) into 
C kp index (xkp is real variable varying from 0 to 9). Using standard
C tables for deriving the 3-hourly ap index from the 3-hourly Kp index
C (e.g., http://www.ngdc.noaa.gov/stp/GEOMAG/kp_ap.shtml) 
C-----------------------------------------------------------------------

        integer		ap,ap_array
        real		kp_array,ap_log_array
        dimension 	ap_array(28),kp_array(28),alap(28)
        data ap_array /0,2,3,4,5,6,7,9,12,15,18,22,27,32,39,48,56,67,
     &                 80,94,111,132,154,179,207,236,300,400/
        
        do 1256 i=2,28
1256       kp_array(i)=(i-1)/3.
		
        if(ap.eq.0) then
        	ckp=0.0
        	return
        	endif
        if(ap.eq.1) then
        	ckp=kp_array(2)/2.
        	return
        	endif
        if(ap.lt.8.and.ap.gt.1) then
        	ckp=kp_array(ap)
        	return
        	endif

        xl_ap=log(ap*1.0)
                	
        i=8
1257    alap(i)=log(ap_array(i)*1.0)
        if(xl_ap.gt.alap(i)) then
                i=i+1
                if(i.le.28) goto 1257
                endif

        slope=(kp_array(i)-kp_array(i-1))/(alap(i)-alap(i-1))
       
		ckp = kp_array(i) + slope * (xl_ap - alap(i))
		
		return
		end
C                   
C        
		subroutine auroral_boundary(xkp,xmlt,cgmlat,ab_mlat)
C-----------------------------------------------------------------------
C Computes equatorward auroral boundary values for givern kp value.
C kp given in units of 0.1 (xkp) for the range from 0.0 to 9.0. Model 
C values are only used for kp=0,1,2,3,4,5,6,7,8,9 and a linear inter-
C polation is applied for intermediate kp values.
C 
C The auroral oval boundary is given as an array for corrected magnetic 
C latitude CGM (ab_mlat). The 48 values correspond to the MLT values 
C of 0.0,0.5,1.0,1.5,2.0 .. 23.5. If the input xmlt is greater than
C -1, then the program determines the CGM latitude, cgmlat, that 
C corresponds to the given MLT value (xmlt).
C
C Y. Zhang and L.J. Paxton, An empirical Kp-dependent global auroral 
C model based on TIMED/GUVI FUV data, Journal of Atmospheric and 
C Solar-Terrestrial Physics 70, 1231–1242, 2008.
C 
C----------------------------------------------------------------------- 

        dimension zp_mlat(48,10),ab_mlat(48),ab_mlt(48)

      data zp_mlat/
     * 66.1,65.9,65.9,66.0,66.2,66.5,66.7,67.0,67.4,67.8,68.2,68.7,
     * 69.4,70.0,70.5,70.8,71.0,71.3,71.9,72.7,73.8,75.4,77.0,77.8,
     * 77.3,76.7,76.4,76.2,75.9,75.4,74.7,74.1,73.5,73.0,72.5,71.9,
     * 71.6,71.1,70.6,70.0,69.2,68.5,67.9,67.5,67.2,66.9,66.7,66.4,
     * 63.0,62.9,63.0,63.0,63.1,63.2,63.3,63.6,64.0,64.5,65.0,65.6,
     * 66.2,66.7,67.1,67.6,68.1,68.6,69.4,70.3,71.5,72.8,74.1,74.9,
     * 74.8,74.7,74.7,74.5,74.0,73.1,72.2,71.2,70.4,69.7,68.9,68.2,
     * 67.5,66.9,66.3,65.8,65.4,65.0,64.5,64.1,63.7,63.3,63.1,62.9,
     * 61.1,61.3,61.5,61.7,61.9,62.1,62.2,62.4,62.7,63.0,63.4,64.0,
     * 64.4,65.0,65.5,65.9,66.4,66.9,67.5,68.3,69.3,70.5,71.6,72.3,
     * 72.7,72.8,72.8,72.4,71.6,70.6,69.7,68.9,68.2,67.5,66.7,65.9,
     * 65.0,64.3,63.7,63.4,63.1,62.8,62.4,62.0,61.5,61.2,61.0,61.0,
     * 59.6,60.0,60.4,60.7,60.7,60.7,60.5,60.4,60.6,61.1,61.8,62.5,
     * 63.0,63.4,63.7,64.1,64.6,65.4,66.2,67.0,67.7,68.5,69.3,70.1,
     * 70.6,70.9,70.9,70.4,69.3,68.0,66.9,66.0,65.2,64.5,63.7,63.0,
     * 62.3,61.7,61.2,60.9,60.6,60.4,60.2,60.1,59.8,59.6,59.4,59.4,
     * 58.5,58.8,59.2,59.4,59.4,59.2,58.9,58.9,59.2,59.7,60.5,61.2,
     * 61.7,62.0,62.3,62.6,63.3,64.1,65.1,65.9,66.4,66.9,67.5,68.4,
     * 69.0,69.2,68.9,68.1,66.7,65.3,64.2,63.4,62.7,62.0,61.1,60.4,
     * 59.7,59.2,58.9,58.6,58.4,58.4,58.4,58.5,58.5,58.4,58.3,58.3,
     * 57.6,57.8,57.9,57.9,57.7,57.5,57.4,57.5,57.9,58.6,59.3,59.8,
     * 60.3,60.5,60.8,61.2,62.0,63.0,64.1,64.9,65.4,65.7,66.2,66.8,
     * 67.2,66.8,66.0,64.8,63.4,62.1,61.3,60.7,60.1,59.3,58.2,57.1,
     * 56.6,56.3,56.1,56.0,56.0,56.2,56.6,56.9,57.1,57.3,57.4,57.5,
     * 54.3,54.9,55.4,55.6,55.6,55.3,55.1,55.0,55.3,55.8,56.5,57.2,
     * 57.8,58.3,58.7,59.4,60.3,61.2,62.0,62.8,63.4,63.9,64.1,64.1,
     * 63.8,63.3,62.5,61.5,60.2,58.9,57.7,56.7,56.0,55.5,55.1,54.8,
     * 54.6,54.3,53.9,53.4,53.2,53.1,53.3,53.4,53.5,53.5,53.6,53.8,
     * 52.9,53.6,54.2,54.6,54.6,54.3,54.0,53.7,53.8,54.1,54.8,55.7,
     * 56.4,56.9,57.4,58.0,58.8,59.6,60.2,60.8,61.7,62.4,62.6,62.2,
     * 61.5,60.7,59.8,58.9,57.9,56.8,55.5,54.5,53.8,53.4,53.3,53.3,
     * 53.5,53.2,52.7,52.1,51.7,51.7,51.7,51.9,51.9,52.0,52.1,52.4,
     * 51.8,52.5,53.2,53.6,53.6,53.3,52.9,52.5,52.4,52.6,53.3,54.2,
     * 55.0,55.6,56.1,56.7,57.3,57.9,58.1,58.6,59.7,60.7,60.9,60.3,
     * 59.2,58.1,57.1,56.3,55.5,54.6,53.5,52.5,51.8,51.4,51.4,51.7,
     * 52.0,51.9,51.4,50.8,50.4,50.3,50.4,50.5,50.6,50.7,50.8,51.2,
     * 50.9,51.7,52.4,52.9,52.9,52.5,52.0,51.5,51.3,51.4,52.1,53.0,
     * 53.9,54.5,55.0,55.5,56.0,56.4,56.4,56.8,58.0,59.3,59.5,58.7,
     * 57.4,56.1,54.9,54.1,53.5,52.8,51.8,50.8,50.1,49.7,49.8,50.4,
     * 50.9,50.9,50.3,49.7,49.3,49.2,49.3,49.4,49.5,49.6,49.8,50.2/

        
        if(xkp.gt.9.0) xkp=9.0
        kp1=int(xkp)+1
        xkp1=int(xkp)*1.0
        kp2=kp1+1
        if(kp2.gt.10) kp2=10

        do i=1,48 
           ab_mlat(i)=zp_mlat(i,kp1)+(xkp-xkp1)*
     &      	(zp_mlat(i,kp2)-zp_mlat(i,kp1))
           enddo
           
        cgmlat=-99.99

        if(xmlt.lt.0.0) return
        
        do i=1,48 
        	ab_mlt(i)=(i-1)*.5
			enddo
        i1=int(xmlt/0.5)+1
        if(i1.ge.48) i1=1
        i2=i1+1
      
        s1=(zp_mlat(i2,kp1)-zp_mlat(i1,kp1))/(ab_mlt(i2)-ab_mlt(i1))
        zmlkp1=zp_mlat(i1,kp1)+(xmlt-ab_mlt(i1))*s1
        s2=(zp_mlat(i2,kp2)-zp_mlat(i1,kp2))/(ab_mlt(i2)-ab_mlt(i1))
        zmlkp2=zp_mlat(i1,kp2)+(xmlt-ab_mlt(i1))*s2
        
        cgmlat=zmlkp1+(xkp-xkp1)*(zmlkp2-zmlkp1)
		return
		end
C
C
