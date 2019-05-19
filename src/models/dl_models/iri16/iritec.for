c iritec.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------        
C
C contains IRIT13, IONCORR, IRI_TEC subroutines to computed the 
C total ionospheric electron content (TEC)
C
c-----------------------------------------------------------------------        
C Corrections
C
C  3/25/96 jmag in IRIT13 as input
C  8/31/97 hu=hr(i+1) i=6 out of bounds condition corrected
C  9/16/98 JF(17) added to input parameters; OUTF(11,50->100)
C  ?/ ?/99 Ne(h) restricted to values smaller than NmF2 for topside        
C 11/15/99 JF(20) instead of JF(17)
C 10/16/00 if hr(i) gt hend then hr(i)=hend
C 12/14/00 jf(30),outf(20,100),oarr(50)
C
C Version-mm/dd/yy-Description (person reporting correction)
C 2000.01 05/07/01 current version
c 2000.02 10/28/02 replace TAB/6 blanks, enforce 72/line (D. Simpson)
c 2000.03 11/08/02 common block1 in iri_tec with F1reg
c 2007.00 05/18/07 Release of IRI-2007
c 2007.02 10/31/08 outf(.,100) -> outf(.,500)
c
C 2012.00 10/05/11 IRI-2012: bottomside B0 B1 model (SHAMDB0D, SHAB1D),
C 2012.00 10/05/11    bottomside Ni model (iriflip.for), auroral foE
C 2012.00 10/05/11    storm model (storme_ap), Te with PF10.7 (elteik),
C 2012.00 10/05/11    oval kp model (auroral_boundary), IGRF-11(igrf.for), 
C 2012.00 10/05/11    NRLMSIS00 (cira.for), CGM coordinates, F10.7 daily
C 2012.00 10/05/11    81-day 365-day indices (apf107.dat), ap->kp (ckp),
C 2012.00 10/05/11    array size change jf(50) outf(20,1000), oarr(100).

c-----------------------------------------------------------------------        
C
C
        subroutine IRIT13(ALATI,ALONG,jmag,jf,iy,md,hour,hbeg,hend,
     &                          tec,tecb,tect)
c-----------------------------------------------------------------------        
c Program for numerical integration of IRI-94 profiles from h=100km
C to h=alth. 
C       
C  INPUT:  ALATI,ALONG  LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C          jmag         =0 geographic   =1 geomagnetic coordinates
C          jf(1:50)     =.true./.false. flags; explained in IRISUB.FOR
C          iy,md        date as yyyy and mmdd (or -ddd)
C          hour         decimal hours LT (or UT+25)
c          hbeg,hend    upper and lower integration limits in km
C 
C  OUTPUT: TEC          Total Electron Content in m-2
C          tecb,tect    percentage of bottomside and topside content
c-----------------------------------------------------------------------        

        dimension       outf(20,1000),oarr(100)
        logical         jf(50)

c
c  select various options and choices for IRI-94
c

        tec = -111.
        tect= -111.
        tecb= -111.

c
c initialize IRI parameter in COMMON blocks
c

        abeg=hbeg
        aend=hend
        astp=hend-hbeg
        call IRI_SUB(JF,JMAG,ALATI,ALONG,IY,MD,HOUR,
     &          abeg,aend,astp,OUTF,OARR)

c
c  calculate total electron content (TEC) in m-2 using the
c  stepsize selection 2 (highest accuracy)
c

        call iri_tec (hbeg,hend,2,tec,tect,tecb)

        return
        end
c
c
        real function ioncorr(tec,f)
c-----------------------------------------------------------------------        
c computes ionospheric correction IONCORR (in m) for given vertical
c ionospheric electron content TEC (in m-2) and frequency f (in Hz)
c-----------------------------------------------------------------------        
        ioncorr = 40.3 * tec / (f*f)
        return
        end
c
c
        subroutine iri_tec (hstart,hend,istep,tectot,tectop,tecbot)
c-----------------------------------------------------------------------        
C subroutine to compute the total ionospheric content
C INPUT:      
C   hstart  altitude (in km) where integration should start
C   hend    altitude (in km) where integration should end
C   istep   =0 [fast, but higher uncertainty <5%]
C           =1 [standard, recommended]
C           =2 [stepsize of 1 km; best TEC, longest CPU time]
C OUTPUT:
C   tectot  total ionospheric content in tec-units (10^16 m^-2)
C   tectop  topside content (in %)
C   tecbot  bottomside content (in %)
C
C The different stepsizes for the numerical integration are 
c defined as follows (h1=100km, h2=hmF2-10km, h3=hmF2+10km, 
c h4=hmF2+150km, h5=hmF2+250km):
C       istep   h1-h2   h2-h3   h3-h4   h4-h5   h5-hend
C       0       2.0km   1.0km   2.5km   exponential approximation
C       1       2.0km   1.0km   2.5km   10.0km  30.0km
C       2       1.0km   0.5km   1.0km   1.0km   1.0km   
C
c-----------------------------------------------------------------------        

        logical         expo
        dimension       step(5),hr(6)
        logical     	f1reg
        common  /block1/hmf2,xnmf2,hmf1,f1reg
     &         /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU
C NEW-GUL------------------------------
c     &         /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,hht,TAU

ctest   
        save

        expo = .false.
        numstep = 5
        xnorm = xnmf2/1000.
C NEW-2003: Half-density: XNETOP at htop in [hmf2,1000 km)
        xxx=xnmf2/2.
        ht1=hmf2
        xne1=xnmf2
        ht2=ht1
        xne2=xne1
        hht=0.0
C NEW-2003: Half-density: XNETOP at htop in [hmf2,1000 km)

        hr(1) = 100.
        hr(2) = hmf2-10.
        hr(3) = hmf2+10.
        hr(4) = hmf2+150.
        hr(5) = hmf2+250.
        hr(6) = hend
        do 2918 i=2,6 
2918            if (hr(i).gt.hend) hr(i)=hend

        if (istep.eq.0) then 
                step(1)=2.0
                step(2)=1.0
                step(3)=2.5
                step(4)=5.
                if (hend.gt.hr(5)) expo=.true.
                endif

        if (istep.eq.1) then
                step(1)=2.0
                step(2)=1.0
                step(3)=2.5
                step(4)=10.0
                step(5)=30.0
                endif

        if (istep.eq.2) then
                step(1)=1.0
                step(2)=0.5
                step(3)=1.0
                step(4)=1.0
                step(5)=1.0
                endif

        sumtop = 0.0
        sumbot = 0.0
C
C find the starting point for the integration
C

        i=0
        ia=1
3       i=i+1
        h=hr(i)
        if(hstart.gt.h) then
                hr(i)=hstart
                ia=i
                goto 3
                endif
C
C start the numerical integration
C
        i=ia
        h=hr(i)
        hu=hr(i+1)
        delx = step(i)
1       h = h + delx
        hh = h
        if (h.ge.hu) then
                delx = hu - h + delx
                hx = hu - delx/2.
                YNE = XE_1(hx)
                if((hx.gt.hmf2).and.(yne.gt.xnmf2)) yne=xnmf2
                yyy = yne * delx / xnorm
                i=i+1
            if(i.lt.6) then
                      h = hr(i)
                      hu = hr(i+1)
                            delx = step(i)
                  endif
        else
              hx = h - delx/2.
                YNE = XE_1(hx)
                if((hx.gt.hmf2).and.(yne.gt.xnmf2)) yne=xnmf2
                yyy = yne * delx / xnorm
        endif
        if (hx.le.hmf2) then
                sumbot = sumbot + yyy
        else
                sumtop = sumtop + yyy

C NEW-GUL: remember xne2 at ht2 :
                ht2=hx
                xne2=yne
C NEW-GUL------------------------------

        endif

C NEW-GUL: interpolate for htop
        if ((hx.le.hmf2).or.(hht.gt.0.0)) goto 4
        if ((xxx.le.xne1).and.(xxx.gt.xne2)) then
           hht=ht1+(ht2-ht1)/(xne2-xne1)*(xxx-xne1)
        else
           ht1=ht2
           xne1=xne2
        endif
C NEW-GUL------------------------------

4       if (expo.and.(hh.ge.hr(4))) goto 5
        if (hh.lt.hend.and.i.lt.6) goto 1

        zzz = sumtop + sumbot
        tectop = sumtop / zzz * 100.
        tecbot = sumbot / zzz * 100.
        tectot = zzz * xnmf2    
        return

5       num_step = 3
        hei_top = hr(4)
        hei_end = hend
        top_end = hei_end - hei_top
        del_hei = top_end / num_step
        xntop = xe_1(hei_end)/xnmf2

        if(xntop.gt.0.9999) then
                ss_t = top_end  
                goto 2345
                endif

        hei_2 = hei_top
        hei_3 = hei_2 + del_hei
        hei_4 = hei_3 + del_hei
        hei_5 = hei_end

        hss = top_end / 4.
C       hss = 360.
        xkk = exp ( - top_end / hss ) - 1.
        x_2 = hei_2
        x_3 =hei_top-hss*alog(xkk*(hei_3 - hei_top)/top_end + 1.) 
        x_4 =hei_top-hss*alog(xkk*(hei_4 - hei_top)/top_end + 1.)
        x_5 = hei_end

        ed_2 = xe_1(x_2)/xnmf2
          if(ed_2.gt.1.) ed_2=1.
        ed_3 = xe_1(x_3)/xnmf2
          if(ed_3.gt.1.) ed_3=1.
        ed_4 = xe_1(x_4)/xnmf2
          if(ed_4.gt.1.) ed_4=1.
        ed_5 = xntop
        if(ed_3.eq.ed_2) then
         ss_2 = ed_3 * (x_3 - x_2)
        else
         ss_2=( ed_3 - ed_2 ) * ( x_3 - x_2 ) / alog ( ed_3 / ed_2 )
        endif
        if(ed_4.eq.ed_3) then
         ss_3 = ed_4 * (x_4 - x_3)
        else
         ss_3=( ed_4 - ed_3 ) * ( x_4 - x_3 ) / alog ( ed_4 / ed_3 )
        endif
        if(ed_5.eq.ed_4) then
         ss_4 = ed_5 * (x_5 - x_4)
        else
         ss_4=( ed_5 - ed_4 ) * ( x_5 - x_4 ) / alog ( ed_5 / ed_4 )
        endif

        ss_t = ss_2 + ss_3 + ss_4 

2345    sumtop = sumtop + ss_t * 1000.
        
        zzz = sumtop + sumbot
        tectop = sumtop / zzz * 100.
        tecbot = sumbot / zzz * 100.
        tectot = zzz * xnmf2

      RETURN
      END

