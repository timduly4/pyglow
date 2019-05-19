c iriorbitmax.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
c Program for computing IRI F2 peak parameters along a satellite orbit
c-----------------------------------------------------------------------
c
c
      INTEGER           pad1(6),jdprof(77)
      DIMENSION         outf(20,1000),oar(100),jfi(6)
      LOGICAL		    jf(50)
      CHARACTER*2       timev(2)
      CHARACTER*3       uni(48),sopt
      CHARACTER*4       IMZ(8),MAP,xtex,coorv(2)
      CHARACTER*5       ITEXT(8)
      CHARACTER*6       dopt,pna(48)
      CHARACTER*8       bopt
      CHARACTER*9       topt,pname(6)
      CHARACTER*10      iopt
      CHARACTER*16      f1opt
      CHARACTER*50  	orbit_input,header,orbit_output

        data jfi/8,9,13,14,15,16/

	  COMMON/const2/icalls,nmono,iyearo,idaynro,rzino,igino,ut0

		icalls=0
		nmono=-1
		iyearo=-1
		idaynro=-1
		rzino=-1
		igino=-1
		ut0=-1
        
        do 6249 i=1,50
6249    oar(i)=-1.0

c open input and output files
c

	    IORB=15
	    type *,'name of file with orbit information'
	    read(5,*) orbit_input
        OPEN(IORB,FILE=orbit_input,STATUS='OLD',ERR=4321,
     &          FORM='FORMATTED')

	    IORB1=16
        type *,'name of output file'
	    read(5,*) orbit_output
        OPEN(IORB1,FILE=orbit_output,STATUS='NEW',ERR=4321,
     &          FORM='FORMATTED')

c read three header files from input file
c
        READ(IORB,1088) header
        write(IORB1,1088) header
        READ(IORB,1088) header
        write(IORB1,1088) header
        READ(IORB,1088) header
        write(IORB1,1088) header
1088    Format(A50)

        jm=0
        iut=1
        htec_max=0.0
        hxin=300.0

c read input file until EOF
c
1       READ(IORB,1087,ERR=4321,END=4322) inyy,indoy,inhh,inmm,inss,
     & 	   xlatin,xlonin
        type*, inyy,indoy,inhh,inmm,inss,xlatin,xlonin
1087    Format(I5,I6,I5,I6,I6,2F11.3)
        jmag=jm
        alati=xlatin
        along=xlonin
        iyyyy=inyy
        mmdd=-indoy
        dhour=inhh*1.0+(inmm+inss/60.)/60.+25.
        HEIBEG=hxin
        HEIEND=hxin
        HEISTP=1.      

        do i=1,30 
              jf(i)=.true.
              enddo
        jf(2)=.false.				  ! no temperatures
        jf(3)=.false.				  ! no ion composition
        jf(5)=.false.               ! URSI foF2 model
c        jf(6)=.false.               ! Newest ion composition model
        jf(12)=.false.              ! no konsol messages
c        jf(21)=.true.			      ! ion drift computed
c        jf(23)=.false.              ! TTS Te model is standard
c        jf(28)=.true.			      ! spread-F computed
        jf(29)=.false.              ! New Topside options
        jf(30)=.false.              ! NeQuick topside
       
        call IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &    HEIBEG,HEIEND,HEISTP,OUTF,OAR)

        WRITE(iorb1,7117) inyy,indoy,inhh,inmm,inss,xlatin,xlonin,
     &     oar(1),oar(2)
7117    Format(I5,I6,I5,I6,I6,2F11.3,2X,E12.6,F9.2)

        goto 1

4321    type*,'ERROR'
        type*,inyy,indoy,inhh,inmm,inss,xlatin,xlonin,hxin
        goto 9876
        
4322    type*,'END'
        type*,inyy,indoy,inhh,inmm,inss,xlatin,xlonin,hxin

9876	stop
	    end
