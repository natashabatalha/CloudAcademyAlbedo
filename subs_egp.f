C-----------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
C 
C THIS FILE CONTAINS SETUP + SUBROUTINES FOR PROGRAM GEOMETRIC 
C IN EGPALB.F
C
C Requires FORTRAN files: egpalb.f, subs_egp.f, raman2.f, trist.f
C 
C Requires input files: 
C *.in, 
C *.cld, 
C GUILLOT.DAT,
C wave_EGP.DAT, 
C H2CIA.DAT,
C full_abs_log_noTiO.736.long (freedman .long to 30u, otherwise 5u), !Old Version
C tpint75.out (for methane in the visible) !Old Version
C
C If Hazes on, files: HAZEINFO, RAGES.HAZ
C 
C Dumps output to both file and screen
C
C......................................................................
C-----------------------------------------------------------------------            
C ATMSETUP (input files here)

      SUBROUTINE ATMSETUP(NLEVEL,Z,TEMP,PRESS,DEN,XMU,CH4,H2,HE,XN2,
     $   CO,H2O,CO2,XNH3,H,Hmin,elec,xNa,xK, RB,CS,FEH,CRH,TIO, VO,
     $   x736, x1060,IPRINT)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (NMAX=201)
      integer pn,tn,pe , cloud_int
      logical do_clouds
      parameter (pn=51,pe=20)
      character*25 header
      CHARACTER*10 :: AA
C INPUTS:  
C NLEVEL    NUMBER OF ALTITUDE LEVELS, J=1 IS AT THE TOP
C Z         ALTITUDE GRID IN KM
C FH2       MIXING RATIO BY NUMBER OF H2

C OUTPUTS AT EACH LEVEL (NOT LAYER AVERAGES):
C TEMP (K), PRESS(BARS), DEN(CM-3), XMU = MEAN MOLECULAR WEIGHT,
C AND CH4, H2, XN2, AR ARE THE NUMBER MIXING RATIOS OF THE GASES

C INTERNAL VARIABLES
C      DIMENSION Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL) !Old Version
C      DIMENSION CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),CO(NLEVEL)    !Old Version
C      DIMENSION HE(NLEVEL),XMU(NLEVEL),elec(NLEVEL)              !Old Version
C     DIMENSION H2O(NLEVEL),XNH3(NLEVEL),H(NLEVEL),Hmin(NLEVEL)  !Old Version
      DIMENSION Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
      DIMENSION CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),CO(NLEVEL),CO2(NLEVEL)
      DIMENSION HE(NLEVEL),XMU(NLEVEL),elec(NLEVEL)
      DIMENSION H2O(NLEVEL),XNH3(NLEVEL),H(NLEVEL),Hmin(NLEVEL)
      DIMENSION XNA(NLEVEL),XK(NLEVEL),H2S(NLEVEL)
      DIMENSION RB(NLEVEL),CS(NLEVEL),FEH(NLEVEL),CRH(NLEVEL)
      DIMENSION TIO(NLEVEL),VO(NLEVEL)   

      common /tristan/ p(pn),t(pn),ro(pn,pn),xnel(pn,pn),ab(pe,pn,pn),code(pe),np,nt

      common /grav/ gin
      common /cloudy/ do_clouds
   
C MEAN MOLECULAR WEIGHT DATA	  
      DATA XMN2/28.0134d0/
      DATA XMH2/2.0158d0/
      DATA XMHE/4.0026d0/
      DATA XMCH4/16.0426d0/
      DATA XCO/28.0104d0/
      DATA XH2O/18.0152d0/
      DATA XMNH3/17.0304d0/
      DATA XMH/1.0d0/
      DATA XCO2/44.01d0/
      DATA XMNA/22.989769d0/
      DATA XMK/39.0983d0/
      DATA XH2S/34.0809d0/
      DATA XRB/85.4678d0/
      DATA XCS/132.90545d0/
      DATA XFEH/56.85d0/
      DATA XCRH/53.00d0/
      DATA XTIO/64.06d0/
      DATA XVO/66.94d0/

c===========READ IN DO CLOUDS==============================================

      OPEN (unit=999,file= 'int.params',status='old')
      READ(999,*) AA, cloud_int 
      CLOSE(999)

C     just for clarity switch to bool 
      if(cloud_int.eq.1) then
        do_clouds=.TRUE.
      else 
        do_clouds=.FALSE.
      endif


C CLOUDS ON/OFF      
c         do_clouds = .true.
c         do_clouds = .true.
      
        if (do_clouds) then
        Print *,' This run does include the Ackerman clouds' 

C SET PATH TO CLOUD FILE
C......................................................................
        open(unit=22, file='input.cld',status='old')        
       read (22,*) header
C......................................................................		
        else
        Print *,' This run does NOT include Ackerman clouds'
	endif
        
C SET PATH TO PT FILE
C......................................................................
      OPEN (55,FILE='input.pt',STATUS='OLD')

C......................................................................
	  
C Gas composition info directly from input P-T profile 
C These are LEVEL quantities
C OPTCV will average them to find LAYER quantities     

C  Read in header information
        read (55,*) header

C Read in the gravity directly from the input file	 
	     READ(55,*) gin  !normal

C Start loop, display J levels
       DO J=1,NLEVEL
c       write(*,*)j


c Very Old File Format
c n, P, T,  CH4, H2O, NH3, CO, N2, H, H2, He, H-, e-
c        READ(55,*) iDUM,press(j),TEMP(J),CH4(j),H2O(j),xNH3(j), CO(j),xN2(j),
c     &    H(j),H2(j),He(j),Hmin(j),elec(j)

C  Old Standard File Format
C  n P T e- H2 H H+ H- VO TiO CO2 He H2O CH4 CO NH3 N2 PH3 H2S  Fe Na  K
c        READ(55,*) iDUM,press(j),TEMP(J),elec(j),H2(j),H(j),xDUM2,Hmin(j),xdum3,xdum4,
c     &    xdum5,He(j),H2O(j),CH4(j),CO(j),xNH3(j),xN2(j),xdum6,xdum7,
c     &         xdum8, xdum9, xdum10

C  Standard File Format
C     n P T e- H2 H H+ H- VO TiO CO2 He H2O CH4 CO NH3 N2 PH3 H2S  Fe Na  K
C     Note: H+, PH3, and Fe not currently used.     
        READ(55,*) iDUM,press(j),TEMP(J),elec(j),H2(j),H(j),xDUM2,
     &    Hmin(j),VO(j),TIO(j),CO2(j),He(j),H2O(j),CH4(j),
     &    CO(j),xNH3(j),xN2(j),XDUM3,H2S(j),xDUM4,xNa(j),xK(j)      
     
C     High Metallicity (Caroline) File Format
C        READ(55,*) iDUM,press(j),TEMP(J),elec(j),H2(j),H(j),xDUM2,
C     &    Hmin(j),VO(j),TIO(j),CO2(j),He(j),H2O(j),CH4(j), xdum8,
C     &    CO(j),xNH3(j),xN2(j),xdum6,H2S(j),
C     &    xNa(j), xK(j),RB(j),CS(j),FEH(j),CRH(j), x736(j), x1060(j)    
        
        
c 1 n 1
c 2 P 0.10000E-05
c 3 T 1029.52
c 4 e- 0.382E-10
c 5 H2 0.838E+00
c 6 H 0.448E-05
c 7 H+ 0.359E-54
c 8 H- 0.185E-19
c 9 VO 0.637E-15
c 10 TiO 0.993E-17
c 11 CO2 0.196E-06
c 12 He 0.162E+00
c 13 H2O 0.265E-03
c 14 CH4 0.189E-13
c 15 CO 0.488E-03
c 16 NH3 0.303E-11
c 17 N2 0.672E-04
c 18 PH3 0.360E-11
c 19 H2S 0.298E-04
c 20 Fe 0.708E-07
c 21 Na 0.390E-05
c 22 K 0.225E-06

C Adhoc change to CH4 abundance
C        CH4(J)=CH4(J)*1.d0

C Adhoc removal of atmospheric species
C        elec(j)=elec(j)*0.d0
C        H2(j)=H2(j)*0.d0
C        H(j)=H(j)*0.d0
C        Hmin(j)=Hmin(j)*0.d0
C        VO(j)=VO(j)*0.d0
C        TIO(j)=TIO(j)*0.d0
C        CO2(j)=CO2(j)*0.d0
C        He(j)=He(j)*0.d0
C        H2O(j)=H2O(j)*0.d0
C        CH4(j)=CH4(j)*0.d0
C        CO(j)=CO(j)*0.d0
C        xNH3(j)=xNH3(j)*0.d0
C        xN2(j)=xN2(j)*0.d0
C        H2S(j)=H2S(j)*0.d0
C        xNa(j)=xNa(j)*0.d0
C        xK(j)=xK(j)*0.d0 
        
        
        
C  Old XMU Calculation	     
C        XMU(J)=XMN2*XN2(J)+XMCH4*CH4(J)+XMHE*HE(J)+XMH2*H2(J)+
C     &         XCO*CO(j)+xh2o*h2o(j)+xmnh3*xnh3(j)+xmh*H(J)
C
c     write(*,*)'XMU=',XMU(J)

C     New XMU Calculation
        XMU(J)=XMN2*XN2(J)+XMCH4*CH4(J)+XMHE*HE(J)+XMH2*H2(J)+
     &         XH2S*H2S(j)+XCO*CO(j)+xh2o*h2o(j)+xmnh3*xnh3(j)+xmh*H(J)+
     &         XCO2*CO2(j)+XMNa*xNa(j)+XMK*xK(j)+XTIO*TIO(j)+XVO*VO(j)
        

C     High Metallicity XMU Calculation
C                XMU(J)=XMN2*XN2(J)+XMCH4*CH4(J)+XMHE*HE(J)+XMH2*H2(J)+XH2S*H2S(j)+
C     &         XCO*CO(j)+xh2o*h2o(j)+xmnh3*xnh3(j)+xmh*H(J)+XCO2*CO2(j)+XMNa*xNa(j)+XMK*xK(j)+
C     &         XRB*RB(j)+XCS*CS(j)+XFEH*FEH(j)+XCRH*CRH(j)+XTIO*TIO(j)+XVO*VO(j)
        
      ENDDO
      close(55)
C......................................................................
C Temporary edit to multiply all pressure levels by x10
C
C      DO J=1,NLEVEL
C         PRESS(J) = 10*PRESS(J)
C      ENDDO
C
C......................................................................      
      

C WRITE TO SCREEN
c	WRITE (6,139)

C CALCULATE DENSITY (k Boltzman 1.38E-16 erg/K)
      DO J=NLEVEL,1,-1
         DEN(J)=PRESS(J)*1E6/(1.381E-16*TEMP(J))
	     Z(J)=0.d0
C AT THIS POINT WE HAVE THE TEMP, PRESS, DEN AND XMU VALUES.
c         WRITE(6,140)J,PRESS(J),DEN(J),TEMP(J),H2(J)*100.
c     &    ,HE(J)*100.,CH4(J)*100.,H2O(J)*100.,XMU(J)

C This is the old methodology that used the premixed tables             
C WE ONLY NEED TO KEEP CH4 from the input file, since Freedman's opacity
C table includes the others. So zero them back out at this point until
C they are read in from Freedman's file.
c     
c	 h2o(j)=0.d0  !zero out H2O, NH3, and CO
c	 xnh3(j)=0.d0 !Freedman's files have opacity for these 3 species
c	 co(j)=0.d0   !at every wavelength (only need methane)
c	 XN2(J)=0.D0

     
      ENDDO 

  139 FORMAT(//'   FIRST GUESS ATMOSPHERE AT LEVELS'/
     &' LVL P(BARS)  # DEN(CM-3) TEMP '
     & , '   %H2    %HE  %CH4       %H2O    MU'  )
  140 FORMAT(1X,I3,1P2E10.3,0PF8.2,2X,2F6.2,3F7.3)

C END INTITAL BACKGROUND ATMOSPHERE SETUP FOR URANUS

C  READ IN HYDROGEN ABUNDANCE DATA, LATER USED IN OPTCV (TRISTAN)
C......................................................................		
      open(unit=35,file='../inputs/GUILLOT.DAT',status='old',
     &     form='formatted')
C......................................................................		
	 
 35   format(1p,8d12.5)
      read(35,*)np,nt,nmol
      read(35,35)(p(ip),ip=1,np)
      read(35,35)(t(it),it=1,nt)
      read(35,35)((ro(it,ip),it=1,nt),ip=1,np)
      read(35,'(6f13.2)')(code(im),im=1,nmol)
      read(35,35)((xnel(it,ip),it=1,nt),ip=1,np)
      do im=1,nmol
       read(35,35)((ab(im,it,ip),it=1,nt),ip=1,np)
      enddo
      close(35)

      RETURN
      END
C----------------------------------------------------------------------- 
C End ATMSETUP
C----------------------------------------------------------------------- 
C-----------------------------------------------------------------------               
C OPTCV (called in GEOM, tweak phase function parameters, scattering here)
C
C THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE VISIBLE
C IT CALCUALTES FOR EACH LAYER J, FOR EACH SPECRAL INTERVAL K IN THE VIS
C
C LAYER: WBAR, DTAU, COSBAR
C LEVEL: TAU

      SUBROUTINE OPTCV(IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      include 'freedman/int.params'
      include 'freedman/limits.com'
      include 'freedman/abs_data.com'
      include 'prog_params'     
      include 'pangle_params'

      PARAMETER (NTEMPS3=198)
      PARAMETER (kline=196) !number of wavelengths in cloud file
      COMMON /ATM/ Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
c      COMMON /GASS/ CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),CO(NLEVEL)
c     & ,HE(NLEVEL),XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER)
c     & ,H2O(NLEVEL),XNH3(NLEVEL),H(NLEVEL),Hmin(NLEVEL),elec(NLEVEL)
      COMMON /GASS/ CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),CO(NLEVEL),
     &  CO2(NLEVEL),HE(NLEVEL),XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER),
     &  H2O(NLEVEL),XNH3(NLEVEL),H(NLEVEL),Hmin(NLEVEL),elec(NLEVEL), 
     &  xNA(NLEVEL), xK(NLEVEL),H2S(NLEVEL),TIO(NLEVEL), VO(NLEVEL),
     &  RB(NLEVEL), CS(NLEVEL), FEH(NLEVEL), CRH(NLEVEL) 
      COMMON /VISGAS/ SOLARF(NSPECV)
      COMMON /AERSOL/ RADIUS(NLAYER), XNUMB(NLAYER), SIGMA(NLAYER)
     & , REALI(NSPECI), XIMGI(NSPECI), REALV(NSPECV), XIMGV(NSPECV)
      COMMON /CLEAR/ CLA(NSPECV),ICLEAR
      COMMON /AERSOLV/ TAERH(NLAYER,NSPECV),QSCTH(NLAYER,NSPECV),
     &         CBARH(NLAYER,NSPECV),XSCT(NLAYER)
      COMMON /CLOUDV/ QEXTC,QSCTC(NSPECV),CBARC(NSPECV)
      COMMON /CLOUD2V/ TAU2(NLAYER),WBAR2(NLAYER),CBAR2(NLAYER)
      COMMON /CLOUD/ RADCLD(NLAYER), XNCLD(NLAYER), XSCTC
     & , RCLDI(NSPECI), XICLDI(NSPECI), RCLDV(NSPECV), XICLDV(NSPECV)
      COMMON /TAUS/ 
     &  TAURV(NSPECV),TAUHV(NSPECV),TAUCV(NSPECV),TAUGV(NSPECV)
     &  ,TAUC2V(NSPECV),TAUCLV(NSPECV)
      COMMON /OPTICV/ DTAUV(NLAYER,NSPECV),TAUV(NLEVEL,NSPECV)
     &                 ,WBARV(NLAYER,NSPECV), COSBV(NLAYER,NSPECV)
     &                 ,NOPTHCK(NSPECV), LAYERBOT(NSPECV)
      COMMON /RAY2/ GCOS2(NLAYER,NSPECV)
      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      COMMON /PLANT/ CSUBP,RSFI,RSFV,F0PI
      COMMON /TEMPS3/ TSHH(NTEMPS3),TSHE(NTEMPS3)
      COMMON /ADJUST/ RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,RCLOUD
      COMMON /CONST/RGAS,RHOP,PI,SIGSB,RTOD 
      COMMON /TESTS/FHZE,FRAY,FCLDC
      COMMON /NEWC/ taucx(nlayer,nspecv),gcx(nlayer,nspecv)
     & 			,wbarcx(nlayer,nspecv)
C      COMMON /XKAPPA/ XKAPCH4(NSPECV,NLAYER),XKAPCO(NSPECV,NLAYER)
C     &     ,XKAPH2O(NSPECV,NLAYER),XKAPNH3(NSPECV,NLAYER)
      COMMON /XKAPPA/ XKAPCH4(NSPECV,NLAYER),XKAPCO(NSPECV,NLAYER)
     &     ,XKAPH2O(NSPECV,NLAYER),XKAPNH3(NSPECV,NLAYER),XKAPCO2(NSPECV,NLAYER),
     & 		XKAPNA(NSPECV,NLAYER), XKAPK(NSPECV,NLAYER), XKAPH2S(NSPECV,NLAYER),
     &	    XKAPRB(NSPECV,NLAYER), XKAPCS(NSPECV,NLAYER), XKAPFEH(NSPECV,NLAYER), 
     & 		XKAPCRH(NSPECV,NLAYER), XKAPTIO(NSPECV,NLAYER),XKAPVO(NSPECV,NLAYER)
      COMMON /XKAP/XKCH4(790)
      common /cloudy/ do_clouds
      COMMON /OPACINDEX/ X736(NLAYER),X1060(NLAYER)

      DIMENSION gasmtx(3),DPOL(3),GNU(2,3),TAUR(3),xxkappa(nspecv)
      dimension xokappa(nspecv), tauhc(nspecv), xfreed(nlayer,nspecv)
      dimension taued(nlayer,kline),ged(nlayer,kline),wbared(nlayer,kline)
      dimension down1(kline),down2(kline),down3(kline),waveEGP(kline),wnumEGP(kline)
      dimension back1(nspecv), back2(nspecv), back3(nspecv)
C      DIMENSION aCH4(NLAYER),aXN2(NLAYER),aH2(NLAYER),aCO(NLAYER)
C      DIMENSION aHE(NLAYER),aXMU(NLAYER),aelec(NLAYER)
C     DIMENSION aH2O(NLAYER),aXNH3(NLAYER),aH(NLAYER),aHmin(NLAYER)
      DIMENSION aCH4(NLAYER),aXN2(NLAYER),aH2(NLAYER),aCO(NLAYER),aCO2(NLAYER)
      DIMENSION aHE(NLAYER),aXMU(NLAYER),aelec(NLAYER)
      DIMENSION aH2O(NLAYER),aXNH3(NLAYER),aH(NLAYER),aHmin(NLAYER)
      DIMENSION axNa(NLAYER),axK(NLAYER), aH2S(NLAYER),aRb(NLAYER)
      DIMENSION aCs(NLAYER), aFeH(NLAYER), aCrH(NLAYER), aTiO(NLAYER), aVO(NLAYER)
      logical do_clouds
      character*25 header
      character*20 charJ736
      character*20 charJ1060
	  

C Planetary phase angle in radians (0 = face on)      
      COMMON /OBSPHA/ PLPHA, PANGLE(NPHAS)
			   
      CHARACTER(80) TFNAME
	  CHARACTER(80) TFNAME2
      CHARACTER(80) TFNAME3
	  
	  
      AMU=1.66057E-24
      xk_boltz=1.381d-16
      dpol(1)=1.022 ! Rayleigh scattering constants
      dpol(2)=1.0
      dpol(3)=1.0
      gnu(1,1)=1.355e-4
      gnu(2,1)=1.235e-6
      gnu(1,2)=3.469e-5
      gnu(2,2)=8.139e-8
      gnu(1,3)=4.318e-4
      gnu(2,3)=3.408e-6
      XN0 = 2.687E19
 
	
c To know what K corresponds to what W#, we must read in a file of the
c W# grid used by the EGP code.
C....................................................................... 
	open(unit=23,file='wave_EGP.dat',status='old') !for BD's wave_180.dat
C....................................................................... 	
	  do K=1,kline
	  read (23,*) dum1, waveEGP(K), wnumEGP(K)
	  enddo
	close(23)	
	
c A bit crude, but now we know at each layer J the tau, g, and
c wbar as a function of wavenumber.  But now we need this on our finer
c grid from 0.35 to 5 microns, at the NSPECV wavelengths.


c Here we'll read in the data from the same *.cld cloud data file, 
c originally input in ATMSETUP, that was generated by the EGP atmosphere code.  
c At NLAYER layers, and 196 wavelength bins, we read in the
c optical depth, assymetry factor g, and single scattering albedo wbar.
c Note: kline is 196 (for EGPs) or 180 (for BDs)

        if (do_clouds) then
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	write(*,*)'cloud interval kline=',kline
	write(*,*)'NLAYER=',NLAYER
	write(*,*)'NLEVEL=',NLEVEL
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	read (22,*) header

	do J=1,NLAYER
	  do K=1,kline                  !optical detph, asymmetry, single scattering
c          read (22,*) dum1, dum2, taued(J,K), ged(J,K), wbared(J,K)
           read (22,*) dum1, dum2, taued(J,K), ged(J,K), wbared(J,K), dum100
 	  enddo
	enddo
	close(22)
	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	Print *,' This run does include the Ackerman clouds ' 
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	do J=1,NLAYER
C Add the cloud	
	  do K=1,kline
	  down1(K)=taued(J,K)
	  down2(K)=ged(J,K)
	  down3(K)=wbared(J,K)
	  enddo
C Interpolate to high resolution grid	  
	  call interp(wnumEGP,down1,kline,back1)
	  call interp(wnumEGP,down2,kline,back2)
	  call interp(wnumEGP,down3,kline,back3)
	
C Over the full nspecv intervals	
	  do K=1,nspecv

c       skip upper layers for low fsed garbage

      if (J.ge.(3)) then
      taucx(J,K)=back1(K) !optical depth
      gcx(J,K)=back2(K) !asymmetry factor
      wbarcx(J,K)=back3(K) !single scattering albedo
      else 
      taucx(J,K)=0.0d0  !zero out the clouds if we don't need them!
      gcx(J,K)=0.0d0
      wbarcx(J,K)=0.0d0       
      endif
c      gcx(J,K)=0.0d0
c	  wbarcx(J,K)=1.0
      enddo
	
	enddo
	close(22)
	
	else
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Print *,' This run does NOT include Ackerman clouds'
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	do J=1,NLAYER
C Over the full nspecv intervals	
	  do K=1,nspecv
	  taucx(J,K)=0.0d0  !zero out the clouds if we don't need them!
	  gcx(J,K)=0.0d0
	  wbarcx(J,K)=0.0d0
	  enddo
	enddo

	endif


C  This is the old opacity method that uses pre-mixed tables
C  READ IN OPACITY DATA		
C Here we'll read in the full Freedman pretabulated x-section data.  We will
C interpolate in it as a function of P and T.  
C......................................................................
C Generic name created to allow for batch processing
C       open(unit=20,file='opacity_data',status='old',form='unformatted')
c No TiO out to 5 um
c      open(unit=20,file='../freedman/full_abs_log_noTiO.736_2',status='old',form='unformatted')
c No TiO out to 30 um      
C      open(unit=20,file='../freedman/full_abs_log_noTiO.736.long',status='old',form='unformatted')
C new kcahoy
C 1x out to 5 micron
C      open(unit=20,file='../freedman/full_abs_log.736',status='old',form='unformatted')
C 3x out to 5 micron
C      open(unit=20,file='../freedman/full_abs_log.+0.5_CO2',status='old',form='unformatted')
c Previous working version below
C 3x out to 5 micron, NO H2O
c      open(unit=20,file='../freedman/full_abs_log_no_h2o',status='old',form='unformatted')
C 3x out to 5 micron, NO H2O version 2
C      open(unit=20,file='../freedman/+0.5_CO2.no.H2O.full_abs_log',status='old',form='unformatted')
C 3x out to 5 micron, NO H2O, NH3, or H2S
C      open(unit=20,file='../freedman/+0.5_CO2.no.H2O.nh3.h2s.full_abs_log',status='old',form='unformatted')
C 3x out to 5 micron, NO H2O, NH3
C      open(unit=20,file='../freedman/+0.5_CO2.no.H2O.nh3.full_abs_log',status='old',form='unformatted') 
C 10x out to 5 micron
c      open(unit=20,file='../freedman/full_abs_log.+1.0_CO2',status='old',form='unformatted')
C 30x out to 5 micron
c      open(unit=20,file='../freedman/full_abs_log.+1.5_CO2',status='old',form='unformatted')
C......................................................................	  

c commenting out all of these lines that deal with the freedman opacities
        
C Read log of absorption coeffs
C      read (20) full_abs ! equiv to xabs
c      write(*,*) 'Full_ABS=',xabs
C      close (20)

C Convert coarse P to log
C      do j=1,NP
C        log_CP(j) = log(Coarse_P(j)) 
C      end do

C     read in file that contains the indicies for the opacity files
C     NKL generates these files from the *.pt files using an IDL code
        open(unit=20, file='input.ind', status='old')
        
        
        
      DO J=1,NLAYER ! mean conditions in each layer

         read(20, *) PBAR, TBAR, x736(J), pdum1, tdum1, x1060(J), pdum2, tdum2
C     TBAR and PBAR are now calculated when the opacity file indicies are calculated
C        TBAR=0.5d0*(TEMP(J)+TEMP(J+1)) 
C        PBAR=DSQRT(PRESS(J)*PRESS(J+1))
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
c        print *,'input',tbar,pbar
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

C The following gets us xsection(T,P) at 2000 frequencies for 1 layer
c     call to do interpolations over all frequencies
c           get_K(t,p) ! call this routine for each t:p combination.
c                      ! T in Kelvins, p in millibars
c                      ! returns an array 'xreturn' [in abs_data.com] of NF[reqs]
c                      ! interpolated points as numbers
C
C	    pmil=pbar*1.d3  !P must be in mbar for get_K
C        call get_K(tbar,pmil) ! interpolate from opacity table to this tbar, pmil

C Old formatting/method        
C We need these abundances for the LAYERS, not LEVELS. 	
C These were read in from the planet/BD's *.in file
C	aCH4(J)=0.5d0*(CH4(J)+CH4(J+1))
C	aXN2(J)=0.5d0*(XN2(J)+XN2(J+1))
C	aH2(J)=0.5d0*(H2(J)+H2(J+1))
C	aCO(J)=0.5d0*(CO(J)+CO(J+1))
C	aHE(J)=0.5d0*(HE(J)+HE(J+1))
C	aXMU(J)=0.5d0*(XMU(J)+XMU(J+1))
C	aelec(J)=0.5d0*(elec(J)+elec(J+1))
C	aH2O(J)=0.5d0*(H2O(J)+H2O(J+1))
C	aXNH3(J)=0.5d0*(XNH3(J)+XNH3(J+1))
C	aH(J)=0.5d0*(H(J)+H(J+1))
C	aHmin(J)=0.5d0*(Hmin(J)+Hmin(J+1))

	aXMU(J)=0.5d0*(XMU(J)+XMU(J+1))
	aCH4(J)= 0.5d0*(CH4(J)+CH4(J+1))
	aXN2(J)= 0.5d0*(XN2(J)+XN2(J+1))
	aH2(J)=  0.5d0*(H2(J)+H2(J+1))
	aCO(J)=  0.5d0*(CO(J)+CO(J+1))
	aCO2(J)= 0.5d0*(CO2(J)+CO2(J+1))
	aHE(J)=  0.5d0*(HE(J)+HE(J+1))
	aelec(J)=0.5d0*(elec(J)+elec(J+1))
	aH2O(J)= 0.5d0*(H2O(J)+H2O(J+1))
	aXNH3(J)=0.5d0*(XNH3(J)+XNH3(J+1))
	aH(J)=   0.5d0*(H(J)+H(J+1))
	aHmin(J)=0.5d0*(Hmin(J)+Hmin(J+1))
	axNa(J)= 0.5d0*(xNa(J)+xNa(J+1))
	axK(J)=  0.5d0*(xK(J)+xK(J+1))
	aH2S(J)= 0.5d0*(H2S(J)+H2S(J+1))
	aRb(J)=  0.5d0*(RB(J)+RB(J+1)) !nb
	aCs(J)=  0.5d0*(CS(J)+CS(J+1)) !nb
	aFeH(J)= 0.5d0*(FEH(J)+FEH(J+1)) !nb
	aCrH(J)= 0.5d0*(CRH(J)+CRH(J+1)) !nb
	aTiO(J)= 0.5d0*(TIO(J)+TIO(J+1)) !nb
	aVO(J)=  0.5d0*(VO(J)+VO(J+1)) !nb
        
C	call getopac(tbar,pbar)

C Old method that uses pre-mixed opacity tables
C This next part is because Freedman's opacities do not include methane in the optical
C Need to include it via this file read
C        OPEN (80,FILE='tpint75.out')
C          DO K=1,NSPECV
C            READ (80,*) DUM, XKAPCH4(K,J)
C            IF (K.LE.790) XKAPCH4(K,J)=XKCH4(K)
C	    IF (K.gt.790) XKAPCH4(K,J)=0.d0
C          ENDDO
C          DO K=1,NSPECV
C            READ (80,*) DUM, XKAPCO(K,J)
C          ENDDO
C          DO K=1,NSPECV
C            READ (80,*) DUM, XKAPH2O(K,J)
C          ENDDO
C          DO K=1,NSPECV
C            READ (80,*) DUM, XKAPNH3(K,J)
C Nspecv goes from low wavelength to high -- freedman is opposite			
C	    xfreed(J,K)=xreturn(NSPECV-K+1)  !layer, freq
C
C          ENDDO
C     CLOSE (80)

C     Read in new opacity files constructed by Caroline Morley (Fall 2015)

        if (x736(j).lt.10.0) then 
			write(charJ736,'(I1)') int(x736(j))
       else if (x736(j).lt.100.0)  then 
			write(charJ736,'(I2)') int(x736(j))
       else if (x736(j).lt.1000.0)  then 
			write(charJ736,'(I3)') int(x736(j))
		endif

       if (x1060(j).lt.10.0) then 
			write(charJ1060,'(I1)') int(x1060(j))
       else if (x1060(j).lt.100.0)  then 
			write(charJ1060,'(I2)') int(x1060(j))
       else if (x1060(j).lt.1000.0)  then 
			write(charJ1060,'(I3)') int(x1060(j))
       else if (x1060(j).lt.1061.0)  then 
			write(charJ1060,'(I4)') int(x1060(j))
		endif

C       print*, 'chars: ',x736(j), charJ1060, charJ736


       OPEN (87,FILE='../opacities/h2o/h2o.'//trim(charJ1060)//'.2000table')
       OPEN (97,FILE='../opacities/karkoshka_ch4_2000.table')
c       OPEN (97,FILE='../opacities/ch4_2014/ch4_2014.'//trim(charJ1060)//'.2000table')
       OPEN (107,FILE='../opacities/co2_2013/co2_2013.'//trim(charJ1060)//'.2000table')
       OPEN (117,FILE='../opacities/nh3/nh3.'//trim(charJ1060)//'.2000table')
       OPEN (127,FILE='../opacities/na/na.'//trim(charJ736)//'.2000table')
       OPEN (137,FILE='../opacities/k/k.'//trim(charJ736)//'.2000table')
       OPEN (147,FILE='../opacities/h2s_2015/h2s_2015.'//trim(charJ1060)//'.2000table')
       OPEN (157,FILE='../opacities/rb/rb.'//trim(charJ736)//'.2000table')
       OPEN (167,FILE='../opacities/cs/cs.'//trim(charJ736)//'.2000table')
       OPEN (177,FILE='../opacities/feh/feh.'//trim(charJ1060)//'.2000table')
       OPEN (187,FILE='../opacities/crh/crh.'//trim(charJ1060)//'.2000table')
       OPEN (197,FILE='../opacities/tio/tio.'//trim(charJ1060)//'.2000table')
       OPEN (207,FILE='../opacities/vo/vo.'//trim(charJ1060)//'.2000table')
       READ (87,*)  header
       READ (87,*)  header
       READ (97,*)  header
       READ (97,*)  header
       READ (107,*)  header
       READ (107,*)  header
       READ (117,*)  header
       READ (117,*)  header
       READ (127,*)  header
       READ (127,*)  header
       READ (137,*)  header
       READ (137,*)  header
       READ (147,*)  header
       READ (147,*)  header
       READ (157,*)  header
       READ (157,*)  header
       READ (167,*)  header
       READ (167,*)  header
       READ (177,*)  header
       READ (177,*)  header
       READ (187,*)  header
       READ (187,*)  header
       READ (197,*)  header
       READ (197,*)  header
       READ (207,*)  header
       READ (207,*)  header
 

       DO K=1,NSPECV
           READ (87,*) DUM, XKAPH2O(K,J) 
           READ (97,*) DUM, XKAPCH4(K,J)
           READ (107,*) DUM, XKAPCO2(K,J)
           READ (117,*) DUM, XKAPNH3(K,J)
           READ (127,*) DUM, XKAPNA(K,J)
           READ (137,*) DUM, XKAPK(K,J)
           READ (147,*) DUM, XKAPH2S(K,J)
           READ (157,*) DUM, XKAPRB(K,J)
           READ (167,*) DUM, XKAPCS(K,J)
           READ (177,*) DUM, XKAPFEH(K,J)
           READ (187,*) DUM, XKAPCRH(K,J)
           READ (197,*) DUM, XKAPTIO(K,J)
           READ (207,*) DUM, XKAPVO(K,J)
c           if (k.eq.1999) print *, 'XKAPH2O(k,j) ', XKAPH2O(k,j)

       ENDDO
       CLOSE (87)
       CLOSE (97)
       CLOSE (107)
       CLOSE (117)
       CLOSE (127)
       CLOSE (137)
       CLOSE (147)
       CLOSE (157)
       CLOSE (167)
       CLOSE (177)
       CLOSE (187)
       CLOSE (197)
       CLOSE (207)
        
	
      ENDDO
      CLOSE (20)
c====================END READ IN LOOPS==================================


C====================NOW START ADDING OPACITY============================
C Start adding up all the layer tau, g, and wbar's
      DO 100 K=1,NSPECV
C ZERO THE COLUMN OPTICAL DEPTHS OF EACH TYPE
C THE OPTICAL DEPTH OF THE TOP OF THE MODEL
C MAY NOT BE ZERO 
      TAURV(K)=0.
      TAUHV(K)=0.
      TAUCV(K)=0.
      TAUCLV(K)=0.
      TAUGV(K)=0.
      TAUC2V(K)=0.
      tauhc(k)=0.
      ICLDFL = 0
      DO 100 J=1,NLAYER

        TBAR=0.5d0*(TEMP(J)+TEMP(J+1))
        PBAR=DSQRT(PRESS(J)*PRESS(J+1))

C      Get Tristan's h-minus results
      call trist(tbar,pbar,hminus,h2minus)

        BMU=0.5d0*(XMU(J+1)+XMU(J))
		
C Old version of COEF1
C       COEF1=RGAS*273.15**2*.5d5* (PRESS(J+1)**2 - PRESS(J)**2)
C    & /(1.01325**2 *EFFG(Z(J))*TBAR*BMU)

C    McKay 1/96 derivation:
       ACOEF = (TBAR/(TEMP(J)*TEMP(J+1)))*
     & (TEMP(J+1)*PRESS(J+1) - TEMP(J)*PRESS(J))/(PRESS(J+1)-PRESS(J))
       BCOEF = (TBAR/(TEMP(J)*TEMP(J+1)))*
     & (TEMP(J) - TEMP(J+1))/(PRESS(J+1)-PRESS(J))
C New COEF1	 
        COEF1=RGAS*273.15**2*.5E5* 
     & (ACOEF* (PRESS(J+1)**2 - PRESS(J)**2) + 
     & BCOEF*(2./3.)*(PRESS(J+1)**3 - PRESS(J)**3) )
     & /(1.01325**2 *EFFG(Z(J))*TBAR*BMU)

C NOW COMPUTE TAUGAS DUE TO THE PIA TERM 
C We arrange for there to always be a PIA term for each wavelength,
C so we don't have to test.

          TAUGAS=0.0d0
C Old version
C           CALL PIA(K,TBAR,PHH,PHE,PH,PHC)
            CALL PIAN(K,TBAR,PHH,PHE,PH,PHC)
			
            TAUGAS=COEF1*
     &       (aH2(J)*aH2(J)*PHH + aH2(J)*aCH4(J)*PHC 
     &         + aH2(J)*aHE(J)*PHE + aH2(J)*aH(J)*PH)


C KCahoy: Deleted the H2-, H- part that used to be below, since
C Freedman's table includes all H and e related opacities. 
C JJF's version may still have this part commented out in it above lines: 
c "Since Richard includes all H-related and e-related opacity in the tables, we can zero
c  out the contributions of these species here.  JJF 3/25/08
c	        taugas = taugas! + tauhmbf + tauhmff + tauh2m"
c !nb
c=========ADD BACK IN H- OPACITY STUFF =================================
              tauhmbf = 0.d0

c   tauhmbfold = 0.d0
             tauhmff = 0.d0
             tauh2m = 0.d0

      if (tbar.gt.600) then
c        It is hot enough to have electrons 
      if (wlnv(k).lt.1.642) then

c        and we are < 1.642 um photodetachemnt threshold
c        we get the the bf H- continuum
         call opa_hmbf(wnov(k),sbf)

c          this uses tristan's abundance table
c        tauhmbfold = sbf*hminus*h2(j)*COLDEN(J)/(xmu(j)*amu)
c          this uses Katharina's abundance table
       tauhmbf = sbf*Hmin(j)*COLDEN(J)/(xmu(j)*amu)

          endif
c         Then we get the H- continuum ff opacity

       call opa_hmff(wnov(k),tbar,sff_hm)

           tauhmff = pbar*1.e6*H(J)*elec(J)*sff_hm*colden(j)/
     &                 (tbar*xmu(j)*xk_boltz*amu)
c         Then we get the H2- ff continuum
c            print *, tbar,wnoi(k),h2min
            call opa_tab_h2mff(tbar,wnov(k),h2min)

            tauh2m = pbar*1.e6*h2(j)*elec(j)*h2min*
     &               colden(j)/(xmu(j)*amu)

        print *, j, k, taugas, tauhmbf, tauhmff, tauh2m
         taugas = taugas + tauhmbf + tauhmff + tauh2m
      endif     

c=======================================================================

C Return values for kappa and create taugas
C           TAUFRE = colden(j)/(axmu(j)*amu) * xfreed(J,K)  ! Freedman opacity
           TAUCH4=COLDEN(J)*aCH4(J)*xkapch4(k,j)*(16.04/axmu(j))
           TAUH2O=COLDEN(J)*aH2O(J)*xkaph2o(k,j)*(18.02/axmu(j))
           TAUNH3=COLDEN(J)*aXNH3(J)*xkapnh3(k,j)*(17.03/axmu(j))
           TAUCO2=COLDEN(J)*aCO2(J)*xkapco2(k,j)*(44.01/axmu(j))
      TAUNA=COLDEN(J)*axNa(J)*xkapna(k,j)*(22.9898/axmu(j))
           TAUK=COLDEN(J)*axK(J)*xkapk(k,j)*(39.0983/axmu(j))
      TAUH2S=COLDEN(J)*aH2S(J)*xkaph2s(k,j)*(34.08/axmu(j))
c     These below can be commented out if doing retrieval studies of cooler
c     Jupiter-like planets 
      TAURB=COLDEN(J)*aRb(J)*xkaprb(k,j)*(85.4678/axmu(j)) !nb
          TAUCS=COLDEN(J)* aCs(J)*xkapcs(k,j)*(132.90545/axmu(j)) !nb
      TAUFEH=COLDEN(J)* aFeH(J)*xkapfeh(k,j)*(56.85/axmu(j)) !nb
      TAUCRH=COLDEN(J)* aCrH(J)*xkapcrh(k,j)*(53.00/axmu(j)) !nb
          TAUTIO=COLDEN(J)* aTiO(J)*xkaptio(k,j)*(64.06/axmu(j)) !nb
      TAUVO=COLDEN(J)* aVO(J)*xkapvo(k,j)*(66.94/axmu(j)) !nb

           TAUGAS=TAUGAS+(TAUCH4+TAUH2O+TAUNH3+TAUCO2+TAUNA+TAUK+TAUH2S+
     &          TAURB+TAUCS+TAUFEH+TAUCRH+
     &          TAUTIO+TAUVO)
C           if (k.eq.1000)  print *, j, k, taugas
C     IF (K.LE.790) THEN
c JJF: making this array only 790 (was 1150) lines long cuts if off for wavelengths
c longer than 1.00 (was 1.616) microns.  At the longer wavelengths we use the CH4 
c opacity that is premixed in R. Freedman's tables.
C			TAUMETH=COLDEN(J)*GAS1(J)*xkapch4(k,j)
C		   ELSE
C			TAUMETH = aCH4(J) * COLDEN(J)/(aXMU(J)*AMU) * XKAPCH4(K,J)
C           ENDIF
C We only include the the *methane* opacity "by hand."  JJF, 3/25/08
C           TAUGAS=TAUGAS+TAUMETH+TAUFRE
c K. Cahoy temporary edit to divide opacities for Mark            
c           TAUGAS=TAUGAS/4.0

C Column (vertical) tau
           TAUGV(K)=TAUGV(K)+TAUGAS
C           if (k.eq.1000)  print *, j, k, TAUGV(K)

C Output at nlayer check (old)
C           if (j.eq.nlayer) then
C             x (84,939) wlnv(k),taugv(k),tauhc(k)
C           endif
C  939 format(1x,f7.4,2e10.3)

C HAZE
C  The HAZE properties are from Rages and the optical calc at the
C  start of this program.  We only use the haze if there is no cloud
C  at that level.  We assume that if there is a cloud the haze particles
C  will have been used as condensation nuclei.  There  is thus also no
C  haze below the cloud either.  Use ICLDFL to keep track of that.
C  But below the cloud we have an absorber in the clear region
C  above the lower cloud.  We only turn on this absorber if the
C  switch ICLEAR has been set equal to 1 in the settings page. 

C       Here we are above the cloud.  FHZE is a factor that changes
C       the haze tau.  It is used as a sensitivity probe.
C       IF ((XNCLD(J).EQ.0.0).AND.(ICLDFL.EQ.0)) THEN
C         fhz = dmin1(fhze,press(j)/0.05)
C         TAEROS=FHZE*TAERH(J,K)
C        QSCTHT = QSCTH(J,K)
C       Here we are in or below the cloud.
C       ELSE
C        ICLDFL = 1
        TAEROS=0.0
        QSCTHT = 0.0
C      ENDIF

C RAYLEIGH
C RAYLEIGH SCATTERING STRAIGHT FROM HANSEN AND TRAVIS, SEE NOTES
C RATIO OF THE LAYER COLUMN NUMBER TO THE TOTAL COLUMN NUMBER ON EARTH. CM-2

C Old version
C     TAURAY=.18*(COLDEN(J)*28.9/(XMU(J)*1013.25))*
C    &(.008569/WLNV(K)**4)*(1.+.0113/WLNV(K)**2+.00013/WLNV(K)**4)

        cfray = 32.*pi**3*1.e21/(3.*2.687e19)
        TAURAY = 0.0
        gasmtx(1) = aH2(J)
        gasmtx(2) = aHE(J)
        gasmtx(3) = aCH4(J)
        cold = COLDEN(J)/(aXMU(J)*AMU)
      do 2940 nn =1,3
        tec = cfray*(dpol(nn)/wlnv(k)**4)*(gnu(1,nn)+gnu(2,nn)/
     1                       wlnv(k)**2)**2
        taur(nn) = COLD*gasmtx(nn)* tec*1.e-5/XN0
        TAURAY = TAURAY + taur(nn)
2940  continue

C Raman scattering switch (turn on)
        wray=dmin1(ramalb(wlnv(k)),0.999999999d0)

C Raman scattering switch (turn off)
C        wray=0.99999999999d0

C CLOUD 
C COMPUTE TAU CLOUD
C Examples for two different layers of cloud, currently off

      TAUCLD=0.0
      tausctc=0.0
      wcld=0.0
      cbarc(k)=0.0

c  Top cloud
c      if (j.eq.13) then
c        taucld=2.0
c        tausctc=2.0
c        wcld=0.999
c        cbarc(k)=0.0
c      endif

c  Bottom cloud
c      if (j.eq.22) then
c        taucld=5.0
c        tausctc=5.0
c        wcld=0.999
c        cbarc(k)=0.0
c      endif

      if (wcld.gt.1.0) wcld=0.99999999

      TAURV(K)=TAURV(K)+TAURAY
      TAUHV(K)=TAUHV(K)+TAEROS
      TAUCV(K)=TAUCV(K)+TAUCLD

C Clear region absorber
C TAUCLV is for the clear region absorber
	  TAUCLR=0.d0
      TAUCLV(K)=TAUCLV(K)+TAUCLR

C H2S cloud
C TAU2 is for the H2S cloud

	  tau2(j)=0.0
      TAUC2V(K)=TAUC2V(K)+TAU2(J)

C TAUGAS
C LOOP OVER THE NTERMS COMPUTE TOTAL TAUGAS
C COMPUTE THE AVERAGE COSBAR AND WBAR

C Weight by the scattering, not total opacity (except for H2S cloud)
C	   tauscth = qscth(j,k)*xnumb(j)*xsct(j)
C	   tausctc = qsctc(k)*xncld(j)*xsctc
C	   tausct2 = tau2(j)
	  tauscth=0.0
C        tausctc=taucld
        tausct2=0.0

C DEFAULT PARAMETERS		
C Default COSBV (asymmetry factor, "g")
        COSBV(J,K)= ( CBARH(J,K)*tauscth + CBARC(K)*tausctc + 
     &     CBAR2(J)*tausct2 + gcx(J,K)*taucx(J,K))
     &    /(tauscth+tausctc+TAURAY+tausct2+TAUCLR+taucx(J,K))
C Default GCOS2 
        GCOS2(j,k) = 0.5*TAURAY/(TAURAY + tauscth +tausctc +tausct2+taucx(J,K))
C Default DTAUV
        DTAUV(J,K)=TAUGAS+TAEROS+TAURAY+TAUCLD+TAU2(J)+TAUCLR+taucx(J,K)
C Default WBARV
        WBARV(J,K)=(
     &     tauray*wray 
     &     + wcld*TAUCLD +TAU2(J)*WBAR2(J) +wbarcx(J,K)*taucx(J,K))
     &     /(TAUGAS+TAEROS+TAURAY+TAUCLD+TAU2(J)+taucx(J,K))

C RAYLEIGH PARAMETERS (use with GINT2 or GINT3 in GEOM)
C Zero COSBV for Rayleigh test case
C        COSBV(J,K)=0
C Constant GCOS2 for Rayleigh test case
C        GCOS2(J,K)=0.5
C DTAUV Rayleigh only
C        DTAUV(J,K)=TAURAY
C WBARV Rayleigh
C        WBARV(J,K)=0.9999
c        WBARV(J,K) = 0.9500
        

C ISOTROPIC PARAMETERS also set RSFV to 0.99 (two places in egpalb.f)		
C Zero COSBV (asymmetry factor, "g") for isotropic case
c        COSBV(J,K)= 0.0
C Zero GCOS2 for isotropic case (used in PUU0 with GINT2, GINT3)
c        GCOS2(j,k) = 0.0
C Default DTAUV
C      DTAUV(J,K)=TAUGAS+TAEROS+TAURAY+TAUCLD+TAU2(J)+TAUCLR+taucx(J,K)
C      DTAUV(J,K) = TAURAY
C      DTAUV(J,K) = 10.d0        
C Constant WBARV (set to 0.99, 0.95, etc?)
c        WBARV(J,K) = 0.9999
c        WBARV(J,K) = 0.9500
        
        
        

C WATER: IN THE CASE OF WATER, uncomment below
C        if (cosbv(j,k).gt.0.70) cosbv(j,k)=0.70

C Print to screen if needed for testing
C        write(*,333) J,K,dtauv(J,K),wbarv(J,K),gcos2(J,K),cosbv(J,K),taugas
C        write(*,333) J,K,TAUGAS,TAEROS,TAURAY,TAUCLD,TAU2(J),TAUCLR,taucx(J,K)
C 333 	FORMAT(2I5,7E13.4)
  909  CONTINUE
  100  CONTINUE


C TOTAL EXTINCTION OPTICAL DEPTHS
         DO 119 K=1,NSPECV
C LOOP OVER NTERMS      
           TAUV(1,K)=0.0
C FLAG VARIABLE
           FLAG=1.0
           banner=1.0
           DO 119 J=1,NLAYER
           TAUV(J+1,K)=TAUV(J,K)+DTAUV(J,K)
             IF ((TAUV(J+1,K).GE.0.5).AND.(FLAG.EQ.1.0)) THEN
               FLAG=0.0
               NOPTHCK(K)=J+1
             ENDIF
             if ((tauv(J+1,K).GE.1.0).AND.(BANNER.EQ.1.0)) THEN
               BANNER=0.0
               LAYERBOT(K)=J+1
             ENDIF
 119     CONTINUE


	  IPRINT = 3
C These are output files, but they need to exist first, can be empty to start
C......................................................................
C   	  OPEN(17, FILE='00_SPEC_INT_DATA.TXT', STATUS='OLD')
C	  OPEN(18, FILE='00_TAUV_NLEVEL.TXT', STATUS='OLD')

	    WRITE(TFNAME, '(a,I3.3,a)') '00_Spec_Int_Data_', idnint(PLPHA*RTOD), '.txt'
        OPEN(17, FILE=TFNAME, STATUS='REPLACE')
C Turned off for speed
C	    WRITE(TFNAME2, '(a,f10.5,a)') '00_TAUV_NLevel', PLPHA, '.txt'
C        OPEN(18, FILE=TFNAME2, STATUS='REPLACE')

C        WRITE(TFNAME3, '(a,f10.5,a)') '00_TAUV_COLDEN', PLPHA, '.txt'
C        OPEN(19, FILE=TFNAME3, STATUS='REPLACE')

C......................................................................	  
      IF (IPRINT .GT. 2) THEN
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           WRITE (6,120)
		   WRITE (17,120)
  120      FORMAT(///'  OPTICAL CONSTANTS IN THE VISIBLE ')
           DO 200 K=1,NSPECV,1
c           WRITE (6,190)
		
	   WRITE (17,190)
c           WRITE (6,210)K,WLNV(K),WNOV(K)
c           WRITE (6,230)REALV(K),XIMGV(K)
	   WRITE (17,210)K,WLNV(K),WNOV(K)
           WRITE (17,230)REALV(K),XIMGV(K)
           DO 195 J=1,NLAYER
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		   

C CONTRIBUTION FUNCTION 
c Compute contribution function (see Griffith, Yelle, & Marley, 1998)
	hq=6.62d-27  !cgs units
	cq=2.998d10
	zkq=1.38d-16
	BB=2.d0*hq*cq**2 / (wlnv(k)/1.d4)**5 * 
     &	 exp(-hq*cq/((wlnv(k)/1.d4)*zkq*temp(J)))
	if(j.eq.nlayer) then
	chap=0.d0
	else
	chap1= -1.d0*(exp(-tauv(j+1,k)/0.5d0)-exp(-tauv(j,k)/0.5d0)) 
	chap2= (dlog10(press(j+1))-dlog10(press(j)))
	chap=chap1/chap2
	endif
	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		   
c           WRITE (6,220)temp(J), press(J), WBARV(J,K),COSBV(J,K)
c     &      ,DTAUV(J,K),TAUV(J,K),colden(J),TAURAY,BB*chap,xfreed(J,K)
	 
           WRITE (17,220)temp(J), press(J), WBARV(J,K),COSBV(J,K)
     &      ,DTAUV(J,K),TAUV(J,K),colden(J),TAURAY,BB*chap,xfreed(J,K)
	 
  195      CONTINUE
  
  
C           WRITE (18,240) TAUV(NLEVEL-1,K)/COLDEN(NLEVEL-1), 
C     &      TAUV(60,K)/COLDEN(60), TAUV(50,K)/COLDEN(50), 
C     &      TAUV(40,K)/COLDEN(40), TAUV(30,K)/COLDEN(30), 
C     &      TAUV(20,K)/COLDEN(20), TAUV(10,K)/COLDEN(20)

C           DO 272 J=1,NLAYER
C272        WRITE (18,240) TAUV(J,K)
C  
C           DO 273 J=1,NLAYER
C273        WRITE (18,240) PRESS(J)         

C           DO 275 J=1,NLAYER
C275        WRITE (19,241) TAUV(J,K)/COLDEN(J)

C           WRITE (18,240) TAUV(70,K), TAUV(69,K), TAUV(68,K), TAUV(67,K), 
C     &      TAUV(66,K), TAUV(65,K), TAUV(64,K), TAUV(63,K), TAUV(62,K), 
C     &      TAUV(61,K), TAUV(60,K), TAUV(59,K), TAUV(58,K), TAUV(57,K), 
C     &      TAUV(56,K), TAUV(55,K), TAUV(54,K), TAUV(53,K), TAUV(52,K),         
C     &      TAUV(51,K), TAUV(50,K), TAUV(49,K), TAUV(48,K), TAUV(47,K),   
C     &      TAUV(46,K), TAUV(45,K), TAUV(44,K), TAUV(43,K), TAUV(42,K),
c     &      TAUV(41,K), TAUV(40,K), TAUV(39,K), TAUV(38,K), TAUV(37,K),  
C     &      TAUV(36,K), TAUV(35,K), TAUV(34,K), TAUV(33,K), TAUV(32,K),
C     &      TAUV(31,K), TAUV(30,K), TAUV(29,K), TAUV(28,K), TAUV(27,K),  
C     &      TAUV(26,K), TAUV(25,K), TAUV(24,K), TAUV(23,K), TAUV(22,K),    
c     &      TAUV(21,K), TAUV(20,K), TAUV(19,K), TAUV(18,K), TAUV(17,K),
C     &      TAUV(16,K), TAUV(15,K), TAUV(14,K), TAUV(13,K), TAUV(12,K), 
C     &      TAUV(11,K), TAUV(10,K), TAUV(09,K), TAUV(08,K), TAUV(07,K), 
C     &      TAUV(06,K), TAUV(05,K), TAUV(04,K), TAUV(03,K), TAUV(02,K),
c     &      TAUV(01,K) 
                  		   
  200      CONTINUE
      END IF
	  
  210 FORMAT(1X,I4,F13.6,F10.2,F10.2,'-',F8.2,F10.3)
  190 FORMAT(1X//'  SNUM  MICRONS   WAVENU   INTERVAL    DELTA-WN')
  230 FORMAT(1X,'NREAL(LAYER)= ',1PE10.3,' NIMG(LAYER)= ',E10.3/
     &' TEMP    PRESSURE   WBAR     COSBAR     DTAU     TAU'
     & ,9X,'COLDEN     GAS    CONTRIB.    XFREED')
  220 FORMAT(F8.1,9(1X,G9.3))
C  240 FORMAT(70(G12.3))
  240 FORMAT(1(G12.3))
  241 FORMAT(1(G12.3))
	  CLOSE(17)
C	  CLOSE(18)
C      CLOSE(19)
	  
      RETURN
      END
C----------------------------------------------------------------------- 
C End OPTCV

C----------------------------------------------------------------------- 
C EFFG (GRAVITY, sort of obsolete since now included in file *.in)

      FUNCTION EFFG(Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /grav/ gin
C      DATA G/50.0/
      EFFG = gin
      RETURN
      END
C----------------------------------------------------------------------- 
C END EFFG

C----------------------------------------------------------------------- 
C HAZESCAT (Called in GEOMETRIC)

      SUBROUTINE HAZESCAT(IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'prog_params'

      COMMON /ATM/ Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
      COMMON /AERSOL/ RADIUS(NLAYER), XNUMB(NLAYER), SIGMA(NLAYER) 
     & , REALI(NSPECI), XIMGI(NSPECI), REALV(NSPECV), XIMGV(NSPECV)
      COMMON /AERSOLV/ TAERH(NLAYER,NSPECV),QSCTH(NLAYER,NSPECV),
     &         CBARH(NLAYER,NSPECV),XSCT(NLAYER)
      COMMON /AERSOLI/ TAERHI(NLAYER,NSPECI)
      COMMON /CLEAR/ CLA(NSPECV),ICLEAR
      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      COMMON /CONST/RGAS,RHOP,PI,SIGSB,RTOD 
      COMMON /ADJUST/ RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,RCLOUD
      REAL*8 NREAL,NIMAG

C     Rmin and Rmax are fractions of Rmean that the integration is
C       carried out to.
      write(*,*)'NLAYER=',NLAYER
      write(*,*)'NLEVEL=',NLEVEL
C Old versions
C	RMIN = 0.25
C	RMAX = 3.0

        rmin=0.0
        rmax=0.0

C     Ar and Br are Rages fits to imag. index, they have uncertainties of
C          +/- 1.1 and 2.2 respectively.  Br has units of micron^(-1)
C Old versions
C Ar = -2.1
C Br = -6.5
        ar=0.0
        br=0.0

C Not currently used		
C     BETAH is the clear region absorber factor, uncert is +/- 13
C     BETAH = 10.
C       print *,'Input Betah Clear region absorber factor (0/10/23)'
C       accept *,BETAH 


C     To choose a logarithmic distribution in FUNC, we need to pick
C       NF = 5
        nf=0.0
C     The parameter before NF is not used in a log dist., so we set  
        B = 0
		
C   Epsq is smallest q considered in integration, Xfunny is ?
C	Xfunny = 2.0
C   Epsq = 1.e-8
        xfunny=0.0
        epsq=0.0
		
C  First we treat the haze in the visible.
	DO 10 J=1,NLAYER
	DO 10 I=1,NSPECV
C        In this version we don't save QEXTH since we calculate
C        all the necessary parameters just once here.
C	 QEXTH(J,I) = 0.0
   	 QSCTH(J,I) = 0.0 
	 CBARH(J,I) = 0.0
10    CONTINUE
C......................................................................	
      OPEN (55,FILE='HAZEINFO.',STATUS='OLD')
C......................................................................		  

      DO 200 J=1,NLAYER
C Old version
C        XSCT(J) = PI * (RADIUS(J)**2) * 1E-8
          xsct(j)=0.0
		  
	IF (XNUMB(J).EQ.0.0) GOTO 200
C Old versions
C 	NREAL = 1.42
C IF (PRESS(J+1).GT.2.4E-3) NREAL = 1.33
C IF (PRESS(J+1).GT.12.0E-3) NREAL = 1.44
        nreal=0.0

	DO 100 I=1,NSPECV
	
C Old version	
C      NIMAG = 100.*DEXP (AR + BR*WLNV(I)) 
C      CALL SCAMIE( WLNV(I),   NREAL,  NIMAG, RMIN, RADIUS(J),
C    &                   RMAX,   B, NF, SIGMA(J),
C    &                   Epsq,   Xfunny,
C    &                   QEXTH, QSCTH(J,I), ALBED, BETA)
C        CBARH(J,I) = BETA
C        TAERH(J,I) = QEXTH*XNUMB(J)*XSCT(J)
C       write (55,7336) j,wlnv(i),qscth(j,i),cbarh(j,i),taerh(j,i)

	   read  (55,*)jdum,wdum,qscth(j,i),cbarh(j,i),taerh(j,i)
       jdum=0.0
       wdum=0.0
       qscth(j,i)=0.0
       cbarh(j,i)=0.0
       taerh(j,i)=0.0
100   CONTINUE
200   CONTINUE 
C Old version
C 7336   FORMAT (2X,I3,4E15.5)
       CLOSE(55)

C Not currently used        
C     Now we do some of the calculation for the absorber in the
C     clear region below the methane cloud.  See Rages et al and
C     notes.
C     print *,'clear region absorber imaginary multiplier (1;10=Tholin)'	
C     accept *, ccm
C     DO 700 I=1,NSPECV
C        X = ccm* DEXP (AR + BR*WLNV(I)) / WLNV(I) 
C        CLA(I) = PI * BETAH * X / 4.46E3
C 700   CONTINUE
      RETURN
      END
C----------------------------------------------------------------------- 
C END HAZESCAT

C----------------------------------------------------------------------- 
C HAZESET (Called in GEOMETRIC)

      SUBROUTINE HAZESET(IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'prog_params'     

      COMMON /ATM/ Z0(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
      COMMON /AERSOL/ RADIUS(NLAYER), XNUMB(NLAYER), SIGMA(NLAYER) 
     & , REALI(NSPECI), XIMGI(NSPECI), REALV(NSPECV), XIMGV(NSPECV)
	DIMENSION PR(100),DN(100),DM(100),Z(100),DD(100),DS(100)
	DIMENSION ZM(100),DDM(100),DSM(100)
C  This subroutine reads in the Haze information and readjusts the
C  haze to fit on the pressure grid used by the program.
      write(*,*)'NLAYER=',NLAYER
      write(*,*)'NLEVEL=',NLEVEL
C......................................................................	
	OPEN (66,FILE='RAGES.HAZ',STATUS='OLD')
C......................................................................	
	READ (66,*) NX
	DO 50 I=2,NX+1
	  READ (66,*) PR(I),PB,Z(I),DN(I),DD(I),DS(I)
C	 mbars to bars:
          PR(I) = PR(I)/1000.D0
50	CONTINUE
	  PR(NX+2) = PB/1000.
	  DN(NX+2) = 0.0
	  Z(NX+2) = Z(NX+1) - (Z(NX) - Z(NX+1))
	  Z(NX+3) = Z(NX+2) - 5.0
          PR(NX+3) = PR(NX+2) + 5.
	  DN(NX+3) = 0.0
	  DN(1) = 0.0
	  PR(1) = 1e-9
          Z(1) = 1000.
          DO 75 I=1,NLEVEL
           DM(I) = 0.0
           DSM(I) = 0.0
           DDM(I) = 0.0
75         ZM(I) = 0.0
 
	CALL HAZES (PR,DN,DM,Z,ZM,DD,DDM,DS,DSM,NX+2)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WRITE(*,105)
105	FORMAT ( /,1X,'LEVEL',4X,'P',5X,'Rbar',6X,
     1           'SIGMA',4X,'N/cm^2')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	DO 100 I=1,NLEVEL-1
 	 XN = 0.
C        IF (PRESS(I+1).LT.1.41) THEN
	 IF (ZM(I).NE.0.0) XN = DM(I)/(ZM(I)*1.E5)
         RADIUS(I) = DDM(I)/2. 
	 SIGMA(I)= DSM(I)
	 XNUMB(I) = DM(I)
C        ENDIF
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 WRITE (*,111) I,PRESS(I),RADIUS(I),SIGMA(I),XNUMB(I)
100	CONTINUE
111	FORMAT (1X,I3,3F9.4,2X,1E11.4)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	END	
C----------------------------------------------------------------------- 
C END HAZESET

C----------------------------------------------------------------------- 
C CLOUDSCAT (Called in subroutine CLD)

      SUBROUTINE CLOUDSCAT 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'prog_params'
      COMMON /ADJUST/ RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,RCLOUD
      COMMON /CLOUDV/ QEXTC,QSCTC(NSPECV),CBARC(NSPECV)
      COMMON /CLOUDI/ QABSCI(NSPECI)
      COMMON /CLOUD/ RADCLD(NLAYER), XNCLD(NLAYER), XSCTC
     & , RCLDI(NSPECI), XICLDI(NSPECI), RCLDV(NSPECV), XICLDV(NSPECV)
C      COMMON /SPECTI/ BWNI(NSPC1I),WNOI(NSPECI),DWNI(NSPECI)
C     &           ,WLNI(NSPECI)
      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      COMMON /CONST/RGAS,RHOP,PI,SIGSB,RTOD
      COMMON /TESTS/FHZE,FRAY,FCLDC
      
C   The total optical depth of the cloud will be constrained
C   to equal that found by Rages et al. 0.66+/-0.18 at 22S, 2.4+/-0.3 65S
C   We use the average of the two values.

C Old version
C       print *,'INPUT TAUFAC CH4 cloud optical depth (.5/1.5/2.7)'
C       accept *, taufac

       TAUFAC = 0.0

C   RCLOUD is cloud particle size in microns

C Old version
C      print *,'Input Rcloud  CH4 cloud size microns (1)'
C      accept *, rcloud

       RCLOUD = 1.0

C Old version	   
C       XSCTC = (RCLOUD**2) * PI * 1.E-8 

         xsctc=0.0

C   Ar and Br are Rages fits to haze imag. index, they have uncertainties of
C      +/- 1.1 and 2.2 respectively.  Br has units of micron^(-1)

C Old version
C      print *,'input Ar (-3.2/-2.1/-1.0), Br (-8.7/-6.5/-4.3)' 
C      accept *,Ar,Br
        Ar = -2.1
        Br = -6.5
		
      PI = 3.141592654
	   
C     Methane cloud properties from Rages et al 1991
      QEXTC = 2.0
      gcloud = 0.34

C Old version
C     fh = 0.72 +/- 0.68 at 22.5 deg S    0.46+/- 0.42 at 65S    
C     FH = 0.72
C      print *, 'Fh (.04/0.72/1.4)'
C      accept *,FH 

C Old version
C      print *,'Cloud contaminant imag index multiplier (1;10=Tholin)'
C      accept *,ccm
	ccm=1.d0

      QRC = 0.05 + ((1.285 - 1.0)**2/(1.285 + 1.)**2)

      DO 100 J=1,NSPECV
	XLAM = WLNV(J)
        XNI =ccm* DEXP(Ar + Br*XLAM)
        XK = PI * XNI / XLAM
	QSCTC(J) = 1.0 + QRC + (1.0 - QRC)*DEXP(-2.0*FH*XK)
        CBARC(J) = fcldc*gcloud
100   CONTINUE

      RETURN
      END
C----------------------------------------------------------------------- 
C END CLOUDSCAT

C----------------------------------------------------------------------- 
C PCH4 (Called in GEOMETRIC, used for info printout)

      FUNCTION PCH4(T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA ASAT,BSAT,PSAT,TSAT/145.99,-0.1259,0.11719,90.68/
      DATA R,XMV/1.9872,16/
C THIS SUBROUTINE RETURNS THE VAPOR PRESSURE OF METHANE
C OVER THE TEMPERATURE INTERVAL (74-97K) IN UNITS OF BARS

C Old version
C IT IS BASED ON DATA FROM THE MATHASEN GAS DATA BOOK P 351
C FITTED BY CP MCKAY 10/85
C     PCH41=3.4543E4 * DEXP(-1145.705/T)

C New version
C    Now we have switched to using Pollack's expression.  It does
C    not differ significantly (but does not have the T limits)
C    See notes.
      XL = ASAT + BSAT*(T-TSAT)
      PCH4 = PSAT * DEXP(-(XMV*XL/R)*((1./T)-(1./TSAT)))
      RETURN
      END
C----------------------------------------------------------------------- 
C END PCH4

C----------------------------------------------------------------------- 
C FUNC (descriptive... used in SCAMIE subroutine, Mie scattering)

      FUNCTION FUNC ( A, RMEAN, B, PO, NF )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Options (50) and (60) added for Venus cloud layers. 3-16-90 --Brad
      DATA    C1  /     10.0      /
      DATA    C2  /     50.0      /
      DATA    PI  /3.14159265358/
      GO TO ( 10,20,30,40 ,50,60) NF
   10     ALPH = (1.0-3.0*B)/B
          FUNC = A**ALPH*DEXP(-A/(RMEAN*B))

      RETURN

   20     Z = (B/PO)*((A/RMEAN)**PO)
          FUNC = (A**B)*DEXP(-Z)

      RETURN

   30   IF ( A .LT. PO ) GO TO 125

          IF ( A .GE. PO .AND. A .LE. C1*RMEAN ) GO TO 135

            IF ( A .GT. C1*RMEAN .AND. A .LE. C2*RMEAN ) GO TO 145

C AT THIS  POINT -A- IS GREATER THAN C2 * RMEAN, SO...
          FUNC = 0.0

      RETURN

  145     FUNC = ((C1*RMEAN/PO)**2)*((PO/A)**(B+2.0))

      RETURN

  135     FUNC = (PO/A)**B

      RETURN

  125     FUNC = 1.0

      RETURN
  40      Z1   = (A/RMEAN)**B
          Z2   = PO*Z1/B

          FUNC = A**(PO)*DEXP(-Z2)
      RETURN

C     Cloud size distribution for Venus.  Mode (1) Use NF = 5
C       SigG  =  width param for distribution.
C       D     =  particle size parameter based on diameter
C       Dg    =  modal particle diameter
  50      Fac1 = (2.0*PI)**(-0.5)
          SigG = PO
          D    = 2.0*A
          Dg   = 2.0*RMEAN
          Fac2 = Fac1 / (D * LOG(SigG))
          FUNC = Fac2 *DEXP(-0.5*DABS((DLOG(D/Dg)/DLOG(SigG)))**2)
      RETURN

C     Cloud size distribution for Venus.  Mode (2) Use NF = 6
   60     SigG = PO
          D    = 2.0*A
          Dg   = 2.0*RMEAN
          Fac1 = 1.0 / (((2.0*PI)**0.5)*SigG)
          FUNC = Fac1 * DEXP(-0.5*((D-Dg)/SigG)**2)
      RETURN
      END
C----------------------------------------------------------------------- 
C END FUNC (distribution used in SCAMIE)

C----------------------------------------------------------------------- 
C DSOLVER (Called in GFLUXV, GINT, GINT2, GINT3, etc. from GEOMETRIC)
C
C THIS SUBROUTINE SOLVES FOR THE COEFFICIENTS OF THE    
C TWO STREAM SOLUTION FOR GENERAL BOUNDARY CONDITIONS   
C NO ASSUMPTION OF THE DEPENDENCE ON OPTICAL DEPTH OF   
C C-PLUS OR C-MINUS HAS BEEN MADE.                      
C NL     = NUMBER OF LAYERS IN THE MODEL                
C CP     = C-PLUS EVALUATED AT TAO=0 (TOP)              
C CM     = C-MINUS EVALUATED AT TAO=0 (TOP)             
C CPM1   = C-PLUS  EVALUATED AT TAOSTAR (BOTTOM)        
C CMM1   = C-MINUS EVALUATED AT TAOSTAR (BOTTOM)        
C EP     = EXP(LAMDA*DTAU)                              
C EM     = 1/EP                                         
C E1     = EP + GAMA *EM                                
C E2     = EP - GAMA *EM                                
C E3     = GAMA*EP + EM                                 
C E4     = GAMA*EP - EM                                 
C BTOP   = THE DIFFUSE RADIATION INTO THE MODEL AT TOP  
C BSURF  = THE DIFFUSE RADIATION INTO THE MODEL AT      
C          THE BOTTOM: INCLUDES EMMISION AND REFLECTION 
C          OF THE UNATTENUATED PORTION OF THE DIRECT    
C          BEAM. BSTAR+RSF*FO*EXP(-TAOSTAR/U0)          
C RSF    = REFLECTIVITY OF THE SURFACE                  
C XK1    = COEFFICIENT OF THE POSITIVE EXP TERM         
C XK2    = COEFFICIENT OF THE NEGATIVE EXP TERM         

      SUBROUTINE DSOLVER(NL,GAMA,CP,CM,CPM1,CMM1
     ,,E1,E2,E3,E4,BTOP,BSURF,RSF,XK1,XK2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C DOUBLE PRECISION VERSION OF SOLVER
      PARAMETER (NMAX=301)	  
C      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION GAMA(NL),CP(NL),CM(NL),
     ,CPM1(NL),CMM1(NL),XK1(NL),XK2(NL)
     ,,E1(NL),E2(NL),E3(NL),E4(NL)
      DIMENSION AF(NMAX),BF(NMAX),CF(NMAX),DF(NMAX),XK(NMAX)
      L=2*NL
	  
C MIXED COEFFICENTS
C THIS VERSION AVOIDS SINGULARITIES ASSOC.
C WITH W0=0 BY SOLVING FOR XK1+XK2, AND XK1-XK2.
      AF(1)=0.0
      BF(1)=GAMA(1)+1.
      CF(1)=GAMA(1)-1.
      DF(1)=BTOP-CMM1(1)
      N=0
      LM2=L-2
C EVEN TERMS
      DO 10 I=2,LM2,2
          N=N+1
          AF(I)=(E1(N)+E3(N))*(GAMA(N+1)-1.)       
          BF(I)=(E2(N)+E4(N))*(GAMA(N+1)-1.)
          CF(I)=2.*(1.-GAMA(N+1)**2)
          DF(I)=(GAMA(N+1)-1.) * (CPM1(N+1) - CP(N))
     &          + (1.-GAMA(N+1))* (CM(N)-CMM1(N+1))
   10 CONTINUE
      N=0
      LM1=L-1
      DO 20 I=3,LM1,2
          N=N+1
          AF(I)=2.*(1.-GAMA(N)**2)
          BF(I)=(E1(N)-E3(N))*(1.+GAMA(N+1))
          CF(I)=(E1(N)+E3(N))*(GAMA(N+1)-1.)
          DF(I)=E3(N)*(CPM1(N+1) - CP(N)) 
     &         + E1(N)*(CM(N) - CMM1(N+1))
   20 CONTINUE
      AF(L)=E1(NL)-RSF*E3(NL)
      BF(L)=E2(NL)-RSF*E4(NL)
      CF(L)=0.0
      DF(L)=BSURF-CP(NL)+RSF*CM(NL)
      CALL DTRIDGL(L,AF,BF,CF,DF,XK)
C UNMIX THE COEFFICIENTS
      DO 28 N=1,NL
      XK1(N)=XK(2*N-1)+XK(2*N)
      XK2(N)=XK(2*N-1)-XK(2*N)
C NOW TEST TO SEE IF XK2 IS REALLY ZERO TO THE LIMIT OF THE
C MACHINE ACCURACY  = 1 .E -30
C XK2 IS THE COEFFICEINT OF THE GROWING EXPONENTIAL AND MUST
C BE TREATED CAREFULLY
      IF (XK2(N) .EQ. 0.0) GO TO 28
      IF (ABS (XK2(N)/XK(2*N-1)) .LT. 1.E-30) XK2(N)=0.0
   28 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------- 
C END DSOLVER 

C----------------------------------------------------------------------- 
C THOLIN (Not currently used anywhere... )

      SUBROUTINE THOLIN(WAVELN,XNR,XNI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(90),XN(90),XK(90)
      DATA W/
     &920.0000,850.0000,774.9000,688.8000,563.5000,387.4000,229.6000,
     &172.2000,140.9000,121.5000, 81.5700, 56.3500, 36.4600, 31.0000,
     & 22.1400, 18.2300, 17.7100, 14.4200, 12.9100, 11.7000, 11.0700,
     & 10.5100,  8.7310,  7.6530,  7.0440,  6.6660,  6.4570,  6.3260,
     &  5.9610,  5.8480,  5.7400,  5.4380,  5.2530,  5.1660,  4.8810,
     &  4.6260,  4.5920,  4.4920,  4.4280,  4.2170,  3.9480,  3.6680,
     &  3.4630,  3.2460,  3.0090,  2.9380,  2.8180,  2.7430,  2.6950,
     &  2.4220,  2.4030,  2.3930,  2.2140,  2.0190,  1.8730,  1.8230,
     &  1.8130,  1.8020,  1.3810,  1.3570,  1.1920,  1.1480,  1.0160,
     &  0.8731,  0.6888,  0.5635,  0.4428,  0.4133,  0.3874,  0.3542,
     &  0.3263,  0.2952,  0.2638,  0.2384,  0.1968,  0.1631,  0.1362,
     &  0.1215,  0.1181,  0.1159,  0.1097,  0.1016,  0.0925,  0.0800,
     &  0.0785,  0.0588,  0.0449,  0.0415,  0.0312,  0.0207/
      DATA XN/
     &2.170,2.170,2.160,2.160,2.150,2.120,2.070,2.040,2.030,2.020,1.930,
     &1.860,1.810,1.810,1.800,1.760,1.740,1.670,1.670,1.640,1.660,1.670,
     &1.710,1.720,1.690,1.690,1.640,1.580,1.430,1.440,1.480,1.550,1.580,
     &1.580,1.610,1.620,1.610,1.610,1.610,1.630,1.640,1.650,1.650,1.650,
     &1.610,1.590,1.580,1.590,1.600,1.620,1.620,1.620,1.630,1.630,1.630,
     &1.640,1.640,1.640,1.640,1.640,1.650,1.650,1.650,1.660,1.680,1.700,
     &1.720,1.690,1.660,1.630,1.640,1.660,1.680,1.680,1.660,1.650,1.700,
     &1.740,1.750,1.750,1.720,1.670,1.580,1.370,1.330,0.963,0.812,0.802,
     &0.850,0.920/
      DATA XK/
     &3.0E-3,1.0E-2,2.9E-2,4.7E-2,7.0E-2,1.0E-1,1.4E-1,1.6E-1,1.6E-1,
     &1.9E-1,2.1E-1,1.9E-1,1.5E-1,1.4E-1,1.8E-1,2.1E-1,2.1E-1,1.7E-1,
     &1.4E-1,9.7E-2,7.9E-2,7.5E-2,9.2E-2,1.3E-1,1.7E-1,2.2E-1,2.6E-1,
     &2.8E-1,1.5E-1,7.0E-2,2.9E-2,1.1E-2,8.7E-3,7.6E-3,1.0E-2,2.7E-2,
     &2.8E-2,1.4E-2,1.1E-2,1.0E-2,1.3E-2,2.1E-2,3.5E-2,5.6E-2,7.5E-2,
     &6.0E-2,2.4E-2,1.1E-2,4.1E-3,1.2E-3,8.5E-4,8.0E-4,8.9E-4,7.2E-4,
     &5.2E-4,4.4E-4,4.2E-4,4.0E-4,4.1E-4,4.2E-4,5.2E-4,6.4E-4,1.0E-3,
     &2.4E-3,8.8E-3,2.3E-2,6.0E-2,7.6E-2,9.1E-2,1.1E-1,1.3E-1,1.5E-1,
     &1.8E-1,2.1E-1,2.2E-1,2.4E-1,2.7E-1,3.7E-1,4.0E-1,4.3E-1,5.0E-1,
     &5.8E-1,6.7E-1,7.7E-1,7.7E-1,6.2E-1,3.8E-1,3.1E-1,1.4E-1,4.9E-2/
      XNR=XN(1)
      XNI=XK(1)
      IF (WAVELN .GT. W(1))  RETURN
      XNR=XN(90)
      XNI=XK(90)
      IF (WAVELN .LT. W(90)) RETURN
      DO 100 I=2,90
      IF (WAVELN .GT. W(I) ) GO TO 101
 100  CONTINUE
 101  CONTINUE
C ALL INTERPOLATION IS IN LOG LAMBDA 
      FACTOR= (dlog(WAVELN) - dlog(W(I)) ) / (dlog(W(I-1)) - dlog(W(I)))
	  
C REAL PART IS LINEARLY INTERPOLATED
      XNR=XN(I) + FACTOR*(XN(I-1) - XN(I))
	  
C IMAGINARY PART IS LOG INTERPOLATED
      XNI=dlog(XK(I)) + FACTOR*(dlog(XK(I-1)) - dlog(XK(I)))
      XNI=exp(XNI)
      RETURN
      END
C----------------------------------------------------------------------- 
C END THOLIN

C----------------------------------------------------------------------- 
C REFLIQ (Not currently used anywhere)

      FUNCTION REFLIQ(W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WAVENO(55),XIMG(55)
      DATA WAVENO/0., 10., 40., 60., 80., 100., 120., 140.,
     &160., 180., 200., 220., 240., 260., 280., 300., 325.,
     &340., 350., 360., 380., 400., 450., 480., 500., 520.,
     &540., 560., 600., 640., 684., 720., 760., 780., 800.,
     &840., 875., 920., 960., 1000., 1040., 1080., 1120., 1160.,
     &1200., 1220., 1240., 1260., 1280., 1300., 1320., 1340.,
     &1400., 1800., 2000./ 
      DATA XIMG/1.0E-5, 5.0E-5, 1.1E-3, 1.5E-3, 1.9E-3, 2.3E-3,
     &2.4E-3, 2.5E-3, 2.4E-3, 2.2E-3, 1.8E-3, 1.3E-3, 1.1E-3,
     &8.6E-4, 6.8E-4, 5.0E-4, 4.0E-4, 2.9E-4, 2.0E-4, 1.8E-4,
     &1.2E-4, 5.0E-5, 2.0E-5, 2.1E-5, 2.3E-5, 2.7E-5, 3.1E-5,
     &3.5E-5, 5.0E-5, 6.2E-5, 9.0E-5, 1.2E-4, 1.7E-4, 2.0E-4,
     &2.4E-4, 3.3E-4, 4.5E-4, 6.2E-4, 9.0E-4, 1.3E-3, 2.0E-3,
     &2.9E-3, 4.4E-3, 6.6E-3, 1.0E-2, 1.3E-2, 1.9E-2, 3.0E-2,
     &7.0E-2, 1.6E-1, 7.0E-2, 6.0E-3, 6.0E-3, 6.0E-3, 6.0E-3/
      DO 100 I=2,55
      IF (W .GT. WAVENO(I)) GO TO 100
      FACTOR= (WAVENO(I) - W )/(WAVENO(I) - WAVENO(I-1))
      REFLIQ=XIMG(I) + FACTOR*(XIMG(I-1) - XIMG(I))
      RETURN
 100  CONTINUE
      REFLIQ=XIMG(55)
      RETURN
      END
C----------------------------------------------------------------------- 
C END REFLIQ

C----------------------------------------------------------------------- 
C HAZES (Not currently used anywhere)

      SUBROUTINE HAZES(PR,DN,DM,Z,ZM,DD,DDM,DS,DSM,NRMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'prog_params'

	DIMENSION PR(100),DN(100),DM(100),Z(100),ZM(100)
	DIMENSION DDM(100),DD(100),DS(100),DSM(100)
        COMMON /ATM/ Z0(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
        
C For each level:
C	PMT = P model top
C	PMB = P model bot
C   	PRT = P Rages top
C	PRB = P Rages bot
C	M counts model layers, N counts haze layers

	DTOT = 0.
	PRMAX = PR(NRMAX)
	M = 1
	NTOP = 1
	NBOT = 2
	PMT = PRESS(M)
	PMB = PRESS(M+1)
        PRT = PR(NTOP)
	PRB = PR(NBOT)

C We have haze as long as we are at pressures less than the
C bottom of the Rages haze layer.

	DO 1000 WHILE (PMT.LT.PRMAX) 
         ZM(M)=0.0
	  DO 100 WHILE (PMT.GT.PR(NTOP+1))
	    NTOP = NTOP + 1
            PRT = PR(NTOP)
100	  CONTINUE
	  DO 200 WHILE (PMB.GT.PRB)
	    NBOT = NBOT + 1
            PRB = PR(NBOT)
200	  CONTINUE
          ZTOP = (Z(NTOP)-Z(NTOP+1))
          ZBOT = (Z(NBOT-1)-Z(NBOT))
	  TN = 0.
C If NBOT - NTOP = 1 then we have a special case of just part
C of 1 Rages layer
	  IF (NBOT-NTOP.EQ.1) THEN
            DZ=(LOG(PRESS(M)/PRESS(M+1))/LOG(PR(NTOP)/PR(NTOP+1)))*ZTOP
	    TN = DN(NTOP) * (DZ/ZTOP)
	    IF (TN.NE.0.0) ZM(M) = DZ
	    DDM(M) = DD(NTOP) 
	    DSM(M) = DS(NTOP) 
C	    PRINT *,'1 LAYER'
C  	    PRINT *,'DN(NTOP), DZ, TN',DN(NTOP),DZ,TN
	  ELSE
C First account for the fractional upper and lower layers
            DZ1 = (LOG(PRESS(M)/PR(NTOP+1))/LOG(PR(NTOP)/PR(NTOP+1))) 
     1              * ZTOP
            DZ2 = (LOG(PRESS(M+1)/PR(NBOT))/LOG(PR(NBOT-1)/PR(NBOT))) 
     1              * ZBOT
            DZ3 = ZBOT - DZ2
	    TN = DN(NTOP) * (DZ1/ZTOP)
            DMT = DD(NTOP) * TN
            DST = DS(NTOP) * TN
	    IF (TN.NE.0.0) ZM(M) = DZ1
C	    print *,'DT TOP',TN
            DT = DN(NBOT-1) * (DZ3/ZBOT)
	    IF (DT.NE.0.0) ZM(M) = ZM(M) + DZ3
C	    PRINT *,'DT BOT',DT
            TN = TN + DT
	    DMT = DMT + DD(NBOT-1) * DT
	    DST = DST + DS(NBOT-1) * DT
C	    PRINT *,'DZ1,DZ2,DZ3',DZ1,DZ2,DZ3
C Now account for the remainder intermediate layers
	    DO 300 I=1,NBOT-NTOP-2
	      DT = DN(NTOP+I)
	      TN = TN + DT
C	      PRINT *,'DT MID',DT
	      DMT = DMT + DD(NTOP+I)*DT
	      DST = DST + DS(NTOP+I)*DT
	      IF (DT.NE.0.0) ZM(M) = ZM(M) + (Z(NTOP+I)-Z(NTOP+I+1))
300	    CONTINUE 
	    DDM(M) = DMT/TN
	    DSM(M) = DST/TN
  	  ENDIF
	DM(M) = TN
	DTOT = DTOT + TN
	M = M + 1
        PMT = PRESS(M)
	PMB = PRESS(M+1)
	PRT = PR(NTOP)
	PRB = PR(NBOT)
1000	CONTINUE

	DTROT = 0
	DO 400 I =1,NRMAX+1
	 DTROT = DTROT +  DN(I)
400	CONTINUE
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	write (*,909) DTROT,DTOT
909	FORMAT(///,'HAZE:',/,' Rages Ntot = ',E12.5,3X
     1       ,'Model Ntot =',E12.5,' (These numbers should agree.)')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	 
	RETURN
	END
C----------------------------------------------------------------------- 
C END HAZES

C----------------------------------------------------------------------- 
C DTRIDGL (Called in DSOLVER)
C
C THIS SUBROUTINE SOLVES A SYSTEM OF TRIDIAGIONAL MATRIX
C EQUATIONS. THE FORM OF THE EQUATIONS ARE:
C A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)
C WHERE I=1,L  LESS THAN 103.
C (REVIEWED -CP)

      SUBROUTINE DTRIDGL(L,AF,BF,CF,DF,XK)
C DOUBLE PRESCISION VERSION OF TRIDGL
      PARAMETER (NMAX=301)
      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION AF(L),BF(L),CF(L),DF(L),XK(L)
      DIMENSION AS(NMAX),DS(NMAX)
      AS(L) = AF(L)/BF(L)
      DS(L) = DF(L)/BF(L)
      DO 10 I=2,L
           X=1./(BF(L+1-I) - CF(L+1-I)*AS(L+2-I))
           AS(L+1-I)=AF(L+1-I)*X
           DS(L+1-I)=(DF(L+1-I)-CF(L+1-I)*DS(L+2-I))*X
   10 CONTINUE
      XK(1)=DS(1)
      DO 20 I=2,L
           XKB=XK(I-1)
           XK(I)=DS(I)-AS(I)*XKB
   20 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------- 
C END DTRIDGL

C----------------------------------------------------------------------- 
C CLD (calls subroutine CLOUDSCAT, called in GEOMETRIC)

      SUBROUTINE CLD(IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C PUT IN A METHANE CLOUD HERE
C THIS ROUTINE SETS UP THE CLOUD DISTRIBUTION     
   
      include 'prog_params'

      COMMON /ATM/ Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
      COMMON /CLOUDV/ QEXTC,QSCTC(NSPECV),CBARC(NSPECV)
C      COMMON /GASS/ CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),CO(NLEVEL)
C     &   ,HE(NLEVEL),XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER)
C     &   ,H2O(NLEVEL),XNH3(NLEVEL),H(NLEVEL),Hmin(NLEVEL),elec(NLEVEL)
      COMMON /GASS/ CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),CO(NLEVEL),
     &  CO2(NLEVEL),HE(NLEVEL),XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER),
     &  H2O(NLEVEL),XNH3(NLEVEL),H(NLEVEL),Hmin(NLEVEL),elec(NLEVEL), 
     &  xNA(NLEVEL), xK(NLEVEL),H2S(NLEVEL),TIO(NLEVEL), VO(NLEVEL),
     &  RB(NLEVEL), CS(NLEVEL), FEH(NLEVEL), CRH(NLEVEL) 
      COMMON /CLOUD/ RADCLD(NLAYER), XNCLD(NLAYER),XSCTC
     & , RCLDI(NSPECI), XICLDI(NSPECI), RCLDV(NSPECV), XICLDV(NSPECV)
      COMMON /PLANT/ CSUBP,RSFI,RSFV,F0PI
      COMMON /ADJUST/ RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,RCLOUD
      COMMON /CONST/RGAS,RHOP,PI,SIGSB,RTOD
      DATA XMN2,XMH2,XMHE,XMCH4/28.0134d0,2.0158d0,4.0026d0,16.0426d0/
      TOTALC=0.0
              
C Critical humidity for cloud to be present:
      XC=.75
      XNC=0.85
      J0 = 0
      DO 190 J=NLAYER,1,-1
        XNCLD(J)=0.
        RADCLD(J)=0.
C J + 1 is level at bottom of layer J
        IF ( CH4(J+1)*PRESS(J+1)/PCH4(TEMP(J+1)) .GT. XC) THEN
C There may be a cloud in this layer, if this is the first cloud layer
C save the level at the base of the cloud.
          IF (J0.EQ.0) THEN
            J0 = J+1
C Save the saturation pressure ratio at this level
            PS0 = PCH4(TEMP(J+1))/PRESS(J+1)
          ENDIF
C We do not put a cloud in if at the base the following 
C condition is not true.  There is always 1 cloud layer.
          IF ((PCH4(TEMP(J+1))/PRESS(J+1))/PS0.GT.XNC) THEN
C We do have a cloud
            RADCLD(J)=RCLOUD
C We use Kathy Rages' well mixed cloud.  Cloud number
C density distribution q = 1 (see notes).     
C We dont worry about the units since we scale the cloud to
C the preset value of TAUFAC in the next loop. 
           XNCLD(J) = (-PRESS(J) + PRESS(J+1))/PRESS(J0) 
C          IF (IPRINT .GT. 1 ) WRITE(6,95) J,XNCLD(J)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           WRITE(6,95) J,XNCLD(J)
95         FORMAT(' CLOUD INSERTED: ',I3,1PE10.3)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           TOTALC=TOTALC+XNCLD(J)
          ELSE
C We are above the cloud top and there is no more cloud.
           RADCLD(J) = 0.0
           XNCLD(J) = 0.0
          ENDIF
        ENDIF
190    CONTINUE

C Old version
C      CTAU=QEXTC*TOTALC*XSCTC
      ctau=0.0

      IF (CTAU.NE.0.0) THEN
C SCALE THE CLOUD DENSITIES TO THE FIXED TAU 
        TOTALC = 0.0
        DO 145 J=1,NLAYER 
         XNCLD(J)=XNCLD(J)*TAUFAC/CTAU 
         TOTALC=TOTALC+XNCLD(J)
145     CONTINUE 
        CTAU=QEXTC*TOTALC*XSCTC
      ELSE
        DO 146 J=1,NLAYER
         XNCLD(J)=0.0
146     CONTINUE
      ENDIF 
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (IPRINT .GT. 0) WRITE(6,98) TOTALC,CTAU
98    FORMAT(' CLOUD COLUMN DENSITY , OPTICAL DEPTH= ' 2E10.2)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C Old version
C      print *,'Input CH4 supersaturation (1)'
C      accept *,csat
      csat=1.d0

      DO 134 J=NLEVEL-1,1,-1
        CH4SAT=PCH4(TEMP(J))/PRESS(J)
        CH4(J)=DMIN1(csat*CH4SAT,CH4(NLEVEL),CH4(J+1))
		
		fstrat = 1d-4
		Proll = 1.d-3
		Pmeso = 1d-6
		fmeso = 1.d-7
	    xm = dlog10(fstrat/fmeso) / dlog10(Proll/Pmeso)
		b = dlog10(fstrat) - xm * dlog10(Proll)
		ch4(j) = dmax1(ch4(j),fstrat)
C Atmosphere above about 1 mbar and below 10-5
C We fit a line in log f/ log P space
	  if (press(j).lt.proll) then
            ch4(j) = 10.d0**(xm*dlog10(press(j)) + b) 
	  endif
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	print *,j,ch4(j)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
134   CONTINUE
      DO 135 J=1,NLEVEL
        XMU(J)=XMN2*XN2(J)+XMCH4*CH4(J)+XMHE*HE(J)+XMH2*H2(J)
135   CONTINUE
      RETURN
      END
C----------------------------------------------------------------------- 
C END CLD	  
	  
C----------------------------------------------------------------------- 
C GET_K (called in OPTCV in subs_egp.f)

      SUBROUTINE get_K (t,p) ! p in mbars : program will take log
                             ! returns xreturn - loop over frequencies

C Include files
      include 'freedman/int.params'
      include 'freedman/abs_data.com'
      include 'freedman/limits.com'

C Parameters
      integer pset,tset,pset_u ! for interpolations
      real Kdum_up,Kdum_low    ! for interpolations
      data tset,pset,pset_u /3*0/     ! first guess for hunt:

C Arrays for interpolated [fine pressure] grid values - save for later

C For linear interpolations in P
      rlint (x1,x2,y1,y2,x) = y1 + (y2 - y1)/(x2 - x1) * (x - x1)

C For interpolations in 1./T - use with T interpolations
      rlint2 (x1,x2,y1,y2,x) = y1 + (y2 - y1)/(1./x2 - 1./x1) * 
     +                         (1./x - 1./x1)

C Take log of P
      plog = log(P)

C Use T to find layer limits - use corrected hunt - replace locate - start with guess of "0"
      call hunt (coarse_T,NT,t,tset)
C Perform checks, can use commented write statements for output	  
c      write (*,*) 't:p tset ',t,p,tset
      tset=min(tset,NT-1)
      tset=max(1,tset)
c      write (*,*) 'After Ck:t:p tset ',t,p,tset,coarse_T(tset),coarse_T(tset+1)

C Find pset at each T and interpolate in P
      NN = upper(tset)
      call hunt(log_CP,NN,plog,pset) ! return pset
C Perform checks, can use commented write statements for output	  
c      write (*,*) 'A t:p:pset ',t,p,pset
      pset=min(pset,Upper(tset)-1)
      pset=max(pset,1)
c      write (*,*) 'A:ck t:p:pset ',t,p,pset,coarse_P(pset),coarse_P(pset+1)


C Now do upper temperature check for same set of Pressures
      if (upper(tset).ne.upper(tset+1)) then
         NN = upper(tset+1)
         call hunt(log_CP,NN,plog,pset_u)
c        write (*,*) 'B t:p:pset ',t,p,pset_u
         pset_u=min(pset_u,Upper(tset+1)-1)
         pset_u=max(pset_u,1)
c        write (*,*) 'B:ck t:p:pset ',t,p,pset_u
      else
         pset_u = pset
      end if

C Now interpolate log K in P at each T - remember to loop over Nu points
      do jf = 1,NF ! loop over Nu
         Kdum_up = rlint(log_CP(pset_u),log_CP(pset_u+1),xabs(pset_u,tset+1,jf),
     +                   xabs(pset_u+1,tset+1,jf),plog)
         Kdum_low = rlint(log_CP(pset),log_CP(pset+1),xabs(pset,tset,jf),
     +                   xabs(pset+1,tset,jf),plog)

C Now interpolate in T as 1/T - convert back to real number
         xreturn(jf) = exp( rlint2(coarse_T(tset),coarse_T(tset+1),Kdum_low,Kdum_up,t) )

      end do ! end loop over Nu points
      return
      end
C----------------------------------------------------------------------- 
C END GET_K	  
	  
C----------------------------------------------------------------------- 
C HUNT (called in GET_K)

      SUBROUTINE hunt(xx,n,x,jlo) ! corrected version
      INTEGER jlo,n
      REAL x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
C----------------------------------------------------------------------- 
C END HUNT
	  
C----------------------------------------------------------------------- 
C BLOCK DATA FREEDMAN INTERPOLATION
C FOR ABSORPTION COEFFS
C USED IN OPTCV and GET_K
      block data

C JJF: This is needed for interpolation within the Freedman pretabulated
C x-sections.  It includes the paramaters for the # of temperatures at
C each pressure.  Pressure ranges from 1d-3 to 3e5 mbar and temperature
C ranges from 75 to 4000 K.  Not every pressure has every temperature.
      include 'freedman/int.params'
      include 'freedman/limits.com'
      include 'freedman/abs_data.com'

C Limits on coarse pressures at each T level
      data lower /42*1/

      data upper / 15,16, 
     +     17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,
     +     18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
     +     18,18,18,18,18,18,18,18,18,18/ ! 42 T levels

      data coarse_P/ 1.E-03, 3.E-03, 1.E-02, 3.E-02, 1.E-01, 3.E-01, 1.E+00,
     +               3.E+00, 1.E+01, 3.E+01, 1.E+02, 3.E+02, 1.E+03, 3.E+03,
     +               1.E+04, 3.E+04, 1.E+05, 3.E+05 /

      data coarse_T/ 75.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0,
     + 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0,
     + 250.0, 275.0, 300.0, 400.0, 500.0, 575.0, 650.0, 725.0,
     + 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0,
     + 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2300.0, 2600.0, 3000.0,
     + 3500.0, 4000.0 /

      end
C----------------------------------------------------------------------- 
C END BLOCK DATA FREEDMAN INTERPOLATION
	  
C----------------------------------------------------------------------- 
C BLOCK DATA GAUSSX (used in POINTS)

       BLOCK DATA GAUSSX
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       COMMON /GAUSS/X(46)

       DATA X/             0.50000000000000,  1.0000000000000,
     *  0.21132486540518,  0.50000000000000,
     *  0.11270166537925,  0.27777777777777,  0.50000000000000,
     *  0.44444444444444,
     *  0.06943184420297,  0.17392742256872,  0.33000947820757,
     *  0.32607257743127,
     *  0.04691007703066,  0.11846344252809,  0.23076534494715,
     *  0.23931433524968,  0.50000000000000,  0.28444444444444,
     *  0.03376524289842,  0.08566224618958,  0.16939530676686,
     *  0.18038078652406,  0.38069040695840,  0.23395696728634,
     *  0.02544604382862,  0.06474248308443,  0.12923440720030,
     *  0.13985269574463,  0.29707742431130,  0.19091502525255,
     *  0.50000000000000,  0.20897959183673,
     *  0.01985507175123,  0.05061426814518,  0.10166676129318,
     *  0.11119051722668,  0.23723379504183,  0.15685332293894,
     *  0.40828267875217,  0.18134189168918,
     *  0.01591988024618,  0.04063719418078,  0.08198444633668,
     *  0.09032408034742,  0.19331428364970,  0.13030534820146/

      END
C----------------------------------------------------------------------- 
C END BLOCK DATA GAUSSX
	  
C----------------------------------------------------------------------- 
C CLOUD2SCAT (called in GEOM when used)

      SUBROUTINE CLOUD2SCAT 
      implicit double precision (a-h,o-z)

      include 'prog_params'
      
      COMMON /ADJUST/ RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,RCLOUD
      COMMON /CLOUD2V/ TAU2(NLAYER),WBAR2(NLAYER),CBAR2(NLAYER)
      COMMON /SPECTV/ WNOV(NSPECV), WLNV(NSPECV)
      COMMON /ATM/ Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
      COMMON /CONST/RGAS,RHOP,PI,SIGSB,RTOD
      COMMON /CLOUDD/ XNCLD2(NLAYER)
       
C This routine sets up the lower H2S cloud.  For now this will be
C permanent, fixed cloud.  The optical properties come from Hammel

      write(*,*)'Set-up fixed lower H2S cloud'

C Lower cloud optical depth (3.0)
       TAU2FAC = 3.0   
       TOT = 0.0

C Lower cloud albedo (W2bar, 0.9985)
 	   ww=0.9985

      DO 100 J=1,NLAYER
       IF ((PRESS(J).GE.2.5).AND.(PRESS(J+1).LE.10.0)) THEN
        tau2(j)=0.0
C Old version
C Make Tau2 homogeneous, proportional to the layer delta_P
C	TAU2(J) = PRESS(J+1)-PRESS(J)
        TOT = TOT + TAU2(J)
        WBAR2(J) = ww
        CBAR2(J) = 0.00
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WRITE(*,99) J
99      FORMAT (1X,'H2S CLOUD INSERTED, LAYER ',I4,'.') 
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       ELSE
        TAU2(J) = 0.0
        WBAR2(J) = 0.0
        CBAR2(J) = 0.0
       ENDIF
100   CONTINUE

C We need an approximate estimate of total column abundance
C of cloud particles.  We assume geometric optics and that
C the lower cloud particle size is the same as the upper
C cloud (RCLOUD).

      NTOT = TAU2FAC / ( 2 * PI * (RCLOUD**2) * 1E-8)
      DO 150 J=1,NLAYER
       XNCLD2(J) = NTOT * TAU2(J)/TOT
	   tau2(j)=0.0
C Old version		 
C       TAU2(J) = TAU2(J) * TAU2FAC/TOT

150   CONTINUE

C   RCLOUD is cloud particle size in microns; we assume it equals
C   CH4 cloud particle size. NTOT is the column abundance
C   of cloud particles.

      RETURN
      END
C----------------------------------------------------------------------- 
C END CLOUD2SCAT
	  
C----------------------------------------------------------------------- 
C SCAMIE (when used, called in OPTCV)

      SUBROUTINE SCAMIE(     LAMB,    NREAL,  NIMAG, RMIN,
     &                        RM,     RMAX,   B, NF, PO,
     &                        Epsq,   Xfunny,
     &                        QEXBAR, QSCBAR, ALBED, BETA)
      implicit double precision (a-h,o-z)

C From an older function definition	  
C    &                        FSH,    SREAL,  SIMAG)


C THIS IS THE MIE SCATTERING SUBROUTINE
C LAMB: WAVELENGTH IN MICRONS
C TO CONVERT MICRONS TO CMS  MULTIPLY BY 10.E-4
C NREAL: REAL INDEX
C NIMAG: IMAGINARY INDEX


      REAL*8        NORMLP,   NREAL,    NIMAG,     LAMB
      REAL*8       Dqmax,    Dq,       Qlast,     FSH
C Old version
C     REAL         SREAL,    SIMAG

      COMPLEX*16   ACAP( 250000 )
C Old version	  
C     COMPLEX   ACAP( 250000 )
C     real*8 alfds,rfrs,rfis,tts,qscatds,qextds,qrprds

      DIMENSION    TT(1),    PIE(3,1), ELTRMX(4,1,2),
     "             TAU(3,1), CSTHT(1), SI2THT(1),
     "             RF(100),   QMI(30),  QMUWT(30)
C Old version
C    "             ,tts(1)

C JOX=0 MEANS NO PRINTOUT       JOX=1 GIVES A PRINTOUT
      DATA      IXMAX  /     50          /
      DATA      M      /      8          /
      DATA      JOX    /      0          /

            TT(1) = 0.0
            JX = 1
            IT = 1
       LL   = 250000
            PI = 3.1415926535897

C SETUP RI & RF VALUES FOR X SPACE INTEGRATION
      XMEAN=2.0*PI*RM/LAMB
      XMIN = XMEAN * RMIN
      XMAX = XMEAN * RMAX

      CALL MKBNDS(XMIN,XMAX,Xfunny,NREAL, NIMAG, IX,  RF)
 162  FORMAT(2X,5(F10.4,2X))
C     write(6,162) (RF(i), i=1,IX)
      RI = XMIN
      NPART = IX
      RFI =  NIMAG
      RFR =  NREAL
      WVNO = 2.*PI/LAMB

      IF (JOX .GE. 1) WRITE(6,165) RI,NPART,(RF(I),I=1,NPART)
  165 FORMAT ('0','RI = ',F10.7,'NPART = ',I10,/,8(1X,F12.7))

      NORMLP = 0.
      QSUMXX  = 0.
      QSUMSS  = 0.
      BETA = 0.0

C CALCULATE MIE SCATTERING
         Dqmax = 0.0
         DO 2000 KI = 1,NPART
         AI = RI
         AF = RF(KI)
         MM = M
         CALL POINTS( AF, AI, MM, QMI, QMUWT)
            Qlast = QSUMXX
            DO 2100 K = 1,MM
            ALF = QMI(K)
            ALFD = QMI(K)
            A = ALF/WVNO
C HERE DAMIE IS CALLED AND WE DESIRE NO PHASE FUNCTION.
C THUS TT(1)=0 AND IT=1 AND JX=1.
C           alfds = sngl(alfd)
C           rfrs = sngl(rfr)
C           rfis = sngl(rfi)
C           tts(1)=sngl(tt(1))
C           CALL DAMIE(  ALFDs,     RFRs,      RFIs,      TTs,
C    "                   JX,       QEXTDs,    QSCATDs,   QRPRDs,
C    "                   ELTRMX,   PIE,      TAU,      CSTHT,
C    "                   SI2THT,   ACAP,     IT,       LL     )
C           qextd=dble(qextds)
C           qscatd = dble(ascatds)
C           qrprd = dble(qrprds)
C Old version 
C            CALL DAMIE(  ALFD,     RFR,      RFI,      TT,
C     "                   JX,       QEXTD,    QSCATD,   QRPRD,
C     "                   ELTRMX,   PIE,      TAU,      CSTHT,
C     "                   SI2THT,   ACAP,     IT,       LL     )
C Old version
C    Here we have substituted DMIESS() for DAMIE() in order to
C    calculate scattering parameters for particles composed of
C    a core material, and a shell, or mantle, of different material.

C     RO  =  ALFD/WVNO
C     RC  =  RO * (1 - FSH)**0.333
C           CALL DMIESS(  RO,     SREAL,     SIMAG,    TT,      JX,
C    &                    QEXTD,  QSCATD,  QRPRD,  ELTRMX,  PIE,
C    &                    TAU,    CSTHT,   SI2THT, ACAP,    QBS, IT,
C    &                    LL,     RC,      RFR,  RFI,   WVNO    )

            QEXT  = QEXTD
            QSCAT = QSCATD

            DIST=FUNC(A,RM,B,PO,NF)

            XXXX = DIST*PI*A*A
            YYYY = QMUWT(K)/WVNO

            QSXX = QSCAT*XXXX*YYYY
            QXXX = QEXT*XXXX*YYYY
            QSUMXX = QSUMXX + QXXX
            QSUMSS = QSUMSS + QSXX
            NORMLP = NORMLP + YYYY*XXXX
            BETA = BETA + QRPRD * XXXX * YYYY

 2100       CONTINUE

       RI = RF(KI)
       Dq = QSUMXX - Qlast
       IF (Dq.LT.Epsq*Dqmax) GOTO 2059
       Dqmax = MAX(Dq, Dqmax)

 2000 CONTINUE


 2059 Radis = RF(KI)*LAMB/(2.0*PI)
      QSCBAR = QSUMSS/NORMLP
C Old version
C     QSCBAR = QSUMSS

      QEXBAR = QSUMXX/NORMLP
C Old version
C     QEXBAR = QSUMXX
      BETA = BETA / QSUMSS
C Item 1, Enclosure A
C          QEXBAR = QEXBAR - QSCBAR * BETA
C          QSCBAR = QSCBAR * ( 1.0 - BETA )

      ALBED = QSCBAR/QEXBAR

      RETURN
      END
C----------------------------------------------------------------------- 
C END SCAMIE
	  
C----------------------------------------------------------------------- 
C DAMIE (called in SCAMIE)

      SUBROUTINE DAMIE(      X,        RFR,      RFI,      THETD,
     &                       JX,       QEXT,     QSCAT,    CTBRQS,
     &                       ELTRMX,   PI,       TAU,      CSTHT,
     &                       SI2THT,   ACAP,     IT,       LL     )
      implicit double precision (a-h,o-z)

      DIMENSION    T(5),     TA(4),    TB(2),
     &             TC(2),    TD(2),    TE(2),
     &             ELTRMX(4,IT,2),     PI(3,IT),
     &             TAU(3,IT),          CSTHT(IT),
     &             SI2THT(IT),         THETD(IT)

      COMPLEX*16    RF,       RRF,      RRFX,     WM1,
     &             FNA,      FNB,      TC1,      TC2,
     &             WFN(2),   FNAP,     FNBP,     ACAP(250000)

       EQUIVALENCE  (FNA,TB(1)),     (FNB,TC(1)),
     &             (FNAP,TD(1)),       (FNBP,TE(1))

C LARGE VERSION
C
C
C      SUBROUTINE FOR COMPUTING THE PARAMETERS OF THE ELECTROMAGNETIC
C      RADIATION SCATTERED BY A SPHERE.
C      THIS SUBROUTINE COMPUTES THE CAPITAL A FUNCTION BY MAKING USE OF
C      THE DOWNWARD RECURRENCE RELATIONSHIP.
C  X: SIZE PARAMETER OF THE SPHERE,( 2 * PI * RADIUS OF THE SPHERE
C      DIVIDED BY THE WAVELENGTH OF THE INCIDENT RADIATION).
C  RF:  REFRACTIVE INDEX OF THE MATERIAL OF THE SPHERE.
C      COMPLEX QUANTITY..FORM: (RFR - I * RFI )
C  THETD(J): ANGLE IN DEGREES BETWEEN THE DIRECTIONS OF THE
C      INCIDENT AND THE SCATTERED RADIATION.  THETD(J) IS< OR= 90.0
C      IF THETD(J) SHOULD HAPPEN TO BE GREATER THAN 90.0, ENTER WITH
C      THE SUPPLEMENTARY VALUE; SEE COMMENTS BELOW ON ELTRMX.
C  JX:  TOTAL NUMBER OF THETD FOR WHICH THE COMPUTATIONS ARE
C      REQUIRED.  JX SHOULD NOT EXCEED **IT** UNLESS THE DIMENSIONS
C      STATEMENTS ARE APPROPRIATEDLY MODIFIED.
C      MAIN PROGRAM SHOULD ALSO HAVE REAL*8 THETD(IT ),ELTRMX(4,IT ,2)
C      THE DEFINITIONS FOR THE FOLLOWING SYMBOLS CAN BE FOUND IN
C         #LIGHT SCATTERING BY SMALL PARTICLES, H.C.VAN DE HULST,
C          JOHN WILEY & SONS, INC., NEW YORK, 1957# .
C    QEXT:  EFFIECIENCY FACTOR FOR EXTINCTION,VAN DE HULST,P.14 & 127.
C    QSCAT: EFFIECINCY FACTOR FOR SCATTERING,V.D. HULST,P.14 & 127.
C    CTBRQS:  AVERAGE(COSINE THETA) * QSCAT,VAN DE HULST,P.128.
C    ELTRMX(I,J,K):  ELEMENTS OF THE TRANSFORMATION MATRIX F,
C         V.D.HULST, P.34,45 & 125.   I=1: ELEMENT M SUB 2..
C         I=2: ELEMENT M SUB 1..      I=3: ELEMENT S SUB 21..
C         I=4: ELEMENT D SUB 21..
C    ELTRMX(I,J,1) REPRESENTS THE ITH ELEMENT OF THE MATRIX FOR
C    THE ANGLE THETD(J).. ELTRMX(I,J,2) REPRESENTS THE ITH ELEMENT
C    OF THE MATRIX FOR THE ANGLE 180.0 - THETD(J) ..
C
C
C      IN THE ORIGINAL PROGRAM THE DIMENSION OF ACAP WAS 7000.
C      FOR CONSERVING SPACE THIS SHOULD BE NOT MUCH HIGHER THAN
C      THE VALUE, N=1.1*(NREAL**2 + NIMAG**2)**.5 * X + 1
C      TA(1): REAL PART OF WFN(1).. TA(2): IMAGINARY PART OF WFN(1)..
C      TA(3): REAL PART OF WFN(2).. TA(4): IMAGINARY PART OF WFN(2).
C      TB(1): REAL PART OF FNA..    TB(2): IMAGINARTY PART OF FNA..
C      TC(1): REAL PART OF FNB..    TC(2): IMAGINARY PART OF FNB..
C      TD(1): REAL PART OF FNAP..   TD(2): IMAGINARY PART OF FNAP.
C      TE(1): REAL PART OF FNBP..   TE(2): IMAGINARY PART OF FNBP.
C      FNAP & FNBP ARE THE PRECEDING VALUES OF FNA & FNB RESPECTIVELY.
C
C
C
      IF(JX.LE. IT)   GO TO 20
      WRITE(6,6)
    6 FORMAT(//'PLEASE READ COMMENTS.'///)
      WRITE(6,7)
    7 FORMAT(//10X,'THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT')
      CALL EXIT
   20 RF =  DCMPLX(RFR,-RFI)
      RRF = 1.0d0/RF
      RX = 1.0d0/X
      RRFX = RRF * RX
      T(1) = (X**2)*(RFR**2 + RFI**2)
      T(1) =  DSQRT(T(1))
      NMX1 = 1.10d0 * T(1)
C     WRITE(6,199) X,RFR,RFI,LL,NMX1
C 199 FORMAT(1H ,*X=*,E14.6,2X,*RFR=*,E14.6,2X,*RFI=*,E14.6,2X,*LL=*,I10
C    1,2X,*NMX1?*,I10)
       IF( NMX1 .LE. LL-1 ) GO TO 21
      WRITE(6,199) X,RFR,RFI,LL,NMX1
  199 FORMAT(' X=',E14.6,'  RFR=',E14.6,'  RFI=',E14.6,'  LL=',I10, '  NMX1=',I10)
       WRITE(6,8)
    8 FORMAT(///10X,'THE UPPER LIMIT OF ACAP IS NOT ENOUGH. SUGGEST YOU GET DETAILED OUTPUT TO MODIFY THE SUBROUTINE.'////)
       STOP
   21  NMX2 = T(1)
C Values of nmx1 and nmx2; changed by Bruce Christofferson for
C longer iterations in calculating large size particle scattering.
       IF( NMX1 .GT. 2100 )  GOTO 22
           NMX1 = 2100
           NMX2 = 1890
   22  ACAP(NMX1 + 1 ) = (0.0d0,0.0d0)
             DO 23 N = 1,NMX1
             NN = NMX1 - N + 1
             ACAP(NN) = (NN+1) * RRFX - 1.0d0/((NN+1)*RRFX + ACAP(NN+1))
   23        CONTINUE
          DO 30 J = 1,JX
          IF ( THETD(J) .LT. 0.0d0 ) THETD(J) =  DABS(THETD(J))
          IF ( THETD(J) .GT. 0.0d0 ) GO TO 24
          CSTHT(J) = 1.0d0
          SI2THT(J) = 0.0d0
            GO TO 30
   24     IF ( THETD(J) .GE. 90.0d0) GO TO 25
          T(1) = ( 3.1415926535898 * THETD(J))/180.0d0
          CSTHT(J) =  DCOS(T(1))
          SI2THT(J) = 1.0d0 - CSTHT(J)**2
            GO TO 30
   25     IF ( THETD(J) .GT. 90.0d0) GO TO 28
          CSTHT(J) = 0.0d0
          SI2THT(J) = 1.0d0
            GO TO 30
   28     WRITE(6,6)
          WRITE(6,5) THETD(J)
    5 FORMAT(10X,'THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN 90
     "DEGREES.  IT IS',E20.10)
          CALL EXIT
   30     CONTINUE
                DO 35 J = 1,JX
                PI(1,J) = 0.0d0
                PI(2,J) = 1.0d0
                TAU(1,J) = 0.0d0
                TAU(2,J) = CSTHT(J)
   35           CONTINUE
      T(1) =  DCOS(X)
      T(2) =  DSIN(X)
      WM1 =  DCMPLX( T(1),-T(2))
      WFN(1) =  DCMPLX(T(2),T(1))
      TA(1) = T(2)
      TA(2) = T(1)
      WFN(2) = RX * WFN(1) - WM1
C     TA(3) = REAL( WFN(2) )
      TA(3) = DBLE( WFN(2) )
C      TA(4) = AIMAG( WFN(2) )
C      TA(4) = DIMAG( WFN(2) )
      TA(4) = IMAG( WFN(2) )
      TC1 = ACAP(1) * RRF + RX
      TC2 = ACAP(1) * RF + RX
      FNA = (TC1*TA(3) - TA(1))/(TC1*WFN(2) - WFN(1))
      FNB = ( TC2*TA(3) - TA(1))/(TC2  * WFN(2) - WFN(1))
      FNAP = FNA
      FNBP = FNB
      T(1) = 1.50d0

C      FROM HERE TO THE STATMENT NUMBER 90, ELTRMX(I,J,K) HAS
C      FOLLOWING MEANING:
C      ELTRMX(1,J,K): REAL PART OF THE FIRST COMPLEX AMPLITUDE.
C      ELTRMX(2,J,K): IMAGINARY PART OF THE FIRST COMPLEX AMPLITUDE.
C      ELTRMX(3,J,K): REAL PART OF THE SECOND COMPLEX AMPLITUDE.
C      ELTRMX(4,J,K): IMAGINARY PART OF THE SECOND COMPLEX AMPLITUDE.
C      K = 1 : FOR THETD(J) AND K = 2 : FOR 180.0 - THETD(J)
C      DEFINITION OF THE COMPLEX AMPLITUDE: VAN DE HULST,P.125.
       TB(1) = T(1) * TB(1)
       TB(2) = T(1) * TB(2)
       TC(1) = T(1) * TC(1)
       TC(2) = T(1) * TC(2)
           DO 60 J = 1,JX
           ELTRMX(1,J,1) = TB(1) * PI(2,J) + TC(1) * TAU(2,J)
           ELTRMX(2,J,1) = TB(2) * PI(2,J) + TC(2) * TAU(2,J)
           ELTRMX(3,J,1) = TC(1) * PI(2,J) + TB(1) * TAU(2,J)
           ELTRMX(4,J,1) = TC(2) * PI(2,J) + TB(2) * TAU(2,J)
           ELTRMX(1,J,2) = TB(1) * PI(2,J) - TC(1) * TAU(2,J)
           ELTRMX(2,J,2) = TB(2) * PI(2,J) - TC(2) * TAU(2,J)
           ELTRMX(3,J,2) = TC(1) * PI(2,J) - TB(1) * TAU(2,J)
           ELTRMX(4,J,2) = TC(2) * PI(2,J) - TB(2) * TAU(2,J)
   60      CONTINUE
      QEXT = 2.0d0 * ( TB(1) + TC(1))
C     WRITE(6,999) TB(1),TC(1)
      QSCAT =(TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2)/0.75d0
      CTBRQS = 0.0d0
      N = 2
   65 T(1) = DBLE(2*N - 1)
      T(2) = DBLE(N - 1)
      T(3) = DBLE(2 * N + 1)
        DO 70 J = 1,JX
        PI(3,J) = (T(1)*PI(2,J)*CSTHT(J)-DBLE(N)*PI(1,J))/T(2)
        TAU(3,J) = CSTHT(J)*(PI(3,J)-PI(1,J))-T(1)*SI2THT(J)*PI(2,J)+
     &  TAU(1,J)
   70   CONTINUE
      WM1 = WFN(1)
      WFN(1) = WFN(2)
      WFN(2) = T(1)*RX*WFN(1) - WM1
C      TA(1) = REAL( WFN(1) )
       TA(1) = DBLE( WFN(1) )
C       TA(2) = AIMAG( WFN(1) )
C       TA(2) = DIMAG( WFN(1) )
      TA(2) = IMAG( WFN(1) )
C      TA(3) = REAL( WFN(2) )
       TA(3) = DBLE( WFN(2) )
C      TA(4) = AIMAG( WFN(2) )
C        TA(4) = DIMAG( WFN(2) )
       TA(4) = IMAG( WFN(2) )
      TC1 = ACAP(N)*RRF + N*RX
      TC2 = ACAP(N)*RF + N*RX
      FNA = (TC1*TA(3)-TA(1))/(TC1*WFN(2) - WFN(1))
      FNB = (TC2*TA(3)-TA(1))/(TC2*WFN(2) - WFN(1))
      T(5) = DBLE(N)
      T(4) = T(1)/(T(5)*T(2))
      T(2) = (T(2)*(T(5) + 1.0d0))/T(5)
      CTBRQS = CTBRQS + T(2)*(TD(1)*TB(1) + TD(2)*TB(2) + TE(1)*TC(1)+
     &TE(2)*TC(2))  + T(4)*(TD(1)*TE(1) + TD(2)*TE(2))
      QEXT = QEXT + T(3)*(TB(1)+TC(1))
      T(4) = TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2
      QSCAT = QSCAT + T(3) *T(4)
      T(2) = DBLE(N*(N+1))
      T(1) = T(3)/T(2)
C     WRITE(6,999) T(1),QEXT
      K = (N/2)*2
       DO 80 J = 1,JX
       ELTRMX(1,J,1)=ELTRMX(1,J,1)+T(1)*(TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,1)=ELTRMX(2,J,1)+T(1)*(TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,1)=ELTRMX(3,J,1)+T(1)*(TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,1)=ELTRMX(4,J,1)+T(1)*(TC(2)*PI(3,J)+TB(2)*TAU(3,J))
           IF ( K .EQ. N) GO TO 75
       ELTRMX(1,J,2)=ELTRMX(1,J,2) +T(1)*(TB(1)*PI(3,J)-TC(1)*TAU(3,J))
       ELTRMX(2,J,2)=ELTRMX(2,J,2)+T(1)*(TB(2)*PI(3,J)-TC(2)*TAU(3,J))
       ELTRMX(3,J,2)=ELTRMX(3,J,2)+T(1)*(TC(1)*PI(3,J)-TB(1)*TAU(3,J))
       ELTRMX(4,J,2)=ELTRMX(4,J,2)+T(1)*(TC(2)*PI(3,J)-TB(2)*TAU(3,J))
          GO TO 80
   75  ELTRMX(1,J,2)=ELTRMX(1,J,2)+T(1)*(-TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,2)=ELTRMX(2,J,2)+T(1)*(-TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,2)=ELTRMX(3,J,2)+T(1)*(-TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,2)=ELTRMX(4,J,2)+T(1)*(-TC(2)*PI(3,J)+TB(2)*TAU(3,J))
   80  CONTINUE
           IF( T(4) .LT.  1.0E-14 ) GO TO 100
       N = N + 1
              DO 90 J = 1,JX
              PI(1,J) = PI(2,J)
              PI(2,J) = PI(3,J)
              TAU(1,J) = TAU(2,J)
              TAU(2,J) = TAU(3,J)
   90         CONTINUE
       FNASV = FNAP
       FNBSV = FNBP
       FNAP = FNA
       FNBP = FNB
         IF ( N .LE. NMX2 )    GO TO 65
       WRITE(6,8)
       WRITE( 6, 105 ) X, RFR, RFI, NMX1, NMX2,
     ;                TB(1), TB(2), TC(1), TC(2), T(4), FNASV,
     ;                FNA, FNBSV, FNB
105   FORMAT(/ 5X, " X = ", E14.6, 2X, "RFR=", E14.6,
     ;        2X, "RFI = ", E14.6, 2X, "MNX1 = ", I3, 2X,
     ;        "NMX2 = ", I3 // 5X, "TB(1) = ", E14.6, 2X,
     ;        "TB(2) = ", E14.6, 2X, " TC(1) = ", E14.6, 2X,
     ;       "TC(2) = ", E14.6 // 5X, "T(4) =", E14.6// 5X,
     ;       "FNASV = ", E14.6, 2X, " FNA = ", 2E14.6,// 5X,
     ;       "FNBSV = " , E14.6, 2X, " FNB = ", 2E14.6  /     )
       CALL EXIT
  100       DO 120 J = 1,JX
            DO 120 K = 1,2
                    DO  115  I= 1,4
                    T(I) = ELTRMX(I,J,K)
  115               CONTINUE
            ELTRMX(2,J,K) = T(1)**2 + T(2)**2
            ELTRMX(1,J,K) = T(3)**2 + T(4)**2
            ELTRMX(3,J,K) = T(1)*T(3) + T(2)*T(4)
            ELTRMX(4,J,K) = T(2)*T(3) - T(4)*T(1)
  120       CONTINUE
      T(1) = 2.0d0 * RX**2
      QEXT = QEXT * T(1)
C      WRITE(6,999) T(1),QEXT
C 999 FORMAT(1H ,3E20.8)
      QSCAT = QSCAT * T(1)
      CTBRQS = 2.d0 * CTBRQS * T(1)

C      THE DETAIL ABOUT THIS SUBROUTINE CAN BE FOUND IN THE FOLLOWING
C      REPORT: # SUBROUTINES FOR COMPUTING THE PARAMETERS OF THE
C      ELECTROMAGNETIC RADIATION SCATTERED BY A SPHERE # J.V. DAVE.
C      I B M SCIENTIFIC CENTER, PALO ALTO , CALIFORNIA.
C      REPORT NO. 320 - 3236 .. MAY 1968 ..

       RETURN
       END
C----------------------------------------------------------------------- 
C END DAMIE
	  
C----------------------------------------------------------------------- 
C MKBNDS (Called by SCAMIE)

      SUBROUTINE MKBNDS(xmin, xmax, Xfunny, n, k, ix, rf)
      implicit double precision (a-h,o-z)
      PARAMETER ( PI = 3.14159265358979323 )

C ICOUNT IS USED TO CALCULATE PLANK FUNCTION.
      PARAMETER (ICOUNT = 2000)
C Control dimensions of parameter arrays for MKBNDS intervals
      PARAMETER(IBNDS=6)
      PARAMETER(INTRVLS=5)
C Boolean values are useful for understanding conditionals in a program
      PARAMETER (IYES = 1)
      PARAMETER (NO  = 0)

C MKBNDS:  Make Bounds                                
C The Make_Bounds subroutine is called by the SCATMIE subroutine to break up a large 
C integration interval into (INTRVLS) smaller intervals; each of these is in turn 
C subdivided into a set of even more closely spaced intervals.        
C
C These spacings may be either linear or logarithmic, and are selected with the 
C intent of optimizing the accuracy of the resultant integral (performed in the 
C calling routine.)                               
C
C The MKBNDS routine returns an array (RF) containing all of the endpoints of all 
C of the intervals.
C                                            
C INPUT:                                              
C xmin, xmax           endpoints of original region  
C Xfunny               scaling factor                
C n, k                 real, imaginary refr. indices 
C
C INTERNAL DATA:                                     
C cb, clinstp, clogfac constants for determining     
C spacings in lin, log regions 
C
C OUTPUT:                                             
C IX                   total number of boundary pts.
C RF()                 the boundary points, from    
C                      RF(1)=xmin to RF(IX)=xmax    
C
C  AUTHOR:               DATE:                        
C  Brad Dalton           March, 1988                 


      REAL*8 xmin, xmax, Xfunny, n, k
      REAL*8 xpeak1, xpeak2, xpeak, Fr, rstart, root
C      REAL*8 xinc, xbdpA, xbdpB
      DOUBLE PRECISION xinc, xbdpA, xbdpB
      DIMENSION cb(IBNDS), clinstp(INTRVLS), clogfac(INTRVLS)
      DIMENSION xbd(IBNDS), step(INTRVLS), rf(1)
	  
      DATA cb      / 1.0d-04, 5.0d-02, 0.5, 5.0, 20.0, 1.0d+04 /
      DATA clinstp / 0.00, 0.1, 0.25, 1.0, 0.00 /
      DATA clogfac / 10.0, 0.0, 0.00, 0.0, 10.0 /

      xpeak1 = MAX(1.0,Xfunny/ABS(n-1.0))
      xpeak2 = 1.0
      Fr     = k/(ABS(n-1.0)+k)
      xpeak  = (1.0-Fr)*xpeak1 + Fr*xpeak2
      DO 111   i = 1,IBNDS
         xbd(i)  = cb(i)*xpeak
 111  CONTINUE
      DO 222   i = 1,INTRVLS
         step(i) = clinstp(i)*xpeak
 222  CONTINUE
 223  FORMAT("  xbd(i): ",6(F12.5,2X))
C     write(6,223) (xbd(i), i=1,IBNDS)

C Now we can calculate the actual limits of the intervals in X-spac

      ix = 0
      DO 666   i = 1,INTRVLS
         if (xmax.LT.xbd(i).OR.xmin.GT.xbd(i+1)) then
            nix     = 0
         else
            xbdpA   = DMAX1(xbd(i),xmin)
            xbdpB   = DMIN1(xbd(i+1),xmax)
            if (clinstp(i).GT.0.0) then
               nix = 1 + (xbdpB-xbdpA)/step(i)
            else
               nix = 1 + clogfac(i)*DLOG10(xbdpB/xbdpA)
            endif
            if (i.EQ.1.OR.ix.EQ.0) then
               rstart = xmin
               idelt = 0
            else
               rstart = rf(ix)
               idelt = ix
            endif
            if (clinstp(i).GT.0.0) then
               xinc = (xbdpB-xbdpA)/nix
               DO 333 j = 1,nix
                  jpr   = j + idelt
                  rf(jpr) = rstart + j*xinc
 333           continue
            else
               root = (xbdpB/xbdpA)**(1.0/nix)
               DO 444 j = 1,nix
                  jpr = j+idelt
                  rf(jpr) = rstart*root**j
 444           continue
            endif
         endif
         ix = ix + nix
 555  FORMAT(2X,"region: ",I2,2X,F10.3," to ",F10.3,
     *       " total ",I3,"values")
 556  FORMAT("  No values calculated for region ",I2)
C     if (nix.NE.0) then
C     write(6,555) i, xbdpA, xbdpB, nix
C     else
C     write(6,556) i
C     endif
 666  continue
      return
      end
C----------------------------------------------------------------------- 
C END MKBNDS
	  
C----------------------------------------------------------------------- 
C POINTS (called in SCAMIE)

       SUBROUTINE POINTS(PTI,PTF,IPTS,QMI,QMUWT)
      implicit double precision (a-h,o-z)

       DIMENSION QMI(1), QMUWT(1)

       COMMON  /GAUSS/ X(46)

       IF( MOD(IPTS,2) .EQ. 0 )    GO TO 100
       I = (IPTS-1)/2
       IPLACE = I*(I+1)*2
       JPTS =(IPTS+1)/2
               DO 1 J = 1,JPTS
               IIP = IPLACE + 2*J
               IIPM = IIP  -   1
               QMI(J) = X(IIPM)
    1          QMUWT(J) = X(IIP)
      GO TO 200

  100 CONTINUE
       I = IPTS/2
       IPLACE = I*(I+1)*2 - IPTS
       JPTS = IPTS/2
               DO 2 J = 1,JPTS
               IIP2 = IPLACE  +  2*J
               IIPM2 = IIP2 - 1
               QMI(J) = X(IIPM2)
    2          QMUWT(J) = X(IIP2)
  200 CONTINUE
       JPTS = (IPTS+1)/2
       IF( MOD(IPTS,2) .EQ. 0 )    JPTS = IPTS/2
               DO 210 I=1,JPTS
               IIPTS = IPTS - I +  1
               QMI(IIPTS) = QMI(I)
               QMI(I) = 1.d0-QMI(I)
  210          QMUWT(IIPTS) = QMUWT(I)
       SCALQ = PTF-PTI
       SCALE = DABS(PTF-PTI)
               DO 211 I = 1,IPTS
               QMUWT(I) = SCALE*QMUWT(I)
  211          QMI(I) = QMI(I)*SCALQ + PTI
       RETURN
       END

C----------------------------------------------------------------------- 
C END POINTS
	  
C----------------------------------------------------------------------- 
C SETPIAI (called in GEOM in egpalb.f)
C JJF: This was taken from the EGP code and dropped in.  The old SETPIAI
C has been removed.

      SUBROUTINE SETPIAI
      implicit double precision (a-h,o-z)
      include 'prog_params'
      
      PARAMETER (NTEMPS3=198,nint=1000)
      
      COMMON /PIAC3/ PIAHH3(NTEMPS3,NSPECV),PIAHE3(NTEMPS3,NSPECV),
     &  PIACH3(NTEMPS3,NSPECV),PIAH3(NTEMPS3,NSPECV)
      COMMON /H2MINUS/ H2M(NTEMPS3,NSPECV)
      COMMON /TEMPS3/ TSHH(NTEMPS3),TSHE(NTEMPS3)
      COMMON /SPECTI/ BWNI(NSPECV),WNOI(NSPECV),DWNI(NSPECV),
     &  WLNI(NSPECV)
      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      common /piaspline/ y2hh(ntemps3,NSPECV),y2he(ntemps3,NSPECV),
     & y2h(ntemps3,NSPECV),y2ch(ntemps3,NSPECV)
      dimension pass(ntemps3),y2x(ntemps3)
      dimension wnc(nint),ph2he(nint,NTEMPS3),ph2h2(nint,NTEMPS3)
      dimension ph2h(nint,NTEMPS3),ph2ch4(nint,NTEMPS3)
      double precision kappa2(NSPECV),kappa3(NSPECV)
  

C Read info from a file
C
      open (8,file='../inputs/H2CIA.DAT')
      read (8,*) idum1, idum2

      if (idum2.ne.ntemps3) then
        print *,'bad input file H2CIA.dat ntemps'
        print *, ntemps3, idum2
		stop
      endif
      if (idum1.ne.nint) then
        print *,'bad input file H2CIA.dat nint'
		stop
      endif
           wtf = 1
      DO kt=1,NTEMPS3
        wtf = wtf+1
	   read (8,*) TSHE(kt)
	   TSHH(kt) = TSHE(kt)
	    do KK = 1, nint
        wtf= wtf+1
		   read (8,*) wnc(kk),ph1,ph2,ph3,ph4 ,ph5

C Leave in logs until we interpolate
C Taking the exp (undo log) puts the units into cm^-1 amagat^-2
C ph2he(kk,kt) = 10.d0**ph2
C ph2h2(kk,kt) = 10.d0**ph1
 		   ph2h2(kk,kt) = ph1

 		  ph2he(kk,kt) = ph2
 		  ph2h(kk,kt)  = ph3
 		  ph2ch4(kk,kt)= ph4
c-nb add fake opacity to see if that is the problem
c      if (wnc(kk).ge.10000) then
c           print *, 'swap',tshe(kt),wnc(kk),ph2h2(kk,kt),ph2h2(500,kt)
c           ph2h2(kk,kt) = ph2h2(500,kt)
c           print *, 'swap',tshe(kt),wnc(kk),ph2h(kk,kt),ph2h(500,kt)
c           ph2h(kk,kt) = ph2h(500,kt)
c      endif
      
      enddo

      enddo
		print *,'T range for H2 opacity ',tshe(1),tshe(ntemps3)
		close(8)

		deltastep = wnc(2) - wnc(1)
		freqmax = wnc(nint)
		freqlo = wnc(1)

 	   do kt=1,ntemps3
		do k=1,NSPECV
	     if (WNOV(k+1).lt.freqmax) then
		  j = 1
		  do while (wnc(j).lt.WNOV(k)) 
			j = j+1
          enddo
		  jlo = j - 1 
		  do while (wnc(j).lt.WNOV(k+1)) 
			j = j+1
          enddo
		  jhi = j 

		   sumhh = 0.0
		   sumhe = 0.0
		   sumh  = 0.0
		   sumch = 0.0

		   count = 0.0

		  do jc = jlo+1, jhi-1

		   sumhh = sumhh + ph2h2(jc,kt)
		   sumhe = sumhe + ph2he(jc,kt)
		   sumh  = sumh  + ph2h(jc,kt)
		   sumch = sumch + ph2ch4(jc,kt)

		   count = count + 1.0
		  enddo

		   xlowhe = ph2he(jlo,kt) + ( (WNOV(k)-wnc(jlo))*
     1		      (- ph2he(jlo,kt) + ph2he(jlo+1,kt))/deltastep )
		   xhihe = ph2he(jhi-1,kt) + ( (WNOV(k+1)-wnc(jhi-1))*
     1		      (- ph2he(jhi-1,kt) + ph2he(jhi,kt))/deltastep )
		   xmeanhe = (sumhe + 0.5*(xlowhe + xhihe))/(count+1.0)	

		   xlowhh = ph2h2(jlo,kt) + ( (WNOV(k)-wnc(jlo))*
     1		      (- ph2h2(jlo,kt) + ph2h2(jlo+1,kt))/deltastep )
		   xhihh = ph2h2(jhi-1,kt) + ( (WNOV(k+1)-wnc(jhi-1))*
     1		      (- ph2h2(jhi-1,kt) + ph2h2(jhi,kt))/deltastep )
		   xmeanhh = (sumhh + 0.5*(xlowhh + xhihh))/(count+1.0)	

		   xlowh = ph2h(jlo,kt) + ( (WNOV(k)-wnc(jlo))*
     1		      (- ph2h(jlo,kt) + ph2h(jlo+1,kt))/deltastep )
		   xhih = ph2h(jhi-1,kt) + ( (WNOV(k+1)-wnc(jhi-1))*
     1		      (- ph2h(jhi-1,kt) + ph2h(jhi,kt))/deltastep )
		   xmeanh = (sumh + 0.5*(xlowh + xhih))/(count+1.0)	

		   xlowch = ph2ch4(jlo,kt) + ( (WNOV(k)-wnc(jlo))*
     1		      (- ph2ch4(jlo,kt) + ph2ch4(jlo+1,kt))/deltastep )
		   xhich = ph2ch4(jhi-1,kt) + ( (WNOV(k+1)-wnc(jhi-1))*
     1		      (- ph2ch4(jhi-1,kt) + ph2ch4(jhi,kt))/deltastep )
		   xmeanch = (sumch + 0.5*(xlowch + xhich))/(count+1.0)	

		 else
		  if (WNOV(k).gt.freqmax) then
			xmeanhh = -33.
			xmeanhe = -33.
			xmeanh = -33.
			xmeanch = -33.
		  else
			xmeanhh = piahh3(kt,k-1)
			xmeanhe = piahe3(kt,k-1)
			xmeanh  = piah3(kt,k-1)
			xmeanch = piach3(kt,k-1)
		  endif
		 endif
		   piahh3(kt,k) = xmeanhh
C                  if ((kt.eq.174).or.(kt.eq.74)) print *,kt,tshh(kt),k,wlni(k),xmeanhh
		   piahe3(kt,k) = xmeanhe
		   piah3(kt,k) = xmeanh
		   piach3(kt,k) = xmeanch
		enddo
	  enddo

C One other thing to do here is to compute the H2- opacity coef
C on the same T grid that H2-H2 opacity uses.
 	   do kt=1,ntemps3
C While we are here, also get Tristan's overtone 3 -> 0 opacity for this T
C at wavelengths where the table has no opacity
        call fit_linsky(tshh(kt),3,wnov,NSPECV,kappa3)
		do k=1,NSPECV
C                 if (piahh3(kt,k).lt.-30.) then
                     piahh3(kt,k) = dlog10(10.d0**piahh3(kt,k)+kappa3(k))
C                    print *,'',kt,k,wlni(k),kappa3(k)
C                 endif
		  if (tshh(kt).lt.600.0) then
                     h2m(kt,k) = 1.d-90
		  else
                     call opa_tab_h2mff(tshh(kt),wnov(k),h2m(kt,k))


		  endif		

          h2m(kt,k) = dlog10(h2m(kt,k))
C      if(tshh(kt).eq.800.)write(*,988)wnoi(k),piahh3(kt,k),piahe3(kt,k),piah3(kt,k),piach3(kt,k),h2m(kt,k),dlog10(kappa3(k))
C 988		format (f9.3,6f14.4)
		enddo
	  enddo

C Now we get set up to compute the PIA terms with splines

C For each spectral interval we need to call
C splint.  Need to compute first derivative
C at end points.
      DO IK=1,NSPECV
        dts1 = tshh(ntemps3)-tshh(ntemps3-1)
        dts0 = tshh(2)-tshh(1)

C H2-H2
        YPN = (piahh3(ntemps3,IK)-piahh3(ntemps3-1,ik))/dts1
        YP1 = (piahh3(2,ik) - piahh3(1,ik))/dts0
	do j=1,ntemps3


	 pass(j)=piahh3(j,ik)

	enddo
	call spline(tshh,pass,ntemps3,yp1,ypn,y2x)
	do j=1,ntemps3
	 y2hh(j,ik)=y2x(j)
	enddo

C H2-He
        YPN = (piahe3(ntemps3,IK)-piahe3(ntemps3-1,ik))/dts1
        YP1 = (piahe3(2,ik) - piahe3(1,ik))/dts0
	do j=1,ntemps3
	 pass(j)=piahe3(j,ik)
	enddo
	call spline(tshh,pass,ntemps3,yp1,ypn,y2x)
	do j=1,ntemps3
	 y2he(j,ik)=y2x(j)
	enddo

C H2-H
        YPN = (piah3(ntemps3,IK)-piah3(ntemps3-1,ik))/dts1
        YP1 = (piah3(2,ik) - piah3(1,ik))/dts0
	do j=1,ntemps3
	 pass(j)=piah3(j,ik)
	enddo
	call spline(tshh,pass,ntemps3,yp1,ypn,y2x)
	do j=1,ntemps3
	 y2h(j,ik)=y2x(j)
	enddo

C H2-CH4
        YPN = (piach3(ntemps3,IK)-piach3(ntemps3-1,ik))/dts1
        YP1 = (piach3(2,ik) - piach3(1,ik))/dts0
	do j=1,ntemps3
	 pass(j)=piach3(j,ik)
	enddo
	call spline(tshh,pass,ntemps3,yp1,ypn,y2x)
	do j=1,ntemps3
	 y2ch(j,ik)=y2x(j)
	enddo

      ENDDO
C This is a bunch of diagnostic stuff for debugging
C the splines
	goto 9999

C	print *,' ',tshh(10),tshh(11)
        do it = 125,7025,100
        ttt = dble (it)
        ik = 174
         call PIAN(IK,ttt,PHH,PHE,PH,PHC)
	print *,ik,ttt,phh,dlog10(phh)
        ik = 74
         call PIAN(IK,ttt,PHH,PHE,PH,PHC)
	print *,ik,ttt,phh,dlog10(phh)
        enddo
C	stop
        do it = 125,2025,100
 	do ik=168,178
         call PIAN(IK,ttt,PHH,PHE,PH,PHC)
         call PIAN(IK+4,ttt,PHH1,PHE1,PH1,PHC1)
         if ((phh.gt.1.).or.(phe.gt.1.)) then
           print *,'ding dong',ik,ttt,phh,phe
         endif
         if ((phc.gt.1.).or.(ph.gt.1.)) then
           print *,'ding dong',ik,ttt,phh,phe
         endif
         xxx = phh1/phh
         rat1 = 1.d5
         rat2 = 1.d-5
         if ((xxx.gt.rat1).or.(xxx.lt.rat2)) then
		print *,'phh',ik,wlni(ik),wlni(ik+4),ttt,phh,phh1
         endif
         xxx = phe1/phe
         if ((xxx.gt.rat1).or.(xxx.lt.rat2)) then
		print *,'phe',ik,ttt,phe,phe1
         endif
         xxx = phc1/phc
         if ((xxx.gt.rat1).or.(xxx.lt.rat2)) then
		print *,'phc',ik,ttt,phc,phc1
         endif
         xxx = ph1/ph
         if ((xxx.gt.rat1).or.(xxx.lt.rat2)) then
		print *,'ph',ik,ttt,ph,ph1
         endif
C        write (*,344)ik,wlni(ik),piahh3(10,ik),phh,piahh3(11,ik)
        enddo
        enddo
9999	continue
C       stop
344	format (i3,f8.3,3f14.5)
      RETURN
      END
C----------------------------------------------------------------------- 
C END SETPIAI
	  
C----------------------------------------------------------------------- 
C FIT_LINSKY (called in SETPIAI)
      
      subroutine fit_linsky(t,va,sigma,nbt,kappa)
	implicit double precision (a-h,o-z)
c Version: 04/01/96
c Fit des bandes v:0->2 et v:0-3 d'apres Linsky (1969) et Lenzuni et al. (1991)
c Auteur: T. Guillot
c Entrees:
c     va: 2 for the first overtone (+dble vibrational), 3 for the second overtone
c     sigma(nbt): wavenumber in cm-1
c Sorties:
c     kappa(nbt): H2-H2 absorption in cm-1.amagat-2

      integer va,nbt
      dimension sigma(nbt)
      integer i,j
	  real*8 kom,kappa(nbt)
      dimension sig0(3),d1(3),d2(3),d3(3),a1(3),a2(3),b1(3),b2(3)
      
c D'apres Lenzuni et al ApJS (1991) and modifications
c (TG -comparison with T=300K experimental data ):
      data sig0/4162.043,8274.650,12017.753/
c      data d1/1.2750d5,1.2750d5,1.2750d5/
      data d1/1.2750d5,1.32d6,1.32d6/
      data d2/2760.,2760.,2760./
      data d3/0.40,0.40,0.40/
c      data d3/0.40,0.60,0.60/
c      data a1/-7.661,-10.05,-11.67/
      data a1/-7.661,-9.70,-11.32/
      data a2/0.5725,0.5725,0.5725/
c      data b1/1.543,1.543,1.543/
c      data b2/0.327,0.327,0.327/
      data b1/0.9376,0.9376,0.9376/
      data b2/0.5616,0.5616,0.5616/
      
      do i=1,nbt
       kappa(i)=0.
      enddo
      
      if ((va.lt.1).or.(va.gt.3)) then
c      write(*,226)va
 226   format('FIT_LINSKY: pas de calcul possible de la transition v:',
     &      '0->0; 0->',i1)
       stop
      endif

c     write(*,116)va
 116  format('Fit de la bande H2-H2 v:0->',i1,' d''apres Linsky',
     &     ' (1969) et Lenzuni (1991)')
      
      j=va
      d=d3(j)*sqrt(d1(j)+d2(j)*t)
      a=10**(a1(j)+a2(j)*log10(t))
      b=10**(b1(j)+b2(j)*log10(t))
      aa=4.d0/13.d0*a/d*exp(1.5d0*d/b)
c     write(*,*)'sig0(j),d,a,b,aa',sig0(j),d,a,b,aa
      do i=1,nbt
       if (sigma(i).lt.sig0(j)) then
        kom=a*d*sigma(i)*exp((sigma(i)-sig0(j))/0.6952/t)/
     &       ((sigma(i)-sig0(j))**2+d*d)
       elseif (sigma(i).lt.sig0(j)+1.5*d) then
        kom=a*d*sigma(i)/((sigma(i)-sig0(j))**2+d*d)
       else
        kom=aa*sigma(i)*exp(-(sigma(i)-sig0(j))/b)
       endif
       kappa(i)=kappa(i)+kom
      enddo
      return
      end

C----------------------------------------------------------------------- 
C END FIT_LINSKY
	  
C----------------------------------------------------------------------- 
C PIAN (called in OPTCV and SETPIAI)
C JJF: This was taken from the EGP code and dropped in here 

      subroutine PIAN(IW,T,PHH,PHE,PH,PHC)
      implicit double precision (a-h,o-z)
      
      include 'prog_params'

      PARAMETER (NTEMPS3=198,nint=1000)
      
      COMMON /PIAC3/ PIAHH3(NTEMPS3,NSPECV),PIAHE3(NTEMPS3,NSPECV),
     &  PIACH3(NTEMPS3,NSPECV),PIAH3(NTEMPS3,NSPECV)
      COMMON /TEMPS3/ TSHH(NTEMPS3),TSHE(NTEMPS3)
      common /piaspline/ y2hh(ntemps3,NSPECV),y2he(ntemps3,NSPECV),
     & y2h(ntemps3,NSPECV),y2ch(ntemps3,NSPECV)
     
 
	if (iw.le.455) then
c This cuts off CIA at wavelengths less than 0.64 microns -- off of table
c 0.5 is actually the end of the table, but there is funny behavior until 0.64
         phh = 0.d0
         phe = 0.d0
         ph = 0.d0
         phc = 0.d0
         return
        endif

c JJF: Let's turn this off:  1/17/07
c      IF ((T.LT.tshh(1)).or.(T.GT.tshh(ntemps3))) THEN
c	   PRINT *,'T ERROR IN PIA',iw,T
c 	   PRINT *,tshh(1),tshh(ntemps3)
c	   T=tshh(1)
c      ENDIF

c JJF: modified so the code will not stop when it hits T<100 K
c IF T<100 K, just use the data for 100 K.
	IF(T.LT.tshh(1)) then
	T=tshh(1)
	endif

      klo = 1
      khi = ntemps3
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(tshh(k).gt.t)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif

      h=tshh(khi)-tshh(klo)
      if (h.eq.0.) pause 'bad t input in pia'
      a=(tshh(khi)-t)/h
      b=(t-tshh(klo))/h

C H2-H2
      phh=a*piahh3(klo,iw)+b*piahh3(khi,iw)+((a**3-a)*y2hh(klo,iw)+(b**3-b)*
     * y2hh(khi,iw))*(h**2)/6.d0

c	if (phh.gt.0.d0) then
c         print *,'ring ring',t,h,a,b,klo,khi,iw
c         print *,'t,h,a,b,klo,khi,iw'
c         print *,piahh3(klo,iw),piahh3(khi,iw),y2hh(klo,iw),y2hh(khi,iw)
c         print *,'piahh3(klo,iw),piahh3(khi,iw),y2hh(klo,iw),y2hh(khi,iw)'
c   endif

C H2-He
      phe=a*piahe3(klo,iw)+b*piahe3(khi,iw)+((a**3-a)*y2he(klo,iw)+(b**3-b)*
     * y2he(khi,iw))*(h**2)/6.d0

C H2-H
      ph=a*piah3(klo,iw)+b*piah3(khi,iw)+((a**3-a)*y2h(klo,iw)+(b**3-b)*
     * y2h(khi,iw))*(h**2)/6.d0

C H2-CH4
      phc=a*piach3(klo,iw)+b*piach3(khi,iw)+((a**3-a)*y2ch(klo,iw)+(b**3-b)*
     * y2ch(khi,iw))*(h**2)/6.d0

c Taking the antilog puts the units into cm^-1 amagat^-2
        PHH = 10.d0 ** PHH
        PHE = 10.d0 ** PHE
        PH  = 10.d0 ** PH
        PHC = 10.d0 ** PHC

      return
      END

C----------------------------------------------------------------------- 
C END PIAN
	  
C----------------------------------------------------------------------- 
C SETCH4 (called in GEOM)
C JJF:  Routine was called SETPIAO for reading in overtone CIA and also
C for reading in CH4.  Now only CH4 is read in--hence the rename.

      SUBROUTINE SETCH4(iprint)
      implicit double precision (a-h,o-z)

      include 'prog_params'

      PARAMETER (NTEMPS3=43,nint=333)
      
      COMMON /ATM/ Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
C      COMMON /TEMPS3/ TSHH(NTEMPS3),TSHE(NTEMPS3)
      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      COMMON /XKAP/ XKCH4(790)
      COMMON /PIAOC/ PIAHHO(NTEMPS3,NSPECV),PIAHEO(NTEMPS3,NSPECV)
      dimension wnc(nint),phho(nint,NTEMPS3),pheo(nint,NTEMPS3)
      COMMON /ERICH/ ech(nspecv)
      dimension wne(3766),ech4(3766)
      double precision pvalue,f(nspecv)


C Now do CH4
c First zero the array
	do k=1,nspecv
	  ech(k) = 0.d0
	enddo

	  print *,' Erich Karkoschka CH4 opacity for use in IR (km-am)^-1'
	  print *,' shortwards of 1 micron.  Freedman table used redward of 1 um'
C Read info from a file
	open (8,file='../inputs/ERICH.DAT')
C        open (9,file='ch4opac.dat')
C        open (10,file='kappch4.dat')
C Number of Karkoshka data points
	ninte = 1750
C Number of Strong data points
c	nints = 662
c Total number of data points
c	nintt = ninte + nints
      	DO  kk=ninte,1,-1
 		read (8,*) wne(kk),x,ech4(kk),x1,x2,x3,x4,x5
C Convert to wavenumber
		wne(kk) = 1.d7/wne(kk) 
C Convert from per km-am to cm2/g
		ech4(kk) = ech4(kk)/71.80d0
		if (wne(kk).lt.1.d4) ech4(kk)=0.d0  !zero out redward of 1 um
C        write(9,12) wne(kk),ech4(kk) 
      	enddo 
 	close(8)

C 3/2008  We stopped using this part of the file.  
C All CH4 redward of the optical
C is now included in the tabulated opacities.
C
c Now read in the Strong CH4 opacity
c	open (8,file='Strong.short')
c      	DO  kk=1,nints
c 	  read (8,*) wne(kk),x,ech4(kk),x1,x2,x3,x4
c Convert from cm^2/molecule to cm^2/g
c       (a 1d-24 has been removed from numerator and denominator)
c	  ech4(kk) = ech4(kk)/(16.d0*1.661d0)
c Now we don't use Strong longward of 6189 wn
c	  if (wne(kk).lt.1.d4) ech4(kk)=0.d0
c      	enddo 
c        close(8)
c        do kk=ninte,1,-1
c        write(9,12) wne(kk),ech4(kk)
c        enddo
c Now read in Dawn's opacities for ch4.
c        open (8,file='ch4n.dat')
c        do kk=ndawn,1,-1
c          read (8,*) wne(kk),ech4(kk)
c Convert from cm^2/molecule to cm^2/g
c        ech4(kk)=ech4(kk)/(16.d0*1.661d0*1d-24)
c          write(9,12) wne(kk),ech4(kk)
c        enddo
c        close(8)
c        close(9)

        call interp(wne,ech4,ninte,f)
        do i=1,790    !1 to 790 is 0.35 to 1.00 microns
          xkch4(i)=f(i)
C          write (10,12) wnov(i),xkch4(i)
        enddo
C        close(10)

c	do kk=1,nintt
c		print *,kk,wne(kk),ech4(kk)
c	enddo
c		freqmax = wne(nintt)
c		freqlo = wne(1)
c		do k=1,nspecv
c	 	 xmeanch=0.0
c		 if ((bwnv(k+1).lt.freqmax).and.(bwnv(k).gt.freqlo)) then
c		  j = 1
c		  do while (wne(j).lt.bwnv(k)) 
c			j = j+1
c                  enddo
c		  jlo = j - 1 
c		  j = j - 1 
c		  do while (wne(j).lt.bwnv(k+1)) 
c			j = j+1
c                  enddo
c		  jhi = j 
c		   sumch = 0.0
c		   count = 0.0
c
c		  do jc = jlo+1, jhi-1
c		   sumch = sumch + ech4(jc)
c		   count = count + 1.0
c		  enddo
c		 ech(k) = sumch /count	
c		endif
c		enddo
c		print *,'CH4 opacity spectrum'
c 		do k=1,nspecv
c  		  print *,k,wlnv(k),ech(k)
c		enddo
c		stop
      RETURN
      END
C----------------------------------------------------------------------- 
C END SETCH4
	  
C----------------------------------------------------------------------- 
C FINDEXACT (does not appear to be currently used)
	
      SUBROUTINE FINDEXACT(xx,values,n,ntemp,f)
      implicit double precision (a-h,o-z)

      include 'prog_params'

      PARAMETER (NTEMPS3=43,nint=466)
      
      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      dimension xx(n),values(n,ntemps3),f(nspecv)
      integer n
      double precision f1,f2,x,x1,x2
      do 10 i=nspecv,1,-1
        x=wnov(i)
        call Locate(xx,n,x,j)
        x1=xx(j)
          if (j.eq.n) then
             f(i)=0.0d0
c            f1=values(j-1,ntemp)
c            f2=values(j,ntemp)
c            x1=xx(j-1)
c            x2=xx(j)
            goto 10
              else
              x2=xx(j+1)
          endif
        f1=values(j,ntemp)
        f2=values(j+1,ntemp)
   20   if ((f1.eq.0.0).and.(f2.eq.0.0)) then
          f(i)=0.0d0
          goto 10
            else
              if (f1.eq.0.0) f1=1.0D-100   
              if (f2.eq.0.0) f2=1.0D-100
            f(i)=((x-x1)/(x2-x1)*log(f2)+(x2-x)/(x2-x1)*log(f1))
     & /log(10.0d0)
            f(i)=10.0d0**f(i)
        endif
   10 continue
      end
C----------------------------------------------------------------------- 
C END FINDEXACT
	  
C----------------------------------------------------------------------- 
C INTERP (called in OPTCV)

      SUBROUTINE INTERP(xx,values,n,f)
      implicit double precision (a-h,o-z)
      
      include 'prog_params'

      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      dimension xx(n),values(n),f(nspecv)
      integer n
      double precision f1,f2,x,x1,x2
      do 10 i=nspecv,1,-1
        x=wnov(i)
        call locate(xx,n,x,j)
        x1=xx(j)
          if (j.eq.n) then
            f1=values(j-1)
            f2=values(j)
            x1=xx(j-1)
            x2=xx(j)
            goto 20
              else
              x2=xx(j+1)
          endif
        f1=values(j)
        f2=values(j+1)
   20     if ((f1.eq.0.0).and.(f2.eq.0.0)) then
            f(i)=0.0d0
            goto 10
              else
                if (f1.eq.0.0) f1=1.0D-100
                if (f2.eq.0.0) f2=1.0D-100
             f(i)=((x-x1)/(x2-x1)*log(f2)+(x2-x)/(x2-x1)*log(f1))
     &  /log(10.0d0)
             f(i)=10.0d0**f(i)
c We can't do log interpolate with negative numbers!  Check for this and fix.	     
	         if ((f1.lt.0.d0).or.(f2.lt.0.d0)) then
		 f(i)=((x-x1)/(x2-x1)*(f2) + (x2-x)/(x2-x1)*(f1))
	         endif
	     
         endif
   10 continue
      end

C----------------------------------------------------------------------- 
C END INTERP
	  
C----------------------------------------------------------------------- 
C GETOPAC (called in OPTCV)

c  This subroutine takes all of the opacity data for ch4, co, h2o, and nh3, for 
c 0.3 to 5 microns and for a specific Temperature and Pressure, gives the 
c opacities at each wavelength.  The output is a two column file.  Column 1 is
c the wavelength, and column 2 is the interpolated opacity for each of the 
c molecules.  ch4 is listed first, co second, h20 third, and nh3 fourth. 

      subroutine getopac(t,p)
      parameter (lines=2000)
      implicit double precision (a-h,o-z)
      dimension filename(4), temp(6), pr(4), op(6,4)
      dimension filen(24) 
      character*9 filename
      character*21 filen
      double precision kappa, lop(6,4)
      double precision lam(lines)

c      print *,'output',t,p
      p=p*1000.		!converts bars to millibars
      filename(1)='ch4.out'
      filename(2)='co.out'
      filename(3)='h2o.out'
      filename(4)='nh3.out'

c  Set up temperature and pressure arrays.

      temp(1)=100.
      temp(2)=200.
      temp(3)=400.
      temp(4)=600.
      temp(5)=800.
      temp(6)=1000.

      pr(1)=10.
      pr(2)=100.
      pr(3)=1000.
      pr(4)=10000.


c Check to see that the point requested is in the range of the data.  If not,
c give opacities to the closest point in the range of the data.

      if (t .lt. temp(1)) then
c       print*,'Temp too small, opacity will be given for t=100K.'
       t=100.
      endif
      if (t .gt. temp(6)) then
c       print*, 'Temp too large, opacity will be given for t=1000K.'
       t=1000.
      endif
      if (p .lt. pr(1)) then
c       print*,'Pressure too small, opacity will be given for p=10mb.'
       p=10.
      endif
      if (p .gt. pr(4)) then
c       print*,'Pressure too large, opacity will be given for p=10000mb.' 
       p=10000.
      endif

c Read in opacity files.  The opacities are read into a 2-D array (6,4).  
c It iterates over the number of lines in the opacity file.  Here the files are
c from 0.3 to 5 microns.  The spacing is every 4 angstroms from 0-1 micron
c and every 25 angstroms for 1-5 microns. 

      open (unit=33, file='tpint75.out')
      do i=1,4
       open (unit=41, file=filename(i)) 
       do if=1,lines
        if (if .gt. 1) rewind (unit=41)
        l=1
        m=1
        ir=1
        is=1
        ku=1
        iu=1  
        do j=1,24
         read (41,*) filen(j)
         open (unit=j+42, file=filen(j))
         if (j .lt. 5) then
          read(j+42,*) lam(if), op(1,l)
          l=l+1
         elseif (j .lt. 9) then
          read(j+42,*) lam(if), op(2,m)
          m=m+1
         elseif (j .lt. 13) then
          read(j+42,*) lam(if), op(3,ir)
          ir=ir+1
         elseif (j .lt. 17) then
          read(j+42,*) lam(if), op(4,is)
          is=is+1
         elseif (j .lt. 21) then
          read(j+42,*) lam(if), op(5,ku)
          ku=ku+1
         else
          read(j+42,*) lam(if), op(6,iu)
          iu=iu+1
         endif
        enddo

        do k=1,6
         do nq=1,4
          if (op(k,nq) .lt. 1D-45) then
           op(k,nq)=1D-45
          endif
          lop(k,nq)=dlog10(op(k,nq))
         enddo
        enddo

c Take the log of the opacities.  If they are too small, set them to a 
c reasonable number.

c Calculate the interpolated opacites for each wavelength and each molecule
c based on the temp and pressure desired.


        call getkappa (temp, pr, t, p, lop, kappa)
        if (kappa .ne. 0.0) then    
         kappa=10.**(kappa)
        else 
         kappa=1.0
        endif
        write (33,*) lam(if), kappa
       enddo
      do j=1,24
       close (j+42)
      enddo
      close (41)
      enddo
      close (33)
      return
      end


C----------------------------------------------------------------------- 
C END GETOPAC
	  
C----------------------------------------------------------------------- 
C GETKAPPA (called in GETOPAC, which is called by OPTCV)

      subroutine getkappa(temp,pr,t,p,lop,kappa)
      implicit double precision (a-h,o-z)
      double precision kappa, lop(6,4)
      dimension temp(6),pr(4),prlog(4),templg(6)
      call locate(temp,6,t,jt)
      call locate(pr,4,p,jp)
      if (jt.eq.6) jt=5
      if (jp.eq.4) jp=3
      plg=dlog10(p)
      tlg=dlog10(t)
      do i=1,6
        templg(i)=dlog10(temp(i))
      enddo
      do i=1,4
        prlog(i)=dlog10(pr(i))
      enddo
      x1=templg(jt)
      x2=templg(jt+1)
      f1=lop(jt,jp)
      f2=lop(jt+1,jp)
      to=(tlg-x1)/(x2-x1)*f2+(x2-tlg)/(x2-x1)*f1
      f1=lop(jt,jp+1)
      f2=lop(jt+1,jp+1)
      tf=(tlg-x1)/(x2-x1)*f2+(x2-tlg)/(x2-x1)*f1
      x1=prlog(jp)
      x2=prlog(jp+1)
      kappa=(plg-x1)/(x2-x1)*tf+(x2-plg)/(x2-x1)*to
      end

C----------------------------------------------------------------------- 
C END GETKAPPA
C-----------------------------------------------------------------------

C----------------------------------------------------------------------- 
C LOCATE (called in INTERP, FINDEXACT, and GETKAPPA)
   
      SUBROUTINE LOCATE(XX,N,X,J) ! modified to correct bugs:Tue May 16 2006
                                  ! from online edition of numerical recipies
C Added by JJF -- same as old LOCATE, but now has correction     
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GE.XX(1)).EQV.(X.GE.XX(JM)))THEN ! corrected
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
C Added code from correction
      if(x.eq.xx(1))then ! Then set the output
         j=1 
      else if(x.eq.xx(n))then
         j=n-1 
      else
         j=jl 
      endif 
      RETURN
      END
C----------------------------------------------------------------------- 
C END LOCATE
	  
C----------------------------------------------------------------------- 
C END OF FILE
