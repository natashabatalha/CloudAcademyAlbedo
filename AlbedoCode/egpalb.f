      PROGRAM GEOMETRIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
C PROGRAM TO CALCULATE THE GEOMETRIC AND BOND ALBEDOS 
C (Originally for Titan)
C
C Requires files: egpalb.f, subs_egp.f, raman2.f, trist.f
C
C Albedo output to file, requires dummy file exists in run dir
C 00_SPECTRAL_ALBEDOS.TXT, 00_DATA_SPEC_VIS.TXT
C
C THE GEOMETRIC ALBEDO IS COMPUTED TWO WAYS: FROM THE SOURCE FUNCTION
C METHOD DEVELOPED BY TOON ET AL. AND FROM THE Q-FIT METHOD IN MCKAYS
C NOTES. 
C
C TGEOM is currently the output Geometric Albedo, egp.f --> GEOMALB
C File I/O in subsegp.f --> ATMSETUP, *.in, *.cld
C Also check GNEFF for effective wavelength for planet model in egp.f
C Also check RSFV, WBAR, COSBV, GCOS2 definitions (subs_egp.f)
C Also check type of GINT (egp.f)
C Phase integral output currently off, egp.f --> GEOM and phase
C 
C
C File OPENs commented with C...........
C Output WRITEs commented with C~~~~~~~~~~~
C Subroutine/function/block data separation C-----------
C
C
C
C                           CODE DEVELOPED JAN 87
C                           C.P. McKAY
C-----------------------------------------------------------------------
C Edited by Jonathan Fortney (November 2006) to
C include pretabulated opacity data from R. Freedman
C Edited by Jonathan Fortney (March 2008) so that
C now pretabulated opacity includes absolutely everything except
C clouds, Rayleigh scattering, and optical CH4 (from Karkoschka)
C-----------------------------------------------------------------------
C Edited by Kerri Cahoy (circa March 2009) to introduce GINT3
C subroutine to include different observer + sun angles (mu0, mu1)
C-----------------------------------------------------------------------
C Edited by Nikole Lewis (circa April 2014) to use the 
C prog_params pangle_params files to simplify operations
C-----------------------------------------------------------------------
C Edited by Nikole Lewis (circa December 2015) to correct numerical 
C issues
      
      include 'freedman/int.params'
      include 'freedman/abs_data.com'
      include 'freedman/limits.com'
      include 'prog_params'
      include 'pangle_params'

      PARAMETER (NTEMPS3=43,nint=466)
 
C Planetary phase angle in radians (0 = face on)      
      COMMON /OBSPHA/ PLPHA, PANGLE(NPHAS)
      COMMON /UBARED/ UBARI,UBARV,UBAR0
      COMMON /UBARED2/ UBAR1, DPHI
      COMMON /UMTXS/ U0MTX,U1MTX,NUMTX,DPMTX
C The UBAR1 is circa March 2009, put in separate common variable	  
      COMMON /ATM/ Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
      COMMON /LAPSE/ DTDP(NLAYER),CONVEQ
C      COMMON /GASS/ CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),CO(NLEVEL)
C     & ,HE(NLEVEL),XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER)
C     & ,H2O(NLEVEL),XNH3(NLEVEL),H(NLEVEL),Hmin(NLEVEL),elec(NLEVEL)
      COMMON /GASS/ CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),CO(NLEVEL),
     &  CO2(NLEVEL),HE(NLEVEL),XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER),
     &  H2O(NLEVEL),XNH3(NLEVEL),H(NLEVEL),Hmin(NLEVEL),elec(NLEVEL), 
     &  xNA(NLEVEL), xK(NLEVEL),H2S(NLEVEL),TIO(NLEVEL), VO(NLEVEL),
     &  RB(NLEVEL), CS(NLEVEL), FEH(NLEVEL), CRH(NLEVEL) 
      COMMON /VISGAS/ SOLARF(NSPECV)
      COMMON /AERSOL/ RADIUS(NLAYER), XNUMB(NLAYER), SIGMA(NLAYER) 
     & ,  REALV(NSPECV), XIMGV(NSPECV)
      COMMON /CLOUD/ RADCLD(NLAYER), XNCLD(NLAYER),XSCTC
     & , RCLDV(NSPECV), XICLDV(NSPECV)
      COMMON /CLOUDD2/ XNCLD2(NLAYER)
      COMMON /TAUS/ 
     &  TAURV(NSPECV),TAUHV(NSPECV),TAUCV(NSPECV),TAUGV(NSPECV)
     &  ,TAUC2V(NSPECV),TAUCLV(NSPECV)
      COMMON /OPTICV/ DTAUV(NLAYER,NSPECV),TAUV(NLEVEL,NSPECV)
     &                 ,WBARV(NLAYER,NSPECV), COSBV(NLAYER,NSPECV)
     &                 ,NOPTHCK(NSPECV), LAYERBOT(NSPECV)
      COMMON /RAY2/ GCOS2(NLAYER,NSPECV)
      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      COMMON /FLUXV/ FNETV(NLEVEL), FUPV(NLEVEL,NSPECV),
     &  FDV(NLEVEL,NSPECV), FMNETV(NLEVEL),FMUPV(NLEVEL),FMDV(NLEVEL)
      COMMON /CLEAR/ CLA(NSPECV),ICLEAR
      COMMON /PLANT/ CSUBP,RSFI,RSFV,F0PI
      COMMON /ADJUST/ RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,RCLOUD
      COMMON /CONST/RGAS,RHOP,PI,SIGSB,RTOD
      COMMON /PLOT1/ VAR(NLEVEL),fink(4)
      COMMON /TESTS/FHZE,FRAY,FCLDC

      DIMENSION WNEFF(702),XNEFF(702)
      DIMENSION BOND(NSPECV),QGEOM(NSPECV),Q1(NSPECV),TGEOM(NSPECV)
      DIMENSION GWGHT(8),GANGL(8),XINT(NLEVEL),rn1(nlayer),phi(nspecv)
c      DIMENSION GWGHT2(10),GANGL2(10),TWGHT2(10),TANGL2(10)  ! for GEOMALB2
      DIMENSION GWGHT3(100),GANGL3(100),TWGHT3(10),TANGL3(10)  ! for GEOMALB3
c      REAL U0MTX(10,10),U1MTX(10,10),NUMTX(10,10)   ! for GEOMALB2
      REAL U0MTX(100,10),U1MTX(100,10),NUMTX(100,10)  ! for GEOMALB3
C     DIMENSION TESTXINT(8) !for use with GEOMALB
c      DIMENSION TESTXINT(100) !for use with GEOMALB2
      DIMENSION TESTXINT(1000) !for use with GEOMALB3
      
c      DIMENSION SLABGEO(NSPECV,10,10) !for GEOMALB2
c      DIMENSION SLABINT(NSPECV,10,10) !for GEOMALB2
      DIMENSION SLABGEO(NSPECV,100,10) !for GEOMALB3
      DIMENSION SLABINT(NSPECV,100,10) !for GEOMALB3

C      DIMENSION OUTXINT(NSPECV,8) !for use with GEOMALB
c      DIMENSION OUTXINT(NSPECV,100) !for use with GEOMALB2
      DIMENSION OUTXINT(NSPECV,1000) !for use with GEOMALB3

      CHARACTER * 15  MESS 
      REAL * 4  RTGM(NSPECV),RWLV(NSPECV)
      REAL * 4  RWNEF(702),RNEFF(702),WIUE(27),XIUE(27)

      CHARACTER(80) TFNAME
      CHARACTER(80) TFNAME2

 
      write(*,*)'NLAYER=',NLAYER
      write(*,*)'NLEVEL=',NLEVEL
      write(*,*)'Number of Phase Divisions:',NPHAS
      write(*,*)'Check:',SIZE(PANGLE)
      
      WRITE(6,1111)
 1111 FORMAT(' GEOMETRIC ALBEDO CALCULATED WITH TWO-STRM SOURCE FUNCTION
     & '/'  BASED ON THE DELTA-EDDINGTON METHOD OF TOON ET AL.
     & ')

      IPRINT=3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO L=1,NPHAS
      PLPHA = PANGLE(L)
      print *,'L= ', L
      print *,'PLPHA= ', PLPHA
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SET UP BACKGROUND ATMOSHPERE 
C Haze column density mult. factor
      fhze = 0.
C Visible surface albedo (also is a data variable later, clean up?)
      rsfv=0.0
C fcldc adjusts cloud g
      fcldc = 1.0

C      Old version of this call
C      CALL ATMSETUP(NLEVEL,Z,TEMP,PRESS,DEN,XMU,
C     & CH4,H2,HE,XN2,CO,H2O,XNH3,H,Hmin,elec,IPRINT)

       CALL ATMSETUP(NLEVEL,Z,TEMP,PRESS,DEN,XMU,
     & CH4,H2,HE,XN2,CO,H2O,CO2,XNH3,H,Hmin,elec,xNa,xK,RB,CS,FEH,CRH,
     & TIO, VO, x736,x1060,IPRINT)
      
        
C NOW CALCULATE THE LAYER AVERAGE GAS MIXING RATIOS.
C OF THE ABSORBING GAS
C AND THE TOTAL LAYER COLUMN MASS (units?)
      DO 159 J=1,NLAYER
      EMU=(XMU(J+1)+XMU(J))*0.5d0 !average mean molecular weight
      COLDEN(J)=RHOP*(PRESS(J+1)-PRESS(J))/EFFG(1.) !column density
      GAS1(J)=(16./EMU)*AVERGE(CH4(J+1),CH4(J)) 
159   CONTINUE

C CALL A ROUTINE THAT SETS UP THE VIS SPECTRAL INTERVALS
      CALL SETSPV(IPRINT)
C Call a routine that sets up pressure-induced absorption
      CALL SETPIAI
C Call a routine that sets up methane	 
      CALL SETCH4

C CALL ROUTINES THAT SET UP THE AEROSOL DIST AND OPTICAL PROP.
c     CALL HAZESET(IPRINT)
c     CALL HAZESCAT(IPRINT)
c     CALL CLOUDSCAT
c     CALL CLD(IPRINT)
c     CALL CLOUD2SCAT

C Density and average molecular mass updates
      DO 9159 J=1,NLAYER
      EMU=(XMU(J+1)+XMU(J))*0.5d0 !layer avg
      COLDEN(J)=RHOP*(PRESS(J+1)-PRESS(J))/EFFG(z(j))
      rn1(j)=averge(ch4(j+1),ch4(j)) 
      GAS1(J)=(16.0426d0/EMU)*rn1(j)
9159   CONTINUE

C CALL A SUBROUTINE THAT SETS UP THE OPTICAL PROPERTIES IN THE 
C VISIBLE. 
      CALL OPTCV(IPRINT)

C PRINT OUT THE AEROSOL AND CLOUD INFORMATION
      IF (IPRINT .GT. 0) THEN
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         WRITE (6,139)
         DO 135 J=1,NLAYER
         WRITE(6,140)J,PRESS(J),TEMP(J),
     &    CH4(J)*PRESS(J)/PCH4(TEMP(J)),CH4(J)*100.,RADCLD(J),XNCLD(J)
     &    ,RADIUS(J),XNUMB(J),RADCLD(J),XNCLD2(J)
  135    CONTINUE
  139 FORMAT(///'   AEROSOLS AND CLOUDS AT LEVELS'/
     &' LVL  P(BARS)  TEMP  RH-CH4'
     & , '  %CH4  CH4 CLOUD:SIZE  COL NUMB   HAZE:SIZE  COL NUM',
     &  '    H2S CLOUD:SIZE  COL NUMB')
  140 FORMAT (1X,I3,1PE10.3,0P2F7.2,F6.2,3(6X,0PF7.3,1PE10.3))
          ENDIF
      IF (IPRINT .GT. 1) THEN
C PRINT OUT SPECTRAL INTERVALS IN THE VISIBLE
        WRITE (6,190) 'VIS '
        WRITE (6,210)(K,WLNV(K),WNOV(K),TAUHV(K),TAUCLV(K),
     &  TAUCV(K),TAUC2V(K),TAURV(K),TAUGV(K),TAUV(NLEVEL,K)
     &  ,XIMGV(K),RCLDV(K), K=1,NSPECV)
 
C PRINT OUT SPECTRAL INTERVALS IN THE VISIBLE

        
        WRITE(TFNAME, '(a,I3.3,a)') '00_Data_Spec_Vis_', idnint(PLPHA*RTOD), '.txt'
        OPEN(7, FILE=TFNAME, STATUS='REPLACE')
        WRITE (7,190) 'VIS '
        WRITE (7,210)(K,WLNV(K),WNOV(K),TAUHV(K),TAUCLV(K),
     &  TAUCV(K),TAUC2V(K),TAURV(K),TAUGV(K),TAUV(NLEVEL,K)
     &  ,XIMGV(K),RCLDV(K), K=1,NSPECV)
        CLOSE(7)
 
      END IF
  210 FORMAT(1X,I4,F10.5,F9.1,1P9E10.3,0P)
  190 FORMAT(///'  SPECTRAL INTERVALS IN ',A4//
     &       '   MICRONS   WAVENU   TAU HAZE',
     & ' TAU CLEAR  TAU CLOUD:CH4    H2S  TAU RAY  TAU GAS   TAU-TOTAL ' 
     &   ,'HAZE:IMG,    CLOUD:REAL')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


C CALL A ROUTINE THAT GETS THE NET FLUXES IN THE VISIBLE
C INTEGRATING OVER ALL SPECTRAL INTERVALS USING THE D-EDDINGTON
      CALL SFLUXV(IPRINT)

C OUTPUT THE SOLAR CONSTANT
      SOLARC=0.0
      DO 552 K=1,NSPECV
      SOLARC=SOLARC+SOLARF(K)
 552  CONTINUE
C Not currently used
      SOLALB=1.0 + 2.* FNETV(1)/SOLARC

       IF (IPRINT .GT. 0 ) THEN 
C OUTPUT THE VISIBLE FLUXES AND PRESENT NET FLUX.
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       WRITE (6,892)
       DO J=1,NLEVEL
       WRITE(6,8922) J,REAL(TEMP(J)),REAL(PRESS(J))
     &  ,FNETV(J),FNETV(J)/(SOLARC*UBAR0)
       enddo
892    FORMAT(/' LEVEL   TEMP    PRESS  ',
     & 'VIS FNET    FNET/SO  ')
8922    FORMAT(2X,I3,2F10.4,2E15.5)
C OUTPUT SELECTED SOLAR FLUX DATA FOR ALL SPECV'S
C NOTE THAT NLEVEL SHOULD NOT DE DIVISIBLE BY 6 OR FORMAT ERROR HERE
       WRITE(96,289)'   UP','SOLAR',(J,J=1,NLEVEL,NLEVEL/8),(K,WLNV(K),
     &  (FUPV(J,K),J=1,NLEVEL,NLEVEL/8),  K=1,NSPECV)
       WRITE(96,289)' DOWN','SOLAR',(J,J=1,NLEVEL,NLEVEL/8),(K,WLNV(K),
     &   (FDV(J,K),J=1,NLEVEL,NLEVEL/8),  K=1,NSPECV)
  289  FORMAT (//A5,'WARD ',A5,' FLUXES',12X,'LEVELS'/
     & ' WAVELENGTH',9I10,/(1X,I4,F9.5,1P9E10.3,0P))
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ENDIF

C COMPUTE THE BOND AND Q-BASED GEOMETRIC ALBEDOS:

        
C THE Q-FIT IS DESCRIBED IN C. McKay's NOTES (but not really used)
        DO K=1,NSPECV
C Not really Bond... not integrated over wavelength
        BOND(K)=2.*FUPV(1,K)/SOLARF(K)
C		print *,'K=', K, '  SOLARF=', SOLARF(K)
C The expression below is an old estimate of the 
C geometric albedo using an expression for its phase
C function (the denominator) that is wlnv-dependent
C and approximately Lambertian (3/2) at short wlnv
C However, TGEOM is reset to 0 in GEOMALB again... 
C So not sure what the point of leaving this in is?

        QGEOM(K)=BOND(K)/(1.5 + 1.27*DEXP(-2.74*WLNV(K) ) )
        TGEOM(K)=QGEOM(K) 
        ENDDO

C NOW COMPUTE THE GEOMETRIC ALBEDO FROM THE TWO STREAM SOURCE
C METHOD OF TOON ET AL.



      
      
      DO 148 K=1,NSPECV
  
c      print *,'K= ', K
C      CALL GEOMALB(NLEVEL,DTAUV(1,K),TAUV(1,K),WBARV(1,K)
C     &   ,COSBV(1,K),gcos2(1,k),TGEOM(K),TESTXINT)


C      CALL GEOMALB2(NLEVEL,DTAUV(1,K),TAUV(1,K),WBARV(1,K)
C     &   ,COSBV(1,K),gcos2(1,k),TGEOM(K),TESTXINT
C     &   ,SLABGEO(K,1:10,1:10),SLABINT(K,1:10,1:10))
     
      CALL GEOMALB3(NLEVEL,DTAUV(1,K),TAUV(1,K),WBARV(1,K)
     &   ,COSBV(1,K),gcos2(1,k),TGEOM(K),TESTXINT
     &   ,SLABGEO(K,1:100,1:10),SLABINT(K,1:100,1:10))

C For geomalb3     
      DO M=1,1000
C For geomalb2     
C      DO M=1,100
C For geomalb      
C      DO M=1,8
c      print *,'M =', M
c      print *,'TESTXINT OUTSIDE = ', TESTXINT(M)
      OUTXINT(K,M) = TESTXINT(M)
c      print *,'OUTXINT = ', OUTXINT(K,M)
      ENDDO
     
  148 CONTINUE

 
C CALL ROUTINES TO ESTIMATE THE PHASE INTEGRAL 
C Set up ability to take second derivatives
       call interpset()
            do 1822 i=1,nspecv
C Phi is the phase integral and phf is the "phase function."
            j=nopthck(i)
            jbot=layerbot(i)
C Check to see if rayleigh?		
            if (cosbv(j,i).eq.0.0) then
              cosbv(j,i)=-1.0
            endif
            call phase(wbarv(j,i),cosbv(j,i),phi(i))

C Subroutine definition returns "phi": phase(ssa,phf,phi)
C wbarv is single scattering albedo visible
C cosbv is "g" or by this acronym, "phase function"
C phi is phase integral

C The output files are currently turned off, apparently...			
C            write (75,*) i,sngl(wlnv(i)),sngl(tgeom(i)),sngl(phi(i))
C            write (750,*) sngl(wlnv(i)),sngl(tgeom(i)),j,jbot
C     &           ,sngl(cosbv(j,i)),sngl(wbarv(j,i)),sngl(phi(i))
1822       continue

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Write albedos to screen
       WRITE(6,782) (K,WLNV(K),GEOMAL(WLNV(K)),TGEOM(K),QGEOM(K),
     & BOND(K)/TGEOM(K),1.5 + 1.27*DEXP(-2.74*WLNV(K) ),
     &      K=1,NSPECV)

       WRITE(TFNAME, '(a,I3.3,a)') '00_Spectral_Albedos_', idnint(PLPHA*RTOD), '.txt'
       OPEN(19, FILE=TFNAME, STATUS='REPLACE')

C Write albedos to file
       WRITE(19,782) (K,WLNV(K),GEOMAL(WLNV(K)),TGEOM(K),QGEOM(K),
     & BOND(K)/TGEOM(K),1.5 + 1.27*DEXP(-2.74*WLNV(K) ),
     &      K=1,NSPECV)	 
 782  FORMAT(//' SPECTRAL ALBEDOS:'/'  K    WAVELN    NEFF   '
     & ,'    GEOMALB   BOND/QFIT  BOND/GEOM  Q-FITTED'/(I4,6(F11.7,2X)))    
	  CLOSE(19)	  
	  
C Write xint to file
	  WRITE(TFNAME, '(a,I3.3,a)') '00_XINT_TGeom_', idnint(PLPHA*RTOD), '.txt'
	  OPEN(25, FILE=TFNAME, STATUS='REPLACE')
C      print *,((OUTXINT(K,M),M=1,8),K=1,NSPECV)	
      write(25,781)
 781  FORMAT('  K, WAVELEN, xint1 for all gangls'/)      
      do 784 K=1,NSPECV
C For geomalb3      
c      write(25,783) K,WLNV(K),(OUTXINT(K,m),m=1,1000)
C For geomalb2      
c      write(25,783) k,WLNV(k),(OUTXINT(k,m),m=1,100)
C For geomalb      
C      write(25,783) k,WLNV(k),(OUTXINT(k,m),m=1,8)
 784  continue
 783  FORMAT(I4,F11.7,100E12.3)
      CLOSE(25)
C-----------------------------------------------------------------------
c       DO NGP=1,100 !with GEOMALB3

c      DO NTP=1,10
c      WRITE(TFNAME, '(a,I2.2,a,I3.3,a,I2.2,a,a)') 'slab_files/SLAB_Albedos', 
c     &         L,'P', NGP,'G', NTP,'T', '.txt' 
c      OPEN(88, FILE=TFNAME, STATUS='REPLACE')

C Write SLAB albedos to file
c       WRITE(88,988) (K,WLNV(K),SLABGEO(K,NGP,NTP),
c     &      K=1,NSPECV) 
c 988  FORMAT(//' SPECTRAL ALBEDOS:'/'  K    WAVELN   GEOM ALB	INT   '
C     & ,/(I4,2(F11.7,2X)))     
C      CLOSE(88)
   
C      ENDDO
C      ENDDO
    
C-----------------------------------------------------------------------      
      
      ENDDO ! on PLPHA

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	  
         STOP 'GEOM DONE'
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------	 
C End GEOMETRIC

C-----------------------------------------------------------------------            
C Block data; define constants
 
      BLOCK DATA  TGMDAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      include 'prog_params'
       
C Planetary phase angle in radians (0 = face on)      
C      COMMON /OBSPHA/ PLPHA, PANGLE(NPHAS)
      COMMON /UBARED/ UBARI,UBARV,UBAR0
      COMMON /UBARED2/ UBAR1, DPHI
      COMMON /LAPSE/ DTDP(NLAYER),CONVEQ
      COMMON /PLANT/ CSUBP,RSFI,RSFV,F0PI
      COMMON /ADJUST/ RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,RCLOUD
      COMMON /CONST/RGAS,RHOP,PI,SIGSB,RTOD
      COMMON /CLEAR/ CLA(NSPECV),ICLEAR
      COMMON /TESTS/ FHZE,FRAY,FCLDC

      DATA PI/3.14159265358979323846/
C RGAS IS THE UNIVERSAL GAS CONSTANT IN UNITS OF: M SEC-2 AMU K-1 KM
      DATA RGAS/8.31432/
C SIGMA IS THE STEFAN-BOLTZMAN CONSTATNT IN CGS UNITS
      DATA SIGSB/5.6677E-5/
C RTOD is Radian to Degrees
      DATA RTOD/57.295780/
C RT CONSTANTS

                     
C INCIDENT AND OBS ANGLES RADIANS
C THESE ARE USED FOR "BOND" CALC IN GEOM
C NEED TO BE NON-ZERO? (Not really useful right now)
      DATA UBAR0,UBAR1/0.5,0.5/

C PLANET SPECIFIC CONSTANTS  
C CONVEQ IS THE ADIABAT DLNT/DLNP
C     DATA CONVEQ/0.0450/
      DATA CONVEQ/0.3250/

C CSUBP IS THE SPECIFIC HEAT AT CONSTANT P OF THE ATMOSPHERE
C IN UNITS OF ERGS K-1 G-1 
      DATA CSUBP/1.05E7/

C RHOP IS THE UNITS CONVERSION FROM TO GET MASS UNITS (G CM-2)
C FROM PRESSURE (BARS) DIVIDED BY GRAVITY (M SEC-2)
C IS EQUAL TO ONE  GM CM-2  BARS-1  M SEC-2
C IF ONE WANTS TO CHANGE UNITS ON PRESSURE THIS
C CONSTANT MUST BE CHANGED
      DATA RHOP/1.E4/    

C RSF IS THE SURFACE REFLECTANCE FOR VIS AND IR
C BLACK = 0, WHITE = 1
      DATA RSFV,RSFI/0.0,0.0/
C      DATA RSFV,RSFI/1.0,0.0/

      
C FOPI IS THE ACTUAL SOLAR FLUX IN ERGS/CM2
C CHECK THAT CORRESPONDS TO PLANET OF INTEREST
C      DATA F0PI/1.5E4/
C This is set to unity later in the code (NKL)

C FHIR IS THE HAZE INFRARED ABSORPTION SCALE FACTOR
      DATA FHIR/0.5/
C FHVIS IS THE HAZE INFRARED ABSORPTION SCALE FACTOR
      DATA FHVIS/1.333333333/
C ICLEAR = 1 means use the absorbing material in the clear region
C between the clouds.
      DATA ICLEAR/1/
C FHZE is a grey multiplicative factor that changes the haze
C optical depth.  It is used for sensitivity tests.  Nominally
C it should be 1.
C     DATA FHZE/1.0/
      END
C-----------------------------------------------------------------------            
C End Define Constants

C-----------------------------------------------------------------------            
C GNEFF 
C These data are from Neff et al. 1984 (K. Cahoy note)
C Geometric albedos of Uranus, Neptune, and Titan
C Not sure what the fitting stuff at the end is exactly yet though
C Probably don't really need this in the output at all.
C Suspect was used to make a comparative plot in Marley et al. 1999
C Also see "GEOMAL" (not to be confused with GEOMALB, GEOMALB2)
C for use of these numbers

      SUBROUTINE GNEFF(WNEFF,XNEFF,NN)
      DIMENSION NEFF(702),XNEFF(702),WNEFF(702)
      DATA (NEFF(I),I=1,240)/ 
     & 69, 65, 58, 83, 55, 80, 74, 64, 56, 61, 92, 76, 66, 60, 77, 71,
     & 73, 66, 75, 76, 73, 67, 76, 72, 82, 72, 76, 73, 74, 74, 78, 76,
     & 80, 70, 83, 77, 79, 76, 77, 80, 80, 83, 81, 81, 81, 83, 85, 80,
     & 85, 83, 88, 85, 86, 86, 89, 88, 88, 88, 89, 90, 87, 89, 91, 91,
     & 91, 92, 93, 93, 94, 95, 95, 96, 97, 97, 97, 97, 98, 98, 99,100,
     & 99, 99,102,103,102,104,105,105,107,106,108,108,109,110,111,112,
     &111,113,113,114,115,116,117,117,119,119,120,121,122,122,123,124,
     &125,126,127,128,128,129,130,131,132,132,133,135,135,137,137,138,
     &140,140,141,141,143,143,144,145,145,145,147,148,149,150,151,152,
     &152,154,154,156,157,157,159,160,161,161,163,165,166,167,168,169,
     &170,171,174,173,174,176,177,177,181,179,181,180,183,183,184,186,
     &188,186,189,190,190,192,194,191,194,194,196,196,196,198,196,194,
     &195,195,197,201,206,208,210,212,214,214,214,215,216,217,218,219,
     &219,221,222,223,222,224,223,226,225,225,227,227,228,228,227,228,
     &228,228,227,228,230,231,235,235,237,237,238,239,239,239,241,243/
      DATA (NEFF(I),I=241,481)/   242,
     &247,247,245,246,247,246,246,248,249,250,250,254,254,256,254,256,
     &254,253,252,253,250,245,237,232,229,227,224,217,215,215,220,225,
     &234,243,251,260,263,273,271,275,272,275,274,275,274,276,276,275,
     &273,274,274,273,273,275,275,273,275,275,272,270,270,274,273,269,
     &268,266,261,262,263,263,263,260,259,259,256,254,253,255,256,261,
     &265,271,277,282,287,288,291,292,291,284,282,278,275,276,279,279,
     &315,284,293,294,291,290,292,290,284,279,268,266,260,250,245,238,
     &240,241,246,253,257,262,265,268,274,275,275,275,270,252,236,222,
     &215,206,198,195,191,187,176,170,166,169,177,179,190,196,197,206,
     &215,224,228,240,249,261,278,284,288,295,296,297,294,298,299,294,
     &299,298,301,301,300,298,294,291,282,361,351,274,306,289,274,273,
     &256,262,260,259,257,256,241,242,244,237,233,229,224,215,207,205,
     &201,198,211,206,208,213,217,214,212,208,200,201,195,198,193,192,
     &194,198,202,203,207,205,211,215,223,227,225,230,237,251,259,274,
     &298,304,301,314,314,316,318,324,339,323,330,334,332,335,342,337/
      DATA (NEFF(I),I=482,702)/ 
     &332,332,323,301,265,249,235,217,209,220,202,208,222,219,232,241,
     &246,248,256,243,243,233,217,211,201,193,178,173,168,160,146,146,
     &148,145,147,157,155,162,171,173,178,181,189,186,189,184,184,178,
     &164,142,123,100, 93,102,106, 99, 82, 88, 94, 94, 93, 97, 94, 88,
     & 91, 93,104,108,106,111,118,128,134,143,149,163,169,172,187,185,
     &179,176,173,180,177,181,186,199,194,207,213,216,216,223,234,245,
     &232,256,261,271,269,286,275,273,289,309,305,299,310,315,313,307,
     &286,262,251,241,240,214,204,211,221,199,181,217,213,199,195,186,
     &188,176,173,163,155,137,131,124,115,122,116,118,109,107,107,115,
     &112,114,113,113,116,106,110,106,106, 95, 87, 99, 86, 93, 72, 72,
     & 75, 71, 82, 81, 84, 71, 79, 82, 75, 75, 84, 88, 82, 83, 72, 72,
     & 72, 73, 56, 58, 59, 71, 54, 65, 67, 68, 90, 74, 70, 75, 87, 75,
     & 95, 51, 79, 66,109,102, 49, 88, 72, 88, 74, 39, 92, 55,104, 98,
     & 71, 87, 97,105, 45,116,188, 99, 99,105,112, 59,162/
       DO I=1,702
       XNEFF(I)=NEFF(I)* (2560./2850.)**2 *0.001
       ENDDO
       WNEFF(1)=3500.E-4
       DO 2 I=2,702,1
       WNEFF(I)=WNEFF(I-1)+10.E-4
 2     CONTINUE
       NN=702
       RETURN
       END     

C-----------------------------------------------------------------------            
C End GNEFF
C-----------------------------------------------------------------------
C GEOMALB3 (Not to be confused with GEOMAL, which scales wvln to planet)
C This one takes different planetary phase angles (see Horak et al 1965)
C GEOMALB only does face-on scenario

       SUBROUTINE GEOMALB3(NLEVEL,DTAUV,TAUV,WBARV,COSBV,gcos2,TGEOM,TESTXINT,SLABGEO,SLABINT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       include 'pangle_params'    
   
       COMMON /PLANT/ CSUBP,RSFI,RSFV,F0PI
       COMMON /OBSPHA/ PLPHA, PANGLE(NPHAS)
       COMMON /UMTXS/ U0MTX,U1MTX,NUMTX,DPMTX 
       DIMENSION DTAUV(1),TAUV(1),WBARV(1),COSBV(1),gcos2(1)
       DIMENSION GWGHT(8),GANGL(8)
       DIMENSION GWGHT3(100),GANGL3(100),TWGHT3(10),TANGL3(10)
       REAL U0MTX(100,10),U1MTX(100,10),NUMTX(100,10),DPMTX(100,10)
       DIMENSION TESTXINT(1000)
       
       DIMENSION SLABGEO(100,10),SLABINT(100,10)
       
C Variables used in this subroutine
       REAL NUTEMP

C THIS SUBROUTINE RETURNS THE GEOMETRIC ALBEDO GIVEN THE LAYER OPACITIES
C FOR EACH SINGLE WAVELENGTH. THE PLANET DISK IS DIVIDED INTO GAUSSIAN
C AND TSCHEBICHEV ANGLES AND WEIGHTS FOR PERFORMING THE INTENSITY
C INTEGRATION AS A FUNCTION OF PLANETARY PHASE ANGLE PLPHA. PLPHA AND
C THE GANGL2 AND TANGL2 ARE USED TO CALCULATE THE UBAR0 AND UBAR1 THAT
C ARE PASSED INTO GINT3. THE GANGL2 AND TANGL2'S ARE FOR N=6 OR 36 
C INTENSITY POINTS ON THE PLANET (18 X 2). USEFUL REFERENCE ON THE
C GAUSSIAN-TSCHEBICHEV QUADRATURE IS CHANDRASEKHAR'S RAD. TRANS. BOOK
C ALSO WIKIPEDIA HAS A GOOD SUMMARY.

C WE COMPUTE THE INTENSITY FROM EACH PLANE PARALLEL SLAB WITH THE 
C TWO STREAM SOURCE METHOD OF TOON ET AL. (1992) (SEE APPENDIX)
C THIS IS EFFECTIVELY EDDINGTION'S SECOND APPROXIMATION 
 
C NLEVEL           IS THE TOTAL NUMBER OF LEVELS, = LAYERS+1
C DTAUV(NLAYER)    IS THE LAYER OPTICAL THICKNESS
C TAUV(NLEVEL)     IS THE COLUMN OPTICAL THICKNESS AT LEVEL
C WBARV(NLAYER)    IS THE LAYER SINGLE SCATTERING ALBEDO
C COSBV(NLAYER)    IS THE LAYER VALUE FOR COSBAR (ASYMMETRY FACTOR)
C TGEOM            IS THE RETURNED GEOMETRIC ALBEDO


      DATA GANGL3/0.0156289844215431, 0.0468716824215916, 0.0780685828134366,
     &            0.109189203580061, 0.140203137236114, 0.171080080538603,
     &            0.201789864095736, 0.232302481844974, 0.262588120371503,    
     &   		  0.292617188038472, 0.322360343900529, 0.351788526372422,
     & 		      0.38087298162463, 0.409585291678302, 0.437897402172031,
     &    		  0.465781649773358, 0.493210789208191, 0.520158019881763,
     &     		  0.546597012065094, 0.572501932621381, 0.597847470247179,
     &   		  0.622608860203708, 0.646761908514129, 0.670283015603141,
     &   		  0.693149199355802, 0.715338117573056, 0.736828089802021,
     &   		  0.757598118519707, 0.777627909649495, 0.796897892390314,
     &   		  0.815389238339176, 0.833083879888401, 0.849964527879591,
     &   		  0.866014688497165, 0.881218679385018, 0.895561644970727,
     &  		  0.90902957098253, 0.921609298145334, 0.93328853504308,
     &  		  0.944055870136256, 0.953900782925492, 0.962813654255816,
     & 		      0.970785775763706, 0.977809358486918, 0.983877540706057,
     &  		  0.988984395242992, 0.993124937037443, 0.996295134733125,
     & 		      0.998491950639596, 0.999713726773441, -0.999713726773441,
     &		      -0.998491950639596, -0.996295134733125, -0.993124937037443,
     & 		      -0.988984395242992, -0.983877540706057, -0.977809358486918,
     & 		      -0.970785775763706, -0.962813654255816, -0.953900782925492,
     & 		      -0.944055870136256, -0.93328853504308, -0.921609298145334,
     & 		      -0.90902957098253, -0.895561644970727, -0.881218679385018,
     & 		      -0.866014688497165, -0.849964527879591, -0.833083879888401,
     & 		      -0.815389238339176, -0.796897892390314, -0.777627909649495,
     & 		      -0.757598118519707, -0.736828089802021, -0.715338117573056,
     & 		      -0.693149199355802, -0.670283015603141, -0.646761908514129,
     & 		      -0.622608860203708, -0.597847470247179, -0.572501932621381,
     & 		   	  -0.546597012065094, -0.520158019881763, -0.493210789208191,
     & 		   	  -0.465781649773358, -0.437897402172031, -0.409585291678302,
     &		      -0.38087298162463, -0.351788526372422, -0.322360343900529,
     &		      -0.292617188038472, -0.262588120371503, -0.232302481844974,
     & 			  -0.201789864095736, -0.171080080538603, -0.140203137236114,
     &		      -0.109189203580061, -0.0780685828134366, -0.0468716824215916,
     & 			  -0.0156289844215431/
     
      DATA GWGHT3/0.0312554234538633, 0.0312248842548494, 0.0311638356962099,
     &            0.0310723374275665, 0.0309504788504911, 0.0307983790311526,
     &            0.0306161865839805, 0.0304040795264549, 0.0301622651051691,   
     &   		  0.0298909795933329, 0.0295904880599126, 0.0292610841106383,
     & 		      0.0289030896011252, 0.0285168543223951, 0.0281027556591013,
     &    		  0.0276611982207923, 0.0271926134465768, 0.0266974591835711,
     &     		  0.0261762192395458, 0.0256294029102083, 0.0250575444815794,
     &   		  0.0244612027079572, 0.023840960265968, 0.0231974231852544,
     &   		  0.0225312202563361, 0.0218430024162473, 0.0211334421125274,
     &   		  0.0204032326462096, 0.0196530874944353, 0.0188837396133747,
     &   		  0.0180959407221285, 0.0172904605683237, 0.0164680861761449,
     &   		  0.0156296210775468, 0.0147758845274413, 0.0139077107037188,
     &  		  0.0130259478929719, 0.0121314576629789, 0.0112251140231866,
     &  		  0.010307802574868, 0.00938041965369416, 0.00844387146966811,
     & 		      0.00749907325546526, 0.00654694845084468, 0.00558842800386471,
     &  		  0.00462445006342216, 0.00365596120132885, 0.00268392537155496,
     & 		      0.00170939265351455, 0.000734634490495073, 0.000734634490495073,
     &		      0.00170939265351455, 0.00268392537155496, 0.00365596120132885,
     & 		      0.00462445006342216, 0.00558842800386471, 0.00654694845084468,
     & 		      0.00749907325546526, 0.00844387146966811, 0.00938041965369416,
     & 		      0.010307802574868, 0.0112251140231866, 0.0121314576629789,
     & 		      0.0130259478929719, 0.0139077107037188, 0.0147758845274413,
     & 		      0.0156296210775468, 0.0164680861761449, 0.0172904605683237,
     & 		      0.0180959407221285, 0.0188837396133747, 0.0196530874944353,
     & 		      0.0204032326462096, 0.0211334421125274, 0.0218430024162473,
     & 		      0.0225312202563361, 0.0231974231852544, 0.023840960265968,
     & 		      0.0244612027079572, 0.0250575444815794, 0.0256294029102083,
     & 		   	  0.0261762192395458, 0.0266974591835711, 0.0271926134465768,
     & 		   	  0.0276611982207923, 0.0281027556591013, 0.0285168543223951,
     &		      0.0289030896011252, 0.0292610841106383, 0.0295904880599126,
     &		      0.0298909795933329, 0.0301622651051691, 0.0304040795264549,
     & 			  0.0306161865839805, 0.0307983790311526, 0.0309504788504911,
     &		      0.0310723374275666, 0.0311638356962099, 0.0312248842548494,
     & 			  0.0312554234538633/
         


      DATA TANGL3/0.1423148382, 0.4154150130, 0.6548607339,
     &            0.8412535328, 0.9594929736, -0.9594929736,
     &            -0.8412535328, -0.6548607339, -0.4154150130,
     &            -0.1423148382/

      DATA TWGHT3/0.2798149423, 0.2363135602, 0.1631221774,
     &            0.0834785409, 0.0226689425, 0.0226689425,
     &            0.0834785409, 0.1631221774, 0.2363135602,
     &            0.2798149423/

      DATA PI/3.1415926/
      F0PI=1. ! WE ONLY SEEK A NORMALIZED ALBEDO
      BTOP=0.0 ! NO LIGHT EMITTED FROM THE SURFACE
C      write(*,*)'Number of Phase Divisions:',NPHAS
C      write(*,*)'Check:',SIZE(PANGLE)

           DO NG=1,100
             DO NT=1,10
               
               NUMTX(NG,NT)=0.0
               U0MTX(NG,NT)=0.0
               U1MTX(NG,NT)=0.0 
               DPMTX(NG,NT)=0.0
               
             ENDDO ! on NT               
           ENDDO ! on NG

           TGEOM=0.0
      
           DO NG=1,100
             DO NT=1,10
               
c              print *,'NG = ', NG
c              print *,'NT = ', NT
c               print *,'PLPHA = ', PLPHA
               NUMTX(NG,NT)=(GANGL3(NG)-(COS(PLPHA)-1)/(COS(PLPHA)+1))/(2/(COS(PLPHA)+1))
               NUTEMP = NUMTX(NG,NT)
c               print *,'NUTEMP = ', NUTEMP
               NUTEMP = ASIN(NUMTX(NG,NT))
c               print *,'ASINNUTEMP = ', NUTEMP          
               U0MTX(NG,NT)=SIN(ACOS(TANGL3(NT)))*COS(NUTEMP-PLPHA)
c               print *, 'U0 = ', U0MTX(NG,NT)
               U1MTX(NG,NT)=SIN(ACOS(TANGL3(NT)))*COS(NUTEMP) 
c               print *, 'U1 = ', U1MTX(NG,NT)
               DPMTX(NG,NT)=ACOS((U0MTX(NG,NT)*U1MTX(NG,NT)-COS(PLPHA))/(SIN(ACOS(U0MTX(NG,NT)))*SIN(ACOS(U1MTX(NG,NT)))+0.0001))  
               
             ENDDO ! on NT               
           ENDDO ! on NG

           TGEOM=0.0
           NUMGT=0           

           DO NG=1,100
             DO NT=1,10
           
              NUMGT=NUMGT+1
              UBAR0=U0MTX(NG,NT)
              UBAR1=U1MTX(NG,NT)

C SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE
      BSURF=0.+ RSFV*UBAR0*F0PI*EXP(-TAUV(NLEVEL)/UBAR0)

C CALL A ROUTINE THAT GETS INTENSITY XINT AT UBAR0 (FACE-ON ONLY) 	  

C THE INPUT PARAMETERS FOR GINT ARE:
C        GINT(NTODO,UBAR0,DTDEL,TDEL,WDEL
C     &     ,CDEL,F0PI,RSF,BTOP,BSURF,XINT,IPRINT)
C
C AND INPUT PARAMETERS FOR GINT2/3 ARE THE SAME EXCEPT
C INCLUDE GCOS2 (GINT2, GINT3) and UBAR1 (GINT3)

C.......................................................................  
C      CALL GINT2(NLEVEL,UBAR0,DTAUV,TAUV,WBARV
C     &   ,COSBV,gcos2,F0PI,RSFV,BTOP,BSURF,XINT1)

       CALL GINT3(NLEVEL,UBAR0,UBAR1,DTAUV,TAUV,WBARV
     &   ,COSBV,gcos2,F0PI,RSFV,BTOP,BSURF,XINT1)
C.......................................................................  
	
	
C SUM OVER THE GAUSS POINTS
c           print *,'XINT1 = ', xint1
c           print *,'NUMGT = ', NUMGT
           TESTXINT(NUMGT)=xint1
C           print *, 'TESTXINT =', TESTXINT(NUMGT)

C-----------------------------------------------------------------------
           SLABGEO(NG,NT)=TWGHT3(NT)*GWGHT3(NG)*XINT1
           SLABGEO(NG,NT)=SLABGEO(NG,NT)/F0PI
           SLABGEO(NG,NT)=0.5*SLABGEO(NG,NT)*(COS(PLPHA)+1.0)
           
c           if (SLABGEO(NG,NT).GT.1) then
c           SLABGEO(NG,NT)=1
c           else if (SLABGEO(NG,NT).LT.0) then
c           SLABGEO(NG,NT)=0
c           end if
           
c           SLABINT(NG,NT)=XINT1
c           if (ISNAN(SLABINT(NG,NT))) then
c           SLABINT(NG,NT)=0
c           end if
           
c
C-----------------------------------------------------------------------

C Write albedos to file
           TGEOM=TGEOM+TWGHT3(NT)*GWGHT3(NG)*XINT1
c           print *,'TGEOM = ', TGEOM
c		   print *,'NG = ', NG
c           print *,'NT = ', NT	                                    
             ENDDO ! on NT
           ENDDO ! on NG           
           
C DIVIDE BY THE LAMBERT DISK RATIO (FOR ALL ALPHA HERE, UNLIKE GEOMALB)
           TGEOM=TGEOM/F0PI
C MULTIPLY BY PLANETARY-PHASE ANGLE SCALING FACTOR (FROM SETTING UP
C THE NUMERICAL INTEGRATION QUADRATURE)
           TGEOM=0.5*TGEOM*(COS(PLPHA)+1.0)		   
		   
        END 

C-----------------------------------------------------------------------            
C END GEOMALB3

C-----------------------------------------------------------------------             
C GEOMALB2 (Not to be confused with GEOMAL, which scales wvln to planet)
C This one takes different planetary phase angles (see Horak et al 1965)
C GEOMALB only does face-on scenario

       SUBROUTINE GEOMALB2(NLEVEL,DTAUV,TAUV,WBARV,COSBV,gcos2,TGEOM,TESTXINT,SLABGEO,SLABINT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       include 'pangle_params'

       COMMON /PLANT/ CSUBP,RSFI,RSFV,F0PI
       COMMON /OBSPHA/ PLPHA, PANGLE(NPHAS)
       COMMON /UMTXS/ U0MTX,U1MTX,NUMTX,DPMTX 
       DIMENSION DTAUV(1),TAUV(1),WBARV(1),COSBV(1),gcos2(1)
       DIMENSION GWGHT(8),GANGL(8)
       DIMENSION GWGHT2(10),GANGL2(10),TWGHT2(10),TANGL2(10)
       REAL U0MTX(10,10),U1MTX(10,10),NUMTX(10,10),DPMTX(10,10)
       DIMENSION TESTXINT(100)
       
       DIMENSION SLABGEO(10,10),SLABINT(10,10)
       
C Variables used in this subroutine
       REAL NUTEMP

C THIS SUBROUTINE RETURNS THE GEOMETRIC ALBEDO GIVEN THE LAYER OPACITIES
C FOR EACH SINGLE WAVELENGTH. THE PLANET DISK IS DIVIDED INTO GAUSSIAN
C AND TSCHEBICHEV ANGLES AND WEIGHTS FOR PERFORMING THE INTENSITY
C INTEGRATION AS A FUNCTION OF PLANETARY PHASE ANGLE PLPHA. PLPHA AND
C THE GANGL2 AND TANGL2 ARE USED TO CALCULATE THE UBAR0 AND UBAR1 THAT
C ARE PASSED INTO GINT3. THE GANGL2 AND TANGL2'S ARE FOR N=6 OR 36 
C INTENSITY POINTS ON THE PLANET (18 X 2). USEFUL REFERENCE ON THE
C GAUSSIAN-TSCHEBICHEV QUADRATURE IS CHANDRASEKHAR'S RAD. TRANS. BOOK
C ALSO WIKIPEDIA HAS A GOOD SUMMARY.

C WE COMPUTE THE INTENSITY FROM EACH PLANE PARALLEL SLAB WITH THE 
C TWO STREAM SOURCE METHOD OF TOON ET AL. (1992) (SEE APPENDIX)
C THIS IS EFFECTIVELY EDDINGTION'S SECOND APPROXIMATION 
 
C NLEVEL           IS THE TOTAL NUMBER OF LEVELS, = LAYERS+1
C DTAUV(NLAYER)    IS THE LAYER OPTICAL THICKNESS
C TAUV(NLEVEL)     IS THE COLUMN OPTICAL THICKNESS AT LEVEL
C WBARV(NLAYER)    IS THE LAYER SINGLE SCATTERING ALBEDO
C COSBV(NLAYER)    IS THE LAYER VALUE FOR COSBAR (ASYMMETRY FACTOR)
C TGEOM            IS THE RETURNED GEOMETRIC ALBEDO

C Following data/code notes for GEOMALB2 for planetary phase angles PLPHA
C other than 0 degrees (face on) which is handled by GEOMALB.
C See Horak et al. 1950 and 1965

C C % Tschebichev weights for n=6
C a1 = 0.426576;
C a2 = 0.274333;
C a3 = 0.0844887;
C a4 = a3;
C a5 = a2;
C a6 = a1;
C C % Tschebichev divisions for n=6
C psi1 = 0.222521;
C psi2 = 0.623490;
C psi3 = 0.900969;
C psi4 = -psi3;
C psi5 = -psi2;
C psi6 = -psi1;
C C % Gaussian weights for n=6
C b1 = 0.467914;
C b2 = 0.360762;
C b3 = 0.171324;
C b4 = b3;
C b5 = b2;
C b6 = b1;
C C % Gaussian divisions for n=6
C xi1 = 0.238619;
C xi2 = 0.661209;
C xi3 = 0.932470;
C xi4 = -xi3;
C xi5 = -xi2;
C xi6 = -xi1;


      DATA GANGL2/0.148874339, 0.4333953941, 0.6794095683,
     &            0.8650633667, 0.9739065285, -0.9739065285,
     &            -0.8650633667, -0.6794095683, -0.4333953941,
     &            -0.148874339/

      DATA GWGHT2/0.2955252247, 0.2692667193, 0.2190863625,
     &            0.1494513492, 0.0666713443, 0.0666713443,
     &            0.1494513492, 0.2190863625, 0.2692667193,
     &            0.2955252247/

      DATA TANGL2/0.1423148382, 0.4154150130, 0.6548607339,
     &            0.8412535328, 0.9594929736, -0.9594929736,
     &            -0.8412535328, -0.6548607339, -0.4154150130,
     &            -0.1423148382/

      DATA TWGHT2/0.2798149423, 0.2363135602, 0.1631221774,
     &            0.0834785409, 0.0226689425, 0.0226689425,
     &            0.0834785409, 0.1631221774, 0.2363135602,
     &            0.2798149423/


      DATA GANGL/0.0446339553, 0.1443662570, 0.2868247571, 0.4548133152, 
     &           0.6280678354, 0.7856915206, 0.9086763921, 0.9822200849/

      DATA GWGHT/0.0032951914, 0.0178429027, 0.0454393195, 0.0791995995, 
     &           0.1060473594, 0.1125057995, 0.0911190236, 0.0445508044/

      DATA PI/3.1415926/
      F0PI=1. ! WE ONLY SEEK A NORMALIZED ALBEDO
      BTOP=0.0 ! NO LIGHT EMITTED FROM THE SURFACE


           DO NG=1,10
             DO NT=1,10
               
               NUMTX(NG,NT)=0.0
               U0MTX(NG,NT)=0.0
               U1MTX(NG,NT)=0.0 
               DPMTX(NG,NT)=0.0
               
             ENDDO ! on NT               
           ENDDO ! on NG

           TGEOM=0.0
      
           DO NG=1,10
             DO NT=1,10
               
               print *,'NG = ', NG
               print *,'NT = ', NT
               print *,'PLPHA = ', PLPHA
               NUMTX(NG,NT)=(GANGL2(NG)-(COS(PLPHA)-1)/(COS(PLPHA)+1))/(2/(COS(PLPHA)+1))
               NUTEMP = NUMTX(NG,NT)
               print *,'NUTEMP = ', NUTEMP
               NUTEMP = ASIN(NUMTX(NG,NT))
               print *,'ASINNUTEMP = ', NUTEMP          
               U0MTX(NG,NT)=SIN(ACOS(TANGL2(NT)))*COS(NUTEMP-PLPHA)
               print *, 'U0 = ', U0MTX(NG,NT)
               U1MTX(NG,NT)=SIN(ACOS(TANGL2(NT)))*COS(NUTEMP) 
               print *, 'U1 = ', U1MTX(NG,NT)
               DPMTX(NG,NT)=ACOS((U0MTX(NG,NT)*U1MTX(NG,NT)-COS(PLPHA))/(SIN(ACOS(U0MTX(NG,NT)))*SIN(ACOS(U1MTX(NG,NT)))+0.0001))
			   
               
             ENDDO ! on NT               
           ENDDO ! on NG

           TGEOM=0.0
           NUMGT=0           

           DO NG=1,10
             DO NT=1,10
           
              NUMGT=NUMGT+1
              UBAR0=U0MTX(NG,NT)
              UBAR1=U1MTX(NG,NT)

C SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE
      BSURF=0.+ RSFV*UBAR0*F0PI*EXP(-TAUV(NLEVEL)/UBAR0)

C CALL A ROUTINE THAT GETS INTENSITY XINT AT UBAR0 (FACE-ON ONLY) 	  

C THE INPUT PARAMETERS FOR GINT ARE:
C        GINT(NTODO,UBAR0,DTDEL,TDEL,WDEL
C     &     ,CDEL,F0PI,RSF,BTOP,BSURF,XINT,IPRINT)
C
C AND INPUT PARAMETERS FOR GINT2/3 ARE THE SAME EXCEPT
C INCLUDE GCOS2 (GINT2, GINT3) and UBAR1 (GINT3)

C.......................................................................  
C      CALL GINT2(NLEVEL,UBAR0,DTAUV,TAUV,WBARV
C     &   ,COSBV,gcos2,F0PI,RSFV,BTOP,BSURF,XINT1)

       CALL GINT3(NLEVEL,UBAR0,UBAR1,DTAUV,TAUV,WBARV
     &   ,COSBV,gcos2,F0PI,RSFV,BTOP,BSURF,XINT1)
C.......................................................................  
	
	
C SUM OVER THE GAUSS POINTS
           print *,'XINT1 = ', xint1
           print *,'NUMGT = ', NUMGT
           TESTXINT(NUMGT)=xint1
C           print *, 'TESTXINT =', TESTXINT(NUMGT)

C-----------------------------------------------------------------------
           SLABGEO(NG,NT)=TWGHT2(NT)*GWGHT2(NG)*XINT1
           SLABGEO(NG,NT)=SLABGEO(NG,NT)/F0PI
           SLABGEO(NG,NT)=0.5*SLABGEO(NG,NT)*(COS(PLPHA)+1.0)	
           SLABINT(NG,NT)=XINT1
c
C-----------------------------------------------------------------------

C Write albedos to file
           TGEOM=TGEOM+TWGHT2(NT)*GWGHT2(NG)*XINT1
           print *,'TGEOM = ', TGEOM
				                                    
             ENDDO ! on NT
           ENDDO ! on NG           
           
C DIVIDE BY THE LAMBERT DISK RATIO (FOR ALL ALPHA HERE, UNLIKE GEOMALB)
           TGEOM=TGEOM/F0PI
C MULTIPLY BY PLANETARY-PHASE ANGLE SCALING FACTOR (FROM SETTING UP
C THE NUMERICAL INTEGRATION QUADRATURE)
           TGEOM=0.5*TGEOM*(COS(PLPHA)+1.0)		   
		   
        END 

C-----------------------------------------------------------------------            
C END GEOMALB2

C-----------------------------------------------------------------------            
C GEOMALB (Not to be confused with GEOMAL, which scales wvln to planet)

       SUBROUTINE GEOMALB(NLEVEL,DTAUV,TAUV,WBARV,COSBV,gcos2,TGEOM,TESTXINT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C       include 'prog_params'

       COMMON /PLANT/ CSUBP,RSFI,RSFV,F0PI
       DIMENSION DTAUV(1),TAUV(1),WBARV(1),COSBV(1),gcos2(1)
       DIMENSION GWGHT(8),GANGL(8)
       DIMENSION TESTXINT(36)

C THIS SUBROUTINE APPLES FOR PLANETARY PHASE = 0 only (face-on)
C THINK "BULL'S EYE" HEMISPHERE/DISK w/PLANE PARALLEL RINGS

C THIS SUBROUTINE RETURNS THE GEOMETRIC ALBEDO GIVEN THE LAYER OPACITIES
C AT A SINGLE WAVELENGTH. THE HEMISPHERE IS DIVIDED INTO EIGHT GAUSSIAN
C POINTS AND THE REFECTED LIGHT IS COMPUTED FROM EACH AND ADDED TO
C FORM THE GEOMETRIC ALBEDO. THE SPHERICAL GEOMETRY OF THE PLANET IS 
C THEREFORE TREATED AS A 16-SIDED POLYGON.

C WE COMPUTE THE INTENSITY FROM EACH PLANE PARALLEL SLAB WITH THE 
C TWO STREAM SOURCE METHOD OF TOON ET AL. (1992) (SEE APPENDIX)
C THIS IS EFFECTIVELY EDDINGTION'S SECOND APPROXIMATION 
 
C NLEVEL           IS THE TOTAL NUMBER OF LEVELS, = LAYERS+1
C DTAUV(NLAYER)    IS THE LAYER OPTICAL THICKNESS
C TAUV(NLEVEL)     IS THE COLUMN OPTICAL THICKNESS AT LEVEL
C WBARV(NLAYER)    IS THE LAYER SINGLE SCATTERING ALBEDO
C COSBV(NLAYER)    IS THE LAYER VALUE FOR COSBAR (ASYMMETRY FACTOR)
C TGEOM            IS THE RETURNED GEOMETRIC ALBEDO

      DATA GANGL/0.0446339553, 0.1443662570, 0.2868247571, 0.4548133152, 
     &           0.6280678354, 0.7856915206, 0.9086763921, 0.9822200849/

      DATA GWGHT/0.0032951914, 0.0178429027, 0.0454393195, 0.0791995995, 
     &           0.1060473594, 0.1125057995, 0.0911190236, 0.0445508044/

      DATA PI/3.1415926/
      F0PI=1. ! WE ONLY SEEK A NORMALIZED ALBEDO
      BTOP=0.0 ! NO LIGHT EMITTED FROM THE SURFACE
 
           TGEOM=0.0
           DO NGAUSS=1,8
           UBAR0=GANGL(NGAUSS)
		   
C.......................................................................  
C Temporary edit for testing GINT3 (K. Cahoy)
C		   UBAR1=UBAR0
C.......................................................................  

C SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE
      BSURF=0.+ RSFV*UBAR0*F0PI*EXP(-TAUV(NLEVEL)/UBAR0)

C CALL A ROUTINE THAT GETS INTENSITY XINT AT UBAR0 (FACE-ON ONLY) 	  

C THE INPUT PARAMETERS FOR GINT ARE:
C        GINT(NTODO,UBAR0,DTDEL,TDEL,WDEL
C     &     ,CDEL,F0PI,RSF,BTOP,BSURF,XINT,IPRINT)
C
C AND INPUT PARAMETERS FOR GINT2/3 ARE THE SAME EXCEPT
C INCLUDE GCOS2 (GINT2, GINT3) and UBAR1 (GINT3)

C.......................................................................  
C      CALL GINT2(NLEVEL,UBAR0,DTAUV,TAUV,WBARV
C     &   ,COSBV,gcos2,F0PI,RSFV,BTOP,BSURF,XINT1)

      CALL GINT3(NLEVEL,UBAR0,UBAR1,DTAUV,TAUV,WBARV
     &   ,COSBV,GCOS2,F0PI,RSFV,BTOP,BSURF,XINT1)
C.......................................................................  
	
	
C SUM OVER THE GAUSS POINTS
           print *,'XINT1 = ', xint1
           TESTXINT(NGAUSS)=xint1
C Write albedos to file

           TGEOM=TGEOM+GWGHT(NGAUSS)*XINT1
				   
		   ENDDO

C DIVIDE BY THE LAMBERT DISK RATIO
           TGEOM=TGEOM*2.*PI/F0PI
		   
		   
        END 

C-----------------------------------------------------------------------            
C END GEOMALB
C-----------------------------------------------------------------------            
C SETSPV
      
      SUBROUTINE SETSPV(IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      include 'prog_params'

      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      COMMON /VISGAS/ SOLARF(NSPECV)
	  
C THIS SUBROUTINE HAS THE WAVELENGTH GRID AND THE ASSOCIATED
C VALUES OF THE SOLAR FLUX AND THE GAS ABSORPTION. 
C It reads in opacities from Freedman's file (except methane in visible)
C freedman *.long from 0.350 to 30 micron
C freedman otherwise from 0.350 to 5 micron 

C      open (71,file='../freedman/wave.coverage.long') ! out to 30 microns
C      open (71,file='../freedman/wave.coverage') ! out to only 5 microns
C      open (71,file='../freedman/newwave.coverage') ! out to only 5 microns

       open (71,file='../opacities/newwave_2000.coverage') ! 0.3-1 micron at R~5000 resolution

      DO I=1,2000
        READ (71,*) dum, WLNV(I)
      ENDDO 
      close (71)  
      xnum=0.00

      data solarf/NSPECV*1./
      DO J=1,NSPECV
		WNOV(J)=1.E4/WLNV(J)
      ENDDO

      RETURN
      END
C-----------------------------------------------------------------------            
C END SETSPV
C-----------------------------------------------------------------------            
C SFLUXV	  
C Note: sfluxv calls gfluxv which uses the differential solver 
C Then gfluxv calculates solar flux at each tau layer for each k wvln
C as well as cumulative, believe the differential solution still
C uses truncated HG as it baseline for gammas and solution parameters 
      SUBROUTINE SFLUXV(IPRINT) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'prog_params'

      COMMON /VISGAS/ SOLARF(NSPECV)
      COMMON /OPTICV/ DTAUV(NLAYER,NSPECV),TAUV(NLEVEL,NSPECV)
     &                 ,WBARV(NLAYER,NSPECV), COSBV(NLAYER,NSPECV)
     &                 ,NOPTHCK(NSPECV), LAYERBOT(NSPECV)
      common /ray2/ gcos2(nlayer,nspecv)
      COMMON /SPECTV/ WNOV(NSPECV),WLNV(NSPECV)
      COMMON /FLUXV/ FNETV(NLEVEL), FUPV(NLEVEL,NSPECV),
     &  FDV(NLEVEL,NSPECV), FMNETV(NLEVEL),FMUPV(NLEVEL),FMDV(NLEVEL)
      COMMON /PLANT/ CSUBP,RSFI,RSFV,F0PI
      COMMON /CONST/RGAS,RHOP,PI,SIGSB,RTOD
      COMMON /UBARED/ UBARI,UBARV,UBAR0
      COMMON /UBARED2/ UBAR1, DPHI
      DIMENSION FUW(NLEVEL),FDW(NLEVEL)

	  write(*,*)'NLAYER=',NLAYER
	  write(*,*)'NLEVEL=',NLEVEL
	  
C ZERO THE NET FLUXES
      DO 212 J=1,NLEVEL
      FNETV(J)=-0.
      FMNETV(J)=-0.
  212 CONTINUE

C WE NOW ENTER A MAJOR LOOP OVER SPECRAL INTERVALS IN THE VISIBLE
C TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

      DO 500 K=1,NSPECV
C ZERO THE SPECTRAL FLUXES IN ANTCIPATION OF SUMMING OVER NTERMS
       DO 214 J=1,NLEVEL
       FUPV(J,K)=0.
       FDV(J,K)=0.
 214   CONTINUE

C SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE
      F0PI=SOLARF(K)
      BTOP=0.0

      BSURF=0.+ RSFV*UBAR0*F0PI*DEXP(-TAUV(NLEVEL,K)/UBAR0)

C WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C WITHIN EACH INTERVAL AT THE MIDPOINT OF THE LAYER
 
C FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO 
C RETURN FLUXES
    
      CALL GFLUXV(NLEVEL,WNOV(K),DTAUV(1,K),TAUV(1,K),WBARV(1,K)
     &   ,COSBV(1,K),F0PI,RSFV,BTOP,BSURF,FUW,FDW,FMUPV,
     &    FMDV,IPRINT)

C NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX 
      DO 300 J=1,NLEVEL
      FMNETV(J)=FMNETV(J)+( FMUPV(J)-FMDV(J) )
      FNETV(J)=FNETV(J)+( FUW(J)-FDW(J) )
      FUPV(J,K)=FUPV(J,K)+FUW(J)
      FDV(J,K)=FDV(J,K)+FDW(J)
  300 CONTINUE
  
  500 CONTINUE
      RETURN 
      END
C-----------------------------------------------------------------------            
C END SETFLUXV
C-----------------------------------------------------------------------            
C AVERGE (gotta love the six-letter Fortranisms)

      FUNCTION AVERGE(X,Y)
C THIS IS USED IN COMPUTING THE GAS MASS LAYER MIXING RATIO
C AND IN THE AEROSOL PROD AND DEN ESTIMATION.
      AVERGE = .5*SQRT(X*Y)+0.25*(X+Y)
      RETURN
      END
C-----------------------------------------------------------------------            
C END AVERGE 
C-----------------------------------------------------------------------            
C GFLUXV

      SUBROUTINE GFLUXV(NTODO,WAVEN,DTDEL,TDEL,WDEL,CDEL
     &          , F0PI,RSF,BTOP,BSURF,FP,FM,FMIDP,FMIDM,IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NL=101)
C  THE VALUE OF NL (101) MUST BE .GE. NTODO
C  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
C  FOR THE VISIBLE  FLUX AT ONE WAVELENGTH K AND SOLVES FOR THE FLUXES
C  AT LEVELS J. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
C  MEASURED FROM THE TOP OF EACH LAEYER (DTAU). TOP OF EACH LAYER HAS  
C  OPTICAL DEPTH TAU(N). LEVEL N IS ABOVE LAYER N: LAYER N
C  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
C  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.

C THIS SUBROUTINE DIFFERS FROM ITS IR CONTERPART IN THAT HERE WE SOLVE FOR 
C THE FLUXES DIRECTLY USING THE GENERALIZED NOTATION OF MEADOR AND WEAVOR
C J.A.S., 37, 630-642, 1980.

      REAL LAMDA
      REAL*8 GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,B81,B82,R81,XK1,XK2
      REAL*8 EM,EP
      COMMON /UBARED/ UBARI,UBARV,UBAR0
	  
      DIMENSION  W0(NL), COSBAR(NL),DTAU(NL),TAU(NL)
      DIMENSION  WDEL(1),CDEL(1),   DTDEL(1),TDEL(1),FM(1),FP(1)
     &    ,FMIDM(1),FMIDP(1)
      DIMENSION ALPHA(NL),LAMDA(NL),XK1(NL),XK2(NL)
      DIMENSION G1(NL),G2(NL),G3(NL)
      DIMENSION GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL)
     &,E1(NL),E2(NL),E3(NL),E4(NL),EXPTRM(NL)
      DATA PI/3.14159265358979323846/
      DATA IDELTA/1/
      NAYER=NTODO-1

C TURN ON THE DELTA-FUNCTION APPROXIMATION IF REQUIRED HERE
      IF (IDELTA .EQ. 0) THEN 
                             DO J=1,NAYER
                             W0(J)=WDEL(J)
                             COSBAR(J)=CDEL(J)
                             DTAU(J)=DTDEL(J)
                             TAU(J)=TDEL(J)
                             ENDDO
                             TAU(NTODO)=TDEL(NTODO)
      ELSE
C THIS IS FOR THE DELTA FUNCTION IF ON, BELOW     
      TAU(1)=TDEL(1)*(1.-WDEL(1)*CDEL(1)**2)
      DO J=1,NAYER
      W0(J)=WDEL(J)*(1.-CDEL(J)**2)/(1.-WDEL(J)*CDEL(J)**2)
      COSBAR(J)=CDEL(J)/(1.+CDEL(J))
      DTAU(J)=DTDEL(J)*(1.-WDEL(J)*CDEL(J)**2)
      TAU(J+1)=TAU(J)+DTAU(J)
      ENDDO
      ENDIF
	  
      DO 20 J=1,NAYER
      ALPHA(J)=SQRT( (1.-W0(J))/(1.-W0(J)*COSBAR(J)) )
C THIS SET OF G'S IS FOR THE TWO STREAM HEMI-CONSTANT
C      G1(J)=(1.-W0(J)*0.5*(1.+COSBAR(J)))/UBARV
C      G2(J)=W0(J)*0.5*(1.-COSBAR(J))/UBARV
C      G3(J)=0.5*(1.-COSBAR(J))
C THIS SET OF G'S (gammas) FOR THE "DOM" APPROACH 
      G1(J)= (SQRT(3.)*0.5)*(2. - W0(J)*(1.+COSBAR(J)))
      G2(J)= (SQRT(3.)*W0(J)*0.5)*(1.-COSBAR(J))
      G3(J)=0.5*(1.-SQRT(3.)*COSBAR(J)*UBAR0)
      LAMDA(J)=SQRT(G1(J)**2 - G2(J)**2)
      GAMA(J)=(G1(J)-LAMDA(J))/G2(J)
   20 CONTINUE
      DO 7 J=1,NAYER
          G4=1.-G3(J)
          DENOM=LAMDA(J)**2 - 1./UBAR0**2
C THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C THE SCATTERING GOES TO ZERO... 
C PREVENT THIS WITH A IF STATEMENT
C      IF ( DENOM .EQ. 0.) THEN
C         DENOM=1.E-10
       IF ( DABS(DENOM) .LT. 1.D-3) THEN
            DENOM = 1.d-3 * (denom/dabs(denom))   
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			   
               WRITE (6,99) 
   99          FORMAT (' DENOM ZERO;  RESET IN GFLUXV, W0=0?')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               END IF
          AM=F0PI*W0(J)*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )/DENOM
          AP=F0PI*W0(J)*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )/DENOM
C CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
C AT THE TOP OF THE LAYER, THAT IS, AT LOWER OPTICAL DEPTH TAU(J)
          CPM1(J)=AP*DEXP(-TAU(J)/UBAR0)
          CMM1(J)=AM*DEXP(-TAU(J)/UBAR0)
C CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C BOTTOM OF THE LAYER, THAT IS, AT HIGHER OPTICAL DEPTH TAU(J+1)
          CP(J)=AP*DEXP(-TAU(J+1)/UBAR0)
          CM(J)=AM*DEXP(-TAU(J+1)/UBAR0)
  7   CONTINUE

C NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
C WARNING IF DTAU(J) IS GREATER THAN ABOUT 35
C WE CLIP IT TO AVOID OVERFLOW.
C EXP (TAU) - EXP(-TAU) WILL BE NONSENSE. THIS IS
C CORRECTED IN THE DSOLVER ROUTINE. 
      DO J=1,NAYER
      EXPTRM(J)=35.
      IF ( LAMDA(J)*DTAU(J) .LT. 35. )  EXPTRM(J)=LAMDA(J)*DTAU(J)
      ENDDO      
C NEED TO CLIP THE EXPONENTIAL HERE.
      DO 8 J=1,NAYER
      EP=DEXP(EXPTRM(J))
      EM=1.0/EP
      E1(J)=EP+GAMA(J)*EM
      E2(J)=EP-GAMA(J)*EM
      E3(J)=GAMA(J)*EP+EM
      E4(J)=GAMA(J)*EP-EM
   8  CONTINUE

      B81=BTOP
      B82=BSURF
      R81=RSF
C CALL THE DSOLVER ROUTINE TO RETURN XK1 and XK2	  
      CALL DSOLVER(NAYER,GAMA,CP,CM,CPM1,CMM1
     &            ,E1,E2,E3,E4,B81,B82,R81,XK1,XK2)

C EVALUATE THE NTODO FLUXES THROUGH THE NAYER LAYERS 
C USE THE TOP (TAU=0) OPTICAL DEPTH EXPRESSIONS TO EVALUATE FP AND FM
C AT THE THE TOP OF EACH LAYER, J = LEVEL J
      DO 46 J=1,NAYER
           FP(J)= XK1(J) + GAMA(J)*XK2(J) + CPM1(J)
           FM(J)=GAMA(J)*XK1(J) + XK2(J) + CMM1(J)
  46  CONTINUE
C USE EXPRESSION FOR BOTTOM FLUX TO GET THE FP AND FM AT NTODO
          J=NAYER
          EP=DEXP(EXPTRM(J))
          EM=1.0/EP
          FP(J+1)=XK1(J)*EP + GAMA(J)*XK2(J)*EM + CP(J)
          FM(J+1)=XK1(J)*EP*GAMA(J) + XK2(J)*EM + CM(J)
C NOTE THAT WE HAVE SOLVED FOR THE FLUXES DIRECTLY AND NO  
C FURTHER INTEGRATION IS NEEDED, THE UBARV TERM IS ABSORBED 
C INTO THE DEFINITION OF G1 THU G4
C UBARV = .5 IS HEMISPHERIC CONSTANT
C UBARV = SQRT(3) IS GAUSS QUADRATURE
C AND OTHER CASES AS PER MEADOR AND WEAVOR, JAS, 37, 630-643,1980.

C ADD THE DIRECT FLUX TERM TO THE DOWNWELLING RADIATION, LIOU 182
          DO 80 J=1,NTODO
          FM(J)=FM(J)+UBAR0*F0PI*DEXP(-TAU(J)/UBAR0)
  80      CONTINUE

C NOW WE CALCULATE THE FLUXES AT THE MIDPOINTS OF THE LAYERS.
C EXACTLY ANALOGOUS TO THE ABOVE COMPUTATION

          DO 1982 J=1,NAYER
          EP=DEXP(0.5*EXPTRM(J))
          EM=1.0/EP
          G4=1.-G3(J)
          DENOM=LAMDA(J)**2 - 1./UBAR0**2
C THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C THE SCATTERING GOES TO ZERO...
C PREVENT THIS WITH A IF STATEMENT
C      IF ( DENOM .EQ. 0.) THEN
C          DENOM=1.E-10
       IF ( DABS(DENOM) .LT. 1.D-3) THEN
            DENOM = 1.d-3 * (denom/dabs(denom))
          
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               WRITE (6,78) 
   78          FORMAT (' DENOM ZERO;  RESET IN GFLUXV, W0=0?')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               END IF
          AM=F0PI*W0(J)*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )/DENOM
          AP=F0PI*W0(J)*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )/DENOM
C CPMID AND CMMID  ARE THE CPLUS AND CMINUS TERMS EVALUATED
C AT THE MIDDLE OF THE LAYER
          TAUMID= (TAU(J)+0.5*DTAU(J) )
          CPMID=AP*DEXP(-TAUMID/UBAR0)
          CMMID=AM*DEXP(-TAUMID/UBAR0)
          FMIDP(J)=XK1(J)*EP + GAMA(J)*XK2(J)*EM + CPMID
          FMIDM(J)=XK1(J)*EP*GAMA(J) + XK2(J)*EM + CMMID
C ADD THE DIRECT FLUX TO THE DOWNWELLING TERM
          FMIDM(J)= FMIDM(J) +UBAR0*F0PI*DEXP(-TAUMID/UBAR0)
 1982 CONTINUE
C NOW PRINTOUT IF NECESSARY
      IF (IPRINT .LT. 9) RETURN
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	  
      WRITE(6,601) F0PI,RSF,BTOP,BSURF
      DO 120 J=1,NAYER
      WRITE (6,301) TAU(J),FP(J),FM(J),DTAU(J),W0(J),COSBAR(J)
     &        ,ALPHA(J), LAMDA(J),G1(J),G2(J),G3(J)
 120  CONTINUE
      J=NTODO
      WRITE (6,301) TAU(J),FP(J),FM(J)
      WRITE (6,602)
      DO 130 J=1,NAYER
      WRITE (6,301) CP(J),CM(J),E1(J),E2(J),E3(J),E4(J),CPM1(J)
     &   ,CMM1(J),GAMA(J)
  130 CONTINUE
  301 FORMAT(1X,1P13E10.3)
  602 FORMAT(' CP(J),CM(J),E1(J),E2(J),EM1(J),EM2(J),CPM1(J),CMM1(J)'
     &    ,',GAMA(J)')
  601 FORMAT(1X,'F0PI,RSF,BTOP,BSURF= ',1P4E10.3,/
     &         ' TAU,FUP,FDOWN,DTAU(J),W0(J),COSBAR(J)'
     &    ,',ALPHA(J),LAMDA(J),G1,G2,G3')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      RETURN
      END
C-----------------------------------------------------------------------            
C END GFLUXV
C-----------------------------------------------------------------------            
C GINT

      SUBROUTINE GINT(NTODO,UBAR0,DTDEL,TDEL,WDEL
     &          ,CDEL,F0PI,RSF,BTOP,BSURF,XINT,IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NL=101)

C THE PARAMETER NL (101) MUST .GE. NLEVEL
C THIS SUBROUTINE IS BASED ON THE GFLUXV SUBROUTINE AND
C USES THE SAME COEFFICIENTS FOR THE TWO STREAM AS IN GFLUXV
C TO COMPUTE THE INTENSITY EMERGENT AT THE BACKSCATTER ANGLE
C FOR USE IN GEOMETRIC ALBEDO CALCULATIONS
C NOTE: 
C GFLUXV FOR DOWNWELLING RADIATION, WHEREAS 
C GINT INCLUDES PHASE FUNCTION FOR SCATTTERING

      REAL LAMDA
      REAL*8 GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,B81,B82,R81,XK1,XK2
      REAL*8 EM,EP
C USED IN ROUTINE
      DIMENSION  W0(NL), COSBAR(NL),DTAU(NL),TAU(NL)
C INPUT PARAMETERS	  
      DIMENSION  WDEL(1),CDEL(1),   DTDEL(1),TDEL(1),XINT(1)
	  
      DIMENSION ALPHA(NL),LAMDA(NL),XK1(NL),XK2(NL)
      DIMENSION G1(NL),G2(NL),G3(NL)
      DIMENSION GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL)
     &,E1(NL),E2(NL),E3(NL),E4(NL),EXPTRM(NL)
      DATA PI/3.14159265358979323846/
      DATA IDELTA/0/
      NAYER=NTODO-1

C TURN ON THE DELTA-FUNCTION APPROXIMATION IF REQUIRED HERE
      IF (IDELTA .EQ. 0) THEN 
                             DO J=1,NAYER
                             W0(J)=WDEL(J)
                             COSBAR(J)=CDEL(J)
                             DTAU(J)=DTDEL(J)
                             TAU(J)=TDEL(J)
                             ENDDO
                             TAU(NTODO)=TDEL(NTODO)
      ELSE
C FOR THE DELTA FUNCTION APPROXIMATION, BELOW     
      TAU(1)=TDEL(1)*(1.-WDEL(1)*CDEL(1)**2)
      DO J=1,NAYER
      W0(J)=WDEL(J)*(1.-CDEL(J)**2)/(1.-WDEL(J)*CDEL(J)**2)
      COSBAR(J)=CDEL(J)/(1.+CDEL(J))
      DTAU(J)=DTDEL(J)*(1.-WDEL(J)*CDEL(J)**2)
      TAU(J+1)=TAU(J)+DTAU(J)
      ENDDO
      ENDIF

      DO 20 J=1,NAYER
      ALPHA(J)=SQRT( (1.-W0(J))/(1.-W0(J)*COSBAR(J)) )

C THIS SET OF G'S IS FOR THE TWO STREAM HEMI-CONSTANT
C      UBARV=0.5
C      G1(J)=(1.-W0(J)*0.5*(1.+COSBAR(J)))/UBARV
C      G2(J)=W0(J)*0.5*(1.-COSBAR(J))/UBARV
C      G3(J)=0.5*(1.-COSBAR(J))

C THIS SET OF G'S IS FOR "DOM" or Gaussian Quadtrature
      G1(J)= (SQRT(3.)*0.5)*(2. - W0(J)*(1.+COSBAR(J)))
      G2(J)= (SQRT(3.)*W0(J)*0.5)*(1.-COSBAR(J))
      G3(J)=0.5*(1.-SQRT(3.)*COSBAR(J)*UBAR0)
      LAMDA(J)=SQRT(G1(J)**2 - G2(J)**2)
      GAMA(J)=(G1(J)-LAMDA(J))/G2(J)
   20 CONTINUE
      DO 7 J=1,NAYER
          G4=1.-G3(J)
          DENOM=LAMDA(J)**2 - 1./UBAR0**2

C THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C THE SCATTERING GOES TO ZERO...
C PREVENT THIS WITH A IF STATEMENT
      IF ( DENOM .EQ. 0.) THEN
               DENOM=1.E-10
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			   
               WRITE (6,99) 
   99          FORMAT (' DENOM ZERO;  RESET IN GFLUXV, W0=0?')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
               END IF
          AM=F0PI*W0(J)*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )/DENOM
          AP=F0PI*W0(J)*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )/DENOM
C CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
C AT THE TOP OF THE LAYER, THAT IS AT LOWER OPTICAL DEPTH TAU(J)
          CPM1(J)=AP*DEXP(-TAU(J)/UBAR0)
          CMM1(J)=AM*DEXP(-TAU(J)/UBAR0)
C CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C BOTTOM OF THE LAYER.  THAT IS AT HIGHER OPTICAL DEPTH TAU(J+1)
          CP(J)=AP*DEXP(-TAU(J+1)/UBAR0)
          CM(J)=AM*DEXP(-TAU(J+1)/UBAR0)
  7   CONTINUE

C NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
C WARNING: IF DTAU(J) IS GREATER THAN ABOUT 35
C EXP (TAU) - EXP(-TAU) WILL BE NONSENSE 
C THIS IS CORRECTED IN THE DSOLVER ROUTINE BUT
C STILL CLIP THE EXPONENT TO AVOID NUMERICAL ISSUES
      DO J=1,NAYER
      EXPTRM(J)=35.
      IF ( LAMDA(J)*DTAU(J) .LT. 35. ) EXPTRM(J)=LAMDA(J)*DTAU(J)
      IF ( LAMDA(J)*DTAU(J) .LT. -35. ) EXPTRM(J)=-35.
      ENDDO      

      DO 8 J=1,NAYER
      EP=DEXP(EXPTRM(J))
      EM=1.0/EP
      E1(J)=EP+GAMA(J)*EM
      E2(J)=EP-GAMA(J)*EM
      E3(J)=GAMA(J)*EP+EM
      E4(J)=GAMA(J)*EP-EM
   8  CONTINUE
   
      B81=BTOP
      B82=BSURF
      R81=RSF

C PASS IN CONSTS AND BOUNDARIES, GET OUT XK1 and XK2	  
      CALL DSOLVER(NAYER,GAMA,CP,CM,CPM1,CMM1
     &            ,E1,E2,E3,E4,B81,B82,R81,XK1,XK2)

C USE EXPRESSION FOR BOTTOM FLUX TO GET THE FP AND FM AT NTODO
          J=NAYER
          EP=DEXP(EXPTRM(J))
          EM=1.0/EP
          FPZERO=XK1(J)*EP + GAMA(J)*XK2(J)*EM + CP(J)

C NOTE THAT WE HAVE SOLVED FOR THE FLUXES DIRECTLY 
C TO GET INTENSITY FOR LAMBERT SURFACE DIVIDE BY PI
          XINT(NTODO) = FPZERO/PI 
          UBARV=1./SQRT(3.)

C REPEAT ABOVE FOR J, NOW THAT HAVE XINT(NTODO)		  
		  
          DO J=NAYER,1,-1
          DENOM=LAMDA(J)**2 - 1./UBAR0**2

C THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C THE SCATTERING GOES TO ZERO...
C PREVENT THIS WITH A IF STATEMENT
      IF ( DENOM .EQ. 0.) THEN
               DENOM=1.E-10
               WRITE (6,99) 
               END IF
          G4=1.-G3(J)
          AM=F0PI*W0(J)*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )/DENOM
          AP=F0PI*W0(J)*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )/DENOM
          AM=AM*DEXP(-TAU(J)/UBAR0)
          AP=AP*DEXP(-TAU(J)/UBAR0)
C TOONS DOM
          G=XK1(J)*(1./UBAR0 - LAMDA(J))
          H=XK2(J)*(1./UBAR0 + LAMDA(J))*GAMA(J)
          A=-AP*(G1(J)-1./UBAR0) + AM*G2(J)
C UPDATED FORM, FROM NOTES   
      G=0.5*W0(J)*XK1(J)*(
     &   (1./UBARV + 3.*COSBAR(J)*UBAR0)
     & + (1./UBARV - 3.*COSBAR(J)*UBAR0)*GAMA(J)   )
      H=0.5*W0(J)*XK2(J)*(
     &   (1./UBARV + 3.*COSBAR(J)*UBAR0)*GAMA(J)
     & + (1./UBARV - 3.*COSBAR(J)*UBAR0) )
      A= 0.5*W0(J)*(1./UBARV + 3.*COSBAR(J)*UBAR0)*AP
     & + 0.5*W0(J)*(1./UBARV - 3.*COSBAR(J)*UBAR0)*AM
C TOONS EDDINGTON
          F1=1+1.5*COSBAR(J)*UBAR0
          F2=1-1.5*COSBAR(J)*UBAR0
          G=W0(J)*XK1(J)*(F1+GAMA(J)*F2)
          H=W0(J)*XK2(J)*(GAMA(J)*F1+F2)
          A=W0(J)*(F1*AP+F2*AM)
		  
          G=G*0.5/PI
          H=H*0.5/PI
          A=A*0.5/PI

C See OPTCV SUBROUTINE, PUU0
C "PHASE FUNCTION" WHERE COSBAR (COSBV) IS RATIO OF
C sum of INDIVIDUAL ASYMs*TAUs / sum of TAUs + TAURAY
C IF TAURAY relatively big, COSBAR small, PUU0 close to 1
C IF TAURAY small, COSBAR closer to 1, PUU0 close to 0
		  
          PUU0=1-COSBAR(J)

          XINT(J)=XINT(J+1)*DEXP(-DTAU(J)/UBAR0) + (W0(J)*F0PI/(8.*PI))
     &   * PUU0* DEXP(-TAU(J)/UBAR0)*(1. - DEXP(-2.*DTAU(J)/UBAR0) )
     &   + A* (1. - DEXP(-2.*DTAU(J)/UBAR0))/2.
     &   + G*(DEXP(EXPTRM(J)-DTAU(J)/UBAR0) - 1.)/(LAMDA(J)*UBAR0 - 1.) 
     &   + H*(1.-DEXP(-EXPTRM(J)-DTAU(J)/UBAR0))/(LAMDA(J)*UBAR0 + 1.) 
         ENDDO
      RETURN
      END

C-----------------------------------------------------------------------            
C END GINT
C-----------------------------------------------------------------------            
C GINT2

      SUBROUTINE GINT2(NTODO,UBAR0,DTDEL,TDEL,WDEL
     &          ,CDEL,GCOS2DEL,F0PI,RSF,BTOP,BSURF,XINT1 )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NL=101)
      
      include 'prog_params'
	  
C THE PARAMETER NL (101) MUST .GE. NLEVEL
C THIS SUBROUTINE IS BASED ON THE GFLUXV SUBROUTINE AND
C USES THE SAME COEFFICIENTS FOR THE TWO STREAM AS IN GFLUXV
C TO COMPUTE THE INTENSITY EMERGENT AT THE BACKSCATTER ANGLE
C FOR USE IN GEOMETRIC ALBEDO CALCULATIONS
C NOTE: 
C GFLUXV FOR DOWNWELLING RADIATION, WHEREAS 
C GINT2 INCLUDES PHASE FUNCTION FOR SCATTTERING

C This version GINT2 was revised by McKay in an attempt 
C to directly include Rayleigh scattering (uses gcos2
C parameter)

      REAL LAMDA
      REAL*8 GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,B81,B82,R81,XK1,XK2
      REAL*8 EM,EP

C USED IN ROUTINE
      DIMENSION  W0(NL), COSBAR(NL),gcos2(NL),DTAU(NL),TAU(NL)
C INPUT PARAMETERS
      DIMENSION  WDEL(1),CDEL(1),gcos2del(1),DTDEL(1),TDEL(1),XINT(NL)

      DIMENSION ALPHA(NL),LAMDA(NL),XK1(NL),XK2(NL)
      DIMENSION G1(NL),G2(NL),G3(NL)
      DIMENSION GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL)
     &,E1(NL),E2(NL),E3(NL),E4(NL),EXPTRM(NL)
      DATA PI/3.14159265358979323846/
      DATA IDELTA/0/
      NAYER=NTODO-1

C TURN ON THE DELTA-FUNCTION IF REQUIRED HERE
      IF (IDELTA .EQ. 0) THEN 
                             DO J=1,NAYER
                             W0(J)=WDEL(J)
                             COSBAR(J)=CDEL(J)
                             GCOS2(J)=GCOS2DEL(J)
                             DTAU(J)=DTDEL(J)
                             TAU(J)=TDEL(J)
                             ENDDO
                             TAU(NTODO)=TDEL(NTODO)
      ELSE
C FOR THE DELTA FUNCTION  HERE...     
      TAU(1)=TDEL(1)*(1.-WDEL(1)*CDEL(1)**2)
      DO J=1,NAYER
      W0(J)=WDEL(J)*(1.-CDEL(J)**2)/(1.-WDEL(J)*CDEL(J)**2)
      COSBAR(J)=CDEL(J)/(1.+CDEL(J))
      DTAU(J)=DTDEL(J)*(1.-WDEL(J)*CDEL(J)**2)
      TAU(J+1)=TAU(J)+DTAU(J)
      ENDDO
      ENDIF

      DO 20 J=1,NAYER
      ALPHA(J)=SQRT( (1.-W0(J))/(1.-W0(J)*COSBAR(J)) )

C THIS SET OF G'S IS FOR THE TWO STREAM HEMI-CONSTANT
C      UBARV=0.5
C      G1(J)=(1.-W0(J)*0.5*(1.+COSBAR(J)))/UBARV
C      G2(J)=W0(J)*0.5*(1.-COSBAR(J))/UBARV
C      G3(J)=0.5*(1.-COSBAR(J))

C SET OF CONSTANTS DETERMINED BY DOM 
      G1(J)= (SQRT(3.)*0.5)*(2. - W0(J)*(1.+COSBAR(J)))
      G2(J)= (SQRT(3.)*W0(J)*0.5)*(1.-COSBAR(J))
      G3(J)=0.5*(1.-SQRT(3.)*COSBAR(J)*UBAR0)
      LAMDA(J)=SQRT(G1(J)**2 - G2(J)**2)
      GAMA(J)=(G1(J)-LAMDA(J))/G2(J)
   20 CONTINUE
      DO 7 J=1,NAYER
          G4=1.-G3(J)
          DENOM=LAMDA(J)**2 - 1./UBAR0**2
C THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C THE SCATTERING GOES TO ZERO...
C PREVENT THIS WITH A IF STATEMENT
      IF ( DENOM .EQ. 0.) THEN
               DENOM=1.E-10
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               WRITE (6,99) 
   99          FORMAT (' DENOM ZERO;  RESET IN GFLUXV, W0=0?')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               END IF
          AM=F0PI*W0(J)*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )/DENOM
          AP=F0PI*W0(J)*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )/DENOM
C  CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
C  AT THE TOP OF THE LAYER, THAT IS LOWER   OPTICAL DEPTH TAU(J)
          CPM1(J)=AP*DEXP(-TAU(J)/UBAR0)
          CMM1(J)=AM*DEXP(-TAU(J)/UBAR0)
C CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C BOTTOM OF THE LAYER.  THAT IS AT HIGHER OPTICAL DEPTH TAU(J+1)
          CP(J)=AP*DEXP(-TAU(J+1)/UBAR0)
          CM(J)=AM*DEXP(-TAU(J+1)/UBAR0)
  7   CONTINUE

C NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
C WARNING: IF DTAU(J) IS GREATER THAN ABOUT 35
C EXP (TAU) - EXP(-TAU) WILL BE NONSENSE 
C THIS IS CORRECTED IN THE DSOLVER ROUTINE BUT
C STILL CLIP THE EXPONENT TO AVOID NUMERICAL ISSUES
      DO J=1,NAYER
      EXPTRM(J)=35.
      IF ( LAMDA(J)*DTAU(J) .LT. 35. ) EXPTRM(J)=LAMDA(J)*DTAU(J)
      IF ( LAMDA(J)*DTAU(J) .LT. -35. ) EXPTRM(J)=-35.
      ENDDO      

      DO 8 J=1,NAYER
      EP=DEXP(EXPTRM(J))
      EM=1.0/EP
      E1(J)=EP+GAMA(J)*EM
      E2(J)=EP-GAMA(J)*EM
      E3(J)=GAMA(J)*EP+EM
      E4(J)=GAMA(J)*EP-EM
   8  CONTINUE

      B81=BTOP
      B82=BSURF
      R81=RSF
      CALL DSOLVER(NAYER,GAMA,CP,CM,CPM1,CMM1
     &            ,E1,E2,E3,E4,B81,B82,R81,XK1,XK2)

C USE EXPRESSION FOR BOTTOM FLUX TO GET THE FP AND FM AT NTODO
          J=NAYER
          EP=DEXP(EXPTRM(J))
          EM=1.0/EP
          FPZERO=XK1(J)*EP + GAMA(J)*XK2(J)*EM + CP(J)

C NOTE THAT WE HAVE SOLVED FOR THE FLUXES DIRECTLY 
C TO GET INTENSITY FOR LAMBERT SURFACE DIVIDE BY PI
          XINT(NTODO) = FPZERO/PI
          UBARV=1./SQRT(3.)

C REPEAT ABOVE FOR J, NOW THAT HAVE XINT(NTODO)

          DO J=NAYER,1,-1
          DENOM=LAMDA(J)**2 - 1./UBAR0**2

		  
C THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C THE SCATTERING GOES TO ZERO...
C PREVENT THIS WITH A IF STATEMENT
      IF ( DENOM .EQ. 0.) THEN
               DENOM=1.E-10
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			   
               WRITE (6,99) 
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			   
               END IF
          G4=1.-G3(J)
          AM=F0PI*W0(J)*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )/DENOM
          AP=F0PI*W0(J)*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )/DENOM
          AM=AM*DEXP(-TAU(J)/UBAR0)
          AP=AP*DEXP(-TAU(J)/UBAR0)
C TOONS DOM
          G=XK1(J)*(1./UBAR0 - LAMDA(J))
          H=XK2(J)*(1./UBAR0 + LAMDA(J))*GAMA(J)
          A=-AP*(G1(J)-1./UBAR0) + AM*G2(J)
C McKay VALUES   
C      G=0.5*W0(J)*XK1(J)*(
C     &   (1./UBARV + 3.*COSBAR(J)*UBAR0)
C     & + (1./UBARV - 3.*COSBAR(J)*UBAR0)*GAMA(J)   )
C      H=0.5*W0(J)*XK2(J)*(
C     &   (1./UBARV + 3.*COSBAR(J)*UBAR0)*GAMA(J)
C     & + (1./UBARV - 3.*COSBAR(J)*UBAR0) )
C      A= 0.5*W0(J)*(1./UBARV + 3.*COSBAR(J)*UBAR0)*AP
C     & + 0.5*W0(J)*(1./UBARV - 3.*COSBAR(J)*UBAR0)*AM

C TOONS EDDINGTON

C TWO STREAM SOURCE FUNTION METHOD WITH G2 TERM (SEE NOTES)
C RAYLEIGH WHEN TAURAY LARGE, COSBAR (COSBV) --> 0 
C AND GCOS2 --> 0.5 (SEE OPTCV WHERE THEY ARE DEFINED)

C UBAR2 adjusted such that when TAURAY large, GEOMALB
C approaches 0.75 (3/4) at short wavelengths as expected
          UBAR2 = 0.767  ! FIT TO EXPECTED RAYLEIGH LIMIT
          F1=1.+1.5*COSBAR(J)*UBAR0
     &  + GCOS2(J)*(3.*UBAR2*UBAR2*UBAR0*UBAR0 - 1.)/2.
          F2=1.-1.5*COSBAR(J)*UBAR0 
     &  + GCOS2(J)*(3.*UBAR2*UBAR2*UBAR0*UBAR0 - 1.)/2.
          G=W0(J)*XK1(J)*(F1+GAMA(J)*F2)
          H=W0(J)*XK2(J)*(GAMA(J)*F1+F2)
          A=W0(J)*(F1*AP+F2*AM)
		  
          G=G*0.5/PI
          H=H*0.5/PI
          A=A*0.5/PI
          PUU0=(1-COSBAR(J))/(1+COSBAR(J)) + GCOS2(J)

          XINT(J)=XINT(J+1)*DEXP(-DTAU(J)/UBAR0) + (W0(J)*F0PI/(8.*PI))
     &   * PUU0* DEXP(-TAU(J)/UBAR0)*(1. - DEXP(-2.*DTAU(J)/UBAR0) )
     &   + A* (1. - DEXP(-2.*DTAU(J)/UBAR0))/2.
     &   + G*(DEXP(EXPTRM(J)-DTAU(J)/UBAR0) - 1.)/(LAMDA(J)*UBAR0 - 1.) 
     &   + H*(1.-DEXP(-EXPTRM(J)-DTAU(J)/UBAR0))/(LAMDA(J)*UBAR0 + 1.) 
         ENDDO
      xint1=xint(1)
      RETURN
      END
	  
C-----------------------------------------------------------------------            
C END GINT2
C-----------------------------------------------------------------------            
C GINT3  

      SUBROUTINE GINT3(NTODO,UBAR0,UBAR1,DTDEL,TDEL,WDEL
     &          ,CDEL,GCOS2DEL,F0PI,RSF,BTOP,BSURF,XINT1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'prog_params'
      include 'pangle_params'

      PARAMETER (NL=101)
      COMMON /OBSPHA/ PLPHA, PANGLE(NPHAS)
C THE PARAMETER NL (101) MUST .GE. NLEVEL
C THIS SUBROUTINE IS BASED ON THE GFLUXV SUBROUTINE AND
C USES THE SAME COEFFICIENTS FOR THE TWO STREAM AS IN GFLUXV
C TO COMPUTE THE INTENSITY EMERGENT AT THE BACKSCATTER ANGLE
C FOR USE IN GEOMETRIC ALBEDO CALCULATIONS
C NOTE: 
C GFLUXV FOR DOWNWELLING RADIATION, WHEREAS 
C GINT2 INCLUDES PHASE FUNCTION FOR SCATTTERING

C This version GINT3 is GINT2 as revised by McKay in an attempt 
C to directly include Rayleigh scattering (uses gcos2
C parameter) as well as revisions to include UBAR1 observer angle
C with contribs Pilorget and Cahoy

      REAL LAMDA
C Here COSTHET used to include aszimuthal asymm      
	  REAL COSTHET
      REAL*8 GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,B81,B82,R81,XK1,XK2
      REAL*8 EM,EP

C USED IN ROUTINE
      DIMENSION  W0(NL), COSBAR(NL),gcos2(NL),DTAU(NL),TAU(NL)
      DIMENSION  SW0(NL), SCOSB(NL), SDTAU(NL), STAU(NL)
C INPUT PARAMETERS	  
      DIMENSION  WDEL(1),CDEL(1),gcos2del(1),DTDEL(1),TDEL(1),XINT(NL)
      DIMENSION ALPHA(NL),LAMDA(NL),XK1(NL),XK2(NL)
      DIMENSION G1(NL),G2(NL),G3(NL)
      DIMENSION GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL)
     &,E1(NL),E2(NL),E3(NL),E4(NL),EXPTRM(NL)
      DATA PI/3.14159265358979323846/


      DATA IDELTA/1/ ! 0 is off, 1 is on
C MAKE SURE YOU SWITCH XINT TO MATCH IDELTA ON / OFF!!!       

C Used with H-G phase function and non azimuthal symm
      DATA PUU0COMPLETE/0/
      REAL DELTAPHI
      REAL DELTAPHI2
      REAL MDELT
      
C      write(*,*)'Number of Phase Divisions:',NPHAS
C      write(*,*)'Check:',SIZE(PANGLE)
            
C The +/-0.0001 is to make sure you don't acos(>1 due to num trunc)


	  DELTAPHI = ACOS((COS(PLPHA)-UBAR1*UBAR0)
     &/(SQRT(1-UBAR1**2)*SQRT(1-UBAR0**2))-0.0001)


      DELTAPHI2 = ACOS((UBAR1*UBAR0-COS(PLPHA))
     &/(SIN(ACOS(UBAR1))*SIN(ACOS(UBAR0))+0.0001))
	 
	 
C      COSTHET = UBAR0*UBAR1+SQRT(1-UBAR0**2)*SQRT(1-UBAR1**2)*COS(DELTAPHI)
      COSTHET = COS(PLPHA)
c      print *,'COSTHET =', COSTHET

c      print *,'DELTAPHI = ', DELTAPHI
c	  print *,'DELTAPHI2 =', DELTAPHI2

      
      
      NAYER=NTODO-1
C TURN ON THE DELTA-FUNCTION IF REQUIRED HERE
C NOTE THE DELTA-FUNCTION ONLY APPLIED IN CALC OF THE
C two stream source function... can't use these as-is
C for all taus here. Ok for within the two-stream but then
C when we insert the t-s F's and integrate them again in each layer, the 
C two-stream taus need to get substituted back to be in terms of the
C other ones... this isn't done here yet so please don't turn
C it on until that gets addressed.
c	  write(*,*)'NLAYER=',NLAYER

      IF (IDELTA .EQ. 0) THEN 
                             DO J=1,NAYER
                             W0(J)=WDEL(J)
                             COSBAR(J)=CDEL(J)
                             GCOS2(J)=GCOS2DEL(J)
                             DTAU(J)=DTDEL(J)
                             TAU(J)=TDEL(J)
                             SW0(J)=WDEL(J)
                             SCOSB(J)=CDEL(J)
                             SDTAU(J)=DTDEL(J)
                             STAU(J)=TDEL(J)
                             ENDDO
                             TAU(NTODO)=TDEL(NTODO)
                             STAU(NTODO)=TDEL(NTODO)
      ELSE
C FOR THE HG (f = g^2) DELTA FUNCTION HERE...     
c      print *, 'DELTA FUNCTION ON'
      TAU(1)=TDEL(1)*(1.-WDEL(1)*CDEL(1)**2)
      DO J=1,NAYER
      W0(J)=WDEL(J)*(1.-CDEL(J)**2)/(1.-WDEL(J)*CDEL(J)**2)
      COSBAR(J)=CDEL(J)/(1.+CDEL(J))
      GCOS2(J)=GCOS2DEL(J)
      DTAU(J)=DTDEL(J)*(1.-WDEL(J)*CDEL(J)**2)
      TAU(J+1)=TAU(J)+DTAU(J)
      SW0(J)=WDEL(J)
      SCOSB(J)=CDEL(J)
      SDTAU(J)=DTDEL(J)
      STAU(J)=TDEL(J)      
      ENDDO
      STAU(NTODO)=TDEL(NTODO)
      ENDIF

      DO 20 J=1,NAYER
      ALPHA(J)=SQRT( (1.-W0(J))/(1.-W0(J)*COSBAR(J)) )

C THIS SET OF G'S IS FOR THE TWO STREAM HEMI-CONSTANT
c      UBARV=0.5
c      G1(J)=(1.-W0(J)*0.5*(1.+COSBAR(J)))/UBARV
c      G2(J)=W0(J)*0.5*(1.-COSBAR(J))/UBARV
c      G3(J)=0.5*(1.-COSBAR(J))

C SET OF CONSTANTS DETERMINED BY DOM aka QUADRATURE
C See Meador and Weaver 1980, uses value 
C "mu1" in M&W = 1/sqrt(3), mu1 is constant
C NOT SAME AS MUBAR1 observer angle used here
C The "G's" correspond to lowercase gammas

C SET OF CONSTANTS DETERMINED BY DOM with 1/sqrt(3)
      G1(J)= (SQRT(3.)*0.5)*(2. - W0(J)*(1.+COSBAR(J)))
      G2(J)= (SQRT(3.)*W0(J)*0.5)*(1.-COSBAR(J))
      G3(J)=0.5*(1.-SQRT(3.)*COSBAR(J)*UBAR0)

C THIS SET OF G's with DOM/Quad + RAYLEIGH-LIKE Phase Function
C Same DOM/Quadrature solving the two-stream. 
C See Meador and Weaver 1980, uses value 
C "mu1" in M&W = 1/sqrt(3), mu1 is constant
C NOT SAME AS MUBAR1 observer angle used here
C The "G's" correspond to lowercase gammas
C SET OF CONSTANTS DETERMINED BY DOM with 1/sqrt(3)
C See K. Cahoy notes for derivation
C      G1(J)= (SQRT(3.)*0.5)*(2. - W0(J)*(1.+COSBAR(J)+GCOS2(J)/3.))
C      G2(J)= (SQRT(3.)*W0(J)*0.5)*(1.-COSBAR(J)-GCOS2(J)/3.)
C      G3(J)=0.5*(1.-SQRT(3.)*COSBAR(J)*UBAR0+GCOS2(J)*0.5*(UBAR0*UBAR0-1))	  

C Same formula for LAMDA and GAMA no matter which G's
	  	  
      LAMDA(J)=SQRT(G1(J)**2 - G2(J)**2)
      GAMA(J)=(G1(J)-LAMDA(J))/G2(J)
   20 CONTINUE
      DO 7 J=1,NAYER
          G4=1.-G3(J)
          DENOM=LAMDA(J)**2 - 1./UBAR0**2
C THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C THE SCATTERING GOES TO ZERO...
C PREVENT THIS WITH A IF STATEMENT
      IF ( DENOM .EQ. 0.) THEN
               DENOM=1.E-10
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			   
               WRITE (6,99) 
   99          FORMAT (' DENOM ZERO;  RESET IN GFLUXV, W0=0?')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               END IF
          AM=F0PI*W0(J)*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )/DENOM
          AP=F0PI*W0(J)*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )/DENOM
C  CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
C  AT THE TOP OF THE LAYER, THAT IS LOWER   OPTICAL DEPTH TAU(J)
          CPM1(J)=AP*DEXP(-TAU(J)/UBAR0)
          CMM1(J)=AM*DEXP(-TAU(J)/UBAR0)
C  CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C  BOTTOM OF THE LAYER.  THAT IS AT HIGHER OPTICAL DEPTH TAU(J+1)
          CP(J)=AP*DEXP(-TAU(J+1)/UBAR0)
          CM(J)=AM*DEXP(-TAU(J+1)/UBAR0)
  7   CONTINUE

C NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
C WARNING: IF DTAU(J) IS GREATER THAN ABOUT 35
C EXP (TAU) - EXP(-TAU) WILL BE NONSENSE 
C THIS IS CORRECTED IN THE DSOLVER ROUTINE BUT
C STILL CLIP THE EXPONENT TO AVOID NUMERICAL ISSUES
      DO J=1,NAYER
      EXPTRM(J)=35.
      IF ( LAMDA(J)*DTAU(J) .LT. 35. ) EXPTRM(J)=LAMDA(J)*DTAU(J)
      IF ( LAMDA(J)*DTAU(J) .LT. -35. ) EXPTRM(J)=-35.
      ENDDO      

      DO 8 J=1,NAYER
      EP=DEXP(EXPTRM(J))
      EM=1.0/EP
      E1(J)=EP+GAMA(J)*EM
      E2(J)=EP-GAMA(J)*EM
      E3(J)=GAMA(J)*EP+EM
      E4(J)=GAMA(J)*EP-EM
   8  CONTINUE

      B81=BTOP
      B82=BSURF
      R81=RSF
      CALL DSOLVER(NAYER,GAMA,CP,CM,CPM1,CMM1
     &            ,E1,E2,E3,E4,B81,B82,R81,XK1,XK2)

C USE EXPRESSION FOR BOTTOM FLUX TO GET THE FP AND FM AT NTODO
          J=NAYER
          EP=DEXP(EXPTRM(J))
          EM=1.0/EP
          FPZERO=XK1(J)*EP + GAMA(J)*XK2(J)*EM + CP(J)
c          print *,'FPZERO= ', FPZERO
c          print *, 'XK1 = ', XK1(J)
c          print *, 'XK2 = ', XK2(J)

C NOTE THAT WE HAVE SOLVED FOR THE FLUXES DIRECTLY 
C TO GET INTENSITY FOR LAMBERT SURFACE DIVIDE BY PI
          XINT(NTODO) = FPZERO/PI
          UBARV=1./SQRT(3.)
          DO J=NAYER,1,-1
          DENOM=LAMDA(J)**2 - 1./UBAR0**2
		  
C THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C THE SCATTERING GOES TO ZERO...
C PREVENT THIS WITH A IF STATEMENT
      IF ( DENOM .EQ. 0.) THEN
               DENOM=1.E-10
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			   
               WRITE (6,99) 
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			   
               END IF
          G4=1.-G3(J)
          AM=F0PI*W0(J)*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )/DENOM
          AP=F0PI*W0(J)*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )/DENOM
          AM=AM*DEXP(-TAU(J)/UBAR0)
          AP=AP*DEXP(-TAU(J)/UBAR0)
C TOONS DOM
C          G=XK1(J)*(1./UBAR0 - LAMDA(J))
C          H=XK2(J)*(1./UBAR0 + LAMDA(J))*GAMA(J)
C          A=-AP*(G1(J)-1./UBAR0) + AM*G2(J)

C MCKAY VALUES   
C      G=0.5*W0(J)*XK1(J)*(f
C     &   (1./UBARV + 3.*COSBAR(J)*UBAR0)
C     & + (1./UBARV - 3.*COSBAR(J)*UBAR0)*GAMA(J)   )
C      H=0.5*W0(J)*XK2(J)*(
C     &   (1./UBARV + 3.*COSBAR(J)*UBAR0)*GAMA(J)
C     & + (1./UBARV - 3.*COSBAR(J)*UBAR0) )
C      A= 0.5*W0(J)*(1./UBARV + 3.*COSBAR(J)*UBAR0)*AP
C     & + 0.5*W0(J)*(1./UBARV - 3.*COSBAR(J)*UBAR0)*AM

C TWO STREAM SOURCE FUNTION METHOD WITH G2 TERM SEE SNOOK AppE NOTES
          UBAR2 = 0.767  ! FIT TO PURE RAYLEIGH LIMIT, ~(1/sqrt(3))^(1/2)
            F1=1.+3*COSBAR(J)*UBAR1
     &  + GCOS2(J)*(3.*UBAR2*UBAR2*UBAR1*UBAR1 - 1.)/2.
          F2=1.-3*COSBAR(J)*UBAR1 
     &  + GCOS2(J)*(3.*UBAR2*UBAR2*UBAR1*UBAR1 - 1.)/2.

C TWO STREAM SOURCE FUNCTION METHOD WITH UPDATED F1, F2
C EDDINGTON WITH Rayleigh-like Phase Function
C F1, F2 see notes by K. Cahoy


          G=W0(J)*XK1(J)*(F1+GAMA(J)*F2)
          H=W0(J)*XK2(J)*(GAMA(J)*F1+F2)
          A=W0(J)*(F1*AP+F2*AM)
CC
          G=G*0.5/PI
          H=H*0.5/PI
          A=A*0.5/PI

C 2-term non-azi Rayleigh phase function option (see notes)
C		  PUU0 = 1+3*COSBAR(J)*UBAR0*UBAR1+0.5*GCOS2(J)*
C     &	  (3*COS(UBAR0*UBAR1)*COS(UBAR0*UBAR1) - 1)


C 2 term representation of Henyey-Greenstein phase functions (see Liou sec 6.5 p313)
C
C         IF (PUU0COMPLETE.EQ.0) THEN
C               PUU0=1-3*COSBAR(J)*UBAR0*UBAR1+5./4.*((COSBAR(J))**2)
C     &         *(3*(UBAR0)**2-1)*(3*(UBAR1)**2-1)
C			   print *,'PUU0 HG AZI SYMM = ', PUU0
C         ELSE
C               COSTHET=-UBAR0*UBAR1+SQRT(1-(UBAR0)**2)*SQRT(1-UBAR1**2)
C     &         *COS(DELTAPHI)
C               print *,'COSTHET = ', COSTHET
C               PUU0=(1-(CDEL(J))**2)/(1+(CDEL(J))**2-CDEL(J)*
C     &         COSTHET)**(1.5)
C               print *,'PUU0HG COMPLETE = ', PUU0
C		 ENDIF

C HG PHASE FUNCTION FOR SINGLE SCATTERING
c          PUU0=(1-SCOSB(J)**2)/SQRT((1+SCOSB(J)**2+2*SCOSB(J)*COSTHET)**3)

C Double HG to do backscatter, also called two-term hg (TTHG)
           PUU0=(1-(SCOSB(J)/2)**2)*(1-SCOSB(J)**2)
     & /SQRT((1+SCOSB(J)**2+2*SCOSB(J)*COSTHET)**3) 
     & +((SCOSB(J)/2)**2)*(1-(-SCOSB(J)/2.)**2)
     & /SQRT((1+(-SCOSB(J)/2.)**2+2*(-SCOSB(J)/2.)*COSTHET)**3)+GCOS2(J)

C This version has higher backscatter by order of magnitude (a bit more, even)
C Double HG to do backscatter, also called two-term hg (TTHG)
c           PUU0=(1-(SCOSB(J)/1.25)**2)*(1-SCOSB(J)**2)
c     & /SQRT((1+SCOSB(J)**2+2*SCOSB(J)*COSTHET)**3) 
c     & +((SCOSB(J)/1.25)**2)*(1-(-SCOSB(J)/1.25)**2)
c     & /SQRT((1+(-SCOSB(J)/1.25)**2+2*(-SCOSB(J)/1.25)*COSTHET)**3)+GCOS2(J)


C Phase function in GINT2 originally, USE THIS ONE FOR RAYLEIGH INCLUSION		  
c          PUU0=(1-COSBAR(J))/(1+COSBAR(J)) + GCOS2(J)

C XINT WITH ONLY UBAR0 (Incident and Observed same angle)
C
C          XINT(J)=XINT(J+1)*DEXP(-DTAU(J)/UBAR0) + (W0(J)*F0PI/(8.*PI))
C     &   * PUU0* DEXP(-TAU(J)/UBAR0)*(1. - DEXP(-2.*DTAU(J)/UBAR0) )
C     &   + A* (1. - DEXP(-2.*DTAU(J)/UBAR0))/2.
C     &   + G*(DEXP(EXPTRM(J)-DTAU(J)/UBAR0) - 1.)/(LAMDA(J)*UBAR0 - 1.) 
C     &   + H*(1.-DEXP(-EXPTRM(J)-DTAU(J)/UBAR0))/(LAMDA(J)*UBAR0 + 1.) 

C XINT WITH UBAR0 and UBAR1 (Incident and Observed different angles)
C Not for use with Delta-Function
C LAST WORKING VERSION (WITHOUT DELTA ON)
c		 XINT(J) = XINT(J+1)*DEXP(-SDTAU(J)/UBAR1) 
c     &   +(SW0(J)*F0PI/(4.*PI))*
c     &   (PUU0)*DEXP(-STAU(J)/UBAR0)*
c     &   (1. - DEXP(-SDTAU(J)*(UBAR0+UBAR1)/(UBAR0*UBAR1)))*
c     &   (UBAR0/(UBAR0+UBAR1))
c     &   +A*(1. - DEXP(-SDTAU(J) *(UBAR0+UBAR1)/(UBAR0*UBAR1)))*
c     &    (UBAR0/(UBAR0+UBAR1))
c     &   +G*(DEXP(EXPTRM(J)-SDTAU(J)/UBAR1) - 1.)/(LAMDA(J)*UBAR1 - 1.)
c     &   +H*(1. - DEXP(-EXPTRM(J)-SDTAU(J)/UBAR1))/(LAMDA(J)*UBAR1 + 1.)
c         ENDDO
c      xint1=xint(1)


C XINT WITH UBAR0 and UBAR1 (Incident and Observed different angles)
C Test with Delta-Function terms (for the inserted source fluxes)
         
         MDELT = (1-SW0(J)*(SCOSB(J)**2))
		 XINT(J) = XINT(J+1)*DEXP(-SDTAU(J)/UBAR1) 
     &   +(SW0(J)*F0PI/(4.*PI))*
     &   (PUU0)*DEXP(-STAU(J)/UBAR0)*
     &   (1. - DEXP(-SDTAU(J)*(UBAR0+UBAR1)/(UBAR0*UBAR1)))*
     &   (UBAR0/(UBAR0+UBAR1))
     &   +A*(1. - DEXP(-SDTAU(J) *(UBAR0+MDELT*UBAR1)/(UBAR0*UBAR1)))*
     &    (UBAR0/(UBAR0+MDELT*UBAR1))
     &   +G*(DEXP(EXPTRM(J)*MDELT-SDTAU(J)/UBAR1) - 1.)/(LAMDA(J)*MDELT*UBAR1 - 1.)
     &   +H*(1. - DEXP(-EXPTRM(J)*MDELT-SDTAU(J)/UBAR1))/(LAMDA(J)*MDELT*UBAR1 + 1.)
         ENDDO
      xint1=xint(1)


      RETURN
      END
C-----------------------------------------------------------------------            
C END GINT3
C-----------------------------------------------------------------------            
C GEOMAL (Not to be confused with GEOMALB)
C Again, this isn't very useful unless GNEFF is updated to read in 
C a relevant file and for a specific planet 
      FUNCTION GEOMAL(WLN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION XNEFF(702),X(702)
      CALL GNEFF(X,XNEFF,NN)
C FINDING THE LOWER NUMBER OF THE CORRECT INTERVAL 
C CONVERT AL TO MICRONS
       AL=WLN*1.E4
       IF (AL .LT. 3500.) THEN
          GEOMAL=XNEFF(1) 
          RETURN          
       ELSEIF (AL .GT. 10510.) THEN
          GEOMAL=XNEFF(NN)
          RETURN
       ENDIF
       NL=10*INT(AL/10.)
       J=(NL-3490 )/10 
       A=XNEFF(J)
       AU=XNEFF(J+1)
       GEOMAL=(((AL-NL)*(AU-A)/10.)+A)
       RETURN
       END
C-----------------------------------------------------------------------            
C END GEOMAL
C-----------------------------------------------------------------------            
C PHASE (Not to be confused with phase angle or PUU0 phase function)
C For phase integral only

      SUBROUTINE phase(ssa,phf,phi)
C Called in program GEOMETRIC as: call phase(wbarv(j,i),cosbv(j,i),phi(i))	  
c     Given the single scattering albedo (ssa) and the phase function
c     (phf) this subroutine will return the phase integral.
c     A negative phf will imply Rayleigh scattering.
c     The subroutine interpset must be called before phase() in order to 
c     calculate and store the second derivatives.
      DOUBLE PRECISION x1a(1:6),x2a(1:9),y2a(1:6,1:9),ya(1:6,1:9)
      DOUBLE PRECISION Ra(1:9),R2a(1:9),ssa,phf,phi,temp,temp2,temp3
      COMMON /PHASEFUNC/ x1a,x2a,y2a,ya,Ra,R2a,m,n

      if (ssa.LT.0.7) ssa=0.7
      if (phf.LT.0) then
         call splint(x2a,Ra,R2a,n,ssa,phi)
      
      else if (phf.GT.0.9) then
         temp = 0.85
         call splin2(x1a,x2a,ya,y2a,m,n,temp,ssa,phi)
         temp2 = phi
         temp = 0.90
         call splin2(x1a,x2a,ya,y2a,m,n,temp,ssa,phi)
         temp3=phi
         phi = temp3 + 20*(temp3-temp2)*(phf-temp)

      else
         call splin2(x1a,x2a,ya,y2a,m,n,phf,ssa,phi)   
      end if
      return
      END
C-----------------------------------------------------------------------            
C END PHASE (phase integral)
C-----------------------------------------------------------------------            
C INTERPSET
C Call before calling phase(), calls splie2, spline
      
      SUBROUTINE interpset()
c     Reads in table and calculates second derivatives.
c     The second derivatives are stored in the arrays y2a and R2a.
      COMMON /PHASEFUNC/ x1a,x2a,y2a,ya,Ra,R2a,m,n
      DOUBLE PRECISION x1a(1:6),x2a(1:9),y2a(1:6,1:9),ya(1:6,1:9)
      DOUBLE PRECISION Ra(1:9),R2a(1:9)
      m=6
      n=9
      open (UNIT=10,FILE='albedo.dat',STATUS='old')

      read(10,*) junk,x2a(9),x2a(8),x2a(7),x2a(6),x2a(5),
     +x2a(4),x2a(3),x2a(2),x2a(1)
      
      read(10,*) junk,Ra(9),Ra(8),Ra(7),Ra(6),Ra(5),Ra(4),
     +Ra(3),Ra(2),Ra(1)

      do i=1,m
         read(10,*) x1a(i),ya(i,9),ya(i,8),ya(i,7),ya(i,6), 
     +   ya(i,5),ya(i,4),ya(i,3),ya(i,2),ya(i,1)
      enddo
      close(10)
      call splie2(x1a,x2a,ya,m,n,y2a) 
      call spline(x2a,Ra,9,1.e30,1.e30,R2a)
      return
      END
C-----------------------------------------------------------------------            
C END INTERPSET
C-----------------------------------------------------------------------            
C SPLIE2 (used in INTERPSET, uses SPLINE)
      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      DOUBLE PRECISION x1a(m),x2a(n),y2a(m,n),ya(m,n)	
      PARAMETER (NN=100)
c	USES spline
c	Given an m by n tabulated function ya(1:m,1:n), and
c	tabulated independent variables x2a(1:n), this routine
c	constructs one-dimensional natural cubic splines or the rows
c	of ya and returns the second derivatives in the array
c	y2a(1:m,1:n). (The array x1a is included in the argument
c	list merely for consistency with routine splin2.)
      INTEGER j,k
      DOUBLE PRECISION y2tmp(NN),ytmp(NN)
      do j=1,m
	do k=1,n
	  ytmp(k)=ya(j,k)
	enddo
	call spline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
	do k=1,n
	  y2a(j,k)=y2tmp(k)
	enddo
      enddo
      return
      END
C-----------------------------------------------------------------------            
C END SPLIE2
C-----------------------------------------------------------------------            
C SPLINE (used in INTERPSET)
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=200)
c	Given arrays x(1:n) and y(1:n) containing a tabulated function,
c	i.e. y(i) = f(x(i)), with x1<x2<...<xN, and given values
c	yp1 and ypn for the first derivative of the interpolating 
c	function at points 1 and n, respectively, this routine returns
c	an array y2(1:n) of length n which contains the second
c  	derivatives of the interpolating function at the tabulated
c	points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger
c	the routine is signaled to set the corresponding boundary 
c	condition for a natural spline, with zero second derivative
c 	on that boundary.
c	Parameter:  NMAX is the largest anticipated value of n.
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
	y2(1)=0.
	u(1)=0.
      else
	y2(1)=-0.5
	u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
	sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	p=sig*y2(i-1)+2.
	y2(i)=(sig-1.)/p
	u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     +/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30) then
      	qn=0.
	un=0.
      else
	qn = 0.5
	un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
	y2(k) = y2(k)*y2(k+1) + u(k)
      enddo
      return
      END
C-----------------------------------------------------------------------            
C END SPLINE
C-----------------------------------------------------------------------            
C SPLIN2 (used in INTERPSET, uses SPLINE and SPLINT)
      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      DOUBLE PRECISION x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=200)
c	USES spline & splint
c	Given x1a,x2a,ya,m,n as described in splie2 and y2a
c	as produced by that routine; and given a desired interpolating
c	point x1,x2; this routine returns an interpolated function
c	value y by bicubic spline interpolation.
      INTEGER j,k
      DOUBLE PRECISION y2tmp(NN),ytmp(NN),yytmp(NN)
      do j=1,m
	do k=1,n
	  ytmp(k)=ya(j,k)
	  y2tmp(k)=y2a(j,k)
	enddo
	call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
      enddo
      call spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END
C-----------------------------------------------------------------------            
C END SPLIN2
C-----------------------------------------------------------------------            
C SPLINT (used in INTERPSET)
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
c	Given the arrays xa(1:n) and ya(1:n) of length n, which
c	tabulate a function (with the xa(i)'x in order), and given
c	the array y2a(1:n), which is the output from spline, and
c	given a value of x, this routine returns a cubic-spline
c	interpolated value y.
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
   1  if (khi-klo.gt.1) then
	k = (khi +klo)/2
	if(xa(k).gt.x) then
	   khi=k
	else
	   klo=k
	endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      !if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)
     +*y2a(khi))*(h**2)/6.
      return
      END
C-----------------------------------------------------------------------            
C END SPLINT
C-----------------------------------------------------------------------            
C End of egp.f file; please also see subs_egp.f
