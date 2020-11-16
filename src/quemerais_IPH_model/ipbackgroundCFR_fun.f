      SUBROUTINE BACKGROUND(fs,xpos,ypos,zpos,n_los,u1,v1,w1,fln)
c     September 2013: MAVEN version of intswan.f
c     Input file: distributions of density and emissivity (multiple scattering)
c     assumed 2.5 D (axi-symmetric)
c     
c     For Hot Model description see Thomas (1978), Lallement et al (1985)
c     
c     V0 Local Insterstellar Cloud Velocity (Solar rest frame)
c     VLON VDEC ecliptic longitude and latitude of Insterstellar wind direction
c     TEMP temperature of insterstellar cloud (LISM)
c     AMU : radiation pressure coefficient (radiation pressure divided by solar flux)
c     TDUR: H lifetime at 1 AU in seconds (1/beta) beta is ionization rate at 1 AU
c     AN: Anisotropy (here always 0)
c     DINF assumed density at infinity: preferred vakue is 0.10 cm-3 (DINF(2))
C
C     compilation seems to require the math library at least, but this
C     is handled by using gfortran to compile on the command line, which
C     automatically includes libraries
C
C     
c     
c     
c     Comments from Mike Chaffin, November 20, 2013: Eric Quemerais wrote
c     this code, and I mostly don't understand it, but it seems like
c     there are three parameters I need to change: the sun lyman alpha
c     brightness, observer position, and the look direction in ecliptic
c     coordinates. These are defined on line 161, lines 206-209, and
c     lines 234-236 of Quemerais' original code.
c     
c     I'll be converting this code to accept command line arguments for
c     these. The first should be the solar lyman alpha line center
c     intensity in photons/cm2/s/Ang, arguments 2-4 the x,y,z, position
c     of Mars, in AU, in ecliptic coordinates; arguments 5-7 the x,y,z
c     components of the look direction in ecliptic coordinates.
c
c     In fortran, all declarations must proceed all executable
c     statements, so the code is harder to read.
c
c      
c     Eric didn't define his variables, which is causing some compiler
c     problems...
      IMPLICIT NONE
      REAL, INTENT(IN) :: fs, xpos, ypos, zpos
      INTEGER, INTENT(IN) :: n_los
      REAL, INTENT(IN), DIMENSION(n_los) :: u1, v1, w1
      REAL, INTENT(OUT), DIMENSION(n_los) ::  fln
      INTEGER :: I_LOS
c--------------------------------------------------------------------------
c
      INTEGER, PARAMETER :: NR=60,NK=19
      INTEGER :: INF
      REAL :: VO, VLON, VDEC, TDUR, AN
      INTEGER :: L, K, ii
      INTEGER :: IJKL, IJKLM, IDB, IDF
      REAL :: BRANCH, XLA, PTF, AM, BOLK, PY, SPI
      REAL :: XNUZ, E2, EMAS, C, DLDN, SIGMAN, SIGMAF, DELNUD
      REAL :: ALPHE, CCC, GZERO, ALAMVENT
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: aa, bb, cc, aa1, bb1, cc1
      REAL :: x1, y1, z1, x2, y2, z2
      REAL :: u2, v2, w2
      REAL :: xot(5), xsn(5)
c--------------------------------------------------------------------------
      REAL :: RAVENT, DECVENT, SO(NR,NK,5), SN(NR,NK,5), SOT(NR,NK)
      REAL :: ALT(NR), DANS(NR,NK)
      REAL :: ANG(NK), ZALT(NR)
      REAL :: TEMP, AMU, DINF(5)
      REAL :: DPI, UA, GRAL, DTAP, DTAU, SIG
      INTEGER :: LMAX, KMAX
      COMMON/DIRVENT/ RAVENT,DECVENT
      COMMON/SOURCE/ SO, SN, SOT
      COMMON/PROF/ ALT, DANS
      COMMON/TRUMO/ ANG, ZALT
      COMMON/BITEM/ TEMP,AMU,DINF
      COMMON/TRAC/ DPI,UA,LMAX,KMAX,GRAL,DTAP,DTAU,SIG
c--------------------------------------------------------------------------
c      DIMENSION VPOS(12)
c      DATA VPOS/1993.,-18.67,32.03,-12.18,1998.,-26.05,39.44,
c     &     -25.71,2003.,-33.07,46.42,-38.79/
c     DATA VPOS/1993.,6.31,45.19,22.42,1998.,4.62,
c     &  61.49,31.16,2003.,2.97,77.25,39.61/  
      integer unit
      REAL*4 DONNE(12)
c--------------------------------------------------------------------------
      DATA DONNE/1.    ,1.21566E-05,0.4162,0.2,
     1     0.8816,1.02500E-05,0.0790,0.1,
     2     1.,5.84E-06,0.276,0.005/
C     IJKL = 1 DONNEES LYMAN-ALPHA
C     IJKL = 2 DONNEES LYMAN-BETA
C     IJKL = 3 DONNEES HELIUM

c.... LECTURE DES MODELES ................................................
      CHARACTER(LEN=*), PARAMETER :: SOURCEFN_FNAME = "/home/mike/"//
     &"Documents/Mars/3D_planetary_RT_model/src/quemerais_IPH_model/"//
     &"fsm99td12v20t80"
      OPEN(UNIT=3,FILE=SOURCEFN_FNAME, 
     1     FORM='FORMATTED',STATUS='OLD')
      READ(3,*) KMAX,LMAX,INF
      READ(3,*) VO,VLON,VDEC,TEMP,AMU,TDUR,AN,DINF(1)
c      print *,inf,lmax,kmax
c      print 2020,VO,VLON,VDEC,TEMP,AMU,TDUR,AN,DINF(1)
      READ(3,*) (ANG(L),L=1,5)
      READ(3,*) (ALT(K),(DANS(K,L),L=1,5),K=1,KMAX)
      READ(3,*) (ANG(L),L=6,10)
      READ(3,*) (ALT(K),(DANS(K,L),L=6,10),K=1,KMAX)
      READ(3,*) (ANG(L),L=11,15)
      READ(3,*) (ALT(K),(DANS(K,L),L=11,15),K=1,KMAX)
      READ(3,*) (ANG(L),L=16,19)
      READ(3,*) (ALT(K),(DANS(K,L),L=16,19),K=1,KMAX)
      READ(3,*) (ANG(L),L=1,5)
      READ(3,*) (ALT(K),(SOT(K,L),L=1,5),K=1,KMAX)
      READ(3,*) (ANG(L),L=6,10)
      READ(3,*) (ALT(K),(SOT(K,L),L=6,10),K=1,KMAX)
      READ(3,*) (ANG(L),L=11,15)
      READ(3,*) (ALT(K),(SOT(K,L),L=11,15),K=1,KMAX)
      READ(3,*) (ANG(L),L=16,19)
      READ(3,*) (ALT(K),(SOT(K,L),L=16,19),K=1,KMAX)
      READ(3,*) (ANG(L),L=1,5)
      READ(3,*) (ALT(K),(SO(K,L,1),L=1,5),K=1,KMAX)
      READ(3,*) (ANG(L),L=6,10)
      READ(3,*) (ALT(K),(SO(K,L,1),L=6,10),K=1,KMAX)
      READ(3,*) (ANG(L),L=11,15)
      READ(3,*) (ALT(K),(SO(K,L,1),L=11,15),K=1,KMAX)
      READ(3,*) (ANG(L),L=16,19)
      READ(3,*) (ALT(K),(SO(K,L,1),L=16,19),K=1,KMAX)
      READ(3,*) (ANG(L),L=1,5)
      READ(3,*) (ALT(K),(SN(K,L,1),L=1,5),K=1,KMAX)
      READ(3,*) (ANG(L),L=6,10)
      READ(3,*) (ALT(K),(SN(K,L,1),L=6,10),K=1,KMAX)
      READ(3,*) (ANG(L),L=11,15)
      READ(3,*) (ALT(K),(SN(K,L,1),L=11,15),K=1,KMAX)
      READ(3,*) (ANG(L),L=16,19)
      READ(3,*) (ALT(K),(SN(K,L,1),L=16,19),K=1,KMAX)
      DO 1212 ii=2,INF
         READ(3,*) VO,VLON,VDEC,TEMP,AMU,TDUR,AN,DINF(ii)
         READ(3,*) (ANG(L),L=1,5)
         READ(3,*) (ZALT(K),(SO(K,L,ii),L=1,5),K=1,KMAX)
         READ(3,*) (ANG(L),L=6,10)
         READ(3,*) (ZALT(K),(SO(K,L,ii),L=6,10),K=1,KMAX)
         READ(3,*) (ANG(L),L=11,15)
         READ(3,*) (ZALT(K),(SO(K,L,ii),L=11,15),K=1,KMAX)
         READ(3,*) (ANG(L),L=16,19)
         READ(3,*) (ZALT(K),(SO(K,L,ii),L=16,19),K=1,KMAX)
         READ(3,*) (ANG(L),L=1,5)
         READ(3,*) (ZALT(K),(SN(K,L,ii),L=1,5),K=1,KMAX)
         READ(3,*) (ANG(L),L=6,10)
         READ(3,*) (ZALT(K),(SN(K,L,ii),L=6,10),K=1,KMAX)
         READ(3,*) (ANG(L),L=11,15)
         READ(3,*) (ZALT(K),(SN(K,L,ii),L=11,15),K=1,KMAX)
         READ(3,*) (ANG(L),L=16,19)
         READ(3,*) (ZALT(K),(SN(K,L,ii),L=16,19),K=1,KMAX)
 1212 continue
      CLOSE(UNIT=3)
 2000 FORMAT(2F8.2,F7.2,F6.0,F6.2,G10.3,F4.1,F10.3)
 2020 FORMAT(2F8.0,F7.2,F6.0,F6.2,G10.3,F4.1,F10.3)
 6001 FORMAT(7X,5(2X,F4.0,3X))
 6501 FORMAT(7X,4(2X,F4.0,3X))
 6012 FORMAT(F8.3,5G12.4)
 6512 FORMAT(F8.3,4G12.4)
 6312 FORMAT(6G12.4)
 2815 FORMAT(F8.3,3G13.6)
C     KMAX = NB. DE VALEURS DE DISTANCES RADIALES DE LA TRAME DE CALCUL
C     (TABLEAUX ZALT-SO-SN ...)
c----------------------------------------------------------------------------
      UA=1.4959E+11
      do ii=1,INF
         dinf(ii)=dinf(ii)*1.E6
      enddo
      do ii=1,kmax
         zalt(ii)=alt(ii)*ua
         alt(ii)=alt(ii)*ua
      enddo
C---  CALCUL DE SIGMAN Section efficace pour la transition LYMAN-ALPHA
c      print *,DINF
C---  integree en frequence (cm2*hertz)
C     IJKL = 1 DONNEES LYMAN-ALPHA
C     IJKL = 2 DONNEES LYMAN-BETA
C     IJKL = 3 DONNEES HELIUM
      IJKL=1
      IJKLM=4*IJKL-3
      BRANCH=DONNE(IJKLM  )
      XLA   =DONNE(IJKLM+1)
      PTF   =DONNE(IJKLM+2)
      DTAU  =DONNE(IJKLM+3)
      AM=1.67333E-27
      BOLK=1.38046E-23
      PY=4*ATAN(1.)
      DPI=PY/180.
      SPI=SQRT(PY)
      XNUZ=1./XLA
      E2=23.0677E-20
      EMAS=9.1084E-28
      C=2.99793E+10
      DLDN=XLA**2*1.E+08/C
      SIGMAN=PY*E2*PTF/(EMAS*C)
C---  SIGMAF idem SIGMAN en cm2*Angstroem
      SIGMAF=SIGMAN*DLDN
C---  Largeur doppler
      DELNUD=XNUZ*SQRT(2.*BOLK*TEMP/AM)*1.E+2
C---- Section efficace au centre de la raie en m2
      SIG=SIGMAN/(SPI*DELNUD)
      SIG=SIG*1.E-4
      DTAP=1./SIG
      C=2.99793E+08
      ALPHE=DLDN*XNUZ*SQRT(2.*TEMP*BOLK/(AM*C**2))
      CCC=SIG/(4.*PY)*BRANCH
c     FS=3.32E+11
      GZERO=FS*SIGMAF
      GRAL=GZERO*1.E-10
c     GRAL passage de photons.s-1.m-2(1/4pi) en rayleigh 
c     -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:
C     PARAMETRES VENT
ccccccccccccccccccccccccccISM direction 
      ALAMVENT = 252.3*dpi
      DECVENT = 8.7*dpi
      A11=SIN(ALAMVENT)
      A12=-COS(ALAMVENT)
      A13=0.
      A21=COS(DECVENT)*COS(ALAMVENT)
      A22=COS(DECVENT)*SIN(ALAMVENT)
      A23=SIN(DECVENT)
      A31=-SIN(DECVENT)*COS(ALAMVENT)
      A32=-SIN(DECVENT)*SIN(ALAMVENT)
      A33=COS(DECVENT)
c     affiche coordonnees ecliptiques et equatoriales du vent.
      aa=COS(DECVENT)*COS(ALAMVENT)
      bb=COS(DECVENT)*SIN(ALAMVENT)
      cc=SIN(DECVENT)
      aa1=aa
      bb1=cos(23.45*dpi)*bb-sin(23.45*dpi)*cc
      cc1=sin(23.45*dpi)*bb+cos(23.45*dpi)*cc
      ravent=atan2(bb1,aa1)/dpi
      if (ravent.lt.0.) ravent=ravent+360.
      decvent=asin(cc1)/dpi
c      print *,'Upwind ecliptic : ',ALAMVENT/dpi,DECVENT/dpi
c      print *,'Upwind equatorial : ',ravent,dcvent
c     
      idb=1
      idf=5
      unit=16
c     do 977 ijk=idb,idf
c     write(unit,'(I4,3F6.2)') ijk,dinf(ijk)/1.E+06,AMU,tdur/1.E+06
c     977        continue

c     OPEN(UNIT=3,FILE='position',FORM='FORMATTED',STATUS='OLD')
c     read(3,*) inum
c     do ic=1,inum
c     read(3,*) numc, idate, tfsk, xpos, ypos, zpos
c     endo

ccccccposition of observer: in AU
c     xpos = cos(252.*dpi)
c     ypos = sin(252.*dpi)
c     zpos = 0.
      
c     coordonnes ecliptiques:
      x1 = xpos
      y1 = ypos
      z1 = zpos

c     coordonnes vent:
      x2 = (a11 * x1 + a12 * y1           )*ua
      y2 = (a21 * x1 + a22 * y1 + a23 * z1)*ua
      z2 = (a31 * x1 + a32 * y1 + a33 * z1)*ua

c      OPEN(UNIT=16,FILE='IPoutput_temp.txt',FORM='FORMATTED',STATUS
c     $     ='UNKNOWN')
c     write(16,2020) VO,VLON,VDEC,TEMP,AMU,TDUR,AN,DINF(2)
c     boucle en latitude
c     
c     do 779 kl = 1, 45
c     
c     de1 = 2. + 4.*float(kl-1) - 90.
c     
c     do 779 kk = 1, 90
c     
c     ra1 = 2. + 4.*float(kk-1) 
c     coordonnees dans le repere ecliptique

c     u1 = cos(ra1*dpi)*cos(de1*dpi)
c     v1 = sin(ra1*dpi)*cos(de1*dpi)
c     w1 = sin(de1*dpi)

      DO I_LOS = 1, n_los
      
      u2=a11*u1(I_LOS)+a12*v1(I_LOS)
      v2=a21*u1(I_LOS)+a22*v1(I_LOS)+a23*w1(I_LOS)
      w2=a31*u1(I_LOS)+a32*v1(I_LOS)+a33*w1(I_LOS)

c     Y axe du vent

      call intensm_ph(x2,y2,z2,u2,v2,w2,xot,xsn,idb,idf)

      IF (MODULO(I_LOS,1000) == 0) THEN
      print *, 'IPH sim done for I_LOS = ', I_LOS, ' OF ', n_los
      END IF
c     c val = 1,5 -> D = 0.05, 0.25
      
c 6784 format(8f8.2)
c 779  continue

c     The second element of the xot (optically thin) and xsn (optically
c     thick) is what we want; this case is for a specified hydrogen
c     density at infinity (in this case 0.05 cm^-3)
c      print *, xot(2), xsn(2)

      fln(I_LOS)=xsn(2)

      END DO

c     close(16)

      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FUNCTION G(TO)

      IMPLICIT NONE
      REAL :: TO, DEPI, DX, G, X, XU, U, VV, DG
      REAL :: GN, DGN, Q
      INTEGER :: KG

c     
      IF (TO.LT.0.) GO TO 29
      IF (TO.LE.2.) GO TO 25
      DEPI=2./SQRT(3.14159265358)
      DX=0.4
      G=0.
      IF (TO.LT.600.) G=DEPI*EXP(-TO)*0.5*DX
      DO 27 KG=1,10
         X=KG*DX
         XU=-X*X
         U=EXP(XU)
         VV=TO*U
         DG=0.
         IF (VV.LT.600.) DG=DEPI*U*U*EXP(-VV)
         G=G+DG*DX
 27   CONTINUE
      GO TO 30
 25   GN=1./SQRT(2.)
      DGN=GN
      Q=1.
    3 DGN=-DGN*TO*SQRT((Q+1)/(Q+2))/Q
      GN=GN+DGN
      Q=Q+1.
      IF (Q.LT.12.) GO TO 3
      G=GN
      GO TO 30
 29   G=0.
 30   RETURN
      END
c//////////////////////////////////////////////////////////////////////////
      FUNCTION T(TO)

      IMPLICIT NONE
      REAL :: TO, DEPI, DX, T, X, XU, U, UU, DT, DTN, TN, Q
      INTEGER :: KT


      IF (TO.LT.0.) GO TO 28
      IF (TO.LE.2.) GO TO 24
      DEPI=2./SQRT(3.14159265358)
      DX=0.4
      T=0.
      IF (TO.LT.600.) T=DEPI*EXP(-TO)*0.5*DX
      DO 26 KT=1,10
         X=KT*DX
         XU=-X*X
         U=EXP(XU)
         UU=TO*U
         DT=0.
         IF (UU.LT.600.) DT=DEPI*U*EXP(-UU)
         T=T+DT*DX
 26   CONTINUE
      GO TO 31
 24   TN=1.
      DTN=TN
      Q=1.
    2 DTN=-DTN*TO/SQRT(Q*(Q+1))
      TN=TN+DTN
      Q=Q+1.
      IF (Q.LT.12.)GO TO 2
      T=TN
      GO TO 31
 28   T=0.
 31   RETURN
      END
c**************************************************************************
      REAL FUNCTION TOP(XF,YF,ZF,XH,YH,ZH,imd)

      IMPLICIT NONE
      REAL, INTENT(IN) :: XF, YF, ZF, XH, YH, ZH
      INTEGER , INTENT(IN):: imd
      REAL :: XA, XB, YA, YB, ZA, ZB
      REAL :: ALTP, RA, RB, DUX, DUY, DUZ, DUR
      REAL :: XAB, YAB, ZAB, TA, DN, DN1
      REAL :: DMA, DSAB, DSA0, SAB, DT
      INTEGER :: KP
      
C     CALCUL DE L'EPAISSEUR OPTIQUE ENTRE DEUX POINTS
      INTEGER, PARAMETER :: NR=60, NK=19
      REAL*4 NORME
c--------------------------------------------------------------------------
      REAL :: SO(NR,NK,5), SN(NR,NK,5), SOT(NR,NK)
      REAL :: ALT(NR), DANS(NR,NK)
      REAL :: ANG(NK), ZALT(NR)
      REAL :: TEMP, AMU, DINF(5)
      REAL :: DPI, UA, GRAL, DTAP, DTAU, SIG
      INTEGER :: LMAX, KMAX
      COMMON/SOURCE/ SO,SN,SOT
      COMMON/PROF/ ALT,DANS
      COMMON/TRUMO/ ANG,ZALT
      COMMON/BITEM/ TEMP,AMU,DINF
      COMMON/TRAC/ DPI,UA,LMAX,KMAX,GRAL,DTAP,DTAU,SIG
c--------------------------------------------------------------------------
      XA=XF/UA
      XB=XH/UA
      YA=YF/UA
      YB=YH/UA
      ZA=ZF/UA
      ZB=ZH/UA
      altp=ALT(1)/UA
      RA=SQRT(XA*XA+YA*YA+ZA*ZA)
      RB=SQRT(XB*XB+YB*YB+ZB*ZB)
      IF (RA.LE.altp.AND.RB.LE.altp) THEN
         GOTO 50
      ENDIF
      IF (RA.GT.RB) THEN
         DUX=XA
         DUY=YA
         DUZ=ZA
         XA=XB
         YA=YB
         ZA=ZB
         XB=DUX
         YB=DUY
         ZB=DUZ
         DUR=RA
         RA=RB
         RB=DUR
      ENDIF
      XAB=XB-XA
      YAB=YB-YA
      ZAB=ZB-ZA
      NORME=SQRT(XAB*XAB+YAB*YAB+ZAB*ZAB)
      IF (NORME.LT..01) GOTO 50
      XAB=XAB/NORME
      YAB=YAB/NORME
      ZAB=ZAB/NORME
      DSA0=NORME/20.
      TA=ACOS(YA/RA)/DPI
      CALL DEN(RA*UA,TA,DN1,KP)
      DN1=DINF(imd)*DN1
      IF (KP.EQ.KMAX) KP=KP-1
      DMA=(ALT(KP+1)-ALT(KP))/3./UA
      DSAB=MIN(DMA,DSA0)
      SAB=0.
      DT=0.
 105  continue
      XA=XA+DSAB*XAB
      YA=YA+DSAB*YAB
      ZA=ZA+DSAB*ZAB
      SAB=SAB+DSAB
      RA=SQRT(XA*XA+YA*YA+ZA*ZA)
      TA=ACOS(YA/RA)/DPI
      CALL DEN(RA*UA,TA,DN,KP)
      DN=DINF(imd)*DN
      DT=DT+(DN+DN1)*.5*DSAB*SIG*UA
c     WRITE(18,*) DN,DT,DSAB,SAB,RA,KP,SIG*UA,NORME
      DN1=DN
      IF (KP.EQ.KMAX) KP=KP-1
      DMA=(ALT(KP+1)-ALT(KP))/3./UA
      DSAB=MIN(DMA,DSA0)
      IF (SAB.LE.NORME) goto 105
      TOP=DT
c     WRITE(18,*) TOP,DN,RA/UA,TA
      RETURN 
 50   TOP=0.
      RETURN 
      END
c########################################################################
      SUBROUTINE DEN(Z,T,DAN,KO)
C     CALCUL DE LA DENSITE EN UN POINT DE DISTANCE GEOCENTRIQUE Z ,T
      IMPLICIT NONE
      REAL :: Z, T, DAN
      INTEGER, PARAMETER :: NR=60, NK=19
      INTEGER :: KO, KK, J, LL, LLP, M, KKP
      REAL :: DT, DU, FL, FLP
c--------------------------------------------------------------------------
      REAL :: SO(NR,NK,5), SN(NR,NK,5), SOT(NR,NK)
      REAL :: ALT(NR), DANS(NR,NK)
      REAL :: ANG(NK), ZALT(NR)
      REAL :: TEMP, AMU, DINF(5)
      REAL :: DPI, UA, GRAL, DTAP, DTAU, SIG
      INTEGER :: LMAX, KMAX
      COMMON/SOURCE/ SO,SN,SOT
      COMMON/PROF/ ALT,DANS
      COMMON/TRUMO/ ANG,ZALT
      COMMON/BITEM/ TEMP,AMU,DINF
      COMMON/TRAC/ DPI,UA,LMAX,KMAX,GRAL,DTAP,DTAU,SIG
c--------------------------------------------------------------------------

      KO = 1
      DAN=0.
      IF (Z.LT.ALT(1)) GO TO 12
      IF (Z.GT.ALT(KMAX)) Z=ALT(KMAX)
      DO 20 J=1,LMAX
         IF (T-ANG(J)) 22,21,20
 21      LL=J
         LLP=J+1
         DT=0.
         GO TO 23
 22      LL=J-1
         LLP=J
         DT=(T-ANG(J-1))/(ANG(J)-ANG(J-1))
         GO TO 23
 20   CONTINUE
 23   DO 24 M=1,KMAX
         IF (Z-ZALT(M)) 26,25,24
 25      KK=M
         KKP=M+1
         DU=0.
         GO TO 27
 26      KK=M-1
         KKP=M
         DU=(Z-ZALT(KK))/(ZALT(KKP)-ZALT(KK))
         GO TO 27
 24   CONTINUE
 27   CONTINUE
      FL=DANS(KK,LL)+DU*(DANS(KKP,LL)-DANS(KK,LL))
      FLP=DANS(KK,LLP)+DU*(DANS(KKP,LLP)-DANS(KK,LLP))
C     INTERPOLATION ENTRE LL ET LLP
      DAN=FL+DT*(FLP-FL)
c     write(12,*) KK,DAN,Z/UA,T
      KO=KK 
 12   RETURN
      END

c**************************************************************************
      SUBROUTINE INTENSM_PH(X,Y,Z,U,V,W,FOT,FLN,idb,idf)
c     imd=1,2,3,4,5 -> dsinf=0.01,0.05,0.1,0.15,0.2
C     CALCUL DE L'INTENSITE LUMINEUSE EN  XYZ,FLN ,DUE A LA
C     DIFFUSION MULTIPLE. (FOT cas optiquement mince)
C     X,Y,Z SONT LES COORDONNEES SATELLITE (Y AXE DU VENT) (UW = 0)
C     U,V,W  DIRECTION DE VISEE .
C     TT  EPAISSEUR OPTIQUE LE LONG DE LA LIGNE DE VISEE
C     SO ET SN FONCTIONS SOURCE  PRIMAIRE ET APRES DIFFUSION MULTIPLE .
c     
      IMPLICIT NONE
      REAL :: X, Y, Z, U, V, W, FOT(5), FLN(5)
      INTEGER :: idb, idf, KO, ii, IMYN
      REAL :: DTT, DFOT, DP, DFLN, S, RR, DUA, YP, R, XAV, YAV, ZAV
      REAL :: TETA, DNA, XP, ZP, RU
      REAL :: cosff, corec, FFNN, DFLNC, TOP, T, TTTII

      INTEGER, PARAMETER :: NR=60, NK=19
      REAL :: DN(5),TT(5),FOTE(5),FN(5),FOO(5)
c--------------------------------------------------------------------------
      REAL :: SO(NR,NK,5), SN(NR,NK,5), SOT(NR,NK)
      REAL :: ALT(NR), DANS(NR,NK)
      REAL :: ANG(NK), ZALT(NR)
      REAL :: TEMP, AMU, DINF(5)
      REAL :: DPI, UA, GRAL, DTAP, DTAU, SIG
      INTEGER :: LMAX, KMAX
      COMMON/SOURCE/ SO,SN,SOT
      COMMON/PROF/ ALT,DANS
      COMMON/TRUMO/ ANG,ZALT
      COMMON/BITEM/ TEMP,AMU,DINF
      COMMON/TRAC/ DPI,UA,LMAX,KMAX,GRAL,DTAP,DTAU,SIG
c--------------------------------------------------------------------------

c     l'axe du vent est l'axe -Oy
      KO=0
      DTT=0.
      DFOT=0.
      DP=0.
      DFLN=0.
      S=0.
      do ii=1,5
         TT(ii)=0.
         FOT(ii)=0.
         FLN(ii)=0.
         DN(ii)=0.
         FOTE(ii)=0.
         FOO(ii)=0.
         FN(ii)=0.
      enddo
      RR=SQRT(X*X+Y*Y+Z*Z)
      IF (RR.GT.ALT(KMAX)) RETURN
      DUA=0.
      YP=Y
      R=RR
      XAV=X
      YAV=Y
      ZAV=Z
      IMYN=0
 52   TETA=ACOS(YP/R)/DPI
      IMYN=IMYN+1
      CALL DEN(R,TETA,DNA,KO)
      do ii=idb,idf
         DN(ii)=DINF(ii)*DNA
      enddo
      IF (DN(idb).EQ.0.) THEN
         DN(idb)=1.
      ENDIF
 62   DP=DTAP*0.05/DN(idb)
      IF (KO.LT.KMAX) THEN
         DUA=(ALT(KO+1)-ALT(KO))/2.
      ELSE
         DUA=(ALT(KMAX)-ALT(KMAX-1))/2.
      ENDIF
      DP=MIN(DP,DUA)
      DP=MAX(DP,UA/10.)
      S=S+DP
      XP=X+S*U
      YP=Y+S*V
      ZP=Z+S*W
      RU=R
      R=SQRT(XP*XP+YP*YP+ZP*ZP)
      IF (R.GT.ALT(KMAX)) GO TO 51
 55   TETA=ACOS(YP/R)/DPI
      CALL IPAL3M(R,TETA,FOTE,FOO,FN,idb,idf)
C     EPAISSEUR OPTIQUE ELEMENTAIRE SUR LE PAS DP
      DTT = TOP(XAV,YAV,ZAV,XP,YP,ZP,idb)
      cosff=(u*xp+v*yp+w*zp)/r
      corec=0.25*cosff*cosff+(11./12.)
c     corec=1.
      do ii=idb,idf
         TT(ii)=TT(ii)+DTT*dinf(ii)/dinf(idb)
         FFNN=FN(ii)+FOO(ii)*(corec-1.)
         TTTII=T(TT(ii))
         DFLNC=FFNN*GRAL*TTTII*DP
         DFLN=FN(ii)*GRAL*TTTII*DP
         FLN(ii)=FLN(ii)+DFLNC
         DFOT=FOTE(ii)*GRAL*DP
         FOT(ii)=FOT(ii)+DFOT*corec
      enddo
      XAV=XP
      YAV=YP
      ZAV=ZP
      GO TO 52
 51   CONTINUE
      do ii=idb,idf
         FOT(ii)=FOT(ii)+GRAL*DN(ii)*UA*UA/R
      enddo
      RETURN
      END
c--------------------------------------------------------------------------
      SUBROUTINE IPAL3M(R,T,CT,FOO,F,idb,idf)
C     INTERPOLATION LINEAIRE A DEUX DIMENSIONS QUI CALCULE UNE VALEUR
C     APPROCHEE DE LA FONCTION SOURCE F EN UN POINT DE LA GEOCOURONNE
C     CONNAISSANT LE TABLEAU SN(K,L,5) DES VALEURS DE LA FONCTION SOURCE
C     EN TOUT LES POINTS DE LA TRAME DE DEFINITION
C     
C     DT=ECART RELATIF A L'ANGLE D'INDICE LL
C     DU=ECART RELATIF A L'ALTITUDE D'INDICE KK
C     

      IMPLICIT NONE
      REAL :: R, T, CT(5), FOO(5), F(5)
      INTEGER :: idb, idf
      INTEGER, PARAMETER :: NR=60, NK=19
      INTEGER :: ii, J, LL, LLP, M, KK, KKP, imd
      REAL :: RZ, DT, DU, FL, FLP, COT
c--------------------------------------------------------------------------
      REAL :: SO(NR,NK,5), SN(NR,NK,5), SOT(NR,NK)
      REAL :: ALT(NR), DANS(NR,NK)
      REAL :: ANG(NK), ZALT(NR)
      REAL :: TEMP, AMU, DINF(5)
      REAL :: DPI, UA, GRAL, DTAP, DTAU, SIG
      INTEGER :: LMAX, KMAX
      COMMON/SOURCE/ SO,SN,SOT
      COMMON/PROF/ ALT,DANS
      COMMON/TRUMO/ ANG,ZALT
      COMMON/BITEM/ TEMP,AMU,DINF
      COMMON/TRAC/ DPI,UA,LMAX,KMAX,GRAL,DTAP,DTAU,SIG
c--------------------------------------------------------------------------
      do ii=1,5
         F(ii)=0.
         CT(ii)=0.
      enddo
      IF (R.LT.ALT(1).OR.R.GE.ALT(KMAX)) RETURN
      RZ=R
C     
C     BOUCLE DETERMINANT LES INDICES D'ANGLES ENCADRANT LE POINT A:LL ;LLP
C     
      DO 20 J=1,LMAX
         IF (T-ANG(J)) 22,21,20
 21      LL=J
         LLP=J+1
         DT=0.
         GO TO 23
 22      LL=J-1
         LLP=J
         DT=(T-ANG(J-1))/(ANG(J)-ANG(J-1))
         GO TO 23
 20   CONTINUE
C     
C     BOUCLE DETERMINANT LES INDICES D'ALT. ENCADRANT LE POINT A:KK ;KKP
C     
 23   DO 24 M=1,KMAX
         IF (RZ-ZALT(M)) 26,25,24
 25      KK=M
         KKP=M+1
         DU=0.
         GO TO 27
 26      KK=M-1
         KKP=M
         DU=(RZ-ZALT(KK))/(ZALT(KKP)-ZALT(KK))
         GO TO 27
 24   CONTINUE
 27   CONTINUE
C     
C     INTERPOLATION DES FONCT.SOURCES ENTRE KK ET KKP SUR LL,PUIS SUR LLP
C     
C     WRITE(8,*) KK,KKP,DU,ALT(KK)/UA,ALT(KKP)/UA
C     WRITE(8,*) LL,LLP,DT
      do imd=idb,idf
         FL=SN(KK,LL,imd)+DU*(SN(KKP,LL,imd)-SN(KK,LL,imd))
         FLP=SN(KK,LLP,imd)+DU*(SN(KKP,LLP,imd)-SN(KK,LLP,imd))
         F(imd)=FL+DT*(FLP-FL)
         FL=SO(KK,LL,imd)+DU*(SO(KKP,LL,imd)-SO(KK,LL,imd))
         FLP=SO(KK,LLP,imd)+DU*(SO(KKP,LLP,imd)-SO(KK,LLP,imd))
         FOO(imd)=FL+DT*(FLP-FL)
         FL=SOT(KK,LL)+DU*(SOT(KKP,LL)-SOT(KK,LL))
         FLP=SOT(KK,LLP)+DU*(SOT(KKP,LLP)-SOT(KK,LLP))
         COT=FL+DT*(FLP-FL)
         CT(imd)=COT*DINF(imd)
      enddo
      RETURN
      END
