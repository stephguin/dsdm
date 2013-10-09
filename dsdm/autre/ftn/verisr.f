C     VERIFICATION DE RSICBA ET DE SCHAST
C 
      SUBROUTINE VERISR(LTYP01)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      LOGICAL LTYP01
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
C 
C     Pour la verification de rsicba et schast
C 
CD    INTEGER  EPSVER , SIGVER , IUNSVE , IUNEVE
CD    INTEGER  SAUVER , SGNVER
CD    INTEGER  ADCOU  , ADSOU  , ADINT  , ADSNT
CD    INTEGER  LONEPS , LONDER , LONDEP , LONSAU
CD    INTEGER  DEPCH0 , EPSCH8
CD    INTEGER  SIGC11 , SAUCH9 , SGNC12 , FTREE7
CD    INTEGER  FTRE10 , ERRDEP , TABNIE(2) , EPSTRA , SAUTRA , DEPTRA
CD    INTEGER  DSIGVE , DEPSVE , FOTDEP  ,  NNOUVA , NDEBEP ,FOTDE1
CD    INTEGER  LONRES , COESIG  , TABNIS(2) ,I , TABNID(2) , FOTDE2
CD    INTEGER  FTECH6 , FTSCH5
C 
CD    DOUBLE PRECISION DIVIS
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      INTEGER     AM2LC      , ADM2LC
      PARAMETER (IDPROG='VERISR')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
CD    IF (NBETLC .EQ. 1 .AND. LTYP01) THEN
C 
C       SEQUENCE DE VERIFICATION DE RSICBA SCHAST
C 
CD      LONEPS = NTETA*NEPS*NGAU1
CD      LONDER = NDDL*NTETA
CD      LONDEP = NDDL*NBMAT
CD      LONSAU = NTETA*NSAU*NGAU2
C 
CD      CALL ADTBDM ('HOO-COUCHE', ADCOU)
CD      CALL ADTBDM ('SOU-COUCHE', ADSOU)
C 
CD      IF (NBINT .GT. 0) THEN
C 
CD        CALL ADTBDM ('HOO-INTERF', ADINT)
CD        CALL ADTBDM ('SOU-INTERF', ADSNT)
C 
CD      ELSE
C 
CD        ADINT = ADCOU
CD        ADSNT = ADSOU
C 
CD      END IF
C 
C       Tableaux des fonctions du temps pour les deltas admissibles
C  
CD      CALL ADTBDM ('TEMPS-SIGM', FTSCH5)
CD      CALL ADTBDM ('TEMPS-EPSI', FTECH6)
C 
C       Tableaux des deformations et contraintes admissibles TOTALES
C 
CD      CALL ADTBDM ('DEPLA-ADMI', DEPCH0)
CD      CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
CD      CALL ADTBDM ('SIG-AD-TOT', SIGC11)
C 
C       POUR LES INTERFACES
C 
CD      IF (NBINT.GT.0) THEN
C 
CD        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
CD        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
C 
CD      ELSE
C 
CD        SAUCH9 = EPSCH8
CD        SGNC12 = SIGC11
C 
CD      ENDIF
C 
CD      CALL ADTBDM ('TEMPS-REEL', FTREE7)
CD      CALL ADTBDM ('TEMP-SI-RE', FTRE10)
C 
CD      CALL POUSMD (NPICET, ERRDEP)
C 
CD      LONRES = LONEPS+LONSAU+LONDER+1
CD      CALL POUSMD (LONRES, COESIG)
CD      EPSTRA = COESIG+1
CD      SAUTRA = EPSTRA+LONEPS
CD      DEPTRA = SAUTRA+LONSAU
C 
CD      TABNIE(1) = CHARAX
CD      TABNIE(2) = LONEPS
C 
CD      TABNID(1) = CHARAX
CD      TABNID(2) = LONDER
C 
CD      CALL EXTRAD (DM(DEPCH0), 2, TABNID(1), 2, 1, DM(DEPTRA), LONDER)
CD      CALL EXTRAD (DM(EPSCH8), 2, TABNIE(1), 2, 1, DM(EPSTRA), LONEPS)
C 
CD      IF (NBINT .GT. 0) THEN
C 
CD        TABNIS(1) = CHARAX
CD        TABNIS(2) = LONSAU
C 
CD        CALL EXTRAD (DM(SAUCH9), 2, TABNIS(1), 2, 1, DM(SAUTRA), LONSAU)
C 
CD      END IF
C 
CD      CALL POUSMD (3*NPICET, FOTDE1)
CD      FOTDE2 = FOTDE1+NPICET
CD      FOTDEP = FOTDE2+NPICET
C 
CD      CALL SCASTE ('sigver', 'sgnver', DM(EPSTRA), DM(SAUTRA),
CD                    1, NPICET,
CD                    DM(FOTDE1), DM(FOTDE2))
C 
CD      CALL IMPTDT ('EVOLUTION CALCULEE DANS SCASTE ',
CD                    DM(FOTDE1), NPICET, 1)
C 
CD      CALL SCHAST (0, 'sigver', 'sgnver', 1,
CD                   DM(EPSTRA), DM(SAUTRA), 1, NPICET,
CD                   DM(FOTDEP))
C 
CD      CALL IMPTDT ('EVOLUTION CALCULEE DANS SCHAST ',
CD                    DM(FOTDEP), NPICET, 1)
C 
C       Calcul de intesp(epssol* K0 * epssol)
C 
CD      CALL NOKGLO (ADCOU, ADINT, 1, 1, DM(EPSTRA), DM(SAUTRA), DIVIS)
C 
C       Calcul de la meilleure evolution des deplacements f1(t)
C 
CD      DO I = 1, NPICET
C 
CD        DM(FOTDEP+I-1) = DM(FOTDEP+I-1 )/DIVIS
C 
CD      END DO
C 
CD      CALL IMPDT ('K0 EPSELA*EPSELA ', DIVIS)
C 
CD      CALL IMPTDT ('EVOLUTION ELSTIQUE RECALCULEE = INITIALE ?? ',
CD                    DM(FOTDEP), NPICET, 1)
C 
CD      TABNIE(1) = CHARAX
CD      TABNIE(2) = NPICET
C 
CD      CALL EXTRAD (DM(FTREE7), 2, TABNIE(1), 2, 1, DM(FOTDEP), NPICET)
C 
CD      CALL IMPTDT ('EVOLUTION ELASTIQUE INITIQLE '//IDPROG,
CD                    DM(FOTDEP), NPICET, 1)
C 
C      CD      CALL ORNEK0 (ADCOU, ADINT, 0, 1, LONEPS, LONSAU,
C      CD                   DM(EPSTRA), DM(SAUTRA), DM(EPSTRA), DM(SAUTRA),
C      CD                   DM(EPSTRA), DM(SAUTRA), NNOUVA, DM(COESIG))
C 
CD      CALL OREDK0 (ADCOU, ADINT, 0, 1,
CD                   LONEPS, LONSAU, LONDER,
CD                   DM(EPSTRA), DM(SAUTRA), DM(DEPTRA),
CD                   DM(EPSTRA), DM(SAUTRA), DM(DEPTRA),
CD                   DM(EPSTRA), DM(SAUTRA), DM(DEPTRA),
CD                   NNOUVA, DM(COESIG))
C 
CD      CALL SCHAST (0, 'sigver', 'sgnver',  1,
CD                   DM(EPSTRA), DM( SAUTRA), 1, NPICET,
CD                   DM(FOTDEP))
C 
CD      CALL IMPTDT ('EVOLUTION CALCULEE DANS SCHAST '//
CD                   'APRES ORTHOGONALISATION '//IDPROG,
CD                    DM(FOTDEP), NPICET, 1)
C 
CD      NBFEPS = 0
C 
CD      CALL RCHAMP (6, NPICET, DM(FOTDEP), DM(FTECH6))
C 
CD      CALL RCHARE (8, LONEPS, DM(EPSTRA), DM(EPSCH8))
C 
CD      IF (NBINT .GT. 0) THEN
C 
CD        CALL RCHARE (9, LONSAU, DM(SAUTRA), DM(SAUCH9))
C 
CD      END IF
C 
C       Appel a la routine remplisssant le fichier scb-numero d'etape locale
C       en y rangeant - delta( sigcha) - K ( delta-tild( epsn) )
C 
CD      NDEBEP = DEADTR-NBFEPS+1
C 
CD      CALL RSICBA (0, ADCOU, ADINT, ADSOU, ADSNT,
CD                   'sigver', 'sgnver',
CD                   1, NBFEPS, NDEBEP, DEADTR,
CD                   EPSCH8, SAUCH9, FTECH6, 1, NPICET)
C 
CD      CALL SCHAST (0, 'sigver', 'sgnver', 1,
CD                   DM(EPSTRA), DM(SAUTRA), 1, NPICET,
CD                   DM(FOTDEP))
C 
CD      CALL IMPTDT
CD               ('EVOLUTION CALCULEE APRES UN PREMIER PASSAGE NUL? '//IDPROG,
CD                 DM(FOTDEP), NPICET, 1)
C 
CD      NBFEPS = 0
CD      DEADTR = 1
CD      SAADTR = 1
C 
CD    ENDIF
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
