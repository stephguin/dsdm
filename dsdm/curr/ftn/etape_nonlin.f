C     Cette routine delete les fichiers directs crees dans
C     le cas d'un calcul sans reprise
C 
      SUBROUTINE DELFIC
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C 
C     Declaration des parametres locaux
C     """""""""""""""""""""""""""""""""
C 
      CHARACTER*20  NOMFIC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DELFIC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      NOMFIC = 'admcou'
C 
      CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
      NOMFIC = 'admint'
C 
      CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
      NOMFIC = 'sigmch'
C 
      CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
      NOMFIC = 'epsich'
C 
      CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
CD    NOMFIC = 'sigver'
C 
CD    CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
CD    NOMFIC = 'epsver'
C 
CD    CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
CD    NOMFIC = 'sinver'
C 
CD    CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
CD    NOMFIC = 'sauver'
C 
CD    CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
      NOMFIC = 'sinoch'
C 
      CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
      NOMFIC = 'sautch'
C 
      CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C 
C     Nouvelle version de PRELIM, routine preliminaire a l'etape locale.
C 
      SUBROUTINE PRELIM 
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  IUNIT1, IUNSIG
C 
      INTEGER  LONDEP, LONDER, LONEPS, LONSAU, LONECF
C 
      CHARACTER*20  NOM, NOMFIC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='PRELIM')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      LONDEP = NDDL*NBMAT
      LONDER = NDDL*NTETA
      LONEPS = NEPS*NGAU1*NTETA
      LONSAU = NSAU*NGAU2*NTETA
C 
C     Choix du nom des fichier dans Q-chapeau : (3)
C     Ouverture du fichier pour les deformations et les contraintes admissibles
C     Modification par rapport aux autres version :
C     on ouvre en direct non formatte => un seul fichier
C 
      NOM = 'admcou'
C 
      NOMFIC = NOM
      LONECF = 12*NPICET
C 
      CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNIT1)
      CALL FERFIC (3, IUNIT1, IDPROG)
C 
      NOM = 'admint'
C 
      NOMFIC = NOM
      LONECF = 6*NPICET
C 
      CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNIT1)
      CALL FERFIC (3, IUNIT1, IDPROG)
C 
C     CREATION du fichier pour les contraintes
C 
      NOM = 'sigmch'
      NOMFIC = NOM
      LONECF = 6*NPICET
      CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNSIG)
      CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     CREATION du fichier pour les deformations
C 
      NOM = 'epsich'
      NOMFIC = NOM
      LONECF = 6*NPICET
      CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNSIG)
      CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     CREATION du fichier pour la verification de rsicba contraintes
C 
CD    NOM = 'sigver'
CD    NOMFIC = NOM
CD    LONECF = 6*NPICET
CD    CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNSIG)
CD    CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     CREATION du fichier pour la verification de rsicba deformations
C 
CD    NOM = 'epsver'
CD    NOMFIC = NOM
CD    LONECF = 6*NPICET
CD    CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNSIG)
CD    CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     CREATION du fichier pour les sauts normaux
C 
      NOM = 'sinoch'
      NOMFIC = NOM
      LONECF = 3*NPICET
      CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNSIG)
      CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     CREATION du fichier pour les sauts aux interfaces
C 
      NOM = 'sautch'
      NOMFIC = NOM
      LONECF = 3*NPICET
      CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNSIG)
      CALL FERFIC (3, IUNSIG, IDPROG)
C 
CD    NOM = 'sinver'
CD    NOMFIC = NOM
CD    LONECF = 3*NPICET
CD    CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONECF, IUNSIG)
CD    CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     CREATION du fichier pour les sauts aux interfaces
C 
CD    NOM = 'sauver'
CD    NOMFIC = NOM
CD    LONECF = 3*NPICET
CD    CALL CFDDNF( 3 ,'q-chapeau' ,NOMFIC , 6 , LONECF ,IUNSIG )
CD    CALL FERFIC( 3 , IUNSIG , IDPROG )
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Subroutine preliminaire de remplissage du fichier des champs
C     admissibles avant la premiere etape locale. On remplit de maniere
C     directe les fichiers :
C 
C     QCADMI : Les accroissement des quantites admissibles pour les couches
C              deformations puis contraintes (12*npicet) par enregistrement
C 
C     QIADMI : Les accroissement des quantites admissibles pourles interfaces
C              sauts puis contraintes normales(6*npicet) par enregistrement
C 
      SUBROUTINE RPCADS
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  NUCOU, PGAU1
      INTEGER  TETA, LONECF, TEMPS
C 
C     Steph.
C 
      INTEGER ZOB, ZOBSIG, DZOBSI, I
C 
C     Pour les interfaces
C 
      INTEGER  NUINT, ADSIGN, ADSAUT
      INTEGER  NUENRS
      INTEGER  SAUAPG, SAPADM, DSAUAP
      INTEGER  AM2LC, ADM2LC
C 
      INTEGER  IUNIT1
C 
C     Pour les nouveaux calculs en temps
C 
      INTEGER  EPSAPG, EPPADM
C 
      INTEGER  SIGC11, FTECH7, EPSCH8, SGNC12, SAUCH9, FTSC10
C 
      INTEGER  ADEPS, ADSIG, DEBQAD, DSIGAP, SIGAPG
C 
      INTEGER  DEPSAP, LONGAU
C 
      CHARACTER*6 IDPROG, NOM, NOMFIC
C 
      PARAMETER (IDPROG='RPCADS')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C -----------------------------------------------------------------------
C 
C     RECHERCHE DES ADRESSES DE TABLEAUX :
C 
C       - contraintes admissibles (nsig, nteta, ngau1, nfotps)
C       - deformations admissibles (neps, nteta, ngau1, nfotps)
C       - accroissements des fonctions du temps donnees (npicet, nbfodo)
C       - intervalles de temps (npicet)
C       - angles des bandes (nteta)
C       - numerateur erreur locale en temps (npicet+1)
C       - denominateur en temps de l'erreur locale en temps (npicet+1)
C 
C     Tableaux des deformations et contraintes admissibles
C 
C     Tableaux d'initialisation  des quantites locales pour le temps zero
C 
      CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
      ADEPS = EPSCH8
      CALL ADTBDM ('SIG-AD-TOT', SIGC11)
      ADSIG = SIGC11
C 
C     POUR LES INTERFACES
C 
      IF (NBINT .GT. 0) THEN
        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
        ADSAUT = SAUCH9
        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
        ADSIGN = SGNC12
      ENDIF
      CALL ADTBDM ('TEMPS-REEL', FTECH7)
      CALL ADTBDM ('TEMP-SI-RE', FTSC10)
C 
C     Creation d'un tableau provisoire pour ranger pour toutes les
C     quantites aux piquets de temps
C     EPSAPG QUI CONTIENT :                                             X6*NPICET
C     SIGAPG QUI CONTIENT :                                             X6*NPICET
C     SIPADM               sigma point admissible                       X6*NPICET
C     EPPADM (SIPADM+6)    epsilon point admissible                     X6*NPICET
C                                                                      __________
C                                                                       24*NPICET
      CALL POUSMD (24*NPICET, EPSAPG)
      SIGAPG = EPSAPG+6*NPICET
C 
C     Quantites stockees en fin de boucle sur le temps
C 
      EPPADM = SIGAPG+6*NPICET
C 
C     Choix du nom des fichier dans Q-chapeau : (3)
C     Ouverture du fichier pour les deformations et les contraintes admissibles.
C     Modification par rapport aux autres versions :
C     on ouvre en direct_non formatte => un seul fichier
C 
      NOM = 'admcou'
C 
      NOMFIC = NOM
      LONECF = 12*NPICET
C 
C     DANS Q-CHAPEAU
C 
      CALL OFDDNF (3, NOMFIC, 6, LONECF, IUNIT1)
C 
      NUENRS = 1
C 
      LONGAU = XINTEG*YINTEG*NBCOL
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES COUCHES
C 
      DO NUCOU = 1, NBCOU
C 
C       BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
        DO PGAU1 = 1, LONGAU
C 
C         BOUCLE ii SUR LES ANGLES
C 
          DO TETA = 1, NTETA
C 
C           On calcule la valeur des champs admissibles pour tous les piquets de
C           temps :  caracterise la position dans le tableau des contraintes
C           ou des deformations admissibles de la 1ere valeur interessante
C 
            CALL ADCTPS (CHARAX, CHARAX, NEPS,
     &                   DM(ADEPS), 1, 1,
     &                   DM(ADSIG), 1, 1,
     &                   DM(FTECH7), 1, 1,
     &                   DM(FTSC10), 1, 1,
C                        on recupere
     &                   DM(EPSAPG), DM(SIGAPG))
C 
C           Pour aller lire le bon point de Gauss dans les tableaux de stockage
C           des deformations et des contraintes.
C 
            ADEPS = ADEPS + CHARAX*NEPS
            ADSIG = ADSIG + CHARAX*NEPS
C 
C           Initialisation des champs avec decalage :
C           EPSAPG <==> EPPADM    EPSILON POINT ADMISSIBLE
C           SIGAPG <==> SIPADM    sigma point admissible
C 
            DEPSAP  = EPSAPG
            DSIGAP  = SIGAPG
            DEBQAD  = EPPADM
C 
C           BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
            DO  TEMPS = 1, NPICET
C 
C             Rangement direct des quantites admissibles provenant
C             de ADCTPS de la facon suivante : - (neps, nsig) X NPICET
C 
              CALL COPITD (6, DM(DEPSAP), DM(DEBQAD))
C 
              DEBQAD = DEBQAD+6
              DEPSAP = DEPSAP+6
C 
              CALL COPITD (6, DM(DSIGAP), DM(DEBQAD))
C 
C             Decalage d'indice des tableaux pour aller lire au bon
C             endroit dans EPSAPG, DSGIAP, DEPPAD, DSIPAD
C 
              DEBQAD = DEBQAD+6
              DSIGAP = DSIGAP+6
C 
C           FIN DE BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
            END DO
C 
            CALL EFDDNF (IUNIT1, DM(EPPADM), EPPADM, LONECF, NUENRS)
CD  	    CALL IMPTDT ('VERIF 4 QADMPR ', DM(EPPADM), 1, LONECF)
            IF ((NUCOU.EQ.1).AND.(PGAU1.EQ.1).AND.(TETA.EQ.1)) THEN
              CALL POUSMD (10, ZOBSIG)
	      ZOB=EPPADM
	      DZOBSI=ZOBSIG
	        DO I = 1, 10
	          DM(ZOBSIG)=DM(ZOB)
	          ZOBSIG=ZOBSIG+1
	          ZOB=ZOB+1
	        END DO
            END IF
            NUENRS = NUENRS +1
C 
C         FIN DE BOUCLE SUR LES ii ANGLES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES COUCHES
C 
      END DO
C 
C     Fermeture de l' unite  admcou
C 
      CALL FERFIC (3, IUNIT1, IDPROG)
CD    CALL IMPTDT ('VERIF 2 CHEPSSIG '//IDPROG, DM(DZOBSI), 1, 10)
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
C     TEST SUR LE NOMBRE D'INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
        NOM = 'admint'
C 
        NOMFIC = NOM
        LONECF = 6*NPICET
C 
C       DANS Q-CHAPEAU
C 
        CALL OFDDNF (3, NOMFIC, 6, LONECF, IUNIT1)
C 
        NUENRS = 1
C 
C       DSACAC    delta saut  point chapeau actuel   XDT           X 3*NPICET
C       DSPCAC    - delta sigma point chapeau actuel XDT           X 3*NPICET
C                                                                 _____________
C                                                                   12*NPICET
        CALL POUSMD (12*NPICET, SAUAPG )
        SIGAPG = SAUAPG+3*NPICET
C 
C       Quantites stockees en fin de boucle sur le temps
C 
        SAPADM = SIGAPG +3*NPICET
C 
        LONGAU = XINTEG*NBCOL
C 
        LONECF = 6*NPICET
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT=1, NBINT
C 
C         BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
          DO PGAU1  = 1, LONGAU
C 
C           BOUCLE ii SUR LES ANGLES
C 
            DO TETA = 1, NTETA
C 
              CALL ADCTPS (CHARAX, CHARAX, NSAU,
     &                     DM(ADSAUT), 1, 1,
     &                     DM(ADSIGN), 1, 1,
     &                     DM(FTECH7), 1, 1,
     &                     DM(FTSC10), 1, 1,
C                          on recupere
     &                     DM(SAUAPG), DM(SIGAPG))
C   
              ADSAUT = ADSAUT + CHARAX*NSAU
              ADSIGN = ADSIGN + CHARAX*NSAU
C 
C             initialisation des champs avec decalage :
C             SIGAPG <==> SAPADM    saut point admissible
C                         SIPADM    sigma point admissible
C 
              DSAUAP  = SAUAPG
              DSIGAP  = SIGAPG
              DEBQAD = SAPADM
C 
C 
C             BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
              DO  TEMPS = 1, NPICET
C 
C               Rangement sequentiel des quantites admissibles provenant
C               de ADITPS de la facon suivante : - ( nsau , nsgn ) X NPICET
C 
                CALL COPITD (3, DM(DSAUAP), DM(DEBQAD))
C 
                DEBQAD = DEBQAD +3
C 
                CALL COPITD (3, DM(DSIGAP), DM(DEBQAD))
C 
C               Decalge d'indice des tableaux pour aller lire au bon endroit dans
C               SAUAPG, DSGIAP, DEPPAD, DSIPAD
C 
                DEBQAD = DEBQAD +3
                DSAUAP   = DSAUAP+3
                DSIGAP   = DSIGAP+3
C 
C             FIN DE BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
              END DO
C 
              CALL EFDDNF (IUNIT1, DM(SAPADM), SAPADM, LONECF, NUENRS)
              NUENRS = NUENRS +1
C 
C           FIN DE BOUCLE SUR ii LES ANGLES
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
          END DO
C 
C       FIN DE BOUCLE SUR LES INTERFACES
C 
        END DO
C 
C     Fermeture de l' unite  admint
C 
      CALL FERFIC (3, IUNIT1, IDPROG)
C 
C     FIN DE TEST SUR LE NOMBRE D'INTERFACES
C 
      END IF
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
      SUBROUTINE TRAREP (NETCHA, AGDCHA, NOUCHR, NOUCHA)
C 
C     Routine de traitement des options de reprise.
C 
C     Pour remettre a zero le nombre des champs admissibles crees
C     apres avoir recalcule la solution admissible par etloca.
C     L'interet est de pouvoir continuer le calcul.
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      LOGICAL NETCHA
C 
C     Pour agrandir les tableaux temps espace
C 
      LOGICAL AGDCHA
C 
C     CHARAX est modifie en  NOUCHR
C 
      INTEGER NOUCHR
C 
C     CHAMAX est modifie en  NOUCHA
C 
      INTEGER NOUCHA
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TRAREP')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (NETCHA) THEN
C 
        CALL TASCHA
C 
      END IF
C 
      IF (AGDCHA) THEN
C 
        CALL MNBCHA (NOUCHR, NOUCHA)
C 
      END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
