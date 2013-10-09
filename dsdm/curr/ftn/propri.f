C     VERSION DU     17 /  10  / 86
C 
      SUBROUTINE PROPRI
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'typcal.h'
      include 'strategie_calcul.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*40    MOT
      LOGICAL         FIN
      INTEGER         LECINT, I, ADEFFO
      INTEGER         ADDEPL, ADDEFA, ADSIGA
      INTEGER         ADSAUT, ADREEL, ADSIGN
      LOGICAL         LECLOG
      LOGICAL         OK
C 
C     REPMAT pour reprendre les matrices dans matricesm;
C     REPGLO reprise a l'etape globale NBETDE avant l'etape locale;
C     REPLOC reprise a l'etape globale NBETDE apres l'etape locale,
C     c'est a dire avant l'etape globale NBETDE+1;
C     NBETDE le numero de l'etape globale;
C     tild2 est le chemin caracterisant le fichier de reprise.
C 
      LOGICAL REPMAT, REPGLO, REPLOC
      INTEGER NBETRE, TYPE, LGCARG
C 
      LOGICAL SUITE
C 
C     Pour remettre a zero le nombre des champs admissibles crees
C     apres avoir recalcule la solution admissible par NETLOC.
C     L'interet est de pouvoir continuer le calcul.
C 
      LOGICAL NETCHA
C 
C     Pour lorsqu'on recalcule la solution admissible par NETLOC mettre a zero,
C     apres l'instabilite, les fonctions du temps correspondant aux accroissements
C     admissibles a zero, a partir d'une erreur trop importante.
C 
      LOGICAL NETINS
C 
C     Pour agrandir les tableaux temps-espace
C 
      LOGICAL AGDCHA
C 
C     CHARAX est modifie en NOUCHR
C 
      INTEGER NOUCHR
C 
C     CHAMAX est modifie en NOUCHA
C 
      INTEGER NOUCHA
C 
      CHARACTER*6     IDPROG, NOMFIC
C 
C     Pour la solution isotrope
C 
      DOUBLE PRECISION EPSILO
      PARAMETER       (IDPROG='PROPRI')
      PARAMETER       (EPSILO=1.D -3)
C 
C     Pour la lecture de la liste d'entier
C 
      INTEGER NLISTE, DEVALT
      PARAMETER (NLISTE=100)
      DOUBLE PRECISION LECDBL
C 
      INTEGER AM2LC, ADM2LC, AM2REE, ADM2RE
      INTEGER LONRES
C 
C     Pour le rangement des champs
C 
      INTEGER DEPCH0, EPSCH8, SIGC11, SAUCH9, SGNC12
C 
C     POUR TEST SUR LES NOUVELLES ROUTINES
C 
      INTEGER FTREE7, FTRE10, APRECI
C 
C     RIGIDITE DE PENALISATION
C 
      LOGICAL PENAL
C 
      DOUBLE PRECISION P1, P2, P3
C 
      INTEGER ADHOII, LONINT, ADSOUI, ADHOIN, ADSOIN, K
      INTEGER TUSED, DT, DTH, DTS, T1, T0
      INTEGER DEBNET, LONECF, IUNIT1
      CHARACTER*3 CARNET
      CHARACTER*6 NOM
C 
C     Pour se faire de la place dans les tableaux M ET DM
C     Pour l'instant, ADM2=0 et AM2=0
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C     On emplile le nom du sous-programme 
C     (~lmtutils/interpret/devel/ftn/a_debug_trac_lib.f)
C 
      CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C     initialisation des adresses des tableaux
C -----------------------------------------------------------------------
C 
      T0 = TUSED ()
C 
C     Initialise : premieres et dernieres adresses libres dans M et DM,
C     longueurs des tableaux contenus, longueurs des tableaux Kon, nombre
C     de termes TSF, coordonnees et poids des points de Gauss, PI, nombre
C     de fichiers ouvert, nombre d'etapes locales effectuees.
C 
      CALL INITIA 
C 
C     Attention : Longueur tild2 inferieure ou egale a 120
C 
      DO I = 1, 120
        tild2(i:i)=' '
      END DO
C 
C     Pour une reprise apres un premier calcul :
C     REPMAT pour reprendre les matrices dans matricesm
C     REPGLO reprise a l'etape globale NBETDE avant l'etape globale
C     REPLOC reprise a l'etape globale NBETDE apres l'etape locale
C     NBETDE le numero de l'etape globale
C 
      CALL LECSEQ ('CARAC_REPRISE ',
     &             'REPRISE DES MATRICES A UNE ETAPE DONNEE')
C 
      CALL MESSAO 
     &         ('NBETDE = 0 POUR UNE REPRISE DE LA SOLUTION ELASTIQUE')
C 
      REPMAT = LECLOG ('REPMAT')
      REPGLO = LECLOG ('REPGLO')
      REPLOC = LECLOG ('REPLOC')
      NBETRE = LECINT ('NBETDE')
C 
      CALL LECCHA ('tild2',tild2)
      LTILD2 = LGCARG (tild2,120)
C 
      IF (.NOT. STANDA) THEN
          CALL LECSEQ ('TRAIT-CHAMP ',
     &              'PERMET LA MODIFICATION DES FONCTIONS TEMPS-ESPACE')
      ELSE
          CALL LFORCS ('POINT-REPR,,,,,,',16,MOT,OK)
      ENDIF
C 
C     Pour nettoyer les champs admissibles "negligeables"
C 
      NETCHA = LECLOG ('NETCHA')
C 
C     Pour, lorsqu'on recalcule la solution admissible par NETLOC, mettre a zero
C     apres l'instabilite, les fonctions du temps correspondant aux accroissements
C     admissibles a zero.
C 
      CALL MESSAO ('L''OPTION NETINS NE FONCTIONNE PAS ')
C 
      NETINS = LECLOG ('NETINS')
C 
C     AGDCHA pour augmenter le nombre de fonctions temps-espace
C     BIGNET pour faire un grand nettoyage
C     NOUCHR nouveau CHARAX, nombre de champs admissibles maxi crees (1+nbiteratio nglobale)
C     NOUCHA nouveau CHAMAX, nombre de champs admissibles a zero +1 maxi crees (mini=2)
C 
      NETINS = .NOT. NETINS
      AGDCHA = LECLOG ('AGDCHA')
      BIGNET = LECLOG ('BIGNET')
      NOUCHR = LECINT ('NOUCHR')
      NOUCHA = LECINT ('NOUCHA')
C 
       CALL LECSEQ ('SYMETRIES', 'ON ENTRE ICI TOUT CE QUI CONCERNE '//
     &               'LES SYMETRIES DU PROBLEME')
C 
C     SYMX indique si la solution est symetrique par rapport a x
C     SYMY indique si la solution est symetrique par rapport a y
C     SYMO indique si la solution est symetrique par rapport a l'origine
C     FORS on force la symetrie (en imposant la symetrie des efforts, meme residuels))
C 
      SYMX   = LECLOG ('SYMX')
      SYMY   = LECLOG ('SYMY')
      SYMO   = LECLOG ('SYMO')
      FORS   = LECLOG ('FORS')
      
      CALL LECSEQ ('STRATEGIE-CALCUL', 'ATTENTION EN CAS DE REPRISE '//
     &             'IL FAUT VERIFIER QU''EN CAS D''ETAPE '//
     &             'PRELIMINAIRE, L''OPTION CHOISIE EST LA MEME QUE '//
     &             'L''OPTION INITIALE')
C 
C     LETAPP indique si on effectue une etape preliminaire
C     UPERIO indique si une seule periode en temps est imposee
C     UNINST indique si on veut decouper uniquement a l'instabilite
C     NBCINM nombre d'iteration maxi pour le probleme cinematique
C     NBSTAM nombre d'iteration maxi pour le probleme en statique
C     NUM type de numerotation : 1 couche, 2 colonne, -1 determination automatique
C     VALPRO indique si on traite l'etape globale par valeurs propres
C 
      LETAPP = LECLOG ('LETAPP')
      UPERIO = LECLOG ('UPERIO')
      UNINST = LECLOG ('UNINST')
      NBCINM = LECINT ('NBCINM')
      NBSTAM = LECINT ('NBSTAM')
      NUM    = LECINT ('NUM')
      VALPRO = LECLOG ('VALPRO')
      NBVALP = LECINT ('NBVALP')
C 
      IF (REPGLO .OR. REPLOC)  GOTO 1515
C 
C     Nettoyage des fichiers crees par les calculs precedents
C 
      CALL DELFIC
C 
C     Evaluation de la place necessaire au calcul
C 
      CALL DETAIL 
C 
C     Initialisation et lecture des donnees geometriques
C 
      CALL DONGEO
C 
      IF (REPMAT) THEN
        CHADIR(1)   = 'matricesm'
        NBFICH(1)   = 1
        NOMFIC      = 'matiso'
        CHAFIC(1,1) = NOMFIC
        CALL REPLAC
        CALL REPLAI
      END IF
C 
C -----------------------------------------------------------------------
C     Sequence de rentree des donnees materiau
C -----------------------------------------------------------------------
C 
      CALL LECSEQ ('COMPORTEMENT','RENTREE DU COMPORTEMENT')
C 
      LENDCO  = LECLOG ('LENDCO')
      LRUPCO  = LECLOG ('LRUPCO')
      LPLACO  = LECLOG ('LPLACO')
      LENINT  = LECLOG ('LENDIN')
      LPLAIN  = LECLOG ('LPLAIN')
C 
C     En cours
C 
CD    LMOSIG  = LECLOG ('LMOSIG')
CD    LRETAR  = LECLOG ('LRETAR')
CD    LRETIN  = LECLOG ('LRETIN')
CD    LPLAEF  = LECLOG ('LPLAEF')
CD    LPLEFF  = LECLOG ('LPLEFF')
C 
C ----------------------------------------------------------------------
C     LECTURE ET INITIALISATION DES CARACTERISTIQUES NON LINEAIRES
C ----------------------------------------------------------------------
      CALL DMATNL
C 
C ----------------------------------------------------------------------
C     DONNEE DE L'EVOLUTION EN TEMPS
C     DONNEE DES EFFORTS ET DES DEPLACEMENTS
C ----------------------------------------------------------------------
      CALL DONTEM
      CALL DONEFD
C 
C ----------------------------------------------------------------------
C     REMPLISSAGE DES TABLEAUX LIES AUX POINTS DE GAUSS
C ----------------------------------------------------------------------
      CALL REMRPG
C 
C ----------------------------------------------------------------------
C     POUR POUVOIR TRAVAILLER AVEC DES ZONES PRE-DELAMINEES
C ----------------------------------------------------------------------
      CALL DONDEL
C 
C ----------------------------------------------------------------------
C     CALCUL DU PROFIL
C ----------------------------------------------------------------------
      CALL PROFIL
C 
C ----------------------------------------------------------------------
C     CALCUL DES MATRICES K0N
C 
C     MATCON <=> Calcul des termes independants des elements
C     CALK0  <=> Calcul des matrices
C ----------------------------------------------------------------------
      CALL MATCON
C 
C ----------------------------------------------------------------------
C     Sequence de calculs des seconds membres pour les efforts imposes :
C 
C     EFFINT <=> Calcul des tableaux constants
C     EFFIMP <=> Calcul et assemblage des efforts sous forme develop
C ----------------------------------------------------------------------
      IF (.NOT. STANDA) THEN
        CALL LECSEQ ('EFFORTS', 'CALCULS DES EFFORTS')
      END IF
      CALL EFFINT
      CALL EFFIMP
C 
      CALL ADTBDM ('MAT-EFFORT', ADEFFO)
C 
C ----------------------------------------------------------------------
C     Sequence de modification des seconds membres pour les
C     deplacements imposes :
C 
C     DDLBLO <=> LISTE DES DDL-BLOQUES
C     DEPIMP <=> TRAITEMENT DES DONNEES
C     ASDEPI <=> MODIFICATION
C ----------------------------------------------------------------------
C 
      CALL DDLBLO 
      CALL DEPBLO
      CALL ASDBLO
C 
C -----------------------------------------------------------------------
C     Sequence de calcul des matrices si il n'y a pas de reprise des
C     matrices (si repmat= .false.)
C -----------------------------------------------------------------------
C 
      CALL GESTDP ('PRECISIONS', NBMAT, APRECI)
      DO I = 0, NBMAT -1
        DM(APRECI+I) = 1.D100
      END DO
C 
      IF (.NOT. REPMAT) THEN
        NOMFIC = 'matiso'
C 
C     DANS Q-CHAPEAU
C 
        CALL DELETF (1, 'matricesm', NOMFIC, 6)
        CALL CALK0
      END IF
C 
C ----------------------------------------------------------------------
C     TERMES DE PENALISATION POUR VOIR PAR GRADIENT CONJUGUE
C ----------------------------------------------------------------------
C 
      CALL LECSEQ ('PENALISATION', 'POUR UNE INTERFACE DE RIGIDITE '//
     &             'INFINIE ?')
C 
      PENAL = LECLOG ('PENAL')
      P1 = LECDBL ('P1')
      P2 = LECDBL ('P2')
      P3 = LECDBL ('P3')
C 
      IF (PENAL) THEN
        CALL ADTBDM ('HOO-INTERF', ADHOII)
        CALL LONGDP ('HOO-INTERF', LONINT)
        CALL ADTBDM ('SOU-INTERF', ADSOUI)
        CALL ADTBDM ('HOIN-ORTHO', ADHOIN)
        CALL ADTBDM ('SOIN-ORTHO', ADSOIN)
        K = 0
        DO I = 1, LONINT/3
          DM(ADHOIN+K) = P1
          DM(ADHOII+K) = (P1+P2)/2.D0
          DM(ADSOUI+K) = (1.D0/P1+1.D0/P2)/2.D0
          DM(ADSOIN+K) = 1.D0/P1
          K = K+1
          DM(ADHOIN+K) = P2
          DM(ADHOII+K) = P3
          DM(ADSOUI+K) = 1.D0/P3
          DM(ADSOIN+K) = 1.D0/P2
          K = K+1
          DM(ADHOIN+K) = P3
          DM(ADHOII+K) = (P1-P2)/2.D0
          DM(ADSOUI+K) = (1.D0/P1-1.D0/P2)/2.D0
          DM(ADSOIN+K) = 1.D0/P3
          K = K+1
        END DO
      END IF
C 
C ----------------------------------------------------------------------
C     REMPLISSAGE DU TABLEAU DES VALEURS DES DEPLACEMENTS PAR DEVELOPPEMENT
C ----------------------------------------------------------------------
      CALL DEPIMP
C 
C ----------------------------------------------------------------------
C     POUR ASSEMBLER LES DEPLACEMENTS IMPOSES NON NULS
C ----------------------------------------------------------------------
      CALL ASDEPI
C 
C ----------------------------------------------------------------------
C     MODIFICATION DES SECONDS MEMBRES EN FONCTION DES DEPLACEMENTS
C     LUS PAR FICHIER (TRACTION ou FLEXION)
C ----------------------------------------------------------------------
      IF (DONFIC .AND. TRACT .AND. (.NOT. DONVOL)) CALL MDEFIT
      IF (DONFIC .AND. (.NOT. TRACT) .AND. (.NOT. DONVOL)) CALL MDEFIF
C      IF (DONFIC .AND. (.NOT. TRACT) .AND. DONVOL) CALL MDEFIV
C 
C -----------------------------------------------------------------------
C     SEQUENCE DE CALCUL DES EFFORTS ISOTROPES
C -----------------------------------------------------------------------
      IF (.NOT. STANDA) THEN
      CALL LECSEQ ('CALISO', 'POUR CALCULS ISOTROPES')
      END IF
C 
C -----------------------------------------------------------------------
C     SEQUENCE DE CALCUL DES DEVELOPPEMENTS SOUS FORME DEVELOPPEE
C     DANS LA SEQUENCE D'INITIALISATION
C     CALISO  <=> CALCUL DE TOUS CES DEVELOPPEMENTS
C -----------------------------------------------------------------------
      CALL ADTBDM ('MAT-DEPLAC', ADDEPL)
      CALL CALISO (DM(ADEFFO), DM(ADDEPL))
C 
C ----------------------------------------------------------------------
C     Calcul de la solution pour le comportement elastique orthotrope
C     CALCUL de la 1ere discretisation en temps et des fonctions du
C     temps associees (appel a DISTEM).
C ----------------------------------------------------------------------
C     Pour ne pas garder les tableaux partiels admissibles
C 
      CALL ENPOUB (AM2REE, ADM2RE)
C 
C     TABLEAUX DES CONTRAINTES ET DEFORMATIONS ADMISSIBLES TOTALES
C 
      CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
      CALL ADTBDM ('SIG-AD-TOT', SIGC11)
C 
C     POUR LES INTERFACES
C 
      IF (NBINT.GT.0)THEN
        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
      ENDIF
C 
      LONRES = 2*NTETA*NFOTPS*(NGAU1*NEPS+NGAU2*NSAU)
C 
      CALL POUSMD (LONRES, ADDEFA)
      ADSIGA = ADDEFA + NTETA*NGAU1*NEPS*NFOTPS
      ADSAUT = ADSIGA + NTETA*NGAU1*NEPS*NFOTPS
      ADSIGN = ADSAUT + NTETA*NGAU2*NSAU*NFOTPS
C 
      CALL SOLORT (DM(ADDEPL),
C                  on recupere
     &             DM(ADDEPL), DM(ADDEFA), DM(ADSIGA),
     &             DM(ADSAUT), DM(ADSIGN))
C 
      CALL RCHARE (8, NTETA*NGAU1*NEPS, DM(ADDEFA), DM(EPSCH8))
      CALL RCHARE (11, NTETA*NGAU1*NEPS, DM(ADSIGA), DM(SIGC11))
C 
      IF (NBINT .GT. 0) THEN
        CALL RCHARE (9, NTETA*NGAU2*NSAU, DM(ADSAUT), DM(SAUCH9))
        CALL RCHARE (12, NTETA*NGAU2*NSAU, DM(ADSIGN), DM(SGNC12))
      END IF
C 
      CALL SOPOUB (AM2REE, ADM2RE)
C 
      IF (.NOT. STANDA) THEN
        CALL LECSEQ ('DEPREE', 'POUR CALCULER LES DEPLACEMENTS REELS ')
      END IF
C 
C ----------------------------------------------------------------------- 
C     CALCUL DES DEPLACEMENTS REELS POUR NE PAS GARDER
C     LES TABLEAUX CORRESPONDANT A DEPREE
C -----------------------------------------------------------------------
C 
      CALL POUSMD (NTETA*NDDL*NFOTPS, ADREEL)
      CALL VRTSM (1, DM(ADDEPL), DM(ADREEL))
      CALL ADTBDM ('DEPLA-ADMI', DEPCH0)
      CALL RCHARE (0, NDDL*NTETA, DM(ADREEL), DM(DEPCH0))
C 
C     calcul des fonctions du temps et rangement pour les evolutions
C     des deformations et contraintes et des erreurs
C 
      CALL NDISTM
      CALL ADTBDM ('VALE-TEMPS', DEVALT)
C 
C ----------------------------------------------------------------------
C     Sauvegarde de l'etape globale (=> type = 0)  numero 0
C     c'est a dire de l'initialisation elastique
C ----------------------------------------------------------------------
C 
      CALL ADTBDM ('TEMPS-REEL', FTREE7)
      CALL ADTBDM ('TEMP-SI-RE', FTRE10)
      CALL ECFTGL (DM(FTREE7), DM(FTRE10))
C 
C     determination des proprietes d'instabilite, etc ...
C 
      CALL DTEXPR (0, 0)
C 
      T1 = TUSED ()
      DT = (T1-T0)/50
      DTH= DT/3600
      DTS= DT-3600*DTH
C       
      CALL MESSAO ('POUR L''ANALYSE ELASTIQUE')
      CALL IMPET  ('TEMPS UTILISE EN HEURES   ', DTH)
      CALL IMPET  ('TEMPS UTILISE EN SECONDES ', DTS)
C 
      CALL SAUCVR (0)
C 
      CALL ADTBDM ('DEPLA-ADMI', DEPCH0)
C 
1515  CONTINUE
C 
      IF (REPGLO) THEN
C 
	CALL REPCVR (0, NBETRE, .FALSE.)
        IF (NBETRE .EQ. 0) THEN
          CALL PRELIM 
        END IF
C 
C       Traitement des options de reprise
C 
        CALL TRAREP (NETCHA, AGDCHA, NOUCHR, NOUCHA)
C       Option grand nettoyage : on vire tous les accroissements construits
C       aux iterations precedentes, on repart sur la solution admissible totale
C       reconstruite a l'etape locale TYPE 1 ci-dessus.
C       A voir : ameliorer DELCHA pour garder certains de ces accroissements
C       (les meilleurs ?)
C 
        IF (BIGNET) THEN
          IF (NBNETT .GT. NBNEMX) THEN
            CALL MESSAO ('TROP DE NETTOYAGES !!')
          END IF
	  IF (NBNETT .EQ. 0) THEN
            CALL GESTEN ('LIS-BIGNET', NBNEMX, DEBNET)
	    M(DEBNET) = NBETGL
	  END IF
          IF (NBNETT .GE. 1) THEN
            CALL ADTBM ('LIS-BIGNET', DEBNET)
	    M(DEBNET+NBNETT) = NBETGL
          END IF
C 
          CALL IDENTI (NBETGL, CARNET)
          NOM = 'ac_'//CARNET
          LONECF = 12*NPICET
          CALL CFDDNF (3, 'q-chapeau', NOM, 6, LONECF, IUNIT1)
          CALL FERFIC (3, IUNIT1, IDPROG)
          NOM = 'ai_'//CARNET
          LONECF = 6*NPICET
          CALL CFDDNF (3, 'q-chapeau', NOM, 6, LONECF, IUNIT1)
          CALL FERFIC (3, IUNIT1, IDPROG)
        END IF
C 
	CALL LOCLOC (1)
	IF (BIGNET) THEN
	  CALL DELCHA
	  NBNETT = NBNETT + 1
	END IF
      END IF
C 
      IF (REPLOC) THEN
        CALL REPCVR (1, NBETRE, .FALSE.)
C 
C       Traitement des options de reprise
C 
C     Steph le 15/12/2000
C     Je vire l'option valpro : sur rubis la compil des routines dspev
C     de chez NETLIB se passe mal
C 
	CALL TRAREP (NETCHA, AGDCHA, NOUCHR, NOUCHA)
CD      IF (VALPRO) THEN
CD 	  CALL NEW_GCNLIN
CD 	ELSE 
	  CALL GCNLIN
CD	END IF
      END IF
C 
1     CALL PROMPT 
      CALL LECART (MOT, FIN)
C 
      IF (FIN) THEN
        CALL ERREUD (0, 'FIN DE FICHIER')
      ENDIF
C 
      IF (MOT(1:12) .EQ. 'ETAPE-LOCALE') THEN
C 
        IF (NBETLC .EQ. 0) THEN
C 
C         Ouverture des fichiers ...
C 
          CALL PRELIM 
C 
C         Remplissage du fichier des accroissements admissibles
C 
          CALL RPCADS
C 
        END IF
C 
        CALL LOCLOC (0)
C 
C       Pour savoir si une nouvelle iteration est necessaire
C 
c
C    emplacement ???? de DEPRPG POUR SAUVGARDER LA VALEUR DES DEPLACEMENT SUR LA MAILLAGE COTENANT
C    LES POINTS DE gAUSS 
C 
        CALL DEPRPG
C
        CALL SUICAL (SUITE)
C 
      ELSE IF (MOT(1:13) .EQ. 'ECRIT_REPRISE') THEN
C 
        TYPE = LECINT ('TYPE')
        CALL SAUCVR (TYPE)
C 
      ELSE IF (MOT(1:13) .EQ. 'ETAPE-GLOBALE') THEN
C 
C     Steph le 15/12/2000
C     Je vire l'option valpro : sur rubis la compil des routines dspev
C     de chez NETLIB se passe mal
C 
CD      IF (VALPRO) THEN
CD	  CALL NEW_GCNLIN 
CD	ELSE 
	  CALL GCNLIN
CD	END IF
C 
        T1 = TUSED()
        DT = (T1-T0)/50
        DTH= DT/3600
        DTS= DT-3600*DTH
        CALL IMPET ('APRES L''ETAPE GLOBALE NUMERO  : ', NBETGL)
        CALL IMPET ('TEMPS TOTAL UTILISE, HEURES   : ', DTH)
        CALL IMPET ('TEMPS TOTAL UTILISE, SECONDES : ', DTS)
C 
        CALL SAUCVR (0)
C 
      ELSE IF (MOT(1:6) .EQ. 'SORTIE') THEN
C 
        GOTO 2
C 
      END IF
C 
      GOTO 1
2     CONTINUE
      CALL RETOUD (IDPROG)
      CALL SOPOUB (AM2LC, ADM2LC)
      RETURN
900   CALL ERREUD (0, 'ERREUR D''OUVERTURE DANS '//IDPROG)
      END
