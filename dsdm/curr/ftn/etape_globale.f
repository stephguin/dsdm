C     Cette routine calcule le nombre d'iterations maxi du
C     gradient conjugue en fonction de l'erreur a l'etape locale.
C 
C     On envoie comme arguments :
C 
C     E ...... ERREUR  erreur a la derniere etape locale
C 
C     On recupere :
C 
C     S ...... NITEMX  nombre d'iterations maxi pour GCADMI
C     S ...... PREGLO  precision de la resolution globale
C 
      SUBROUTINE CNBITE (ERREUR, NITEMX, PREGLO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      INTEGER           NITEMX
C 
      DOUBLE PRECISION  ERREUR, PREGLO
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CNBITE')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      IF (ERREUR .GT. 1.D0) THEN
        NITEMX = 3
      ELSE IF (ERREUR .LE. 1.D0 .AND. ERREUR .GT. .9D0) THEN
        NITEMX = 6
      ELSE IF (ERREUR .GT. .8D0 .AND. ERREUR .LE. .9D0) THEN
        NITEMX = 9
      ELSE IF (ERREUR .GT. .7D0 .AND. ERREUR .LE. .8D0) THEN
        NITEMX = 12
      ELSE IF (ERREUR .GT. .6D0 .AND. ERREUR .LE. .7D0) THEN
        NITEMX = 15
      ELSE IF (ERREUR .GT. .5D0 .AND. ERREUR .LE. .6D0) THEN
        NITEMX = 18
      ELSE IF (ERREUR .GT. .4D0 .AND. ERREUR .LE. .5D0) THEN
        NITEMX = 21
      ELSE IF (ERREUR .GT. .3D0 .AND. ERREUR .LE. .4D0) THEN
        NITEMX = 24
      ELSE IF (ERREUR .GT. .2D0 .AND. ERREUR .LE. .3D0) THEN
        NITEMX = 27
      ELSE IF (ERREUR .LE. .2D0 .AND. ERREUR .GE. 0.D0) THEN
        NITEMX = 30
      ELSE
	CALL ERREUD (0, 'ERREUR NEGATIVE DANS '//idprog)
      END IF
C 
      PREGLO = ERREUR/2.D0
C 
      CALL IMPET ('NOMBRE MAXI D''ITERATIONS CALCULE DANS '// IDPROG,
     &             NITEMX)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule les differentes fonctions ALPHA admissibles
C     cinematiquement a zero pour les NFTGLO fonctions du temps a
C     l'iteration concernee de l'etape globale. Celles-ci sont rangees
C     (nddl, nbmat, nftglo). Ces fonctions sont calculees apres l'etape
C     d'optimisation preliminaire sur les fonctions admissibles a zero
C     obtenues grace aux iterations precedentes.
C 
C     On envoie comme arguments et on recupere :
C 
C     E ...... FOTEMP les fonctions du temps associees
C 
      SUBROUTINE GCNLIN
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'strategie_calcul.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  AM2LC , ADM2LC, ADREEL
      INTEGER  DEPSOL, EPSSOL, SAUSOL
C 
      LOGICAL  LTRACN, LTRACP
C 
      DOUBLE PRECISION NORME, INTERV
C 
C     Pour l'optimisation en temps
C 
C     Nombre d'iterations en temps max
C 
      INTEGER NITEMX, ADSMEG, ADINF
      INTEGER SGNZER, SIGZER, NBTOUR
      INTEGER ERRDEP, ERRSIG, ANUMER
      INTEGER FOTEST, FTREE7, NUMELO
      INTEGER EPSCH8, SGNC12, SIGC11, SAUCH9, FTRE10
      INTEGER LONEPSCH8, LONSIGC11
      INTEGER FTECH6, FTSCH5, DEPCH0
C 
C     Pour le calcul du nombre d'iterations dans GCADMI
C 
      INTEGER NITGCA, NBTTST
      DOUBLE PRECISION ERFINL
C 
C     Pour essayer de construire une contrainte rigoureusement sa a 0
C 
      INTEGER SIGRES, SGNRES
C 
C     Nombre d'terations en espace max
C 
      INTEGER NISPMX, NTOUTE
C 
      INTEGER AFEPSI, AFCONT, DFCONT, INFOGL
      INTEGER AERREU, ADINTE, ADPICE, DEBFT, FOCONS, I
      INTEGER FOTDEP, FCONPA, FOTDPP
C 
      DOUBLE PRECISION NORIN, NORAC, NORFT, NORINP, NORACP
      DOUBLE PRECISION PRECIS, PREGLO
C 
      LOGICAL REPRIS
C 
C     Longueurs des tableaux locaux
C 
      INTEGER LONEPS, LONSAU, LONDEP, LONDER, LONRES
      INTEGER LFTA0C, LFTA0D
C 
      DOUBLE PRECISION DBLLOC, NORM1
C 
C     Pour ne pas calculer les termes nuls
C 
      INTEGER   NBDEV, TNUDEV
C 
      CHARACTER*3 CARETG
C 
      INTEGER     NPICDE, NPICFI, PERIOD, NPERIO, APERIO
      INTEGER     TUSED, DT, DTH, DTS, T1, T0, T2
C 
C     Precision en temps
C 
      PARAMETER (PRECIS = 2.D-1)
      INTEGER LONTAB, AMATPR, BMATPR, LWORK, LIWORK, IWORK, WORK, DIMP
      INTEGER VCDEB, VPDEB, AMACUS, AMACUR
      INTEGER MCDEB, MCCURR, MCCURS, VPCURR, VCCURR, VPCURS
      DOUBLE PRECISION X
      CHARACTER*2 CHETLC
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='GCNLIN')
C 
      CHARACTER*3 BARATI
      INTEGER MIN, IFAIL, INFOIN, J, K, NBVAL
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL IDENTI (NBETGL, CARETG)
C 
      T0 = TUSED()
C 
C     Pour UPERIO = .TRUE. on redefinit le decoupage en temps
C     suivant l'erreur, le taux de variation de l'erreur, le passage
C     par une instabilite. A essayer ?
C 
      CALL DPERIO
C 
C -----------------------------------------------------------------------
C     PRELIMINAIRE I ON VA CHERCHER LES TABLEAUX QU'IL NOUS FAUT
C -----------------------------------------------------------------------
C 
      CALL ADTBM ('PERIOD -'//CARETG, APERIO)
      NPERIO = M(APERIO)
      NBETGL = NBETGL+1
      CALL IDENTI (NBETGL, CARETG)
C 
      CALL MESSAO (
     &  '******************************************************
     &  \
     &  \DEBUT DU TRAITEMENT DE L''ETAPE GLOBALE NUMERO '//CARETG //'
     &  \
     &  \******************************************************')
C 
C     Appel du tableau des valeurs precedentes de l'erreur,
C     le nombre de fonctions du temps utilisees pour les evolutions totales
C 
      CALL ADTBDM ('INF-ERREUR', ADINF)
      DM(ADINF+NPICET+1) = NBDPTR
      DM(ADINF+NPICET+2) = EVCOTR
      LONEPS = NTETA*NEPS*NGAU1
      LONDER = NDDL*NTETA
      LONDEP = NDDL*NBMAT
      LONSAU = NTETA*NSAU*NGAU2
C 
C     Appel des tableaux pour les valeurs des piquets, des intervalles,
C     eventuellement affectees par un redecoupage en temps DPERIO.
C 
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
      CALL ADTBDM ('VALE-TEMPS', ADPICE)
      INTERV = DM(ADINTE+1)
C 
C     Appel des caracteristiques des tableaux des fonctions du temps
C     pour les deltas admissibles
C  
      CALL INFODP ('TEMPS-SIGM', FTSCH5, LFTA0C)
      CALL INFODP ('TEMPS-EPSI', FTECH6, LFTA0D)
C 
C     Tableaux des fonctions du temps pour la solution totale admissible en deformations,
C     en contraintes
C 
      CALL ADTBDM ('TEMPS-REEL', FTREE7)
      CALL ADTBDM ('TEMP-SI-RE', FTRE10)
      NFTGLO = 1
C 
C     Tableau des numerateurs de l'erreur aux piquets a l'issue de l'etape locale
C 
      CALL ADTBDM ('NUMERA-LOC', NUMELO)
      CALL POUSMD (NPICET, ANUMER)
      CALL COPID (DM(NUMELO+1), DM(ANUMER), NPICET, 1)
C 
C -----------------------------------------------------------------------
C     PRELIMINAIRE II RESERVATION DE PLACE POUR LES TABLEAUX PROVISOIRES
C -----------------------------------------------------------------------
C 
C     On cree le tableau inf-gl (info etape globale) dans lequel on trouvera
C     pour l'etape globale courante :
C 
C        le nombre maxi d'iterations en espace
C        le nombre maxi d'iterations en temps
C        le nombre maxi de decoupages en temps
C        le nombre d'iterations en espace reellement effectue
C        le nombre de fonctions du temps creees par iteration
C 
C     Il y en a par iteration en espace (max 5) 3 fois le nombre d'iterations
C     en temps +3 => 15x5.
C     Pour les containtes on en a 3 par iteration en espace.
C 
      CALL GESTEN ('INF_GL_'//CARETG, 8, INFOGL)
      M(INFOGL  ) = 25
      M(INFOGL+1) = 4
      M(INFOGL+2) = 3
C 
C     Tableaux des fonctions servant aux test de convergence en temps.
C     Pour chaque iteration en temps on stocke 3 nouvelles fonctions :
C     celle calculee, celle normee, celle regularisee.
C 
C     195 et 45 ??????????????? a modifier (coherence avec NBDPTR)
C 
      CALL GESTDP ('FT_EPS_'//CARETG, 195*NPICET, AFEPSI)
      CALL GESTDP ('FT_CON_'//CARETG, 45*NPICET, AFCONT)
      DEBFT  = AFEPSI
      DFCONT = AFCONT
      NITEMX = 15
      NISPMX = M(INFOGL)
C 
C     Dans fotest on stocke deux evolutions successives en temps pour arreter le jeu temps-espace.
C 
      LONRES  = 5*NPICET
      CALL POUSMD (LONRES, FOTEST)
C 
C     Pour le calcul des normes en temps et espace de la contrainte delta(-sigba)
C 
      FOCONS = FOTEST+NPICET
      FOTDEP = FOCONS+NPICET
      FCONPA = FOTDEP+NPICET
      FOTDPP = FCONPA+NPICET
      DO I = 1, NPICET
        DM(FOCONS+I-1) = 1.D0
      END DO
C 
C -----------------------------------------------------------------------
C     DEBUT DE L'ETAPE GLOBALE
C -----------------------------------------------------------------------
C 
C     On determine l'evolution initiale; l'erreur de resolution n'ayant pas de
C     sens en temps, on choisit comme fonction de depart en temps pour toutes
C     les iterations l'erreur a l'etape globale normee en temps.
C 
      CALL ADTBDM ('ERREUR-LOC', AERREU)
      CALL COPID (DM(AERREU+1), DM(FOTDEP), NPICET, 1)
      CALL NORMFT (ADINTE, ADPICE, DM(FOTDEP), DM(FOTDEP), NORFT)
C 
C     Reservation de place pour les nouvelles contraintes
C 
      LONRES  = LONSAU+LONEPS+LONDEP
      CALL POUSMD (LONRES, ADSMEG)
      SIGZER = ADSMEG+LONDEP
      SGNZER = SIGZER+LONEPS
C 
C     Reservation de place et mise a zero du tableau des efforts globaux
C 
C     Tableau des developpements non nuls
C 
      CALL POUSME (NBMAT, TNUDEV)
C 
C     Calcul de l'effort correspondant a intemps(-delta(sigmch))
C     obtenu a l'etape locale
C 
      CALL EFFPRE (0, ADSMEG, DM(FOCONS), DM(SIGZER), DM(SGNZER),
     &             NBDEV, M(TNUDEV), NORIN)
      NORAC = NORIN
      CALL IMPDT ('NORME INITIALE CORRESPONDANT A sigmch ', NORIN)
C 
C     Calcul du nombre d'iterations maxi du gradient conjugue
C     en fonction de l'erreur a l'etape locale
C 
      ERFINL =   DM(AERREU+NPICET)
      CALL CNBITE (ERFINL, NITGCA, PREGLO)
      LONRES = 2*NPICET
      CALL POUSMD (LONRES, ERRDEP)
      ERRSIG = ERRDEP+NPICET
      IF (LETAPP .AND. (NBETLC .GT. 1) .AND. (.NOT. BIGNET)) THEN
C 
C       On passe par l'etape preliminaire
C 
        CALL ETAPRE
        T1 = TUSED()
        DT = (T1-T0)/50
        DTH= DT/3600
        DTS= DT-3600*DTH
C 
C       Calcul de l'effort correspondant a intemps(-delta(sigmch))
C       apres optimisation des fonctions de l'espace calculees aux
C       etapes globales precedentes.
C 
        CALL EFFPRE (0, ADSMEG, DM(FOCONS), DM(SIGZER), DM(SGNZER),
     &               NBDEV, M(TNUDEV), NORAC)
        CALL IMPDT ('NORME APRES ETAPRE CORRESPONDANT A SIGMCH ', NORIN)
      END IF
C 
      IF ((.NOT. LETAPP) .AND. (NBETLC .GT. 1)) THEN
C 
        NBFEPS =0
        NBFSIG =0
C 
        CALL BALAID (LFTA0C, DM(FTSCH5))
        CALL BALAID (LFTA0D, DM(FTECH6))
C 
      END IF
C 
      CALL MESSAO (
     &  '******************************************************
     &  \
     &  \DEBUT DU TRAITEMENT DU PROBLEME CINEMATIQUE
     &  \
     &  \******************************************************')
C 
      LONRES = LONDEP + LONDER + 2*(LONEPS+LONSAU)
C 
      CALL POUSMD (LONRES, DEPSOL)
C 
      EPSSOL = DEPSOL + LONDEP
      SAUSOL = EPSSOL + LONEPS
      ADREEL = SAUSOL + LONSAU
      SIGRES = ADREEL + LONDER
      SGNRES = SIGRES + LONEPS
C 
      CALL MESSAO (
     &  '******************************************************
     &  \
     &  \  DEBUT DE LA BOUCLE SUR LES DIFFERENTES PERIODES EN TEMPS
     &  \
     &  \******************************************************')
C 
      DO PERIOD = 1, NPERIO
C 
C       Pour la periode de temps consideree les piquets sont
C 
        NPICDE = M(APERIO+PERIOD)
        NPICFI = M(APERIO+PERIOD +1)-1
        PREGLO = DM(AERREU+NPICFI)/2.D0
C 
        DO I = 1, NPICDE-1
          DM(FCONPA+I-1) = 0.D0
        END DO
        DO I = NPICDE, NPICFI
          DM(FCONPA+I-1) = 1.D0
        END DO
        DO I = NPICFI+1, NPICET
          DM(FCONPA+I-1) = 0.D0
        END DO
        DO I = 1, NPICDE-1
          DM(FOTDPP+I-1) = 0.D0
        END DO
        DO I = NPICDE, NPICFI
          DM(FOTDPP+I-1) = DM(FOTDEP+I-1)
        END DO
        DO I = NPICFI+1, NPICET
          DM(FOTDPP+I-1) = 0.D0
        END DO
C 
        CALL NORMFT (ADINTE, ADPICE, DM(FOTDPP), DM(FOTDPP), NORFT)
        CALL COPID (DM(FOTDPP), DM(DEBFT), NPICET, 1)
C 
C       Calcul de l'effort correspondant a intemps(-delta(sigmch))
C       apres l'etape preliminaire, et mesure de sa norme pour
C       l'intervalle de temps considere.
C 
        CALL EFFPRE (0, ADSMEG, DM(FCONPA), DM(SIGZER), DM(SGNZER),
     &               NBDEV, M(TNUDEV), NORINP)
        CALL IMPDT ('NORME EN DEBUT ETAGLO CORRESPONDANT A sigmch ',
     &               NORINP)
        NORACP = NORINP
        NBTOUR = 1
C 
C       La version actuelle consiste a obtenir la precision souhaitee sur
C       chaque periode en temps.
C 
        CALL MESSAO (
     &  '******************************************************
     &  \
     &  \    DEBUT DES ITERATIONS EN ESPACE
     &  \
     &  \******************************************************')
C 
        DO WHILE ((NBTOUR .LE. NBCINM) .AND.
     &            (NORACP/NORINP .GT. PREGLO))
C 
	  CALL IMPET ('ON BOUCLE EN ESPACE, NBTOUR ', NBTOUR)
          NTOUTE = 0
          NORME  = 1.D0
          NBTOUR = NBTOUR+1
          REPRIS = .FALSE.
          LONRES = LONDEP +LONEPS+LONSAU
          CALL MENADM (DEPSOL, LONRES)
C 
          CALL MESSAO (
     &  '******************************************************
     &  \
     &  \      DEBUT DE L''OPTIMISATION EN TEMPS
     &  \
     &  \******************************************************')
C 
C         On passe au moins une fois dans tous les cas
C 
          DO WHILE ((NTOUTE .EQ. 0) .OR. ((NTOUTE .LE. NITEMX)
     &               .AND. (NORME .GT. PRECIS)))
	    NTOUTE =  NTOUTE+1
	    CALL IMPET ('ON BOUCLE EN TEMPS, NTOUTE ', NTOUTE)
C 
C           Calcul de  la solution DEPSOL cinematiquement admissible
C           a zero du probleme :
C 
C           Intemps (K0)[ DEPSOL ] = Intemps (fotemp*-*DELTA(sigchap))
C 
C           DM(SIGZER) est pour le moment la contrainte Intemps (fotemp*-*DELTA(sigchap))
C 
            CALL EFFPRE (0, ADSMEG, DM(DEBFT), DM(SIGZER), DM(SGNZER),
     &                   NBDEV, M(TNUDEV), NORM1)
            CALL IMPDT ('NORME DANS ITE ETAGLO CORRESPONDANT A sigmch ',
     &                   NORM1)
            LONRES = LONEPS+LONSAU
            CALL MENADM (SIGRES, LONRES)
            T2 = TUSED()
            CALL GCADMI (NITGCA, REPRIS, ADSMEG, NBDEV, M(TNUDEV),
     &                   DM(DEPSOL), DM(EPSSOL), DM(SAUSOL), DM(SIGRES),
     &                   DM(SGNRES))
            T1 = TUSED()
            DT = (T1-T2)/50
            DTH= DT/3600
            DTS= DT-3600*DTH
C 
C           Verification que les deplacements obtenus sont bien admisssibles a zero
C 
D           CALL VDPVZE (DM(DEPSOL))
C 
C           On determine la meilleure evolution pour la solution en deplacement
C 
            CALL EVDPOP (DM(EPSSOL), DM(SAUSOL), NPICDE, NPICFI,
C                        et on recupere :
     &                   DM(DEBFT+NPICET))
C 
C           On regularise la fonction du temps obtenue
C 
CD         CALL WVISC1 (NPICET, DBLLOC, DM(ADINTE), DM(DEBFT+NPICET),
C                        et on recupere
CD   &                   DM(DEBFT+2*NPICET))
C 
C           On norme la fonction du temps obtenue
C 
            CALL COPITD (NPICET, DM(DEBFT+NPICET), DM(DEBFT+2*NPICET))
	    CALL NORMFT (ADINTE , ADPICE, DM(DEBFT+2*NPICET),
C                        et on recupere
     &                   DM(DEBFT+2*NPICET), NORFT)
C 
C           On calcule la difference des fonctions test normees precedentes
C           et actuelles. Pour savoir si le processus en temps a converge,
C           on calcule la norme de cette difference.
C 
            CALL SOU (DM(DEBFT), DM(DEBFT+2*NPICET),
C                     et on recupere
     &                DM(FOTEST), NPICET, 1)
            CALL SCALFT (ADINTE, ADPICE, DM(FOTEST), DM(FOTEST), NORME)
            CALL IMPDT ('NORME**2 DES DIFFERENCES EN TEMPS '//IDPROG,
     &                   NORME)
C 
            IF (NORME .LT. PRECIS .OR. NTOUTE .EQ. NITEMX) THEN
              CALL IMPDT  ('POUR UNE VALEUR DE PRECISION ', PRECIS)
              CALL IMPET  ('NOMBRE DE TOURS EN TEMPS     ', NTOUTE)
              CALL IMPDT  ('NORME AU CARRE DES DIFFERENCES DES '//
     &                    'EVOLUTIONS EN TEMPS < PRECIS  ', NORME)
              CALL IMPTDT ('MEILLEURE EVOLUTION DES DEPLACEMENTS ',
     &                      DM(DEBFT+NPICET), NPICET, 1)
            END IF
C 
C           On est deja passe par gcadmi pour la determination
C           de la deformation solution
C 
            REPRIS = .TRUE.
C 
C           On prend comme nouvelle fonction du temps la fonction du temps
C           optimale precedente normee
C 
            DEBFT = DEBFT+2*NPICET
C 
          END DO
C 
          CALL MESSAO (
     &  '******************************************************
     &  \
     &  \      FIN DE L''OPTIMISATION EN TEMPS
     &  \
     &  \******************************************************')
C 
          CALL COPID (DM(FOTDPP), DM(DEBFT), NPICET, 1)
C 
C         On retire a la contrainte SIGZER ayant servi au calcul de
C         la solution residuelle la contrainte residuelle de GCAMI
C         qui existe lorsque l'on n'a pas exactement converge.
C 
C         Steph. La suite est passee en commentaire. En effet, vu qu'on modifie
C         SIGZER par la suite, ce passage est inutile.
C 
          CALL SOU (DM(SIGZER), DM(SIGRES), DM(SIGZER), LONEPS, 1)
          IF (NBINT .GT. 0) THEN
            CALL SOU (DM(SGNZER), DM(SGNRES), DM(SGNZER), LONSAU, 1)
          END IF
C 
C         Calcul de la contrainte associee au deplacement DEPSOL calcule precedemment 
C         par GCADMI. Contrainte 'admissible a zero' au sens du probleme residuel :
C         K(epssol)-delta-sigchap
C 
          CALL CSIAD0 (DM(DEPSOL), DM(EPSSOL), DM(SAUSOL),
C                      on modifie
     &                 DM(SIGZER), DM(SGNZER))
C 
C         Orthogonalisation de la nouvelle fonction de l'espace, calcul de son
C         evolution optimale et remplissage de sigmch par la contrainte residuelle
C 
          T2 = TUSED()
C 
C         S'il n'y a qu'un decoupage en temps et que l'on passe par une etape
C         preliminaire
C 
          IF (UPERIO .AND. LETAPP) THEN
            CALL ETORTH (DM(DEBFT-NPICET), DM(DEPSOL), DM(EPSSOL),
     &                   DM(SAUSOL), DM(SIGZER), DM(SGNZER),
     &                   NPICDE, NPICFI)
D           CALL VORTEP (0)
          ELSE
            CALL ETORDE (DM(DEPSOL), DM(EPSSOL),
     &                   DM(SAUSOL), NPICDE, NPICFI)
            CALL ETORCO (DM(SIGZER), DM(SGNZER), NPICDE, NPICFI)
D	    CALL VORTEP (1)
          END IF
C 
C 
CD        T1 = TUSED()
CD        DT = (T1-T2)/50
CD        DTH= DT/3600
CD        DTS= DT-3600*DTH
CD        CALL MESSAO ('TEMPS UTILISE DANS ETORTH DANS '//IDPROG)
CD        CALL IMPET  ('HEURES   : ', DTH)
CD        CALL IMPET  ('SECONDES : ', DTS)
C 
C         Calcul de l'effort residuel correspondant a intemps(-delta sigmba)
C 
          CALL IMPTDT ('VERIF3 SIGZER ', DM(SIGZER), 1, 100)
          CALL EFFPRE (0, ADSMEG, DM(FCONPA), DM(SIGZER), DM(SGNZER),
     &                 NBDEV, M(TNUDEV), NORACP)
          CALL IMPDT ('NORME DANS ITE ETAGLO CORRESPONDANT A SIGCH ',
     &                 NORACP)
          CALL IMPDT ('VALEUR DE NORACP/NORINP DANS GCNLIN ',
     &                 NORACP/NORINP)
          CALL IMPDT ('VALEUR DE PREGLO DANS GCNLIN        ',
     &                 PREGLO)
        END DO
C 
        CALL MESSAO (
     &  '******************************************************
     &  \
     &  \    FIN DES ITERATIONS EN ESPACE
     &  \
     &  \******************************************************')
C 
      END DO
C 
      CALL MESSAO (
     &  '******************************************************
     &  \
     &  \  FIN DE LA BOUCLE SUR LES DIFFERENTES PERIODES EN TEMPS
     &  \
     &  \******************************************************')
C 
      CALL MESSAO (
     &  '******************************************************
     &  \
     &  \FIN DU TRAITEMENT DU PROBLEME CINEMATIQUE
     &  \
     &  \******************************************************')
C 
C     Calcul de l'effort residuel total correspondant a intemps(-delta sigmba)
C 
      CALL EFFPRE (0, ADSMEG, DM(FOCONS), DM(SIGZER), DM(SGNZER),
     &             NBDEV, M(TNUDEV), NORAC)
      CALL IMPDT ('VALEUR DE NORAC/NORIN FINAL EN DEFORMATION ',
     &              NORAC/NORIN)
C 
C     Test pour un eventuel traitement du probleme statique
C 
      IF (NBSTAM .GE. 1) THEN
C 
        CALL MESSAO (
     &  '******************************************************
     &  \
     &  \DEBUT DU TRAITEMENT DU PROBLEME STATIQUE
     &  \
     &  \******************************************************')
C 
C       Calcul de l'effort correspondant a intemps(-delta(sigmch))
C       obtenu a l'etape locale.
C 
        CALL EFFPRE (0, ADSMEG, DM(FOCONS), DM(SIGZER), DM(SGNZER),
     &               NBDEV, M(TNUDEV), NORIN)
        NORAC = NORIN
C 
C       Comme on ne comprend pas le probleme de convergence en temps
C       en statique on s'inspire du cinematique
C 
        NBTTST  = NTOUTE
C 
        CALL MESSAO (
     &  '******************************************************
     &  \
     &  \  DEBUT DE LA BOUCLE SUR LES DIFFERENTES PERIODES EN TEMPS
     &  \
     &  \******************************************************')
C  
        DO PERIOD = 1, NPERIO
          NPICDE = M(APERIO+PERIOD)
          NPICFI = M(APERIO+PERIOD +1)-1
          PREGLO = DM(AERREU+NPICFI)/2.D0
          DO I = 1, NPICDE-1
            DM(FCONPA+I-1) = 0.D0
          END DO
          DO I = NPICDE, NPICFI
            DM(FCONPA+I-1) = 1.D0
          END DO
          DO I = NPICFI+1, NPICET
            DM(FCONPA+I-1) = 0.D0
          END DO
          DO I = 1, NPICDE-1
            DM(FOTDPP+I-1) = 0.D0
          END DO
          DO I = NPICDE, NPICFI
            DM(FOTDPP+I-1) = DM(FOTDEP+I-1)
          END DO
          DO I = NPICFI+1, NPICET
            DM(FOTDPP+I-1) = 0.D0
          END DO
          CALL NORMFT (ADINTE, ADPICE, DM(FOTDPP), DM(FOTDPP), NORFT)
          CALL COPID (DM(FOTDPP), DM(DEBFT), NPICET, 1)
C 
C         Calcul de l'effort correspondant a intemps(-delta(sigmch))
C         apres l'etape preliminaire, et mesure de sa norme pour
C         l'intervalle de temps considere
C   
          CALL EFFPRE (0, ADSMEG, DM(FCONPA), DM(SIGZER), DM(SGNZER),
     &                 NBDEV, M(TNUDEV), NORINP)
          NORACP = NORINP
          NBTOUR = 1
C 
C         La version actuelle consiste a obtenir la precision souhaitee sur
C         chaque periode en temps.
C 
          CALL MESSAO (
     &  '******************************************************
     &  \
     &  \    DEBUT DES ITERATIONS EN ESPACE
     &  \
     &  \******************************************************')

          DO WHILE ((NBTOUR .LE. NBSTAM) .AND.
     &              (NORACP/NORINP .GT. PREGLO))
	    CALL IMPET ('ON BOUCLE EN ESPACE, NBTOUR = ', NBTOUR)
	    NTOUTE = 0
            NORME  = 1.D0
            NBTOUR = NBTOUR+1
            REPRIS = .FALSE.
            LONRES = LONDEP +LONEPS+LONSAU
            CALL MENADM (DEPSOL, LONRES)
C 
            CALL MESSAO (
     &  '******************************************************
     &  \
     &  \      DEBUT DE L''OPTIMISATION EN TEMPS
     &  \
     &  \******************************************************')
C 
C           On passe quand meme un fois
C 
            DO WHILE ((NTOUTE .EQ. 0) .OR. ((NTOUTE .LE. NBTTST)
     &                 .AND. (NORME .GT. PRECIS)))
              NTOUTE =  NTOUTE+1
	      CALL MESSAO (' ')
	      CALL IMPET ('ON BOUCLE EN TEMPS, NTOUTE = ', NTOUTE)
              CALL MESSAO (' ')
C 
C             Calcul de la solution DEPSOL cinematiquement admissible a zero du probleme :
C 
C             Intemps (K0)[DEPSOL] = Intemps (fotemp*-*DELTA(sigchap))
C 
C             DM(SIGZER) est pour le moment la contrainte Intemps (fotemp*-*DELTA(sigchap))
C 
              CALL EFFPRE (0, ADSMEG, DM(DEBFT), DM(SIGZER), DM(SGNZER),
     &                     NBDEV, M(TNUDEV), NORM1)
              LONRES = LONEPS+LONSAU
              CALL MENADM (SIGRES, LONRES)
              T2 = TUSED()
              CALL GCADMI (NITGCA, REPRIS, ADSMEG, NBDEV, M(TNUDEV),
     &                     DM(DEPSOL), DM(EPSSOL), DM(SAUSOL),
     &                     DM(SIGRES), DM(SGNRES))
              T1 = TUSED()
              DT = (T1-T2)/50
              DTH= DT/3600
              DTS= DT-3600*DTH
C 
C             Verification que les deplacements obtenus sont bien admisssible a zero
C 
D             CALL VDPVZE (DM(DEPSOL))
C 
C             Calcul  de la contrainte admissible a zero au sens du probleme RESIDUEL :
C             K(epssol)-delta-sigchap
C 
              CALL CSIAD0 (DM(DEPSOL), DM(EPSSOL), DM(SAUSOL),
C                          on modifie
     &                     DM(SIGZER), DM(SGNZER))
C 
C             On determine la meilleure evolution pour la solution en contrainte
C 
              CALL EVSIOP (DM(SIGZER), DM(SGNZER), NPICDE, NPICFI,
C                          on recupere
     &                     DM(DEBFT+NPICET))
C 
C             On modifie les fonctions du temps en fonction 
C 
CD            CALL WVISC1 (NPICET, DBLLOC, DM(ADINTE), DM(DEBFT+NPICET),
C                          on recupere
CD   &                     DM(DEBFT+2*NPICET))
C 
              CALL COPITD (NPICET, DM(DEBFT+NPICET), DM(DEBFT+2*NPICET))
	      CALL NORMFT (ADINTE, ADPICE, DM(DEBFT+2*NPICET),
C                          on recupere
     &                     DM(DEBFT+2*NPICET), NORFT)
C 
C             On calcule la difference des fonctions test normees precedentes
C             et actuelles. Pour savoir si le processus en temps a converge,
C             on calcule la norme de cette difference.
C 
C              CALL IMPTDT ('FONCTION DEPART POUR LE TEST '//IDPROG,
C     &                       DM(DEBFT),NPICET,1)
C 
C              CALL IMPTDT ('FONCTION ARRIVEE  POUR LE TEST '//IDPROG,
C    &                       DM(DEBFT+2*NPICET), NPICET, 1)
C 
              CALL SOU (DM(DEBFT), DM(DEBFT+2*NPICET),
C                       on recupere
     &                  DM(FOTEST), NPICET, 1)
              CALL SCALFT (ADINTE, ADPICE, DM(FOTEST),
     &                     DM(FOTEST), NORME)
              CALL IMPDT ('NORME AU CARRE DES DIFFERENCES EN TEMPS ',
     &                     NORME)
              IF (NORME .LT. PRECIS .OR . NTOUTE .EQ. NITEMX) THEN
                CALL IMPDT ('POUR UNE VALEUR DE PRECISION ', PRECIS)
                CALL IMPET ('NOMBRE DE TOURS EN TEMPS     ', NTOUTE)
                CALL IMPDT ('NORME AU CARRE DES DIFFERENCES DES '//
     &                       'EVOLUTIONS EN TEMPS < PRECIS ', NORME)
                CALL IMPTDT ('MEILLEURE EVOLUTION DES DEPLACEMENTS ',
     &                         DM(DEBFT+NPICET), NPICET, 1)
              END IF
C 
C             On est deja passe par gcadmi pour la determination
C             de la deformation solution.
C 
              REPRIS = .TRUE.
C 
C             On prend comme nouvelle fonction du temps la fonction du temps
C             optimale precedente normee.
C 
              DEBFT = DEBFT+2*NPICET
C 
            END DO
C 
            CALL MESSAO (
     &  '******************************************************
     &  \
     &  \      FIN DE L''OPTIMISATION EN TEMPS
     &  \
     &  \******************************************************')
C   
            CALL COPID (DM(FOTDPP), DM(DEBFT), NPICET, 1)
C 
C           Orthogonalisation de la nouvelle fonction de l'espace.
C           Calcul de son evolution optimale et remplissage de sigmch
C           par la contrainte residuelle.
C 
            T2 = TUSED()
C 
            CALL ETORDE (DM(DEPSOL), DM(EPSSOL),
     &                   DM(SAUSOL), NPICDE, NPICFI)
            CALL ETORCO (DM(SIGZER), DM(SGNZER), NPICDE, NPICFI)
C 
            T1 = TUSED()
            DT = (T1-T2)/50
            DTH= DT/3600
            DTS= DT-3600*DTH
C 
C           Calcul de l'effort residuel correspondant a intemps(-delta sigmba )
C 
            CALL EFFPRE (0, ADSMEG, DM(FCONPA), DM(SIGZER), DM(SGNZER),
     &                   NBDEV, M(TNUDEV), NORACP)
            CALL IMPDT ('VALEUR DE NORACP/NORINP DANS GCNLIN ',
     &                    NORACP/NORINP )
            CALL IMPDT ('VALEUR DE PREGLO DANS GCNLIN ',
     &                    PREGLO )
          END DO
C 
          CALL MESSAO (
     &  '******************************************************
     &  \
     &  \    FIN DES ITERATIONS EN ESPACE
     &  \
     &  \******************************************************')
C 
        END DO
C 
        CALL MESSAO (
     &  '******************************************************
     &  \
     &  \  FIN DE LA BOUCLE SUR LES DIFFERENTES PERIODES EN TEMPS
     &  \
     &  \******************************************************')
C 
C       Calcul de l'effort residuel total correspondant a intemps(-delta sigmba)
C 
        CALL EFFPRE (0, ADSMEG, DM(FOCONS), DM(SIGZER), DM(SGNZER),
     &               NBDEV, M(TNUDEV), NORAC)
        CALL IMPDT ('VALEUR DE NORAC/NORIN FINAL EN CONTRAINTE ',
     &               NORAC/NORIN)
C 
        CALL MESSAO (
     &  '******************************************************
     &  \
     &  \FIN DU TRAITEMENT DU PROBLEME STATIQUE
     &  \
     &  \******************************************************')
C 
C     Fin du test sur un eventuel traitement du probleme statique
C 
      END IF
C 
C     Ecriture des fonctions du temps admissibles de l'etape globale
C 
      CALL ECFTGL (DM(FTREE7), DM(FTRE10))
      CALL DTEXPR (0, NBETGL)
C 
      CALL MESSAO (
     &  '******************************************************
     &  \
     &  \FIN DU TRAITEMENT DE L''ETAPE GLOBALE NUMERO '//CARETG //'
     &  \
     &  \******************************************************')
C 
      T1 = TUSED()
      DT = (T1-T0)/50
      DTH= DT/3600
      DTS= DT-3600*DTH
      CALL MESSAO ('TEMPS UTILISE EN TOUT DANS '//IDPROG)
      CALL IMPET  ('TEMPS UTILISE EN HEURES    ', DTH )
      CALL IMPET  ('TEMPS UTILISE EN SECONDES  ', DTS )
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine determine les efforts globaux qui correspondent a
C     l'integration en temps du contenu de sigmch, c'est a dire de
C     [delta_sigma_chapeau - fi(t)* sigi] (construit a l'etape globale
C     actuelle) integre en temps avec FOTEMP; elle lit le fichier direct
C     sigmch rempli a l'etape locale, donc elle doit etre utilisee apres
C     l'etape locale.
C 
C     On procede en au moins deux etapes :
C 
C       - determination des differents efforts pour toutes les valeurs de teta :
C                               _         .
C       - valeurs de teta  <=> -B delta(sign)
C       - developpement en serie de Fourier de ces efforts
C 
C     On envoie comme arguments :
C 
C     E ...... TYPE    = 0 => Efforts associes a -delta(sigmch)
C                      = 1 => Efforts associes a -delta(epsich)
C     E ...... ADSMEG  l'adresse de depart des seconds membres pour l'etape globale
C                      ceux-ci sont assembles dans DM a partir de ADSMEG (NDDL, NBMAT)
C     E....... FOTEMP  valeurs aux differents piquets de temps de l'accroissement de
C                      fonction du temps telle que :
C                      F = intemp(fotemp*-B delta(sign))
C     E....... SIGSSM  le residu en contraintes dans les couches
C     E....... SGNSSM  le residu en contraintes dans les interfaces
C 
C     Et on recupere :
C 
C     S ...... NORME   effort residuel total correspondant a intemps(-delta sigmba)
C 
C     Les seconds membres sont mis a zero en entree <=> dm(adsmeg, ...) = 0.d0  !!!!
C 
      SUBROUTINE EFFPRE (TYPE, ADSMEG, FOTEMP, SIGSSM, SGNSSM,
     &                   NBDEV, TNUDEV, NORME)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER          TYPE, ADSMEG, NBDEV, TNUDEV(NBMAT)
C 
      DOUBLE PRECISION FOTEMP(NPICET), NORME
      DOUBLE PRECISION SIGSSM(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION SGNSSM(NEPS*NTETA*NGAU1)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TLOCN1, DEBGAU, DBGAUI, ADINTE, DBDDLU, DBDDLV, DBDDLW
      INTEGER  DBDDIU, DBDDIV, DBDDIW, BSPCAC, ASIGIN, TIGINT, LONRES 
      INTEGER  DSPCAC, DSPCAI, LONECF, BSPCAI, ASININ, TININT, NFT
      INTEGER  NUCOU, NUCOL, ADRGAU, X, Y, SIGINT
      INTEGER  NUINT, SININT, TETA, H, K, ADTETA, I, ADCOU, ADINT
      INTEGER  AM2LC, ADM2LC
C 
      LOGICAL PASSP
C 
      DOUBLE PRECISION  A, B, R, RAYONC, POIDG, WI
      DOUBLE PRECISION  KLOC(17), COUISO(10)
      DOUBLE PRECISION  KLOCI(9), INTISO(3)
      DOUBLE PRECISION  TETORT, TETCAL
      DOUBLE PRECISION  EPSLOC (6) , SAULOC(3)
C 
      LOGICAL LTYP1
C 
C     POUR LE FICHIER SIGMCH
C 
      INTEGER       NUENSI, LONSCH, IUNSIG, DEBSIG
      CHARACTER*6   IDPROG
      CHARACTER*3   BARATE
      CHARACTER*20  NOMFIC, NOM
C 
      PARAMETER (IDPROG='EFFPRE')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      LONRES = NDDL*NBMAT
      CALL MENADM (ADSMEG, LONRES)
      NFT = NFTGLO
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
C 
      K=XINTEG*(XINTEG-1)/2
      H=YINTEG*(YINTEG-1)/2
C 
C     NOMFIC est le nom de fichier dans lequel sont stockees les quantites barre
C     de l'etape locale precedente. Ouverture du fichier pour les contraintes.
C 
      LTYP1 = .FALSE.
C 
      IF (TYPE .EQ. 0) THEN
	NOM = 'sigmch'
      ELSE IF (TYPE .EQ. 1) THEN
        LTYP1 = .TRUE.
        NOM = 'epsich'
        CALL ADTBDM ('HOO-COUCHE', ADCOU)
        IF (NBINT .GT. 0) THEN
          CALL ADTBDM ('HOO-INTERF', ADINT)
        ELSE
          ADINT = ADCOU
        END IF
      END IF
C 
      NOMFIC = NOM
      LONSCH = 6*NPICET
      NUENSI  = 1
      CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C 
C     Recherche des differents tableaux utiles
C     TLOCN1      pour le rangement des numeros de noeuds
C     TAB-GAUSS   pour les fonctions elementaires
C     INTE-TEMPS  pour les intervalles de temps
C 
      CALL ADTBM  ('TLOCN1    ', TLOCN1)
      CALL ADTBDM ('TAB-GAUSS ', DEBGAU)
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
C 
C     Ouverture d'un tableau partiel pour ranger les numeros
C     des ddl des deplacements ranges calcul (u, v, w) par developpement croissant
C 
      CALL POUSME (NDDLEL, DBDDLU)
      DBDDLV = DBDDLU+12
      DBDDLW = DBDDLV+12
C 
C     Pour les interfaces
C 
      DBDDIU = DBDDLU
      DBDDIV = DBDDIU+8
      DBDDIW = DBDDIV+8
C 
C     Ouverture d'un tableau partiel pour ranger les
C     quantites stockees en fin de boucle sur le temps (etape-locale)
C     et les contraintes * GT integrees sur le temps pour tous les teta
C 
C     DSPCAC delta sigma point chapeau actuel XDT                  X 6*NPICET
C     ASIGIN sigma*GT integre sur le temps puis sigma developpe    X 6*(NTETA+1)
C     TIGINT transpose de sigma*gt                                 X 6*NTETA
C                                                              ____________________
C                                                               6*NPICET+12*NTETA+6
      LONRES = 6*NPICET+12*NTETA+6
C 
      CALL POUSMD (LONRES, DSPCAC)
      BSPCAC = DSPCAC
C 
C     Pour les quantites de l'etape locale pour les couches
C 
      ASIGIN = DSPCAC+6*NPICET
      TIGINT = ASIGIN+6*(NTETA+1)
      LONECF = 6*NPICET
      DEBSIG = 1
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES COUCHES
C 
      DO NUCOU= 1, NBCOU
C 
C     Recherche des caracteristiques geometriques de la couche
C 
        CALL VALEP (NUCOU, B)
        IF (LTYP1) THEN
          CALL ANGCOU (NUCOU, TETORT)
          CALL COMISO (ADCOU, NUCOU, COUISO)
        END IF
C 
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1, NBCOL
C 
          CALL VALRAY (NUCOL, RAYONC, A)
C 
C         Recherche des numeros de ddl ranges calcul pour l'element
C 
          CALL DDLCAL (1, NUCOU, NUCOL, TLOCN1, M(DBDDLU))
          CALL DDLCAL (2, NUCOU, NUCOL, TLOCN1, M(DBDDLV))
          CALL DDLCAL (3, NUCOU, NUCOL, TLOCN1, M(DBDDLW))

C         DEBGAU est l'adresse pour x=1, y=1 dans TAB-GAUSS
C 
          ADRGAU = DEBGAU
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT X
C 
          DO X = 1, XINTEG
C 
            WI=POIDS(K+X)
C 
C           Calcul du rayon au point de gauss
C 
            R = RAYONC + A*GAUSS(XINTEG*(XINTEG-1)/2+X)
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Y
C 
            DO Y = 1, YINTEG
C 
              POIDG = WI*POIDS(H+Y)
C 
C             BOUCLE iv SUR LES ANGLES
C 
              SIGINT = ASIGIN
C 
              DO TETA = 1, NTETA
C 
C               On recupere les quantites stockees dans le directory Q-CHAPEAU
C               pour le point de gauss considere pour tous les temps c'est a dire :
C 
C               DSPCAC delta sigma point chapeau actuel XDT    X 6*NPICET
C 
C               Toutes ces quantites sont multiplies par interv
C 
C               Lecture directe du fichier (dans q-chapeau) des accroissements
C               admissibles de l'etape locale precedente.
C 
                CALL LFDDNF (DM(DSPCAC), DSPCAC, LONSCH, IUNSIG, NUENSI)
		NUENSI = NUENSI +1
C 
C               Equivalent a une boucle sur les piquets de temps pour integrer (sigma DT)*GT(t)
C               LTYP1 pour une sauvegarde en deformations a l'etape locale
C                     pour une sauvegarde en contraintes sinon
C 
                IF (LTYP1) THEN
                  TETCAL = DM (ADTETA+TETA-1)
C 
C                 RICLOC : Rigidite couche dans la base locale
C 
                  CALL RICLOC (COUISO, TETORT, TETCAL, KLOC)
                  CALL STSIPA (DM(ADINTE), 6, DM(DSPCAC),
     &                         FOTEMP, EPSLOC)
                  DO I = 1, 6
                    EPSLOC(I) = -EPSLOC(I)
                  END DO
                  CALL MULORT (KLOC(1), KLOC(10), KLOC(14), KLOC(17),
     &                         EPSLOC, DM(SIGINT))
                ELSE
                  CALL STSIPA (DM(ADINTE), 6, DM(DSPCAC),
     &                         FOTEMP, DM(SIGINT))
                END IF
                SIGINT = SIGINT+6

C             FIN DE BOUCLE iv SUR LES ANGLES
C 
              END DO
C 
              CALL COPITD (6*NTETA, DM(ASIGIN), SIGSSM(DEBSIG))
              DEBSIG = DEBSIG+ 6*NTETA
C 
C             Rangement des contraintes integrees sur le temps
C             sous la forme (nteta, 6)
C 
              CALL TRANSP (DM(ASIGIN), DM(TIGINT), 6, NTETA)
C 
C             Calcul des valeurs developpes rangees calcul des contraintes
C             integrees sur le temps rangees dans ASIGIN
C 
              CALL SIGDEV (DM(TIGINT), DM(ASIGIN))
              CALL SMITEC (POIDG, NFT, DM(ADRGAU), A, B, R,
     &                     M(DBDDLU), DM(ASIGIN), ADSMEG)
C 
C             Pour aller lire dans TAB_GAUSS au bon point de GAUSS
C 
              ADRGAU      = ADRGAU+36
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Z
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES COLONNES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES COUCHES
C 
      END DO
C 
C     Fermeture du fichier
C 
      CALL FERFIC (3, IUNSIG, IDPROG)
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
C     TEST SUR SUR LES INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
C 
C       Ouverture d'un tableau partiel pour ranger les
C       quantites stockees en fin de boucle sur le temps (etape-locale)
C       et les contrainte * GT integrees sur le temps pour tous les teta
C 
C       DSPCAI delta sigma point chapeau actuel XDT                 X 3*NPICET
C       ASININ sigma*GT integre sur le temps puis sigma developpe   X 3*(NTETA+1)
C       TININT transpose de sigma*gt                                X 3*NTETA
C                                                                ____________________
C                                                                 3*NPICET+6*NTETA+3
C       Pour les quantites de l'etape locale pour les interfaces
C 
        LONRES =  3*NPICET+6*NTETA+3
        CALL POUSMD (LONRES, DSPCAI)
        BSPCAI = DSPCAI
        ASININ = BSPCAI+3*NPICET
        TININT = ASININ+3*(NTETA+1)
C 
C       OUVERTURE du fichier pour les contraintes normales
C 
        NOM = 'sinoch'
        IF (TYPE  .EQ. 0) THEN
          NOM = 'sinoch'
        ELSE IF (TYPE  .EQ. 1) THEN
          NOM = 'sautch'
        END IF
C 
        NOMFIC = NOM
        NUENSI = 1
        LONSCH = 3*NPICET
        CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
        CALL ADTBDM ('GAUSSINTER', DBGAUI)
C 
        DEBSIG = 1
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT = 1, NBINT
C 
C         SYMPAR : .TRUE. si le plan de symetrie est une interface
C 
          PASSP = .FALSE.
          IF (NUINT .EQ. 1 .AND. SYMPAR) PASSP = .TRUE.
          IF (LTYP1) THEN
            CALL IOMISO (ADINT, NUINT, INTISO)
            CALL ANGINT (NUINT, TETORT)
          END IF
C 
C         BOUCLE i SUR LES COLONNES
C 
          DO NUCOL = 1, NBCOL
C 
            CALL VALRAY (NUCOL, RAYONC, A)
C 
C           Recherche des numeros de ddl ranges calcul pour l'element
C 
            CALL DDLICA (NUINT, NUCOL, TLOCN1, M(DBDDIU),
     &                   M(DBDDIV), M(DBDDIW))

C           DEBGAU est l'adresse pour x=1, y=1 dans TAB-GAUSS
C 
            ADRGAU = DBGAUI
C 
C           BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
C 
            DO X = 1, XINTEG
C 
              POIDG =POIDS(K+X)
C 
C             Calcul du rayon au point de gauss
C 
              R  = RAYONC + A*GAUSS( XINTEG*(XINTEG-1)/2 +X )
C 
C             BOUCLE iii SUR LES ANGLES
C 
              SININT = ASININ
C 
              DO TETA = 1, NTETA
C 
C               On recupere les quantites stockees dans la directory  Q-CHAPEAU
C               pour le point de gauss considere pour tout les temps c'est a dire :
C 
C               DEPCAI delta epsilon point chapeau actuel XDT    X 3*NPICET
C               DSPCAI - delta sigma point chapeau actuel XDT    X 3*NPICET
C               CMTGAI comportement tangent actuel        XDT    X 9*NPICET
C 
C               Toutes ces quantites sont multipliees par interv
C               Lecture sequentielle du fichier dans Q-chapeau :
C 
C               Lecture directe du fichier (dans q-chapeau) des accroissements
C               admissibles de l'etape locale precedente.
C 
                CALL LFDDNF (DM(DSPCAI), DSPCAI, LONSCH, IUNSIG, NUENSI)
                NUENSI = NUENSI +1
C 
C               Equivalent a une boucle sur les piquets de temps pour integrer (sigma DT) * GT(t)
C 
                IF (LTYP1) THEN
C 
C               Recherche de l'angle de la bande (tetcal) correspondant a teta
C 
                  TETCAL = DM(ADTETA+TETA-1)
                  CALL INTELA (INTISO, TETORT, TETCAL, KLOCI)
                  CALL STSIPA (DM(ADINTE), 3,
     &                         DM(DSPCAI), FOTEMP, SAULOC)
                  DO I = 1, 3
                    SAULOC(I) = -SAULOC(I)
                  END DO
                  CALL IULORT (KLOCI(1), SAULOC, DM(SININT))
                ELSE
                  CALL STSIPA (DM(ADINTE), 3,
     &                         DM(DSPCAI), FOTEMP, DM(SININT))
                END IF
C 
                SININT = SININT+3
C 
C             FIN DE BOUCLE iii SUR LES ANGLES
C 
              END DO
C 
              CALL COPITD (3*NTETA, DM(ASININ), SGNSSM(DEBSIG))
              DEBSIG = DEBSIG+ 3*NTETA
C 
C             Rangement des contraintes normales integrees sur le temps (nteta, 3)
C 
              CALL TRANSP (DM(ASININ), DM(TININT), 3, NTETA)

C 
C             Calcul des valeurs developpes rangees calcul des contraintes
C             normales integrees sur le temps rangees dans ASIGIN
C 
              CALL SINDEV (DM(TININT), DM(ASININ))
              IF (PASSP) THEN
                CALL SMITEP (POIDG, NFT, DM(ADRGAU), A,
     &                       R, M(DBDDIU), DM(ASININ), ADSMEG)
              ELSE
                CALL SMITEI (POIDG, NFT, DM(ADRGAU), A,
     &                       R, M(DBDDIU), DM(ASININ), ADSMEG)
              END IF
C 
C             Pour aller lire dans TAB_GAUSS au bon point de GAUSS
C 
              ADRGAU = ADRGAU+8
C 
C           FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES COLONNES
C 
          END DO
C 
C       FIN DE BOUCLE SUR LES INTERFACES
C 
        END DO
C 
C       Fermeture du fichier
C 
        CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     FIN DU TEST SUR LES INTERFACES
C 
      END IF
C 
C     Mise a zero des termes non significatifs des efforts
C 
      CALL MZBLOC (1, DM(ADSMEG), NBDEV, TNUDEV, NORME)
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule la meilleure evolution des deplacements
C     et la meilleure evolution des contraintes sur les quantites
C     admissibles a zero construites precedemment.
C 
C     On modifie en fin de routine le fichier des DELTA SIGMA BARRE.
C 
      SUBROUTINE ETAPRE
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
      INTEGER AM2LC, ADM2LC
      INTEGER FTSCH5, FTECH6
      INTEGER EPSCH8, SGNC12, SIGC11, SAUCH9, FTRE10
      INTEGER DEPCH0, EPSTRA, SAUTRA, DEPTRA
      INTEGER NUEPS, NDEBEP, NDEBSI, NUSIG
      INTEGER DEPSTR, DSAUTR, DDEPTR, DSIGTR, DSGNTR
      INTEGER FOTDEP, I, FOTSIG
      INTEGER DEBTDE, DEBTSI, LONDER, NBFONC
      INTEGER LONEPS, LONSAU
      INTEGER SIGTRA, SGNTRA, FTREE7, COEFF, COESIG
      INTEGER ADINT, ADCOU, ADSOU, ADSNT
      INTEGER TABNIE(2), TABNIS(2), TABNID(2)
C 
C     LONGUEURS DES TABLEAUX DE TRAVAIL
C 
      INTEGER  TRAV1, TRAV2, ADINTE
C 
      CHARACTER*6 IDPROG
      CHARACTER*3 BARATE
      PARAMETER (IDPROG='ETAPRE')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C 
C     Appel des tableaux pour les rigidite et souplesse des couches
C 
      CALL ADTBDM ('HOO-COUCHE', ADCOU)
      CALL ADTBDM ('SOU-COUCHE', ADSOU)
C 
C     Appel des tableaux pour les rigidite et souplesse des interfaces
C 
      IF (NBINT .GT. 0) THEN
        CALL ADTBDM ('HOO-INTERF', ADINT)
        CALL ADTBDM ('SOU-INTERF', ADSNT)
      ELSE
        ADINT = ADCOU
        ADSNT = ADSOU
      END IF
C 
C     Appel du tableau des intervalles entre piquets de temps
C 
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
C 
C     Appel des tableaux des deplacements, deformations et contraintes
C     admissibles totales pour les couches
C  
      CALL ADTBDM ('DEPLA-ADMI', DEPCH0)
      CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
      CALL ADTBDM ('SIG-AD-TOT', SIGC11)
C 
C     Appel des tableaux des deformations et contraintes admissibles totales
C     pour les interfaces
C 
      IF (NBINT .GT. 0) THEN
        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
      ELSE
        SAUCH9 = EPSCH8
        SGNC12 = SIGC11
      ENDIF
C 
C     Appel des tableaux des fonctions du temps pour l'evolution des valeurs
C     admissibles en contraintes, deformations; destines a etre remis a zero
C     puis remplis a l'etape preliminaire.
C 
      CALL IMPET ('NOMBRE DE F(TPS) DELTA CAD0 DANS '//IDPROG, NBFEPS)
      CALL IMPET ('NOMBRE DE F(TPS) DELTA SAD0 DANS '//IDPROG, NBFSIG)
      CALL ADTBDM ('TEMPS-SIGM', FTSCH5)
      CALL ADTBDM ('TEMPS-EPSI', FTECH6)
C 
C     Appel des tableaux des fonctions du temps pour les valeurs admissibles
C     en contraintes, deformations.
C 
      CALL IMPET ('NOMBRE DE F(TPS) TOTALE CA DANS '//IDPROG, DEADTR)
      CALL IMPET ('NOMBRE DE F(TPS) TOTALE SA DANS '//IDPROG, COTORE)           
      CALL ADTBDM ('TEMPS-REEL', FTREE7)
      CALL ADTBDM ('TEMP-SI-RE', FTRE10)
C 
      LONDER = NDDL*NTETA
      LONEPS = NEPS*NGAU1*NTETA
      LONSAU = NSAU*NGAU2*NTETA
C 
C -----------------------------------------------------------------------
C 
C     STOCKAGE SEQUENTIEL DES NBFEPS DEFORMATIONS
C 
C -----------------------------------------------------------------------
C 
C     Ceci est necessaire pour calculer les meilleures evolutions
C     par l'intermediaire de schast.
C     On utilise toutes les fonctions de l'espace admissibles a zero et
C     stockees de 2 a nbfeps < ou = a charax.
C     NDEBEP est le numero du premier champ servant au calcul de
C     l'accroissement admissible a zero en deformation calcule a l'etape
C     globale precedente.
C 
C     On extrait du tableau des deformations totales stockees
C     (charax, loneps) les deformations allant du numero NDEBEP a DEADTR;
C     normalement on extrait NFEPS deformations "CA a 0".
C 
      NDEBEP = DEADTR-NBFEPS+1
      TRAV1  = NBFEPS*(LONEPS+LONSAU+LONDER)+ NBFEPS
      CALL POUSMD (TRAV1, COEFF)
      EPSTRA    = COEFF  + NBFEPS
      SAUTRA    = EPSTRA + NBFEPS*LONEPS
      DEPTRA    = SAUTRA + NBFEPS*LONSAU
      DEPSTR    = EPSTRA
      DSAUTR    = SAUTRA
      DDEPTR    = DEPTRA
      TABNIE(1) = CHARAX
      TABNIE(2) = LONEPS
      TABNIS(1) = CHARAX
      TABNIS(2) = LONSAU
      TABNID(1) = CHARAX
      TABNID(2) = LONDER
C 
      DO NUEPS = NDEBEP, DEADTR
        CALL EXTRAD (DM(EPSCH8), 2, TABNIE(1), 2,
     &               NUEPS, DM(DEPSTR), LONEPS)
        DEPSTR = DEPSTR+LONEPS
        CALL EXTRAD (DM(DEPCH0), 2, TABNID(1), 2,
     &               NUEPS, DM(DDEPTR), LONDER)
        DDEPTR = DDEPTR+LONDER
        IF (NBINT .GT. 0) THEN
          CALL EXTRAD (DM(SAUCH9), 2, TABNIS(1), 2,
     &                 NUEPS, DM(DSAUTR), LONSAU)
          DSAUTR = DSAUTR+LONSAU
        END IF
      END DO
C 
C -----------------------------------------------------------------------
C 
C     Calcul de intesp(-delta(sigcha* EPSILONI) pour la valeur
C     de delta(sigcha) correspondant a l'etape locale precedente.
C     Cette evolution est aussi, puisque les deformations sont
C     orthonormalisees au sens de K0, la meilleure evolution des
C     deplacements fi(t) pour le delta admissible de l'etape globale actuelle.
C 
C -----------------------------------------------------------------------
C 
      CALL POUSMD (NBFEPS*NPICET, FOTDEP)
      CALL SCHAST (0, 'sigmch', 'sinoch',  NBFEPS,
     &             DM(EPSTRA), DM(SAUTRA), 1, NPICET, DM(FOTDEP))
C 
C     Regularisation des evolutions des deplacements
C 
      DEBTDE = FOTDEP
D     CALL IMPTDT ('MEILLEURE EVOLUTION DES DEPLACEMENTS ',
D    &              DM(FOTDEP), NPICET, NBFEPS)
C 
C -----------------------------------------------------------------------
C 
C     Modification des fonctions du temps champs pour tenir compte
C     de l'etape preliminaire pour pouvoir calculer la valeur reelle
C     qui est la somme de la valeur reelle precedente et de la
C     contribution de l'etape preliminaire
C 
C -----------------------------------------------------------------------
C 
      CALL ADDCHA (NPICET, NBDPTR, NBFEPS, DM(FOTDEP), DM(FTREE7))
C 
C     Pour aller ranger les evolutions au bon endroit on remet
C     NBFEPS a zero puis on range les nouvelles evolutions des
C     deplacements et deformations.
C 
      NBFONC = NBFEPS
      NBFEPS = 0
      DEBTDE = FOTDEP
      DO I = 1, NBFONC
        CALL RCHAMP (6, NPICET, DM(DEBTDE), DM(FTECH6))
        DEBTDE = DEBTDE+NPICET
      END DO
C 
C     Appel a la routine remplisssant le fichier scb-numero d'etape locale
C     en y rangeant - delta(sigcha) - K (delta-tild(epsn))
C 
C     Les arguments de RSICBA sont :
C 
C     type 0 ou 1
C     adcou, adint   => debut du tableau K0
C     adsou, adsnt   => debut du tableau K0-1
C     1              => numero de la 1ere evolution                (chamax)
C     NBFEPS         => numero de la derniere evolution            (chamax)
C     NDEBEP         => numero de la 1ere deformation "CA a 0"     (charax)
C     DEADTR         => numero de la derniere deformation "CA a 0" (charax)
C     EPSCH8, SAUCH9 => EPS-AD-TOT
C     FTECH6         => TEMPS-EPSI
C     NPICDE         => piquet de temps de debut d'integration
C     NPICFI         => piquet de temps de fin d'integration
C 
C 
      CALL RSICBA (0, ADCOU, ADINT, ADSOU, ADSNT,
     &             'sigmch', 'sinoch',
     &             1, NBFEPS, NDEBEP, DEADTR,
     &             EPSCH8, SAUCH9, FTECH6, 1, NPICET)
C 
C     TEST SUR LA COHERENCE DE SCHAST ET DE RSICBA
C 
D     CALL SCHAST (0, 'sigmch', 'sinoch', NBFEPS,
D    &             DM(EPSTRA), DM(SAUTRA), 1, NPICET,
D    &             DM(FOTDEP))
D     CALL IMPTDT ('MEILLEURES EVOLUTIONS DES DEFORMATIONS '//
D    &             'APRES UN PREMIER PASSAGE NULLES? '//IDPROG,
D    &              DM(FOTDEP), NPICET, NBFEPS)
C 
C -----------------------------------------------------------------------
C 
C     STOCKAGE SEQUENTIEL DES NBFSIG CONTRAINTES
C 
C -----------------------------------------------------------------------
C 
C     NDEBSI est le numero du premier champ servant au calcul de
C     l'accroissement admissible a zero en contrainte calcule a l'etape
C     globale precedente.
C 
C     On orthogonalise le dernier champ cree par rapport a tous les
C     champs servant a l'orthogonalisation qui sont au nombre de
C     NBFEPS -1 pour les deformations.
C 
C     On extrait du tableau des contraintes totales stockes (nbchare, loneps)
C     les deformations allant du numero NDEBSI a DEADTR
C 
      NDEBSI = COTORE-NBFSIG+1
      TRAV2  =  NBFSIG*(LONEPS+LONSAU)+ NBFSIG
      CALL POUSMD (TRAV2, COESIG)
      SIGTRA = COESIG+NBFSIG
      DSIGTR = SIGTRA
      SGNTRA = SIGTRA+NBFSIG*LONEPS
      DSGNTR = SGNTRA
C 
      DO NUSIG = NDEBSI, COTORE
C 
        CALL EXTRAD (DM(SIGC11), 2, TABNIE(1),
     &               2, NUSIG, DM(DSIGTR), LONEPS)
        DSIGTR = DSIGTR+LONEPS
        IF (NBINT .GT. 0) THEN
          CALL EXTRAD (DM(SGNC12), 2, TABNIS(1),
     &                 2, NUSIG, DM(DSGNTR), LONSAU)
          DSGNTR = DSGNTR+LONSAU
        END IF
C 
      END DO
C 
      TRAV1  = NBFSIG*NPICET
      CALL POUSMD (TRAV1, FOTSIG)
C 
C     Calcul de INTESP(DELTA(EPSADMI)* EPSCHA) pour la valeur
C     de delta(epscha) correspondant a l'etape locale
C 
      CALL SCHAST (2, 'sigmch', 'sinoch', NBFSIG,
     &             DM(SIGTRA), DM(SGNTRA),
     &             1, NPICET, DM(FOTSIG))
      CALL MUMARE (-1.D0, NPICET*NBFSIG, DM(FOTSIG), DM(FOTSIG))
C 
C     Regularisation des evolutions des contraintes
C 
D     CALL IMPTDT ('MEILLEURE EVOLUTION DES CONTRAINTES '//IDPROG,
D    &              DM(FOTSIG), NPICET, NBFSIG)
C 
C -----------------------------------------------------------------------
C 
C     Modification des fonctions du temps champs pour tenir compte
C     de l'etape preliminaire pour pouvoir calculer la valeur reelle
C     qui est la somme de la valeur reelle precedente et de la
C     contribution de l'etape preliminaire.
C 
C -----------------------------------------------------------------------
C 
      CALL ADDCHA (NPICET, EVCOTR, NBFSIG, DM(FOTSIG), DM(FTRE10))
C 
C -----------------------------------------------------------------------
C 
C     Pour aller ranger les evolutions au bon endroit on remet NBFSIG
C     a un, puis on range les nouvelles evolutions des contraintes.
C 
C -----------------------------------------------------------------------
C 
      NBFONC = NBFSIG
      NBFSIG = 0
      DEBTSI = FOTSIG
      DO I = 1, NBFONC
        CALL RCHAMP (5, NPICET, DM(DEBTSI), DM(FTSCH5))
        DEBTSI = DEBTSI+NPICET
      END DO
C 
C     Appel a la routine remplissant le fichier ecbarr
C     en y rangeant (K-1)(delta(sigcha)-(delta-tild(sign))
C 
      CALL RSICBA (2, ADSOU, ADSNT, ADCOU, ADINT,
     &             'sigmch', 'sinoch',
     &             1, NBFSIG, NDEBSI, COTORE,
     &             SIGC11, SGNC12, FTSCH5, 1, NPICET)
C 
C     TEST SUR LA COHERENCE DE SCHAST ET DE RSICBA
C 
D     CALL SCHAST (2, 'sigmch', 'sinoch', NBFSIG,
D    &             DM(SIGTRA), DM(SGNTRA), 1, NPICET,
D    &             DM(FOTSIG))
D     CALL IMPTDT ('MEILLEURES EVOLUTIONS DES CONTRAINTES '//
D    &             'APRES UN 1er PASSAGE NULLES? ',
D    &              DM(FOTSIG), NPICET, NBFEPS)
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
C 
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule la meilleure evolution des deplacements
C     et la meilleure evolution des contraintes sur les quantites
C     admissibles a zero construites precedemment .
C 
C     On modifie en fin de routine le fichier des DELTA SIGMA BARRE .
C 
C     On envoie comme arguments :
C 
C     E ...... FTDEOP   evolution optimale associee a la derniere
C                       deformation "CAD0" construite
C     E ...... DEPSOL   solution developpee Fourier en deplacement
C                       de ce probleme
C     E ...... EPSSOL   solution en deformation de ce probleme
C     E ...... SAUSOL   solution en saut de ce probleme
C     E ...... SIGZER   contraintes rigoureusement admissibles
C                       a zero elements finis associees en sortie
C     E ...... SGNZER   contraintes normales rigoureusement
C                       admissibles a zero associees en sortie
C     E ...... NPICDE   piquet de temps de debut d'integration
C     E ...... NPICFI   piquet de temps de fin d'integration
C 
      SUBROUTINE ETORTH (FTDEOP, DEPSOL, EPSSOL, SAUSOL,
     &                   SIGZER, SGNZER, NPICDE, NPICFI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NPICDE, NPICFI
C 
      DOUBLE PRECISION  FTDEOP(NPICET)
      DOUBLE PRECISION  DEPSOL(NDDL*NBMAT), EPSSOL(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION  SIGZER(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION  SGNZER(NSAU*NTETA*NGAU2)
      DOUBLE PRECISION  SAUSOL(NSAU*NTETA*NGAU2)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER AM2LC, ADM2LC, FTSCH5, FTECH6
      INTEGER EPSCH8, SGNC12, SIGC11, SAUCH9, FTRE10
      INTEGER DEPCH0, LONDEP
      INTEGER NDEBSI, NBCPOR
      INTEGER NUSIG, DSIGTR, DSGNTR
      INTEGER FOTSIG, LONDER
      INTEGER NNOUV, LONEPS, LONSAU, NNOUVA
      INTEGER SIGTRA, SGNTRA, FTREE7, COESIG
      INTEGER ADINT, ADCOU, ADSOU, ADSNT
      INTEGER TABNIE(2), TABNIS(2), TABNID(2)
C 
C     LONGUEURS DE TABLEAUX DE TRAVAIL
C 
      INTEGER TRAV1, TRAV2, ADINTE
C 
C     ADRESSE DE TRAVAIL
C 
      INTEGER  ADTEP1, ADTEP2
      INTEGER ADREEL
      DOUBLE PRECISION DIVIS, NORMF
C 
C     POUR LES VERIFS
C 
      DOUBLE PRECISION NORVER
      INTEGER  DBCHEP, COEFF, EPSTRA, SAUTRA, DEPTRA, DEPSTR, DSAUTR
      INTEGER  DDEPTR, NUEPS, ADEPS1, ADSAU1, ADDEP1, ADEPS2, ADSAU2
      INTEGER  ADDEP2, FOTDEP
      LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ETORTH')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL MESSAO ('ON ENTRE DANS '// IDPROG)
C 
      LONDER = NDDL*NTETA
      LONDEP = NDDL*NBMAT
      LONEPS = NEPS*NGAU1*NTETA
      LONSAU = NSAU*NGAU2*NTETA
C 
C     Appel des tableaux pour la rigidite et la souplesse couches
C 
      CALL ADTBDM ('HOO-COUCHE', ADCOU)
      CALL ADTBDM ('SOU-COUCHE', ADSOU)
C 
C     Appel des tableaux pour la rigidite et la souplesse interfaces
C 
      IF (NBINT .GT. 0) THEN
        CALL ADTBDM ('HOO-INTERF', ADINT)
        CALL ADTBDM ('SOU-INTERF', ADSNT)
      ELSE
        ADINT = ADCOU
        ADSNT = ADSOU
      END IF
C 
C     Appel du tableau pour les valeurs des intervalles de temps
C 
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
C 
C     Appel des tableaux des deplacements, deformations et contraintes
C     admissibles totales couches
C  
      CALL ADTBDM ('DEPLA-ADMI', DEPCH0)
      CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
      CALL ADTBDM ('SIG-AD-TOT', SIGC11)
C 
C     Appel des tableaux des sauts et contraintes admissibles totales
C     interfaces
C  
      IF (NBINT .GT. 0) THEN
        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
      ELSE
        SAUCH9 = EPSCH8
        SGNC12 = SIGC11
      ENDIF
C 
C     Appel des tableaux pour les fonctions du temps des evolutions
C     en contraintes et deformations admissibles
C 
      CALL ADTBDM ('TEMPS-SIGM', FTSCH5)
      CALL ADTBDM ('TEMPS-EPSI', FTECH6)
C 
C     Appel des tableaux pour les fonctions du temps des evolutions
C     en contraintes et deformations admissibles totales
C 
      CALL ADTBDM ('TEMPS-REEL', FTREE7)
      CALL ADTBDM ('TEMP-SI-RE', FTRE10)
      TABNIE(1) = CHARAX
      TABNIE(2) = LONEPS
      TABNIS(1) = CHARAX
      TABNIS(2) = LONSAU
      TABNID(1) = CHARAX
      TABNID(2) = LONDER
C 
C ***********************************************************************
C 
C     Les champs admissibles sont construits un a un dans l'etape
C     cinematique de l'etape globale. Comme il y a eu une etape
C     preliminaire, ce nouveau champ est theoriquement orthogonal
C     aux champs utilises pour l'optimisation dans etapre et aux
C     champs precedents construits pour satisfaire au critere de
C     precision de l'etape globale; il suffit de le normer et de
C     modifier en consequence l'evolution optimale correspondante.
C 
C ***********************************************************************
C 
      CALL NOKGLO (ADCOU, ADINT, 1, 1, EPSSOL, SAUSOL, NORMF)
      NORMF = DSQRT(NORMF)
      DIVIS = 1.D0/NORMF
C 
C     Normalisation des quantites attachees au deplacement
C     solution par rapport a Ke
C 
      CALL MUMARE (DIVIS, LONDEP, DEPSOL, DEPSOL)
      CALL MUMARE (DIVIS, LONEPS, EPSSOL, EPSSOL)
      CALL MUMARE (DIVIS, LONSAU, SAUSOL, SAUSOL)
      CALL POUSMD (LONDER, ADREEL)
      CALL VRTSM (1, DEPSOL, DM(ADREEL))
C 
C     Modification de la fonction du temps optimale correspondante
C 
      CALL MUMARE (NORMF, NPICET, FTDEOP, FTDEOP)
D     CALL IMPTDT ('EVOLUTION DE LA CORRECTION EN DEFORMATION DANS '
D    &             //IDPROG//' VERSION NORMALE ', FTDEOP, NPICET, 1)
C 
C     Rangement des quantites admissibles a zero
C 
      CALL RCHARE (0, LONDER, DM(ADREEL), DM(DEPCH0))
      CALL RCHARE (8, LONEPS, EPSSOL, DM(EPSCH8))
      IF (NBINT .GT. 0) THEN
         CALL RCHARE (9, LONSAU, SAUSOL, DM(SAUCH9))
      END IF
C 
C     Fonction du temps deformations admissibles a zero
C 
      CALL RCHAMP (6, NPICET, FTDEOP, DM(FTECH6))
C 
C     Fonction du temps deformations totales
C 
      CALL RCHARE (7, NPICET, FTDEOP, DM(FTREE7))
C 
      GOTO 1000
C 
C     Fin du travail normal sur les deformations :
C     Le champ cree est orthogonal par construction (le residu est lui-meme
C     orthogonal). Ce qui suit est ecrit a titre de verificationseulement.
C 
C ***********************************************************************
C 
C     Sequence de verification de l'orthogonalite d'orthonormalisation
C 
C ***********************************************************************
C 
C     Un seul champ admissible a ete construit a l'etape de correction de l'etape
C     globale, on verifie son independance et on le retire si neccessaire
C 
      NNOUV  = 1
C 
C     MODIFICATION DES CHAMPS POUR LES DEFORMATIONS
C 
C     DBCHEP est le numero du premier champ servant au calcul de
C     l'accroissement admissible a zero en deformation calcule a l'etape
C     globale precedente
C 
C     DBCHEP = DEADTR-NBFEPS+1
C 
C     On orthogonalise le dernier champ cree par rapport a tous les
C     champs servant a l'orthogonalisation qui sont au nombre de
C     NBFEPS pour les deformations.
C 
C     On extrait du tableau des deformations totales stockees
C     (charax, loneps) les deformations allant du numero DBCHEP a DEADTR
C 
      DBCHEP = DEADTR-NBFEPS+1
C 
C     NBCPOR est le nombre de champ servant a l'orthogonalisation
C 
      NBCPOR =  NBFEPS
      TRAV1  =  NBCPOR*(LONEPS+LONSAU+LONDER)+ NBCPOR
      CALL POUSMD (TRAV1, COEFF)
      EPSTRA  = COEFF  + NBCPOR
      SAUTRA  = EPSTRA + NBCPOR*LONEPS
      DEPTRA  = SAUTRA + NBCPOR*LONSAU
      DEPSTR = EPSTRA
      DSAUTR = SAUTRA
      DDEPTR = DEPTRA
C 
      DO NUEPS = DBCHEP, DEADTR
C 
        CALL EXTRAD (DM(EPSCH8), 2, TABNIE(1),
     &               2, NUEPS, DM(DEPSTR), LONEPS)
        DEPSTR = DEPSTR+LONEPS
        CALL EXTRAD (DM(DEPCH0), 2, TABNID(1),
     &               2, NUEPS, DM(DDEPTR), LONDER)
        DDEPTR = DDEPTR+LONDER
        IF (NBINT .GT. 0) THEN
          CALL EXTRAD (DM(SAUCH9), 2, TABNIS(1),
     &                 2, NUEPS, DM(DSAUTR), LONSAU)
          DSAUTR = DSAUTR+LONSAU
        END IF
C 
      END DO
C 
C     Orthogonalisation du nouveau deplacement developpe admissible a zero, au sens de K0
C 
      ADEPS1 = EPSTRA
      ADSAU1 = SAUTRA
      ADDEP1 = DEPTRA
      ADEPS2 = EPSTRA+(NBCPOR-1)*LONEPS
      ADSAU2 = SAUTRA+(NBCPOR-1)*LONSAU
      ADDEP2 = DEPTRA+(NBCPOR-1)*LONDER
      CALL OREDK0 (ADCOU, ADINT, NBCPOR-1, NNOUV,
     &             LONEPS, LONSAU, LONDER,
     &             DM(ADEPS1), DM(ADSAU1), DM(ADDEP1),
     &             DM(ADEPS2), DM(ADSAU2), DM(ADDEP2),
     &             DM(ADEPS2), DM(ADSAU2), DM(ADDEP2),
     &             NNOUVA, DM(COEFF))
      NNOUV = NNOUVA
      NBDEPR = NBDEPR-1
      DEADTR = DEADTR-1
      SAADTR = SAADTR-1
      CALL IMPET ('NOMBRE DE NOUVEAUX VECTEURS ', NNOUVA)
      IF (NNOUVA .EQ. 0) THEN
         CALL IMPET ('PAS DE NOUVELLE EVOLUTION DES DEPLACEMENTS '
     &               //IDPROG, NNOUVA)
C 
C       On decremente le nombre de champ servant au calcul des delta admissibles
C 
       NBCPOR  = NBCPOR -1
      ELSE
C 
C       Si le nouveau vecteur n'etait pas combinaison lineaire des precedents
C       on le stocke a la place de son ancienne valeur, sinon il est inutile.
C       Decrementation du nombre de vecteurs stockes, car on ne passe pas dans
C       le IF (qui incremente par appel a RCHARE) dans le else on decremente
C       le nonbre d'evolutions reelles qui autrement ont deja ete stockees.
C 
        CALL RCHARE (0, LONDER, DM(ADDEP2), DM(DEPCH0))
        CALL RCHARE (8, LONEPS, DM(ADEPS2), DM(EPSCH8))
        IF (NBINT .GT. 0) THEN
          CALL RCHARE (9, LONSAU, DM(ADSAU2), DM(SAUCH9))
        END IF
C 
C       Verification de la norme
C 
D       CALL NOKGLO (ADCOU, ADINT, 1, 1, DM(ADEPS2), DM(ADSAU2), NORVER)
D       CALL IMPDT ('NORME = 1? DANS '//IDPROG, NORVER)
        CALL POUSMD (NPICET, FOTDEP)
C 
C       A priori on a calcule l'evolution optimale sur la derniere
C       fonction de l'espace uniquement. Calcul de
C       intesp(delta(-sigcha* ADEPS2) pour la valeur de delta(sigcha)
C       correspondant au residu en contrainte.
C 
C       Cette evolution est aussi, puisque les deformations sont
C       orthonormalisees au sens de K0, la meilleure evolution des
C       deplacements fi(t) pour le delta admissible de l'etape globale actuelle.
C 
        CALL SCHAST (0, 'sigmch', 'sinoch', NNOUV,
     &               DM(ADEPS2), DM(ADSAU2), NPICDE, NPICFI, DM(FOTDEP))
        CALL IMPTDT ('MEILLEURE EVOLUTION DE LA CORRECTION EN '//
     &              'DEFORMATION DANS ' //IDPROG, DM(FOTDEP), NPICET, 1)
C 
C ***********************************************************************
C 
C       Pour aller ranger les evolutions au bon endroit on remet
C       NBFEPS a zero puis on range les nouvelles evolutions des
C       deplacements et deformations.
C 
C ***********************************************************************
C 
C       Fonction du temps pour les deformations admissibles a zero
C 
        NBFEPS = NBFEPS-1
        NBDPTR = NBDPTR-1
        CALL RCHAMP (6, NPICET, DM(FOTDEP), DM(FTECH6))
C 
C       Fonction du temps pour les deformations totales
C 
        CALL RCHARE (7, NPICET, DM(FOTDEP), DM(FTREE7))
C 
C       Appel a la routine remplisssant le fichier scb-numero d'etape locale
C       en y rangeant - delta(sigcba)-K(delta-tild(epsn))
C 
        CALL RSICBA (0, ADCOU, ADINT, ADSOU, ADSNT, 'sigmch', 'sinoch',
     &                NBFEPS, NBFEPS, DEADTR, DEADTR,
     &                EPSCH8, SAUCH9, FTECH6, NPICDE, NPICFI)
C 
C      TEST normalement la nouvelle evolution optimale
C      devrait etre nulle puisque l'on a deja retire
C      l'evolution optimale correspondant a adeps2
C 
D       CALL SCHAST (0, 'sigmch', 'sinoch', NNOUV, DM(ADEPS2),
D    &               DM(ADSAU2), NPICDE, NPICFI, DM(FOTDEP))
D       CALL IMPTDT ('MEILLEURE EVOLUTION DE LA CORRECTION EN '//
D    &         'DEFORMATION APRES UNE PREMIERE OPTIMISATION => 0? DANS '
D    &          //IDPROG, DM(FOTDEP), NPICET, 1)
C 
      END IF
C 
C     Si on est en version Debug le travail a deja ete fait
C     Appel a la routine remplissant le fichier scb-numero d'etape locale
C     en y rangeant -delta(sigcba)- K(delta-tild(epsn))
C 
1000  CONTINUE
      CALL RSICBA (0, ADCOU, ADINT, ADSOU, ADSNT,
     &            'sigmch', 'sinoch',
     &             NBFEPS, NBFEPS, DEADTR, DEADTR,
     &             EPSCH8, SAUCH9, FTECH6, NPICDE, NPICFI)
C 
C     NDEBSI est le numero du premier champ servant au calcul de
C     l'accroissement admissible a zero en contrainte calcule a l'etape
C     globale precedente.
C 
C     Pour comprendre => [~/alliant/delami/renseignement/nom_champ]
C 
C     NDEBSI = COTORE-NBFSIG
C 
C     On orthogonalise le dernier champ cree par rapport a tous les
C     champs servant a l'orthogonalisation qui sont au nombre de NBFSIG.
C 
      NDEBSI = COTORE-NBFSIG+1
      NBCPOR =  NBFSIG+1
      TRAV2  =  NBFSIG*(LONEPS+LONSAU)+ NBCPOR
      CALL POUSMD (TRAV2, COESIG)
      SIGTRA = COESIG+NBCPOR
      DSIGTR = SIGTRA
      SGNTRA = SIGTRA+NBFSIG*LONEPS
      DSGNTR = SGNTRA
      DO NUSIG = NDEBSI, COTORE
C 
        CALL EXTRAD (DM(SIGC11), 2, TABNIE(1),
     &               2, NUSIG, DM(DSIGTR), LONEPS)
        DSIGTR = DSIGTR+LONEPS
        IF (NBINT .GT. 0) THEN
          CALL EXTRAD (DM(SGNC12), 2, TABNIS(1),
     &                 2, NUSIG, DM(DSGNTR), LONSAU)
          DSGNTR = DSGNTR+LONSAU
        END IF
C 
      END DO
C 
      NNOUV = 1
      ADTEP1 = SIGTRA
      ADTEP2 = SGNTRA
C 
      CALL ORNEK0 (ADSOU, ADSNT, NBCPOR-1, NNOUV, LONEPS, LONSAU,
     &             DM(ADTEP1), DM(ADTEP2), SIGZER, SGNZER,
     &             SIGZER, SGNZER, NNOUVA, DM(COESIG))
C 
      NNOUV = NNOUVA
      IF (NNOUVA .EQ. 0) THEN
        CALL IMPET 
     &      ('PAS DE NOUVELLE EVOLUTION EN CONTRAINTES '
     &       //IDPROG, NNOUVA)
      ELSE
        CALL RCHARE (11, LONEPS, SIGZER, DM(SIGC11))
        IF (NBINT .GT. 0) THEN
          CALL RCHARE (12, LONSAU, SGNZER, DM(SGNC12))
        END IF
C 
C       Calcul de la meilleure evolution des contraintes fi(t)
C       pour le delta admissible de l'etape globale actuelle
C 
        TRAV1  = NBCPOR*NPICET
        CALL POUSMD (TRAV1, FOTSIG)
C 
C       Calcul de INTESP(DELTA(EPSADMI)* EPSCHA) pour la valeur
C       de delta(epscha) correspondant a l'etape locale
C 
        CALL SCHAST (2, 'sigmch', 'sinoch', NNOUV,
     &               SIGZER, SGNZER, NPICDE, NPICFI, DM(FOTSIG))
        CALL MUMARE (-1.D0, NPICET, DM(FOTSIG), DM(FOTSIG))
D       CALL IMPTDT ('MEILLEURE EVOLUTION DE LA CORRECTION EN '//
D    &              'CONTRAINTES DANS '// IDPROG, DM(FOTSIG), NPICET, 1)
C 
C       On range la nouvelle evolution des contraintes.
C       contraintes admissibles a zero
C 
        CALL RCHAMP (5, NPICET, DM(FOTSIG), DM(FTSCH5))
C 
C       contraintes totales
C 
        CALL RCHARE (10, NPICET, DM(FOTSIG), DM(FTRE10))
C 
C       Appel a la routine remplissant le fichier ecbarr
C       en y rangeant (K-1)(delta(sigcha) -  (delta-tild(sign))
C 
        CALL RSICBA (2, ADSOU, ADSNT, ADCOU, ADINT,
     &               'sigmch', 'sinoch',
     &               NBFSIG, NBFSIG, COTORE, COTORE,
     &               SIGC11, SGNC12, FTSCH5, NPICDE, NPICFI)
C 
C       TEST normalement la nouvelle evolution optimale
C       devrait etre nulle puisqu'on a deja retire
C       l'evolution optimale correspondant a ADTEP3
C 
D       CALL SCHAST (2, 'sigmch', 'sinoch', NNOUV,
D    &               SIGZER, SGNZER, NPICDE, NPICFI, DM(FOTSIG))
D       CALL IMPTDT ('MEILLEURE EVOLUTION DE LA CORRECTION EN '//
D    &               'CONTRAINTE APRES UNE PREMIERE OPTIMISATION '//
D    &               '(=> 0?) DANS '// IDPROG, DM(FOTSIG), NPICET, 1)
C 
      END IF
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD     CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule la meilleure evolution des deplacements
C     sur les quantites admissibles a zero construites precedemment.
C 
C     On modifie en fin de routine le fichier des DELTA SIGMA BARRE.
C 
C     On envoie comme arguments :
C 
C     E ...... DEPSOL  La solution developpee fourier en deplacement de ce probleme
C     E ...... EPSSOL  La solution en deformation de ce probleme
C     E ...... SAUSOL  La solution en saut de ce probleme
C     E ...... NPICDE  piquet de temps de debut d'integration
C     E ...... NPICFI  piquet de temps de fin d'integration
C 
      SUBROUTINE ETORDE (DEPSOL, EPSSOL, SAUSOL, NPICDE, NPICFI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'strategie_calcul.h'
C 
      DOUBLE PRECISION  DEPSOL(NDDL*NBMAT), EPSSOL(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION  SAUSOL(NSAU*NTETA*NGAU2)
      INTEGER NPICDE, NPICFI
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  AM2LC, ADM2LC, FTDEOP
      INTEGER  FTECH6, EPSCH8, SAUCH9
      INTEGER  DEPCH0, LONDEP, EPSTRA, SAUTRA, DEPTRA
      INTEGER  NUEPS, NBCPOR
      INTEGER  DEPSTR, DSAUTR, DDEPTR
      INTEGER  LONDER, NNOUV, LONEPS, LONSAU, NNOUVA
      INTEGER  FTREE7, COEFF
      INTEGER  ADINT, ADCOU, ADSOU, ADSNT
      INTEGER  TABNIE(2), TABNIS(2), TABNID(2), FOTDEP
C 
C     LONGUEURS DE TABLEAUX DE TRAVAIL
C 
      INTEGER  TRAV1, ADINTE
C 
C     ADRESSES DE TRAVAIL
C 
      INTEGER  ADREEL, DBCHEP
C 
      DOUBLE PRECISION DIVIS
C 
C     POUR LES VERIFS
C 
CD    DOUBLE PRECISION NORVER
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      CHARACTER*3 BARATE
      PARAMETER (IDPROG='ETORDE')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL MESSAO ('ON ENTRE DANS '//IDPROG)
C 
      LONDER = NDDL*NTETA
      LONDEP = NDDL*NBMAT
      LONEPS = NEPS*NGAU1*NTETA
      LONSAU = NSAU*NGAU2*NTETA
C 
C     Appel des tableaux pour les valeurs des rigidites et
C     souplesses couches
C 
      CALL ADTBDM ('HOO-COUCHE', ADCOU)
      CALL ADTBDM ('SOU-COUCHE', ADSOU)
C 
C     Appel des tableaux pour les valeurs des rigidites et
C     souplesses interfaces
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL ADTBDM ('HOO-INTERF', ADINT)
        CALL ADTBDM ('SOU-INTERF', ADSNT)
C 
      ELSE
C 
        ADINT = ADCOU
        ADSNT = ADSOU
C 
      END IF
C 
C     Appel du tableau des valeurs des intervalles de temps
C 
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
C 
C     Appel des tableaux des accroissements en deplacements et
C     deformations admissibles couches
C  
      CALL ADTBDM ('DEPLA-ADMI', DEPCH0)
      CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
C 
      IF (NBINT .GT. 0) THEN
C 
C       Appel des tableaux des accroissements en sauts de
C       deplacement interfaces
C  
        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
C 
      ELSE
C 
        SAUCH9 = EPSCH8
C 
      ENDIF
C 
C     Reservation de place pour une fonction du temps
C 
      TRAV1  = NPICET
      CALL POUSMD (TRAV1, FTDEOP)
C 
C     Calcul de la norme au sens de K0 du champ EPSSOL solution du GCADMI
C 
      CALL NOKGLO (ADCOU, ADINT, 1, 1, EPSSOL, SAUSOL, DIVIS)
C 
      IF (DIVIS .LT. 1.D -10) THEN
C 
         CALL IMPDT ('NORME DE LA NOUVELLE DEFORMATION ADMISSIBLE A '//
     &               'ZERO TRES PETITE ', DIVIS)
C 
         GOTO 2000
C 
      END IF
C 
C     Normalisation de EPSSOL et SAUSOL
C 
      DIVIS = 1.D0/ DSQRT(DIVIS)
      CALL MUMARE (DIVIS, LONEPS, EPSSOL, EPSSOL)
      IF (NBINT .GT. 0) THEN
        CALL MUMARE (DIVIS, LONSAU, SAUSOL, SAUSOL)
      END IF
C 
C     Calcul de INTESP(DELTA(EPSADMI)* EPSCHA) pour la valeur
C     de delta(epscha) correspondant a l'etape locale
C 
C    CALL IDENTI (NBETLC, BARATE)
      CALL SCHAST (0, 'sigmch', 'sinoch', 1, EPSSOL,
     &             SAUSOL, NPICDE, NPICFI, DM(FTDEOP))
C 
      CALL IMPTDT ('MEILLEURE EVOLUTION DE LA CORRECTION EN '//
     &             'DEFORMATION DANS ' //IDPROG, DM(FTDEOP), NPICET, 1)
C 
C     Modif Steph. La regularisation semble mettre un sacre bordel.
C     Pour voir, je vire dans le cas VALPRO = .TRUE.
C 
CD   IF (.NOT. VALPRO) THEN
CD     CALL WVISC1 (NPICET, .9D0, DM(ADINTE), DM(FTDEOP), DM(FTDEOP))
CD    END IF
C 
C     Fonctions du temps pour les champs des accroissements de deformations
C     admissibles a zero
C 
      CALL ADTBDM ('TEMPS-EPSI', FTECH6)
C 
C     Fonctions du temps pour les champs de deformations admissibles
C 
      CALL ADTBDM ('TEMPS-REEL', FTREE7)
C 
      TABNIE(1) = CHARAX
      TABNIE(2) = LONEPS
C 
      TABNIS(1) = CHARAX
      TABNIS(2) = LONSAU
C 
      TABNID(1) = CHARAX
      TABNID(2) = LONDER
C 
C     Rangement des fonctions du temps pour les champs des accroissements de 
C     deformations admissibles a zero
C     NBFEPS = MIN0(NBFEPS+1, CHAMAX)
C 
      CALL RCHAMP (6, NPICET, DM(FTDEOP), DM(FTECH6))
C 
C     Rangement des fonctions du temps pour les champs de deplacements admissibles
C     NBDPTR = NBDPTR+1
C 
      CALL RCHARE (7, NPICET, DM(FTDEOP), DM(FTREE7))
C 
      CALL POUSMD (LONDER, ADREEL)
      CALL VRTSM  (1, DEPSOL, DM(ADREEL))
C 
C     Rangement du nouveau champ des deplacements admissibles
C     NBDEPR = NBDEPR+1
C 
      CALL RCHARE (0, LONDER, DM(ADREEL), DM(DEPCH0))
C 
C     Rangement du nouveau champ des deformations admissibles couches
C     DEADTR = DEADTR+1
C 
      CALL RCHARE (8, LONEPS, EPSSOL, DM(EPSCH8))
C 
      IF (NBINT .GT. 0) THEN
C 
C       Rangement du nouveau champ des sauts admissibles interfaces
C       SAADTR = SAADTR+1
C 
        CALL RCHARE (9, LONSAU, SAUSOL, DM(SAUCH9))
C 
      END IF
C 
C     Actualisation du residu correspondant au GCADMI dans sigmch (sinoch)
C 
      CALL RSICBA (0, ADCOU, ADINT, ADSOU, ADSNT, 'sigmch',
     &             'sinoch', NBFEPS, NBFEPS, DEADTR, DEADTR,
     &             EPSCH8, SAUCH9, FTECH6, NPICDE, NPICFI)
C 
      IF ((.NOT. LETAPP) .OR. BIGNET) GOTO 2000
C 
C     Stockage du champ des deformations admissibles totales
C     avant orthogonalisation du nouveau champ (pour verification)
C 
      CALL DEORMO (0)
C 
C     Stockage du champ des accroissements des deformations admissibles a zero
C     avant orthogonalisation du nouveau champ (pour verification)
C 
      CALL DEORMO (1)
C 
C     Un seul champ (??) admissible a ete construit a l'etape de
C     correction  de l'etape globale, on verifie son independance
C     et on le retire si neccessaire.
C 
      NNOUV  = 1
C 
C     ORTHOGONALISATION DES CHAMPS DES DEFORMATIONS
C 
C     DBCHEP est le numero du premier champ des accroissements des deformations
C     admissibles a zero en deformation calcule a l'etape globale precedente.
C 
      DBCHEP = DEADTR-NBFEPS+1
C 
C     On orthogonalise le dernier champ cree par rapport a tous les champs servant
C     a l'orthogonalisation qui sont au nombre de NBFEPS pour les deformations.
C 
C     On extrait du tableau des deformations totales stockees (charax, loneps) les
C     deformations allant du numero DBCHEP a DEADTR.
C 
C     NBCPOR est le nombre de champs servant a l'orthogonalisation.
C 
      NBCPOR =  NBFEPS
      TRAV1  =  NBCPOR*(LONEPS+LONSAU+LONDER)+ NBCPOR
C 
      CALL POUSMD (TRAV1, COEFF)
C 
      EPSTRA = COEFF  + NBCPOR
      SAUTRA = EPSTRA + NBFEPS*LONEPS
      DEPTRA = SAUTRA + NBFEPS*LONSAU
      DEPSTR = EPSTRA
      DSAUTR = SAUTRA
      DDEPTR = DEPTRA
C 
      DO NUEPS = DBCHEP, DEADTR
C 
        CALL EXTRAD (DM(EPSCH8), 2, TABNIE(1),
     &               2, NUEPS, DM(DEPSTR), LONEPS)
C 
        DEPSTR = DEPSTR+LONEPS
C 
        CALL EXTRAD (DM(DEPCH0), 2, TABNID(1),
     &               2, NUEPS, DM(DDEPTR), LONDER)
C 
        DDEPTR = DDEPTR+LONDER
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL EXTRAD (DM(SAUCH9), 2, TABNIS(1),
     &                 2, NUEPS, DM(DSAUTR), LONSAU)
C 
          DSAUTR = DSAUTR+LONSAU
C 
        END IF
C 
      END DO
C 
C     Orthogonalisation des champs des accroissements en deplacement, deformation
C     couches et interface admissible a zero, au sens de K0
C 
      CALL OREDK0 (ADCOU, ADINT, NBCPOR-1, 1,
     &             LONEPS, LONSAU, LONDER,
     &             DM(EPSTRA), DM(SAUTRA), DM(DEPTRA),
     &             EPSSOL, SAUSOL, DEPSOL,
     &             EPSSOL, SAUSOL, DEPSOL,
     &             NNOUVA, DM(COEFF))
C 
      NNOUV = NNOUVA
C 
C     Modification des fonctions du temps reelles pour les champs de deformations
C     admissibles pour tenir compte des orthogonalisations
C 
      CALL MOCHAR (NPICET, NBDPTR, DM(COEFF), NBCPOR, DM(FTREE7))
C 
C     Modification des fonctions du temps pour les champs des accroissements de
C     deformations pour tenir compte des orthogonalisations.
C 
      CALL MOCHAM (NPICET, NBFEPS, DM(COEFF), NBCPOR, DM(FTECH6))
C  
      IF (NNOUVA .EQ. 0) THEN
C 
        NBFEPS = NBFEPS -1
        NBDPTR = NBDPTR -1
        DEADTR = DEADTR -1
        SAADTR = SAADTR -1
        NBDEPR = NBDEPR -1
C       NBCPOR = NBCPOR -1
C 
        CALL IMPET ('PAS DE NOUVELLE EVOLUTION DES DEFORMATIONS DANS '
     &              //IDPROG, NNOUVA)
C 
      ELSE
C 
        DEADTR = DEADTR -1
        SAADTR = SAADTR -1
        NBDEPR = NBDEPR -1
C 
        CALL VRTSM (1, DEPSOL, DM(ADREEL))
C 
        CALL RCHARE (0, LONDER, DM(ADREEL), DM(DEPCH0))
        CALL RCHARE (8, LONEPS, EPSSOL, DM(EPSCH8))
C 
        IF (NBINT .GT. 0) THEN
C 
           CALL RCHARE (9, LONSAU, SAUSOL, DM(SAUCH9))
C 
        END IF
C 
      END IF
C 
C 
C     Verification du champ des deformations admissibles totales
C     apres orthogonalisation du nouveau champ
C 
      CALL FIORMO(0)
C 
C     Verification du champ des accroissements de deformations admissibles a zero
C     apres orthogonalisation du nouveau champ
C 
      CALL FIORMO(1)
C 
      IF (NNOUV .NE. 0) THEN
C 
C       TEST normalement la nouvelle evolution optimale devrait etre nulle
C       puisqu'on a deja retire l'evolution optimale
C 
        CALL POUSMD (NPICET, FOTDEP)
C 
        CALL SCHAST (0, 'sigmch', 'sinoch', NNOUV,
     &               EPSSOL, SAUSOL, NPICDE, NPICFI, DM(FOTDEP))
C 
        CALL IMPTDT ('MEILLEURE EVOLUTION DE LA CORRECTION EN '//
     &               'DEFORMATION APRES UNE PREMIERE OPTIMISATION '//
     &               'NULLE DANS '//IDPROG//' MAIS REGULARISATION ',
     &                DM(FOTDEP), NPICET, 1)
C 
      END IF
C 
2000  CONTINUE
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule la meilleure evolution  des contraintes
C     sur les quantites admissibles a zero construites precedemment.
C 
C     On modifie en fin de routine le fichier des DELTA EPSI BARRE.
C 
C     Etape d'orthogonalisation des contraintes et calcul de la fonction
C     du temps associee.
C 
C     On envoie comme arguments :
C 
C     E ...... SIGZER les contraintes rigoureusement admissibles a zero elements finis
C     E ...... SGNZER les contraintes normales rigoureusement admissibles a zero
C     E ...... NPICDE  piquet de temps de debut d'integration
C     E ...... NPICFI  piquet de temps de fin d'integration
C 
      SUBROUTINE ETORCO (SIGZER, SGNZER, NPICDE, NPICFI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'strategie_calcul.h'
C 
      DOUBLE PRECISION  SIGZER(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION  SGNZER(NSAU*NTETA*NGAU2)
C 
      INTEGER NPICDE, NPICFI
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  AM2LC, ADM2LC, FTSCH5
      INTEGER  SGNC12, SIGC11, FTRE10
      INTEGER  NDEBSI, NBCPOR, NUSIG
      INTEGER  DSIGTR, DSGNTR, FOTSIG
      INTEGER  NNOUV, LONEPS, LONSAU, NNOUVA
      INTEGER  SIGTRA, SGNTRA, COESIG
      INTEGER  ADINT, ADCOU, ADSOU, ADSNT
      INTEGER  TABNIE(2), TABNIS(2)
C 
C     LONGUEURS DE TABLEAUX DE TRAVAIL
C 
      INTEGER  TRAV1, TRAV2, ADINTE
C 
C     ADRESSE DE TRAVAIL
C 
      DOUBLE PRECISION DIVIS
C 
C     POUR LES VERIFS
C 
      DOUBLE PRECISION NORVER
C 
      LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      CHARACTER*3 BARATE
      PARAMETER (IDPROG='ETORC0')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL MESSAO ('ON ENTRE DANS '// IDPROG)
C 
      LONEPS = NEPS*NGAU1*NTETA
      LONSAU = NSAU*NGAU2*NTETA
C 
      CALL ADTBDM ('HOO-COUCHE', ADCOU)
      CALL ADTBDM ('SOU-COUCHE', ADSOU)
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL ADTBDM ('HOO-INTERF', ADINT)
        CALL ADTBDM ('SOU-INTERF', ADSNT)
C 
      ELSE
C 
         ADINT = ADCOU
         ADSNT = ADSOU
C 
      END IF
C 
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
      CALL ADTBDM ('SIG-AD-TOT', SIGC11)
C 
C     POUR LES INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
C 
      ELSE
C 
        SGNC12 = SIGC11
C 
      ENDIF
C 
C     Evolutions admissibles a zero
C 
      CALL ADTBDM ('TEMPS-SIGM', FTSCH5)
C 
C     Evolutions admissibles
C 
      CALL ADTBDM ('TEMP-SI-RE', FTRE10)
C 
      TABNIE(1) = CHARAX
      TABNIE(2) = LONEPS
C 
      TABNIS(1) = CHARAX
      TABNIS(2) = LONSAU
C 
C ***********************************************************************
C 
C     TRAITEMENT  DES CONTRAINTES
C 
C ***********************************************************************
C 
C     Calcul de la meilleure evolution des contraintes fi(t)
C     pour le delta admissible de l'etape globale actuelle
C 
      TRAV1 = NPICET
C 
      CALL POUSMD (TRAV1, FOTSIG)
C 
      CALL NOKGLO (ADSOU, ADSNT, 1, 1, SIGZER, SGNZER, DIVIS)
C 
      IF (DIVIS .LT. 1.D-10) THEN
C 
         CALL IMPDT (
     &        'NORME DE LA NOUVELLE CONTRAINTE ADMISSIBLE A ZERO '//
     &        'TRES PETITE ', DIVIS)
C 
         GOTO 2000
C 
      END IF
C 
      DIVIS = 1.D0/ DSQRT(DIVIS)
C 
      CALL MUMARE (DIVIS, LONEPS, SIGZER, SIGZER)
C 
      IF (NBINT .GT. 0) THEN
        CALL MUMARE (DIVIS, LONSAU, SGNZER, SGNZER)
      END IF
C 
C     Calcul de INTESP(DELTA(EPSADMI)* EPSCHA) pour la valeur
C     de delta(epscha) correspondant a l'etape locale
C 
C     CALL IDENTI (NBETLC, BARATE)
      CALL SCHAST (2, 'sigmch', 'sinoch', 1,
     &             SIGZER, SGNZER, NPICDE, NPICFI, DM(FOTSIG))
C 
      CALL MUMARE (-1.D0, NPICET, DM(FOTSIG), DM(FOTSIG))
C 
C     Modif Steph. La regularisation semble mettre un sacre bordel.
C     Pour voir, je vire dans le cas VALPRO = .TRUE.
C 
CD    IF (.NOT. VALPRO) THEN
CD      CALL WVISC1 (NPICET, .9D0, DM(ADINTE), DM(FOTSIG), DM(FOTSIG))
CD    END IF
C 
      CALL IMPTDT ('MEILLEURE EVOLUTION DE LA CORRECTION EN '//
     &             'CONTRAINTE DANS ' //IDPROG, DM(FOTSIG), NPICET, 1)
C 
C     Rangement des fonctions du temps pour les champs des accroissements de 
C     contraintes admissibles a zero
C     NBFSIG = MIN0(NBFSIG+1, CHAMAX)
C 
      CALL RCHAMP (5, NPICET, DM(FOTSIG), DM(FTSCH5))
C 
C     Rangement des fonctions du temps pour les champs de contraintes admissibles
C     EVCOTR = EVCOTR+1
C 
      CALL RCHARE (10, NPICET, DM(FOTSIG), DM(FTRE10))
C 
C     Rangement du nouveau champ des contraintes admissibles couches
C     COTORE = COTORE+1
C 
      CALL RCHARE (11, LONEPS, SIGZER, DM(SIGC11))
C 
      IF (NBINT .GT. 0) THEN
C 
C       Rangement du nouveau champ des contraintes admissibles interfaces
C       CNTORE = CNTORE+1
C 
        CALL RCHARE (12, LONSAU, SGNZER, DM(SGNC12))
C 
      END IF
C 
C     Actualisation du residu correspondant au GCADMI dans sigmch (sinoch)
C 
      CALL RSICBA (2, ADSOU, ADSNT, ADCOU, ADINT,
     &             'sigmch', 'sinoch',
     &             NBFSIG, NBFSIG, COTORE, COTORE,
     &             SIGC11, SGNC12, FTSCH5, NPICDE, NPICFI)
C 
C     Si il n'y a pas d'etape preliminaire pas besoin d'orthogonaliser.
C 
      IF ((.NOT. LETAPP) .OR. BIGNET) GOTO 2000
C 
C     Stockage du champ des contraintes admissibles totales
C     avant orthogonalisation du nouveau champ (pour verification)
C 
      CALL DEORMO (2)
C 
C     Stockage du champ des accroissements des contraintes admissibles a zero
C     avant orthogonalisation du nouveau champ (pour verification)
C 
      CALL DEORMO (3)
C 
C     NDEBSI est le numero du premier champ servant au calcul de l'accroissement
C     admissible a zero en contrainte calcule a l'etape globale precedente
C 
C     On orthogonalise le dernier champ cree par rapport a tous les
C     champs servant a l'orthogonalisation qui sont au nombre de NBFSIG
C 
      NDEBSI = COTORE-NBFSIG+1
C 
      NBCPOR =  NBFSIG
      TRAV2  =  NBFSIG*(LONEPS+LONSAU)+ NBCPOR
C 
      CALL POUSMD (TRAV2, COESIG)
C 
      SIGTRA = COESIG+NBCPOR
      DSIGTR = SIGTRA
      SGNTRA = SIGTRA+NBFSIG*LONEPS
      DSGNTR = SGNTRA
C 
      DO NUSIG = NDEBSI, COTORE
C 
        CALL EXTRAD (DM(SIGC11), 2, TABNIE(1),
     &               2, NUSIG, DM(DSIGTR), LONEPS)
C 
        DSIGTR = DSIGTR+LONEPS
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL EXTRAD (DM(SGNC12), 2, TABNIS(1),
     &                 2, NUSIG, DM(DSGNTR), LONSAU)
C 
          DSGNTR = DSGNTR+LONSAU
C 
        END IF
C 
      END DO
C 
      NNOUV = 1
C 
C     Orthogonalisation des nouveaux vecteurs
C 
      CALL ORNEK0 (ADSOU, ADSNT, NBCPOR-1, NNOUV, LONEPS, LONSAU,
     &             DM(SIGTRA), DM(SGNTRA), SIGZER, SGNZER,
     &             SIGZER, SGNZER, NNOUVA, DM(COESIG))
C 
      NNOUV = NNOUVA
C 
C     Modification des fonctions du temps reelles pour les champs de contraintes
C     admissibles pour tenir compte des orthogonalisations
C 
      CALL MOCHAR (NPICET, EVCOTR, DM(COESIG), NBCPOR, DM(FTRE10))
C 
C     Modification des fonctions du temps pour les champs des accroissements de
C     contraintes pour tenir compte des orthogonalisations.
C 
      CALL MOCHAM (NPICET, NBFSIG, DM(COESIG), NBCPOR, DM(FTSCH5))
C 
      IF (NNOUVA .EQ. 0) THEN
C 
        NBFSIG = NBFSIG -1
        EVCOTR = EVCOTR -1
        NBCPOR = NBCPOR -1
        COTORE = COTORE -1
        CNTORE = CNTORE -1
C 
        CALL IMPET ('PAS DE NOUVELLE EVOLUTION DES CONTRAINTES DANS '
     &               //IDPROG, NNOUVA)
C 
      ELSE
C 
      COTORE = COTORE-1
      CNTORE = CNTORE-1
C 
      CALL RCHARE (11, LONEPS, SIGZER, DM(SIGC11))
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL RCHARE (12, LONSAU, SGNZER, DM(SGNC12))
C 
        END IF
C 
      END IF
C 
C     Verification du champ des contraintes admissibles totales
C     apres orthogonalisation du nouveau champ
C 
CD    CALL FIORMO (0)
D     CALL FIORMO (2)
C 
C     Verification du champ des accroissements  de contraintes admissibles a zero
C     apres orthogonalisation du nouveau champ
C 
D     CALL FIORMO (3)
C 
      IF (NNOUVA .NE. 0) THEN
C 
C       TEST normalement la nouvelle evolution optimale
C       devrait etre nulle puisqu'on a deja retiree.
C 
        CALL SCHAST (2, 'sigmch', 'sinoch', NNOUV,
     &               SIGZER, SGNZER, NPICDE, NPICFI, DM(FOTSIG))
C 
        CALL IMPTDT ('MEILLEURE EVOLUTION DE LA CORRECTION EN '//
     &               'CONTRAINTE APRES UNE PREMIERE OPTIMISATION '//
     &               '(=> 0?) DANS '//IDPROG//' MAIS REGULARISATION ',
     &                DM(FOTSIG), NPICET, 1)
C 
      END IF
C 
2000  CONTINUE
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine remplit le fichier direct sigmch
C     en y rangeant delta(sigcha)-K(delta-tild(epsn))
C 
C     On envoie comme arguments :
C 
C     TYPE  = 0 => On envoie a ranger des deformations
C                  entre npicd et npicfi sinon on ne modifie rien
C     TYPE  = 1 => On envoie a ranger des contraintes
C 
C     Le probleme a resoudre est sign+(-sigc)=K * espsn
C 
C     si type = 0 on envoie les deformations et il s'agit
C     de retrancher a (-sigc) K * epsn
C     si type = 1 on envoie les contraintes et il s'agit
C     d'ajouter a sign(-sigc)
C 
C     ADCOU EST L'ADRESSE DE DEPART DU TABLEAU :
C 
C     HOO-COUCHE => NORME AU SENS DE K0
C     SOU-COUCHE => NORME AU SENS DE K0-1
C 
C     ADINT EST L'ADRESSE DE DEPART DU TABLEAU :
C 
C     HOO-INTERF => NORME AU SENS DE K0
C     SOU-INTERF => NORME AU SENS DE K0-1
C 
C     ADSOU EST L'ADRESSE "COMPLEMENTAIRE"
C 
C     DSOU = AD('SOU-COUCHE') SI ADCOU = AD('HOO-COUCHE')
C     et vice versa, meme chose pour les interfaces
C 
C     On envoie comme arguments :
C 
C     E ...... CARESP   Le caractere * 20 nom du fichier du type deformation
C     E ...... CARSAU   Le caractere * 20 nom du fichier du type saut
C     E ...... NPICDE   piquet de temps de debut d'integration
C     E ...... NPICFI   piquet de temps de fin d'integration
C 
C     Hors de NPICDE et NPICFI la fonction du temps en sortie doit etre nulle
C 
C     Et on recupere :
C 
C     S ...... EPSCH2 et SAUCH4 adresses de depart des tableaux 
C                     des deformations ou contraintes admissibles
C     S ...... FTECH6 Les evolutions en temps correspondantes
C 
      SUBROUTINE RSICBA (TYPE, ADCOU, ADINT, ADSOU, ADSNT,
     &                   CAREPS, CARSAU, NDEBFT, NFINFT, NDEBCH, NFINCH,
     &                   EPSCH2, SAUCH4, FTECH6, NPICDE, NPICFI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'strategie_calcul.h'
C 
      INTEGER TYPE
C 
      INTEGER ADCOU, ADINT, ADSOU, ADSNT, NPICDE, NPICFI
      INTEGER NDEBFT, NFINFT, NDEBCH, NFINCH
C 
      CHARACTER*20  NOM, NOMFIC, CAREPS, CARSAU
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      LOGICAL  LTYP0, LTYP2
C 
      INTEGER  IUNSIG, NUENSI, ADEPS, NUETCO, PGAU1
      INTEGER  DSPCAC, SIGCBA, NUCOU, LONGAU, DSIGCB
      INTEGER  TETA, TEMPS, BSPCAC, LONSCH
C 
C     Pour les interfaces
C 
      INTEGER  NUINT, ADSAUT
      INTEGER  AM2LC, ADM2LC
C 
C     Pour les nouveaux calculs en temps
C 
      INTEGER EPSAPG, EPPADM
      INTEGER FTECH6, EPSCH2, SAUCH4
C 
      DOUBLE PRECISION KLOC(17) , COUISO(10)
      DOUBLE PRECISION INTISO(3)
      DOUBLE PRECISION TETORT, TETCAL
      INTEGER          ADTETA
C 
CD    LOGICAL  LTRACP
C 
      CHARACTER*6 IDPROG
C 
      PARAMETER (IDPROG='RSICBA')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
D     CALL IMPET ('DEBUT DES FONCTIONS DU TEMPS  '//IDPROG, NDEBFT)
D     CALL IMPET ('  FIN DES FONCTIONS DU TEMPS  '//IDPROG, NFINFT)
D     CALL IMPET ('DEBUT DES FONCTIONS EN ESPACE '//IDPROG, NDEBCH)
D     CALL IMPET ('  FIN DES FONCTIONS EN ESPACE '//IDPROG, NFINCH)
C 
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
C 
      LTYP0 = .FALSE.
      IF (TYPE .EQ. 0) LTYP0 = .TRUE.
C 
C     Calcul du numero d'etape locale correspondant pour test
C 
      NUETCO = 1
C 
C     RECHERCHE DES ADRESSES DE TABLEAUX :
C 
C     contraintes admissibles (nsig, nteta, ngau1, nfotps)
C     deformations admissibles (neps, nteta, ngau1, nfotps)
C 
C     accroissements des fonctions du temps donnees (npicet, nbfodo)
C 
C     intervalles de temps (npicet)
C 
      ADEPS = EPSCH2
C 
C      POUR LES INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
        ADSAUT = SAUCH4
C 
       ENDIF
C 
C      Creation d'un tableau provisoire pour ranger toutes les
C      quantites au piquets de temps
C      EPSAPG delta-tild (epsn)                    XDT               X 6*NPICET
C      DSPCAC - (delta sigma point chapeau actuel) XDT               X 6*NPICET
C      SIGCBA -  delta( sigma-point-chapeau) -K( delta-tild( epsn )) X 6*NPICET
C                                                                   _____________
C                                                                    18*NPICET
      CALL POUSMD (18*NPICET, EPSAPG)
C 
      BSPCAC = EPSAPG+6*NPICET
      SIGCBA = BSPCAC+6*NPICET
C 
C     NOMFIC est le nom de fichier dans lequel sont stockees les
C     quantites barre de l'etape locale precedente;
C     ouverture du fichier pour les contraintes
C 
      NOM = CAREPS
      NOMFIC = NOM
      LONSCH = 6*NPICET
      NUENSI  = 1
      CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C 
C     Nombre de points de Gauss dans une couche
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
        IF (LTYP0) THEN
C 
          CALL ANGCOU (NUCOU, TETORT)
          CALL COMISO (ADCOU, NUCOU, COUISO)
C 
        END IF
C 
C       BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
        DO PGAU1  = 1, LONGAU
C 
C         BOUCLE ii SUR LES ANGLES
C 
          DO TETA = 1, NTETA
C 
C           recherche de l'angle de la bande (tetcal) correspondant a teta
C 
            TETCAL = DM(ADTETA+TETA-1)
C 
            IF (LTYP0) CALL RICLOC (COUISO, TETORT, TETCAL, KLOC)
C 
C           On calcule la valeur des deformations admissibles pour tous les piquets
C           de temps pour les deformations provenant de l'etape preliminaire.
C 
            CALL VCPGTE
C 
C           Valeur des Champs au Point de Gauss pour tous les TEmps
C                   on envoie
     &                (CHAMAX, NEPS, NDEBFT, NFINFT, NDEBCH, NFINCH,
     &                 DM(ADEPS), DM(FTECH6),
C                   on recupere
     &                 DM(EPSAPG))
C 
C           Pour aller lire le bon point de Gauss dans les tableaux de stockage
C           des deformations et des contraintes.
C 
            ADEPS = ADEPS + CHARAX*NEPS
C 
C           Lecture des accroissements admissibles de l'etape locale precedente .
C 
            CALL LFDDNF (DM(BSPCAC), BSPCAC, LONSCH, IUNSIG, NUENSI)
C 
C           initialisation des champs avec decalage :
C           EPSAPG <==> EPPADM EPSILON POINT ADMISSIBLE
C           DSPCAC      delta epsilon point chapeau actuel XDT
C           DSIGCB      delta sigma barre actuel XDT
C 
            EPPADM  = EPSAPG
            DSPCAC  = BSPCAC
            DSIGCB  = SIGCBA
C 
C           BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
	    DO  TEMPS = 1, NPICDE-1
C 
C             Decalage d'indice des tableaux pour aller lire au bon endroit dans
C             EPSADM, SIPADM, DEPCAC, DSPCAC, CMTGAC :
C 
              CALL COPITD (6, DM(DSPCAC), DM(DSIGCB))
C 
              EPPADM = EPPADM+6
              DSPCAC = DSPCAC+6
              DSIGCB = DSIGCB+6
C 
            END DO
C 
            DO  TEMPS = NPICDE, NPICFI
C 
              IF (LTYP0) THEN
C 
                CALL CASIBA (DM(EPPADM), KLOC, DM(DSPCAC),
     &                       DM(DSIGCB))
C 
              ELSE
C 
                CALL ADDMAD (6, DM(EPPADM), DM(DSPCAC), DM(DSIGCB))
C 
              END IF
C 
C             Decalage d'indice des tableaux pour aller lire au bon endroit dans
C             EPSADM, SIPADM, DEPCAC, DSPCAC, CMTGAC :
C 
              EPPADM   = EPPADM+6
              DSPCAC   = DSPCAC+6
              DSIGCB   = DSIGCB+6
C 
            END DO
C 
            DO  TEMPS = NPICFI+1, NPICET
C 
C             Decalage d'indice des tableaux pour aller lire au bon endroit dans
C             EPSADM, SIPADM, DEPCAC, DSPCAC, CMTGAC :
C 
              CALL COPITD (6, DM(DSPCAC), DM(DSIGCB))
C 
              EPPADM   = EPPADM+6
              DSPCAC   = DSPCAC+6
              DSIGCB   = DSIGCB+6
C 
C           FIN DE BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
            END DO
C 
C           Pour avoir directement les contraintes permettant le calcul
C           des seconds membres pour l'etape globale on stocke :
C           -delta(sigma-point-chapeau)-K(delta-tild(epsn)).
C 
C           Ecriture des accroissements admissibles de l'etape globale actuelle.
C 
            CALL EFDDNF (IUNSIG, DM(SIGCBA), SIGCBA, LONSCH, NUENSI)
            NUENSI = NUENSI +1
C 
C         FIN DE BOUCLE ii SUR LES ANGLES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES POINTS DE GAUSS
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES COUCHES
C 
      END DO
C 
      CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     CD    CALL IMPTDT ('DENOMINATEUR POUR LES COUCHES '//IDPROG,
C     CD                  DM(DDENOM+1), 1, NPICET)
C 
C     CD    CALL IMPTDT ('NUMERATEUR POUR LES COUCHES '//IDPROG,
C     CD                  DM(DNUMER+1), 1, NPICET)
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
C     DEBUT DE TEST SUR LES INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
C       NOMFIC est le nom de fichier dans lequel sont stockees les
C       quantites barre de l'etape locale precedente;
C       ouverture du fichier pour les contraintes
C 
        NOM    = CARSAU
        NOMFIC = NOM
        LONSCH = 3*NPICET
        NUENSI = 1
C 
        CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C 
C       Remise des poubelles a zero ==> gain de place et FAUTE ?
C 
C       CALL SOPOUB (AM2LC, ADM2LC)
C 
C       creation d'un tableau provisoire pour ranger pour toutes les
C       quantites au piquets de temps
C       EPSAPG delta-tild( saun )                   XDT               X 3*NPICET
C       DSPCAC - (delta sign point chapeau actuel)  XDT               X 3*NPICET
C       SIGCBA -  delta(sign-point-chapeau) - K (delta-tild(saun))    X 3*NPICET
C                                                                   _____________
C                                                                     9*NPICET
        CALL POUSMD (9*NPICET, EPSAPG)
C 
        BSPCAC = EPSAPG+3*NPICET
        SIGCBA = BSPCAC+3*NPICET
        LONGAU = XINTEG*NBCOL
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT = 1, NBINT
C 
          IF (LTYP0) THEN
C 
            CALL IOMISO (ADINT, NUINT, INTISO)
            CALL ANGINT (NUINT, TETORT)
C 
          END IF
C 
C         BOUCLE i SUR LES POINTS DE GAUSS DANS L'INTERFACE
C 
          DO PGAU1  = 1, LONGAU
C 
C           BOUCLE ii SUR LES ANGLES
C 
            DO TETA = 1, NTETA
C 
C             Recherche de l'angle de la bande (tetcal) correspondant a teta
C 
              IF ( LTYP0 ) THEN
C 
                TETCAL = DM(ADTETA+TETA-1)
C 
                CALL INTELA (INTISO, TETORT, TETCAL, KLOC)
C 
              END IF
C 
C             on calcule la valeur des deformations  admissibles pour tous les piquets
C             de temps pour les deformations provenant de l'etape preliminaire.
C 
              CALL VCPGTE
C                     Valeur des Champs au Point de Gauss pour tout les TEmps
C                     on envoie
     &                   (CHAMAX, NSAU, NDEBFT, NFINFT, NDEBCH, NFINCH,
     &                    DM(ADSAUT), DM(FTECH6),
C                     on recupere
     &                    DM(EPSAPG))
C 
C             Pour aller lire le bon point de Gauss dans les tableaux de stockage
C             des deformations et des contraintes.
C 
              ADSAUT = ADSAUT + CHARAX*NSAU
C 
C             Lecture des accroissements admissibles de l'etape locale precedente.
C 
              CALL LFDDNF (DM(BSPCAC), BSPCAC, LONSCH, IUNSIG, NUENSI)
C 
C             Initialisation des champs avec decalage :
C             EPSAPG <==> EPPADM    EPSILON POINT ADMISSIBLE
C             DSPCAC    - delta sign  point chapeau actuel XDT
C             DSIGCB    - delta sign barre          actuel XDT
C 
              EPPADM  = EPSAPG
              DSPCAC  = BSPCAC
              DSIGCB  = SIGCBA
C 
C             BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
              DO  TEMPS = 1, NPICDE-1
C 
C 
C               Decalage d'indice des tableaux pour aller lire au bon endroit dans
C               EPSADM, SIPADM, DEPCAC, DSPCAC, CMTGAC :
C 
                CALL COPITD (3, DM(DSPCAC), DM(DSIGCB))
C 
                EPPADM   = EPPADM+3
                DSPCAC   = DSPCAC+3
                DSIGCB   = DSIGCB+3
C 
              END DO
C 
              DO  TEMPS = NPICDE, NPICFI
C 
                IF (LTYP0) THEN
C 
                  CALL CASNBA (DM(EPPADM), KLOC, DM(DSPCAC), DM(DSIGCB))
C 
                ELSE
C 
                  CALL ADDMAD (3, DM(EPPADM), DM(DSPCAC), DM(DSIGCB))
C 
                END IF
C 
C               Decalage d'indice des tableaux pour aller lire au bon endroit dans
C               EPSADM, SIPADM, DEPCAC, DSPCAC, CMTGAC :
C 
                EPPADM   = EPPADM+3
                DSPCAC   = DSPCAC+3
                DSIGCB   = DSIGCB+3
C 
              END DO
C 
              DO TEMPS = NPICFI+1, NPICET
C 
C               Decalage d'indice des tableaux pour aller lire au bon endroit dans
C               EPSADM, SIPADM, DEPCAC, DSPCAC, CMTGAC :
C 
                CALL COPITD (3,  DM(DSPCAC), DM(DSIGCB))
C 
                EPPADM   = EPPADM+3
                DSPCAC   = DSPCAC+3
                DSIGCB   = DSIGCB+3
C 
C             FIN DE BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
              END DO
C 
C             Pour avoir directement les contraintes permettant le calcul
C             des seconds membres pour l'etape globale on stocke :
C             delta(sgn-point-chapeau)-K(delta-tild(saun );
C             ecriture des accroissements admissibles de l'etape locale actuelle.
C 
              CALL EFDDNF (IUNSIG, DM(SIGCBA), SIGCBA, LONSCH, NUENSI)
              NUENSI = NUENSI +1
C 
C           FIN DE BOUCLE ii SUR LES ANGLES
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES POINTS DE GAUSS DANS L'INTERFACE
C 
          END DO
C 
C       FIN DE BOUCLE SUR LES INTERFACES
C 
        END DO
C 
        CALL FERFIC (3, IUNSIG, IDPROG)
C 
      END IF
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine renvoie la valeur du champ sigmch ou epsich
C     modifiee par l'etape preliminaire soit :
C 
C     sigmch -K0(eppadm) ou epsich -K0-1(eppadm)
C 
C     On envoie comme arguments :
C 
C     E ...... EPPADM    epsilon point admissible (nsig) X DT
C     E ...... CMTGAC    comportement elastique stocke K, B, A, C <=> 9+4+3+1
C     E ...... EPSCAC    -2 X delta epsilon point chapeau actuel
C 
C     Et on recupere :
C 
C     S ...... DSIGCB    -delta sigma point barre actuel
C                        c-a-d : -2delta( sigma-point-chapeau) -K(delta-tild(epsn)).
C 
      SUBROUTINE CASIBA (EPPADM, CMTGAC, EPSCAC, DSIGCB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  EPPADM(6), EPSCAC(6), DSIGCB(6)
      DOUBLE PRECISION  CMTGAC(17)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         I
C 
      DOUBLE PRECISION  EPSLOC( 6 )
C 
      CHARACTER*6 IDPROG
C 
CD    LOGICAL     LTRACN , LTRACP
C 
      PARAMETER (IDPROG='CASIBA')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      CALL MULORT (CMTGAC(1), CMTGAC(10), CMTGAC(14), CMTGAC(17),
     &             EPPADM(1), EPSLOC(1))
      DO I = 1, 6
        DSIGCB(I) =  EPSCAC(I)-EPSLOC(I)
      END DO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... SAPADM    saut point admissible ( nsig ) X DT
C     E ...... CMTGAC    comportement elastique stocke K(9)
C     E ...... DSACAC    -2delta saut point chapeau actuel
C     E ...... CMTGAC    comportement tangent
C 
C     Et on recupere :
C 
C     S ...... DSGNCB    -2delta sigma point barre actuel
C                         c-a-d :  -2delta(sigma-point-chapeau) -K (delta-tild(epsn))
C 
      SUBROUTINE CASNBA (SAPADM, CMTGAC, SAUCAC, DSGNCB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION SAPADM(3) , SAUCAC(3) , DSGNCB(3)
C 
      DOUBLE PRECISION  CMTGAC(9)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         I
C 
      DOUBLE PRECISION SAULOC( 3)
C 
      CHARACTER*6 IDPROG
C 
CD    LOGICAL     LTRACN , LTRACP
C 
       PARAMETER (IDPROG='CASNBA')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACN(1)) THEN
C 
CD      CALL OMPTDN ('MATRICE TANGENTE EN ENTREE ', CMTGAC, 9, 1)
CD      CALL OMPTDN ('SAUT POINT ADMISSIBLE EN ENTREE ',
CD                    SAPADM(1), 3, 1)
CD      CALL OMPTDN ('-2 DELTA SAUT POINT CHAPEAU EN ENTREE ',
CD                    SAUCAC(1), 3, 1)
C 
CD    END IF
C 
      CALL IULORT (CMTGAC(1), SAPADM(1), SAULOC(1))
C 
      DO I = 1, 3
C 
        DSGNCB(I) =  SAUCAC(I)-SAULOC(I)
C 
      END DO
C 
CD    IF (LTRACN(1)) THEN
C 
CD      CALL OMPTDN ( 'DELTA SIGMA POINT BARRE EN SORTIE ',
CD                     DSGNCB(1), 3, 1)
C 
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule pour un point de gauss considere et pour
C     un jeu de valeurs de sigma corespondant les efforts elementaires
C     ranges calcul qui correspondent => pour un element donne il faut
C     sommer sur l'element.
C 
C     On envoie comme arguments :
C 
C     E ...... POIDG   le poids au point de gauss de l'element
C     E ...... TABGAU  tableaux N pour les points de
C                      Gauss concernes (provient de ADR-TGAUSS)
C     E ...... N       le numero de developpement concerne
C     E ...... A       la demi-longeur de l'element
C     E ...... B       la demi-hauteur de l'element
C     E ...... R       le rayon au point de Gauss de l'element
C     E ...... VALSIG  la valeur des contraintes sur l'element
C                      !!!!! ranges calcul
C     Et on recupere :
C 
C     S ...... UN(12)  Valeur de la composante N de l'effort
C                      !!!!!!  U rangee calcul
C     S ...... VN(12)  Valeur de la composante N de l'effort
C                      !!!!!!  V rangee calcul
C     S ...... WN(12)  Valeur de la composante N de l'effort
C                      !!!!!!  W rangee calcul
C                      N > OU = 0  ==> ( COSN , SINN , COSN )
C                      N < 0       ==> ( SINN , COSN , SINN)
C 
      SUBROUTINE BNSIGM (POIDG, TABGAU, N, A, B, R, VALSIG,
     &                   UN, VN, WN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER            N
C 
      DOUBLE PRECISION   TABGAU(36) , A , B , R , VALSIG(6) , UN(12)
      DOUBLE PRECISION   VN(12), WN(12) ,  POIDG
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    LOGICAL           LTRACN , LTRACP
C 
      INTEGER           AD2, DEP
      INTEGER           AM2LC, ADM2LC
C 
      DOUBLE PRECISION  NLOC, R2, OM1, OM2, OM3
C 
      DOUBLE PRECISION  SI1, SI2
C 
      CHARACTER*6 IDPROG
       PARAMETER (IDPROG='BNSIGM')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACP(1) )THEN
C 
CD      CALL IMPDP ( ' POIDG ' , POIDG )
CD      CALL IMPDP ( ' R     ' , R )
CD      CALL IMPDP ( ' A     ' , A )
CD      CALL IMPDP ( ' B     ' , B )
CD      CALL OMPTDP( ' SIGDEV' , VALSIG, 6 , 1 )
C 
CD    END IF
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
      NLOC       = DBLE(N)
      R2         = DSQRT(2.D0)/2.D0
C 
C     Constantes de multiplication
C 
      OM1 = POIDG*R*B
      OM2 = POIDG*A*B
      OM3 = POIDG*R*A
      SI1 = VALSIG(2)-NLOC*R2*VALSIG(3)
      SI2 = NLOC*VALSIG(2)-R2*VALSIG(3)

C 
      CALL POUSMD (36+36 , AD2)
      DEP = AD2+36
C 
      CALL LINLOC( A , 12 , TABGAU(1 ) , DM(AD2)   )
      CALL LINLOC( A , 12 , TABGAU(13) , DM(AD2+12))
      CALL LINLOC( A , 12 , TABGAU(25) , DM(AD2+24))
C 
C     Calcul de FUN
C  
C     CALL MUMARE( VALSIG(1)*OM1       , 12, DM(AD2+12) , DM(DEP   ) )
C     CALL MUMARE( SI1*OM2             , 12, DM(AD2)    , DM(DEP+12) )
C     CALL MUMARE( R2*VALSIG(5)*OM3    , 12, DM(AD2+24) , DM(DEP+24) )
C 
      CALL HOMAT( VALSIG(1)*OM1    , DM(AD2+12) , DM(DEP   ) ,12 , 1 )
      CALL HOMAT( SI1*OM2          , DM(AD2)    , DM(DEP+12) ,12 , 1 )
      CALL HOMAT( R2*VALSIG(5)*OM3 , DM(AD2+24) , DM(DEP+24) ,12 , 1 )
C 
C 
C     CALL ADDMAD( 12, DM(DEP)    , DM(DEP+12)  , DM(DEP+12) )
C     CALL ADDMAD( 12, DM(DEP+12) , DM(DEP+24) , UN(1)       )
C 
      CALL ADD( DM(DEP)    , DM(DEP+12)  , DM(DEP+12) , 12 , 1 )
      CALL ADD( DM(DEP+12) , DM(DEP+24) , UN(1)       , 12 , 1 )
C  
C     Fin du calcul de FUN, debut du calcul de FVN
C  
C     CALL MUMARE( OM2*SI2              , 12, DM(AD2)    , DM(DEP )   )
C     CALL MUMARE( R2*OM1*VALSIG(3)     , 12, DM(AD2+12) , DM(DEP+12) )
C     CALL MUMARE( R2*OM3*VALSIG(4)     , 12, DM(AD2+24) , DM(DEP+24) )
C 
      CALL HOMAT( OM2*SI2          , DM(AD2)    , DM(DEP )   , 12 , 1 )
      CALL HOMAT( R2*OM1*VALSIG(3) , DM(AD2+12) , DM(DEP+12) , 12 , 1 )
      CALL HOMAT( R2*OM3*VALSIG(4) , DM(AD2+24) , DM(DEP+24) , 12 , 1 )
C  
C     CALL ADDMAD( 12, DM(DEP)     , DM(DEP+12)  , DM(DEP+12) )
C     CALL ADDMAD( 12, DM(DEP+12) , DM(DEP+24) ,  VN(1)      )
C 
      CALL ADD( DM(DEP)    , DM(DEP+12) , DM(DEP+12) , 12 , 1 )
      CALL ADD( DM(DEP+12) , DM(DEP+24) , VN(1)      , 12 , 1 )
C  
C     Fin du calcul de FVN, debut du calcul de FWN
C  
C      CALL MUMARE(-R2*NLOC*OM2*VALSIG(4) , 12, DM(AD2   ) , DM(DEP   ) )
C      CALL MUMARE( R2*VALSIG(5)*OM1      , 12, DM(AD2+12) , DM(DEP+12) )
C      CALL MUMARE( VALSIG(6)*OM3         , 12, DM(AD2+24) , DM(DEP+24) )
C 
      CALL HOMAT(-R2*NLOC*OM2*VALSIG(4) , DM(AD2) , DM(DEP) , 12 , 1 )
      CALL HOMAT( R2*VALSIG(5)*OM1 , DM(AD2+12) , DM(DEP+12) , 12 , 1 )
      CALL HOMAT( VALSIG(6)*OM3    , DM(AD2+24) , DM(DEP+24) , 12 , 1 )
C 
C     CALL ADDMAD( 12, DM(DEP ) , DM(DEP +12) , DM(DEP +12) )
C     CALL ADDMAD( 12, DM(DEP+12) , DM(DEP+24) , WN(1)      )
C 
      CALL ADD( DM(DEP ) , DM(DEP +12) , DM(DEP +12) , 12 , 1 )
      CALL ADD( DM(DEP+12) , DM(DEP+24) , WN(1) , 12 , 1 )
C 
CD    IF( LTRACN(1) ) THEN
CD      CALL IMPTDN ('VALEUR DE FUN ', UN(1) , 1 , 12 )
CD      CALL IMPTDN ('VALEUR DE FVN ', VN(1) , 1 , 12 )
CD      CALL IMPTDN ('VALEUR DE FWN ', WN(1) , 1 , 12 )
CD    END IF
C 
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Pour les elements de type couche.
C 
C     Cette routine calcule pour tous les termes du developpement
C     en serie de Fourier pour le point de GAUSS correspondant a
C     TABGAUS, les deplacements elementaires de U, V, W  qui
C     correspondent a NDDLU, NDDLV, NDDLW : ceux-ci  sont directement
C     assembles dans DM a partir de ADSMEG de la facon suivante :
C     (NFTGLO, NDDL, NBMAT) ==> on peut directement utiliser la
C     routine de descente remontee pour plusieurs seconds membres,
C     les seconds membres correspondant au meme rang de Fourier pour
C     des fonctions du temps differentes etant ranges sequentiellement.
C 
C     On envoie comme arguments :
C 
C     E ...... POIDG     le poids au point de gauss de l'element
C     E ...... NFT       le numero de la fonction du temps
C     E....... TABGAU    les tableaux N du point de gauss
C     E ...... A         la demi-longeur de l'element
C     E ...... B         la demi-hauteur de l'element
C     E ...... R         le rayon au point de gauss de l'element
C     E ...... NUDDL     numero des ddl du deplacement ranges
C     E                  calcul pour l'element auquel appartient
C     E                  le point de GAUSS c'est a dire
C     E                  NDDLU, NDDLV, NDDLW calcul
C     E ...... SIGDEV    developpement serie de fourier reel des
C     E                  contraintes integrees sur le temps
C     E                  rangees  composantes des contraintes
C     E                  apres composantes  (6, nbmat)
C     E ...... ADSMEG    l'adresse de depart des seconds membres
C     E                  pour l'etape globale ceux-ci sont
C     E                  assembles dans DM a partir de ADSMEG de
C     E                  la facon suivante (NDDL, NFTGLO, NBMAT)
C 
      SUBROUTINE SMITEC (POIDG, NFT, TABGAU, A, B, R, NUDDL,
     &                   SIGDEV, ADSMEG)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER  NFT, NUDDL(36), ADSMEG
C 
      DOUBLE PRECISION  TABGAU(36), A, B, R, SIGDEV(6*NBMAT), POIDG
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   I, NUDEV, ADEFFO, DEBSIG, DECAL
      DOUBLE PRECISION DEPLA(36)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SMITEC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      DEBSIG = 1
C 
C     Les seconds membres correspondant a NFT sont ranges a partir de
C     ADEFFO.
C 
      ADEFFO = ADSMEG+(NFT-1)*NDDL
C 
C     Calcul du decalage entre deux developpements correspondant a un
C     meme effort puisque l'on range les efforts correspondant a un
C     meme developpement sequentiellement.
C 
      DECAL  = NDDL*NFTGLO
C 
      DO I = 1, NBMAT
C 
        NUDEV = I-NTDSFG-1
C 
        CALL BNSIGM (POIDG, TABGAU, NUDEV, A, B, R,
     &               SIGDEV(DEBSIG), DEPLA(1), DEPLA(13), DEPLA(25))
C 
C       Le decalage du a NUDEV dans ASVEFI est de (NUDEV-NTDSFG)*NDDL;
C       on le met a zero en imposant NUDEV = NTDSFG et on impose le bon
C       decalage a l'aide de ADEFFO.
C 
        CALL ASVEFI (ADEFFO, -NTDSFG, NUDDL, 36, DEPLA(1))
        DEBSIG   =  DEBSIG+6
        ADEFFO   =  ADEFFO+DECAL
C 
      END DO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Pour les elements de type interface
C 
C     Cette routine calcule pour tous les termes du developpement
C     en serie de Fourier pour le point de GAUSS correspondant a
C     TABGAUS les deplacements elementaires de U, V, W qui
C     correspondent a NDDLU, NDDLV, NDDLW; ceux-ci  sont directement
C     assembles dans DM a partir de ADSMEG de la facon suivante :
C     (NFTGLO, NDDL, NBMAT  ==> on peut directement utiliser la
C     routine de descente remontee pour plusieurs seconds membres;
C     les seconds membres correspondant au meme rang de Fourier pour
C     des fonctions du temps differentes etant ranges sequentiellement.
C 
C     On envoie comme arguments :
C 
C     E ...... POIDG   Le poids au point de gauss de l'element
C     E ...... NFT     le numero de la fonction du temps
C     E ...... TABGAU  les tableaux N du point de gauss
C     E ......  A      La demi-longeur de l'element
C     E ...... R       Le rayon au point de gauss de l'element
C     E ...... NUDDL   numero des ddl du deplacement ranges
C                      calcul pour l'element auquel appartient
C                      le point de GAUSS c'est a dire NDDLU, NDDLV, NDDLW calcul
C     E ...... SINDEV  developpement serie de Fourier reel des
C                      contraintes normales integrees sur le temps
C                      rangees  composantes des contraintes
C                      apres composantes (6, nbmat)
C     E ...... ADSMEG  l'adresse de depart des seconds membres
C                      pour l'etape globale ceux-ci sont
C                      assembles dans DM a partir de ADSMEG de
C                      la facon suivante (NDDL, NFTGLO, NBMAT)
C 
      SUBROUTINE SMITEI (POIDG, NFT, TABGAU, A, R, NUDDL,
     &                   SINDEV, ADSMEG)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER  NFT , NUDDL(24) , ADSMEG
C 
      DOUBLE PRECISION  TABGAU(8) , A ,  R , SINDEV( 3*NBMAT)
C 
      DOUBLE PRECISION POIDG
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   I , NUDEV ,  ADEFFO , DEBSIG  , DECAL
C 
      DOUBLE PRECISION DEPLA(24)
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SMITEI')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      DEBSIG = 1
C 
C     Les seconds membres correspondant a NFT sont ranges a partir de
C     ADEFFO.
C 
      ADEFFO = ADSMEG+(NFT-1)*NDDL
C 
C     Calcul du decalage entre deux developpements correspondant a un
C     meme effort puisque l'on range les efforts correspondant a un
C     meme developpement sequentiellement.
C 
      DECAL  = NDDL*NFTGLO
C 
      DO I =1 , NBMAT
C 
        NUDEV = I-NTDSFG-1
C 
CD      CALL IMPEN( ' Pour le developpement ' , NUDEV )
C 
        CALL BNSIGN (POIDG, TABGAU, NUDEV, A, R, SINDEV (DEBSIG),
     &               DEPLA(1), DEPLA(9), DEPLA(17))
C 
C     Le decalage du a NUDEV dans ASVEFI est de (NUDEV+NTDSFG)*NDDL.
C     On le met a zero en imposant NUDEV = NTDSFG et on impose le bon
C     decalage a l'aide de ADEFFO.
C 
        CALL ASVEFI( ADEFFO, -NTDSFG , NUDDL , 24 , DEPLA(1) )
        DEBSIG   =  DEBSIG+3
        ADEFFO   =  ADEFFO+DECAL
C 
      END DO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule pour un point de gauss considere et pour
C     un jeu de valeurs de sigma correspondant les deplacements elementaires
C     ranges calcul qui correspondent => pour un element donne il faut
C     sommer sur l'element.
C 
C     On envoie comme arguments :
C 
C     E ...... POIDG   Le poids au point de gauss de l'element
C     E ...... TABGAU  Tableaux N pour les points de
C                      Gauss concernes (provient de ADR-TGAUSS)
C     E ......  N      Le numero de developpement concerne
C     E ...... A       La demi-longeur de l'element
C     E ...... R       Le rayon au point de gauss de l'element
C     E ...... VALSIN  La valeur des contraintes sur l'element
C                      !!!!!    ranges calcul (developpees Fourier)
C 
C     Et on recupere :
C 
C     S ...... UN(8)   Valeur de la composante N de l'effort
C                      !!!!!!  U rangee calcul
C     S ...... VN(8)   Valeur de la composante N de l'effort
C                      !!!!!!  V rangee calcul
C     S ...... WN(8)   Valeur de la composante N de l'effort
C                      !!!!!!  W rangee calcul
C                      N > OU = 0  ==> ( COSN , SINN , COSN )
C                      N < 0       ==> ( SINN , COSN , SINN)
C 
      SUBROUTINE BNSIGN (POIDG, TABGAU, N, A, R,  VALSIN, UN, VN, WN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER            N
C 
      DOUBLE PRECISION   TABGAU(8), A, R, VALSIN(3), UN(8)
      DOUBLE PRECISION   VN(8), WN(8), POIDG
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    LOGICAL LTRACN , LTRACP
C 
      INTEGER     AD2
      INTEGER     AM2LC, ADM2LC
C 
      DOUBLE PRECISION   MULT
C 
      CHARACTER*6 IDPROG
       PARAMETER (IDPROG='BNSIGN')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL POUSMD (8 , AD2)
      CALL LINLOC( A , 8 , TABGAU(1 ) , DM(AD2)   )
      MULT = A*R*POIDG
C 
C     TCD CALL MUMARE( VALSIN(1)*MULT  , 8, DM(AD2) , UN(1) )
C     TCD CALL MUMARE( VALSIN(2)*MULT  , 8, DM(AD2) , VN(1) )
C     TCD CALL MUMARE( VALSIN(3)*MULT  , 8, DM(AD2) , WN(1) )
C 
      CALL HOMAT( VALSIN(1)*MULT  , DM(AD2) , UN(1) , 8 , 1 )
      CALL HOMAT( VALSIN(2)*MULT  , DM(AD2) , VN(1) , 8 , 1 )
      CALL HOMAT( VALSIN(3)*MULT  , DM(AD2) , WN(1) , 8 , 1 )
C  
CD    IF( LTRACN(1) ) THEN
CD      CALL IMPTDN ('VALEUR DE FUN ', UN(1) , 1 , 8 )
CD      CALL IMPTDN ('VALEUR DE FVN ', VN(1) , 1 , 8 )
CD      CALL IMPTDN ('VALEUR DE FWN ', WN(1) , 1 , 8 )
CD    END IF
C  
CD    IF( LTRACP(1) ) THEN
CD      CALL IMPTDP ('VALEUR DE SIGMA en entree ', VALSIN(1) , 3 , 1 )
CD    END IF
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C  
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine orthogonalise les nnouv nouveaux vecteurs
C     par rapport aux anciens (orthogonalises) et a eux meme.
C     On modifie nnouv si les vecteurss sont trop proches.
C 
C     ADCOU EST L'ADRESSE DE DEPART DU TABLEAU :
C 
C       HOO-COUCHE => NORME AU SENS DE K0
C       SOU-COUCHE => NORME AU SENS DE K0-1
C 
C     ADINT EST L'ADRESSE DE DEPART DU TABLEAU :
C 
C       HOO-INTERF => NORME AU SENS DE K0
C       SOU-INTERF => NORME AU SENS DE K0-1
C 
C     L'OTHOGONALISATION EST A PRENDRE AU SENS DE K0 OU DE K0-1
C 
C     On envoie comme arguments :
C 
C     E ...... NBANC    nombre de vecteur des iterations
C     E                          PRECEDENTE
C     E ...... NNOUV    nombre de nouveaux vecteurs
C     E ...... LONEPS   la longueur d'une deformation
C     E ...... LONSAU   la longueur d'un saut
C     ES...... VECANC   les anciens vecteurs orthogonalises
C     ES...... VECNOU   les nouveaux vecteurs a orthogonaliser
C 
C     Et on recupere :
C 
C     ES...... NNOUVA   nombre de nouveaux vecteurs apres orthogonalisation
C     ES...... VECORT   NNOUVA (sortie) vecteurs orthogonalises
C     S ...... COEFF    Coefficient d'orthogonalisation
C 
      SUBROUTINE ORNEK0 (ADCOU, ADINT, NBANC, NNOUV, LONEPS, LONSAU,
     &                   EPSANC, SAUANC, EPSNOU, SAUNOU,
     &                   EPSORT, SAUORT, NNOUVA, COEFF)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER       ADCOU, ADINT
      INTEGER       LONEPS, LONSAU, NBANC, NNOUV, NNOUVA
C 
      DOUBLE PRECISION  EPSANC(LONEPS*NBANC),  EPSNOU(LONEPS*NNOUV)
      DOUBLE PRECISION  SAUANC(LONSAU*NBANC),  SAUNOU(LONSAU*NNOUV)
      DOUBLE PRECISION  EPSORT(LONEPS*NNOUV), COEFF(NBANC+1)
      DOUBLE PRECISION  SAUORT(LONSAU*NNOUV)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I, K
      INTEGER  J, LONREE, LONRES
      INTEGER  DEBUE, DEBUS, JDEBUE, KDEBUE, JDEBUS, KDEBUS
      INTEGER  PLAC, DBEPS, DBSAU
C 
      DOUBLE PRECISION NORMF, PSCALV, SEUIL, MULT, DIVIS
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ORNEK0')
      PARAMETER (SEUIL = 1.D -10)
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      DBEPS = 1
      DBSAU = 1
C 
      DO I = 1, NNOUV
C 
        CALL NOKGLO (ADCOU, ADINT, 1, 1,
     &               EPSNOU(DBEPS), SAUNOU(DBSAU), NORMF)
C 
        NORMF = DSQRT(NORMF)
        DIVIS = 1.D0/NORMF
C 
        CALL MUMARE (DIVIS, LONEPS, EPSNOU(DBEPS), EPSORT(DBEPS))
C 
        CALL MUMARE (DIVIS, LONSAU, SAUNOU(DBSAU), SAUORT(DBSAU))
C 
        DBSAU = DBSAU+LONSAU
        DBEPS = DBEPS+LONEPS
C 
      END DO
C 
      MULT = NORMF
C 
      NNOUVA =  NNOUV
C 
      J = 1
C 
C     BOUCLE POUR PASSER EN REVUE TOUS LES VECTEURS
C 
      DO WHILE (J .LE.  NNOUVA)
C 
C       Orthogonalisation du (j+1)ieme vecteur
C 
        JDEBUE = 1 + (J-1)*LONEPS
        JDEBUS = 1 + (J-1)*LONSAU
        KDEBUE = 1
        KDEBUS = 1
C 
C       BOUCLE i POUR
C       Remise de debut au debut du jieme vecteur du temps <=> jdebut
C 
        DO K = 1, NBANC
C 
          DEBUE = JDEBUE
          DEBUS = JDEBUS
C 
C         Calcul du produit scalaire de fnouv(j) avec ft(k)
C         dont l'adresse de depart est KDEBU.
C 
          CALL SCKGLO (ADCOU, ADINT, 1, 1,
     &                 EPSORT(JDEBUE), SAUORT(JDEBUS), 1,
     &                 EPSANC(KDEBUE), SAUANC(KDEBUS), PSCALV)
C 
          COEFF(K) = PSCALV
C 
C         Modification du (j+1)ieme vecteur comme suit
C         vect(j+1) = vect(j+1) - scalt (vect(j+1), vect(k))* vect(k)
C 
          DO PLAC = 1, LONEPS
C 
            EPSORT(DEBUE) = EPSORT(DEBUE) - PSCALV*EPSANC(KDEBUE)
            DEBUE         = DEBUE+1
            KDEBUE        = KDEBUE+1
C 
          END DO
C 
          DO PLAC = 1, LONSAU
C 
            SAUORT(DEBUS) = SAUORT(DEBUS) - PSCALV*SAUANC(KDEBUS)
            DEBUS         = DEBUS+1
            KDEBUS        = KDEBUS+1
C 
          END DO
C 
C       FIN DE BOUCLE i POUR
C       Remise de debut au debut du jieme vecteur du temps <=> jdebut
C 
        END DO
C 
        DO K = 1, J-1
C 
C       BOUCLE i POUR
C       Remise de debut au debut du jieme vecteur du temps <=> jdebut
C 
          DEBUE  = JDEBUE
          KDEBUE = 1
C 
          DEBUS  = JDEBUS
          KDEBUS = 1
C 
C         Calcul du produit scalaire de ft(j+1)  avec ft(k)
C         dont l'adresse de depart est KDEBUT
C 
          CALL SCKGLO (ADCOU, ADINT, 1, 1,
     &                 EPSORT(JDEBUE), SAUORT(JDEBUS), 1,
     &                 EPSORT(KDEBUE), SAUORT(KDEBUS), PSCALV)
C 
C         Modification du (j+1)ieme vecteur comme suit
C         vect(j+1) = vect(j+1) - scalt ( vect(j+1) , vect(k) )* vect(k)
C 
          DO PLAC = 1, LONEPS
C 
            EPSORT(DEBUE) = EPSORT(DEBUE) - PSCALV*EPSORT(KDEBUE)
            DEBUE         = DEBUE+1
            KDEBUE        = KDEBUE+1
C 
          END DO
C 
          DO PLAC = 1, LONSAU
C 
            SAUORT(DEBUS) = SAUORT(DEBUS) - PSCALV*SAUORT(KDEBUS)
            DEBUS         = DEBUS+1
            KDEBUS        = KDEBUS+1
C 
          END DO
C 
C       FIN DE BOUCLE i POUR
C       Remise de debut au debut du jieme vecteur du temps <=> jdebut
C 
        END DO
C 
        IF (NBANC .GT. 0) THEN
C 
          CALL NOKGLO (ADCOU, ADINT, 1, 1,
     &                 EPSORT(JDEBUE), SAUORT(JDEBUS),
     &                 NORMF)
C 
          NORMF = DSQRT(NORMF)
          IF (NORMF .LT. SEUIL)	goto 1001  
C 
          DIVIS = 1.D0/NORMF
C 
          CALL MUMARE (DIVIS, LONEPS, EPSORT(JDEBUE), EPSORT(JDEBUE))
C 
          CALL MUMARE (DIVIS, LONSAU, SAUORT(JDEBUS), SAUORT(JDEBUS))
C 
        ELSE
C 
          NORMF = 1.D0
C 
        END IF
C 
1001    CONTINUE
        IF (NORMF .LT. SEUIL) THEN
C 
          CALL IMPDT
     &    ('NORME TROP PETITE POUR CONSERVER LE VECTEUR : NORME '
     &    , NORMF)
C 
C         On retourne a la meme valeur de j apres avoir recopie les valeurs
C         restantes sur la longueur lonrec des fonctions a partir de la jieme fonction
C 
          LONREE = (NNOUVA-J)*LONEPS
          LONRES = (NNOUVA-J)*LONSAU
C 
          IF (LONREE .GT. 0) THEN
            CALL COPITD (LONREE, EPSORT(DEBUE), EPSORT(JDEBUE))
          END IF
C 
          IF (LONRES .GT. 0) THEN
            CALL COPITD (LONRES, SAUORT(DEBUS), SAUORT(JDEBUS))
          END IF
C 
C         Decrementation du nombre de vecteur
C 
          NNOUVA  = NNOUVA-1
          J       = J-1
C 
C         SEQUENCE DE MODIFICATION ???????????????
C 
C         DO I = 1, NBANC +1
C 
C           COEFF(I) = 0.D0
C 
C         END DO
C 
C         ELSE
C 
C           MULT = MULT*NORMF
C           COEFF( NBANC+1) = MULT
C 
        END IF
C 
        COEFF( NBANC+1) = MULT*NORMF
C 
        DO  I =1, NBANC
C 
          COEFF( I) = MULT*COEFF(I)
C 
        END DO
C 
C     END IF
C 
        J = J+1
C 
C     FIN DE LA BOUCLE POUR PASSER EN REVUE TOUS LES VECTEURS
C 
      END DO
C 
CD      CALL MESSAO ('SEQUENCE DE VERIFICATION DANS '//IDPROG )
C 
CD      CALL IMPET ('NOMBRE DE NOUVEAUX VECTEURS INITIAUX ', NNOUV)
CD      CALL IMPET ('NOMBRE DE NOUVEAUX VECTEURS FINAUX   ', NNOUVA)
C 
CD      DO  J = 1, NNOUVA
C 
CD        CALL IMPET ('POUR LE NOUVEAU VECTEUR NUMERO ', J)
C 
CD        JDEBUE = 1+(J-1)*LONEPS
CD        JDEBUS = 1+(J-1)*LONSAU
C 
CD        CALL NOKGLO (ADCOU, ADINT, 1, 1,
CD                     EPSORT(JDEBUE), SAUORT(JDEBUS), NORMF)
C 
CD        CALL IMPDT ('LA NORME VAUT ', NORMF)
C 
CD        DO I = 1, NBANC
C 
CD          CALL IMPET ('POUR L''ANCIEN VECTEUR NUMERO ', I)
C 
CD          KDEBUE = 1+(I-1)*LONEPS
CD          KDEBUS = 1+(I-1)*LONSAU
C 
CD          CALL SCKGLO (ADCOU, ADINT, 1, 1,
CD                       EPSORT(JDEBUE), SAUORT(JDEBUS), 1,
CD                       EPSANC(KDEBUE), SAUANC(KDEBUS), PSCALV)
C 
CD          CALL IMPDT ('LE PRODUIT SCALAIRE VAUT ', PSCALV)
C 
CD        END DO
C 
CD     END DO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine orthogonalise les NNOUV nouveaux vecteurs
C     par rapport aux anciens (orthogonalises) et a eux-memes.
C     On modifie NNOUV si les vecteurs sont trop proches.
C 
C     L'OTHOGONALISATION EST A PRENDRE AU SENS DE K0 OU DE K0-1
C 
C     SAUF QUE :
C 
C     LES DEPLACEMENTS SONT AUSSI ORTHOGONALISES AU SENS DE K0
C 
C     ADCOU EST L'ADRESSE DE DEPART DU TABLEAU :
C 
C     HOO-COUCHE => NORME AU SENS DE K0
C     SOU-COUCHE => NORME AU SENS DE K0-1
C 
C     ADINT EST L'ADRESSE DE DEPART DU TABLEAU :
C 
C     HOO-INTERF => NORME AU SENS DE K0
C     SOU-INTERF => NORME AU SENS DE K0-1
C 
C     On envoie comme arguments :
C 
C     E ...... NBANC    nombre de vecteurs des iterations precedentes
C     E....... NNOUV    nombre de nouveaux vecteurs
C     E....... LONEPS   la longueur d'une deformation
C     E....... LONSAU   la longueur d'un saut
C     ES...... VECANC   les anciens vecteurs orthogonalises
C     ES...... VECNOU   les nouveaux vecteurs a orthogonaliser
C 
C     Et on recupere :
C 
C     ES...... NNOUVA   nombre de nouveaux vecteurs apres orthogonalisation
C     ES...... VECORT   NNOUVA vecteurs orthogonalises
C     S ...... COEFF    Coefficient d'orthogonalisation
C 
      SUBROUTINE OREDK0 (ADCOU, ADINT, NBANC, NNOUV,
     &                   LONEPS, LONSAU, LONDEP,
     &                   EPSANC, SAUANC, DEPANC, EPSNOU, SAUNOU, DEPNOU,
     &                   EPSORT, SAUORT, DEPORT, NNOUVA, COEFF)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER  ADCOU, ADINT
      INTEGER  LONEPS, LONSAU, LONDEP, NBANC, NNOUV, NNOUVA
C 
      DOUBLE PRECISION  EPSANC(LONEPS*NBANC), EPSNOU(LONEPS*NNOUV)
      DOUBLE PRECISION  SAUANC(LONSAU*NBANC), SAUNOU(LONSAU*NNOUV)
      DOUBLE PRECISION  DEPANC(LONDEP*NBANC), DEPNOU(LONDEP*NNOUV)
C 
      DOUBLE PRECISION  EPSORT(LONEPS*NNOUV), COEFF(NBANC+1)
      DOUBLE PRECISION  SAUORT(LONSAU*NNOUV)
      DOUBLE PRECISION  DEPORT(LONDEP*NNOUV)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I, K
      INTEGER  J, LONREE, LONRES, LONRED
      INTEGER  DEBUE, DEBUS, JDEBUE, KDEBUE, JDEBUS, KDEBUS
      INTEGER  DEBUD, JDEBUD, KDEBUD
      INTEGER  PLAC, DBEPS, DBSAU, DBDEP
C 
      DOUBLE PRECISION NORMF, PSCALV, SEUIL, MULT, DIVIS
      PARAMETER (SEUIL = 1.D -10)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='OREDK0')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      DBEPS = 1
      DBSAU = 1
      DBDEP = 1
C 
      DO I = 1, NNOUV
C 
        CALL NOKGLO (ADCOU, ADINT, 1, 1,
     &               EPSNOU(DBEPS), SAUNOU(DBSAU), NORMF)
C 
        NORMF = DSQRT(NORMF)
C 
        DIVIS = 1.D0/NORMF
C 
        CALL MUMARE (DIVIS, LONEPS, EPSNOU(DBEPS), EPSORT(DBEPS))
        CALL MUMARE (DIVIS, LONSAU, SAUNOU(DBSAU), SAUORT(DBSAU))
        CALL MUMARE (DIVIS, LONDEP, DEPNOU(DBDEP), DEPORT(DBDEP))
C 
        DBEPS = DBEPS+LONEPS
        DBSAU = DBSAU+LONSAU
        DBDEP = DBDEP+LONDEP
C 
      END DO
C 
      MULT = NORMF
C 
      NNOUVA =  NNOUV
C 
      J = 1
C 
C     BOUCLE POUR PASSER EN REVUE TOUS LES VECTEURS
C 
      DO WHILE (J .LE. NNOUVA)
C 
C       Orthogonalisation du (j+1)ieme vecteur
C 
        JDEBUE  = 1 + (J-1)*LONEPS
        JDEBUS  = 1 + (J-1)*LONSAU
        JDEBUD  = 1 + (J-1)*LONDEP
        KDEBUE  = 1
        KDEBUS  = 1
        KDEBUD  = 1
C 
C       remise de debut au debut du jieme vecteur du temps <=> jdebut
C 
        DO K = 1, NBANC
C 
          DEBUE = JDEBUE
          DEBUS = JDEBUS
          DEBUD = JDEBUD
C 
C         Calcul du produit scalaire de fnouv(j) avec ft(k)
C         dont l'adresse de depart est KDEBU
C 
          CALL SCKGLO (ADCOU, ADINT, 1, 1,
     &                 EPSORT(JDEBUE), SAUORT(JDEBUS), 1,
     &                 EPSANC(KDEBUE), SAUANC(KDEBUS), PSCALV)
C 
          COEFF( K) = PSCALV
C 
C         Modification du (j+1)ieme vecteur comme suit
C         vect(j+1) = vect(j+1) - scalt (vect(j+1), vect(k))* vect(k)
C 
          DO PLAC = 1, LONEPS
C 
            EPSORT(DEBUE) = EPSORT(DEBUE) - PSCALV*EPSANC(KDEBUE)
            DEBUE         = DEBUE+1
            KDEBUE        = KDEBUE+1
C 
          END DO
C 
          DO PLAC = 1, LONSAU
C 
            SAUORT(DEBUS) = SAUORT(DEBUS) - PSCALV*SAUANC(KDEBUS)
            DEBUS         = DEBUS+1
            KDEBUS        = KDEBUS+1
C 
          END DO
C 
          DO PLAC = 1, LONDEP
C 
            DEPORT(DEBUD) = DEPORT(DEBUD) - PSCALV*DEPANC(KDEBUD)
            DEBUD = DEBUD +1
            KDEBUD = KDEBUD +1
C 
          END DO
C 
C       remise de debut au debut du jieme vecteur du temps <=> jdebut
C 
        END DO
C 
C       remise de debut au debut du jieme vecteur du temps <=> jdebut
C 
        DO K = 1 , J-1
C 
          DEBUE  = JDEBUE
          KDEBUE = 1
C 
          DEBUS  = JDEBUS
          KDEBUS = 1
C 
          DEBUD = JDEBUD
          KDEBUD = 1
C 
C         Calcul du produit scalaire de ft(j+1) avec ft(k)
C         dont l'adresse de depart est KDEBUT
C 
          CALL SCKGLO (ADCOU, ADINT, 1, 1,
     &                 EPSORT(JDEBUE), SAUORT(JDEBUS), 1,
     &                 EPSORT(KDEBUE), SAUORT(KDEBUS), PSCALV)
C 
C          Modification du (j+1)ieme vecteur comme suit
C          vect(j+1) = vect(j+1) - scalt ( vect(j+1) , vect(k) )* vect(k)
C 
          DO PLAC = 1, LONEPS
C 
            EPSORT(DEBUE) = EPSORT(DEBUE) - PSCALV*EPSORT(KDEBUE)
            DEBUE         = DEBUE+1
            KDEBUE        = KDEBUE+1
C 
          END DO
C 
          DO PLAC = 1, LONSAU
C 
            SAUORT(DEBUS) = SAUORT(DEBUS) - PSCALV*SAUORT(KDEBUS)
            DEBUS         = DEBUS+1
            KDEBUS        = KDEBUS+1
C 
          END DO
C 
          DO PLAC = 1, LONDEP
C 
            DEPORT(DEBUD) = DEPORT(DEBUD) - PSCALV*DEPORT(KDEBUD)
            DEBUD = DEBUD +1
            KDEBUD = KDEBUD +1
C 
          END DO
C 
C       remise de debut au debut du jieme vecteur du temps <=> jdebut
C 
        END DO
C 
C       NORMALISATION DU VECTEUR
C 
        IF (NBANC .GT. 0) THEN
C 
          CALL NOKGLO (ADCOU, ADINT, 1, 1,
     &                 EPSORT(JDEBUE), SAUORT(JDEBUS), NORMF)
C 
          NORMF = DSQRT( NORMF )
          DIVIS = 1.D0/NORMF
C 
          CALL MUMARE (DIVIS, LONEPS, EPSORT(JDEBUE), EPSORT(JDEBUE))
          CALL MUMARE (DIVIS, LONSAU, SAUORT(JDEBUS), SAUORT(JDEBUS))
          CALL MUMARE (DIVIS, LONDEP, DEPORT(JDEBUD), DEPORT(JDEBUD))
C 
        ELSE
C 
          NORMF = 1.D0
C 
        END IF
C 
CD      CALL IMPET ('POUR LE NOUVEAU VECTEUR CONSERVE ', J)
CD      CALL IMPDT ('VALEUR DE LA NORME DU NOUVEAU VECTEUR ', NORMF)
C 
        IF (NORMF .LT. SEUIL) THEN
C 
C 
          CALL IMPDT ('NORME TROP PETITE POUR CONSERVER LE VECTEUR : '//
     &                'NORME ', NORMF)
C 
C         On retourne a la meme valeur de j apres avoir recopie les valeurs
C         restantes sur la longueur lonrec des fonctions a partir de la jieme fonction
C 
          LONREE = (NNOUVA-J)*LONEPS
          LONRES = (NNOUVA-J)*LONSAU
          LONRED = (NNOUVA-J)*LONDEP
C 
          IF (LONREE .GT. 0) THEN
            CALL COPITD (LONREE, EPSORT(DEBUE), EPSORT(JDEBUE))
            CALL COPITD (LONRED, DEPORT(DEBUD), DEPORT(JDEBUD))
          END IF
C 
          IF (LONRES .GT. 0) THEN
            CALL COPITD (LONRES, SAUORT(DEBUS), SAUORT(JDEBUS))
          END IF
C 
C         Decrementation du nombre de vecteurs
C 
          NNOUVA  = NNOUVA-1
          J       = J-1
C 
        END IF
C 
        COEFF(NBANC+1) = MULT*NORMF
C 
        DO  I =1, NBANC
C 
          COEFF(I) = MULT*COEFF(I)
C 
        END DO
C 
        J = J+1
C 
C     FIN DE LA BOUCLE POUR PASSER EN REVUE TOUS LES VECTEURS
C 
      END DO
C 
D        CALL MESSAO ('SEQUENCE DE VERIFICATION DANS        '//IDPROG)
D        CALL IMPET  ('NOMBRE DE NOUVEAUX VECTEURS INITIAUX ', NNOUV)
D        CALL IMPET  ('NOMBRE DE NOUVEAUX VECTEURS FINAUX   ', NNOUVA)
C 
D        DO  J = 1, NNOUVA
C 
D          CALL IMPET ('POUR LE NOUVEAU VECTEUR NUMERO ', J)
C 
D          JDEBUE = 1+(J-1)*LONEPS
D          JDEBUS = 1+(J-1)*LONSAU
C 
D          CALL NOKGLO (ADCOU, ADINT, 1, 1,
D    &                  EPSORT(JDEBUE), SAUORT(JDEBUS), NORMF)
C 
D          CALL IMPDT ('LA NORME VAUT ', NORMF)
C 
D          DO I = 1, NBANC
C 
D            CALL IMPET ('POUR L''ANCIEN VECTEUR NUMERO ', I)
C 
D            KDEBUE = 1+(I-1)*LONEPS
D            KDEBUS = 1+(I-1)*LONSAU
C 
D           CALL SCKGLO (ADCOU, ADINT, 1, 1,
D    &                   EPSORT(JDEBUE), SAUORT(JDEBUS), 1,
D    &                   EPSANC(KDEBUE), SAUANC(KDEBUS), PSCALV)
C 
D            CALL IMPDT ('LE PRODUIT SCALAIRE VAUT ', PSCALV)
C 
D         END DO
C 
D       END DO
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine effectue le produit scalaire en espace des NFONCT
C     deformations admissibles a zero ALPHAI envoyees stockees de 1 a
C     nfonct par les contraintes STOCKEE DANS SID.
C 
C     Le resultat de cette routine est un tableau de NFONCT
C     fonctions rangees sequentiellement.
C 
C     On envoie comme arguments :
C 
C     E ...... TYPE     <=> 0 on fait le produit normal
C     E ...... TYPE     <=> 1 on fait le produit par Ke
C     E ...... TYPE     <=> 2 on fait le produit par Ke-1
C     E ...... SIDEPS   character*6 = epsich ou sigmch suivant les besoins
C     E ...... SIDSAU   character*6 = sautch ou sinoch suivant les besoins
C     E ...... NFONCT   le nombre de fonctions du temps a orthogonaliser
C     E ...... ALPHAI   Deformations admissibles a zero rangees (NEPS, NTETA, NGAU1, NFONCT)
C     E ...... SAUTI    Sauts admissibles a zero ranges (NSAU, NTETA, NGAU2, NFONCT)
C     E....... NPICDE   piquet de temps de debut d'integration
C     E....... NPICFI   piquet de temps de fin d'integration
C 
C     Hors de NPICDE et NPICFI la fonction du temps en sortie doit etre nulle.
C 
C     Et on recupere :
C 
C     S ...... SACTEP   les fonctions du temps INTESP(-delta( sign )*ALPHAI)
C                       rangees (NPICET, NFONCT)
C 
      SUBROUTINE SCHAST (TYPE, SIDEPS, SIDSAU, NFONCT,
     &                   ALPHAI, SAUTI, NPICDE, NPICFI, SACTEP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      CHARACTER*10 SIDEPS, SIDSAU
C 
      INTEGER  NFONCT, TYPE
C 
      DOUBLE PRECISION  SACTEP(NPICET*NFONCT)
      DOUBLE PRECISION  ALPHAI(NEPS*NTETA*NGAU1*NFONCT)
      DOUBLE PRECISION  SAUTI(NSAU*NTETA*NGAU2*NFONCT)
      DOUBLE PRECISION  RAYONC, A, JLOC, RLOC
C 
      INTEGER NPICDE, NPICFI
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      LOGICAL LTYP0, LTYP12
C 
C     Pour le traitement du cas 12
C 
      DOUBLE PRECISION KLOC(17), COUISO(10)
      DOUBLE PRECISION KLOCI(9), INTISO(3)
      DOUBLE PRECISION TETORT, TETCAL
      INTEGER          ADTETA
C 
      INTEGER ADCOU, ADINT
C 
      DOUBLE PRECISION TRACE, B
      DOUBLE PRECISION MULTI, MULTJ
C 
      INTEGER     EPSILO, SAUTLO, SAUI, EPSI
      INTEGER     LONEPS, LONSAU
      INTEGER     SIGMAC, SIGNAC, NUCOU
      INTEGER     TETA, NUCOL, BMTGAC
      INTEGER     DEBUT, I, TEMPS, NUINT
      INTEGER     AM2LC, ADM2LC, HINT, KINT
      INTEGER     IUNSIG, LONSCH, NUENSI
      INTEGER     FONTET, TETTOT, LONTOT, ADTET
      INTEGER     AVTET, X, Y
      INTEGER     EPSLOC, DBCPG, LONLIF
C 
      INTEGER ADPOPG, DBPOPG
C 
C     CD    INTEGER     ADTEI , AVTEI , SACTEI
C 
CD    LOGICAL     LTRACN, LTRACP
C 
      CHARACTER*6  IDPROG
      CHARACTER*10 NOMFIC
      PARAMETER (IDPROG='SCHAST')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACP(1)) THEN
C 
CD      CALL IMPTDP ('ALPHAI ', ALPHAI(1), NEPS*NTETA*NGAU1, NFONCT)
CD      CALL IMPTDP ('SAUTI ', SAUTI(1), NSAU*NTETA*NGAU2, NFONCT)
C 
CD    END IF
C 
      LTYP0  = .FALSE.
      LTYP12 = .FALSE.
C 
      IF (TYPE .EQ. 0) LTYP0 = .TRUE.
C 
      IF (TYPE .EQ. 1) THEN
         LTYP12 = .TRUE.
         CALL ADTBDM ('HOO-COUCHE', ADCOU)
         IF (NBINT .GT. 0) THEN
           CALL ADTBDM ('HOO-INTERF', ADINT)
         ELSE
           ADINT = ADCOU
         END IF
      END IF
      IF (TYPE .EQ. 2) THEN
         LTYP12 = .TRUE.
         CALL ADTBDM ('SOU-COUCHE', ADCOU)
         IF (NBINT .GT. 0) THEN
           CALL ADTBDM ('SOU-INTERF', ADINT)
         ELSE
           ADINT = ADCOU
         END IF
      END IF
C 
      IF (LTYP12) THEN
        CALL ADTBDM ('ANGLES-GEO', ADTETA )
      ENDIF
C 
C     Pour aller lire les points de Gauss
C 
      HINT = XINTEG*(XINTEG-1)/2
      KINT = YINTEG*(YINTEG-1)/2
C 
      NOMFIC = SIDEPS
C 
C     NOMFIC est le nom de fichier dans lequel sont stockees les quantites barre
C     de l'etape locale precedente. Ouverture du fichier pour les contraintes.
C 
      LONSCH = 6*NPICET
      NUENSI  = 1
      CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C 
C     Longueur du stockage pour tous les angles d'une fonction
C 
      FONTET = NTETA*NPICET
C 
C     Longueur totale du tableau partiel en teta
C 
      TETTOT = NFONCT*FONTET
C 
C     Longueur totale du tableau des fonctions
C 
      LONTOT = NFONCT*NPICET
C 
      CALL BALAI (SACTEP, LONTOT, 1)
C 
C     SIGMAC    contrainte adm actuelle         X 6*NPICET
C     ADTET     tableau en teta                 X TETTOT
C 
      LONLIF   = 6*NPICET
C 
      CALL POUSMD (LONLIF+NFONCT*NEPS, SIGMAC)
      BMTGAC   = SIGMAC
      EPSLOC   = BMTGAC+LONLIF
C 
      CALL GSPOUD (TETTOT, ADTET)
      AVTET    = ADTET-1
C 
      EPSI   = 1
      LONEPS = NEPS*NTETA*NGAU1
C 
C     Pour recuperer les poids aux points de gauss
C 
      CALL ADTBDM ('POIDS-PGCR', ADPOPG)
      DBPOPG   = ADPOPG
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES COUCHES
C 
      DO NUCOU= 1, NBCOU
C 
        CALL VALEP (NUCOU, B)
C 
        IF (LTYP12) THEN
          CALL ANGCOU (NUCOU, TETORT)
          CALL COMISO (ADCOU, NUCOU, COUISO)
        END IF
C 
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1, NBCOL
C 
         CALL VALRAY (NUCOL, RAYONC, A)
         JLOC    = A*B
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT EN R
C 
          DO X = 1, XINTEG
C 
            RLOC   = RAYONC + A*GAUSS( HINT + X )
            MULTI  = POIDS(HINT+X)*JLOC*RLOC
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT EN Z
C 
            DO Y = 1, YINTEG
C 
C             MULTJ  = POIDS(KINT+Y)*MULTI
C 
              MULTJ = DM(DBPOPG)
              DBPOPG = DBPOPG+1
C 
C             BOUCLE iv SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
                IF (LTYP12) THEN
                  TETCAL = DM(ADTETA+TETA-1)
                  CALL RICLOC (COUISO, TETORT, TETCAL, KLOC)
                END IF
C 
C               Lecture directe du fichier (dans q-chapeau) des accroissements
C               admissibles de l'etape locale precedente.
C 
                CALL LFDDNF (DM(BMTGAC), BMTGAC, LONSCH, IUNSIG, NUENSI)
                NUENSI = NUENSI +1
C 
C               BOUCLE v SUR LES PREMIERES FONCTIONS DU TEMPS
C 
                DO I = 1, NFONCT
C 
C                 EPSILO est l'adresse dans ALPHAI de la deformation
C                 du point de Gauss actuel pour la ieme fonction.
C 
                  EPSILO = EPSI+(I-1)*LONEPS
C 
C                 Deformation de la fonction suivante au meme pt de gauss
C 
                  DEBUT = AVTET+(I-1)*FONTET+TETA
                  SIGMAC     = BMTGAC
C 
C                 BOUCLE vi SUR LES PIQUETS DE TEMPS
C 
                  DO  TEMPS = 1, NPICDE-1
C 
                     DEBUT      = DEBUT+NTETA
                     SIGMAC     = SIGMAC+NEPS
C 
                  END DO
C 
                  DO  TEMPS = NPICDE, NPICFI
C 
                    IF(LTYP0)
     &              CALL SCAVEC (NEPS, DM(SIGMAC),
     &                           ALPHAI(EPSILO), TRACE)
                    IF (LTYP12)
     &              CALL TRAC12 (KLOC, DM(SIGMAC),
     &                           ALPHAI(EPSILO), TRACE)
C 
                    DM(DEBUT) = DM(DEBUT)+MULTJ*TRACE
                    DEBUT     = DEBUT+NTETA
C 
                    SIGMAC    = SIGMAC+NEPS
C 
                  END DO
C 
                  DO  TEMPS = NPICFI+1, NPICET
C 
                    DEBUT      = DEBUT+NTETA
                    SIGMAC     = SIGMAC+NEPS
C 
C                 FIN DE BOUCLE vi SUR LES PIQUETS DE TEMPS
C 
                  END DO
C 
C               FIN DE BOUCLE v SUR LES DEFORMATIONS
C 
                END DO
C 
                EPSI = EPSI+NEPS
C 
C            FIN DE BOUCLE iv SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT EN Z
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT EN R
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES COLONNES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES COUCHES
C 
      END DO
C 
C     Fermeture de l'unite ou est stocke SIGMA
C 
      CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     Sequence d'integration en teta et calcul de la contribution au point de Gauss
C 
C     CALL ITBTET (LONTOT, DM(ADTET), SACTEP)
C 
C     CD    CALL IMPTDT ('Contribution des couches '//IDPROG,
C     CD                 SACTEP, NPICET, NFONCT)
C  
C **********************************************************************
C *
C *     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C *
C **********************************************************************
C 
C      DEBUT DU TEST SUR LES INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
        NOMFIC = SIDSAU
C 
C       NOMFIC est le nom de fichier dans lequel sont stockees les
C       quantites barre de l'etape locale precedente
C       ouverture du fichier pour les contraintes
C 
        LONSCH = 3*NPICET
        NUENSI  = 1
        CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C 
C       Debut  des valeurs du saut au point de Gauss
C 
        DBCPG = 1
C 
C       SAUTGAC   -sign admissible  actuel      X 3*NPICET
C 
        LONLIF   = NSAU*NPICET
C 
        CALL POUSMD (LONLIF+NFONCT*NSAU, SIGNAC)
        BMTGAC   = SIGNAC
        EPSLOC   = BMTGAC+LONLIF
C 
C       ADTET     tableau en teta              X TETTOT
C 
C       ADTEI = ADTET
C       AVTEI = AVTET
C 
C       CD    CALL GSPOUD (TETTOT, ADTEI)
C       CD    AVTEI    = ADTEI-1
C 
        SAUI   = 1
        LONSAU = NSAU*NTETA*NGAU2
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT = 1, NBINT
          IF (LTYP12) THEN
            CALL IOMISO (ADINT, NUINT, INTISO)
            CALL ANGINT (NUINT, TETORT)
          END IF
C 
C         BOUCLE i SUR LES COLONNES
C 
          DO NUCOL = 1, NBCOL
C 
C            CALL VALRAY (NUCOL, RAYONC, A)
C 
C           BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
C 
            DO X = 1, XINTEG
C 
C             RLOC  = RAYONC + A*GAUSS(HINT + X)
C             MULTI  = POIDS(HINT+X)*A*RLOC
C 
              MULTI = DM(DBPOPG)
              DBPOPG = DBPOPG+1
C 
C             BOUCLE iii SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
                IF (LTYP12) THEN
                  TETCAL = DM(ADTETA+TETA-1)
                  CALL INTELA (INTISO, TETORT, TETCAL, KLOCI)
                END IF
C 
C               Lecture directe du fichier (dans q-chapeau) des accroissements
C               admissibles de l'etape locale precedente.
C 
                CALL LFDDNF (DM(BMTGAC), BMTGAC, LONSCH,
     &                       IUNSIG, NUENSI)
C 
                NUENSI = NUENSI +1
C 
C               BOUCLE iv SUR LES PREMIERES FONCTIONS DU TEMPS
C 
                DO I = 1, NFONCT
C 
C                 Definition des debuts de tableaux des
C                 deformations de la 1ere fonction.
C 
                  SAUTLO = SAUI+(I-1)*LONSAU
C 
C                 Deformation de la fonction suivante au meme pt de gauss
C 
                  DEBUT = AVTET+(I-1)*FONTET+TETA
C 
C                 Definition des debuts de tableaux des
C                 sauts de la 1ere fonction
C 
                  SIGNAC = BMTGAC
C 
C                 BOUCLE v SUR LES PIQUETS DE TEMPS
C 
                  DO  TEMPS = 1, NPICDE-1
C 
                    DEBUT  = DEBUT+NTETA
C 
                    SIGNAC   = SIGNAC+3
C 
                  END DO
C 
                  DO  TEMPS = NPICDE, NPICFI
C 
                    IF (LTYP0)
     &              CALL SCAVEC (NSAU, DM(SIGNAC), SAUTI(SAUTLO),
     &                           TRACE)
C 
                    IF (LTYP12)
     &              CALL PROI12 (KLOCI, DM(SIGNAC), SAUTI(SAUTLO),
     &                           TRACE)
C 
                    DM(DEBUT) = DM(DEBUT)+ MULTI*TRACE
C 
                    DEBUT = DEBUT+NTETA
C 
                    SIGNAC = SIGNAC+3
C 
                  END DO
C 
                  DO  TEMPS = NPICFI+1, NPICET
C 
                    DEBUT = DEBUT+NTETA
C 
                    SIGNAC = SIGNAC+3
C 
C                 FIN DE BOUCLE v SUR LES PIQUETS DE TEMPS
C 
                  END DO
C 
C               FIN DE BOUCLE iv SUR LES FONCTIONS DE TEMPS
C 
                END DO
C 
                SAUI = SAUI+NSAU
C 
C             FIN DE BOUCLE iii SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES COLONNES
C 
          END DO
C 
C       FIN DE BOUCLE SUR LES INTERFACES
C 
        END DO
C 
C       Fermeture de l'unite ou est stocke SIGMA
C 
        CALL FERFIC (3, IUNSIG, IDPROG)
C 
C       sequence d'integration en teta et calcul de
C       la contribution au point de Gauss.
C 
C       CD    CALL GSPOUD (LONTOT, SACTEI)
C       CD    CALL IMPTDT ('Contribution de l''interface '//IDPROG,
C       CD                  DM(SACTEI), NPICET, NFONCT)
C       CD    CALL ADDMAD (NPICET*NFONCT, SACTEP, DM(SACTEI), SACTEP)
C 
C     FIN DE TEST SUR LES INTERFACES
C 
      END IF
C 
      CALL ITBTET (LONTOT, DM(ADTET), SACTEP)
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1)) THEN
C 
CD      CALL OMPTDN ('FONCTION DU TEMPS ',
CD                    SACTEP, NPICET, NFONCT)
C 
CD    END IF
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
