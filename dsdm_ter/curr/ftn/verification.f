C     Cette routine sert a verifier que les deplacements des corrections admissisbles 
C     a zero sont nuls la ou on impose des deplacements.
C 
C     (verification de l'admissibilite a zero en deplacement)
C 
C      Si le deplacement est superieur a 10E-10 on le signale
C  
C     On envoie comme arguments :
C                                 
C     ES ...... DEPDEV valeur des efforts en entree ranges (nddl, nbmat)
C 
      SUBROUTINE VDPVZE (DEPDEV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION DEPDEV(NDDL*NBMAT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   DDLBLO, NBTDDL, MAT, I, DEBUT  
      INTEGER   ADDBLO                    
C 
      DOUBLE PRECISION  SUPTAB, COMPAR, VAL
C 
      LOGICAL LOGVAL
CD    LOGICAL LTRACP, LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VDPVZE')
C 
CD     CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBM ('DDL-BLOQUE', DDLBLO)            
      CALL LONGEN ('DDL-BLOQUE', NBTDDL)
C 
      DDLBLO= DDLBLO-1
      DEBUT = 0
      LOGVAL = .FALSE.
C 
      DO MAT = 1, NBMAT
C 
        CALL NORSUP (NDDL, DEPDEV(DEBUT+1), SUPTAB)
        COMPAR = 1.D-10*SUPTAB
C 
        DO I =  1, NBTDDL
C 
          VAL = DABS (DEPDEV(DEBUT + M(DDLBLO+I)))
C 
          IF (VAL .GT. COMPAR) THEN
C 
            LOGVAL = .TRUE.
	    CALL IMPET ('POUR LE DEVELOPPEMENT NUMERO      '//IDPROG,
     &                   MAT-NTDSFG+1  ) 
            CALL IMPET ('POUR LE DDL A DEP IMPOSE  NUMERO  ',
     &                   M(DDLBLO+I))
            CALL IMPDT ('10-10*VALEUR > SUP DU DEPLACEMENT ', COMPAR)
            CALL IMPDT ('VALEUR DU DEPLACEMENT IMPOSE =0?  ', VAL)
C 
          END IF
C 
        END DO
C 
        DEBUT = DEBUT + NDDL 
C 
      END DO
C 
      IF (.NOT. LOGVAL) THEN
        CALL MESSAO ('VERIFICATION DEPLACEMENTS ADMISSIBLES OK')
      ENDIF
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END                         
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Verification des orthogonalisations et des modifications des champs
C     accroissements en deformations ou contraintes.
C  
C     TYPE = 0, 1, 2 ou 3
C  
C     SI TYPE = 0 alors on verifie les valeurs totales en deformation
C     SI TYPE = 1 alors on verifie les accroissements en deformation
C     SI TYPE = 2 alors on verifie les valeurs totales en contrainte
C     SI TYPE = 3 alors on verifie les accroissements en contrainte
C  
C     On cree dans q-chapeau un fichier verif(0-3) detruit par FIORMO.
C 
      SUBROUTINE DEORMO (TYPE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'    
      include 'cominc_visu.h'    
C 
      INTEGER TYPE
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
C     LOGIQUE SUIVANT TYPE
C 
      LOGICAL LTYP0, LTYP1, LTYP2, LTYP3, LTYP02, LTYP13
C 
      CHARACTER*1 IDETYP
C 
      INTEGER NBFTMX, DBFT, FIFT, DBCHA, FICHA, FONTEM 
      INTEGER CHAESP, CHAES2, LONPRO, CHAAPG
      INTEGER LONCHA, NUENCA, IUADCH
      INTEGER ADCHA, DQADMP, CHAADM 
      INTEGER NUCOL, NUCOU, X, Y, TETA, TEMPS
      INTEGER NUINT, DEBGAU, FINGAU, PGAU1
C 
      CHARACTER*6 NOM, NOMFIC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DEORMO')
C 
      INTEGER    AM2LC, ADM2LC     
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)     
C 
C -----------------------------------------------------------------------
C 
C     Logique de type
C 
      LTYP0  = .FALSE.
      LTYP1  = .FALSE.
      LTYP2  = .FALSE.
      LTYP3  = .FALSE.
      LTYP02 = .FALSE.
      LTYP13 = .FALSE.
C 
      IF (TYPE .EQ. 0) LTYP0 = .TRUE.
      IF (TYPE .EQ. 1) LTYP1 = .TRUE.
      IF (TYPE .EQ. 2) LTYP2 = .TRUE.
      IF (TYPE .EQ. 3) LTYP3 = .TRUE.
      IF ((TYPE .EQ. 0) .OR .(TYPE .EQ. 2)) LTYP02 = .TRUE.
      IF ((TYPE .EQ. 1) .OR .(TYPE .EQ. 3)) LTYP13 = .TRUE.
C 
C 
C     RECHERCHE DES ADRESSES DES TABLEAUX :
C 
C     contraintes admissibles (nsig, nteta, ngau1, nfotps)
C     deformations admissibles (neps, nteta, ngau1, nfotps)
C     accroissements des fonctions du temps donnees (npicet, nbfodo)
C     initialisation des quantites locales pour le temps zero
C    
C     NDEBCA est le numero du premier champ servant au calcul  
C   
      CALL IDENT1 (TYPE, IDETYP)
C 
      IF (LTYP0) THEN
C 
        NBFTMX =  CHARAX
        CALL ADTBDM ('TEMPS-REEL', FONTEM)
C 
        DBFT  = 1
        FIFT  = NBDPTR
C 
        DBCHA =  1
        FICHA =  DEADTR
C 
D       CALL IMPET ('RECALCUL : TYPE DANS '//IDPROG, TYPE)
C 
        CALL ADTBDM ('EPS-AD-TOT', CHAESP)
C 
C       POUR LES INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
           CALL ADTBDM ('SAU-AD-TOT', CHAES2)
C 
        ENDIF
C 
      ELSE IF (LTYP1) THEN
C 
        NBFTMX = CHAMAX
        DBFT   = 1
        FIFT   = NBFEPS
C 
        DBCHA  =  DEADTR-NBFEPS+1
        FICHA  =  DEADTR
C 
        CALL ADTBDM ('TEMPS-EPSI', FONTEM)
        CALL ADTBDM ('EPS-AD-TOT', CHAESP)
C 
C       POUR LES INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
           CALL ADTBDM ('SAU-AD-TOT', CHAES2)
C 
        ENDIF
C 
D       CALL IMPET ('RECALCUL : TYPE DANS '//IDPROG, TYPE)
C 
      ELSE IF (LTYP2) THEN
C 
        NBFTMX =  CHARAX
        CALL ADTBDM ('TEMP-SI-RE', FONTEM)
C 
        DBFT  = 1
        FIFT  = EVCOTR
C 
        DBCHA =  1
        FICHA =  COTORE
C 
        CALL ADTBDM ('SIG-AD-TOT', CHAESP)
C 
C       POUR LES INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL ADTBDM ('SGN-AD-TOT', CHAES2)
C 
        ENDIF
C 
D       CALL IMPET ('RECALCUL : TYPE DANS '//IDPROG, TYPE)
C 
        ELSE IF (LTYP3) THEN
C 
        NBFTMX =  CHAMAX
        DBFT   = 1
        FIFT   = NBFSIG
C 
        DBCHA  =  COTORE-NBFSIG+1
        FICHA  =  COTORE
C 
        CALL ADTBDM ('TEMPS-SIGM', FONTEM) 
        CALL ADTBDM ('SIG-AD-TOT', CHAESP)
C 
C       POUR LES INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL ADTBDM ('SGN-AD-TOT', CHAES2)
C 
        ENDIF
C 
D       CALL IMPET ('RECALCUL : TYPE DANS '//IDPROG, TYPE)
C 
      END IF 
C 
      LONPRO = 6*NPICET
C 
      CALL POUSMD (LONPRO, CHAAPG)
C 
C     Choix du nom des fichier dans Q-chapeau : (3)
C 
      NOM    = 'vech_'//IDETYP
      NOMFIC = NOM                          
      LONCHA = 6*NPICET
      NUENCA = 1
      CALL OFDDNF (3, NOMFIC, 6, LONCHA, IUADCH)
C 
      ADCHA  = CHAESP
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
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL  = 1, NBCOL   
C 
C       BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT X
C 
          DO X  = 1, XINTEG
C 
C           BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Y
C 
	    DO  Y =  1, YINTEG
C 
C             BOUCLE iii SUR LES ANGLES
C  
              DO TETA = 1, NTETA 
C 
C               ATTENTION TOUTE LES QUANTITES SONT MULTIPLIES PAR INTERV
C               on calcule la valeur des deltas des champs admissibles pour tout les 
C               piquets de temps : debut caracterise la position dans le tableau des 
C               contraintes ou des deformations admissibles de la 1ere valeur interessante.
C 
                CALL VCPGTE (
C               Valeur des Champs aux Points de Gauss pour tous les TEmps
C               on envoie
     &          NBFTMX, NEPS, DBFT, FIFT, DBCHA, FICHA, DM(ADCHA),
     &          DM(FONTEM),
C               on recupere
     &          DM(CHAAPG))
C 
C               Pour aller lire le bon point de Gauss dans les tableaux de stockage 
C 
                ADCHA = ADCHA + CHARAX*NEPS
                CALL EFDDNF (IUADCH, DM(CHAAPG), CHAAPG, LONCHA, NUENCA)
                NUENCA = NUENCA + 1
C 
C             FIN DE BOUCLE iii SUR LES ANGLES
C  
              END DO
C 
C           FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Y
C 
            END DO                                  
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT X
C 
          END DO
C 
C         FIN DE BOUCLE i SUR LES COLONNES
C 
        END DO
C 
C     BOUCLE SUR LES COUCHES
C 
      END DO
C 
      CALL FERFIC (3, IUADCH, IDPROG)
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
C       Ouverture du fichier pour les contraintes normales.
C 
        LONCHA = 3*NPICET
        NOM    = 'vsch_'//IDETYP
        NOMFIC = NOM                          
        LONCHA = 3*NPICET
        NUENCA = 1
        CALL OFDDNF (3, NOMFIC, 6, LONCHA, IUADCH)
        ADCHA  = CHAES2
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT = 1 , NBINT
C 
C         Recherche des caracteristiques du comportement de l'interface.
C 
          DEBGAU   =  1+(NUINT-1)*XINTEG*NBCOL
          FINGAU   =  NUINT*XINTEG*NBCOL
C 
C         BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
          DO PGAU1 = DEBGAU, FINGAU
C 
C           BOUCLE iii SUR LES ANGLES
C 
            DO TETA = 1 , NTETA
C 
              CALL VCPGTE (
C             Valeur des Champs au Point de Gauss pour tout les TEmps
C             on envoie
     &        NBFTMX, NSAU, DBFT, FIFT, DBCHA, FICHA,
     &                DM(ADCHA), DM(FONTEM),
C             on recupere
     &        DM(CHAAPG))
C 
C             Pour aller lire le bon point de Gauss dans les tableaux de stockage 
C 
              ADCHA = ADCHA + CHARAX*NSAU
C 
              CALL EFDDNF (IUADCH, DM(CHAAPG), CHAAPG, LONCHA, NUENCA)
              NUENCA = NUENCA +1
C 
C           FIN DE BOUCLE iii SUR LES ANGLES
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
        CALL FERFIC (3, IUADCH, IDPROG)
C 
C     FIN DE TEST SUR LES INTERFACES
C 
      END IF 
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
C     Verification des orthogonalisations et des modifications 
C     des fonctions du temps 
C  
C     TYPE = 0, 1, 2 ou 3
C  
C     SI TYPE = 0 alors on verifie les valeurs totales en deformation
C     SI TYPE = 1 alors on verifie les accroissements  en deformation
C     SI TYPE = 2 alors on verifie les valeurs totales en contrainte
C     SI TYPE = 3 alors on verifie les accroissements  en contrainte
C  
C     On cree dans q chapeau un fichier verif(0-3) detruit par FIORMO
C 
      SUBROUTINE FIORMO (TYPE)
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'    
      include 'cominc_visu.h'    
C 
      INTEGER TYPE
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
C     LOGIQUE SUIVANT TYPE
C 
      LOGICAL LTYP0, LTYP1, LTYP2, LTYP3, LTYP02, LTYP13, ZOBI
C 
      CHARACTER*1 IDETYP
C 
      INTEGER NBFTMX, DBFT, FIFT, DBCHA, FICHA, FONTEM 
      INTEGER CHAESP, CHAES2, LONPRO, CHAAPG
      INTEGER LONCHA, NUENCA, IUADCH
      INTEGER ADCHA, DQADMP, CHAADM, DBCHAD
      INTEGER NUCOL, NUCOU, X, Y, TETA, TEMPS
      INTEGER NUINT, DEBGAU, FINGAU, PGAU1
      INTEGER I, DBADM, DBPAG
C 
      DOUBLE PRECISION  CHADIF(6), TEST
C 
      CHARACTER*6 NOM, NOMFIC
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='FIORMO')
C 
      INTEGER    AM2LC, ADM2LC     
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)     
C 
C -----------------------------------------------------------------------
CD    CALL IMPET ('LONGUEUR RESERVEE DANS DM '//IDPROG, LDM)
C 
C     Logique de type
C 
      LTYP0  = .FALSE.
      LTYP1  = .FALSE.
      LTYP2  = .FALSE.
      LTYP3  = .FALSE.
      LTYP02 = .FALSE.
      LTYP13 = .FALSE.
C 
      IF (TYPE .EQ. 0) LTYP0 = .TRUE.
      IF (TYPE .EQ. 1) LTYP1 = .TRUE.
      IF (TYPE .EQ. 2) LTYP2 = .TRUE.
      IF (TYPE .EQ. 3) LTYP3 = .TRUE.
      IF ((TYPE .EQ. 0) .OR .(TYPE .EQ. 2)) LTYP02 = .TRUE.
      IF ((TYPE .EQ. 1) .OR .(TYPE .EQ. 3)) LTYP13 = .TRUE.
C 
C     RECHERCHE DES ADRESSES DE TABLEAUX :
C 
C     contraintes  admissibles (nsig, nteta, ngau1, nfotps)
C     deformations admissibles (neps, nteta, ngau1, nfotps)
C     accroissements des fonctions du temps donnees (npicet, nbfodo)
C     initialisation  des quantites locales pour le temps zero
C 
C     NDEBCA est le numero du premier champ servant au calcul  
C 
C     Pour le calcul du champ la premiere adresse a servir 
C     dans le tableau est ADCHA le decalage de NDEBCA se faisant dans ADCTPS
C 
      CALL IDENT1 (TYPE, IDETYP)
C 
      IF (LTYP0) THEN
C 
        NBFTMX =  CHARAX
        CALL ADTBDM ('TEMPS-REEL', FONTEM)
C 
        DBFT  = 1
        FIFT  = NBDPTR
C 
        DBCHA =  1
        FICHA =  DEADTR
C 
D       CALL IMPET ('RECALCUL : TYPE DANS '//IDPROG, TYPE)
C 
        CALL ADTBDM ('EPS-AD-TOT', CHAESP)
C 
C       POUR LES INTERFACES
C 
        IF (NBINT. GT. 0) THEN
C 
           CALL ADTBDM ('SAU-AD-TOT', CHAES2)
C 
        ENDIF
C 
      ELSE IF (LTYP1) THEN
C 
        NBFTMX =  CHAMAX
        DBFT   = 1
        FIFT   = NBFEPS
C 
        DBCHA  =  DEADTR-NBFEPS+1
        FICHA  =  DEADTR
C 
        CALL ADTBDM ('TEMPS-EPSI', FONTEM)
        CALL ADTBDM ('EPS-AD-TOT', CHAESP)
C 
C       POUR LES INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL ADTBDM ('SAU-AD-TOT', CHAES2)
C 
        ENDIF
C 
D       CALL IMPET ('RECALCUL : TYPE DANS '//IDPROG, TYPE)
C 
      ELSE IF (LTYP2) THEN
C 
        NBFTMX =  CHARAX
        CALL ADTBDM ('TEMP-SI-RE', FONTEM)
C 
        DBFT  = 1
        FIFT  = EVCOTR
C 
        DBCHA =  1
        FICHA =  COTORE
C 
        CALL ADTBDM ('SIG-AD-TOT', CHAESP)
C 
C       POUR LES INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL ADTBDM ('SGN-AD-TOT', CHAES2)
C 
        ENDIF
C 
D       CALL IMPET ('RECALCUL : TYPE DANS '//IDPROG, TYPE)
C 
        ELSE IF (LTYP3) THEN
C 
        NBFTMX =  CHAMAX
        DBFT   = 1
        FIFT   = NBFSIG
C 
        DBCHA  =  COTORE-NBFSIG+1
        FICHA  =  COTORE
C 
        CALL ADTBDM ('TEMPS-SIGM', FONTEM) 
        CALL ADTBDM ('SIG-AD-TOT', CHAESP)
C 
C       POUR LES INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL ADTBDM ('SGN-AD-TOT', CHAES2)
C 
        END IF
C 
D       CALL IMPET ('RECALCUL : TYPE DANS '//IDPROG, TYPE)
C 
      END IF 
C 
      LONPRO = 12*NPICET
C 
      CALL POUSMD (LONPRO, CHAAPG)
C 
      DBCHAD = CHAAPG+6*NPICET
C 
C     Choix du nom des fichier dans Q-chapeau :(3)
C 
      NOM    = 'vech_'//IDETYP
      NOMFIC = NOM                          
      LONCHA = 6*NPICET
      NUENCA = 1
      CALL OFDDNF (3, NOMFIC, 6, LONCHA, IUADCH)
C 
      ADCHA  = CHAESP
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
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1, NBCOL   
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT X
C 
          DO X = 1, XINTEG
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Y
C 
            DO  Y = 1, YINTEG
C 
C             BOUCLE iv SUR LES ANGLES
C 
              DO TETA = 1, NTETA 
C 
C               ATTENTION TOUTE LES QUANTITES SONT MULTIPLIES PAR INTERV
C               on calcule la valeur des deltas des champs admissibles pour tout les 
C               piquets de temps : debut caracterise la position dans le tableau des 
C               contraintes ou des deformations admissibles de la 1ere valeur interessante.
C                                                                
                CALL VCPGTE (
C               Valeur des Champs au Point de Gauss pour tous les TEmps
C               on envoie
     &          NBFTMX, NEPS, DBFT, FIFT, DBCHA, FICHA,
     &          DM(ADCHA), DM(FONTEM),
C               on recupere
     &          DM(CHAAPG))
C 
C               Pour aller lire le bon point de Gauss dans les tableaux de stockage 
C 
                ADCHA = ADCHA + CHARAX*NEPS
                CALL LFDDNF (DM(DBCHAD), DBCHAD, LONCHA, IUADCH, NUENCA)
                NUENCA = NUENCA +1
                DBADM = DBCHAD
                DBPAG = CHAAPG
C 
C               BOUCLE v SUR LES PIQUETS DE TEMPS
C 
                ZOBI = .FALSE.
		DO TEMPS = 1, NPICET
                  CALL SOUMAD (6, DM(DBPAG), DM(DBADM), CHADIF)
                  DO I = 1, 6 
                    TEST = DABS(CHADIF(I))
                    IF (TEST .GT. 1.D-5*DABS(DM(DBPAG+I-1)) 
     &                 .AND. TEST .GT. 1.D-15) THEN
C 
                     ZOBI = .TRUE.
                    END IF
                  END DO
C 
C               FIN DE BOUCLE v SUR LES PIQUETS DE TEMPS
C 
                END DO
C 
C             FIN DE BOUCLE iv SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Y
C 
            END DO                                  
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT X
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
      CALL FERFIC (3, IUADCH, IDPROG)          
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
C       Ouverture du fichier pour les contraintes normales.
C 
        LONCHA = 3*NPICET
C 
C       Choix du nom des fichier dans Q-chapeau : (3)
C 
        NOM    = 'vsch_'//IDETYP
        NOMFIC = NOM                          
        LONCHA = 3*NPICET
        NUENCA = 1
        CALL OFDDNF (3, NOMFIC, 6, LONCHA, IUADCH)
        ADCHA  = CHAES2
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT = 1, NBINT
C 
          DEBGAU =  1+(NUINT-1)*XINTEG*NBCOL
          FINGAU =  NUINT*XINTEG*NBCOL
C 
C         BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
          DO PGAU1  = DEBGAU , FINGAU
C 
C           BOUCLE iii SUR LES ANGLES
C 
            DO TETA = 1, NTETA
C 
              CALL VCPGTE (
C             Valeur des Champs au Point de Gauss pour tous les TEmps
C             on envoie
     &        NBFTMX, NSAU, DBFT, FIFT, DBCHA, FICHA,
     &        DM(ADCHA), DM(FONTEM),
C             on recupere
     &        DM(CHAAPG))
C 
C             Pour aller lire le bon point de Gauss dans les tableaux de stockage 
C 
              ADCHA = ADCHA + CHARAX*NSAU
              CALL LFDDNF (DM(DBCHAD), DBCHAD, LONCHA, IUADCH, NUENCA)
              NUENCA = NUENCA + 1                        
C 
              DBPAG = CHAAPG
              DBADM = DBCHAD
C 
C             BOUCLE iv SUR LES PIQUETS DE TEMPS
C 
              DO TEMPS = 1, NPICET
C 
                CALL SOUMAD (3, DM(DBPAG), DM(DBADM), CHADIF)
                DO I = 1 , 3 
                  TEST = DABS(CHADIF(I))
                  IF (TEST .GT. 1.D-5*DABS(DM(DBPAG+I-1)) 
     &               .AND. TEST .GT. 1.D-10) THEN
                   ZOBI = .TRUE.
                  END IF
                END DO
C 
C             FIN DE BOUCLE iv SUR LES PIQUETS DE TEMPS
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES ANGLES
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
        CALL FERFIC (3, IUADCH, IDPROG)
C 
C       FIN DE TEST SUR LES INTERFACES
C 
        IF (.NOT. ZOBI) THEN
          CALL MESSAO ('VERIFICATION DEPLACEMENTS TOTAUX OK')
        ELSE
          CALL MESSAO ('PROBLEME DEPLACEMENTS TOTAUX')
	  CALL IMPET ('APPEL DEORMO-FIORMO TYPE ', TYPE)
        END IF
C 
      END IF 
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
