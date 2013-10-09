C     Cette routine remplit les tableaux des points et des rayons
C     aux points de gauss dans l'ordre adapte au produit scalaire en espace
C 
C     On envoie comme arguments:
C 
      SUBROUTINE REMRPG
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
      DOUBLE PRECISION  A, B, JLOC
      DOUBLE PRECISION  RLOC, RAYONC, MULTI
C 
      INTEGER     NUCOL, NUCOU, NUINT
      INTEGER     AM2LC, ADM2LC, HINT, KINT
      INTEGER     X, Y, ADPOPG, ADRAPG, LONPOP, LONRAP
      INTEGER     DBPOPG, DBRAPG
C 
CD    LOGICAL     LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='REMRPG')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
      LONPOP = NBCOL*XINTEG*(NBCOU*YINTEG+NBINT)
      LONRAP = NBCOL*XINTEG*(NBCOU+NBINT)
C 
      CALL GESTDP('POIDS-PGCR', LONPOP , ADPOPG )
      CALL GESTDP('RAYON-PGCR', LONRAP , ADRAPG )
C 
C     Pour aller lire les points de Gauss
C 
      HINT = XINTEG*(XINTEG-1)/2
      KINT = YINTEG*(YINTEG-1)/2
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
      DBPOPG = ADPOPG
      DBRAPG = ADRAPG
C 
C     BOUCLE SUR LES COUCHES
C 
      DO NUCOU = 1 , NBCOU
C 
        CALL VALEP( NUCOU , B )
C 
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1 , NBCOL
C 
          CALL VALRAY( NUCOL , RAYONC , A )
          JLOC    = A*B
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS EN r
C 
          DO X = 1 ,XINTEG
C 
            RLOC         = RAYONC + A*GAUSS( HINT + X )
            DM(DBRAPG)   = RLOC
               DBRAPG    = DBRAPG+1
            MULTI        = POIDS(HINT+X)*JLOC*RLOC
C 
C          BOUCLE iii SUR LES POINTS DE GAUSS EN z
C 
            DO Y = 1 ,YINTEG
C 
              DM(DBPOPG)  = POIDS(KINT+Y)*MULTI
CD 	      CALL IMPDT ('VERIF 1 MULPPG ', DM(DBPOPG))
              DBPOPG   = DBPOPG+1
C 
C          FIN DE BOUCLE iii SUR LES POINTS DE GAUSS EN z
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS EN r
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
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES INTERFACES
C 
      DO NUINT = 1 , NBINT
C 
C        BOUCLE i SUR LES COLONNES
C 
         DO NUCOL = 1 , NBCOL
C 
               CALL VALRAY( NUCOL , RAYONC , A )
C 
C             BOUCLE ii SUR LES POINTS DE GAUSS EN r
C 
               DO X = 1 ,XINTEG
C 
               RLOC         = RAYONC + A*GAUSS( HINT + X )
               DM(DBRAPG)   = RLOC
               DBRAPG       = DBRAPG+1
C 
               DM(DBPOPG)  = POIDS(HINT+X)*A*RLOC
               DBPOPG   = DBPOPG+1
C 
C             FIN DE BOUCLE ii SUR LES POINTS DE GAUSS r
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES COLONNES
C 
           END DO
C 
C     FIN DE BOUCLE SUR LES INTERFACES
C 
      END DO
C 
C -----------------------------------------------------------------------
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C -----------------------------------------------------------------------
C 
CD    IF (LTRACN(1))THEN
C 
CD    CALL IMPTDT ('POIDS RANGES CROISSANT ', DM(ADPOPG), LONPOP, 1)
C 
CD      CALL OMPTDN( ' RAYONS RANGES CROISSANT' ,
CD                     DM(ADRAPG)  , LONRAP , 1 )

CD    END IF
C 
      CALL SOPOUB(AM2LC,ADM2LC)
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Calcul de la solution orthotrope associee a la
C     solution isotrope transverse en deplacement
C 
C     Cette solution sera prise comme solution admissible initiale
C     (dans le cas de plusieurs histoires il faudrait faire plusieurs
C     calculs)
C 
C     On envoie comme arguments :
C 
C     ES...... DEISOD  Valeur de la solution isotrope
C                      transverse developpee.
C 
C     Et on recupere :
C 
C     ES...... DEORTD  Valeur des deplacements orthotropes solution
C                      developpes
C     S ...... EPSORT  Valeur des deformations orthotrope solution
C     S ...... SIGORT  Valeur des contraintes orthotropes solution
C     S ...... SAUORT  Valeur des sauts orthotropes solution
C     S ...... SGNORT  Valeur des contraintes normales orthotropes
C                      solution

      SUBROUTINE SOLORT (DEISOD, DEORTD, EPSORT, SIGORT, SAUORT, SGNORT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C     ENTREE
C 
      DOUBLE PRECISION DEISOD(NBMAT*NDDL)
C 
C     SORTIE
C 
      DOUBLE PRECISION DEORTD(NBMAT*NDDL)
      DOUBLE PRECISION EPSORT(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION SIGORT(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION SAUORT(NSAU*NTETA*NGAU2)
      DOUBLE PRECISION SGNORT(NSAU*NTETA*NGAU2)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
C     POUR LES COUCHES :
C 
      INTEGER          ADISOS, ADTETA, NUCOU, ADRORT
      INTEGER          DEBGAU, FINGAU, PGAU1, TETA
      INTEGER          EPSISO, SIGISO, DEPCHA, DSICHA, DBSICH
C 
C     POUR LES INTERFACES : 
C 
      INTEGER          NUINT, SAUISO, SINISO, DSACHA, DSNCHA
      INTEGER          DBSNCH, DBSCHA
      INTEGER          NBDEV, TNUDEV
C 
C     POUR LES SOLUTIONS ADMISSIBLES A ZERO :
C 
      INTEGER          ADSMEG, EPSSOL, DEPSOL, SIGRES, SGNRES
      INTEGER          SAUSOL, LONRES
      INTEGER          AM2LC, ADM2LC, CARPOI
      DOUBLE PRECISION CARNLI(12), SOUORT( 17 )
      DOUBLE PRECISION YMAXAC, YMAXPR, EPCORT(6), SACORT(3)
      INTEGER          ADCOU, ADSOU, ADINT, ADSNT
      DOUBLE PRECISION KLOC(17), COUISO(10), SOUISO(10), SLOC(17)
      DOUBLE PRECISION KLOCI(9), INTISO(3), INSISO(3), SLOCI(9)
      DOUBLE PRECISION TETORT, ANGORT, TETCAL
      DOUBLE PRECISION NORME
      LOGICAL          REPRIS
CD    LOGICAL          LTRACP
      CHARACTER*6      IDPROG
      PARAMETER        (IDPROG='SOLORT')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL IMPET ('ADM2 EN ENTREE DANS '//IDPROG, ADM2)
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
      CALL ADTBDM ('HOO-COUCHE', ADCOU)
      CALL ADTBDM ('SOU-COUCHE', ADSOU)
C 
      IF (NBINT.GT. 0) THEN
        CALL ADTBDM ('HOO-INTERF', ADINT)
        CALL ADTBDM ('SOU-INTERF', ADSNT)
      ELSE
          ADINT = ADCOU
          ADSNT = ADSOU
      END IF
C 
C     Calcul des deformations associees a la solution isotrope
C 
      CALL TOUEPS (1, DEISOD, EPSORT)
C 
C     Calcul des contraintes associees a la solution isotrope
C 
      CALL TOUSIG (1, EPSORT, SIGORT)
C 
      IF (NBINT .GT. 0) THEN
C 
C       Calcul des sauts associes a la solution isotrope
C 
        CALL TOUSAU (1, DEISOD, SAUORT)
C 
C       Calcul des contraintes normales associes a la solution isotrope
C 
        CALL TOUSGN (1, SAUORT, SGNORT)
C 
      END IF
C 
      EPSISO= 1
      SIGISO= 1
C 
      CALL POUSMD (NTETA*NGAU1*12, DEPCHA)
C 
      DSICHA = DEPCHA+NTETA*NGAU1*6
      DBSICH = DSICHA
C 
C -----------------------------------------------------------------------
C 
C     DEBUT DU CALCUL AVEC LA SOLUTION ISOTROPE TRANSVERSE
C 
C -----------------------------------------------------------------------
C 
C     Calcul des delta entre isotrope et orthotrope
C -----------------------------------------------------------------------
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES COUCHES
C 
      DO NUCOU = 1 , NBCOU
C 
C     Recherche des caracteristiques du comportement de la couche
C 
        CALL COMISO (ADCOU, NUCOU, COUISO)
        CALL COMISO (ADSOU, NUCOU, SOUISO)
C 
        CALL ANGCOU( NUCOU , TETORT )
C 
C       BOUCLE i SUR LES POINTS DE GAUSS DANS LA COUCHE
C 
        DEBGAU   =  1+(NUCOU-1)*XINTEG*YINTEG*NBCOL
        FINGAU   =  NUCOU*XINTEG*YINTEG*NBCOL
C 
        CALL CNLINC (NUCOU, CARNLI, ADISOS, ADRORT)

        DO PGAU1  = DEBGAU, FINGAU
C 
C         BOUCLE ii SUR LES ANGLES
C 
          DO TETA = 1 , NTETA
C 
C           Recherche de l'angle de la bande (tetcal) correspondant a teta
C 
            TETCAL = DM (ADTETA+TETA-1)
C 
C           Initialisation des quantites de l'iteration precedente
C 
            CALL RICLOC (COUISO, TETORT, TETCAL, KLOC)
C 
C           Souplesse tangente
C 
            CALL RICLOC (SOUISO, TETORT, TETCAL, SLOC)
C 
            ANGORT = TETCAL - TETORT
C 
            CALL CHAPOR
     &            (ANGORT, CARNLI, EPSORT(EPSISO), SIGORT(SIGISO),
     &             SLOC, KLOC,
C 
C           Et on recupere :
C 
C           S ...... DEPCHA   delta epsilon point chapeau
C           S ...... DSICHA   -delta sigma point chapeau
C           S ...... YMAXAC   valeur maximale de Ybcou actuelle
C 
     &              EPCORT, DM(DEPCHA), DM(DSICHA), YMAXAC)
C 
C     <=> on calcule -2delta (epsilon-point-chapeau)
C 
            CALL MUMARE (-2.D0, 6, DM(DEPCHA), DM(DEPCHA))
C  
C     Pour avoir directement les contraintes permettant le calcul
C     des seconds membres pour l'etape globale on multiplie par 2
C     <=> on stocke -2delta( sigma-point-chapeau)
C 
            CALL ADDMAD (6, DM(DSICHA), DM(DSICHA), DM(DSICHA))
C 
            EPSISO = SIGISO+6
            SIGISO = SIGISO+6
            DEPCHA = DEPCHA+6
            DSICHA = DSICHA+6
C 
C         FIN DE BOUCLE ii SUR LES ANGLES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES POINTS DE GAUSS DANS LA COUCHE
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES COUCHES
C 
      END DO
C 
C     SEQUENCE POUR LES INTERFACES
C 
C     Debut de test sur le nombre d'interface
C     Pour le cas ou il n'y a pas d'interface
C 
      DBSNCH = DBSICH
C 
      IF (NBINT .GT. 0) THEN
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
        YMAXAC   = 0.D0
        SINISO   = 1
        SAUISO   = 1
C 
        CALL POUSMD(NTETA*NGAU2*6 , DSACHA )
C 
        DBSCHA = DSACHA
        DSNCHA = DSACHA+NTETA*NGAU2*3
        DBSNCH = DSNCHA
C 
C     BOUCLE SUR LES INTERFACES
C 
      DO NUINT = 1 , NBINT
C 
C     Recherche des caracteristiques du comportement de l'interface
C 
        CALL IOMISO (ADINT, NUINT, INTISO)
        CALL IOMISO (ADSNT, NUINT, INSISO)
        CALL ANGINT (NUINT, TETORT)
C 
        CALL CNLINI (NUINT, CARNLI, SOUORT)
C 
C       BOUCLE i SUR LES POINTS DE GAUSS POUR L'INTERFACE
C  
        DEBGAU   =  1+(NUINT-1)*XINTEG*NBCOL
        FINGAU   =  NUINT*XINTEG*NBCOL
C 
        DO PGAU1  = DEBGAU , FINGAU
C 
C         BOUCLE ii SUR LES ANGLES
C 
          DO TETA = 1 , NTETA
C 
C           Recherche de l' angle de la bande (tetcal) correspondant a teta
C 
            TETCAL = DM (ADTETA+TETA-1)
            ANGORT = TETCAL-TETORT
C 
            CALL INTELA (INTISO, TETORT, TETCAL, KLOCI)
            CALL INTELA (INSISO, TETORT, TETCAL, SLOCI)
C 
C           Matrice elastique comportement tangent
C 
            CALL  IHAPOR
     &            (ANGORT, CARNLI, SAUORT(SAUISO), SGNORT(SINISO),
     &             SLOCI, KLOCI, SACORT, DM(DSACHA), DM(DSNCHA), YMAXAC)
C 
C           Pour avoir directement les sauts permettant le calcul
C           des seconds membres pour l'etape globale on multiplie par 2
C           <=> on stocke -2delta( saut-point-chapeau)
C 
            CALL MUMARE (-2.D0, 3, DM(DSACHA), DM(DSACHA))
C 
C           Pour avoir directement les contraintes permettant le calcul
C           des seconds membres pour l'etape globale on multiplie par 2
C           <=> on calcule -2delta( sigma-point-chapeau)
C 
            CALL ADDMAD (3, DM(DSNCHA), DM(DSNCHA), DM(DSNCHA))
C 
            SINISO   = SINISO+3
            SAUISO   = SAUISO+3
            DSNCHA   = DSNCHA+3
            DSACHA   = DSACHA+3
C 
C         FIN DE BOUCLE ii SUR LES ANGLES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES POINTS DE GAUSS DANS L'INTERFACE
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES INTERFACES
C 
      END DO
C 
C     fin de test sur le nombre d'interface
C 
      END IF
C 
C -----------------------------------------------------------------------
C 
C     FIN DU CALCUL AVEC LA SOLUTION ISOTROPE TRANSVERSE
C 
C -----------------------------------------------------------------------
C 
C     Calcul des efforts globaux associes aux delta
C -----------------------------------------------------------------------
C     GOTO 1001
C 
      CALL GSPOUD( NDDL*NBMAT , ADSMEG )
C 
      CALL POUSME (NBMAT , TNUDEV )
C 
C     CALL NOUEFF (1, DM( DBSICH), DM(DBSNCH))
C 
      CALL EFFORT (DM(DBSICH), DM(DBSNCH), ADSMEG, NBDEV,
     &             M(TNUDEV),NORME)
C 
C -----------------------------------------------------------------------
C 
C     Calcul de la solution admissible a zero orthotrope
C -----------------------------------------------------------------------
      LONRES = NDDL*NBMAT+NTETA*(NGAU1*NEPS+NGAU2*NSAU)
      CALL GSPOUD( LONRES , DEPSOL )
C 
      EPSSOL = DEPSOL + NDDL*NBMAT
      SAUSOL = EPSSOL + NTETA*NGAU1*NEPS
C 
C     a priori on ne se sert pas du residu
C 
      LONRES = NTETA*(NGAU1*NEPS+NGAU2*NSAU)
      CALL POUSMD (LONRES, SIGRES)
      SGNRES = SIGRES + NTETA*NGAU1*NEPS
C 
      REPRIS = .FALSE.
C 
      CALL GCADMI (30, REPRIS, ADSMEG, NBDEV, M(TNUDEV),
     &             DM(DEPSOL), DM(EPSSOL), DM(SAUSOL),
     &             DM(SIGRES), DM(SGNRES))
C 
CD         CALL VDPVZE( DM(DEPSOL) )
C 
C     Calcul de la contrainte admissible
C     a zero au sens du probleme apres l'etape preliminaire.
C 
        CALL CSIAD0 (DM(DEPSOL), DM(EPSSOL), DM(SAUSOL),
     &               DM(DBSICH), DM(DBSNCH))
C 
C -----------------------------------------------------------------------
C 
C     Calcul de la solution admissible orthotrope
C -----------------------------------------------------------------------
      CALL ADD (DEISOD, DM(DEPSOL), DEORTD, NDDL , NBMAT)
      CALL ADD (EPSORT, DM(EPSSOL), EPSORT, NTETA*NGAU1, NEPS)
      CALL ADD (SIGORT, DM(DBSICH), SIGORT, NTETA*NGAU1, NEPS)
C 
      IF (NBINT .GT.0) THEN
C 
        CALL ADD (SAUORT, DM(SAUSOL), SAUORT, NTETA*NGAU2, NSAU)
C 
C       MODIF CALL ADD (SGNORT, DM(SGNZER), SGNORT, NTETA*NGAU2, NSAU)
C 
        CALL ADD (SGNORT, DM(DBSNCH), SGNORT, NTETA*NGAU2, NSAU)
C 
      END IF
C 
C -----------------------------------------------------------------------
C 
C     Determination du point le plus endommage
C -----------------------------------------------------------------------
      YMAXPR = 0.D0
C 
C     Reservation de place pour les caracteristiques du point le
C     plus endommage (au vu d'une prevision elastique), ce
C     qui devrait permettre de trouver une discretisation en temps
C     pas trop debile (couche+interface)
C 
      CALL GESTEN ('PT-CHARGES', 6, CARPOI)
C 
C -----------------------------------------------------------------------
C 
C     Pour determiner le point de passge elastique-endommagement
C     on fait un tour a vide avec la solution chapeau calculee a partir
C     de la solution orthotrope obtenue par gadmor
C -----------------------------------------------------------------------
      EPSISO = 1
      SIGISO = 1
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
        CALL POUSMD (NTETA*NGAU1*12, DEPCHA)
C 
        DSICHA = DEPCHA+NTETA*NGAU1*6
        DBSICH = DSICHA
C 
C       recherche des caracteristiques du comportement de la couche
C 
        CALL ANGCOU (NUCOU, TETORT)
C 
        CALL COMISO (ADCOU, NUCOU, COUISO)
        CALL COMISO (ADSOU, NUCOU, SOUISO)
C 
C       BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
        DEBGAU   =  1+(NUCOU-1)*XINTEG*YINTEG*NBCOL
        FINGAU   =  NUCOU*XINTEG*YINTEG*NBCOL
C 
        CALL CNLINC (NUCOU, CARNLI, ADISOS, ADRORT)

        DO PGAU1  = DEBGAU , FINGAU
C 
C         BOUCLE ii SUR LES ANGLES
C 
          DO TETA = 1, NTETA
C 
C         RECHERCHE DE L'ANGLE DE LA BANDE (tetcal) CORRESPONDANT A teta
C 
            TETCAL = DM(ADTETA+TETA-1)
C 
C           INITIALISATION DES QUANTITES DE L'ITERATION PRECEDENTE
C 
            CALL RICLOC (COUISO, TETORT, TETCAL, KLOC)
C 
C           SOUPLESSE TANGENTE
C 
            CALL RICLOC (SOUISO, TETORT, TETCAL, SLOC)
C 
            ANGORT = TETCAL - TETORT
C 
C           CALCULD'UN COUPLE CHAPEAU AVEC UNE DIRECTION DE DESCENTE
C           -K ORTHOTROPE
C 
            CALL CHAPOR
     &         (ANGORT, CARNLI, EPSORT(EPSISO), SIGORT(SIGISO),
     &          SLOC, KLOC, EPCORT, DM(DEPCHA), DM(DSICHA), YMAXAC)
C 
C           Si le point est plus 'endommage' que les points precedents on
C           stocke ses caracteristiques
C 
            IF (YMAXAC .GT. YMAXPR) THEN
C 
              M(CARPOI)    = NUCOU
              M(CARPOI+1)  = PGAU1
              M(CARPOI+2)  = TETA
              YMAXPR       = YMAXAC
C 
CD            CALL IMPEN( ' PGAU1' ,PGAU1 )
CD            CALL IMPEN( ' TETA '  ,TETA )
CD            CALL IMPDN( ' YMAXAC ' ,YMAXAC )
CD            CALL OMPTDN( 'EPSORT ' ,EPSORT(EPSISO ) , 6, 1 )
CD            CALL OMPTDN( 'SIGORT ' ,SIGORT(EPSISO ) , 6, 1 )
CD            CALL OMPTDN( 'EPCORT ' ,EPCORT , 6 , 1 )
C 
            END IF
C 
C          <=> on stocke -2delta(epsilon-point-chapeau)
C 
            CALL MUMARE (-2.D0, 6, DM(DEPCHA), DM(DEPCHA))
C 
C           Pour avoir directement les contraintes permettant le calcul
C           des seconds membres pour l'etape globale, on multiplie par 2
C           <=> on stocke -2delta sigma-point-chapeau)
C 
            CALL ADDMAD (6, DM(DSICHA), DM(DSICHA), DM(DSICHA))
C 
            EPSISO = SIGISO+6
            SIGISO = SIGISO+6
            DEPCHA = DEPCHA+6
            DSICHA = DSICHA+6
C 
C         FIN DE BOUCLE ii SUR LES ANGLES
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
C     Test sur le nombre d'interfaces.
C     Pour le cas ou il n'y a pas d'interface.
C 
       DBSNCH = DBSICH
C 
      IF (NBINT .GT. 0) THEN
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
        YMAXPR   = 0.D0
        SINISO   = 1
        SAUISO   = 1
        CARPOI   = CARPOI+3
C 
        CALL POUSMD (NTETA*NGAU2*6, DSACHA)
C 
        DBSCHA = DSACHA
        DSNCHA = DSACHA+NTETA*NGAU2*3
        DBSNCH = DSNCHA
C 
C     BOUCLE SUR LES INTERFACES
C 
      DO NUINT = 1 , NBINT
C 
C       RECHERCHE DES CARACTERISTIQUES DU COMPORTEMENT DE LA COUCHE
C 
        CALL IOMISO (ADINT, NUINT, INTISO)
        CALL IOMISO (ADSNT, NUINT, INSISO)
        CALL ANGINT (NUINT, TETORT)
C 
        CALL CNLINI (NUINT, CARNLI, SOUORT)
C 
C       BOUCLE i SUR LES POINTS DE GAUSS DANS LA COUCHE
C 
        DEBGAU   =  1+(NUINT-1)*XINTEG*NBCOL
        FINGAU   =  NUINT*XINTEG*NBCOL
C 
        DO PGAU1  = DEBGAU , FINGAU
C 
C         BOUCLE ii SUR LES ANGLES
C 
          DO TETA = 1 , NTETA
C 
C         RECHERCHE DE L'ANGLE DE LA BANDE (tetcal) CORRESPONDANT A teta
C 
            TETCAL = DM(ADTETA+TETA-1)
            ANGORT = TETCAL-TETORT
C 
            CALL INTELA (INSISO, TETORT, TETCAL, SLOCI)
            CALL INTELA (INTISO, TETORT, TETCAL, KLOCI)
C 
C           MATRICE ELASTIQUE COMPORTEMENT TANGENT
C 
            CALL  IHAPOR
     &        (ANGORT, CARNLI, SAUORT(SAUISO), SGNORT(SINISO),
     &         SLOCI, KLOCI, SACORT, DM(DSACHA), DM(DSNCHA), YMAXAC)
C 
C     Pour avoir directement les sauts permettant le calcul
C     des seconds membres pour l'etape globale, on multiplie par 2
C     <=> on stocke -2delta(saut-point-chapeau)
C 
           CALL MUMARE (-2.D0, 3, DM(DSACHA), DM(DSACHA))
C 
C     Pour avoir directement les contraintes permettant le calcul
C     des seconds membres pour l'etape globale, on multiplie par 2
C     <=> on stocke -2delta(sigma-point-chapeau)
C 
          CALL ADDMAD (3, DM(DSNCHA), DM(DSNCHA), DM(DSNCHA))
C 
            IF (YMAXAC .GT. YMAXPR) THEN
C 
              M( CARPOI )    = NUINT
              M( CARPOI+1 )  = PGAU1
              M( CARPOI+2 )  = TETA
              YMAXPR         = YMAXAC
C 
CD            CALL IMPEN( ' NUINT' ,NUINT )
CD            CALL IMPEN( ' PGAU1' ,PGAU1 )
CD            CALL IMPEN( ' TETA '  ,TETA )
CD            CALL IMPDN( ' YMAXAC ' ,YMAXAC )
CD            CALL OMPTDN( 'SAUORT ' ,SAUORT(SAUISO ) , 3, 1 )
CD            CALL OMPTDN( 'SINORT ' ,SGNORT(SINISO ) , 3, 1 )
CD            CALL OMPTDN( 'SACORT ' ,SACORT , 3 , 1 )
C 
            END IF
C 
            SINISO   = SINISO+3
            SAUISO   = SAUISO+3
            DSNCHA   = DSNCHA+3
            DSACHA   = DSACHA+3
C 
C         FIN DE BOUCLE ii SUR LES ANGLES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES POINTS DE GAUSS DANS L'INTERFACE
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES INTERFACES
C 
      END DO
C 
C     FIN DE TEST SUR LE NOMBRE D'INTERFACES
C 
      END IF

C     1001   CONTINUE
C 
      CALL SOPOUB (AM2LC, ADM2LC)
      CALL IMPET ('ADM2 EN SORTIE DANS'//IDPROG, ADM2)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments:
C 
C     E ...... ANGORT   angle entre la base locale et base d'orthotropie
C     E ...... CARNLI   tableau des caracteristiques non lineaires de la couche
C     E ...... EPSISO   epsilon isotrope admissible (nsig) 
C     E ...... SIGISO   sigma isotrope  admissible (neps)
C 
C     (PROVIENT DE L'ITERATION PRECEDENTE)
C 
C     E ...... SOUORT   souplesse elastique orthotrope locale
C                       stocke K , B , A , C <=> 9 + 4 + 3+ 1
C     E ...... KORT     rigidite elastique orthotrope locale
C                       stocke K , B , A , C <=> 9 + 4 + 3+ 1
C     E ...... EPCORT   epsilon orthotrope chapeau
C 
C     Et on recupere :
C 
C     S ...... DEPCHA   delta epsilon point chapeau
C     S ...... DSICHA   -delta sigma point chapeau
C     S ...... YMAXAC   valeur maximale de Ybcou actuelle

      SUBROUTINE CHAPOR
     &            (ANGORT, CARNLI, EPSISO, SIGISO, SOUORT, KORT,
     &             EPCORT, DEPCHA, DSICHA, YMAXAC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C     Caracteristiques de la couche et du point
C 
      DOUBLE PRECISION ANGORT, CARNLI(12)
C 
C     Provient des champs admissibles de l'etape locale precedente
C 
C              Provient des champs admissibles
C 
      DOUBLE PRECISION EPSISO(6), SIGISO(6)
C 
C              Caracteristiques de la couche et du point
C 
      DOUBLE PRECISION SOUORT(17), KORT(17), EPCORT(6)
C 
      DOUBLE PRECISION DEPCHA(6)
      DOUBLE PRECISION DSICHA(6), YMAXAC
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER          I
      INTEGER          AM2LC, ADM2LC
      DOUBLE PRECISION EPSLOC(6)
C 
      CHARACTER*6 IDPROG
C 
CD    LOGICAL     LTRACN , LTRACP
C 
      PARAMETER (IDPROG='CHAPOR')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
CD    IF ( LTRACN(1) ) THEN
C 
CD        CALL OMPTDN('souplesse orthotrope ', SOUORT ,17 , 1 )
CD        CALL OMPTDN('rigidite orthotrope ', KORT ,17 , 1 )
C 
CD        CALL OMPTDN( ' sigma  admissible isotrope '
CD         ,SIGISO(1), 6 ,  1 )
CD        CALL OMPTDN( ' epsilon  admissible isotrope '
CD         ,EPSISO(1), 6 ,  1 )
C 
CD    END IF
C 
      CALL MULORT (SOUORT(1), SOUORT(10), SOUORT(14), SOUORT(17),
     &             SIGISO, EPSLOC)
C 
C     DEPCAC = ( eps-point-chapeau - eps-point-n) X DT
C 
      DO I = 1 , 6
        DEPCHA(I)  = .5D0*( EPSLOC(I) - EPSISO(I) )
      END DO
C   
C     DSPCAC =  -(sig-point-chapeau - sig-point-n) XDT
C 
      CALL MULORT (KORT(1), KORT(10), KORT(14), KORT(17),
     &             DEPCHA(1), DSICHA(1))
C 
      DO I = 1 , 6
        EPSLOC(I)  = .5D0*( EPSLOC(I) + EPSISO(I) )
      END DO
C 
      CALL QBORTH (ANGORT, EPSLOC, EPCORT)
C 
      CALL YELAST (CARNLI, EPCORT, YMAXAC)
C 
CD     IF (LTRACN(1)) THEN
C 
CD      CALL OMPTDN( 'Delta epsilon  chapeau en sortie '
CD       ,DEPCHA(1), 6 ,  1 )
CD      CALL OMPTDN( 'Delta sigma chapeau en sortie '
CD       ,DSICHA(1), 6 ,  1 )
C 
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
C     On envoie comme arguments :
C 
C     E ...... ANGORT   angle entre la base locale et base d'orthotropie
C     E ...... CARNLI   tableau des caracteristiques non lineaires de la couche .
C     E ...... SAUISO   saut isotrope admissible ( nsig )
C     E ...... SINISO   sigma normal isotrope  admissible ( neps )
C 
C     (PROVIENT DE L'ITERATION PRECEDENTE)
C 
C     E ...... SOUORT   souplesse elastique orthotrope locale
C     E ...... KORT     rigidite elastique orthotrope locale
C 
C     Et on recupere :
C 
C     S ...... SACORT   epsilon orthotrope chapeau
C     S .....  DSACHA   delta epsilon point chapeau
C     S ...... DSICHA   -delta sigma point chapeau
C 
      SUBROUTINE IHAPOR
     &            (ANGORT, CARNLI, SAUISO, SINISO, SOUORT, KORT,
     &             SACORT, DSACHA, DSICHA, YMAXAC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C     Provient des champs admissibles de l'etape locale precedente
C 
C              Caracteristiques de la couche et du point 
C 
      DOUBLE PRECISION SOUORT(9), KORT(9), ANGORT, CARNLI(6)
C 
C              Provient des champs admissibles
C 
      DOUBLE PRECISION SAUISO(3), SINISO(3), SACORT(3), YMAXAC
C 
      DOUBLE PRECISION DSACHA(3)
      DOUBLE PRECISION DSICHA(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER          I
      DOUBLE PRECISION SAULOC(3)
C 
      CHARACTER*6      IDPROG
C 
CD    LOGICAL     LTRACN, LTRACP
C 
      PARAMETER (IDPROG='IHAPOR')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF ( LTRACP(1) ) THEN
C 
CD      CALL OMPTDP( ' saut  admissible en entree '
CD       ,SAUISO(1), 3 ,  1 )
CD      CALL OMPTDP( ' sigma  admissible en entree '
CD       ,SINISO(1), 3 ,  1 )
CD      CALL OMPTDP('matrice tangente en entree ', KORT ,3 , 3 )
CD      CALL OMPTDP('SOUPLESSE TANGENTE EN ENTREE ',SOUORT  ,3 , 3 )
C 
CD    END IF
C 
      CALL IULORT (SOUORT, SINISO, SAULOC)
C 
C     SOL = INTERV *(2 saut-point-chapeau - saut-point-n)
C 
C     ===> MODIFICATION POUR AVOIR :
C 
C     SOL = INTERV * (2 saut-point-chapeau - saut-point-n)
C 
C     DSACAC = INTERV*( saut-point-chapeau - saut-point-n)X DT
C 
      DO I = 1 , 3
        DSACHA(I)  = .5D0*( SAULOC(I) - SAUISO(I) )
      END DO
C 
C     DSPCAC =  -(sig-point-chapeau - sig-point-n)
C 
      CALL IULORT (KORT, DSACHA(1), DSICHA(1))

C 
      DO I = 1 , 3
        SAULOC(I)  = .5D0*( SAULOC(I) + SAUISO(I) )
      END DO
C 
      CALL QIBORT (ANGORT, SAULOC, SACORT)
C 
      CALL YELASI (CARNLI, SACORT, YMAXAC)
C 
CD    IF (LTRACN(1)) THEN
C 
CD      CALL OMPTDN( 'Delta saut chapeau en sortie '
CD       ,DSACHA(1), 3 ,  1 )
CD      CALL OMPTDN( 'Delta sigma chapeau en sortie '
CD       ,DSICHA(1), 3 ,  1 )
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
C     Cette routine determine les (efforts globaux)/2 qui proviennent de
C     l'etape locale.
C 
C     On procede en au moins deux etapes :
C        - determination des differents efforts pour toutes les
C          valeurs de teta :
C        - valeurs de teta  <=> -2B delta( sign )
C        - developpement en series de Fourier de ces efforts
C 
C     On envoie comme arguments et on modifie : dm(admeg)
C 
C     E ...... DSICHA  delta sigma point chapeau
C     E ...... DSNCHA  delta sigma normale point chapeau
C     ES...... ADSMEG  l'adresse de depart des seconds membres
C                      pour l'etape globale ceux-ci sont
C                      assembles dans DM a partir de ADSMEG de
C                      la facon suivante (NDDL, NFTGLO, NBMAT)
C 
C    !!!!  seconds membres doivent etre mis a zero en entree< = >
C          dm (adsmeg, ...) = 0.d0                              !!!!
C 
C     Et on recupere :
C 
C     S ...... NBDEV  nombre de developpemts non nuls
C     S ...... TNUDEV tablleau de developpemts non nuls
C 
      SUBROUTINE EFFORT (DSICHA, DSNCHA, ADSMEG, NBDEV,TNUDEV,NORME)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION DSICHA(NEPS*NTETA*NGAU1)
C 
      DOUBLE PRECISION DSNCHA(NSAU*NTETA*NGAU2), NORME
C 
      INTEGER          ADSMEG, NBDEV, TNUDEV(NBMAT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TLOCN1, DEBGAU, DBGAUI
      INTEGER  DBDDLU, DBDDLV, DBDDLW, DBDDIU, DBDDIV, DBDDIW
      INTEGER  ASIGIN, TIGINT, DBSIGC, DBSINC
      INTEGER  ASININ, TININT
      INTEGER  NUCOU, NUCOL, ADRGAU, X, Y
      INTEGER  NUINT, H, K
      INTEGER  AM2LC, ADM2LC
C 
      LOGICAL PASSP
C 
      DOUBLE PRECISION  A, B, R, RAYONC, POIDG, WI
C 
C     Pour nettoyer la partie des efforts sans signification
C 
      DOUBLE PRECISION SCAL, TEST
C 
      INTEGER          RESU, ARESU, I, DEBUT
C 
      INTEGER          NUDEV
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EFFORT')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL TESTAD ( ADSMEG , NDDL*NBMAT , IDPROG )
C 
      K=XINTEG*(XINTEG-1)/2
      H=YINTEG*(YINTEG-1)/2
C 
C     Recherche des differents tableaux utiles
C 
      CALL ADTBM ( 'TLOCN1    ' , TLOCN1 )
      CALL ADTBDM( 'TAB-GAUSS ' , DEBGAU )
C 
C     Ouverture d'un tableau partiel pour ranger les numeros
C     des ddl des deplacements calcules (u, v, w) par developpement croissant
C 
      CALL POUSME(NDDLEL, DBDDLU)
      DBDDLV = DBDDLU+12
      DBDDLW = DBDDLV+12
C 
C     Pour les interfaces
C 
      DBDDIU = DBDDLU
      DBDDIV = DBDDIU+8
      DBDDIW = DBDDIV+8
C 
C     ASIGIN    sigma
C               puis sigma developpe ====> X 6*(NTETA+1)
C     TIGINT    transpose de sigma   ====> X 6*NTETA
C 
      CALL POUSMD( 12*NTETA+6 , ASIGIN )
C 
      TIGINT = ASIGIN+6*(NTETA+1)
C 
C     Pour les quantites de l'etape locale pour les interfaces
C 
      ASININ = ASIGIN
      TININT = TIGINT
C 
C     1ere adresse dans DSICHA et DSNCHA
C 
      DBSIGC =  1
      DBSINC =  1
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES COUCHES
C 
      DO NUCOU = 1 , NBCOU
C 
C     Recherche des caracteristiques geometriques de la couche
C 
CD        CALL IMPEN('POUR NUCOU = ' , NUCOU)
C 
          CALL VALEP( NUCOU , B )
C 
CD        CALL IMPDP( 'VALEUR DE B ', B  )
C 
C         BOUCLE i SUR LES COLONNES
C 
          DO NUCOL = 1,NBCOL
C 
CD          CALL IMPEN('POUR NUCOL = ' , NUCOL)
C 
            CALL VALRAY( NUCOL , RAYONC , A )
C 
CD          CALL IMPDP( 'VALEUR DE A ', A  )
C 
C           Recherche  des numeros de ddl ranges tableau pour l'element
C 
            CALL DDLCAL (1, NUCOU, NUCOL, TLOCN1, M(DBDDLU))
            CALL DDLCAL (2, NUCOU, NUCOL, TLOCN1, M(DBDDLV))
            CALL DDLCAL (3, NUCOU, NUCOL, TLOCN1, M(DBDDLW))

C           DEBGAU est l'adresse pour x=1, y=1 dans TAB-GAUSS
C 
            ADRGAU = DEBGAU
C 
C           BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            DO X = 1 , XINTEG
C 
              WI=POIDS(K+X)
C 
C             Calcul du rayon au point de gauss
C 
              R  = RAYONC + A*GAUSS( XINTEG*(XINTEG-1)/2 +X )
C 
C             BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT z
C 
              DO Y = 1 , YINTEG
C 
CD              CALL IMPEP ( 'VALEUR DE Y', Y )
C 
                POIDG = WI*POIDS(H+Y)
C 
C               Rangement des contraintes
C               sous la forme ( nteta , 6 )
C 
                CALL TRANSP (DSICHA (DBSIGC), DM (TIGINT), 6, NTETA)
C 
C               Calcul des valeurs developpees rangees calcul des contraintes
C               integrees sur le temps rangees dans ASIGIN
C 
                CALL SIGDEV (DM(TIGINT), DM(ASIGIN))
                CALL SMITEC (POIDG, 1, DM(ADRGAU), A, B, R,
     &                       M (DBDDLU), DM( ASIGIN), ADSMEG)
C 
C     pour aller lire dans TAB_GAUSS au bon point de GAUSS
C 
                DBSIGC      = DBSIGC+6*NTETA
C 
                ADRGAU      = ADRGAU+36
C 
C             FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT z
C 
              END DO
C 
C           FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            END DO
C 
C        FIN DE BOUCLE i SUR LES COLONNES
C 
          END DO
C 
C     FIN DE BOUCLE SUR LES COUCHES
C 
      END DO
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES INTERFACES
C 
      IF( NBINT.GT.0) THEN
C 
        CALL ADTBDM( 'GAUSSINTER' , DBGAUI)
C 
        DO NUINT = 1 , NBINT
C 
          PASSP = .FALSE.
          IF(SYMPAR .AND. NUINT .EQ. 1 ) PASSP = .TRUE.
C 
C         BOUCLE i SUR LES COLONNES
C 
          DO NUCOL = 1,NBCOL
C 
            CALL VALRAY( NUCOL , RAYONC , A )
C 
CD          CALL IMPDP( 'VALEUR DE A ', A  )
C 
C           Recherche des numeros de ddl ranges calcul pour l'element
C 
            CALL DDLICA (NUINT, NUCOL, TLOCN1, M(DBDDIU),
     &                   M(DBDDIV), M(DBDDIW))

C           DEBGAU est l'adresse pour x=1, y=1 dans TAB-GAUSS
C 
            ADRGAU = DBGAUI
C 
C           BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            DO X = 1 , XINTEG
C 
CD            CALL IMPEP ( 'VALEUR DE X ', X )
C 
              POIDG =POIDS(K+X)
C 
C             Calcul du rayon au point de gauss
C 
              R  = RAYONC + A*GAUSS( XINTEG*(XINTEG-1)/2 +X )
C 
CD            CALL IMPDP ( 'VALEUR DE R ', R )
C 
C             Rangement des contraintes normales integrees sur le temps
C             sous la forme (nteta, 3)
C  
              CALL TRANSP (DSNCHA (DBSINC), DM (TININT), 3, NTETA)
C 
C             Calcul des valeurs developpees rangees calcul des contraintes
C             normales integrees sur le temps rangees dans ASIGIN
C 
              CALL SINDEV( DM(TININT) , DM(ASININ) )
C 
              IF ( PASSP ) THEN
C 
                CALL SMITEP( POIDG , 1 , DM(ADRGAU) , A , R ,
     &                       M( DBDDIU) , DM( ASININ) , ADSMEG )
C 
              ELSE
C 
                CALL SMITEI( POIDG , 1 , DM(ADRGAU) , A , R ,
     &                       M( DBDDIU) , DM( ASININ) , ADSMEG )
C 
              END IF
C 
C             Pour aller lire dans TAB_GAUSS au bon point de GAUSS
C 
              DBSINC      = DBSINC+3*NTETA
              ADRGAU      = ADRGAU+8
C 
C           FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES COLONNES
C  
          END DO
C 
C       BOUCLE SUR LES INTERFACES
C 
        END DO
C 
      END IF
C 
C     Pour eviter les ennuis numeriques
C 
      SCAL =  0.D0
C 
      DEBUT = ADSMEG
      CALL POUSMD( NBMAT , ARESU)
CD    CALL IMPET ('ARESU dans '//IDPROG, ARESU)
      RESU = ARESU
C 
      DO I = 1 , NTDSFG
C 
        CALL SCAVEC( NDDL , DM(DEBUT) , DM(DEBUT) , DM(RESU) )
CD      CALL IMPET ('POUR I', I)
CD      CALL IMPDT ( ' DM(RESU) dans '//IDPROG, DM(RESU) )
C 
        DM(RESU) = PI*DM(RESU)
        SCAL = SCAL+ DM(RESU)
        RESU = RESU+1
        DEBUT = DEBUT + NDDL
C 
      END DO
C  
      CALL SCAVEC (NDDL, DM(DEBUT), DM(DEBUT), DM(RESU))
C 
CD    CALL IMPET ('POUR I', I)
CD    CALL IMPDT ('DM(RESU) DANS '//IDPROG, DM(RESU))
C 
      DM(RESU) = 2.D0*PI*DM(RESU)
      SCAL  = SCAL+ DM(RESU)
      RESU  = RESU+1
      DEBUT = DEBUT + NDDL
C 
      DO I = NTDSFG+2 , NBMAT
C 
        CALL SCAVEC (NDDL, DM(DEBUT), DM(DEBUT), DM(RESU))
C 
CD      CALL IMPET ('POUR I', I)
CD      CALL IMPDT ('DM(RESU) DANS '//IDPROG, DM(RESU))
C 
        DM(RESU) = PI*DM(RESU)
        SCAL  = SCAL+ DM(RESU)
        RESU  = RESU+1
        DEBUT = DEBUT + NDDL
C 
      END DO
C 
C     MODIF SCAL  = SCAL/DBLE(NBMAT)
C 
      TEST = 1.D -8*SCAL
C 
      DEBUT = ADSMEG
      NORME = 0.D0
      NBDEV = 0
C 
CD    CALL IMPET ('2 ARESU '//IDPROG, ARESU)
C 
      DO I  = ARESU , ARESU+NBMAT-1
C 
       IF(DM(I).LT.TEST) THEN
C 
           NUDEV = I-ARESU-NTDSFG
C 
CD         CALL IMPET ('On balaye effort NUDEV0 '//idprog , NUDEV)
C 
           CALL BALAID (NDDL, DM(DEBUT) )
C 
        ELSE
C 
        NBDEV           = NBDEV+1
        TNUDEV(NBDEV)   = I-ARESU+1
C 
        NORME = NORME+DM(I)
C 
        END IF
C 
CD      CALL IMPDT('Valeur relative du developpement'//
CD   &     ' numero dans '//idprog , DM(I)/SCAL)
C 
        DEBUT = DEBUT+NDDL
C 
      END DO
C 
      NORME = DSQRT( NORME)
C 
      CALL IMPTDT ('VALEUR DES DEVELOPPEMENTS ', DM(ARESU), NBMAT, 1)
C 
CD    IF (LTRACN(1))THEN
C 
CD      CALL OMPTDN ('EFFORT GLOBAL ',
CD                    DM(ADSMEG), NDDL, NBMAT)
C 
CD    END IF
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
