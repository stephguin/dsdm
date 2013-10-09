C     Cette routine calcule la contribution au point de Gauss de la
C     deformation NUEPS => DEPS est la nouvelle contrainte.
C     La matrice KT etant symetrique, on ne fait varier j que de
C     1 a i, ce qui correpond au rangement de KIJPG
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions du temps
C     E ...... NUEPS    numero de la deformation
C     E ...... EPSILO   valeur de la deformation au point de Gauss
C     E ...... K        la valeur des diferentes matrices tangentes
C                       au point de gauss
C     E ...... DBNSIG   l'adresse du debut des contraintes calculees
C                       SIGNOU
C 
C     Et on modifie et recupere :
C 
C     S ...... LONEC    la derniere adresse ecrite du tableau
C                       des nouvelles contraintes
C     ES...... SIGNOU   le tableau des contraintes au point de gauss
C 
      SUBROUTINE CTSIGC (NFONCT, NUEPS, EPSILO, K,
     &                   DBNSIG, SIGNOU, LONEC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NFONCT, NUEPS, DBNSIG, LONEC
      DOUBLE PRECISION  K(17*NFONCT*(NFONCT+1)/2)
      DOUBLE PRECISION  EPSILO(NEPS)
      DOUBLE PRECISION  SIGNOU(NGAU1*NEPS*NTETA*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  SIGMA(6)
C 
      INTEGER     I, DIJ, DECALSI, DEBSIG
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CTSIGC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
C     DCALSI decalage de SIGMA dans SIGNOU
C 
      DECALSI = NGAU1*NEPS*NTETA
      DEBSIG  = DBNSIG
C 
C     DEBUT DE LA MATRICE KIJ (J = NUEPS)
C 
      DIJ = 17*(NUEPS-1)+1
C 
C     Partie inferieure de la diagonale
C 
C 
      DO  I = 1, NUEPS-1
        CALL MULORT (K(DIJ), K(DIJ+9), K(DIJ+13), K(DIJ+16),
     &               EPSILO, SIGMA)
        CALL ADDMAD (NEPS, SIGMA, SIGNOU(DEBSIG), SIGNOU(DEBSIG))
        DIJ     = DIJ   +17*(NFONCT-I)
        DEBSIG  = DEBSIG+DECALSI
      END DO
C 
C    Partie superieure de la diagonale
C 
      DO  I =  NUEPS , NFONCT
        CALL MULORT (K(DIJ), K(DIJ+9), K(DIJ+13), K(DIJ+16),
     &               EPSILO, SIGMA)
        CALL ADDMAD (NEPS, SIGMA, SIGNOU(DEBSIG), SIGNOU(DEBSIG))
        DIJ     = DIJ   +17
        DEBSIG  = DEBSIG+DECALSI
      END DO
C 
      LONEC = DEBSIG-DECALSI+5
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule la contribution au point de Gauss de la
C     deformation NUEPS => DEPS est la nouvelle contrainte.
C     La matrice KT etant symetrique, on ne fait varier j que de
C     1 a i, ce qui correpond au rangement de KIJPG.
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions du temps
C     E ...... NUSAU    numero de la deformation
C     E ...... SAUT     valeur de la dformation au point de Gauss
C     E ...... K        la valeur des differentes matrices tangentes
C                       au point de gauss
C     E ...... DBNSIG   l'adresse du debut des contraintes calculees SIGNOU
C     E ...... LONEC    la derniere adresse  ecrite du tableau
C 
C     Et on recupere :
C 
C     ES...... SGNNOU   le tableau des contraintes au point de gauss

      SUBROUTINE CTSIGI (NFONCT, NUSAU, SAUT,
     &                   K, DBNSIG, SGNNOU, LONEC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NFONCT , NUSAU , DBNSIG , LONEC
      DOUBLE PRECISION  K(9*NFONCT*(NFONCT+1)/2)
      DOUBLE PRECISION  SAUT(NSAU)
      DOUBLE PRECISION  SGNNOU(NGAU2*NSAU*NTETA*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  SIGMA(3)
C 
      INTEGER     I , DIJ , DECALSI  , DEBSIG
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CTSIGI')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
C     DCALSI decalage de SIGMA dans SIGNOU
C 
      DECALSI = NGAU2*NSAU*NTETA
      DEBSIG  = DBNSIG
C 
C     DEBUT DE LA MATRICE KIJ (J = NUEPS)
C 
      DIJ = 9*(NUSAU-1)+1
C 
C     PARTIE INFERIEURE A LA DIAGONALE
C 
      DO  I = 1, NUSAU-1
C 
        CALL IULORT (K(DIJ), SAUT, SIGMA)
C 
        CALL ADDMAD (NSAU, SIGMA, SGNNOU(DEBSIG), SGNNOU(DEBSIG))
C 
        DIJ     = DIJ   +9*(NFONCT-I)
        DEBSIG  = DEBSIG+DECALSI
C 
      END DO
C 
C     PARTIE SUPERIEURE A LA DIAGONALE
C 
      DO  I =  NUSAU , NFONCT
C 
        CALL IULORT (  K(DIJ) , SAUT , SIGMA )
C 
        CALL ADDMAD( NSAU , SIGMA , SGNNOU(DEBSIG ) ,SGNNOU(DEBSIG ) )
C 
        DIJ     = DIJ   +9
        DEBSIG  = DEBSIG+DECALSI
C 
      END DO
C 
      LONEC = DEBSIG-DECALSI+2
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACP(1))THEN
C 
CD      CALL OMPTDP( ' SGNNOU modifie', SGNNOU(1) ,DECALSI ,NFONCT )
C 
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Attention !!! Cette routine renvoie dans dm (adsmeg) la
C     difference avec dm (adsmeg) de depart
C 
C     Cette routine determine les nouveaux efforts globaux qui proviennent de
C     l'iteration BFGS. Il y en a autant que de fonctions du temps a
C     l'etape globale. On procede en (au moins!) deux etapes :
C 
C       - determination des differents efforts pour toutes les
C         valeurs de teta  <=> B ( NOUSIG )
C       - developpement en serie de Fourier de ces efforts
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT  le nombre de fonctions du temps de l'iteration
C     E ...... SIGNOU  tableau des nouvelles contraintes
C     E ...... SGNNOU  tableau des nouvelles contraintes normales
C 
C     Et on recupere le tableau dans DM ecrit a partir de ADSMEG :
C 
C     ES...... ADSMEG  l'adresse de depart des seconds membres
C                      pour le nouveau calcul global; ceux-ci sont
C                      assembles dans DM a partir de ADSMEG de
C                      la facon suivante : (NDDL, NFONCT, NBMAT)
C     S....... NORME   la norme de ces efforts

      SUBROUTINE NOUEFF (NFONCT, SIGNOU, SGNNOU, ADSMEG,
     &                   NBDEV, TNUDEV, NORME)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER          ADSMEG, NFONCT, NBDEV, TNUDEV(NBMAT)
      DOUBLE PRECISION SIGNOU(NEPS*NTETA*NGAU1*NFONCT), NORME
      DOUBLE PRECISION SGNNOU(NSAU*NTETA*NGAU2*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TLOCN1, DEBGAU, DBGAUI, AFTGLO, DFTGLO
      INTEGER  DBDDLU, DBDDLV, DBDDLW, DBDDIU, DBDDIV, DBDDIW
      INTEGER  NFT, NUCOU, NUCOL, ADRGAU, X, Y
      INTEGER  NUINT, LONTEC, LONTEI, DSIDEV
      INTEGER  DSGNTR, DSNDEV, DBSIGM, DBSIGN, DSIGTR
      INTEGER  AM2LC, ADM2LC, H, K
C 
      LOGICAL PASSP
C 
      DOUBLE PRECISION  A, B, R, RAYONC, POIDG, WI
C 
CD    LOGICAL  LTRACN, LTRACP
C 
C     pour nettoyer la partie des efforts sans signification
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NOUEFF')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL TESTAD ( ADSMEG , NDDL*NFTGLO*NBMAT , IDPROG )
C 
      K=XINTEG*(XINTEG-1)/2
      H=YINTEG*(YINTEG-1)/2
C 
C 
C     RECHERCHE DES DIFFERENTS TABLEAUX UTILES
C 
      CALL ADTBM ( 'TLOCN1    ' , TLOCN1 )
      CALL ADTBDM( 'TAB-GAUSS ' , DEBGAU )
C 
      DFTGLO      = AFTGLO
C 
C     Ouverture d'un tableau partiel pour ranger les numeros
C     des ddl des deplacements ranges calcul (u,v,w) par developpement croissant
C 
      CALL POUSME(NDDLEL, DBDDLU)
      DBDDLV = DBDDLU+12
      DBDDLW = DBDDLV+12
C 
C     pour les interfaces
C  
      DBDDIU = DBDDLU
      DBDDIV = DBDDIU+8
      DBDDIW = DBDDIV+8
C 
C     Ouverture d'un tableau partiel pour ranger les valeurs transposees
C     puis developpees des nouvelles contrainte au point de gauss pour
C     tout teta pour les couches puis pour les interfaces
C 
      LONTEC   = NTETA*NSIG
      LONTEI   = NTETA*NSAU
      CALL POUSMD( LONTEC+NTDSFG*(NSIG+NSAU)+LONTEI , DSIGTR  )
C 
      DSIDEV   = DSIGTR+LONTEC
C 
      DSGNTR   = DSIDEV+NTDSFG*NSIG
C 
      DSNDEV   = DSGNTR+LONTEI
C 
C     BOUCLES SUR LES FONCTIONS DE TEMPS GLOBALES
C 
      DBSIGM = 1
      DBSIGN = 1
C 
      DO NFT = 1, NFONCT
C 
C       adresse de depart de la fonction en temps
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C       BOUCLE i SUR LES COUCHES
C 
        DO NUCOU = 1 , NBCOU
C 
C         RECHERCHE DES CARACTERISTIQUES GEOMETRIQUES DE LA COUCHE
C 
C         CALL IMPEN ('POUR NUCOU = ', NUCOU)
C 
          CALL VALEP (NUCOU, B)
C 
C         CALL IMPDP ('VALEUR DE B', B)
C 
C           BOUCLE ii SUR LES COLONNES
C 
            DO NUCOL = 1,NBCOL
C 
C           CALL IMPEN ('POUR NUCOL = ', NUCOL)
C 
            CALL VALRAY (NUCOL, RAYONC, A)
C 
C           CALL IMPDP( 'VALEUR DE A ', A  )
C 
C           RECHERCHE DES NUMEROS DE DDL RANGES CALCUL POUR L'ELEMENT
C 
            CALL DDLCAL (1, NUCOU, NUCOL, TLOCN1, M(DBDDLU))
            CALL DDLCAL (2, NUCOU, NUCOL, TLOCN1, M(DBDDLV))
            CALL DDLCAL (3, NUCOU, NUCOL, TLOCN1, M(DBDDLW))

C           DEBGAU EST L'ADRESSE POUR x=1, y=1 DANS TAB-GAUSS
C 
            ADRGAU = DEBGAU
C 
C             BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
              DO X = 1 , XINTEG
C 
C             CALL IMPEP ( 'VALEUR DE X ', X )
C 
              WI=POIDS(K+X)
C 
C             CALCUL DU RAYON AU POINT DE GAUSS
C 
              R  = RAYONC + A*GAUSS( XINTEG*(XINTEG-1)/2 +X )
C 
C             CALL IMPDP ('VALEUR DE R ', R)
C 
C               BOUCLE iv SUR LES POINTS DE GAUSS SUIVANT z
C 
                DO Y = 1 , YINTEG
C 
C               CALL IMPEP ( 'VALEUR DE Y', Y )
C 
                POIDG = WI*POIDS(H+Y)
C 
C               RANGEMENT DES CONTRAINTES SOUS LA FORME (nteta, 6)
C 
                CALL TRANSP (SIGNOU(DBSIGM), DM(DSIGTR), 6, NTETA)
                CALL HOMAT (-1.D0, DM(DSIGTR), DM(DSIGTR),
     &                       6, NTETA)
C 
                DBSIGM = DBSIGM+LONTEC
C 
C               CALCUL DES VALEURS DEVELOPPEES RANGEES CALCUL DES CONTRAINTES
C               INTEGREES SUR LE TEMPS RANGEES DANS DSIDEV
C 
                CALL SIGDEV (DM(DSIGTR), DM(DSIDEV))
                CALL SMITEC (POIDG, NFT, DM(ADRGAU), A, B, R,
     &                       M(DBDDLU), DM(DSIDEV), ADSMEG)
C 
C               POUR ALLER LIRE DANS TAB_GAUSS AU BON POINT DE GAUSS
C 
                ADRGAU      = ADRGAU+36
C 
C               FIN DE BOUCLE iv SUR LES POINTS DE GAUSS SUIVANT z
C 
                END DO
C 
C             FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
              END DO
C 
C           FIN DE BOUCLE ii SUR LES COLONNES
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES COUCHES
C 
          END DO
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
C       BOUCLE i SUR LES INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
CD      CALL IMPMN (' ***** POUR LES INTERFACES ***** ' )
C 
        CALL ADTBDM ('GAUSSINTER', DBGAUI)
C 
        DO NUINT = 1 , NBINT
C 
        PASSP = .FALSE.
        IF (SYMPAR .AND. NUINT .EQ. 1) PASSP = .TRUE.
C 
C         BOUCLE ii SUR LES COLONNES
C 
          DO NUCOL = 1,NBCOL
C 
          CALL VALRAY (NUCOL, RAYONC, A)
C 
C         CALL IMPDP ('VALEUR DE A ', A)
C 
C         RECHERCHE DES UMEROS DE DDL RANGES CALCUL POUR L'ELEMENT
C 
          CALL DDLICA (NUINT, NUCOL, TLOCN1, M(DBDDIU),
     &                   M(DBDDIV), M(DBDDIW))
C 
C         DEBGAU EST L'ADRESSE POUR x=1, y=1 DANS TAB-GAUSS
C 
          ADRGAU = DBGAUI
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C  
            DO X = 1 , XINTEG
C 
CD          CALL IMPEP ( 'VALEUR DE X ', X )
C 
            POIDG=POIDS(K+X)
C 
C           CALCUL DU RAYON AU POINT DE GAUSS
C 
            R  = RAYONC + A*GAUSS( XINTEG*(XINTEG-1)/2 +X )
C 
C           CALL IMPDP ('VALEUR DE R ', R)
C 
C           RANGEMENT DES CONTRAINTES NORMALES INTEGREES SUR LE TEMPS
C           SOUS LA FORME (nteta, 3)
C 
            CALL TRANSP (SGNNOU(DBSIGN), DM(DSGNTR), 3, NTETA)
            CALL HOMAT (-1.D0, DM(DSGNTR), DM(DSGNTR),
     &                     3, NTETA)
            DBSIGN = DBSIGN+LONTEI
C 
C           CALCUL DES VALEURS DEVELOPPEES RANGEES CALCUL DES CONTRAINTES
C           NORMALES RANGEES DANS DSNDEV
C 
            CALL SINDEV (DM(DSGNTR), DM(DSNDEV))
C 
            IF (PASSP) THEN
              CALL SMITEP (POIDG, NFT, DM(ADRGAU), A, R,
     &                       M(DBDDIU), DM(DSNDEV), ADSMEG)
            ELSE
              CALL SMITEI (POIDG, NFT, DM(ADRGAU), A, R,
     &                       M(DBDDIU), DM(DSNDEV), ADSMEG)
            ENDIF
C 
C           POUR ALLER LIRE DANS TAB_GAUSS AU BON POINT DE GAUSS
C 
            ADRGAU      = ADRGAU+8
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            END DO
C  
C         FIN DE BOUCLE ii SUR LES COLONNES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES INTERFACES
C 
        END DO
C 
        END IF
C 
C     FIN DE BOUCLE SUR LES FONCTIONS DE TEMPS GLOBALES
C 
      END DO
C 
C     MISE A ZERO DES TERMES NON SIGNIFICATIFS EN EFFORT
C 
      CALL MZBLOC (NFONCT, DM(ADSMEG), NBDEV, TNUDEV(1), NORME)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
