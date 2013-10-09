C     Cette routine va chercher l'adresse de debut de tableau pour la
C     matrice de souplesse orthotrope du pli N0

C     On envoie comme arguments :
C 
C     E ...... N0 numero de couche
C 
C     Et on recupere :
C 
C     S ...... M0 1ere adresse du tableau de la souplesse
C                 orthotrope de la couche
C 
      SUBROUTINE SOCORT (N0, KORT)
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
      INTEGER N0,M0,M2,P,R
C 
      DOUBLE PRECISION KORT(17)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='COMCOU')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('SOUP-ORTHO', R)
      CALL ADTBM ('TYP-COUCHE', M2)
C 
C     P caracterise le type de comportement
C 
      P=M(M2+N0-1)
      M0=R+17*(P-1)
C 
      CALL COPITD (17, DM(M0), KORT)
C 
CD    CALL IMPEP  ('VALEUR DU TYPE DE COMPORTEMENT ',P)
CD    CALL IMPEP  ('VALEUR DE LA IERE ADRESSE M0   ',M0)
CD    CALL IMPTDP ('VALEUR DE LA SOUPLESSE         ', KORT, 1, 17)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine va chercher la direction
C     d'orthotropie de la couche consideree
C 
C     On envoie comme arguments :
C 
C     E ...... N0 numero de couche
C 
C     Et on recupere :
C 
C     S ...... TETA angle de la couche
C 
      SUBROUTINE ANGCOU (N0, TETA)
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
      INTEGER N0, R
C 
      DOUBLE PRECISION TETA
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ANGCOU')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ADTBDM ('ANGLES-COU', R)
      TETA=DM(R+N0-1)
C 
CD    CALL IMPDP('VALEUR DE L''ANGLE',TETA)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     MULTIPILICATION PAR UNE MATRICE SYMETRIQUE (3,3)
C 
C     On envoie comme arguments :
C 
C     E ...... K(9)     comportement travaillant avec les sauts
C     E ...... SAUT     Valeur des sauts
C 
C     Et on recupere :
C 
C     S ...... SIGMAN   Valeur des contraintes normales
C 
      SUBROUTINE MUMS33 (K, SAUT, SIGMAN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION  K(3,3), SAUT(3), SIGMAN(3)
C 
C -----------------------------------------------------------------------
      SIGMAN(1) = K(1,1) * SAUT(1) + K(1,2) * SAUT(2) + K(1,3) * SAUT(3)
      SIGMAN(2) = K(1,2) * SAUT(1) + K(2,2) * SAUT(2) + K(2,3) * SAUT(3)
      SIGMAN(3) = K(1,3) * SAUT(1) + K(2,3) * SAUT(2) + K(3,3) * SAUT(3)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule TRACE = TR (EPS2 K(EPS1))
C 
C     On envoie comme arguments :
C 
C     E ...... K     la matrice stockee (17)
C     E ...... EPS1  la 1ere deformation (6)
C     E ...... EPS2  la 2eme deformation (6)
C 
C     Et on recupere :
C 
C     S ...... TRACE

      SUBROUTINE TRAC12 (K, EPS1, EPS2, TRACE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  TRACE, K(17), EPS1(NEPS), EPS2(NEPS)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION SIGMA(6)
C 
      INTEGER    I
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TRAC12')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL MULORT (K(1), K(10), K(14), K(17), EPS1(1), SIGMA(1))
C 
C     CALCUL DE LA TRACE
C 
      TRACE = 0.D0
C 
      DO I = 1, 6
C 
         TRACE = TRACE + EPS2(I)*SIGMA(I)
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
C     Cette routine calcule les quantites (contraintes ou deformations)
C     dans la base d'orthotropie a partir des quantites connues dans la base
C     locale; ANGORT est l'angle entre la base locale et la base d'orthotropie
C 
C     On envoie comme arguments :
C 
C     E ...... ANGORT  angle entre la base d'orthotropie et la base locale
C     E ...... QPLABL  quantite chapeau dans la base locale
C 
C     Et on recupere :
C 
C     S ...... QPLABO  quantite chapeau dans la base d'orthotropie
C 
      SUBROUTINE QBORTH (ANGORT, QPLABL, QPLABO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION ANGORT, QPLABL(6), QPLABO(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION COS, SIN, C2, S2, SC, R2
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='QBORTH')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      COS   = DCOS(ANGORT)
      SIN   = DSIN(ANGORT)
      R2    = DSQRT( 2.D0 )
      C2    = COS*COS
      S2    = SIN*SIN
C 
C     attention sc * r2
C 
      SC    = R2*SIN*COS
C 
C     MODIF DU 25/07/96
C 
      QPLABO( 1 ) =  C2*QPLABL( 1 )+S2*QPLABL( 2 )-SC*QPLABL( 3 )
      QPLABO( 2 ) =  S2*QPLABL( 1 )+C2*QPLABL( 2 )+SC*QPLABL( 3 )
      QPLABO( 3 ) =  SC*( QPLABL(1)-QPLABL(2) )+(C2-S2)*QPLABL(3)
      QPLABO( 4 ) =  COS* QPLABL(4)+SIN*QPLABL(5)
      QPLABO( 5 ) =  COS* QPLABL(5)-SIN*QPLABL(4)
      QPLABO( 6 ) =  QPLABL(6)
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
C     E ...... DIM   dimension des matrices
C     E ...... MAT1  1ere matrice
C     E ...... MAT2  2ene matrice (D pour double precision)
C 
C     Et on recupere :
C 
C     S ...... MATSOU=MAT1-MAT2
C 
      SUBROUTINE SOUMAD (DIM, MAT1, MAT2, MATSOU)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  DIM,I
C 
      DOUBLE PRECISION MAT1(DIM),MAT2(DIM),MATSOU(DIM)
C 
C -----------------------------------------------------------------------
      DO I=1, DIM
        MATSOU(I) = MAT1(I) - MAT2(I)
      ENDDO
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine additionne deux matrices. C'est quand meme balaise ce
C     quon fait avec dsdm.
C 
C     On envoie comme arguments :
C 
C     E ...... DIM dimensions des matrices
C     E ...... MAT1 1ere matrice
C     E ...... MAT2 2ene matrice (D pour double precision)
C 
C     Et on recupere :
C 
C     S ...... MATADD = MAT1+MAT2
C 
      SUBROUTINE ADDMAD (DIM, MAT1, MAT2, MATADD)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  DIM,I
C 
      DOUBLE PRECISION MAT1(DIM),MAT2(DIM),MATADD(DIM)
C 
C -----------------------------------------------------------------------
      DO I=1,DIM
        MATADD(I)=MAT1(I)+MAT2(I)
      ENDDO
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine multiplie une matrice par un reel
C 
C     On envoie comme arguments :
C 
C     E ...... REEL   terme qui multiplie la matrice,
C     E ...... DIM    dimension de la matrice et MAT, 1ere adresse de la matrice
C 
C     Et on recupere :
C 
C     S ...... MATMU  nouvelle matrice
C 
      SUBROUTINE MUMARE (REEL, DIM, MAT, MATMU)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      INTEGER  DIM, I
C 
      DOUBLE PRECISION  MAT(DIM), MATMU(DIM), REEL
C 
C -----------------------------------------------------------------------
      DO I = 1, DIM
        MATMU(I) = REEL*MAT(I)
      ENDDO
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Multiplication Couche Orthotrope dans la base D'ORthotropie
C 
C     On envoie comme arguments :
C 
C     E ...... K(3,3)   Partie du comportement travaillant sur
C                       la partie plane des deformations
C     E ...... B(2,2)   Partie du comportement travaillant sur
C                       le cisaillement normal
C     E ...... A(3)     Partie du comportement reliant la partie plane
C                       des deformations a la deformation normale
C     E ...... C        Relie la contrainte a la deformation normale
C     E ...... EPSILO   Valeur des deformations
C 
C     Et on recupere :
C 
C     S ...... SIGMA    Valeur des contraintes
C 
C     EPSILO DIFFERENT DE SIGMA

      SUBROUTINE MCOROR (K, EPSILO, SIGMA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  K(17), EPSILO(6), SIGMA(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    LOGICAL LTRACP
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='MCOROR')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      SIGMA(1) = K(1) * EPSILO(1) + K(2) * EPSILO(2)
     &           + K(14) * EPSILO(6)
      SIGMA(2) = K(2) * EPSILO(1) + K(5) * EPSILO(2)
     &           + K(15) * EPSILO(6)
      SIGMA(3) =  K(9) * EPSILO(3)
C 
      SIGMA(4) = K(10) * EPSILO(4)
      SIGMA(5) = K(11) * EPSILO(5)
C 
      SIGMA(6) = K(14)*EPSILO(1) + K(15)*EPSILO(2)+ K(17)*EPSILO(6)
C 
CD    IF ( LTRACP(1) ) THEN
CD      CALL OMPTDP ('K EN ENTREE ' , K(1), 3, 3)
CD      CALL OMPTDP ('B EN ENTREE ' , K(10), 2, 2)
CD      CALL OMPTDP ('A EN ENTREE ' , K(14), 3, 1)
CD      CALL IMPDP  ('C EN ENTREE ' , K(17))
CD      CALL OMPTDP ('VECTEUR EN ENTREE ', EPSILO(1), 6, 1)
CD      CALL OMPTDP ('VECTEUR EN SORTIE ', SIGMA (1), 6, 1)
CD    END IF
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
C     E ...... K(3,3)   Partie du comportement travaillant sur
C                       la partie plane des deformations
C     E ...... B(2,2)   Partie du comportement travaillant sur
C                       le cisaillement normal
C     E ...... A(3)     Partie du comportement reliant
C                       la partie plane des deformations a la
C                       deformation normale
C     E ...... C        Relie la contrainte a la deformation normale
C     E ...... EPSILO   Valeur des deformations
C 
C     Et on recupere :
C 
C     S ...... SIGMA    Valeur des contraintes
C 
C     EPSILO DIFFERENT DE SIGMA
C 
      SUBROUTINE MULORT (K, B, A, C, EPSILO, SIGMA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  K(3,3), B(2,2), A(3), C, EPSILO(6), SIGMA(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    LOGICAL LTRACP
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='MULORT')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      SIGMA(1) = K(1,1) * EPSILO(1) + K(2,1) * EPSILO(2) +
     &           K(3,1) * EPSILO(3) + A(1)   * EPSILO(6)
      SIGMA(2) = K(1,2) * EPSILO(1) + K(2,2) * EPSILO(2) +
     &           K(3,2) * EPSILO(3) + A(2)   * EPSILO(6)
      SIGMA(3) = K(1,3) * EPSILO(1) + K(2,3) * EPSILO(2) +
     &           K(3,3) * EPSILO(3) + A(3)   * EPSILO(6)
C 
      SIGMA(4) = B(1,1) * EPSILO(4) + B(2,1) * EPSILO(5)
      SIGMA(5) = B(1,2) * EPSILO(4) + B(2,2) * EPSILO(5)
C 
      SIGMA(6) = A(1)*EPSILO(1) + A(2)*EPSILO(2)+A(3)*EPSILO(3)
     &           + C*EPSILO(6)
C 
CD     IF ( LTRACP(1) ) THEN
CD       CALL OMPTDP( ' K en entree ' , K(1,1) , 3 ,3 )
CD       CALL OMPTDP( ' B en entree ' , B(1,1) , 2 , 2)
CD       CALL OMPTDP( ' A en entree ' , A(1)   , 3 , 1)
CD       CALL IMPDP(  ' C en entree ' , C )
CD       CALL OMPTDP( ' vecteur en entree ', EPSILO(1) , 6 , 1)
CD       CALL OMPTDP( ' vecteur en sortie ', SIGMA (1) , 6 , 1)
CD     END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Inversion d'une matrice symetrique 3,3
C 
C     On envoie comme arguments :
C 
C     E................ K(3,3)
C 
C     Et on recupere :
C 
C     S................ INVK(3,3)

      SUBROUTINE INVMS3 (K, INVK)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C  
      DOUBLE PRECISION K(3,3), INVK(3,3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      DOUBLE PRECISION DETK, INVDET, COF(3)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='INVMS3')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      COF(1)          =  ( K(2,2)*K(3,3)-K(3,2)*K(3,2) )
      COF(2)          =  -( K(2,1)*K(3,3)-K(3,2)*K(3,1) )
      COF(3)          =  ( K(2,1)*K(3,2)-K(2,2)*K(3,1) )
C 
      DETK            = K(1,1)*COF(1)+ K(2,1)*COF(2)
     &                  +K(3,1)*COF(3)
      INVDET          = 1.D0/DETK
C 
      INVK (1,1)      = INVDET*COF(1)
      INVK (2,1)      = INVDET*COF(2)
      INVK (3,1)      = INVDET*COF(3)
      INVK (2,2)      = INVDET*( -K(3,1)*K(3,1)+K(1,1)*K(3,3))
      INVK (3,2)      = -INVDET*( K(1,1)*K(3,2)-K(2,1)*K(3,1))
      INVK (3,3)      = -INVDET*( K(1,2)*K(1,2)-K(2,2)*K(1,1))
      INVK (2,3)      = INVK(3,2)
      INVK (1,2)      = INVK(2,1)
      INVK (1,3)      = INVK(3,1)
C 
C 
CD     IF(LTRACP(1))THEN
C 
CD      CALL OMPTDP( ' MATRICE 3,3 EN ENTREE', K(1,1),3,3   )
CD      CALL OMPTDP( ' MATRICE 3,3 INVERSEE ', INVK(1,1),3,3)
CD      CALL MATM02( K , INVK , ID , 3 , 3 ,3 )
CD      CALL IMPTDT( ' MATRICE IDENTITE(3,3)? ', ID(1,1),3,3)
C 
CD     END IF
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine effectue la trace deformations ALPHAI par les
C     contraintes SIGMA.
C 
C     On envoie comme arguments :
C 
C     E ...... ALPHAI    deformations rangees (NEPS, NTETA, NGAU1)
C     E ...... SAUTI     sauts ranges (NSAU, NTETA, NGAU2)
C      
C     E ...... SIGMA     contraintes rangees (NEPS, NTETA, NGAU1)
C     E ...... SIGNOR    contraintes normales rangees (NSAU, NTETA, NGAU2)
C 
C     Et on recupere :
C 
C     S ...... TRACE    la trace
C     S ...... TRALOC   le sup des traces par element

      SUBROUTINE TRAGLO (ALPHAI, SAUTI, SIGMA, SIGNOR,
     &                   TRACE, TRALOC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  TRACE, TRALOC
C 
      DOUBLE PRECISION  ALPHAI(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION  SAUTI (NSAU*NTETA*NGAU2)
C 
      DOUBLE PRECISION  SIGMA(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION  SIGNOR(NSAU*NTETA*NGAU2)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION SCAL, A, B, JLOC
      DOUBLE PRECISION RLOC, RAYONC, MULTI, MULTJ
C 
      INTEGER     DEBPG, DECAlP, NUCOU
      INTEGER     TETA , NUCOL
      INTEGER     DEBUT, NUINT
      INTEGER     AM2LC, ADM2LC, HINT, KINT
      INTEGER     ADTET
      INTEGER     AVTET, X, Y
C 
CD    LOGICAL            LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TRAGLO')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
CD     IF (LTRACP(1)) THEN
C 
CD       CALL OMPTDP ('ALPHAI', ALPHAI(1), NEPS*NTETA, NGAU1)
CD       CALL OMPTDP ('SAUTI ', SAUTI(1), NSAU*NTETA, NGAU2)
CD       CALL OMPTDP ('SIGMA ', SIGMA(1), NEPS*NTETA, NGAU1)
CD       CALL OMPTDP ('SIGNOR', SIGNOR(1), NSAU*NTETA, NGAU2)
C 
CD     END IF
C 
C     Pour aller lire les points de Gauss
C 
      HINT = XINTEG*(XINTEG-1)/2
      KINT = YINTEG*(YINTEG-1)/2
C 
C     decalage pour passer d'une fonction ALPHAI a la suivante
C 
      DECALP = NEPS*NTETA*NGAU1
C 
C     Debut pour la 1ere fonction ALPHAI des valeurs de epsilon au point de Gauss
C 
      DEBPG = 1
C 
C     ADTET tableau en teta X NTETA
C 
      CALL GSPOUD (NTETA, ADTET)
C 
      AVTET    = ADTET-1
C 
      TRALOC = 0.D0
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
        CALL VALEP( NUCOU , B )
C 
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1, NBCOL
C 
          CALL VALRAY (NUCOL, RAYONC, A)
          JLOC    = A*B
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
          DO X = 1, XINTEG
C 
            RLOC   = RAYONC + A*GAUSS(HINT + X)
            MULTI  = POIDS(HINT+X)*JLOC*RLOC
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT z
C 
            DO Y = 1, YINTEG
C 
              MULTJ  = POIDS(KINT+Y)*MULTI
C 
C 
C             BOUCLE iv SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
C               DEFINITION DES DEBUTS DE TABLEAUX DES DEFORMATIONS
C 
                DEBUT = AVTET+TETA
C 
C               DEBUT DE BOUCLE v SUR LES DEUXIEMES FONCTIONS DU TEMPS
C 
                CALL SCAVEC (NEPS, ALPHAI(DEBPG), SIGMA(DEBPG), SCAL)
C 
                DM(DEBUT) = DM(DEBUT)+ MULTJ*SCAL
C 
                TRALOC = MAX (TRALOC, SCAL)
C 
                DEBUT  = DEBUT+NTETA
C 
C               PASSAGE AUX DEFORMATIONS DU POINT DE GAUSS (TETA compris) SUIVANT
C 
                DEBPG  = DEBPG+NEPS
C 
C             FIN DE BOUCLE iv SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT z
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
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
C     DECLAGE POUR PASSER D'UNE FONCTION SAUTI A LA SUIVANTE
C 
      DECALP = NSAU*NTETA*NGAU2
C 
C     DEBUT POUR LA 1ere FONCTION SAUTI DES VALEURS DU SAUT AU POINT DE GAUSS
      DEBPG = 1
C 
C     ADTET tableau en teta X TETTOT
C 
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT = 1, NBINT
C 
C         BOUCLE i SUR LES COLONNES
C 
          DO NUCOL = 1, NBCOL
C 
            CALL VALRAY (NUCOL, RAYONC, A)
C 
C           BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            DO X = 1, XINTEG
C 
              RLOC  = RAYONC + A*GAUSS(HINT + X)
              MULTI  = POIDS(HINT+X)*A*RLOC
C 
C             BOUCLE iii SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
                DEBUT = AVTET+TETA
C 
                CALL SCAVEC (NSAU, SAUTI(DEBPG), SIGNOR(DEBPG), SCAL)
C 
                DM(DEBUT) = DM(DEBUT)+ MULTI*SCAL
C 
                 DEBUT  = DEBUT+NTETA
C 
C                PASSAGE AU SAUT DU POINT DE GAUSS (TETA COMPRIS) SUIVANT
C 
                DEBPG  = DEBPG+NSAU
C 
C             FIN DE BOUCLE iii SUR LES ANGLES
C 
              END DO
C 
C          FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
           END DO
C 
C       FIN DE BOUCLE i SUR LES COLONNES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES INTERFACES
C 
      END DO
C 
C     SEQUENCE D'INTEGRATION EN TETA ET CALCUL DE LA CONTRIBUTION DU POINT DE GAUSS
C 
      CALL ITBTET (1, DM(ADTET), TRACE)
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    CALL IMPDN ('TRACE                  ', TRACE)
CD    CALL IMPDN ('SUP DES TRACES LOCALES ', TRALOC)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C 
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C    ...on envoie comme arguments le numero de couches NC
C    et on recupere...EPAIS, l'epaisseur de la couche divisee par 2
C 
      SUBROUTINE VALEP(NC, EPAIS )
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
       INTEGER NC ,A1 ,AR1
C  
       DOUBLE PRECISION EPAIS
C  
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VALEP')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('EPAISSEURS',A1)
      AR1=A1-1
      EPAIS=DM(AR1+NC)
C 
CD    CALL IMPEP('POUR LA COUCHE NUMERO :',NC)
CD    CALL IMPDP('VALEUR DE L''EPAISSEUR',EPAIS)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C    On envoie comme arguments :
C 
C    E ...... NUCOL    numero de la colonne
C 
C    Et on recupere :
C 
C    S ...... RAYONC   Le rayon au centre de l'element
C    S ...... LONG     La demi-longueur de l'element
C 
      SUBROUTINE VALRAY (NUCOL, RAYONC, LONG)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      INTEGER  NUCOL
C 
      DOUBLE PRECISION  RAYONC, LONG
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   ADRAY, ADLON
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VALRAY')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('RAYON     ', ADRAY)
C 
CD    CALL IMPEP('1ERE ADRESSE DU TABLEAU RAYON ', ADRAY)
CD    CALL IMPDP('1ERE VALEUR DU TABLEAU RAYON ', DM(ADRAY))
C 
      RAYONC   = DM(ADRAY+NUCOL-1)
C 
CD    CALL IMPEP('POUR LA COLONNE NUMERO :',NUCOL)
CD    CALL IMPEP('VALEUR NUCOL-1',(NUCOL-1))
CD    CALL IMPEP('VALEUR NUCOL-1',(NUCOL-1))
CD    CALL IMPDP('VALEUR DU RAYON AU CENTRE DANS :',RAYONC)
C 
      CALL ADTBDM ('LONGUEURS ', ADLON)
      LONG     = DM(ADLON+NUCOL-1)
C 
CD    CALL IMPDP(
CD    'VALEUR DE LA LONGUEUR DANS :'//IDPROG,LONG)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine effectue le produit scalaire des 2 vecteurs vect1
C     et vect2 de longueur LONG => = scal.
C 
C     On envoie comme arguments :
C 
C     E ...... LONG   la longueur des vecteurs
C     E ...... VECT1  le 1er vecteur
C     E ...... VECT2  le second vecteur
C 
C     Et on recupere :
C 
C     S ...... SCAL   le produit scalaire des 2 vecteurs
C 
      SUBROUTINE SCAVEC (LONG, VECT1, VECT2, SCAL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      INTEGER             LONG
      DOUBLE PRECISION    VECT1(LONG), VECT2(LONG), SCAL
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         I
C 
      SCAL = 0.D0
C 
      DO I = 1, LONG
C 
        SCAL = SCAL + VECT1(I)*VECT2(I)
C 
      END DO
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule l'integrale en teta d'un tableau range (LONFON, NTETA)
C 
C     On envoie comme arguments :
C 
C     E ...... LONFON  la 1ere dimension du tableau
C     E ...... FIJTET  le tableau (NTETA, LONFON)
C 
C     Et on recupere :
C 
C     S ...... TLKJI   le tableau integre range (LONFON)
C 
C     Attention : normalement tableau a 4 indices !!!
C 
      SUBROUTINE ITBTET (LONFON, FIJTET, TLKJI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER      LONFON
C 
      DOUBLE PRECISION  FIJTET(LONFON * NTETA)
      DOUBLE PRECISION  TLKJI(LONFON)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         DEBUT, I
C 
CD    LOGICAL         LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ITBTET')
C 
CD    CALL WLKBCD(IDPROG)
C 
C     CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     SEQUENCE D'IMPRESSION EN TRACE PROFONDE
C 
CD    IF (LTRACN(1))THEN
C 
CD      CALL OMPTDN( ' FIJTET ',
CD                     FIJTET(1) , NTETA ,LONFON  )
C 
CD    END IF
C 
      DEBUT = 1
C 
      DO I = 1 , LONFON
C 
C       CALL IMPEP ('TERME NUMERO ', I)
C 
C       CALL OMPTDP ('VECTEUR EXTRAIT ', DM(VECT), LONG, 1)

        CALL INTETA (FIJTET(DEBUT), TLKJI(I))
C 
        DEBUT = DEBUT+NTETA
C 
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1)) THEN
C 
CD      CALL OMPTDN ('TABLEAU INTEGRE ',
CD                    TLKJI(1), LONFON, 1)
C 
CD    END IF
C 
C     CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
         SUBROUTINE INTETA( FTETA , INTF )
C -
C -    Cette routine calcule le'integrale en teta d'une fonction
C -  en supposant un developpement en serie de Fourier
C -
C On envoie comme arguments:
C E................ FTETA    (nteta  ) la valeur de la fonction
C                            aux nteta angles
C !!!!!  ATTENTION !!! FTETA est modifie
C Et on recupere:
C S................ INTF     valeur de l'integrale
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
C
      include 'cominc.h'
C
      DOUBLE PRECISION   FTETA( NTETA ) , INTF
C
C**********************************************************************
C
C *   Declaration des parametres locaux
C     """""""""""""""""""""""""""""""""
C
      INTEGER   INV  , FIMA
      INTEGER   AM2LC , ADM2LC
CD    LOGICAL  LTRACN , LTRACP
C
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='INTETA')
C
C
CD    CALL WLKBCD(IDPROG)
C
C***********************************************************************
      CALL ENPOUB(AM2LC,ADM2LC)
C -
CD    IF (LTRACP(1))THEN
C -
CD      CALL OMPTDP(' VALEURS REELLES AUX ANGLES '
CD                    , FTETA(1) , NTETA , 1 )
C -
CD    END IF
C utilisation de la transformation de fourier inverse qui a N
C valeur reelles associe N coeefficients du developpement en serie de
C Fourier complexe .
      INV = 1
C -
C - on recupere dans SINREE la partie reelle et dans FIMA la partie
C - complexe qui en entree est " mise "a zero
C -
      CALL GSPOUD ( NTETA , FIMA )
C -
      CALL FT01AD( NTETA ,INV , FTETA  , DM( FIMA ) )
C -
      INTF = 2.D0*PI*FTETA(1)
C -
CD    CALL IMPDN( ' VALEUR DE L''INTEGRALE EN TETA', INTF)
C***********************************************************************

      CALL SOPOUB(AM2LC,ADM2LC)
CD    CALL RETOUD(IDPROG)
      RETURN
      END
