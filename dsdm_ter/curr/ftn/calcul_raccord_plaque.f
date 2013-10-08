C     Ce programme calcul les caracteristiques plaque :
C 
C     < Kcp > ,< Kcpz > , <kcpz2 >
C 
C     Rappel:
C 
C     Caracteristiques du comportement de la
C     couche dans sa base d'orthotropie:
C 
C     |K11  K12   0   0   0  A1|  E11
C     |     K22   0   0   0  A2|  E22
C     |          K66  0   0   0|  R2*E12
C     |               B1  0   0|  R2*E23      =====>DEFORMATIONS
C     |                   B2  0|  R2*E13
C     |   SYM                 C|  E33
C 
C     On envoie
C     E............. cote1   la cote inferieure par rapporta la surface moyenne
C     E............. cote2   la cotesuperieure par rapport a la surface moyenne
C 
C     Et on recupere:        !!!!!! Dans la base globale !!!!!!
C 
C     S............. KCP     somme( cote1 ,cote2) KCP
C 
C     S............. KCPZ    somme( cote1 ,cote2) KCPz
C                                                     2
C     S............. KCPZ2   somme( cote1 ,cote2) KCPz

      SUBROUTINE CARPLA (COTE1, COTE2, KCP, KCPZ, KCPZ2)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
      DOUBLE PRECISION  COTE1, COTE2
C 
      DOUBLE PRECISION  KCP(3,3), KCPZ(3,3), KCPZ2(3,3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      DOUBLE PRECISION  TETORT, TETCAL
C 
      DOUBLE PRECISION  R2
C 
      DOUBLE PRECISION  HAUT1, HAUT2, COTEIN, EPASUP, MUL
C 
      LOGICAL APP
C 
      INTEGER  I, J, ADHAU, NCOINF, NCOSUP, NBCOCO
C 
      INTEGER  ADANG, ADCOIS, ADCOOR, DEBCOO
C 
      INTEGER  AM2LC, ADM2LC
C 
CD    LOGICAL LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CARPLA')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C 
      R2 = DSQRT(2.D0)
C 
C     On range cote1 et cote2 par ordre croissant
C 
      HAUT1 = DMIN1(COTE1,COTE2)
      HAUT2 = DMAX1(COTE1,COTE2)
      NCOINF = 0
      NCOSUP = NBCOU
      IF(SYM) HAUT1 =0.D0
C 
C 
C     Le tableau HAU-CE-TOT contient les distances au centre
C     de tous les plans en partant du plan inferieur du stratifie
C 
      CALL ADTBDM('HAU-CE-TOT', ADHAU)
C 
      DO I = 1 , NBCOU
C 
         CALL APPART( DM(ADHAU+I-1), DM(ADHAU+I), HAUT1, APP )
C 
         IF(APP) THEN
C 
CD         CALL IMPDT(' La cote inferieure '//IDPROG , HAUT1)
C 
CD         CALL IMPET(' appartient a la couche '//IDPROG , I )
C 
           NCOINF = I
           GOTO 1
C 
         END IF
C 
      END DO
C 
1     CONTINUE
C 
      DO I = NCOINF , NBCOU
C 
         CALL APPART( DM(ADHAU+I-1), DM(ADHAU+I), HAUT2, APP )
C 
         IF(APP) THEN
C 
CD         CALL IMPDT(' La cote superieure '//IDPROG , HAUT2)
C 
CD         CALL IMPET(' appartient a la couche '//IDPROG , I )
C 
           NCOSUP = I
           GOTO 2
C 
         END IF
C 
      END DO
C 
2     CONTINUE
C 
C     NORMALEMENT NCOINF est plus petit que NCOSUP
C 
CD    IF( NCOINF .GT. NCOSUP )THEN
C 
CD      CALL ERREUD ( 0 ,'NCOINF .GT. NCOSUP DANS'//IDPROG )
C 
CD    END IF
C 
C     NOMBRE DE COUCHES CONCERNEES PAR LE CALCUL
C 
      NBCOCO = NCOSUP-NCOINF+1
C 
C     Reservation de place pour les comportements concernes
C     et pour les angles calcul/base d'orthotropie concernes
C 
      CALL GSPOUD( (2*17+1) , ADANG )
C 
C     Adresse de depart du comportement contraintes planes
C     range type orthotrope
C 
      DEBCOO = ADANG + 1
      ADCOOR = DEBCOO+17
C 
C     A PRIORI ON FAIT LE CALCUL DANS LA BASE GLOBALE === >
C 
      TETCAL = 0.D0
C 
C     on fait le menage
C 
      CALL BALAID ( 9 , KCP(1,1) )
      CALL BALAID ( 9 , KCPZ(1,1) )
      CALL BALAID ( 9 , KCPZ2(1,1) )
C 
      DO I = NCOINF , NCOSUP-1
c      
        EPASUP = DM(ADHAU+I)
	COTEIN = DM(ADHAU+I-1)
C 
        CALL IMPDT('COTEIN',COTEIN)
        CALL IMPDT('EPASUP',EPASUP)
C 
        CALL IMPET( ' POUR LA COUCHE DANS '//IDPROG , I )

        CALL ANGCOU( I , DM(ADANG+I-1)  )
C 
        TETORT =  DM(ADANG+I-1)
C 
CD        CALL IMPDT ( 'VALEUR DE L''ORIENTATION DE LA BASE'
CD      &//'D''ORTHOTROPIE '//IDPROG,TETORT)
C 
C     Adresse de depart ( dans dm !! ) du comportement range type isotrope
C 
        CALL COMCOU( I , ADCOIS  )
C 
        CALL COMPLA( DM(ADCOIS  )  , TETORT , TETCAL ,
     &       DM(ADCOOR) ,DM(ADCOOR+9) ,DM(ADCOOR+13) ,DM(ADCOOR+16) )
C 
C       MUL = (h2-h1)
C 
        MUL = ( EPASUP-COTEIN )
CD      CALL IMPDT( 'MUL', MUL )
C 
        CALL MUMARE( MUL, 9, DM(ADCOOR), DM(DEBCOO) )
        CALL ADDMAD( 9 , KCP(1,1)   , DM(DEBCOO)    , KCP(1,1)    )
C 
CD         CALL IMPTDT(' KCPLOCAL   ',DM(DEBCOO), 3,3 )
CD         CALL IMPTDT(' KCP         ',KCP(1,1), 3,3 )
C 
        MUL = MUL*( EPASUP+COTEIN )/2.D0
C 
CD      CALL IMPDT( 'MULZ', MUL )
C 
        CALL MUMARE( MUL , 9 ,  DM(DEBCOO)  , DM(DEBCOO)  )
C 
CD         CALL IMPTDT(' KCPZLOCAL   ',DM(DEBCOO), 3,3 )
CD         CALL IMPTDT(' KCPZ        ',KCPZ(1,1), 3,3 )
C 
        CALL ADDMAD( 9 , KCPZ(1,1)  , DM(DEBCOO)    , KCPZ(1,1)   )
C 
C       MUL =  ( h2**3 -h1**3)/3 (ANCIENNE VERSION)
C 
        MUL = ( EPASUP*EPASUP*EPASUP-COTEIN*COTEIN*COTEIN)
     &        / ( 3.D0 )
CD      CALL IMPDT( 'MULZ2', MUL )
C 
        CALL MUMARE( MUL , 9 ,  DM(ADCOOR)  , DM(ADCOOR)  )
C 
CD         CALL IMPTDT(' KCPZ2LOCAL   ',DM(ADCOOR), 3,3 )
CD         CALL IMPTDT(' KCPZ2        ',KCPZ2(1,1), 3,3 )
C 
        CALL ADDMAD( 9 , KCPZ2(1,1) , DM(ADCOOR)    , KCPZ2(1,1)  )
C 
      END DO
C 
      I = NCOSUP
C 
      CALL IMPET( ' POUR LA COUCHE DANS '//IDPROG , I )
C 
      CALL ANGCOU( I , DM(ADANG+I-1)  )
C 
      TETORT =  DM(ADANG+I-1)
C 
CD     CALL IMPDT ( 'VALEUR DE L''ORIENTATION DE LA BASE'
CD      //'D''ORTHOTROPIE '//IDPROG,TETORT)
C 
C     Adresse de depart ( dans dm !! ) du comportement range type isotrope
C 
      CALL COMCOU( I , ADCOIS  )
C 
      CALL COMPLA( DM(ADCOIS  ), TETORT, TETCAL,
     &         DM(ADCOOR), DM(ADCOOR+9), DM(ADCOOR+13), DM(ADCOOR+16))
C 
      COTEIN = DM(ADHAU+NCOSUP-1)
      EPASUP = DM(ADHAU+NCOSUP)
C 
C     MUL = (h2-h1)
C 
      MUL = ( EPASUP-COTEIN )
C 
      CALL MUMARE( MUL , 9 ,  DM(ADCOOR)  , DM(DEBCOO)  )
      CALL ADDMAD( 9 , KCP(1,1)   , DM(DEBCOO)    , KCP(1,1)    )
C 
CD       CALL IMPTDT(' KCPLOCAL   ',DM(DEBCOO), 3,3 )
CD       CALL IMPTDT(' KCP        ',KCP(1,1), 3,3 )
C  
C     MUL = (h2+h1)/2 => mul*mulpre = ( h2**2 -h1**2)/2
C 
      MUL = MUL*( EPASUP+COTEIN )/2.D0
C 
      CALL MUMARE( MUL , 9 ,  DM(DEBCOO)  , DM(DEBCOO)  )
      CALL ADDMAD( 9 , KCPZ(1,1)  , DM(DEBCOO)    , KCPZ(1,1)   )      
C 
      CALL IMPTDT(' KCPZLOCAL   ',DM(DEBCOO), 3,3 )
      CALL IMPTDT(' KCPZ        ',KCPZ(1,1), 3,3 )
C
C     MUL =  ( h2**3 -h1**3)/3
C 
      MUL = ( EPASUP*EPASUP*EPASUP-COTEIN*COTEIN*COTEIN)
     &        / ( 3.D0 )
C 
      CALL MUMARE( MUL , 9 ,  DM(ADCOOR)  , DM(ADCOOR)  )
      CALL ADDMAD( 9 , KCPZ2(1,1) , DM(ADCOOR)    , KCPZ2(1,1)  )
C 
CD       CALL IMPTDT(' KCPZ2LOCAL   ',DM(ADCOOR), 3,3 )
CD       CALL IMPTDT(' KCPZ2       ',KCPZ2(1,1), 3,3 )
C 
      IF (SYM ) THEN
        DO I = 1 , 3
          DO J = 1 , 3
             KCPZ(I,J)  = 0.D0
             KCP (I,J)  = 2.D0*KCP(I,J)
             KCPZ2(I,J) = 2.D0*KCPZ2(I,J)
          END DO
        END DO
      END IF
C 
C     IF(LTRACP(1))THEN        ?????????????????IF INACTIF ???????????????
C 
CD      CALL IMPTDT(' KCP   ',KCP  (1,1),3,3 )
C 
CD      CALL IMPTDT(' KCPZ  ',KCPZ (1,1),3,3 )
C 
CD      CALL IMPTDT(' KCPZ2 ',KCPZ2(1,1),3,3 )
C 
C     END IF
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments:
C     E........... ISODEV   Caracteristiques du comportement de la
C                            couche dans sa base d'orthotropie
C 
C     |K11  K12   0   0   0  A1|  E11
C     |     K22   0   0   0  A2|  E22
C     |          K66  0   0   0|  R2*E12
C     |               B1  0   0|  R2*E23  =====>DEFORMATIONS
C     |                   B2  0|  R2*E13
C     |   SYM                 C|  E33
C 
C     E........... TETORT   Valeur de l'angle que fait cette base avec
C                            le repere global
C     E........... TETCAL   Valeur de l'angle pour lequel on calcule
C                            le comportement
C     Et on recupere:
C     S........... KCP(3,3) Partie du comportement contraintes planes
C                           travaillant sur la partie plane des deformations
C     S........... B(2,2)   Partie du comportement travaillant sur
C                           le cisaillement normal
C     S........... A(3)     Partie du comportement reliant
C                           la partie plane des deformations a la
C                           deformation normale
C     S........... C        Relie la contrainte a la deformation
C                           normale

      SUBROUTINE COMPLA( ISODEV , TETORT , TETCAL , KCP , B , A , C )
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C
      DOUBLE PRECISION  ISODEV(10) , TETORT , TETCAL , KCP(3,3), B(2,2)
      DOUBLE PRECISION  A(3) , C
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      DOUBLE PRECISION  R2 , COS2 , COS4 , SIN2 , SIN4 , TETLOC
C 
CD    LOGICAL LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='COMPLA')
C 
CD    CALL WLKBCD(IDPROG)
C 
C ----------------------------------------------------------------------
C     CALL IMPTDT( 'ISODEV '//IDPROG, ISODEV(1),10,1)
C 
      R2 = DSQRT(2.D0)
C 
      TETLOC = TETCAL-TETORT
C 
      COS2   = DCOS(2.D0*TETLOC)
      SIN2   = DSIN(2.D0*TETLOC)
      COS4   = DCOS(4.D0*TETLOC)
      SIN4   = DSIN(4.D0*TETLOC)
C 
C     CALCUL DE A
C 
      A(1)    = ISODEV(4)+COS2*ISODEV(8)
      A(2)    = ISODEV(4)-COS2*ISODEV(8)
      A(3)    = -SIN2*R2*ISODEV(8)
C 
C     CALL IMPTDT( 'A '//IDPROG, A(1),3,1)
C 
C     CALCUL DE C
C 
      C       = ISODEV(6)
C     CALL IMPDP('C DANS'//IDPROG,C)
C 
C 
C     CALCUL DE KCP
C 
      KCP(1,1)  = ISODEV(1)+COS2*ISODEV(7)+COS4*ISODEV(10)
     &            - A(1)*A(1)/C
      KCP(2,2)  = ISODEV(1)-COS2*ISODEV(7)+COS4*ISODEV(10)
     &            - A(2)*A(2)/C
      KCP(1,2)  = ISODEV(2)-COS4*ISODEV(10)
     &            - A(1)*A(2)/C
      KCP(2,1)  = KCP(1,2)
      KCP(1,3)  = -R2*(SIN2*ISODEV(7)/2.D0+SIN4*ISODEV(10))
     &            - A(1)*A(3)/C
      KCP(3,1)  = KCP(1,3)
      KCP(3,3)  = ISODEV(3)-2.D0*COS4*ISODEV(10)
     &            - A(3)*A(3)/C
      KCP(2,3)  = -R2*(SIN2*ISODEV(7)/2.D0-SIN4*ISODEV(10))
     &            - A(2)*A(3)/C
      KCP(3,2)  = KCP(2,3)
C 
CD    CALL IMPTDT(' KCP '//IDPROG,KCP(1,1),3,3 )
C 
C     CALCUL DE B
C 
      B(1,1)  = ISODEV(5)+COS2*ISODEV(9)
      B(2,2)  = ISODEV(5)-COS2*ISODEV(9)
      B(1,2)  = SIN2*ISODEV(9)
      B(2,1)  = B(1,2)
C 
CD    CALL IMPTDT(' B '//IDPROG,B(1,1),2,2 )
C 
CD    CALL IMPDT (
CD    'VALEUR DE L''ORIENTATION DE LA BASE D''ORTHOTROPIE   ',TETORT)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C ROUTINE QUI TOURNE KCP : KCP2 = P KCP TP, ROTATION DE TETLOC
C
      SUBROUTINE RTCPLA( TETLOC , KCP ,KCP2 )
C 
      DOUBLE PRECISION   KCP(3,3), KCP2(3,3) , TETLOC 
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      DOUBLE PRECISION  R2 , C2 , S2 , SC
      DOUBLE PRECISION  R(3,3), RT(3,3), LOC(3,3)
      INTEGER I , J , K
C 
CD    LOGICAL LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='RTCPLA')
C 
CD    CALL WLKBCD(IDPROG)
C 
C ----------------------------------------------------------------------
C 
      R2 = DSQRT(2.D0)
C 
      C2   = DCOS(TETLOC)*DCOS(TETLOC)
      S2   = DSIN(TETLOC)*DSIN(TETLOC)   
      SC     = SIN(TETLOC)*COS(TETLOC)
C ROTATION 
      R(1,1) = C2
      R(1,2) = S2
      R(1,3) = R2*SC  
      R(2,1) = S2
      R(2,2) = C2
      R(2,3) = -R2*SC
      R(3,1) = -R2*SC
      R(3,2) = R2*SC
      R(3,3) = C2-S2
C ROTATION TRANSPOSEE
      RT(1,1) = C2
      RT(1,2) = S2
      RT(1,3) = -R2*SC  
      RT(2,1) = S2
      RT(2,2) = C2
      RT(2,3) = R2*SC
      RT(3,1) = R2*SC
      RT(3,2) = -R2*SC
      RT(3,3) = C2-S2   
C
      DO I = 1, 3
        DO J = 1, 3
	   LOC ( J, I ) = 0.D0
	   KCP2 ( J, I) = 0.D0
	END DO
      END DO       
C 
C     CALCUL DE R KCP
C
      DO I = 1, 3
        DO J = 1, 3
	   DO K =1 ,3
	     LOC ( J, I ) = LOC (J ,I)+ R( J , K) *KCP (K , I)
	   END DO 
	END DO
      END DO 
C 
C     CALCUL DE R KCP
C
      DO I = 1, 3
        DO J = 1, 3
	   DO K =1 ,3
	     KCP2 ( J, I ) = KCP2 (J ,I)+ LOC( J , K) *RT (K , I)
	   END DO 
	END DO
      END DO 
C 
CD    CALL IMPTDT (' KCP      '//IDPROG, KCP(1,1), 3, 3)
CD    CALL IMPTDT (' P KCP TP '//IDPROG, KCP2(1,1), 3, 3)
C       
      RETURN
      END
