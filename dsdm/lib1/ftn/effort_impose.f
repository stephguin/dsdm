C     On envoie comme arguments :
C 
C     E........... ADEFFO  1ere adresse de MAT-EFFORT
C     E........... NUDEV   numero du developpement concerne (-NTDSFG , NTDSFG )
C     E........... NUDDL   numero des ddl concernes
C     E........... NBDDL   nombre de ces ddl
C     E........... VALEUR  valeur de ces ddl
C 
C     Assemble alors la matrice des seconds membres
C 
      SUBROUTINE ASVEFI (ADEFFO, NUDEV, NUDDL, NBDDL, VALEUR)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER          NBDDL, NUDDL(NBDDL), ADEFFO, NUDEV
      DOUBLE PRECISION VALEUR(NBDDL)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   I, ADEPPA, j
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='ASVEFI')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      ADEPPA  = ADEFFO + (NTDSFG + NUDEV)*NDDL
      DO I=1 , NBDDL
        J = ADEPPA+NUDDL(I)-1
        DM(J)  = DM(J) + VALEUR(I)
      ENDDO
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C 
C     Donnees des efforts et des deplacements :
C 
C     On caracterise les zones par une suite de 6-uplets indiquant :
C     (numero de bord, 1er numero de colonne ou de couche, dernier numero
C     de colonne ou de couche, CARACU, CARACV, CARACW).
C 
C     CONVENTION :
C 
C     CARAC = 0  <==> effort impose nul
C     CARAC = 1  <==> blocage
C     CARAC = 2  <==> blocage noeud par noeud
C     CARAC = 3  <==> blocage noeud par noeud par developpement
C 
C     A priori blocage des deplacements de solide autrement c'est 
C     un peu bizzare.
C 
C     CARAC = -1 <==> autre type d'effort impose
C     CARAC = -2 <==> autre type de deplacement impose
C     CARAC = -3 <==>  deplacement imposes noeud par noeud
C 
C     (numero de bord, 1er numero (colonne ou couche), numero de
C     noeud => 1 a 4  si nubord = 1 ou 3 !( 2,4 derivees)
C           => 1 a 6  si nubord = 2 ou 4 !( 2,4,6 derivees)
C 
      SUBROUTINE EFFIMP
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
      INTEGER  ADLOC1, ADZONE, DEBZON
      INTEGER  I, NUCO, TYPDEP, NBTDDL, NBDEV, INDDEP, DBDGJ
      INTEGER  NUBORD, NBDEGJ, INCDEP
      INTEGER  ADPDRE, ADNDEV, ADPDEV, ADCOEF, PLADEP
      INTEGER  J, INDPDG, NUDEVJ
      INTEGER  ADEFFO, ADEFIN, DDL(12), ADEPPA, EFFORT, NFT
C 
      DOUBLE PRECISION POSIC, DIMENS
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EFFIMP')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('EFFORT-INT', ADEFIN)
      CALL ADTBDM ('MAT-EFFORT', ADEFFO)
      CALL ADTBM  ('P-DEGRES  ', ADPDRE)
      CALL ADTBM  ('NUM-DEVELO', ADNDEV)
      CALL ADTBM  ('P-DEVELOPP', ADPDEV)
      CALL ADTBDM ('COEFF-POLY', ADCOEF)
C 
      CALL ADTBM  ('TLOCN1    ', ADLOC1)
C 
CD    CALL IMPEP ('ADRESSE DE DEPART DANS M DE TLOCN1 '//IDPROG, ADLOC1)
C 
      CALL ADTBM  ('ZONE-CARAC', ADZONE)
C 
C -----------------------------------------------------------------------
C     Sequence de modification des efforts pour les efforts imposes
C -----------------------------------------------------------------------
CD    CALL IMPMN ('RENTREE DANS LA SEQUENCE DE MODIF DES EFFORTS')
C 
      DEBZON   = ADZONE
      INDPDG   = ADPDRE
      INCDEP   = ADPDEV
      NBTDDL   = 0
      NFT      = NBFODO
C 
      EFFORT=ADEFFO+NDDL*NBMAT*(NFT-1)
C 
      DO I= 1, NBZONE(NFT)
C 
        NUBORD = M(DEBZON)
C 
        DO TYPDEP= 1,3
C 
          INDDEP   = INCDEP+TYPDEP
C 
          CALL IMPET ('POUR TYPDEP ', TYPDEP)
C 
          PLADEP  = DEBZON+2+TYPDEP
C 
          IF (M(PLADEP) .EQ. -2) THEN
C 
            NBDEV   = M(INDDEP)-M(INDDEP-1)
            INDPDG  = INDPDG+NBDEV
C 
CD           CALL IMPMP ('ON EST EN DEPLACEMENTS IMPOSES')
CD           CALL IMPEP ('NOMBRE DE DEVELOPPEMENT',NBDEV)
CD           CALL IMPEP ('PLACE DANS POINTEUR DE DEGRES',INDPDG)
C 
          ENDIF
C 
          IF (M(PLADEP) .EQ. -1) THEN
C 
C            CARAC = -1 <==> autre type d'effort impose
C            CALL IMPMT('ON EST EN EFFORTS IMPOSES')
C 
            NBDEV   = M(INDDEP)-M(INDDEP-1)
C 
            DO NUCO = M(DEBZON+1),M(DEBZON+2)
C 
              CALL CAGEOB (NUBORD, NUCO, POSIC, DIMENS)
C 
CD            CALL IMPEP ('NUBORD APRES CAGEOB ', NUBORD)
CD            CALL IMPEP ('NBDEV  APRES CAGEOB ', NBDEV)
C 
              IF (NUBORD .EQ. 1) THEN
C 
                CALL NDDLC (TYPDEP, 1, NUCO, ADLOC1 , DDL)
C 
              ELSE IF (NUBORD .EQ. 3) THEN
C 
                CALL NDDLC (TYPDEP, NBCOU, NUCO, ADLOC1, DDL)
C 
              ELSE IF (NUBORD .EQ. 2) THEN
C 
                CALL NDDLC (TYPDEP,  NUCO, 1, ADLOC1, DDL)
C 
              ELSE IF (NUBORD .EQ. 4) THEN
C 
                CALL NDDLC (TYPDEP, NUCO, NBCOL, ADLOC1, DDL)
C 
              ENDIF
C 
CD            CALL IMPEP ('NUBORD APRES NDDLC ', NUBORD)
CD            CALL IMPEP ('NBDEV  APRES NDDLC ', NBDEV)
C 
              DO J =1, NBDEV
C 
                DBDGJ    = M(INDPDG+J-1)
                NUDEVJ   = M(ADNDEV+M(INDDEP-1)+J-1)
                ADEPPA   = EFFORT+NDDL*(NTDSFG+NUDEVJ)
                NBDEGJ   = M(INDPDG+J)-DBDGJ
C 
                CALL IMPET ('POUR LE DEVELOPPEMENT ', NUDEVJ)
C 
CD              CALL IMPEP( 'NUBORD AVANT CALEFF', NUBORD)
CD              CALL IMPDP(' VALEUR AU CENTRE ', POSIC )
CD              CALL IMPDP(' DEMIE DIMENSION  ', DIMENS )
C 
                CALL IMPTET ('NUM-DDDL-CROI ', DDL(1), 1, 12)
C 
                CALL CALEFF (ADEPPA-1, DM(ADEFIN), NUBORD, NBDEGJ,
     &                       DM(ADCOEF+DBDGJ), DDL, POSIC, DIMENS)
C 
              ENDDO
C 
            ENDDO
C 
            INDPDG = INDPDG+NBDEV
C 
          ENDIF
C 
        ENDDO
C 
        DEBZON   = DEBZON+6
        INCDEP   = INCDEP+3
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
C     QUE FAIT CETTE ROUTINE :
C 
C     On calcule les tableaux neccessaires au calcul des
C     termes d'efforts (MAT-EFFORT).
C     Attention les tableaux ne sont pas multiplies par MATLOC,
C     ils doivent l'etre dans les routines ASSBFI.
C 
      SUBROUTINE EFFINT
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  H1, H2, H3, H4, X
      DOUBLE PRECISION  W1, W2, W3, W4
      DOUBLE PRECISION  L1, L2, L3, Y
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I, ADL1, ADR1, H, K
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EFFINT')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     Creation du tableau des termes integres pour les EFFORTS :
C 
C            - nom: EFFORT-INT
C            - 1ere adresse libre: ADL1
C            - numero :N1
C 
C     le 1er tableau est le tableau integre pour les interfaces
C     le 2eme tableau est le tableau integre*x pour les interfaces
C -----------------------------------------------------------------------
        CALL GESTDP('EFFORT-INT',25,ADL1)
        CALL MENADM(ADL1,25)
        ADR1=ADL1-1
C 
C -----------------------------------------------------------------------
C     calcul du tableau des efforts pour les bords inferieurs et superieurs
C     puis interieur
C 
C     Le principe de rangement est le suivant :
C 
C     Pour les bords inferieurs et superieurs :
C 
C     D'abord H1, H1*X, H1*X*X, H1*X*X*X
C        puis H3, H3*X, H3*X*X, H3*X*X*X
C        puis H2, H2*X, H2*X*X, H2*X*X*X
C        puis H4, H4*X, H4*X*X, H4*X*X*X
C 
C     Pour le bord interieur :
C 
C     D'abord L1, L1*Y, L1*Y*Y
C        puis L2, L2*Y, L2*Y*Y
C        puis L3, L3*Y, L3*Y*Y
C -----------------------------------------------------------------------
C 
C     POUR LES BORDS INFERIEURS ET SUPERIEURS
C 
      K=XINTEG*(XINTEG-1)/2
      DO I=1 , XINTEG
         X            = GAUSS(K+I)
         W1           = POIDS(K+I)
         W2           = W1*X
         W3           = W2*X
         W4           = W3*X
C 
C        Pour le noeud gauche pour le deplacement
C 
         DM(ADR1+1)   = DM(ADR1+1)  + W1*H1(X)
         DM(ADR1+2)   = DM(ADR1+2)  + W2*H1(X)
         DM(ADR1+3)   = DM(ADR1+3)  + W3*H1(X)
         DM(ADR1+4)   = DM(ADR1+4)  + W4*H1(X)
C 
C        Pour le noeud droit pour le deplacement
C 
         DM(ADR1+5)   = DM(ADR1+5)  + W1*H3(X)
         DM(ADR1+6)   = DM(ADR1+6)  + W2*H3(X)
         DM(ADR1+7)   = DM(ADR1+7)  + W3*H3(X)
         DM(ADR1+8)   = DM(ADR1+8)  + W4*H3(X)
C 
C        Pour le noeud gauche pour la derivee du  deplacement
C 
         DM(ADR1+9)   = DM(ADR1+9)  + W1*H2(X)
         DM(ADR1+10)  = DM(ADR1+10) + W2*H2(X)
         DM(ADR1+11)  = DM(ADR1+11) + W3*H2(X)
         DM(ADR1+12)  = DM(ADR1+12) + W4*H2(X)
C 
C        Pour le noeud droit pour la derivee du  deplacement
C 
         DM(ADR1+13)  = DM(ADR1+13) + W1*H4(X)
         DM(ADR1+14)  = DM(ADR1+14) + W2*H4(X)
         DM(ADR1+15)  = DM(ADR1+15) + W3*H4(X)
         DM(ADR1+16)  = DM(ADR1+16) + W4*H4(X)
C 
       ENDDO
C 
C     POUR LE BORD INTERIEUR
C 
      ADR1=ADR1+16
      H=YINTEG*(YINTEG-1)/2
      DO I=1 , YINTEG
         W1           = POIDS(H+I)
         Y            = GAUSS(H+I)
         W2           = W1*Y
         W3           = W2*Y
C 
C        Pour le noeud inferieur
C 
         DM(ADR1+1)   = DM(ADR1+1)  + W1*L1(Y)
         DM(ADR1+2)   = DM(ADR1+2)  + W2*L1(Y)
         DM(ADR1+3)   = DM(ADR1+3)  + W3*L1(Y)
C 
C        Pour le noeud milieu
C 
         DM(ADR1+4)   = DM(ADR1+4)  + W1*L3(Y)
         DM(ADR1+5)   = DM(ADR1+5)  + W2*L3(Y)
         DM(ADR1+6)   = DM(ADR1+6)  + W3*L3(Y)
C 
C        Pour le noeud superieur
C 
         DM(ADR1+7)   = DM(ADR1+7)  + W1*L2(Y)
         DM(ADR1+8)   = DM(ADR1+8)  + W2*L2(Y)
         DM(ADR1+9)   = DM(ADR1+9)  + W3*L2(Y)
       ENDDO
C 
CD    CALL IMPEP ('POUR XINTEG ',XINTEG)
CD    CALL IMPTDP(' Efforts integres pour les fonctions HI dans:'
CD                //IDPROG ,DM(ADL1),1 ,16)
CD    CALL IMPEP ('POUR YINTEG ',YINTEG)
CD    CALL IMPTDP(' Efforts integres pour les fonctions LI dans:'
CD                 //IDPROG,DM(ADR1+1),1 ,9)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine remplit les tableaux des efforts
C 
C     On envoie comme arguments :
C 
C     E ...... ADDEPA  Adresse de depart -1 du tableau des efforts concernes
C                      ========> resultat de CALL ADEFF0
C                               (Est calculee en fonction de la
C                                valeur du rang du developpement en series
C                                et donc du type d'effort (1 ou 2))
C     E ...... EFFINT  Tableau des efforts integres
C     E ...... NUBORD  numero du bord concerne si 1 ====> bord inferieur
C                                              si 2 ====> bord interieur
C                                              si 3 ====> bord superieur
C     E ...... NBCONS=DEGRES+1  le nombre de constantes definissant
C                               l'effort concerne
C     E ...... VALEFF(NBCONS)   la valeur de ces efforts
C 
C     E ...... NUDDL            valeur des numero de ddl ranges croissants
C                               correspondant au deplacement
C     E ...... RC               position du centre de l'element :
C                               rayon ou hauteur
C     E ...... A                demi-longueur de l'element
C                               demi-longueur ou epaisseur
C 
      SUBROUTINE CALEFF(ADDEPA, EFFINT, NUBORD, NBCONS,
     $     VALEFF, NUDDL, RC, A)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NUBORD , NBCONS , NUDDL(12) , ADDEPA
      DOUBLE PRECISION  VALEFF(NBCONS) , EFFINT(25) , RC , A
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER    TEST , NDDLCO , EFFLOC , DDL(6) , ADLOC , NBTINT
      INTEGER    NTCAL , DECAL , I , J
C 
      DOUBLE PRECISION   C(8)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALEFF')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (NBCONS.EQ.0)THEN
        GOTO 1
      ENDIF
      IF(NUBORD.EQ.1.OR.NUBORD.EQ.3)THEN
        TEST=NBCOU
      ELSE IF(NUBORD.EQ.2)THEN
        TEST=NBCOL
      ELSE
        CALL ERREUD(0,
     $  ' MAUVAIS INDICATEUR DE BORD DANS :'//IDPROG)
      ENDIF
C -----------------------------------------------------------------------
C     Correspondance entre la numerotation elementaire et la numerotation
C     reelle croissantE issue des routines NDDL(U,V,W)C.
C     WARNING !!!!! La correspondance est faite a partir des bords
C     et les ddl correspondant aux deplacements sont ranges avant.
C 
C     NTCAL est le nombre de termes de calcul pour un degre donne.
C     EFFLOC indique a partir d'ou il faut aller chercher les termes
C     dans EFFORT-INT.
C     NBTINT indique le nombre de termes correspondant a un deplacement
C     dans EFFORT-INT
C -----------------------------------------------------------------------
CD    CALL IMPTDP('VALEUR DES EFFORTS INTEGRES DANS: '//IDPROG,
CD    EFFINT(1),1,25)
CD    CALL IMPDP ('VALEUR DE LA LONGUEUR DANS :'//IDPROG,A)
CD    CALL IMPDP ('VALEUR AU CENTRE DE L''ELEMENT DANS :'//IDPROG,RC)
CD    CALL IMPEP('TYPE DE BORD ',NUBORD)
CD    CALL IMPEP ('VALEUR DU DEGRES+1 ',NBCONS)
CD     DO I=1,NBCONS
CD       CALL IMPDP('VALEUR DES EFFORTS ',VALEFF(I))
CD     ENDDO
C 
C -----------------------------------------------------------------------
      IF(NUM.EQ.1)THEN
        IF(NUBORD.EQ.1) THEN
          DDL(1)        = NUDDL( 1 )
          DDL(2)        = NUDDL( 3 )
          DDL(3)        = NUDDL( 2 )
          DDL(4)        = NUDDL( 4 )
          NDDLCO        = 2
          EFFLOC        = 0
          TEST          = NBCOU
          NTCAL         = NBCONS+1
          NBTINT        = 4
        ELSE IF (NUBORD.EQ.3) THEN
          DDL(1)        = NUDDL( 9 )
          DDL(2)        = NUDDL( 11 )
          DDL(3)        = NUDDL( 10 )
          DDL(4)        = NUDDL( 12 )
          NDDLCO        = 2
          EFFLOC        = 0
          TEST          = NBCOU
          NTCAL         = NBCONS+1
          NBTINT        = 4
        ELSE IF (NUBORD.EQ .2) THEN
          DDL(1)        = NUDDL( 1 )
          DDL(2)        = NUDDL( 5 )
          DDL(3)        = NUDDL( 9 )
          NDDLCO        = 3
          TEST          = NBCOL
          EFFLOC        = 16
          NTCAL         = NBCONS
          NBTINT        = 3
        ELSE
          CALL ERREUD(0,
     $    ' MAUVAIS INDICATEUR DE BORD DANS :'//IDPROG)
        ENDIF
      ELSE IF(NUM.EQ.2) THEN
        IF(NUBORD.EQ.1) THEN
          DDL(1)        = NUDDL( 1 )
          DDL(2)        = NUDDL( 7 )
          DDL(3)        = NUDDL( 2 )
          DDL(4)        = NUDDL( 8 )
          NDDLCO        = 2
          EFFLOC        = 0
          TEST          = NBCOU
          NTCAL         = NBCONS+1
          NBTINT        = 4
        ELSE IF(NUBORD.EQ.3) THEN
          DDL(1)        = NUDDL( 5 )
          DDL(2)        = NUDDL( 11 )
          DDL(3)        = NUDDL( 6 )
          DDL(4)        = NUDDL( 12 )
          NDDLCO        = 2
          EFFLOC        = 0
          TEST          = NBCOU
          NTCAL         = NBCONS+1
          NBTINT        = 4
        ELSE IF(NUBORD.EQ.2) THEN
          DDL(1)        = NUDDL( 1 )
          DDL(2)        = NUDDL( 3 )
          DDL(3)        = NUDDL( 5 )
          NDDLCO        = 3
          TEST          = NBCOL
          EFFLOC        = 16
          NTCAL         = NBCONS
          NBTINT        = 3
        ELSE
          CALL ERREUD(0,
     $    'VALEUR ILLICITE DE NUM DANS : '//IDPROG)
        ENDIF
      ENDIF
C -----------------------------------------------------------------------
C     Pour les bords 1 et 3 :
C 
C     Les derniers termes de C sont multiplies par A
C     car ils correspondent a H2 et H4 (donc aux derivees du
C     deplacement ).
C 
C     Pour le bord 2 :
C 
C     Seuls les deplacements sont concernes
C -----------------------------------------------------------------------
C 
C     SEQUENCE DE CALCUL DES EFFORTS 
C 
C -----------------------------------------------------------------------
C 
C     DEBUT DE LA SEQUENCE POUR LES BORDS INF ET SUP
C 
C -----------------------------------------------------------------------
      IF(NUBORD.EQ.1.OR.NUBORD.EQ.3) THEN
        IF (NBCONS.EQ.1) THEN
          C(1)        = A*RC*VALEFF(1)
          C(2)        = A*A*VALEFF(1)
          C(5)        = A*C(1)
          C(6)        = A*C(2)
        ELSE IF (NBCONS.EQ.2) THEN
          C(1)        = A*RC*(VALEFF(1)+RC*VALEFF(2))
          C(2)        = A*A*(VALEFF(1)+RC*2*VALEFF(2))
          C(3)        = A*A*A*VALEFF(2)
          C(5)        = A*C(1)
          C(6)        = A*C(2)
          C(7)        = A*C(3)
        ELSE IF (NBCONS.EQ.3) THEN
          C(1)        = A*RC*(VALEFF(1)+RC*(VALEFF(2)+RC*VALEFF(3)))
          C(2)        = A*A*(VALEFF(1)+RC*(2*VALEFF(2)+3*RC*VALEFF(3)))
          C(3)        = A*A*A*(VALEFF(2)+3*RC*VALEFF(3))
          C(4)        = A*(RC**3)*VALEFF(3)
          C(5)        = A*C(1)
          C(6)        = A*C(2)
          C(7)        = A*C(3)
          C(8)        = A*C(4)
        ELSE
         CALL ERREUD(0,
     $  ' MAUVAIS PASSAGE DU NBCONS DANS :'//IDPROG)
        ENDIF
CD      CALL IMPTDP('VALEUR DES CONSTANTES C(I) DANS: '//IDPROG,
CD      C(1),1,8)
C 
C     POUR LES DDL DU DEPLACEMENT
C 
        DO I= 1 , NDDLCO
          ADLOC    = ADDEPA+DDL(I)
          DECAL  = EFFLOC+NBTINT*(I-1)
          DO J= 1 , NTCAL
            DM(ADLOC)= DM(ADLOC)+C(J)*EFFINT(DECAL+J)
          ENDDO
        ENDDO
C 
C     POUR LES DDL DES DERIVEES DU DEPLACEMENT
C 
        DO I= NDDLCO+1 , 2*NDDLCO
          ADLOC    = ADDEPA+DDL(I)
CD        CALL IMPEP('POUR LE DDL NUMERO DANS '//IDPROG,DDL(I))
          DECAL  = EFFLOC+NBTINT*(I-1)
          DO J= 1 , NTCAL
            DM(ADLOC)= DM(ADLOC)+C(J+NBTINT)*EFFINT(DECAL+J)
          ENDDO
        ENDDO
      ELSE
C 
C -----------------------------------------------------------------------
C 
C     DEBUT  DE LA SEQUENCE POUR LE BORD INTERIEUR
C 
C -----------------------------------------------------------------------
        IF (NBCONS.EQ.1) THEN
          C(1)        = A*VALEFF(1)
        ELSE IF (NBCONS.EQ.2) THEN
          C(1)        = A*(VALEFF(1)+RC*VALEFF(2))
          C(2)        = A*A*VALEFF(2)
        ELSE IF (NBCONS.EQ.3) THEN
          C(1)        = A*(VALEFF(1)+RC*(VALEFF(2)+RC*VALEFF(3)))
          C(2)        = A*A*(VALEFF(2)+2*RC*VALEFF(3))
          C(3)        = A*A*A*VALEFF(3)
        ELSE
         CALL ERREUD(0,
     $  ' MAUVAIS PASSAGE DU NBCONS DANS :'//IDPROG)
        ENDIF
        DO I= 1 , NDDLCO
          ADLOC    = ADDEPA+DDL(I)
          DECAL  = EFFLOC+NBTINT*(I-1)
          DO J= 1 , NTCAL
            DM(ADLOC)= DM(ADLOC)+C(J)*EFFINT(DECAL+J)
          ENDDO
        ENDDO
      ENDIF
C 
  1   CONTINUE
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
