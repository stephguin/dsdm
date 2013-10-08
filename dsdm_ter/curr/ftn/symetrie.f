C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine assemble pour un element donne les termes
C     provenant de CALKI pour toutes les matrices KOn
C 
C     On envoie comme arguments :
C 
C     ASSEMBLAGE PARTICULIER DU AU CAS DE SYMETRIE SUR
C     L'INTERFACE NUM 1
C 
C     E ...... N      valeur du rang de Fourier
C     E ...... TK0    tableau des adresses de depart des matrices KON
C     E ...... UCAL   tableau des adresses elementaires UCAL
C     E ...... TAB1   tableau du terme independant de n provenant de
C                     CALKI1 pour l'element considere (le meme pour U et V)
C     E ...... TAB3   tableau du terme independant de n provenant de
C                     CALKI3 pour l'element considere
C     E ...... NDDLU  tableau des numeros reels des ddl du deplacement
C                     U de l'element d'interface considere
C     E ...... NDDLV  tableau des numeros reels des ddl du deplacement
C                     V de l'element d'interface considere
C     E ...... NDDLW  tableau des numeros reels des ddl du deplacement
C                     U de l'element d'interface considere
C     E ...... P1     1ere adresse du tableau profil
C 
      SUBROUTINE ASSBIP (N, TK0N, UCAL, TAB1, TAB3, NDDLU, NDDLV,
     &                   NDDLW, P1)
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
      INTEGER UCAL(8), NDDLU(8), P1, NUI
      INTEGER UCALI, N, I, J, NUJ, ADU1
      INTEGER NDDLV(8), NDDLW(8), NVI, NWI
      INTEGER NVJ, NWJ, ADV1, ADW1
      INTEGER TK0N, PROFU(8), PROFV(8), PROFW(8)
      DOUBLE PRECISION  TAB1(8,8), TAB3(8,8), TCAL1, TCAL3
C 
      CHARACTER*6 IDPROG, BARATI
      PARAMETER (IDPROG='ASSBIP')
C 
CD    CALL WLKBCD(IDPROG)
C 
CD    CALL IMPTDP('VALEUR DE TAB1 DANS '//IDPROG,TAB1(1,1),8,8)
CD    CALL IMPTDP('VALEUR DE TAB3 DANS '//IDPROG,TAB3(1,1),8,8)
CD    CALL IMPTEN('RANG U CROISSANT '//IDPROG,NDDLU(1),1,8)
CD    CALL IMPTEN('RANG V CROISSANT '//IDPROG,NDDLV(1),1,8)
CD    CALL IMPTEN('RANG W CROISSANT '//IDPROG,NDDLW(1),1,8)
CD    CALL IMPTEN('CORRESPONDANCE CALCUL'//IDPROG,UCAL(1),1,8)
CD    CALL IMPEN('VALEUR DE NTSFG DANS ASSBI',NTDSFG)
C 
      DO I=5,8
         PROFU(I)    =M(P1+NDDLU(I))
         PROFV(I)    =M(P1+NDDLV(I))
         PROFW(I)    =M(P1+NDDLW(I))
      ENDDO
C 
CD    CALL IMPTEN('PROF DE U',PROFU(1),1,8)
CD    CALL IMPTEN('PROF DE V',PROFV(1),1,8)
CD    CALL IMPTEN('PROF DE W',PROFW(1),1,8)
C 
        DO I=5,8
          NUI    = NDDLU(I)
          NVI    = NDDLV(I)
          NWI    = NDDLW(I)
          UCALI = UCAL(I)
          DO J=I,8
C 
CD          CALL IMPEP('POUR J',J)
C 
            NUJ=NDDLU(J)
            NVJ=NDDLV(J)
            NWJ=NDDLW(J)
            ADU1=TK0N-1+NUI-NUJ+M(P1+NUJ)
            ADV1=TK0N-1+NVI-NVJ+M(P1+NVJ)
            ADW1=TK0N-1+NWI-NWJ+M(P1+NWJ)
            TCAL1=TAB1(UCALI,UCAL(J))
            TCAL3=TAB3(UCALI,UCAL(J))
            DM(ADU1)=DM(ADU1)+TCAL1
            DM(ADV1)=DM(ADV1)+TCAL1
            DM(ADW1)=DM(ADW1)+TCAL3
          ENDDO
        ENDDO
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
        WRITE(BARATI(1:6),'(I6)')N
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Pour les elements de type interface dans le cas particulier
C     d'une symetrie a l'interface centarle.
C 
C     Cette routine calcule pour tous les termes du developpement
C     en serie de Fourier pour le point de Gauss correspondant a
C     TABGAUS les deplacements elementaires de U, V, W  qui
C     correspondent a NDDLU, NDDLV, NDDLW; ceux-ci  sont directement
C     assembles dans DM a partir de ADSMEG de la facon suivante :
C     (NFTGLO, NDDL, NBMAT) ==> on peut directement utiliser la
C     routine de descente remontee pour plusieurs seconds membres.
C     Les seconds membres correspondant au meme rang de Fourier pour
C     des fonctions du temps differentes etant rangees sequentiellement
C 
C     On envoie comme arguments :
C 
C     E ......  POIDG   Le poids au point de gauss de l'element
C     E ......  NFT     le numero de la fonction du temps
C     E ...... TABGAU   les tableaux N du point de gauss
C     E ...... A        La demi-longeur de l'element
C     E ......  R       Le rayon au point de gauss de l'element
C     E ...... NUDDL    numero des ddl du deplacement ranges
C     E                 calcul pour l'element auquel appartient
C     E                 le point de GAUSS c'est a dire NDDLU, NDDLV, NDDLW calcul
C     E ...... SINDEV   developpement serie de Fourier reel des
C     E                 contraintes normales integrees sur le temps
C     E                 rangees composantes des contraintes
C     E                 apres composantes (6, nbmat)
C     E ...... ADSMEG   l'adresse de depart des seconds membres
C     E                 pour l'etape globale ceux-ci sont
C     E                 assembles dans DM a partir de ADSMEG de
C     E                 la facon suivante (NDDL, NFTGLO, NBMAT)
C 
      SUBROUTINE SMITEP (POIDG, NFT, TABGAU, A, R, NUDDL,
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
      PARAMETER (IDPROG='SMITEP')
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
      DO I = 1 , NBMAT
C 
        NUDEV = I-NTDSFG-1
C 
CD      CALL IMPEN( ' Pour le developpement ' , NUDEV )
C 
        CALL BNSIGN( POIDG , TABGAU , NUDEV , A , R , SINDEV( DEBSIG) ,
     &               DEPLA(1) , DEPLA(9) , DEPLA(17) )

        CALL DDLPAR ( NUDDL ,DEPLA)
C 
C     Le decalage du a NUDEV dans ASVEFI est de (NUDEV+NTDSFG)*NDDL.
C     On le met a zero en imposant NUDEV = NTDSFG et on impose le bon
C     decalage a l'aide de ADEFFO.
C 
        CALL ASVEFI( ADEFFO, -NTDSFG , NUDDL , 12 , DEPLA(1) )
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
C     Pour les elements de type interface dans le cas particulier
C     d'une symetrie a l'interface centrale.
C 
C     On envoie comme arguments :
C 
C     ES ...... NUDDL   le poids au point de gauss de l'element
C     ES ...... DEPLA   le numero de la fonction du temps

      SUBROUTINE DDLPAR (NUDDL, DEPLA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER   NUDDL(24)
C 
      DOUBLE PRECISION DEPLA(24)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DDLPAR')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
C     Pour n'assembler que 12 valeurs et non pas 24
C 
      IF ( NUM . EQ. 1) THEN
C 
C     Pour U
C 
        NUDDL(1)  = NUDDL(3)
        NUDDL(2)  = NUDDL(4)
        NUDDL(3)  = NUDDL(7)
        NUDDL(4)  = NUDDL(8)
C 
C     Pour v
C 
        NUDDL(5)  = NUDDL(11)
        NUDDL(6)  = NUDDL(12)
        NUDDL(7)  = NUDDL(15)
        NUDDL(8)  = NUDDL(16)
C 
C     Pour w
C 
        NUDDL(9)  = NUDDL(19)
        NUDDL(10) = NUDDL(20)
        NUDDL(11) = NUDDL(23)
        NUDDL(12) = NUDDL(24)
C 
CD     DO I = 1, 12
C 
CD       IF ( NUDDL(I) .LT.1 ) THEN
C 
CD         CALL ERREUD ( 0, ' Mauvais numerO de '//
CD          'ddl dans '//IDPROG )
C 
CD       END IF
C 
CD    END DO
C 
C 
C     Pour U
C 
        DEPLA(1)  = DEPLA(3)
        DEPLA(2)  = DEPLA(4)
        DEPLA(3)  = DEPLA(7)
        DEPLA(4)  = DEPLA(8)
C 
C     Pour v
C 
        DEPLA(5)  = DEPLA(11)
        DEPLA(6)  = DEPLA(12)
        DEPLA(7)  = DEPLA(15)
        DEPLA(8)  = DEPLA(16)
C 
C     Pour w
C 
        DEPLA(9)  = DEPLA(19)
        DEPLA(10) = DEPLA(20)
        DEPLA(11) = DEPLA(23)
        DEPLA(12) = DEPLA(24)
C 
      ELSE
C 
C     Pour U
C 
        NUDDL(1)  = NUDDL(5)
        NUDDL(2)  = NUDDL(6)
        NUDDL(3)  = NUDDL(7)
        NUDDL(4)  = NUDDL(8)
C 
C     Pour v
C 
        NUDDL(5)  = NUDDL(13)
        NUDDL(6)  = NUDDL(14)
        NUDDL(7)  = NUDDL(15)
        NUDDL(8)  = NUDDL(16)
C 
C     Pour w
C 
        NUDDL(9)  = NUDDL(21)
        NUDDL(10) = NUDDL(22)
        NUDDL(11) = NUDDL(23)
        NUDDL(12) = NUDDL(24)
C 
CD     DO I = 1, 12
C 
CD       IF ( NUDDL(I) .LT.1 ) THEN
C 
CD         CALL ERREUD ( 0, ' Mauvais numerO de '//
CD          'ddl dans '//IDPROG )
C 
CD       END IF
C 
CD    END DO
C 
C 
C     Pour U
C 
        DEPLA(1)  = DEPLA(5)
        DEPLA(2)  = DEPLA(6)
        DEPLA(3)  = DEPLA(7)
        DEPLA(4)  = DEPLA(8)
C 
C     Pour v
C 
        DEPLA(5)  = DEPLA(13)
        DEPLA(6)  = DEPLA(14)
        DEPLA(7)  = DEPLA(15)
        DEPLA(8)  = DEPLA(16)
C 
C     Pour w
C 
        DEPLA(9)  = DEPLA(19)
        DEPLA(10) = DEPLA(20)
        DEPLA(11) = DEPLA(23)
        DEPLA(12) = DEPLA(24)
C 
      END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cas particulier de l'interface centrale calcul symetrique
C 
C     On envoie comme arguments :
C 
C     E ...... TOUDEP   Tableau MAT-DEPLA (DM(ADEPPA))
C     E ...... RANINT   Tableau RANINV-INT (M(ARINCO))
C     E ...... NDDLU    Numero des deplacements de U ranges croissant
C     E ...... NDDLV    Numero des deplacements de V ranges croissant
C     E ...... NDDLW    Numero des deplacements de W ranges croissant
C 
C     Et on recupere :
C 
C     S ...... DEPCAN   Tableau des valeurs des deplacements pour
C                       l'element ranges calcul par ordre
C                       croissant de numero de developpement
C                       c-a-D -NTDSFG ------- > NTDSFG

         SUBROUTINE SADVCP (TOUDEP, RANINT, NDDLU, NDDLV, NDDLW
     &                      , DEPCAN )
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  TOUDEP(NBMAT*NDDL) , DEPCAN(24*NBMAT)
      INTEGER           RANINT(8) , NDDLU(8) , NDDLV(8) , NDDLW(8)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SADVCP')
      INTEGER     DECALN , K , U , V , W  , NUDEV
C  
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      DECALN = 0
      K      = 1
C 
CD    CALL IMPEP(' NDDL DANS '//IDPROG,NDDL)
C 
      DO NUDEV  =  1 , NBMAT
C 
CD       CALL IMPEP(' POUR NUDEV ',NUDEV-NTDSFG-1)
C 
         DO U =  1 , 8
                DEPCAN(K)   = 0.D0
                K           = K+1
         ENDDO
C 
         DO V =  1 , 8
                DEPCAN(K)   = 0.D0
                K           = K+1
         ENDDO
C 
         IF ( NUM .EQ. 1 ) THEN
C 
           DO W =  1 , 2
                DEPCAN(K)   = 0.D0
                K           = K+1
           ENDDO
C 
           DO W =  3 , 4
C 
C         CALL IMPEP('POWR W ', W)
C         CALL IMPEP('NDDLW(RANINT(W)) ', NDDLW(RANINT(W)))
C 
                DEPCAN(K)   = TOUDEP( DECALN+NDDLW ( RANINT( W ) ))
                K           = K+1
           ENDDO

           DO W =  5 , 6
                DEPCAN(K)   = 0.D0
                K           = K+1
           ENDDO
C 
           DO W =  7 , 8
C 
C         CALL IMPEP('POWR W ', W)
C         CALL IMPEP('NDDLW(RANINT(W)) ', NDDLW(RANINT(W)))
C 
                DEPCAN(K)   = TOUDEP( DECALN+NDDLW ( RANINT( W ) ))
                K           = K+1
           ENDDO
C 
         ELSE
C 
           DEPCAN(K)   = 0.D0
           K           = K+1
           DEPCAN(K)   = 0.D0
           K           = K+1
           DEPCAN(K)   = TOUDEP( DECALN+ NDDLW (5) )
           K           = K+1
           DEPCAN(K)   = TOUDEP( DECALN+ NDDLW (7) )
           K           = K+1
           DEPCAN(K)   = 0.D0
           K           = K+1
           DEPCAN(K)   = 0.D0
           K           = K+1
           DEPCAN(K)   = TOUDEP( DECALN+ NDDLW (6) )
           K           = K+1
           DEPCAN(K)   = TOUDEP( DECALN+ NDDLW (8) )
           K           = K+1
C 
         END IF
C 
         DECALN = DECALN+NDDL
C 
      ENDDO
C 
CD    CALL OMPTDN ('DEPLACEMENTS RANGES CALCULS PAR DEVELOPPEMENT ',
CD                  DEPCAN(1) , 24 , NBMAT )
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On cree le tableau des numeros de matrices K0N a calculer
C     en fonction des symetries
C 
C     sym/x : n<0 non calculees
C     sym/y : n>=0 et n impaires non calculees et n<0 et n paires non calculees
C     sym/o : n impaires non calculees

      SUBROUTINE TRASYM
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
      INTEGER I , MOD , NMAT , NONMAT , J
      INTEGER DEBNUM , ADNUM , ADEPPA , ADEPP1
      INTEGER DECAL , NMATP
      INTEGER AM2LC, ADM2LC
C 
      LOGICAL TEST
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TRASYM')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
      NMAT = 0
      CALL POUSME( NBMAT , ADNUM)
      DEBNUM = ADNUM
C 
      IF (SYMX .OR . ( SYMY .AND. SYMO) ) GOTO 1
C 
      DO I = -NTDSFG , -1
C 
         IF (SYMY .AND. MOD(I,2).EQ.0) GOTO 11
         IF (SYMO .AND. MOD(I,2).NE.0) GOTO 11
         M(DEBNUM) = I
         DEBNUM = DEBNUM+1
         NMAT   = NMAT+1
C 
11       CONTINUE
      END DO
C 
1     CONTINUE
C 
      DO I = 0 , NTDSFG
C 
         IF (  ( SYMY .OR. SYMO) .AND. MOD(I,2).NE.0) GOTO 12
         M(DEBNUM) = I
         DEBNUM = DEBNUM+1
         NMAT   = NMAT+1
C 
12       CONTINUE
      END DO
C 
2     CONTINUE
C 
      DECAL = 0
      NMATP = NMAT
C 
C       IF ( M(ADNUM) .EQ . -NTDSFG )  THEN
C          NMAT  = NMAT-1
C          DECAL = 1
C       END IF
C 
C       IF ( M(ADNUM+NMATP-1) .EQ . NTDSFG)  NMAT = NMAT-1
C 
      CALL GESTEN ('NUMDEV-MAT', NMAT , ADEPPA )
C 
      CALL COPITE (NMAT, M(ADNUM+DECAL), M(ADEPPA))
C 
      CALL GESTEN ('NUMERO-MAT', NMAT, ADEPP1)
C 
      DO I = 1 , NMAT
       M(ADEPP1+I-1) = M(ADEPPA+I-1)+NTDSFG+1
      END DO
C 
      CALL IMPTET ('NUMERO DES DEV DES MATRICES STOCKEES '//IDPROG,
     &              M(ADEPPA), 1, NMAT)
C 
      CALL GESTEN ('MAT-NULLES', NBMAT-NMAT, ADEPPA)
C 
      NONMAT = 0
C 
      DO I = 1 , NBMAT
C 
         TEST = .TRUE.
         DO J =1 ,NMAT
           IF( I .EQ. M(ADEPP1+J-1) ) TEST = .FALSE.
         END DO
         IF(TEST) THEN
            M(ADEPPA+NONMAT) = I
            NONMAT           = NONMAT+1
         END IF
      END DO
C 
      IF( NONMAT .NE. NBMAT-NMAT) 
     &  CALL ERREUD (0, 'Incoherence entre matrice nulle et non nulle '
     &               //IDPROG)
C 
      CALL IMPTET ('NUMERO DES MATRICES STOCKEES '//IDPROG,
     &              M(ADEPP1), 1, NMAT)
      CALL IMPTET ('NUMERO DES MATRICES NULLES '//IDPROG,
     &              M(ADEPPA), 1, NONMAT)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
C      END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Renvoie les numeros d'enregistremet des matrices
C     correspondant au numero de developpemet TNUDEV
C 
      SUBROUTINE NUMACA (NBDEV, TNUDEV, TNUENR)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NBDEV, TNUDEV(NBDEV), TNUENR(NBDEV)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER I, NMAT, J
      INTEGER ADEPPA
      INTEGER AM2LC, ADM2LC
C 
      LOGICAL  PROB
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NUMACA')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
CD    CALL IMPTET ('NUMERO DES EFFORTS NON NULS '//IDPROG,
CD                   TNUDEV, 1, NBDEV)
C 
      CALL INFOEN ('NUMERO-MAT', ADEPPA, NMAT)
CD    CALL IMPTET ('NUMERO DES MATRICES         '//IDPROG,
CD                  M(ADEPPA), 1, NMAT)

C 
      IF (NMAT .LT. NBDEV) THEN
       CALL IMPET ('NBDEV DANS '//IDPROG, NBDEV)
       CALL IMPET ('NMAT DNAS '//IDPROG, NMAT)
       CALL ERREUD (0, 'Pas assez de matrice dans        ' //IDPROG)
      END IF

C 
      DO I = 1, NBDEV
C 
         PROB = .TRUE.
         DO J = 1 , NMAT
           IF (TNUDEV(I) .EQ. M(ADEPPA+J -1)) PROB = .FALSE.
         END DO
C 
        IF ( PROB )
     &  CALL ERREUD (0 , 'probleme de coherence dans ' //IDPROG)
C 
      END DO
C 
      DO I = 1, NBDEV
C 
         DO J = 1 , NMAT
           IF (TNUDEV(I) .EQ. M(ADEPPA+J-1)) TNUENR(I) = J
         END DO
C 
      END DO
C 
CD    CALL IMPTET ('NUMERO D''ENREGISTREMENT DES MATRICES  '//IDPROG,
CD                  TNUENR, 1, NBDEV)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Mise a zero des efforts correspondant aux symetries
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT le nombre de fonctions du temps
C     ES...... EFFORT valeur des efforts en entree ranges
C     ES              (nddl, nfonct, nbmat)
C 
C     Et on recupere :
C 
C     ES...... EFFORT valeur des efforts en sortie
C 
      SUBROUTINE MZSYM (NFONCT, EFFORT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION EFFORT(NDDL*NBMAT*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   MAT, I, DEBUT, NFONCT
      INTEGER   ADEPPA, NONMAT
C 
CD    LOGICAL LTRACP, LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MZSYM')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL INFOEN('MAT-NULLES', ADEPPA , NONMAT )
C 
CD    IF ( LTRACP(1) ) THEN
C 
CD      DO I = 1 , NBMAT
C 
CD        CALL OMPTDP(' EFFORT EN ENTREE ',EFFORT(DEBUT),NDDL,NFONCT)
C 
CD      END DO
C 
CD      DEBUT = DEBUT+NDDL*NFONCT
C 
CD    END IF
C 
      DO MAT = 1 , NONMAT
C 
        DEBUT =  (M( ADEPPA +MAT-1) -1)*NDDL
C 
        DO I =  1 , NDDL
C 
          EFFORT( DEBUT + I ) = 0.D0
C 
        END DO
C 
      END DO
C 
CD    IF ( LTRACP(1) ) THEN
C 
CD      CALL IMPTDP(' EFFORT EN SORTIE ',EFFORT(1),NDDL,NBMAT)
C 
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END

