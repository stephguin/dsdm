C     On envoie comme arguments :
C 
C     E ....... NP le numero de comportement
C 
C     Et on recupere :
C 
C     S ....... NC le nombre d'interfaces concernees NC
C     S ....... C  leurs numeros dans le tableau C

      SUBROUTINE NUINT(NP, NC, C)
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
      INTEGER NP , NC ,C(NBINT),I, K ,  M1, R1
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NUINT')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBM('TYP-COUCHE',M1)
      R1=M1-1+NBCOU
C 
      K=0
C 
      DO I=1,NBINT
        IF(M(R1+I).EQ.NP)THEN
          K=K+1
          C(K)=I
        ENDIF
      ENDDO
C 
      NC=K
C 
CD    CALL IMPEN('POUR LE COMPORTEMENT DE L''INTERFACE NP',NP)
CD    CALL IMPEN('NOMBRE D''INTERFACES CONCERNEES',NC)
CD    CALL IMPTEN('NUMEROS DES INTERFACES CONCERNEES',C(1),1,NC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine remplit le tableau de replacage pour les interfaces :
C 
C                - le  tableau RANCAL-INT qui donne la numerotation
C                  de calcul U,V,W rangee dans un ordre tel que les
C                  numeros reels correspondants dans les matrices
C                  stockees profil soient en ordre croissant.
C                  I > J ==== > ddlU (TAB(I)) > ddlU( TAB(J))
C                - le  tableau RANINV-INT qui a partir de la
C                  numerotation U croissant donne LA NUMEROTATION
C                  calcul:
C                  U(1)  , U(2)  , U(4)  , U(3)  puis:
C                  U'(1) , U'(2) , U'(4) , U'(3) ou
C                  1 , 2 , 3 , 4 correspond a la numerotation
C                  de l'interface
C                  4-------3
C                  1-------2
C  
      SUBROUTINE REPLAI
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER APM1, APM2, APM1P, APM2P
      INTEGER APM3, APM3P
      INTEGER AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='REPLAI')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     creation d'un tableau provisoire pour ranger les numeros
C     des noeuds dans l'ordre correspondant au rangement
C     croissant des numeros reels
C           -1ere adresse:APM1
C -----------------------------------------------------------------------
C 
      CALL GESTEN('RANCAL-INT',8,APM1)
      CALL GESTEN('RANINV-INT',8,APM2)
      CALL GESTEN('RANTAB-INT',8,APM3)
C 
C -----------------------------------------------------------------------
C     POUR LE CAS SYMPAR
C -----------------------------------------------------------------------
C 
      CALL GESTEN('RANCAL-INP',8,APM1P)
      CALL GESTEN('RANINV-INP',8,APM2P)
      CALL GESTEN('RANTAB-INP',8,APM3P)
C 
C -----------------------------------------------------------------------
C 
C     ce tableau depend du type de numerotation choisie
C 
      IF(NUM.EQ.1)THEN
C  
C                           3,7 -------- 4,8
C                                                ordre des ddl
C     rangement calcul      1,5 -------- 2,6
C 
C                           3,7 -------- 4,8
C                                                toujours les meme numero
C     rangement tableau     1,5 -------- 2,6
C 
C                           5,6 -------- 7,8
C     rangement croissant                         ordre des ddl
C                           1,2 -------- 3,4
C 
C 
C     POUR RANCAL-INT
C     ddl(rancal ) => ddl range croissant a partir de ddl range calcul
C     ou tableau croissant
C 
      M(APM1)    = 1
      M(APM1+1)  = 5
      M(APM1+2)  = 2
      M(APM1+3)  = 6
      M(APM1+4)  = 3
      M(APM1+5)  = 7
      M(APM1+6)  = 4
      M(APM1+7)  = 8
C 
C     POUR RANINV-INT
C     ddl(raninv ) => ddl range calcul a partir de ddl croissant
C 
      M(APM2)    = 1
      M(APM2+1)  = 3
      M(APM2+2)  = 5
      M(APM2+3)  = 7
      M(APM2+4)  = 2
      M(APM2+5)  = 4
      M(APM2+6)  = 6
      M(APM2+7)  = 8
C 
C     POUR RANTAB-INT
C     ddl(raninv ) => ddl range tableau a partir de ddl croissant
C 
      M(APM3)    = 1
      M(APM3+1)  = 3
      M(APM3+2)  = 5
      M(APM3+3)  = 7
      M(APM3+4)  = 2
      M(APM3+5)  = 4
      M(APM3+6)  = 6
      M(APM3+7)  = 8
C 
C                           5,7 -------- 6,8
C                                                ordre des ddl
C     rangement calcul      X,X -------- X,X
C 
C                           3,7 -------- 4,8
C                                                toujours les meme numero
C     rangement tableau     1,5 -------- 2,6
C 
C                           5,6 -------- 7,8
C     rangement croissant                         ordre des ddl
C                           X,X -------- X,X
C 
C     ddl(rancal ) => ddl range croissant a partir de ddl range calcul
C     ou tablaeau croissant
C 
      M(APM1P)    = 1
      M(APM1P+1)  = 5
      M(APM1P+2)  = 2
      M(APM1P+3)  = 6
      M(APM1P+4)  = 3
      M(APM1P+5)  = 7
      M(APM1P+6)  = 4
      M(APM1P+7)  = 8
C 
C   ddl(raninv ) => ddl range calcul a partir de ddl croissant
C 
      M(APM2P)    = 1
      M(APM2P+1)  = 3
      M(APM2P+2)  = 5
      M(APM2P+3)  = 7
      M(APM2P+4)  = 2
      M(APM2P+5)  = 4
      M(APM2P+6)  = 6
      M(APM2P+7)  = 8

C     POUR RANTAB-INP
C     ddl(raninv ) => ddl range tableau a partir de ddl croissant
      M(APM3P)    = 1
      M(APM3P+1)  = 3
      M(APM3P+2)  = 5
      M(APM3P+3)  = 7
      M(APM3P+4)  = 2
      M(APM3P+5)  = 4
      M(APM3P+6)  = 6
      M(APM3P+7)  = 8
C 
      ENDIF
C 
      IF(NUM.EQ.2)THEN
C 
C                           3,7 -------- 4,8
C 
C     rangement calcul      1,5 -------- 2,6
C 
C                           3,7 -------- 4,8
C                                                 toujours les meme numero
C     rangement tableau     1,5 -------- 2,6
C 
C                           3,4 -------- 7,8
C     rangement croissant
C                           1,2 -------- 5,6
C 
C     POUR RANCAL-INT
C     ddl(rancal ) => ddl range croissant a partir de ddl range calcul
C 
      M(APM1)    = 1
      M(APM1+1)  = 5
      M(APM1+2)  = 3
      M(APM1+3)  = 7
C 
C     modif        M(APM1+3)  = 8
C 
      M(APM1+4)  = 2
      M(APM1+5)  = 6
      M(APM1+6)  = 4
      M(APM1+7)  = 8
C 
C     modif        M(APM1+7)  = 7
C     POUR RANINV-INT
C     ddl(raninv ) => ddl range calcul a partir de ddl croissant
C 
      M(APM2)    = 1
      M(APM2+1)  = 5
      M(APM2+2)  = 3
      M(APM2+3)  = 7
      M(APM2+4)  = 2
      M(APM2+5)  = 6
C  
C     modif de aout 92
C  
C     M(APM2+6)  = 8
C     M(APM2+7)  = 4
      M(APM2+6)  = 4
      M(APM2+7)  = 8
C  
C     POUR RANTAB-INT
C     ddl(raninv ) => ddl range tableau a partir de ddl croissant
C 
      M(APM3)    = 1
      M(APM3+1)  = 5
      M(APM3+2)  = 3
      M(APM3+3)  = 7
      M(APM3+4)  = 2
      M(APM3+5)  = 6
      M(APM3+6)  = 4
      M(APM3+7)  = 8
C 
C     modif        M(APM3+6)  = 8
C     modif        M(APM3+7)  = 4
C 
C                           5,7 -------- 6,8
C 
C     rangement calcul      x,x -------- x,x
C 
C 
C                           3,7 -------- 4,8
C 
C     rangement tableau     1,5 -------- 2,6
C 
C                           5,6      -------- 7,8
C     rangement croissant
C                         -100,-100 -------- -100,-100
C 
C 
C     POUR RANCAL-INP
C     ddl(rancal ) => ddl range croissant a partir de ddl calcul
C     ou tableau croissant
C 
      M(APM1P)    = 1
      M(APM1P+1)  = 5
      M(APM1P+2)  = 3
      M(APM1P+3)  = 7
      M(APM1P+4)  = 3
      M(APM1P+5)  = 7
      M(APM1P+6)  = 4
      M(APM1P+7)  = 8
C 
C     POUR RANINV-INP
C     ddl(raninv ) => ddl range calcul a partir de ddl croissant
C 
      M(APM2P)    = 1
      M(APM2P+1)  = 5
      M(APM2P+2)  = 3
      M(APM2P+3)  = 7
      M(APM2P+4)  = 5
      M(APM2P+5)  = 7
      M(APM2P+6)  = 6
      M(APM2P+7)  = 8
C 
C     POUR RANTAB-INP
C     ddl(raninv ) => ddl range tableau a partir de ddl croissant
C 
      M(APM3P)    = 1
      M(APM3P+1)  = 5
      M(APM3P+2)  = 3
      M(APM3P+3)  = 7
      M(APM3P+4)  = 3
      M(APM3P+5)  = 7
      M(APM3P+6)  = 4
      M(APM3P+7)  = 8
C 
      ENDIF
C 
CD    CALL IMPTET('TABLEAU RANCAL-INT',M(APM1),1,8)
CD    CALL IMPTET('TABLEAU RANINV-INT',M(APM2),1,8)
CD    CALL IMPTET('TABLEAU RANTAB-INT',M(APM3),1,8)
CD    CALL IMPTET('TABLEAU RANCAL-INP',M(APM1P),1,8)
CD    CALL IMPTET('TABLEAU RANINV-INP',M(APM2P),1,8)
CD    CALL IMPTET('TABLEAU RANTAB-INP',M(APM3P),1,8)
C 
C     remise des adresses des tableaux partiels a "zero"
C 
      CALL SOPOUB(AM2LC,ADM2LC)
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... A   demi-longueur de l'element
C     E ...... NB  dimension de la matrice a multiplier
C     E ...... LIN ligne d'entree (necdelamirement 12 OU 8?)
C 
C     Et on recupere :
C 
C     S ...... LINA(1..6)=LIN(1..6) ET LINA(7..12) = A*LIN( 7..12)
C 
      SUBROUTINE LINLOC (A, NB, LIN, LINA)
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
      INTEGER I, NB, DEMI
C 
      DOUBLE PRECISION  A, LIN(NB), LINA(NB)
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='LINLOC')
C 
CD    CALL WLKBCD(IDPROG)
C 
CD    IF(MOD(NB,2).NE.0)THEN
CD       CALL ERREUD(0,
CD      'La dimension de la ligne dans calk0 est impaire ')
CD    ENDIF
C 
      DEMI=NB/2
      DO I=1,DEMI
        LINA(I)=LIN(I)
      ENDDO
      DO I=DEMI+1,NB
        LINA(I)=A*LIN(I)
      ENDDO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0I
C 
C     VERSION DU     14 /  3  / 87
C 
C     On envoie comme argument :
C 
C     E ........ A    demi-largeur de l'element
C     E ........ RC   rayon de l'element
C     E ........ ISO  comportement isotrope transverse de l'interface
C     E ........ TI   adresse de depart des tableaux integres pour
C                     l'interface
C 
C     Et on recupere :
C 
C     S ........ TAB1, TAB3 le tableau utile a l'assemblage
C 
      SUBROUTINE CALKI (A, RC, ISO, TI, TAB1, TAB3)
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
      INTEGER  TI,I,ZI,K,L,AD0,AP0,ADR0,APR0,ADR1,NUMERO
      INTEGER     AM2LC      , ADM2LC
C 
      DOUBLE PRECISION  ISO(3),TAB1(64),TAB3(64),A,RC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALKI ')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C    Creation d'un tableau provisoire de premiere adresse ad0 pour
C    ranger les adresses de depart des tableaux integres qui nous interessent
C -----------------------------------------------------------------------
CD    CALL IMPDN('VALEUR DE A DANS CALKI',A)
CD    CALL IMPDN('VALEUR DE K1ISO DANS CALKI',ISO(1))
CD    CALL IMPDN('VALEUR DE K3ISO DANS CALKI',ISO(2))
CD    CALL IMPDN('VALEUR DE RC  DANS CALKI',RC)
C 
      CALL GSPOUE(2, AD0)
      AP0=AD0-1
      K=AD0
      M(K)=TI
      K=K+1
      M(K)=TI+64
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR0 pour
C     ranger les constantes issues du comportement qui nous interessent
C -----------------------------------------------------------------------
      DO NUMERO=1,2
        CALL GSPOUD(2*65,ADR0)
        ADR1=ADR0+2
        APR0=ADR0-1
        L=ADR0
        DM(L)=ISO(NUMERO)*A*RC
C 
CD      CALL IMPDN('VALEUR DE K1ISO*A*RC',DM(L))
C 
        DM(L+1)=ISO(NUMERO)*A*A
C 
CD      CALL IMPDN('VALEUR DE K1IS0*A*A',DM(L+1))
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR1 pour
C     ranger les ... issues du comportement  qui nous interessent
C -----------------------------------------------------------------------
        DO I=1,2
C 
C     TCD CALL MUMARE(DM(APR0+I),64,DM(M(AP0+I)),DM(ADR1+64*(I-1)))
C 
          CALL HOMAT(DM(APR0+I),DM(M(AP0+I)),DM(ADR1+64*(I-1)),64,1)
        ENDDO
        ZI=ADR1+64
C 
C     TCD CALL ADDMAD(64,DM(ZI-64),DM(ZI),DM(ZI))
C 
        CALL ADD(DM(ZI-64),DM(ZI),DM(ZI),64,1)
C 
C     RANGEMENT DANS LES TABLEAUX DE SORTIE DE TAB1 et TAB3
C 
        IF(NUMERO.EQ.1)THEN
          CALL MATLOC(8,A,DM(ADR1+64),TAB1)
        ELSE
          CALL MATLOC(8,A,DM(ADR1+64),TAB3)
        ENDIF
      ENDDO
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD      CALL IMPTDN('VALEUR DE TAB1'//IDPROG,TAB1,8,8)
CD      CALL IMPTDN('VALEUR DE TAB3'//IDPROG,TAB3,8,8)
C 
C     REMISE DES TABLEAUX PARTIELS A ZERO
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
C     CECI DOIT SE TROUVER DANS CALK0I
C 
C     VERSION DU     14 /  2  / 87
C 
C     On envoie comme arguments ;
C 
C     E ....... NUINT  le numero d'interface
C     E ....... NUCOL  le numero de colonne
C     E ....... TLOCN1 l'adresse de depart du tableau TLOCN1
C 
C     Et on recupere :
C 
C     S ....... NDDLU,NDDLV,NDDLW
C               tableau des numeros reels des ddl correspondant a
C               U ,V, W, stockes croissants
C 
      SUBROUTINE NDDLIU(NUINT, NUCOL, TLOCN1, NDDLU, NDDLV, NDDLW)
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
      INTEGER NDDLU(8),NDDLV(8),NDDLW(8),I,NUINT,NUCOL,NI,NELAV,TLOCN1
      INTEGER INDICI
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NDDLIU')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
       NELAV=TLOCN1+6*NBCOL*NBCOU-1
     & +4*(NBCOL*(NUINT-1)+NUCOL-1)
C 
CD    CALL IMPEN('POUR L''INTERFACE NUMERO ',NUINT)
CD    CALL IMPEN('ET LA COLONNE NUMERO ',NUCOL)
CD    CALL IMPEN('LE NOMBRE D''ELEMENTS AVANT EST',NELAV)
C 
      DO I=1,4
C 
CD      CALL IMPEP('POUR I ',I)
C 
C     NI EST LE PREMIER NUMERO DE DDL DE L'INTERFAC NINT 
C     POUR LE NOEUD NUMERO I
C 
        NI        = 6*(M(NELAV+I)-1)+1
        INDICI    = 2*(I-1)+1
C 
CD      CALL IMPEP('PREMIER NUMERO DE DDL DE L''INTERFACE ',NI)
C 
        NDDLU(INDICI)           = NI
        NDDLU(INDICI+1)         = NI+1
        NDDLV(INDICI)           = NI+2
        NDDLV(INDICI+1)         = NI+3
        NDDLW(INDICI)           = NI+4
        NDDLW(INDICI+1)         = NI+5
C 
      ENDDO
C 
CD    CALL IMPEN('POUR NUCOU='//IDPROG,NUINT)
CD    CALL IMPEN('POUR NUCOL='//IDPROG,NUCOL)
CD    CALL IMPTEN(
CD    'TABLEAU DES NUMEROS DE DDL STOCKES CROISSANTS U',NDDLU,1,8)
CD    CALL IMPTEN(
CD    'TABLEAU DES NUMEROS DE DDL STOCKES CROISSANTS V',NDDLV,1,8)
CD    CALL IMPTEN(
CD    'TABLEAU DES NUMEROS DE DDL STOCKES CROISSANTS W',NDDLW,1,8)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0I
C 
C     VERSION DU     14 /  3  / 87
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine assemble pour un element donne
C     les termes provenant de CALKI pour toutes les matrices KOn
C 
C     On envoie comme arguments :
C 
C     E ...... N        valeur du rang de Fourier
C     E ...... TK0      tableau des adresses de depart des matrices KON
C     E ...... UCAL     tableau des adresses elementaires UCAL
C     E ...... TAB1     tableau du terme independant de n provenant de
C                       CALKI1 pour l'element considere (le meme pour U et V)
C     E ...... TAB3     tableau du terme independant de n provenant de
C                       CALKI3 pour l'element considere
C     E ...... NDDLU    tableau des numeros reels des ddl du deplacement
C                       U de l'element d'interface considere
C     E ...... NDDLV    tableau des numeros reels des ddl du deplacement
C                       V de l'element d'interface considere
C     E ...... NDDLW    tableau des numeros reels des ddl du deplacement
C                       U de l'element d'interface considere
C     E ...... P1       1ere adresse du tableau profil

      SUBROUTINE ASSBI (N, TK0N, UCAL, TAB1, TAB3, NDDLU,
     &                  NDDLV, NDDLW, P1)
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
      INTEGER NDDLV(8), NDDLW(8), NVI, NWI, NVJ, NWJ, ADV1, ADW1
      INTEGER TK0N, PROFU(8), PROFV(8), PROFW(8)
      DOUBLE PRECISION  TAB1(8,8), TAB3(8,8), TCAL1,TCAL3
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ASSBI')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPTDP('VALEUR DE TAB1 DANS '//IDPROG,TAB1(1,1),8,8)
CD    CALL IMPTDP('VALEUR DE TAB3 DANS '//IDPROG,TAB3(1,1),8,8)
CD    CALL IMPTEN('RANG U CROISSANT '//IDPROG,NDDLU(1),1,8)
CD    CALL IMPTEN('RANG V CROISSANT '//IDPROG,NDDLV(1),1,8)
CD    CALL IMPTEN('RANG W CROISSANT '//IDPROG,NDDLW(1),1,8)
CD    CALL IMPTEN('CORRESPONDANCE CALCUL'//IDPROG,UCAL(1),1,8)
CD    CALL IMPEN('VALEUR DE NTSFG DANS ASSBI',NTDSFG)
C 
      DO I=1,8
         PROFU(I)    =M(P1+NDDLU(I))
         PROFV(I)    =M(P1+NDDLV(I))
         PROFW(I)    =M(P1+NDDLW(I))
      ENDDO
C 
CD    CALL IMPTEN('PROF DE U',PROFU(1),1,8)
CD    CALL IMPTEN('PROF DE V',PROFV(1),1,8)
CD    CALL IMPTEN('PROF DE W',PROFW(1),1,8)
C 
      DO I=1,8
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
C     WRITE (BARATI(1:6),'(I6)')N
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
C     E ...... NUINT  le numero d'interface
C     E ...... NUCOL  le numero de colonne
C     E ...... TLOCN1 l'adresse de depart du tableau TLOCN1
C 
C     Et on recupere :
C 
C     S ...... NDDLU,NDDLV,NDDLW tableau des numeros reels des ddl 
C                                correspondants a U, V, W,stockes tableau
C 
      SUBROUTINE DDLICA (NUINT, NUCOL, TLOCN1, NDDLUC, NDDLVC, NDDLWC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NDDLUC(8), NDDLVC(8), NDDLWC(8), NUINT, NUCOL, TLOCN1
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER NDDLU(8), NDDLV(8), NDDLW(8), I
      INTEGER PLACE, RANINV
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DDLICA')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPEN('POUR numint='//IDPROG,NUINT)
C 
CD    CALL IMPEN('POUR NUCOL='//IDPROG,NUCOL)
C 
      CALL ADTBM('RANINV-INT',RANINV)
      IF(SYMPAR .AND.NUINT.EQ.1)CALL ADTBM('RANINV-INP',RANINV)
C 
      RANINV = RANINV -1
C 
      CALL NDDLIU( NUINT,NUCOL,TLOCN1, NDDLU,NDDLV,NDDLW)

      DO I = 1 , 8
C 
        PLACE     = M(RANINV+I)
        NDDLUC(I) = NDDLU( PLACE )
        NDDLVC(I) = NDDLV( PLACE )
        NDDLWC(I) = NDDLW( PLACE )
C 
      END DO
C 
CD    CALL IMPTEN(
CD    'tableau des no de ddl de stockes calcul U',NDDLUC,1,8)
CD    CALL IMPTEN(
CD    'tableau des no de ddl de stockes calcul V',NDDLVC,1,8)
CD    CALL IMPTEN(
CD    'tableau des no de ddl de stockes calcul W',NDDLWC,1,8)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule la trace au point de Gauss des NFONCT
C     fonctions de l'espace WK par les NFONCT fonction VL.
C     Le resultat de cette routine est le tableau TLKJI a 4 indices
C     tel que, en rangement informatique :
C 
C     I = 1 , NFONCT
C       J = 1 , I
C         K = 1 , NK
C           L = 1 , NL
C                                        _
C     Pour les interfaces
C 
C             TLKJI ( l , k , j , i ) =  Wk Kij Vl
C 
C 
C     La matrice KT etant symetrique, on ne fait varier J que de
C     1 a I, ce qui correpond au rangement de KIJPG.
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions
C     E ...... KIJPG    le tableau des matrices au point de
C                       Gauss range (9,[ i , j])
C     E ...... WK       le 1er tableau de sauts range (NSAU,NK)
C     E ...... VL       le 2eme tableau de sauts range (NSAU,NL)
C 
C     Et on recupere :
C 
C     ES...... TLKJIE  le tableau des traces au point de gauss
C                      en entree range comme indique plus haut
C     ES...... TLKJIS  le tableau des traces au point de gauss
C                      en sortie range comme indique plus haut
C 
C     REM : inutile de mettre TLKJI a 0. On peut faire TLKJIE + TLKJIS.

         SUBROUTINE SIBFPG (NFONCT, MULT, KIJPG, NK, WK, NL, VL,
     &                      TLKJIE, TLKJIS)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NFONCT  , NK , NL
C 
      DOUBLE PRECISION  KIJPG(9*NFONCT*(NFONCT+1)/2)
      DOUBLE PRECISION  WK (NSAU*NK)   , MULT
      DOUBLE PRECISION  VL (NSAU*NL)
      DOUBLE PRECISION  TLKJIE( NL , NK , NFONCT , NFONCT )
      DOUBLE PRECISION  TLKJIS( NL , NK , NFONCT , NFONCT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION TRACE
C 
      INTEGER     I , J , K , L , DEBKIJ , DEPSK  , DEPSL
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SIBFPG')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     DEBUT DE LA MATRICE KIJ DANS KIJPG
C 
      DEBKIJ = 1
C 
C     DEBSK DEBUT DU SAUT WK
C     DEBSL DEBUT DU SAUT VL
C 
      DO  I = 1, NFONCT
C 
        DO   J = 1, I
C 
          DEPSK = 1
C 
          DO   K = 1, NK
C 
            DEPSL = 1
C 
            DO   L = 1, NL
C 
              CALL PROI12 (KIJPG(DEBKIJ), VL(DEPSL), WK(DEPSK),
     &                     TRACE)
C 
              TLKJIS (L, K, J, I) = TLKJIE (L, K, J, I)
     &                                      + TRACE*MULT
C 
              DEPSL = DEPSL+NSAU
C 
            END DO
C 
            DEPSK = DEPSK+NSAU
C 
          END DO
C 
          DEBKIJ = DEBKIJ+9
C 
        END DO
C 
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
C 
CD      CALL OMPTDN( ' TLKJIS SOUS LA FORME (NL*NK,N2) ',
CD                     TLKJIS(1,1,1,1) , NL*NK ,NFONCT*NFONCT )
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
C     Cette routine calcule TRACE = TR( SAUT2 K(SAUT1) )
C 
C     On envoie comme arguments :
C 
C     E ...... K      la matrice stockee (9)
C     E ...... SAUT1  la 1ere deformation (3)
C     E ...... SAUT2  la 2eme deformation (3)
C 
C     Et on recupere :

C     S ...... TRACE

      SUBROUTINE PROI12 (K, SAUT1, SAUT2, TRACE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  TRACE , K(9) , SAUT1(NSAU) , SAUT2(NSAU)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION SIGMA( 3 )
C 
      INTEGER    I
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='PROI12')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL IULORT (K(1), SAUT1(1), SIGMA(1))
C 
C     Calcul de la trace
C 
      TRACE = 0.D0
C 
      DO I = 1 , 3
C 
        TRACE = TRACE + SAUT2(I)*SIGMA(I)
C 
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
C 
CD     CALL IMPDN(' VALEUR DE LA TRACE ', TRACE )
C 
CD    END IF
C 
C     SEQUENCE D'IMPRESSION EN TRACE PROFONDE
C 
CD    IF (LTRACP(1))THEN
C 
CD      CALL OMPTDP(' K EN ENTREE     ', K(1)      ,9 , 1)
CD      CALL OMPTDP(' SAUT1 EN ENTREE ', SAUT1(1)  ,3 , 1)
CD      CALL OMPTDP(' K ( SAUT1 )     ', SIGMA(1)  ,3 , 1)
CD      CALL OMPTDP(' SAUT2 EN ENTREE ', SAUT2(1)  ,3 , 1)
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
C     Cette routine calcule la trace au point de Gauss des NFONCT
C     fonctions de l'espace WK par les NFONCT fonctions VL.
C     Le resultat de cette routine est le tableau TLKJI a 4 indices
C     tel que, en rangement informatique :
C 
C     I = 1 , NFONCT
C       J = 1 , I
C         K = 1 , NK
C           L = 1 , K
C 
C     Pour les interfaces
C 
C             TLKJI ( l , k , j , i ) =  Wk Kij Vl
C 
C 
C     La matrice KT etant symetrique, on ne fait varier J que de
C     1 a I, ce qui correpond au rangement de KIJPG.
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions
C     E ...... MULT     le multiplicateur d'integration
C     E ...... KIJPG    le tableau des matrices au point de
C                       Gauss range (9, [i,j])
C     E ...... NK       le nombre de vecteurs WK
C     E ...... WK       le 1er tableau de sauts range (NSAU,NK)
C 
C            ATTENTION !! WK = VL NK = NL
C 
C     E ...... NL       le nombre de vecteurs VL
C     E ...... VL       le 2eme tableau de sauts range (NSAU,NK)
C 
C     Et on recupere :
C 
C     ES ..... TLKJIE    le tableau des traces au point de Gauss
C                        en entree ranges comme indique plus haut
C     ES ..... TLKJIS    le tableau des traces au point de Gauss
C                        en sortie ranges comme indique plus haut
C 
C     REM : inutile de mettre TLKJI a 0.

      SUBROUTINE NIBFPG (NFONCT, MULT, KIJPG, NK, WK, NL, VL,
     &                   TLKJIE, TLKJIS)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NFONCT , NK , NL
C 
      DOUBLE PRECISION  KIJPG(9*NFONCT*(NFONCT+1)/2)
      DOUBLE PRECISION  WK (NSAU*NK ) , MULT
      DOUBLE PRECISION  VL (NSAU*NL )
      DOUBLE PRECISION  TLKJIE( NL , NK , NFONCT , NFONCT )
      DOUBLE PRECISION  TLKJIS( NL , NK , NFONCT , NFONCT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION TRACE
C 
      INTEGER     I , J , K , L , DEBKIJ , DEPSK  , DEPSL
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NIBFPG')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF ( NL . NE . NK  ) THEN
C 
        CALL ERREUD( 0, ' NL diff  de NK mauvaise utilisation de '
     &                  //IDPROG)
C 
      END IF
C 
C     Debut de la matrice KIJ dans KIJPG
C 
      DEBKIJ = 1
C 
C     DEBSK debut du  saut WK
C     DEBSL debut du  saut VL
C 
      DO  I = 1, NFONCT
C 
        DO   J = 1 , I
C 
          DEPSK = 1
C 
          DO   K = 1 , NK
C 
            DEPSL = 1
C 
            DO   L = 1 , K
C 
              CALL PROI12 (KIJPG(DEBKIJ), VL(DEPSL), WK(DEPSK),
     &                     TRACE)
C 
              TLKJIS( L , K , J , I ) = TLKJIE( L , K , J , I )
     &                                  + MULT*TRACE
C 
              DEPSL = DEPSL+NSAU
C 
            END DO
C 
            DEPSK = DEPSK+NSAU
C 
          END DO
C 
          DEBKIJ = DEBKIJ+9
C 
        END DO
C 
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
C 
CD      CALL OMPTDN( ' TLKJIS SOUS LA FORME (NK*NL,N2) ',
CD                     TLKJIS(1,1,1,1) , NK*NL ,NFONCT*NFONCT )
C 
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END

