C     TAB S7 2; () 1:72
C ----------------------------------------------------------------------
C 
C      PREVISION DU DELAMINAGE D'UNE PLAQUE TROUEE
C 
C      L.M.T.-------------------CACHAN----------------------------------
C 
C      VERSION DU     17 /  10  / 86
C 
C      MODIFIEE LES :
C 
C      ALLIX OLIVIER
C 
C ----------------------------------------------------------------------
C      QUE FAIT CETTE SUBROUTINE :
C 
C         Routine calculant aux points de gauss les termes
C         independants de la geometrie de la matrice E.F. des
C         deformations. Chacune des fonctions de base est calculee
C         au points de gauss necdelamires puis rangees dans un tableau
C         de nom TAB-GAUSS.
C 
C         Les tableaux du type NN sont ensuite calcules pour tous
C         les points de gauss et ranges dans le tableau MAT-CONST.
C         Le meme travail est effectue pour les termes d'interface.
C         On calcule d'abord les termes N(x,1),N(x,-1) et on les
C         range dans le tableau GAUSSINTER. Les produits interessants
C         sont ranges dans le tableau MAT-INTERF
C 
C         ATTENTION MODIF :
C 
C         Les fonctions sont rangees dans l'ordre de calcul,
C         c'est a dire :
C 
C         U=N1*u1+N3*u2+N5*u3+N7*u4+N9*u5+N11*u6
C          +N2*u'1+N4*u'2+N6*u'3+N8*u'4+N10*u'5+N12*u'6
C 
C ----------------------------------------------------------------------
C 
      SUBROUTINE MATCON
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
      INTEGER P, Q, I, J, ADLC1, ADLC2
      INTEGER VAL, VALI, H, K, DEBMAT
      INTEGER LONG, INDIC, INDIC2, ADLCI1, ADLCI2, LONGI
C 
      DOUBLE PRECISION N1, N2, N3, N4, N5, N6
      DOUBLE PRECISION N7, N8, N9, N10, N11, N12
      DOUBLE PRECISION NY1, NY2, NY3, NY4, NY5, NY6
      DOUBLE PRECISION NY7, NY8, NY9, NY10, NY11
      DOUBLE PRECISION NX1, NX2, NX3, NX4, NX5, NX6
      DOUBLE PRECISION NX7, NX8, NX9, NX10, NX11
      DOUBLE PRECISION NX12, NY12, H1, H2, H3, H4
      DOUBLE PRECISION X, Y
      CHARACTER*6 IDPROG
C 
C     POUR TEST
C 
      DOUBLE PRECISION WI, W
      INTEGER TEST1, TEST2
      PARAMETER (IDPROG='MATCON')
C 
CD    CALL WLKBCD(IDPROG)
C 
C ----------------------------------------------------------------------
C     POUR TEST
C 
      P=XINTEG
      Q=YINTEG
C 
C ----------------------------------------------------------------------
C     Creation d'un tableau  pour ranger les valeurs
C     des fonctions elementaires aux differents points de gauss
C     1ere adresse:ADLC1
C ----------------------------------------------------------------------
      CALL GESTDP('TAB-GAUSS ',36*P*Q,ADLC1)
      CALL MENADM(ADLC1,36*P*Q)
CD    CALL IMPEN('T1',T1)
      VAL=ADLC1-1
C 
      CALL DEBUDP(DEBMAT)
      ADLC2 = DEBMAT
C 
      INDIC=0
      INDIC2=0
      K=P*(P-1)/2
      H=Q*(Q-1)/2
C 
      DO I=1,P
        X=GAUSS(K+I)
        WI=POIDS(K+I)
C 
        DO J=1,Q
C 
        Y=GAUSS(H+J)
        W=WI*POIDS(H+J)
C 
        DM(VAL+1)=N1(X,Y)
C 
        DM(VAL+2)=N3(X,Y)
C 
        DM(VAL+3)=N5(X,Y)
C 
        DM(VAL+4)=N7(X,Y)
C 
        DM(VAL+5)=N9(X,Y)
C 
        DM(VAL+6)=N11(X,Y)
C 
        DM(VAL+7)=N2(X,Y)
C 
        DM(VAL+8)=N4(X,Y)
C 
        DM(VAL+9)=N6(X,Y)
C 
        DM(VAL+10)=N8(X,Y)
C 
        DM(VAL+11)=N10(X,Y)
C 
        DM(VAL+12)=N12(X,Y)
C 
        DM(VAL+13)=NX1(X,Y)
C 
        DM(VAL+14)=NX3(X,Y)
C 
        DM(VAL+15)=NX5(X,Y)
C 
        DM(VAL+16)=NX7(X,Y)
C 
        DM(VAL+17)=NX9(X,Y)
C 
        DM(VAL+18)=NX11(X,Y)
C 
        DM(VAL+19)=NX2(X,Y)
C 
        DM(VAL+20)=NX4(X,Y)
C 
        DM(VAL+21)=NX6(X,Y)
C 
        DM(VAL+22)=NX8(X,Y)
C 
        DM(VAL+23)=NX10(X,Y)
C 
        DM(VAL+24)=NX12(X,Y)
C 
        DM(VAL+25)=NY1(X,Y)
C 
        DM(VAL+26)=NY3(X,Y)
C 
        DM(VAL+27)=NY5(X,Y)
C 
        DM(VAL+28)=NY7(X,Y)
C 
        DM(VAL+29)=NY9(X,Y)
C 
        DM(VAL+30)=NY11(X,Y)
C 
        DM(VAL+31)=NY2(X,Y)
C 
        DM(VAL+32)=NY4(X,Y)
C 
        DM(VAL+33)=NY6(X,Y)
C 
        DM(VAL+34)=NY8(X,Y)
C 
        DM(VAL+35)=NY10(X,Y)
C 
        DM(VAL+36)=NY12(X,Y)
C 
CD      CALL IMPDN(' POUR XINTEG',X)
CD      CALL IMPDN('POUR YINTEG',Y)
CD      CALL IMPTDN('FONCTIONS N AUX POINTS DE GAUSS '
CD      ,DM(VAL+1),36,1)
C 
        INDIC=INDIC+1
C 
C     CALCUL DE NN
C 
CD      CALL IMPEN('AVANT APPEL A VECVEC',INDIC)
        CALL MATM02(DM(VAL+1),DM(VAL+1),DM(ADLC2),12,1,12)
C 
C     TCD CALL VVECVE(DM(VAL+1),DM(VAL+1),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDN('TABLEAUX CONSTANTS NN',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
C                _ 
C     CALCUL DE NNx
C 
        CALL MATM02(DM(VAL+1),DM(VAL+13),DM(ADLC2),12,1,12)
C 
C     TCD CALL VVECVE(DM(VAL+1),DM(VAL+13),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NNx',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
C               _
C     CALCUL DE NxN
C 
        CALL MATM02(DM(VAL+13),DM(VAL+1),DM(ADLC2),12,1,12)
C 
C     TCD CALL VVECVE(DM(VAL+13),DM(VAL+1),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NxN',DM(ADLC2),12,12)
        ADLC2=ADLC2+144
C               _
C     CALCUL DE NNy
C 
        CALL MATM02(DM(VAL+1),DM(VAL+25),DM(ADLC2),12,1,12)
C 	
C     TCD CALL VVECVE(DM(VAL+1),DM(VAL+25),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NNy',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
C               _
C     CALCUL DE NyN
C 
        CALL MATM02(DM(VAL+25),DM(VAL+1),DM(ADLC2),12,1,12)
C 
C     TCD CALL VVECVE(DM(VAL+25),DM(VAL+1),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NyN',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
C               _
C     CALCUL DE NxNx
C 
        CALL MATM02(DM(VAL+13),DM(VAL+13),DM(ADLC2),12,1,12)
C 
C     TCD CALL VVECVE(DM(VAL+13),DM(VAL+13),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NxNx',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
C               _
C     CALCUL DE NxNy
C 
        CALL MATM02(DM(VAL+13),DM(VAL+25),DM(ADLC2),12,1,12)
C 
C     TCD CALL VVECVE(DM(VAL+13),DM(VAL+25),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NxNy',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
C               _
C     CALCUL DE NyNx
C 
        CALL MATM02(DM(VAL+25),DM(VAL+13),DM(ADLC2),12,1,12)
C 
C     TCD CALL VVECVE(DM(VAL+25),DM(VAL+13),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NyNx',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
C               _
C     CALCUL DE NyNy
C 
        CALL MATM02(DM(VAL+25),DM(VAL+25),DM(ADLC2),12,1,12)
C 
C     TCD CALL VVECVE(DM(VAL+25),DM(VAL+25),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NyNy',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
C               _   _
C     CALCUL DE NNx+NxN
C 
        CALL ADD(DM(ADLC2-8*144),DM(ADLC2-7*144),DM(ADLC2),12,12)
C 
CD      CALL IMPDP(' POUR XINTEG',X)
CD      CALL IMPDP('POUR YINTEG',Y)
CD      CALL IMPTDP('TABLEAUX CONSTANTS NNx+NxN',DM(ADLC2),12,12)
C 
        ADLC2=ADLC2+144
        VAL=VAL+36
        ENDDO
      ENDDO
C 
      LONG=(ADLC2-DEBMAT)
CD    CALL IMPEN('LONGUEUR DE MAT-CONST',LONG)
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau de nom MAT-CONST, de 1ere adresse debmat
C     pour ranger les valeurs des tableaux des fonctions elementaires 
C     aux differents points de gauss
C -----------------------------------------------------------------------
      CALL GESTDP('MAT-CONST ',LONG,TEST1)
      CALL TESTEN (TEST1, DEBMAT,IDPROG)
C 
C -----------------------------------------------------------------------
C 
C     SEQUENCE DE CALCUL POUR LES INTERFACES
C -----------------------------------------------------------------------
      IF(NBINT.GT.0)THEN
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau des fonctions elementaires aux differents
C     points de gauss pour ranger les valeurs pour les interfaces, de
C     1ere adresse:ADLCI1
C -----------------------------------------------------------------------
        CALL GESTDP('GAUSSINTER',8*P,ADLCI1)
        CALL MENADM(ADLCI1,8*P)
        VALI=ADLCI1-1
        INDIC=0
C 
C     M CALL IMPEN('ADLCI1',ADLCI1)
C 
      CALL DEBUDP(DEBMAT)
      ADLCI2 = DEBMAT
      DO I=1,P
          X=GAUSS(K+I)
          DM(VALI+1)   = -1.D0*H1(X)
          DM(VALI+2)   = -1.D0*H3(X)
          DM(VALI+3)   = H1(X)
          DM(VALI+4)   = H3(X)
          DM(VALI+5)   = -1.D0*H2(X)
          DM(VALI+6)   = -1.D0*H4(X)
          DM(VALI+7)   = H2(X)
          DM(VALI+8)   = H4(X)
C 
CD        CALL IMPDN('POUR X',X)
CD        CALL IMPTDN(
CD        'FONCTIONS AUX POINTS DE GAUSS POUR LES INTERFACES'
CD        ,DM(VALI+1),8,1)
C 
           INDIC=INDIC+1
C             _
C     CALCUL DE H(x)*H(x)
C 
CD        CALL IMPDP(' POUR XINTEG',X)
          CALL MATM02(DM(VALI+1),DM(VALI+1),DM(ADLCI2),8,1,8)
C 
C     TCD CALL VECVEC(DM(VALI+1),DM(VALI+1),DM(ADLCI2),8,8)
C 
CD        CALL IMPTDP('TABLEAUX CONSTANTS H(x)*H(x)',DM(ADLCI2),8,8)
C 
          ADLCI2=ADLCI2+64
          VALI=VALI+8
C 
        ENDDO
C 
        LONGI=(ADLCI2-DEBMAT)
C 
CD      CALL IMPEN('LONGUEUR DE MAT-INTERF',LONGI)
C 
C -----------------------------------------------------------------------
C       Creation d'un tableau de 1ere adresse ADLCI2 pour ranger les
C       valeurs des tableaux des fonctions elementaires aux differents
C       points de gauss
C -----------------------------------------------------------------------
        CALL GESTDP('MAT-INTERF',LONGI,TEST2)
        CALL TESTEN (TEST2,DEBMAT,IDPROG)
      ENDIF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
