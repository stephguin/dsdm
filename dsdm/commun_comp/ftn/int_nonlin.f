C     Integration a saut imposee de la plasticite de l'interface
C 
C     On envoie comme arguments :
C 
C     E ...... ANGORT                      angle entre la base locale et la base d'orthotropie
C 
C     E ...... R0, AL, BE, A1, A2, K1, K2  les caracteristiques plastiques de l'interface
C     E ...... RPRE                        ancienne valeur du seuil
C     E ...... PPRE                        ancienne valeur de la plasticite cumulee
C     E....... SAUPOP                      saut plastique 1 et 2 precedent
C     E....... SAUEOP                      saut elastique precedent
C 
C     PROVIENT DES CHAMPS ADMISSIBLES
C 
C     E ...... SAPADM                      saut point admissible (nsig) X DT
C     E ...... SGPADM                      sgn point admissible (nsig) X DT
C     E ...... QADMPR                      quantites admissible de l'etape
C                                          locale precedente (neps, nsig) X DT
C     Et on recupere :
C 
C     POUR LA BASE D'ORTHOTROPIE
C 
C     S ...... RNOU                        nouvelle valeur du seuil
C     S ...... PNOU                        nouvelle valeur de la plasticite cumulee
C     S ...... SAUPOR                      saut plastique 1 et 2 actuel
C     S ...... SAUEOR                      saut elastique actuel
C 
C     POUR LA BASE LOCALE
C 
C     S ...... SAADMP                      saut admissible precedent
C     S ...... SAADM                       saut admissible actuel
C 
      SUBROUTINE INPLAI
C                 Pour la base d'orthotropie
     &           (ANGORT, R0, AL, BE, A1, A2, K1, K2,
     &            RPRE, PPRE, RAC, PAC, SAUPOP, SAUPOR, SAUEOP,
     &            SAUEOR,
C                 Pour la base locale
     &            QADMPR, SAPADM, SGPADM, SAADMP, SAADM)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'typcal.h'
C 
C     Caracteristiques de l'interface et du point
C 
      DOUBLE PRECISION ANGORT, R0, AL, BE, A1, A2, K1, K2
C 
C     Provient des champs admissibles de l'etape locale precedente
C 
      DOUBLE PRECISION QADMPR(6)
C 
C     PROVIENT DES CHAMPS ADMISSIBLES
C 
      DOUBLE PRECISION SAPADM(3), SAADMP(3)
C 
C     PROVIENT DE L'ITERATION PRECEDENTE
C 
      DOUBLE PRECISION  RPRE, PPRE, SAUPOP(2), SAUEOP(3)
C 
C     QUANTITES ACTUELLES
C 
      DOUBLE PRECISION SAADM(3), RAC, PAC, SAUEOR(3), SAUPOR(2)
      DOUBLE PRECISION SGPADM(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER          I, J
C 
      DOUBLE PRECISION CEFF(2), CEQ
      DOUBLE PRECISION TEST, DELP, SAUIMP(3)
C 
      CHARACTER*6 IDPROG
C 
C     Pour la verification
C 
CD    DOUBLE PRECISION CP1, CP2, SAUVE(2)
CD    LOGICAL LTRACP
C 
CD    PARAMETER (IDPROG='INPLAI')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     Calcul des accroissements des quantites admissibles
C     de 1 a 3  => sauts
C     de 4 a 6  => contraintes
C 
CD    CALL IMPTDT ('VERIF SAPADM ', SAPADM, 1, 6)
CD    CALL IMPTDT ('VERIF QADMPR ', QADMPR, 1, 6)
C 
      DO I = 1, 3
C 
C       Calcul du taux de saut admissible
C 
        QADMPR(I) = SAPADM(I) + QADMPR(I)
        SAPADM(I) = QADMPR(I)
        SAADM (I) = SAADMP(I) + SAPADM(I)
C 
C       Calcul du taux de contrainte admissible
C 
        J = I+3
C 
        QADMPR(J) = SGPADM(I) + QADMPR(J)
        SGPADM(I) = QADMPR(J)
C 
      END DO
C 
C     Calcul des sauts chapeau dans la base d'orthotropie
C     saut chapeau = saut admissible
C 
      CALL QIBORT (ANGORT, SAPADM, SAUEOR)
      CALL ADDMAD (3, SAUEOR, SAUEOP, SAUEOR)
C 
      CALL QIBORT (ANGORT, SAADM, SAUIMP)
C 
C     PREDICTION ELASTIQUE : CALCUL DE LA DEFORMATION ELASTIQUE EN
C     L'ABSENCE D'ACCROISSEMENT PLASTIQUE ET D'EVOLUTION DES MODULES
C 
C     Calcul des pseudo-contraintes effectives, en realite liees aux sauts elastiques :
C 
C     CEFF(I, 1) = S12EFF = G120*EPSILON12EL
C     CEFF(I, 2) = S22EFF = E20 *<E22EL>
C                                                                        +
C     ATTENTION E12 = R2 EPSILON 12
C 
C 
C     Calcul des pseudo-contraintes effectives, en realite liees aux sauts elastiques :
C 
      CEFF(1)  =  K1*SAUEOR(1)
      CEFF(2)  =  K2*SAUEOR(2)
C 
C     Pour verification de la qualite et de la validite de l'integration
C 
CD    CP1 = CEFF(1)
CD    CP2 = CEFF(2)
C 
      CEQ = DSQRT (A1*CEFF(1)*CEFF(1) + A2*CEFF(2)*CEFF(2))
C 
C     TEST POUR SAVOIR SI LE PREDICTEUR ELASTIQUE EST DANS LA SURFACE DE CHARGE
C 
      TEST = CEQ-RPRE
C 
      IF (.NOT. LPLAIN) THEN
        TEST = -1.D0
      END IF
C 
C     SI LE PREDICTEUR EST DANS LA SURFACE DE CHARGE ALORS LA SOLUTION
C     ELASTIQUE ENDOMMAGEABLE EST BONNE SINON :
C 
      IF ((TEST .GT. 0.D0)) THEN
C 
C       VALEUR DU SEUIL INITIAL
C 
        RAC = RPRE
C 
C       Resolution de (DABS (CEQ))+ R(p+dp)  =0
C       On calcule en realite le saut elastique
C 
        CALL VALDPI (R0, AL, BE, A1, A2, K1, K2,
     &               PPRE, CEFF(1), RAC, DELP)
C 
        PAC       = PPRE +DELP
C 
        SAUEOR(1)  = CEFF(1)/K1
        SAUEOR(2)  = CEFF(2)/K2
C 
        SAUPOR(1)  = SAUIMP(1)-SAUEOR(1)
        SAUPOR(2)  = SAUIMP(2)-SAUEOR(2)
C 
C        Pour verifier la qualite et de la validite de l'integration
C 
CD       SAUVE(1)  = -(CEFF(1)-CP1)/K1 + SAUPOP(1)
CD       SAUVE(2)  = -(CEFF(2)-CP2)/K2 + SAUPOP(2)
C 
CD       IF (DABS(SAUVE(1) - SAUPOR(1)) .GT. 1D-6) THEN
CD         CALL IMPDT ('SAUVE(1)'//IDPROG, SAUVE(1))
CD         CALL IMPDT ('EPSPOR(1)'//IDPROG, SAUPOR(1))
CD       END IF
C 
CD       IF (DABS(SAUVE(2) - SAUPOR(2)) .GT. 1D-6) THEN
CD         CALL IMPDT( 'EPSPVE(2)'//IDPROG, SAUVE(2))
CD         CALL IMPDT ('EPSPOR(2)'//IDPROG, SAUPOR(2))
CD       END IF
C 
      ELSE
C 
        SAUPOR(1)  = SAUPOP(1)
        SAUPOR(2)  = SAUPOP(2)
C 
        PAC        = PPRE
        RAC        = RPRE
C 
      END IF
C 
CD    IF (LTRACP(1)) THEN
CD      CALL OMPTDN ('QUANTITES ADMISSIBLES ACTUELLES '//IDPROG,
CD                    QADMPR, 6, 2)
CD      CALL OMPTDN ('SAUTS ADMISSIBLES ACTUELS  '//IDPROG,
CD                    QADMPR, 6, 2)
CD      CALL OMPTDN ('SAUTS ELASTIQUES (ORTHOTROPIE) 1 ET 2
CD                    ACTUELS ', SAUEOR, 2, 1)
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     MODIFIEE LE 25/07/96 PAR OA + DL
C     DL  Ajout du seuil YR le 3/09/96
C 
C     On envoie l'accroissement de temps deltat
C     On envoie les valeurs precedentes des endommagements
C 
C     Principe
C 
C     .               n
C     D = < Y' - d >       si d < 1  et  YD < YR
C                  +
C     Comme d1=d2=d, on ne calcule que d par la procedure 'inrprf'
C 
      SUBROUTINE VAENDI
C                       On envoie
     &                  (K, N, Y0, YC, YR, YD, MINT, DELTAT,
     &                  TESTD, Y, ENDOP,
C                       On recupere
     &                  ENDO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'typcal.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  K, N, Y0, YC, YR, YD, MINT, MULINT
      DOUBLE PRECISION  DELTAT ,  Y , ENDOP(3) ,ENDO(3)
      LOGICAL  TESTD
C 
      DOUBLE PRECISION  PAR
      INTEGER I
C 
C -----------------------------------------------------------------------
      IF (.NOT. LENINT) GOTO 5
C 
      IF (.NOT. TESTD)  GOTO 5
C 
C     Coefficient multiplicateur de (Y-Y0) (voir modele dans etaloc.f)
C 
      MULINT = (MINT + 1.D0) * (YC - Y0)
      MULINT = MINT / MULINT
C 
C     INITIALISATION DES VALEURS DE D ET DP
C 
      DO I = 1, 3
        ENDO(I) = ENDOP(I)
      END DO
C 
C     ON CALCULE L'EVOLUTION DE D
C     Calcul du Y'
C 
      PAR = MULINT * (Y - Y0)
      IF (PAR .LT. 0.D0) GOTO 5
      PAR = PAR ** MINT
C 
      CALL INRPRF (DELTAT, K, N, PAR, YD, YR, ENDOP(3), ENDO(3), 1.D0)
C 
      ENDO(3) = DMIN1 (ENDO(3), 1.D0)
      ENDO(2) = ENDO(3)
      ENDO(1) = ENDO(3)
C 
5     CONTINUE
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Integration de la plasticite pour l'interface
C 
C     On envoie comme arguments :
C 
C     E ...... R0, AL, BE  caracteristiques de plasticite  G  = G/1-D
C     E ...... PRE         valeur precedente de la deformation plastique equivalente
C     E ...... C           la valeur de la contrainte effective
C     E ...... RAC         la valeur du seuil de plasticite
C 
C     On recupere :
C 
C     S ...... DELP       l'accroissement de deformation plastique equivalente
C 
C 
      SUBROUTINE VALDPI (R0, AL, BE, A1, A2, K1, K2, PRE, C, RAC, DELP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'typcal.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  PRE, C(3), RAC
C 
      DOUBLE PRECISION  R0, AL, BE, A1, A2, K1, K2
C 
      DOUBLE PRECISION  CEQ, ACCP, CP1, CP2
C 
      DOUBLE PRECISION  F, P(500), DIVIS
      DOUBLE PRECISION  DELP, MUL1, MUL2
C 
      LOGICAL TEST
C 
      INTEGER I, J
C 
C -----------------------------------------------------------------------
      P(1)    = PRE
      CEQ     = DSQRT (A1*C(1)*C(1) + A2*C(2)*C(2))
C 
       F       = CEQ - RAC
C 
      IF (P(1) .EQ. 0.D0) THEN
        ACCP =  F /K1
      ELSE
        DIVIS   = A1*A1*K1*C(1)*C(1)+A2*A2*K2*C(2)*C(2)
        DIVIS   = DIVIS/( CEQ*CEQ)+ BE*AL/( P(1)**(1.D0-AL) )
        ACCP    = F/DIVIS
      END IF
C 
      CP1 = C(1)
      CP2 = C(2)
C 
      J    = 0
C 
1     CONTINUE
C 
      P(1)  = PRE +ACCP
C 
      IF (P(1) .LT. 0.D0) THEN
        ACCP = ACCP/2.D0
        GOTO1
      END IF
C 
      CALL  SEUIL (R0, AL, BE, P(1), RAC)
C 
      MUL1 = (1.D0 -A1*K1*ACCP/CEQ )
      MUL2 = (1.D0 -K2*A2*ACCP/CEQ )
C 
      IF (MUL1 .LT. 0.D0 .OR. MUL2 .LT. 0.D0) THEN
        ACCP = ACCP /2.D0
        GOTO  1
      END IF
C 
      C(1)  = CP1*MUL1
      C(2)  = CP2*MUL2
C 
      CEQ   =  DSQRT (A1*C(1)*C(1) + A2*C(2)*C(2))
C 
      F = CEQ - RAC
C 
      IF (F .LT. 0.D0) THEN
        J = J+1
        ACCP = ACCP/10.D0
        GOTO 1
      END IF
C 
      TEST = .FALSE.
C 
      IF (DABS(F) .LT. 1.D -6) TEST = .TRUE.
C 
C     Calcul de l'accroissement de deformation effective
C 
      CP1 = C(1)
      CP2 = C(2)
C 
      J    = 0
C 
      I = 1
C 
      DO WHILE (.NOT. TEST .AND. (I .LT. 500))
C 
        I    = I+1
C 
        DIVIS   = A1*A1*K1*C(1)*C(1)+A2*A2*K2*C(2)*C(2)
        DIVIS   = DIVIS/( CEQ*CEQ)+ BE*AL/( P(I-1)**(1.D0-AL) )
        ACCP    = F/DIVIS
C 
2       CONTINUE
C 
        P(I)  = P(I-1) + ACCP
        DELP  = P(I)   - PRE
C 
        CALL  SEUIL (R0, AL, BE, P(I), RAC)
C 
        MUL1 = (1.D0 -A1*K1*ACCP/CEQ)
        MUL2 = (1.D0 -K2*A2*ACCP/CEQ)
C 
        IF (MUL1 .LT. 0.D0 .OR. MUL2 .LT. 0.D0) THEN
          ACCP = ACCP /2.D0
          GOTO  2
        END IF
C 
        C(1)  =  CP1*MUL1
        C(2)  =  CP2*MUL2
C 
        CEQ   =  DSQRT (A1*C(1)*C(1) + A2*C(2)*C(2))
C 
        F = CEQ - RAC
C 
        IF (F .LT. 0.D0) THEN
          J = J+1
          ACCP = ACCP/2.D0
          GOTO 2
        END IF
C 
        IF (DABS(F) .LT. 1.D -6) TEST = .TRUE.
C 
        CP1 = C(1)
        CP2 = C(2)
C 
      END DO
C 
      DELP = P(I) - PRE
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... N0 numero d'interface
C 
C     On recupere :
C 
C     S ...... Y0      valeur du seuil
C     S ...... YC      equivalent a GIC
C     S ...... YR      seuil fragile a l'arrachement (mode I)
C     S ...... ALP     puissance pour le couplage
C     S ...... MINT    puissance de la loi d'evolution
C     S ...... GAM1    valeur du couplage d'endommagement 1
C     S ...... GAM2    valeur du couplage d'endommagement 2
C     S ...... A       valeur du couplage energie de cisaillement et energie de compression
C     S ...... E10     valeur initiale du module E1
C     S ...... E20     valeur initiale du module E2
C     S ...... E30     valeur initiale du module E3
C     S ...... SORT    valeur du tableau des souplesses dans la base d'orthotropie
C     S ...... SRI13   valeur limite dans le critere quadratique elastique
C     S ...... SRI23   valeur limite dans le critere quadratique elastique
C 
      SUBROUTINE CONLII (N0, R0, AL, BE, A1, A2,
     &                   K, N, Y0, YC, YR, ALP, MINT, G1, G2,
     &                   A, K1, K2, K0, SORT, SRI13, SRI23)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER N0
C 
      DOUBLE PRECISION    R0 , AL , BE , A1 , A2 , K1 , K2 , K0
      DOUBLE PRECISION    Y0, ALP , MINT, SRI13, SRI23
      DOUBLE PRECISION    K , N , YC, YR, G1, G2, A, SORT(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER M2, P, R, MNLIN, ADCAEN, M0, ADSOUP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CONLII')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('SOIN-ORTHO', ADSOUP)
      CALL ADTBDM ('HOIN-ORTHO', R)
      CALL ADTBM  ('TYP-COUCHE', M2)
      CALL ADTBDM ('CARAI-NONL', ADCAEN)
C 
C     P caracterise le type de comportement
C 
      P=M(M2+NBCOU+N0-1)
C 
CD    CALL IMPEP ('VALEUR DU TYPE DE COMPORTEMENT ', P)
C 
      M0 = R+3*(P-1)
      K1    = DM(M0)
      K2    = DM(M0+1)
      K0    = DM(M0+2)
C 
      M0      = ADSOUP+3*(P-1)
      SORT(1) = DM(M0)
      SORT(2) = DM(M0+1)
      SORT(3) = DM(M0+2)
C 
      MNLIN = ADCAEN+17*(P-1)
C 
      R0    =  DM(MNLIN)
      AL    =  DM(MNLIN+1)
      BE    =  DM(MNLIN+2)
      A1    =  DM(MNLIN+3)
      A2    =  DM(MNLIN+4)
      K     =  DM(MNLIN+5)
      N     =  DM(MNLIN+6)
      Y0    =  DM(MNLIN+7)
      YC    =  DM(MNLIN+8)
      YR    =  DM(MNLIN+9)
      ALP   =  DM(MNLIN+10)
      MINT  =  DM(MNLIN+11)
      G1    =  DM(MNLIN+12)
      G2    =  DM(MNLIN+13)
      A     =  DM(MNLIN+14)
      SRI13 =  DM(MNLIN+15)
      SRI23 =  DM(MNLIN+16)
C 
CD    CALL IMPET ('VALEUR DU NUMERO D''INTERFACE N0 ', N0)
CD    CALL IMPET ('VALEUR DE LA IERE ADRESSE M0     ', M0)
CD    CALL IMPET ('VALEUR DE LA IERE ADRESSE MNLIN  ', MNLIN)
C 
CD    CALL IMPDT ('VALEUR DE K1                     ', K1)
CD    CALL IMPDT ('VALEUR DE K2                     ', K2)
CD    CALL IMPDT ('VALEUR DE K0                     ', K0)
CD    CALL IMPDT ('VALEUR DE R0                     ', R0)      
CD    CALL IMPDT ('VALEUR DE AL                     ', AL)
CD    CALL IMPDT ('VALEUR DE BE                     ', BE)      
CD    CALL IMPDT ('VALEUR DE A1                     ', A1)           
CD    CALL IMPDT ('VALEUR DE A2                     ', A2)           
CD    CALL IMPDT ('VALEUR DE K                      ', K)           
CD    CALL IMPDT ('VALEUR DE N                      ', N)           
CD    CALL IMPDT ('VALEUR DE Y0                     ', Y0)           
CD    CALL IMPDT ('VALEUR DE YC                     ', YC)
CD    CALL IMPDT ('VALEUR DE YR                     ', YR)
CD    CALL IMPDT ('VALEUR DE ALP                    ', ALP)
CD    CALL IMPDT ('VALEUR DE MINT                   ', MINT)
CD    CALL IMPDT ('VALEUR DE G1                     ', G1)
CD    CALL IMPDT ('VALEUR DE G2                     ', G2)
CD    CALL IMPDT ('VALEUR DE A                      ', A)      
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     On envoie comme arguments :
C 
C     E ...... K(9)     comportement travaillant avec les sauts
C     E ...... SAUT     Valeur des sauts
C 
C     Et on recupere :
C 
C     S ...... SIGMAN   Valeur des contraintes normales
C 
      SUBROUTINE MCORIN (K, SAUT, SIGMAN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  K(3), SAUT(3), SIGMAN(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MCORIN')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      SIGMAN(1) = K(1) * SAUT(1)
      SIGMAN(2) = K(2) * SAUT(2)
      SIGMAN(3) = K(3) * SAUT(3)
C 
CD    CALL OMPTDP ('SIGMA NORMAL ', SIGMAN(1), 1, 3)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
