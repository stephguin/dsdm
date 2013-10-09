C     Cette routine permet de definir les caracteristiques en temps du
C     chargement, les caracteristiques pour un decoupage en plusieurs
C     periodes, ...
C 
      SUBROUTINE DONTEM
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
      INTEGER  ACARTE, NBZONT, IMPOSE
      INTEGER  NLU, ADFOTE, NMAX, I
      INTEGER  AM2LC, ADM2LC
      INTEGER  ADINTE, DEVALT, FOTEMP
      INTEGER  LECINT, NFT
C 
      DOUBLE PRECISION LECDBL
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      LOGICAL  OK
C 
      CHARACTER*40 MOT
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DONTEM')
C 
CD    CALL WLKBCD (IDPROG)
C 
C     Dans DONTEM, lecture des caracterisiques en temps du chargement
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
      CALL LECSEQ ('DONTEM', 'CARACTERISTIQUES EN TEMPS ')
C 
      DUREE  = LECDBL ('DUREE')
      NBFODO = LECINT ('NBFODO')
C 
C     Dans CALTEM, lecture des caracteristiques pour le decoupage en temps
C 
      NFOTPS = NBFODO
      NFTGLO = NBFODO
C 
      DUREAC = DUREE
C 
      NPICMX = NPICET
C 
      IF (.NOT. STANDA) THEN
        CALL LECSEQ ('CALTEM', 'POUR PLUSIEURS DECOUPAGES EN TEMPS')
      ELSE
        CALL LFORCS ('CALTEM,,,,,,', 12, MOT, OK)
      END IF
C 
      IMPOSE = LECINT ('IMPOSE')
      NBZONT = LECINT ('NBZONT')
C 
C -----------------------------------------------------------------------
C     Dans le tableau CARA_TEMPS sont ranges :
C     si le decoupage en temps est impose => IMPOSE = 1
C     alors le nombre de decoupage en temps
C     le numero decoupage actuel depuis la reprise
C     le numero de decoupage actuel
C     le nombre de piquets de temps des differents decoupage
C -----------------------------------------------------------------------
      CALL GESTEN ('CARA_TEMPS', 5+NBZONT, ACARTE)
      M(ACARTE)   = IMPOSE
      M(ACARTE+1) = NBZONT
      M(ACARTE+2) = 1
      M(ACARTE+3) = 1
      M(ACARTE+4) = NPICET
C 
      NMAX = 4
      CALL GESTDP ('TEMPS-DONS', NMAX*NBFODO, ADFOTE)
C 
C -----------------------------------------------------------------------
C     Pour boucler sur les fonctions du temps
C -----------------------------------------------------------------------
      CALL MESSAO (
     &' ON CARACTERISE LES FONCTIONS DU TEMPS PAR 4 DONNEES DOUBLES:
     &\    1ERE DONNEE (CORRESPOND AU TYPE) :
     &\    0.=> Polynome du 3eme degre en fonction de (t/duree)
     &\         nul a t= 0  ==> 3 DONNEES
     &\    1.=> montee descente  ==> 3 DONNEES
     &\         valeur maxi, temps correspondant, temps d''annulation')
C  
12    CONTINUE
C 
      CALL LECLDP (DM(ADFOTE), NMAX, NLU)
C  
      IF (NLU .NE. NMAX) THEN
C 
        CALL MESSAO (
     &   'VOUS N''AVEZ PAS RENTRE LE BON NOMBRE DE VALEURS 
     &   \POUR LES FONCTIONS DU TEMPS RECOMMENCEZ')
          GOTO 12
      END IF
C  
      CALL GESTDP ('LONG_INTER', NBZONT+1, ACARTE)
C 
      DM(ACARTE)   = 0.D0
      DM(ACARTE+1) = DUREE
C 
C -----------------------------------------------------------------------
C     Pour boucler sur les periodes en temps
C -----------------------------------------------------------------------
1     CONTINUE
C 
      IF (.NOT. STANDA) THEN
        CALL MESSAO ('ON DONNE LES VALEURS DES DIFFERENTS TEMPS '//
     &               'EN FIN D''INTERVALLE')
        CALL LECLDP (DM(ACARTE+2), NBZONT, NLU)
      END IF
C  
      IF (NLU .NE. (NBZONT-1)) THEN
         CALL IMPET ('NOMBRE DE ZONES EN TEMPS ', NBZONT)
         CALL IMPET ('INCOHERENCE : NLU+1 DIFFERENT DE NBZONT, '//
     &               'RECOMMENCEZ ', NLU+1)
         GOTO 1
      ENDIF
C  
      DO  I = 2, NLU +1
        IF (DM(ACARTE+I) .LE. DM(ACARTE+I-1)) THEN
          CALL ERREUD (0, 'TEMPS2 < TEMPS1 '//IDPROG)
        END IF
      END DO
C 
      CALL GESTDP ('INTE-TEMPS', NPICET, ADINTE)
      CALL GESTDP ('VALE-TEMPS', NPICET+1, DEVALT)
      CALL GESTDP ('F-TEMPS-DO', NPICET, FOTEMP)
C 
      DO NFT = 1, NBFODO
        CALL IMPET ('FONCTION DU TEMPS NUMERO ', NFT)
        CALL IMPTDT ('FONCTIONS DU TEMPS ', DM(ADFOTE), 4, 1)
        ADFOTE = ADFOTE+4
      END DO
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
C     Cette routine calcule la 1ere discretisation en temps
C 
      SUBROUTINE NDISTM
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
C     Pour les couches
C 
      INTEGER          I, AM2LC, ADM2LC
C 
      DOUBLE PRECISION  INTERV
C 
C     Pour les temps
C 
      INTEGER ADINTE, DEINTE, DEVALT, ADFOTE, DEFOTE
      INTEGER DETEMP, FOTEMP
C 
      DOUBLE PRECISION VALFIN
C 
C     Adresse des caracteristiques dans la base d'orthotropie
C     Pour les interfaces
C 
      INTEGER  VAFTLO, NPICLO
C 
C     Pour ne pas faire doublon avec les tableaux
C 
      INTEGER  ACARTE, AINTER, NUCONV, ADTAU2, LONINT
      INTEGER  ADERR, LONVAL, LONFOT, LONCAR, LONGIN
C 
C     Pour le dessin
C 
      CHARACTER*6   IDPROG
C 
      PARAMETER (IDPROG='NDISTM')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL INFODP ('INTE-TEMPS', ADINTE, LONINT)
      DEINTE = ADINTE-1
C 
      CALL INFODP ('VALE-TEMPS', DEVALT, LONVAL)
C 
      CALL INFODP ('TEMPS-DONS', ADFOTE, LONFOT)
      DEFOTE    = ADFOTE
C 
C     Definition et remplissage du tableau des accroissements des
C     fonctions du temps elastique. Stockage : npicet
C 
      CALL INFODP ('F-TEMPS-DO', FOTEMP, LONFOT)
      DETEMP    = FOTEMP-1
C 
C     Pour obtenir la discretisation voulue pour la premiere etape
C 
      CALL INFOEN ('CARA_TEMPS', ACARTE, LONCAR)
C 
      CALL IMPTET ('CARACTERISTIQUES DU TEMPS DANS '//IDPROG ,
     &            M(ACARTE), 1, LONCAR)
C 
      CALL INFODP ('LONG_INTER', AINTER , LONGIN)
C 
      CALL IMPTDN ('LONG-INTER DANS '//IDPROG ,
     &              DM(ACARTE), 1, LONCAR)
C 
        IF (M(ACARTE) .EQ . 1) THEN
          NUCONV = M(ACARTE+2)
          TEMDI1 = DM (AINTER+NUCONV)-DM(AINTER+NUCONV-1)
          CALL IMPDN ('TEMDI1 DANS '//IDPROG, TEMDI1)
        ELSE
          TEMDI1 = 1.5D0*DMIN1(DM(DEVALT+NPICET), DM(VAFTLO+NPICLO))
          TEMDI1 = DMIN1 (DUREE, TEMDI1)
        END IF
C 
        NPICET = NPICMX
C 
        CALL IMPMN ('ON IMPOSE LA DISCRETISATION ')
        INTERV = TEMDI1/NPICMX
C 
C       REMPLISSAGE DES FONCTIONS DU TEMPS
C 
        DM(DEVALT) = 0.D0
C 
        DO I = 1, NPICMX
          DM(ADINTE+I-1)   = INTERV
          DM(DEVALT+I)     = DM( DEVALT+I-1)+INTERV
        END DO
C 
      CALL GESTDP ('PARAM-CHAR', NPICET+1, ADERR)
      DM(ADERR) = 0.D0
C 
        VALFIN = 0.D0
C 
        DO I = 1, NPICMX
          CALL VAFTEL (DM(ADFOTE), DM(DEVALT+I-1), DM(DEVALT+I),
     &                 INTERV, DM(FOTEMP+I-1))
          DM(ADERR+I) = DM(ADERR+I-1)+DM(FOTEMP+I-1)
          VALFIN = VALFIN + DM(FOTEMP+I-1)
        END DO
C 
      CALL IMPTDT
     & ('Parametre de charge dans '//IDPROG,
     &   DM(ADERR), NPICET+1, 1)
C 
      CALL IMPDT
     &     ('**** VALEUR FINALE DE LA FONCTION DU TEMPS *** ', VALFIN)
C 
      CALL GESTDP ('TAUX-CHA-2', NPICET, ADTAU2)
C 
      FOTEMP = FOTEMP-1
      ADTAU2 = ADTAU2-1
C 
      DO I = 1 , NPICET
        DM(ADTAU2+I) =  DM(FOTEMP+I)*DM(FOTEMP+I)
      END DO
C 
      CALL IMPTDT 
     & ('TABLEAU TAUX-CHA-2 '//IDPROG,
     &   DM(ADTAU2+1), NPICET, 1)
C 
C --------------------------------------------------------------
C     Fin de la sequence due aux interfaces
C --------------------------------------------------------------
C 
C     On definit la duree actuelle prise en compte dans le calcul
C 
      DUREAC = DM(DEVALT+NPICET)
      CALL IMPDT ('DUREAC DANS '//IDPROG, DUREAC)
C 
      NPICAC = NPICET
C 
C     Duree de l'intervalle de temps pris en compte
C 
      LTEMPS = DM(DEVALT+NPICET)-DM(DEVALT)
      CALL IMPDT ('LTEMPS DANS '//IDPROG, DUREAC)
C 
C     Piquets de temps de reprise mis a npicet
C 
      PICREP = NPICET
C 
      CALL IMPTDT ('VALEURS DU TEMPS DANS '//IDPROG,
     &              DM(DEVALT), NPICET+1, 1)
C 
      CALL IMPTDT ('INTERVALLES DE TEMPS DANS '//IDPROG,
     &              DM(ADINTE), NPICET, 1)
C 
      CALL IMPTDT ('ACCROISSEMENTS FONCTIONS DU '//
     &             'TEMPS DONNEES (NPICET) DANS '//IDPROG,
     &              DM(FOTEMP+1), NPICET, 1)
C 
      CALL PRETEM
C 
C     Pour le calcul de l'erreur totale
C     TR((siga-sigc)*(epsa-epsc)) + ((signa-signc)*(saua-sauc))
C     On construit un denominateur pour les couches (DE-COU-LOC) et un denominateur
C     pour les interfaces, afin de pouvoir visualiser les contributions relatives
C     couches/interfaces locales a l'indicateur d'erreur.
C 
      CALL GESTDP ('NUMERA-LOC', NPICET+1, ADERR)
      DM(ADERR) = 0.D0
C 
      CALL GESTDP ('DENOMI-LOC', NPICET+1, ADERR)
      DM(ADERR) = 0.D0
C 
      CALL GESTDP ('DE-COU-LOC', NPICET+1, ADERR)
      DM(ADERR) = 0.D0
C 
      CALL GESTDP ('DE-INT-LOC', NPICET+1, ADERR)
      DM(ADERR) = 0.D0
C 
      CALL GESTDP ('ERREUR-LOC', NPICET+1, ADERR)
      DM(ADERR) = 0.D0
C 
C     Pour la 1ere etape locale ou aucune erreur n'a ete
C     calculee on met la valeur precedente a "l'infini"
C 
      DO I = 1, NPICET
       DM(ADERR+I) = 10000.D0
      END DO
C 
C     Dans info-erreur on garde la valeur precedente de l'erreur
C     + le nombre de fonctions du temps utilisee=> 2
C 
      CALL GESTDP ('INF-ERREUR', NPICET+3, ADERR)
      DM(ADERR) = 0.D0
C 
      DO I = 1, NPICET
       DM(ADERR+I) = 10000.D0
      END DO
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Valeur de l'Accroissement de la Fonction du Temps ELastique
C 
C     Cette routine calcule l'accroissement de l'evolution elastique
C     entre les temps T1 et T2 (T2 > T1).
C 
C     On envoie comme arguments :
C 
C     E ...... T1       valeur du temps initial
C     E ...... INTERV   T2-T1
C     E ...... T2       valeur du temps final
C     E ...... NBVALE   nombre de valeurs a partir desquelles sont estimees
C     E ...... CARTEM   constantes caracteristiques de l'evolution
C     E                 f(t) = C1*t +C2*(t**2)
C     Et on recupere :
C 
C     S ...... ACCFEL   fel(t2) - fel(t1)
C 
      SUBROUTINE VAFTEL (CARTEM, T1, T2, INTERV, ACCFEL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  T1 , T2 , CARTEM(4) , INTERV
      DOUBLE PRECISION  ACCFEL
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  T1LOC , T2LOC
C 
CD    LOGICAL           LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VAFTEL')
C 
CD    CALL WLKBCD(IDPROG)
C 
C     CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
CD    IF ( T1 .GT . T2 ) THEN
CD        CALL IMPDT( ' T1 ' , T1 )
CD        CALL IMPDT( ' T2 ' , T2 )
CD       CALL ERREUD(0, 'T1 > T2 DANS '//IDPROG)
CD    END IF
C 
C     RAPPEL : On caracterise les fonctions du temps par 4 DONNEES DOUBLES :
C 
C     - 1er donnee (correspond au type) :
C     - 0.=> Polynome du 3eme degre !! en fonction de (t/duree)
C     - nul t= 0  ==> 3 DONNEES
C     - 1.=> montee descente  ==> 3 DONNEES
C     - valeur maxi, temps correspondant, temps d'annulation ')
C 
C     On calcule les acroissements de la solution elastique
C 
      T1LOC = T1/DUREE
      T2LOC = T2/DUREE
C 
      IF (CARTEM(1) .EQ. 0.D0) THEN
C 
C     TYPE = POLYNOME EN TEMPS/DUREE
C 
       ACCFEL = (T2LOC-T1LOC)*(CARTEM(2)
     &           + CARTEM(3)*(T2LOC+T1LOC)
     &           + CARTEM(4)*(T2LOC*T2LOC+T1LOC*T1LOC
     &           + T2LOC*T1LOC))
C 
      ELSE IF (CARTEM(1) .EQ .1.D0) THEN
C 
C     TYPE = MONTEE_DESCENTE
C 
        IF (T2 .LE. CARTEM(3)) THEN
C 
C       ON EST DANS LA PHASE MONTEE
C 
          ACCFEL = CARTEM(2)*(T2 -T1)/CARTEM(3)
C 
        ELSE IF( T1 .GE. CARTEM(3) )THEN
C 
C       ON EST DANS LA PHASE DESCENTE
C 
          ACCFEL = CARTEM(2)*(T2 -T1)/(CARTEM(3)-CARTEM(4) )
C 
        ELSE
C 
          CALL MESSAO ('ON EST ENTRE MONTEE ET DESCENTE '//idprog)
C 
          CALL IMPDT ('T1           ', T1)
          CALL IMPDT ('T2           ', T2)
          CALL IMPDT ('T VALEUR MAXI', CARTEM(3))
C 
          ACCFEL =CARTEM(2)*(T2-CARTEM(4))/( CARTEM(3)-CARTEM(4) )
     &            -CARTEM(2)*T1/CARTEM(3)
C 
        END IF
C 
      ELSE
C 
        CALL IMPDT( ' CARTEM(1) = TYPE DIFF DE )> OU 1.' , CARTEM(1) )
        CALL ERREUD(0,' MAUVAIS PASSAGE DE TYPE DANS '//IDPROG )
C 
      END IF
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1)) THEN
CD     CALL IMPDN ('T1       ', T1)
CD     CALL IMPDN ('T2       ', T2)
CD     CALL OMPTDN ('CARTEM   ', CARTEM, 4, 1)
CD     CALL IMPDN ('ACCFEL   ', ACCFEL)
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
      SUBROUTINE APPART( B1 , B2 , VAL , APP )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      DOUBLE PRECISION  B1 , B2 , VAL
C 
      LOGICAL APP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='APPART')
C 
      DOUBLE PRECISION  INFB , SUPB
C 
CD    CALL WLKBCD(IDPROG)
C 
      INFB = MIN( B1 , B2 )
      SUPB = MAX( B1 , B2 )
C 
      APP = .TRUE .
C 
      IF( VAL . GT. SUPB .OR. VAL .LT.INFB ) THEN
        APP = .FALSE.
      END IF
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
