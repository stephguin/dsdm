C     Cette routine calcule les variables d'endommagement.
C 
C     Steph le 05/11/99 : on travaille desormais en contraintes.
C     YD1 et YD sont toujours explicites en deformations, donc en temps.
C     YDP n'est plus explicite en temps. On itere en precision sur les valeurs
C     de d1, d, dp dans la boucle j de REPELA, sous-iteration de la boucle i de
C     partition endommagement/plasticite dans INCOCO.
C 
C     L'activation de la rupture en traction transverse n'est pas explicite.
C 
C     On envoie comme arguments :
C 
C     E ...... YC, YCP, K, KP   |
C     E ...... Y0, Y0P, B, BP   | parametres d'endommagement
C     E ...... N, NP, YCS, YTS  |
C     E ...... D1N              | variables d'endommagement
C     E ...... DN               | au pas de temps precedent
C     E ...... DPN              | 
C     E ...... YD1              | forces d'endommagement
C     E ...... YD               | au pas de temps courant
C     E ...... YDP              | 
C     E ...... TESTD1           | logiques vrais si la rupture est
C     E ...... TESTD            | activee au pas de temps precedent
C     E ...... TESTDP           |
C 
C     Et en Entree-Sortie :
C 
C     ES ..... LOGRUP1            logique vrai si rup. activee a l'iter. precedente
C     ES ..... LOGRUP2            logique vrai si rup. activee a l'iter. courante
C                             
C     Et on recupere :
C                             
C     E ...... D(4)             | variables d'endommagement au pas de temps courant
C 
C 
C     Formulation de YD et YDP pour une evolution lente:
C     -------------------------------------------------
C      (t+1)    0                            2     0                           2
C     Y  =    G    <<EPILON(1,2) ELASTIQUE>>  +  G   <<EPILON(1,3) ELASTIQUE>>
C      d       12                                 13
C      (t+1)                     2             2                    2             2
C     Y   =        <<SIGMA(2,2)>> / (2 E (1-d') )  +  <<SIGMA(3,3)>> / (2 E (1-d') )
C      d'                               2                                  3
C 
C     Lois d'evolutions lentes pour d et d', resolution de l'equadiff :
C     ----------------------------------------------------------------
C     .(t+1)              (t+1)     (t+1)  1/2    1/2     1/2    (t+1)   N
C     d      = k   < ( ( Y     + b Y      )   - Y   ) / Y     - d      >       si d < 1
C                         d         d'          0       c               +
C 
C     .(t+1)            (t+1)       (t+1)  1/2   1/2     1/2    (t+1)   N'
C     d'     = k'  < ( Y       + b'Y      )   - Y'  )/ Y'    - d'      >       si d'< 1
C                       d'           d          0      c                +
C 
C     Formulation de YD1 et YDP pour une evolution rapide:
C     ---------------------------------------------------
C      (t+1)
C     Y  =  10  [<<EPILON(1,1) ELASTIQUE>>  -  EPTLIM] / EPTLIM     en traction
C      d1
C 
C      (t+1)
C     Y  =  10  [<<EPILON(1,1) ELASTIQUE>>  -  EPCLIM] / EPCLIM     en compression
C      d1
C 
C      (t+1)      (t+1)
C     Y  =  10  [Y  - Y' ] / Y'
C      d'         d'   ts     ts
C 
C     Lois d'evolutions rapides pour d1 et d', resolution de l'equadiff :
C     ------------------------------------------------------------------
C     .(t+1)         (t+1)   (t+1)
C     d1    = kf  < Y     - d1    >   si d < 1
C                    d1             +      1
C     .(t+1)         (t+1)   (t+1)
C     d'    = kf  < Y     - d'    >   si d'< 1
C                    d'             +
C 
      SUBROUTINE VAENDO (DELTAT, COENDO, DN, YD1, YD, YDP,
     &                   NOCVD, D, RUPCOU)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'typcal.h'
      LOGICAL  RUPCOU(2)
      LOGICAL  NOCVD(3)
C 
      DOUBLE PRECISION  COENDO(12)
      DOUBLE PRECISION  DELTAT, DN(3), YD1, YD, YDP, D(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  B, K, Y0, YC, YCS, N, BP, KP, Y0P, YCP, YTS, NP
      DOUBLE PRECISION  YEQ, YEQP, YFL, YDPLOC
      DOUBLE PRECISION  KF
C 
C     KF est la 'vitesse' d'evolution pour les variables rapides
C 
      PARAMETER         (KF = 10000.D0)
      CHARACTER*6        IDPROG
      PARAMETER         (IDPROG = 'VAENDO')
C 
C -----------------------------------------------------------------------
      IF ((.NOT. LENDCO) .AND. (.NOT. LRUPCO)) THEN
        D(1) = 0.D0
	D(2) = 0.D0
	D(3) = 0.D0
	GOTO 4
      END IF
C 
      B   = COENDO(1)
      K   = COENDO(2)
      Y0  = COENDO(3)
      YC  = COENDO(4)
      YCS = COENDO(5)
      N   = COENDO(6)
      BP  = COENDO(7)
      KP  = COENDO(8)
      Y0P = COENDO(9)
      YCP = COENDO(10)
      YTS = COENDO(11)
      NP  = COENDO(12)
C 
C     INITIALISATION DES VALEURS DE D1, D ET DP
C 
      D(1) = DN(1)
      D(2) = DN(2)
      D(3) = DN(3)
C 
C     CALCUL DE Df (endommagement des fibres)
C 
      IF (NOCVD(1) .OR. (.NOT. LRUPCO)) GOTO 2
C 
      CALL INTRET (DELTAT, KF, 1.D0, YD1, DN(1), D(1), 1.D0)
2     CONTINUE
C 
C     CALCUL DE Dps
C 
      IF (NOCVD(2)) GOTO 3
C 
      YEQ = DMAX1((DSQRT(YD + B*YDP)-Y0), 0.D0)/YC
      CALL INTRET (DELTAT, K, N, YEQ, DN(2), D(2), 1.D0)
3     CONTINUE
C 
C     CALCUL DE Dpt
C 
      IF (NOCVD(3)) GOTO 4
C 
      IF ((.NOT. RUPCOU(1)) .AND. (.NOT. RUPCOU(2))) THEN
        YFL  = DSQRT(YDP + BP*YD)-YTS
	IF (YFL .GT. 0.D0) THEN
	  RUPCOU(1) = .TRUE.
	  YDPLOC  = YTS*YTS
	  D(3) = DMAX1((DSQRT(YDPLOC + BP*YD)-Y0P), 0.D0)/YCP
          D(3) = DMIN1 (D(3), 1.D0)
        ELSE
	  YEQP = DMAX1((DSQRT(YDP + BP*YD)-Y0P), 0.D0)/YCP
   	  YEQP = DMIN1(YEQP, 1.D0) 
	  CALL INTRET (DELTAT, KP, NP, YEQP, DN(3), D(3), 1.D0)
	END IF
	GOTO 4
      END IF
C 
      IF (RUPCOU(1) .AND. (.NOT. RUPCOU(2))) THEN
        YFL  = DSQRT(YDP + BP*YD)-YTS
        IF (YFL .GT. 0.D0) THEN
	  YDPLOC = 10.D0*DMAX1(YDP-YTS, 0.D0)/YTS
          YDPLOC = DMIN1 (YDP, 1.D0)
	  RUPCOU(2)=.TRUE.
   	  CALL INTRET (DELTAT, KF, 1.D0, YDPLOC, DN(3), D(3), 1.D0)
	ELSE
	  RUPCOU(1) = .FALSE.
	  RUPCOU(2)= .TRUE.
          YEQP = DMAX1((DSQRT(YDP + BP*YD)-Y0P), 0.D0)/YCP
          YEQP = DMIN1 (YEQP, 1.D0)
   	  CALL INTRET (DELTAT, KP, 1.D0, YEQP, DN(3), D(3), 1.D0)
	END IF
	GOTO 4
      END IF
C 
      IF (RUPCOU(1) .AND. RUPCOU(2)) THEN
	YDPLOC = 10.D0*DMAX1(YDP-YTS, 0.D0)/YTS
        YDPLOC = DMIN1 (YDP, 1.D0)
   	CALL INTRET (DELTAT, KF, 1.D0, YDPLOC, DN(3), D(3), 1.D0)
	GOTO 4
      END IF
C 
      IF ((.NOT. RUPCOU(1)) .AND. RUPCOU(2)) THEN
CD 	CALL MESSAO ('CAS A LA CON ACTIVE SUR D(3) !!!!!!!!')
	YEQP = DMAX1((DSQRT(YDP + BP*YD)-Y0P), 0.D0)/YCP
        YEQP = DMIN1 (YEQP, 1.D0)
   	CALL INTRET (DELTAT, KP, NP, YEQP, DN(3), D(3), 1.D0)
	GOTO 4
      END IF
C 
4     CONTINUE
C 
      D(1) = DMAX1(D(1), DN(1))
      D(2) = DMAX1(D(2), DN(2))
      D(3) = DMAX1(D(3), DN(3))
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CASITR calcule le champ de contraintes associe au champ de deformation
C     envoye en entree. Desormais, les modules elastiques sont actualises ici, pour
C     traiter les problemes lies aux retours en compression. Avec la formulation en
C     contraintes, ces calculs ne sont plus explicites.
C 
C     On envoie comme arguments :
C 
C     E ...... D2D3       option calcul contraintes planes ou contraintes HP
C     E ...... ELAIN(6)   coefficients elastiques vierges
C     E ...... COFIBR(3)  coefficients pour le comportement non
C     E ...... ENDOMS(3)  endommagement d1, d, d' et d''
C     E ...... LOGIN      logiques ranges E11, E22, E33 : l'endommagement est actif
C                         en traction, inactif en compression
C     E ...... DE         deformations elastiques dans la base d'orthotropie
C 
C     Et on recupere :
C 
C     S ...... C          contraintes chapeau dans la base dans la base d'orthotropie 
C     S ...... EAC        modules ranges E11, E22, E33 en sortie
C     S ...... LOGOUT     logiques ranges E11, E22, E33 en sortie
C 
      SUBROUTINE CASITR (D2D3, ELAIN, COFIBR, ENDOMS,
     &                   LOGIN, DE, C, EAC, LOGOUT, LOGSIT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
C     S Souplesse initiale
C 
      DOUBLE PRECISION D2D3, ELAIN(9), COFIBR(3)
      DOUBLE PRECISION DE(6), ENDOMS(3), EAC(3), C(6)
      LOGICAL          LOGSIT, LOGIN(3), LOGOUT(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION C2(6), C3(6), C4(6), C5(6), C6(6)
      DOUBLE PRECISION C7(6), C8(6), ELAOUT(9)
      LOGICAL          TEST00, TEST01, TEST02, TEST03, TEST04, LOG1(3)
      LOGICAL          TEST05, TEST06, TEST07, TEST08
      CHARACTER*6      IDPROG
      PARAMETER       (IDPROG='CASITR')
C 
C -----------------------------------------------------------------------
      TEST00 = .FALSE.
      TEST01 = .FALSE.
      TEST02 = .FALSE.
      TEST03 = .FALSE.
      TEST04 = .FALSE.
      TEST05 = .FALSE.
      TEST06 = .FALSE.
      TEST07 = .FALSE.
      TEST08 = .FALSE.
C 
C -----------------------------------------------------------------------
C     Cas 00 : tout se passe comme au pas de temps precedent
C     Sinon, on envisagera tous les cas possibles (cas 1 a 8)
C -----------------------------------------------------------------------
C 
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN,
     &             ELAOUT, C)
      CALL TESCAS (C, LOGOUT)
      IF ((LOGIN(1) .EQV. LOGOUT(1)) .AND. (LOGIN(2) .EQV. LOGOUT(2))
     &    .AND. (LOGIN(3) .EQV. LOGOUT(3))) THEN
        TEST00 = .TRUE.
        IF (LOGSIT) THEN
          CALL IMPTDT ('VERIF 1 CASITR ', C, 1, 6)
        END IF
	GOTO 1
      END IF
C 
C -----------------------------------------------------------------------
C     Cas 01 : E1=E10(1-D1), E2=E20(1-DP), E3=E30(1-DP)
C              SIG11>0,      SIG22>0,      SIG33>0
C -----------------------------------------------------------------------
C 
      LOG1(1) = .TRUE.
      LOG1(2) = .TRUE.
      LOG1(3) = .TRUE.
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN,
     &             ELAOUT, C)
      CALL TESCAS (C, LOGOUT)
      IF ((LOG1(1) .EQV. LOGOUT(1)) .AND. (LOG1(2) .EQV. LOGOUT(2))
     &    .AND. (LOG1(3) .EQV. LOGOUT(3))) THEN
        TEST01 = .TRUE.
	IF (LOGSIT) THEN
	  CALL IMPTDT ('VERIF 0 CASITR ', C, 1, 6)
	END IF
	GOTO 1
      END IF
C 
C -----------------------------------------------------------------------
C     Steph. 15/09/99. Je vire le cote compression fibres qui fait chier
C     Cas 02 : E1=E10 , E2=E20(1-DP), E3=E30(1-DP),
C              SIG11<0, SIG22>0,      SIG33>0
C -----------------------------------------------------------------------
C 
      LOG1(1) = .FALSE.
      LOG1(2) = .TRUE.
      LOG1(3) = .TRUE.
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN,
     &             ELAOUT, C2)
      CALL TESCAS (C2, LOGOUT)
      IF ((LOG1(1) .EQV. LOGOUT(1)) .AND. (LOG1(2) .EQV. LOGOUT(2))
     &    .AND. (LOG1(3) .EQV. LOGOUT(3))) TEST02 = .TRUE.
      IF (TEST02 .AND. LOGSIT) THEN
        CALL IMPTDT ('VERIF 2 CASITR ', C2, 1, 6)
      END IF
C 
C -----------------------------------------------------------------------
C     Cas 03 : E1=E10(1-D1), E2=E20   , E3=E30(1-DP)
C              SIG11>0,      SIG22<0  , SIG33>0
C -----------------------------------------------------------------------
C 
      LOG1(1) = .TRUE.
      LOG1(2) = .FALSE.
      LOG1(3) = .TRUE.
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN,
     &             ELAOUT, C3)
      CALL TESCAS (C3, LOGOUT)
      IF ((LOG1(1) .EQV. LOGOUT(1)) .AND. (LOG1(2) .EQV. LOGOUT(2))
     &    .AND. (LOG1(3) .EQV. LOGOUT(3))) TEST03 = .TRUE.
      IF (TEST03 .AND. LOGSIT) THEN
        CALL IMPTDT ('VERIF 3 CASITR ', C3, 1, 6)
      END IF
C 
C -----------------------------------------------------------------------
C     Cas 04 : E1=E10(1-D1), E2=E20(1-DP), E3=E30
C              SIG11>0,      SIG22>0,      SIG33<0
C -----------------------------------------------------------------------
C 
      LOG1(1) = .TRUE.
      LOG1(2) = .TRUE.
      LOG1(3) = .FALSE.
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN, ELAOUT, C4)
      CALL TESCAS (C4, LOGOUT)
      IF ((LOG1(1) .EQV. LOGOUT(1)) .AND. (LOG1(2) .EQV. LOGOUT(2))
     &    .AND. (LOG1(3) .EQV. LOGOUT(3))) TEST04 = .TRUE.
      IF (TEST04 .AND. LOGSIT) THEN
        CALL IMPTDT ('VERIF 4 CASITR ', C4, 1, 6)
      END IF
C 
      IF (TEST02 .AND. (.NOT. TEST03) .AND. (.NOT. TEST04)) THEN
	CALL COPITD (6, C2, C)
	GOTO 1
      ELSE IF ((.NOT. TEST02) .AND. TEST03 .AND. (.NOT. TEST04)) THEN
        CALL COPITD (6, C3, C)
	GOTO 1
      ELSE IF ((.NOT. TEST02) .AND. (.NOT. TEST03) .AND. TEST04) THEN
        CALL COPITD (6, C4, C)
	GOTO 1
      END IF
C 
C -----------------------------------------------------------------------
C     Steph. 15/09/99. Je vire le cote compression fibres qui fait chier
C     Cas 05 : E1=E10 , E2=E20, E3=E30(1-DP),
C              SIG11<0, SIG22<0 et SIG33>0
C -----------------------------------------------------------------------
C 
      LOG1(1) = .FALSE.
      LOG1(2) = .FALSE.
      LOG1(3) = .TRUE.
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN, ELAOUT, C5)
      CALL TESCAS (C5, LOGOUT)
      IF (TEST05 .AND. LOGSIT) THEN
        CALL IMPTDT ('VERIF 5 CASITR ', C5, 1, 6)
      END IF
      IF ((LOG1(1) .EQV. LOGOUT(1)) .AND. (LOG1(2) .EQV. LOGOUT(2))
     &    .AND. (LOG1(3) .EQV. LOGOUT(3))) TEST05 = .TRUE.
C 	      
C -----------------------------------------------------------------------
C     Steph. 15/09/99. Je vire le cote compression fibres qui fait chier
C     Cas 06  : E1=E10 , E2=E20(1-DP), E3=E30,
C               SIG11<0, SIG22>0 et SIG33<0
C -----------------------------------------------------------------------
C 
      LOG1(1) = .FALSE.
      LOG1(2) = .TRUE.
      LOG1(3) = .FALSE.
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN, ELAOUT, C6)
      CALL TESCAS (C6, LOGOUT)
      IF ((LOG1(1) .EQV. LOGOUT(1)) .AND. (LOG1(2) .EQV. LOGOUT(2))
     &    .AND. (LOG1(3) .EQV. LOGOUT(3))) TEST06 = .TRUE.
      IF (TEST06 .AND. LOGSIT) THEN
        CALL IMPTDT ('VERIF 6 CASITR ', C6, 1, 6)
      END IF
C 	      
C -----------------------------------------------------------------------
C       Cas 07 : E1=E10(1-D1), E2=E20, E3=E30, SIG11>0, SIG22<0 et SIG33<0
C -----------------------------------------------------------------------
C 
      LOG1(1) = .TRUE.
      LOG1(2) = .FALSE.
      LOG1(3) = .FALSE.
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN, ELAOUT, C7)
      CALL TESCAS (C7, LOGOUT)
      IF ((LOG1(1) .EQV. LOGOUT(1)) .AND. (LOG1(2) .EQV. LOGOUT(2))
     &    .AND. (LOG1(3) .EQV. LOGOUT(3))) TEST07 = .TRUE.
      IF (TEST07 .AND. LOGSIT) THEN
        CALL IMPTDT ('VERIF 7 CASITR ', C7, 1, 6)
      END IF
C 
      IF (TEST05 .AND. (.NOT. TEST06) .AND. (.NOT. TEST07)) THEN
	CALL COPITD (6, C5, C)
	GOTO 1
      ELSE IF ((.NOT. TEST05) .AND. TEST06 .AND. (.NOT. TEST07)) THEN
        CALL COPITD (6, C6, C)
	GOTO 1
      ELSE IF ((.NOT. TEST05) .AND. (.NOT. TEST06) .AND. TEST07) THEN
        CALL COPITD (6, C7, C)
	GOTO 1
      END IF
C 
C -----------------------------------------------------------------------
C     Steph. 15/09/99. Je vire le cote compression fibres qui fait chier
C     Cas 08 : E1=DMAX1(E10(1+GAM*SIG11/E10), 0), E2=E20(1-DP), E3=E30,
C              SIG11<0, SIG22>0 et SIG33<0
C -----------------------------------------------------------------------
C 
      LOG1(1) = .FALSE.
      LOG1(2) = .FALSE.
      LOG1(3) = .FALSE.
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN, ELAOUT, C8)
      CALL TESCAS (C8, LOGOUT)
      IF ((LOG1(1) .EQV. LOGOUT(1)) .AND. (LOG1(2) .EQV. LOGOUT(2))
     &    .AND. (LOG1(3) .EQV. LOGOUT(3))) THEN
        TEST08 = .TRUE.
        CALL COPITD (6, C8, C)
      END IF
      IF (TEST08 .AND. LOGSIT) THEN
        CALL IMPTDT ('VERIF 8 CASITR ', C8, 1, 6)
      END IF
C 
1     CONTINUE
C 
      IF ((.NOT. TEST01) .AND. (.NOT. TEST02) .AND. (.NOT. TEST03) .AND.
     &    (.NOT. TEST04) .AND. (.NOT. TEST05) .AND. (.NOT. TEST06) .AND.
     &    (.NOT. TEST07) .AND. (.NOT. TEST08) .AND. (.NOT. TEST00)) THEN
        CALL MESSAO ('HALTE LA : PAS DE SOLUTION ??')
	STOP
      END IF
      EAC(1) = ELAOUT(1)
      EAC(2) = ELAOUT(2)
      EAC(3) = ELAOUT(3)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     ACTURI est une routine qui calcule l'operateur de rigidite endommage
C     (9 termes).
C 
C     On envoie comme arguments :
C 
C     E ...... D2D3       option de calcul : 2D-SIG ou 3D-SIG
C     E ...... ELAIN(9)   coefficients elastiques initiaux
C                         sequence (E11, E22, E33, G12, G23, G13, NU12, NU23, NU13)
C     E ...... COFIBR(3)  coefficients pour le comportement NL sens fibres
C     E ...... DEFIN      deformations
C     E ...... ENDOIN(3)  endommagement d1, d et d'
C     E ...... LOGIN(3)   tableau de logiques ranges E11, E22, E33
C                        (.TRUE., .TRUE., .TRUE., ) si D1 est active (-----> SIG11 > 0)
C                                                   si D  est active (-----> SIG22 > 0)
C                                                   si DP est active (-----> SIG33 > 0)
C     Et on recupere :
C 
C     S ...... ELAOUT(9)  coefficients elastiques endommages
C     S ...... SIOUT(6)   nouvelles contraintes
C 
      SUBROUTINE ACTURI (D2D3, ELAIN, COFIBR, DEFIN, ENDOIN,
     &                   LOGIN, ELAOUT, SIOUT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION D2D3, ELAIN(9), COFIBR(3), DEFIN(6)
      DOUBLE PRECISION ENDOIN(3)
      LOGICAL          LOGIN(3)
      DOUBLE PRECISION ELAOUT(9), SIOUT(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION RIG(9), DIVIS, D1IN, DIN, DPIN
      DOUBLE PRECISION EPTLIM, EPCLIM, GAM
C 
C -----------------------------------------------------------------------
      EPTLIM = COFIBR(1)
      EPCLIM = COFIBR(2)
      GAM    = COFIBR(3)
C 
C     Steph. 15/09/99. Je vire le cote non-lineaire sous compression
C     sens fibres qui fait chier un maximum.
C 
      D1IN   = ENDOIN(1)
      DIN    = ENDOIN(2)
      DPIN   = ENDOIN(3)
      CALL COPITD (9, ELAIN, ELAOUT)
C 
      IF (D2D3 .EQ. 2.) THEN
        IF (LOGIN(1)) THEN
	  ELAOUT(1) = ELAIN(1)*(1.D0-D1IN)
          ELAOUT(7) = ELAIN(7)*(1.D0-D1IN)
          ELAOUT(9) = ELAIN(9)*(1.D0-D1IN)
	ELSE
CD 	  ELAOUT(1) = ELAIN(1)*DMAX1(1+GAM*DMIN1(DEFIN(1), 0.D0), 0.D0)
	  ELAOUT(1) = ELAIN(1)
	END IF
        IF (LOGIN(2)) THEN
	  ELAOUT(2) = ELAIN(2)*(1.D0-DPIN)
          ELAOUT(8) = ELAIN(8)*(1.D0-DPIN)
	END IF
        ELAOUT(4) = ELAIN(4)*(1.D0-DIN)
C 
      ELSE IF (D2D3 .EQ. 3.) THEN
        IF (LOGIN(1)) THEN
          ELAOUT(1) = ELAIN(1)*(1.D0-D1IN)
          ELAOUT(7) = ELAIN(7)*(1.D0-D1IN)
          ELAOUT(9) = ELAIN(9)*(1.D0-D1IN)
        ELSE
CD 	  ELAOUT(1) = ELAIN(1)*DMAX1(1+GAM*DMIN1(DEFIN(1), 0.D0), 0.D0)
	  ELAOUT(1) = ELAIN(1)
	END IF
	IF (LOGIN(2)) THEN
	  ELAOUT(2) = ELAIN(2)*(1.D0-DPIN)
          ELAOUT(8) = ELAIN(8)*(1.D0-DPIN)
	END IF
CD 	IF (LOGIN(3)) THEN
CD        ELAOUT(3) = ELAIN(3)*(1.D0-DPIN)
CD 	END IF
	ELAOUT(4) = ELAIN(4)*(1.D0-DIN)
        ELAOUT(6) = ELAIN(6)*(1.D0-DIN)
CD      ELAOUT(5) = ELAIN(5)*(1.D0-DIN)
      END IF
C 
      DIVIS = 1.D0-(ELAIN(7)/ELAIN(1))*(ELAOUT(7)*ELAOUT(2))
     &            -(ELAIN(9)/ELAIN(1))*(ELAOUT(9)*ELAOUT(3))
     &            -(ELAIN(8)/ELAIN(2))*(ELAOUT(8)*ELAOUT(3))
     &            -(2*(ELAIN(7)/ELAIN(1))*ELAOUT(2)*ELAOUT(9)*
     &             (ELAIN(8)/ELAIN(2))*ELAOUT(3))
C 
      RIG(1) = ELAOUT(1)*(1.D0-(ELAIN(8)/ELAIN(2))*ELAOUT(8)*ELAOUT(3))
      RIG(1) = RIG(1)/DIVIS
      RIG(2) = ELAOUT(2)*(1.D0-(ELAIN(9)/ELAIN(1))*ELAOUT(9)*ELAOUT(3))
      RIG(2) = RIG(2)/DIVIS
      RIG(3) = ELAOUT(3)*(1.D0-(ELAIN(7)/ELAIN(1))*ELAOUT(7)*ELAOUT(2))
      RIG(3) = RIG(3)/DIVIS
      RIG(4) = ELAOUT(5)
      RIG(5) = ELAOUT(6)
      RIG(6) = ELAOUT(4)
      RIG(7) = ELAOUT(2)*((ELAIN(8)/ELAIN(2))*ELAOUT(3)+
     &                     ELAOUT(7)*(ELAIN(9)/ELAIN(1))*ELAOUT(3))
      RIG(7) = RIG(7)/DIVIS
      RIG(8) = ELAOUT(3)*(ELAOUT(9)+ELAOUT(7)*ELAOUT(8))
      RIG(8) = RIG(8)/DIVIS
      RIG(9) = ELAOUT(2)*(ELAOUT(7)+ELAOUT(9)*
     &                   (ELAIN(8)/ELAIN(2))*ELAOUT(3))
      RIG(9) = RIG(9)/DIVIS
C 
      SIOUT(1)= RIG(1)*DEFIN(1)+RIG(9)*DEFIN(2)+RIG(8)*DEFIN(6)
      SIOUT(2)= RIG(9)*DEFIN(1)+RIG(2)*DEFIN(2)+RIG(7)*DEFIN(6)
      SIOUT(3)= RIG(6)*DEFIN(3)
      SIOUT(4)= RIG(4)*DEFIN(4)
      SIOUT(5)= RIG(5)*DEFIN(5)
      SIOUT(6)= RIG(8)*DEFIN(1)+RIG(7)*DEFIN(2)+RIG(3)*DEFIN(6)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     CASIEX calcule le champ de contraintes associe au champ de deformations
C     envoye en entree. Ici, les modules elastiques sont imposes.
C     (EX pour explicite en deformations et modules)
C 
C     On envoie comme arguments :
C 
C     E ...... D2D3       option calcul contraintes planes ou contraintes HP
C     E ...... ELAIN(6)   coefficients elastiques vierges
C     E ...... COFIBR(3)  coefficients pour le comportement non lineaire
C     E ...... ENDOMS(3)  endommagement d1, d, d' et d''
C     E ...... LOGIN      logiques ranges E11, E22, E33 : l'endommagement est actif
C                         en traction, inactif en compression
C     E ...... DE         deformations elastiques dans la base d'orthotropie
C 
C     Et on recupere :
C 
C     S ...... C          contraintes chapeau dans la base dans la base d'orthotropie 
C 
      SUBROUTINE CASIEX (D2D3, ELAIN, COFIBR, ENDOMS, LOGIN, DE, C)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
C     S Souplesse initiale
C 
      DOUBLE PRECISION D2D3, ELAIN(9), COFIBR(3)
      DOUBLE PRECISION DE(6), ENDOMS(3), C(6)
      LOGICAL          LOGIN(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION ELAOUT(9)
      CHARACTER*6      IDPROG
      PARAMETER       (IDPROG='CASIEX')
C 
C -----------------------------------------------------------------------
      CALL ACTURI (D2D3, ELAIN, COFIBR, DE, ENDOMS, LOGIN,
     &             ELAOUT, C)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine recupere les coefficients des lois de comportement
C     de la couche.
C 
C     On envoie comme argument :
C 
C     E ...... D2D3       option de calcul : 2D-SIG ou 3D-SIG
C     E ...... N0         numero de couche
C 
C     Et on recupere :
C 
C     S ...... COENDO(12) coefficients l'endo (b,k,Y0,Yc,Ycs,n,b',k',Yo',Yc',Yts,n')
C     S ...... COPLAS(4)  coefficients pour la loi plastique (R0, beta, alpha, a2)
C     S ...... COFIBR(3)  coefficients pour le comp. NL fibres (EPCLIM, EPTLIM, GAMMA)
C     S ...... ELATRI(6)  coefficients elastiques initiaux
C     S ...... CRITER(5)  limites (22 trac., 22 comp., 12, 13, 23)
C 
      SUBROUTINE CONLIC (D2D3, N0,
C                        et on recupere...
     &                   COFIBR, COPLAS, COENDO, M0, ELATRI, CRITER)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER            N0
      DOUBLE PRECISION   D2D3
      DOUBLE PRECISION   COFIBR(3), COPLAS(4), COENDO(12)
      DOUBLE PRECISION   ELATRI(9), CRITER(5)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER M2, P, R, MNLIN, ADCAEN, R2, R1, M0, M1
C 
CD    LOGICAL LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CONLIC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('SOUP-ORTHO', R)
      CALL ADTBDM ('ELAS-ORTHO', R2)
      CALL ADTBDM ('HOOK-ORTHO', R1)
      CALL ADTBM  ('TYP-COUCHE', M2)
      CALL ADTBDM ('CARAC-NONL', ADCAEN)
C 
C     P caracterise le type de comportement
C 
      P = M(M2+N0-1)
      M0 = R+17*(P-1)
      M1 = R2+9*(P-1)
C 
C     Remplissage de ELATRI
C 
      ELATRI(1)  = DM(M1)
      ELATRI(2)  = DM(M1+1)
      ELATRI(3)  = DM(M1+2)
      ELATRI(4)  = DM(M1+3)
      ELATRI(5)  = DM(M1+4)
      ELATRI(6)  = DM(M1+5)
      ELATRI(7)  = DM(M1+6)
      ELATRI(8)  = DM(M1+7)
      ELATRI(9)  = DM(M1+8)
C 
      MNLIN = ADCAEN+25*(P-1)
C 
      COFIBR(1)  = DM(MNLIN)
      COFIBR(2)  = DM(MNLIN + 1)
      COFIBR(3)  = DM(MNLIN + 2)
      COPLAS(1)  = DM(MNLIN + 3)
      COPLAS(2)  = DM(MNLIN + 4)
      COPLAS(3)  = DM(MNLIN + 5)
      COPLAS(4)  = DM(MNLIN + 6)
      COENDO(1)  = DM(MNLIN + 7)
      COENDO(2)  = DM(MNLIN + 8)
      COENDO(3)  = DM(MNLIN + 9)
      COENDO(4)  = DM(MNLIN + 10)
      COENDO(5)  = DM(MNLIN + 11)
      COENDO(6)  = DM(MNLIN + 12)
      COENDO(7)  = DM(MNLIN + 13)
      COENDO(8)  = DM(MNLIN + 14)
      COENDO(9)  = DM(MNLIN + 15)
      COENDO(10) = DM(MNLIN + 16)
      COENDO(11) = DM(MNLIN + 17)
      COENDO(12) = DM(MNLIN + 18)
C 
      CRITER(1)  = DM(MNLIN + 19)
      CRITER(2)  = DM(MNLIN + 20)
      CRITER(3)  = DM(MNLIN + 21)
      CRITER(4)  = DM(MNLIN + 22)
      CRITER(5)  = DM(MNLIN + 23) 
      D2D3       = DM(MNLIN + 24)
C 
C     Plus simple pour ecrire l'endommagement : YC  = DSQRT(YC)
C                                               YCP = DSQRT(YCP)
C                                               Y0  = DSQRT(Y0)/YC
      COENDO(4)  = DSQRT(COENDO(4))
      COENDO(10) = DSQRT(COENDO(10))
      COENDO(3)  = DSQRT(COENDO(3))/COENDO(4)
      COENDO(9)  = DSQRT(COENDO(9))/COENDO(10)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
      SUBROUTINE SEUIL (R0, AL, BE, P, R)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  AL, BE
      DOUBLE PRECISION  P, R, R0
C 
C -----------------------------------------------------------------------
      R = R0 + BE*(P**AL)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
      SUBROUTINE TESCAS (C, LOG)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION C(6)
      LOGICAL          LOG(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """""""""""""""""""""""""""""""""
C 
      CHARACTER*6      IDPROG
      PARAMETER       (IDPROG='TESCAS')
C 
C -----------------------------------------------------------------------
      LOG(1) = .FALSE.
      LOG(2) = .FALSE.
      LOG(3) = .FALSE.
C 
      IF (C(1) .GE. 0.D0) LOG(1) = .TRUE.
      IF (C(2) .GE. 0.D0) LOG(2) = .TRUE.
      IF (C(6) .GE. 0.D0) LOG(3) = .TRUE.
C 
      RETURN
      END
