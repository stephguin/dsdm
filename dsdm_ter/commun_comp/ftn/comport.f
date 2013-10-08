C     Routine d'integration du comportement couche appelee par etaloc.
C     Endommagement et plasticite. Ecriture en contrainte.
C 
C     E ...... ELAIN(6)   coefficients elastiques vierges (E1,E2,E3,G12,G23,G13)
C     E ...... COENDO(12) coefficients pour l'endo (b,k,Y0,Yc,Ycs,n,b',k',Yo',Yc',Yts,n')
C     E ...... COPLAS(4)  coefficients pour la plasticite (R0, beta, alpha, a2)
C     E ...... COFIBR(3)  coefficients pour le comp. NL fibres (EPCLIM, EPTLIM, GAMMA)
C     E ...... CRITER(5)  limites (22 trac., 22 comp., 12, 13, 23)
C 
      SUBROUTINE INCOCO (INTERV, MULPPG,
C                        Parametres materiau
     &		         D2D3, COFIBR, COPLAS, COENDO, SORT, ELAIN,
C                        Taux admissibles dans la base ortho. en Entree     
     &                   EPPAO, SIPAO,
C                        Valeurs admissibles dans la base ortho. (en Entree-Sortie)     
     &                   EPSAOR, SIGAOR,
C                        Quantites chapeau dans la base ortho. (en Entree-Sortie)     
     &                   EPSCOR, SIGCOR, RPRE, PPRE, EPSEOR,
C                        Quantites chapeau dans la base ortho. stockees (en Entree) 
     &                   ENDOME, MODULE, CASIN, YDPIN, CRICOE, EPLAOE,
C                        Quantites chapeau dans la base ortho. stockees (en Sortie) 
     &                   ENDOMS, MODULS, CAOUT, YDPOUT, CRICOS, EPLAOS,
C                        Quantites erreur dans la base ortho. (en Sortie)
     &                   ERRVIS, VEPPOR, VSIPOR, DSIGOR, NUMERR, DENOMI,
C                        Logiques a passer a NETLOC
     &                   RUPPRE, LOGIMP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'typcal.h'
C         
      DOUBLE PRECISION   INTERV, MULPPG(YINTEG), ELAIN(9), D2D3
      DOUBLE PRECISION   COFIBR(3), COPLAS(4), COENDO(12)
      DOUBLE PRECISION   EPPAO(6*YINTEG), SIPAO(6*YINTEG)
      DOUBLE PRECISION   EPSAOR(6*YINTEG), SIGAOR(6*YINTEG)
      DOUBLE PRECISION   EPSCOR(6*YINTEG)
      DOUBLE PRECISION   SIGCOR(6*YINTEG), EPSEOR(6*YINTEG)
      DOUBLE PRECISION   RPRE(YINTEG), PPRE(YINTEG)
      DOUBLE PRECISION   ENDOME(3), MODULE(3*YINTEG), CRICOE
      DOUBLE PRECISION   ENDOMS(3), MODULS(3*YINTEG), CRICOS
      DOUBLE PRECISION   EPLAOE(6*YINTEG), YDPIN, EPLAOS(6*YINTEG)
      DOUBLE PRECISION   ERRVIS(YINTEG), VEPPOR(6*YINTEG), YDPOUT
      DOUBLE PRECISION   VSIPOR(6*YINTEG), DSIGOR(6*YINTEG) 
      DOUBLE PRECISION   NUMERR, DENOMI
      LOGICAL            LOGIMP, RUPPRE
      LOGICAL            CASIN(3*YINTEG), CAOUT(3*YINTEG)
      INTEGER            SORT
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER            Y, DEB, DEC, DBSIG, NBTOUR, CPT3, I
      INTEGER            CPT2(YINTEG), CPT1(YINTEG), CPT
      INTEGER            CP1(YINTEG), CP2(YINTEG)
      DOUBLE PRECISION   SIGADD(6), DENLOC, DEPLAO(6*YINTEG)
      DOUBLE PRECISION   TEST, ENDOMP(3), ENDODO(3)
      DOUBLE PRECISION   RPREIN(3), PPREIN(3), EPSELO(6), NORRES(YINTEG)
      DOUBLE PRECISION   NORRIN(YINTEG), RPREDE(3), PPREDE(3)
      DOUBLE PRECISION   EPSEIN(6*YINTEG), EPSEOUT(6*YINTEG)
      DOUBLE PRECISION   EPSEOUTP(6*YINTEG)
      INTEGER            AM2LC, ADM2LC, SIGCHA, LONPRO
      LOGICAL            LOGPLA, RUPCOU(2), LOGZOB
      LOGICAL            CASDEB(3), CASBID(3)
      DOUBLE PRECISION   DEFBID(6), UNITE(4), SIGBID(6), MODBID(9)
C 
      CHARACTER*6  IDPROG
      PARAMETER   (IDPROG='INCOCO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C 
C     On calcule epsilon_adm et sigma_adm
C 
      LONPRO = 6*YINTEG
      CALL ADDMAD (LONPRO, EPPAO, EPSAOR, EPSAOR)
      CALL ADDMAD (LONPRO, SIPAO, SIGAOR, SIGAOR)
      CALL COPITD (LONPRO, EPLAOE, EPLAOS)
      CALL COPITD (LONPRO, EPSAOR, EPSCOR)
      CALL COPITD (3, ENDOME, ENDOMS)
C 
C     Initialisation elastique : on considere que tout l'increment
C     de deformation admissible est elastique.
C 
      CALL COPITD (LONPRO, EPSEOR, EPSEIN)
      CALL ADDMAD (LONPRO, EPPAO, EPSEIN, EPSEIN)
      CALL COPITD (LONPRO, EPSEIN, EPSEOUT)
      CALL COPITD (LONPRO, EPSEIN, EPSEOUTP)
      CALL POUSMD (LONPRO, SIGCHA)
C 
      YDPOUT = YDPIN
      DEB = 1
      DO Y = 1, YINTEG
        RPREDE(Y) = RPRE(Y)
	PPREDE(Y) = PPRE(Y)
	CALL SCAVEC (6, EPSEIN(DEB), EPSEIN(DEB), NORRIN(Y))
	DEB = DEB + 6
      END DO
C 
C     Initialisations pour la boucle i 
C 
      CPT = 0
      NBTOUR = 0
      TEST   = 1.D0
      DO I = 1, 2
        RUPCOU(I) = .FALSE.
      END DO
C 
C     Boucle i : en precision sur la repartition endo/plasto
C 
      DO WHILE (TEST .GT. 1.D-3)
C       
	NBTOUR = NBTOUR+1
CD 	IF (LOGIMP) THEN
CD 	  CALL IMPET ('BOUCLE I ', NBTOUR)
CD 	  CALL IMPTDT ('VERIF ENDO   ', ENDOMS, 1, 3)
CD 	  CALL IMPTDT ('VERIF PLASTO ', DEPLAO, 1, 6)
CD 	END IF
	CALL COPITD (LONPRO, EPSEOUT, EPSEOUTP)
        CALL COPITD (3, ENDOMS, ENDOMP)
C 
C       Integration de l'endommagement (constant en Y)
C 
	CALL REPELA (D2D3, ELAIN, INTERV, COFIBR, COENDO,
     &               EPSEOUTP, ENDOME, ENDOMS, YDPIN, YDPOUT,
     &               CASIN, CPT3, RUPPRE, RUPCOU, LOGIMP)
C 
C       Initialisations pour la boucle en Y
C 
	DEB = 1
	DEC = 1
C 
C       Integration de la plasticite (variable en Y)
C 
	LOGPLA = .FALSE.
	DO Y  = 1, YINTEG
	  IF (LOGIMP .AND. (Y .EQ. 3)) LOGPLA = .TRUE.
	  RPREIN(Y) = RPREDE(Y)
	  PPREIN(Y) = PPREDE(Y)
	  CALL NINPLC (D2D3, ELAIN, COFIBR, COPLAS, ENDOMS,
     &                 CASIN(DEC), RPREIN(Y), PPREIN(Y), DEPLAO(DEB),
     &                 EPSEIN(DEB), EPSEOUT(DEB), CP1(Y), CP2(Y),
     &                 LOGPLA, RUPPRE, RUPCOU, LOGIMP)
	  CPT     = CPT+CP1(Y)
	  CPT     = CPT+CP2(Y)
	  DEC     = DEC+3
	  DEB     = DEB+6
	  LOGPLA = .FALSE.
	END DO
C 
        CPT  = CPT+CPT3
	DEB  = 1
	CALL SOUMAD (3, ENDOMS, ENDOMP, ENDODO)
	CALL SCAVEC (3, ENDODO, ENDODO, TEST)
   	TEST = DSQRT(TEST)
CD 	DO Y = 1, YINTEG
CD 	  CALL SOUMAD (6, EPSEOUTP(DEB), EPSEOUT(DEB), EPSELO)
CD  	  CALL SCAVEC (6, EPSELO, EPSELO, NORRES(Y))	
CD 	  NORRES(Y) = NORRES(Y)/NORRIN(Y)
CD        TEST = DMAX1 (TEST, DSQRT(NORRES(Y)))
CD        DEB  = DEB+6
CD      END DO
CD 	CALL IMPTDT ('VERIF ENDO ', ENDOMS, 1, 3)
CD 	CALL IMPTDT ('VERIF PLAS ', RPREIN, 1, 3)
CD 	CALL IMPDT ('VERIF TEST I ', TEST)
C 
C     Fin de la boucle i
C 
      END DO
C 
CD    CALL MESSAO (' ')
      DEB = 1
      DO Y = 1, YINTEG
	RPRE(Y) = RPREIN(Y)
	PPRE(Y) = PPREIN(Y)
	CALL ADDMAD (6, EPLAOS(DEB), DEPLAO(DEB), EPLAOS(DEB))
	DEB = DEB + 6
	LOGPLA = .FALSE.
      END DO
C 
      IF (RUPCOU(1) .AND. RUPCOU(2)) RUPPRE = .TRUE.
      CALL COPITD (LONPRO, EPSEOUT, EPSEOR)
C 
      DEB   = 1
      DEC   = 1
      DBSIG = SIGCHA
C 
      DO Y = 1, YINTEG
	CALL CASITR (D2D3, ELAIN, COFIBR, ENDOMS,
     &               CASIN(DEC), EPSEOUT(DEB),
     &               DM(DBSIG), MODULS(DEC), CAOUT(DEC),.FALSE.)
	DEC   = DEC   + 3
	DEB   = DEB   + 6
	DBSIG = DBSIG + 6
      END DO
C 
      DEB  = 1
      NUMERR = 0.D0
      DENOMI = 0.D0
      DO Y = 1, YINTEG
C 
C       Calcul de sigma_point_chapeau X DELTA(T) = sigma_chapeau(t+1) -
C       sigma_chapeau(t) stocke dans vsipor
C 
        CALL SOUMAD (6, DM(SIGCHA), SIGCOR(DEB), VSIPOR(DEB))
        CALL COPITD (6, DM(SIGCHA), SIGCOR(DEB))
C 
C       Calcul de sigma_point_admissible - sigma_point_chapeau stocke dans VSIPOR
C 
        CALL SOUMAD (6, SIPAO(DEB), VSIPOR(DEB), VSIPOR(DEB))
C 
C       Calcul de sigma_chapeau - sigma_admissible stocke dans DSIGOR
C       et contribution du point de Gauss courant au numerateur de l'erreur
C       exploitee en visu pour les carto d'ecart chapeau/admissible
C 
        CALL SOUMAD (6, SIGCOR(DEB), SIGAOR(DEB), DSIGOR(DEB))
        CALL TRAC12 (SORT, DSIGOR(DEB), DSIGOR(DEB), ERRVIS(Y))
        ERRVIS(Y) = ERRVIS(Y)*MULPPG(Y)
C 
C       Calcul de sigma_chapeau + sigma_admissible stocke dans SIGADD
C       et contribution du point de Gauss courant au denominateur de l'erreur
C       On calcule NUMERR = int_espace((sig_cha-sig_adm)*(1/K0)*(sig_cha-sig_adm))
C               et DENOMI = int_espace((sig_cha-sig_adm)*(1/K0)*(sig_cha-sig_adm))
C 
        CALL ADDMAD (6, SIGCOR(DEB), SIGAOR(DEB), SIGADD)
        CALL TRAC12 (SORT, SIGADD, SIGADD, DENLOC) 
        DENOMI = DENOMI + MULPPG(Y)*DENLOC
        NUMERR = NUMERR + ERRVIS(Y)    
        DEB    = DEB + 6
        SIGCHA = SIGCHA + 6
C 
      END DO
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
C     Cette routine calcule l'endommagement pour une deformation elastique
C     imposee. Boucle J en precision sur l'endommagement, sous-iteration de
C     la boucle I de partition elasticite/plasticite dans INCOCO.
C 
C     On envoie comme arguments :
C 
C     E ...... D2D3       option de calcul : 2D-SIG ou 3D-SIG
C     E ...... ELAIN(9)  coefficients elatiques initiaux
C     E ...... INTERV     longueur de l'intervalle de temps courant
C     S ...... COENDO(12) coefficients l'endo (b,k,Y0,Yc,Ycs,n,b',k',Yo',Yc',Yts,n')
C     S ...... COPLAS(4)  coefficients pour la loi plastique (R0, beta, alpha, a2)
C     S ...... COFIBR(3)  coefficients pour le comp. NL fibres (EPCLIM, EPTLIM, GAMMA)
C     S ...... ELATRI(6)  coefficients elastiques initiaux
C     S ...... CRITER(5)  limites (22 trac., 22 comp., 12, 13, 23)
C     E ...... EPSELA(6)  deformations elastiques moyennees dans l'epaisseur du pli
C     E ...... DIN(3)     endommagements d1, d, d' et d'' au pas de temps precedent
C     E ...... CASPRE(3)  tableau de 3 logiques caracterisant l'etat de contraintes
C                         au pas de temps precedent : traction (TRUE) ou compression
C                         (FALSE) pour SIG11, SIG22 et SIG33
C 
C     Et en entree-sortie :
C 
C     ES ..... RUPPRE     tableau de logiques pour activation des mecanismes rapides
C     ES ..... RUPCOU(2)  couple de logiques pour la rupture en traction au pas courant
C 
C     Et on recupere :
C 
C     S ...... DOUT(3)    endommagements d1, d, d' et d'' au pas de temps courant
C 
      SUBROUTINE REPELA (D2D3, ELAIN, INTERV, COFIBR, COENDO, 
     &                   EPSELA, DIN, DOUT, YDPIN, YDPOUT,
     &                   CASPRE, CPT3, RUPPRE, RUPCOU, LOGIMP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'typcal.h'
C 
      DOUBLE PRECISION   ELAIN(9), INTERV
      DOUBLE PRECISION   COFIBR(3), COENDO(12)
      DOUBLE PRECISION   EPSELA(6*YINTEG)
      DOUBLE PRECISION   D2D3
      DOUBLE PRECISION   DIN(3), DOUT(3), YDPIN
      LOGICAL            LOGIMP, RUPPRE, RUPCOU(2)
      LOGICAL            CASPRE(3*YINTEG)
      INTEGER            CPT3
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION   EPS, C(6), DEF(6), EPSIN(6), YDPOUT
      DOUBLE PRECISION   EPTLIM, EPCLIM, YTS
      DOUBLE PRECISION   YD1, YD, MODOUT(3), DIVISY, SIGOUT(6)
      DOUBLE PRECISION   DOUTP(3)
      INTEGER            Y, I, DEB, DEC
      LOGICAL            CASCOU(3*YINTEG), NOCVD(3)
      INTEGER            AM2LC, ADM2LC
      CHARACTER*6        IDPROG
      PARAMETER         (IDPROG='REPELA')
      DOUBLE PRECISION   PRECISION
      PARAMETER         (PRECISION = 1.D-6)
      INTEGER            NTRMAX
      PARAMETER         (NTRMAX = 300)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     Initialisations au pas de temps precedent
C 
      DO I = 1, 3
        DOUT(I)  = DIN(I)
        NOCVD(I) = .FALSE.
      END DO
C 
      IF (.NOT. LENDCO) GOTO 2
C 
      DIVISY = DBLE(YINTEG)
      EPTLIM = COFIBR(1)
      EPCLIM = COFIBR(2)
      YTS    = COENDO(11)
C 
      IF (DIN(1) .EQ. 1.D0) NOCVD(1) = .TRUE.
      IF (DIN(2) .EQ. 1.D0) NOCVD(2) = .TRUE.
      IF (DIN(3) .EQ. 1.D0) NOCVD(3) = .TRUE.
C 
C     Initialisations de la boucle en Y.
C 
      DEB = 1 
      DO I = 1, 6
	DEF(I) = 0.D0
      END DO
C 
C     BOUCLE EN Y : Calcul de la moyenne des deformations dans l'epaisseur du pli
C 
      DO Y = 1, YINTEG
	CALL COPITD (6, EPSELA(DEB), EPSIN)
	DEF(1) = DEF(1)+ EPSIN(1)/DIVISY
	DEF(2) = DEF(2)+ DMAX1(EPSIN(2), 0.D0)/DIVISY
	DEF(3) = DEF(3)+ EPSIN(3)/DIVISY
	DEF(4) = DEF(4)+ EPSIN(4)/DIVISY
	DEF(5) = DEF(5)+ EPSIN(5)/DIVISY
	DEF(6) = DEF(6)+ DMAX1(EPSIN(6), 0.D0)/DIVISY
	DEB  = DEB + 6
      END DO
C 
C     FIN DE LA BOUCLE EN Y
C 
C     Calcul de Yd1 et Yd explicite en moyennes des deformations
C 
      YD1 = 0.D0
      IF (DEF(1) .GT. EPTLIM) THEN
        YD1 = 10.D0*(DEF(1)-EPTLIM)/EPTLIM
      ELSE IF (DEF(1) .LT. EPCLIM) THEN
        YD1 = 10.D0*DABS((DEF(1)-EPCLIM))/EPCLIM
      END IF
      YD    = ELAIN(4)*(DABS(DEF(3)))**2.D0
     &      + ELAIN(6)*(DABS(DEF(5)))**2.D0
CD   &      + ELAIN(5)*(DABS(DEF(4)))**2.D0
C 
C     Initialisations pour la boucle en precision
C 
      CPT3 = 0
      EPS  = 1.D0
C 
C     BOUCLE EN PRECISION J : calcul de Yd', d, d' (implicites en contraintes)
C 
      DO WHILE ((EPS .GT. PRECISION) .AND. (CPT3 .NE. NTRMAX))
C 
C       Initialisations sur l'iteration precedente
C 
	DOUTP(1) = DOUT(1)
	DOUTP(2) = DOUT(2)
	DOUTP(3) = DOUT(3)
	CPT3     = CPT3+1
CD 	IF (LOGIMP) THEN
CD 	  CALL IMPET ('BOUCLE J ', CPT3)
CD 	  CALL IMPTDT ('VERIF ENDO ', DOUT, 1, 3)
CD 	  CALL IMPTLT ('VERIF RUPT ', RUPCOU, 1, 2)
CD 	  CALL IMPTLT ('VERIF TRAC ', CASCOU, 1, 9)
CD 	  CALL IMPDT ('VERIF YD  ', YD)
CD 	  CALL IMPTDT ('VERIF DEF ', DEF, 1, 6)
CD 	  CALL IMPTDT ('VERIF ELAIN ', ELAIN, 1, 9)
CD 	END IF
C 
C       Initialisations pour la boucle en Y
C 
	DEB  = 1
	DEC  = 1
	DO I = 1, 6
	  C(I) = 0.D0
	END DO
C  
C       BOUCLE EN Y : Calcul de la moyenne des contraintes dans l'epaisseur du pli
C 
	DO Y = 1, YINTEG
	  CALL COPITD (6, EPSELA(DEB), EPSIN)
	  CALL CASITR (D2D3, ELAIN, COFIBR, DOUTP,
     &                 CASPRE(DEC), EPSIN,
     &                 SIGOUT, MODOUT, CASCOU(DEC), .FALSE.)
	  C(1) = C(1)+(SIGOUT(1)*SIGOUT(1)/DIVISY)
          C(2) = C(2)+(DMAX1(SIGOUT(2), 0.D0)*
     &                 DMAX1(SIGOUT(2), 0.D0)/DIVISY)
          C(3) = C(3)+(SIGOUT(3)*SIGOUT(3)/DIVISY)
	  C(4) = C(4)+(SIGOUT(4)*SIGOUT(4)/DIVISY)
	  C(5) = C(5)+(SIGOUT(5)*SIGOUT(5)/DIVISY)
	  C(6) = C(6)+(DMAX1(SIGOUT(6), 0.D0)*
     &                 DMAX1(SIGOUT(6), 0.D0)/ DIVISY)
	  DEB  = DEB + 6
	  DEC  = DEC + 3
        END DO
C 	
C       FIN DE LA BOUCLE EN Y
C 
	IF ((DOUTP(3) .NE. 1.D0) .AND.
     &      (.NOT. (RUPCOU(1) .AND. (.NOT. RUPCOU(2)))) .AND.
     &      (.NOT. (RUPCOU(1) .AND. RUPCOU(2)))) THEN
          YDPOUT=  C(2)/(2.D0*ELAIN(2)*(1.D0-DOUTP(3))*(1.D0-DOUTP(3)))
CD   &           + C(6)/(2.D0*ELAIN(3)*(1.D0-DOUTP(3))*(1.D0-DOUTP(3)))
	END IF
	CALL VAENDO (INTERV, COENDO, DIN, YD1, YD, YDPOUT,
     &               NOCVD, DOUT, RUPCOU)
C 
C       Calcul de la precision atteinte
C 
	EPS = DABS  (DOUT(1)-DOUTP(1))
	EPS = DMAX1 (EPS, DABS (DOUT(2)-DOUTP(2)))
	EPS = DMAX1 (EPS, DABS (DOUT(3)-DOUTP(3)))
C 
      END DO
C 
C     FIN DE LA BOUCLE EN PRECISION
C 
      IF (CPT3 .EQ. NTRMAX) THEN
	CALL MESSAO ('PROBLEME DANS LA BOUCLE ENDO '//IDPROG)
        CALL IMPDT ('PRECISION ATTEINTE ', EPS)
	STOP
      END IF
C 
2     CONTINUE
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine calcule la correction plastique pour une deformation totale
C     imposee. Boucle K en precision sur la plasticite, sous-iteration de la
C     boucle I de partition elasticite/plasticite dans INCOCO. Pour eviter les 
C     problemes, on bloque la plasticite dans les cas :
C 
C            -----> rupture fibre au pas precedent
C            -----> rupture en traction au pas precedent
C            -----> rupture en traction confirmee au pas courant (RUPCOU = (TRUE, TRUE))
C            -----> endom. en cisaillement au pas precedent > 0.99
C     
C 
C     On envoie comme arguments :
C 
C     E ...... D2D3             option de calcul : 2D-SIG ou 3D-SIG
C     E ...... ELAIN            coefficients elatiques vierges
C     E ...... GAM, R0, BE  |   caracterisiques non-lineaires
C     E ...... EPTLIM       |   pour la plasticite
C     E ...... ENDOIN(3)        endommagements d1, d, d' et d'' au pas de temps precedent
C     E ...... CASPRE
C     E ...... EPSEOR           deformations elastiques au pas precedent
C     E ...... EPSEIN           deformations elastiques a l'iteration I precedente
C     E ...... CASPRE           tableau de 3 logiques caracterisant l'etat de contraintes
C                               au pas de temps precedent : traction (TRUE) ou compression
C                               (FALSE) pour SIG11, SIG22 et SIG33
C 
C     Et en Entree-Sortie :
C                             
C     ES ..... RPRE             seuil plastique
C     ES ..... PPRE             def. plastique cumulee
C     ES ..... EPSPO            increment de deformations plastiques
C     ES ..... RUPPRE       |   tableau de logiques pour activation des mecanismes rapides
C     ES ..... RUPCOU(2)    |   couple de logiques pour la rupture en traction au pas courant
C 
C 
C     Et on recupere :
C                             
C     S ...... EPSEOUT          deformations elastiques a l'iteration I+1
C 
      SUBROUTINE NINPLC (D2D3, ELAIN, COFIBR, COPLAS, ENDOIN,
     &                   CASPRE, RPRE, PPRE, EPSPO, EPSEIN, EPSEOUT,
     &                   CPTNEG, CPTBLB, LOGIMP, RUPPRE, RUPCOU,
     &                   LOGZOB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'typcal.h'
C 
      DOUBLE PRECISION D2D3
      DOUBLE PRECISION ELAIN(9), COFIBR(3), COPLAS(4)
      DOUBLE PRECISION ENDOIN(3)
      LOGICAL          CASPRE(3)
      DOUBLE PRECISION RPRE, PPRE, EPSPO(6)
      DOUBLE PRECISION EPSEIN(6), EPSEOUT(6)
      INTEGER          CPTNEG, CPTBLB
      LOGICAL          LOGIMP, RUPPRE, RUPCOU(2), LOGZOB
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      DOUBLE PRECISION R0, BE, AL, A2, ACP1, ACP2
      DOUBLE PRECISION NORSIE, KEDHSI(6), RLO, MUL, ACCP, DELP
      DOUBLE PRECISION EPSEOUTP(6), DIVIS, VALF, VALFP, RP, D(3)
      DOUBLE PRECISION HSIGMA(6), HSIGTI(6), DEPSPO(6), SIGMOR(6)
      DOUBLE PRECISION DEPSEO(6), EAC(3), SIGMAI(6), SIGMAP(6)
      DOUBLE PRECISION SITILD(6), NORSII, NORSIP, PPREI, EPSELO(6)
      DOUBLE PRECISION DELP1, DELP2
      DOUBLE PRECISION PRECISION, ACCP1, ZBLB
      PARAMETER        (PRECISION = 1.D-6)
      LOGICAL          CASCOU(3), CPTCP1, CPTCP2
      INTEGER          I, COMPT, KK
      CHARACTER*6      IDPROG
      PARAMETER       (IDPROG='NINPLC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      D(1)  = ENDOIN(1)
      D(2)  = ENDOIN(2)
      D(3)  = ENDOIN(3)
      R0    = COPLAS(1)
      BE    = COPLAS(2)
      AL    = COPLAS(3)
      A2    = COPLAS(4)
      ACCP = 0.D0
      DO I = 1, 6
        EPSPO(I)  = 0.D0
	DEPSPO(I) = 0.D0
      END DO
C 
C     Initialisation des compteurs
C 
      CPTNEG = 0
      CPTBLB = 0
C 
      CALL COPITD (6, EPSEIN, EPSEOUT)
      IF (.NOT. LPLACO) GOTO 5
C 
      CALL CASITR (D2D3, ELAIN, COFIBR, D, CASPRE, EPSEIN,
     &             SIGMOR, EAC, CASCOU, .FALSE.)
CD    IF (LOGIMP) CALL IMPTDT ('VERIF SIGMOR ', SIGMOR, 1, 6)
C 
C     Blocage de la plasticite : endo fragile detecte au pas precedent
C 
CD    IF ((RUPPRE .AND. (CASCOU(1) .OR. CASCOU(2) .OR. CASCOU(3))) .OR.
CD   &    (D(2) .EQ. 1.D0)) GOTO 5
      IF ((RUPPRE) .OR. (D(2) .EQ. 1.D0)) GOTO 5
C 
C     Blocage de la plasticite : endo fragile detecte a l'iteration i courante
C 
      IF (RUPCOU(1) .OR. RUPCOU(2)) GOTO 5
C 
C     Calcul de la fonction seuil pour un increment de deformation elastique
C     initialise sur l'increment de deformation total admissible.
C 
      CALL SIGTIL (D2D3, CASCOU, D, SIGMOR, SIGMOR)
      CALL NSIGEQ (D2D3, A2, SIGMOR, NORSIE)
      VALF = NORSIE-RPRE
      IF (VALF .LT. 0.D0) THEN
        GOTO 5
      END IF
C 
CD    CALL MESSAO ('YA PLASTO')
      PPREI  = PPRE
      NORSII = NORSIE
      CPTCP1 = .FALSE.
      CPTCP2 = .FALSE.
C 
C -----------------------------------------------------------------------
C     Cas p=0 au pas de temps precedent : p**(1.D0-AL); correction grossiere
C -----------------------------------------------------------------------
C 
      KK = 0
      IF (PPRE .EQ. 0.D0) THEN
C 
	ACCP = PRECISION
1       CONTINUE
   	KK = KK + 1
CD 	IF (LOGIMP) CALL IMPET ('BOUCLE KK  ', KK)
CD 	IF (LOGIMP) CALL IMPDT ('VERIF VALF ', VALF)
CD 	VALFP = VALF
        CALL SEUIL (R0, AL, BE, ACCP, RLO)
        CALL NSIGEQ (D2D3, A2, SIGMOR, NORSII)
        CALL CAHSIP (D2D3, A2, SIGMOR, HSIGTI)
        MUL = ACCP/NORSII
        CALL MUMARE (MUL, 6, HSIGTI, DEPSPO)
        CALL MUMARE (-1.D0, 6, DEPSPO, DEPSEO)
   	CALL ADDMAD (6, EPSEIN, DEPSEO, EPSEOUTP)
   	CALL CASIEX (D2D3, ELAIN, COFIBR, D, CASCOU, EPSEOUTP, SIGMOR)
        CALL SIGTIL (D2D3, CASCOU, D, SIGMOR, SIGMOR)
        VALF = NORSII-RLO
CD 	IF (VALF .GT. VALFP) THEN
CD 	  CALL MESSAO ('OH OH')
CD 	END IF
	CALL COPITD (6, EPSEOUTP, EPSEOUT)
	IF ((VALF .LT. 0.D0) .AND. (.NOT. CPTCP1)) THEN
	  PPRE = 0.D0
	  GOTO 5
	END IF
	IF ((VALF .LT. 0.D0) .AND. CPTCP1) THEN
	  ACCP = ACCP1
	  GOTO 4
	END IF
	IF (VALF .GT. 0.D0) THEN
	  CPTCP1 = .TRUE.
	  ACCP1 = ACCP
	  ACCP = ACCP+PRECISION
	  GOTO 1
	END IF
C 
      END IF
C 
C -----------------------------------------------------------------------
C     Cas p<>0 au pas de temps precedent, p initialise dans ce qui precede
C -----------------------------------------------------------------------
C 
2     CONTINUE
      PPREI = PPRE + ACCP
      CALL SEUILP (AL, BE, PPREI, RP)
      CALL SEUIL (R0, AL, BE, PPREI, RLO)
      CALL CASIEX (D2D3, ELAIN, COFIBR, D, CASCOU, EPSEOUT, SIGMAI)
      CALL SIGTIL (D2D3, CASCOU, D, SIGMAI, SIGMAP)
      CALL NSIGEQ (D2D3, A2, SIGMAP, NORSII)
      VALF = NORSII-RLO
      CALL CAHSIP (D2D3, A2, SIGMAP, HSIGTI)
      CALL CASIEX (D2D3, ELAIN, COFIBR, D, CASCOU, HSIGTI, KEDHSI)
      CALL SCAVEC (6, HSIGTI, KEDHSI, DIVIS)
      DIVIS = RP+(DIVIS/(NORSII*NORSII))
C 
      COMPT = 0
      ZBLB  = 1.D0
      DO WHILE (VALF .GT. 0.D0)
        COMPT = COMPT + 1
CD 	IF (LOGIMP) CALL IMPDT ('VERIF ZBLB ', ZBLB)
	DELP  = VALF/DIVIS
	DELP2 = 0.D0
   	ACCP  = DELP
	PPREI = PPRE + ACCP
3       CONTINUE
	IF (ZBLB .LT. PRECISION) THEN
	  DELP = DELP2
	  GOTO 4
	END IF
	CPTBLB = CPTBLB+1
	CALL SEUIL (R0, AL, BE, PPREI, RLO)
	MUL  = DELP/NORSII
        CALL MUMARE (MUL, 6, HSIGTI, DEPSPO)
        CALL MUMARE (-1.D0, 6, DEPSPO, DEPSEO)
	CALL ADDMAD (6, EPSEOUT, DEPSEO, EPSEOUTP)
	CALL SOUMAD (6, EPSEOUT, EPSEOUTP, EPSELO)
	CALL SCAVEC (6, EPSELO, EPSELO, ZBLB)
	ZBLB = DSQRT(ZBLB)
   	CALL CASIEX (D2D3, ELAIN, COFIBR, D, CASCOU, EPSEOUTP, SIGMAI)
        CALL SIGTIL (D2D3, CASCOU, D, SIGMAI, SIGMAP)
        CALL NSIGEQ (D2D3, A2, SIGMAP, NORSIP)
        VALF = NORSIP-RLO
	IF ((VALF .LT. 0.D0) .AND. (.NOT. CPTCP2)) THEN
	  CPTCP2 = .TRUE.
	  DELP1 = DELP
	  ACCP = DELP
	  DELP = DELP * .5D0
	  ACCP = DELP
	  PPREI = PPRE + ACCP
	  GOTO 3
	END IF
	IF ((VALF .LT. 0.D0) .AND. CPTCP2) THEN
	  DELP2 = DELP
	  ACCP = DELP
	  DELP = DELP * .5D0
	  ACCP = DELP
	  PPREI = PPRE + ACCP
	  GOTO 3
	END IF
	CALL SEUILP (AL, BE, PPREI, RP)
	CALL CAHSIP (D2D3, A2, SIGMAP, HSIGTI)
        CALL CASIEX (D2D3, ELAIN, COFIBR, D, CASCOU, HSIGTI, KEDHSI)
        CALL SCAVEC (6, HSIGTI, KEDHSI, DIVIS)
	DIVIS = RP+(DIVIS/(NORSIP*NORSIP))
	CALL COPITD (6, EPSEOUTP, EPSEOUT)
      END DO
C 
4     CONTINUE
C 
      PPRE = PPRE + ACCP
      CALL SEUILP (AL, BE, PPRE, RP)
      CALL SEUIL (R0, AL, BE, PPRE, RPRE)
      CALL CASIEX (D2D3, ELAIN, COFIBR, D, CASCOU, EPSEOUT, SIGMAI)
      CALL SIGTIL (D2D3, CASCOU, D, SIGMAI, SIGMAP)
      CALL NSIGEQ (D2D3, A2, SIGMAP, NORSII)
      VALF = NORSII-RPRE
      CALL CAHSIP (D2D3, A2, SIGMAP, HSIGTI)
      CALL CASIEX (D2D3, ELAIN, COFIBR, D, CASCOU, HSIGTI, KEDHSI)
      CALL SCAVEC (6, HSIGTI, KEDHSI, DIVIS)
      DIVIS = RP+(DIVIS/(NORSII*NORSII))
      DELP  = VALF/DIVIS
      MUL  = DELP/NORSII
      CALL MUMARE (MUL, 6, HSIGTI, DEPSPO)
CD    CALL IMPTDT ('VERIF PLASTO OUT ', DEPSPO, 1, 6)
      CALL COPITD (6, DEPSPO, EPSPO)
C 
5     CONTINUE
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Calcul de HSIGMA pour la plasticite
C 
      SUBROUTINE CAHSIP (D2D3, A2, SIGMA, HSIGMA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION   D2D3, A2
      DOUBLE PRECISION   SIGMA(6), HSIGMA(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      DOUBLE PRECISION   ALOC
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CAHSIP')
C 
C -----------------------------------------------------------------------
C 
      ALOC      = A2
      ALOC      = DSQRT(A2)
      HSIGMA(1) = 0.D0
CD    HSIGMA(2) = A2*SIGMA(2)
      HSIGMA(2) = A2*DMAX1(0.D0,SIGMA(2))
      HSIGMA(3) = SIGMA(3)
      HSIGMA(4) = 0.D0
CD    HSIGMA(5) = SIGMA(5)
      HSIGMA(5) = 0.D0
CD    HSIGMA(6) = A2*SIGMA(6)
      HSIGMA(6) = A2*DMAX1(0.D0,SIGMA(6))
CD    HSIGMA(6) = 0.D0
C 
C       d2d3 = 2  ---->  version ancienne (SIG12 et SIG22 seules c)
C       d2d3 = 3  ---->  nouvelle version (SIG12, SIG23, SIG33 et SIG22)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Calcul de la norme de SIGEQ en plasticite
C 
      SUBROUTINE NSIGEQ (D2D3, A2, SIGMA, NORSIE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION   NORSIE, A2, ALOC, D2D3
      DOUBLE PRECISION   SIGMA(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      CHARACTER*6        IDPROG
      PARAMETER         (IDPROG='NSIGEQ')
      DOUBLE PRECISION   HSIGMA(6)
C 
C 
CD    CALL WLKBCD(IDPROG)
C -----------------------------------------------------------------------
C 
      ALOC      = A2
      ALOC      = DSQRT(A2)
      HSIGMA(1) = 0.D0
CD    HSIGMA(2) = A2*SIGMA(2)
      HSIGMA(2) = A2*DMAX1(0.D0,SIGMA(2))
      HSIGMA(3) = SIGMA(3)
      HSIGMA(4) = 0.D0
CD    HSIGMA(5) = SIGMA(5)
      HSIGMA(5) = 0.D0
CD    HSIGMA(6) = A2*SIGMA(6)
      HSIGMA(6) = A2*DMAX1(0.D0,SIGMA(6))
CD    HSIGMA(6) = 0.D0
      CALL SCAVEC (6, HSIGMA, HSIGMA, NORSIE)
      NORSIE = DSQRT(NORSIE)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Derivee du seuil par rapport a la plasticite
C 
      SUBROUTINE SEUILP (AL, BE, P, RP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION  AL, BE
      DOUBLE PRECISION  P, RP
C 
C -----------------------------------------------------------------------
C 
      RP = BE*AL*(P**(AL-1.D0))
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Calcul de la contrainte equivalente pour la plasticite
C 
      SUBROUTINE SIGTIL (D2D3, CASCOU, D, SIGMA, STILDE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION   D(3)
      LOGICAL            CASCOU(3)
      DOUBLE PRECISION   D2D3
      DOUBLE PRECISION   SIGMA(6), STILDE(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      DOUBLE PRECISION   DPS, DPT
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CAHSIP')
C 
CD    CALL WLKBCD(IDPROG)
C -----------------------------------------------------------------------
C 
      DPS   = D(2)
      DPT   = D(3)
C 
CD    STILDE(1) = DMAX1(0.D0,SIGMA(1))
      STILDE(1) = SIGMA(1)
CD    IF (CASCOU(2)) THEN
CD      STILDE(2) = SIGMA(2)/(1.D0-DPT)
CD    ELSE
CD      STILDE(2) = DMAX1(0.D0,SIGMA(2))
        STILDE(2) = SIGMA(2)
CD    END IF
CD    STILDE(3) = SIGMA(3)/(1.D0-DPS)
      STILDE(3) = SIGMA(3)
      STILDE(4) = SIGMA(4)
CD    STILDE(5) = SIGMA(5)/(1.D0-DPS)
      STILDE(5) = SIGMA(5)
CD    IF (CASCOU(3)) THEN
CD      STILDE(6) = SIGMA(6)/(1.D0-DPS)
CD    ELSE
CD      STILDE(6) = DMAX1(0.D0,SIGMA(6))
        STILDE(6) = SIGMA(6)
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
