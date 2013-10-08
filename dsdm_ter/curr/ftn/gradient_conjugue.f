C     Cette routine verifie l'orthogonalisation des nb deformations
C 
C     La verification  est a prendre au sens de K0OU de K0-1
C 
C     Les deplacements sont aussi verifies au sens de K0
C 
C     ADCOU est l'adresse de depart du tableau :
C 
C       HOO-COUCHE => NORME AU SENS DE K0
C       SOU-COUCHE => NORME AU SENS DE K0-1
C 
C     ADINT est l'adresse de depart du tableau :
C 
C       HOO-INTERF => NORME AU SENS DE K0
C       SOU-INTERF => NORME AU SENS DE K0-1
C 
C     On envoie comme arguments :
C 
C     E ...... NB       nombre de vecteurs 
C     E....... LONEPS   longueur d'une deformation
C     E....... LONSAU   longueur d'un saut
C     E....... VEC...   vecteurs a verifier
C 
      SUBROUTINE VEORES (ADCOU, ADINT, NB, LONEPS, LONSAU, EPS, SAU)
C 
C ----------------------------------------------------------------------- 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER       ADCOU, ADINT
      INTEGER       LONEPS, LONSAU, NB
C 
      DOUBLE PRECISION  EPS(LONEPS*NB)
      DOUBLE PRECISION  SAU(LONSAU*NB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  J, K
      INTEGER  DBEPS, DBSAU
      INTEGER  KDEBUE, KDEBUS
C  
      DOUBLE PRECISION NORMF, PSCALV
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VEORES')
C 
C -----------------------------------------------------------------------
      CALL IMPET ('NOMBRE DE VECTEURS A VERIFIER '//IDPROG, NB)
C 
      DBEPS = 1
      DBSAU = 1
C 
      DO J = 1, NB 
        CALL NOKGLO (ADCOU, ADINT, 1, 1, EPS(DBEPS), SAU(DBSAU), NORMF)
        CALL IMPET ('POUR LA DEFORMATION  '//IDPROG, J) 
        NORMF = DSQRT(NORMF) 
        CALL IMPDT ('VALEUR DE SA NORME      '//IDPROG, NORMF)
        KDEBUE    = 1
        KDEBUS    = 1	
	DO K = 1, J-1
           CALL  SCKGLO (ADCOU, ADINT, 1, 1,
     &                   EPS(DBEPS), SAU(DBSAU), 1,
     &                   EPS(KDEBUE), SAU(KDEBUS), PSCALV)
C 
C         CALL IMPET ('PRODUIT SCALAIRE AVEC L''ANCIEN VECTEUR ', K) 
C         CALL IMPDT ('VALEUR DU PRODUIT SCALAIRE ', PSCALV)
C 
          KDEBUE    = KDEBUE + LONEPS
          KDEBUS    = KDEBUS + LONSAU	
	END DO 
        DBEPS = DBEPS+LONEPS
        DBSAU = DBSAU+LONSAU	
      END DO
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine Verifie l'ORthogonalisation de Toutes les deformations(EP) 
C     admissibles a zero
C 
      SUBROUTINE VORTEP (TYPE)
C 
C     TYPE = 0 VERIFICATION PAR RAPPORT A K0
C     TYPE = 1 VERIFICATION PAR RAPPORT A K0-1       
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'strategie_calcul.h'
      integer TYPE
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  NBCHA
      INTEGER  AM2LC, ADM2LC, FTDEOP
      INTEGER  FTECH6, EPSCH8, SAUCH9
      INTEGER  DEPCH0, LONDEP, EPSTRA, SAUTRA, DEPTRA
      INTEGER  NUEPS, NBCPOR
      INTEGER  DEPSTR, DSAUTR, DDEPTR
      INTEGER  LONDER, NNOUV, LONEPS, LONSAU, NNOUVA
      INTEGER  FTREE7, COEFF
      INTEGER  ADINT, ADCOU
      INTEGER  TABNIE(2), TABNIS(2), TABNID(2), FOTDEP
C 
C     LONGUEURS DE TABLEAUX DE TRAVAIL
C 
      INTEGER  TRAV1, ADINTE
C 
C     ADRESSE DE TRAVAIL
C 
      INTEGER ADREEL, DBCHEP
C 
      DOUBLE PRECISION DIVIS
C 
C     POUR LES VERIFS
C 
CD    DOUBLE PRECISION NORVER
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VORTEP')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
D       IF( TYPE .NE. 0 .AND. TYPE .NE. 1) THEN
D         CALL IMPET ('ERREUR TYPE = '//IDPROG , TYPE)
D         STOP
D       END IF
C 
      LONEPS = NEPS*NGAU1*NTETA
      LONSAU = NSAU*NGAU2*NTETA
C 
      IF (TYPE .EQ. 0) THEN
        CALL MESSAO ('VERIFICATION POUR TOUS LES CHAMPS CA 0')
        CALL ADTBDM ('HOO-COUCHE', ADCOU)
C 
C       Tableaux des deplacements, deformations et contraintes admissibles
C 	 
	CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
        IF (NBINT .GT. 0) THEN
          CALL ADTBDM ('HOO-INTERF', ADINT)
	  CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
        ELSE
          ADINT = ADCOU
	  SAUCH9 = EPSCH8
        END IF	 
      ENDIF 
C 
      IF (TYPE .EQ. 1) THEN
        CALL MESSAO ('VERIFICATION POUR TOUS LES CHAMPS SA A 0')      
        CALL ADTBDM ('SOU-COUCHE', ADCOU)
C 
C       Tableaux des deplacements deformations et contraintes admissibles
C 	 
	CALL ADTBDM ('SIG-AD-TOT', EPSCH8)
        IF (NBINT .GT. 0) THEN
          CALL ADTBDM ('SOU-INTERF', ADINT)
	  CALL ADTBDM ('SGN-AD-TOT', SAUCH9)
        ELSE
          ADINT = ADCOU
	  SAUCH9 = EPSCH8
        END IF	 
      ENDIF 
C 
      TABNIE(1) = CHARAX
      TABNIE(2) = LONEPS
C 
      TABNIS(1) = CHARAX
      TABNIS(2) = LONSAU
C 
C     On extrait du tableau des deformations totales stockees (charax, loneps) les
C     deformations allant du numero 2 a DEADTR.
C 
      NBCHA = DEADTR-1
      TRAV1  =  NBCHA*(LONEPS+LONSAU)
C 
      CALL POUSMD (TRAV1, EPSTRA)
C 
      SAUTRA  = EPSTRA + NBCHA*LONEPS
C 
      DEPSTR = EPSTRA
      DSAUTR = SAUTRA
C 
      DO NUEPS = 2 , DEADTR
C 
        CALL EXTRAD (DM(EPSCH8), 2, TABNIE(1),
     &               2, NUEPS, DM(DEPSTR), LONEPS)
C 
        DEPSTR = DEPSTR+LONEPS
C 
       IF (NBINT .GT. 0) THEN
C 
          CALL EXTRAD (DM(SAUCH9), 2, TABNIS(1),
     &                 2, NUEPS, DM(DSAUTR), LONSAU)
C 
          DSAUTR = DSAUTR+LONSAU
C 
        END IF
C 
      END DO
C 
C     verification de l'Orthogonalisation
C 
      CALL VEORES (ADCOU, ADINT, NBCHA, 
     &             LONEPS, LONSAU,
     &             DM(EPSTRA), DM(SAUTRA))
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
      RETURN
      END      
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine a pour but de determiner la meilleure contrainte
C     admissible elements finis tout en etant proche de la resolution
C     avec la matrice tangente par gradient conjugue.
C 
C     On calcule la solution DEPSOL cinematiquement admissible
C     a zero pour l'iteration en espace concernee de l'etape globale a savoir :
C 
C     Intemps (fotemp*K0*fotemp)[DEPSOL] = Intemps (fotemp*-*DELTA(sigchap))
C 
C     On envoie comme arguments :
C 
C     E ...... NITEMX  nombre d'iterations maxi du gradient conjugue
C     E ...... REPRIS  logique indiquant si il y a eu une iteration en temps
C     ES...... ADSMEG  l'adresse de depart des seconds membres;
C                      pour l'etape globale ceux-ci sont assembles dans DM
C                      a partir de ADSMEG de la facon suivante (NDDL, NFTGLO, NBMAT)
C 
C     Les seconds membres doivent etre mis a zero en entree dm(adsmeg, ...) = 0.d0
C 
C     E ...... NBDEV   nombre de developpements non nuls
C     E ...... TNUDEV  tableau des developpements non nuls

C     Et on recupere :
C 
C     S ...... DEPSOL  solution developpee Fourier en deplacement de ce probleme
C     S ...... EPSSOL  solution en deformation de ce probleme
C     S ...... SAUSOL  solution en saut de ce probleme
C     S ...... SIGRES  contrainte residuelle couche calculee en isotrope
C     S ...... SGNRES  contrainte residuelle interface calculee en isotrope
C 
C     Attention epssol, depsol, sausol, sigres, sgnres doivent etre mise a zero
C 
      SUBROUTINE GCADMI (NITEMX, REPRIS, ADSMEG, NBDEV, TNUDEV,
     &                   DEPSOL, EPSSOL, SAUSOL, SIGRES, SGNRES)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER          NITEMX, NBDEV, TNUDEV(NBMAT)
C 
      DOUBLE PRECISION DEPSOL(NDDL*NBMAT)
      DOUBLE PRECISION EPSSOL(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION SAUSOL(NSAU*NTETA*NGAU2)
      DOUBLE PRECISION SIGRES(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION SGNRES(NSAU*NTETA*NGAU2)
C 
      LOGICAL REPRIS
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER ADSMEG, AM2LC, ADM2LC
      INTEGER NUITER, DEPOPT, EPSOPT
      INTEGER AMPLC, ADMPLC, SAUOPT
      INTEGER DEPPLA, NDEPPL, EPSILO
      INTEGER ADSEN, NEPSIL, SAUT, NSAUT
      INTEGER SIGITE, SINITE
      INTEGER ADLADN, EPLADN, SALADN
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      DOUBLE PRECISION NORACT, PRECIS, SCAL
      DOUBLE PRECISION NORIN, DNKDN, DNRN, MUN, DNKSN, LANDAN
C 
C     Pour les longueurs des tableaux provisoires
C 
      INTEGER LONEPS, LONSAU, LONDEP, LONDER, LONRES
      INTEGER ADITER, ADNORM, DBITER, DBNORM
      INTEGER ADCHOO, ADIHOO
      INTEGER AEPSDN, ASAUDN
C 
      DOUBLE PRECISION VZERO
C 
C     Pour test
C 
      INTEGER RESU, DEBUT, I, ARESU, RESUD, ARESUD
C 
C     Pour les verifications sur la solution
C 
CD    INTEGER  SIGSOL, SGNSOL
CD    DOUBLE PRECISION TRAVE1, TRAVE2, TRALOV
C 
      PARAMETER (PRECIS = 1.D -6)
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='GCADMI')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C 
C     Pour verification si non-convergence
C 
      SCAL =  0.D0
C 
      DEBUT = ADSMEG
      CALL POUSMD (NBMAT, ARESUD)
      RESUD = ARESUD
C 
      DO I = 1, NTDSFG
C 
        CALL SCAVEC (NDDL, DM(DEBUT), DM(DEBUT), DM(RESUD))
        DM(RESUD) = PI*DM(RESUD)
        SCAL  = SCAL + DM(RESUD)
        RESUD = RESUD +1
        DEBUT = DEBUT + NDDL
C 
      END DO
C 
      CALL SCAVEC (NDDL, DM(DEBUT), DM(DEBUT), DM(RESUD))
      DM(RESUD) = 2.D0*PI*DM(RESUD)
      SCAL  = SCAL+ DM(RESUD)
      RESUD = RESUD +1
      DEBUT = DEBUT + NDDL
C 
      DO I = NTDSFG+2, NBMAT
C 
        CALL SCAVEC (NDDL, DM(DEBUT), DM(DEBUT), DM(RESUD))
        DM(RESUD) = PI*DM(RESUD)
        SCAL  = SCAL+ DM(RESUD)
        RESUD = RESUD +1
        DEBUT = DEBUT + NDDL
C 
      END DO
C 
      DO I = ARESUD, ARESUD+NBMAT-1
C 
        DM(I) = DM(I)/SCAL
C 
      END DO
C 
C     fin de verification si non-convergence
C 
      CALL GSPOUD (2*(NIBFEM+1), ADITER)
C 
      CALL ADTBDM ('HOO-COUCHE', ADCHOO)
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL ADTBDM ('HOO-INTERF', ADIHOO)
C 
      ELSE
C 
        ADIHOO = ADCHOO
C 
      END IF
C 
      ADNORM = ADITER+NIBFEM+1
      DBITER = ADITER
      DBNORM = ADNORM
C 
      LONEPS = NTETA*NEPS*NGAU1
      LONDER = NDDL*NTETA
      LONDEP = NDDL*NBMAT
      LONSAU = NTETA*NSAU*NGAU2
C 
C     RESERVATION DE PLACE POUR LES NOUVEAUX DEPLACEMENTS ET DEFORMATIONS
C 
      LONRES = 2*(LONDEP + LONEPS +LONSAU)
      CALL POUSMD (LONRES, DEPPLA)
C 
      NDEPPL = DEPPLA
      EPSILO = DEPPLA+LONDEP
      NEPSIL = EPSILO
      SAUT   = EPSILO+LONEPS
      NSAUT  = SAUT
C 
C     RESERVATION DE PLACE POUR LES TERMES CALCULES
C 
      DEPOPT = SAUT   + LONSAU
      EPSOPT = DEPOPT + LONDEP
      SAUOPT = EPSOPT + LONEPS
C 
      CALL SCADEV (DM(ADSMEG), DM(ADSMEG), NORIN)
C 
      NORIN  = DSQRT(NORIN)
C 
      DM(DBNORM)  = NORIN
      DBNORM      = DBNORM+1
      DM(DBITER)  = 0.D0
      DBITER      = DBITER+1
C 
      CALL IMPDT ('NORME INITIALE DE L''EFFORT ', NORIN)
      CALL POUSMD (LONDEP, ADSEN)
      CALL COPID (DM(ADSMEG), DM(ADSEN), LONDEP, 1)
C 
C     Calcul de la direction de recherche initiale, en deplacement, deformation et saut
C 
      IF (.NOT. REPRIS) THEN
C 
C       Si c'est le 1er passage pour l'iteration  actuelle on
C       calcule la direction de recherche initiale par Ko
C 
        CALL CALIS1 (DM(ADSMEG), NBDEV, TNUDEV(1), DM(NDEPPL))
C 
D       CALL VDPVZE (DM(NDEPPL))
C 
CD      CALL IMPTDT ('DEPLACEMENT : DIRECTION INITIALE',
CD   &                DM(NDEPPL), NDDL, NBMAT)
C 
        CALL TOUEPS (NFTGLO, DM(NDEPPL), DM(NEPSIL))
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL TOUSAU (NFTGLO, DM(NDEPPL), DM(NSAUT))
C 
        END IF
C 
      ELSE
C 
C       Sinon on choisit comme direction de recherche initiale
C       la solution de l'iteration en espace precedente, puis
C       on remet a zero les solutions en deplacements, contraintes
C       et deformations.
C 
        CALL COPID (DEPSOL, DM(NDEPPL), LONDEP, 1)
        CALL BALAI (DEPSOL, LONDEP, 1)
        CALL COPID (EPSSOL, DM(NEPSIL), LONEPS, 1)
        CALL BALAI (EPSSOL, LONEPS, 1)
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL COPID (SAUSOL, DM(NSAUT), LONSAU, 1)
          CALL BALAI (SAUSOL, LONSAU, 1)
C 
        END IF
C 
      END IF
C 
      NUITER = 1
C 
C ***********************************************************************
C 
C     DEBUT DU DO WHILE DES ITERATIONS
C 
C ***********************************************************************
C 
C     On modifie car il faut avoir une solution correspondant
C     vraiment a la direction de recherche
C 
      DO WHILE (NUITER .LT. 30)
C 
C     DO WHILE (NUITER .LT. NITEMX)
C 
1001    CONTINUE
C 
        CALL ENPOUB (AMPLC, ADMPLC)
C                     _
C       CALCUL DE dnKdn
C 
        CALL NOKGLO (ADCHOO, ADIHOO, NFTGLO, NFTGLO,
     &               DM(NEPSIL), DM(NSAUT), DNKDN)
C 
C       CALCUL DE dnRn
C 
        CALL SCADEV (DM(ADSMEG), DM(NDEPPL), DNRN)
CD      CALL IMPDT ('dN*RN ', DNRN)
C 
C       CALCUL DE MUN
C 
        MUN = DNRN/DNKDN
C 
C       CALCUL DE LA DIRECTION OPTIMALE Dn TELLE QUE Dn =  MUn*dn
C       ROUTINE HOMAT INTERPRETEUR-GESDYN : MULTIPLICATION DE MATRICE
C 
        CALL HOMAT (MUN, DM(NDEPPL), DM(NDEPPL), LONDEP, 1)
        CALL HOMAT (MUN, DM(NEPSIL), DM(NEPSIL), LONEPS, 1)
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL HOMAT (MUN, DM(NSAUT), DM(NSAUT), LONSAU, 1)
C 
        END IF
C                     _
C       CALCUL DE DnKDn
C 
        DNKDN = MUN*MUN*DNKDN
C 
C       CALCUL DE LA NOUVELLE SOLUTION Un+1 TELLE QUE :
C       Un+1 = Un+ MUn*dn = Un+Dn
C       ROUTINE ADD INTERPRETEUR-GESDYN : ADDITION DE DEUX MATRICES
C 
        CALL ADD (DM(NDEPPL), DEPSOL, DEPSOL, LONDEP, 1)
        CALL ADD (DM(NEPSIL), EPSSOL, EPSSOL, LONEPS, 1)
C 
        IF (NBINT .GT. 0) THEN
          CALL ADD (DM(NSAUT), SAUSOL, SAUSOL, LONSAU, 1)
        END IF
C 
C       CALCUL DU NOUVEAU RESIDU Rn+1 TEL QUE : Rn+1 = F-K(Un+1) = Rn-K(Dn)
C 
        LONRES = LONEPS+LONSAU
C 
        CALL GSPOUD (LONRES, SIGITE)
C 
        SINITE = SIGITE + LONEPS
C 
        CALL NOUSIG (NFTGLO, DM(NEPSIL), DM(NSAUT),
     &               DM(SIGITE), DM(SINITE))
C 
        CALL NOUEFF (NFTGLO, DM(SIGITE), DM(SINITE), ADSMEG,
     &               NBDEV, TNUDEV(1), NORACT)
C 
CD      CALL IMPDT ('NORACT/NORIN ', NORACT/NORIN)
C 
        DM(DBNORM) = NORACT
        DBNORM     = DBNORM+1
        DM(DBITER) = DBLE(NUITER)
        DBITER     = DBITER+1
C 
CD      CALL IMPDT ('VALEUR DE LA NORME DU RESIDU ', NORACT)
C 
C       Test de convergence sur le nouveau residu
C 
        IF  ((NORACT/NORIN .LT. PRECIS) .OR. (NUITER .EQ. NIBFEM)) THEN
C 
C         Verification de la coherence entre nokglo et traglo et nousig
C 
          CALL IMPDT ('VALEUR DE PRECISION DANS '//IDPROG, PRECIS)
          CALL IMPDT ('VALEUR DE  NORACT/NORIN  '//IDPROG,
     &                 NORACT/NORIN)
C 
          CALL IMPET ('VALEUR DE NUITER '//IDPROG, NUITER)
C 
CD         LONRES = LONEPS+LONSAU
C 
CD         CALL GSPOUD(LONRES , SIGSOL )
CD         SGNSOL = SIGSOL+LONEPS
C 
CD         CALL NOUSIG (NFTGLO, EPSSOL, SAUSOL, DM(SIGSOL), DM(SGNSOL))
C 
CD         CALL TRAGLO (EPSSOL, SAUSOL, DM(SIGSOL), DM(SGNSOL),
CD   &                  TRAVE1, TRALOV)
C 
CD         CALL IMPDT ('Norme en Energie de la solution par TRAGLO '
CD   &                  // IDPROG, DSQRT(TRAVE1))
C 
CD         CALL NOKGLO (ADCHOO, ADIHOO, NFTGLO, NFTGLO, EPSSOL, SAUSOL,
CD   &                  TRAVE2)
C 
CD         CALL IMPDT ('Norme en Energie de la solution par NOKGLO '
CD   &                  //IDPROG, DSQRT(TRAVE2) )
C 
CD         CALL SCADEV (DM(ADSEN), DEPSOL, TRALOV)
C 
CD         CALL IMPDT ('Norme en Energie de la solution par F U '
CD   &                  //IDPROG, DSQRT(TRALOV))
C 
C ***********************************************************************
C 
C          ON FORCE LA FIN DU DO WHILE
C 
C ***********************************************************************
           GOTO 1000
C 
C ***********************************************************************
C 
C       FIN DU TRAITEMENT POUR LA CONVERGENCE
C 
C ***********************************************************************
        END IF
C 
C       Calcul de la nouvelle solution sn+1 telle que : Ko (sn+1) = Rn+1
C 
        CALL CALIS1 (DM(ADSMEG), NBDEV, TNUDEV(1), DM(DEPOPT))
C 
CD      CALL MESSAO ('2EME APPEL A VDPVZE DANS '//IDPROG)
CD      CALL VDPVZE (DM(DEPOPT))
C 
C       Calcul des deformations admissibles a zero rangees (NEPS, NTETA, NGAU1)
C 
        CALL TOUEPS (NFTGLO, DM(DEPOPT), DM(EPSOPT))
C 
C       Calcul des sauts admissibles a zero  ranges (NSAU, NTETA, NGAU2)
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL TOUSAU (NFTGLO, DM(DEPOPT), DM(SAUOPT))
C 
        END IF
C 
C       Calcul de DnKsn=1
C 
        CALL SCKGLO (ADCHOO, ADIHOO, NFTGLO, NFTGLO,
     &               DM(NEPSIL), DM(NSAUT), NFTGLO,
     &               DM(EPSOPT), DM(SAUOPT), DNKSN)
C 
C       Calcul de landan = - DnKsn/DnKDn
C 
        LANDAN = - DNKSN/DNKDN
CD      CALL IMPDT ('LANDAN = - DN*K*SN/DNKDN ', LANDAN)
CD 
CD      Orthogonalisation de la nouvelle direction sn+1
CD      par rapport a l' ancienne direction pour K ****
CD 
CD      1)    ***** CALCUL DE Dn+1 *****
CD 
CD            ***** En deplacement ****
CD 
CD        POUR SAUVEGARDER DN POUR LA VERIFICATION
CD 
CD      LONRES = LONEPS+LONSAU
CD      CALL POUSMD( LONRES , AEPSDN )
CD      ASAUDN = AEPSDN+LONEPS
CD      CALL COPITD(LONEPS ,DM(NEPSIL),DM(AEPSDN))
CD      CALL COPITD(LONSAU ,DM(NSAUT),DM(ASAUDN))
CD 
        LONRES = LONDEP+LONEPS+LONSAU
        CALL POUSMD (LONRES, ADLADN)
        EPLADN = ADLADN + LONDEP
        SALADN = EPLADN + LONEPS
C 
C             ***** En deplacement ****
C 
        CALL HOMAT (LANDAN, DM(NDEPPL), DM(ADLADN), LONDEP, 1)
        CALL ADD (DM(ADLADN), DM(DEPOPT), DM(NDEPPL), NDDL, NBMAT)
C 
C             ***** En deformation ****
C 
        CALL HOMAT (LANDAN, DM(NEPSIL), DM(EPLADN), LONEPS, 1)
        CALL ADD (DM(EPLADN), DM(EPSOPT), DM(NEPSIL), LONEPS, 1)
C 
C             ***** En saut ****
C 
        IF (NBINT .GT. 0) THEN
C 
          CALL HOMAT (LANDAN, DM(NSAUT), DM(SALADN), LONSAU, 1)
          CALL ADD (DM(SALADN), DM(SAUOPT), DM(NSAUT), LONSAU, 1)
C 
        END IF
C 
C       Verification de l'orthogonalite de DN et de DN+1
C       Steph. Le calcul plante 'segmentation violation' 
C       quand j'active SCKGLO, dans gcnlin/gcadmi/sckglo/copitd
C 
CD      CALL SCKGLO (ADCHOO, ADIHOO, NFTGLO, NFTGLO, DM(NEPSIL ),
CD   &               DM(NSAUT), NFTGLO, DM(AEPSDN), DM(ASAUDN), VZERO)
CD      CALL IMPDT ('VALEUR DU PRODUIT SCALAIRE APRES ORTHOGONALISATION'
CD   &               //' NULLE DANS '//IDPROG, VZERO)
C 
C       Remise a zero des tableaux partiels de l'iteration
C 
        NUITER = NUITER +1
        CALL SOPOUB(AMPLC,ADMPLC)
C 
C ***********************************************************************
C 
C     FIN DU DO WHILE DES ITERATIONS
C 
C ***********************************************************************
C 
      END DO
C 
C     Si on sort par depassement du nombre d'iterations
C 
      CALL IMPET ('SORTIE PAR NOMBRE D''ITERATIONS DANS '// IDPROG,
     &             NUITER)
      CALL IMPDT ('NORACT DANS '//IDPROG, NORACT)
      CALL IMPDT ('NORIN  DANS '//IDPROG, NORIN)
      CALL IMPDT ('PRECIS dDANS '//IDPROG, PRECIS)
C 
      SCAL =  0.D0
C 
      DEBUT = ADSMEG
      CALL POUSMD (NBMAT, ARESU)
      RESU = ARESU
C 
      DO I = 1 , NTDSFG
C 
        CALL SCAVEC (NDDL, DM(DEBUT), DM(DEBUT), DM(RESU))
C 
        DM(RESU) = PI*DM(RESU)
        SCAL = SCAL+ DM(RESU)
        RESU = RESU+1
        DEBUT = DEBUT + NDDL
C 
      END DO
C 
      CALL SCAVEC (NDDL, DM(DEBUT), DM(DEBUT), DM(RESU))
      DM(RESU) = 2.D0*PI*DM(RESU)
      SCAL  = SCAL+ DM(RESU)
      RESU  = RESU+1
      DEBUT = DEBUT + NDDL
C 
      DO I = NTDSFG+2 , NBMAT
C 
        CALL SCAVEC (NDDL, DM(DEBUT), DM(DEBUT), DM(RESU))
        DM(RESU) = PI*DM(RESU)
        SCAL  = SCAL+ DM(RESU)
        RESU  = RESU+1
        DEBUT = DEBUT + NDDL
C 
      END DO
C 
      DO I = ARESU, ARESU+NBMAT-1
C 
        DM(I) = DM(I)/SCAL
C 
      END DO
C 
      CALL IMPTDT ('NORME DE L''EFFORT PAR DEVELOPPEMENT ',
     &              DM(ARESUD), 1, NBMAT)
C 
CD    CALL IMPTDT ('NORME DU RESIDU PAR DEVELOPPEMENT ',
CD   &              DM(ARESU), 1, NBMAT)
C 
C     Sortie du do while si le residu est satisfaisant
C 
1000  CONTINUE
C 
      CALL IMPET ('ON SORT A  L''ITERATION ', NUITER)
      CALL IMPDT ('POUR UNE VALEUR PRECIS '// IDPROG, PRECIS)
C 
      CALL IMPDT ('POUR UNE VALEUR NORACT/NORIN<PRECIS? ', NORACT/NORIN)
C 
      CALL OMPTDN ('DEPLACEMENT SOLUTION DEVELOPPEE ', DEPSOL, NDDL,
     &              NBMAT)
C 
C     Calcul de la contrainte associee au residu pour construire
C     une contrainte rigoureusement SA a 0
C 
      CALL CALIS1 (DM(ADSMEG), NBDEV, TNUDEV(1), DM(DEPOPT))
C 
CD    CALL MESSAO ('3EME APPEL A VDPVZE DANS '//IDPROG)
CD    CALL VDPVZE (DM(DEPOPT))
C 
C     Calcul des deformations admissibles a zero ranges (NEPS, NTETA, NGAU1)
C 
      CALL TOUEPS (NFTGLO, DM(DEPOPT), DM(EPSOPT))
C 
C     Calcul des sauts admissibles a zero ranges (NSAU, NTETA, NGAU2)
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL TOUSAU (NFTGLO, DM(DEPOPT), DM(SAUOPT))
C 
      END IF
C 
      CALL SIGISO (NFTGLO, DM(EPSOPT), DM(SAUOPT), SIGRES, SGNRES)
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
C     ANCIENNEMENT AOUSIG => VOIR NOUVELLE DONNEE
C 
C     !!!!! ATTENTION !!!! necessite de mettre les contraintes a zero
C 
C     Cette routine calcule la contribution de la nouvelle deformation NUEPS
C     point de Gauss par point de Gauss a la nouvelle contrainte SIGNOU
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions du temps de l'iteration
C     E ...... TOUEPS   deformation rangee (NEPS*NTETA*NGAU1*NFONCT)
C     E ...... TOUSAU   les sauts ranges (NSAU*NTETA*NGAU1*NFONCT)
C 
C     Et on recupere :
C 
C     S ...... SIGNOU   nouvelles contraintes rangees
C                       (NEPS,NTETA,NGAU1 ,NFONCT)
C     S ...... SGNNOU   nouvelles contraintes normales rangees
C                       (NSAU,NTETA,NGAU2,NFONCT)

      SUBROUTINE NOUSIG (NFONCT, TOUEPS, TOUSAU, SIGNOU, SGNNOU)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NFONCT
C 
      DOUBLE PRECISION TOUEPS(NEPS*NTETA*NGAU1*NFONCT)
      DOUBLE PRECISION SIGNOU(NEPS*NTETA*NGAU1*NFONCT)
      DOUBLE PRECISION TOUSAU(NSAU*NTETA*NGAU2*NFONCT)
      DOUBLE PRECISION SGNNOU(NSAU*NTETA*NGAU2*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  NUCOU, NUCOL, TETA, X, Y
      INTEGER  NUINT, NUSAU
      INTEGER  DEPS, NUEPS, DEBSIG, LONEC
      INTEGER  DSAU, DEBSGN
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      DOUBLE PRECISION KLOC(17), COUISO(10)
      DOUBLE PRECISION KLOCI(9), INTISO(3)
      DOUBLE PRECISION TETORT, TETCAL
C 
      INTEGER          ADTETA
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NOUSIG')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
C 
      DEPS = 1
C 
C     BOUCLE SUR LES FONCTIONS DU TEMPS
C 
      DO NUEPS = 1, NFONCT
C 
        DEBSIG = 1
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C       BOUCLE i SUR LES COUCHES
C 
        DO NUCOU = 1 , NBCOU
C 
          CALL ANGCOU (NUCOU, TETORT)
          CALL COCISO (NUCOU, COUISO)
C 
C         BOUCLE ii SUR LES COLONNES
C 
          DO NUCOL = 1, NBCOL
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT X
C 
            DO X = 1, XINTEG
C 
C             BOUCLE iv SUR LES POINTS DE GAUSS SUIVANT Y
C 
              DO Y = 1, YINTEG
C 
C               BOUCLE v SUR LES ANGLES
C 
                DO TETA = 1, NTETA
C 
C                 Recherche de l'angle de la bande correspondant a teta
C 
                  TETCAL = DM (ADTETA+TETA-1)
C 
C                 On recupere les quantites stockees dans Q-CHAPEAU
C                 pour le point de gauss, pour toutes les matrices
C 
                  CALL RICLOC (COUISO, TETORT, TETCAL, KLOC)
                  CALL CTSIGC (NFONCT, NUEPS, TOUEPS(DEPS),
     &                         KLOC, DEBSIG, SIGNOU, LONEC)
                  DEPS   = DEPS  +NEPS
                  DEBSIG = DEBSIG+NSIG
C 
C               FIN DE BOUCLE v SUR LES ANGLES
C 
                END DO
C 
C             FIN DE BOUCLE iv SUR LES POINTS DE GAUSS SUIVANT z
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES COLONNES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES COUCHES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES FONCTIONS DU TEMPS
C 
      END DO
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
      DSAU   = 1
C 
C     BOUCLE SUR LES FONCTIONS DU TEMPS
C 
      DO NUSAU = 1, NFONCT
C 
        DEBSGN = 1
C 
C       BOUCLE i SUR LES INTERFACES
C 
        DO NUINT = 1 , NBINT
C 
        CALL COIISO (NUINT, INTISO)
        CALL ANGINT (NUINT, TETORT)
C 
C         BOUCLE ii SUR LES COLONNES
C 
          DO NUCOL = 1, NBCOL
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            DO X = 1, XINTEG
C 
C             BOUCLE iv SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
C               RECHERCHE DE L'ANGLE DE LA BANDE (tetcal) CORRESPONDANT A teta
C 
                TETCAL = DM(ADTETA+TETA-1)
C 
C               OM RECUPERE LES QUANTITES STOCKEES DANS LE DIRECTORY Q-CHAPEAU
C               POUR LE POINT DE GAUSS CONSIDERE POUR TOUS LES TEMPS
C 
                CALL INTELA (INTISO, TETORT, TETCAL, KLOCI)
                CALL CTSIGI  (NFONCT, NUSAU, TOUSAU(DSAU),
     &                        KLOCI, DEBSGN, SGNNOU, LONEC)
                DSAU    = DSAU   + NSAU
                DEBSGN  = DEBSGN + NSAU
C 
C             FIN DE BOUCLE iv SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES COLONNES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES INTERFACES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES FONCTIONS DU TEMPS
C 
      END DO
C 
C 
      IF (NBINT .GE. 1) THEN
        CALL TESTEN (DEBSGN-1, NSAU*NTETA*NGAU2, IDPROG)
        CALL TESTEN (DSAU-1, NSAU*NTETA*NGAU2*NFONCT, IDPROG)
        CALL TESTEN (LONEC, NSAU*NTETA*NGAU2*NFONCT, IDPROG)
      END IF
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
C     E ...... ISODEV   Caracteristiques du comportement de la
C                       couche dans sa base d'orthotropie
C 
C           |K11  K12   0   0   0  A1|  E11
C           |     K22   0   0   0  A2|  E22
C           |          K66  0   0   0|  R2*E12
C           |               B1  0   0|  R2*E23  =====>DEFORMATIONS
C           |                   B2  0|  R2*E13
C           |   SYM                 C|  E33
C 
C     E ...... TETORT   Valeur de l'angle que fait cette base avec
C                       le repere global
C     E ...... TETCAL   Valeur de l'angle pour lequel on calcule
C                       le comportement
C 
C     Et on recupere :  K(17) STOCKE
C 
C     S ...... K(3,3)   Partie du comportement travaillant sur
C                       la partie plane des deformations
C     S ...... B(2,2)   Partie du comportement travaillant sur
C                       le cisaillement normal
C     S ...... A(3)     Partie du comportement reliant
C                       la partie plane des deformations a la
C                       deformation normale
C     S ...... C        Relie la contrainte a la deformation normale
C 
      SUBROUTINE RICLOC (ISODEV, TETORT, TETCAL, K)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  ISODEV(10) , TETORT , TETCAL , K(17)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  R2, COS2, COS4, SIN2, SIN4, TETLOC
C 
CD    LOGICAL LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='RICLOC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      R2 = DSQRT(2.D0)
C 
      TETLOC = TETCAL-TETORT
      COS2   = DCOS(2.D0*TETLOC)
      SIN2   = DSIN(2.D0*TETLOC)
      COS4   = DCOS(4.D0*TETLOC)
      SIN4   = DSIN(4.D0*TETLOC)
C 
C     CALCUL DE K
C 
      K(1)  = ISODEV(1)+COS2*ISODEV(7)+COS4*ISODEV(10)
      K(2)  = ISODEV(2)-COS4*ISODEV(10)
      K(3)  = -R2*(SIN2*ISODEV(7)/2.D0+SIN4*ISODEV(10))
      K(4)  = K(2)
      K(5)  = ISODEV(1)-COS2*ISODEV(7)+COS4*ISODEV(10)
      K(6)  = -R2*(SIN2*ISODEV(7)/2.D0-SIN4*ISODEV(10))
      K(7)  = K(3)
      K(8)  = K(6)
      K(9)  = ISODEV(3)-2.D0*COS4*ISODEV(10)
C 
C     CALCUL DE B
C 
      K(10)  = ISODEV(5)+COS2*ISODEV(9)
      K(11)  = SIN2*ISODEV(9)
      K(12)  = K(11)
      K(13)  = ISODEV(5)-COS2*ISODEV(9)
C 
C     CALCUL DE A
C 
      K(14)    = ISODEV(4)+COS2*ISODEV(8)
      K(15)    = ISODEV(4)-COS2*ISODEV(8)
      K(16)    = -SIN2*R2*ISODEV(8)
C 
      K(17)     = ISODEV(6)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     !!!!! ATTENTION !!!! necessite de mettre les contraintes a zero
C 
C     Cette routine calcule la contribution de la nouvelle deformation NUEPS
C     point de Gauss par point de Gauss a la nouvelle contrainte SIGNOU.
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions du temps de l'iteration
C     E ...... TOUEPS   deformation rangee (NEPS*NTETA*NGAU1*NFONCT)
C     E ...... TOUSAU   les sauts ranges (NSAU*NTETA*NGAU1*NFONCT)
C 
C     Et on recupere :
C 
C     S ...... SIGNOU   nouvelles contraintes rangees
C                       (NEPS, NTETA, NGAU1, NFONCT)
C     S ...... SGNNOU   nouvelles contraintes normales rangees
C                       (NSAU, NTETA, NGAU2, NFONCT)
C 
      SUBROUTINE SIGISO (NFONCT, TOUEPS, TOUSAU, SIGNOU, SGNNOU)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NFONCT
C 
      DOUBLE PRECISION TOUEPS(NEPS*NTETA*NGAU1*NFONCT)
      DOUBLE PRECISION SIGNOU(NEPS*NTETA*NGAU1*NFONCT)
      DOUBLE PRECISION TOUSAU(NSAU*NTETA*NGAU2*NFONCT)
      DOUBLE PRECISION SGNNOU(NSAU*NTETA*NGAU2*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  NUCOU, NUCOL, TETA, X, Y
      INTEGER  NUINT, NUSAU
      INTEGER  DEPS, NUEPS, DEBSIG, LONEC
      INTEGER  DSAU, DEBSGN
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      DOUBLE PRECISION KLOC(17), COUISO(10)
      DOUBLE PRECISION KLOCI(9), INISO(3)
C 
      INTEGER          ADTETA
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SIGISO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
      DEPS = 1
C 
C     BOUCLE SUR LES FONCTIONS DU TEMPS
C 
      DO NUEPS = 1 , NFONCT
C 
        DEBSIG = 1
C 
C **********************************************************************
C *
C *     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C *
C **********************************************************************
C 
C       BOUCLE i SUR LES COUCHES
C 
        DO NUCOU = 1, NBCOU
C 
          CALL COCISO (NUCOU, COUISO)
          CALL RICISO (COUISO, KLOC)
C 
C         POUR TEST
C 
C         BOUCLE ii SUR LES COLONNES
C 
          DO NUCOL = 1, NBCOL
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            DO X = 1, XINTEG
C 
C             BOUCLE iv SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT z
C 
              DO Y = 1, YINTEG
C 
C               BOUCLE v SUR LES ANGLES
C 
                DO TETA = 1, NTETA
C 
C                 ON RECUPERE LES QUANTITES STOCKEES DANS LE DIRECTORY Q-CHAPEAU
C                 POUR LE POINT DE GAUSS CONSIDERE POUR TOUTES LES MATRICES :
C 
                  CALL CTSIGC (NFONCT, NUEPS, TOUEPS(DEPS),
     &                         KLOC, DEBSIG, SIGNOU, LONEC)
C 
                  DEPS   = DEPS  +NEPS
                  DEBSIG = DEBSIG+NSIG
C 
C               FIN DE BOUCLE v SUR LES ANGLES
C 
                END DO
C 
C             FIN DE BOUCLE iv SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT z
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES COLONNES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES COUCHES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES FONCTIONS DE TEMPS
C 
      END DO
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
      DSAU   = 1
C 
C     BOUCLE SUR LES FONCTIONS DE TEMPS
C 
      DO NUSAU = 1 , NFONCT
C 
        DEBSGN = 1
C 
C       BOUCLE i SUR LES INTERFACES
C 
        DO NUINT = 1 , NBINT
C 
          CALL COIISO (NUINT, INISO)
          CALL INTISO (INISO,  KLOCI)
C 
C         BOUCLE ii SUR LES COLONNES
C  
          DO NUCOL = 1, NBCOL
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            DO X = 1, XINTEG
C 
C             BOUCLE iv SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
C 
C               ON RECUPERE LES QUANTITES STOCKEES DANS LE DIRECTORY Q-CHAPEAU
C               POUR LE POINT DE GAUSS CONSIDERE POUR TOUS LES TEMPS
C 
                CALL CTSIGI (NFONCT, NUSAU, TOUSAU(DSAU),
     &                       KLOCI, DEBSGN, SGNNOU, LONEC)
                DSAU    = DSAU   + NSAU
                DEBSGN  = DEBSGN + NSAU
C 
C             FIN DE BOUCLE iv SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES COLONNES
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES INTERFACES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES FONCTIONS DE TEMPS
C 
      END DO
C 
      IF (NBINT .GE. 1) THEN
        CALL TESTEN (DEBSGN-1, NSAU*NTETA*NGAU2, IDPROG)
        CALL TESTEN (DSAU-1, NSAU*NTETA*NGAU2*NFONCT, IDPROG)
        CALL TESTEN (LONEC, NSAU*NTETA*NGAU2*NFONCT, IDPROG)
      END IF
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
C     E ...... ISODEV   Caracteristiques du comportement de la
C                       couche dans sa base d'orthotropie
C 
C           |K11  K12   0   0   0  A1|  E11
C           |     K22   0   0   0  A2|  E22
C           |          K66  0   0   0|  R2*E12
C           |               B1  0   0|  R2*E23  =====>DEFORMATIONS
C           |                   B2  0|  R2*E13
C           |   SYM                 C|  E33
C 
C     E ...... TETORT   Valeur de l'angle que fait cette base avec
C                       le repere global
C     E ...... TETCAL   Valeur de l'angle pour lequel on calcule
C                       le comportement
C 
C     Et on recupere :
C 
C     S ...... K(17)    comporetment isotrope transverse
C     S ...... K(3,3)   Partie du comportement travaillant sur
C                       la partie plane des deformations
C     S ...... B(2,2)   Partie du comportement travaillant sur
C                       le cisaillement normal
C     S ...... A(3)     Partie du comportement reliant la partie plane
C                       des deformations a la deformation normale
C     S ...... C        Relie la contrainte a la deformation normale
C 
      SUBROUTINE RICISO (ISODEV, K)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  ISODEV(10) , K(17)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD     LOGICAL LTRACP
C 
      CHARACTER*6 IDPROG
       PARAMETER (IDPROG='RICISO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     CALCUL DE K
C 
      K(1)  = ISODEV(1)
      K(2)  = ISODEV(2)
      K(3)  = 0.D0
      K(4)  = K(2)
      K(5)  = ISODEV(1)
      K(6)  = 0.D0
      K(7)  = K(3)
      K(8)  = K(6)
      K(9)  = ISODEV(3)
C 
C     CALCUL DE B
C 
      K(10)  = ISODEV(5)
      K(11)  = 0.D0
      K(12)  = 0.D0
      K(13)  = ISODEV(5)
C 
C     CALCUL DE A
C 
      K(14)    = ISODEV(4)
      K(15)    = ISODEV(4)
      K(16)    = 0.D0
C 
      K(17)     = ISODEV(6)
C 
CD    IF(LTRACP(1))THEN
CD      CALL IMPTDP(' K ',K(1),3,3 )
CD      CALL IMPTDP(' B ',K(10),2,2 )
CD       CALL OMPTDP( 'A', K(14),3,1)
CD       CALL IMPDP('C',K(17))
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
C     E ...... ISODEV   Caracteristiques du comportement de
C                       l'interface dans sa base d'orthotropie
C     E ...... TETORT   Valeur de l'angle que fait cette base avec
C                       le repere global
C     E ...... TETCAL   Valeur de l'angle pour lequel on calcule
C                       le comportement
C 
C     Et on recupere :
C 
C     S ...... K(9)     comportement travaillant sur les sauts
C                       range k11, k22, k12, k33

      SUBROUTINE INTELA (ISODEV, TETORT, TETCAL, K)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  ISODEV(3) , TETORT , TETCAL , K(9)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  COS2 , SIN2 , TETLOC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='INTELA')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      TETLOC = TETCAL-TETORT
C 
      COS2   = DCOS(2.D0*TETLOC)
      SIN2   = DSIN(2.D0*TETLOC)
C 
C     CALCUL DE K
C 
      K(1)    =  ISODEV(1)+ISODEV(3)*COS2
      K(2)    =  -ISODEV(3)*SIN2
      K(3)    =  0.D0
      K(4)    =  K(2)
      K(5)    =  ISODEV(1)-ISODEV(3)*COS2
      K(6)    =  0.D0
      K(7)    =  0.D0
      K(8)    =  0.D0
      K(9)    =    ISODEV(2)
C 
CD    CALL IMPDP (
CD     'VALEUR DE L''ORIENTATION PAR RAPPORT AU REPERE GLOBAL',TETCAL)
CD    CALL IMPDP (
CD     'VALEUR DE L''ORIENTATION DE LA BASE D''ORTHOTROPIE   ',TETORT)
CD    CALL OMPTDP('VALEUR DU COMPORTEMENT D''INTERFACE',
CD     K(1) , 9 , 1)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
         SUBROUTINE INTISO( ISODEV ,  K  )
C On envoie comme arguments:
C E................ ISODEV   Caracteristiques du comportement de
C                            l''interface dans sa base d'orthotropie
C Et on recupere:
C S................ K(9)     comportement isotropetravaillant sur les sauts
C                            range , k11,k22,k12,k33
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
C
      include 'cominc.h'
C
      DOUBLE PRECISION  ISODEV(3) ,  K(9)
C**********************************************************************
C
C
      CHARACTER*6 IDPROG
       PARAMETER (IDPROG='INTISO')
C
C
CD    CALL WLKBCD(IDPROG)
C
C***********************************************************************
C CALCUL DE K
      K(1)    =  ISODEV(1)
      K(2)    =  0.D0
      K(3)    =  0.D0
      K(4)    =  0.D0
      K(5)    =  ISODEV(1)
      K(6)    =  0.D0
      K(7)    =  0.D0
      K(8)    =  0.D0
      K(9)    =  ISODEV(2)
C -
CD    CALL OMPTDP('VALEUR DU COMPORTEMENT D''INTERFACE',
CD     K(1) , 9 , 1)
C -
C***********************************************************************
CD    CALL RETOUD(IDPROG)
      RETURN
      END

C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     ADCOU est l'adresse de depart du tableau :
C 
C     HOO-COUCHE => NORME AU SENS DE K0
C     SOU-COUCHE => NORME AU SENS DE K0-1
C 
C     ADINT est l'adresse de depart du tableau :
C 
C     HOO-INTERF => NORME AU SENS DE K0
C     SOU-INTERF => NORME AU SENS DE K0-1
C 
C     Cette routine calcule l'integrale en espace des NFONCT fonctions de
C     l'espace WK par les NFONCT fonctions WL pour les matrices K0.
C     Le resultat de cette routine est le tableau TLKJI a 4 indices
C     tel que, en rangement informatique :
C 
C     I = 1 , NFONCT
C       J = 1 , I
C         K = 1 , NK
C           L = 1 , K
C                                 _
C             TLKJI ( l , k , j , i ) = int ( TR( Wk Kij WL )  )
C 
C 
C     La matrice KT etant symetrique on ne fait varier J que de
C     1 a I ce qui correspond au rangement de KIJPG
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions
C     E ...... LONK     la longueur de l'enregistrement de KIJPG
C     E                 pour le point de Gauss I, J fixes
C     E ...... EPSWK    le r tableau de deformations range
C     E                 (NEPS*NTETA*NGAU1*NK)
C     E ...... SAUTK    le 1er tableau de saut range
C     E                 (NSAU*NTETA*NGAU2*NK)
C 
C     Et on recupere :
C 
C     S ...... TLKJI    le tableau des traces au point de gauss
C                       ranges comme indique plus haut
C 
      SUBROUTINE NOKGLO (ADCOU, ADINT, NFONCT, NK, EPSWK, SAUTK, TLKJI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           ADCOU, ADINT
      INTEGER           NFONCT, NK
      DOUBLE PRECISION  EPSWK(NEPS*NTETA*NGAU1*NK)
      DOUBLE PRECISION  SAUTK(NSAU*NTETA*NGAU2*NK)
      DOUBLE PRECISION  TLKJI(NK, NK, NFONCT, NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION   A, B
      DOUBLE PRECISION   RAYONC, MULTI, MULTJ
C 
C     La contribution de chaque point de gauss a TLKJI
C 
      INTEGER     DEBPG, DECALP, NUCOU
      INTEGER     TETA, NUCOL, EPSI
      INTEGER     NUINT
      INTEGER     AM2LC, ADM2LC, HINT, KINT
      INTEGER     LONTET, LONEPS, DEBWK, WK
      INTEGER     LONTLK, DBLCLK, LCLKJI
      INTEGER     X, Y, I
      INTEGER     ADPOPG, DBPOPG
C 
      DOUBLE PRECISION KLOC(17), COUISO(10)
      DOUBLE PRECISION KLOCI(9), INTISO(3)
C 
      DOUBLE PRECISION TETORT, TETCAL
      INTEGER          ADTETA
C 
CD    LOGICAL            LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NOKGLO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
C 
C     Pour aller lire les points de Gauss
C 
      HINT = XINTEG*(XINTEG-1)/2
      KINT = YINTEG*(YINTEG-1)/2
C  
C     longueur du stockage pour tout un angle de la valeur des NK deformations
C 
      LONEPS = NK*NEPS
C 
      CALL POUSMD (LONEPS, DEBWK)
C 
C     longueur totale du tableau TLKJI
C 
      LONTLK = NFONCT*NFONCT*NK*NK
C 
C     longueur du stockage pour tous les angles d'une fonction
C     representant la contribution de chaque bande a TLKJI
C 
      LONTET = NTETA*LONTLK
C 
      CALL GSPOUD (LONTET, DBLCLK)
C 
C     decalage pour passer d'une fonction EPS a la suivante
C 
      DECALP = NEPS*NTETA*NGAU1
C 
C     debut pour les fonctions EPS des valeurs de epsilon au point de Gauss
C 
      DEBPG = 1
C 
C     Pour recuperer les poids au point de gauss
C 
      CALL ADTBDM ('POIDS-PGCR',  ADPOPG)
      DBPOPG   = ADPOPG
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES COUCHES
C 
      DO NUCOU = 1, NBCOU
C 
        CALL VALEP (NUCOU, B)
        CALL ANGCOU (NUCOU, TETORT)
        CALL COMISO (ADCOU, NUCOU, COUISO)
C 
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1, NBCOL
C 
C          CALL VALRAY( NUCOL , RAYONC , A )
C          JLOC    = A*B
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS LA COUCHE SUIVANT r
C 
          DO X = 1, XINTEG
C 
C            RLOC   = RAYONC + A*GAUSS( HINT + X )
C            MULTI  = POIDS(HINT+X)*JLOC*RLOC
C 
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS LA COUCHE SUIVANT y
C 
            DO Y = 1, YINTEG
C 
              MULTJ = DM(DBPOPG)
              DBPOPG = DBPOPG+1
C 
C             MULTJ  = POIDS(KINT+Y)*MULTI
C 
C             BOUCLE iv SUR LES ANGLES
C 
              LCLKJI  = DBLCLK
C 
              DO TETA = 1, NTETA
C 
C               Recherche de l'angle de la bande correspondant a teta
C 
                TETCAL = DM( ADTETA+TETA-1)
C 
                CALL RICLOC( COUISO , TETORT , TETCAL , KLOC  )
C 
C               Definition des debuts de tableaux des deformations de la 1ere
C               fonction
C 
                EPSI  = DEBPG
C 
C               Debut de boucle sur les fonctions
C 
                WK = DEBWK
C 
                DO I = 1 , NK
C 
                  CALL COPITD (6, EPSWK(EPSI), DM(WK))
                  WK = WK+NEPS
C 
C                 Lecture des valeurs au meme point de GAUSS des
C                 deformations suivantes dans EPSWK
C 
                  EPSI = EPSI+DECALP
C 
                END DO
C 
C               Passage aux deformations du point de Gauss (TETA compris) suivant
C  
                DEBPG  = DEBPG+NEPS
C 
C               Calcul de la contribution du point de gauss a TLKJI
C 
                CALL NOBFPG (NFONCT, MULTJ, KLOC, NK, DM(DEBWK),
     &                       NK, DM(DEBWK), DM(LCLKJI), DM(LCLKJI))
C 
                LCLKJI  = LCLKJI+LONTLK
C 
C             FIN DE BOUCLE iv SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS LA COUCHE SUIVANT y
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS LA COUCHE SUIVANT r
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
CD    CALL IMPMP( ' ******SEQUENCE POUR LES INTERFACES*****')
C 
C     Decalage pour passer d'une fonction SAUTI a la suivante
C 
      DECALP = NSAU*NTETA*NGAU2
C 
C     Debut pour la 1ere fonction SAUTI des valeurs du saut au
C     point de Gauss
C 
      DEBPG = 1
C 
C     Longueur du stockage pour tout un angle de la valeur des
C     NK deformations
C 
      LONEPS = NK*NSAU
C 
      CALL POUSMD( LONEPS , DEBWK )

      IF ( NBINT . GT . 0 ) THEN
C 
C     BOUCLE SUR LES INTERFACES
C 
      DO NUINT = 1 , NBINT
C 
        CALL IOMISO( ADINT , NUINT , INTISO )
        CALL ANGINT( NUINT , TETORT )
C 
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1 , NBCOL
C 
          CALL VALRAY( NUCOL , RAYONC , A )
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'INTERFACE EN r
C 
          DO X = 1 ,XINTEG
C 
C           RLOC  = RAYONC + A*GAUSS( HINT + X )
C           MULTI  = POIDS(HINT+X)*A*RLOC
C 
            MULTI = DM(DBPOPG)
            DBPOPG = DBPOPG+1
C 
C           BOUCLE iii SUR LES ANGLES
C 
            LCLKJI  = DBLCLK
C 
            DO TETA = 1 , NTETA
C 
C           Recherche de l' angle de la bande (tetcal)
C           correspondant a teta
C 
              TETCAL = DM( ADTETA+TETA-1)
C 
              CALL INTELA( INTISO , TETORT , TETCAL , KLOCI  )
C 
C             Definition des debuts de tableaux des sauts de la 1ere
C             fonction
C 
              EPSI  = DEBPG
C 
C             Debut de boucle sur les  fonctions
C 
              WK = DEBWK
C 
              DO I = 1 , NK
C 
                CALL COPITD (3, SAUTK(EPSI), DM(WK))
                WK = WK+NSAU

C 
C               Lecture des valeurs au meme point de GAUSS des
C               sauts suivants dans SAUTK et SAUTL
C  
                EPSI = EPSI+DECALP
C 
              END DO
C 
C             passage aux sauts du point de Gauss ( TETA compris )
C             suivant
C 
              DEBPG  = DEBPG+NSAU
C 
C             Calcul de la contribution du point de gauss a TLKJI
C 
              CALL NIBFPG (NFONCT, MULTI, KLOCI, NK, DM(DEBWK),
     &                     NK, DM(DEBWK), DM(LCLKJI), DM(LCLKJI))
C 
               LCLKJI  = LCLKJI+LONTLK
C 
C           FIN DE BOUCLE iii SUR LES ANGLES
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'INTERFACE EN r
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
      END IF
C 
C     sequence d'integration en teta et calcul de la contribution
C     de tout les angles au tableau  => remplissage de TLKJI
C 
      CALL ITLKJI (NFONCT, NK, NK, DM(DBLCLK), TLKJI(1,1,1,1))
C 
C     exploitation des symetries
C 
      CALL SYMENO (NFONCT, NK, TLKJI(1,1,1,1))
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
C 
CD      CALL OMPTDN( ' TLKJI SOUS LA FORME (NK*NK , N2 )',
CD                     TLKJI(1,1,1,1) , NK*NK ,NFONCT*NFONCT )
C 
CD    END IF
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
C     ADCOU EST L'ADRESSE DE DEPART DU TABLEAU :
C 
C       'HOO-COUCHE' => PRODUIT SCALAIRE AU SENS DE K0
C       'SOU-COUCHE' => PRODUIT SCALAIRE AU SENS DE K0-1
C 
C     ADINT EST L'ADRESSE DE DEPART DU TABLEAU :
C 
C       'HOO-INTERF' => PRODUIT SCALAIRE AU SENS DE K0
C       'SOU-INTERF' => PRODUIT SCALAIRE AU SENS DE K0-1
C 
C 
C     Cette routine calcule l'integrale en espace des NK
C     fonctions de l'espace WK par les NL fonction VL pour les
C     matrices Kij stokees dans MATGLO.
C     Le resultat de cette routine est le tableau TLKJI a 4 indices
C     tel que, en rangement informatique :
C 
C     I = 1 , NFONCT
C       J = 1 , I
C         K = 1 , NK
C           L = 1 , NL
C                                 
C             TLKJI ( l , k , j , i ) = int ( TR( Wk Kij Vl )  )
C 
C 
C     La matrice KT etant symetrique, on ne fait varier J que de
C     1 a I ce qui correpond au rangement de KIJPG.
C 
C     On envoie comme arguments :
C 
C      E ...... NFONCT   le nombre de fonctions du temps
C      E ...... NK       le nombre de vecteurs WK
C      E ...... LONK     la longueur de l'enregistrement de KIJPG
C                        pour le point de Gauss I, J fixes
C      E ...... EPSWK    le 1er tableau de deformations range
C                        (NEPS*NTETA*NGAU1*NFONCT)
C      E ...... SAUTK    le 1er tableau de sauts range
C                        (NSAU*NTETA*NGAU2*NFONCT)
C      E ...... NL       le nombre de vecteur VL
C      E ...... EPSVL    le 2eme tableau de deformations range
C                        (NEPS*NTETA*NGAU1*NFONCT)
C      E ...... SAUTL    le 2eme tableau de saut range
C                        (NSAU*NTETA*NGAU2*NFONCT)
C      Et on recupere :
C 
C      S ...... TLKJI    le tableau des traces au point de Gauss
C                        range comme indique plus haut

      SUBROUTINE SCKGLO (ADCOU, ADINT, NFONCT, NK, EPSWK,SAUTK, NL,
     &                   EPSVL, SAUTL, TLKJI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           ADCOU, ADINT
      INTEGER           NFONCT, NK, NL
C 
      DOUBLE PRECISION  EPSWK(NEPS*NTETA*NGAU1*NK)
      DOUBLE PRECISION  EPSVL(NEPS*NTETA*NGAU1*NL)
      DOUBLE PRECISION  SAUTK(NSAU*NTETA*NGAU2*NK)
      DOUBLE PRECISION  SAUTL(NSAU*NTETA*NGAU2*NL)
      DOUBLE PRECISION  TLKJI(NL,NK,NFONCT,NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION   B
      DOUBLE PRECISION   MULTI, MULTJ
C 
C     La contribution de chaque point de gauss a TLKJI
C 
      INTEGER     DEBPGK, DEBPGL, DECALP, NUCOU
      INTEGER     TETA, NUCOL,  EPSIK, EPSIL
      INTEGER     NUINT
      INTEGER     AM2LC, ADM2LC, HINT, KINT
      INTEGER     LONLIF
      INTEGER     LONTET, LONEPK, LONEPL, DEBWK, WK, DEBVL, VL
      INTEGER     LONTLK, DBLCLK, LCLKJI
      INTEGER     X, Y, I
      INTEGER     ADPOPG, DBPOPG
C 
      DOUBLE PRECISION KLOC(17), COUISO(10)
      DOUBLE PRECISION KLOCI(9), INTISO(3)
      DOUBLE PRECISION TETORT, TETCAL
      INTEGER          ADTETA
C 
CD    LOGICAL            LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SCKGLO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('ANGLES-GEO', ADTETA )
C 
C     Pour aller lire les points de Gauss
C 
      HINT = XINTEG*(XINTEG-1)/2
      KINT = YINTEG*(YINTEG-1)/2
C 
C     Longueur du stockage pour tout un angle de la valeur des
C     NK deformations EPSWK
C 
      LONEPK = NK*NEPS
C 
C     Longueur du stockage pour tout un angle de la valeur des
C     NL deformations EPSWL
C 
      LONEPL = NL*NEPS
C 
      CALL POUSMD( LONEPK+LONEPL , DEBWK )
      DEBVL = DEBWK+LONEPK
C 
C     Longueur totale du tableau TLKJI
C 
      LONTLK = NFONCT*NFONCT*NK*NL
C 
C     Longueur du stockage pour tous les angles d'une fonction
C     representant la contribution de chaque bande a TLKJI
C 
      LONTET = NTETA*LONTLK
C 
C     Longueur du stockage pour tous les angles d'une fonction
C     representant la contribution de chaque bande a TLKJI
C     ==> mise a zero  <=> DBLCLK
C 
      CALL GSPOUD( LONTET , DBLCLK )
C 
C     Decalage pour passer d'une fonction EPSK ou EPSL a la suivante
C 
      DECALP = NEPS*NTETA*NGAU1
C 
C     Debut pour les fonction EPSK des valeurs de epsilon au point de Gauss
C 
      DEBPGK = 1
C 
C     Debut pour les fonctions EPSL des valeurs de epsilon au point de Gauss
C 
      DEBPGL = 1
C 
      CALL ADTBDM ('POIDS-PGCR', ADPOPG)
      DBPOPG   = ADPOPG
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES COUCHES
C 
      DO NUCOU = 1, NBCOU
C 
        CALL VALEP (NUCOU, B)
C 
        CALL ANGCOU (NUCOU, TETORT)
C 
        CALL COMISO (ADCOU, NUCOU, COUISO)
C 
C       Pour test
C 
C       COUISO(7) = 0.D0
C       COUISO(8) = 0.D0
C       COUISO(9) = 0.D0
C       COUISO(10) = 0.D0
C 
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1, NBCOL
C 
C         CALL VALRAY (NUCOL, RAYONC, A)
C         JLOC    = A*B
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT r
C 
          DO X = 1 ,XINTEG
C 
C           RLOC   = RAYONC + A*GAUSS( HINT + X )
C           MULTI  = POIDS(HINT+X)*JLOC*RLOC
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT z
C 
            DO Y = 1 ,YINTEG
C 
              MULTJ = DM(DBPOPG)
              DBPOPG = DBPOPG+1
C 
C             MULTJ  = POIDS(KINT+Y)*MULTI
C 
C             BOUCLE iv SUR LES ANGLES
C 
              LCLKJI  = DBLCLK
C 
              DO TETA = 1 , NTETA
C 
C               RECHERCHE DE L'ANGLE DE LA BANDE (tetcal) CORRESPONDANT A teta
C 
                TETCAL = DM( ADTETA+TETA-1)
C 
                CALL RICLOC( COUISO , TETORT , TETCAL , KLOC  )
C 
C               DEFINITION DES DEBUTS DE TABLEAUX DE LA PREMIERE FONCTION
C 
                EPSIK  = DEBPGK
                EPSIL  = DEBPGL
C 
C               BOUCLE v SUR LES FONCTIONS
C 
                WK = DEBWK
                VL = DEBVL
C 
                DO I = 1 , NK
C 
                  CALL COPITD (6, EPSWK(EPSIK), DM(WK))
                  WK = WK+NEPS
C 
C                 LECTURE DES VALEURS AU MEME POINT DE GAUSS DES DEFORMATIONS
C                 SUIVANTES DANS EPSWK
C 
                  EPSIK = EPSIK+DECALP
C 
                END DO
C 
                DO I = 1 , NL
C 
                  CALL COPITD (6, EPSVL(EPSIL), DM(VL))
                  VL = VL+NEPS
C 
C                 LECTURE DES VALEURS AU MEME POINT DE GAUSS DES DEFORMATIONS
C                 SUIVANTES DANS EPSVL
C 
                  EPSIL = EPSIL+DECALP
C 
C               FIN DE BOUCLE v SUR LES FONCTIONS
C 
                END DO
C 
C               PASSAGE AUX DEFORMATIONS DU POINT DE GAUSS SUIVANT (TETA compris)
C 
                DEBPGK  = DEBPGK+NEPS
                DEBPGL  = DEBPGL+NEPS
C 
C               CALCUL DE LA CONTRIBUTION DU POINT DE GAUSS A TLKJI
C 
                CALL SCBFPG (NFONCT, MULTJ, KLOC, NK, DM(DEBWK),
     &	                     NL, DM(DEBVL), DM(LCLKJI), DM(LCLKJI))
C 
                LCLKJI  = LCLKJI+LONTLK
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
      IF ( NBINT . GT .0 ) THEN
C 
CD    CALL IMPMN( ' ******SEQUENCE POUR LES INTERFACES*****')
C 
C     DECALAGE POUR PASSER D'UNE FONCTION SAUTI A LA SUIVANTE
C 
      DECALP = NSAU*NTETA*NGAU2
C 
C     DEBUT POUR LA 1ere FONCTION SAUTK DES VALEURS DU SAUT AU POINT DE GAUSS
C 
      DEBPGK = 1
C 
C     DEBUT POUR LA 1ere FONCTION SAUTL DES VALEURS DU SAUT AU POINT DE GAUSS
C 
      DEBPGL = 1
C 
      LONLIF   = 9*NFONCT*(NFONCT+1)/2
C 
C     BOUCLE SUR LES INTERFACES
C 
      DO NUINT= 1, NBINT
C 
CD      CALL IMPET ('INTERFACE ', NUINT)
	CALL IOMISO (ADINT, NUINT, INTISO)
        CALL ANGINT (NUINT, TETORT)
C 
C       BOUCLE i SUR LES COLONNES
C 
       DO NUCOL = 1, NBCOL
C 
CD        CALL IMPET ('COLONNE ', NUCOL)
C 
C       CALL VALRAY (NUCOL, RAYONC, A)
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMET SUIVANT r
C 
          DO X = 1, XINTEG
C 
CD          CALL IMPET ('XGAUSS ', X)
C 
C           RLOC  = RAYONC + A*GAUSS( HINT + X )
C           MULTI  = POIDS(HINT+X)*A*RLOC
C 
            MULTI = DM(DBPOPG)
            DBPOPG = DBPOPG+1
C 
C           BOUCLE iii SUR LES ANGLES
C 
            LCLKJI  = DBLCLK
C 
            DO TETA = 1, NTETA
C 
CD            CALL IMPET ('ANGLE ', TETA)
C 
C             RECHERCHE DE L'ANGLE DE LA BANDE (tetcal) CORRESPONDANT A teta
C 
              TETCAL = DM(ADTETA+TETA-1)
C 
              CALL INTELA (INTISO, TETORT, TETCAL, KLOCI)
C 
C             DEFINITION DES DEBUTS DE TABLEAUX DES SAUTS DE LA 1ere FONCTION
C 
              EPSIK  = DEBPGK
              EPSIL  = DEBPGL
C 
C             BOUCLE iv SUR LES FONCTIONS
C 
              WK = DEBWK
              VL = DEBVL
C 
              DO I = 1, NK
C 
C 
                CALL COPITD (3, SAUTK(EPSIK), DM(WK))
                WK = WK+NSAU
C 
C               LECTURE DES VALEURS AU MEME POINT DE GAUSS DES
C               SAUTS SUIVANTS DANS SAUTK
C 
                EPSIK = EPSIK+DECALP
C 
              END DO
C 
              DO I = 1 , NL
C 
                CALL COPITD( 3  , SAUTL(EPSIL) , DM(VL) )
                VL = VL+NSAU
C 
C               LECTURE DES VALEURS AU MEME POINT DE GAUSS DES
C               SAUTS SUIVANTS DANS SAUTL
C 
                EPSIL = EPSIL+DECALP
C 
              END DO
C 
C             PASSAGE AUX SAUTS DU POINT DE GAUSS (TETA compris) suivant
C 
               DEBPGK  = DEBPGK+NSAU
               DEBPGL  = DEBPGL+NSAU
C 
C             CALCUL DE LA CONTRIBUTION DU POINT DE GAUSS A TLKJI
C 
              CALL SIBFPG( NFONCT , MULTI , KLOCI , NK ,
     &                     DM(DEBWK), NL , DM( DEBVL) ,  DM(LCLKJI) ,
     &                     DM(LCLKJI) )
C 
               LCLKJI  = LCLKJI+LONTLK
C 
C           FIN DE BOUCLE iii SUR LES ANGLES
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMET SUIVANT r
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
      END IF
C 
C     sequence d'integration en teta et calcul de la contribution
C     de tout les angles au tableau  => remplissage de TLKJI
C 
      CALL ITLKJI (NFONCT, NK, NL, DM(DBLCLK), TLKJI)
C 
      CALL SYMESC (NFONCT, NK, NL, TLKJI(1,1,1,1))
C 
CD    IF (LTRACP(1))THEN
C 
CD     CALL OMPTDP( 'TLKJI apres symetrie diagonale  ( NK*NL  , N2 )',
CD                   TLKJI(1,1,1,1)  , NK*NL , NFONCT*NFONCT )
C 
CD    END IF
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
C 
CD      CALL OMPTDN( ' TLKJI SOUS LA FORME N2,(NK*NL) ',
CD                     TLKJI(1,1,1,1) , NL*NK ,NFONCT*NFONCT )
C 
CD    END IF
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine a pour but de determiner la contrainte a zero
C     elements finis associee a la solution en deplacement
C 
C     On envoie comme arguments :
C 
C     E ...... DEPSOL La solution developpee fourier en deplacement
C                     de ce probleme
C     E ...... EPSSOL La solution en deformation de ce probleme
C     E ...... SAUSOL La solution en saut de ce probleme
C 
C     Et on recupere :
C 
C     S ...... SIGZER les contraintes rigoureusement admissibles
C     S               a zero elements finis associees en sortie
C     S ...... SGNZER les contraintes normales rigoureusement
C     S               admissibles a zero associees en sortie
C 
C     RAPPEL
C 
C     Intemps(fotemp*K0*fotemp)[DEPSOL] = Intemps(fotemp*-*DELTA(sigchap))

      SUBROUTINE CSIAD0 (DEPSOL, EPSSOL, SAUSOL, SIGZER, SGNZER)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION DEPSOL(NDDL*NBMAT), EPSSOL(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION SIGZER(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION SGNZER(NSAU*NTETA*NGAU2)
      DOUBLE PRECISION SAUSOL(NSAU*NTETA*NGAU2)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  AM2LC, ADM2LC
C 
CD    LOGICAL  LTRACN, LTRACP
C 
C     Pour les longueurs provisoires
C 
      INTEGER LONEPS, LONSAU, LONDEP, LONDER, LONRES
      INTEGER SIGITE, SINITE
C 
C     Pour test
C 
CD    DOUBLE PRECISION TDENOM, DENOM, TRACE, TRALOC, RAP
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='CSIAD0')
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     CALCUL DE LA CONTRAINTE ISSUE DES ITERATIONS
C 
      LONEPS = NTETA*NEPS*NGAU1
      LONDER = NDDL*NTETA
      LONDEP = NDDL*NBMAT
      LONSAU = NTETA*NSAU*NGAU2
C 
      LONRES = LONEPS+LONSAU
      CALL GSPOUD (LONRES, SIGITE)
C 
      SINITE = SIGITE + LONEPS
C 
      CALL NOUSIG (1, EPSSOL, SAUSOL, DM(SIGITE), DM(SINITE))
C 
      CALL SOU (DM(SIGITE), SIGZER, SIGZER, LONEPS, 1)
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL SOU (DM(SINITE), SGNZER, SGNZER, LONSAU, 1)
C 
      END IF
C 
C     TEST TRALOC = SUP DES TRACES LOCALES
C 
C     CD    CALL TRAGLO (EPSSOL, SAUSOL, DM(SIGITE), DM(SINITE),
C     CD                 DENOM, TDENOM)
C     CD    CALL TRAGLO (EPSSOL, SAUSOL, SIGZER, SGNZER,
C     CD                 TRACE, TRALOC)
C     CD    RAP = TRACE/DENOM
C     CD    CALL IMPDN ('EPSSOL*SIGZER/EPSOL*K*EPSOL = 0?', RAP)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
