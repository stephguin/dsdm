C     Preliminaire a l'etape locale.
C 
      SUBROUTINE PRETEM
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
      INTEGER  FTECH6, FOTEMP, FTREE7, FTRE10, FTSCH5
      INTEGER  LONRES
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='PRETEM')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C     Ouverture des tableaux des accroissements admissibles qu'on 
C     remplira a l'etape globale. Ils sont stockes (chamax, npicet)
C 
      LONRES = NPICET*CHAMAX
      CALL GESTDP ('TEMPS-SIGM', LONRES, FTSCH5)
      CALL GESTDP ('TEMPS-EPSI', LONRES, FTECH6)
C 
C     Ouverture des tableaux des solutions admissibles.
C     On recopie la premiere correspondant a l'elasticite;
C     ils sont stockes (charax, npicet).
C 
      LONRES = NPICMX*CHARAX
      CALL GESTDP ('TEMPS-REEL', LONRES, FTREE7)
      CALL GESTDP ('TEMP-SI-RE', LONRES, FTRE10)
      CALL ADTBDM ('F-TEMPS-DO', FOTEMP)
      CALL RCHARE (7, NPICET, DM(FOTEMP), DM(FTREE7))
      CALL RCHARE (10, NPICET, DM(FOTEMP), DM(FTRE10))
C 
CD    CALL IMPTDT ('TEMPS-REEL', DM(FTREE7), CHARAX, NPICET)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Regularisation des fonctions du temps
C 
C     On envoie comme arguments :
C 
C     E ...... NPOINT   nombre de points
C     E ...... COEFA    coefficient de viscosite
C     E ...... INTERV   pas de la discretisation
C     E ...... F        fonction du temps a modifier
C 
C     Et on recupere :
C 
C     S ...... Y        nouvelle fonction du temps
C 
      SUBROUTINE WVISC1 (NPOINT, COEFA, INTERV, F, Y)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NPOINT
C 
      DOUBLE PRECISION COEFA
      DOUBLE PRECISION INTERV(NPOINT)
      DOUBLE PRECISION F(NPOINT)
      DOUBLE PRECISION Y(NPOINT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  A0
C 
C     Valeurs intermediaires
C 
      DOUBLE PRECISION COEF1, COEF2
C 
C     Indices
C 
      INTEGER IPOINT, AD, AM2LC, ADM2LC, I
C 
      CHARACTER*6 IDPROG
      PARAMETER  (IDPROG = 'WVISC1')
C 
CD    CALL WLKBCD (IDPROG)
C -----------------------------------------------------------------------
      CALL ENPOUB (AM2LC, ADM2LC)
C 
      CALL POUSMD (NPOINT, AD)
      A0  = F(1)/INTERV(1)
      DM(AD) = A0
C 
      Y(1) = F(1)
C 
      COEF2 = .9D0
      COEF1 = 1.D0 / (1.D0 + COEF2)
C 
      AD = AD -1
C 
      DO IPOINT = 2, NPOINT
C 
        DM(AD +IPOINT) = COEF1*(F(IPOINT)/INTERV(IPOINT)+ COEF2*A0 )
        A0             = DM(AD+IPOINT)
        Y(IPOINT)      = A0*INTERV(IPOINT)
C 
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
C     Cette routine calcule le produit scalaire L2 de deux fonctions
C     du temps en faisant une interpolation lineaire.
C      (1/duree) INT(f1*f2)
C 
C     On envoie comme arguments :
C 
C     E ...... ADINTE  adresse du tableau des intervalles en temps
C     E ...... ADPICE  adresse du tableau des piquets de temps
C     ES...... FOTEM1  valeur de la 1ere fonction du temps pour les npicet
C     ES...... FOTEM2  valeur de la 2eme fonction du temps pour les npicet
C 
C     Et on recupere :
C 
C     S ...... PSCALT  valeur du produit scalaire
C 
      SUBROUTINE SCALFT (ADINTE, ADPICE, FOTEM1, FOTEM2, PSCALT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER ADINTE, ADPICE
C 
      DOUBLE PRECISION FOTEM1(NPICET), FOTEM2(NPICET), PSCALT
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TEMPS, DBVLT, I
CD    LOGICAL  LTRACN, LTRACP
      DOUBLE PRECISION  VALPR1, MULT, FI
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SCALFT')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF ( LTRACP(1) ) THEN
CD      CALL OMPTDP ('1ERE FONCTION DU TEMPS ', FOTEM1(1), NPICET, 1)
CD      CALL OMPTDP ('2EME FONCTION DU TEMPS ', FOTEM2(1), NPICET, 1)
CD    END IF
C 
      PSCALT = 0.D0
      VALPR1 = FOTEM1(1)*FOTEM2(1)
      DBVLT  = ADINTE-1
C 
      DO TEMPS = 1 , NPICET
C 
          MULT   = DM(DBVLT+TEMPS)
C 
          FI     = FOTEM1(TEMPS)*FOTEM2(TEMPS)
C 
          PSCALT = PSCALT+MULT*(FI+ VALPR1)*.5D0
C 
          VALPR1 = FI
C 
      END DO
C 
      PSCALT = PSCALT/LTEMPS
C 
CD    CALL IMPDN ('VALEUR DU PRODUIT SCALAIRE ', PSCALT)
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine norme les fonctions du temps.
C      Norme = [(1/duree) INT(f*f) ]**0.5
C 
C     On envoie comme arguments :
C 
C     E ...... ADINTE  adresse du tableau des intervalles en temps
C     E ...... ADPICE  adresse du tableau des piquets de temps
C     ES...... FOTEMP  valeur de la fonction du temps pour les npicet
C 
C     Et on recupere :
C 
C     S ...... FOTENR  valeurs normees
C     S ...... NORME   valeur de la norme
C 
      SUBROUTINE NORMFT (ADINTE, ADPICE, FOTEMP, FOTENR, NORME)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER   ADINTE, ADPICE
C 
      DOUBLE PRECISION FOTEMP(NPICET), FOTENR(NPICET), NORME
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NORMFT')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      CALL SCALFT (ADINTE, ADPICE, FOTEMP, FOTEMP, NORME)
      NORME = DSQRT(NORME)
      DO I = 1, NPICET
        FOTENR(I) = FOTEMP(I)/NORME
      END DO
C 
CD    IF (LTRACP(1)) THEN
CD      CALL IMPTDT ('FONCTION DU TEMPS NORMEE ', FOTENR(1), NPICET, 1)
CD      CALL IMPDT ('NORME EN TEMPS ', NORME)
CD    END IF
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule les sommes sur le temps de NFONCT quantites
C     type contraintes (ou deformations...) dt * une fonction du temps
C     (NFONCT=6 pour des contraintes) * (1/duree)
C 
C     On suppose que les quantites sont multipliees par l'intervalle de
C     temps t(i)-t(i-1)
C 
C     On envoie comme arguments :
C 
C     E ...... INTERV  tableau des intervalles de temps
C     E ...... NFONCT  le nombre de fonctions a sommer avec GT
C     E ...... SIGMA   valeur des NFONCT quantites pour tous les npicets * DT
C                      stockees (NFONCT, npicet)
C     E ...... GT      valeur de la fonction du temps pour les npicet
C 
C     Et on recupere :
C 
C     S ...... SSIGGT  valeur du produit scalaire des NFONCT quantites par GT
C 
      SUBROUTINE STSIPA (INTERV, NFONCT, SIGMA, GT, SSIGGT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      INTEGER           NFONCT
C 
C     MODIF NPICAC --> NPICET
C 
      DOUBLE PRECISION  SIGMA(NFONCT*NPICET), GT(NPICET)
      DOUBLE PRECISION  INTERV(NPICET)
      DOUBLE PRECISION  SSIGGT(NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   I, VALP1
      INTEGER   DVALP1, TEMPS
      INTEGER   AM2LC, ADM2LC, ADSIG
C 
      DOUBLE PRECISION  MULT, VALAC2,  T1
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='STSIPA')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL POUSMD (NFONCT, VALP1)
      DVALP1    = VALP1
C 
      MULT      = INTERV(1)
      VALAC2    = GT(1)
C 
C     Boucle sur les NFONCT quantites
C 
      ADSIG = 1
C 
      DO I = 1, NFONCT
C 
        T1           = SIGMA(ADSIG)*VALAC2
        ADSIG        = ADSIG+1
        SSIGGT(I)    = MULT*(T1)
        DM(VALP1)    = T1
        VALP1        = VALP1+1
C 
      END DO
C 
C     MODIF      DO TEMPS = 1, NPICAC
C 
      DO TEMPS = 2, NPICET
C 
        MULT      = .5D0*INTERV(TEMPS)
        VALAC2    = GT(TEMPS)
C 
C       Boucle sur les NFONCT quantites
C 
        VALP1  = DVALP1
C 
        DO I = 1, NFONCT
C 
          T1           = SIGMA(ADSIG)*VALAC2
          ADSIG        = ADSIG+1
C  
          SSIGGT(I)    = SSIGGT(I) +MULT*(T1 + DM(VALP1))
C 
          DM(VALP1)    = T1
          VALP1        = VALP1+1
C 
        END DO
C 
      END DO
C 
      DO I = 1, NFONCT
C 
        SSIGGT(I)    = SSIGGT(I)/LTEMPS
C 
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
C     Cette routine calcule la meilleure evolution des deplacements
C     et la meilleure evolution des contraintes au sens de l'etape globale.
C 
C     On envoie comme arguments :
C 
C     E ...... SIGZER  les contraintes developpees admissibles a zero
C     E ...... SGNZER  les contraintes normales developpees
C     E ...... NPICDE  piquet de temps de debut d'integration
C     E ...... NPICFI  piquet de temps de fin d'integration

C     Et on recupere :
C 
C     S ...... FTSIGS  fonction du temps des contraintes en sortie
C 
C     Attention SIGZER et SGNZER sont mis a zero en entree.
C 
      SUBROUTINE EVSIOP (SIGZER, SGNZER, NPICDE, NPICFI, FTSIGS)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NPICDE, NPICFI
C 
      DOUBLE PRECISION  FTSIGS( NPICET )
      DOUBLE PRECISION  SIGZER(NTETA*NEPS*NGAU1)
      DOUBLE PRECISION  SGNZER(NTETA*NSAU*NGAU2)
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
      INTEGER  I
C 
C     Pour l'optimisation en temps
C 
      INTEGER DEKTEM , SCKTEM
      INTEGER DESTEM , SCSTEM
      INTEGER ERRTP1 , ERRTP2
      INTEGER ADCHOO , ADIHOO
C 
      DOUBLE PRECISION DIVIS
C 
      CHARACTER*6 IDPROG
      CHARACTER*3 BARATE
      PARAMETER (IDPROG='EVSIOP')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('SOU-COUCHE', ADCHOO)
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL ADTBDM ('SOU-INTERF', ADIHOO)
C 
      ELSE
C 
        ADIHOO = ADCHOO
C 
      END IF
C 
      CALL GSPOUD (4*NPICET, SCKTEM)
C 
      SCSTEM = SCKTEM + NPICET
      ERRTP1 = SCSTEM + NPICET
      ERRTP2 = ERRTP1 + NPICET
C 
C     Calcul de intesp(sigadm* INVKT * sigadm)
C 
      CALL NOKGLO (ADCHOO, ADIHOO, NFTGLO,
     &             NFTGLO, SIGZER, SGNZER, DIVIS)
C 
C     Calcul de intesp(delta(sig_point_chapeau)*K0-1*sigzer)
C 
C     CALL SCASTE ('epsich', 'sautch', SIGZER, SGNZER, NPICDE, NPICFI,
C     &             DM(SCSTEM), DM(ERRTP2))
C 
C     CALL IDENTI (NBETLC, BARATE)
      CALL SCHAST (2, 'sigmch', 'sinoch', 1, SIGZER,
     &             SGNZER, NPICDE, NPICFI, DM(SCSTEM))
C 
      CALL MUMARE (-1.D0, NPICET, DM(SCSTEM), DM(SCSTEM))
C 
CD    CALL IMPTDN ('-sigzer*delta(epscha)) ',
CD                  DM(SCSTEM), NPICET, 1)
C 
      DEKTEM = SCKTEM
      DESTEM = SCSTEM
C 
      DO I = 1 , NPICET
C 
        FTSIGS(I)  = DM(DESTEM)/DIVIS
C 
        DESTEM     = DESTEM+1
C 
      END DO
C 
CD    CALL IMPTDN ('MEILLEURE EVOLUTION DES CONTRAINTES ',
CD                    FTSIGS(1), NPICET, 1)
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
C     Cette routine calcule la meilleure evolution des deplacements
C     au sens de l'etape globale.
C 
C     On envoie comme arguments :
C 
C     E ...... EPSSOL   la solution en deformation
C     E ...... SAUSOL   la solution en saut
C     E ...... NPICDE   piquet de temps de debut d'integration
C     E ...... NPICFI   piquet de temps de fin d'integration
C 
C     Et on recupere :

C     S ...... FTDEPS   fonction du temps des deplacements en sortie
C 
      SUBROUTINE EVDPOP (EPSSOL, SAUSOL, NPICDE, NPICFI, FTDEPS)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NPICDE, NPICFI
C 
      DOUBLE PRECISION  FTDEPS(NPICET)
      DOUBLE PRECISION  EPSSOL(NTETA*NEPS*NGAU1)
      DOUBLE PRECISION  SAUSOL(NTETA*NSAU*NGAU2)
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
      INTEGER  I
C 
C     Pour l'optimisation en temps
C 
      INTEGER DEKTEM, SCKTEM, LONRES
      INTEGER DESTEM, SCSTEM
      INTEGER ERRTP1, ERRTP2
      INTEGER  ADIHOO, ADCHOO
C 
      DOUBLE PRECISION DIVIS
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EVDPOP')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('HOO-COUCHE', ADCHOO)
      IF (NBINT .GT. 0) THEN
        CALL ADTBDM ('HOO-INTERF', ADIHOO)
      ELSE
        ADIHOO = ADCHOO
      END IF
C 
      LONRES = 4*NPICET
      CALL POUSMD (LONRES, SCKTEM)
      SCSTEM = SCKTEM + NPICET
      ERRTP1 = SCSTEM + NPICET
      ERRTP2 = ERRTP1 + NPICET
C 
C     Calcul de intesp(epssol*K0*epssol)
C 
      CALL NOKGLO (ADCHOO, ADIHOO, NFTGLO, NFTGLO, EPSSOL, SAUSOL,
     &             DIVIS)
C 
C     Calcul de intesp(-delta(sigba)* epssol)
C 
      CALL SCASTE ('sigmch', 'sinoch', EPSSOL, SAUSOL,
     &              NPICDE, NPICFI, DM(SCSTEM), DM(ERRTP2))
C 
C     Calcul de la meilleure evolution des deplacements f1(t)
C 
      DESTEM = SCSTEM
C 
      DO I = 1, NPICET
        FTDEPS(I) = DM(DESTEM)/DIVIS
        DEKTEM = DEKTEM+1
        DESTEM = DESTEM+1
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
C     Cette routine effectue le produit scalaire en espace de -delta(sign)
C     calculee a l'etape locale par la deformation alphai.
C     Le resultat de cette routine est la fonction du temps SCATEP.
C 
C     On envoie comme arguments :
C 
C     E ...... CARESP   Le caractere * 20 nom du fichier du type deformation
C     E ...... CARSAU   Le caractere * 20 nom du fichier du type saut
C     E ...... ALPHAI   adresse des deformations admissibles a zero rangees (NEPS, NTETA, NGAU1)
C     E ...... SAUTI    adressse des sauts admissibles a zero ranges(NSAU, NTETA, NGAU2)
C     E ...... NPICDE   piquet de temps de debut d'integration
C     E ...... NPICFI   piquet de temps de fin d'integration
C 
C     Hors de NPICDE et NPICFI la fonction du temps en sortie doit etre nulle
C 
C     Et on recupere :
C 
C     S ...... SACTEP   la fonction du  temps INTESP(-delta(sig-point-chapeau)*ALPHAI)
C     S ...... SUPESP   la fonction du temps correspondant au sup sur les elements

      SUBROUTINE SCASTE (CAREPS, CARSAU, ALPHAI, SAUTI,
     &                   NPICDE, NPICFI, SACTEP, SUPESP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  SACTEP(NPICET)
      DOUBLE PRECISION  ALPHAI(NEPS*NTETA*NGAU1)
      DOUBLE PRECISION  SAUTI(NSAU*NTETA*NGAU2)
      DOUBLE PRECISION  SUPESP(NPICET)
C 
      INTEGER NPICDE, NPICFI
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION EPSILI(6),TRACE, A, B, JLOC
      DOUBLE PRECISION RLOC, RAYONC, MULTI, MULTJ
C 
      INTEGER     SIGMAC, SIGNAC, DEBPG, DECAlP, NUCOU
      INTEGER     TETA, NUCOL, BMTGAC, EPSI
      INTEGER     DEBUT, TEMPS, NUINT
      INTEGER     AM2LC, ADM2LC, HINT, KINT
      INTEGER     LONLIF
      INTEGER     NBFON, FONTET, TETTOT, LONTOT, ADTET
      INTEGER     AVTET, X, Y
      INTEGER     LONSCH, IUNSIG, NUENSI
C 
CD    LOGICAL            LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      CHARACTER*20  NOMFIC, NOM, CAREPS, CARSAU
      PARAMETER (IDPROG='SCASTE')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACP(1)) THEN
C 
CD      CALL OMPTDP ('ALPHAI ',ALPHAI(1), NEPS*NTETA, NGAU1)
CD      CALL OMPTDP ('SAUTI  ',SAUTI(1), NSAU*NTETA, NGAU2)
C 
CD    END IF
C 
C     Pour aller lire les points de Gauss
C 
      HINT = XINTEG*(XINTEG-1)/2
      KINT = YINTEG*(YINTEG-1)/2
C 
C     Recherche du nom du fichier ou est stocke sigma barre
C 
      NOM = CAREPS
      NOMFIC = NOM
      LONSCH = 6*NPICET
      NUENSI  = 1
      CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C 
C     Nombre de fonctions du temps  calculees
C 
      NBFON  = 1
C 
C     Longueur du stockage pour tous les angles d'une fonction
C 
      FONTET = NTETA*NPICET
C 
C     Longueur totale du tableau partiel en teta
C 
      TETTOT = NBFON*FONTET
C 
C     Longueur totale du tableau des fonctions
C 
      LONTOT = NBFON*NPICET
C 
C     Decalage pour passer d'une fonction ALPHAI a la suivante
C 
      DECALP = NEPS*NTETA*NGAU1
C 
C     Debut pour la 1ere fonction ALPHAI des valeurs de epsilon au point de Gauss
C 
      DEBPG = 1
C 
C     MODIF => BALAID DE SCATEP
C 
      CALL BALAI (SACTEP, LONTOT, 1)
      CALL BALAI (SUPESP, LONTOT, 1)
C 
C     SIGMAC    -2*contrainte adm actuelle      X 6*NPICET
C     ADTET     tableau en teta                 X TETTOT
C     DEBFON    tableau en Yinteg               X LONTOT
C 
      LONLIF   = 6*NPICET
C 
      CALL POUSMD (LONLIF, SIGMAC)
      BMTGAC = SIGMAC
C 
      CALL GSPOUD (TETTOT, ADTET)
      AVTET  = ADTET-1
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
C       BOUCLE i SUR LES COLONNES
C 
        DO NUCOL = 1, NBCOL
C 
          CALL VALRAY (NUCOL, RAYONC, A)
          JLOC = A*B
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
C 
          DO X = 1, XINTEG
C 
            RLOC   = RAYONC + A*GAUSS( HINT + X )
            MULTI  = POIDS(HINT+X)*JLOC*RLOC
C 
C           BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Z
C 
            DO Y = 1, YINTEG
C 
             MULTJ = POIDS(KINT+Y)*MULTI
C 
C             BOUCLE iv SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
C             Lecture directe du fichier (dans q-chapeau) des accroissements
C             admissibles de l'etape locale precedente.
C 
                CALL LFDDNF (DM(BMTGAC), BMTGAC, LONSCH,
     &                       IUNSIG, NUENSI)
C 
                NUENSI = NUENSI +1
C 
C               Definition des debuts de tableaux des deformations
C 
                EPSI  = DEBPG
C 
                CALL COPITD( 6  , ALPHAI(EPSI) , EPSILI(1) )
C 
                DEBUT = AVTET+TETA
C 
                SIGMAC     = BMTGAC
C 
C               DEBUT DE BOUCLE v SUR LES PIQUETS DE TEMPS
C 
                DO  TEMPS = 1, NPICDE-1
C 
                  DEBUT  = DEBUT+NTETA
C 
                  SIGMAC   = SIGMAC+NEPS
C 
                END DO
C 
                DO  TEMPS = NPICDE, NPICFI
C 
                  CALL SCAVEC (NEPS, DM(SIGMAC), EPSILI(1), TRACE)
C 
                  DM(DEBUT) = DM(DEBUT)+ MULTJ*TRACE
                  DEBUT  = DEBUT+NTETA
C 
                  SUPESP(TEMPS) = MAX(SUPESP(TEMPS), TRACE)
C 
                  SIGMAC   = SIGMAC+NEPS
C 
                END DO
C 
                DO  TEMPS = NPICFI+1 , NPICET
C 
                   DEBUT  = DEBUT+NTETA
C 
                   SIGMAC   = SIGMAC+NEPS
C 
C               FIN DE BOUCLE v SUR LES PIQUETS DE TEMPS
C 
                END DO
C 
C               Passage aux deformations du point de Gauss (TETA compris) suivant
C 
                DEBPG  = DEBPG+NEPS
C 
C             FIN DE BOUCLE iv SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Z
C 
            END DO
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
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
C     Fermeture de l'unite ou est stocke SIGMA
C 
      CALL FERFIC (3, IUNSIG, IDPROG)
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
C     TEST SUR LE NOMBRE D'INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
C       Debut des valeurs du saut au point de Gauss
C 
        DEBPG = 1
C 
C       SAUTGAC   -*sign admissible actuel        X 3*NPICET
C       ADTET     tableau en teta                 X TETTOT
C       DEBFON    tableau en xinteg               X LONTOT
C 
C       Ouverture du fichier pour les contraintes normales
C 
        NOM = CARSAU
        NOMFIC = NOM
        LONSCH = 3*NPICET
        NUENSI = 1
        CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT = 1, NBINT
C 
C         BOUCLE i SUR LES COLONNES
C 
          DO NUCOL = 1, NBCOL
C 
            CALL VALRAY (NUCOL, RAYONC, A)
C 
C           BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
C 
            DO X = 1, XINTEG
C 
              RLOC   = RAYONC + A*GAUSS( HINT + X )
              MULTI  = POIDS(HINT+X)*A*RLOC
C 
C             BOUCLE iii SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
C               Lecture directe du fichier (dans q-chapeau) des accroissements
C               admissibles de l'etape locale precedente.
C 
                CALL LFDDNF (DM(BMTGAC), BMTGAC, LONSCH,
     &                       IUNSIG, NUENSI)
C 
                NUENSI = NUENSI +1
C 
C               Definition des debuts de tableaux des sauts de la 1ere fonction
C 
                EPSI  = DEBPG
C 
                CALL COPITD (3, SAUTI(EPSI), EPSILI(1))
C 
                DEBUT = AVTET+TETA
C 
                SIGNAC     = BMTGAC
C 
C               BOUCLE iv SUR LES PIQUERS DE TEMPS
C 
                DO TEMPS = 1, NPICDE-1
C 
                  DEBUT  = DEBUT+NTETA
C 
                  SIGNAC   = SIGNAC+3
C 
                END DO
C 
                DO TEMPS = NPICDE, NPICFI
C 
                  CALL SCAVEC( NSAU, DM(SIGNAC), EPSILI(1), TRACE)
                  DM(DEBUT) = DM(DEBUT)+ MULTI*TRACE
C 
                  SUPESP(TEMPS) = MAX(SUPESP(TEMPS), TRACE)
C 
                  DEBUT  = DEBUT+NTETA
C 
                  SIGNAC = SIGNAC+3
C 
                END DO
C 
                DO TEMPS = NPICFI+1, NPICET
C 
                  DEBUT  = DEBUT+NTETA
C 
                  SIGNAC = SIGNAC+3
C 
C               FIN DE BOUCLE iv SUR LES PIQUETS DE TEMPS
C 
                END DO
C 
C               Passage aux sauts du point de Gauss (TETA compris) suivant
C 
                DEBPG  = DEBPG+NSAU
C 
C             FIN DE BOUCLE iii SUR LES ANGLES
C 
              END DO
C 
C           FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT R
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES COLONNES
C 
          END DO
C 
C       FIN DE BOUCLE SUR LES INTERFACES
C 
        END DO
C 
C 
C       Fermeture de l'unite ou est stocke SIGMA
C 
        CALL FERFIC (3, IUNSIG, IDPROG)
C 
C     FIN DE TEST SUR LE NOMBRE D'INTERFACES
C 
      END IF
C 
C     Sequence d'integration en teta et calcul de la contribution au point de Gauss
C 
      CALL ITBTET (LONTOT, DM(ADTET), SACTEP)
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1)) THEN
C 
CD      CALL OMPTDN ('FONCTION DU TEMPS ', SACTEP, NPICET, 1)
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
C     Modification des fonctions du temps reelles pour tenir compte des
C     orthogonalisations. A priori nterme = npicet.
C 
C     On envoie comme arguments :
C 
C     E ...... NTERME    nombre de composantes du champ
C     E ...... NBCHAR    le nombre de champs ecrits
C     E ...... RCOMPO    les coefficients dus a l'orthogonalisation
C     E ...... NBCHMO    le nombre de champs a modifier
C 
C     Et on recupere :
C 
C     S ...... CHAMP     Le tableau = TBSTOC(NUCHAM, LONCHA)
C 
      SUBROUTINE MOCHAR (NTERME, NBCHAR, RCOMPO, NBCHMO, CHAMP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NTERME, NBCHAR, NBCHMO
C 
      DOUBLE PRECISION CHAMP(CHARAX, NTERME), RCOMPO(NBCHMO)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER TERME, NCHAMP, K
C 
      DOUBLE PRECISION COEFF
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MOCHAR')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C     Initialisation
C 
C     BOUCLE SUR LES FONCTIONS
C 
C     On ne modifie pas la 1ere qui correspond a l'initialisation
C 
      K = 1
C 
C     Addition de la contribution de la fonction orthonormee aux vecteurs precedents.
C     Si ????, on ne rentre pas dans la boucle.
C 
      IF ((NBCHAR-NBCHMO+1) .GT. (NBCHAR-1)) GOTO 1000
C 
      DO NCHAMP = NBCHAR-NBCHMO+1, NBCHAR-1
C 
        COEFF =  RCOMPO(K)
        K     =  K+1
C 
C       BOUCLE i SUR TOUS LES TERMES
C 
        DO TERME = 1, NTERME
C 
C       Modif CHAMP(NCHAMP,TERME) = COEFF*CHAMP(NCHAMP,TERME)
C 
        CHAMP(NCHAMP, TERME) = CHAMP(NCHAMP, TERME)
     &                       + COEFF*CHAMP(NBCHAR, TERME)
C 
C       FIN DE BOUCLE i SUR TOUS LES TERMES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES FONCTIONS
C 
      END DO
C 
C     Calcul de la contribution de la fonction orthonormee
C 
1000  CONTINUE
C 
      COEFF =  RCOMPO(K)
C 
      DO TERME = 1, NTERME
C 
        CHAMP(NBCHAR, TERME) =   COEFF * CHAMP(NBCHAR, TERME)
C 
      END DO
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Modification des fonctions du temps pour tenir compte des orthogonalisations
C     A priori nterme = npicet
C 
C     On envoie comme arguments :
C 
C     E ...... NTERME    la longueur du  champ
C     E ...... NBCHAR    le nombre de champs ecrits
C     E ...... RCOMPO    les coefficients dus a l'orthogonalisation
C     E ...... NBCHMO    le nombre de champs a modifier
C 
C     Et on recupere :
C 
C     S ...... CHAMP     Le tableau  = TBSTOC (NUCHAM, LONCHA)
C 
      SUBROUTINE MOCHAM (NTERME, NBCHAR, RCOMPO, NBCHMO, CHAMP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NTERME
C 
      INTEGER NBCHAR
C 
      INTEGER NBCHMO
C 
      DOUBLE PRECISION CHAMP(CHAMAX, NTERME)
C 
      DOUBLE PRECISION  RCOMPO(NBCHMO)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER TERME, NCHAMP, K
C 
      DOUBLE PRECISION COEFF
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MOCHAM')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C     Initialisation
C 
CD    CALL IMPTDP ('CHAMP  ', CHAMP(1,1), CHAMAX, NTERM)
CD    CALL IMPTDP ('RCOMPO ', RCOMPO(1), NBCHMO, 1)
CD    CALL IMPEN  ('NTERNE ', NTERME)
CD    CALL IMPEN  ('NBCHAR ', NBCHAR)
CD    CALL IMPEN  ('NBCHMO ', NBCHMO)
C 
C     IF (NBETLT .EQ. 66) THEN
C 
C       CALL IMPTDT ('CHAMP  ', CHAMP(1,1), CHAMAX, NTERME)
C       CALL IMPTDT ('RCOMPO ', RCOMPO(1), NBCHMO, 1)
C       CALL IMPET  ('NTERNE ', NTERME)
C       CALL IMPET  ('NBCHAR ', NBCHAR)
C       CALL IMPET  ('NBCHMO ', NBCHMO)
C 
C     END IF
C 
C     BOUCLE SUR LES FONCTIONS
C 
C     On ne modifie pas la 1ere qui correspond a l'initialisation
C 
      K = 1
C 
C     Addition de la contribution de la fonction orthonormee aux vecteurs precedents
C 
C     Si ????, on ne rentre pas dans la boucle
C 
      IF ((NBCHAR-NBCHMO+1) .GT. (NBCHAR-1)) GOTO 1000
C 
        DO NCHAMP = NBCHAR-NBCHMO+1, NBCHAR-1
C 
C       MODIF COEFF =  (1.D0+RCOMPO(K))
C 
        COEFF =  RCOMPO(K)
        K     =  K+1
C 
C       BOUCLE i SUR TOUS LES TERMES
C 
        DO TERME = 1, NTERME
C 
C       MODIF CHAMP (NCHAMP, TERME) = COEFF*CHAMP (NCHAMP, TERME)
C 
          CHAMP (NCHAMP, TERME) = CHAMP (NCHAMP, TERME)
     &                          + COEFF*CHAMP (NBCHAR, TERME)
C 
        END DO
C 
C       FIN DE BOUCLE i SUR TOUS LES TERMES
C 
C     FIN DE BOUCLE SUR LES FONCTIONS
C 
      END DO
C 
C     Calcul la contribution de la fonction orthonormee
C 
1000  CONTINUE
      COEFF =  RCOMPO(K)
C 
      DO TERME = 1, NTERME
C 
        CHAMP (NBCHAR, TERME) =   COEFF * CHAMP (NBCHAR, TERME)
C 
      END DO
C 
CD    CALL OMPTDP ('VALEUR DU CHAMP ', CHAMP, CHAMAX, NTERME)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Modification des fonctions du temps champs pour tenir compte
C     de l'etape preliminaire  !!!A priori nterme = npicet  !!!
C     Cette modification concerne les fonctions du temps reelles.
C 
C     On envoie comme arguments :
C 
C     E ...... NTERME    la longueur d'un terme
C     E ...... NBCHAR    nombre de champs reels
C     E ...... NBCMOD    Nombre de champs modifies
C     E ...... CHAADD    champ a additioner
C 
C     Et on recupere :
C 
C     ES ...... CHAMP    champ reel en sortie
C 
      SUBROUTINE ADDCHA (NTERME, NBCHAR, NBCMOD, CHAADD, CHAMP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'

C     Nombre de composantes du champ
C 
      INTEGER NTERME
C 
C     Nombre de champs reels deja remplis
C 
      INTEGER NBCHAR
C 
C     Nombre de champs reels a modifier
C 
      INTEGER NBCMOD
C 
C     Tableau des champs
C 
      DOUBLE PRECISION CHAMP(CHARAX, NTERME)
C 
C     Tableau de recomposition du vecteur precedemment orthogonalise
C 
      DOUBLE PRECISION  CHAADD(NTERME, NBCMOD)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER TERME, NCHAMP, NMODIF
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ADDCHA')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPTDP ('CHAMP  ', CHAMP(1,1), CHARAX, NTERME)
CD    CALL IMPTDP ('CHAADD ', CHAADD(1,1), NTERME, NBCMOD)
CD    CALL IMPEP  ('NBCHAR ', NBCHAR)
CD    CALL IMPEP  ('NBCMOD ', NBCMOD)
C 
C     On ne modifie pas la 1ere qui correspond a l'initialisation;
C     addition de la contribution de la fonction orthonormee aux vecteurs precedents
C 
        NMODIF = 0
C 
        DO NCHAMP = NBCHAR-NBCMOD+1, NBCHAR
C 
          NMODIF = NMODIF+1
C 
C         boucle sur tous les termes
C 
          DO TERME = 1, NTERME
C 
            CHAMP(NCHAMP, TERME) = CHAMP(NCHAMP, TERME)  +
     &                             CHAADD(TERME, NMODIF)
C 
          END DO
C 
        END DO
C 
CD    CALL OMPTDP ('VALEUR DU CHAMP ', CHAMP(1,1), CHARAX, NTERME)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C    Cette routine permet de ne garder dans les acroissements admissibles
C    totaux que ceux qui sont significatif i.e. dont la fonction du temps
C    correspondante a une norme relative non negligeable (en effet toutes
C    les fonctions de l'espace sont normees par rapport rapport a K0 ou K0-1)
C 
C    Rappel :
C 
C     5 < = > evolution des contraintes partielles
C     6 < = > evolution des deformations partielles
C     0 < = > deplacement totaux reels
C     7 < = > evolution des deplacements totaux reels
C     8 < = > deformations admissibles totales reelles
C     9 < = > sauts admissibles totaux reels
C    10 < = > evolution des quantites de type contraintes totales reelles
C    11 < = > contraintes totales reelles
C    12 < = > contraintes normales totales reelles
C 
      SUBROUTINE TASCHA
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
CD    LOGICAL  LTRACN , LTRACP
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER EPSCH8, SGNC12, SIGC11, SAUCH9, FTRE10
      INTEGER FTECH6, FTSCH5, FTREE7, DEPCH0
      INTEGER LONEPS, LONDER, LONSAU
      INTEGER DEBNOR, LONRES, DEBFT, NUFT, ADNOR
      INTEGER DEBDPR, DEBEPR, DEBSAR
      INTEGER ADPICE, ADINTE, FINBOU
      INTEGER TABFT(2), TABDEP(2), TABEPS(2), TABSAU(2)
C 
      DOUBLE PRECISION INTERV, NORCOM, TEST
C 
      INTEGER   POSI
      INTEGER   AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TASCHA')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL MESSAO ('ON ENTRE DANS '//IDPROG)
C 
C     Tableaux des fonctions du temps pour les deltas admissibles
C  
      CALL ADTBDM ('TEMPS-SIGM', FTSCH5)
      CALL ADTBDM ('TEMPS-EPSI', FTECH6)
C 
C     (TYPE .EQ. 5) = > NBFSIG
C     (TYPE .EQ. 6) = > NBFEPS
C 
        POSI   = NBFSIG+1
        POSI   = NBFEPS+1
C 
      CALL IMPET ('NBFSIG EN ENTREE DANS '//IDPROG, NBFSIG)
      CALL IMPET ('NBFEPS EN ENTREE DANS '//IDPROG, NBFEPS)
C 
C     Tableaux des deformations et contraintes admissibles TOTALES
C 
      CALL ADTBDM ('DEPLA-ADMI', DEPCH0)
      CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
      CALL ADTBDM ('SIG-AD-TOT', SIGC11)
C 
C     POUR LES INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
C 
      ELSE
C 
        SAUCH9 = EPSCH8
        SGNC12 = SIGC11
C 
      ENDIF
C 
      CALL ADTBDM ('TEMPS-REEL', FTREE7)
      CALL ADTBDM ('TEMP-SI-RE', FTRE10)
C 
      LONEPS = NTETA*NEPS*NGAU1
      LONDER = NDDL*NTETA
      LONSAU = NTETA*NSAU*NGAU2
C 
C    (TYPE .EQ. 0) = > NBDEPR
C    (TYPE .EQ. 7) = > NBDPTR
C    (TYPE .EQ. 8) = > DEADTR
C    (TYPE .EQ. 9) = > SAADTR
C    (TYPE .EQ. 10)= > EVCOTR
C    (TYPE .EQ. 11)= > COTORE
C    (TYPE .EQ. 12)= > CNTORE
C 
      CALL IMPET ('NBDPTR EN ENTREE DANS '//IDPROG, NBDPTR)
      CALL IMPET ('EVCOTR EN ENTREE DANS '//IDPROG, EVCOTR)
C 
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
      CALL ADTBDM ('VALE-TEMPS', ADPICE)
C 
      INTERV = DM(ADINTE+1)
C 
      TABFT(1)  = CHARAX
      TABFT(2)  = NPICET
C 
      TABDEP(1) = CHARAX
      TABDEP(2) = LONDER
C 
      TABEPS(1) = CHARAX
      TABEPS(2) = LONEPS
C 
      TABSAU(1) = CHARAX
      TABSAU(2) = LONSAU
C 
C     Traitement des deformations admissibles
C 
      NORCOM     = 0.D0
C 
      LONRES = NBDPTR-1+NPICET
      CALL POUSMD (LONRES, DEBNOR)
      DEBFT  = DEBNOR+NBDPTR-1
      DEBNOR = DEBNOR -2
C 
C     Calcul de toutes les normes au carre en temps des deformations admissibles
C     a zero et de la norme au carre
C 
CD    CALL IMPTDN ('FONCTION DU TEMPS EN DEFORMATION EN ENTREE '
CD                 //IDPROG, DM(FTREE7), CHARAX, NPICET)
C 
      DO NUFT = 2, NBDPTR
C 
CD      CALL IMPEN ('POUR LA DEFORMATION NUMERO '//idprog, NUFT)
C 
        ADNOR  =  DEBNOR+NUFT
C 
        CALL EXTRAD (DM(FTREE7), 2, TABFT(1),
     &               2, NUFT, DM(DEBFT), NPICET)
C 
CD      CALL IMPTDN ('FONCTION DU TEMPS EXTRAITE '
CD                   //IDPROG, DM(DEBFT), 1, NPICET)
C 
        CALL SCALFT (ADINTE, ADPICE, DM(DEBFT), DM(DEBFT), DM(ADNOR))
C 
CD      CALL IMPDN ('LA NORME EN TEMPS VAUT '//idprog, DM(ADNOR))
C 
        NORCOM = DMAX1(DM(ADNOR), NORCOM)
C 
      END DO
C 
CD    CALL IMPDN ('NORME MAXI EN DEFORMATION ', NORCOM)
C 
C     Test pour savoir si la deformation admissible est a garder
C     si oui on la recopie des qu'il y a de la place.
C     Comme on s'interesse au champ reel on conserve de toute facon
C     le premier champ qui correspond a la solution elastique.
C 
      FINBOU  = NBDPTR
C 
      NBDPTR = 1
      NBDEPR = 1
      DEADTR = 1
      SAADTR = 1
C 
C     Comme on s'interesse au champ reel on ne conserve pas le premier champ.
C 
      NBFEPS = 0
C 
      TEST = 1.D -4*NORCOM
C 
CD    CALL IMPDN ('DEFORAMATIONS ELIMINEES EN DESSOUS DE ', TEST)
C 
      LONRES = LONEPS + LONDER + LONSAU
C 
      CALL POUSMD (LONRES, DEBDPR)
C 
      DEBEPR = DEBDPR+LONDER
      DEBSAR = DEBEPR+LONEPS
C 
      DO NUFT = 2 , FINBOU
C 
        ADNOR = NUFT+DEBNOR
C 
        IF (DM(ADNOR) .LT. TEST) THEN
C 
CD        CALL IMPDP ('VALEUR DE LA NORME ', DM(ADNOR))
C 
CD        CALL IMPET
CD          ('ELIMINATION DE LA DEFORMATION REELLE NUMERO '//IDPROG, NUFT)
C 
        ELSE
C 
C     La valeur des fonctions du temps partielles n'a aucune importance;
C     ici on en calcule uniquement le nombre residuel.
C 
CD        CALL IMPET
CD          ('COPIE DE LA DEFORMATION REELLE NUMERO '//IDPROG, NUFT)
C 
          NBFEPS = NBFEPS +1
          NBFEPS = MIN0(CHAMAX, NBFEPS)
C 
          CALL EXTRAD (DM(FTREE7), 2, TABFT(1),
     &                 2, NUFT, DM(DEBFT), NPICET)
C 
          CALL RCHARE (7, NPICET, DM(DEBFT), DM(FTREE7))
C 
          CALL EXTRAD (DM(DEPCH0), 2, TABDEP(1),
     &                 2, NUFT, DM(DEBDPR), LONDER)
C 
          CALL RCHARE (0, LONDER, DM(DEBDPR), DM(DEPCH0))
C 
          CALL EXTRAD (DM(EPSCH8), 2, TABEPS(1),
     &                 2, NUFT, DM(DEBEPR), LONEPS)
C 
          CALL RCHARE (8, LONEPS, DM(DEBEPR), DM(EPSCH8))
C 
          CALL EXTRAD (DM(SAUCH9), 2, TABSAU(1),
     &                 2, NUFT, DM(DEBSAR), LONSAU)
C 
          CALL RCHARE (9, LONSAU, DM(DEBSAR), DM(SAUCH9))
C 
        END IF
C 
      END DO
C 
C     Traitement des contraintes admissibles
C 
      NORCOM     = 0.D0
C 
      LONRES = EVCOTR-1+NPICET
      CALL POUSMD (LONRES, DEBNOR)
      DEBFT  = DEBNOR+EVCOTR-1
      DEBNOR = DEBNOR -2
C 
C     Calcul de toutes les normes au carre en temps des deformations admissibles
C     a zero et de la norme au carre.
C 
      DO NUFT = 2 , EVCOTR
C 
        ADNOR  =  DEBNOR+NUFT
C 
        CALL EXTRAD (DM(FTRE10), 2, TABFT(1),
     &               2, NUFT, DM(DEBFT), NPICET)
C 
        CALL SCALFT (ADINTE, ADPICE, DM(DEBFT), DM(DEBFT), DM(ADNOR))
C 
CD      CALL IMPEN ('POUR LA CONTRAINTE NUMERO '//IDPROG, NUFT)
CD      CALL IMPDN ('LA NORME EN TEMS VAUT     '//IDPROG, DM(ADNOR))
C 
        NORCOM = DMAX1( DM(ADNOR) , NORCOM)
C 
      END DO
C 
CD    CALL IMPDN ('NORME MAXI EN CONTRAINTE ', NORCOM)
C 
C     Test pour savoir si la contrainte admissible est a garder;
C     si oui on la recopie des qu'il y a de la place.
C     Comme on s'interesse au champ reel on conserve de toute facon
C     le premier champ qui correspond a la solution elastique.
C 
      FINBOU = EVCOTR
C 
      EVCOTR = 1
      COTORE = 1
      CNTORE = 1
C 
C     Comme on s'interesse au champ reel on ne conserve pas le premier champ.
C 
      NBFSIG = 0
C 
C     Comme on s'interesse au champ partiel on ne conserve pas le premier champ
C 
      TEST = 1.D -4*NORCOM
C 
CD    CALL IMPDN ('CONTRAINTES ELIMINEES EN DESSOUS DE ', TEST)
C 
      DO NUFT = 2 , FINBOU
C 
        ADNOR = NUFT+DEBNOR
C 
        IF (DM(ADNOR) .LT. TEST) THEN
C 
CD       CALL IMPDP ('VALEUR DE LA NORME ', DM(ADNOR))
C 
         CALL IMPET
     &  ('ELIMINATION DE LA CONTRAINTE REELLE NUMERO '//IDPROG, NUFT)
C 
        ELSE
C 
C     La valeur des fonctions du temps partielles n'a aucune importance;
C     ici on en calcule uniquement le nombre residuel.
C 
CD       CALL IMPEN
CD      ('COPIE DE LA DEFORMATION REELLE NUMERO '//IDPROG, NUFT)
C 
          NBFSIG = NBFSIG+1
          NBFSIG = MIN0(CHAMAX, NBFSIG)
C 
          CALL EXTRAD (DM(FTRE10), 2, TABFT(1),
     &                 2, NUFT, DM(DEBFT), NPICET)
C 
          CALL RCHARE (10, NPICET, DM(DEBFT), DM(FTRE10))
C 
          CALL EXTRAD (DM(SIGC11), 2, TABEPS(1),
     &                 2, NUFT, DM(DEBEPR), LONEPS)
C 
          CALL RCHARE (11, LONEPS, DM(DEBEPR), DM(SIGC11))
C 
          CALL EXTRAD (DM(SGNC12), 2, TABSAU(1),
     &                 2, NUFT, DM(DEBSAR), LONSAU)
C 
          CALL RCHARE (12, LONSAU, DM(DEBSAR), DM(SGNC12))
C 
        END IF
C 
      END DO
C 
      CALL IMPET ('NBFEPS EN SORTIE DANS '//IDPROG ,NBFEPS)
      CALL IMPET ('NBFSIG EN SORTIE DANS '//IDPROG ,NBFSIG)
C 
      CALL IMPET ('NBDPTR EN SORTIE DANS '//IDPROG ,NBDPTR)
      CALL IMPET ('EVCOTR EN SORTIE DANS '//IDPROG, EVCOTR)
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
C    Cette routine est appelee dans le cas de l'option de reprise BIGNET = .TRUE.
C    Alors on fait un grand nettoyage pour pouvoir continuer le calcul, c'est a dire
C    qu'on reinitialise les tableaux de rangement des champ d'accroissements des
C    contraintes et deformations ainsi que des fonctions du temps associees. On repart
C    de la solution admissible totale a l'etape globale precedente.
C 
C    Rappel :
C 
C     5  < = > evolution des contraintes partielles
C     6  < = > evolution des deformations partielles
C     0  < = > deplacement totaux reels
C     7  < = > evolution des deplacements totaux reels
C     8  < = > deformations admissibles totales reelles
C     9  < = > sauts admissibles totaux reels
C     10 < = > evolution des quantites de type contraintes totales reelles
C     11 < = > contraintes totales reelles
C     12 < = > contraintes normales totales reelles
C 
      SUBROUTINE DELCHA
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
      INTEGER   AM2LC, ADM2LC, FTSCH5, FTECH6, LONRES, DEBUT
      INTEGER   LONEPS, LONSAD, DEBNET, LONECF, IUNIT1
      CHARACTER*10 NOM
      CHARACTER*3  CARNET
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DELCHA')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NBFEPS = 0
      NBDPTR = 1
      NBDEPR = 1
      DEADTR = 1
      SAADTR = 1
      NBFSIG = 0
      EVCOTR = 1
      COTORE = 1
      CNTORE = 1
C 
C     Nettoyage des tableaux des fonctions du temps
C 
      LONRES = NBFEPS * NPICET
      CALL ADTBDM ('TEMP-SI-RE', FTSCH5)
      CALL MENADM (FTSCH5, LONRES)
      CALL ADTBDM ('TEMPS-REEL', FTECH6)
      CALL MENADM (FTECH6, LONRES)
      CALL ADTBDM ('TEMPS-SIGM', FTSCH5)
      CALL MENADM (FTSCH5, LONRES)
      CALL ADTBDM ('TEMPS-EPSI', FTECH6)
      CALL MENADM (FTECH6, LONRES)
C  
C     Nettoyage des tableaux des champs
C 
      LONEPS = NEPS*NTETA*NGAU1*CHAMAX
      LONSAD = NTETA*NSAU*NBINT*XINTEG*NBCOL*CHARAX
      CALL ADTBDM ('EPS-AD-TOT', DEBUT)
      CALL MENADM (DEBUT, LONEPS)
      CALL ADTBDM ('SIG-AD-TOT', DEBUT)
      CALL MENADM (DEBUT, LONEPS)
      CALL ADTBDM ('SAU-AD-TOT', DEBUT)
      CALL MENADM (DEBUT, LONSAD)
      CALL ADTBDM ('SGN-AD-TOT', DEBUT)
      CALL MENADM (DEBUT, LONSAD)
C 
C     En sortant, on retablit BIGNET a .FALSE.
C 
      BIGNET = .FALSE.
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
