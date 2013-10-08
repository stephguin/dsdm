C     Cette routine remplit le second membre pour l'optimisation
C     par BFGS en espace. Soient Fi les efforts globaux et
C     Vi les deplacements associes par le calcul isotrope transverse,
C     le second membre correspond a fik = INT (VK . Fi), k variant d'abord.
C 
C     On envoie comme arguments :
C 
C     E ...... DEPDV1 le developpement des 1er deplacements
C                     stocke (NDDL, NBMAT)
C     E ...... DEPDV2 le developpement des seconds deplacements
C                     stocke (NDDL, NBMAT)
C 
C     Et on recupere :
C 
C     S ...... SCAL

      SUBROUTINE SCADEV (DEPDV1, DEPDV2, SCAL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION    DEPDV1(NDDL*NBMAT)
      DOUBLE PRECISION    DEPDV2(NDDL*NBMAT)
      DOUBLE PRECISION    SCAL
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      DOUBLE PRECISION RESU
C 
      INTEGER DEBUT, NTOT
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SCADEV')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     SEQUENCE D'IMPRESSION EN TRACE PROFONDE
C 
CD    IF ( LTRACP(1) )THEN
C 
CD      CALL OMPTDP ('DEPDV1 (NDDL,NBMAT) ', DEPDV1(1), NDDL, NBMAT)
CD      CALL OMPTDP ('DEPDV2 (NDDL*NBMAT) ', DEPDV2(1), NDDL, NBMAT)
C 
CD    END IF
C 
      DEBUT = 1
      SCAL =  0.D0
C 
      NTOT = NDDL*NTDSFG
      CALL SCAVEC (NTOT, DEPDV1(DEBUT), DEPDV2(DEBUT), RESU)
C 
C     TCD CALL VPSCAL (DEPDV1(DEBUT), DEPDV2(DEBUT), RESU, NTOT)
C 
      SCAL = SCAL+ PI*RESU
      DEBUT = DEBUT + NTOT
C 
      CALL SCAVEC (NDDL, DEPDV1(DEBUT), DEPDV2(DEBUT), RESU)
C 
C     TCD CALL VPSCAL (DEPDV1(DEBUT), DEPDV2(DEBUT), RESU, NDDL)
C 
      SCAL = SCAL+ 2.D0*PI*RESU
      DEBUT = DEBUT + NDDL
C 
      CALL SCAVEC (NTOT, DEPDV1(DEBUT), DEPDV2(DEBUT), RESU)
C 
C     TCD CALL VPSCAL (DEPDV1(DEBUT), DEPDV2(DEBUT), RESU, NTOT)
C 
      SCAL = SCAL+ PI*RESU
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1)) THEN
C 
CD       CALL IMPDN ('produit scalaire ', SCAL)
C 
CD    END IF
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
