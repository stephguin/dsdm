C     On peut confondre  EFFDEV et DEPDEV en entree
C 
C     Calcul en isotrope transverse avec un second membre
C     ranges (NDDL, NBMAT).
C  
C     On envoie comme arguments :
C 
C     ES...... EFFDEV  Valeur des efforts developpes
C     S ...... TNUDEV  tableau des numeros de developpements non nuls
C     S ...... NBDVV   nombre de developpements non nuls
C     S ...... EFFDEV  Valeur des efforts developpes
C 
C     Et on recupere :
C 
C     ES...... DEPDEV  Valeur des deplacements solution developpee rangee :
C                      (NDDL, NBMAT)
C 
      SUBROUTINE CALIS1 (EFFDEV, NBDEV, TNUDEV, DEPDEV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'strategie_calcul.h'
C 
      DOUBLE PRECISION   EFFDEV(NBMAT*NDDL)
C 
      INTEGER            NBDEV
      INTEGER            TNUDEV(NBDEV)
C 
      DOUBLE PRECISION   DEPDEV(NBMAT*NDDL)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER          TNUENR
      INTEGER          I, DEBEFF, DEBMAT, ADPRO
      INTEGER          IUNIT, DEPLA
      INTEGER          DECAL
      INTEGER          AM2LC, ADM2LC, JCE
      INTEGER          NUENRS, PAS
C 
      INTEGER          NTRESO,  NIMAT, NBRES
      INTEGER          DEBMA(NBCE), DEBEF(NBCE), DEBDE(NBCE)
C 
CD    LOGICAL          LTRACP, LTRACN
C 
      CHARACTER*6 IDPROG  , NOM
      PARAMETER (IDPROG='CALIS1')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL POUSME (NBDEV, TNUENR)
      CALL NUMACA (NBDEV, TNUDEV, M(TNUENR))
      TNUENR = TNUENR-1
C 
C     Mise a zero de depdev; seuls les termes non nuls seront calcules
C  
      NTRESO = NDDL*NBMAT
      CALL BALAID (NTRESO, DEPDEV(1))
C 
C     RESERVATION DE PLACE POUR K0N SANS REMISE A ZERO
C 
      DECAL =NDDL
C 
      NTRESO = NBCE*NTMAT
      CALL POUSMD (NTRESO, DEBMAT)
C 
      CALL ADTBM ('PRODL     ', ADPRO)
      NOM      = 'matiso'
C 
C     ouverture du fichier des matrices
C 
      CALL OFDDNF (1, NOM, 6, NTMAT, IUNIT)
C 
      NIMAT = (NBDEV/NBCE) * NBCE
C 
      PAS = 0
C 
      IF (NIMAT .GT. 0) THEN
C 
C         DO I = 1, NIMAT, NBCE
C         SG PB BOUCLE DO
C 
        DO I = 1, NIMAT
          DEBEFF   = 1+(TNUDEV(I)-1)*NDDL
          DEBEF(1) = DEBEFF
          DEBMA(1) = DEBMAT
          DEBDE(1) = DEBEFF
C 
          DO JCE = 2, NBCE
C 
CD          CALL IMPET('1 NUM MAT DANS CALIS1 JCE; ', JCE)
C 
            DEBEF(JCE) = 1+(TNUDEV( I+JCE-1) -1)*NDDL
            DEBMA(JCE) = DEBMA(JCE-1) + NTMAT
            DEBDE(JCE) = DEBEF(JCE)
          ENDDO
C 
CD        N  = I-1-NTDSFG
C 
CD        CALL IMPEN ('POUR N ', N)
C 
C         LECTURE a partir de DM(DEBMAT) du fichier K0N
C         hypothese : les fichiers K0N ont ete les premiers a etre ecrits
C 
          DO JCE = 1, NBCE
C 
C           Lecture de la matrice correspondant au numero
C 
            PAS    = PAS+1
C   
C           MODIF          NUENRS = TNUDEV(PAS)
C 
            NUENRS = M(TNUENR+PAS)
C 
CD          CALL IMPET('2 NUM MAT DANS CALIS1 JCE;',JCE)
CD          CALL IMPET('3 DEBMA(JCE) ',DEBMA(JCE) )
CD          CALL IMPET('4 NUENRS ',NUENRS )
CD          CALL IMPET('5 NTMAT ', NTMAT )
CD          CALL IMPDT('6 DM(DEBMA(JCE)) ', DM(DEBMA(JCE)) )
C 
            CALL LFDDNF (DM(DEBMA(JCE)), DEBMA(JCE),
     &                   NTMAT, IUNIT, NUENRS)
C 
          ENDDO
C 
C         Resolution par K0N pour toutes les fonctions du temps :
C 
C         vd$l concur
C         vd$l cncall
C 
          DO JCE = 1, NBCE
            CALL DEREPP (NDDL, M(ADPRO), DM(DEBMA(JCE)), NFTGLO,
     &                 EFFDEV(DEBEF(JCE)), DEPDEV(DEBDE(JCE)))
          ENDDO
C 
         DEPLA  = DEPLA + NBCE*DECAL
C 
        END DO
C 
      END IF
C 
C     LE RESTE
C 
      NBRES = NBDEV - NIMAT
C 
      IF (NBRES .GT. 0) THEN
C 
C       OA        DO I = NIMAT+1 , NBDEV , NBRES
C 
C       OA           DEBEFF = 1+(TNUDEV( I) -1 )*NDDL
C       OA           DEBEF(1) = DEBEFF
C       OA           DEBMA(1) = DEBMAT
C       OA           DEBDE(1) = DEBEFF
C 
C       OA           DO JCE = 2,NBRES
C       OA             DEBEF(JCE) = 1+( TNUDEV( I+JCE-1) -1 )*NDDL
C       OA             DEBMA(JCE) = DEBMA(JCE-1) + NTMAT
C       OA             DEBDE(JCE) = DEBEF(JCE)
C       OA           ENDDO
C 
        DEBEF(1) = 1+( TNUDEV(NIMAT+1)-1)*NDDL
        DEBMA(1)   = DEBMAT
        DEBDE(1)   = DEBEF(1)
C 
        DO JCE = 2, NBRES
C 
          DEBEF(JCE) = 1+( TNUDEV(NIMAT+JCE) -1 )*NDDL
          DEBMA(JCE) = DEBMA(JCE-1) + NTMAT
          DEBDE(JCE) = DEBEF(JCE)
C 
        END DO

C 
C     LECTURE a partir de DM(DEBMAT) du fichier K0N
C     hypothese : les fichiers K0N ont ete les premiers a etre ecrits
C 
        DO JCE = 1, NBRES
          PAS    = PAS +1
C     MODIF          NUENRS = TNUDEV(PAS)
          NUENRS = M(TNUENR+PAS)
C 
CD         CALL IMPET('APRES RESTE NUM MAT DANS CALIS1 JCE;',JCE)
C 
          CALL LFDDNF (DM(DEBMA(JCE)), DEBMA(JCE), NTMAT, IUNIT, NUENRS)
C 
        ENDDO
C 
C     Resolution par K0N pour toutes les fonctions du temps :
C 
C     vd$l concur
C     vd$l cncall
C 
        DO JCE = 1,NBRES
          CALL DEREPP( NDDL , M(ADPRO), DM(DEBMA(JCE)), NFTGLO ,
     &                 EFFDEV(DEBEF(JCE)) , DEPDEV(DEBDE(JCE)) )
        ENDDO
C 
      END IF
C 
      CALL FERFIC (1, IUNIT, IDPROG)
C 
C     Mise a zero stricte des ddl bloques
C 
      CALL MZDDLB (NFTGLO, DEPDEV)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
