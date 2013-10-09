C     On envoie comme arguments :
C 
C     E ...... TABGAU  Tableaux N pour le point de gauss concerne
C     E ...... N       Le numero de developpement concerne
C     E ...... A       La demi-longeur de l'element
C     E ...... B       La demi-hauteur de l'element
C     E ...... R       Le rayon au point de gauss de l'element
C     E ...... VALDEP  La valeur des deplacements sur l'element
C                      !!!!! ranges calcul ===== > PASSAGE PAR RAN-INV
C     Et on recupere :
C 
C     S ...... EPSILN  Valeur de la composante N des deformations
C                      rangee :
C                      (e11,e22,r2*e12,r2*e23,r2*e13,e33)
C                      attention  !!!!!!!! le developpement est du
C                      type suivant:
C                      N > OU = 0  ==> (COSN,COSN,COSN,SINN,SINN,CONS)
C                      N <         ==> (SINN,SINN,SINN,COSN,COSN,SINN)
C 
      SUBROUTINE EPSDEV (TABGAU, N, A, B, R, VALDEP, EPSILN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER          N
      DOUBLE PRECISION TABGAU(36) , A , B , R , VALDEP(36) , EPSILN(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER           AD2  ,ADD
      INTEGER           AM2LC , ADM2LC , I
      DOUBLE PRECISION  NLOC , R2 , EPSLOC , EPLOC(2)
C 
CD     CHARACTER*6 IDPROG
CD     PARAMETER (IDPROG='EPSDEV')
C 
CD    CALL WLKBCD(IDPROG)
C 
C     VALEURS EN ENTREE
C 
CD    CALL OMPTDP(' VALEURS DE TABGAUSS ' , TABGAU(1) , 36 , 1 )
CD    CALL OMPTDP(' VALEURS DES DEPLACEMENTS RANGES CALCULS   '
CD                   , VALDEP(1) , 36 , 1 )
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
CD    CALL IMPDP( ' VALEUR DE A ',A)
CD    CALL IMPDP( ' VALEUR DE B ',B)
CD    CALL IMPDP( ' VALEUR DE R ',R)
C 
      NLOC       = DBLE(N)
      R2         = DSQRT(2.D0)/2.D0
C 
      CALL POUSMD (36 , AD2)
      AD2 = AD2 -1
C 
      DO I = 1  , 6
        EPSILN(I) = 0.D0
      END DO

      DO I = 1  , 6
        DM (AD2+I) = TABGAU(I)
      END DO
      DO I = 7  , 12
        DM (AD2+I) = A*TABGAU(I)
      END DO
      DO I = 13  , 18
        DM (AD2+I) = TABGAU(I)
      END DO
      DO I = 19  , 24
        DM (AD2+I) = A*TABGAU(I)
      END DO
      DO I = 25  , 30
        DM (AD2+I) = TABGAU(I)
      END DO
      DO I = 31  , 36
        DM (AD2+I) = A*TABGAU(I)
      END DO
C 
      ADD = AD2+12
      DO I = 1 , 12
        EPSILN(1)  = EPSILN(1)+DM(ADD +I)*VALDEP(I)
      END DO
      EPSILN(1)  = EPSILN(1)/A
C 
      EPLOC(1) =0.D0
      EPLOC(2) =0.D0
C 
      DO I = 1 , 12
        EPLOC(1)  = EPLOC(1)+DM(AD2+I)*VALDEP(I)
      END DO
      DO I = 1 ,12
        EPLOC(2)  = EPLOC(2)+DM(AD2+I)*VALDEP(I+12)
      END DO
C 
      EPSILN(2)  = (NLOC*EPLOC(2)+EPLOC(1))/R
C 
      DO I = 13 ,24
        EPSILN(3)  = EPSILN(3)+DM(AD2+I)*VALDEP(I)
      END DO
C 
      EPSILN(3)  = R2*( (-NLOC*EPLOC(1)-EPLOC(2))/R + EPSILN(3)/A )
C 
      EPSLOC = 0.D0
C 
      DO I = 13 ,24
        EPSLOC  = EPSLOC+DM(ADD +I)*VALDEP(I)
      END DO
C 
      DO I = 1 ,12
        EPSILN(4)  = EPSILN(4)+DM(AD2+I)*VALDEP(24+I)
      END DO
C 
      EPSILN(4)  = R2*( EPSLOC/B -NLOC*EPSILN(4)/R)
C 
      DO I = 1 ,12
        EPSILN(5)  = EPSILN(5)+DM(ADD +I)*VALDEP(I+24)
      END DO
C 
      ADD = AD2+24
C 
      EPSLOC = 0.D0
C 
      DO I = 1 ,12
        EPSLOC  = EPSLOC+DM(ADD +I)*VALDEP(I)
      END DO
C 
      EPSILN(5)  = R2*( EPSLOC/B +EPSILN(5)/A)
C 
      DO I = 25 ,36
        EPSILN(6) = EPSILN(6)+ DM(AD2+I)*VALDEP(I)
      END DO
      EPSILN(6)  = EPSILN(6)/B
C 
CD    CALL IMPTDN ('VALEUR DE EPSILON(N) ', EPSILN(1) , 1 , 6 )
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
C     POUR L'INITIALISATION ELASTIQUE
C 
C     CALCUL EN ISOTROPE TRANSVERSE : !!!! meme s'il y a
C     plusieurs seconds membres, on fait le calcul fonction du
C     temps par fonction du temps !!!!
C 
C     On envoie comme arguments :
C 
C     E ...... EFFDEV  Valeur des efforts developpes
C 
C     Et on recupere :
C 
C     S ...... DEPDEV  Valeur des deplacements solutions developpes

      SUBROUTINE CALISO (EFFDEV, DEPDEV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION   EFFDEV(NBMAT*NDDL*NFOTPS)
      DOUBLE PRECISION   DEPDEV(NBMAT*NDDL*NFOTPS)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I, DEBEFF, DEBMAT, ADPRO, EFFORT
      INTEGER  IN
      INTEGER  DEPLA, N, NFT, IUNIT, ADEPPA, NMAT
      INTEGER  AM2LC, ADM2LC
C 
CD    INTEGER  VERDEP, TEST1, NBTDDL, VERTRA
CD    LOGICAL  LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER  (IDPROG='CALISO')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     RESERVATION DE PLACE POUR K0N SANS REMISE A ZERO
C -----------------------------------------------------------------------
C 
      CALL POUSMD (NTMAT, DEBMAT)
      CALL ADTBM ('PRODL     ', ADPRO)
C 
      DEBEFF=1
      CALL GSPOUD (2*NDDL, DEPLA)
      EFFORT=DEPLA+NDDL
C 
      CALL OFDDNF (1, 'matiso', 6, NTMAT, IUNIT)
C 
      CALL INFOEN ('NUMDEV-MAT', ADEPPA, NMAT)
C 
      DO IN = 1, NMAT
C 
        N = M(ADEPPA +IN -1)
        I = N+NTDSFG+1
        DEBEFF = 1+NDDL*(I-1)
C 
C -----------------------------------------------------------------------
C     ECRITURE A PARTIR DE DM(DEBMAT) DU FICHIER K0N
C     HYPOTHESE : LES FICHIERS K0N ONT ETES LES PREMIERS A ETRE ECRITS
C -----------------------------------------------------------------------
C 
        CALL LFDDNF (DM(DEBMAT), DEBMAT, NTMAT, IUNIT, IN)
C 
C -----------------------------------------------------------------------
C     RESOLUTION PAR K0N POUR TOUTES LES FONCTIONS DU TEMPS :
C -----------------------------------------------------------------------
C 
        DO NFT =  1 , NFOTPS
          CALL DESREM (NDDL, M(ADPRO), DM(DEBMAT), EFFDEV(DEBEFF),
     &                 DEPDEV(DEBEFF) )
C 
CD         CALL PRODLV (NDDL, M(ADPRO), DM(DEBMAT), DEPDEV(DEBEFF),
CD                      DM(DEPLA))
C 
          DEBEFF = DEBEFF+NDDL*NBMAT
        END DO
C 
      END DO
C 
      IF (DONFIC .AND. TRACT) THEN
CD       CALL MESSAO ('APPEL A VEDIFT DANS '//IDPROG)
        CALL VEDIFT (0, DEPDEV(1), NDDL*NBMAT)
      END IF
C 
      CALL FERFIC (1, IUNIT, IDPROG)
C 
CD    IF (DONFIC) THEN
C 
CD      CALL MESSAO ('VERIFICATION DES DEP DANS '//IDPROG)
C 
CD      CALL GSPOUD (NDDL*NTETA, VERDEP)
C 
CD      CALL ADTBM ('DDL-BLOQUE', TEST1)
CD      TEST1 = TEST1 -1
CD      CALL LONGEN ('DDL-BLOQUE', NBTDDL)
C 
CD      CALL VRTSM (1, DEPDEV(1), DM(VERDEP))
C 
CD      IF (TRACT) THEN
CD        CALL VEDIFT (1, DM(VERDEP), NDDL*NTETA)
CD      END IF
C 
CD    END IF
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
C     Pour une fonction du temps :
C     Calcule les deformations admissibles pour tous les points de Gauss;
C     elles sont rangees (neps, nteta, ngaus1, nfotps)
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT Nombre de fonctions du temps
C     E ...... DEPDEV Tableau des deplacements developpes
C 
C     Et on recupere :
C 
C     S ...... EPSADM Tableau des deformations reelles (admissibles)
C                     (neps, nteta, ngau1, nfonct)
C 
      SUBROUTINE TOUEPS (NFONCT, DEPDEV, EPSADM)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER            NFONCT
      DOUBLE PRECISION   DEPDEV(NBMAT*NDDL*NFONCT)
      DOUBLE PRECISION   EPSADM(NEPS*NGAU1*NTETA*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  DEPINV, ADRGAU, TABNIV(2), LONGAU
      INTEGER  ADM2DP, ADM2EP, ADM2VE, NUCOU, NUCOL, TLOCN1
      INTEGER  DBDDLU, DBDDLV, DBDDLW, DEBGAU, X, Y, NUDEV
      INTEGER  DBDEPN, DBEPSN, DECAL, DBPREE, NFT
      INTEGER  AM2LC, ADM2LC
      DOUBLE PRECISION  RAYONC, A, B, R
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TOUEPS')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     RECHERCHE DES DIFFERENTS TABLEAUX UTILES
C 
      CALL ADTBM ('DEPINV-COU', DEPINV)
C 
CD    CALL IMPTEP('RANINV ', M(DEPINV) , 1 ,12)
C 
      CALL ADTBM  ('TLOCN1    ', TLOCN1)
      CALL ADTBDM ('TAB-GAUSS ', DEBGAU)
C 
C     POUR UTILISER EXTRAD
C 
      TABNIV(1)   = 6
      TABNIV(2)   = NBMAT
C 
C     Ouverture d'un tableau partiel pour ranger les numeros
C     des ddl des deplacements croissants (u, v, w) par developpement croissant
C 
      CALL GSPOUE (NDDLEL, DBDDLU)
      DBDDLV = DBDDLU+12
      DBDDLW = DBDDLV+12
C 
C 
C     Ouverture d'un tableau partiel pour ranger les deplacements
C     calcules (u,v,w) par developpement croissant, puis les deformations
C 
C     MODIF DU 21-12   CALL GSPOUD (NBMAT*(NDDLEL+NEPS+1)*NFONCT, ADM2DP)
C 
      CALL POUSMD( NBMAT*(NDDLEL+NEPS+1)*NFONCT, ADM2DP )
      ADM2EP = ADM2DP+NBMAT*NDDLEL
      ADM2VE = ADM2EP+NBMAT*NEPS
C 
C     Decalage d'indice pour epsree
C 
      DECAL       = NEPS*NTETA
      DBPREE      = 1
C 
      DO NFT = 1, NFONCT
C 
CD      CALL IMPEN(' NUMERO DE FONCTION DU TEMPS ' , NFT )
C 
        DO NUCOU = 1, NBCOU
C 
CD        CALL IMPEN('POUR NUCOU = ' , NUCOU)
C 
          CALL VALEP (NUCOU, B)
C 
CD        CALL IMPDP ('VALEUR DE B ', B)
C 
          DO NUCOL = 1, NBCOL
            CALL VALRAY (NUCOL, RAYONC, A)
C 
CD          CALL IMPDP ('VALEUR DE A ', A)
C 
C     RECHERCHE  DES NUMEROS DE DDL RANGES CROISSANT POUR L'ELEMENT
C 
            CALL NDDLC (1 , NUCOU, NUCOL, TLOCN1, M(DBDDLU))
            CALL NDDLC (2 , NUCOU, NUCOL, TLOCN1, M(DBDDLV))
            CALL NDDLC (3 , NUCOU, NUCOL, TLOCN1, M(DBDDLW))
C 
C     REMPLISSAGE DE  DEPCAN : Tableau des valeurs des deplacements pour
C     l'element ranges calcul (U, V, W) par ordre croissant de numero de
C     developpement, c-a-d :  - NTDSFG ------- > NTDSFG
C 
            CALL DPDVCA (DEPDEV(1+(NFT-1)*NDDL*NBMAT), M(DEPINV),
     &               M(DBDDLU), M(DBDDLV), M(DBDDLW), DM(ADM2DP))
C 
C     DEBGAU est l'adresse pour x=1 , y=1 dans TAB-GAUSS
C 
            ADRGAU = DEBGAU
C 
C     Boucle sur les points de gauss
C 
            DO X = 1 , XINTEG
C 
CD            CALL IMPEP ( 'VALEUR DE X ', X )
C 
C     Calcul du rayon au point de gauss
C 
              R  = RAYONC + A*GAUSS( XINTEG*(XINTEG-1)/2 + X )
              DO Y = 1 , YINTEG
C 
CD            CALL IMPEP ( 'VALEUR DE Y', Y )
C 
C     Pour le point de gauss considere (==> tabgau considere), on
C     calcule la valeur de epsilon (n)
C  
C     DBDEPN est l'adresse pour n du tableau des deplacements provenant
C     de DPDVCA
C 
                DBDEPN = ADM2DP
C 
C     DBEPSN est l'adresse pour n du tableau des deformations
C 
                DBEPSN = ADM2EP
                DO NUDEV = -NTDSFG , NTDSFG
C 
CD                CALL IMPEP ( 'VALEUR DE NUDEV ', NUDEV )
C 
                  CALL EPSDEV (DM(ADRGAU), NUDEV, A, B, R,
     &                         DM(DBDEPN), DM(DBEPSN))
                  DBDEPN = DBDEPN+36
                  DBEPSN = DBEPSN+6
                ENDDO
                ADRGAU = ADRGAU+36
C 
C     On appelle epsree qui pour un point de gauss donne et a partir d'un
C     rangement calcule des deplacements et fournit la valeur des deformations
C     rangees (neps, nteta, npgau1)
C 
                CALL EPSREE (DM(ADM2EP), EPSADM(DBPREE))
                DBPREE = DBPREE+DECAL
              ENDDO
            ENDDO
          END DO
        ENDDO
      ENDDO
C 
C     POUR VERIFICATION QU'ON NE DEBORDE PAS DANS TABGAUSS
C 
      CALL LONGDP ('TAB-GAUSS', LONGAU)
      CALL TESTEN (ADRGAU-DEBGAU, LONGAU, 'TOUEPS')
      IF ( (DBPREE-1-DECAL*NGAU1*NFONCT).NE.0)THEN
        CALL ERREUD
     $(0,'MAUVAIS CALCUL DE LONGUEUR DE TABLEAU DANS '//IDPROG)
      ENDIF
C 
CD    DO NFT = 1 , NFONCT
CD      CALL IMPEN( ' NUMERO DE FONCTION DU TEMPS ', NFT )
CD      CALL IMPMN(' N > OU = 0  ==> (COSN,COSN,COSN,COSN,COSN,CONS)')
CD      CALL IMPMN('  N <  0     ==> (SINN,SINN,SINN,SINN,SINN,SINN)')
CD      CALL IMPMN(' (e11,e22,r2*e12,r2*e23,r2*e13,e33) ')
CD      CALL OMPTDN(' DEFORMATIONS ADMISSIBLES RANGEES 6*(NTETA *NGAU1)'
CD                   , EPSADM(1+(NFT-1)* NEPS*NTETA*NGAU1)
CD                   , NEPS , NTETA*NGAU1 )
CD    END DO
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
C     On envoie comme arguments :
C 
C     E ...... COMISO(6) Caracteristiques du comportement isotrope
C                        transverse de la couche
C     E ...... EPSILO(6) Valeurs des deformations
C                        (e11, e22, r2*e12, r2*e23, r2*e13, e33)
C     Et on recupere :
C 
C     S ...... SIGMA(6)  Valeur des contraintes "ISOTROPES"
C 
      SUBROUTINE ISOSIG (COMISO, EPSILO, SIGMA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  COMISO(6), EPSILO(6), SIGMA(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ISOSIG')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPTDP('VALEUR DE EPSILON ISO ', EPSILO(1) , 1 , 6 )
CD    CALL IMPTDP('VALEUR DU COMPORT ISO ', COMISO(1) , 1 , 6 )
C 
      SIGMA(1) = COMISO(1)*EPSILO(1)+COMISO(2)*EPSILO(2)
     &           +COMISO(4)*EPSILO(6)
      SIGMA(2) = COMISO(2)*EPSILO(1)+COMISO(1)*EPSILO(2)
     &           +COMISO(4)*EPSILO(6)
      SIGMA(3) = COMISO(3)*EPSILO(3)
      SIGMA(4) = COMISO(5)*EPSILO(4)
      SIGMA(5) = COMISO(5)*EPSILO(5)
      SIGMA(6) = COMISO(4)*(EPSILO(1)+EPSILO(2))+COMISO(6)*EPSILO(6)
C 
CD    CALL IMPTDN('VALEUR DE SIGMA ISO ', SIGMA(1) , 1 , 6 )
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT  Nombre de fonctions du temps
C     E ...... TOUEPS  Deformations (isotropes)
C                      (neps, nteta, ngau1, nfonct)
C 
C     Et on recupere :
C 
C     S ...... SIGADM  Contraintes (isotropes)
C                      (neps, nteta, ngau1, nfonct)

      SUBROUTINE TOUSIG (NFONCT, TOUEPS, SIGADM)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER            NFONCT
      DOUBLE PRECISION   TOUEPS(NEPS*NTETA*NGAU1*NFONCT)
      DOUBLE PRECISION   SIGADM(NTETA*NEPS*NGAU1*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     ADCOM , NEPSCO , NUCOU , DEBEPS  , EPS , NFT
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TOUSIG')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     ADCOM est l'adresse de depart du comportement isotrope transverse
C     de la couche
C     NEPSCO est le nombre de deformation dans une couche (teta
C     compris)
C 
      NEPSCO = NTETA*XINTEG*YINTEG*NBCOL
      DEBEPS  = 1
C 
      DO NFT = 1 , NFONCT
C 
CD      CALL IMPEN (' POUR LA FONCTION DU TEMPS NUMERO', NFT )
C 
        DO NUCOU = 1 , NBCOU
          CALL COMCOU( NUCOU , ADCOM )
          DO EPS = 1 , NEPSCO
            CALL ISOSIG (DM(ADCOM), TOUEPS(DEBEPS), SIGADM(DEBEPS))
            DEBEPS = DEBEPS +6
          END DO
        END DO
      END DO
C 
CD    DO NFT = 1 , NFONCT
CD      CALL IMPEN (' POUR LA FONCTION DU TEMPS NUMERO', NFT )
CD      CALL OMPTDN(' CONTRAINTES ADMISSIBLES RANGEES 6*NTETA*NGAU1 '
CD         , SIGADM(1+(NFT-1)*6*NTETA*NGAU1) , 6 * NTETA*NGAU1 , 1)
CD    END DO
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
C           L = 1 , NL
C                                           _
C             TLKJI ( l , k , j , i ) = TR( Wk Kij Vl )
C 
C     La matrice KT etant symetrique, on ne fait varier J que de
C     1 a I, ce qui correpond au rangement de KIJPG.
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions du temps
C     E ...... MULT     le multiplicateur d'integration du point de Gauss
C     E ...... KIJPG    le tableau des matrices au point de
C                       Gauss range (17,[ i , j])
C     E ...... Nk       le nombre de fonctions WK
C     E ...... WK       le 1er tableau de deformations range (NEPS,NFONCT)
C     E ...... NL       le nombre de fonctions VL
C     E ...... VL       le 2eme tableau de deformations range (NEPS,NFONCT)
C 
C     Et on recupere :
C 
C     ES...... TLKJIE   le tableau des traces au point de gauss
C                       en entree range comme indique plus haut
C     ES...... TLKJIS   le tableau des traces au point de gauss
C                       en sortie range comme indique plus haut
C                       => + TLKJIE
C 
C     REM : inutile de mettre TLKJI a 0. On peut faire TLKJIE = TLKJIS !!!!
C 
      SUBROUTINE SCBFPG (NFONCT, MULT, KIJPG, NK,
     &                   WK, NL, VL, TLKJIE, TLKJIS)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NFONCT , NK , NL
      DOUBLE PRECISION  KIJPG(17*NFONCT*(NFONCT+1)/2)
      DOUBLE PRECISION  WK (NEPS*NK) , MULT
      DOUBLE PRECISION  VL (NEPS*NL)
      DOUBLE PRECISION  TLKJIE( NL , NK , NFONCT , NFONCT )
      DOUBLE PRECISION  TLKJIS( NL , NK , NFONCT , NFONCT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  TRACE
      INTEGER           I , J , K , L , DEBKIJ , DEPSK  , DEPSL
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SCBFPG')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     DEBUT DE LA MATRICE KIJ DANS KIJPG
C 
      DEBKIJ = 1
C 
C     DEBSK DEBUT DE LA DEFORMATION WK
C     DEBSL DEBUT DE LA DEFORMATION VL
C 
      DO  I = 1, NFONCT
        DO   J = 1 , I
          DEPSK = 1
          DO   K = 1 , NK
            DEPSL = 1
            DO   L = 1 , NL
              CALL TRAC12 (KIJPG(DEBKIJ), VL(DEPSL), WK(DEPSK),
     &                     TRACE)
              TLKJIS (L, K, J, I) = TLKJIE (L, K, J, I)
     &                                  +TRACE*MULT
              DEPSL = DEPSL+NEPS
            END DO
            DEPSK = DEPSK+NEPS
          END DO
          DEBKIJ = DEBKIJ+17
        END DO
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
CD      CALL OMPTDN( ' TLKJIS SOUS LA FORME [ (NL*NK),N2 ]',
CD                     TLKJIS(1,1,1,1) , NL*NK,NFONCT*NFONCT )
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
      INTEGER FUNCTION BIJNNK( N ,NK , P , Q  )
C
C - permet de faire la correspondance entre la notation tensorielle
C - et matricielle  : bijection de [ 1, NK ] X [ 1, N ] => [ 1 , NK*N ]
C -                                (   P     ,    Q  ) =>     K
C On envoie comme arguments:
C E................ N        la dimension
C E................ P        1er element du couple
C E................ Q        2eme element du couple
C E                          (NEPS,NFONCT)
C Et on recupere:  BIJN2 ( N , P , Q ) = K tq :
C S................ K        ( Q-1)*NK +P
C -
C*      DECLARATION DES PARAMETRES GLOBAUX
C       """"""""""""""""""""""""""""""""""
      INTEGER N , P , Q , NK
C
      BIJNNK = (Q-1)*NK+P
C
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule la valeur de TLKJI comme integrale en teta
C     de la fonction FIJTET rangee (NTETA, NFONCT**4).
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT le nombre de fonctions
C 
C         attention : normalement tableau a 5 indices !!
C 
C     E ...... FIJTET les fonctions rangees (NFONCT**4, NTETA)
C 
C     Et on recupere:
C 
C         attention : normalement tableau a 4 indices !!
C 
C     S ...... TLKJI les fonctions integrees rangees (NFONCT**4)
C 
      SUBROUTINE ITLKJI (NFONCT, NK, NL, FIJTET, TLKJI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NFONCT  , NK , NL
      DOUBLE PRECISION  FIJTET( NFONCT*NFONCT*NK*NL* NTETA )
      DOUBLE PRECISION  TLKJI ( NFONCT*NFONCT*NK*NL )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         LONFON  , TABNIV(2)  , I
      INTEGER         AM2LC   , ADM2LC , VECT , LONG
C 
CD    LOGICAL         LTRACN  , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ITLKJI')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     SEQUENCE D'IMPRESSION EN TRACE PROFONDE
C 
CD    IF (LTRACP(1))THEN
CD      CALL OMPTDP( ' FIJTET ',
CD                     FIJTET(1) , NFONCT*NFONCT*NK*NL ,NTETA )
CD    END IF
C 
      LONFON    = NFONCT*NFONCT*NK*NL
      TABNIV(1) = LONFON
      TABNIV(2) = NTETA
C 
      CALL POUSMD( NTETA , VECT )
C 
      DO I = 1 , LONFON
C 
CD       CALL IMPEP (' TERME NUMERO ', I)
C 
        CALL EXTRAD (FIJTET, 2, TABNIV, 2, I, DM(VECT), LONG)
C 
CD       CALL OMPTDP(' VECTEUR EXTRAIT ', DM(VECT), LONG, 1)
C 
        CALL INTETA( DM(VECT) , TLKJI(I) )
C 
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
CD      CALL OMPTDN( ' TLKJI SOUS LA FORME (NL*NK , N2 ) ',
CD                     TLKJI(1) , NL*NK ,NFONCT*NFONCT )
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
C     Cette routine calcule la trace au point de Gauss des NFONCT
C     fonctions de l'espace WK par les NFONCT fonction  VL.
C     Le resultat de cette routine est le tableau TLKJI a 4 indices
C     tel que, en rangement informatique :
C 
C     I = 1 , NFONCT
C       J = 1 , I
C         K = 1 , NK
C           L = 1 , K
C                                           _
C             TLKJI ( l , k , j , i ) = TR( Wk Kij Vl )
C 
C     La matrice KT etant symetrique on ne fait varier J que de
C     1 a I  ce qui correpond au rangement de KIJPG.
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT   le nombre de fonctions
C     E ...... MULT     le multiplicateur d'integration
C     E ...... KIJPG    le tableau des matrices au point de
C                       Gauss range (17, [ i , j]
C     E ...... WK       le 1er tableau de deformations range (NEPS, NFONCT)
C 
C                       ATTENTION WK = VL
C 
C     E....... VL       le 2eme tableau de deformations range (NEPS, NFONCT)
C 
C     Et on recupere :
C 
C     ES...... TLKJIE   le tableau des traces au point de gauss
C     ES                en entre ranges comme indique plus haut
C     ES...... TLKJIS   le tableau des traces au point de gauss
C     ES                en sortie ranges comme indique plus haut
C 
C     REM : inutile de mettre TLKJI a 0.
C     REM : on peut faire TLKJIe = TLKJIS
C 
      SUBROUTINE NOBFPG (NFONCT, MULT, KIJPG, NK, WK, NL, VL,
     &                   TLKJIE, TLKJIS)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NFONCT ,  NK , NL
      DOUBLE PRECISION  KIJPG(17*NFONCT*(NFONCT+1)/2)
      DOUBLE PRECISION  WK (NEPS*NK)
      DOUBLE PRECISION  VL (NEPS*NL)  , MULT
      DOUBLE PRECISION  TLKJIE( NL , NK , NFONCT , NFONCT )
      DOUBLE PRECISION  TLKJIS( NL , NK , NFONCT , NFONCT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION TRACE
      INTEGER          I , J , K , L , DEBKIJ , DEPSK  , DEPSL
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NOBFPG')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF ( NL . NE . NK  ) THEN
        CALL ERREUD( 0, ' NL diff  de NK mauvaise utilisation de '
     &                //IDPROG)
      END IF
C 
C     SEQUENCE D'IMPRESSION EN TRACE PROFONDE
C 
CD    IF (LTRACP(1))THEN
CD      CALL OMPTDP( ' WK en entree ',WK(1) , NEPS ,NK )
CD      CALL OMPTDP( ' CMTGIJ EN ENTREE ',KIJPG(1)
CD                 ,17,NFONCT*(NFONCT+1)/2)
CD      CALL OMPTDP( ' TLKJIE SOUS LA FORME (NK*NK , N2 )',
CD                     TLKJIE(1,1,1,1) , NK*NK ,NFONCT*NFONCT )
CD    END IF
C 
C     Debut de la matrice KIJ dans KIJPG
C 
      DEBKIJ = 1
C 
C     DEBSK debut de la deformation WK
C     DEBSL debut de la deformation VL
C 
      DO  I = 1, NFONCT
        DO   J = 1 , I
          DEPSK = 1
          DO   K = 1 , NK
            DEPSL = 1
            DO   L = 1 , K
              CALL TRAC12( KIJPG( DEBKIJ) , VL( DEPSL) , WK(DEPSK) ,
     &                     TRACE )
              TLKJIS( L , K , J , I ) = TLKJIE( L , K , J , I )
     &                                  + MULT*TRACE
              DEPSL = DEPSL+NEPS
            END DO
            DEPSK = DEPSK+NEPS
          END DO
          DEBKIJ = DEBKIJ+17
        END DO
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
CD      CALL OMPTDN( ' TLKJIS SOUS LA FORME (NK*NK , N2 )',
CD                     TLKJIS(1,1,1,1) , NK*NK ,NFONCT*NFONCT )
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine exploite les symetries du tableau TIJKL dans
C     le cas de la routine NOKGLO.
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT  le nombre de fonctions
C     E ...... NK      le nombre de vecteurs WK
C     ES...... TLKJI   RANGEE (NK, NK, NFONCT, NFONCT)

      SUBROUTINE SYMENO (NFONCT, NK, TLKJI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER             NFONCT , NK
      DOUBLE PRECISION    TLKJI( NK , NK ,NFONCT, NFONCT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER    I , J , K , L
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SYMENO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACP(1))THEN
CD     CALL OMPTDP( 'TLKJI en  entree NK*NL , N2  ', TLKJI(1,1,1,1) ,
CD                  NK*NK , NFONCT*NFONCT )
CD    END IF
C 
C     Symetrie par rapport a la diagonale de chaque bloc;
C     (I > OU = J ), ( L > K )
C 
      DO   I = 1 , NFONCT
        DO   J = 1 , I
          DO   K = 1 , NK
             DO  L = K+1 , NK
                TLKJI( L , K , J , I ) = TLKJI( K , L , J , I )
             END DO
          END DO
        END DO
      END DO
C 
CD    IF (LTRACP(1))THEN
CD      CALL OMPTDP( 'TLKJI apres symetrie diagonale  (NK*NK , N2)  ',
CD                   TLKJI(1,1,1,1) , NK*NK , NFONCT*NFONCT )
CD    END IF
C 
C     Symetrie des blocs non diagonaux
C     ( I < J ), ( L , K  quelconques )
C 
      DO   I = 1 , NFONCT
        DO   J = I+1  , NFONCT
          DO   K = 1 , NK
             DO  L = 1 , NK
                TLKJI( L , K , J , I ) = TLKJI( L , K , I , J )
             END DO
          END DO
        END DO
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE

CD    IF (LTRACN(1))THEN
CD      CALL OMPTDN( 'TLKJI en  sortie (NK*NK , N2)  ', TLKJI(1,1,1,1) ,
CD                  NK*NK , NFONCT*NFONCT )
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine exploite les symetries du tableau TIJKL dans
C     le cas de la routine SCKGLO.
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT  le nombre de fonctions
C     E ...... NK      le nombre de vecteurs WK
C     E ...... NL      le nombre de vecteurs VL
C     ES...... TLKJI   RANGEE (NL, NK, NFONCT, NFONCT)

      SUBROUTINE SYMESC (NFONCT, NK, NL, TLKJI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER             NFONCT , NK , NL
      DOUBLE PRECISION    TLKJI( NL , NL ,NFONCT, NFONCT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I , J , K , L
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SYMESC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACP(1))THEN
CD      CALL OMPTDP( 'TLKJI en  entree ( NK*NL , N2)  ',TLKJI(1,1,1,1),
CD                  NK*NL , NFONCT*NFONCT )
CD    END IF
C 
C     Symetrie des blocs non diagonaux
C     (I < J), (L, K quelconques)
C 
      DO   I = 1, NFONCT
        DO   J = I+1  , NFONCT
          DO   K = 1 , NK
             DO  L = 1 , NL
                TLKJI (L, K, J, I) = TLKJI (L, K, I, J)
             END DO
          END DO
        END DO
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
CD      CALL OMPTDN( 'TLKJI en  sortie (NL*NK , N2 ) ', TLKJI(1,1,1,1) ,
CD                  NL*NK , NFONCT*NFONCT )
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
