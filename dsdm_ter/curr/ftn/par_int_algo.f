C     On envoie comme arguments :
C 
C     E ...... TABGAU  Tableaux N pour le point de
C                      gauss concerne
C     E ...... A       La demi-longeur de l'element
C     E ...... VALDEP  La valeur des deplacements sur l'element
C                      !!!!! ranges calcul ===== > PASSAGE PAR RAN-INV
C                      (VALDEP( 3*8 ))
C 
C     Et on recupere :
C 
C     S ...... SAUTN  Valeur de la composante N des sauts rangee :
C                           (U, V, W)
C                           attention  !!!!!!!! LE DEVELOPPEMENT EST DU
C                           type suivant:
C                           N > OU = 0  ==> ( COSN , SINN , COSN)
C                           N <      0  ==> ( SINN , COSN , SINN)
C 
      SUBROUTINE SAUDEV (TABGAU, A,  VALDEP, SAUTN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION   TABGAU(8) , A , VALDEP(24) , SAUTN(3)
      INTEGER AD2
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SAUDEV')
      INTEGER     AM2LC      , ADM2LC
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------      
CD    CALL OMPTDP('DEPDEV POUR L''INTERFACE', VALDEP(1),8,3)
CD    CALL IMPDP('DEMIE LONGUEUR ', A)
C 
      CALL GSPOUD( 8 , AD2 )
      CALL LINLOC( A , 8 , TABGAU(1 ) , DM(AD2)   )
C 
      CALL MATM02 (DM(AD2), VALDEP(1), SAUTN(1), 1, 8, 1)
      CALL MATM02 (DM(AD2), VALDEP(9), SAUTN(2), 1, 8, 1)
      CALL MATM02 (DM(AD2), VALDEP(17), SAUTN(3), 1, 8, 1)
C 
C     TCD CALL PSCAL (DM(AD2), VALDEP(1), SAUTN(1), 8)
C     TCD CALL PSCAL (DM(AD2), VALDEP(9), SAUTN(2), 8)
C     TCD CALL PSCAL (DM(AD2), VALDEP(17), SAUTN(3), 8)
C 
CD    CALL IMPTDN('VALEUR DU SAUT DE U , V , W ', SAUTN(1) , 1 , 3 )
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
         SUBROUTINE SADVCA( TOUDEP , RANINT , NDDLU , NDDLV , NDDLW
     &                      , DEPCAN )
C On envoie comme arguments:
C E................  TOUDEP   Tableau MAT-DEPLA (DM(ADDEPA))
C E................  RANINT   Tableau RANINV-INT (M(ARINCO))
C E................  NDDLU    Numero des deplacements de U ranges
C                             Croissant
C E................  NDDLV    Numero des deplacements de V ranges
C                             Croissant
C E................  NDDLW    Numero des deplacements de W ranges
C                             Croissant
C
C Et on recupere:
C S................  DEPCAN   Tableau des valeurs des deplacements pour
C                             l'element ranges calcul par ordre
C                             croissant de numero de developpement
C                             c-a-D -NTDSFG ------- > NTDSFG
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
C
      include 'cominc.h'
C
      DOUBLE PRECISION  TOUDEP(NBMAT*NDDL) , DEPCAN(24*NBMAT)
      INTEGER           RANINT(8) , NDDLU(8) , NDDLV(8) , NDDLW(8)
C
C**********************************************************************
C *     Declaration des parametres locaux
C  -
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SADVCA')
C
      INTEGER     DECALN , K , U , V , W  , NUDEV
C  -
CD    CALL WLKBCD(IDPROG)
C
C***********************************************************************
      DECALN = 0
      K      = 1
CD    CALL IMPEP(' NDDL DANS '//IDPROG,NDDL)
      DO NUDEV  =  1 , NBMAT
CD       CALL IMPEP(' POUR NUDEV ',NUDEV-NTDSFG-1)
         DO U =  1 , 8
C                CALL IMPEP('POUR U ', U)
C                CALL IMPEP('NDDLU(RANINT(U)) ' ,NDDLU(RANINT(U)))
                DEPCAN(K)   = TOUDEP( DECALN+NDDLU ( RANINT( U ) ))
                K           = K+1
         ENDDO
         DO V =  1 , 8
C                CALL IMPEP('POUR V ', V)
C                CALL IMPEP('NDDLV(RANINT(V)) ' ,NDDLV(RANINT(V)))
                DEPCAN(K)   = TOUDEP( DECALN+NDDLV ( RANINT( V ) ))
                K           = K+1
         ENDDO
         DO W =  1 , 8
C                CALL IMPEP('POWR W ', W)
C                CALL IMPEP('NDDLW(RANINT(W)) ' ,NDDLW(RANINT(W)))
                DEPCAN(K)   = TOUDEP( DECALN+NDDLW ( RANINT( W ) ))
                K           = K+1
         ENDDO
         DECALN = DECALN+NDDL
      ENDDO
CD    CALL OMPTDN ('DEPLACEMENTS RANGES CALCULS PAR DEVELOPPEMENT ',
CD                  DEPCAN(1) , 24 , NBMAT )
C***********************************************************************
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Pour une fonction du temps :
C     Calcule les sauts admissibles pour tout les points de gauss
C     ils sont ranges (nsau, nteta, ngaus2)
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT Nombre de fonctions du temps
C     E ...... DEPDEV Tableau des deplacements developpes
C 
C     Et on recupere :
C 
C     S ...... SAUADM Tableau des deformations reelles (admissibles)
C                     (nsau, nteta, ngau2, nfonct)
C 
      SUBROUTINE TOUSAU (NFONCT, DEPDEV, SAUADM)
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
      DOUBLE PRECISION   SAUADM(NTETA*NSAU*NGAU2*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  RANINV, ADRGAU, LONGAU, RANINP
      INTEGER  ADM2DP, ADM2EP, ADM2VE, NUINT, NUCOL, TLOCN1
      INTEGER  ADDLU, ADDLV, ADDLW, DEBGAU, X, NUDEV
      INTEGER  DBDEPN, DBSAUN, DECAL, DBPREE, NFT
C 
      DOUBLE PRECISION   A, RAYONC
C 
      LOGICAL PASSP
C 
C     POUR TEST
C 
      INTEGER     AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TOUSAU')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     RECHERCHE DES DIFFERENTS TABLEAUX UTILES
C 
      CALL ADTBM ('RANINV-INT', RANINV)
      CALL ADTBM ('RANINV-INP', RANINP)
      CALL ADTBM ('TLOCN1    ', TLOCN1)
      CALL ADTBDM ('GAUSSINTER', DEBGAU)
C 
C     Ouverture d'un tableau partiel pour ranger les numeros
C     des ddl des deplacements calcules (u, v, w) par developpement croissant
C 
      CALL GSPOUE (24, ADDLU)
      ADDLV = ADDLU+8
      ADDLW = ADDLV+8
C 
C     Ouverture d'un tableau partiel pour ranger les deplacements
C     croissant (u, v, w) par developpement croissant, puis les deformations
C 
      CALL GSPOUD( NBMAT*(24+NSAU+1),ADM2DP )
      ADM2EP = ADM2DP+NBMAT*24
      ADM2VE = ADM2EP+NBMAT*NSAU
C 
C     Decalage d'indice pour epsree
C 
      DECAL       = NSAU*NTETA
      DBPREE      = 1
C 
      DO NFT = 1, NFONCT
C 
        DO NUINT = 1, NBINT
C 
        PASSP = .FALSE.
        IF (SYMPAR .AND. NUINT .EQ. 1) PASSP = .TRUE.
C 
          DO NUCOL = 1,NBCOL
C 
            CALL VALRAY( NUCOL , RAYONC , A )
C 
C     RECHERCHE  DES NUMEROS DE DDL RANGES CROISSANT POUR L'ELEMENT
C 
            CALL NDDLIU (NUINT, NUCOL, TLOCN1, M(ADDLU), M(ADDLV),
     &                   M(ADDLW) )
C 
C     REMPLISSAGE DE  DEPCAN : Tableau des valeurs des deplacements pour
C     l'element ranges calcul(U , V , W)  par ordre croissant de numero de
C     developpement c-a-d :  - NTDSFG ------- > NTDSFG
C 
            IF (PASSP) THEN
C 
              CALL SADVCP (DEPDEV(1+(NFT-1)*NDDL*NBMAT), M(RANINP)
     &        , M(ADDLU), M(ADDLV), M(ADDLW), DM(ADM2DP))
C 
            ELSE
C 
              CALL SADVCA (DEPDEV(1+(NFT-1)*NDDL*NBMAT), M(RANINV)
     &        , M(ADDLU), M(ADDLV), M(ADDLW), DM(ADM2DP))
C 
            END IF
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
C     Pour le point de gauss considere ( ==> tabgau considere) on
C     calcule la valeur de epsilon( n )
C 
C     DBDEPN est l'adresse pour n du tableau des deplacements provenant
C     de DPDVCA
C 
              DBDEPN = ADM2DP
C 
C     DBSAUN est l'adresse pour n du tableau des deformations
C 
              DBSAUN = ADM2EP
C 
              DO NUDEV = -NTDSFG , NTDSFG
C 
                  CALL SAUDEV( DM(ADRGAU) , A , DM(DBDEPN)
     &                       ,DM(DBSAUN) )
C 
                  DBDEPN = DBDEPN+24
                  DBSAUN = DBSAUN+3
C 
              ENDDO
C 
              ADRGAU = ADRGAU+8
C 
C     On appelle sauree qui pour un point de gauss donne et a partir d'un
C     rangement calcul des deplacements fourni la valeur des deformations
C     rangees (NSAU , nteta , npgau2 )
C 
              CALL SAUREE( DM(ADM2EP) , SAUADM(DBPREE) )
C 
              DBPREE = DBPREE+DECAL
C 
            END DO
C 
          END DO
C 
        END DO
C 
      END DO
C 
      CALL LONGDP('GAUSSINTER',LONGAU )
C 
      CALL TESTEN( LONGAU , ADRGAU-DEBGAU , 'TOUSAU' )
C 
CD    DO NFT = 1, NFONCT
CD      CALL IMPEN (' NUMERO DE FONCTION DU TEMPS ', NFT )
CD      CALL OMPTDN(' SAUTS ADMISSIBLES '
CD       , SAUADM(1+(NFT-1)*NSAU*NTETA*XINTEG*NBINT)
CD       , NSAU , NTETA*XINTEG*NBINT)
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
C     E ...... NFONCT Nombre de fonctions du temps
C     E ...... TOUSAU  sauts reels (nsau, nteta, ngau2, nfonct)
C 
C     Et on recupere :
C 
C     S ...... SGNORM  Contraintes normales (isotropes)
C                     (nsau, nteta, ngau2, nfonct)
C 
         SUBROUTINE TOUSGN (NFONCT, TOUSAU, SGNORM)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      INTEGER            NFONCT
      DOUBLE PRECISION   TOUSAU(NSAU*NTETA*NGAU2*NFONCT)
      DOUBLE PRECISION   SGNORM(NSAU*NTETA*NGAU2*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     ADCOM, NSAUCO, NUINT, DEBSAU
      INTEGER     NSAUT, NFT
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TOUSGN')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     ADCOM est l'adresse de depart du comportement isotrope transverse
C     de la couche
C     NSAUCO est le nombre de sauts par interface (teta compris)
C 
      NSAUCO = NTETA*XINTEG*NBCOL
      DEBSAU  = 1
C 
      DO  NFT = 1 , NFONCT
        DO NUINT = 1 , NBINT
          CALL COMINT( NUINT , ADCOM )
          DO NSAUT = 1 , NSAUCO
            CALL SIGNIS( DM(ADCOM) , TOUSAU(DEBSAU) , SGNORM(DEBSAU) )
            DEBSAU    = DEBSAU +3
          END DO
        END DO
      END DO
C 
CD    DEB = 1
CD    DO NFT = 1 , NFONCT
CD      CALL IMPEN (' FONCTION DU TEMPS NUMERO', NFT )
CD      CALL OMPTDN(' CONTRAINTES ADMISSIBLES RANGEES 3*NTETA*NGAU2 '
CD                   , SGNORM(DEB) , 3 , NTETA*NGAU2 )
CD      DEB = DEB+NSAU*NTETA*NGAU2
C 
CD    END DO
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
C     E ...... COMISO(2) Caracteristiques du comportement isotrope
C                        transverse de la l'interface
C     E ...... SAUT(3)   Valeurs des sauts sr, s0, sz
C 
C     Et on recupere :
C 
C     S ...... SIGMAN(3) Valeur des contraintes normales "ISOTROPES"
C 
      SUBROUTINE SIGNIS (COMISO, SAUT, SIGMAN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  COMISO(2), SAUT(3), SIGMAN(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SIGNIS')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPTDP('VALEUR DE SAUT ISO ', SAUT(1) , 1 , 3 )
CD    CALL IMPTDP(' COMPORT ISO DE L''INTERFACE', COMISO(1) , 1 , 2 )
C 
      SIGMAN(1) = COMISO(1)*SAUT(1)
C 
      SIGMAN(2) = COMISO(1)*SAUT(2)
C 
      SIGMAN(3) = COMISO(2)*SAUT(3)
C 
CD    CALL IMPMN (' n > = 0 ==> cosn , sinn , cosn  ')
CD    CALL IMPMN (' n <   0 ==> sinn , cosn , sinn  ')
CD    CALL IMPMN (' s13 , s23, s33 ')
CD    CALL IMPTDN('VALEUR DE SIGMA ISO ', SIGMAN(1) , 1 , 3 )
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
