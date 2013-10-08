C     Cette routine calcule le profil des matrices K0n stockees profil
C 
      SUBROUTINE PROFIL
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
      INTEGER T1, TR1, T2, TR2, A1
      INTEGER R1, AP1, AR1, IE, K, J
      INTEGER X, Y, NUMIN, NUACT, I
      INTEGER MAX0, NNN, LARBAN
      INTEGER AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='PROFIL')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     RECHERCHE DE LA PREMIERE ADRESSE DANS TABELE
C 
      CALL ADTBM ('TABELE     ', T1)
C 
C     POUR LES ELEMENTS DE TYPE 1
C 
CD    CALL IMPEP ('1ERE ADRESSE DE TABELE', T1)
C 
      TR1=T1-1
C 
C     POUR LES ELEMENTS DE TYPE 2
C 
      T2=T1+9*NEL1
      TR2=T2-1
C 
C -----------------------------------------------------------------------
C     Ouverture d'un tableau de nom PRODL, de 1ere adresse A1
C     pour ranger le profil par ddl
C -----------------------------------------------------------------------
      CALL GESTEN ('PRODL     ', NDDL+1, A1)
      CALL MENAM (A1, NDDL+1)
      R1=A1-1
C 
C -----------------------------------------------------------------------
C     Ouverture d'un tableau partiel de 1ere adresse AP1
C     pour ranger le profil par noeud
C -----------------------------------------------------------------------
      CALL GSPOUE (NNOEUD, AP1)
      AR1=AP1-1
C 
C -----------------------------------------------------------------------
C     Calcul du profil par noeud
C -----------------------------------------------------------------------
C 
C     CONTRIBUTION DES ELEMENTS DE COUCHE
C 
      DO IE=1, NEL1
C 
C     X est la position du plus petit numero de noeud de l'element
C 
        X=T1-6+9*IE
CD      CALL IMPEP('X',X)
        NUMIN=M(X)
CD      CALL IMPEP('NUMIN',NUMIN)
        DO I=0,5
          NUACT=M(X+I)
          Y=AR1+NUACT
CD        CALL IMPEP('VALEUR DE Y',Y)
          M(Y)=MAX0(M(Y),NUACT-NUMIN)
CD        CALL IMPEP('VALEUR DE M(Y)',M(Y))
        ENDDO
      ENDDO
C 
C     CONTRIBUTION DES ELEMENTS D'INTERFACE
C 
      IF (SYMPAR) THEN
        DO IE = 1, NEL2
C 
C     X est la position du plus petit numero de noeud de l'element
C 
          X=T2-4+7*IE
          NUMIN=M(X)
          IF ( NUMIN .LE. 0 ) NUMIN=M(X+2)
          DO I=0,3
            NUACT=M(X+I)
            IF ( NUACT .GT. 0 ) THEN
              Y=AR1+NUACT
              M(Y)=MAX0(M(Y),NUACT-NUMIN)
            END IF
          ENDDO
        ENDDO
      ELSE
        DO IE=1,NEL2
C 
C     X est la position du plus petit numero de noeud de l'element
C 
          X=T2-4+7*IE
          NUMIN=M(X)
          DO I=0,3
            NUACT=M(X+I)
            Y=AR1+NUACT
            M(Y)=MAX0(M(Y),NUACT-NUMIN)
          ENDDO
        ENDDO
      ENDIF
CD    CALL IMPTET('PROFIL PAR NOEUD',M(AP1),1,NNOEUD)
C 
C -----------------------------------------------------------------------
C    Modification du profil en profil par ddl
C -----------------------------------------------------------------------
      K=1
      DO I=1,NNOEUD
        NNN=6*M(AR1+I)
          DO J=1,6
            K=K+1
            M(R1+K)=NNN+J+M(R1+K-1)
          ENDDO
      ENDDO
      NTMAT=M(R1+K)
      CALL IMPET ('NOMBRE DE TERMES PAR MATRICE ', NTMAT)
C 
C -----------------------------------------------------------------------
C     POUR DIMENSIONNER LES FICHIERS
C 
C     impression de NTMAT
C     impression du profil par ddl
C -----------------------------------------------------------------------
CD    CALL IMPTET ('PROFIL PAR DDL ', M(A1), 1, NDDL+1)
C   
C -----------------------------------------------------------------------
C     CALCUL DE LA DEMI-LARGEUR DE BANDE
C -----------------------------------------------------------------------
      DO I = 1 , NDDL+1
         LARBAN = MAX0 (M(A1+I)-M(A1+I-1), LARBAN)
      END DO
      CALL IMPET ('DEMI-LARGEUR DE BANDE ', LARBAN)
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
C     Cette routine ouvre les tableaux de rangement des matrices
C     KOn ainsi que le tableau des adresses de depart de ces
C     tableaux de nom ADRESSES.
C 
C     Cette routine ouvre aussi les tableaux de rangement des
C     efforts (MAT-EFFORT).
C 
C     Cette routine ouvre aussi les tableaux de rangement des
C     deplacements(MAT-DEPLA).
C 
C     Les efforts seront ranges par valeurs de N croissantes
C     de -NTDSFG a NTDSFG

      SUBROUTINE OUVTAB
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
      INTEGER  M1, I, NMAT, R1, LONEFF
      INTEGER  NTERME, LONSAD
      INTEGER  ADEFF, ADDEPL, ADEPS, ADSIG, ADSAU, ADSIGN
C 
C     Pour la solution isotrope
C 
      INTEGER LONEPS
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='OUVTAB')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      LONEFF = NBMAT*NDDL*NBFODO
      NMAT   = NBMAT
C 
C -----------------------------------------------------------------------
C     Adresses du tableau des efforts initiaux
C -----------------------------------------------------------------------
C 
      CALL GESTDP ('MAT-EFFORT', LONEFF, ADEFF)
      CALL MENADM (ADEFF, LONEFF)
      CALL GESTEN ('ADR-EFFORT', NBMAT, M1)
      R1    = M1-1
      NTERME = ADEFF
      DO I = 1, NBFODO
        M(R1+I)  = NTERME
        NTERME   = NTERME+NDDL*NBMAT
      ENDDO
C 
C -----------------------------------------------------------------------
C     Adresses des tableaux admissibles a l'iteration
C -----------------------------------------------------------------------
C 
      CALL GESTDP ('DEPLA-ADMI', CHARAX*NTETA*NDDL, ADDEPL)
      CALL GESTDP ('MAT-DEPLAC', LONEFF, ADDEPL)
C  
      LONEPS = NEPS*NTETA*NGAU1*CHARAX
      LONEFF = NBMAT*NDDL*CHARAX
      LONSAD = NTETA*NSAU*NBINT*XINTEG*NBCOL*CHARAX
C 
C -----------------------------------------------------------------------
C     Adresses des tableaux des deformations et contraintes admissibles
C -----------------------------------------------------------------------
C 
      CALL GESTDP ('EPS-AD-TOT', LONEPS, ADEPS)
      CALL GESTDP ('SIG-AD-TOT', LONEPS, ADSIG)
C 
C     Pour les interfaces
C 
      IF (NBINT .GT. 0) THEN
        CALL GESTDP ('SAU-AD-TOT', LONSAD, ADSAU)
	CALL GESTDP ('SGN-AD-TOT', LONSAD, ADSIGN)
      ENDIF
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C  
C     Cette routine calcule la place necessaire a l'execution
C 
      SUBROUTINE DETAIL 
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
      INTEGER TAILDM, TAISTM, TAIMAT, TAIMIN, LECINT
      INTEGER IREST
      LOGICAL LECLOG
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DETAIL')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL LECSEQ ('DETAIL', 'DECLARATION DES PARAMETRES DE TAILLE ')
C 
C     Nombre de champs admissibles a zero maxi crees
C 
      CHAMAX = LECINT ('CHAMAX')
C 
C     Nombre de champs admissibles maxi crees
C 
      CHARAX = LECINT ('CHARAX')
C 
C     SYM indique si la solution est symetrique par rapport au plan moyen
C     
      SYM = LECLOG ('SYM')
C 
C     Donnees maillage :
C     NOMBRE-DE-COUCHES total meme en cas de symetrie
C     NOMBRE-DE-COLONNES est le nombre d'elements finis dans une couche
C     NTETA est le nombre NTETA d'angles qui vaut 2x(le nombre de termes
C     de la serie de fourier )
C 
      NBCOU= LECINT ('NOMBRE-DE-COUCHES')
      NBCOL  = LECINT ('NOMBRE-DE-COLONNES')
      NTETA  = LECINT ('NTETA')
C 
C     Ccalcul du nombre de matrice NBMAT
C 
      NTDSFG = NTETA/2
      NBMAT = 2*NTDSFG+1
C 
      NPICET = LECINT ('NB-PAS-TEMPS')
C 
      XINTEG=LECINT ('XINTEG')
      YINTEG=LECINT ('YINTEG')
C 
C     Donnees a-priori : nombres de composantes en contraintes, de composantes
C     en deformations, de composantes en sauts, de ddl par element
C 
      NSIG    = 6
      NEPS    = 6
      NSAU    = 3
      NDDLEL  = 36
C 
C     SYMPAR : logique indiquant si le plan de symetrie est une interface
C 
      SYMPAR = .FALSE.
      NBINT    = NBCOU-1
C 
      IF (SYM) THEN
        IREST = MOD(NBCOU,2)
        NBCOU = NBCOU/2+IREST
C 
        CALL IMPET (IDPROG// 'SYMETRIE => NOMBRE DE COUCHES ', NBCOU)
        NBINT  = NBCOU
C 
        SYMPAR = .TRUE.
        IF (IREST .EQ. 1) THEN
          NBINT   = NBINT-1
          SYMPAR = .FALSE.
        END IF
C 
        CALL IMPET (IDPROG// 'SYMETRIE => NOMBRE D''INTERFACES ', NBINT)
C 
      END IF
C 
C     DONNEES DEPENDANTES DE LA GEOMETRIE
C 
      NNOEUD = 3*NBCOU*(NBCOL+1)
      NDDL    = 6*NNOEUD
      NEL1    = NBCOU*NBCOL
      NEL2    = NBINT*NBCOL
      NBANGL  = NBINT+NBCOU
C 
      CALL IMPET (IDPROG// ' NOMBRE DE DDL PAR BANDE                   '
     &            //'                 ', NDDL)
C 
      CALL IMPET (IDPROG// ' NOMBRE DE DDL TOTAL                       '
     &            //'                 ', NDDL*NTETA)
C 
C 
C     CALCUL DU NOMBRE DE POINTS DE GAUSS TOTAL DANS LES COUCHES
C 
      NGAU1  = XINTEG*YINTEG*NBCOU*NBCOL
C 
C     CALCUL DU NOMBRE DE POINTS DE GAUSS TOTAL DANS LES INTERFACES
C 
      NGAU2  = XINTEG*NBINT*NBCOL
C 
      TAILDM = 2*CHAMAX*XINTEG*NBCOL*NTETA*(6*YINTEG*NBCOU+3*NBINT)
C 
      CALL IMPET (IDPROG// ' LDM NE DOIT PAS ETRE INFERIEUR A          '
     &            //'               ', TAILDM)
C 
      TAISTM =   12*NTETA*NPICET*(2*NGAU1+NGAU2)
C 
      CALL IMPET (IDPROG// ' QUANTITES STOCKEES A L''ETAPE LOCALE      '
     &            //'               ', TAISTM)
C 
      IF (2*NBCOL .LE. (3*NBCOU-1)) THEN
C 
        TAIMAT =  NBMAT* NDDL*12*(NBCOL+2)
C 
        CALL IMPET (IDPROG// ' LARGEUR DE BANDE                        '
     &              //'               ',12*(NBCOL+2))
C 
      ELSE
C 
        TAIMAT =  NBMAT* NDDL*18*(NBCOU+1)
C 
        CALL IMPET (IDPROG// ' LARGEUR DE BANDE                        '         
     &	            //'               ',18*(NBCOU+1))
C 
      END IF
C 
      CALL IMPET (IDPROG// ' TAILLE ESTIMEE DES MATRICES AU MOYEN D''UN'
     &            //' STOCKAGE BANDE ', TAIMAT)
C 
      TAIMIN = 2*TAILDM+TAISTM+TAIMAT
C 
      CALL IMPET (IDPROG// ' PLACE MINI SUR LE DISQUE POUR TOURNER ',
     &            TAIMIN)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
