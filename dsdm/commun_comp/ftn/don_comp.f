C     QUE FAIT CETTE ROUTINE?:
C     DMATNL est une routine de lecture et d'initialisation
C     des caracteristiques non lineaires .
C 
C     on y rempli des tableaux (donnees de comportement)
C     par type de comportement :
C 
C     Caracteristiques en endommagement : CARAC-ENDO
C       
C     Caracteristiques en plasticite    : CARAC-PLAS
C 
      SUBROUTINE DMATNL
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
      INTEGER         I
      INTEGER         NP1, NTTEI, NARG1, NARG2, RLOC, NP2
      INTEGER         AD3, AC3
      INTEGER         NCOULU, M1 , R1 , NI
      INTEGER         APM2, APR2 , NP , NT
C 
C     Pour SYMPAR=.true.
C 
      LOGICAL          INTSYM, DIM2, DIM3
      INTEGER          ISYM    
      INTEGER          KISYM   
      INTEGER          AM2LC, ADM2LC
      CHARACTER*4      NUMERO
      CHARACTER*6      IDPROG
      DOUBLE PRECISION LECDBL, D2D3
C 
C     pour la rentree du comportement des couches
C 
      DOUBLE PRECISION     E11 , E22 , V12  , V23 , V13 
      DOUBLE PRECISION     G12 , G23 , G13 , E33
      DOUBLE PRECISION     S11 , S22 , S12 , S66 , A1 , A2 , B1 , B2 , C
      INTEGER SOUPLI, ADSOUP, ADELAS
      INTEGER NSOUOR, NELAOR
C 
C     pour la rentree du comportement des interfaces
C 
      DOUBLE PRECISION     S1 , S2 , S3
      DOUBLE PRECISION     E1 , E2 , E3 , G1 , G2 , A
      INTEGER ADIAEN , NSOIOR
      INTEGER ADSOUI
      INTEGER NTCEND , ADCAEN , NNONL                              
C 
C     NOMBRE DE CARACTERISTIQUES NON LINEAIRES PAR COUCHE (15)
C                                                     
      INTEGER NBCNLC, NBCNLI, P, CAENLO, SOUPLO, ELASLO, CAENLI
      DOUBLE PRECISION   EPTLIM , EPCLIM , GAM , R0 , BETA , ALPHA 
      DOUBLE PRECISION   B ,K , Y0 , YC , YR , YCS , N
      DOUBLE PRECISION   PB , PK , PY0 , PYC , YTS , NPP
      DOUBLE PRECISION   SR22T, SR22C, SR12, SR13, SR23
      DOUBLE PRECISION   SRI13, SRI23 
      DOUBLE PRECISION   MINT , ALP
C 
      PARAMETER (IDPROG='DMATNL')
C 
      CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C 
C     rentree du comportement des couches
C 
C     creation d'un tableau pour les angles
C     -de nom ANGLES-COU
C     -il comprend:
C        les angles des couches (de la couche1---->couche nbcou
C        les angles des INTERFACES (
C                   interface 1 adresse (nbcou+1)---->
C                   interface nbcou-1(adresse 2*nbcou-1))
C 
C     NBANGL est la longueur du tableau des angles
C 
        CALL GESTDP ('ANGLES-COU', NBANGL,AD3)
        AC3=AD3-1
C 
         CALL MESSAO (
     &      'Rentree de la sequence d''empilement en donnant les
     &      \angles en degres a partir de la couche inferieure ou
     &      \du milieu si le probleme est symetrique.')
C                                         
         CALL LECLDP (DM(AD3), NBCOU, NCOULU)
C 
         DO I  = 1, NBCOU
           DM(AC3+I) = DM(AC3+I)*PI/180.
         ENDDO
C 
C -----------------------------------------------------------------------
C     rangement dans ANGLES-COU de l'angle des  interfaces
C -----------------------------------------------------------------------
C 
        IF (NBINT.GT.0) THEN
          IF (.NOT. SYMPAR) THEN     
            DO I = 0, NBINT-1
              DM(AD3+NBCOU+I)=(DM(AD3+I+1)+DM(AD3+I))/2.D0      
            ENDDO
          ELSE                                                  
C 
C           Sinon l'interface milieu a un angle de 0
C 
            DM(AD3+NBCOU)= DM(AD3)
            DO I=1,NBINT-1
              DM(AD3+NBCOU+I)=(DM(AC3+I+1)+DM(AC3+I))/2.D0      
            ENDDO
          ENDIF
        ENDIF
C 
C -----------------------------------------------------------------------
C     creation d'un tableau pour ranger les types de comportement
C     concernant les couches(soit NP) on retrouvera la 1ere constante
C     d'elasticite de la couche a l'adresse 10*(NP-1)+1 dans HOO-COUCHE,
C     de nom TYP-COUCHE
C 
C     Steph le 08/09/99. On stocke un operateur de rigidite pour eviter
C     des problemes numeriques poses sur la matrice de souplesse quand
C     l'endommagement atteint 1. -----> 12 termes en plus des 17.
C -----------------------------------------------------------------------
C    
      NARG1 = 17*NBCOU
      CALL GSPOUD (NARG1, ADSOUP)
      NARG2 = 9*NBCOU
      CALL GSPOUD (NARG2, ADELAS)
      CALL GESTEN ('TYP-COUCHE', NBCOU+NBINT, M1)
      R1=M1-1
C 
C     DL du 26/07/96 : Ajout des parametres YCS, PB, PY0, YTS
C     SG du 23/03/98 : Ajout des parametres criteres elastiques 
C     SG du 06/08/98 : Ajout de l'option d2d3 pour ydp(2D) ou ydp(3D)
C    
      NBCNLC =  25
      NTCEND = NBCOU*NBCNLC
C 
      CALL GSPOUD (NTCEND, ADCAEN)
C 
C -----------------------------------------------------------------------
C     creation d'un tableau provisoire pour ranger les numeros
C     des couches concernes par le comportement NP; 1ere adresse : APM2
C -----------------------------------------------------------------------
C 
      CALL GSPOUE (NBCOU, APM2)
      APR2=APM2-1
      CALL MENAM (APM2, NBCOU)
C 
C -----------------------------------------------------------------------
C     tant que toutes les couches n'ont pas de comportement,
C     NP est le nombre de comportements deja rentres
C -----------------------------------------------------------------------
4     NP=0
C 
C     NP1 est le nombre de couches ayant deja un comportement
C 
      NP1=0
C  
      DO WHILE (NP1 .NE. NBCOU)
C  
            CALL MESSAO (
     $     '\donnee de Eij, Gig, vij tels que dans la base 
     $      \d''orthotropie (s = contraintes, e = deformation) :
     $      \   eii = sii/Ei - ( vij/Ej) sjj   (i diff de j)
     $      \   eij = sij / (2Gig)             (i diff de j)')
C  
C  
C     Attention ! les valeurs des Gij sont multiplies par 2 
C     pour avoir le rapport direct entre sij et eij
C  
            CALL LECSEQ ('RIGIDITE',
     $      'MATRICE DE HOOKE DANS LA BASE D''ORTHOTROPIE')
C 
C -----------------------------------------------------------------------
C     NT est la 1ere adresse libre dans SOU-ORTHO
C -----------------------------------------------------------------------
            NSOUOR = ADSOUP+17*NP
	    NELAOR = ADELAS+12*NP
C 
CD          CALL IMPEP('NT=',NT)
C 
C -----------------------------------------------------------------------
C     affectation des valeurs lues 
C -----------------------------------------------------------------------
            E11  = LECDBL ('E11')
            V12  = LECDBL ('V12')
            E22  = LECDBL ('E22')
            G12  = 2.D0*LECDBL ('G12')
            V13  = LECDBL ('V13')
            V23  = LECDBL ('V23')
            G13  = 2.D0*LECDBL ('G13')
            G23  = 2.D0*LECDBL ('G23')
            E33  = LECDBL ('E33')
C 
C ----------------------------------------------------------------------- 
C     Calcul des souplesses 
C ----------------------------------------------------------------------- 
            S11  = 1.D0/E11
            S12  = -V12/E11
            S22  = 1.D0/E22
            S66  = 1.D0/G12
            A1   = -V13/E11
            A2   = -V23/E22
            B2   = 1.D0/G13
            B1   = 1.D0/G23
            C    = 1.D0/E33
C 
            DM(NSOUOR)    = S11
            DM(NSOUOR+1)  = S12
            DM(NSOUOR+2)  = 0.D0
            DM(NSOUOR+3)  = S12
            DM(NSOUOR+4)  = S22
            DM(NSOUOR+5)  = 0.D0
            DM(NSOUOR+6)  = 0.D0
            DM(NSOUOR+7)  = 0.D0
            DM(NSOUOR+8)  = S66
            DM(NSOUOR+9)  = B1
            DM(NSOUOR+10) = 0.D0
            DM(NSOUOR+11) = 0.D0
            DM(NSOUOR+12) = B2
            DM(NSOUOR+13) = A1
            DM(NSOUOR+14) = A2
            DM(NSOUOR+15) = 0.D0
            DM(NSOUOR+16) = C
C 
            DM(NELAOR)    = E11
            DM(NELAOR+1)  = E22
            DM(NELAOR+2)  = E33
            DM(NELAOR+3)  = G12
            DM(NELAOR+4)  = G23
            DM(NELAOR+5)  = G13
            DM(NELAOR+6)  = V12
            DM(NELAOR+7)  = V23
            DM(NELAOR+8)  = V13
C 
C -----------------------------------------------------------------------
C     SG 6/08/98
C     Je rajoute en option d2d3 le choix de yd(3D) ou yd(2D)
C     Il y a une CARAC-NLIN de plus donc 25 par couche
C     On decale pour ranger dans CARAC-NLIN
C -----------------------------------------------------------------------                      
            NNONL = NBCNLC*NP+ADCAEN
C 
C ----------------------------------------------------------------------- 
C     CARACTERISTIQUES DANS LA DIRECTION DES FIBRES
C ----------------------------------------------------------------------- 
            CALL MESSAO (

     $     '\EPTLIM et EPCLIM limites en deformation en traction  
     $      \et compression  dans la direction des fibres, puis 
     $      \donnee de GAM caracterisant le module secant en 
     $      \compression : E11= E110*(1+gam*inf(EPS11, 0))')
C 
            CALL LECSEQ ('LIMSIG ', 'DONNEE DE EPTLIM, EPCLIM, GAM ')
C 
            EPTLIM       =  LECDBL ('EPTLIM')
            EPCLIM       =  LECDBL ('EPCLIM')
            GAM          =  LECDBL ('GAM')
C 
            DM(NNONL)    =  EPTLIM
            DM(NNONL+1)  =  EPCLIM
            DM(NNONL+2)  =  GAM
C  
C -----------------------------------------------------------------------
C     CARACTERISTIQUES DE PLASTICITE
C -----------------------------------------------------------------------
            CALL MESSAO (
     $     '\Donnee de Ro (seuil de plasticite en MPA), de BETA (MPA)
     $      \et de ALPHA tels que : R(p) = Ro + BETA*(p ** ALPHA)  
     $      \(ou p est la plasticite cumulee)
     $      \Donnee de A2, terme de couplage entre 12 et 22')
C  
            CALL LECSEQ ('CPLAST', 'donnee de RO, BETA, ALPHA, A2')
C 
            R0       =  LECDBL ('R0')
            BETA     =  LECDBL ('BETA')
            ALPHA    =  LECDBL ('ALPHA')
            A2       =  LECDBL ('A2')
C 
            DM(NNONL+3)     =  R0
            DM(NNONL+4)     =  BETA
            DM(NNONL+5)     =  ALPHA
            DM(NNONL+6)     =  A2
C 
C -----------------------------------------------------------------------
C    Rentree des caracteristiques d'endommagement des couches
C -----------------------------------------------------------------------
            CALL MESSAO (
     $      '
     $      \                                    2
     $      \  Yd  = 2G120.( epsilon12_elastique)   
     $      \                                         2
     $      \  Ydp = E220.( sup(epsilon22_elastique,0))   
     $      \                        1/2    1/2      1/2
     $      \  Y   = ( ( Yd + b.Ydp )   - Y0   ) / Yc
     $      \  .                          n
     $      \  d   = k. sup[ 0 , (Y - d) ]
     $      \                         1/2     1/2       1/2
     $      \  Yp  = ( ( Ydp + bp.Yd )   - Y0p   ) / Ycp
     $      \  .                             np
     $      \  dp  = kp. sup[ 0 , (Yp - dp) ]     tels que :
     $      \  
     $      \  G12 = G120 (1 - d)  et  E2 = E20 (1 - dp)' )
C  
            CALL LECSEQ ('CENDOM',
     $      ' donnee de b , k ,Y0 , Yc ,n , kp , Ycp , np
     $      \ sr22t, sr22c, sr12, sr13, sr23, d2d3  ')
C 
            B   =  LECDBL ('B')
            K   =  LECDBL ('K')
            Y0  =  LECDBL ('Y0')
            YC  =  LECDBL ('YC')
            YCS =  LECDBL ('YS')
            N   =  LECDBL ('N')
            PB  =  LECDBL ('PB')
            PK  =  LECDBL ('PK')
            PY0 =  LECDBL ('PY0')
            PYC =  LECDBL ('PYC')
            YTS =  LECDBL ('PYS')
            NPP =  LECDBL ('PN')
            SR22T  =  LECDBL ('SR22T')
            SR22C  =  LECDBL ('SR22C')
            SR12 =  LECDBL ('SR12')
            SR13 =  LECDBL ('SR13')
            SR23 =  LECDBL ('SR23')
            D2D3 =  LECDBL ('D2D3')
	    DIM2 = .TRUE.
	    DIM3 = .TRUE.
	    IF (D2D3 .EQ. 2.) DIM2 = .FALSE.
	    IF (D2D3 .EQ. 3.) DIM3 = .FALSE.
C 
            IF (DIM2 .AND. DIM3) THEN
C           
	      CALL MESSAO ('MAUVAISE SAISIE DE L''OPTION 2D OU 3D POUR 
     &                     \LE CALCUL DE L''ENDOMMAGEMENT INTRALAMINAIRE
     &                     \        OPTION 2D  : D2D3 = 2.
     &                     \        OPTION 3D  : D2D3 = 3.')
              GOTO 4
            END IF
C 
            DM( NNONL+7)  =  B
            DM( NNONL+8)  =  K
            DM( NNONL+9)  =  Y0
            DM( NNONL+10) =  YC
            DM( NNONL+11) =  YCS
            DM( NNONL+12) =  N
            DM( NNONL+13) =  PB
            DM( NNONL+14) =  PK
            DM( NNONL+15) =  PY0
            DM( NNONL+16) =  PYC
            DM( NNONL+17) =  YTS
            DM( NNONL+18) =  NPP
            DM( NNONL+19) =  SR22T
            DM( NNONL+20) =  SR22C
            DM( NNONL+21) =  SR12
            DM( NNONL+22) =  SR13
            DM( NNONL+23) =  SR23           
            DM( NNONL+24)  =  D2D3
C 
C -----------------------------------------------------------------------
C     rangement dans TYP-COUCHE de NP pour les couches concernees 
C -----------------------------------------------------------------------
            NP=NP+1
C 
CD           CALL IMPEP('NP=',NP)
C 
          CALL MESSAO ('Numero des couches concernees ( -1 => TOUTE )')
C 
          CALL LECLEN (M(APM2), NBCOU, NCOULU)
C 
          IF (M(APM2) .EQ. -1) THEN
            DO I = 1, NBCOU
                M(R1+I)=NP
            ENDDO
            NP1=NP1+NBCOU
          ELSE 
            DO I=1,NCOULU
              M(R1+M(APR2+I))=NP
            ENDDO
            NP1=NP1+NCOULU
          END IF
C 
          IF (NP1.GT.NBCOU) THEN
            CALL MESSAO ('trop grand nombre de couche
     $                   \recommencez la rentree de donnee ')
            GOTO 4
          ENDIF
C  
      ENDDO
C 
C -----------------------------------------------------------------------
C     creation d'un tableau pour les constantes d'elasticite
C -----------------------------------------------------------------------
C 
      NKCOU = NP  
      CALL IMPET ('NOMBRE DE TYPES DE COMPORTEMENT DE COUCHE ', NKCOU)
      NARG1 = 17*NKCOU
      NARG2 = 9*NKCOU
C 
      CALL GESTDP ('SOUP-ORTHO', NARG1, SOUPLO)
      CALL COPITD (NARG1, DM(ADSOUP), DM(SOUPLO))
      ADSOUP = SOUPLO
C  
      CALL GESTDP ('ELAS-ORTHO', NARG2, ELASLO)
      CALL COPITD (NARG2, DM(ADELAS), DM(ELASLO))
      ADELAS = ELASLO
C 
      NTCEND  = NKCOU*NBCNLC
      CALL GESTDP ('CARAC-NONL', NTCEND, CAENLO)
      CALL COPITD (NTCEND, DM(ADCAEN), DM(CAENLO))
      ADCAEN = CAENLO
C -----------------------------------------------------------------------
C     creation d'un tableau pour les constantes d'elasticite des interfaces
C     NTTEI est le nombre maxi de termes d'elasticite d'interface
C -----------------------------------------------------------------------
C 
      NTTEI=3*NBINT
C 
      CALL GSPOUD (NTTEI, ADSOUI)     
C 
C -----------------------------------------------------------------------
C     DR Ajout le 22/01/96 : Parametre A en plus en non lineaire
C     DL Ajout le 25/07/96 : Parametre y0, alp, m en plus en non lineaire
C     DL Ajout le 03/09/96 : Parametre yr (seuil fragile en arrachement)
C     SG Ajout le 23/03/98 : Parametres CRITERE ELASTIQUES
C -----------------------------------------------------------------------
      NBCNLI = 17
C 
      KISYM = 0        
C 
      CALL GSPOUD (NBCNLI*NBINT, ADIAEN)
C 
C -----------------------------------------------------------------------
C     NI est le nombre de comportement d'interface deja rentre
C -----------------------------------------------------------------------
      NI=0            
      RLOC=R1+NBCOU
C 
C -----------------------------------------------------------------------
C     tant que toute les interfaces n'ont pas de comportement
C     NP2 est le nombre d'interface ayant deja un comportement
C -----------------------------------------------------------------------
5     NP2=0
C 
      DO WHILE(NP2.NE.NBINT)
C 
CD          CALL IMPMP('AVANT APPEL A LECSEQ A RIGINT')
C 
	    CALL MESSAO(
     $     '\  donnee de Ei tels que : 
     $      \   si3 =  Eix[|U3|]   ')
C  
            CALL LECSEQ('RIGINT','valeur des termes de la matrice')
C 
C -----------------------------------------------------------------------
C     NT est la 1ere adresse libre dans SOIN-ORTHO
C -----------------------------------------------------------------------
            NT=3*NI+ADSOUI
C 
C -----------------------------------------------------------------------
C     affectation des valeurs lues 
C -----------------------------------------------------------------------
            E1 = LECDBL('E1')
            E2 = LECDBL('E2')
            E3 = LECDBL('E3')
C 
            S1 = 1.D0/E1
            S2 = 1.D0/E2
            S3 = 1.D0/E3
C 
C -----------------------------------------------------------------------
C     REMPLISSAGE DE SOIN-ORTHO
C -----------------------------------------------------------------------
            NSOIOR       = ADSOUI+3*NI
            DM(NSOIOR)   = S1
            DM(NSOIOR+1) = S2
            DM(NSOIOR+2) = S3
C 
C -----------------------------------------------------------------------
C     NNONL est la 1ere adresse libre dans CARAI-NONL
C -----------------------------------------------------------------------
            NNONL=NBCNLI*NI+ADIAEN
C 
C -----------------------------------------------------------------------
C     Rentree des caracteristiques de plasticite      
C -----------------------------------------------------------------------
            CALL MESSAO (
     $     '\   Soient :
     $      \   R( p ) = Ro + BETA*( p ** ALPHA )  
     $      \   ou le seuil est defini par :
     $      \   f = seq - s(3,3) -R < ou = 0
     $      \                       2      2  1/2
     $      \   avec : seq =( A1.s13 +A2s23 )      '  )
C  
            CALL LECSEQ ('IPLAST',
     $      ' donnee de RO , BETA  , ALPHA , A1 , A2 ')
C 
            R0       =  LECDBL('R0')
            ALPHA    =  LECDBL('ALPHA')
            BETA     =  LECDBL('BETA')
            A1       =  LECDBL('A1')
            A2       =  LECDBL('A2')
C 
            DM( NNONL)       =  R0
            DM( NNONL+1)     =  ALPHA
            DM( NNONL+2)     =  BETA
            DM( NNONL+3)     =  A1
            DM( NNONL+4)     =  A2
C  
            CALL MESSAO(
     $      '  Soient :
     $      \           2                    2
     $      \  Yd1 = s13 / 2.k10 , Yd2  = s23 / 2.k20
     $      \                 2  
     $      \  Yd = sup(0,s33) / 2.k0
     $      \                                               
     $      \                 alpha             alpha    alpha
     $      \  Y  = (GAM1.Yd1)      + (GAM2.Yd2)     + Yd     
     $      \ 
     $      \       (1/alpha)
     $      \  Y = Y
     $      \                        2
     $      \  Y  =  Y - a.inf(s33,0) / 2.k0
     $      \                
     $      \  Y  = sup(0,Y)    
     $      \  .                                 m       n
     $      \  d  = k.(  ((m/m+1)*(Y-Y0)/(Yc-Y0))  -  d )
     $      \                        
     $      \  d1 = d2 = d
     $      \                        
     $      \  d = 1  si  Yd > Yr         ' )
C  
            CALL LECSEQ ('IENDOM',
     $  'Donnee de GAM1, GAM2, a, k, Y0, Yc, Yr, ALP, M, n, sr13, sr23')
C 
            K     =  LECDBL('K')
            N     =  LECDBL('N')
            Y0    =  LECDBL('Y0')
            YC    =  LECDBL('YC')
            YR    =  LECDBL('YR')
            ALP   =  LECDBL('ALP')
            MINT  =  LECDBL('M')
            G1    =  LECDBL('GAM1')
            G2    =  LECDBL('GAM2')
            A     =  LECDBL('A')
            SRI13 =  LECDBL('SRI13')
            SRI23 =  LECDBL('SRI23')
	    
C 
            DM( NNONL+5)  =  K
            DM( NNONL+6)  =  N
            DM( NNONL+7)  =  Y0
            DM( NNONL+8)  =  YC
            DM( NNONL+9)  =  YR
            DM( NNONL+10) =  ALP
            DM( NNONL+11) =  MINT
            DM( NNONL+12) =  G1
            DM( NNONL+13) =  G2
            DM( NNONL+14) =  A
            DM( NNONL+15) =  SRI13
            DM( NNONL+16) =  SRI23           	    
C 
C -----------------------------------------------------------------------
C     rangement dans TYP-COUCHE de NI pour les couches concernees
C -----------------------------------------------------------------------
            NI = NI+1
C 
            CALL MESSAD(
     &      'Numero des interfaces concernees (-1 => toute)')
C 
            CALL LECLEN(M(APM2),NBINT,NCOULU)
C                                      
            IF (M(APM2) .EQ. -1)  THEN
C 
              NP2=NP2+NBINT
              INTSYM   = .FALSE. 
C 
              IF(SYMPAR)THEN 
                INTSYM = .TRUE.
                IF ( NBINT .GT.1) THEN 
C 
C -----------------------------------------------------------------------
C     on attribue au comportement de l'interface centrale 
C     le numero 2
C -----------------------------------------------------------------------
                  M(RLOC+1)= NI+1
                  ISYM     = 1
                  KISYM    = NI+1
                  DO I=2,NBINT
                    M(RLOC+I)=NI
                  ENDDO
                ELSE
C 
C -----------------------------------------------------------------------
C     Il y a une seule interface
C -----------------------------------------------------------------------
                  KISYM     = 1
                  NI        = NI-1  
                  M(RLOC+1) =NI+1               
                END IF
              ELSE 
                DO I = 1, NBINT
                  M(RLOC+I)=NI
                END DO 
              END IF
C 
            ELSE 
C                          
              INTSYM = .FALSE.
C 
              DO I=1,NCOULU
C                                                 
                IF( SYMPAR .AND. M(APM2+I-1) .EQ.1 ) THEN
                  INTSYM = .TRUE.
                  IF ( NCOULU .GT. 1) THEN 
                    ISYM   = I            
                    M(RLOC+M(APM2+I-1))=NI+1
                    KISYM              =NI+1 
                  ELSE 
                    KISYM              =NI              
                    M(RLOC+M(APM2+I-1))=NI
                    NI                 = NI-1   
                  END IF
                ELSE 
                  M(RLOC+M(APM2+I-1))=NI
                END IF 
C 
              END DO
C 
              NP2=NP2+NCOULU
C 
            ENDIF
C 
            IF (INTSYM ) THEN 
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau de comportement particulier
C     pour l'interface centrale
C -----------------------------------------------------------------------
              NSOIOR  = ADSOUI+3*NI
              NNONL   = NBCNLI*NI+ADIAEN
C 
C -----------------------------------------------------------------------
C     Comme on utilise DM(NSOIOR) uniquemeNt pour des produits scalaire
C     on impose egalement S1=S2=0 => le saut plan est nul
C     Le saut de w etant le double du deplacememnt normal de 
C     l'interface on impose egalement S3 = .5D0*S3
C -----------------------------------------------------------------------
              DM(NSOIOR)   = 0.D0
              DM(NSOIOR+1) = 0.D0
              DM(NSOIOR+2) = 0.5D0*S3
C 
C -----------------------------------------------------------------------
C     NNONL est la 1ere adresse libre dans CARAI-NONL
C -----------------------------------------------------------------------
              DM( NNONL)      =  R0
              DM( NNONL+1)    =  ALPHA
              DM( NNONL+2)    =  BETA
              DM( NNONL+3)    =  A1
              DM( NNONL+4)    =  A2
C  
              DM( NNONL+5)  =  K
              DM( NNONL+6)  =  N
              DM( NNONL+7)  =  Y0
              DM( NNONL+8)  =  YC
              DM( NNONL+9)  =  YR
              DM( NNONL+10) =  ALP
              DM( NNONL+11) =  MINT
              DM( NNONL+12) =  G1
              DM( NNONL+13) =  G2
              DM( NNONL+14) =  A
              DM( NNONL+15) =  SRI13
              DM( NNONL+16) =  SRI23           
C 
              NI = NI+1  
C 
          END IF
C 
          IF(NP2.GT.NBINT)THEN
              CALL MESSAO('trop grand nombre d''interfaces
     $                    \recommencez la rentree de donnee ')
              GOTO 5
          ENDIF
C 
      ENDDO
C 
      NKINT  = NI
C 
      CALL IMPET ( 'nb de type de comportement d''interface', NKINT )
C 
C -----------------------------------------------------------------------            
C     NTTEI est le nombre total de termes d'elasticite d'interface
C -----------------------------------------------------------------------
      NTTEI=3*NKINT
C 
      CALL GESTDP ('SOIN-ORTHO',NTTEI,SOUPLI)     
      CALL COPITD (NTTEI, DM(ADSOUI), DM(SOUPLI))
      ADSOUI =SOUPLI
C 
      CALL GESTDP ('CARAI-NONL', NBCNLI*NKINT, CAENLI)
      CALL COPITD (NBCNLI*NKINT, DM(ADIAEN), DM(CAENLI))
      ADIAEN = CAENLI
C 
C -----------------------------------------------------------------------
C     SEQUENCE D'IMPRESSION DES RENTREES POUR LES COUCHES
C -----------------------------------------------------------------------
      DO  I = 1 , NBCOU
C 
        WRITE(NUMERO(1:4),'(I4)')I
C 
        CALL MESSAO ('Pour la couche numero : '//numero )
        CALL IMPDT ('L''angle est en degres ' , DM(AC3+I)*180./PI )
C 
        P=M(M1+I-1)
C 
        CALL IMPET  ('Le numero de comportement est ',P)
        CALL IMPTDT ('La souplesse elastique est caracterisee par ',
     &                DM(SOUPLO+17*(P-1)), 1, 17)
        CALL IMPTDT ('les coefficients elastiques sont ',
     &                DM(ELASLO+9*(P-1)), 1, 9)
        CALL IMPTDT ('Les caracteristiques non lineaires sont ',
     &                DM(CAENLO+NBCNLC*(P-1)), 1, NBCNLC)
      END DO
C 
C -----------------------------------------------------------------------
C     POUR LES INTERFACES
C -----------------------------------------------------------------------                              
      DO  I = 1 , NBINT
C 
        WRITE(NUMERO(1:4),'(I4)')I
C 
        CALL MESSAO ('Pour l''interface : '//numero )
        CALL IMPDT  ('L''angle est en degres ' , 
     &                DM(AC3+NBCOU+I)*180./PI )
C 
        P=M(M1+NBCOU+I-1)
C 
        CALL IMPET ( ' Le numero de comportement est' ,P )
        CALL IMPTDT ( ' La souplesse elastique est caracterisee par ',
     &    DM(SOUPLI+3*(P-1)), 1, 3 )
        CALL IMPTDT ( ' Les caracteristiques non lineaires sont ',
     &    DM(CAENLI+NBCNLI*(P-1)), 1, NBCNLI )

      END DO
C 
C -----------------------------------------------------------------------
C     TRAITEMENT DES INFOS PRECEDENTES
C -----------------------------------------------------------------------
      CALL TRAINL(KISYM)
C 
      CALL SOPOUB(AM2LC,ADM2LC)
      CALL RETOUD(IDPROG)
C 
      RETURN
      END
