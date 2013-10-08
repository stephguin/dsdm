      SUBROUTINE LFIDON
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'cominc_visu.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  U, V, I, NBORLU, ABORD, NUMDIR, IUNIT, J
      INTEGER  DEVREE, DEVIMA, DEVDRE, DEVDIM, K
      INTEGER  AURREE, AVRREE, ADTETA, ADHAU , ADCOTE
      INTEGER  DEVFU, DEVFV, DEVFW, DEVFRR, DEVFR0
      INTEGER  W, RR, R0
      INTEGER  LGCARG, LNOMFI
C 
      DOUBLE PRECISION DHAUT
      DOUBLE PRECISION UX, UY, RX, RY
      DOUBLE PRECISION C, S
C 
      CHARACTER*9   NOMDIR
      CHARACTER*6   IDPROG
      CHARACTER*20  NOMFIC
      CHARACTER*120 NOM
C 
      INTEGER  AM2LC, ADM2LC
C 
      PARAMETER (IDPROG='LFIDON')
C 
C     include 'identi.h'
C 
C -----------------------------------------------------------------------
      CHARACTER*6      CPLAST(6) , IPLAST(2), CRITER(2)
      CHARACTER*3      CENDOM(3) , ENDOMI(3)
      CHARACTER*10     ENDDIF(3)
      CHARACTER*9      BORD(4)
      CHARACTER*3      CONTRA(6), DEFORM(6)
      CHARACTER*4      DCONTR(6)
      CHARACTER*4      DDEFOR(6)
      CHARACTER*1      DEPLA(3)
C 
      BORD(1)    = 'inferieur'
      BORD(2)    = 'interieur'
      BORD(3)    = 'superieur'
      BORD(4)    = 'exterieur'
C 
      DEPLA(1)   = 'u'
      DEPLA(2)   = 'v'
      DEPLA(3)   = 'w'
C 
C  SGU 05/09/13, VISORT pas clair
C
      VISORT=.TRUE.
      VISUXY=.FALSE.
      IF (VISORT) THEN
C 
         CONTRA(1)  = 'C11'
         CONTRA(2)  = 'C22'
         CONTRA(3)  = 'C12'
         CONTRA(4)  = 'C23'
         CONTRA(5)  = 'C13'
         CONTRA(6)  = 'C33'
C 
         DEFORM(1)  = 'E11'
         DEFORM(2)  = 'E22'
         DEFORM(3)  = 'E12'
         DEFORM(4)  = 'E23'
         DEFORM(5)  = 'E13'
         DEFORM(6)  = 'E33'
C 
         DCONTR(1)  = 'DC11'
         DCONTR(2)  = 'DC22'
         DCONTR(3)  = 'DC12'
         DCONTR(4)  = 'DC23'
         DCONTR(5)  = 'DC13'
         DCONTR(6)  = 'DC33'
C 
         DDEFOR(1)  = 'DE11'
         DDEFOR(2)  = 'DE22'
         DDEFOR(3)  = 'DE12'
         DDEFOR(4)  = 'DE23'
         DDEFOR(5)  = 'DE13'
         DDEFOR(6)  = 'DE33'
C 
      ELSE IF (VISUXY) THEN
C 
         CONTRA(1)  = 'Cxx'
         CONTRA(2)  = 'Cyy'
         CONTRA(3)  = 'Cxy'
         CONTRA(4)  = 'Cyz'
         CONTRA(5)  = 'Cxz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Exx'
         DEFORM(2)  = 'Eyy'
         DEFORM(3)  = 'Exy'
         DEFORM(4)  = 'Eyz'
         DEFORM(5)  = 'Exz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCxx'
         DCONTR(2)  = 'DCyy'
         DCONTR(3)  = 'DCxy'
         DCONTR(4)  = 'DCyz'
         DCONTR(5)  = 'DCxz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DExx'
         DDEFOR(2)  = 'DEyy'
         DDEFOR(3)  = 'DExy'
         DDEFOR(4)  = 'DEyz'
         DDEFOR(5)  = 'DExz'
         DDEFOR(6)  = 'DEzz'
C 
      ELSE
C 
         CONTRA(1)  = 'Crr'
         CONTRA(2)  = 'C00'
         CONTRA(3)  = 'Cr0'
         CONTRA(4)  = 'C0z'
         CONTRA(5)  = 'Crz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Err'
         DEFORM(2)  = 'E00'
         DEFORM(3)  = 'Er0'
         DEFORM(4)  = 'E0z'
         DEFORM(5)  = 'Erz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCrr'
         DCONTR(2)  = 'DC00'
         DCONTR(3)  = 'DCr0'
         DCONTR(4)  = 'DC0z'
         DCONTR(5)  = 'DCrz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DErr'
         DDEFOR(2)  = 'DE00'
         DDEFOR(3)  = 'DEr0'
         DDEFOR(4)  = 'DE0z'
         DDEFOR(5)  = 'DErz'
         DDEFOR(6)  = 'DEzz'
C 
      ENDIF
C 
      CPLAST(1)  = 'epsp11' 
      CPLAST(2)  = 'epsp22' 
      CPLAST(3)  = 'epsp12' 
      CPLAST(4)  = 'epsp23' 
      CPLAST(5)  = 'epsp13' 
      CPLAST(6)  = 'epsp33' 
C 
      CENDOM(1)   = 'dfi'
      CENDOM(2)   = 'dps' 
      CENDOM(3)   = 'dpt'
C 
      IPLAST(1)  = 'SAUTP1' 
      IPLAST(2)  = 'SAUTP2' 
C 
      ENDOMI(1)   = 'di1' 
      ENDOMI(2)   = 'di2'
      ENDOMI(3)   = 'di3'
C 
      ENDDIF(1)   = 'di1-di1ini'
      ENDDIF(1)   = 'di2-di2ini'
      ENDDIF(1)   = 'di3-di3ini'
c -
      CRITER(1)  = 'CRIT-S'
      CRITER(2)  = 'CRIT-N'
C -----------------------------------------------------------------------
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------      
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
C 
C     Le tableau HAU-CE-TOT contient les distances au centre
C     de tous les plans en partant du plan inferieur du stratifie
C 
      CALL ADTBDM ('HAU-CE-TOT', ADHAU)
C 
      CALL GESTDP ('HAUT-NOEUD', 3*NBCOU, ADCOTE)
C 
C     On remplit le tableau des hauteurs des noeuds
C 
      K = 0
      DO I = 1, NBCOU
        DHAUT = DABS((DM(ADHAU+I)-DM(ADHAU+I-1)) /2.D0)
        DM(ADCOTE+K) = DM(ADHAU+I-1)
        K = K+1
        DO J = 1, 2
          DM(ADCOTE+K) = DHAUT + DM(ADCOTE+K-1)
          K = K+1
        END DO
      END DO
C 
CD    CALL IMPDT (' demi-epaisseur totale ', THICK )
CD    CALL IMPTDT ('HAUT-NOEUD DANS '//IDPROG, DM(ADCOTE), 3*NBCOU, 1)
C 
C -----------------------------------------------------------------------
C     Entree dans la sequence de lecture des fichiers de donnees
C -----------------------------------------------------------------------
C 
      IF (DONFIC .AND. TRACT) THEN
C 
         CALL MESSAO ('FICHIER DE DONNEES DE TYPE TENSION')
         CALL MESSAO ('LECTURE DES FICHIERS DE DONNEES DANS '//IDPROG)
C 
         CALL ADTBM ('TYP-FICHIE',  ABORD )
         CALL LONGEN('TYP-FICHIE', NBORLU)
C 
         NUMDIR = 12
         NOMDIR = 'donneefic'
C 
         DO I = 1 , NBORLU
           NOMFIC = CHAFIC( 12 , I )
           LNOMFI = LGCARG ( CHAFIC( 12 , I ) ,20)
           CALL MESSAO (
     &          'POUR LE BORD'//bord(M(ABORD +I-1))//
     &          '\LE FICHIER DE DONNEES DANS LE REPERTOIRE'//
     &          '\DONNEEFIC A POUR NOM : '//NOMFIC)
         END DO
C 
C -----------------------------------------------------------------------
C     Tableau des valeurs reelles de Ur et U0 en fonction des angles
C -----------------------------------------------------------------------
         CALL GESTDP('U-R-IMP-FI', NTETA , U  )
         CALL GESTDP('V-R-IMP-FI', NTETA , V  )
C 
	 CALL OUVFCD( NUMDIR , NOMFIC  ,LNOMFI , 'F' , NTETA ,IUNIT )
         CALL IMPET('IUNIT', IUNIT)
         INQUIRE( UNIT= IUNIT , NAME=NOM)
CD       CALL IMPCT('POUR PATHNAME DANS '//IDPROG,   NOM )
C 
CD       NOM = tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LNOMFI)
CD       CALL IMPCT('POUR PATHNAME DANS '//IDPROG,   NOM )
C 
         DO I = 0 , NTETA-1
	   READ(IUNIT ,*) DM(U+I), DM(V+I)
         END DO
C 
         CALL IMPTDT ('VALEURS DE U LUES '//IDPROG , DM(U) , 1 , NTETA )
         CALL IMPTDT ('VALEURS DE V LUES '//IDPROG , DM(V) , 1 , NTETA )
C 
         IF ( .NOT. POLAR ) THEN
           CALL MESSAO ('TRANSFORMATION EN COORDONNEES POLAIRES')
           DO I = 0 , NTETA-1
             C       = DCOS(DM(ADTETA+I))
             S       = DSIN(DM(ADTETA+I))
             UX      = DM(U+I)
             UY      = DM(V+I)
             DM(U+I) = C*UX+S*UY
             DM(V+I) = C*UY-S*UX
           END DO
         END IF
C 
C -----------------------------------------------------------------------
C     Traitement des donnees :
C 
C     Il faut :  a ) calculer la valeur des deplacements sur tout le bord
C                b ) developper les deplacements en series de Fourier
C                c ) modifier les seconds membres en consequence
C                d ) modifier les diagonales des matrices K0n en consequence
C 
C     Realisation du a : calcul de la valeur des deplacements sur tout le bord
C 
C     Pour le moment en tension, principe :
C 
C         u et v sont identiques sur tout le bord, w depend de la couche
C         de facon a assurer que la contrainte normale est nulle
C         ( remarque grad w est nul )
C                        m
C          comme pour le moment on n'a pas la valeur de la deformation de tension
C          on ne peut pas faire mieux qu'imposer w nul sur la ligne moyenne
C 
C          Donc on developpe u et v en series de Fourier (w etant nul)
C 
C -----------------------------------------------------------------------
C     Adresses des debuts des tableaux provisoires :
C -----------------------------------------------------------------------
C 
C     Tableau des valeurs developpees en series de Fourier
C     de Ur et U0 en fonction des angles
C 
      CALL GESTDP('U-T-IMP-FI', NBMAT , DEVREE )
      CALL GESTDP('V-T-IMP-FI', NBMAT , DEVIMA )
C 
C     Tableau des valeurs developpees en series de Fourier
C     des derivees par rapport a teta de Ur et U0 en fonction
C     des angles, de facon a calculer eventuellement les gradients
C 
      CALL GESTDP('U0-T-IM-FI', NBMAT , DEVDRE )
      CALL GESTDP('V0-T-IM-FI', NBMAT , DEVDIM )
C 
CD    CALL MESSAO('CALCUL DU DEVELOPPEMENT DE U')
      CALL DEVSFO( NTETA , DM(U) , 0 , DM(DEVREE) )
C 
C ----------------------------------------------------------------------- 
C     Calcul du developpement en series de Fourier des derivees
C     en teta de u puis de leurs valeurs reelles
C -----------------------------------------------------------------------
      CALL POUSMD( 2*NTETA ,  AURREE )
      AVRREE = AURREE+NTETA
C 
CD    CALL MESSAO('CALCUL DES VALEURS REELLES DES DERIVEES DE U')
      CALL DERTET( 0, DM(DEVREE), DM(DEVDRE) )
C 
      CALL VAREFO( NTETA, DM(DEVDRE), 0, DM(AURREE) )
C 
CD    CALL MESSAO('CALCUL DU DEVELOPPEMENT DE V')
      CALL DEVSFO( NTETA, DM(V), 1, DM(DEVIMA) )
C 
C ----------------------------------------------------------------------- 
C     Calcul du developpement en serieS de Fourier des derivees
C     en teta de v
C -----------------------------------------------------------------------
CD    CALL MESSAO ('CALCUL DES VALEURS REELLES DES DERIVEES DE V')
      CALL DERTET ( 1, DM(DEVIMA), DM(DEVDIM) )
      CALL VAREFO ( NTETA, DM(DEVDIM), 1 ,DM(AVRREE) )
C 
      CALL TDEFIT
C 
      ELSE IF ( DONFIC .AND. (.NOT. TRACT) ) THEN
C 
         CALL MESSAO ('FICHIER DE DONNEES DE TYPE FLEXION')
         CALL MESSAO ('LECTURE DES FICHIERS DE DONNEES DANS '
     &                //IDPROG )
C 
         CALL ADTBM ('TYP-FICHIE', ABORD)
         CALL LONGEN ('TYP-FICHIE', NBORLU)
C 
         NUMDIR = 12
         NOMDIR = 'donneefic'
C 
         DO I = 1 , NBORLU
           NOMFIC = CHAFIC(12, I)
           LNOMFI = LGCARG (CHAFIC(12, I), 20)
           CALL MESSAO (
     &          'POUR LE BORD '// bord(M(ABORD +I-1))//
     &          ' LE FICHIER DE DONNEES DANS LE REPERTOIRE '//
     &          '\donneefic A POUR NOM : '// NOMFIC )
         END DO
C 
C        Tableau des valeurs reelles de Ur, U0 et W, Rer, Re0
C        en fonction des angles
C 
        CALL GESTDP ('U-R-IMP-FI', NTETA , U)
        CALL GESTDP ('V-R-IMP-FI', NTETA , V)
        CALL GESTDP ('W-R-IMP-FI', NTETA , W)
        CALL GESTDP ('RRR-IMP-FI', NTETA , RR)
        CALL GESTDP ('R0R-IMP-FI', NTETA , R0)
C 
        CALL OUVFCD (NUMDIR, NOMFIC, LNOMFI, 'F', NTETA, IUNIT)
C 
        DO I = 0, NTETA-1
          READ (IUNIT,*) DM(U+I), DM(V+I), DM(W+I), DM(RR+I), DM(R0+I)
        END DO
C 
        CALL IMPTDT('VALEUR DE U LUES  '//IDPROG, DM(U) , 1 , NTETA )
        CALL IMPTDT('VALEUR DE V LUES  '//IDPROG, DM(V) , 1 , NTETA )
        CALL IMPTDT('VALEUR DE W LUES  '//IDPROG, DM(W) , 1 , NTETA )
        CALL IMPTDT('VALEUR DE RR LUES '//IDPROG, DM(RR), 1 , NTETA )
        CALL IMPTDT('VALEUR DE R0 LUES '//IDPROG, DM(R0), 1 , NTETA )
C 
        IF (.NOT. POLAR) THEN
C 
          CALL MESSAO ('TRANSFORMATION EN COORDONNEES POLAIRES')
C 	  
          DO I = 0, NTETA-1
            C        = DCOS(DM(ADTETA+I))
            S        = DSIN(DM(ADTETA+I))
            UX       = DM(U+I)
            UY       = DM(V+I)
            RX       = DM(RR+I)
            RY       = DM(R0+I)
            DM(U+I)  = C*UX+S*UY
            DM(V+I)  = C*UY-S*UX
            DM(RR+I) = C*RX+S*RY
            DM(R0+I) = C*RY-S*RX
          END DO
C 
	  CALL IMPTDT ('VALEUR DE U  '//IDPROG, DM(U) , 1, NTETA)
          CALL IMPTDT ('VALEUR DE V  '//IDPROG, DM(V) , 1, NTETA)
          CALL IMPTDT ('VALEUR DE W  '//IDPROG, DM(W) , 1, NTETA)
          CALL IMPTDT ('VALEUR DE Rx '//IDPROG, DM(RR), 1, NTETA)
          CALL IMPTDT ('VALEUR DE Ry '//IDPROG, DM(R0), 1, NTETA)
C 
        END IF
C 
C       Traitement des donnees :
C 
C       Il faut :  a ) calculer la valeur des deplacements sur tout le bord
C                  b ) developper les deplacements en series de Fourier
C                  c ) modifier les seconds membres en consequence
C                  d ) modifier les diagonales des matrices K0n en consequence
C 
C       Realisation du a : calcul de la valeur des deplacements sur tout le bord
C 
C       En flexion principe :
C 
C       u et v, w rer et re0 sont calcules sur la ligne moyenne
C       Donc on developpe u, v, w, rer, re0 en series de Fourier
C       Le calcul des deplacements dans l'epaisseur sera effectue dans
C 
C       Tableau des valeurs developpees en serie de Fourier
C       de Ur et U0 W RR R0 en fonction des angles
C 
        CALL GESTDP ('U-T-IMP-FI', NBMAT , DEVFU )
        CALL GESTDP ('V-T-IMP-FI', NBMAT , DEVFV )
        CALL GESTDP ('W-T-IMP-FI', NBMAT , DEVFW )
        CALL GESTDP ('RRT-IMP-FI', NBMAT , DEVFRR )
        CALL GESTDP ('R0T-IMP-FI', NBMAT , DEVFR0 )
C 
CD      CALL MESSAO ('Calcul du developpement de U')
C 
        CALL DEVSFO ( NTETA , DM(U) , 0 , DM(DEVFU) )
C 
CD      CALL MESSAO('Calcul du developpement de V')
        CALL DEVSFO( NTETA , DM(V) , 1 , DM(DEVFV) )
C 
CD      CALL MESSAO('Calcul du developpement de W')
        CALL DEVSFO( NTETA , DM(W) , 0 , DM(DEVFW) )
C 
CD      CALL MESSAO('Calcul du developpement de RR')
        CALL DEVSFO( NTETA , DM(RR) , 1 , DM(DEVFRR) )
C 
CD      CALL MESSAO('Calcul du developpement de R0')
        CALL DEVSFO( NTETA , DM(R0) , 0 , DM(DEVFR0) )
C 
        CALL MESSAO ('DEVELOPPEMENT EN SERIES DE FOURIER')
C
  	  CALL IMPTDT('DEV DE U   '//IDPROG,
     &                    DM(DEVFU) , 1 , NBMAT)
          CALL IMPTDT('DEV DE V   '//IDPROG,
     &                    DM(DEVFV) , 1 , NBMAT)
          CALL IMPTDT('DEV DE W   '//IDPROG,
     &                    DM(DEVFW) , 1 , NBMAT)
          CALL IMPTDT('DEV DE RR  '//IDPROG,
     &                    DM(DEVFRR), 1 , NBMAT)
          CALL IMPTDT('DEV DE R0  '//IDPROG,
     &                    DM(DEVF00), 1 , NBMAT)
C
	CALL TDEFIF
C 
      END IF
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
C     Traitement des deplacements donnes par fichier en traction
C     Raccord a une solution plaque en deplacement en traction

      SUBROUTINE TDEFIT
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'cominc_visu.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  ADLOC1, ADDGDP
      INTEGER  NUCO, TYPDEP, NBDDL, NBTDDL
      INTEGER  NUBORD, TEST1
      INTEGER  AM2LC, ADM2LC
C 
C     Pour les deplacements imposes pour U et V
C 
      INTEGER NBTDDU(2), DDLU(2), ADDEP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TDEFIT')
C 
C     include 'identi.h'
C 
C -----------------------------------------------------------------------
      CHARACTER*6      CPLAST(6) , IPLAST(2), CRITER(2)
      CHARACTER*3      CENDOM(3) , ENDOMI(3)
      CHARACTER*10     ENDDIF(3)
      CHARACTER*9      BORD(4)
      CHARACTER*3      CONTRA(6), DEFORM(6)
      CHARACTER*4      DCONTR(6)
      CHARACTER*4      DDEFOR(6)
      CHARACTER*1      DEPLA(3)
C 
      BORD(1)    = 'inferieur'
      BORD(2)    = 'interieur'
      BORD(3)    = 'superieur'
      BORD(4)    = 'exterieur'
C 
      DEPLA(1)   = 'u'
      DEPLA(2)   = 'v'
      DEPLA(3)   = 'w'
C 
C  SGU 05/09/13, VISORT pas clair
C
      VISORT=.TRUE.
      VISUXY=.FALSE.
      IF (VISORT) THEN
C 
         CONTRA(1)  = 'C11'
         CONTRA(2)  = 'C22'
         CONTRA(3)  = 'C12'
         CONTRA(4)  = 'C23'
         CONTRA(5)  = 'C13'
         CONTRA(6)  = 'C33'
C 
         DEFORM(1)  = 'E11'
         DEFORM(2)  = 'E22'
         DEFORM(3)  = 'E12'
         DEFORM(4)  = 'E23'
         DEFORM(5)  = 'E13'
         DEFORM(6)  = 'E33'
C 
         DCONTR(1)  = 'DC11'
         DCONTR(2)  = 'DC22'
         DCONTR(3)  = 'DC12'
         DCONTR(4)  = 'DC23'
         DCONTR(5)  = 'DC13'
         DCONTR(6)  = 'DC33'
C 
         DDEFOR(1)  = 'DE11'
         DDEFOR(2)  = 'DE22'
         DDEFOR(3)  = 'DE12'
         DDEFOR(4)  = 'DE23'
         DDEFOR(5)  = 'DE13'
         DDEFOR(6)  = 'DE33'
C 
      ELSE IF (VISUXY) THEN
C 
         CONTRA(1)  = 'Cxx'
         CONTRA(2)  = 'Cyy'
         CONTRA(3)  = 'Cxy'
         CONTRA(4)  = 'Cyz'
         CONTRA(5)  = 'Cxz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Exx'
         DEFORM(2)  = 'Eyy'
         DEFORM(3)  = 'Exy'
         DEFORM(4)  = 'Eyz'
         DEFORM(5)  = 'Exz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCxx'
         DCONTR(2)  = 'DCyy'
         DCONTR(3)  = 'DCxy'
         DCONTR(4)  = 'DCyz'
         DCONTR(5)  = 'DCxz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DExx'
         DDEFOR(2)  = 'DEyy'
         DDEFOR(3)  = 'DExy'
         DDEFOR(4)  = 'DEyz'
         DDEFOR(5)  = 'DExz'
         DDEFOR(6)  = 'DEzz'
C 
      ELSE
C 
         CONTRA(1)  = 'Crr'
         CONTRA(2)  = 'C00'
         CONTRA(3)  = 'Cr0'
         CONTRA(4)  = 'C0z'
         CONTRA(5)  = 'Crz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Err'
         DEFORM(2)  = 'E00'
         DEFORM(3)  = 'Er0'
         DEFORM(4)  = 'E0z'
         DEFORM(5)  = 'Erz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCrr'
         DCONTR(2)  = 'DC00'
         DCONTR(3)  = 'DCr0'
         DCONTR(4)  = 'DC0z'
         DCONTR(5)  = 'DCrz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DErr'
         DDEFOR(2)  = 'DE00'
         DDEFOR(3)  = 'DEr0'
         DDEFOR(4)  = 'DE0z'
         DDEFOR(5)  = 'DErz'
         DDEFOR(6)  = 'DEzz'
C 
      ENDIF
C 
      CPLAST(1)  = 'epsp11' 
      CPLAST(2)  = 'epsp22' 
      CPLAST(3)  = 'epsp12' 
      CPLAST(4)  = 'epsp23' 
      CPLAST(5)  = 'epsp13' 
      CPLAST(6)  = 'epsp33' 
C 
      CENDOM(1)   = 'dfi'
      CENDOM(2)   = 'dps' 
      CENDOM(3)   = 'dpt'
C 
      IPLAST(1)  = 'SAUTP1' 
      IPLAST(2)  = 'SAUTP2' 
C 
      ENDOMI(1)   = 'di1' 
      ENDOMI(2)   = 'di2'
      ENDOMI(3)   = 'di3'
C 
      ENDDIF(1)   = 'di1-di1ini'
      ENDDIF(1)   = 'di2-di2ini'
      ENDDIF(1)   = 'di3-di3ini'
c -
      CRITER(1)  = 'CRIT-S'
      CRITER(2)  = 'CRIT-N'
C -----------------------------------------------------------------------
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     adresse du tableau de localisation des numeros
C     de noeuds ranges par ordre croissant par elements
C     ranges couche par couche
C -----------------------------------------------------------------------
      CALL ADTBM ('TLOCN1    ', ADLOC1)
C 
C     POUR REMPLIR LE TABLEAU DDL-DPI-FI
C 
      CALL DEBUEN( ADDGDP )
C 
C     POUR REMPLIR LES TABLEAUX 'DDL-UI-FI4' , 'DDL-VI-FI4'
C 
      CALL POUSME( 2*NDDL ,  DDLU(1) )
C 
      DDLU(2) = DDLU(1)+NDDL
C 
C     Interet : on impose la forme du deplacement sur tout un bord
C     exemple : raccord avec une solution plaque
C 
      NBTDDL    = 0
      NBTDDU(1) = 0
      NBTDDU(2) = 0
C 
C     Les deplacements sont imposes sur le bord exterieur
C 
      NUBORD = 4
C 
      DO TYPDEP = 1 , 2
        DO NUCO = 1 , NBCOU
C 
CD        CALL IMPET('NUCO '//IDPROG ,NUCO)
CD        CALL IMPET('POUR '//DEPLA(TYPDEP)//' '//IDPROG ,TYPDEP)
C 
          CALL NDDLCD(TYPDEP,NUBORD,NUCO,ADLOC1,M(ADDGDP+NBTDDL),NBDDL)
          ADDEP =  DDLU(TYPDEP)+NBTDDU(TYPDEP)
          CALL COPITE( NBDDL , M(ADDGDP+NBTDDL) , M( ADDEP ) )
C 
CD        CALL IMPTET('VALEUR DES DDL A DEPLACEMETS IMPOSES'//
CD               ' DANS ' //IDPROG   ,M(ADDGDP+NBTDDL),1,NBDDL)
C 
          NBTDDU(TYPDEP)        = NBTDDU(TYPDEP)+NBDDL
          NBTDDL                = NBTDDL+NBDDL
C 
        END DO
      END DO
C 
C     Creation du tableau des developpements des ddl a
C     deplacement imposes et verification de l'adresse de depart
C   
      CALL GESTEN('DDL-DPI-FI',NBTDDL,TEST1)
C 
D     CALL IMPTET(' DDL POUR LES DEPLACEMENTS DANS '//IDPROG
D    &           ,M(ADDGDP),1,NBTDDL)
C 
C     Creation du tableau des developpements des numeros de ddl correspondant
C     a U imposes par fichier sur le bord exterieur
C 
      CALL GESTEN('DDL-UI-FI4',NBTDDU(1),TEST1)
      CALL COPITE( NBTDDU(1) , M(DDLU(1)) , M( TEST1) )
C 
D      CALL IMPTET(' DDL U IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
D    &           ,M(TEST1),1,NBTDDU(1)  )
C 
C     Creation du tableau des developpements des numeros de ddl correspondant
C     a V imposes par fichier sur le bord exterieur
C 
      CALL GESTEN('DDL-VI-FI4',NBTDDU(2),TEST1)
      CALL COPITE( NBTDDU(2) , M(DDLU(2)) , M( TEST1) )
C 
D     CALL IMPTET(' DDL V IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
D    &           ,M(TEST1),1,NBTDDU(2)  )
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Traitement des deplacements donnes par fichier en flexion
C     Determination des numeros de ddl concernes sur le bord interieur
C     Raccord a une solution plaque en deplacement en flexion
C 
      SUBROUTINE TDEFIF
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'cominc_visu.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  ADLOC1, ADDGDP
      INTEGER  NUCO, TYPDEP, NBDDL, NBTDDL
      INTEGER  NUBORD, TEST1, REST, MOD
      INTEGER  AM2LC, ADM2LC
C 
C     Pour les deplacaments imposes pour U et V
C 
      INTEGER NBTDDU(5), DDLU(5), ADDEP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TDEFIF')
C 
C     include 'identi.h'
C 
C -----------------------------------------------------------------------
      CHARACTER*6      CPLAST(6) , IPLAST(2), CRITER(2)
      CHARACTER*3      CENDOM(3) , ENDOMI(3)
      CHARACTER*10     ENDDIF(3)
      CHARACTER*9      BORD(4)
      CHARACTER*3      CONTRA(6), DEFORM(6)
      CHARACTER*4      DCONTR(6)
      CHARACTER*4      DDEFOR(6)
      CHARACTER*1      DEPLA(3)
C 
      BORD(1)    = 'inferieur'
      BORD(2)    = 'interieur'
      BORD(3)    = 'superieur'
      BORD(4)    = 'exterieur'
C 
      DEPLA(1)   = 'u'
      DEPLA(2)   = 'v'
      DEPLA(3)   = 'w'
C 
C  SGU 05/09/13, VISORT pas clair
C
      VISORT=.TRUE.
      VISUXY=.FALSE.
      IF (VISORT) THEN
C 
         CONTRA(1)  = 'C11'
         CONTRA(2)  = 'C22'
         CONTRA(3)  = 'C12'
         CONTRA(4)  = 'C23'
         CONTRA(5)  = 'C13'
         CONTRA(6)  = 'C33'
C 
         DEFORM(1)  = 'E11'
         DEFORM(2)  = 'E22'
         DEFORM(3)  = 'E12'
         DEFORM(4)  = 'E23'
         DEFORM(5)  = 'E13'
         DEFORM(6)  = 'E33'
C 
         DCONTR(1)  = 'DC11'
         DCONTR(2)  = 'DC22'
         DCONTR(3)  = 'DC12'
         DCONTR(4)  = 'DC23'
         DCONTR(5)  = 'DC13'
         DCONTR(6)  = 'DC33'
C 
         DDEFOR(1)  = 'DE11'
         DDEFOR(2)  = 'DE22'
         DDEFOR(3)  = 'DE12'
         DDEFOR(4)  = 'DE23'
         DDEFOR(5)  = 'DE13'
         DDEFOR(6)  = 'DE33'
C 
      ELSE IF (VISUXY) THEN
C 
         CONTRA(1)  = 'Cxx'
         CONTRA(2)  = 'Cyy'
         CONTRA(3)  = 'Cxy'
         CONTRA(4)  = 'Cyz'
         CONTRA(5)  = 'Cxz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Exx'
         DEFORM(2)  = 'Eyy'
         DEFORM(3)  = 'Exy'
         DEFORM(4)  = 'Eyz'
         DEFORM(5)  = 'Exz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCxx'
         DCONTR(2)  = 'DCyy'
         DCONTR(3)  = 'DCxy'
         DCONTR(4)  = 'DCyz'
         DCONTR(5)  = 'DCxz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DExx'
         DDEFOR(2)  = 'DEyy'
         DDEFOR(3)  = 'DExy'
         DDEFOR(4)  = 'DEyz'
         DDEFOR(5)  = 'DExz'
         DDEFOR(6)  = 'DEzz'
C 
      ELSE
C 
         CONTRA(1)  = 'Crr'
         CONTRA(2)  = 'C00'
         CONTRA(3)  = 'Cr0'
         CONTRA(4)  = 'C0z'
         CONTRA(5)  = 'Crz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Err'
         DEFORM(2)  = 'E00'
         DEFORM(3)  = 'Er0'
         DEFORM(4)  = 'E0z'
         DEFORM(5)  = 'Erz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCrr'
         DCONTR(2)  = 'DC00'
         DCONTR(3)  = 'DCr0'
         DCONTR(4)  = 'DC0z'
         DCONTR(5)  = 'DCrz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DErr'
         DDEFOR(2)  = 'DE00'
         DDEFOR(3)  = 'DEr0'
         DDEFOR(4)  = 'DE0z'
         DDEFOR(5)  = 'DErz'
         DDEFOR(6)  = 'DEzz'
C 
      ENDIF
C 
      CPLAST(1)  = 'epsp11' 
      CPLAST(2)  = 'epsp22' 
      CPLAST(3)  = 'epsp12' 
      CPLAST(4)  = 'epsp23' 
      CPLAST(5)  = 'epsp13' 
      CPLAST(6)  = 'epsp33' 
C 
      CENDOM(1)   = 'dfi'
      CENDOM(2)   = 'dps' 
      CENDOM(3)   = 'dpt'
C 
      IPLAST(1)  = 'SAUTP1' 
      IPLAST(2)  = 'SAUTP2' 
C 
      ENDOMI(1)   = 'di1' 
      ENDOMI(2)   = 'di2'
      ENDOMI(3)   = 'di3'
C 
      ENDDIF(1)   = 'di1-di1ini'
      ENDDIF(1)   = 'di2-di2ini'
      ENDDIF(1)   = 'di3-di3ini'
c -
      CRITER(1)  = 'CRIT-S'
      CRITER(2)  = 'CRIT-N'
C -----------------------------------------------------------------------
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Adresse du tableau de localisation des numeros
C     de noeuds ranges par ordre croissant par elements
C     ranges couche par couche
C -----------------------------------------------------------------------
      CALL ADTBM ('TLOCN1    ', ADLOC1)
C 
C     POUR REMPLIR LE TABLEAU DDL-DPI-FI
C 
      CALL DEBUEN( ADDGDP )
C 
C     POUR REMPLIR LES TABLEAUX 'DDL-UI-FI4' , 'DDL-VI-FI4'
C 
      CALL POUSME( 5*NDDL ,  DDLU(1) )
C 
      DDLU(2) = DDLU(1)+NDDL
      DDLU(3) = DDLU(2)+NDDL
      DDLU(4) = DDLU(3)+NDDL
      DDLU(5) = DDLU(4)+NDDL
C 
C     Interet : on impose la forme du deplacement sur tout un bord
C     Exemple : raccord avec une solution plaque
C 
      NBTDDL    = 0
      NBTDDU(1) = 0
      NBTDDU(2) = 0
C 
C     w est l'objet d'un traitement particulier puisque on ne l'impose
C     que sur la ligne moyenne
C 
      NBTDDU(3) = 0
      NBTDDU(4) = 0
      NBTDDU(5) = 0
C 
C     Les deplacements sont imposes sur le bord exterieur
C 
      NUBORD = 4
C 
      DO TYPDEP = 1 , 2
         DO NUCO = 1 , NBCOU
C 
CD         CALL IMPET('NUCO '//IDPROG ,NUCO)
CD         CALL IMPET('POUR '//DEPLA(TYPDEP)//' '//IDPROG ,NUCO)
C 
	   CALL NDDLCD(
     &           TYPDEP,NUBORD,NUCO,ADLOC1,M(ADDGDP+NBTDDL),NBDDL)
           ADDEP =  DDLU(TYPDEP)+NBTDDU(TYPDEP)
           CALL COPITE( NBDDL , M(ADDGDP+NBTDDL) , M( ADDEP ) )
C 
CD         CALL IMPTET('VALEUR DES DDL A DEPLACEMETS IMPOSES'//
CD   &           ' DANS ' //IDPROG   ,M(ADDGDP+NBTDDL),1,NBDDL)
C 
           NBTDDU(TYPDEP)        = NBTDDU(TYPDEP)+NBDDL
           NBTDDL                = NBTDDL+NBDDL
C 
         END DO
      END DO
C 
C     Pour w
      REST   = MOD( NBCOU , 2)
      TYPDEP = 3
      NUCO = NBCOU/2
C 
CD    CALL IMPET('NUCO '//IDPROG ,NUCO)
CD    CALL IMPET('POUR '//DEPLA(TYPDEP)//' '//IDPROG ,TYPDEP)
C 
C      CALL NDDLCD(TYPDEP,NUBORD,NUCO,ADLOC1,M(ADDGDP+NBTDDL),NBDDL)
      IF (NUCO .GE. 1) THEN
        CALL NDDLCD(TYPDEP,NUBORD,NUCO,ADLOC1,M(ADDGDP+NBTDDL),NBDDL)
C 
        ADDEP =  DDLU(TYPDEP)+NBTDDU(TYPDEP)
C 
        M( ADDEP )       = M(ADDGDP+NBTDDL+4)
        M(ADDGDP+NBTDDL) = M(ADDGDP+NBTDDL+4)
C 
CD    CALL IMPET('VALEUR DU DDL DE W A DEPLACEMETS IMPOSES'//
CD                 ' DANS ' //IDPROG   ,M(ADDEP))
C 
        NBTDDU(TYPDEP)        = NBTDDU(TYPDEP)+1
        NBTDDL                = NBTDDL+1
      END IF
C 
      IF ( REST .EQ. 1 ) THEN
         NUCO = NBCOU/2+1
C 
CD       CALL IMPET('NUCO '//IDPROG ,NUCO)
CD       CALL IMPET('POUR '//DEPLA(TYPDEP)//' '//IDPROG ,NUCO)
C 
        CALL NDDLCD(TYPDEP,NUBORD,NUCO,ADLOC1,M(ADDGDP+NBTDDL),NBDDL)
        ADDEP =  DDLU(TYPDEP)+NBTDDU(TYPDEP)
        M( ADDEP )        = M(ADDGDP+NBTDDL+4)
        M(ADDGDP+NBTDDL)  = M(ADDGDP+NBTDDL+4)
C 
CD       CALL IMPET('VALEUR DU DDL DE W A DEPLACEMETS IMPOSES'//
CD                  ' DANS ' //IDPROG   ,M(ADDEP))
C 
        NBTDDU(TYPDEP)        = NBTDDU(TYPDEP)+1
        NBTDDL                = NBTDDL+1
      END IF
C    
C  
C     Creation du tableau des developpements des ddl a
C     deplacements imposes et verification de l'adresse de depart
C  
      CALL GESTEN('DDL-DPI-FI',NBTDDL,TEST1)
C 
       CALL IMPTET(' DDL POUR LES DEPLACEMENTS DANS '//IDPROG
     &          ,M(ADDGDP),1,NBTDDL)
C 
C     Creation du tableau des developpements des numeros de ddl correspondant
C     a U imposes par fichier sur le bord exterieur
C 
      CALL GESTEN('DDL-UI-FI4',NBTDDU(1),TEST1)
      CALL COPITE( NBTDDU(1) , M(DDLU(1)) , M( TEST1) )
C 
       CALL IMPTET(' DDL U IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
     &       ,M(ADDGDP),1,NBTDDU(1)  )
C 
C     Creation du tableau des developpements des numeros de ddl correspondant
C     a V imposes par fichier sur le bord exterieur
C 
      CALL GESTEN('DDL-VI-FI4',NBTDDU(2),TEST1)
      CALL COPITE( NBTDDU(2) , M(DDLU(2)) , M( TEST1) )
C 
      CALL IMPTET(' DDL V IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
     &             ,M(TEST1),1,NBTDDU(2)  )
C 
      CALL GESTEN('DDL-WI-FI4',NBTDDU(3),TEST1)
      CALL COPITE( NBTDDU(3) , M(DDLU(3)) , M( TEST1) )
C 
      CALL IMPTET(' DDL W IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
     &           ,M(TEST1),1,NBTDDU(3)  )
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Modification des deplacements donnes par fichier en FLEXION.
C  
C     Raccord a une solution plaque en deplacement en FLEXION.

      SUBROUTINE MDEFIF
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER U, V, W, RR, R0, TYP(5)
      INTEGER DDLU(3), LONGU(3), VERLON(3)
      INTEGER MOD, DBCOTE
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  ADCOTE, ADCOTL
      INTEGER  ADLOC1, I, NUCO, TYPDEP, NBDDL
      INTEGER  NUBORD, J
      INTEGER  NFT, NUDDL(6), REST
      INTEGER  APRECI, ADEFFO, EFFOLC, EFFORT, DBDDL
      INTEGER  AM2LC, ADM2LC
C 
      DOUBLE PRECISION VALDEP
      DOUBLE PRECISION MULDIA, SUP, INF
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MDEFIF')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB( AM2LC , ADM2LC )
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('HAUT-NOEUD', ADCOTE)
C 
C     'U-T-IMP-FI' Tableau NTETA des valeurs developpees du deplacement
C                       radial
C 
C     'V-T-IMP-FI' Tableau NTETA des valeurs developpees du deplacement
C                       orthoradial
C 
C     'W-T-IMP-FI' Tableau NTETA des valeurs developpees du deplacement
C                       normal
C 
C     'RRT-IM-FI' Tableau NTETA des valeurs developpees de la rotation
C                       radiale
C 
C     'R0T-T-IM-FI' Tableau NTETA des valeurs developpees de la rotation
C                       orthoradiale
C 
      CALL ADTBDM('U-T-IMP-FI',  U  )
      TYP(1) = U  + NTDSFG
      CALL IMPTDT ('TYP(1)', DM(TYP(1)), 1, NTDSFG)
C 
      CALL ADTBDM('V-T-IMP-FI',  V  )
      TYP(2) = V  + NTDSFG
      CALL IMPTDT ('TYP(2)', DM(TYP(2)), 1, NTDSFG)
C 
      CALL ADTBDM('W-T-IMP-FI',  W  )
      TYP(3) = W  + NTDSFG
C 
      CALL ADTBDM('RRT-IMP-FI',  RR )
      TYP(4) = RR + NTDSFG
C 
      CALL ADTBDM('R0T-IMP-FI',  R0 )
      TYP(5) = R0 + NTDSFG
C 
C     Adresse du tableau de localisation des numeros
C     de noeuds ranges par ordre croissant par elements
C     ranges couche par couche
C 
      CALL ADTBM ('TLOCN1    ', ADLOC1)
C 
C     Tableau des developpements des numeros de ddl correspondant
C     a U imposes par fichier sur le bord exterieur
C 
      CALL INFOEN('DDL-UI-FI4', DDLU(1) ,LONGU(1) )
C 
      CALL IMPTET(' DDL U IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
     &             ,M(DDLU(1) ),1,LONGU(1)  )
C 
C     Tableau des developpements des numeros de ddl correspondant
C     a V imposes par fichier sur le bord exterieur
C 
      CALL INFOEN('DDL-VI-FI4',DDLU(2),LONGU(2))
C 
      CALL IMPTET(' DDL v IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
     &             ,M(DDLU(2) ),1,LONGU(2)  )
C 
C     Tableau des developpements des numeros de ddl correspondant
C     a W imposes par fichier sur le bord exterieur
C 
      CALL INFOEN('DDL-WI-FI4',DDLU(3),LONGU(3))
C 
      CALL IMPTET(' DDL W IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
     &             ,M(DDLU(3) ),1,LONGU(3)  )
C 
C     POUR VERIFIER LES LONGUEURS LUES
C 
      VERLON(1) = 0
      VERLON(2) = 0
      VERLON(3) = 0
C 
C -----------------------------------------------------------------------
C     Sequence de modification des efforts pour les deplacements imposes
C -----------------------------------------------------------------------
C 
C     Pour multiplier les deplacements par les termes
C     diagonaux bloques
C 
      CALL ADTBDM('PRECISIONS',APRECI)
      CALL IMPTDT ('PRECISIONS, ', DM(APRECI), 1, NBMAT)
C 
      APRECI = APRECI+NTDSFG
C 
      CALL GSPOUD( 100 , EFFOLC)
C 
C     adeffo est l'adresse de depart du tableau des efforts pour
C     la nft-ieme fonction du temps
C 
      CALL ADTBDM('MAT-EFFORT',ADEFFO)
C 
      NFT    = 1
      EFFORT = ADEFFO+(NFT-1)*NDDL*NBMAT
      NBDDL  = 3
C 
C     POUR U
C 
      TYPDEP= 1
      DBDDL = 0
      ADCOTL = ADCOTE
C 
      DO NUCO =  1 , NBCOU
        DO J = -NTDSFG , NTDSFG
          MULDIA = DM( APRECI+J)
          DBCOTE = ADCOTL
          DO I = 1 , NBDDL
            VALDEP = DM( TYP(TYPDEP)+J)+DM( TYP(5)+J)*DM(DBCOTE)
            DM( EFFOLC+I-1) = MULDIA*VALDEP
            DBCOTE = DBCOTE+1
          END DO
	  DO I =1, NBDDL
CD          print*, 'DM(EFFOLC+I-1)', DM(EFFOLC+I-1)
CD	    print*, 'BORD 1, DDLU 1 ', I, ', ',
CD   &                M(DDLU(1)+DBDDL+I-1)
          END DO
CD	  print*, ' '
          CALL ASVEFI( EFFORT, J, M(DDLU(TYPDEP)+DBDDL ), NBDDL,
     &                 DM(EFFOLC) )
        END DO
        VERLON(TYPDEP) = VERLON(TYPDEP) + NBDDL
        DBDDL  = DBDDL  + NBDDL
        ADCOTL = ADCOTL + NBDDL
      END DO
C 
CD    CALL TESTEN ( LONGU(TYPDEP) , VERLON(TYPDEP) , IDPROG )
C 
C     POUR V 
C 
      TYPDEP= 2
      DBDDL = 0
      ADCOTL = ADCOTE
C 
      DO NUCO =  1 , NBCOU
        DO J = -NTDSFG , NTDSFG
          MULDIA = DM( APRECI+J)
          DBCOTE = ADCOTL
          DO I = 1 , NBDDL
            VALDEP = DM( TYP(TYPDEP)+J)-DM( TYP(4)+J)*DM(DBCOTE)
            DM( EFFOLC+I-1) = MULDIA*VALDEP
            DBCOTE = DBCOTE+1
          END DO
	  DO I =1, NBDDL
CD          print*, 'DM(EFFOLC+I-1)', DM(EFFOLC+I-1)
CD          print*, 'BORD 1, DDLU 2 ', I, ', ',
CD   &                M(DDLU(2)+DBDDL+I-1)
          END DO
CD	  print*, ' '
          CALL ASVEFI( EFFORT, J, M(DDLU(TYPDEP)+DBDDL), NBDDL,
     &                 DM(EFFOLC) )
        END DO
        VERLON(TYPDEP) = VERLON(TYPDEP) + NBDDL
        DBDDL  = DBDDL  + NBDDL
        ADCOTL = ADCOTL + NBDDL
      END DO
C 
CD    CALL TESTEN ( LONGU(TYPDEP) , VERLON(TYPDEP) , IDPROG )
C 
C     POUR W
C 
      TYPDEP= 3
      DBDDL = 0
      REST = MOD( NBCOU ,2 )
      NUCO = NBCOU/2
C 
      DO J = -NTDSFG , NTDSFG
        MULDIA = DM( APRECI+J)
        VALDEP = DM( TYP(TYPDEP)+J)
        DM( EFFOLC) = MULDIA*VALDEP
        CALL ASVEFI( EFFORT, J, M(DDLU(TYPDEP)+DBDDL ), 1, DM(EFFOLC) )
      END DO
C 
      VERLON(TYPDEP) = VERLON(TYPDEP) + 1
C 
      DBDDL  = DBDDL+1
C 
      IF ( REST . EQ. 1 ) THEN
        NUCO = NBCOU/2+1
        DO J = -NTDSFG , NTDSFG
          MULDIA = DM( APRECI+J)
          VALDEP = DM( TYP(TYPDEP)+J)
          DM( EFFOLC) = MULDIA*VALDEP
          CALL ASVEFI( EFFORT, J, M(DDLU(TYPDEP)+DBDDL ), 1, DM(EFFOLC))
        END DO
        VERLON(TYPDEP) = VERLON(TYPDEP) + 1
        DBDDL  = DBDDL+1
      END IF
C 
CD    CALL TESTEN ( LONGU(TYPDEP) , VERLON(TYPDEP) , IDPROG )
C 
C    DR  AJOUT le 16/2/96 : Hypothese de Kirchoff sur le bord de la plaque
C 
C    POUR W, r
C 
      TYPDEP = 3
      NUBORD = 4
C 
      DO NUCO =  1 , NBCOU
        CALL NDDLCR(TYPDEP, NUBORD, NUCO, ADLOC1, NUDDL, NBDDL)
        IF (NBDDL .NE. 3) THEN
          CALL ERREUD(0, ' ERREUR FATALE DANS '//IDPROG)
        END IF
        DO J = -NTDSFG , NTDSFG
          MULDIA = DM( APRECI+J)
          DO I = 1 , NBDDL
            DM( EFFOLC+I-1) = MULDIA*DM( TYP(5)+J )
          END DO
C 
C     DR          CALL ASVEFI( EFFORT, J, NUDDL, NBDDL, DM(EFFOLC) )
C 
        END DO
      END DO
C 
      CALL SOPOUB( AM2LC , ADM2LC )
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Modification des deplacements donnes par fichier en TRACTION.
C     Deplacements traites par penalisation.
C 
C     Raccord a une solution plaque en deplacement en TRACTION.
C 
      SUBROUTINE MDEFIT
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER U, V, TYP(2)
      INTEGER DDLU(2), LONGU(2), VERLON(2)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  ADLOC1
      INTEGER  I, NUCO, TYPDEP, NBDDL
      INTEGER  J, NFT
      INTEGER  APRECI, ADEFFO, EFFOLC, EFFORT, DBDDL
      INTEGER  AM2LC, ADM2LC
C 
      DOUBLE PRECISION MULDIA
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MDEFIT')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     'U-T-IMP-FI' Tableau NTETA des valeurs developpees du deplacement radial
C 
      CALL ADTBDM ('U-T-IMP-FI', U)
      TYP(1) = U+NTDSFG
C 
C     'V-T-IMP-FI' Tableau NTETA des valeurs developpees du deplacement orthoradial
C 
      CALL ADTBDM ('V-T-IMP-FI', V)
      TYP(2) = V+NTDSFG
C 
C     Adresse du tableau de localisation des numeros
C     de noeuds ranges par ordre croissant par elements
C     ranges couche par couche
C 
      CALL ADTBM ('TLOCN1    ', ADLOC1)
C 
C     Tableau des developpements des numeros de ddl correspondant
C     a U imposes par fichier sur le bord exterieur
C 
      CALL ADTBM ('DDL-UI-FI4',DDLU(1) )
      CALL LONGEN ('DDL-UI-FI4',LONGU(1))
C 
CD    CALL IMPTET(' DDL U IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
CD                ,M(DDLU(1) ),1,LONGU(1)  )
C 
C     Tableau des developpements des numeros de ddl correspondant
C     a V imposes par fichier sur le bord exterieur
C 
      CALL ADTBM ('DDL-VI-FI4',DDLU(2))
      CALL LONGEN ('DDL-VI-FI4',LONGU(2))
C 
CD    CALL IMPTET(' DDL v IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG
CD                ,M(DDLU(2) ),1,LONGU(2)  )
C 
C     POUR VERIFIER LES LONGUEURS LUES
C 
      VERLON(1) = 0
      VERLON(2) = 0
C 
C -----------------------------------------------------------------------
C     Sequence de modification des efforts pour les deplacements imposes
C -----------------------------------------------------------------------
C 
C     Pour multiplier les deplacements par les termes
C     diagonaux bloques
C 
      CALL ADTBDM ('PRECISIONS', APRECI)
C 
      APRECI = APRECI+NTDSFG
C 
      CALL GSPOUD (10, EFFOLC)
C 
C     adeffo est l'adresse de depart du tableau des efforts pour
C     la nft-ieme fonction du temps
C 
      CALL ADTBDM ('MAT-EFFORT', ADEFFO)
C 
      NFT   = 1
      EFFORT= ADEFFO+(NFT-1)*NDDL*NBMAT
      NBDDL = 3
C 
      DO TYPDEP= 1 , 2
        DBDDL = 0
        DO NUCO =  1 , NBCOU
          DO J = -NTDSFG , NTDSFG
            MULDIA = DM(APRECI+J)
            DO I = 1 , NBDDL
              DM(EFFOLC+I-1) = MULDIA*DM(TYP(TYPDEP)+J)
            END DO
            CALL ASVEFI (EFFORT, J, M(DDLU(TYPDEP)+DBDDL ), NBDDL,
     &                   DM(EFFOLC))
          END DO
          VERLON(TYPDEP) = VERLON(TYPDEP) + NBDDL
          DBDDL = DBDDL+NBDDL
        END DO
C 
CD      CALL TESTEN ( LONGU(TYPDEP) , VERLON(TYPDEP) , IDPROG )
C 
      END DO
C 
CD    CALL GSPOUD( NDDL*(NBMAT+NTETA) , VEREFF )
C 
CD    VERDEP = VEREFF+NDDL*NBMAT
C 
CD    DO I = -NTDSFG , NTDSFG
CD      DEBEFF = ADEFFO+(I+NTDSFG)*NDDL
CD      CALL MUMARE( 1.D0/DM(APRECI+I), NDDL, DM(DEBEFF), DM(VEREFF) )
CD      VEREFF = VEREFF+NDDL
CD    END DO
C 
CD    VEREFF =VEREFF-NBMAT*NDDL
C 
CD    CALL VRTSM ( 1 , DM(VEREFF) , DM(VERDEP)  )
C 
CD    CALL MESSAO(' APPEL A VEDIFT DANS '//IDPROG )
CD    CALL VEDIFT ( 1 , DM(VERDEP) , NTETA*NDDL )
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
C     VERification des Deplacements Imposes  par Fichier en Traction
C 
C     Pour vefifier si la solution soit developpee soit reelle
C     est conforme aux valeurs imposees.
C 
C     On envoie comme arguments :
C 
C     E ...... TYPE    0 = DEVELOPPEE
C     E ...... TYPE    1 = REELLE
C     E ...... DEPSOL  Solution testee
C     E ...... LONG    Longueur de DEPSOL => (nddl,nbmat) si 0
C                                         => (nteta,nddl) si 1
C 
      SUBROUTINE VEDIFT (TYPE, DEPSOL, LONG)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'cominc_visu.h'
C 
      INTEGER           TYPE, LONG
      DOUBLE PRECISION  DEPSOL(LONG)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
C     Pour les deplacaments imposes  pour U et V
C 
      INTEGER  DDLU(2), U, V, TYP(2), LONGU(2)
C 
      INTEGER TYPDEP, NUDDL, NUMDDL
      INTEGER PLAC, J
      INTEGER AM2LC, ADM2LC
      INTEGER VALDEP, DEPIMP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VEDIFT')
C 
C     include 'identi.h'
C 
C -----------------------------------------------------------------------
      CHARACTER*6      CPLAST(6) , IPLAST(2), CRITER(2)
      CHARACTER*3      CENDOM(3) , ENDOMI(3)
      CHARACTER*10     ENDDIF(3)
      CHARACTER*9      BORD(4)
      CHARACTER*3      CONTRA(6), DEFORM(6)
      CHARACTER*4      DCONTR(6)
      CHARACTER*4      DDEFOR(6)
      CHARACTER*1      DEPLA(3)
C 
      BORD(1)    = 'inferieur'
      BORD(2)    = 'interieur'
      BORD(3)    = 'superieur'
      BORD(4)    = 'exterieur'
C 
      DEPLA(1)   = 'u'
      DEPLA(2)   = 'v'
      DEPLA(3)   = 'w'
C 
C  SGU 05/09/13, VISORT pas clair
C
      VISORT=.TRUE.
      VISUXY=.FALSE.
      IF (VISORT) THEN
C 
         CONTRA(1)  = 'C11'
         CONTRA(2)  = 'C22'
         CONTRA(3)  = 'C12'
         CONTRA(4)  = 'C23'
         CONTRA(5)  = 'C13'
         CONTRA(6)  = 'C33'
C 
         DEFORM(1)  = 'E11'
         DEFORM(2)  = 'E22'
         DEFORM(3)  = 'E12'
         DEFORM(4)  = 'E23'
         DEFORM(5)  = 'E13'
         DEFORM(6)  = 'E33'
C 
         DCONTR(1)  = 'DC11'
         DCONTR(2)  = 'DC22'
         DCONTR(3)  = 'DC12'
         DCONTR(4)  = 'DC23'
         DCONTR(5)  = 'DC13'
         DCONTR(6)  = 'DC33'
C 
         DDEFOR(1)  = 'DE11'
         DDEFOR(2)  = 'DE22'
         DDEFOR(3)  = 'DE12'
         DDEFOR(4)  = 'DE23'
         DDEFOR(5)  = 'DE13'
         DDEFOR(6)  = 'DE33'
C 
      ELSE IF (VISUXY) THEN
C 
         CONTRA(1)  = 'Cxx'
         CONTRA(2)  = 'Cyy'
         CONTRA(3)  = 'Cxy'
         CONTRA(4)  = 'Cyz'
         CONTRA(5)  = 'Cxz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Exx'
         DEFORM(2)  = 'Eyy'
         DEFORM(3)  = 'Exy'
         DEFORM(4)  = 'Eyz'
         DEFORM(5)  = 'Exz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCxx'
         DCONTR(2)  = 'DCyy'
         DCONTR(3)  = 'DCxy'
         DCONTR(4)  = 'DCyz'
         DCONTR(5)  = 'DCxz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DExx'
         DDEFOR(2)  = 'DEyy'
         DDEFOR(3)  = 'DExy'
         DDEFOR(4)  = 'DEyz'
         DDEFOR(5)  = 'DExz'
         DDEFOR(6)  = 'DEzz'
C 
      ELSE
C 
         CONTRA(1)  = 'Crr'
         CONTRA(2)  = 'C00'
         CONTRA(3)  = 'Cr0'
         CONTRA(4)  = 'C0z'
         CONTRA(5)  = 'Crz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Err'
         DEFORM(2)  = 'E00'
         DEFORM(3)  = 'Er0'
         DEFORM(4)  = 'E0z'
         DEFORM(5)  = 'Erz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCrr'
         DCONTR(2)  = 'DC00'
         DCONTR(3)  = 'DCr0'
         DCONTR(4)  = 'DC0z'
         DCONTR(5)  = 'DCrz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DErr'
         DDEFOR(2)  = 'DE00'
         DDEFOR(3)  = 'DEr0'
         DDEFOR(4)  = 'DE0z'
         DDEFOR(5)  = 'DErz'
         DDEFOR(6)  = 'DEzz'
C 
      ENDIF
C 
      CPLAST(1)  = 'epsp11' 
      CPLAST(2)  = 'epsp22' 
      CPLAST(3)  = 'epsp12' 
      CPLAST(4)  = 'epsp23' 
      CPLAST(5)  = 'epsp13' 
      CPLAST(6)  = 'epsp33' 
C 
      CENDOM(1)   = 'dfi'
      CENDOM(2)   = 'dps' 
      CENDOM(3)   = 'dpt'
C 
      IPLAST(1)  = 'SAUTP1' 
      IPLAST(2)  = 'SAUTP2' 
C 
      ENDOMI(1)   = 'di1' 
      ENDOMI(2)   = 'di2'
      ENDOMI(3)   = 'di3'
C 
      ENDDIF(1)   = 'di1-di1ini'
      ENDDIF(1)   = 'di2-di2ini'
      ENDDIF(1)   = 'di3-di3ini'
c -
      CRITER(1)  = 'CRIT-S'
      CRITER(2)  = 'CRIT-N'
C -----------------------------------------------------------------------
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL GSPOUD (2*LONG, VALDEP)
      VALDEP = VALDEP-1
      DEPIMP = VALDEP+LONG
C 
      IF (TYPE .EQ. 0) THEN
C 
        CALL MESSAO
     &    ('VERIFICATION DES DEPLACEMENTS DEVELOPPES DANS '//IDPROG)
C 
        CALL ADTBM ('DDL-UI-FI4', DDLU(1))
        CALL LONGEN ('DDL-UI-FI4', LONGU(1))
C 
CD      CALL IMPTEN ('DDL U IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG,
CD                    M(DDLU(1)), 1, LONGU(1))
C 
C       Creation du tableau des developpements des numeros de ddl correspondant
C       a V imposes par fichier sur le bord exterieur.
C 
        CALL ADTBM ('DDL-VI-FI4', DDLU(2))
        CALL LONGEN ('DDL-VI-FI4', LONGU(2))
C 
CD      CALL IMPTEN ('DDL V IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG,
CD                    M(DDLU(2)), 1, LONGU(2) )
C 
C       U-T-IMP-FI Tableau NTETA des valeurs developpees du deplacement radial
C 
        CALL ADTBDM ('U-T-IMP-FI', U)
        TYP(1) = U+NTDSFG
C 
C       V-T-IMP-FI Tableau NTETA des valeurs developpes du deplacement orthoradial
C 
        CALL ADTBDM ('V-T-IMP-FI', V)
        TYP(2) = V+NTDSFG
C 
C       Sequence de verification des deplacements imposes par fichier dans le cas developpe
C 
        IF (LONG .NE. NDDL*NBMAT) THEN
          CALL IMPET ('POUR LE CAS DEVELOPPE : TYPE = ', TYPE )
          CALL ERREUD (0 , 'MAUVAISE LONGUEUR DE DEPSOL DANS '//IDPROG)
        END IF
C 
        DO TYPDEP= 1, 2
C 
          DO NUMDDL = 0, LONGU(TYPDEP)-1
C 
            NUDDL =  M( DDLU(TYPDEP) + NUMDDL )
C 
            DO J = -NTDSFG , NTDSFG
              PLAC = (J+NTDSFG)*NDDL +NUDDL
              DM(DEPIMP+J) = DM(TYP(TYPDEP)+J)
              DM(VALDEP+J) = DEPSOL(PLAC)
              IF (DABS(DM(DEPIMP+J)-DEPSOL(PLAC)) .GT. 1.D -8) THEN
                CALL MESSAO ('!!!!!!!! PROBLEME DANS '//IDPROG )
                CALL IMPET ('POUR LE NUMERO DE DDL DE TYPE '
     &                      //DEPLA(TYPDEP)//' '//BORD(4)//
     &                      ' DANS '//IDPROG, NUDDL )
                CALL IMPET ('POUR LE DEVELOPPEMENT   ', J)
                CALL IMPDT ('LA VALEUR IMPOSEE EST   ', DM(DEPIMP+J))
                CALL IMPDT ('LA VALEUR DE DEPSOL EST ', DM(VALDEP+J))
              END IF
            END DO
C 
CD            CALL IMPEN ('POUR LE NUMERO DE DDL DE TYPE '
CD                       //DEPLA(TYPDEP)//' '//BORD(4)//
CD                        ' DANS '//IDPROG, NUDDL)
CD            CALL IMPTDN ('DEP DEVELOPPE CALCULE '//IDPROG,
CD                          DM(VALDEP+1), NTDSFG, 1)
CD            CALL IMPTDN ('DEP DEVELOPPE IMPOSE  '//IDPROG ,
CD                          DM(DEPIMP+1), NTDSFG, 1)
C 
          END DO
C 
        END DO
C 
      END IF
C 
      IF (TYPE .EQ. 1) THEN
C 
        CALL MESSAO
     &  ('VERIFICATION DES DEPLACEMENTS REELS DANS'//IDPROG)
C 
        CALL ADTBM ('DDL-UI-FI4', DDLU(1))
        CALL LONGEN ('DDL-UI-FI4', LONGU(1))
C 
CD      CALL IMPTEN ('DDL U IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG,
CD                    M(DDLU(1)), 1, LONGU(1))
C 
C       Creation du tableau des developpements des numeros de ddl correspondant
C       a V imposes par fichier sur le bord exterieur
C 
        CALL ADTBM ('DDL-VI-FI4', DDLU(2))
        CALL LONGEN ('DDL-VI-FI4', LONGU(2))
C 
CD      CALL IMPTEN ('DDL V IMPOSE SUR LE BORD EXTERIEUR DANS '//IDPROG,
CD                    M(DDLU(2)), 1, LONGU(2))
C 
C       U-R-IMP-FI Tableau NTETA des valeurs developpees du deplacement radial
C 
        CALL ADTBDM ('U-R-IMP-FI', U)
        TYP(1) = U-1
C 
C       V-R-IMP-FI Tableau NTETA des valeurs developpes du deplacement orthoradial
C 
        CALL ADTBDM('V-R-IMP-FI',  V )
        TYP(2) = V-1
C 
C        Sequence de verification des  deplacements imposes  par fichier dans le cas developpe
C 
        IF (LONG .NE. NDDL*NTETA) THEN
          CALL IMPET ('POUR LE CAS DEVELOPPE : TYPE =' , TYPE)
          CALL ERREUD (0, 'MAUVAISE LONGUEUR DE DEPSOL DANS '//IDPROG)
        END IF
C 
        DO TYPDEP= 1, 2
C 
          DO NUMDDL = 0, LONGU(TYPDEP)-1
C 
            NUDDL =  M( DDLU(TYPDEP) + NUMDDL )
            PLAC   = (NUDDL-1)*NTETA
C 
            DO J = 1, NTETA
              DM(DEPIMP+J) = DM(TYP(TYPDEP)+J)
              DM(VALDEP+J)  = DEPSOL(PLAC+J)
              IF (DABS(DM(VALDEP+J) - DM(DEPIMP+J)) .GT. 1.D -6) THEN
                CALL MESSAO ('!!!!!!!! PROBLEME DANS '//IDPROG )
                CALL IMPET ('POUR LE NUMERO DE DDL DE TYPE '
     &                      //DEPLA(TYPDEP)//' '//BORD(4)//
     &                      ' DANS '//IDPROG, NUDDL )
                CALL IMPET ('NUMERO D''ANGLE '//IDPROG, J)
                CALL IMPDT ('LA VALEUR REELLE IMPOSEE EST ',
     &                       DM(DEPIMP+J))
                CALL IMPDT ('LA VALEUR DE REELE CALCULEE EST ',
     &                       DM(VALDEP+J))
              END IF
            END DO
C 
CD          CALL IMPEN('POUR LE NUMERO DE DDL DE TYPE '
CD          //DEPLA(TYPDEP)//' '//BORD(4)//' DANS '//IDPROG
CD          , NUDDL )
CD          CALL IMPTDN( 'DEP REEL CALCULE '//IDPROG ,
CD                 DM(VALDEP+1),   NTETA , 1 )
CD          CALL IMPTDN( 'DEP REEL IMPOSE '//IDPROG  ,
CD                   DM(DEPIMP+1),NTETA , 1 )
C 
          END DO
C 
        END DO
C 
      END IF
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
