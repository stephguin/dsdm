C     On envoie comme arguments :
C 
C     E ........ NC le numero de couches
C 
C     Et on recupere :
C 
C     S ........ EPAIS l'epaisseur de l'element
C     S ........ HAUTEC la hauteur au centre de l'element
C 
      SUBROUTINE VALHAU (NC, EPAIS, HAUTEC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER   NC
C  
      DOUBLE PRECISION EPAIS, HAUTEC
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   A1
C  
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VALHAU')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('HAUTEURS  ',A1)
      HAUTEC   = DM(A1-1+NC)
      CALL ADTBDM('EPAISSEURS',A1)
      EPAIS    = DM(A1-1+NC)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     Le type de deplacement  TYPDEP
C 
C           si TYPDEP=1 ==========> U
C           si TYPDEP=2 ==========> V
C           si TYPDEP=3 ==========> W
C 
C     Le numero de couche
C     Le numero de colonne
C     L'adresse de depart du tableau TLOCN1
C 
C     Et on recupere :
C 
C     NDDLDE    tableau des numeros reels des ddl correspondant aux
C               deplacements stockes croissants.
C 
      SUBROUTINE NDDLC (TYPDEP, NUCOU, NUCOL, TLOCN1, NDDLDE)
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
      INTEGER NDDLDE(12),I,NUCOU,NUCOL,NI,NELAV,TLOCN1,TYPDEP,DECAL
      CHARACTER*1 IDDEPL
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NDDLC ')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF ((NUCOU .GT. NBCOU) .OR. (NUCOU .LT. 1)) THEN
        CALL IMPET ('VALEUR DE NUCOU ', NUCOU)
        CALL ERREUD (0, 'MAUVAIS PASSAGE DE NUMERO DE COUCHE '//IDPROG)
      ENDIF
      IF ((NUCOL .GT. NBCOL) .OR. (NUCOL .LT. 1)) THEN
        CALL IMPET ('VALEUR DE NUCOL ', NUCOL)
        CALL ERREUD (0, 'MAUVAIS PASSAGE DE NUMERO DE COLONNE '//IDPROG)
      ENDIF
      IF (TYPDEP .EQ. 1) THEN
        DECAL   =  1
        IDDEPL  = 'U'
      ELSE IF (TYPDEP .EQ. 2) THEN
        DECAL   =  3
        IDDEPL  = 'V'
      ELSE IF (TYPDEP .EQ. 3) THEN
        DECAL   =  5
        IDDEPL  = 'W'
      ELSE
        CALL ERREUD (0, 'MAUVAIS PASSAGE DE TYPE DE DEPL DANS '//IDPROG)
      ENDIF
      NELAV     = 6*(NBCOL*(NUCOU-1)+(NUCOL-1))
      NELAV     = NELAV+TLOCN1-1
C 
      DO I=1,6
        NI=6*(M(NELAV+I)-1)+DECAL
        NDDLDE(2*(I-1)+1)=NI
        NDDLDE(2*I)=NI+1
      ENDDO
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... NUBORD
C     E ...... NUCO
C 
C     Et on recupere :
C 
C     S ...... RC     rayon au centre
C     S ...... A      demi-longueur
C     S ...... ZC     rayon au centre
C     S ...... B      demi-longueur
C 
      SUBROUTINE CAGEO (NUBORD, NUCO, RC, A, ZC, B)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NUBORD, NUCO
      DOUBLE PRECISION  A, RC, B, ZC
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CAGEO ')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (NUBORD.EQ.1)THEN
        CALL VALRAY( NUCO,RC,A)
        CALL VALHAU( 1 ,B,ZC)
      ELSE IF (NUBORD.EQ.3)THEN
        CALL VALRAY( NUCO,RC,A)
        CALL VALHAU(NBCOU, B , ZC)
      ELSE IF (NUBORD.EQ.2)THEN
        CALL VALRAY( 1 ,RC,A)
        CALL VALHAU(NUCO, B , ZC)
      ELSE IF (NUBORD.EQ.4)THEN
        CALL VALRAY( NBCOL ,RC,A)
        CALL VALHAU(NUCO, B , ZC)
      ENDIF
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E............ NUBORD
C     E............ NUCO
C 
C     Et on recupere :
C 
C     S............ POSIC    position du centre
C     S............ DIMENS   demi-longueur ou hauteur

      SUBROUTINE CAGEOB (NUBORD, NUCO, POSIC, DIMENS)
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
      INTEGER           NUBORD, NUCO
      DOUBLE PRECISION  POSIC , DIMENS
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CAGEOB')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (NUBORD.EQ.1)THEN
        CALL VALRAY( NUCO,POSIC,DIMENS)
      ELSE IF (NUBORD.EQ.3)THEN
        CALL VALRAY( NUCO,POSIC,DIMENS)
      ELSE IF (NUBORD.EQ.2)THEN
        CALL VALHAU(NUCO, DIMENS , POSIC)
      ELSE IF (NUBORD.EQ.4)THEN
        CALL VALHAU(NUCO, DIMENS , POSIC)
      ENDIF
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine remplit le tableau de replacage pour les couches.
C 
C     Le tableau RANCAL-COU qui donne la numerotation de calcul (U, V, W)
C     rangee dans un ordre tel que les numeros reels correspondant dans
C     les matrices stockees profil soient en ordre croissant.
C 
      SUBROUTINE REPLAC
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
C ----------------------------------------------------------------------- 
C     DECLARATION DES PARAMETRES LOCAUX
C     """""""""""""""""""""""""""""""""
C 
      INTEGER M2, R2, X, APM1, ARM1, I, K
      INTEGER ARM2, APM2, DECAL, TYPDEP, APM3, ARM3
      INTEGER AM2LC, ADM2LC
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='REPLAC')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C ----------------------------------------------------------------------- 
C     Creation d'un tableau provisoire pour ranger les numeros des noeuds
C     dans l'ordre correspondant au rangement croissant des numeros reels.
C     1ere adresse : APM1
C 
      CALL GSPOUE (18, APM1)
      ARM1=APM1-1
      APM2 =APM1+6
      ARM2 =APM2-1
      APM3 =APM2+6
      ARM3 =APM3-1
C  
C     Ce tableau depend du type de numerotation choisie
C 
      IF (NUM .EQ. 1) THEN
C  
        M(APM1)      =1
        M(APM1+1)    =2
        M(APM1+2)    =5
        M(APM1+3)    =6
        M(APM1+4)    =4
        M(APM1+5)    =3
C  
        M(APM2)      =1
        M(APM2+1)    =7
        M(APM2+2)    =31
        M(APM2+3)    =25
        M(APM2+4)    =13
        M(APM2+5)    =19
C  
        M(APM3)      =1
        M(APM3+1)    =3
        M(APM3+2)    =11
        M(APM3+3)    =9
        M(APM3+4)    =5
        M(APM3+5)    =7
C  
      ENDIF
C  
      IF (NUM .EQ. 2) THEN
C  
        M(APM1)      =1
        M(APM1+1)    =5
        M(APM1+2)    =4
        M(APM1+3)    =2
        M(APM1+4)    =6
        M(APM1+5)    =3
C  
        M(APM2)      =1
        M(APM2+1)    =19
        M(APM2+2)    =31
        M(APM2+3)    =13
        M(APM2+4)    =7
        M(APM2+5)    =25
C  
        M(APM3)      =1
        M(APM3+1)    =7
        M(APM3+2)    =11
        M(APM3+3)    =5
        M(APM3+4)    =3
        M(APM3+5)    =9
C  
      ENDIF
C 
C     Creation du tableau de rangement des numeros de ddl des couches
C     de la maniere suivante : les ddl correspondant aux deplacements
C     sont ranges d'abord, puis ceux correspondant aux derivees dans
C     l'ordre de calcul.
C 
C          - nom : RANCAL-COU
C          - 1ere adresse libre : M2
C          - numero : N2
C 
C     Modification : on remplit RANCAL-COU 
C     TAB(RANCAL(I), RANCAL(J))=K0N(NDDLC(I), NDDLC(J))
C 
      CALL GESTEN ('RANCAL-COU', 12, M2)
      R2    = M2-1
      K     = M2
C 
      DO I=1, 6
        X       = M(ARM1+I)
        M(K)    = X
        M(K+1)  = X+6
        K       = K+2
      ENDDO
C 
CD    CALL IMPTEN ('TABLEAU RANCAL-COU', M(M2), 1, 12)
C 
C     Creation du tableau de rangement des numeros de ddl des couches
C     de la maniere a ce que ddl(tab(i)) corresponde au rangement
C     calcul
C          - de nom RANINV-COU
C          - 1ere adresse libre: M2
C          - numero :N3
C 
      CALL GESTEN ('RANINV-COU', 36, M2)
      R2    = M2-1
      K     = M2
C 
      DECAL =0
      DO TYPDEP = 1 , 3
        DO I=1,6
          X     = M(ARM2+I)+DECAL
          M(K)  = X
          M(K+6)= X+1
          K     = K+1
        ENDDO
        K     = K+6
        DECAL = DECAL+2
      ENDDO
C 
CD     CALL IMPTEN('TABLEAU RANINV-COU',M(M2),1,36)
C  
C -----------------------------------------------------------------------
C       creation du tableau de rangement des numeros de ddl des couches
C       de la maniere suivante:
C       les ddl correspondant aux deplacements  sont ranges d'abord
C       puis ceux correspondants aux derivees dans l'ordre calcul :
C 
C            - DE NOM DEPINV-COU
C            - 1ere adresse libre: M2
C            - numero :N4
C -----------------------------------------------------------------------
      CALL GESTEN ('DEPINV-COU', 12, M2)
      R2    = M2-1
      K     = M2
C 
      DO I=1,6
        X       = M(ARM3+I)
        M(K)    = X
        M(K+6)  = X+1
        K       = K+1
      ENDDO
C 
CD    CALL IMPTEN ('TABLEAU DEPINV-COU', M(M2), 1, 12)
C 
C     remise des adresses des tableaux partiels a "zero"
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
C     On envoie comme arguments :
C 
C     E ...... TOUDEP   Tableau MAT-DEPLA developpes (DM(ADDEPA))
C     E ...... DINVCO   Tableau DEPINV-COU (M(ARINCO))
C     E ...... NDDLU    Numero des deplacements de U ranges croissant
C     E ...... NDDLV    Numero des deplacements de V ranges croissant
C     E ...... NDDLW    Numero des deplacements de W ranges croissant
C 
C     Et on recupere :
C 
C     S ...... DEPCAN   Tableau des valeurs des deplacements pour
C                       l'element ranges calcul par ordre
C                       croissant de numero de developpement
C                       c-a-D -NTDSFG ------- > NTDSFG
C 
      SUBROUTINE DPDVCA (TOUDEP, DINVCO, NDDLU, NDDLV, NDDLW, DEPCAN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  TOUDEP(NDDL*NBMAT) , DEPCAN(36*NBMAT)
      INTEGER           DINVCO(36) , NDDLU(12) , NDDLV(12) , NDDLW(12)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DPDVCA')
      INTEGER    NG, DECALN, K, U, V, W, NUDEV
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NG = 2*NTDSFG+1
      DECALN = 0
      K      = 1
      DO NUDEV = 1, NG
         DO U = 1, 12
           DEPCAN(K) = TOUDEP(DECALN+NDDLU(DINVCO(U)))
           K         = K+1
         ENDDO
         DO V = 1, 12
           DEPCAN(K) = TOUDEP(DECALN+NDDLV(DINVCO(V)))
           K         = K+1
         ENDDO
         DO W = 1, 12
           DEPCAN(K) = TOUDEP(DECALN+NDDLW(DINVCO(W)))
           K         = K+1
         ENDDO
         DECALN = DECALN+NDDL
      ENDDO
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine va chercher les noeuds ou on impose des conditions
C     aux limites.
C 
C     On envoie comme arguments :
C 
C     Le type de deplacement  TYPDEP
C 
C            si TYPDEP=1 ==========> U
C            si TYPDEP=2 ==========> V
C            si TYPDEP=1 ==========> W
C 
C     Le numero de bord
C     Le numero de couche ou de colonne suivant le type de bord
C     L'adresse de depart du tableau TLOCN1
C 
C     Et on recupere :
C 
C     DDL      tableau des numeros reels des ddl (correspondant aux
C              deplacements) aux bords et a nuco stockes croissants.
C     NBDDL    nombre de ddl concernes
C 
      SUBROUTINE NDDLCB (TYPDEP, NUBORD, NUCO, TLOCN1, DDL, NBDDL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      INTEGER       NUCO, NUBORD, NBDDL, TLOCN1, TYPDEP, DDL(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER       NDDLDE(12), DECAL
      CHARACTER*1   IDDEPL
      CHARACTER*6   IDPROG
      PARAMETER    (IDPROG='NDDLCB')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (TYPDEP .EQ. 1) THEN
        DECAL   =  1
        IDDEPL  = 'U'
      ELSE IF (TYPDEP .EQ. 2) THEN
        DECAL   =  3
        IDDEPL  = 'V'
      ELSE IF (TYPDEP .EQ. 3) THEN
        DECAL   =  5
        IDDEPL  = 'W'
      ELSE
        CALL ERREUD (0, 'MAUVAIS PASSAGE DE TYPE DE DEP DANS '//IDPROG)
      ENDIF
C 
      IF (NUM .EQ. 1) THEN
        IF (NUBORD .EQ. 1) THEN
          CALL NDDLC (TYPDEP, 1, NUCO, TLOCN1, NDDLDE)
          DDL(1)        = NDDLDE(1)
          DDL(2)        = NDDLDE(2)
          DDL(3)        = NDDLDE(3)
          DDL(4)        = NDDLDE(4)
          NBDDL         = 4
        ELSE IF (NUBORD.EQ.3) THEN
          CALL NDDLC (TYPDEP, NBCOU, NUCO, TLOCN1, NDDLDE)
          DDL(1)        = NDDLDE(9)
          DDL(2)        = NDDLDE(10)
          DDL(3)        = NDDLDE(11)
          DDL(4)        = NDDLDE(12)
          NBDDL         = 4
        ELSE IF (NUBORD.EQ .2) THEN
          CALL NDDLC (TYPDEP, NUCO, 1, TLOCN1, NDDLDE)
          DDL(1)        = NDDLDE(1)
          DDL(2)        = NDDLDE(2)
          DDL(3)        = NDDLDE(5)
          DDL(4)        = NDDLDE(6)
          DDL(5)        = NDDLDE(9)
          DDL(6)        = NDDLDE(10)
          NBDDL         = 6
        ELSE IF (NUBORD.EQ .4) THEN
          CALL NDDLC (TYPDEP, NUCO, NBCOL, TLOCN1, NDDLDE)
          DDL(1)        = NDDLDE(3)
          DDL(2)        = NDDLDE(4)
          DDL(3)        = NDDLDE(7)
          DDL(4)        = NDDLDE(8)
          DDL(5)        = NDDLDE(11)
          DDL(6)        = NDDLDE(12)
          NBDDL         = 6
        ELSE
          CALL ERREUD (0, ' MAUVAIS INDICATEUR DE BORD DANS :'//IDPROG)
        ENDIF
      ELSE IF (NUM .EQ. 2) THEN
        IF (NUBORD .EQ. 1) THEN
          CALL NDDLC (TYPDEP, 1, NUCO, TLOCN1, NDDLDE)
          DDL(1)        = NDDLDE(1)
          DDL(2)        = NDDLDE(2)
          DDL(3)        = NDDLDE(7)
          DDL(4)        = NDDLDE(8)
          NBDDL         = 4
        ELSE IF (NUBORD .EQ. 3) THEN
          CALL NDDLC (TYPDEP, NBCOU, NUCO, TLOCN1, NDDLDE)
          DDL(1)        = NDDLDE(5)
          DDL(2)        = NDDLDE(6)
          DDL(3)        = NDDLDE(11)
          DDL(4)        = NDDLDE(12)
          NBDDL         = 4
        ELSE IF (NUBORD .EQ. 2) THEN
          CALL NDDLC (TYPDEP, NUCO, 1, TLOCN1, NDDLDE)
          DDL(1)        = NDDLDE(1)
          DDL(2)        = NDDLDE(2)
          DDL(3)        = NDDLDE(3)
          DDL(4)        = NDDLDE(4)
          DDL(5)        = NDDLDE(5)
          DDL(6)        = NDDLDE(6)
          NBDDL         = 6
        ELSE IF (NUBORD .EQ. 4) THEN
          CALL NDDLC (TYPDEP, NUCO, NBCOL, TLOCN1, NDDLDE)
          DDL(1)        = NDDLDE(7)
          DDL(2)        = NDDLDE(8)
          DDL(3)        = NDDLDE(9)
          DDL(4)        = NDDLDE(10)
          DDL(5)        = NDDLDE(11)
          DDL(6)        = NDDLDE(12)
          NBDDL         = 6
        ELSE
          CALL ERREUD (0, 'VALEUR ILLICITE DE NUM DANS : '//IDPROG)
        ENDIF
      ENDIF
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
C     E ...... TABLE     Tableau double precision stocke
C                        (TABNIV(1), ...... , TABNIV(NBNIV))
C     E ...... NBNIV     Entier caracterisant le nombre d'indices
C                        du tableau d'entree
C     E ...... TABNIV    Tableau des valeurs des differents niveaux
C     E ...... NUNIV     Numero du niveau dont on veut extraire un vecteur
C     E ...... DEBUT     Adresse de debut dans table du vecteur
C                        que l'on veut extraire
C     Et on recupere :
C 
C     S ...... VECT      Vecteur extrait (de longueur TABNIV (NUNIV))
C     S ...... LONG      longueur du vecteur extrait
C 
C            Le vecteur obtenu correspond a tous les indices fixes sauf
C            l'indice nuniv que l'on fait varier de 1 a TABNIV(NUNIV).

      SUBROUTINE EXTRAD (TABLE, NBNIV, TABNIV, NUNIV, DEBUT, VECT, LONG)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER            NBNIV, TABNIV(NBNIV), NUNIV, DEBUT, LONG
      DOUBLE PRECISION   TABLE(*) , VECT(*)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     I , DECAL , TOT , K
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EXTRAD')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      DECAL = 1
      DO  I = 1, NUNIV-1
        DECAL = DECAL*TABNIV(I)
      ENDDO
      TOT     = DEBUT
      LONG    = TABNIV(NUNIV)
      K       = 1
      DO  I = 1 , LONG
        VECT(K) = TABLE(TOT)
        K       = K+1
        TOT     = TOT+DECAL
      ENDDO
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... TYPDEP le type de deplacement  TYPDEP
C                     si TYPDEP=1 ==========> U
C                     si TYPDEP=2 ==========> V
C                     si TYPDEP=3 ==========> W
C     E ...... NUCOU  le numero de couche
C     E ...... NUCOL  le numero de colonne
C     E ...... TLOCN1 l'adresse de depart du tableau TLOCN1
C 
C     Et on recupere :
C 
C     S ...... NDDLCA tableau des numeros reels des ddl correspondant au
C                     deplacements stockes calcul

      SUBROUTINE DDLCAL (TYPDEP, NUCOU, NUCOL, TLOCN1, NDDLCA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER NDDLCA(12)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER NDDLDE(12),I,NUCOU,NUCOL,TLOCN1,TYPDEP
      INTEGER RANINV
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DDLCAL ')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL NDDLC( TYPDEP, NUCOU, NUCOL, TLOCN1, NDDLDE)
      CALL ADTBM('DEPINV-COU',RANINV)
      RANINV = RANINV-1
      DO I = 1 , 12
        NDDLCA(I) = NDDLDE( M(RANINV+I) )
      END DO
CD    CALL RETOUD(IDPROG)
      RETURN
      END
