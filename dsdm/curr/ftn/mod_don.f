C     Cette routine lit et initialise les nombres de couches,
C     de colonnes, et les tableaux pour les donnees geometriques.
C 
      SUBROUTINE DONGEO
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
C     """""""""""""""""""""""""""""""""
C  
      INTEGER         LECINT, I
      INTEGER         NEPLU, NP, APD1, NLLU, ADHAU
      INTEGER         AD0, AC0, AD1, AC1, APC1
      INTEGER         J, ADRAY, ACRAY
      INTEGER         AM2LC, ADM2LC
      INTEGER         APM0, ARM0, NBORLU, ABORD
      INTEGER         NUMDIR
C 
      LOGICAL         LECLOG
C 
      DOUBLE PRECISION   LONGPRE, HAUT, LONBAN
      DOUBLE PRECISION   EPSILON, LOU, LECDBL
C 
      CHARACTER*20  FICDON
      CHARACTER*9   NOMDIR
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DONGEO')
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
C     Initialisation de la precision
C 
      EPSILON=1.D -12
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL LECSEQ ('DONNEES-PAR-FICHIER','DONNEES PAR FICHIER')
C 
      DONFIC    = LECLOG ('DONFIC')
      DONHER    = LECLOG ('DONHER')
      TRACT     = LECLOG ('TRACT')
      POLAR     = LECLOG ('POLAR')
C
      IF (DONHER) THEN
        CALL MESSAO ('DONNEZ LE NOMBRE D''ELEMENTS NBELZC DE LA ZONE DE 
     &               \CONTACT. SI NBELZC=0, ON RECALCULE LA REPARTITION 
     &               \DE PRESSION SUR CHAQUE ELEMENT A PARTIR DE LA FOR
     &               \CE DE CONTACT.')
         print*, DONFIC, DONHER, TRACT, POLAR, DONVOL
	 NBELZC = LECINT ('NBELZC')
      END IF
C 
      IF (DONFIC) THEN
C 
        CALL MESSAO  ('INDIQUEZ LES TYPES DE BORD CONCERNES PAR DES 
     $                \DONNEES PAR FICHIER')
C 
C       Gestion du tableau des entiers
C 
        CALL GSPOUE(4 , APM0)
        ARM0=APM0-1
        CALL MENAM(APM0,4)
        CALL LECLEN (M(APM0), 4, NBORLU)
        CALL GESTEN ('TYP-FICHIE', NBORLU, ABORD)
        CALL COPITE (NBORLU, M(APM0), M(ABORD))
C 
C       Il n'y pas besoin d'ouvrir le fichier qui existe
C       deja mais il faut declarer la directorie
C 
        NUMDIR     = 12
        NOMDIR     = 'donneefic'
        CHADIR( 12) = NOMDIR
C 
        DO I = 1 , NBORLU
          CALL LECSEQ ('NOM-FICHIER-DONNEE',
     &                 'POUR LE BORD '// bord(M(ARM0+I))//' DONNER '//
     &                 'LE NOM DU FICHIER DE DONNEES (6 CARACTERES) '//
     &                 'DANS LE REPERTOIRE donneefic')
          CALL LECCHA ('FICDON', FICDON)
          DO J =1, 20
            CHAFIC(12,I)(J:J) = ' '
          END DO
          CHAFIC(12, I) =  FICDON
        END DO
C 
      END IF
C 
      CALL LECSEQ ('GEOMETRIE', 'DONNEES POUR LA GEOMETRIE ET LE '//
     &             'MAILLAGE')
C 
      PREGRA = 1.D -6
C 
C     Donnees geometrie
C 
      RAYON  = LECDBL ('RAYON')
      RAYEXT =  LECDBL ('RAYEXT')
C 
      IF ( RAYEXT .LE. RAYON ) THEN
        CALL ERREUD (0, 'RAYEXT < RAYON DANS '// IDPROG)
      END IF
C 
C     Calcul de ces angles
C 
      CALL ANGGEO
C 
      CALL CATELE
C 
C     Creation du tableau des longueurs :
C            - nom : longueurs
C            - 1ere adresse libre: AD0
C            - numero :N0
C 
      CALL GESTDP ('LONGUEURS ', NBCOL, AD0)
      AC0=AD0-1
C 
10    CALL MESSAO ('ENTREE DES DIFFERENTES LONGUEURS DES ELEMENTS DU
     &              \MAILLAGE DANS L''ORDRE A PARTIR DU TROU. SI UNE
     &              \SEULE LONGUEUR EST LUE, CE SERA LA 1ere D''UNE
     &              \SERIE ARITHMETIQUE.')
C 
      CALL LECLDP (DM(AD0), NBCOL, NLLU)
      CALL IMPET ('NOMBRE DE LONGUEURS LUES ', NLLU)
C 
C    WARNING modif on divise les longueurs par deux pour avoir directement
C    l'element d'integration
C 
      DO I=1, NLLU
        DM(AC0+I)=DM(AC0+I)/2.D0
      ENDDO
C 
C    Creation du tableau des rayons :
C           - nom : RAYON
C           - 1ere adresse libre : ADRAY
C           - numero : NRAY
C 
      CALL GESTDP ('RAYON     ', NBCOL, ADRAY)
      ACRAY=ADRAY-1
      LONGPRE=RAYON
C 
C     Pour calculer la longueur de la bande (remplissage
C     des tableaux LONGUEURS et RAYON)
C 
      LONBAN = 0.D0
C 
      IF (NLLU .EQ. 1) THEN
C 
C     LOU est l'unique longueur lue /2.D0*NBCOL
C     LOU=DM(AD0)/DBLE(NBCOL)
C 
        LOU = 0.D0
C 
        IF (NBCOL .GT. 1) THEN
C 
C         LOU est la raison multiplee par deux
C 
          LOU =  (RAYEXT-RAYON) /DBLE (NBCOL) -2.D0*DM(AD0)
          LOU = LOU  / DBLE (NBCOL-1)
        END IF
C 
        DO I=1, NBCOL-1
          DM(AC0+I)       = DM(AD0)+DBLE(I-1)*LOU
          LONBAN          = LONBAN+2.D0*DM(AC0+I)
          DM(ACRAY+I)     = LONGPRE+DM(AC0+I)
          LONGPRE         = DM(ACRAY+I)+DM(AC0+I)
        ENDDO
C 
       DM(AC0+NBCOL)      = ( RAYEXT-RAYON - LONBAN )/2.D0
       LONBAN             =  LONBAN+2.D0*DM(AC0+NBCOL)
       DM(ACRAY+NBCOL)    =  LONGPRE+DM(AC0+NBCOL)
C 
       CALL MESSAO ('DANS LE CAS DE LA RAISON ARITHMETIQUE ')
C 
       CALL IMPTDT ('TABLEAU DES DEMI-LONGUEURS ', DM(AD0), 1, NBCOL)
C 
       CALL IMPTDT ('TABLEAU DES RAYONS AU CENTRE ',
     &               DM(ADRAY), 1, NBCOL)
C 
      ELSE
C 
        DO I=1, NBCOL
          LONBAN          = LONBAN+2.D0*DM(AC0+I)
          DM(ACRAY+I)     = DM(AC0+I)+LONGPRE
          LONGPRE         = DM(ACRAY+I)+DM(AC0+I)
        ENDDO
      ENDIF
C 
C     Test pour verifier que toutes les longueurs sont positives
C 
      DO I = 1, NBCOL
        IF (DM(AC0+I) .LT. 0.D0) THEN
          CALL IMPET ('LONGUEUR DE L''ELEMENT NEGATIVE OU NULLE ', I)
          CALL IMPDT ('LONGUEUR ', DM(AC0+I))
          CALL ERREUD (0, 'ERREUR DE DONNEE DANS '// IDPROG)
        END IF
      END DO
C 
C     Test pour savoir si le nombre total de colonnes lu NLLU=NBCOL
C 
      IF ((NLLU .NE. NBCOL) .AND. (NLLU.NE.1)) THEN
        CALL MESSAO ('VOUS N''AVEZ PAS RENTRE LE NOMBRE 
     &               \TOTAL DE COUCHES INDIQUE : RECOMMENCEZ ')
        GOTO 10
      ELSE IF (RAYON+LONBAN .NE. RAYEXT ) THEN
        CALL IMPDT ('RAYEXT ', RAYEXT)
        CALL IMPDT ('LONGUEUR DE BANDE ', LONBAN)
        CALL IMPDT ('RAYON ', RAYON)
        CALL MESSAO ('RAYON+LONBAN .NE. RAYEXT RECOMMENCEZ LA DONNEE')
        GOTO 10
      ENDIF
C 
C     creation du tableau des epaisseurs :
C            - nom: epaisseurs
C            - 1ere adresse libre AD1
C 
      CALL GESTDP ('EPAISSEURS', NBCOU, AD1)
C 
      AC1=AD1-1
C 
C     creation d'un tableau provisoire pour ranger les epaisseurs :
C            - 1ere adresse libre APD1
C 
      CALL GSPOUD (NBCOU, APD1)
      APC1=APD1-1
      CALL MENADM (APD1, NBCOU)
C 
C     Rangement dans DM a partir de APD1 des epaisseurs lues
C     Les epaisseurs sont divisees par 2 pour avoir directement
C     l'element d'integration
C 
3     CALL MESSAO ('ENTREE DES DIFFERENTES EPAISSEURS : DE LA COUCHE
     &             \INFERIEURE A LA COUCHE SUPERIEURE. SI UNE SEULE
     &             \COUCHE EST LUE, TOUTES LES COUCHES ONT LA MEME
     &             \EPAISSEUR. SI LE STRATIFIE EST SYMETRIQUE, ON RENTRE
     &             \LES EPAISSEURS DES COUCHES DE LA COUCHE MOYENNE A
     &             \LA COUCHE SUPERIEURE. SI LE NOMBRE DE COUCHES EST
     &             \IMPAIR, ON DONNE L''EPAISSEUR TOTALE DE LA COUCHE
     &             \CENTRALE.')
C 
C     Lit une liste de variables double precision, et remplit
C     le tableau provisoire cree precedemment.

      CALL LECLDP (DM(APD1), NBCOU, NEPLU)
C 
      IF (NEPLU .NE. 1 .AND. NEPLU .NE. NBCOU) THEN
        CALL IMPET ('NOMBRE D''EPAISSEURS DONNEES DIFFERENT DE NBCOU 
     &               ET DE 1 ', NEPLU)
        CALL IMPET ('NOMBRE DE COUCHES POUR LE CALCUL NBCOU ',NBCOU)
        CALL MESSAD ('VOUS N''AVEZ PAS LE NOMBRE TOTAL DE COUCHES 
     &               INDIQUE : RECOMMENCEZ')        
      GOTO 3
      END IF
C 
C     S'il n'y a qu'une epaisseur (NEPLU=1) : rangement dans epaisseurs
C 
      IF (NEPLU .EQ. 1) THEN
        DO I=1, NBCOU
          DM(AC1+I)=DM(APD1)/2.D0
        ENDDO
        NP=NBCOU
      ELSE
        DO I=1, NEPLU
          DM(AC1+I)=DM(APD1+I-1)/2.D0
        ENDDO
      END IF
C 
      IF (SYM .AND. (.NOT. SYMPAR)) THEN
        DM(AD1) = DM(AD1)/2.D0
      ENDIF
C 
C     Test pour savoir si une epaisseur a ete affectee par couche :
C     on regarde s'il n'y a pas d'epaisseur nulle (<EPSILON)
C 
      DO I=0, (NBCOU-1)
        IF (DABS(DM(AD1+I)) .LT. EPSILON) THEN
           CALL MESSAO ('UNE DES COUCHES A UNE EPAISSEUR NULLE 
     $                  \RECOMMENCEZ LA SAISIE DES EPAISSEURS')
           GOTO 3
        ENDIF
      ENDDO
C 
C     Creation du tableau des hauteurs au centre :
C            - nom : HAUTEURS
C            - 1ere adresse libre : ADHAU
C            - numero : NHAU
C 
      CALL GESTDP ('HAUTEURS  ', NBCOU, ADHAU)
      DM(ADHAU) = DM(AD1)
      HAUT = 2.D0*DM(ADHAU)
C 
C     Calcul de la mi-epaisseur
C 
      THICK = DM(ADHAU)
C 
      DO I=1, NBCOU-1
         DM(ADHAU+I)  = HAUT+ DM(AD1+I)
         HAUT         = DM(ADHAU+I)+DM(AD1+I)
         THICK        = THICK + DM(AD1+I)
      ENDDO
C 
      IF (SYM) THICK = 2.D0*THICK
C 
CD    CALL IMPDT ('DEMI-EPAISSEUR ', THICK)
C 
C     Cherche le numero du tableau double precision LONGUEURS
C 
      CALL ADTBDM ('LONGUEURS                  ', AD0)
      CALL IMPTDT ('TAB DES DEMI-LONGUEURS     ', DM(AD0), 1, NBCOL)
      CALL IMPTDT ('TAB DES DEMI-EPAISSEURS    ', DM(AD1), 1, NBCOU)
      CALL IMPTDT ('TAB DES HAUTEURS-AU-CENTRE ', DM(ADHAU), 1, NBCOU)
      CALL IMPTDT ('TAB DES RAYONS-AU-CENTRE   ', DM(ADRAY), 1, NBCOL)
C 
C     Pour le calcul du raccord plaque
C     Le tableau HAU-CE-TOT contient les distances au plan median
C     de tous les plans du stratifie a partir du plan inferieur.
C 
      CALL GESTDP ('HAU-CE-TOT', NBCOU+1, ADHAU)
C 
      HAUT            = -THICK
      IF (SYM) HAUT   = 0.D0
      DM(ADHAU)       =  HAUT
C 
      DO I=1,NBCOU
         DM(ADHAU+I)  = HAUT+ 2.D0*DM(AD1+I-1)
         HAUT         = DM(ADHAU+I)
      ENDDO
C 
CD    CALL IMPTDT ('TABLEAU DES HAUTEURS DES PLANS ',
CD   &              DM(ADHAU),1,NBCOU+1)
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
C     Cette routine remplit le tableau des angles des bandes en foncion
C     de la valeur de ntdsfg
 
      SUBROUTINE ANGGEO
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER  I, ADTETA
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ANGGEO')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      CALL GESTDP ('ANGLES-GEO', NTETA, ADTETA)
      IF (NTETA .GT. 0)THEN
        DO I = 0, NTETA-1
          DM(ADTETA+I) = 2.D0*PI*DBLE(I)/DBLE(NTETA)
        ENDDO
      END IF
CD    CALL IMPTDN ('ANGLES GEOMETRIQUES', DM(ADTETA), 1, NTETA)
CD    CALL RETOUD (IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... N0      numero de couche
C     E ...... YO      valeur Y initiale
C     E ...... YC      valeur Y de critique
C     E ...... B       valeur du couplage d'endommagement
C     E ...... E220    valeur initiale du module E22
C     E ...... G120    valeur initiale du module 2*G12
C     E ...... R0      valeur du seuil de plasticite
C     E ...... BETA    valeurs telle que p =  beta * R
C     E ...... ALPHA
C     E ...... A2      valeur du couplage de plasticite

C     On recupere dans le tableau COMNLI stocke sequentiellement comme suit :
C 
C     S ...... S11INF  valeur mini de la contrainte 11
C     S ...... S11SUP  valeur maxi de la contrainte 11
C     S ...... S22SUP  valeur maxi de la contrainte 22
C     S ...... M0      valeur de l'adresse de depart du tableau des souplesses
C                      elastiques dans la base d'orthotropie
C     S ...... ADHONP  valeur de l'adresse de depart du tableau des rigidites
C                      elastiques dans la base d'orthotropie
C 
      SUBROUTINE CNLINC (N0, COMNLI, M0, ADHONP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER N0, M0, ADHONP
C 
      DOUBLE PRECISION COMNLI(12)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER M2, P, R, MNLIN, ADCAEN, R1
C 
CD     LOGICAL LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CNLINC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('SOUP-ORTHO', R)
      CALL ADTBDM ('HOOK-ORTHO', R1)
      CALL ADTBM ('TYP-COUCHE', M2)
      CALL ADTBDM ('CARAC-NONL', ADCAEN)
C 
C     P caracterise le type de comportement
C 
      P=M(M2+N0-1)
C 
CD    CALL IMPEP('VALEUR DU TYPE DE COMPORTEMENT',P)
C 
      M0        = R+17*(P-1)
      ADHONP    = R1+17*(P-1)
C 
      MNLIN = ADCAEN+10*(P-1)
C 
      COMNLI( 1 )  = DM( MNLIN )
      COMNLI( 2 )  = DM( MNLIN + 1 )
      COMNLI( 3 )  = DM( MNLIN + 2 )
      COMNLI( 4 )  = 1.D0/DM( M0+4 )
      COMNLI( 5 )  = 1.D0/DM( M0+8 )
      COMNLI( 6 )  = DM( MNLIN + 3 )
      COMNLI( 7 )  = DM( MNLIN + 4 )
      COMNLI( 8 )  = DM( MNLIN + 5 )
      COMNLI( 9 )  = DM( MNLIN + 6 )
      COMNLI( 10 ) = DM( MNLIN + 7 )
      COMNLI( 11 ) = DM( MNLIN + 8 )
      COMNLI( 12 ) = DM( MNLIN + 9 )
C 
CD    IF( LTRACN(1) ) THEN
CD      CALL IMPEN('VALEUR DUNUMERO DE COUCHE N0',N0)
CD      CALL IMPEN('VALEUR DE LA IERE ADRESSE M0',M0)
CD      CALL OMPTDN(' Valeur des souplesses orthotropes'
CD                  , DM(M0),17,1)
CD      CALL IMPEN('VALEUR DE LA IERE ADRESSE MNLIN',MNLIN)
CD      CALL OMPTDN(' Valeur des caracteristiques non lineaires '
CD                  , COMNLI(1) , 12 , 1 )
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
C     E ...... N0      numero d'interface .
C 
C     On recupere dans le tableau COMNLI stocke sequentiellement comme suit :
C 
C     S ...... YC      valeur de Y critique
C     S ...... GAM1    valeur du couplage d'endommagement 1
C     S ...... GAM2    valeur du couplage d'endommagement 2
C     S ...... E10     valeur initiale du module E1
C     S ...... E20     valeur initiale du module E2
C     S ...... E30     valeur initiale du module E3
C     S ...... SOUORT  valeur du tableau des souplesses dans la base d'orthotropie
C 
      SUBROUTINE CNLINI (N0, COMNLI, SOUORT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER N0
      DOUBLE PRECISION COMNLI(6) , SOUORT(9)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER M2 , P , R , MNLIN , ADSOUP , ADCAEN , M0
C 
CD    LOGICAL LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CNLINI')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('SOIN-ORTHO' , ADSOUP )
      CALL ADTBDM('HOIN-ORTHO' , R )
      CALL ADTBM( 'TYP-COUCHE' , M2 )
      CALL ADTBDM('CARAI-NONL' , ADCAEN )
C 
C     P caracterise le type de comportement
C 
      P=M(M2+NBCOU+N0-1)
C 
CD    CALL IMPEP('VALEUR DU TYPE DE COMPORTEMENT',P)
C 
      M0 = r+3*(P-1)
C 
      MNLIN = ADCAEN+3*(P-1)
C 
      COMNLI( 1 )  = DM( MNLIN )
      COMNLI( 2 )  = DM( MNLIN + 1 )
      COMNLI( 3 )  = DM( MNLIN + 2 )
      COMNLI( 4 )  = DM( M0)
      COMNLI( 5 )  = DM( M0+1 )
      COMNLI( 6 )  = DM( M0+2 )
C 
      M0 = ADSOUP+3*(P-1)
C 
      SOUORT( 1 )  = DM(M0)
      SOUORT( 2 )  = 0.D0
      SOUORT( 3 )  = 0.D0
      SOUORT( 4 )  = 0.D0
      SOUORT( 5 )  = DM(M0+1)
      SOUORT( 6 )  = 0.D0
      SOUORT( 7 )  = 0.D0
      SOUORT( 8 )  = 0.D0
      SOUORT( 9 )  = DM(M0+2)
C 
CD    IF( LTRACN(1) ) THEN
CD      CALL IMPEN('VALEUR DU NUMERO D''INTERFACE N0',N0)
CD      CALL IMPEN('VALEUR DE LA IERE ADRESSE M0',M0)
CD      CALL IMPEN('VALEUR DE LA IERE ADRESSE ADSOUP', ADSOUP )
CD      CALL OMPTDN('Valeur des souplesses dans la base d''orthotropie'
CD                , DM(ADSOUP) , 3 , 1 )
CD      CALL IMPEN('VALEUR DE LA IERE ADRESSE MNLIN',MNLIN)
CD      CALL OMPTDN(' Valeur des caracteristiques non lineaires '
CD                , COMNLI(1) , 6 , 1 )
CD      CALL OMPTDN(' Valeur des souplesses ( orthotropie )  '
CD                , SOUORT(1) , 9 , 1 )
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
