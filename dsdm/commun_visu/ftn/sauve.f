C     Pour sauvegarder tout les tableaux definitifs des commons en vue
C     d'etre reaffectes apres point-reprise
C 
      SUBROUTINE SAUCVR (TYPE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'typcal.h'
      include 'strategie_calcul.h'
C 
      INTEGER TYPE
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I, IERNAM, NBETDE
      INTEGER  IUNIT, LONENR, RESTE, NBENR, DEBUT
C 
      CHARACTER*6   IDPROG
      CHARACTER*13  NOMFIC
      CHARACTER*3   BARATI

      PARAMETER (IDPROG='SAUCVR')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL IDENTI (NBETGL, BARATI)
      NOMFIC = 'sau_etglo_'//BARATI
C 
      CALL IMPCT ('NOM DU FICHIER DE SAUVEGARDE ',
     &             tild1(1:LTILD1)//'/'//NOMFIC)
C 
      IUNIT = 7
C 
CD    IF ((NBETGL .GT. 0) .AND. ((NBETGL/3))*3 .NE. NBETGL) THEN
      CALL IDENTI (NBETGL-1, BARATI)
      NOMFIC = 'sau_etglo_'//BARATI
      OPEN (UNIT=IUNIT, FILE=tild1(1:LTILD1)//'/'//NOMFIC,
     &      FORM='UNFORMATTED', ACCESS='SEQUENTIAL',
     &      IOSTAT=IERNAM, ERR=900)
      CLOSE (IUNIT, ERR=1100, STATUS='DELETE')
CD    END IF
C 
C     10 pages de doubles
C 
      LONENR = 2560
C 
      IF (TYPE .EQ . 0) THEN
        CALL IDENTI (NBETGL, BARATI)
        NOMFIC = 'sau_etglo_'//BARATI
      ELSE IF (TYPE .EQ . 1) THEN
        CALL IDENTI (NBETLC, BARATI)
        NOMFIC = 'sau_etloc_'//BARATI
      ELSE
        CALL IMPET ('Type de sauvegarde inconnu ', TYPE)
        CALL ERREUD (0, 'Mauvais passage de type dans '//IDPROG )
      END IF
C 
      OPEN (UNIT=IUNIT, FILE=tild1(1:LTILD1)//'/'//NOMFIC,
     &      FORM='UNFORMATTED', ACCESS='SEQUENTIAL',
     &      IOSTAT=IERNAM, ERR=900)
C 
      WRITE (IUNIT) ADM
      WRITE (IUNIT) AM
      WRITE (IUNIT) LONGDM
      WRITE (IUNIT) AM2EN
      WRITE (IUNIT) LONGM
      WRITE (IUNIT) ADM1
      WRITE (IUNIT) ADM2
      WRITE (IUNIT) AM1
      CALL IMPET ('LONGUEUR DES REELS ', ADM2-ADM1)
      WRITE (IUNIT) AM2
      WRITE (IUNIT) ADM2EN
      WRITE (IUNIT) NBTADM
      WRITE (IUNIT) NBTAM
      WRITE (IUNIT) VERAM
      WRITE (IUNIT) VERADM
C 
      CALL ECFICE (IUNIT, M(1), 1, AM1-1)
C 
      DEBUT = 1
      NBENR = (ADM1-1)/LONENR
      RESTE = (ADM1-1) - LONENR*NBENR
C 
      DO I = 1, NBENR
        CALL ECFICD (IUNIT, DM(DEBUT), DEBUT, LONENR)
        DEBUT = DEBUT+LONENR
      END DO
C 
      IF (RESTE .GT. 0) THEN
        CALL ECFICD (IUNIT, DM(DEBUT), DEBUT, RESTE)
      ENDIF
C 
      WRITE (IUNIT) NBCOU
      WRITE (IUNIT) NBCOL
      WRITE (IUNIT) NNOEUD
      WRITE (IUNIT) NDDL
      WRITE (IUNIT) NEL1
      WRITE (IUNIT) NEL2
      WRITE (IUNIT) NUM
      WRITE (IUNIT) NBZONE
      WRITE (IUNIT) NCDPIM
      WRITE (IUNIT) NCBPIM
      WRITE (IUNIT) NCEFIM
      WRITE (IUNIT) NCBLIM
      WRITE (IUNIT) NBFODO
C 
      WRITE (IUNIT) RAYON
      WRITE (IUNIT) RAYEXT
      WRITE (IUNIT) THICK
C 
      WRITE (IUNIT) NKCOU
      WRITE (IUNIT) NKINT
      WRITE (IUNIT) NBINT
      WRITE (IUNIT) NBANGL
C 
      WRITE (IUNIT) XINTEG
      WRITE (IUNIT) YINTEG
      WRITE (IUNIT) NTDSF
      WRITE (IUNIT) NTDSFG
      WRITE (IUNIT) NBMAT
      WRITE (IUNIT) NTMAT
C 
      WRITE (IUNIT) GAUSS
      WRITE (IUNIT) POIDS
      WRITE (IUNIT) PI
C 
      WRITE (IUNIT) CHARM
      WRITE (IUNIT) CHARDM
C 
      WRITE (IUNIT) NSIG
      WRITE (IUNIT) NSAU
      WRITE (IUNIT) NEPS
      WRITE (IUNIT) NGAU1
      WRITE (IUNIT) NDDLEL
      WRITE (IUNIT) NTETA
      WRITE (IUNIT) NGAU2
C 
      WRITE (IUNIT) NBFICH
      WRITE (IUNIT) NUMFIC
      WRITE (IUNIT) NBDIR
      WRITE (IUNIT) LONNOM
C 
      WRITE (IUNIT) CHADIR
      WRITE (IUNIT) CHAFIC
C 
      WRITE (IUNIT) NOTAXE
      WRITE (IUNIT) NOMSAU
C 
      WRITE (IUNIT) ELAS
C 
      WRITE (IUNIT) NPICET
      WRITE (IUNIT) NFOTPS
      WRITE (IUNIT) NBNETT
CD    WRITE (IUNIT) NBNEMX
      WRITE (IUNIT) NBETLC
      WRITE (IUNIT) NPICMX
      WRITE (IUNIT) NBETLT
      WRITE (IUNIT) NPICAC
      WRITE (IUNIT) PICREP
      WRITE (IUNIT) NBETGL
C 
      WRITE (IUNIT) DUREE
      WRITE (IUNIT) ERRLOC
      WRITE (IUNIT) DUREAC
      WRITE (IUNIT) DCONV1
      WRITE (IUNIT) LTEMPS
C 
      WRITE (IUNIT) NFTGLO
C 
      WRITE (IUNIT) CHAMAX
      WRITE (IUNIT) NBFSIG
      WRITE (IUNIT) NBFEPS
      WRITE (IUNIT) NBDEPR
      WRITE (IUNIT) NBDEVR
      WRITE (IUNIT) NBDPTR
      WRITE (IUNIT) CHARAX
      WRITE (IUNIT) DEADTR
      WRITE (IUNIT) SAADTR
      WRITE (IUNIT) EVCOTR
      WRITE (IUNIT) COTORE
      WRITE (IUNIT) CNTORE
C 
      WRITE (IUNIT) PREGRA
C 
      WRITE (IUNIT) LENDCO
      WRITE (IUNIT) LRUPCO
      WRITE (IUNIT) LPLACO
      WRITE (IUNIT) LENINT
      WRITE (IUNIT) LPLAIN
C 
C     encours
C 
CD    WRITE (IUNIT) LMOSIG
CD    WRITE (IUNIT) LRETAR
CD    WRITE (IUNIT) LRETIN
CD    WRITE (IUNIT) LPLAEF
CD    WRITE (IUNIT) LPLEFF
C 
C     fincour
C 
      WRITE (IUNIT) SYM
      WRITE (IUNIT) SYMPAR
C 
      WRITE (IUNIT) SYMX
      WRITE (IUNIT) SYMY
      WRITE (IUNIT) SYMO
      WRITE (IUNIT) FORS
C 
CD    CALL IMPLT ('ENDOMMAGEMENT COUCHE    ', LENDCO)
CD    CALL IMPLT ('RUPTURE FIBRE           ', LRUPCO)
CD    CALL IMPLT ('PLASTICITE COUCHE       ', LPLACO)
CD    CALL IMPLT ('ENDOMMAGEMENT INTERFACE ', LENINT)
CD    CALL IMPLT ('PLASTICITE INTERFACE    ', LPLAIN)
C 
C     STEPHANE LE 22/10/98
C     ON SAUVE NUPATE POUR LE RECUPERER LORS D'UNE REPRISE
C 
      WRITE (IUNIT) NUPATE
C 
      CLOSE (IUNIT)
C  
C     On va essayer de reconstruire les fichiers necessaires a la
C     realisation des etapes locales et globales => pas reecfi
C  
C     => la reprise s'effectue apres une etape globale.
C  
C     CALL REECFI
C  
CD    CALL RETOUD(IDPROG)
C 
      RETURN
C 
1100  CALL IMPET (' NUMERO D''ERREUR ', IERNAM)
      CALL IMPET (' A l''etape globale ', NBETDE)
      CALL IMPCT (' FICHIER', tild1(1:LTILD1)//'/'//NOMFIC)
      CALL ERREUD 
     &(0, 'PB POUR DELETER LE FICHIER APRES OPEN DANS '//IDPROG)
C 
900   CALL IMPCT (' FICHIER', tild1(1:LTILD1)//'/'//NOMFIC)
      CALL IMPET (' A l''etape globale ', NBETDE )
      CALL ERREUD (0,'ERREUR D''OUVERTURE DANS '//IDPROG)
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Pour reprendre tous les tableaux definitifs des commons a une etape
C     globale donnee, avant (type = 0) ou apres (type = 1) etape locale en
C     vue de les reaffecter apres point-reprise.
C 
C     On envoie comme arguments :
C    
C     E ...... TYPE 0   etape globale
C              TYPE 1   etape locale
C     E ...... NBETDE   numero d'etape globale de reprise
C     E ...... COMPAC   logique indiquant si l'on doit compacter les champs
C                       de l'espace au nombre necessaire pour la visu

      SUBROUTINE REPCVR (TYPE, NBETDE, COMPAC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C    
      include 'cominc.h'
      include 'typcal.h'
      include 'strategie_calcul.h'
C 
      LOGICAL COMPAC
C 
      INTEGER NBETDE, TYPE
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER I, IERNAM
      INTEGER IUNIT, LONENR, RESTE, NBENR, DEBUT
      INTEGER NOUCHR, NOUCHA
      INTEGER ADANG
      INTEGER ETLOLO, ETGLLO
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='REPCVR')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C    10 pages de doubles
C 
      LONENR = 2560
C 
      IF (TYPE .EQ. 0) THEN
        ETGLLO = NBETDE
        ETLOLO = NBETDE
      ELSE IF (TYPE .EQ. 1) THEN
        ETGLLO = NBETDE
        ETLOLO = NBETDE+1
      END IF
C 
      IUNIT = 7
C 
      OPEN (UNIT=IUNIT, FILE=tild2(1:LTILD2),
     &      FORM='UNFORMATTED', STATUS='OLD',
     &      ACCESS='SEQUENTIAL', IOSTAT=IERNAM , ERR=900)
C 
      print*, 'FILE ', tild2(1:LTILD2)
      READ(IUNIT) ADM
      READ(IUNIT) AM
      READ(IUNIT) LONGDM
      READ(IUNIT) AM2EN
      READ(IUNIT) LONGM
      READ(IUNIT) ADM1
      READ(IUNIT) ADM2
      READ(IUNIT) AM1
      READ(IUNIT) AM2
      READ(IUNIT) ADM2EN
      READ(IUNIT) NBTADM
      READ(IUNIT) NBTAM
      READ(IUNIT) VERAM
      READ(IUNIT) VERADM
C 
      CALL IMPET ('LONGUEUR DIPONIBLE LDMEFF DANS DM', LDMEFF)
      CALL IMPET ('LONGUEUR A LIRE DANS DM', ADM1)
C 
      CALL LIFCEU (M(1), 1, AM1-1, IUNIT)
C 
      DEBUT = 1
      NBENR = (ADM1-1)/LONENR
      RESTE = ( ADM1-1 )- LONENR*NBENR
C 
      DO I = 1, NBENR
        CALL LIFCDU (DM(DEBUT), DEBUT, LONENR, IUNIT)
        DEBUT = DEBUT+LONENR
      END DO
C 
      IF (RESTE .GT. 0) THEN
        CALL LIFCDU (DM(DEBUT), DEBUT, RESTE, IUNIT)
      ENDIF
C 
      READ(IUNIT) NBCOU
      READ(IUNIT) NBCOL
      READ(IUNIT) NNOEUD
      READ(IUNIT) NDDL
      READ(IUNIT) NEL1
      READ(IUNIT) NEL2
      READ(IUNIT) NUM
      READ(IUNIT) NBZONE
      READ(IUNIT) NCDPIM
      READ(IUNIT) NCBPIM
      READ(IUNIT) NCEFIM
      READ(IUNIT) NCBLIM
      READ(IUNIT) NBFODO
C 
      READ(IUNIT) RAYON
      READ(IUNIT) RAYEXT
      READ(IUNIT) THICK
C 
      READ(IUNIT) NKCOU
      READ(IUNIT) NKINT
      READ(IUNIT) NBINT
      READ(IUNIT) NBANGL
C 
      READ(IUNIT) XINTEG
      READ(IUNIT) YINTEG
      READ(IUNIT) NTDSF
      READ(IUNIT) NTDSFG
      READ(IUNIT) NBMAT
      READ(IUNIT) NTMAT
C 
      READ(IUNIT) GAUSS
      READ(IUNIT) POIDS
      READ(IUNIT) PI
C 
      READ(IUNIT) CHARM
      READ(IUNIT) CHARDM
C 
      READ(IUNIT) NSIG
      READ(IUNIT) NSAU
      READ(IUNIT) NEPS
      READ(IUNIT) NGAU1
      READ(IUNIT) NDDLEL
      READ(IUNIT) NTETA
      READ(IUNIT) NGAU2
C 
      READ(IUNIT) NBFICH
C 
C    MODIF pour ne pas limiter betement
C 
      CALL BALAIE (NBDIRM, NBFICH)
C 
      READ (IUNIT) NUMFIC
C 
C    MODIF pour ne pas limiter betement
C 
      CALL BALAIE (NBFIMX*NBDIRM, NUMFIC)
C 
      READ(IUNIT) NBDIR
      READ(IUNIT) LONNOM
C 
      READ(IUNIT) CHADIR
      READ(IUNIT) CHAFIC
C 
      READ(IUNIT) NOTAXE
      READ(IUNIT) NOMSAU
C 
      READ(IUNIT) ELAS
C 
      READ(IUNIT) NPICET
      READ(IUNIT) NFOTPS
      READ(IUNIT) NBNETT
CD    READ(IUNIT) NBNEMX
      READ(IUNIT) NBETLC
C 
      IF (ETLOLO .NE. NBETLC) THEN
        CALL IMPET ('REPRISE DEMANDEE A L''ETAPE '//
     &              'GLOBALE NUMERO ', NBETDE)
        CALL IMPET ('POUR LE TYPE ', TYPE)
        CALL IMPET ('NUMERO D''ETAPE LOCALE THEORIQUE ', ETLOLO)
        CALL IMPET ('NUMERO D''ETAPE LOCALE LU EN REPRISE ', NBETLC)
        CALL MESSAO ('WHAT''S THAT MESS !!!! '// IDPROG)
      END IF
C 
      READ(IUNIT) NPICMX
      READ(IUNIT) NBETLT
      READ(IUNIT) NPICAC
      READ(IUNIT) PICREP
      READ(IUNIT) NBETGL
C 
      IF (NBETDE .GT. NBETGL) THEN
        CALL IMPET ('REPRISE DEMANDEE A L''ETAPE GLOBALE NUMERO     ',
     &               NBETDE)
        CALL IMPET ('SAUVEGARDE EFFECTUEE A L''ETAPE GLOBALE NUMERO ',
     &               NBETGL)
        CALL MESSAO ('WHAT''S THAT MESS !!!! '// IDPROG)
      END IF

      READ(IUNIT) DUREE
      READ(IUNIT) ERRLOC
      READ(IUNIT) DUREAC
      READ(IUNIT) DCONV1
      READ(IUNIT) LTEMPS
C 
      READ(IUNIT) NFTGLO
C 
      READ(IUNIT) CHAMAX
      READ(IUNIT) NBFSIG
      READ(IUNIT) NBFEPS
      READ(IUNIT) NBDEPR
      READ(IUNIT) NBDEVR
      READ(IUNIT) NBDPTR
      READ(IUNIT) CHARAX
      READ(IUNIT) DEADTR
      READ(IUNIT) SAADTR
      READ(IUNIT) EVCOTR
      READ(IUNIT) COTORE
      READ(IUNIT) CNTORE
C 
      READ(IUNIT) PREGRA
      READ(IUNIT) LENDCO
      READ(IUNIT) LRUPCO
      READ(IUNIT) LPLACO
      READ(IUNIT) LENINT
      READ(IUNIT) LPLAIN
C 
C    encours
C 
C    CD     READ(IUNIT) LMOSIG
C    CD     READ(IUNIT) LRETAR
C    CD     READ(IUNIT) LRETIN
C    CD     READ(IUNIT) LPLAEF
C    CD     READ(IUNIT) LPLEFF
C 
C    fincours
C 
      READ(IUNIT) SYM
      READ(IUNIT) SYMPAR
C 
      READ(IUNIT) SYMX
      READ(IUNIT) SYMY
      READ(IUNIT) SYMO
      READ(IUNIT) FORS
C 
C     STEPHANE LE 22/10/98
C     ON RECUPERE LE NUPATE
C 
      READ(IUNIT) NUPATE
C 
      CLOSE(IUNIT)
C 
C    derniere adresse libre dans M
C 
      AM2    =LM
      AM2EN  =LM
C 
C    derniere adresse libre dans DM
C 
      ADM2EN    =LDMEFF
      ADM2      =LDMEFF
C 
C    COMPAC en mode visualisation, pour appeler seulement les champs dont on a besoin ?
C 
      IF (COMPAC) THEN
        CALL IMPET ('NOMBRE TOTAL DE CHAMPS RESERVES     ', CHARAX)
        CALL IMPET ('NOMBRE D''EVOLUTIONS EN DEFORMATION ', NBDPTR)
        CALL IMPET ('NOMBRE D''EVOLUTIONS EN CONTRAINTES ', EVCOTR)
        CALL IMPET ('NOMBRE DE DEFORMATIONS CREEES       ', DEADTR)
        CALL IMPET ('NOMBRE DE CONTRAINTES CREEES        ', COTORE)
        IF( (DEADTR .NE. CHARAX) .AND. (COTORE .NE. CHARAX)) THEN
          NOUCHR = MAX0 (DEADTR, COTORE)
          NOUCHA = MAX0 (NBFSIG, NBFEPS)
          CALL MNBCHA (NOUCHR, NOUCHA)
          CALL SAUCVR (TYPE)
          COMPAC = .FALSE.
        END IF
      END IF
C 
C     On va essayer de reconstruire les fichiers necessaires a la
C     realisation des etapes locales et globales => pas reecfi
C  
C     => la reprise s'effectue apres une etape globale.
C  
C     CALL RELIFI(NUREP)
C 
C     INFORMATION GENERALE NOTAMMENT POUR LA VISU
C 
      CALL IMPET ('NOMBRE DE PAS DE TEMPS        ', NPICET)
      CALL IMPET ('NOMBRE DE COUCHES             ', NBCOU)
      CALL IMPET ('NOMBRE DE COLONNES            ', NBCOL)
      CALL IMPET ('NOMBRE D''ANGLES              ', NTETA)
      CALL IMPET ('PAS DE REPRISE AVEC PREEND    ', NUPATE)
      CALL ADTBDM ('ANGLES-COU', ADANG)
      CALL IMPTDT ('VALEUR DES ANGLES DE COUCHES ', DM(ADANG), 1, NBCOU)
      CALL IMPLT ('ENDOMMAGEMENT COUCHE          ', LENDCO)
      CALL IMPLT ('RUPTURE FIBRE                 ', LRUPCO)
      CALL IMPLT ('PLASTICITE COUCHE             ', LPLACO)
      CALL IMPLT ('ENDOMMAGEMENT INTERFACE       ', LENINT)
      CALL IMPLT ('PLASTICITE INTERFACE          ', LPLAIN)
C 
C    encours
C 
C    CD     CALL IMPLT ('modele en contrainte                 ', LMOSIG)
C    CD     CALL IMPLT ('modele avec retard                   ', LRETAR)
C    CD     CALL IMPLT ('modele avec retard apres instabilite ', LRETIN)
C    CD     CALL IMPLT ('deformation plastique effective      ', LPLAEF)
C    CD     CALL IMPLT ('deformation plastique finie          ', LPLEFF)
C 
C    fincours
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
900   CALL ERREUD (0, 'ERREUR D''OUVERTURE DANS '//IDPROG)
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine permet de modifier la taille des tableaux de stockage
C     des fonctions du temps et fonctions de l'espace pour le moment
C     uniquement en nombre de champs mais pas en nombre de piquets de temps.
C     Ceci permet de continuer une resolution lorsque l'on a  depasse le
C     nombre maxi de champs autorises. Cela permet egalement de ne garder
C     que les champs calcules pour la visu.
C 
C     On envoie comme arguments :
C 
C     E ...... NOUCHR   nouveau nombre de champs reels desires
C     E ...... NOUCHA   nouveau nombre de champs partiels desires
C 
C     Et on recupere :
C 
C     S ...... CHARAX  ancien nombre de champs reels desires
C     S ...... CHAMAX  ancien nombre de champs partiels desires
C 
C     ON SAUVEGARDE DETRUIT PUIS ON AGRANDIT LES TABLEAUX
C 
C     Pour les fonctions du temps definissant les accroissements
C     admissibles a zero
C 
C     E ...... NBFSIG  5 < = > evolution des contraintes
C     E ...... NBFEPS  6 < = > evolution des deformations
C 
C     Pour les fonctions du temps definissant la solution admissible
C 
C     E ...... NBDEPR  0 < = > deplacements totaux reels
C     E ...... NBDPTR  7 < = > evolution des deplacements totaux reels
C     E ...... DEADTR  8 < = > deformations admissibles totales reelles
C     E ...... SAADTR  9 < = > sauts admissibles totaux reels
C     E ...... EVCOTR 10 < = > evolution des quantites de type contraintes
C     E                        totales reelles
C     E ...... COTORE 11 < = > contraintes totales reelles
C     E ...... CNTORE 12 < = > contraintes normales totales reelles
C 
      SUBROUTINE MNBCHA (NOUCHR, NOUCHA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER   NOUCHR , NOUCHA
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER EPSCH8, SGNC12, SIGC11, SAUCH9, FTRE10
      INTEGER FTECH6, FTSCH5, FTREE7, DEPCH0
      INTEGER LONEPS, LONDER, LONSAU
      INTEGER LONRES, LONPRE, LONPLU, LONDIS
      INTEGER LONG(9,4) , NBNOM
      INTEGER LONRE2, DBSAUV, ADSAUV, POUB, ADEPF(9)
      INTEGER LONAVT(9), ADEPPA(9), LONTAB, TAB(2)
      INTEGER ACHARX, ACHAMX
      INTEGER I, J
      INTEGER AM2LC, ADM2LC
C 
      CHARACTER*10 LISNOM(9)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MNBCHA')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C    Tableaux des fonctions du temps pour les deltas admissibles
C 
      CALL ADTBDM ('TEMPS-SIGM', FTSCH5)
      CALL ADTBDM ('TEMPS-EPSI', FTECH6)
C 
C    (TYPE . EQ . 5) = > NBFSIG
C    (TYPE . EQ . 6) = > NBFEPS
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
      IF (NBCOU.GT.1) THEN
        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
      ELSE
        SAUCH9 = EPSCH8
        SGNC12 = SIGC11
      ENDIF
C 
      CALL ADTBDM ('TEMPS-REEL', FTREE7)
      CALL ADTBDM ('TEMP-SI-RE', FTRE10)
C 
      LONEPS = NTETA*NEPS*NGAU1
      LONDER = NDDL*NTETA
      LONSAU = NTETA*NSAU*NGAU2
C 
C    (TYPE .EQ.  0)= > NBDEPR
C    (TYPE .EQ.  7)= > NBDPTR
C    (TYPE .EQ.  8)= > DEADTR
C    (TYPE .EQ.  9)= > SAADTR
C    (TYPE .EQ. 10)= > EVCOTR
C    (TYPE .EQ. 11)= > COTORE
C    (TYPE .EQ. 12)= > CNTORE
C 
      CALL IMPET ('NBDPTR EN ENTREE DANS '//IDPROG, NBDPTR)
      CALL IMPET ('EVCOTR EN ENTREE DANS '//IDPROG, EVCOTR)
C 
C    Calcul de la place a reserver si celle ci est trop grande
C    par rapport a la place disponible; on cherche a segmenter.
C 
      LONRES = (NOUCHR+NOUCHR)*(NPICET+LONEPS+LONSAU)+NOUCHR*LONDER
     &        +2*(NOUCHA+NOUCHA)*NPICET
C 
C    Calcul de la place totale a prendre
C 
      LONPRE = 2*CHARAX*(NPICET+LONEPS+LONSAU)+CHARAX*LONDER
     &        +2*CHAMAX*NPICET
C 
      LONPLU = LONRES-LONPRE
C 
C     Calcul de la place disponible
C 
      LONDIS = ADM2-ADM1
C 
      IF (LONPLU .GT. LONDIS) THEN
        CALL IMPET ('LONGUEUR NECESSAIRE POUR LE REDIMENSIONNEMENT '//
     &              'DES FONCTIONS TEMPS-ESPACE ', LONPLU)
        CALL IMPET ('LONGUEUR DISPONIBLE DANS DM ', LONDIS)
        CALL ERREUD (0, 'REDIMENSIONNEMENT IMPOSSIBLE  '// IDPROG)
      END IF
C 
C     test pour savoir s'il faut segmenter
C 
      IF ((LONPRE+LONEPS) .GT. LONDIS) THEN
        CALL IMPET ('LONGUEUR NECESSAIRE POUR LE REDIMENSIONNEMENT '//
     &              'DIRECT DES FONCTIONS TEMS-ESPACE ', LONPRE+LONEPS)
        CALL IMPET ('LONGUEUR DISPONIBLE DANS DM ', LONDIS )
        CALL ERREUD (0, 'SEGMENTATION NON PROGRQMMEE : IMPOSSIBLE '
     &               //IDPROG)
      END IF
C 
      LISNOM(1) = 'TEMPS-SIGM'
      LONG(1,1) = NOUCHA
      LONG(1,2) = NPICET
      LONG(1,3) = NBFSIG
      LONG(1,4) = 5
      LISNOM(2) = 'TEMPS-EPSI'
      LONG(2,1) = NOUCHA
      LONG(2,2) = NPICET
      LONG(2,3) = NBFEPS
      LONG(2,4) = 6
      LISNOM(3) = 'TEMP-SI-RE'
      LONG(3,1) = NOUCHR
      LONG(3,2) = NPICET
      LONG(3,3) = EVCOTR
      LONG(3,4) = 10
      LISNOM(4) = 'TEMPS-REEL'
      LONG(4,1) = NOUCHR
      LONG(4,2) = NPICET
      LONG(4,3) = NBDPTR
      LONG(4,4) = 7
      LISNOM(5) = 'DEPLA-ADMI'
      LONG(5,1) = NOUCHR
      LONG(5,2) = LONDER
      LONG(5,3) = NBDEPR
      LONG(5,4) = 0
      LISNOM(6) = 'EPS-AD-TOT'
      LONG(6,1) = NOUCHR
      LONG(6,2) = LONEPS
      LONG(6,3) = DEADTR
      LONG(6,4) = 8
      LISNOM(7)  = 'SIG-AD-TOT'
      LONG(7,1) = NOUCHR
      LONG(7,2) = LONEPS
      LONG(7,3) = COTORE
      LONG(7,4) = 11
C 
      NBNOM     =7
C 
      IF (NBCOU .GT. 1) THEN
        LISNOM(8) = 'SAU-AD-TOT'
        LONG(8,1) = NOUCHR
        LONG(8,2) = LONSAU
        LONG(8,3) = SAADTR
        LONG(8,4) = 9
        LISNOM(9) = 'SGN-AD-TOT'
        LONG(9,1) = NOUCHR
        LONG(9,2) = LONSAU
        LONG(9,3) = CNTORE
        LONG(9,4) = 12
        NBNOM = 9
      ENDIF
C 
C     Remise a zero des nombres de fonctions du temps-espace
C 
      NBFSIG  = 0
      NBFEPS  = 0
      NBDEPR  = 0
      NBDPTR  = 0
      DEADTR  = 0
      SAADTR  = 0
      EVCOTR  = 0
      COTORE  = 0
      CNTORE  = 0
C 
      LONRE2 = LONPRE+LONEPS
C 
      IF((ADM2-LONRE2).LT.ADM1)THEN
         CALL IMPET ('ADM1 ', ADM1)
         CALL IMPET ('ADM2 ', ADM2)
         CALL IMPET ('LONGUEUR DU TABLEAU PROVISOIRE ', LONRE2)
         CALL ERREUD (0,
     &   'ECRITURE SUR DM AVANT LA 1ere ADRESSE LIBRE '//IDPROG)
      ENDIF
C 
      CALL POUSMD (LONRE2, DBSAUV)
      ADSAUV = DBSAUV
      POUB   = DBSAUV+LONPRE
C 
      DO I = 1, NBNOM
        CALL INFODP (LISNOM(I), ADEPPA(I), LONAVT(I))
        CALL COPITD (LONAVT(I), DM(ADEPPA(I)), DM(ADSAUV))
        ADSAUV = ADSAUV+LONAVT(I)
      END DO
C 
      CALL RETIDP (LISNOM, NBNOM)
C 
C     Dimensionnement aux nouvelles longueurs
C 
      ACHAMX = CHAMAX
      ACHARX = CHARAX
C 
      CHAMAX = NOUCHA
      CHARAX = NOUCHR
C 
      ADSAUV = DBSAUV
C 
C     Dimensionnement aux nouvelles longueurs
C 
      DO I = 1, 2
C 
        LONTAB = LONG(I,1)*LONG(I,2)
        CALL GESTDP (LISNOM(I), LONTAB, ADEPF(I))
C 
        TAB(1) = ACHAMX
        TAB(2) = LONG(I,2)
C 
        CALL IMPCT ('POUR LE TABLEAU         ', LISNOM(I))
        CALL IMPET ('ANCIENNE 1ere DIMENSION ', TAB(1))
        CALL IMPET ('NOUVELLE 1ere DIMENSION ', LONG(I,1))
        CALL IMPET ('NOUVELLE 2eme DIMENSION ', LONG(I,2))
        DO J = 1, LONG(I,3)
          CALL EXTRAD (DM(ADSAUV), 2, TAB(1),
     &                 2, J, DM(POUB), LONG(I,2))
          CALL RCHAMP (LONG(I,4), LONG(I,2), DM(POUB), DM(ADEPF(I)))
        END DO
C 
        ADSAUV = ADSAUV+LONAVT(I)
C 
      END DO
C 
      DO I = 3, NBNOM
C 
        LONTAB = LONG(I,1)*LONG(I,2)
        CALL GESTDP (LISNOM(I), LONTAB, ADEPF(I))
C 
        TAB(1) = ACHARX
        TAB(2) = LONG(I,2)
C 
        CALL IMPCT ('POUR LE TABLEAU         ', LISNOM(I))
        CALL IMPET ('ANCIENNE 1ere DIMENSION ', TAB(1))
        CALL IMPET ('NOUVELLE 1ere DIMENSION ', LONG(I,1))
        CALL IMPET ('NOUVELLE 2eme DIMENSION ', LONG(I,2))
C 
        DO J = 1, LONG(I,3)
          CALL EXTRAD (DM(ADSAUV), 2, TAB(1),
     &                 2, J, DM(POUB), LONG(I,2))
          CALL RCHARE (LONG(I,4), LONG(I,2), DM(POUB), DM(ADEPF(I)))
        END DO
C 
        ADSAUV = ADSAUV+LONAVT(I)
C 
      END DO
C 
      CALL IMPET ('NBFSIG EN SORTIE DANS '//IDPROG, NBFSIG)
      CALL IMPET ('NBFEPS EN SORTIE DANS '//IDPROG, NBFEPS)
      CALL IMPET ('NBDPTR EN SORTIE DANS '//IDPROG, NBDPTR)
      CALL IMPET ('EVCOTR EN SORTIE DANS '//IDPROG, EVCOTR)
C 
CD     DO I = 1, 4
CD       CALL IMPCN  ('POUR LE TABLEAU  ', LISNOM(I))
CD       CALL IMPEN  ('1ere DIMENSION   ', LONG(I,1))
CD       CALL IMPEN  ('2eme DIMENSION   ', LONG(I,2))
CD       CALL IMPEN  ('NOMBRE UTILISE   ', LONG(I,3))
CD       CALL IMPTDN ('VALEUR EN SORTIE ', DM(ADEPF(I)),
CD                     LONG(I,1), LONG(I,2))
CD    END DO
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
