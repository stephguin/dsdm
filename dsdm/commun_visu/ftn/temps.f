C     Cette routine doit servir pour le calcul des champs admissibles
C     et des deltas des champs admissibles.
C       Pour les champs admissibles             => nbftmx = charax
C       Pour les deltas des champs admissibles  => nbftmx = chamax
C 
C     Calcul des contraintes ou des deformations a un pas de temps donne
C
C  A LIRE ET VERIFIER
C
C ****************
C OA 10/01/2013
C     pour une  composante pour tout point de GAUSS 
C     
C     COORDONNES POLAIRES QUANTITES STOCKEES (NSIG OU NSAUT , NTETA, le reste)
C      Pour une quantite au NOEUD a verifier ATTENTION *** ( NTETA, Nbre ddl par C       noeuds = 6, le reste) 
C  
C 
C  Pour une version coordonnees cartesienne il existe VCHXYT ne semble marcher
C  que pour les quantités aux points de gauss stockes (NSIG, NTETA)
C ********************************************
C 
C     CHARAX est le nombre maximal de vecteurs crees au cours de l'execution,
C     initialisation comprise STOCKE (CHARAX, NTERME)
C 
C     Le resultat de cette routine est une fonction de l'espace rangee (nterme, npicet)
C 
C         VAL( t , m ) =   TEMPSi (NDEBFT, t) * CHAMPi (NDEBCH, m)
C                         - - - - - - - - - - - - - - - - - - -
C                        + TEMPSi (NFINFT, t) * CHAMPi (NFINCH, m)
C 
C         Une condition sine-qua-non est donc que
C         NFINCH-NDEBCH = NFINFT -NDEBFT
C 
C     ATTENTION !!! Dans champ, la valeur des quantites au point de Gauss
C 
C     Valeur des Champs REELS <=> CHARAX !!!!pour  Un Pas de TEmps
C 
      SUBROUTINE VCHAPT (NBFTMX, NTERME, NDEBFT, NFINFT, NDEBCH,
     &                   NFINCH, CHAMP, FTEMPS, NUPICE, VAL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C     Nombre de composantes des deformations ou contraintes
C 
      INTEGER NTERME
C 
C     Numero du piquet de temps pour lequel on veut faire  le calcul
C 
      INTEGER NUPICE
C 
C     Nombre de composantes des deformations ou contraintes
C 
      INTEGER NBFTMX
C 
C     Nombres de champs deja remplis
C 
      INTEGER NDEBFT , NFINFT , NDEBCH , NFINCH
C 
C     Tableau des variables du temps
C 
      DOUBLE PRECISION FTEMPS( NBFTMX * NPICET )
C 
C     Tableau des variables de l'espace
C 
      DOUBLE PRECISION CHAMP(CHARAX*NTERME)
C 
C     Resultat
C 
      DOUBLE PRECISION  VAL(NTERME)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  NBCHAR
      INTEGER  AM2LC, ADM2LC
C 
C     Indices de boucle et de ligne de vecteur
C 
      INTEGER TERME, NCHAMP, K, DEBUT, DEBTEM, VALTEM
      INTEGER DEBK, DCOEFF, KCOEFF, PICET, DEFOTE, DEVALT
C 
C     Coefficients
C 
      DOUBLE PRECISION COEFF, TDEB, VALACT
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VCHAPT')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Initialisation
C 
      CALL BALAID (NTERME, VAL)
C 
CD    CALL IMPEN ('NUPICE', NUPICE)
CD    IF ((NFINFT-NDEBFT) .NE. (NFINCH-NDEBCH)) THEN
CD       CALL IMPET ('NFINFT-NDEBFT DANS '//IDPROG, NFINFT-NDEBFT)
CD       CALL IMPET ('NFINCH-NFINFT DANS '//IDPROG, NFINCH-NDEBCH)
CD       CALL ERREUD ( 0 , 'ABSURDE VOYONS !!!!!')
CD    END IF
C 
C     DEBTEM est l'adresse precedant le debut des valeurs
C     interessantes dans le tableau des temps.
C 
C     La fonction du temps correspond a l'acroissement
C     des fonctions admissibles =>
C     Pour avoir leurs valeurs au temps correspondant
C     a nupice il faut les sommer jusqu'a  nupice -1.
C     De plus si le temps de depart de l'intervalle
C     actuellement pris en compte n'est pas zero il faut
C     retirer a la 1ere fonction du temps correspondant
C     a la solution elastique sa valeur au temps T(0).
C 
C     Determination des fonctions du temps donnees;
C     recherche des constantes de definition
C 
C     Reservation de place pour les valeurs des fonctions du temps
C 
      NBCHAR = NFINFT-NDEBFT+1
      CALL GSPOUD (NBCHAR, DCOEFF)
C 
C     Calcul de la valeur au temps correspondant a NUPICE
C     des evolutions admissibles
C 
C     Boucle sur tous les piquets de temps jusqu'a nupice
C 
      DO PICET = 1, NUPICE
        DEBTEM = NBFTMX*(PICET-1)+1
        KCOEFF = DCOEFF
        DO NCHAMP = NDEBFT, NFINFT
          DM(KCOEFF) = DM(KCOEFF)+FTEMPS(DEBTEM)
          DEBTEM     = DEBTEM+1
          KCOEFF     = KCOEFF+1
        END DO
      END DO
C 
C     Boucle sur tous les termes
C 
      DEBUT = 1
      DEBK  = NDEBCH
      DO TERME = 1, NTERME
        KCOEFF = DCOEFF
        K      = DEBK
C 
C       Boucle sur toutes les fonctions remplies
C 
        DO NCHAMP = NDEBCH, NFINCH
          COEFF      = DM(KCOEFF)
          KCOEFF     = KCOEFF+1
          VAL(DEBUT) = VAL(DEBUT)+COEFF*CHAMP(K)
          K          = K+1
        END DO
        DEBUT = DEBUT+1
        DEBK  = DEBK+CHARAX
      END DO
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine calcule la valeur des champs ou des accroissements des champs
C     admissibles pour la couche.
C 
C     NTERME           =>  nombre d'accroissements calcules en general NTERME=NEPS ou NSAU
C     CHAMPS           =>  NBFTMX = CHARAX
C     DELTA DES CHAMPS =>  NBFTMX = CHAMAX
C 
C     CACULS POUR UN POINT DE GAUSS DES TAUX DE CONTRAINTES OU DES
C     TAUX DE DEFORMATION A TOUS LES INSTANTS STOCKES X INTERV : (NEPS, NPICET)
C 
C     On envoie comme arguments :
C 
C     E ...... CHAEPS   partie du tableau des deformations admissibles
C                       correspondant au point de gauss stocke (CHARAX, NEPS)
C     E ...... CHASIG   partie du tableau des contraintes admissibles
C                       correspondant au point de gauss stocke (CHARAX, NSIG)
C     E ...... FTEPS    tableau des fonctions du temps pour les champs de deformations
C                       admissibles (FTEPMX, NPICET)
C     E ...... FTSIG    tableau des fonctions du temps pour les champs de contraintes
C                       admissibles (FTSIMX, NPICET)
C     E ...... NFEPS    nombre de fonctions du temps pour les champs de deformations admissibles
C     E ...... NFSIG    nombre de fonctions du temps pour les champs de contraintes admissibles
C 
C     Et on recupere :
C 
C     S ...... EPSTEM   accroissements des deformations admissibles pour tous
C                       les piquets de temps stockes (NEPS, NPICET)
C     S ...... SIGTEM   accroissements des contraintes admissibles pour tous
C                       les piquets de temps stockes (NSIG, NPICET)
C 
      SUBROUTINE ADCTPS (FTEPMX, FTSIMX, NTERME, CHAEPS, NDEBEP, NFINEP,
     &                   CHASIG, NDEBSI, NFINSI, FTEPS, NDEBFE, NFINFE,
     &                   FTSIG, NDEBFS, NFINFS,
C                        Et on recupere :
     &                   EPSTEM, SIGTEM)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER  FTEPMX, FTSIMX, NTERME,  NDEBEP, NFINEP
      INTEGER  NDEBSI, NFINSI
      INTEGER  NDEBFE, NFINFE
      INTEGER  NDEBFS, NFINFS
      DOUBLE PRECISION  CHAEPS(CHARAX*NEPS), CHASIG(CHARAX*NEPS)
      DOUBLE PRECISION  EPSTEM(NEPS*NPICET), SIGTEM(NEPS*NPICET)
      DOUBLE PRECISION  FTEPS(FTEPMX*NPICET), FTSIG(FTSIMX*NPICET)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ADCTPS')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      CALL VCPGTE (
C               Valeur des Champs au Point de Gauss pour tous les TEmps
C               On envoie :
     &            FTEPMX, NTERME, NDEBFE, NFINFE, NDEBEP, NFINEP,
     &            CHAEPS, FTEPS,
C               On recupere :
     &            EPSTEM)
C 
      CALL VCPGTE (
C               Valeur des Champs au Point de Gauss pour tous les TEmps
C               On envoie :
     &            FTSIMX, NTERME, NDEBFS, NFINFS, NDEBSI, NFINSI,
     &            CHASIG, FTSIG,
C               On recupere :
     &            SIGTEM)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine doit servir pour le calcul des champs admissibles
C     et des deltas des champs admissibles
C           pour les les champs admissibles         => nbftmx = charax
C           pour les deltas des champs admissibles  => nbftmx = chamax
C 
C     Calcul des contraintes ou des deformations pour tous les piquets donnes
C     pour NC composantes au point de GAUSS.
C 
C     CHARAX est le nombre maximal de vecteurs crees au cours de l'execution
C     initialisation comprise STOCKE (CHARAX, NTERME)
C 
C     Le resultat de cette routine est une fonction du temps et de l'espace
C     rangee (nterme, npicet)
C 
C     Somme de
C 
C           VAL(t, m) = TEMPSi(NDEBFT, t) * CHAMPi(NDEBCH, m)
C                   + TEMPSi(NFINFT, t) * CHAMPi(NFINCH, m)
C 
C     Une condition sinequanon est donc que NFINCH-NDEBCH = NFINFT -NDEBFT
C 
C     ATTENTION il faut envoyer dans champ la valeur des quantites au point
C     de Gauss. Si le tableau total en espace est stocke (NCMC, LONG) et que
C     l'on veut calculer l'evolution a partir de NUTERM < = LONG, il faut envoyer
C     dans champ le tableau total NCM*(NUTERM-1)+1.
C 
      SUBROUTINE VCPGTE (
C                    On envoie
     &                  NBFTMX, NTERME, NDEBFT, NFINFT, NDEBCH, NFINCH,
     &                  CHAMP, FTEMPS,
C                    On recupere
     &                  VAL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C     Nombre de composantes des deformations ou contraintes
C 
      INTEGER NTERME , NBFTMX
C 
C     Nombre de champs deja remplis
C 
      INTEGER NDEBFT , NFINFT , NDEBCH , NFINCH
C 
C     Tableau des variables du temps
C 
      DOUBLE PRECISION FTEMPS( NBFTMX * NPICET )
C 
C     Tableau des variables de l'espace
C 
      DOUBLE PRECISION CHAMP(CHARAX,NTERME)
C 
C     Resultat
C 
      DOUBLE PRECISION  VAL(NTERME*NPICET)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
C     Indice de boucle et de ligne de vecteur
C 
      INTEGER PICET, TERME, NCHAMP, K, DEBUT, DEBUK
C 
C     Coefficient
C 
      DOUBLE PRECISION COEFF
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VCPGTE')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C 
C     Initialisation
C 
      CALL BALAID (NTERME*NPICET, VAL)
C 
      DEBUK = NDEBFT
      DEBUT = 0
C 
C     BOUCLE SUR LES PIQUETS DE TEMPS
C 
      DO  PICET = 1 , NPICET
C 
C      BOUCLE i SUR LES FONCTIONS
C 
        K = DEBUK
C 
        DO NCHAMP = NDEBCH, NFINCH
C 
          COEFF =  FTEMPS(K)
          K     =  K+1
C 
C        BOUCLE iii SUR TOUS LES TERMES
C 
          DO TERME = 1, NTERME
C 
            VAL(DEBUT+TERME) = VAL(DEBUT+TERME)
     &                         + COEFF * CHAMP(NCHAMP, TERME)
C 
C        FIN DE BOUCLE iii SUR TOUS LES TERMES
C 
          END DO
C 
C      FIN DE BOUCLE ii SUR LES FONCTIONS
C 
        END DO
C 
        DEBUT = DEBUT + NTERME
        DEBUK = DEBUK + NBFTMX
C 
C     FIN DE BOUCLE SUR LES PIQUETS DE TEMPS
C 
      END DO
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Valeur des Champs REELS <=> CHARAX !!!!pour  Un Pas de TEmps
C 
C     On envoie comme arguments :
C 
C     E ...... NTERME  nombre de composantes des deformations ou contraintes
C     E ...... NBCHAR  nombre de champs deja remplis
C     E ...... NBFTMX  nombre de fonctions du temps remplies
C     E ...... CHAMP   fonctions de l'espace stockees (1, charax)
C     E ...... NUPICE  numero du piquet de temps pour lequel on veut faire le calcul
C 
C     Et on recupere :
C 
C     S ...... VAL     tableau des variables de l'espace resultat
C 
      SUBROUTINE VCUPTE (NTERME, NBCHAR, NBFTMX, CHAMP,
     &                   FTEMPS, NUPICE, VAL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C     nombre de composantes des deformations ou contraintes
C 
      INTEGER NTERME
C 
C     nombre de champs deja remplis
C 
      INTEGER NBCHAR
C 
C     nombre de fonctions du temps remplies
C 
      INTEGER NBFTMX
C 
C     numero du piquet de temps pour lequel on veut faire le calcul
C 
      INTEGER NUPICE
C 
C     tableau des variables du temps
C 
      DOUBLE PRECISION FTEMPS(CHARAX * NPICET)
C 
C     tableau des variables de l'espace
C 
      DOUBLE PRECISION CHAMP(CHARAX*NTERME)
C 
C     resultat
C 
      DOUBLE PRECISION  VAL(NTERME)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
C     indice de boucle et de ligne de vecteur
C 
      INTEGER TERME, NCHAMP, K, DEBUT, DEBTEM
      INTEGER DEBK, DCOEFF, KCOEFF, PICET, DEFOTE, DEVALT
C 
C     coefficient
C 
      DOUBLE PRECISION COEFF, TDEB, VALACT
C 
      INTEGER  AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VCUPTE')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL BALAID (NTERME, VAL)
C 
CD    CALL IMPEN  ('NUPICE', NUPICE)
CD    CALL IMPEN  ('NBCHAR', NBCHAR)
C 
C     CALL IMPTDP ('TEMPS', FTEMPS(1), CHARAX, NPICET)
C     CALL IMPTDP ('CHAMP', CHAMP(1), CHARAX, NTERME)
C 
C     DEBTEM est l'adresse precedant le debut des valeurs
C     interessantes dans le tableau des temps.
C 
C     La fonction du temps correspond a l'acroissement
C     des fonctions admissibles =>
C     Pour avoir leurs valeurs au temps correspondant
C     a nupice il faut les sommer jusqu'a nupice -1.
C     De plus si le temps de depart de l'intervalle
C     actuellement pris en compte n'est pas zero il faut
C     retirer a la 1ere fonction du temps correspondant
C     a la solution elastique sa valeur au temps T(0).
C 
C     Determination des fonctions du temps donnees
C     recherche des constantes de definition
C 
C     Calcul de la valeur au temps correspondant a NUPICE
C     de l'evolution elastique. => valact
C 
      CALL ADTBDM ('PARAM-CHAR', DEFOTE)
      CALL ADTBDM ('VALE-TEMPS', DEVALT)
C 
      TDEB    = DM(DEVALT+NUPICE) - DM(DEVALT)
C 
CD    CALL IMPDN ('TDEB',TDEB)
C 
      VALACT  = DM(DEFOTE+NUPICE)
C 
CD    CALL IMPDN ('VALACT', VALACT)
C 
C     Reservation de place pour les valeurs des fonctions du temps
C 
      CALL GSPOUD (NBCHAR, DCOEFF)
      DM(DCOEFF) = VALACT
C 
C     Calcul de la valeur au temps correspondant a NUPICE
C     des evolutions admissibles a zero.
C 
C     Boucle sur tous les piquets de temps jusqu'a nupice
C 
      DO PICET = 1, NUPICE
C 
        DEBTEM = NBFTMX*(PICET-1)+2
        KCOEFF = DCOEFF+1
C 
        DO NCHAMP = 2, NBCHAR
C 
          DM(KCOEFF)     = DM(KCOEFF)+FTEMPS(DEBTEM)
C 
          DEBTEM         = DEBTEM+1
          KCOEFF         = KCOEFF+1
C 
        END DO
C 
      END DO
C 
CD    CALL IMPTDN ('Valeurs des fonctions du temps pour  nupice ',
CD                  DM(DCOEFF), NBCHAR, 1)
C 
C     boucle sur tous les termes
C 
      DEBUT = 1
      DEBK  = 1
C 
      DO TERME = 1, NTERME
C 
        KCOEFF = DCOEFF
        K      = DEBK
C 
C       boucle sur toutes les fonctions remplies
C 
        DO NCHAMP = 1, NBCHAR
C 
          COEFF            = DM (KCOEFF)
          KCOEFF           = KCOEFF+1
C 
          VAL(DEBUT)       = VAL(DEBUT) +  COEFF * CHAMP(K)
C 
          K                = K+1
C 
        END DO
C 
        DEBUT = DEBUT+1
        DEBK  = DEBK+CHARAX
C 
      END DO
C 
C     CALL OMPTDP ('Valeur de la fonction au piquet de temps ',
C    &              VAL, NTERME, 1)
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
C     Cette routine range le nouveau champ calcule sous la forme (nbqt+1, loncha)
C 
C     On envoie comme arguments :
C 
C     E ...... TYPE      Type de quantite stockee [5, 6]
C 
C     A priori ne sert plus que pour les fonctions du temps
C     permettant de calculer les accroissements admissibles.
C 
C              TYPE=5 < = > evolution des contraintes
C              TYPE=6 < = > evolution des deformations
C      
C     E ...... LONCHA    longueur du terme a stocker
C     E ...... CHAMP     valeur du terme a stocker
C 
      SUBROUTINE RCHAMP (TYPE, LONCHA, CHAMP, TBSTOC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      DOUBLE PRECISION CHAMP(LONCHA), TBSTOC(CHAMAX, LONCHA)
C 
      INTEGER          TYPE, LONCHA
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   I, POSI, J
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='RCHAMP')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      IF (TYPE .EQ. 5) THEN
        POSI   = NBFSIG+1
        NBFSIG = MIN0(POSI, CHAMAX)
      ELSE IF (TYPE .EQ. 6) THEN
        POSI   = NBFEPS+1
        NBFEPS = MIN0(POSI, CHAMAX)
      ELSE
        CALL IMPET ('VALEUR DU TYPE ', TYPE)
        CALL ERREUD (0, 'MAUVAIS PASSAGE DE TYPE DANS '//IDPROG)
      END IF
C 
C     si posi > CHAMAX nombre max des delta on recopie dans tbstoc
C     les posi-1 deniers champs
C 
      IF (POSI .GT. CHAMAX) THEN
C 
        POSI = CHAMAX
C 
        DO J = 1, LONCHA
C 
          DO  I = 1, CHAMAX-1
C 
            TBSTOC(I, J) = TBSTOC(I+1, J)
C 
          END DO
C 
        END DO
C 
      END IF
C 
      DO I = 1, LONCHA
C 
        TBSTOC(POSI, I) = CHAMP(I)
C 
      END DO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine range le nouveau champ calcule sous la forme (nbqt+1, loncha)
C 
C     On envoie comme arguments :
C 
C     E ...... TYPE      Type de quantite stockee [0, 7]
C 
C                        0  < = > deplacement admissibles totaux reels
C                        7  < = > accroissement des deformations totales reelles
C                        8  < = > deformations admissibles totales reelles
C                        9  < = > sauts admissibles totaux reels
C                        10 < = > accroissement des quantites de type contraintes
C                                 totales reelles
C                        11 < = > contraintes totales reelles
C                        12 < = > contraintes normales totales reelles
C      
C     E ...... LONCHA    longueur du terme a stocker
C     E ...... CHAMP     valeur du terme a stocker
C 
      SUBROUTINE RCHARE (TYPECH, LONCHA, CHAMP, TBSTOC)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  CHAMP(LONCHA), TBSTOC(CHARAX, LONCHA)
      INTEGER           TYPECH, LONCHA
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     I, POSI
      CHARACTER*6 IDPROG
      PARAMETER  (IDPROG='RCHARE')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      IF (TYPECH .EQ. 0) THEN
        NBDEPR = NBDEPR+1
        POSI   = NBDEPR
        ELSE IF (TYPECH .EQ. 7) THEN
        NBDPTR = NBDPTR+1
        POSI   = NBDPTR
        ELSE IF (TYPECH .EQ. 8) THEN
        DEADTR = DEADTR+1
        POSI   = DEADTR
        ELSE IF (TYPECH .EQ. 9) THEN
        SAADTR = SAADTR+1
        POSI   = SAADTR
        ELSE IF (TYPECH .EQ. 10) THEN
        EVCOTR = EVCOTR+1
        POSI   = EVCOTR
        ELSE IF (TYPECH .EQ. 11) THEN
        COTORE = COTORE+1
        POSI   = COTORE
        ELSE IF (TYPECH .EQ. 12) THEN
        CNTORE = CNTORE+1
        POSI   = CNTORE
        ELSE
        CALL IMPET ('VALEUR DU TYPE ', TYPECH)
        CALL ERREUD (0, 'MAUVAIS PASSAGE DE TYPE DANS '//IDPROG)
      END IF
C 
      IF (POSI .GT. CHARAX) THEN
        CALL IMPET ('VALEUR DU TYPE ', TYPECH)
        CALL IMPET ('NOMBRE DE CHAMPS REELS MAXI ', CHARAX)
        CALL IMPET ('POSI > CHARAX ', POSI )
        CALL ERREUD (0, 'TROP GRAND NOMBRE DE '//IDPROG)
      END IF
C 
      DO I = 1, LONCHA
        TBSTOC (POSI, I) = CHAMP (I)
      END DO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine range les fonctions du temps admissibles
C     pour l'etape globale dans un tableau :
C 
C        'FT-EPS-'//NUETGL (NBDPTR, NPICET)
C        'FT-CON-'//NUETGL (EVCOTR, NPICET)
C 
C     On envoie comme arguments :
C 
C     E ...... FTEPS rangee (CHARAX, NPICET)
C     E ...... FTSIG rangee (CHARAX, NPICET)
C 
      SUBROUTINE ECFTGL (FTEPS, FTSIG)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      DOUBLE PRECISION FTEPS(CHARAX, NPICET)
      DOUBLE PRECISION FTSIG(CHARAX, NPICET)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   I, J, ADEPS, ADSIG, ADFEPS, ADFSIG
C 
      CHARACTER*3 CARETG
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ECFTGL')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL MESSAO ('ON ENTRE DANS '// IDPROG)
C 
C     NBDPTR est le nombre d'evolutions des deplacements totaux deja stockes
C 
      CALL IDENTI (NBETGL, CARETG)
      CALL GESTDP ('FT-EPS-'//CARETG, NBDPTR*NPICET, ADFEPS)
C 
      ADEPS = ADFEPS
      DO J = 1, NPICET
        DO I = 1, NBDPTR
          DM(ADEPS) = FTEPS(I, J)
          ADEPS     = ADEPS+1
       END DO
      END DO
C 
C     EVCOTR est le nombre d'evolutions des contraintes totales deja stockees
C 
      CALL GESTDP ('FT-CON-'//CARETG, EVCOTR*NPICET, ADFSIG)
C 
      ADSIG = ADFSIG
      DO J = 1, NPICET
        DO I = 1, EVCOTR
          DM(ADSIG) = FTSIG (I, J)
          ADSIG     = ADSIG+1
       END DO
      END DO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Valeur des derivees des champs REELS <=>
C     La derivee est estimee a partir de la moyenne de deux pas
C     de temps consecutifs en prolongeant a la fin
C 
C     On envoie comme arguments :
C 
C     E ...... NTERME  nombre de composantes des deformations ou contraintes
C     E ...... NBCHAR  nombre de champs deja remplis
C     E ...... NBFTMX  nombre de fonctions du temps remplies
C     E ...... CHAMP   fonctions de l'espace stockees (1, charax)
C     E ...... NUPICE  numero du piquet de temps pour lequel on veut faire le calcul
C 
C     Et on recupere :
C 
C     S ...... VAL     tableau des variables de l'espace resultat
C 
      SUBROUTINE VDERTE (NTERME, NBCHAR, NBFTMX, CHAMP,
     &                   FTEMPS, NUPICE, VAL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C     nombre de composantes des deformations ou contraintes
C 
      INTEGER NTERME
C 
C     nombre de champs deja remplis
C 
      INTEGER NBCHAR
C 
C     nombre de fonctions du temps remplies
C 
      INTEGER NBFTMX
C 
C     numero du piquet de temps pour lequel on veut faire  le calcul
C 
      INTEGER NUPICE
C 
C     tableau des variables du temps
C 
      DOUBLE PRECISION FTEMPS(CHARAX * NPICET)
C 
C     tableau des variables de l'espace
C 
      DOUBLE PRECISION CHAMP(CHARAX*NTERME)
C 
C     resultat
C 
      DOUBLE PRECISION  VAL(NTERME)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
C     indice de boucle et de ligne de vecteur
C 
      INTEGER TERME, NCHAMP, K, DEBUT
      INTEGER DEBK, DCOEFF, KCOEFF, PICET, DEVALT
      INTEGER DEBTE1
C 
C     coefficient
C 
      DOUBLE PRECISION COEFF, TDEB
C 
      INTEGER  AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VDERTE')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL BALAID (NTERME, VAL)
C 
C     DEBTEM est l'adresse precedant le debut des valeurs
C     interessantes dans le tableau des temps.
C 
C     La fonctions du temps correspond a l'acroissement
C     des fonctions admissibles =>
C     Pour avoir leurs valeurs au temps correspondant
C     a nupice il faut les sommer jusqu'a  nupice -1.
C     De plus si le temps de depart de l'intervalle
C     actuellement pris en compte n'est pas zero il faut
C     retirer a la 1ere fonction du temps correspondant
C     a la solution elastique sa valeur au temps T(0).
C 
C     Determination des fonctions du temps donnees
C     recherche des constantes de definition
C 
      CALL ADTBDM ('VALE-TEMPS', DEVALT)
C 
      TDEB    = DM(DEVALT+NUPICE) - DM(DEVALT)
C 
C     Reservation de place pour les valeurs des fonctions du temps
C 
      CALL GSPOUD (NBCHAR, DCOEFF)
C 
C     Recherche de la valeur au temps correspondant a NUPICE
C     des evolutions admissibles
C 
C 
      PICET =  NUPICE
C 
      KCOEFF = DCOEFF
C 
C     IF( NUPICE .NE. NPICET ) THEN
C 
C       DEBTE1 = NBFTMX*(PICET-1) + 1
C       DEBTE2 = NBFTMX*PICET     + 1
C 
C       DO NCHAMP = 1 , NBCHAR
C 
C          DM( KCOEFF )=.5D0* (FTEMPS(DEBTE1) + FTEMPS(DEBTE2))
C 
C          DEBTE1      = DEBTE1+1
C          DEBTE2      = DEBTE2+1
C          KCOEFF      = KCOEFF+1
C 
C        END DO
C 
C      ELSE
C 
       DEBTE1 = NBFTMX*(PICET-1)+1
C 
       DO NCHAMP = 1, NBCHAR
C 
          DM(KCOEFF) = FTEMPS(DEBTE1)
          DEBTE1     = DEBTE1+1
          KCOEFF     = KCOEFF+1
C 
       END DO
C 
C      END IF
C 
CD    CALL IMPTDN ('Valeurs de la derivee du temps pour nupice ',
CD                 DM(DCOEFF), NBCHAR, 1)
C 
C     boucle sur tous les termes
C 
      DEBUT = 1
      DEBK  = 1
C 
      DO TERME = 1, NTERME
C 
        KCOEFF = DCOEFF
        K      = DEBK
C 
C       boucle sur toutes les fonctions remplies
C 
        DO NCHAMP = 1, NBCHAR
C 
          COEFF            = DM ( KCOEFF )
          KCOEFF           = KCOEFF+1
          VAL(DEBUT)       = VAL(DEBUT) + COEFF * CHAMP(K)
          K                = K+1
C 
        END DO
C 
        DEBUT = DEBUT+1
        DEBK  = DEBK+CHARAX
C 
      END DO
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
