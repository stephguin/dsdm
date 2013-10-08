C     Nouvelle Etape locale avec resolution du comportement facon Newton modifiee.
C     Ecrit en contrainte tild mais le modele en deformation plastique est classique
C 
C     TYPE = 0, 1, 2, 3
C 
C     SI TYPE = 0 alors on est dans une etape locale classique
C                 ou les quantites admissibles sont calculees
C                 comme somme des quantites admissibles de l'etape
C                 locale precedente et des delta admissibles de l'etape
C                 globale precedente.
C 
C     SI TYPE = 1 alors on est dans une etape locale apres un appel a repcvr :
C                 la seule modification par rapport a une etape locale classique
C                 est que les champs admissibles sont recalcules completement
C                 a partir de ADCTPS
C 
C     SI TYPE = 2 alors on est dans une etape locale apres un appel a repcvr
C                 + visualisation
C 
C     SI TYPE = 3 alors on est dans une etape locale apres un appel a repcvr
C                 + visualisation uniquemement des points demandes; on ne
C                 sauvegarde rien; pour les couches il faut iterer aussi sur
C                 les points avec le meme Y.
C 
C     On calcule les quantites  :
C 
C         QCADMI : Les accroissement des quantites admissibles pour
C                  les couches, en deformations puis contraintes.
C                  (12*npicet) par enregistrement .
C 
C         QIADMI : Les accroissement des quantites admissibles pour
C                  les interfaces, en sauts puis contraintes normales
C                  (6*npicet) enregistrement.
C 
C     C'est egalement un nouveau type d'integration de facon a avoir un
C     endommagement constant dans l'epaisseur couche et retarde; on integre
C     egalement la plasticite.
C 
      SUBROUTINE LOCLOC (TYPE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'cominc_visu.h'
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
C     LOGIQUE SUIVANT TYPE
C 
      LOGICAL LTYP0, LTYP1, LTYP2, LTYP3
      LOGICAL LTYP01, LTYP12, LTYP23, LTY123
C 
C     UNITE
C 
      INTEGER IUNSIG
      INTEGER IUADMC, LONADM, LONSCH, IUADMN, IUADMI, IUACOU
      INTEGER LONPLA, LONEND, LONENP, NUENSU
      INTEGER NUENRS, NUENSI, NUENEN, NUENEP, NUENPL, NUENRN, NUENRI
      CHARACTER*20 NOM, NOMFIC
C 
C     TABLEAUX ADMISSIBLES
C 
      INTEGER EPSCH2, SIGCH1, SAUCH4, SGNCH3
C 
C     ADRESSES DE DEPART CORRESPONDANTES
C 
      INTEGER  ADEPS, ADSIG
C 
C     TABLEAUX DU TEMPS
C 
      INTEGER  FTSCH5, FTECH6, ADFOTE, ADINTE, DEVALT
      DOUBLE PRECISION INTERV, INTEIP
C 
C     TABLEAUX DES ANGLES
C 
      INTEGER  ADTETA
      DOUBLE PRECISION TETCAL, TETORT, ANGORT
C 
C     TABLEAUX DES POIDS
C 
      INTEGER  ADPOPG, DBPOPG
      DOUBLE PRECISION MULPPC(YINTEG)
      DOUBLE PRECISION MULPPI
C 
C     INDICES DE BOUCLE
C 
      INTEGER NUCOU, NUCOL, X, Y, TETA, TEMPS
C 
C     Pour le calcul de l'erreur locale
C 
      INTEGER  ERRTOT, NERRLC, DERRLC
      INTEGER  DCOULO, DINTLO
      DOUBLE PRECISION NUMERR, DENOMI, SIGCAO(6)
C 
C     Pour le passage en temps de Yd'
C 
      DOUBLE PRECISION YDPIN, YDPOUT
C 
C     Pour la determination du prochain nombre de picets
C 
      INTEGER ERRANC
C 
C     Pour les fichiers des champs admissibles
C 
      INTEGER QADMPR, QADMPS, QADMPN, DQADMP, DQADMN, QADMPI, DQADMI
C 
C     INDICATEUR COURANT
C 
      INTEGER LONPRO
      INTEGER RPRE, PPRE, EPSAPG
      INTEGER SIGAPG, DSIGAP
C 
C     INDICATEUR INITIAL
C 
      INTEGER DRPRE, DPPRE,  DPSPOR
      INTEGER DPSEOR
C 
C     Pour l'integration de la plasticite (base locale)
C 
      INTEGER DEPADM, EPADM
C 
C     Pour l'integration de l'endommagement (base locale)
C 
      INTEGER DMODUL, MODUL, DENDOM, ENDOM, LONRES, ENDCIN, ENDCSU
      INTEGER DENDCI, DENDCS
      LOGICAL NOCVD
      DOUBLE PRECISION YD, YDP, DIVISY
C 
C     Pour la visu
C 
      INTEGER NUCRI, NUCRC
C 
C     Pour l'initialisation
C 
C     Pour le non-lineaire
C 
      DOUBLE PRECISION   D2D3
      DOUBLE PRECISION   COFIBR(3), COPLAS(4), COENDO(12), ELATRI(9)
      DOUBLE PRECISION   CRITER(5)
      DOUBLE PRECISION   TR232, TR132, SRI13, SRI23, SAUIMP(3)
      DOUBLE PRECISION   TRC22, TRC12
      INTEGER            ADSORT 
C 
C     Pour le calcul des contraintes
C 
      INTEGER            SICAC, DSICAC, DSPCAC, BSPCAC
      INTEGER            ESPCAC, DEPCAC, DEPSAP
      DOUBLE PRECISION   SICOR(6), DSICAO(6)
      DOUBLE PRECISION   SORT(17), DEPCAO(6), EPSPVI(3)
C 
C     Pour la reprise
C 
C     Pour la visualisation
C 
      INTEGER NBPTIN, NUPDES
C 
C     Pour la visu de plusieurs etapes globales
C 
      INTEGER     NLUET, NBLUET, ADNUET, NBPDES
      CHARACTER*3 CARETG
      LOGICAL     LVISV
      DOUBLE PRECISION EPCLO(6), SIGCLO(6)
      DOUBLE PRECISION SIGALO(6)
C 
C     Pour la visualisation directe
C 
      INTEGER INDI13, INDI23,  NLU, INDI33
      INTEGER DBDEBU, DBCARX, DBDIRC, DBTABD
      LOGICAL LPASS3
      CHARACTER*3   CARCOU, CARCOL, CARGAX, CARTET, CARGAY, CARITE
C 
C     Pour la visu
C 
      INTEGER TABDES, IUNEND, IUNENP, IUNPLA
      INTEGER ADCARA, DBCARA, I, ADEBUT
      INTEGER INDIC
      INTEGER IUNMOD, NUENMO, LONMOD
      DOUBLE PRECISION ZERO(6), EPPCOR(6), EPCOR(6)
      DOUBLE PRECISION VMODA(2)
      LOGICAL LVISU
C 
C     Pour la visu des erreurs et des differences chapeaux admissibles
C 
      INTEGER DERRVI, ERRVIS, DDSGOR, DSIGOR, LONERR, LONDSG
      INTEGER NUERR , IUNERR, IUNDSG, NUDSG
      INTEGER LONCRI, IUNCRI, DDCRII, DCRIIN
      INTEGER LONCRC, IUNCRC, DDCRIC, DCRICO
C 
C     Pour les tests
C 
CD    DOUBLE PRECISION   ANSCSC, ADSCSC
C 
C     Pour l'interface
C 
      INTEGER NUGAU, NUINT, DEBGAU, FINGAU, PGAU1
      INTEGER PREEND, DPREEN
C 
      LOGICAL NOCVD1, NOCVD2
C 
      INTEGER ADSAUT, ADSIGN
      INTEGER SAUAPG, SGNAPG, SAUADM, SAUPOR
      INTEGER SGNCAC, SAUEOR,  DSGPAC, DSAPAC
      INTEGER DSGNAP, DSAUAP, DAUPOR, DSAADM
C 
      INTEGER DSGCAC, DSAEOR,  BSGPAC, BSAPAC
C 
      DOUBLE PRECISION A1, G1, G2, A, K0, K1, K2
      DOUBLE PRECISION R0, A2, K, N, YC
      DOUBLE PRECISION BE, AL, YAC, YACA, YACMIN, YACAMAX, YACMAX
C 
C     Pour la nouvelle version (voir david)
C 
      DOUBLE PRECISION  ALPINT, MINT, YOINT, YRINT
      DOUBLE PRECISION  YD1, YD2, YD3, INVALP
C 
      INTEGER           AM2LC, ADM2LC
C 
      DOUBLE PRECISION  SINOR(3), SINPOR(3), SINOP(3)
      DOUBLE PRECISION  SGAPOR(3), SGNAOR(3)
      DOUBLE PRECISION  SINSOM(3), DSGPOR(3)
      DOUBLE PRECISION  DSACAO(3), SAUCOR(3)
C 
C     Nouvelles declarations
C 
      INTEGER  EPSLAO, DEPLAO, VEPPXY, DVEPPX, VSIPXY, DVSIPX
      INTEGER  VSIPOR, DSIPOR, SVIPOR
      INTEGER  EPSAOR, DEPAOR, SIGAOR, DSIAOR
      INTEGER  EPSCOR, DEPCOR, SIGCOR, DSICOR
      INTEGER  EPPAOR, DEPPAO, SIPAOR, DSIPAO
      INTEGER  EPLAO, DADEND, EPSEOR
C      
C     Variable locale pour la visu LVIS3 en base polaire
C 
C     DOUBLE PRECISION  EPSAR0(6), SIGAR0(6), EPSCR0(6), SIGCR0(6)
      DOUBLE PRECISION  DPINF, DPSUP
      INTEGER           EPSAR0, SIGAR0, EPSCR0, SIGCR0
C  
C     Steph 27/05/99, pb avec QBORTH (travaille sur des tab de longueur 6)
C 
      INTEGER DEBPRO, BEBPRO, DSIPXY, LONLON
      INTEGER IUAVER, DEBUT, DEB, LONVER, VARVER
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='NETLOC')
      LOGICAL  LOGIERR
C 
C     Steph 03/11/99, pour le calcul implicite des contraintes : on garde
C     des informations du pas de temps precedent (en l'occurence, les retours en
C     compression)
C 
      LOGICAL     CASIN(3*YINTEG), CASOUT(3*YINTEG)
C 
      INTEGER     TUSED, DT, DTH, DTM, DTS, T1, T0, DEBIM1, DEBIM2
      LOGICAL     LOGIMP, RUPPRE
      INTEGER     DEBNET, SUPNET, IUADLE, NUENLE
      INTEGER     IUADEC, NUENEC, IUADIN, NUENIN
      INTEGER     DDEPSA, DDSIGA, DDEPSC, DDSIGC, DDEPLA, DDMODU
      CHARACTER*3 CARNET
C 
CD    CALL WLKBCD(IDPROG)
C 
       DO I =1, 6
         zero(i)=0.d0
       END DO
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     T0 = TUSED()
C 
C     Logique de type
C 
      LTYP0  = .FALSE.
      LTYP1  = .FALSE.
      LTYP2  = .FALSE.
      LTYP3  = .FALSE.
      LTYP01 = .FALSE.
      LTYP12 = .FALSE.
      LTYP23 = .FALSE.
      LTY123 = .FALSE.
C 
      IF (TYPE .EQ. 0) LTYP0 = .TRUE.
      IF (TYPE .EQ. 1) LTYP1 = .TRUE.
      IF (TYPE .EQ. 2) LTYP2 = .TRUE.
      IF (TYPE .EQ. 3) LTYP3 = .TRUE.
      IF ((TYPE .EQ. 0) .OR .(TYPE .EQ. 1)) LTYP01 = .TRUE.
      IF ((TYPE .EQ. 1) .OR .(TYPE .EQ. 2)) LTYP12 = .TRUE.
      IF ((TYPE .EQ. 2) .OR .(TYPE .EQ. 3)) LTYP23 = .TRUE.
      IF (LTYP12 .OR. LTYP3) LTY123 = .TRUE.
C 
C     Recherche des adresses des tableaux :
C 
C       - deformations admissibles (neps, nteta, ngau1, nfotps)
C       - contraintes admissibles (nsig, nteta, ngau1, nfotps)
C       - accroissements des fonctions du temps donnees (npicet, nbfodo)
C       - intervalles de temps (npicet)
C       - angles des bandes (nteta)
C       - numerateur erreur locale en temps (npicet+1)
C       - denominateur en temps de l'erreur locale en  temps (npicet+1)
C 
C     Tableaux d'initialisation des quantites locales pour le temps zero
C 
C     On distingue pour l'emploi de ADCTPS le cas de repise ou non
C 
C     NDEBEP est le numero du premier champ servant au calcul de
C     l'accroissement admissible a zero en deformation calcule a l'etape
C     globale precedente.
C     NDEBSI est le numero du premier champ servant au calcul de
C     l'accroissement admissible a zero en contrainte calcule a l'etape
C     globale precedente.
C 
C     Pour le calcul de l'accroissement en deformation, la premiere
C     adresse a servir dans le tableau des deformations reelles est
C     ADEPS, le decalage de NDEBEP se faisant dans ADCTPS.
C 
C     on recherche l'adresse de depart du tableau des fonctions
C     du temps servant au calcul des deltas admissibles a zero
C     => chamax fonctions maxis ou charax.
C 
      CALL ADTBDM ('EPS-AD-TOT', EPSCH2)
      CALL ADTBDM ('SIG-AD-TOT', SIGCH1)
C 
C     Pour les interfaces
C 
      IF (NBINT .GT. 0) THEN
         CALL ADTBDM ('SAU-AD-TOT', SAUCH4)
         CALL ADTBDM ('SGN-AD-TOT', SGNCH3)
      ENDIF
C 
      IF (LTYP0) THEN
C 
C     Cas 0 : on poursuit le calcul
C 
        CALL MESSAO ('ETAPE LOCALE TYPE 0')
	FTEPMX =  CHAMAX
        FTSIMX =  CHAMAX
        DBFTEP  = 1
        DBFTSI  = 1
        FIFTEP  = NBFEPS
        FIFTSI  = NBFSIG
        DBCHEP =  DEADTR-NBFEPS+1
        DBCHSI =  COTORE-NBFSIG+1
        FICHEP =  DEADTR
        FICHSI =  COTORE
C 
        CALL ADTBDM ('TEMPS-SIGM', FTSCH5)
        CALL ADTBDM ('TEMPS-EPSI', FTECH6)
C 
	DEBIM1 = FTSCH5
	DEBIM2 = FTECH6
C 
      ELSE IF (LTYP1) THEN
C 
C     Cas 1 : on reprend le calcul
C 
        CALL MESSAO ('ETAPE LOCALE TYPE 1')
        FTEPMX =  CHARAX
        FTSIMX =  CHARAX
        CALL ADTBDM ('TEMP-SI-RE', FTSCH5)
        CALL ADTBDM ('TEMPS-REEL', FTECH6)
C 
	DEBIM1 = FTSCH5
	DEBIM2 = FTECH6
C 
        DBFTEP  = 1
        DBFTSI  = 1
        FIFTEP  = NBDPTR
        FIFTSI  = EVCOTR
        DBCHEP =  1
        DBCHSI =  1
        FICHEP =  DEADTR
        FICHSI =  COTORE
C 
C     Cas 2 : on effectue une reprise dans le cadre d'une visualisation
C 
      ELSE IF (LTYP2) THEN
C 
        CALL MESSAO ('DONNEE DU NUMERO D''ETAPE GLOBALE A VISUALISER')
        CALL POUSME (1, ADNUET)
        CALL LECLEN (M(ADNUET), 1, NLUET)
	NUETGL = M(ADNUET)
        CALL IDENTI (M(ADNUET), CARETG)
        CALL INFODP ('FT-EPS-'//CARETG, FTECH6, FTEPMX)
        FTEPMX = FTEPMX/NPICET
        CALL INFODP ('FT-CON-'//CARETG, FTSCH5, FTSIMX)
        FTSIMX = FTSIMX/NPICET
C 
        DBFTEP  = 1
        DBFTSI  = 1
        FIFTEP  = FTEPMX
        FIFTSI  = FTSIMX
        DBCHEP =  1
        DBCHSI =  1
        FICHEP =  FTEPMX
        FICHSI =  FTSIMX
C 
C     Cas 2 : on effectue une reprise dans le cadre d'une visualisation
C 
      ELSE IF (LTYP2 .AND. (NBNETT .GT. 1)) THEN
C 
        CALL MESSAO ('DONNEE DU NUMERO D''ETAPE GLOBALE A VISUALISER')
        CALL POUSME (1, ADNUET)
        CALL LECLEN (M(ADNUET), 1, NLUET)
	NUETGL = M(ADNUET)
        CALL IDENTI (M(ADNUET), CARETG)
        FTEPMX =  CHARAX
        FTSIMX =  CHARAX
        CALL ADTBDM ('TEMP-SI-RE', FTSCH5)
        CALL ADTBDM ('TEMPS-REEL', FTECH6)
C 
        DBFTEP  = 1
        DBFTSI  = 1
        FIFTEP  = NBDPTR
        FIFTSI  = EVCOTR
        DBCHEP =  1
        DBCHSI =  1
        FICHEP =  DEADTR
        FICHSI =  COTORE
C 
      END IF
C 
      CALL ADTBDM ('F-TEMPS-DO', ADFOTE)
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
      CALL ADTBDM ('VALE-TEMPS', DEVALT)
C 
C     Pour le calcul de l'erreur totale
C     TR((siga-sigc) K-1(siga-sigc)) + ((signa-signc).K-1 ((signa-signc)
C     On construit un denominateur pour les couches (DE-COU-LOC) et un denominateur
C     pour les interfaces, afin de pouvoir visualiser les contributions relatives
C     couches/interfaces locales a l'indicateur d'erreur.
C 
      CALL ADTBDM ('NUMERA-LOC', NERRLC)
      CALL MENADM (NERRLC+1, NPICET)
      CALL ADTBDM ('DE-COU-LOC', DCOULO)
      CALL MENADM (DCOULO+1, NPICET)
      CALL ADTBDM ('DE-INT-LOC', DINTLO)
      CALL MENADM (DINTLO+1, NPICET)
      CALL ADTBDM ('DENOMI-LOC', DERRLC)
      CALL MENADM (DERRLC+1, NPICET)
      CALL ADTBDM ('ERREUR-LOC', ERRTOT)
C 
C     On va chercher les poids aux points de Gauss ranges (cou, col, x, y, teta)
C 
      CALL ADTBDM ('POIDS-PGCR', ADPOPG)
C 
C     Pour determiner le nombre de picet de la prochaine
C     etape locale, on conserve l'ancienne valeur de l'erreur.
C     Dans info -erreur on garde la valeur precedente de l'erreur
C     + le nombre de fonctionS du temps utiliseeS => 2
C 
      CALL ADTBDM ('INF-ERREUR', ERRANC)
      CALL COPITD (NPICET+1, DM(ERRTOT), DM(ERRANC))
      CALL MENADM (ERRTOT, NPICET+1)
C       
C -----------------------------------------------------------------------
C 
C     Creation de tableaux provisoires pour ranger toutes les quantites aux
C     piquets de temps.
C 
C     Accroissements des taux admissibles de l'etape globale precedente, X INTERV
C     EPSAPG en deformation            X      6*NPICET*YINTEG
C     SIGAPG en contrainte             X      6*NPICET*YINTEG
C 
      LONPRO = YINTEG*12*NPICET 
      CALL POUSMD (LONPRO, DEPSAP)
      DSIGAP = DEPSAP + 6*NPICET*YINTEG
C 
C     On commence par les quantites sauvegardees dans les fichiers de Q-CHAPEAU
C     admcou et sigmch (TYPE 0-1),
C     plasti, endomm, cricou, modul, erreur, dsigor (TYPE 2)
C 
C     Par rapport a l'ancienne version de etloca on recupere tous les points dans l'epaisseur.
C 
C     Les quantites du type n*(npicet+1) correspondent
C     a des valeurs a chaque piquet de temps et il s'agit de faire un decalage
C     de n apres chaque piquet de temps pour recopier la valeur qui vient d'etre
C     calculee comme nouvelle valeur precedente.
C 
C     Quantites sauvegardees dans le cas TYPE 0, 1 et 2
C     ______________________________________________
C 
C     DQADMP deltas admissibles precedents glob            X      12*NPICET*YINTEG
C     DQADMN copie de DQADMP                               X      12*NPICET*YINTEG
C     modif on en met 6 au lieu de 5 plus pratique
C     EPLAO  epsilon plastique ortho                       X   6*(NPICET+1)*YINTEG
C     VEPPXY delta epsilon point chapeau XDT               X       6*NPICET*YINTEG
C            eps-point-chapeau - eps-point-adm glo
C     VSIPXY delta sigma point chapeau XDT                 X       6*NPICET*YINTEG
C            sig-point-chapeau - sig-point-adm glo
C     ERRVIS denominateur de l'erreur au point             X         NPICET*YINTEG
C     DSIGOR sigchap-sigadm ortho                          X       6*NPICET*YINTEG
C     DDCRIC critere de couche                             X         NPICET*YINTEG
C     ENDOM  endommagement df, dps,dpt                     X          3*(NPICET+1)
C     MODUL  modules E11AC(Y) et E22AC(Y) et E33AC(Y)      X   3*(NPICET+1)*YINTEG
C                                                             __________________________
C 
C                                                         YINTEG*(53*NPICET+9)+3*(NPICET+1)
C 
      LONPRO = YINTEG*(53*NPICET+9)+3*(NPICET+1) 
C 
      CALL POUSMD (LONPRO, DQADMP)
C       
C     On calcule les adresses des debuts de tableau
C 
      DQADMN = DQADMP + 12*NPICET*YINTEG
      DEPLAO = DQADMN + 12*NPICET*YINTEG
      DVEPPX = DEPLAO + 6*(NPICET+1)*YINTEG
      DVSIPX = DVEPPX + 6*NPICET*YINTEG
      DERRVI = DVSIPX + 6*NPICET*YINTEG
      DDSGOR = DERRVI + NPICET*YINTEG
      DDCRIC = DDSGOR + 6*NPICET*YINTEG
      DENDOM = DDCRIC + NPICET*YINTEG
      DMODUL = DENDOM + 3*(NPICET+1)
C 
C     Pour les quantites a sauvegarder dans des fichiers pour le TYPE 3 seulement, on stocke la valeur
C     au pas de temps courant. Le stockage pour la visu est effectue dans sauvic (dans visu_dsdm,
C     dans delami c'est une routine bidon compilee a vide).
C 
C     Les valeurs admissible et chapeau sont a priori a initialiser a zero.
C 
C     Quantites sauvegardees dans le cas TYPE 3 seulement
C     ________________________________________________
C 
C     EPSAOR    epsilon admissible ortho                   X    6*YINTEG
C     SIGAOR    sigma admissible ortho                     X    6*YINTEG
C     EPSCOR    epsilon chapeau ortho                      X    6*YINTEG
C     SIGCOR    sigma chapeau ortho                        X    6*YINTEG
C                                                       ___________________
C 
C                                                              24*YINTEG
      LONPRO = 24*YINTEG
C 
C     On calcule les adresses des debuts de tableau
C 
      CALL POUSMD (LONPRO, DEPAOR)
      DSIAOR = DEPAOR+6*YINTEG
      DEPCOR = DSIAOR+6*YINTEG
      DSICOR = DEPCOR+6*YINTEG
C 
C     Reservation de place pour les quantites a passer d'un pas de temps a l'autre pour
C     l'integration de la loi de comportement, ce sont des variables d'entree-sortie pour INCOCO.
C 
C     RPRE    seuil plastique                              X      YINTEG
C     PPRE    deformation plastique cumulee                X      YINTEG
C     EPPAOR  epsilon point chapeau admissible orthotropie X    6*YINTEG 
C     SIPAOR  sigma point chapeau admissible orthotropie   X    6*YINTEG
C     EPSEOR  epsilon elastique orthotropie                X    6*YINTEG 
C     DSIPOR  sigma point orthotropie                      X    6*YINTEG 
C                                                       ___________________
C 
C                                                              26*YINTEG
      LONPRO = 26*YINTEG
C 
C     On calcule les adresses des debuts de tableau
C 
      CALL POUSMD (LONPRO, DRPRE)
      DPPRE  = DRPRE+YINTEG
      DEPPAO = DPPRE+YINTEG
      DSIPAO = DEPPAO+6*YINTEG
      EPSEOR = DSIPAO+6*YINTEG
      DSIPOR = EPSEOR+6*YINTEG
C 
C     Si en sequence de visualisation on ne veut pas
C     visualiser de quantite sur les couches
C 
      IF (LTYP2 .AND. (.NOT. LVISC)) GOTO 2000
C 
      IF (LTYP01) THEN
C 
        NOM = 'endomm'
        NOMFIC = NOM
        CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
        LONEND = 3*(NPICET+1)
        NUENEN = 1
        CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONEND, IUNEND)
C 
C       Ouverture du fichier pour lecture des quantites admissibles en poursuite
C 
        NOM = 'admcou'
        NOMFIC = NOM
        LONADM = 12*NPICET
        CALL OFDDNF (3, NOMFIC, 6, LONADM, IUADMC)
        NUENRS = 1
C  
C       Ouverture du fichier pour lecture des quantites admissibles en reprise
C 
	IF (LTYP1 .AND. (NBNETT .GE. 1)) THEN
	  CALL ADTBM ('LIS-BIGNET', DEBNET)
	  SUPNET = M(DEBNET+NBNETT-1)
	  CALL IDENTI (SUPNET, CARNET)
	  NOM = 'ac_'//CARNET
	  CALL OFDDNF (3, NOM, 6, LONADM, IUADLE)
	END IF
C 
C       Ouverture du fichier pour ecriture des quantites admissibles en reprise
C 
	IF (LTYP1 .AND. BIGNET) THEN
	  CALL ADTBM ('LIS-BIGNET', DEBNET)
	  SUPNET = M(DEBNET+NBNETT)
	  CALL IDENTI (SUPNET, CARNET)
	  NOM = 'ac_'//CARNET
	  CALL OFDDNF (3, NOM, 6, LONADM, IUADEC)
	END IF
C 
C 
C       Ouverture du fichier pour les contraintes delta_sigma_chapeau
C 
        NOM = 'sigmch'
        NOMFIC = NOM
        LONSCH = 6*NPICET
        NUENSI  = 1
        CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C 
      END IF
C 
      IF (LTYP2) THEN
C 
C       Ouverture du fichier pour lecture des quantites admissibles en poursuite
C 
        NOM = 'admcou'
        NOMFIC = NOM
        LONADM = 12*NPICET
        CALL OFDDNF (3, NOMFIC, 6, LONADM, IUADMC)
        NUENRS = 1
C  
C       Ouverture du fichier pour lecture des quantites admissibles en reprise
C 
	IF (NBNETT .GE. 1) THEN
	  CALL ADTBM ('LIS-BIGNET', DEBNET)
	  SUPNET = M(DEBNET+NBNETT-1)
	  CALL IDENTI (SUPNET, CARNET)
	  NOM = 'ac_'//CARNET
	  CALL OFDDNF (3, NOM, 6, LONADM, IUADLE)
	END IF
C 
        IF (LPLAST) THEN
          NOM = 'plasti'
          NOMFIC = NOM
          CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C        Longueur d'un enregistrement
C 
          NUENPL  = 1
          LONPLA  = 6*(NPICET+1)
          CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONPLA, IUNPLA)
        END IF
C 
        IF (LENDOM) THEN
          NOM = 'endomm'
          NOMFIC = NOM
          CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C        Longueur d'un enregistrement
C 
          LONEND = 3*(NPICET+1)
          NUENEN = 1
          CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONEND, IUNEND)
        END IF
C 
          IF (CRICOU) THEN
            NOM = 'cricou'
            NOMFIC = NOM
            CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C           Longueur d'un enregistrement
C 
            LONCRC= NPICET
            NUCRC = 1
            CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONCRC, IUNCRC)
          END IF
C 
        IF (LMODUL) THEN
          NOM = 'module'
          NOMFIC = NOM
          CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C        Longueur d'un enregistrement
C 
          LONMOD = 2*(NPICET+1)
          NUENMO  = 1
          CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONMOD, IUNMOD)
        END IF
C 
        IF (LVERRE) THEN
          NOM = 'erreur'
          NOMFIC = NOM
          CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C        Longueur d'un enregistrement
C 
          LONERR= NPICET
          NUERR  = 1
          CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONERR, IUNERR)
        END IF
C 
        IF (LDSIGO) THEN
          NOM = 'dsigor'
          NOMFIC = NOM
          CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C        Longueur d'un enregistrement
C 
          LONDSG= 6*NPICET
          NUDSG  = 1
          CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONDSG, IUNDSG)
        END IF
C 
      ENDIF
C 
C -----------------------------------------------------------------------
C 
C     SEQUENCE POUR LA VISU DIRECTE
C 
C -----------------------------------------------------------------------
C 
      IF (LTYP3) THEN
C 
        IF ((.NOT. LCPLAN) .AND. (.NOT. LCNORM) .AND. (.NOT. LCNONL))
     &  GOTO 2000
	LONPRO  = 7+YINTEG
        CALL POUSME (LONPRO, DBDEBU)
        DBCARA = DBDEBU+5
        DBCARX = DBCARA+1
        DBDIRC = DBCARX+1
C 
C       Pour visualiser les points suivants
C 
1000    CONTINUE
C 
        ADEBUT  = DBDEBU
        ADCARA  = DBCARA
C 
        CALL MESSAO 
     &             ('POINT A VISUALISER DANS LES COUCHES
     &              \ON PRECISE POUR LE POINT 5 DONNEES :
     &              \  1er   = NUMERO DE COUCHE
     &              \  2eme  = NUMERO DE COLONNE
     &              \  3eme  = NUMERO DE POINT DE GAUSS EN X
     &              \  4eme  = NUMERO DE L''ANGLE CONCERNE
     &              \  5eme  = NUMERO DE POINT DE GAUSS EN Y')
C 
        CALL LECLEN (M(ADEBUT), 5, NLU)
C 
C       S'IL N'Y A PAS DE POINT LU ON SORT
C 
        IF (NLU .EQ. 0) GOTO 2000
C 
        ADEBUT = ADEBUT-1
C 
        M(ADCARA)  =   (M(ADEBUT+1)-1)*NBCOL*XINTEG*NTETA*YINTEG
     &               + (M(ADEBUT+2)-1)*XINTEG*NTETA*YINTEG
     &               + (M(ADEBUT+3)-1)*NTETA*YINTEG
     &               + (M(ADEBUT+4)-1)*YINTEG
     &               +  M(ADEBUT+5)
C 
        CALL IMPET ('CARACTERISTIQUES DU POINT A VISUALISER ',
     &               M(DBCARA))
        CALL MESSAO ('NUMERO DES ETAPES GLOBALES A VISUALISER :
     &               \ON DONNE POUR LE POINT PRECEDENT DANS UN
     &               \ORDRE CROISSANT LES NUMEROS D''ETAPES
     &               \GLOBALES POUR ANALYSER LA CONVERGENCE ')
        CALL POUSME (NBETGL+1, ADNUET)
        CALL LECLEN (M(ADNUET), NBETGL+1, NLUET)
C 
C       TABDES adresse de depart du tableau provisoire des quantites a visualiser 
C       au point de Gauss courant sur les piquets de temps. Normalement appele 
C       uniquement dans le cas visu LTYP3. A VERIFIER.

        NBLUET = 0
        NUPDES = 0
        NBPDES = NLUET
        CALL GSPOUD (30*(NPICET+1)*NLUET, TABDES)
        DBTABD = TABDES
C 
        ADNUET = ADNUET-1
1110    CONTINUE
C 
        CALL MENADM (DEPAOR, 6*YINTEG)
        CALL MENADM (DSIAOR, 6*YINTEG)
        CALL MENADM (DEPCOR, 6*YINTEG)
        CALL MENADM (DSICOR, 6*YINTEG)
	CALL MENADM (DQADMP, 12*YINTEG)
	TABDES = DBTABD
        INDI13 = 0
        INDI23 = 0
        INDI33 = 0
C 
        ADCARA  = DBCARA
        ADNUET = ADNUET+1
C 
        CALL IDENTI (M(ADNUET), CARETG)
        CALL INFODP ('FT-EPS-'//CARETG, FTECH6, FTEPMX)
        FTEPMX = FTEPMX/NPICET
        CALL INFODP ('FT-CON-'//CARETG, FTSCH5, FTSIMX)
C 
        FTSIMX =  FTSIMX/NPICET
        DBFTEP =  1
        DBFTSI =  1
        FIFTEP =  FTEPMX
        FIFTSI =  FTSIMX
        DBCHEP =  1
        DBCHSI =  1
        FICHEP =  FTEPMX
        FICHSI =  FTSIMX
C 
C -----------------------------------------------------------------------
C 
C       FIN DE SEQUENCE POUR LA VISU DIRECTE
C 
C -----------------------------------------------------------------------
C 
      END IF
C 
C     remise a leurs valeurs initiales des adresses pour la visu
C 
      ADEPS  = EPSCH2
      ADSIG  = SIGCH1
      DBPOPG = ADPOPG
C 
      LONPRO = 6*YINTEG
      CALL POUSMD (LONPRO, DEBPRO)
C 
C ***********************************************************************
C 
C     1ERE BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE COUCHE
C 
C ***********************************************************************
C 
C     BOUCLE SUR LES COUCHES
C 
      DIVISY = DBLE(YINTEG)
      INDIC = 0
C 
      DO NUCOU= 1, NBCOU
C 
C       Recherche des caracteristiques du comportement de la couche
C 
        CALL CONLIC (D2D3, NUCOU,
     &               COFIBR, COPLAS, COENDO, ADSORT, ELATRI, CRITER)
C 
	TRC22=CRITER(1)*CRITER(1)
	TRC12=CRITER(1)*CRITER(1)
C 
C       Recherche de la souplesse dans la base d'orthotropie, 
C       de l'orientation de la couche
C 
        CALL SOCORT (NUCOU, SORT)
        CALL ANGCOU (NUCOU, TETORT)
C 
C       BOUCLE SUR LES COLONNES
C 
        DO NUCOL  = 1, NBCOL
C 
C         BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT X
C 
          DO X = 1, XINTEG
C 
	    DO Y = 1 , YINTEG
	      MULPPC(Y) = DM(DBPOPG)
	      DBPOPG = DBPOPG + 1
	    END DO 
C 	    
C           BOUCLE SUR LES ANGLES
C 
	    DO TETA = 1, NTETA
C 
C             Recherche de l'angle de la bande (tetcal) correspondant a teta
C 
              LOGIMP = .FALSE.
	      TETCAL = DM(ADTETA+TETA-1)
              ANGORT = TETCAL-TETORT
C 
C -----------------------------------------------------------------------
C 
C             SEQUENCE POUR LA VISU DIRECTE
C 
C -----------------------------------------------------------------------
C 
 	      IF (LTYP3) THEN
                LVISU  = .FALSE.
                LVISV  = .FALSE.
		INDI33 = 0
		INDI13 = 0
		DO Y = 1, YINTEG
		  INDI33 = INDI33+1
		  INDI23 = INDI23+1
		  IF (M(ADCARA) .EQ. INDI23) THEN
C 
C                   LE POINT EST A VISUALISER
C 
                    INDI13 = INDI13+1
                    NBLUET = NBLUET+1
                    NUPDES = NUPDES+1
                    LVISU = .TRUE.
                    IF (NBLUET .EQ. NLUET)  LVISV = .TRUE.
                    CALL IDENTI (NUCOU, CARCOU)
                    CALL IDENTI (NUCOL, CARCOL)
                    CALL IDENTI (X, CARGAX)
                    CALL IDENTI (TETA, CARTET)
                    CALL IDENTI (Y, CARGAY)
                    NOMSAU(1) = '   POINT DE GAUSS : '//
     &              CARCOU//', '//CARCOL//', '//CARGAX//', '//
     &              CARTET//', '//CARGAY
                    GOTO 1210
C 
                  END IF
C 
		END DO
C 
1210            CONTINUE
C 
		IF (INDI13 .EQ. 0) THEN
C 
C                 ON SAUTE LE POINT
C 
                  GOTO 1200
C 
                END IF
C 
              END IF
C 
              IF (LVISU .AND. (LTYP3)) THEN
C 
C               Sauvegarde pour la visu des quantites initiales
C 
                CALL BALAID (6, EPCOR)
                VMODA(1) = ELATRI(1)
                VMODA(2) = ELATRI(2)
		CALL SAUVIC (0, ZERO, ZERO, ZERO, ZERO, ZERO,
     &                       VMODA(1), ZERO, TABDES, NUPDES)
C 
              END IF
C -----------------------------------------------------------------------
C 
C               FIN DE SEQUENCE POUR LA VISU DIRECTE
C 
C -----------------------------------------------------------------------
C 
C             Pour relire les bonnes quantites admissibles
C             et pour ecrire tous les fichiers dans le meme sens
C             numeros d'enregistrement pour les quantites dependant de y
C 
              NUENRS = (NUCOU-1)*NBCOL*XINTEG*YINTEG*NTETA +
     &                 (NUCOL-1)*XINTEG*YINTEG*NTETA +
     &                 (X-1)*YINTEG*NTETA+TETA 
              NUENLE = NUENRS
	      NUENRN = NUENRS
	      NUENSI = NUENRS
              NUENPL = NUENRS
              NUERR  = NUENRS
              NUDSG  = NUENRS
C 
C             Numeros d'enregistrement pour les quantites ne dependant pas de y
C 
              NUCRC  = (NUCOU-1)*NBCOL*XINTEG*NTETA +
     &                 (NUCOL-1)*XINTEG*NTETA +
     &                 (X-1)*NTETA+TETA 
	      NUENMO = NUCRC
	      NUENEN = NUCRC
C 
C             Initialisation des adresses de depart pour le calcul des quantites
C             admissibles
C 	       
              ADEPS  = EPSCH2+CHARAX*NEPS*(NUENRS-1)
              ADSIG  = SIGCH1+CHARAX*NEPS*(NUENRS-1)
              QADMPR = DQADMP
              QADMPN = DQADMN
              EPSAPG = DEPSAP
              SIGAPG = DSIGAP
C 
C             BOUCLE SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Y
C             pour le calcul des quantites admissibles
C 
              DO Y = 1, YINTEG
C 
C               Reconstruction des accroissements des taux admissibles de l'etape globale
C               precedente, X INTERV
C 
                CALL ADCTPS (FTEPMX, FTSIMX, NEPS,
     &                       DM(ADEPS), DBCHEP, FICHEP,
     &                       DM(ADSIG), DBCHSI, FICHSI,
     &                       DM(FTECH6), DBFTEP, FIFTEP,
     &                       DM(FTSCH5), DBFTSI, FIFTSI,
C                            Et on recupere :
     &                       DM(EPSAPG), DM(SIGAPG))
C 
   		ADEPS = ADEPS+NTETA*NEPS*CHARAX
                ADSIG = ADSIG+NTETA*NSIG*CHARAX
C 
   		IF (LTYP0) THEN
		  CALL LFDDNF (DM(QADMPR), QADMPR, LONADM,
     &                         IUADMC, NUENRS)
                END IF
C 
                IF (LTYP1 .AND. (NBNETT .EQ. 0)) THEN
                  CALL MENADM (QADMPR, 12*NPICET)
                END IF
C   
  		IF (LTYP1 .AND. (NBNETT .GE. 1)) THEN
                  CALL LFDDNF (DM(QADMPR), QADMPR, LONADM,
     &                         IUADLE, NUENRS)
		END IF
C 
                IF (LTYP2 .OR. LTYP3) THEN
                  CALL MENADM (QADMPR, 12*NPICET)
                END IF
C   
   		IF (LTYP2 .AND. (NBNETT .GE. 1)) THEN
                  CALL LFDDNF (DM(QADMPR), QADMPR, LONADM,
     &                         IUADLE, NUENRS)
                  NUENRS = NUENRS + NTETA
   		END IF
C 
C               Calcul des accroissements de quantites admissibles
C               de 1 a 6  => deformations
C               de 7 a 12 => contraintes
C 
C               BOUCLE SUR LES PIQUETS DE TEMPS
C 
                DO TEMPS = 1, NPICET
C 
                  DO I = 1, 6
C 
C                   Calcul des taux de deformation admissible
C 
                    DM(QADMPR) = DM(EPSAPG) + DM(QADMPR)
                    QADMPR = QADMPR+1
		    EPSAPG = EPSAPG+1
C 
                  END DO
C 
C                 Calcul des taux de contrainte admissible
C 
                  DO I = 1, 6
                    DM(QADMPN) = DM(QADMPR)
                    DM(QADMPR) = DM(SIGAPG) + DM(QADMPR)
                    QADMPR = QADMPR+1
		    SIGAPG = SIGAPG+1
                  END DO
C 
C               FIN DE BOUCLE SUR LES PIQUETS DE TEMPS
C 	      
		END DO 
C 		 
		QADMPR = QADMPR-12*NPICET
C 
                IF (LTYP0) THEN
                  CALL EFDDNF (IUADMC, DM(QADMPR),
     &                         QADMPR, LONADM, NUENRS)
                END IF
                IF (LTYP1 .AND. BIGNET) THEN
                  CALL EFDDNF (IUADEC, DM(QADMPR),
     &                         QADMPR, LONADM, NUENRS)
                  CALL EFDDNF (IUADMC, DM(QADMPR),
     &                         QADMPR, LONADM, NUENRS)
                END IF
                IF (LTYP1 .AND. (.NOT. BIGNET)) THEN
                  CALL EFDDNF (IUADMC, DM(QADMPR),
     &                         QADMPR, LONADM, NUENRS)
                END IF
C 
                NUENRS = NUENRS+NTETA
	        QADMPR = QADMPR+12*NPICET
		QADMPN = QADMPN+12*NPICET
C 
C             FIN DE BOUCLE SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT Y
C 
	      END DO 
C 
C             CHANGEMENT DE QADMPR (12, NPICET, YINTEG) EN QADMPR (12, YINTEG, NPICET)
C             Necessaire pour lire point de gauss par point de gauss dans la boucle
C             en temps
C 
   	      CALL TAB132 (12, NPICET, YINTEG, DM(DQADMP))
C 
C             Initialisation des adresses de depart
C 
              EPLAO  = DEPLAO
	      VEPPXY = DVEPPX
	      VSIPXY = DVSIPX
              ERRVIS = DERRVI
              DSIGOR = DDSGOR
	      DCRICO = DDCRIC
	      ENDOM  = DENDOM
	      MODUL  = DMODUL
	      EPSAOR = DEPAOR
	      SIGAOR = DSIAOR
	      EPSCOR = DEPCOR
	      SIGCOR = DSICOR	      
              RPRE   = DRPRE
              PPRE   = DPPRE
              QADMPR = DQADMP
C 
C             Initialisation des champs avec decalage et des valeurs entree sortie 
C             sauf eppaor et sipaor calculees a partir de QADMPR
C 
	      DM(ENDOM)   = 0.D0
              DM(ENDOM+1) = 0.D0
              DM(ENDOM+2) = 0.D0
	      DEB = 0
	      DO Y = 1, YINTEG
	        DM(MODUL)   = ELATRI(1)
		MODUL = MODUL+1
                DM(MODUL) = ELATRI(2)
		MODUL = MODUL+1
	        DM(MODUL) = ELATRI(3)
		MODUL = MODUL+1
                DM(RPRE)    = COPLAS(1)
		RPRE = RPRE+1
                DM(PPRE)    = 0.D0
		PPRE = PPRE+1
		DEB = DEB+3
	      END DO
	      MODUL  = DMODUL
C 
C             Initialisations de la boucle en temps
C 
              CALL MENADM (EPLAO , 6*YINTEG)
              CALL MENADM (EPSAOR, 6*YINTEG)
              CALL MENADM (SIGAOR, 6*YINTEG)
	      CALL MENADM (EPSCOR, 6*YINTEG)
	      CALL MENADM (SIGCOR, 6*YINTEG)
	      CALL MENADM (EPSEOR, 6*YINTEG)
C             TEST A RETIRER A PRIORI
  	      CALL MENADM (VSIPXY, 6*YINTEG*NPICET)
	      CALL MENADM (DSIGOR, 6*YINTEG)
	      NUMERR = 0.D0
	      DENOMI = 0.D0
C 
	      RUPPRE = .FALSE.
	      DO I = 1, (3*YINTEG)
	        CASIN(I)  = .TRUE.
	        CASOUT(I) = .TRUE.
	      END DO
   	      YDPIN = 0.D0
C 
C             BOUCLE SUR LES PIQUETS DE TEMPS
C 
	      LOGIMP = .FALSE.
CD 	      IF (LOGIMP) THEN
CD 	        CALL IMPET ('COUCHE  ', NUCOU)
CD              CALL IMPET ('COLONNE ', NUCOL)
CD              CALL IMPET ('X       ', X)
CD              CALL IMPET ('TETA    ', TETA)
CD 	      END IF
C 
	      DO  TEMPS = 1, NPICET
C 
CD              IF ((NUCOU .EQ. 8) .AND. (NUCOL .EQ. 1) .AND. (X .EQ. 1)
CD   &               .AND. (TETA .EQ. 2) .AND. (TEMPS .EQ. 10)) THEN
CD   	          LOGIMP = .TRUE.
CD 		END IF
		LOGIERR = .FALSE.
		INTERV  = DM(ADINTE+TEMPS-1)
	        VSIPOR = DSIPOR
	        CALL MENADM (VSIPOR, 6*YINTEG)
	        EPPAOR = DEPPAO
                EPSAOR = DEPAOR
		SIPAOR = DSIPAO
		SIGAOR = DSIAOR
C 
C               Dans QADMPR il y a les quantites admissibles dans la base globale
C               rangees (12, YINTEG, NPICET)
C 
                DO Y = 1, YINTEG  
                  CALL QBORTH (ANGORT, DM(QADMPR), DM(EPPAOR))
                  QADMPR = QADMPR+6
		  EPPAOR = EPPAOR+6
                  CALL QBORTH (ANGORT, DM(QADMPR), DM(SIPAOR))
		  QADMPR = QADMPR+6
		  SIPAOR = SIPAOR+6
		END DO
C 
CD 		CALL IMPET ('TEMPS   ', TEMPS)
		CALL INCOCO (INTERV, MULPPC,
C Parametres materiau
     &		     D2D3, COFIBR, COPLAS, COENDO, DM(ADSORT), ELATRI,
C Taux admissibles dans la base d'orthotropie en Entree (X DELTA T)
     &               DM(DEPPAO), DM(DSIPAO),
C Valeurs admissibles dans la base d'orthotropie (non stockees sont en Entree-Sortie)     
     &               DM(EPSAOR), DM(SIGAOR), 
C Quantites chapeau dans la base d'orthotropie (non stockees sont en Entree-Sortie)     
     &               DM(EPSCOR), DM(SIGCOR), DM(DRPRE), DM(DPPRE),
     &               DM(EPSEOR),
C Quantites chapeau dans la base d'orthotropie stockee (en Entree) 
     &               DM(ENDOM), DM(MODUL), CASIN, YDPIN,
     &               DM(DCRICO), DM(EPLAO),
C Quantites chapeau dans la base d'orthotropie stockee (en Sortie) 
     &               DM(ENDOM+3), DM(MODUL+3*YINTEG), CASOUT, YDPOUT,
     &               DM(DCRICO+1), DM(EPLAO+6*YINTEG),
C Quantites erreur dans la base d'orthotropie (en Sortie) seront tournees apres incoco     
     &               DM(ERRVIS), DM(VEPPXY), DM(VSIPOR), DM(DSIGOR),
     &               NUMERR, DENOMI, RUPPRE, LOGIMP)
C 
		DM(NERRLC+TEMPS) = DM(NERRLC+TEMPS) + NUMERR
                DM(DERRLC+TEMPS) = DM(DERRLC+TEMPS) + DENOMI
                DM(DCOULO+TEMPS) = DM(DCOULO+TEMPS) + DENOMI
C 
CD 		CALL IMPTDT ('VERIF PLAST ', DM(EPLAO+6*YINTEG), 1, 6)
		DO I = 1, (3*YINTEG)
CD 		  CALL IMPET ('ZBLB ', I)
CD 		  IF (CASIN(I) .NE. CASOUT(I)) THEN
CD 		    IF (CASIN(I) .AND. (.NOT. CASOUT(I)) .AND.
CD   &                  (TEMPS .NE. 1)) THEN
CD 		      CALL MESSAO ('RETOUR EN COMPRESSION')
CD 		    END IF
CD 		    IF (CASIN(I) .AND. (.NOT. CASOUT(I)) .AND.
CD   &                  (TEMPS .EQ. 1)) THEN
CD 		      CALL MESSAO ('INIT EN COMPRESSION')
CD 		    END IF
CD 		    IF ((.NOT. CASIN(I)) .AND. CASOUT(I) .AND.
CD   &                  (TEMPS .NE. 1)) THEN
CD 		      CALL MESSAO ('SUPER ZARBI')
CD 		    END IF
CD 		  END IF
		  CASIN(I) = CASOUT(I)
		END DO
   		YDPIN = YDPOUT
C 
		IF (LOGIERR) THEN
		  CALL MESSAO ('PB CONVERGENCE '//IDPROG)
		  CALL IMPET ('COUCHE  ', NUCOU)
		  CALL IMPET ('COLONNE ', NUCOL)
		  CALL IMPET ('GAUSS X ', X)
		  CALL IMPET ('TETA    ', TETA)
		  CALL IMPET ('TEMPS   ', TEMPS)
		END IF
C 
C               Mise dans les quantites liees a l'iteration precedente des resultats
C               de l'iteration actuelle
C 
		DO Y = 1, YINTEG
		  CALL QBORTH (-ANGORT, DM(VSIPOR), DM(VSIPXY))
		  VSIPOR = VSIPOR + 6
		  VSIPXY = VSIPXY + 6
		END DO
C 
                IF (LVISU .AND. LTYP3) THEN
                  IF (LOGR0) THEN
                    CALL QBORTH (-ANGORT, DM(EPSAOR), DM(EPSAR0))
                    CALL QBORTH (-ANGORT, DM(SIGAOR), DM(SIGAR0))
                    CALL QBORTH (-ANGORT, DM(EPSCOR), DM(EPSCR0))
                    CALL QBORTH (-ANGORT, DM(SIGCOR), DM(SIGCR0))
         		    EPSAR0 = EPSAR0+((INDI33-1)*6)
        		    SIGAR0 = SIGAR0+((INDI33-1)*6)
        		    EPSCR0 = EPSCR0+((INDI33-1)*6)
        		    SIGCR0 = SIGCR0+((INDI33-1)*6)
        		    EPLAO  = EPLAO+((INDI33-1)*6)
                 	    DDMODU = MODUL+((INDI33-1)*3)
		    CALL SAUVIC (TEMPS, DM(EPSAR0), DM(SIGAR0),
     &	                         DM(EPSCR0), DM(SIGCR0),
     &                           DM(ENDOM), DM(DDMODU), DM(EPLAO),
     &                           TABDES, NUPDES)
                  ELSE
   		    DDEPSA = DEPAOR+((INDI33-1)*6)
   		    DDSIGA = DSIAOR+((INDI33-1)*6)
   		    DDEPSC = DEPCOR+((INDI33-1)*6)
   		    DDSIGC = DSICOR+((INDI33-1)*6)
   		    DDEPLA = EPLAO+((INDI33-1)*6)
		    DDMODU = MODUL+((INDI33-1)*3)
                    CALL SAUVIC (TEMPS, DM(DDEPSA), DM(DDSIGA),
     &		                 DM(DDEPSC), DM(DDSIGC), DM(ENDOM),
     &		                 DM(DDMODU), DM(DDEPLA), TABDES, NUPDES)
                  END IF
                END IF
C s
                ENDOM  = ENDOM  + 3
                MODUL  = MODUL  + 3*YINTEG
                DCRICO = DCRICO + 1
                EPLAO  = EPLAO  + 6*YINTEG
                ERRVIS = ERRVIS + YINTEG
                VEPPXY = VEPPXY + 6*YINTEG
                DSIGOR = DSIGOR + 6*YINTEG
C 
C             FIN DE BOUCLE SUR LES PIQUETS DE TEMPS
C 
	      END DO
C 
              LOGIMP = .FALSE.
	      IF (LTYP01) THEN
		CALL EFDDNF (IUNEND, DM(DENDOM),
     &                       DENDOM, LONEND, NUENEN)
              ENDIF
C 
              IF (LTYP2) THEN
C 
                IF (LMODUL) THEN
C 
C                 Ecriture des modules 11 et 22
C 
                  CALL EFDDNF (IUNMOD, DM(DMODUL),
     &                         DMODUL, LONMOD, NUENMO)
                END IF
C 
                IF (LENDOM) THEN
C 
C                 Ecriture des endommagements
C 
CD                CALL IMPTDT ('VERIF COUCHE ', DM(DENDOM), 1, LONEND)
		  CALL EFDDNF (IUNEND, DM(DENDOM),
     &                         DENDOM, LONEND, NUENEN)
                ENDIF
C 
		IF (CRICOU) THEN
C 
C                 Ecriture du critere elastique
C 
		  CALL EFDDNF (IUNCRC, DM(DDCRIC),
     &                         DDCRIC, LONCRC, NUCRC)
                ENDIF
C 
		IF (LPLAST) THEN
C 
C                 Ecriture des deformations plastiques
C 
		  CALL TAB132 (6, YINTEG, NPICET+1, DM(DEPLAO))
		  CALL POUSMD (LONPLA, DEBUT)
		  DEB = DEBUT
		  EPLAO = DEPLAO
                  DO Y = 1, YINTEG
		    DEBUT = DEB		   	
   		    DO I = 1, NPICET+1
		      DM(DEBUT) = DM(EPLAO)
		      DEBUT     = DEBUT+1
		      DM(DEBUT) = DM(EPLAO+1)
		      DEBUT     = DEBUT+1
  		      DM(DEBUT) = DM(EPLAO+2)
  		      DEBUT     = DEBUT+1
  		      DM(DEBUT) = DM(EPLAO+3)
  		      DEBUT     = DEBUT+1
  		      DM(DEBUT) = DM(EPLAO+4)
  		      DEBUT     = DEBUT+1
  		      DM(DEBUT) = DM(EPLAO+5)
  		      DEBUT     = DEBUT+1
   		      EPLAO     = EPLAO+6
		    END DO
                    CALL EFDDNF (IUNPLA, DM(DEB), DEB, LONPLA, NUENPL)
                    EPLAO  = DEPLAO + 6*(NPICET+1)
                    NUENPL = NUENPL + NTETA
		  END DO 
                END IF
C 
                IF (LVERRE) THEN
C 
C                 Ecriture de la contribution locale de l'erreur
C 
                  CALL TAB132 (1, YINTEG, NPICET, DM(DERRVI))
		  ERRVIS = DERRVI 
                  DO Y = 1, YINTEG 
                    CALL EFDDNF (IUNERR, DM(ERRVIS),
     &                           ERRVIS, LONERR, NUERR)
CD                  CALL IMPTDT ('VERIF ECRITURE ', DM(ERRVIS),
CD   &                            1, LONERR)
		    ERRVIS = ERRVIS+NPICET
                    NUERR = NUERR + NTETA
		  END DO 		    
                END IF
C 
                IF (LDSIGO) THEN
C 
C                 Ecriture de l'ecart chapeau-admissible
C 
                  CALL TAB132 (6, YINTEG, NPICET, DM(DDSGOR))
		  DSIGOR = DDSGOR
                  DO Y = 1, YINTEG 
                    CALL EFDDNF (IUNDSG, DM(DSIGOR),
     &                           DSIGOR, LONDSG, NUDSG)
                    NUDSG  = NUDSG + NTETA
		    DSIGOR = DSIGOR+ 6*NPICET
		  END DO 
                END IF
C 
              END IF
C 
1200          CONTINUE
C 
C             CHOIXD est compilee a vide dans delami, mais appelle des routines
C             graphiques dans visu_dsdm (---> repertoire faux_commun_visu)
C 
              IF (LVISV .AND. LTYP3) THEN
                CALL CHOIXD (3, NBPDES, DM(TABDES), DM(DEVALT))
              END IF
C 
C             Pour avoir directement les contraintes permettant le calcul
C             des seconds membres pour l'etape globale
C 
C             <=> on stocke - delta (sigma-point-chapeau)
C  PAS FAIT   <=> on stocke   delta (eps-point-chapeau)
C 
              IF (LTYP01) THEN
                CALL TAB132 (6, YINTEG, NPICET, DM(DVSIPX))
		VSIPXY = DVSIPX
       		DO Y = 1, YINTEG
                  CALL EFDDNF (IUNSIG, DM(VSIPXY),
     &                         VSIPXY, LONSCH, NUENSI)
                  VSIPXY = VSIPXY + 6*NPICET
                  NUENSI = NUENSI + NTETA
		END DO 
              END IF
C 
C           FIN DE BOUCLE SUR LES ANGLES 
C 
	    END DO
C 
1300        CONTINUE
C 
C         FIN DE BOUCLE ii SUR LES POINTS DE GAUSS DANS L'ELEMENT SUIVANT X
C 
          END DO
C 
C       FIN DE BOUCLE i SUR LES COLONNES
C 
        END DO
C 
C     FIN DE BOUCLE SUR LES COUCHES
C 
      END DO
C 
      IF (LTYP2) THEN
        IF (LPLAST) CALL FERFIC (3, IUNPLA, IDPROG)
        IF (LENDOM) CALL FERFIC (3, IUNEND, IDPROG)
	IF (CRICOU) CALL FERFIC (3, IUNCRC, IDPROG)	  
        IF (LMODUL) CALL FERFIC (3, IUNMOD, IDPROG)
        IF (LVERRE) CALL FERFIC (3, IUNERR, IDPROG)
        IF (LDSIGO) CALL FERFIC (3, IUNDSG, IDPROG)
      END IF
C 
      IF (LTYP01) THEN
        CALL FERFIC (3, IUADMC, IDPROG)
        CALL FERFIC (3, IUNSIG, IDPROG)
	CALL FERFIC (3, IUNEND, IDPROG)
      END IF
      IF (LTYP1 .AND. (NBNETT .GE. 1)) THEN
	CALL FERFIC (3, IUADLE, IDPROG)
      END IF
      IF (LTYP1 .AND. BIGNET) THEN
	CALL FERFIC (3, IUADEC, IDPROG)
      END IF
C 
      IF (LTYP2) THEN
        CALL FERFIC (3, IUADMC, IDPROG)
      END IF
      IF (LTYP2 .AND. (NBNETT .GE. 1)) THEN
        CALL FERFIC (3, IUADLE, IDPROG)
      END IF
C 
C     Pour l'etape suivante en visu directe
C     Pour un autre point de Gauss en visu directe
C     Sinon, fin de la visu des points de couche.
C 
      IF (LTYP3) THEN
        IF (NBLUET .NE. NLUET) THEN
	  GOTO 1110
	ELSE
	  GOTO 1000
	END IF
      END IF
C 
2000  CONTINUE
C 
C     Si en sequence de visualisation on ne veut pas
C     visualiser de quantites sur les interfaces
C 
      IF (LTYP2 .AND. (.NOT. LVISI)) GOTO 2001
C 
C ***********************************************************************
C 
C     2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C 
C ***********************************************************************
C 
C     INITIALISATION DES YAC (RECHERCHE COEFF COUPLAGE CIS-COMP A)
C 
      YACMIN  = 0.D0
      YACAMAX = 0.D0
      YACMAX  = 0.D0
C 
C      DEBUT DU TEST SUR LES INTERFACES
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL ADTBDM ('PREDELAINT', PREEND)
C 
        DPREEN = PREEND
C 
C       DEBUT DU TEST SUR LE TYPE
C  
        IF (LTYP01) THEN
C 
C         Ouverture du fichier pour les quantites admissibles en poursuite
C 
          NOM = 'admint'
          NOMFIC = NOM
          LONADM = 6*NPICET
          CALL OFDDNF (3, NOMFIC, 6, LONADM, IUADMC)
          NUENRS = 1
C 
          NOM = 'endomm'
          NOMFIC = NOM
          LONADM = 3*(NPICET+1)
          CALL OFDDNF (3, NOMFIC, 6, LONADM, IUACOU)
          NUENIN = 1
          NUENSU = 1
C 
C         Ouverture du fichier pour lecture des quantites admissibles en reprise
C 
          LONADM = 6*NPICET
	  IF (LTYP1 .AND. (NBNETT .GE. 1)) THEN
	    CALL ADTBM ('LIS-BIGNET', DEBNET)
	    SUPNET = M(DEBNET+NBNETT-1)
	    CALL IDENTI (SUPNET, CARNET)
	    NOM = 'ai_'//CARNET
	    CALL OFDDNF (3, NOM, 6, LONADM, IUADMI)
	    NUENRI = 1
	  END IF
C 
C         Ouverture du fichier pour ecriture des quantites admissibles en reprise
C 
	  IF (LTYP1 .AND. BIGNET) THEN
	    CALL ADTBM ('LIS-BIGNET', DEBNET)
	    SUPNET = M(DEBNET+NBNETT)
	    CALL IDENTI (SUPNET, CARNET)
	    NOM = 'ai_'//CARNET
	    CALL OFDDNF (3, NOM, 6, LONADM, IUADIN)
	    NUENIN = 1
	  END IF
C 
C         Ouverture du fichier pour les contraintes normales
C 
	  LONSCH = 3*NPICET
          NOM = 'sinoch'
          NOMFIC = NOM
          NUENSI  = 1
          CALL OFDDNF (3, NOMFIC, 6, LONSCH, IUNSIG)
C  
        END IF
C 
        IF (LTYP2) THEN
C 
          NOM = 'endomm'
          NOMFIC = NOM
          LONADM = 3*(NPICET+1)
          CALL OFDDNF (3, NOMFIC, 6, LONADM, IUACOU)
          NUENIN = 1
          NUENSU = 1
C 
C         Ouverture du fichier pour les quantites admissibles en poursuite
C 
          NOM = 'admint'
          NOMFIC = NOM
          LONADM = 6*NPICET
          CALL OFDDNF (3, NOMFIC, 6, LONADM, IUADMC)
          NUENRS = 1
C 
C         Ouverture du fichier pour lecture des quantites admissibles en reprise
C 
   	  IF (NBNETT .GE. 1) THEN
   	    CALL ADTBM ('LIS-BIGNET', DEBNET)
   	    SUPNET = M(DEBNET+NBNETT-1)
   	    CALL IDENTI (SUPNET, CARNET)
   	    NOM = 'ai_'//CARNET
   	    CALL OFDDNF (3, NOM, 6, LONADM, IUADMI)
   	    NUENRI = 1
   	  END IF
C 
          IF (LPLASI) THEN
            NOM = 'plasin'
            NOMFIC = NOM
            CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C           Longueur d'un enregistrement
C 
            NUENPL  = 1
            LONPLA  = 2*(NPICET+1)
            CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONPLA, IUNPLA)
          END IF
C 
          IF (LENDIN) THEN
            NOM = 'endint'
            NOMFIC = NOM
            CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C           Longueur d'un enregistrement
C 
            LONEND = 3*(NPICET+1)
            NUENEN  = 1
            CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONEND, IUNEND)
          END IF
C 
          IF (LENDIP) THEN
            NOM = 'enddif'
            NOMFIC = NOM
            CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C           Longueur d'un enregistrement
C 
            LONENP = 3*(NPICET+1)
            NUENEP  = 1
            CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONENP, IUNENP)
          END IF
C 
          IF (LVERRI) THEN
            NOM = 'errint'
            NOMFIC = NOM
            CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C           Longueur d'un enregistrement
C 
            LONERR= NPICET
            NUERR  = 1
            CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONERR, IUNERR)
          END IF
C 
          IF (LDSIGN) THEN
            NOM = 'dsgnor'
            NOMFIC = NOM
            CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C           Longueur d'un enregistrement
C 
            LONDSG= 3*NPICET
            NUDSG = 1
            CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONDSG, IUNDSG)
          END IF
C 
          IF (CRIINT) THEN
            NOM = 'criint'
            NOMFIC = NOM
            CALL DELETF (3, 'q-chapeau', NOMFIC, 6)
C 
C           Longueur d'un enregistrement
C 
            LONCRI= NPICET
            NUCRI = 1
            CALL CFDDNF (3, 'q-chapeau', NOMFIC, 6, LONCRI, IUNCRI)
          END IF
C 
        ENDIF
C 
C     Declaration de tableaux locaux ceux du type n*(npicet+1) corespondent
C     a des valeurs a chaque piquet de temps et il s'agit de faire un decalage
C     de n apres chaque piquet de temps pour recopier la valeur qui vient d'etre
C     calculee comme nouvelle valeur precedente.
C 
C     Creation d'un tableau provisoire pour ranger pour toutes les
C     quantites au piquets de temps
C 
C     SAUAPG QUI CONTIENT : SAUADM ( SGNADM+3) epsilon point admissible X 3*NPICET
C 
C     SIGAPT QUI CONTIENT : SIPADM POUR TOUT (NTETA, YINTEG)            X 3*NPICET
C 
C     QADMPR delta  admissibles de l'etape locale precedente            X 6*NPICET
C     QADMPI copie de QADMPR                                            X 6*NPICET
C 
C     RPRE seuil                                                        X   (NPICET+1)
C     PPRE saut plastique equivalent                                    X   (NPICET+1)
C     SAUPO sauts plastiques dans la base d'orthotropie                 X 2*(NPICET+1)
C 
C     SAUADM saut admissible actuel                                     X 3*(NPICET+1)
C 
C     SGNCAC sign   chapeau actuel                                      X 3*(NPICET+1)
C     SAUEO saut   elastique dans la base d'orthotropie                 X 3*(NPICET+1)
C     ENDOM endommagement d, dp, ds                                     X 3*(NPICET+1)
C     DENDCI endommagement dps, dpt                                     X 3*(NPICET+1)
C     DENDCS endommagement dps, dpt                                     X 3*(NPICET+1)
C     DSGPAC delta sign point chapeau actuel XDT                        X 3*NPICET
C     DSAPAC delta saut point chapeau actuel XDT                        X 3*NPICET
C     ERRVIS denominateur de l'erreur au point de Gauss                 X NPICET
C     DSIGOR sigchap-sigadm orthotropie                                 X 3*NPICET
C     DCRIN CRITERE D'INTERFACE                                         X   NPICET
C                                                                   _________________
C 
C                                                                      51*NPICET+26
        LONPRO = 51*NPICET+26
C 
        CALL POUSMD (LONPRO, DSAUAP)
C 
        DSGNAP = DSAUAP + 3*NPICET
        QADMPR = DSGNAP + 3*NPICET
        QADMPI = QADMPR + 6*NPICET
        DRPRE  = QADMPI + 6*NPICET
        DPPRE  = DRPRE  + NPICET+1
        DAUPOR = DPPRE  + NPICET+1
        DSAADM = DAUPOR + 2*(NPICET+1)
C 
        DSGCAC = DSAADM + 3*(NPICET+1)
        DSAEOR = DSGCAC + 3*(NPICET+1)
        DENDOM = DSAEOR + 3*(NPICET+1)
	DENDCI = DENDOM + 3*(NPICET+1)
	DENDCS = DENDCI + 3*(NPICET+1)
        BSGPAC = DENDCS + 3*(NPICET+1)
        BSAPAC = BSGPAC + 3*NPICET
C 
        DERRVI = BSAPAC + 3*NPICET
        DDSGOR = DERRVI + NPICET
        DDCRII = DDSGOR + 3*NPICET
C 
C       Pour la verification de rsicba et schast
C 
CD      LONPRO = 6*NPICET
CD      CALL POUSMD (LONPRO, SAUVER)
CD      SGNVER = SAUVER+3*NPICET
C 
C -----------------------------------------------------------------------
C 
C       SEQUENCE POUR LA VISU DIRECTE
C 
C -----------------------------------------------------------------------
C 
        IF (LTYP3) THEN
C 
          NOM = 'endomm'
          NOMFIC = NOM
          LONADM = 3*(NPICET+1)
          CALL OFDDNF (3, NOMFIC, 6, LONADM, IUACOU)
          NUENIN = 1
          NUENSU = 1
C 
          CALL POUSME (5, DBCARA)
          DBDEBU = DBCARA+1
C 
C         Pour visualiser les points suivants
C 
1001      CONTINUE
C 
          ADCARA  = DBCARA
          ADEBUT =  DBDEBU
C 
          CALL MESSAO 
     &      ('POINT A VISUALISER DANS LES INTERFACES
     &      \
     &      \ON PRECISE POUR LE POINT 4 DONNEES :
     &      \ 
     &      \  1er    = NUMERO D''INTERFACE
     &      \  2eme   = NUMERO DE COLONNE
     &      \  3eme   = NUMERO DE POINT DE GAUSS EN X
     &      \  4eme   = NUMERO DE L''ANGLE CONCERNE ')
C 
          CALL LECLEN (M(ADEBUT), 4, NLU)
C 
C         S' IL N'Y A PAS DE POINT LU ON SORT
C 
          IF (NLU .EQ. 0) GOTO 2001
C 
          ADEBUT = ADEBUT-1
          M(ADCARA) = (M(ADEBUT+1)-1)*NBCOL*NTETA*XINTEG
     &               +(M(ADEBUT+2)-1)*NTETA*XINTEG
     &               +(M(ADEBUT+3)-1)*NTETA
     &               + M(ADEBUT+4)
C 
          CALL IMPET ('CARACTERISTIQUES DU POINT A VISUALISER ',
     &                 M(DBCARA))
          CALL MESSAO ('NUMERO DES ETAPES GLOBALES A VISUALISER :
     &                 \ON DONNE POUR LE POINT PRECEDENT DANS UN
     &                 \ORDRE CROISSANT LES NUMEROS D''ETAPES
     &                 \GLOBALES POUR ANALYSER LA CONVERGENCE ')
C 
          CALL POUSME (NBETGL+1, ADNUET)
          CALL LECLEN (M(ADNUET), NBETGL+1, NLUET)
C 
C         Nombre d'etapes globale deja considerees
C 
          NBLUET = 0
          NUPDES = 0
          NBPDES = NLUET
          CALL GSPOUD (18*(NPICET+1)*NLUET, TABDES)
          DBTABD =TABDES
C 
          ADNUET = ADNUET-1
2110      CONTINUE
          PREEND =  DPREEN
          ADCARA  = DBCARA
          ADEBUT =  DBDEBU
          TABDES = DBTABD
          INDIC = 0
          ADNUET = ADNUET+1
          CALL IDENTI (M(ADNUET), CARETG)
          CALL INFODP ('FT-EPS-'//CARETG, FTECH6, FTEPMX)
          FTEPMX = FTEPMX/NPICET
          CALL INFODP ('FT-CON-'//CARETG, FTSCH5, FTSIMX)
          FTSIMX = FTSIMX/NPICET
C 
          DBFTEP  = 1
          DBFTSI  = 1
          FIFTEP  = FTEPMX
          FIFTSI  = FTSIMX
C 
          DBCHEP =  1
          DBCHSI =  1
          FICHEP =  FTEPMX
          FICHSI =  FTSIMX
C 
        END IF
C 
C -----------------------------------------------------------------------
C 
C       FIN DE SEQUENCE POUR LA VISU DIRECTE
C 
C -----------------------------------------------------------------------
C 
        INDIC = 0
        ADSAUT = SAUCH4
        ADSIGN = SGNCH3
        DBPOPG = ADPOPG+NGAU1
C 
C       BOUCLE SUR LES INTERFACES
C 
        DO NUINT= 1, NBINT
C 
          CALL ANGINT (NUINT, TETORT)
C 
C         recherche des caracteristiques du comportement de l'interface
C 
          CALL CONLII (NUINT, R0, AL, BE, A1, A2, K, N, YOINT, YC,
     &                 YRINT, ALPINT, MINT, G1, G2, A, K1, K2, K0,
     &                 SORT, SRI13, SRI23)
C 
          TR132  = SRI13*SRI13
	  TR232  = SRI23*SRI23
          DEBGAU =  1+(NUINT-1)*XINTEG*NBCOL
          FINGAU =  NUINT*XINTEG*NBCOL
C 
C         BOUCLE i SUR LES POINTS DE GAUSS DANS L'ELEMENT
C 
          DO PGAU1  = DEBGAU, FINGAU
C 
            MULPPI = DM(DBPOPG)
            DBPOPG = DBPOPG+1
C 
C           BOUCLE ii SUR LES ANGLES
C 
            DO TETA = 1, NTETA
C    
	      ENDCIN  = DENDCI
	      ENDCSU  = DENDCS
              INDIC = INDIC+1
              NUGAU = MOD(PGAU1-DEBGAU+1, XINTEG)
              IF (NUGAU .EQ. 0) NUGAU = XINTEG
              NUCOL  = (PGAU1-DEBGAU-NUGAU+1)/XINTEG+1
	      LONADM = 3*(NPICET+1)
	      NUENIN = (NUINT-1)*NBCOL*XINTEG*NTETA
     &                +(NUCOL-1)*XINTEG*NTETA
     &                +(NUGAU-1)*NTETA+TETA
	      CALL LFDDNF (DM(ENDCIN), ENDCIN, LONADM, IUACOU, NUENIN)
	      NUENSU = (NUINT)*NBCOL*XINTEG*NTETA
     &                +(NUCOL-1)*XINTEG*NTETA
     &                +(NUGAU-1)*NTETA+TETA
	      CALL LFDDNF (DM(ENDCSU), ENDCSU, LONADM, IUACOU, NUENSU)
C 
              IF (LTYP3) THEN
C 
                LVISU  = .FALSE.
                LVISV  = .FALSE.
                LPASS3 = .FALSE.
C 
C               DEBUT DE TEST POUR SAVOIR SI ON CONSERVE LE POINT DE GAUSS COURANT
C 
                IF (INDIC .EQ. M(DBCARA)) THEN
C 
C                 Le point est a visualiser
C 
                  LVISU = .TRUE.
                  NUPDES = NUPDES+1
                  NBLUET = NBLUET+1
                  IF (NBLUET .EQ. NLUET)  LVISV = .TRUE.
C 
                  CALL IDENTI (NUINT, CARCOU)
                  CALL IDENTI (NUCOL, CARCOL)
                  CALL IDENTI (NUGAU, CARGAX)
                  CALL IDENTI (TETA, CARTET)
                  CALL IDENTI (NBETLC, CARITE)
C 
                  NOMSAU(1) =
     &            '                      '//CARITE//' , '//CARCOU//
     &            ' , '//CARCOL//' , '//CARGAX//' , '//CARTET
C 
                END IF
C 
C               FIN DE TEST POUR SAVOIR SI ON CONSERVE LE POINT DE GAUSS COURANT
C  
                IF (LTYP3 .AND. .NOT. LVISU) THEN
                  ADSAUT = ADSAUT + CHARAX*NSAU
                  ADSIGN = ADSIGN + CHARAX*NSAU
                  PREEND = PREEND + 3
                  GOTO 1101
                END IF
C 
              END IF
C 
C             ATTENTION TOUTE LES QUANTITES SONT MULTIPLIES PAR INTERV
C             On calcule la valeur des champs admissibles pour tous les piquets de
C             temps : debut caracterise la position dans le tableau des contraintes
C             normales ou des sauts admissibles de la 1ere valeur interessante
C 
              SAUAPG = DSAUAP
              SGNAPG = DSGNAP
C 
              CALL ADCTPS (FTEPMX, FTSIMX, NSAU,
     &                     DM(ADSAUT), DBCHEP, FICHEP,
     &                     DM(ADSIGN), DBCHSI, FICHSI,
     &                     DM(FTECH6), DBFTEP, FIFTEP,
     &                     DM(FTSCH5), DBFTSI, FIFTSI,
C                          Et on recupere:
     &                     DM(SAUAPG), DM(SGNAPG))

C 
C             Pour aller lire le bon point de Gauss dans les tableaux
C             de stockage des deformations et des contraintes.
C 
	      ADSAUT = ADSAUT + CHARAX*NSAU
              ADSIGN = ADSIGN + CHARAX*NSAU
C 
              LONADM =6*NPICET
	      IF (LTYP0) THEN
                CALL LFDDNF (DM(QADMPR), QADMPR, LONADM, IUADMC, NUENRS)
              END IF
              IF (LTYP1 .AND. (NBNETT .EQ. 0)) THEN
		CALL MENADM (QADMPR, LONADM)
              END IF
              IF (LTYP1 .AND. (NBNETT .GE. 1)) THEN
                CALL LFDDNF (DM(QADMPR), QADMPR, LONADM, IUADMI, NUENRS)
              END IF
              IF (LTYP2 .OR. LTYP3) THEN
                CALL MENADM (QADMPR, 6*NPICET)
              END IF
              IF (LTYP2 .AND. (NBNETT .GE. 1)) THEN
                CALL LFDDNF (DM(QADMPR), QADMPR, LONADM, IUADMI, NUENRS)
                NUENRS = NUENRS + 1
              END IF
C 
	      DQADMP = QADMPR
              DQADMI = QADMPI
C 
C             Recherche de l'angle de la bande (tetcal) correspondant a teta
C 
              TETCAL = DM(ADTETA+TETA-1)
              ANGORT = TETCAL-TETORT
C 
C -----------------------------------------------------------------------
C 
C             Debut d'initialisation des quantites chapeaux
C 
C -----------------------------------------------------------------------
C 
C             initialisation des champs avec decalage :
C             SIGAPG, SIPADM sigma point admissible
C             DSAPAC delta saut point chapeau actuel XDT
C             DSPCAC - delta sigma point chapeau actuel XDT
C 
              RPRE   = DRPRE
              PPRE   = DPPRE
              SAUPOR = DAUPOR
              SAUADM = DSAADM
C 
              SGNCAC  = DSGCAC
              SAUEOR  = DSAEOR
              ENDOM   = DENDOM
              DSGPAC  = BSGPAC
              DSAPAC  = BSAPAC
C 
              ERRVIS = DERRVI
              DSIGOR = DDSGOR
              DCRIIN = DDCRII
C 
              DM(RPRE) = R0
              DM(PPRE) = 0.D0
C 
              CALL MENADM (SAUPOR, 2)
              CALL MENADM (SAUEOR, 3)
              CALL MENADM (SAUADM, 3)
C 
              IF (NUPATE .EQ. 0) THEN      
                DM(ENDOM)   = DM(PREEND)
                PREEND = PREEND +1
                DM(ENDOM+1) = DM(PREEND)
                PREEND = PREEND +1
                DM(ENDOM+2) = DM(PREEND)
                PREEND = PREEND +1
              ELSE
	        CALL MENADM (DENDOM, LONEND)
	      END IF
C 
              CALL MENADM (SGNCAC, 3)
C 
C -----------------------------------------------------------------------
C 
C             Fin d'initialisation des quantites chapeaux
C 
C -----------------------------------------------------------------------
C 
C             Sauvegarde pour la visu des quantites initiales
C             SAUVII est compilee a vide dans delami, mais appelle des routines
C             graphiques dans visu_dsdm (---> repertoire faux_commun_visu)
C 
              IF (LVISU .AND. LTYP3) THEN
                CALL SAUVII (0, ZERO, ZERO, ZERO, ZERO, DM(ENDOM),
     &                       ZERO, ZERO, TABDES, NUPDES)
              END IF
C 
C               Indicateurs pour savoir si les endommagements precedents valent 1
C 
              NOCVD  = .TRUE.
              NOCVD1 = .TRUE.
              NOCVD2 = .TRUE.
C 
              CALL BALAID (3, SINOP)
              CALL BALAID (3, SGNAOR)
C 
C               DEBUT DE BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
              DO  TEMPS = 1, NPICET
C 
                IF (NUPATE .EQ. TEMPS) THEN      
                  DM(ENDOM)   = DM(PREEND)
                  PREEND = PREEND +1
                  DM(ENDOM+1) = DM(PREEND)
                  PREEND = PREEND +1
                  DM(ENDOM+2) = DM(PREEND)
                  PREEND = PREEND +1
		END IF

                INTERV = DM(ADINTE+TEMPS-1)
C 
                IF (DM(ENDOM) .EQ. 1.D0) THEN
                  DM(ENDOM+3) = DM(ENDOM)
		  GOTO 1111
                END IF
C 
	        IF (BIGNET) CALL COPITD (6, DQADMP, DQADMI)
C 
C               Integration a saut impose de la plasticite de l'interface
C 
		CALL INPLAI
C                      Pour la base d'orthotropie
     &                (ANGORT, R0, AL, BE, A1, A2, K1, K2,
     &                 DM(RPRE), DM(PPRE), DM(RPRE+1), DM(PPRE+1),
     &                 DM(SAUPOR), DM(SAUPOR+2),
     &                 DM(SAUEOR), DM(SAUEOR+3),
C                      Pour la base locale
     &                 DM(DQADMP), DM(SAUAPG), DM(SGNAPG),
     &                 DM(SAUADM), DM(SAUADM+3))
C 
	        SAUEOR = SAUEOR +3
                SAUPOR = SAUPOR +2
C 
C               Modele original LMT 
C 
                INVALP = 1.D0/ALPINT
                YD1 = 0.D0
                YD2 = 0.D0
                YD3 = 0.5D0*K0*DMAX1(DM(SAUEOR+2), 0.D0)
     &                        *DMAX1(DM(SAUEOR+2), 0.D0)
                YAC = (((G1*YD1)**ALPINT) + ((G2*YD2)**ALPINT)
     &                + (YD3 **ALPINT))**INVALP
C 
C              DR Ajout le 17/1/96 : On tient compte de la pression hydrostatique
C                                     au niveau de l'interface
C                 A : Parametre materiau tel que  A = Ecisail / Ecompress en elasticite
C 
                YACA    = K0*DMIN1(DM(SAUEOR+2), 0.D0)
     &                      *DMIN1(DM(SAUEOR+2), 0.D0)/2.D0
                YAC     = YAC - (A*YACA)
                YAC     = DMAX1 (YAC, 0.D0)
                YACAMAX = DMAX1 (YACA, YACAMAX)
   		YACMAX  = DMAX1 (YAC, YACMAX)
                IF (YACAMAX .EQ. YACA) YACMIN = YAC
C 
                CALL VAENDI (K, N, YOINT, YC, YRINT, YD3, MINT,
     &                       INTERV, NOCVD, YAC, DM(ENDOM),
     &                       DM(ENDOM+3))
C 
C               Calcul de la contrainte chapeau dans la base d'orthotropie
C 
                SINOR(1)  =  K1*(1.D0-DM(ENDOM+3))*DM(SAUEOR)
                SINOR(2)  =  K2*(1.D0-DM(ENDOM+4))*DM(SAUEOR+1)
C 
C               On introduit ici le caractere unilateral en 33
C 
                SINOR(3)  =
     &             K0 * (1.D0-DM(ENDOM+5)) * DMAX1 (DM(SAUEOR+2), 0.D0)
     &            +K0 * DMIN1 (DM(SAUEOR+2), 0.D0)
C 
C               Modele pheno. ecrit en endommagement des plis adjacents ;
C               reproduction du schéma en double hélice, zone de traction 
C               interfaciale, bande de compression interfaciale...
C 
   		DPINF = DMAX1(DM(ENDCIN+4), DM(ENDCIN+5))
   		DPSUP = DMAX1(DM(ENDCSU+4), DM(ENDCSU+5))
   		DM(ENDOM+3) = DMAX1(0.D0, DPINF-DPSUP)
   		DM(ENDOM+3) = DMAX1(DM(ENDOM), DM(ENDOM+3))
1111            CONTINUE
                DM(ENDOM+4) = DM(ENDOM+3)
                DM(ENDOM+5) = DM(ENDOM+3)
C 
C               Calcul de la contrainte chapeau dans la base d'orthotropie
C 
C               SINOR(1)  =  K1*(1.D0-DM(ENDOM+3))*DM(SAUEOR)
C               SINOR(2)  =  K2*(1.D0-DM(ENDOM+4))*DM(SAUEOR+1)
C 
C               On introduit ici le caractere unilateral en 33
C 
C               SINOR(3)  =
C    &             K0 * (1.D0-DM(ENDOM+5)) * DMAX1(DM(SAUEOR+2), 0.D0)
C    &            +K0 * DMIN1(DM(SAUEOR+2), 0.D0)
C 
C               CALCUL 'VALABLE' A L'ETAPE LOCALE 1 (SOLUTION ELASTIQUE)
C 
                CALL QIBORT (ANGORT, DM(SAUADM+3), SAUIMP)
C 		  
                DM(DCRIIN)=(((K1*SAUIMP(1))*(K1*SAUIMP(1)))/TR132)
     & 		          +(((K2*SAUIMP(2))*(K2*SAUIMP(2)))/TR232)   
C       
C               Calcul de  sgn_point_chapeau dans la base d'orthotropie
C 
                SINPOR(1) = SINOR(1) -SINOP(1)
                SINPOR(2) = SINOR(2) -SINOP(2)
                SINPOR(3) = SINOR(3) -SINOP(3)
C 
C               Decalage d'indice des tableaux pour aller lire au bon endroit
C 
                DQADMP   = DQADMP + 6
                DQADMI   = DQADMI + 6
                RPRE     = RPRE   + 1
                PPRE     = PPRE   + 1
                SAUADM   = SAUADM + 3
                SGNCAC   = SGNCAC + 3
                ENDOM    = ENDOM  + 3
		ENDCIN   = ENDCIN + 3
		ENDCSU   = ENDCSU + 3
C 
C               Calcul de sgn_point_n dans la base d'orthotropie
C 
		CALL QIBORT (ANGORT, DM(SGNAPG), SGAPOR)
C 
C               Calcul de sigma_adm dans la base d'orthotropie
C 
                CALL ADDMAD (3, SGNAOR, SGAPOR, SGNAOR)
C 
C               Calcul de sig-cha -sigad dans la base d'orthotropie
C 
                CALL SOUMAD (3, SINOR, SGNAOR, DM(DSIGOR))
C 
C               Calcul de la contribution du point au numerateur de l'erreur
C 
                CALL PROIOR (SORT, DM(DSIGOR), DM(DSIGOR), NUMERR)
C 
                IF (NUMERR .LT. 0.D0) THEN
                  CALL IMPET  ('NUM INTERFACE     ', NUINT)
                  CALL IMPET  ('NUM COLONNE       ', NUCOL)
                  CALL IMPET  ('NUM X             ', X)
                  CALL IMPET  ('NUM ANGLE         ', TETA)
                  CALL IMPET  ('PAS DE TEMPS      ', TEMPS)
                  CALL IMPTDT ('K-1 ORTHOTROPIE   ', SORT, 1, 3)
                  CALL IMPTDT ('DSIGOR            ', DM(DSIGOR), 1, 3)
                  CALL IMPDT  ('K-I*DSIGOR**2=NUMERR ', NUMERR)
                END IF
C 
                DM(ERRVIS)       = MULPPI*NUMERR
                DM(NERRLC+TEMPS) = DM(NERRLC+TEMPS) + DM(ERRVIS)
C 
C               Calcul de sig-cha +sigad dans la base d'orthotropie
C 
                CALL ADDMAD (3, SGNAOR, SINOR, SINSOM)
C 
C               Calcul de la contribution du point au denominateur de l'erreur
C 
                CALL PROIOR (SORT, SINSOM, SINSOM, DENOMI)
                DM(DERRLC+TEMPS) = DM(DERRLC+TEMPS) + MULPPI*DENOMI
                DM(DINTLO+TEMPS) = DM(DINTLO+TEMPS) + MULPPI*DENOMI
C 
C               Calcul de dspcac = ( sigma point chapeau ) - (sigma point adm )
C               dans la base d'orthotropie
C 
                CALL SOUMAD (3, SINPOR, SGAPOR, DSGPOR)
                CALL MCORIN (SORT, DSGPOR, DSACAO)
C 
C               Calcul de depcac = K0-1 (( sigma chapeau ) - (sigma adm ))
C               dans la base locale
C 
                CALL QIBORT (-ANGORT, DSACAO, DM(DSAPAC))
C 
C               Calcul de -dspcac = - (( sigma chapeau ) - (sigma adm ))
C               dans la base locale
C 
                CALL QIBORT (-ANGORT, DSGPOR, DM(DSGPAC))
                CALL MUMARE (-1.D0, 3, DM(DSGPAC), DM(DSGPAC))
C 
                IF (LVISU .AND. LTYP3) THEN
C 
C                 On remplit le tableau DM a partir de advisu
C                 Calcul de saut-chapeau dans la base d'orthotropie
C 
                  CALL ADDMAD (2, DM(SAUEOR), DM(SAUPOR), SAUCOR)
                  SAUCOR(3) = DM(SAUEOR+2)
C 
                  IF (LOGR0) THEN
                    CALL QIBORT (-ANGORT, SAUCOR, EPCLO)
                    CALL QIBORT (-ANGORT, SGNAOR, SIGALO)
                    CALL QIBORT (-ANGORT, SINOR, SIGCLO)
                    CALL SAUVII (TEMPS, EPCLO, SIGALO,
     &                           EPCLO, SIGCLO, DM(ENDOM), DM(SAUPOR),
     &                           DM(ERRVIS), TABDES, NUPDES)
                  ELSE
		    CALL SAUVII (TEMPS, SAUCOR, SGNAOR,
     &                           SAUCOR, SINOR, DM(ENDOM), DM(SAUPOR),
     &                           DM(ERRVIS), TABDES, NUPDES)
                  END IF
C 
                END IF
C 
C                 Decalage d'indice des tableaux pour aller lire au bon endroit
C 
                DSAPAC   = DSAPAC + 3
                DSGPAC   = DSGPAC + 3
                SAUAPG   = SAUAPG + 3
                SGNAPG   = SGNAPG + 3
                DSIGOR   = DSIGOR + 3
                ERRVIS   = ERRVIS + 1
                DCRIIN   = DCRIIN + 1
C 
                SINOP(1) = SINOR(1)
                SINOP(2) = SINOR(2)
                SINOP(3) = SINOR(3)
C 
C               FIN DE BOUCLE iii SUR LES PIQUETS DE TEMPS
C 
              END DO
C 
C             CHOIXD est compilee a vide dans delami, mais appelle des routines
C             graphiques dans visu_dsdm (---> repertoire faux_commun_visu)
C 
              IF (LVISV .AND. LTYP3) THEN
                CALL CHOIXD (4, NBPDES, DM(TABDES), DM(DEVALT))
              END IF
C 
              IF (LTYP2) THEN
C 
C               Ecriture des sauts plastiques
C 
                IF (LPLASI) THEN
                  CALL EFDDNF (IUNPLA, DM(DAUPOR),
     &                         DAUPOR, LONPLA, NUENPL)
                  NUENPL = NUENPL +1
                END IF
C 
C               Ecriture des endommagements
C 
                IF (LENDIN) THEN
                  CALL EFDDNF (IUNEND, DM(DENDOM),
     &                         DENDOM, LONEND, NUENEN)
                  NUENEN = NUENEN + 1
                END IF
C 
                IF (LENDIP) THEN
                  CALL EFDDNF (IUNENP, DM(DENDOM),
     &                         DENDOM, LONENP, NUENEP)
                  NUENEP = NUENEP + 1
                END IF
C 
              END IF
C 
C             Pour avoir directement les contraintes permettant
C             le calcul des seconds membres pour l'etape globale
C 
C               <=> on stocke -delta(sgn-point-chapeau)
C               <=> on stocke  delta(sau-point-chapeau)
C 
              IF (LTYP01) THEN
C 
C               Ecriture de la solution admissible de l'etape locale actuelle.
C 
		CALL EFDDNF (IUNSIG, DM(BSGPAC), BSGPAC, LONSCH, NUENSI)
                IF (LTYP0) THEN
                  CALL EFDDNF (IUADMC, DM(QADMPR),
     &                         QADMPR, LONADM, NUENRS)
                END IF
                IF (LTYP1 .AND. BIGNET) THEN
                  CALL EFDDNF (IUADIN, DM(QADMPR),
     &                         QADMPR, LONADM, NUENRS)
                  CALL EFDDNF (IUADMC, DM(QADMPR),
     &                         QADMPR, LONADM, NUENRS)
                END IF
                IF (LTYP1 .AND. (.NOT. BIGNET)) THEN
                  CALL EFDDNF (IUADMC, DM(QADMPR),
     &                         QADMPR, LONADM, NUENRS)
                END IF
		NUENSI = NUENSI +1
                NUENRS = NUENRS +1
              END IF
C 
              IF (LTYP2) THEN
C 
                IF (LVERRI) THEN
                  CALL EFDDNF (IUNERR, DM(DERRVI),
     &                         DERRVI, LONERR, NUERR)
                  NUERR = NUERR+1
                END IF
C 
                IF (LDSIGN) THEN
                  CALL EFDDNF (IUNDSG, DM(DDSGOR),
     &                         DDSGOR, LONDSG, NUDSG)
                  NUDSG = NUDSG +1
                END IF
C 
	        IF (CRIINT) THEN
                  CALL EFDDNF (IUNCRI, DM(DDCRII),
     &                         DDCRII, LONCRI, NUCRI)
                  NUCRI = NUCRI +1
                END IF
C 
              END IF
C 
C             SORTIE POUR TYPE = 3 SI LE POINT N'EST PAS A VISUALISER
C 
1101          CONTINUE
C 
C           FIN DE BOUCLE ii SUR LES ANGLES
C 
            END DO
C 
C         FIN DE BOUCLE i SUR LES POINTS DE GAUSS DANS L'INTERFACE
C 
          END DO
C 
C       FIN DE BOUCLE SUR LES INTERFACES
C 
        END DO
C 
	IF (LENDIN) CALL FERFIC (3, IUNEND, IDPROG)
	IF (LENDIP) CALL FERFIC (3, IUNENP, IDPROG)
        IF (LPLASI) CALL FERFIC (3, IUNPLA, IDPROG)
        IF (LVERRI) CALL FERFIC (3, IUNERR, IDPROG)
        IF (LDSIGN) CALL FERFIC (3, IUNDSG, IDPROG)
        IF (CRIINT) CALL FERFIC (3, IUNCRI, IDPROG)
C 
C       Pour l'etape suivante en visu directe
C       Pour un autre point de Gauss en visu directe
C       Sinon, fin de la visu des points de couche.
C 
        IF (LTYP3) THEN
          IF (NBLUET .NE. NLUET) THEN
	    GOTO 2110
	  ELSE
	    GOTO 1001
	  END IF
        END IF
C 
C     FIN DU TEST SUR LES INTERFACES
C 
      END IF
C 
      IF (LTYP01) THEN
        CALL FERFIC (3, IUADMC, IDPROG)
        CALL FERFIC (3, IUNSIG, IDPROG)
        CALL FERFIC (3, IUACOU, IDPROG)
      END IF
      IF (LTYP3) THEN
        CALL FERFIC (3, IUACOU, IDPROG)
      END IF
      IF (LTYP1 .AND. (NBNETT .GE. 1)) THEN
	CALL FERFIC (3, IUADMI, IDPROG)
      END IF
      IF (LTYP1 .AND. BIGNET) THEN
	CALL FERFIC (3, IUADIN, IDPROG)
      END IF
      IF (LTYP2) THEN
	CALL FERFIC (3, IUADMC, IDPROG)
      END IF
      IF (LTYP2 .AND. (NBNETT .GE. 1)) THEN
	CALL FERFIC (3, IUADMI, IDPROG)
      END IF
C 
C **********************************************************************
C *
C *   FIN DE 2EME BOUCLE SUR L'ESPACE <===> ELEMENT DE TYPE INTERFACE
C *
C **********************************************************************
C 
C     CALCUL DE L'ERREUR
C 
      CALL IMPTDT ('DENOMI AV DE L''ERREUR ',
     &              DM(DERRLC+1), 1, NPICET)
      CALL intgft (1, ADINTE, DM(DERRLC+1), DM(DERRLC+1))
      CALL IMPTDT ('DENOMI AP DE L''ERREUR ',
     &              DM(DERRLC+1), 1, NPICET)
      CALL IMPTDT ('NUMERR AV DE L''ERREUR ', DM(NERRLC+1), 1, NPICET)
      CALL intgft (1, ADINTE, DM(NERRLC+1), DM(NERRLC+1))
      CALL IMPTDT ('NUMERR AP DE L''ERREUR ', DM(NERRLC+1), 1, NPICET)
C 
      DO TEMPS  = 1, NPICET
        DM(ERRTOT+TEMPS) = DMAX1(DM(ERRTOT+TEMPS-1),
     &                     DSQRT(DM(NERRLC+TEMPS)/DM(DERRLC+TEMPS)))
      END DO
C 
      IF (ERTOTA .AND. LTYP2) THEN
        CALL DESERR (NPICET, DM(DEVALT), DM(NERRLC), DM(DERRLC),
     &               DM(ERRTOT))
      END IF
C 
C     Determination du prochain nombre de piquets par comparaison de
C     l'erreur actuelle et de l'ancienne erreur par une valeur maxi
C     de l'erreur actuellement .5
C 
      NPICAC = NPICET
C 
      CALL IMPDT  ('DUREE TOTALE DE L''HISTOIRE ', DUREE)
      CALL IMPTDT ('ERREUR GLOBALE ', DM(ERRTOT), 1, NPICET+1)
C 
C     Incrementation du nombre d'etapes locales effectuees
C 
      CALL ECERGL (DM(NERRLC), DM(DERRLC), DM(ERRTOT))
C 
      CALL IMPET ('NOMBRE PARTIEL D''ETAPES LOCALES ', NBETLC)
C 
      NBETLC = NBETLC+1
C 
C     T1 = TUSED()
      DT = (T1-T0)/50
      DTH= DT/3600
      DTM= (DT-3600*DTH)/60
      DTS= (DT-3600*DTH-60*DTM)
      CALL MESSAO ('TEMPS UTILISE EN TOUT DANS '//IDPROG)
      CALL IMPET  ('TEMPS UTILISE EN HEURES   : ', DTH)
      CALL IMPET  ('TEMPS UTILISE EN MINUTES  : ', DTM)
      CALL IMPET  ('TEMPS UTILISE EN SECONDES : ', DTS)
C 
C     Verification de rsicba et schast
C 
      CALL VERISR (LTYP01)
C 
2001  CONTINUE
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
