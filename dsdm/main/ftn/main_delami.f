!> Bonjour je suis DSDM et je calcule des delaminages
C
!> @author
!> -# Olivier Allix LMT-Cachan
!> -# Stephane Guinard EADS IW
!
! DESCRIPTION: 
C ----------------------------------------------------------------------
C  
C>     QUE FAIT CETTE SUBROUTINE :
C> 
C>           Dans ce programme principal on initialise interpret.
C>           Les resultats errones seront dans le fichier 'erreur-lecture'.
C>           Les resultats qui sont imprimes avec l'aide des routines
C>           imp... seront imprimes sur le fichier 'delami-resu', si la
C>           trace est activee dans les routines correspondantes.
C> 
C>           Elles se terminent par :
C> 
C>                 - T si on veut toujours imprimer
C>                 - N si on veut imprimer en trace normale
C>                 - P si on veut imprimer en trace profonde
C> 
C>           Elle appelle ensuite les routines :
C> 
C>          CALCUL pour rentrer les donnees utilisateur et dans CALCUL
C>                 - inigeo (initialisation du common geometrie)
C>                   apres que les VALeurs aient ete affectees.
C>                 - si apres cette sequence d'initialisation on  veut
C>                   sortir on tape le mot cle exit
C>                 - si on desire executer on tape le mot cle
C>                   execution
C 
C 
C ----------------------------------------------------------------------
C     CETTE ROUTINE APPELLE:
C 
C                    -CALCUL
C                    -CATELE
C 
C ----------------------------------------------------------------------
C 
      PROGRAM MAIN_DELAMI
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
      INTEGER         ITAB(4),I, NUECR, NUFICH, LECINT, IUECRA
      LOGICAL         LECLOG, STAT, RECURS
      LOGICAL         OK
      CHARACTER*40    MOT
      CHARACTER*30    FICRESU
C 
      INTEGER         LGCARG, IUCLAV
      INTEGER         LONGGE
C 
      INTEGER         TUSED, DT, DTH, DTS, T1, T0
C 
      DATA (ITAB(I), I=1, 4) / 6, 6, 5, 2/
C 
C -----------------------------------------------------------------------
C     Initialisation de l'interpreteur de mots-clefs
C     (cf /home/sg/DDD/DSDM/interpret/mac_dep/ftn/insint.f)
C 
      CALL INSINT 
C 
C     Definition du prompt
C 
      CALL SPROMP ('GIVE')
C 
C     Initialisation de la debugg-tools-lib : affecte un numero d'unite
C     aux differents fichiers accessibles par les routines de la
C     "debugg-tools-lib"; permet de demander un controle de recursivite;
C     permet de demander un chronometrage et l'edition de statistiques.
C     (cf ~lmtutils/interpret/devel/ftn/a_debug_trac_lib.f)
C 
      CALL IDEBTO (ITAB, 3, .TRUE., .FALSE., 'ERREUR-LECTURE')
C 
C     Definition du fichier de traces : lit le nom du fichier de definition
C     des modules tracables, et des blocs-de-modules.
C     (cf ~lmtutils/interpret/devel/ftn/a_debug_trac_lib.f)
C 
      CALL LECTRA ('lhes')
C 
C     Definition du fichier contenant les mots-cles :
C     lit et transcrit les informations contenues dans le fichier help
C     (cf ~lmtutils/interpret/devel/ftn/lechel.f)
C 
      CALL LECHEL ('/home/sg/dsdm/delami/lhdo')
C 
C     On met tild a blanc. Attention il faut que la longueur du tild
C     ne depasse pas 120.
C 
      DO I = 1, 120
C 
        tild(i:i)=' '
C 
      END DO
C 
C     Pour une utilisation non-standard avec tous les commentaires internes
C 
      CALL LECSEQ ('UTIL', 'TYPE D''UTILISATION')
C 
      STANDA  = LECLOG ('STANDARD')
C
C     SGU le 04/08/2008
C     Appel a SETDIR pour une relecture en reprise de rep/rep_dir
C     Sans redonner la main à l'utilisateur
C
C      CALL SETDIR
C      STANDA  = .FALSE.
C 
      IF (.NOT. STANDA ) THEN
C 
        CALL LECSEQ ('INITIALISATION-INTERPRETEUR', ' ')
C 
        ELSE
C 
C       Traite la carte contenue dans la chaine de caractere CHAINE.
C       Imagine les formats elementaires utiles a la relecture,
C       stocke dans COMMON/FORCAR/ et COMMON/FORBLA/. Relit la carte du
C       COMMON/CARTE/, suivant le format imagine par GETFOS. Le mot clef
C       est renvoye complet dans la variable MOT.
C       Ignore les cartes commentaire commencant par '**'.
C       (cf ~lmtutils/interpret/devel/ftn/a_interpret_1.f)
C 
        CALL LFORCS ('INIT-INTER,,,,,,', 16, MOT, OK)
C 
      ENDIF
C 
C     Lit la valeur d'un argument de type CHARACTER
C     (cf ~lmtutils/interpret/devel/ftn/a_interpret_5.f)
C 
      CALL LECCHA ('tild', tild)
      LTILD = LGCARG (tild, 120)
C 
C     Attention il faut que la longueur du tild ne depasse pas 120.
C     On met tild1 a blanc.
C 
      DO I = 1, 120
C 
      tild1(i:i)=' '
C 
      END DO
C 
      CALL LECCHA ('tild1', tild1)
      LTILD1 = LGCARG (tild1, 120)
C 
      CALL LECCHA ('FICRESU', FICRESU)
      NUFICH  = LECINT ('NUFICH')
      RECURS = LECLOG ('RECURSIVITE')
      STAT   = LECLOG ('STATISTIQUE-CHRONO')
C 
      NUECR   = IUECRA ()
      ITAB(1) = NUECR
C 
      IF (NUFICH .NE. NUECR) THEN
C 
C       Ouverture du fichier
C 
        OPEN (UNIT=NUFICH, FILE=tild(1:LTILD)// '/'// FICRESU)
C 
      END IF
C 
C     On reaffecte le numero d'unite pour le fichier d'impression.
C 
      ITAB(2) = NUFICH
C 
C     On recupere le numero d'unite qui avait ete affecte pour
C     le fichier reprise. (cf ~lmtutils/interpret/devel/ftn/a_debug_trac_lib.f)
C 
      ITAB(3) = IUCLAV ()
C 
      CALL IDEBTO (ITAB, 3, RECURS, STAT, 'ERREUR_LECTURE.ERR')
C 
      CALL INFOME (LONGGE)
      LDMEFF = MIN0 (LONGGE, LDM)
      CALL IMPET ('LONGUEUR RESERVEE          : '//'MAIN', LDM)
      CALL IMPET ('LONGUEUR PAR INFOME        : '//'MAIN', LONGGE)
      CALL IMPET ('LONGUEUR DE TRAVAIL LDMEFF : '//'MAIN', LDMEFF)
C 
C     Comptage du temps user (cf ~lmtutils/divers/ftn/mac_dep.f)
C 
      T0 = TUSED ()
C 
      CALL PROPRI
C 
      T1 = TUSED ()
      DT = (T1-T0)/50
      DTH= DT/3600
      DTS= DT-3600*DTH
      CALL IMPET ('TEMPS UTILISE, HEURES   : ', DTH)
      CALL IMPET ('TEMPS UTILISE, SECONDES : ', DTS)
C 
C     Edition des statistiques d'appel
C 
      CALL STAWLD 
C 
C     Relache le terminal graphique
C 
      CALL XXIT 
      END
