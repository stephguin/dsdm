C -
C -   Nombre de fonction du temps en deformation pour l'etape
C -   globale consideree
C -
      INTEGER FTEPMX
C -
C -   Nombre de fonction du temps en contrainte pour l'etape
C -   globale consideree
C -
      INTEGER FTSIMX
C -   
C -   Numero de debut de fonction du temps en deformation pour l'etape
C -   globale consideree
C -
      INTEGER DBFTEP
C -   
C -   Numero de debut de fonction du temps en contrainte pour l'etape
C -   globale consideree
C -
      INTEGER DBFTSI
C -   
C -   Numero de debut de fonction de l'espace en deformation pour l'etape
C -   globale consideree
C -
      INTEGER DBCHEP
C -   
C -   Numero de debut de fonction de l'espace en contrainte pour l'etape
C -   globale consideree
C -
      INTEGER DBCHSI
C -   
C -   Numero de fin de fonction du temps en deformation pour l'etape
C -   globale consideree
C -
      INTEGER FIFTEP
C -   
C -   Numero de fin de fonction du temps en contrainte pour l'etape
C -   globale consideree
C -
      INTEGER FIFTSI
C -   
C -   Numero de fin de fonction de l'espace en deformation pour l'etape
C -   globale consideree
C -
      INTEGER FICHEP
C -   
C -   Numero de fin de fonction de l'espace en contrainte pour l'etape
C -   globale consideree
C -
      INTEGER FICHSI
C -
C -   Numero d'etape globale demandee
C -
      INTEGER NUETGL
C -
C   si on veut tracer en postscript
C -
       LOGICAL LPOSCR
C -
C   si on veut tracer en postscript couleur ou non
C -
       LOGICAL LPSCOL
C -
C   commentaires en anglais
C -
       LOGICAL ANGLAI 
C -
C   pour visualiser la zone bord au sein d'une plaque rectangulaire
C   vises pour la plaque rectangulaire
C - visest autrement
C -
       LOGICAL VISES , VISEST
C -          
C -
C - de coins
C -     
       DOUBLE PRECISION XMINPL, XMAXPL , YMINPL , YMAXPL
C -
C -    POUR IMPOSER LES SEUILS A 0 ET 1  DANS LECHAM
C -    SI ON VOIT DE L'ENDOMMAGEMENT
C -
       LOGICAL VIENDO             
C -
C -    Pour normaliser les sorties plasto-plan
C -    pour comparer les differentes versions du modele
C -
       LOGICAL VIPLAS
C -
C= en r                   NORMAR
C= en teta                NORMAT
C= en z                   NORMAZ
C= en contrainte          NORMAC
C= en deformation         NORMAD
C= en deplacement         NORMAU
C -
      DOUBLE PRECISION NORMAR , NORMAT ,  NORMAZ , NORMAC
      DOUBLE PRECISION NORMAD , NORMAU ,  NORABS , NORORD
C -
C - de longueur reel de dessin
C -     
      REAL*4 LONGDE , LARGDE
C -
C -   Pour creer par etloca lorsque TYPE =2 les fichiers
C -   de visu pour les couchea
C -   Pour visualiser les champs chapeau dans les couches
      LOGICAL LVISC
      LOGICAL LPLAST , LENDOM , LMODUL , LVERRE , LDSIGO
C - 
C -   Pour visualiser les champs chapeau dans les toutes les 
C -   couches à tous les pas de temps
C -
      LOGICAL LTOTAL
C - 
C -   Pour creer par etloca lorsque TYPE = 2 les fichiers
C -   de visu pour les interfaces
C -   Pour visualiser les champs chapeau dans les interfaces
C -   Pour calculer la surface totale delaminee (d>0.95),
C -   la surface totale endommagee (d>0.05)
      LOGICAL LVISI
      LOGICAL LPLASI, LENDIN, LENDIP, LSUDEL, LVERRI, LDSIGN
C -
C -   Pour creer par etloca lorsque TYPE =2 les fichiers
C -   de visu pour l'erreur totale
      LOGICAL ERTOTA
C -   Pour visualiser les criteres elastiques dans les interfaces
      LOGICAL CRIINT
C 
C     Pour visualiser les champs chapeau dans les couches
       LOGICAL CRICOU  
c -
c -   Pour visualiser les differences en contrainte chapeau admissible
C -
      LOGICAL VDIFFE
c - 
c -   Pour indiquer si le fichiers de visaulisation ont deja ete construits
c -
      LOGICAL DEJPAS
C -
C -   Logique pour la Visualisation en coordonnees cartesiennes
C -
      LOGICAL  VISUXY
C -
C -   Logique pour la Visualisation en coordonnees cartesiennes
C -   dans la base d'orthotropie
C -
      LOGICAL  VISORT
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C -
C -    POUR LES DESSINS  viens de cominc.h le 26-07-96
C -
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C -
C *     titre  et sa longueur    
        INTEGER LTITRE , ILTIT , LONTIT 
        CHARACTER*100 TITRE
C -
C -   Epaisseur du trait pour de beaux dessins
C -
      INTEGER EPAISS
      INTEGER ZOBZOB
C -
C   si on est sur ecran noir et blanc
C   si on veut tracer les erreurs
C   si on veut tracer les fonctions simples de controle ( appel a DESFON )
C   si on veut tracer les fonction du temps en contraintes
C   si on veut tracer les fonction du temps en deformation
C   si on veut tracer l'ancienne version de dessin 
C   si on veut tracer le paramatre de charge
C   si on veut tracer le nombre d'iterations
C   si on veut tracer toute la courbe parametre de charge
C   si on le fichier paramc a deja ete rempli  pour cette exemple
       LOGICAL LERREU 
       LOGICAL LDFONC 
       LOGICAL LDFTCO 
       LOGICAL LDFTEP 
       LOGICAL ANCDES 
       LOGICAL LDPARC 
       LOGICAL LDTAUX 
       LOGICAL LDTPAR
       LOGICAL LITPAC
       LOGICAL COUTEM
       LOGICAL INTTEM         
C -POUR CHAPPG
c LOGR0=.t. => base er, e0, ez sinon , e1,e2, ez (base d'orthotropie)
c COUCHES : LCPLAN ( S11,S22,S12 ) LNORM( S23,S13,S33 ) LCNONL (d,d', EPSp)
c INTERFS : LIORM( S23,S13,S33 )   LINONL (d1,d2,d3', Up)
c LDIFCO indique si on veut une visualisation des differences chapeau admissible 
c si DETAIL on trace alors etape globale par etape globale
c
       LOGICAL LOGR0, LCPLAN, LCNORM, LCNONL, LINORM, LINONL 
       LOGICAL LDIFCO, DETAIL
C -
C * DESSIN SUR BENSON SI TBENS
      LOGICAL TBENS 
C * PAS DE DESSIN SUR ECRAN SI BATCH
      LOGICAL BATCH
C * SI NOIR ET BLANC
      LOGICAL BLACKW
C -
C -
      COMMON /VISUTEMPS/ FTEPMX , FTSIMX , DBFTEP , DBFTSI ,
     &                   DBCHEP , DBCHSI , FIFTEP , FIFTSI ,
     &                   FICHEP , FICHSI , NUETGL
C -
      COMMON /TYPTRACE/  LPOSCR , LPSCOL , ANGLAI 
C -
      COMMON /VISUESSAI/XMINPL , XMAXPL  , YMINPL , YMAXPL 
      COMMON /VISUDES  /LONGDE , LARGDE  , VISES , VISEST
      COMMON /VICOURBE  / NORMAR , NORMAT , NORMAZ , NORMAC,
     &                    NORMAD , NORMAU , NORABS , NORORD
C -
      COMMON /TYPVIS/LVISC, LPLAST, LENDOM, LMODUL, LVERRE, LDSIGO,
     &               LVISI, LPLASI, LENDIN, LENDIP, LSUDEL, LVERRI,
     &               LDSIGN, VDIFFE, DEJPAS, ERTOTA, VIENDO, CRIINT,
     &               CRICOU, VIPLAS, LTOTAL
C -
      COMMON /CHAPPG/  LOGR0, LCPLAN, LCNORM, LCNONL, LINORM, LINONL, 
     &                 LDIFCO , DETAIL
c
      COMMON /VISUXY/ VISUXY , VISORT
      COMMON /ETITRE/    ILTIT , LONTIT 
      COMMON /CTITRE/    TITRE
      COMMON / LOGDES  /BLACKW , LERREU , 
     &                  LDFONC , LDFTCO  , LDFTEP , ANCDES , LDPARC ,
     &                  LDTAUX , LDTPAR   , LITPAC , 
     &                  COUTEM , INTTEM  , TBENS , BATCH 

