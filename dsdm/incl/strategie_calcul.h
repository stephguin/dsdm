C -
C   pour passer par l'etape preliminaire
C -
      LOGICAL LETAPP
C -
C   pour passer a une resolution par valeurs propres dans l'etape globale
C -
      LOGICAL VALPRO
C -
C   pour passer a une resolution par valeurs propres dans l'etape globale
C -
      INTEGER NBVALP
C -
C   pour imposer le nombre de fonctions du temps dans le cas VALPRO = .T.
C -
      LOGICAL UPERIO
C -
C   pour imposer ne diviser qu'a l'instabilite
C -
      LOGICAL UNINST
C -
C   nombre d'iterations max du probleme statique
C -
      INTEGER NBSTAM
C -
C   nombre d'iterations max du probleme cinematique
C -
      INTEGER NBCINM
c -
C   logical pour une reprise avec grand nettoyage
C -
      LOGICAL BIGNET
c -
C   SYMX indique si la solution est symetrique par rapport a x
C -
      LOGICAL SYMX
C -
C   SYMY indique si la solution est symetrique par rapport a y
C -
      LOGICAL SYMY
C -
C   SYMO indique si la solution est symetrique par rapport a l'origine
C -
      LOGICAL SYMO                                            
c -
C   FORS on force la symetrie ( en imposant la symetrie des efforts ( meme residuels))
C -
      LOGICAL FORS
c -                                            
      COMMON /STRCAL/ LETAPP, VALPRO, UPERIO, UNINST, BIGNET
c -
      COMMON /DIMSTR/ NBCINM, NBSTAM, NBVALP
C *
      COMMON/SYMCAL/ SYMX, SYMY, SYMO, FORS
