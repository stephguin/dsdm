C -
C * LOGIQUE INDIQUANT LES OPTIONS DE COMPORTEMENT
C -   LENDCO      On tient compte de l'endommagement des couches
C -   LRUPCO      On tient compte de la rupture des fibres
C -   LPLACO      On tient compte de la plasticite des couches
C -   LENDIN      On tient compte de l'endommagement des interfaces
C -   LPLAIN      On tient compte de la plasticite des interfaces
C -   D2D3        = 2. calcul "plaque"
C -   D2D3        = 3. calcul "3D", participation des contraintes HP
C -                    aux forces d'endommagement
C -
      LOGICAL    LENDCO , LRUPCO , LPLACO , LENINT , LPLAIN
C -   
C * LOGIQUE INDIQUANT LE TYPE DE COMPORTEMENT
C -   LMOSIG      modele en contrainte
C -   LRETAR      modele avec retard 
C -   LRETIN      modele avec retard apres instabilite
C -   LPLAEF      deformation plastique effective
C -   LPLEFF      deformation plastique finie
C *                                                        
      LOGICAL   LMOSIG , LRETAR , LRETIN , LPLAEF , LPLEFF
c -
      COMMON/ TYPCAL / LENDCO, LRUPCO, LPLACO, LENINT, LPLAIN
      COMMON/ MODELE / LMOSIG, LRETAR, LRETIN, LPLAEF, LPLEFF
C *
