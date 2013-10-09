C     Caracteristiques non lineaires pour la couche .
C 
C     On envoie comme arguments :
C 
C     E ...... YO          valeur de Y initiale
C     E ...... YC          valeur de Y critique
C     E ...... B           valeur du couplage d'endommagement
C     E ...... E22         valeur initiale du module E22
C     E ...... G12         valeur initiale du module 2*G12
C     E ...... EPSORT (3)  valeur des deformations planes dans la
C                          base d'orthotropie (actuelle)
C     Et on recupere :
C 
C     S ...... YBNOU       nouvelle valeur de Y barre
C 
      SUBROUTINE YELAST (CARNLI, EPSORT, YBNOU)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION  CARNLI( 5 ) , EPSORT(3)
      DOUBLE PRECISION   YBNOU
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  PPOS , INTER
C 
CD    LOGICAL LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='YELAST')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF( LTRACN(1) ) THEN
CD      CALL IMPDN(' Valeur de s22 '  ,  EPSORT(2)  )
CD      CALL IMPDN(' Valeur de s12 '  ,  EPSORT(3)  )
CD      CALL IMPDN(' Valeur de E22 '  ,  CARNLI(4)  )
CD      CALL IMPDN(' Valeur de 2G12 '  , CARNLI(5)  )
CD    END IF
C 
      PPOS    = DDIM( EPSORT(2) , 0.D0 )
C 
      INTER   = EPSORT(3)*EPSORT(3)*CARNLI(5)
     &        + CARNLI(3)*PPOS*PPOS*CARNLI(4)
      YBNOU   = DSQRT( INTER )
C 
CD    IF( LTRACN(1) ) THEN
CD      CALL IMPDN(' Valeur de ybnou  ' , YBNOU  )
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Caracteristiques non lineaires pour l'interface.
C 
C     Le tableau CARNLI contient sequentiellement :
C 
C     E ...... YC         valeur  de Y  critique
C     E ...... GAM1       valeur du couplage d'endommagement 1
C     E ...... GAM2       valeur du couplage d'endommagement 2
C     E ...... E10        valeur initiale du module E1
C     E ...... E20        valeur initiale du module E2
C     E ...... E30        valeur initiale du module E3
C     E ...... SAUORT (3) valeur des sauts plan dans la
C                          base d'orthotropie (actuelle)
C 
C     Et on recupere :
C 
C     S ...... YBNOU      nouvelle valeur de Y barre

      SUBROUTINE YELASI (CARNLI, SAUORT, YBNOU)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION  CARNLI( 6 ) , SAUORT(3)
      DOUBLE PRECISION  YBNOU
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  PPOS , INTER
C 
CD    LOGICAL LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='YELASI')
C 
CD    CALL WLKBCD(IDPROG)
C 
C ----------------------------------------------------------------------- 
CD    IF( LTRACN(1) ) THEN
CD      CALL IMPDN(' Valeur de s11 '  ,  SAUORT(1)  )
CD      CALL IMPDN(' Valeur de s22 '  ,  SAUORT(2)  )
CD      CALL IMPDN(' Valeur de s33 '  ,  SAUORT(3)  )
CD    END IF
C 
      PPOS    = DDIM( SAUORT(3) , 0.D0 )
C 
      INTER   = CARNLI(4)*SAUORT(1)*CARNLI(2)*SAUORT(1)
     &        + CARNLI(5)*SAUORT(2)*CARNLI(3)*SAUORT(2)
     &        + CARNLI(6)*PPOS*PPOS
C 
      YBNOU   = .5D0*INTER
C 
CD    IF( LTRACN(1) ) THEN
CD      CALL IMPDN(' Valeur de ybnou  '  , YBNOU  )
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
