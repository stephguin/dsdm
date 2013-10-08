      SUBROUTINE TABSUP( LONG , TAB , SUPTAB )
C -
C - Cette routine renvoie le taaleau sup du tableau double TAB 
C -
C -
C On envoie comme arguments:
C E................LONG    longueur de TAB 
C E................TAB     tableau double 
C
C Et on recupere:
C S................ SUPTABL tableau sup du tableau TAB 
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
C
C
      INTEGER       LONG
C
      DOUBLE PRECISION  TAB( LONG) , SUPTAB( LONG )
C
C**********************************************************************
C
C *   Declaration des parametres locaux
C     """""""""""""""""""""""""""""""""
C
      INTEGER I
C
      DOUBLE PRECISION  MAXVAL
C
C______________________________________________________________________
C -
      MAXVAL = 0.D0
C -
      DO I = 1, LONG
C -
        MAXVAL    = MAX( MAXVAL , TAB(I) )
        SUPTAB(I) = MAXVAL
C -
      END DO
C -
      RETURN
      END     
C
      SUBROUTINE NORSUP( LONG , TAB , SUPTAB )
C -
C - Cette routine la valeur sup du dabs( TAB  )
C -
C -
C On envoie comme arguments:
C E................LONG    longueur de TAB 
C E................TAB     tableau double 
C
C Et on recupere:
C S................ SUPTABL tableau sup du tableau TAB 
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
C
C
      INTEGER       LONG
C
      DOUBLE PRECISION  TAB( LONG) , SUPTAB
C
C**********************************************************************
C
C *   Declaration des parametres locaux
C     """""""""""""""""""""""""""""""""
C
      INTEGER I
C
      DOUBLE PRECISION  MAXVAL
C
C______________________________________________________________________
C -
      SUPTAB = 0.D0
C -
      DO I = 1, LONG
C -
        SUPTAB    = MAX( SUPTAB , DABS(TAB(I)) )
C -
      END DO
C -
      RETURN
      END     
C -                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C -
      SUBROUTINE SUPTAB( LONG , TAB , MAXVAL )
C -
C - Cette routine renvoie le sup du tableau double TAB 
C -
C -
C On envoie comme arguments:
C E................LONG    longueur de TAB 
C E................TAB     tableau double 
C
C Et on recupere:
C S................ MAXVAL sup du tableau double TAB 
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
      INTEGER       LONG
C
      DOUBLE PRECISION  TAB( LONG) , MAXVAL
C
C**********************************************************************
C
C *   Declaration des parametres locaux
C     """""""""""""""""""""""""""""""""
C
      INTEGER I
C
C______________________________________________________________________
C -
      MAXVAL = 0.D0
C -
      DO I = 1, LONG
C -
        MAXVAL = MAX( MAXVAL , DABS( TAB(I) ) )
C -
      END DO
C -
      RETURN
      END
C
      SUBROUTINE MUMACV( DIM , MAT , VECTEN , VECSOR )
C multiplication d'une matrice carre symetrique par un vecteur
C On envoie comme arguments:
C E................DIM     Dimension de la matrice carre
C E................MAT     Matrice carre
C E................VECTEN  Vecteur en entree
C
C Et on recupere:
C S................VECSOR  Vecteur en sortie
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
C
      INTEGER       DIM
C
      DOUBLE PRECISION MAT(DIM,DIM), VECTEN(DIM) , VECSOR(DIM)
C
C**********************************************************************
C
C *   Declaration des parametres locaux
C     """""""""""""""""""""""""""""""""
C
      INTEGER         I , J 
      DOUBLE PRECISION VEC
C
C
C***********************************************************************
Cvd$:r noconcur
Cvd$:r novector
      DO J =1, DIM
        VEC =0.D0
        DO I =1 , DIM
          VEC = VEC + MAT(I,J)*VECTEN(I)
        END DO
        VECSOR(J) = VEC
      END DO
C***********************************************************************
      RETURN
      END
