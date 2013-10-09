C     On envoie comme arguments :
C 
C     E ...... MAT   Matrice stockee profil dont on cherche le + gd
C                    terme diagonal
C     E ...... DIM   Dimension de cette matrice
C     E ...... NP    Tableau profil
C 
C     Et on recupere :
C 
C     S ...... VALSUP Valeur superieure des termes de la diagonale
C 
      SUBROUTINE TERSUP (MAT, DIM, NP, VALSUP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER   DIM , NP(DIM+1)
      DOUBLE PRECISION  MAT(NTMAT), VALSUP
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TERSUP')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      VALSUP=0.D0
      DO I = 2, DIM+1
        VALSUP = DMAX1( VALSUP, MAT(NP(I)))
      ENDDO
C 
CD    CALL IMPDP('VALEUR DU TERME + GD TERME DIAG '//IDPROG,VALSUP)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
