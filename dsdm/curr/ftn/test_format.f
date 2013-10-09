C===============================================
C    Test si l'exposant du format des donnees
C       n'est pas trop petit (<1.e-99)
C              ou trop grand (>1.e+99)
C===============================================
C
       SUBROUTINE   TEST_FORMAT (DIM,TAB)
C
       INTEGER            DIM, I
       DOUBLE  PRECISION  TAB(DIM)
C
       DO  I=1, DIM
          IF (DABS(TAB(I)) .LT. .1D -99)  TAB(I) = 0.D0
          IF (DABS(TAB(I)) .GT. .1D +99)  TAB(I) = .1D99
       ENDDO
C
       RETURN
       END
