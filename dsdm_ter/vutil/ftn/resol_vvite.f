C     Les seconds membres sont stockes sous la forme (NBSM, N).
C     En entree on peut confondre SM et SOL.
C 
C     On envoie comme arguments :
C 
C     E ...... N     Taille de A decomposee CHOLESKY A=LU
C     E ...... PRO   Pointeur diagonal
C     E ...... A     Termes de A stockes profil
C     E ...... NBSM  Nombre de seconds membres
C 
C     Et on recupere :
C 
C     ES...... SM    Second membre; on peut faire F=U en entree
C     ES...... SOL   Solution
C 
      SUBROUTINE DEREPP (N, PRO, A, NBSM, SM, SOL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      INTEGER N,NBSM
      INTEGER PRO(N+1)
      DOUBLE PRECISION A(*),SM(NBSM,N),SOL(NBSM,N)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER I,J,K
      INTEGER KDEB,KFIN
      INTEGER ADRL,ADR,IPRO,NBT
      DOUBLE PRECISION ZZ,UNSURA
C 
C -----------------------------------------------------------------------
C     vd$r noconcur
C 
C     DESCENTE
C 
      DO I=1,N
        ADRL = PRO(I)
        IPRO = PRO(I+1)
        NBT =  IPRO - ADRL -1
C 
C       WRITE( 6 , * ) 'A(IPRO)' , A(IPRO)
C 
        UNSURA = 1.D0/A(IPRO)
C 
        KDEB = I-NBT
        KFIN = I-1
        DO J=1,NBSM
          ADR = ADRL
          ZZ = 0.D0
          DO K = KDEB , KFIN
            ADR = ADR + 1
            ZZ = ZZ + A(ADR)*SOL(J,K)
          ENDDO
          SOL(J,I) = (SM(J,I) - ZZ) * UNSURA
        ENDDO
      ENDDO
C 
C     REMONTEE
C 
      DO I=N,1,-1
        ADRL = PRO(I)
        IPRO = PRO(I+1)
        NBT =  IPRO - ADRL -1
        UNSURA = 1.D0/A(IPRO)
C 
C     vd$l novector
C 
        DO J=1,NBSM
          SOL(J,I) = SOL(J,I) * UNSURA
        ENDDO
C 
        KDEB = I-NBT
        KFIN = I-1
        DO J=1,NBSM
          ADR = ADRL
          ZZ = SOL(J,I)
          DO K =  KDEB ,KFIN
            ADR = ADR + 1
            SOL(J,K) = SOL(J,K) - A(ADR)*ZZ
          ENDDO
        ENDDO
      ENDDO
C 
      RETURN
      END
