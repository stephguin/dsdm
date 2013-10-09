C     On envoie comme arguments :
C 
C     E ...... N        Taille de A decomposee CHOLESKY A=LU
C     E ...... NP       Pointeur diagonal
C     E ...... A        Termes de A stockes profil
C 
C     ES ...... F       Second membre
C                       On peut faire F=U en entree
C     ES ...... U       Solution

      SUBROUTINE DESREM (N, NP, A, F, U)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER         N, NP(N+1)
C 
      DOUBLE PRECISION  A(NTMAT), F(N), U(N)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     NBTI, JDEB, NADRI, I
      INTEGER     P, JFIN, IR, NBTIR, NADRIR, IRDEB, IRFIN, IL
C 
      DOUBLE PRECISION   X
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DESREM')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     DESCENTE
C 
CD    CALL IMPTDP ('SECOND MEMBRE DANS : '//IDPROG, F(1), 1, N)
C 
      DO I = 1, N
        X      = F(I)
        NBTI   = NP(I+1)-NP(I)
        JDEB   = I-NBTI+1
        NADRI  = NP(I+1)-I
        JFIN   = I-1
        IF (JDEB .LE. JFIN) THEN
          DO P = JDEB, JFIN
            X  = X-A(NADRI+P)*U(P)
          ENDDO
        ENDIF
        U(I)   = X/A(NADRI+I)
      ENDDO
C 
C     REMONTEE
C 
      DO I = 1, N
        IR     = N-I+1
        NBTIR  = NP(IR+1)-NP(IR)
        NADRIR = NP(IR+1)-IR
        IRDEB  = IR-NBTIR+1
        U(IR)  = U(IR)/A(NADRIR+IR)
        X      = U(IR)
        IF (IRDEB .LT. IR) THEN
          IRFIN    = IR-1
          DO IL = IRDEB, IRFIN
C 
C     MODIF APRES = NADRI???==>NADRIR
C 
            U(IL)     = U(IL)-A(NADRIR+IL)*X
          ENDDO
        ENDIF
      ENDDO
C 
CD    CALL IMPTDP ('SOLUTION DANS : '//IDPROG, U(1), 1, N)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
