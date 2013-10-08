C     On envoie comme arguments :
C 
C     E ...... n      dimension de la matrice symetrique a decomposer profil
C     E ...... np     pointeur diagonal de a (tableau du profil par ddl)
C     E ...... EPS    precision test provient du tableau DP
C 
C     Argument entree-sortie (eventuellement) :
C 
C     ES ..... A==>L  A matrice a decomposer profil on recupere en
C                     sortie L a la place de A c-a-d la matrice
C                     triangulaire inferieure
C 
C     Et on recupere :
C 
C     S ...... NPIV   Nombre de pivots inferieurs au test EPS
C 
      SUBROUTINE OCHOLS (N, NP, EPS, A, L, NPIV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER  N, NP(N+1), NPIV
C 
      DOUBLE PRECISION  A(NTMAT), L(NTMAT), EPS
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     NBTI, JDEB, NADRI, NADRJ, NBTJ, PDEB, PFIN, I
      INTEGER     J, P
C 
      DOUBLE PRECISION X
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='OCHOLS')
C 
      CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C     CALL MESSAO ('ON ENTRE DANS '//IDPROG)
C 
C     Pour verification : impression de A et NP en entree
C 
CD    CALL IMPTET ('TABLEAU PROFIL EN ENTREE DANS '//IDPROG, NP(1),
CD   &              1, N+1)
CD    CALL IMPTDT ('MATRICE A STOCKEE PROFIL EN ENTREE
CD   &             DANS '//IDPROG, A(1), 1, NP(N+1))
C 
C     Calcul du 1er terme de la matrice triangulaire inferieure
C     Verification du pivot apres initialisation du nb de pivot
C 
      NPIV  = 0
C 
      CALL TESPIV (A(1), 1, EPS, NPIV)
      CALL IMPET ('NOMBRE DE PIVOTS OU PIV < TEST', NPIV)
C 
      L(1)   = DSQRT(A(1))
      IF (N .EQ. 1) THEN
        GOTO 1
      ENDIF
C 
C     Calcul des termes de 2 au nombre total de termes stockes
C 
      DO I = 2, N
C 
C      if (i .eq. 56) stop
C       CALL IMPET ('I dans '// IDPROG, I)
C 
        NBTI = NP(I+1)-NP(I)
C 
C       CALL IMPET ('NBTI dans '// IDPROG, NBTI)
C 
        JDEB = I-NBTI+1
C 
C 	CALL IMPET ('JDEB dans '// IDPROG, JDEB)
C 
        NADRI = NP(I+1)-I
C 
C 	CALL IMPET ('NADRI dans '// IDPROG, NADRI)
C 
        DO J = JDEB, I
C 
C         CALL IMPET ('J dans '// IDPROG, J)
C 
          X = A(NADRI+J)
C 
C         CALL IMPDT ('X dans '// IDPROG, X)
C 
          NADRJ = NP(J+1)-J
C 
C         CALL IMPET ('NADRJ dans '// IDPROG, NADRJ)
C 
          IF (J .GT. JDEB) THEN
C 
C           CALL IMPET ('J dans '// IDPROG, J)
C 
            NBTJ  = NP(J+1)-NP(J)
C 
C           CALL IMPET ('PDEB dans '// IDPROG, PDEB)
C 
            PDEB  = MAX0 (JDEB,J-NBTJ+1)
C 
C           CALL IMPET ('I dans '// IDPROG, I)
C 
            PFIN  = J-1
C 
C           CALL IMPET ('PFIN dans '// IDPROG, PFIN)
C 
            DO P = PDEB, PFIN
C 
C             CALL IMPDT ('X dans '// IDPROG, X)
C             CALL IMPDT ('L(NADRI+P) dans '// IDPROG, L(NADRI+P))
C             CALL IMPDT ('L(NADRJ+P) dans '// IDPROG, L(NADRJ+P))
C 
                 IF ( DABS( L(NADRI+P)) .LT. 1.D-50) GOTO 33
                X = X-L(NADRI+P)*L(NADRJ+P)
33              CONTINUE       		
            ENDDO
          ENDIF
          IF (J .LT. I) THEN
C 
C           MODIF APRES = A===>L
C           CALL IMPDT ('X dans '// IDPROG, X)
C           CALL IMPDT ('L(NADRI+J) dans '// IDPROG, L(NADRI+J))
C 
            IF ( X . LT . 1.D -50) GOTO 44
             L(NADRI+J) = X/L(NADRJ+J)
44	     L(nadri+j) = 0.D0
          ELSE
C 
c            CALL TESPIV (X, I, EPS, NPIV)
C 
C           CALL IMPDT ('L(NADRI+I) dans '// IDPROG, L(NADRI+I))
C 
            IF ( X . LT . 1.D -50) GOTO 22
            L(NADRI+I) = DSQRT(X)
22          L(NADRI+I) = 0. D0	    
C 
          ENDIF
C 
        ENDDO
C 
      ENDDO
C 
      IF (NPIV .GT. 0) THEN
C 
        CALL IMPEP ('LE NOMBRE DE PIVOTS INFERIEURS A LA PRECISION
     &               EST ', NPIV)
      ENDIF
C 
      CALL IMPTDP ('MATRICE TRIANGULAIRE EN SORTIE : '//IDPROG,
     &              L(1), 1, NP(N+1))
C 
      CALL RETOUD (IDPROG)
C 
1     RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     On envoie comme arguments :
C 
C     E ...... PIV     Valeur du pivot teste
C     E ...... NUM     Numero de ce pivot dans la matrice
C     E ...... TEST    Valeur test du pivot
C
C     Et on recupere :
C 
C     S ...... NPIV    Nombre de pivots inferieurs au test
C 
      SUBROUTINE TESPIV (PIV, NUM, TEST, NPIV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      INTEGER   NUM, NPIV
C 
      DOUBLE PRECISION   PIV, TEST
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TESPIV')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (PIV .LT. TEST) THEN
        NPIV  = NPIV+1
CD      CALL IMPET ('POUR LE PIVOT NUMERO DANS '//IDPROG, NUM)
CD      CALL IMPDT ('LA VALEUR DE CE PIVOT EST ', PIV)
      ENDIF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     On envoie comme arguments :
C 
C     E ...... n      dimension de la matrice symetrique a decomposer profil
C     E ...... np     pointeur diagonal de a (tableau du profil par ddl)
C     E ...... EPS    precision test provient du tableau DP
C 
C     Argument entree-sortie (eventuellement) :
C 
C     ES ..... A==>L  A matrice a decomposer profil on recupere en
C                     sortie L a la place de A c-a-d la matrice
C                     triangulaire inferieure
C 
C     Et on recupere :
C 
C     S ...... NPIV   Nombre de pivots inferieurs au test EPS
C 
      SUBROUTINE CHOLSP (N, NP, EPS, A, L, NPIV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER  N, NP(N+1), NPIV
C 
      DOUBLE PRECISION  A(NTMAT), L(NTMAT), EPS
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     NBTI, JDEB, NADRI, NADRJ, NBTJ, PDEB, PFIN, I
      INTEGER     J, P
C 
      DOUBLE PRECISION X
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CHOLSP')
C 
      CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C     CALL MESSAO ('ON ENTRE DANS '//IDPROG)
C 
C     Calcul du 1er terme de la matrice triangulaire inferieure
C     Verification du pivot apres initialisation du nb de pivot
C 
      NPIV  = 0
C 
      CALL TESPIV (A(1), 1, EPS, NPIV)
C     CALL IMPET ('NOMBRE DE PIVOTS OU PIV < TEST', NPIV)
C 
      L(1)   = DSQRT(A(1))
      IF (N .EQ. 1)  GOTO 1
C 
C     Calcul des termes de 2 au nombre total de termes stockes
C 
      DO I = 2, N
        NBTI = NP(I+1)-NP(I)
        JDEB = I-NBTI+1
        NADRI = NP(I+1)-I
        DO J = JDEB, I
          X = A(NADRI+J)
          NADRJ = NP(J+1)-J
          IF (J .GT. JDEB) THEN
            NBTJ  = NP(J+1)-NP(J)
            PDEB  = MAX0 (JDEB,J-NBTJ+1)
            PFIN  = J-1
            DO P = PDEB, PFIN
                 IF ( DABS( L(NADRI+P)) .LT. 1.D-50) GOTO 33
                X = X-L(NADRI+P)*L(NADRJ+P)
33              CONTINUE       		
            ENDDO
          ENDIF
          IF (J .LT. I) THEN
             L(NADRI+J) = X/L(NADRJ+J)
          ELSE
            L(NADRI+I) = DSQRT(X)	    
          ENDIF
        ENDDO
      ENDDO
C 
C     IF (NPIV .GT. 0) THEN
C 
C       CALL IMPEP ('LE NOMBRE DE PIVOTS INFERIEURS A LA PRECISION
C    &               EST ', NPIV)
C     ENDIF
C 
C     CALL IMPTDP ('MATRICE TRIANGULAIRE EN SORTIE : '//IDPROG,
C    &              L(1), 1, NP(N+1))
C 
      CALL RETOUD (IDPROG)
C 
1     RETURN
      END
C      
