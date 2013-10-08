C                    .               n
C     Integration de D = k < Y - d >
C     Implicite la valeur de Y est donnee > On utilise Newton
C 
      SUBROUTINE INTRET (DELTAT, K, N, Y, DN, DAC, DMAX)
C      
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION  DELTAT, K, N, Y, DN, DAC, DMAX
      LOGICAL           LOGFRA
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  TEST, KLOC, DELTY, GAC, GINI, GPAC
      DOUBLE PRECISION  DELTD
      INTEGER           NBFOIS
C -----------------------------------------------------------------------
C 
C     INITIALISATION DES VALEURS DE D
C 
      DAC   =  DN
      TEST  =  Y-DAC
      DELTD =  0.D0
C 
      IF ((TEST .LE. 0.D0) .OR. (DN .EQ. DMAX)) GOTO 1	
C 
      KLOC  = K*DELTAT
C 
      IF (N .EQ. 1.D0) THEN
        DAC = (DN+KLOC*Y)/(1.D0+KLOC)
        DAC =  DMIN1(DAC, DMAX)
        GOTO 1
      END IF
C 
      IF (TEST .GT. 0.D0 .AND. (DAC .LT. DMAX)) THEN
C 
        KLOC  = K*DELTAT
        DELTY = Y-DN
C 
C       ACCROISSEMENT DE D
C 
        GINI   =  -KLOC*(DELTY**N)
        GAC    =   GINI
        NBFOIS =   0
C 
        DO WHILE ((DABS(GAC/GINI). GT. 1.D-5)
     &             .AND. NBFOIS .LT. 50)
C 
          GPAC   = 1.D0 + KLOC*N*((DELTY-DELTD)**(N-1.D0))
          DELTD  = DELTD - GAC/GPAC
          GAC    = DELTD - KLOC*((DELTY-DELTD)**N)
          NBFOIS = NBFOIS+1
C 
        END DO
C 
      END IF
C 
      DAC = DN + DELTD
      IF (DAC .GT. DMAX) DAC = DMAX
      IF (DAC .LT. DN) DAC = DN
C 
1     RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Integration de D = k < Y - d >
C     Implicite : la valeur de Y est donnee. On utilise Newton.
C 
      SUBROUTINE INRPRF (DELTAT, K, N, YRP, YRF, YS, DN, DAC, DMAX)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  DELTAT
      DOUBLE PRECISION  K, N, YRP, YRF, YS, DN, DAC, DMAX
      DOUBLE PRECISION  YFL, ACC1, ACC2, KF
C -----------------------------------------------------------------------
C 
C     INITIALISATION DES VALEURS DE D
C 
      CALL INTRET (DELTAT, K, N, YRP, DN, ACC2, DMAX)
      ACC2 = ACC2 - DN
      DAC = DMIN1 (DN+ACC2, DMAX)
C 
      RETURN
      END
