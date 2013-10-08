C     Fonction representant les forces exterieures de Hertz
C 
       DOUBLE PRECISION FUNCTION FHERTZ (R, MILG, RM, RAYCONT, FBASE)
C 
       external  FBASE
C 
       DOUBLE PRECISION R, MILG, RM, K, RAYCONT, FBASE
C 
       K = ((R*MILG + RM)/RAYCONT)**2
       IF (K .LT. 1.) THEN
          FHERTZ = DSQRT(1.D0-K) * FBASE(R) * (R*MILG+RM)
       ELSE
          FHERTZ = 0.D0
       ENDIF
C 
       RETURN
C 
       END
