C 
C    Integration numerique par les trapezes de la fonction de Hertz
C 
      DOUBLE PRECISION FUNCTION INT_FHERTZ (MILG, RM, RAYCONT, FBASE)
C 
      external  FBASE
C 
      DOUBLE PRECISION SOM, MILG, RM, DX, RAYCONT, FHERTZ, FBASE
      INTEGER          PAS,I
C 
      PAS = 100
      DX = 2.D0/PAS
C 
      SOM = FHERTZ (-1.D0, MILG, RM, RAYCONT, FBASE)/2.D0
      DO I = 1, PAS-1
        SOM = SOM + FHERTZ (I*DX-1.D0, MILG, RM, RAYCONT, FBASE)
      ENDDO
      SOM = SOM + FHERTZ (1.D0, MILG, RM, RAYCONT, FBASE)/2.D0
C 
      INT_FHERTZ = SOM*DX*MILG
C 
      RETURN
C 
      END
