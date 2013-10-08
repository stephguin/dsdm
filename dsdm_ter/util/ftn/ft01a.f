C     FT01AD        06/02/79
C     NAME FT01AD(R)                 CHECK
C 
      SUBROUTINE FT01AD(IT,INV,TR,TI)                                          1
C 
C     STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)                        2
C 
      DOUBLE PRECISION  TH,TI,TR,UI,UM,UR,WI,WR,WS,Z,ZI,ZR                     3
C 
C     PHIL AJOUT DE I I0 ......
C 
      INTEGER i ,  I0 , i1 , i2 , i3 , kk , k  , kk1 , l , l1
      INTEGER II , J , J0 , J1 , J2
      INTEGER IT , INV , KJUMP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='FT01AD')
C 
C     THIS ROUTINE CALCULATES THE FOURIER TRANSFORM OF EQUALLY SPACE           5
C     F(N)  N=0,1,...,IT-1                                                     6
C     THE DATA IS TAKEN TO BE PERIODIC IE.   F(N+IT) = F(N)                    7
C     ++++++ ARGUMENTS SET BY THE CALLING PROGRAM ++++++                       8
C     IT IS THE PROBLEM SIZE AND MUST BE A POWER OF 2                          9
C     INV = 2 FOR DIRECT TRANSFORM IE.                                        10
C     G(M) = SUM OVER N=0,1,..,IT-1 OF F(N)*EXP(2PI*SQRT(-1)*N*M/IT)          11
C     FOR M=0,1,...,IT-1                                                      12
C     INV = 1 FOR INVERSE TRANSFORM IE.                                       13
C     F(N) = (1./IT)*(SUM OVER M=0,1,..,IT-1 OF G(M)*EXP(-2PI*SQRT(-1)*       14
C     FOR N =0,1,...,IT-1                                                     15
C     TR(I)    I=1,2,..,IT MUST CONTAIN REAL PART OF DATA                     16
C     TI(I)    I=1,2,..,IT MUST CONTAIN THE IMAGINARY PART OF DATA            17
C     ++++++ ARGUMENTS SET BY ROUTINE ++++++                                  18
C     IF IT IS NOT A POWER OF 2 INV IS SET TO -1 FOR ERROR RETURN             19
C     TR(I)    I=1,2,..,IT IS SET TO REAL PART OF TRANSFORM                   20
C     TI(I)    I=1,2,..,IT IS SET TO THE IMAGINARY PART OF TRANSFORM          21
C     THE METHOD USED IN THIS ROUTINE IS DISCRIBED IN                         22
C     (GENTLEMAN AND SANDE, PROC. FALL JOINT COMPUTER CONFER. 1966)           23
C 
C     MODIF
C     DIMENSION TR(4),TI(4),UR(15),UI(15) 
C                                                                             24
      DIMENSION TR(*),TI(*),UR(15),UI(15)                                     24
C     TD
      DOUBLE PRECISION PIOLIV,DEUXPI
C     TD     DATA KJUMP/1/                                                    25
C     TD
C     vd$r noconcur
      KJUMP = 1
C     TD
      GO TO (100,200),KJUMP                                                   26
C     TD  100 UM=0.5D0                                                        27
  100 CONTINUE
      PIOLIV = DACOS(-1.D0)
      DEUXPI = 2.D0 * PIOLIV
      UM=0.5D0                                                                27
      DO 50 I=1,15                                                            28
        UM=0.5D0*UM                                                           29
C     MODIF OLIVIER      TH=6.283185307178D0*UM                               30
C     TD     TH=2.D0*DACOS(-1.D0)*UM                                          30
        TH= DEUXPI * UM                                                       30
        UR(I)=DCOS(TH)                                                        31
        UI(I)=DSIN(TH)                                                        32
   50 CONTINUE                                                                32
      KJUMP=2                                                                 33
  200 UM=1.0D0                                                                34
      GO TO(1,2),INV                                                          35
    1 UM=-1.0D0                                                               36
    2 I0=2                                                                    37
      DO 3 I=2,16                                                             38
      I0=I0+I0                                                                39
      IF(I0-IT)3,4,5                                                          40
    3 CONTINUE                                                                41
C     ERROR IN IT - SET INV=-1 AND RETURN                                     42
    5 INV=-1                                                                  43
      RETURN                                                                  44
C     IT= 2**I - INITIALISE OUTER LOOP                                        45
    4 I0=I                                                                    46
      II=I0                                                                   47
      I1=IT/2                                                                 48
      I3=1                                                                    49
C     START MIDDLE LOOP                                                       50
   10 K=0                                                                     51
      I2=I1+I1                                                                52
C     CALCULATE TWIDDLE FACTOR E(K/I2)                                        53
   11 WR=1.                                                                   54
      WI=0.                                                                   55
      KK=K                                                                    56
      J0=I0                                                                   57
   24 IF(KK)21,22,21                                                          58
   21 J0=J0-1                                                                 59
      KK1=KK                                                                  60
      KK=KK/2                                                                 61
      IF(KK1-2*KK)23,21,23                                                    62
   23 WS=WR*UR(J0)-WI*UI(J0)                                                  63
      WI=WR*UI(J0)+WI*UR(J0)                                                  64
      WR=WS                                                                   65
      GO TO 24                                                                66
   22 WI=WI*UM                                                                67
C     START INNER LOOP                                                        68
      J=0                                                                     69
C     DO 2*2 TRANSFORM                                                        70
   31 L=J*I2+K                                                                71
      L1=L+I1                                                                 72
      ZR=TR(L+1)+TR(L1+1)                                                     73
      ZI=TI(L+1)+TI(L1+1)                                                     74
      Z=WR*(TR(L+1)-TR(L1+1))-WI*(TI(L+1)-TI(L1+1))                           75
      TI(L1+1)=WR*(TI(L+1)-TI(L1+1))+WI*(TR(L+1)-TR(L1+1))                    76
      TR(L+1)=ZR                                                              77
      TR(L1+1)=Z                                                              78
      TI(L+1)=ZI                                                              79
C     INDEX J LOOP                                                            80
      J=J+1                                                                   81
      IF(J-I3)31,12,12                                                        82
C     INDEX K LOOP                                                            83
   12 K=K+1                                                                   84
      IF(K-I1)11,6,6                                                          85
C     INDEX OUTER LOOP                                                        86
    6 I3=I3+I3                                                                87
      I0=I0-1                                                                 88
      I1=I1/2                                                                 89
      IF(I1)51,51,10                                                          90
C     UNSCRAMBLE                                                              91
   51 J=1                                                                     92
      UM=1.                                                                   93
      GO TO(61,52),INV                                                        94
   61 UM=1./FLOAT(IT)                                                         95
   52 K=0                                                                     96
      J1=J                                                                    97
      DO 53 I=1,II                                                            98
      J2=J1/2                                                                 99
      K=2*(K-J2)+J1                                                          100
   53 J1=J2                                                                  101
   54 IF(K-J)66,56,55                                                        102
   56 TR(J+1)=TR(J+1)*UM                                                     103
      TI(J+1)=TI(J+1)*UM                                                     104
      GO TO 66                                                               105
   55 ZR=TR(J+1)                                                             106
      ZI=TI(J+1)                                                             107
      TR(J+1)=TR(K+1)*UM                                                     108
      TI(J+1)=TI(K+1)*UM                                                     109
      TR(K+1)=ZR*UM                                                          110
      TI(K+1)=ZI*UM                                                          111
   66 J=J+1                                                                  112
      IF(J-IT+1)52,57,57                                                     113
   57 TR(1)=TR(1)*UM                                                         114
      TI(1)=TI(1)*UM                                                         115
      TR(IT)=TR(IT)*UM                                                       116
      TI(IT)=TI(IT)*UM                                                       117
      RETURN                                                                 118
      END                                                                    119
