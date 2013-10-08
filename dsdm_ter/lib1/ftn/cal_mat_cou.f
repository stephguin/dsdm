C     On envoie comme arguments :
C 
C     E ....... A   demi-longueur de l'element
C     E ....... NB  dimension de la matrice a multiplier
C     E ....... MAT matrice d'entree( necdelamirement 12*12)
C 
C     Et on recupere :
C                    | 1   A   |
C     S ....... MATA=|         | * MAT
C                    | A   A*A |
C 
      SUBROUTINE MATLOC (NB, A, MAT, MATA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER I, J, NB, DEMI
C 
      DOUBLE PRECISION  A, MAT(NB, NB), MATA(NB, NB)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MATLOC')
C 
CD    CALL WLKBCD(IDPROG)
C 
CD    IF(MOD(NB,2).NE.0)THEN
CD      CALL ERREUD(0,
CD      'LA DIMENSION DE LA MATRICE DANS MATLOC EST IMPAIRE ')
CD    ENDIF
C 
      DEMI=NB/2
      DO I=1,DEMI
        DO J=1,DEMI
        MATA(I,J)=MAT(I,J)
        ENDDO
      ENDDO
      DO I=DEMI+1,NB
        DO J=1,DEMI
        MATA(I,J)=A*MAT(I,J)
        ENDDO
      ENDDO
      DO I=1,DEMI
        DO J=DEMI+1,NB
        MATA(I,J)=A*MAT(I,J)
        ENDDO
      ENDDO
      DO I=DEMI+1,NB
        DO J=DEMI+1,NB
        MATA(I,J)=A*A*MAT(I,J)
        ENDDO
      ENDDO
C 
C    CD CALL IMPTDP('RESULTAT DE LA MULTIPLICATION', MATA(1,1), NB, NB)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine calcule a partir des tableaux constants
C     contenus dans MAT-CONST les tableaux integres qui sont
C     necdelamires au calcul des matrices de rigidites.
C 
C     Les 10 premiers tableaux sont les sommes ponderees
C     (sur les points de gauss) des tableaux correspondant de
C     MAT-CONST 
C                                      _    _    _    _
C     Les 4 derniers sont les tableaux NxNx, NxNy, NyNx, NyNy,
C     multiplies par x sur chaque point de gauss et sommes
C     ensuite sur tous les points de gauss.
C 
C     Attention les tableaux ne sont pas multiplies par MATLOC,
C     ils doivent l'etre dans les routines CALKIJ.
C 
C     On calcule aussi les tableaux necdelamires au calcul des
C     matrices de rigidite des interfaces (MAT-INTERF).
C 
C     Attention les tableaux ne sont pas multiplies par MATLOC,
C     ils doivent l'etre dans les routines CALKI.
C 
C     On calcule aussi les tableaux necdelamires au calcul des
C     termes d'efforts (MAT-EFFORT).
C 
C     Attention les tableaux ne sont pas multiplies par MATLOC,
C     ils doivent l'etre dans les routines ASSBFI.

      SUBROUTINE TABINT
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER D0, L0, LONGTA, NTABL, NTERME
      INTEGER D1, R1 ,K, H, ADLC2
      INTEGER D0I, D1I, R1I, ADLC2I, N0, ADR0
      INTEGER I, J, L, R, NTERAV, ADR0I
      INTEGER AM2LC, ADM2LC
C 
      DOUBLE PRECISION XSI, WI, WJ, W, MULT
C 
CD    LOGICAL LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TABINT')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL NUTBDM('MAT-CONST ',N0)
C 
C     D0 EST LA PREMIERE ADRESSSE DE MAT-CONST
C 
      D0=ADM(N0)
      L0=LONGDM(N0)
C 
CD    CALL IMPEN('LONGUEUR DE MAT-CONST',L0)
C 
C     LONGTA EST LE NOMBRE DE TERMES STOCKES POUR LE CALCUL D'UN TABLEAU
C 
      LONGTA=XINTEG*YINTEG*144
C 
C     NTABL EST LE NOMBRE DE TABLEAUX STOCKES PAR POINT DE GAUSS
C 
      NTABL=L0/LONGTA
C 
CD    CALL IMPEN('NOMBRE DE TABLEAUX STOCKES PAR POINT DE GAUSS',NTABL)
C 
C     NTERME EST LE NOMBRE DE TERMES STOCKES PAR POINT DE GAUSS
C 
      NTERME=144*NTABL
C 
CD    CALL IMPEN(
CD    'NOMBRE DE TERMES STOCKES PAR TABLEAU ET PAR POINT DE GAUSS',NTERME)
C 
C -----------------------------------------------------------------------
C     Creation du tableau des termes integres
C 
C             nom : TAB-INTEGC
C             1ere adresse libre : D1
C             numero : N1
C -----------------------------------------------------------------------
      CALL GESTDP('TAB-INTEGC',(NTABL+4)*144,D1)
      CALL MENADM(D1,(NTABL+4)*144)
      R1=D1-1
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR0 pour
C     ranger le tableau des termes multiplies par les poids (ou par w*x)
C -----------------------------------------------------------------------
      CALL GSPOUD(144, ADR0)
C 
C     INTEGRATION DES TABLEAUX A L'AIDE DES POIDS
C 
      K=XINTEG*(XINTEG-1)/2
      H=YINTEG*(YINTEG-1)/2
C 
      ADLC2=D0
      DO I=1,XINTEG
        XSI=GAUSS(K+I)
        WI=POIDS(K+I)
        DO J=1,YINTEG
          WJ=POIDS(H+J)
          W=WI*WJ
C 
CD        CALL IMPDN('VALEUR DE W '//IDPROG,W)
C 
          MULT=W*XSI
C 
CD        CALL IMPDN('VALEUR DE MULT '//IDPROG,MULT)
C 
          NTERAV=((I-1)*YINTEG+J-1)*NTERME
          DO L=1,NTABL
C 
C     TCD CALL MUMARE(W,144,DM(ADLC2),DM(ADR0))
C 
            CALL HOMAT(W,DM(ADLC2),DM(ADR0),144,1)
C 
C     TCD CALL ADDMAD(144,DM(D1+(L-1)*144),DM(ADR0),DM(D1+(L-1)*144))
C 
            CALL ADD(DM(D1+(L-1)*144) , DM(ADR0) , DM(D1+(L-1)*144) ,
     &                144 , 1 )
            ADLC2=ADLC2+144
          ENDDO
          DO R=1,4
C 
C     TCD CALL MUMARE(MULT,144,DM(D0+(4+R)*144+NTERAV),DM(ADR0))
C     TCD CALL ADDMAD(144,DM(D1+(NTABL-1+R)*144),DM(ADR0),DM(D1+(NTABL-1+R)*144))
C 
            CALL HOMAT( MULT , DM(D0+(4+R)*144+NTERAV) , DM(ADR0) ,
     &                  144 , 1 )
            CALL ADD( DM(D1+(NTABL-1+R)*144) , DM(ADR0) ,
     &                 DM(D1+(NTABL-1+R)*144) , 144 , 1 )
          ENDDO
        ENDDO
      ENDDO
C 
      IF(NBINT.GT.0)THEN
C -----------------------------------------------------------------------
C 
C     S'IL Y A DES INTERFACES
C 
C     Creation du tableau des termes integres pour les interfaces
C             nom : TAB-INTEGI
C             1ere adresse libre : D1I
C             numero : N1I
C -----------------------------------------------------------------------
        CALL ADTBDM('MAT-INTERF',D0I)
C 
C     D0I EST LA PREMIERE ADRESSE DE MAT-INTERF
C     LE PREMIER TABLEAU EST LE TABLEAU INTEGRE POUR LES INTERFACES
C     LE DEUXIEME TABLEAU EST LE TABLEAU INTEGRE*x POUR LES INTERFACES
C 
        CALL GESTDP('TAB-INTEGI',2*64,D1I)
        R1I=D1I-1
        CALL MENADM(D1I,2*64)
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR0I pour
C     ranger le tableau integre des termes multiplies par les poids (ou par w*x)
C -----------------------------------------------------------------------
        CALL GSPOUD(64, ADR0I)
        ADLC2I=D0I
C 
        DO I=1,XINTEG
          XSI=GAUSS(K+I)
          WI =POIDS(K+I)
          MULT=WI*XSI
C 
CD        CALL IMPDN('VALEUR DE MULT '//IDPROG,MULT)
C 
C     TCD CALL MUMARE(WI,64,DM(ADLC2I),DM(ADR0I))
C     TCD CALL ADDMAD(64,DM(D1I),DM(ADR0I),DM(D1I))
C     TCD CALL MUMARE(MULT,64,DM(ADLC2I),DM(ADR0I))
C     TCD CALL ADDMAD(64,DM(D1I+64),DM(ADR0I),DM(D1I+64))
C 
          CALL HOMAT(WI,DM(ADLC2I),DM(ADR0I),64,1)
          CALL ADD(DM(D1I),DM(ADR0I),DM(D1I),64,1)
          CALL HOMAT(MULT,DM(ADLC2I),DM(ADR0I),64,1)
          CALL ADD(DM(D1I+64),DM(ADR0I),DM(D1I+64),64,1)
          ADLC2I=ADLC2I+64
        ENDDO
      ENDIF
C 
CD    IF(NBINT.GT.0)THEN
C 
CD      CALL IMPTDN('VALEUR DU TABLEAU H*H',DM(D1I),8,8)
CD      CALL IMPTDN('VALEUR DU TABLEAU x*H*H',DM(D1I+64),8,8)
C 
CD    END IF
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine calcule NN/R integre sur l'element

C     On envoie comme arguments :
C 
C     E ....... A la demi-longueur de l'element
C     E ....... RC distance du centre de l'element au centre du trou
C 
C     Et on recupere :
C 
C     S ....... NN/R integre

      SUBROUTINE CANNSR (A, RC, NNSURR)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER D0, L0, LONGTA, NTABL, NTERME, K, H
      INTEGER I, J, NTERAV, N0, ADR0
      INTEGER AM2LC, ADM2LC
      DOUBLE PRECISION XSI, WI, WJ, W, MULT, RI
C 
      DOUBLE PRECISION A,RC,NNSURR(144)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CANNSR')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ENPOUB(AM2LC,ADM2LC)
      CALL NUTBDM('MAT-CONST ',N0)
C 
C     D0 EST LA PREMIERE ADRESSE DE MAT-CONST
C 
      D0=ADM(N0)
      L0=LONGDM(N0)
C 
CD    CALL IMPEP('LONGUEUR DE MAT-CONST', L0)
C 
C     LONGTA EST LE NOMBRE DE TERMES STOCKES POUR LE CALCUL D'UN TABLEAU
C 
      LONGTA=XINTEG*YINTEG*144
C 
C     NTABL EST LE NOMBRE DE TABLEAUX STOCKES PAR POINT DE GAUSS
C 
      NTABL=L0/LONGTA
C 
CD    CALL IMPEP('NOMBRE DE TABLEAUX STOCKES PAR POINT DE GAUSS',NTABL)
C 
C     NTERME EST LE NOMBRE DE TERMES STOCKES PAR POINT DE GAUSS
C 
      NTERME=144*NTABL
C 
C -----------------------------------------------------------------------
C     creation d'un tableau provisoire de 1ere adresse ADR0 pour 
C     ranger le tableau des termes multiplies par les RI*W
C -----------------------------------------------------------------------
      CALL GSPOUD(144, ADR0)
C 
      K=XINTEG*(XINTEG-1)/2
      H=YINTEG*(YINTEG-1)/2
C 
      DO I=1,XINTEG
        XSI=GAUSS(K+I)
        WI=POIDS(K+I)
        RI=1.D0/(A*XSI+RC)
C 
CD      CALL IMPDP('VALEUR DE RI',RI)
C 
        DO J=1,YINTEG
          WJ=POIDS(H+J)
          W=WI*WJ
C 
CD        CALL IMPDP('VALEUR DE W',W)
C 
          MULT=W*RI
C 
CD        CALL IMPDP('VALEUR DE MULT',MULT)
C 
          NTERAV=((I-1)*YINTEG+J-1)*NTERME
C 
C     TCD CALL MUMARE(MULT,144,DM(D0+NTERAV),DM(ADR0))
C     TCD CALL ADDMAD(144,NNSURR,DM(ADR0),NNSURR)
C 
          CALL HOMAT(MULT,DM(D0+NTERAV),DM(ADR0),144,1)
          CALL ADD(NNSURR,DM(ADR0),NNSURR,144,1)
        ENDDO
      ENDDO
C 
CD    CALL IMPTDN('VALEUR DU TABLEAU NNSURR INTEGRE',NNSURR,12,12)
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme argument :
C 
C     E ...... A    demi-largeur de l'element
C     E ...... B    demi-hauteur de l'element
C     E ...... RC   rayon de l'element
C     E ...... ISO  comportement isotrope transverse de la couche
C     E ...... TI   adresse de depart des tableaux integres
C     E ...... NR   adresse de depart du tableau NN/R
C 
C     Et on recupere :
C 
C     S ...... TAB1, TABN les deux tableaux utiles a l'assemblage
C 
      SUBROUTINE CALK11 (A, B, RC,I SO, TI, NR, TAB1, TABN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TI,NR,I,ZI,K,L,AD0,AP0,ADR0,APR0,ADR1
      INTEGER     AM2LC      , ADM2LC
C 
      DOUBLE PRECISION  ISO(6),TAB1(144),TABN(144),A,B,RC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALK11')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ad0 pour 
C     ranger les adresses de depart des tableaux integres qui nous interessent
C -----------------------------------------------------------------------
C 
CD    CALL IMPDN('VALEUR DE A DANS CALK11     ',A     )
CD    CALL IMPDN('VALEUR DE K11ISO DANS CALK11',ISO(1))
CD    CALL IMPDN('VALEUR DE RC  DANS CALK11   ',RC    )
CD    CALL IMPDN('VALEUR DE B  DANS CALK11    ',B     )
C  
      CALL GSPOUE(7, AD0)
      AP0   = AD0-1
      K     = AD0
      M(K)  = TI+5*144
      K     = K+1
      M(K)  = TI+10*144
      K     = K+1
      M(K)  = NR
      K     = K+1
      M(K)  = TI+9*144
      K     = K+1
      M(K)  = TI+8*144
      K     = K+1
      M(K)  = TI+13*144
      K     = K+1
      M(K)  = NR
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR0 pour 
C     ranger les constantes issues du comportement qui nous interessent
C -----------------------------------------------------------------------
      CALL GSPOUD(7*145,ADR0)
      ADR1     = ADR0+7
      APR0     = ADR0-1
      L        = ADR0
      DM(L)    = ISO(1)*B*RC/A
      DM(L+1)  = ISO(1)*B
      DM(L+2)  = ISO(1)*B*A
      DM(L+3)  = ISO(2)*B
      DM(L+4)  = ISO(5)*A*RC/(2.D0*B)
      DM(L+5)  = ISO(5)*A*A/(2.D0*B)
      DM(L+6)  = ISO(3)*A*B/2.D0
C 
CD    CALL IMPDN('VALEUR DE (K11+B/2)BRC/A   ',DM(L)  )
CD    CALL IMPDN('VALEUR DE (K11+B/2)B       ',DM(L+1))
CD    CALL IMPDN('VALEUR DE K11ISO*A*B       ',DM(L+2))
CD    CALL IMPDN('VALEUR DE K12ISO*B         ',DM(L+3))
CD    CALL IMPDN('VALEUR DE BISO*ARC/2B      ',DM(L+4))
CD    CALL IMPDN('VALEUR DE BISO*AA/2B       ',DM(L+5))
CD    CALL IMPDN('VALEUR DE K66ISO*A*B/2     ',DM(L+6))
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR1 pour
C      ranger les ... issues du comportement  qui nous interessent
C -----------------------------------------------------------------------
      DO I=1,7
C 
C     TCD CALL MUMARE(DM(APR0+I),144,DM(M(AP0+I)),DM(ADR1+144*(I-1)))
C 
        CALL HOMAT (DM(APR0+I), DM(M(AP0+I)),
     &	            DM(ADR1+144*(I-1)), 144, 1)
      END DO
      DO I=1,5
        ZI=ADR1+144*I
C 
C     TCD CALL ADDMAD(144,DM(ZI-144),DM(ZI),DM(ZI))
C 
        CALL ADD(DM(ZI-144),DM(ZI),DM(ZI),144,1)
      ENDDO
C 
C     RANGEMENT DANS LE PREMIER TABLEAU DE SORTIE DE LA PARTIE
C     INDEPENDANTE DE N DE K11
C 
      CALL MATLOC (12, A, DM(ADR1+5*144), TAB1)
C 
C     RANGEMENT DANS LE DEUXIEME TABLEAU DE SORTIE DE LA PARTIE
C     DEPENDANT DE N DE K11 
C 
      CALL MATLOC (12, A, DM(ADR1+6*144), TABN)
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD    CALL IMPTDN('VALEUR DE TAB1'//IDPROG,TAB1,12,12)
C 
C     REMISE DES ADRESSES DES TABLEAUX PARTIELS A ZERO
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 

      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme argument :
C 
C     E ...... A    demi-largeur de l'element
C     E ...... B    demi-hauteur de l'element
C     E ...... RC   rayon de l'element
C     E ...... ISO  comportement isotrope transverse de la couche
C     E ...... TI   adresse de depart des tableaux integres
C     E ...... NR   adresse de depart du tableau NN/R
C 
C     Et on recupere :
C 
C     S ...... TABK12 le tableau utile a l'assemblage
C 
         SUBROUTINE CALK12 (A, B, RC, ISO, TI, NR, TABK12)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  ISO(3),TABK12(144),A,B,RC
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TI,NR,I,ZI,K,L,AD0,AP0,ADR0,APR0,ADR1
C 
      INTEGER     AM2LC      , ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALK12')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ad0 pour ranger
C     les adresses de depart des tableaux integres qui nous interessent
C -----------------------------------------------------------------------
CD    CALL IMPDN('VALEUR DE A DANS CALK12     ',A     )
CD    CALL IMPDN('VALEUR DE K11ISO DANS CALK12',ISO(1))
CD    CALL IMPDN('VALEUR DE RC  DANS CALK12   ',RC    )
CD    CALL IMPDN('VALEUR DE B  DANS CALK12    ',B     )
C 
      CALL GSPOUE(3,AD0)
      AP0       = AD0-1
      K         = AD0
      M(K)      = TI+2*144
      K         = K+1
      M(K)      = TI+144
      K         = K+1
      M(K)      = NR
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR0 pour ranger
C     les constantes issues du comportement qui nous interessent
C -----------------------------------------------------------------------
      CALL GSPOUD(3*145,ADR0)
      ADR1      = ADR0+3
      APR0      = ADR0-1
      L         = ADR0
      DM(L)     = ISO(2)*B
      DM(L+1)   = -ISO(3)*B/2.D0
      DM(L+2)   = (ISO(1)+ISO(3)/2.D0)*A*B
C 
CD    CALL IMPDN('VALEUR DE K12ISO*B      ',DM(L))
CD    CALL IMPDN('VALEUR DE -K66ISO*B/2   ',DM(L+1))
CD    CALL IMPDN('VALEUR DE (K11+K66/2)AB ',DM(L+2)  )
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR1 pour ranger
C     les ... issu du comportement qui nous qui nous interesse
C -----------------------------------------------------------------------
      DO I=1,3
C 
C     TCD CALL MUMARE(DM(APR0+I),144,DM(M(AP0+I)),DM(ADR1+144*(I-1)))
C 
        CALL HOMAT(DM(APR0+I),DM(M(AP0+I)),DM(ADR1+144*(I-1)),144,1)
      ENDDO
      DO I=1,2
        ZI=ADR1+144*I
C 
C     TCD CALL ADDMAD(144,DM(ZI-144),DM(ZI),DM(ZI))
C 
        CALL ADD(DM(ZI-144),DM(ZI),DM(ZI),144,1)
      ENDDO
C 
C     RANGEMENT DANS LE TABLEAU DE SORTIE DE TABK12
C     LE TABLEAU K12 SERA A MULTIPLIER PAR N DANS ASSBUW
C 
      CALL MATLOC(12,A,DM(ADR1+2*144),TABK12)
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD    CALL IMPTDN('VALEUR DE TABK12'//IDPROG,TABK12,12,12)
C 
C     REMISE DES ADRESSES DES TABLEAUX PARTIELS A ZERO
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme argument :
C 
C     E ...... A    demi-largeur de l'element
C     E ...... B    demi-hauteur de l'element
C     E ...... RC   rayon de l'element
C     E ...... ISO  comportement isotrope transverse de la couche
C     E ...... TI   adresse de depart des tableaux integres
C     E ...... NR   adresse de depart du tableau NN/R
C 
C     Et on recupere :
C 
C     S ...... TABK23 le tableau utile a l'assemblage
C 
      SUBROUTINE CALK23 (A, B, RC, ISO, TI, NR, TABK23)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TI,NR,I,ZI,K,L,AD0,AP0,ADR0,APR0,ADR1
      INTEGER     AM2LC      , ADM2LC
C 
      DOUBLE PRECISION  ISO(6),TABK23(144),A,B,RC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALK23')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ad0 pour ranger
C     les adresses de depart des tableaux integres qui nous interessent
C -----------------------------------------------------------------------
      CALL GSPOUE(2,AD0)
      AP0       = AD0-1
      K         = AD0
      M(K)      = TI+3*144
      K         = K+1
      M(K)      = TI+4*144
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR0 pour ranger
C     les constantes issues du comportement qui nous interessent
C -----------------------------------------------------------------------
      CALL GSPOUD(2*145, ADR0)
      ADR1      = ADR0+2
      APR0      = ADR0-1
      L         = ADR0
      DM(L)     = ISO(4)*A
      DM(L+1)   = -ISO(5)*A/2.D0
C 
CD    CALL IMPDN('VALEUR DE AISO*A  ',DM(L)  )
CD    CALL IMPDN('VALEUR DE -BISO*A/2',DM(L+1))
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR1 pour ranger
C     les ... issu du comportement qui nous interesse
C -----------------------------------------------------------------------
      DO I=1,2
C 
C     TCD CALL MUMARE(DM(APR0+I),144,DM(M(AP0+I)),DM(ADR1+144*(I-1)))
C 
        CALL HOMAT(DM(APR0+I),DM(M(AP0+I)),DM(ADR1+144*(I-1)),144,1)
      ENDDO
      ZI=ADR1+144
C 
C     TCD CALL ADDMAD(144,DM(ZI-144),DM(ZI),DM(ZI))
C 
      CALL ADD(DM(ZI-144),DM(ZI),DM(ZI),144,1)
C 
C     RANGEMENT DANS LE TABLEAU DE SORTIE DE TABK23
C     LE TABLEAU TABK23 SERA A MULTIPLIER PAR N DANS ASSBVW
C 
      CALL MATLOC(12,A,DM(ADR1+144),TABK23)
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD      CALL IMPTDN('VALEUR DE TABK23'//IDPROG,TABK23,12,12)
C 
C     REMISE DES ADRESSES DES TABLEAUX PARTIELS A ZERO
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     01 /  3  / 87
C 
C     On envoie comme argument :
C 
C     E ...... A    demi-largeur de l'element
C     E ...... B    demi-hauteur de l'element
C     E ...... RC   rayon de l'element
C     E ...... ISO  comportement isotrope transverse de la couche
C     E ...... TI   adresse de depart des tableaux integres
C     E ...... NR   adresse de depart du tableau NN/R
C 
C     Et on recupere :
C 
C     S ...... TAB13 le tableau utile a l'assemblage
C 
      SUBROUTINE CALK13 (A, B, RC, ISO, TI, NR, TAB13)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TI,NR,I,ZI,K,L,AD0,AP0,ADR0,APR0,ADR1
C 
      INTEGER     AM2LC      , ADM2LC
      DOUBLE PRECISION  ISO(5),TAB13(144),A,B,RC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALK13')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ad0 pour ranger
C     les adresses de depart des tableaux integres qui nous interessent
C -----------------------------------------------------------------------
CD    CALL IMPDN('VALEUR DE A DANS CALK13   ',A     )
CD    CALL IMPDN('VALEUR DE AISO DANS CALK13',ISO(4))
CD    CALL IMPDN('VALEUR DE BISO DANS CALK13',ISO(5))
CD    CALL IMPDN('VALEUR DE RC  DANS CALK13 ',RC    )
CD    CALL IMPDN('VALEUR DE B  DANS CALK13  ',B     )
C 
      CALL GSPOUD(5*145, ADR0)
      ADR1      = ADR0+5
      APR0      = ADR0-1
      L         = ADR0
      DM(L)     = ISO(4)*RC
      DM(L+1)   = ISO(4)*A
      DM(L+2)   = ISO(4)*A
      DM(L+3)   = ISO(5)*RC/2.D0
      DM(L+4)   = ISO(5)*A/2.D0
C 
CD    CALL IMPDN('VALEUR DE (A)*RC  ',DM(L)  )
CD    CALL IMPDN('VALEUR DE (A2)*A   ',DM(L+1))
CD    CALL IMPDN('VALEUR DE AISO*A      ',DM(L+2))
C 
      CALL GSPOUE(5,AD0)
      AP0     = AD0-1
      K       = AD0
      M(K)    = TI+6*144
      K       = K+1
      M(K)    = TI+11*144
      K       = K+1
      M(K)    = TI+3*144
      K       = K+1
      M(K)    = TI+7*144
      K       = K+1
      M(K)    = TI+12*144
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR0 pour ranger
C     les constantes issues du comportement qui nous interessent
C -----------------------------------------------------------------------
      DO I=1,5
C 
C     TCD CALL MUMARE(DM(APR0+I),144,DM(M(AP0+I)),DM(ADR1+144*(I-1)))
C 
        CALL HOMAT(DM(APR0+I),DM(M(AP0+I)),DM(ADR1+144*(I-1)),144,1)
      ENDDO
      DO I=1,4
        ZI=ADR1+144*I
C 
C     TCD CALL ADDMAD(144,DM(ZI-144),DM(ZI),DM(ZI))
C 
        CALL ADD(DM(ZI-144),DM(ZI),DM(ZI),144,1)
      ENDDO
C 
C     RANGEMENT DANS LE TABLEAU DE SORTIE DE TABK13
C 
      CALL MATLOC(12,A,DM(ADR1+4*144),TAB13)
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD    CALL IMPTDN('VALEUR DE TAB13'//IDPROG,TAB13,12,12)
C 
C     REMISE DES ADRESSES DES TABLEAUX PARTIELS A ZERO
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     5 /  3  / 87
C 
C     On envoie comme argument :
C 
C     E ...... A    demi-largeur de l'element
C     E ...... B    hauteur de l'element
C     E ...... RC   rayon de l'element
C     E ...... ISO  comportement isotrope transverse de la couche
C     E ...... TI   adresse de depart des tableaux integres
C     E ...... NR   adresse de depart du tableau NN/R
C 
C     Et on recupere :
C 
C     S ...... TAB1, TABN les deux tableaux utiles a l'assemblage
C 
      SUBROUTINE CALK22 (A, B, RC, ISO, TI, NR, TAB1, TABN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TI,NR,I,ZI,K,L,AD0,AP0,ADR0,APR0,ADR1
      INTEGER     AM2LC      , ADM2LC
C 
      DOUBLE PRECISION  ISO(6),TAB1(144),TABN(144),A,B,RC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALK22')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ad0 pour ranger
C     les adresses de depart des tableaux integres qui nous interessent
C -----------------------------------------------------------------------
CD    CALL IMPDN('VALEUR DE A DANS CALK22    ',A     )
CD    CALL IMPDN('VALEUR DE K11ISO DANS CALK2',ISO(1))
CD    CALL IMPDN('VALEUR DE RC  DANS CALK22  ',RC    )
CD    CALL IMPDN('VALEUR DE B  DANS CALK22   ',B     )
C 
      CALL GSPOUE(7, AD0)
      AP0       = AD0-1
      K         = AD0
      M(K)      = NR
      K         = K+1
      M(K)      = TI+5*144
      K         = K+1
      M(K)      = TI+10*144
      K         = K+1
      M(K)      = TI+9*144
      K         = K+1
      M(K)      = TI+8*144
      K         = K+1
      M(K)      = TI+13*144
      K         = K+1
      M(K)      = NR
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere ADR0 adresse pour ranger
C     les constantes issues du comportement qui nous interessent
C -----------------------------------------------------------------------
      CALL GSPOUD(7*145,ADR0)
      ADR1       = ADR0+7
      APR0       = ADR0-1
      L          = ADR0
      DM(L)      = ISO(3)*A*B/2.D0
      DM(L+1)    = ISO(3)*RC*B/(2.D0*A)
      DM(L+2)    = ISO(3)*B/2.D0
      DM(L+3)    = -ISO(3)*B/2.D0
      DM(L+4)    = ISO(5)*A*RC/(2.D0*B)
      DM(L+5)    = ISO(5)*A*A/(2.D0*B)
      DM(L+6)    = ISO(1)*A*B
C 
CD    CALL IMPDN('VALEUR DE K66*A*B/2      ',DM(L)    )
CD    CALL IMPDN('VALEUR DE K66*RC*B/2A    ',DM(L+1)  )
CD    CALL IMPDN('VALEUR DE K66IS0*B/2     ',DM(L+2)  )
CD    CALL IMPDN('VALEUR DE -K66IS0*B/2    ',DM(L+3)  )
CD    CALL IMPDN('VALEUR DE K11AB          ',DM(L+6)  )
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere ADR1 pour ranger
C     les ... issues du comportement qui nous interessent
C -----------------------------------------------------------------------
      DO I=1,7
C 
C     TCD CALL MUMARE(DM(APR0+I),144,DM(M(AP0+I)),DM(ADR1+144*(I-1)))
C 
        CALL HOMAT(DM(APR0+I),DM(M(AP0+I)),DM(ADR1+144*(I-1)),144,1)
      ENDDO
      DO I=1,5
        ZI=ADR1+144*I
C 
C     TCD CALL ADDMAD(144,DM(ZI-144),DM(ZI),DM(ZI))
C 
        CALL ADD(DM(ZI-144),DM(ZI),DM(ZI),144,1)
      ENDDO
C 
C     RANGEMENT DANS LE PREMIER TABLEAU DE SORTIE DE LA PARTIE
C     INDEPENDANTE DE N DE K22
C 
      CALL MATLOC(12,A,DM(ADR1+5*144),TAB1)
C 
C     RANGEMENT DANS LE DEUXIEME TABLEAU DE SORTIE DE LA PARTIE
C     DEPENDANTE DE N K22
C 
      CALL MATLOC(12,A,DM(ADR1+6*144),TABN)
C 
C     SEQUENCE D'IMPRESION EVENTUELLE
C 
CD    CALL IMPTDN('VALEUR DE TAB1'//IDPROG,TAB1,12,12)
CD    CALL IMPTDN('VALEUR DE TABN'//IDPROG,TABN,12,12)
C 
C    REMISE DES ADRESSES DES TABLEAUX PARTIELS A ZERO
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 

      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     5 /  3  / 87
C 
C     On envoie comme argument :
C 
C     E ...... A    demi-largeur de l'element
C     E ...... B    hauteur de l'element
C     E ...... RC   rayon de l'element
C     E ...... ISO  comportement isotrope transverse de la couche
C     E ...... TI   adresse de depart des tableaux integres
C     E ...... NR   adresse de depart du tableau NN/R
C 
C     Et on recupere :
C 
C     S ...... TAB1, TABN les deux tableaux utiles a l'assemblage
C 
      SUBROUTINE CALK33 (A, B, RC, ISO, TI, NR, TAB1, TABN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TI,NR,I,ZI,K,L,AD0,AP0,ADR0,APR0,ADR1
      INTEGER     AM2LC      , ADM2LC
C 
      DOUBLE PRECISION  ISO(6),TAB1(144),TABN(144),A,B,RC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALK33')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ad0 pour ranger
C     les adresses de depart des tableaux integres qui nous interessent
C -----------------------------------------------------------------------
CD    CALL IMPDN('VALEUR DE A DANS CALK33     ',A     )
CD    CALL IMPDN('VALEUR DE K11ISO DANS CALK2 ',ISO(1))
CD    CALL IMPDN('VALEUR DE RC  DANS CALK33   ',RC    )
CD    CALL IMPDN('VALEUR DE B  DANS CALK33    ',B     )
C 
      CALL GSPOUE(5, AD0)
      AP0       = AD0-1
      K         = AD0
      M(K)      = TI+8*144
      K         = K+1
      M(K)      = TI+13*144
      K         = K+1
      M(K)      = TI+5*144
      K         = K+1
      M(K)      = TI+10*144
      K         = K+1
      M(K)      = NR
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR0 pour ranger
C     les constantes issues du comportement qui nous interessent
C -----------------------------------------------------------------------
      CALL GSPOUD(5*145,ADR0)
      ADR1       = ADR0+5
      APR0       = ADR0-1
      L          = ADR0
      DM(L)      = ISO(6)*A*RC/B
      DM(L+1)    = ISO(6)*A*A/B
      DM(L+2)    = ISO(5)*B*RC/(2.D0*A)
      DM(L+3)    = ISO(5)*B/2.D0
      DM(L+4)    = ISO(5)*A*B/2.D0
C 
CD    CALL IMPDN('VALEUR DE (C+B/2)*ARC/B    ',DM(L)  )
CD    CALL IMPDN('VALEUR DE (C+B/2)*AA/B     ',DM(L+1))
CD    CALL IMPDN('VALEUR DE BISO*AB/2        ',DM(L+2))
C 
C -----------------------------------------------------------------------
C     Creation d'un tableau provisoire de premiere adresse ADR1 pour ranger
C     les ... issues du comportement qui nous qui nous interessent
C -----------------------------------------------------------------------
      DO I=1,5
C 
C     TCD CALL MUMARE(DM(APR0+I),144,DM(M(AP0+I)),DM(ADR1+144*(I-1)))
C 
        CALL HOMAT(DM(APR0+I),DM(M(AP0+I)),DM(ADR1+144*(I-1)),144,1)
      ENDDO
      DO I= 1,3
        ZI=ADR1+144*I
C 
C     TCD CALL ADDMAD(144,DM(ZI-144),DM(ZI),DM(ZI))
C 
        CALL ADD(DM(ZI-144),DM(ZI),DM(ZI),144,1)
      ENDDO
C 
C     RANGEMENT DANS LE PREMIER TABLEAU DE SORTIE DE LA PARTIE
C     INDEPENDANTE DE N DE K11
C 
      CALL MATLOC(12,A,DM(ADR1+3*144),TAB1)
C 
C     RANGEMENT DANS LE DEUXIEME TABLEAU DE SORTIE DE LA PARTIE
C     DEPENDANTE DE N*N DE K33
C 
      CALL MATLOC(12,A,DM(ADR1+4*144),TABN)
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD    CALL IMPTDN('VALEUR DE TAB1'//IDPROG,TAB1,12,12)
CD    CALL IMPTDN('VALEUR DE TABN'//IDPROG,TABN,12,12)
C 
C     REMISE DES ADRESSES DES TABLEAUX PARTIELS A ZERO
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     14 /  2  / 87
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine assemble pour un element donne les termes
C     provenant de CALK11 pour toutes les matrices KOn
C 
C     On envoie comme arguments :
C 
C     E ...... N        valeur du rang de Fourier
C     E ...... TK0N     tableau des adresses de depart de la matrice KON
C     E ...... UCAL     tableau des adresses elementaires UCAL
C     E ...... TAB1     tableau du terme independant de n provenant de
C                       CALK11 pour l'element considere
C     E ...... TABN     tableau du terme dependant de n provenant de
C                       CALK11 pour l'element considere
C     E ...... NDDLU    tableau des numeros reels des ddl du deplacement
C                       U de l'element considere
C     E ...... P1       1ere adresse du tableau profil

      SUBROUTINE ASSBU (N, TK0N, UCAL, TAB1, TABN, NDDLU, P1)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER UCAL(12), NDDLU(12), P1, NI, UCALI, N, I
      INTEGER J, NJ, AD1
      INTEGER UCALJ, TK0N
C 
      DOUBLE PRECISION  TAB1(12,12), TABN(12,12), TCAL, DBLE
C 
      CHARACTER*6 IDPROG, BARATI
      PARAMETER (IDPROG='ASSBU')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
        WRITE (BARATI(1:6),'(I6)')N
C 
CD      CALL IMPTDP('VALEUR DU TABLEAU
CD      K0('//BARATI(1:6),DM(TK0N),1,NTMAT)
C 
        DO I=1,12
          NI    = NDDLU(I)
          UCALI = UCAL(I)
C 
CD        CALL IMPEP('POUR I',I)
CD        CALL IMPEP('VALEUR DE NI DANS ASSBU',NI)
CD        CALL IMPEP('VALEUR DE UCALI DANS ASSBU',UCALI)
C 
          DO J=I,12
C 
CD          CALL IMPEP('POUR J',J)
C 
            NJ       = NDDLU(J)
            UCALJ    = UCAL(J)
            AD1      = TK0N-1+NI-NJ+M(P1+NJ)
            TCAL     = TAB1(UCALI,UCALJ)+DBLE(N*N)*TABN(UCALI,UCALJ)
            DM(AD1)  = DM(AD1)+TCAL
          END DO
        END DO
C 
C -----------------------------------------------------------------------
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD      WRITE(BARATI(1:6),'(I6)')N
CD      CALL IMPTDP('VALEUR DU TABLEAU
CD      K0('//BARATI(1:6),DM(TK0N),1,NTMAT)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     14 /  2  / 87
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine assemble pour un element donne
C     les termes provenant de CALK22 pour toutes les matrices KOn
C 
C     On envoie comme arguments :
C 
C     E ...... N        valeur du rang de Fourier
C     E ...... TK0N     tableau des adresses de depart des matrices KON
C     E ...... UCAL     tableau des adresses elementaires UCAL
C     E ...... TAB1     tableau du terme independant de n provenant de
C                       CALK22 pour l'element considere
C     E ...... TABN     tableau du terme dependant de n provenant de
C                       CALK22 pour l'element considere
C     E ...... NDDLV    tableau des numeros reels des ddl du deplacement
C                       U de l'element considere
C     E ...... P1       1ere adresse du tableau profil
C 
      SUBROUTINE ASSBV (N, TK0N, UCAL, TAB1, TABN, NDDLV, P1)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER UCAL(12), NDDLV(12), P1, NI, UCALI
      INTEGER N, I, J, NJ, AD1
      INTEGER UCALJ, TK0N
C 
      DOUBLE PRECISION  TAB1(12,12), TABN(12,12), TCAL, DBLE
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ASSBV')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPTEN('RANGEMENT VCAL DANS '//IDPROG,UCAL(1),1,12)
CD    CALL IMPTEN('RANGEMENT  CROISSANT DANS '//IDPROG,NDDLV(1),1,12)
C 
C    X CALL IMPEN('VALEUR DE NTSFG DANS ASSBV',NTDSFG)
C 
        DO I=1,12
          NI    = NDDLV(I)
          UCALI = UCAL(I)
C 
CD        CALL IMPEP('POUR I',I)
CD        CALL IMPEP('VALEUR DE NI DANS ASSBV',NI)
CD        CALL IMPEP('VALEUR DE UCALI DANS ASSBV',UCALI)
C 
          DO J=I,12
C 
CD          CALL IMPEP('POUR J',J)
C 
            NJ    = NDDLV(J)
            UCALJ = UCAL(J)
            AD1   = TK0N-1+NI-NJ+M(P1+NJ)
            TCAL  = TAB1(UCALI,UCALJ)+DBLE(N*N)*TABN(UCALI,UCALJ)
            DM(AD1)=DM(AD1)+TCAL
C 
CD          CALL IMPDN('VALEUR DE DM(AD1) APRES=TCAL?',DM(AD1))
C 
          ENDDO
        ENDDO
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD      WRITE (BARATI(1:6),'(I6)')N
CD      CALL IMPTDP ('VALEUR DU TABLEAU
CD      K0('//BARATI(1:6),DM(TK0N),1,NTMAT)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     14 /  2  / 87
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine assemble pour un element donne
C     les termes provenant de CALK33 pour toutes les matrices KOn
C 
C     On envoie comme arguments :
C 
C     E ...... N        valeur du rang de Fourier
C     E ...... TK0      tableau des adresses de depart des matrices KON
C     E ...... UCAL     tableau des adresses elementaires UCAL
C     E ...... TAB1     tableau du terme independant de n provenant de
C                       CALK33 pour l'element considere
C     E ...... TABN     tableau du terme dependant de n provenant de
C                       CALK33 pour l'element considere
C     E ...... NDDLW    tableau des numeros reels des ddl du deplacement
C                       W de l'element considere
C     E ...... P1       1ere adresse du tableau profil
C 
      SUBROUTINE ASSBW (N, TK0N, UCAL, TAB1, TABN, NDDLW, P1)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER UCAL(12), NDDLW(12), P1, NI, UCALI
      INTEGER N, I, J, NJ, AD1
      INTEGER UCALJ, TK0N
C 
      DOUBLE PRECISION  TAB1(12,12), TABN(12,12), TCAL, DBLE
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ASSBW')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPTEN('RANGEMENT UCAL DANS '//IDPROG,UCAL(1),1,12)
CD    CALL IMPTEN('RANGEMENT  CROISSANT DANS '//IDPROG,NDDLW(1),1,12)
C 
        DO I=1,12
          NI    = NDDLW(I)
          UCALI = UCAL(I)
C 
CD        CALL IMPEP('POUR I',I)
CD        CALL IMPEP('VALEUR DE NI DANS ASSBW',NI)
CD        CALL IMPEP('VALEUR DE UCALI DANS ASSBW',UCALI)
C 
          DO J=I,12
C 
CD          CALL IMPEP('POUR J',J)
C 
            UCALJ   = UCAL(J)
            NJ      = NDDLW(J)
C 
CD          CALL IMPEP('VALEUR DE NJ DANS ASSBW',NJ)
C 
            AD1     = TK0N-1+NI-NJ+M(P1+NJ)
C 
CD          CALL IMPEP('VALEUR DE M(P1+NJ) DANS ASSBW',M(P1+NJ))
CD          CALL IMPEP('VALEUR DE AD1 DANS ASSBW',AD1)
C 
            TCAL    = TAB1(UCALI,UCALJ)+DBLE(N*N)*TABN(UCALI,UCALJ)
            DM(AD1) = DM(AD1)+TCAL
          ENDDO
        ENDDO
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD      WRITE(BARATI(1:6),'(I6)')N
CD      CALL IMPTDP('VALEUR DU TABLEAU
CD      K0('//BARATI(1:6),DM(TK0N),1,NTMAT)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     02 /  3  / 87
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine assemble pour un element donne les termes provenant
C     de CALK12 pour toutes les matrices KOn
C 
C     On envoie comme arguments :
C 
C     E ...... N        valeur du rang de Fourier
C     E ...... TK0      tableau des adresses de depart des matrices KON
C     E ...... UCAL     tableau des adresses elementaires UCAL
C     E ...... TABK12   tableau du terme dependant de n provenant de
C                       CALK12 pour l'element considere
C     E ...... NDDLU    tableau des numeros reels des ddl du deplacement
C                       V de l'element considere
C     E ...... NDDLV    tableau des numeros reels des ddl du deplacement
C                       V de l'element considere
C     E ...... P1       1ere adresse du tableau profil

         SUBROUTINE ASSBUV (N, TK0N, UCAL, TABK12, NDDLU, NDDLV, P1)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER UCAL(12), NDDLU(12), NDDLV(12)
      INTEGER N, I, J, NJ, AD1, P1, NI, TK0N, MOD
      INTEGER ADI1, ADJ1, UCALI, UCALJ
C 
      DOUBLE PRECISION  TABK12(12,12), DBLE
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ASSBUV')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPEN('VALEUR DE NTSFG DANS ASSBUV',NTDSFG)
CD    CALL IMPTEN('RANGEMENT UCAL DANS '//IDPROG,UCAL(1),1,12)
CD    CALL IMPTEN('RANG U  CROISSANT DANS '//IDPROG,NDDLU(1),1,12)
CD    CALL IMPTEN('RANG V  CROISSANT DANS '//IDPROG,NDDLV(1),1,12)
C 
        DO J = 1, 12
          NJ       = NDDLV(J)
          UCALJ    = UCAL(J)
          ADJ1     = -NJ+M(P1+NJ)-1+TK0N
C 
CD        CALL IMPEP('VALEUR DE ADJ1 DANS ASSBUV',ADJ1)
CD        CALL IMPEP('POUR J',J)
CD        CALL IMPEP('VALEUR DE NJ DANS ASSBUV',NJ)
C 
          DO I = 1, J+MOD(J,2)
C 
CD          CALL IMPEP('POUR I',I)
C 
            NI     =NDDLU(I)
C 
CD          CALL IMPEP('VALEUR DE NI DANS ASSBUV',NI)
C 
            AD1    =NI+ADJ1
            DM(AD1)=DM(AD1)+DBLE(N)*TABK12(UCAL(I),UCALJ)
          ENDDO
        ENDDO
        DO I = 3,12
C 
CD        CALL IMPEP('POUR I',I)
C 
          NI=NDDLU(I)
C 
CD        CALL IMPEP('VALEUR DE NI DANS ASSBUV',NI)
C 
          ADI1     = -NI+M(P1+NI)-1+TK0N
          UCALI    = UCAL(I)
          DO J = 1, I-1-MOD(I-1,2)
            NJ       =  NDDLV(J)
            AD1      =  NJ+ADI1
            UCALJ    =  UCAL(J)
            DM(AD1)  =DM(AD1)+DBLE(N)*TABK12(UCALI,UCALJ)
C 
CD          CALL IMPDN('VALEUR DE DM(AD1) APRES=TCAL?',DM(AD1))
C 
          ENDDO
        ENDDO
C 
C     SEQUENCE D'IMPRESSIOM EVENTUELLE
C 
CD      WRITE (BARATI(1:6),'(I6)')N
CD      CALL IMPTDP ('VALEUR DU TABLEAU
CD      K0('//BARATI(1:6),DM(TK0N),1,NTMAT)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     02 /  3  / 87
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine assemble pour un element donne
C     les termes provenant de CALK23 pour toutes les matrices KOn
C 
C     On envoie comme arguments :
C 
C     E ...... N        valeur du rang de Fourier
C     E ...... TK0      tableau des adresses de depart des matrices KON
C     E ...... UCAL     tableau des adresses  elementaires UCAL
C     E ...... TABK23   tableau du terme dependant de n provenant de
C                       CALK12 pour l'element considere
C     E ...... NDDLV    tableau des numeros reels des ddl du deplacement
C                       V de l'element considere
C     E ...... NDDLW    tableau des numeros reels des ddl du deplacement
C                       W de l'element considere
C     E ...... P1       1ere adresse du tableau profil

         SUBROUTINE ASSBVW (N, TK0N, UCAL, TABK23, NDDLV, NDDLW, P1)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER UCAL(12), NDDLV(12), NDDLW(12)
      INTEGER N, I, J, NJ, AD1, P1, NI, MOD
      INTEGER ADI1, ADJ1, TK0N, UCALI, UCALJ
C 
      DOUBLE PRECISION  TABK23(12,12), DBLE
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ASSBVW')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPTEN('RANGEMENT UCAL DANS '//IDPROG,UCAL(1),1,12)
CD    CALL IMPTEN('RANGEMENT  CROISSANT DE V DANS'
CD                //IDPROG,NDDLV(1),1,12)
CD    CALL IMPTEN('RANG DE W  CROISSANT DANS '//IDPROG,NDDLW(1),1,12)
CD    CALL IMPEN('VALEUR DE NTSFG DANS ASSBVW',NTDSFG)
C 
        DO J=1,12
          NJ       = NDDLW(J)
          ADJ1     = -NJ+M(P1+NJ)-1+TK0N
C 
CD        CALL IMPEP('POUR J',J)
CD        CALL IMPEP('VALEUR DE NJ DANS ASSBVW',NJ)
C 
          UCALJ    = UCAL(J)
          DO I = 1,J+MOD(J,2)
C 
CD          CALL IMPEP('POUR I',I)
C 
            NI     =NDDLV(I)
C 
CD          CALL IMPEP('VALEUR DE NI DANS ASSBVW',NI)
C 
            AD1    =NI+ADJ1
C 
CD          CALL IMPEP('POUR AD1',AD1)
CD          CALL IMPDN('VALEUR DE DM(AD1) AVANT',DM(AD1))
C 
            DM(AD1)=DM(AD1)+DBLE(N)*TABK23(UCAL(I),UCALJ)
C 
CD          CALL IMPDN('VALEUR DE DM(AD1) APRES=TCAL?',DM(AD1))
C 
          ENDDO
        ENDDO
        DO I  = 3, 12
C 
CD        CALL IMPEP('POUR I',I)
C 
          UCALI    = UCAL(I)
          NI       = NDDLV(I)
          DO J = 1, I-1-MOD(I-1,2)
C 
CD          CALL IMPEP('VALEUR DE NI DANS ASSBUV',NI)
C 
            NJ       =  NDDLW(J)
            ADI1     = -NI+M(P1+NI)-1+TK0N
            AD1      =  NJ+ADI1
C 
CD          CALL IMPEP('VALEUR DE AD1 DANS ASSBUV',AD1)
CD          CALL IMPDN('VALEUR DE DM(AD1) AVANT',DM(AD1))
C 
            DM(AD1)=DM(AD1)+DBLE(N)*TABK23(UCALI,UCAL(J))
C 
CD          CALL IMPDN('VALEUR DE DM(AD1) APRES=TCAL?',DM(AD1))
C 
          ENDDO
        ENDDO
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD      WRITE(BARATI(1:6),'(I6)')N
CD      CALL IMPTDP(
CD      'K0('//BARATI(1:6),DM(TK0N),1,NTMAT)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     02 /  3  / 87
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine assemble pour un element donne
C     les termes provenant de CALK13 pour toutes les matrices KOn
C 
C     On envoie comme arguments :
C 
C     E ...... N        valeur du rang de Fourier
C     E ...... TK0N     tableau des adresses de depart des matrices KON
C     E ...... UCAL     tableau des adresses  elementaires UCAL
C     E ...... TAB13    tableau du terme independant de n provenant de
C                       CALK13 pour l'element considere
C     E ...... NDDLU    tableau des numeros reels des ddl du deplacement
C                       V de l'element considere
C     E ...... NDDLW    tableau des numeros reels des ddl du deplacement
C                       V de l'element considere
C     E ...... P1       1ere adresse du tableau profil

      SUBROUTINE ASSBUW (N, TK0N, UCAL, TAB13, NDDLU, NDDLW, P1)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER UCAL(12), NDDLU(12), NDDLW(12)
      INTEGER N, I, J, NJ, AD1, P1, NI, MOD, TK0N
      INTEGER ADI1, ADJ1, UCALI, UCALJ
C 
      DOUBLE PRECISION  TAB13(12,12)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ASSBUW')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    CALL IMPEN('VALEUR DE NTSFG DANS ASSBUW',NTDSFG)
CD    CALL IMPTEN('RANGEMENT UCAL DANS '//IDPROG,UCAL(1),1,12)
CD    CALL IMPTEN('RANG DE U  CROISSANT DANS '//IDPROG,NDDLU(1),1,12)
CD    CALL IMPTEN('RANG DE W  CROISSANT DANS '//IDPROG,NDDLW(1),1,12)
C 
        DO J=1,12
          NJ       = NDDLW(J)
          ADJ1     = -NJ+M(P1+NJ)-1+TK0N
          UCALJ    = UCAL(J)
C 
CD          CALL IMPEP('VALEUR DE ADJ1 DANS ASSBUW',ADJ1)
CD          CALL IMPEP('POUR J',J)
CD          CALL IMPEP('VALEUR DE NJ DANS ASSBUW',NJ)
C 
          DO I = 1, J+MOD(J,2)
C 
CD          CALL IMPEP('POUR I',I)
C 
            NI     = NDDLU(I)
            UCALI  = UCAL(I)
C 
CD          CALL IMPEP('VALEUR DE NI DANS ASSBUW',NI)
C 
            AD1    =NI+ADJ1
C 
CD            CALL IMPEP('POUR AD1',AD1)
CD            CALL IMPDN('VALEUR DE DM(AD1) AVANT',DM(AD1))
C 
              DM(AD1)=DM(AD1)+TAB13(UCALI,UCALJ)
C 
CD            CALL IMPDN('VALEUR DE DM(AD1) APRES=TCAL?',DM(AD1))
C 
          ENDDO
        ENDDO
        DO I =3,12
          NI       = NDDLU(I)
          UCALI    = UCAL(I)
          DO J= 1 , I-1-MOD(I-1,2)
            NJ         = NDDLW(J)
            UCALJ      = UCAL(J)
C 
CD            CALL IMPEP('POUR I',I)
CD            CALL IMPEP('VALEUR DE NI DANS ASSBUW',NI)
C 
            ADI1     = -NI+M(P1+NI)-1+TK0N
            AD1      =  NJ+ADI1
C 
CD            CALL IMPEP('VALEUR DE AD1 DANS ASSBUW',AD1)
CD            CALL IMPDN('VALEUR DE DM(AD1) AVANT',DM(AD1))
C 
            DM(AD1)=DM(AD1)+TAB13(UCALI,UCALJ)
C 
CD            CALL IMPDN('VALEUR DE DM(AD1) APRES=TCAL?',DM(AD1))
C 
          ENDDO
        ENDDO
C 
C     SEQUENCE D'IMPRESSION EVENTUELLE
C 
CD        WRITE(BARATI(1:6),'(I6)')N
CD        CALL IMPTDP('VALEUR DU TABLEAU
CD      K0('//BARATI(1:6),DM(TK0N),1,NTMAT)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     CECI DOIT SE TROUVER DANS CALK0
C 
C     VERSION DU     17 /  2  / 87
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine calcule les tableaux d'integration pour le
C     calcul des termes en 1/r; on ouvre le tableau MATCOL
C     de dimension nbcol*144 pour ranger le resultat des
C     integrations sur les dfferents elements
 
      SUBROUTINE TABCOL
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER AD0,AC0,I,ACRAY,ACLONG
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TABCOL')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
C     Creation du tableau des valeurs integrees sur les colonnes
C 
C             nom : MAT-COL
C             1ere adresse libre: AD0
C 
C -----------------------------------------------------------------------
      CALL GESTDP('MATCOL    ',144*NBCOL,AD0)
      AC0=AD0-1
      CALL MENADM(AD0,144*NBCOL)
C 
C    CALL IMPEP('VALEUR DU NB DE TABLEAUX DANS '//IDPROG,NBTADM)
C 
      CALL ADTBDM('LONGUEURS ',ACLONG)
      ACLONG=ACLONG-1
      CALL ADTBDM('RAYON     ',ACRAY)
      ACRAY=ACRAY-1
C 
      DO I=1,NBCOL
        CALL CANNSR(DM(ACLONG+I),DM(ACRAY+I),DM(AD0+(I-1)*144))
      ENDDO
C 
CD    CALL IMPTDP('TABLEAU MATCOL DES VALEURS INTEGREES PAR
CD                 COLONNES',DM(AD0),12,12*NBCOL)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     QUE FAIT CETTE ROUTINE :
C 
C     On envoie comme arguments :
C     Et on recupere :

C     VERSION DU     16 /  2  / 87
C 
      SUBROUTINE CALK0
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER ADLONG, ADRAYO, ADINTE, NUMINT
      INTEGER ADMATC, ADCOMP, MNBCOU, MNUCOU
      INTEGER ADPROD, TYPDEP, NDDLU, ADTAB1, I, J, K, NCOU
      INTEGER ADRANU, ADTLOC, NDDLV, NDDLW, N, ADRUI, ADRUIP
      INTEGER AINTGI, MNBINT, MNUINT, ADTAB3
      INTEGER K11, K11N, K12, K22, K22N, K13, K23, K33, K33N
      INTEGER K0N, NPIV, NDDLUI, NDDLVI, NDDLWI, ADCOMI
C  
      DOUBLE PRECISION EPAIS, VALSUP, EPSILO, EPS, GROS
C 
C     POUR LA SYMETRIE
C 
      LOGICAL  PASSP
C 
      INTEGER  ADMASO, NMAT, IN
C 
C     POUR PASSER AU DIRECT
C 
      INTEGER  NUENRS
C  
CD    LOGICAL  LTRACP, LTRACN
C  
      INTEGER  IUNIT, APRECI
      INTEGER  AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      CHARACTER*9 NOMDIR
      PARAMETER   (IDPROG='CALK0')
      PARAMETER   (GROS=1.D15)
      PARAMETER   (EPSILO=1.D -10)
C 
      CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
      CALL IMPET ('ADM2 EN ENTREE DE '//IDPROG, ADM2)
C 
C -----------------------------------------------------------------------
C 
C     REMPLISSAGE DU TABLEAU ADRESSES
C     REMPLISSAGE DU TABLEAU MATCOL
C 
      CALL TABCOL
C 
C     REMPLISSAGE DU TABLEAU TAB-ING
C     REMPLISSAGE DES TABLEAUX TAB-INTEGC ET TAB-INTEGI
C 
      CALL TABINT
C 
C     TRAITEMENT DES SYMETRIES
C 
      CALL TRASYM
C 
      CALL INFOEN('NUMDEV-MAT', ADMASO , NMAT  )
C 
C     REMPLISSAGE DU TABLEAU RANCAL-COU
C 
      CALL REPLAC
      CALL ADTBM ('RANCAL-COU', ADRANU)
C  
      IF (NBINT .GT. 0) THEN
C  
C     REMPLISSAGE DU TABLEAU RANCAL-INT
C 
        CALL REPLAI
C 
C     CALCUL ELEMENTAIRE CORRESPONDANT A U, V, W RANGES DANS L'ORDRE REEL
C 
        CALL ADTBM ('RANCAL-INT', ADRUI)
        CALL ADTBM ('RANCAL-INP', ADRUIP)
        CALL ADTBDM ('TAB-INTEGI', AINTGI)
C  
      END IF
C 
C     ADRAN U,V,W EST LA PREMIERE ADRESSE DU TABLEAU DES NUMEROS DE DDL
C     CALCUL ELEMENTAIRE CORRESPONDANT A U, V, W RANGES DANS L'ORDRE REEL
C 
      CALL ADTBM ('PRODL     ', ADPROD)
      CALL ADTBDM ('LONGUEURS ', ADLONG)
      CALL ADTBDM ('RAYON     ', ADRAYO)
      CALL ADTBDM ('TAB-INTEGC', ADINTE)
      CALL ADTBDM ('MATCOL    ', ADMATC)
      CALL ADTBM ('TLOCN1    ', ADTLOC)
C 
C     OUVERTURE DU TABLEAU PRECISION
C 
      CALL ADTBDM ('PRECISIONS', APRECI)
C  
C     PLACE RESERVEE POUR LA MATRICE K0N => UNE SEULE PLACE EST NECESSAIRE
C     ON EFFECTUE LA REMISE A ZERO NECESSAIRE DANS LA BOUCLE
C 
      CALL POUSMD (NTMAT, K0N)
C  
C     PLACE RESERVEE POUR LES TABLEAUX CALKIJ
C 
            CALL GSPOUD (9*144, K11)
            K11N = K11 +144
            K12  = K11N+144
            K22  = K12 +144
            K22N = K22 +144
            K13  = K22N+144
            K23  = K13 +144
            K33  = K23 +144
            K33N = K33+144
C 
C -----------------------------------------------------------------------
C     Ouverture d'un tableau partiel de longueur 64
C     pour le tableau TAB1 calcule pour chaque element par CALKI
C -----------------------------------------------------------------------
            CALL GSPOUD (2*64 , ADTAB1)
            ADTAB3 = ADTAB1+64
C 
C     PLACE RESERVEE POUR LES DEPLACEMENTS DES COUCHES
C 
            CALL GSPOUE (36, NDDLU)
            NDDLV  = NDDLU+12
            NDDLW  = NDDLV+12
C 
C     PLACE RESERVEE POUR LES DEPLACEMENTS DES INTERFACES
C 
            CALL GSPOUE (24, NDDLUI)
            NDDLVI = NDDLUI+8
            NDDLWI = NDDLVI+8
C 
C     PLACE RESERVEE POUR LE COMPORTEMENT DES COUCHES
C 
            CALL GSPOUD (10, ADCOMP)
C 
C     PLACE RESERVEE POUR LE COMPORTEMENT DES INTERFACES
C 
            CALL GSPOUD (3, ADCOMI)
C 
C     PLACE RESERVEE POUR LES NUMEROS DE COUCHES
C 
          CALL GSPOUE (NBCOU+1, MNBCOU)
          MNUCOU= MNBCOU+1
C 
C     PLACE RESERVEE POUR LES NUMEROS D'INTERFACES
C 
          CALL GSPOUE (NBINT+1, MNBINT)
          MNUINT=MNBINT+1
C 
C -----------------------------------------------------------------------
C     Ouverture du fichier ou sont stockees toutes les matrices
C -----------------------------------------------------------------------
C 
C     OUVERTURE DU FICHIER POUR LES DEFORMATIONS
C 
      NOMDIR = 'matricesm'
C 
      CALL CFDDNF (1, NOMDIR, 'matiso', 6, NTMAT, IUNIT)
C 
      NUENRS = 1
C 
C     MODIF CALL CREFIC (1, NOMDIR, 'MATISO', 6 ,'U', IUNIT)
C 
C     MODIF DO N = -NTDSFG , NTDSFG
C 
C     DEBUT DE BOUCLE SUR LE NOMBRE DE TERMES DE LA DECOMPOSITION EN SERIES DE FOURIER
C 
      DO IN = 1, NMAT
C 
        N =  M(ADMASO+IN-1)
C  
C       Remise a zero de la matrice de KON
C 
        CALL MENADM (K0N, NTMAT)
C 
C       DEBUT DE TEST SUR LE NOMBRE D'INTERFACES
C 
        IF (NBINT .GT. 0) THEN
C 
C         DEBUT DE BOUCLE i SUR LES INTERFACES
C 
          DO I = 1, NKINT
C 
C           On retrouve le comportement isotrope transverse correspondant a NKINT
C 
            CALL VALCOI (I, DM(ADCOMI))
C 
C           Recherche des numeros des interfaces concernees par le comportement precedent
C  
            CALL NUINT (I, M(MNBINT), M(MNUINT))
C  
C           DEBUT DE BOUCLE ii SUR LES INTERFACES CONCERNEES
C 
	    DO J = 1, M(MNBINT)
C 
              NUMINT=M(MNUINT+J-1)
              PASSP = .FALSE.
C 
              IF (SYMPAR .AND. NUMINT .EQ. 1) PASSP = .TRUE.
C  
C             DEBUT DE BOUCLE iii SUR LES COLONNES
C 
              DO K = 1, NBCOL
C 
                CALL NDDLIU (NUMINT, K, ADTLOC, M(NDDLUI),
     &                       M(NDDLVI), M(NDDLWI))
                CALL CALKI (DM(ADLONG+K-1), DM(ADRAYO+K-1), DM(ADCOMI),
     &                      AINTGI,DM(ADTAB1),DM(ADTAB3))
C 
                IF (PASSP) THEN
C 
                  CALL ASSBIP (N, K0N, M(ADRUIP), DM(ADTAB1),
     &                         DM(ADTAB3), M(NDDLUI), M(NDDLVI),
     &                         M(NDDLWI), ADPROD)
C 
                ELSE
C 
                  CALL ASSBI (N, K0N, M(ADRUI), DM(ADTAB1), DM(ADTAB3),
     &                        M(NDDLUI), M(NDDLVI), M(NDDLWI), ADPROD)
C 
                END IF
C 
CD              IF( N .EQ. -NTDSFG)  THEN
C 
CD                CALL IMPEN ( ' NUMINT ' , NUMINT )
CD                CALL IMPEN ( ' NUCOL  ' , K )
CD                CALL MESSAO ( ' IMPRESSION APRES LES INTERFACES ')
C 
CD                 CALL IMPDIA( K0N ,N)
C 
CD              END IF
C  
C             FIN DE BOUCLE iii SUR LES COLONNES
C 
              ENDDO
C  
C           FIN DE BOUCLE ii SUR LES INTERFACES CONCERNEES
C 
            ENDDO
C 
C         DEBUT DE BOUCLE i SUR LES INTERFACES
C 
          ENDDO
C 
CD        IF (LTRACN(1) .AND. N .EQ. 0) THEN
C 
CD          CALL IMPMAP (M(ADPROD), NDDL, DM(K0N))
C 
CD        ELSE IF (LTRACN(1) .AND. N .EQ. 2) THEN
C 
CD          CALL IMPMAP (M(ADPROD), NDDL, DM(K0N))
C 
CD        ELSE IF (LTRACN(1) .AND. N .EQ. 4) THEN
C 
CD         CALL IMPMAP (M(ADPROD), NDDL, DM(K0N))
C 
CD        END IF
C 
C       FIN DE TEST SUR LE NOMBRE D'INTERFACES
C 
        ENDIF
C 
C       DEBUT DE BOUCLE i SUR LE NOMBRE DE COUCHES
C 
        DO I = 1, NKCOU
C  
C         On retrouve le comportement isotrope transverse correspondant a NKCOU
C 
          CALL VALCOC (I, DM(ADCOMP))
C 
C         Recherche des numeros de couches concernees par le comportement precedent
C 
          CALL NUCOU (I, M(MNBCOU), M(MNUCOU))
C 
C         DEBUT DE BOUCLE ii SUR LES COUCHES CONCERNEES
C 
          DO J = 1, M(MNBCOU)
C 
            NCOU=M(MNUCOU+J-1)
            CALL VALEP (NCOU, EPAIS)
C 
C           DEBUT DE BOUCLE iii SUR LES COLONNES
C 
            DO K = 1, NBCOL
C 
C 
              TYPDEP    = 1
              CALL NDDLC (TYPDEP, NCOU, K, ADTLOC, M(NDDLU))
              TYPDEP    = 2
              CALL NDDLC (TYPDEP, NCOU, K, ADTLOC, M(NDDLV))
              TYPDEP    = 3
              CALL NDDLC (TYPDEP, NCOU, K, ADTLOC, M(NDDLW))
              CALL CALK11 (DM(ADLONG+K-1), EPAIS, DM(ADRAYO+K-1),
     &                     DM(ADCOMP), ADINTE, ADMATC+144*(K-1),
     &                     DM(K11), DM(K11N))
              CALL ASSBU (N, K0N, M(ADRANU), DM(K11), DM(K11N),
     &                    M(NDDLU), ADPROD)
              CALL CALK12 (DM(ADLONG+K-1), EPAIS, DM(ADRAYO+K-1),
     &                     DM(ADCOMP), ADINTE, ADMATC+144*(K-1),
     &                     DM(K12))
              CALL ASSBUV (N, K0N,M(ADRANU), DM(K12), M(NDDLU),
     &                     M(NDDLV), ADPROD)
              CALL CALK22 (DM(ADLONG+K-1), EPAIS, DM(ADRAYO+K-1),
     &                     DM(ADCOMP), ADINTE, ADMATC+144*(K-1),
     &                     DM(K22), DM(K22N))
              CALL ASSBV (N, K0N, M(ADRANU), DM(K22), DM(K22N),
     &                    M(NDDLV), ADPROD)
C 
CD            CALL IMPEP ('NCOU AV CALK13 ', NCOU)
C 
              CALL CALK13 (DM(ADLONG+K-1), EPAIS, DM(ADRAYO+K-1),
     &                     DM(ADCOMP), ADINTE, ADMATC+144*(K-1),
     &                     DM(K13))
              CALL ASSBUW (N, K0N, M(ADRANU), DM(K13), M(NDDLU),
     &                     M(NDDLW), ADPROD)
              CALL CALK23 (DM(ADLONG+K-1), EPAIS, DM(ADRAYO+K-1),
     &                     DM(ADCOMP),ADINTE, ADMATC+144*(K-1), DM(K23))
              CALL ASSBVW (N, K0N, M(ADRANU), DM(K23), M(NDDLV),
     &                     M(NDDLW), ADPROD)
              CALL CALK33 (DM(ADLONG+K-1), EPAIS, DM(ADRAYO+K-1),
     &                     DM(ADCOMP), ADINTE, ADMATC+144*(K-1),
     &                     DM(K33), DM(K33N))
              CALL ASSBW (N, K0N, M(ADRANU), DM(K33),
     &                    DM(K33N), M(NDDLW), ADPROD)
C 
C           FIN DE BOUCLE iii SUR LES COLONNES
C 
            ENDDO
C 
C         FIN DE BOUCLE ii SUR LES COUCHES CONCERNEES
C 
          ENDDO
C 
C       FIN DE BOUCLE i SUR LE NOMBRE DE COUCHES
C 
        ENDDO
C 
C       Determination de la precision test pour KON
C 
        CALL TERSUP( DM (K0N), NDDL , M(ADPROD) , VALSUP)
        EPS  = VALSUP * EPSILO
C       CALL IMPDT ('PRECISION EPS', EPS)
C  
C       MODIFICATION DE LA DIAGONALE EN FONCTION DES BLOCAGES ET DES
C       DEPLACEMENTS IMPOSES
C 
C       CALL IMPET ('POUR N DANS CALK0 AVANT BLOCAGE ', N)
C       CALL IMPDT ('VALSUP ', VALSUP)
C       CALL IMPDIA ( K0N ,N-NTDSFG-1)
C 
        CALL MOPRDI (N, DM(K0N), K0N, VALSUP)
        DM(APRECI+N+NTDSFG)  = VALSUP*GROS
C 
C       CALL IMPET('POUR N DANS CALK0 APRES BLOCAGE', N)
C       CALL IMPDIA( K0N ,N-NTDSFG-1)
C 
C       Decomposition profil de KON
C 
        CALL CHOLSP (NDDL, M(ADPROD), EPS, DM(K0N), DM(K0N), NPIV)
C 
C       Sequence d'impression eventuelle
C 
CD      IF (LTRACP(1)) THEN
CD        CALL IMPMAP( M(ADPROD) , NDDL ,DM(K0N)  )
CD        CALL IMPDIA( K0N ,N-NTDSFG-1)
CD        CALL IMPMAP( M(ADPROD) , NDDL ,DM(K0N)  )
CD      END IF
C 
C       Ouverture d'un fichier de nom K0//IDENTI pour ranger K0
C 
C       CALL IDENTI (N, NUMERO)
C 
C       NOM = 'MAT'//NUMERO
C 
C       NOMFIC = NOM
C       NOMDIR = 'matricesm'
C 
C       Ecriture de l'enregistrement correspondant a la matrice actuelle
C 
        CALL EFDDNF (IUNIT, DM(K0N), K0N, NTMAT, NUENRS)
        NUENRS = NUENRS +1
C 
C       CALL ECFICD( IUNIT , DM(K0N) , K0N ,NTMAT )
C  
C       CALL FERFIC( IUNIT , IDPROG )
C 
C     FIN DE BOUCLE SUR LES SERIES DE FOURIER
C 
      END DO
C 
      CALL FERFIC (1, IUNIT, IDPROG)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
      CALL IMPET ('ADM2 EN SORTIE DE '//IDPROG, ADM2)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
