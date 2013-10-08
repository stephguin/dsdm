C     VERSION DU 02/10/86
C 
C     ALLIX OLIVIER
C 
C     QUE FAIT CETTE SUBROUTINE?:
C 
C         cette routine remplit le tableau TABELE suivant
C 
C         la numerotation ligne si num=1;<=====>2*(nbcol+2).le.3*(nbcou+1)
C                          '''''                2*nbcol     .le.3*nbcou-1
C         la numerotation colonne si num=2;<===>2*(nbcol+2).gt.3*(nbcou+1)
C                          ''''''               2*nbcol    .gt.3*nbcou-1
C -----------------------------------------------------------------------
C         dans la numerotation ligne on a:
C 
C         l'element de la couche q de la colonne p a pour numero:
C                     (Q-1)NBCOL+P
C         et pour numero de premier noeud:
C                      3(NBCOL+1)(Q-1)+P
C         l'element de l'interface r de la colonne p a pour numero:
C                      NEL1 +(R-1)NBCOL+P
C         et pour numero de premier noeud:
C                      (3R-1)(NBCOL+1)+P
C -----------------------------------------------------------------------
C         dans la numerotation colonne on a:
C 
C         l'element de la couche q de la colonne p a pour numero:
C                     (P-1)NBCOU+Q
C         et pour numero de premier noeud:
C                      3(NBCOU)(P-1)+3Q-2
C         l'element de l'interface r de la colonne p a pour numero:
C                      NEL1 +(P-1)(NBCOU-1)+R
C         et pour numero de premier noeud:
C                      3(NBCOU(P-1)+R)
C -----------------------------------------------------------------------
C         on range sequentiellement les informations concernant
C         les elements comme suit:
C 
C         element couche(type 1)
C 
C         - type
C         - numero de couche
C         - numero de colonne
C         - numero de premier noeud    numerotation 4--------------3
C         - numero de second noeud                  5--------------6
C         - numero de troisieme noeud  choisie      1--------------2
C         - numero de quatrieme noeud
C         - numero de cinquieme noeud
C         - numero de sixieme noeud
C 
C         element interface(type 2)
C 
C         - type
C         - numero d' interface
C         - numero de colonne
C         - numero de premier noeud
C         - numero de second noeuCD     numerotation 4------- 3
C         - numero de troisieme noeud   choisie      1--------2
C         - numero de quatrieme noeud
C -----------------------------------------------------------------------
C       cette routine remplit le tableau de localisation des numeros
C       de noeuds ranges par ordre croissant par elements ranges couche
C       par couche ===> TLOCN1
C       a la lecture les numeros reels des noeuds sont ranges de
C       facon croissante par element
C       ensuite les numeros de noeud des interfaces sont ranges de
C       maniere identique
C 
      SUBROUTINE CATELE
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER C, E, M1
C 
C     L EST LA LONGUEUR DU TABLEAU TABELE
C 
      INTEGER L
C 
C     P EST UN INDICE DE BOUCLE SUR LES LIGNES
C 
      INTEGER P
C 
C     Q EST UN INDICE DE BOUCLE SUR LES COLONNES
C 
      INTEGER Q
C 
C     R EST UN INDICE DE BOUCLE SUR LES INTERFACES
C 
      INTEGER R, R1
C 
      INTEGER X, M2, K
C 
      INTEGER A
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CATELE')
C 
      E=NBCOL
      C=NBCOU
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     Calcul du type de numerotation
C 
      CALL CALNUM
C 
C     Ouverture d'un tableau d'entiers :
C        de nom  TABELE
C        de 1ere adresse M1
C        de numero N1
C 
      L=9*NEL1+7*NEL2
      CALL GESTEN ('TABELE    ', L, M1)
C 
C     On fait les valeurs de tabele a zero
C 
      CALL MENAM (M1, L)
C 
C     Ouverture d'un tableau d'entiers :
C        de nom  TLOCN1
C        de 1ere adresse M2
C        de numero N2
C 
      CALL GESTEN ('TLOCN1    ', NEL1*6+4*NEL2, M2)
      K=M2
C 
      IF (NUM .EQ. 1) THEN
C 
C     On est en numerotation ligne
C 
CD    CALL IMPMN ('RENTREE DANS LA NUMEROTATION LIGNE')
C 
        DO Q=1, C
          DO P=1, E
            A=M1+9*((Q-1)*E+P-1)
            M(A)=1
            M(A+1)=Q
            M(A+2)=P
            M(A+3)=3*(E+1)*(Q-1)+P
            M(K)=M(A+3)
            K=K+1
            M(A+4)=M(A+3)+1
            M(K)=M(A+4)
            K=K+1
            M(A+7)=M(A+3)+E+1
            M(K)=M(A+7)
            K=K+1
            M(A+8)=M(A+4)+E+1
            M(K)=M(A+8)
            K=K+1
            M(A+6)=M(A+4)+2*(E+1)-1
            M(K)=M(A+6)
            K=K+1
            M(A+5)=M(A+6)+1
            M(K)=M(A+5)
            K=K+1
          ENDDO
        ENDDO
C 
C       La sequence suivante n'est a activer que si il y a des interfaces
C 
        IF (NBINT .GT. 0) THEN
C 
CD        CALL IMPMN ('RENTREE DANS LA NUMEROTATION COLONNE POUR LES
CD                     INTERFACES')
C 
          DO R1=1, NBINT
            R = R1
C 
C           SYMPAR indique si le plan de symetrie est une interface
C 	  
            IF (SYMPAR) R = R1-1
            DO P=1,E
              X=M1+9*C*E+7*((R1-1)*E+P-1)
              M(X)      =2
              M(X+1)    =R1
              M(X+2)    =P
              M(X+3)    =(3*R-1)*(E+1)+P
              M(K)      =M(X+3)
              K         =K+1
              M(X+4)    =M(X+3)+1
              M(K)      =M(X+4)
              K         =K+1
              M(X+6)    =M(X+4)+E
              M(K)      =M(X+6)
              K         =K+1
              M(X+5)    =M(X+6)+1
              M(K)      =M(X+5)
              K         =K+1
            ENDDO
          ENDDO
        ENDIF
      ELSE
C 
C     On est en numerotation colonne
C  
        CALL IMPMN ('RENTREE DANS LA NUMEROTATION COLONNE')
        DO Q=1, C
          DO P=1, E
            A=M1+9*((P-1)*C+Q-1)
            M(A)=1
            M(A+1)=Q
            M(A+2)=P
            M(A+3)=3*C*(P-1)+3*Q-2
            M(K)=M(A+3)
            K=K+1
            M(A+7)=M(A+3)+1
            M(K)=M(A+7)
            K=K+1
            M(A+6)=M(A+3)+2
            M(K)=M(A+6)
            K=K+1
            M(A+4)=M(A+3)+3*C
            M(K)=M(A+4)
            K=K+1
            M(A+8)=M(A+4)+1
            M(K)=M(A+8)
            K=K+1
            M(A+5)=M(A+4)+2
            M(K)=M(A+5)
            K=K+1
          ENDDO
        ENDDO
C 
C       La sequence suivante n'est a activer que si il y a des interfaces
C 
        IF (NBINT .GE. 0) THEN
C 
CD        CALL IMPMN ('RENTREE DANS LA NUMEROTATION COLONNE POUR LES
CD                     INTERFACES')
C 
          IF (.NOT. SYMPAR) THEN
            DO R1=1, NBINT
              R = R1
              DO P=1, E
                X=M1+9*C*E+7*((P-1)*NBINT+R1-1)
                M(X)      = 2
                M(X+1)    = R1
                M(X+2)    = P
                M(X+3)    = 3*(C*(P-1)+R1)
                M(K)      = M(X+3)
                K         = K+1
                M(X+6)    = M(X+3)+1
                M(K)      = M(X+6)
                K         = K+1
                M(X+4)    = M(X+3)+3*C
                M(K)      = M(X+4)
                K         = K+1
                M(X+5)    = M(X+4)+1
                M(K)      = M(X+5)
                K         = K+1
              ENDDO
            ENDDO
          ELSE
C 
C           Pour l'interface 1
C 
            DO P=1, E
              A=M1+9*C*E+7*((P-1)*NBINT)
              M(A)=2
              M(A+1)=1
              M(A+2)=P
              M(A+3)=-100
              M(K)=M(A+3)
              K=K+1
              M(A+4)=-100
              M(K)=M(A+4)
              K=K+1
              M(A+6)= 3*C*(P-1)+1
              M(K)=M(A+6)
              K=K+1
              M(A+5)= M(A+6)+3*C
              M(K)=M(A+5)
              K=K+1
            ENDDO
C 
C           Pour les interfaces de 2 a NBINT
C 
            DO R1=1, NBINT-1
              DO P=1, E
                X=M1+9*C*E+7*((P-1)*NBINT+R1)
                M(X)      = 2
                M(X+1)    = R1+1
                M(X+2)    = P
                M(X+3)    = 3*(C*(P-1)+R1)
                M(K)      = M(X+3)
                K         = K+1
                M(X+6)    = M(X+3)+1
                M(K)      = M(X+6)
                K         = K+1
                M(X+4)    = M(X+3)+3*C
                M(K)      = M(X+4)
                K         = K+1
                M(X+5)    = M(X+4)+1
                M(K)      = M(X+5)
                K         = K+1
              ENDDO
            ENDDO
          END IF
        ENDIF
      ENDIF
C 
C     S'il y a des interfaces 
C 
CD    IF (NBINT .GE. 1) THEN
CD      CALL IMPTET ('TABLEAU DES ELEMENTS DE TYPE 2 ',
CD                    M(M1+9*NEL1), 7, NEL2)
CD      CALL IMPTET ('TABLEAU TLOCN1 POUR LES INTERFACES ',
CD                    M(M2+6*NEL1),4,NEL2)
CD    ENDIF
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule num :
C 
C           - num=1 si la numerotation est de type ligne
C           - num=2 si la numerotation est de type colonne
C 
      SUBROUTINE CALNUM
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
C     """""""""""""""""""""""""""""""""

      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALNUM')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      IF (NUM .EQ. -1) THEN
C 
        IF ((2*NBCOL) .LE. (3*NBCOU-1)) THEN
          NUM=1
        ELSE
          NUM=2
        ENDIF
C 
      END IF
C 
      CALL IMPET ('VALEUR DE NUM (TYPE DE NUMEROTATION) ',NUM)
      CALL MESSAO ('NUM=1 POUR NUMEROTATION LIGNE')
      CALL MESSAO ('NUM=2 POUR NUMEROTATION COLONNE')
C 
CD    IF ((NUM .NE. 1) .AND. (NUM .NE. 2))
CD    CALL ERREUD (0, 'PROBLEME DE NUMEROTATION '//IDPROG)
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
