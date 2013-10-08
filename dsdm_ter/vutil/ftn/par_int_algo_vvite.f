C     On envoie comme arguments :
C 
C     E ...... SAUDEV  Tableau des sauts de deplacement pour un point
C                      de gauss range (3, NBMAT)
C 
C     Et on recupere :
C 
C     S ...... SAUCOS  Tableau des deformations correspondant aux
C                           COSINUS ( NTDSFG+1 ,3 )
C     S ...... SAUSIN  Tableau des deformations correspondant aux
C                           SINUS ( NTDSFG ,3 )
C 
      SUBROUTINE SASICO (SAUDEV, SAUCOS, SAUSIN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION     SAUDEV(3,NBMAT)
      DOUBLE PRECISION     SAUSIN(3*(NTDSFG+1)), SAUCOS(3*(NTDSFG+1))
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   N, KS, KC, NTDP1, NTDP2
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='SASICO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     REMPLISSAGE DANS L'ORDRE COS(0),...COS(N) <=> SAUCOS ,KC
C     REMPLISSAGE DANS L'ORDRE SIN(1),...SIN(N) <=> SAUSIN ,KS
C 
      KC = 1
      KS = 1
C 
CD    CALL IMPEP( 'NTDSFG DS SASICO POUR I=1,2 ',NTDSFG )
C 
      NTDP1 = NTDSFG + 1
      NTDP2 = NTDSFG + 2
      DO N = 1 , NTDSFG
        SAUSIN(KS) = SAUDEV( 1 , NTDP1-N)
        KS         = KS+1
      ENDDO
      DO N = NTDP1 , NBMAT
        SAUCOS(KC) = SAUDEV( 1 , N)
        KC         = KC+1
      ENDDO
C 
      DO N = 1 , NTDP1
        SAUCOS(KC) = SAUDEV( 2 , NTDP2-N )
        KC         = KC+1
      ENDDO
      DO N = NTDP2 , NBMAT
        SAUSIN(KS) = SAUDEV( 2 , N )
        KS         = KS+1
      ENDDO
C 
      DO N = 1 , NTDSFG
        SAUSIN(KS) = SAUDEV( 3 , NTDP1-N)
        KS         = KS+1
      ENDDO
      DO N = NTDP1 , NBMAT
        SAUCOS(KC) = SAUDEV( 3 , N)
        KC         = KC+1
      ENDDO
C 
CD    CALL OMPTDN('TABLEAU DES SAUCOS', SAUCOS(1) , NTDSFG+1 , 3 )
CD    CALL OMPTDN('TABLEAU DES SAUSIN', SAUSIN(1) , NTDSFG , 3 )
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
