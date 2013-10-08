C     On envoie comme arguments :
C  
C     E ...... K(9)     comportement travaillant avec les sauts
C     E ...... SAUT     Valeur des sauts
C 
C     Et on recupere :
C 
C     S ...... SIGMAN   Valeur des contraintes normales
C 
      SUBROUTINE IULORT (K, SAUT, SIGMAN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  K(3,3), SAUT(3), SIGMAN(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='IULORT')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      SIGMAN(1) = K(1,1) * SAUT(1) + K(2,1) * SAUT(2) + K(3,1) * SAUT(3)
      SIGMAN(2) = K(1,2) * SAUT(1) + K(2,2) * SAUT(2) + K(3,2) * SAUT(3)
      SIGMAN(3) = K(1,3) * SAUT(1) + K(2,3) * SAUT(2) + K(3,3) * SAUT(3)
C 
CD    CALL OMPTDP(' SIGMA NORMAL ', SIGMAN(1) , 1 , 3 )
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
