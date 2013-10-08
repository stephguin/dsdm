C     Cette routine calcule a partir des valeurs developpees rangees
C     de 1 a nteta/2 :
C               N <  0      ==> (SINN,SINN,SINN,SINN,SINN,SINN)
C     de nteta/2+1 a nteta  :
C               N > OU = 0  ==> (COSN,COSN,COSN,COSN,COSN,COSN)
C     les valeurs developpes rangess
C     de 1 a nteta/2 :
C               N <  0      ==> (SINN,SINN,COSN,COSN,SINN,SINN)
C     de nteta/2+1 a nteta  :
C               N > OU = 0  ==> (COSN,COSN,SINN,SINN,COSN,COSN)
C 
C     On envoie comme arguments Et on recupere :
C 
C     ES ....... SIGDEV  developpement en serie de fourier reel
C                        ranges matriciel composante des contraintes
C              .         apres composante range (6, nbmat)
C 
      SUBROUTINE PERSIG (SIGDEV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  SIGDEV( 6,NBMAT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  I , J
C 
      DOUBLE PRECISION RESER1 , RESER2
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='PERSIG')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      DO  I = 3 , 4
C 
        DO J = 1 , NTDSFG
C 
          RESER1                  = SIGDEV( I , J)
          RESER2                  = SIGDEV( I , NBMAT+1-J )
          SIGDEV( I , J)          = RESER2
          SIGDEV( I , NBMAT+1-J ) = RESER1
C 
        END DO
C 
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
C 
CD     CALL OMPTDN( ' SIGDEV en sortie ', SIGDEV(1,1) , 6 , NBMAT )
C 
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END

