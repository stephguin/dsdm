C     !!! MODIF : Pour calculer directement la contribution a l'energie
C     on nultiplie par pi les termes different de 0 et par
C     2pi le terme 0 !!!
C 
C     Cette routine calcule a partir des valeurs complexes  obtenues
C     en utilisant la transformation de fourier inverse des contraintes
C     reelles les coefficients du developpement en serie de fourier de
C     NTDSFG termes ranges de 1 a nteta/2 :
C               N <  0      ==> (SINN,SINN,SINN,SINN,SINN,SINN)
C                          de nteta/2+1 a nteta  :
C               N > OU = 0  ==> (COSN,COSN,COSN,COSN,COSN,COSN)
C 
C     On envoie comme arguments :
C 
C     E ...... REECOM  partie reelle du developpement complexe
C                      range (6, NTDSFG+1)
C     E ...... IMACOM  partie imaginaire du developpement complexe
C                      range (6, NTDSFG)
C 
C     Et on recupere :
C 
C     S................ SIGDEV  developpement serie de Fourier reel
C                               ranges matriciel composante des contraintes
C                               apres composante range (6, nbmat)
C 
      SUBROUTINE COEDEV (REECOM, IMACOM, SIGDEV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  REECOM(6*(NTDSFG+1)), IMACOM(6*NTDSFG)
      DOUBLE PRECISION  SIGDEV(6*NBMAT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  FIN, DEBUT, K, I, J
      INTEGER  AM2LC, ADM2LC
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='COEDEV')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     SEQUENCE D'IMPRESSION EN TRACE PROFONDE
C 
CD    IF (LTRACP(1))THEN
C 
CD     CALL OMPTDP('partie reelle en entree', REECOM(1) , 6 , NTDSFG+1)
CD     CALL OMPTDP('partie imaginaire en entree ',IMACOM(1),6, NTDSFG )
C 
CD    END IF
C 
C     Multiplication par 2 de la partie reelle autre que ceux
C     correspondant a TETA = 0 et teta = NTDSFG
C 
C     TCD CALL MUMARE( 2.D0 , 6*(NTDSFG-1) , REECOM(7) , REECOM(7) )
C 
      CALL HOMAT( 2.D0 , REECOM(7) , REECOM(7) , 6 , (NTDSFG-1) )
C 
C     Remplissage du tableau SIGDEV dans l'ordre contraintes puis
C               sin N <------------> cos N
C     On permutera ensuite les 3eme et 4eme  lignes de SIGDEV
C     correspondant a sigma 12 et sigma 23
C 
C     Pour ranger les coefficients des sinus => les 6 1er sont nuls
C     car ils correspondent a sinus(  NTDSFG )
C 
      CALL BALAID( 6 , SIGDEV(1) )
C 
C     Les 6*(NTDSFG-1) coefficients restants sont stockes dans IMATRA
C     jusqu'a IMATRA+6*(NTDSFG-1)-1 qui correspond au terme
C     sigma33 (sin(ntdsfg)); on veut ranger ce tableau a partir de SIGDEV
C     de 7 dans un ordre decroissant en TETA et croissant en sigma.
C     On lit dans l'ordre inverse IMATRA => les contraintes sont
C     lues dans l'ordre inverse  et multipliee par -2
C 
      FIN     = 6*NTDSFG
      DEBUT   = 13
      K       = 0
C 
      DO  I =1 , NTDSFG -1
C 
        DO J = 1 , 6
C 
          SIGDEV( DEBUT -J )  = -2.D0*IMACOM( FIN -K )
          K                   = K+1
C 
        END DO
C 
        DEBUT  = DEBUT+6
C 
      END DO
C 
C     Rangement  des coefficients des cosinus dans SIGDEV
C 
C     TCD CALL COPITD (6*(NTDSFG+1), REECOM(1), SIGDEV(6*NTDSFG+1))
C 
      CALL COPID (REECOM(1), SIGDEV(6*NTDSFG+1), 6, (NTDSFG+1))
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1))THEN
C 
CD    CALL OMPTDN( ' SIGDEV en sortie ', SIGDEV(1) , 6 , NBMAT )
C 
CD    END IF
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
