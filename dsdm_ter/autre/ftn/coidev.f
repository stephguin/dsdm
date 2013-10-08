C     !!! MODIF : Pour calculer directement la contribution a l'energie
C     on multiplie par pi les termes differents de 0 et par
C     2pi le terme 0 !!!
C 
C     Cette routine calcule a partir des valeurs complexes obtenues
C     en utilisant la transformation de Furier inverse des contraintes
C     reelles les coefficients du developpement en serie de Fourier de
C     NTDSFG termes ranges
C     - de 1 a nteta/2 :
C               N <  0      ==> (SINN,COSN,SINN)
C     - de nteta/2+1 a nteta  :
C               N > OU = 0  ==> (COSN,SINN,COSN)
C 
C     On envoie comme arguments :
C 
C     E ...... REECOM  partie reelle du developpement complexe
C                      range (3, NTDSFG+1)
C     E ...... IMACOM  partie imaginaire du developpement complexe
C                      range (3, NTDSFG)
C 
C     Et on recupere :
C 
C     S ...... SINDEV  developpement en serie de Fourier reel
C                      ranges matriciel composante des contraintes
C                      apres composante range (3, nbmat)
C 
      SUBROUTINE COIDEV (REECOM, IMACOM, SINDEV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  REECOM(3*(NTDSFG+1)), IMACOM(3*NTDSFG)
      DOUBLE PRECISION  SINDEV(3*NBMAT )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  S22, SIG22
      INTEGER  FIN, DEBUT, K, I, J, TABNIV(2), LONG
      INTEGER  AM2LC, ADM2LC
C 
CD    LOGICAL  LTRACN , LTRACP
C  
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='COIDEV')
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
CD     CALL OMPTDP('partie reelle en entree', REECOM(1) , 3 , NTDSFG+1)
CD     CALL OMPTDP('partie imaginaire en entree ',IMACOM(1),3, NTDSFG )
C 
CD    END IF
C 
C     Pour stocker les contraintes sigma 12 et igma 23 extraites
C 
      CALL POUSMD( NBMAT , S22 )
C 
C     Multiplication par 2 de la partie reelle autre que ceux
C     correspondant a TETA = 0 et teta= N
C 
      CALL HOMAT( 2.D0 , REECOM(4) , REECOM(4) , 3 , (NTDSFG-1) )
C 
C     Remplissage du tableau SINDEV dans l'ordre contraintes puis
C     sin N <------------> cos N.
C     O permutera ensuite la 2eme lignes de SINDEV correspondant a sigma 2.
C 
C     Pour ranger les coefficients des sinus => les 3 1er sont nuls
C     car ils correspondent a sinus (NTDSFG).
C 
      CALL BALAID( 3 , SINDEV(1) )
C 
C     Les 3*(NTDSFG-1) coefficients restant sont stockes dans IMATRA
C     jusqu'a IMATRA+3*(NTDSFG-1)-1 qui correspond au terme
C     sigma33( sin(ntdsfg); on veut ranger ce tableau a partir de SINDEV
C     de 7 dans un ordre decroissant en TETA et croissant en sigma.
C     On lit dans l'ordre inverse IMATRA => les contraintes sont
C     lues dans l'ordre inverse  et multipliee par -2.
C 
      FIN     = 3*NTDSFG
      DEBUT   = 7
      K       = 0
C 
      DO  I =1 , NTDSFG -1
C 
        DO J = 1 , 3
C 
          SINDEV( DEBUT -J )  = -2.D0*IMACOM( FIN -K )
          K                   = K+1
C 
        END DO
C 
        DEBUT  = DEBUT+3
C 
      END DO
C 
C     Rangement des coefficients des cosinus dans SINDEV
C 
C     TCD CALL COPITD( 3*(NTDSFG+1) , REECOM(1) , SINDEV(3*NTDSFG+1) )
C 
      CALL COPID( REECOM(1) , SINDEV(3*NTDSFG+1) , 3 , (NTDSFG+1) )
C 
C     Impression intermediaire de sigdev
C 
CD    IF (LTRACP(1))THEN
C 
CD     CALL OMPTDP( ' SINDEV avant rearrangement de sigma 22  '
CD                   , SINDEV(1) , 3 , NBMAT )
C 
CD    END IF
C 
C     SINDEV est rempli dans l'ordre contraintes puis
C     sin N <------------> cos N.
C     Permutation de 2eme  ligne de SINDEV
C     correspondant a sigma 22  pour la ranger 
C     cos N <------------> sin N
C 
      TABNIV(1) = 3
      TABNIV(2) = NBMAT
C 
C     Extraction de sigma 22
C 
      CALL EXTRAD (SINDEV, 2, TABNIV, 2, 2, DM(S22), LONG)
C 
      SIG22   = 2
C 
      DO I = 1 , NBMAT
C 
        SINDEV( SIG22)    = DM( S22 + NBMAT -I )
        SIG22             = SIG22+3
C 
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE

CD    IF (LTRACN(1))THEN
C 
CD     CALL OMPTDN( ' SINDEV en sortie ', SINDEV(1) , 3 , NBMAT )
C 
CD    END IF
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
