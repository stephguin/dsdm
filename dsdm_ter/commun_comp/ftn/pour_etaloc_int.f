C     On envoie comme arguments :
C 
C     E ...... N0 numero de couche
C 
C     Et on recupere :
C 
C     S ...... TETA angle de la couche
C 
      SUBROUTINE ANGINT (N0, TETA)
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
      INTEGER N0, R
C 
      DOUBLE PRECISION TETA
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ANGCOU')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('ANGLES-COU',R)
      TETA=DM(R+N0+NBCOU-1)
CD    CALL IMPDP('VALEUR DE L''ANGLE',TETA)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule les quantites (contraintes ou deformations)
C     dans la base d'orthotropie a partir des quantites connues dans la base
C     locale; TETCAL est l'angle entre la base d'orthotropie et la base locale
C 
C     On envoie comme arguments :
C 
C     E ...... TETCAL  angle entre la base d'orthotropie et la base locale
C     E ...... QPLABL  quantite chapau dans la base locale
C 
C     Et on recupere :
C 
C     E ...... QPLABO  quantite chapau dans la base d'orthoropie

      SUBROUTINE QIBORT (TETCAL, QPLABL, QPLABO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION TETCAL , QPLABL(3) , QPLABO(3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION COS , SIN
C 
          CHARACTER*6 IDPROG
      PARAMETER (IDPROG='QIBORT')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      COS   = DCOS(TETCAL)
      SIN   = DSIN(TETCAL)
C 
      QPLABO( 1 ) =   COS*QPLABL( 1 )-SIN*QPLABL( 2 )
      QPLABO( 2 ) =   SIN*QPLABL( 1 )+COS*QPLABL( 2 )
      QPLABO( 3 ) =   QPLABL( 3 )
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule TRACE = TR(SAUT2)K(SAUT1))
C 
C     On envoie comme arguments :
C 
C     E ...... K      la matrice stockee (3)
C     E ...... SAUT1  la 1ere deformation (3)
C     E ...... SAUT2  la 2eme deformation (3)
C 
C     Et on recupere :
C 
C     S ...... TRACE
C 
      SUBROUTINE PROIOR (K, SAUT1, SAUT2, TRACE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
      DOUBLE PRECISION  TRACE, K(3), SAUT1(NSAU), SAUT2(NSAU)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER    I
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='PROIOR')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      TRACE = 0.D0
C 
      DO I = 1, 3
C 
        TRACE = TRACE + SAUT2(I)*K(I)*SAUT1(I)
C 
      END DO
C 
C     SEQUENCE D'IMPRESSION EN TRACE NORMALE
C 
CD    IF (LTRACN(1)) THEN
C 
CD      CALL IMPDN ('VALEUR DE LA TRACE ', TRACE)
C 
CD    END IF
C 
C     SEQUENCE D'IMPRESSION EN TRACE PROFONDE
C 
CD    IF (LTRACP(1)) THEN
C 
CD      CALL OMPTDP ('K EN ENTREE     ', K(1), 3, 1)
CD      CALL OMPTDP ('SAUT1 EN ENTREE ', SAUT1(1), 3, 1)
CD      CALL OMPTDP ('SAUT2 EN ENTREE ', SAUT2(1), 3, 1)
C 
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
