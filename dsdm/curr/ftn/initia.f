C      VERSION DU 02/10/86
C 
C 
       SUBROUTINE INITIA
C 
C ----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
C ----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES LOCAUX
C     """""""""""""""""""""""""""""""""
C 
      INTEGER I
C 
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='INITIA')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
C ----------------------------------------------------------------------
C 
C     INITIALISATION 
C 
C ----------------------------------------------------------------------
C 
C     Premiere adresse libre dans M
C 
      AM1=1
C 
C     Premiere adresse libre dans DM
C 
      ADM1=1
C 
C     Derniere adresse libre dans M
C 
      AM2  =LM
      AM2EN=LM
C 
C     Derniere adresse libre dans DM
C 
      ADM2EN=LDM
      ADM2  =LDM
C 
C     Nombre de tableaux dans M
C 
      NBTAM=0
C 
C     Nombre de tableaux dans DM
C 
      NBTADM=0
C  
C     Initialisation des longueurs des tableaux contenus dans M et DM
C 
C 
      DO I=1, NBTTM
        LONGM(I)=0
      ENDDO
      DO I=1, NBTTDM
        LONGDM(I)=0
      ENDDO
C 
C     Initialisation des adresses des tableaux contenus dans M et DM
C 
      DO I=1, NBTTM
        AM(I)=0
      ENDDO
      DO I=1, NBTTDM
        ADM(I)=0
      ENDDO
C 
C     Initialisation des longueurs des tableaux Kon
C 
      NTMAT=0
C 
C     Initialisation du nombre de termes du developpement
C     en serie de Fourier
C 
      NTDSF=0
C 
C     RENTREE DE PI, PI1 = DACOS(-1.D0)
C 
      PI  = DACOS(-1.D0)
C 
C     Coordonnees des points de gauss et valeurs des poids correspondants
C 
      GAUSS(1)   = 0.D0
      GAUSS(2)   = -1.D0/DSQRT(3.D0)
      GAUSS(3)   = -GAUSS(2)
      GAUSS(4)   = -DSQRT(3.D0/5.D0)
      GAUSS(5)   = 0.D0
      GAUSS(6)   = -GAUSS(4)
      GAUSS(7)   = -DSQRT((3.D0+2.D0*DSQRT(6.D0/5.D0))/7.D0)
      GAUSS(8)   = -DSQRT(( 3.D0-2.D0*DSQRT(6.D0/5.D0))/7.D0 )
      GAUSS(9)   = -GAUSS(8)
      GAUSS(10)  = -GAUSS(7)
C 
      POIDS(1) =2.D0
      POIDS(2) =1.D0
      POIDS(3) =1.D0
      POIDS(4) =5.D0/9.D0
      POIDS(5) =8.D0/9.D0
      POIDS(6) =5.D0/9.D0
      POIDS(7) =0.5D0-1.D0/(6.D0*DSQRT(6.D0/5.D0))
      POIDS(8) =0.5D0+1.D0/(6.D0*DSQRT(6.D0/5.D0))
      POIDS(9) =POIDS(8)
      POIDS(10)=POIDS(7)
C 
CD    CALL IMPTDN ('POINTS DE GAUSS ', GAUSS(1), 1, 10)
CD    CALL IMPTDN ('POIDS  ASSOCIES ', POIDS(1), 1, 10)
C 
C     NOMBRE DE FICHIERS OUVERTS
C 
      NBFIOU = 0
C 
C     NOMBRE TOTAL D'ETAPE LOCALE EFFECTUEE +1
C 
      NBETLT   = 1
C 
C     NOMBRE D'ETAPE LOCALE ACTUELLE EFFECTUEE +1
C 
      NBETLC   = 0
      NBETGL   = 0
      NBNETT   = 0
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
