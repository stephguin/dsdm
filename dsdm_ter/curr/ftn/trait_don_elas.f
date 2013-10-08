C     QUE FAIT CETTE ROUTINE?:
C     TRAINL traite les donnees de comportement pour remplir les
C     tableaux utiles
C 
C     ...ON ENVOIE COMME ARGUMENTS
C      KISYM = numero du comportement de l'interface 1
C     quand on est dans le cas SYMPAR

      SUBROUTINE TRAINL (KISYM)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
      INTEGER KISYM
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      INTEGER     I, NTTEI, AD2, AD4
      INTEGER     M1, R1, NT
      INTEGER     AM2LC, ADM2LC
C  
C     pour la rentree du comportement des couches
C 
      DOUBLE PRECISION     G12, G23, G13
      DOUBLE PRECISION     S11, S22, S12, S66, A1, A2, B1, B2, C
      DOUBLE PRECISION     K11, K22, K12, K66
      INTEGER              SOUPLI, RIGID, ASOUP, ADSOUP, NR
      INTEGER              ADHOOR, ADSOOR, NHOOKR, NSOUOR
C 
C     pour la rentree du comportement des interfaces
C 
      DOUBLE PRECISION     S1, S2, S3
      DOUBLE PRECISION     E1, E2, E3 , MAT33(3,3) , MATROT(3,3), TETLOC
      DOUBLE PRECISION     EX, EY, VXY, VYX, G , TEST1 , TEST2, TETA
      INTEGER              ADTETA
      INTEGER              NRI, ADSOUI, NHOIOR
      INTEGER              ADHOIN, ADSOIN
      INTEGER              LONSOUP
      INTEGER              AKCP, AKCPZ, AKCPZ2, AAZ
      LOGICAL              INDIC
C 
      CHARACTER*6 IDPROG
      PARAMETER  (IDPROG='TRAINL')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     Tableau pour ranger les types de comportement
C     concernant les couches(soit NP) on retrouvera la 1ere constante
C     d'elasticite de la couche a l'adresse 17*(NP-1)+1 dans HOOK-ORTHO
C -----------------------------------------------------------------------
C 
      CALL ADTBM ('TYP-COUCHE', M1)
      R1=M1-1
C 
      CALL INFODP ('SOUP-ORTHO', ADSOOR, LONSOUP)
      CALL IMPEP ('LONSOUP ', LONSOUP)
C 
      NKCOU = LONSOUP/17
      ASOUP = ADSOUP-1
C 
      CALL GESTDP ('HOO-COUCHE', 10*NKCOU, AD2)
      CALL GESTDP ('HOOK-ORTHO', 17*NKCOU, ADHOOR)
      CALL GESTDP ('SOU-COUCHE', 10*NKCOU, ADSOUP)
      ASOUP = ADSOUP-1
C 
      CALL POUSMD (18, SOUPLI)
      RIGID = SOUPLI+9
C 
      DO I = 0, NKCOU-1
C 
         NSOUOR        = ADSOOR+17*I
C 
         S11  =   DM(NSOUOR)
         S12  =   DM(NSOUOR+1)
         S22  =   DM(NSOUOR+4)
         S66  =   DM(NSOUOR+8)
         A1   =   DM(NSOUOR+13)
         A2   =   DM(NSOUOR+14)
         B2   =   DM(NSOUOR+12)
         B1   =   DM(NSOUOR+9)
         C    =   DM(NSOUOR+16)
C 
         G12   = 1.D0/S66
         G23   = 1.D0/B1
         G13   = 1.D0/B2
C 
C     Remplissage de sous-couches
C 
         NT   = 10*I+ADSOUP
C 
C -----------------------------------------------------------------------
C     transformation en valeurs interessantes
C     definition de la matrice isotrope transverse
C -----------------------------------------------------------------------
C 
C     TERME S11ISO
C 
         DM(NT)   = (3.D0*( S11+S22)+2.D0*(S12+S66))/8.D0
C 
C     TERME S12ISO
C 
         DM(NT+1) = (6.D0*S12+S11+S22-2.D0*S66)/8.D0
C 
C     TERMES S66ISO
C 
         DM(NT+2) = ( S11+S22+2.D0*(S66-S12))/4.D0
         DM(NT+3) = (A1+A2)/2.D0
         DM(NT+4) = (B1+B2)/2.D0
         DM(NT+5) =  C
C 
C     TERMES DU TYPE 2*TETA
C 
         DM(NT+6) = (S11-S22)/2.D0
         DM(NT+7) = (A1-A2)/2.D0
         DM(NT+8) = (B1-B2)/2.D0
C 
C     TERMES DU TYPE 4*TETA
C 
         DM(NT+9) = (S11+S22-2.D0*(S12+S66))/8.D0
C 
C -----------------------------------------------------------------------
C     Calcul des rigidites
C -----------------------------------------------------------------------
         DM( SOUPLI)     =  S11
         DM( SOUPLI+1)   =  S12
         DM( SOUPLI+2)   =  A1
         DM( SOUPLI+3)   =  S12
         DM( SOUPLI+4)   =  S22
         DM( SOUPLI+5)   =  A2
         DM( SOUPLI+6)   =  A1
         DM( SOUPLI+7)   =  A2
         DM( SOUPLI+8)   =  C
C 
         CALL INVMS3 (DM(SOUPLI), DM(RIGID))
C 
         K11  = DM( RIGID)
         K12  = DM( RIGID +1)
         K22  = DM( RIGID +4)
         K66  = G12
         A1   = DM( RIGID +2)
         A2   = DM( RIGID +5)
         B1   = G23
         B2   = G13
         C    = DM( RIGID +8)
C 
C -----------------------------------------------------------------------
C        REMPLISSAGE DE HOOK-ORTHO
C -----------------------------------------------------------------------
         NHOOKR       = ADHOOR+17*I
C 
         DM(NHOOKR)   = K11
         DM(NHOOKR+1) = K12
         DM(NHOOKR+2) = 0.D0
         DM(NHOOKR+3) = K12
         DM(NHOOKR+4) = K22
         DM(NHOOKR+5) = 0.D0
         DM(NHOOKR+6) = 0.D0
         DM(NHOOKR+7) = 0.D0
         DM(NHOOKR+8) = K66
C 
         DM(NHOOKR+9)  = B1
         DM(NHOOKR+10) = 0.D0
         DM(NHOOKR+11) = 0.D0
         DM(NHOOKR+12) = B2
C 
         DM(NHOOKR+13) = A1
         DM(NHOOKR+14) = A2
         DM(NHOOKR+15) = 0.D0
C 
         DM(NHOOKR+16) = C
C 
C -----------------------------------------------------------------------
C     Nr est la 1ere adresse libre dans HOO-COUCHE
C     taux des valeurs rangees isotrope transerse
C     + termes complemetaires
C 
C -----------------------------------------------------------------------
         NR=10*I+AD2
C 
CD       CALL IMPEP('NR=',NR)
C 
C -----------------------------------------------------------------------
C     definition de la matrice isotrope transverse
C -----------------------------------------------------------------------
C 
C     TER11ISO
C 
         DM(NR)=(3.D0*( K11+K22)+2.D0*(K12+K66))/8.D0
C 
C     TER12ISO
C 
         DM(NR+1)=(6.D0*K12+K11+K22-2.D0*K66)/8.D0
C 
C     TER66ISO
C 
         DM(NR+2)=( K11+K22+2.D0*(K66-K12))/4.D0
         DM(NR+3)=(A1+A2)/2.D0
         DM(NR+4)=(B1+B2)/2.D0
         DM(NR+5)=C
C 
C     TERMU TYPE 2*TETA
C 
         DM(NR+6)=(K11-K22)/2.D0
         DM(NR+7)=(A1-A2)/2.D0
         DM(NR+8)=(B1-B2)/2.D0
C 
C     TERMU TYPE 4*TETA
C 
         DM(NR+9)=(K11+K22-2.D0*(K12+K66))/8.D0
C 
      END DO
C 
C -----------------------------------------------------------------------
C     s'il y a  des interfaces-
C -----------------------------------------------------------------------
      CALL INFODP('SOIN-ORTHO',ADSOIN,NTTEI)
      NKINT = NTTEI/3
C 
      CALL GESTDP('SOU-INTERF', NTTEI , ADSOUI)
      CALL GESTDP('HOIN-ORTHO', NTTEI , ADHOIN)
      CALL GESTDP('HOO-INTERF',NTTEI,AD4)
C 
      DO I = 0 , NKINT-1
C 
         NT=3*I+ADSOIN
C 
         S1 = DM(NT)
         S2 = DM(NT+1)
         S3 = DM(NT+2)
C 
         NT=3*I+ADSOUI
C 
C -----------------------------------------------------------------------
C     1ere valeur pour la matrice isotrope transverse
C -----------------------------------------------------------------------
         DM(NT)  =(S1+S2)/2.D0
         DM(NT+2)=(S1-S2)/2.D0
C 
C -----------------------------------------------------------------------
C     2eme valeur pour la matrice isotrope transverse
C -----------------------------------------------------------------------
         DM(NT+1)=S3
C 
         IF (SYMPAR .AND. (KISYM .EQ. (I+1)) ) THEN
C 
           E1      = 0.D0
           E2      = 0.D0
           E3      = 1.D0/S3
           DM(NT)  = 0.D0
           DM(NT+2)= 0.D0
C 
         ELSE
C 
           E1 = 1.D0/S1
           E2 = 1.D0/S2
           E3 = 1.D0/S3
C 
         ENDIF
C 
C -----------------------------------------------------------------------
C     REMPLISSAGE DE HOIN-ORTHO
C -----------------------------------------------------------------------
         NHOIOR       = ADHOIN+3*I
C 
         DM(NHOIOR)   = E1
         DM(NHOIOR+1) = E2
         DM(NHOIOR+2) = E3
 
C     NRRI est la 1ere adresse libre dans HOO-INTERF
C 
         NRI=3*I+AD4
C 
C -----------------------------------------------------------------------
C     transformation en valeurs interessantes
C 
C     1ere valeur pour la matrice isotrope transverse
C 
         DM(NRI)  =(E1+E2)/2.D0
C 
C     2eme valeur pour la matrice isotrope transverse
C 
         DM(NRI+1)=E3
         DM(NRI+2)=(E1-E2)/2.D0
C 
      END DO
C 
CD    CALL IMPTET(' TYPES DE COMP DE COUCHE',M(M1),1,NBCOU)
C 
CD    CALL IMPTDT('TABLEAU SOU-COUCHE ',DM(ADSOUP),10,NKCOU)
CD    CALL IMPTDT('TABLEAU HOO-COUCHE ',DM(AD2)   ,10,NKCOU)
CD    CALL IMPTDT('TABLEAU SOUP-ORTHO ',DM(ADSOOR)   ,17,NKCOU)
C 
CD    CALL IMPTET(' TYPES DE COMP D''INTERFACES',M(M1+NBCOU),1,NBINT)
CD    CALL IMPTDT('TABLEAU SOU-INTERF ',DM(ADSOUI),3,NKINT)
CD    CALL IMPTDT('TABLEAU SOIN-ORTHO ',DM(ADSOIN),3,NKINT)
CD    CALL IMPTDT('TABLEAU HOO-INTERF ',DM(AD4)   ,3,NKINT)
CD    CALL IMPTDT('TABLEAU HOIN-ORTHO ',DM(ADHOIN),3,NKINT)
C 
C 
C     Pour test sur le calcul de la solution plaque
C 
      CALL GSPOUD( 30 , AKCP )
      AKCPZ   = AKCP  +9
      AKCPZ2  = AKCPZ +9
      AAZ     = AKCPZ2+9
C 
      CALL CARPLA( -THICK ,THICK ,
     &                DM(AKCP) , DM(AKCPZ) , DM(AKCPZ2)  )
C      
      CALL IMPTDT ( ' <KCP>   empilement ' , DM(AKCP), 3 ,3)
C  
      INDIC = .FALSE.  
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
      ADTETA = ADTETA -1
      DO I = 1 , NTETA
        TETLOC = DM(ADTETA+I)
	CALL RTCPLA( TETLOC , DM(AKCP) , MATROT(1,1) ) 
C 	
C 	TEST1 = DABS( MATROT(1,3)/MATROT(1,1) )
C 	TEST2 = DABS( MATROT(2,3)/MATROT(1,1) )
C 	IF ( (TEST1. LT. 1D-6) .AND. (TEST2 .LT. 1D-6) )THEN
C 	  INDIC = .TRUE.
C           CALL INVMS3( DM(AKCP) , MAT33(1,1) )
C 
C           EX  = 1.D0/(2.D0*THICK*MAT33(1,1) )
C           EY  = 1.D0/(2.D0*THICK*MAT33(2,2) )
C           VXY = -MAT33(1,2)/MAT33(1,1)
C           VYX = -MAT33(2,1)/MAT33(2,2)
C           G   = 1.D0/(2.D0*THICK*MAT33(3,3))
C 
C           CALL IMPDT ( ' DEMI EPAISSEUR ' , THICK)
C           CALL IMPDT ('ANGLE D''ORTHOTROPIE EQUIVALENT ', TETLOC)
C           CALL IMPDT( 'EX  POUR KCP '//IDPROG,EX)
C           CALL IMPDT( 'EY  POUR KCP '//IDPROG,EY)
C           CALL IMPDT( 'VXY POUR KCP '//IDPROG,VXY)
C           CALL IMPDT( 'VYX  POUR KCP '//IDPROG,VYX)
C           CALL IMPDT( 'VYX  POUR KCP '//IDPROG,VYX)	  
C 	  GOTO 1001
C 	   
C 	END IF 
      END DO  
1001  CONTINUE
      IF (.NOT. INDIC) THEN
        CALL MESSAO( ' ON NE TROUVE PAS DE BASE D''ORTHOTROPIE '//
     &    '\ POUR KCP DANS ' //IDPROG)	
      END IF
C
      CALL IMPTDT ( ' <KCPZ>  empilement ', DM(AKCPZ), 3 ,3)
      CALL IMPTDT ( ' <KCPZ2> empilement ', DM(AKCPZ2), 3 ,3)
C 
      CALL INVMS3( DM(AKCP) , MAT33(1,1) )
C 
      EX  = 1.D0/(2.D0*THICK*MAT33(1,1) )
      EY  = 1.D0/(2.D0*THICK*MAT33(2,2) )
      VXY = -MAT33(1,2)/MAT33(1,1)
      VYX = -MAT33(2,1)/MAT33(2,2)
      G   = 1.D0/(2.D0*THICK*MAT33(3,3))
C 
      CALL IMPDT ( ' DEMI EPAISSEUR ' , THICK)
      CALL IMPDT('ANGLE D''ORTHOTROPIE EQUIVALENT ', TETA)
      CALL IMPDT( 'EX  POUR KCP '//IDPROG,EX)
      CALL IMPDT( 'EY  POUR KCP '//IDPROG,EY)
      CALL IMPDT( 'VXY POUR KCP '//IDPROG,VXY)
      CALL IMPDT( 'VYX  POUR KCP '//IDPROG,VYX)
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
