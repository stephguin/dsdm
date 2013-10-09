C     On envoie comme arguments :
C 
C     E ...... AD L'adresse de depart du tableau 'HOO-INTERF'
C                 ou  du tableau 'SOU-INTERF'
C     E ...... NP numero du comportement de l'interface
C 
C     Et on recupere :
C 
C     S ...... ISO(3) comportement isotrope transverse de l'interface
C                     range comme dans HOO-INTERF
C                     (E1+E2)/2  , E3 , (E1-E2)2
C 
      SUBROUTINE IOMISO (AD, NUINT, ISO)
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
      INTEGER NUINT, AD, AVANT, P, M2
      DOUBLE PRECISION ISO(3)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='IOMISO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBM('TYP-COUCHE',M2)
      P=M(M2+NBCOU+NUINT-1)
C 
C     P caracterise le type de comportement
C 
      AVANT=AD +3*(P-1)-1
      ISO(1)=DM(AVANT+1)
      ISO(2)=DM(AVANT+2)
      ISO(3)=DM(AVANT+3)
C 
C     IF( SYMPAR .AND. (NUINT .EQ.1) ) THEN
C 
C     Pour assurer qu'il n'y a pas de saut a l'interface centrale
C     on impose K1=K2=0 => la contrainte de cisaillement est nulle
C 
C     ISO(1) = 0.D0
C     ISO(2) = 2.D0*ISO(2)
C 
C     Le saut de w etant le double du deplacememnt normal de
C     l'interface on impose egalement
C 
C     ISO(3) = 2.D0*ISO(3)
C 
C     END IF
C 
CD    CALL IMPTDN('TERMES DE ISO',ISO(1),1,3)
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
C     E ...... NP numero du comportement de l'interface
C 
C     Et on recupere :
C 
C     S ...... ISO(3) comportement isotrope transverse de l'interface
C 
      SUBROUTINE COIISO (NUINT, ISO)
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
      INTEGER NUINT, AD, AVANT, P, M2
      DOUBLE PRECISION ISO(3)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='COIISO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('HOO-INTERF', AD)
      CALL ADTBM ('TYP-COUCHE', M2)
      P=M(M2+NBCOU+NUINT-1)
C 
C     P CARACTERISE LE TYPE DE COMPORTEMENT
C 
      AVANT=AD +3*(P-1)-1
      ISO(1)=DM(AVANT+1)
      ISO(2)=DM(AVANT+2)
      ISO(3)=DM(AVANT+3)
C 
CD    CALL IMPTDN('TERMES DE ISO',ISO(1),1,3)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... N0 numero d'interface
C 
C     Et on recupere :
C 
C     S ...... M0 1ere adresse du tableau de comportement
C                 isotrope de l'interface
C 
      SUBROUTINE COMINT(N0,M0)
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
      INTEGER N0,M0,M2,P,R
C 
      DOUBLE PRECISION ISO(3)
C 
      LOGICAL LTRACP
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='COMINT')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('HOO-INTERF',R)
      CALL ADTBM('TYP-COUCHE',M2)
      P=M(M2+NBCOU+N0-1)
C 
C     P caracterise le type de comportement
C 
      M0=R+3*(P-1)
C 
      IF ( LTRACP(1) )THEN
C  
CD      CALL IMPEP('VALEUR DU TYPE DE COMPORTEMENT',P)
CD      CALL IMPEP('VALEUR DE LA IERE ADRESSE M0',M0)
        CALL VALCOI (P, ISO)
CD      CALL OMPTDP(' ISO ', ISO(1),3,1)
C 
      END IF
CD    CALL OMPTDN(' COMPORTEMENT ISOTROPE ', DM(M0),3,1)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C    On envoie comme arguments :
C 
C    E ....... NP   numero du comportement de l'interface
C 
C    Et on recupere :
C 
C    S ....... ISO(3) comportement isotrope transverse de l'interface
C 
      SUBROUTINE VALCOI (NP, ISO)
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
      INTEGER NP, AD, DR1, AVANT
      DOUBLE PRECISION ISO(3)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VALCOI')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('HOO-INTERF',AD)
C 
CD    CALL IMPTDP('TABLEAU HOO-INTERF DANS '//IDPROG,DM(AD),1,3)
C 
      DR1=AD -1
      AVANT=DR1+3*(NP-1)
      ISO(1)=DM(AVANT+1)
      ISO(2)=DM(AVANT+2)
      ISO(3)=DM(AVANT+3)
C 
CD    CALL IMPTDN('TERMES DE ISO',ISO(1),1,3)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END

