C     On envoie comme arguments :
C 
C     Et on recupere :
C 
      SUBROUTINE COMISO (DR1, NUCOU, ISO)
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
      INTEGER NUCOU, DR1, AVANT, I, M2, P
      DOUBLE PRECISION ISO(10)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='COMISO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBM ('TYP-COUCHE', M2)
C  
C     P caracterise le type de comportement
C 
      P=M(M2+NUCOU-1)
C 
      AVANT=DR1+10*(P-1)-1
C 
      DO I=1,10
        ISO(I)=DM(AVANT+I)
      ENDDO
C 
CD    CALL IMPTDN('TERMES DE ISO',ISO(1),1,10)
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
C     E ...... NUCOU nombre de couches
C 
C     Et on recupere :
C 
C     S ...... ISO tableau de 10 termes extraits de DM ??
C 
      SUBROUTINE COCISO (NUCOU, ISO)
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
      INTEGER NUCOU, DR1, AVANT, I, M2, P
      DOUBLE PRECISION ISO(10)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='COCISO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('HOO-COUCHE', DR1)
      CALL ADTBM ('TYP-COUCHE', M2)
C 
C     P CARACTERISE LE TYPE DE COMPORTEMENT
C 
      P=M(M2+NUCOU-1)
C 
      AVANT=DR1+10*(P-1)-1
C 
      DO I=1,10
        ISO(I)=DM(AVANT+I)
      ENDDO
C 
CD    CALL IMPTDN('TERMES DE ISO',ISO(1),1,10)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     QUE FAIT CETTE SUBROUTINE?:
C     cette subroutine retrouve l'adresse de depart du
C     comportement d'une couche:
C 
C     ...on envoie comme argument N0, numero de couche
C     et on recupere...M0, 1ere adresse du tableau de comportement
C     de la couche
C 
      SUBROUTINE COMCOU(N0,M0)
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
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='COMCOU')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM('HOO-COUCHE',R)
      CALL ADTBM('TYP-COUCHE',M2)
C 
C     P caracterise le type de comportement
C 
      P=M(M2+N0-1)
      M0=R+10*(P-1)
C 
CD    CALL IMPEP('VALEUR DU TYPE DE COMPORTEMENT',P)
CD    CALL IMPEP('VALEUR DE LA IERE ADRESSE M0',M0)
CD    CALL IMPTDP(' VALEUR DU COMPORTEMENT ', DM(M0),1,10)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... NP le numero de comportement
C 
C     Et on recupere :
C 
C     S ...... NC le nombre de couches concernees
C     S ...... C  les numeros de ces couches dans le tableau

         SUBROUTINE NUCOU (NP, NC, C)
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
      INTEGER NP, NC, C(NBCOU), I, K, M1, R1
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NUCOU')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ADTBM('TYP-COUCHE',M1)
      R1=M1-1
C 
      K=0
C 
      DO I=1,NBCOU
        IF(M(R1+I).EQ.NP)THEN
          K=K+1
          C(K)=I
        ENDIF
      ENDDO
C 
      NC=K
C 
CD     CALL IMPEP('POUR LE COMPORTEMENT DE COUCHENP',NP)
CD     CALL IMPEP('NOMBRE DE COUCHES CONCERNEES',NC)
CD     CALL IMPTEN('NUMEROS DES COUCHES CONCERNEES',C(1),1,NC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
         SUBROUTINE VALCOC(NP,ISO)
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
      INTEGER NP, DR1, AVANT, I
      DOUBLE PRECISION ISO(10)
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VALCOC')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ADTBDM('HOO-COUCHE',DR1)
      DR1   = DR1-1
      AVANT=DR1+10*(NP-1)
      DO I=1,10
        ISO(I)=DM(AVANT+I)
      ENDDO
C 
CD    CALL IMPTDN('TERMES DE ISO',ISO(1),1,10)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END

