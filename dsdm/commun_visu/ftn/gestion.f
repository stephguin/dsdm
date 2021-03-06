C     Cette routine change un tableau double (lon1,lon2,lon3)
C     en un tableau (lon1,lon3,lon2) le tableau est en Entree-Sortie.
C 
C     On envoie comme arguments :
C 
C     E ...... LONG Longueur du tableau a ranger
C     E ...... TAB  Tableau double a ranger
C 
      SUBROUTINE TAB132 (LON1, LON2, LON3, TAB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
      INTEGER           LON1, LON2, LON3
C 
      DOUBLE PRECISION  TAB(LON1, LON2, LON3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      INTEGER           LONG, DEBUT, I, J, K, DEB
C 
      CHARACTER*6       IDPROG
      PARAMETER        (IDPROG='TAB132')
C 
      INTEGER           AM2LC, ADM2LC, LDE, LFI
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      LONG = LON1*LON2*LON3
      CALL POUSMD (LONG, DEBUT)
      DEB = DEBUT
C 
      DO J= 1, LON2
        DO K = 1, LON3
  	  DO I = 1, LON1
            DM(DEB) = TAB(I, J, K)
  	    DEB = DEB+1
  	  END DO
         END DO	    
      END DO
C 
      CALL COPITD (LONG, DM(DEBUT), TAB)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
      
C     Cette routine change un tableau double (lon1,lon2,lon3)
C     en un tableau (lon3,lon2,lon1) le tableau est en Entree-Sortie.
C 
C     On envoie comme arguments :
C 
C     E ...... LONG Longueur du tableau a ranger
C     E ...... TAB  Tableau double a ranger
C 
      SUBROUTINE TAB321 (LON1, LON2, LON3, TAB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      include 'cominc.h'
C 
      INTEGER           LON1, LON2, LON3
C 
      DOUBLE PRECISION  TAB(LON1, LON2, LON3)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      INTEGER           LONG, DEBUT, I, J, K, DEB
C 
      CHARACTER*6       IDPROG
      PARAMETER        (IDPROG='TAB321')
C 
      INTEGER           AM2LC, ADM2LC, LDE, LFI
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      LONG = LON1*LON2*LON3
      CALL POUSMD (LONG, DEBUT)
      DEB = DEBUT
C 
      DO I= 1, LON1
        DO J = 1, LON2
  	  DO K = 1, LON3
            DM(DEB) = TAB(I, J, K)
  	    DEB = DEB+1
  	  END DO
         END DO	    
      END DO
C 
      CALL COPITD (LONG, DM(DEBUT), TAB)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
      
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     On envoie comme arguments :
C 
C     E ...... NUMERO   un entier (0 , 9)
C 
C     Et on recupere :
C 
C     S ...... COMPNO   le character*1 correspondant
C 
      SUBROUTINE IDENT1 (NUMERO, COMPNO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      CHARACTER*1  COMPNO
C 
      INTEGER      NUMERO
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='IDENT1')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      COMPNO = '0'
      IF (NUMERO .GE. 0 .AND. NUMERO .LT. 10) THEN
        WRITE (COMPNO(1:1), '(I1)') NUMERO
      END IF
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
      SUBROUTINE INFOME (LONGGE)
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
      INTEGER LONGGE
C 
C -----------------------------------------------------------------------
      LONGGE =LDM
      CALL IMPET ('LDM EN GESTION DANS INFOME ', LONGGE)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Recopie le tableau double en entree dans tabsor
C 
C     On envoie comme arguments :
C 
C     E ...... LONG    longueur totale du tableau
C     E ...... TABENT  tableau entree(long)
C 
C     Et on recupere :
C 
C     S ...... TABSOR  tableau sortie(long)
C 
      SUBROUTINE COPITD (LONG, TABENT, TABSOR)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER       LONG
C  
      DOUBLE PRECISION   TABENT(LONG), TABSOR(LONG)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER I
C 
C -----------------------------------------------------------------------
C 
      DO I =1, LONG
        TABSOR(I) = TABENT(I)
      END DO
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Recopie le tableau ENTIER en entree dans tabsor
C     On envoie comme arguments:
C 
C     On envoie comme arguments :
C 
C     E ...... LONG    longueur totale du tableau
C     E ...... TABENT  tableau en entree (long)
C 
C     Et on recupere :
C 
C     S ...... TABSOR  tableau en sortie (long)
C 
      SUBROUTINE COPITE (LONG, TABENT, TABSOR)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      INTEGER  LONG
C 
      INTEGER  TABENT(LONG), TABSOR(LONG)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
       INTEGER I
C 
      DO I =1, LONG
        TABSOR(I) = TABENT(I)
      END DO
C  
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine donne un petit coup de balai dans le tableau TAB.
C 
C     On envoie comme arguments :
C 
C     E ...... LONG Longueur du tableau a nettoyer
C     E ...... TAB  Tableau double a nettoyer
C 
      SUBROUTINE BALAID (LONG, TAB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER           LONG
C 
      DOUBLE PRECISION  TAB(LONG)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C  
      INTEGER         I
C 
C -----------------------------------------------------------------------
      DO I= 1 , LONG
        TAB(I) = 0.D0
      END DO
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
      SUBROUTINE BALAIE ( LONG , TAB )
C On envoie comme arguments:
C E................ LONG Longueur du tableau a nettoyer
C E................ TAB  Tableau entier a nettoyer
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
C
      INTEGER           LONG , TAB(LONG)
C
C
C**********************************************************************
C
C *   Declaration des parametres locaux
C     """""""""""""""""""""""""""""""""
C
      INTEGER         I
C
C***********************************************************************
      DO I= 1 , LONG
        TAB(I) = 0
      END DO
C -----------------------------------------------------------------------
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine nettoie le tableau M de ADRESS a ADRESS+LONG-1

      SUBROUTINE MENAM (ADRESS, LONG)
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER I, ADRESS, LONG
C 
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MENAM')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      DO I=ADRESS, ADRESS+LONG-1
         M(I)=0
      ENDDO
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine nettoie le tableau DM de ADRESS a ADRESS+LONG-1
C 
      SUBROUTINE MENADM (ADRESS, LONG)
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER I,ADRESS,LONG
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MENADM')
C 
CD    CALL WLKBCD(IDPROG)
CD    IF((LONG).LT.0)THEN
CD       CALL IMPET(' LONG ',LONG)
CD       CALL ERREUD(0,
CD     'LONGUEUR DE TABLEAU NEGATIVE  '//IDPROG)
CD    ENDIF
CD    IF((ADRESS+LONG-1).GT.LDMEFF)THEN
CD       CALL IMPET(' ADRESS+LONG-1> LDMEFF ',ADRESS+LONG-1)
CD       CALL ERREUD(0,
CD     'ECRITURE AU DELA DE LDMEFF  '//IDPROG)
CD    ENDIF
C 
      DO I = ADRESS, ADRESS+LONG-1
         DM(I)=0.D0
      END DO
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
C     E ....... NOM

C     Et on recupere :
C 
C     S ....... N
 
      SUBROUTINE NUTBM (NOM, N)
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
      INTEGER N,I
C 
      CHARACTER*(*)   NOM
      CHARACTER*10    NOML
C 
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG=' NUTBM')
C 
CD    CALL WLKBCD(IDPROG)
C 
      NOML = NOM
      DO I=1,NBTTM
        IF(NOML.EQ.CHARM(I))THEN
          N=I
C 
CD        CALL IMPCP('NOM DU TABLEAU',NOM)
CD        CALL IMPCP('NOM DU TABLEAU',NOML)
CD        CALL IMPEP('NUMERO CORRESPONDANT',N)
C 
          GOTO 1
        ENDIF
      ENDDO
          CALL IMPCT('NOM DU TABLEAU',NOM)
          CALL IMPCT('NOM DU TABLEAU',NOML)
          CALL ERREUD(0,'LE NOM DE TABLEAU DONT ON CHERCHE LE NUMERO '
     $                //'N''EXISTE PAS')
C 
1     CONTINUE
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine cherche le numero du tableau double precision
C     dont le nom est en entree
C 
      SUBROUTINE NUTBDM (NOM, N)
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
      INTEGER N,I
C 
      CHARACTER*(*) NOM
      CHARACTER*10 NOML
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NUTBDM')
C 
CD    CALL WLKBCD(IDPROG)
C 
      NOML = NOM
      DO I=1,NBTTDM
        IF(NOML.EQ.CHARDM(I))THEN
          N=I
C 
CD        CALL IMPCP('NOM DU TABLEAU',NOM)
CD        CALL IMPCP('NOM DU TABLEAU',NOML)
CD        CALL IMPEP('NUMERO DU TABLEAU',N)
C 
          GOTO 1
        ENDIF
      ENDDO
      CALL IMPCT('NOM DU TABLEAU', NOM)
      CALL IMPCT ('NOM DU TABLEAU', NOML)
      CALL ERREUD (0, 'LE NOM DE TABLEAU DONT ON CHERCHE LE NUMERO '
     $             //'N''EXISTE PAS')
C
1     CONTINUE
C
CD    CALL RETOUD(IDPROG)
C
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine cherche le numero du tableau double precision
C     dont le nom est en entree.
C 
C     On envoie comme arguments :
C 
C     E ...... NOM  le nom du tableau D.P
C 
C     Et on recupere :
C 
C     S ...... ADDEPA  adresse de depart du tableau
C 
      SUBROUTINE ADTBDM (NOM, ADDEPA)
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
      INTEGER ADDEPA, I ,N
C 
      CHARACTER*(*) NOM
      CHARACTER*10 NOML
C 
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ADTBDM')
C 
CD    CALL WLKBCD(IDPROG)
C 
      NOML = NOM
      DO I=1, NBTTDM
        IF (NOML .EQ. CHARDM(I)) THEN
          N       = I
          ADDEPA  = ADM(N)
CD        CALL IMPCP ('NOM DU TABLEAU', NOM)
CD        CALL IMPCP ('NOM DU TABLEAU', NOML)
CD        CALL IMPEP ('NUMERO DU TABLEAU', N)
CD        CALL IMPEP ('ADRESSE DE DEPART DU TABLEAU DANS DM', ADDEPA)
          GOTO 1
        ENDIF
      ENDDO
          CALL IMPCT ('NOM DU TABLEAU DOUBLE PRECISION', NOM)
          CALL IMPCT ('NOM DU TABLEAU DOUBLE PRECISION', NOML)
          CALL ERREUD (0,
     $        'LE NOM DE TABLEAU DONT ON CHERCHE L'' ADRESSE DE DEPART '
     $       //'N''EXISTE PAS')
C 
1     CONTINUE
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine cherche le numero du tableau d'entiers
C     dont le nom est en entree.
C 
C     Oon envoie comme argument :
C 
C     E ...... NOM le nom du tableau ENTIER
C 
C     Et on recupere :
C 
C     S ...... ADDEPA  adresse de depart du tableau
C 
      SUBROUTINE ADTBM (NOM,ADDEPA)
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
      INTEGER ADDEPA, I, N
C 
      CHARACTER*(*) NOM
      CHARACTER*10 NOML
C 
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='ADTBM ')
C 
CD    CALL WLKBCD(IDPROG)
C 
      NOML = NOM
      DO I=1,NBTTM
        IF(NOML.EQ.CHARM(I))THEN
          N       =I
          ADDEPA  =AM(N)
C 
CD        CALL IMPCP('NOM DU TABLEAU',NOM)
CD        CALL IMPCP('NOM DU TABLEAU',NOML)
CD        CALL IMPEP('NUMERO DU TABLEAU',N)
CD        CALL IMPEP('ADRESSE DE DEPART DU TABLEAU DANS M',ADDEPA)
C 
          GOTO 1
        ENDIF
      ENDDO
      CALL IMPCT('NOM DU TABLEAU ENTIER ', NOM)
      CALL IMPCT('NOM DU TABLEAU ENTIER ', NOML)
      CALL ERREUD(0,
     $ 'LE NOM DE TABLEAU DONT ON CHERCHE L'' ADRESSE DE DEPART '
     $ //'N''EXISTE PAS')
C 
1     CONTINUE
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C    Cette routine va chercher la premiere adresse libre dans le
C    tableau des entiers
C 
      SUBROUTINE DEBUEN (ADDEPA)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER ADDEPA
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='DEBUEN')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      ADDEPA  = AM1
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
      SUBROUTINE DEBUDP (ADDEPA)
C 
C -----------------------------------------------------------------------
C     Cette routine va chercher la premiere adresse libre dans le
C     tableau des double-precision
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER ADDEPA
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='DEBUDP')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      ADDEPA  = ADM1
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C 
C     Cette routine affecte de nouvelles valeurs aux dernieres adresses
C     libres dans les tableaux M et DM. On ecrase le reste.
C 
C     On envoie comme arguments :
C 
C     E ...... AM2LC   nouvelle derniere adresse libre dans le tableau M
C     E ...... ADM2LC  nouvelle derniere adresse libre dans le tableau DM
C 
      SUBROUTINE ENPOUB (AM2LC, ADM2LC)
C 
C ----------------------------------------------------------------------
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
      INTEGER     AM2LC, ADM2LC
C  
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='ENPOUB')
C 
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      ADM2LC  = ADM2
      AM2LC   = AM2
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     E ........ ENT1 ET ENT2 SONT COMPARES DANS LA ROUTINE NOMSUB
C 
      SUBROUTINE TESTEN (ENT1, ENT2, NOMSUB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER      ENT1 , ENT2
      CHARACTER*6 NOMSUB
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='TESTEN')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (ENT1 .NE. ENT2 )THEN
C 
        CALL IMPET ('VALEUR DE ENT1 ', ENT1)
        CALL IMPET ('VALEUR DE ENT2 ', ENT2)
        CALL ERREUD (0, 'LES ENTIERS SONT DIFFERENTS DANS '//NOMSUB)
      END IF
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine affecte de nouvelles valeurs aux dernieres adresses
C     libres dans les tableaux M et DM.
C 
C     On envoie comme arguments :
C 
C     E ...... AM2LC   ancienne derniere adresse libre dans le tableau M
C     E ...... ADM2LC  ancienne derniere adresse libre dans le tableau DM
C 
C     On recupere :
C 
C     S ...... AM2     nouvelle derniere adresse libre dans le tableau M
C     S ...... ADM2    nouvelle derniere adresse libre dans le tableau DM
C 
      SUBROUTINE SOPOUB (AM2LC, ADM2LC)
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER     AM2LC, ADM2LC
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='SOPOUB')
CD    INTEGER DIFMAX
C 
CD    CALL WLKBCD (IDPROG)
C 
C ----------------------------------------------------------------------- 
      ADM2 = ADM2LC
      AM2  = AM2LC
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     GEStion des Tableaux ENtiers
C 
C     On envoie comme arguments :
C 
C     E ...... L       Longueur du tableau
C     E ...... NOM     Nom du tableau d'entier
C 
C     Et on recupere :
C 
C     S ...... ADDEPA  Adresse de depart
 
      SUBROUTINE GESTEN (NOM, L, ADDEPA)
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER L, ADDEPA
C 
      CHARACTER*10  NOM
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='GESTEN')
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      ADDEPA          = AM1
      NBTAM           = NBTAM+1
      AM(NBTAM)       = AM1
      LONGM(NBTAM)    = L
      AM1             = AM1+L
C 
      CHARM(NBTAM)    = NOM
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     GEStion des Tableaux Double Precision
C 
C     On envoie comme arguments :
C 
C     E ...... NOM     Nom du tableau double precision
C     E ...... L       Longueur du tableau
C 
C     Et on recupere :
C 
C     S ...... ADDEPA  Adresse de depart
C 

      SUBROUTINE GESTDP (NOM, L, ADDEPA)
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER L, ADDEPA
C 
      CHARACTER*10  NOM
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='GESTDP')
C  
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
      IF ((ADM2-ADM1) .LT. L) THEN
        CALL IMPCT ('NOM DU TABLEAU ', NOM)
	CALL IMPET ('LONGUEUR DISPO POUR LE TABLEAU ', ADM2-ADM1)
        CALL IMPET ('LONGUEUR LIBRE POUR LE TABLEAU ', L)
        CALL IMPET ('ADM1                           ', ADM1)
        CALL IMPET ('ADM2                           ', ADM2)
        CALL ERREUD (0, 'LONGUEUR DE TABLEAU TROP GRANDE')
      ENDIF
C 
      ADDEPA         = ADM1
      NBTADM         = NBTADM+1
      ADM(NBTADM)    = ADM1
      LONGDM(NBTADM) = L
      ADM1           = ADM1+L
      CHARDM(NBTADM) = NOM
C 
      RETURN
      END 
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Creation du tableau provisoire d'entiers et nettoyage
C 
C     On envoie comme argument :
C 
C     E ...... L        la longueur du tableau a creer
C 
C     Et on recupere :
C 
C     S ...... ADDEPA   son adresse de depart
C 
C     + nettoyage
C 
      SUBROUTINE GSPOUE (L, ADDEPA)
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
C      """""""""""""""""""""""""""""""""
C 
      INTEGER L, ADDEPA
C 
C     L est la longueur du tableau provisoire
C 
      INTEGER I
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='GSPOUE')
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C 
CD    IF ((L) .LT. 0) THEN
CD      CALL IMPET ('LONG ', L)
CD      CALL ERREUD (0, 'LONGUEUR DE TABLEAU NEGATIVE '// IDPROG)
CD    ENDIF
C 
CD    IF ((AM2-L) .LT. AM1) THEN
CD      CALL IMPET ('AM1 ', AM1)
CD      CALL IMPET ('AM2 ', AM2)
CD      CALL IMPET ('LONGUEUR DU TABLEAU PROVISOIRE ', L)
CD      CALL ERREUD (0, 'ECRITURE SUR M AVANT LA 1 ADRESSE LIBRE '
CD   &              // IDPROG)
CD    ENDIF
CD    CALL IMPEP ('AM2 AVANT ', AM2)
CD    CALL IMPEP ('LONGUEUR DU TABLEAU PROVISOIRE ', L)
C 
      AM2   =AM2-L
C 
C     Nettoyage de M entre AM2 en sortie et AM2 en entree
C 
      DO I=0, L-1
         M(AM2+I)=0.D0
      ENDDO
C 
      ADDEPA =AM2
C 
CD    CALL IMPEP ('VALEUR ADDEPA ', ADDEPA)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Gestion du tableau provisoire de double precision et nettoyage
C 
C     On envoie comme argument :
C 
C     E ...... L        la longueur du tableau a creer
C 
C     Et on recupere :
C 
C     S ...... ADDEPA   son adresse de depart
C 
      SUBROUTINE GSPOUD (L, ADDEPA)

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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER L , ADDEPA
C 
C     L est la longueur du tableau provisoire
C 
      INTEGER I
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='GSPOUD')
CD    CALL WLKBCD (IDPROG)
C 
CD    IF ((L) .LT. 0) THEN
CD      CALL IMPET ('LONG ', L)
CD      CALL ERREUD (0, 'LONGUEUR DE TABLEAU NEGATIVE '// IDPROG)
CD    ENDIF
C 
CD    IF ((ADM2-L) .LT. ADM1) THEN
CD      CALL IMPET ('ADM1 ', ADM1)
CD      CALL IMPET ('ADM2 ', ADM2)
CD      CAL IMPET ('LONGUEUR DU TABLEAU PROVISOIRE ', L)
CD      CALL ERREUD(0,
CD                  'ECRITURE SUR DM AVANT LA 1 ADRESSE LIBRE '//IDPROG)
CD    ENDIF
CD    CALL IMPEP ('ADM2 AVANT ', ADM2)
CD    CALL IMPEP ('LONGUEUR DU TABLEAU PROVISOIRE ', L)
C 
      ADM2=ADM2-L
C 
C     Nettoyage de DM entre ADM2 en sortie et ADM2 en entree
C 
        DO I=0,L-1
          DM(ADM2+I)=0.D0
        ENDDO
C 
      ADDEPA =ADM2
C 
CD    CALL IMPEP ('VALEUR ADDEPA ', ADDEPA)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine cherche la longueur du tableau d'entiers
C     dont le nom est en entree
C 
C     ...ON ENVOIE COMME ARGUMENTS  NOM
C     ET ON RECUPERE...            LONG
C 
      SUBROUTINE LONGEN(NOM,LONG)
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
C     """""""""""""""""""""""""""""""""
C 
      INTEGER I , LONG
C 
      CHARACTER*(*) NOM
      CHARACTER*10  NOML
C 
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='LONGEN')
C 
CD    CALL WLKBCD(IDPROG)
C 
      NOML = NOM
      DO I=1,NBTTM
        IF(NOML.EQ.CHARM(I))THEN
          LONG = LONGM(I)
C 
CD        CALL IMPCP('NOM DU TABLEAU',NOM)
CD        CALL IMPCP('NOM DU TABLEAU',NOML)
CD        CALL IMPEP('NUMERO CORRESPONDANT',I)
C 
          GOTO 1
        ENDIF
      ENDDO

CD     IF((LONG).LT.0)THEN
CD       CALL IMPET(' LONG ',LONG)
CD       CALL ERREUD(0,
CD     'LONGUEUR DE TABLEAU NEGATIVE'//IDPROG)
CD     ENDIF
C 
      CALL IMPCT('NOM DU TABLEAU',NOM)
      CALL IMPCT('NOM DU TABLEAU',NOML)
      CALL ERREUD(0,
     &'LE NOM DE TABLEAU CHERCHE DANS LONGEN N''EXISTE PAS')
C 
1     CONTINUE
C 
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C ----------------------------------------------------------------------- 
C -----------------------------------------------------------------------
C 
C     Elle cherche la longueur du tableau double precision
C     dont le nom est en entree
C 
C     On envoie comme arguments :
C 
C     E ...... NOM
C 
C     Et on recupere :
C 
C     S ...... LONG
C 
      SUBROUTINE LONGDP(NOM,LONG)
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
      INTEGER I , LONG
C 
      CHARACTER*(*) NOM
      CHARACTER*10  NOML
C 
CD     CHARACTER*6 IDPROG
CD     PARAMETER (IDPROG='LONGDP')
C 
CD     CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD     IF((LONG).LT.0)THEN
CD       CALL IMPET(' LONG ',LONG)
CD       CALL ERREUD(0,
CD     'LONGUEUR DE TABLEAU NEGATIVE  '//IDPROG)
CD     ENDIF
C 
      NOML = NOM
      DO I=1,NBTTDM
        IF(NOML.EQ.CHARDM(I))THEN
          LONG = LONGDM(I)
C 
CD        CALL IMPCP('NOM DU TABLEAU',NOM)
CD        CALL IMPCP('NOM DU TABLEAU',NOML)
CD        CALL IMPEP('NUMERO CORRESPONDANT',I)
C 
          GOTO 1
        ENDIF
      ENDDO
      CALL IMPCT('NOM DU TABLEAU',NOM)
      CALL IMPCT('NOM DU TABLEAU',NOML)
      CALL ERREUD(0,
     &'LE NOM DE TABLEAU CHERCHE DANS LONGDP N''EXISTE PAS')
C 
1     CONTINUE
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     ...On envoie comme arguments:
C 
C     E ..... NOM
C 
C     Et on recupere :
C 
C     S ..... ADDEPA ,   LONG
C 
      SUBROUTINE INFOEN (NOM, ADDEPA, LONG)
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
      INTEGER I, LONG, ADDEPA
C 
      CHARACTER*(*) NOM
      CHARACTER*10  NOML
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='LONGEN')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NOML = NOM
      DO I = 1, NBTTM
        IF (NOML .EQ. CHARM(I)) THEN
          LONG   = LONGM(I)
          ADDEPA = AM(I)
          GOTO 1
        END IF
      END DO
      CALL IMPCT ('NOM DU TABLEAU', NOM)
      CALL IMPCT ('NOM DU TABLEAU', NOML)
      CALL ERREUD(0,
     &'LE NOM DE TABLEAU CHERCHE DANS LONGEN N''EXISTE PAS')
C 
1     CONTINUE
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     elle cherche la longueur du tableau double precision
C     dont le nom est en entree
C 
C     On envoie comme argument :
C 
C     E ......  NOM
C     E ......  ADDEPA adesse de depart dans DM
C 
C     Et on recupere :
C 
C     S ......  LONG
C 
      SUBROUTINE INFODP (NOM, ADDEPA, LONG)
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
      INTEGER I , LONG , ADDEPA
C 
      CHARACTER*(*) NOM
      CHARACTER*10  NOML
C 
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='INFODP')
C 
CD    CALL WLKBCD(IDPROG)
C 
CD     IF((LONG).LT.0)THEN
CD       CALL IMPET(' LONG ',LONG)
CD       CALL ERREUD(0,
CD     'LONGUEUR DE TABLEAU NEGATIVE  '//IDPROG)
CD     ENDIF
C  
      NOML = NOM
      DO I=1,NBTTDM
        IF(NOML.EQ.CHARDM(I))THEN
          LONG = LONGDM(I)
          ADDEPA  =ADM(I)
CD        CALL IMPCP('NOM DU TABLEAU',NOM)
CD        CALL IMPCP('NOM DU TABLEAU',NOML)
CD        CALL IMPEP('NUMERO CORRESPONDANT',I)
          GOTO 1
        ENDIF
      ENDDO
      CALL IMPCT('NOM DU TABLEAU',NOM)
      CALL IMPCT('NOM DU TABLEAU',NOML)
      CALL ERREUD(0,
     &'LE NOM DE TABLEAU CHERCHE DANS LONGDP N''EXISTE PAS')
C 
1     CONTINUE
C 
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... ADRESS 
C     E ...... + LONG < LDMEFF DANS LA
C     E ...... LA ROUTINE NOMSUB ???
C 
      SUBROUTINE TESTAD (ADRESS, LONG, NOMSUB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C
      INTEGER      ADRESS , LONG
      CHARACTER*6  NOMSUB
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='TESTAD')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD     IF((LONG).LT.0)THEN
CD       CALL IMPET(' LONG ',LONG)
CD       CALL ERREUD(0,
CD     'LONGUEUR DE TABLEAU NEGATIVE  '//IDPROG)
CD     ENDIF
C 
      IF (ADRESS+ LONG. GT . LDMEFF )THEN
        CALL ERREUD(0,
     &  ' ECRASEMENT DE DM DANS '//NOMSUB )
      END IF
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     ...ON ENVOIE COMME ARGUMENTS  LONG LA LONGUEUR DU FICHIER
C     ==> PAS DE REMISE A ZERO DE DM
C 
      SUBROUTINE POUSMD (LONG, ADDEPA)
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
      INTEGER  LONG, ADDEPA
C 
C 
C 
CD     CHARACTER*6 IDPROG
CD     PARAMETER (IDPROG='POUSMD')
C 
CD     CALL WLKBCD(IDPROG)
C 
CD      IF((ADM2-LONG).LT.ADM1)THEN
CD        CALL IMPET(' ADM1 ',ADM1)
CD        CALL IMPET(' ADM2 ',ADM2)
CD        CALL IMPET(' LONGUEUR DU TABLEAU PROVISOIRE',LONG)
CD        CALL ERREUD(0,
CD     'ECRITURE SUR DM AVANT LA 1 ADRESSE LIBRE '//IDPROG)
CD     ENDIF
C 
CD     IF((ADM2).GT.LDMEFF)THEN
CD       CALL IMPET(' LDMEFF ',LDMEFF)
CD       CALL IMPET(' ADM2 ',ADM2)
CD       CALL ERREUD(0,
CD     'ECRITURE SUR DM APRES LDMEFF '//IDPROG)
CD     ENDIF
CD     IF((LONG).LT.0)THEN
CD       CALL IMPET(' LONG ',LONG)
CD       CALL ERREUD(0,
CD     'LONGUEUR DE TABLEAU NEGATIVE  '//IDPROG)
CD     ENDIF
CD    CALL IMPEP(' ADM2 avant ',ADM2)
CD    CALL IMPEP(' LONGUEUR DU TABLEAU PROVISOIRE',LONG)
      ADM2   = ADM2 - LONG
      ADDEPA = ADM2
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C     ...ON ENVOIE COMME ARGUMENTS  LONG LA LONGUEUR DU FICHIER
C     ==> PAS DE REMISE A ZERO DE M
 
      SUBROUTINE POUSME( LONG , ADDEPA )
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
      INTEGER  LONG  , ADDEPA
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='POUSME')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF((AM2-LONG).LT.AM1)THEN
C 
CD        CALL IMPET(' AM1 ',AM1)
CD        CALL IMPET(' AM2 ',AM2)
CD        CALL IMPET(
CD        'LONGUEUR DU TABLEAU PROVISOIRE ENTIER',LONG)
CD        CALL ERREUD(0,
CD        'ECRITURE SUR M AVANT LA 1 ADRESSE LIBRE '//IDPROG)
CD    ENDIF
C 
CD    CALL IMPEP(' AM2 AVANT ',AM2)
CD    CALL IMPEP(' LONGUEUR DU TABLEAU PROVISOIRE',LONG)
CD     IF((AM2).GT.LM)THEN
CD       CALL IMPET(' LM ',LM)
CD       CALL IMPET(' AM2 ',AM2)
CD       CALL ERREUD(0,
CD     'ECRITURE SUR M APRES LM '//IDPROG)
CD     ENDIF
CD     IF((LONG).LT.0)THEN
CD       CALL IMPET(' LONG ',LONG)
CD       CALL ERREUD(0,
CD     'LONGUEUR DE TABLEAU NEGATIVE  '//IDPROG)
CD     ENDIF
C 
      AM2   = AM2 -LONG
      ADDEPA = AM2
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... NUMERO   un entier (1, 99)
C 
C     Et on recupere :
C 
C     S ...... COMPNO le character*2 correspondant
C 
      SUBROUTINE IDENT2 (NUMERO, COMPNO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      CHARACTER*2  COMPNO
C 
      INTEGER      NUMERO
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='IDENT2')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      COMPNO = '00'
      IF (NUMERO .GE. 0 .AND. NUMERO .LT.10) THEN
        WRITE (COMPNO(2:2), '(I1)') NUMERO
      ELSE IF (NUMERO .GE. 10) THEN
        WRITE (COMPNO(1:2), '(I2)') NUMERO
      END IF
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C 
C    On envoie comme arguments :
C 
C    E............ LONG    Longueur du tableau d'entiers correspondant
C                          a un nombre de termes a lire
C    E............ TABENT  Ce tableau
C    E............ TYPE    Un character*8 servant a l'identifier
C 
      SUBROUTINE POINTE(LONG, TABENT, TYPE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER          LONG, TABENT(LONG)
      CHARACTER*8      TYPE
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER           I, AR1, AP1
      CHARACTER*10      NOM
C 
CD     CHARACTER*6 IDPROG
CD     PARAMETER (IDPROG='POINTE')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NOM      = 'P-        '
      WRITE(NOM(3:10),'(A8)')TYPE
      CALL GESTEN( NOM , LONG+1, AP1 )
      AR1         = AP1-1
      M(AP1)      = 0
      DO I= 1, LONG
        M(AP1+I)      = M(AR1+I)+TABENT(I)
      ENDDO
C 
CD    CALL IMPTEN(' TABLEAU POINTEUR '//NOM,M(AP1),1,LONG+1)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... LISNOM  Liste des noms de tableaux a retirer
C     E ...... NBTAB   nombre de tableaux a retirer
C 
      SUBROUTINE RETIDP (LISNOM, NBTAB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER       NBTAB
C 
      CHARACTER*10  LISNOM(NBTAB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER    I, J, K
      INTEGER    LONRES
      INTEGER    DENUAD, DEBNUM, ADNUAD
      INTEGER    NBTBDB, NBTABA
      INTEGER    DEBLIS, ADLIS, CONT
      INTEGER    NUTAB, LONTAB, ADEPPA, ADTANC
C 
      CHARACTER*10  NOMTAB
C 
CD    CHARACTER*6 IDPROG
CD    PARAMETER (IDPROG='RETIDP')
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      LONRES = 2*NBTAB
C 
      CALL POUSME (LONRES ,DENUAD)
C 
      DEBNUM = DENUAD +NBTAB
C 
C     Recuperation des adresses de depart et des longueurs des
C     tableaux a supprimer
C 
      DO I = 1, NBTAB
C 
        ADNUAD = DENUAD +I-1
C 
        CALL NUTBDM (LISNOM(I), M(ADNUAD))
C 
      END DO
C 
      M(DEBNUM) = M(DENUAD)
C 
      DEBNUM = DEBNUM-1
      DENUAD = DENUAD -1
C 
      DO I = 2, NBTAB
C 
        ADNUAD = DENUAD +I
C 
C      Classement des numeros de ces tableaux par ordre croissant
C 
        DO J = 1, I-1
C 
          IF (M(ADNUAD) .LE. M(DEBNUM+J)) THEN
C 
            DO K = I-1, J , -1
C 
              M(DEBNUM+K+1) = M(DEBNUM+K)
C 
            END DO
C 
            M(DEBNUM+J) = M( ADNUAD)
C 
C          On ne passe qu'une fois dans le if la liste etant construite
C 
            GOTO 1000
C 
          END IF
C 
        END DO
C 
        M(DEBNUM+I) = M( ADNUAD)
C 
1000    CONTINUE
C 
      END DO
C 
CD    CALL IMPEN  (' NOMBRE TOTAL DE TABLEAUX AU DEPART ', NBTADM)
CD    CALL IMPTEN ('LISTE DES NUMEROS DE DEPART         ', M(DENUAD +1), 1, NBTAB)
CD    CALL IMPTEN ('LISTE DES NUMEROS CROISSANTS        ',M(DEBNUM+1) ,1, NBTAB)
C 
C     A partir du premier numero de tableau supprime on tasse =>
C     on recopie sur les tableaux supprimes.
C 
C     NBTBDB est le nombre de tableaux avant suppression.
C 
      NBTBDB = NBTADM
C 
C     NBTADM est remis au nombre de tableaux non rearranges.
C 
      NBTADM = M(DEBNUM+1)-1
      ADM1   = ADM(NBTADM)+LONGDM(NBTADM)
C 
C     NBTABA est le nombre de tableaux rearranges
C     Or on rearrange a partir de m(debnum+1) + 1
C     => dans une liste de  NBTBDB- NBTDAM-1 tableaux
C     tous les tableaux sauf (NBTAB-1) tableaux sont a supprimer.
C 
      NBTABA = NBTBDB-NBTADM-NBTAB
C 
      CALL POUSME (NBTABA, DEBLIS)
C 
      DEBLIS = DEBLIS-1
      ADLIS  = DEBLIS
      CONT   = 0
C 
C     Le tableau de numero M(DEBNUM+1) est a retirer. Boucle
C     quand meme pour les impressions et la comprehension.
C 
      DO I = M(DEBNUM+1), NBTBDB
C 
        DO K = 1, NBTAB
C 
          IF (I .EQ. M(DEBNUM+K)) THEN
C 
C         Le tableau est a supprimer. On ne recopie pas son
C         numero dans les numeros des tableaux a recopier.
C 
            CALL IMPET ('TABLEAU A SUPPRIMER NUMERO ', I)
            GOTO 2000

          END IF
C 
        END DO
C 
        CONT     = CONT+1
        ADLIS    = ADLIS+1
        M(ADLIS) = I
C 
2000    CONTINUE
C 
      END DO
C 
      CALL IMPET ('NOMBRE DE TABLEAUX REARRANGES PREVU   ', NBTABA)
      CALL IMPET ('NOMBRE DE TABLEAUX REARRANGES TROUVES ', CONT)
C 
CD    CALL IMPTET('LISTE DES NUMEROS DE TABLEAUX RECOPIES ',
CD                 M(DEBLIS+1), 1, CONT)
C 
      DO I = 1, NBTABA
C 
        NUTAB  = M(DEBLIS+I)
        LONTAB = LONGDM(NUTAB)
        NOMTAB = CHARDM(NUTAB)
C 
C      aAncienne adresse de depart dans DM
C 
        ADTANC = ADM(NUTAB)
C 
CD      CALL IMPET ('NUMERO DU TABLEAU RECOPIE ', NUTAB)
CD      CALL IMPCT ('NOM DU TABLEAU RECOPIE    ', NOMTAB)
C 
        CALL GESTDP (NOMTAB, LONTAB, ADEPPA)
        CALL COPITD (LONTAB, DM(ADTANC), DM(ADEPPA))
C 
CD      CALL IMPET ('NOUVEU NOMBRE DE TABLEAUX  ', NBTADM)
C 
      END DO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     Cette routine construit la chaine de caracteres pour l'objet
C     donne en argument.
C 
C     On envoie comme argument :
C 
C     E ...... NUMERO   un entier (-99, 999)
C 
C     Et on recupere :
C 
C     S ...... COMPNO   le character*3 correspondant
C 
      SUBROUTINE IDENTI (NUMERO, COMPNO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      CHARACTER*3  COMPNO
C 
      INTEGER   NUMERO
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='IDENTI')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      IF ((NUMERO .LT.-99) .OR. (NUMERO.GE. 1000)) THEN
C 
         CALL IMPET ('ERREUR NUMERO < 0 ou > 999 dans '//IDPROG,
     &   NUMERO )
C 
      END IF
C 
      COMPNO = '000'
      IF (NUMERO .LT. 0 .AND. NUMERO .GT.-10) THEN
        WRITE (COMPNO(2:3), '(I2)') NUMERO
      ELSE IF (NUMERO .LE.-10) THEN
        WRITE (COMPNO(1:3), '(I3)') NUMERO
      ELSE IF (NUMERO .GE. 0 .AND. NUMERO .LT.10) THEN
        WRITE (COMPNO(3:3), '(I1)') NUMERO
      ELSE IF (NUMERO .GE. 10) THEN
        WRITE (COMPNO(2:3), '(I2)') NUMERO
      ELSE IF (NUMERO .GE. 100) THEN
        WRITE (COMPNO(1:3), '(I3)') NUMERO
      END IF
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
  
