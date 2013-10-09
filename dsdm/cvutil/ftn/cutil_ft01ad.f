C     Cette routine calcule les valeurs developpees des contraintes
C     (ou deformations) en fonction des valeurs reelles.
C 
C     On envoie comme arguments :
C 
C     E ...... SIGREE  (nteta, 6) la valeur de ces contraintes
C 
C     Et on recupere :
C 
C     S....... SIDEVF  (6, nbmat) la valeur des contraintes
C                      developpees en series de FOURIER rangees
C                      - de 1 a nteta/2 :
C                      N <  0 ==> (SINN, SINN, SIN, SIN, SINN, SINN)
C                      - de nteta/2+1 a nteta  :
C                      N > OU = 0 ==> (COSN, COSN, COSNN, COSN, COSN, COSN)
C 
      SUBROUTINE SIGDEV (SIGREE, SIDEVF)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION SIGREE( 6*NTETA ) , SIDEVF( 6*( NBMAT ) )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER INV  , SIGIMA , DEBREE(6) , DEBIMA(6) , I
      INTEGER     AM2LC      , ADM2LC   , REETRA , IMATRA
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SIGDEV')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Utilisation de la transformation de Fourier inverse qui a N
C     valeurs reelles associe N coeefficients du developpement en serie de
C     Fourier complexe.
C 
      INV = 1
C 
C     On recupere dans SIGREE la partie reelle et dans SIGCOM la partie
C     complexe qui en entree est " mise "a zero
C 
      CALL GSPOUD (6*NTETA, SIGIMA)
C 
C      DEBREE = 1
C      DEBIMA = SIGIMA
C 
C      write(6,*) ' --- sigdev avant appel ft01aD --- '
C      PAUSE ' SIGDEV AVANT FT01AD '
C      vd$l concur
C      vd$l cncall
C 
      DO  I = 1 , 6
C 
        DEBREE(I) = 1 + (I-1)*NTETA
        DEBIMA(I) = SIGIMA + (I-1)*NTETA
        CALL FT01AD( NTETA , INV , SIGREE(DEBREE(I)) , DM(DEBIMA(I)) )
C 
      END DO
C 
C      PAUSE ' SIGDEV APRES FT01AD '
C 
C      transposition du tableau REECOM <=> REETRA   ( 6 , NTDSFG+1 )
C      transposition du tableau IMACOM <=> IMATRA   ( 6 , NTDSFG )
C 
      CALL POUSMD( 12*NTETA , REETRA )
      IMATRA     = REETRA + 6*NTETA
C 
      CALL TRANSP( SIGREE(1)  , DM(REETRA) , NTETA , 6 )
      CALL TRANSP( DM(SIGIMA) , DM(IMATRA) , NTETA , 6 )
C 
      CALL COEDEV ( DM(REETRA) , DM (IMATRA)  , SIDEVF )
C 
      CALL PERSIG(  SIDEVF )
C 
      CALL SOPOUB (AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine calcule les valeurs developpees des contraintes
C     (ou deformations) en fonction des valeurs reelles.
C 
C     On envoie comme arguments :
C 
C     E ...... SIGREE  (nteta, 3) la valeur de ces contraintes
C 
C     Et on recupere :
C 
C     S ...... SIDEVF  (3, nbmat) la valeur des contraintes
C                      developpees en serie de Fourier rangees
C                      - de 1 a nteta/2 :
C                      N <  0      ==> (SINN,COSN,SINN)
C                      - de nteta/2+1 a nteta  :
C                      N > OU = 0  ==> (COSN,SINN,COSN)
C 
      SUBROUTINE SINDEV (SIGREE, SIDEVF)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION   SIGREE( 6*NTETA ) , SIDEVF( 6*( NBMAT ) )
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER INV  , SIGIMA , DEBREE(3) , DEBIMA(3) , I
      INTEGER     AM2LC      , ADM2LC   , REETRA , IMATRA
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SINDEV')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
C     Utilisation de la transformation de Fourier inverse qui a N
C     valeur reelles associe N coeefficients du developpement en serie de
C     Fourier complexe.
C 
      INV = 1
C 
C     On recupere dans SIGREE la partie reelle et dans SIGCOM la partie
C     complexe qui en entree est " mise " a zero
C 
      CALL GSPOUD ( 3*NTETA , SIGIMA )
C 
C     DEBREE = 1
C     DEBIMA = SIGIMA
C 
C     write(6,*) ' --- sindev avant appel ft01aD --- '
C     PAUSE ' SINDEV AVANT FT01AD '
C     vd$l concur
C     vd$l cncall
C 
      DO  I = 1 , 3
C 
        DEBREE(I) = 1 + (I-1)*NTETA
        DEBIMA(I) = SIGIMA + (I-1)*NTETA
        CALL FT01AD( NTETA , INV , SIGREE(DEBREE(I)) , DM(DEBIMA(I)) )
C 
      END DO
C 
C     write(6,*) ' --- sindev apres appel ft01aD --- '
C     PAUSE ' SINDEV APRES FT01AD '
C 
C 
C     Transposition du tableau REECOM <=> REETRA   ( 3 , NTDSFG+1 )
C     Transposition du tableau IMACOM <=> IMATRA   ( 3 , NTDSFG )
C 
      CALL POUSMD( 6*NTETA , REETRA )
      IMATRA     = REETRA + 3*NTETA
C 
      CALL TRANSP( SIGREE(1)  , DM(REETRA) , NTETA , 3 )
      CALL TRANSP( DM(SIGIMA) , DM(IMATRA) , NTETA , 3 )
C 
      CALL COIDEV ( DM(REETRA) , DM (IMATRA)  , SIDEVF )
C 
      CALL SOPOUB(AM2LC,ADM2LC)
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
C     E ...... SAUDEV  Tableau des sauts pour un point range 
C                      (3, NBMAT)
C 
C     Et on recupere :
C 
C     S ...... SAREEL  Tableau des sauts reels stockes
C                      (3, NTETA)
C 
      SUBROUTINE SAUREE (SAUDEV, SAREEL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION  SAREEL(3*NTETA), SAUDEV(3*NBMAT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
CD    LOGICAL LTRACP , LTRACN
C 
      INTEGER   ADCOS, ADSIN, ADREE, ADIMA
      INTEGER   DEBCOS(3), DEBSIN(3), DEBREE(3), DEBIMA(3)
      INTEGER   I, INV
      INTEGER   AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SAUREE')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C -----------------------------------------------------------------------
C     IMPRESSION EN ENTREE
C 
CD    CALL OMPTDN (
CD    'SAUTS DEVELOPPES EN ENTREE POUR LE POINT DE GAUSS',
CD    SAUDEV(1) ,  3 , NBMAT)
C 
      CALL GSPOUD( 6*NTDSFG+6   + 6*NTETA , ADCOS)
C  
C     Adresse des debuts de tableaux provisoires
C 
      ADSIN   = ADCOS+3*(NTDSFG+1)
      ADIMA   = ADSIN+3*NTDSFG
      ADREE   = ADIMA+3*NTETA
C 
      DEBCOS(1)  = ADCOS
      DEBSIN(1)  = ADSIN
      DEBIMA(1)  = ADIMA
      DEBREE(1)  = ADREE
C 
C     vd$l novector
C     vd$l noconcur
C 
      DO I=2,3
        DEBCOS(I) = DEBCOS(I-1) + NTDSFG+1
        DEBSIN(I) = DEBSIN(I-1) + NTDSFG
        DEBREE(I) = DEBREE(I-1) + NTETA
        DEBIMA(I) = DEBIMA(I-1) + NTETA
      ENDDO
C 
      CALL SASICO( SAUDEV , DM(ADCOS) , DM(ADSIN) )
C 
      INV = 2
C 
C     PAUSE ' SAUREE AVANT FT01AD '
C 
C    vd$l concur
C    vd$l cncall 
C 
      DO I = 1 , 3
C 
CD      CALL IMPEP( 'POUR LE SAUT NUMERO ', I )
C 
        CALL TABTDF(NTDSFG, DM(DEBCOS(I)), DM(DEBSIN(I)),
     &               DM(DEBREE(I)), DM(DEBIMA(I)))
C 
C     On recupere le tableau de longueur 2n des coefficients reels du
C     developpement en series d'exponentielles puis le tableau de longueur
C     2n des coefficients imaginaires du developpement en series
C     d'exponentielles .
C 
CD      CALL IMPEP( 'NTDSFG AV FT01AD ',NTDSFG )
C 
        CALL FT01AD( NTETA , INV , DM(DEBREE(I)) , DM(DEBIMA(I)) )
C 
CD      CALL IMPEP( 'NTDSFG AP FT01AD ',NTDSFG )
C 
C     Incrementation pour aller lire les tableaux des cosinus et des sinus
C 
CD      CALL OMPTDP('VALEURS REELLES DES SAUTS ',
CD                   DM(DEBREE(I)) , NTETA ,1)
CD      CALL OMPTDP('VALEURS IMAGINAIRES DES SAUTS ',
CD                   DM(DEBIMA(I)) , NTETA ,1)
C 
      ENDDO
C 
C     PAUSE ' SAUREE APRES FT01AD '
C 
CD    IF (LTRACP(1) ) THEN
C 
CD      CALL OMPTDP( ' TABLEAU DES SAUTS REELS DANS ' // IDPROG
CD       , DM(ADREE) , NTETA , 3 )
CD      CALL OMPTDP( ' TABLEAU DES SAUTS IMAGINAIRES DANS ' // IDPROG
CD       , DM(ADIMA) , NTETA , 3 )
C 
C     Pour test sur l'operation inverse
C 
CD    CALL POUSMD( 3*NTETA , ATEST )
CD    CALL SINDEV( DM(ADREE) , DM(ATEST) )
CD    CALL OMPTDN( '  SAUT DEV = ENTREE?? ' // IDPROG
CD     , DM(ATEST),  3 , NBMAT )
CD    END IF
C 
      CALL TRANSP( DM(ADREE) , SAREEL(1) , NTETA , NSAU )
C 
CD    CALL OMPTDN( ' TABLEAU DES SAUTS REELS DANS ' // IDPROG
CD     , SAREEL(1) ,  NSAU , NTETA )
C 
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On peut confondre SM et VRTSM en entree SM(NDDL, NBMAT), VRSM(NTETA, NDDL)
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT    Nombre de fonctions du temps nombre de seconds membres
C     ES...... SM        Valeurs developpees de tous les seconds
C                        membres ranges calcul
C  
C     Et on recupere :
C  
C     ES...... VRSM      Valeurs reelles de tous les seconds membres
C                        ranges ddl par ddl
C 
      SUBROUTINE VRTSM (NFONCT, SM, VRSM)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER   NFONCT
      DOUBLE PRECISION  SM(NBMAT*NDDL*NFONCT), VRSM(NTETA*NDDL*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER    NG, NG2, LONG, AD2, AD22, DEBSIN, DEBCOS
      INTEGER    DEBREE, DEBIMA, INV, I, NFT, DEBUT
      INTEGER    AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER  (IDPROG='VRTSM ')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL IMPET ('ADM2 EN ENTREE DANS '//IDPROG, ADM2)
C 
      NG     = 2*NTDSFG+1
      NG2    = NTETA
C 
C     MODIF 21-12      CALL GSPOUD (LONG*NFONCT, AD2)
C 
      LONG   = NG2*NDDL
      CALL POUSMD (NG2*NDDL+NDDL*NG*NFONCT, AD2)
      AD22   = AD2+NDDL*NG*NFONCT
      DEBSIN = AD2
      DEBCOS = AD2+NTDSFG
      DEBREE = 1
      DEBIMA = AD22
      INV    = 2
C 
C     MODIF 25
C 
      DEBUT  = 1
C 
      DO NFT = 1, NFONCT
C 
C       M(1) est ici un argument inutile
C 
        CALL DESICO (M(1), NDDL, SM(DEBUT), DM(AD2))
C  
C       modif
C 
        DEBSIN = AD2
        DEBCOS = DEBSIN+NTDSFG
        DEBIMA = AD22
C 
C       fin de modif
C 
        DO I= 1,NDDL
C 
          CALL TABTDF (NTDSFG, DM(DEBCOS), DM(DEBSIN), VRSM(DEBREE),
     &                 DM(DEBIMA))
          CALL FT01AD (2*NTDSFG, INV, VRSM(DEBREE), DM(DEBIMA))
C 
          DEBSIN  = DEBSIN+NG
          DEBCOS  = DEBCOS+NG
          DEBREE  = DEBREE+NG2
          DEBIMA  = DEBIMA+NG2
        END DO
C 
C       MODIF 25
C 
        DEBUT = NBMAT*NDDL+DEBUT
C 
      END DO
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
