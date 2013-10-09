C     On envoie comme arguments :
C 
C     E ...... N        numero de la matrice ( -ntdsfg , ntdsfg )
C     E ...... MATRIC   matrice dont on modifie la diagonale
C                       pour tenir compte a la fois
C                       des blocages et deplacements imposes
C     E ...... ADRESS   adresse de depart de la matrice dans DM
C     E ...... TERSUP   + grand terme de la diagonale
C 
      SUBROUTINE MOPRDI (N, MATRIC, ADRESS, TERSUP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER ADRESS   , N
      DOUBLE PRECISION MATRIC(NTMAT)  , TERSUP
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER    ADPRO, ADBLO, LONBLO
C 
      INTEGER    ADNBLO, LDVBLO, DEBEFF, FINEFF, DBBLO
      INTEGER    NTDIAG, TERDIA, DBDIA, I, LONDIA
      INTEGER    AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MOPRDI')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL TESTAD( ADRESS , NTMAT , IDPROG)
C 
C     SEQUENCE DE MODIFICATION DU PROFIL DES MATRICES K0N
C 
      CALL ADTBM('PRODL     ',ADPRO)
C 
C     POUR LES BLOCAGES => doit etre valable pour tout
C                          les numeros de developpement
C 
      CALL ADTBM ('DDL-BLOQUE',ADBLO)
      CALL LONGEN('DDL-BLOQUE',LONBLO)
C 
C     ADEP-DVBLO EST LE TABLEAU DES ADRESSES DE DEPART
C     DES DEVELOPPEMENTS DES EFFORTS CONCERNES PAR LES BLOCAGES
C 
      CALL ADTBM( 'ADEP-DVBLO'  , ADNBLO )
C 
      CALL LONGEN( 'ADEP-DVBLO', LDVBLO )
C 
C     POUR TESTER SI UNE ADRESSE DE DEPART DANS ADEP-DVBLO
C     CORRESPOND AU DEVELOPPEMENT TESTE
C 
      DEBEFF = ( N+NTDSFG)*NDDL
      FINEFF =  (N+NTDSFG+1)*NDDL
C 
      NTDIAG = 0
C 
      CALL POUSME(NDDL , TERDIA )
      DBDIA = TERDIA
C 
       DO I = 0 , LDVBLO-1
C 
C     ON TESTE S'IL EXISTE DES NUMEROS DE DEVELOPPEMENT
C 
         IF( ( M( ADNBLO+I) .GT. DEBEFF)
     &       .AND.  ( M( ADNBLO+I) .LE. FINEFF) ) THEN
C 
CD         CALL IMPET( 'POUR DEBEFF ' ,DEBEFF  )
CD         CALL IMPET( 'POUR FINEFF ' ,FINEFF  )
D          CALL IMPET( 'POUR LA MATRICE N '//IDPROG , N )
CD         CALL IMPET( 'POUR LE TERME ADEP-DVBLO' , M(ADNBLO+I) )
C 
           M(DBDIA)  = M(ADNBLO+I)-DEBEFF
C 
D          CALL IMPET (' ON MODIFIE LE TERME DIAGONAL ',
D    &             M(DBDIA) )
C 
           DBDIA     = DBDIA+1
C 
           NTDIAG = NTDIAG+1
C 
         END IF
C 
      END DO
C 
      IF ( NTDIAG . GT . 0 ) THEN
C 
D      CALL IMPET( 'POUR LA MATRICE N ' , N )
D      CALL IMPET( 'NOMBRE DE TERME DIAGONAUX BLOQUES',
D    &   NTDIAG )
C 
        LONDIA = LONBLO+NTDIAG
        CALL POUSME( LONDIA , DBBLO )
        CALL COPITE( LONBLO , M( ADBLO) , M( DBBLO) )
        CALL COPITE( NTDIAG , M( TERDIA ) , M( DBBLO+LONBLO) )
C 
        CALL MODIMB( M(DBBLO) , LONDIA , M(ADPRO)
     &              , MATRIC , ADRESS, TERSUP)
C 
       ELSE
C 
         CALL MODIMB( M(ADBLO) , LONBLO , M(ADPRO)
     &               , MATRIC , ADRESS, TERSUP)
C 
       END IF
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
C     Modification de la diagonale d'une matrice symetrique stockee
C     profil en vue des deplacements imposes
C 
C     On envoie comme arguments :
C 
C     E ...... NUDDL    Numero des ddl concernes par les dep imposes
C     E ...... NBDDL    Nombre de ddl concernes
C     E ...... PROF     Tableau profil associe a MATPRO
C     E ...... MATRIC   La matrice  dont on veut modifier la diagonale
C     E ...... ADRESS   Adresse dans DM du debut de la matrice
C     E ...... TERSUP   Valeur sup de la diagonale

      SUBROUTINE MODIMB (NUDDL, NBDDL, PROF, MATRIC, ADRESS, TERSUP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER     PROF(NDDL+1), NUDDL(NBDDL), NBDDL, PLAMAT
      INTEGER     ADRESS
      DOUBLE PRECISION  MATRIC(NTMAT)  , TERSUP
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER                  I
      DOUBLE PRECISION         GROS  , DIAGO
      PARAMETER                (GROS=1.D15)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MODIMB')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     CALL IMPET('NTMAT DANS'//IDPROG, NTMAT)
C 
      IF( NTMAT+ADRESS . GT . LDMEFF) THEN
        CALL IMPET(' NTMAT+ADRESS ',NTMAT+ADRESS)
        CALL IMPET(' LDMEFF       ',LDMEFF)
        CALL ERREUD(0,' ECRASEMENT DE DM DANS'//IDPROG)
      ENDIF
C 
C      CALL IMPTET('NUMERO DES PIVOTS CONCERNES DANS '//IDPROG,
C    &    NUDDL(1) ,1 , NBDDL  )
C      CALL IMPTET('PROFIL DE LA MATRICE DANS '//IDPROG,
C     &    PROF(1) ,1 , (NBDDL+1)  )
C 
      DIAGO = TERSUP*GROS
C      CALL IMPDT('DIAGO DANS '//IDPROG, DIAGO)      
C 
      DO I= 1,NBDDL
C 
        PLAMAT            = PROF(NUDDL(I)+1)
        MATRIC( PLAMAT)   = DIAGO
C 
      ENDDO
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Mise a zero des efforts correspondant aux ddl a deplacement impose.
C     ===>  verification de l'admissibilite a zero en deplacement
C     Mise a zero des efforts non significatifs.
C     + 1er et dernier terme des developpements ?!
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT le nombre de fonctions du temps
C     ES...... EFFORT valeur des efforts en entree ranges (nddl, nfonct,nbmat)
C 
C     Et on recupere :
C 
C     ES...... EFFORT  valeur des efforts en sortie
C     S ...... NBDEV   nombre de developpements non mis a zero
C     S ...... TNUDEV  tableau des numeros de developpements non mis a zero
C     S ...... NORME   la norme de ce developpement

      SUBROUTINE MZBLOC (NFONCT, EFFORT, NBDEV, TNUDEV, NORME)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION EFFORT(NDDL*NBMAT), NORME
C 
      INTEGER          NBDEV, TNUDEV(NBMAT)
      INTEGER NFONCT
C 
C     pour nettoyer la partie des efforts sans signification
C 
      DOUBLE PRECISION SCAL, TEST
C 
      INTEGER          RESU, ARESU, I, DEBUT
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MZBLOC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL MZDDLB (1, EFFORT)
      CALL MZSYM (1, EFFORT)
      SCAL =  0.D0
C 
      DEBUT = NDDL+1
      CALL POUSMD (NBMAT, ARESU)
      RESU = ARESU
C 
C     mise a zero ds nudev  -ntdsfg
C 
      DM(RESU) = 0.D0
      RESU     = RESU+1
C 
      DO I = 2, NTDSFG
C 
        CALL SCAVEC (NDDL, EFFORT(DEBUT), EFFORT(DEBUT), DM(RESU))
C 
        DM(RESU) = PI*DM(RESU)
        SCAL = SCAL+ DM(RESU)
        RESU = RESU+1
        DEBUT = DEBUT + NDDL
C 
      END DO
C 
      CALL SCAVEC (NDDL, EFFORT(DEBUT), EFFORT(DEBUT), DM(RESU))
      DM(RESU) = 2.D0*PI*DM(RESU)
      SCAL  = SCAL+ DM(RESU)
      RESU  = RESU+1
      DEBUT = DEBUT + NDDL
C 
      DO I = NTDSFG+2, NBMAT-1
C 
        CALL SCAVEC (NDDL, EFFORT(DEBUT), EFFORT(DEBUT), DM(RESU))
        DM(RESU) = PI*DM(RESU)
        SCAL  = SCAL+ DM(RESU)
        RESU  = RESU+1
        DEBUT = DEBUT + NDDL
C 
      END DO
C 
C     mise a zero ds nudev  ntdsfg
C 
      DM(RESU) = 0.D0
C 
      TEST = 1.D-8*SCAL
C 
      DEBUT = 1
      NORME = 0.D0
C 
C     A EXPLIQUER CORRESPONDANCE AVEC LA MISE A ZERO DE -NTDSFG
C 
      NBDEV = 0
C 
      DO I  = ARESU, ARESU+NBMAT-1
C 
        IF (DM(I) .LT. TEST) THEN
C 
           CALL BALAID (NDDL, EFFORT(DEBUT))
C 
CD         NUDEV = I-ARESU-NTDSFG
C 
CD         CALL IMPET ('On a balaye le developpemet de l''effort'//
CD        ' numero dans '//idprog , NBDEV)
C 
CD         CALL IMPDT( ' Valeur relative ', DSQRT( DM(I) /SCAL) )
C 
        ELSE
C 
CD         NUDEV = I-ARESU-NTDSFG
C 
CD         CALL IMPET('Effort'//
CD         ' numero dans '//idprog , NUDEV)
C 
CD         CALL IMPDT( ' Valeur relative ', DSQRT( DM(I) /SCAL) )
C 
           NBDEV = NBDEV+1
C 
           TNUDEV(NBDEV) = I-ARESU+1
C 
           NORME = NORME+DM(I)
C 
        END IF
C 
        DEBUT = DEBUT+NDDL
C 
      END DO
C 
C 
C     Indication sur les residus extremes
C 
C     MODIFCD    CALL SCAVEC( NDDL , EFFORT(1) , EFFORT(1) , DM(ARESU) )
C     MODIFCD    DEBUT = 1+(NBMAT-1)*NDDL
C     MODIFCD    CALL SCAVEC( NDDL , EFFORT(DEBUT) , EFFORT(DEBUT) ,
C     MODIFCD                 DM(ARESU+NBMAT-1) )
C     MODIFCD    NORVER = NORME +DM(ARESU)+DM(ARESU+NBMAT-1)
C 
CD    RAP = DM(ARESU)/NORME
C 
CD    IF( RAP .GT. 0.2) THEN
C 
CD      CALL IMPDN('rapport du developpement -ntdsfg au total',
C 
C     MODIFCD     DM(ARESU)/NORVER )
C 
CD      RAP )
C 
CD      CALL IMPDT('valeur du developpement', DM(ARESU) )
C 
CD    END IF
C 
CD    RAP = DM(ARESU+NBMAT-1)/NORME
C 
CD    IF( RAP .GT. 0.2) THEN
C 
CD      CALL IMPDN('rapport du developpement ntdsfg au total',
C 
C     MODIFCD     DM(ARESU+NBMAT-1)/NORVER )
C 
CD      RAP )
C 
CD      CALL IMPDN('valeur du developpement', DM(ARESU+NBMAT-1) )
C 
CD    END IF
C 
CD     CALL IMPTDN ('VALEUR DES DEVELOPPEMENTS ',
CD     DM(ARESU),1,NBMAT)
C 
C     A EXPLIQUER  ET A AMELIORER
C 
C     MODIF       CALL BALAID(NDDL,EFFORT(1))
C     MODIF       DEBUT = 1+NDDL*(NBMAT-1)
C     MODIF       NBDEV = NBDEV+1
C     MODIF       TNUDEV(NBDEV)=NBMAT
C     MODIF       CALL BALAID(NDDL,EFFORT(DEBUT))
C 
      NORME = DSQRT( NORME)       
C 
C     CALL IMPTET ('Numero des efforts restants dans '//IDPROG,
C     &             TNUDEV(1), 1 , NBDEV)
C 
CD    IF ( LTRACP(1) ) THEN
C 
CD      CALL IMPTDP ('EFFORT EN SORTIE ',EFFORT(1), NDDL, NBMAT)
C 
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Mise a zero des efforts correspondant aux ddl a deplacements imposes
C     ===>  verification de l'admissibilite a zero en deplacement
C 
C     On envoie comme arguments :
C 
C     E ...... NFONCT le nombre de fonctions du temps
C     ES...... EFFORT valeur des efforts en entree ranges
C                     (nddl, nfonct, nbmat)
C 
C     Et on recupere :
C 
C     ES...... EFFORT valeur des efforts en sortie
C 
      SUBROUTINE MZDDLB( NFONCT , EFFORT )
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'strategie_calcul.h'
C 
      DOUBLE PRECISION EFFORT(NDDL*NBMAT*NFONCT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER   DDLBLO , NBTDDL , MAT , I , DEBUT  , NFONCT
C 
CD    LOGICAL LTRACP , LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='MZDDLB')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBM( 'DDL-BLOQUE',DDLBLO )
C 
      CALL LONGEN('DDL-BLOQUE',NBTDDL)
C 
CD    IF ( LTRACP(1) ) THEN
C 
CD      CALL IMPTEP('LISTE DES DDL BLOQUES ',M(DDLBLO),NBTDDL,1)
C 
CD      DEBUT = 1
C 
CD      DO I = 1 , NBMAT
C 
CD        CALL OMPTDP(' EFFORT EN ENTREE ',EFFORT(DEBUT),NDDL,NFONCT)
C 
CD      END DO
C 
CD      DEBUT = DEBUT+NDDL*NFONCT
C 
CD    END IF
C 
      DDLBLO = DDLBLO - 1
      DEBUT  = 0
C 
      DO MAT = 1 , NBMAT*NFONCT
C 
        DO I =  1 , NBTDDL
C 
          EFFORT( DEBUT + M( DDLBLO+I ) ) = 0.D0
C 
        END DO
C 
          DEBUT = DEBUT + NDDL
C 
      END DO
C 
      IF(FORS) CALL MZSYM(NFONCT, EFFORT)
C 
CD    IF ( LTRACP(1) ) THEN
C 
CD      CALL IMPTDP(' EFFORT EN SORTIE ', EFFORT(1), NDDL, NBMAT)
C 
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
