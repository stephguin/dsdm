C     Delete un fichier double precision non formate .
C 
C     On envoie comme arguments :
C 
C     E ...... NUMDIR          numero de la directory
C     E ...... NOMDIR          nom de la directory
C     E ...... NOMFIC(1:LG)    nom du fichier dans NOMDIR
C     E ...... LG              longueur de la chaine de caracteres 'nomdir'

      SUBROUTINE DELETF (NUMDIR, NOMDIR, NOMFIC, LG)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
C 
      INTEGER        NUMDIR, LG
      CHARACTER*40   NOMFIC
      CHARACTER*9    NOMDIR
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER      IERNAM
C 
      CHARACTER*6  IDPROG
      CHARACTER*12 FORM, DIR
      CHARACTER*20 RESU1, RESU2
      PARAMETER    (IDPROG='DELETF')
C 
      LOGICAL      EXT, OPE
      INTEGER      NUMFI, LRESU1, LRESU2, RECLEN
C 
C -----------------------------------------------------------------------
      PRINT*, 'NOM DE FICHIER ', NOMFIC(1:LG)
      PRINT*, 'NOM DE DIRECTORIE ', NOMDIR
C 
      INQUIRE (FILE=tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LG),
     &         EXIST=EXT)
C 
      IF (EXT) THEN
C 
        INQUIRE (FILE=tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LG),
     &           OPENED = OPE, NUMBER=NUMFI)
C 
        IF (OPE) THEN
C 
          CLOSE (NUMFI, ERR=900, STATUS='DELETE')
C 
          ELSE
C 
          INQUIRE (FILE=tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LG),
     &             DIRECT=DIR, FORMATTED =FORM, RECL= RECLEN)
C 
          IF (RECLEN .EQ. 0) THEN
C 
            RECLEN = 8
C 
          END IF
C 
          IF( DIR(1:1) .EQ. 'N') THEN
C 
            RESU1 = 'SEQUENTIAL'
            LRESU1 = 10
C 
          ELSE
C 
          RESU1 = 'DIRECT'
          LRESU1 = 6
C 
          END IF
C 
          IF (FORM(1:1) .EQ. 'Y') THEN
C 
            RESU2  = 'FORMATTED'
            LRESU2 = 9
C 
          ELSE
C 
            RESU2 = 'UNFORMATTED'
            LRESU2 = 11
C 
          END IF
C 
          OPEN (UNIT=99, FILE=tild(1:LTILD)//'/'//NOMDIR//'/'
     &         //NOMFIC(1:LG), ERR=1000)
C 
          CLOSE (99, ERR=1100, STATUS='DELETE')
C 
        END IF
C 
      END IF
C 
      RETURN
C 
1100  CONTINUE
      CALL IMPET ('NUMERO D''ERREUR ', IERNAM)
      CALL ERREUD (0, 'PB POUR DELETER LE FICHIER APRES OPEN'//IDPROG)
C 
1000  CONTINUE
      CALL IMPET ('NUMERO D''ERREUR ', IERNAM)
      CALL ERREUD (0, 'PB POUR DELETER LE FICHIER APRES OPEN '//IDPROG)
C 
900   CONTINUE
      CALL IMPET ('NUMERO D''ERREUR ', IERNAM)
      CALL ERREUD (0, 'PB POUR DELETER LE FICHIER '//IDPROG)
C 
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Creation de fichier double precision en direct non formatte.
C 
C     On envoie comme arguments :
C 
C     E........... NUMDIR    numero de la directory
C     E........... NOMDIR    nom de la directory
C     E........... NOMFIC    nom du fichier dans NOMDIR
C     E........... NBENRT    nombre total de doubles de l'enregistrement
C     E                      du fichier dans NOMDIR
C 
C     Et on recupere :
C 
C     S........... IUNIT    numero d'unite correspondant

      SUBROUTINE CFDDNF (NUMDIR, NOMDIR, NOMFIC, LGNMFI, NBENRT, IUNIT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER        NUMDIR, IUNIT, NBENRT, LGNMFI
      CHARACTER*40   NOMFIC
      CHARACTER*9    NOMDIR
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER        IERNAM
C 
      CHARACTER*6    IDPROG
      PARAMETER      (IDPROG='CFDDNF')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NBFICH( NUMDIR) = NBFICH( NUMDIR)+1
      CHAFIC( NBFICH(NUMDIR) , NUMDIR) = NOMFIC
C 
CD    CALL IMPCP(' NOM DE FICHIER ', NOMFIC )
C 
      CHADIR( NUMDIR) = NOMDIR
C 
CD    CALL IMPCP(' NOM DE DIRECTORY ', NOMDIR )
C 
      IF (NBFICH(NUMDIR) . GE . NBFIMX)THEN
        CALL ERREUD(0, 'NUFICH > NBFIMX DANS '//IDPROG)
      END IF
C 
CD    CALL IMPEP(' NOMBRE DE FICHIER DE NOMDIR',NBFICH(NUMDIR) )
C 
      CALL NUNFOU( IUNIT )
C 
CD    CALL IMPEP(' NUMERO D''UNITE',IUNIT )
C 
        OPEN (UNIT=IUNIT
     &  ,FILE=tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LGNMFI)
     &            ,FORM='UNFORMATTED' ,RECL = NBENRT*8
     &            ,ACCESS='DIRECT',IOSTAT=IERNAM , ERR=900 )
C 
CD     CALL IMPCN('NOM DE FICHIER    '//idprog, NOMFIC(1:LGNMFI))
CD     CALL IMPCN(' NOM DE DIRECTORY  '//idprog, NOMDIR )
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
900   CONTINUE
       CALL IMPET('POUR LE NUMERO DE FICHIER ', NBFICH(NUMDIR) )
       CALL IMPCT('NOM DE FICHIER    '//idprog, NOMFIC(1:LGNMFI))
       CALL IMPCT(' NOM DE DIRECTORY  '//idprog, NOMDIR )
       CALL IMPCT(' file '//idprog,
     &  tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LGNMFI) )
       CALL ERROR_$PRINT(IERNAM)
       CALL IMPET(' NUMERO D''ERREUR ', IERNAM )
       CALL ERREUD (0,'PB D''OUVERTURE DANS'//IDPROG)
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Ouverture de fichier double precision en direct non formatte
C 
C     On envoie comme arguments :
C 
C     E ...... NUMDIR    numero de la directory
C     E ...... NOMFIC    nom du fichier dans NOMDIR
C     E ...... NBENRT    nombre total d'enregistrements
C                        du fichier dans NOMDIR
C     Et on recupere :
C 
C     S ...... IUNIT     numero d'unite correspondant

      SUBROUTINE OFDDNF (NUMDIR, NOMFIC, LGNMFI, NBENRT, IUNIT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER        NUMDIR, IUNIT, NBENRT, LGNMFI
      CHARACTER*40   NOMFIC
      CHARACTER*9    NOMDIR
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER        IERNAM
C 
      CHARACTER*6    IDPROG
      PARAMETER      (IDPROG='OFDDNF')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NOMDIR = CHADIR(NUMDIR)
C 
CD    CALL IMPCP(' NOM DE DIRECTORIE', NOMDIR)
C 
C     CALL NUFICH (NUMDIR, NOMFIC, NUMEFI)
C 
      IF (NBFICH(NUMDIR) . GE . NBFIMX)THEN
        CALL ERREUD(0, 'NUFICH > NBFIMX DANS '//IDPROG)
      END IF
C 
      CALL NUNFOU( IUNIT )
C 
CD    CALL IMPEP(' NUMERO D''UNITE',IUNIT)
C 
C     On ferme avant d'ouvrir pour supprimer probleme dans OPEN
C 
        CLOSE(UNIT=IUNIT)
        OPEN (UNIT=IUNIT
     &  ,FILE=tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LGNMFI)
     &            ,FORM='UNFORMATTED' ,RECL = 8*NBENRT
     &            ,ACCESS='DIRECT',IOSTAT=IERNAM , ERR=900 )
C 
CD     CALL IMPCN('NOM DE FICHIER    '//idprog, NOMFIC(1:LGNMFI))
CD     CALL IMPCN(' NOM DE DIRECTORY  '//idprog, NOMDIR )
CD     CALL IMPEN(' NUMERO D''UNITE '//idprog, IUNIT )
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
900   CONTINUE
       CALL IMPCT('NOM DE FICHIER    '//idprog, NOMFIC(1:LGNMFI))
       CALL IMPCT(' NOM DE DIRECTORY  '//idprog, NOMDIR )
       CALL IMPET(' NUMERO D''UNITE '//idprog, IUNIT )
       CALL ERROR_$PRINT(IERNAM)
       CALL IMPET(' NUMERO D''ERREUR ', IERNAM )
       CALL ERREUD (0,'PB D''OUVERTURE DANS'//IDPROG)
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Ecriture de fichier double precision en direct non formatte
C 
C     On envoie comme arguments :
C 
C     E ...... IUNIT    LE NUMERO D'UNITE DU FICHIER
C     E ...... TAB      LE TABLEAU DOUBLE PRECISION
C     E ...... ADRESS   L''ADRESSE DU DEBUT DU TABLEAU DANS DM
C     E ...... LONTAB   SA LONGUEUR
C     E ...... NUENRS   Le numero d'enregistrement
C 
      SUBROUTINE EFDDNF (IUNIT, TAB, ADRESS, LONTAB, NUENRS)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           LONTAB, IUNIT, ADRESS, NUENRS
C 
      DOUBLE PRECISION  TAB(LONTAB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EFDDNF')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF( ADRESS+LONTAB . GT . LDMEFF   ) THEN
CD      CALL ERREUD(0,'LECTURE EN DEHORS DE DM DANS '//IDPROG)
CD    END IF
C 
      WRITE(IUNIT,REC=NUENRS)TAB
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Lecture de fichier double precision en direct non formatte
C 
C     Permet d'affecter dans tab le contenu de l'unite IUNIT prealablement
C     ouverte et correspondant a un stockage direct
C 
C     On envoie comme arguments :
C 
C     E ...... TAB      LE TABLEAU DOUBLE PRECISION QUE L'ON REMPLIT
C     E ...... ADRESS   L''ADRESSE DU DEBUT DU TABLEAU DANS DM
C     E ...... LONTAB   SA LONGUEUR
C     E ...... IUNIT    Le numero d'unite
C     E ...... NUENRS   Le numero d'enregistrement
C 
      SUBROUTINE LFDDNF (TAB, ADRESS, LONTAB, IUNIT, NUENRS)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER             LONTAB, IUNIT, ADRESS, NUENRS
      DOUBLE PRECISION    TAB(LONTAB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='LFDDNF')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      IF (ADRESS+LONTAB .GT. LDMEFF) THEN
        CALL ERREUD (0, 'ECRASEMENT DE TABLEAU DANS '//IDPROG)
      END IF
C 
      READ (IUNIT, REC=NUENRS) TAB
C  
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
