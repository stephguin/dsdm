C     Ecriture des fichiers pour visumail. Les fichiers sont formattes.
C  
C     On envoie comme arguments :
C 
C     E ...... IUNIT    LE NUMERO D'UNITE DU FICHIER
C     E ...... TAB      LE TABLEAU DOUBLE PRECISION
C     E ...... TABENT   LE TABLEAU ENTIER
C     E ...... ADRESS   L'ADRESSE DU DEBUT DU TABLEAU DANS DM
C     E ...... LONTAB   SA LONGUEUR
C 
      SUBROUTINE EFVISM (IUNIT, TAB, TABENT, ADRESS, LONTAB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER       LONTAB, IUNIT, I, ADRESS
      INTEGER       STATUS, TABENT(LONTAB)
C 
      DOUBLE PRECISION    TAB(LONTAB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EFVISM')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (ADRESS+LONTAB .GT. LDMEFF) THEN
CD      CALL ERREUD (0,'LECTURE EN DEHORS DE DM DANS '// IDPROG)
CD    END IF
CD    CALL IMPEN ('NUMERO D''UNITE', IUNIT)
C 
      CALL TEST_FORMAT (LONTAB, TAB)
      WRITE (IUNIT, FMT=15, ERR=99000, IOSTAT=STATUS)
     &      (TABENT(I), TAB(I), I=1, LONTAB)
15    FORMAT(I4,',',E12.6)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
C 
99000 CONTINUE
      CALL IMPCT ('NOM DE ROUTINE ', IDPROG)
      CALL IMPET ('NUMERO D''UNITE', IUNIT)
      CALL IMPET ('STATUS ', STATUS)
      CALL ERROR_$PRINT (STATUS)
      WRITE (6, 99100) IUNIT, STATUS
99100 FORMAT (/,'PROBLEME ECRITURE SUR UNITE : ', I10,/,
     &         ' STATUS : ', I10,/)
C 
      CALL ERREUD (0,' BOFFFFFFFF ???')
C 
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... ACCES     sequentiel TYPE unformatted
C     E ...... NUMDIR    numero de la directory
C     E ...... NOMDIR    nom de la directory
C     E ...... NOMFIC    nom du fichier dans NOMDIR
C     E ...... LGNMFI    longueur du nom du fichier dans NOMDIR
C     E ...... TYPE      (F ou U) Formatted ou Unformatted
C     E ...... NBENRT    nombre total de double du fichier dans NOMDIR
C 
C     Et on recupere :
C 
C     S....... IUNIT     numero d'unite correspondant
C 
         SUBROUTINE CREFIC (NUMDIR, NOMDIR, NOMFIC,
     &                       LGNMFI, TYPE, NBENRT, IUNIT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER        NUMDIR  , IUNIT , LGNMFI , NBENRT
      CHARACTER*1    TYPE
      CHARACTER*40   NOMFIC
      CHARACTER*9    NOMDIR
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         IERNAM
C 
      CHARACTER*6 IDPROG
C      CHARACTER*30  TILDLO
      PARAMETER (IDPROG='CREFIC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NBFICH(NUMDIR) = NBFICH (NUMDIR)+1
      LONNOM(NBFICH(NUMDIR), NUMDIR) = LGNMFI
      CHAFIC(NBFICH(NUMDIR), NUMDIR) = NOMFIC
C 
CD    CALL IMPCN('NOM DE FICHIER DANS '//IDPROG, NOMFIC(1:LGNMFI))
C 
      CHADIR(NUMDIR) = NOMDIR
C 
CD    CALL IMPEN('NUMERO DE DIRECTORY '//IDPROG, NUMDIR)
CD    CALL IMPCN('   NOM DE DIRECTORY '//IDPROG, NOMDIR)
C 
      IF (NBFICH(NUMDIR) .GE. NBFIMX) THEN
        CALL ERREUD (0, 'NUFICH > NBFIMX DANS '//IDPROG)
      END IF
C 
CD    CALL IMPEP('NOMBRE DE FICHIER DE NOMDIR ', NBFICH(NUMDIR) )
C 
      CALL NUNFOU (IUNIT)
CD     CALL IMPEN('NUMERO D''UNITE ',IUNIT)
C 
       IF (TYPE .EQ. 'U') THEN
         OPEN (UNIT=IUNIT,
     &         FILE=tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LGNMFI),
     &         FORM='UNFORMATTED' , RECL = 8*NBENRT,
     &         ACCESS='SEQUENTIAL', IOSTAT=IERNAM , ERR=900 )
       ELSE IF (TYPE .EQ. 'F') THEN
         OPEN (UNIT=IUNIT,
     &         FILE=tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LGNMFI),
     &         FORM='FORMATTED' , RECL = 8*NBENRT,
     &         ACCESS='SEQUENTIAL', IOSTAT=IERNAM, ERR=900)
       ELSE
         CALL IMPET ('POUR LE NUMERO DE FICHIER ', NBFICH(NUMDIR))
         CALL IMPCT ('POUR LE TYPE  ', TYPE)
         CALL IMPCT ('POUR LE NOM DE FICHIER ', NOMFIC(1:LGNMFI))
         CALL IMPCT ('POUR LE NOM DE DIRECTORY ', NOMDIR )
         CALL ERREUD (0, 'MAUVAIS PASSAGE DE TYPE DANS '//IDPROG)
       END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
900   CONTINUE
      CALL IMPCT ('NOM DE ROUTINE' , IDPROG )
      CALL IMPET ('POUR LE NUMERO DE FICHIER ', NBFICH(NUMDIR) )
      CALL IMPET ('POUR LE NUMERO D''UNITE ', IUNIT )
      CALL IMPET ('POUR LE NUMERO DE FICHIER ', NBFICH(NUMDIR) )
      CALL IMPTET ('NUMERO DE FICHIERS OUVERTS ', NUFIOU(1), NBFIOU, 1)
      CALL ERROR_$PRINT (IERNAM)
C 
      CALL IMPET ('NUMERO D''ERREUR ', IERNAM )
      CALL IMPCT ('NOM COMPLET DE FICHIER ',
     &             tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LGNMFI))
      CALL ERREUD (0,'PROBLEME D''OUVERTURE DANS '//IDPROG)
C 
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C     On envoie comme arguments:
C 
C     ACCES sequentiel TYPE unformatted
C     E................ NUMDIR    numero de la directory
C     E................ NOMFIC    nom du fichier dans NOMDIR
C     E................ LGNMFI    longueur du nom du fichier dans NOMDIR
C     E................ TYPE      (F ou U) Formatted ou Unformatted
C     E................ NBENRT    nombre d'enregistrement
C 
C     On recupere :
C 
C     S...............  IUNIT    numero d'unite correspondant
C 
      SUBROUTINE OUVFCD(NUMDIR, NOMFIC, LGNMFI, TYPE, NBENRT,
     &                  IUNIT )
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER        NUMDIR, IUNIT, LGNMFI, NBENRT
      CHARACTER*1    TYPE
      CHARACTER*40   NOMFIC
      CHARACTER*120  NOM
      CHARACTER*9    NOMDIR
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER        IERNAM,LNOM
C 
      CHARACTER*6    IDPROG
      PARAMETER      (IDPROG='OUVFCD')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NOMDIR = CHADIR(NUMDIR)
C 
CD    CALL IMPEN('NUMERO DE DIRECTORY DANS '//IDPROG, NUMDIR)
C 
      IF (NBFICH(NUMDIR) . GE . NBFIMX) THEN
        CALL ERREUD(0, 'NUFICH > NBFIMX DANS '//IDPROG)
      END IF
C 
      CALL NUNFOU( IUNIT )
C 
CD    CALL IMPEN('NUMERO D''UNITE ', IUNIT)
      NOM = tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LGNMFI)
      LNOM = LTILD +LGNMFI+9+1+1
C 
CD    CALL IMPCT('POUR PATHNAME ', NOM )
C 
      IF ( TYPE .EQ. 'U' ) THEN
         OPEN (UNIT=IUNIT,
     &         FILE=NOM(1:LNOM)
     &         ,FORM='FORMATTED' ,
C          RECL = 8*NBENRT
     &          ACCESS='SEQUENTIAL',IOSTAT=IERNAM , ERR=900 )
       ELSE IF ( TYPE . EQ . 'F') THEN
         OPEN (UNIT=IUNIT,
     &         FILE=NOM(1:LNOM)
     &         ,FORM='FORMATTED' ,RECL = 8*NBENRT
     &         ,ACCESS='SEQUENTIAL',IOSTAT=IERNAM , ERR=900 )
            ELSE
              CALL IMPET('POUR LE NUMERO DE FICHIER ', NBFICH(NUMDIR) )
              CALL IMPCT('POUR LE TYPE  ', TYPE  )
              CALL IMPCT('POUR LE NOM DE FICHIER ', NOMFIC(1:LGNMFI) )
              CALL IMPCT('POUR LE NOM DE DIRECTORY ', NOMDIR )
              CALL ERREUD(0, 'MAUVAIS PASSAGE DE TYPE DANS '//IDPROG)
      END IF
C 
CD    INQUIRE (UNIT=IUNIT,NAME=NOM)
CD    CALL IMPCT(' NOM PAR INUIRE '//IDPROG, NOM)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
900   CONTINUE
      CALL IMPCT('NOM DE ROUTINE ', IDPROG)
      CALL IMPCT('POUR LE NOM DE FICHIER ', NOMFIC(1:LGNMFI))
      CALL IMPCT('POUR LE NOM DE DIRECTORY ', NOMDIR)
      CALL IMPCT('POUR PATHNAME ',
     &            tild(1:LTILD)//'/'//NOMDIR//'/'//NOMFIC(1:LGNMFI) )
      CALL IMPET(' NUMERO D''ERREUR ', IERNAM )
      CALL ERROR_$PRINT(IERNAM)
      CALL ERREUD (0,'PB D''OUVERTURE DANS'//IDPROG)
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... IUNIT    LE NUMERO D'UNITE DU FICHIER
C     E ...... TAB      LE TABLEAU DOUBLE PRECISION
C     E ...... ADRESS   L''ADRESSE DU DEBUT DU TABLEAU DANS DM
C     E ...... LONTAB   SA LONGUEUR

      SUBROUTINE ECFICD (IUNIT, TAB, ADRESS, LONTAB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER       LONTAB, IUNIT, ADRESS
      INTEGER       STATUS
C 
      DOUBLE PRECISION    TAB(LONTAB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='ECFICD')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (ADRESS+LONTAB .GT. LDMEFF) THEN
CD      CALL ERREUD (0, 'LECTURE EN DEHORS DE DM DANS '//IDPROG)
CD    END IF
CD    CALL IMPEN ('numero d''unite ', IUNIT)
C 
      WRITE (IUNIT, ERR=99000, IOSTAT=STATUS) TAB
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
99000 CONTINUE
      CALL IMPCT ('NOM DE ROUTINE ', IDPROG)
      CALL IMPET ('numero d''unite ', IUNIT)
      CALL IMPET ('Status ', STATUS)
      CALL ERROR_$PRINT (STATUS)
      WRITE (6,99100) IUNIT, STATUS
99100 FORMAT(/,'PROBLEME ECRITURE SUR UNITE : ',I10,/,
     &         'STATUS : ', I10,/)
C 
      CALL ERREUD (0, ' BOFFFFFFFF ???')
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     On envoie comme arguments :
C 
C     E ...... IUNIT    LE NUMERO D'UNITE DU FICHIER
C     E ...... TAB      LE TABLEAU DOUBLE PRECISION
C     E ...... ADRESS   L'ADRESSE DU DEBUT DU TABLEAU DANS M
C     E ...... LONTAB   SA LONGUEUR

      SUBROUTINE ECFICE (IUNIT, TAB, ADRESS, LONTAB)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER       LONTAB, IUNIT, ADRESS
      INTEGER       STATUS
C 
      INTEGER       TAB(LONTAB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='ECFICE')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (ADRESS+LONTAB .GT. LM) THEN
        CALL ERREUD (0,'LECTURE EN DEHORS DE M DANS '//IDPROG)
      END IF
C 
CD    CALL IMPEN ('numero d''unite ', IUNIT)
C 
      WRITE (IUNIT, ERR=99000, IOSTAT=STATUS) TAB
C 
C     CLOSE (IUNIT)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
99000 CONTINUE
      CALL IMPCT ('NOM DE ROUTINE ', IDPROG)
      CALL IMPET ('numero d''unite ', IUNIT)
      CALL IMPET ('Status ', STATUS)
      CALL ERROR_$PRINT (STATUS)
      WRITE (6,99100) IUNIT, STATUS
99100 FORMAT (/,'PROBLEME ECRITURE SUR UNITE : ', I10,/,
     &          ' STATUS : ', I10,/)
C 
      CALL ERREUD (0, 'BOFFFFFFFF ???')
C 
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Permet d'affecter dans tab le contenu de l'unite IUNIT prealablement
C     ouverte et correspondant a un stockage sequentiel.
C 
C     On envoie comme arguments :
C 
C     E ...... TAB      LE TABLEAU DOUBLE PRECISION QUE L'ON REMPLIT
C     E ...... ADRESS   L'ADRESSE DU DEBUT DU TABLEAU DANS DM
C     E ...... LONTAB   SA LONGUEUR
C     E ...... IUNIT    le numero d'unite
C 
      SUBROUTINE LIFCDU (TAB, ADRESS, LONTAB, IUNIT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER       LONTAB, IUNIT,  ADRESS
      INTEGER       STATUS
C 
      DOUBLE PRECISION    TAB(LONTAB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='LIFCDU')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (ADRESS+LONTAB .GT. LDMEFF) THEN
        CALL ERREUD (0, 'ECRASEMENT DE TABLEAU DANS '//IDPROG)
      END IF
C 
CD    CALL IMPEN ('NUMERO D''UNITE ' ,IUNIT)
C 
      READ (IUNIT, ERR=99000, IOSTAT=STATUS, END =400) TAB
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
99000 CONTINUE
      CALL IMPCT ('NOM DE ROUTINE ', IDPROG)
      CALL PERROR ('LIFCDU')
      WRITE (6,99100) IUNIT, STATUS
99100 FORMAT(/,'PROBLEME ECRITURE SUR UNITE : ', I10, /,
     &         'STATUS : ', I10, /)
C 
      CALL ERREUD (0, 'PAS DE FIN DE FICHIER DANS '//IDPROG)
C 
400   CONTINUE
      CALL ERREUD (0, 'FIN DE FICHIER DANS '//IDPROG)
C 
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Permet d'affecter dans tab le contenu de l'unite IUNIT prealablement
C     ouverte et correspondant a un stockage sequentiel.
C 
C     On envoie comme arguments :
C 
C     E ...... TAB      LE TABLEAU DOUBLE PRECISION QUE L'ON REMPLIT
C     E ...... ADRESS   L'ADRESSE DU DEBUT DU TABLEAU DANS M
C     E ...... LONTAB   SA LONGUEUR
C     E ...... IUNIT    le numero d'unite
C 
         SUBROUTINE LIFCEU (TAB, ADRESS, LONTAB, IUNIT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER   LONTAB, IUNIT, ADRESS
C 
      INTEGER   TAB(LONTAB)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='LIFCEU')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (ADRESS+LONTAB .GT. LM) THEN
        CALL ERREUD (0, 'ECRASEMENT DE M DANS '//IDPROG)
      END IF
C 
CD    CALL IMPEN ('NUMERO D''UNITE' ,IUNIT)
C 
      READ (IUNIT, ERR=800, END =400) TAB
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
800   CONTINUE
      CALL IMPCT ('NOM DE ROUTINE ', IDPROG)
      CALL ERREUD (0, 'PB DE LECTURE DANS '//IDPROG)
C 
400   CONTINUE
      CALL IMPCT ('NOM DE ROUTINE ', IDPROG)
      CALL ERREUD (0, 'FIN DE FICHIER DANS '//IDPROG)
C 
      END
