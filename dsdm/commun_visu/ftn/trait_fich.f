C    Cette routine renvoie un numero d'unite non ouverte pour toute
C    ouverture de nouveau fichier
C 
C    ...ON ENVOIE COMME ARGUMENTS
C    ET ON RECUPERE... N

      SUBROUTINE NUNFOU(N)
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
      INTEGER N
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NUNFOU')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF ( NBFIOU . GT . 0 ) THEN
C 
        NBFIOU         = NBFIOU + 1
        NUFIOU(NBFIOU) = NUFIOU(NBFIOU-1)+1
        N              = NUFIOU(NBFIOU)
C 
      ELSE
C 
C     5 = >  clavier
C     6 = >  ecran
C     7 = >  sauvegarde
C 
        NBFIOU         = 1
        NUFIOU(NBFIOU) = 8
        N              = NUFIOU(NBFIOU)
C 
      END IF
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine retire de la liste des unites ouvertes le numero
C     d'unite N, decremente le nombre d'unites ouvertes et actualise
C     les numeros d'unites ouvertes
C 
C     On envoie comme arguments :
C 
C     Et on recupere :
C 
C     S ...... NUMDIR
C     S ...... N
C     S ...... NOMSUB
C 
      SUBROUTINE FERFIC( NUMDIR  , N , NOMSUB )
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
      INTEGER N , NUMDIR
C 
      CHARACTER*6 NOMSUB
C 
      INTEGER I , PLAC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='FERFIC')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     Recherche de l'adresse du numero d'unite dans NUFIOU
C     dans l'ordre inverse (a priori plus rapide ?)
C 
      DO I =  NBFIOU , 1 , -1
C 
        IF (NUFIOU(I) . EQ . N) THEN
C 
          PLAC = I
          GOTO 1
C 
        END IF
C 
      END DO
C 
C     Si le numero d'unite n'existe pas dans NUFIOU
C 
      CALL IMPET ('LE NUMERO D''UNITE ', N)
      CALL ERREUD (0, 'N''EST PAS OUVERT DANS '//NOMSUB)
C 
1     CONTINUE
C  
      CLOSE(N)
C 
CD    CALL IMPEN ('FERFIC : NUMERO D''UNITE FERMEE DANS '//NOMSUB,N )
C 
      IF (PLAC .EQ. NBFIOU) THEN
C 
        NUFIOU( NBFIOU) = 0
C 
      ELSE
C 
        DO I = PLAC+1 , NBFIOU
C 
          NUFIOU(I-1) = NUFIOU(I)
C 
        END DO
C 
      END IF
C 
      NBFIOU = NBFIOU - 1
C 
C     Decrementation du nombre de fichiers ouverts dans la
C     directory de numero NUMDIR
C 
      NBFICH(NUMDIR) = NBFICH(NUMDIR)-1
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
         SUBROUTINE FERFID( NUMDIR  , N , NOMSUB )
C -
C - Cette routine retire de la liste des unites  ouvertes le numero
C - d'unite N et decremente le nombre d'unite ouverte etactualise
C - les numeros d'unites ouvertese
C -
C       ...ON ENVOIE COMME ARGUMENTS
C          ET ON RECUPERE... NUMDIR , N , NOMSUB
C
C[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]]][]
C
C*      DECLARATION DES PARAMETRES GLOBAUX
C       """"""""""""""""""""""""""""""""""
C
      include 'cominc.h'
C
C**********************************************************************
C
C*    DECLARATION DES PARAMETRES LOCAUX
C     """""""""""""""""""""""""""""""""
C
      INTEGER N , NUMDIR
C
      CHARACTER*6 NOMSUB
C
      INTEGER I , PLAC , IERNAM
C
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='FERFID')
C
CD    CALL WLKBCD(IDPROG)
C
C     Recherche de l'adresse du numero d'unite dans NUFIOU
C     dans l'ordre inverse ( a priori plus rapide ?! )
C
      DO I =  NBFIOU , 1 , -1
C -
        IF (  NUFIOU(I) . EQ . N ) THEN
C -
          PLAC = I
          GOTO 1
C -
        END IF
C -
      END DO
C -
C     Si le numero d'unite n'existe pas dans NUFIOU
C -
      CALL IMPET ( ' LE NUMERO D''UNITE ' , N )
      CALL ERREUD ( 0 , ' N''EST PAS OUVERT DANS '//NOMSUB )
C -
1     CONTINUE
C -
      CLOSE(N, ERR=1100,STATUS='DELETE')
C
CD    CALL IMPEN('FERFIC :NUMERO D''UNITE FERMEE DANS '//NOMSUB,N )
C -
      IF ( PLAC .EQ. NBFIOU ) THEN
C -
        NUFIOU( NBFIOU) = 0
C -
      ELSE
C -
        DO I = PLAC+1 , NBFIOU
C -
          NUFIOU( I-1) = NUFIOU(I)
C -
        END DO
C -
      END IF
C -
      NBFIOU = NBFIOU - 1
C -
C -     Decrementation du nombre de fichiers ouverts dans la
C -   directory de numero NUMDIR .
C -
      NBFICH( NUMDIR ) = NBFICH( NUMDIR)-1
C -
CD    CALL RETOUD(IDPROG)
C
      RETURN
1100  CALL IMPET(' NUMERO D''ERREUR dans '//IDPROG, IERNAM )
      CALL ERREUD(0,
     &  'PB POUR DELETER LE FICHIER  DANS'//IDPROG)
C -
      END

