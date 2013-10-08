C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
C
C
C     L.M.T.------------CACHAN------------------------------
C
C     version du 23 / 09 / 1986
C
C     Jean-Yves COGNARD, Michel POSS
C
C     ensembles de sous programmes pour les indications  des traces
C     normales
C
C     traces s'effectuant sur NFIMPR et NFECRA
C     ESSAI D'IMPRESSION OLIVIER
C     ..................................................................
C
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
      SUBROUTINE TRAAON
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
C
C     ne sert a rien
C
C     ..................................................................
C
C *   baratin @ imprimer
      CHARACTER NOM *(*)
C *
      include 'comm_trace_inter.h'
c
      include 'files_debint_inc.h'
c
C *   entier , r{el , tableau @ imprimer (nombre de lignes et colonnes)
      INTEGER NL , NC
      DOUBLE PRECISION  TABLO(NL,NC)
      INTEGER ITABLO(NL,NC)
      INTEGER NB , NN , I
      CHARACTER*16 FORMAT
C
      INTEGER II , JJ
C
      RETURN
C
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
      ENTRY OMPTDN ( NOM , TABLO , NL , NC )
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
C
C     FONCTION :
C     impression transpose d'un tableau r{el double precision
C
C
C     indique le numero de la ligne (ou de la colonne suivant le type
C     d'impression) si le nombre de lignes du tableau est sup{rieur
C     @ 6
C
C     traces normales
C
C     Exemple : CALL IMPTDN ( 'contraintes',CONTRA,1,6)
C
C     R{sultat: - TDN - , TABLEAU : contraintes  1 ligne(s) 6 colonne(s)
C
C     Exemple : CALL IMPTDN ( 'XXXXX',XXXXXX,10,6)
C
C     R{sultat:
C        - TDN - , TABLEAU : XXXXXX (transpos{) 10 ligne(s) 6 colonne(s)
C
C     ..................................................................
C
C ... E    CHARACT  NOM    : Nom du tableau @ imprimer
C ... E    DOU-PRE  TABLO  : Tableau @ imprimer
C ... E    INTEGER  NL     : Nombre de lignes du tableau
C ... E    INTEGER  NC     : Nombre de colonnes du tableau
C
C     ..................................................................
C -
C -   test d'impression
      IF (LN) THEN
C -
C -     test de duplication des traces sur l'ecran
        IF (LS) THEN
C -
           WRITE(NFIMPR,21)NOM,NL,NC
           WRITE(NFECRA,21)NOM,NL,NC
           DO II=1,NC
             IF (NC.GT.2) THEN
               WRITE (NFIMPR,11) II
               WRITE (NFECRA,11) II
             END IF
             WRITE (NFIMPR,22) (TABLO(JJ,II),JJ=1,NL)
             WRITE (NFECRA,22) (TABLO(JJ,II),JJ=1,NL)
           ENDDO
         ELSE
C -
           WRITE(NFIMPR,21)NOM,NL,NC
           DO II=1,NC
             IF (NC.GT.2) THEN
               WRITE (NFIMPR,11) II
             END IF
             WRITE (NFIMPR,22) (TABLO(JJ,II),JJ=1,NL)
           ENDDO
         END IF
      END IF
11    FORMAT(' ligne numero : ',I5)
21    FORMAT(/,' * TDN * , TABLEAU :',A ,
     &       I5,' ligne(s) ',I5,' colonne(s) ')
22    FORMAT(10D13.5)
      RETURN
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
      ENTRY OMPTDT ( NOM , TABLO , NL , NC )
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
C
C     FONCTION :
C     impression transpose d'un tableau r{el double precision
C
C
C     indique le numero de la ligne (ou de la colonne suivant le type
C     d'impression) si le nombre de lignes du tableau est sup{rieur
C     @ 6
C
C     traces normales
C
C     Exemple : CALL IMPTDN ( 'contraintes',CONTRA,1,6)
C
C     R{sultat: - TDN - , TABLEAU : contraintes  1 ligne(s) 6 colonne(s)
C
C     Exemple : CALL IMPTDN ( 'XXXXX',XXXXXX,10,6)
C
C     R{sultat:
C        - TDN - , TABLEAU : XXXXXX (transpos{) 10 ligne(s) 6 colonne(s)
C
C     ..................................................................
C
C ... E    CHARACT  NOM    : Nom du tableau @ imprimer
C ... E    DOU-PRE  TABLO  : Tableau @ imprimer
C ... E    INTEGER  NL     : Nombre de lignes du tableau
C ... E    INTEGER  NC     : Nombre de colonnes du tableau
C
C     ..................................................................
C -
C -   test d'impression
C -     test de duplication des traces sur l'ecran
        IF (LS) THEN
C -
           WRITE(NFIMPR,21)NOM,NL,NC
           WRITE(NFECRA,21)NOM,NL,NC
           DO II=1,NC
             IF (NC.GT.2) THEN
               WRITE (NFIMPR,11) II
               WRITE (NFECRA,11) II
             END IF
             WRITE (NFIMPR,22) (TABLO(JJ,II),JJ=1,NL)
             WRITE (NFECRA,22) (TABLO(JJ,II),JJ=1,NL)
           ENDDO
         ELSE
C -
           WRITE(NFIMPR,21)NOM,NL,NC
           DO II=1,NC
             IF (NC.GT.2) THEN
               WRITE (NFIMPR,11) II
             END IF
             WRITE (NFIMPR,22) (TABLO(JJ,II),JJ=1,NL)
           ENDDO
         END IF
c11    FORMAT(' ligne numero : ',I5)
c21    FORMAT(/,' * TDN * , TABLEAU :',A ,
c     &       I5,' ligne(s) ',I5,' colonne(s) ')
c22    FORMAT(10D13.5)
      RETURN

C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
      ENTRY OMPTEN ( NOM , ITABLO , NL , NC )
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
C
C     FONCTION :
C     impression d'un tableau entier
C
C     meme r}gle d'impression que pour IMPTDN
C
C     traces normales
C
C     Exemple : CALL IMPTEN ( 'indicateurs',INDICA,1,6)
C
C     R{sultat: - TEN - , TABLEAU : indicateurs 1 ligne(s) 6 colonne(s)
C
C     ..................................................................
C
C ... E    CHARACT  NOM    : Nom du tableau @ imprimer
C ... E    INTEGER  ITABLO : Tableau @ imprimer
C ... E    INTEGER  NL     : Nombre de lignes du tableau
C ... E    INTEGER  NC     : Nombre de colonnes du tableau
C
C     ..................................................................
C -
C -   test d'impression
      IF (LN) THEN
C -
C -     test de duplication des traces sur l'ecran
        IF (LS) THEN
C -
C -       recherche du maxi entre le nombre de lignes et de colonnes
            WRITE(NFIMPR,41)NOM,NL,NC
            WRITE(NFECRA,41)NOM,NL,NC
            DO II=1,NC
              WRITE (NFIMPR,42) (ITABLO(JJ,II),JJ=1,NL)
              WRITE (NFECRA,42) (ITABLO(JJ,II),JJ=1,NL)
            ENDDO
        ELSE
            WRITE(NFIMPR,41)NOM,NL,NC
            DO II=1,NC
              WRITE (NFIMPR,42) (ITABLO(JJ,II),JJ=1,NL)
            ENDDO
        ENDIF
      ENDIF
C -
41    FORMAT(/,' * TEN * , TABLEAU :',A ,
     &       I5,' ligne(s) ',I5,' colonne(s) ')
42    FORMAT(4(5 I6 , 2X))
      RETURN
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     FONCTION : impression transpose d'un tableau reel double precision
C 
C     Indique le numero de la ligne (ou de la colonne suivant le type
C     d'impression) si le nombre de lignes du tableau est superieur a 6
C 
C     Traces normales
C 
C     Exemple : CALL IMPTDN ('contraintes', CONTRA, 1, 6)
C 
C     Resultat: - TDN - , TABLEAU : contraintes  1 ligne(s) 6 colonne(s)
C 
C     Exemple : CALL IMPTDP ('XXXXX', XXXXXX, 10, 6)
C 
C     Resultat: - TDN - , TABLEAU : XXXXXX (transpose 10 ligne(s) 6 colonne(s))
C 
C     On envoie comme arguments :
C 
C     E ...... CHARACT  NOM    : Nom du tableau a imprimer
C     E ...... DOU-PRE  TABLO  : Tableau a imprimer
C     E ...... INTEGER  NL     : Nombre de lignes du tableau
C     E ...... INTEGER  NC     : Nombre de colonnes du tableau
C 
      ENTRY OMPTDP (NOM, TABLO, NL, NC)
C 
C -----------------------------------------------------------------------
C     Test d'impression
C 
      IF (LP) THEN
C 
           WRITE (NFIMPR, 211) NOM, NL, NC
           DO II=1,NC
             IF (NC .GE. 2) THEN
               WRITE (NFIMPR, 111) II
             ENDIF
             WRITE (NFIMPR, 221) (TABLO(JJ, II), JJ=1, NL)
           ENDDO
C 
      END IF
C 
111   FORMAT ('LIGNE NUMERO : ', I5)
211   FORMAT (/, ' * TDN *, TABLEAU : ', A ,
     &        I5, 'LIGNE(S) ', I5, 'COLONNE(S)')
221   FORMAT(10D13.5)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
      SUBROUTINE OTD4IP ( NOM , TAB4 , N1 , N2 , N3 , N4  )
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
C
C     FONCTION :
C     impression transpose d'un tableau r{el double precision  plusieurs
C     indices
C
C     indique le numero de la ligne (ou de la colonne suivant le type
C     d'impression) si le nombre de lignes du tableau est sup{rieur
C     @ 6
C
C     traces normales
C
C     Exemple : CALL IMPTDN ( 'contraintes',CONTRA,1,6)
C
C     R{sultat: - TDN - , TABLEAU : contraintes  1 ligne(s) 6 colonne(s)
C
C     Exemple : CALL IMPTDN ( 'XXXXX',XXXXXX,10,6)
C
C     R{sultat:
C        - TDN - , TABLEAU : XXXXXX (transpos{) 10 ligne(s) 6 colonne(s)
C
C     ..................................................................
C
C ... E    CHARACT  NOM    : Nom du tableau @ imprimer
C ... E    DOU-PRE  TABLO  : Tableau @ imprimer
C ... E    INTEGER  NL     : Nombre de lignes du tableau
C ... E    INTEGER  NC     : Nombre de colonnes du tableau
C
      CHARACTER NOM *(*)
      INTEGER N1 , N2 , N3 , N4  , DEB1 , DEB2 , I , J
      DOUBLE PRECISION  TAB4( N1 , N2 , N3 , N4 )
C     ..................................................................
      DEB1 = 1
      DEB2 = 1
C -
      DO  I = 1 , N4
C -
        DEB1 = 1
C -
CD      CALL IMPEP(' 4eme colonne numero ', I )
C -
        DO J = 1 , N3
C -
CD        CALL IMPEP(' 3eme colonne numero ', J )
C -
CD        CALL OMPTDP(NOM , TAB4(1,1,DEB1,DEB2) , N1 , N2 )
          DEB1 = DEB1+1
C -
        END DO
C -
        DEB2 = DEB2+1
C -
      END DO
C -
      RETURN
C -
      END
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
      SUBROUTINE OTD4IN ( NOM , TAB4 , N1 , N2 , N3 , N4  )
C [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]
C
C     FONCTION :
C     impression transpose d'un tableau r{el double precision  plusieurs
C     indices
C
C     indique le numero de la ligne (ou de la colonne suivant le type
C     d'impression) si le nombre de lignes du tableau est sup{rieur
C     @ 6
C
C     traces normales
C
C     Exemple : CALL IMPTDN ( 'contraintes',CONTRA,1,6)
C
C     R{sultat: - TDN - , TABLEAU : contraintes  1 ligne(s) 6 colonne(s)
C
C     Exemple : CALL IMPTDN ( 'XXXXX',XXXXXX,10,6)
C
C     R{sultat:
C        - TDN - , TABLEAU : XXXXXX (transpos{) 10 ligne(s) 6 colonne(s)
C
C     ..................................................................
C
C ... E    CHARACT  NOM    : Nom du tableau @ imprimer
C ... E    DOU-PRE  TABLO  : Tableau @ imprimer
C ... E    INTEGER  NL     : Nombre de lignes du tableau
C ... E    INTEGER  NC     : Nombre de colonnes du tableau
C
C -
      CHARACTER NOM *(*)
      INTEGER N1 , N2 , N3 , N4  , DEB1 , DEB2 , I , J
      DOUBLE PRECISION  TAB4( N1 , N2 , N3 , N4 )
C     ..................................................................
      DEB2 = 1
C -
      DO  I = 1 , N4
C -
CD      CALL IMPEN(' 4eme colonne numero ', I )
C -
        DEB1 = 1
C -
        DO J = 1 , N3
C -
CD        CALL IMPEN(' 3eme colonne numero ', J )
C -
CD        CALL OMPTDN(NOM, TAB4(1,1,DEB1,DEB2) , N1 , N2 )
          DEB1 = DEB1+1
C -
        END DO
C -
        DEB2 = DEB2+1
C -
      END DO
C -
      RETURN
C -
      END
