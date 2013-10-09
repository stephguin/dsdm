C     Lecture de l'effort maximal de contact de type Hertz
C     ON RECAINTEGER           NBELZC , NUDDL(6)LCULE LA REPARTITION DE PRESSION DE TYPE HERTZ
C     SUR CHAQUE ELEMENT A PARTIR DE CETTE EFFORT DE CONTACT
C     => CONCERNE W ET W'  SUR LE BORD SUP ( NUBORD = 3 )
C     LES EFFORTS SONT SUPPOSES INDEPENDANT DE TETA => NUDEV = 0
C     LA SURFACE DE CONTACT RESTE IDENTIQUE AU COURS DU TEMPS
C 
      SUBROUTINE EHERTZ
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      external H1, H2, H3, H4
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*9  NOMDIR
      CHARACTER*20 NOMFIC
      INTEGER      LNOMFI, IUNIT, NUMDIR, LGCARG
C 
      INTEGER ADEFFO, ADLOC1, ADLONG, ADRAY
      INTEGER NUCO, TYPDEP, NUDEV, NUBORD, NBDDL
C 
C     Nombre d'elements maximum de la zone de contact
C 
      INTEGER           NUDDL(6)
      DOUBLE PRECISION  TMAX, FORCONT, RAYCONT
      DOUBLE PRECISION  FACT, RM, MILG, VALEUR (2)
      DOUBLE PRECISION  INT_FHERTZ, H1, H2, H3, H4
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EHERTZ')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      NUMDIR = 12
      NOMDIR = 'donneefic'
      NOMFIC = 'effort_hertz'
      LNOMFI = LGCARG (NOMFIC, 20)
      CALL OUVFCD (NUMDIR, NOMFIC, LNOMFI, 'F', 1, IUNIT)
      READ (IUNIT, *) TMAX, FORCONT, RAYCONT
      CLOSE (IUNIT)
C 
      CALL ADTBDM ('MAT-EFFORT', ADEFFO)
      CALL ADTBDM ('LONGUEURS ', ADLONG)
      ADLONG = ADLONG-1
      CALL ADTBDM ('RAYON     ', ADRAY)
      ADRAY = ADRAY-1
      CALL ADTBM ('TLOCN1    ', ADLOC1)
C 
C     Hypotheses :
C 
      TYPDEP = 3
      NUDEV  = 0
      NUBORD = 3
C 
      FACT = -3.D0*FORCONT/(2.D0*PI*RAYCONT**2)
      VALEUR(1) = 0.D0
      VALEUR(2) = 0.D0
C 
      DO NUCO = 1, NBCOL
        CALL NDDLCB (TYPDEP, NUBORD, NUCO, ADLOC1, NUDDL, NBDDL)
C 
C       Test pour voir si on est sur le rayon de contact
C 
        IF (DM(ADRAY+NUCO)-DM(ADLONG+NUCO) .LT. RAYCONT) THEN
          NBELZC = NBELZC+1
          RM = DM(ADRAY+NUCO)
          MILG = DM(ADLONG+NUCO)
C 
C         Assemblage des efforts (nb de ddl par noeud 2)
C 
          VALEUR(1) = VALEUR(1)+FACT*INT_FHERTZ (MILG, RM, RAYCONT, H1)
          VALEUR(2) = VALEUR(2)+
     &                FACT*MILG*INT_FHERTZ (MILG, RM, RAYCONT, H2)
          CALL ASVEFI (ADEFFO, NUDEV, NUDDL, 2, VALEUR)
          VALEUR(1) = FACT*INT_FHERTZ (MILG, RM, RAYCONT, H3)
          VALEUR(2) = FACT*MILG*INT_FHERTZ (MILG, RM, RAYCONT, H4)
        ELSE
          GO TO 100
        ENDIF
      END DO
C 
C     Assemblage pour le dernier noeud
C 
100   CONTINUE
C 
      CALL ASVEFI (ADEFFO, NUDEV, NUDDL(3), 2, VALEUR)
C 
      CALL IMPET ('NBELZC', NBELZC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
