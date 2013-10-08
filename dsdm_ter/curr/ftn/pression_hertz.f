C     LECTURE DES EFFORTS DE TYPE PRESSION DE HERTZ CALCULES
C     PAR LA THEORIE DES PLAQUES  (D. ROMAND)
C     => CONCERNE W ET W' SUR LE BORD SUP ( NUBORD = 3 )
C     LES EFFORTS SONT SUPPOSES INDEPENDANTS DE TETA => NUDEV = 0
C     CE SONT DES DONNEES NOEUD PAR NOEUD DE TYPE E.F.
C     LES MAILLAGES 2D ET 3D SONT SUPPOSES IDENTIQUES SUR LA
C     SURFACE DE CONTACT
C     LA SURFACE DE CONTACT RESTE IDENTIQUE AU COURS DU TEMPS
C 
      SUBROUTINE PHERTZ
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
      CHARACTER*9  NOMDIR
      CHARACTER*20 NOMFIC
      INTEGER      LNOMFI, IUNIT, NUMDIR, LGCARG
C 
      INTEGER ADEFFO, ADLOC1
      INTEGER NUCO, TYPDEP, NUDEV, NUBORD, NBDDL
C 
C     Nombre d'elements maximum de la zone de contact
C 
      INTEGER           NUDDL(6)
      DOUBLE PRECISION  VALEUR (2)
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='PHERTZ')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      NUMDIR = 12
      NOMDIR = 'donneefic'
      NOMFIC = 'pression_hertz'
      LNOMFI = LGCARG (NOMFIC, 20)
      CALL OUVFCD (NUMDIR, NOMFIC, LNOMFI, 'F', NBELZC+1, IUNIT)
      CALL ADTBDM ('MAT-EFFORT', ADEFFO)
      CALL ADTBM ('TLOCN1    ', ADLOC1)
C 
C     Hypotheses :
C 
      TYPDEP = 3
      NUDEV  = 0
      NUBORD = 3
C 
      DO NUCO = 1, NBELZC
        CALL NDDLCB (TYPDEP, NUBORD, NUCO, ADLOC1, NUDDL, NBDDL)
        IF (NBDDL .NE. 4) THEN
          CALL ERREUD (0, 'ERREUR FATALE DANS '// IDPROG)
        END IF
C 
C     Assemblage des efforts nb de ddl concernes 2 car on assemble noeud par noeud
C 
        READ (IUNIT , *) VALEUR(1), VALEUR(2)
        CALL ASVEFI (ADEFFO, NUDEV, NUDDL, 2, VALEUR)
      END DO
C 
C     Pour le dernier noeud
C 
      READ (IUNIT , *) VALEUR(1), VALEUR(2)
      CALL ASVEFI (ADEFFO, NUDEV, NUDDL(3), 2, VALEUR)
C 
      CLOSE (IUNIT)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
