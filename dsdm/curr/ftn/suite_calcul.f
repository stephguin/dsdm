C     Cette routine determine si on peut passer a la suite du calcul
C     en temps; dans ce cas on remplit les tableaux d'initialisation
C     pour les couches et les interfaces pour la discretisation suivante.
C 
      SUBROUTINE SUICAL
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
      INTEGER AERREU
      DOUBLE PRECISION PRECIS, FERLOC
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
C 
      PARAMETER (PRECIS=1.D -30)
      PARAMETER (IDPROG='SUICAL')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBDM ('ERREUR-LOC', AERREU)
C 
      FERLOC = DM(AERREU+NPICET)
C 
      IF  (FERLOC .LT. PRECIS) THEN
C 
        CALL IMPET  ('ON A CONVERGE A L''ETAPE GLOBALE : ', NBETLC-1)
        CALL IMPDT  ('ERREUR GLOBALE FINALE            : ', FERLOC)
        CALL IMPDT  ('POUR UNE PRECISION DE            : ', PRECIS)
        CALL MESSAO ('STOP FIN NORMALE DANS            : '//IDPROG)
        STOP
C 
      ENDIF
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
