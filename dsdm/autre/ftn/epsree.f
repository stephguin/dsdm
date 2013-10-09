C     On envoie comme arguments :
C 
C     E ...... EPSDEV  Tableau des deformations pour un ddl range
C                      (6, NBMAT)
C     Et on recupere :
C 
C     S ...... EPREEL  Tableau des deformations correspondant aux
C                      EPREEL (6, NTETA)
C 
      SUBROUTINE EPSREE (EPSDEV, EPREEL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION    EPREEL(6*NTETA), EPSDEV(6*NBMAT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         ADCOS, ADSIN, DEBCOS(6), ADIMA, DEBSIN(6)
      INTEGER         ADREE, DEBIMA(6), DEBREE(6), K, INV
      INTEGER         AM2LC, ADM2LC
C 
C     Pour test
C 
CD    LOGICAL         LTRACP, LTRACN
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='EPSREE')
C 
CD    CALL WLKBCD(IDPROG)
C 
CD    IF  (LTRACP(1)) THEN
C 
CD       CALL OMPTDP ('EPSDEV EN ENTREE', EPSDEV(1), 6, NBMAT)
C 
CD    ENDIF
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL POUSMD (12*NTDSFG+6 + 12*NTETA, ADCOS)
C 
C     Adresse des debuts de tableaux provisoires
C 
      ADSIN   = ADCOS+6*(NTDSFG+1)
      ADIMA   = ADSIN+6*NTDSFG
      ADREE   = ADIMA+6*NTETA
C 
C 
      DEBCOS(1)  = ADCOS
      DEBSIN(1)  = ADSIN
      DEBIMA(1)  = ADIMA
      DEBREE(1)  = ADREE
C 
C     vd$l noconcur
C 
      DO K = 2 , 6
        DEBCOS(K) = DEBCOS(K-1) + NTDSFG+1
        DEBSIN(K) = DEBSIN(K-1) + NTDSFG
        DEBREE(K) = DEBREE(K-1) + NTETA
        DEBIMA(K) = DEBIMA(K-1) + NTETA
      ENDDO
C 
      CALL EPSICO (EPSDEV, DM(ADCOS), DM(ADSIN))
C 
      INV = 2
C 
C     PAUSE ' EPSREE AVANT FT01AD '
C 
C     modif vd$l concur
C     modif vd$l cncall
C 
      DO K = 1 , 6
C 
        CALL TABTDF (NTDSFG, DM(DEBCOS(K)), DM(DEBSIN(K)),
     &               DM(DEBREE(K)), DM(DEBIMA(K)))
C 
C     On recupere le tableau de longueur 2n des coefficients reels du
C     developpement en series d'exponentielles puis le tableau de longueur
C     2n des coefficients imaginaires du developpement en series
C     d'exponentielles.
C 
        CALL FT01AD (NTETA, INV, DM(DEBREE(K)), DM(DEBIMA(K)))
C 
C     CALL OMPTDN ('partie reelle, 6', DM(DEBREE),  NTETA, 1)
C 
      ENDDO
C 
C     PAUSE ' EPSREE APRES FT01AD '
C 
CD     IF ( LTRACP(1)) THEN
C 
CD        CALL OMPTDP(' TABLEAU DES EPSILON IMAGINAIRES DANS ' // IDPROG
CD    &   , DM(ADIMA) , NTETA , 6 )
C 
C     POUR TEST SUR L'OPERATION INVERSE
C 
CD     CALL POUSMD( 6*NTETA , ATEST )
CD     CALL SIGDEV( DM(ADREE) , DM(ATEST) )
CD     CALL OMPTDN( '  EPSILON DEV = ENTREE?? ' // IDPROG
CD    & , DM(ATEST),  6 , NBMAT )
C 
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C 
C     chardm ecrase dans toueps...epsree...
C     Au cours de l'execution de la procedure transp...
C     Ben pouqoui???
C 
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C 
CD     END IF
C 
      CALL TRANSP (DM(ADREE), EPREEL(1), NTETA, NEPS)
C 
CD     CALL OMPTDN ('EPSILON REELLES NTETA, 6' // IDPROG,
CD    &              EPREEL(1), 6, NTETA)
C 
      CALL SOPOUB (AM2LC,ADM2LC)
C 
CD     CALL RETOUD(IDPROG)
C 
      RETURN
      END
