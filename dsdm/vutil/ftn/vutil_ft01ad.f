C     On envoie comme arguments :
C 
C     E............ N Le nombre de termes du developpement tel que :
C                   U = Uo +...+ U(n)cos(n)===> ici n+1 termes
C 
C     ATTENTION 2N doit etre une puissance de 2
C 
C     E............ C    Tableau de longueur n+1 des coefficients des
C                        cosinus, on recupere un tableau de longueur N
C     E............ S    Tableau de longueur n des coefficients des
C                        sinus, on recupere un tableau de longueur N
C 
C     S............ REEL Tableau de longueur 2n des coefficients reels
C                        du developpement en series d'exponentielles
C     S............IMAG  Tableau de longueur 2n des coefficients
C                        imaginaires du developpement en series
C                        d'exponentielles

      SUBROUTINE TABTDF (N, C, S, REEL, IMAG)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      INTEGER           N, N2
      DOUBLE PRECISION  C(N+1), S(N), REEL(2*N), IMAG(2*N)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER           I
C 
      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='TABTDF')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C     CALL IMPTDT('TABLEAU DES COSINUS DANS :'//IDPROG,
C    &             C(1),1,N+1)
C     CALL IMPTDT('TABLEAU DES SINUS DANS :'//IDPROG,
C    &             S(1),1,N)
      IF(MOD(N,2).NE.0) THEN
         CALL IMPET( 'VALEUR DE N DANS '//IDPROG,N)
         CALL ERREUD(0,'N N''EST PAS UNE PUISSANCE DE 2')
      ENDIF
C 
      REEL(1)        = C(1)
      IMAG(1)        = 0.D0
C 
      DO I = 2 , N
        REEL(I)       = C(I)*.5D0
        IMAG(I)       = -S(I-1)*.5D0
      ENDDO
C 
      N2 = 2*N
C 
      DO I = N+2 , N2
        REEL(I)       = C(N2-I+2)*.5D0
        IMAG(I)       = S(N2-I+1)*.5D0
      ENDDO
C 
C     TERMES INTRODUITS POUR ASSURER UNE PUISSANCE DE 2
C 
      REEL(N+1)      = C(N+1)
      IMAG(N+1)      = 0.D0
C 
C     CALL IMPTDT('TABLEAU DES COEFFS REELS DE L''ECRITURE'
C    &        // ' EXPONENTIELLE :'//IDPROG, REEL(1),2*N,1)
C     CALL IMPTDT('TABLEAU DES COEFFS IMAGINAIRES DE L''ECRITURE '
C    &          //'EXPONENTIELLE:'//IDPROG, IMAG(1),2*N,1)

C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     !!!! Attention valable pour les deplacements uniquement !!!!
C 
C     On envoie comme arguments :
C 
C     E ...... NUDDL   Tableau des ddl
C     E ...... NBDDL   Longueur du tableau des ddl
C     E ...... VADEVO  Valeurs des termes developpes
C 
C     Attention : si nbddl = ndll, le rangement en entree est du type
C     (nddl, nbmat), sinon : il est du type (nbmat, ndll)
C 
C     Et on recupere :
C 
C     S ...... VASICO  Valeurs ou les sinus et les cosinus sont regroupes
C 
         SUBROUTINE DESICO (NUDDL, NBDDL, VADEVO, VASICO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           NBDDL, NUDDL(NBDDL)
      DOUBLE PRECISION  VADEVO(NBDDL*NBMAT), VASICO(NBDDL*NBMAT)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     I, DECALI, NUDEV
      INTEGER     NG, REST
      INTEGER     AD21, AR21, AD22, AR22, K, N
      INTEGER     AM2LC, ADM2LC

      CHARACTER*6 IDPROG
      PARAMETER   (IDPROG='DESICO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C     WARNING !!!!!! SI NBDDL = NDDL LES VALEURS SONT RANGEES SOUS LA FORME
C     (NDDL,2*NTDSFG+1)
C     WARNING !!!!!! Pour utiliser FT01AD on est oblige de sortir
C     sequentiellemnet les valeurs correspondant a un meme ddl 
C     =========> en sortie les valeurs sont rangees ddl par ddl
C 
       NG = 2*NTDSFG+1
C 
CD     CALL OMPTDP ('VALEURS DEVELOPPEES ', VADEVO(1), NBDDL, NG)
C 
       IF (NBDDL .EQ. NDDL) THEN
C 
          CALL POUSME (2*NG, AD21)
          AR21  = AD21-1
          AD22  = AD21+NG
          AR22  = AD22-1
C 
C         POUR LES TERMES EN U ET W  <===> 1, 2, 5, 6
C 
          DO I = 1, NTDSFG
            M(AR21+I) = NDDL*(NTDSFG-I)
          ENDDO
          DO I = NTDSFG, NG-1
            M(AD21+I) = NDDL*I
          ENDDO
C 
CD        CALL IMPTEP ('REPLACAGE CAS 1 ', M(AD21), 1, NG)
C 
C         POUR LES TERMES EN V  <===> 3, 4
C 
          DO I = 1, NTDSFG
            M(AR22+I) = NDDL*(NTDSFG+I)
          ENDDO
          DO I = 0, NTDSFG
            M(AD22+NTDSFG+I) = NDDL*(NTDSFG-I)
          ENDDO
C 
CD        CALL IMPTEP ('REPLACAGE CAS 2 ', M(AD22), 1, NG)
C 
          K       = 1
          DECALI  = 0
          DO I = 1, NNOEUD
C 
CD          CALL IMPEP ('Pour le dll ', DECALI+1)
C 
            DO N  = 1, NG
              VASICO(K)      = VADEVO(DECALI+M(AR21+N)+1)
              K              = K+1
            ENDDO
C 
CD          CALL IMPTDP ('VALEURS SIN-COS ', VASICO(K-NG), 1, NG)
CD          CALL IMPEP ('Pour le dll ', DECALI+2)
C 
            DO N  = 1, NG
              VASICO(K)      = VADEVO(DECALI+M(AR21+N)+2)
              K              = K+1
            ENDDO
C 
CD          CALL IMPTDP ('VALEURS SIN-COS ', VASICO(K-NG),1,NG)
CD          CALL IMPEP ('Pour le dll ', DECALI+3)
C 
            DO N  = 1, NG
              VASICO(K)      = VADEVO(DECALI+M(AR22+N)+3)
              K              = K+1
            ENDDO
C 
CD          CALL IMPTDP ('VALEURS SIN-COS ', VASICO(K-NG), 1, NG)
CD          CALL IMPEP ('Pour le dll ', DECALI+4)
            DO N  = 1, NG
              VASICO(K)      = VADEVO(DECALI+M(AR22+N)+4)
              K              = K+1
            ENDDO
C 
CD          CALL IMPTDP ('VALEURS SIN-COS ', VASICO(K-NG), 1, NG)
CD          CALL IMPEP ('Pour le dll ', DECALI+5)
C 
            DO N  = 1, NG
              VASICO(K)      = VADEVO(DECALI+M(AR21+N)+5)
              K              = K+1
            ENDDO
C 
CD          CALL IMPTDP ('VALEURS SIN-COS ', VASICO(K-NG), 1, NG)
CD          CALL IMPEP ('Pour le dll ', DECALI+6)
C 
            DO N  = 1, NG
              VASICO(K)      = VADEVO(DECALI+M(AR21+N)+6)
              K              = K+1
            ENDDO
C 
CD          CALL IMPTDP ('VALEURS SIN-COS ', VASICO(K-NG), 1, NG)
C 
            DECALI  = DECALI + 6
          ENDDO
        ELSE
C 
C       WARNING !!!!!! SI NBDDL DIFFERENT DE NDDL LES VALEURS SONT RANGEES
C       SOUS LA FORME (2*NTDSFG+1,NBDDL)
C 
          DO I=1, NBDDL
            REST = MOD(NUDDL(I),6)
            IF (REST.EQ.3.OR.REST.EQ.4)THEN
              DECALI = NG*I-NTDSFG
C 
C       WARNING !!!!!! LES SINUS SONT AUSSI RANGES PAR ORDRE CROISSANT
C 
              VASICO(DECALI)                = VADEVO(DECALI)
C 
C       vd$l nodepchk
C 
              DO NUDEV =1, NTDSFG
                VASICO(DECALI-NTDSFG-1+NUDEV)  = VADEVO(DECALI+NUDEV)
                VASICO(DECALI+NUDEV)           = VADEVO(DECALI-NUDEV)
              ENDDO
            ELSE
              DECALI = NG*I-NTDSFG
              VASICO(DECALI) = VADEVO(DECALI)
C 
C       vd$l nodepchk
C 
             DO NUDEV =1, NTDSFG
                VASICO(DECALI-NTDSFG-1+NUDEV)  = VADEVO(DECALI-NUDEV)
                VASICO(DECALI+NUDEV)           = VADEVO(DECALI+NUDEV)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C 
CD     CALL IMPTDN(' VALEURS DEVELOPPEES SIN-COS',VASICO(1),NG,NBDDL)
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
C     E ...... EPSDEV  Tableau des deformations pour un ddl range
C                      (6, NBMAT)
C 
C     Et on recupere :
C 
C     E ...... EPSCOS  Tableau des deformations correspondant aux
C                           COSINUS ( NTDSFG+1 ,6 )
C     E ...... EPSSIN  Tableau des deformations correspondant aux
C                           SINUS ( NTDSFG ,6 )
C 
         SUBROUTINE EPSICO (EPSDEV, EPSCOS, EPSSIN)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION     EPSDEV(6*NBMAT) , EPSSIN(6*(NTDSFG))
      DOUBLE PRECISION      EPSCOS(6*(NTDSFG+1))
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  AM2LC , ADM2LC , ADSIN
C 
      INTEGER ADE12C , ADE12S , ADE23C , ADE23S
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='EPSICO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC , ADM2LC)
C 
C -----------------------------------------------------------------------
C     REMPLISSAGE DANS L'ORDRE SIN(1),...SIN(N) <=> EPSSIN ,KS
C  
      CALL POUSMD( (NEPS+4)*NTDSFG , ADSIN )
C 
      ADE12C = ADSIN +NEPS*NTDSFG
      ADE12S = ADE12C+NTDSFG
      ADE23C = ADE12S+NTDSFG
      ADE23S = ADE23C+NTDSFG
C 
      CALL SINCRO ( NEPS ,  EPSDEV(1)  , DM(ADSIN) )
      CALL TRANSP ( DM(ADSIN), EPSSIN(1) ,NEPS ,NTDSFG )
C 
      CALL POUSMD( NEPS*NTDSFG , ADSIN )
      CALL SINCRO ( NEPS ,  EPSDEV(1)  , DM(ADSIN) )
      CALL TRANSP ( DM(ADSIN), EPSSIN(1) ,NEPS ,NTDSFG )
C 
C     REMPLISSAGE DANS L'ORDRE SIN(1),...SIN(N) <=> EPSSIN ,KS
C 
      CALL TRANSP ( EPSDEV(6*NTDSFG+1) ,EPSCOS(1) ,NEPS ,NTDSFG+1 )
C 
      CALL COPID( EPSSIN( 2*NTDSFG+1)     , DM(ADE12C) , NTDSFG , 1)
      CALL COPID( EPSSIN( 3*NTDSFG+1)     , DM(ADE23C) , NTDSFG , 1)
      CALL COPID( EPSCOS( 2*(NTDSFG+1)+2) , DM(ADE12S) , NTDSFG , 1)
      CALL COPID( EPSCOS( 3*(NTDSFG+1)+2) , DM(ADE23S) , NTDSFG , 1)
C 
      CALL COPID( DM(ADE12S) , EPSSIN( 2*NTDSFG+1  )   , NTDSFG , 1 )
      CALL COPID( DM(ADE23S) , EPSSIN( 3*NTDSFG+1  )   , NTDSFG , 1 )
      CALL COPID( DM(ADE12C) , EPSCOS( 2*(NTDSFG+1)+2) , NTDSFG , 1 )
      CALL COPID( DM(ADE23C) , EPSCOS( 3*(NTDSFG+1)+2) , NTDSFG , 1 )
C 
CD    CALL OMPTDN('TABLEAU DES EPSCOS', EPSCOS(1) , NTDSFG+1 , 6 )
CD    CALL OMPTDN('TABLEAU DES EPSSIN', EPSSIN(1) , NTDSFG , 6 )
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
