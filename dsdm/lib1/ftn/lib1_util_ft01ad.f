C     Cette routine permet de ranger de maniere croissante un tableau de
C     sinus range prealablemt (LONG, NTDSFG), ou les N varient
C     de NTDSFG a 1; on les range (LONG, NTDSFG) N variant de 1 a NTDSFG
C 
C     On envoie comme arguments :
C 
C     E................ LONG     1ere dimension du tableau des sinus
C     E................ RANDEV   tableau des sinus range
C                                NTDSFG ---> -1
C                                     1 ---> LONG
C 
C     Et on recupere :
C 
C     E................ RANCRO   tableau des sinus range
C                                1 ---> NTDSFG
C     E                           1---> LONG

      SUBROUTINE SINCRO (LONG ,  RANDEV  ,  RANCRO )
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           LONG
      DOUBLE PRECISION  RANDEV( LONG*NTDSFG) , RANCRO(LONG*NTDSFG)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     I     , J  , DEBUT , FIN
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='SINCRO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACP(1))THEN
CD      CALL OMPTDP('RANDEV (LONG, NTDSFG)', RANDEV(1)
CD                  ,LONG, NTDSFG)
CD    END IF
C 
      FIN   = NTDSFG*LONG
      DEBUT = LONG+1
C 
      DO  I =1 , NTDSFG
        DO J = 1 , LONG
          RANCRO( DEBUT -J )  = RANDEV( FIN  )
          FIN                   = FIN-1
        END DO
        DEBUT  = DEBUT+LONG
      END DO
C 
CD    IF (LTRACN(1))THEN
CD      CALL OMPTDN(' RANCRO ( LONG , NTDSFG)', RANCRO(1)
CD                  , LONG , NTDSFG )
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine permet de ranger de maniere croissante un tableau de
C     sinus range prealablement (LONG, N) ou les N varient
C     de NMAX a 1; on les range (LONG, NMAX) N variant de 1 a NMAX
C 
C     On envoie comme arguments :
C 
C     E............ LONG     1ere dimension du tableau des sinus
C     E............ RANDEV   tableau des sinus range
C                            NTDSFG ---> -1
C                                 1 ---> LONG
C 
C     Et on recupere :
C 
C     E............ RANCRO   tableau des sinus range
C                            1 ---> NTDSFG
C                            1 ---> LONG

      SUBROUTINE INVRAN( LONG , NMAX , RANDEV , RANCRO )
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER           LONG , NMAX
      DOUBLE PRECISION  RANDEV( LONG*NMAX) , RANCRO(LONG*NMAX)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER     I     , J  , DEBUT , FIN
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='INVRAN')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACP(1))THEN
CD      CALL OMPTDP(' RANDEV ( LONG , NTDSFG)', RANDEV(1)
CD                  , LONG , NMAX )
CD    END IF
C 
      FIN   = NMAX*LONG
      DEBUT = LONG+1
C 
      DO  I =1 , NMAX
        DO J = 1 , LONG
          RANCRO( DEBUT -J )  = RANDEV( FIN  )
          FIN                   = FIN-1
        END DO
        DEBUT  = DEBUT+LONG
      END DO
C 
C 
CD    IF (LTRACN(1))THEN
CD      CALL OMPTDN('RANCRO (LONG, NMAX)', RANCRO(1),
CD                  LONG, NMAX)
CD    END IF
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Subroutine de developpement en series de Fourier
C     correspondant a un ensemble de ntermes valeurs egalement
C     reparties tout les pi/NTERME
C 
C     REMARQUE :  VALEUR N'EST PAS MODIFIE
C 
C     On envoie comme arguments:
C 
C     E........... NTERME  Nombre de termes ou est connue la fonction
C                          (!!!!!   MULTIPLE DE 4   !!!!!)
C     E........... VALEUR  La valeur de ces termes (aux points
C                          2k*pi/NTERME) rangement :
C                          (sin1 ,.,sinn, 1 , cos1,..., cosn)
C     Et on recupere:
C 
C     S........... VALDEV  Les composantes du developpement reel en
C     S...........         series de Fourier suivant le rangement type
C                          (tableau de longueur NTERME+1)

      SUBROUTINE REEDEV( NTERME , VALEUR ,  VALDEV )
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C
      INTEGER            NTERME
      DOUBLE PRECISION   VALEUR(NTERME) , VALDEV(NTERME+1)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         LONG , I , INV , AD1 , AD2
      INTEGER         AM2LC , ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='REEDEV')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------      
      LONG  = NTERME/2
C 
      CALL POUSMD( NTERME , AD1 )
C 
      CALL COPITD( NTERME , VALEUR(1) , DM(AD1) )
C 
      CALL GSPOUD(NTERME , AD2)
C 
      INV = 1
      CALL FT01AD( NTERME , INV , DM(AD1) , DM(AD2) )
C  
CD    CALL IMPTDT(' PARTIE REELLE     '//IDPROG, DM(AD1) , 1, LONG+1 )
CD    CALL IMPTDT(' PARTIE IMAGINAIRE '//IDPROG, DM(AD2) , 1, LONG )
C  
      DO I = 1 , LONG-1
         VALDEV(I)      = -2.D0*DM( AD2+I)
      ENDDO
C  
      VALDEV(LONG)    = 0.D0
      VALDEV(LONG+1) = DM(AD1)
C  
      DO I = 2 ,LONG
         VALDEV(I+LONG) = 2.D0*DM(AD1+I-1)
      ENDDO
C 
      VALDEV(NTERME+1) = DM(AD1+LONG)
C  
CD    CALL IMPTDN(' VALEURS DEVELOPPEES '//IDPROG
CD                , VALDEV(1) , 1 , NTERME+1)
C  
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Subroutine de derivee d'une fonction de teta
C 
C     REMARQUE :  VALDEV N'EST PAS MODIFIE
C 
C     On envoie comme arguments:
C 
C     E............ TYPE    Indique le type de rangement en entree
C                           et en sortie
C                           0 ==> (sinn ,.,sin1,1,cos1,..., cosn)
C                           1 ==> (cosn ,.,cos1,1,sin1,..., sinn)
C     E............ VALDEV  Les composantes du developpement reel en
C     E............         series de Fourier suivant le rangement type
C                           (tableau de longueur 2*NTETA+1)
C     Et on recupere :
C 
C     S............. VDEVDR  Les composantes du developpement reel de la derivee
C     S.............         series de Fourier suivant le rangement type
C                           (tableau de longueur 2*NTETA+1)
C 
      SUBROUTINE DERTET( TYPE , VALDEV , VDEVDR )
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      DOUBLE PRECISION   VALDEV(NBMAT) , VDEVDR(NBMAT+1)
      INTEGER            TYPE
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER          I ,J
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DERTET')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF ( TYPE .EQ. 0 ) THEN
C 
        DO I = 1 , NTDSFG
          VDEVDR(I) = DBLE(NTDSFG-I+1)*VALDEV(I)
        END DO
        DO I = NTDSFG+1 , NBMAT
          J = I -NTDSFG-1
          VDEVDR(I) = -DBLE(J)*VALDEV(I)
        END DO
C 
      ELSE IF ( TYPE .EQ. 1 ) THEN
C 
        DO I = 1 , NTDSFG+1
          VDEVDR(I) = -DBLE(NTDSFG-I+1)*VALDEV(I)
        END DO
        DO I = NTDSFG+2 , NBMAT
          J = I -NTDSFG-1
          VDEVDR(I) = DBLE(J)*VALDEV(I)
        END DO
C 
      END IF
C 
CD    CALL IMPTDT('VALEURS DEVELOPPEES DE LA DERIVEE DANS '//IDPROG
CD                ,VDEVDR(1),1,NBMAT )
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     REMARQUE : VALEUR N'EST PAS MODIFIE
C 
C     On envoie comme arguments:
C 
C     E............. NTERME  Nombre de termes ou est connue la fonction
C                            !!!!!   MULTIPLE DE 4   !!!!!
C     E............. VALEUR  La valeur de ces termes (aux points
C                             2k*pi/NTERME)
C     E............. TYPE    Indique le type de rangement souhaite
C                            0 ==> (sinn ,.,sin1,1,cos1,..., cosn)
C                            1 ==> (cosn ,.,cos1,1,sin1,..., sinn)
C     Et on recupere:
C 
C     S............. VALDEV  Les composantes du developpement reel en
C     S.............         series de Fourier suivant le rangement type
C                            (tableau de longueur 2*TERME+1)
C 
C     VALDEV N'A PAS A ETRE MIS A ZERO EN ENTREE
C 
      SUBROUTINE DEVSFO (NTERME, VALEUR, TYPE, VALDEV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER            NTERME, TYPE
      DOUBLE PRECISION   VALEUR(NTERME), VALDEV(NTERME+1)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         LONG, I, AD1
      INTEGER         AM2LC, ADM2LC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DEVSFO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C -----------------------------------------------------------------------
      LONG = NTERME / 2
C 
      IF (TYPE.EQ.0)THEN
C  
      CALL POUSMD( NTERME+1 , AD1 )
C  
      CALL REEDEV( NTERME , VALEUR , DM(AD1) )
C  
C     Rangement des cosinus
C  
      CALL COPITD( LONG+1 , DM(AD1+LONG) , VALDEV(LONG+1) )
C 
      DO I = 1 , LONG
        VALDEV(I) = DM(AD1+LONG-I)
      END DO
C  
C     Sequence de verification : on recalcule par varefo
C  
CD    CALL POUSMD( 3*NTERME+1 , ADIMA)
C  
C     Adresses des debuts des tableaux provisoires
C  
CD    ADREE    = ADIMA + NTERME
CD    DEBREE   = ADREE
CD    DEBIMA   = ADIMA
CD    ADSIN    = ADREE + NTERME
CD    ADCOS    = ADSIN + LONG
C  
CD    CALL VAREFO( NTERME , VALDEV , TYPE , DM(DEBREE) )
C 
CD    CALL IMPTDT(' Valeurs developpees sinus n a 1 '//IDPROG,
CD                  VALDEV(1) , 1 , LONG)
CD    CALL IMPTDT(' Valeurs developpees cosinus 0 a n '//IDPROG
CD                , VALDEV(1+LONG) , 1 , LONG+1)
C  
      ELSE IF (TYPE.EQ.1)THEN
C  
        CALL POUSMD( NTERME+1 , AD1 )
C 
C     Fournit les valeurs sin de 1 a n cosinus de 0 a n
C 
        CALL REEDEV( NTERME , VALEUR , DM(AD1) )
C  
        CALL COPITD( LONG , DM(AD1) , VALDEV(LONG+2) )
C  
        DO I = 1 ,LONG+1
          VALDEV(I) = DM( AD1 + NTERME+1-I)
        ENDDO
C  
C     Sequence de verification
C  
C     AdresseS des debuts de tableaux provisoires
C  
CD    CALL POUSMD( 2*NTERME , ADIMA)
C  
CD    ADREE    = ADIMA+NTERME
CD    DEBREE   = ADREE
CD    DEBIMA   = ADIMA
CD    ADSIN    = ADREE + NTERME
CD    ADCOS    = ADSIN + LONG
C 
CD    CALL VAREFO( NTERME , VALDEV , TYPE , VALEUR )
C  
CD    CALL IMPTDT(' Valeurs developpees sinus 1 a n '//IDPROG,
CD                  VALDEV(1+LONG) , 1 , LONG)
CD    CALL IMPTDT(' Valeurs developpees cosinus n a 0'//IDPROG
CD                , VALDEV(1) , 1 , LONG+1)
C  
      ELSE
        CALL ERREUD(0,' MAUVAIS PASSAGE DE TYPE DANS : '//IDPROG)
      ENDIF
C  
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     REMARQUE :  VALDEV N'EST PAS MODIFIE
C 
C     On envoie comme arguments:
C 
C     E........... NTERME  Nombre de termes du developpement
C                          en series de Fourier -1 = NBANGLE
C     E........... VALDEV  La valeur de ces termes ( NTERME+1)
C     E........... TYPE    Indique le type
C                           0 ==> (sinn ,.,sin1,1,cos1,..., cosn)
C                           1 ==> (cosn ,.,cos1,1,sin1,..., sinn)
C     Et on recupere:
C 
C     E........... VALEUR  La valeur de ces termes aux
C                          2k*pi/(NTERME) angles
C 
      SUBROUTINE VAREFO (NTERME, VALDEV, TYPE, VALREE)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER            NTERME, TYPE
      DOUBLE PRECISION   VALREE(NTERME), VALDEV(NTERME+1)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER         LONG, INV
      INTEGER         AM2LC, ADM2LC
      INTEGER         ADREE, ADIMA, DEBREE, DEBIMA
      INTEGER         ADCOS, ADSIN
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='VAREFO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
CD    CALL IMPET ( 'TYPE DANS ' //IDPROG , TYPE )
C 
      LONG = NTERME / 2
C 
      IF (TYPE.EQ.0)THEN
C  
C     CALL IMPTDT(' Valeurs developpees  '//IDPROG
C    &            , VALDEV(1) , 1 , NTERME+1)
C 
C ----------------------------------------------------------------------- 
C     Sequence de verification : on recalcule par tabtdf les coefficients
C     reels du developpement sous forme exponentielle
C -----------------------------------------------------------------------
      CALL POUSMD( 3*NTERME+1 , ADIMA)
C  
C     Adresse des debuts de tableaux provisoires
C  
      ADREE    = ADIMA + NTERME
      DEBREE   = ADREE
      DEBIMA   = ADIMA
      ADSIN    = ADREE + NTERME
      ADCOS    = ADSIN + LONG
C  
      CALL SINCRO( 1 , VALDEV(1) , DM(ADSIN) )
C  
      CALL COPITD( LONG+1 , VALDEV(1+LONG) , DM(ADCOS) )
C  
      CALL TABTDF( LONG , DM(ADCOS) , DM(ADSIN) ,
     &               VALREE(1) , DM(DEBIMA) )
C  
C     On recupere le tableau de longueur 2n des coefficients reels du
C     developpement en series d'exponentielles puis le tableau de longueur
C     2n des coefficients imaginaires du developpement en series
C     d'exponentielles.
C 
      INV = 2
      CALL FT01AD( NTERME , INV , VALREE(1) , DM(DEBIMA) )
C 
C     Incrementation pour aller lire les tableaux des cosinus et des sinus
C 
CD     CALL IMPTDT('VALEURS REELLES OBTENUES DANS '//IDPROG,
CD                  VALREE(1) , NTERME ,1)
C 
C      CALL IMPTDT('VALEURS IMAGINAIRES DES SAUTS '//IDPROG,
C    &              DM(DEBIMA) , NTERME ,1)
C  
C      CALL IMPTDT(' Valeurs developpees sinus n a 1 '//IDPROG,
C    &               VALDEV(1), 1, LONG)
C      CALL IMPTDT(' Valeurs developpees cosinus 0 a n '//IDPROG,
C    &               VALDEV(1+LONG), 1, LONG+1)
C  
      ELSE IF (TYPE.EQ.1)THEN
C  
C     Adresses des debuts des tableaux provisoires
C  
        CALL POUSMD( 2*NTERME , ADIMA)
C  
        ADREE    = ADIMA+NTERME
        DEBREE   = ADREE
        DEBIMA   = ADIMA
        ADSIN    = ADREE + NTERME
        ADCOS    = ADSIN + LONG
C 
        CALL INVRAN( 1 , LONG+1 , VALDEV(1)     , DM(ADCOS) )
C 
        CALL COPITD( LONG , VALDEV(2+LONG) , DM(ADSIN) )
C  
CD      CALL IMPTDT('VALEURS DEVELOPPEES SINUS 1 a N '//IDPROG,
CD                   DM(ADSIN), 1, LONG)
C  
        CALL TABTDF( LONG ,DM(ADCOS), DM(ADSIN) ,
     &               VALREE(1), DM(DEBIMA))
C  
C     On recupere le tableau de longueur 2n des coefficients reels du
C     developpement en series d'exponentielles puis le tableau de longueur
C     2n des coefficients imaginaires du developpement en series
C     d'exponentielles .
C 
        INV = 2
        CALL FT01AD( NTERME, INV, VALREE(1), DM(DEBIMA) )
C 
C     Incrementation pour aller lire les tableaux des cosinus et des sinus
C 
CD      CALL IMPTDT('VALEURS REELLES OBTENUES DANS '//IDPROG,
CD                   VALREE(1) , NTERME ,1)
C 
C       CALL IMPTDT('VALEURS IMAGINAIRES DES SAUTS '//IDPROG,
C    &               DM(DEBIMA), NTERME, 1)
C       CALL IMPTDT('VALEURS DEVELOPPEES SINUS 1 a n '//IDPROG,
C    &               VALDEV(1+LONG), 1, LONG)
C       CALL IMPTDT('VALEURS DEVELOPPEES COSINUS n a 0'//IDPROG,
C    &               VALDEV(1), 1, LONG+1)
C  
      ELSE
        CALL ERREUD(0,' MAUVAIS PASSAGE DE TYPE DANS : '//IDPROG)
      ENDIF
C  
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
