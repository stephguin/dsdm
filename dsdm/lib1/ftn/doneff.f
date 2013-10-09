C     Donnees des efforts et des deplacements :
C 
C     On caracterise les zones par une suite de 6-uplets indiquant :
C 
C     1er cas
C     -------
C 
C     (numero de bord, 1er numero de colonne ou de couche, dernier numero
C     de colonne ou de couche, CARACU, CARACV, CARACW)
C 
C     CONVENTION :
C 
C     CARAC = 0  <==> effort impose nul
C     CARAC = 1  <==> blocage
C     CARAC = 2  <==> blocage noeud par noeud
C     CARAC = 3  <==> blocage noeud par noeud par developpement
C 
C     A priori, blocage des deplacements de solide autrement
C     c'est un peu bizzare
C 
C     2eme cas
C     --------
C 
C     (numero de bord,1er numero (colonne ou couche), numero de
C     noeud => 1 a 4  si nubord = 1 ou 3 !( 2,4 derivees)
C           => 1 a 6  si nuborD = 2 ou 4 !( 2,4,6 derivees)
C 
C     CARAC = 4  <==> on impose les deplacements de solide
C     CARAC = -1 <==> autre type d''effort imposes
C     CARAC = -2 <==> autre type de deplacement imposes
C     CARAC = -3 <==>  deplacement imposes noeud par noeud')
C 
      SUBROUTINE DONEFD
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'cominc_visu.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  NMAX, NLU, TEST, A1, NMAX1, TEST1, TEST2, TEST3
      INTEGER  A1SYM
      INTEGER  NMAX2, R1, R2, DEBDEG, INDDEV, INDDEG, NUBORD, D1
      INTEGER  I, J, K, L
      INTEGER  DEBUTI,  NLU1,  NLU2, LONDEV, LONDEG, A2
      INTEGER  INDZON, TYPDEP, DEBDEV, NUCO, MOD, LECINT, NFT
      INTEGER  NBZOTO
      INTEGER  AM2LC , ADM2LC, COMPT, COMBDD
C 
C     Longueur des donnees pour les deplacements de solide
C 
      INTEGER LOSOLI
C 
      CHARACTER*3 CARNU1, CARNU2
      CHARACTER*4 NUMERO
      CHARACTER*2 CCOMPT
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DONEFD')
C  
C     include 'identi.h'
C 
C -----------------------------------------------------------------------
      CHARACTER*6      CPLAST(6) , IPLAST(2), CRITER(2)
      CHARACTER*3      CENDOM(3) , ENDOMI(3)
      CHARACTER*10     ENDDIF(3)
      CHARACTER*9      BORD(4)
      CHARACTER*3      CONTRA(6), DEFORM(6)
      CHARACTER*4      DCONTR(6)
      CHARACTER*4      DDEFOR(6)
      CHARACTER*1      DEPLA(3)
C 
      BORD(1)    = 'inferieur'
      BORD(2)    = 'interieur'
      BORD(3)    = 'superieur'
      BORD(4)    = 'exterieur'
C 
      DEPLA(1)   = 'u'
      DEPLA(2)   = 'v'
      DEPLA(3)   = 'w'
C 
      CPLAST(1)  = 'epsp11' 
      CPLAST(2)  = 'epsp22' 
      CPLAST(3)  = 'epsp12' 
      CPLAST(4)  = 'epsp23' 
      CPLAST(5)  = 'epsp13' 
      CPLAST(6)  = 'epsp33' 
C 
      CENDOM(1)   = 'dfi'
      CENDOM(2)   = 'dps' 
      CENDOM(3)   = 'dpt'
C 
      IPLAST(1)  = 'SAUTP1' 
      IPLAST(2)  = 'SAUTP2' 
C 
      ENDOMI(1)   = 'di1' 
      ENDOMI(2)   = 'di2'
      ENDOMI(3)   = 'di3'
C 
      ENDDIF(1)   = 'di1-di1ini'
      ENDDIF(1)   = 'di2-di2ini'
      ENDDIF(1)   = 'di3-di3ini'
c -
      CRITER(1)  = 'CRIT-S'
      CRITER(2)  = 'CRIT-N'
C -----------------------------------------------------------------------
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------     
      CALL LECSEQ ('DONEFD', 'DONNEE DES EFFORTS ET DEPLACEMENTS ')
      NBZOTO = LECINT ('NBZOTO')
C 
      IF (SYM .AND. (.NOT. SYMPAR)) THEN
        CALL MESSAO ('PAR RAISON DE SYMETRIE, ON CREE UNE ZONE EN PLUS')
        NBZOTO = NBZOTO+1
        CALL IMPET ('NBZOTO MODIFIE', NBZOTO)
      END IF
C 
      CALL LFIDON
C 
C     OUVERTURE DES TABLEAUX DES EFFORTS
C 
C     DR  Ajout le 14/12/95 :
C     Steph : DONHER si donnees en efforts par fichiers
C             EHERTZ si donnees en effort maxi au centre sous le poincon
C             PHERTZ si donnees en pression sur la surface poinconnee ??
C    ouverture des tableaux pour ranger efforts imposes, deplacements
C    imposes (conditions aux limites)
C 
      CALL OUVTAB
C 
      IF (DONHER) THEN
C 
        IF (NBELZC .EQ. 0) THEN
C 
C       Le fichier de donnees comporte la force maximale
C 
          CALL EHERTZ
C 
        ELSE
C 
C       Le fichier de donnees comporte la repartition de pression maximale
C 
          CALL PHERTZ
C 
        ENDIF
C 
      ENDIF
C 
      NCDPIM      = 0
      NCEFIM      = 0
      NCBLIM      = 0
      CALL GESTEN ('ZONE-CARAC', 6*NBZOTO, TEST1)
      A1    = TEST1
      A1SYM = A1
      INDDEV      = 0
      INDDEG      = 0
C 
C     POUR REMPLIR LE TABLEAU 'NUM-DEVELO'
C 
      CALL DEBUEN (A2)
      DEBDEV = A2
C 
C     POUR REMPLIR LE TABLEAU 'COEFF-POLY'
C 
      CALL DEBUDP (D1)
      DEBDEG = D1
      NMAX1  = NBZOTO*NBMAT*3
      NMAX2  = 3*NMAX1
C 
C     R1 est l'adresse de depart du tableau poubelle qui doit etre
C     transforme en tableau pointeur dans par appel a pointe
C 
      CALL GSPOUE (NMAX1+NMAX2, R1)
      R2  = R1+NMAX1
C 
      IF (NBZOTO .EQ. 0) THEN
        NLU = 0
        TEST= 0
        NMAX= 0
        GOTO 2
      END IF
C  
      CALL MESSAO (
     &' On caracterise les zones par une suite de 6-uplets indiquant :
     &\ (numero de bord, 1er numero de colonne ou de couche, dernier numero
     &\ de colonne ou de couche, CARACU, CARACV, CARACW ); CONVENTION :
     &\ 1=inferieur, 2=interieur, 3=superieur, 4=exterieur
     &\ CARAC =  0 <==> effort impose nul
     &\ CARAC =  1 <==> blocage
     &\ CARAC =  4 <==> on impose les deplacements de solide
     &\ CARAC = -1 <==> autre type d''effort impose
     &\ CARAC = -2 <==> autre type de deplacement impose ')
C 
      CALL MESSAO (
     &' Pour les deplacements imposes noeud par noeud, on caracterise
     &\ le noeud et le deplacement imposes par des 6-uplets indiquant :
     &\ (numero de bord, 1er numero (colonne ou couche), numero de ddl
     &\ pour la composante => 1 a 4 (nubord = 1, 3) !( 2,4 derivees)
     &\                    => 1 a 6 (nubord = 2, 4) !( 2,4,6 derivees)
     &\ puis   (CARACU, CARACV,CARACW )
     &\ avec :
     &\ CARAC = 2  <==> blocage noeud par noeud
     &\       = 3  <==> blocage noeud par noeud par developpement
     &\       =-3  <==> deplacement imposes noeud par noeud')
C 
      NMAX     = 6*NBZOTO
C 
      IF (SYM .AND. (.NOT. SYMPAR)) THEN
        M(A1)   = 1
        M(A1+1) = 1
        M(A1+2) = NBCOL
        M(A1+3) = 0
        M(A1+4) = 0
        M(A1+5) = 1
        A1SYM   = A1+6
        NMAX  = NMAX -6
      END IF
C 
1     CALL LECLEN (M(A1SYM), NMAX, NLU)
      IF (SYM .AND. (.NOT. SYMPAR)) NLU = NLU+6
C 
      TEST = MOD(NLU,6)
      IF (TEST .NE. 0) THEN
        CALL MESSAO ('VOUS N''AVEZ PAS RENTRE DE 6-UPLETS RECOMMENCEZ')
        GOTO 1
      ENDIF
C 
2     CONTINUE
C 
      NFT              = NBFODO
      NBZONE(NFT)      = NLU/6
C 
C     POUR LE TABLEAU ZONE-CARAC
C 
      INDZON      = A1+2
C 
C     Longueur des donnees pour les deplacements de solide
C 
      LOSOLI = 0
C 
C     Fait pour un nombre de fonctions du temps = 1
C 
C     Calcul du nombre de zone a deplacements imposes
C     Bizzare il manque carac = 2
C     NCBPIM
C 
      DO I = 1, NBZONE(1)
         DO TYPDEP =1, 3
           IF (M(INDZON+TYPDEP) .EQ. -2)THEN
             DO NUCO = M(INDZON-1),M(INDZON)
               NCDPIM   = NCDPIM+1
             ENDDO
           ELSE  IF (M(INDZON+TYPDEP) .EQ. -3) THEN
             NCDPIM   = NCDPIM+1
           ELSE  IF (M(INDZON+TYPDEP) .EQ. 3) THEN
             NCBPIM   = NCBPIM+1
           ELSE  IF (M(INDZON+TYPDEP) .EQ. 4) THEN
             NCBPIM   = NCBPIM+1
           ENDIF
         ENDDO
         INDZON = INDZON+6
      ENDDO
C  
C     COMBDD : Uniquement pour le cas 3 (le blocage noeud par noeud par
C     developpement)
C 
      COMBDD = 0
      COMPT  = 0
C  
      DO I = 1, NBZONE(NFT)
C 
C       On boucle sur les zones et leurs caracteristiques
C 
        DEBUTI    = 6*(I-1)+A1
        NUBORD = M(DEBUTI)
C 
        CALL IDENTI (M(DEBUTI+1), CARNU1)
        CALL IDENTI (M(DEBUTI+2), CARNU2)
C 
        DO J = 1, 3
C 
C         On boucle sur les zones et leurs caracteristiques, en effort puis en
C         deplacement impose
C 
          IF (M(DEBUTI+2+J) .EQ. -1) THEN
C 
C           Cas des efforts imposes non nuls
C 
            CALL MESSAO (
     &      'POUR LE BORD '//bord(nubord)//' ET POUR '//depla(j)//',
     &      \POUR LA ZONE COMPRISE ENTRE LES ELEMENTS NUMERO '//CARNU1
     &      //' ET '//CARNU2//'
     &      \ LES NUMEROS DES DEVELOPPEMENTS CONCERNES SONT :')
C 
C           Stockage des numeros de developpement en series
C           correspondant a des efforts imposes non nuls
C 
            CALL LECLEN (M(DEBDEV), NMAX1, NLU1)
            M(R1+INDDEV)   = NLU1
            INDDEV         = INDDEV+1
C 
            DO K= 1, NLU1
C 
C           Stockage des caracteristiques en r ou z suivant les bords
C           associes aux  numeros de developpement en series
C           correspondant a des efforts imposes non nuls
C 
              WRITE (NUMERO(1:4), '(I4)') M(DEBDEV-1+K)
C 
              CALL MESSAO (
     &      'POUR LE BORD '//bord(nubord)//', POUR '//depla(j)//',
     &     \POUR LA ZONE COMPRISE ENTRE LES ELEMENTS NUMERO '//CARNU1
     &     //' ET '//CARNU2//' et
     &     \POUR LE DEVELOPPEMENT' //numero//' DONNER DANS L''ORDRE LES 
     &     \COEFFICIENTS ao ...an (MPa) CARACTERISANT LE POLYNOME DES 
     &     \EFFORTS SURFACIQUES (en r ou z)
     &     \(DEGRE MAX= 3 POUR r, 2 POUR z)')
C 
              CALL LECLDP (DM(DEBDEG), NMAX2, NLU2)
C 
C             Les donnees sont des efforts surfaciques (homogenes a des pressions)
C             l'unite est le MPA. Lorsque l'on developpe en series de Fourier il faut
C             multiplier par 2PI pour le terme du developpement 0 et Pi pour les autres
C 
              IF (M(DEBDEV-1+K) .EQ. 0) THEN
                DO L = 1, NLU2
                  DM(DEBDEG+L-1)= 2.D0*PI*DM(DEBDEG+L-1)
                END DO
              ELSE
                DO L = 1, NLU2
                  DM(DEBDEG+L-1)= PI*DM(DEBDEG+L-1)
                END DO
              END IF
C 
              DEBDEG         = DEBDEG+NLU2
              M(R2+INDDEG)   = NLU2
              INDDEG         = INDDEG+1
C 
            ENDDO
C 
            DEBDEV      = DEBDEV+NLU1
C 
          ELSE IF (M(DEBUTI+2+J) .EQ. -2) THEN
C 
C           Cas des deplacements imposes autres que bloques ou
C           imposes noeud par noeud
C           Interet : on impose la forme du deplacement sur tout un bord
C           exemple : raccord avec une solution plaque
C 
            CALL MESSAO (
     &      'POUR LE BORD '//bord(nubord)//' ET POUR '//depla(j)//',
     &      \POUR LA ZONE COMPRISE ENTRE LES ELEMENTS NUMERO '//CARNU1
     &      //' ET '//CARNU2//' LES
     &      \NUMEROS DE DEVELOPPEMENT CONCERNES SONT :')
C 
            CALL LECLEN (M(DEBDEV), NMAX1, NLU1)
C 
C     Stockage des numeros de developpements en series
C     correspondant a des deplacements imposes autres que bloques ou
C     imposes noeud par noeud
C 
            M (R1+INDDEV) = NLU1
C 
CD          CALL IMPEP ('NOMBRE DE DEVELOPPEMENTS LUS ', NLU1)
CD          CALL IMPTEP ('VALEURS ', M(DEBDEV), 1, NLU1)
C 
            INDDEV         = INDDEV+1
C 
C     On boucle sur les numeros de developpement
C 
            DO K= 1, NLU1
C 
              WRITE (NUMERO(1:4), '(I4)') M(DEBDEV-1+K)
C 
              IF (NUBORD .EQ. 1 .OR. NUBORD .EQ. 3) THEN
C 
C     bord(1)    = 'inferieur'
C     bord(3)    = 'superieur'
C 
C     => Polynome en r
C 
                CALL MESSAO (
     &      'POUR LE BORD '//bord(nubord)//', POUR '//depla(j)//'
     &     \POUR LA ZONE COMPRISE ENTRE LES ELEMENTS NUMERO '//CARNU1
     &      //' ET '//CARNU2//' ET,
     &     \POUR LE DEVELOPPEMENT' //numero//' DONNER DANS L''ORDRE LES
     &     \COEFFICIENTS ao ...an*(r**n) CARACTERISANT LE POLYMOME DES
     &     \DEPLACEMENTS.')
C 
              ELSE
C 
C      bord(2)    = 'interieur'
C      bord(4)    = 'exterieur'
C 
C      => Polynome en z
C 
                CALL MESSAO (
     &      'POUR LE BORD '//bord(nubord)//', POUR '//depla(j)//',
     &     \ POUR LA ZONE COMPRISE ENTRE LES ELEMENTS NUMERO '//CARNU1
     &      //' ET '//CARNU2//' et,
     &     \POUR LE DEVELOPPEMENT' //numero//' DONNER DANS L''ORDRE :
     &     \LE NOMBRE DE COEFFICIENTS DU POLYNOME EN z(EN REEL), LA
     &     \VALEUR DE CES COEFFICIENTS, LA VALEUR DES COEFFICIENTS POUR
     &     \LES POLYNOMES EN r (DEGRE MAX= 3 POUR r, 2 POUR z)')
C 
              END IF
C 
C     On stocke le degres des polynome determine par le nombre
C     constantes lues
C 
              CALL LECLDP (DM(DEBDEG), NMAX2, NLU2)
              DEBDEG         = DEBDEG+NLU2
              M(R2+INDDEG)   = NLU2
              INDDEG         = INDDEG+1
C 
            ENDDO
C 
            DEBDEV      = DEBDEV+NLU1
C 
          ELSE IF (M(DEBUTI+2+J) .EQ. -3) THEN
C 
C -----------------------------------------------------------------------
C     Cas des deplacements imposes noeud par noeud
C -----------------------------------------------------------------------
            COMPT = COMPT+1
            CALL IDENT2( COMPT , CCOMPT )
C 
            CALL MESSAO (
     &      'POUR LE BORD '//bord(nubord)//' ET POUR '//depla(j)//',
     &      \POUR L''ELEMENT NUMERO '//CARNU1
     &      //' ET LE NOEUD NUMERO '//CARNU2//' et,
     &      \ POUR LE DEPLACEMENT IMPOSE NUMERO ' //CCOMPT//' LES
     &      \NUMEROS DES DEVELOPPEMENTS CONCERNES SONT :')
C 
            CALL LECLEN (M(DEBDEV),NMAX1,NLU1)
            M(R1+INDDEV)   = NLU1
C 
CD          CALL IMPET('NOMBRE DE DEVELOPPEMENTS LUS',NLU1)
CD          CALL IMPTET('VALEURS',M(DEBDEV),1,NLU1)
C 
            INDDEV         = INDDEV+1
C 
C     On boucle sur les numeros de developpement
C 
            DO K = 1, NLU1
C 
              WRITE (NUMERO(1:4), '(I4)') M(DEBDEV-1+K)
C 
              IF (NUBORD .EQ. 1 .OR. NUBORD .EQ. 3)THEN
C 
C     bord(1)    = 'inferieur'
C     bord(3)    = 'superieur'
C     => Polynome en r
C 
                CALL MESSAO (
     &     'POUR LE BORD '//bord(nubord)//', POUR '//depla(j)//', ET
     &     \POUR L''ELEMENT NUMERO '//CARNU1
     &      //' ET LE NOEUD NUMERO '//CARNU2//' , ET
     &     \POUR LE DEVELOPPEMENT' //numero//' DONNER DANS L''ORDRE LES
     &     \COEFFICIENTS ao ...an*(r**n) CARACTERISANT LE POLYNOME DES
     &     \DEPLACEMENTS.')
C 
              ELSE
C 
C     bord(2)    = 'interieur'
C     bord(4)    = 'exterieur'
C     => Polynome en z
C 
                CALL MESSAO (
     &      'POUR LE BORD '//bord(nubord)//', POUR '//depla(j)//', ET
     &      \POUR L''ELEMENT NUMERO '//CARNU1
     &      //' ET LE NOEUD NUMERO '//CARNU2//', ET
     &     \POUR LE DEVELOPPEMENT' //numero//' DONNER DANS L''ORDRE :
     &     \LE NOMBRE DE COEFFICIENTS DU POLYNOME EN z(EN REEL), LA
     &     \VALEUR DE CES COEFFICIENTS, LA VALEUR DES COEFFICIENTS POUR
     &     \LES POLYNOMES EN r (DEGRE MAX= 3 POUR r, 2 POUR z)')
C 
              ENDIF
C 
              CALL LECLDP (DM(DEBDEG), NMAX2, NLU2)
C 
C     On stocke le degre des polynomes determine par le nombre de
C     constantes lues.
C 
              DEBDEG         = DEBDEG+NLU2
              M(R2+INDDEG)   = NLU2
              INDDEG         = INDDEG+1
C 
            ENDDO
C 
            DEBDEV = DEBDEV+NLU1
C 
          ELSE IF (M(DEBUTI+2+J) .EQ. 3) THEN
C 
C -----------------------------------------------------------------------
C     Cas des blocages imposes noeud par noeud
C     On demande les numeros de developpement
C     On n'a pas a preciser les polynomes.
C -----------------------------------------------------------------------
            COMBDD = COMBDD +1
            CALL IDENT2 (COMBDD, CCOMPT)
C 
            CALL MESSAO (
     &      'POUR LE BORD '//bord(nubord)//' ET POUR '//depla(j)//',
     &      \POUR LE BLOCAGE IMPOSE NUMERO ' //CCOMPT//'
     &      \POUR L''ELEMENT NUMERO '//CARNU1
     &      //' ET LE NOEUD NUMERO '//CARNU2//' ET LES
     &      \NUMEROS DES DEVELOPPEMENTS CONCERNES SONT :')
C 
            CALL LECLEN (M(DEBDEV),NMAX1,NLU1)
            M(R1+INDDEV)   = NLU1
C 
CD          CALL IMPET ('NOMBRE DE DEVELOPPEMENTS DE BLOCAGE LUS ', NLU1)
CD          CALL IMPTET ('VALEURS ', M(DEBDEV), 1, NLU1)
C 
            INDDEV         = INDDEV+1
            DEBDEV         = DEBDEV+NLU1
C 
          ELSE
            INDDEV      = INDDEV+1
          ENDIF
C 
        END DO
C 
      END DO
C 
      LONDEV = DEBDEV-A2
      LONDEG = DEBDEG-D1
C 
C -----------------------------------------------------------------------
C    Creation du tableau  'NUM-DEVELO'  de tous les developpements
C    pour tous les types de conditons aux limites
C    Normalement TEST2 = A2 => le tabelau est deja rempli
C -----------------------------------------------------------------------
      CALL GESTEN ('NUM-DEVELO',LONDEV,TEST2)
      CALL TESTEN (TEST2,A2 , IDPROG)
C 
C -----------------------------------------------------------------------
C    Creation du tableau 'COEFF-POLY' de tous les developpements
C    pour tous les types de conditons aux limites
C -----------------------------------------------------------------------
      CALL GESTDP ('COEFF-POLY', LONDEG, TEST3)
      CALL TESTEN ( TEST3 , D1, IDPROG)
C 
C -----------------------------------------------------------------------
C    Creation du tableau pointeur 'P-DEVELOPP'
C -----------------------------------------------------------------------
      CALL POINTE (INDDEV, M(R1), 'DEVELOPP')
C 
C -----------------------------------------------------------------------
C    Creation du tableau pointeur 'P-DEGRES  '
C -----------------------------------------------------------------------
      CALL POINTE ( INDDEG,M(R2),'DEGRES  ')
C 
      CALL IMPTET ('NOMBRE DE  DEVELOPPEMENTS',M(R1),1,INDDEV)
      CALL IMPTET ('NOMBRE DE  DEGRES',M(R2),1,INDDEG)
      CALL IMPTET ('NUMEROS DES DEVELOPPEMENTS',M(A2),1,LONDEV)
      CALL IMPTDT ('COEFFICIENTS DES POLYNOMES',DM(D1),1,LONDEG)
C  
CD    DO NFT = 1, NBFODO
CD      CALL IMPEN ( ' FONCTION DU TEMPS NUMERO ' , NFT )
CD      CALL IMPTET(' CARACTERISTIQUES DES ZONES',M(A1),6,NBZONE(NFT) )
CD    END DO
C  
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Cette routine remplit le tableau 'VAL-DEP-DV'
C     des valeurs des deplacements par developpement.
C 
C     On l'utilise a l'aide du tableau :
C     pointeur des valeurs des deplacements par developpements'
C 
      SUBROUTINE DEPIMP
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
      INTEGER  ADLOC1, ADZONE, DEBZON, ADDGDP, INDPD, ADPVDD, INDPVA
      INTEGER  I, NUCO, TYPDEP, NBDDL, NBTDDL, NBDEV, INDDEP, DBDGJ
      INTEGER  NUBORD, A2, NBDEGJ, LONPVD, LONVDD ,INCDEP
      INTEGER  ADPDRE, ADNDEV, ADPDEV, ADCOEF, PLADEP
      INTEGER  J, INDPDG, INVDPD, ADVDPD
      INTEGER  TEST1, TEST2, TEST3, NFT, NUDDL(6)
      DOUBLE PRECISION RC, A, ZC, B
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DEPIMP')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
C 
C     Adresse du tableau pointeur 'P-DEGRES  '
C 
      CALL ADTBM ('P-DEGRES  ',ADPDRE)
C 
C     Adresse du tableau 'NUM-DEVELO' de tous les developpements
C     pour tous les types de conditons limites
C 
      CALL ADTBM ('NUM-DEVELO',ADNDEV)
C 
C     Adresse du tableau pointeur 'P-DEVELOPP'
C 
      CALL ADTBM ('P-DEVELOPP',ADPDEV)
C 
C     Adresse du tableau 'COEFF-POLY' de tous les developpements
C     pour tous les types de conditons limites
C 
      CALL ADTBDM('COEFF-POLY',ADCOEF)
C 
C     Adresse du tableau de localisation des numeros
C     de noeuds ranges par ordre croissant par elements
C     ranges couche par couche
C 
      CALL ADTBM ('TLOCN1    ',ADLOC1)
C 
C     Adresse du tableau de caracteristiques des zones
C 
      CALL ADTBM ('ZONE-CARAC',ADZONE)
C 
C -----------------------------------------------------------------------
C     Pour les deplacements imposes
C -----------------------------------------------------------------------
      DEBZON  = ADZONE
      NBTDDL  = 0
C 
C     Creation du tableau des pointeurs des groupes de developpements
C     pour les ddl a deplacements imposes
C  
      CALL GESTEN('P-D -GDV-DP',NCDPIM+1,A2)
      INDPD = A2+1
C 
      CALL MENAM(A2,NCDPIM+1)
C 
C     POUR REMPLIR LE TABLEAU DDL-GDV-DP
C 
      CALL DEBUEN( ADDGDP )
C 
      NFT = NBFODO
C 
      DO I= 1, NBZONE(NFT)
C 
CD      CALL IMPET ('nuzone',i)
C 
        NUBORD = M(DEBZON)
C 
        DO TYPDEP= 1, 3
C 
C     CARAC = -3 <==>  deplacement impose noeud par noeud
C 
          IF (M(DEBZON+2+TYPDEP).EQ.-3)THEN
            NUCO =  M(DEBZON+1)
            CALL NDDLCB(TYPDEP, NUBORD, NUCO, ADLOC1,NUDDL,NBDDL)
            CALL IMPTET ('VALEUR DES DDL',NUDDL ,1,NBDDL)
            M(ADDGDP+1)  = NUDDL(  M(DEBZON+2) )
            CALL IMPET ('VALEUR DU DDL BLOQUE ', M(ADDGDP+1) )
            NBTDDL       = NBTDDL+1
            M(INDPD)     = NBTDDL
            INDPD = INDPD +1
          ELSE IF (M(DEBZON+2+TYPDEP) .EQ. -2) THEN
C 
C     CARAC = -2 <==> autre type de deplacement impose
C 
C     Cas des deplacements imposes autres que bloques ou
C     imposes noeud par noeud
C 
C     Interet : on impose la forme du deplacement sur tout un bord
C     exemple : raccord avec une solution plaque
C
            DO NUCO = M(DEBZON+1) ,M(DEBZON+2)
              CALL NDDLCB(
     &              TYPDEP,NUBORD,NUCO ,ADLOC1,M(ADDGDP+NBTDDL),NBDDL)
              NBTDDL      = NBTDDL+NBDDL
              M(INDPD)    = NBTDDL
              INDPD = INDPD +1
            END DO
C 
          ENDIF
C 
        ENDDO
C 
        DEBZON   = DEBZON+6
C 
      END DO
C 
C     Creation du tableau des developpements des ddl a deplacement imposes
C 
      CALL GESTEN('DDL-GDV-DP',NBTDDL,TEST1)
      CALL TESTEN (TEST1 , ADDGDP , IDPROG )
C 
CD    CALL IMPTEN(
CD    ' POINTEUR-DES DDL PAR GROUPE DE DEVELOPPEMENTS POUR LES
CD      DEPLACEMENTS' ,M(A2),1,NCDPIM+1)
C 
CD    CALL IMPTEN(' DDL PAR DEVELOPPEMENT POUR LES DEPLACEMENTS'
CD                ,M(ADDGDP),1,NBTDDL)
C 
C -----------------------------------------------------------------------
C     Sequence de modification des efforts pour les deplacements imposes
C -----------------------------------------------------------------------
      CALL MESSAO ('RENTREE DANS LA SEQUENCE DE MODIF DES DEPLACEMENTS')
      DEBZON   = ADZONE
C 
C     POUR LE TABLEAU ENTIER : P-V-DP-GDV
C 
      CALL DEBUEN( ADPVDD )
      INDPVA   = ADPVDD
C 
C     POUR LE TABLEAU DOUBLE : VAL-DEP-DEV
C 
      CALL DEBUDP( ADVDPD )
      INVDPD = ADVDPD
C 
      INDPDG   = ADPDRE
      INCDEP   = ADPDEV
      M(ADPVDD)= 0
      NBTDDL =0
C 
CD    CALL IMPEN( 'FONCTION DU TEMPS NUMERO' , NFT )
C 
      DO I= 1, NBZONE(NFT)
C 
        NUBORD = M(DEBZON)
C 
CD      CALL IMPEP( 'DEBZON -ADZONE +1 ', DEBZON-ADZONE +1)
C 
        DO TYPDEP= 1,3
C 
          INDDEP   = INCDEP+TYPDEP
          PLADEP  = DEBZON+2+TYPDEP
C 
          IF (M(PLADEP).EQ.-1)THEN
            NBDEV   = M(INDDEP)-M(INDDEP-1)
            INDPDG  = INDPDG+NBDEV
C 
CD          CALL IMPEP('PLACE DANS POINTEUR DE DEGRES',INDPDG)
C 
          ENDIF
C 
          IF (M(PLADEP).EQ.-2)THEN
C 
C     ON EST EN DEPLACEMENTS IMPOSES
C 
            NBDEV   = M(INDDEP)-M(INDDEP-1)
C 
            DO NUCO = M(DEBZON+1),M(DEBZON+2)
              CALL CAGEO(NUBORD,NUCO,RC,A,ZC,B)
              DO J =1, NBDEV
                 DBDGJ    = M(INDPDG+J-1)
                 NBDEGJ   = M(INDPDG+J)-DBDGJ
                 CALL CALDEP( NUBORD, NBDEGJ, DM(ADCOEF+DBDGJ),RC,A,
     &                        ZC, B , NBDDL, DM(INVDPD))
                 NBTDDL       = NBTDDL +NBDDL
                 INDPVA       = INDPVA+1
                 M(INDPVA)    = NBTDDL
                 INVDPD = INVDPD +NBDDL
              ENDDO
            ENDDO
C 
            INDPDG = INDPDG+NBDEV
C 
          ENDIF
C 
          IF (M(PLADEP).EQ.-3)THEN
C 
            CALL IMPMT 
     &       ('ON EST EN DEPLACEMENTS IMPOSES NOEUD PAR NOEUD ')
C 
            NBDEV   = M(INDDEP)-M(INDDEP-1)
C 
            CALL IMPEP ('NOMBRE DE DEVELOPPEMENTS',NBDEV)
C 
            NUCO = M(DEBZON+1)
C 
            CALL IMPET (' NUMERO DE COUCHE OU DE COLONNE', NUCO)
C 
            DO J =1, NBDEV
              DBDGJ    = M(INDPDG+J-1)
              CALL IMPET ('POUR LE DEVELOPPEMENT'
     &                   ,M(ADNDEV+M(INDDEP-1)+J-1))
              NBDEGJ   = M(INDPDG+J)-DBDGJ
              CALL IMPET (' NOMBRE DE DEGRES',NBDEGJ)
              DM( INVDPD) = DM(ADCOEF+DBDGJ)
              NBTDDL       = NBTDDL +1
              INDPVA       = INDPVA+1
              M(INDPVA)    = NBTDDL
              INVDPD = INVDPD +1
            ENDDO
C 
            INDPDG = INDPDG+NBDEV
C 
          ENDIF
C 
        END DO
C 
        INCDEP   = INCDEP+3
        DEBZON   = DEBZON+6
C 
      END DO
C 
      LONPVD = INDPVA-ADPVDD +1
      LONVDD = INVDPD -ADVDPD
C 
C     Creation du tableau des valeurs par terme de developpement
C     des deplacements imposes
C 
      CALL GESTEN('P-V-DP-GDV',LONPVD,TEST2)
      CALL TESTEN(TEST2 , ADPVDD , IDPROG )
C 
      CALL GESTDP('VAL-DEP-DV',LONVDD,TEST3)
      CALL TESTEN( TEST3 , ADVDPD ,IDPROG)
C 
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Remplit le tableau : 'NUM-DDLBLO'
C     DES DDL PAR DEVELOPPEMENT POUR LES BLOCAGES '//IDPROG
C 
C     On caracterise les zones par une suite de 6-uplets indiquant :
C     (numero de bord, 1er numero de colonne ou de couche, dernier numero
C     de colonne ou de couche, CARACU, CARACV,CARACW)
C 
C     CONVENTION :
C 
C     CARAC = 0  <==> effort impose nul
C     CARAC = 1  <==> blocage
C     CARAC = 2  <==> blocage noeud par noeud
C     CARAC = 3  <==> blocage noeud par noeud par developpement
C 
C     A priori blocage des deplacements de solide
C     autrement c'est un peu bizzare
C 
C     CARAC = -1 <==> autre type d'effort impose
C     CARAC = -2 <==> autre type de deplacement impose
C     CARAC = -3 <==> deplacement impose noeud par noeud
C 
C     (numero de bord,1er numero (colonne ou couche), numero de
C     noeuD => 1 a 4  si nubord = 1 ou 3 !( 2,4 derivees)
C         => 1 a 6  si nuborD = 2 ou 4 !( 2,4,6 derivees)
C 
      SUBROUTINE DEPBLO
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
      INTEGER  ADLOC1, ADZONE, DEBZON, ADDGDP
      INTEGER  I, NUCO, TYPDEP, NBDDL
      INTEGER  NUBORD
      INTEGER  ADPDRE, ADNDEV, ADPDEV, ADCOEF
      INTEGER  NFT, NUDDL(6), DBDGDP
C 
      CHARACTER*6 IDPROG
      PARAMETER  (IDPROG='DEPBLO')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL ADTBM ('P-DEGRES  ',ADPDRE)
      CALL ADTBM ('NUM-DEVELO',ADNDEV)
      CALL ADTBM ('P-DEVELOPP',ADPDEV)
      CALL ADTBDM('COEFF-POLY',ADCOEF)
      CALL ADTBM ('TLOCN1    ',ADLOC1)
      CALL ADTBM ('ZONE-CARAC',ADZONE)
C 
C -----------------------------------------------------------------------
C     POUR LES DEPLACEMENTS IMPOSES
C -----------------------------------------------------------------------
      DEBZON  = ADZONE
      CALL GESTEN('NUM-DDLBLO',NCBPIM,ADDGDP)
      DBDGDP = ADDGDP
C 
      NFT = NBFODO
C 
      DO I= 1, NBZONE(NFT)
C 
        NUBORD = M(DEBZON)
C 
        DO TYPDEP= 1, 3
C 
C    CARAC = 3  <==> blocage noeud par noeud par developpement
C                 => A priori blocage pour ?
C                    c'est un peu bizzare
C 
         IF (M(DEBZON+2+TYPDEP).EQ.3)THEN
           NUCO= M(DEBZON+1)
           CALL NDDLCB(TYPDEP, NUBORD, NUCO, ADLOC1, NUDDL, NBDDL)
           M(ADDGDP)  = NUDDL( M(DEBZON+2) )
           ADDGDP     = ADDGDP+1
         ENDIF
C 
        ENDDO
C 
        DEBZON   = DEBZON+6
C 
      END DO
C 
      CALL IMPTET(' DDL PAR DEVELOPPEMENT POUR LES BLOCAGES '//IDPROG
     &            ,M(DBDGDP), NCBPIM , 1)
C 
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Principe : tous les deplacements bloques (pour tous les
C     numeros de developpement) seront stockes dans ddl-bloque.
C 
      SUBROUTINE DDLBLO
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
      INTEGER     ADPRO, ADLOC1, DEBZON, DEBUT, NUBORD, NUCO, K
      INTEGER     TYPDEP, NBTDDL, NBDDL, ADZONE, I, J, DEBLOC
      INTEGER     TEST1, DDSEUL, NUDDL
      INTEGER     AM2LC, ADM2LC
      INTEGER     ADGDP2, ADGDP4, LONDFI2, LONDFI4, LONGDFI
C 
      CHARACTER*6 IDPROG
      PARAMETER  (IDPROG='DDLBLO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL ADTBM('PRODL     ',ADPRO)
      CALL ADTBM('ZONE-CARAC',ADZONE)
      CALL ADTBM('TLOCN1    ',ADLOC1)
C 
      CALL GSPOUE(12,DEBLOC)
      DEBZON   = ADZONE
C 
C -----------------------------------------------------------------------
C     POUR LE TABLEAU ENTIER : DDL-BLOQUE
C 
C     Principe : tous les deplacements bloques
C     (pour tous les numeros de developpement) seront stockes
C     dans ddl-bloque ==> on recupere l'adresse de depart du tableau
C     des developpements des ddl a deplacements imposes
C     et verification de l'adresse de depart.
C -----------------------------------------------------------------------
C 
      CALL DEBUEN (DEBUT)
C 
C     NOMBRE DE DDL BLOQUES RENTRE ACTUELLEMENT
C 
      NBTDDL = 0
C 
C     S'IL Y A EU DES DONNEES PAR FICHIER
C 
      IF (DONFIC) THEN
        CALL ADTBM ('DDL-DPI-FI2', ADGDP2)
	print*, 'ADGDP2 ', ADGDP2
        CALL LONGEN ('DDL-DPI-FI2', LONDFI2)
        CALL ADTBM ('DDL-DPI-FI4', ADGDP4)
	print*, 'ADGDP4 ', ADGDP4
        CALL LONGEN ('DDL-DPI-FI4', LONDFI4)
	LONGDFI=LONDFI2+LONDFI4
        CALL COPITE (LONDFI2, M(ADGDP2), M(DEBUT))
        CALL COPITE (LONDFI4, M(ADGDP4), M(DEBUT+LONDFI2))
        NBTDDL = LONGDFI
      END IF
       CALL IMPTET(' DDL POUR LES DEPLACEMENTS DANS '//IDPROG
     &          ,M(ADGDP2),1,LONDFI2)
       CALL IMPTET(' DDL POUR LES DEPLACEMENTS DANS '//IDPROG
     &          ,M(ADGDP4),1,LONDFI4)
C 
C     !!!!!!!!!  ON SUPPOSE QUE LES ZONES OU SONT IMPOSES LES
C     !!!!!!!!!  DEPLACEMENTS SONT LES MEMES POUR TOUTES LES
C     !!!!!!!!!  FONCTIONS DU TEMPS
C 
      DO I= 1, NBZONE(1)
C 
        NUBORD = M(DEBZON)
C 
        DO TYPDEP= 1, 3
C 
C     POUR LES BLOCAGES DE TOUTE UNE ZONE
C 
          IF (M(DEBZON+2+TYPDEP).EQ.1.OR.M(DEBZON+2+TYPDEP).EQ.-2) THEN
C 
C     CAS DES DPLACEMENTS BLOQUES OU (1)
C     DE RACCORD AVEC UNE SOLUTION PLAQUE (-2)
C 
            DO NUCO= M(DEBZON+1), M(DEBZON+2)
C 
              CALL NDDLCB (TYPDEP, NUBORD, NUCO, ADLOC1,M(DEBLOC),NBDDL)
C 
              DO J=1,NBDDL
C 
                DO K=1,NBTDDL
C 
C     EVITE DE STOCKER PLUSIEURS FOIS UN DDL DEJA BLOQUE
C 
                  IF (M(DEBLOC+J-1).EQ.M(DEBUT+K-1))  GOTO 1
                ENDDO
C 
                M(DEBUT+NBTDDL) =M(DEBLOC+J-1)
                NBTDDL    = NBTDDL+1
C 
1               CONTINUE
C 
              ENDDO
C 
            ENDDO
C 
C    NOEUDS BLOQUES
C 
          ELSE IF ( M(DEBZON+2+TYPDEP).EQ.2)THEN
C 
            NUCO = M(DEBZON+1)
            NUDDL= M(DEBZON+2)
C 
            CALL IMPET ('POUR LE BORD '//IDPROG , NUBORD)
            CALL IMPET ('POUR LE DEPLACMENT '  , TYPDEP)
            CALL IMPET ('POUR LE COCOL      '  , NUCO)
            CALL IMPET ('POUR LE DDL        '  , NUDDL)
C 
            CALL NDDLCB (TYPDEP, NUBORD, NUCO, ADLOC1,M(DEBLOC),NBDDL)
C 
            DDSEUL = M( DEBLOC+NUDDL-1)
            CALL IMPET( ' DDL SEUL'//IDPROG , DDSEUL )
C 
            DO K=1,NBTDDL
              IF (DDSEUL.EQ.M(DEBUT+K-1))  GOTO 10
            ENDDO
C 
            M(DEBUT+NBTDDL) = DDSEUL
            NBTDDL    = NBTDDL+1
C 
10          CONTINUE
C 
          ENDIF
C 
        ENDDO
C 
        DEBZON   = DEBZON+6
C 
      ENDDO
C 
C     TABLEAU DES DDL A IMPOSER POUR TOUT DEVELOPPEMENT
C 
      CALL GESTEN('DDL-BLOQUE ',NBTDDL,TEST1)
      CALL TESTEN( TEST1 ,DEBUT , IDPROG )
      CALL IMPET ('NB DE DDL BLOQUES',NBTDDL)
      CALL IMPTET ('DDL BLOQUES POUR TOUT DEVELOPPEMENTS PAR DONFIC'
     &             //IDPROG , M(DEBUT) , 1 , LONGDFI)
      CALL IMPTET ('DDL BLOQUES POUR TOUT DEVELOPPEMENTS AUTREMENT'
     &             //IDPROG , M(DEBUT+LONGDFI) , 1 , (NBTDDL-LONGDFI))     
C 
      CALL SOPOUB(AM2LC,ADM2LC)
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C    TRAITE LE CAS CARAC = 3 POUR STOCKER LES NUMEROS
C 
C      <==>  blocage noeud par noeud par developpement
C        =>  A priori blocage des deplacements de solide
C            autrement c'est un peu bizzare
C 
C    Ce programme remplit le tableau :
C 
C    ADEP-DVBLO est le tableau des adresses de depart dans M
C    correspondant aux numeros de ddl concernes par les blocages
C    ne jouant pas pour tous les developpements dans la matrice 
C    de tous les developpements des efforts

      SUBROUTINE ASDBLO
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
      INTEGER   ADZONE, ADPDEV, PLADEP
      INTEGER   DEBDEV, NBDEV, NUDEVJ
      INTEGER   I, DEBZON, J, TYPDEP
      INTEGER   N1, N4, ADNDEV, N6, NUCO, NFT
      INTEGER   DBPDEV, INDDEV
      INTEGER   NUDDL, AM2LC, ADM2LC
      INTEGER   ADDBLO, PADBLO, INDBLO, DBABLO
      INTEGER   DDLBLO, IDDBLO
C 
      CHARACTER*6 IDPROG
      DOUBLE PRECISION GROS
C 
      PARAMETER (GROS=1.D15)
      PARAMETER (IDPROG='ASDBLO')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL ADTBM ('ZONE-CARAC',ADZONE)
      CALL ADTBM ('P-DEVELOPP',ADPDEV)
      CALL ADTBM ('NUM-DDLBLO',DDLBLO)
      CALL ADTBM ('NUM-DEVELO',ADNDEV)
C 
      CALL NUTBM ('P-DEVELOPP', N1 )
      CALL IMPTET('P-DEVELOPP '//IDPROG , M(ADPDEV) ,1 ,LONGM(N1))
C 
      CALL NUTBM ('NUM-DDLBLO', N4 )
      CALL IMPTET('NUM-DDLBLO', M(DDLBLO) ,1 ,LONGM(N4))
C 
      CALL NUTBM('NUM-DEVELO', N6 )
      CALL IMPTET('NUM-DEVELO', M(ADNDEV) ,1 ,LONGM(N6))
C 
      DEBZON = ADZONE
C 
C    DEBUT DES NUMEROS DE NOEUDS BLOQUES : POUR LA ZONE 3
C 
      IDDBLO = DDLBLO
C 
C    DEBUT DU POINTEUR DES DEVELOPPEMENTS : TOUTE ZONES
C 
      DBPDEV = ADPDEV
C 
      NFT    = NBFODO
C 
      CALL POUSME( NCBLIM*NTETA , PADBLO )
      INDBLO = 0
      DBABLO = PADBLO
C 
      DO I= 1, NBZONE(NFT)
C 
        CALL IMPET( ' NUZONE ' , I )
C 
        DO TYPDEP= 1,3
C 
          INDDEV  = DBPDEV+TYPDEP
          PLADEP  = DEBZON+2+TYPDEP
C 
          IF (M(PLADEP).EQ.3)THEN
C 
C     CARAC = 3  <==> blocage noeud par noeud par developpement
C 
C     A priori blocage des deplacements de solide autrement c'est
C     un peu bizzare
C 
CD          CALL IMPMT (' BLOCAGES IMPOSES ' )
C 
            NUDDL    = M(IDDBLO)
            IDDBLO   = IDDBLO+1
C 
CD          CALL IMPET('  NUDDL BLOQUES ' , NUDDL )
C 
            DEBDEV  = M(INDDEV-1)
C 
CD          CALL IMPET('VALEUR AV DANS LE TAB DES POINTS  ', DEBDEV )
C 
            NBDEV   = M(INDDEV)-DEBDEV
C 
CD          CALL IMPET (' NOMBRE DE DEVELOPPEMENTS', NBDEV)
C 
            NUCO = M(DEBZON+1)
C 
CD          CALL IMPET (' POUR LE NUMERO DE COU OU COL ', NUCO)
CD          CALL IMPET (' POUR LE NUMERO LOCAL DE DDL  ', M(DEBZON+2) )
C 
            DO J =1, NBDEV
              CALL IMPET('POUR LE DEVELOPPEMENT DE NUMERO   ',J)
              NUDEVJ    = M(ADNDEV+DEBDEV+J-1)
C 
CD            CALL IMPET('NUMERO DU DEVELOPPEMENT',NUDEVJ)
C 
              M(DBABLO) = (NUDEVJ+NTDSFG)*NDDL + NUDDL
C 
CD            CALL IMPET('NUMERO DU DDL DANS LE TB TOTAL ', M(DBABLO) )
C 
              DBABLO    =  DBABLO+1
              INDBLO    =  INDBLO+1
            ENDDO
C 
          ENDIF
C 
        ENDDO
C 
        DBPDEV   = DBPDEV+3
        DEBZON   = DEBZON+6
C 
      ENDDO
C 
C     ADEP-DVBLO est le tableau des adresses de depart dans M
C     correspondant aux numeros de ddl concernes par les blocages
C     ne jouant pas pour tous les developpements dans la matrice de tous
C     les developpements des efforts
C 
      CALL GESTEN( 'ADEP-DVBLO', INDBLO , ADDBLO )
      CALL COPITE( INDBLO , M(PADBLO) , M(ADDBLO) )
      CALL IMPTET('ADEP-DVBLO' , M(ADDBLO) , INDBLO , 1 )
C 
      CALL SOPOUB(AM2LC,ADM2LC)
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Pour assembler les deplacements imposes non-nuls :
C 
C     CARAC = 1  <==> blocage
C     CARAC = 2  <==> blocage noeud par noeud
C     CARAC = 3  <==> blocage noeud par noeud par developpement
C                (numero de bord,1er numero (colonne ou couche), numero de
C                noeud => 1 a 4  si nubord = 1 ou 3 !( 2,4 derivees)
C                      => 1 a 6  si nubord = 2 ou 4 !( 2,4,6 derivees)
C     CARAC = -1 <==> autre type d''effort impose
C     CARAC = -2 <==> autre type de deplacement impose
C     CARAC = -3 <==> deplacement impose noeud par noeud')
C 
      SUBROUTINE ASDEPI
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
      INTEGER     ADEFFO, ADZONE , ADPDGD ,ADPDEV , ADPVDD
      INTEGER     ADVDPD, INDDEP, INCDEP,  INDVAD, PLADEP
      INTEGER     DEBDEV, NBDEV, DBDDL, NBDDL, NUDEVJ, DEBVDD
      INTEGER     NBVDD, I , DEBZON, J, TYPDEP, INDPVD, ADDGDP
      INTEGER     N1, N2, N3, N4, N5, ADNDEV, N6, NUCO  , NFT
      INTEGER     EFFORT , APRECI , EFFOLC , K
      INTEGER     AM2LC , ADM2LC
C 
      DOUBLE PRECISION   MULDIA
C 
      CHARACTER*6 IDPROG
      DOUBLE PRECISION GROS
C 
      PARAMETER (GROS=1.D15)
      PARAMETER (IDPROG='ASDEPI')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB(AM2LC,ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL ADTBM ('ZONE-CARAC',ADZONE)
      CALL ADTBDM('MAT-EFFORT',ADEFFO)
      CALL ADTBM ('P-DEVELOPP',ADPDEV)
      CALL ADTBM ('P-V-DP-GDV',ADPVDD)
      CALL ADTBDM('VAL-DEP-DV',ADVDPD)
      CALL ADTBM ('P-D -GDV-DP',ADPDGD)
      CALL ADTBM ('DDL-GDV-DP',ADDGDP)
      CALL ADTBM ('NUM-DEVELO',ADNDEV)
C 
      CALL NUTBM('P-DEVELOPP',N1)
C 
CD    CALL IMPTEP('P-DEVELOPP',M(ADPDEV),1,LONGM(N1))
C 
      CALL NUTBM('P-V-DP-GDV',N2)
C 
CD    CALL IMPTEP('P-V-DP-GDV',M(ADPVDD),1,LONGM(N2))
C 
      CALL NUTBM('P-D -GDV-DP',N3)
C 
CD    CALL IMPTEP('P-D -GDV-DP',M(ADPDGD),1,LONGM(N3))
C 
      CALL NUTBM('DDL-GDV-DP',N4)
C 
CD    CALL IMPTEP('DDL-GDV-DP',M(ADDGDP),1,LONGM(N4))
C 
      CALL NUTBDM('VAL-DEP-DV',N5)
C 
CD    CALL IMPTDP('VAL-DEP-DV',DM(ADVDPD),1,LONGDM(N5))
C 
      CALL NUTBM('NUM-DEVELO',N6)
C 
CD    CALL IMPTEP('NUM-DEVELO',M(ADNDEV),1,LONGM(N6))
C 
C     Modification pour multiplier les deplacements par les termes
C     diagonaux bloques
C 
      CALL ADTBDM('PRECISIONS',APRECI)
      CALL GSPOUD( 10 , EFFOLC)
C 
      DEBZON = ADZONE
      INDPVD = ADPDGD +1
      INCDEP = ADPDEV
      INDVAD = ADPVDD
      NFT    = NBFODO
C 
C     Le cas de blocage ne necessite pas de mofifier les seconds
C     membres mais uniquement la matrice => pas de traitement de :
C 
C     CARAC = 1  <==> blocage
C     CARAC = 2  <==> blocage noeud par noeud
C 
CD    CALL IMPEN( 'FONCTION DU TEMPS NUMERO' , NFT )
C 
C     adeffo est l'adresse de depart du tableau des efforts pour
C     la nft-ieme fonction du temps
C 
      EFFORT = ADEFFO+(NFT-1)*NDDL*NBMAT
C 
      DO I= 1, NBZONE(NFT)
C 
        DO TYPDEP= 1,3
C 
          INDDEP   = INCDEP+TYPDEP
          PLADEP  = DEBZON+2+TYPDEP
C 
          IF (M(PLADEP).EQ.-2) THEN
C 
C        Cas des deplacements imposes autres que bloques ou
C        imposes noeud par noeud
C 
C        Interet : on impose la forme du deplacement sur tout un bord
C        exemple : raccord avec une solution plaque
C 
            DEBDEV  = M(INDDEP-1)
            NBDEV   = M(INDDEP)-M(INDDEP-1)
C 
            DO NUCO = M(DEBZON+1), M(DEBZON+2)
              DBDDL   = M(INDPVD -1)
              NBDDL   = M(INDPVD)-DBDDL
              DO J =1, NBDEV
                NUDEVJ  = M(ADNDEV+DEBDEV+J-1)
                INDVAD = INDVAD +1
                DEBVDD = M(INDVAD -1)
                NBVDD = M(INDVAD)-DEBVDD
                MULDIA = DM( APRECI+NUDEVJ+NTDSFG)
                DO K = 1 ,NBDDL
                  DM( EFFOLC+K-1) = MULDIA*DM(ADVDPD +DEBVDD+K-1)
                END DO
                CALL ASVEFI( EFFORT, NUDEVJ, M(ADDGDP+DBDDL), NBDDL,
     &                       DM(EFFOLC))
              ENDDO
              INDPVD = INDPVD +1
            ENDDO
C 
          ENDIF
C 
          IF (M(PLADEP).EQ.-3)THEN
C 
C        Cas des deplacements imposes noeud par noeud
C 
            DEBDEV  = M(INDDEP-1)
            NBDEV   = M(INDDEP)-M(INDDEP-1)
C 
            NUCO = M(DEBZON+1)
            DBDDL   = M(INDPVD -1)
            NBDDL   = M(INDPVD)-DBDDL
C 
            DO J =1, NBDEV
              NUDEVJ  = M(ADNDEV+DEBDEV+J-1)
              INDVAD = INDVAD +1
              DEBVDD = M(INDVAD -1)
              NBVDD = M(INDVAD)-DEBVDD
              MULDIA = DM( APRECI+NUDEVJ+NTDSFG)
              DO K = 1 ,NBDDL
                CALL IMPDT (' VALEUR DU DEP IMPOSE '//IDPROG
     &                      , DM(ADVDPD +DEBVDD+K-1) )
                DM( EFFOLC+K-1) = MULDIA*DM(ADVDPD +DEBVDD+K-1)
                CALL IMPDT (' VALEUR DU TERME D''EFFORTS CORRESPONDANT '
     &                      //IDPROG, DM( EFFOLC+K-1)  )
              END DO
              CALL ASVEFI( EFFORT, NUDEVJ, M(ADDGDP+DBDDL), NBDDL,
     &                     DM(EFFOLC))
            ENDDO
C 
            INDPVD = INDPVD +1
C 
          ENDIF

        ENDDO
C 
        DEBZON   = DEBZON+6
        INCDEP   = INCDEP+3
C 
      ENDDO
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
C     QUE FAIT CETTE ROUTINE :
C 
C     Cette routine fait la correspondance entre deplacement
C     donne sous forme polynomiale et la valeur correspondante
C     pour les ddl
C 
C     On envoie comme arguments :
C 
C     E ...... NUBORD  numero du bord concerne si 1 ====> bord inferieur
C                                              si 2 ====> bord interieur
C                                              si 3 ====> bord superieur
C                                              si 4 ====> bord exterieur
C     E ...... NBCONS=DEGRES  le nombre de constantes definissant
C                             l'effort concerne
C     E ...... COEFF (NBCONS) la valeur de ces efforts
C     E ...... RC             rayon au centre de l'element
C     E ...... A              demi-longueur de l'element
C     E ...... ZC             hauteur au centre de l'element
C     E ...... B              demi-hauteur de l'element
C 
C     Et on recupere :
C 
C     S ...... NBDDL          nombre des ddl affectes
C     S ...... DEPCAL         valeur correspondante
C 
      SUBROUTINE CALDEP(NUBORD, NBCONS, COEFF, RC, A, ZC,
     &                  B, NBDDL, DEPCAL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      INTEGER           NUBORD , NBCONS, NBDDL
      DOUBLE PRECISION  COEFF(4), RC, A, ZC, B, DEPCAL(6)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER           I, NBCR, NBCZ
      DOUBLE PRECISION  X1, X12, X2, X22, PZ(3), PR, PDR, RLOC
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALDEP')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      IF (NBCONS.EQ.0)THEN
        GOTO 1
      ENDIF
      IF(NUBORD.NE.1.AND.NUBORD.NE.3
     &.AND.NUBORD.NE.2.AND.NUBORD.NE.4)THEN
        CALL ERREUD(0,
     $  ' MAUVAIS INDICATEUR DE BORD DANS :'//IDPROG)
      ENDIF
C 
CD    CALL IMPDP ('VALEUR DE RC  :'//IDPROG,RC)
CD    CALL IMPDP ('VALEUR DE ZC  :'//IDPROG,ZC)
CD    CALL IMPDP ('VALEUR DE  A  :'//IDPROG,A )
CD    CALL IMPDP ('VALEUR DE  B  :'//IDPROG,B )
CD    CALL IMPEP('TYPE DE BORD        ',NUBORD)
CD    DO I=1,NBCONS
CD       CALL IMPDP('VALEUR DES COEFFS ',COEFF(I))
CD    ENDDO
C 
C -----------------------------------------------------------------------
C     Debut de la sequence pour les bords inf et sup
C -----------------------------------------------------------------------
      IF(NUBORD.EQ.1.OR.NUBORD.EQ.3) THEN
        IF (NBCONS.GT.4)THEN
C 
CD        CALL IMPEP (' POUR NUBORD', NUBORD)
C 
          CALL ERREUD(0, 'DEGRES DES DEPLACEMENTS >3 DANS ' //IDPROG)
        ENDIF
C 
CD      CALL IMPEP ('VALEUR DU DEGRES+1 ',NBCONS)
C 
        X1     = RC-A
        X12    = X1*X1
        X2     = RC+A
        X22    = X2*X2
        NBDDL  = 4
        IF (NBCONS.EQ.1) THEN
          DEPCAL(1)        = COEFF(1)
          DEPCAL(3)        = COEFF(1)
        ELSE IF (NBCONS.EQ.2) THEN
          DEPCAL(1)        = COEFF(1)+COEFF(2)*X1
          DEPCAL(2)        = COEFF(2)
          DEPCAL(3)        = COEFF(1)+COEFF(2)*X2
          DEPCAL(4)        = COEFF(2)
        ELSE IF (NBCONS.EQ.3) THEN
          DEPCAL(1)        = COEFF(1)+COEFF(2)*X1+COEFF(3)*X12
          DEPCAL(2)        = COEFF(2)+2*COEFF(3)*X1
          DEPCAL(3)        = COEFF(1)+COEFF(2)*X2+COEFF(3)*X22
          DEPCAL(4)        = COEFF(2)+2*COEFF(3)*X2
        ELSE
          CALL ERREUD(0,' MAUVAIS PASSAGE DU NBCONS DANS :'//IDPROG)
        ENDIF
      ENDIF
C 
C -----------------------------------------------------------------------
C     Debut de la sequence pour les bords interieurs et exterieurs
C -----------------------------------------------------------------------
        IF(NUBORD.EQ.2.OR.NUBORD.EQ.4)THEN
          NBCZ=IDNINT(COEFF(1))
          NBCR=NBCONS-NBCZ-1
C 
CD        CALL IMPEP('VALEUR DE NBCR',NBCR)
CD        CALL IMPEP('VALEUR DE NBCZ',NBCZ)
C 
          IF (NBCZ.GT.3.OR.NBCZ.LT.0 )THEN
            CALL ERREUD(0, 'DEGRES DES DEPLACEMENTS EN z > 2 ')
          ENDIF
          IF (NBCR.GT.4.OR.NBCR.LT.0)THEN
            CALL ERREUD(0, 'DEGRES DES DEPLACEMENTS EN r > 3 ')
          ENDIF
          NBDDL =6
          IF(NUBORD.EQ.2) RLOC=RC-A
          IF(NUBORD.EQ.4) RLOC=RC+A
          PZ(1)  = COEFF(2)
          PZ(2)  = COEFF(2)
          PZ(3)  = COEFF(2)
          DO I=2,NBCZ
            PZ(1)  = COEFF(1+I)*((-B+ZC)**(I-1))+PZ(1)
            PZ(2)  = COEFF(1+I)*((ZC)**(I-1))+PZ(2)
            PZ(3)  = COEFF(1+I)*((B+ZC)**(I-1))+PZ(3)
          ENDDO
C 
C     SI NBCR=0 => LA CONSTANTE DE MULTIPLICATION DES DEPLACEMENTS EST 1.D0
C 
          IF ( NBCR .EQ. 0) THEN
            PR     = 1.D0
          ELSE
            PR     = 0.D0
          END IF
C 
          DO I = 1, NBCR
            PR     = COEFF(1+NBCZ+I)*((RLOC)**(I-1)) +PR
          ENDDO
          PDR    = 0.D0
          DO I=2,NBCR
            PDR     = DBLE(I-1)*COEFF(1+NBCZ+I)*((RLOC)**(I-2))+PDR
          ENDDO
          DEPCAL(1)  = PZ(1)*PR
          DEPCAL(2)  = PZ(1)*PDR
          DEPCAL(3)  = PZ(2)*PR
          DEPCAL(4)  = PZ(2)*PDR
          DEPCAL(5)  = PZ(3)*PR
          DEPCAL(6)  = PZ(3)*PDR
        ENDIF
C 
1     CONTINUE
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
C     Le type de deplacement  TYPDEP
C 
C            si TYPDEP=1 ==========> U
C            si TYPDEP=2 ==========> V
C            si TYPDEP=3 ==========> W
C 
C     Le numero de bord
C     Le numero de couche ou de colonne suivant le type de bord
C     L'adresse de depart du tableau TLOCN1
C 
C     On recupere :
C 
C     DDL       tableau des numeros reels des ddl (correspondant
C               aux deplacements) et aux bords stockes croissants.
C               Cette routine se sert de NDDLCB mais ne renvoie pas
C               pour les bords 2 et 4 les ddl des derivees du
C               deplacement
C     NBDDL     nombre de ddl concernes
C 
      SUBROUTINE NDDLCD( TYPDEP, NUBORD, NUCO, TLOCN1, DDL, NBDDL)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      INTEGER         NBDDL , DDL(6), TYPDEP, NUBORD, NUCO, TLOCN1
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NDDLCD')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
      CALL NDDLCB( TYPDEP, NUBORD, NUCO, TLOCN1, DDL, NBDDL)
      IF(NUBORD.EQ.2.OR.NUBORD.EQ.4)THEN
        DDL(2)      = DDL(3)
        DDL(3)      = DDL(5)
        NBDDL       = 3
      ENDIF
C 
CD    CALL RETOUD(IDPROG)
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
         SUBROUTINE NDDLCR( TYPDEP, NUBORD, NUCO, TLOCN1, DDL, NBDDL)
C 
C       ...on envoie comme arguments
C                   le type de deplacement  TYPDEP
C                   si TYPDEP=1 ==========> U
C                   si TYPDEP=2 ==========> V
C                   si TYPDEP=3 ==========> W
C                   le numero de bord
C                   le numero de couche ou de colonne suivant le type de
C                   bord
C                   l'adresse de depart du tableau TLOCN1
C          et on recupere...
C          DDL       tableau des numeros reels des ddl (correspondants
C                    au rotations)  et au bord stockes croissants.
C                    cette routine ce sert de NDDLCB mais ne renvoie pas
C                    pour les bords 2 et 4 les ddl des deplacements
C          NBDDL     nombre de ddl concernes
C======================================================================
C *     Declaration des parametres globaux
C       """"""""""""""""""""""""""""""""""
C
      INTEGER         NBDDL, DDL(6), TYPDEP, NUBORD, NUCO, TLOCN1
C
C**********************************************************************
C *   Declaration des parametres locaux
C     """""""""""""""""""""""""""""""""
C
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='NDDLCR')
C
CD    CALL WLKBCD(IDPROG)
C
C***********************************************************************
      CALL NDDLCB( TYPDEP, NUBORD, NUCO, TLOCN1, DDL, NBDDL)
      IF(NUBORD.EQ.2.OR.NUBORD.EQ.4)THEN
        DDL(1)      = DDL(2)
        DDL(2)      = DDL(4)
        DDL(3)      = DDL(6)
        NBDDL       = 3
      ENDIF
C***********************************************************************
CD    CALL RETOUD(IDPROG)
      RETURN
      END
