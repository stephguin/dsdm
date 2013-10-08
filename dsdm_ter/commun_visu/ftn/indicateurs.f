
C     Cette routine calcule l'integrale d'une fonction
C 
C     On envoie comme arguments :
C 
C     E ...... TYPE    0 = Integrale d'une fonction non nulle
C                          a l'origine normale (cas rare)
C     E ...... TYPE    1 = Integrale d'une fonction non nulle
C                          a l'origine de type accroissement
C     E ...... TYPE    2 = Integrale d'une fonction non nulle
C                          a l'origine de type produit de 2  accroissements
C     E ...... TYPE    3 = Integrale d'une fonction nulle
C                          a l'origine de type acroissement
C     E ...... TYPE    4 = Integrale d'une fonction nulle a l'origine normale
C                          mais vu le stockage c'est le cas le moins courant
C 
C     E ...... ADINTE  adresse du tableau des intervalles en temps
C     ES...... FOTEM1  valeur de la 1ere fonction du temps pour les npicet
C 
C     Et on recupere :
C 
C     ES...... FOTEM2  valeur de la fonction du temps pour integrer
C 
      SUBROUTINE intgft (TYPE, ADINTE, FOTEM1, FOTEM2)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER ADINTE, TYPE
C 
      DOUBLE PRECISION FOTEM1(NPICET), FOTEM2(NPICET)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  TEMPS, DBVLT
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      DOUBLE PRECISION VALPR1, MULT, FI
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='intgft')
C 
CD    CALL WLKBCD(IDPROG)
C 
C -----------------------------------------------------------------------
CD    IF (LTRACP(1)) THEN
CD      CALL OMPTDP('1ERE FONCTION DU TEMPS ',
CD                   FOTEM1(1), NPICET, 1)
CD    END IF
C 
      IF (TYPE .EQ. 0) THEN
C 
C     On integre une fonction nulle a l'origine
C     sans accroissement (derivee pure)
C 
        MULT      = DM(ADINTE)
        DBVLT     = ADINTE-1
        FI        = FOTEM1(1)
        VALPR1    = FOTEM1(1)
        FOTEM2(1) = MULT*(VALPR1+FI)*.5D0
C 
        DO TEMPS = 2, NPICET
C 
          MULT   = DM(DBVLT+TEMPS)
          FI     = FOTEM1(TEMPS)
C 
          FOTEM2(TEMPS) = FOTEM2(TEMPS-1) +MULT*(FI+ VALPR1)*.5D0
C 
          VALPR1 = FI
C 
        END DO
C 
C     TYPE 1 = Integrale d'une fonction non nulle a l'origine de type accroissement
C 
      ELSE IF (TYPE  .EQ. 1 )THEN
C 
        VALPR1    = FOTEM1(1)
        FI        = FOTEM1(1)
        FOTEM2(1) = (FI+ VALPR1)*.5D0
C 
        DO TEMPS = 2 , NPICET
C 
          FI     = FOTEM1(TEMPS)
C 
          FOTEM2(TEMPS) = FOTEM2(TEMPS-1) +(FI+ VALPR1)*.5D0
C 
          VALPR1 = FI
C 
        END DO
C 
       ELSE IF (TYPE .EQ. 2) THEN
C 
C      TYPE 2 = Integrale d'une fonction non nulle a l'origine de type produit
C      de 2 accroissements
C 
        FI        = FOTEM1(1)
        VALPR1    = FOTEM1(1)
        MULT      = 1.D0/DM(ADINTE)
        FOTEM2(1) = MULT*(FI+ VALPR1)*.5D0
C 
        DBVLT     = ADINTE-1
C 
        DO TEMPS = 2, NPICET
C 
          MULT   = 1.D0/DM(DBVLT+TEMPS)
C 
          FI     = FOTEM1(TEMPS)
C 
          FOTEM2(TEMPS) = FOTEM2(TEMPS-1) +MULT*(FI+ VALPR1)*.5D0
C 
          VALPR1 = FI
C 
        END DO
C 
       ELSE IF (TYPE  .EQ. 3 )THEN
C 
C      TYPE 3 = Integrale d'une fonction nulle a l'origine de type accroissement
C 
          VALPR1    = FOTEM1(1)
          FOTEM2(1) =  VALPR1*.5D0
C 
          DO TEMPS = 2, NPICET
C 
            FI     = FOTEM1(TEMPS)
C 
            FOTEM2(TEMPS) = FOTEM2(TEMPS-1) +(FI+ VALPR1)*.5D0
C 
            VALPR1 = FI
C 
          END DO
C 
       ELSE IF (TYPE .EQ. 4) THEN
C 
C      TYPE 4 = Integrale d'une fonction normale et  nulle a l'origine, vu le
C      stockage  des taux. C'est le cas le moins courant.
C 
        VALPR1    = FOTEM1(1)
        MULT      = DM(ADINTE)
        FOTEM2(1) = MULT*VALPR1*.5D0
C 
        DBVLT     = ADINTE-1
C 
        DO TEMPS = 2, NPICET
C 
          MULT   = DM(DBVLT+TEMPS)
C 
          FI     = FOTEM1(TEMPS)
C 
          FOTEM2(TEMPS) = FOTEM2(TEMPS-1) +MULT*(FI+ VALPR1)*.5D0
C 
          VALPR1 = FI
C 
        END DO
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
C     Cette routine effectue le dessin du travail des efforts
C     exterieurs en fonction du parametre de charge.
C 
C     On envoie comme arguments :
C 
C     E ...... TYPE    SI TYPE = 0 alors on ne visualise pas
C                      SI TYPE = 2 ou 3 on visualise
C     E ...... NETGLO  numero d'etape globale
C 
      SUBROUTINE DTEXPR (TYPE, NETGLO)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'cominc_visu.h'
C 
      INTEGER TYPE, NETGLO
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER  AERREU, I
C 
      INTEGER EPSCH8, SGNC12, SIGC11, SAUCH9
      INTEGER DEPCH0, FTECH6, FTSCH5
C 
      INTEGER LONDEP, LONDER, LONEPS, LONSAU
      INTEGER DEVALT, VFTTOT, ADINTE
C 
      INTEGER EPSSOL, SAUSOL, SIGSOL, SGNSOL
      INTEGER TRAVEX, FOTEMP, LONRES, PARAMC
C 
      INTEGER PUISEX
      INTEGER DEPUIS
      INTEGER INDPUI
      INTEGER FCTINS
      INTEGER DSIGSO, DSGNSO
C 
      CHARACTER*50     CARABS, CARORD
C 
      DOUBLE PRECISION TRALOC
C 
C     DCONV1 est un double qui indique s'il y a deja eu convergence
C     au moins une fois sur tout un intervalle de temps. <=> 1.d0
C     On l'initie a .FALSE. et on le SAUVEGARDE.
C 
CD    LOGICAL  LTRACN , LTRACP
C 
      CHARACTER*3 CARETG
      CHARACTER*6 IDPROG
      INTEGER     AM2LC      , ADM2LC
      PARAMETER (IDPROG='DTEXPR')
C 
CD    CALL WLKBCD(IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL IDENTI (NETGLO, CARETG)
      CALL INFODP ('FT-EPS-'//CARETG, FTECH6, FTEPMX)
C 
      FTEPMX = FTEPMX/NPICET
C 
      CALL INFODP ('FT-CON-'//CARETG, FTSCH5, FTSIMX)
C 
      FTSIMX = FTSIMX/NPICET
C 
C     Creation pour l'etape globale consideree des tableaux
C     Travail des efforts exterieurs
C     Puissance des efforts exterieurs
C     Puissance des taux des efforts exterieurs
C     Integration de la puissance des taux des efforts exterieurs
C 
      CALL GESTDP ('TRAVEX-'//CARETG, NPICET+1, TRAVEX)
      DM(TRAVEX) = 0.D0
      CALL GESTDP ('PUISEX-'//CARETG, NPICET, PUISEX)
      CALL GESTDP ('DEPUIS-'//CARETG, NPICET, DEPUIS)
      CALL GESTDP ('INDPUI-'//CARETG, NPICET+1, INDPUI)
      DM(INDPUI) = 0.D0
      CALL GESTDP ('FCTINS-'//CARETG, NPICET+1, FCTINS)
      DM(FCTINS) = 0.D0
C 
      CALL ADTBDM ('PARAM-CHAR', PARAMC)
      CALL ADTBDM ('ERREUR-LOC', AERREU)
      CALL ADTBDM ('INTE-TEMPS', ADINTE)
      CALL ADTBDM ('VALE-TEMPS', DEVALT)
      CALL ADTBDM ('F-TEMPS-DO', FOTEMP)
C 
C     Tableaux des deformations et contraintes admissibles TOTALES
C  
      CALL ADTBDM ('EPS-AD-TOT', EPSCH8)
      CALL ADTBDM ('SIG-AD-TOT', SIGC11)
      CALL ADTBDM ('DEPLA-ADMI', DEPCH0)
C 
C     POUR LES INTERFACES, tableaux des sauts et contraintes admissibles totales
C 
      IF (NBINT .GT. 0) THEN
C 
        CALL ADTBDM ('SAU-AD-TOT', SAUCH9)
        CALL ADTBDM ('SGN-AD-TOT', SGNC12)
C 
      ELSE
C 
        SAUCH9 = EPSCH8
        SGNC12 = SIGC11
C 
      ENDIF
C 
      LONDEP = NDDL*NBMAT
      LONDER = NDDL*NTETA
      LONEPS = NEPS*NGAU1*NTETA
      LONSAU = NSAU*NGAU2*NTETA
C 
C     On calcule le travail des efforts exterieurs:
C     = travail des (efforts donnes + efforts de liaisons)
C 
      LONRES = 3*(LONEPS+LONSAU) + NPICET + 1
      CALL POUSMD (LONRES, EPSSOL)
C 
      SAUSOL = EPSSOL+LONEPS
      DSIGSO = SAUSOL+LONSAU
      DSGNSO = DSIGSO+LONEPS
      SIGSOL = DSGNSO+LONSAU
      SGNSOL = SIGSOL+LONEPS
C 
C     VFTTOT est la valeur du parametre de charge
C 
      VFTTOT = SGNSOL+LONSAU
C 
C     On calcule au moyen de vcupte les valeur reelles au piquet
C     de temps I des :
C       - deformations
C       - sauts
C       - contraintes
C       - contraintes normales
C     puis on calcule l'oppose du travail des efforts interieurs,
C     que l'on divise par le parametre de charge.
C 
      DO I = 1, NPICET
C 
         CALL VDERTE (LONEPS, FTEPMX, FTEPMX, DM(EPSCH8),
     &                DM(FTECH6), I, DM(EPSSOL))
         CALL VCUPTE (LONEPS, FTSIMX, FTSIMX, DM(SIGC11),
     &                DM(FTSCH5), I, DM(SIGSOL))
         CALL VDERTE (LONEPS, FTSIMX, FTSIMX, DM(SIGC11),
     &                DM(FTSCH5), I, DM(DSIGSO))
C 
          IF (NBINT .GT. 0) THEN
C 
            CALL VDERTE (LONSAU, FTEPMX, FTEPMX, DM(SAUCH9),
     &                   DM(FTECH6), I, DM(SAUSOL))
            CALL VCUPTE (LONSAU, FTSIMX, FTSIMX, DM(SGNC12),
     &                   DM(FTSCH5), I, DM(SGNSOL))
            CALL VDERTE (LONSAU, FTSIMX, FTSIMX, DM(SGNC12),
     &                   DM(FTSCH5), I, DM(DSGNSO))
C 
          END IF
C 
          CALL TRAGLO (DM(EPSSOL), DM(SAUSOL), DM(SIGSOL),
     &                 DM(SGNSOL), DM(PUISEX+I-1), TRALOC)
          CALL TRAGLO (DM(EPSSOL), DM(SAUSOL), DM(DSIGSO),
     &                 DM(DSGNSO), DM(DEPUIS+I-1), TRALOC)
C 
        END DO
C 
        CALL intgft (2, ADINTE, DM(DEPUIS), DM(INDPUI+1))
        CALL intgft (3, ADINTE, DM(PUISEX), DM(TRAVEX+1))
C 
        CALL IMPTDT 
     &     ('VALEUR DE LA PUISSANCE DES EFFORTS EXTERIEURS          ',
     &       DM(PUISEX), NPICET, 1)
        CALL IMPTDT 
     &     ('VALEUR DU TRAVAIL DES EFFORTS EXTERIEURS               ',
     &       DM(TRAVEX), NPICET+1, 1)
        CALL IMPTDT  
     &     ('VALEUR DE LA PUISSANCE DES TAUX DES EFFORTS EXTERIEURS ',
     &       DM(DEPUIS), NPICET, 1)
        CALL IMPTDT 
     &     ('VALEUR DE L''INTEGRALE DE LA PUISSANCE DU TAUX DES'
     &    //'EFFORTS EXTERIEURS                                     ',                                    
     &       DM(INDPUI), NPICET+1, 1)
C 
C     Calcul de la fonction indicatrice de la stabilite,
C     c'est a dire de INT (spoint* epoint/ char_point**2)
C     qui doit etre une droite quel que soit le chargement
C     impose par une fonction du temps
C 
      CALL CALFIN (ADINTE, DEPUIS, DM(FCTINS))
C 
      CALL IMPTDT 
     &     ('VALEUR DE LA FONCTION D''INSTABILITE                   ',
     &       DM(FCTINS), NPICET+1, 1)
C 
C     Dessin du travail des efforts exterieurs en fonction du parametre de charge
C 
      IF ((TYPE .EQ. 2) .OR. (TYPE .EQ. 3)) THEN
C 
        CARABS = 'temps'
        CARORD = 'Puissance'
C 
C       DESFON est compilee a vide dans delami, mais appelle des routines
C       graphiques dans visu_dsdm (---> repertoire faux_commun_visu)
C 
        CALL DESFON (NPICET, DM(DEVALT+1), CARABS,
     &               DM(PUISEX), CARORD)
C 
        CARORD = 'Travail des efforts exterieurs'
C 
        CALL DESFON (NPICET+1, DM(DEVALT), CARABS,
     &               DM(TRAVEX), CARORD)
C 
        CARORD = 'Puissance taux efforts exterieurs'
C 
        CALL DESFON (NPICET, DM(DEVALT+1), CARABS,
     &               DM(DEPUIS), CARORD)
C 
        CARORD = 'Integrale puissance taux efforts exterieurs'
C 
        CALL DESFON(  NPICET+1 , DM(DEVALT) , CARABS
     &             , DM(INDPUI) , CARORD  )
C 
        CARORD = 'Fonction d''instabilite'
C 
        CALL DESFON (NPICET+1, DM(DEVALT), CARABS,
     &               DM(FCTINS), CARORD)
C 
      END IF
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
C     Cette routine calcule la fonction de detection de l'instabilite
C     dans le cas d'un effort impose par une fonction du temps.
C 
C     On envoie comme arguments :
C 
C     E ...... ADINTE  adresse du tableau des intervalles en temps
C     E ...... DEPUIS  adresse du tableau de la puissance du taux des efforts exterieurs
C 
C     Et on recupere :
C 
C     ES...... FOTEMP  valeur de la fonction d'instabilite
C 
      SUBROUTINE CALFIN (ADINTE, DEPUIS, FOTEMP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
C 
      INTEGER ADINTE, DEPUIS
C 
      DOUBLE PRECISION FOTEMP(NPICET+1)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER PTAUDC, ITACH2, TYPE, I, J, AFOTEM
C 
      LOGICAL TROUVE
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='CALFIN')
C 
CD    CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
      CALL POUSMD (NPICET, PTAUDC)
      CALL POUSMD (NPICET, ITACH2)
      CALL ADTBDM ('F-TEMPS-DO', AFOTEM)
C 
      DO I = 0, NPICET-1
C 
        DM(ITACH2+I) =  DM (AFOTEM+I)*DM(AFOTEM+I)
C 
      END DO
C 
CD    CALL IMPTDT
CD     ('Tableau TAUX-CHA-2 '//IDPROG, DM(ITACH2+1), NPICET, 1)
C 
C     Determination du type d'integration. TYPE = 0 : Integrale d'une fonction non nulle
C     a l'origine normale (cas rare).
C 
      TYPE = 0
C 
      DO I = 1, NPICET
C 
        IF (DABS(DM(ITACH2+I-1)) .LT. 1.D -14) THEN
C 
          CALL IMPET ('POUR LE PAS DE TEMPS ', I)
C 
          CALL IMPDT
     &     ('CAS SPECIAL OU LE TAUX DE CHARGEMENT AU CARRE VAUT ',
     &     DM(ITACH2+I-1))
C 
          IF( DABS(DM(DEPUIS+I-1)). LT. 1.D -10 )  THEN
C 
            CALL IMPDT 
     &       ('CAS ENCORE PLUS SPECIAL LA PUISSANCE DU TAUX VAUT ',
     &       DM(DEPUIS+I-1))
C 
            TROUVE = .FALSE.
C 
C           On recherche le pas de temps significatif le plus proche de I :
C                par valeur inferieure si I > 1 ;
C                par valeur superieure si J > 1 ;
C 
            IF (I . GT . 1) THEN
C 
              J = I
C 
              TROUVE = .FALSE.
C 
              DO WHILE (.NOT. TROUVE .AND. J .NE. 0)
C 
                J = J-1
                IF (DABS(DM(DEPUIS+J-1)). GT. 1.D -5 ) TROUVE= .TRUE.
C 
              END DO
C 
            ELSE
C 
C           TYPE = 4 : Integrale d'une fonction nulle a l'origine normale
C           mais vu le stockage c'est le cas le moins courant
C 
              TYPE = 4
C 
            END IF
C 
            IF (I .EQ. 1 .OR. (.NOT. TROUVE)) THEN
C 
              TROUVE = .FALSE.
              J = I
C 
              DO WHILE (.NOT. TROUVE .AND. J .NE. NPICET+1)
C 
                J = J+1
                IF (DABS(DM(DEPUIS+J-1)) .GT. 1.D-5) TROUVE= .TRUE.
C 
              END DO
C 
            END IF
C 
              IF (TROUVE) THEN
C 
                IF (DABS(DM(ITACH2+J-1)) .LT. 1.D-14) THEN
C 
                 CALL IMPET ('POUR LE PAS DE TEMPS ', I)
                 CALL IMPET 
     &                     ('POUR LE 1er PAS DE TEMPS SIGNIFICATIF ', J)
                 CALL IMPDT ('LA VALEUR DE LA PUISSANCE DU TAUX VAUT',
     &                        DABS(DM(DEPUIS+J-1)))
                 CALL IMPDT 
     &                     ('MAIS LE TAUX DE CHARGEMENT AU CARRE VAUT ',
     &                       DM(ITACH2+J-1))
C 
                 DM(PTAUDC+I-1) = DM(DEPUIS+J-1)/1.D -14
C 
               ELSE
C 
                 CALL IMPET ('ON REMPLACE LE PAS DE TEMPS ', I)
                 CALL IMPET ('PAR LE PAS DE TEMPS         ', J)

                 DM(PTAUDC+I-1) = DM(DEPUIS+J-1)/DM(ITACH2+J-1)
C 
               END IF
C 
             ELSE
C 
               CALL MESSAO ( 'C''EST A DEVENIR FOU, AUCUN '//
     &         'PAS DE TEMPS N''EST SIGNIFICATIF ')
C 
               CALL ERREUD (0, ' DANS '//IDPROG)
C 
             END IF
C 
           ELSE
C 
             CALL IMPET ('POUR LE PAS DE TEMPS ', I)
             CALL IMPDT ('LA VALEUR DE LA PUISSANCE DU TAUX VAUT ',
     &                    DABS(DM(DEPUIS+I-1)))
             CALL IMPDT ('MAIS LE TAUX DE CHARGEMENT AU CARRE VAUT ',
     &                    DM(ITACH2+J-1))
C 
             DM(PTAUDC+I-1) = DM(DEPUIS+I-1)/1.D -14
C 
           END IF
C 
         ELSE
C 
C          ON EST DANS UN CAS NORMAL
C 
           DM(PTAUDC+I-1) = DM(DEPUIS+I-1)/DM(ITACH2+I-1)
C 
         END IF
C 
      END DO
C 
      CALL IMPTDT ('INDICATEUR D''INSTABILITE AVANT INTEGRATION ',
     &              DM(PTAUDC), 1, NPICET)
C 
      CALL intgft (TYPE, ADINTE, DM(PTAUDC), FOTEMP(2))
C 
CD    CALL RETOUD(IDPROG)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Determination des periodes en temps.
C     actuellement quatre maxi, dont trois comptees.
C 
C     une non comptee lorsque l'erreur globale est inferieure a 2%
C     une lorsque l'erreur globale est comprise entre 2% et 25%
C     une pour les piquets compris entre 25% et l'instabilite
C     une au dela de l'instabilite
C 
C     le principe est de reperer uniquement les debuts de periode,
C     les fins de periode etant obtenue par -1.
C     Si l'instabilite est comprise entre 2% et 25% on s'arrete a l'instabilite
C     puis on cree une periode au dela de l'instabilite
C 
C     S'il n'y a pas d'instabilite on cree au maxi 3 periodes
C     < 2% , de 2% a  25% , au dela de 25%

      SUBROUTINE DPERIO
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'strategie_calcul.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      INTEGER DEPUIS, AERREU, I, J
C 
      INTEGER INTINS, INT2, APERIO, TABPER, ADSAUT, NSAUT
      INTEGER INSP(2), NPICSP
C 
      DOUBLE PRECISION VALCOM
C 
      LOGICAL LOGINS
C 
CD    LOGICAL  LTRACN, LTRACP
C 
      CHARACTER*3 CARETG
      CHARACTER*6 IDPROG
      INTEGER     AM2LC, ADM2LC
      PARAMETER  (IDPROG='DPERIO')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
      CALL MESSAO ('ON ENTRE DANS '// IDPROG)
      CALL IDENTI (NBETGL, CARETG)
C 
C     Puissance des taux des efforts exterieurs
C 
      CALL ADTBDM ('DEPUIS-'//CARETG, DEPUIS)
C 
C     Erreur globale a l'etape locale
C 
      CALL ADTBDM ('ERREUR-LOC', AERREU)
C 
      CALL POUSME (NPICET+2, APERIO)
C 
      LOGINS = .FALSE.
      I = 0
C 
      DO WHILE ((.NOT. LOGINS) .AND. (I .LT. NPICET))
        IF (DM(DEPUIS+I) .LT. 0.D0) THEN
          LOGINS = .TRUE.
C 
C     intins est l'intervalle ou se produit l'instabilite
C 
          INTINS =  I+1
        END IF
        I = I+1
      END DO
C 
      INT2  = 1
      I     = 1
C 
      DO WHILE ((DM(AERREU+I) .LE. 0.02) .AND. (I .LE. NPICET))
C 
C      int2 est le premier intervalle ou l'erreur est >  0.02
C 
        INT2 = I+1
        I    = I+1
      END DO
C 
      IF (INT2 .EQ. NPICET+1) THEN
C 
        M(APERIO)   = 1
        M(APERIO+1) = 1
        M(APERIO+2) = NPICET+1
C 
C       Le decoupage en temps est termine
C 
      END IF
C 
C     Pour eviter de subdiviser, une idee serait de prendre la racine
C     carree de DM(AERREU+NPICET)
C 
      VALCOM = DM(AERREU+NPICET)/DBLE(NPICET)
C 
      NSAUT = 0
      CALL POUSME (NPICET, ADSAUT)
C 
      DO I = 1, NPICET
        IF ((DM(AERREU+I)-DM(AERREU+I-1)) .GT. 2.D0*VALCOM) THEN
C 
C        i est l'intervalle pour lequel l'acroissement d'erreur
C        est 2x superieur a l'accroissement moyen
C 
          M(ADSAUT+NSAUT) = I
          NSAUT           = NSAUT+1
        END IF
      END DO
C 
      IF (NSAUT .NE. 0)  M(ADSAUT+NSAUT) = NPICET +1
C 
C    Determination du nombre de piquets speciaux et rangement par ordre croissant
C 
      IF (LOGINS) THEN
        IF (INTINS .LT. INT2) THEN
          NPICSP      = 2
          INSP(1)     = INTINS
          INSP(2)     = INT2
          M(APERIO  ) = 2
          M(APERIO+1) = INSP(1)
          M(APERIO+2) = INSP(2)
        ELSE IF (INTINS .GT. INT2) THEN
          NPICSP      = 2
          INSP(1)     = INT2
          INSP(2)     = INTINS
          M(APERIO  ) = 2
          M(APERIO+1) = INSP(1)
          M(APERIO+2) = INSP(2)
        ELSE IF (INTINS .EQ. INT2) THEN
          NPICSP      = 1
          INSP(1)     = INT2
          M(APERIO  ) = 1
          M(APERIO+1) = INSP(1)
        END IF
      ELSE
          NPICSP      = 1
          INSP(1)     = INT2
          M(APERIO  ) = 1
          M(APERIO+1) = INSP(1)
      END IF
C 
      J = NPICSP+1
C 
      IF (NSAUT .NE. 0) THEN
C 
        I = 0
        J = 1
        DO WHILE ((M(ADSAUT+I) .LT. INSP(1)) .AND. (I .LT. NSAUT))
            M(APERIO+J) = M(ADSAUT+I)
                J       = J +1
                I       = I +1
            M(APERIO) = M(APERIO)+1
        END DO
        M(APERIO+J) = INSP(1)
            J       = J +1
C 
        IF (M(ADSAUT+I) .EQ. INSP(1)) I = I+1
C 
        IF (NPICSP .EQ. 2) THEN
C 
          DO WHILE (M(ADSAUT+I) .GT. INSP(1)
     &             .AND. M(ADSAUT+I) .LT. INSP(2)
     &             .AND. (I .LT. NSAUT))
              M(APERIO+J) = M(ADSAUT+I)
                  J       = J +1
                  I       = I +1
              M(APERIO) = M(APERIO)+1
          END DO
          M(APERIO+J) = INSP(2)
              J       = J +1
C 
C         Pour le test qui suit
C 
          INSP(1) = INSP(2)
          IF( M(ADSAUT+I) .EQ. INSP(1) ) I = I+1
        END IF
C 
        DO WHILE (I .LT. NSAUT)
              M(APERIO+J) = M(ADSAUT+I)
                  J       = J +1
                  I       = I +1
                M(APERIO) = M(APERIO)+1
        END DO
C 
      END IF
C 
      M(APERIO+J) = NPICET+1
C 
C     Le decoupage en temps est termine
C 
      IF (LOGINS) THEN
C 
        IF (NSAUT . EQ. 0) THEN
          M(APERIO)   = 2
          M(APERIO+1) = INT2
          M(APERIO+2) = INTINS
          M(APERIO+3) = NPICET+1
        ELSE
C 
          M(APERIO)   = 2+NSAUT
          M(APERIO+1) = INT2
C 
          J = 2
C 
          IF (INTINS .LT. M(ADSAUT)) THEN
            M(APERIO+J) = INTINS
            J           = J +1
            DO I = 0 , NSAUT-1
              M(APERIO+J) = M(ADSAUT+I)
              J = J+1
            END DO
            M(APERIO+J) = NPICET+1
C 
          ELSE
C 
            DO I = 0, NSAUT-1
              IF( M(ADSAUT+I) .LT. INTINS .AND.
     &                M(ADSAUT+I+1) .GE. INTINS) THEN
                M(APERIO+J) = INTINS
                J           = J +1
                M(APERIO+J) = M(ADSAUT+I)
              ELSE
                M(APERIO+J) = M(ADSAUT+I)
              END IF
              J = J+1
            END DO
            M(APERIO+J) = NPICET+1
C 
          END IF
C 
        END IF
C 
C     Si pas d'instabilite
C 
      ELSE
C 
        IF (NSAUT .EQ. 0) THEN
          M(APERIO)   = 1
          M(APERIO+1) = INT2
          M(APERIO+2) = NPICET+1
C 
        ELSE
          M(APERIO)     = 1+NSAUT
          M(APERIO+1)   = INT2
C 
          J = 2
          DO I = 0 , NSAUT-1
            M(APERIO+J) = M(ADSAUT+I)
            J = J+1
          END DO
          M(APERIO+J) = NPICET+1
C 
        END IF
C 
      END IF
C 
C      pour tester la procedure avec tous les pas de temps
C 
C      M(APERIO)   = NPICET
C 
C      DO I =  1, NPICET+1
C 
C        M(APERIO+I) = I
C 
C      END DO
C 
C      pour alliant tant que la procedure n'est pas au point
C 
      IF (UPERIO) THEN
C 
        M(APERIO)   = 1
        M(APERIO+1) = 1
        M(APERIO+2) = NPICET+1
C 
      END IF
C 
      IF (UNINST) THEN
        IF (LOGINS) THEN
           IF (INTINS .EQ.1 .OR. INTINS.EQ. NPICET) THEN
             M(APERIO)   = 1
             M(APERIO+1) = 1
             M(APERIO+2) = NPICET+1
           ELSE
             IF (INTINS .GT. INT2) THEN
               M(APERIO)   = 2
               M(APERIO+1) = INT2
               M(APERIO+2) = INTINS
               M(APERIO+3) = NPICET+1
             ELSE
              M(APERIO)   = 1
              M(APERIO+1) = INT2
              M(APERIO+2) = NPICET+1
             ENDIF
           END IF
        ELSE
           M(APERIO)   = 1
           M(APERIO+1) = INT2
           M(APERIO+2) = NPICET+1
        END IF
      END IF
C 
      CALL GESTEN ('PERIOD -'//CARETG, M(APERIO)+2, TABPER)
C 
      CALL COPITE (M(APERIO)+2, M(APERIO), M(TABPER))
C 
      CALL IMPTET ('PERIOD -'// CARETG, M(TABPER), 1, M(TABPER)+2)
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
