C     Caracterisation des eventuelles zones delaminees
C 
C     On utilise comme arguments :
C 
C     E ...... DINDEL logique : calcul avec ou sans delaminage initial
C     E ...... TYPDEV entier  : 0 pour descro del ini par Fourier
C                               1 pour descro del ini par angle, rayon, d=0 ou 1
C                               2 pour descro del ini par angle, rayon, d=f(rayon)
C     On obtient en sortie :
C 
C     S ...... PREDELAINT tableau des endommagements initiaux
C 
      SUBROUTINE DONDEL
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
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
      INTEGER         AM2LC, ADM2LC, I, J
C 
      INTEGER         ADINT, NINLU, ADCOEF, NBCOEF
      INTEGER         DBCOEF, ADVCOE, DBVCOE
      INTEGER         ADNBCO, DBNBCO, ADDEPA, ADDEPB, ADD0A, ADD0B
C 
      CHARACTER*3     CARINT, CARNUM
C 
      LOGICAL         LECLOG, DINDEL
      INTEGER         TYPDEV, LECINT
C 
C     NOM DES FICHIERS DE DONNEES EVENTUELS
C 
      INTEGER  APM0
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='DONDEL')
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
      IF (VISORT) THEN
C 
         CONTRA(1)  = 'C11'
         CONTRA(2)  = 'C22'
         CONTRA(3)  = 'C12'
         CONTRA(4)  = 'C23'
         CONTRA(5)  = 'C13'
         CONTRA(6)  = 'C33'
C 
         DEFORM(1)  = 'E11'
         DEFORM(2)  = 'E22'
         DEFORM(3)  = 'E12'
         DEFORM(4)  = 'E23'
         DEFORM(5)  = 'E13'
         DEFORM(6)  = 'E33'
C 
         DCONTR(1)  = 'DC11'
         DCONTR(2)  = 'DC22'
         DCONTR(3)  = 'DC12'
         DCONTR(4)  = 'DC23'
         DCONTR(5)  = 'DC13'
         DCONTR(6)  = 'DC33'
C 
         DDEFOR(1)  = 'DE11'
         DDEFOR(2)  = 'DE22'
         DDEFOR(3)  = 'DE12'
         DDEFOR(4)  = 'DE23'
         DDEFOR(5)  = 'DE13'
         DDEFOR(6)  = 'DE33'
C 
      ELSE IF (VISUXY) THEN
C 
         CONTRA(1)  = 'Cxx'
         CONTRA(2)  = 'Cyy'
         CONTRA(3)  = 'Cxy'
         CONTRA(4)  = 'Cyz'
         CONTRA(5)  = 'Cxz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Exx'
         DEFORM(2)  = 'Eyy'
         DEFORM(3)  = 'Exy'
         DEFORM(4)  = 'Eyz'
         DEFORM(5)  = 'Exz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCxx'
         DCONTR(2)  = 'DCyy'
         DCONTR(3)  = 'DCxy'
         DCONTR(4)  = 'DCyz'
         DCONTR(5)  = 'DCxz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DExx'
         DDEFOR(2)  = 'DEyy'
         DDEFOR(3)  = 'DExy'
         DDEFOR(4)  = 'DEyz'
         DDEFOR(5)  = 'DExz'
         DDEFOR(6)  = 'DEzz'
C 
      ELSE
C 
         CONTRA(1)  = 'Crr'
         CONTRA(2)  = 'C00'
         CONTRA(3)  = 'Cr0'
         CONTRA(4)  = 'C0z'
         CONTRA(5)  = 'Crz'
         CONTRA(6)  = 'Czz'
C 
         DEFORM(1)  = 'Err'
         DEFORM(2)  = 'E00'
         DEFORM(3)  = 'Er0'
         DEFORM(4)  = 'E0z'
         DEFORM(5)  = 'Erz'
         DEFORM(6)  = 'Ezz'
C 
         DCONTR(1)  = 'DCrr'
         DCONTR(2)  = 'DC00'
         DCONTR(3)  = 'DCr0'
         DCONTR(4)  = 'DC0z'
         DCONTR(5)  = 'DCrz'
         DCONTR(6)  = 'DCzz'
C 
         DDEFOR(1)  = 'DErr'
         DDEFOR(2)  = 'DE00'
         DDEFOR(3)  = 'DEr0'
         DDEFOR(4)  = 'DE0z'
         DDEFOR(5)  = 'DErz'
         DDEFOR(6)  = 'DEzz'
C 
      ENDIF
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
CD    CALL WLKBCD (IDPROG)
C 
C     INITIALISATION DE LA PRECISION
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C 
C     ENTREE DANS LA SEQUENCE D'INITIALISATION DES DONNEES
C 
      CALL LECSEQ ('ZONE-DELAMI',
     &             'POUR CARACTERISER LES ZONES DELAMINEES')
C 
      DINDEL = LECLOG ('DINDEL')
      TYPDEV = LECINT ('TYPDEV')
      NUPATE = LECINT ('NUPATE')
C 
      IF (DINDEL) THEN
C             
        CALL MESSAO ('INDIQUER LES NUMEROS D''INTERFACES PREDELAMINEES :
     &               \ATTENTION IL EST IMPERATIF DE DONNER LES NUMEROS 
     &               \D''INTERFACE PAR ORDRE CROISSANT 
     &               \(i.e. A PARTIR DU BAS)')
C 
        CALL POUSME (NBINT, APM0)
C 
        CALL LECLEN (M(APM0), NBINT, NINLU)
        CALL GESTEN ('NUMINT-DEL', NINLU, ADINT)
        CALL COPITE (NINLU, M(APM0), M(ADINT))
C 	
	IF (TYPDEV .EQ. 0) THEN
C 
C         TABLEAU DES NOMBRES DE COEFFICIENTS PAR INTERFACE
C 
          CALL GESTEN ('NBCOINTDEL', NINLU, ADNBCO)
          DBNBCO = ADNBCO
C 
          CALL POUSME (NINLU*NBMAT, ADCOEF)
          DBCOEF = ADCOEF
C 
          CALL POUSMD (NINLU*NTETA, ADVCOE)
          DBVCOE = ADVCOE
C 
          DO I =  0 , NINLU-1
C 
            CALL IDENTI (M(ADINT+I), CARINT)
C 
            CALL MESSAO ('POUR L''INTERFACE '//carint//' DONNER LES 
     &                   \NUMEROS DE DEVELOPPEMENTS EN SERIES DE FOURIER
     &                   \ PRINCIPE : +n => r =  cos(n0)
     &                   \            -n => r =  sin(n0)        ' )
C 
            CALL LECLEN (M(ADCOEF), NTDSFG, M(ADNBCO))
C 
            DO J = 0, M(ADNBCO)-1
C 
              CALL IDENTI (M(ADCOEF+J), CARNUM)
              CALL MESSAO ('DONNEE, POUR L''INTERFACE '//CARINT//
     &                     ' ET POUR LE NUMERO DE DEVELOPPEMENT, '
     &                     //CARNUM// ' DU COEFFICIENT.' )
C 
              CALL LECLDP (DM(ADVCOE), NTETA, NBCOEF)
              ADVCOE     = ADVCOE+NBCOEF
C 
            END DO
C 
            ADCOEF = ADCOEF+M(ADNBCO)
            ADNBCO    = ADNBCO+1
C 
          END DO
C 
C 
CD        CALL IMPTET ('NUMEROS DES INTERFACES PREDELAMINEES ',
CD                      M(ADINT), 1, NINLU)
C 
CD        CALL IMPTET ('NOMBRE DE COEF PAR INTERFACE PREDELAMINEE ',
CD                      M(DBNBCO), 1, NINLU)
C 
          CALL GESTEN ('NUCOINTDEL', ADCOEF-DBCOEF, ADDEPA)
          CALL COPITE (ADCOEF-DBCOEF, M(DBCOEF), M(ADDEPA))
C 
CD        CALL IMPTET ( NUMEROS DES COEF PAR INTERFACE PREDELAMINEE',
CD                    M(ADDEPA), 1, ADCOEF-DBCOEF)
C 
          CALL GESTDP ('VLCOINTDEL', ADVCOE-DBVCOE, ADDEPA)
          CALL COPITD (ADVCOE-DBVCOE, DM(DBVCOE), DM(ADDEPA))
C 
CD        CALL IMPTDT ('VALEUR DES COEF ', DM(ADDEPA), 1, ADVCOE-DBVCOE)
C 
          ADDEPA = ADDEPA + NTETA
C 
C       Donnee angle par angle, cas d'une evolution d=0 ou 1
C  
        ELSE IF (TYPDEV .EQ. 1) THEN
C 	
	  CALL GESTDP ('RAYONSDELA', NINLU*NTETA, ADDEPA)
C 
          ADD0A = ADDEPA
C 
	  DO I =  0, NINLU-1
	  
            CALL IDENTI (M(ADINT+I), CARINT)
C 
1002        CONTINUE        
            CALL MESSAO ('POUR L''INTERFACE '//carint//' DONNEZ ANGLE 
     &                   \PAR ANGLE LES NTETA VALEURS DES RAYONS DE LA 
     &                   \COURBE DELIMITANT LE DELAMINAGE ')	  
C 
            CALL LECLDP (DM(ADDEPA), NTETA, NBCOEF)
C 	    
	    IF (NBCOEF .NE. NTETA) THEN
               CALL MESSAO ('VOUS N''AVEZ PAS DONNE LE BON
     &                      \NTETA COEFF; RECOMMENCEZ ')
      
               GOTO 1002 
C 	      
	    END IF  
C 
            ADDEPA = ADDEPA + NTETA
C 
          END DO
C 
C       Donnee angle par angle, cas d'une evolution adoucie de d en fonction du rayon
C 
        ELSE IF (TYPDEV .EQ. 2) THEN
C 
	  CALL GESTDP ('RAYONSDELA', NINLU*NTETA, ADDEPA)
	  CALL GESTDP ('RAYONSENDO', NINLU*NTETA, ADDEPB)
C 
          ADD0A = ADDEPA
	  ADD0B = ADDEPB
C 
	  DO I =  0, NINLU-1
	  
            CALL IDENTI (M(ADINT+I), CARINT)
C 
1003        CONTINUE        
            CALL MESSAO ('POUR L''INTERFACE '//carint//' DONNEZ ANGLE
     &                   \PAR ANGLE LES NTETA VALEURS DES RAYONS DE LA
     &                   \COURBE DELIMITANT LE DELAMINAGE')	  
C 
            CALL LECLDP (DM(ADDEPA), NTETA, NBCOEF)
C 
            CALL MESSAO ('POUR L''INTERFACE '//carint//' DONNEZ ANGLE
     &                   \PAR ANGLE LES NTETA VALEURS DES RAYONS DE LA
     &                   \COURBE DELIMITANT LA ZONE ENDOMMAGEE')	  
C 	    
            CALL LECLDP (DM(ADDEPB), NTETA, NBCOEF)
C 
	    IF (NBCOEF .NE. NTETA) THEN
C 
              CALL MESSAO ('VOUS N''AVEZ PAS DONNE LE BON
     &                     \NTETA COEFF; RECOMMENCEZ ')
      
              GOTO 1003 
C 	      
	    END IF  
C 
	  ADDEPA = ADDEPA + NTETA
	  ADDEPB = ADDEPB + NTETA
C   
          END DO
C 
        END IF 
C 
      END IF
C 
      IF (.NOT. DINDEL) THEN
C 
        CALL GESTEN ('NUMINT-DEL', 0, ADINT) 
        CALL GESTEN ('NBCOINTDEL', 0, ADNBCO)
        CALL GESTEN ('NUCOINTDEL', 0, ADDEPA)
        CALL GESTDP ('VLCOINTDEL', 0, ADDEPA)
	CALL GESTDP ('RAYONSDELA', 0, ADDEPA)
	CALL GESTDP ('RAYONSENDO', 0, ADDEPA)
C 
      END IF 
C 
      IF (TYPDEV .EQ. 1) THEN
C 
        CALL IMPTDT ('RAYONS ZONE DELAMINEE  ', DM(ADD0A), 1,
     &                NINLU*NTETA)
      END IF
C 
      IF (TYPDEV .EQ. 2) THEN
C 
CD      CALL IMPTDT ('RAYONS ZONE DELAMINEE  ', DM(ADD0A), 1,
CD   &                NINLU*NTETA)
CD      CALL IMPTDT ('RAYONS ZONE ENDOMMAGEE ', DM(ADD0B), 1,
CD   &                NINLU*NTETA)
      END IF
C 
      CALL TRADEL (TYPDEV)
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
C     Cette routine traite le predelaminage
C 
C     On envoie comme argument :
C 
C     E ...... TYDEV  0 pour une description en series de Fourier du delaminage initial
C                     1 pour une description en TETA du rayon de la zone delaminee, d=0 ou 1
C                     2 pour une description en TETA du rayon de la zone delaminee + raccordement
C 
      SUBROUTINE TRADEL (TYPDEV)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      INTEGER  TYPDEV
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  A, DANGL, MULT, RAYDEL, RAYEND
      DOUBLE PRECISION  RLOC, RLOCI, RAYONC
C 
      INTEGER    TETA, NUCOL
      INTEGER    NUINT, L
      INTEGER    AM2LC, ADM2LC, HINT, KINT
      INTEGER    X, ADPOPG
      INTEGER    DBPOPG
C 
      INTEGER    ADTETA, NINLU, ADINT, ADNBCO
      INTEGER    NBNUM , ADNUM, NBVAL, ADVAL, DBVAL
      INTEGER    LONRES, ENDOM, ENDODO, DBNUM
      INTEGER    NBCOE, COEFF, NUDEV
C 
      LOGICAL    LOGDEL
      INTEGER    LONRAY, LONEND, ADDRAY , ADDEND, DEBRAY, DEBEND
C 
CD    LOGICAL    LTRACN, LTRACP
C 
      CHARACTER*6 IDPROG
      PARAMETER (IDPROG='TRADEL')
C 
CD    CALL WLKBCD (IDPROG)
C 
      CALL ENPOUB (AM2LC, ADM2LC)
C 
C -----------------------------------------------------------------------
C 
      CALL ADTBDM ('POIDS-PGCR', ADPOPG )
C 
C     POUR ALLER LIRE LES POINTS DE GAUSS
C 
      HINT = XINTEG*(XINTEG-1)/2
      KINT = YINTEG*(YINTEG-1)/2
C 
      CALL ADTBDM ('ANGLES-GEO', ADTETA)
      CALL INFOEN ('NUMINT-DEL', ADINT, NINLU)
C 
      IF (TYPDEV .EQ. 0) THEN
C 
        CALL ADTBM ('NBCOINTDEL', ADNBCO)
        CALL INFOEN ('NUCOINTDEL', ADNUM, NBNUM)
        DBNUM = ADNUM
C 
        CALL INFODP ('VLCOINTDEL', DBVAL, NBVAL)
C 		
      END IF	
C 
      DBPOPG = ADPOPG+NGAU1
C 
      LONRES = 3*NGAU2*NTETA
	
      IF (TYPDEV .EQ. 1) THEN
C 
        CALL INFODP ('RAYONSDELA', DEBRAY, LONRAY)
        ADDRAY = DEBRAY	
C 
      END IF				
C 	
      IF (TYPDEV .EQ. 2) THEN
C 
        CALL INFODP ('RAYONSDELA', DEBRAY, LONRAY)
	CALL INFODP ('RAYONSENDO', DEBEND, LONEND)
        ADDRAY = DEBRAY
	ADDEND = DEBEND
C 
      END IF
C 				
      CALL GESTDP ('PREDELAINT', LONRES, ENDOM)
      ENDODO = ENDOM
C 
C     BOUCLE SUR LES INTERFACES
C 
      DO NUINT= 1, NBINT
C 
        LOGDEL = .FALSE.
C 
        DO L = 0, NINLU-1
C 
          IF (NUINT .EQ. M(ADINT+L)) THEN
C 
            LOGDEL = .TRUE.
C 
CD          CALL IMPET ('NUMERO D''INTERFACE DANS '//IDPROG, NUINT)
C 
            IF (TYPDEV .EQ. 0) THEN
              NBCOE  = M(ADNBCO)
              ADNBCO = ADNBCO+1
C 
C            ELSE 
C 
C 	         GOTO 1000
C 
            END IF
C 
          END IF
C 
C 1000      CONTINUE
C 
        END DO
C 
C       BOUCLE TEST (INTERFACE DELAMINEE ?)
C 
	IF (LOGDEL) THEN
C 
C         BOUCLE i SUR LES COLONNES
C 
          DO NUCOL = 1, NBCOL
C 
            CALL VALRAY (NUCOL, RAYONC, A)
C 
C           BOUCLE ii SUR LES POINTS DE GAUSS EN X
C 
            DO X = 1, XINTEG
C 
              RLOC  = RAYONC + A*GAUSS(HINT+X)
C 
C             BOUCLE iii SUR LES ANGLES
C 
              DO TETA = 1, NTETA
C 
                DANGL = DM(ADTETA+TETA-1)
C 
C               BOUCLE iv SUR LES TYPES DE DESCRO DU PREDEL
C 
		IF (TYPDEV .EQ. 0) THEN
C 
		  RAYDEL = 0.D0
                  ADVAL = DBVAL
                  ADNUM = DBNUM
C 
                  DO COEFF = 1, NBCOE
C 
                    NUDEV = M(ADNUM)
                    ADNUM = ADNUM+1
C 
                    IF (NUDEV .LT. 0) THEN
                      MULT = DSIN(DBLE(-NUDEV)*DANGL)
                    ELSE
                      MULT = DCOS(DBLE(NUDEV)*DANGL)
                    ENDIF
C 
                    RAYDEL = MULT*DM(ADVAL)+RAYDEL
                    ADVAL  = ADVAL+1
C 
                  END DO
C 
                  IF (RAYDEL .LT. RLOC) THEN
C 
		    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
C 
                  ELSE
C 
		    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
C 
		  END IF
C 
                ELSE IF (TYPDEV .EQ. 1) THEN
C 
                  RAYDEL = DM(ADDRAY)
		  ADDRAY = ADDRAY + 1
C 
                  IF (RAYDEL .LT. RLOC) THEN
C 
		    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
C 
                  ELSE
C 
		    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
C 
		  END IF
C 
                ELSE IF (TYPDEV .EQ. 2) THEN
C 
                  RAYDEL = DM(ADDRAY)
		  RAYEND = DM(ADDEND)
		  ADDRAY = ADDRAY + 1
		  ADDEND = ADDEND + 1
C 
                  IF (RAYEND .LT. RLOC) THEN
C 
		    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 0.D0
                    ENDOM     = ENDOM+1
C 
                  ELSE IF (RLOC .LT. RAYDEL) THEN
C 
		    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
                    DM(ENDOM) = 1.D0
                    ENDOM     = ENDOM+1
C 
                  ELSE
C 
		    RLOCI     = RLOC - RAYDEL
               	    DM(ENDOM) = .5D0*(1.D0+(DCOS((PI*RLOCI)/
     &                          (RAYEND - RAYDEL))))
		    ENDOM     = ENDOM+1
          	    DM(ENDOM) = .5D0*(1.D0+(DCOS((PI*RLOCI)/
     &                          (RAYEND - RAYDEL))))
		    ENDOM     = ENDOM+1
         	    DM(ENDOM) = .5D0*(1.D0+(DCOS((PI*RLOCI)/
     &                          (RAYEND - RAYDEL))))
 		    ENDOM     = ENDOM+1

                  END IF
C 		  
C               FIN DE BOUCLE v SUR LES TYPES DE DESCRO DU PREDEL
C 	
                END IF
C 
C             FIN DE BOUCLE iv SUR LES ANGLES
C 
	      END DO
C 
C           FIN DE BOUCLE iii SUR LES POINTS DE GAUSS EN X
C 
              ADDRAY = DEBRAY
	      ADDEND = DEBEND
	    END DO
C 
C         FIN DE BOUCLE ii SUR LES COLONNES
C 
            ADDRAY = DEBRAY
	    ADDEND = DEBEND
	  END DO
C 
C       FIN DE BOUCLE TEST (INTERFACE DELAMINEE?)
C 
        END IF
C 
        DBNUM = DBNUM+NBCOE
        DBVAL = DBVAL+NBCOE
C 
C     FIN DE BOUCLE SUR LES INTERFACES
C 
        DEBRAY = DEBRAY + NTETA
	DEBEND = DEBEND + NTETA
        ADDRAY = DEBRAY
	ADDEND = DEBEND
      END DO
C 
      CALL SOPOUB (AM2LC, ADM2LC)
C 
CD    CALL RETOUD (IDPROG)
C 
      RETURN
      END
