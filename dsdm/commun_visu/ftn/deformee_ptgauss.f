        Subroutine triten (lontab, tab, tabtri)
C SUBROUTINE TRIANT LES VALEURS DU TABLEAU ENTIER TAB ET LES METTANT PAR ORDRE c CROISSANT DaNS TABTRI       
        integer lontab, tab(lontab), tabtri(lontab)
c
        do  i=1,lontab-1
        do  j=i+1,lontab
           if (TAB(j).LT.TAB(i)) then
            TABTRI(i)= TAB(j) 
           TABTRI(j)= TAB(I)
          endif
        end do
       end do
        end
C
C ROUTINES UTILE POUR RETROUVER LES CHOSES DANS VISUAL DEPPLA MAIPAR DTEXPR
C DANS DELAMI VCHAPT VCUPTE, VDERTE ,VRTSM, VPCGTE RCHARE CATELE LIB1/NUMEROTATION

C     Ecriture dans le repertoire resultats des fichiers exploitables
C     par PARAVIEW en terme de deplacement aux les points de Gauss fictifs
C     associés au maillage créé par MAIPAR 
C c'est à dire au maillage POUR LE MOMENT NE NUMEROTATION COUCHE
C de toutes les points de Gauss dans les bandes du calcul DSDM
C range (UX, nteta, ---rangement ligne ncolpg, ncouchpg )
C(UY, nteta, ncolpg, ncouchpg )
C(UZ, nteta, ncolpg, ncouchpg )
C npicet 
C 
C le principe on recupere par VCHAPT les deplacements en coordonnees polaire stockes
C (nteta,numoretation couche ou colonne, pas de temps par pas de temps)
C on inverse nteta et numerotation couche ou colonne pour avoir tous les deplacements
C dans une bande en coordonnees polaires. 
C a partir des deplacements aux noeuds de chaque elements on interpole au moyen des
C fonctions de base elementaires (attention Ã  la numeroration particuliere de ces 
C fonctiosn de base pour pouvoir ensuite calculer en coordonnees polaires
C         U=N1*u1+N3*u2+N5*u3+N7*u4+N9*u5+N11*u6
C          +N2*u'1+N4*u'2+N6*u'3+N8*u'4+N10*u'5+N12*u'6
C On effectue la rotation pour se ramener en coordonnees cartesiennes
C on ne stocke que Ux,UY,UZ on adonc un tableau
C (UX,UY,UZ, nbnoeud bande, nteta) par piquet de temps
C  on extrait pour avoir (Ux, nteta, nbnoued bande)
C puis                   (Uy, nteta, nbnoued bande)
C puis                   (Uz, nteta, nbnoued bande)
C 
C
      SUBROUTINE DEPRPG
C 
C Subroutine DEPRPG : crée les deplacements reels au points de Gauss associés à maipar
C                      pour tout les points de Gauss 
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      include 'cominc.h'
      include 'incl_visu_mail.h'
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C NNOERF NOMBRE DE NOEUDS SUR L'ELEMENT DE REFERENCE DISCRETISE POUR CONTENIR TOUT LES POINTS DE GAUSS
C
      INTEGER           NNOERF, DEBX, DEBY, VAL,LVFBNG,ADFBNG,DEBUT,P,Q
      DOUBLE PRECISION  X, Y , LONGX, LONGY
      DOUBLE PRECISION  N1,N2, N3,N4, N5, N6, N7,N8,N9,N10,N11, N12

      INTEGER           LONDEN , ADRXYT
      INTEGER           LONDEG , AGRXYT      
C    
C   nuetgl pour visu ? 11/01/2013  en resolution NBETGL
       CHARACTER*3 CARETG
C
C      Pour les fcts du temps
      INTEGER FTREE7, FTEPMX, FIFTEP, FICHEP, DBFTEP, DBCHEP
      INTEGER DEPREE, DEPCH0
D     INTEGER LONDEP, NBDEP      
C       
      INTEGER ADANGL, ADER, ADTETA, I
      INTEGER LONDEL, AVANT , DEB
      INTEGER  ADRU, ADRUP,ADRV,ADRVP,ADRW,ADRWP
      INTEGER L, K, J , TEST , DECOU, DECOL 
      INTEGER DECALE , TABNIV(2),LONG, NUY, NUX, NBX, COOR, FNPG, POSTG    
      INTEGER LONXYC, AXYC, ACXY , ADEPPA, POSI, DEP, TETA,POSTDC     
      INTEGER LXYC1, LXYC2 , VALDEP, NBNBAG , DEBUR
C 28/08/2013 DEBUT pr DEBDEP et  par DEBUN  pour les coordonnees    
      INTEGER DEBUN , DEBDEP, INDIC
C 29/08/2013 test sur la coherence des deplacemnts imposes sur le bord exterieur a teta fixe en polaire
D      INTEGER ADURZ, ADVTZ, ADWZ, DEBURZ, DEBVT,DEBW, DECDEP, TESTDEP 
D      INTEGER  ADUYGA, ADUXGA, ADWGA, DADUX, DADUY, DADW,LONEPG

C       
D     INTEGER LON         
C
      INTEGER            AM2LC, ADM2LC
C Pour ecriture du fichier deplacement au maillage pt de Gauss
C
      INTEGER IUNIC
      CHARACTER*30  NOMFI1
      INTEGER       STATUS       
C 
      CHARACTER*9 NOMDIR 
      CHARACTER*6 IDPROG 
      PARAMETER (IDPROG='DEPRPG')
C 
D      CALL WLKBCD (IDPROG)
C 
C -----------------------------------------------------------------------
C
      CALL ENPOUB (AM2LC,ADM2LC)
C 
      NOMDIR = 'resultats'
C
      CALL GSPOUD(XINTEG+1, DEBX)
      CALL GSPOUD(YINTEG+1, DEBY)
C 
C 
      LONGX=2.D0/DBLE(XINTEG)
      LONGY=2.D0/DBLE(YINTEG)      
C Creation des coordonnees des noeuds fictifs de l'elements de reference aux point de Gauss en x et y
C par exemple enx pour 2 pt de Gauss -1, -1/3,1/3 et 1 on pourrait mettre la vrai position des pts de Gauss
CLe maillage cree encadre les points de Gauss si on a xinteg pt de gauss il faut XINTEG +1 point pour les encadrer Cen X de mÃªme en y 
C
      DEBUT= DEBX
      X= -1.D0
      DO I=0,XINTEG
        DM(DEBUT)= X
        X=  X+LONGX
        DEBUT = DEBUT+1
        END DO
      DEBUT = DEBY  
      Y= -1.D0    
      DO J=0, YINTEG
           DM(DEBUT)= Y
           Y=Y+LONGY
           DEBUT = DEBUT +1
      END DO 
C 
C ----------------------------------------------------------------------
C     Creation d'un tableau  pour ranger les valeurs
C     des fonctions elementaires ASSOCIÉS AUX DEPLACEMENTS aux differents points de gauss
C     1ere adresse:ADLC1
C ----------------------------------------------------------------------

C RAPPELS
C 
C         Dans le tableau appele TABELE 
C 
C         la numerotation ligne si num=1;<=====>2*(nbcol+2).le.3*(nbcou+1)
C                          '''''                2*nbcol     .le.3*nbcou-1
C         la numerotation colonne si num=2;<===>2*(nbcol+2).gt.3*(nbcou+1)
C                          ''''''               2*nbcol    .gt.3*nbcou-1
C -----------------------------------------------------------------------
C         dans la numerotation ligne on a:
C 
C         l'element de la couche q de la colonne p a pour numero:
C                     (Q-1)NBCOL+P
C         et pour numero de premier noeud:
C                      3(NBCOL+1)(Q-1)+P
C         l'element de l'interface r de la colonne p a pour numero:
C                      NEL1 +(R-1)NBCOL+P
C         et pour numero de premier noeud:
C                      (3R-1)(NBCOL+1)+P
C -----------------------------------------------------------------------
C         dans la numerotation colonne on a:
C 
C         l'element de la couche q de la colonne p a pour numero:
C                     (P-1)NBCOU+Q
C         et pour numero de premier noeud:
C                      3(NBCOU)(P-1)+3Q-2
C         l'element de l'interface r de la colonne p a pour numero:
C                      NEL1 +(P-1)(NBCOU-1)+R
C         et pour numero de premier noeud:
C                      3(NBCOU(P-1)+R)
C -----------------------------------------------------------------------
C         on range sequentiellement les informations concernant
C         les elements comme suit:
C 
C         element couche(type 1)
C 
C         - type
C         - numero de couche
C         - numero de colonne
C         - numero de premier noeud    numerotation 4--------------3
C         - numero de second noeud                  5--------------6
C         - numero de troisieme noeud  choisie      1--------------2
C         - numero de quatrieme noeud
C         - numero de cinquieme noeud
C         - numero de sixieme noeud
C 
C         element interface(type 2)
C 
C         - type
C         - numero d' interface
C         - numero de colonne
C         - numero de premier noeud
C         - numero de second noeued     numerotation 4------- 3
C         - numero de troisieme noeud   choisie      1--------2
C         - numero de quatrieme noeud
C ----------------------------------------------------------------------
C
C      P=XINTEG
C      Q=YINTEG
      NNOERF= (XINTEG+1)*(YINTEG+1)    
C
C    Nombre de valeurs des fonctions de base sur au noeud de l'element de reference  du maillage point de Gauss 
      LVFBNG = 12*NNOERF
C     TABLEAU DES VALEURS DES FONCTIONS DE BASE AUX NOEUDS DU MAILLAGE GAUSS SUR L'ELEMENT DE REFERENCE
      CALL GSPOUD(LVFBNG,ADFBNG)
      VAL=ADFBNG-1
C Pour pouvoir ensuite calculer
C         U=N1*u1+N3*u2+N5*u3+N7*u4+N9*u5+N11*u6
C          +N2*u'1+N4*u'2+N6*u'3+N8*u'4+N10*u'5+N12*u'6
C
C on determine les valeurs des fonctions de base aux noeuds du maillage point de gauss
C
      DO J=0,YINTEG
C    cordonnee en Y du point considere sur l'element de reference  
        Y = DM(DEBY+J)
        DO I=0,XINTEG
C         cordonnee en X du point considere sur l'element de reference  
          X        = DM(DEBX+I)
C
          DM(VAL+1)=N1(X,Y)
C 
          DM(VAL+2)=N3(X,Y)
C 
          DM(VAL+3)=N5(X,Y)
C 
          DM(VAL+4)=N7(X,Y)
C 
          DM(VAL+5)=N9(X,Y)
C 
          DM(VAL+6)=N11(X,Y)
C
          DM(VAL+7)=N2(X,Y)
C 
          DM(VAL+8)=N4(X,Y)
C 
          DM(VAL+9)=N6(X,Y)
C 
          DM(VAL+10)=N8(X,Y)
C 
          DM(VAL+11)=N10(X,Y)
C 
          DM(VAL+12)=N12(X,Y)
C 

c
D         CALL IMPDT(' POUR X', X)
D         CALL IMPDT(' POUR Y', Y)
D         CALL IMPTDT('6 FCTS 6 DERIV AU POINT DE GAUSS',DM(VAL+1),1,12)         
C 
          VAL       =VAL +12 
        END DO      
C   
      END DO
C    nuetgl pour visu ? 11/01/2013  CALL IDENTI( NUETGL ,  CARETG)
      CALL IDENTI( NBETGL ,  CARETG)
      
      CALL INFODP('FT-EPS-'//CARETG,FTREE7,FTEPMX)
C -
      FTEPMX = FTEPMX/NPICET
C -   
      DBFTEP  = 1
      FIFTEP  = FTEPMX
C -
C -
D      CALL IMPET(' DEBUT DES FONCTIONS DU TEMPS'//
D    &     ' EN DEPLACEMENT '//IDPROG, DBFTEP )
D     CALL IMPET('   FIN DES FONCTIONS DU TEMPS'//
D    &     ' EN DEPLACEMENT '//IDPROG, FIFTEP )
C -
C -                                  
      CALL ADTBDM('DEPLA-ADMI' , DEPCH0)
D     CALL INFODP('DEPLA-ADMI' , DEPCH0, LONDEP) 
D     NBDEP= LONDEP/(NDDL*NTETA)  
D     CALL IMPET('LONGUEUR DEP_ADMI'//IDPROG, LONDEP)
D     CALL IMPET('Nb Depla ADMI maximum'//IDPROG, NBDEP)
D     CALL IMPET('Nb Depla ADMI NBDEPR'//IDPROG, NBDEPR)       
c -
      CALL ADTBDM ('ANGLES-GEO', ADANGL)
C 27/08/2913      
      ADANGL= ADANGL-1
C -      
      CALL GSPOUD ( 2*NTETA , ADER)
      ADER   = ADER-1
      ADTETA = ADER+NTETA
C       
      DO I = 1 ,NTETA
         DM( ADER+I  ) = DCOS( DM( ADANGL+I) )
         DM( ADTETA+I) = DSIN( DM( ADANGL+I) )
      END DO
C - 
D     CALL IMPTDT('valeur des cosinus',DM(ADER+1),1,NTETA)
D     CALL IMPTDT('valeur des sinus', DM(ADTETA+1),1,NTETA)
      TABNIV(1)  = NTETA
      TABNIV(2)  = NDDL

C-          --------->
C POUR LA NUMEROTATION COUCHE (en admettant que l'on note les ddl croissant    
C Principe pour sauter une ligne ajouter 6*(NBCOL+1)*NTETA
C          pour passer au noeud de la ligne d'au dessus idem
C          pour passer au noeud ˆ droite 6*NTETA
C 
C POUR PASSER DU PREMIER DDL :
C             DU PREMIER NOEUD AU NOEUD 2  DECALAGE 6 (NBDDL PAR NOEUD) : 6
C POUR PASSER DU PREMIER AU NOEUD 3  DECALAGE  6* (NBCOL+1) *2  +  6
C POUR PASSER DU PREMIER AU NOEUD 4  DECALAGE  6* (NBCOL+1) *2 
C POUR PASSER DU PREMIER AU NOEUD 5  DECALAGE  6*(NBCOL+1)
C POUR PASSER DU PREMIER AU NOEUD 6  DECALAGE  6*(NBCOL+1)+ 6
C 
      CALL GSPOUE(5, DECALE) 
      M(DECALE   )= 6 
      M(DECALE+1 )= 12* (NBCOL+1) +  6   
      M(DECALE+2 )= 12* (NBCOL+1)  
      M(DECALE+3 )= 6*(NBCOL+1)  
      M(DECALE+4 )= 6*(NBCOL+1)+ 6
C
D     CALL IMPTET('decalage pour la notation elementaire',M(DECALE),1,5)      
C      
C RAPPELS
C 
C         Dans le tableau appele TABELE 
C 
C         la numerotation ligne si num=1;<=====>2*(nbcol+2).le.3*(nbcou+1)
C                          '''''                2*nbcol     .le.3*nbcou-1
C         la numerotation colonne si num=2;<===>2*(nbcol+2).gt.3*(nbcou+1)
C                          ''''''               2*nbcol    .gt.3*nbcou-1
C -----------------------------------------------------------------------
C         dans la numerotation ligne on a:
C 
C         l'element de la couche q de la colonne p a pour numero:
C                     (Q-1)NBCOL+P
C         et pour numero de premier noeud:
C                      3(NBCOL+1)(Q-1)+P
C         l'element de l'interface r de la colonne p a pour numero:
C                      NEL1 +(R-1)NBCOL+P
C         et pour numero de premier noeud:
C                      (3R-1)(NBCOL+1)+P
C -----------------------------------------------------------------------
C         dans la numerotation colonne on a:
C 
C         l'element de la couche q de la colonne p a pour numero:
C                     (P-1)NBCOU+Q
C         et pour numero de premier noeud:
C                      3(NBCOU)(P-1)+3Q-2
C         l'element de l'interface r de la colonne p a pour numero:
C                      NEL1 +(P-1)(NBCOU-1)+R
C         et pour numero de premier noeud:
C                      3(NBCOU(P-1)+R)
C -----------------------------------------------------------------------
C         on range sequentiellement les informations concernant
C         les elements comme suit:
C 
C         element couche(type 1)
C 
C         - type
C         - numero de couche
C         - numero de colonne
C         - numero de premier noeud    numerotation 4--------------3
C         - numero de second noeud                  5--------------6
C         - numero de troisieme noeud  choisie      1--------------2
C         - numero de quatrieme noeud
C         - numero de cinquieme noeud
C         - numero de sixieme noeud
C 
C         element interface(type 2)
C 
C         - type
C         - numero d' interface
C         - numero de colonne
C         - numero de premier noeud
C         - numero de second noeued     numerotation 4------- 3
C         - numero de troisieme noeud   choisie      1--------2
C         - numero de quatrieme noeud
C ----------------------------------------------------------------------
C

C
C TAILLE DU TABLEAU DES DEPLACEMENT AUX NOEUD DU MAILLAGE REEL A TOUT PIQUET DE TEMPS 
C 
      LONDEN= NTETA*NDDL*NPICET
C
      CALL GSPOUD(LONDEN , ADRXYT)
C      
C Valeur des deplacement aux noeuds du maillage "reel" en coordonnees cartesiennes 
CNOMBRE DE VALEURS ASSOCIES AU DEPLACEMENT ( suivant X, Y, Z) POUR LE MAILLAGE DEFINIT CPOUR CONTENIR TOUS LES POINTS DE GAUSS *NPICET
C Attention On ne calcule que les valeurs du deplacement aux noeuds => 3 par les valeurs des derivees en r 
C
      NBNBAG= (NBCOL*XINTEG+1)*NBCOU*(YINTEG+1)
      LONDEG=NTETA*3*NBNBAG*NPICET
C
D      CALL GSPOUD(3*NBNBAG, ADUXGA)
D      ADUYGA= ADUXGA+NBNBAG
D      ADWGA= ADUYGA+NBNBAG
C
      CALL GSPOUD(LONDEG , AGRXYT)  
C-          
      POSI  =AGRXYT
      DEPREE= ADRXYT     
C
C RESERVATION DE PLACE ET MISE A ZERO DU TABLEAU DES VALEURS DES DEPLACEMENTS AUX
C NOEUDS DU MAILLAGE
C  CONTENANT LES POINTS DE GAUSS POUR TOUT LES PIQUETS DE TEMPS
C
C   Nombre de ddl sans les derivees dans l'element de reference avec les points de Gauss      
       LONXYC =3* NNOERF
       CALL GSPOUD(LONXYC,AXYC)
       CALL GSPOUD(LONXYC,ACXY) 
       CALL POUSMD(NDDL, ADEPPA)
       CALL POUSMD (12, VALDEP)  
C- 28/08 2013 pour tester le nombre de valeur construite pour les dp maillage pt de gauss
      INDIC=0
C  
c 
C 29/08/2013 test sur la coherence des deplacemnts imposes sur le bord exterieur a teta fixe en polaire
D      CALL GSPOUD(9*NBCOU,ADURZ)
D      ADVTZ=ADURZ+3*NBCOU
D      ADWZ=ADVTZ+3*NBCOU    
C
    
      DO I=1,  NPICET
C -
C 28/08/2013         POSTDC=POSI
C
         POSTDC=AGRXYT+(I-1)*(NBCOL*XINTEG+1)*NBCOU*(YINTEG+1)*3*NTETA
C
         CALL VCHAPT(
C -
C -        Valeur des champs admissible au pas de tepms ici les deplacements en
C           coordonnées polaires
C -                                             
C      .... on envoie
     &    FTEPMX , NDDL*NTETA , DBFTEP , FIFTEP , DBFTEP , FIFTEP ,
     &    DM(DEPCH0) , DM(FTREE7)  , I ,
C      .... on recupere
     &      DM(DEPREE)  )
D          CALL IMPET('au pas de temps'//IDPROG,I)
C          CALL IMPTDT('DEP-polaire', DM(DEPREE), 1, NTETA*NDDL)           
C    
        DO TETA= 1, NTETA 
C
D          CALL IMPET('au pas de temps'//IDPROG,I)
D          CALL IMPET('numero angle'//IDPROG,TETA)
C
          CALL EXTRAD(DM(DEPREE),2, TABNIV,2, TETA , DM(ADEPPA), LONG)
CD          CALL IMPTDT('DEP-polaire num=1', DM(ADEPPA), 1, NDDL)           	  
c 
C 29/08/2013 test sur la coherence des deplacemnts imposes sur le bord exterieur
D            
D 
C 29/08/2013 test sur la coherence des deplacemnts imposes sur le bord exterieur a teta fixe en polaire

D          DEBURZ=ADEPPA+ 6*NBCOL
D          DEBVT= DEBURZ +2
D          DEBW = DEBVT+2  
D          DECDEP= 6*(NBCOL+1)
C
D         DO TESTDEP= 0, 3*NBCOU-1
D           DM(ADURZ+TESTDEP)= DM(DEBURZ+TESTDEP*DECDEP)
D           DM(ADVTZ+TESTDEP)= DM(DEBVT+TESTDEP*DECDEP)
D           DM(ADWZ+TESTDEP)= DM(DEBW+TESTDEP*DECDEP)
D         END DO
C
D          CALL IMPTDT('UR dans epaisseur', DM(ADURZ), 1, 3*NBCOU)
D          CALL IMPTDT('VT dans epaisseur', DM(ADVTZ), 1, 3*NBCOU) 
D          CALL IMPTDT('W  dans epaisseur', DM(ADWZ) , 1, 3*NBCOU)    
C 
          DO L= 0, NBCOU-1
CD         CALL IMPET('POUR LA COUCHE '//idprog, L+1) 	  
C   
          DO K = 0 , NBCOL-1
CD         CALL IMPET('POUR LA COLONNE '//idprog, K+1) 	  
C  
C adresse du premier ddl U  
C         dans la numerotation ligne on a:
C 
C         l'element de la couche q de la colonne p              
C         a pour numero de premier noeud:
C                      6*(3*(NBCOL+1)(Q-1)+(P-1))pour le premier numero de DDL
C 28/08/2013 DEBUT remplace par DEBDEP
          DEBDEP =ADEPPA  +6*(3*L*(NBCOL+1)+ K) 
CD         CALL IMPET('ADRESSE DE PREMIER DDL DE LA BANDE '//idprog,
CD    &     DEBDEP-ADEPPA+1) 	     
          NBX = XINTEG-1
          IF (K .EQ. (NBCOL-1)) NBX=XINTEG 
C     
               CALL MENADM(AXYC,LONXYC)
               POSTG= AXYC
C
C Debut du calcul en tout noeud de Gauss de l'element pour tooute les coordonnees:
C Ur,Ut,W les valeurs des dep par interpolation aus noeuds de Gauss stocke(NNOERF, NCOOR=3)	
C	  
               DO COOR = 0, 2
CD        CALL IMPET( 'POUR LA COORDONNEE', COOR+1)	       
C                                                       
C - DEBUT DU TABLEAU DES VALEURS AUX NOEUDS
C   DEBUT U v OU W du noeud 1
C   DEBUR POUR LA DERIVEE DU DEPLACEMENT du noeud BAS DROITE DE L'ELEMENT TRAITE c(NUMEROTATION COUCHE)    
C 28/08/2013 DEBUT remplace DEBDEP par DEBUN
                     DEBUN=  DEBDEP+2*COOR      
                     DEBUR = DEBUN+1  
C
C     CALL IMPTET('decalage pour la notation elementaire',M(DECALE),1,5)      
C            
C
                     DM (VALDEP    )= DM(DEBUN)
                     DM (VALDEP+1  )= DM(DEBUN+M(DECALE))
                     DM (VALDEP+2  )= DM(DEBUN+M(DECALE+1))
                     DM (VALDEP+3  )= DM(DEBUN+M(DECALE+2))                       
                     DM (VALDEP+4  )= DM(DEBUN+M(DECALE+3))
                     DM (VALDEP+5  )= DM(DEBUN+M(DECALE+4))        
                     DM (VALDEP+6  )= DM(DEBUR)
                     DM (VALDEP+7  )= DM(DEBUR+M(DECALE))
                     DM (VALDEP+8  )= DM(DEBUR+M(DECALE+1))
                     DM (VALDEP+9  )= DM(DEBUR+M(DECALE+2))                       
                     DM (VALDEP+10 )= DM(DEBUR+M(DECALE+3))
                     DM (VALDEP+11 )= DM(DEBUR+M(DECALE+4))  
C		     
C            CALL IMPTDT('VALEURS DEP RANG ELEMENT',DM(VALDEP),1,12)   		C    
C		     
C   ADRESSE DE DEBUT DES VALEURS DES FONCTIONS DE BASE AUX NOEUDS DU MAILLAGE GAUSS
C   SUR L'ELEMENT DE REFERENCE = ADFBNG
C
                 FNPG=ADFBNG
C               
                 DO NUY = 0 , YINTEG             
C 
                   DO NUX = 0, NBX 
c    
                     DO DEP = 0, 11
                              DM(POSTG) = DM(POSTG)
     &                        +DM(VALDEP+DEP)*DM(FNPG)
                              FNPG       = FNPG+1   
                     END DO
                   POSTG = POSTG+1           
c         
                   END DO    
C  
                 END DO
c                
               END DO
c                 
C FIN BOUCLE EN COOR  
C  A REMETTRE SI PB                END DO
C
C  Fin du calcul en tout noeud de Gauss de l'element pour toute les coordonnees:
C Ur,Ut,W les valeurs des dep par interpolation aux noeuds de Gauss
C
C          CALL IMPTDT('VALEURS DEPL AU NOEUDS GAUSS ELT pol', 
C    &                        DM(AXYC), (YINTEG+1)*(NBX+1), 3)              
C		      
C On les tranpose pour les avoir (ncoor+3, NNOERF sans redondance entre les elements d'ou nbx) 
              LXYC1= (NBX+1)*(YINTEG+1)
              LXYC2 =  3
              CALL TRANSP( DM(AXYC)  , DM(ACXY) , LXYC1 , LXYC2 )     
C
              ADRU=ACXY
              ADRV=ACXY+1
              ADRW=ACXY+2
               DO NUY = 0 , YINTEG 
                 DO NUX = 0, NBX 
C
C 28/08/2013
C Decalage pour passer d'une numerotation d'un elt point de gauss a la numerotation globale point de Gauss 
C dans une bande (la gestion de teta est fait ensuite par une transposition le pb se pose par bande
C On rappelle qu'on est numerotation ligne 
C     pour une colonne donnee et un nuy donne (decalage en x uniquement) le decalage correspondant est
C
C      (nucol-1)*XINTEG * 3 (3 deplacement)
C 
C     pour passer a un  nuy donne le decalage est
C     (nuy-1)*(NBCOL*XINTEG+1)*3
C 
C     pour une couche donnee le decalage correspondant est
C     (NUCOU-1)*(NBCOL*XINTEG+1)*(YINTEG+1)*3
C
C pour un teta donne le decalage est
C
C      (TETA-1)*(NBCOL*XINTEG+1)*NBCOU*(YINTEG+1)*3
C      
C MODI 20/08/2013        POSI= POSTDC +(I-1)*NTETA*(NBCOL*XINTEG+1)*NBCOU*(YINTEG+1)*3	
C
C
      POSI=POSTDC+(TETA-1)*(NBCOL*XINTEG+1)*NBCOU*(YINTEG+1)*3
     &	  +  L*(NBCOL*XINTEG+1)*(YINTEG+1)*3
     &	  +  K*XINTEG * 3
     &    +  NUY*(NBCOL*XINTEG+1)*3
     &    +  3*NUX
                      DM(POSI)  = DM(ADER+TETA)*DM(ADRU)
     &                            -DM(ADTETA+TETA)*DM(ADRV)
                      POSI=POSI+1
                      DM(POSI)= DM(ADTETA+TETA)*DM(ADRU)
     &                            +DM(ADER+TETA)*DM(ADRV)
                      POSI=POSI+1
                      DM(POSI)= DM(ADRW)
		      POSI=POSI+1 
                      ADRU   =ADRU+3
                      ADRV   =ADRV+3
                      ADRW   =ADRW+3
		      INDIC =INDIC +3
C                    
                  END DO
                END DO   
C  MODIF 28/08/2013        CALL IMPET('POSI /FIN negatif? '//IDPROG, POSI-AGRXYT-LONDEG)
CD          CALL IMPET('LONGUEUR TAB DEP ELT GAUSS',(YINTEG+1)*(NBX+1)*3) 		
CD          CALL IMPTDT('VALEURS DEPL AU NOEUDS GAUSS ELT xyz', 
CD    &          DM(POSI-(YINTEG+1)*(NBX+1)*3-1), 1,(YINTEG+1)*(NBX+1)*3)  
C                            
C FIN BOUCLE EN COLONNE                 
            END DO
C FIN BOUCLE EN COUCHE      
          END DO
c fin de boucle EN TETA  
C MODIF 28/08/2013A RETIRER SI PB ET REMETTRE LIGNE 507
          END DO
C        
      DEPREE= DEPREE+NDDL*NTETA
C
C     TAB321 change un tableau double (lon1,lon2,lon3)
C     en un tableau (lon3,lon2,lon1) le tableau est en Entree-Sortie.
C 
C      NBNBAG= (NBCOL*XINTEG+1)*NBCOU*(YINTEG+1)
C      INTEGER  ADUYGA, ADUXGA, ADWGA
C
C
D      DO TETA = 1, NTETA
D          CALL IMPET('pour teta =', TETA)
D          DADUX=POSTDC+ (TETA-1)*3*NBNBAG+3*NBCOL*XINTEG
D          DADUY= DADUX+1
D          DADW = DADUY+1  
D          DECDEP= 3*(NBCOL*XINTEG+1)
D          LONEPG= NBCOU*(YINTEG+1)
C
D         DO TESTDEP= 0, LONEPG-1
D           DM(ADUXGA+TESTDEP)= DM(DADUX+TESTDEP*DECDEP)
D           DM(ADUYGA+TESTDEP)= DM(DADUY+TESTDEP*DECDEP)
D           DM(ADWGA+TESTDEP)= DM(DADW+TESTDEP*DECDEP)
D         END DO
C
D          CALL IMPTDT('UX GAUSS dans epaisseur', DM(ADUXGA), 1,LONEPG )
D          CALL IMPTDT('UY GAUSS dans epaisseur', DM(ADUYGA), 1, LONEPG) 
D          CALL IMPTDT('W GAUSS  dans epaisseur', DM(ADWGA) , 1, LONEPG)  
C
D       END DO  

      LONDEG=NTETA*3*NBNBAG*NPICET
C	 
         CALL TAB321 (3, NBNBAG, NTETA, DM(POSTDC))
C FIN BOUCLE EN TEMPS
      END DO 
c      
C MODIF 28/08/2013 D     CALL IMPET('TEST SUR LA LONGUEUR DU TABLEAU',POSI-AGRXYT-LONDEG)
D     CALL IMPET('TEST SUR LA LONGUEUR DU TABLEAU', INDIC-LONDEG)
C
C TAILLE DU TABLEAU DES DEPLACEMENT AUX NOEUD DU MAILLAGE REEL A TOUT PIQUET DE TEMPS 
C sans les derivees d'ou divise par 2
C
       LONDEN=LONDEN/2
D      NOMFI1 = 'dep_reel_polaire_sans_der'
D      CALL NUNFOU(IUNIC)
D      CALL CREFIC (2, NOMDIR, NOMFI1, 15, 'F', LONDEN, IUNIC)
D      WRITE (IUNIC, FMT=16, ERR=99000, IOSTAT=STATUS)
D    &      (DM(ADRXYT-1+I), I=1, LONDEN)
D16    FORMAT(E12.5)    
D      CALL FERFIC (2, IUNIC, IDPROG)     
C 
      NOMFI1 = 'dep_xyz_noeud_gauss'
      CALL NUNFOU(IUNIC)
      CALL CREFIC (2, NOMDIR, NOMFI1, 15, 'F', LONDEG, IUNIC)
      WRITE (IUNIC, FMT=15, ERR=99000, IOSTAT=STATUS)
     &      (DM(AGRXYT-1+I), I=1, LONDEG)
15    FORMAT(E12.5)
C   
      CALL FERFIC (2, IUNIC, IDPROG)
C
      CALL SOPOUB (AM2LC, ADM2LC)
C 
D     CALL RETOUD (IDPROG)
      RETURN
C 
99000 CONTINUE
      CALL IMPCT ('NOM DE ROUTINE ', IDPROG)
      CALL IMPET ('NUMERO D''UNITE', IUNIT)
      CALL IMPET ('STATUS ', STATUS)
      CALL ERROR_$PRINT (STATUS)
      WRITE (6, 99100) IUNIT, STATUS
99100 FORMAT (/,'PROBLEME ECRITURE SUR UNITE : ', I10,/,
     &         ' STATUS : ', I10,/)
C 
      CALL ERREUD (0,' BOFFFFFFFF ???')

      END  
