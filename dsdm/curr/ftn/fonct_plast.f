C     On envoie comme arguments :
C 
C     E ...... S12CP             valeur precedente de s12 chapeau
C     E ...... S22CP             valeur precedente de s22 chapeau
C     E ...... AS12CP            accroissement de s12 chapeau
C     E ...... AS22CP            accroissement de s22 chapeau
C     E ...... CD                ancienne valeur de d
C     E ...... DP                ancienne valeur de d'
C 
C     Le tableau carnli stocke sequentiellement comme suit :
C 
C                YO                valeur Yinitiale
C                YC                valeur  de Ycritique
C                B                 valeur du couplage d'endommagement
C                E220              valeur initiale du module E22
C                G120              valeur initiale du module 2*G12
C                Ro                valeur du seuil de plasticite
C                BETA              valeurs telle que p =  beta * R
C                ALPHA
C                A2                valeur du couplage de plasticite
C                S11INF            valeur mini de la contrainte 11
C                S11SUP            valeur maxi de la contrainte 11
C                S22SUP            valeur maxi de la contrainte 22
C 
C     Les contraintes sont les contraintes chapau dans la base d'orthotropie
C 
      DOUBLE PRECISION FUNCTION AINTG (S12CP, S22P, AS12CP, AS22P,
     &                                 D, DP, CARNLI)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION  S12CP, S22P, AS12CP, AS22P, D, DP
      DOUBLE PRECISION  CARNLI(12)
C 
C -----------------------------------------------------------------------
C  
C     DECLARATION DES PARAMETRES LOCAUX
C     """"""""""""""""""""""""""""""""""
C 
      DOUBLE PRECISION  TERM, D2, DS122, DS222
C 
      TERM   = S12CP+AS12CP
      D2     = 1.D0-DP
      DS122  = (TERM*TERM - S12CP*S12CP)/( D2*D2)
C 
      TERM   = S22P+AS22P
      D2     = 1.D0-DP
      DS222  = 2.D0*CARNLI(9)*(TERM*TERM - S22P*S22P)/( D2*D2)
C 
      AINTG = DS122 + DS222
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Les contraintes sont les contraintes chapau dans la base d'orthotropie
C 
      DOUBLE PRECISION FUNCTION TOTG (GINTPR, ACGINT)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION GINTPR , ACGINT
C 
      TOTG = DDIM (GINTPR+ACGINT, 0.D0)
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Calcul du seuil de plasticite :
C 
C     On envoie comme arguments:
C     E ...... RP      valeur du seuil plastique precedent locale
C     E ...... GAC     valeur actuelle de la fonction TOTG
C 
C     Et on recupere :
C 
C     E ...... ACPLAS  Logique indiquant s'il y a plastification
C     E ...... DELTR   acroissement du seuil
C 
      DOUBLE PRECISION FUNCTION RSEUIL (RP, GAC, ACPLAS, DELTR)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION RP, GAC, DELTR
      LOGICAL          ACPLAS
C 
      DELTR = GAC - RP

      IF (DELTR .GT. 0.D0) THEN
C 
        RSEUIL = GAC
        ACPLAS = .TRUE.
C 
      ELSE
C 
        RSEUIL = RP
        ACPLAS = .FALSE.
        DELTR  = 0.D0
C 
      END IF
C 
      RETURN
      END
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C 
C     Calcul du seuil de plasticite .
C 
C     On envoie comme arguments :
C 
C     E ...... PANC    ancienne valeur de p deformation plastique cumulee
C     E ...... ALPHA   coeff de la loi d'ecrouissage
C                                              alpha
C     E ...... BETA           R = Ro + beta * p
C     E....... DELTR   acroissement du seuil
C 
C     Et on recupere :
C 
C     S ...... DELTP   acroissement de la deformation plastique cumulee
C 
      DOUBLE PRECISION FUNCTION PNOUV (PANC, BETA, ALPHA, DELTR, DELTP)
C 
C -----------------------------------------------------------------------
C 
C     DECLARATION DES PARAMETRES GLOBAUX
C     """"""""""""""""""""""""""""""""""
C   
      DOUBLE PRECISION PANC, DELTR, BETA, ALPHA, DELTP, PUIS
C 
      PUIS = 1.D0/ALPHA
C 
      DELTP = ((DELTR/BETA + (PANC**ALPHA)) ** PUIS)
      PNOUV = PANC + DELTP
C 
      RETURN
      END
