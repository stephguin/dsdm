c 
c     ------------------------------
c     incl_graphics
c     ------------------------------
c
c
c                       * LONGUEUR D'UN TEXTE * 
c                       * UTILISE PAR PAD_$CREATE_WINDOW *
      integer*2 name_length 
c
c                       * TYPE DE PAD *
c                       * UTILISE PAR PAD_$CREATE_WINDOW *
      integer*2 pad_type 
c
c                       * NUMERO D'UNITE *
c                       * UTILISE PAR PAD_$CREATE_WINDOW *
      integer*2 node_unit
c
c                       * POSITION ET DIMENSIONS DE LA FENETRE *
c                       * UTILISE PAR PAD_$CREATE_WINDOW *
      integer*2 window(4)
c                       * NUMERO DE LA FENETRE *
c                       * RENVOYE PAR PAD_$INQ_WINDOWS *
c                       * UTILISE PAR PAD_$SET_AUTO_CLOSE *
      integer*2 window_no
c
c                       * NUMERO DU STREAM *
c                       * RENVOYE PAR PAD_$CREATE_WINDOW *  
c                       * UTILISE PAR GPR_$INIT          *
      integer*2 stream_id 
c
c                       * STATUS *   
c                       * RENVOYE PAR TOUS LES APPELS *
      integer status
c
c                       * MODE D'OPERATION *
c                       * UTILISE PAR GPR_$INIT *
      integer*2 op_mode 
c                        
c                       * TAILLE DE LA BITMAP *
c                       * UTILISE DANS GPR_$INIT *
      integer*2 init_bitmap_size(2)
c
c                       * NUMERO DU PLAN DE DESSIN *
c                       * UTILISE PAR GPR_$INIT *
      integer*2 hi_plane_id
c                       * DESCRIPTEUR DE LA BITMAP *
c                       * RENVOYE PAR GPR_$INIT *
      integer *4 init_bitmap_desc 
c
c                       * FAUT-IL SUPPRIMER LE FICHIER/FENETRE *
c                       * UTILISE PAR GPR_$TERMINATE * 
      logical delete_display
c
c                       * L'ECRAN EST NOIR *
c                       * RENVOYE PAR GPR_$ACQUIRE_DISPLAY 
      logical unobsc
c
c                       * TYPE D'INTERRUPTION PERMISE *
c                       * UTILISE PAR GPR_$ENABLE_INPUT *
      integer*2 event_type
c
c                       * CARACTERES ADMIS PAR RETICU *
c                       * UTILISE PAR GPR_$ENABLE_INPUT *
      integer*4 key_set(8)
      save key_set
c
c                       * CURSEUR ACTIF *
c                       * UTILISE PAR GPR_$SET_CURSOR_ACTIVE *
      logical cursor
c
c                       * CARACTERE RECU * 
c                       * UTILISE PAR GPR_$EVENT_WAIT *
      character*1 event_data
c
c                       * POSITION DU CURSEUR *
c                       * UTILISE PAR GPR_$EVENT_WAIT *
      integer*2 curs_position(2)
c
c
c
c                       * DIMENSIONS DE LA FENETRE EN COMMON *
      integer*2 pyhau 
      real hau,lar
c
c
c                       * FONCTIONS DE CONVERSION CM/PIXELS *
      integer*2 cmpi,cmpiy
      real picm
c
c     * ----------------------------------------------------
c     declaration external fait chier dans le block-data.
c
ccccc      external cmpi,cmpiy,picm 
c     * ----------------------------------------------------
c
c
c                       * ON EST EN MODE DIRECT OU EN MODE FRAME *
      logical direct_mode 
c
c                       * EMULATION NORSK-DATA *
      logical emul_norsk
c
c
C
C                       * TRAITEMENT DE LA COULEUR *
C *
C *   type de terminal (couleur ou non) ,affecte par trunic
      LOGICAL TERCOL
C *                   
C *   anciene couleur et couleur courante
      INTEGER ILASCO , ICOLOR
C
      COMMON /COULEU / TERCOL , ILASCO , ICOLOR
      SAVE   /COULEU /
C
C                       * MEMORISATION STREAM_ID *
      INTEGER SSTRMD
C
C                       * FONTE *
      INTEGER LONFNT
      CHARACTER*256 FONTE
      INTEGER*2 HORIZONTAL_SPACING,CHARACTER_WIDTH
      REAL LARCAR
      INTEGER PLARCAR
      common /cpyhau/pyhau,hau,lar,direct_mode,emul_norsk,stream_id 
     &,SSTRMD,LONFNT,FONTE,LARCAR,PLARCAR
C
C
C     ----------------------------------------------
C     - status des appels gpr$  , logique appel ok
C
      LOGICAL   LOKGPR
C
      COMMON / CSTATU / STATUS  ,  LOKGPR
C
c
      COMMON / GRAPHG / init_bitmap_desc 
c
C
