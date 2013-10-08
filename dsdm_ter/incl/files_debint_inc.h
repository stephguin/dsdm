cs-->   files_debint_inc.ftn
C * fichiers debug-tools et interpreteur.
C
      INTEGER NFECRA,NFIMPR,NFCLAV,NFERR,NFHIST,
     &        NFDONU,NFHELN,NFHELF,
     &        NFREPR,NFARGU,NFMOCL,
     &        NFCLA1,NFCLA2,NFTRAC,NFRAB2  ,
     &        nfecrs,nfimps,nfclas
C
      COMMON /FUNITI/ NFECRA,NFIMPR,NFCLAV,NFERR,NFHIST,
     &                NFDONU,NFHELN,NFHELF,
     &                NFREPR,NFARGU,NFMOCL,
     &                NFCLA1,NFCLA2,NFTRAC,NFRAB2,
     &                nfecrs,nfimps,nfclas
C
      CHARACTER*80    FIECRA,FIIMPR,FICLAV,FIERR,FIHIST
      CHARACTER*80    FIDONU,FIHELN,FIHELF
      CHARACTER*80    FIREPR,FIARGU,FIMOCL
      CHARACTER*80    FICLA1,FICLA2,FITRAC,FIRAB2
C
      COMMON /FNAMEI/ FIECRA,FIIMPR,FICLAV,FIERR,FIHIST,
     &                FIDONU,FIHELN,FIHELF,
     &                FIREPR,FIARGU,FIMOCL,
     &                FICLA1,FICLA2,FITRAC,FIRAB2
C
      LOGICAL         LOECRA,LOIMPR,LOCLAV,LOERR,LOHIST
      LOGICAL         LODONU,LOHELN,LOHELF
      LOGICAL         LOREPR,LOARGU,LOMOCL
      LOGICAL         LOCLA1,LOCLA2,LOTRAC,LORAB2
C
      COMMON /LOPENI/ LOECRA,LOIMPR,LOCLAV,LOERR,LOHIST,
     &                LODONU,LOHELN,LOHELF,
     &                LOREPR,LOARGU,LOMOCL,
     &                LOCLA1,LOCLA2,LOTRAC,LORAB2
C
      save /FUNITI/
      save /FNAMEI/
      save /LOPENI/
C
