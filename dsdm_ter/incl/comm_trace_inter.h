cs-->   comm_trace_inter.ftn
C
C         FICHIER : COMM-TRACE-INTER
C
C     * Variables de la trace conditionnelle.
C
C     * Pile des routines a tracer.
      INTEGER      NROTRA,NROTRM,IRT
      PARAMETER   (NROTRM=40)
      CHARACTER*6  ROTRAC(NROTRM)
C
C     L.ROTR sont les logiques des routines dont la trace a ete activee.
C     L.     sont les logiques de la routine courante.
C     L.ALL  sont affectes si on sohaite tout tracer.
C     L.LU   sont les logiques lus.
C
C     * Logiques de trace associ{s aux routines actives.
      LOGICAL      LSROTR(NROTRM),LS,LSALL
      LOGICAL      LNROTR(NROTRM),LN,LNALL
      LOGICAL      LPROTR(NROTRM),LP,LPALL
C
C
      COMMON / LTRACE / NROTRA,LSALL,LNALL,LPALL,
     &                         LS,LN,LP,LSROTR,LNROTR,LPROTR
      COMMON / CTRACE / ROTRAC
C
