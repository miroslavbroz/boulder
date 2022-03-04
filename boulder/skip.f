c =====================================================================
      SUBROUTINE SKIP (IUNIT)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ROUTINE TO SKIP LINES BEGINNING WITH THE STRING /*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      CHARACTER  STRING*80,LEADER*2
c
      DATA LEADER/'/*'/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   10 READ(IUNIT,1000,END=20) STRING
      call testchar (NCHAR, STRING)
      II = NCHAR
      IF(II.EQ.0) GOTO 10
c
      IF (INDEX(STRING,LEADER).EQ.1) THEN
         GOTO 10
      ELSE
         BACKSPACE (IUNIT)
      ENDIF
c
 1000 FORMAT(A80)
c
   20 RETURN
      END
C
