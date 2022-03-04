      SUBROUTINE TESTCHAR (NCHAR, STRING)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  DETERMINES THE MINIMUM NONBLANK LENGTH OF A STRING
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      CHARACTER*(*) STRING
      CHARACTER BLANK
c
      DATA BLANK/' '/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      NMAX=LEN(STRING)
      NCHAR=0
c
      DO 10 I=1,NMAX
         ITEST=NMAX-I+1
         IF(STRING(ITEST:ITEST).NE.BLANK) THEN
           NCHAR=ITEST
           RETURN
         ENDIF
10    CONTINUE
c
      RETURN
      END

