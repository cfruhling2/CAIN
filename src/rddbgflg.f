	SUBROUTINE RDDBGFLG(LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
      IMPLICIT NONE
      INTEGER LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
	INCLUDE 'include/alloccm.h'
      INCLUDE 'include/ctrlcm.h'
C      INCLUDE 'include/readcm.h'
C
      INTEGER MOP
      PARAMETER (MOP=1)
      CHARACTER*12 OP(MOP)/'ID'/
      INTEGER NFF(MOP)/1/
	INTEGER ID1
      INTEGER NC,K,J,NF,I,N
C
      IRTN=0
      IF(LN(1,2).EQ.0) RETURN
      CALL CMDBLK('DEBUGFLAG',MOP,OP,0,NFF,MBL,NBL,KOP,KBL,REL,
     %    LNKW,LNBL,LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.1) GOTO 990
      IF(IRTN.EQ.1.OR.NBL.LE.0) RETURN
      J=1
      DO 120 I=1,MOP
        ID(I)=J
        J=J+MAX(0,NFF(I))
 120  CONTINUE
      DO 180 I=1,J-1
        PAR(I)=UNDEF
 180  CONTINUE
	NGSTRRD=0
	ID1=0

      DO 300 J=1,NBL
        I=KBL(J)
        IF(NFF(I).GE.1) THEN
          CALL BLKREC(LNBL(1,1,J),LINE,NCHAR,TEXT,NC,IRTN,MSGFL)
          IF(IRTN.NE.0) GOTO 990
          DO 220 K=1,NFF(I)
	      FFF(K)=UNDEF
 220      CONTINUE
          CALL BLKRHS(TEXT(1:NC),FFF,MFFF,NF,GSTRRD,NGSTRRD,IRTN)
          IF(IRTN.NE.0) GOTO 990
          IF(NF.GT.NFF(I)) GOTO 940
          IF(NF.GE.1) THEN
            DO 240 K=1,NF
              IF(FFF(K).NE.UNDEF) PAR(ID(I)+K-1)=FFF(K)
 240        CONTINUE
          ENDIF
	  ELSE
        ENDIF
 300  CONTINUE
	DO 400 I=1,MOP
	  IF(OP(I).EQ.'ID') THEN
	    IF(PAR(ID(I)).NE.UNDEF) ID1=NINT(PAR(ID(I))%X)
	  ENDIF
 400  CONTINUE
	IF(ID1.GE.1.AND.ID1.LE.8) IDBGFLG(ID1)=1
	RETURN
 940  IRTN=1004
      WRITE(MSGFL,945) OP(J)
 945  FORMAT(' (SUBR.RDDBGFLG) Too many numbers for ',
     %  'operand "',A,'".')
      GOTO 995
 990  IRTN=1009
      GOTO 995
 995  RETURN
      END
