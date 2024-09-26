MODULE XASDET_MODULE
CHARACTER(256) :: FILE
LOGICAL(4) :: TINIT=.FALSE.
LOGICAL(4) :: TWAKE=.FALSE.
END MODULE XASDET_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE XASDET$SETL4(ID,VAL)
!      *************************************************************************
       USE XASDET_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       LOGICAL(4),INTENT(IN) :: VAL
!      *************************************************************************
       IF(ID.EQ.'INIT') THEN
         TINIT=VAL
       ELSE IF(ID.EQ.'WAKE') THEN
         TWAKE=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('XASDET$SETL4')
       END IF
       RETURN
       END SUBROUTINE XASDET$SETL4
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE XASDET$GETL4(ID,VAL)
!      *************************************************************************
       USE XASDET_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       LOGICAL(4),INTENT(OUT) :: VAL
!      *************************************************************************
       IF(ID.EQ.'INIT') THEN
         VAL=TINIT
       ELSE IF(ID.EQ.'WAKE') THEN
         VAL=TWAKE
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('XASDET$GETL4')
       END IF
       RETURN
       END SUBROUTINE XASDET$GETL4
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE XASDET$SETCH(ID,VAL)
!      *************************************************************************
       USE XASDET_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       CHARACTER(*),INTENT(IN) :: VAL
!      *************************************************************************
       IF(ID.EQ.'FILE') THEN
         FILE=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('XASDET$SETCH')
       END IF
       RETURN
       END SUBROUTINE XASDET$SETCH
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE XASDET$GETCH(ID,VAL)
!      *************************************************************************
       USE XASDET_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       CHARACTER(*),INTENT(OUT) :: VAL
!      *************************************************************************
       IF(ID.EQ.'FILE') THEN
         VAL=FILE
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('XASDET$GETCH')
       END IF
       RETURN
       END SUBROUTINE XASDET$GETCH
! !
! !      ..1.........2.........3.........4.........5.........6.........7.........8
!        SUBROUTINE XASDET$WRITE
! !      *************************************************************************
!        USE XASDET_MODULE
!        IMPLICIT NONE
!        INTEGER(4) :: NFIL
! !      *************************************************************************
!                           CALL TRACE$PUSH('XASDET$WRITE')
!        IF(.NOT.TINIT) RETURN
!        IF(.NOT.TWAKE) RETURN

!        CALL FILEHANDLER$SETFILE('XASDET',.FALSE.,FILE)
!        CALL FILEHANDLER$SETSPECIFICATION('XASDET','STATUS','UNKNOWN')
!        CALL FILEHANDLER$SETSPECIFICATION('XASDET','POSITION','REWIND')
!        CALL FILEHANDLER$SETSPECIFICATION('XASDET','ACTION','READWRITE')
!        CALL FILEHANDLER$SETSPECIFICATION('XASDET','FORM','UNFORMATTED')
!        CALL FILEHANDLER$UNIT('XASDET',NFIL)

!        CALL XASDET_PHI(NFIL)

!                           CALL TRACE$POP
!        END SUBROUTINE XASDET$WRITE
! !
!        SUBROUTINE XASDET_PHI(NFIL)  ! MARK: XASDET_PSPHI
! !      *************************************************************************
!        USE SETUP_MODULE, ONLY: NSP,THIS
!        IMPLICIT NONE
!        INTEGER(4), INTENT(IN) :: NFIL
!        INTEGER(4) :: ISP
!        INTEGER(4) :: NR
!        REAL(8) :: DEX
!        REAL(8) :: R1
!        CHARACTER(8) :: GRIDTYPE
! !      *************************************************************************
!                           CALL TRACE$PUSH('XASDET_PHI')
! !      NUMBER OF SETUPS/SPECIES
!        WRITE(NFIL)NSP
! !      LOOP THROUGH SETUPS
!        DO ISP=1,NSP
!          CALL SETUP$ISELECT(ISP)
! !        GET GRID DETAILS
!          CALL RADIAL$GETI4(THIS%GID,'NR',NR)
!          CALL RADIAL$GETR8(THIS%GID,'DEX',DEX)
!          CALL RADIAL$GETR8(THIS%GID,'R1',R1)
!          CALL RADIAL$GETCH(THIS%GID,'TYPE',GRIDTYPE)
! !        WRITE GRID DETAILS
!          WRITE(NFIL)GRIDTYPE,NR,DEX,R1,THIS%LNX
! !        WRITEAEPHI AND PSPHI
!          WRITE(NFIL)THIS%AEPHI
!          WRITE(NFIL)THIS%PSPHI
!          CALL SETUP$UNSELECT
!        END DO
!                           CALL TRACE$POP
!        END SUBROUTINE XASDET_PHI

!        SUBROUTINE XASDET_PDOS(NFIL)  ! MARK: XASDET_PDOS
! !      *************************************************************************
!        USE PDOS_MODULE
!        IMPLICIT NONE
!        INTEGER(4), INTENT(IN) :: NFIL
!        INTEGER(4) :: ISP,IKPT
!        INTEGER(4) :: LNX1
!        INTEGER(4) :: ILOGICAL
!        INTEGER(4) :: IB1
! !      *************************************************************************
!                           CALL TRACE$PUSH('XASDET_PDOS')
!        WRITE(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX
!        WRITE(NFIL)LNX(:),LOX(:,:),ISPECIES(:)
!        ILOGICAL=1
!        IF(.NOT.TINV)ILOGICAL=0
!        WRITE(NFIL)NKDIV(:),ISHIFT(:),RNTOT,NEL,ILOGICAL
!        ILOGICAL=1
!        IF(.NOT.TSHIFT)ILOGICAL=0
!        WRITE(NFIL)ILOGICAL

!        WRITE(NFIL)RBAS,R,ATOMID

!        WRITE(NFIL)XK,NB,WKPT
!        DO IB1=1,NB
!                           CALL TRACE$POP
!        RETURN
!        END SUBROUTINE XASDET_PDOS


      
