! CONTROL FILE
!    ALLOW FOR SHIFT OF AMPLITUDES TO A DIFFERENT EQUIVALENT ATOM
! PARALLELIZATION NECESSARY?
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE XCCNTL_MODULE  ! MARK: XCNTL_MODULE
!     **************************************************************************
!     ** XCNTL MODULE FOR PAW XRAY TOOL                                       **
!     **************************************************************************
        USE LINKEDLIST_MODULE, ONLY: LL_TYPE
        TYPE(LL_TYPE) :: LL_CNTL
        LOGICAL(4) :: TOVERLAP=.FALSE. ! OVERLAP INPUT FILE PRESENT
        CHARACTER(256) :: OVERLAPFILE ! FILENAME FOR OVERLAP INPUT FILE
        SAVE
      END MODULE XCCNTL_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE RIXSAMPL_MODULE

      TYPE ONEKPT_TYPE
      INTEGER(4) :: NB1
      INTEGER(4) :: NOCC
      COMPLEX(8), ALLOCATABLE :: XY(:,:) ! (NB1-NOCC,NOCC)
      REAL(8), ALLOCATABLE :: EIG(:) ! (NB1)
      END TYPE ONEKPT_TYPE

      TYPE TWOKPT_TYPE
      INTEGER(4) :: NB1
      INTEGER(4) :: NOCC
      COMPLEX(8), ALLOCATABLE :: X(:) ! (NB1-NOCC)
      COMPLEX(8), ALLOCATABLE :: Y(:) ! (NOCC)
      REAL(8), ALLOCATABLE :: EIG(:) ! (NB1)
      END TYPE TWOKPT_TYPE

      TYPE RESULT_TYPE
      INTEGER(4) :: NB1
      INTEGER(4) :: NOCC
      COMPLEX(8), ALLOCATABLE :: XY(:,:) ! (NB1-NOCC,NOCC)
      REAL(8), ALLOCATABLE :: ABSXY(:,:) ! (NB1-NOCC,NOCC)
      REAL(8), ALLOCATABLE :: EIG(:)
      END TYPE RESULT_TYPE

      TYPE SPEC_TYPE
      CHARACTER(256) :: FILENAME ! FILENAME FOR THE SPEC
      REAL(8) :: RHOLE(3) ! POSITION OF THE COREHOLE
      TYPE(ONEKPT_TYPE), ALLOCATABLE :: ONEKPT(:,:) ! (NKPTTOT,NSPIN)
      TYPE(TWOKPT_TYPE), ALLOCATABLE :: TWOKPT(:,:) ! (NKPTTOT,NSPIN)
      END TYPE SPEC_TYPE

      ! THE SAME FOR EVERYTHING
      INTEGER(4) :: NKPT
      INTEGER(4) :: NKPTTOT
      INTEGER(4) :: NSPIN
      INTEGER(4) :: NKDIV(3)
      INTEGER(4), ALLOCATABLE :: TOIRRKPT(:) ! (NKPTOT)
      LOGICAL(4), ALLOCATABLE :: TIRRKPT(:) ! (NKPTTOT)
      CHARACTER(6) :: FLAG
      REAL(8) :: EMIN
      REAL(8) :: EMAX
      REAL(8) :: ELIGHT
      REAL(8) :: NORMAL(3)
      REAL(8) :: KDIRI(3)
      COMPLEX(8) :: POLI(2)
      COMPLEX(8) :: POLXYZI(3)
      REAL(8) :: KI(3)
      REAL(8) :: XKI(3)
      REAL(8) :: KDIRO(3)
      COMPLEX(8) :: POLO(2)
      COMPLEX(8) :: POLXYZO(3)
      REAL(8) :: KO(3)
      REAL(8) :: XKO(3)
      REAL(8) :: Q(3)
      REAL(8) :: XQ(3)
      REAL(8) :: QAPPROX(3)
      REAL(8) :: XQAPPROX(3)
      REAL(8) :: QERROR
      LOGICAL(4) :: TKPTSHIFT


      LOGICAL(4) :: TFIRST=.TRUE.

      INTEGER(4) :: NSPEC=-1 ! NUMBER OF FCH CALCULATIONS

      TYPE(SPEC_TYPE), ALLOCATABLE, TARGET :: SPEC(:) ! (NSPEC)

      TYPE(SPEC_TYPE), POINTER :: THIS ! POINTER TO THE CURRENT SPEC
      LOGICAL(4) :: SELECTED=.FALSE. ! SELECTED SPEC FOR THE CURRENT TASK

      TYPE(RESULT_TYPE), ALLOCATABLE :: RES(:,:) !(NKPTTOT,NSPIN)
      END MODULE RIXSAMPL_MODULE

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM PAW_RIXSCOMBINE ! MARK: PAW_RIXSCOMBINE
      USE CLOCK_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NFIL
      CHARACTER(32) :: DATIME
      CALL MPE$INIT
                          CALL TRACE$PUSH('PAW_RIXSCOMBINE')
                          CALL TIMING$START
      ! INITIALIZE FILES
      CALL INITIALIZEFILEHANDLER
      CALL FILEHANDLER$UNIT('PROT',NFIL)
      CALL CPPAW_WRITEVERSION(NFIL)
      CALL CLOCK$NOW(DATIME)
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='(80("="),T15,"  PROGRAM STARTED ",A32,"  ")')DATIME
      WRITE(NFIL,FMT='(80("="))')

      CALL FILEHANDLER$UNIT('XCCNTL',NFIL)
      CALL XCCNTL$READ(NFIL)
      CALL XCCNTL$AMPLITUDE

      CALL INPUT

      CALL FILEHANDLER$UNIT('PROT',NFIL)

      CALL RIXSAMPL$REPORT(NFIL)

      CALL RIXSAMPL$CALCULATE

                          CALL TIMING$PRINT('~',NFIL)
      CALL MPE$CLOCKREPORT(NFIL)
      CALL CLOCK$NOW(DATIME)
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='(80("="),T15,"  PROGRAM FINISHED ",A32,"  ")')DATIME
      WRITE(NFIL,FMT='(80("="))')

                          CALL TRACE$POP
      CALL ERROR$NORMALSTOP
      END PROGRAM PAW_RIXSCOMBINE

      ! ONEKPT
      ! ATOMS NEED TO BE COMBINED FIRST
      ! STORE THE AMPLITUDES PARALLELIZED OVER KPT AND SPIN
      ! PERFORM SUM OVER ATOMS WITHIN EACH TASK
      ! FORM ABSOLUTE SQUARE OF THE AMPLITUDE WITHIN EACH TASK
      ! SUM OVER KPT, SPIN, IOCC, IEMP WITHIN EACH TASK
      ! COMBINE OVER TASKS

      ! TWO KPT REQUIRES MULTIPLYING THE COMBINATIONS OF IOCC AND IEMP FIRST
      ! THEN THE ABOVE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INPUT  ! MARK: INPUT
!     **************************************************************************
!     ** INPUT: READS THE INPUT FILE AND INITIALIZES RIXSAMPL MODULE          **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NFIL
      INTEGER(4) :: NSPEC
      INTEGER(4) :: ISPEC
      CHARACTER(256) :: FILENAME
                          CALL TRACE$PUSH('INPUT')

      CALL RIXSAMPL$GETI4('NSPEC',NSPEC)

      CALL FILEHANDLER$SETFILE('RIXSAMPL',.TRUE.,-'.RIXS')
      CALL FILEHANDLER$SETSPECIFICATION('RIXSAMPL','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('RIXSAMPL','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('RIXSAMPL','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('RIXSAMPL','FORM','UNFORMATTED')
      

      DO ISPEC=1,NSPEC
        CALL RIXSAMPL$ISELECT(ISPEC)

        CALL RIXSAMPL$GETCH('FILENAME',FILENAME)

        CALL FILEHANDLER$SETFILE('RIXSAMPL',.FALSE.,TRIM(ADJUSTL(FILENAME)))
        CALL FILEHANDLER$UNIT('RIXSAMPL',NFIL)

        CALL INPUTHEADER(NFIL)

        CALL INPUTDATA(NFIL)

        CALL RIXSAMPL$UNSELECT
      ENDDO
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE INPUT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INPUTHEADER(NFIL)  ! MARK: INPUTHEADER
!     **************************************************************************
!     ** INPUTHEADER: READS THE HEADER OF THE INPUT FILE                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NFIL
      INTEGER(4) :: NKPT
      INTEGER(4) :: NKPTTOT
      INTEGER(4) :: NSPIN
      INTEGER(4) :: NKDIV(3)
      REAL(8) :: RHOLE(3)
      INTEGER(4), ALLOCATABLE :: IARR(:)
      LOGICAL(4), ALLOCATABLE :: TIRRKPT(:)
      INTEGER(4), ALLOCATABLE :: TOIRRKPT(:) ! (NKPTTOT)
      CHARACTER(6) :: FLAG
      INTEGER(4) :: IKPT
      INTEGER(4) :: IKPTTOT
      REAL(8) :: EMIN,EMAX,ELIGHT
      REAL(8) :: NORMAL(3)
      REAL(8) :: KDIRI(3)
      COMPLEX(8) :: POLI(2)
      COMPLEX(8) :: POLXYZI(3)
      REAL(8) :: KI(3)
      REAL(8) :: XKI(3)
      REAL(8) :: KDIRO(3)
      COMPLEX(8) :: POLO(2)
      COMPLEX(8) :: POLXYZO(3)
      REAL(8) :: KO(3)
      REAL(8) :: XKO(3)
      REAL(8) :: Q(3)
      REAL(8) :: XQ(3)
      REAL(8) :: QAPPROX(3)
      REAL(8) :: XQAPPROX(3)
      REAL(8) :: QERROR
      INTEGER(4) :: ILOGICAL
      LOGICAL(4) :: TKPTSHIFT
                          CALL TRACE$PUSH('INPUTHEADER')
      ! NKPT,NKPTTOT,NSPIN,NKDIV,RHOLE
      READ(NFIL)NKPT,NKPTTOT,NSPIN,NKDIV,RHOLE
      CALL RIXSAMPL$SETI4('NKPT',NKPT)
      CALL RIXSAMPL$SETI4('NKPTTOT',NKPTTOT)
      CALL RIXSAMPL$SETI4('NSPIN',NSPIN)
      CALL RIXSAMPL$SETI4A('NKDIV',3,NKDIV)
      CALL RIXSAMPL$SETR8A('RHOLE',3,RHOLE)
      ALLOCATE(IARR(NKPTTOT))
      ALLOCATE(TIRRKPT(NKPTTOT))
      ALLOCATE(TOIRRKPT(NKPTTOT))
      ! TIRRKPT(IARR)
      READ(NFIL)TOIRRKPT,IARR
      DO IKPTTOT=1,NKPTTOT
        TIRRKPT(IKPTTOT)=IARR(IKPTTOT).EQ.1
      END DO
      DEALLOCATE(IARR)
      CALL RIXSAMPL$SETL4A('TIRRKPT',NKPTTOT,TIRRKPT)
      CALL RIXSAMPL$SETI4A('TOIRRKPT',NKPTTOT,TOIRRKPT)
      DEALLOCATE(TIRRKPT)
      DEALLOCATE(TOIRRKPT)

      READ(NFIL)FLAG,EMIN,EMAX,ELIGHT,NORMAL,KDIRI,POLI,POLXYZI,KI,XKI,KDIRO, &
     &          POLO,POLXYZO,KO,XKO,Q,XQ,QAPPROX,XQAPPROX,QERROR,ILOGICAL
      TKPTSHIFT=ILOGICAL.EQ.1
      CALL RIXSAMPL$SETCH('FLAG',FLAG)
      CALL RIXSAMPL$SETR8('EMIN',EMIN)
      CALL RIXSAMPL$SETR8('EMAX',EMAX)
      CALL RIXSAMPL$SETR8('ELIGHT',ELIGHT)
      CALL RIXSAMPL$SETR8A('NORMAL',3,NORMAL)
      CALL RIXSAMPL$SETR8A('KDIRI',3,KDIRI)
      CALL RIXSAMPL$SETC8A('POLI',2,POLI)
      CALL RIXSAMPL$SETC8A('POLXYZI',3,POLXYZI)
      CALL RIXSAMPL$SETR8A('KI',3,KI)
      CALL RIXSAMPL$SETR8A('XKI',3,XKI)
      CALL RIXSAMPL$SETR8A('KDIRO',3,KDIRO)
      CALL RIXSAMPL$SETC8A('POLO',2,POLO)
      CALL RIXSAMPL$SETC8A('POLXYZO',3,POLXYZO)
      CALL RIXSAMPL$SETR8A('KO',3,KO)
      CALL RIXSAMPL$SETR8A('XKO',3,XKO)
      CALL RIXSAMPL$SETR8A('Q',3,Q)
      CALL RIXSAMPL$SETR8A('XQ',3,XQ)
      CALL RIXSAMPL$SETR8A('QAPPROX',3,QAPPROX)
      CALL RIXSAMPL$SETR8A('XQAPPROX',3,XQAPPROX)
      CALL RIXSAMPL$SETR8('QERROR',QERROR)
      CALL RIXSAMPL$SETL4('TKPTSHIFT',TKPTSHIFT)

      CALL RIXSAMPL$SETL4('TFIRST',.FALSE.)
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE INPUTHEADER
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INPUTDATA(NFIL)  ! MARK: INPUTDATA
!     **************************************************************************
!     ** READS THE INPUT FILE DATA                                            **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NFIL
      INTEGER(4) :: NKPTTOT
      INTEGER(4) :: IKPT
      INTEGER(4) :: NKPT,NSPIN
      INTEGER(4) :: ISPIN
      INTEGER(4) :: NB1,NOCC
      INTEGER(4) :: RKPT,RSPIN ! READ KPT AND SPIN
      CHARACTER(6) :: FLAG
      REAL(8), ALLOCATABLE :: EIG(:) ! (NB1)
      COMPLEX(8), ALLOCATABLE :: XY(:,:) ! (NB1-NOCC,NOCC)
      COMPLEX(8), ALLOCATABLE :: X(:) ! (NB1-NOCC)
      COMPLEX(8), ALLOCATABLE :: Y(:) ! (NOCC)
                          CALL TRACE$PUSH('INPUTDATA')
      CALL RIXSAMPL$INIT
      CALL RIXSAMPL$GETI4('NKPT',NKPT)
      CALL RIXSAMPL$GETI4('NSPIN',NSPIN)
      CALL RIXSAMPL$GETCH('FLAG',FLAG)

      READ(NFIL)NKPTTOT
      CALL RIXSAMPL$SETI4('NKPTTOT',NKPTTOT)

      
      DO IKPT=1,NKPTTOT
        DO ISPIN=1,NSPIN
          ! IKPTTOT,ISPIN,NB1,NOCC
          READ(NFIL)RKPT,RSPIN,NB1,NOCC
          IF(RKPT.GT.NKPTTOT.OR.RKPT.LT.1) THEN
            CALL ERROR$MSG('READ KPT OUT OF RANGE')
            CALL ERROR$I4VAL('RKPT',RKPT)
            CALL ERROR$I4VAL('NKPTTOT',NKPTTOT)
            CALL ERROR$STOP('INPUTDATA')
          END IF
          IF(RSPIN.GT.NSPIN.OR.RSPIN.LT.1) THEN
            CALL ERROR$MSG('READ SPIN OUT OF RANGE')
            CALL ERROR$I4VAL('RSPIN',RSPIN)
            CALL ERROR$I4VAL('NSPIN',NSPIN)
            CALL ERROR$STOP('INPUTDATA')
          END IF
          ALLOCATE(EIG(NB1))
          ! EIG(NB1)
          READ(NFIL)EIG
          ! AMPLITUDES FOR TWOKPT
          IF(FLAG.EQ.+'TWOKPT') THEN
            ALLOCATE(X(NB1-NOCC))
            ALLOCATE(Y(NOCC))
            ! X(NB1-NOCC)
            READ(NFIL) X
            ! Y(NOCC)
            READ(NFIL) Y
            ! STORE DATA IN THE SPEC
            CALL RIXSAMPL$SETTWOKPT(RKPT,RSPIN,NB1,NOCC,EIG,X,Y)
            DEALLOCATE(X)
            DEALLOCATE(Y)
          ELSE IF(FLAG.EQ.+'ONEKPT') THEN
            ALLOCATE(XY(NB1-NOCC,NOCC))
            ! XY(NB1-NOCC,NOCC)
            READ(NFIL) XY
            ! STORE DATA IN THE SPEC
            CALL RIXSAMPL$SETONEKPT(RKPT,RSPIN,NB1,NOCC,EIG,XY)
            DEALLOCATE(XY)
          END IF
          DEALLOCATE(EIG)
        ENDDO ! END ISPIN
      ENDDO ! END IKPTTOT
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE INPUTDATA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$INIT  ! MARK: RIXSAMPL$INIT
!     **************************************************************************
!     ** INITIALIZE RIXSAMPL STORAGE                                          **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: SPEC,NSPEC,TFIRST,FLAG,NKPTTOT,NSPIN
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4) :: ISPEC
      ! DO NOT INITIALIZE BEFORE NSPEC,NKPTTOT,NSPIN ARE SET
      IF(TFIRST) RETURN
      IF(.NOT.ALLOCATED(SPEC)) THEN
        CALL ERROR$MSG('RIXSAMPL SPEC NOT ALLOCATED')
        CALL ERROR$STOP('RIXSAMPL$INIT')
      END IF
      DO ISPEC=1,NSPEC
        IF(FLAG.EQ.+'ONEKPT') THEN
          IF(.NOT.ALLOCATED(SPEC(ISPEC)%ONEKPT)) THEN
            ALLOCATE(SPEC(ISPEC)%ONEKPT(NKPTTOT,NSPIN))
          END IF
        ELSE IF(FLAG.EQ.+'TWOKPT') THEN
          IF(.NOT.ALLOCATED(SPEC(ISPEC)%TWOKPT)) THEN
            ALLOCATE(SPEC(ISPEC)%TWOKPT(NKPTTOT,NSPIN))
          END IF
        ELSE
          CALL ERROR$MSG('RIXSAMPL FLAG NOT RECOGNIZED')
          CALL ERROR$CHVAL('FLAG: ',FLAG)
          CALL ERROR$STOP('RIXSAMPL$INIT')
        END IF        
      ENDDO
      RETURN
      END SUBROUTINE RIXSAMPL$INIT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$REPORT(NFIL)  ! MARK: RIXSAMPL$GETI4
!     **************************************************************************
!     ** REPORT DATA                                                          **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: NKPT,NKPTTOT,NSPIN,NKDIV,FLAG,EMIN,EMAX, &
     &      ELIGHT,NORMAL,KDIRI,POLI,POLXYZI,KI,XKI,KDIRO,POLO,POLXYZO,KO,XKO, &
     &      Q,XQ,QAPPROX,XQAPPROX,QERROR,TKPTSHIFT,TFIRST
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NFIL
      REAL(8) :: EV
      REAL(8) :: ANGSTROM
                          CALL TRACE$PUSH('RIXSAMPL$REPORT')
      IF(TFIRST) THEN
        CALL ERROR$MSG('RIXSAMPL NOT INITIALIZED')
        CALL ERROR$STOP('RIXSAMPL$REPORT')
      END IF
      CALL CONSTANTS('EV',EV)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      WRITE(NFIL,'(80("#"))')
      WRITE(NFIL,FMT='(A25)')'RIXSAMPL MODULE'
      WRITE(NFIL,'(80("#"))')
      WRITE(NFIL,FMT='(A12,I10)')'NKPT:',NKPT
      WRITE(NFIL,FMT='(A12,I10)')'NKPTTOT:',NKPTTOT
      WRITE(NFIL,FMT='(A12,I10)')'NSPIN:',NSPIN
      WRITE(NFIL,FMT='(A12,3I4)')'NKDIV:',NKDIV(:)
      WRITE(NFIL,FMT='(A12,X,A)')'FLAG:',FLAG
      ! TODO: CHECK IF THIS HOLDS
      IF(EMIN.NE.-HUGE(1.D0)) WRITE(NFIL,FMT='(A12,F10.3)')'EMIN:',EMIN/EV
      IF(EMAX.NE.HUGE(1.D0)) WRITE(NFIL,FMT='(A12,F10.3)')'EMAX:',EMAX/EV
      WRITE(NFIL,FMT='(A12,F10.4)')'ELIGHT[EV]:',ELIGHT/EV
      WRITE(NFIL,FMT='(A12,3F10.4)')'NORMAL:',NORMAL(:)
      ! INCOMING LIGHT
      WRITE(NFIL,FMT='(A)') 'INCOMING LIGHT'
      WRITE(NFIL,FMT='(A12,3F10.4)')'KDIRI:',KDIRI(:)
      WRITE(NFIL,FMT='(A12,3F10.4)')'KI[1/AA]:',KI(:)*ANGSTROM
      WRITE(NFIL,FMT='(A12,3F10.4)')'XKI:',XKI(:)
      WRITE(NFIL,FMT=-'(A12,2(F8.5,SP,F8.5,"I ",S))')'POLI:',POLI(:)
      WRITE(NFIL,FMT=-'(A12,3(F8.5,SP,F8.5,"I ",S))')'POLXYZI:',POLXYZI(:)
      ! OUTGOING LIGHT
      WRITE(NFIL,FMT='(A)') 'OUTGOING LIGHT'
      WRITE(NFIL,FMT='(A12,3F10.4)')'KDIRO:',KDIRO(:)
      WRITE(NFIL,FMT='(A12,3F10.4)')'KO[1/AA]:',KO(:)*ANGSTROM
      WRITE(NFIL,FMT='(A12,3F10.4)')'XKO:',XKO(:)
      WRITE(NFIL,FMT=-'(A12,2(F8.5,SP,F8.5,"I ",S))')'POLO:',POLO(:)
      WRITE(NFIL,FMT=-'(A12,3(F8.5,SP,F8.5,"I ",S))')'POLXYZO:',POLXYZO(:)
      ! MOMENTUM TRANSFER
      WRITE(NFIL,FMT='(A,L2)')'MOMENTUM TRANSFER ACTIVE: ',TKPTSHIFT
      WRITE(NFIL,FMT='(A12,3F10.4)')'Q[1/AA]:',Q(:)*ANGSTROM
      WRITE(NFIL,FMT='(A12,3F10.4)')'QAPPROX:',QAPPROX(:)*ANGSTROM
      WRITE(NFIL,FMT='(A12,F10.4)')'QERROR:',QERROR*ANGSTROM
      WRITE(NFIL,FMT='(A12,3F10.4)')'XQ:',XQ(:)
      WRITE(NFIL,FMT='(A12,3F10.4)')'XQAPPROX:',XQAPPROX(:)

      CALL RIXSAMPL_REPORTSPEC(NFIL)
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE RIXSAMPL$REPORT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL_REPORTSPEC(NFIL)  ! MARK: RIXSAMPL_REPORTSPEC
!     **************************************************************************
!     ** REPORT SPECTRA                                                       **        
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: NKPTTOT,NSPIN,NSPEC,SPEC,FLAG
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NFIL
      INTEGER(4) :: ISPEC
      INTEGER(4) :: IKPT
      INTEGER(4) :: ISPIN
      REAL(8) :: ANGSTROM

      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      DO ISPEC=1,NSPEC
        WRITE(NFIL,'(80("-"))')
        WRITE(NFIL,FMT='(A12,I10)')'AMPLITUDE:',ISPEC
        WRITE(NFIL,FMT='(A12,X,A)')'FILE:',TRIM(SPEC(ISPEC)%FILENAME)
        WRITE(NFIL,FMT='(A12,3F14.6)')'R HOLE:',SPEC(ISPEC)%RHOLE(:)/ANGSTROM
        WRITE(NFIL,FMT='(A)') 'WARNING: COUNTING OF KPT MIGHT VARY (NKPTTOT VS. NKPT)'
        DO IKPT=1,NKPTTOT
          DO ISPIN=1,NSPIN
            IF(FLAG.EQ.+'ONEKPT') THEN
              WRITE(NFIL,FMT='(4I12)')IKPT,ISPIN, &
     &                                SPEC(ISPEC)%ONEKPT(IKPT,ISPIN)%NB1, &
     &                                SPEC(ISPEC)%ONEKPT(IKPT,ISPIN)%NOCC
            ELSE IF(FLAG.EQ.+'TWOKPT') THEN
              WRITE(NFIL,FMT='(4I12)')IKPT,ISPIN, &
     &                                SPEC(ISPEC)%TWOKPT(IKPT,ISPIN)%NB1, &
     &                                SPEC(ISPEC)%TWOKPT(IKPT,ISPIN)%NOCC
            END IF
          ENDDO ! END ISPIN
        ENDDO ! END IKPT
      ENDDO ! END ISPEC
      END SUBROUTINE RIXSAMPL_REPORTSPEC
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$GETI4(ID,VAL)  ! MARK: RIXSAMPL$GETI4
!     **************************************************************************
!     ** GET INTEGER VALUE FROM RIXSAMPL MODULE                               **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: NKPT,NSPIN,NKPTTOT,NSPEC,TFIRST
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(OUT) :: VAL
      IF(ID.EQ.'NKPT') THEN
        IF(TFIRST) THEN
          CALL ERROR$MSG('RIXSAMPL GETI4 CALLED BEFORE INITIALIZATION')
          CALL ERROR$STOP('RIXSAMPL$GETI4')
        END IF
        VAL=NKPT
      ELSE IF(ID.EQ.'NKPTTOT') THEN
        ! TODO: TFIRST DOES NOT PROTECT FROM CALLING THIS TOO EARLY
        IF(TFIRST) THEN
          CALL ERROR$MSG('RIXSAMPL GETI4 CALLED BEFORE INITIALIZATION')
          CALL ERROR$STOP('RIXSAMPL$GETI4')
        END IF
        VAL=NKPTTOT
      ELSE IF(ID.EQ.'NSPIN') THEN
        IF(TFIRST) THEN
          CALL ERROR$MSG('RIXSAMPL GETI4 CALLED BEFORE INITIALIZATION')
          CALL ERROR$STOP('RIXSAMPL$GETI4')
        END IF
        VAL=NSPIN
      ELSE IF(ID.EQ.'NSPEC') THEN
        VAL=NSPEC
      ELSE
        CALL ERROR$MSG('RIXSAMPL GETI4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$GETI4')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$GETI4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETI4(ID,VAL)  ! MARK: RIXSAMPL$SETI4
!     **************************************************************************
!     ** SET INTEGER VALUE IN RIXSAMPL MODULE                                 **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: NKPT,NSPIN,NKPTTOT,NSPEC,TFIRST,SPEC
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: VAL
      IF(ID.EQ.'NKPT') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.NKPT) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$I4VAL('NKPT',NKPT)
            CALL ERROR$I4VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETI4')
          END IF
        END IF
        NKPT=VAL
      ELSE IF(ID.EQ.'NKPTTOT') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.NKPTTOT) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$I4VAL('NKPTTOT',NKPTTOT)
            CALL ERROR$I4VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETI4')
          END IF
        END IF
        NKPTTOT=VAL
      ELSE IF(ID.EQ.'NSPIN') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.NSPIN) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$I4VAL('NSPIN',NSPIN)
            CALL ERROR$I4VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETI4')
          END IF
        END IF
        NSPIN=VAL
      ELSE IF(ID.EQ.'NSPEC') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.NSPEC) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$I4VAL('NSPEC',NSPEC)
            CALL ERROR$I4VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETI4')
          END IF
        END IF
        IF(VAL.LE.0) THEN
          CALL ERROR$MSG('NSPEC SHOULD BE POSITIVE')
          CALL ERROR$I4VAL('NSPEC',NSPEC)
          CALL ERROR$I4VAL('VAL',VAL)
          CALL ERROR$STOP('RIXSAMPL$SETI4')
        END IF
        IF(.NOT.ALLOCATED(SPEC)) THEN
          ALLOCATE(SPEC(VAL))
        END IF
        NSPEC=VAL
      ELSE
        CALL ERROR$MSG('RIXSAMPL SETI4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$SETI4')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$SETI4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETI4A(ID,LEN,VAL)  ! MARK: RIXSAMPL$SETI4A
!     **************************************************************************
!     ** SET INTEGER ARRAY VALUE IN RIXSAMPL MODULE                           **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: NKDIV,TFIRST,TOIRRKPT,NKPTTOT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      INTEGER(4), INTENT(IN) :: VAL(LEN)
      IF(ID.EQ.'NKDIV') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('NKDIV LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETI4A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(NKDIV.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETI4A')
          END IF
        END IF
        NKDIV=VAL
      ELSE IF(ID.EQ.'TOIRRKPT') THEN
        IF(LEN.NE.NKPTTOT) THEN
          CALL ERROR$MSG('TOIRRKPT LENGTH SHOULD BE NKPTTOT')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NKPTTOT',NKPTTOT)
          CALL ERROR$STOP('RIXSAMPL$SETI4A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(TOIRRKPT.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETI4A')
          END IF
        END IF
        TOIRRKPT=VAL
      ELSE
        CALL ERROR$MSG('RIXSAMPL SETI4A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$SETI4A')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$SETI4A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETR8(ID,VAL)  ! MARK: RIXSAMPL$SETR8
!     **************************************************************************
!     ** SET REAL(8) VALUE IN RIXSAMPL MODULE                                 **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: EMIN,EMAX,ELIGHT,QERROR,TFIRST
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      REAL(8), INTENT(IN) :: VAL
      IF(ID.EQ.'EMIN') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.EMIN) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$R8VAL('EMIN',EMIN)
            CALL ERROR$R8VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETR8')
          END IF
        END IF
        EMIN=VAL
      ELSE IF(ID.EQ.'EMAX') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.EMAX) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$R8VAL('EMAX',EMAX)
            CALL ERROR$R8VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETR8')
          END IF
        END IF
        EMAX=VAL
      ELSE IF(ID.EQ.'ELIGHT') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.ELIGHT) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$R8VAL('ELIGHT',ELIGHT)
            CALL ERROR$R8VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETR8')
          END IF
        END IF
        ELIGHT=VAL
      ELSE IF(ID.EQ.'QERROR') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.QERROR) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$R8VAL('QERROR',QERROR)
            CALL ERROR$R8VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETR8')
          END IF
        END IF
        QERROR=VAL
      ELSE
        CALL ERROR$MSG('RIXSAMPL SETR8 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$SETR8')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$SETR8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETR8A(ID,LEN,VAL)  ! MARK: RIXSAMPL$SETR8A
!     **************************************************************************
!     ** SET REAL(8) ARRAY VALUE IN RIXSAMPL MODULE                           **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: THIS,SELECTED,NORMAL,KDIRI,KI,XKI,KDIRO,KO,XKO, &
     &                           Q,XQ,QAPPROX,XQAPPROX,TFIRST
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      REAL(8), INTENT(IN) :: VAL(LEN)
      IF(ID.EQ.'NORMAL') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('NORMAL LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(NORMAL.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        NORMAL=VAL
      ELSE IF(ID.EQ.'KDIRI') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('KDIRI LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(KDIRI.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        KDIRI=VAL
      ELSE IF(ID.EQ.'KI') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('KI LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(KI.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        KI=VAL
      ELSE IF(ID.EQ.'XKI') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('XKI LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(XKI.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        XKI=VAL
      ELSE IF(ID.EQ.'KDIRO') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('KDIRO LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(KDIRO.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        KDIRO=VAL
      ELSE IF(ID.EQ.'KO') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('KO LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(KO.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        KO=VAL
      ELSE IF(ID.EQ.'XKO') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('XKO LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(XKO.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        XKO=VAL
      ELSE IF(ID.EQ.'Q') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('Q LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(Q.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        Q=VAL
      ELSE IF(ID.EQ.'XQ') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('XQ LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(XQ.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        XQ=VAL
      ELSE IF(ID.EQ.'QAPPROX') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('QAPPROX LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(QAPPROX.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        QAPPROX=VAL
      ELSE IF(ID.EQ.'XQAPPROX') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('XQAPPROX LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(XQAPPROX.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETR8A')
          END IF
        END IF
        XQAPPROX=VAL
      ELSE IF(ID.EQ.'RHOLE') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('RHOLE LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        IF(.NOT.SELECTED) THEN
          CALL ERROR$MSG('RIXSAMPL SETR8A CALLED WITHOUT SELECTED SPEC')
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('RIXSAMPL$SETR8A')
        END IF
        THIS%RHOLE=VAL
      ELSE
        CALL ERROR$MSG('RIXSAMPL SETR8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$SETR8A')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$SETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETC8A(ID,LEN,VAL)  ! MARK: RIXSAMPL$SETC8A
!     **************************************************************************
!     ** SET COMPLEX(8) ARRAY VALUE IN RIXSAMPL MODULE                         **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: POLI,POLXYZI,POLO,POLXYZO,TFIRST
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      COMPLEX(8), INTENT(IN) :: VAL(LEN)
      IF(ID.EQ.'POLI') THEN
        IF(LEN.NE.2) THEN
          CALL ERROR$MSG('POLI LENGTH SHOULD BE 2')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETC8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(POLI.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETC8A')
          END IF
        END IF
        POLI=VAL
      ELSE IF(ID.EQ.'POLXYZI') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('POLXYZI LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETC8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(POLXYZI.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETC8A')
          END IF
        END IF
        POLXYZI=VAL
      ELSE IF(ID.EQ.'POLO') THEN
        IF(LEN.NE.2) THEN
          CALL ERROR$MSG('POLO LENGTH SHOULD BE 2')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETC8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(POLO.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETC8A')
          END IF
        END IF
        POLO=VAL
      ELSE IF(ID.EQ.'POLXYZO') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('POLXYZO LENGTH SHOULD BE 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('RIXSAMPL$SETC8A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(POLXYZO.NE.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETC8A')
          END IF
        END IF
        POLXYZO=VAL
      ELSE
        CALL ERROR$MSG('RIXSAMPL SETC8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$SETC8A')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$SETC8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$GETCH(ID,VAL)  ! MARK: RIXSAMPL$GETCH
!     **************************************************************************
!     ** GET CHARACTER VALUE FROM RIXSAMPL MODULE                             **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: FLAG,TFIRST,THIS,SELECTED
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      CHARACTER(*), INTENT(OUT) :: VAL
      IF(ID.EQ.'FLAG') THEN
        VAL=FLAG
        IF(TFIRST) THEN
          CALL ERROR$MSG('RIXSAMPL GETCH CALLED BEFORE INITIALIZATION')
          CALL ERROR$STOP('RIXSAMPL$GETCH')
        END IF
      ELSE IF(ID.EQ.'FILENAME') THEN
        IF(.NOT.SELECTED) THEN
          CALL ERROR$MSG('RIXSAMPL GETCH CALLED WITHOUT SELECTED SPEC')
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('RIXSAMPL$GETCH')
        END IF
        VAL=TRIM(ADJUSTL(THIS%FILENAME))
      ELSE
        CALL ERROR$MSG('RIXSAMPL GETCH ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$GETCH')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$GETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETCH(ID,VAL)  ! MARK: RIXSAMPL$SETCH
!     **************************************************************************
!     ** SET CHARACTER VALUE IN RIXSAMPL MODULE                               **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: FLAG,TFIRST,THIS,SELECTED
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      CHARACTER(*), INTENT(IN) :: VAL
      IF(ID.EQ.'FLAG') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NE.FLAG) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('FLAG',FLAG)
            CALL ERROR$CHVAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETCH')
          END IF
        END IF
        FLAG=VAL
      ELSE IF(ID.EQ.'FILENAME') THEN
        IF(.NOT.SELECTED) THEN
          CALL ERROR$MSG('RIXSAMPL SETCH CALLED WITHOUT SELECTED SPEC')
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('RIXSAMPL$SETCH')
        END IF
        THIS%FILENAME=TRIM(ADJUSTL(VAL))
      ELSE
        CALL ERROR$MSG('RIXSAMPL SETCH ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$SETCH')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$SETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETL4(ID,VAL)  ! MARK: RIXSAMPL$SETL4
!     **************************************************************************
!     ** SET LOGICAL VALUE IN RIXSAMPL MODULE                                 **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: TFIRST,TKPTSHIFT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      LOGICAL(4), INTENT(IN) :: VAL
      IF(ID.EQ.'TFIRST') THEN
        TFIRST=VAL
      ELSE IF(ID.EQ.'TKPTSHIFT') THEN
        IF(.NOT.TFIRST) THEN
          IF(VAL.NEQV.TKPTSHIFT) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$L4VAL('TKPTSHIFT',TKPTSHIFT)
            CALL ERROR$L4VAL('VAL',VAL)
            CALL ERROR$STOP('RIXSAMPL$SETL4')
          END IF
        END IF
        TKPTSHIFT=VAL
      ELSE
        CALL ERROR$MSG('RIXSAMPL SETL4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$SETL4')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$SETL4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETL4A(ID,LEN,VAL)  ! MARK: RIXSAMPL$SETL4A
!     **************************************************************************
!     ** SET LOGICAL ARRAY VALUE IN RIXSAMPL MODULE                           **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: NKPT,NKPTTOT,TIRRKPT,TFIRST
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      LOGICAL(4), INTENT(IN) :: VAL(LEN)
      IF(ID.EQ.'TIRRKPT') THEN
        IF(TFIRST) ALLOCATE(TIRRKPT(NKPTTOT))
        IF(LEN.NE.NKPTTOT) THEN
          CALL ERROR$MSG('TIRRKPT LENGTH SHOULD BE NKPTTOT')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NKPTTOT',NKPTTOT)
          CALL ERROR$STOP('RIXSAMPL$SETL4A')
        END IF
        IF(.NOT.TFIRST) THEN
          IF(ANY(TIRRKPT.NEQV.VAL)) THEN
            CALL ERROR$MSG('VALUES SHOULD NOT DIFFER')
            CALL ERROR$CHVAL('ID: ',ID)
            CALL ERROR$STOP('RIXSAMPL$SETL4A')
          END IF
        END IF
        TIRRKPT=VAL
      ELSE
        CALL ERROR$MSG('RIXSAMPL SETL4A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('RIXSAMPL$SETL4A')
      END IF
      RETURN
      END SUBROUTINE RIXSAMPL$SETL4A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETONEKPT(IKPT,ISPIN,NB1,NOCC,EIG,XY)  ! MARK: RIXSAMPL$SETONEKPT
!     **************************************************************************
!     ** SET ONEKPT DATA IN RIXSAMPL MODULE                                   **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: SELECTED,THIS,NKPTTOT,NSPIN
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(IN) :: NB1
      INTEGER(4), INTENT(IN) :: NOCC
      REAL(8), INTENT(IN) :: EIG(NB1)
      COMPLEX(8), INTENT(IN) :: XY(NB1-NOCC,NOCC)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('RIXSAMPL NOT SELECTED')
        CALL ERROR$STOP('RIXSAMPL$SETONEKPT')
      END IF
      IF(.NOT.ALLOCATED(THIS%ONEKPT)) THEN
        CALL ERROR$MSG('RIXSAMPL ONEKPT NOT ALLOCATED')
        CALL ERROR$STOP('RIXSAMPL$SETONEKPT')
      END IF
      IF(IKPT.GT.NKPTTOT.OR.IKPT.LT.1) THEN
        CALL ERROR$MSG('RIXSAMPL IKPT OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT',IKPT)
        CALL ERROR$I4VAL('NKPTTOT',NKPTTOT)
        CALL ERROR$STOP('RIXSAMPL$SETONEKPT')
      END IF
      IF(ISPIN.GT.NSPIN.OR.ISPIN.LT.1) THEN
        CALL ERROR$MSG('RIXSAMPL ISPIN OUT OF RANGE')
        CALL ERROR$I4VAL('ISPIN',ISPIN)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('RIXSAMPL$SETONEKPT')
      END IF
      IF(ALLOCATED(THIS%ONEKPT(IKPT,ISPIN)%EIG)) THEN
        CALL ERROR$MSG('EIG IN ONEKPT ALREADY INITIALIZED')
        CALL ERROR$I4VAL('IKPT',IKPT)
        CALL ERROR$I4VAL('ISPIN',ISPIN)
        CALL ERROR$STOP('RIXSAMPL$SETONEKPT')
      END IF
      IF(.NOT.ALLOCATED(THIS%ONEKPT(IKPT,ISPIN)%EIG)) THEN
        ALLOCATE(THIS%ONEKPT(IKPT,ISPIN)%EIG(NB1))
      END IF
      IF(ALLOCATED(THIS%ONEKPT(IKPT,ISPIN)%XY)) THEN
        CALL ERROR$MSG('XY IN ONEKPT ALREADY INITIALIZED')
        CALL ERROR$I4VAL('IKPT',IKPT)
        CALL ERROR$I4VAL('ISPIN',ISPIN)
        CALL ERROR$STOP('RIXSAMPL$SETONEKPT')
      END IF
      IF(.NOT.ALLOCATED(THIS%ONEKPT(IKPT,ISPIN)%XY)) THEN
        ALLOCATE(THIS%ONEKPT(IKPT,ISPIN)%XY(NB1-NOCC,NOCC))
      END IF
      THIS%ONEKPT(IKPT,ISPIN)%NB1=NB1
      THIS%ONEKPT(IKPT,ISPIN)%NOCC=NOCC
      THIS%ONEKPT(IKPT,ISPIN)%EIG=EIG
      THIS%ONEKPT(IKPT,ISPIN)%XY=XY
      RETURN
      END SUBROUTINE RIXSAMPL$SETONEKPT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$SETTWOKPT(IKPT,ISPIN,NB1,NOCC,EIG,X,Y)  ! MARK: RIXSAMPL$SETTWOKPT
!     **************************************************************************
!     ** SET TWOKPT DATA IN RIXSAMPL MODULE                                    **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: SELECTED,THIS,NKPTTOT,NSPIN
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(IN) :: NB1
      INTEGER(4), INTENT(IN) :: NOCC
      REAL(8), INTENT(IN) :: EIG(NB1)
      COMPLEX(8), INTENT(IN) :: X(NB1-NOCC)
      COMPLEX(8), INTENT(IN) :: Y(NOCC)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('RIXSAMPL NOT SELECTED')
        CALL ERROR$STOP('RIXSAMPL$SETTWOKPT')
      END IF
      IF(.NOT.ALLOCATED(THIS%TWOKPT)) THEN
        CALL ERROR$MSG('RIXSAMPL TWOKPT NOT ALLOCATED')
        CALL ERROR$STOP('RIXSAMPL$SETTWOKPT')
      END IF
      IF(IKPT.GT.NKPTTOT.OR.IKPT.LT.1) THEN
        CALL ERROR$MSG('RIXSAMPL IKPT OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT',IKPT)
        CALL ERROR$I4VAL('NKPTTOT',NKPTTOT)
        CALL ERROR$STOP('RIXSAMPL$SETTWOKPT')
      END IF
      IF(ISPIN.GT.NSPIN.OR.ISPIN.LT.1) THEN
        CALL ERROR$MSG('RIXSAMPL ISPIN OUT OF RANGE')
        CALL ERROR$I4VAL('ISPIN',ISPIN)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('RIXSAMPL$SETTWOKPT')
      END IF
      IF(.NOT.ALLOCATED(THIS%TWOKPT(IKPT,ISPIN)%EIG)) THEN
        ALLOCATE(THIS%TWOKPT(IKPT,ISPIN)%EIG(NB1))
      END IF
      IF(.NOT.ALLOCATED(THIS%TWOKPT(IKPT,ISPIN)%X)) THEN
        ALLOCATE(THIS%TWOKPT(IKPT,ISPIN)%X(NB1-NOCC))
      END IF
      IF(.NOT.ALLOCATED(THIS%TWOKPT(IKPT,ISPIN)%Y)) THEN
        ALLOCATE(THIS%TWOKPT(IKPT,ISPIN)%Y(NOCC))
      END IF
      THIS%TWOKPT(IKPT,ISPIN)%NB1=NB1
      THIS%TWOKPT(IKPT,ISPIN)%NOCC=NOCC
      THIS%TWOKPT(IKPT,ISPIN)%EIG=EIG
      THIS%TWOKPT(IKPT,ISPIN)%X=X
      THIS%TWOKPT(IKPT,ISPIN)%Y=Y
      RETURN
      END SUBROUTINE RIXSAMPL$SETTWOKPT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$CALCULATE  ! MARK: RIXSAMPL$CALCULATE
!     **************************************************************************
!     ** CALCULATE CROSS-SECTION FROM AMPLITUDES                              **
!     **************************************************************************
      IMPLICIT NONE
                          CALL TRACE$PUSH('RIXSAMPL$CALCULATE')
      
      CALL RIXSAMPL_CALCULATEONEKPT

      CALL RIXSAMPL_CALCULATETWOKPT

                          CALL TRACE$POP
      RETURN
      END SUBROUTINE RIXSAMPL$CALCULATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL_CALCULATEONEKPT  ! MARK: RIXSAMPL_CALCULATEONEKPT
!     **************************************************************************
!     ** CALCULATE CROSS-SECTION FROM AMPLITUDES FOR ONEKPT                   **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: FLAG,NSPEC,RES,NKPTTOT,NSPIN,Q,THIS
      USE STRINGS_MODULE
      IMPLICIT NONE
      COMPLEX(8), PARAMETER :: CI=(0.D0,1.D0)
      INTEGER(4) :: ISPEC
      INTEGER(4) :: IKPT
      INTEGER(4) :: ISPIN
      COMPLEX(8) :: EXPFAC
      INTEGER(4) :: NB1,NOCC
! WARNING: REQUIRES ALL (IKPT,ISPIN) TO HAVE THE SAME NB1,NOCC FOR ALL SPEC
!          SHOULD BE GIVEN AS ALL REFERENCE THE SAME GROUND STATE SIMULATION
      IF(FLAG.NE.+'ONEKPT') RETURN
                          CALL TRACE$PUSH('RIXSAMPL_CALCULATEONEKPT')
      IF(ALLOCATED(RES)) THEN
        CALL ERROR$MSG('RES ARRAY ALREADY ALLOCATED')
        CALL ERROR$STOP('RIXSAMPL_CALCULATEONEKPT')
      END IF
      ALLOCATE(RES(NKPTTOT,NSPIN))

      ! ALLOCATE AND INITIALIZE DATA
      CALL RIXSAMPL$ISELECT(1)
      DO IKPT=1,NKPTTOT
        DO ISPIN=1,NSPIN
          NB1=THIS%ONEKPT(IKPT,ISPIN)%NB1
          NOCC=THIS%ONEKPT(IKPT,ISPIN)%NOCC
          RES(IKPT,ISPIN)%NB1=NB1
          RES(IKPT,ISPIN)%NOCC=NOCC
          ALLOCATE(RES(IKPT,ISPIN)%EIG(NB1))
          RES(IKPT,ISPIN)%EIG(:)=THIS%ONEKPT(IKPT,ISPIN)%EIG(:)
          ALLOCATE(RES(IKPT,ISPIN)%XY(NB1-NOCC,NOCC))
          RES(IKPT,ISPIN)%XY(:,:)=(0.D0,0.D0)
          ALLOCATE(RES(IKPT,ISPIN)%ABSXY(NB1-NOCC,NOCC))
        ENDDO
      ENDDO
      CALL RIXSAMPL$UNSELECT


      DO ISPEC=1,NSPEC
        CALL RIXSAMPL$ISELECT(ISPEC)
        ! EXP(I*Q*D)
        EXPFAC=EXP(CI*DOT_PRODUCT(Q,THIS%RHOLE))
        DO IKPT=1,NKPTTOT
          DO ISPIN=1,NSPIN
            RES(IKPT,ISPIN)%XY(:,:)=RES(IKPT,ISPIN)%XY(:,:)+ &
     &                              EXPFAC*THIS%ONEKPT(IKPT,ISPIN)%XY(:,:)
          ENDDO
        ENDDO
        CALL RIXSAMPL$UNSELECT
      ENDDO

      DO IKPT=1,NKPTTOT
        DO ISPIN=1,NSPIN
          RES(IKPT,ISPIN)%ABSXY(:,:)=REAL(RES(IKPT,ISPIN)%XY,KIND=8)**2+AIMAG(RES(IKPT,ISPIN)%XY)**2
        ENDDO
      ENDDO
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE RIXSAMPL_CALCULATEONEKPT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL_CALCULATETWOKPT  ! MARK: RIXSAMPL_CALCULATETWOKPT
!     **************************************************************************
!     ** CALCULATE CROSS-SECTION FROM AMPLITUDES FOR TWOKPT                   **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: FLAG,RES,NKPTTOT,NSPIN,NSPEC,THIS,Q
      USE STRINGS_MODULE
      IMPLICIT NONE
      COMPLEX(8), PARAMETER :: CI=(0.D0,1.D0)
      INTEGER(4) :: ISPEC
      INTEGER(4) :: IKPT
      INTEGER(4) :: ISPIN
      COMPLEX(8) :: EXPFAC
      COMPLEX(8) :: CVAR
      INTEGER(4) :: NB1,NOCC
      INTEGER(4) :: IEMP,IOCC
! WARNING: REQUIRES ALL (IKPT,ISPIN) TO HAVE THE SAME NB1,NOCC FOR ALL SPEC
!          SHOULD BE GIVEN AS ALL REFERENCE THE SAME GROUND STATE SIMULATION
      IF(FLAG.NE.+'TWOKPT') RETURN
                          CALL TRACE$PUSH('RIXSAMPL_CALCULATETWOKPT')

      IF(ALLOCATED(RES)) THEN
        CALL ERROR$MSG('RES ARRAY ALREADY ALLOCATED')
        CALL ERROR$STOP('RIXSAMPL_CALCULATEONEKPT')
      END IF
      ALLOCATE(RES(NKPTTOT,NSPIN))

      ! ALLOCATE AND INITIALIZE DATA
      CALL RIXSAMPL$ISELECT(1)
      DO IKPT=1,NKPTTOT
        DO ISPIN=1,NSPIN
          NB1=THIS%TWOKPT(IKPT,ISPIN)%NB1
          NOCC=THIS%TWOKPT(IKPT,ISPIN)%NOCC
          RES(IKPT,ISPIN)%NB1=NB1
          RES(IKPT,ISPIN)%NOCC=NOCC
          ALLOCATE(RES(IKPT,ISPIN)%EIG(NB1))
          RES(IKPT,ISPIN)%EIG(:)=THIS%TWOKPT(IKPT,ISPIN)%EIG(:)
          ALLOCATE(RES(IKPT,ISPIN)%XY(NB1-NOCC,NOCC))
          RES(IKPT,ISPIN)%XY(:,:)=(0.D0,0.D0)
          ALLOCATE(RES(IKPT,ISPIN)%ABSXY(NB1-NOCC,NOCC))
        ENDDO
      ENDDO
      CALL RIXSAMPL$UNSELECT

      DO ISPEC=1,NSPEC
        CALL RIXSAMPL$ISELECT(ISPEC)
        ! EXP(I*Q*D)
        EXPFAC=EXP(CI*DOT_PRODUCT(Q,THIS%RHOLE))
        DO IKPT=1,NKPTTOT
          call trace$i4val('ikpt',ikpt)
          DO ISPIN=1,NSPIN
            call trace$i4val('ispin',ISPIN)
            NB1=RES(IKPT,ISPIN)%NB1
            NOCC=RES(IKPT,ISPIN)%NOCC
            DO IOCC=1,NOCC
              DO IEMP=1,NB1-NOCC
                CVAR=EXPFAC*THIS%TWOKPT(IKPT,ISPIN)%X(IEMP)*THIS%TWOKPT(IKPT,ISPIN)%Y(IOCC)
                RES(IKPT,ISPIN)%XY(IEMP,IOCC)=RES(IKPT,ISPIN)%XY(IEMP,IOCC)+CVAR
              ENDDO ! END IEMP
            ENDDO ! END IOCC
          ENDDO ! END ISPIN
        ENDDO ! END IKPT
        CALL RIXSAMPL$UNSELECT
      ENDDO ! END ISPEC

      DO IKPT=1,NKPTTOT
        DO ISPIN=1,NSPIN
          RES(IKPT,ISPIN)%ABSXY(:,:)=REAL(RES(IKPT,ISPIN)%XY,KIND=8)**2+AIMAG(RES(IKPT,ISPIN)%XY)**2
        ENDDO
      ENDDO
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE RIXSAMPL_CALCULATETWOKPT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$ISELECT(ISPEC)  ! MARK: RIXSAMPL$ISELECT
!     **************************************************************************
!     ** SELECT SPEC FOR RIXSAMPL DATA                                       **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: SELECTED,THIS,SPEC,NSPEC
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: ISPEC
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('RIXSAMPL ALREADY SELECTED')
        CALL ERROR$STOP('RIXSAMPL$ISELECT')
      END IF
      IF(.NOT.ALLOCATED(SPEC)) THEN
        CALL ERROR$MSG('RIXSAMPL SPEC NOT ALLOCATED')
        CALL ERROR$STOP('RIXSAMPL$ISELECT')
      END IF
      IF(NSPEC.EQ.-1) THEN
        CALL ERROR$MSG('RIXSAMPL NSPEC NOT SET')
        CALL ERROR$STOP('RIXSAMPL$ISELECT')
      END IF
      IF(ISPEC.GT.NSPEC.OR.ISPEC.LT.1) THEN
        CALL ERROR$MSG('RIXSAMPL ISPEC OUT OF RANGE')
        CALL ERROR$I4VAL('ISPEC',NSPEC)
        CALL ERROR$I4VAL('NSPEC',NSPEC)
        CALL ERROR$STOP('RIXSAMPL$ISELECT')
      END IF
      THIS=>SPEC(ISPEC)
      SELECTED=.TRUE.
      RETURN
      END SUBROUTINE RIXSAMPL$ISELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RIXSAMPL$UNSELECT  ! MARK: RIXSAMPL$UNSELECT
!     **************************************************************************
!     ** UNSELECT SPEC FOR RIXSAMPL DATA                                     **
!     **************************************************************************
      USE RIXSAMPL_MODULE, ONLY: SELECTED,THIS,SPEC
      IMPLICIT NONE
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT UNSELECT SPECTRUM IF NONE SELECTED')
        CALL ERROR$STOP('RIXSAMPL$UNSELECT')
      END IF
      NULLIFY(THIS)
      SELECTED=.FALSE.
      RETURN
      END SUBROUTINE RIXSAMPL$UNSELECT
!      
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INITIALIZEFILEHANDLER
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: CNTLNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: NARGS
                          CALL TRACE$PUSH('INITIALIZEFILEHANDLER')
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE XRAY TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GET_COMMAND_ARGUMENT(1,CNTLNAME)
      ISVAR=INDEX(CNTLNAME,-'.XCCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=CNTLNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('XCCNTL',.FALSE.,CNTLNAME)
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE INITIALIZEFILEHANDLER
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STANDARDFILES
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(32)        :: ID
      INTEGER(4)          :: THISTASK
      INTEGER(4)          :: NTASKS
                                   CALL TRACE$PUSH('STANDARDFILES')
      ! ========================================================================
      ! == SET STANDARD FILENAMES                                             ==
      ! ========================================================================
      ! ==  ERROR FILE =========================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XCERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
      ! ==  PROTOCOL FILE ======================================================
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      ID=+'PROT'
      IF(THISTASK.GT.1) THEN
        CALL FILEHANDLER$SETFILE(ID,.FALSE.,-'/DEV/NULL')
      ELSE
        CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XCPROT')
      END IF
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
      ! ==  CONTROL FILE  == ===================================================
      ID=+'XCCNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XCCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XCCNTL$READ(NFIL)  ! MARK: XCCNTL$READ
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XCCNTL_MODULE, ONLY: LL_CNTL
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NFIL
      LOGICAL(4) :: TCHK
                          CALL TRACE$PUSH('XCCNTL$READ')
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$MARK(LL_CNTL,1)
      CALL LINKEDLIST$EXISTL(LL_CNTL,'XCCNTL',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('CONTROL FILE DOES NOT CONTAIN !XCCNTL BLOCK')
        CALL ERROR$STOP('XCCNTL$READ')
      END IF
                          CALL TRACE$POP
      END SUBROUTINE XCCNTL$READ
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XCCNTL$AMPLITUDE  ! MARK: XCCNTL$AMPLITUDE
!     **************************************************************************
!     ** READ !XCCNTL!AMPLITUDE BLOCKS FROM CONTROL FILE                      **
!     **************************************************************************
      USE XCCNTL_MODULE, ONLY: LL_CNTL
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      LOGICAL(4) :: TCHK
      INTEGER(4) :: NSPEC
      INTEGER(4) :: ISPEC
      CHARACTER(256) :: FILENAME
                          CALL TRACE$PUSH('XCCNTL$AMPLITUDE')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'XCCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'AMPLITUDE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('CONTROL FILE DOES NOT CONTAIN !XCCNTL!AMPLITUDE BLOCKS')
        CALL ERROR$STOP('XCCNTL$AMPLITUDE')
      END IF

      CALL LINKEDLIST$NLISTS(LL_CNTL,'AMPLITUDE',NSPEC)
      CALL RIXSAMPL$SETI4('NSPEC',NSPEC)

      DO ISPEC=1,NSPEC
        CALL RIXSAMPL$ISELECT(ISPEC)
        CALL LINKEDLIST$SELECT(LL_CNTL,'~')
        CALL LINKEDLIST$SELECT(LL_CNTL,'XCCNTL')
        CALL LINKEDLIST$SELECT(LL_CNTL,'AMPLITUDE',ISPEC)

        ! READ FILE NAME
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!AMPLITUDE:FILE NOT FOUND')
          CALL ERROR$I4VAL('ISPEC',ISPEC)
          CALL ERROR$STOP('XCCNTL$AMPLITUDE')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FILENAME)
        CALL RIXSAMPL$SETCH('FILENAME',FILENAME)
        CALL RIXSAMPL$UNSELECT
      ENDDO
    
                          CALL TRACE$POP
      END SUBROUTINE XCCNTL$AMPLITUDE