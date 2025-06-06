!     ==  OPEN QUESTIONS  ======================================================
! WHAT HAPPENS FOR DIFFERENT AMOUNT OF BANDS AT DIFFERENT K-POINTS?

!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE XCNTL_MODULE  ! MARK: XCNTL_MODULE
!     **************************************************************************
!     ** XCNTL MODULE FOR PAW XRAY TOOL                                       **
!     **************************************************************************
        USE LINKEDLIST_MODULE, ONLY: LL_TYPE
        TYPE(LL_TYPE) :: LL_CNTL
        SAVE
      END MODULE XCNTL_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE KSMAP_MODULE  ! MARK: KSMAP_MODULE
!     **************************************************************************
!     ** KSMAP MODULE FOR PAW XRAY TOOL                                       **
!     ** CONTROLS THE PARALLELIZATION OF K-POINTS AND SPINS                   **
!     **************************************************************************
      LOGICAL(4) :: INITIALIZED=.FALSE. ! INITIALIZED
      INTEGER(4) :: NKPTG  ! #(KPOINTS TOTAL)
      INTEGER(4) :: NSPING ! #(SPINS TOTAL)
      INTEGER(4) :: RTASK  ! TASK RESPONSIBLE FOR READ/WRITE
      INTEGER(4), ALLOCATABLE :: KSMAP(:,:) ! (NKPT,NSPIN) K-POINT AND SPIN RESPONSIBILITY OF TASKS
      END MODULE KSMAP_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE SIMULATION_MODULE  ! MARK: SIMULATION_MODULE
!     **************************************************************************
!     ** SIMULATION MODULE FOR PAW XRAY TOOL                                  **
!     **************************************************************************
      TYPE SIMULATION_TYPE
      CHARACTER(6) :: ID ! SIMULATION ID
      LOGICAL(4) :: INITIALIZED=.FALSE. ! INITIALIZED
      CHARACTER(256) :: FILE ! FILENAME
      INTEGER(4) :: NAT ! #(ATOMS)
      INTEGER(4) :: NSP ! #(SETUPS)
      INTEGER(4) :: NKPT ! #(KPOINTS)
      INTEGER(4) :: NSPIN ! #(SPINS)
      INTEGER(4) :: NDIM ! #(DIMENSIONS)
      INTEGER(4) :: NPRO ! #(PROJECTIONS)
      INTEGER(4) :: LNXX
      CHARACTER(6) :: FLAG
      INTEGER(4), ALLOCATABLE :: LMNX(:) ! (NSP)
      INTEGER(4), ALLOCATABLE :: LNX(:) ! (NSP)
      INTEGER(4), ALLOCATABLE :: LOX(:,:) ! (LNXX,NSP)
      INTEGER(4), ALLOCATABLE :: MAP(:,:) ! (NAT,LNXX) INDEX-1 OF ATOM AND LN
      LOGICAL(4) :: TINV ! INVERSION SYMMETRY
      INTEGER(4) :: NKDIV(3) ! K-POINT DIVISIONS
      INTEGER(4) :: ISHIFT(3) ! K-POINT SHIFTS
      REAL(8) :: RNTOT
      REAL(8) :: NEL
      REAL(8) :: ETOT=HUGE(1.D0)  ! TOTAL ENERGY
      REAL(8) :: EDFT=HUGE(1.D0)  ! TOTAL DFT ENERGY
      REAL(8) :: ECORE=HUGE(1.D0)  ! TOTAL CORE ENERGY
      INTEGER(4) :: SPACEGROUP
      LOGICAL(4) :: TSHIFT
      REAL(8) :: RBAS(3,3) ! BASIS VECTORS
      REAL(8) :: VCELL=-HUGE(1.D0) ! CELL VOLUME
      REAL(8), ALLOCATABLE :: R(:,:) ! (3,NAT) ATOMIC POSITIONS
      CHARACTER(16), ALLOCATABLE :: ATOMID(:) ! (NAT) ATOM IDENTIFIERS
      INTEGER(4), ALLOCATABLE :: ISPECIES(:) ! (NAT) SPECIES INDEX
      REAL(8), ALLOCATABLE :: XK(:,:) ! (3,NKPT) K-POINTS IN REL. COORD.
      REAL(8), ALLOCATABLE :: WKPT(:) ! (NKPT) K-POINT WEIGHTS
      COMPLEX(8), ALLOCATABLE :: DENMAT(:,:,:,:) ! (LMNXX,LMNXX,NDIM,NAT) DENSITY MATRIX
      END TYPE SIMULATION_TYPE
      
      TYPE(SIMULATION_TYPE) :: GROUND ! GROUND STATE SIMULATION (1)
      TYPE(SIMULATION_TYPE) :: EXCITE ! EXCITED STATE SIMULATION (2)
      TYPE(SIMULATION_TYPE), POINTER :: THIS ! POINTER TO CURRENT SIMULATION
      LOGICAL(4), SAVE :: SELECTED=.FALSE. ! SIMULATION SELECTED
      END MODULE SIMULATION_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE STATE_MODULE  ! MARK: STATE_MODULE
!     **************************************************************************
!     ** STATE MODULE FOR PAW XRAY TOOL                                       **
!     **************************************************************************
      ! WARNING: WHY ARE POINTERS USED HERE?
      TYPE STATE_TYPE
      ! IF OCCUPATION > OCCPERCENT*MAXOCC THEN OCCUPIED 
      ! DEACTIVATED IF EFERMI GIVEN, MAXOCC=2.0/1.0 FOR NSPIN=1/2
      REAL(8) :: OCCPERCENT=0.5D0 
      INTEGER(4)          :: NB=-HUGE(1) ! #(STATES)
      REAL(8), POINTER    :: EIG(:) ! (NB) EIGENVALUES
      REAL(8), POINTER    :: OCC(:) ! (NB) OCCUPATIONS
      INTEGER(4)          :: NOCC=-HUGE(1) ! #(OCCUPIED STATES)
      COMPLEX(8), POINTER :: PROJ(:,:,:) ! (NDIM,NB,NPRO) PROJECTIONS
      END TYPE STATE_TYPE

      INTEGER(4) :: NKPTG
      INTEGER(4) :: NSPING
      INTEGER(4) :: NDIM
      INTEGER(4) :: NPRO
      LOGICAL(4) :: INITIALIZED=.FALSE. ! INITIALIZED
      TYPE(STATE_TYPE), ALLOCATABLE :: GROUND(:,:) ! (NKPT,NSPIN) GROUND STATE
      TYPE(STATE_TYPE), ALLOCATABLE :: EXCITE(:,:) ! (NKPT,NSPIN) EXCITED STATE
      TYPE(STATE_TYPE), POINTER :: THIS ! POINTER TO CURRENT STATE ARRAY
      LOGICAL(4), SAVE :: SELECTED=.FALSE. ! STATE SELECTED
      END MODULE STATE_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE XSETUP_MODULE  ! MARK: XSETUP_MODULE
!     **************************************************************************
!     ** XSETUP MODULE FOR PAW XRAY TOOL                                      **
!     ** NOT AVAILABLE IF START FROM RESTART FILE                             **
!     **************************************************************************
      TYPE XSETUP_TYPE
      INTEGER(4) :: GID ! GRID ID FOR RADIAL MESH
      REAL(8) :: ECORE ! CORE ENERGY
      INTEGER(4) :: LNX
      REAL(8), ALLOCATABLE :: PSPHI(:,:) ! PSEUDO PARTIAL WAVE
      REAL(8), ALLOCATABLE :: AEPHI(:,:) ! AE PARTIAL WAVE
      INTEGER(4) :: NBATOM ! #(ATOMIC WAVE FUNCTIONS)
      INTEGER(4), ALLOCATABLE :: LATOM(:) ! (NBATOM) L QUANTUM NUMBER
      REAL(8), ALLOCATABLE :: AEPSI(:,:) ! (NR,NBATOM) AE ATOMIC WAVE FUNCTION
      END TYPE XSETUP_TYPE

      TYPE(XSETUP_TYPE), ALLOCATABLE :: GROUND(:) ! GROUND STATE SETUPS
      TYPE(XSETUP_TYPE), ALLOCATABLE :: EXCITE(:) ! EXCITED STATE SETUPS
      TYPE(XSETUP_TYPE), POINTER :: THIS ! POINTER TO CURRENT SETUP ARRAY
      LOGICAL(4), SAVE :: SELECTED=.FALSE. ! SETUP SELECTED
      LOGICAL(4), SAVE :: TRSTRT=.FALSE. ! START FROM RESTART FILE
      END MODULE XSETUP_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE OVERLAP_MODULE  ! MARK: OVERLAP_MODULE
!     **************************************************************************
!     ** OVERLAP MODULE FOR PAW XRAY TOOL                                     **
!     **************************************************************************
      ! TODO: CHECK WHAT IS ACTUALLY NECESSARY
      TYPE OVERLAP_TYPE
      INTEGER(4) :: NB1 ! #(BANDS IN SIMULATION 1)
      INTEGER(4) :: NB2 ! #(BANDS IN SIMULATION 2)
      INTEGER(4) :: NOCC ! #(OCCUPIED BANDS)
      COMPLEX(8), ALLOCATABLE :: PW(:,:) ! (NB2,NB1) PLANE WAVE OVERLAP MATRIX
      COMPLEX(8), ALLOCATABLE :: AUG(:,:) ! (NB2,NB1) AUGMENTATION OVERLAP MATRIX
      COMPLEX(8), ALLOCATABLE :: OV(:,:,:) ! (NB2,NB1) OVERLAP MATRIX
      ! IT IS CHECKED THAT NOCC1=NOCC2
      ! OVERLAP MATRIX
      ! | A | B |
      ! |---+---|
      ! | C | D |
      COMPLEX(8), ALLOCATABLE :: A(:,:) ! (NOCC,NOCC)
      COMPLEX(8), ALLOCATABLE :: B(:,:) ! (NOCC,NB1-NOCC)
      COMPLEX(8), ALLOCATABLE :: C(:,:) ! (NB2-NOCC,NOCC)
      COMPLEX(8), ALLOCATABLE :: D(:,:) ! (NB2-NOCC,NB1-NOCC)
      ! AINV=A^(-1)
      COMPLEX(8), ALLOCATABLE :: AINV(:,:) ! (NOCC,NOCC)
      COMPLEX(8) :: ADET ! DETERMINANT OF A
      COMPLEX(8), ALLOCATABLE :: DIPOLE(:,:) ! (3,NB2) DIPOLE MATRIX ELEMENTS
      END TYPE OVERLAP_TYPE
      
      TYPE(OVERLAP_TYPE), ALLOCATABLE, TARGET :: OVLARR(:,:) ! (NKPT,NSPIN) OVERLAP ARRAYS
      TYPE(OVERLAP_TYPE), POINTER :: THIS ! POINTER TO CURRENT OVERLAP
      INTEGER(4) :: NKPTG  ! #(KPOINTS TOTAL)
      INTEGER(4) :: NSPING ! #(SPINS TOTAL)
      LOGICAL(4), SAVE :: SELECTED=.FALSE. ! OVERLAP SELECTED
      LOGICAL(4) :: INITIALIZED=.FALSE. ! INITIALIZED
      END MODULE OVERLAP_MODULE


!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE PAW_XRAY_MODULE  ! MARK: PAW_XRAY_MODULE
      TYPE SETTINGS_TYPE
        CHARACTER(32) :: COREHOLE ! ATOM ID OF CORE HOLE
        INTEGER(4) :: IATOM ! ATOM INDEX WITH CORE HOLE IN SIM1 (FOR SIM2 USE ATOMMAP)
        INTEGER(4) :: NCORE ! N QUANTUM NUMBER OF CORE HOLE
        INTEGER(4) :: LCORE ! L QUANTUM NUMBER OF CORE HOLE
      END TYPE SETTINGS_TYPE

      TYPE(SETTINGS_TYPE) :: SETTINGS

      REAL(8), ALLOCATABLE :: S(:,:,:) ! (NAT,LNXX1,LNXX2) ATOMIC OVERLAP MATRIX
      INTEGER(4), ALLOCATABLE :: ATOMMAP(:) ! (NAT) ATOM INDEX MAP BRA -> KET

      INTEGER(4), ALLOCATABLE :: KSMAP(:,:) ! (NKPT,NSPIN) K AND S RESPONSIBILITY OF TASKS
      INTEGER(4) :: NKPTG  ! #(KPOINTS TOTAL)
      INTEGER(4) :: NSPING ! #(SPINS TOTAL)
      INTEGER(4) :: RTASK  ! TASK RESPONSIBLE FOR READ/WRITE
      END MODULE PAW_XRAY_MODULE

!     ==========================================================================
!     ==                    PROGRAM PAW_XRAY                                  ==
!     ==========================================================================
      PROGRAM PAW_XRAY  ! MARK: PAW_XRAY
      USE LINKEDLIST_MODULE
      USE CLOCK_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NFIL
      CHARACTER(32) :: DATIME
      INTEGER(4) :: NTASKS
      INTEGER(4) :: THISTASK
!     **************************************************************************
      CALL MPE$INIT
                          CALL TRACE$PUSH('PAW_XRAY')
                          CALL TIMING$START
      CALL MPE$QUERY('~',NTASKS,THISTASK)
!     INITIALIZE FILES
      CALL INITIALIZEFILEHANDLER
      CALL FILEHANDLER$UNIT('PROT',NFIL)
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL,FMT='()')
        CALL CPPAW_WRITEVERSION(NFIL)
        CALL CLOCK$NOW(DATIME)        
        WRITE(NFIL,FMT='(80("="))')
        WRITE(NFIL,FMT='(80("="),T15,"  PROGRAM STARTED ",A32,"  ")')DATIME
        WRITE(NFIL,FMT='(80("="))')
        WRITE(NFIL,FMT='(A,I3)')'NUMBER OF TASKS: ',NTASKS
      END IF
!     ==========================================================================
!     == READ XRAY CONTROL FILE                                               ==
!     ==========================================================================
                          CALL TIMING$CLOCKON('XCNTL')
      CALL FILEHANDLER$UNIT('XCNTL',NFIL)
      CALL XCNTL$READ(NFIL)
      ! FINISH READING CONTROL FILE
                          CALL TIMING$CLOCKOFF('XCNTL')
      ! READ RESTART OR INITIALIZE SIMULATIONS

      ! CALCULATE OVERLAP

      ! CALCULATE DIPOLE MATRIX ELEMENTS

      ! ...












                          CALL TIMING$PRINT('~',NFIL)
      CALL MPE$CLOCKREPORT(NFIL)
      IF(THISTASK.EQ.1) THEN
        CALL CLOCK$NOW(DATIME)
        WRITE(NFIL,FMT='(80("="))')
        WRITE(NFIL,FMT='(80("="),T15,"  PROGRAM FINISHED ",A32,"  ")')DATIME
        WRITE(NFIL,FMT='(80("="))')
      END IF
                          CALL TRACE$POP
      CALL ERROR$NORMALSTOP
      STOP

                          CALL TRACE$POP
      END PROGRAM PAW_XRAY

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
      ! **************************************************************************
                          CALL TRACE$PUSH('INITIALIZEFILEHANDLER')
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE XRAY TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GET_COMMAND_ARGUMENT(1,CNTLNAME)
      ISVAR=INDEX(CNTLNAME,-'.XCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=CNTLNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('XCNTL',.FALSE.,CNTLNAME)
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
      ! **************************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
      ! ========================================================================
      ! == SET STANDARD FILENAMES                                             ==
      ! ========================================================================
      ! ==  ERROR FILE =========================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
      ! ==  PROTOCOL FILE ======================================================
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      ID=+'PROT'
      ! TODO: THIS ASSUMES THAT TASKS 1 IS ALWAYS READTASK
      IF(THISTASK.GT.1) THEN
        CALL FILEHANDLER$SETFILE(ID,.FALSE.,-'/DEV/NULL')
      ELSE
        CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XPROT')
      END IF
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
      ! ==  CONTROL FILE  == ===================================================
      ID=+'XCNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES
!

!     ==========================================================================
!     ==========================================================================
!     ==                    XCNTL MODULE FUNCTIONS                            ==
!     ==========================================================================
!     ==========================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XCNTL$READ(NFIL)  ! MARK: XCNTL$READ
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XCNTL_MODULE, ONLY: LL_CNTL
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NFIL
      LOGICAL(4) :: TCHK
                          CALL TRACE$PUSH('XCNTL$READ')
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$MARK(LL_CNTL,1)
      CALL LINKEDLIST$EXISTL(LL_CNTL,'XCNTL',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('CONTROL FILE DOES NOT CONTAIN !XCNTL BLOCK')
        CALL ERROR$STOP('XCNTL$READ')
      END IF
                          CALL TRACE$POP
      END SUBROUTINE XCNTL$READ
!
!     ==========================================================================
!     ==========================================================================
!     ==                 SIMULATION MODULE FUNCTIONS                          ==
!     ==========================================================================
!     ==========================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$CONSISTENCY  ! MARK: SIMULATION$CONSISTENCY
!     **************************************************************************
!     ** CHECK CONSISTENCY OF THE TWO SIMULATIONS                             **
!     **************************************************************************
! TODO: IMPLEMENT CHECK IF BOTH SIMULATIONS HAVE THE SAME PARAMETERS FOR CERTAIN
!       PROPERTIES, SUCH AS THE NUMBER OF ATOMS, SPINS, K-POINTS, ETC.
      END SUBROUTINE SIMULATION$CONSISTENCY
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$SELECT(ID)  ! MARK: SIMULATION$SELECT
!     **************************************************************************
!     ** SELECT SIMULATION IN SIMULATION MODULE                               **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT SELECT SIMULATION WHEN ALREADY SELECTED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%FLAG)
        CALL ERROR$CHVAL('ID TO SELECT: ',ID)
        CALL ERROR$STOP('SIMULATION$SELECT')
      END IF
      IF(ID.EQ.'GROUND') THEN
        THIS=>GROUND
      ELSE IF(ID.EQ.'EXCITE') THEN
        THIS=>EXCITE
      ELSE
        CALL ERROR$MSG('SIMULATION ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID TO SELECT: ',ID)
        CALL ERROR$STOP('SIMULATION$SELECT')
      END IF
      SELECTED=.TRUE.
      RETURN
      END SUBROUTINE SIMULATION$SELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$UNSELECT()  ! MARK: SIMULATION$UNSELECT
!     **************************************************************************
!     ** UNSELECT SIMULATION IN SIMULATION MODULE                             **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT UNSELECT SIMULATION WHEN NOT SELECTED')
        CALL ERROR$STOP('SIMULATION$UNSELECT')
      END IF
      SELECTED=.FALSE.
      NULLIFY(THIS)
      RETURN
      END SUBROUTINE SIMULATION$UNSELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$ISELECT(I)  ! MARK: SIMULATION$ISELECT
!     **************************************************************************
!     ** SELECT SIMULATION IN SIMULATION MODULE BY INDEX                      **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: I
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT SELECT SIMULATION WHEN ALREADY SELECTED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%FLAG)
        CALL ERROR$CHVAL('INDEX TO SELECT: ',I)
        CALL ERROR$STOP('SIMULATION$ISELECT')
      END IF
      IF(I.EQ.1) THEN
        THIS=>GROUND
      ELSE IF(I.EQ.2) THEN
        THIS=>EXCITE
      ELSE
        CALL ERROR$MSG('SIMULATION INDEX NOT RECOGNIZED')
        CALL ERROR$CHVAL('INDEX TO SELECT: ',I)
        CALL ERROR$STOP('SIMULATION$ISELECT')
      END IF
      SELECTED=.TRUE.
      RETURN
      END SUBROUTINE SIMULATION$ISELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$INIT(NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG)  ! MARK: SIMULATION$INIT
!     **************************************************************************
!     ** INITIALIZE SIMULATION MODULE                                         **
!     ** NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG ARE SIMULATION PARAMETERS     **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NAT
      INTEGER(4), INTENT(IN) :: NSP
      INTEGER(4), INTENT(IN) :: NKPT
      INTEGER(4), INTENT(IN) :: NSPIN
      INTEGER(4), INTENT(IN) :: NDIM
      INTEGER(4), INTENT(IN) :: NPRO
      INTEGER(4), INTENT(IN) :: LNXX
      CHARACTER(6), INTENT(IN) :: FLAG
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$INIT')
      END IF
      IF(THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('SIMULATION ALREADY INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP
      END IF
      IF(ID.NE.'GROUND'.AND.ID.NE.'EXCITE') THEN
        CALL ERROR$MSG('SIMULATION ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID TO INITIALIZE: ',ID)
        CALL ERROR$STOP('SIMULATION$INIT')
      END IF
      THIS%NAT=NAT
      THIS%NSP=NSP
      THIS%NKPT=NKPT
      THIS%NSPIN=NSPIN
      THIS%NDIM=NDIM
      THIS%NPRO=NPRO
      THIS%LNXX=LNXX
      THIS%FLAG=FLAG
      ! ALLOCATE(THIS%LMNX(NSP))
      ALLOCATE(THIS%LNX(NSP))
      ALLOCATE(THIS%LOX(LNXX,NSP))
      ! ALLOCATE(THIS%MAP(NAT,LNXX))
      ALLOCATE(THIS%NKDIV(3))
      ALLOCATE(THIS%ISHIFT(3))
      ALLOCATE(THIS%RBAS(3,3))
      ALLOCATE(THIS%R(3,NAT))
      ALLOCATE(THIS%ATOMID(NAT))
      ALLOCATE(THIS%ISPECIES(NAT))
      ALLOCATE(THIS%XK(3,NKPT))
      ALLOCATE(THIS%WKPT(NKPT))
      ALLOCATE(THIS%DENMAT(LMNXX,LMNXX,NDIM,NAT))
      THIS%INITIALIZED=.TRUE.
      RETURN
      END SUBROUTINE SIMULATION$INIT

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$SETI4A(ID,LEN,VAL)  ! MARK: SIMULATION$SETI4A
!     **************************************************************************
!     ** SET INTEGER VALUE IN SIMULATION MODULE                               **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      INTEGER(4), INTENT(IN) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$SETI4')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$SETI4')
      END IF
      IF(ID.EQ.'LMNX') THEN
        IF(LEN.NE.THIS%NSP) THEN
          CALL ERROR$MSG('LENGTH OF LMNX ARRAY MUST BE EQUAL TO NSP')
          CALL ERROR$CHVAL('LENGTH OF LMNX: ',LEN)
          CALL ERROR$CHVAL('NSP: ',THIS%NSP)
          CALL ERROR$STOP('SIMULATION$SETI4')
        END IF
        THIS%LMNX=VAL
      ELSE IF(ID.EQ.'LNX') THEN
        IF(LEN.NE.THIS%NSP) THEN
          CALL ERROR$MSG('LENGTH OF LNX ARRAY MUST BE EQUAL TO NSP')
          CALL ERROR$CHVAL('LENGTH OF LNX: ',LEN)
          CALL ERROR$CHVAL('NSP: ',THIS%NSP)
          CALL ERROR$STOP('SIMULATION$SETI4')
        END IF
        THIS%LNX=VAL
      ELSE IF(ID.EQ.'LOX') THEN
        IF(LEN.NE.THIS%LNXX*THIS%NSP) THEN
          CALL ERROR$MSG('LENGTH OF LOX ARRAY MUST BE EQUAL TO LNXX*NSP')
          CALL ERROR$CHVAL('LENGTH OF LOX: ',LEN)
          CALL ERROR$CHVAL('LNXX: ',THIS%LNXX)
          CALL ERROR$CHVAL('NSP: ',THIS%NSP)
          CALL ERROR$STOP('SIMULATION$SETI4')
        END IF
        THIS%LOX(:,:)=RESHAPE(VAL,(/THIS%LNXX,THIS%NSP/))
      ELSE IF(ID.EQ.'MAP') THEN
        IF(LEN.NE.THIS%NAT*THIS%LNXX) THEN
          CALL ERROR$MSG('LENGTH OF MAP ARRAY MUST BE EQUAL TO NAT*LNXX')
          CALL ERROR$CHVAL('LENGTH OF MAP: ',LEN)
          CALL ERROR$CHVAL('NAT: ',THIS%NAT)
          CALL ERROR$CHVAL('LNXX: ',THIS%LNXX)
          CALL ERROR$STOP('SIMULATION$SETI4')
        END IF
        THIS%MAP(:,:)=RESHAPE(VAL,(/THIS%NAT,THIS%LNXX/))
      ELSE IF(ID.EQ.'NKDIV') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('LENGTH OF NKDIV ARRAY MUST BE 3')
          CALL ERROR$CHVAL('LENGTH OF NKDIV: ',LEN)
          CALL ERROR$STOP('SIMULATION$SETI4')
        END IF
        THIS%NKDIV=VAL
      ELSE IF(ID.EQ.'ISHIFT') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('LENGTH OF ISHIFT ARRAY MUST BE 3')
          CALL ERROR$CHVAL('LENGTH OF ISHIFT: ',LEN)
          CALL ERROR$STOP('SIMULATION$SETI4')
        END IF
        THIS%ISHIFT=VAL
      ELSE IF(ID.EQ.'ISPECIES') THEN
        IF(LEN.NE.THIS%NAT) THEN
          CALL ERROR$MSG('LENGTH OF ISPECIES ARRAY MUST BE EQUAL TO NAT')
          CALL ERROR$CHVAL('LENGTH OF ISPECIES: ',LEN)
          CALL ERROR$CHVAL('NAT: ',THIS%NAT)
          CALL ERROR$STOP('SIMULATION$SETI4')
        END IF
        THIS%ISPECIES=VAL
      ELSE
        CALL ERROR$MSG('SIMULATION SETI4A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$SETI4')
      END IF 
      RETURN
      END SUBROUTINE SIMULATION$SETI4A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$GETI4A(ID,LEN,VAL)  ! MARK: SIMULATION$GETI4A
!     **************************************************************************
!     ** GET INTEGER VALUE ARRAY FROM SIMULATION MODULE                       **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      INTEGER(4), INTENT(OUT) :: VAL(LEN)
      INTEGER(4) :: ISP
      INTEGER(4) :: LN
      INTEGER(4) :: IPRO
      INTEGER(4) :: IAT
      INTEGER(4) :: L
      INTEGER(4) :: M
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$GETI4A')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$GETI4A')
      END IF
      IF(ID.EQ.'LMNX') THEN
        IF(LEN.NE.THIS%NSP) THEN
          CALL ERROR$MSG('LENGTH OF LMNX ARRAY MUST BE EQUAL TO NSP')
          CALL ERROR$CHVAL('LENGTH OF LMNX: ',LEN)
          CALL ERROR$CHVAL('NSP: ',THIS%NSP)
          CALL ERROR$STOP('SIMULATION$GETI4A')
        END IF
        IF(.NOT.ALLOCATED(THIS%LMNX)) THEN
          ALLOCATE(THIS%LMNX(THIS%NSP))
          THIS%LMNX=0
          DO ISP=1,THIS%NSP
            DO LN=1,THIS%LNX(ISP)
              THIS%LMNX(ISP)=THIS%LMNX(ISP)+2*THIS%LOX(LN,ISP)+1
            ENDDO
          ENDDO
        END IF
        VAL=THIS%LMNX
      ELSE IF(ID.EQ.'LNX') THEN
        IF(LEN.NE.THIS%NSP) THEN
          CALL ERROR$MSG('LENGTH OF LNX ARRAY MUST BE EQUAL TO NSP')
          CALL ERROR$CHVAL('LENGTH OF LNX: ',LEN)
          CALL ERROR$CHVAL('NSP: ',THIS%NSP)
          CALL ERROR$STOP('SIMULATION$GETI4A')
        END IF
        VAL=THIS%LNX
      ELSE IF(ID.EQ.'LOX') THEN
        IF(LEN.NE.THIS%LNXX*THIS%NSP) THEN
          CALL ERROR$MSG('LENGTH OF LOX ARRAY MUST BE EQUAL TO LNXX*NSP')
          CALL ERROR$CHVAL('LENGTH OF LOX: ',LEN)
          CALL ERROR$CHVAL('LNXX: ',THIS%LNXX)
          CALL ERROR$CHVAL('NSP: ',THIS%NSP)
          CALL ERROR$STOP('SIMULATION$GETI4A')
        END IF
        VAL(:)=RESHAPE(THIS%LOX,(/THIS%LNXX*THIS%NSP/))
      ELSE IF(ID.EQ.'MAP') THEN
        IF(LEN.NE.THIS%NAT*THIS%LNXX) THEN
          CALL ERROR$MSG('LENGTH OF MAP ARRAY MUST BE EQUAL TO NAT*LNXX')
          CALL ERROR$CHVAL('LENGTH OF MAP: ',LEN)
          CALL ERROR$CHVAL('NAT: ',THIS%NAT)
          CALL ERROR$CHVAL('LNXX: ',THIS%LNXX)
          CALL ERROR$STOP('SIMULATION$GETI4A')
        END IF
        IF(.NOT.ALLOCATED(THIS%MAP)) THEN
          ALLOCATE(THIS%MAP(THIS%NAT,THIS%LNXX))
          THIS%MAP(:,:)=0
          IPRO=0
          DO IAT=1,THIS%NAT
            ISP=THIS%ISPECIES(IAT)
            DO LN=1,THIS%LNX(ISP)
              THIS%MAP(IAT,LN)=IPRO
              L=THIS%LOX(LN,ISP)
              DO M=1,2*L+1
                IPRO=IPRO+1
              ENDDO ! END M
            ENDDO ! END LN
          ENDDO ! END IAT
        END IF
        VAL(:)=RESHAPE(THIS%MAP,(/THIS%NAT*THIS%LNXX/))
      ELSE IF(ID.EQ.'NKDIV') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('LENGTH OF NKDIV ARRAY MUST BE 3')
          CALL ERROR$CHVAL('LENGTH OF NKDIV: ',LEN)
          CALL ERROR$STOP('SIMULATION$GETI4A')
        END IF
        VAL=THIS%NKDIV
      ELSE IF(ID.EQ.'ISHIFT') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('LENGTH OF ISHIFT ARRAY MUST BE 3')
          CALL ERROR$CHVAL('LENGTH OF ISHIFT: ',LEN)
          CALL ERROR$STOP('SIMULATION$GETI4A')
        END IF
        VAL=THIS%ISHIFT
      ELSE IF(ID.EQ.'ISPECIES') THEN
        IF(LEN.NE.THIS%NAT) THEN
          CALL ERROR$MSG('LENGTH OF ISPECIES ARRAY MUST BE EQUAL TO NAT')
          CALL ERROR$CHVAL('LENGTH OF ISPECIES: ',LEN)
          CALL ERROR$CHVAL('NAT: ',THIS%NAT)
          CALL ERROR$STOP('SIMULATION$GETI4A')
        END IF
        VAL=THIS%ISPECIES
      ELSE
        CALL ERROR$MSG('SIMULATION GETI4A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$GETI4A')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$GETI4A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$SETI4(ID,VAL)  ! MARK: SIMULATION$SETI4
!     **************************************************************************
!     ** SET INTEGER VALUE IN SIMULATION MODULE                               **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$SETI4')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$SETI4')
      END IF
      IF(ID.EQ.'SPACEGROUP') THEN
        THIS%SPACEGROUP=VAL
      ELSE
        CALL ERROR$MSG('SIMULATION SETI4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$SETI4')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$SETI4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$GETI4(ID,VAL)  ! MARK: SIMULATION$GETI4
!     **************************************************************************
!     ** GET INTEGER VALUE FROM SIMULATION MODULE                             **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(OUT) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$GETI4')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$GETI4')
      END IF
      IF(ID.EQ.'NAT') THEN
        VAL=THIS%NAT
      ELSE IF(ID.EQ.'NSP') THEN
        VAL=THIS%NSP
      ELSE IF(ID.EQ.'NKPT') THEN
        VAL=THIS%NKPT
      ELSE IF(ID.EQ.'NSPIN') THEN
        VAL=THIS%NSPIN
      ELSE IF(ID.EQ.'NDIM') THEN
        VAL=THIS%NDIM
      ELSE IF(ID.EQ.'NPRO') THEN
        VAL=THIS%NPRO
      ELSE IF(ID.EQ.'LNXX') THEN
        VAL=THIS%LNXX
      ELSE IF(ID.EQ.'SPACEGROUP') THEN
        VAL=THIS%SPACEGROUP
      ELSE
        CALL ERROR$MSG('SIMULATION GETI4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$GETI4')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$GETI4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$SETR8A(ID,LEN,VAL)  ! MARK: SIMULATION$SETR8A
!     **************************************************************************
!     ** SET REAL VALUE ARRAY IN SIMULATION MODULE                            **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      REAL(8), INTENT(IN) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$SETR8A')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$SETR8A')
      END IF
      IF(ID.EQ.'RBAS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('LENGTH OF RBAS ARRAY MUST BE 9')
          CALL ERROR$CHVAL('LENGTH OF RBAS: ',LEN)
          CALL ERROR$STOP('SIMULATION$SETR8A')
        END IF
        THIS%RBAS=RESHAPE(VAL,(/3,3/))
        THIS%VCELL=THIS%RBAS(1,1)*(THIS%RBAS(2,2)*THIS%RBAS(3,3)-THIS%RBAS(2,3)*THIS%RBAS(3,2)) &
     &         +THIS%RBAS(2,1)*(THIS%RBAS(3,2)*THIS%RBAS(1,3)-THIS%RBAS(3,3)*THIS%RBAS(1,2)) &
     &         +THIS%RBAS(3,1)*(THIS%RBAS(1,2)*THIS%RBAS(2,3)-THIS%RBAS(1,3)*THIS%RBAS(2,2)) 
      ELSE IF(ID.EQ.'R') THEN
        IF(LEN.NE.THIS%NAT*3) THEN
          CALL ERROR$MSG('LENGTH OF R ARRAY MUST BE EQUAL TO NAT*3')
          CALL ERROR$CHVAL('LENGTH OF R: ',LEN)
          CALL ERROR$CHVAL('NAT: ',THIS%NAT)
          CALL ERROR$STOP('SIMULATION$SETR8A')
        END IF
        THIS%R(:,:)=RESHAPE(VAL,(/3,THIS%NAT/))
      ELSE IF(ID.EQ.'XK') THEN
        IF(LEN.NE.THIS%NKPT*3) THEN
          CALL ERROR$MSG('LENGTH OF XK ARRAY MUST BE EQUAL TO NKPT*3')
          CALL ERROR$CHVAL('LENGTH OF XK: ',LEN)
          CALL ERROR$CHVAL('NKPT: ',THIS%NKPT)
          CALL ERROR$STOP('SIMULATION$SETR8A')
        END IF
        THIS%XK(:,:)=RESHAPE(VAL,(/3,THIS%NKPT/))
      ELSE IF(ID.EQ.'WKPT') THEN
        IF(LEN.NE.THIS%NKPT) THEN
          CALL ERROR$MSG('LENGTH OF WKPT ARRAY MUST BE EQUAL TO NKPT')
          CALL ERROR$CHVAL('LENGTH OF WKPT: ',LEN)
          CALL ERROR$CHVAL('NKPT: ',THIS%NKPT)
          CALL ERROR$STOP('SIMULATION$SETR8A')
        END IF
        THIS%WKPT(:)=VAL
      ELSE
        CALL ERROR$MSG('SIMULATION SETR8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$SETR8A')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$SETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$GETR8A(ID,LEN,VAL)  ! MARK: SIMULATION$GETR8A
!     **************************************************************************
!     ** GET REAL VALUE ARRAY FROM SIMULATION MODULE                          **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      REAL(8), INTENT(OUT) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$GETR8A')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$GETR8A')
      END IF
      IF(ID.EQ.'RBAS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('LENGTH OF RBAS ARRAY MUST BE 9')
          CALL ERROR$CHVAL('LENGTH OF RBAS: ',LEN)
          CALL ERROR$STOP('SIMULATION$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS%RBAS,(/9/))
      ELSE IF(ID.EQ.'R') THEN
        IF(LEN.NE.THIS%NAT*3) THEN
          CALL ERROR$MSG('LENGTH OF R ARRAY MUST BE EQUAL TO NAT*3')
          CALL ERROR$CHVAL('LENGTH OF R: ',LEN)
          CALL ERROR$CHVAL('NAT: ',THIS%NAT)
          CALL ERROR$STOP('SIMULATION$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS%R,(/3*THIS%NAT/))
      ELSE IF(ID.EQ.'XK') THEN
        IF(LEN.NE.THIS%NKPT*3) THEN
          CALL ERROR$MSG('LENGTH OF XK ARRAY MUST BE EQUAL TO NKPT*3')
          CALL ERROR$CHVAL('LENGTH OF XK: ',LEN)
          CALL ERROR$CHVAL('NKPT: ',THIS%NKPT)
          CALL ERROR$STOP('SIMULATION$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS%XK,(/3*THIS%NKPT/))
      ELSE IF(ID.EQ.'WKPT') THEN
        IF(LEN.NE.THIS%NKPT) THEN
          CALL ERROR$MSG('LENGTH OF WKPT ARRAY MUST BE EQUAL TO NKPT')
          CALL ERROR$CHVAL('LENGTH OF WKPT: ',LEN)
          CALL ERROR$CHVAL('NKPT: ',THIS%NKPT)
          CALL ERROR$STOP('SIMULATION$GETR8A')
        END IF
        VAL(:)=THIS%WKPT
      ELSE
        CALL ERROR$MSG('SIMULATION GETR8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$GETR8A')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$GETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$SETR8(ID,VAL)  ! MARK: SIMULATION$SETR8
!     **************************************************************************
!     ** SET REAL VALUE IN SIMULATION MODULE                                  **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      REAL(8), INTENT(IN) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$SETR8')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$SETR8')
      END IF
      IF(ID.EQ.'RNTOT') THEN
        THIS%RNTOT=VAL
      ELSE IF(ID.EQ.'NEL') THEN
        THIS%NEL=VAL
      ELSE IF(ID.EQ.'ETOT') THEN
        THIS%ETOT=VAL
      ELSE IF(ID.EQ.'EDFT') THEN
        THIS%EDFT=VAL
      ELSE IF(ID.EQ.'ECORE') THEN
        THIS%ECORE=VAL
      ELSE
        CALL ERROR$MSG('SIMULATION SETR8 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$SETR8')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$SETR8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$GETR8(ID,VAL)  ! MARK: SIMULATION$GETR8
!     **************************************************************************
!     ** GET REAL VALUE FROM SIMULATION MODULE                                **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      REAL(8), INTENT(OUT) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$GETR8')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$GETR8')
      END IF
      IF(ID.EQ.'RNTOT') THEN
        VAL=THIS%RNTOT
      ELSE IF(ID.EQ.'NEL') THEN
        VAL=THIS%NEL
      ELSE IF(ID.EQ.'ETOT') THEN
        ! IF ETOT IS HUGE, IT MEANS IT HAS NOT BEEN CALCULATED YET
        IF(THIS%ETOT.EQ.HUGE(1.D0)) THEN
          ! IF ECORE IS HUGE, IT MEANS IT HAS NOT BEEN CALCULATED YET
          IF(THIS%ECORE.EQ.HUGE(1.D0)) CALL SIMULATION_ECORE
          THIS%ETOT=THIS%EDFT+THIS%ECORE
        END IF
        VAL=THIS%ETOT
      ELSE IF(ID.EQ.'EDFT') THEN
        VAL=THIS%EDFT
      ELSE IF(ID.EQ.'ECORE') THEN
        ! IF ECORE IS HUGE, IT MEANS IT HAS NOT BEEN CALCULATED YET
        IF(THIS%ECORE.EQ.HUGE(1.D0)) CALL SIMULATION_ECORE
        VAL=THIS%ECORE
      ELSE IF(ID.EQ.'VCELL') THEN
        ! IF VCELL IS NEGATIVE HUGE, IT MEANS IT HAS NOT BEEN CALCULATED YET
        IF(THIS%VCELL.EQ.-HUGE(1.D0)) THEN
          CALL ERROR$MSG('VCELL NOT CALCULATED YET, NEEDS RBAS')
          CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
          CALL ERROR$STOP('SIMULATION$GETR8')
        END IF
        VAL=THIS%VCELL  
      ELSE
        CALL ERROR$MSG('SIMULATION GETR8 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$GETR8')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$GETR8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$SETL4(ID,VAL)  ! MARK: SIMULATION$SETL4
!     **************************************************************************
!     ** SET LOGICAL VALUE IN SIMULATION MODULE                               **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      LOGICAL(4), INTENT(IN) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$SETL4')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$SETL4')
      END IF
      IF(ID.EQ.'TINV') THEN
        THIS%TINV=VAL
      ELSE IF(ID.EQ.'TSHIFT') THEN
        THIS%TSHIFT=VAL
      ELSE
        CALL ERROR$MSG('SIMULATION SETL4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$SETL4')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$SETL4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$GETL4(ID,VAL)  ! MARK: SIMULATION$GETL4
!     **************************************************************************
!     ** GET LOGICAL VALUE FROM SIMULATION MODULE                             **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      LOGICAL(4), INTENT(OUT) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$GETL4')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$GETL4')
      END IF
      IF(ID.EQ.'TINV') THEN
        VAL=THIS%TINV
      ELSE IF(ID.EQ.'TSHIFT') THEN
        VAL=THIS%TSHIFT
      ELSE
        CALL ERROR$MSG('SIMULATION GETL4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$GETL4')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$GETL4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$SETCHA(ID,LEN,VAL)  ! MARK: SIMULATION$SETCHA
!     **************************************************************************
!     ** SET CHARACTER VALUE ARRAY IN SIMULATION MODULE                       **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      CHARACTER(*), INTENT(IN) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$SETCHA')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$SETCHA')
      END IF
      IF(ID.EQ.'ATOMID') THEN
        IF(LEN.NE.THIS%NAT) THEN
          CALL ERROR$MSG('LENGTH OF ATOMID ARRAY MUST BE EQUAL TO NAT')
          CALL ERROR$CHVAL('LENGTH OF ATOMID: ',LEN)
          CALL ERROR$CHVAL('NAT: ',THIS%NAT)
          CALL ERROR$STOP('SIMULATION$SETCHA')
        END IF
        THIS%ATOMID=VAL
      ELSE
        CALL ERROR$MSG('SIMULATION SETCHA ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$SETCHA')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$SETCHA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$GETCHA(ID,LEN,VAL)  ! MARK: SIMULATION$GETCHA
!     **************************************************************************
!     ** GET CHARACTER VALUE ARRAY FROM SIMULATION MODULE                     **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: LEN
      CHARACTER(*), INTENT(OUT) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$GETCHA')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$GETCHA')
      END IF
      IF(ID.EQ.'ATOMID') THEN
        IF(LEN.NE.THIS%NAT) THEN
          CALL ERROR$MSG('LENGTH OF ATOMID ARRAY MUST BE EQUAL TO NAT')
          CALL ERROR$CHVAL('LENGTH OF ATOMID: ',LEN)
          CALL ERROR$CHVAL('NAT: ',THIS%NAT)
          CALL ERROR$STOP('SIMULATION$GETCHA')
        END IF
        VAL=THIS%ATOMID
      ELSE
        CALL ERROR$MSG('SIMULATION GETCHA ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$GETCHA')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$GETCHA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$SETCH(ID,VAL)  ! MARK: SIMULATION$SETCH
!     **************************************************************************
!     ** SET CHARACTER VALUE IN SIMULATION MODULE                             **
!     ** CAN SET VALUES BEFORE INITIALIZATION                                 **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      CHARACTER(*), INTENT(IN) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$SETCH')
      END IF
      IF(ID.EQ.'ID') THEN
        THIS%ID=VAL
      ELSE IF(ID.EQ.'FILE') THEN
        THIS%FILE=VAL
      ELSE
        CALL ERROR$MSG('SIMULATION SETCH ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$SETCH')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$SETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION$GETCH(ID,VAL)  ! MARK: SIMULATION$GETCH
!     **************************************************************************
!     ** GET CHARACTER VALUE FROM SIMULATION MODULE                           **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      CHARACTER(*), INTENT(OUT) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION$GETCH')
      END IF
      IF(.NOT.THIS%INITIALIZED) THEN
        CALL ERROR$MSG('SIMULATION NOT INITIALIZED')
        CALL ERROR$CHVAL('SELECTED SIMULATION ID: ',THIS%ID)
        CALL ERROR$STOP('SIMULATION$GETCH')
      END IF
      IF(ID.EQ.'ID') THEN
        VAL=THIS%ID
      ELSE IF(ID.EQ.'FILE') THEN
        VAL=THIS%FILE
      ELSE IF(ID.EQ.'FLAG') THEN
        VAL=THIS%FLAG
      ELSE
        CALL ERROR$MSG('SIMULATION GETCH ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('SIMULATION$GETCH')
      END IF
      RETURN
      END SUBROUTINE SIMULATION$GETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMULATION_ECORE
!     **************************************************************************
!     ** CALCULATE ECORE FOR SELECTED SIMULATION                              **
!     ** REQUIRES SELECTED SIMULATION                                         **
!     **************************************************************************
      USE SIMULATION_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      REAL(8) :: E
      INTEGER(4) :: IAT
      INTEGER(4) :: ISP
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('SIMULATION_ECORE')
      END IF
      THIS%ECORE=0.D0
      CALL XSETUP$SELECT(THIS%ID)
      DO IAT=1,THIS%NAT
        DO ISP=1,THIS%NSP
          CALL XSETUP$GETR8('ECORE',ISP,E)
          THIS%ECORE=THIS%ECORE+E
        ENDDO
      ENDDO
      CALL XSETUP$UNSELECT
      RETURN
      END SUBROUTINE SIMULATION_ECORE
!
!     ==========================================================================
!     ==========================================================================
!     ==                   XSETUP MODULE FUNCTIONS                            ==
!     ==========================================================================
!     ==========================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$SELECT(ID)  ! MARK: XSETUP$SELECT
!     **************************************************************************
!     ** SELECT XSETUP ARRAY                                                  **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: GROUND,EXCITE,SELECTED,THIS
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT SELECT SETUP WHEN ALREADY SELECTED')
        CALL ERROR$STOP('XSETUP$SELECT')
      END IF
      IF(ID.EQ.'GROUND') THEN
        IF(.NOT.ALLOCATED(GROUND)) THEN
          CALL ERROR$MSG('GROUND SETUP ARRAY NOT ALLOCATED')
          CALL ERROR$STOP('XSETUP$SELECT')
        END IF
        THIS=>GROUND
      ELSE IF(ID.EQ.'EXCITE') THEN
        IF(.NOT.ALLOCATED(EXCITE)) THEN
          CALL ERROR$MSG('EXCITE SETUP ARRAY NOT ALLOCATED')
          CALL ERROR$STOP('XSETUP$SELECT')
        END IF
        THIS=>EXCITE
      ELSE
        CALL ERROR$MSG('XSETUP SELECT ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('XSETUP$SELECT')
      END IF
      SELECTED=.TRUE.
      RETURN
      END SUBROUTINE XSETUP$SELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$UNSELECT()  ! MARK: XSETUP$UNSELECT
!     **************************************************************************
!     ** UNSELECT XSETUP ARRAY                                                **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED, THIS
      IMPLICIT NONE
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT UNSELECT SETUP ARRAY WHEN NOT SELECTED')
        CALL ERROR$STOP('XSETUP$UNSELECT')
      END IF
      SELECTED=.FALSE.
      NULLIFY(THIS)
      RETURN
      END SUBROUTINE XSETUP$UNSELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$ISELECT(I)  ! MARK: XSETUP$ISELECT
!     **************************************************************************
!     ** SELECT XSETUP ARRAY BY INTEGER INDEX                                 **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,GROUND,EXCITE,THIS
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: I
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT SELECT SETUP WHEN ALREADY SELECTED')
        CALL ERROR$STOP('XSETUP$ISELECT')
      END IF
      IF(I.EQ.1) THEN
        IF(.NOT.ALLOCATED(GROUND)) THEN
          CALL ERROR$MSG('GROUND SETUP ARRAY NOT ALLOCATED')
          CALL ERROR$STOP('XSETUP$ISELECT')
        END IF
        THIS=>GROUND
      ELSE IF(I.EQ.2) THEN
        IF(.NOT.ALLOCATED(EXCITE)) THEN
          CALL ERROR$MSG('EXCITE SETUP ARRAY NOT ALLOCATED')
          CALL ERROR$STOP('XSETUP$ISELECT')
        END IF
        THIS=>EXCITE
      ELSE
        CALL ERROR$MSG('XSETUP ISELECT INDEX NOT RECOGNIZED')
        CALL ERROR$CHVAL('INDEX: ',I)
        CALL ERROR$STOP('XSETUP$ISELECT')
      END IF
      SELECTED=.TRUE.
      RETURN
      END SUBROUTINE XSETUP$ISELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$INIT(NSP)  ! MARK: XSETUP$INIT
!     **************************************************************************
!     ** INITIALIZE XSETUP ARRAY                                              **  
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: THIS,SELECTED
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NSP
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$INIT')
      END IF
      IF(NSP.LE.0) THEN
        CALL ERROR$MSG('INVALID NUMBER OF SPECIES')
        CALL ERROR$CHVAL('NUMBER OF SPECIES: ',NSP)
        CALL ERROR$STOP('XSETUP$INIT')
      END IF
      IF(ALLOCATED(THIS)) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT INITIALIZE SETUP ARRAY WHEN ALREADY ALLOCATED')
        CALL ERROR$STOP('XSETUP$INIT')
      END IF
      ALLOCATE(THIS(NSP))
      RETURN
      END SUBROUTINE XSETUP$INIT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$NEW(ISP,GID,ECORE,LNX,NBATOM)  ! MARK: XSETUP$NEW
!     **************************************************************************
!     ** ADD NEW SETUP TO SETUP ARRAY                                         **
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: ISP
      INTEGER(4), INTENT(IN) :: GID
      REAL(8), INTENT(IN) :: ECORE
      INTEGER(4), INTENT(IN) :: LNX
      INTEGER(4), INTENT(IN) :: NBATOM
      INTEGER(4) :: NR
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$INIT')
      END IF
      IF(ISP.LE.0.OR.ISP.GT.SIZE(THIS)) THEN
        CALL ERROR$MSG('INVALID SPECIES INDEX')
        CALL ERROR$CHVAL('SPECIES INDEX: ',ISP)
        CALL ERROR$STOP('XSETUP$NEW')
      END IF
      CALL RADIAL$GETI4(GID,'NR',NR)
      THIS(ISP)%GID=GID
      THIS(ISP)%ECORE=ECORE
      THIS(ISP)%LNX=LNX
      ALLOCATE(THIS%PSPHI(NR,LNX))
      ALLOCATE(THIS%AEPHI(NR,LNX))
      THIS(ISP)%NBATOM=NBATOM
      ALLOCATE(THIS%LATOM(NBATOM))
      ALLOCATE(THIS%AEPSI(NR,NBATOM))
      RETURN
      END SUBROUTINE XSETUP$NEW
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$GETI4A(ID,ISP,LEN,VAL)  ! MARK: XSETUP$GETI4A
!     **************************************************************************
!     ** GET INTEGER VALUE ARRAY FROM SELECTED SETUP                          **
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: ISP
      INTEGER(4), INTENT(IN) :: LEN
      INTEGER(4), INTENT(OUT) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$GETI4A')
      END IF
      IF(ISP.LE.0.OR.ISP.GT.SIZE(THIS)) THEN
        CALL ERROR$MSG('INVALID SPECIES INDEX')
        CALL ERROR$CHVAL('SPECIES INDEX: ',ISP)
        CALL ERROR$STOP('XSETUP$GETI4A')
      END IF
      IF(ID.EQ.'LATOM') THEN
        IF(LEN.NE.THIS(ISP)%NBATOM) THEN
          CALL ERROR$MSG('LENGTH OF LATOM ARRAY MUST BE EQUAL TO NBATOM')
          CALL ERROR$CHVAL('LENGTH OF LATOM: ',LEN)
          CALL ERROR$CHVAL('NBATOM: ',THIS(ISP)%NBATOM)
          CALL ERROR$STOP('XSETUP$GETI4A')
        END IF
        VAL(:)=THIS(ISP)%LATOM(:)
      ELSE
        CALL ERROR$MSG('XSETUP GETI4A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('XSETUP$GETI4A')
      END IF
      RETURN
      END SUBROUTINE XSETUP$GETI4A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$SETI4A(ID,ISP,LEN,VAL)  ! MARK: XSETUP$SETI4A
!     **************************************************************************
!     ** SET INTEGER VALUE ARRAY IN SELECTED SETUP                            **
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: ISP
      INTEGER(4), INTENT(IN) :: LEN
      INTEGER(4), INTENT(IN) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$SETI4A')
      END IF
      IF(ISP.LE.0.OR.ISP.GT.SIZE(THIS)) THEN
        CALL ERROR$MSG('INVALID SPECIES INDEX')
        CALL ERROR$CHVAL('SPECIES INDEX: ',ISP)
        CALL ERROR$STOP('XSETUP$SETI4A')
      END IF
      IF(ID.EQ.'LATOM') THEN
        IF(LEN.NE.THIS(ISP)%NBATOM) THEN
          CALL ERROR$MSG('LENGTH OF LATOM ARRAY MUST BE EQUAL TO NBATOM')
          CALL ERROR$CHVAL('LENGTH OF LATOM: ',LEN)
          CALL ERROR$CHVAL('NBATOM: ',THIS(ISP)%NBATOM)
          CALL ERROR$STOP('XSETUP$SETI4A')
        END IF
        THIS(ISP)%LATOM(:)=VAL(:)
      ELSE
        CALL ERROR$MSG('XSETUP SETI4A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('XSETUP$SETI4A')
      END IF
      RETURN
      END SUBROUTINE XSETUP$SETI4A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$GETI4(ID,ISP,VAL)  ! MARK: XSETUP$GETI4
!     **************************************************************************
!     ** GET INTEGER VALUE FROM SELECTED SETUP                                **
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: ISP
      INTEGER(4), INTENT(OUT) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$GETI4')
      END IF
      IF(ISP.LE.0.OR.ISP.GT.SIZE(THIS)) THEN
        CALL ERROR$MSG('INVALID SPECIES INDEX')
        CALL ERROR$CHVAL('SPECIES INDEX: ',ISP)
        CALL ERROR$STOP('XSETUP$GETI4')
      END IF
      IF(ID.EQ.'GID') THEN
        VAL=THIS(ISP)%GID
      ELSE IF(ID.EQ.'NBATOM') THEN
        VAL=THIS(ISP)%NBATOM
      ELSE
        CALL ERROR$MSG('XSETUP GETI4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('XSETUP$GETI4')
      END IF
      RETURN
      END SUBROUTINE XSETUP$GETI4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$GETR8A(ID,ISP,LEN,VAL)  ! MARK: XSETUP$GETR8A
!     **************************************************************************
!     ** GET REAL VALUE ARRAY FROM SELECTED SETUP                             **
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: ISP
      INTEGER(4), INTENT(IN) :: LEN
      REAL(8), INTENT(OUT) :: VAL(LEN)
      REAL(8) :: NR
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$GETR8A')
      END IF
      IF(ISP.LE.0.OR.ISP.GT.SIZE(THIS)) THEN
        CALL ERROR$MSG('INVALID SPECIES INDEX')
        CALL ERROR$CHVAL('SPECIES INDEX: ',ISP)
        CALL ERROR$STOP('XSETUP$GETR8A')
      END IF
      CALL RADIAL$GETI4(THIS(ISP)%GID,'NR',NR)
      IF(ID.EQ.'PSPHI') THEN
        IF(LEN.NE.NR*THIS(ISP)%LNX) THEN
          CALL ERROR$MSG('LENGTH MUST BE EQUAL TO NR*LNX')
          CALL ERROR$CHVAL('LENGTH: ',LEN)
          CALL ERROR$CHVAL('NR*LNX: ',NR*THIS(ISP)%LNX)
          CALL ERROR$STOP('XSETUP$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS(ISP)%PSPHI,(/LEN/))
      ELSE IF(ID.EQ.'AEPHI') THEN
        IF(LEN.NE.NR*THIS(ISP)%LNX) THEN
          CALL ERROR$MSG('LENGTH MUST BE EQUAL TO NR*LNX')
          CALL ERROR$CHVAL('LENGTH: ',LEN)
          CALL ERROR$CHVAL('NR*LNX: ',NR*THIS(ISP)%LNX)
          CALL ERROR$STOP('XSETUP$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS(ISP)%AEPHI,(/LEN/))
      ELSE IF(ID.EQ.'AEPSI') THEN
        IF(LEN.NE.NR*THIS(ISP)%NBATOM) THEN
          CALL ERROR$MSG('LENGTH MUST BE EQUAL TO NR*NBATOM')
          CALL ERROR$CHVAL('LENGTH: ',LEN)
          CALL ERROR$CHVAL('NR*NBATOM: ',NR*THIS(ISP)%NBATOM)
          CALL ERROR$STOP('XSETUP$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS(ISP)%AEPSI,(/LEN/))
      ELSE
        CALL ERROR$MSG('XSETUP GETR8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('XSETUP$GETR8A')
      END IF
      RETURN
      END SUBROUTINE XSETUP$GETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$SETR8A(ID,ISP,LEN,VAL)  ! MARK: XSETUP$SETR8A
!     **************************************************************************
!     ** SET REAL VALUE ARRAY IN SELECTED SETUP                               **
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: ISP
      INTEGER(4), INTENT(IN) :: LEN
      REAL(8), INTENT(IN) :: VAL(LEN)
      REAL(8) :: NR
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$SETR8A')
      END IF
      IF(ISP.LE.0.OR.ISP.GT.SIZE(THIS)) THEN
        CALL ERROR$MSG('INVALID SPECIES INDEX')
        CALL ERROR$CHVAL('SPECIES INDEX: ',ISP)
        CALL ERROR$STOP('XSETUP$SETR8A')
      END IF
      CALL RADIAL$GETI4(THIS(ISP)%GID,'NR',NR)
      IF(ID.EQ.'PSPHI') THEN
        IF(LEN.NE.NR*THIS(ISP)%LNX) THEN
          CALL ERROR$MSG('LENGTH MUST BE EQUAL TO NR*LNX')
          CALL ERROR$CHVAL('LENGTH: ',LEN)
          CALL ERROR$CHVAL('NR*LNX: ',NR*THIS(ISP)%LNX)
          CALL ERROR$STOP('XSETUP$SETR8A')
        END IF
        THIS(ISP)%PSPHI=RESHAPE(VAL,(/NR,THIS(ISP)%LNX/))
      ELSE IF(ID.EQ.'AEPHI') THEN
        IF(LEN.NE.NR*THIS(ISP)%LNX) THEN
          CALL ERROR$MSG('LENGTH MUST BE EQUAL TO NR*LNX')
          CALL ERROR$CHVAL('LENGTH: ',LEN)
          CALL ERROR$CHVAL('NR*LNX: ',NR*THIS(ISP)%LNX)
          CALL ERROR$STOP('XSETUP$SETR8A')
        END IF
        THIS(ISP)%AEPHI=RESHAPE(VAL,(/NR,THIS(ISP)%LNX/))
      ELSE IF(ID.EQ.'AEPSI') THEN
        IF(LEN.NE.NR*THIS(ISP)%NBATOM) THEN
          CALL ERROR$MSG('LENGTH MUST BE EQUAL TO NR*NBATOM')
          CALL ERROR$CHVAL('LENGTH: ',LEN)
          CALL ERROR$CHVAL('NR*NBATOM: ',NR*THIS(ISP)%NBATOM)
          CALL ERROR$STOP('XSETUP$SETR8A')
        END IF
        THIS(ISP)%AEPSI=RESHAPE(VAL,(/NR,THIS(ISP)%NBATOM/))
      ELSE
        CALL ERROR$MSG('XSETUP SETR8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('XSETUP$SETR8A')
      END IF
      RETURN
      END SUBROUTINE XSETUP$SETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$GETR8(ID,ISP,VAL)  ! MARK: XSETUP$GETR8
!     **************************************************************************
!     ** GET REAL VALUE FROM SELECTED SETUP                                   **
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: ISP
      REAL(8), INTENT(OUT) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$GETR8')
      END IF
      IF(ISP.LE.0.OR.ISP.GT.SIZE(THIS)) THEN
        CALL ERROR$MSG('INVALID SPECIES INDEX')
        CALL ERROR$CHVAL('SPECIES INDEX: ',ISP)
        CALL ERROR$STOP('XSETUP$GETR8')
      END IF
      IF(ID.EQ.'ECORE') THEN
        VAL=THIS(ISP)%ECORE
      ELSE
        CALL ERROR$MSG('XSETUP GETR8 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('XSETUP$GETR8')
      END IF
      RETURN
      END SUBROUTINE XSETUP$GETR8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XSETUP$SETR8(ID,ISP,VAL)  ! MARK: XSETUP$SETR8
!     **************************************************************************
!     ** SET REAL VALUE IN SELECTED SETUP                                     **
!     ** REQUIRES SELECTED SETUP ARRAY                                        **
!     **************************************************************************
      USE XSETUP_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: ISP
      REAL(8), INTENT(IN) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('XSETUP$SETR8')
      END IF
      IF(ISP.LE.0.OR.ISP.GT.SIZE(THIS)) THEN
        CALL ERROR$MSG('INVALID SPECIES INDEX')
        CALL ERROR$CHVAL('SPECIES INDEX: ',ISP)
        CALL ERROR$STOP('XSETUP$SETR8')
      END IF
      IF(ID.EQ.'ECORE') THEN
        THIS(ISP)%ECORE=VAL
      ELSE
        CALL ERROR$MSG('XSETUP SETR8 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('XSETUP$SETR8')
      END IF
      RETURN
      END SUBROUTINE XSETUP$SETR8
!
!     ==========================================================================
!     ==========================================================================
!     ==                    KSMAP MODULE FUNCTIONS                            ==
!     ==========================================================================
!     ==========================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KSMAP$INIT(NKPTG_,NSPING_)
!     **************************************************************************
!     ** INITIALIZE KSMAP FOR TASKS                                           **
!     ** KSMAP(NKPTG,NSPING) IS THE MAP OF KPOINTS AND SPINS TO TASKS         **
!     ** RTASK IS THE TASK RESPONSIBLE FOR READING AND WRITING                **
!     **************************************************************************
      USE KSMAP_MODULE, ONLY: NKPTG,NSPING,RTASK,KSMAP,INITIALIZED
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NKPTG_
      INTEGER(4), INTENT(IN) :: NSPING_
      INTEGER(4) :: NTASKS,THISTASK
      INTEGER(4) :: IKPT,ISPIN,IND
      IF(INITIALIZED) THEN
        ! IF ALREADY INITIALIZED, CHECK IF SAME VALUES
        IF(NPKTG_.NE.NKPTG.OR.NSPING_.NE.NSPING) THEN
          CALL ERROR$MSG('KSMAP ALREADY INITIALIZED WITH DIFFERENT VALUES')
          CALL ERROR$I4VAL('OLD NKPTG: ',NKPTG)
          CALL ERROR$I4VAL('OLD NSPING: ',NSPING)
          CALL ERROR$I4VAL('NEW NKPTG: ',NKPTG_)
          CALL ERROR$I4VAL('NEW NSPING: ',NSPING)
          CALL ERROR$STOP('KSMAP$INIT')
        END IF
        RETURN
      END IF
                          CALL TRACE$PUSH('KSMAP$INIT')
      NSPING=NSPING_
      NKPTG=NKPTG_
      ALLOCATE(KSMAP(NKPTG,NSPING))
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      ! ONLY ONE TASK: 
      !   BOTH READ TASK AND RESPONSIBLE FOR EVERY KPOINT
      IF(NTASKS.EQ.1) THEN
        RTASK=1
        KSMAP=1
      ! MORE THAN ONE TASK:
      !   FIRST IS READ TASK, OTHERS ARE RESPONSIBLE FOR DIFFERENT KPOINTS
      ELSE
        RTASK=1
        DO IKPT=1,NKPTG
          DO ISPIN=1,NSPING
            IND=NSPING*(IKPT-1)+ISPIN
            IND=MOD(IND-1,NTASKS-1)+2
            KSMAP(IKPT,ISPIN)=IND
          ENDDO
        ENDDO
      END IF
      INITIALIZED=.TRUE.
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE KSMAP$INIT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KSMAP$REPORT(NFIL)
!     **************************************************************************
!     ** REPORT KSMAP TO FILE                                                 **
!     **************************************************************************
      USE KSMAP_MODULE, ONLY: NKPTG,NSPING,RTASK,KSMAP,INITIALIZED
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NFIL
      INTEGER(4), PARAMETER :: PERLINE=14
      INTEGER(4) :: IKPT,ISPIN,IND
      INTEGER(4) :: BEGINKPT,ENDKPT
      INTEGER(4) :: NTASKS,THISTASK
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      IF(THISTASK.NE.RTASK) RETURN
                          CALL TRACE$PUSH('KSMAP$REPORT')
      IF(.NOT.INITIALIZED) THEN
        CALL ERROR$MSG('KSMAP NOT INITIALIZED')
        CALL ERROR$STOP('KSMAP$REPORT')
      END IF        
      WRITE(NFIL,FMT='(A8,I5)')'RW TASK:',RTASK
      BEGINKPT=1
      DO WHILE(.TRUE.)
        ! KPT IN THIS LINE (HALF FOR NSPIN=2)
        ENDKPT=MIN(BEGINKPT+PERLINE/NSPING-1,NKPTG)
        IF(NSPING.EQ.1) THEN
          WRITE(NFIL,FMT='(A8,*(I4,"|"))')'KPOINT:',(IKPT, IKPT=BEGINKPT,ENDKPT)
        ELSE
          WRITE(NFIL,FMT='(A8,*(I9,"|"))')'KPOINT:',(IKPT, IKPT=BEGINKPT,ENDKPT)
        END IF
        WRITE(NFIL,FMT='(A8,*(I4,"|"))')'SPIN:',(MOD(IND-1,NSPING)+1, IND=BEGINKPT*NSPING-1,ENDKPT*NSPING)
        WRITE(NFIL,FMT='(A8,*(I4,"|"))')'KSMAP:',(KSMAP(IKPT,:), IKPT=BEGINKPT,ENDKPT)
        IF(ENDKPT.NE.NKPTG) THEN
          BEGINKPT=ENDKPT+1
        ELSE
          EXIT
        END IF
      ENDDO
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE KSMAP$REPORT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KSMAP$READTASK(RTASK_)  ! MARK: KSMAP$READTASK
!     **************************************************************************
!     ** GET READ TASK                                                        **
!     ** REQUIRES KSMAP INITIALIZED                                           **
!     **************************************************************************
      USE KSMAP_MODULE, ONLY: RTASK,INITIALIZED
      IMPLICIT NONE
      INTEGER(4), INTENT(OUT) :: RTASK_
      IF(.NOT.INITIALIZED) THEN
        CALL ERROR$MSG('KSMAP NOT INITIALIZED')
        CALL ERROR$STOP('KSMAP$READTASK')
      END IF
      RTASK_=RTASK
      RETURN
      END SUBROUTINE KSMAP$READTASK
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KSMAP$WORKTASK(IKPT,ISPIN,TASK_)  ! MARK: KSMAP$WORKTASK
!     **************************************************************************  
!     ** GET WORK TASK FOR KPOINT AND SPIN                                    **
!     ** REQUIRES KSMAP INITIALIZED                                           **
!     **************************************************************************
      USE KSMAP_MODULE, ONLY: KSMAP,INITIALIZED,NKPTG,NSPING
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(OUT) :: TASK_
      IF(.NOT.INITIALIZED) THEN
        CALL ERROR$MSG('KSMAP NOT INITIALIZED')
        CALL ERROR$STOP('KSMAP$WORKTASK')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$CHVAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('KSMAP$WORKTASK')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPING) THEN
        CALL ERROR$MSG('INVALID SPIN INDEX')
        CALL ERROR$CHVAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('KSMAP$WORKTASK')
      END IF
      TASK_=KSMAP(IKPT,ISPIN)
      RETURN
      END SUBROUTINE KSMAP$WORKTASK
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KSMAP$WORKTASKSKPT(IKPT,TASKS_)  ! MARK: KSMAP$WORKTASKSKPT
!     **************************************************************************  
!     ** GET WORK TASKS FOR KPOINT                                            **
!     ** REQUIRES KSMAP INITIALIZED                                           **
!     **************************************************************************
      USE KSMAP_MODULE, ONLY: KSMAP,INITIALIZED,NKPTG,NSPING
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(OUT) :: TASKS_(NSPING)
      IF(.NOT.INITIALIZED) THEN
        CALL ERROR$MSG('KSMAP NOT INITIALIZED')
        CALL ERROR$STOP('KSMAP$WORKTASK')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$CHVAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('KSMAP$WORKTASK')
      END IF
      TASKS_(:)=KSMAP(IKPT,:)
      RETURN
      END SUBROUTINE KSMAP$WORKTASKSKPT
!
!     ==========================================================================
!     ==========================================================================
!     ==                    STATE MODULE FUNCTIONS                            ==
!     ==========================================================================
!     ==========================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$INIT(NKPT_,NSPIN_,NDIM_,NPRO_)  ! MARK: STATE$INIT
!     **************************************************************************
!     ** INITIALIZE STATE MODULE                                              **
!     **************************************************************************
      USE STATE_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NKPT_
      INTEGER(4), INTENT(IN) :: NSPIN_
      INTEGER(4), INTENT(IN) :: NDIM_
      INTEGER(4), INTENT(IN) :: NPRO_
                          CALL TRACE$PUSH('STATE$INIT')
      IF(INITIALIZED) THEN
        ! IF ALREADY INITIALIZED, CHECK IF SAME VALUES
        IF(NKPT_.NE.NKPTG.OR.NSPIN_.NE.NSPING.OR.NDIM_.NE.NDIM.OR.NPRO_.NE.NPRO) THEN
          CALL ERROR$MSG('STATE ALREADY INITIALIZED WITH DIFFERENT VALUES')
          CALL ERROR$I4VAL('OLD NKPTG: ',NKPTG)
          CALL ERROR$I4VAL('OLD NSPING: ',NSPING)
          CALL ERROR$I4VAL('OLD NDIM: ',NDIMG)
          CALL ERROR$I4VAL('OLD NPRO: ',NPRO)
          CALL ERROR$I4VAL('NEW NKPTG: ',NKPT_)
          CALL ERROR$I4VAL('NEW NSPING: ',NSPIN_)
          CALL ERROR$I4VAL('NEW NDIM: ',NDIM_)
          CALL ERROR$I4VAL('NEW NPRO: ',NPRO_)
          CALL ERROR$STOP('STATE$INIT')
        END IF
        RETURN
      END IF
      NKPTG=NKPT_
      NSPING=NSPIN_
      NDIM=NDIM_
      NPRO=NPRO_
      ALLOCATE(GROUND(NKPTG,NSPING))
      ALLOCATE(EXCITE(NKPTG,NSPING))
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE STATE$INIT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$SELECT(ID)  ! MARK: STATE$SELECT
!     **************************************************************************
!     ** SELECT STATE ARRAY                                                   **
!     **************************************************************************
      USE STATE_MODULE, ONLY: GROUND,EXCITE,THIS,SELECTED
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT SELECT STATE ARRAY WHEN ALREADY SELECTED')
        CALL ERROR$STOP('STATE$SELECT')
      END IF
      IF(ID.EQ.'GROUND') THEN
        IF(.NOT.ALLOCATED(GROUND)) THEN
          CALL ERROR$MSG('GROUND STATE ARRAY NOT ALLOCATED')
          CALL ERROR$STOP('STATE$SELECT')
        END IF
        THIS=>GROUND
      ELSE IF(ID.EQ.'EXCITE') THEN
        IF(.NOT.ALLOCATED(EXCITE)) THEN
          CALL ERROR$MSG('EXCITE STATE ARRAY NOT ALLOCATED')
          CALL ERROR$STOP('STATE$SELECT')
        END IF
        THIS=>EXCITE
      ELSE
        CALL ERROR$MSG('STATE SELECT ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('STATE$SELECT')
      END IF
      SELECTED=.TRUE.
      RETURN
      END SUBROUTINE STATE$SELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$UNSELECT()  ! MARK: STATE$UNSELECT
!     **************************************************************************
!     ** UNSELECT STATE ARRAY                                                 **
!     **************************************************************************
      USE STATE_MODULE, ONLY: THIS,SELECTED
      IMPLICIT NONE
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT UNSELECT STATE ARRAY WHEN NOT SELECTED')
        CALL ERROR$STOP('STATE$UNSELECT')
      END IF
      SELECTED=.FALSE.
      DEALLOCATE(THIS)
      RETURN
      END SUBROUTINE STATE$UNSELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$ISELECT(I)  ! MARK: STATE$ISELECT
!     **************************************************************************
!     ** SELECT STATE ARRAY BY INDEX                                          **
!     **************************************************************************
      USE STATE_MODULE, ONLY: GROUND,EXCITE,THIS,SELECTED
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: I
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT SELECT STATE ARRAY WHEN ALREADY SELECTED')
        CALL ERROR$STOP('STATE$ISELECT')
      END IF
      IF(I.EQ.1) THEN
        IF(.NOT.ALLOCATED(GROUND)) THEN
          CALL ERROR$MSG('GROUND STATE ARRAY NOT ALLOCATED')
          CALL ERROR$STOP('STATE$ISELECT')
        END IF
        THIS=>GROUND
      ELSE IF(I.EQ.2) THEN
        IF(.NOT.ALLOCATED(EXCITE)) THEN
          CALL ERROR$MSG('EXCITE STATE ARRAY NOT ALLOCATED')
          CALL ERROR$STOP('STATE$ISELECT')
        END IF
        THIS=>EXCITE
      ELSE
        CALL ERROR$MSG('STATE ISELECT INDEX NOT RECOGNIZED')
        CALL ERROR$I4VAL('INDEX: ',I)
        CALL ERROR$STOP('STATE$ISELECT')
      END IF
      SELECTED=.TRUE.
      RETURN
      END SUBROUTINE STATE$ISELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$NEW(IKPT,ISPIN,NB)  ! MARK: STATE$NEW
!     **************************************************************************
!     ** CREATE NEW STATE IN ARRAY                                            **
!     ** REQUIRES SELECTED STATE ARRAY                                        **
!     **************************************************************************
      USE STATE_MODULE, ONLY: THIS,SELECTED,NKPTG,NSPING
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(IN) :: NB
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO STATE ARRAY SELECTED')
        CALL ERROR$STOP('STATE$NEW')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$I4VAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('STATE$NEW')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPING) THEN
        CALL ERROR$MSG('INVALID SPIN INDEX')
        CALL ERROR$I4VAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('STATE$NEW')
      END IF
      IF(NB.LE.0) THEN
        CALL ERROR$MSG('NUMBER OF BANDS MUST BE POSITIVE')
        CALL ERROR$I4VAL('NUMBER OF BANDS: ',NB)
        CALL ERROR$STOP('STATE$NEW')
      END IF
      IF(THIS(IKPT,ISPIN)%NB.NE.-HUGE(1)) THEN
        CALL ERROR$MSG('STATE ALREADY EXISTS')
        CALL ERROR$I4VAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$I4VAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('STATE$NEW')
      END IF
      ! INITIALIZE STATE
      THIS(IKPT,ISPIN)%NB=NB
      ALLOCATE(THIS(IKPT,ISPIN)%EIG(NB))
      ALLOCATE(THIS(IKPT,ISPIN)%OCC(NB))
      ! TODO: PROJ IS NOT ALLOCATED HERE AS IT MIGHT NOT BE USED WITH A RESTART FILE
      RETURN
      END SUBROUTINE STATE$NEW
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$GETI4(ID,IKPT,ISPIN,VAL)  ! MARK: STATE$GETI4
!     **************************************************************************
!     ** GET INTEGER VALUE FROM SELECTED STATE                                **
!     ** REQUIRES SELECTED STATE ARRAY                                        **
!     **************************************************************************
      USE STATE_MODULE, ONLY: THIS,SELECTED,NKPTG,NSPING
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(OUT) :: VAL
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO STATE ARRAY SELECTED')
        CALL ERROR$STOP('STATE$GETI4')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$I4VAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('STATE$GETI4')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPING) THEN
        CALL ERROR$MSG('INVALID SPIN INDEX')
        CALL ERROR$I4VAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('STATE$GETI4')
      END IF
      IF(ID.EQ.'NB') THEN
        VAL=THIS(IKPT,ISPIN)%NB
      ELSE IF(ID.EQ.'NOCC') THEN
        ! TODO: TRIGGER CALCULATION OF NOCC HERE
        CALL ERROR$MSG('NOCC NOT IMPLEMENTED YET')
        CALL ERROR$STOP('STATE$GETI4')
        VAL=THIS(IKPT,ISPIN)%NOCC
      ELSE
        CALL ERROR$MSG('STATE GETI4 ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('STATE$GETI4')
      END IF
      RETURN
      END SUBROUTINE STATE$GETI4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$SETR8A(ID,IKPT,ISPIN,LEN,VAL)  ! MARK: STATE$SETR8A
!     **************************************************************************
!     ** SET REAL VALUE ARRAY IN SELECTED STATE                               **
!     ** REQUIRES SELECTED STATE ARRAY                                        **
!     **************************************************************************
      USE STATE_MODULE, ONLY: THIS,SELECTED,NKPTG,NSPING
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(IN) :: LEN
      REAL(8), INTENT(IN) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO STATE ARRAY SELECTED')
        CALL ERROR$STOP('STATE$SETR8A')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$I4VAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('STATE$SETR8A')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPING) THEN
        CALL ERROR$MSG('INVALID SPIN INDEX')
        CALL ERROR$I4VAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('STATE$SETR8A')
      END IF
      IF(ID.EQ.'EIG') THEN
        IF(LEN.NE.THIS(IKPT,ISPIN)%NB) THEN
          CALL ERROR$MSG('LEN MUST BE EQUAL TO NUMBER OF BANDS')
          CALL ERROR$I4VAL('LEN: ',LEN)
          CALL ERROR$I4VAL('NB: ',THIS(IKPT,ISPIN)%NB)
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('STATE$SETR8A')
        END IF
        THIS(IKPT,ISPIN)%EIG(:)=VAL(:)
      ELSE IF(ID.EQ.'OCC') THEN
        IF(LEN.NE.THIS(IKPT,ISPIN)%NB) THEN
          CALL ERROR$MSG('LENGTH MUST BE EQUAL TO NUMBER OF BANDS')
          CALL ERROR$I4VAL('LEN: ',LEN)
          CALL ERROR$I4VAL('NB: ',THIS(IKPT,ISPIN)%NB)
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('STATE$SETR8A')
        END IF
        THIS(IKPT,ISPIN)%OCC=VAL(:)
      ELSE
        CALL ERROR$MSG('STATE SETR8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('STATE$SETR8A')
      END IF
      RETURN
      END SUBROUTINE STATE$SETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$GETR8A(ID,IKPT,ISPIN,LEN,VAL)  ! MARK: STATE$GETR8A
!     **************************************************************************
!     ** GET REAL VALUE ARRAY FROM SELECTED STATE                             **
!     ** REQUIRES SELECTED STATE ARRAY                                        **
!     **************************************************************************
      USE STATE_MODULE, ONLY: THIS,SELECTED,NKPTG,NSPING
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(IN) :: LEN
      REAL(8), INTENT(OUT) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO STATE ARRAY SELECTED')
        CALL ERROR$STOP('STATE$GETR8A')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$I4VAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('STATE$GETR8A')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPING) THEN
        CALL ERROR$MSG('INVALID SPIN INDEX')
        CALL ERROR$I4VAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('STATE$GETR8A')
      END IF
      IF(ID.EQ.'EIG') THEN
        IF(LEN.NE.THIS(IKPT,ISPIN)%NB) THEN
          CALL ERROR$MSG('LEN MUST BE EQUAL TO NUMBER OF BANDS')
          CALL ERROR$I4VAL('LEN: ',LEN)
          CALL ERROR$I4VAL('NB: ',THIS(IKPT,ISPIN)%NB)
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('STATE$GETR8A')
        END IF
        VAL(:)=THIS(IKPT,ISPIN)%EIG(:)
      ELSE IF(ID.EQ.'OCC') THEN
        IF(LEN.NE.THIS(IKPT,ISPIN)%NB) THEN
          CALL ERROR$MSG('LENGTH MUST BE EQUAL TO NUMBER OF BANDS')
          CALL ERROR$I4VAL('LEN: ',LEN)
          CALL ERROR$I4VAL('NB: ',THIS(IKPT,ISPIN)%NB)
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('STATE$GETR8A')
        END IF
        VAL(:)=THIS(IKPT,ISPIN)%OCC(:)
      ELSE
        CALL ERROR$MSG('STATE GETR8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('STATE$GETR8A')
      END IF
      RETURN
      END SUBROUTINE STATE$GETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$SETC8A(ID,IKPT,ISPIN,LEN,VAL)  ! MARK: STATE$SETC8A
!     **************************************************************************
!     ** SET COMPLEX VALUE ARRAY IN SELECTED STATE                            **
!     ** REQUIRES SELECTED STATE ARRAY                                        **
!     **************************************************************************
      USE STATE_MODULE, ONLY: THIS,SELECTED,NKPTG,NSPING,NDIM,NPRO
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(IN) :: LEN
      COMPLEX(8), INTENT(IN) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO STATE ARRAY SELECTED')
        CALL ERROR$STOP('STATE$SETC8A')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$I4VAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('STATE$SETC8A')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPING) THEN
        CALL ERROR$MSG('INVALID SPIN INDEX')
        CALL ERROR$I4VAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('STATE$SETC8A')
      END IF
      IF(ID.EQ.'PROJ') THEN
        IF(LEN.NE.NDIM*THIS(IKPT,ISPIN)%NB*NPRO) THEN
          CALL ERROR$MSG('LEN MUST BE EQUAL TO NDIM*NB*NPRO')
          CALL ERROR$I4VAL('LEN: ',LEN)
          CALL ERROR$I4VAL('NDIM*NB*NPRO: ',NDIM*THIS(IKPT,ISPIN)%NB*NPRO)
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('STATE$SETC8A')
        END IF
        IF(.NOT.ALLOCATED(THIS(IKPT,ISPIN)%PROJ)) THEN
          ALLOCATE(THIS(IKPT,ISPIN)%PROJ(NDIM,THIS(IKPT,ISPIN)%NB,NPRO))
        END IF
        THIS(IKPT,ISPIN)%PROJ=RESHAPE(VAL,(/NDIM,THIS(IKPT,ISPIN)%NB,NPRO/))
      ELSE
        CALL ERROR$MSG('STATE SETC8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('STATE$SETC8A')
      END IF
      RETURN
      END SUBROUTINE STATE$SETC8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STATE$GETC8A(ID,IKPT,ISPIN,LEN,VAL)  ! MARK: STATE$GETC8A
!     **************************************************************************
!     ** GET COMPLEX VALUE ARRAY FROM SELECTED STATE                          **
!     ** REQUIRES SELECTED STATE ARRAY                                        **
!     **************************************************************************
      USE STATE_MODULE, ONLY: THIS,SELECTED,NKPTG,NSPING,NDIM,NPRO
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      INTEGER(4), INTENT(IN) :: LEN
      COMPLEX(8), INTENT(OUT) :: VAL(LEN)
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO STATE ARRAY SELECTED')
        CALL ERROR$STOP('STATE$GETC8A')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$I4VAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('STATE$GETC8A')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPING) THEN
        CALL ERROR$MSG('INVALID SPIN INDEX')
        CALL ERROR$I4VAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('STATE$GETC8A')
      END IF
      IF(ID.EQ.'PROJ') THEN
        IF(LEN.NE.NDIM*THIS(IKPT,ISPIN)%NB*NPRO) THEN
          CALL ERROR$MSG('LEN MUST BE EQUAL TO NDIM*NB*NPRO')
          CALL ERROR$I4VAL('LEN: ',LEN)
          CALL ERROR$I4VAL('NDIM*NB*NPRO: ',NDIM*THIS(IKPT,ISPIN)%NB*NPRO)
          CALL ERROR$CHVAL('ID: ',ID)
          CALL ERROR$STOP('STATE$GETC8A')
        END IF
        IF(.NOT.ALLOCATED(THIS(IKPT,ISPIN)%PROJ)) THEN
          CALL ERROR$MSG('PROJ NOT ALLOCATED')
          CALL ERROR$I4VAL('IKPT: ',IKPT)
          CALL ERROR$I4VAL('ISPIN: ',ISPIN)
          CALL ERROR$STOP('STATE$GETC8A')
        END IF
        VAL(:)=RESHAPE(THIS(IKPT,ISPIN)%PROJ,(/NDIM*THIS(IKPT,ISPIN)%NB*NPRO/))
      ELSE
        CALL ERROR$MSG('STATE GETC8A ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID: ',ID)
        CALL ERROR$STOP('STATE$GETC8A')
      END IF
      RETURN
      END SUBROUTINE STATE$GETC8A
!
!     ==========================================================================
!     ==========================================================================
!     ==                    OVERLAP MODULE FUNCTIONS                          ==
!     ==========================================================================
!     ==========================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE OVERLAP$INIT(NKPT_,NSPIN_)  ! MARK: OVERLAP$INIT
!     **************************************************************************
!     ** INITIALIZE OVERLAP MODULE                                            **
!     **************************************************************************
      USE OVERLAP_MODULE, ONLY: NKPTG,NSPING,OVLARR,INITIALIZED
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NKPT_
      INTEGER(4), INTENT(IN) :: NSPIN_
                          CALL TRACE$PUSH('OVERLAP$INIT')
      IF(INITIALIZED) THEN
        CALL ERROR$MSG('OVERLAP ALREADY INITIALIZED')
        CALL ERROR$STOP('OVERLAP$INIT')
      END IF
      NKPTG=NKPT_
      NSPING=NSPIN_
      ALLOCATE(OVLARR(NKPTG,NSPING))
      INITIALIZED=.TRUE.
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE OVERLAP$INIT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE OVERLAP$SELECT(IKPT,ISPIN)  ! MARK: OVERLAP$SELECT
!     **************************************************************************
!     ** SELECT OVERLAP FOR KPOINT AND SPIN                                   **
!     **************************************************************************
      USE OVERLAP_MODULE, ONLY: OVLARR,THIS,SELECTED,NKPTG,NSPING
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: IKPT
      INTEGER(4), INTENT(IN) :: ISPIN
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT SELECT OVERLAP WHEN ALREADY SELECTED')
        CALL ERROR$STOP('OVERLAP$SELECT')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPTG) THEN
        CALL ERROR$MSG('INVALID KPOINT INDEX')
        CALL ERROR$I4VAL('KPOINT INDEX: ',IKPT)
        CALL ERROR$STOP('OVERLAP$SELECT')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPING) THEN
        CALL ERROR$MSG('INVALID SPIN INDEX')
        CALL ERROR$I4VAL('SPIN INDEX: ',ISPIN)
        CALL ERROR$STOP('OVERLAP$SELECT')
      END IF
      THIS=>OVLARR(IKPT,ISPIN)
      SELECTED=.TRUE.
      RETURN
      END SUBROUTINE OVERLAP$SELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE OVERLAP$UNSELECT()  ! MARK: OVERLAP$UNSELECT
!     **************************************************************************
!     ** UNSELECT OVERLAP                                                    **
!     **************************************************************************
      USE OVERLAP_MODULE, ONLY: THIS,SELECTED
      IMPLICIT NONE
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT UNSELECT OVERLAP WHEN NOT SELECTED')
        CALL ERROR$STOP('OVERLAP$UNSELECT')
      END IF
      SELECTED=.FALSE.
      NULLIFY(THIS)
      RETURN
      END SUBROUTINE OVERLAP$UNSELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XRAY$READ  ! MARK: XRAY$READ
!     **************************************************************************
!     ** READ SIMULATION FILES IN PARALLEL AND CALCULATE PLANE WAVE OVERLAP   **
!     ** ON THE FLY.                                                          **
!     **************************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      ! CHANGING THIS WILL LIKELY BREAK THINGS
      INTEGER(4), PARAMETER :: NSIM=2
      LOGICAL(4), PARAMETER :: TPR
      INTEGER(4) :: NFIL(2)
      CHARACTER(6) :: ID
      INTEGER(4) :: NAT,NSP,NKPT,NSPIN,NDIM,LNXX
      INTEGER(4) :: NPRO(NSIM)
      CHARACTER(6) :: FLAG
      INTEGER(4), ALLOCATABLE :: LNX(:)  ! (NSP)
      INTEGER(4), ALLOCATABLE :: LOX(:,:)  ! (LNXX,NSP)
      INTEGER(4) :: NKDIV(3)
      INTEGER(4) :: ISHIFT(3)
      REAL(8) :: RNTOT
      REAL(8) :: NEL
      REAL(8) :: EDFT
      INTEGER(4) :: ILOGICAL
      LOGICAL(4) :: TINV
      INTEGER(4) :: SPACEGROUP
      LOGICAL(4) :: TSHIFT
      REAL(8), ALLOCATABLE :: R(:,:)  ! (3,NAT)
      CHARACTER(16), ALLOCATABLE :: ATOMID(:)  ! (NAT)
      INTEGER(4), ALLOCATABLE :: ISPECIES(:)  ! (NAT)
      REAL(8) :: RBAS(3,3)
      REAL(8) :: VCELL
      INTEGER(4) :: GID
      CHARACTER(8) :: GRIDTYPE
      INTEGER(4) :: NR
      REAL(8) :: DEX
      REAL(8) :: R1
      REAL(8) :: ECORE
      REAL(8), ALLOCATABLE :: PSPHI(:,:)  ! (NR,LNX)
      REAL(8), ALLOCATABLE :: AEPHI(:,:)  ! (NR,LNX)
      INTEGER(4) :: NBATOM
      INTEGER(4), ALLOCATABLE :: LATOM(:)  ! (NBATOM)
      REAL(8), ALLOCATABLE :: AEPSI(:,:)  ! (NR,NBATOM)
      INTEGER(4) :: NBX
      INTEGER(4) :: NKPT_
      INTEGER(4) :: NSPIN_
      REAL(8), ALLOCATABLE :: OCC(:,:,:)  ! (NMX,NKPT_,NSPIN_)
      REAL(8), ALLOCATABLE :: XK(:,:)  ! (3,NKPT_)
      REAL(8), ALLOCATABLE :: WKPT(:)  ! (NKPT_)
      CHARACTER(8) :: KEY(NSIM)
      INTEGER(4) :: NGG(NSIM)
      INTEGER(4) :: NDIM_(NSIM)
      INTEGER(4) :: NB(NSIM)
      LOGICAL(4) :: TSUPER(NSIM)

      INTEGER(4), ALLOCATABLE :: IGVEC(:,:,:) ! (NSIM,3,NGG)
      REAL(8) :: XK_(NSIM,3)
      COMPLEX(8), ALLOCATABLE :: PSIK1(:,:,:) ! (NGG,NDIM,NB)
      COMPLEX(8), ALLOCATABLE :: PSIK2(:,:,:) ! (NGG,NDIM,NB)
      COMPLEX(8), ALLOCATABLE :: PROJ1(:,:,:) ! (NDIM,NB,NPRO)
      COMPLEX(8), ALLOCATABLE :: PROJ2(:,:,:) ! (NDIM,NB,NPRO)
      REAL(8), ALLOCATABLE :: EIG1(:) ! (NB)
      REAL(8), ALLOCATABLE :: EIG2(:) ! (NB)

      COMPLEX(8), ALLOCATABLE :: PW(:,:) ! (NB2,NB1)

      INTEGER(4) :: ISIM
      INTEGER(4) :: ISP
      INTEGER(4) :: IKPT
      INTEGER(4) :: ISPIN
      INTEGER(4) :: IB1
      INTEGER(4) :: IB2
      INTEGER(4) :: NTASKS,THISTASK,RTASK,WTASK
      INTEGER(4), ALLOCATABLE :: WTASKSKPT(:) ! (NSPING)

                          CALL TRACE$PUSH('XRAY$READ')
                          CALL TIMING$CLOCKON('XRAY$READ')
      CALL MPE$QUERY('~',NTASKS,THISTASK)

      DO ISIM=1,NSIM
        CALL SIMULATION$ISELECT(ISIM)
        ! ERROR: ID AND FILE MUST HAVE BEEN SET BEFORE
        CALL SIMULATION$GETCH('ID',ID)
        IF(THISTASK.EQ.1) THEN
          IF(TPR) CALL TRACE$PASS('READING GENERAL QUANTITIES')
          ! ERROR: FILEHANDLING FOR SIMULATION FILES MISSING
          CALL FILEHANDLER$UNIT(ID,NFIL(ISIM))
          REWIND(NFIL(ISIM))
          ! ======================================================================
          ! == READ GENERAL QUANTATIES                                          ==
          ! ======================================================================
          ! NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG
          READ(NFIL(ISIM))NAT,NSP,NKPT,NSPIN,NDIM,NPRO(ISIM),LNXX,FLAG
        END IF
        CALL MPE$BROADCAST('~',1,NAT)
        CALL MPE$BROADCAST('~',1,NSP)
        CALL MPE$BROADCAST('~',1,NKPT)
        CALL MPE$BROADCAST('~',1,NSPIN)
        CALL MPE$BROADCAST('~',1,NDIM)
        CALL MPE$BROADCAST('~',1,NPRO(ISIM))
        CALL MPE$BROADCAST('~',1,LNXX)
        CALL MPE$BROADCAST('~',1,FLAG)
        CALL SIMULATION$INIT(NAT,NSP,NKPT,NSPIN,NDIM,NPRO(ISIM),LNXX,FLAG)
        ! DISTRIBUTE WORKLOAD FOR K-POINT AND SPIN LOOPS WHILE READING
        CALL KSMAP$INIT(NKPT,NSPIN)
        CALL KSMAP$READTASK(RTASK)

        ALLOCATE(LNX(NSP))
        ALLOCATE(LOX(LNXX,NSP))
        IF(THISTASK.EQ.RTASK) THEN
          ! LNX(NSP),LOX(LNXX,NSP)
          READ(NFIL(ISIM))LNX,LOX
          ! NKDIV(3),ISHIFT(3),RNTOT,NEL,EDFT,ILOGICAL
          READ(NFIL(ISIM))NKDIV,ISHIFT,RNTOT,NEL,EDFT,ILOGICAL
          TINV=.FALSE.
          IF(ILOGICAL.EQ.1) TINV=.TRUE.
          ! SPACEGROUP,TSHIFT
          READ(NFIL(ISIM))SPACEGROUP,ILOGICAL
          TSHIFT=.FALSE.
          IF(ILOGICAL.EQ.1) TSHIFT=.TRUE.
        END IF
        CALL MPE$BROADCAST('~',RTASK,LNX)
        CALL MPE$BROADCAST('~',RTASK,LOX)
        CALL MPE$BROADCAST('~',RTASK,NKDIV)
        CALL MPE$BROADCAST('~',RTASK,ISHIFT)
        CALL MPE$BROADCAST('~',RTASK,RNTOT)
        CALL MPE$BROADCAST('~',RTASK,NEL)
        CALL MPE$BROADCAST('~',RTASK,EDFT)
        CALL MPE$BROADCAST('~',RTASK,TINV)
        CALL MPE$BROADCAST('~',RTASK,SPACEGROUP)
        CALL MPE$BROADCAST('~',RTASK,TSHIFT)
        CALL SIMULATION$SETI4A('LNX',NSP,LNX)
        CALL SIMULATION$SETI4A('LOX',LNXX*NSP,LOX)
        CALL SIMULATION$SETI4A('NKDIV',3,NKDIV)
        CALL SIMULATION$SETI4A('ISHIFT',3,ISHIFT)
        CALL SIMULATION$SETR8('RNTOT',RNTOT)
        CALL SIMULATION$SETR8('NEL',NEL)
        CALL SIMULATION$SETR8('EDFT',EDFT)
        CALL SIMULATION$SETL4('TINV',TINV)
        CALL SIMULATION$SETI4('SPACEGROUP',SPACEGROUP)
        CALL SIMULATION$SETL4('TSHIFT',TSHIFT)
        ! ======================================================================
        ! ==  READ ATOMIC STRUCTURE                                           ==
        ! ======================================================================
        ALLOCATE(R(3,NAT))
        ALLOCATE(ATOMID(NAT))
        ALLOCATE(ISPECIES(NAT))
        IF(THISTASK.EQ.RTASK) THEN
          IF(TPR) CALL TRACE$PASS('READING ATOMIC STRUCTURE')
          ! RBAS(3,3),R(3,NAT),ATOMID(NAT),ISPECIES(NAT)
          READ(NFIL(ISIM))RBAS,R,ATOMID,ISPECIES
        END IF
        CALL MPE$BROADCAST('~',RTASK,RBAS)
        CALL MPE$BROADCAST('~',RTASK,R)
        CALL MPE$BROADCAST('~',RTASK,ATOMID)
        CALL MPE$BROADCAST('~',RTASK,ISPECIES)
        CALL SIMULATION$SETR8A('RBAS',9,RBAS)
        CALL SIMULATION$SETR8A('R',3*NAT,R)
        CALL SIMULATION$SETCH('ATOMID',NAT,ATOMID)
        CALL SIMULATION$SETI4A('ISPECIES',NAT,ISPECIES)
        CALL SIMULATION$GETR8('VCELL',VCELL)
        ! ======================================================================
        ! ==  ELEMENT SPECIFIC QUANTITIES                                     ==
        ! ======================================================================
        CALL XSETUP$ISELECT(ISIM)
        CALL XSETUP$INIT(NSP)
        DO ISP=1,NSP
          IF(THISTASK.EQ.RTASK) THEN
            IF(TPR) CALL TRACE$PASS('READING ELEMENT SPECIFIC QUANTITIES')
            ! GRIDTYPE,NR,DEX,R1
            READ(NFIL(ISIM))GRIDTYPE,NR,DEX,R1
            ! ECORE
            READ(NFIL(ISIM))ECORE
          END IF
          CALL MPE$BROADCAST('~',RTASK,GRIDTYPE)
          CALL MPE$BROADCAST('~',RTASK,NR)
          CALL MPE$BROADCAST('~',RTASK,DEX)
          CALL MPE$BROADCAST('~',RTASK,R1)
          CALL MPE$BROADCAST('~',RTASK,ECORE)
          CALL RADIAL$NEW(GRIDTYPE,GID)
          CALL RADIAL$SETI4(GID,'NR',NR)
          CALL RADIAL$SETR8(GID,'DEX',DEX)
          CALL RADIAL$SETR8(GID,'R1',R1)
          ALLOCATE(PSPHI(NR,LNX(ISP)))
          ALLOCATE(AEPHI(NR,LNX(ISP)))
          IF(THISTASK.EQ.RTASK) THEN
            ! PSPHI(NR,LNX(ISP))
            READ(NFIL(ISIM))PSPHI
            ! AEPHI(NR,LNX(ISP))
            READ(NFIL(ISIM))AEPHI
            ! NBATOM
            READ(NFIL(ISIM))NBATOM
          END IF
          CALL MPE$BROADCAST('~',RTASK,PSPHI)
          CALL MPE$BROADCAST('~',RTASK,AEPHI)
          CALL MPE$BROADCAST('~',RTASK,NBATOM)
          ALLOCATE(LATOM(NBATOM))
          ALLOCATE(AEPSI(NR,NBATOM))
          IF(THISTASK.EQ.RTASK) THEN
            ! LATOM(NBATOM)
            READ(NFIL(ISIM))LATOM
            ! AEPSI(NR,NBATOM)
            READ(NFIL(ISIM))AEPSI
          END IF
          CALL MPE$BROADCAST('~',RTASK,LATOM)
          CALL MPE$BROADCAST('~',RTASK,AEPSI)
          CALL XSETUP$NEW(ISP,GID,ECORE,LNX(ISP),NBATOM)
          CALL XSETUP$SETR8A('PSPHI',ISP,NR*LNX(ISP),PSPHI)
          CALL XSETUP$SETR8A('AEPHI',ISP,NR*LNX(ISP),AEPHI)
          CALL XSETUP$SETI4A('LATOM',ISP,NBATOM,LATOM)
          CALL XSETUP$SETR8A('AEPSI',ISP,NR*NBATOM,AEPSI)
          DEALLOCATE(PSPHI)
          DEALLOCATE(AEPHI)
          DEALLOCATE(LATOM)
          DEALLOCATE(AEPSI)
        ENDDO ! END ISP
        ! ======================================================================
        ! ==  OCCUPATIONS AND K-POINTS AND THEIR WEIGHTS                      ==
        ! ======================================================================
        IF(THISTASK.EQ.RTASK) THEN
          IF(TPR) CALL TRACE$PASS('READING K-POINTS AND OCCUPATIONS')
          ! NBX,NKPT_,NSPIN_
          READ(NFIL(ISIM))NBX,NKPT_,NSPIN_
          IF(NKPT_.NE.NKPT.OR.NSPIN_.NE.NSPIN) THEN
            CALL ERROR$MSG('NKPT OR NSPIN NOT CONSISTENT WITH INITIALIZATION')
            CALL ERROR$I4VAL('NKPT: ',NKPT_)
            CALL ERROR$I4VAL('NSPIN: ',NSPIN_)
            CALL ERROR$I4VAL('INITIALIZED NKPT: ',NKPT)
            CALL ERROR$I4VAL('INITIALIZED NSPIN: ',NSPIN)
            CALL ERROR$STOP('XRAY$READ')
          END IF
        END IF
        CALL MPE$BROADCAST('~',RTASK,NBX)
        CALL MPE$BROADCAST('~',RTASK,NKPT)
        CALL MPE$BROADCAST('~',RTASK,NSPIN)

        ALLOCATE(OCC(NBX,NKPT,NSPIN))
        ALLOCATE(XK(3,NKPT))
        ALLOCATE(WKPT(NKPT))
        ! OCC(NBX,NKPT,NSPIN), XK(3,NKPT),WKPT(NKPT)
        IF(THISTASK.EQ.RTASK) THEN
          READ(NFIL(ISIM))OCC,XK,WKPT
        END IF
        CALL MPE$BROADCAST('~',RTASK,OCC)
        CALL MPE$BROADCAST('~',RTASK,XK)
        CALL MPE$BROADCAST('~',RTASK,WKPT)
        CALL SIMULATION$SETR8A('XK',3*NKPT,XK)
        CALL SIMULATION$SETR8A('WKPT',NKPT,WKPT)
        DEALLOCATE(XK)
        DEALLOCATE(WKPT)

        CALL STATE$INIT(NKPT,NSPIN,NDIM,NPRO(ISIM))
        CALL STATE$ISELECT(ISIM)
        DO IKPT=1,NKPT
          DO ISPIN=1,NSPIN
            CALL KSMAP$WORKTASK(IKPT,ISPIN,WTASK)
            IF(THISTASK.NE.WTASK.AND.THISTASK.NE.RTASK) CYCLE
            CALL STATE$NEW(IKPT,ISPIN,NBX)
            CALL STATE$SETR8A('OCC',IKPT,ISPIN,NBX,OCC(:,IKPT,ISPIN))
          ENDDO ! END ISPIN
        ENDDO ! END IKPT
        DEALLOCATE(OCC)
        DEALLOCATE(LNX)
        DEALLOCATE(LOX)
        DEALLOCATE(R)
        DEALLOCATE(ATOMID)
        DEALLOCATE(ISPECIES)
        CALL XSETUP$UNSELECT
        CALL SIMULATION$UNSELECT
        CALL STATE$UNSELECT
      ENDDO ! END ISIM

      ! ========================================================================
      ! ==  WAVE FUNCTIONS AND PROJECTIONS                                    ==
      ! ========================================================================
      CALL OVERLAP$INIT(NKPT,NSPIN)
      ALLOCATE(WTASKSKPT(NSPIN))
      DO IKPT=1,NKPT
        CALL KSMAP$WORKTASKSKPT(IKPT,WTASKSKPT)
        ! CYCLE EARLY IF THIS TASK IS NOT WORKING ON THIS KPOINT
        IF(ALL(THISTASK.NE.WTASKSKPT(:)).AND.THISTASK.NE.RTASK) CYCLE

        IF(THISTASK.EQ.RTASK) THEN
          IF(TPR) CALL TRACE$I4VAL(' READ/SEND IKPT: ',IKPT)
          ! READ GENERAL INFORMATION ABOUT K POINT
          DO ISIM=1,NSIM
            ! KEY,NGG,NDIM,NB,TSUPER
            READ(NFIL(ISIM))KEY(ISIM),NGG(ISIM),NDIM_(ISIM),NB(ISIM),ILOGICAL
            TSUPER(ISIM)=.FALSE.
            IF(ILOGICAL.EQ.1) TSUPER(ISIM)=.TRUE.
            IF(KEY(ISIM).NE.'PSI') THEN
              CALL ERROR$MSG('KEY NOT "PSI"')
              CALL ERROR$MSG('FILE IS CORRUPTED')
              CALL ERROR$CHVAL('KEY: ',KEY(ISIM))
              CALL ERROR$I4VAL('IKPT: ',IKPT)
              CALL ERROR$STOP('XRAY$READ')
            END IF
          ENDDO ! END ISIM
          ! TODO: CHECK IF NGG IS THE SAME? IS THAT REQUIRED?
          ! READ K POINT AND G VECTORS 
          ALLOCATE(IGVEC(NSIM,3,NGG(1)))
          DO ISIM=1,NSIM
            ! XK_(3),IGVEC(ISIM,3,NGG)
            READ(NFIL(ISIM))XK_(ISIM,:),IGVEC(ISIM,:,:)
          ENDDO ! END ISIM
          CALL TESTKPT!(NSIM,KEY,NGG,NDIM_,TSUPER,XK_,IGVEC)
          DEALLOCATE(IGVEC)
        END IF
        DO ISPIN=1,NSPIN
          CALL KSMAP$WTASK(IKPT,ISPIN,WTASK)
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,NGG)
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,NDIM_)
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,NB)
        ENDDO ! END ISPIN
        ! WARNING: REQUIRES SUPER WAVE FUNCTIONS TO BE UNRAVELLED
        ALLOCATE(PSIK1(NGG(1),NDIM_(1),NB(1)))
        ALLOCATE(PSIK2(NGG(2),NDIM_(2),NB(2)))
        ALLOCATE(PROJ1(NDIM_(1),NB(1),NPRO(1)))
        ALLOCATE(PROJ2(NDIM_(2),NB(2),NPRO(2)))
        ALLOCATE(EIG1(NB(1)))
        ALLOCATE(EIG2(NB(2)))
        ALLOCATE(PW(NB(2),NB(1)))

        DO ISPIN=1,NSPIN
          CALL KSMAP$WTASK(IKPT,ISPIN,WTASK)
          IF(THISTASK.NE.WTASK.AND.THISTASK.NE.RTASK) CYCLE

          IF(THISTASK.EQ.RTASK) THEN
            ! READ PLANE WAVE BASIS
            READ(NFIL(1))PSIK1
            READ(NFIL(2))PSIK2
            ! READ PROJECTIONS
            READ(NFIL(1))PROJ1
            READ(NFIL(2))PROJ2
            ! READ EIGENVALUES
            READ(NFIL(1))EIG1
            READ(NFIL(2))EIG2
          END IF
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,PSIK1)
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,PSIK2)
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,PROJ1)
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,PROJ2)
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,EIG1)
          CALL MPE$SENDRECEIVE('~',RTASK,WTASK,EIG2)
          
          ! TODO: BREAKS IF NBX != NB(1) OR NB(2)
          CALL STATE$ISELECT(1)
          CALL STATE$SETR8A('EIG',IKPT,ISPIN,NB(1),EIG1)
          CALL STATE$SETC8A('PROJ',IKPT,ISPIN,NDIM_(1)*NB(1)*NPRO(1),PROJ1)
          CALL STATE$UNSELECT
          CALL STATE$ISELECT(2)
          CALL STATE$SETR8A('EIG',IKPT,ISPIN,NB(2),EIG2)
          CALL STATE$SETC8A('PROJ',IKPT,ISPIN,NDIM_(2)*NB(2)*NPRO(2),PROJ2)
          CALL STATE$UNSELECT

          CALL TIMING$CLOCKON('XRAY$READ_SCALARPRODUCT')
          IF(THISTASK.EQ.WTASK) THEN
            CALL TRACE$I4VAL(' CALCULATING IKPT: ',IKPT)
            DO IB2=1,NB(2)
              DO IB1=1,NB(1)
                ! NO DIM LOOP AS NDIM=1
                ! SCALARPRODUCT (SUM OVER G VECTORS)
                ! PW(I,J)=<PSI1(J)|PSI2(I)>
                ! WARNING: CHECK IF SCALARPRODUCT IS CORRECT ALSO WITH CONJG
                !          SHOULD BE THE CASE AS CALL OF ZGEMM IS DONE WITH 'C' OPTION
                ! WARNING: CHECK WILL NOT WORK IF ONLY USING REAL WAVE FUNCTIONS AT GAMMA POINT
                CALL LIB$SCALARPRODUCTC8(.FALSE.,NGG(1),1,PSIK1(:,1,IB1),1,PSIK2(:,1,IB2),PW(IB2,IB1))
              ENDDO ! END IB1
            ENDDO ! END IB2
            PW=PW*VCELL
          END IF
          CALL TIMING$CLOCKOFF('XRAY$READ_SCALARPRODUCT')
          ! MISSING: DETECTION OF OCCUPIED BANDS BASED ON FERMI ENERGY
        ENDDO ! END ISPIN
        DEALLOCATE(PW)
        DEALLOCATE(PSIK1)
        DEALLOCATE(PSIK2)
        DEALLOCATE(PROJ1)
        DEALLOCATE(PROJ2)
        DEALLOCATE(EIG1)
        DEALLOCATE(EIG2)
      ENDDO ! END IKPT
                          CALL TIMING$CLOCKOFF('XRAY$READ')
                          CALL TRACE$POP
      RETURN

      CONTAINS 
        SUBROUTINE TESTKPT!(NSIM,KEY,NGG,NDIM,TSUPER,XK,IGVEC)
        ! INTEGER(4), INTENT(IN) :: NSIM
        ! CHARACTER(8), INTENT(IN) :: KEY(NSIM)
        ! INTEGER(4), INTENT(IN) :: NGG(NSIM)
        ! INTEGER(4), INTENT(IN) :: NDIM(NSIM)
        ! LOGICAL(4), INTENT(IN) :: TSUPER(NSIM)
        ! REAL(8), INTENT(IN) :: XK(NSIM,3)
        ! INTEGER(4), INTENT(IN) :: IGVEC(NSIM,3,NGG(1))
        IF(KEY(1).NE.KEY(2)) THEN
          CALL ERROR$MSG('KEYS NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$CHVAL('KEY1',KEY(1))
          CALL ERROR$CHVAL('KEY2',KEY(2))
          CALL ERROR$STOP('RIXS$READ')
        END IF
        IF(NGG(1).NE.NGG(2)) THEN
          CALL ERROR$MSG('NGG NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('NGG1',NGG(1))
          CALL ERROR$I4VAL('NGG2',NGG(2))
          CALL ERROR$STOP('RIXS$READ')
        END IF
        IF(NDIM_(1).NE.NDIM_(2)) THEN
          CALL ERROR$MSG('NDIM NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('NDIM1',NDIM(1))
          CALL ERROR$I4VAL('NDIM2',NDIM(2))
          CALL ERROR$STOP('RIXS$READ')
        END IF
        ! TODO: CHECK NDIM AGAINST PREVIOUSLY READ ONE
        IF(TSUPER(1).NEQV.TSUPER(2)) THEN
          CALL ERROR$MSG('TSUPER NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$L4VAL('TSUPER1',TSUPER(1))
          CALL ERROR$L4VAL('TSUPER2',TSUPER(2))
          CALL ERROR$STOP('RIXS$READ')
        END IF
        ! CHECK IF XK SAME IN BOTH SIMULATIONS
        IF(SUM(ABS(XK_(1,:)-XK_(2,:)))>TOL) THEN
          CALL ERROR$MSG('XK NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$R8VAL('XK1',XK(1,:))
          CALL ERROR$R8VAL('XK2',XK(2,:))
          CALL ERROR$STOP('RIXS$READ')
        END IF
        ! CHECK IF IGVEC SAME IN BOTH SIMULATIONS
        IF(ANY(IGVEC(1,:,:).NE.IGVEC(2,:,:))) THEN
          CALL ERROR$MSG('IGVEC NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$STOP('RIXS$READ')
        END IF
        END SUBROUTINE TESTKPT
      END SUBROUTINE XRAY$READ


      ! TODO: ATOMMAPPING
      ! TODO: DATACONSISTENCY
