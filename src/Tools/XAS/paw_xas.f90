!........1.........2.........3.........4.........5.........6.........7.........8
      MODULE XAS_MODULE
      TYPE SETUP_TYPE
        INTEGER(4) :: GID ! GRID ID FOR RADIAL MESH
        REAL(8), ALLOCATABLE :: PSPHI(:,:) ! PSEUDO PARTIAL WAVE
        REAL(8), ALLOCATABLE :: AEPHI(:,:) ! AE PARTIAL WAVE
      END TYPE SETUP_TYPE

      TYPE STATE_TYPE
        INTEGER(4)          :: NB ! #(STATES)
        REAL(8), POINTER    :: EIG(:) ! (NB) EIGENVALUES
        REAL(8), POINTER    :: OCC(:)
        COMPLEX(8), POINTER :: PROJ(:,:,:) ! (NDIM,NB,NPRO) PROJECTIONS
      END TYPE STATE_TYPE

      TYPE OVERLAP_TYPE
        COMPLEX(8), POINTER :: PW(:,:) ! (NB1,NB2) PLANE WAVE OVERLAP
        REAL(8), POINTER    :: S(:,:) ! (NPRO1,NPRO2) ATOMIC OVERLAP
      END TYPE OVERLAP_TYPE

!     GENERAL SETTINGS OF A SIMULATION 
      TYPE SIMULATION_TYPE
        CHARACTER(11) :: ID ! SIMULATION IDENTIFIER (GROUNDSTATE,EXCITESTATE)
        CHARACTER(256) :: FILE ! FILENAME
        INTEGER(4) :: NAT ! #(ATOMS)
        INTEGER(4) :: NSP ! #(SETUPS)
        INTEGER(4) :: NKPT ! #(KPOINTS)
        INTEGER(4) :: NSPIN ! #(SPINS)
        INTEGER(4) :: NDIM ! #(DIMENSIONS)
        INTEGER(4) :: NPRO ! #(PROJECTIONS)
        INTEGER(4) :: LNXX
        CHARACTER(6) :: FLAG
        INTEGER(4), ALLOCATABLE :: LNX(:) ! (NSP)
        INTEGER(4), ALLOCATABLE :: LOX(:,:) ! (LNXX,NSP)
        LOGICAL(4) :: TINV ! INVERSION SYMMETRY
        INTEGER(4) :: NKDIV(3) ! K-POINT DIVISIONS
        INTEGER(4) :: ISHIFT(3) ! K-POINT SHIFTS
        REAL(8) :: RNTOT
        REAL(8) :: NEL
        INTEGER(4) :: SPACEGROUP
        LOGICAL(4) :: TSHIFT
        REAL(8) :: RBAS(3,3) ! BASIS VECTORS
        REAL(8), ALLOCATABLE :: R(:,:) ! (3,NAT) ATOMIC POSITIONS
        CHARACTER(16), ALLOCATABLE :: ATOMID(:) ! (NAT) ATOM IDENTIFIERS
        INTEGER(4), ALLOCATABLE :: ISPECIES(:) ! (NAT) SPECIES INDEX
        REAL(8), ALLOCATABLE :: XK(:,:) ! (3,NKPT) K-POINTS IN REL. COORD.
        REAL(8), ALLOCATABLE :: WKPT(:) ! (NKPT) K-POINT WEIGHTS
        TYPE(SETUP_TYPE), ALLOCATABLE :: SETUP(:) ! (NSP) ARRAY OF SETUPS
        TYPE(STATE_TYPE), ALLOCATABLE :: STATEARR(:,:) ! (NKPT,NSPIN) ARRAY OF STATES
        TYPE(STATE_TYPE), POINTER :: STATE ! CURRENT STATE
      END TYPE SIMULATION_TYPE

! WARNING: CHANGING THIS NUMBER MIGHT BREAK SOME CODE
!     INDEX 1 IS GROUNDSTATE AND 2 EXCITED STATE
      INTEGER(4), PARAMETER :: NSIM=2  ! #(SIMULATIONS)
      LOGICAL,SAVE :: SELECTED=.FALSE.
      TYPE(SIMULATION_TYPE), TARGET :: SIM(NSIM) ! (NSIM) ARRAY OF SIMULATIONS
      TYPE(SIMULATION_TYPE), POINTER :: THIS ! CURRENT SIMULATION
      TYPE(OVERLAP_TYPE), ALLOCATABLE, TARGET :: OVERLAPARR(:,:) ! (NKPT,NSPIN)
      TYPE(OVERLAP_TYPE), POINTER :: OVERLAP ! CURRENT OVERLAP
      END MODULE XAS_MODULE

      MODULE XASCNTL_MODULE  ! MARK: XASCNTL_MODULE
        USE LINKEDLIST_MODULE, ONLY: LL_TYPE
        TYPE(LL_TYPE) :: LL_CNTL
        SAVE
      END MODULE XASCNTL_MODULE

      ! MODULE XASSTP_MODULE  ! MARK: XASSTP_MODULE
      !   USE XASPDOS_MODULE, ONLY: NPD
      !   USE LINKEDLIST_MODULE, ONLY: LL_TYPE
      !   TYPE STP_TYPE
      !     INTEGER(4) :: NFIL=0
      !     TYPE(LL_TYPE) :: LL
      !     INTEGER(4) :: GID=0
      !     INTEGER(4) :: LNX=0
      !     INTEGER(4) :: NR=0
      !     INTEGER(4) :: NB ! #(STATES IN AE ATOMIC CALCULATION)
      !     INTEGER(4) :: NC ! #(CORE STATES)
      !   END TYPE STP_TYPE
      !   TYPE(STP_TYPE) :: STP(NPD)
      ! END MODULE XASSTP_MODULE

      PROGRAM PAW_XAS  ! MARK: PAW_XAS
      ! USE XASPDOS_MODULE, ONLY: NPD,PD
      USE XASCNTL_MODULE, ONLY: LL_CNTL
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NFIL
      INTEGER(4) :: IPD
      INTEGER(4) :: I
      INTEGER(4) :: NUM_ARGS
!     **************************************************************************
                          CALL TRACE$PUSH('MAIN')
      CALL MPE$INIT
!     INITIALIZE FILES
      CALL INITIALIZEFILEHANDLER

      CALL FILEHANDLER$UNIT('PROT',NFIL)
      CALL CPPAW_WRITEVERSION(NFIL)
!     READ XCNTL FILE
      CALL FILEHANDLER$UNIT('XCNTL',NFIL)
      CALL XASCNTL$READ(NFIL)

      CALL XAS$SELECT('GROUNDSTATE')
      CALL XASCNTL$FILES('GROUNDSTATE')
      CALL XAS$UNSELECT

      CALL XAS$SELECT('EXCITESTATE')
      CALL XASCNTL$FILES('EXCITESTATE')
      CALL XAS$UNSELECT

      CALL XAS$READ

      CALL XAS$SELECT('GROUNDSTATE')
      CALL XAS$REPORT
      CALL XAS$UNSELECT
      CALL XAS$SELECT('EXCITESTATE')
      CALL XAS$REPORT
      CALL XAS$UNSELECT

      CALL DATACONSISTENCY


!     REPORT UNUSED LINKEDLISTS
      CALL FILEHANDLER$UNIT('PROT',NFIL)
      CALL LINKEDLIST$REPORT_UNUSED(LL_CNTL,NFIL)


     
      ! CALL XAS$ISELECT(1)
      ! OPEN(UNIT=NFIL,FILE='gs_1_psphi.dat',STATUS='REPLACE')
      ! CALL XAS$WRITEPHI(NFIL,1,'PSPHI')
      ! CLOSE(NFIL)
      ! OPEN(UNIT=NFIL,FILE='gs_2_psphi.dat',STATUS='REPLACE')
      ! CALL XAS$WRITEPHI(NFIL,2,'PSPHI')
      ! CLOSE(NFIL)
      ! OPEN(UNIT=NFIL,FILE='gs_1_aephi.dat',STATUS='REPLACE')
      ! CALL XAS$WRITEPHI(NFIL,1,'AEPHI')
      ! CLOSE(NFIL)
      ! OPEN(UNIT=NFIL,FILE='gs_2_aephi.dat',STATUS='REPLACE')
      ! CALL XAS$WRITEPHI(NFIL,2,'AEPHI')
      ! CLOSE(NFIL)
      ! CALL XAS$UNSELECT

      ! CALL XAS$ISELECT(2)
      ! OPEN(UNIT=NFIL,FILE='exc_1_psphi.dat',STATUS='REPLACE')
      ! CALL XAS$WRITEPHI(NFIL,1,'PSPHI')
      ! CLOSE(NFIL)
      ! OPEN(UNIT=NFIL,FILE='exc_2_psphi.dat',STATUS='REPLACE')
      ! CALL XAS$WRITEPHI(NFIL,2,'PSPHI')
      ! CLOSE(NFIL)
      ! OPEN(UNIT=NFIL,FILE='exc_1_aephi.dat',STATUS='REPLACE')
      ! CALL XAS$WRITEPHI(NFIL,1,'AEPHI')
      ! CLOSE(NFIL)
      ! OPEN(UNIT=NFIL,FILE='exc_2_aephi.dat',STATUS='REPLACE')
      ! CALL XAS$WRITEPHI(NFIL,2,'AEPHI')
      ! CLOSE(NFIL)
      ! CALL XAS$UNSELECT

      CALL FILEHANDLER$CLOSEALL
                          CALL TRACE$POP
      CALL ERROR$NORMALSTOP
      STOP
      END PROGRAM PAW_XAS

      SUBROUTINE XAS$ISELECT(I)  ! MARK: XAS$ISELECT
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XAS_MODULE, ONLY: NSIM,SELECTED,THIS,SIM
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: I
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$ISELECT')
      IF(I.GT.NSIM.OR.I.LT.0) THEN
        CALL ERROR$MSG('I NOT IN RANGE')
        CALL ERROR$I4VAL('I',I)
        CALL ERROR$I4VAL('NSIM',NSIM)
        CALL ERROR$STOP('XAS$ISELECT')
      END IF
      IF(I.EQ.0) THEN
        IF(.NOT.SELECTED) THEN
          CALL ERROR$MSG('CANNOT UNSELECT A SIMULATION THAT IS NOT SELECTED')
          CALL ERROR$STOP('XAS$ISELECT')
        END IF
        SELECTED=.FALSE.
        NULLIFY(THIS)
      ELSE
        IF(SELECTED) THEN
          CALL ERROR$MSG('ANOTHER SIMULATION IS ALREADY SELECTED')
          CALL ERROR$STOP('XAS$ISELECT')
        END IF
        THIS=>SIM(I)
        SELECTED=.TRUE.
      ENDIF
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE XAS$ISELECT


      SUBROUTINE XAS$SELECT(ID)  ! MARK: XAS$SELECT
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XAS_MODULE, ONLY: SELECTED,THIS,SIM
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$SELECT')
      IF(SELECTED) THEN
        CALL ERROR$MSG('ANOTHER SIMULATION IS ALREADY SELECTED')
        CALL ERROR$CHVAL('SELECTED ID',THIS%ID)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('XAS$SELECT')
      END IF
      SELECTED=.TRUE.
      IF(ID.EQ.'GROUNDSTATE') THEN
        THIS=>SIM(1)
      ELSE IF(ID.EQ.'EXCITESTATE') THEN
        THIS=>SIM(2)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$MSG('ID MUST BE GROUNDSTATE OR EXCITESTATE')
        CALL ERROR$STOP('XAS$SELECT')
      END IF
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE XAS$SELECT


      SUBROUTINE XAS$UNSELECT()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XAS_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$UNSELECT')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CANNOT UNSELECT A SIMULATION THAT IS NOT SELECTED')
        CALL ERROR$STOP('XAS$UNSELECT')
      END IF
      SELECTED=.FALSE.
      NULLIFY(THIS)
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE XAS$UNSELECT
! !
! !     ..................................................................
!       SUBROUTINE XAS$READPDOS(NFIL,IPD)  ! MARK: XAS$READPDOS
! !     ******************************************************************
! !     ** READS PDOS FILE FROM NFIL INTO THE XAS MODULE                **
! !     ** INPUT: NFIL - UNIT NUMBER OF PDOS FILE                       **
! !     **        IPD  - INDEX OF PDOS TO BE READ                       **
! !     ** ADAPTED FROM PDOS$READ                                       **
! !     ******************************************************************
!       USE XASPDOS_MODULE, ONLY: PD,NPD
!       USE LINKEDLIST_MODULE
!       IMPLICIT NONE
!       INTEGER(4),INTENT(IN)  :: NFIL
!       INTEGER(4),INTENT(IN)  :: IPD
!       INTEGER(4)             :: ISP,IKPT,ISPIN,IB
!       INTEGER(4)             :: LNX1,NB
!       INTEGER(4)             :: IOS
!       CHARACTER(82)          :: IOSTATMSG
!       LOGICAL(4)             :: TCHK
!       REAL(8)                :: OCCSUM
!       INTEGER(4)             :: ILOGICAL

!       INTEGER(4)             :: NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,SPACEGROUP
!       LOGICAL(4)             :: TINV,TSHIFT
!       INTEGER(4)                               :: NKDIV(3)
!       INTEGER(4)                               :: ISHIFT(3)
!       REAL(8)                                  :: RNTOT,NEL
!       REAL(8)                                  :: RBAS(3,3)
!       CHARACTER(6)           :: FLAG
! !     ******************************************************************
!                              CALL TRACE$PUSH('XAS$PDOSREAD')
! !     CHECK SELECTION OF PDOS INDEX
!       IF(IPD.LT.1.OR.IPD.GT.NPD)THEN
!         CALL ERROR$MSG('PDOS INDEX NOT IN RANGE')
!         CALL ERROR$I4VAL('IPD',IPD)
!         CALL ERROR$I4VAL('NPD',NPD)
!         CALL ERROR$STOP('XAS$READPDOS')
!       END IF
! !
! !     ==================================================================
! !     == GENERAL QUANTITIES                                           ==
! !     ==================================================================
!       TCHK=.FALSE.
!       READ(NFIL,ERR=100)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG
!       TCHK=.TRUE.
!  100  CONTINUE
!       IF(.NOT.TCHK) THEN
!         PRINT*,'WARNING: NO OCCUPATIONS PRESENT IN PDOS FILE'
!         PRINT*,'            OCCUPATIONS WILL BE SET TO 0'
!         FLAG='LEGACY'
!         REWIND(NFIL)
!         READ(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX
!       END IF
      
!       ALLOCATE(PD(IPD)%LNX(NSP))
!       ALLOCATE(PD(IPD)%LOX(LNXX,NSP))
!       ALLOCATE(PD(IPD)%ISPECIES(NAT))
!       READ(NFIL)PD(IPD)%LNX(:),PD(IPD)%LOX(:,:),PD(IPD)%ISPECIES(:)
      
!       IF(FLAG.EQ.'181213')THEN
! !         == GFORTRAN LOGICAL REPRESENTATION DEFINED WITH TRUE=1, FALSE=0     ==
! !         https://gcc.gnu.org/onlinedocs/gfortran/compiler-characteristics/
! !         internal-representation-of-logical-variables.html
! !         == IFORT LOGICAL REPRESENTATION DEFINED WITH VALUE OF LAST BIT      ==
! !         https://www.intel.com/content/www/us/en/docs/fortran-compiler/
! !         developer-guide-reference/2024-2/logical-data-representations.html
! !         == BOTH SHARE MEANING OF LAST BIT 1=TRUE, 0=FALSE                   ==
! !         == ENSURES BACKWARDS COMPATIBILITY WITH OLD PDOS FILES              ==
!         READ(NFIL)NKDIV(:),ISHIFT(:),RNTOT,NEL,ILOGICAL
!         TINV=BTEST(ILOGICAL,0)
!         READ(NFIL)SPACEGROUP,ILOGICAL
!         TSHIFT=BTEST(ILOGICAL,0)
!       ENDIF
! !
! !     ==================================================================
! !     == ATOMIC STRUCTURE                                             ==
! !     ==================================================================
!       ALLOCATE(PD(IPD)%R(3,NAT))
!       ALLOCATE(PD(IPD)%ATOMID(NAT))
!       READ(NFIL)RBAS(:,:),PD(IPD)%R(:,:),PD(IPD)%ATOMID(:)
! !
! !     ==================================================================
! !     == ELEMENT SPECIFIC QUANTITIES                                  ==
! !     ==================================================================
!       ALLOCATE(PD(IPD)%IZ(NSP))
!       ALLOCATE(PD(IPD)%RAD(NSP)); PD(IPD)%RAD=0.D0
!       ALLOCATE(PD(IPD)%PHIOFR(LNXX,NSP)); PD(IPD)%PHIOFR=0.D0
!       ALLOCATE(PD(IPD)%DPHIDR(LNXX,NSP)); PD(IPD)%DPHIDR=0.D0
!       ALLOCATE(PD(IPD)%OV(LNXX,LNXX,NSP)); PD(IPD)%OV=0.D0
!       DO ISP=1,NSP
!         LNX1=PD(IPD)%LNX(ISP)
!         READ(NFIL)PD(IPD)%IZ(ISP),PD(IPD)%RAD(ISP),PD(IPD)%PHIOFR(1:LNX1,ISP) &
!      &            ,PD(IPD)%DPHIDR(1:LNX1,ISP),PD(IPD)%OV(1:LNX1,1:LNX1,ISP)
!       ENDDO
! !
! !     ==================================================================
! !     ==  NOW READ PROJECTIONS                                       ==
! !     ==================================================================
!       OCCSUM=0.0D0
!       ALLOCATE(PD(IPD)%XK(3,NKPT))
!       ALLOCATE(PD(IPD)%WKPT(NKPT))
!       ALLOCATE(PD(IPD)%STATEARR(NKPT,NSPIN))
!       DO IKPT=1,NKPT
!         DO ISPIN=1,NSPIN
!           PD(IPD)%STATE=>PD(IPD)%STATEARR(IKPT,ISPIN)
!           READ(NFIL,ERR=9998,END=9998)PD(IPD)%XK(:,IKPT),NB,PD(IPD)%WKPT(IKPT)
!           PD(IPD)%STATE%NB=NB
!           ALLOCATE(PD(IPD)%STATE%EIG(NB))
!           ALLOCATE(PD(IPD)%STATE%VEC(NDIM,NPRO,NB))
!           ALLOCATE(PD(IPD)%STATE%OCC(NB))
!           DO IB=1,NB
!             IF(FLAG.EQ.'LEGACY') THEN
!               PD(IPD)%STATE%OCC(:)=0.D0
!               READ(NFIL,ERR=9999,IOSTAT=IOS)PD(IPD)%STATE%EIG(IB),PD(IPD)%STATE%VEC(:,:,IB)
!             ELSE
!               READ(NFIL,ERR=9999,IOSTAT=IOS)PD(IPD)%STATE%EIG(IB) &
!     &                          ,PD(IPD)%STATE%OCC(IB),PD(IPD)%STATE%VEC(:,:,IB)
!               OCCSUM=OCCSUM+PD(IPD)%STATE%OCC(IB)
!             END IF
!           ENDDO
!         ENDDO
!       ENDDO
! PRINT*,"OCCSUM",OCCSUM
! !     SET FIXED SIZE VARIABLES IN PDOS STRUCTURE
!       PD(IPD)%FLAG=FLAG
!       PD(IPD)%NAT=NAT
!       PD(IPD)%NSP=NSP
!       PD(IPD)%NKPT=NKPT
!       PD(IPD)%NSPIN=NSPIN
!       PD(IPD)%NDIM=NDIM
!       PD(IPD)%NPRO=NPRO
!       PD(IPD)%NKDIV(:)=NKDIV(:)
!       PD(IPD)%ISHIFT(:)=ISHIFT(:)
!       PD(IPD)%RNTOT=RNTOT
!       PD(IPD)%NEL=NEL
!       PD(IPD)%TINV=TINV
!       PD(IPD)%LNXX=LNXX
!       PD(IPD)%RBAS(:,:)=RBAS(:,:)
!       PD(IPD)%SPACEGROUP=SPACEGROUP
!       PD(IPD)%TSHIFT=TSHIFT
!                              CALL TRACE$POP
!       RETURN
!  9999 CONTINUE
!       CALL FILEHANDLER$IOSTATMESSAGE(IOS,IOSTATMSG)
!       CALL ERROR$MSG('ERROR READING PDOS FILE')
!       CALL ERROR$I4VAL('IOS',IOS)
!       CALL ERROR$CHVAL('IOSTATMSG',IOSTATMSG)
!       CALL ERROR$I4VAL('IB',IB)
!       CALL ERROR$I4VAL('IKPT',IKPT)
!       CALL ERROR$I4VAL('ISPIN',ISPIN)
!       CALL ERROR$I4VAL('NPRO',NPRO)
!       CALL ERROR$STOP('XAS$READPDOS')
!       STOP
!  9998 CONTINUE
!       CALL ERROR$MSG('ERROR READING PDOS FILE')
!       CALL ERROR$MSG('OLD VERSION: VARIABLE WKPT IS NOT PRESENT')
!       CALL ERROR$MSG('PRODUCE NEW PDOS FILE')
!       CALL ERROR$STOP('XAS$READPDOS')
!       STOP
!       END SUBROUTINE XAS$READPDOS
!
!     ..................................................................
      SUBROUTINE DATACONSISTENCY  ! MARK: DATACONSISTENCY
!     **************************************************************************
!     ** CHECKS CONSISTENCY OF PDOS DATA                                      **
! TODO: IMPLEMENT THIS SUBROUTINE FURTHER, TEST THE CHECKS, ADD ADDITIONAL CHECKS
! TODO: COMBINE WITH TEST SUBROUTINE IN XAS$READ
!     **************************************************************************
      USE XAS_MODULE, ONLY: THIS,SIMULATION_TYPE
      IMPLICIT NONE
      TYPE(SIMULATION_TYPE), POINTER :: THIS1,THIS2
      INTEGER(4) :: I
      REAL(8),PARAMETER :: TOL=1.D-8
!     **************************************************************************
                          CALL TRACE$PUSH('DATACONSISTENCY')
      CALL XAS$SELECT('GROUNDSTATE')
      THIS1=>THIS
      CALL XAS$UNSELECT
      CALL XAS$SELECT('EXCITESTATE')
      THIS2=>THIS
      CALL XAS$UNSELECT
      
      IF(THIS1%NAT.NE.THIS2%NAT)THEN
        CALL ERROR$MSG('NAT INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$I4VAL('GROUNDSTATE%NAT',THIS1%NAT)
        CALL ERROR$I4VAL('EXCITESTATE%NAT',THIS2%NAT)
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
      IF(THIS1%NKPT.NE.THIS2%NKPT)THEN
        CALL ERROR$MSG('NKPT INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$I4VAL('GROUNDSTATE%NKPT',THIS1%NKPT)
        CALL ERROR$I4VAL('EXCITESTATE%NKPT',THIS2%NKPT)
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
      IF(THIS1%NSPIN.NE.THIS2%NSPIN)THEN
        CALL ERROR$MSG('NSPIN INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$I4VAL('GROUNDSTATE%NSPIN',THIS1%NSPIN)
        CALL ERROR$I4VAL('EXCITESTATE%NSPIN',THIS2%NSPIN)
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
      IF(THIS1%NDIM.NE.1.OR.THIS2%NDIM.NE.1)THEN
        CALL ERROR$MSG('ONLY IMPLEMENTED FOR NDIM=1')
        CALL ERROR$I4VAL('GROUNDSTATE%NDIM',THIS1%NDIM)
        CALL ERROR$I4VAL('EXCITESTATE%NDIM',THIS2%NDIM)
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
! TODO: CHECK IF NECESARRY
      IF(THIS1%TINV.NEQV.THIS2%TINV)THEN
        CALL ERROR$MSG('TINV INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$L4VAL('GROUNDSTATE%TINV',THIS1%TINV)
        CALL ERROR$L4VAL('EXCITESTATE%TINV',THIS2%TINV)
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
      IF(ANY(ABS(THIS1%NKDIV-THIS2%NKDIV).NE.0))THEN
        CALL ERROR$MSG('NKDIV INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
      IF(ANY(ABS(THIS1%ISHIFT-THIS2%ISHIFT).NE.0))THEN
        CALL ERROR$MSG('ISHIFT INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
      IF(THIS1%TSHIFT.NEQV.THIS2%TSHIFT)THEN
        CALL ERROR$MSG('TSHIFT INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$L4VAL('GROUNDSTATE%TSHIFT',THIS1%TSHIFT)
        CALL ERROR$L4VAL('EXCITESTATE%TSHIFT',THIS2%TSHIFT)
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
      IF(SUM(ABS(THIS1%RBAS-THIS2%RBAS)).GT.TOL) THEN
        CALL ERROR$MSG('RBAS INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
! TODO: CHECK FOR ATOMIC POSITIONS
! TODO: CHECK FOR XK
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE DATACONSISTENCY

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XASCNTL$READ(NFIL)  ! MARK: XASCNTL$READ
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XASCNTL_MODULE, ONLY: LL_CNTL
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NFIL
!     **************************************************************************
                          CALL TRACE$PUSH('XASCNTL$READ')
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$MARK(LL_CNTL,1)
                          CALL TRACE$POP
      END SUBROUTINE XASCNTL$READ
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XASCNTL$FILES(ID)  ! MARK: XASCNTL$FILES
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XASCNTL_MODULE, ONLY: LL_CNTL
      USE XAS_MODULE, ONLY: THIS,SELECTED
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      CHARACTER(11), INTENT(IN) :: ID ! GROUNDSTATE OR EXCITESTATE
      LOGICAL(4) :: TCHK
      INTEGER(4) :: NFIL
      CHARACTER(256) :: FILENAME
!     **************************************************************************
                          CALL TRACE$PUSH('XASCNTL$FILES')
!     CHECK IF A SIMULATION IS SELECTED
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('XASCNTL$FILES')
      END IF
!     CHECK IF FLAG IS RECOGNIZED
      IF(ID.NE.'GROUNDSTATE'.AND.ID.NE.'EXCITESTATE') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$MSG('ID MUST BE GROUNDSTATE OR EXCITESTATE')
        CALL ERROR$STOP('XASCNTL$FILES')
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'XCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,ID,1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!'//ID//' NOT FOUND IN XCNTL')
        CALL ERROR$STOP('XASCNTL$FILES')
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,ID)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('FILE NOT FOUND IN !'//ID)
        CALL ERROR$STOP('XASCNTL$FILES')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'FILE',0,FILENAME)
    
      CALL XAS$SETCH('FILE',FILENAME)
      CALL XAS$SETCH('ID',ID)
      CALL FILEHANDLER$SETFILE(ID,.FALSE.,TRIM(FILENAME))
                          CALL TRACE$POP
      END SUBROUTINE XASCNTL$FILES
!      
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INITIALIZEFILEHANDLER
!     **************************************************************************
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: XASINNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: NARGS
!     **************************************************************************
                          CALL TRACE$PUSH('INITIALIZEFILEHANDLER')
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE XAS TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GET_COMMAND_ARGUMENT(1,XASINNAME)
      ISVAR=INDEX(XASINNAME,-'.XCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=XASINNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('XCNTL',.FALSE.,XASINNAME)
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
!     **************************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
!     ==========================================================================
!     == SET STANDARD FILENAMES                                               ==
!     ==========================================================================
!     ==  ERROR FILE ===========================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!     ==  PROTOCOL FILE ========================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!     ==  CONTROL FILE  == =====================================================
      ID=+'XCNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!     ==  GROUNDSTATE XAS FILE   ===============================================
      ID=+'GROUNDSTATE'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.GROUND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!     ==  EXCITED STATE XAS FILE   =============================================
      ID=+'EXCITESTATE'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.EXCITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$READ  ! MARK: XAS$READ
!     **************************************************************************
!     **  READ FILE PRODUCED BY SIMULATION CODE                               **
!     **************************************************************************
      USE XAS_MODULE, ONLY: THIS,SELECTED,SIMULATION_TYPE,SIM,OVERLAP,OVERLAPARR
      USE RADIAL_MODULE, ONLY: NGID
      IMPLICIT NONE
      INTEGER(4), PARAMETER :: NSIM=2
      REAL(8), PARAMETER :: TOL=1.D-8
      INTEGER(4) :: NFIL(NSIM)
      INTEGER(4) :: IS
      INTEGER(4) :: ILOGICAL
      INTEGER(4) :: ISP
      INTEGER(4) :: GID
      INTEGER(4) :: GID_
      INTEGER(4) :: IGID
      CHARACTER(8) :: GRIDTYPE
      CHARACTER(8) :: GRIDTYPE_
      INTEGER(4) :: NR
      INTEGER(4) :: NR_
      REAL(8) :: DEX
      REAL(8) :: DEX_
      REAL(8) :: R1
      REAL(8) :: R1_
      INTEGER(4) :: LNX
      INTEGER(4) :: IB
      INTEGER(4) :: IB1,IB2
      INTEGER(4) :: IKPT
      INTEGER(4) :: IG
      INTEGER(4) :: NKPT
      INTEGER(4) :: NKPT_
      INTEGER(4) :: NB_
      INTEGER(4) :: NSPIN_
      INTEGER(4) :: ISPIN
      LOGICAL(4) :: TCHK
      CHARACTER(8) :: KEY(NSIM)
      INTEGER(4) :: NGG(NSIM)
      INTEGER(4) :: NDIM(NSIM)
      INTEGER(4) :: NB(NSIM)
      INTEGER(4) :: NBH(NSIM)
      LOGICAL(4) :: TSUPER(NSIM)
      INTEGER(4), ALLOCATABLE :: IGVEC(:,:,:) ! (NSIM,3,NGG)
      REAL(8) :: XK(NSIM,3)
      COMPLEX(8), ALLOCATABLE :: PSIK1(:,:,:) ! (NGG,NDIM,NB)
      COMPLEX(8), ALLOCATABLE :: PSIK2(:,:,:) ! (NGG,NDIM,NB)
      COMPLEX(8), ALLOCATABLE :: PROJ1(:,:,:)  ! (NDIM,NB,NPRO)
      COMPLEX(8), ALLOCATABLE :: PROJ2(:,:,:)  ! (NDIM,NB,NPRO)
      REAL(8), ALLOCATABLE :: OCC(:,:,:) ! (NB,NKPT,NSPIN)
      REAL(8), ALLOCATABLE :: EIG1(:) ! (NB)
      REAL(8), ALLOCATABLE :: EIG2(:) ! (NB)
      COMPLEX(8), ALLOCATABLE :: OVLAP(:,:) ! (NB,NB)
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$READ')
!     ==========================================================================
!     == LOOP OVER BOTH SIMULATIONS                                           ==
!     == REQUIRED TO CALCULATE OVERLAP ON READ AND NOT STORE PLANE WAVE BASIS ==
!     ==========================================================================
      DO IS=1,NSIM
        CALL XAS$ISELECT(IS)
        CALL FILEHANDLER$UNIT(THIS%ID,NFIL(IS))
        REWIND(NFIL(IS))
!       ========================================================================
!       == READ GENERAL QUANTATIES                                            ==
!       ========================================================================
!       NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG
        READ(NFIL(IS))THIS%NAT,THIS%NSP,THIS%NKPT,THIS%NSPIN,THIS%NDIM,THIS%NPRO, &
       &          THIS%LNXX,THIS%FLAG
        ALLOCATE(THIS%LNX(THIS%NSP))
        ALLOCATE(THIS%LOX(THIS%LNXX,THIS%NSP))
!       LNX(NSP),LOX(LNXX,NSP)
        READ(NFIL(IS))THIS%LNX,THIS%LOX
!       NKDIV(3),ISHIFT(3),RNTOT,NEL,TINV
        READ(NFIL(IS))THIS%NKDIV,THIS%ISHIFT,THIS%RNTOT,THIS%NEL,ILOGICAL
        THIS%TINV=.FALSE.
        IF(ILOGICAL.EQ.1) THIS%TINV=.TRUE.
!       SPACEGROUP,TSHIFT
        READ(NFIL(IS))THIS%SPACEGROUP,ILOGICAL
        THIS%TSHIFT=.FALSE.
        IF(ILOGICAL.EQ.1) THIS%TSHIFT=.TRUE.
!       ========================================================================
!       == READ ATOMIC STRUCTURE                                              ==
!       ========================================================================
        ALLOCATE(THIS%R(3,THIS%NAT))
        ALLOCATE(THIS%ATOMID(THIS%NAT))
        ALLOCATE(THIS%ISPECIES(THIS%NAT))
!       RBAS(3,3),R(3,NAT),ATOMID(NAT),ISPECIES(NAT)
        READ(NFIL(IS))THIS%RBAS,THIS%R,THIS%ATOMID,THIS%ISPECIES
!       ==========================================================================
!       == ELEMENT SPECIFIC QUANTITIES                                          ==
!       ==========================================================================
        ALLOCATE(THIS%SETUP(THIS%NSP))
!       LOOP THROUGH SPECIES
        DO ISP=1,THIS%NSP
!         GRIDTYPE,NR,DEX,R1
          READ(NFIL(IS))GRIDTYPE,NR,DEX,R1
!         CHECK IF RADIAL GRID WITH SAME PROPERTIES ALREADY EXISTS
          TCHK=.FALSE.
          IF(NGID.GT.0) THEN
            DO IGID=1,NGID
              CALL RADIAL$GETCH(IGID,'TYPE',GRIDTYPE_)
              CALL RADIAL$GETI4(IGID,'NR',NR_)
              CALL RADIAL$GETR8(IGID,'DEX',DEX_)
              CALL RADIAL$GETR8(IGID,'R1',R1_)
              IF(GRIDTYPE.EQ.GRIDTYPE_.AND.NR.EQ.NR_.AND.DEX.EQ.DEX_.AND.R1.EQ.R1_) THEN
                GID=IGID
                TCHK=.TRUE.
                EXIT
              END IF
            ENDDO
          ENDIF
!         IF GRID ALREADY EXISTS, USE ITS GID. ELSE CREATE NEW GRID
          IF(TCHK) THEN
            THIS%SETUP(ISP)%GID=GID
          ELSE
            CALL RADIAL$NEW(GRIDTYPE,THIS%SETUP(ISP)%GID)
            CALL RADIAL$SETI4(THIS%SETUP(ISP)%GID,'NR',NR)
            CALL RADIAL$SETR8(THIS%SETUP(ISP)%GID,'DEX',DEX)
            CALL RADIAL$SETR8(THIS%SETUP(ISP)%GID,'R1',R1)
          ENDIF
!         ALLOCATE AND READ (AUXILIARY) PARTIAL WAVES
          LNX=THIS%LNX(ISP)
          ALLOCATE(THIS%SETUP(ISP)%PSPHI(NR,LNX))
          ALLOCATE(THIS%SETUP(ISP)%AEPHI(NR,LNX))
!         PSPHI(NR,LNX)
          READ(NFIL(IS))THIS%SETUP(ISP)%PSPHI
!         AEPHI(NR,LNX)
          READ(NFIL(IS))THIS%SETUP(ISP)%AEPHI
        ENDDO
!       ==================================================================
!       == OCCUPATIONS AND K-POINTS AND THEIR WEIGHTS                   ==
!       ==================================================================
        READ(NFIL(IS))NB_,NKPT_,NSPIN_
        ALLOCATE(THIS%STATEARR(NKPT_,NSPIN_))
        ALLOCATE(THIS%XK(3,NKPT_))
        ALLOCATE(THIS%WKPT(NKPT_))
        ALLOCATE(OCC(NB_,NKPT_,NSPIN_))
!       OCC(NB,NKPT,NSPIN),XK(3,NKPT),WKPT(NKPT)
        READ(NFIL(IS))OCC,THIS%XK,THIS%WKPT
        DO IKPT=1,NKPT_
          DO ISPIN=1,NSPIN_
            THIS%STATE=>THIS%STATEARR(IKPT,ISPIN)
            THIS%STATE%NB=NB_
            ALLOCATE(THIS%STATE%OCC(NB_))
            THIS%STATE%OCC(:)=OCC(:,IKPT,ISPIN)
          ENDDO
        ENDDO
        DEALLOCATE(OCC)
        CALL XAS$UNSELECT
      ENDDO

!     ==================================================================
!     == DATA CONSISTENCY CHECKS BETWEEN BOTH SIMULATIONS             ==
!     ==================================================================
      CALL DATACONSISTENCY

!     ALLOCATE STORAGE ARRAYS FOR OVERLAP_TYPES AND STATE_TYPES
      ALLOCATE(OVERLAPARR(SIM(1)%NKPT,SIM(1)%NSPIN))
!
!     ==================================================================
!     == WAVE FUNCTIONS AND PROJECTIONS                               ==
!     ==================================================================
      NKPT=SIM(1)%NKPT
!     LOOP OVER K POINTS
      DO IKPT=1,NKPT
!       READ GENERAL INFORMATION ABOUT K POINT
        DO IS=1,NSIM
!         KEY,NGG,NDIM,NB,NBH,TSUPER
          READ(NFIL(IS))KEY(IS),NGG(IS),NDIM(IS),NB(IS),ILOGICAL
          TSUPER(IS)=.FALSE.
          IF(ILOGICAL.EQ.1) TSUPER(IS)=.TRUE.
        ENDDO
!       READ K POINT AND G VECTORS (PREVIOUSLY CHECKED IF #(NGG1,NGG2) SAME)
        ALLOCATE(IGVEC(NSIM,3,NGG(1)))
        DO IS=1,NSIM
!
!         XK(3),IGVEC(3,NGG)
          READ(NFIL(IS))XK(IS,:),IGVEC(IS,:,:)
        ENDDO
!       ==================================================================
!       == DATA CHECKS                                                  ==
!       ==================================================================
! NOTE: DATA CHECKS USE VARIABLES FROM THIS SUBROUTINE, POSITION OF CALL MATTERS        
        CALL TEST
        DEALLOCATE(IGVEC)

! WARNING: REQUIRES SUPER WAVE FUNCTIONS TO BE UNRAVELED
! TODO: STORAGE SOLUTION
        ALLOCATE(PSIK1(NGG(1),NDIM(1),NB(1)))
        ALLOCATE(PSIK2(NGG(2),NDIM(2),NB(2)))
        ALLOCATE(EIG1(NB(1)))
        ALLOCATE(EIG2(NB(2)))
!       LOOP OVER SPIN
        DO ISPIN=1,SIM(1)%NSPIN
!         SET STATE POINTER
          DO IS=1,NSIM
!           SET POINTER TO SPECIFIC STATE(IKPT,ISPIN)
            SIM(IS)%STATE=>SIM(IS)%STATEARR(IKPT,ISPIN)
!           CHECK IF NB CONSISTENT FOR OCCUPATIONS AND EIGENVALUES/PROJECTIONS
            IF(NB(IS).NE.SIM(IS)%STATE%NB) THEN
              CALL ERROR$MSG('NB FOR PROJECTIONS AND EIGENVALUES NOT THE SAME')
              CALL ERROR$MSG('FOR OCCUPATIONS AND EIGENVALUES/PROJECTIONS')
              CALL ERROR$I4VAL('NB_OCC',SIM(IS)%STATE%NB)
              CALL ERROR$I4VAL('NB_PROJ',NB(IS))
              CALL ERROR$STOP('XAS$READ')
            END IF
!           ALLOCATE STORAGE FOR PROJECTIONS
            ALLOCATE(SIM(IS)%STATE%PROJ(NDIM(IS),NB(IS),SIM(IS)%NPRO))
!           ALLOCATE STORAGE FOR EIGENVALUES
            ALLOCATE(SIM(IS)%STATE%EIG(NB(IS)))
          ENDDO
!         SET POINTER TO SPECIFIC OVERLAP(IKPT,ISPIN)
          OVERLAP=>OVERLAPARR(IKPT,ISPIN)
! WARNING: FIND OUT REQUIRED ORDER OF OVERLAP
!         ALLOCATE STORAGE FOR OVERLAP
          ALLOCATE(OVERLAP%PW(NB(1),NB(2)))

!         READ PLANE WAVE BASIS
          READ(NFIL(1))PSIK1
          READ(NFIL(2))PSIK2
!         READ PROJECTIONS
          READ(NFIL(1))SIM(1)%STATE%PROJ
          READ(NFIL(2))SIM(2)%STATE%PROJ
!         READ EIGENVALUES
          READ(NFIL(1))SIM(1)%STATE%EIG
          READ(NFIL(2))SIM(2)%STATE%EIG

          DO IB1=1,NB(1) ! LOOP OVER BANDS
            DO IB2=1,NB(2) ! LOOP OVER BANDS
!             NO NDIM LOOP AS NDIM=1
!             SUM OVER G VECTORS
              OVERLAP%PW(IB1,IB2)=SUM(CONJG(PSIK1(:,1,IB1))*PSIK2(:,1,IB2))
            ENDDO ! END LOOP OVER BANDS
          ENDDO ! END LOOP OVER BANDS
        ENDDO ! END LOOP OVER SPIN
        DEALLOCATE(EIG2)
        DEALLOCATE(EIG1)
        DEALLOCATE(PSIK1)
        DEALLOCATE(PSIK2)
      ENDDO ! END LOOP OVER K POINTS
                          CALL TRACE$POP
      
      CONTAINS
        SUBROUTINE TEST!(NSIM,KEY,NGG,NDIM,TSUPER,XK,IGVEC)
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
          CALL ERROR$STOP('XAS$READ')
        END IF
        IF(NGG(1).NE.NGG(2)) THEN
          CALL ERROR$MSG('NGG NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('NGG1',NGG(1))
          CALL ERROR$I4VAL('NGG2',NGG(2))
          CALL ERROR$STOP('XAS$READ')
        END IF
        IF(NDIM(1).NE.NDIM(2)) THEN
          CALL ERROR$MSG('NDIM NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('NDIM1',NDIM(1))
          CALL ERROR$I4VAL('NDIM2',NDIM(2))
          CALL ERROR$STOP('XAS$READ')
        END IF
! TODO: CHECK NDIM AGAINST PREVIOUSLY READ ONE
        IF(TSUPER(1).NEQV.TSUPER(2)) THEN
          CALL ERROR$MSG('TSUPER NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$L4VAL('TSUPER1',TSUPER(1))
          CALL ERROR$L4VAL('TSUPER2',TSUPER(2))
          CALL ERROR$STOP('XAS$READ')
        END IF
!       CHECK IF XK SAME IN BOTH SIMULATIONS
        IF(SUM(ABS(XK(1,:)-XK(2,:)))>TOL) THEN
          CALL ERROR$MSG('XK NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$R8VAL('XK1',XK(1,:))
          CALL ERROR$R8VAL('XK2',XK(2,:))
          CALL ERROR$STOP('XAS$READ')
        END IF
!       CHECK IF IGVEC SAME IN BOTH SIMULATIONS
        IF(ANY(IGVEC(1,:,:).NE.IGVEC(2,:,:))) THEN
          CALL ERROR$MSG('IGVEC NOT THE SAME IN BOTH SIMULATIONS')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$STOP('XAS$READ')
        END IF
        END SUBROUTINE TEST
      END SUBROUTINE XAS$READ
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$REPORT  ! MARK: XAS$REPORT
!     **************************************************************************
!     ** REPORT DATA FOR A SELECTED SIMULATION                                **
!     **************************************************************************
      USE XAS_MODULE, ONLY: THIS, SELECTED
      IMPLICIT NONE
      INTEGER(4) :: NFIL
      INTEGER(4) :: IAT,ISP,ISPIN,IKPT
      CHARACTER(256) :: FORMAT
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$REPORT')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('XAS$REPORT')
      END IF
      CALL FILEHANDLER$UNIT('PROT',NFIL)
      WRITE(NFIL,'(80("#"))')
      WRITE(NFIL,FMT='(A14,A14)')'SIMULATION',THIS%ID
      WRITE(NFIL,'(80("#"))')
      WRITE(NFIL,FMT='(A10,A)')'FILE:',TRIM(THIS%FILE)
      WRITE(NFIL,FMT='(A10,I10)')'NAT:',THIS%NAT
      WRITE(NFIL,FMT='(A10,I10)')'NSP:',THIS%NSP
      WRITE(NFIL,FMT='(A10,I10)')'NKPT:',THIS%NKPT
      WRITE(NFIL,FMT='(A10,I10)')'NSPIN:',THIS%NSPIN
      WRITE(NFIL,FMT='(A10,I10)')'NDIM:',THIS%NDIM
      WRITE(NFIL,FMT='(A10,I10)')'NPRO:',THIS%NPRO
      WRITE(NFIL,FMT='(A10,I10)')'LNXX:',THIS%LNXX
      WRITE(NFIL,FMT='(A10,L10)')'TINV:',THIS%TINV
      WRITE(NFIL,FMT='(A10,3I10)')'NKDIV:',THIS%NKDIV(:)
      WRITE(NFIL,FMT='(A10,3I10)')'ISHIFT:',THIS%ISHIFT(:)
      WRITE(NFIL,FMT='(A10,F10.4)')'RNTOT:',THIS%RNTOT
      WRITE(NFIL,FMT='(A10,F10.4)')'NEL:',THIS%NEL
      WRITE(NFIL,FMT='(A10,I10)')'SPACEGR.:',THIS%SPACEGROUP
      WRITE(NFIL,FMT='(A10,L10)')'TSHIFT:',THIS%TSHIFT
      WRITE(NFIL,FMT='(A10,3F10.4)')'RBAS:',THIS%RBAS(:,1)
      DO IAT=2,3
        WRITE(NFIL,FMT='(A10,3F10.4)')' ',THIS%RBAS(:,IAT)
      ENDDO
      WRITE(NFIL,FMT='(5A10)')'ATOM','X','Y','Z','SPECIES'
      DO IAT=1,THIS%NAT
        WRITE(NFIL,FMT='(A10,3F10.4,I10)')THIS%ATOMID(IAT),THIS%R(:,IAT),THIS%ISPECIES(IAT)
      ENDDO
      WRITE(NFIL,FMT='(4A10)')'SETUP','GID','LNX','LOX'
      DO ISP=1,THIS%NSP
        IF(THIS%LNXX.GT.9) THEN
          WRITE(FORMAT,'("(3I10,",I2,"I10)")')THIS%LNXX
        ELSE
          WRITE(FORMAT,'("(3I10,",I1,"I10)")')THIS%LNXX
        END IF
        WRITE(NFIL,FMT=FORMAT)ISP,THIS%SETUP(ISP)%GID,THIS%LNX(ISP),THIS%LOX(:,ISP)
      ENDDO
      CALL RADIAL$REPORT(NFIL)
      WRITE(NFIL,FMT='(3A10)')'IKPT','ISPIN','NB'
      DO IKPT=1,THIS%NKPT
        DO ISPIN=1,THIS%NSPIN
          WRITE(NFIL,FMT='(3I10)')IKPT,ISPIN,THIS%STATEARR(IKPT,ISPIN)%NB
        ENDDO
      ENDDO
                          CALL TRACE$POP
      END SUBROUTINE XAS$REPORT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$GETCH(ID,VAL)  ! MARK: XAS$GETCH
!     **************************************************************************
!     ** GET CHARACTER IN XAS MODULE                                          **
!     **************************************************************************
      USE XAS_MODULE, ONLY: THIS, SELECTED
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: ID
      CHARACTER(*),INTENT(OUT):: VAL
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('XAS$GETCH')
      END IF
      IF(ID.EQ.'FILE') THEN
        VAL=THIS%FILE
      ELSE IF(ID.EQ.'FLAG') THEN
        VAL=THIS%FLAG
      ELSE IF(ID.EQ.'ID') THEN
        VAL=THIS%ID
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('XAS$GETCH')
      END IF
      RETURN
      END SUBROUTINE XAS$GETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$SETCH(ID,VAL)  ! MARK: XAS$SETCH
!     **************************************************************************
!     ** SET CHARACTER IN XAS MODULE                                          **
!     **************************************************************************
      USE XAS_MODULE, ONLY: THIS, SELECTED
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: ID
      CHARACTER(*) ,INTENT(IN):: VAL
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('XAS$SETCH')
      END IF
      IF(ID.EQ.'FILE') THEN
        THIS%FILE=VAL
      ELSE IF(ID.EQ.'FLAG') THEN
        THIS%FLAG=VAL
      ELSE IF(ID.EQ.'ID') THEN
        THIS%ID=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('XAS$SETCH')
      END IF
      RETURN
      END SUBROUTINE XAS$SETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$WRITEPHI(NFIL,ISP,ID)  ! MARK: XAS$WRITEPHI
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XAS_MODULE, ONLY: THIS  ! SELECTS SIMULATION
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4), INTENT(IN) :: ISP ! SELECT SETUP
      CHARACTER(*) ,INTENT(IN) :: ID ! SELECT AEPHI OR PSPHI
      INTEGER(4) :: NR
      INTEGER(4) :: IR
      REAL(8), ALLOCATABLE :: R(:) ! (NR)
      REAL(8), ALLOCATABLE :: PSI(:,:) ! (NR,LNX)
!     **************************************************************************
      CALL RADIAL$GETI4(THIS%SETUP(ISP)%GID,'NR',NR)
      ALLOCATE(R(NR))
      ALLOCATE(PSI(NR,THIS%LNX(ISP)))
      CALL RADIAL$R(THIS%SETUP(ISP)%GID,NR,R)
      IF(ID.EQ.'AEPHI') THEN
        PSI=THIS%SETUP(ISP)%AEPHI
      ELSE IF(ID.EQ.'PSPHI') THEN
        PSI=THIS%SETUP(ISP)%PSPHI
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('XAS$WRITEPHI')
      END IF
      DO IR=1,NR
        WRITE(NFIL,*) R(IR),PSI(IR,:)
      ENDDO
      DEALLOCATE(PSI)
      DEALLOCATE(R)
      END SUBROUTINE XAS$WRITEPHI


