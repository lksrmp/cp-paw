!........1.........2.........3.........4.........5.........6.........7.........8
      MODULE XAS_MODULE
      TYPE SETUP_TYPE
        INTEGER(4) :: GID
        REAL(8), ALLOCATABLE :: PSPHI(:,:)
        REAL(8), ALLOCATABLE :: AEPHI(:,:)
      END TYPE SETUP_TYPE
      TYPE STATE_TYPE
        INTEGER(4)          :: NB
        REAL(8), POINTER    :: EIG(:)
        REAL(8), POINTER    :: OCC(:)
        COMPLEX(8), POINTER :: PRO(:,:,:)  
      END TYPE STATE_TYPE
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
        TYPE(SETUP_TYPE), ALLOCATABLE :: SETUP(:)
        TYPE(STATE_TYPE), ALLOCATABLE :: STATEARR(:,:)
        TYPE(STATE_TYPE), POINTER :: STATE
      END TYPE SIMULATION_TYPE
! WARNING: CHANGING THIS NUMBER MIGHT BREAK SOME CODE
!     INDEX 1 IS GROUNDSTATE AND 2 EXCITED STATE
      INTEGER(4), PARAMETER :: NSIM=2  ! #(SIMULATIONS)
      LOGICAL,SAVE :: SELECTED=.FALSE.
      TYPE(SIMULATION_TYPE), TARGET :: SIM(NSIM)
      TYPE(SIMULATION_TYPE), POINTER :: THIS
      END MODULE XAS_MODULE


!       MODULE XASPDOS_MODULE  ! MARK: XASPDOS_MODULE
!         TYPE STATE_TYPE
!           INTEGER(4)          :: NB
!           REAL(8), POINTER    :: EIG(:)
!           REAL(8), POINTER    :: OCC(:)
!           COMPLEX(8), POINTER :: VEC(:,:,:)  
!         END TYPE STATE_TYPE
!         TYPE PDOS_TYPE
!           CHARACTER(6)                 :: FLAG
!           INTEGER(4)                   :: NAT
!           INTEGER(4)                   :: NSP
!           INTEGER(4)                   :: NKPT
!           INTEGER(4)                   :: NSPIN
!           INTEGER(4)                   :: NDIM
!           INTEGER(4)                   :: NPRO
!           INTEGER(4)                   :: NKDIV(3)
!           INTEGER(4)                   :: ISHIFT(3)
!           REAL(8)                      :: RNTOT
!           REAL(8)                      :: NEL
!           LOGICAL(4)                   :: TINV
!           INTEGER(4)                   :: LNXX
!           REAL(8)                      :: RBAS(3,3)
!           REAL(8)   ,ALLOCATABLE       :: R(:,:)
!           INTEGER(4),ALLOCATABLE       :: LNX(:)
!           INTEGER(4),ALLOCATABLE       :: LOX(:,:)
!           INTEGER(4),ALLOCATABLE       :: ISPECIES(:)
!           CHARACTER(16),ALLOCATABLE    :: ATOMID(:)
!           INTEGER(4),ALLOCATABLE       :: IZ(:)
!           REAL(8)   ,ALLOCATABLE       :: RAD(:)
!           REAL(8)   ,ALLOCATABLE       :: PHIOFR(:,:)
!           REAL(8)   ,ALLOCATABLE       :: DPHIDR(:,:)
!           REAL(8)   ,ALLOCATABLE       :: OV(:,:,:)
!           REAL(8)   ,ALLOCATABLE       :: XK(:,:)
!           REAL(8)   ,ALLOCATABLE       :: WKPT(:)
!           TYPE(STATE_TYPE),ALLOCATABLE :: STATEARR(:,:)
!           TYPE(STATE_TYPE),POINTER     :: STATE
!           INTEGER(4)                   :: SPACEGROUP
!           LOGICAL(4)                   :: TSHIFT
!         END TYPE PDOS_TYPE
! ! WARNING: CHANGING THIS NUMBER MIGHT BREAK SOME CODE
!         INTEGER(4), PARAMETER :: NPD=2  ! NUMBER OF PDOS FILES
!         TYPE(PDOS_TYPE), TARGET :: PD(NPD) 
!       END MODULE XASPDOS_MODULE

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
!     READ XCNTL FILE
      CALL FILEHANDLER$UNIT('XCNTL',NFIL)
      CALL XASCNTL$READ(NFIL)

      CALL XAS$SELECT('GROUNDSTATE')
      CALL XASCNTL$FILES('GROUNDSTATE')
      CALL XAS$READ
      CALL XAS$REPORT
      CALL XAS$UNSELECT

      CALL XAS$SELECT('EXCITESTATE')
      CALL XASCNTL$FILES('EXCITESTATE')
      CALL XAS$READ
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
      IF(ANY(ABS(THIS1%RBAS-THIS2%RBAS).GT.TOL)) THEN
        CALL ERROR$MSG('RBAS INCONSISTENT BETWEEN SIMULATIONS')
        CALL ERROR$STOP('DATACONSISTENCY')
      END IF
! TODO: CHECK FOR ATOMIC POSITIONS
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
      USE XAS_MODULE, ONLY: THIS,SELECTED
      USE RADIAL_MODULE, ONLY: NGID
      IMPLICIT NONE
      INTEGER(4) :: NFIL
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
      INTEGER(4) :: NB
      INTEGER(4) :: IB
      INTEGER(4) :: IKPT
      INTEGER(4) :: ISPIN
      LOGICAL(4) :: TCHK
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$READ')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('XAS$READ')
      END IF
      CALL FILEHANDLER$UNIT(THIS%ID,NFIL)
      REWIND(NFIL)
!     ==========================================================================
!     == READ GENERAL QUANTATIES                                              ==
!     ==========================================================================
      READ(NFIL)THIS%NAT,THIS%NSP,THIS%NKPT,THIS%NSPIN,THIS%NDIM,THIS%NPRO, &
     &          THIS%LNXX,THIS%FLAG
      ALLOCATE(THIS%LNX(THIS%NSP))
      ALLOCATE(THIS%LOX(THIS%LNXX,THIS%NSP))
      READ(NFIL)THIS%LNX(:),THIS%LOX(:,:)
      READ(NFIL)THIS%NKDIV(:),THIS%ISHIFT(:),THIS%RNTOT,THIS%NEL,ILOGICAL
      THIS%TINV=.FALSE.
      IF(ILOGICAL.EQ.1) THIS%TINV=.TRUE.
      READ(NFIL)THIS%SPACEGROUP,ILOGICAL
      THIS%TSHIFT=.FALSE.
      IF(ILOGICAL.EQ.1) THIS%TSHIFT=.TRUE.
!     ==========================================================================
!     == READ ATOMIC STRUCTURE                                                ==
!     ==========================================================================
      ALLOCATE(THIS%R(3,THIS%NAT))
      ALLOCATE(THIS%ATOMID(THIS%NAT))
      ALLOCATE(THIS%ISPECIES(THIS%NAT))
      READ(NFIL)THIS%RBAS(:,:),THIS%R(:,:),THIS%ATOMID(:),THIS%ISPECIES(:)
!     ==========================================================================
!     == ELEMENT SPECIFIC QUANTITIES                                          ==
!     ==========================================================================
      ALLOCATE(THIS%SETUP(THIS%NSP))
!     LOOP THROUGH SPECIES
      DO ISP=1,THIS%NSP
        READ(NFIL)GRIDTYPE,NR,DEX,R1
!       CHECK IF RADIAL GRID WITH SAME PROPERTIES ALREADY EXISTS
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
!       IF GRID ALREADY EXISTS, USE ITS GID. ELSE CREATE NEW GRID
        IF(TCHK) THEN
          THIS%SETUP(ISP)%GID=GID
        ELSE
          CALL RADIAL$NEW(GRIDTYPE,THIS%SETUP(ISP)%GID)
          CALL RADIAL$SETI4(THIS%SETUP(ISP)%GID,'NR',NR)
          CALL RADIAL$SETR8(THIS%SETUP(ISP)%GID,'DEX',DEX)
          CALL RADIAL$SETR8(THIS%SETUP(ISP)%GID,'R1',R1)
        ENDIF
!       ALLOCATE AND READ (AUXILIARY) PARTIAL WAVES
        LNX=THIS%LNX(ISP)
        ALLOCATE(THIS%SETUP(ISP)%PSPHI(NR,LNX))
        ALLOCATE(THIS%SETUP(ISP)%AEPHI(NR,LNX))
        READ(NFIL)THIS%SETUP(ISP)%PSPHI(:,:)
        READ(NFIL)THIS%SETUP(ISP)%AEPHI(:,:)
      ENDDO
!     ==========================================================================
!     == READ PROJECTIONS                                                     ==
!     ==========================================================================
      ALLOCATE(THIS%XK(3,THIS%NKPT))
      ALLOCATE(THIS%WKPT(THIS%NKPT))
      ALLOCATE(THIS%STATEARR(THIS%NKPT,THIS%NSPIN))
      DO IKPT=1,THIS%NKPT
        DO ISPIN=1,THIS%NSPIN
          THIS%STATE=>THIS%STATEARR(IKPT,ISPIN)
          READ(NFIL)THIS%XK(:,IKPT),NB,THIS%WKPT(IKPT)
          THIS%STATE%NB=NB
          ALLOCATE(THIS%STATE%EIG(NB))
          ALLOCATE(THIS%STATE%PRO(THIS%NDIM,THIS%NPRO,NB))
          ALLOCATE(THIS%STATE%OCC(NB))
          DO IB=1,NB
            READ(NFIL)THIS%STATE%EIG(IB)
            READ(NFIL)THIS%STATE%OCC(IB)
            READ(NFIL)THIS%STATE%PRO(:,:,IB)
          ENDDO
        ENDDO ! END LOOP OVER SPIN
      ENDDO ! END LOOP OVER KPOINTS
                          CALL TRACE$POP
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
      INTEGER(4) :: IAT,ISP
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
      INTEGER(4) :: IR
      REAL(8), ALLOCATABLE :: R(:) ! (NR)
      REAL(8), ALLOCATABLE :: PSI(:,:) ! (NR,LNX)
!     **************************************************************************
      ALLOCATE(R(THIS%SETUP(ISP)%NR))
      ALLOCATE(PSI(THIS%SETUP(ISP)%NR,THIS%SETUP(ISP)%LNX))
      CALL RADIAL$R(THIS%SETUP(ISP)%GID,THIS%SETUP(ISP)%NR,R)
      IF(ID.EQ.'AEPHI') THEN
        PSI=THIS%SETUP(ISP)%AEPHI
      ELSE IF(ID.EQ.'PSPHI') THEN
        PSI=THIS%SETUP(ISP)%PSPHI
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('XAS$WRITEPHI')
      END IF
      DO IR=1,THIS%SETUP(ISP)%NR
        WRITE(NFIL,*) R(IR),PSI(IR,:)
      ENDDO
      DEALLOCATE(PSI)
      DEALLOCATE(R)
      END SUBROUTINE XAS$WRITEPHI


