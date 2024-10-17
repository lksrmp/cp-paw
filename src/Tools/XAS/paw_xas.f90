!........1.........2.........3.........4.........5.........6.........7.........8
      MODULE XAS_MODULE
      TYPE SETUP_TYPE
        INTEGER(4) :: GID ! GRID ID FOR RADIAL MESH
        REAL(8), ALLOCATABLE :: PSPHI(:,:) ! PSEUDO PARTIAL WAVE
        REAL(8), ALLOCATABLE :: AEPHI(:,:) ! AE PARTIAL WAVE
        INTEGER(4) :: NBATOM ! #(ATOMIC WAVE FUNCTIONS)
        INTEGER(4), ALLOCATABLE :: LATOM(:) ! (NBATOM) L QUANTUM NUMBER
        REAL(8), ALLOCATABLE :: AEPSI(:,:) ! (NR,NBATOM) AE ATOMIC WAVE FUNCTION
      END TYPE SETUP_TYPE

      TYPE STATE_TYPE
        REAL(8) :: OCCTOL=1.D-3 ! OCCUPATION > OCCTOL IS OCCUPIED
        INTEGER(4)          :: NB ! #(STATES)
        REAL(8), POINTER    :: EIG(:) ! (NB) EIGENVALUES
        REAL(8), POINTER    :: OCC(:) ! (NB) OCCUPATIONS
        INTEGER(4)          :: NOCC ! #(OCCUPIED STATES)
        COMPLEX(8), POINTER :: PROJ(:,:,:) ! (NDIM,NB,NPRO) PROJECTIONS
      END TYPE STATE_TYPE

      TYPE OVERLAP_TYPE
        COMPLEX(8), ALLOCATABLE :: PW(:,:) ! (NB2,NB1) PLANE WAVE OVERLAP
        COMPLEX(8), ALLOCATABLE :: AUG(:,:) ! (NB2,NB1) AUGMENTATION OVERLAP
        COMPLEX(8), ALLOCATABLE :: OV(:,:) ! (NB2,NB1) OVERLAP
!       IT IS CHECKED THAT NOCC1=NOCC2
        COMPLEX(8), ALLOCATABLE :: OVOCC(:,:) ! (NOCC2,NOCC1) OCCUPIED OVERLAP
        COMPLEX(8), ALLOCATABLE :: OVEMP(:,:) ! (NB2-NOCC2,NOCC1) EMPTY OVERLAP
        COMPLEX(8), ALLOCATABLE :: KMAT(:,:) ! (NB2-NOCC2,NOCC1) K MATRIX ELEMENTS
        COMPLEX(8), ALLOCATABLE :: DIPOLE(:,:) ! (3,NB2) DIPOLE MATRIX ELEMENTS
        COMPLEX(8) :: ADET ! DETERMINANT OF OCCUPIED OVERLAP
       END TYPE OVERLAP_TYPE

      TYPE SETTINGS_TYPE
        INTEGER(4) :: NSPEC ! #(SPECTRA)
        CHARACTER(256) :: ATOM ! ATOM ID WITH CORE HOLE
        INTEGER(4) :: IATOM ! ATOM INDEX WITH CORE HOLE
        INTEGER(4) :: NCORE ! N QUANTUM NUMBER OF CORE HOLE
        INTEGER(4) :: LCORE ! L QUANTUM NUMBER OF CORE HOLE
      END TYPE SETTINGS_TYPE

      TYPE SPECTRUM_TYPE
        REAL(8) :: NORMAL(3) ! NORMAL VECTOR OF COORDINATES
        REAL(8) :: KDIR(3) ! DIRECTION OF K-VECTOR
        COMPLEX(8) :: POL(2) ! POLARIZATION VECTOR
        COMPLEX(8) :: POLXYZ(3) ! POLARIZATION VECTOR IN CARTESIAN COORDINATES
      END TYPE SPECTRUM_TYPE

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
        INTEGER(4), ALLOCATABLE :: MAP(:,:) ! (NAT,LNXX) INDEX-1 OF ATOM AND LN
        LOGICAL(4) :: TINV ! INVERSION SYMMETRY
        INTEGER(4) :: NKDIV(3) ! K-POINT DIVISIONS
        INTEGER(4) :: ISHIFT(3) ! K-POINT SHIFTS
        REAL(8) :: RNTOT
        REAL(8) :: NEL
        REAL(8) :: ETOT
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

      TYPE(SETTINGS_TYPE) :: SETTINGS

      TYPE(SPECTRUM_TYPE), ALLOCATABLE, TARGET :: SPECTRUMARR(:) ! (NSPEC)
      TYPE(SPECTRUM_TYPE), POINTER :: SPECTRUM ! CURRENT SPECTRUM

      REAL(8), ALLOCATABLE :: S(:,:,:) ! (NAT,LNXX1,LNXX2) ATOMIC OVERLAP MATRIX

      LOGICAL(4) :: TSIM=.FALSE.
      LOGICAL(4) :: TOVERLAP=.FALSE.
      LOGICAL(4) :: TSPECTRUM=.FALSE.
      LOGICAL(4) :: TSETTINGS=.FALSE.
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
!     ==========================================================================
!     == READ XCNTL FILE                                                      ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('XCNTL',NFIL)
      CALL XASCNTL$READ(NFIL)
!     READ FILES FOR DFT CALCULATION DATA
      CALL XAS$SELECT('GROUNDSTATE')
      CALL XASCNTL$FILES('GROUNDSTATE')
      CALL XAS$UNSELECT
      CALL XAS$SELECT('EXCITESTATE')
      CALL XASCNTL$FILES('EXCITESTATE')
      CALL XAS$UNSELECT
!     READ SETTINGS FOR SPECTRA
      CALL XASCNTL$SPECTRUM
!     CALCULATE POLARISATION IN CARTESIAN COORDINATES FROM NORMAL, K, AND POL
      CALL XAS$POLARISATION
!     READ DFT CALCULATION DATA
      CALL XAS$READ
!     REPORT DFT CALCULATION DATA
      CALL XAS$SELECT('GROUNDSTATE')
      CALL XAS$REPORTSIMULATION
      CALL XAS$UNSELECT
      CALL XAS$SELECT('EXCITESTATE')
      CALL XAS$REPORTSIMULATION
      CALL XAS$UNSELECT

      CALL DATACONSISTENCY

!     READ SETTINGS FOR CORE HOLE
      CALL XASCNTL$COREHOLE

      CALL XAS$OVERLAP

      CALL XAS$DIPOLEMATRIX

      CALL XAS$ADET
      CALL XAS$KMAT

      CALL XAS$REPORTSIMULATION

      CALL XAS$REPORTSETTINGS

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
      RETURN
      END SUBROUTINE XAS$SELECT


      SUBROUTINE XAS$UNSELECT()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XAS_MODULE, ONLY: SELECTED,THIS
      IMPLICIT NONE
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CANNOT UNSELECT A SIMULATION THAT IS NOT SELECTED')
        CALL ERROR$STOP('XAS$UNSELECT')
      END IF
      SELECTED=.FALSE.
      NULLIFY(THIS)
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
! TODO: IMPLEMENT VARIABLE ORDER OF ATOMS
      DO I=1,THIS1%NAT
        IF(SUM(ABS(THIS1%R(:,I)-THIS2%R(:,I))).GT.TOL) THEN
          CALL ERROR$MSG('ATOMIC POSITIONS INCONSISTENT BETWEEN SIMULATIONS')
          CALL ERROR$I4VAL('I',I)
          CALL ERROR$STOP('DATACONSISTENCY')
        END IF
      ENDDO
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
      SUBROUTINE XASCNTL$SPECTRUM  ! MARK: XASCNTL$SPECTRUM
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XASCNTL_MODULE, ONLY: LL_CNTL
      USE XAS_MODULE, ONLY: SPECTRUM,SPECTRUMARR,TSPECTRUM,SETTINGS
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      LOGICAL(4) :: TCHK
      INTEGER(4) :: NSPEC
      INTEGER(4) :: ISPEC
      REAL(8) :: REALPOL(2)
!     *************************************************************************
! WARNING: READING OF COMPLEX POLARIZATION IS NOT IMPLEMENTED
! TODO: IMPLEMENT READING OF COMPLEX POLARIZATION
                          CALL TRACE$PUSH('XASCNTL$SPECTRUM')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'XCNTL')
!     CHECK IF AT LEAST ONE SPECTRUM IS DEFINED
      CALL LINKEDLIST$EXISTL(LL_CNTL,'SPECTRUM',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('SPECTRUM NOT FOUND IN XCNTL')
        CALL ERROR$STOP('XASCNTL$SPECTRUM')
      END IF
!     ALLOCATE CORRECT NUMBER OF SPECTRA
      CALL LINKEDLIST$NLISTS(LL_CNTL,'SPECTRUM',NSPEC)
      ALLOCATE(SPECTRUMARR(NSPEC))
      SETTINGS%NSPEC=NSPEC
!     LOOP OVER SPECTRA AND READ DATA
      DO ISPEC=1,NSPEC
        SPECTRUM=>SPECTRUMARR(ISPEC)
        CALL LINKEDLIST$SELECT(LL_CNTL,'~')
        CALL LINKEDLIST$SELECT(LL_CNTL,'XCNTL')
        CALL LINKEDLIST$SELECT(LL_CNTL,'SPECTRUM',ISPEC)
!       NORMAL OF THE SYSTEM
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NORMAL',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('NORMAL NOT FOUND IN !SPECTRUM')
          CALL ERROR$I4VAL('ISPEC',ISPEC)
          CALL ERROR$STOP('XASCNTL$SPECTRUM')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'NORMAL',1,SPECTRUM%NORMAL)
!       DIRECTION OF THE K-VECTOR
        CALL LINKEDLIST$EXISTD(LL_CNTL,'KDIR',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('KDIR NOT FOUND IN !SPECTRUM')
          CALL ERROR$I4VAL('ISPEC',ISPEC)
          CALL ERROR$STOP('XASCNTL$SPECTRUM')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'KDIR',1,SPECTRUM%KDIR)
!       POLARIZATION OF THE LIGHT
        CALL LINKEDLIST$EXISTD(LL_CNTL,'POL',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('POL NOT FOUND IN !SPECTRUM')
          CALL ERROR$I4VAL('ISPEC',ISPEC)
          CALL ERROR$STOP('XASCNTL$SPECTRUM')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'POL',1,REALPOL)
        SPECTRUM%POL=CMPLX(REALPOL,KIND=8)
      ENDDO
      TSPECTRUM=.TRUE.
                          CALL TRACE$POP
      END SUBROUTINE XASCNTL$SPECTRUM
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XASCNTL$COREHOLE  ! MARK: XASCNTL$COREHOLE
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XASCNTL_MODULE, ONLY: LL_CNTL
      USE XAS_MODULE, ONLY: SETTINGS,SIM,TSIM,TSETTINGS
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      LOGICAL(4) :: TCHK
      INTEGER(4) :: NCORE
      INTEGER(4) :: LCORE
      CHARACTER(256) :: ATOM
      INTEGER(4) :: IAT
                          CALL TRACE$PUSH('XASCNTL$COREHOLE')
      IF(.NOT.TSIM) THEN
        CALL ERROR$MSG('SIMULATION DATA NOT ALLOCATED')
        CALL ERROR$MSG('MUST BE CALLED AFTER XAS$READ TO ACCESS ATOM NAMES')
        CALL ERROR$STOP('XASCNTL$COREHOLE')
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'XCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'COREHOLE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!COREHOLE NOT FOUND IN !XCNTL')
        CALL ERROR$STOP('XASCNTL$COREHOLE')
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'COREHOLE')

      CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('ATOM NOT FOUND IN !COREHOLE')
        CALL ERROR$STOP('XASCNTL$COREHOLE')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,ATOM)
      SETTINGS%ATOM=TRIM(ATOM)
      TCHK=.FALSE.
      DO IAT=1,SIM(2)%NAT
        IF(TRIM(ATOM).EQ.TRIM(SIM(2)%ATOMID(IAT))) THEN
          TCHK=.TRUE.
          EXIT
        END IF
      ENDDO
      IF(TCHK) THEN
        SETTINGS%IATOM=IAT
      ELSE
        CALL ERROR$MSG('ATOM NOT FOUND IN EXCITESTATE')
        CALL ERROR$CHVAL('ATOM',ATOM)
        CALL ERROR$STOP('XASCNTL$COREHOLE')
      END IF

      CALL LINKEDLIST$EXISTD(LL_CNTL,'NCORE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('NCORE NOT FOUND IN !COREHOLE')
        CALL ERROR$STOP('XASCNTL$COREHOLE')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'NCORE',1,NCORE)
      SETTINGS%NCORE=NCORE
! TODO: CHECK IF LCORE IS VALID FOR GIVEN NCORE
      CALL LINKEDLIST$EXISTD(LL_CNTL,'LCORE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('LCORE NOT FOUND IN !COREHOLE')
        CALL ERROR$STOP('XASCNTL$COREHOLE')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'LCORE',1,LCORE)
      IF(LCORE.NE.0) THEN
        CALL ERROR$MSG('IMPLEMENTATION ONLY FOR S COREHOLES (LCORE=0)')
        CALL ERROR$I4VAL('LCORE',LCORE)
        CALL ERROR$STOP('XASCNTL$COREHOLE')
      END IF
      SETTINGS%LCORE=LCORE

      TSETTINGS=.TRUE.
                          CALL TRACE$POP
      END SUBROUTINE XASCNTL$COREHOLE
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
      USE XAS_MODULE, ONLY: THIS,SELECTED,SIMULATION_TYPE,SIM,TSIM, &
     &                      OVERLAP,OVERLAPARR
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
      REAL(8) :: V
      INTEGER(4) :: LNX
      INTEGER(4) :: IB
      INTEGER(4) :: IAT
      INTEGER(4) :: LN
      INTEGER(4) :: L
      INTEGER(4) :: M
      INTEGER(4) :: IPRO
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
        READ(NFIL(IS))THIS%NKDIV,THIS%ISHIFT,THIS%RNTOT,THIS%NEL,THIS%ETOT,ILOGICAL
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
        V=THIS%RBAS(1,1)*(THIS%RBAS(2,2)*THIS%RBAS(3,3)-THIS%RBAS(2,3)*THIS%RBAS(3,2)) &
     &  +THIS%RBAS(2,1)*(THIS%RBAS(3,2)*THIS%RBAS(1,3)-THIS%RBAS(3,3)*THIS%RBAS(1,2)) &
     &  +THIS%RBAS(3,1)*(THIS%RBAS(1,2)*THIS%RBAS(2,3)-THIS%RBAS(1,3)*THIS%RBAS(2,2)) 
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
!         NBATOM
          READ(NFIL(IS))THIS%SETUP(ISP)%NBATOM
!         ALLOCATE AND READ ATOMIC WAVE FUNCTIONS
          ALLOCATE(THIS%SETUP(ISP)%LATOM(THIS%SETUP(ISP)%NBATOM))
          ALLOCATE(THIS%SETUP(ISP)%AEPSI(NR,THIS%SETUP(ISP)%NBATOM))
!         LATOM(NBATOM)
          READ(NFIL(IS))THIS%SETUP(ISP)%LATOM       
!         AEPSI(NR,NBATOM)
          READ(NFIL(IS))THIS%SETUP(ISP)%AEPSI
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
!           COUNT NUMBER OF OCCUPIED STATES
            THIS%STATE%NOCC=-1
            DO IB=1,NB_
              IF(THIS%STATE%OCC(IB).LT.THIS%STATE%OCCTOL) THEN
                THIS%STATE%NOCC=IB-1
                EXIT
              ENDIF
            ENDDO
            IF(THIS%STATE%NOCC.EQ.-1) THEN
              THIS%STATE%NOCC=NB_
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(OCC)
!       ==================================================================
!       == MAPPING OF PROJECTION INDICES                                ==
!       ==================================================================
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
            ENDDO
          ENDDO
        ENDDO
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
          IF(KEY(IS).NE.'PSI') THEN
            CALL ERROR$MSG('KEY NOT "PSI"')
            CALL ERROR$MSG('FILE IS CORRUPTED')
            CALL ERROR$CHVAL('KEY',KEY(IS))
            CALL ERROR$I4VAL('IKPT',IKPT)
            CALL ERROR$STOP('XAS$READ')
          END IF
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
          ALLOCATE(OVERLAP%PW(NB(2),NB(1)))

!         READ PLANE WAVE BASIS
          READ(NFIL(1))PSIK1
          READ(NFIL(2))PSIK2
!         READ PROJECTIONS
          READ(NFIL(1))SIM(1)%STATE%PROJ
          READ(NFIL(2))SIM(2)%STATE%PROJ
!         READ EIGENVALUES
          READ(NFIL(1))SIM(1)%STATE%EIG
          READ(NFIL(2))SIM(2)%STATE%EIG

          DO IB2=1,NB(2) ! LOOP OVER BANDS
            DO IB1=1,NB(1) ! LOOP OVER BANDS
!             NO NDIM LOOP AS NDIM=1
!             SCALARPRODUCT (SUM OVER G VECTORS)
!             PW(I,J)=<PSI1(J)|PSI2(I)>
! WARNING: CHECK IF SCALARPRODUCT IS CORRECT ALSO WITH COMJG
              CALL LIB$SCALARPRODUCTC8(.FALSE.,NGG(1),1,PSIK1(:,1,IB1),1,PSIK2(:,1,IB2),OVERLAP%PW(IB2,IB1))
            ENDDO ! END LOOP OVER BANDS
          ENDDO ! END LOOP OVER BANDS
          OVERLAP%PW=OVERLAP%PW*V
!           DO IB1=1,NB(1) ! LOOP OVER BANDS
!             DO IB2=1,NB(2) ! LOOP OVER BANDS
! !             NO NDIM LOOP AS NDIM=1
! !             SUM OVER G VECTORS
! !             PW(I,J)=<PSI1(I)|PSI2(J)>
!               OVERLAP%PW(IB1,IB2)=V*SUM(CONJG(PSIK1(:,1,IB1))*PSIK2(:,1,IB2))
!             ENDDO ! END LOOP OVER BANDS
!           ENDDO ! END LOOP OVER BANDS

        ENDDO ! END LOOP OVER SPIN
        DEALLOCATE(EIG2)
        DEALLOCATE(EIG1)
        DEALLOCATE(PSIK1)
        DEALLOCATE(PSIK2)
      ENDDO ! END LOOP OVER K POINTS
      TSIM=.TRUE.
                          CALL TRACE$POP
      RETURN

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
      SUBROUTINE XAS$REPORTSIMULATION  ! MARK: XAS$REPORTSIMULATION
!     **************************************************************************
!     ** REPORT DATA FOR A SELECTED SIMULATION                                **
!     **************************************************************************
      USE XAS_MODULE, ONLY: THIS,SELECTED,S,OVERLAP,OVERLAPARR
      IMPLICIT NONE
      INTEGER(4) :: NFIL
! TODO: REMOVE HARDCODED OUTPUT
      integer(4) :: nfilo=11
      integer(4) :: nfilc=12
      INTEGER(4) :: NB1,NB2
      INTEGER(4) :: IAT,ISP,ISPIN,IKPT,IPRO,IB
      CHARACTER(256) :: FORMAT
      REAL(8) :: EV
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$REPORTSIMULATION')
      CALL CONSTANTS('EV',EV)
      CALL FILEHANDLER$UNIT('PROT',NFIL)
      IF(SELECTED) THEN
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
        WRITE(NFIL,FMT='(A10,F10.4)')'ETOT [EV]:',THIS%ETOT/EV
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
        IF(THIS%LNXX.GT.9) THEN
          WRITE(FORMAT,'("(3I10,",I2,"I10)")')THIS%LNXX
        ELSE
          WRITE(FORMAT,'("(3I10,",I1,"I10)")')THIS%LNXX
        END IF
        DO ISP=1,THIS%NSP
          WRITE(NFIL,FMT=FORMAT)ISP,THIS%SETUP(ISP)%GID,THIS%LNX(ISP),THIS%LOX(:,ISP)
        ENDDO
        WRITE(NFIL,FMT='(A)')'MAPPING OF PROJECTIONS (INDEX-1)'
        WRITE(NFIL,FMT='(A10,A10)')'ATOM','MAP'
        IF(THIS%LNXX.GT.9) THEN
          WRITE(FORMAT,'("(I10,",I2,"I10)")')THIS%LNXX
        ELSE
          WRITE(FORMAT,'("(I10,",I1,"I10)")')THIS%LNXX
        END IF
        DO IAT=1,THIS%NAT
          WRITE(NFIL,FMT=FORMAT)IAT,THIS%MAP(IAT,:)
        ENDDO
        CALL RADIAL$REPORT(NFIL)
        WRITE(NFIL,FMT='(3A10)')'IKPT','ISPIN','NB'
        DO IKPT=1,THIS%NKPT
          DO ISPIN=1,THIS%NSPIN
            WRITE(NFIL,FMT='(3I10)')IKPT,ISPIN,THIS%STATEARR(IKPT,ISPIN)%NB
          ENDDO
        ENDDO
      ELSE
        WRITE(NFIL,'(80("#"))')
        WRITE(NFIL,FMT='(A19)')'GENERAL INFORMATION'
        WRITE(NFIL,'(80("#"))')
        CALL XAS$SELECT('EXCITESTATE')
        NB2=THIS%STATE%NB
        IF(THIS%LNXX.GT.9) THEN
          WRITE(FORMAT,'("(",I2,"F10.4)")')THIS%LNXX
        ELSE
          WRITE(FORMAT,'("(",I1,"F10.4)")')THIS%LNXX
        END IF
        CALL XAS$UNSELECT
        CALL XAS$SELECT('GROUNDSTATE')
        NB1=THIS%STATE%NB
        WRITE(NFIL,FMT='(A)')'ATOMIC OVERLAP MATRIX S'
        DO IAT=1,THIS%NAT
          WRITE(NFIL,FMT='(A10,I10)')'ATOM',IAT
          DO IPRO=1,THIS%LNXX
            WRITE(NFIL,FMT=FORMAT)S(IAT,IPRO,:)
          ENDDO
        ENDDO
! TODO: PROPER OUTPUT FOR OVERLAP MATRIX
        OVERLAP=>OVERLAPARR(1,1)
        IF(NB1.GT.99) THEN
          WRITE(FORMAT,'("(",I3,"F10.6)")')NB1
        ELSE IF(NB1.GT.9.AND.NB1.LT.100) THEN
          WRITE(FORMAT,'("(",I2,"F10.6)")')NB1
        ELSE
          WRITE(FORMAT,'("(",I1,"F10.6)")')NB1
        END IF
        open(nfilo,file='pw.dat')
        open(nfilc,file='pw_c.dat')
        WRITE(NFIL,FMT='(A)')'PLANE WAVE OVERLAP MATRIX OF FIRST K POINT/SPIN'
        DO IB=1,NB2
          WRITE(NFIL,FMT=FORMAT)CDABS(OVERLAP%PW(IB,:))
          WRITE(nfilo,FMT=FORMAT)CDABS(OVERLAP%PW(IB,:))
          write(nfilc,*)OVERLAP%PW(IB,:)
        ENDDO
        close(nfilo)
        close(nfilc)
        open(nfilo,file='aug.dat')
        open(nfilc,file='aug_c.dat')
        WRITE(NFIL,FMT='(A)')'AUGMENTATION OVERLAP MATRIX OF FIRST K POINT/SPIN'
        DO IB=1,NB2
          WRITE(NFIL,FMT=FORMAT)CDABS(OVERLAP%AUG(IB,:))
          WRITE(nfilo,FMT=FORMAT)CDABS(OVERLAP%AUG(IB,:))
          write(nfilc,*)OVERLAP%AUG(IB,:)
        ENDDO
        close(nfilo)
        close(nfilc)
        WRITE(NFIL,FMT='(A)')'OVERLAP MATRIX OF FIRST K POINT/SPIN'
        open(nfilo,file='ov.dat')
        open(nfilc,file='ov_c.dat')
        DO IB=1,NB2
          WRITE(NFIL,FMT=FORMAT)CDABS(OVERLAP%OV(IB,:))
          WRITE(nfilo,FMT=FORMAT)CDABS(OVERLAP%OV(IB,:))
          write(nfilc,*)OVERLAP%OV(IB,:)
        ENDDO
        close(nfilo)
        close(nfilc)
        WRITE(NFIL,FMT='(A)')'OCCUPIED OVERLAP OF FIRST K POINT/SPIN'
        DO IB=1,THIS%STATEARR(1,1)%NOCC
          WRITE(NFIL,FMT=FORMAT)CDABS(OVERLAP%OVOCC(IB,:))
        ENDDO
        WRITE(NFIL,FMT='(A)')'EMPTY OVERLAP OF FIRST K POINT/SPIN'
        DO IB=1,NB2-THIS%STATEARR(1,1)%NOCC
          WRITE(NFIL,FMT=FORMAT)CDABS(OVERLAP%OVEMP(IB,:))
        ENDDO
        WRITE(NFIL,FMT='(A)')'KMAT OF FIRST K POINT/SPIN'
        DO IB=1,NB2-THIS%STATEARR(1,1)%NOCC
          WRITE(NFIL,FMT=FORMAT)CDABS(OVERLAP%KMAT(IB,:))
        ENDDO
! TODO: OUTPUT FOR DIPOLE ELEMENTS
        CALL XAS$UNSELECT
      ENDIF
                          CALL TRACE$POP
      END SUBROUTINE XAS$REPORTSIMULATION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$REPORTSETTINGS  ! MARK: XAS$REPORTSETTINGS
!     **************************************************************************
!     ** REPORT SETTINGS FOR XAS CALCULATION                                  **
!     **************************************************************************
      USE XAS_MODULE, ONLY: SETTINGS,SPECTRUM,SPECTRUMARR
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NFIL
      INTEGER(4) :: ISPEC
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$REPORTSETTINGS')
      CALL FILEHANDLER$UNIT('PROT',NFIL)
      WRITE(NFIL,'(80("#"))')
      WRITE(NFIL,FMT='(A19)')'SETTINGS'
      WRITE(NFIL,'(80("#"))')
      WRITE(NFIL,FMT='(A10,I10)')'NSPECTRA:',SETTINGS%NSPEC
      WRITE(NFIL,FMT='(A10,A)')'HOLE ATOM:',TRIM(SETTINGS%ATOM)
      WRITE(NFIL,FMT='(A10,I10)')'IND ATOM:',SETTINGS%IATOM
      WRITE(NFIL,FMT='(A10,I10)')'NCORE:',SETTINGS%NCORE
      WRITE(NFIL,FMT='(A10,I10)')'LCORE:',SETTINGS%LCORE
      DO ISPEC=1,SETTINGS%NSPEC
        SPECTRUM=>SPECTRUMARR(ISPEC)
        WRITE(NFIL,'(80("#"))')
        WRITE(NFIL,FMT='(A10,I10)')'SPECTRUM:',ISPEC
        WRITE(NFIL,FMT='(A10,3F10.4)')'NORMAL:',SPECTRUM%NORMAL(:)
        WRITE(NFIL,FMT='(A10,3F10.4)')'KDIR:',SPECTRUM%KDIR(:)
        WRITE(NFIL,FMT=-'(A10,2(F8.5,SP,F8.5,"I ",S))')'POL:',SPECTRUM%POL(:)
        WRITE(NFIL,FMT=-'(A10,3(F8.5,SP,F8.5,"I ",S))')'POLXYZ:',SPECTRUM%POLXYZ(:)
      ENDDO
                          CALL TRACE$POP
      END SUBROUTINE XAS$REPORTSETTINGS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$OVERLAP  ! MARK: XAS$OVERLAP
!     **************************************************************************
!     ** CALCULATE OVERLAP MATRIX                                             **
!     **************************************************************************
      USE XAS_MODULE, ONLY: S,OVERLAP,OVERLAPARR,SIM,STATE_TYPE,TOVERLAP
      IMPLICIT NONE
      LOGICAL(4), PARAMETER :: TTEST=.TRUE. ! TEST EIGENVALUES OF OVERLAP
      TYPE(STATE_TYPE), POINTER :: STATE1,STATE2
      INTEGER(4) :: NKPT,NSPIN
      INTEGER(4) :: IKPT,ISPIN
      INTEGER(4) :: NOCC
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$OVERLAP')
      IF(ALLOCATED(S)) THEN
        CALL ERROR$MSG('ATOMIC OVERLAP MATRIX ALREADY CALCULATED')
        CALL ERROR$STOP('XAS$OVERLAP')
      END IF
      IF(.NOT.ALLOCATED(OVERLAPARR)) THEN
        CALL ERROR$MSG('PLANE WAVE OVERLAP MATRIX NOT AVAILABLE')
        CALL ERROR$MSG('SHOULD BE CALCULATED ON READ')
        CALL ERROR$MSG('USE XAS$READ TO READ DATA')
        CALL ERROR$STOP('XAS$OVERLAP')
      END IF
      CALL XAS$OVERLAPATOMMATRIX

      CALL XAS$OVERLAPAUGMENTATION

      NKPT=SIM(1)%NKPT
      NSPIN=SIM(1)%NSPIN
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          OVERLAP=>OVERLAPARR(IKPT,ISPIN)
          STATE1=>SIM(1)%STATEARR(IKPT,ISPIN)
          STATE2=>SIM(2)%STATEARR(IKPT,ISPIN)
          ALLOCATE(OVERLAP%OV(STATE2%NB,STATE1%NB))
          OVERLAP%OV(:,:)=OVERLAP%PW(:,:)+OVERLAP%AUG(:,:)
!         CHECK OF NUMBER OF OCCUPIED STATES IS THE SAME TO PRODUCE SQUARE MATRIX
          IF(STATE1%NOCC.NE.STATE2%NOCC) THEN
            CALL ERROR$MSG('NUMBER OF OCCUPIED STATES NOT THE SAME')
            CALL ERROR$MSG('HAS NUMBER OF ELECTRONS CHANGED OR VARIABLE OCCUPATIONS?')
            CALL ERROR$I4VAL('NOCC1',STATE1%NOCC)
            CALL ERROR$I4VAL('NOCC2',STATE2%NOCC)
            CALL ERROR$I4VAL('IKPT',IKPT)
            CALL ERROR$I4VAL('ISPIN',ISPIN)
            CALL ERROR$STOP('XAS$OVERLAP')
          END IF
          NOCC=STATE1%NOCC
          ALLOCATE(OVERLAP%OVOCC(NOCC,NOCC))
          ALLOCATE(OVERLAP%OVEMP(STATE2%NB-NOCC,NOCC))
          OVERLAP%OVOCC(:,:)=OVERLAP%OV(1:NOCC,1:NOCC)
          OVERLAP%OVEMP(:,:)=OVERLAP%OV(NOCC+1:STATE2%NB,1:NOCC)
        ENDDO
      ENDDO  
      TOVERLAP=.TRUE.  
                          CALL TRACE$POP
      END SUBROUTINE XAS$OVERLAP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$OVERLAPAUGMENTATION  ! MARK: XAS$OVERLAPAUGMENTATION
!     **************************************************************************
!     ** OVERLAP CONTRIBUTION FROM PAW AUGMENTATION                           **
!     **************************************************************************
      USE XAS_MODULE, ONLY: S,OVERLAP,OVERLAPARR,SIM,STATE_TYPE
      IMPLICIT NONE
      TYPE(STATE_TYPE), POINTER :: STATE1,STATE2
      INTEGER(4) :: NKPT,NSPIN
      INTEGER(4) :: IKPT,ISPIN
      INTEGER(4) :: IB1,IB2
      INTEGER(4) :: IPRO1,IPRO2
      COMPLEX(8) :: CVAR
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$OVERLAPAUGMENTATION')
      NKPT=SIM(1)%NKPT
      NSPIN=SIM(1)%NSPIN
  !   LOOP OVER K POINTS
      DO IKPT=1,NKPT
  !     LOOP OVER SPIN
        DO ISPIN=1,NSPIN
          STATE1=>SIM(1)%STATEARR(IKPT,ISPIN)
          STATE2=>SIM(2)%STATEARR(IKPT,ISPIN)
          OVERLAP=>OVERLAPARR(IKPT,ISPIN)
          ALLOCATE(OVERLAP%AUG(STATE2%NB,STATE1%NB))
          OVERLAP%AUG(:,:)=(0.D0,0.D0)
  !       LOOP OVER BANDS OF FIRST SIMULATION
          DO IB2=1,STATE2%NB
  !         LOOP OVER BANDS OF SECOND SIMULATION
            DO IB1=1,STATE1%NB
              CALL XAS$OVERLAPSTATE(STATE1,IB1,STATE2,IB2,CVAR)
!             AUG(I,J)=<PSI1(J)|PSI2(I)>
              OVERLAP%AUG(IB2,IB1)=CVAR
            ENDDO ! END LOOP OVER BANDS OF SECOND SIMULATION
          ENDDO ! END LOOP OVER BANDS OF FIRST SIMULATION
        ENDDO ! END LOOP OVER SPIN
      ENDDO ! END LOOP OVER K POINTS
                          CALL TRACE$POP
      END SUBROUTINE XAS$OVERLAPAUGMENTATION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$OVERLAPSTATE(STATE1,NB1,STATE2,NB2,OVLAP)  ! MARK: XAS$OVERLAPSTATE
!     **************************************************************************
!     ** CALCULATE ATOMIC OVERLAP BETWEEN TWO STATES                          **
!     ** <STATE1|STATE2>                                                      **
!     **************************************************************************
      USE XAS_MODULE, ONLY: STATE_TYPE,S,SIM
      IMPLICIT NONE
      TYPE(STATE_TYPE), INTENT(IN) :: STATE1
      INTEGER(4), INTENT(IN) :: NB1
      TYPE(STATE_TYPE), INTENT(IN) :: STATE2
      INTEGER(4), INTENT(IN) :: NB2
      COMPLEX(8), INTENT(OUT) :: OVLAP
      INTEGER(4) :: IAT
      INTEGER(4) :: ISP1,ISP2
      INTEGER(4) :: LN1,LN2
      INTEGER(4) :: IPRO1,IPRO2
      INTEGER(4) :: L1,L2
      INTEGER(4) :: M1,M2
      COMPLEX(8) :: CVAR
!     **************************************************************************
      OVLAP=(0.D0,0.D0)
      DO IAT=1,SIM(1)%NAT
        ISP1=SIM(1)%ISPECIES(IAT)
! WARNING: REQUIRES THE SAME ORDER OF ATOMS IN BOTH SIMULATIONS
        ISP2=SIM(2)%ISPECIES(IAT)
        DO LN1=1,SIM(1)%LNX(ISP1)
          DO LN2=1,SIM(2)%LNX(ISP2)
            L1=SIM(1)%LOX(LN1,ISP1)
            L2=SIM(2)%LOX(LN2,ISP2)
            IF(L1.NE.L2) CYCLE
            IPRO1=SIM(1)%MAP(IAT,LN1)
            DO M1=1,2*L1+1
              IPRO1=IPRO1+1
              IPRO2=SIM(2)%MAP(IAT,LN2)
              DO M2=1,2*L2+1
                IPRO2=IPRO2+1
                IF(M1.NE.M2) CYCLE
                CVAR=CONJG(STATE1%PROJ(1,NB1,IPRO1))*STATE2%PROJ(1,NB2,IPRO2)
                OVLAP=OVLAP+CVAR*S(IAT,LN1,LN2)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END DO
      END SUBROUTINE XAS$OVERLAPSTATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$OVERLAPATOMMATRIX  ! MARK: XAS$OVERLAPATOMMATRIX
!     **************************************************************************
!     ** CALCULATE ATOMIC OVERLAP MATRIX                                      **
!     ** REAL SPHERICAL HARMONICS GIVE L1=L2 AND M1=M2                        **
!     ** ONLY CALCULATE FOR ONE ARBITRARY M1=M2 AS RESULT IS THE SAME WITHIN  **
!     ** THE SAME L-SHELL; M1=M2 MUST BE ENSURED ELSEWHERE                    **
!     **************************************************************************
! WARNING: REQUIRES THE SAME ORDER OF ATOMS IN BOTH SIMULATIONS
      USE XAS_MODULE, ONLY: SIM,S
      IMPLICIT NONE
      INTEGER(4) :: IAT1,IAT2
      INTEGER(4) :: ISP1,ISP2
      INTEGER(4) :: LN1,LN2
      INTEGER(4) :: L1,L2
      REAL(8) :: SVAL
      INTEGER(4),ALLOCATABLE :: IATMAP1(:),IATMAP2(:)
      REAL(8), ALLOCATABLE :: PSPHI1(:),PSPHI2(:)
      REAL(8), ALLOCATABLE :: AEPHI1(:),AEPHI2(:)
      REAL(8), ALLOCATABLE :: AUX(:)
      REAL(8), ALLOCATABLE :: R(:)
      INTEGER(4) :: NR1,NR2
      INTEGER(4) :: GID1,GID2
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$OVERLAPATOMMATRIX')
      ALLOCATE(S(SIM(1)%NAT,SIM(1)%LNXX,SIM(2)%LNXX))
      S=0.D0
      DO IAT1=1,SIM(1)%NAT
        ISP1=SIM(1)%ISPECIES(IAT1)
        GID1=SIM(1)%SETUP(ISP1)%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        ALLOCATE(PSPHI1(NR1))
        ALLOCATE(AEPHI1(NR1))
        ALLOCATE(PSPHI2(NR1))
        ALLOCATE(AEPHI2(NR1))
        ALLOCATE(AUX(NR1))
        ALLOCATE(R(NR1))
        CALL RADIAL$R(GID1,NR1,R)

        DO LN1=1,SIM(1)%LNX(ISP1)
          L1=SIM(1)%LOX(LN1,ISP1)
          IAT2=IAT1
          ISP2=SIM(2)%ISPECIES(IAT2)
          GID2=SIM(2)%SETUP(ISP2)%GID
          CALL RADIAL$GETI4(GID2,'NR',NR2)

          DO LN2=1,SIM(2)%LNX(ISP2)
            L2=SIM(2)%LOX(LN2,ISP2)
            IF(L1.NE.L2) CYCLE
!           GET RADIAL FUNCTIONS
            PSPHI1=SIM(1)%SETUP(ISP1)%PSPHI(:,LN1)
            AEPHI1=SIM(1)%SETUP(ISP1)%AEPHI(:,LN1)  
!           IF GRIDS ARE DIFFERENT, MAP ONTO FIRST GRID   
            IF(GID1.NE.GID2) THEN
              CALL RADIAL$CHANGEGRID(GID2,NR2,SIM(2)%SETUP(ISP2)%PSPHI(:,LN2), &
     &                               GID1,NR1,PSPHI2)
              CALL RADIAL$CHANGEGRID(GID2,NR2,SIM(2)%SETUP(ISP2)%AEPHI(:,LN2), &
     &                               GID1,NR1,AEPHI2)
            ELSE
              PSPHI2=SIM(2)%SETUP(ISP2)%PSPHI(:,LN2)
              AEPHI2=SIM(2)%SETUP(ISP2)%AEPHI(:,LN2)
            END IF
! TODO: CHECK IF ORDER IS CORRECT
            AUX=AEPHI1*AEPHI2-PSPHI1*PSPHI2
            AUX=AUX*R**2
! TODO: CHECK IF AUX IS ZERO OUTSIDE OF R_AUG
            CALL RADIAL$INTEGRAL(GID1,NR1,AUX,SVAL)
            S(IAT1,LN1,LN2)=SVAL
          ENDDO
        ENDDO
        DEALLOCATE(PSPHI1)
        DEALLOCATE(AEPHI1)
        DEALLOCATE(PSPHI2)
        DEALLOCATE(AEPHI2)
        DEALLOCATE(AUX)
        DEALLOCATE(R)
      ENDDO
      ! IPRO=0
      ! WRITE(*,FMT='(6A10)')'IAT','ISP','LN','L','M','IPRO'
      ! DO IAT1=1,SIM(1)%NAT
      !   ISP1=SIM(1)%ISPECIES(IAT1)
      !   DO LN1=1,SIM(1)%LNX(ISP1)
      !     L1=SIM(1)%LOX(LN1,ISP1)
      !     DO M1=1,2*L1+1
      !       IPRO=IPRO+1
      !       WRITE(*,FMT='(6I10)')IAT1,ISP1,LN1,L1,M1,IPRO
      !     ENDDO
      !   ENDDO
      ! ENDDO
                          CALL TRACE$POP
      END SUBROUTINE XAS$OVERLAPATOMMATRIX
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$POLARISATION  ! MARK: XAS$POLARISATION
!     **************************************************************************
!     **  CALCULATE POLARISATION FROM K VECTOR, SURFACE NORMAL,               **
!     **  AND POLARISATION OF INCIDENT LIGHT                                  **
!     **************************************************************************
      USE XAS_MODULE, ONLY: SPECTRUM,SPECTRUMARR,SETTINGS
      IMPLICIT NONE
      REAL(8), PARAMETER :: TOL=1.D-10
      INTEGER(4) :: ISPEC
      REAL(8) :: KVEC(3)
      REAL(8) :: WORK(3)
      REAL(8) :: WORK2(3)
      REAL(8) :: SVAR
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$POLARISATION')
      IF(.NOT.ALLOCATED(SPECTRUMARR)) THEN
        CALL ERROR$MSG('NO SPECTRUM AVAILABLE')
        CALL ERROR$STOP('XAS$POLARISATION')
      END IF
!     LOOP OVER SPECTRA
      DO ISPEC=1,SETTINGS%NSPEC
        SPECTRUM=>SPECTRUMARR(ISPEC)
        KVEC=SPECTRUM%KDIR
        CALL CROSS_PROD(KVEC,SPECTRUM%NORMAL,WORK)
        SVAR=NORM2(WORK)
        IF(SVAR.LT.TOL) THEN
          CALL ERROR$MSG('K VECTOR AND SURFACE NORMAL ARE PARALLEL')
          CALL ERROR$MSG('POLARISATION VECTOR SET TO ARBITRARY ORTHOGONAL VECTOR')
          CALL VEC_ORTHO(KVEC,WORK)
          SVAR=NORM2(WORK)
        END IF
        WORK=WORK/SVAR
        CALL CROSS_PROD(WORK,KVEC,WORK2)
        WORK2=WORK2/NORM2(WORK2)
        SPECTRUM%POLXYZ=SPECTRUM%POL(1)*WORK+ &
      &                 SPECTRUM%POL(2)*WORK2
        SVAR=SQRT(SUM(ABS(SPECTRUM%POLXYZ)**2))
        SPECTRUM%POLXYZ=SPECTRUM%POLXYZ/SVAR
      ENDDO
                          CALL TRACE$POP
      END SUBROUTINE XAS$POLARISATION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$DIPOLEMATRIX ! MARK: XAS$DIPOLEMATRIX
!     **************************************************************************
!     ** CALCULATE DIPOLE MATRIX ELEMENTS IN EXCITED ORBITAL BASIS            **
!     ** SIMULATION 2 IS THE EXCITED STATE, SIMULATION 1 IS THE GROUND STATE  **
!     **************************************************************************
      USE XAS_MODULE, ONLY: SETTINGS,SIM,SETUP_TYPE,STATE_TYPE,OVERLAP,OVERLAPARR
      IMPLICIT NONE
      LOGICAL(4), PARAMETER :: TTEST=.FALSE.
      REAL(8), PARAMETER :: PI=4.D0*ATAN(1.D0)
      INTEGER(4) :: NFIL
      INTEGER(4) :: IATOM ! ATOM WITH CORE HOLE
      INTEGER(4) :: ISP   ! SETUP INDEX OF ATOM WITH CORE HOLE
      TYPE(SETUP_TYPE), POINTER :: STP ! POINTER TO SETUP
      TYPE(STATE_TYPE), POINTER :: STATE
      INTEGER(4) :: N
      INTEGER(4) :: IB
      INTEGER(4) :: LN
      INTEGER(4) :: NR
      REAL(8), ALLOCATABLE :: AEPSI(:) ! (NR) ATOMIC RADIAL PART
      REAL(8) :: GAUNT(3)
      INTEGER(4) :: LLCORE ! LL OF CORE HOLE
      INTEGER(4) :: LLVAL ! LL OF VALENCE ORBITAL
      INTEGER(4) :: L
      INTEGER(4) :: M
      INTEGER(4) :: IKPT
      INTEGER(4) :: ISPIN
      INTEGER(4) :: IPRO
      REAL(8), ALLOCATABLE :: RADINT(:) ! (LNX(ISP)) RADIAL INTEGRAL VALUES
      REAL(8), ALLOCATABLE :: R(:)
      REAL(8), ALLOCATABLE :: WORK(:)
      COMPLEX(8) :: CVAR(3)
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$DIPOLEMATRIX')
      IF(TTEST) THEN
        CALL FILEHANDLER$UNIT('PROT',NFIL)
        WRITE(NFIL,'(80("#"))')
        WRITE(NFIL,FMT='(A)')'DIPOLE MATRIX CALCULATION'
        WRITE(NFIL,'(80("#"))')
      ENDIF
      IATOM=SETTINGS%IATOM
      ISP=SIM(2)%ISPECIES(IATOM)
      STP=>SIM(2)%SETUP(ISP)
      CALL RADIAL$GETI4(STP%GID,'NR',NR)
      ALLOCATE(AEPSI(NR))
      ALLOCATE(R(NR))
      CALL RADIAL$R(STP%GID,NR,R)
      ALLOCATE(WORK(NR))
!     ==========================================================================
!     == SELECT CORRECT RADIAL PART FOR CORE ORBITAL                          ==
!     == NOTE: THIS ASSUME THE FOLLOWING STRUCTURE                            ==
!     ==       |IB | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |10 |                  ==
!     ==       |---|---|---|---|---|---|---|---|---|---|---|                  ==
!     ==       | N | 1 | 2 | 2 | 3 | 3 | 3 | 4 | 4 | 4 | 4 |                  ==
!     ==       | L | 0 | 0 | 1 | 0 | 1 | 2 | 0 | 1 | 2 | 3 |                  ==
!     ==========================================================================
      N=0
      DO IB=1,STP%NBATOM
        IF(STP%LATOM(IB).EQ.0) N=N+1
        IF(N.EQ.SETTINGS%NCORE.AND.STP%LATOM(IB).EQ.SETTINGS%LCORE) THEN
          AEPSI(:)=STP%AEPSI(:,IB)
          IF(TTEST) WRITE(NFIL,FMT='(A10,I10)')'INDEX RAD:',IB
        ENDIF
      ENDDO
! WARING: THIS ASSUMES THAT THE CORE ORBITAL IS AN S ORBITAL
! TODO: GENERALISE TO ARBITRARY ORBITALS
      CALL LLOFLM(SETTINGS%LCORE,0,LLCORE)
      IF(TTEST) WRITE(NFIL,FMT='(A10,I10)')'LLCORE:',LLCORE
!     PRE-CALCULATE RADIAL INTEGRALS
      ALLOCATE(RADINT(SIM(2)%LNX(ISP)))
      DO LN=1,SIM(2)%LNX(ISP)
        WORK=STP%AEPHI(:,LN)*R**3*AEPSI
        CALL RADIAL$INTEGRAL(STP%GID,NR,WORK,RADINT(LN))
      ENDDO
!     LOOP OVER K POINTS
      DO IKPT=1,SIM(2)%NKPT
!       LOOP OVER SPIN
        DO ISPIN=1,SIM(2)%NSPIN
          STATE=>SIM(2)%STATEARR(IKPT,ISPIN)
          OVERLAP=>OVERLAPARR(IKPT,ISPIN)
          ALLOCATE(OVERLAP%DIPOLE(3,STATE%NB))
          IF(TTEST) WRITE(NFIL,FMT='(A10,I10,A10,I10)')'IKPT:',IKPT,'ISPIN:',ISPIN
!         LOOP OVER BANDS
          DO IB=1,STATE%NB
            IF(TTEST) WRITE(NFIL,FMT='(A10,I10)')'BAND:',IB
            CVAR=(0.D0,0.D0)
            DO LN=1,SIM(2)%LNX(ISP)
              L=SIM(2)%LOX(LN,ISP)
              IPRO=SIM(2)%MAP(IATOM,LN)
              DO M=-L,L
                IPRO=IPRO+1
                CALL LLOFLM(L,M,LLVAL)
!               CALCULATE GAUNT COEFFICIENT VECTOR
!               X COMPONENT
                CALL SPHERICAL$GAUNT(LLVAL,2,LLCORE,GAUNT(1))
!               Y COMPONENT
                CALL SPHERICAL$GAUNT(LLVAL,4,LLCORE,GAUNT(2))
!               Z COMPONENT
                CALL SPHERICAL$GAUNT(LLVAL,3,LLCORE,GAUNT(3))
                CVAR=CVAR+CONJG(STATE%PROJ(1,IB,IPRO))*RADINT(LN)*GAUNT
              ENDDO
            ENDDO
            CVAR=SQRT(4.D0*PI/3.D0)*CVAR
            IF(TTEST) WRITE(NFIL,FMT='(A10,3(F8.5,SP,F8.5,"I ",S))')'DIPOLE:',CVAR(:)
            OVERLAP%DIPOLE(:,IB)=CVAR
          ENDDO ! END LOOP OVER BANDS
        ENDDO ! END LOOP OVER SPINS
      ENDDO ! END LOOP OVER K POINTS          
      DEALLOCATE(AEPSI)
      DEALLOCATE(R)
      DEALLOCATE(WORK)
      DEALLOCATE(RADINT)
                          CALL TRACE$POP
      END SUBROUTINE XAS$DIPOLEMATRIX
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$ADET ! MARK: XAS$ADET
!     **************************************************************************
!     ** CALCULATE DETERMINANT ADET FOR OCCUPIED OVERLAP MATRIX               **
!     **************************************************************************
      USE XAS_MODULE, ONLY: OVERLAP,OVERLAPARR,SIM
      IMPLICIT NONE
      INTEGER(4) :: IKPT
      INTEGER(4) :: ISPIN
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$ADET')
      DO IKPT=1,SIM(1)%NKPT
        DO ISPIN=1,SIM(1)%NSPIN
          OVERLAP=>OVERLAPARR(IKPT,ISPIN)
          SIM(1)%STATE=>SIM(1)%STATEARR(IKPT,ISPIN)
          CALL LIB$DETC8(SIM(1)%STATE%NOCC,OVERLAP%OVOCC,OVERLAP%ADET)
          WRITE(*,*)'ADET:',OVERLAP%ADET
        ENDDO
      ENDDO
                          CALL TRACE$POP
      END SUBROUTINE XAS$ADET
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS$KMAT  ! MARK: XAS$KMAT
!     **************************************************************************
!     ** CALCULATE K MATRIX FOR XAS CALCULATION                               **
!     **************************************************************************
      USE XAS_MODULE, ONLY: OVERLAP,OVERLAPARR,SIM
      IMPLICIT NONE
      INTEGER(4) :: IKPT
      INTEGER(4) :: ISPIN
      INTEGER(4) :: NOCC
      INTEGER(4) :: NB1,NB2
      COMPLEX(8), ALLOCATABLE :: WORK(:,:)
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$KMAT')
      DO IKPT=1,SIM(1)%NKPT
        DO ISPIN=1,SIM(1)%NSPIN
          OVERLAP=>OVERLAPARR(IKPT,ISPIN)
          SIM(1)%STATE=>SIM(1)%STATEARR(IKPT,ISPIN)
          SIM(2)%STATE=>SIM(2)%STATEARR(IKPT,ISPIN)
          NOCC=SIM(1)%STATE%NOCC
          NB1=SIM(1)%STATE%NB
          NB2=SIM(2)%STATE%NB
          ALLOCATE(OVERLAP%KMAT(NB2-NOCC,NB1))
          ALLOCATE(WORK(NOCC,NOCC))
          CALL LIB$INVERTC8(NOCC,OVERLAP%OVOCC,WORK)
          CALL LIB$MATMULC8(NB2-NOCC,NOCC,NOCC,OVERLAP%OVEMP,WORK,OVERLAP%KMAT)
          DEALLOCATE(WORK)
        ENDDO
      ENDDO
                          CALL TRACE$POP
      END SUBROUTINE XAS$KMAT
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LLOFLM(L,M,LL) ! MARK: LLOFLM
!     **************************************************************************
!     ** CALCULATE LL=L*L+L-M+1                                               **
!     **************************************************************************
      INTEGER(4), INTENT(IN) :: L
      INTEGER(4), INTENT(IN) :: M
      INTEGER(4), INTENT(OUT) :: LL
!     **************************************************************************
      LL=L*L+L-M+1
      RETURN
      END SUBROUTINE LLOFLM
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CROSS_PROD(A,B,C)  ! MARK: CROSS_PROD
!     **************************************************************************
!     ** CALCULATE CROSS PRODUCT OF TWO VECTORS                               **
!     **************************************************************************
      REAL(8), INTENT(IN) :: A(3)
      REAL(8), INTENT(IN) :: B(3)
      REAL(8), INTENT(OUT) :: C(3)
!     **************************************************************************
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END SUBROUTINE CROSS_PROD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VEC_ORTHO(A,B)  ! MARK: VEC_ORTHO
!     **************************************************************************
!     ** CALCULATE ORTHOGONAL VECTOR TO A                                      **
!     **************************************************************************
! TODO: UNDERSTAND
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: A(3)
      REAL(8), INTENT(OUT) :: B(3)
      REAL(8) :: VECVAR(3)
!     **************************************************************************
      VECVAR = (/1.D0,0.D0,0.D0/)
      IF(DOT_PRODUCT(A,VECVAR).EQ.NORM2(A)) THEN
        VECVAR = (/0.D0,1.D0,0.D0/)
      ENDIF
      CALL CROSS_PROD(A,VECVAR,B)
      END SUBROUTINE VEC_ORTHO