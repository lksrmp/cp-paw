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
        COMPLEX(8), POINTER :: VEC(:,:,:)  
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
        REAL(8) :: NKDIV(3) ! K-POINT DIVISIONS
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
      CALL XASCNTL$FILES('GROUNDSTATE')
      CALL XASCNTL$FILES('EXCITESTATE')

!     READ FIRST FILE
      CALL XAS$SELECT('GROUNDSTATE')
      CALL XAS$READ
      CALL XAS$UNSELECT

!     READ SECOND FILE
      CALL XAS$SELECT('EXCITESTATE')
      CALL XAS$READ
      CALL XAS$UNSELECT

!     REPORT UNUSED LINKEDLISTS
      CALL FILEHANDLER$UNIT('PROT',NFIL)
      CALL LINKEDLIST$REPORT_UNUSED(LL_CNTL,NFIL)
      CALL RADIAL$REPORT(NFIL)

     
      CALL XAS$ISELECT(1)
      OPEN(UNIT=NFIL,FILE='gs_1_psphi.dat',STATUS='REPLACE')
      CALL XAS$WRITEPHI(NFIL,1,'PSPHI')
      CLOSE(NFIL)
      OPEN(UNIT=NFIL,FILE='gs_2_psphi.dat',STATUS='REPLACE')
      CALL XAS$WRITEPHI(NFIL,2,'PSPHI')
      CLOSE(NFIL)
      OPEN(UNIT=NFIL,FILE='gs_1_aephi.dat',STATUS='REPLACE')
      CALL XAS$WRITEPHI(NFIL,1,'AEPHI')
      CLOSE(NFIL)
      OPEN(UNIT=NFIL,FILE='gs_2_aephi.dat',STATUS='REPLACE')
      CALL XAS$WRITEPHI(NFIL,2,'AEPHI')
      CLOSE(NFIL)
      CALL XAS$UNSELECT

      CALL XAS$ISELECT(2)
      OPEN(UNIT=NFIL,FILE='exc_1_psphi.dat',STATUS='REPLACE')
      CALL XAS$WRITEPHI(NFIL,1,'PSPHI')
      CLOSE(NFIL)
      OPEN(UNIT=NFIL,FILE='exc_2_psphi.dat',STATUS='REPLACE')
      CALL XAS$WRITEPHI(NFIL,2,'PSPHI')
      CLOSE(NFIL)
      OPEN(UNIT=NFIL,FILE='exc_1_aephi.dat',STATUS='REPLACE')
      CALL XAS$WRITEPHI(NFIL,1,'AEPHI')
      CLOSE(NFIL)
      OPEN(UNIT=NFIL,FILE='exc_2_aephi.dat',STATUS='REPLACE')
      CALL XAS$WRITEPHI(NFIL,2,'AEPHI')
      CLOSE(NFIL)
      CALL XAS$UNSELECT

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

! !     ..................................................................
!       SUBROUTINE DATACONSISTENCY  ! MARK: DATACONSISTENCY
! !     **************************************************************************
! !     ** CHECKS CONSISTENCY OF PDOS DATA                                      **
! ! TODO: IMPLEMENT THIS SUBROUTINE FURTHER, TEST THE CHECKS, ADD ADDITIONAL CHECKS
! !     **************************************************************************
!       USE XASPDOS_MODULE, ONLY: PD,NPD
!       IMPLICIT NONE
!       INTEGER(4) :: IPD
!       INTEGER(4) :: I,J
!       DO IPD=2,NPD
!         IF(PD(IPD)%FLAG.NE.PD(1)%FLAG)THEN
!           CALL ERROR$MSG('FLAG INCONSISTENT BETWEEN PDOS FILES')
!           CALL ERROR$I4VAL('IPD',IPD)
!           CALL ERROR$CHVAL('PD(IPD)%FLAG',PD(IPD)%FLAG)
!           CALL ERROR$CHVAL('PD(1)%FLAG',PD(1)%FLAG)
!           CALL ERROR$STOP('DATACONSISTENCY')
!         END IF
!         IF(PD(IPD)%NAT.NE.PD(1)%NAT)THEN
!           CALL ERROR$MSG('NAT INCONSISTENT BETWEEN PDOS FILES')
!           CALL ERROR$I4VAL('IPD',IPD)
!           CALL ERROR$I4VAL('PD(IPD)%NAT',PD(IPD)%NAT)
!           CALL ERROR$I4VAL('PD(1)%NAT',PD(1)%NAT)
!           CALL ERROR$STOP('DATACONSISTENCY')
!         END IF
!         IF(PD(IPD)%NKPT.NE.PD(1)%NKPT)THEN
!           CALL ERROR$MSG('NKPT INCONSISTENT BETWEEN PDOS FILES')
!           CALL ERROR$I4VAL('IPD',IPD)
!           CALL ERROR$I4VAL('PD(IPD)%NKPT',PD(IPD)%NKPT)
!           CALL ERROR$I4VAL('PD(1)%NKPT',PD(1)%NKPT)
!           CALL ERROR$STOP('DATACONSISTENCY')
!         END IF
!         IF(PD(IPD)%NSPIN.NE.PD(1)%NSPIN)THEN
!           CALL ERROR$MSG('NSPIN INCONSISTENT BETWEEN PDOS FILES')
!           CALL ERROR$I4VAL('IPD',IPD)
!           CALL ERROR$I4VAL('PD(IPD)%NSPIN',PD(IPD)%NSPIN)
!           CALL ERROR$I4VAL('PD(1)%NSPIN',PD(1)%NSPIN)
!           CALL ERROR$STOP('DATACONSISTENCY')
!         END IF
!         IF(PD(IPD)%NDIM.NE.PD(1)%NDIM)THEN
!           CALL ERROR$MSG('NDIM INCONSISTENT BETWEEN PDOS FILES')
!           CALL ERROR$I4VAL('IPD',IPD)
!           CALL ERROR$I4VAL('PD(IPD)%NDIM',PD(IPD)%NDIM)
!           CALL ERROR$I4VAL('PD(1)%NDIM',PD(1)%NDIM)
!           CALL ERROR$STOP('DATACONSISTENCY')
!         END IF
!         DO I=1,3
!           IF(PD(IPD)%NKDIV(I).NE.PD(1)%NKDIV(I))THEN
!             CALL ERROR$MSG('NKDIV INCONSISTENT BETWEEN PDOS FILES')
!             CALL ERROR$I4VAL('IPD',IPD)
!             CALL ERROR$I4VAL('PD(IPD)%NKDIV',PD(IPD)%NKDIV(I))
!             CALL ERROR$I4VAL('PD(1)%NKDIV',PD(1)%NKDIV(I))
!             CALL ERROR$STOP('DATACONSISTENCY')
!           END IF
!         ENDDO
!         DO I=1,3
!           DO J=1,3
!             IF(PD(IPD)%RBAS(I,J).NE.PD(1)%RBAS(I,J))THEN
!               CALL ERROR$MSG('RBAS INCONSISTENT BETWEEN PDOS FILES')
!               CALL ERROR$I4VAL('IPD',IPD)
!               CALL ERROR$R8VAL('PD(IPD)%RBAS',PD(IPD)%RBAS(I,J))
!               CALL ERROR$R8VAL('PD(1)%RBAS',PD(1)%RBAS(I,J))
!               CALL ERROR$STOP('DATACONSISTENCY')
!             END IF
!           ENDDO
!         ENDDO
!       ENDDO
!       END SUBROUTINE DATACONSISTENCY

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
! !
! !     ...1.........2.........3.........4.........5.........6.........7.........8
!       SUBROUTINE XASCNTL$FILES  ! MARK: XASCNTL$FILES
! !     **************************************************************************
!       USE XASCNTL_MODULE, ONLY: LL_CNTL
!       USE LINKEDLIST_MODULE
!       USE STRINGS_MODULE
!       IMPLICIT NONE
!       LOGICAL(4) :: TCHK
!       LOGICAL(4) :: PDOS1
!       LOGICAL(4) :: PDOS2
!       INTEGER(4) :: NUM
!       INTEGER(4) :: I
!       CHARACTER(32) :: ID
!       CHARACTER(256) :: FILENAME
!       LOGICAL(4) :: EXT
! !     **************************************************************************
!                           CALL TRACE$PUSH('XASCNTL$FILES')
!       CALL LINKEDLIST$SELECT(LL_CNTL,'~')
!       CALL LINKEDLIST$SELECT(LL_CNTL,'XCNTL')
!       CALL LINKEDLIST$EXISTL(LL_CNTL,'FILES',1,TCHK)
!       IF(.NOT.TCHK) THEN
!         CALL ERROR$MSG('!FILES NOT FOUND IN XCNTL')
!         CALL ERROR$STOP('XASCNTL$FILES')
!       END IF
!       CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
!       CALL LINKEDLIST$NLISTS(LL_CNTL,'FILE',NUM)
! print *,'NUM',NUM
!       IF(NUM.EQ.0) THEN
!         CALL ERROR$MSG('NO FILES DEFINED IN XCNTL')
!         CALL ERROR$STOP('XASCNTL$FILES')
!       END IF
!       DO I=1,NUM
! print *,'I',I
!         CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',I)
!         CALL LINKEDLIST$GET(LL_CNTL,'ID',0,ID)
!         ID=+ID
!         CALL LINKEDLIST$GET(LL_CNTL,'NAME',0,FILENAME)
!         CALL LINKEDLIST$GET(LL_CNTL,'EXT',0,EXT)
!         CALL FILEHANDLER$SETFILE(ID,EXT,FILENAME)
!         IF(ID.EQ.'PDOS1') THEN
!           PDOS1=.TRUE.
!         ELSE IF(ID.EQ.'PDOS2') THEN
!           PDOS2=.TRUE.
!         END IF
!         CALL LINKEDLIST$SELECT(LL_CNTL,'..')
!       ENDDO
!       IF((.NOT.PDOS1).AND.(.NOT.PDOS2)) THEN
!         CALL ERROR$MSG('PDOS FILES NOT DEFINED IN !XCNTL!FILES')
!         CALL ERROR$STOP('XASCNTL$FILES')
!       ELSE IF(.NOT.PDOS1) THEN
!         CALL ERROR$MSG('PDOS1 FILE NOT DEFINED IN !XCNTL!FILES')
!         CALL ERROR$STOP('XASCNTL$FILES')
!       ELSE IF(.NOT.PDOS2) THEN
!         CALL ERROR$MSG('PDOS2 FILE NOT DEFINED IN !XCNTL!FILES')
!         CALL ERROR$STOP('XASCNTL$FILES')
!       END IF
!                           CALL TRACE$POP
!       END SUBROUTINE XASCNTL$FILES

      SUBROUTINE XASCNTL$FILES(FLAG)  ! MARK: XASCNTL$FILES
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE XASCNTL_MODULE, ONLY: LL_CNTL
      USE XAS_MODULE, ONLY: THIS,SELECTED
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      CHARACTER(11), INTENT(IN) :: FLAG ! GROUNDSTATE OR EXCITESTATE
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
      IF(FLAG.NE.'GROUNDSTATE'.AND.FLAG.NE.'EXCITESTATE') THEN
        CALL ERROR$MSG('FLAG NOT RECOGNIZED')
        CALL ERROR$CHVAL('FLAG',FLAG)
        CALL ERROR$MSG('FLAG MUST BE GROUNDSTATE OR EXCITESTATE')
        CALL ERROR$STOP('XASCNTL$FILES')
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'XCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,FLAG,1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!'//FLAG//' NOT FOUND IN XCNTL')
        CALL ERROR$STOP('XASCNTL$FILES')
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,FLAG)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('FILE NOT FOUND IN !'//FLAG)
        CALL ERROR$STOP('XASCNTL$FILES')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'FILE',0,FILENAME)
    
      CALL XAS$SETCH('FILE',FILENAME)
      CALL XAS$SETCH('FLAG',FLAG)
      CALL FILEHANDLER$SETFILE(FLAG,.FALSE.,TRIM(FILENAME))
                          CALL TRACE$POP
      END SUBROUTINE XASCNTL$FILES




! !
! !     ...1.........2.........3.........4.........5.........6.........7.........8
!       SUBROUTINE XASSTP$NEW(NFIL_,ISTP,TCHK)  ! MARK: XASSTP$NEW
! !     **************************************************************************
! !     **                                                                      **
! !     **************************************************************************
!       USE XASSTP_MODULE, ONLY: STP
!       USE LINKEDLIST_MODULE
!       IMPLICIT NONE
!       INTEGER(4)  ,INTENT(IN) :: NFIL_
!       INTEGER(4) ,INTENT(IN) :: ISTP
!       LOGICAL(4)  ,INTENT(OUT):: TCHK
!       TYPE(LL_TYPE)           :: LL_STP_
!       CHARACTER(6)            :: STRING
!       CHARACTER(16)           :: TYPE
!       REAL(8)                 :: R1
!       REAL(8)                 :: DEX
!       INTEGER(4)              :: LNX
! !     **********************************************************************
!       LL_STP_=STP(ISTP)%LL
!       TCHK=.FALSE.
! !     == FIRST REMOVE OLD LIST IF NEEDED
!       IF(STP(ISTP)%NFIL.NE.0) THEN
!         CALL LINKEDLIST$SELECT(LL_STP_,'~')
! !       CALL LINKEDLIST$DELETE(LL_STP_)
!         CALL LINKEDLIST$RMLIST(LL_STP_,'SETUP')
!         STP(ISTP)%GID=0
!         STP(ISTP)%LNX=0
!         STP(ISTP)%NR=0
!         STP(ISTP)%NFIL=0
!       END IF
!       IF(NFIL_.EQ.0) RETURN
! !
! !     ======================================================================  
! !     == CHECK IF THE FILE IS A LINKEDLIST                                ==  
! !     ======================================================================  
! !PRINT*,'CHECK'
!       TCHK=.TRUE.
!       REWIND NFIL_
!       READ(NFIL_,*)STRING
!       REWIND NFIL_
!       IF(STRING.NE.'!SETUP') THEN
!         TCHK=.FALSE.  
!         RETURN
!       END IF
!       STP(ISTP)%NFIL=NFIL_
! !
! !     ======================================================================  
! !     == CREATE NEW LINKED LIST AND READ                                  ==  
! !     ======================================================================  
! !PRINT*,'CREATE LINKEDLIST',STRING
!       CALL LINKEDLIST$NEW(LL_STP_)
! !PRINT*,'READ',NFIL
!       CALL LINKEDLIST$READ(LL_STP_,STP(ISTP)%NFIL,'MONOMER')
!       TCHK=.TRUE.
! !
! !     ======================================================================  
! !     == READ SOME DATA                                                   ==  
! !     ======================================================================  
! !     == LNX
!       CALL LINKEDLIST$SELECT(LL_STP_,'~')
!       CALL LINKEDLIST$SELECT(LL_STP_,'SETUP')
!       CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!       CALL LINKEDLIST$GET(LL_STP_,'LNX',0,LNX)
!       STP(ISTP)%LNX=LNX
! !     == GRID
!       CALL LINKEDLIST$SELECT(LL_STP_,'~')
!       CALL LINKEDLIST$SELECT(LL_STP_,'SETUP')
!       CALL LINKEDLIST$SELECT(LL_STP_,'GRID')
!       CALL LINKEDLIST$GET(LL_STP_,'TYPE',0,TYPE)
!       CALL LINKEDLIST$GET(LL_STP_,'R1',0,R1)
!       CALL LINKEDLIST$GET(LL_STP_,'DEX',0,DEX)
!       CALL LINKEDLIST$GET(LL_STP_,'NR',0,STP(ISTP)%NR)
!       CALL RADIAL$NEW(TRIM(TYPE),STP(ISTP)%GID)
!       CALL RADIAL$SETR8(STP(ISTP)%GID,'R1',R1)
!       CALL RADIAL$SETR8(STP(ISTP)%GID,'DEX',DEX)
!       CALL RADIAL$SETI4(STP(ISTP)%GID,'NR',STP(ISTP)%NR)
! !     == NB,NC
!       CALL LINKEDLIST$SELECT(LL_STP_,'~')
!       CALL LINKEDLIST$SELECT(LL_STP_,'SETUP')
!       CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!       CALL LINKEDLIST$GET(LL_STP_,'NB',0,STP(ISTP)%NB)
!       CALL LINKEDLIST$GET(LL_STP_,'NC',0,STP(ISTP)%NC)
!       STP(ISTP)%LL=LL_STP_
!       RETURN
!       END

! !
! !     ...1.........2.........3.........4.........5.........6.........7.........8
!       SUBROUTINE XASSTP$GETI4(ID,ISTP,VAL)  ! MARK: XASSTP$GETI4
! !     **************************************************************************
! !     **                                                                      **
! !     **************************************************************************
!       USE XASSTP_MODULE, ONLY: STP
!       USE LINKEDLIST_MODULE
!       IMPLICIT NONE
!       CHARACTER(*),INTENT(IN) :: ID
!       INTEGER(4)  ,INTENT(IN) :: ISTP
!       INTEGER(4)  ,INTENT(OUT):: VAL
!       TYPE(LL_TYPE)           :: LL_STP_
! !     **************************************************************************
!       LL_STP_=STP(ISTP)%LL
!       CALL LINKEDLIST$SELECT(LL_STP_,'~')
!       CALL LINKEDLIST$SELECT(LL_STP_,'SETUP')
! !
! !     ==========================================================================
! !     == #(PARTIAL WAVES)                                                     ==
! !     ==========================================================================
!       IF(ID.EQ.'LNX') THEN
!         VAL=STP(ISTP)%LNX
! !
! !     ==========================================================================
! !     == #(RADIAL GRID POINTS)                                                ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'NR') THEN
!         VAL=STP(ISTP)%NR
! !
! !     ==========================================================================
! !     == GRID ID                                                              ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'GID') THEN
!         VAL=STP(ISTP)%GID
! !
! !     ==========================================================================
! !     == NUMBER OF CORE STATES                                                ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'NC') THEN
!         CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!         CALL LINKEDLIST$GET(LL_STP_,'NC',0,VAL)
! !
! !     ==========================================================================
! !     == NUMBER OF ATOMIC STATES AVAILABLE                                    ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'NB') THEN
!         CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!         CALL LINKEDLIST$GET(LL_STP_,'NB',0,VAL)
! !
! !     ==========================================================================
! !     ==  WRONG ID                                                            ==
! !     ==========================================================================
!       ELSE
!         CALL ERROR$MSG('ID NOT RECOGNIZED')
!         CALL ERROR$CHVAL('ID',ID)
!         CALL ERROR$STOP('XASSTP$GETI4')
!       END IF
!       RETURN
!       END

! !
! !     ...1.........2.........3.........4.........5.........6.........7.........8
!       SUBROUTINE XASSTP$GETR8(ID,ISTP,VAL)  ! MARK: XASSTP$GETR8
! !     **************************************************************************
! !     **                                                                      **
! !     **************************************************************************
!       USE XASSTP_MODULE, ONLY: STP
!       USE LINKEDLIST_MODULE
!       IMPLICIT NONE
!       CHARACTER(*),INTENT(IN) :: ID
!       INTEGER(4)  ,INTENT(IN) :: ISTP
!       REAL(8)     ,INTENT(OUT):: VAL
!       TYPE(LL_TYPE)           :: LL_STP_
! !     **********************************************************************
!       LL_STP_=STP(ISTP)%LL
!       CALL LINKEDLIST$SELECT(LL_STP_,'~')
!       CALL LINKEDLIST$SELECT(LL_STP_,'SETUP')
!       IF(ID.EQ.'AEZ') THEN
!         CALL LINKEDLIST$SELECT(LL_STP_,'GENERIC')
!         CALL LINKEDLIST$GET(LL_STP_,'AEZ',0,VAL)
!       ELSE IF(ID.EQ.'PSZ') THEN
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'PSZ',0,VAL)
!       ELSE IF(ID.EQ.'RCSM') THEN
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'RCSM',0,VAL)
!       ELSE
!         CALL ERROR$MSG('ID NOT RECOGNIZED')
!         CALL ERROR$CHVAL('ID',ID)
!         CALL ERROR$STOP('XASSTP$GETR8')
!       END IF
!       RETURN
!       END
! !
! !     ...1.........2.........3.........4.........5.........6.........7.........8
!       SUBROUTINE XASSTP$GETI4A(ID,ISTP,LENG,VAL)  ! MARK: XASSTP$GETI4A
! !     **************************************************************************
! !     **                                                                      **
! !     **************************************************************************
!       USE XASSTP_MODULE, ONLY: STP
!       USE LINKEDLIST_MODULE
!       IMPLICIT NONE
!       CHARACTER(*),INTENT(IN) :: ID
!       INTEGER(4)  ,INTENT(IN) :: ISTP
!       INTEGER(4)  ,INTENT(IN) :: LENG
!       INTEGER(4)  ,INTENT(OUT):: VAL(LENG)
!       TYPE(LL_TYPE)           :: LL_STP_
!       INTEGER(4)  ,ALLOCATABLE:: IWORK(:)
! !     **********************************************************************
!       LL_STP_=STP(ISTP)%LL
!       CALL LINKEDLIST$SELECT(LL_STP_,'~')
!       CALL LINKEDLIST$SELECT(LL_STP_,'SETUP')
! !
! !     ==  ANGULAR MOMENTA OF THE PARTIAL WAVES  ==========================
!       IF(ID.EQ.'LOX') THEN
!         IF(LENG.NE.STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.LNX)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETI4A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'LOX',0,VAL)
! !
! !     ==  ANGULAR MOMENTUM OF THE CORE STATES =============================
!       ELSE IF(ID.EQ.'LOFC') THEN
!         IF(LENG.NE.STP(ISTP)%NC) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NC)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NC',STP(ISTP)%NC)
!           CALL ERROR$STOP('XASSTP$GETI4A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!         ALLOCATE(IWORK(STP(ISTP)%NB))
!         CALL LINKEDLIST$GET(LL_STP_,'L',0,IWORK)
!         VAL(:)=IWORK(1:STP(ISTP)%NC)
!         DEALLOCATE(IWORK)
! !
! !     == #(NODES) OF THE CORE STATES =====================================
!       ELSE IF(ID.EQ.'NNOFC') THEN
!         IF(LENG.NE.STP(ISTP)%NC) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NC)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NC',STP(ISTP)%NC)
!           CALL ERROR$STOP('XASSTP$GETI4A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!         ALLOCATE(IWORK(STP(ISTP)%NB))
!         CALL LINKEDLIST$GET(LL_STP_,'NN',0,IWORK)
!         VAL(:)=IWORK(1:STP(ISTP)%NC)
!         DEALLOCATE(IWORK)
!       ELSE
!         CALL ERROR$MSG('ID NOT RECOGNIZED')
!         CALL ERROR$CHVAL('ID',ID)
!         CALL ERROR$STOP('XASSTP$GETI4A')
!       END IF
!       RETURN
!       END
! !
! !     ...1.........2.........3.........4.........5.........6.........7.........8
!       SUBROUTINE XASSTP$GETR8A(ID,ISTP,LENG,VAL)  ! MARK: XASSTP$GETR8A
! !     **************************************************************************
! !     **                                                                      **
! !     **************************************************************************
!       USE XASSTP_MODULE, ONLY: STP
!       USE LINKEDLIST_MODULE
!       IMPLICIT NONE
!       CHARACTER(*),INTENT(IN) :: ID
!       INTEGER(4)  ,INTENT(IN) :: ISTP
!       INTEGER(4)  ,INTENT(IN) :: LENG
!       REAL(8)     ,INTENT(OUT):: VAL(LENG)
!       TYPE(LL_TYPE)           :: LL_STP_
!       REAL(8)     ,ALLOCATABLE:: WORK(:)
!       INTEGER(4)              :: I,I1,I2
!       INTEGER(4)  ,ALLOCATABLE:: LOFi(:)
! !     **************************************************************************
!       LL_STP_=STP(ISTP)%LL
!       CALL LINKEDLIST$SELECT(LL_STP_,'~')
!       CALL LINKEDLIST$SELECT(LL_STP_,'SETUP')
! !
! !     ==========================================================================
! !     == PROJECTOR FUNCTIONS                                                  ==
! !     ==========================================================================
!       IF(ID.EQ.'PRO') THEN
!         IF(LENG.NE.STP(ISTP)%NR*STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         DO I=1,STP(ISTP)%LNX
!           I1=STP(ISTP)%NR*(I-1)+1
!           I2=STP(ISTP)%NR*I
!           CALL LINKEDLIST$GET(LL_STP_,'PRO',I,VAL(I1:I2))
!         ENDDO
! !
! !     ==========================================================================
! !     == ALL-ELECTRON PARTIAL WAVES                                           ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'AEPHI') THEN
!         IF(LENG.NE.STP(ISTP)%NR*STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         DO I=1,STP(ISTP)%LNX
!           I1=STP(ISTP)%NR*(I-1)+1
!           I2=STP(ISTP)%NR*I
!           CALL LINKEDLIST$GET(LL_STP_,'AEPHI',I,VAL(I1:I2))
!         ENDDO
! !
! !     ==========================================================================
! !     == PSEUDO PARTIAL WAVES                                                 ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'PSPHI') THEN
!         IF(LENG.NE.STP(ISTP)%NR*STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         DO I=1,STP(ISTP)%LNX
!           I1=STP(ISTP)%NR*(I-1)+1
!           I2=STP(ISTP)%NR*I
!           CALL LINKEDLIST$GET(LL_STP_,'PSPHI',I,VAL(I1:I2))
!         ENDDO
! !
! !     ==========================================================================
! !     == ALL-ELECTRON CORE DENSITY                                            ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'AECORE') THEN
!         IF(LENG.NE.STP(ISTP)%NR) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'AECORE',0,VAL)
! !
! !     ==========================================================================
! !     == PSEUDO CORE DENSITY                                                  ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'PSCORE') THEN
!         IF(LENG.NE.STP(ISTP)%NR) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'PSCORE',0,VAL)
! !
! !     ==========================================================================
! !     == POTENTIAL VADD                                                       ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'VADD') THEN
!         IF(LENG.NE.STP(ISTP)%NR) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'VADD',0,VAL)
! !
! !     ==========================================================================
! !     == ONE-CENTER DIFFERENCE KINETIC ENERGY MATRIX                          ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'DT') THEN
!         IF(LENG.NE.STP(ISTP)%LNX*STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.LNX**2)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'DT',0,VAL)
! !
! !     ==========================================================================
! !     == ONE-CENTER DIFFERENCE OVERLAP MATRIX                                 ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'DO') THEN
!         IF(LENG.NE.STP(ISTP)%LNX*STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.LNX**2)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'DO',0,VAL)
! !
! !     ==========================================================================
! !     == ONE-CENTER DIFFERENCE HAMILTONIAN OF THE ATOM                        ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'DH') THEN
!         IF(LENG.NE.STP(ISTP)%LNX*STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.LNX**2)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         CALL LINKEDLIST$GET(LL_STP_,'DH',0,VAL)
! !
! !     ==========================================================================
! !     == ALL-ELECTRON ATOMIC POTENTIAL                                        ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'AEPOT') THEN
!         IF(LENG.NE.STP(ISTP)%NR) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!         CALL LINKEDLIST$GET(LL_STP_,'AEPOT',0,VAL)
! !
! !     ==========================================================================
! !     == ENERGIES OF THE CORE STATES                                          ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'EOFC') THEN
!         IF(LENG.NE.STP(ISTP)%NC) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NC)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NC',STP(ISTP)%NC)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!         ALLOCATE(WORK(STP(ISTP)%NB))
!         CALL LINKEDLIST$GET(LL_STP_,'E',0,WORK)
!         VAL(:)=WORK(1:STP(ISTP)%NC)
!         DEALLOCATE(WORK)
! !
! !     ==========================================================================
! !     == OCCUPATIONS OF THE CORE STATES                                       ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'FOFC') THEN
!         IF(LENG.NE.STP(ISTP)%NC) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NC)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NC',STP(ISTP)%NC)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!         ALLOCATE(WORK(STP(ISTP)%NB))
!         CALL LINKEDLIST$GET(LL_STP_,'OCC',0,WORK)
!         VAL(:)=WORK(1:STP(ISTP)%NC)
!         DEALLOCATE(WORK)
! !
! !     ==========================================================================
! !     == NODELESS PARTIAL WAVES                                               ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'NDLSPHI') THEN
!         IF(LENG.NE.STP(ISTP)%NR*STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         DO I=1,STP(ISTP)%LNX
!           I1=STP(ISTP)%NR*(I-1)+1
!           I2=STP(ISTP)%NR*I
!           CALL LINKEDLIST$GET(LL_STP_,'NDLSPHI',I,VAL(I1:I2))
!         ENDDO
! !
! !     ==========================================================================
! !     == KINETIC ENERGY OF NODELESS PARTIAL WAVES                             ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'NDLSTPHI') THEN
!         IF(LENG.NE.STP(ISTP)%NR*STP(ISTP)%LNX) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$I4VAL('LNX',STP(ISTP)%LNX)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         CALL LINKEDLIST$SELECT(LL_STP_,'AUGMENTATION')
!         DO I=1,STP(ISTP)%LNX
!           I1=STP(ISTP)%NR*(I-1)+1
!           I2=STP(ISTP)%NR*I
!           CALL LINKEDLIST$GET(LL_STP_,'NDLSTPHI',I,VAL(I1:I2))
!         ENDDO
! !
! !     ==========================================================================
! !     == READ SET OF NODELESS ATOMIC CORE STATES                              ==
! !     ==========================================================================
!       ELSE IF(ID.EQ.'AEPSICORE') THEN
!         CALL LINKEDLIST$SELECT(LL_STP_,'AESCF')
!         IF(LENG.NE.STP(ISTP)%NR*STP(ISTP)%NC) THEN
!           CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*NC)')
!           CALL ERROR$CHVAL('ID',ID)
!           CALL ERROR$I4VAL('LENG',LENG)
!           CALL ERROR$I4VAL('NR',STP(ISTP)%NR)
!           CALL ERROR$I4VAL('Nc',STP(ISTP)%NC)
!           CALL ERROR$STOP('XASSTP$GETR8A')
!         END IF
!         DO I=1,STP(ISTP)%NC
!           I1=1+(I-1)*STP(ISTP)%NR
!           I2=I1-1+STP(ISTP)%NR
!           CALL LINKEDLIST$GET(LL_STP_,'PHI-NODELESS',I,VAL(I1:I2))
!         ENDDO
!         ALLOCATE(LOFi(STP(ISTP)%NC))
!         CALL XASSTP$GETI4A('LOFC',ISTP,STP(ISTP)%NC,LOFi)
!         CALL XASSTP_NDLSTONODAL(STP(ISTP)%GID,STP(ISTP)%NR,STP(ISTP)%NC,LOFI,VAL)
!         DEALLOCATE(LOFi)
! !
! !     ==========================================================================
! !     == WRONG ID                                                             ==
! !     ==========================================================================
!       ELSE
!         CALL ERROR$MSG('ID NOT RECOGNIZED')
!         CALL ERROR$CHVAL('ID',ID)
!         CALL ERROR$STOP('XASSTP$GETR8A')
!       END IF
!       RETURN
!       END

! !
! !     ...1.........2.........3.........4.........5.........6.........7.........8
!       SUBROUTINE XASSTP_NDLSTONODAL(GID,NR,NC,LOFI,PHI)  ! MARK: XASSTP_NDLSTONODAL
! !     **************************************************************************
! !     ** TRANSFORM NODELESS CORES STATES INTO CORRECT ALL-ELECTRON CORE STATES**
! !     **************************************************************************
!       IMPLICIT NONE
!       INTEGER(4)  ,INTENT(IN)    :: GID
!       INTEGER(4)  ,INTENT(IN)    :: NR
!       INTEGER(4)  ,INTENT(IN)    :: NC
!       INTEGER(4)  ,INTENT(IN)    :: LOFI(NC)
!       REAL(8)     ,INTENT(INOUT) :: PHI(NR,NC)
!       REAL(8)                    :: R(NR)
!       REAL(8)                    :: AUX(NR)
!       REAL(8)                    :: VAL
!       INTEGER(4)                 :: IB1,IB2
! !     **************************************************************************
!       call radial$r(gid,nr,r)
!       DO IB1=1,NC
! !
! !       ========================================================================
! !       == NORMALIZE                                                          ==
! !       ========================================================================
!         AUX(:)=R(:)**2*PHI(:,IB1)**2
!         CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!         PHI(:,IB1)=PHI(:,IB1)/SQRT(VAL)
! !
! !       ========================================================================
! !       == ORTHOGONALIZE HIGHER STATES                                        ==
! !       ========================================================================
!         DO IB2=IB1+1,NC
!           IF(LOFI(IB2).NE.LOFI(IB1)) CYCLE
!           AUX(:)=R(:)**2*PHI(:,IB1)*PHI(:,IB2)
!           CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!           PHI(:,IB2)=PHI(:,IB2)-PHI(:,IB1)*VAL
!         ENDDO
!       ENDDO
!       RETURN
!       END
!      
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INITIALIZEFILEHANDLER
!     **************************************************************************
!     **************************************************************************
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: XASINNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: NARGS
!     **************************************************************************
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
      RETURN
      END SUBROUTINE INITIALIZEFILEHANDLER
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STANDARDFILES
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(32)        :: ID
!     **************************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
!  
!     ==========================================================================
!     == SET STANDARD FILENAMES                                               ==
!     ==========================================================================
!
!     ==  ERROR FILE ===========================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOL FILE ========================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =====================================================
      ID=+'XCNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PDOS FILE   ==========================================================
      ID=+'GROUNDSTATE'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.GROUND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     ==  PDOS FILE   ==========================================================
      ID=+'EXCITESTATE'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.EXCITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
! !
! !     ==  DENSITY OF STATES FILE PRODUCES AS OUTPUT ============================
! !     ==  WILL BE ATTACHED TO DIFFERENT FILES DURING EXECUTION =================
!       ID=+'PDOSOUT'
!       CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PDOSOUT')
!       CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
!       CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
!       CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
!       CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
! !
! !     ==  NUMBER OF STATES FILE PRODUCES AS OUTPUT =============================
! !     ==  WILL BE ATTACHED TO DIFFERENT FILES DURING EXECUTION =================
!       ID=+'PNOSOUT'
!       CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PNOSOUT')
!       CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
!       CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
!       CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
!       CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
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
      LOGICAL(4) :: TCHK
!     **************************************************************************
                          CALL TRACE$PUSH('XAS$READ')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO SIMULATION SELECTED')
        CALL ERROR$STOP('XAS$READ')
      END IF
      CALL FILEHANDLER$UNIT(THIS%FLAG,NFIL)
      REWIND(NFIL)
!     ==========================================================================
!     == READ GENERAL QUANTATIES                                              ==
!     ==========================================================================
      READ(NFIL)THIS%NAT,THIS%NSP,THIS%NKPT,THIS%NSPIN,THIS%NDIM,THIS%NPRO, &
     &          THIS%LNXX,THIS%FLAG
      ALLOCATE(THIS%LNX(THIS%NSP))
      ALLOCATE(LOX(THIS%LNXX,THIS%NSP))
      READ(NFIL)THIS%LNX,THIS%LOX
      READ(NFIL)THIS%NKDIV,THIS%ISHIFT,THIS%RNTOT,THIS%NEL,ILOGICAL
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
      READ(NFIL)THIS%RBAS,THIS%R,THIS%ATOMID,THIS%ISPECIES
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
        READ(NFIL)THIS%SETUP(ISP)%PSPHI
        READ(NFIL)THIS%SETUP(ISP)%AEPHI
      ENDDO
!     ==========================================================================
!     == READ PROJECTIONS                                                     ==
!     ==========================================================================
! TODO: CONTINUE HERE
                          CALL TRACE$POP
      END SUBROUTINE XAS$READ
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XAS_READPHI(NFIL)  ! MARK: XAS_READPHI
!     **************************************************************************
!     ** READ RADIAL GRIDS AND INITIALIZE RADIAL MODULE                       **
!     **************************************************************************
      USE XAS_MODULE, ONLY: THIS
      USE RADIAL_MODULE, ONLY: NGID
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4) :: NSP
      INTEGER(4) :: ISP
      INTEGER(4) :: IGID
      INTEGER(4) :: IGID_
      CHARACTER(8) :: GRIDTYPE
      CHARACTER(8) :: GRIDTYPE_
      INTEGER(4) :: NR
      INTEGER(4) :: NR_
      REAL(8) :: DEX
      REAL(8) :: DEX_
      REAL(8) :: R1
      REAL(8) :: R1_
      LOGICAL(4) :: TCHK
!     **************************************************************************
                          CALL TRACE$PUSH('XAS_READPHI')
!     NUMBER OF SETUPS/SPECIES
      READ(NFIL) THIS%NSP
!     CHECK IF SETUP ALREADY ALLOCATED
      IF(ALLOCATED(THIS%SETUP)) THEN
        CALL ERROR$MSG('SETUP ALREADY ALLOCATED')
        CALL ERROR$STOP('XAS_READPHI')
      END IF
      ALLOCATE(THIS%SETUP(THIS%NSP))
!     LOOP THROUGH SETUPS
      DO ISP=1,THIS%NSP
        READ(NFIL)GRIDTYPE,NR,DEX,R1,THIS%SETUP(ISP)%LNX
        THIS%SETUP(ISP)%NR=NR
!       CHECK IF RADIAL GRID WITH SAME PROPERTIES ALREADY EXISTS
        TCHK=.FALSE.
        IF(NGID.GT.0) THEN
          DO IGID=1,NGID
            CALL RADIAL$GETCH(IGID,'TYPE',GRIDTYPE_)
            CALL RADIAL$GETI4(IGID,'NR',NR_)
            CALL RADIAL$GETR8(IGID,'DEX',DEX_)
            CALL RADIAL$GETR8(IGID,'R1',R1_)
            IF(GRIDTYPE.EQ.GRIDTYPE_.AND.NR.EQ.NR_.AND.DEX.EQ.DEX_.AND.R1.EQ.R1_) THEN
              IGID_=IGID
              TCHK=.TRUE.
              EXIT
            END IF
          ENDDO
        ENDIF
!       IF GRID ALREADY EXISTS, USE ITS GID. ELSE CREATE NEW GRID
        IF(TCHK) THEN
          THIS%SETUP(ISP)%GID=IGID_
        ELSE
          CALL RADIAL$NEW(GRIDTYPE,THIS%SETUP(ISP)%GID)
          CALL RADIAL$SETI4(THIS%SETUP(ISP)%GID,'NR',NR)
          CALL RADIAL$SETR8(THIS%SETUP(ISP)%GID,'DEX',DEX)
          CALL RADIAL$SETR8(THIS%SETUP(ISP)%GID,'R1',R1)
        ENDIF

!       ALLOCATE AND READ (AUXILIARY) PARTIAL WAVES
        ALLOCATE(THIS%SETUP(ISP)%AEPHI(NR,THIS%SETUP(ISP)%LNX))
        READ(NFIL) THIS%SETUP(ISP)%AEPHI
        ALLOCATE(THIS%SETUP(ISP)%PSPHI(NR,THIS%SETUP(ISP)%LNX))
        READ(NFIL) THIS%SETUP(ISP)%PSPHI
      ENDDO
                          CALL TRACE$POP
      END SUBROUTINE XAS_READPHI

      SUBROUTINE XAS$GETCH(ID,VAL)  ! MARK: XAS$GETCH
!     **************************************************************************
!     **                                                                      **
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
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('XAS$GETCH')
      END IF
      RETURN
      END SUBROUTINE XAS$GETCH

      SUBROUTINE XAS$SETCH(ID,VAL)  ! MARK: XAS$SETCH
!     **************************************************************************
!     **                                                                      **
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
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('XAS$SETCH')
      END IF
      RETURN
      END SUBROUTINE XAS$SETCH

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
      REAL(8), ALLOCATABLE :: R(:)
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


