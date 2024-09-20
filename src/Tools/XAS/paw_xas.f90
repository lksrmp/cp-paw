!........1.........2.........3.........4.........5.........6.........7.........8
      MODULE XAS_MODULE
        TYPE STATE_TYPE
          INTEGER(4)          :: NB
          REAL(8), POINTER    :: EIG(:)
          REAL(8), POINTER    :: OCC(:)
          COMPLEX(8), POINTER :: VEC(:,:,:)  
        END TYPE STATE_TYPE
        TYPE PDOS_TYPE
          CHARACTER(6)                 :: FLAG
          INTEGER(4)                   :: NAT
          INTEGER(4)                   :: NSP
          INTEGER(4)                   :: NKPT
          INTEGER(4)                   :: NSPIN
          INTEGER(4)                   :: NDIM
          INTEGER(4)                   :: NPRO
          INTEGER(4)                   :: NKDIV(3)
          INTEGER(4)                   :: ISHIFT(3)
          REAL(8)                      :: RNTOT
          REAL(8)                      :: NEL
          LOGICAL(4)                   :: TINV
          INTEGER(4)                   :: LNXX
          REAL(8)                      :: RBAS(3,3)
          REAL(8)   ,ALLOCATABLE       :: R(:,:)
          INTEGER(4),ALLOCATABLE       :: LNX(:)
          INTEGER(4),ALLOCATABLE       :: LOX(:,:)
          INTEGER(4),ALLOCATABLE       :: ISPECIES(:)
          CHARACTER(16),ALLOCATABLE    :: ATOMID(:)
          INTEGER(4),ALLOCATABLE       :: IZ(:)
          REAL(8)   ,ALLOCATABLE       :: RAD(:)
          REAL(8)   ,ALLOCATABLE       :: PHIOFR(:,:)
          REAL(8)   ,ALLOCATABLE       :: DPHIDR(:,:)
          REAL(8)   ,ALLOCATABLE       :: OV(:,:,:)
          REAL(8)   ,ALLOCATABLE       :: XK(:,:)
          REAL(8)   ,ALLOCATABLE       :: WKPT(:)
          TYPE(STATE_TYPE),ALLOCATABLE :: STATEARR(:,:)
          TYPE(STATE_TYPE),POINTER     :: STATE
          INTEGER(4)                   :: SPACEGROUP
          LOGICAL(4)                   :: TSHIFT
        END TYPE PDOS_TYPE
        TYPE(PDOS_TYPE), TARGET :: PD(2) 
      END MODULE XAS_MODULE

      PROGRAM PAW_XAS
      use XAS_MODULE, only: PD
      IMPLICIT NONE
      INTEGER(4) :: NFIL=10
      INTEGER(4) :: IPD
      INTEGER(4) :: I
      CHARACTER(256) :: PDOS_FILES(2)
      DO I=1,2
        CALL GET_COMMAND_ARGUMENT(I, PDOS_FILES(I))
        ! open the files
        OPEN(UNIT=NFIL,FILE=PDOS_FILES(I),STATUS='OLD',FORM='UNFORMATTED')
        CALL XAS$READ(NFIL,I)
      ENDDO
      DO I=1,2
        WRITE(*,*) '*********************************************'
        WRITE(*,*) 'FILENAME = ', TRIM(ADJUSTL(PDOS_FILES(I)))
        WRITE(*,*) 'PDOS INDEX = ', I
        WRITE(*,*) 'NAT = ',        PD(I)%NAT
        WRITE(*,*) 'NSP = ',        PD(I)%NSP
        WRITE(*,*) 'NKPT = ',       PD(I)%NKPT
        WRITE(*,*) 'NSPIN = ',      PD(I)%NSPIN
        WRITE(*,*) 'NDIM = ',       PD(I)%NDIM
        WRITE(*,*) 'NPRO = ',       PD(I)%NPRO
        WRITE(*,*) 'NKDIV = ',      PD(I)%NKDIV
        WRITE(*,*) 'ISHIFT = ',     PD(I)%ISHIFT
        WRITE(*,*) 'RNTOT = ',      PD(I)%RNTOT
        WRITE(*,*) 'NEL = ',        PD(I)%NEL
        WRITE(*,*) 'TINV = ',       PD(I)%TINV
        WRITE(*,*) 'LNXX = ',       PD(I)%LNXX
        WRITE(*,*) 'RBAS = ',       PD(I)%RBAS
        WRITE(*,*) 'SPACEGROUP = ', PD(I)%SPACEGROUP
        WRITE(*,*) 'TSHIFT = ',     PD(I)%TSHIFT
      ENDDO
      
      END PROGRAM PAW_XAS
!
!     ..................................................................
      SUBROUTINE XAS$READ(NFIL,IPD)
!     ******************************************************************
!     ** READS PDOS FILE FROM NFIL INTO THE XAS MODULE                **
!     ** INPUT: NFIL - UNIT NUMBER OF PDOS FILE                       **
!     **        IPD  - INDEX OF PDOS TO BE READ                       **
!     ******************************************************************
      USE XAS_MODULE, ONLY: PD
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4),INTENT(IN)  :: IPD
      INTEGER(4)             :: ISP,IKPT,ISPIN,IB
      INTEGER(4)             :: LNX1,NB
      INTEGER(4)             :: IOS
      CHARACTER(82)          :: IOSTATMSG
      LOGICAL(4)             :: TCHK
      REAL(8)                :: OCCSUM
      INTEGER(4)             :: ILOGICAL

      INTEGER(4)             :: NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,SPACEGROUP
      LOGICAL(4)             :: TINV,TSHIFT
      INTEGER(4)                               :: NKDIV(3)
      INTEGER(4)                               :: ISHIFT(3)
      REAL(8)                                  :: RNTOT,NEL
      REAL(8)                                  :: RBAS(3,3)
      CHARACTER(6)           :: FLAG
!     ******************************************************************
                             CALL TRACE$PUSH('XAS$READ')
!     CHECK SELECTION OF PDOS INDEX
      IF(IPD.NE.1.AND.IPD.NE.2)THEN
        CALL ERROR$MSG('PDOS INDEX MUST BE 1 OR 2')
        CALL ERROR$I4VAL('IPD',IPD)
        CALL ERROR$STOP('XAS$READ')
      END IF
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      TCHK=.FALSE.
      READ(NFIL,ERR=100)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG
      TCHK=.TRUE.
 100  CONTINUE
      IF(.NOT.TCHK) THEN
        PRINT*,'WARNING: NO OCCUPATIONS PRESENT IN PDOS FILE'
        PRINT*,'            OCCUPATIONS WILL BE SET TO 0'
        FLAG='LEGACY'
        REWIND(NFIL)
        READ(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX
      END IF
      
      ALLOCATE(PD(IPD)%LNX(NSP))
      ALLOCATE(PD(IPD)%LOX(LNXX,NSP))
      ALLOCATE(PD(IPD)%ISPECIES(NAT))
      READ(NFIL)PD(IPD)%LNX(:),PD(IPD)%LOX(:,:),PD(IPD)%ISPECIES(:)
      
      IF(FLAG.EQ.'181213')THEN
!         == GFORTRAN LOGICAL REPRESENTATION DEFINED WITH TRUE=1, FALSE=0     ==
!         https://gcc.gnu.org/onlinedocs/gfortran/compiler-characteristics/
!         internal-representation-of-logical-variables.html
!         == IFORT LOGICAL REPRESENTATION DEFINED WITH VALUE OF LAST BIT      ==
!         https://www.intel.com/content/www/us/en/docs/fortran-compiler/
!         developer-guide-reference/2024-2/logical-data-representations.html
!         == BOTH SHARE MEANING OF LAST BIT 1=TRUE, 0=FALSE                   ==
!         == ENSURES BACKWARDS COMPATIBILITY WITH OLD PDOS FILES              ==
        READ(NFIL)NKDIV(:),ISHIFT(:),RNTOT,NEL,ILOGICAL
        TINV=BTEST(ILOGICAL,0)
        READ(NFIL)SPACEGROUP,ILOGICAL
        TSHIFT=BTEST(ILOGICAL,0)
      ENDIF
!
!     ==================================================================
!     == ATOMIC STRUCTURE                                             ==
!     ==================================================================
      ALLOCATE(PD(IPD)%R(3,NAT))
      ALLOCATE(PD(IPD)%ATOMID(NAT))
      READ(NFIL)RBAS(:,:),PD(IPD)%R(:,:),PD(IPD)%ATOMID(:)
!
!     ==================================================================
!     == ELEMENT SPECIFIC QUANTITIES                                  ==
!     ==================================================================
      ALLOCATE(PD(IPD)%IZ(NSP))
      ALLOCATE(PD(IPD)%RAD(NSP)); PD(IPD)%RAD=0.D0
      ALLOCATE(PD(IPD)%PHIOFR(LNXX,NSP)); PD(IPD)%PHIOFR=0.D0
      ALLOCATE(PD(IPD)%DPHIDR(LNXX,NSP)); PD(IPD)%DPHIDR=0.D0
      ALLOCATE(PD(IPD)%OV(LNXX,LNXX,NSP)); PD(IPD)%OV=0.D0
      DO ISP=1,NSP
        LNX1=PD(IPD)%LNX(ISP)
        READ(NFIL)PD(IPD)%IZ(ISP),PD(IPD)%RAD(ISP),PD(IPD)%PHIOFR(1:LNX1,ISP) &
     &            ,PD(IPD)%DPHIDR(1:LNX1,ISP),PD(IPD)%OV(1:LNX1,1:LNX1,ISP)
      ENDDO
!
!     ==================================================================
!     ==  NOW READ PROJECTIONS                                       ==
!     ==================================================================
      OCCSUM=0.0D0
      ALLOCATE(PD(IPD)%XK(3,NKPT))
      ALLOCATE(PD(IPD)%WKPT(NKPT))
      ALLOCATE(PD(IPD)%STATEARR(NKPT,NSPIN))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          PD(IPD)%STATE=>PD(IPD)%STATEARR(IKPT,ISPIN)
          READ(NFIL,ERR=9998,END=9998)PD(IPD)%XK(:,IKPT),NB,PD(IPD)%WKPT(IKPT)
          PD(IPD)%STATE%NB=NB
          ALLOCATE(PD(IPD)%STATE%EIG(NB))
          ALLOCATE(PD(IPD)%STATE%VEC(NDIM,NPRO,NB))
          ALLOCATE(PD(IPD)%STATE%OCC(NB))
          DO IB=1,NB
            IF(FLAG.EQ.'LEGACY') THEN
              PD(IPD)%STATE%OCC(:)=0.D0
              READ(NFIL,ERR=9999,IOSTAT=IOS)PD(IPD)%STATE%EIG(IB),PD(IPD)%STATE%VEC(:,:,IB)
            ELSE
              READ(NFIL,ERR=9999,IOSTAT=IOS)PD(IPD)%STATE%EIG(IB) &
    &                          ,PD(IPD)%STATE%OCC(IB),PD(IPD)%STATE%VEC(:,:,IB)
              OCCSUM=OCCSUM+PD(IPD)%STATE%OCC(IB)
            END IF
          ENDDO
        ENDDO
      ENDDO
PRINT*,"OCCSUM",OCCSUM
!     SET FIXED SIZE VARIABLES IN PDOS STRUCTURE
      PD(IPD)%FLAG=FLAG
      PD(IPD)%NAT=NAT
      PD(IPD)%NSP=NSP
      PD(IPD)%NKPT=NKPT
      PD(IPD)%NSPIN=NSPIN
      PD(IPD)%NDIM=NDIM
      PD(IPD)%NPRO=NPRO
      PD(IPD)%NKDIV(:)=NKDIV(:)
      PD(IPD)%ISHIFT(:)=ISHIFT(:)
      PD(IPD)%RNTOT=RNTOT
      PD(IPD)%NEL=NEL
      PD(IPD)%TINV=TINV
      PD(IPD)%LNXX=LNXX
      PD(IPD)%RBAS(:,:)=RBAS(:,:)
      PD(IPD)%SPACEGROUP=SPACEGROUP
      PD(IPD)%TSHIFT=TSHIFT
                             CALL TRACE$POP
      RETURN
 9999 CONTINUE
      CALL FILEHANDLER$IOSTATMESSAGE(IOS,IOSTATMSG)
      CALL ERROR$MSG('ERROR READING PDOS FILE')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL ERROR$CHVAL('IOSTATMSG',IOSTATMSG)
      CALL ERROR$I4VAL('IB',IB)
      CALL ERROR$I4VAL('IKPT',IKPT)
      CALL ERROR$I4VAL('ISPIN',ISPIN)
      CALL ERROR$I4VAL('NPRO',NPRO)
      CALL ERROR$STOP('XAS$READ')
      STOP
 9998 CONTINUE
      CALL ERROR$MSG('ERROR READING PDOS FILE')
      CALL ERROR$MSG('OLD VERSION: VARIABLE WKPT IS NOT PRESENT')
      CALL ERROR$MSG('PRODUCE NEW PDOS FILE')
      CALL ERROR$STOP('XAS$READ')
      STOP
      END SUBROUTINE XAS$READ