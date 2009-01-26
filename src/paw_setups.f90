!===============================================================================
!===============================================================================
!===============================================================================
!====      FORMFACTORS                                                      ====
!===============================================================================
!===============================================================================
!===============================================================================
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SETUP_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: SETUP                                                              **
!**                                                                           **
!**  PURPOSE:                                                                 **
!**                                                                           **
!**  READS PSEUDOPOTENTIALS AND ATOMIC PSEUDO-WAVEFUNCTIONS                   **
!**  GIVEN ON A LINEAR, RADIAL MESH                                           **
!**              AND                                                          **
!**  CALCULATES THE FORMFACTORS FOR LOCAL AND NONLOCAL                        **
!**  CONTRIBUTIONS TO THE PSEUDOPOTENTIAL                                     **
!**  IN PLANE WAVE REPRESENTATION                                             **
!**                                                                           **
!**  THE LOCAL COMPONENT CONTAINS ALSO THE POTENTIAL OF A                     **
!**  (GAUSSIAN SHAPED) CHARGEDENSITY, WHICH COMPENSATES                       **
!**  THE IONIC CHARGE.                                                        **
!**                                                                           **
!**                                                                           **
!**  RCBG (R) = -ZATOM/(SQRT(PI)*RCRHO)**3 * EXP(-(R/RCRHO)**2)               **
!**                                                                           **
!**  OUTPUT:                                                                  **
!**    RCRHO     RADIUS OF GAUSSIAN IN RHOPS                                  **
!**    RHOPS     COMPENSATION CHARGE DENSITY                                  **
!**    VLOC      LOCAL CONTRIBUTION TO PSEUDOPOTENTIAL                        **
!**    WNL       NONLOCAL CONTRIBUTION (PROJECTOR)                            **
!**                                                                           **
!**                                                                           **
!**                                                                           **
!**  LOCAL ORBITALS:                                                          **
!**    LOCAL ORBITALS ARE APPROXIMATIONS FOR WANNIER FUNCTIONS FOR TREATING   **
!**    LOCAL CORRELATION EFFECTS AND FOR CHEMICAL BOND ANALYSIS               **
!**    THEY ARE PARTIAL WAVES SUPERIMPOSED SO THAT THEY HAVE A NODE           **
!**    AT A GIVEN RADIUS.                                                     **
!**                                                                           **
!**      |CHI_I>=\SUM_J |PHI_J> AMAT(J,I)                                     **
!**      |PSI>=\SUM_J |CHI_J> BMAT(J,I) <PRO|PSITILDE>                        **
!**                                                                           **
!**    SETR4('RADCHI',RAD)                                                    **
!**    SETI4A('NOFLCHI',4,NOFL)                                               **
!**    GETI4('LNXCHI',LNXCHI)                                                 **
!**    GETI4A('LOXCHI',LNXCHI,LOXCHI)                                         **
!**    GETR8A('AMATCHI',LNXPHI*LNXCHI,AMAT)                                   **
!**    GETR8A('BMATCHI',LNXCHI*LNXPHI,BMAT)                                   **
!**                                                                           **
!**                                                                           **
!**  AN UNSELECT FUNCTION HAS BEEN INCLUDED BUT NOT FULLY IMPLEMENTED.        **
!**  BY FORCING UNSELECT BEFORE SELECT, ONE CAN SAFEGUARD THAT A SUBROUTINE   **
!**  CHANGES THE SETTING OF A PARENT ROUTINE.                                 **
!**                                                                           **
!**                                                                           **
!**                                              P.E. BLOECHL, (1991-2008)    **
!*******************************************************************************
TYPE SETUPPARMS_TYPE
  CHARACTER(128)  :: ID 
  REAL(8)         :: POW_POT=0.D0
  logical(4)      :: tval0_pot       
  REAL(8)         :: VAL0_POT
  REAL(8)         :: RC_POT
  REAL(8)         :: POW_CORE
  logical(4)      :: tval0_core       
  REAL(8)         :: VAL0_CORE
  REAL(8)         :: RC_CORE
  CHARACTER(32)   :: TYPE
  REAL(8),POINTER :: RCL(:)
  REAL(8),POINTER :: LAMBDA(:)
END TYPE SETUPPARMS_TYPE
TYPE ATOMWAVES_TYPE
  INTEGER(4)         :: NB=-1
  INTEGER(4)         :: NC=-1
  INTEGER(4),POINTER :: LOFI(:)
  REAL(8)   ,POINTER :: EOFI(:)
  REAL(8)   ,POINTER :: FOFI(:)
  REAL(8)   ,POINTER :: AEPSI(:,:)
  REAL(8)   ,POINTER :: AEPOT(:)
END TYPE ATOMWAVES_TYPE
TYPE THIS_TYPE
INTEGER(4)             :: I            ! INTEGER IDENTIFIER (ISPECIES)
CHARACTER(32)          :: ID           ! IDENTIFIER (SPECIES-NAME)
INTEGER(4)             :: GID          ! GRID ID FOR R-SPACE GRID
INTEGER(4)             :: GIDG         ! GRID ID FOR G-SPACE GRID
REAL(8)                :: AEZ          ! ATOMIC NUMBER
REAL(8)                :: RCBG
REAL(8)                :: RCSM         ! GAUSSIAN DECAY FOR COMPENSATION CHARGE
INTEGER(4)             :: LX           ! HIGHEST ANGULAR MOMENTUM
INTEGER(4)             :: LNX          ! #(ORBITAL SHELLS)
INTEGER(4)             :: LMNX
INTEGER(4)             :: LMRX         ! #(ANGULAR MOMENTA FOR 1C-DENSITY)
INTEGER(4),POINTER     :: LOX(:)       !(LNX) MAIN ANGULAR MOMENTA 
INTEGER(4),POINTER     :: ISCATT(:)    !(LNX) =-1 FOR SEMI-CORE STATE
                                       !      = 0 FOR VALENCE STATE   (PHI)
                                       !      = 1 FOR 1ST SCATTERING STATE (PHIDOT)
REAL(8)   ,POINTER     :: VADD(:)      !(NR)
REAL(8)   ,POINTER     :: PSPOT(:)     !(NR)
REAL(8)   ,POINTER     :: AECORE(:)    !(NR)  CORE ELECTRON DENSITY
REAL(8)   ,POINTER     :: PSCORE(:)    !(NR)  PSEUDIZED ELECTRON DENSITY
REAL(8)   ,POINTER     :: PRO(:,:)     !(NR,LNX)  PROJECTOR FUNCTIONS
REAL(8)                :: RBOX         ! PARTIAL WAVES HAVE OUTER NODE AT RBOX
REAL(8)   ,POINTER     :: EOFLN(:)     !(LNX)  ENERGIES FOR PARTIAL WAVE CONSTRUCTION
REAL(8)   ,POINTER     :: ESCATT(:)    !(LNX)  ENERGIES FOR SCATTERING WAVE FUNCTIONS
REAL(8)   ,POINTER     :: AEPHI(:,:)   !(NR,LNX)  AE PARTIAL WAVES
REAL(8)   ,POINTER     :: PSPHI(:,:)   !(NR,LNX)  PS PARTIAL WAVES
REAL(8)   ,POINTER     :: UPHI(:,:)    !(NR,LNX)  NODELESS PARTIAL WAVES
REAL(8)   ,POINTER     :: TUPHI(:,:)   !(NR,LNX)  KINETIC ENERGY OF ABOVE
REAL(8)   ,POINTER     :: NLPHIDOT(:,:)!(NR,LNX)  NL SCATTERING PARTIAL WAVES
REAL(8)   ,POINTER     :: PSPHIDOT(:,:)!(NR,LNX)  PS SCATTERING PARTIAL WAVES
REAL(8)   ,POINTER     :: AEPHIDOT(:,:)!(NR,LNX)  AE SCATTERING PARTIAL WAVES
REAL(8)   ,POINTER     :: DTKIN(:,:)   !(LNX,LNX) 1C-KIN. EN. MATRIX ELEMENTS
REAL(8)   ,POINTER     :: DOVER(:,:)   !(LNX,LNX) 1C-OVERLAP MATRIX ELEMENTS
REAL(8)   ,POINTER     :: VADDOFG(:)   !(NGX)
REAL(8)   ,POINTER     :: PSCOREOFG(:) !(NGX)
REAL(8)   ,POINTER     :: VHATOFG(:)   !(NGX)
REAL(8)   ,POINTER     :: NHATPRIMEOFG(:)  !(NGX)
REAL(8)   ,POINTER     :: PROOFG(:,:)  !(NR,LNX)  
LOGICAL(4)             :: LOCORBINI=.FALSE. !LOCAL ORBS ARE INITIALIZED IF TRUE
REAL(8)                :: LOCORBRAD(4)=5.D0  ! RADIUS OF LOCAL ORBITALS
INTEGER(4)             :: LOCORBNOFL(4)=0 ! #(LOCAL ORBITALS PER L)
INTEGER(4)             :: LOCORBLNX=0     ! #(LOCAL ORBITAL-SHELLS)
INTEGER(4),POINTER     :: LOCORBLOX(:)    ! L FOR EACH LOCAL ORBITAL-SHELL
REAL(8)   ,POINTER     :: LOCORBAMAT(:,:) ! |CHI>=|PHI>*AMAT
REAL(8)   ,POINTER     :: LOCORBBMAT(:,:) ! |PSI>=|CHI>BMAT<PTILDE|PSITILDE>
REAL(8)                :: M
REAL(8)                :: ZV
REAL(8)                :: PSG2
REAL(8)                :: PSG4
CHARACTER(32)          :: SOFTCORETYPE
CHARACTER(16)          :: FILEID
TYPE(SETUPPARMS_TYPE)  :: PARMS
TYPE(ATOMWAVES_TYPE)   :: ATOM
TYPE(THIS_TYPE),POINTER:: NEXT
END TYPE THIS_TYPE

!
!INTEGER(4)             :: NR    !=250
!INTEGER(4)             :: NRX=250
!REAL(8)   ,PARAMETER   :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
!INTEGER(4),PARAMETER   :: NG=256
!REAL(8)                :: G1          ! FIRST POINT ON THE RADIAL G-GRID
INTEGER(4)              :: NSP=0
LOGICAL,SAVE            :: SELECTED=.FALSE. ! CONTAINER FOR ACTUAL SETTING
INTEGER(4)              :: LMRXX=0
INTEGER(4)              :: LMNXX=0
INTEGER(4)              :: LNXX=0
TYPE(THIS_TYPE),POINTER :: FIRST
TYPE(THIS_TYPE),POINTER :: THIS
TYPE FASTACCESS_TYPE  
  TYPE(THIS_TYPE),POINTER :: THIS
END TYPE FASTACCESS_TYPE
TYPE(FASTACCESS_TYPE),ALLOCATABLE :: FASTACCESS(:)
LOGICAL(4)             :: TINTERNALSETUPS=.TRUE.
END MODULE SETUP_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$ISELECT(I)
!     **************************************************************************
!     **  SELECTS A SETUP PER INTEGER INDEX                                   **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: I
      INTEGER(4)            :: J
!     **************************************************************************
!
!     ==========================================================================
!     ==  DEFINE FASTACCESS LOOKUP ARRAY                                      ==
!     ==========================================================================
      IF(.NOT.ALLOCATED(FASTACCESS)) THEN
        IF(.NOT.ASSOCIATED(FIRST)) THEN
          CALL ERROR$MSG('NO SETUPS DEFINED')
          CALL ERROR$STOP('SETUP$ISELECT')
        END IF
        THIS=>FIRST
        NSP=1
        DO WHILE(ASSOCIATED(THIS%NEXT))
          NSP=NSP+1
          THIS=>THIS%NEXT
        ENDDO
        ALLOCATE(FASTACCESS(NSP))
        FASTACCESS(1)%THIS=>FIRST
        DO J=2,NSP
          FASTACCESS(J)%THIS=>FASTACCESS(J-1)%THIS%NEXT 
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CHECK IF I LIES OUT OF RANGE                                        ==
!     ==========================================================================
      IF(I.GT.NSP.OR.I.LT.0) THEN
        CALL ERROR$MSG('INDEX I OUT OF RANGE')
        CALL ERROR$I4VAL('I',I)
        CALL ERROR$I4VAL('NSP',NSP)
        CALL ERROR$STOP('SETUP$ISELECT')
      END IF
!
!     ==========================================================================
!     ==  SELECT OR UNSELECT                                                  ==
!     ==========================================================================
      IF(I.EQ.0) THEN
!       ==  UNSELECT ===========================================================
        IF(.NOT.SELECTED) THEN
          CALL ERROR$MSG('SAFEGUARD FUNCTION:')
          CALL ERROR$MSG('CANNOT UNSELECT A SETUP THAT IS NOT SELECTED')
          CALL ERROR$STOP('SETUP$ISELECT')
        END IF
        SELECTED=.FALSE.
        NULLIFY(THIS)
      ELSE 
!       ==  SELECT =============================================================
!        IF(SELECTED) THEN
!          CALL ERROR$MSG('SAFEGUARD FUNCTION:')
!          CALL ERROR$MSG('ANOTHER SETUP IS ALREADY SELECTED:')
!          CALL ERROR$MSG('UNSELECT BEFORE SELECT')
!          CALL ERROR$STOP('SETUP$ISELECT')
!        END IF
        THIS=>FASTACCESS(I)%THIS
        SELECTED=.TRUE.
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$UNSELECT()
!     **************************************************************************
!     **  UNSELECTS                                                           **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT UNSELECT A SETUP THAT IS NOT SELECTED')
        CALL ERROR$STOP('SETUP$UNSELECT')
      END IF
      SELECTED=.FALSE.
      NULLIFY(THIS)
      RETURN
      END
!!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$SELECT(ID)
!     **************************************************************************
!     **  SELECTS A SETUP PER ID                                              **
!     **  AND CREATES A NEW, IF IT DOES NOT EXIST                             **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
!     **************************************************************************
!
!     == CHECK IF ALREADY SELECTED
      IF(ASSOCIATED(THIS)) THEN
        IF(THIS%ID.EQ.ID) RETURN
      END IF
!
!     == CHECK IF PRESENT ===
      IF(ASSOCIATED(FIRST)) THEN
        THIS=>FIRST
        DO 
          IF(THIS%ID.EQ.ID) RETURN
          IF(.NOT.ASSOCIATED(THIS%NEXT))EXIT
          THIS=>THIS%NEXT
        ENDDO
        ALLOCATE(THIS%NEXT)
        THIS%NEXT%I=THIS%I+1
        THIS=>THIS%NEXT
      ELSE
        ALLOCATE(FIRST)
        THIS=>FIRST
        THIS%I=1
      END IF
!
!     == CREATE NEW
      IF(ALLOCATED(FASTACCESS)) DEALLOCATE(FASTACCESS)
      THIS%ID    =ID
      THIS%AEZ   =0.D0
      THIS%RCBG  =0.D0
      THIS%RCSM  =0.D0
      THIS%LNX   =0
      THIS%LMNX  =0
      THIS%LMRX  =0
      THIS%PSG2  =0.D0
      THIS%PSG4  =0.D0
      THIS%SOFTCORETYPE='NONE'
      NULLIFY(THIS%LOX)     !(LNX)
      NULLIFY(THIS%VADD)    !(NRX)
      NULLIFY(THIS%AECORE)  !(NRX)
      NULLIFY(THIS%PSCORE)  !(NRX)
      NULLIFY(THIS%PRO)     !(NRX,LNX)
      NULLIFY(THIS%AEPHI)   !(NRX,LNX)
      NULLIFY(THIS%PSPHI)   !(NRX,LNX)
      NULLIFY(THIS%DTKIN)   !(LNXX,LNX)
      NULLIFY(THIS%DOVER)   !(LNXX,LNX)
      NULLIFY(THIS%VADDOFG) !(NGX)
      NULLIFY(THIS%PSCOREOFG) !(NGX)
      NULLIFY(THIS%VHATOFG) !(NGX)
      NULLIFY(THIS%NHATPRIMEOFG) !(NGX)
      NULLIFY(THIS%PROOFG)  !(NGX,LNX)
      THIS%LOCORBINI=.FALSE. ! INITIALIZED?
      THIS%LOCORBRAD(:)=5.D0  ! RADIUS OF LOCAL ORBITALS
      THIS%LOCORBNOFL(:)=0 ! (4) #(LOCAL ORBITALS PER L)
      THIS%LOCORBLNX=0     ! #(LOCAL ORBITAL-SHELLS)
      NULLIFY(THIS%LOCORBLOX)  !(LOCORBLNX)  L FOR EACH LOCAL ORBITAL-SHELL
      NULLIFY(THIS%LOCORBAMAT) !(LNX,LOCORBLNX) |CHI>=|PHI>*AMAT
      NULLIFY(THIS%LOCORBBMAT) !(LOCORBLNX,LNX) |PSI>=|CHI>BMAT<PTILDE|PSITILDE>

      NULLIFY(THIS%NEXT)
      WRITE(THIS%FILEID,*)THIS%I
      THIS%FILEID='ATOM'//ADJUSTL(THIS%FILEID)
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETCH(ID,VAL)
!     **************************************************************************
!     **  COLLECTS INTERNAL DATA                                              **
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      CHARACTER(*),INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ID') THEN
        VAL=THIS%ID
      ELSE IF(ID.EQ.'SOFTCORETYPE') THEN
        VAL=THIS%SOFTCORETYPE
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('SETUP$GETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETL4(ID,VAL)
!     **************************************************************************
!     **  COLLECTS INTERNAL DATA                                              **
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'INTERNALSETUPS') THEN
        VAL=TINTERNALSETUPS
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('SETUP$GETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETI4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'NSP') THEN
        VAL=NSP
      ELSE IF(ID.EQ.'ISP') THEN
        VAL=THIS%I
      ELSE IF(ID.EQ.'LNXX') THEN
        VAL=LNXX
      ELSE IF(ID.EQ.'LMNXX') THEN
        VAL=LMNXX
      ELSE IF(ID.EQ.'LMRXX') THEN
        VAL=LMRXX
      ELSE IF(ID.EQ.'GID') THEN
        VAL=THIS%GID
      ELSE IF(ID.EQ.'GIDG') THEN
        VAL=THIS%GIDG
      ELSE IF(ID.EQ.'NR') THEN
        CALL RADIAL$GETI4(THIS%GID,'NR',VAL)
      ELSE IF(ID.EQ.'LNX') THEN
        VAL=THIS%LNX
      ELSE IF(ID.EQ.'LMNX') THEN
        VAL=THIS%LMNX
      ELSE IF(ID.EQ.'LMRX') THEN
        VAL=THIS%LMRX
      ELSE IF(ID.EQ.'LNXCHI') THEN
        CALL ERROR$MSG('ID LNXCHI IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$GETI4')
        CALL SETUP_RENEWLOCORB()
        VAL=THIS%LOCORBLNX    ! #(LOCAL ORBITALS)
      ELSE IF(ID.EQ.'NC') THEN
        VAL=THIS%ATOM%NC
      ELSE IF(ID.EQ.'NB') THEN
        VAL=THIS%ATOM%NB
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$SETI4A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(IN)  :: VAL(LEN)
      INTEGER(4)               :: I
!     **************************************************************************
      IF(ID.EQ.'NOFLCHI') THEN
        IF(LEN.NE.4) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$SETI4A')
        END IF
!       == CHECK IF VALUE HAS CHANGED ==========================================        
        DO I=1,4
          THIS%LOCORBINI=THIS%LOCORBINI.AND.(THIS%LOCORBNOFL(I).EQ.VAL(I))
        ENDDO
        THIS%LOCORBNOFL=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETI4A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(OUT) :: VAL(LEN)
!     **************************************************************************
!     ==========================================================================
!     == ANGULAR MOMENTA OF THE PARTIAL WAVES                                 ==
!     ==========================================================================
      IF(ID.EQ.'LOX') THEN
        IF(LEN.NE.THIS%LNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%LOX
!
!     ==========================================================================
!     == COUNTS PARTIAL WAVES RELATIVE TO THE HIGHEST VALENCE STATE           ==
!     ==========================================================================
      ELSE IF(ID.EQ.'ISCATT') THEN
        IF((LEN.NE.THIS%LNX).OR.(THIS%LNX.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        IF(.NOT.ASSOCIATED(THIS%ISCATT)) THEN
          CALL ERROR$MSG('DATA NOT AVAILABLE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%ISCATT
!
!     ==========================================================================
!     == MAIN ANGULAR MOMENTA OF ATOMIC STATES                                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'LB') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%ATOM%LOFI
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$SETR8(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(IN)  :: VAL
!     **************************************************************************
      IF(ID.EQ.'RADCHI') THEN
        CALL ERROR$MSG('ID RADCHI IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$SETR8')
!!$!       == EXTENT OF LOCAL ORBITALS
!!$        THIS%LOCORBINI=.FALSE.
!!$        THIS%LOCORBRAD(:)=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$SETR8A(ID,LENG,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LENG
      REAL(8)     ,INTENT(IN)  :: VAL(LENG)
!     **************************************************************************
      IF(ID.EQ.'RADCHI') THEN
        CALL ERROR$MSG('ID RADCHI IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$SETR8A')
!!$!       == EXTENT OF LOCAL ORBITALS
!!$        THIS%LOCORBINI=.FALSE.
!!$        IF(LENG.GT.4) THEN
!!$          CALL ERROR$MSG('INCONSISTENT SIZE')
!!$          CALL ERROR$MSG('DIMENSION OF ARRAY MUST NOT BE LARGER THAN 4')
!!$          CALL ERROR$I4VAL('LENG',LENG)
!!$          CALL ERROR$CHVAL('ID',ID)
!!$          CALL ERROR$STOP('SETUP$SETR8')
!!$        END IF
!!$        THIS%LOCORBRAD(:)=0.D0
!!$        THIS%LOCORBRAD(:)=VAL(1:LENG)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$GETR8(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VAL
!     ******************************************************************
      IF(ID.EQ.'AEZ') THEN
        VAL=THIS%AEZ
      ELSE IF(ID.EQ.'RCSM') THEN
        VAL=THIS%RCSM
      ELSE IF(ID.EQ.'RCBG') THEN
        VAL=THIS%RCBG
      ELSE IF(ID.EQ.'ZV') THEN
        VAL=THIS%ZV
      ELSE IF(ID.EQ.'M') THEN
        VAL=THIS%M
      ELSE IF(ID.EQ.'<G2>') THEN
        VAL=THIS%PSG2
      ELSE IF(ID.EQ.'<G4>') THEN
        VAL=THIS%PSG4
      ELSE IF(ID.EQ.'RBOX') THEN
        VAL=THIS%RBOX  ! OUTER NODE OF PARTIAL WAVES
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETR8A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      REAL(8)     ,INTENT(OUT) :: VAL(LEN)
      INTEGER(4)               :: NR
      INTEGER(4)               :: NC
      INTEGER(4)               :: LRHOX
!     **************************************************************************
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
!
!     ==========================================================================
!     == PROJECTOR FUNCTIONS                                                  ==
!     ==========================================================================
      IF(ID.EQ.'PRO') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PRO,(/LEN/))
! 
!     ==========================================================================
!     ==  ENERGIES FOR PARTIAL WAVE CONSTRUCTION                              ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EOFLN') THEN
        IF(LEN.NE.THIS%LNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%EOFLN
! 
!     ==========================================================================
!     ==  ENERGIES FOR SCATTERING PARTIAL WAVE CONSTRUCTION                   ==
!     ==========================================================================
      ELSE IF(ID.EQ.'ESCATT') THEN
        IF(LEN.NE.THIS%LNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ESCATT
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%AEPHI,(/LEN/))
!
!     ==========================================================================
!     ==  PSEUDO PARTIAL WAVES                                                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PSPHI,(/LEN/))
!
!     ==========================================================================
!     ==  NODELESS PARTIAL WAVES                                              ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NLPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%UPHI,(/LEN/))
!
!     ==========================================================================
!     ==  KINETIC ENERGY OPERATOR APPLIED TO NODELESS PARTIAL WAVES           ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NLTPHI') THEN
!       == DO NOT FORGET TO CLEAN UP THIS%TUPHI!!!
        CALL ERROR$MSG('OPTION IS MARKED FOR DELETION')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8A')
!
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%TUPHI,(/LEN/))
!
!     ==========================================================================
!     ==  NODELESS SCATTERING PARTIAL WAVES                                   ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NLPHIDOT') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%NLPHIDOT,(/LEN/))
!
!     ==========================================================================
!     ==  SCATTERING PSEUDO PARTIAL WAVES                                     ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSPHIDOT') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PSPHIDOT,(/LEN/))
!
!     ==========================================================================
!     ==  SCATTERING ALL-ELECTRON PARTIAL WAVES                               ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPHIDOT') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%AEPHIDOT,(/LEN/))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AECORE') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%AECORE
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSCORE') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%PSCORE
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'VADD') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%VADD
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'DEKIN') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%DTKIN,(/LEN/))
!
!     ==========================================================================
!     ==  OVERLAP DIFFERENCE MATRIX ELEMENTS <AEPSI|AEPSI>-<PSPSI|PSPSI>      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'DO') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%DOVER,(/LEN/))
!
!     ==========================================================================
!     ==  ATOMIC WAVE FUNCTIONS                                               ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPSI') THEN
        IF(LEN.NE.NR*THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%ATOM%AEPSI,(/LEN/))
!
!     ==========================================================================
!     ==  ATOMIC OCCUPATIONS                                                  ==
!     ==========================================================================
      ELSE IF(ID.EQ.'FOFI') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ATOM%FOFI
!
!     ==========================================================================
!     ==  ATOMIC ONE-PARTICLE ENERGIES                                        ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EOFI') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ATOM%EOFI
!
!     ==========================================================================
!     ==  OCCUPATION OF ATOMIC ORBITALS                                       ==
!     ==========================================================================
      ELSE IF(ID.EQ.'FB') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ATOM%FOFI
!
!     ==========================================================================
!     ==  ATOMIC ENERGY LEVELS                                                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EB') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ATOM%EOFI
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPOT') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL(:)=THIS%ATOM%AEPOT(:)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NUCPOT') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        CALL RADIAL$NUCPOT(THIS%GID,NR,THIS%AEZ,VAL)
!
!     ==========================================================================
!     ==  WRONG ID                                                            ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8A')
      END IF
      RETURN
      END  
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETFOFG(ID,TDER,IND,NG_,G2,CELLVOL,F)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID      !IDENTIFIER
      LOGICAL(4)  ,INTENT(IN)  :: TDER    ! CALCULATE RADIAL DERVATIVE
      INTEGER(4)  ,INTENT(IN)  :: IND     ! SELECTOR (USED ONLY FOR ID=PRO)
      INTEGER(4)  ,INTENT(IN)  :: NG_     ! #(PLANE WAVES)
      REAL(8)     ,INTENT(IN)  :: G2(NG_) ! G**2
      REAL(8)     ,INTENT(OUT) :: F(NG_)  
      REAL(8)     ,ALLOCATABLE :: FOFG(:)  !(NG)
      REAL(8)                  :: PI
      INTEGER(4)               :: IG
      INTEGER(4)               :: NG
      INTEGER(4)               :: NGAMMA
      REAL(8)                  :: CELLVOL
      REAL(8)                  :: G
      INTEGER(4)               :: GIDG
!     ******************************************************************
      IF(NG_.EQ.0) RETURN
      CALL RADIAL$GETI4(THIS%GIDG,'NR',NG)
      ALLOCATE(FOFG(NG))
      IF(ID.EQ.'PRO') THEN
        IF(IND.LT.1.OR.IND.GT.THIS%LNX) THEN
          CALL ERROR$MSG('LN OUT OF RANGE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('IND',IND)
          CALL ERROR$I4VAL('LNX',THIS%LNX)
          CALL ERROR$STOP('SETUP$GETFOFG')
        END IF
        FOFG(:)=THIS%PROOFG(:,IND)
      ELSE IF(ID.EQ.'PSCORE') THEN
        FOFG(:)=THIS%PSCOREOFG(:)
      ELSE IF(ID.EQ.'VADD') THEN
        FOFG(:)=THIS%VADDOFG(:)
      ELSE IF(ID.EQ.'V0') THEN
        FOFG(:)=THIS%VHATOFG(:)
      ELSE IF(ID.EQ.'G0') THEN
        FOFG(:)=THIS%NHATPRIMEOFG(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETFOFG')
      END IF
!
!     ==================================================================
!     == INTERPOLATE VALUES FROM RADIAL GRID
!     ==================================================================
      GIDG=THIS%GIDG
      NGAMMA=0
      IF(TDER) THEN
        G=SQRT(G2(1))
        IF(G.LT.1.D-6) NGAMMA=1
        CALL RADIAL$DERIVATIVE(GIDG,NG,FOFG,G,F(1))
        F(1)=G*F(1)
        DO IG=2,NG_
          IF(ABS(G2(IG)-G2(IG-1)).LT.1.D-6) THEN
            F(IG) =F(IG-1)
          ELSE
            G=SQRT(G2(IG))
            IF(G.LT.1.D-6) NGAMMA=IG
            CALL RADIAL$DERIVATIVE(GIDG,NG,FOFG,G,F(IG))
            F(IG)=G*F(IG)
          END IF
        ENDDO
      ELSE
        G=SQRT(G2(1))
        IF(G.LT.1.D-6) NGAMMA=1
        CALL RADIAL$VALUE(GIDG,NG,FOFG,G,F(1))
        DO IG=2,NG_
          IF(ABS(G2(IG)-G2(IG-1)).LT.1.D-6) THEN
            F(IG) =F(IG-1)
          ELSE
            G=SQRT(G2(IG))
            IF(G.LT.1.D-6) NGAMMA=IG
            CALL RADIAL$VALUE(GIDG,NG,FOFG,G,F(IG))
          END IF
        ENDDO
      END IF
      DEALLOCATE(FOFG)
!
!     ==================================================================
!     == CORRECT EXTRAPOLATION TO THE GAMMA POINT                     ==
!     ==================================================================
      IF(NGAMMA.NE.0) THEN
        PI=4.D0*ATAN(1.D0)
        IF(TDER) THEN
          NGAMMA=0.D0
        ELSE
          IF(ID.EQ.'G0') THEN 
            F(NGAMMA)=4.D0*PI
          ELSE IF(ID.EQ.'V0') THEN
            F(NGAMMA)=PI*(THIS%RCBG**2-THIS%RCSM**2)*4.D0*PI
          END IF
        END IF
      END IF
!
!     ==================================================================
!     == DIVIDE BY CELLVOL                                            ==
!     ==================================================================
      F=F/CELLVOL
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_READ
!     ******************************************************************
!     **  READ SELECTED SETUP                                         **
!     **  REQUIRES INFORMATION FROM ATOMTYPELIST                      **
!     **    NAME; LRHOX                                               **
!     **  REQUIRES THE FILEHANDLER TO KNOW THE SETUP FILE             **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      REAL(8)   ,PARAMETER  :: TOL=1.D-6
      INTEGER(4)            :: LMRXCUT
      INTEGER(4)            :: GID
      INTEGER(4)            :: GIDG
      INTEGER(4)            :: NFIL
      REAL(8)               :: G1
      REAL(8)   ,PARAMETER  :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
      INTEGER(4),PARAMETER  :: NG=250
      REAL(8)               :: R1,DEX
      INTEGER(4)            :: NR,NRX
      INTEGER(4)            :: IR
      INTEGER(4)            :: LN
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: IRMAX
      INTEGER(4)            :: L,LX,ISVAR,LNOLD,LNX
      INTEGER(4)            :: NC !SANTOS040616
      INTEGER(4)            :: LN1,LN2,LN1A,LN2A
      INTEGER(4),ALLOCATABLE:: NPRO(:)
      INTEGER(4),ALLOCATABLE:: IWORK(:)
      REAL(8)   ,ALLOCATABLE:: DWORK(:,:,:) 
      REAL(8)               :: PI,FOURPI,Y0
      LOGICAL(4)            :: TNEWFORMAT
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)               :: SVAR
      REAL(8)               :: PSZ   ! LEGACY ONLY
!     ******************************************************************
                            CALL TRACE$PUSH('SETUP_READ')
!     == CERTAIN OPTIONS ONLY WORK IF ALL SETUPS ARE CALCULATED  =======
!     == INTERNALLY USING SETUP_READ_NEW. THIS SWITCH ALLOWS TO ========
!     == LOCK OPTIONS THAT DO NOT WORK =================================
      TINTERNALSETUPS=.FALSE.
!
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('SETUP_READ')
      END IF
!
      CALL ATOMTYPELIST$NAME(THIS%I,THIS%ID)
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETR8('M',THIS%M)
      CALL ATOMTYPELIST$GETR8('ZV',THIS%ZV)
      CALL ATOMTYPELIST$GETR8('PS<G2>',THIS%PSG2)
      CALL ATOMTYPELIST$GETR8('PS<G4>',THIS%PSG4)
      CALL ATOMTYPELIST$GETCH('SOFTCORETYPE',THIS%SOFTCORETYPE)
!
      CALL FILEHANDLER$UNIT(THIS%FILEID,NFIL)
      CALL SETUPREAD$NEW(NFIL,TNEWFORMAT)
!
      IF(TNEWFORMAT) THEN
        CALL SETUPREAD$GETI4('LNX',LNX)
        CALL SETUPREAD$GETI4('GID',THIS%GID)
        GID=THIS%GID
        CALL SETUPREAD$GETI4('NR',NR)
        NRX=NR
        THIS%LNX=LNX
      ELSE
        CALL INPOT$GRID(NFIL,R1,DEX,NR)
        CALL RADIAL$NEW('LOG',GID)
        THIS%GID=GID
        NRX=NR
        CALL RADIAL$SETI4(GID,'NR',NR)
        CALL RADIAL$SETR8(GID,'DEX',DEX)
        CALL RADIAL$SETR8(GID,'R1',R1)
        CALL INPOT$LNX(NFIL,LNX)
        THIS%LNX=LNX
      END IF
      ALLOCATE(THIS%LOX(LNX))
      ALLOCATE(THIS%VADD(NRX))
      ALLOCATE(THIS%AECORE(NRX))
      ALLOCATE(THIS%PSCORE(NRX))
      ALLOCATE(THIS%PRO(NRX,LNX))
      ALLOCATE(THIS%AEPHI(NRX,LNX))
      ALLOCATE(THIS%PSPHI(NRX,LNX))
      ALLOCATE(THIS%UPHI(NRX,LNX))
      ALLOCATE(THIS%TUPHI(NRX,LNX))
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%DOVER(LNX,LNX))
      THIS%VADD=0.D0
      THIS%AECORE=0.D0
      THIS%PSCORE=0.D0
      THIS%PRO=0.D0
      THIS%AEPHI=0.D0
      THIS%PSPHI=0.D0
!     
!     ==================================================================
!     ==  READ PSEUDOPOTENTIALS AND PSEUDO WAVE FUNCTIONS             ==
!     ==================================================================
                            CALL TRACE$PASS('READ SETUP FILES')
      THIS%RCBG=1.D0/SQRT(0.218D0)
      
!      CALL INPOT$READALL(NFIL,NRX,R1,DEX,NR,THIS%LNX,THIS%LOX &
!     &         ,THIS%AEZ,THIS%PSZ,THIS%PSPHI,THIS%AEPHI &
!     &         ,THIS%VADD,THIS%RCSM,THIS%DTKIN,THIS%DOVER &
!     &         ,IRCCOR,THIS%AECORE,THIS%PSCORE,THIS%PRO)
      IF(TNEWFORMAT) THEN
        CALL SETUPREAD$GETI4A('LOX',LNX,THIS%LOX)
        CALL SETUPREAD$GETR8('AEZ',THIS%AEZ)
        CALL SETUPREAD$GETR8A('PSPHI',NR*LNX,THIS%PSPHI)
        CALL SETUPREAD$GETR8A('AEPHI',NR*LNX,THIS%AEPHI)
        CALL SETUPREAD$GETR8A('PRO',NR*LNX,THIS%PRO)
        CALL SETUPREAD$GETR8A('NDLSPHI',NR*LNX,THIS%UPHI)
        CALL SETUPREAD$GETR8A('NDLSTPHI',NR*LNX,THIS%TUPHI)
        CALL SETUPREAD$GETR8A('VADD',NR,THIS%VADD)
        CALL SETUPREAD$GETR8('RCSM',THIS%RCSM)
        CALL SETUPREAD$GETR8A('DT',LNX*LNX,THIS%DTKIN)
        CALL SETUPREAD$GETR8A('DO',LNX*LNX,THIS%DOVER)
        CALL SETUPREAD$GETR8A('PSCORE',NR,THIS%PSCORE)
        CALL SETUPREAD$GETR8A('AECORE',NR,THIS%AECORE)
        CALL SETUPREAD$GETR8A('AEPOT',NR,THIS%ATOM%AEPOT)
      ELSE
        CALL INPOT$READALL(NFIL,NRX,R1,DEX,NR,THIS%LNX,THIS%LOX &
     &         ,THIS%AEZ,PSZ,THIS%PSPHI,THIS%AEPHI &
     &         ,THIS%VADD,THIS%RCSM,THIS%DTKIN,THIS%DOVER &
     &         ,THIS%AECORE,THIS%PSCORE,THIS%PRO)
      END IF
      CALL FILEHANDLER$CLOSE(THIS%FILEID)
PRINT*,'NEW FORMAT?',TNEWFORMAT
ALLOCATE(R(NR))
CALL RADIAL$R(GID,NR,R)
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*THIS%AECORE*Y0*R**2,SVAR)
PRINT*,'INT AECORE ',SVAR
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*THIS%PSCORE*Y0*R**2,SVAR)
PRINT*,'INT PSCORE ',SVAR
PRINT*,'AEZ ',THIS%AEZ
PRINT*,'RCSM ',THIS%RCSM
!STOP
!     
!     ==================================================================
!     == LIMIT NUMBER OF PROJECTORS FOR EACH L                        ==
!     ==================================================================
      LX=0
      DO LN=1,THIS%LNX
        LX=MAX(LX,THIS%LOX(LN))
      ENDDO
      ALLOCATE(NPRO(LX+1)) 
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETI4A('NPRO',LX+1,NPRO)
      DO L=0,LX
        ISVAR=0
        DO LN=1,THIS%LNX
          IF(THIS%LOX(LN).NE.L)CYCLE
          ISVAR=ISVAR+1
          IF(ISVAR.GT.NPRO(L+1)) THEN
            THIS%LOX(LN)=-1     ! MARK PROJECTORS TO BE DELETED BY LOX=-1
            ISVAR=ISVAR-1
          END IF
        ENDDO
        NPRO(L+1)=ISVAR
      ENDDO
      CALL ATOMTYPELIST$SETI4A('NPRO',LX+1,NPRO)
      DEALLOCATE(NPRO)
!
      LNOLD=THIS%LNX
      LNX=0
      DO LN=1,LNOLD
        IF(THIS%LOX(LN).NE.-1) LNX=LNX+1
      ENDDO
      THIS%LNX=LNX
!
!     == FOLD DOWN ARRAYS FOR PROJECTORS AND PARTIALWAVES, LOX =========
      ALLOCATE(DWORK(NRX,LNOLD,5))
      ALLOCATE(IWORK(LNOLD))
      DWORK(:,:,1)=THIS%PRO(:,:)
      DWORK(:,:,2)=THIS%AEPHI(:,:)
      DWORK(:,:,3)=THIS%PSPHI(:,:)
      DWORK(:,:,4)=THIS%UPHI(:,:)
      DWORK(:,:,5)=THIS%TUPHI(:,:)
      IWORK(:)=THIS%LOX(:)
      DEALLOCATE(THIS%PRO)
      DEALLOCATE(THIS%AEPHI)
      DEALLOCATE(THIS%PSPHI)
      DEALLOCATE(THIS%UPHI)
      DEALLOCATE(THIS%TUPHI)
      DEALLOCATE(THIS%LOX)
      ALLOCATE(THIS%PRO(NRX,THIS%LNX))
      ALLOCATE(THIS%AEPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%PSPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%UPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%TUPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%LOX(THIS%LNX))
      ISVAR=0
      DO LN=1,LNOLD
        IF(IWORK(LN).EQ.-1) CYCLE
        ISVAR=ISVAR+1
        THIS%PRO(:,ISVAR)=DWORK(:,LN,1)
        THIS%AEPHI(:,ISVAR)=DWORK(:,LN,2)
        THIS%PSPHI(:,ISVAR)=DWORK(:,LN,3)
        THIS%UPHI(:,ISVAR)=DWORK(:,LN,4)
        THIS%TUPHI(:,ISVAR)=DWORK(:,LN,5)
        THIS%LOX(ISVAR)=IWORK(LN)
      ENDDO
      DEALLOCATE(DWORK)
!
!     == FOLD DOWN ARRAYS FOR DTKIN AND DOVER ==========================
      ALLOCATE(DWORK(LNOLD,LNOLD,2))
      DWORK(:,:,1)=THIS%DTKIN(:,:)
      DWORK(:,:,2)=THIS%DOVER(:,:)
      DEALLOCATE(THIS%DTKIN)
      DEALLOCATE(THIS%DOVER)
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%DOVER(LNX,LNX))
      LN1A=0
      DO LN1=1,LNOLD
        IF(IWORK(LN1).EQ.-1) CYCLE
        LN1A=LN1A+1
        LN2A=0
        DO LN2=1,LNOLD
          IF(IWORK(LN2).EQ.-1) CYCLE
          LN2A=LN2A+1
          THIS%DTKIN(LN1A,LN2A)=DWORK(LN1,LN2,1)
          THIS%DOVER(LN1A,LN2A)=DWORK(LN1,LN2,2)
        ENDDO
      ENDDO
      DEALLOCATE(DWORK)
!
      DEALLOCATE(IWORK)

!!$PRINT*,'=============================================================='
!!$PRINT*,'AEZ  ',THIS%AEZ
!!$PRINT*,'RCSM ',THIS%RCSM
!!$PRINT*,'LNX  ',THIS%LNX
!!$PRINT*,'LOX  ',THIS%LOX
!!$PRINT*,'PRO  ',THIS%PRO
!!$PRINT*,'AEPHI',THIS%AEPHI
!!$PRINT*,'PSPHI',THIS%PSPHI
!!$PRINT*,'VADD ',THIS%VADD
!!$PRINT*,'PSCORE',THIS%PSCORE
!!$PRINT*,'AECORE',THIS%AECORE
!!$PRINT*,'DOVER',THIS%DOVER
!!$PRINT*,'DTKIN',THIS%DTKIN
!CALL ERROR$STOP('FORCED STOP IN SETUP')
!     
!     ==================================================================
!     == SET VALUES BEYOND A CERTAIN RADIUS EXACTLY TO ZERO           ==
!     ==================================================================
                            CALL TRACE$PASS('CHECK MAX. RADIUS')
      CALL RADIAL$R(GID,NR,R)
      IRMAX=0
      DO IR=1,NR
        TCHK=(ABS(THIS%VADD(IR)).LT.TOL)
        TCHK=TCHK.AND.(ABS(THIS%PSCORE(IR)-THIS%AECORE(IR)).LT.TOL)
        DO LN=1,THIS%LNX
          TCHK=TCHK.AND. &
     &           (ABS(THIS%AEPHI(IR,LN)-THIS%PSPHI(IR,LN)).LT.TOL)
        ENDDO
!       == LDAPLUSU REQUIRES A SOMEWHAT LARGER RADIUS ==================
        TCHK=TCHK.AND.(R(IR).GE.6.D0)  
        IF(.NOT.TCHK) IRMAX=IR
      ENDDO
      DO IR=IRMAX+1,NR
        THIS%VADD(IR)=0.D0
        DO LN=1,THIS%LNX
          THIS%AEPHI(IR,LN)=0.D0
          THIS%PSPHI(IR,LN)=0.D0
          THIS%UPHI(IR,LN)=0.D0
          THIS%TUPHI(IR,LN)=0.D0
        ENDDO
      ENDDO
!     
!     ================================================================
!     ==  DEFINE ARRAYS                                             ==
!     ================================================================
                            CALL TRACE$PASS('DEFINE ARRAYS')
!
!     == SELECT NATURAL VALUES =======================================
      THIS%LMNX=0
      THIS%LMRX=0
      DO LN=1,THIS%LNX
        THIS%LMNX=THIS%LMNX+2*THIS%LOX(LN)+1
        THIS%LMRX=MAX(THIS%LMRX,(2*THIS%LOX(LN)+1)**2)
      ENDDO
!
!     == LIMIT MAX ANGULAR MOMENTUM FOR THE DENSITY TO MAX VALUE =======
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETI4('LRHOX',LMRXCUT)
      LMRXCUT=(LMRXCUT+1)**2
      THIS%LMRX=MIN(THIS%LMRX,LMRXCUT)
      CALL ATOMTYPELIST$UNSELECT
!     
!     ==================================================================
!     ==  UPDATE GLOBAL VARIABLES                                     ==
!     ==================================================================
      LMNXX=MAX(LMNXX,THIS%LMNX)
      LMRXX=MAX(LMRXX,THIS%LMRX)
      LNXX=MAX(LNXX,THIS%LNX)
!     
!     ==================================================================
!     ==  PERFORM BESSELTRANSFORMS                                    ==
!     ==================================================================
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      CALL RADIAL$GETR8(GID,'DEX',DEX)
      G1=GMAX*EXP(-DEX*DBLE(NG-1))
      CALL RADIAL$NEW('LOG',GIDG)
      THIS%GIDG=GIDG
      CALL RADIAL$SETI4(GIDG,'NR',NG)
      CALL RADIAL$SETR8(GIDG,'R1',G1)
      CALL RADIAL$SETR8(GIDG,'DEX',DEX)
PRINT*,'GIDG ',GIDG,G1,DEX,NG
!       
!     == VADD (VBAR) ===================================================
      ALLOCATE(THIS%VADDOFG(NG))
      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%VADD,GIDG,NG,THIS%VADDOFG)
      THIS%VADDOFG(:)=FOURPI*THIS%VADDOFG(:)
!     == PSCORE (VBAR) =================================================
      ALLOCATE(THIS%PSCOREOFG(NG))
      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%PSCORE,GIDG,NG,THIS%PSCOREOFG)
      THIS%PSCOREOFG(:)=FOURPI*THIS%PSCOREOFG(:)
!     == PROJECTORS ====================================================
      ALLOCATE(THIS%PROOFG(NG,LNX))
      DO LN=1,LNX
        L=THIS%LOX(LN)
        CALL RADIAL$BESSELTRANSFORM(L,GID,NR,THIS%PRO(:,LN),GIDG,NG,THIS%PROOFG(:,LN))
        THIS%PROOFG(:,LN)=FOURPI*THIS%PROOFG(:,LN)
      ENDDO
!     == COMPENSATION GAUSSIAN =========================================
      ALLOCATE(THIS%NHATPRIMEOFG(NG))
      ALLOCATE(THIS%VHATOFG(NG))
      CALL SETUP_COMPOFG(THIS%RCBG,THIS%RCSM,GIDG,NG &
     &                  ,THIS%NHATPRIMEOFG,THIS%VHATOFG)
!      
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_READ_NEW
!     **************************************************************************
!     **  READ SELECTED SETUP                                                 **
!     **  REQUIRES INFORMATION FROM ATOMTYPELIST                              **
!     **    NAME; LRHOX                                                       **
!     **  REQUIRES THE FILEHANDLER TO KNOW THE SETUP FILE                     **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: NBX=19
      INTEGER(4)            :: LOFI(NBX)
      INTEGER(4)            :: SOFI(NBX)
      REAL(8)               :: FOFI(NBX)
      INTEGER(4)            :: NNOFI(NBX)
      REAL(8)               :: EOFI(NBX)
      INTEGER(4)            :: GID
      INTEGER(4)            :: GIDG
      REAL(8)               :: DEX
      REAL(8)               :: G1
      REAL(8)   ,PARAMETER  :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
      INTEGER(4),PARAMETER  :: NG=250
      INTEGER(4)            :: NR
      INTEGER(4)            :: LX,LNX
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC 
      INTEGER(4),ALLOCATABLE:: NPRO(:)
      REAL(8)   ,ALLOCATABLE:: R(:)
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)   ,ALLOCATABLE:: RC(:)
      REAL(8)   ,ALLOCATABLE:: LAMBDA(:)  ! SECOND PARAMETER FOR PSEUDIZATION
      REAL(8)               :: RBOX
      REAL(8)               :: ROUT
      REAL(8)               :: AEZ
      REAL(8)               :: ZV
      REAL(8)               :: ETOT
      CHARACTER(32)         :: KEY
      REAL(8)   ,ALLOCATABLE:: PSI(:,:)
      LOGICAL(4)            :: TCHK
      REAL(8)               :: PI,FOURPI,Y0,C0LL
      REAL(8)               :: SVAR
      INTEGER(4)            :: IR,IB,LN,L
      INTEGER(4)            :: LRHOX
      CHARACTER(64)         :: STRING
      CHARACTER(64)         :: PSEUDIZATION
      INTEGER(4)            :: NFIL
!     **************************************************************************
                            CALL TRACE$PUSH('SETUP_READ_NEW')
!
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      C0LL=1.D0/SQRT(FOURPI)
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('SETUP_READ')
      END IF
!
      CALL ATOMTYPELIST$NAME(THIS%I,THIS%ID)
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETR8('M',THIS%M)
      CALL ATOMTYPELIST$GETR8('PS<G2>',THIS%PSG2)
      CALL ATOMTYPELIST$GETR8('PS<G4>',THIS%PSG4)
      CALL ATOMTYPELIST$GETCH('SOFTCORETYPE',THIS%SOFTCORETYPE)
!
!     ==========================================================================
!     == EXTRACT LNX,LOX,LX FROM NPRO AS DEFINED IN STRC INPUT FILE           ==
!     ==========================================================================
      ALLOCATE(NPRO(10)) 
      CALL ATOMTYPELIST$GETI4A('NPRO',10,NPRO)
      LNX=SUM(NPRO)
      ALLOCATE(LOX(LNX))
      LN=0
      DO L=0,SIZE(NPRO)-1
        DO WHILE (NPRO(L+1).GT.0)
          LN=LN+1
          NPRO(L+1)=NPRO(L+1)-1
          LOX(LN)=L
        ENDDO
      ENDDO
      LX=MAXVAL(LOX)
      DEALLOCATE(NPRO)
      THIS%LMNX=SUM(2*LOX(:)+1)
      THIS%LNX=LNX
      ALLOCATE(THIS%LOX(LNX))
      THIS%LOX(:)=LOX(:LNX)
!
!     == LIMIT LRHOX TO THE MAXIMUM CONSISTENT WITH WAVE FUNCTION CUTOFF
      CALL ATOMTYPELIST$GETI4('LRHOX',LRHOX)
      LRHOX=MIN(2*LX,LRHOX)
      CALL ATOMTYPELIST$SETI4('LRHOX',LRHOX)
      THIS%LMRX=(LRHOX+1)**2
!     
!     ==  UPDATE GLOBAL VARIABLES ==============================================
      LMNXX=MAX(LMNXX,THIS%LMNX)
      LMRXX=MAX(LMRXX,THIS%LMRX)
      LNXX=MAX(LNXX,THIS%LNX)
!
!     ==========================================================================
!     == READ SETUP INFORMATION FROM PARAMETER FILE                           ==
!     ==========================================================================
      ALLOCATE(THIS%PARMS%RCL(LX+1))
      ALLOCATE(THIS%PARMS%LAMBDA(LX+1))
      CALL ATOMTYPELIST$GETCH('ID',THIS%PARMS%ID)
      CALL ATOMLIB$SCNTLLOOKUP(THIS%PARMS%ID,GID,AEZ,ZV,RBOX,LX &
     &             ,THIS%PARMS%TYPE,THIS%PARMS%RCL,THIS%PARMS%LAMBDA &
     &             ,THIS%RCSM &
     &             ,THIS%PARMS%POW_POT,THIS%PARMS%RC_POT,this%parms%tval0_pot &
     &             ,THIS%PARMS%VAL0_POT &
     &             ,THIS%PARMS%POW_CORE,THIS%PARMS%RC_CORE &
     &             ,this%parms%tval0_core,THIS%PARMS%VAL0_CORE)
      THIS%GID=GID
      THIS%AEZ=AEZ
      THIS%ZV=ZV
      THIS%RCBG=1.D0/SQRT(0.218D0)
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      ROUT=R(NR-3)
      DEALLOCATE(R)
!     
!     ==========================================================================
!     ==  PERFORM ALL-ELECTRON CALCULATION FOR THE ATOM IN A BOX              ==
!     ==========================================================================
      ALLOCATE(THIS%ATOM%AEPOT(NR))
      ALLOCATE(PSI(NR,NBX))
      KEY='START,REL,NONSO'
      CALL ATOMLIB$AESCF(GID,NR,KEY,ROUT,AEZ,NBX,NB,LOFI,SOFI,FOFI,NNOFI &
    &                   ,ETOT,THIS%ATOM%AEPOT,EOFI,PSI)
!
!     == MAP ATOMIC DATA ON GRID ===============================================
      THIS%ATOM%NB=NB
      ALLOCATE(THIS%ATOM%LOFI(NB))
      ALLOCATE(THIS%ATOM%FOFI(NB))
      ALLOCATE(THIS%ATOM%EOFI(NB))
      ALLOCATE(THIS%ATOM%AEPSI(NR,NB))
      THIS%ATOM%LOFI=LOFI(1:NB)
      THIS%ATOM%FOFI=FOFI(1:NB)
      THIS%ATOM%EOFI=EOFI(1:NB)
      THIS%ATOM%AEPSI=PSI(:,1:NB)
!
!     == DETERMINE NUMBER OF CORE SHELLS =======================================
      SVAR=AEZ-ZV
      NC=0
      IF(SVAR.GT.1.D-3) THEN
        DO IB=1,NB
          SVAR=SVAR-FOFI(IB)
          NC=IB
          IF(ABS(SVAR).LT.1.D-3) EXIT
          IF(SVAR.LT.-1.D-3) THEN
            CALL ERROR$MSG('ZV INCONSISTENT WITH COMPLETE ANGULAR MOMENTUM SHELLS')
            CALL ERROR$R8VAL('AEZ',AEZ)
            CALL ERROR$R8VAL('ZV',ZV)
            CALL ERROR$R8VAL('SVAR',SVAR)
            CALL ERROR$STOP('SETUP_READ_NEW')
          END IF
        ENDDO
      END IF
      THIS%ATOM%NC=NC
!
!     ==========================================================================
      ALLOCATE(THIS%ISCATT(LNX))
      CALL SETUP_MAKEISCATT(AEZ,NB,NC,LOFI(1:NB),NNOFI(1:NB),LNX,LOX &
     &                     ,THIS%ISCATT)
!
!     ==========================================================================
!     == CALCULATE AND PSEUDIZE CORE DENSITY                                  ==
!     ==========================================================================
      ALLOCATE(THIS%AECORE(NR))
      ALLOCATE(THIS%PSCORE(NR))
      THIS%AECORE(:)=0.D0
      DO IB=1,NC
        THIS%AECORE(:)=THIS%AECORE(:)+FOFI(IB)*PSI(:,IB)**2*C0LL
      ENDDO
      CALL ATOMIC_PSEUDIZE(GID,NR,THIS%PARMS%POW_CORE,this%parms%tval0_core &
     &         ,THIS%PARMS%VAL0_CORE,THIS%PARMS%RC_CORE,THIS%AECORE,THIS%PSCORE)
      DEALLOCATE(PSI)
!
!     ==========================================================================
!     == CONSTRUCT PARTIAL WAVES                                              ==
!     ==========================================================================
      ALLOCATE(THIS%VADD(NR))
      ALLOCATE(THIS%PSPOT(NR))
      ALLOCATE(THIS%PRO(NR,LNX))
      ALLOCATE(THIS%AEPHI(NR,LNX))
      ALLOCATE(THIS%PSPHI(NR,LNX))
      ALLOCATE(THIS%UPHI(NR,LNX))
      ALLOCATE(THIS%TUPHI(NR,LNX))
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%DOVER(LNX,LNX))
      ALLOCATE(THIS%EOFLN(LNX))
      ALLOCATE(THIS%ESCATT(LNX))
      ALLOCATE(THIS%NLPHIDOT(NR,LNX))
      ALLOCATE(THIS%PSPHIDOT(NR,LNX))
      ALLOCATE(THIS%AEPHIDOT(NR,LNX))
      THIS%VADD=0.D0
      THIS%PRO=0.D0
      THIS%AEPHI=0.D0
      THIS%PSPHI=0.D0
      THIS%RBOX=RBOX  ! OUTER NODE OF PARTIAL WAVES
!
      ALLOCATE(RC(LNX))
      ALLOCATE(LAMBDA(LNX))
      DO LN=1,LNX
        L=LOX(LN)
        RC(LN)=THIS%PARMS%RCL(L+1)
        LAMBDA(LN)=THIS%PARMS%LAMBDA(L+1)
      ENDDO
!
      CALL ATOMIC_MAKEPARTIALWAVES(GID,NR,KEY,AEZ,THIS%ATOM%AEPOT,NB,NC &
     &           ,LOFI(1:NB),SOFI(1:NB),NNOFI(1:NB),EOFI(1:NB),FOFI(1:NB) &
     &           ,RBOX,ROUT,LNX,LOX,THIS%PARMS%TYPE,RC,LAMBDA &
     &           ,THIS%ISCATT,THIS%EOFLN,THIS%ESCATT &
     &           ,THIS%AEPHI,THIS%PSPHI,THIS%UPHI,THIS%PRO,THIS%DTKIN,THIS%DOVER &
     &           ,THIS%PSPOT,THIS%PARMS%POW_POT,this%parms%tval0_pot &
     &           ,THIS%PARMS%VAL0_POT,THIS%PARMS%RC_POT &
     &           ,THIS%RCSM,THIS%VADD,THIS%NLPHIDOT,THIS%AEPHIDOT,THIS%PSPHIDOT)
!     
!     ==========================================================================
!     ==  calculate and print scattering properties                           ==
!     ==========================================================================
      WRITE(STRING,FMT='(F3.0)')AEZ
      STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
      CALL FILEHANDLER$SETFILE('TMP',.FALSE.,-'SCATT'//TRIM(STRING))
      CALL FILEHANDLER$SETSPECIFICATION('TMP','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('TMP','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('TMP','FORM','FORMATTED')
      CALL FILEHANDLER$UNIT('TMP',NFIL)
      CALL SETUP_TESTSCATTERING(NFIL)
      CALL FILEHANDLER$CLOSE('TMP')
!     
!     ==========================================================================
!     ==  PERFORM BESSELTRANSFORMS                                            ==
!     ==========================================================================
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      CALL RADIAL$GETR8(GID,'DEX',DEX)
      G1=GMAX*EXP(-DEX*DBLE(NG-1))
G1=1.175D-4     
DEX=0.05D0
      CALL RADIAL$NEW('LOG',GIDG)
      THIS%GIDG=GIDG
      CALL RADIAL$SETI4(GIDG,'NR',NG)
      CALL RADIAL$SETR8(GIDG,'R1',G1)
      CALL RADIAL$SETR8(GIDG,'DEX',DEX)
PRINT*,'GIDG ',GIDG,G1,DEX,NG
!       
!     == VADD (VBAR) ===========================================================
      ALLOCATE(THIS%VADDOFG(NG))

      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%VADD,GIDG,NG,THIS%VADDOFG)
      THIS%VADDOFG(:)=FOURPI*THIS%VADDOFG(:)
!     == PSCORE (VBAR) =========================================================
      ALLOCATE(THIS%PSCOREOFG(NG))
      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%PSCORE,GIDG,NG,THIS%PSCOREOFG)
      THIS%PSCOREOFG(:)=FOURPI*THIS%PSCOREOFG(:)
!     == PROJECTORS ============================================================
      ALLOCATE(THIS%PROOFG(NG,LNX))
      DO LN=1,LNX
        L=THIS%LOX(LN)
        CALL RADIAL$BESSELTRANSFORM(L,GID,NR,THIS%PRO(:,LN),GIDG,NG &
     &                             ,THIS%PROOFG(:,LN))
        THIS%PROOFG(:,LN)=FOURPI*THIS%PROOFG(:,LN)
      ENDDO
!     == COMPENSATION GAUSSIAN =================================================
      ALLOCATE(THIS%NHATPRIMEOFG(NG))
      ALLOCATE(THIS%VHATOFG(NG))
      CALL SETUP_COMPOFG(THIS%RCBG,THIS%RCSM,GIDG,NG &
     &                  ,THIS%NHATPRIMEOFG,THIS%VHATOFG)
!      
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$REPORT(NFIL)
!     **************************************************************************
!     **  REPORT SETUP INFORMATION                                            **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NFIL
      INTEGER(4)                 :: L,NPRO,LN,NPROSUM
      TYPE(THIS_TYPE),POINTER    :: THIS1
      CHARACTER(64)              :: STRING
      REAL(8)                    :: U,PI,Y0
      INTEGER(4)                 :: THISTASK,NTASKS
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS$GET('U',U)
      THIS1=>FIRST
      DO 
        WRITE(NFIL,*)
        CALL REPORT$TITLE(NFIL,'ATOMIC SETUP: '//TRIM(THIS1%ID))
        CALL REPORT$R8VAL(NFIL,'ATOMIC NUMBER',THIS1%AEZ,' ')
        CALL REPORT$R8VAL(NFIL,'ATOMIC MASS  ',THIS1%M/U,'U')
        CALL REPORT$R8VAL(NFIL,'#(VALENCE ELECTRONS)',THIS1%ZV,' ')
        L=0
        NPROSUM=0
        DO WHILE (NPROSUM.LT.THIS1%LNX)
          NPRO=0
          DO LN=1,THIS1%LNX
            IF(THIS1%LOX(LN).EQ.L) NPRO=NPRO+1
          ENDDO
          IF(NPRO.NE.0) THEN
             WRITE(STRING,*)L
             STRING='#(PROJECTORS FOR L='//ADJUSTL(STRING)
             STRING=TRIM(STRING)//')'
             CALL REPORT$I4VAL(NFIL,TRIM(STRING),NPRO,' ')
          END IF
          L=L+1
          NPROSUM=NPROSUM+NPRO
        ENDDO        

        CALL REPORT$I4VAL(NFIL,'MAX #(ANGULAR MOMENTA (L,M) FOR 1C-DENSITY)' &
     &                        ,THIS1%LMRX,' ')
        CALL REPORT$R8VAL(NFIL,'GAUSSIAN DECAY FOR COMPENSATION DENSITY ' &
     &                        ,THIS1%RCSM,'ABOHR ')
        CALL REPORT$R8VAL(NFIL,'GAUSSIAN DECAY FOR EXTENDED COMPENSATION DENSITY' &
     &                        ,THIS1%RCBG,'ABOHR ')
        CALL REPORT$R8VAL(NFIL,'PARAMETER PS<G2> FOR MASS RENORMALIZATION' &
     &                        ,THIS1%PSG2,'H')
        CALL REPORT$R8VAL(NFIL,'PARAMETER PS<G4> FOR MASS RENORMALIZATION' &
     &                        ,THIS1%PSG4,'H')
!
!       == REPORT SETTINGS FOR THE SETUP CONSTRUCTION ==========================
        IF(THIS1%PARMS%POW_POT.NE.0.D0) THEN
          CALL REPORT$CHVAL(NFIL,'SETUP ID',THIS1%PARMS%ID)
          IF(THIS1%PARMS%TYPE.EQ.'KERKER') THEN
            CALL REPORT$CHVAL(NFIL,'CONSTRUCTION METHOD','KERKER TYPE')
          ELSE IF(THIS1%PARMS%TYPE.EQ.'HBS') THEN
            CALL REPORT$CHVAL(NFIL,'CONSTRUCTION METHOD','HAMANN-BACHELET-SCHLUETER TYPE')
          ELSE
            CALL REPORT$CHVAL(NFIL,'CONSTRUCTION METHOD','UNKNOWN')
          END IF
          DO L=0,MAXVAL(THIS1%LOX)
            WRITE(STRING,*)L
            STRING='PARTIAL WAVE PSEUDIZATION PARAMETER RC FOR L='//ADJUSTL(STRING)
            CALL REPORT$R8VAL(NFIL,TRIM(STRING),THIS1%PARMS%RCL(L+1),'ABOHR')
            IF(THIS1%PARMS%TYPE.EQ.'HBS') THEN
              WRITE(STRING,*)L
              STRING='PARTIAL WAVE PSEUDIZATION PARAMETER LAMBDA FOR L='//ADJUSTL(STRING)
              CALL REPORT$R8VAL(NFIL,TRIM(STRING),THIS1%PARMS%LAMBDA(L+1),'ABOHR')
            END IF 
         ENDDO
          if(this%parms%tval0_pot) then
            CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: VAL(0)' &
    &                            ,THIS1%PARMS%VAL0_POT*Y0,'H')
          else
            CALL REPORT$CHVAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: VAL(0)' &
    &                            ,'NOT DETERMINED')
          end if
          CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: RC' &
    &                            ,THIS1%PARMS%RC_POT,'ABOHR')
          CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: POWER' &
    &                            ,THIS1%PARMS%POW_POT,'')
          if(this%parms%tval0_core) then
            CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: VAL(0)' &
    &                            ,THIS1%PARMS%VAL0_CORE*Y0,'H')
          else
            CALL REPORT$CHVAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: VAL(0)' &
    &                            ,'NOT DETERMINED')
          end if
          CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: RC' &
    &                            ,THIS1%PARMS%RC_CORE,'ABOHR')
          CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: POWER' &
    &                            ,THIS1%PARMS%POW_CORE,'')
        END IF
        IF(.NOT.ASSOCIATED(THIS1%NEXT)) EXIT
        THIS1=>THIS1%NEXT
      ENDDO
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$READ()
!     ******************************************************************
!     **  READ SETUP                                                  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)               :: ISP,NSP1
      CHARACTER(32)            :: NAME
      LOGICAL     ,ALLOCATABLE :: TNEW(:) ! SWITCH BETWEEN INTERNAL ATOM 
                                          ! PROGRAM AND READING THE FILE
!     ******************************************************************
                            CALL TRACE$PUSH('SETUP$READ')
!
!     ==================================================================
!     ==  CREATE SETUPS                                               ==
!     ==================================================================
      CALL ATOMTYPELIST$LENGTH(NSP1)
      ALLOCATE(TNEW(NSP1))
      DO ISP=1,NSP1
        CALL ATOMTYPELIST$NAME(ISP,NAME)
        CALL ATOMTYPELIST$SELECT(NAME)
        CALL ATOMTYPELIST$GETL4('TID',TNEW(ISP))
        CALL ATOMTYPELIST$UNSELECT
        CALL SETUP$SELECT(NAME)
      ENDDO
!
!     ==================================================================
!     ==  READ SETUP FILES                                            ==
!     ==================================================================
      DO ISP=1,NSP1
        CALL SETUP$ISELECT(ISP)
        IF(TNEW(ISP)) THEN
          CALL SETUP_READ_NEW
        ELSE
          CALL SETUP_READ
        END IF
      ENDDO
      DEALLOCATE(TNEW)
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$AEPARTIALWAVES(ISP_,NRX_,LNX_,AEPHI_)
!LEGACY CODE! -> SETUP$GETR8A('AEPHI'
!     ******************************************************************
!     **  RETURN AE PARTIAL WAVES ON THE RADIAL GRID                  **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      INTEGER(4),INTENT(IN) :: LNX_
      REAL(8)   ,INTENT(OUT):: AEPHI_(NRX_,LNX_)
      INTEGER(4)            :: NR
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      IF(NRX_.NE.NR) THEN
        CALL ERROR$MSG('INCONSISTENT GRID SIZE')
        CALL ERROR$STOP('SETUP$AEPARTIALWAVES')
      END IF
      IF(LNX_.NE.THIS%LNX) THEN
        CALL ERROR$MSG('INCONSISTENT #(PARTIAL WAVES)')
        CALL ERROR$STOP('SETUP$AEPARTIALWAVES')
      END IF
      AEPHI_(:,:)=THIS%AEPHI(:,:)
      RETURN  
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$PSPARTIALWAVES(ISP_,NRX_,LNX_,PSPHI_)
!LEGACY CODE! -> SETUP$GETR8A('PSPHI'
!     ******************************************************************
!     **  RETURN PS PARTIAL WAVE ON A RADIAL GRID                     **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      INTEGER(4),INTENT(IN) :: LNX_
      REAL(8)   ,INTENT(OUT):: PSPHI_(NRX_,LNX_)
      INTEGER(4)            :: NR
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      IF(NRX_.NE.NR) THEN
        CALL ERROR$MSG('INCONSISTENT GRID SIZE')
        CALL ERROR$STOP('SETUP$AEPARTIALWAVES')
      END IF
      IF(LNX_.NE.THIS%LNX) THEN
        CALL ERROR$MSG('INCONSISTENT #(PARTIAL WAVES)')
        CALL ERROR$STOP('SETUP$AEPARTIALWAVES')
      END IF
      PSPHI_(:,:)=THIS%PSPHI(:,:)
      RETURN  
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$1COVERLAP(ISP_,LNXX_,DOVER_)
!LEGACY CODE! -> SETUP$GETR8A('DO'
!     ******************************************************************
!     **  RETURN 1-C- OVERLAP OF THE PARTIAL WAVES                    **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNXX_
      REAL(8)   ,INTENT(OUT):: DOVER_(LNXX_,LNXX_)
      INTEGER(4)            :: LN1,LN2
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO LN1=1,THIS%LNX
        DO LN2=1,THIS%LNX
          DOVER_(LN1,LN2)=THIS%DOVER(LN1,LN2)
        ENDDO
      ENDDO
      RETURN  
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$1CKINETIC(ISP_,LNXX_,DTKIN_)
!LEGACY CODE! -> SETUP$GETR8A('DEKIN'
!     ******************************************************************
!     **  RETURN 1-C- KINETIC ENERGY OVERLAP OF THE PARTIAL WAVES     **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNXX_
      REAL(8)   ,INTENT(OUT):: DTKIN_(LNXX_,LNXX_)
      INTEGER(4)            :: LN1,LN2
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO LN1=1,THIS%LNX
        DO LN2=1,THIS%LNX
          DTKIN_(LN1,LN2)=THIS%DTKIN(LN1,LN2)
        ENDDO
      ENDDO
      RETURN  
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$RCSM(ISP_,RCSM_)
!LEGACY CODE! -> SETUP$GETR8('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: RCSM_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      RCSM_=THIS%RCSM
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$RCBG(ISP_,RCBG_)
!LEGACY CODE! -> SETUP$GETR8('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: RCBG_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      RCBG_=THIS%RCBG
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LNX(ISP_,LNX_)
!LEGACY CODE! -> SETUP$GETI4('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LNX_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      LNX_=THIS%LNX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMNX(ISP_,LMNX_)
!LEGACY CODE! -> SETUP$GETI4('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LMNX_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      LMNX_=THIS%LMNX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMNXX(LMNXX_)
!LEGACY CODE! -> SETUP$GETI4('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT):: LMNXX_
!     ******************************************************************
      LMNXX_=LMNXX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LOFLN(ISP_,LNX_,LOX_)
!LEGACY CODE! -> SETUP$GETI4A('LOX'
!     ******************************************************************
!     **  RETURN NUMBER MAIN ANGULAR MOMENTUM OF PARTIAL WAVES        **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNX_
      INTEGER(4),INTENT(OUT):: LOX_(LNX_)
      INTEGER(4)            :: LN
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      IF(LNX_.GT.THIS%LNX) THEN
        CALL ERROR$MSG('LNX ON INPUT TOO SMALL')
        CALL ERROR$OVERFLOW('LNX_',LNX_,THIS%LNX)
        CALL ERROR$STOP('SETUP$LOFLN')
      END IF
      DO LN=1,THIS%LNX
        LOX_(LN)=THIS%LOX(LN)
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$NSPECIES(NSP_)
!LEGACY CODE! -> SETUP$GETI4
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NSP_
!     ******************************************************************
      NSP_=NSP
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMRXX(LMRXX_)
!LEGACY CODE! -> SETUP$GETI4
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ****************************************************************** 
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: LMRXX_
!     ******************************************************************
      LMRXX_=LMRXX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LNXX(LNXX_)
!LEGACY CODE! -> SETUP$GETI4
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: LNXX_
!     ******************************************************************
      LNXX_=LNXX
      RETURN  
      END
!
!     ..........................................................................
      SUBROUTINE SETUP$LMRX(ISP_,LMRX_)
!LEGACY CODE! -> SETUP$GETI4
!     **************************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY                  **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LMRX_
!     **************************************************************************
      CALL SETUP$ISELECT(ISP_)
      LMRX_=THIS%LMRX
      RETURN  
      END
!
!     ..........................................................................
      SUBROUTINE SETUP$AEZ(ISP_,AEZ_)
!LEGACY CODE! ->SETUP$GETR8
!     **************************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY                  **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: AEZ_
!     **************************************************************************
      CALL SETUP$ISELECT(ISP_)
      AEZ_=THIS%AEZ
      RETURN  
      END
!
!     ..........................................................................
      SUBROUTINE SETUP_COMPOFG(RCBG,RCSM,GIDG,NG,G0,V0)
!     **************************************************************************
!     **                                                                      **
!     **  COMPENSATION DENSITY AND POTENTIAL ON A RADIAL GRID                 **
!     **  IN G-SPACE                                                          **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RCBG
      REAL(8)   ,INTENT(IN) :: RCSM
      INTEGER(4),INTENT(IN) :: GIDG
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(OUT):: G0(NG)
      REAL(8)   ,INTENT(OUT):: V0(NG)
      REAL(8)   ,PARAMETER  :: EPSILONGAMMA=1.D-7
      REAL(8)               :: PI
      REAL(8)               :: SVAR1,SVAR2,SVAR3,SVAR4
      REAL(8)               :: BGGAUSS,SMGAUSS
      INTEGER(4)            :: IG
      REAL(8)               :: GARR(NG)
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      SVAR1=-0.25D0*RCBG**2
      SVAR2=-0.25D0*RCSM**2
      SVAR3=4.D0*PI
      SVAR4=-4.D0*PI*SVAR3
      CALL RADIAL$R(GIDG,NG,GARR)
      DO IG=1,NG
        BGGAUSS=EXP(SVAR1*GARR(IG)**2)
        SMGAUSS=EXP(SVAR2*GARR(IG)**2)
        G0(IG)=SVAR3*BGGAUSS
        V0(IG)=SVAR4*(BGGAUSS-SMGAUSS)/GARR(IG)**2
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP_CHECKGAUSS(CELLVOL,RC,TOL,GMAX,RMAX,TCHKR,TCHKG)
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN) :: CELLVOL
      REAL(8)    ,INTENT(IN) :: RC
      REAL(8)    ,INTENT(IN) :: TOL
      REAL(8)    ,INTENT(IN) :: GMAX
      REAL(8)    ,INTENT(IN) :: RMAX
      LOGICAL(4) ,INTENT(OUT):: TCHKR
      LOGICAL(4) ,INTENT(OUT):: TCHKG
      REAL(8)                :: CHECK
!     ******************************************************************
      TCHKG=.TRUE.
      TCHKR=.TRUE.
      CHECK=-1.D0/CELLVOL*EXP(-0.25D0*(RC*GMAX)**2)
      IF(ABS(CHECK).GT.TOL) TCHKG=.FALSE. 
      CHECK=-EXP(-(RMAX/RC)**2)
      IF(ABS(CHECK).GT.TOL) TCHKR=.FALSE. 
      RETURN
      END
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$LNX(NFIL,LNX)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4),INTENT(OUT) :: LNX
      REAL(8)                :: R1,DEX
      INTEGER(4)             :: NR
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$LNX')
      REWIND NFIL
      READ(NFIL,FMT='(F15.10,F10.5,2I4)')R1,DEX,NR,LNX
                              CALL TRACE$POP
      RETURN
      END
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$GRID(NFIL,R1,DEX,NR)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4),INTENT(OUT) :: NR
      REAL(8)   ,INTENT(OUT) :: R1
      REAL(8)   ,INTENT(OUT) :: DEX
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$GRID')
      REWIND NFIL
      READ(NFIL,FMT='(F15.10,F10.5,2I4)')R1,DEX,NR
                              CALL TRACE$POP
      RETURN
      END
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$NC(NFIL,NC)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4),INTENT(OUT) :: NC
      REAL(8)                :: R1,DEX,PSZ,AEZ,RCSM
      INTEGER(4)             :: NR,LNX
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$LNX')
      REWIND NFIL
      NC=0
      READ(NFIL,ERR=1000,FMT='(F15.10,F10.5,2I4,2F5.2,F15.12,I5)') &
     &               R1,DEX,NR,LNX,PSZ,AEZ,RCSM,NC
1000 CONTINUE
                              CALL TRACE$POP
      RETURN
      END
!
!     ...........................................INPOT..........................
      SUBROUTINE INPOT$READALL(NFIL,NRX,R1,DEX,NR,LNX,LOX,AEZ,PSZ &
     &         ,PSPHI,AEPHI,VADD,RCSM &
     &         ,DTKIN,DOVER,RHOCOR,PSCORR,PRO)
!     **************************************************************************
!     **                                                                      **
!     ******************************************* P.E. BLOECHL, 1991 ***********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)  :: NFIL
      INTEGER(4) ,INTENT(IN)  :: NRX
      REAL(8)    ,INTENT(IN)  :: R1
      REAL(8)    ,INTENT(IN)  :: DEX
      INTEGER(4) ,INTENT(IN)  :: NR
      INTEGER(4) ,INTENT(IN)  :: LNX
      REAL(8)    ,INTENT(OUT) :: AEZ
      REAL(8)    ,INTENT(OUT) :: PSZ
      REAL(8)    ,INTENT(OUT) :: RCSM
      REAL(8)    ,INTENT(OUT) :: VADD(NRX)
      REAL(8)    ,INTENT(OUT) :: PRO(NRX,LNX)
      INTEGER(4) ,INTENT(OUT) :: LOX(LNX)
!INTEGER(4) ,INTENT(OUT) :: IRCCOR
      REAL(8)    ,INTENT(OUT) :: DTKIN(LNX,LNX)
      REAL(8)    ,INTENT(OUT) :: DOVER(LNX,LNX)
      REAL(8)    ,INTENT(OUT) :: AEPHI(NRX,LNX)
      REAL(8)    ,INTENT(OUT) :: PSPHI(NRX,LNX)
      REAL(8)    ,INTENT(OUT) :: RHOCOR(NRX)
      REAL(8)    ,INTENT(OUT) :: PSCORR(NRX)
      REAL(8)                 :: R11,DEX1
      INTEGER(4)              :: NR1,LNX1,I,IR,LN1,LN2,LN
!     **************************************************************************
                              CALL TRACE$PUSH('INPOT$READALL')
      REWIND NFIL
!     == OPTIONAL VALUE NC FOLLOWING RCSM IS NOT READ 
      READ(NFIL,6000)R11,DEX1,NR1,LNX1,PSZ,AEZ,RCSM
!      READ(NFIL,6000)R11,DEX1,NR1,LNX1,PSZ,AEZ,RCSM,IRCCOR
6000  FORMAT(F15.10,F10.5,2I4,2F5.2,F15.12,I5)
      IF(R11.NE.R1.OR.DEX1.NE.DEX.OR.NR1.NE.NR) THEN
        CALL ERROR$MSG('ONLY ONE TYPE OF RADIAL GRID ALLOWED')
        CALL ERROR$STOP('INPOT')
      END IF


!      IF(IRCCOR.LE.0.OR.IRCCOR.GT.NR) THEN
!      PRINT*,'WARNING! NO MT-RADIUS SPECIFIED FOR ATOM WITH Z=',AEZ
!       IRCCOR=NR-2
!      END IF
      IF(LNX.NE.LNX1) THEN
        CALL ERROR$MSG('LNX OUT OF RANGE')
        CALL ERROR$I4VAL('LNX ON INPUT',LNX)
        CALL ERROR$I4VAL('LNX ON FILE',LNX1)
        CALL ERROR$STOP('INPOT')
      END IF
                              CALL TRACE$PASS('BEFORE LOX')
      READ(NFIL,6020)(LOX(I),I=1,LNX)
6020  FORMAT(14I5)
                              CALL TRACE$PASS('BEFORE VADD')
      READ(NFIL,6100)(VADD(IR),IR=1,NR)
!     ====  RHOCOR = CORE CHARGE DENSITY  ======================================
                              CALL TRACE$PASS('BEFORE RHOCOR')
      READ(NFIL,6100)(RHOCOR(IR),IR=1,NR)
!     ====  PSCORR = PSEUDO CORE CHARGE DENSITY ================================
                              CALL TRACE$PASS('BEFORE PSCORR')
      READ(NFIL,6100)(PSCORR(IR),IR=1,NR)
!     ====  DTKIN = <AEPHI|-DELTA/2|AEPHI> - <PSPHI|-DELTA/2|PSPHI> ============
                              CALL TRACE$PASS('BEFORE DTKIN')
      READ(NFIL,6100)((DTKIN(LN1,LN2),LN1=1,LNX),LN2=1,LNX)
!     PRINT*,'DTKIN ',(DTKIN(LN,LN),LN=1,LNX)
!     ====  DOVER = <AEPHI|AEPHI> - <PSPHI|PSPHI> ==============================
                              CALL TRACE$PASS('BEFORE DOVER')
      READ(NFIL,6100)((DOVER(LN1,LN2),LN1=1,LNX),LN2=1,LNX)
                              CALL TRACE$PASS('BEFORE PROJECTORS')
      DO LN=1,LNX
        READ(NFIL,6100)(PRO(IR,LN),IR=1,NR)
        READ(NFIL,6100)(AEPHI(IR,LN),IR=1,NR)
        READ(NFIL,6100)(PSPHI(IR,LN),IR=1,NR)
6100    FORMAT(SP,5E14.8)
      ENDDO
                              CALL TRACE$POP
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE ATOMLIB$READSCNTL(LL_SCNTL)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE) ,INTENT(OUT) :: LL_SCNTL
      CHARACTER(128)             :: SCNTLNAME='STP.CNTL'
      TYPE(LL_TYPE) ,SAVE        :: LL_SCNTL_SAVE
      INTEGER(4)                 :: NFIL
      LOGICAL(4)    ,SAVE        :: TINI=.FALSE.
!     **************************************************************************
      IF(TINI) THEN
        LL_SCNTL=LL_SCNTL_SAVE
        RETURN
      END IF
      TINI=.TRUE.
!!$      CALL FILEHANDLER$SETFILE(+'SCNTL',.FALSE.,-SCNTLNAME)
!!$      CALL FILEHANDLER$SETSPECIFICATION(+'SCNTL','STATUS','OLD')
!!$      CALL FILEHANDLER$SETSPECIFICATION(+'SCNTL','POSITION','REWIND')
!!$      CALL FILEHANDLER$SETSPECIFICATION(+'SCNTL','ACTION','READ')
!!$      CALL FILEHANDLER$SETSPECIFICATION(+'SCNTL','FORM','FORMATTED')
!!$      CALL FILEHANDLER$UNIT('SCNTL',NFIL)
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_SCNTL)
      CALL FILEHANDLER$UNIT('PARMS_STP',NFIL)
      CALL LINKEDLIST$READ(LL_SCNTL,NFIL,'~')
      LL_SCNTL_SAVE=LL_SCNTL
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$SCNTLLOOKUP(SETUPID,GID,AEZ,ZV,RBOX,LX,TYPE,RCL,LAMBDA &
     &                                ,RCSM,POW_POT,RC_POT,tval0_pot,VAL0_POT &
     &                                 ,POW_CORE,RC_CORE,tval0_core,VAL0_CORE)
!     **************************************************************************
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: SETUPID
      INTEGER(4)  ,INTENT(OUT):: GID
      REAL(8)     ,INTENT(OUT):: AEZ
      REAL(8)     ,INTENT(OUT):: ZV
      REAL(8)     ,INTENT(OUT):: RBOX
      INTEGER(4)  ,INTENT(IN) :: LX
      CHARACTER(*),INTENT(OUT):: TYPE
      REAL(8)     ,INTENT(OUT):: RCL(LX+1)
      REAL(8)     ,INTENT(OUT):: LAMBDA(LX+1)
      REAL(8)     ,INTENT(OUT):: RCSM
      REAL(8)     ,INTENT(OUT):: POW_POT
      REAL(8)     ,INTENT(OUT):: RC_POT
      logical(4)  ,intent(out):: tval0_pot
      REAL(8)     ,INTENT(OUT):: VAL0_POT
      REAL(8)     ,INTENT(OUT):: POW_CORE
      REAL(8)     ,INTENT(OUT):: RC_CORE
      logical(4)  ,intent(out):: tval0_core
      REAL(8)     ,INTENT(OUT):: VAL0_CORE
      INTEGER(4)              :: NFIL
      TYPE(LL_TYPE)           :: LL_SCNTL
      INTEGER(4)              :: ITH
      LOGICAL(4)              :: TCHK,TCHK1,TCHK2,TCHK3
      CHARACTER(128)          :: ID
      CHARACTER(2)            :: EL
      REAL(8)                 :: R1
      REAL(8)                 :: DEX
      REAL(8)                 :: SVAR,SVAR1,SVAR2
      INTEGER(4)              :: NR
      INTEGER(4)              :: LENG
      INTEGER(4)              :: L,LN
      INTEGER(4)              :: IZ
      REAL(8)     ,ALLOCATABLE:: R(:)
      REAL(8)     ,ALLOCATABLE:: RCL1(:)
      REAL(8)     ,ALLOCATABLE:: lambda1(:)
      REAL(8)                 :: RCOV
      REAL(8)                 :: DMIN,DMAX,RX
      REAL(8)                 :: PI,Y0
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL ATOMLIB$READSCNTL(LL_SCNTL)
      CALL LINKEDLIST$SELECT(LL_SCNTL,'~')
      CALL LINKEDLIST$SELECT(LL_SCNTL,'SCNTL')
      ITH=0
      DO 
        ITH=ITH+1
        CALL LINKEDLIST$EXISTL(LL_SCNTL,'SETUP',ITH,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('SETUP DATA NOT FOUND')
          CALL ERROR$I4VAL('ITH',ITH)
          CALL ERROR$CHVAL('SETUPID',SETUPID)
          CALL ERROR$STOP('ATOMLIB$READCNTL')
        END IF
        CALL LINKEDLIST$SELECT(LL_SCNTL,'SETUP',ITH)
        CALL LINKEDLIST$EXISTD(LL_SCNTL,'ID',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!SCNTL!SETUP:ID NOT SPECIFIED')
          CALL ERROR$STOP('ATOMLIB$READCNTL')
        END IF
        CALL LINKEDLIST$GET(LL_SCNTL,'ID',1,ID)
        IF(ID.EQ.SETUPID) THEN
!
!         ======================================================================
!         == COLLECT RADIAL GRID                                              ==
!         ======================================================================
          CALL LINKEDLIST$EXISTL(LL_SCNTL,'GRID',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP!GRID NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$SELECT(LL_SCNTL,'GRID')

          CALL LINKEDLIST$EXISTD(LL_SCNTL,'DMIN',1,TCHK)
          IF(TCHK) THEN
!           == R1 ==============================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'DMIN',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:DMIN NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'DMIN',1,DMIN)
!           == R1 ==============================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'DMAX',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:DMAX NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'DMAX',1,DMAX)
!           == R1 ==============================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'RMAX',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:RMAX NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'RMAX',1,RX)
!
            CALL RADIAL$GRIDPARAMETERS(DMIN,DMAX,RX,R1,DEX,NR)
!
!           == CHECK FOR CONFLICT ==============================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'R1',1,TCHK1)
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'DEX',1,TCHK2)
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'NR',1,TCHK3)
            IF(TCHK1.OR.TCHK2.OR.TCHK3) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:R1,DEX,NR AND DMIN,DMAX,RX')
              CALL ERROR$MSG('ARE NOT COMPATIBLE')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
          ELSE
!           == R1 ==============================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'R1',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:R1 NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'R1',1,R1)
!           == DEX =============================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'DEX',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:DEX NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'DEX',1,DEX)
!           === NR =============================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'NR',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:R NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'NR',1,NR)
          END IF
!
!         == CREATE GRID =======================================================
          CALL RADIAL$NEW('SHLOG',GID)
          CALL RADIAL$SETR8(GID,'R1',R1)
          CALL RADIAL$SETR8(GID,'DEX',DEX)
          CALL RADIAL$SETI4(GID,'NR',NR)
!
          CALL LINKEDLIST$SELECT(LL_SCNTL,'..')
!
!         ======================================================================
!         == IDENTIFY ATOM                                                    ==
!         ======================================================================
!         == COLLECT ATOMIC NUMBER =============================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'EL',1,TCHK)
          EL='XX'
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'EL',1,EL)
!
!         == COLLECT ATOMIC NUMBER =============================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'Z',1,TCHK)
          IF(TCHK) THEN
            IF(EL.EQ.'XX') THEN 
              CALL ERROR$MSG('!SCNTL!SETUP:AEZ AND EL="XX"')
              CALL ERROR$MSG('MUST NOT BE SPECIFIED SIMULTANEOUSLY')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'Z',1,AEZ)
          ELSE
            IF(EL.EQ.'XX') THEN
              CALL ERROR$MSG('!SCNTL!SETUP:AEZ NOT SPECIFIED')
              CALL ERROR$MSG('AND ELEMENT NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            ELSE
              CALL PERIODICTABLE$GET(EL,'Z',AEZ)
            END IF
          END IF
!
!         == COLLECT #VALENCE ELECTRONS ========================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'ZV',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP:ZV NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$GET(LL_SCNTL,'ZV',1,ZV)
!
!         == COLLECT BOX RADIUS ================================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RBOX/RCOV',1,TCHK)
          ALLOCATE(R(NR))
          CALL RADIAL$R(GID,NR,R)
          IF(.NOT.TCHK) THEN
            RBOX=R(NR-3)
          ELSE
            CALL LINKEDLIST$GET(LL_SCNTL,'RBOX/RCOV',1,RBOX)
            CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
            RBOX=RBOX*RCOV 
            IF(RBOX.GT.R(NR-3)) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:RBOX EXCEEDS GRID DIMENSIONS')
              CALL ERROR$R8VAL('AEZ ',AEZ)
              CALL ERROR$CHVAL('SETUPID ',SETUPID)
              CALL ERROR$R8VAL('RBOX ',RBOX)
              CALL ERROR$R8VAL('Rcov ',Rcov)
              CALL ERROR$R8VAL('R(NR-3)',R(NR-3))
              CALL ERROR$R8VAL('R(NR)',R(NR))
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
          END IF
          DEALLOCATE(R)
!
!         ======================================================================
!         == PARAMETER FOR PARTIALWAVE CONSTRUCTION                           ==
!         ======================================================================
          TYPE='KERKER'
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'TYPE',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'TYPE',1,TYPE)
!
!         == matching radius ===================================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RCL/RCOV',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP:RCL NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$SIZE(LL_SCNTL,'RCL/RCOV',1,LENG)
          ALLOCATE(RCL1(LENG))
          CALL LINKEDLIST$GET(LL_SCNTL,'RCL/RCOV',1,RCL1)
          IF(LENG.LT.LX+1) THEN
            CALL ERROR$MSG('!SCNTL!SETUP:RCL ARRAY TOO SHORT')
            CALL ERROR$I4VAL('MAX L ON FILE ',LENG-1)
            CALL ERROR$I4VAL('MAX L REQUESTED',LX)
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
          LENG=MIN(SIZE(RCL),LENG)
          RCL=RCL1(1:LENG)*RCOV
          deallocate(rcl1)
!
!         == additional lambda parameter required for 'hbs' type construction ==
          LAMBDA(:)=0.D0
          IF(TYPE.EQ.'HBS') THEN
            LAMBDA(:)=6.D0    !DEFAULT FOR HBS CONSTRUCTION
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'LAMBDA',1,TCHK)
            IF(TCHK) THEN
              CALL LINKEDLIST$SIZE(LL_SCNTL,'LAMBDA',1,LENG)
              ALLOCATE(LAMBDA1(LENG))
              CALL LINKEDLIST$GET(LL_SCNTL,'LAMBDA',1,LAMBDA1)
              IF(LENG.LT.LX+1) THEN
                CALL ERROR$MSG('!SCNTL!SETUP:LAMBDA ARRAY TOO SHORT')
                CALL ERROR$I4VAL('MAX L ON FILE ',LENG-1)
                CALL ERROR$I4VAL('MAX L REQUESTED',LX)
                CALL ERROR$STOP('ATOMLIB$READCNTL')
              END IF
              LENG=MIN(SIZE(LAMBDA),LENG)
              LAMBDA=LAMBDA1(1:LENG)
              deallocate(lambda1)
            END IF
          END IF
!
!         == decay constant for short-ranged compensation charge density ======       
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RCSM/RCOV',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_SCNTL,'RCSM/RCOV',1,RCSM)
            RCSM=RCSM*RCOV
          ELSE
            CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
            RCSM=0.25D0*RCOV
          END IF
!
!         ======================================================================
!         == DEFINE POTENTIAL PSEUDIZATION                                    ==
!         ======================================================================
          CALL LINKEDLIST$EXISTL(LL_SCNTL,'POT',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP!POT NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$SELECT(LL_SCNTL,'POT')
!  
          POW_POT=2.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'POW',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'POW',1,POW_POT)
!
          VAL0_POT=0.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'VAL0',1,TCHK)
          tval0_pot=tchk
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'VAL0',1,VAL0_POT)
          VAL0_POT=VAL0_POT/Y0
!
          RC_POT=1.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RC/RCOV',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'RC/RCOV',1,RC_POT)
          CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
          RC_POT=RC_POT*RCOV

          CALL LINKEDLIST$SELECT(LL_SCNTL,'..')
!
!         ======================================================================
!         == DEFINE CORE PSEUDIZATION                                         ==
!         ======================================================================
          CALL LINKEDLIST$EXISTL(LL_SCNTL,'CORE',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP!CORE NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$SELECT(LL_SCNTL,'CORE')
!  
          POW_CORE=2.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'POW',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'POW',1,POW_CORE)
!
          VAL0_CORE=0.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'VAL0',1,TCHK)
          TVAL0_CORE=TCHK
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'VAL0',1,VAL0_CORE)
          VAL0_CORE=VAL0_CORE/Y0
!
          RC_CORE=1.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RC/RCOV',1,TCHK)
          IF(TCHK) CALL LINKEDLIST$GET(LL_SCNTL,'RC/RCOV',1,RC_CORE)
          CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
          RC_CORE=RC_CORE*RCOV
!
          CALL LINKEDLIST$SELECT(LL_SCNTL,'..')
!
!         ======================================================================
!         == DONE                                                             ==
!         ======================================================================
          EXIT        
        END IF
        CALL LINKEDLIST$SELECT(LL_SCNTL,'..')
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMIC_MAKEPARTIALWAVES(GID,NR,KEY,AEZ,AEPOT &
     &                  ,NB,NC,LOFI,SOFI,NNOFI,EOFI,FOFI &
     &           ,RBOX,ROUT,LNX,LOX,TYPE,RC,LAMBDA,ISCATT,EOFLN,ESCATT &
     &                  ,AEPHI,PSPHI,NLPHI,PRO,DT,DOVER,PSPOT &
     &                  ,POW_POT,tval0_pot,VAL0_POT,RC_POT,RCSM,VADD &
     &                  ,NLPHIDOT,AEPHIDOT,PSPHIDOT)
!     **************************************************************************
!     **  CONSTRUCTS  THE SETUP                                               **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTEnT(IN) :: TYPE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      CHARACTER(*),INTENT(IN) :: KEY
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)   ! ALL ELECTRON POTENTIAL
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NC
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      INTEGER(4),INTENT(IN) :: SOFI(NB)
      INTEGER(4),INTENT(IN) :: NNOFI(NB)
      REAL(8)   ,INTENT(IN) :: EOFI(NB)
      REAL(8)   ,INTENT(IN) :: ROUT    ! HARD SPHERE FOR ATOM HAS RADIUS RBOX
      REAL(8)   ,INTENT(IN) :: RBOX    ! RADIUS FOR BOUNDARY CONDITION OF PHI
      REAL(8)   ,INTENT(IN) :: FOFI(NB)
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: ISCATT(LNX)
      REAL(8)   ,INTENT(IN) :: RC(LNX)
      REAL(8)   ,INTENT(IN) :: LAMBDA(LNX)
      REAL(8)   ,INTENT(IN) :: POW_POT
      logical(4),intent(in) :: tval0_pot
      REAL(8)   ,INTENT(IN) :: VAL0_POT
      REAL(8)   ,INTENT(IN) :: RC_POT
      REAL(8)   ,INTENT(IN) :: RCSM
      REAL(8)   ,INTENT(OUT):: VADD(NR)
      REAL(8)   ,INTENT(OUT):: PSPOT(NR)
      REAL(8)   ,INTENT(OUT):: EOFLN(LNX)   ! ENERGY OF PARTIAL WAVE
      REAL(8)   ,INTENT(OUT):: ESCATT(LNX)  ! ENERGY OF SCATTERING PARTIAL WAVE
      REAL(8)   ,INTENT(OUT):: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: NLPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PRO(NR,LNX)
      REAL(8)   ,INTENT(OUT):: DT(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: DOVER(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: NLPHIDOT(NR,LNX)  
      REAL(8)   ,INTENT(OUT):: PSPHIDOT(NR,LNX)  
      REAL(8)   ,INTENT(OUT):: AEPHIDOT(NR,LNX)  
      INTEGER(4),ALLOCATABLE:: NPROL(:)
      INTEGER(4),ALLOCATABLE:: NCL(:)
      REAL(8)               :: DH(LNX,LNX)
      REAL(8)               :: PSPHIPROBARE(LNX,LNX)
      REAL(8)               :: TRANSPHI(LNX,LNX)
      REAL(8)               :: TRANSPHIINV(LNX,LNX)
      REAL(8)               :: TRANSPRO(LNX,LNX)
      REAL(8)               :: TRANSU(LNX,LNX)
      REAL(8)               :: TRANSUINV(LNX,LNX)
      REAL(8)               :: EOFI1(NB)
      REAL(8)               :: EOFICOMP(2,NB-NC)
      REAL(8)               :: UOFI(NR,NB)
      REAL(8)               :: TUOFI(NR,NB)
      REAL(8)               :: AEPSI(NR,NB)
      REAL(8)               :: TNLPHI(NR,LNX)
      REAL(8)               :: TAEPHI(NR,LNX)
      REAL(8)               :: TPSPHI(NR,LNX)
      REAL(8)               :: QN(NR,LNX)
      REAL(8)               :: TQN(NR,LNX)
      REAL(8)               :: QNP(NR,LNX)
      REAL(8)               :: PSPHIP(NR,LNX)
      REAL(8)               :: BAREPRO(NR,LNX)
      REAL(8)               :: PHITEST1(NR,LNX)
      REAL(8)   ,ALLOCATABLE:: PHITEST(:,:)
      REAL(8)   ,ALLOCATABLE:: TPHITEST(:,:)
      REAL(8)   ,ALLOCATABLE:: PRO1(:,:)
      REAL(8)   ,ALLOCATABLE:: DH1(:,:)
      REAL(8)   ,ALLOCATABLE:: DO1(:,:)
      REAL(8)               :: AERHO(NR),PSRHO(NR)
      REAL(8)               :: G(NR),DREL(NR),G1(NR),PHI(NR),PHI1(NR)
      REAL(8)               :: E
      REAL(8)               :: RC1
      REAL(8)               :: EHOMO
      INTEGER(4)            :: LX
      INTEGER(4)            :: L,IB,LN,IR,IB1,IB2,LN1,LN2,I,ISO
      INTEGER(4)            :: NV,NPRO,IV,IPRO,IPRO1,IPRO2
      REAL(8)               :: PI,Y0
      REAL(8)               :: X0,Z0,DX,XM,ZM
      INTEGER(4)            :: ISTART,IBI
      INTEGER(4)            :: NITER=100
      REAL(8)               :: PHIPHASE
      REAL(8)               :: R(NR)
      REAL(8)   ,PARAMETER  :: TOL=1.D-7
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)   ,ALLOCATABLE:: AUXARR(:,:)
      REAL(8)               :: VAL,DER,JVAL,JDER,KVAL,KDER
      REAL(8)               :: SVAR,SVAR1,SVAR2
      LOGICAL(4)            :: TREL,TSO,TCHK
      REAL(8)   ,ALLOCATABLE:: A(:,:),AINV(:,:)
      REAL(8)   ,ALLOCATABLE:: PROJ(:)
      REAL(8)               :: AEPSIF(NR,NB-NC)
      REAL(8)               :: PSPSIF(NR,NB-NC)
      REAL(8)               :: EH,EXC
      INTEGER(4)            :: NN,NN0
      INTEGER(4)            :: NFIL
      CHARACTER(64)         :: STRING
      REAL(8)               :: RCOV    !COVALENT RADIUS
      REAL(8)               :: RASA    !COVALENT RADIUS*APPROX 1.15
      REAL(8)               :: RNORM   !NORMALIZATIONS ARE DONE WITHIN RNORM
      REAL(8)               :: RBND   !RADIUS FOR BOUNDARY CONDITIONS
      REAL(8)               :: DFACP,DFACM    !DOUBLE FACTORIAL OF 2L+1; 2L-1
      LOGICAL   ,PARAMETER  :: TTEST=.TRUE.
      LOGICAL   ,PARAMETER  :: TWRITE=.TRUE.
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      LX=MAX(MAXVAL(LOX),MAXVAL(LOFI))
      CALL RADIAL$R(GID,NR,R)
      CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RCOV)
      CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RASA)
      RNORM=RBOX
      RBND=RNORM
      PHIPHASE=1.D0    !NODE AT R=RBND!
!
!     ==========================================================================
!     == RESOLVE KEY                                                          ==
!     ==========================================================================
      TREL=INDEX(KEY,'NONREL').EQ.0
      IF(TREL.AND.INDEX(KEY,'REL').EQ.0) THEN
        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
      END IF
      TSO=INDEX(KEY,'NONSO').EQ.0
      IF(TSO.AND.INDEX(KEY,'SO').EQ.0) THEN
        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
      END IF
!
!     == DETERMINE HIGHEST CORE STATE FOR EACH ANGULAR MOMENTUM ================
      ALLOCATE(NCL(0:LX))
      NCL(:)=0
      DO IB=1,NC
        L=LOFI(IB)
        NCL(L)=MAX(NCL(L),IB)
      ENDDO
!
!     == DETERMINE NUMBER OF PROJECTORS FOR EACH ANGULAR MOMENTUM ==============
      ALLOCATE(NPROL(0:LX))
      NPROL(:)=0
      DO LN=1,LNX
        L=LOX(LN)
        NPROL(L)=NPROL(L)+1
      ENDDO
!
!     ==========================================================================
!     == CONSTRUCT PSEUDO POTENTIAL                                           ==
!     ==========================================================================
      CALL ATOMIC_PSEUDIZE(GID,NR,POW_POT,tval0_pot,VAL0_POT,RC_POT,AEPOT,PSPOT)
!
!     ==========================================================================
!     == CONSTRUCT NODELESS WAVE FUNCTIONS                                    ==
!     ==========================================================================
      DO L=0,LX
        G(:)=0.D0
        DO IB=1,NB
          IF(LOFI(IB).NE.L) CYCLE
          E=EOFI(IB)
          DREL(:)=0.D0
          IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,ROUT,DREL,G,0,AEPOT,E,UOFI(:,IB))
          EOFI1(IB)=E
          TUOFI(:,IB)=G+(E-AEPOT(:)*Y0)*UOFI(:,IB)
          G(:)=UOFI(:,IB)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REPORT SETTINGS ON WAVE FUNCTIONS                                    ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20," Z=",F5.0," ")')AEZ
        WRITE(6,FMT='(82("="),T20," ENERGIES FOR ATOMIC WAVE FUNCTIONS ")')
        WRITE(6,FMT='(82("="),T20," OLD: AE SCHRODINGER EQUATION           ")')
        WRITE(6,FMT='(82("="),T20," NEW: NODELESS EQUATION                 ")')
        WRITE(6,FMT='(82("="),T20," DIFFERENCE DUE TO RELATIVISTIC EFFECTS ")')
        DO IB=1,NB
          WRITE(6,FMT='("IB=",I3," L=",I2," F=",F10.5," E[NEW]=",F15.5 &
     &                                               ," E[OLD]=",F15.5)') &
     &                  IB,LOFI(IB),FOFI(IB),EOFI1(IB),EOFI(IB)
        ENDDO
!       CALL SETUP_WRITEPHI('UOFI.DAT',GID,NR,NB,UOFI)
      END IF
!
!     ==========================================================================
!     == CONSTRUCT NODELESS PARTIAL WAVES                                     ==
!     ==========================================================================
      DFACP=1.D0
      DO L=0,LX
        ISO=0
        DFACM=DFACP
        DFACP=DFACP*REAL(2*L+1,KIND=8)
        E=0.D0
        DO IB=NC+1,NB
          IF(LOFI(IB).NE.L) CYCLE
          E=EOFI(IB)
          EXIT
        ENDDO
        G(:)=0.D0
        IF(NCL(L).NE.0)G(:)=UOFI(:,NCL(L))
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
!         == SEARCH NODELESS STATE =============================================
          ISTART=1
          X0=E
          DX=1.D-2
          CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
          DO I=1,NITER
            E=X0
            IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
            CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,ISO,G,L,E,1,PHI)
            CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RBND,Z0)
            Z0=Z0-PHIPHASE
            IF(ABS(2.D0*DX).LE.TOL) EXIT
            CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
          ENDDO
          IF(ABS(DX).GT.TOL) THEN
            CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
            CALL ERROR$MSG('BOUND STATE NOT FOUND')
            CALL ERROR$STOP('ATOMIC_MAKEPARTIALWAVES')
          END IF
          EOFLN(LN)=E
          NLPHI(:,LN)=PHI(:)
          TNLPHI(:,LN)=G(:)+(E-AEPOT(:)*Y0)*PHI(:)
          G(:)=NLPHI(:,LN)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REPORT SETTINGS ON PARTIAL WAVES                                     ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20," ENERGIES FOR PARTIAL-WAVE CONSTRUCTION")')
        WRITE(6,FMT='("RBOX=",F9.5)')RBOX
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," E=",F10.5," RC=",F6.3)') &
     &                      LN,LOX(LN),EOFLN(LN),RC(LN)
        ENDDO
      END IF
!
!     ==========================================================================
!     == NORMALIZE EACH ANGULAR MOMENTUM SO THAT FIRST PARTIAL WAVE IS NORMAL ==
!     ==========================================================================
      DO L=0,LX
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          PHI(:)=NLPHI(:,LN)
!         == ORTHOGONALIZE TO CORE STATES  =====================================
          DO IB=NC,1,-1
            IF(LOFI(IB).NE.L) CYCLE
            AUX(:)=R(:)**2*UOFI(:,IB)*PHI(:)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR1)
            AUX(:)=R(:)**2*UOFI(:,IB)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR2)
            VAL=SVAR1/SVAR2
            PHI(:)=PHI(:)-UOFI(:,IB)*VAL
          ENDDO
!         == NORMALIZATION FACTOR  =============================================
          AUX(:)=R(:)**2*PHI(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RNORM,VAL)
          VAL=1.D0/SQRT(VAL)
          CALL RADIAL$VALUE(GID,NR,PHI,MAXVAL(RC),SVAR)
          VAL=SIGN(VAL,SVAR)
!         == SCALE ALL WAVE FUNCTIONS ==========================================
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            UOFI(:,IB) = UOFI(:,IB)*VAL
            TUOFI(:,IB)=TUOFI(:,IB)*VAL
          ENDDO
!         == SCALE ALL PARTIAL WAVES  ==========================================
          DO LN1=1,LNX
            IF(LOX(LN1).NE.L) CYCLE
            NLPHI(:,LN1)   =NLPHI(:,LN1)*VAL
            TNLPHI(:,LN1)  =TNLPHI(:,LN1)*VAL
          ENDDO
          EXIT ! SCALE ONLY ONCE PER L
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CONSTRUCT QN FUNCTIONS        (H-E)|QN>=|UC>                         ==
!     ==========================================================================
      TRANSU(:,:)=0.D0
      DO L=0,LX
        IPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          IPRO=IPRO+1
          SVAR=1.D0
          DO LN1=1,LN
            IF(LOX(LN1).NE.L) CYCLE
            TRANSU(LN1,LN)=TRANSU(LN1,LN)+SVAR
            SVAR=SVAR*(EOFLN(LN)-EOFLN(LN1))
          ENDDO
        ENDDO
      ENDDO
      CALL LIB$INVERTR8(LNX,TRANSU,TRANSUINV)
!
      QN=MATMUL(NLPHI,TRANSU)
      TQN=MATMUL(TNLPHI,TRANSU)

      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  TRANSU ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F15.10)')TRANSU(LN1,:)
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  TRANSU^(-1) ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F15.10)')TRANSUINV(LN1,:)
        ENDDO
      END IF

CALL SETUP_WRITEPHI('UN.DAT',GID,NR,LNX,NLPHI)
CALL SETUP_WRITEPHI('QN.DAT',GID,NR,LNX,QN)
!
!     ==========================================================================
!     == TEST EQUATION FOR QN                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  TEST QN EQ.  ")')
        DO L=0,LX
          IPRO=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IPRO=IPRO+1
            PRO(:,LN)=TQN(:,LN)+(AEPOT*Y0-EOFLN(LN))*QN(:,LN)
            IF(NCL(L).NE.0) PRO(:,LN)=PRO(:,LN)-UOFI(:,NCL(L))
            WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|UC>=",F20.15)') &
     &                     LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == CORE-ORTHOGONALIZE QN TO OBTAIN NODAL AE PARTIAL WAVES               ==
!     == USE LADDER OF NODELESS WAVE FUNCTIONS                                ==
!     == DUE TO NEGLECT OF THE SMALL COMPONENT THE ORTHOGONALIZATION IS NOT   ==
!     == EXACT. HOWEVER THIS CHOICE ENSURES THAT AEPHI FULFILLS THE           ==
!     == SCHRODINGER EQUATION                                                 ==
!     ==========================================================================
      AEPHI(:,:) =QN(:,:)   
      TAEPHI(:,:)=TQN(:,:)  
      DO LN=1,LNX
        SVAR=1.D0
        DO IB=NC,1,-1
          IF(LOFI(IB).NE.LOX(LN)) CYCLE
          SVAR=SVAR/(EOFLN(LN)-EOFI1(IB))
          AEPHI(:,LN) = AEPHI(:,LN)+ UOFI(:,IB)*SVAR
          TAEPHI(:,LN)=TAEPHI(:,LN)+TUOFI(:,IB)*SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == TEST ORTHOGONALITY OF AEPHI TO CORE STATES                           ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(A(LNX,NC))
        A(:,:)=0.D0
        DO LN=1,LNX
          DO IB=1,NC
            IF(LOX(LN).NE.LOFI(IB)) CYCLE
            AUX(:)=R(:)**2*UOFI(:,IB)*AEPHI(:,LN)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
            A(LN,IB)=VAL
            AUX(:)=R(:)**2*UOFI(:,IB)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
            A(LN,IB)=A(LN,IB)/SQRT(VAL)
          ENDDO
        ENDDO
        WRITE(6,FMT='(82("="),T20," <UC|AEPHI>/SQRT(<UC|UC>  ")')
        WRITE(6,FMT='(82("="),T10," DEVIATION DUE NEGLECT OF SMALL COMPONENT")')
        WRITE(6,FMT='(82("="),T10," ORTHOGONALIZATION DONE BASED ON ENERGIES")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F10.5)')A(LN1,:)
        ENDDO
        DEALLOCATE(A)
      END IF
!
!     ==========================================================================
!     == TEST EQUATION FOR ALL-ELECTRON PARTIAL WAVES AEPHI                   ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  TEST AEPHI EQ.  ")')
        DO LN=1,LNX
          PRO(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
          WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|AEPHI_N> =",F20.15)') &
      &                 LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
        ENDDO
      END IF
!
!     ==========================================================================
!     == CONSTRUCT PSEUDO PARTIAL WAVES                                       ==
!     ==========================================================================
      PSPHI=QN
      TPSPHI=TQN
      IF(TYPE.EQ.'KERKER') THEN
        DO L=0,LX
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            RC1=RC(LN)
          ENDDO
          NPRO=NPROL(L)
          IF(NPRO.EQ.0) CYCLE
          ALLOCATE(PHITEST(NR,NPRO))
          ALLOCATE(TPHITEST(NR,NPRO))
          IPRO=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IPRO=IPRO+1
            PHITEST(:,IPRO)=QN(:,LN)
            TPHITEST(:,IPRO)=TQN(:,LN)
          ENDDO
          CALL ATOMIC_MAKEPSPHI(GID,NR,RC1,L,NPRO,PHITEST,TPHITEST)
          IPRO=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IPRO=IPRO+1
            PSPHI(:,LN)=PHITEST(:,IPRO)
            TPSPHI(:,LN)=TPHITEST(:,IPRO)
          ENDDO
          DEALLOCATE(PHITEST)
          DEALLOCATE(TPHITEST)
        ENDDO
      ELSE IF(TYPE.EQ.'HBS') THEN
        CALL ATOMIC_MAKEPSPHI_HBS(GID,NR,LNX,LOX,EOFLN,RC,LAMBDA,PSPOT &
     &                               ,rbnd,PSPHI,TPSPHI)
      ELSE
        CALL ERROR$MSG('PSEUDIZATION TYPE IS UNKNOWN')
        CALL ERROR$MSG('CAN BE "BESSEL" OR "HBS"')
        CALL ERROR$STOP('ATOMIC_MAKEPARTIALWAVES')
      END IF
!
!     ==========================================================================
!     == CONSTRUCT PROJECTOR FUNCTIONS                                        ==
!     ==========================================================================
      DO LN=1,LNX
        BAREPRO(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
      ENDDO
!     ==  CLEANUP OF NUMERICAL ERRORS INDUCED BY A TOO COARSE GRID AT LARGE ====
!     ==  RADII. 2.D0*RCOV IS CHOSEN ARBITRARILY ===============================
      DO IR=1,NR
        IF(R(IR).LT.RNORM) CYCLE
        BAREPRO(IR:,:)=0.D0
        EXIT
      ENDDO

      IF(TTEST) THEN
        CALL SETUP_WRITEPHI('PRO-BARE.DAT',GID,NR,LNX,BAREPRO)
      END IF
!
!     ==========================================================================
!     == CHECK PAW EQUATION FOR PSEUDO PARTIALWAVES                           ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(PHITEST(NR,LNX))
        ALLOCATE(TPHITEST(NR,LNX))
        DO LN=1,LNX
          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN) &
     &                  -BAREPRO(:,LN)
          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
        ENDDO
!
        WRITE(6,FMT='(82("="),T20,"  TEST RAW PAW EQUATION  ")')
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," RAW PAW EQ.",F10.5 &
     &                                 ," SCHR. EQ.",F10.5)') &
     &         LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))),MAXVAL(ABS(PHITEST(:,LN)))
        ENDDO
        DEALLOCATE(PHITEST)
        DEALLOCATE(TPHITEST)
      ENDIF
!
!     ==========================================================================
!     == ENFORCE BIORTHOGONALIZATION                                          ==
!     ==========================================================================
      ALLOCATE(A(LNX,LNX))
      CALL BIORTHOMATRICES(GID,NR,RBOX,LNX,LOX,PSPHI,BAREPRO,TRANSPHI,TRANSPRO)
      A=MATMUL(TRANSPRO,TRANSPOSE(TRANSPHI))
      PRO=MATMUL(BAREPRO,A)
      CALL LIB$INVERTR8(LNX,TRANSPHI,TRANSPHIINV)
      DEALLOCATE(A)
!CALL SETUP_WRITEPHI(-'PRO.DAT',GID,NR,LNX,MATMUL(BAREPRO,TRANSPRO))
!CALL SETUP_WRITEPHI(-'PSPHI.DAT',GID,NR,LNX,MATMUL(PSPHI,TRANSPHI))
!CALL SETUP_WRITEPHI(-'AEPHI.DAT',GID,NR,LNX,MATMUL(AEPHI,TRANSPHI))
!
!     ==========================================================================
!     == CHECK BIORTHOGONALIZATION                                            ==
!     ==========================================================================
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=R(:)**2*PSPHI(:,LN1)*PRO(:,LN2)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          IF(LN1.EQ.LN2)VAL=VAL-1.D0
          IF(ABS(VAL).GT.1.D-5) THEN
            CALL ERROR$MSG('BIORTHOGONALIZATION FAILED')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('LN1',LN1)
            CALL ERROR$I4VAL('LN2',LN2)
            CALL ERROR$R8VAL('DEVIATION',VAL)
            CALL ERROR$STOP('ATOMIC_MAKEPARTIALWAVES')
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DT,DO                                                                ==
!     == ATTENTION!! DO NOT MAKE DT AND DH SYMMETRIC. THEY ARE NOT HERMITEAN! ==
!     ==========================================================================
      DT=0.D0
      DOVER=0.D0
      DH=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*TAEPHI(:,LN2)-PSPHI(:,LN1)*TPSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          DT(LN1,LN2)=VAL
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          DOVER(LN1,LN2)=VAL
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPOT(:)*Y0*AEPHI(:,LN2) &
      &                  -PSPHI(:,LN1)*PSPOT(:)*Y0*PSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          DH(LN1,LN2)=DT(LN1,LN2)+VAL
        ENDDO
      ENDDO
      IF(TTEST) THEN
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," DTKIN-TRANSPOSE(DTKIN)=",10F10.5)') &
     &                                    LN,(DT(LN,LN2)-DT(LN2,LN),LN2=1,LNX)
        ENDDO
      END IF
!
!     ==========================================================================
!     == CHECK PAW EQUATION FOR PSEUDO PARTIALWAVES                           ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(TPHITEST(NR,LNX))    ! HOLDS TEST FOR PSEUDO
        ALLOCATE(PHITEST(NR,LNX))     ! HOLDS TEST FOR ALL-ELECTRON 
        ALLOCATE(PROJ(LNX))
        DO LN=1,LNX
!         == DETERMINE PROJECTIONS
          PROJ(:)=0.D0
          DO LN1=1,LNX
            IF(LOX(LN1).NE.LOX(LN)) CYCLE
            AUX(:)=R(:)**2*PRO(:,LN1)*PSPHI(:,LN)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
            PROJ(LN1)=VAL
          ENDDO
          WRITE(6,FMT='("LN=",I2," <P|PSPHI>=",10F10.5)')LN,PROJ
!
          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
!TPHITEST(:,LN)=TPHITEST(:,LN)-BAREPRO(:,LN)
          PHITEST1(:,LN)=0.D0          
          DO LN1=1,LNX
            IF(LOX(LN1).NE.LOX(LN)) CYCLE
            SVAR=0.D0
            DO LN2=1,LNX
              IF(LOX(LN2).NE.LOX(LN)) CYCLE
              SVAR=SVAR+(DH(LN1,LN2)-EOFLN(LN2)*DOVER(LN1,LN2))*PROJ(LN2)
            ENDDO
            TPHITEST(:,LN)=TPHITEST(:,LN)+PRO(:,LN1)*SVAR
            PHITEST1(:,LN)=PHITEST1(:,LN)+PRO(:,LN1)*SVAR
          ENDDO
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  TEST PAW EQUATION  ")')
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," PAW EQ.",F20.5 &
     &                                 ," SCHR. EQ.",F20.5," DPRO ",F20.5)') &
     &          LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))) &  
     &                    ,MAXVAL(ABS(PHITEST(:,LN))) &   ! ALL-ELECTRON EQ
     &                    ,MAXVAL(ABS(PHITEST1(:,LN)+BAREPRO(:,LN)))
        ENDDO
!        CALL SETUP_WRITEPHI('PRO-TEST.DAT',GID,NR,LNX,-BAREPRO)
        DEALLOCATE(PROJ)
        DEALLOCATE(PHITEST)
        DEALLOCATE(TPHITEST)
      END IF
!
!     ==========================================================================
!     == CONSTRUCT PHIDOT FUNCTIONS                                           ==
!     ==========================================================================
!     == TRANSFORMATION TO NODELESS REPRESENTATION =============================
      QNP(:,:)   =MATMUL(QN,TRANSUINV)
      PSPHIP(:,:)=MATMUL(PSPHI,TRANSUINV)
      DO L=0,LX
        NPRO=NPROL(L)
        IF(NPRO.EQ.0) CYCLE
        ALLOCATE(DH1(NPRO,NPRO))
        ALLOCATE(DO1(NPRO,NPRO))
        ALLOCATE(PRO1(NR,NPRO))
        ALLOCATE(PROJ(NPRO))
        IPRO1=0
        DO LN1=1,LNX
          IF(LOX(LN1).NE.L) CYCLE
          IPRO1=IPRO1+1
          PRO1(:,IPRO1)=PRO(:,LN1)
          IPRO2=0
          DO LN2=1,LNX
            IF(LOX(LN2).NE.L) CYCLE
            IPRO2=IPRO2+1
            DH1(IPRO1,IPRO2)=DH(LN1,LN2)
            DO1(IPRO1,IPRO2)=DOVER(LN1,LN2)
          ENDDO
        ENDDO
        G(:)=0.D0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          ESCATT(LN)=MIN(-0.0D0,EOFLN(LN))
          E=ESCATT(LN)
          IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
!
!         == CALCULATE QN  FUNCTION ===========================================
          G(:)=0.D0
          IF(NCL(L).GT.0)G(:)=UOFI(:,NCL(L))
          CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,E,1,PHI)
          IF(NCL(L).EQ.0) THEN
            CALL RADIAL$VALUE(GID,NR,PHI,1.D-2,SVAR1)
            CALL RADIAL$VALUE(GID,NR,QN(:,LN),1.D-2,SVAR2)
            PHI(:)=PHI(:)*SVAR2/SVAR1
          END IF
          NLPHIDOT(:,LN)=PHI(:)
!
!         == CALCULATE PSEUDO WAVE FUNCTIONS ===================================
          G(:)=0.D0
          CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO1,DH1,DO1,G,PHI)
          CALL RADIAL$VALUE(GID,NR,PHI,MAXVAL(RC),SVAR)
          CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LN),MAXVAL(RC),VAL)
          PSPHIDOT(:,LN)=PHI(:)*VAL/SVAR
!
          IF(ABS(E-EOFLN(LN)).GT.1.D-3) THEN
            DO LN1=1,LNX
              IF(LOX(LN1).NE.L) CYCLE
              IF(ISCATT(LN1).GT.ISCATT(LN)) CYCLE
              NLPHIDOT(:,LN)   =(   NLPHIDOT(:,LN)-   QNP(:,LN1))/(E-EOFLN(LN1))
              PSPHIDOT(:,LN)=(PSPHIDOT(:,LN)-PSPHIP(:,LN1))/(E-EOFLN(LN1))     
            ENDDO
          ELSE
!           == ENERGY DERIVATIVE OF QN ========================================
            G(:)=NLPHIDOT(:,LN)
            CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,E,1,PHI)
            NLPHIDOT(:,LN)=PHI
!           == ENERGY DERIVATIVE OF PSPHI ====================================
            PHI=PSPHIDOT(:,LN)    
            DO IPRO=1,NPRO
              AUX(:)=R(:)**2*PRO1(:,IPRO)*PHI(:)
              CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
              CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,PROJ(IPRO))
            ENDDO
            G(:)=PHI(:)+MATMUL(PRO1,MATMUL(DO1,PROJ))
            CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO1,DH1,DO1,G,PHI1)
            CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LN)-PHI1,MAXVAL(RC),VAL)
            CALL RADIAL$VALUE(GID,NR,PHI,MAXVAL(RC),SVAR)
            PSPHIDOT(:,LN)=PHI1(:)+PHI*VAL/SVAR
            DO LN1=1,LNX
              IF(LOX(LN1).NE.L) CYCLE
              IF(ISCATT(LN1).GE.ISCATT(LN)) CYCLE
              NLPHIDOT(:,LN)   =(   NLPHIDOT(:,LN)-   QNP(:,LN1))/(E-EOFLN(LN1))
              PSPHIDOT(:,LN)=(PSPHIDOT(:,LN)-PSPHIP(:,LN1))/(E-EOFLN(LN1))          
            ENDDO
          END IF
!
!         == CONSTRUCT AEPHIDOT BY CORE-ORTHOGONALIZATION ======================
          AEPHIDOT(:,LN)=NLPHIDOT(:,LN)
          DO IB=NC,1,-1
            IF(LOFI(IB).NE.L) CYCLE
            AUX=R(:)**2*UOFI(:,IB)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR1)
            AUX=R(:)**2*UOFI(:,IB)*AEPHIDOT(:,LN)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR2)
            AEPHIDOT(:,LN)=AEPHIDOT(:,LN)-UOFI(:,IB)*SVAR2/SVAR1
          ENDDO
        ENDDO
        DEALLOCATE(DH1)
        DEALLOCATE(DO1)
        DEALLOCATE(PRO1)
        DEALLOCATE(PROJ)
      ENDDO
!
!     ==========================================================================
!     == BACK TRANSFORM                                                       ==
!     ==========================================================================
GOTO 100
      QN=MATMUL(QN,TRANSUINV)
      TQN=MATMUL(TQN,TRANSUINV)
      PSPHI=MATMUL(PSPHI,TRANSUINV)
      TPSPHI=MATMUL(TPSPHI,TRANSUINV)
      AEPHI=MATMUL(AEPHI,TRANSUINV)
      TAEPHI=MATMUL(TAEPHI,TRANSUINV)
      PRO=MATMUL(PRO,TRANSPOSE(TRANSU))
      DT=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DT,TRANSUINV))
      DOVER=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DOVER,TRANSUINV))
      DH=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DH,TRANSUINV))
!
!     == TEST IF BACK TRANSFORM WAS SUCCESSFUL ================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  TEST BACK TRANSFORM  ")')
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," DIFF. NDLSS PHI",F10.5 &
     &                                 ," DIFF. KIN.OP NDLSS. PHI ",F10.5)') &
     &          LN,LOX(LN),MAXVAL(ABS(QN(:,LN)-NLPHI(:,LN))) &
     &                    ,MAXVAL(ABS(TQN(:,LN)-TNLPHI(:,LN)))
        ENDDO
      END IF
100 CONTINUE
!
!     ==========================================================================
!     == RENORMALIZE WAVE FUNCTIONS AND PROJECTOR FUNCTIONS                   ==
!     ==========================================================================
GOTO 10001
      DO L=0,LX
        IPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          IPRO=IPRO+1
!         == NORMALIZE PS PARTIAL WAVE =========================================
!          IF(IPRO.EQ.1) THEN
          IF(ISCATT(LN).LE.0) THEN    ! NORMALIZE VALENCE AND SEMI-CORE STATES
            AUX(:)=R(:)**2*PSPHI(:,LN)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RNORM,VAL)
            VAL=VAL+DOVER(LN,LN)
            VAL=1.D0/SQRT(VAL)
          END IF
          DO LN2=LN,LNX
            IF(LOX(LN2).NE.L) CYCLE
            PSPHI(:,LN2)   =   PSPHI(:,LN2)*VAL
            TPSPHI(:,LN2)  =  TPSPHI(:,LN2)*VAL
            AEPHI(:,LN2)   =   AEPHI(:,LN2)*VAL
            TAEPHI(:,LN2)  =  TAEPHI(:,LN2)*VAL
            NLPHI(:,LN2)   =   NLPHI(:,LN2)*VAL
            TNLPHI(:,LN2)  =  TNLPHI(:,LN2)*VAL
            PSPHIDOT(:,LN2)=PSPHIDOT(:,LN2)*VAL
            AEPHIDOT(:,LN2)=AEPHIDOT(:,LN2)*VAL
            NLPHIDOT(:,LN2)=NLPHIDOT(:,LN2)*VAL
            QN(:,LN2)      =      QN(:,LN2)*VAL
            TQN(:,LN2)     =     TQN(:,LN2)*VAL
            PRO(:,LN2)     =     PRO(:,LN2)/VAL
            DH(LN2,:)      =      DH(LN2,:)*VAL
            DT(LN2,:)      =      DT(LN2,:)*VAL
            DOVER(LN2,:)   =   DOVER(LN2,:)*VAL
            DH(:,LN2)      =      DH(:,LN2)*VAL
            DT(:,LN2)      =      DT(:,LN2)*VAL
            DOVER(:,LN2)   =   DOVER(:,LN2)*VAL
          ENDDO
        ENDDO
      ENDDO
10001 CONTINUE
!
!     ==========================================================================
!     == CALCULATE DENSITY FOR UNSCREENING                                    ==
!     ==========================================================================
      AERHO(:)=0.D0
      PSRHO(:)=0.D0
      EOFICOMP(:,:)=0.D0
      DO L=0,LX
        NPRO=NPROL(L)
        IF(NPRO.EQ.0) CYCLE
        ALLOCATE(DH1(NPRO,NPRO))
        ALLOCATE(DO1(NPRO,NPRO))
        ALLOCATE(PRO1(NR,NPRO))
        ALLOCATE(PROJ(NPRO))
        IPRO1=0
        DO LN1=1,LNX
          IF(LOX(LN1).NE.L) CYCLE
          IPRO1=IPRO1+1
          PRO1(:,IPRO1)=PRO(:,LN1)
          IPRO2=0
          DO LN2=1,LNX
            IF(LOX(LN2).NE.L) CYCLE
            IPRO2=IPRO2+1
            DH1(IPRO1,IPRO2)=DH(LN1,LN2)
            DO1(IPRO1,IPRO2)=DOVER(LN1,LN2)
          ENDDO
        ENDDO
        NN0=-1
        G(:)=0.D0
        DO IB=NC+1,NB
          IF(LOFI(IB).NE.L) CYCLE
          IF(NN0.EQ.-1)NN0=NNOFI(IB)
           E=EOFI1(IB)
!
          G(:)=0.D0
          IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,ROUT,DREL,G,NNOFI(IB),AEPOT &
       &                             ,E,AEPSIF(:,IB-NC))
          SVAR1=E
          EOFICOMP(1,IB-NC)=E
!
          NN=NNOFI(IB)-NN0
          G(:)=0.D0
          CALL ATOMLIB$PAWBOUNDSTATE(GID,NR,L,NN,ROUT,PSPOT,NPRO,PRO1,DH1,DO1 &
     &                              ,G,E,PSPSIF(:,IB-NC))
          SVAR2=E
          EOFICOMP(2,IB-NC)=E
          IF(ABS(SVAR2-EOFI1(IB)).GT.1.D-2) THEN
            CALL ERROR$MSG('INACCURACY WHILE UNSCREENING PS POTENTIAL')
            CALL ERROR$MSG('ONE-PARTICLE ENERGIES OBTAINED FROM PAW ')
            CALL ERROR$MSG('DISAGREE WITH THOSE FROM THE AE CALCULATION')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('IB',IB)
            CALL ERROR$R8VAL('TARGET: E[EV]',EOFI1(IB)*27.211D0)
            CALL ERROR$R8VAL('TARGET: E[EV]',EOFI(IB)*27.211D0)
            CALL ERROR$R8VAL('AE:     E[EV]',SVAR1*27.211D0)
            CALL ERROR$R8VAL('PAW:    E[EV]',SVAR2*27.211D0)
            CALL ERROR$R8VAL('( PAW-AE)[EV]',(SVAR2-SVAR1)*27.211D0)
            CALL ERROR$R8VAL('(PAW-REF)[EV]',(SVAR2-EOFI1(IB))*27.211D0)
            CALL ERROR$R8VAL('( AE-REF)[EV]',(SVAR1-EOFI1(IB))*27.211D0)
            CALL ERROR$STOP('ATOMLIB_MAKEPARTIALWAVES')
          END IF
          IF(TTEST) THEN
             WRITE(6,FMT='("DEVIATION OF THE ATOMIC ENERGY LEVELS IN EV")')
             WRITE(6,FMT='("OBTAINED ONCE WITH PAW AND THE AE CALCULATION")')
             WRITE(6,FMT='("L",I2," PAW-REF ",F10.5,"EV; AE-REF ",F10.5," EV")') &
         &          L,(SVAR2-EOFI1(IB))*27.211D0,(SVAR1-EOFI1(IB))*27.211D0
          END IF
!
          DO IR=1,NR-2
            IF(R(IR).LT.ROUT) CYCLE
            PSPSIF(IR+2:,IB-NC)=0.D0
            AEPSIF(IR+2:,IB-NC)=0.D0
            IF(PSPSIF(IR-2,IB-NC)*AEPSIF(IR-2,IB-NC).LT.0.D0) &
     &                                          AEPSIF(:,IB-NC)=-AEPSIF(:,IB-NC)
            EXIT
          ENDDO
!
!         == CALCULATE PROJECTIONS PROJ=<P|PS-PSI>  ============================
          DO IPRO=1,NPRO
            AUX(:)=R(:)**2*PSPSIF(:,IB-NC)*PRO1(:,IPRO)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,PROJ(IPRO))
          ENDDO
!
!         == NORMALIZE PS WAVE FUNCTION=========================================
          AUX(:)=R(:)**2*PSPSIF(:,IB-NC)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
          DO IPRO1=1,NPRO
            DO IPRO2=1,NPRO
              VAL=VAL+PROJ(IPRO1)*DO1(IPRO1,IPRO2)*PROJ(IPRO2)
            ENDDO
          ENDDO
          PSPSIF(:,IB-NC)=PSPSIF(:,IB-NC)/SQRT(VAL)
          PROJ=PROJ/SQRT(VAL)
!
    PRINT*,'PROJ ',L,PROJ
!
!         == NORMALIZE AE WAVE FUNCTION=========================================
          AUX(:)=R(:)**2*AEPSIF(:,IB-NC)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
          AEPSIF(:,IB-NC)=AEPSIF(:,IB-NC)/SQRT(VAL)
!
          PSRHO(:)=PSRHO(:)+FOFI(IB)*PSPSIF(:,IB-NC)**2*Y0
          AERHO(:)=AERHO(:)+FOFI(IB)*AEPSIF(:,IB-NC)**2*Y0
        ENDDO
        DEALLOCATE(DH1)
        DEALLOCATE(DO1)
        DEALLOCATE(PRO1)
        DEALLOCATE(PROJ)
      ENDDO      

!     == REPORT ===============================================================
      WRITE(6,FMT='(82("="),T20,"  AE- AND PS-ENERGY LEVELS ")')
      DO IB=NC+1,NB
        WRITE(6,FMT='("IB=",I2," L=",I2," AE-ENERGY:",F10.5," EV;" &
     &                                 ," PS-ENERGY:",F10.5," EV")') &
     &          IB,LOFI(IB),EOFICOMP(:,IB-NC)*27.211D0
      ENDDO
!
!!$      IF(TTEST) THEN
!!$        CALL SETUP_WRITEPHI('PSPSIF.DAT',GID,NR,NB-NC,PSPSIF)
!!$        CALL SETUP_WRITEPHI('AEPSIF.DAT',GID,NR,NB-NC,AEPSIF)
!!$        ALLOCATE(AUXARR(NR,2))
!!$        AUXARR(:,1)=AERHO
!!$        AUXARR(:,2)=PSRHO
!!$        CALL SETUP_WRITEPHI('RHO.DAT',GID,NR,2,AUXARR)
!!$        DEALLOCATE(AUXARR)
!!$      END IF
!
!     ==========================================================================
!     == UNSCREENING                                                          ==
!     ==========================================================================
      SVAR=AEZ
      DO IB=1,NC
        SVAR=SVAR-FOFI(IB)
      ENDDO
      CALL ATOMIC_UNSCREEN(GID,NR,Rout,SVAR,AERHO,PSRHO,PSPOT,RCSM,VADD)
!!$      IF(TTEST) THEN
!!$        ALLOCATE(AUXARR(NR,2))
!!$        AUXARR(:,1)=VADD
!!$        AUXARR(:,2)=PSPOT
!!$        CALL SETUP_WRITEPHI('VADD.DAT',GID,NR,2,AUXARR)
!!$        DEALLOCATE(AUXARR)
!!$      END IF
!
!     ==========================================================================
!     == CUT OFF THE EXPONENTIALLY GROWING TAIL OF THE PARTIALWAVES
!     ==========================================================================
      DO IR=1,NR
        IF(R(IR).GT.MAX(2.D0*RCOV,RNORM)) THEN
          I=IR+1
          EXIT
        END IF
      ENDDO
      IR=I
      AEPHI(IR:,:)=0.D0
      PSPHI(IR:,:)=0.D0
      NLPHI(IR:,:)=0.D0
      aePHIDOT(IR:,:)=0.D0
      NLPHIDOT(IR:,:)=0.D0
      PSPHIDOT(IR:,:)=0.D0
      QN(IR:,:)=0.D0
      PRO(IR:,:)=0.D0
!
!     ==========================================================================
!     == TRANSFORM ONTO SEQUENTIAL RESPRESENTATION                            ==
!     ==========================================================================
!!$PSPHI=MATMUL(PSPHI,TRANSPHI)
!!$AEPHI=MATMUL(AEPHI,TRANSPHI)
!!$NLPHI=MATMUL(NLPHI,TRANSPHI)
!!$PRO=MATMUL(PRO,TRANSPOSE(TRANSPHIINV))
!!$DT=MATMUL(TRANSPOSE(TRANSPHI),MATMUL(DT,TRANSPHI))
!!$DOVER=MATMUL(TRANSPOSE(TRANSPHI),MATMUL(DOVER,TRANSPHI))
!
!     ==========================================================================
!     == WRITE INFORMATION TO FILE                                            ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  DTKIN  ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DT(LN1,:)
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  DOVERLAP  ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DOVER(LN1,:)
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  DHAMILTONIAN ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DH(LN1,:)
        ENDDO
!
        WRITE(6,FMT='(82("="),T20,"  AE-OVERLAP ")')
        ALLOCATE(A(LNX,LNX))
        A(:,:)=0.D0
        DO LN1=1,LNX
          DO LN2=1,LNX
            IF(LOX(LN2).NE.LOX(LN1)) CYCLE
            AUX(:)=R(:)**2*AEPHI(:,LN1)*AEPHI(:,LN2)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RNORM,VAL)
            A(LN1,LN2)=VAL
          ENDDO
          WRITE(6,FMT='(20F12.3)')A(LN1,:)
        ENDDO
        DEALLOCATE(A)
!
        WRITE(6,FMT='(82("="),T20," <PRO|PSPHI> ")')
        ALLOCATE(A(LNX,LNX))
        A(:,:)=0.D0
        DO LN1=1,LNX
          DO LN2=1,LNX
            IF(LOX(LN2).NE.LOX(LN1)) CYCLE
            AUX(:)=R(:)**2*PRO(:,LN1)*PSPHI(:,LN2)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
            A(LN1,LN2)=VAL
          ENDDO
          WRITE(6,FMT='(20F12.3)')A(LN1,:)
        ENDDO
        DEALLOCATE(A)
      END IF
      IF(TWRITE) THEN
        WRITE(STRING,FMT='(F3.0)')AEZ
        STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
!       == AE PARTIAL WAVES
        CALL SETUP_WRITEPHI(-'AEPHI'//TRIM(STRING),GID,NR,LNX,AEPHI)
!
!       == PS PARTIAL WAVES ====================================================
        CALL SETUP_WRITEPHI(-'PSPHI'//TRIM(STRING),GID,NR,LNX,PSPHI)
!
!       == NODELESS PARTIAL WAVES ==============================================
        CALL SETUP_WRITEPHI(-'NLPHI'//TRIM(STRING),GID,NR,LNX,NLPHI)
!
!       == NODELESS PARTIAL WAVES ==============================================
        CALL SETUP_WRITEPHI(-'QPHI'//TRIM(STRING),GID,NR,LNX,QN)
!
!       == PROJECTOR FUNCTIONS =================================================
        CALL SETUP_WRITEPHI(-'PRO'//TRIM(STRING),GID,NR,LNX,PRO)
!
!       == NODELESS SCATTERING PARTIAL WAVES ===================================
        CALL SETUP_WRITEPHI(-'NLPHIDOT'//TRIM(STRING),GID,NR,LNX,NLPHIDOT)
!
!       == SCATTERING PSEUDO PARTIAL WAVES =====================================
        CALL SETUP_WRITEPHI(-'PSPHIDOT'//TRIM(STRING),GID,NR,LNX,PSPHIDOT)
!
!       == SCATTERING ALL-ELECTRON PARTIAL WAVES ===============================
        CALL SETUP_WRITEPHI(-'AEPHIDOT'//TRIM(STRING),GID,NR,LNX,AEPHIDOT)
!
!       == POTENTIALS  =========================================================
        ALLOCATE(AUXARR(NR,4))
        AUXARR(:,1)=AEPOT
        AUXARR(:,2)=PSPOT
        AUXARR(:,3)=PSPOT-VADD
        AUXARR(:,4)=VADD
        CALL SETUP_WRITEPHI(-'POT'//TRIM(STRING),GID,NR,4,AUXARR)
        DEALLOCATE(AUXARR)
!
!       == QBAR ================================================================
        DO LN=1,LNX
          WRITE(*,FMT='("LN= ",I3," L=",I1," ISCATT=",I3," QBAR ",F10.5)') &
     &                 LN,LOX(LN),ISCATT(LN)
        ENDDO
      END IF
!
!STOP 'FORCED: IN MAKEPARTIALWAVES'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_TESTSCATTERING(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: NC
      INTEGER(4)             :: NB
      REAL(8)                :: EMIN
      REAL(8)                :: EMAX
      INTEGER(4),PARAMETER   :: NE=200
      REAL(8)                :: DE
      REAL(8)                :: E
      INTEGER(4)             :: LX
      INTEGER(4)             :: GID ! GRID ID
      INTEGER(4)             :: NR  ! #(GRID POINTS)
      REAL(8)   ,ALLOCATABLE :: R(:)  !(NR) RADIAL GRID
      INTEGER(4)             :: LNX ! #(PARTIAL WAVES)
      INTEGER(4),ALLOCATABLE :: LOX(:) !(LNX) ANGULAR MOMENTA OF PARTIAL WAVES
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      REAL(8)   ,ALLOCATABLE :: G(:)
      REAL(8)   ,ALLOCATABLE :: PHI(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: AEPHASE(:,:)
      REAL(8)   ,ALLOCATABLE :: PAWPHASE(:,:)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:)  !(NR,LNX) PROJECTOR FUNCTIONS
      REAL(8)   ,ALLOCATABLE :: DH(:,:)  !(LNX,LNX) HAMILTON-DIFFERENCE
      REAL(8)   ,ALLOCATABLE :: DO(:,:)  !(LNX,LNX) OVERLAP DIFFERENCE
      REAL(8)                :: VAL
      REAL(8)                :: PI,Y0
      REAL(8)                :: RCOV
      REAL(8)                :: DPHASE
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: IE,LN1,LN2,IPRO1,IPRO2,L,LN,IB
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      NC=THIS%ATOM%NC
      NB=THIS%ATOM%NB
      EMIN=MINVAL(THIS%ATOM%EOFI(NC+1:NB))-0.2D0
      EMAX=MAX(MAXVAL(THIS%ATOM%EOFI(NC+1:NB)),1.D0)
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      ALLOCATE(DREL(NR))
      ALLOCATE(G(NR))
      ALLOCATE(PHI(NR))
      ALLOCATE(AUX(NR))
      LNX=THIS%LNX
      ALLOCATE(LOX(LNX))
      LOX=THIS%LOX
      LX=MAX(MAXVAL(THIS%LOX),MAXVAL(THIS%ATOM%LOFI(:NB)))
      CALL PERIODICTABLE$GET(THIS%AEZ,'R(COV)',RCOV)
!
!     ==========================================================================
!     ==  LOOP OVER ANGULAR MOMENTA                                           ==
!     ==========================================================================
      ALLOCATE(AEPHASE(LX+1,NE))
      ALLOCATE(PAWPHASE(LX+1,NE))
      DO L=0,LX
        NPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          NPRO=NPRO+1
        ENDDO
        DPHASE=0
        DO IB=1,NC
          IF(THIS%ATOM%LOFI(IB).NE.L) CYCLE
          DPHASE=DPHASE+1.D0
        ENDDO
        ALLOCATE(DH(NPRO,NPRO))
        ALLOCATE(DO(NPRO,NPRO))
        ALLOCATE(PRO(NR,NPRO))
        IPRO1=0
        DO LN1=1,LNX
          IF(LOX(LN1).NE.L) CYCLE
          IPRO1=IPRO1+1
          PRO(:,IPRO1)=THIS%PRO(:,LN1)
          IPRO2=0
          DO LN2=1,LNX
            IF(LOX(LN2).NE.L) CYCLE
            IPRO2=IPRO2+1
            AUX(:)=THIS%AEPHI(:,LN1)*THIS%AEPHI(:,LN2)*THIS%ATOM%AEPOT*Y0 &
     &            -THIS%PSPHI(:,LN1)*THIS%PSPHI(:,LN2)*THIS%PSPOT*Y0
            AUX(:)=R(:)**2*AUX(:)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
            DH(IPRO1,IPRO2)=THIS%DTKIN(LN1,LN2)+VAL
            DO(IPRO1,IPRO2)=THIS%DOVER(LN1,LN2)
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  LOOP OVER ENERGIES                                                ==
!       ========================================================================
        DO IE=1,NE
          E=EMIN+DE*REAL(IE-1,KIND=8)
          CALL SCHROEDINGER$DREL(GID,NR,THIS%ATOM%AEPOT,E,DREL)
          G(:)=0.D0          
          CALL SCHROEDINGER$SPHERICAL(GID,NR,THIS%ATOM%AEPOT,DREL,0,G,L,E,1,PHI)
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RCOV,AEPHASE(L+1,IE))
          IF(NPRO.GT.0) THEN
            G(:)=0.D0          
            CALL ATOMLIB_PAWDER(GID,NR,L,E,THIS%PSPOT,NPRO,PRO,DH,DO,G,PHI)
          END IF
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RCOV,PAWPHASE(L+1,IE))
          PAWPHASE(L+1,IE)=PAWPHASE(L+1,IE)+DPHASE
        ENDDO
        DEALLOCATE(DH)
        DEALLOCATE(DO)
        DEALLOCATE(PRO)
      ENDDO
!
!     ==========================================================================
!     ==  WRITE RESULT TO FILE                                                ==
!     ==========================================================================
      rewind(nfil)
      DO IE=1,NE
        E=EMIN+DE*REAL(IE-1,KIND=8)
        WRITE(NFIL,FMT='(20F10.5)')E,AEPHASE(:,IE),PAWPHASE(:,IE)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_MAKEISCATT(AEZ,NB,NC,LOFI,NNOFI,LNX,LOX,ISCATT)
!     **************************************************************************
!     **                                                                      **
!     ** THE PARAMETER ISCATT DETERMINES WHETHER A PARTIAL WAVE IS ATTRIBUTED **
!     ** TO A SEMI-CORE, VALENCE, OR SCATTERING STATE                         **
!     **   ISCATT=-2    STATE BELOW HIGHEST CORE STATE                        **
!     **   ISCATT=-1    SEMI CORE STATE                                       **
!     **   ISCATT= 0    VALENCE STATE (PHI)                                   **
!     **   ISCATT= 1    FIRST SCATTERING STATE (PHIDOT)                       **
!     **   ISCATT= 2    SECOND SCATTERING STATE (PHIDOTDOT)                   **
!     **                                                                      **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: AEZ
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: NC
      INTEGER(4),INTENT(IN)  :: LOFI(NB)
      INTEGER(4),INTENT(IN)  :: NNOFI(NB)
      INTEGER(4),INTENT(IN)  :: LNX
      INTEGER(4),INTENT(IN)  :: LOX(LNX)
      INTEGER(4),INTENT(OUT) :: ISCATT(LNX)
      INTEGER(4),ALLOCATABLE :: NCL(:)
      INTEGER(4)             :: LX
      INTEGER(4)             :: IB,L,LN
      INTEGER(4)             :: ISVAR
!     **************************************************************************
      LX=MAXVAL(LOX)
!
!     == DETERMINE HIGHEST CORE STATE FOR EACH ANGULAR MOMENTUM ================
      ALLOCATE(NCL(0:LX))
      NCL(:)=0
      DO IB=1,NC
        L=LOFI(IB)
        NCL(L)=MAX(NCL(L),IB)
      ENDDO
!
!     == DIVIDE PARTIAL WAVES TO VALENCE AND SCATTERING TYPES ==================
      DO L=0,LX
        IF(L.EQ.0) THEN
          CALL PERIODICTABLE$GET(NINT(AEZ),'#NODES(S)',ISVAR)
        ELSE IF(L.EQ.1) THEN
          CALL PERIODICTABLE$GET(NINT(AEZ),'#NODES(P)',ISVAR)
        ELSE IF(L.EQ.2) THEN
          CALL PERIODICTABLE$GET(NINT(AEZ),'#NODES(D)',ISVAR)
        ELSE IF(L.EQ.3) THEN
          CALL PERIODICTABLE$GET(NINT(AEZ),'#NODES(F)',ISVAR)
        ELSE
          CALL ERROR$STOP('ATOMIC_MAKEISCATT')
        END IF
        IF(NCL(L).GT.0) THEN
          ISVAR=ISVAR-NNOFI(NCL(L))
        ELSE
          ISVAR=ISVAR+1
        END IF
!       == ISVAR IS NOW THE NUMBER OF VALENCE SHELLS
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          ISCATT(LN)=-ISVAR+1
          ISVAR=ISVAR-1
        ENDDO
      ENDDO
      DEALLOCATE(NCL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_WRITEPHI(FILE,GID,NR,NPHI,PHI)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        IF(R(IR).GT.3.D0.AND.MAXVAL(ABS(PHI(IR,:))).GT.1.D+3) EXIT
        WRITE(100,FMT='(F15.10,2X,20(F25.15,2X))')R(IR),PHI(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE ATOMIC_UNSCREEN(GID,NR,RBOX,AEZ,AERHO,PSRHO,PSPOT,RCSM,VADD)
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      REAL(8)     ,INTENT(IN) :: RBOX
      REAL(8)     ,INTENT(IN) :: AEZ
      REAL(8)     ,INTENT(IN) :: AERHO(NR)
      REAL(8)     ,INTENT(IN) :: PSRHO(NR)
      REAL(8)     ,INTENT(IN) :: PSPOT(NR)
      REAL(8)     ,INTENT(IN) :: RCSM
      REAL(8)     ,INTENT(OUT):: VADD(NR)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: AUX(NR),AUX1(NR),SVAR
      REAL(8)                 :: POT(NR),potxc(nr),poth(nr)
      REAL(8)                 :: PI,Y0
      REAL(8)                 :: QLM
      REAL(8)                 :: ALPHA,CL
      REAL(8)                 :: GRHO(NR)
      INTEGER(4)              :: IR
      REAL(8)                 :: RH,GRHO2,VXC,VGXC,EXC,DUMMY1,DUMMY2,DUMMY3
real(8):: auxarr(nr,5)
!     ************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
auxarr(:,1)=aerho
auxarr(:,2)=psrho
CALL setup_WRITEPHI('TEST',GID,NR,2,auxarr(:,:2))
!
AUX(:)=PSRHO(:)*R(:)**2
CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
PRINT*,'PSRHO CHARGE ',QLM*Y0*4.D0*PI
 AUX(:)=AERHO(:)*R(:)**2
CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
PRINT*,'AERHO CHARGE ',QLM*Y0*4.D0*PI
CALL RADIAL$INTEGRAl(GID,NR,AUX,qlm)
PRINT*,'AERHO CHARGE ',QLM*Y0*4.D0*PI,rbox,r(nr)
!
!     ========================================================================
!     == MOMENT OF DIFFERENCE CHARGE DENSITY                                ==
!     ========================================================================
      AUX(:)=(AERHO(:)-PSRHO(:))*R(:)**2
!     == the integration must be performed outward to the end, because =======
!     == small deviations lead to long-range tail in vadd ====================
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
! CALL RADIAL$INTEGRAl(GID,NR,AUX,qlm) 
      QLM=QLM-AEZ*Y0    ! CHARGE =QM/Y0
PRINT*,'CHARGE ',QLM/Y0,AEZ
!
!     ========================================================================
!     == ADD COMPENSATION DENSITY AND DETERMINE ELECTROSTATIC POTENTIA      ==
!     ========================================================================
      ALPHA=1.D0/RCSM**2
      CALL GAUSSN(0,ALPHA,CL)
      SVAR=QLM*CL
      AUX(:)=PSRHO(:)+SVAR*EXP(-ALPHA*R(:)**2)
 CALL RADIAL$INTEGRATE(GID,NR,AUX*R(:)**2,AUX1)
 CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
PRINT*,'QLM ',QLM*Y0*4.D0*PI
      CALL RADIAL$POISSON(GID,NR,0,AUX,POTh)
!
!     ========================================================================
!     == EXCHANGE AND CORRELATION                                           ==
!     ========================================================================
      CALL RADIAL$DERIVE(GID,NR,PSRHO(:),GRHO)
      DO IR=1,NR
        RH=PSRHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        POTxc(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      POTXC(:)=POTXC(:)-AUX(:)
      IF(R(1).GT.1.D+10) THEN
         POTXC(:)=POTXC(:)-2.D0/R(:)*GRHO(:)
      ELSE
        POTXC(2:)=POTXC(2:)-2.D0/R(2:)*GRHO(2:)
        POTXC(1)=POTXC(1)-2.D0/R(2)*GRHO(2)
      END IF
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOR LOW DENSITIES                                  ==
!     == ATTENTION: THIS CUTOFF MUST BE COMPLETELY CONSISTENT WITH            ==
!     ==            ATOMLIB$BOXVOFRHO                                         ==
!     ==========================================================================
      DO IR=1,NR
        IF(AERHO(IR)*Y0.LT.1.D-6) POTXC(IR)=0.D0
      ENDDO
!
!     ========================================================================
!     == VADD                                                               ==
!     ========================================================================
      pot=poth+potxc
      VADD(:)=PSPOT(:)-POT(:)
auxarr(:,1)=pspot
auxarr(:,2)=pot
auxarr(:,3)=vadd
CALL SETUP_WRITEPHI('TESTPOT',GID,NR,3,AUXARR(:,:3))
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE BIORTHOMATRICES(GID,NR,RBOX,LNX,LOX,PSPHI,PRO,TRANSPHI,TRANSPRO)
!     **                                                                  **
!     ** DETERMINES THE MATRICES TRANSPHI AND TRANSPRO SUCH THAT          **
!     **     |PHI-BAR>:=|PHI>TRANSSPHI                                    **
!     **     |PRO-BAR>:=|PRO>TRANSSPRO                                    **
!     **  OBEY  <PHIBAR(I)|PROBAR(J)>=DELTA(I,J)    (KRONECKER DELTA)     **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      REAL(8)    ,INTENT(IN)     :: RBOX
      INTEGER(4) ,INTENT(IN)     :: LNX
      INTEGER(4) ,INTENT(IN)     :: LOX(LNX)
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,LNX)
      REAL(8)    ,INTENT(IN)     :: PRO(NR,LNX)
      REAL(8)    ,INTENT(OUT)    :: TRANSPHI(LNX,LNX)
      REAL(8)    ,INTENT(OUT)    :: TRANSPRO(LNX,LNX)
      INTEGER(4)                 :: LN,LN1,LN2
      REAL(8)                    :: AUX(NR),AUX1(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: SVAR
      LOGICAL(4),PARAMETER       :: TTEST=.TRUE.
      REAL(8)                    :: PROPSI(LNX,LNX)
      REAL(8)                    :: TRANSPROPSI(LNX,LNX)
!     *********************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     =====================================================================
!     == CALCULATE INITIAL VIOLATION OF BIORTHOGONALITY                  ==
!     =====================================================================
      PROPSI(:,:)=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=PRO(:,LN1)*PSPHI(:,LN2)*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)       
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)       
          PROPSI(LN1,LN2)=SVAR
        ENDDO
      ENDDO
!
!     =====================================================================
!     == COLLECT TRANSFORMATION MATRIX BETWEEN NEW AND OLD               ==
!     =====================================================================
      TRANSPRO(:,:)=0.D0
      TRANSPHI(:,:)=0.D0
      DO LN=1,LNX
        TRANSPRO(LN,LN)=1.D0
        TRANSPHI(LN,LN)=1.D0
      ENDDO
      TRANSPROPSI(:,:)=PROPSI(:,:)
!
      DO LN1=1,LNX
        DO LN2=1,LN1-1
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!         == ORTHOGONALIZE PARTIAL WAVES ===================================
          SVAR=TRANSPROPSI(LN2,LN1)
          TRANSPHI(:,LN1)   =TRANSPHI(:,LN1)-TRANSPHI(:,LN2)*SVAR
          TRANSPROPSI(:,LN1)=TRANSPROPSI(:,LN1)-TRANSPROPSI(:,LN2)*SVAR          
!         == ORTHOGONALIZE PROJECTOR=====================================
          SVAR=TRANSPROPSI(LN1,LN2)
          TRANSPRO(:,LN1)   =TRANSPRO(:,LN1)-TRANSPRO(:,LN2)*SVAR
          TRANSPROPSI(LN1,:)=TRANSPROPSI(LN1,:)-TRANSPROPSI(LN2,:)*SVAR          
        ENDDO
        SVAR=TRANSPROPSI(LN1,LN1)
        TRANSPRO(:,LN1)=TRANSPRO(:,LN1)/SVAR
        TRANSPROPSI(LN1,:)=TRANSPROPSI(LN1,:)/SVAR
      ENDDO
!
!     =====================================================================
!     == CHECK RESULT                                                    ==
!     =====================================================================
      TRANSPROPSI(:,:)=MATMUL(TRANSPOSE(TRANSPRO),MATMUL(PROPSI,TRANSPHI))
      DO LN=1,LNX
        TRANSPROPSI(LN,LN)=TRANSPROPSI(LN,LN)-1.D0
      ENDDO
      SVAR=MAXVAL(TRANSPROPSI)
      IF(SVAR.GT.1.D-5) THEN
DO LN=1,LNX
WRITE(*,FMT='(20E10.2)')TRANSPROPSI(:,LN)
ENDDO
        CALL ERROR$MSG('BIORTHOGONALIZATION INACCURATE')
        CALL ERROR$R8VAL('MAX. DEV.',SVAR)
        CALL ERROR$STOP('BIORTHOMATRICES')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMIC_MAKEPSPHI_HBS(GID,NR,LNX,LOX,EOFLN,RC,LAMBDA,PSPOT &
     &                               ,rbnd,PSPHI,TPSPHI)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
use strings_module
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)     :: GID
      INTEGER(4),INTENT(IN)     :: NR
      INTEGER(4),INTENT(IN)     :: LNX
      INTEGER(4),INTENT(IN)     :: LOX(LNX)
      REAL(8)   ,INTENT(IN)     :: EOFLN(LNX)
      REAL(8)   ,INTENT(IN)     :: LAMBDA(LNX)
      REAL(8)   ,INTENT(IN)     :: RC(LNX)
      REAL(8)   ,INTENT(IN)     :: PSPOT(NR)
      REAL(8)   ,INTENT(IN)     :: rbnd
      REAL(8)   ,INTENT(INOUT)  :: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(INOUT)  :: TPSPHI(NR,LNX)
      INTEGER(4),PARAMETER      :: NITER=100
      INTEGER(4),PARAMETER      :: ISO=0
      REAL(8)   ,PARAMETER      :: TOL=1.D-6
      REAL(8)   ,PARAMETER      :: cmin=1.D-8
      REAL(8)                   :: E
      REAL(8)                   :: PI,Y0
      REAL(8)                   :: AUX(NR),DREL(NR),G(NR),POT(NR),PHI(NR)
      REAL(8)                   :: C(NR)
      REAL(8)                   :: R(NR)
      REAL(8)                   :: X0,XM,Z0,ZM,DX
      REAL(8)                   :: PHIPHASE
      INTEGER(4)                :: IRBND
      INTEGER(4)                :: ISTART,IBI
      INTEGER(4)                :: L
      INTEGER(4)                :: LN,ITER,IR
      logical(4)                :: convg
      real(8)                   :: svar
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      IRBND=0
      DO IR=1,NR-2
        IRBND=IR
        IF(r(ir).GT.RBND) EXIT
      ENDDO
!
!     ==========================================================================
!     ==  CONSTRUCT                                                           ==
!     ==========================================================================
      DO LN=1,LNX
        L=LOX(LN)
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSPHI(:,LN),RBND,PHIPHASE)
        E=EOFLN(LN)
        C(:)=EXP(-(R(:)/RC(LN))**lambda(ln))
!       == cut off c, if it falls below minimum ================================
        svar=rc(ln)*(-log(cmin))**(1.d0/lambda(ln))
        irbnd=0
        do ir=1,nr
          irbnd=ir
          if(r(ir).gt.svar) exit
        enddo
        C(:)=(C(:)-C(IRBND))/(1.D0-C(IRBND))
        C(IRBND:)=0.D0
!
!       == find outermost point considering the potential ======================
        do ir=nr-2,1,-1
          if(ir.lt.irbnd) exit
          svar=tpsphi(ir,ln)+(pspot(ir)*y0-e)*psphi(ir,ln)
          svar=abs(svar/psphi(ir,ln))
          if(svar.gt.1.d-5) then
            irbnd=ir
            exit
          end if
        enddo
!
!       == loop to find pseudo partial wave ====================================
        ISTART=1
        X0=0.d0
        DX=1.D-2
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        DO ITER=1,NITER
          POT(:)=PSPOT(:)+X0*C(:)
          DREL(:)=0.D0
          G(:)=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,ISO,G,L,E,1,PHI)
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RBND,Z0)
          Z0=PHIPHASE-z0
          convg=(ABS(2.D0*DX).LE.TOL)
          IF(convg) EXIT
          CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        ENDDO
        IF(.NOT.CONVG) THEN
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$STOP('ATOMIC_MAKEPSPHI_HBS')
        END IF
        phi(:)=phi(:)/phi(irbnd)*psphi(irbnd,ln)
        psphi(:,ln)=phi(:)
        TPSPHI(:,LN)=(E-POT(:)*Y0)*phi(:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMIC_MAKEPSPHI(GID,NR,RC,L,NB,PSPHI,TPSPHI)
!     **                                                                      **
!     **  PSEUDIZES A WAVE FUNCTION BY MATCHING SPHERICAL BESSEL FUNCTIONS    **
!     **                                                                      **
!     **  FIRST WE FIND THE VALUES K_I SO THAT THE SPHERICAL BESSEL           **
!     **  FUNCTION J_L(K_I*R) MATCHES THE LOGARITHMIC DERIVATIVE              **
!     **  OF THE FUNCTION TO BE PSEUDIZED.                                    **
!     **                                                                      **
!     **  THE FIRST TWO ARE THEN MATCHED SO THAT THE PHI AND TPHI             **
!     **  ARE CONTINUOUS. THE DERIVATIVE OF PHI IS AUTOMATICALLY              **
!     **  CONTINUOUS THROUGH THE CORRECT CHOICE OF THE K_I.                   **
!     **                                                                      **
!     **  ATTENTION! SINCE WE MATCH TO TAYLOR COEFFICIENTS, ONE SHOULD        **
!     **  ALSO MATCH ENERGY DERIVATIVES OF THE SPHERICAL BESSEL FUNCTIONS.    **
!     **                                                                      **
!     **  NOTE THAT THE TAYLOR EXPANSION OF THE WAVE FUNCTIONS IN ENERGY      **
!     **  BEHAVE DIFFERENT FROM THE ENERGY DEPENDENT SOLUTION.                **
!     **  THE TAYLOR EXPANSION CAN CREATE TWO NODES AT ONCE, WHICH            **
!     **  CAUSES DIFFICULTIES.                                                **
!     **                                                                      **
!     **                                                                      ** 
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)     :: GID
      INTEGER(4),INTENT(IN)     :: NR
      INTEGER(4),INTENT(IN)     :: L
      INTEGER(4),INTENT(IN)     :: NB
      REAL(8)   ,INTENT(IN)     :: RC
      REAL(8)   ,INTENT(INOUT)  :: PSPHI(NR,NB)
      REAL(8)   ,INTENT(INOUT)  :: TPSPHI(NR,NB)
      INTEGER(4),PARAMETER      :: NZEROS=2
      INTEGER(4),PARAMETER      :: NITER=100
      INTEGER(4),PARAMETER      :: NX=1000
      LOGICAL(4),PARAMETER      :: TDIFFERENTIABLELAPLACIAN=.FALSE.
      REAL(8)   ,PARAMETER      :: XMAX=30.D0
      REAL(8)                   :: R(NR)
      INTEGER(4)                :: IR,I,J,ITER
      REAL(8)                   :: TOL=1.D-12
      REAL(8)                   :: PI
      REAL(8)                   :: VAL,VALT,DERT,VALJ,VALTJ
      REAL(8)                   :: X(NX)
      REAL(8)                   :: JL(NX),DJLDR(NX),D2JLDR2(NX),PHASE(NX)
      REAL(8)                   :: TJL(NX)
      REAL(8)                   :: JLOFKR(NR,NZEROS),TJLOFKR(NR,NZEROS)
      REAL(8)                   :: DJLOFKR(NR,NZEROS),D2JLOFKR(NR,NZEROS),RDUMMY(NR,NZEROS)
      REAL(8)                   :: AUX1(NR),AUX2(NR)
      REAL(8)                   :: NN
      INTEGER(4)                :: ISTART,IBI
      REAL(8)                   :: X0,Y0,DX,XM,YM
      REAL(8)                   :: KI(NZEROS)
      REAL(8)                   :: A11,A12,A21,A22,DET
      REAL(8)                   :: AINV11,AINV12,AINV21,AINV22
      REAL(8)                   :: V1,V2,C1,C2,C3,C1A,C2A,C1B,C2B
      REAL(8)                   :: DTJ1,DTJ2,DTJ3
      REAL(8)                   :: SVAR1,SVAR2,FAC,SVAR
      INTEGER(4)                :: IB
      REAL(8)                   :: SUPPHI(NR,NB)
      REAL(8)                   :: SUPTPHI(NR,NB)
      REAL(8)                   :: SHIFT
      REAL(8)                   :: X1,XDEX
      REAL(8)                   :: OFFSETNN
      INTEGER(4)                :: GIDX
      INTEGER(4)                :: IRMATCH
!     **************************************************************************
!CALL SETUP_WRITEPHI('TESTCONSTRUCTPSPHI1',GID,NR,NB,PSPHI)
      PI=4.D0*ATAN(1.D0)
      CALL RADIAL$NEW('SHLOG',GIDX)
      XDEX=1.D-5
      X1=XMAX/XDEX/REAL(NX)
      CALL RADIAL$SETR8(GIDX,'R1',X1)
      CALL RADIAL$SETR8(GIDX,'DEX',XDEX)
      CALL RADIAL$SETI4(GIDX,'NR',NX)
      CALL RADIAL$R(GIDX,NX,X)
!
!     ==========================================================================
!     ==  MAKE A COPY OF THE INPUT WAVE FUNCTIONS                             ==
!     ==========================================================================
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$XOFR(GID,RC,VAL)
      IRMATCH=INT(VAL)
      SUPPHI(:,:)=PSPHI(:,:)
      SUPTPHI(:,:)=TPSPHI(:,:)
!
!     =======================================================================
!     == CALCULATE BESSEL FUNCTION                                         ==
!     =======================================================================
      DO IR=1,NX
        CALL SPECIALFUNCTION$BESSEL(L,X(IR),JL(IR)) 
      ENDDO
!
!     ==========================================================================
!     == CALCULATE PHASESHIFT OF BESSEL FUNCTION                              ==
!     ==========================================================================
      CALL RADIAL$VERLETD1(GIDX,NX,JL,DJLDR)
      NN=0.D0
      DO IR=2,NX
        IF(JL(IR)*JL(IR-1).LT.0.D0)NN=NN+1.D0
        PHASE(IR)=0.5D0-ATAN(X(IR)/RC*DJLDR(IR)/JL(IR))/PI+NN
      ENDDO
!     == AVOID GLITCH
      PHASE(1)=PHASE(2)
      PHASE(NX)=PHASE(NX-1)
!
!     ==========================================================================
!     == DETERMINE OFFSET IN THE NUMBER OF NODES                              ==
!     ==========================================================================
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,SUPPHI(:,1),RC,SHIFT)
      OFFSETNN=REAL(INT(SHIFT),KIND=8)
!     
!     ==========================================================================
!     == PSEUDIZE EVERY FUNCTION ON THE SUPPORT ENERGY GRID                   ==
!     ==========================================================================
      DO IB=1,NB
!       ========================================================================
!       == DETERMINE K-VALUES FOR WHICH THE LOGARITHIC DERIVATIVE MATCHES     ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,SUPPHI(:,IB),RC,SHIFT)
        SHIFT=SHIFT-OFFSETNN
        IF(SHIFT.LT.PHASE(1)) THEN
          CALL ERROR$MSG('MATCHING INCLUDING FIRST BESSEL FUNCTION NOT POSSIBLE')
          CALL ERROR$MSG('TRY TO INCREASE MATCHING RADIUS')
          CALL ERROR$I4VAL('IB',IB)
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('SHIFT',SHIFT)
          CALL ERROR$R8VAL('PHASE(1)',PHASE(1))
          CALL ERROR$R8VAL('OFFSETNN',OFFSETNN)
          CALL ERROR$R8VAL('RC',RC)
          CALL ERROR$STOP('CONSTRUCTPSPHI')
PRINT*,'WARNING SHIFT WAS REQUIRED ',L,IB,PHASE(1),SHIFT
           SHIFT=SHIFT+REAL(INT(PHASE(1)-SHIFT+1.D0))
        END IF
!
!       ========================================================================
!       == SEARCH BESSEL FUNCTIONS JL(K*R) WHICH HAVE IDENTICAL LOGARITHMIC   ==
!       == DERIVATIVES TO THE NODELESS FUNCTION AT THE RADIUS RCS             ==
!       ========================================================================
        DO I=1,NZEROS
          ISTART=1
          X0=0.4D0+1.D-3
          DX=0.1D0
          CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
          DO ITER=1,NITER
            CALL RADIAL$VALUE(GIDX,NX,PHASE,X0,Y0)
            Y0=Y0-SHIFT
            CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
            IF(X0.GT.XMAX) THEN
              CALL ERROR$MSG('RESULTING K OUT OF RANGE')
              CALL ERROR$STOP('CONSTRUCTPSPHI')
            END IF
            IF(ABS(DX).LT.TOL) EXIT
          ENDDO
          IF(ABS(DX).GT.TOL) THEN
            CALL ERROR$MSG('LOOP NOT CONVERGED')
            CALL ERROR$STOP('CONSTRUCTPSPHI')
          END IF
          KI(I)=X0/RC
          IF(KI(I).LT.0.D0) THEN
            CALL ERROR$MSG('KI NEGATIVE')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('IB',IB)
            CALL ERROR$I4VAL('I',I)
            CALL ERROR$R8VAL('KI(I)',KI(I))
            CALL ERROR$STOP('CONSTRUCTPSPHI')
          END IF
          SHIFT=SHIFT+1.D0
        ENDDO
PRINT*,'KI ',KI
!
!       ========================================================================
!       == PLACE BESSELFUNCTIONS ON RADIAL GRID                               ==
!       ========================================================================
        DO I=1,NZEROS
          DO IR=1,NR
            CALL SPECIALFUNCTION$BESSEL(L,KI(I)*R(IR),JLOFKR(IR,I)) 
          ENDDO
          TJLOFKR(:,I)=0.5D0*KI(I)**2*JLOFKR(:,I)
        ENDDO
!
!       ========================================================================
!       == CONSTRUCT NEW PSEUDO WAVE FUNCTION                                 ==
!       ========================================================================
        DO I=1,NZEROS
!         == SELECT AND NORMALIZE THIS AND HIGHER BESSEL FUNCTIONS WITH ========
!         == RESPECT TO PROPERTY ===============================================
          DO J=I,NZEROS
            IF(I.EQ.1) THEN
              CALL RADIAL$VALUE(GID,NR,JLOFKR(:,J),RC,SVAR)
            ELSE IF(I.EQ.2) THEN
              CALL RADIAL$VALUE(GID,NR,TJLOFKR(:,J),RC,SVAR)
            ELSE IF(I.EQ.3) THEN
              CALL RADIAL$DERIVATIVE(GID,NR,TJLOFKR(:,J),RC,SVAR)     
              SVAR=SVAR*1.D+3
            END IF
            IF(ABS(SVAR).LT.1.D-7) THEN
              CALL ERROR$MSG('PROBLEM WITH DIVIDE-BY ZERO')
              CALL ERROR$MSG('CHANGE MATCHING RADIUS FOR PS-PARTIALWAVE CONSTRUCTION')
              CALL ERROR$I4VAL('L',L)
              CALL ERROR$I4VAL('IB',IB)
              CALL ERROR$R8VAL('RC',RC)
              CALL ERROR$I4VAL('I',I)
              CALL ERROR$I4VAL('J',J)
              CALL ERROR$STOP('CONSTRUCTPSPHI')
            END IF
            JLOFKR(:,J) =JLOFKR(:,J)/SVAR
            TJLOFKR(:,J)=TJLOFKR(:,J)/SVAR
          ENDDO
!
!         == PROJECT PROPERTY FROM HIGHER BESSEL FUNCTIONS =====================
          DO J=I+1,NZEROS
            JLOFKR(:,J) =JLOFKR(:,J)-JLOFKR(:,I)
            TJLOFKR(:,J)=TJLOFKR(:,J)-TJLOFKR(:,I)
          ENDDO
!
!         == EXTRACT MISSING PART OF THIS PROPERTY IN THE PARTIAL WAVE =========
          AUX1(:)=SUPPHI(:,IB)
          AUX2(:)=SUPTPHI(:,IB)
          DO J=1,I-1
            AUX1(:)=AUX1(:)-JLOFKR(:,J)
            AUX2(:)=AUX2(:)-TJLOFKR(:,J)
          ENDDO
          IF(I.EQ.1) THEN
            CALL RADIAL$VALUE(GID,NR,AUX1,RC,SVAR)
          ELSE IF(I.EQ.2) THEN
            CALL RADIAL$VALUE(GID,NR,AUX2,RC,SVAR)
          ELSE IF(I.EQ.3) THEN
            CALL RADIAL$DERIVATIVE(GID,NR,AUX2,RC,SVAR)     
            SVAR=SVAR*1.D+3
          END IF
!         == SCALE THIS BESSEL FUNCTION TO MATCH THE MISSING PART ==============
          JLOFKR(:,I) =JLOFKR(:,I)*SVAR
          TJLOFKR(:,I)=TJLOFKR(:,I)*SVAR
        ENDDO
!
!       ========================================================================
!       ==                                                                    ==
!       ========================================================================
        DO I=1,NZEROS
          DO J=1,I-1
            JLOFKR(:,I)=JLOFKR(:,I)+JLOFKR(:,J)
            TJLOFKR(:,I)=TJLOFKR(:,I)+TJLOFKR(:,J)
          ENDDO
        ENDDO
        SUPPHI(:IRMATCH,IB) =JLOFKR(:IRMATCH,NZEROS)
        SUPTPHI(:IRMATCH,IB)=TJLOFKR(:IRMATCH,NZEROS)
!
!!$IF(IB.EQ.1)OPEN(100,FILE='JLSD1.DAT')
!!$IF(IB.EQ.2)OPEN(100,FILE='JLSD2.DAT')
!!$IF(IB.EQ.3)OPEN(100,FILE='JLSD3.DAT')
!!$REWIND(100)
!!$DO IR=1,NR
!!$  WRITE(100,FMT='(20F20.8)')R(IR),PSPHI(IR,IB),JLOFKR(IR,:) &
!!$ &                               ,TPSPHI(IR,IB),TJLOFKR(IR,:)
!!$ENDDO
!!$CLOSE(100)
      ENDDO
!CALL WRITEPHI('TESTCONSTRUCTPSPHI2',GID,NR,NB,SUPPHI)
!STOP 'FORCED IN CONSTRUCTPSPHI'
      PSPHI(:IRMATCH,:)=SUPPHI(:IRMATCH,:)
      TPSPHI(:IRMATCH,:)=SUPTPHI(:IRMATCH,:)
      RETURN
      END

!
!     ......................................................................
      SUBROUTINE ATOMIC_PSEUDIZE(GID,NR,POW,tval,VAL0_,RC,AEF,PSF)
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      REAL(8)    ,INTENT(IN)     :: POW   !POWER
      logical(4) ,intent(in)     :: tval
      REAL(8)    ,INTENT(IN)     :: VAL0_
      REAL(8)    ,INTENT(IN)     :: RC
      REAL(8)    ,INTENT(IN)     :: AEF(NR)
      REAL(8)    ,INTENT(OUT)    :: PSF(NR)
      REAL(8)                    :: VAL0
      REAL(8)                    :: VAL,DER
      REAL(8)                    :: V1,V2
      REAL(8)                    :: A11,A12,A21,A22
      REAL(8)                    :: DET
      REAL(8)                    :: AINV11,AINV12,AINV21,AINV22
      REAL(8)                    :: C0,C1,C2
      INTEGER(4)                 :: IR,IR0
      REAL(8)                    :: R(NR)
      LOGICAL(4),PARAMETER       :: TTEST=.FALSE.
      REAL(8)                    :: SVAR
      REAL(8)                    :: AUX(NR)
!     *********************************************************************
      CALL RADIAL$VALUE(GID,NR,AEF,RC,VAL)
      CALL RADIAL$DERIVATIVE(GID,NR,AEF,RC,DER)
      if(tval) then
        val0=val0_
      else
        val0=val-rc*der/pow
      end if
!     == EQUATION SYSTEM A*C=V      
      V1=VAL-VAL0
      V2=DER
      A11=RC**POW
      A12=RC**(POW+2.D0)
      A21=POW       *RC**(POW-1.D0)
      A22=(POW+2.D0)*RC**(POW+1.D0)
!     == INVERT MATRIX
      DET=A11*A22-A12*A21
      AINV11=A22/DET
      AINV22=A11/DET
      AINV12=-A12/DET
      AINV21=-A21/DET
!     == SOLVE EQUATION SYSTEM
      C0=VAL0
      C1=AINV11*V1+AINV12*V2
      C2=AINV21*V1+AINV22*V2
if(.not.tval)print*,'--c2-- ',c2,val0
!     == CALCULATE PS FUNCTION
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IF(R(IR).GT.RC) THEN
          IR0=IR
          EXIT
        END IF
      ENDDO
      PSF(:IR0-1)=C0+R(:IR0-1)**POW*(C1+C2*R(:IR0-1)**2)
      PSF(IR0:)=AEF(IR0:)
      IF(TTEST) THEN
        SVAR=C0+C1*RC**POW+C2*RC**(POW+2.D0)
        PRINT*,'TEST PESUDIZE 1 ',SVAR,VAL,SVAR-VAL
        SVAR=POW*C1*RC**(POW-1.D0)+(POW+2.d0)*C2*RC**(POW+1.D0)
        PRINT*,'TEST PESUDIZE 2 ',SVAR,DER,SVAR-DER
        SVAR=C0
        PRINT*,'TEST PESUDIZE 3 ',SVAR,VAL0,SVAR-VAL0
        AUX(:)=C0+R(:)**POW*(C1+C2*R(:)**2)
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$AECORE(ISP_,NRX_,AECORE_)
!LEGACY CODE! -> SETUP$GETR8A('AECORE'
!     ******************************************************************
!     **  RETURN AE CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: AECORE_(NRX_)
!     ******************************************************************
      CALL ERROR$MSG('CODE IS MARKED FOR DELETION')
      CALL ERROR$MSG('USE SETUP$GETR8A INSTEAD')
      CALL ERROR$STOP('SETUP$AECORE')
      CALL SETUP$ISELECT(ISP_)
      AECORE_(:)=THIS%AECORE(:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$PSCORE(ISP_,NRX_,PSCORE_)
!LEGACY CODE! -> SETUP$GETR8A('PSCORE'
!     ******************************************************************
!     **  RETURN PS CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: PSCORE_(NRX_)
!     ******************************************************************
      CALL ERROR$MSG('CODE IS MARKED FOR DELETION')
      CALL ERROR$MSG('USE SETUP$GETR8A INSTEAD')
      CALL ERROR$STOP('SETUP$PSCORE')
      CALL SETUP$ISELECT(ISP_)
      PSCORE_(:)=THIS%PSCORE(:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$VBAR(ISP_,NRX_,VBAR_)
!LEGACY CODE! -> SETUP$GETR8A('VADD'
!     ******************************************************************
!     **  RETURN PS CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: VBAR_(NRX_)
!     ******************************************************************
      CALL ERROR$MSG('CODE IS MARKED FOR DELETION')
      CALL ERROR$MSG('USE SETUP$GETR8A INSTEAD')
      CALL ERROR$STOP('SETUP$VBAR')
      CALL SETUP$ISELECT(ISP_)
      VBAR_(:)=THIS%VADD(:)
      RETURN  
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_RENEWLOCORB()
!     **************************************************************************
!     **                                                                      **
!     **  DRIVER ROUTINE FOR SETUP_CHIFROMPHI                                 **
!     **                                                                      **
!     **  CHECKS WHETHER AN UPDATE IS NECESSARY; REALLOCATES ARRAYS AND       **
!     **  COMPUTES NEW TRANSFORMATION MATRICES                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ORBRAD
      USE PERIODICTABLE_MODULE
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LN,L,I,NR,LXX
      REAL(8)   ,ALLOCATABLE :: RAD(:)
      CHARACTER(3)           :: ID ! CAN BE OLD OR NEW
      REAL(8)                :: RASA
!     **************************************************************************
      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
      CALL ERROR$STOP('SETUP_RENEWLOCORB')
      IF(THIS%LOCORBINI) RETURN
      THIS%LOCORBINI=.TRUE.
                            CALL TRACE$PUSH('SETUP_RENEWLOCORB')
!     == LOCORBINI IS RESET TO FALSE, WHEVEVER THE DATA ARE CHANGED VIA SETUP$SET
PRINT*,'A1',THIS%LOCORBNOFL
PRINT*,'A2',THIS%LOCORBLNX
!
!     == REALLOCATE DATA =======================================================
      IF(SUM(THIS%LOCORBNOFL).NE.THIS%LOCORBLNX) THEN
PRINT*,'B'
        THIS%LOCORBLNX=SUM(THIS%LOCORBNOFL)
        IF(ASSOCIATED(THIS%LOCORBLOX)) DEALLOCATE(THIS%LOCORBLOX)
        IF(ASSOCIATED(THIS%LOCORBAMAT)) DEALLOCATE(THIS%LOCORBAMAT)
        IF(ASSOCIATED(THIS%LOCORBBMAT)) DEALLOCATE(THIS%LOCORBBMAT)
        ALLOCATE(THIS%LOCORBLOX(THIS%LOCORBLNX))
        ALLOCATE(THIS%LOCORBAMAT(THIS%LNX,THIS%LOCORBLNX))
        ALLOCATE(THIS%LOCORBBMAT(THIS%LOCORBLNX,THIS%LNX))
        LN=0
        DO L=0,3
          DO I=1,THIS%LOCORBNOFL(L+1)
            LN=LN+1
            THIS%LOCORBLOX(LN)=L
          ENDDO
        ENDDO
      END IF
PRINT*,'C'
!
!     == RECALCULATE LOCAL ORBITALS ============================================
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      ALLOCATE(RAD(THIS%LOCORBLNX))
      DO LN=1,THIS%LOCORBLNX
        L=THIS%LOCORBLOX(LN)
        RAD(LN)=THIS%LOCORBRAD(L+1)
      ENDDO
      CALL SETUP_CHIFROMPHI(THIS%GID,NR,THIS%LNX,THIS%LOX,THIS%AEPHI &
     &               ,THIS%LOCORBNOFL,RAD,THIS%LOCORBLNX &
     &               ,THIS%LNX,THIS%LOCORBLOX,THIS%LOCORBAMAT,THIS%LOCORBBMAT)
      DEALLOCATE(RAD)
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_CHIFROMPHI(GID,NR,LNX,LOX,PHI,NORB,RCUT &
     &                           ,LNXCHI,LNXPHI,LOXCHI,AMAT,BMAT)
!     **                                                                      **
!     **  DEFINES PROJECTOR FUNCTIONS FOR LOCAL ORBITALS IN TERMS OF          **
!     **  THE CONVENTIONAL PARTIAL WAVE PROJECTOR FUNCTIONS                   **
!     **                                                                      **
!     **     <P_CHI_I|= SUM_J BMAT(I,J) <P_PHI_J|                             **
!     **     |CHI_I>  = SUM_J |PHI_J> AMAT(J,I)                               **
!     **                                                                      **
!     **  THE LOCAL ORBITALS ARE SELECTED BY NORB WHICH SPECIFIES THE NUMBER  **
!     **  OF LOCAL ORBITALS PER ANGULAR MOMENTUM CHANNEL.                     **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID     ! GRID ID
      INTEGER(4),INTENT(IN) :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LNX     ! #(PARTIAL WAVES) 
      INTEGER(4),INTENT(IN) :: LOX(LNX)! ANGULAR MOMENTUM PER PARTIAL WAVE
      REAL(8)   ,INTENT(IN) :: PHI(NR,LNX)    ! PARTIAL WAVES
      INTEGER(4),INTENT(IN) :: NORB(4) ! X#(HEAD ORBITALS PER L)
      REAL(8)   ,INTENT(IN) :: RCUT(LNX)  ! EXTENT OF HEAD FUNCTIONS
      INTEGER(4),INTENT(IN) :: LNXCHI  ! #(LOCAL ORBITALS)
      INTEGER(4),INTENT(IN) :: LNXPHI  ! #(PARTIAL WAVES)
      INTEGER(4),INTENT(OUT):: LOXCHI(LNXCHI)  ! ANGULAR MOMENTUM PER LOCAL ORB.
      REAL(8)   ,INTENT(OUT):: BMAT(LNXCHI,LNXPHI) ! MAPPING OF PROJECTORS
      REAL(8)   ,INTENT(OUT):: AMAT(LNXPHI,LNXCHI)     ! MAPPING OF WAVE FUNCTIONS
      REAL(8)   ,ALLOCATABLE:: CHI(:,:)
      REAL(8)   ,ALLOCATABLE:: CHI1(:,:)
      REAL(8)   ,ALLOCATABLE:: A(:,:)
      REAL(8)   ,ALLOCATABLE:: B(:,:)
      REAL(8)   ,ALLOCATABLE:: A1(:,:)
      REAL(8)   ,ALLOCATABLE:: MAT(:,:)
      REAL(8)   ,ALLOCATABLE:: MATINV(:,:)
      REAL(8)   ,ALLOCATABLE:: R(:)        
      REAL(8)   ,ALLOCATABLE:: G(:)        
      REAL(8)   ,PARAMETER  :: RCG=1.D-2
      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:)
      REAL(8)               :: SVAR1,SVAR2
      INTEGER(4)            :: IR
      INTEGER(4)            :: NX,N,LX,L,LN,LNCHI,NOFL,NOFL2,ISVAR
      INTEGER(4)            :: N1,N2,LN1,LN2,L1,L2
      CHARACTER(64)         :: FILE
!     **************************************************************************
                            CALL TRACE$PUSH('SETUP_CHIFROMPHI')
      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
      CALL ERROR$STOP('SETUP_CHIFROMPHI')
!
!     ==========================================================================
!     ==  COLLECT DATA FROM SETUP OBJECT                                      ==
!     ==========================================================================
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  CONSISTENCY CHECKS                                                  ==
!     ==========================================================================
!     == CHECK IF #(LOCAL ORBITALS) IS CONSISTENT WITH ARRAY DIMENSIONS ========
      IF(SUM(NORB).NE.LNXCHI) THEN
        CALL ERROR$MSG('NORB AND LNXCHI INCONSISTENT')
        CALL ERROR$I4VAL('NORB',LNX)
        CALL ERROR$I4VAL('LNXCHI',LNXCHI)
        CALL ERROR$STOP('SETUP_CHIFROMPHI')
      END IF
!
!     == CHECK IF #(PARTIAL WAVES) FROM SETUP OBJECT IS CONSISTENT WITH ========
!     == ARRAY DIMENSIONS ======================================================
      IF(LNX.NE.LNXPHI) THEN
        CALL ERROR$MSG('LNX AND LNXPHI INCONSISTENT')
        CALL ERROR$I4VAL('LNX',LNX)
        CALL ERROR$I4VAL('LNXPHI',LNXPHI)
        CALL ERROR$STOP('SETUP_CHIFROMPHI')
      END IF
!
!     == CHECK IF NUMBER OF LOCAL ORBITALS DOES NOT EXCEED NUMBER OF PARTIAL WAVES
      LX=MAXVAL(LOX)
      DO L=0,LX
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE !CONSIDER ONLY PARTIAL WAVES WITH THE CORRECT L
          NOFL=NOFL+1
        ENDDO
        IF(NORB(L+1).GT.NOFL) THEN
          CALL ERROR$MSG('#(CORRELATED ORBITALS) EXCEEDS #(PARTIAL WAVES)')
          CALL ERROR$STOP('SETUP_CHIFROMPHI')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  START WITH PARTIAL WAVES AS LOCAL ORBITALS                          ==
!     == |CHI_I>=SUM_I |PHI_J>A(J,I)
!     ==========================================================================
      ALLOCATE(CHI(NR,LNX))        
      CHI(:,:)=0.D0
      ALLOCATE(A(LNX,LNX))        
      A(:,:)=0.D0
      DO LN=1,LNX
        L=LOX(LN)
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        A(LN,LN)=1.D0
        CHI(:,LN)=PHI(:,LN)
      ENDDO
!
!     ==========================================================================
!     == MAKE HEAD FUNCTION ANTIBONDING WITH NODE AT RCUT ======================
!     ==========================================================================
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          NOFL=NOFL+1
          IF(NOFL.GT.NORB(L+1)) CYCLE  ! CONSIDER ONLY HEAD FUNCTIONS
!         == SEARCH FOR NEXT PARTIAL WAVE WITH THIS L =========================
          LN1=0
          DO LN2=LN+1,LNX
            IF(LOX(LN2).NE.L) CYCLE  
            LN1=LN2
            EXIT
          ENDDO
          IF(LN1.EQ.0) THEN
            CALL ERROR$MSG('CAN NOT LOCALIZE')
            CALL ERROR$MSG('NO TAIL FUNCTION AVAILABLE')
            CALL ERROR$STOP('SETUP_CHIFROMPHI')
          END IF
!
!         == NOW CONTINUE; LN1 POINTS TO THE NEXT PHI WITH THIS L ==============
!         == IMPOSE NODE CONDITION==============================================
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RCUT(LN),SVAR1)
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN1),RCUT(LN),SVAR2)
          IF(SVAR1.EQ.0.D0.AND.SVAR2.EQ.0.D0) THEN
            CALL ERROR$MSG('PARTIAL WAVES ARE TRUNCATED INSIDE OF RCUT')
            CALL ERROR$MSG('THIS IS A FLAW OF THE IMPLEMENTATION')
            CALL ERROR$MSG('CHOOSE SMALLER RCUT')
            CALL ERROR$STOP('LDAPLUSU_CHIFROMPHI')
          END IF
          CHI(:,LN)=CHI(:,LN)*SVAR2-CHI(:,LN1)*SVAR1
          A(:,LN)=A(:,LN)*SVAR2-A(:,LN1)*SVAR1
!         == ORTHOGONALIZE TO THE LOWER HEAD FUNCTIONS =========================
          DO LN1=1,LN-1
            IF(LOX(LN1).NE.L) CYCLE  
            AUX(:)=CHI(:,LN)*CHI(:,LN1)*R(:)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,MIN(RCUT(LN),RCUT(LN1)),SVAR1)
            CHI(:,LN)=CHI(:,LN)-CHI(:,LN1)*SVAR1
            A(:,LN)=A(:,LN)-A(:,LN1)*SVAR1
          ENDDO
!         == NORMALIZE HEAD FUNCTION =============================================
          AUX(:)=CHI(:,LN)**2*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RCUT(LN),SVAR1)
          SVAR1=1.D0/SQRT(SVAR1)
          CHI(:,LN)=CHI(:,LN)*SVAR1
          A(:,LN)=A(:,LN)*SVAR1
        ENDDO
      END DO          
!
!     ==========================================================================
!     == DELOCALIZE TAIL FUNCTIONS                                            ==
!     ==========================================================================
!     == DEFINE WEIGHT FUNCTION       
      ALLOCATE(G(NR))
      G(:)=EXP(-(R(:)/RCG)**2)
!     ==  MINIMIZE SUCCESSIVELY THE WEIGHT OF <CHI2+X*CHI1|W|CHI2+X*CHI1> ======
!     ==  I.E. <CHI1|W|CHI2>+X<CHI1|W|CHI1> ====================================
      ALLOCATE(CHI1(NR,LNX))
      CHI1(:,:)=CHI(:,:)
      ALLOCATE(A1(LNX,LNX))        
      A1(:,:)=A(:,:)
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          NOFL=NOFL+1
!         == MINIMIZE CONTRIBUTION NEAR THE CENTER 
          NOFL2=0
          DO LN2=1,LN-1
            IF(LOX(LN2).NE.L) CYCLE
            NOFL2=NOFL2+1
            AUX(:)=G(:)*CHI1(:,LN)*CHI1(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
            AUX(:)=G(:)*CHI1(:,LN2)**2*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
            SVAR1=SVAR1/SVAR2
            CHI1(:,LN)=CHI1(:,LN)-CHI1(:,LN2)*SVAR1
            A1(:,LN)=A1(:,LN)-A1(:,LN2)*SVAR1
          ENDDO
        ENDDO
      ENDDO
!
!     = NOW MAP ONLY TAIL FUNCTIONS BACK AND LEAVE HEAD FUNCTIONS UNTOUCHED        
      DO L=0,LX
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          NOFL=NOFL+1
          IF(NOFL.LE.NORB(L+1)) CYCLE ! LEAVE HEAD FUNCTIONS ALONE
          CHI(:,LN)=CHI1(:,LN)
          A(:,LN)=A1(:,LN)
        ENDDO
      ENDDO
      DEALLOCATE(G)
      DEALLOCATE(CHI1)
      DEALLOCATE(A1)
      DEALLOCATE(AUX)
      DEALLOCATE(AUX1)
!
!     ==========================================================================
!     == NOW INVERT THE MATRIX A: B:=A^{-1}                                   ==
!     == |CHI_I>=SUM_J |PHI_J>A(J,I)                                          ==
!     == 1=|PHI><P|=|PHI>A B <P|=|CHI> ( B<P| )                               ==
!     == <P_CHI_I|=SUM_J B(I,J)<P_PHI_J|                                      ==
!     ==========================================================================
      ALLOCATE(B(LNX,LNX))
      B(:,:)=0.D0
      LX=MAXVAL(LOX)
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
!
        NX=0
        DO LN=1,LNX
          IF(LOX(LN).EQ.L) NX=NX+1
        ENDDO
        IF(NX.EQ.0) CYCLE   

        ALLOCATE(MAT(NX,NX))
        ALLOCATE(MATINV(NX,NX))
!        
!       == MAP A INTO MATRIX MAT ===============================================
        N1=0
        DO LN1=1,LNX
          L1=LOX(LN1)
          IF(L1.NE.L) CYCLE
          N1=N1+1
          N2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            IF(L2.NE.L) CYCLE
            N2=N2+1
            MAT(N2,N1)=A(LN2,LN1)
          ENDDO
        ENDDO
!
!       == INVERT MATRIX =======================================================
        CALL LIB$INVERTR8(NX,MAT,MATINV)
!
!       == MAP MATINV INTO B ===================================================
        N1=0
        DO LN1=1,LNX
          L1=LOX(LN1)
          IF(L1.NE.L) CYCLE
          N1=N1+1
          N2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            IF(L2.NE.L) CYCLE
            N2=N2+1
            B(LN1,LN2)=MATINV(N1,N2) 
          ENDDO
        ENDDO
        DEALLOCATE(MAT)
        DEALLOCATE(MATINV)
      ENDDO
!
!     ==========================================================================
!     === REMOVE TAIL FUNCTIONS                                               ==
!     ==========================================================================
      LNCHI=0
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        NOFL=0
        DO LN=1,LNX
          IF(L.NE.LOX(LN)) CYCLE
          NOFL=NOFL+1
          IF(NOFL.GT.NORB(L+1)) EXIT
          LNCHI=LNCHI+1
          IF(LNCHI.GT.LNXCHI) THEN
            CALL ERROR$MSG('LNCHI OUT OF RANGE')
            CALL ERROR$I4VAL('LNCHI',LNCHI)
            CALL ERROR$I4VAL('LNXCHI',LNXCHI)
            CALL ERROR$STOP('SETUP_CHIFROMPHI')
          END IF
          LOXCHI(LNCHI)=L
          BMAT(LNCHI,:)=B(LN,:)
          AMAT(:,LNCHI)=A(:,LN)
        END DO
      ENDDO
      IF(LNCHI.NE.LNXCHI) THEN
        CALL ERROR$MSG('LNCHI AND LNXCHI ARE INCONSISTENT')
        CALL ERROR$I4VAL('LNCHI',LNCHI)
        CALL ERROR$I4VAL('LNXCHI',LNXCHI)
        CALL ERROR$STOP('SETUP_CHIFROMPHI')
      END IF
!
!     ==========================================================================
!     ==  CLEAN UP                                                            ==
!     ==========================================================================
      DEALLOCATE(CHI)
      DEALLOCATE(A)
      DEALLOCATE(B)
      DEALLOCATE(R)
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$BUILDPARMS()
!     **************************************************************************
!     ** GENERATES AUTOMATICALLY A PAERAMETER FILE FOR THE SETUP CONSTRUCTION **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NFIL=1001
      INTEGER(4),PARAMETER :: NFIL2=1002
      CHARACTER(2)     :: EL
      CHARACTER(12)    :: ID
      CHARACTER(14)    :: SY
      CHARACTER(5)     :: TYPE="'HBS'"
      REAL(8)          :: RMAX=7.D0
      INTEGER(4)       :: I
      REAL(8)          :: AEZ,ZCORE,ZV,RCOV,ANGSTROM
      CHARACTER(1)     :: ATOMTYPE
!     **************************************************************************
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      OPEN(NFIL,FILE='STP.CNTL_HBSAUTOMATIC')
      OPEN(NFIL2,FILE='TEST.STRC_HBSAUTOMATIC')
      rewind nfil
      rewind nfil2
      WRITE(NFIL,FMT='("!SCNTL")')
      DO I=1,105
        CALL PERIODICTABLE$GET(I,'SYMBOL',EL)
        CALL PERIODICTABLE$GET(I,'Z',AEZ)
        CALL PERIODICTABLE$GET(I,'ZCORE',ZCORE)
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
        ZV=AEZ-ZCORE
        ID="'"//TRIM(EL)//'_HBS'//"'"
        SY="'"//EL//"'"
        IF(SY(3:3).EQ.' ')SY(3:3)='_'
        IF(I.GE.1.AND.I.LE.2)   ATOMTYPE='A'   ! ALKALI AND EARTH ALKALI
        IF(I.GE.3.AND.I.LE.4)   ATOMTYPE='A'   ! ALKALI AND EARTH ALKALI
        IF(I.GE.5.AND.I.LE.10)  ATOMTYPE='M'   ! MAIN GROUP
        IF(I.GE.11.AND.I.LE.12) ATOMTYPE='A'   
        IF(I.GE.13.AND.I.LE.18) ATOMTYPE='M'   
        IF(I.GE.19.AND.I.LE.20) ATOMTYPE='A'   
        IF(I.GE.21.AND.I.LE.30) ATOMTYPE='T'   !TRANSITION METAL AND POOR METALS  
        IF(I.GE.31.AND.I.LE.36) ATOMTYPE='M'   
        IF(I.GE.37.AND.I.LE.38) ATOMTYPE='A'   
        IF(I.GE.39.AND.I.LE.48) ATOMTYPE='T'   
        IF(I.GE.49.AND.I.LE.54) ATOMTYPE='M'   
        IF(I.GE.55.AND.I.LE.56) ATOMTYPE='A'   
        IF(I.GE.57.AND.I.LE.71) ATOMTYPE='F'   
        IF(I.GE.72.AND.I.LE.80) ATOMTYPE='T'   
        IF(I.GE.81.AND.I.LE.86) ATOMTYPE='M'   
        IF(I.GE.87.AND.I.LE.88) ATOMTYPE='A'   
        IF(I.GE.89.AND.I.LE.103) ATOMTYPE='F'   
        IF(I.GE.104.AND.I.LE.112) ATOMTYPE='T'   
        IF(I.GE.113.AND.I.LE.118) ATOMTYPE='M'   
!
!       == RMAX MUST BE LARGER THAN THE BOX RADIUS AND SHOULD AT LEAST BE AS ===
!       == LARGE AS THE BOND DISTANCE TO OXYGEN ================================
        RMAX=MAX(2.D0*RCOV+0.5d0,RCOV+0.73d0*ANGSTROM)
!
!       ========================================================================
!       == VALENCE-ONLY SETUPS                                                ==
!       ========================================================================
!       == EXCLUDE D AND F-SHELLS FOR MAINGROUP ELEMENTS (S-P)
        ZV=AEZ-ZCORE
        IF(I.GE.31.AND.I.LE.36) ZV=ZV-10.D0
        IF(I.GE.49.AND.I.LE.54) ZV=ZV-10.D0
        IF(I.GE.72.AND.I.LE.80) ZV=ZV-14.D0
        IF(I.GE.81.AND.I.LE.86) ZV=ZV-24.D0
        IF(I.GE.104) ZV=ZV-14.D0
        WRITE(NFIL,FMT='(T2,"!SETUP  ID=",A," EL=",A," ZV=",f3.0)') &
     &     TRIM(ID),TRIM(SY),ZV
        WRITE(NFIL,FMT='(T10,"TYPE=",A," RBOX/RCOV=2.0 RCSM/RCOV=0.25")') &
           TYPE
        WRITE(NFIL,FMT='(T10,"RCL/RCOV= 0.75 0.75 0.75 0.75")') 
        WRITE(NFIL,FMT='(T10,"LAMBDA= 6. 6. 6. 6.")') 
        WRITE(NFIL,FMT='(T4,"!GRID DMIN=5.E-6 DMAX=0.1 RMAX=",F4.1," !END")')RMAX 
        WRITE(NFIL,FMT='(T4,"!POT  POW=3. VAL=-1.2 RC/RCOV=0.67 !END")') 
        WRITE(NFIL,FMT='(T4,"!CORE POW=3. VAL=0.1  RC/RCOV=0.67 !END")') 
        WRITE(NFIL,FMT='(T2,"!END")') 
        WRITE(NFIL,*) 
!
!       == NOW PREPARE INPUT FOR TEST STRUCTURE FILE
        WRITE(NFIL2,FMT='(T2,"!SPECIES NAME=",A," ID=",A)')TRIM(SY),TRIM(ID)
        IF(I.LE.2) THEN
          WRITE(NFIL2,FMT='(T4,"NPRO=1  LRHOX=2")')
        ELSE IF(I.LE.18) THEN
          WRITE(NFIL2,FMT='(T4,"NPRO=1 1 0  LRHOX=2")')
        ELSE IF(I.LT.56) THEN
          WRITE(NFIL2,FMT='(T4,"NPRO=1 1 1  LRHOX=2")')
        ELSE      
          WRITE(NFIL2,FMT='(T4,"NPRO=1 1 1 1 LRHOX=2")')
        END IF
        IF(ATOMTYPE.EQ.'A') THEN
          WRITE(NFIL2,FMT='(T4,"!HYBRID_X NCORROFL= 1 0 0 0 CV=T HFWEIGHT=0.25 !END")')
        ELSE IF(ATOMTYPE.EQ.'M') THEN
          WRITE(NFIL2,FMT='(T4,"!HYBRID_X NCORROFL= 1 1 0 0 CV=T HFWEIGHT=0.25 !END")')
        ELSE IF(ATOMTYPE.EQ.'T') THEN
          WRITE(NFIL2,FMT='(T4,"!HYBRID_X NCORROFL= 1 0 1 0 CV=T HFWEIGHT=0.25 !END")')
        ELSE IF(ATOMTYPE.EQ.'F') THEN
          WRITE(NFIL2,FMT='(T4,"!HYBRID_X NCORROFL= 1 0 1 1 CV=T HFWEIGHT=0.25 !END")')
        END IF
        WRITE(NFIL2,FMT='(T2,"!END")')
!
!       ========================================================================
!       == SEMI-CORE SETUPS                                                   ==
!       ========================================================================
!       == INCLUDE CORE-S-P SHELL INTO THE VALENCE SO THAT ALL STATES 
!       == WITH A MINIMUM MAIN QUANTUM NUMBER ARE INCLUDED
        IF(AEZ.GE.11) THEN
          ZV=AEZ-ZCORE+8.d0
!       == F-STATES ARE INCLUDED ONLY FOR F-ELEMENTS
!          IF(I.GE.72.AND.I.LE.86) ZV=ZV-14.d0
!          IF(I.GE.104.AND.I.LE.118) ZV=ZV-14.d0
          ID="'"//TRIM(EL)//'_HBS_SC'//"'"
          SY="'"//EL//"'"
          IF(SY(3:3).EQ.' ')SY(3:3)='_'
          WRITE(NFIL,FMT='(T2,"!SETUP  ID=",A," EL=",A," ZV=",f3.0)') &
       &     TRIM(ID),TRIM(SY),ZV
          WRITE(NFIL,FMT='(T10,"TYPE=",A," RBOX/RCOV=1.2 RCSM/RCOV=0.25")') &
               TYPE
          WRITE(NFIL,FMT='(T10,"RCL/RCOV= 0.5 0.5 0.5 0.5")') 
          WRITE(NFIL,FMT='(T10,"LAMBDA= 6. 6. 6. 6.")') 
          WRITE(NFIL,FMT='(T4,"!GRID DMIN=5.E-6 DMAX=0.1 RMAX=",F4.1," !END")')RMAX 
          WRITE(NFIL,FMT='(T4,"!POT  POW=3. VAL=-2.2 RC/RCOV=0.5  !END")') 
          WRITE(NFIL,FMT='(T4,"!CORE POW=3. VAL= 0.1 RC/RCOV=0.5  !END")') 
          WRITE(NFIL,FMT='(T2,"!END")') 
          WRITE(NFIL,*) 
!
!         == NOW PREPARE INPUT FOR TEST STRUCTURE FILE
          SY="'"//EL//"_SC'"
          WRITE(NFIL2,FMT='(T2,"!SPECIES NAME=",A," ID=",A)')TRIM(SY),TRIM(ID)
          IF(I.LE.2) THEN
            WRITE(NFIL2,FMT='(T4,"NPRO=1  LRHOX=2")')
          ELSE IF(I.LE.18) THEN
            WRITE(NFIL2,FMT='(T4,"NPRO=2 2 0  LRHOX=2")')
          ELSE IF(I.LT.56) THEN
            WRITE(NFIL2,FMT='(T4,"NPRO=2 2 1  LRHOX=2")')
          ELSE      
            WRITE(NFIL2,FMT='(T4,"NPRO=2 2 1 1 LRHOX=2")')
          END IF
          IF(ATOMTYPE.EQ.'A') THEN
            WRITE(NFIL2,FMT='(T4,"!HYBRID_X NCORROFL= 2 1 0 0 CV=T HFWEIGHT=0.25 !END")')
          ELSE IF(ATOMTYPE.EQ.'M') THEN
            WRITE(NFIL2,FMT='(T4,"!HYBRID_X NCORROFL= 2 2 0 0 CV=T HFWEIGHT=0.25 !END")')
          ELSE IF(ATOMTYPE.EQ.'T') THEN
            WRITE(NFIL2,FMT='(T4,"!HYBRID_X NCORROFL= 2 1 1 0 CV=T HFWEIGHT=0.25 !END")')
          ELSE IF(ATOMTYPE.EQ.'F') THEN
            WRITE(NFIL2,FMT='(T4,"!HYBRID_X NCORROFL= 2 1 1 1 CV=T HFWEIGHT=0.25 !END")')
          END IF
          WRITE(NFIL2,FMT='(T2,"!END")')
        END IF
      ENDDO
      WRITE(NFIL,FMT='("!END")') 
      WRITE(NFIL,FMT='("!EOB")') 
      WRITE(NFIL,*)
      CLOSE(NFIL)
      RETURN
      END
