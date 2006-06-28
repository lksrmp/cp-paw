!
!.......................................................................
MODULE MIXPULAY_MODULE
  INTEGER(4)           ,SAVE :: NHISTORY=3
  REAL(8)  ,ALLOCATABLE,SAVE :: RHOOLD(:,:)  ! CHARGE DENSITY
  REAL(8)  ,ALLOCATABLE,SAVE :: RESOLD(:,:) ! 1CENTER DENSITY MATRIX
  REAL(8)  ,ALLOCATABLE,SAVE :: RHOOUT(:) ! 1CENTER DENSITY MATRIX
  REAL(8)  ,ALLOCATABLE,SAVE :: A(:,:)
  REAL(8)  ,ALLOCATABLE,SAVE :: AINV(:,:)
  REAL(8)  ,ALLOCATABLE,SAVE :: TEST(:,:)
  REAL(8)  ,ALLOCATABLE,SAVE :: B(:)
  REAL(8)  ,ALLOCATABLE,SAVE :: ALPHABAR(:)
  LOGICAL(4)           ,SAVE :: TFIRST=.TRUE.
  INTEGER(4)           ,SAVE :: IITER
END MODULE MIXPULAY_MODULE

!     ..................................................................
      SUBROUTINE WAVES_MIXRHO_PUL(RHODIM,DENDIM,RHO,DENMAT,TCONV,&
           CONVPSI) ! USING R AND ONLY A
!     ******************************************************************
!     ** ALLPY A PULAY MIXING SCHEME,                                 **
!     ** REF: KRESSE, FURTMUELLER COM. MATER. SCI. 6, 15 (1996)       **
!     **     EQUATIONS 83 TO 87                                       **
!     ******************************************************************
      ! TODO: DEALLOCATE-MECHANISM
      USE MIXPULAY_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: RHODIM
      INTEGER(4) ,INTENT(IN)    :: DENDIM
      REAL(8)    ,INTENT(INOUT) :: RHO(RHODIM) !RHO(NRL,NDIMD)
      real(8)    ,INTENT(INOUT) :: DENMAT(DENDIM) !(LMNXX,LMNXX,NDIMD,NAT)
      LOGICAL(4) ,INTENT(OUT)   :: TCONV
      REAL(8)    ,INTENT(OUT)   :: CONVPSI
      REAL(8)                   :: CONV=1.D-8
      REAL(8)                   :: G=0.5D0 !MIXING FACTOR
      REAL(8)                   :: SVAR,c(1,1)
      REAL(8)                   :: R22,R12,R11
      INTEGER(4)                :: IH1,IH2,NBACK
      INTEGER(4)                :: RHOX
!     ******************************************************************
      !FIRST STEP
      IF(TFIRST) THEN
         IITER=1
         IF(NHISTORY.LT.1) THEN
            CALL ERROR$MSG('NHISTORY MUST BE AT LEAST 1 (FOR LINEAR MIXING)')
            CALL ERROR$STOP('WAVES_MIXRHO')
         END IF
         ! ALLOCATE ARRAYS
         ALLOCATE(RHOOLD(RHODIM+DENDIM,NHISTORY))
         RHOOLD=0.D0
         ALLOCATE(RESOLD(RHODIM+DENDIM,NHISTORY))
         RESOLD=0.D0
         PRINT*,'MIXER, ALLOCATING ',NHISTORY*(RHODIM+DENDIM)*2.D0* &
              8.D0/1024.D0,' KB'
         ALLOCATE(A(NHISTORY,NHISTORY))
         A=0.D0
         ALLOCATE(AINV(NHISTORY,NHISTORY))
         AINV=0.D0
         ALLOCATE(ALPHABAR(NHISTORY))
         ALPHABAR=0.D0
         ! ASSIGN FIRST DENSITY
         RHOOLD(1:RHODIM,1)=RHO(:)
         RHOOLD(RHODIM+1:RHODIM+DENDIM,1)=DENMAT(:)
         TFIRST=.FALSE.
         TCONV=.FALSE.
         CONVPSI=1.D-3 ! 3 AUCH GUT
         RETURN
      END IF
      ! AFTER FIRST STEP: DO MIXING

!!$! DO NOTHING BUT RETURNING THE START DENSITY
!!$RHO(:)=RHOOLD(1:RHODIM,1)
!!$DENMAT(:)=RHOOLD(RHODIM+1:RHODIM+DENDIM,1)
!!$CONVPSI=1.D-6 
!!$RETURN
      ALLOCATE(RHOOUT(RHODIM+DENDIM))
      IITER=IITER+1

      PRINT*,' -------------- MIXER ITER: ',IITER,' ----------'
      CALL LIB$SCALARPRODUCTR8(.TRUE.,RHODIM,1,&
           RHO(1:RHODIM)-RHOOLD(1:RHODIM,1),1,&
           RHO(1:RHODIM)-RHOOLD(1:RHODIM,1),c)
      svar=c(1,1)
      !SVAR=MAXVAL(ABS(RHO(1:RHODIM)-RHOOLD(1:RHODIM,1)))
      CALL MPE$COMBINE('NONE','+',SVAR)
WRITE(*,"('CG-MIXER: RHO   :',3F10.7)") RHO(1:3)
WRITE(*,"('CG-MIXER: RHOOLD:',3F10.7)") RHOOLD(1:3,1)
PRINT*,'CG-MIXER ITER ',IITER,' RHO-DIFF (RES) =',SVAR
      IF(SVAR.GT.1.D-1) THEN
         CONVPSI=1.D-4 ! 4 AUCH GUT
      ELSE
         CONVPSI=1.D-6
      END IF
      IF(SVAR.LT.1.D-4) CONVPSI=1.D-8
      TCONV=.FALSE. !(SVAR.LT.CONV) ! DO NOT USE DESITY CRITERION

      ! RESOLD(:,1) CONTAINS THE NEWEST 
      RESOLD(1:RHODIM,1)=RHO(:)-RHOOLD(1:RHODIM,1)
      RESOLD(RHODIM+1:RHODIM+DENDIM,1)=DENMAT(:) &
           -RHOOLD(RHODIM+1:RHODIM+DENDIM,1)
      !RHOOUT=RHOOLD(:,1)+G*RESOLD(:,1)  ! LINEAR PART OF EQ. 92
      NBACK=MIN(NHISTORY,IITER-1)
PRINT*,'MIXER NBACK', NBACK

DO IH1=1,NHISTORY 
   WRITE(*,"('    MIX RESOLD(',I1,')    ',4F10.5)") IH1,RESOLD(1:4,IH1)
END DO
PRINT*,'  MIX      ======='
DO IH1=1,NHISTORY 
   WRITE(*,"('    MIX RHOOLD(',I1,')    ',4F10.5)") IH1,RHOOLD(1:4,IH1)
END DO
PRINT*,'  MIX      ======='

      RHOX=RHODIM ! USE ONLY DENSITY FOR ADJUSTING FITTING
      !RHOX=RHODIM+DENDIM ! ALSO USE DENAMT FOR FITTING
!!$      DO IH1=2,NBACK ! EQ. 91
!!$         DO IH2=IH1,NBACK
!!$            CALL LIB$SCALARPRODUCTR8(.FALSE.,RHOX,1,&
!!$                 RESOLD(1:RHOX,IH1),1,&
!!$                 RESOLD(1:RHOX,IH2),&
!!$                 A(IH1,IH2))
!!$            A(IH2,IH1)=A(IH1,IH2)
!!$         END DO
!!$      END DO
      ! EQ. 87
      CALL LIB$SCALARPRODUCTR8(.TRUE.,RHOX,NBACK,&
           RESOLD(1:RHOX,1:NBACK),NBACK,RESOLD(1:RHOX,1:NBACK),&
           A(1:NBACK,1:NBACK))

DO IH1=1,NHISTORY 
  WRITE(*,"('MIX A        ',5F10.5)") A(:,IH1)
END DO
      IF(NBACK.GE.1) THEN
         CALL LIB$INVERTR8(NBACK,A(1:NBACK,1:NBACK),AINV(1:NBACK,1:NBACK))
      END IF
DO IH1=1,NHISTORY 
  WRITE(*,"('MIX AINV     ',5F15.5)") AINV(IH1,:)
END DO
!!$CALL LIB$MATMULR8(NBACK-1,NBACK-1,NBACK-1,A,AINV,TEST)
!!$DO IH1=2,NHISTORY 
!!$  WRITE(*,"('MIX AINV*A   ',4F15.5)") TEST(IH1,:)
!!$END DO
      SVAR=0.D0
      DO IH1=1,NBACK ! EQ 87
         DO IH2=1,NBACK
            SVAR=SVAR+AINV(IH1,IH2)
         END DO
      END DO
PRINT*,'MIX SUM AINV:',SVAR
      ALPHABAR=0.D0
      DO IH1=1,NBACK ! EQ 87
         DO IH2=1,NBACK
            ALPHABAR(IH1)=ALPHABAR(IH1)+AINV(IH1,IH2)
         END DO
      END DO
      ALPHABAR=ALPHABAR/SVAR
WRITE(*,"('MIX ALPHABAR ',5F10.5)") ALPHABAR(:)
      RHOOUT=0.D0
      DO IH1=1,NBACK ! EQ 83
         RHOOUT=RHOOUT+ALPHABAR(IH1)*(RHOOLD(:,IH1)+G*RESOLD(:,IH1))
      END DO
WRITE(*,"('CG-MIXER: RHOOUT:',3F10.7)") RHOOUT(1:3)
      ! MAP RHOOLD ON RHO AND DENMAT
      RHO(:)=RHOOUT(1:RHODIM)
      DENMAT(:)=RHOOUT(RHODIM+1:RHODIM+DENDIM)
      ! PROPAGATE RHOOLD AND RESOLD
      DO IH1=NHISTORY,2,-1
         RHOOLD(:,IH1)=RHOOLD(:,IH1-1)
         RESOLD(:,IH1)=RESOLD(:,IH1-1)
      END DO
      RHOOLD(:,1)=RHOOUT
      DEALLOCATE(RHOOUT)
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE WAVES_MIXRHO_PUL_DR(RHODIM,DENDIM,RHO,DENMAT,TCONV,&
           CONVPSI) ! USING DR FIRST LARGE, THEN SMALL
!     ******************************************************************
!     ** ALLPY A PULAY MIXING SCHEME,                                 **
!     ** REF: KRESSE, FURTMUELLER COM. MATER. SCI. 6, 15 (1996)       **
!     ** USING THE DELTA-R SCHEME OF EQUATIONS 88 TO 92               **
!     ******************************************************************
      ! TODO: DEALLOCATE-MECHANISM
      USE MIXPULAY_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: RHODIM
      INTEGER(4) ,INTENT(IN)    :: DENDIM
      REAL(8)    ,INTENT(INOUT) :: RHO(RHODIM) !RHO(NRL,NDIMD)
      REAL(8)    ,INTENT(INOUT) :: DENMAT(DENDIM) !(LMNXX,LMNXX,NDIMD,NAT)
      LOGICAL(4) ,INTENT(OUT)   :: TCONV
      REAL(8)    ,INTENT(OUT)   :: CONVPSI
      REAL(8)                   :: CONV=1.D-6
      REAL(8)                   :: G=0.5D0 !MIXING FACTOR
      REAL(8)                   :: SVAR
      REAL(8)                   :: R22,R12,R11
      INTEGER(4)                :: IH1,IH2,NBACK
      INTEGER(4)                :: RHOX
!     ******************************************************************
      !FIRST STEP
      IF(TFIRST) THEN
         IITER=1
         IF(NHISTORY.LT.1) THEN
            CALL ERROR$MSG('NHISTORY MUST BE AT LEAST 1 (FOR LINEAR MIXING)')
            CALL ERROR$STOP('WAVES_MIXRHO')
         END IF
         ALLOCATE(RHOOLD(RHODIM+DENDIM,NHISTORY))
         ALLOCATE(RESOLD(RHODIM+DENDIM,NHISTORY))
         PRINT*,'MIXER, ALLOCATING ',NHISTORY*(RHODIM+DENDIM)*2.D0* &
              8.D0/1024.D0,' KB'
         ALLOCATE(A(2:NHISTORY,2:NHISTORY))
         A=0.D0
         ALLOCATE(AINV(2:NHISTORY,2:NHISTORY))
         AINV=0.D0
         ALLOCATE(TEST(2:NHISTORY,2:NHISTORY))
         TEST=0.D0
         ALLOCATE(B(2:NHISTORY))
         B=0.D0
         ALLOCATE(ALPHABAR(2:NHISTORY))
         ALPHABAR=0.D0
         RHOOLD=0.D0
         RESOLD=0.D0
         RHOOLD(1:RHODIM,1)=RHO(:)
         RHOOLD(RHODIM+1:RHODIM+DENDIM,1)=DENMAT(:)
         TFIRST=.FALSE.
         TCONV=.FALSE.
         CONVPSI=1.D-3 ! 3 AUCH GUT
         RETURN
      END IF
      ! STEPS AFTER FIRST STEP: DO MIXING
      ALLOCATE(RHOOUT(RHODIM+DENDIM))
      IITER=IITER+1

PRINT*,' -------------- MIXER ITER: ',IITER,' ----------'
SVAR=MAXVAL(ABS(RHO(1:RHODIM)-RHOOLD(1:RHODIM,1)))
WRITE(*,"('CG-MIXER: RHO   :',3F10.7)") RHO(1:3)
WRITE(*,"('CG-MIXER: RHOOLD:',3F10.7)") RHOOLD(1:3,1)
PRINT*,'CG-MIXER ITER ',IITER,' RHO-DIFF (RES) =',SVAR
      IF(SVAR.GT.1.D-1) THEN
         CONVPSI=1.D-4 ! 4 AUCH GUT
      ELSE
         CONVPSI=1.D-6
      END IF
      IF(SVAR.LT.1.D-4) CONVPSI=1.D-8
      TCONV=(SVAR.LT.CONV)

      ! RESOLD(:,1) CONTAINS THE NEWEST 
      RESOLD(1:RHODIM,1)=RHO(:)-RHOOLD(1:RHODIM,1)
      RESOLD(RHODIM+1:RHODIM+DENDIM,1)=DENMAT(:) &
           -RHOOLD(RHODIM+1:RHODIM+DENDIM,1)
      IF(NHISTORY.GE.2) THEN
         RESOLD(:,2)=RESOLD(:,1)-RESOLD(:,2) ! STORE DELTA R IN RESOLD
      END IF
      RHOOUT=RHOOLD(:,1)+G*RESOLD(:,1)  ! LINEAR PART OF EQ. 92
      NBACK=MIN(NHISTORY,IITER-1)
PRINT*,'MIXER NBACK', NBACK

DO IH1=1,NHISTORY 
   WRITE(*,"('    MIX RESOLD(',I1,')    ',4F10.5)") IH1,RESOLD(1:4,IH1)
END DO
PRINT*,'        ======='
DO IH1=1,NHISTORY 
   WRITE(*,"('    MIX RHOOLD(',I1,')    ',4F10.5)") IH1,RHOOLD(1:4,IH1)
END DO
PRINT*,'        ======='

      RHOX=RHODIM ! USE ONLY DENSITY FOR ADJUSTING FITTING
      !RHOX=RHODIM+DENDIM ! ALSO USE DENAMT FOR FITTING
      ! THE FOLLOWING LOOP MAY BE DONE BY ONLY ONE CALL

!!$      DO IH1=2,NBACK ! EQ. 91
!!$         DO IH2=IH1,NBACK
!!$            CALL LIB$SCALARPRODUCTR8(.FALSE.,RHOX,1,&
!!$                 RESOLD(1:RHOX,IH1),1,&
!!$                 RESOLD(1:RHOX,IH2),&
!!$                 A(IH1,IH2))
!!$            A(IH2,IH1)=A(IH1,IH2)
!!$         END DO
!!$      END DO
      CALL LIB$SCALARPRODUCTR8(.TRUE.,RHOX,NBACK-1,&
           RESOLD(1:RHOX,2:NBACK),NBACK-1,RESOLD(1:RHOX,2:NBACK),&
           A(2:NBACK,2:NBACK))

DO IH1=2,NHISTORY ! EQ. 91
  WRITE(*,"('MIX A        ',4F10.5)") A(:,IH1)
END DO
      DO IH1=2,NBACK ! EQ 90 B=<DR_J|R1>
         CALL LIB$SCALARPRODUCTR8(.FALSE.,RHOX,1, &
              RESOLD(1:RHOX,IH1),1,RESOLD(1:RHOX,1),B(IH1))
      END DO
WRITE(*,"('MIX B        ',4F10.5)") B
      IF(NBACK.GE.2) THEN
         CALL LIB$INVERTR8(NBACK-1,A(2:NBACK,2:NBACK),AINV(2:NBACK,2:NBACK))
      END IF
DO IH1=2,NHISTORY 
  WRITE(*,"('MIX AINV     ',4F15.5)") AINV(IH1,:)
END DO
!!$CALL LIB$MATMULR8(NBACK-1,NBACK-1,NBACK-1,A,AINV,TEST)
!!$DO IH1=2,NHISTORY 
!!$  WRITE(*,"('MIX AINV*A   ',4F15.5)") TEST(IH1,:)
!!$END DO
      ALPHABAR=0.D0
      DO IH1=2,NBACK ! EQ 90 
         DO IH2=IH1,NBACK
            ALPHABAR(IH1)=ALPHABAR(IH1)-AINV(IH1,IH2)*B(IH2)
         END DO
      END DO
WRITE(*,"('MIX ALPHABAR ',4F10.5)") ALPHABAR(:)
      DO IH1=2,NBACK ! SUM IN EQ 92  
         RHOOUT=RHOOUT+ALPHABAR(IH1)*(RHOOLD(:,IH1-1)-RHOOLD(:,IH1)+ &
              G*(RESOLD(:,IH1)))
      END DO
      ! MAP RHOOLD ON RHO AND DENMAT
!WRITE(*,"('CG-MIXER: RHOOUT:',3F10.7)") RHOOUT(1:3)
      RHO(:)=RHOOUT(1:RHODIM)
      DENMAT(:)=RHOOUT(RHODIM+1:RHODIM+DENDIM)
      ! PROPAGATE RHOOLD AND RESOLD
      DO IH1=NHISTORY,2,-1
         RHOOLD(:,IH1)=RHOOLD(:,IH1-1)
         RESOLD(:,IH1)=RESOLD(:,IH1-1)
      END DO
      RHOOLD(:,1)=RHOOUT
      DEALLOCATE(RHOOUT)
      RETURN
      END 


