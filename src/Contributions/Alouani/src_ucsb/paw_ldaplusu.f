      MODULE LDAPLUSU_MODULE
      LOGICAL(4)   :: TDO=.FALSE.  
      LOGICAL(4)   :: TFIRST=.TRUE.  
      INTEGER(4)   :: POT_TYPE=3
      INTEGER(4)   :: PRINTOCC=1
      REAL(8)      :: XU=8.D0 
      REAL(8)      :: XJ=0.95D0 
      REAL(8)      :: RC_OFF=-1.0D0 
      REAL(8),DIMENSION(50)                      :: MU_B  
      REAL(8),DIMENSION(50)                      :: CHARGED  
      REAL(8),DIMENSION(-2:2,-2:2)               :: UM 
      REAL(8),DIMENSION(-2:2,-2:2)               :: JM 
      COMPLEX(8),DIMENSION(-2:2,-2:2,-2:2,-2:2)  :: VEMMM 
      END MODULE LDAPLUSU_MODULE
!
!     ..................................................................
      SUBROUTINE LDAPLUSU__SET
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      TDO=.TRUE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LDAPLUSU__SETTING(TDO_)
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(OUT) :: TDO_
      TDO_=TDO
      RETURN
      END
!
!     ...................................................AUTOPI.........
      SUBROUTINE LDAPLUSU(IAT,NRX,LNX,LMNXX,NSPIN,LOX &
     &                ,R1,DEX,AEZ,AEPHI,DENMAT,DETOT,DH)
!     **                                                              **
!     **                                                              **
      USE LDAPLUSU_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      INTERFACE 
        SUBROUTINE READ_PROJO(PRO) 
         IMPLICIT NONE 
         REAL(8),DIMENSION(:,:,:),POINTER :: PRO 
        END SUBROUTINE 
      END INTERFACE 


      INTEGER(4),INTENT(IN)    :: IAT 
      INTEGER(4),INTENT(IN)    :: NRX
      INTEGER(4),INTENT(IN)    :: LMNXX
      INTEGER(4),INTENT(IN)    :: NSPIN
      INTEGER(4),INTENT(IN)    :: LNX
      INTEGER(4),INTENT(IN)    :: LOX(LNX)
      REAL(8)   ,INTENT(IN)    :: R1
      REAL(8)   ,INTENT(IN)    :: DEX
!     INTEGER(4),INTENT(IN)    :: NR
      REAL(8)   ,INTENT(IN)    :: AEZ
      REAL(8)   ,INTENT(IN)    :: AEPHI(NRX,LNX)
      REAL(8)   ,INTENT(IN)    :: DENMAT(LMNXX,LMNXX,NSPIN)
      REAL(8)   ,INTENT(INOUT) :: DH(LMNXX,LMNXX,NSPIN)
      REAL(8)   ,ALLOCATABLE   :: OCC(:,:)
      REAL(8)   ,ALLOCATABLE   :: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE   :: POT(:,:,:)
      COMPLEX(8),ALLOCATABLE   :: POT_MAT(:,:,:)
      REAL(8)                  :: X(LNX,LNX)
      REAL(8)                  :: DH1(LMNXX,LMNXX,NSPIN)
      INTEGER(4)               :: LN
      REAL(8)                  :: RJ(0:3)
      REAL(8)                  :: RU(0:3)
      REAL(8)                  :: WORK(NRX)
!     == POINTER ARRAYS


!***********************************************************************
! BEN MODIF 

      REAL(8)                          :: tolerance 
      INTEGER(4)                       ::IMAX 
      REAL(8),DIMENSION(:,:,:),POINTER :: PRO   
      REAL(8)                          :: F1(NRX)
      REAL(8)                          :: F2(NRX)
      REAL(8),DIMENSION(LNX,LNX)       :: R_MT 
      REAL(8)                          :: HAMSPIN(LMNXX,LMNXX,NSPIN)
      INTEGER(4)                       :: alloc_error

      REAL(8)                          :: n_dwn,n_up,n_tot,Vc_up,Vc_dwn   
      REAL(8)                          :: tempo1,tempo2,tempo3,Vc_moy,xu1    
      REAL(8)                          ::  e_lda,e_lda_u,elda_u     
      REAL(8)                          ::  nb_total,nb_spin(2)      
      COMPLEX(8)                       :: temp1,temp2 


      INTEGER(4),ALLOCATABLE  :: LN_OF_LMN(:)
      INTEGER(4),ALLOCATABLE  :: L_OF_LMN(:)
      INTEGER(4),ALLOCATABLE  :: M_OF_LMN(:)
      INTEGER(4),ALLOCATABLE  :: LM_OF_LMN(:)

      INTEGER(4) :: lmn1,lmn2,ispin1,ispin2,lm1,lm2,lm3,lm4   
      INTEGER(4) :: ilmn1,ilmn2,ilma,ilmb,im1,im2,ima,imb,im3,im4  
      INTEGER(4) :: lmn,lm,m,l 
 
!***********************************************************************

      IF(.NOT.TDO) RETURN
      WRITE(*,*) 
      WRITE(*,*) 'ENTER IN LDA+U SUBROUTINE IAT= ',IAT  


      XEXP=DEXP(DEX)

!     ==================== ALLOCATE ARRAYS =============================

      ALLOCATE(OCC(LMNXX,NSPIN))
      ALLOCATE(RHO(5,5,NSPIN))
      ALLOCATE(POT(LMNXX,LMNXX,NSPIN))
      ALLOCATE(POT_MAT(LMNXX,LMNXX,NSPIN))
      ALLOCATE(LN_OF_LMN(LMNXX))
      ALLOCATE(L_OF_LMN(LMNXX))
      ALLOCATE(M_OF_LMN(LMNXX))
      ALLOCATE(LM_OF_LMN(LMNXX))



!     ========================================================
      LMN=0
      DO LN=1,LNX
         L=LOX(LN)
         DO M=1,2*L+1
          LMN=LMN+1
          LN_OF_LMN(LMN)=LN
          L_OF_LMN(LMN)=L
          M_OF_LMN(LMN)=M
          LM_OF_LMN(LMN)=L*L+M
         ENDDO !(M)
      ENDDO !(LN)

      MAX_LMN=LMNXX 
      MAX_LM=LM_OF_LMN(LMNXX) 


!     ==================================================================
!     ==  OVERLAP                                                     ==
!     ==================================================================

!     ---------------------------------------------------
!     -EVALUATE RASA FOR EACH LM AND FOR THE  ATOM (IAT)- 
!     ---------------------------------------------------


!!    CALL READ_PROJO(PRO) 

!!     R_MT(:,:) = 0. 
!!     tolerance=.01 
!!     IMAX=0 

!!     DO LN1=1,LNX
!!     L1=LOX(LN1)
!!     DO LN2=LN1,LNX
!!        L2=LOX(LN2)
!!        IF (L1.EQ.L2) THEN
!!          F1=PRO(:,LN1,IAT)   
!!          F2=PRO(:,LN2,IAT)   
!!          IR=50 
!!          DO IR=50,NRX  
!!             IF (((ABS(F1(IR)).gt.ABS(F1(IR+1)))          &
!!                .AND. (ABS(F1(IR+1)).gt.ABS(F1(IR+2)))    &
!!                .AND. (ABS(F1(IR+2)).gt.ABS(F1(IR+3)))    &
!!                .AND. (ABS(F1(IR)).lt.tolerance))         &
!!             .OR. ((ABS(F2(IR)).gt.ABS(F2(IR+1)))         &
!!                .AND. (ABS(F2(IR+1)).gt.ABS(F2(IR+2)))    &
!!                .AND. (ABS(F2(IR+2)).gt.ABS(F2(IR+3)))    &
!!                .AND. (ABS(F2(IR)).lt. tolerance))) EXIT
!!          ENDDO 
!!          IMAX=IR  
!!          R_MT(LN1,LN2)=R1 * EXP( (IMAX - 1) * DEX)
!!        ENDIF !(L1=L2)  
!!     ENDDO !(LN2) 
!!     ENDDO !(LN1) 


!----------------
! fix cutoff radius 
!----------------

!!     IF ((L_OF_LMN(MAX_LMN) .GE. 2).AND.    &  
!!         (RC_OFF .GE. 0.D0)) R_MT(:,:)=RC_OFF

       IF (LNX .GT. 2) R_MT(:,:)=2.1D0  
       IF (LNX == 2)   R_MT(:,:)=1.7D0  

!----------------

!!     IF (ASSOCIATED(PRO)) THEN 
!!     DEALLOCATE(PRO,STAT=alloc_error) 
!!     IF (alloc_error /= 0) THEN 
!!        WRITE(6,*) 'ldaplusu: Cannot deallocate PRO' 
!!        STOP 
!!     ENDIF 
!!     ENDIF 

!-------------------------------------------------------------------
!                EVALUATE OVERLAP X(I,J)=<PHI_I|PHI_J> 
!-------------------------------------------------------------------

      WRITE(*,*) 
      WRITE(*,*) 'OVERLAP AEPHI'

      X(:,:)=0.D0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO LN2=LN1,LNX
          L2=LOX(LN2)
          X(LN1,LN2)=0.D0
          X(LN2,LN1)=0.D0
          IF(L1.EQ.L2) THEN
            RI=R1/XEXP
            DO IR=1,NRX
              RI=RI*XEXP
              WORK(IR)=AEPHI(IR,LN1)*AEPHI(IR,LN2)*RI**2
              IF(RI.GT. R_MT(LN1,LN2)) WORK(IR)=0.D0
            ENDDO
            CALL RADIAL__INTEGRAL(R1,DEX,NRX,WORK,X(LN1,LN2))
            X(LN2,LN1)=X(LN1,LN2)
            WRITE(*,*) 'OVER ', LN1,LN2,' = ', X(LN1,LN2) 
          END IF
        ENDDO
      ENDDO
      WRITE(*,*) 

!     ==================================================================
!     ==  CALCULATE OCCUPATIONS                                       ==
!     ==================================================================

!.......................................................................
!    TRANSFORM DENMAT IN SPIN UP AND DOWN FROM TOTAL AND SPIN 
!.......................................................................
     
      HAMSPIN(:,:,:) = 0. 
      HAMSPIN(:,:,1) =0.5*(DENMAT(:,:,1) + DENMAT(:,:,2)) 
      HAMSPIN(:,:,2) =0.5*(DENMAT(:,:,1) - DENMAT(:,:,2))


      MAX_LM=LM_OF_LMN(LMNXX) 

      CALL OCCUPATION(IAT,MAX_LM,NSPIN,LMNXX,LNX,LOX,HAMSPIN,X,OCC,RHO) 

!*****************************
      IF (LOX(LNX) .GE. 2) THEN 
!*****************************

!     ==================================================================
!     == CALCULATE POTENTIAL CORRECTION AND TOTAL ENERGY CONTRIBUTION == 
!     ==================================================================


      DO ISPIN=1,NSPIN
        DO LM=1,LMNXX
          DO LM2=1,LMNXX 
             POT(LM,LM2,ISPIN)=0.D0
             POT_MAT(LM,LM2,ISPIN)=0.D0
          ENDDO 
        ENDDO
      ENDDO

      IF (POT_TYPE == 1) THEN 
         CALL POT_TYPE1(LMNXX,LNX,NSPIN,LM_OF_LMN,   & 
                        LN_OF_LMN,OCC,X,POT_MAT,DETOT) 
      ELSEIF (POT_TYPE == 2) THEN  
         CALL POT_TYPE2(LMNXX,LNX,NSPIN,LM_OF_LMN,M_OF_LMN,   & 
                        L_OF_LMN,LN_OF_LMN,OCC,X,POT_MAT,DETOT) 
      ELSEIF (POT_TYPE == 3) THEN  
         CALL POT_TYPE3(LMNXX,LNX,NSPIN,M_OF_LMN,   & 
                        LN_OF_LMN,OCC,RHO,X,POT_MAT,DETOT) 
      ELSEIF (POT_TYPE == 4) THEN  
         CALL POT_TYPE4(LMNXX,LNX,NSPIN,LM_OF_LMN,M_OF_LMN,   & 
                        L_OF_LMN,LN_OF_LMN,OCC,RHO,X,POT_MAT,DETOT) 
      ENDIF 


!.......................................................................
!    TRANSFORM HAMITONIAN DH IN SPIN UP AND DOWN FROM TOTAL AND SPIN 
!.......................................................................


      HAMSPIN(:,:,1)  = DH(:,:,1) + DH(:,:,2)
      HAMSPIN(:,:,2)  = DH(:,:,1) - DH(:,:,2)

      IF (PRINTOCC /= 0) THEN 
        PRINT*,'   '
        PRINT*,' HAMILTONIAN  MATRIX BEFORE '
        DO ISPIN=1,NSPIN
         DO ILMN1=5,14
          WRITE(*,FMT='(10F10.3)')   & 
                   (27.2D0*HAMSPIN(ILMN1,ILMN2,ISPIN),ILMN2=5,14)
         ENDDO
        PRINT*,'   '
        ENDDO
      ENDIF 

!=================================================================
!======  ADD HUBBARD CORRECTION TO HAMILTONIAN DH(I,J)      ======
!=================================================================

      ELDA_U = 0.D0 
!     DETOT  = 0.D0 

      DO ISPIN=1,NSPIN 
       DO ILMN1=1,MAX_LMN 
        DO ILMN2=1,MAX_LMN 
           HAMSPIN(ILMN1,ILMN2,ISPIN) = HAMSPIN(ILMN1,ILMN2,ISPIN)  & 
                                 + REAL(POT_MAT(ILMN1,ILMN2,ISPIN)) 
!          ELDA_U = ELDA_U + HAMSPIN(LMN1,LMN2,ISPIN)           & 
!                 * REAL(POT_MAT(ILMN1,ILMN2,ISPIN)) 
        ENDDO !(ILMN2)  
       ENDDO !(ILMN1)  
      ENDDO !(ISPIN)  

!     WRITE(*,*) 
!     WRITE(*,*) 'COMPARAISON ENERGY' 
!     WRITE(*,FMT='("DETOT =",F8.5," ELDA_U=",F8.5)') DETOT, ELDA_U 
!     WRITE(*,*) 
!     DETOT = ELDA_U 

      IF (PRINTOCC /= 0) THEN 
        PRINT*,'   '
        PRINT*,' HAMILTONIAN  MATRIX AFTER '
        DO ISPIN=1,NSPIN
         DO ILMN1=5,14
          WRITE(*,FMT='(10F10.3)')   & 
                   (27.2D0*HAMSPIN(ILMN1,ILMN2,ISPIN),ILMN2=5,14)
         ENDDO
        PRINT*,'   '
        ENDDO
      ENDIF 


!.......................................................................
!    TRANSFORM BACK HAMITONIAN DH IN TOTAL AND SPIN FROM SPIN UP AND DWN  
!.......................................................................

      DH(:,:,1)=0.5D0*(HAMSPIN(:,:,1)+HAMSPIN(:,:,2))
      DH(:,:,2)=0.5D0*(HAMSPIN(:,:,1)-HAMSPIN(:,:,2))

!=====================================================================



!************************
      ENDIF !(LNX=2) 
!************************
 
!************************
!************************

      DEALLOCATE(POT)
      DEALLOCATE(POT_MAT)
      DEALLOCATE(RHO)
      DEALLOCATE(OCC)
      DEALLOCATE(LN_OF_LMN)
      DEALLOCATE(L_OF_LMN)
      DEALLOCATE(M_OF_LMN)
      DEALLOCATE(LM_OF_LMN)

      WRITE(*,*) 'GET OUT LDA+U SUBROUTINE ' 
      WRITE(*,*) 
      RETURN
      END

!....................................................
!*****************************************************

      SUBROUTINE COULOMB_EXCHANGE2 
      USE LDAPLUSU_MODULE
      IMPLICIT NONE

        INTERFACE
         SUBROUTINE GAUNT(CGMAT)
          REAL(8),DIMENSION(0:4,-2:2,-2:2) :: CGMAT
         END SUBROUTINE GAUNT
        END INTERFACE



      REAL(8)                          :: F0,F2,F4
      INTEGER,PARAMETER                :: l=2
      INTEGER                          :: m1,m2,m3,m4,k,ma,mb,mc,md
      REAL(8)                          :: U0,U2,U4,J0,J2,J4,coeff
      REAL(8)                          :: Pi,pi_4,Z1,XU1,XJ1  
      COMPLEX(8)                       :: s2,si2
      COMPLEX(8),DIMENSION(-2:2,-2:2)  :: Dm,Am,Bm
      COMPLEX(8)                       :: U00,U22,U44,J00,J22,J44,prefac
      REAL(8),DIMENSION(0:4,-2:2,-2:2) :: CGMAT



          CALL GAUNT(CGMAT)



       Pi=3.141592654
       pi_4=4*Pi
       XU1=XU/27.2D0 
       XJ1=XJ/27.2D0 
       F0=XU1  
       F2=(14/1.625)*XJ1 
       F4=(70./13.)*XJ1  

      Z1=1.D0/SQRT(2.d0)
      S2=(0.70710678D0,0.D0)
      SI2=(0.D0,0.70710678D0)
      DM(:,:)=0.0D0

!=========================================================
!        DEFINE MATRIX TRANSFORMATION FROM CUBIC         =
!               TO SPHERICAL HAROMIC                     =
!                 Y_cub= Dm. Y_sphe                      =
!=========================================================


      Dm(-2,-2)=s2 ;                                Dm(-2,2)=s2
                   Dm(-1,-1)=s2;          Dm(-1,1)=-s2
                               Dm(0,0)=1.
                   Dm(1,-1)=si2 ;         Dm(1,1)=si2
      Dm(2,-2)=si2;                                 Dm(2,2)=-si2

      Am(:,:)=0.d0
      Bm(:,:)=0.d0

      DO m1=-2,2
       DO m2=-2,2
          U0=0.d0 ; U2=0.d0 ; U4=0.d0

          DO ma=-2,2
           DO mb=-2,2
            DO mc=-2,2
             DO md=-2,2
                IF ((ma-mc)/=(md-mb)) CYCLE
                prefac = CONJG(Dm(ma,m1))*CONJG(Dm(mb,m2)) &
                       * Dm(mc,m1) * Dm(md,m2)
                k=0
                U0=U0 + prefac * F0*CGMAT(k,ma,mc)*CGMAT(k,md,mb)

                k=2
                U2=U2 + prefac * F2*CGMAT(k,ma,mc)*CGMAT(k,md,mb)

                k=4
                U4=U4 + prefac * F4*CGMAT(k,ma,mc)*CGMAT(k,md,mb)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          Am(m1,m2)=U0+U2+U4
         ENDDO
      ENDDO


      DO m1=-2,2
       DO m2=-2,2
          J0=0.d0 ; J2=0.d0 ; J4=0.d0

          DO ma=-2,2
           DO mb=-2,2
            DO mc=-2,2
             DO md=-2,2
                IF ((ma-mc)/=(md-mb)) CYCLE
                prefac = CONJG(Dm(ma,m1))*CONJG(Dm(mb,m2)) &
                       * Dm(mc,m2) * Dm(md,m1)
                k=0
                J0=J0 + prefac * F0*CGMAT(k,ma,mc)*CGMAT(k,md,mb)
    
                k=2
                J2=J2 + prefac * F2*CGMAT(k,ma,mc)*CGMAT(k,md,mb)

                k=4
                J4=J4 + prefac * F4*CGMAT(k,ma,mc)*CGMAT(k,md,mb)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          Bm(m1,m2)=J0+J2+J4
         ENDDO
      ENDDO

      UM(:,:)=REAL(Am(:,:))
      JM(:,:)=REAL(Bm(:,:))


! DO m1=-2,2
!       write(6,"('Amm(',i3,';:)=',(3x,5f6.2))")   &
!                 m1,(REAL(Am(m1,m2))*27.2,m2=-2,2)
! ENDDO
!  write(6,*)
!  write(6,*)
!  write(6,*)
! DO m1=-2,2
!       write(6,"('Bmm(',i3,';:)=',(3x,5f6.2))")   &
!                 m1,(REAL(Bm(m1,m2))*27.2,m2=-2,2)
! ENDDO

      END SUBROUTINE COULOMB_EXCHANGE2


!...........................................................

        SUBROUTINE COULOMB_EXCHANGE3
         USE LDAPLUSU_MODULE 
         IMPLICIT NONE

        INTERFACE
         SUBROUTINE GAUNT(CGMAT) 
          IMPLICIT NONE 
          REAL(8),DIMENSION(0:4,-2:2,-2:2) :: CGMAT
         END SUBROUTINE  GAUNT
        END INTERFACE


      REAL(8)                          :: F0,F2,F4
      INTEGER,PARAMETER                :: L=2
      INTEGER                          :: M1,M2,M3,M4,K,MA,MB,MC,MD
      REAL(8)                          :: U0,U2,U4,J0,J2,J4,COEFF
      REAL(8)                          :: PI,PI_4,Z1,XU1,XJ1  
      COMPLEX(8)                       :: S2,SI2
      COMPLEX(8),DIMENSION(-2:2,-2:2)  :: DM,AM,BM
      COMPLEX(8)                       :: U00,U22,U44,J00,J22,J44,PREFAC
      REAL(8),DIMENSION(0:4,-2:2,-2:2) :: CGMAT


          PI=3.141592654
          PI_4=4*PI

          XU1=XU/27.2D0 
          XJ1=XJ/27.2D0 

          F0=XU1
          F2=(14/1.625)*XJ1
          F4=(70./13.)*XJ1

          CALL GAUNT(CGMAT) 


      Z1=1.D0/SQRT(2.d0) 
      S2=(0.70710678D0,0.D0)
      SI2=(0.D0,0.70710678D0)
      DM(:,:)=0.0D0

!=========================================================
!        DEFINE MATRIX TRANSFORMATION FROM CUBIC         =
!               TO SPHERICAL HAROMIC                     =
!                 Y_cub= Dm. Y_sphe                      =
!=========================================================



      DM(-2,-2)=S2 ;                                DM(-2,2)=S2
                   DM(-1,-1)=S2;          DM(-1,1)=-S2
                               DM(0,0)=1.
                   DM(1,-1)=si2 ;         DM(1,1)=SI2
      DM(2,-2)=SI2;                                 DM(2,2)=-SI2

!=========================================================
!  EVALUATE 4-CENTER INTEGRALS <i,j|Vee|k,l>             =
!=========================================================

         VEMMM(:,:,:,:)=(0.d0,0.d0)

        DO M1=-2,2
         DO M3=-2,2
          DO M2=-2,2
           DO M4=-2,2
            U0=0.d0 ; U2=0.d0 ; U4=0.d0

            DO MA=-2,2
             DO MB=-2,2
              DO MC=-2,2
               DO MD=-2,2
                  IF ((MA-MB)/=(MD-MC)) CYCLE
                  PREFAC = CONJG(DM(MA,M1))*CONJG(DM(MC,M3)) &
                         * DM(MB,M2) * DM(MD,M4)
                  K=0
                  U0=U0 + PREFAC * F0*CGMAT(K,MA,MB)*CGMAT(K,MD,MC)
      
                  K=2
                  U2=U2 + PREFAC * F2*CGMAT(K,MA,MB)*CGMAT(K,MD,MC)

                  K=4
                  U4=U4 + PREFAC * F4*CGMAT(K,MA,MB)*CGMAT(K,MD,MC)
               ENDDO !(MA)
              ENDDO !(MB)
             ENDDO !(MC)
            ENDDO !(MD)
            VEMMM(M1,M3,M2,M4)=U0+U2+U4
           ENDDO !(M4)
          ENDDO !(M3)
         ENDDO !(M2)
        ENDDO !(M1)

!-----------------------------------------------


!  DO m1=-2,2
!        write(6,"('Umm(',i3,';:)=',(3x,5f6.2))") m1,   &
!                 (REAL(Um(m1,m2))*27.2,m2=-2,2)
!  ENDDO
!
!  DO m1=-2,2
!        write(6,"('Jmm(',i3,';:)=',(3x,5f6.2))") m1,   &
!                 (REAL(Jm(m1,m2))*27.2,m2=-2,2)
!  ENDDO
 
!...............................

       END SUBROUTINE COULOMB_EXCHANGE3
!=======================================
      SUBROUTINE GAUNT(CGMAT)
      IMPLICIT NONE

      REAL(8),DIMENSION(0:4,-2:2,-2:2)    :: CGMAT


      CGMAT(:,:,:) = 0.D0 


! k=0

      CGMAT(0,+2,+2)=1.
      CGMAT(0,-2,-2)=1.
      CGMAT(0,+1,+1)=1.
      CGMAT(0,-1,-1)=1.
      CGMAT(0,-0,-0)=1.


! k=2

      CGMAT(2,-2,-2)=-.28571428571428571428
      CGMAT(2,-2,-1)=.34992710611188258545
      CGMAT(2,-2,-0)=-.28571428571428571428
      CGMAT(2,-2,+1)=0.
      CGMAT(2,-2,+2)=0.

      CGMAT(2,-1,-2)=-.34992710611188258545
      CGMAT(2,-1,-1)=.14285714285714285714
      CGMAT(2,-1,-0)=.14285714285714285714
      CGMAT(2,-1,+1)=-.34992710611188258545
      CGMAT(2,-1,+2)=0.


      CGMAT(2,-0,-2)=-.28571428571428571428
      CGMAT(2,-0,-1)=-.14285714285714285714
      CGMAT(2,-0,-0)=.28571428571428571428
      CGMAT(2,-0,+1)=-.14285714285714285714
      CGMAT(2,-0,+2)=-.28571428571428571428

      CGMAT(2,+1,-2)=0.
      CGMAT(2,+1,-1)=-.34992710611188258545
      CGMAT(2,+1,-0)=.14285714285714285714
      CGMAT(2,+1,+1)=.14285714285714285714
      CGMAT(2,+1,+2)=-.34992710611188258545


      CGMAT(2,+2,-2)=0.
      CGMAT(2,+2,-1)=0.
      CGMAT(2,+2,-0)=-.28571428571428571428
      CGMAT(2,+2,+1)=.34992710611188258545
      CGMAT(2,+2,+2)=-.28571428571428571428


! k=4

      CGMAT(4,-2,-2)=.04761904761904761904
      CGMAT(4,-2,-1)=-.10647942749998998554
      CGMAT(4,-2,-0)=.18442777839082937548
      CGMAT(4,-2,+1)=-.28171808490950552583
      CGMAT(4,-2,+2)=.39840953644479787998

      CGMAT(4,-1,-2)=.10647942749998998554
      CGMAT(4,-1,-1)=-.19047619047619047619
      CGMAT(4,-1,-0)=.26082026547865053021
      CGMAT(4,-1,+1)=-.30116930096841707923
      CGMAT(4,-1,+2)=.28171808490950552583


      CGMAT(4,-0,-2)=.18442777839082937548
      CGMAT(4,-0,-1)=-.26082026547865053021
      CGMAT(4,-0,-0)=.28571428571428571428
      CGMAT(4,-0,+1)=-.26082026547865053021
      CGMAT(4,-0,+2)=.18442777839082937548

      CGMAT(4,+1,-2)=.28171808490950552583
      CGMAT(4,+1,-1)=-.30116930096841707923
      CGMAT(4,+1,-0)=.26082026547865053021
      CGMAT(4,+1,+1)=-.19047619047619047619
      CGMAT(4,+1,+2)=.10647942749998998554


      CGMAT(4,+2,-2)=.39840953644479787998
      CGMAT(4,+2,-1)=-.28171808490950552583
      CGMAT(4,+2,-0)=.18442777839082937548
      CGMAT(4,+2,+1)=-.10647942749998998554
      CGMAT(4,+2,+2)=.04761904761904761904


      RETURN
      END SUBROUTINE


!======================================
      SUBROUTINE READ_PROJO(PRO) 
       IMPLICIT NONE 
       REAL(8),DIMENSION(:,:,:),POINTER :: pro 

       INTEGER(4)                       :: NRX 
       INTEGER(4)                       :: NAT  
       INTEGER(4)                       :: NSP  
       INTEGER(4)                       :: LNXX  
       INTEGER(4)                       :: IR  
       INTEGER(4)                       :: LN  
       INTEGER(4)                       :: IAT   
       INTEGER(4)                       :: ISP   
       INTEGER(4)                       :: alloc_error    
       INTEGER(4),ALLOCATABLE           :: LNX(:) 
       INTEGER(4),ALLOCATABLE           :: ISPECIES(:) 
       REAL(8),ALLOCATABLE              :: tab(:) 
       NULLIFY(PRO) 

       OPEN(UNIT=402,FILE='PROJO',FORM='FORMATTED',STATUS='OLD')

      READ(402,*) NRX
      READ(402,*) NAT
      READ(402,*) NSP
      READ(402,*) LNXX 



      ALLOCATE(LNX(NSP),STAT=ALLOC_ERROR) 
      IF (ALLOC_ERROR /= 0) THEN
      WRITE(6,*) 'READ_PROJO: CANNOT ALLOCATE LNX'
      STOP
      ENDIF

      ALLOCATE(ISPECIES(NAT),STAT=ALLOC_ERROR) 
      IF (ALLOC_ERROR /= 0) THEN
         WRITE(6,*) 'READ_PROJO: CANNOT ALLOCATE ISPECIES'
         STOP
      ENDIF

      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES) 

      READ(402,*) (LNX(ISP),ISP=1,NSP)

     ALLOCATE(PRO(NRX,LNXX,NAT),STAT=ALLOC_ERROR) 
      IF (ALLOC_ERROR /= 0) THEN 
         WRITE(6,*) 'READ_PROJO: CANNOT ALLOCATE PROJO'
         STOP 
      ENDIF 
      PRO(:,:,:)=0. 


       DO IAT=1,NAT
          ISP=ISPECIES(IAT) 

          DO IR=1,NRX
             READ(402,*)(PRO(IR,LN,IAT),LN=1,LNX(ISP))
          ENDDO
       ENDDO
       CLOSE(402)
    
       DEALLOCATE(LNX,STAT=ALLOC_ERROR) 
       IF (ALLOC_ERROR /= 0) THEN
          WRITE(6,*) 'READ_PROJO: cANNOT DEALLOCATE LNX'
          STOP
       ENDIF
       DEALLOCATE(ISPECIES,STAT=ALLOC_ERROR)  
       IF (ALLOC_ERROR /= 0) THEN
          WRITE(6,*) 'READ_PROJO: cANNOT DEALLOCATE LNX'
          STOP
       ENDIF

       RETURN
      END SUBROUTINE READ_PROJO  

!     ..................................................................
      SUBROUTINE LDAPLUSU$REPORT(NFIL)
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IPOT
      CHARACTER(24)         :: TYPE 


      
      IF(.NOT.TDO) RETURN
      WRITE(NFIL,*)
      WRITE(NFIL,FMT='("LDA+U REPORT"/12("="))')
      WRITE(NFIL,*) 
      WRITE(NFIL,FMT='("HUBBARD CORRECTION APPLY ON NI-D ORBITALS")') 
!     IF (POT_TYPE==3) THEN 
!      TYPE=" ROTATIONNALLY INVARIANT"
!      WRITE(NFIL,FMT='("LDA+U POTENTIAL TYPE:",34("."),A24)') TYPE  
!     ENDIF  
      WRITE(NFIL,FMT='(" POT_TYPE = ",45("."),I3)') POT_TYPE 
      WRITE(NFIL,FMT='(" U        = ",45("."),F5.1,"  [EV]")') XU 
      WRITE(NFIL,FMT='(" J        = ",45("."),F6.2," [EV]")') XJ 
      WRITE(NFIL,FMT='(" RC_OFF   = ",45("."),F5.1," [A.U]")') RC_OFF 
      WRITE(NFIL,*)

      RETURN
      END


!     ..................................................................
     SUBROUTINE OCCUPATION(IAT,MAX_LM,NSPIN,LMNXX,LNX,LOX,HAMSPIN,X,OCC,RHO)  
      USE LDAPLUSU_MODULE  
      IMPLICIT NONE 
      INTEGER(4),INTENT(IN)    :: IAT  
      INTEGER(4),INTENT(IN)    :: MAX_LM 
      INTEGER(4),INTENT(IN)    :: NSPIN 
      INTEGER(4),INTENT(IN)    :: LMNXX 
      INTEGER(4),INTENT(IN)    :: LNX
      INTEGER(4),INTENT(IN)    :: LOX(LNX)
      REAL(8)                  :: HAMSPIN(LMNXX,LMNXX,NSPIN)
      REAL(8)                  :: X(LNX,LNX)
      REAL(8)                  :: OCC(LMNXX,NSPIN)
      REAL(8)                  :: RHO(5,5,NSPIN)

      INTEGER(4) :: ISPIN,LM,LM2,LM1,L1,LMN1,M1,M2,L2,LMN2,LN1,LN2   
      REAL(8)    :: NB_SPIN(2),NB_TOTAL 

      DO ISPIN=1,NSPIN
        DO LM=1,LMNXX
          OCC(LM,ISPIN)=0.D0
        ENDDO
      ENDDO

      DO LM1=1,5 
       DO LM2=1,5 
          RHO(LM1,LM2,1)=0.D0 
          RHO(LM1,LM2,2)=0.D0 
       ENDDO 
      ENDDO 
 
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO M1=1,2*L1+1
          LMN1=LMN1+1
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO M2=1,2*L2+1
              LMN2=LMN2+1
   
              IF ((L1 == 2).AND.(L2 == 2)) THEN 
               DO ISPIN=1,NSPIN  
                 RHO(M1,M2,ISPIN)=RHO(M1,M2,ISPIN) & 
                    +HAMSPIN(LMN1,LMN2,ISPIN)*X(LN1,LN2)  
               ENDDO 
              ENDIF 

              IF(L1.EQ.L2.AND.M1.EQ.M2) THEN
                LM=L1**2+M1
                DO ISPIN=1,NSPIN
                  OCC(LM,ISPIN)=OCC(LM,ISPIN) &
     &             +HAMSPIN(LMN1,LMN2,ISPIN)*X(LN1,LN2)
                ENDDO
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO




!------------------------------------------------------------
!== EVALUATE  CHARGE AND MAGNETIC MOMENT                =====
!------------------------------------------------------------

      NB_SPIN(1)=0.D0 
      NB_SPIN(2)=0.D0 

      DO ISPIN=1,NSPIN  
         DO LM1=5,MAX_LM 
            NB_SPIN(ISPIN) = NB_SPIN(ISPIN) & 
                           + OCC(LM1,ISPIN) 
         ENDDO
      ENDDO 

      NB_TOTAL=NB_SPIN(1)+NB_SPIN(2) 
      
      IF (LNX ==2) THEN 
         NB_TOTAL=0.D0 
         DO LM1=2,MAX_LM 
            NB_TOTAL=NB_TOTAL+OCC(LM1,1)+OCC(LM1,2) 
         ENDDO 
      ENDIF 

      MU_B(IAT)   = NB_SPIN(1)-NB_SPIN(2)
      CHARGED(IAT) = NB_TOTAL 

!-----------------------------------------------------------

  IF (PRINTOCC /= 0) THEN 
      PRINT*,'OCC ',NSPIN,LMNXX
      DO ISPIN=1,NSPIN
        WRITE(*,FMT='(9F8.5)')(OCC(LM,ISPIN),LM=1,MAX_LM)
      ENDDO

      DO ISPIN=1,NSPIN
        WRITE(*,FMT='("NB_SPIN(",I1,")= ",F6.3)') ISPIN,NB_SPIN(ISPIN)
      ENDDO
      WRITE(*,FMT='("MAGNETIC MOMENT= ",F6.3)') MU_B(IAT) 
      WRITE(*,FMT='("TOTAL    CHARGE= ",F6.3)') CHARGED(IAT) 

      PRINT*,'   ' 
      PRINT*,' DENSITY MATRIX ' 
      DO ISPIN=1,NSPIN
       DO LM1=1,5 
        WRITE(*,FMT='(9F8.5)')(RHO(LM1,LM2,ISPIN),LM2=1,5)
       ENDDO 
      PRINT*,'   ' 
      ENDDO
  ENDIF 


      END SUBROUTINE OCCUPATION 


!     ..................................................................

     SUBROUTINE POT_TYPE1(LMNXX,LNX,NSPIN,LM_OF_LMN,LN_OF_LMN, & 
                                            OCC,X,POT_MAT,DETOT)
      USE LDAPLUSU_MODULE  
      IMPLICIT NONE 

      INTEGER(4)          :: LMNXX 
      INTEGER(4)          :: LNX  
      INTEGER(4)          :: NSPIN  
      INTEGER(4)          :: LM_OF_LMN(LMNXX)
      INTEGER(4)          :: LN_OF_LMN(LMNXX)
      REAL(8)             :: OCC(LMNXX,NSPIN)
      REAL(8)             :: X(LNX,LNX)
      COMPLEX(8)          :: POT_MAT(LMNXX,LMNXX,NSPIN)
      REAL(8),INTENT(OUT) :: DETOT 

      INTEGER(4) :: LMN1,LN1,LM1,ISPIN  
      REAL(8)    :: XU1,TEMP1,TEMP2  

      XU1=XU/27.2
      POT_MAT(:,:,:)=0.D0

      DO ISPIN=1,NSPIN 
       DO LMN1=1,LMNXX 
         LM1=LM_OF_LMN(LMN1) 
         LN1=LN_OF_LMN(LMN1) 
         POT_MAT(LMN1,LMN1,ISPIN)=(XU1*(0.5D0-OCC(LM1,ISPIN)))*X(LN1,LN1)  
       ENDDO !(LMN1) 
      ENDDO !(ISPIN)  

     END SUBROUTINE POT_TYPE1


!     ..................................................................

     SUBROUTINE POT_TYPE2(LMNXX,LNX,NSPIN,LM_OF_LMN,M_OF_LMN,    & 
                           L_OF_LMN,LN_OF_LMN,OCC,X,POT_MAT,DETOT)
      USE LDAPLUSU_MODULE  
      IMPLICIT NONE 

      INTEGER(4)          :: LMNXX 
      INTEGER(4)          :: LNX  
      INTEGER(4)          :: NSPIN  
      INTEGER(4)          :: LM_OF_LMN(LMNXX)
      INTEGER(4)          :: M_OF_LMN(LMNXX)
      INTEGER(4)          :: L_OF_LMN(LMNXX)
      INTEGER(4)          :: LN_OF_LMN(LMNXX)
      REAL(8)             :: OCC(LMNXX,NSPIN)
      REAL(8)             :: X(LNX,LNX)
      COMPLEX(8)          :: POT_MAT(LMNXX,LMNXX,NSPIN)
      REAL(8),INTENT(OUT) :: DETOT 

      INTEGER(4) :: LMN1,LN1,LM1,ISPIN,MAX_LM,MAX_LMN,ISPIN2    
      INTEGER(4) :: ILMN1,ILMN2,IM1,IM2,ILN1,ILM1,ILM2,L1,M1 
      REAL(8)    :: XU1,XJ1,UEFF,TEMPO1,TEMPO2,TEMPO3,OVER,V1    
      REAL(8)    :: NB_SPIN(2),NB_TOTAL 
      XU1=XU/27.2
      XJ1=XJ/27.2
      UEFF=XU1-XJ1/2.D0 
      POT_MAT(:,:,:)=0.D0
      MAX_LMN=LMNXX 
      MAX_LM=LM_OF_LMN(LMNXX) 

!------------------------------------
!== EVALUATE UM AND JM MATRICES    ==  
!------------------------------------

   IF (TFIRST .EQV. .TRUE.) THEN
      TFIRST=.FALSE.
      CALL COULOMB_EXCHANGE2 
   ENDIF

!--------------------------------------------


!       DO ISPIN=1,2
!          ISPIN2=1
!          IF (ISPIN == 1) ISPIN2=2

!          DO ILMN1=1,MAX_LMN
!             ILM1=LM_OF_LMN(ILMN1)
!             IF (L_OF_LMN(ILMN1) == 2) THEN
!                TEMPO1=0
!                TEMPO2=0
!                DO ILMN2=1,MAX_LM
!                   IF (L_OF_LMN(ILMN2) == 2) THEN
!                      IM1=M_OF_LMN(ILMN1)-3
!                      IM2=M_OF_LMN(ILMN2)-3

!                      TEMPO1=TEMPO1+(UM(IM1,IM2)-UEFF)*   &
!                         OCC(ILMN2,ISPIN2)

!                      IF (IM1 /= IM2) THEN
!                         TEMPO2=TEMPO2+(UM(IM1,IM2)-JM(IM1,IM2)-UEFF)*   &
!                           OCC(ILMN2,ISPIN)
!                      ENDIF !(M1/=M2)
!                   ENDIF !(L2=2)
!                ENDDO !(ILMN2)
!                TEMPO3=UEFF*(1/2.-OCC(ILM1,ISPIN))-XJ1/4.
!                ILN1=LN_OF_LMN(ILMN1)
!                OVER=X(ILN1,ILN1)
!                V1=TEMPO1+TEMPO2+TEMPO3
!                POT_MAT(ILMN1,ILMN1,ISPIN)=V1*OVER
!             ENDIF !(L1=2)
!          ENDDO !(ILMN1)
!       ENDDO !(ISPIN)

! PRINT*,'POT_TYPE2   '
!     PRINT*,' POTENTIAL MATRIX '
!     DO ISPIN=1,NSPIN
!      DO ILMN1=5,14
!       WRITE(*,FMT='(10F8.5)')   & 
!                (27.2D0*REAL(POT_MAT(ILMN1,ILMN2,ISPIN)),ILMN2=5,14)
!      ENDDO
!     PRINT*,'   '
!     ENDDO

!---------------------------------------
! NEW VERSION 
!---------------------------------------


      NB_SPIN(1)=0.D0 
      NB_SPIN(2)=0.D0 

      L1=2 
      DO ISPIN=1,NSPIN  
         DO M1=1,2*L1+1
            LM1=L1**2+M1 
            NB_SPIN(ISPIN) = NB_SPIN(ISPIN) & 
                           + OCC(LM1,ISPIN) 
         ENDDO
      ENDDO 

      NB_TOTAL=NB_SPIN(1)+NB_SPIN(2) 



        DO ISPIN=1,2
           ISPIN2=1
           IF (ISPIN == 1) ISPIN2=2

           DO ILMN1=1,MAX_LMN
              ILM1=LM_OF_LMN(ILMN1)
              IF (L_OF_LMN(ILMN1) == 2) THEN
                 TEMPO1=0
                 TEMPO2=0
                 DO ILMN2=1,MAX_LM
                    IF (L_OF_LMN(ILMN2) == 2) THEN
                       IM1=M_OF_LMN(ILMN1)-3
                       IM2=M_OF_LMN(ILMN2)-3

                       TEMPO1=TEMPO1+(UM(IM1,IM2)+UM(IM2,IM1))*   &
                          OCC(ILMN2,ISPIN2)

                       IF (IM1 /= IM2) THEN
                          TEMPO2=TEMPO2+ OCC(ILMN2,ISPIN) *  & 
                            (UM(IM1,IM2)-JM(IM1,IM2)+UM(IM2,IM1)-JM(IM2,IM1))   
                       ENDIF !(M1/=M2)
                    ENDIF !(L2=2)
                 ENDDO !(ILMN2)
                 TEMPO3= - XU1 * (NB_TOTAL-0.5D0)          & 
                         + XJ1 * (NB_SPIN(ISPIN)-0.5D0)
                 ILN1=LN_OF_LMN(ILMN1)
                 OVER=X(ILN1,ILN1)
                 V1=0.5D0 * (TEMPO1 + TEMPO2) + TEMPO3
                 POT_MAT(ILMN1,ILMN1,ISPIN)=V1*OVER
              ENDIF !(L1=2)
           ENDDO !(ILMN1)
        ENDDO !(ISPIN)

  IF (PRINTOCC /= 0) THEN 
  PRINT*,'   '
      PRINT*,' POTENTIAL MATRIX '
      DO ISPIN=1,NSPIN
       DO ILMN1=5,14
        WRITE(*,FMT='(10F8.5)')   & 
                 (27.2D0*REAL(POT_MAT(ILMN1,ILMN2,ISPIN)),ILMN2=5,14)
       ENDDO
      PRINT*,'   '
      ENDDO
  ENDIF 


!================================
!== EVALUATE TOTAL ENERGY      ==
!================================



   DETOT= 0.D0
   TEMPO1= 0.D0 
   TEMPO2= 0.D0 

  DO ISPIN=1,NSPIN 
     ISPIN2=1
     IF (ISPIN == 1) ISPIN2=2

   DO ILM1=1,MAX_LM 
   DO ILM2=1,MAX_LM 
     IF (L_OF_LMN(ILM1) /= 2) CYCLE 
     IF (L_OF_LMN(ILM2) /= 2) CYCLE 
       IM1=M_OF_LMN(ILM1)-3
       IM2=M_OF_LMN(ILM2)-3
 
      TEMPO1=TEMPO1 + UM(IM1,IM2) * OCC(ILM1,ISPIN)   & 
                  * OCC(ILM2,ISPIN2)  

      IF (IM1 /= IM2) THEN 
         TEMPO2= TEMPO2 + (UM(IM1,IM2)-JM(IM1,IM2))      & 
                      * OCC(ILM1,ISPIN)*OCC(ILM2,ISPIN) 
      ENDIF 
   ENDDO !(ILM2) 
   ENDDO !(ILM1) 
  ENDDO !(ISPIN) 

   TEMPO3 = XU1 * NB_TOTAL*(NB_TOTAL-1)           &
          - XJ1 * (NB_SPIN(1)*(NB_SPIN(1)-1)      & 
                    +  NB_SPIN(2)*(NB_SPIN(2)-1))
  
   DETOT=0.5D0*(TEMPO1 + TEMPO2 - TEMPO3)  

     END SUBROUTINE POT_TYPE2



!     ..................................................................
      SUBROUTINE POT_TYPE3(LMNXX,LNX,NSPIN,M_OF_LMN,LN_OF_LMN,  & 
                                        OCC,RHO,X,POT_MAT,DETOT)
      USE LDAPLUSU_MODULE  
      IMPLICIT NONE 

      INTEGER(4),INTENT(IN)    :: LMNXX 
      INTEGER(4),INTENT(IN)    :: LNX  
      INTEGER(4),INTENT(IN)    :: NSPIN  
      INTEGER(4),INTENT(IN)    :: M_OF_LMN(LMNXX)
      INTEGER(4),INTENT(IN)    :: LN_OF_LMN(LMNXX)
      REAL(8),INTENT(IN)       :: OCC(LMNXX,NSPIN)
!     REAL(8),INTENT(IN)       :: RHO(5,5,NSPIN)
      REAL(8),INTENT(INOUT)    :: RHO(5,5,NSPIN)
      REAL(8),INTENT(IN)       :: X(LNX,LNX)
      COMPLEX(8),INTENT(OUT)   :: POT_MAT(LMNXX,LMNXX,NSPIN)
      REAL(8),INTENT(OUT)      :: DETOT 


      INTERFACE 
       SUBROUTINE COULOMB_EXCHANGE3
        USE LDAPLUSU_MODULE 
        IMPLICIT NONE
       END SUBROUTINE 
      END INTERFACE 

      INTEGER(4) :: ILMN1,MAX_LMN,ILMN2,ILMA,ILMB,MAX_LM 
      INTEGER(4) :: IM1,IM2,IMA,IMB,ISPIN,L1,M1,LM1,ISPIN2  
      INTEGER(4) :: ILN1,ILN2,LM3,LM4   
      INTEGER(4) :: IM3,IM4,LM2,ISPIN1,LMAX  

      REAL(8)    :: XU1,XJ1,OVER,V1,TEMPO1,TEMPO2,TEMPO3  
      REAL(8)    :: NB_TOTAL,NB_SPIN(2) 
      COMPLEX(8) :: TEMP1,TEMP2,TEMP3  

      
      XU1=XU/27.2
      XJ1=XJ/27.2

!     TEMP
      WRITE(*,*)"XU IS EQUAL TO", XU
      WRITE(*,*)"XJ IS EQUAL TO", XJ
!     END TEMP


      MAX_LMN=LMNXX 
      MAX_LM=9 

      POT_MAT(:,:,:)=(0.D0,0.D0) 

!------------------------------------
!== EVALUATE FOUR-CENTER INTEGRALS ==  
!------------------------------------
 
      IF (TFIRST .EQV. .TRUE.) THEN 
         TFIRST=.FALSE. 
         CALL COULOMB_EXCHANGE3
      ENDIF 



!------------------------------------------------------------
!== EVALUATE NB OF ELECTRONS IN FOR BOTH SPIN FOR D-ATOMS  ==  
!------------------------------------------------------------

      NB_SPIN(1)=0.D0 
      NB_SPIN(2)=0.D0 

      L1=2 
      DO ISPIN=1,NSPIN  
         DO M1=1,2*L1+1
            LM1=L1**2+M1 
            NB_SPIN(ISPIN) = NB_SPIN(ISPIN) & 
                           + OCC(LM1,ISPIN) 
         ENDDO
      ENDDO 

      NB_TOTAL=NB_SPIN(1)+NB_SPIN(2) 

!---------------------------------------------------
!  EVALUATION OF POTENTIAL 
!---------------------------------------------------


      DO ISPIN=1,NSPIN 
       DO ILMN1=5,MAX_LMN
       DO ILMN2=5,MAX_LMN

           TEMP1=0.D0
           TEMP2=0.D0
           TEMP3=0.d0

          DO ISPIN2=1,NSPIN 
           DO ILMA=5,MAX_LM
           DO ILMB=5,MAX_LM


              IM1=M_OF_LMN(ILMN1)-3
              IM2=M_OF_LMN(ILMN2)-3
              IMA=M_OF_LMN(ILMA)-3
              IMB=M_OF_LMN(ILMB)-3

              TEMP1 = TEMP1 + VEMMM(IM1,IMA,IM2,IMB) *   &
                      RHO(ILMA-4,ILMB-4,ISPIN2)

              IF (ISPIN == ISPIN2) THEN 
                 TEMP1 = TEMP1 - VEMMM(IM1,IMA,IMB,IM2) *  & 
                       RHO(ILMA-4,ILMB-4,ISPIN2) 
              ENDIF 
           ENDDO !(ILMB)
           ENDDO !(ILMA)
          ENDDO !(ISPIN2)
  
        IF (M_OF_LMN(ILMN1) == M_OF_LMN(ILMN2)) THEN
           TEMP3 = - XU1 *(NB_TOTAL-1.D0/2.D0)    &
                   + XJ1 *(NB_SPIN(ISPIN)-1.D0/2.D0)
        ENDIF

        ILN1 = LN_OF_LMN(ILMN1)
        ILN2 = LN_OF_LMN(ILMN2)
        OVER = X(ILN1,ILN2)

        POT_MAT(ILMN1,ILMN2,ISPIN)=(TEMP1 + TEMP3) * OVER

       ENDDO !(ILMN2)
       ENDDO !(ILMN1)
      ENDDO !(ISPIN)

      IF (PRINTOCC /= 0) THEN 
        PRINT*,' POTENTIAL MATRIX '
        DO ISPIN=1,NSPIN
         DO ILMN1=5,14
          WRITE(*,FMT='(10F8.5)')   & 
               (27.2D0*REAL(POT_MAT(ILMN1,ILMN2,ISPIN)),ILMN2=5,14)
         ENDDO
        PRINT*,'   '
        ENDDO
      ENDIF  
!----------------------------------------------------------

!================================
!== EVALUATE TOTAL ENERGY      ==  
!================================

      DETOT= 0.D0 
      TEMP1=(0.D0,0.D0) 
      TEMP2=(0.D0,0.D0) 

      DO ISPIN1=1,NSPIN 
      DO ISPIN2=1,NSPIN 

      DO LM1=5,MAX_LM 
      DO LM2=5,MAX_LM 
      DO LM3=5,MAX_LM 
      DO LM4=5,MAX_LM 

         IM1=M_OF_LMN(LM1)-3  
         IM2=M_OF_LMN(LM2)-3  
         IM3=M_OF_LMN(LM3)-3  
         IM4=M_OF_LMN(LM4)-3  

         TEMP1 = TEMP1 + VEMMM(IM1,IM3,IM2,IM4)     & 
               * RHO(LM1-4,LM2-4,ISPIN1)            & 
               * RHO(LM3-4,LM4-4,ISPIN2) 

         IF (ISPIN1 == ISPIN2) THEN 
            TEMP1 = TEMP1 - (VEMMM(IM1,IM3,IM4,IM2)    & 
                  * RHO(LM1-4,LM2-4,ISPIN1)            & 
                  * RHO(LM3-4,LM4-4,ISPIN2))  
         ENDIF 
  
      ENDDO !(LM1)
      ENDDO !(LM2)
      ENDDO !(LM3)
      ENDDO !(LM4)

      ENDDO !(ISPIN2)  
      ENDDO !(ISPIN1)  

      TEMP1= 0.5D0 * TEMP1 


      TEMP2 = 0.5D0 * XU1 *  NB_TOTAL *(NB_TOTAL-1)       & 
            - 0.5D0 * XJ1 * (NB_SPIN(1)*(NB_SPIN(1)-1) +  &
                          NB_SPIN(2)*(NB_SPIN(2)-1))


      DETOT = REAL(TEMP1  - TEMP2) 


     END SUBROUTINE POT_TYPE3 
!     ..................................................................
!     ..................................................................

     SUBROUTINE POT_TYPE4(LMNXX,LNX,NSPIN,LM_OF_LMN,M_OF_LMN,    & 
                           L_OF_LMN,LN_OF_LMN,OCC,RHO,X,POT_MAT,DETOT)
      USE LDAPLUSU_MODULE  
      IMPLICIT NONE 

      INTEGER(4)          :: LMNXX 
      INTEGER(4)          :: LNX  
      INTEGER(4)          :: NSPIN  
      INTEGER(4)          :: LM_OF_LMN(LMNXX)
      INTEGER(4)          :: M_OF_LMN(LMNXX)
      INTEGER(4)          :: L_OF_LMN(LMNXX)
      INTEGER(4)          :: LN_OF_LMN(LMNXX)
      REAL(8)             :: OCC(LMNXX,NSPIN)
      REAL(8)             :: RHO(5,5,NSPIN)
      REAL(8)             :: X(LNX,LNX)
      COMPLEX(8)          :: ZC1 
      COMPLEX(8)          :: POT_MAT(LMNXX,LMNXX,NSPIN)
      REAL(8),INTENT(OUT) :: DETOT 

      INTEGER(4) :: LMN1,LN1,LM1,ISPIN,MAX_LM,MAX_LMN,ISPIN2    
      INTEGER(4) :: ILMN1,ILMN2,IM1,IM2,ILN1,ILM1,ILM2,L1,M1,ILN2 
      REAL(8)    :: XU1,XJ1,UEFF,TEMPO1,TEMPO2,TEMPO3,OVER,V1    
      REAL(8)    :: NB_SPIN(2),NB_TOTAL 


!-------------------------
!== INITIALIZE OBJECTS  ==
!-------------------------

      XU1=XU/27.2D0
      XJ1=XJ/27.2D0
      ZC1 = (1.D0,0.D0) 
      POT_MAT(:,:,:)= (0.D0,0.D0) 
      MAX_LMN=LMNXX 
      MAX_LM=LM_OF_LMN(LMNXX) 


!=======================================================================
!             EVALUATION OF LDA+U POTENTIAL
!=======================================================================

        DO ISPIN=1,2

         DO ILMN1=1,MAX_LMN
            IF (L_OF_LMN(ILMN1) /= 2) CYCLE 
            IM1 = M_OF_LMN(ILMN1) 

              DO ILMN2=1,MAX_LMN
                IF (L_OF_LMN(ILMN2) /= 2) CYCLE 
                IM2 = M_OF_LMN(ILMN2) 

                  TEMPO1 = -1.D0 * (XU1 - XJ1) * RHO(IM1,IM2,ISPIN) 

                  IF (IM1 == IM2) THEN
                     TEMPO1 = TEMPO1 + (XU1 - XJ1) *0.5D0 
                  ENDIF !(M1/=M2)

                 ILN1 = LN_OF_LMN(ILMN1)
                 ILN2 = LN_OF_LMN(ILMN2)
                 OVER=X(ILN1,ILN2)

                 POT_MAT(ILMN1,ILMN2,ISPIN)= ZC1 * TEMPO1 * OVER
              ENDDO !(ILMN2)

           ENDDO !(ILMN1)
        ENDDO !(ISPIN)

!=======================================================================
!                         TEST PRINTOUT
!=======================================================================

      IF (PRINTOCC /= 0) THEN 
        PRINT*,'   '
        PRINT*,' POTENTIAL MATRIX '
        DO ISPIN=1,NSPIN
         DO ILMN1=5,14
          WRITE(*,FMT='(10F8.5)')   & 
                   (27.2D0*REAL(POT_MAT(ILMN1,ILMN2,ISPIN)),ILMN2=5,14)
         ENDDO
        PRINT*,'   '
        ENDDO
      ENDIF 

!=======================================================================
! TOTAL ENERGY LDA+U DOUBLE COUNTING CORRECTION
! EQ.5, DUDAREV ET AL PRB 57,1505 (1998)
!=======================================================================

!-----------------------------------------------------------------------
!     HARTREE (LDA+U - LSDA) CONTRIBUTION
!-----------------------------------------------------------------------


      DETOT  = 0.D0
      TEMPO1 = 0.D0 

      DO ISPIN=1,NSPIN 

       DO ILM1=1,MAX_LM 
       DO ILM2=1,MAX_LM 
         IF (L_OF_LMN(ILM1) /= 2) CYCLE 
         IF (L_OF_LMN(ILM2) /= 2) CYCLE 
           IM1=M_OF_LMN(ILM1)
           IM2=M_OF_LMN(ILM2)
 
          TEMPO1 = TEMPO1 - (XU1 -XJ1) * 0.5D0  & 
                          * RHO(IM1,IM2,ISPIN) * RHO(IM2,IM1,ISPIN) 

         IF (ILM1 == ILM2) THEN 
            TEMPO1 = TEMPO1 + (XU1 -XJ1) * 0.5D0 * RHO(IM1,IM1,ISPIN) 
         ENDIF 
       ENDDO !(ILM2) 
       ENDDO !(ILM1) 
      ENDDO !(ISPIN) 

  
      DETOT= TEMPO1
      IF (PRINTOCC /= 0) THEN 
        WRITE(*,*) '   '
        WRITE(*,FMT='("LDA+U DOUBLE COUNTING ")')  
        WRITE(*,FMT='("==========================")') 
        WRITE(*,FMT='("E(e-e)-E(LSDA) =",F10.5)') DETOT
      ENDIF 

      END SUBROUTINE POT_TYPE4

