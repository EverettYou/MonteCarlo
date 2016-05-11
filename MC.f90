! TO COMPILE -----------
! 1. create signature file by:
!    f2py MC.f90 -m MC -h MC.pyf
! 2. build extension by:
!    f2py -c MC.pyf MC.f90
! +-----------------------------------+
! |   Fortran core for Monte Carlo    |
! +-----------------------------------+
MODULE CORE
    IMPLICIT NONE
    ! data pool linked to python
    INTEGER :: DOF ! onsite degrees of freedom
    INTEGER :: NSITE ! number of sites
    INTEGER :: NA, NB, NC ! number of sublattice sites
    INTEGER, ALLOCATABLE :: CONFIG(:) ! configuration
    INTEGER, ALLOCATABLE :: CHI(:,:)  ! chi table
    INTEGER, ALLOCATABLE :: IRNG(:)   ! adjacency range list
    INTEGER, ALLOCATABLE :: JLST(:)   ! adjacent site index list
    REAL(8), ALLOCATABLE :: KLST(:)   ! energy coefficient list
    ! private workspace
    INTEGER :: NBLK ! number of bulk sites
    REAL(8), ALLOCATABLE :: BIAS(:,:)   ! local bias field
    REAL(8), ALLOCATABLE :: WEIGHT(:,:) ! local weights
    REAL(8), ALLOCATABLE :: RND(:)      ! random numbers
CONTAINS
    ! initialization
    SUBROUTINE INIT()
        ! allocate workspace
        ! BIAS and WEIGHT are of the shape [DOF, NSITE]
        IF (ALLOCATED(BIAS)) THEN
            IF (ANY(SHAPE(BIAS) /= [DOF, NSITE])) THEN
                DEALLOCATE(BIAS)
                ALLOCATE(BIAS(DOF,NSITE))
            END IF
        ELSE
            ALLOCATE(BIAS(DOF,NSITE))
        END IF
        IF (ALLOCATED(WEIGHT)) THEN
            IF (ANY(SHAPE(WEIGHT) /= [DOF, NSITE])) THEN
                DEALLOCATE(WEIGHT)
                ALLOCATE(WEIGHT(DOF,NSITE))
            END IF
        ELSE
            ALLOCATE(WEIGHT(DOF,NSITE))
        END IF
        ! bulk sites = A + B sites
        NBLK = NA + NB
        ! RND will only use bulk sites
        IF (ALLOCATED(RND)) THEN
            IF (ANY(SHAPE(RND) /= [NBLK])) THEN
                DEALLOCATE(RND)
                ALLOCATE(RND(NBLK))
            END IF
        ELSE
            ALLOCATE(RND(NBLK))
        END IF
        ! random seed
        CALL INIT_RAND_SEED()
    END SUBROUTINE INIT
    ! random seed
    SUBROUTINE INIT_RAND_SEED()
        INTEGER :: I, N, CLOCK
        INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
        CALL RANDOM_SEED(SIZE = N)
        ALLOCATE(SEED(N))
        CALL SYSTEM_CLOCK(COUNT = CLOCK)
        SEED = CLOCK + 37*[(I-1,I=1,N)]
        CALL RANDOM_SEED(PUT = SEED)
        DEALLOCATE(SEED)
    END SUBROUTINE INIT_RAND_SEED
    ! one Monte Carlo step goes over odd and even sites
    SUBROUTINE MCSTEP()
        ! prepare RND
        CALL RANDOM_NUMBER(RND)
        CALL SAMPLE(1, NA)
        CALL SAMPLE(NA+1, NA+NB)
    END SUBROUTINE MCSTEP
    ! block Gibbs sampling
    SUBROUTINE SAMPLE(LB, UB)
        INTEGER, INTENT(IN) :: LB, UB
        INTEGER :: I
        ! set BIAS, WEIGHT, and adjust RND for random choice
        CALL SET_BLOCK(LB, UB)
        ! make random choice by binary search
        FORALL (I = LB:UB)
            CONFIG(I) = CHOOSE(WEIGHT(:, I),RND(I))
        END FORALL
    END SUBROUTINE SAMPLE
    ! calculate accumulated weights
    SUBROUTINE SET_BLOCK(LB, UB)
        INTEGER, INTENT(IN) :: LB, UB
        INTEGER :: I
        ! calculate BIAS and WEIGHT
        FORALL (I = LB:UB)
            BIAS(:,I) = MATMUL(CHI(:,CONFIG(JLST(IRNG(I):IRNG(I+1)-1))),&
                        KLST(IRNG(I):IRNG(I+1)-1))
            WEIGHT(:,I) = EXP(BIAS(:,I))
        END FORALL
        ! in-place accumulate WEIGHT to CDF (unnormalized)
        DO I = 2, DOF
            WEIGHT(I, LB:UB) = WEIGHT(I, LB:UB) + WEIGHT(I-1, LB:UB)
        END DO
        ! multiply RND by the total weight
        RND(LB:UB) = RND(LB:UB) * WEIGHT(DOF, LB:UB)
    END SUBROUTINE SET_BLOCK
    ! random choice
    ! given CDF array W, and random number X
    ! return a choice in [1:DOF]
    PURE FUNCTION CHOOSE(W, X) RESULT (R)
        REAL(8), INTENT(IN) :: W(DOF) ! weight array
        REAL(8), INTENT(IN) :: X ! a random number
        INTEGER :: L, R, M
        ! binary search
        L = 0
        R = DOF
        ! if R > L+1 do the search
        DO WHILE (R - L > 1)
            M = (L + R)/2 ! set mid point
            IF (W(M) >= X) THEN
                R = M
            ELSE
                L = M
            END IF
        END DO
        ! now X is determined in the bin (L,R=L+1)
        ! L will be the result of the choice
    END FUNCTION CHOOSE
    ! core dump
    SUBROUTINE DUMP()
        INTEGER :: NLST
        NLST = SIZE(JLST)
        
        PRINT *, "MC.core dump to MC.core.dat"
        OPEN (UNIT=99, FILE="MC.core.dat", STATUS="REPLACE", ACCESS="STREAM")
        WRITE(99) DOF, NSITE, NA, NB, NC, NLST
        WRITE(99) CONFIG, CHI, IRNG, JLST, KLST
        WRITE(99) NBLK
        WRITE(99) BIAS, WEIGHT, RND
        CLOSE(99)
    END SUBROUTINE DUMP
    ! core load
    SUBROUTINE LOAD()
        INTEGER :: NLST
        
        PRINT *, "MC.core load from MC.core.dat"
        OPEN (UNIT=99, FILE="MC.core.dat", STATUS="UNKNOWN", ACCESS="STREAM")
        READ(99) DOF, NSITE, NA, NB, NC, NLST
        IF (ALLOCATED(CONFIG)) DEALLOCATE(CONFIG)
        IF (ALLOCATED(CHI)) DEALLOCATE(CHI)
        IF (ALLOCATED(IRNG)) DEALLOCATE(IRNG)
        IF (ALLOCATED(JLST)) DEALLOCATE(JLST)
        IF (ALLOCATED(KLST)) DEALLOCATE(KLST)
        ALLOCATE(CONFIG(NSITE), CHI(DOF,DOF), IRNG(NSITE+1), JLST(NLST), KLST(NLST))
        READ(99) CONFIG, CHI, IRNG, JLST, KLST
        READ(99) NBLK
        IF (ALLOCATED(BIAS)) DEALLOCATE(BIAS)
        IF (ALLOCATED(WEIGHT)) DEALLOCATE(WEIGHT)
        IF (ALLOCATED(RND)) DEALLOCATE(RND)
        ALLOCATE(BIAS(DOF,NSITE), WEIGHT(DOF,NSITE), RND(NBLK))
        READ(99) BIAS, WEIGHT, RND
        CLOSE(99)
    END SUBROUTINE LOAD

END MODULE CORE

! ========== for debug use ==========
! SUBROUTINE TEST_CHOOSE()
!     USE CORE
!     INTEGER :: L
!     DOF = 6
!     L = CHOOSE([2._8,3._8,7._8,15._8,20._8,23._8],22.99_8)
!     PRINT *, L
! END SUBROUTINE TEST_CHOOSE
! SUBROUTINE TEST_MCSTEP()
!     USE CORE
!     CALL MCSTEP()
!     PRINT *, CONFIG
! END SUBROUTINE TEST_MCSTEP
! 
! PROGRAM MAIN
!     USE CORE
!     
!     CALL LOAD()
!     CALL TEST_MCSTEP()
! !     CALL TEST_CHOOSE()
! END PROGRAM MAIN