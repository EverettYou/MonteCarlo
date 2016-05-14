! TO COMPILE -----------
! 1. create signature file by:
!    f2py MC.f90 -m MC -h MC.pyf --overwrite-signature
! 2. remove private variables and procedures in MC.pyf
! 3. build extension by:
!    f2py -c MC.pyf MC.f90 --opt='-O3' --build-dir build
! +-----------------------------------+
! |   Fortran core for Monte Carlo    |
! +-----------------------------------+
MODULE CORE
    USE IEEE_ARITHMETIC
    IMPLICIT NONE
    ! data pool linked to python
    ! system variables
    INTEGER :: DOF    ! onsite degrees of freedom
    INTEGER :: NSITE  ! number of sites
    INTEGER :: NA, NB, NC ! number of sublattice sites
    INTEGER :: NLST   ! size of JLST and KLST
    INTEGER, ALLOCATABLE :: CHI(:,:)  ! chi table
    INTEGER, ALLOCATABLE :: IRNG(:)   ! adjacency range list
    INTEGER, ALLOCATABLE :: JLST(:)   ! adjacent site index list
    REAL(8), ALLOCATABLE :: KLST(:)   ! energy coefficient list
    ! state variables
    REAL(8) :: ENERGY ! energy
    INTEGER, ALLOCATABLE :: CONFIG(:) ! configuration    
    INTEGER, ALLOCATABLE :: HIST(:)   ! histogram of spin (bin count)
    ! private workspace
    INTEGER :: NBLK ! number of bulk sites
    REAL(8), ALLOCATABLE :: BIAS(:,:)   ! local bias field
    REAL(8), ALLOCATABLE :: WEIGHT(:,:) ! local weights
    REAL(8), ALLOCATABLE :: RND(:)      ! random numbers
CONTAINS
    ! initialization
    SUBROUTINE INIT()
        ! allocate workspace
        ! bulk sites = A + B sites
        NBLK = NA + NB
        IF (ALLOCATED(BIAS)) THEN
            IF (ANY(SHAPE(BIAS) /= [DOF, NBLK])) THEN
                DEALLOCATE(BIAS)
                ALLOCATE(BIAS(DOF, NBLK))
            END IF
        ELSE
            ALLOCATE(BIAS(DOF, NBLK))
        END IF
        IF (ALLOCATED(WEIGHT)) THEN
            IF (ANY(SHAPE(WEIGHT) /= [DOF, NBLK])) THEN
                DEALLOCATE(WEIGHT)
                ALLOCATE(WEIGHT(DOF, NBLK))
            END IF
        ELSE
            ALLOCATE(WEIGHT(DOF, NBLK))
        END IF
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
    ! k steps of MC update (each update runs over all sites)
    ! MODE = 0: update without physical observables
    ! MODE = 1: update with physical observables
    SUBROUTINE RUN(STEPS, MODE)
        INTEGER, INTENT(IN) :: STEPS, MODE
        INTEGER :: STEP
        SELECT CASE (MODE)
        CASE(0)
            DO STEP = 1, STEPS
                CALL RANDOM_NUMBER(RND) ! prepare RND
                CALL SAMPLE0(1, NA)
                CALL SAMPLE0(NA+1, NA+NB)
            END DO
            ! set energy = NaN indicates: energy is not up to date
            ENERGY = IEEE_VALUE(1.D0, IEEE_QUIET_NAN)
            HIST(1) = -1 ! set HIST(1)=-1 indicates: HIST is not up to date
        CASE(1)
            IF (IEEE_IS_NAN(ENERGY)) CALL GET_ENERGY()
            IF (HIST(1) < 0) CALL GET_HIST()
            DO STEP = 1, STEPS
                CALL RANDOM_NUMBER(RND) ! prepare RND
                CALL SAMPLE1(1, NA)
                CALL SAMPLE1(NA+1, NA+NB)            
            END DO
        END SELECT
    END SUBROUTINE RUN
    ! block Gibbs sampling (mode 0: no update physical variables)
    SUBROUTINE SAMPLE0(LB, UB)
        INTEGER, INTENT(IN) :: LB, UB
        INTEGER :: I
        ! set BIAS, WEIGHT, and adjust RND for random choice
        CALL SET_BLOCK(LB, UB)
        ! make random choice by binary search
        FORALL (I = LB:UB)
            CONFIG(I) = CHOOSE(WEIGHT(:, I),RND(I))
        END FORALL
    END SUBROUTINE SAMPLE0
    ! block Gibbs sampling (mode 1: update physical variables)
    SUBROUTINE SAMPLE1(LB, UB)
        INTEGER, INTENT(IN) :: LB, UB
        INTEGER :: I, G0, G1
        ! set BIAS, WEIGHT, and adjust RND for random choice
        CALL SET_BLOCK(LB, UB)
        ! make random choice by binary search
        DO I = LB, UB
            G0 = CONFIG(I)
            G1 = CHOOSE(WEIGHT(:, I),RND(I))
            CONFIG(I) = G1
            ENERGY = ENERGY + BIAS(G0, I) - BIAS(G1, I)
            HIST(G0) = HIST(G0) - 1
            HIST(G1) = HIST(G1) + 1
        END DO
    END SUBROUTINE SAMPLE1
    ! calculate accumulated weights
    SUBROUTINE SET_BLOCK(LB, UB)
        INTEGER, INTENT(IN) :: LB, UB
        INTEGER :: I
        ! set BIAS, WEIGHT block, and renormalize RND
        FORALL (I = LB:UB)
            BIAS(:,I) = SET_BIAS(I) ! faster than MATMUL
            WEIGHT(:,I) = SET_WEIGHT(BIAS(:,I)) ! combined
            RND(I) = RND(I) * WEIGHT(DOF, I)            
        END FORALL
    END SUBROUTINE SET_BLOCK
    ! calculate bias field (as a vector of DOF components)
    ! given the site index I
    PURE FUNCTION SET_BIAS(I) RESULT (B)
        INTEGER, INTENT(IN) :: I
        REAL(8) :: B(DOF)
        INTEGER :: J
        B = 0.D0
        DO J = IRNG(I), IRNG(I+1)-1
            B = B + CHI(:,CONFIG(JLST(J)))*KLST(J)
        END DO
        ! this is several times faster then the MATMUL implementation
        ! because on matrix is constructed before multiplication
    END FUNCTION SET_BIAS
    ! calculate accumulated weight
    ! given the site index I
    PURE FUNCTION SET_WEIGHT(B) RESULT (W)
        REAL(8), INTENT(IN) :: B(DOF)
        REAL(8) :: W(DOF)
        INTEGER :: K
        W(1) = EXP(B(1))
        ! directly accumulate CDF, not constructing PDF first
        DO K = 2, DOF
            W(K) = W(K-1) + EXP(B(K))
        END DO
    END FUNCTION SET_WEIGHT
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
    ! calculate energy
    SUBROUTINE GET_ENERGY()
        INTEGER :: I
        ENERGY = 0.D0
        DO I = 1, NSITE
            ENERGY = ENERGY - DOT_PRODUCT(CHI(CONFIG(I), &
                     CONFIG(JLST(IRNG(I):IRNG(I+1)-1))), &
                     KLST(IRNG(I):IRNG(I+1)-1))
        END DO
    END SUBROUTINE GET_ENERGY 
    ! collect total spin (histogram)
    SUBROUTINE GET_HIST()
        INTEGER :: I, G
        HIST = 0
        DO I = 1, NBLK
            G = CONFIG(I)
            HIST(G) = HIST(G) + 1
        END DO
    END SUBROUTINE GET_HIST
    ! core dump
    SUBROUTINE DUMP()        
        PRINT *, "MC.core dump to MC.core.dat"
        OPEN (UNIT=99, FILE="MC.core.dat", STATUS="REPLACE", ACCESS="STREAM")
        WRITE(99) DOF, NSITE, NA, NB, NC, NLST, NBLK
        WRITE(99) CHI, IRNG, JLST, KLST
        WRITE(99) ENERGY, CONFIG, HIST
        WRITE(99) BIAS, WEIGHT, RND
        CLOSE(99)
    END SUBROUTINE DUMP
    ! core load
    SUBROUTINE LOAD() 
        PRINT *, "MC.core load from MC.core.dat"
        OPEN (UNIT=99, FILE="MC.core.dat", STATUS="UNKNOWN", ACCESS="STREAM")
        READ(99) DOF, NSITE, NA, NB, NC, NLST, NBLK
        IF (ALLOCATED(CHI)) DEALLOCATE(CHI)
        IF (ALLOCATED(IRNG)) DEALLOCATE(IRNG)
        IF (ALLOCATED(JLST)) DEALLOCATE(JLST)
        IF (ALLOCATED(KLST)) DEALLOCATE(KLST)
        IF (ALLOCATED(CONFIG)) DEALLOCATE(CONFIG)
        IF (ALLOCATED(HIST)) DEALLOCATE(HIST)
        IF (ALLOCATED(BIAS)) DEALLOCATE(BIAS)
        IF (ALLOCATED(WEIGHT)) DEALLOCATE(WEIGHT)
        IF (ALLOCATED(RND)) DEALLOCATE(RND)
        ALLOCATE(CHI(DOF,DOF), IRNG(NSITE+1), JLST(NLST), KLST(NLST), &
                 CONFIG(NSITE), HIST(DOF), &
                 BIAS(DOF,NBLK), WEIGHT(DOF,NBLK), RND(NBLK))
        READ(99) CHI, IRNG, JLST, KLST
        READ(99) ENERGY, CONFIG, HIST
        READ(99) BIAS, WEIGHT, RND
        CLOSE(99)
    END SUBROUTINE LOAD
END MODULE CORE
! Measurement module
MODULE PHYSICS
    USE CORE
    IMPLICIT NONE
    ! data variables
    INTEGER :: NSPIN ! number of spins to monitor
    INTEGER, ALLOCATABLE :: MONITOR(:) ! sites for monitoring
    REAL(8) :: ENERGY1, ENERGY2
    REAL(8), ALLOCATABLE :: MAGNET1(:), MAGNET2(:,:) ! (DOF, DOF)
    REAL(8), ALLOCATABLE :: SPINS(:,:) ! (DOF, NSPIN)
CONTAINS
    ! launch measurement environment 
    SUBROUTINE LAUNCH()
        ! magnetization 1st moment (vector)
        IF (ALLOCATED(MAGNET1)) THEN
            IF (ANY(SHAPE(MAGNET1) /= [DOF])) THEN
                DEALLOCATE(MAGNET1)
                ALLOCATE(MAGNET1(DOF))
            END IF
        ELSE
            ALLOCATE(MAGNET1(DOF))
        END IF
        ! magnetization 2nd moment (matrix)
        IF (ALLOCATED(MAGNET2)) THEN
            IF (ANY(SHAPE(MAGNET2) /= [DOF, DOF])) THEN
                DEALLOCATE(MAGNET2)
                ALLOCATE(MAGNET2(DOF, DOF))
            END IF
        ELSE
            ALLOCATE(MAGNET2(DOF, DOF))
        END IF
        ! spin 1st moment (array of vectors)
        ! all higher order moments are the same as 1st moment
        IF (ALLOCATED(SPINS)) THEN
            IF (ANY(SHAPE(SPINS) /= [DOF, NSPIN])) THEN
                DEALLOCATE(SPINS)
                ALLOCATE(SPINS(DOF, NSPIN))
            END IF
        ELSE
            ALLOCATE(SPINS(DOF, NSPIN))
        END IF
        ! get energy and hist ready
        IF (IEEE_IS_NAN(ENERGY)) CALL GET_ENERGY()
        IF (HIST(1) < 0) CALL GET_HIST()
    END SUBROUTINE LAUNCH
    ! measure energy and spin data over steps
    SUBROUTINE MEASURE(STEPS)
        INTEGER, INTENT(IN) :: STEPS
        INTEGER :: STEP, I, J
        REAL(8) :: MAGNET(DOF)
        
        CALL LAUNCH() ! launch measurement environment
        ! clear data pool
        ENERGY1 = 0.D0
        ENERGY2 = 0.D0
        MAGNET1 = 0.D0
        MAGNET2 = 0.D0
        SPINS = 0.D0
        ! run MC with updates, and collect data
        DO STEP = 1, STEPS
            CALL RANDOM_NUMBER(RND) ! prepare RND
            CALL SAMPLE1(1, NA)
            CALL SAMPLE1(NA+1, NA+NB)
            ! energy measurement
            ENERGY1 = ENERGY1 + ENERGY
            ENERGY2 = ENERGY2 + ENERGY**2
            ! magnetization measurement
            FORALL (I = 1:DOF)
                MAGNET(I) = REAL(HIST(I),8)/NBLK
            END FORALL
            MAGNET1 = MAGNET1 + MAGNET
            FORALL (I = 1:DOF, J = 1:DOF)
                MAGNET2(I, J) = MAGNET2(I, J) + MAGNET(I)*MAGNET(J)
            END FORALL
            ! spin measurement
            FORALL (I = 1:NSPIN)
                SPINS(CONFIG(MONITOR(I)), I) = SPINS(CONFIG(MONITOR(I)), I) + 1.D0
            END FORALL
        END DO
        ! normalized by STEPS
        ENERGY1 = ENERGY1/STEPS
        ENERGY2 = ENERGY2/STEPS
        MAGNET1 = MAGNET1/STEPS
        MAGNET2 = MAGNET2/STEPS
        SPINS = SPINS/STEPS
    END SUBROUTINE MEASURE
END MODULE PHYSICS

! ========== for debug use ==========
! SUBROUTINE TEST_CHOOSE()
!     USE CORE
!     INTEGER :: L
!     DOF = 6
!     L = CHOOSE([2._8,3._8,7._8,15._8,20._8,23._8],22.99_8)
!     PRINT *, L
! END SUBROUTINE TEST_CHOOSE
! SUBROUTINE TEST_RUN()
!     USE CORE
!     CALL RUN(1,1)
!     PRINT *, ENERGY
! END SUBROUTINE TEST_RUN
! SUBROUTINE TEST_MEASURE()
!     USE PHYSICS
!     CALL INIT_RAND_SEED()
!     CALL MEASURE(10)
!     PRINT *, NSPIN
!     PRINT *, ENERGY1, ENERGY2
!     PRINT *, MAGNET1
!     PRINT *, CONFIG
! END SUBROUTINE TEST_MEASURE
! 
! PROGRAM MAIN
!     USE CORE
!     
!     CALL LOAD()
! !     CALL TEST_RUN()
! !     CALL TEST_MEASURE()
! !     CALL TEST_CHOOSE()
! END PROGRAM MAIN