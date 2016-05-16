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
    REAL,    ALLOCATABLE :: KLST(:)   ! action coefficient list
    ! state variables
    REAL :: ACTION ! action
    REAL :: BETA   ! inverse temperature
    INTEGER, ALLOCATABLE :: CONFIG(:) ! configuration    
    INTEGER, ALLOCATABLE :: HIST(:)   ! histogram of spin (bin count)
    ! private workspace
    INTEGER :: NBLK ! number of bulk sites
    REAL, ALLOCATABLE :: BIAS(:,:)   ! local bias field
    REAL, ALLOCATABLE :: WEIGHT(:,:) ! local weights
    REAL, ALLOCATABLE :: RND(:)      ! random numbers
    REAL, ALLOCATABLE :: KLST1(:)    ! an initial copy KLST for tempering
CONTAINS
    ! ------ Initialization ------
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
        ! make an initial copy of KLST
        ! it is better to update KLST from its initial copy
        ! to avoid error accumulation after several updates
        IF (NLST > 0) THEN
            KLST1 = KLST
        END IF
    END SUBROUTINE INIT
    ! automatic random seed
    SUBROUTINE INIT_RAND_SEED()
        INTEGER :: I, N, CLOCK
        INTEGER, DIMENSION(:), ALLOCATABLE :: SD
        CALL RANDOM_SEED(SIZE = N)
        ALLOCATE(SD(N))
        CALL SYSTEM_CLOCK(COUNT = CLOCK)
        SD = CLOCK + 37*[(I-1,I=1,N)]
        CALL RANDOM_SEED(PUT = SD)
        DEALLOCATE(SD)
    END SUBROUTINE INIT_RAND_SEED
    ! manually set random seed
    SUBROUTINE SEED(S)
        INTEGER, INTENT(IN) :: S
        INTEGER :: I, N
        INTEGER, DIMENSION(:), ALLOCATABLE :: SD
        CALL RANDOM_SEED(SIZE = N)
        ALLOCATE(SD(N))
        SD = S + 37*[(I-1,I=1,N)]
        CALL RANDOM_SEED(PUT = SD)
        DEALLOCATE(SD)
    END SUBROUTINE SEED
    ! ------ Sampling ------
    ! k steps of MC update (each update runs over all sites)
    SUBROUTINE RUN(STEPS)
        INTEGER, INTENT(IN) :: STEPS
        INTEGER :: STEP
        ! run a simple sampling without keeping track of states
        DO STEP = 1, STEPS
            CALL RANDOM_NUMBER(RND) ! prepare RND
            CALL SAMPLE0(1, NA)
            CALL SAMPLE0(NA+1, NA+NB)
        END DO
        ! recalculate the action and hist after running
        CALL GET_ACTION()
        CALL GET_HIST()
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
            ACTION = ACTION + BIAS(G0, I) - BIAS(G1, I)
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
        REAL :: B(DOF)
        INTEGER :: J
        B = 0.
        DO J = IRNG(I), IRNG(I+1)-1
            B = B + CHI(:,CONFIG(JLST(J)))*KLST(J)
        END DO
        ! this is several times faster then the MATMUL implementation
        ! because on matrix is constructed before multiplication
    END FUNCTION SET_BIAS
    ! calculate accumulated weight
    ! given the site index I
    PURE FUNCTION SET_WEIGHT(B) RESULT (W)
        REAL, INTENT(IN) :: B(DOF)
        REAL :: W(DOF)
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
        REAL, INTENT(IN) :: W(DOF) ! weight array
        REAL, INTENT(IN) :: X ! a random number
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
    ! ------ State Modification ------
    ! set beta and update action coefficients
    SUBROUTINE SET_BETA(NEW_BETA)
        REAL, INTENT(IN) :: NEW_BETA 
        BETA = NEW_BETA ! set new beta
        ! update action coefficients
        IF (NLST > 0) THEN
            KLST = KLST1 * BETA ! from the initial copy
        END IF
        ! action is also changed by resetting temperature
        CALL GET_ACTION() ! recalculate action
    END SUBROUTINE SET_BETA
    ! calculate action
    SUBROUTINE GET_ACTION()
        INTEGER :: I
        ACTION = 0.
        DO I = 1, NSITE ! for all sites
            ACTION = ACTION - DOT_PRODUCT(CHI(CONFIG(I), &
                     CONFIG(JLST(IRNG(I):IRNG(I+1)-1))), &
                     KLST(IRNG(I):IRNG(I+1)-1))
        END DO
        ! each bond is doubled counted, so 1/2
        ACTION = ACTION/2
    END SUBROUTINE GET_ACTION 
    ! collect total spin (histogram)
    SUBROUTINE GET_HIST()
        INTEGER :: I, G
        HIST = 0
        DO I = 1, NBLK ! for bulk sites only
            G = CONFIG(I)
            HIST(G) = HIST(G) + 1
        END DO
    END SUBROUTINE GET_HIST
    ! ------ I/O System ------
    ! core dump
    SUBROUTINE DUMP()        
        PRINT *, "MC.core dump to MC.core.dat"
        OPEN (UNIT=99, FILE="MC.core.dat", STATUS="REPLACE", ACCESS="STREAM")
        WRITE(99) DOF, NSITE, NA, NB, NC, NLST, NBLK
        WRITE(99) CHI, IRNG, JLST, KLST
        WRITE(99) ACTION, BETA, CONFIG, HIST
        WRITE(99) BIAS, WEIGHT, RND, KLST1
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
        IF (ALLOCATED(KLST1)) DEALLOCATE(KLST1)
        ALLOCATE(CHI(DOF,DOF), IRNG(NSITE+1), JLST(NLST), KLST(NLST), &
                 CONFIG(NSITE), HIST(DOF), &
                 BIAS(DOF,NBLK), WEIGHT(DOF,NBLK), RND(NBLK), &
                 KLST1(NLST))
        READ(99) CHI, IRNG, JLST, KLST
        READ(99) ACTION, BETA, CONFIG, HIST
        READ(99) BIAS, WEIGHT, RND, KLST1
        CLOSE(99)
    END SUBROUTINE LOAD
END MODULE CORE
! Measurement module
MODULE DATA
    USE CORE
    IMPLICIT NONE
    ! data variables
    INTEGER :: NSPIN ! number of spins to monitor
    INTEGER, ALLOCATABLE :: MONITOR(:) ! sites for monitoring
    REAL :: ENERGY1, ENERGY2
    REAL, ALLOCATABLE :: MAGNET1(:), MAGNET2(:,:) ! (DOF, DOF)
    REAL, ALLOCATABLE :: SPINS(:,:) ! (DOF, NSPIN)
    ! time series
    INTEGER :: NTS ! length of the time series
    REAL, ALLOCATABLE :: ETS(:)    ! energy density series
    REAL, ALLOCATABLE :: MTS(:,:) ! magnetization density series
CONTAINS
    ! measure energy and spin data over steps
    SUBROUTINE MEASURE(STEPS)
        INTEGER, INTENT(IN) :: STEPS
        INTEGER :: STEP, I, J
        REAL :: ENERGY, MAGNET(DOF)
        
        ! launch measurement environment
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
        CALL GET_ACTION()
        CALL GET_HIST()
        ! clear data pool
        ENERGY1 = 0.
        ENERGY2 = 0.
        MAGNET1 = 0.
        MAGNET2 = 0.
        SPINS = 0.
        ! run MC with updates, and collect data
        DO STEP = 1, STEPS
            CALL RANDOM_NUMBER(RND) ! prepare RND
            CALL SAMPLE1(1, NA)
            CALL SAMPLE1(NA+1, NA+NB)
            ! energy measurement
            ENERGY = ACTION/BETA/NSITE
            ENERGY1 = ENERGY1 + ENERGY
            ENERGY2 = ENERGY2 + ENERGY**2
            ! magnetization measurement
            MAGNET = REAL(HIST)/NBLK
            MAGNET1 = MAGNET1 + MAGNET
            FORALL (I = 1:DOF, J = 1:DOF)
                MAGNET2(I, J) = MAGNET2(I, J) + MAGNET(I)*MAGNET(J)
            END FORALL
            ! spin measurement
            FORALL (I = 1:NSPIN)
                SPINS(CONFIG(MONITOR(I)), I) = SPINS(CONFIG(MONITOR(I)), I) + 1.
            END FORALL
        END DO
        ! normalized by STEPS
        ENERGY1 = ENERGY1/STEPS
        ENERGY2 = ENERGY2/STEPS
        MAGNET1 = MAGNET1/STEPS
        MAGNET2 = MAGNET2/STEPS
        SPINS = SPINS/STEPS
    END SUBROUTINE MEASURE
    ! collect time series
    SUBROUTINE COLLECT(STEPS)
        INTEGER, INTENT(IN) :: STEPS
        INTEGER :: T
        NTS = STEPS
        ! reallocation for ETS, MTS
        IF (ALLOCATED(ETS).OR.ALLOCATED(MTS)) THEN
            IF (ANY(SHAPE(ETS) /= [NTS]) .OR. &
                ANY(SHAPE(MTS) /= [DOF,NTS])) THEN
                IF (ALLOCATED(ETS)) DEALLOCATE(ETS)
                IF (ALLOCATED(MTS)) DEALLOCATE(MTS)
                ALLOCATE(ETS(NTS), MTS(DOF, NTS))
            END IF
        ELSE
            ALLOCATE(ETS(NTS), MTS(DOF, NTS))
        END IF
        ! get energy and hist ready
        CALL GET_ACTION()
        CALL GET_HIST()
        ! start collect time series
        DO T = 1, NTS
            CALL RANDOM_NUMBER(RND) ! prepare RND
            CALL SAMPLE1(1, NA)
            CALL SAMPLE1(NA+1, NA+NB)
            ETS(T) = ACTION/BETA/NSITE
            MTS(:,T) = REAL(HIST)/NBLK
        END DO
    END SUBROUTINE COLLECT
END MODULE DATA

! ========== for debug use ==========
! SUBROUTINE TEST_CHOOSE()
!     USE CORE
!     INTEGER :: L
!     DOF = 6
!     L = CHOOSE([2.,3.,7.,15.,20.,23.],22.99)
!     PRINT *, L
! END SUBROUTINE TEST_CHOOSE
! SUBROUTINE TEST_RUN()
!     USE CORE
!     !CALL RUN(1,1)
!     PRINT *, CONFIG
!     CALL GET_ACTION()
!     PRINT *, ACTION
! END SUBROUTINE TEST_RUN
! SUBROUTINE TEST_MEASURE()
!     USE DATA
!     CALL INIT_RAND_SEED()
!     CALL MEASURE(100)
!     PRINT *, NSPIN
!     PRINT *, ENERGY1, ENERGY2
!     PRINT *, MAGNET1
!     PRINT *, CONFIG
! END SUBROUTINE TEST_MEASURE
! SUBROUTINE TEST_SETBETA()
!     USE CORE
!     PRINT *, BETA, ACTION
!     PRINT *, KLST(1:3),'...'
!     PRINT *, KLST1(1:3),'...'
!     CALL SET_BETA(1.)
!     PRINT *, BETA, ACTION
!     PRINT *, KLST(1:3),'...'
!     PRINT *, KLST1(1:3),'...'
! END SUBROUTINE TEST_SETBETA
! 
! PROGRAM MAIN
!     USE CORE
!     
!     CALL LOAD()
!     CALL TEST_SETBETA()
! !     CALL TEST_RUN()
! !     CALL TEST_MEASURE()
! !     CALL TEST_CHOOSE()
! END PROGRAM MAIN