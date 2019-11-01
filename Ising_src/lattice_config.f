        MODULE lattice_config
        USE simulation
!============================================================================
!       Variable dictionary
!============================================================================
!       SIMULATION PARAMETERS
        REAL,PARAMETER          :: mag_h = 0.01     !magnetic field h
        REAL,PARAMETER          :: kb    = 1.0      !Boltzmann constant
!============================================================================
!       LATTICE DIMENSIONS
        INTEGER                 :: lattice_size 
        INTEGER,ALLOCATABLE     :: lattice (:,:)    !2D Ising model lattice
!============================================================================
!       ENERGIES & MAGNETIZATION
        REAL                    :: previous_energy = 0.0
        REAL                    :: energy_diff
        REAL                    :: total_energy
        REAL                    :: pot_energy
        REAL                    :: int_energy
        REAL                    :: magnetization
        INTEGER                 :: total_spin
!============================================================================
!       MONTE CARLO & PROBABILITIES
        REAL                    :: prob_ene  !prob as a function of energy
        REAL                    :: prob_diff !prob as a function of energy difference
        INTEGER                 :: row, col           ! row and column indices
        LOGICAL                 :: isaccept = .false. ! true=accept
!============================================================================
!       CONTINUE FROM PREVIOUS SIMULATION WHEN NECESSARY
        INTEGER                 :: state_count = 1      !number of configuration
        LOGICAL                 :: exist_count =.false.
        LOGICAL                 :: exist_state =.false.
        LOGICAL                 :: iscontinue = .false.
        CONTAINS
!============================================================================

        SUBROUTINE read_lattice
! reads lattice size and allocates lattice arrays
        OPEN(40,file='lattice_size.dat',status='old')
        READ(40,*) lattice_size
        CLOSE(40)
        ALLOCATE(lattice(lattice_size,lattice_size))
        END SUBROUTINE read_lattice

        SUBROUTINE continue_previous
! allows to continue simulation from previous saved state if old files
! exist
        IMPLICIT NONE
        INTEGER         :: iter,jter
        INTEGER         :: dummy
        INQUIRE(file="final_state_count.dat", exist=exist_count)
        INQUIRE(file="final_state.dat", exist=exist_state)
        IF(exist_count .AND. exist_state) THEN
          iscontinue = .true.
          OPEN(90,file='final_state_count.dat',status='old')
          READ(90,*) state_count
          CLOSE(90)
          OPEN(91,file='final_state.dat',status='old')
          DO iter=1,lattice_size
            DO jter=1,lattice_size
               READ(91,*) dummy, dummy, lattice(iter,jter) 
            END DO
          END DO
          CLOSE(91)
        END IF
        END SUBROUTINE continue_previous

        SUBROUTINE initiate_lattice
! create an NxN 2D lattice and assign spins to each lattice point randomly
        IMPLICIT NONE
        INTEGER :: iter,jter
        REAL    :: rnum
        INTEGER :: spin                 ! spin value, +1 or -1
        IF ( .not.iscontinue ) THEN
          DO iter=1,lattice_size
            DO jter =1,lattice_size
               CALL RANDOM_NUMBER(rnum)
               IF (rnum > 0.5) THEN
                 spin = 1
               ELSE
                 spin = -1
               END IF
               lattice(iter,jter) = spin 
            END DO
          END DO
        END IF
        CALL calc_total_energy
        previous_energy = total_energy
        CALL calc_magnet
        CALL calc_prob(total_energy,prob_ene)
        OPEN(20,file='results.dat',access='append')
        IF (.not.iscontinue)  CALL save_results
        END SUBROUTINE initiate_lattice

        SUBROUTINE stop_simulation
        IMPLICIT NONE
        CLOSE(20)
        CALL save_final_state
        END SUBROUTINE stop_simulation

        SUBROUTINE save_final_state
! save final state or configuration at the end of simulation
        IMPLICIT NONE
        INTEGER         :: iter, jter
        OPEN(25,file='final_state.dat')
          DO iter = 1, lattice_size
            DO jter = 1, lattice_size
              WRITE(25,*) iter, jter, lattice(iter,jter)
            END DO
          END DO
        CLOSE(25)

        OPEN(26,file='final_state_count.dat')
        WRITE(26,*) state_count 
        CLOSE(26)
        END SUBROUTINE save_final_state

        SUBROUTINE pick_random_lattice(row,col)
! randomly pick a lattice point to be flipped
! called by SUBROUTINE flip_lattice
        IMPLICIT NONE
        INTEGER,INTENT(OUT) :: row, col    ! row and column indices
        REAL                :: rnum,cnum   ! row and column number generated randomly
        CALL RANDOM_NUMBER(rnum)
        CALL RANDOM_NUMBER(cnum)
        rnum=rnum * lattice_size 
        row=NINT(rnum)
        cnum=cnum * lattice_size 
        col=NINT(cnum)
! apply periodic boundary conditions
        IF (row == 0) THEN
          row = lattice_size
        ELSE IF (row > lattice_size) THEN
          row = 1
        END IF
        IF (col == 0) THEN
          col = lattice_size
        ELSE IF (col > lattice_size) THEN
          col = 1
        END IF
!        WRITE(*,*) "row&col are:", row,col
        END SUBROUTINE pick_random_lattice


        SUBROUTINE flip_lattice      
! action to flip lattice spin
        IMPLICIT NONE
        CALL pick_random_lattice(row,col)
        lattice(row,col) = -lattice(row,col)
        END SUBROUTINE flip_lattice

        SUBROUTINE monte_carlo
! monte carlo method, deciding when to accept or reject new state
        IMPLICIT NONE
        REAL    :: rnum
        isaccept=.false.
        IF (energy_diff < 0.0) THEN
          ! accept
          isaccept=.true.
        ELSE
          energy_diff = total_energy - previous_energy
          CALL calc_prob(energy_diff,prob_diff)
          CALL RANDOM_NUMBER(rnum)
          IF ( prob_diff > rnum ) THEN
          ! accept
          isaccept=.true.
          ELSE
          ! reject, flip back
          lattice(row,col) = -lattice(row,col)
          END IF
        END IF
        
        IF (isaccept) THEN 
          state_count = state_count + 1
!          WRITE(*,*) 'accept'
!          WRITE(*,*) 'state_count', state_count
          previous_energy = total_energy
          CALL calc_magnet
          CALL calc_prob(total_energy,prob_ene)
          CALL save_results
!          WRITE(*,*) 'total_energy', total_energy
        ELSE
!          WRITE(*,*) 'reject'
        END IF
        END SUBROUTINE monte_carlo

        SUBROUTINE calc_prob(energy,prob)         
! calculate probability as a function of difference in energy or energy
! of a state
        IMPLICIT NONE
        REAL,INTENT(IN)          :: energy
        REAL,INTENT(OUT)         :: prob
        prob = - energy /(kb*tempera)
        prob = EXP(prob) 
!        WRITE(*,*) "Probability", prob
        END SUBROUTINE calc_prob

        SUBROUTINE calc_total_spin
! calculate total spin value, considering all lattice points
! called by SUBROUTINE calc_magnet & calc_pot_energy
        IMPLICIT NONE
        INTEGER         :: iter,jter
        total_spin = 0
        DO iter=1,lattice_size
          DO jter=1,lattice_size
             total_spin = total_spin + lattice(iter,jter)
          END DO
        END DO
        END SUBROUTINE calc_total_spin

        SUBROUTINE calc_magnet  
! calculate magnetization
        IMPLICIT NONE
        INTEGER         :: iter,jter
!        CALL calc_total_spin
        magnetization = REAL(total_spin) / REAL(lattice_size**2)
!        WRITE(*,*) "Magnetization",magnetization
        END SUBROUTINE calc_magnet



        SUBROUTINE save_results  
! save physical quantities and probability
        IMPLICIT NONE
        WRITE(20,*) state_count,total_energy,magnetization
        !WRITE(21,*) state_count,energy_diff
        !WRITE(22,*) state_count,int_energy,pot_energy
        END SUBROUTINE save_results

        SUBROUTINE calc_pot_energy(energy)
! calculate potential energy of the contributing spins 
! called by SUBROUTINE calc_total_energy
        IMPLICIT NONE
        REAL,INTENT(OUT)  :: energy
        INTEGER           :: iter,jter
        CALL calc_total_spin
        energy = -mag_h * total_spin        ! negative sign introduced
!        WRITE(*,*) 'Energy is:', energy
        END SUBROUTINE calc_pot_energy


        SUBROUTINE calc_int_energy(energy)  
! calculate the interaction energy with neighbours in a 2D square
! lattice
! called by SUBROUTINE calc_total_energy
        IMPLICIT NONE
        REAL,INTENT(OUT)       :: energy
        INTEGER                :: iter,jter
        INTEGER                :: north,south,east,west ! neighbour index
        REAL,PARAMETER         :: eps = 1.0   ! epsilon
        REAL                   :: neigh       ! neighbouring energy
        energy = 0.0
        DO iter=1,lattice_size
          DO jter=1,lattice_size
             neigh = 0.0
             north = iter-1
             south = iter+1
             east  = jter+1      
             west  = jter-1      
             IF (north < 1)            north = lattice_size
             IF (south > lattice_size) south = 1
             IF (east > lattice_size)  east = 1
             IF (west < 1)             west = lattice_size
             neigh = neigh + lattice(iter,jter)*lattice(north,jter)*eps
             neigh = neigh + lattice(iter,jter)*lattice(south,jter)*eps
             neigh = neigh + lattice(iter,jter)*lattice(iter,east)*eps
             neigh = neigh + lattice(iter,jter)*lattice(iter,west)*eps
             energy = energy + neigh
          END DO
        END DO
        energy = - energy / 2.0    ! negative sign introduced & double counting fixed
!        WRITE(*,*) "Interaction energy", energy
        END SUBROUTINE calc_int_energy


        SUBROUTINE calc_total_energy
! calculate total energy, which is
! sum of potential energy & interaction energy
        IMPLICIT NONE
        CALL calc_pot_energy(pot_energy)
        CALL calc_int_energy(int_energy)
        total_energy = pot_energy + int_energy
        energy_diff = total_energy - previous_energy
!        WRITE(*,*) "Potential energy", pot_energy
!        WRITE(*,*) "Interaction energy", int_energy
!        WRITE(*,*) "Total energy", total_energy
!        WRITE(*,*) "Energy diff", energy_diff

        END SUBROUTINE calc_total_energy

        END MODULE lattice_config
