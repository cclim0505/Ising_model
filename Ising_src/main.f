        PROGRAM MAIN
        USE lattice_config
        USE simulation
        USE average
        IMPLICIT NONE
        INTEGER :: iter,jter
        CALL read_parameters   ! read number of loops & temperature
        CALL read_lattice      ! read lattice_size
        CALL continue_previous ! continue from previous state if old files exist
        CALL initiate_lattice
!========================================
! Main simulation loop
!========================================
        DO iter =1,loop
          CALL flip_lattice
          CALL calc_total_energy
          CALL monte_carlo               ! accept or reject 
        END DO
        CALL stop_simulation
!========================================
! Calculate average values
!========================================
        CALL read_in_val
        CALL calc_averages
        END PROGRAM MAIN

