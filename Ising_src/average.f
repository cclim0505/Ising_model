        MODULE average
        USE simulation
        USE lattice_config 
        REAL,DIMENSION(:),ALLOCATABLE   :: ene           !energy array
        REAL,DIMENSION(:),ALLOCATABLE   :: ene_sq        !energy squared array
        REAL,DIMENSION(:),ALLOCATABLE   :: mag           !magnetization array
        REAL                            :: ave_ene       !average energy
        REAL                            :: ave_ene_sq    !average energy squared
        REAL                            :: ave_mag       !average magnetization
        REAL                            :: heat_capacity ! heat capacity
        CONTAINS

        SUBROUTINE read_in_val
        IMPLICIT NONE
        INTEGER         :: iter
        INTEGER         :: dummy
        ALLOCATE(ene(state_count))
        ALLOCATE(ene_sq(state_count))
        ALLOCATE(mag(state_count))
        OPEN(60,file='results.dat',status='old')
        DO iter=1,state_count
          READ(60,*) dummy,ene(iter),mag(iter)
        END DO
        CLOSE(60)
        END SUBROUTINE read_in_val

        
        SUBROUTINE calc_averages
        IMPLICIT NONE
        CALL calc_ave_ene
        CALL calc_ave_mag
        CALL calc_heat_capacity
        CALL save_averages
        END SUBROUTINE calc_averages

        SUBROUTINE save_averages
        IMPLICIT NONE
        OPEN(61,file='average.dat')
        WRITE(61,*) tempera, ave_ene, ave_mag, heat_capacity
        CLOSE(61)
        END SUBROUTINE save_averages

        SUBROUTINE calc_ave_ene
        IMPLICIT NONE
        INTEGER         :: iter
        ave_ene = 0.0
        DO iter = 1,state_count 
          ave_ene = ave_ene + ene(iter) 
        END DO
        ave_ene = ave_ene / REAL(state_count)
        END SUBROUTINE calc_ave_ene


        SUBROUTINE calc_ave_mag
        IMPLICIT NONE
        INTEGER         :: iter
        ave_mag = 0.0
        DO iter =1, state_count
          ave_mag = ave_mag + mag(iter)
        END DO
        ave_mag = ave_mag / REAL(state_count)
        END SUBROUTINE calc_ave_mag


        SUBROUTINE calc_heat_capacity
        IMPLICIT NONE
        INTEGER         :: iter
        REAL            :: beta
        beta = 1.0 / (kb*(tempera**2))
        ene_sq = ene**2
        ave_ene_sq = 0.0
        DO iter = 1,state_count 
          ave_ene_sq = ave_ene_sq + ene_sq(iter) 
        END DO
        ave_ene_sq = ave_ene_sq / REAL(state_count)
        heat_capacity = ave_ene_sq - (ave_ene**2)
        heat_capacity = heat_capacity * beta
        END SUBROUTINE calc_heat_capacity

        END MODULE average

