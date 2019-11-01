        MODULE simulation
        INTEGER         :: loop
        REAL            :: tempera
        CONTAINS
        SUBROUTINE read_parameters
        OPEN(50,file='loop.dat',status='old')
        OPEN(51,file='temperature.dat',status='old')
        READ(50,*) loop
        READ(51,*) tempera
        CLOSE(50)
        CLOSE(51)
        END SUBROUTINE read_parameters
        END MODULE simulation
