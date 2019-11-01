        PROGRAM MAIN
        IMPLICIT NONE
        INTEGER        :: iter
        INTEGER        :: start, loops
        REAL           :: incre = 0.02
        REAL           :: tempera
        CHARACTER(LEN=30)      :: frame
        CHARACTER(LEN=30)   :: file_out
        start =  0
        loops = 300
        DO iter=start, loops
          tempera = REAL(start)+incre*REAL(iter)
          WRITE(*,*) tempera
          WRITE(*,'(F5.2)') tempera
!          WRITE(*,*) REAL(start)+incre*REAL(iter)
!          WRITE(*,*) iter
!          WRITE(frame,"(I3)") INT(tempera*100)
          WRITE(frame,"(F4.2)") tempera
!          frame = adjustl(frame)
          WRITE(*,*) frame
          file_out = trim(frame)//"_kelvins.dat"
          OPEN(21,file=trim(file_out),status='replace')
          WRITE(21,"(F4.2)") tempera
          CLOSE(21)
        END DO
        END PROGRAM MAIN
