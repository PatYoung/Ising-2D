SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, ALLOCATABLE :: seed(:)
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    CALL SYSTEM_CLOCK(COUNT = clock)
    seed = clock + 37 * (/(i - 1, i = 1, n)/)
    CALL RANDOM_SEED(put = seed)
    DEALLOCATE(seed)
END SUBROUTINE init_random_seed



PROGRAM Ising
    IMPLICIT NONE
    INTEGER i,j,T,n,m,x,y,clock
    INTEGER,PARAMETER :: p = 60
    INTEGER,PARAMETER :: q = 60
    REAL(KIND = 8) E,dE,starttime,endtime,a,b,c,Ei,Ef,stay
    REAL(KIND = 8) :: S(p + 2,q + 2)
    OPEN(UNIT = 11, POSITION = 'APPEND', FILE = '1.dat')    
    
    CALL init_random_seed()
    !WRITE(*,*) "n ->(n < 10000000)"
    !READ(*,*) n
    n = 100000
    T = 3.5
    E = 0
    stay = 1
    DO i = 1,p
       DO j = 1,q
            S(i,j) = 1
        END DO
    END DO

    DO i = 2,p
        S(i,1) = S(i,q + 1)
        S(i,q + 2) = S(i,2)
    END DO

    DO j = 2,q
        S(1,j) = S(p + 1,j)
        S(p + 2,j) = S(2,j)
    END DO

    CALL SYSTEM_CLOCK(COUNT = clock)
    starttime = clock

    !DO i = 1,p
    !    DO j = 1,q
    !        WRITE(*,*) S(i,j)
    !    END DO
    !END DO
    DO i = 1,n
        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(b)
        x = INT(p * a + 2)
        y = INT(q * b + 2)
        m = S(x,y)
        S(x,y) = - m
        Ei = S(x - 1,y) * (- S(x,y)) + S(x + 1,y) * (- S(x,y)) + S(x,y - 1) * (- S(x,y)) + S(x,y + 1) * (- S(x,y))
        Ef = S(x - 1,y) * S(x,y) + S(x + 1,y) * S(x,y) + S(x,y - 1) * S(x,y) + S(x,y + 1) * S(x,y)
        dE = Ef - Ei
        !E = E + dE
        !WRITE(*,*) x,y,S(x,y),dE,m
    
        stay = exp(-dE/T)
        CALL RANDOM_NUMBER(c)
        IF (c < stay) THEN
            S(x,y) = - m
            E = E + dE
        ELSE
            S(x,y) = m
            E = E
        END IF
        IF (i > 0) THEN
            !E = E / 10000
            WRITE(11,*) i,E
        END IF
    END DO
    !WRITE(*,*) stay,S(x,y)
    !DO i = 1,p
    !    DO j = 1,q
    !        WRITE(*,*) S(i,j)
    !    END DO
    !END DO
    CALL SYSTEM_CLOCK(COUNT = clock)
    endtime = clock

    WRITE(*,*) "Time:",endtime - starttime, "ms"

END PROGRAM
