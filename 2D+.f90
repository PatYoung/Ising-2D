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
    INTEGER i,j,n,m,x,y,clock
    INTEGER,PARAMETER :: p = 60
    INTEGER,PARAMETER :: q = 60
    REAL(KIND = 8) E,dE,starttime,endtime,a,b,c,d,Ei,Ef,stay,T,He,He2,PFC,ACC,Cv
    REAL(KIND = 8) :: S(p + 2,q + 2)
    OPEN(UNIT = 11, STATUS = 'REPLACE', POSITION = 'APPEND', FILE = '1.dat')
    OPEN(UNIT = 12, POSITION = 'APPEND', FILE = '2.dat')
    CALL init_random_seed()
    !WRITE(*,*) "n ->(n < 10000000)"
    !READ(*,*) n
    WRITE(*,*) "T ->"
    READ(*,*) T
    n = 500000
    !T = 1.5
    E = 0
    stay = 1
    DO i = 2,p + 1
        DO j = 2,q + 1
            CALL RANDOM_NUMBER(d)
            
            IF (d < 0.5) THEN
                S(i,j) = - 1
            ELSE
                S(i,j) = 1
            END IF
            !WRITE(*,*) "S(",i,j,")",S(i,j) 
        END DO
    END DO

    DO i = 2,p + 1
        S(i,1) = S(i,q + 1)
        S(i,q + 2) = S(i,2)
    END DO

    DO j = 2,q + 1
        S(1,j) = S(p + 1,j)
        S(p + 2,j) = S(2,j)
    END DO

    DO x = 2,p + 1
        DO y = 2,q + 1
            E = E + S(x - 1,y) * S(x,y) + S(x + 1,y) * S(x,y) + S(x,y - 1) * S(x,y) + S(x,y + 1) * S(x,y)
        END DO
    END DO
    E = E/2
    WRITE(*,*) "00000",E
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
        
        IF (stay > 1) THEN
            stay = 1
        END IF
        
        IF (i > 300000) THEN
            ACC = ACC + stay
            He = He + E
            He2 = He2 + E * E
        END IF
        !IF (i > (n - 1000)) THEN
        !    WRITE(*,*) i,E,ACC,He
        !END IF
    END DO
    
    He = He / n
    He2 = He2 /n
    PFC = He2 - He * He
    ACC = ACC / n
    Cv = PFC / (T * T)
    WRITE(*,*) "He = ",He,"He2 = ",He2
    WRITE(*,*) "ACC = ",ACC,"PFC = ",PFC 
    WRITE(*,*) "Cv = ",Cv
    WRITE(12,*) T,Cv
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
