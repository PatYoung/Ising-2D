!二维伊辛模型能量最小 
!利用重要性抽样方法
!
!杨航
!
!最后修改：2018.05.10


!随机数种子初始化子程序，时随机数与系统时间有关
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
    INTEGER,PARAMETER :: p = 1
    INTEGER,PARAMETER :: q = 60
    REAL(KIND = 8) E,dE,starttime,endtime,a,b,c,d,Ei,Ef,stay
    REAL(KIND = 8) :: S(p + 2,q + 2)                !用一个62×62维的数组表示二维伊辛模型，多出的两个为考虑周期边界条件添加的
    OPEN(UNIT = 11, STATUS = 'REPLACE', POSITION = 'APPEND', FILE = '1.dat')
    
    CALL init_random_seed()
    WRITE(*,*) "n ->(n < 10000000)"
    READ(*,*) n
    !n = 10000                                       !n为迭代次数，可以修改，但是应该是有上限的
    T = 3.5                                         !温度T,可以修改
    E = 0                                           !总能量E
    stay = 1                                        !自旋改变的概率，与温度T和ΔE有关
    !DO i = 2,p + 1                                  !初始自旋赋值
    DO j = 2,q + 1
        CALL RANDOM_NUMBER(d)                   !取零到一的随机数
            
        IF (d < 1) THEN                         !这里随机取初始值，d < 1时，初始值均为-1，与均为1一样，修改d小于的数，影响自旋初始取值概率
            S(2,j) = - 1                        !满足条件，自旋初始取-1
        ELSE
            S(2,j) = 1                          !不满足条件，自旋初始值取1
        END IF
        !WRITE(*,*) "S(",i,j,")",S(i,j)         !这句格式写的有问题，写对了的话，应该可以看初始的取值情况，懒得改了
    END DO
    !END DO

    !DO i = 2,p + 1                                  !周期性边界条件
    S(2,1) = S(2,q + 1)
    S(2,q + 2) = S(2,2)
    !END DO

    !DO j = 2,q + 1
    !    S(1,j) = S(p + 1,j)                         !周期性边界条件
    !    S(p + 2,j) = S(2,j)
    !END DO
    
    !WRITE(*,*) S(i,j)

    !DO x = 2,p + 1                                  !求初始能量值，有没有问题我还没细想
    DO y = 2,q + 1
        E = E + S(2,y - 1) * S(2,y) + S(2,y + 1) * S(2,y)! + S(x,y - 1) * S(x,y) + S(x,y + 1) * S(x,y)
    END DO
    !END DO
    E = E/2                                         !这里除了个2，是我感觉上面重复计算了一次，只是感觉，也没细想
    WRITE(*,*) "00000",(E/(p * q))                  !输出初始能量
    WRITE(11,*) 0,(E/(p * q))
    CALL SYSTEM_CLOCK(COUNT = clock)
    starttime = clock                               !记录初始模拟时间

    !DO i = 2,p + 1
    !    DO j = 2,q + 1
    !        WRITE(*,*) S(i,j)
    !    END DO
    !END DO
    DO i = 1,n
        !CALL RANDOM_NUMBER(a)                       !用随机数a和b表示随意的在（2，61）取一个点，1和62是边界
        CALL RANDOM_NUMBER(b)
        !x = INT(p * a + 2)                          !把（0，1）对应到（2，61）
        y = INT(q * b + 2)
        x = 2
        m = S(x,y)
        S(x,y) = - m
        Ei = S(2,y - 1) * m + S(2,y + 1) * m! + S(x,y - 1) * m + S(x,y + 1) * m
        Ef = S(2,y - 1) * S(x,y) + S(2,y + 1) * S(x,y)! + S(x,y - 1) * S(x,y) + S(x,y + 1) * S(x,y)
        dE = Ef - Ei                                !计算自旋反向后能量变化，以计算自旋改变概率
        !IF (i < 10) THEN
        !    WRITE(*,*) x,y,S(x,y),dE,m
        !END IF
        stay = exp(-dE/T)                           !自旋改变概率
        
        CALL RANDOM_NUMBER(c)                       !利用随机数来满足这样的概率取值
        
        IF (c < stay) THEN                          !满足时，到下个态，自旋反向
            S(x,y) = - m
            E = E + dE
        ELSE                                        !不满足时，保持不变，能量不变
            S(x,y) = m
            E = E
        END IF
        
        IF (i > 0) THEN
            WRITE(11,*) i,(E/(p * q))               !输出到1.dat文件中，将11改为*是输出到命令行，这里IF不用写，我只想改下条件，取个截断
        END IF
    END DO
    !WRITE(*,*) stay,S(x,y)
    !DO i = 1,p
    !    DO j = 1,q
    !        WRITE(*,*) S(i,j)
    !    END DO
    !END DO
    CALL SYSTEM_CLOCK(COUNT = clock)
    endtime = clock                                 !记录结束模拟时间

    WRITE(*,*) "Time:",endtime - starttime, "ms"    !计算模拟时间，并输出

END PROGRAM
