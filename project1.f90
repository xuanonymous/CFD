PROGRAM project1
    USE SUBS
    USE TYPES_VARS
    ! Problem 1 REAL(DP), PARAMETER :: A = 0, B = 4, DX=0.1
    ! Problem 2 
    REAL(DP), PARAMETER :: A=0,B=4,dx=0.1,t1=0,t2=1,dt=0.1,c=2
    REAL(DP) :: H,L,tH,TL,V,N=1
    INTEGER :: NPTSX,NPTST, I
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: X, T
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: UIC, U
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI
    CHARACTER(LEN=20) :: pf,pf1
        CALL gridgen(L,A,B,DX,X,nx)
        CALL gridgen(L,T1,T2,DT,T,NT)
        CALL WAVE_IC(x,nx,uic,2,1)
        CALL LAX_WAVE(UIC,U,DX,DT,C,NX,NT)
        PF = 'soln_out.dat'
        pf1 = '1dlax.dat'
    OPEN(1,file=pf1,status='unknown')
    DO i=1,nx
        write(1,*) x(i), ' ', u(i,int(Nt*0.5)) 
    END DO
END PROGRAM project1
