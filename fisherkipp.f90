PROGRAM FisherKpp
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    ! Time / Space Variables -----------------------------------------
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1,Q
    REAL(DP) :: DX,DT
    REAL(DP) :: V,R
    REAL(DP) :: L,tH,trange,tc,C
    INTEGER :: NPTSX,NPTST,I,J,NEQN,IC_CASE,N, PLOTT, PLOTT2
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: X, T, XX, QQ, AA, BB, CC
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: UIC, U, LWU, UEXACT, UBM, E, KTW
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI, M, D, KPPU
	CHARACTER(LEN=20) :: pf
	pf = 'kpp.dat'
    ! Allow user to pick which time step to plot-----------------------
    print*, "Enter a Graph Time to plot"
    read(*,*) tc
    T1=0
    T2=tc
    ALPHA=0.8
    KAPPA=0.1
    ! Domain Range , from A to B
    A=0
    B=1
    C=1
    ! KPP / IC Coeffiecients
    B1=0.5                    ! KPP coeffient
    N1=1                      ! sin(n1*pi*x)
    DX=0.1                  ! space step
    DT=0.01                   ! time step
    q = 1.0                   ! wave IC parameter
    n1 = 1
	CALL gridgen(L,A,B,DX,X,NPTSX)
    CALL gridgen(trange,T1,T2,DT,T,NPTST)
    CALL WAVE_IC(x,nptsx,UIC,N1,q)
    CALL KPP_FDE(alpha,dt,dx,kappa,nptsx,nptst,UIC,KPPU,B1)
	!SUBROUTINE KPP_FDE(a,dt,dx,kp,NPTSX,NPTST,UIC,ku,b)
    ALLOCATE(uexact(nptst,nptsx))
	uexact(1,:) = 0
    DO i=1,nptst
		DO j=2, nptsx
			uexact(i,j) = exp(t(i)*(b1-(kappa*n1*n1*pi*pi)))*(sin(n1*pi*x(j)))
        END DO
    END DO
    write(*,*)'x grid pts:',nptsx,'t grid pts',nptst
 	OPEN(1,file=pf,status='unknown')
 	DO i=1,nptsx 
		write(1,*) x(i), uexact(int(tc),i), KPPU(int(tc),i)
 	END DO
    A=-100
    B=100
    DX=1
    T1=0
    T2=tc
    dt = 0.001
	CALL gridgen(L,A,B,DX,X,NPTSX)
	CALL gridgen(trange,T1,T2,DT,T,NPTST)
    alpha=0.4
    kappa=1
    r = KAPPA*DT/(DX*DX)

    b1=2
    CALL KPP_TW(r,alpha,dt,NPTSX,NPTST,ktw,b1,C,x,t)
	
END PROGRAM FisherKpp
