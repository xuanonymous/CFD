PROGRAM PS5
    USE SUBS
    USE TYPES_VARS
	USE MATRIXSOLV
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DT
    REAL(DP) :: V,plotT,r,plotT2
    REAL(DP) :: H,L,tH,TL,tc,C
    INTEGER :: NPTSX,NPTST,I,J,NEQN,IC_CASE,N
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: X, T, XX, QQ, AA, BB, CC, Y
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: UIC, U, LWU, UEXACT, UBM, E, KTW
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI, M, D, KPPU
	! Initialization variables for Conjugate Gradient Method
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Am, Am_t
	REAL(dp), ALLOCATABLE, DIMENSION(:) :: X_v, B_v, h_pos,fd
	! Poisson Variables
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: A_P
	REAL(DP) :: BETA, OMEGA, tol
	INTEGER :: N_cg,rows,cols
	CHARACTER(LEN=20) :: pf
	A = -1
	B = 1
	to1 = 0.0001
	omega = 1.3
	dx = 0.1
	dy = dx
	beta = dy/dx
	CALL gridgen(L,A,B,DX,X,nx)
	CALL gridgen(L,A,B,DX,Y,ny)
	rows = nx
	cols = ny
	h = DX
	write(*,*) rows, cols
	ALLOCATE(X_v(rows),B_v(rows),U(cols,rows),FD(rows))
	CALL P_MATRIX(Am,cols,rows,NX)
	
	x_v(:) = 1

	! #########################################
	b_v(:) = 1
	CALL SOR(Am,X_V,B_V,u,h,ROWS,COLS,beta,omega,tol)
	pf = 'poisson.dat'
	OPEN(1,file=pf,status='unknown')
	DO i=1,rows
		write(1,*) x(i), y(i), u(i,:)
	END DO
	!CALL JI(Am,X_V,B_V)
   	!CALL space_initial(L,A,B,DX,H,NPTSX)
    !CALL time_initial(TL,T1,T2,DT,tH,NPTST)
    !CALL xgridgen(A,H,NPTSX,X)
    !CALL tgridgen(T1,TH,NPTST,T)
END PROGRAM PS5
