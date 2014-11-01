PROGRAM PS6
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2

    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DT
    REAL(DP) :: V,plotT,r,plotT2
    REAL(DP) :: H,L,tH,TL,tc,C
    INTEGER(DP) :: nx,nt,I,J,NEQN,IC_CASE,N,ux,uy,x1
    ! Mesh Variables ###################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: X, T
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: UIC, U
	! Thomas Algorithm Variables ###################################
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: XX, QQ, AA, BB, CC
    ! Fisher Variables #################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) ::  LWU, UEXACT, UBM, E, KTW
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI, M, D, KPPU
    ! Conjugate Gradient Method Variables ##########################
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Am, Am_t
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: X_v, B_v
    ! Problem 4.49 Storage Variables with underscore  
    ! Representing v=1.0,0.6,0.3 respectively
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: U_1,U_6,U_3,UBR_1,UBR_6,UBR_3
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: UBL,UBR,UBM_1,UBM_6,UBM_3
    INTEGER(DP) :: N_cg, d_vis, x_1,x_2,t_1,t_2
    ! Burger's Equation Variables ###################################
    ! Where: F = b * u^2/2 + c * u
    	REAL(DP) :: B_be,C_be, mu
    	CHARACTER(LEN=20) :: p1, p2, tec
	d_vis=0
   	ux = 0
	uy = 1
	x1 = 2
	B_be = 1
	C_be = 0
	DT = 0.06
	DX = 0.1
	A = 0
	B = 5
	T1 = 0
	T2 = 1
	mu = 10**(-6)
    ! CALL BE_LINEAR(U,X,T,mu,nt,nx,DX,DT,c)
    ! v= 1
    DT = 0.1
    DX = 0.1
    CALL gridgen(L,A,B,DX,X,nx)
    CALL gridgen(L,T1,T2,DT,T,nt)
    ALLOCATE(U_1(nx,nt),UIC(nx,nt),UBM_1(nx,nt))
    CALL RIEMANN_IC(UIC,Ux,Uy,x1,nx)
    CALL BURGERS_LAX(U_1,UIC,nx,nt,DX,DT,b_be,c_be)
    CALL BURGERS_MAC(UBM_1,UIC,nx,nt,DX,DT,b_be,c_be)
    OPEN(UNIT=12,FILE="Lax11.TXT",ACTION="WRITE")
    DO i=1,nx
        WRITE(12,*) X(i),',',U_1(i,2)
    END DO
    OPEN(UNIT=12,FILE="Mac1.TXT",ACTION="WRITE")
    DO i=1,nx
    	WRITE(12,*) X(i),',',UBM_1(i,2)
    END DO
    ! v= 0.3
    DT = 0.03
    DX = 0.1
    CALL gridgen(L,A,B,DX,X,nx)
    CALL gridgen(L,T1,T2,DT,T,nt)
    DEALLOCATE(UBM_1,U_1,UIC)
    ALLOCATE(U_3(nx,nt),UIC(nx,nt),UBM_3(nx,nt))
    CALL RIEMANN_IC(UIC,Ux,Uy,x1,nx)
    CALL BURGERS_LAX(U_3,UIC,nx,nt,DX,DT,b_be,c_be)
    CALL BURGERS_MAC(UBM_3,UIC,nx,nt,DX,DT,b_be,c_be)
    OPEN(UNIT=12,FILE="Lax13.TXT",ACTION="WRITE")
    DO i=1,nx
        WRITE(12,*) X(i),',',U_3(i,10)
    END DO
    OPEN(UNIT=12,FILE="Mac13.TXT",ACTION="WRITE")
    DO i=1,nx
        WRITE(12,*) X(i),',',UBM_3(i,10)
    END DO
    ! v= 0.6
    DT = 0.06
    DX = 0.1
    CALL gridgen(L,A,B,DX,X,nx)
    CALL gridgen(L,T1,T2,DT,T,nt)
    CALL RIEMANN_IC(UIC,Ux,Uy,x1,nx)
    CALL BURGERS_LAX(U_6,UIC,nx,nt,DX,DT,b_be,c_be)
    CALL BURGERS_MAC(UBM_6,UIC,nx,nt,DX,DT,b_be,c_be)
    !PRINT these----------------------------------------------------
    OPEN(UNIT=12,FILE="Lax16.TXT",ACTION="WRITE")
        DO i=1,nx
            WRITE(12,*) X(i),',',U_6(i,int(nt*0.1))
    END DO
    OPEN(UNIT=12,FILE="Mac16.TXT",ACTION="WRITE")
        DO i=1,nx
            WRITE(12,*) X(i),',',UBM_6(i,int(nt*0.1))
    END DO
    ! CALL BE_LINEAR(U,X,T,mu,nt,nx,DX,DT,c)
    ! v= 1 ############################# SIN IC #####################
    DT = 0.001
    DX = 0.001
    A = -4*pi
    B = 4*pi
    T1 = 0
    T2 = 1
    DEALLOCATE(UIC,X,T)
    CALL gridgen(L,A,B,DX,X,nx)
    CALL gridgen(L,T1,T2,DT,T,nt)
    ALLOCATE(UIC(nx,nt),UBR_1(nx,nt),U_1(nx,nt),UBM_1(nx,nt))
    CALL BE_SIN_IC(UIC,nx,nt,x)
    CALL BURGERS_ROE(UBR_1,UIC,nx,nt,DX,DT,b_be,c_be)
    CALL BURGERS_LAX(U_1,UIC,nx,nt,DX,DT,b_be,c_be)
    CALL BURGERS_MAC(UBM_1,UIC,nx,nt,DX,DT,b_be,c_be)
	nt = int(nt*0.1)
    OPEN(UNIT=12,FILE="Roe.TXT",ACTION="WRITE")
    DO i=1,nx
      WRITE(12,*) X(i),',',UBR_1(i,int(0.1*nt)),',',U_1(i,int(0.1*nt)),',',UBM_1(i,int(0.1*nt))
    END DO
    tec = 'roe_tec.dat'
    t_1 = 1
    t_2 = nt-4
    x_1 = 1
    x_2 = nx-4
    CALL tecplot2(UBR_1,tec,t,x,nt,nx,t_1,t_2,x_1,x_2,t_2)
END PROGRAM PS6
