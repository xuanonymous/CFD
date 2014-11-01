PROGRAM PS7
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT
    REAL(DP) :: V,plotT,r,plotT2
    REAL(DP) :: H,L,tH,TL,tc,C
    INTEGER :: nx,ny,nt,I,J,NEQN,IC_CASE,N
    ! Mesh Variables ####################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: X, Y, T
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: UIC, U
    ! 3d Varaibles (x,y,t) ##############################################
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: IC_2D, U_2d
	! Thomas Algorithm Variables ########################################
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: XX, QQ, AA, BB, CC
    ! Fisher Variables ##################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) ::  LWU, UEXACT, UBM, E, KTW
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI, M, D, KPPU
	! Conjugate Gradient Method Variables ###############################
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Am, Am_t
	REAL(dp), ALLOCATABLE, DIMENSION(:) :: X_v, B_v
	! Representing v=1.0,0.6,0.3 respectively
	INTEGER :: N_cg
	! 2D Adv/Diff Variables #############################################
	! Where: F = b * u^2/2 + c * u
	! Ut + u Ux + v Uy = v(Uxx + Uyy)
	! u1 = u ; u2 = v
	REAL(DP) :: NU, U1, U2
	REAL(DP) :: B_be,C_be,D_be,E_be, mu
	CHARACTER(LEN=20) :: plot_file
	DX = 0.02
	DY = 0.02
	A = -PI
	B =  PI
	! H = 0.025, 0.0125
	DT = 0.01
	T1 = 0
	T2 = 1
	nu = 10**(-3)
	H = DT/DX
	WRITE(12,*) H
	U1 = 1.0
	U2 = 2.0
	B_be = 0
    C_be = U1
    D_be = 0
    E_be = U2
    ! v= 1
    CALL xgridgen(L,A,B,DX,X,nx)
    CALL xgridgen(L,A,B,DY,Y,ny)
    CALL tgridgen(T1,T2,DT,T,nt)
    ALLOCATE(IC_2D(nx,ny,1))
    DO j=1,nx
        DO i=1,ny
            IC_2D(j,i,1) = sin(2*pi*x(j))*sin(2*pi*y(i))
        END DO
    END DO 
    !CALL RIEMANN_IC(UIC,1,0,11,nt,nx,DX,DT)
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_be,C_be,D_be,E_be,nu,1)
    plot_file = 'soln_out.dat'
	!CALL AD_MAC_2D(U_2D,IC_2D,nu,nx,ny,nt,DX,DY,DT,B_be,c_be)
	call tecplot(U_2D,plot_file,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    !CALL surf(U_2D,pause=-1.0,terminal='windows',contour='both')
END PROGRAM PS7
