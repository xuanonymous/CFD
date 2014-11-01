PROGRAM PS8
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT
    REAL(DP) :: plotT,r,plotT2
    REAL(DP) :: H,L,tH,TL,tc,C
    INTEGER  :: nx,ny,nt,I,J,K,NEQN,IC_CASE,N
    ! Mesh Variables ####################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:)    :: X, Y, T
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: UIC, U
    ! 2d Varaibles (x,y,t) ##############################################
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:):: IC_2D, U_2d
	! Thomas Algorithm Variables ########################################
    REAL(DP), ALLOCATABLE, DIMENSION(:)    :: XX, QQ, AA, BB, CC
    ! Fisher Variables ##################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)  ::  LWU, UEXACT, UBM, E, KTW
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: PHI, M, D, KPPU
	! Conjugate Gradient Method Variables ###############################
	REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: Am, Am_t
	REAL(dp), ALLOCATABLE, DIMENSION(:)    :: X_v, B_v
	integer :: N_cg
	! 2D Adv/Diff Variables #############################################
	! Where: F = b * u^2/2 + c * u
	! Ut + u Ux + v Uy = v(Uxx + Uyy)
	! u1 = u ; u2 = v
	REAL(DP) :: NU, U1, U2
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: B_x,C_x,B_y,C_y
	CHARACTER(LEN=20) :: pf
	! SOD Shock Tube : ###################################################
	! Rho = density
	! Mach = U / A : mach #
	! u_scalar, v_scalar, scalar components of velocity
	! Tp : temperature
	! gamma : 1.4
	REAL(DP) :: A_s
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: Mach,lambda,Tp,He,gamma
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: P,rho,u_scalar,v_scalar
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: F,Fp,Fn, VARS
	DX = 0.01
	DY = 0
	A = 0
	B = 1
	T1 = 0
	T2 = 1
	nu = 0
	CALL gridgen(L,A,B,DX,X,nx)
    ALLOCATE(IC_2D(nx,1,1),B_x(nx),C_x(nx),B_y(ny),C_y(ny),F(nx,3))
	ALLOCATE(P(NX),rho(nx),u_scalar(nx),v_scalar(nx),gamma(nx),mach(nx),he(nx))
	ALLOCATE(U(nx,3),Tp(nx))
	! Shock Tube Initial Conditions########################################
	! Speed of sound
	gamma(:) = 1.4
	A_s = SQRT(MAXVAL(GAMMA)*8.3145*293.15)
    DO j=1,nx
        u_scalar(j) = 0
        v_scalar(j) = 0
		Mach(j) = u_scalar(j)/A_s
		He(j) = (gamma(j)/(gamma(j)-1))+0.5*(u_scalar(j)**2+v_scalar(j)**2)
    END DO 
	DO i=1,int(nx*0.5)
		P(i) = 1
		rho(i) = 1.0 
	END DO
	DO i=int(nx*0.5),nx
		P(i)   = 0.1
		rho(i) = 0.125
	END DO
	! Temperature, Kelvins correpsonding with universal
	! gas constant of 8.3145 J/mol K
	Tp(:) = 305
	DT = DX/A_s
    CALL gridgen(L,T1,T2,DT,T,nt)
	CALL PKG(rho,u_scalar,v_scalar,P,Tp,He,nx,gamma,Vars)
	CALL EULER1D_MAC(U,nx,nt,DX,DT,vars)
    !#####################################################################
    !DEALLOCATE(U_2D)
    
    !pf = 'AD_MAC_2D_1.dat'
    !CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    !call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    !DEALLOCATE(U_2D)
    
    !pf = 'AD_ROE_2D_1.dat'
	!CALL AD_ROE_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
	!call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
	!DEALLOCATE(U_2D,IC_2D)
END PROGRAM PS8
