PROGRAM lid
    USE subs
    USE types_vars
    USE matrixsolv
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT,R
    REAL(DP) :: H,L,tH,TL,tc,C
    INTEGER(DP)  :: I,J,K,N,counter
    INTEGER(DP) :: nx,ny,nt,n_1
    ! Mesh Variables ####################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:)    :: X, Y, T, Omega_magnitude,CC,DD
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: F,Fp,Fn, VARS
	REAL(DP) :: DELP,cfl
 	! Poisson Variables
 	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: A_P,X_v, B_v,fd
 	REAL(DP) :: BETA, Omega_SOR, tol
 	INTEGER(DP) :: N_cg,rows,cols
 	CHARACTER(LEN=20) :: pf,pf1,pf2,pf3
 	! VTE Variables
 	REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: Omega,Psi,V, U
 	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Psi_s,Omega_s,Am,Omega_n
 	REAL(DP) :: Re_L, nu
 	pf  = 'lid_ReL_50_Str.dat'
	pf1 = 'lid_ReL_50_Vrt.dat'
	pf2 = 'lid_ReL_50_U.dat'
	pf3 = 'lid_ReL_50_V.dat'
 	! Lid Driven Cavity Problem
 	! Advection Diffusion Equation solved with Roe Scheme
 	! Poisson Equation Solved with GS-SOR
 	Re_L = 1000
 	A = 0
 	B = 1
 	T1 = 0
 	T2 = 0.1
	n_1 = 1
 	to1 = 0.01
 	Omega_SOR = 1
 	dx = 0.01
 	dy = 0.01
 	dt = 0.001
 	h = DX
 	nu = 0.0001
 	beta = dy/dx
 	CALL gridgen(L,T1,T2,DT,T,nt)
 	CALL gridgen(L,A,B,DX,X,nx)
 	CALL gridgen(L,A,B,DY,Y,ny)
 	rows = nx
 	cols = ny
 	ALLOCATE(X_v(rows,cols),B_v(rows,cols),U(rows,cols,nt),FD(rows,cols))
 	ALLOCATE(Omega(nx,ny,nt),Psi(nx,ny,nt),Omega_s(nx,ny),V(rows,cols,nt))
	ALLOCATE(Omega_n(nx,ny),CC(NX),DD(NX))
 	! Build the 5-pt Scheme banded matrix to be used throughout
 	CALL P_MATRIX(Am,nx,ny,NX)
 	! Initialize the Initial Conditions
 	! Stream Function, Vorticity, Velocity 
 	! Set arrays to = 0 for Omega and looping without
 	! touching initial conditions 
 	! avoids the need to update them
 	! -- will need to renew omega / voriticty at top plate
 	Psi(:,:,:) = 0
 	Omega(:,:,:) = 0
 	U(:,:,:) = 0
 	U(:,1,:) = Re_L*nu/L
 	Psi(:,1,:) = Re_L*nu/L
	V(:,:,:) = 0
 	n = 0
	!Vorticity at t=0 @ the top plate boundary
	Omega(:,1,1) = -(Re_L*nu/L)/dy	
	CC(:) =0
	DD(:) =0
	! Stream function @ t=0 @ the top plate
	!Psi(i,2,1) = -Psi(i,1,1) + U(i,1,1)*dy
	write(*,*) 'check'
 	DO n = 1,nt
    	!Omega_s(:,:) = Omega(:,:,n)
    	CALL VTE_ROE_2D(Omega,nx,ny,DX,DY,DT,U,cc,V,dd,nu,n)
		b_v(:,:) = Omega(:,:,n+1)
		write(*,*) 'u :', U(int(nx*0.5),int(ny*0.5),n+1)
		write(*,*) 'Vorticity :', Omega(int(nx*0.5),2,n+1)
		write(*,*) 'Stream :', Omega(int(nx*0.5),2,n+1)
		! Call GS-SOR to solve stream fn elliptic PDE
		! Appearing as Poisson Equation 
		CALL SOR(Am,b_v,psi_s,h,nx,ny,x,y,beta,tol)
		! Store Poisson Solution In the Stream(Psi) Array
		Psi(:,:,n) = Psi_s(:,:)
		DO i=2,nx-1
			DO j=2,ny-1
			! calculate velocities from stream function
			U(i,j,n+1) = 0.5*(Psi(i+1,j,n)-Psi(i-1,j,n))/dy
			V(i,j,n+1) = -0.5*(Psi(i,j+1,n)-Psi(i,j-1,n))/dx
			END DO
		END DO
		! Top Wall - Vorticity BC
		Omega(:,1,n+1) = -((Re_L*nu/L)-U(:,1,n))/dy
		! Left wall - Vorticity BC
		Omega(1,:,n+1) = 2*(Psi(1,:,n)-Psi(2,:,n))/(dx*dx)
		! Right wall - Vorticity BC
		Omega(nx,:,n+1) = 2*(Psi(nx,:,n)-Psi(nx-1,:,n))/(dx*dx)
		! Bottom Wall - Vorticity BC
		Omega(:,ny,n+1) = 2*(Psi(:,ny,n)-Psi(:,ny-1,n))/(dy*dy)
		! convergence criteria
		! Check for convergence to steady state off the boundary
		! This should eventually become constant 
		DELP = 0
		Omega_n(:,:) = Omega(:,:,n+1)
		DELP = ABS(SUM(Omega_n(:,:)-Omega_s(:,:)))
		IF ( DELP .LE. 0.001 ) THEN
		write(*,*) 'time step: ', n, nt, 'total'
		EXIT
		END IF
	END DO
 T2 = n*dt
 A=A+1
 CALL gridgen(L,T1,T2,DT,T,nt)
 CALL tecplot(Psi,pf,x,y,t,nx,ny,nt,n_1,nx,n_1,ny,n_1,nt,n_1)
 CALL tecplot(Omega,pf1,x,y,t,nx,ny,nt,n_1,nx,n_1,ny,n_1,nt,n_1)
 CALL tecplot(U,pf2,x,y,t,nx,ny,nt,n_1,nx,n_1,ny,n_1,nt,n_1)
 CALL tecplot(V,pf3,x,y,t,nx,ny,nt,n_1,nx,n_1,ny,n_1,nt,n_1)
DEALLOCATE(X,Y,T,Omega,Psi,U)
END PROGRAM lid
