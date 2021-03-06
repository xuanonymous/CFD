PROGRAM lid
    USE subs
    USE types_vars
    USE matrixsolv
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT,R
    REAL(DP) :: H,L,tH,TL,tc,C
    INTEGER  :: nx,ny,nt,I,J,K,N,counter
    ! Mesh Variables ####################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:)    :: X, Y, T, Omega_magnitude
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: F,Fp,Fn, VARS
	REAL(DP) :: DELP,cfl
 	! Poisson Variables
 	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: A_P,X_v, B_v,fd
 	REAL(DP) :: BETA, Omega_SOR, tol
 	INTEGER :: N_cg,rows,cols
 	CHARACTER(LEN=20) :: pf,pf1,pf2,pf3
 	! VTE Variables
 	REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: Omega,Psi,V, U
 	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Psi_s,Omega_s,Am,Omega_n
 	REAL(DP) :: Re_L, nu
 	pf = 'lid_ReL_50_Str.dat'
	pf1 = 'lid_ReL_50_Vrt.dat'
	pf2 = 'lid_ReL_50_U.dat'
	pf3 = 'lid_ReL_50_V.dat'
 	! Lid Driven Cavity Problem
 	! Advection Diffusion Equation solved with FTCS
 	! Poisson Equation Solved with GS-SOR
 	Re_L = 1000
 	A = 0
 	B = 8
 	T1 = 0
 	T2 = 1
 	to1 = 0.01
 	Omega_SOR = 1
 	dx = 0.01
 	dy = dx
 	dt = 0.01
 	h = DX
 	nu = 0.001
 	beta = dy/dx
 	CALL gridgen(L,T1,T2,DT,T,nt)
 	CALL gridgen(L,A,B,DX,X,nx)
 	CALL gridgen(L,A,B,DY,Y,ny)
 	rows = nx
 	cols = ny
 	ALLOCATE(X_v(rows,cols),B_v(rows,cols),U(rows,cols,nt),FD(rows,cols))
 	ALLOCATE(Omega(nx,ny,nt),Psi(nx,ny,nt),Omega_s(nx,ny),V(rows,cols,nt))
	ALLOCATE(Omega_n(nx,ny))
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
 	U(:,1,1) = Re_L*nu/L
 	U(1,:,1) = Re_L*nu/L
 	Psi(:,1,:) = Re_L*nu/L
	V(:,:,:) = 0
 	n = 0
	!Vorticity at t=0 @ the top plate boundary
	Omega(:,1,1) = -(Re_L*nu/L)/dy	

	! Stream function @ t=0 @ the top plate
	!Psi(i,2,1) = -Psi(i,1,1) + U(i,1,1)*dy

 	DO n=1,nt-1
    Omega_s(:,:) = Omega(:,:,n)
		DO i=2,nx-1
			DO j=2,ny-1
			! solve 2d, adv-diff VTE 
		    	Omega(i,j,n+1)  =	&	 
   				-dt*0.5*u(i,j,n)*(Omega(i+1,j,n)-Omega(i-1,j,n))/Dx    	&
				-dt*0.5*v(i,j,n)*(Omega(i,j+1,n)-Omega(i,j-1,n))/Dy    	&
				+ nu*DT*(Omega(i,j+1,n)					                &
				-2*Omega(i,j,n)+Omega(i-1,j,n))/(dy*dy)                 &
				+ nu*DT*(Omega(i+1,j,n)					                &
				-2*Omega(i,j,n)+Omega(i,j-1,n))/(dx*dx)                 &
				+ Omega(i,j,n) 
				! store vorticity data in b_v to send to
				! the gs-sor poisson solver
				b_v(i,j) = Omega(i,j,n+1)
		    END DO
		END DO
		write(*,*) 'u :', U(int(nx*0.5),int(ny*0.5),n+1)
		write(*,*) 'Vorticity :', Omega(int(nx*0.5),2,n+1)
		write(*,*) 'Stream :', Omega(int(nx*0.5),2,n+1)
		! Call GS-SOR to solve stream fn elliptic PDE
		! Appearing as Poisson Equation 
		CALL SOR(Am,b_v,psi_s,h,nx,ny,beta,Omega_SOR,tol)
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
		Omega(:,1,n+1) = -((Re_L*nu/L)-U(:,1,n+1))/dy
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
		DELP = ABS(SUM(Omega_n(:,:)-Omega_s(:,:)))/dt
		IF ( DELP .LE. 0.001 ) EXIT
	END DO
 T2 = n*dt
 CALL gridgen(L,T1,T2,DT,T,nt)
 write(*,*) 'Total Iterations :', n
 OPEN(1,file=pf,status='unknown')
 	DO i=1,rows 
		write(1,*) Psi(i,:,nt-1) 
 	END DO
 OPEN(1,file=pf1,status='unknown')
 	DO i=1,rows 
		write(1,*) Omega(i,:,nt-1) 
 	END DO
 OPEN(1,file=pf2,status='unknown')
 	DO i=1,rows 
		write(1,*) U(i,:,nt-1) 
 	END DO
 OPEN(1,file=pf3,status='unknown')
 	DO i=1,rows 
		write(1,*) V(i,:,nt-1) 
 	END DO
 DEALLOCATE(X,Y,T,Omega,Psi,U)
 END PROGRAM lid
