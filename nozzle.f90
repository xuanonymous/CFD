PROGRAM nozzle
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT,R
    REAL(DP) :: H,L,tH,TL,tc,C
    INTEGER  :: nx,ny,nt,I,J,K,N,counter
    ! Mesh Variables ####################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:)    :: X, Y, T
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: UIC, U
	CHARACTER(LEN=20) :: pf
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: Mach,lambda,Tp,He,gamma,Rg
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: P,rho,u_scalar,v_scalar,sound
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: F,Fp,Fn, VARS
	! Nozzle Conditions ########################################
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: AREA,Et,P_N,P_N1,e_int,DELP1
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: H_source,dAdX
	REAL(DP) :: DELP,cfl
	DY = 0
	A =  0
	B =  10
	nx = 10000
	DX = (B-A)/nx
	nu = 0
	CFL = 0.00005
	CALL gridgen(L,A,B,DX,X,nx)
	ALLOCATE(u_scalar(nx),v_scalar(nx),He(nx),gamma(nx),P(nx),rho(nx),f(nx,3))
	ALLOCATE(Fp(nx,3),Fn(nx,3),u(nx,3),Tp(nx),Rg(nx),sound(nx),Mach(nx),DELP1(NX))
	ALLOCATE(AREA(NX),dAdX(NX,1),Et(nx),H_source(nx,3),e_int(nx),P_N1(nx),P_N(Nx))
	!################## Initial Conditions
	!##### Temperature 
	Tp(:) = 70
	!##### Universal Gas Constant = 1716 ft - lbf/slugs - R
	Rg(:) = 1716
	gamma(:) = 1.4
	!################## Initial Conditions at the Inlet ########################
	Mach(:) = 1.5
	He(:) = (gamma(:)/(gamma(:)-1))+0.5*(u_scalar(:)**2)
	! 1000 lbf/square foot
	P(:) = 1000
	! initial density = 0.00237 slugs / ft^3
	rho(:) = 0.00237
	gamma(:) = 1.4
		P(:) = 1000
		! initial density = 0.00237 slugs / ft^3
		rho(:) = 0.00237
	    SOUND(:) = SQRT(GAMMA(:)*P(:)/RHO(:))
	    u_scalar(:) = 1.5*sound(:)
	!MACH(2:nx) = U_SCALAR(2:nx)/SOUND(2:nx)
	! 1000 lbf/square foot
	! Boundary Conditions at the Outlet
	! Extrapolation 
	! rho(nx) = 2*rho(nx-1) - rho(nx-2) 
	! u_scalar(nx) = 2*u_scalar(nx-1) - u_scalar(nx-2) 
	! et(nx) = 2*et(nx-1) - et(nx-2) 
	v = (DT/DX)
	!DT = 0.4*DX/(ABS(MAXVAL(U_SCALAR))+MAXVAL(sound))
	WRITE(*,*) 'DT:', DT
	WRITE(*,*) 'DX:', DX
	WRITE(*,*) 'Van Leer FVS CFL: ', v	
	! Allocate the Nozzle Area Array
	CALL nozzle_area(AREA,X,nx)
	! Calculate dA/dX from O(2) first derivative
	CALL SOTPF(AREA,DX,dAdX,nx)
	!dAdX(:,1) = (1/(exp(0.8*x(:)-4)+exp(-0.8*x(:)+4)))*0.347
	! Specify subsonice velocity at the exit of the domain 390 ft/s
    CALL solvec(u,rho,u_scalar,P,gamma,nx)
    CALL fluxvec(u,f,gamma,nx)
    U(:,1) = U(:,1)*Area(:)
    U(:,2) = U(:,2)*Area(:)
    U(:,3) = U(:,3)*Area(:)
    f(:,1) = f(:,1)*Area(:)
    f(:,2) = f(:,2)*Area(:)
    f(:,3) = f(:,3)*Area(:)
    MACH(:) = U_SCALAR(:)/SOUND(:)
	e_int(:) = P(:)*(1-rho(:))/((gamma(:)-1)*rho(:))
	Et(:) = e_int(:) + 0.5*(u_scalar(:)**2)
	He(:) = (e_int(:) + P(:))/rho(:)
    P_N1(:) = P(:)
    CALL AUSM(f,u_scalar,Fp,Fn,rho,He,P,e_int,mach,gamma,sound,nx)
	DO
		dt = CFL*DX/MAXVAL(u_scalar+sound)
		v = dt/dx
		! Update the Source term then multiply it by dAdX(i)
		CALL source_vec(H_source,p,nx)
		DO counter=1,nx
			H_source(counter,1) = dAdX(counter,1)*H_source(counter,1)*dt
			H_source(counter,2) = dAdX(counter,1)*H_source(counter,1)*dt
			H_source(counter,3) = dAdX(counter,1)*H_source(counter,1)*dt
		END DO
		DO n=1,3
			DO i=2,nx-1
           		u(i,n) = u(i,n)-v*(Fp(i,n)+Fn(i+1,n)-(Fp(i-1,n)+Fn(i,n))) &
           				 +H_source(i,n)
			END DO
    	END DO
    	u(nx,:) = 2*u(nx-1,:) - u(nx-2,:)
    	u(1,:) = u(2,:)
		P_N(:) = P(:)
    	CALL DECODE_SOL(u,u_scalar,p,rho,sound,gamma,nx)
    	CALL fluxvec(u,f,gamma,nx)
    	! Extrapolate the other boundary conditions at the prims
    	rho(nx) = 2*rho(nx-1) - rho(nx-2)
        u_scalar(nx) = 2*u_scalar(nx-1) - u_scalar(nx-2)
        et(nx) = 2*et(nx-1) - et(nx-2)
    	MACH(:) = U_SCALAR(:)/SOUND(:)
		e_int(:) = P(:)*(1-rho(:))/((gamma(:)-1)*rho(:))
		Et(:) = e_int(:) + 0.5*(u_scalar(:)**2)
		He(:) = (e_int(:) + P(:))/rho(:)
    	! Check for convergence to steady state
    	delp = 0
	    v = (DT/DX)
    	CALL AUSM(f,u_scalar,Fp,Fn,rho,He,P,e_int,mach,gamma,sound,nx)
		DELP1(:) = P(:) - P_N(:)
		DELP = (SUM(DELP1))
		IF ( DELP .LE. 0.1 ) EXIT
	END DO
	pf = '1d_nozzle1.dat'
    OPEN(1,file=pf,status='unknown')
	DO i=1,nx
		write(1,*) x(i), ' ', rho(i), ' ', u_scalar(i), ' ', P(i), ' ', Mach(I)
	END DO
	END PROGRAM nozzle
