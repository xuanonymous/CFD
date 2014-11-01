PROGRAM nozzle2
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT,R
    REAL(DP) :: H,L,tH,TL,tc,C
    integer(dp)  :: nx,ny,nt,I,J,K,N,counter
    ! Mesh Variables ####################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:)    :: X, Y, T
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: UIC, U, U_new
	CHARACTER(LEN=20) :: pf
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: Mc,lambda,Tp,He,gma
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: P,rho,ux,v_scalar,sound,Rg
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: F,Fp,Fn, VARS
	! Nozzle Conditions #################################################
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: AREA,Et,P_N,P_N1,e_int,DELP1
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: H_source,dAdX
	REAL(DP) :: DELP,cfl,cv
	DY = 0
	A =  0
	B =  10
	nx = 100
	DX = (B-A)/nx
	nu = 0
	CFL = 0.4	
	CALL gridgen(L,A,B,DX,X,nx)
	ALLOCATE(ux(nx),v_scalar(nx),He(nx),gma(nx),P(nx))
	ALLOCATE(rho(nx),f(nx,3),Fp(nx,3),Fn(nx,3),u(nx,3),Tp(nx),Rg(nx))
	ALLOCATE(sound(nx),Mc(nx),DELP1(NX),H_source(nx,3),e_int(nx))
	ALLOCATE(AREA(NX),dAdX(NX,1),Et(nx),P_N1(nx),P_N(Nx),U_new(nx,3))
	!################## Initial Conditions ###############################
	gma(:) = 1.4
	Rg(:)=1716
	cv= (1/(gma(1)-1))*Rg(1)
	DO i=1,int(2.8/dx)
		ux(i) = 1153.0
		P(i) = 1000
		! initial density = 0.00237 slugs / ft^3
		rho(i) = 0.00237
		U(i,1) = 0.00237
		U(i,2) = 2.7323
		U(i,3) = 4075
    END DO  
	DO i=int(2.8/dx)+1,nx
		ux(i) = 390.75
		rho(i) = 0.00237
		P(i) = 1000
		U(i,1) = 0.00237
		U(i,2) = 0.92608
		U(i,3) = 2680.93
	END DO 	
	Tp(:) = 527.67/(1+(gma(:)-1)*0.5*Mc(:)**2)
	v = (DT/DX)
	! Allocate the Nozzle Area Array
	CALL nozzle_area(AREA,X,nx)
	! A(i) = 1.398 + 0.347*TANH(0.8*(X(i)-4.0))
	do i = 1,nx
		DADX(i,1) = 0.347*(1 - 0.8*(tanh(0.8*x(i)-4.0)**2))
	end do 
    e_int(:) = P(:)/(rho(:)*(gma(:)-1))
	Et(:) = e_int(:) + 0.5*(ux(:)**2)
   	He(:) = (Et(:) + P(:))/rho(:)
    P_N1(:) = P(:)
   	SOUND(:) = abs(SQRT(gma(:)*P(:)/RHO(:)))
	Mc(:) = ux(:)/sound(:)
	Mc(nx) = 2*Mc(nx-1)-Mc(nx-2)
	f(:,1) = rho(:)*ux(:)
	f(:,2) = rho(:)*ux(:)*ux(:) + P(:)
	f(:,3) = ux(:)*rho(:)*(Et(:)+P(:))
    CALL AUSM(f,ux,Fp,Fn,rho,He,P,e_int,Mc,gma,sound,nx)
	!CALL van_leer(f,ux,Mc,Fp,Fn,rho,sound,gma,nx)
	k=0
	! Time marching 1-D Euler FTCS Flux Splitting FVM #########
	DO
		k=k+1
		write(*,*) 'delp ', delp, &
					'T ', Tp(int(0.7*nx)), &
					'p ', p(int(0.7*nx)), &
					'M ', Mc(int(0.7*nx))
		! Update the Source term then multiply it by dAdX(i)
		CALL source_vec(H_source,p,nx)
		H_source(:,2) = dAdX(:,1)*H_source(:,2)*dt
		dt = CFL*DX/minval(abs(ux)+sound(:))
		v = dt/dx
		DO n=1,3
			DO i=2,nx-1
           		u(i,n) = u(i,n)-v*(Fp(i,n)+Fn(i+1,n)  &
           				-(Fp(i-1,n)+Fn(i,n)))+H_source(i,n)/Area(i)
			END DO
    	END DO
		P_N(:) = P(:)
    	CALL DECODE_SOL(u,ux,p,rho,sound,gma,nx)
	! Enforce back pressure boundary condition where pressure
    ! is enforced as a boundary condition on the outlet and all other
    ! information is taken from the interior
	! Inflow Boundary Conditions:
        P(1) = 1000
		rho(1) = 0.00237	
		ux(nx) = 1153.0
	! Outflow Boundary Conditions: 
		ux(nx) = 390.75
		P(nx) = 2459.6
        rho(nx) = 2*rho(nx-1) - rho(nx-2) 
    ! Update Primitive Variables
        SOUND(:) = abs(SQRT(gma(:)*P(:)/RHO(:)))
		Mc(:) = ux(:)/sound(:)
		Tp(:) = 527.67/(1+(gma(:)-1)*0.5*Mc(:)**2)
        e_int(:) = P(:)/(rho(:)*(gma(:)-1))
        e_int(nx) = 2*e_int(nx-1) - e_int(nx-2)
		Et(:) = e_int(:) + 0.5*(ux(:)**2)
   		He(:) = (Et(:) + P(:))/rho(:)
	! Check for convergence to steady state###############################
    	delp = 0
    	DELP1(:) = P(:)-P_N(:)
		DELP = abs(SUM(DELP1))
		IF (DELP .LE. 0.1) EXIT
	! ######################################################################
    	CALL solvec(u,rho,ux,P,gma,nx)
    	CALL fluxvec(u,f,gma,nx)
    	CALL AUSM(f,ux,Fp,Fn,rho,He,P,Et,Mc,gma,sound,nx)
		!CALL van_leer(f,ux,Mc,Fp,Fn,rho,sound,gma,nx)
	END DO
	
	pf = 'n2.dat'
    OPEN(1,file=pf,status='unknown')
	DO i=1,nx
		write(1,*) x(i),' ',rho(i),' ', ux(i),' ', P(i), ' ', Mc(I), &
					' ', dAdX(i,1),' ',Tp(i)
	END DO
    DEALLOCATE(x,rho,ux,p,Mc,area,tp)
END PROGRAM nozzle2

