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
	CFL = 0.01
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
    END DO 
	DO i=int(2.8/dx)+1,nx
		ux(i) = 390.75
		rho(i) = 0.00237
		P(i) = 1000
	END DO 	
	Tp(:) = 527.67/(1+(gma(:)-1)*0.5*Mc(:)**2)
	v = (DT/DX)
	! Allocate the Nozzle Area Array
	CALL nozzle_area(AREA,X,nx)
	! Calculate dA/dX from O(2) first derivative
	CALL SOTPF(AREA,DX,dAdX,nx)
	e_int(:) = P(:)/(rho(:)*(gma(:)-1))
	Et(:) = rho(:)*(e_int(:) + 0.5*(ux(:)**2))
    He(:) = (rho(:)*e_int(:) + P(:))/rho(:)
    P_N1(:) = P(:)
   	SOUND(:) = abs(SQRT(gma(:)*P(:)/RHO(:)))
	Mc(:) = ux(:)/sound(:)
    Mc(1) = 1.5
    CALL solvec(u,rho,ux,P,gma,nx)
    CALL fluxvec(u,f,gma,nx)
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
		DO n=1,3
			DO i=2,nx-1
				dt = CFL*DX/(ux(i)+sound(i))
				v = dt/dx
				H_source(i,2) = dAdX(i,1)*H_source(i,2)*dt
           		u(i,n) = u(i,n)-v*(Fp(i,n)+Fn(i+1,n)  &
           				-(Fp(i-1,n)+Fn(i,n)))/Area(i)+H_source(i,n)/Area(i)
			END DO
    	END DO
		P_N(:) = P(:)
    	CALL DECODE_SOL(u,ux,p,rho,sound,gma,nx)
	! Enforce back pressure boundary condition where pressure
    ! is enforced as a boundary condition on the outlet and all other
    ! information is taken from the interior
        P(1) = 1000
        ux(1) = 1153.0
		ux(nx) = 390.75
        P(nx) = 2459.6
    ! Extrapolate boundary conditions
        rho(nx) = 2*rho(nx-1) - rho(nx-2) 
        SOUND(:) = abs(SQRT(gma(:)*P(:)/RHO(:)))
		Mc(:) = ux(:)/sound(:)
		Mc(1) = 1.5
		Tp(:) = 527.67/(1+(gma(:)-1)*0.5*Mc(:)**2)
        e_int(:) = P(:)/(rho(:)*(gma(:)-1))
        e_int(nx) = 2*e_int(nx-1) - e_int(nx-2)
   		Et(:) = rho(:)*e_int(:) + 0.5*rho(:)*ux(:)*ux(:)
        !Et(nx) = 2*Et(nx-1) - Et(nx-2)
   		He(:) = (rho(:)*e_int(:) + P(:))/rho(:)
        !He(nx) = 2*He(nx-1) - He(nx-2)
	! Check for convergence to steady state###############################
    	delp = 0
    	DELP1(:) = (P_N(:)-P(:))
		DELP = (SUM(ABS(DELP1)))
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
		write(1,*) x(i), ' ', rho(i), &
					' ', ux(i), &
					' ', P(i), ' ', Mc(I), &
					' ', Area(i),' ',Tp(i)
	END DO
    DEALLOCATE(x,rho,ux,p,Mc,area,tp)
END PROGRAM nozzle2

