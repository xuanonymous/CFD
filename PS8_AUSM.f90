PROGRAM PS8_AUSM
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT,R
    REAL(DP) :: H,L,tH,TL,tc,C
    INTEGER  :: nx,ny,nt,I,J,K,N
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
	! Re = ideal gas constant
	! Tp = temperature
	! Mach = U / A : mach #
	! u_scalar, v_scalar, scalar components of velocity
	! Tp : temperature
	! gamma : 1.4
	REAL(DP) :: A_s
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: Mach,lambda,Tp,He,gamma,Rg,e_int,Et
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: P,rho,u_scalar,v_scalar,sound
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: F,Fp,Fn,VARS,uo
	! Shock Tube Initial Conditions########################################
	DY = 0
	A = 0
	B = 1
	nx = 100
	CFL = 0.00001
	DX = (B-A)/nx
	T1 = 0
	T2 = 0.1452
	nu = 0
	CALL gridgen(L,A,B,DX,X,nx)
	ALLOCATE(u_scalar(nx),v_scalar(nx),He(nx),gamma(nx),P(nx),rho(nx),f(nx,3))
	ALLOCATE(Fp(nx,3),Fn(nx,3),u(nx,3),Tp(nx),Rg(nx),sound(nx),Mach(nx),uo(Nx,3))
	ALLOCATE(Et(nx),e_int(nx))
    DO j=1,nx
        u_scalar(j) = 0
        v_scalar(j) = 0
    END DO 
	DO i=int(nx*0.5)+1,nx
		P(i) = 0.1
		rho(i) = 0.125
	END DO
	DO i=1,int(nx*0.5)
		P(i)   = 1.0
		rho(i) = 1.0
	END DO
	gamma(:) = 1.4
	SOUND(:) = SQRT(GAMMA(:)*P(:)/RHO(:))
	MACH(:) = U_SCALAR(:)/SOUND(:)
	e_int(:) = P(:)*(1-rho(:))/((gamma(:)-1)*rho(:))
	Et(:) = e_int(:) + 0.5*(u_scalar(:)**2)
	He(:) = (e_int(:) + P(:))/rho(:)
    CALL solvec(u,rho,u_scalar,P,gamma,nx)
    CALL fluxvec(u,f,gamma,nx)
    CALL DECODE_SOL(u,u_scalar,p,rho,sound,gamma,nx)
	MACH(:) = U_SCALAR(:)/SOUND(:)
	dt=CFL*dx
	v = (DT/DX)
	e_int(:) = P(:)*(1-rho(:))/((gamma(:)-1)*rho(:))
    CALL gridgen(L,T1,T2,DT,T,nt)
    	WRITE(*,*) 'AUSM FVS CFL: ', v
		WRITE(*,*) 'iteration :', k
		write(*,*) 'total iterations: ', nt
	DO k=1,nt
		CALL AUSM(f,u_scalar,Fp,Fn,rho,He,P,e_int,mach,gamma,sound,nx)
		DO n=1,3
			DO i=2,nx-1
           		u(i,n) = u(i,n)-v*(Fp(i,n)+Fn(i+1,n)-(Fp(i-1,n)+Fn(i,n)))
			END DO
    	END DO
    	U(nx,:) = U(nx-1,:)
		U(1,:) = U(2,:)
    	CALL DECODE_SOL(u,u_scalar,p,rho,sound,gamma,nx)
    	CALL fluxvec(u,f,gamma,nx)
		MACH(:) = U_SCALAR(:)/SOUND(:)
		e_int(:) = P(:)*(1-rho(:))/((gamma(:)-1)*rho(:))
		Et(:) = e_int(:) + 0.5*(u_scalar(:)**2)
		He(:) = (e_int(:) + P(:))/rho(:)
		write(*,*) 'ux center:', u_scalar(int(nx*0.5))
		write(*,*) 'ux boundary:', u_scalar(int(nx-10))	
		write(*,*) 'Mach #', mach(int(nx*0.5))
		write(*,*) 'Pressure ', P(int(nx*0.5))
		write(*,*) 'rho ', rho(int(nx*0.5))
		write(*,*) 'sound ', sound(int(nx*0.5))
		WRITE(*,*) 'DT:', DT
		WRITE(*,*) 'DX:', DX
		DEALLOCATE(fp,fn)
	END DO
	pf = '1d_sod_ausm.dat'
    OPEN(1,file=pf,status='unknown')
	DO i=1,nx
		write(1,*) x(i), ' ', rho(i), ' ', u_scalar(i), ' ', P(i), ' ' 
	END DO
    !#####################################################################
    !DEALLOCATE(U_2D)
    !CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    !call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    !DEALLOCATE(U_2D)
END PROGRAM PS8_AUSM
