PROGRAM PS7
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT
    REAL(DP) :: V,plotT,r,plotT2
    REAL(DP) :: H,L,tH,TL,tc,C
    integer :: nx,ny,nt,I,J,NEQN,IC_CASE,N
    ! Mesh Variables ####################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: X, Y, T
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: UIC, U
    ! 2d Varaibles (x,y,t) ##############################################
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: IC_2D, U_2d
	! Thomas Algorithm Variables ########################################
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: XX, QQ, AA, BB, CC
    ! Fisher Variables ##################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) ::  LWU, UEXACT, UBM, E, KTW
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI, M, D, KPPU
	! Conjugate Gradient Method Variables ###############################
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Am, Am_t
	REAL(dp), ALLOCATABLE, DIMENSION(:) :: X_v, B_v
	integer :: N_cg
	! 2D Adv/Diff Variables #############################################
	! Where: F = b * u^2/2 + c * u
	! Ut + u Ux + v Uy = v(Uxx + Uyy)
	! u1 = u ; u2 = v
	REAL(DP) :: NU, U1, U2
	REAL(DP) :: B_x,C_x,B_y,C_y, mu 
	CHARACTER(LEN=20) :: pf
	! Problem 1 : ########################################################
	H = 0.05
	DT = 0.4*H
	DX = 0.05
	DY = 0.05
	A = -PI
	B =  PI
	! H = 0.025, 0.0125
    H = 0.05
	T1 = 0
	T2 = 1
	nu = 10**(-3)
	U1 = 1.0
    U2 = 2.0
    B_x = 0
    C_x = U1
    B_y = 0
    C_y = U2
    CALL xgridgen(L,A,B,DX,X,nx)
    CALL xgridgen(L,A,B,DY,Y,ny)
    CALL tgridgen(T1,T2,DT,T,nt)
    ALLOCATE(IC_2D(nx,ny,1))
    DO j=1,nx
        DO i=1,ny
            IC_2D(j,i,1) = sin(x(j))*sin(y(i))
        END DO
    END DO 
    ! case 1 #############################################################
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'AD_LAX_2D_1.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_MAC_2D_1.dat'
    CALL AD_MAC_2D(U_2d,IC_2D,nu,nx,ny,nt,DX,DY,DT,B_x,C_x,B_y,C_y)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_ROE_2D_1.dat'
	CALL AD_ROE_2D(U_2d,IC_2D,nu,nx,ny,nt,DX,DY,DT,B_x,C_x,B_y,C_y)
	call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
	DEALLOCATE(U_2D,IC_2D)
	
	! case 2 ##############################################################
	H = 0.0025
    DT = 0.4*H
    CALL xgridgen(L,A,B,DX,X,nx)
    CALL xgridgen(L,A,B,DY,Y,ny)
    CALL tgridgen(T1,T2,DT,T,nt)
    ALLOCATE(IC_2D(nx,ny,1))
    DO j=1,nx
        DO i=1,ny
            IC_2D(j,i,1) = sin(x(j))*sin(y(i))
        END DO
    END DO 
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'AD_LAX_2D_2.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_MAC_2D_2.dat'
    CALL AD_MAC_2D(U_2d,IC_2D,nu,nx,ny,nt,DX,DY,DT,B_x,C_x,B_y,C_y)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_ROE_2D_2.dat'
    CALL AD_ROE_2D(U_2d,IC_2D,nu,nx,ny,nt,DX,DY,DT,B_x,C_x,B_y,C_y)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D,IC_2D)
    
    ! case 3 ################################################################
    H = 0.0125
    DT = 0.4*H
    CALL xgridgen(L,A,B,DX,X,nx)
    CALL xgridgen(L,A,B,DY,Y,ny)
    CALL tgridgen(T1,T2,DT,T,nt)
    ALLOCATE(IC_2D(nx,ny,1))
    DO j=1,nx
        DO i=1,ny
            IC_2D(j,i,1) = sin(x(j))*sin(y(i))
        END DO
    END DO 
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'AD_LAX_2D_3.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_MAC_2D_3.dat'
    CALL AD_MAC_2D(U_2d,IC_2D,nu,nx,ny,nt,DX,DY,DT,B_x,C_x,B_y,C_y)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_ROE_2D_3.dat'
    CALL AD_ROE_2D(U_2d,IC_2D,nu,nx,ny,nt,DX,DY,DT,B_x,C_x,B_y,C_y)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D,IC_2D)
END PROGRAM PS7
