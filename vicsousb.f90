PROGRAM sbr
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
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: U_U, U_V,B_x,C_x,B_y,C_y
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: RD
	REAL(DP):: mu,Ax,Bx,Ay,By,Ro
	CHARACTER(LEN=20) :: pf
	! viscous burgersN#################################################
    nt = 100
    nx = 100
    ny = 100
    
    A = -PI
    B = PI
    T1 = 0
    T2 = 1.5
    T_1 = 1.5
    nu = 10**(-3)
    
    DT = (T2-T1)/nt
    DX = (B-A)/nx
    DY = (B-A)/ny

    Ax = 0.1
    Bx = 0.6
    Ay = 0.1
    By = 0.6
    CALL xgridgen(L,A,B,DX,X,nx)
    CALL xgridgen(L,A,B,DY,Y,ny)
    CALL tgridgen(T1,T2,DT,T,nt)
    ALLOCATE(U_U(nx),U_V(ny),C_x(nx),C_y(ny),RD(nx,ny),B_x(nx),B_y(ny))
    ALLOCATE(IC_2D(nx,ny,1))
    DO j=1,nx
        DO i=1,ny
            IC_2D(j,i,1) = 0   
            C_x(j) = 1
            C_y(i) = 1 
        END DO
    END DO
    DO j=int((Ax-A)/DX),int((Bx-A)/DX)
        DO i=int((Ay-A)/DY),int((By-A)/DY)
            IC_2D(j,i,1) = 1
        END DO
    END DO 
    B_x = 0
    B_y = 0
    ! case 1 ######################################################
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'BELAX_1.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'BEMAC_1.dat'
    CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'BEROE_1.dat'
	CALL AD_ROE_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
	call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
	DEALLOCATE(U_2D)

	! case 2 ######################################################
    H = 0.025
    DT = 0.4*H
    CALL tgridgen(T1,T2,DT,T,nt)
	CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'BELAX_2.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'BEMAC_2.dat'
    CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'BEROE_2.dat'
    CALL AD_ROE_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    
	! case 3 ################################################################
    H = 0.0125
    DT = 0.4*H
    CALL tgridgen(T1,T2,DT,T,nt)
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'BELAX_3.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'BEMAC_3.dat'
    CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'BEROE_3.dat'
    CALL AD_ROE_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
END PROGRAM sbr
