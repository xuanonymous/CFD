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
	! Error Variables
	REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: E1,E2,E3,E4,E5,E6
	CHARACTER(LEN=20) :: pf
	! SOLID BODY ROTATION#################################################
    H = 0.05
    DT = 0.4*H
    DX = 0.04
    DY = 0.04
    A = -PI
    B = PI
    T1 = 0
    T2 = 0.25
    nu = 0
   ! Initial Conditions for Solid Body Rotation
    Ax = 0
    Bx = 0.5
    Ay = 1
    By = 1.25
    Ro = 0.2
    CALL xgridgen(L,A,B,DX,X,nx)
    CALL xgridgen(L,A,B,DY,Y,ny)
    CALL tgridgen(T1,T2,DT,T,nt)
    ALLOCATE(U_U(nx),U_V(ny),C_x(nx),C_y(ny),RD(nx,ny),B_x(nx),B_y(ny))
    ALLOCATE(E1(nx,ny,nt),E2(nx,ny,nt),E3(nx,ny,nt),E4(nx,ny,nt),E5(nx,ny,nt),E6(nx,ny,nt))
    ALLOCATE(IC_2D(nx,ny,1))
    DO j=1,nx
        DO i=1,ny
            RD(j,i) = 0
            IC_2D(j,i,1) = 0   
            U_U(j) = -(y(j)-0.5)
            U_V(i) = +(x(i)-0.5)
            C_x(j) = U_U(j)
            C_y(i) = U_V(i) 
            RD(j,i) = MIN(SQRT((x(j)-0.5)**2 + (y(i)-1.25)**2),Ro)/Ro
            IC_2D(j,i,1) = (1 + COS(RD(j,i)*PI))*0.25  
        END DO
    END DO

    B_x = 0
    B_y = 0
    ! case 1 ######################################################
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'AD_LAX_SBR_1.dat'
    call tecplot(E1,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_MAC_SBR_1.dat'
    CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_ROE_SBR_1.dat'
	CALL AD_ROE_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
	call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
	DEALLOCATE(U_2D)

	! case 2 ######################################################
    H = 0.025
    DT = 0.4*H
    CALL tgridgen(T1,T2,DT,T,nt)
	CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'AD_LAX_SBR_2.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_MAC_SBR_2.dat'
    CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_ROE_SBR_2.dat'
    CALL AD_ROE_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    
	! case 3 ################################################################
    H = 0.0125
    DT = 0.4*H
    CALL tgridgen(T1,T2,DT,T,nt)
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'AD_LAX_SBR_3.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_MAC_SBR_3.dat'
    CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_ROE_SBR_3.dat'
    CALL AD_ROE_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
END PROGRAM sbr
