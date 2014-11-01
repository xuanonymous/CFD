PROGRAM swirl
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
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: U_U, U_V,B_x,C_x,B_y,C_y,g
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: RD
    REAL(DP):: mu,Ax,Bx,Ay,By,Ro, T_1
    CHARACTER(LEN=20) :: pf
    
    ! SWIRL 1#################################################
    H = 0.05
    DT = 0.4*H
    DX = 0.04
    DY = 0.04
    A = -pi
    B = pi
    T1 = 0
    T2 = 1.5
    T_1 = 1.5
    nu = 0
    
    Ax = 0
    Bx = 0.5
    Ay = 1
    By = 1.5
    Ro = 0.2
    CALL xgridgen(L,A,B,DX,X,nx)
    CALL xgridgen(L,A,B,DY,Y,ny)
    CALL tgridgen(T1,T2,DT,T,nt)
    allocate(g(nt))
    ALLOCATE(U_U(nx),U_V(ny),C_x(nx),C_y(ny),RD(nx,ny),B_x(nx),B_y(ny))
    ALLOCATE(IC_2D(nx,ny,1))
    DO n=1,nt
        g(n) = cos(pi*t(n)/T_1) 
            DO j=1,nx
                DO i=1,ny
                RD(j,i) = MIN(sqrt((x(j)-Bx)**2 + (y(i)-1.25)**2),Ro)/Ro
                IC_2D(j,i,1) = (1 + cos(pi*RD(j,i)))*0.25  
                U_U(j) = +((sin(pi*x(j))**2)*sin(2*pi*y(i))*g(n))
                U_V(i) = -((sin(pi*y(i))**2)*sin(2*pi*x(j))*g(n))
                C_x(j) = U_U(j)
                C_y(i) = U_V(i) 
                END DO
            END DO
    END DO
    B_x = 0
    B_y = 0
    ! case 1 #############################################################
    CALL AD_LAX_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    pf = 'AD_LAX_SWIRL_1.dat'
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_MAC_SWIRL_1.dat'
    CALL AD_MAC_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
    
    pf = 'AD_ROE_SWIRL_1.dat'
    CALL AD_ROE_2D(U_2D,IC_2d,nx,nt,ny,dx,dy,dt,B_x,C_x,B_y,C_y,nu)
    call tecplot(U_2D,pf,x,y,t,nx,ny,nt,1,nx,1,ny,1,nt,40)
    DEALLOCATE(U_2D)
END PROGRAM swirl
