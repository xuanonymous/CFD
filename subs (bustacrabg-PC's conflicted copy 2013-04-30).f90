MODULE SUBS
	IMPLICIT none
	CONTAINS
	!CALL surf(U_2D,pause=-1.0,terminal='windows',contour='both')
	! #############################################################

SUBROUTINE tecplot2(phi,c_file,xx,yy,nx,ny,ix1,ix2,iy1,iy2,io)
    USE types_vars
    !-----------------------------------------------
    ! Example of output for contour plotting with tecplot plotting program
    !
    ! input:
    ! phi = variable to be plotted
    ! c_file = name of file for output
    ! nx, ny = size of xx and yy arrays
    ! xx = values of x
    ! yy = values of y
    ! ix1, ix2 = beginning and ending index values for plotting in x-direction
    ! iy1, iy2 = beginning and ending index values for plotting in y-direction
    !
    ! Refer to dimension statements below for more information
    !
    !-----------------------------------------------
    IMPLICIT NONE
    INTEGER(DP), INTENT(IN) :: nx, ny
    INTEGER(DP) :: io,i,j,k,ix1,ix2,iy1,iy2,it1,it2 ! Local values
    REAL(DP), INTENT(IN) :: phi(ix1:ix2,iy1:iy2),xx(ix1:ix2),yy(iy1:iy2)
    CHARACTER(LEN=20) :: c_file
    write(*,*) nx, ny
    open(io,file=c_file,status='unknown')
    write(io,*) 'VARIABLES = "x", "y", "Value"'
    write(io,*) 'ZONE I=', nx, ', J=', ny, 'DATAPACKING=POINT'
        do j = iy1,iy2
            do i = ix1, ix2
                write(io,*) xx(i),yy(j),phi(i,j)
            end do
        end do
    close(io)
END SUBROUTINE tecplot2
SUBROUTINE tecplot(phi,c_file,xx,yy,tt,nx,ny,nt,ix1,ix2,iy1,iy2,it1,it2,io)
    USE types_vars
    !-----------------------------------------------
    ! Example of output for contour plotting with tecplot plotting program
    !
    ! input:
    ! phi = variable to be plotted
    ! c_file = name of file for output
    ! nx, ny = size of xx and yy arrays
    ! xx = values of x
    ! yy = values of y
    ! ix1, ix2 = beginning and ending index values for plotting in x-direction
    ! iy1, iy2 = beginning and ending index values for plotting in y-direction
    !
    ! Refer to dimension statements below for more information
    !
    !-----------------------------------------------
    IMPLICIT NONE
    INTEGER(DP), INTENT(IN) :: nx, ny, nt
    INTEGER(DP) :: io,i,j,k,ix1,ix2,iy1,iy2,it1,it2 ! Local values
    REAL(DP), INTENT(IN) :: phi(ix1:ix2,iy1:iy2,it1:it2),xx(ix1:ix2)
    REAL(DP), INTENT(IN) :: yy(iy1:iy2),tt(it1:it2)
    CHARACTER(LEN=20) :: c_file
    write(*,*) nx, ny, nt
    open(io,file=c_file,status='unknown')
    write(io,*) 'VARIABLES = "x", "y", "t", "Value"'
    write(io,*) 'ZONE I=', nx, ', J=', ny, ', K=', nt, 'DATAPACKING=POINT'
    do k = it1,it2
        do j = iy1,iy2
            do i = ix1, ix2
                write(io,*) xx(i),yy(j),tt(k),phi(i,j,k)
            end do
        end do
    end do
    close(io)
END SUBROUTINE tecplot
! Written by: Phillip Langsdon 2013-----Except above this line
! --------------------------------------------------------
! Grid Generation-----------------------------------------
SUBROUTINE gridgen(L,A,B,DX,X,nx)
    USE TYPES_VARS
    REAL(DP), INTENT(OUT) :: L
	REAL(DP), INTENT(OUT) , DIMENSION(:) :: X
    REAL(DP), INTENT(IN) :: DX,A,B   
	INTEGER(DP), INTENT(OUT) :: nx	
    INTEGER(DP) :: I
	L = B-A
    nx = L/DX
	ALLOCATE(X(nx))
	DO i=1,nx
		x(i)=A+i*dx
	END DO
END SUBROUTINE gridgen

! ########Initial Conditions / Boundary Conditions ##################
! ###################################################################
SUBROUTINE WAVE_IC(X,npts,UIC,N1,q)
        USE types_vars
        INTEGER(DP), INTENT(IN) :: q,n1  ! 'I' is the n value in the IC
        REAL(DP), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: X
        REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: UIC
        INTEGER(DP), INTENT(IN) :: NPTS
        INTEGER(DP) :: J
        ALLOCATE (UIC(npts,1))
        DO J=1,npts
                UIC(j,1) = sin(n1*pi*x(j)/q)
        END DO
END SUBROUTINE WAVE_IC

! Dirichlet Boundary Conditions:------------------------------------
SUBROUTINE DCBC(DX,DT,V,C,X,nx,UIC,N,DCV)
        USE types_vars
        REAL(DP), INTENT(IN) :: DX,DT,V,C,DCV 
        REAL(DP), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: X
        REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: UIC
        INTEGER(DP), INTENT(IN) :: nx,N
        INTEGER(DP) :: J
        ALLOCATE (UIC(nx,1))
        DO J=1,nx
                UIC(j,1) = DCV
        END DO
END SUBROUTINE DCBC

! Setup Riemann Discontinuous IC ##################################
SUBROUTINE RIEMANN_IC(UIC,U1,U2,x1,nx)
    USE types_vars
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: UIC
    INTEGER(DP), INTENT(IN) :: U1, U2, X1, nx
    INTEGER(DP) :: J, X3, X4
    ALLOCATE(UIC(nx,nx))
    DO j=1,x1
        UIC(j,1) = U1
    END DO
    DO j=x1+1,nx-1
        UIC(j,1) = U2
    END DO
END SUBROUTINE RIEMANN_IC
!subroutine surf(x(:), y(:), z(:,:), pause, palette, terminal, 
! filename, pm3d, contour, persist, input)
! IC : U(x,0) = sin(x) + 0.5*sin(0.5x) on 0 < x < 2pi ############
SUBROUTINE BE_SIN_IC(UIC,nx,nt,x)
    USE types_vars
	REAL(DP), INTENT(IN), DIMENSION(:) :: X
    INTEGER(DP), INTENT(IN) :: nx, nt
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: UIC
    INTEGER(DP) :: J,I
    ALLOCATE(UIC(nx,nt))
    DO J=1,nx
        		uic(j,1) = sin(x(j)) + 0.5*sin(0.5*x(j))
    END DO
END SUBROUTINE BE_SIN_IC

!######################--Derivatives--############################
! Forward Difference First Derivative-----------------------------
! Approximation---------------------------------------------------
SUBROUTINE deriv1(u,H,ux,nx)
        USE TYPES_VARS
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: u
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: ux  
        REAL(DP), INTENT(IN) :: H
        INTEGER(DP), INTENT(IN) :: nx
        INTEGER(DP) :: i
        
        ALLOCATE(UX(nx))
        DO i=1,nx
        IF (i-1>1) THEN
          ux(i)=(u(i) - u(i-1) ) / H
        END IF
        IF (i-1<1) THEN
          ux(i)=(u(i)-u(nx-1)) / H
        END IF
        END DO
END SUBROUTINE deriv1

! First Derivative O(2)------------------------------------
SUBROUTINE deriv12(U,dx,UX12,NPTS)
        USE TYPES_VARS
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: U
        REAL(DP), INTENT(IN) :: dx
        INTEGER(DP), INTENT(IN) :: NPTS
        INTEGER(DP) :: i
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: UX12
        ALLOCATE(UX12(NPTS))
        ! First order du/dx finite difference centered scheme 
        ! Second Order Truncation Accuracy
        ! Applying Periodic Boundary conditions through
        ! 1-dimensional grid
        DO i=1,NPTS,1
        IF (i-1 .LT. 1) THEN
                UX12(I)=(U(I+1)-(U(NPTS-1)))*0.5*dx
        ELSE IF (i+1 > npts) THEN
                UX12(I)=(U(2)-(U(I-1)))*0.5*dx
        ELSE   
                UX12(I)=(U(I+1)-(U(I-1)))*0.5*dx
        END IF
        END DO 
END SUBROUTINE deriv12

! First Derivative O(2) Three point formula
SUBROUTINE SOTPF(U,dx,UX123,NPTS)
        USE TYPES_VARS
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: U
        REAL(DP), INTENT(IN) :: DX
        INTEGER(DP), INTENT(IN) :: NPTS
        INTEGER(DP) :: i
        REAL(DP), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: UX123
        ALLOCATE(UX123(NPTS,1))
        ! Second Order Three Point formula
        ! For first derivative
        DO i=1,NPTS
        	IF (i+2 > NPTS .AND. i+1 .EQ. NPTS) THEN
				UX123(I,1) = (-3*U(i)+4*U(2)-U(2))*0.5*dx
			ELSE IF (i+1 < NPTS .AND. i+2 > NPTS) THEN
				UX123(I,1) = (-3*U(i)+4*U(i+1)-U(2))*0.5*dx
			ELSE IF(i+2 < NPTS .AND. i+1 < NPTS) THEN
				UX123(I,1) = (-3*U(i)+4*U(i+1)-U(i+2))*0.5*dx
			END IF
        END DO
END SUBROUTINE SOTPF

SUBROUTINE TDERIV1(f2u,H,tux,npts)
        USE TYPES_VARS
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(IN):: F2U
        REAL(DP), INTENT(IN) :: H
        INTEGER(DP), INTENT(IN) :: NPTS
        INTEGER(DP) :: I
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: TUX
        ALLOCATE(TUX(NPTS))
        DO I=1, NPTS
                TUX(I) = (f2U(I+1)-f2U(I))/H
        END DO
END SUBROUTINE TDERIV1

SUBROUTINE TDERIV12(f2u,H,TUX12,NPTS)
        USE TYPES_VARS
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: F2U
        REAL(DP), INTENT(IN) :: H
        INTEGER(DP), INTENT(IN):: NPTS
        INTEGER(DP) :: I
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: TUX12
        ALLOCATE(TUX12(NPTS))
        DO I = 1 , NPTS-1
                TUX12(I) = f2U(I+1)-2*f2U(I)+f2U(I-1)
        END DO
END SUBROUTINE TDERIV12
 
! First Derivative O(4)----------------------------------------------
SUBROUTINE deriv14(U,H,UX14,NPTS)
        USE TYPES_VARS
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: U
        REAL(DP), INTENT(IN) :: H
        INTEGER(DP), INTENT(IN) :: NPTS
        INTEGER(DP) :: i
        REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: UX14
        ALLOCATE(UX14(NPTS))
        ! First order du/dx finite difference cenetered scheme
        ! Fourth order truncation accuracy
        ! PBC Applied
        ! 1-dimensional grid
        DO i=1,NPTS,1
            IF (I-1<1) THEN
                UX14(I)=(-U(I+2)+8*U(I+1)-8*U(NPTS-1)+U(NPTS-2))/(12*H)
            ELSE IF (I-2<1) THEN
                UX14(I)=(-U(I+2)+8*U(I+1)-8*U(i-1)+U(NPTS-2))/(12*H)
            ELSE IF (I+2>NPTS .and. I+1<NPTS) THEN
                UX14(I)=(-U(2)+8*U(I+1)-8*U(I-1)+U(I-2))/(12*H)
            ELSE IF (I+2>NPTS .and. I+1>NPTS) THEN
                UX14(I)=(-U(3)+8*U(2)-8*U(I-1)+U(I-2))/(12*H)
            ELSE IF (I+1<NPTS .AND. I-2>0) THEN
                UX14(I)=(-U(I+2)+8*U(I+1)-8*U(i-1)+U(i-2))/(12*H)
            END IF
        END DO        
END SUBROUTINE deriv14

! KPP Equation Numerical Solver f(UIC)
SUBROUTINE KPP_FDE(a,dt,dx,kp,nx,nt,UIC,ku,b)
	USE types_vars
	REAL(DP), INTENT(IN) :: a,dt,B,kp,dx
	REAL(DP), INTENT(IN), ALLOCATABLE, DIMENSION(:,:)::uic
	REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:)::ku
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: del,del_2
	REAL(DP) :: a1,a2,a3,a4,a5
	INTEGER(DP), INTENT(IN) :: nx,nt
	INTEGER(DP) :: k,n

	ALLOCATE (ku(nt,nx),del(nt,nx),del_2(nt,nx))
	! initiliaze the initial conditions
	DO n=1, nx
		KU(1,n)=UIC(n,1)
	END DO
	! initilaize the boundary conditions u(1,:) = 0 , u(nx,:) = 0
	DO k=1, nt
		KU(k,1)=0
		KU(k,nx)=0
	END DO

	a1 = a/(2*dt)
	a2 = kp/(dx**2)
	a3 = b/2
	a4 = (1-a)/dt
	a5 = a1+a4+a2

	DO k=2,nt-1
        DO n=2,nx
		    IF (n-1 .lt. 0 ) THEN
				ku(k+1,n) = (a1*ku(nx-1,n)+a2*(ku(k,n+1) &
				-ku(nx-1,n)+ku(k,n-1)) + a4*ku(k,n) &
				+a3*(ku(k,n+1)*(1-ku(k,n+1)) &
				+ku(k,n-1)*(1-ku(k,n-1))))/a5
			ELSE IF (n+1 .gt. nx)THEN
				ku(k+1,n) = (a1*ku(k-1,n)+a2*(ku(k,2) &
				-ku(k-1,n)+ku(k,n-1)) + a4*ku(k,n) &
				+a3*(ku(k,2)*(1-ku(k,2)) &
				+ku(k,n-1)*(1-ku(k,n-1))))/a5
			ELSE
				ku(k+1,n) = (a1*ku(k-1,n)+a2*(ku(k,n+1) &
				-ku(k-1,n)+ku(k,n-1)) + a4*ku(k,n) &
				+a3*(ku(k,n+1)*(1-ku(k,n+1)) &
				+ku(k,n-1)*(1-ku(k,n-1))))/a5
			END IF
		END DO
	END DO
END SUBROUTINE KPP_FDE

! KPP Equation - Traveling wave solution-------------------------]
SUBROUTINE KPP_TW(r,alpha,dt,nx,nt,ktw,B1,C,x,t)
        USE types_vars
        REAL(DP), INTENT(IN) :: ALPHA, B1,r,dt
        REAL(DP), INTENT(IN), DIMENSION(:) :: T, X
        REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: ktw
        INTEGER(DP), INTENT(IN) :: nx,nt    
        REAL(DP) :: C   
        INTEGER(DP) :: K,N
        ALLOCATE (ktw(nx,nt))
        DO n=1, nt
                ktw(1,n)=1
                ktw(nx,n)=0
        END DO
        C=4
        DO N=1,nt
                DO K=1,nx-1
                  ktw(k,n) = 1/((1+C*exp((-5/6)*b1*t(n) &
							 +SQRT(B1))*x(k)/6)**2)
                END DO
        END DO
END SUBROUTINE KPP_TW
! Lax Wave Method------------------------------------------------
SUBROUTINE LAX_WAVE(UIC,U,DX,DT,C,NX,NT)
        USE types_vars
        REAL(DP), INTENT(IN) :: DX,DT,C
        REAL(DP), INTENT(IN), ALLOCATABLE, DIMENSION(:,:) :: UIC
        REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: U
        INTEGER(DP) :: j,n
        INTEGER(DP), INTENT(IN) :: nx,nt
        REAL(DP) :: V
        ALLOCATE(U(nx,nt))
        V = C*DT/DX
        do j=1,nx
                U(j,1) = UIC(j,1)
        end do
        do j=1,nx
			do n=2,nt-1
             if (j-1<1) then
                u(j,n+1) = u(j,n) + 0.5*(u(j+1,n)-2*u(j,n)+u(nx,n))     &
                -V*0.5*(u(j+1,n)-u(nx,n))
             else if (j+1>nx) then
                u(j,n+1) = u(j,n) + 0.5*(u(2,n)-2*u(j,n)+u(j-1,n))      &
                -V*0.5*(u(2,n)-u(j-1,n))
             else 
                 u(j,n+1) = u(j,n) + 0.5*(u(j+1,n)-2*u(j,n)+u(j-1,n))   &
                 -V*0.5*(u(j+1,n)-u(j-1,n))
             end if
           end do
        end do
END SUBROUTINE LAX_WAVE

! Lax-Wendroff Method---###########################################
! Wave Equation Solve
SUBROUTINE LW_WAVE(DX,DT,V,C,X,T,UIC,LWU,nx,nt)
        USE types_vars
        REAL(DP), INTENT(IN) :: DX,DT,C,V
        REAL(DP), INTENT(IN), ALLOCATABLE, DIMENSION(:,:) :: UIC
        REAL(DP), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: X, T
        REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: LWU
        INTEGER(DP) :: j,n
        INTEGER(DP), INTENT(IN) :: nx,nt
        ALLOCATE(LWU(nx,nt)) 
        do j=1,nx
                LWU(j,1) = UIC(j,1)
        end do
        do j=1,nx
              do n=1,nt-1
                if (j-1<1) then
                    lwu(j,n+1) = &
				    lwu(j,n)-c*dt*((lwu(j+1,n)-lwu(nx,n))/(2*dx) &
                    -0.5*c*dt*dt*(lwu(j+1,n)-2*lwu(j,n)+lwu(nx,n))/(dx*dx))
                else if (j+1>nx) then
                    lwu(j,n+1) = & 
					lwu(j,n)-c*dt*((lwu(2,n)-lwu(j-1,n))/(2*dx) &
                    -0.5*c*c*dt*dt*(lwu(2,n)-2*lwu(j,n)+lwu(j-1,n))/(dx*dx))
                else if (j>0 .and. j+1<nx) then
                    lwu(j,n+1) = & 
					lwu(j,n)-c*dt*((lwu(j+1,n)-lwu(j-1,n))/2*dx) &
                    -0.5*c*c*dt*dt*(lwu(j+1,n)-2*lwu(j,n)+lwu(j-1,n))/dx*dx
                end if
            end do
        end do
END SUBROUTINE LW_WAVE

! 2D Advection Diffusion Lax Method ######################################
! Ut + u * Ux + v * Uy = *nu(Uxx + Uyy)
! Uxx and Uyy Central Differenced O(2) 
! Lax Method here is O(1) so O(2) for viscous terms is acceptable
SUBROUTINE AD_LAX_2D(U,UIC,nx,nt,ny,DX,DY,DT,b,c,d,e,nu)
    USE types_vars
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:,:,:) :: UIC
    REAL(DP), INTENT(OUT), DIMENSION(:,:,:), ALLOCATABLE :: U
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: F,G
    REAL(DP), INTENT(IN), DIMENSION(:) :: b, c, d, e
    REAL(DP), INTENT(IN) :: DT, DX, DY, nu
    REAL(DP) :: vy, vx, rx, ry, r, v
    INTEGER(DP), INTENT(IN) :: nx, ny, nt
    INTEGER(DP) :: j,i,n
    ! B=0, C=1 -> 2D Linear Advection-Diffusion Equation
    ! F = cu + b*u^2/2
    ! G = cu + b*u^2/2 
    !  : used to switch viscous terms on and off
    !  = 1 -> viscous terms remain
    !  = 0 -> becomes inviscid equation
    ALLOCATE(U(nx,ny,nt),F(nx,ny,nt),G(nx,ny,nt))
    ! Initialize the IC's ############################################
    DO j=1,nx
        DO i=1,ny
            U(j,i,1) = UIC(j,i,1)
            ! Initialize F = cu + b*u^2/2
            F(j,i,1) = (U(j,i,1)**2)*b(j)*0.5 + U(j,i,1)*c(j)
            G(j,i,1) = (U(j,i,1)**2)*d(i)*0.5 + U(j,i,1)*e(i)
            END DO
    END DO
    vy =  0.5*DT/DY
    vx =  0.5*DT/DX
    v = (maxval(uic))*dt/dx
    rx = nu*DT/(DX*DX)
    ry = nu*DT/(DY*DY) 
    r = rx + ry
    write(*,*) '2D-LAX Stability r =:', rx, 'cfl:', v
    DO J=1,nx
      DO I=1,ny
        DO N=1,nt-1
            IF( i + 1 .GT. ny .and. J + 1 .GT. nx) THEN
            U(j,i,n+1) = (u(2,2,n)+u(j-1,i-1,n))*0.5                        &
                        -(F(2,i,n)-F(j-1,i,n))*vx                           &
                        -(G(j,2,n)-G(j,i-1,n))*vy                           &
                          +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx       &
                          +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry         
            end if
            
            if ( i + 1 .gt. ny .and. j - 1 .LT. 1 ) THEN
                U(j,i,n+1) = (u(j+1,2,n)+u(nx-1,i-1,n))*0.5                 &
                           -(F(j+1,i,n)-F(nx-1,i,n))*vx                     &
                           -(G(j,2,n)-G(j,i-1,n))*vy                        &
                          +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx    &
                          +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry       
            end if
            IF ( i + 1 .gt. ny .and. j + 1 .lt. ny .and. j - 1 .GT. 1 ) THEN 
            U(j,i,n+1) = (u(j+1,2,n)+u(j-1,i-1,n))*0.5                      &
                         -(F(j+1,i,n)-F(j-1,i,n))*vx                        &
                         -(G(j,2,n)-G(j,i-1,n))*vy                          &
                          +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx     &
                          +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry       
            end if
            IF ( i - 1 .LT. 1 .and. J + 1 .GT. nx ) THEN
            U(j,i,n+1) = (u(2,i+1,n)+u(j-1,ny-1,n))*0.5                     &
                         -(F(2,i,n)-F(j-1,i,n))*vx                          &
                         -(G(j,i+1,n)-G(j,ny-1,n))*vy                       &
                          +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx       &
                          +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry
			end if
            if ( i - 1 .lt. 1 .and. j - 1 .LT. 1 ) THEN  
            U(j,i,n+1) = (u(j+1,i+1,n)+u(nx-1,ny-1,n))*0.5                  &
                         -(F(j+1,i,n)-F(nx-1,i,n))*vx                       &
                         -(G(j,i+1,n)-G(j,ny-1,n))*vy                       &
                          +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx    &
                          +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry          
            END IF
            IF ( i - 1 .lt. 1 .and. j + 1 .lt. ny .and. j - 1 .GT. 1 ) THEN
            U(j,i,n+1) = (u(j+1,i+1,n)+u(j-1,ny-1,n))*0.5                   &
                         -(F(j+1,i,n)-F(j-1,i,n))*vx                        &
                         -(G(j,i+1,n)-G(j,ny-1,n))*vy                       &
                          +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx     &
                          +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry 
			END IF
            IF ( i + 1 .lt. ny .and. i - 1 .GT. 1 .AND. J + 1 .GT. nx ) THEN           
            U(j,i,n+1) = (u(2,i+1,n)+u(j-1,i-1,n))*0.5                      &
                         -(F(2,i,n)-F(j-1,i,n))*vx                          &
                         -(G(j,i+1,n)-G(j,i-1,n))*vy                        &
                          +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx       &
                          +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry 
			END IF
            if ( i + 1 .lt. ny .and. i - 1 .GT. 1 .AND. j - 1 .LT. 1 ) THEN 
            U(j,i,n+1) = (u(j+1,i+1,n)+u(nx-1,i-1,n))*0.5                   &
                         -(F(j+1,i,n)-F(nx-1,i,n))*vx                       &
                         -(G(j,i+1,n)-G(j,i-1,n))*vy                        &
                          +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx    &
                          +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry       
            END IF
            IF (i+1.lt.ny.and.i-1.GT.1.AND.j+1.lt.nx.and. j-1 .GT. 1) THEN 
            U(j,i,n+1) = (u(j+1,i+1,n)+u(j-1,i-1,n))*0.5                    &
                         -(F(j+1,i,n)-F(j-1,i,n))*vx                        &
                         -(G(j,i+1,n)-G(j,i-1,n))*vy                        &
                          +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx     &
                          +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry 
            END IF
            F(j,i,n+1) = (U(j,i,n+1)**2)*b(j)*0.5 + U(j,i,n+1)*c(j)
            G(j,i,n+1) = (U(j,i,n+1)**2)*d(i)*0.5 + U(j,i,n+1)*e(i)   
        END DO
      END DO
   END DO
END SUBROUTINE AD_LAX_2D

! 2D ADE MacCormack Method ##########################################
! MacCormack (1969) Predictor-Corrector version of the LW Scheme
! Jacobian does not appear
! O(2) central differenced viscous 
SUBROUTINE AD_MAC_2D(U,UIC,nx,nt,ny,DX,DY,DT,b,c,d,e,nu)
    USE types_vars
    REAL(DP), INTENT(IN), DIMENSION(:,:,:) :: UIC
    REAL(DP), INTENT(OUT), DIMENSION(:,:,:), ALLOCATABLE :: U
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: F, Fp, G, Gp, Up
    REAL(DP), INTENT(IN), DIMENSION(:) :: B,C,D,E
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: vx,vy
    ! Up : U predicted
    INTEGER(DP), INTENT(IN) :: nx, nt, ny
    REAL(DP), INTENT(IN) :: DT,DX,DY,nu
    REAL(DP) :: rx, ry, r
    INTEGER(DP) :: n,j,i
    ALLOCATE(U(nx,ny,nt),F(nx,ny,nt),G(nx,ny,nt),vx(nx),Vy(ny))
	ALLOCATE(Fp(nx,ny,nt),Gp(nx,ny,nt),UP(nx,ny,nt))
    ! Initialize the IC's ############################################
    ! Initialize F = b*u2/2 + cu
    DO j=1,nx
		DO i=1,ny
        		U(j,i,1) = UIC(j,i,1)
        		F(j,i,1) = (U(j,i,1)**2)*b(j)*0.5 + U(j,i,1)*c(j)
        		G(j,i,1) = (U(j,i,1)**2)*d(i)*0.5 + U(j,i,1)*e(i)
		END DO
    END DO
    DO J=1,NX
        DO I=1,NY
            vX(J) = DT/DX
            vY(I) = DT/DY   
        END DO
    END DO
	rx = nu*DT/(DX*DX)
	ry = nu*DT/(DY*DY) 
	r = (rx + ry)
    write(*,*) '2D-MAC r =:', rx, '<= 0.5', 'cfl:', maxval(vx)*maxval(uic), '<=1'
    DO N=1,nt-1
      DO J=1,nx
          DO I=1,ny
            IF ( i + 1 .GT. ny .and. J + 1 .GT. nx ) THEN 
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(2,2,n)-F(j,i,n))                        &
                            - vy(i)*(G(2,2,n)-G(j,i,n))                        &
                            +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx*0.5      &
                            +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(j-1,i-1,n+1))                   &
                            -vy(i)*(Gp(j,i,n+1)-Gp(j-1,i-1,n+1)))*0.5              &
                            +(Up(2,i,n+1)-2*up(j,i,n+1) +Up(j-1,i,n+1))*rx*0.5   &
                            +(Up(j,2,n+1)-2*up(j,i,n+1) +Up(j,i-1,n+1))*ry*0.5
            end if
            if ( i + 1 .gt. ny .and. j - 1 .LT. 1 ) THEN
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(j+1,2,n)-F(j,i,n))                        &
                            - vy(i)*(G(j+1,2,n)-G(j,i,n))                        &
                            +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx*0.5      &
                            +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(nx-1,i-1,n+1))                   &
                            -vy(i)*(Gp(j,i,n+1)-Gp(nx-1,i-1,n+1)))*0.5              &
                            +(Up(j+1,i,n+1)-2*up(j,i,n+1) +Up(nx-1,i,n+1))*rx*0.5   &
                            +(Up(j,2,n+1)-2*up(j,i,n+1) +Up(j,i-1,n+1))*ry*0.5
            end if
            IF ( i + 1 .gt. ny .and. j + 1 .lt. nx .and. j - 1 .GT. 1 ) THEN
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(j+1,2,n)-F(j,i,n))                        &
                            - vy(i)*(G(j+1,2,n)-G(j,i,n))                        &
                            +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx*0.5      &
                            +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(j-1,i-1,n+1))                   &
                            -vy(i)*(Gp(j,i,n+1)-Gp(j-1,i-1,n+1)))*0.5              &
                            +(Up(j+1,i,n+1)-2*up(j,i,n+1) +Up(j-1,i,n+1))*rx*0.5   &
                            +(Up(j,2,n+1)-2*up(j,i,n+1) +Up(j,i-1,n+1))*ry*0.5
            end if
            IF ( i - 1 .LT. 1 .and. J + 1 .GT. nx ) THEN
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(2,i+1,n)-F(j,i,n))                        &
                            - vy(i)*(G(2,i+1,n)-G(j,i,n))                        &
                            +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx*0.5      &
                            +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(j-1,ny-1,n+1))                   &
                            -vy(i)*(Gp(j,i,n+1)-Gp(j-1,ny-1,n+1)))*0.5              &
                            +(Up(2,i,n+1)-2*up(j,i,n+1) +Up(j-1,i,n+1))*rx*0.5   &
                            +(Up(j,i+1,n+1)-2*up(j,i,n+1) +Up(j,ny-1,n+1))*ry*0.5
            end if
            if ( i - 1 .LT. 1  .and. j - 1 .LT. 1 ) THEN
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(j+1,i+1,n)-F(j,i,n))                        &
                            - vy(i)*(G(j+1,i+1,n)-G(j,i,n))                        &
                            +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx*0.5      &
                            +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(nx-1,ny-1,n+1))                   &
                            -vy(i)*(Gp(j,i,n+1)-Gp(nx-1,ny-1,n+1)))*0.5              &
                            +(Up(j+1,i,n+1)-2*up(j,i,n+1) +Up(nx-1,i,n+1))*rx*0.5   &
                            +(Up(j,i+1,n+1)-2*up(j,i,n+1) +Up(j,ny-1,n+1))*ry*0.5
            end if
            IF ( i - 1 .LT. 1  .and. j + 1 .lt. nx .and. j - 1 .GT. 1 ) THEN
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(j+1,i+1,n)-F(j,i,n))                        &
                            - vy(i)*(G(j+1,i+1,n)-G(j,i,n))                        &
                            +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx*0.5      &
                            +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(j-1,ny-1,n+1))                   &
                            -vy(i)*(Gp(j,i,n+1)-Gp(j-1,ny-1,n+1)))*0.5              &
                            +(Up(j+1,i,n+1)-2*up(j,i,n+1) +Up(j-1,i,n+1))*rx*0.5   &
                            +(Up(j,i+1,n+1)-2*up(j,i,n+1) +Up(j,ny-1,n+1))*ry*0.5
            end if
            IF ( i + 1 .lt. ny .and. i - 1 .GT. 1 .AND. J + 1 .GT. nx ) THEN
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(2,i+1,n)-F(j,i,n))                        &
                            - vy(i)*(G(2,i+1,n)-G(j,i,n))                        &
                            +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx*0.5      &
                            +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(j-1,i-1,n+1))                   &
                            -vy(i)*(Gp(j,i,n+1)-Gp(j-1,i-1,n+1)))*0.5              &
                            +(Up(2,i,n+1)-2*up(j,i,n+1) +Up(j-1,i,n+1))*rx*0.5     &
                            +(Up(j,i+1,n+1)-2*up(j,i,n+1) +Up(j,i-1,n+1))*ry*0.5
            end if
            if ( i + 1 .lt. ny .and. i - 1 .GT. 1 .AND. j - 1 .LT. 1 ) THEN
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(j+1,i+1,n)-F(j,i,n))                        &
                            - vy(i)*(G(j+1,i+1,n)-G(j,i,n))                        &
                            +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx*0.5     &
                            +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(nx-1,i-1,n+1))                  &
                            -vy(i)*(Gp(j,i,n+1)-Gp(nx-1,i-1,n+1)))*0.5             &
                            +(Up(j+1,i,n+1)-2*up(j,i,n+1) +Up(nx-1,i,n+1))*rx*0.5  &
                            +(Up(j,i+1,n+1)-2*up(j,i,n+1) +Up(j,i-1,n+1))*ry*0.5
            end if
            IF (i+1.lt.ny.and.i-1.GT.1.AND.j+1.lt.nx.and. j-1 .GT. 1) THEN
                ! Predictor Step
                Up(j,i,n+1) = u(j,i,n)                                             &
                            - vx(j)*(F(j+1,i+1,n)-F(j,i,n))                        &
                            - vy(i)*(G(j+1,i+1,n)-G(j,i,n))                        &
                            +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx*0.5      &
                            +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry*0.5
                ! Update F and G predicted from Up            
                Fp(j,i,n+1) = (Up(j,i,n+1)**2)*b(j)*0.5 + Up(j,i,n+1)*c(j)
                Gp(j,i,n+1) = (Up(j,i,n+1)**2)*d(i)*0.5 + Up(j,i,n+1)*e(i)
                ! Correction Step
                U(j,i,n+1) = (u(j,i,n)+Up(j,i,n+1)                                 &
                            -vx(j)*(Fp(j,i,n+1)-Fp(j-1,i-1,n+1))                   &
                            -vy(i)*(Gp(j,i,n+1)-Gp(j-1,i-1,n+1)))*0.5              &
                            +(Up(j+1,i,n+1)-2*up(j,i,n+1) +Up(j-1,i,n+1))*rx*0.5   &
                            +(Up(j,i+1,n+1)-2*up(j,i,n+1) +Up(j,i-1,n+1))*ry*0.5
            end if
                F(j,i,n+1) = (U(j,i,n+1)**2)*0.5*b(j) + c(j)*u(j,i,n+1)
                G(j,i,n+1) = (U(j,i,n+1)**2)*0.5*d(i) + e(i)*u(j,i,n+1)
        END DO
      END DO
   END DO
END SUBROUTINE AD_MAC_2D

SUBROUTINE AD_ROE_2D(U,UIC,nx,nt,ny,DX,DY,DT,b,c,d,e,nu)
    USE types_vars
    REAL(DP), INTENT(IN), DIMENSION(:,:,:) :: UIC
    REAL(DP), INTENT(OUT), DIMENSION(:,:,:), ALLOCATABLE :: U
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: F, Fp, G, Gp, Up
    ! Up : U predicted
    INTEGER(DP), INTENT(IN) :: nx, nt, ny
    REAL(DP), INTENT(IN) :: DT, DX, DY, nu
    REAL(DP) :: v, rx, ry, r
    REAL(DP), INTENT(IN), DIMENSION(:) :: B,C,D,E
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: vx, vy
    INTEGER(DP) :: n,j,i
    ALLOCATE(U(nx,ny,nt),F(nx,ny,nt),G(nx,ny,nt))
    ! Initialize the IC's ############################################
    ! Initialize F = b*u2/2 + cu
    ALLOCATE(VX(NX),VY(NY))
    DO j=1,nx
		DO i=1,ny
        		U(j,i,1) = UIC(j,i,1)
        		F(j,i,1) = (U(j,i,1)**2)*b(J)*0.5 + U(j,i,1)*c(J)
        		G(j,i,1) = (U(j,i,1)**2)*d(i)*0.5 + U(j,i,1)*e(i)
		END DO
    END DO
    DO J=1,NX
        DO I=1,NY
            vX(J) =  DT/DX
            vY(I) =  DT/DY   
        END DO
    END DO
    rx = nu*DT/(DX*DX)
    ry = nu*DT/(DY*DY) 
    r = rx + ry
    write(*,*) '2D-ROE Stability r =:', rx, 'cfl :', vx(1)*maxval(uic)
    DO J=1,nx
      DO I=1,ny
        DO N=1,nt-1
            IF ( i + 1 .GT. ny .and. J + 1 .GT. nx ) THEN 
                    u(j,i,n+1) = u(j,i,n) &
                    -((F(2,2,n)-F(j-1,i-1,n))  &
                    -abs(c(j)+(u(2,2,n)+u(j,i,n))*b(j)*0.5)*(u(2,2,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(j-1,i-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(j-1,i-1,n)))*(DT/DX)*0.5 &
                    -((G(2,2,n)-G(j-1,i-1,n))  &
                    -abs(e(i)+(u(2,2,n)+u(j,i,n))*d(i)*0.5)*(u(2,2,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(j-1,i-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(j-1,i-1,n)))*(DT/DY)*0.5 &
                    +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx*0.5     &
                    +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry*0.5  
            end if
            if ( i + 1 .gt. ny .and. j - 1 .LT. 1 ) THEN
                    u(j,i,n+1) = u(j,i,n) & 
                    -((F(j+1,2,n)-F(nx-1,i-1,n))  &
                    -abs(c(j)+(u(j+1,2,n)+u(j,i,n))*b(j)*0.5)*(u(j+1,2,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(nx-1,i-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(nx-1,i-1,n)))*(DT/DX)*0.5 &
                    -((G(j+1,2,n)-G(nx-1,i-1,n))  &
                    -abs(e(i)+(u(j+1,2,n)+u(j,i,n))*d(i)*0.5)*(u(j+1,2,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(nx-1,i-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(nx-1,i-1,n)))*(DT/DY)*0.5 &
                    +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx    &
                    +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry
            end if
            IF ( i + 1 .gt. ny .and. j + 1 .lt. nx .and. j - 1 .GT. 1 ) THEN
                    u(j,i,n+1) = u(j,i,n) &
                    -((F(j+1,2,n)-F(j-1,i-1,n))  &
                    -abs(c(j)+(u(j+1,2,n)+u(j,i,n))*b(j)*0.5)*(u(j+1,2,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(j-1,i-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(j-1,i-1,n)))*(DT/DX)*0.5 &
                    -((G(j+1,2,n)-G(j-1,i-1,n))  &
                    -abs(e(i)+(u(j+1,2,n)+u(j,i,n))*d(i)*0.5)*(u(j+1,2,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(j-1,i-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(j-1,i-1,n)))*(DT/DY)*0.5 &
                    +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx     &
                    +(U(j,2,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry  
            end if
            IF ( i - 1 .LT. 1 .and. J + 1 .GT. nx ) THEN
                    u(j,i,n+1) = u(j,i,n) & 
                    -((F(2,i+1,n)-F(j-1,ny-1,n))  &
                    -abs(c(j)+(u(2,i+1,n)+u(j,i,n))*b(j)*0.5) &
                    *(u(2,i+1,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(j-1,ny-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(j-1,ny-1,n)))*(DT/DX)*0.5 &
                    -((G(2,i+1,n)-G(j-1,ny-1,n))  &
                    -abs(e(i)+(u(2,i+1,n)+u(j,i,n))*d(i)*0.5)*(u(2,i+1,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(j-1,ny-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(j-1,ny-1,n)))*(DT/DY)*0.5 &
                    +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx     &
                    +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry 
            end if
            if ( i - 1 .LT. 1  .and. j - 1 .LT. 1 ) THEN
                    u(j,i,n+1) = u(j,i,n) &
                    -((F(j+1,i+1,n)-F(nx-1,ny-1,n))  &
                    -abs(c(j)+(u(j+1,i+1,n)+u(j,i,n))*b(j)*0.5)*(u(j+1,i+1,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(nx-1,ny-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(nx-1,ny-1,n)))*(DT/DX)*0.5 &
                    -((G(j+1,i+1,n)-G(nx-1,ny-1,n))  &
                    -abs((e(i)+(u(j+1,i+1,n)+u(j,i,n))*d(i)*0.5)*(u(j+1,i+1,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(nx-1,ny-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(nx-1,ny-1,n))))*(DT/DY)*0.5 &
                    +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx     &
                    +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry 
            end if
            IF ( i - 1 .LT. 1  .and. j + 1 .lt. nx .and. j - 1 .GT. 1 ) THEN
                    u(j,i,n+1) = u(j,i,n) &
                    -((F(j+1,i+1,n)-F(j-1,ny-1,n))  &
                    -abs(c(j)+(u(j+1,i+1,n)+u(j,i,n))*b(j)*0.5)*(u(j+1,i+1,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(j-1,ny-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(j-1,ny-1,n)))*(DT/DX)*0.5 &
                    -((G(j+1,i+1,n)-G(j-1,ny-1,n))  &
                    -abs(e(i)+(u(j+1,i+1,n)+u(j,i,n))*d(i)*0.5)*(u(j+1,i+1,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(j-1,ny-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(j-1,ny-1,n)))*(DT/DY)*0.5 &
                    +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx     &
                    +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,ny-1,n+1))*ry       
            end if
            IF ( i + 1 .lt. ny .and. i - 1 .GT. 1 .AND. J + 1 .GT. nx ) THEN
                    u(j,i,n+1) = u(j,i,n) &
                    -((F(2,i+1,n)-F(j-1,i-1,n))  &
                    -abs(c(j)+(u(2,i+1,n)+u(j,i,n))*b(j)*0.5)*(u(2,i+1,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(j-1,i-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(j-1,i-1,n)))*(DT/DX)*0.5 &
                    -((G(2,i+1,n)-G(j-1,i-1,n))  &
                    -abs(e(i)+(u(2,i+1,n)+u(j,i,n))*d(i)*0.5)*(u(2,i+1,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(j-1,i-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(j-1,i-1,n)))*(DT/DY)*0.5 &
                    +(U(2,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx     &
                    +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry 
            end if
            if ( i + 1 .lt. ny .and. i - 1 .GT. 1 .AND. j - 1 .LT. 1 ) THEN
                    u(j,i,n+1) = u(j,i,n) &       
                    -((F(j+1,i+1,n)-F(nx-1,i-1,n))  &
                    -abs(c(j)+(u(j+1,i+1,n)+u(j,i,n))*b(j)*0.5)*(u(j+1,i+1,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(nx-1,i-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(nx-1,i-1,n)))*(DT/DX)*0.5 &
                    -((G(j+1,i+1,n)-G(nx-1,i-1,n))  &
                    -abs(e(i)+(u(j+1,i+1,n)+u(j,i,n))*d(i)*0.5)*(u(j+1,i+1,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(nx-1,i-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(nx-1,i-1,n)))*(DT/DY)*0.5 &
                    +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(nx-1,i,n+1))*rx     &
                    +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry 
            end if
            IF (i+1.lt.ny.and.i-1.GT.1.AND.j+1.lt.nx.and.j-1.GT.1) THEN
                    u(j,i,n+1) = u(j,i,n) &
                    -((F(j+1,i+1,n)-F(j-1,i-1,n))  &
                    -abs(c(j)+(u(j+1,i+1,n)+u(j,i,n))*b(j)*0.5)*(u(j+1,i+1,n)-u(j,i,n)) &
                    +abs(c(j) + (u(j,i,n) + u(j-1,i-1,n))*b(j)*0.5) &
                    *(u(j,i,n)-u(j-1,i-1,n)))*(DT/DX)*0.5 &
                    -((G(j+1,i+1,n)-G(j-1,i-1,n))  &
                    -abs(e(i)+(u(j+1,i+1,n)+u(j,i,n))*d(i)*0.5)*(u(j+1,i+1,n)-u(j,i,n)) &
                    +abs(e(i) + (u(j,i,n) + u(j-1,i-1,n))*d(i)*0.5) &
                    *(u(j,i,n)-u(j-1,i-1,n)))*(DT/DY)*0.5 &
                    +(U(j+1,i,n+1)-2*u(j,i,n+1) +U(j-1,i,n+1))*rx     &
                    +(U(j,i+1,n+1)-2*u(j,i,n+1) +U(j,i-1,n+1))*ry 
            end if
                F(j,i,n+1) = (U(j,i,n+1)**2)*0.5*b(j) + c(j)*u(j,i,n+1)
                G(j,i,n+1) = (U(j,i,n+1)**2)*0.5*d(i) + e(i)*u(j,i,n+1)
        END DO
      END DO
   END DO
END SUBROUTINE AD_ROE_2D

SUBROUTINE VTE_ROE_2D(U,nx,ny,DX,DY,DT,b,c,d,e,nu,q)
    USE types_vars
    REAL(DP), INTENT(INOUT), DIMENSION(:,:,:), ALLOCATABLE :: U
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: F, Fp, G, Gp, Up
    ! Up : U predicted
    INTEGER(DP), INTENT(IN) :: nx, ny,q
    REAL(DP), INTENT(IN) :: DX, DY, nu,DT
    REAL(DP), INTENT(IN), DIMENSION(:,:,:) :: B,D
    REAL(DP) :: v, rx, ry, r
    REAL(DP), INTENT(IN), DIMENSION(:) :: C,E
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: vx, vy
    INTEGER(DP) :: n,j,i
    ALLOCATE(F(nx,ny,q),G(nx,ny,q))
    ! Initialize the IC's ############################################
    ! Initialize F = b*u2/2 + cu
    ALLOCATE(VX(NX),VY(NY))
    DO i=1,nx
		DO j=1,ny
        		F(i,j,q) = (U(i,j,q)**2)*b(i,j,q)*0.5 + U(i,j,q)*c(i)
        		G(i,j,q) = (U(i,j,q)**2)*d(i,j,q)*0.5 + U(i,j,q)*e(j)
		END DO
    END DO
    vX(:) =  DT/DX
    vY(:) =  DT/DY   
    rx = nu*DT/(DX*DX)
    ry = nu*DT/(DY*DY) 
    r = rx + ry
    write(*,*) '2D-ROE Stability r =:', rx, 'cfl :', vx(1)*maxval(u)
    DO J=2,nx-1
      DO I=2,ny-1
              u(i,j,q+1) = u(i,j,q) &
              -((F(i+1,j,q)-F(i-1,j,q))  &
               -abs(c(i)+(u(i+1,j,q)+u(i,j,q))*b(i,j,q)*0.5)*(u(i+1,j,q)-u(i,j,q)) &
               +abs(c(i) + (u(i,j,q) + u(i-1,j,q))*b(i,j,q)*0.5) &
               *(u(i,j,q)-u(i-1,j,q)))*(DT/DX)*0.5 &
               -((G(i,j+1,q)-G(i,j-1,q))  &
               -abs(e(j)+(u(i,j+1,q)+u(i,j,q))*d(i,j,q)*0.5)*(u(i,j+1,q)-u(i,j,q)) &
               +abs(e(j) + (u(i,j,q) + u(i,j-1,q))*d(i,j,q)*0.5) &
               *(u(i,j,q)-u(i-1,j-1,q)))*(DT/DY)*0.5 &
               +(U(i+1,j,q)-2*u(i,j,q) +U(i-1,j,q))*rx     &
               +(U(i,j+1,q)-2*u(i,j,q) +U(i,j-1,q))*ry 
        END DO
      END DO
END SUBROUTINE VTE_ROE_2D

! Burgers Lax Method ###############################################
SUBROUTINE BURGERS_LAX(U,UIC,nx,nt,DX,DT,b,c)
    USE types_vars
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: UIC
    REAL(DP), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: U
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: F
    INTEGER(DP), INTENT(IN) :: nx, nt
    REAL(DP), INTENT(IN) :: DT, DX, b, c
    INTEGER(DP) :: n,j
    
    ALLOCATE(U(nx,nt),F(nx,nt))
    ! Initialize the IC's ############################################
    U(1:nx,1) = UIC(1:nx,1)
    ! Initialize F = u^2
    F(1:nx,1) = (U(1:nx,1)**2)*b*0.5 + U(1:nx,1)*c
    
    DO J=1,nx
        DO N=1,nt-1
            IF (j+1 .GT. nx) THEN
                U(j,n+1) = (U(2,n)+U(j-1,n)-(DT/DX)*(F(2,n)-F(j-1,n)))*0.5
            END IF
            IF (j-1 .LT. 1 ) THEN
                U(j,n+1) = (U(j+1,n)+U(nx-1,n)-(DT/DX)*(F(j+1,n)-F(nx-1,n)))*0.5
            END IF
            IF (j+1 .LT. nx .AND. j-1 .GT. 1 ) THEN
                U(j,n+1) = (U(j+1,n)+U(j-1,n)-(DT/DX)*(F(j+1,n)-F(j-1,n)))*0.5
            END IF
                F(j,n+1) = (U(j,n+1)**2)*b*0.5 + U(j,n+1)*c
        END DO
   END DO
END SUBROUTINE BURGERS_LAX

! Burgers MacCormack Method ##########################################
! MacCormack (1969) Predictor-Corrector version of the LW Scheme
! Jacobian does not appear
! Invisicid Burger's Equation Solver-Good resolution for discontinuities 
SUBROUTINE BURGERS_MAC(U,UIC,nx,nt,DX,DT,b,c)
    USE types_vars
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: UIC
    REAL(DP), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: U
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: F, Fp, Up
    INTEGER(DP), INTENT(IN) :: nx, nt
    REAL(DP), INTENT(IN) :: DT, DX, b, c
    REAL(DP) :: v, vv
    INTEGER(DP) :: n,j
    ALLOCATE(U(nx,nt),F(nx,nt),Fp(nx,nt),UP(nx,nt))
    ! Initialize the IC's ############################################
    ! Initialize F = b*u2/2 + cu
    DO j=1,nx
        U(j,1) = UIC(j,1)
        F(j,1) = (U(j,1)**2)*b*0.5 + U(j,1)*c
    END DO
    v = (DT/DX)
	vv = v*abs(MAXVAL(UIC))
	WRITE(*,*) 'MacCormack Method (Burgers) Stability ', v, '.LEQ. 1'
    DO J=1,nx
        DO N=1,nt-1
            IF (j+1 .GT. nx) THEN
                ! Predictor Step
                Up(j,n+1) = u(j,n) - (F(2,n)-F(j,n))*v
                ! Update Fp
                Fp(j,n+1) = (Up(j,n+1)**2)*b*0.5 + Up(j,n+1)*c
                ! Correction Step
                U(j,n+1) = (u(j,n)+Up(j,n+1)-(Fp(j,n+1)-Fp(j-1,n+1)*v))*0.5
            END IF
            IF (j-1 .LT. 1) THEN
                ! Predictor Step
                Up(j,n+1) = u(j,n) - v*(F(j+1,n)-F(j,n))*0.5
                ! Update Fp
                Fp(j,n+1) = (Up(j,n+1)**2)*b*0.5 + Up(j,n+1)*c
                ! Correction Step
                U(j,n+1) = (u(j,n)+Up(j,n+1)-v*(Fp(j,n+1)-Fp(nx-1,n+1)))*0.5
            END IF
            IF (j+1 .LT. nx .AND. j-1 .GT. 1 ) THEN
                ! Predictor Step
                Up(j,n+1) = u(j,n) - (F(j+1,n)-F(j,n))*0.5*v
                ! Update Fp
                Fp(j,n+1) = (Up(j,n+1)**2)*b*0.5 + Up(j,n+1)*c
                ! Correction Step
                U(j,n+1) = (u(j,n)+Up(j,n+1)-v*(Fp(j,n+1)-Fp(j-1,n+1)))*0.5
            END IF
            F(j,N+1) = (U(j,n+1)**2)*b*0.5 + U(j,n+1)*c
        END DO
    END DO
END SUBROUTINE BURGERS_MAC

!END SUBROUTINE NEWTONSMETHOD(F,FP,Pr,nx)
!REAL(DP), INTENT(INOUT), DIMENSION(:) PR
!REAL(DP) :: pLR,aL,aR,aRL,p_R

!PLR = Pr(1)/Pr(NX) 
!aL = (gamma(1)+1)/(gamma(1)-1)
!aR = (gamma(nx)+1)/(gamma(nx)-1)
!aRL = aR/aL

!DO WHILE (tol .LEQ. 0.001)
!p_R(n+1)= p_R(n) -  
!tol = Pr-
!END DO
!END SUBROUTINE NEWTONSMETHOD

SUBROUTINE RANKINE(M,UX,RHO,A,PR,GAMMA,NX)
    USE types_vars
    REAL(DP), INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: M,gamma,A
	REAL(DP), INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: Pr,rho,ux
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: Mach, p21
    INTEGER(DP), INTENT(IN) :: nx
    ALLOCATE(MACH(NX),p21(Nx))
    	Mach(:) = ux(:)/a(:)
    	IF (maxval(MACH) .LE. -1) THEN
			M(:) = 0
		END IF
		IF ( maxval(MACH) .GE. -1 .AND. maxval(MACH) .LE. 1 ) THEN 
			M(:) = ((Mach(:)+1)/2)**2 
		END IF
		IF ( maxval(MACH) .GE. 1 ) THEN 
			M(:) = Mach(:) 
		END IF
		! in 2d, n(slope) * the free stream velocity 
		M(int(nx*0.5)) = 													   &
				SQRT(1/(2*gamma(1))*(Pr(int(nx*0.5)+1)/(Pr(1)))+(gamma(1)-1))
		Ux(int(0.5*nx)) = a(1)*M(int(nx*0.5)) - ux(1)
		Pr(int(0.5*nx)) = Ux(int(0.5*nx))
		!Pr(int(0.5*nx)+1) = 
    	   !Pr(int(0.5*nx)-1) =
    	
		Ux(int(0.5*nx)+1) = & 
				(2*a(1)*M(int(nx*0.5))**2)/((gamma(1)+1)*M(int(nx*0.5)+1))     &
				 - Ux(1)
		!Ux(int(0.5*nx)-1) = 
		
		rho(int(nx*0.5)+1) = 												   &
				rho(1)*((Pr(int(nx*0.5)+1)/Pr(1))+((gamma(1)-1)/(gamma(1)+1))) &
				*(1/(1+(gamma(1)-1)*(Pr(int(nx*0.5)+1)/Pr(1))/(gamma(1)+1)))
		!rho(int(nx*0.5)-1) = 												   &
				
		write(*,*) rho, ux, m
END SUBROUTINE RANKINE

SUBROUTINE EULER1D_MAC(U,Pr,rho,ux,uy,T,H,gamma,R,x,nx,nt,DX,DT)
    USE types_vars
    REAL(DP), INTENT(INOUT), DIMENSION(:,:), ALLOCATABLE :: U
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: F, Fp, Up, Un, F_p, F_n,F_pp,F_np
	REAL(DP), INTENT(IN), DIMENSION (:) :: x
    REAL(DP), INTENT(IN) :: DX
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: M,T,H,gamma,A_s,c,Mach
	REAL(DP), INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: Pr,rho,ux,uy,R
    REAL(DP) :: v,dt
    INTEGER(DP), INTENT(IN) :: nx, nt
    INTEGER(DP) :: n,i,k
	CHARACTER(LEN=20) :: pf
	ALLOCATE(F(nx,3),Fp(nx,3),Up(nx,3),Un(nx,3),F_pp(nx,3),F_np(nx,3))
    ALLOCATE(F_p(nx,3),F_n(nx,3),m(nx),C(NX),A_s(nx),Mach(nx))
    c(:) = SQRT(GAMMA*R*T)
    !CALL RANKINE(M,UX,RHO,C,PR,GAMMA,NX)
    CALL solvec(u,rho,ux,Pr,gamma,nx)
	CALL fluxvec(u,f,gamma,nx)
	CALL van_leer(f,Ux,M,F_p,F_n,rho,c,gamma,nx)
	DO k=1,3
    	v = (DT/DX)
		WRITE(*,*) 'DT:', DT
		WRITE(*,*) 'DX:', DX
		WRITE(*,*) 'MacCormack 1D Euler Stability ', v
		!Predictor Step
		DO n=1,3
			DO i=2,nx-1
           		Up(i,n) = u(i,n)-v*(F_p(i,n)+F_n(i+1,n)-(F_p(i-1,n)+F_n(i,n)))
			END DO
    	END DO
    	Up(1,:) = Up(2,:)
    	Up(nx,:) = Up(nx-1,:)
    	! Decode the solution into prims
		CALL DECODE_SOL(up,ux,pr,rho,c,gamma,nx)
		write(*,*) 'decode sol uX', UX
		CALL fluxvec(up,fp,gamma,nx)
		!CALL van_leer(Fp,Ux,M,Pr,F_pp,F_np,rho,c,gamma,ux,nx)
		!DO n=1,3
		!	DO i=2,nx-1
		!		!Correction Step
       	!		Up(i,n)= 0.5*(U(i,n)+Up(i,n) &
		!		-v*0.5*(F_np(i,n)+F_pp(i+1,n)-(F_np(i-1,n)+F_pp(i,n))))
		!	END DO
    	    !END DO
		! Decode the solution into prims
		CALL DECODE_SOL(up,ux,pr,rho,c,gamma,nx)
		! update the flux vector
		!CALL fluxvec(up,fp,gamma,nx)
		! split the flux
		!CALL van_leer(Fp,Ux,M,Pr,F_p,F_n,rho,c,gamma,ux,nx)
	END DO
	pf = '1d_euler_mac.dat'
    OPEN(1,file=pf,status='unknown')
	DO i=1,nx
		write(1,*) x(i), ' ', rho(i), ' ', ux(i), ' ', Pr(i), ' ' 
	END DO
END SUBROUTINE EULER1D_MAC

SUBROUTINE solvec(u,rho,ux,p,gamma,nx)
	USE types_vars
	REAL(DP), ALLOCATABLE,  INTENT(OUT), DIMENSION(:,:) :: u
	REAL(DP), INTENT(IN),DIMENSION(:) :: rho,ux,p,gamma
	INTEGER(DP), INTENT(IN) :: nx
	INTEGER(DP) :: i
	ALLOCATE(U(nx,3))
    ! Initialize the Solution Vector
		U(:,1) = rho(:)
		U(:,2) = rho(:)*ux(:)
		U(:,3) = P(:)/(gamma(:)-1.0) + 0.5*(rho(:)*ux(:)**2)
END SUBROUTINE solvec

SUBROUTINE fluxvec(u,f,gamma,nx)
	USE types_vars
	REAL(DP), ALLOCATABLE, INTENT(OUT), DIMENSION(:,:) :: F
	REAL(DP), ALLOCATABLE, INTENT(IN), DIMENSION(:,:) :: U
	REAL(DP), INTENT(IN), DIMENSION(:) :: gamma
	INTEGER(DP), INTENT(IN) :: nx
	INTEGER(DP) :: i
	ALLOCATE(F(nx,3))
    ! Initialize the Flux Vector 
	F(:,1) = U(:,2)
	F(:,2) = (gamma(:)-1.0)*u(:,3)+(3.0-gamma(:))*U(:,2)**2/(2.0*U(:,1))
	F(:,3) = gamma(:)*U(:,3)*U(:,2)/U(:,1) - (gamma(:)-1.0)*U(:,2)**3/(2.0*U(:,1)**2)
END SUBROUTINE fluxvec

SUBROUTINE source_vec(h,p,nx)
	USE types_vars
	REAL(DP), ALLOCATABLE, INTENT(OUT), DIMENSION(:,:) :: H
	REAL(DP), INTENT(IN), DIMENSION(:) :: P
	INTEGER(DP), INTENT(IN) :: nx
	INTEGER(DP) :: i
	ALLOCATE(H(nx,3))
    ! Initialize the Flux Vector 
	H(:,1) = 0
	H(:,2) = P(:)
	H(:,3) = 0
END SUBROUTINE source_vec

SUBROUTINE nozzle_area(a,x,nx)
	USE types_vars
	REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: A
	REAL(DP), INTENT(IN), DIMENSION(:) :: X
	INTEGER(DP), INTENT(IN) :: NX
	INTEGER(DP) :: I
	ALLOCATE(A(NX))
	DO i=1,nx
		A(i) = 1.398 + 0.347*TANH(0.8*X(i)-0.4)
	END DO
END SUBROUTINE nozzle_area

! Flux VECTOR Split - Van Leer #########################
SUBROUTINE VAN_LEER(F,ux,Mach,Fp,Fn,rho,a,gamma,nx)
    USE types_vars
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: F
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: Fp, Fn
    REAL(DP), ALLOCATABLE,  DIMENSION(:,:) :: Fe
    REAL(DP), INTENT(IN), DIMENSION(:) :: rho,ux,gamma,mach,a
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: M,Mp,Mn
    INTEGER(DP), INTENT(IN) :: nx
    INTEGER(DP) :: i,j,k
    ALLOCATE(Fe(nx,3),Fp(nx,3),Fn(nx,3),M(nx),Mp(nx),Mn(nx)) 
    M(:) = Mach(:)
	DO I=1,NX
		IF ( abs(M(I)) .LE. 1 ) THEN
		! Mach Split
		Mp(i) = 0.25*((M(i)+1)**2)
		Mn(i) = -0.25*((M(i)-1)**2)
		ELSE 
		! Mach split
		Mp(i) = 0.5*(M(i)+abs(M(i)))
		Mn(i) = 0.5*(M(i)-abs(M(i)))
		END IF 
	END DO
    !IF (maxval(MACH) .LE. -1) THEN
	!	M(:) = 0 
	!END IF
	!IF ( maxval(MACH) .GE. -1 .AND. maxval(MACH) .LE. 1 ) THEN 
	!	M(:) = ((Mach(:)+1)/2)**2 
	!END IF
	!IF ( maxval(MACH) .GE. 1 ) THEN 
	!	M(:) = Mach(:) 
	!END IF
	Fe(:,1) = rho(:)*a(:)*((M(:)+1)**2)*0.25
	Fe(:,2) = rho(:)*a(:)*((M(:)+1)**2)*0.25
	Fe(:,3) = rho(:)*a(:)*((M(:)+1)**2)*0.25
	
	Fp(:,1) = Fe(:,1)
	Fp(:,2) = Fe(:,2)*2*((gamma(:)-1)*Mp(:)*a(:)+2*a(:))/gamma(:)
	Fp(:,3) = Fe(:,3)*((gamma(:)-1)*Mp(:)*a(:)+2*a(:))**2/(2*(gamma(:)**2-1))
	Fn(:,1) = F(:,1) - Fp(:,1)
	Fn(:,2) = F(:,2) - Fp(:,2)
	Fn(:,3) = F(:,3) - Fp(:,3)
END SUBROUTINE van_leer

!##############Liou and Steffen (1993) Flux Vector Splitting #################
!#The flux is split into convective and pressure terms
!#The L or R primitive variables are designated by M_half .GEQ. 0 or 
!############### M_half < 0
!#############################################################################
SUBROUTINE AUSM(F,ux,Fc,Fp,rho,H,Pr,e,mach,gamma,a,nx)
    USE types_vars
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: F
    REAL(DP), INTENT(IN), DIMENSION(:) :: gamma
	REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: Fc,Fp
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Am, A_hat
	REAL(DP), INTENT(IN), DIMENSION(:) :: ux,mach,e,rho,H,Pr,a
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: M,Mp,Mn,M_half,Pp,Pn
	INTEGER(DP), INTENT(IN) :: nx
    INTEGER(DP) :: i,j,k
    ALLOCATE(M(nx),Fc(nx,3),Fp(nx,3),Pp(nx),Pn(nx))
	ALLOCATE(Mp(nx),Mn(nx),M_half(nx),Am(nx,3),A_hat(nx,3))
	M(:) = Mach(:) 
	DO I=1,NX
		IF ( abs(M(I)) .LE. 1 ) THEN
		! Pressure Split
		Pp(i) = 0.5*Pr(i)*((1+M(i))**2)
		Pn(i) = 0.5*Pr(i)*((1-M(i))**2)
		! Mach Split
		Mp(i) = 0.25*((M(i)+1)**2)
		Mn(i) = -0.25*((M(i)-1)**2)
		ELSE 
	    ! Pressure split
		Pp(i) = 0.5*Pr(i)*(M(i)+abs(M(i)))/M(i)
		Pn(i) = 0.5*Pr(i)*(M(i)-abs(M(i)))/M(i)
		! Mach split
		Mp(i) = 0.5*(M(i)+abs(M(i)))
		Mn(i) = 0.5*(M(i)-abs(M(i)))
		END IF 
		M_half(i) = Mn(i) + Mp(i)
		IF ( M_half(I) .GE. 0 ) THEN
			A_hat(i,1) = rho(i)*a(i)
			A_hat(i,2) = rho(i)*a(i)*ux(i)
			A_hat(i,3) = rho(i)*a(i)*H(i)		
		ELSE
		IF ( I+1 .LT. NX ) THEN
			A_hat(i,1) = rho(i+1)*a(i+1)
			A_hat(i,2) = rho(i+1)*a(i+1)*ux(i+1)
			A_hat(i,3) = rho(i+1)*a(i+1)*H(i+1)	
		ELSE IF ( I .EQ. NX ) THEN
			A_hat(i,1) = A_hat(nx-1,1)
			A_hat(i,2) = A_hat(nx-1,2)
			A_hat(i,3) = A_hat(nx-1,3)
		END IF
		END IF
	END DO
	DO j=1,nx
		Am(j,1) = M_half(j)*A_hat(j,1)
		Am(j,2) = M_half(j)*A_hat(j,2)
		Am(j,3) = M_half(j)*A_hat(j,3)
	END DO
	Fc(:,1) = Am(:,1)
	Fc(:,2) = Am(:,2)
	Fc(:,3) = Am(:,3)
	Fp(:,1) = 0
	Fp(:,2) = Pp(:) + Pn(:)
	Fp(:,3) = 0
END SUBROUTINE AUSM

! Flux Vector Decode ################################# 
SUBROUTINE DECODE_VL(F,nx,u,v,pr,rho,t,h,m,gamma)
	USE types_vars
	REAL(DP), INTENT(IN), DIMENSION(:,:) :: F
	REAL(DP), ALLOCATABLE, INTENT(OUT), DIMENSION(:):: u,v,rho,t,h,m,pr
	REAL(DP), INTENT(IN), DIMENSION(:) :: gamma
	INTEGER(DP), INTENT(IN) :: nx
	REAL(DP) :: r
	! Flux vector decoded into primitive variables
	! R - universal gas constant 8.3145 J/mol K
	! R - 1545.4 ft*lbf/lbmol*R 
	! Air : R = 287 m^2/(s^2 K) and Gamma = 1.4
	ALLOCATE(pr(nx),u(nx),v(nx),rho(nx),h(nx),m(nx),t(nx))
	v(:) = F(:,3)/F(:,1)
	u(:) = (gamma(:)/(1+gamma(:)))*(F(:,2)/F(:,1))  &
			+SQRT((gamma(:)/(gamma(:)+1))**2 &
		   - ((gamma(:)-1)/(gamma(:)+1))*(2*H(:)-v(:)**2))
	Pr(:) = F(:,2) - rho(:)*(u(:)**2)
	rho(:) = F(:,1)/u(:)
END SUBROUTINE DECODE_VL

! Flux Vector Decode ################################# 
SUBROUTINE DECODE_FLUX(F,nx,u,v,rho,t,h,m,r,gamma,units)
	USE types_vars
	REAL(DP), INTENT(IN), DIMENSION(:,:) ::F
	REAL(DP), ALLOCATABLE, INTENT(OUT), DIMENSION(:):: u,v,rho,t,h,m
	REAL(DP), INTENT(IN), DIMENSION(:) :: gamma
	INTEGER(DP), INTENT(IN) :: nx
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: a,b,c,p
	REAL(DP), INTENT(IN) :: r
	INTEGER(DP) :: i, units
	! Flux vector decoded into primitive variables
	! R - universal gas constant 8.3145 J/mol K
	! R - 1545.4 ft*lbf/lbmol*R 
	! Air : R = 287 m^2/(s^2 K) and Gamma = 1.4
	ALLOCATE(a(nx),b(nx),c(nx),p(nx),u(nx),v(nx),rho(nx),h(nx),m(nx),t(nx))
		! polynomial coefficients for calculating density
		A(:) = (F(:,3)**2)/(2*F(:,1))
		B(:) = F(:,1)*F(:,2)*gamma(:)/(1-gamma(:))
		C(:) = (F(:,1)**3)*(gamma(:)+1)/(2*(gamma(:)-1))
		rho(:) = (-b(:) + sqrt(b(:)**2 - 4*a(:)*c(:)))/2*a(:)
		! velocity components
		u(:) = F(:,1) / rho(:)
		v(:) = F(:,3) / F(:,1)
		! Pressure
		P(:) = F(:,2) - F(:,1)*u(:)
		! Temperature
		T(:) = P(:) / rho(:)*R
		! Total Enthalpy --
		! Steady, isoenergetic flow, constant H:
		! H = (gamma/(gamma-1))*P/rho + (u**2 + v**2)*0.5
		H(:) = SQRT(F(:,1)**2 + F(:,2)**2 + F(:,3)**2) + (P(:)/rho(:)) 
		M(:) = SQRT(u(:)**2+v(:)**2)/a
END SUBROUTINE DECODE_FLUX

SUBROUTINE WALL_PRESSURE(P,Omega,nu,dx,dy,dt,nx,ny,nt)
	USE types_vars
	REAL(DP), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: P
	REAL(DP), INTENT(IN), DIMENSION(:,:,:) :: Omega
	REAL(DP), INTENT(IN) :: nu,dy,dx,dt
	INTEGER(DP), INTENT(IN) :: ny,nx,nt
	INTEGER(DP) :: i,j,n
	do n=1,nt
		do i=2,nx-1
			do j=2,ny-1
				p(i+1,1,n) = -nu*2*dx*(-3*omega(i,1,n)+4*omega(i,2,n)       &
							 -omega(i,3,n))/(2*dy)
			end do
		end do
	end do
END SUBROUTINE

! Solution Vector Decode ################################# 
SUBROUTINE DECODE_SOL(u,ux,pr,rho,c,gamma,nx)
	USE types_vars
	REAL(DP), INTENT(IN), DIMENSION(:,:) :: U
	REAL(DP), ALLOCATABLE, INTENT(OUT), DIMENSION(:):: ux,pr,rho,c
	INTEGER(DP), INTENT(IN) :: nx
	REAL(DP), INTENT(IN),DIMENSION(:) :: gamma
	INTEGER(DP) :: i, units
    ALLOCATE(ux(nx),pr(nx),rho(nx),c(nx))
	! Solution Vector Decode ################################# 
	rho(:) = u(:,1)
	ux(:) = u(:,2)/u(:,1)
	pr(:) = (gamma-1.0)*(u(:,3)-0.5*(u(:,2)**2/u(:,1)))
	c(:) = SQRT(gamma*abs(PR(:))/rho(:))
END SUBROUTINE DECODE_SOL

! Roe Scheme ###############################################
SUBROUTINE BURGERS_ROE(U,UIC,nx,nt,DX,DT,b,c)
    USE types_vars
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: UIC
    REAL(DP), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: U
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: F
    INTEGER(DP), INTENT(IN) :: nx, nt 
    REAL(DP), INTENT(IN) :: DT, DX, b, c
    REAL(DP) :: v, r
    INTEGER(DP) :: n,j
    ALLOCATE(U(nx,nt),F(nx,nt))
    ! Initialize the IC's ############################################
    ! Initialize F = u2/2
    U(:,1) = UIC(:,1)
    F(:,1) = (U(:,1)**2)*b*0.5 + U(:,1)*c
    v = DT/DX
    r = 1/(DX*DX)
    write(*,*) 'Roes Method (Burgers) Stability: ' , v 
       DO J=1,nx
        DO N=1,nt-1
            IF (j+1 .GT. nx) THEN
                u(j,n+1) = u(j,n) - ((F(2,n)-F(j-1,n)) &
                - (c + (u(2,n) + u(j,n))*b*0.5)*(u(2,n)-u(j,n)) &
                + (c + (u(j,n) + u(j-1,n))*b*0.5)*(u(j,n)+u(j-1,n)))*v*0.5 &
                + r*(u(2,n)-u(j,n)*2 + u(j-1,n))
            END IF
            IF (j-1 .LT. 1 ) THEN
                u(j,n+1) = u(j,n) - ((F(j+1,n)-F(nx-1,n)) &
                - (c + (u(j+1,n) + u(j,n))*b*0.5)*(u(j+1,n)-u(j,n)) &
                + (c + (u(j,n) + u(nx-1,n))*b*0.5)*(u(j,n)+u(nx-1,n)))*v*0.5 &
                + r*(u(j+1,n)-u(j,n)*2 + u(nx-1,n)) 
            END IF
            IF (j+1 .LT. nx .AND. j-1 .GT. 1 ) THEN
                u(j,n+1) = &
                u(j,n) - ((F(j+1,n)-F(j-1,n)) &
                - (c + (u(j+1,n) + u(j,n))*b*0.5)*(u(j+1,n)-u(j,n)) &
                + (c + (u(j,n) + u(j-1,n))*b*0.5)*(u(j,n)+u(j-1,n)))*v*0.5 &
                + r*(u(j+1,n)-u(j,n)*2 + u(j-1,n)) 
            END IF
             F(j,n+1) = (U(j,n+1)**2)*b*0.5 + U(j,n+1)*c
        END DO
   END DO
END SUBROUTINE BURGERS_ROE

! Linearized Burger's Exact Solution
SUBROUTINE BE_LINEAR(U,X,T,mu,nt,nx,DX,DT,c)
    USE types_vars
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: U
    REAL(DP), INTENT(IN) :: DT,DX,MU,C
    REAL(DP), INTENT(IN), DIMENSION(:) :: X, T
    INTEGER(DP) :: nt, nx
    INTEGER(DP) :: J,N
    ALLOCATE(U(nx,nt))
    DO j=1,nx
        DO n=1,nt
            u(j,n) = exp(-4*pi*pi*mu*t(n))*sin(2*pi*(x(j)-c*t(n)))
        END DO
    END DO
END SUBROUTINE BE_LINEAR


! For solving non-linear Poisson Equation using Newton's method
! The solution vector needs to be mapped into a 1D array :
SUBROUTINE MAP_SOLN(A,U,FD,NX,NY)
	USE types_vars
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: U
	REAL(DP), INTENT(IN), DIMENSION(:) :: FD
	REAL(DP), ALLOCATABLE, INTENT(OUT), DIMENSION(:,:) :: A
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: UV
	REAL(DP)  :: H
    INTEGER(DP), intent(in) :: nx,ny
    INTEGER(DP) :: I,J,K,nm2,n2
	ALLOCATE(UV(NX),A(NX,NY))
	nm2 = nx-1	
	h = nx*ny
	k=1
	n2 =1 
	do j=2,nx-1
		do i=2,nx-1		
		uv(k) = u(i,j) 
		k=k+1
		end do 
	end do
	do j = 1, n2
		do i = 1, n2
			if ( j .eq. i ) then
				a(i,j) = -4.0-(h**2)*fd(i)
			else if (j.eq.(i+1).and.mod(i,nm2)/=0.or.j.eq.(i-1).and.mod(j,nm2)/=0) then
				a(i,j) = 1.0
			else if (j.eq.(i-nm2).or.j.eq.(i+nm2)) then
				a(i,j) = 1.0
			else
				a(i,j) = 0.0
			end if 
		end do
	end do
END SUBROUTINE MAP_SOLN

! 5-Point Scheme FDE Poisson Matrix Assembler 
SUBROUTINE P_MATRIX(A,m,n,NX)
	USE types_vars
	! A = [M,N]
	REAL(DP), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: A
	! Matrix Size
	INTEGER(DP), INTENT(IN) :: NX,m,n
	INTEGER(DP) :: I,J,K,S
	s = 0.5*(M)
	ALLOCATE(A(M,M))
	A(:,:) = 0
	DO J=2,M-1
		A(j,j) = 4
		IF (J+1 < M-1) THEN
			A(j+1,j) = -1
			A(j,j+1) = -1
		END IF
	END DO
	DO j=2,M-S-1
		A(j+s,j) = -1
		A(j,j+s) = -1
	END DO
END SUBROUTINE P_MATRIX
END module subs
