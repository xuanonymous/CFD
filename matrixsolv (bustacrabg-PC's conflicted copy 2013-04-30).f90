! Matrix Solvers
MODULE matrixsolv
IMPLICIT none
CONTAINS

! 'Cyclic Thomas Algorithm' :
SUBROUTINE thomas_per(neqn,bb,aa,cc,xx,qq)
USE types_vars
USE subs
    REAL(DP), INTENT(IN), DIMENSION(:) :: qq, aa, bb, cc
    INTEGER(DP), INTENT(IN) :: neqn
    INTEGER(DP) :: I
    REAL(DP), INTENT(OUT), DIMENSION(:) :: xx
    REAL(DP), DIMENSION(NEQN) :: QQ2, X1, X2
    ! AA(1) CC(1) 0 ----------0 BB(1) QQ(1)
    ! BB(1) AA(2) CC(2) 0 ---------0  QQ(2)
    ! 0     BB(2) AA(3) CC(3) ---- 0  QQ(3)
    !                 AA(N-1) CC(n-1) QQ(N-1)
    !                 BB(n-1)   AA(N) QQ(N)
    qq2 = 0.0d0
    qq2(1) = -bb(neqn)
    qq2(neqn-1) = -cc(neqn-1)

    call thomas(neqn-1,bb(1:neqn-1),aa(1:neqn-1),cc(1:neqn-1), &
            x1(1:neqn-1),qq(1:neqn-1))
    call thomas(neqn-1,bb(1:neqn-1),aa(1:neqn-1),cc(1:neqn-1), &
            x2(1:neqn-1),qq2(1:neqn-1))
    xx(neqn) = (qq(neqn)-cc(neqn)*x1(1)-bb(neqn-1)*x1(neqn-1))/ &
            (aa(neqn)+cc(neqn)*x2(1)+bb(neqn-1)*x2(neqn-1))
    do i=1,neqn-1
            xx(i) = x1(i)+x2(i)*xx(neqn)
    end do
END SUBROUTINE thomas_per

SUBROUTINE thomas(neqn,bb,aa,cc,xx,qq)
    	USE types_vars
		USE subs
        REAL(DP), INTENT(IN), DIMENSION(neqn) :: qq, aa, bb, cc
        INTEGER(DP), INTENT(IN) :: neqn
        INTEGER(DP) :: i
        REAL(DP), INTENT(OUT), DIMENSION(neqn) :: xx
        REAL(DP), DIMENSION(neqn) :: l, u, d, yy
        !
        ! aa(1) cc(1) 0 ----------0
        ! bb(1) aa(2) cc(2) 0 -----0
        ! 0 ! bb(2) aa(3) cc(3)---0
        !  aa(n-1)cc(n-1)
        !  bb(n-1)aa(n)
        ! LU decomposition
        d(1) = aa(1)
        u(1) = cc(1)
        do i = 1, neqn-2
                l(i)   = bb(i)/d(i)
                d(i+1) = aa(i+1) - l(i)*u(i)
                u(i+ 1) = cc(i+1)
        end do
        l(neqn-1) = bb(neqn-1)/d(neqn-1)
        d(neqn)   = aa(neqn) - l(neqn-1)*u(neqn-1)
        ! Forward substituion Lyy=qq
        yy(1) = qq(1)
        do i = 2, neqn
                yy(i) = qq(i) - l(i-1)*yy(i-1)
        end do
        ! Backward substitution Ux=yy
        xx(neqn) = yy(neqn)/d(neqn)
        do i = neqn-1,1,-1
                 xx(i) = (yy(i) - u(i)*xx(i+1))/d(i)
        end do
END SUBROUTINE thomas

SUBROUTINE boxmethod(DX,DT,V,C,X,NPTSX,NPTST,NEQN,UBM,UIC)
  USE types_vars
  USE subs
  REAL(DP), INTENT(IN) :: DX,DT,C,V
  REAL(DP), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: X
  REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: UBM
  REAL(DP), INTENT(IN), ALLOCATABLE, DIMENSION(:,:) :: UIC
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: BB,AA,CC,XX,QQ
  INTEGER(DP), INTENT(IN) :: NPTSX, NPTST
  INTEGER(DP)  :: T,J,NEQN,NP
  ! Diagonal
  ALLOCATE(AA(NPTSx))
  DO j=1,nptsx
    AA(j)=1-V
  END DO
   ! Lower Diagonal
  ALLOCATE(BB(NPTSx))
  DO j=1,nptsx
   BB(J)=0
  END DO
   ! Upper Diagonal
  ALLOCATE(CC(NPTSx))
  DO J=1,nptsx
    CC(J)=1+V
  END DO
  NEQN=NPTSx
  ALLOCATE(UBM(NPTSX,NPTST),XX(NPTSx),QQ(NPTSx))
  DO J=1,NPTSX
    UBM(J,1)=UIC(J,1)       
    QQ(J) = UIC(J,1)
  END DO
  DO J=2, NPTSX-1
    DO t=1,NPTST-1
      CALL thomas_per(neqn,bb,aa,cc,xx,qq)
        UBM(J,T) = XX(J)
        QQ(J)=QQ(J)-(bb(J)/cc(j))*qq(j-1)
    END DO
  END DO
END SUBROUTINE boxmethod

! Conjugate Gradient Method(rougher) ####################################
! Method of Steepest Descent
SUBROUTINE CG(A,x,b,n)
    USE types_vars
    USE subs
    INTEGER(DP),  INTENT(IN) :: N
    REAL(DP), INTENT(IN), DIMENSION(:), ALLOCATABLE :: b
    REAL(DP), INTENT(INOUT), DIMENSION(:), ALLOCATABLE :: x
	REAL(DP), DIMENSION(:), ALLOCATABLE :: r,d,q,z,w,lm,c,c_1,alpha
    REAL(DP), INTENT(IN), DIMENSION(:,:), ALLOCATABLE :: A
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dd,rr,xx
    REAL(DP) :: tol, h, temp
    INTEGER(DP) :: i,j,m,y
    ALLOCATE(r(n),d(n),q(n),c(n),c_1(n),z(n),w(n),alpha(n),lm(n))
	ALLOCATE(rr(1,n),dd(1,n),xx(1,n)) 
    ! Check to see if matrix is positive definite or not
	! if it is, run the algorithm
    ! x.'Ax > 0 x .NEQ. 0
	! use m to keep track of number of iterations
	! to convergence
	m = 0
	C_1 = DOT_PRODUCT(MATMUL(A,x),x)
	h = maxval(SQRT(C_1))
	r = MATMUL(A,x) - b

	IF (h .GT. 0) THEN
    	DO 	
			m = m + 1
			y = 0 
			d(1:n) = -r(1:n)
			lm(1:n) = -DOT_PRODUCT(R,DD(1,1:n))/DOT_PRODUCT(dd(1,1:n),MATMUL(A,d))
			!lm(1:n) = -(dd(1,1:n)*r(1:n))/(dd(1,1:n)*MATMUL(A(1:n,1:n),d(1:n)))
			x = x + lm*d
			r = r + lm*MATMUL(A,d)
			! Check for convergence IF sqrt(r_t*r) < tolerance then
            ! algorithm has reached convergence
            temp = SQRT(DOT_PRODUCT(RR(1,1:n),R(1:n)))
            tol = 0.00000000000000000000000000000000000000000001
            IF (temp <= tol) THEN
                write(*,*) x
                write(*,*) "no of iterations: ", m
                RETURN
            END IF
            ! Exit if NaN's are encountered
            DO I=1,n
                IF ( X(i)/X(i) /= 1) THEN 
                    Y = Y + 1
                END IF
                IF ( Y .eq. 3 ) RETURN 
            END DO
        END DO
	END IF
END SUBROUTINE CG

! SOR ##########################################################
SUBROUTINE SOR(A,b,u,h,nx,ny,beta,omega,tol)
    USE types_vars
    USE subs
	REAL(dp), ALLOCATABLE, INTENT(out), DIMENSION(:,:) :: u
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: A,B
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: U1, Residual,r
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Tr
    INTEGER(DP), INTENT(IN) :: nx, ny
    REAL(DP), INTENT(IN) :: tol, beta,h
	REAL(dp), INTENT(out) :: omega
	REAL(dp) :: mu

    INTEGER(DP) :: m, n, i, j, k, num1
    N = SIZE(A(:,1))
    M = SIZE(A(1,:))
    ! U1 temporary storage of values
    ALLOCATE(U1(m,n),Residual(m,n),Tr(m),R(m,n),u(m,n))
    ! Calculate the Point-Jacobi Eigenvalue :
	mu = 0.25*((cos(pi/M) + cos(pi/N))**2)
	omega = ABS(2/(1+sqrt(1-mu**2)))
	write(*,*) 'omega :', omega, 'mu ', mu
	U(:,:) = A(:,:)
	U1(:,:) = U(:,:)
	R(:,:) = b(:,:)
    k = 1
    num1 = 0
    DO
    	k = k + 1
        DO i=2,m-1
            DO j=2,n-1
                U1(i,j) = ((u1(i-1,j)+u(i+1,j)+u1(i,j-1)+u(i,j+1))+h*h*r(i,j))*0.25
                U1(i,j) = (1-omega)*U(i,j) + omega*U1(i,j)
				Residual(i,j) = U(i,j) - U1(i,j)
				U(i,j) = U1(i,j)
				IF ( ABS(Residual(i,j)) .LE. tol ) THEN
					write(*,*) "GS-SOR Iterations: ", k
				RETURN
	    END IF
		  END DO
        END DO

    END DO
END SUBROUTINE SOR

! Jacobi Iterative Method(smoother) ######################################
SUBROUTINE JI(A,x,b)
    USE types_vars
    USE subs
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: A
    REAL(DP), INTENT(INOUT), DIMENSION(:) :: X,b
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: U, Un, lm, f
    INTEGER(DP) :: m, n, k, num1
    INTEGER(DP) :: i,j
    k=0
    num1 = 0
    N = SIZE(A(:,1))
    M = SIZE(A(1,:))
    
    ALLOCATE(Un(m,n),U(m,n),lm(m,n),f(m,n))
    DO I=1,M
        DO J=1,N
            U(j,i) = A(j,i) 
        END DO
    END DO
    DO i=1,M
        DO j=1,N
            lm(j,i) =  0.5*(cos(i*pi/M)+cos(j*pi/N))
        END DO
    END DO
    DO
        k = k + 1
        DO i=1,m
            DO j=1,n-1
				IF ( i-1 .LT. 1 .AND. j-1 .GT. 1 ) THEN
                		Un(j,i)= (U(j,M-1)+U(j,i+1)+U(j-1,i)+U(j+1,i)-f(j,1))*0.25
				END IF
				IF ( j-1 .LT. 1 .AND. i-1 .GT. 1 .AND. i+1 .LT. N ) THEN
                		Un(j,i)= (U(j,i-1)+U(j,i+1)+U(N-1,i)+U(j+1,i)-f(j,1))*0.25
				END IF
				IF ( i-1 .LT. 1 .AND. j-1 .LT. 1 ) THEN
                		Un(j,i)= (U(j,N-1)+U(j,i+1)+U(M-1,i)+U(j+1,i)-f(j,1))*0.25
				END IF
				IF (i-1 .GT. 1 .AND. i+1 .LT. N .AND. j-1 .GT. 1) THEN
					Un(j,i)= (U(j,i-1)+U(j,i+1)+U(j-1,i)+U(j+1,i)-f(j,1))*0.25
				END IF
            END DO
        END DO
    END DO 
END SUBROUTINE JI
END MODULE MATRIXSOLV
