! Matrix Solvers
MODULE matrixsolv
IMPLICIT none
CONTAINS

! 'Cyclic Thomas Algorithm' :
SUBROUTINE thomas_per(neqn,bb,aa,cc,xx,qq)
USE types_vars
USE subs
    REAL(DP), INTENT(IN), DIMENSION(:) :: qq, aa, bb, cc
    INTEGER, INTENT(IN) :: neqn
    INTEGER :: I
    REAL(DP), INTENT(OUT), DIMENSION(:) :: xx
    REAL(DP), DIMENSION(NEQN) :: QQ2, X1, X2
    ! AA(1) CC(1) 0 ----------0 BB(1) QQ(1)
    ! BB(1) AA(2) CC(2) 0 ---------0  QQ(2)
    ! 0     BB(2) AA(3) CC(3) ---- 0  QQ(3)
    !                 AA(N-1) CC(n-1) QQ(N-1)
    !                 BB(n-1)   AA(N) QQ(N)
    !
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
        INTEGER, INTENT(IN) :: neqn
        INTEGER :: i
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
        INTEGER, INTENT(IN) :: NPTSX, NPTST
        INTEGER  :: T,J,NEQN,NP
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

! Enter your matrix in row major order as usual
! Plug it into this subroutine
! Continue your operations as usual
SUBROUTINE MATRIX(A)
USE types_vars
	USE subs
	REAL(DP),INTENT(inout),DIMENSION(:,:),ALLOCATABLE :: A
	A = TRANSPOSE(RESHAPE(A,SHAPE(A)))
END SUBROUTINE

! Conjugate Gradient Method ####################################
SUBROUTINE CG(A,x,b,r,n,alpha)
USE types_vars
USE subs
INTEGER, INTENT(IN) :: N
REAL(DP),INTENT(IN),DIMENSION(:,:), ALLOCATABLE :: b
REAL(DP), INTENT(INOUT), DIMENSION(:,:), ALLOCATABLE :: x, a
REAL(DP) :: alpha, r
INTEGER :: i
!CALL matrix(A)
r = b - MATMUL(A,x)
alpha = &
MATMUL(TRANSPOSE(r),r)/(TRANSPOSE(r)*A*r)
i = 0
DO
    IF (abs(sum(r)) .LT. 10**(-3)) EXIT
    i = i + 1
    X=X+alpha*r
    R=R-alpha*r
    alpha = &
    MATMUL(TRANSPOSE(r),r)/(MATMUL(TRANSPOSE(r),MATMUL(A,r)))
END DO
END SUBROUTINE CG
END MODULE matrixsolv


