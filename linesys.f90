PROGRAM linesys
        USE linear
        USE types_vars
        USE subs
        INTEGER :: i,n,p
        REAL(DP),ALLOCATABLE :: a(:,:),j(:,:),b(:),x1(:),x2(:),x3(:),c(:,:),v(:),x(:)
        LOGICAL :: ok
        print *, 'Matrix Order?'
        read *,n
        allocate (a(n,n),j(n,n),b(n),x1(n),x2(n),x3(n), c(n,n), v(n),x(n))
        do i=1,n
                write(*, '(a,i2,a)') ' Elements of row ', i, ' of a?'
                read *, a(i,:)
                write(*, '(a,i2,a)') 'Component ', i,' of b?'
                read *, b(i)
        end do
        ! Analytical Solution to Jacobian of given 
        ! Vector valued functions
        c(1,1) = x1(:)**3
        c(1,2) = -2*x2(:)
        c(1,3) = -2
        c(2,1) = x1(:)**3
        c(2,2) = -6*x3(:)**2
        c(2,3) = 7
        c(3,1) = x2(:)*x3(:)*x3(:)
        c(3,2) = -1
        c(3,3) = 0
        j(1,1) = 3*x1(:)**2
        j(1,2) = -2
        j(1,3) = 0
        j(2,1) = 3*x1(:)
        j(2,2) = -10*x3(:)
        j(2,3) = 0
        j(3,1) = 0
        j(3,2) = x3(:)**2
        j(3,3) = 2*x2(:)*x3(:)
        call solve(a, maxval(abs(a))*1.0e-10,b,ok)
        if (ok) then
                write(*, '(/,a,/,(5f12.4))') ' Solution is', b
        else
                print*, 'The matrix is singular'
        end if
        ! Newton's Method to Solve Non-Linear
        ! System of Equations
        DO i=1,n
                IF (abs(maxval(v))<5e-04) THEN
                DO p=1,n
                v(i)=-c(x(i),p)      
                END DO
                END IF
        END DO
        DO i=1,n
                x(I+1)=x(i)+v(i)
        END DO    
end program linesys 
