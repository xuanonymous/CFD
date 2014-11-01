MODULE linear
    integer, parameter, public :: kind=selected_real_kind(10)
    public :: solve

CONTAINS
    SUBROUTINE solve(a, piv_tol, b, ok)
    !arguments
    real(kind), intent(inout), dimension(:,:) :: a
    
    !The matrix a.
    real(kind), intent(in) :: piv_tol
    
    !Smallest acceptable pivot.
    real(kind), intent(inout), dimension(:) :: b
    
    !The rhs vector on entry. Overwritten by solution.
    logical, intent(out) :: ok
    !True after a successful entry and false otherwise.
    ! Local variables
    integer :: i !Row index
    integer :: j !Column index
    integer :: n !Matrix order
    real(kind), dimension(size(b)) :: row
    
    !Automatic array needed for workspace.
    real(kind) :: element !Workspace variable.
    
    n = size(b)
    ok = size(a, 1) == n .and. size(a,2) == n
    if (.not.ok) then
        return
    end if
    
    do j = 1, n
    ! Update elements in column j.
        do i = 1, j-1
            a(i+1:n, j) = a(i+1:n, j) - a(i,j)*a(i+1:n,i)
        end do
    ! Find pivot and check its size (using maxval just to obtain
    ! a scalar).
    	i = maxval(maxloc(abs(a(j:n, j)))) + j - 1
    !maxval and maxloc
        if (abs(a(i,j)) < piv_tol) then
    		ok=.false.
        return
    end if
    
    ! If necessary apply row interchange
        if (i/=j) then
            row = a(j,:); a(j,:) = a(i,:); a(i,:) = row
            element = b(j); b(j) = b(i); b(i) = element
        end if

    !Compute elements j+1 : n of j-th column.
        a(j+1:n,j) = a(j+1:n, j) / a(j,j)
    end do

    !Forward substitution
    do i = 1, n-1
        b(i+1:n) = b(i+1:n) - b(i)*a(i+1:n,i)
    end do
    
    !Back substitution
    do j = n, 1, -1
        b(j) = b(j)/a(j,j)
        b(1:j-1) = b(1:j-1) - b(j)*a(1:j-1,j)
    end do
    end subroutine solve
END MODULE linear


