PROGRAM nozzle2
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    REAL(DP) :: T1,T2,ALPHA,KAPPA,A,B,B1,N1 
    REAL(DP) :: DX,DY,DT
    REAL(DP) :: H,L,tH,TL,tc,C
    integer(dp)  :: nx,ny,nt,I,J,K,N,counter
    ! Mesh Variables ####################################################
    REAL(DP), ALLOCATABLE, DIMENSION(:)    :: X, Y, T
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: UIC, U, R
	CHARACTER(LEN=20) :: pf
	! Nozzle Conditions #################################################
	DY = 0
	A =  0
	B =  10
	nx = 100
	DX = (B-A)/nx
	nu = 0
	CFL = 0.05
	CALL gridgen(L,A,B,DX,X,nx)
	ALLOCATE(R(nx,ny))
	!################## Initial Conditions ###############################
	gma(:) = 1.4
	Rg(:)=1716
	cv= (1/(gma(1)-1))*Rg(1)
	DO i=1,nx
		DO j=1,ny
			R(i,j) = SQRT((0.5+i)**2 + (0.5+j)**2)
    END DO 

	pf = 'n2.dat'
    OPEN(1,file=pf,status='unknown')
	DO i=1,nx
		write(1,*) x(i), ' ', rho(i), &
					' ', ux(i), &
					' ', P(i), ' ', Mc(I), &
					' ', Area(i),' ',Tp(i)
	END DO
    DEALLOCATE(x,rho,ux,p,Mc,area,tp)
END PROGRAM nozzle2

