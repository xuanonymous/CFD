PROGRAM uvdis
    USE subs
    USE types_vars
    USE matrixsolv
    USE gnufor2
    REAL(DP)                                :: T1,T2
    REAL(DP)		                    :: DX,DY,DZ,DT,aq,aw,rq,tq,p
    REAL(DP)     	                    :: Lx,Ly,Lz,x1,x2,y1,y2,z1,z2
    INTEGER(dp)  		            :: nx,ny,nt,nz,I,J,K
    REAL(DP), ALLOCATABLE, DIMENSION(:)     :: X,Y,Z,T,N,G,RHO
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: R
    CHARACTER(LEN=20) :: pf
	P = 100
	aw = -log(0.70)
	aq = -log(0.95)/0.2
	rq = 1.25
	tq = 0.2
	z1 = 0
	z2 = 5
	y1 = 0
	y2 = 5
	x1 = 0
	x2 = 5
	dX = (x2-x1)/10
	dY = (y2-y1)/10
	dZ = (z2-z1)/10
	CALL gridgen(Lx,x1,x2,DX,X,nx)
	CALL gridgen(Ly,y1,y2,DY,Y,ny)
	CALL gridgen(Lz,z1,z2,DZ,Z,nz)
	ALLOCATE(RHO(nx),g(3),n(100))
    	g(:) = (/(k, k=1,3)/)
	n(:) = (/(k, k=-75,75)/)
	DO i=1,nx
	DO j=1,ny
	DO k=1,nz
		rho(i) = SQRT(x(i)**2+y(j)**2+z(k)**2)
	END DO
	END DO
	END DO
	pf = 'uv_data.dat'
   	OPEN(1,file=pf,status='unknown')
	DO i=1,nx
	  DO j=1,ny
            DO k=1,nz
	      write(1,*) x(i),',',  y(j),',',  z(k),',', rho(i)
	    END DO
	  END DO
	END DO
    DEALLOCATE(x,rho)
END PROGRAM uvdis

