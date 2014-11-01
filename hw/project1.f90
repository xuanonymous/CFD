PROGRAM project1
    USE SUBS
    USE TYPES_VARS
    ! Problem 1 REAL(DP), PARAMETER :: A = 0, B = 4, DX=0.1
    ! Problem 2 
    REAL(DP), PARAMETER :: A=0,B=40,C=1,V=1,dx=1,dt=c*dx*V,t1=0,t2=18
    REAL(DP) :: H,L,tH,TL,N=1
    INTEGER :: NPTSX,NPTST,I,J,NEQN,G
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: X, T, XX, QQ, AA, BB, CC
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: UIC, U, LWU, UEXACT, UBM
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI, M, D
    CHARACTER(LEN=20) :: plot_file
        ! Select Case for V=? that will increment
        ! file names to increment each time
        ! throw in an IF statement to see if
        ! File name xists, if it does, increment,
        ! If it doesnt, create it
        CALL space_initial(L,A,B,DX,H,NPTSX)
        CALL time_initial(TL,T1,T2,DT,tH,NPTST)
        CALL xgridgen(A,H,NPTSX,X)
        CALL tgridgen(T1,TH,NPTST,T)
        CALL WAVE_IC(DX,DT,V,C,X,NPTSX,UIC,N)
        CALL LAX_WAVE(DX,DT,V,C,X,T,UIC,U,NPTSX,NPTST)
        CALL LW_WAVE(DX,DT,V,C,X,T,UIC,LWU,NPTSX,NPTST)
        CALL tecplot(u,plot_file,X,U,nptsx,nptsx,1,nptsx,1,nptsx,40)
        CALL tecplot(lwu,plot_file,X,LWU,nptsx,nptsx,1,nptsx,1,nptsx,40)
        ALLOCATE(UEXACT(NPTSX,NPTST))
        DO j=1,NPTSX
                DO i=1, NPTST
                        uexact(j,i) = sin(2*n*3.14*x(j)/40)*t(i) 
                END DO
        END DO
        CALL tecplot(uexact,plot_file,X,uexact,nptsx,nptsx,1,nptsx,1,nptsx,40)
        CALL BOXMETHOD(DX,DT,V,C,X,NPTSX,NPTST,NEQN,UBM,N)
        CALL tecplot(UBM,plot_file,X,UBM,NPTSX,NPTSX,1,NPTSX,1,NPTSX,40)
        !WRITE(*,*) 'x values', X(:)
        !WRITE(*,*) UX(:)
        !WRITE(*,*) UX12(:)
        !WRITE(*,*) UX123(:)
END PROGRAM project1
