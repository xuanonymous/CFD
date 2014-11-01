!-------------------------------------------------------------------------------
! - Simple program for calculating formation free energies (configurational entropy 
! not accounted for) of intrinsic point defects in Si, using the Bazant empirical 
! potential. 
!
! - Originally designed for the "Computational Physics" class taught by Stefan Goedecker
! at the Department of Physics, University of Basel, Switzerland
!
! - For background on the free-energy calculations of intrinsic point defects in silicon
! see: O. K. Al-Mushadani and R. J. Needs, Phys. Rev. B 68, 235205 (2003).
!
!
! (c)2004  E. Penev
!
!-------------------------------------------------------------------------------
program si_point_def

  implicit none

  integer               :: i, ierr, j, n, n0, natoms
  integer, parameter    :: nsd_max   = 1000            ! max. num of SD steps
  integer, parameter    :: nt        = 21              ! num of T sampling points

  real(8), parameter    :: a0        = 5.42981d0       ! lattice constant of bulk Si in A
  real(8), parameter    :: m_si      = 28.0855d0       ! atomic mass of Si [amu]
  real(8), parameter    :: feps      = 1.d-5           ! conv. crit. for forces

  ! some constants and conversion factors
  real(8), parameter    :: k_to_ev     = 8.617343d-5         !    K --> eV
  real(8), parameter    :: invcm_to_ev = 1.23984191d-4       ! 1/cm --> eV
  real(8), parameter    :: ev_to_kj    = 4.184d0 * 23.061d0  ! eV -> kcal/mol -> kJ/mol
  real(8), parameter    :: pi          = 3.141592653589793d0 ! pi
  real(8), parameter    :: c           = 299792458.0d0       ! speed of light
  real(8), parameter    :: to_wavenumbers = 1.0d+11 / (2.0d0 * pi * c)

  real(8), dimension(3) :: alat = 3.d0 * (/a0, a0, a0/) ! 3x3x3 bulk block section
  real(8)               :: t1, t2, fnrm, rand, count
  real(8)               :: e0, e1, e2
  real(8)               :: fvib0(nt), fvib1(nt), fvib2(nt)
  real(8)               :: t                     ! temperature
  character(20)         :: file1 = 'si_ffcd'     ! input atomic coordinates for ffcd
  character(20)         :: file2 = 'si_vacancy'  ! input atomic coordinates for vacancy

  ! global allocatable arrays
  real(8), allocatable  :: atmcrd(:,:)           ! atomic coordinates
  real(8), allocatable  :: f(:,:), hess(:,:)     ! forces and hessian matrices
  real(8), allocatable  :: eig(:)                ! normal modes
  !-------------------------------------------------------------------------------

  call cpu_time ( t1 )

  ! initialize bazant() call counter
  count = 0.0d0

  !-------------------------------------------------------------------------------
  ! 0.  Perfect bulk Si
  !-------------------------------------------------------------------------------
  write (6,'(72("-"))')
  write (6,'(A)') '>>> Doing perfect bulk Si'
  write (6,'(72("-"))')

  call xyz_read ( 'si_bulk216' )
  allocate ( f(1:3,1:natoms) )

  ! Save the number of atoms in the perfect Si supercell
  n0 = natoms

  call gradient ( e0, f(:,:), atmcrd(:,:), print = .true. )

  allocate ( eig(1:3*natoms), stat=ierr )
  if (ierr /= 0) stop '** main: allocation failed...'

  ! Calculate normal modes of perfect Si
  call get_modes ( eig(:), atmcrd(:,:), 1.0d-4, alat(:), print=.true. )

  ! Set initial temperature in [K]
  t = 1.0d0
  
  ! Calculate energy for different temperatures
  do i=1,nt
     call get_fvib ( fvib0(i), t, eig(:) )
     print *, t, fvib0(i)
     t = t + 100.d0
  end do

  !-------------------------------------------------------------------------------
  ! I. Fourfold coordinated point defect (FFCD)
  !-------------------------------------------------------------------------------
  write (6,*) ''
  write (6,'(72("-"))')
  write (6,'(A)') '>>> Doing Si fourfold coordinated point defect (FFCD)'
  write (6,'(72("-"))')

  call xyz_read ( trim (file1) )

  if ( allocated (f) ) deallocate ( f )
  allocate ( f(1:3,1:natoms) )

  if ( allocated (eig) ) deallocate ( eig )
  allocate ( eig(1:3*natoms) )

  call gradient ( e1, f(:,:), atmcrd(:,:), print = .true. )

  ! Steepest descent optimization ---------------------------
  write (6,'(A)') '>>> Doing random shift of the coordinates'

  ! First do small random shift to avoid symmetry lock-in
  call random_seed ()
  do i=1,natoms
     do j=1,3
        call random_number ( rand )
        atmcrd(j,i) = atmcrd(j,i) + ( rand - 0.5d0 )/2.0d0
     end do
  end do

  write (6,'(72("-"))')
  write (6,'(A)') '>>> Steepest descent optimization: step = 1/50'
  write (6,'(72("-"))')
  write (6,'(A)') '>>>    it      Etot         rmsf'

  do i=1,nsd_max
     
     ! commit SD step
     atmcrd(:,:) = atmcrd(:,:) + 2.0d-2 * f(:,:)

     ! calculate new gradient
     call gradient ( e1, f(:,:), atmcrd(:,:) )

     ! calculate rms force
     fnrm   = sqrt ( sum(f(:,:)*f(:,:)) / (3*natoms) )

     ! control print
     if ( mod (i,20) == 0 .or. i==1 ) write (6,'(A,I5,X,2F12.6)') '>>> ', i, e1, fnrm

     ! check for convergence
     if ( fnrm < feps ) exit

  end do

  write (6,'(A,I5,2(A,F12.6))') '>>> Convergence reached: it= ', i, &
                                ', Etot= ', e1, ', rmsf= ', fnrm

  write (6,'(72("-"))')
  write (6,*) ''

  ! Write Vsim format
  call write_vsim   ( trim (file1)//'_relaxed', alat(:), atmcrd(:,:) )
  !----------------------------------------------------------

  ! Calculate normal modes of "defective" Si
  call get_modes ( eig(:), atmcrd(:,:), 1.0d-4, alat(:), print=.true. )

  ! Set initial temperature in [K]
  t = 1.0d0
  
  ! Calculate energy for different temperatures
  do i=1,nt
     call get_fvib ( fvib1(i), t, eig(:) )
     print *, t, fvib1(i)
     t = t + 100.d0
  end do

  !-------------------------------------------------------------------------------
  ! II. Vacancy defect
  !-------------------------------------------------------------------------------
  write (6,*) ''
  write (6,'(72("-"))')
  write (6,'(A)') '>>> Doing Si vacancy'
  write (6,'(72("-"))')

  call xyz_read ( trim (file2) )

  if ( allocated (f) ) deallocate ( f )
  allocate ( f(1:3,1:natoms) )

  if ( allocated (eig) ) deallocate ( eig )
  allocate ( eig(1:3*natoms) )

  call gradient ( e2, f(:,:), atmcrd(:,:), print = .true. )

  ! Steepest descent optimization ---------------------------
  write (6,'(A)') '>>> Doing small random shift of the coordinates'

  ! First do small random shift to avoid symmetry lock-in
  do i=1,natoms
     do j=1,3
        call random_number ( rand )
        atmcrd(j,i) = atmcrd(j,i) + ( rand - 0.5d0 )/2.0d0
     end do
  end do

  write (6,'(72("-"))')
  write (6,'(A)') '>>> Steepest descent optimization: step = 1/50'
  write (6,'(72("-"))')
  write (6,'(A)') '>>>    it      Etot         rmsf'

  do i=1,nsd_max

     ! commit SD step
     atmcrd(:,:) = atmcrd(:,:) + 2.0d-2 * f(:,:)

     ! calculate new gradient
     call gradient ( e2, f(:,:), atmcrd(:,:) )

     ! calculate rms force
     fnrm   = sqrt ( sum (f(:,:)*f(:,:)) / (3*natoms) )

     ! control print
     if ( mod (i,20) == 0 .or. i==1 ) &
          write (6,'(A,I5,X,2F12.6)') '>>> ', i, e2, fnrm

     ! check for convergence
     if ( fnrm < feps ) exit

  end do

  write (6,'(A,I5,2(A,F12.6))') '>>> Convergence reached: it= ', i, &
                                ', Etot= ', e2, ', rmsf= ', fnrm

  write (6,'(72("-"))')
  write (6,*) ''

  ! Write Vsim format
  call write_vsim   ( trim (file2)//'_relaxed', alat(:), atmcrd(:,:) )
  !----------------------------------------------------------

  ! Calculate normal modes
  call get_modes ( eig(:), atmcrd(:,:), 1.0d-4, alat(:), print=.true. )

  ! Set initial temperature in [K]
  t = 1.0d0
  
  ! Calculate free energy for different temperatures
  do i=1,nt
     call get_fvib ( fvib2(i), t, eig(:) )
     print *, t, fvib2(i)
     t = t + 100.d0
  end do


  !-------------------------------------------------------------------------------
  ! Postprocessing
  !-------------------------------------------------------------------------------
  write (6,*) ''  
  write (6,'(72("-"))')
  write (6,'(A)') '>>> Formation free energy [eV]'
  write (6,'(72("-"))')
  write (6,'(A)') ' Temperature    FFCD      Vacancy'

  t = 1.0d0
  do i=1,nt
     write (6,'(3f12.6)')  t, &
          ! formation free energy of a Si ffcd
          ( e1 + fvib1(i) ) - ( e0 + fvib0(i) ) , &
          ! formation free energy of Si vacancy
           ( e2 + fvib2(i) ) - float (natoms) * ( e0 + fvib0(i) ) / n0  
     t = t + 100.d0
  end do

  !-------------------------------------------------------------------------------
  ! Finalize 
  !-------------------------------------------------------------------------------
  call cpu_time ( t2 )

  write (6, '(72("-"))')
  write (*, *) '               CPU Time= ', t2 - t1, 'sec.'
  write (6, '(72("-"))')

  if ( allocated ( atmcrd ) ) deallocate ( atmcrd)
  if ( allocated ( hess   ) ) deallocate ( hess  )

  deallocate ( f, eig )



!===============================================================================
contains
!===============================================================================


  !----------------------------------------------
  subroutine gradient ( energy, ff, coords, print )
  !----------------------------------------------
  !
  ! This is a simple driver routine to calculate energy and forces
  ! for silicon. 
  !
    implicit none

    real(8), intent(out) :: energy
    real(8), intent(out) :: ff(:,:)
    real(8), intent(in)  :: coords(:,:)
    logical, optional    :: print

    real(8)              :: coord, ener_var, coord_var
    logical              :: qprint
    !-------------------------------------------------

    qprint = .false.
    if ( present (print) ) qprint = print

    call bazant ( natoms, alat, coords, ff, energy, coord, ener_var, coord_var, count )
  
    if ( qprint ) then
       write (6,'(72("-"))')
       write (6,'(A,f15.6)') '          total energy [eV]= ', energy
       write (6,'(A,f15.6)') '           rms force [eV/A]= ', &
                             sqrt ( sum ( ff(:,:) * ff(:,:) ) / (3*natoms) )
       write (6,'(A,f15.6)') 'variance of the energy/atom= ', ener_var
       write (6,'(A,f15.6)') 'average coordination number= ', coord 
       write (6,'(72("-"))')
       write (6,*) ''
    end if

  end subroutine gradient


  !----------------------------------------------
  subroutine hessian ( coords, delta )
  !----------------------------------------------

    implicit none
 
    ! array arguments
    real(8), intent(in) :: coords(:,:) ! atomic coordinates

    ! scalar arguments.
    real(8), intent(in) :: delta       ! step size in [A]
    
    ! local scalars.
    integer             :: i, iatom, ierr, index, j
    real(8)             :: ener, temp
    
    ! local arrays.
    real(8), allocatable, dimension(:,:) :: gm1, gm2, gp2
    !-----------------------------------------------------------

   if ( allocated ( hess ) ) deallocate ( hess )
   allocate ( hess(1:3*natoms,1:3*natoms), stat=ierr )
   if (ierr /= 0 ) stop '** hessian: allocation failed...'

   hess(:,:) = 0.0d0

   allocate ( gm1(1:3,1:natoms), gm2(1:3,1:natoms), gp2(1:3,1:natoms), stat=ierr )
   if (ierr /= 0 ) stop '** hessian: allocation failed...'

   do iatom = 1,natoms
      do i = 1,3         
         temp = atmcrd(i,iatom) ! save initial coordinates

         !---  get function at -2h
         atmcrd(i,iatom) = temp - 2.d0 * delta
         call gradient ( ener, f, coords )
         gm2 = -f

         !---  get function at +2h
         atmcrd(i,iatom) = temp + 2.d0 * delta
         call gradient ( ener, f, coords )
         gp2 = -f

         !--- function at -h
         atmcrd(i,iatom) = temp - delta
         call gradient ( ener, f, coords )
         gm1 = -f

         !--- function at +h
         atmcrd(i,iatom) = temp + delta
         call gradient ( ener, f, coords )

         atmcrd(i,iatom) = temp ! restore coordinates

         index = 3 * ( iatom - 1 ) + i

         hess(:,index) = pack ( -8.d0 * (f + gm1) - gp2 + gm2, .true. )
      end do
   end do

   ! [eV / Angstr.^2]
   hess(:,:) = hess(:,:) / ( 12.d0 * delta )

   ! Mass-weighted Hessian matrix
   hess(:,:) = ev_to_kj * hess(:,:) / m_si  ! complies with DYNAMO f90 library

   deallocate ( gm1, gm2, gp2 )

 end subroutine hessian
 
 
 !---------------------------------------------------------
 subroutine get_modes ( modes, coords, delta, alat, print )
 !---------------------------------------------------------
   
   implicit none
   
   ! Arguments
   real(8), intent(out) :: modes(:)
   real(8), intent(in)  :: coords(:,:)
   real(8), intent(in)  :: delta
   real(8), intent(in)  :: alat(3)
   logical, optional    :: print
   
   ! Local variables
   integer              :: i, info
   real(8)              :: e0, fac
   real(8), allocatable :: work(:)
   logical              :: qprint
   !----------------------------------------------------------
   
   qprint = .false.
   if ( present (print) ) qprint = print
   
   ! Get hessian
   call hessian ( coords, delta )
   
   ! Get number of degrees of freedom: 3N
   n = size ( atmcrd )
   
   ! Allocate workspace
   allocate ( work(1:3*n-1), stat=ierr )
   if (ierr /= 0) stop '** get_modes: allocation failed...'
   
   ! Checksums
   if ( qprint ) write (6,'(2(a,e16.9))') '>>> Checksums: ', &
        sum ( abs ( hess(:,:) ) ), ', ', sum ( hess(:,:) - transpose ( hess(:,:) ) )

   ! Diagonalize symmetric matrix using LAPACK's dsyev() driver routine
   call dsyev ('N', 'U', n, hess, n, modes,  work, 3*n-1, info)
   if ( info /= 0 ) stop '** get_modes: dsyev failed... '

    ! normal mode frequencies in [1/cm]
   modes = to_wavenumbers * sign ( sqrt ( abs ( modes(:) ) ), modes(:) )
   
   ! Optional printout of the eigenfrequencies
   if ( qprint ) then
      write (6, "(/21('-'),a,21('-'))" ) " Harmonic Frequencies [1/cm] "
      write (6, '(5(f10.3,x))'), modes(1:30)
      write (6, '(5x,32("."))')
      write (6, '(5(f10.3,x))'), modes(n-29:n)
      write (6, '(72("-"))')
   end if

   ! Deallocate workspace
   deallocate ( work )

 end subroutine get_modes


 !---------------------------------------------- 
 subroutine get_fvib ( fvib, t, modes )
 !----------------------------------------------

   implicit none

   real(8), intent(out) :: fvib      ! vibrational free energy
   real(8), intent(in)  :: t         ! temperature
   real(8), intent(in)  :: modes(:)  ! normal modes

   integer              :: i
   real(8)              :: fac
   !----------------------------------------------------------

   ! initialize fvib
   fvib = 0.d0

   ! Sum over the 3N-3 non-zero modes
   do i=4,size ( modes )
      ! convert modes from [1/cm] --> [eV], and T from [K] --> [eV]
      fac =  - ( modes(i) * invcm_to_ev ) / ( t * k_to_ev )
      fvib = fvib - log ( exp ( 0.5d0 * fac ) / ( 1.d0 - exp ( fac ) ) )
   enddo

   ! Convert the vibrational free energy in eV
   fvib = t * k_to_ev * fvib
   
 end subroutine get_fvib
 

 !----------------------------------------------
 subroutine xyz_read ( file, print )
 !----------------------------------------------

   implicit none
   
   ! Arguments
   character(len=*), intent(in)  :: file
   logical, optional             :: print 
   
   ! Local variables
   integer       :: i, iu = 5
   character( 4) :: ext='.xyz' ! default input file extension
   character( 2) :: atmlab     ! atomic lable
   character(62) :: dummy
   logical       :: qprint
   !-----------------------------------------------------------
   
   qprint = .false.
   if ( present ( print ) ) qprint = print
   
   write (6,'(2A)') '>>> Reading data from ', file//ext
   open (iu, file=file//ext, status='old')
   
   read (iu,*) natoms
   write(6,'(A,I5)') '>>> natoms= ', natoms
  
   if ( allocated ( atmcrd ) ) deallocate ( atmcrd )
   allocate ( atmcrd(1:3,1:natoms) )
   atmcrd(:,:) = 0.d0
   
   read (iu,*) dummy
   
   ! Read in coordinate record
   do i=1,natoms
      read (iu, *) atmlab, atmcrd(1:3,i)
   end do

   ! Optional printout
   if ( qprint ) then
      do i=1,natoms
         write (6,'(A,A2,X,3F12.6)'), '>>> ', atmlab, atmcrd(:,i)
      end do
   end if
   
   close(iu)

 end subroutine xyz_read


  !-----------------------------------------------
  subroutine write_vsim ( vsim_file, box, coords )
  !-----------------------------------------------
  
    implicit none
    
    ! subroutine arguments
    character(len=*), intent(in) :: vsim_file
    real(8), intent(in)          :: coords(:,:), box(3)
    
    ! local vars & parameters
    integer                      :: i, nat
    real(8)                      :: bbox(3)
    real(8), allocatable         :: x(:,:)
    
    ! vsim format
    character(11) :: vsim_form = '(3f15.9,A4)'
    character(6)  ::        ext= '.ascii'
    character(60) ::     title = '! Realxed geometry'
    !-----------------------------------------------

    nat = size (coords, 2)

    allocate ( x(3,nat) )

    do i=1,3
       x(i,:)  = coords(i,:) - minval ( coords(i,:) )
       bbox(i) = maxval ( x(i,:) )
    end do
    
    open (unit=9, file=trim ( vsim_file )//ext, status='unknown')
    !------------------------------------------------
    ! writing generic format suitable for the vsim software
    !------------------------------------------------
    write (9, '(a)') trim ( title )
    write (9, '(3f15.9)') bbox(1), 0.d0, bbox(2)
    write (9, '(3f15.9)') 0.d0,    0.d0, bbox(3)
    
    do i=1,nat
       write (9, vsim_form) x(1:3,i), ' Si'
    enddo
    !------------------------------------------------
    close (9)
    deallocate ( x )
    return
  end subroutine write_vsim

  
end program si_point_def
