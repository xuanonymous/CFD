  MODULE types_vars

! This module defines the KIND types of all the variables used in the code: 
! I4B, I2B and I1B for integer variables, SP and DP for real variables (and
! SPC and DPC for corresponding complex cases), and LGT for the default 
! logical type. This follows the convention used the Numerical Recipes for 
! Fortran 90 types module 'nrtype', pp. 1361
!
! Symbolic names for kind types of 4-, 2- and 1-byte integers:   

  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)

! Symbolic names for kind types of single- and double-precison reals

  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  
! Symbolic names for kind types of single- and double-precison complex

  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))

! Symbolic name for kind type of default logical

  INTEGER, PARAMETER :: LOGIC = KIND(.true.)

! Frequently used mathematical constants (with precision to spare)

  REAL(SP), PARAMETER :: one   = 1.0_sp
  REAL(SP), PARAMETER :: zero  = 0.0_sp
  REAL(SP), PARAMETER :: pi    = 3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: pio2  = 1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: twopi = 6.283185307179586476925286766559005768394_sp

  REAL(DP), PARAMETER :: one_d   = 1.0_dp
  REAL(DP), PARAMETER :: zero_d  = 0.0_dp
  REAL(DP), PARAMETER :: pi_d    = 3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: pio2_d  = 1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: twopi_d = 6.283185307179586476925286766559005768394_dp

  END MODULE types_vars
  
!-------------------------------------------------------

  PROGRAM precise

  USE types_vars
  IMPLICIT NONE


  REAL :: a, b
  REAL(SP) :: c, d
  REAL(DP) :: e, f

  a = 4.5678
  b = 9.1234
  c = 5.3214
  d = 7.8890
  e = 1.2345

  print*, a, b, c, d, e
  print*, kind(a), kind(b), kind(c), kind(d), kind(e)

  END PROGRAM precise
