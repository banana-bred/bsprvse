! ================================================================================================================================ !
module bsprvse__constants

  use iso_fortran_env, only: real64

  implicit none

  private

  integer,     parameter, public :: wp    = real64
    !! The working precision of real and complex types
  real(wp),    parameter, public :: zero  = 0.0_wp
  real(wp),    parameter, public :: one   = 1.0_wp
  real(wp),    parameter, public :: two   = 2.0_wp
  real(wp),    parameter, public :: three = 3.0_wp
  real(wp),    parameter, public :: five  = 5.0_wp
  real(wp),    parameter, public :: six   = 6.0_wp
  real(wp),    parameter, public :: Ryd   = 219474.6313710_wp
    !! multiplication factor to convert atomic units -> inverse centimeters
  real(wp),    parameter, public :: au2ev = 27.2113834e0_wp
    !! multiplication factor to convert atomic units -> eV
  real(wp),    parameter, public :: au2amu = 5.4858010860603975e-4_wp
    !! multiplication factor to convert atomic units of mass (electron mass = 1) to atomic mass units (Daltons)

  complex(wp), parameter, public :: ci = (zero, one)
    !! $\sqrt{-1}$

  character(13), parameter, public :: numeric = "0123456789.+-"
    !! characters considered "numeric"

  ! -- ASCII constants
  integer, parameter, public :: uppercase_a = ichar('A')
  integer, parameter, public :: uppercase_z = ichar('Z')
  integer, parameter, public :: lowercase_a = ichar('a')
  integer, parameter, public :: lowercase_z = ichar('z')

! ================================================================================================================================ !
end module
! ================================================================================================================================ !
