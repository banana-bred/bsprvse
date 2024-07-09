! ================================================================================================================================ !
module bsprvse__utilities

  implicit none

  private

  public :: i2char
  public :: ndigits
  public :: die
  public :: swapvals
  public :: swapcols
  public :: append

  interface die
    module procedure :: die_1
    module procedure :: die_2
  end interface die

  interface swapvals
    !! Swap two values of an array
    module procedure :: swapvals_i
    module procedure :: swapvals_r
    module procedure :: swapvals_c
  end interface swapvals

  interface swapcols
    !! Swap two columnss of an array
    module procedure :: swapcols_i
    module procedure :: swapcols_r
    module procedure :: swapcols_c
  end interface swapcols

  interface append
    !! Append an element to an array. If not allocated, allocate the array to contain that element
    module procedure append_i
    module procedure append_r
    module procedure append_c
  end interface append


! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function ndigits(n) result(num)
    !! Returns number of characters an integer will occupy

    use bsprvse__constants, only: one

    integer, intent(in) :: n
    integer :: num

    num = 1

    if(n .eq. 0) return

    num = floor(log10(abs(n) * one)) + 1

    ! -- account for minus sign
    if(n .lt. 1) num = num + 1

  end function ndigits

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  recursive function i2char(i) result(output)
    !! Convert the integer i to a character(n). This function determines the value to use such that the output character array
    !! 'output' is exactly the size needed to fit the integer i, e.g.,
    !!   int2char0(1) -> '1'
    !!   int2char0(19) -> '19'

    integer, intent(in) :: i
    character(:), allocatable :: output

    character(:), allocatable :: errmsg
    character(:), allocatable :: frmt

    integer :: nn

    nn = ndigits(i)

    ! -- allocate characters so that they're wide enough to write to
    allocate(character(nn) :: output)

    write(output, '(I0)') i

  end function i2char

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine die_1(message)
    !! Stop program execution with a message
    use iso_fortran_env, only: stderr => error_unit
    character(*), intent(in), optional :: message
    write(stderr,*)
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@                                                           @@@@")')
    write(stderr,'("@@@@                          ERROR                            @@@@")')
    write(stderr,'("@@@@                                                           @@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,*)
    if(.not.present(message)) error stop
    write(stderr,'("STOP",1X,"::",1X,A)') message
    write(stderr,*)
    error stop
  end subroutine die_1
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine die_2(message1,message2)
    !! Stop program execution with two messages
    use iso_fortran_env, only: stderr => error_unit
    character(*), intent(in) :: message1
    character(*), intent(in) :: message2
    write(stderr,*)
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@                                                           @@@@")')
    write(stderr,'("@@@@                          ERROR                            @@@@")')
    write(stderr,'("@@@@                                                           @@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,*)
    write(stderr,'("STOP",1X,"::",1X,A,/,A)') message1, message2
    write(stderr,*)
    error stop
  end subroutine die_2

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine swapvals_i(arr, i, j)
    !! Swap the values indexd by `i` and `j` of the array `arr`
    use bsprvse__constants, only: wp
    integer, intent(inout) :: arr(:)
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer :: tmp
    tmp = arr(i)
    arr(i) = arr(j)
    arr(j) = tmp
  end subroutine swapvals_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine swapvals_r(arr, i, j)
    !! Swap the values indexd by `i` and `j` of the array `arr`
    use bsprvse__constants, only: wp
    real(wp), intent(inout) :: arr(:)
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(wp) :: tmp
    tmp = arr(i)
    arr(i) = arr(j)
    arr(j) = tmp
  end subroutine swapvals_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine swapvals_c(arr, i, j)
    !! Swap the values indexd by `i` and `j` of the array `arr`
    use bsprvse__constants, only: wp
    complex(wp), intent(inout) :: arr(:)
    integer, intent(in) :: i
    integer, intent(in) :: j
    complex(wp) :: tmp
    tmp = arr(i)
    arr(i) = arr(j)
    arr(j) = tmp
  end subroutine swapvals_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine swapcols_i(arr, i, j)
    !! Swap the values indexd by `i` and `j` of the array `arr`
    use bsprvse__constants, only: wp
    integer, intent(inout) :: arr(:,:)
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, allocatable :: tmp(:)
    tmp = arr(:, i)
    arr(:, i) = arr(:, j)
    arr(:, j) = tmp
  end subroutine swapcols_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine swapcols_r(arr, i, j)
    !! Swap the values indexd by `i` and `j` of the array `arr`
    use bsprvse__constants, only: wp
    real(wp), intent(inout) :: arr(:,:)
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(wp), allocatable :: tmp(:)
    tmp = arr(:, i)
    arr(:, i) = arr(:, j)
    arr(:, j) = tmp
  end subroutine swapcols_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine swapcols_c(arr, i, j)
    !! Swap the values indexd by `i` and `j` of the array `arr`
    use bsprvse__constants, only: wp
    complex(wp), intent(inout) :: arr(:,:)
    integer, intent(in) :: i
    integer, intent(in) :: j
    complex(wp), allocatable :: tmp(:)
    tmp = arr(:, i)
    arr(:, i) = arr(:, j)
    arr(:, j) = tmp
  end subroutine swapcols_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! APPEND INTERFACE
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine append_i(arr, new)
    !! Append element "new" to array "arr"

    implicit none

    integer, intent(in)                 :: new
    integer, intent(inout), allocatable :: arr(:)

    select case(allocated(arr))
    case(.true.)
      arr = [arr, new]
    case(.false.)
      arr = [new]
    end select

  end subroutine append_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine append_r(arr, new)
    !! Append element "new" to array "arr"

    use bsprvse__constants, only: wp

    implicit none

    real(wp), intent(in)                 :: new
    real(wp), intent(inout), allocatable :: arr(:)

    select case(allocated(arr))
    case(.true.)
      arr = [arr, new]
    case(.false.)
      arr = [new]
    end select

  end subroutine append_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine append_c(arr, new)
    !! Append element "new" to array "arr"

    use bsprvse__constants, only: wp

    implicit none

    complex(wp), intent(in)                 :: new
    complex(wp), intent(inout), allocatable :: arr(:)

    select case(allocated(arr))
    case(.true.)
      arr = [arr, new]
    case(.false.)
      arr = [new]
    end select

  end subroutine append_c

! ================================================================================================================================ !
end module bsprvse__utilities
! ================================================================================================================================ !
