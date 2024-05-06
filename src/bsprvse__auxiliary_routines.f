! ================================================================================================================================ !
module bsprvse__auxiliary_routines

  use bsprvse__constants, only: wp

  implicit none

  private

  public :: spline
  public :: spl

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine SPLINE(N, X, Y, CM)
    use bsprvse__constants, only: zero, three, six
    integer, intent(in) :: N
    real(wp), intent(in) :: X(N)
    real(wp), intent(in) :: Y(N)
    real(wp), intent(out) :: CM(N)
    integer :: i
    real(wp) :: A
    real(wp) :: C
    real(wp) :: ALPHA(N)
    real(wp) :: BETA(N)
    real(wp) :: GAMM(N)
    real(wp) :: B(N)
    CM(1) = zero
    CM(N) = zero
    do I =  3,N
     A =  X(I - 1) - X(I - 2)
     C =  X(I) - X(I - 1)
     ALPHA(I - 2) =  (A + C) / three
     BETA(I - 2) =  C / six
     GAMM(I - 2) =  BETA(I - 2)
     B(I - 2) =  (Y(I) - Y(I - 1)) / C - (Y(I - 1) - Y(I - 2)) / A
    enddo
    call TRIDIA(ALPHA, BETA, GAMM, B, CM(2:N), N - 2)
  end subroutine SPLINE

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine TRIDIA(ALPHA, BETA, GAMM, B, X, N)
    integer, intent(in) :: N
    real(wp), intent(inout) :: ALPHA(:)
    real(wp), intent(in) :: BETA(:)
    real(wp), intent(in) :: GAMM(:)
    real(wp), intent(inout) :: B(:)
    real(wp), intent(inout) :: X(:)

    integer :: i, j
    real(wp) :: RAP

    do i = 2, N
     RAP = BETA(i - 1) / ALPHA(i - 1)
     ALPHA(i) = ALPHA(i) - RAP * GAMM(i - 1)
     B(i) = B(i) - RAP * B(i - 1)
    enddo

    X(N) = B(N) / ALPHA(N)

    do J = 2,N
      i = N - J + 1
     X(i) = (B(i) - GAMM(i) * X(i + 1)) / ALPHA(i)
    enddo

  end subroutine TRIDIA

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function SPL(N, X, Y, M, T) result(res)
    integer, intent(in) :: N
    real(wp), intent(in) :: X(N)
    real(wp), intent(in) :: Y(N)
    real(wp), intent(in) :: M(N)
    real(wp), intent(in) :: T
    real(wp) :: res

    integer :: i, k
    real(wp) :: G, E, F

    IF(T.LE.X(1)) GO TO 30
    IF(T.GE.X(N)) GO TO 40
    K=2
    10   IF(T.LE.X(K)) GO TO 20
    K=K+1
    GO TO 10
    20   E=X(K)-X(K-1)
    F=X(K)-T
    G=T-X(K-1)
    res = (M(K-1)*F*F*F+M(K)*G*G*G+(6.*Y(K)-M(K)*E*E)*G+(6.*Y(K-1)- M(K-1)*E*E)*F)/(6.*E)
    RETURN
    30   E=X(2)-X(1)
    res = ((Y(2)-Y(1))/E-M(2)*E/6.)*(T-X(1))+Y(1)
    RETURN
    40   E=X(N)-X(N-1)
    res = ((Y(N)-Y(N-1))/E+M(N-1)*E/6.)*(T-X(N))+Y(N)

  end function SPL

! ================================================================================================================================ !
end module bsprvse__auxiliary_routines
! ================================================================================================================================ !
