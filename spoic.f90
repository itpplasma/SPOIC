program spoic

  implicit none



  contains

  subroutine afield(x, A)
    use allglobal, only : ivol
    real(8), intent(in) :: x      ! Position s, theta, zeta
    real(8), intent(out) :: A(2)  ! Covariant A_theta and A_zeta

    real(8) :: TT(0:Lrad(ivol),0:1) ! this is almost identical to cheby; 17 Dec 15;
  end subroutine afield

end program spoic
