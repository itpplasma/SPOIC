module field_mod

    implicit none

    type :: Field
      real(8) :: epsmn = 0.001d0  ! Perturbation strength
      real(8) :: Ath, Aph
      real(8), dimension(3) :: dAth, dAph
      ! second derivatives: drdr, drdth, drdph, dthdth, dthdph, dphdph
      real(8), dimension(6) :: d2Ath, d2Aph
    end type Field

    contains

    ! for testing -> circular tokamak
    subroutine eval_field(f, r, th, ph, mode_secders)
    !
      type(Field), intent(inout) :: f
      real(8), intent(in) :: r, th, ph
      integer, intent(in) :: mode_secders

      real(8) :: B0, iota0, a, R0, cth, sth

       B0 = 1.0    ! magnetic field modulus normalization
       iota0 = 1.0 ! constant part of rotational transform
       a = 0.5     ! (equivalent) minor radius
       R0 = 1.0    ! (equivalent) major radius

       cth = cos(th)
       sth = sin(th)

       f%Ath = B0*(r**2/2d0 - r**3/(3d0*R0)*cth)
       f%Aph = -B0*iota0*(r**2/2d0-r**4/(4d0*a**2))   + f%epsmn*cos(2d0*th + ph)

       f%dAth(1) = B0*(r - r**2/R0*cth)
       f%dAth(2) = B0*r**3*sth/(3.0*R0)
       f%dAth(3) = 0d0

       f%dAph(1) = -B0*iota0*(r-r**3/a**2)
       f%dAph(2) = 0d0    - f%epsmn*2d0*sin(2d0*th + ph)
       f%dAph(3) = 0d0    - f%epsmn*sin(2d0*th + ph)

       if (mode_secders <= 0) return

       f%d2Ath(1) = B0*(1d0 - 2d0*r/R0*cth)
       f%d2Ath(4) = B0*r**3*cth/(3d0*R0)
       f%d2Ath(3) = 0d0
       f%d2Ath(2) = B0*r**2/R0*sth
       f%d2Ath(5) = 0d0
       f%d2Ath(6) = 0d0

       f%d2Aph(1) = -B0*iota0*(1d0-3d0*r**2/a**2)
       f%d2Aph(4) = 0d0
       f%d2Aph(3) = 0d0
       f%d2Aph(2) = 0d0
       f%d2Aph(5) = 0d0
       f%d2Aph(6) = 0d0

    end subroutine eval_field

end module field_mod
