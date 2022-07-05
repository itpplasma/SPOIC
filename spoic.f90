  program test_read_spec
    use read_spec
    implicit none

    type(SpecOutput) :: s
    character(*), parameter :: filename = 'test.sp.h5'

    logical :: Lcoordinatesingularity

    integer, parameter :: ivol = 1
    real(8), allocatable :: cheby(:,:), zernike(:,:,:), TT(:,:)

    real(8), parameter :: zero=0.0d0, one=1.0d0, half=0.5d0
    real(8), parameter :: machprec=1.11d-16, small = 10000*machprec
    real(8) :: st(2), zeta, lss, teta, sbar, arg, carg, sarg

    integer :: lvol, ii, ll, mi, ni

    real(8) :: Az, At

    zeta = 0.3
    st(1) = 0.3
    st(2) = 0.5

    lss = st(1) ; teta = st(2)
    lvol = ivol

    !write(*,*) "reading '",filename,"'..."
    call loadSpec(s, filename)
    !write(*,*) "done"

    !write(*,"(A,F4.2)") "SPEC version: ", s%version
    !write(*,"(A,99I4)") "Lrad:", s%input%physics%Lrad

    associate(Lrad => s%input%physics%Lrad, &
      Mpol => s%input%physics%Mpol, mn => s%output%mn, &
      im => s%output%im, in => s%output%in, &
      Aze => s%vector_potential%Aze, Ate => s%vector_potential%Ate, &
      Azo => s%vector_potential%Azo, Ato => s%vector_potential%Ato)

      allocate(cheby(0:Lrad(ivol),0:1), zernike(0:Lrad(1),0:Mpol,0:1), &
        TT(0:Lrad(ivol),0:1))

      if (lvol < 1) then
        Lcoordinatesingularity=.True.
      else
        Lcoordinatesingularity=.False.
      end if

      lss = st(1) ; teta = st(2)

      if( abs(lss).gt.one ) goto 9999 ! out of domain;

      if( Lcoordinatesingularity ) sbar = max( ( lss + one ) * half, small )

      if (Lcoordinatesingularity) then
        !print *, 'Zernike', sbar, Lrad(ivol), Mpol
        call get_zernike(sbar, Lrad(ivol), Mpol, zernike(:,:,0:1))
        !print *, zernike(:,:,0:1)
      else
        !print *, 'Cheby'
        call get_cheby(lss, Lrad(ivol), cheby(0:Lrad(ivol),0:1))
        !print *, cheby(0:Lrad(ivol),0:1)
      end if

      Az = 0d0
      At = 0d0

      do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ; arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg) ! shorthand; 20 Apr 13;

        !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        if( Lcoordinatesingularity ) then ! regularization factor depends on mi; 17 Dec 15;

        do ll = 0, Lrad(lvol) ; TT(ll,0:1) = (/        zernike(ll,mi,0),        zernike(ll,mi,1)*half                      /)
        enddo

        else

        do ll = 0, Lrad(lvol) ; TT(ll,0:1) = (/             cheby(ll,0),             cheby(ll,1)                           /)
        enddo

        endif ! end of if( Lcoordinatesingularity ) ; 16 Jan 15;

    !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        do ll = 0, Lrad(lvol) ! loop over Chebyshev summation; 20 Feb 13;
          Az = Az + Aze(lvol,ii) * TT(ll,0) * carg + Azo(lvol,ii) * TT(ll,0) * sarg
          At = At + Ate(lvol,ii) * TT(ll,0) * carg + Ato(lvol,ii) * TT(ll,0) * sarg
        enddo ! end of do ll; 10 Dec 15;

    !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      enddo ! end of do ii = 1, mn;

      !print *, Az, At

      call trace

    end associate

9999 call finalize

  contains

    subroutine trace
      use field_mod, only: Field, eval_field

      type(Field) :: f
      integer, parameter :: nphi = 32, nturn = 20
      integer :: i, j
      real(8) :: z(3), r, Athold, h, fun, dfun

      z(1) = 0.3d0; r = z(1)
      z(2) = 0.1d0
      z(3) = 0.2d0

      h = 2d0*3.14159265358979323846d0/nphi

      call eval_field(f, z(1), z(2), z(3), 1)
      Athold = f%Ath

      do i = 1,nturn*nphi
        do j = 1,10
          call eval_field(f, r, z(2), z(3), 1)
          fun = f%dAth(1)*(f%Ath - Athold) &
              + h*(f%dAph(2)*f%dAth(1) - f%dAph(1)*f%dAth(2))
          dfun = f%d2Ath(1)*(f%Ath - Athold) + f%dAth(1)**2 &
               + h*(f%d2Aph(2)*f%dAth(1) + f%dAph(2)*f%d2Ath(1) &
                  - f%d2Aph(1)*f%dAth(2) - f%dAph(1)*f%d2Ath(2))

          r = r - fun/dfun
        end do
        call eval_field(f, r, z(2), z(3), 1)
        Athold = f%Ath
        z(1) = r; z(2) = z(2) + h*f%dAph(1)/f%dAth(1); z(3) = z(3) + h
        print *, z
      end do
    end subroutine trace

    subroutine initialize
      ! TODO: Allocate zernike, cheby and TT based on maximum shape
    end subroutine initialize

    subroutine finalize
      deallocate(cheby, zernike)
      call freeSpec(s)
    end subroutine finalize

  end program test_read_spec
