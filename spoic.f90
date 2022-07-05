  program test_read_spec
    use read_spec
    implicit none

    type(SpecOutput) :: s
    character(*), parameter :: filename = 'test.sp.h5'

    real(8) :: lvol
    logical :: Lcoordinatesingularity

    integer, parameter :: ivol = 1
    real(8), allocatable :: cheby(:,:), zernike(:,:,:)

    real(8), parameter :: zero=0.0d0, one=1.0d0, half=0.5d0
    real(8), parameter :: machprec=1.11d-16, small = 10000*machprec
    real(8) :: st(2), zeta, lss, teta, sbar

    double precision :: x_dp
    real(8) :: x_real
    real :: x_real4

    print *, kind(x_dp), kind(x_real), kind(x_real4), kind(1.0e0)

    zeta = 0.3
    st(1) = 0.3
    st(2) = 0.5

    lss = st(1) ; teta = st(2)

    write(*,*) "reading '",filename,"'..."
    call loadSpec(s, filename)
    write(*,*) "done"

    write(*,"(A,F4.2)") "SPEC version: ", s%version
    write(*,"(A,99I4)") "Lrad:", s%input%physics%Lrad

    associate(Lrad => s%input%physics%Lrad, &
      Mpol => s%input%physics%Mpol)

      allocate(cheby(0:Lrad(ivol),0:1), zernike(0:Lrad(1),0:Mpol,0:1))

      if (lvol < 1) then
        Lcoordinatesingularity=.True.
      else
        Lcoordinatesingularity=.False.
      end if

      lss = st(1) ; teta = st(2)

      if( abs(lss).gt.one ) goto 9999 ! out of domain;

      if( Lcoordinatesingularity ) sbar = max( ( lss + one ) * half, small )

      if (Lcoordinatesingularity) then
        print *, 'Zernike', sbar, Lrad(ivol), Mpol
        call get_zernike(sbar, Lrad(ivol), Mpol, zernike(:,:,0:1))
        print *, zernike(:,:,0:1)
      else
        print *, 'Cheby'
        call get_cheby(lss, Lrad(ivol), cheby(0:Lrad(ivol),0:1))
        print *, cheby(0:Lrad(ivol),0:1)
      end if


9999  deallocate(cheby, zernike)

    end associate

    call freeSpec(s)
  end program test_read_spec
