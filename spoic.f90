  program test_read_spec
    use read_spec
    implicit none

    type(SpecOutput) :: s
    character(*), parameter :: filename = 'spec_out.h5'

    double precision :: lvol
    logical :: Lcoordinatesingularity

    integer, parameter :: ivol = 1
    double precision, allocatable :: cheby(:,:), zernike(:,:,:)

    double precision, parameter :: zero=0.0d0, one=1.0d0, half=0.5d0
    double precision, parameter :: machprec = 1.11d-16, small = 10000*machprec
    double precision st(2), zeta, lss, teta, sbar


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
        call get_zernike(sbar, Lrad(ivol), Mpol, zernike(:,:,0:1))
      else
        call get_cheby(lss, Lrad(ivol), cheby(0:Lrad(ivol),0:1))
      end if

9999  deallocate(cheby, zernike)

    end associate

    call freeSpec(s)
  end program test_read_spec
