  program test_read_spec
    use read_spec
    implicit none

    type(SpecOutput) :: s
    character(*), parameter :: filename = 'spec_out.h5'

    double precision :: lvol
    logical :: Lcoordinatesingularity

    integer, parameter :: ivol = 1
    double precision :: cheby(0:Lrad(ivol),0:1), zernike(0:Lrad(1),0:Mpol,0:1)

    write(*,*) "reading '",filename,"'..."
    call loadSpec(s, filename)
    write(*,*) "done"

    write(*,"(A,F4.2)") "SPEC version: ", s%version
    write(*,"(A,99I4)") "Lrad:", s%input%physics%Lrad

    associate(Lrad => s%input%physics%Lrad, &
      Mpol => s%input%physics%Mpol)

    if (lvol < 1) then
      Lcoordinatesingularity=.True.
    else
      Lcoordinatesingularity=.False.
    end if

    if (Lcoordinatesingularity) then
      call get_zernike(sbar, Lrad(ivol), Mpol, zernike(:,:,0:1))
    else
      call get_cheby(lss, Lrad(ivol), cheby(0:Lrad(ivol),0:1))
    end if

    end associate

    call freeSpec(s)
  end program test_read_spec
