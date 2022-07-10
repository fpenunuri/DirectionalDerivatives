  !Example of Sec. 5.4 (Dual number method)
  
  include "../ModulesDual/dualz1_mod.f90"
  include "../ModulesDual/dualz2_mod.f90"
  include "../ModulesDual/dualz3_mod.f90"
  include "../ModulesDual/dualz4_mod.f90"

  include "../ModulesDual/DMLF_dual_mod.f90"

  !function test module
  module ftestd_mod
    use dualz4_mod
    implicit none
    private

    public :: F5
  contains

    function F5(r) result(fr)
      type(dualz4), intent(in), dimension(:) :: r
      type(dualz4) :: fr
      type(dualz4) :: x, y, z, w, u

      x=r(1); y=r(2); z=r(3); w=r(4); u=r(5)

      if(size(r) /= 5) stop 'vector argument dimension is not 5'
      fr = cos(x*y/u)*z/w + 3*sin(x*u)*sin(y/u)*log(x/u*y*z/w)
    end function F5

  end module ftestd_mod

  !main program
  program main
    !use dualz4_mod
    use DMLF_dual_mod
    use ftestd_mod
    implicit none

    real(8),dimension(5) :: r5v
    real(8) :: d1,d2,d3,d4

    r5v=[1.1d0,2.2d0,3.3d0,4.4d0,5.5d0]

    d1 = df1(F5,1,r5v)
    d2 = df2(F5,[3,3],r5v)
    d3 = df3(F5,[5,4,2],r5v)
    d4 = df4(F5,[5,3,4,1],r5v)
    print*,d1
    print*,d2
    print*,d3
    print*,d4
  end program main
  
