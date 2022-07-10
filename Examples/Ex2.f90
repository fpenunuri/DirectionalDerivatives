  !Example of Sec. 5.2 (Dual number method)
  
  include "../ModulesDual/dualz1_mod.f90"
  include "../ModulesDual/dualz2_mod.f90"
  include "../ModulesDual/dualz3_mod.f90"
  include "../ModulesDual/dualz4_mod.f90"

  include "../ModulesDual/DMLF_dual_mod.f90"

  !function test module
  module ftest_mod
    implicit none
    private

    public :: f1appx

  contains
    function f1appx(q) result(fr)
      use dualz4_mod
      type(dualz4), intent(in), dimension(:) :: q
      type(dualz4), allocatable, dimension(:) :: fr
      type(dualz4) :: q1, q2, q3

      q1 = q(1); q2 = q(2); q3 = q(3)
      allocate (fr(2))
      fr = [sin(q1 * q2**2 * q3**3), cos(q1 * q2**2 * q3**3)]
    end function f1appx
  end module ftest_mod

  !main program
  program Eappx
    use DMLF_dual_mod
    use ftest_mod
    implicit none

    real(8), dimension(3), parameter :: v = [1,2,3]
    real(8), dimension(3), parameter :: q = [1.1d0,2.2d0,3.3d0]
    complex(8), dimension(3), parameter :: qc = q
    real(8), dimension(2) :: JFqv
    integer :: k

    JFqv = d1mlf(f1appx,v,q,2)

    do k = 1, 2
       print*,JFqv(k)
    end do
  end program Eappx
