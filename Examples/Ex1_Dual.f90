  !Example of Sec. 5.1 (Dual number method)

  include "../ModulesDual/dualz1_mod.f90"
  include "../ModulesDual/dualz2_mod.f90"
  include "../ModulesDual/dualz3_mod.f90"
  include "../ModulesDual/dualz4_mod.f90"

  include "../ModulesDual/DMLF_dual_mod.f90"

  !function test module
  module f_mod
    implicit none
    private

    public ::  InvCosWF

  contains
    function  InvCosWF(q) result(fr)
      use dualz4_mod
      type(dualz4), intent(in), dimension(:) :: q
      type(dualz4) :: fr
      type(dualz4), dimension(size(q)-1) :: xcut1, xcut2, num, disc

      xcut1 = q(1:size(q)-1)
      xcut2 = q(2:size(q))
      num = -(xcut1**2 + xcut2**2 + 0.5*xcut1*xcut2)
      disc = xcut1**2 + xcut2**2 + 0.5*xcut1*xcut2
      fr = -sum(exp(num/8.0)*cos(4*sqrt(disc)))
    end function InvCosWF
  end module f_mod

  !main program
  program Ea1
    use DMLF_dual_mod
    use f_mod
    implicit none

    integer, parameter :: m = 20
    real(8), dimension(m) :: x, y, z, w, q
    real(8) :: derdir4
    integer :: k
    real(8) :: t1,t2,t

    do k=1,m
       x(k)=k
    end do

    y = sin(x); z = cos(x); w = sqrt(x); q = log(x)

    call cpu_time(t1)
    derdir4 = d4mlf(InvCosWF,x,y,z,w,q)
    call cpu_time(t2)
    t = t2-t1

    write(*,*) "m=", m
    write(*,*) "d4=", derdir4
    write(*,*) "t(s)=",t
  end program Ea1
