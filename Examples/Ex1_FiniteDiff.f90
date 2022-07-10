  !Example of Sec. 5.1 (FMFD method)
  include '../ModulesFinDiff/finite_diff4_mod.f90'

  !function test module
  module ftest_mod    
    implicit none
    private

    public ::  InvCosWF
  contains
    function  InvCosWF(q) result(fr)
      use iso_fortran_env,  only: prc => real128
      real(prc), intent(in), dimension(:) :: q
      real(prc) :: fr
      real(prc), dimension(size(q)-1) :: xcut1, xcut2, num, disc

      xcut1 = q(1:size(q)-1)
      xcut2 = q(2:size(q))
      num = -(xcut1**2 + xcut2**2 + 0.5*xcut1*xcut2)
      disc = xcut1**2 + xcut2**2 + 0.5*xcut1*xcut2
      fr = -sum(exp(num/8.0)*cos(4*sqrt(disc)))
    end function InvCosWF
  end module ftest_mod

  !main program
  program Ea1D
    use iso_fortran_env,  only: prc => real128
    use finite_diff4_mod
    use ftest_mod
    implicit none
    integer, parameter :: m = 5, N=4
    real(prc), dimension(m) :: x, y, z, w, q
    real(prc) :: derdir4, Aijkl
    integer :: i,j,k,l
    real :: t1,t2,t
    real(prc), parameter :: h=1d-5

    do k=1,m
       x(k)=k
    end do

    y = sin(x); z = cos(x); w = sqrt(x); q = log(x)

    derdir4=0
    call cpu_time(t1)
    do i=1,m
       do j=1,m
          do k=1,m
             do l=1,m
                Aijkl = df4_diff(InvCosWF,[i,j,k,l],q,h,N)
                derdir4 = derdir4 + Aijkl*x(i)*y(j)*z(k)*w(l)
             end do
          end do
       end do
    end do
    call cpu_time(t2)
    t = t2-t1
    write(*,*) "m=", m
    write(*,*) "d4=", derdir4
    write(*,*) "t(s)=",t    
  end program Ea1D

