  !Example of Sec. 5.3 (Dual number method)
  
  include "../ModulesDual/dualz1_mod.f90"
  include "../ModulesDual/dualz2_mod.f90"
  include "../ModulesDual/dualz3_mod.f90"
  include "../ModulesDual/dualz4_mod.f90"

  include "../ModulesDual/DMLF_dual_mod.f90"

  !function test module
  module ftest_mod   
    implicit none
    private

    public :: fE3artic

  contains
    function fE3artic(q) result(fr)
      use dualz4_mod
      type(dualz4), intent(in), dimension(:) :: q
      type(dualz4), allocatable, dimension(:) :: fr
      type(dualz4) :: q1, q2
      integer :: k, dimf, dimq

      dimf = 3
      allocate(fr(dimf))

      dimq = 2
      if(dimq /= size(q)) stop 'dim q is not 2'

      q1 = q(1); q2 = q(2)
      fr(1) = log(q1* q2**2)     
      fr(2) = sqrt(q2)/q1

      fr(3) = sin(q1*q2)
      do k=1,500-1
         fr(3) = sin(fr(3))
      end do
    end function fE3artic
  end module ftest_mod

  !main program
  program E3artic
    use DMLF_dual_mod
    use ftest_mod
    use dualz4_mod
    implicit none

    real(8), dimension(2), parameter :: v1 = [0.1d0,0.2d0], v2 = sin(v1)
    real(8), dimension(2), parameter :: q = [1.1d0,2.2d0]
    complex(8), dimension(2), parameter :: v1c = v1, v2c = v2
    complex(8), dimension(2), parameter :: qc = q + (0,0.1d0)
    type(dualz4) :: e1
    type(dualz4), dimension(3) :: dmlfvals
    complex(8), dimension(3) :: dml3, dml4
    integer :: k
    real(8), dimension(2) ::  q1p, q2p
    real(8), dimension(3) :: accel

    e1 = dualz4(0,1,0,0,0)
    dmlfvals = fE3artic(qc + e1*v1c)

    !derivative MLF
    dml3 =  dmlfvals%f3
    dml4 =  dmlfvals%f4

    print*,"der3MLF=", dml3
    print*,'------'
    print*,"der4MLF=", dml4
    print*,'------'
    

    !acceleration
    q1p = [0.5d0,-2.7d0]
    q2p = [-0.1d0,0.7d0]

    accel = d1mlf(fE3artic,q2p,q,3) + d2mlf(fE3artic,q1p,q1p,q,3)
    print*,'acceleration components'
    do k=1,3
       write(*,'(*(f0.11,2x))')  accel(k)
    end do
  end program E3artic
  
