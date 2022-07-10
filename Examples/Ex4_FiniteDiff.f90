  !Example of Sec. 5.4 (FMFD method)

  include '../ModulesFinDiff/finite_diff4_mod.f90'

  !function test module
  module ftest_mod
    implicit none
    private

    public :: F5_4, F5_8, F5_16
  contains

    function F5_4(r) result(fr)
      real(4), intent(in), dimension(:) :: r
      real(4) :: fr
      real(4) :: x, y, z, w, u

      x=r(1); y=r(2); z=r(3); w=r(4); u=r(5)

      if(size(r) /= 5) stop 'vector argument dimension is not 5'
      fr = cos(x*y/u)*z/w + 3*sin(x*u)*sin(y/u)*log(x*y*z/(u*w))
    end function F5_4

    function F5_8(r) result(fr)
      real(8), intent(in), dimension(:) :: r
      real(8) :: fr
      real(8) :: x, y, z, w, u

      x=r(1); y=r(2); z=r(3); w=r(4); u=r(5)

      if(size(r) /= 5) stop 'vector argument dimension is not 5'
      fr = cos(x*y/u)*z/w + 3*sin(x*u)*sin(y/u)*log(x/u*y*z/w)
    end function F5_8

    function F5_16(r) result(fr)
      real(16), intent(in), dimension(:) :: r
      real(16) :: fr
      real(16) :: x, y, z, w, u

      x=r(1); y=r(2); z=r(3); w=r(4); u=r(5)

      if(size(r) /= 5) stop 'vector argument dimension is not 5'
      fr = cos(x*y/u)*z/w + 3*sin(x*u)*sin(y/u)*log(x/u*y*z/w)
    end function F5_16
  end module ftest_mod

  !main program
  program main
    use finite_diff4_mod
    use ftest_mod
    implicit none
    real(16) :: h16, d1_16, d2_16, d3_16, d4_16
    real(8)  :: h8, d1_8, d2_8, d3_8, d4_8
    real(4)  :: h4, d1_4, d2_4, d3_4, d4_4
    real(16), dimension(5) :: r0v_16
    real(8),  dimension(5) :: r0v_8
    real(4),  dimension(5) :: r0v_4

    !the evaluating point [1.1,2.2,3.3,4.4,5.5]
    r0v_4  = [1.1,2.2,3.3,4.4,5.5]
    r0v_8  = [1.1_8,2.2_8,3.3_8,4.4_8,5.5_8]
    r0v_16 = [1.1_16,2.2_16,3.3_16,4.4_16,5.5_16]

    h4  = 1.0/1e5    
    h8  = 1.0_8/1e5
    h16 = 1.0_16/1e5

    d1_4  = df1_diff(F5_4,[1],r0v_4,h4,8)    
    d1_8  = df1_diff(F5_8,[1],r0v_8,h8,8) 
    d1_16 = df1_diff(F5_16,[1],r0v_16,h16,8) 

    print*,d1_4
    print*,d1_8
    print*,d1_16   
    print*,'----'
    d2_16 = df2_diff(F5_16,[3,3],r0v_16,h16,8)
    d2_8  = df2_diff(F5_8,[3,3],r0v_8,h8,8)   
    d2_4  = df2_diff(F5_4,[3,3],r0v_4,h4,8)   
    print*,d2_4
    print*,d2_8
    print*,d2_16   
    print*,'----'

    d3_16 = df3_diff(F5_16,[5,4,2],r0v_16,h16,8)
    d3_8  = df3_diff(F5_8,[5,4,2],r0v_8,h8,8)    
    d3_4  = df3_diff(F5_4,[5,4,2],r0v_4,h4,8)   
    print*,d3_4
    print*,d3_8
    print*,d3_16   
    print*,'----'

    d4_16 = df4_diff(F5_16,[5,3,4,1],r0v_16,h16,8) 
    d4_8  = df4_diff(F5_8 ,[5,3,4,1],r0v_8,h8,8)  
    d4_4  = df4_diff(F5_4 ,[5,3,4,1],r0v_4,h4,8)  
    print*,d4_4
    print*,d4_8
    print*,d4_16   
  end program main
