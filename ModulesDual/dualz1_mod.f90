!Dual numbers of complex components
!This can be used to compute first order derivatives
!F. Peñuñuri
!UADY, Merida Yucatan Mexico
!2021
!
module dualz1_mod
  implicit none
  private

  !Dual number type
  !Notice we use complex(8) components
  !dualz1 number
  !dualz1 ---> [f0,f1]
  type, public :: dualz1
     complex(8) :: f0, f1    
  end type dualz1

  public :: assignment (=)

  public :: operator(==), operator(/=), operator (+),  operator (-),   &
       operator (*), operator(**), operator (/)

  public :: sin, cos, exp, log, tan, sqrt, abs, asin, acos, atan, cosh,&
       sinh, tanh, atanh, asinh, acosh, atan2, matmul, sum, product,   &
       conjg, dot_product

  public :: f0part_dz1, f1part_dz1
  public :: todualz1 

  !INTERFACES
  !equal assignment
  interface assignment (=)
     module procedure igual1z_dz1   !dualz1 <--- complex8
     module procedure igual1z4_dz1  !dualz1 <--- complex4
     module procedure igual1r_dz1   !dualz1 <--- real8
     module procedure igual1r4_dz1  !dualz1 <--- real4
     module procedure igual1i8_dz1  !dualz1 <--- integer8
     module procedure igual1i4_dz1  !dualz1 <--- integer4
  end interface assignment (=)

  !Logical equal operator 
  interface operator (==) 
     module procedure eq_dz1     !dualz1 == dualz1
  end interface operator (==)

  !Logical not equal operator 
  interface operator (/=) 
     module procedure noteq_dz1  !dualz1 /= dualz1
  end interface operator (/=)

  !Plus operator
  interface operator (+)
     module procedure mas_dz1      !unary 
     module procedure suma_dz1     !dualz1 + dualz1
     module procedure sumadz8_dz1  !dualz1 + complex8
     module procedure sumadz4_dz1  !dualz1 + complex4
     module procedure sumadzr8_dz1 !dualz1 + real8
     module procedure sumadzr4_dz1 !dualz1 + real4
     module procedure sumadzi8_dz1 !dualz1 + integer8
     module procedure sumadzi4_dz1 !dualz1 + integer4          
     module procedure sumaz8d_dz1  !complex8 + dualz1
     module procedure sumaz4d_dz1  !complex4 + dualz1
     module procedure sumar8dz_dz1 !real8 + dualz1
     module procedure sumar4dz_dz1 !real4 + dualz1
     module procedure sumai8dz_dz1 !integer8 + dualz1
     module procedure sumai4dz_dz1 !integer4 + dualz1
  end interface operator (+)

  !Minus operator
  interface operator (-)
     module procedure menos_dz1     !unary
     module procedure resta_dz1     !dualz1 - dualz1
     module procedure restadz8_dz1  !dualz1 - complex8
     module procedure restadz4_dz1  !dualz1 - complex4
     module procedure restadzr8_dz1 !dualz1 - real8
     module procedure restadzr4_dz1 !dualz1 - real4
     module procedure restadzi8_dz1 !dualz1 - integer8
     module procedure restadzi4_dz1 !dualz1 - integer4          
     module procedure restaz8d_dz1  !complex8 - dualz1
     module procedure restaz4d_dz1  !complex4 - dualz1
     module procedure restar8dz_dz1 !real8 - dualz1
     module procedure restar4dz_dz1 !real4 - dualz1
     module procedure restai8dz_dz1 !integer8 - dualz1
     module procedure restai4dz_dz1 !integer4 - dualz1
  end interface operator (-)

  !Times operator
  interface operator (*)
     module procedure times_dz1     !dualz1*dualz1
     module procedure timesdz8_dz1  !dualz1 * complex8
     module procedure timesdz4_dz1  !dualz1 * complex4
     module procedure timesdzr8_dz1 !dualz1 * real8
     module procedure timesdzr4_dz1 !dualz1 * real4
     module procedure timesdzi8_dz1 !dualz1 * integer8
     module procedure timesdzi4_dz1 !dualz1 * integer4          
     module procedure timesz8d_dz1  !complex8 * dualz1
     module procedure timesz4d_dz1  !complex4 * dualz1
     module procedure timesr8dz_dz1 !real8 * dualz1
     module procedure timesr4dz_dz1 !real4 * dualz1
     module procedure timesi8dz_dz1 !integer8 * dualz1
     module procedure timesi4dz_dz1 !integer4 * dualz1
  end interface operator (*)

  !Power operator
  interface operator (**)
     module procedure power_dz1     !dualz1**dualz1
     module procedure powerdz8_dz1  !dualz1 ** complex8
     module procedure powerdz4_dz1  !dualz1 ** complex4
     module procedure powerdzr8_dz1 !dualz1 ** real8
     module procedure powerdzr4_dz1 !dualz1 ** real4
     module procedure powerdzi8_dz1 !dualz1 ** integer8
     module procedure powerdzi4_dz1 !dualz1 ** integer4          
     module procedure powerz8d_dz1  !complex8 ** dualz1
     module procedure powerz4d_dz1  !complex4 ** dualz1
     module procedure powerr8dz_dz1 !real8 ** dualz1
     module procedure powerr4dz_dz1 !real4 ** dualz1
     module procedure poweri8dz_dz1 !integer8 ** dualz1
     module procedure poweri4dz_dz1 !integer4 ** dualz1
  end interface operator (**)

  !Division operator
  interface operator (/)
     module procedure div_dz1     !dualz1 / dualz1
     module procedure divdz8_dz1  !dualz1 / complex8
     module procedure divdz4_dz1  !dualz1 / complex4
     module procedure divdzr8_dz1 !dualz1 / real8
     module procedure divdzr4_dz1 !dualz1 / real4
     module procedure divdzi8_dz1 !dualz1 / integer8
     module procedure divdzi4_dz1 !dualz1 / integer4          
     module procedure divz8d_dz1  !complex8 / dualz1
     module procedure divz4d_dz1  !complex4 / dualz1
     module procedure divr8dz_dz1 !real8 / dualz1
     module procedure divr4dz_dz1 !real4 / dualz1
     module procedure divi8dz_dz1 !integer8 / dualz1
     module procedure divi4dz_dz1 !integer4 / dualz1
  end interface operator (/)

  !Mathematical Functions
  interface sin
     module procedure sin_dz1
  end interface sin

  interface cos 
     module procedure cos_dz1
  end interface cos

  interface exp 
     module procedure exp_dz1
  end interface exp

  interface log
     module procedure log_dz1
  end interface log

  !fortran90 does not have tan (among others) for complex arguments
  interface tan
     module procedure tan_dz1
     module procedure tan_z 
  end interface tan

  interface sqrt
     module procedure sqrt_dz1
  end interface sqrt

  interface abs
     module procedure abs_dz1
  end interface abs

  interface asin
     module procedure asin_dz1
     module procedure asin_z
  end interface asin

  interface acos
     module procedure acos_dz1
     module procedure acos_z
  end interface acos

  interface atan
     module procedure atan_dz1
     module procedure atan_z 
  end interface atan

  interface sinh
     module procedure sinh_dz1
     module procedure sinh_z
  end interface sinh

  interface cosh
     module procedure cosh_dz1
     module procedure cosh_z
  end interface cosh

  interface tanh
     module procedure tanh_dz1
     module procedure tanh_z
  end interface tanh

  interface atanh
     module procedure atanh_dz1
     module procedure atanh_z
  end interface atanh

  interface asinh
     module procedure asinh_dz1
     module procedure asinh_z
  end interface asinh

  interface acosh
     module procedure acosh_dz1
     module procedure acosh_z
  end interface acosh

  interface atan2
     module procedure atan2_dz1
     module procedure atan2_z
  end interface atan2

  interface conjg
     module procedure conjg_dz1
  end interface conjg

  !some vector and matrix functions
  !Matrix multiplication
  interface matmul
     module procedure matmul22_dz1    !dual_rank2*dual_rank2     
     module procedure matmul_dr21_dz1 !dual_rank2*dual_rank1
     module procedure matmul_dr12_dz1 !dual_rank1*dual_rank2
     !
     module procedure matmul22_dz1C
     module procedure matmul22_dz1R
     module procedure matmul22_Cdz1
     module procedure matmul22_Rdz1
     !
     module procedure matmul_dr21_dz1C
     module procedure matmul_dr21_dz1R
     module procedure matmul_dr21_Cdz1
     module procedure matmul_dr21_Rdz1
     !
     module procedure matmul_dr12_dz1C
     module procedure matmul_dr12_dz1R
     module procedure matmul_dr12_Cdz1
     module procedure matmul_dr12_Rdz1
  end interface matmul

  interface sum
     module procedure sum2_dz1  !sum(dual_rank2,k)
     module procedure sum0_dz1  !sum(dual_rank2)
     module procedure sum10_dz1 !sum(dual_rank1)
  end interface sum

  !product for dual2 vectors x(k)
  interface product
     module procedure prod1_dz1
     module procedure prod0_dz1
     module procedure prod2_dz1 !product(dual_rank2,k)
  end interface product

  !<A|B>
  interface dot_product
     module procedure dot_product_dz1
  end interface dot_product

  !auxiliar function  
  interface todualz1
     module procedure dz1todualz1
     module procedure z8todualz1
     module procedure z4todualz1
     module procedure r8todualz1
     module procedure r4todualz1
     module procedure i8todualz1
     module procedure i4todualz1
  end interface todualz1
  !terminan interfaces

  !Functions and subroutines associated to the above interfaces
contains

  !Assignment, equal operator
  !dualz1 <--- complex8
  elemental subroutine igual1z_dz1(A, z)
    type(dualz1), intent(out) :: A
    complex(8), intent(in) :: z

    A%f0 = z
    A%f1 = 0
  end subroutine igual1z_dz1

  !Assignment, equal operator
  !dualz1 <--- complex4
  elemental subroutine igual1z4_dz1(A, z)
    type(dualz1), intent(out) :: A
    complex(4), intent(in) :: z    

    A%f0 = z
    A%f1 = 0
  end subroutine igual1z4_dz1

  !Assignment, equal operator
  !dualz1 <--- real8
  elemental subroutine igual1r_dz1(A, x)
    type(dualz1), intent(out) :: A
    real(8), intent(in) :: x

    A%f0 = x
    A%f1 = 0
  end subroutine igual1r_dz1

  !Assignment, equal operator
  !dualz1 <--- real4
  elemental subroutine igual1r4_dz1(A, x)
    type(dualz1), intent(out) :: A
    real(4), intent(in) :: x

    A%f0 = x
    A%f1 = 0
  end subroutine igual1r4_dz1

  !Assignment, equal operator
  !dualz1 <--- integer8
  elemental subroutine  igual1i8_dz1(A, i)
    type(dualz1), intent(out) :: A
    integer(8), intent(in) :: i

    A%f0 = i
    A%f1 = 0
  end subroutine igual1i8_dz1

  !Assignment, equal operator
  !dualz1 <--- integer4
  elemental subroutine  igual1i4_dz1(A, i)
    type(dualz1), intent(out) :: A
    integer(4), intent(in) :: i

    A%f0 = i
    A%f1 = 0
  end subroutine igual1i4_dz1

  !Logical equal operator
  elemental function eq_dz1(lhs, rhs) result(f_res)
    type (dualz1), intent(in) :: lhs, rhs
    logical :: f_res
    logical :: eqf0, eqf1

    eqf0 = lhs%f0 == rhs%f0
    eqf1 = lhs%f1 == rhs%f1

    f_res = all([eqf0,eqf1])
  end function eq_dz1

  !Logical not equal operator
  elemental function noteq_dz1(lhs, rhs) result(f_res)
    type (dualz1), intent(in) :: lhs, rhs
    logical :: f_res

    f_res = .not.(lhs == rhs)
  end function noteq_dz1

  !+dualz1 (unary)
  elemental function mas_dz1(A) result(f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res

    f_res = A
  end function mas_dz1

  !dualz1 + dualz1
  elemental function suma_dz1(A,B) result(f_s)
    type(dualz1), intent(in) :: A, B
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + B%f0
    f_s%f1 = A%f1 + B%f1    
  end function suma_dz1

  !dualz1 + complex8
  elemental function sumadz8_dz1(A,z8) result(f_s)
    type(dualz1), intent(in) :: A
    complex(8), intent(in) :: z8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + z8
    f_s%f1 = A%f1     
  end function sumadz8_dz1

  !dualz1 + complex4
  elemental function sumadz4_dz1(A,z4) result(f_s)
    type(dualz1), intent(in) :: A
    complex(4), intent(in) :: z4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + z4
    f_s%f1 = A%f1     
  end function sumadz4_dz1

  !dualz1 + real8
  elemental function sumadzr8_dz1(A,x8) result(f_s)
    type(dualz1), intent(in) :: A
    real(8), intent(in) :: x8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + x8
    f_s%f1 = A%f1     
  end function sumadzr8_dz1

  !dualz1 + real4
  elemental function sumadzr4_dz1(A,x4) result(f_s)
    type(dualz1), intent(in) :: A
    real(4), intent(in) :: x4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + x4
    f_s%f1 = A%f1     
  end function sumadzr4_dz1

  !dualz1 + integer8
  elemental function sumadzi8_dz1(A,i8) result(f_s)
    type(dualz1), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + i8
    f_s%f1 = A%f1     
  end function sumadzi8_dz1

  !dualz1 + integer4
  elemental function sumadzi4_dz1(A,i4) result(f_s)
    type(dualz1), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + i4
    f_s%f1 = A%f1     
  end function sumadzi4_dz1

  !complex8 + dualz1
  elemental function sumaz8d_dz1(zz8,A) result(f_s)
    type(dualz1), intent(in) :: A
    complex(8), intent(in) :: zz8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + zz8
    f_s%f1 = A%f1
  end function sumaz8d_dz1

  !complex4 + dualz1 
  elemental function sumaz4d_dz1(zz4,A) result(f_s)
    type(dualz1), intent(in) :: A
    complex(4), intent(in) :: zz4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + zz4
    f_s%f1 = A%f1     
  end function sumaz4d_dz1

  !real8 + dualz1
  elemental function sumar8dz_dz1(xx8,A) result(f_s)
    type(dualz1), intent(in) :: A
    real(8), intent(in) :: xx8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + xx8
    f_s%f1 = A%f1     
  end function sumar8dz_dz1

  !real4 + dualz1
  elemental function sumar4dz_dz1(xx4,A) result(f_s)
    type(dualz1), intent(in) :: A
    real(4), intent(in) :: xx4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + xx4
    f_s%f1 = A%f1     
  end function sumar4dz_dz1

  !integer8 + dualz1
  elemental function sumai8dz_dz1(ii8,A) result(f_s)
    type(dualz1), intent(in) :: A
    integer(8), intent(in) :: ii8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + ii8
    f_s%f1 = A%f1     
  end function sumai8dz_dz1

  !integer4 + dualz1
  elemental function sumai4dz_dz1(ii4,A) result(f_s)
    type(dualz1), intent(in) :: A
    integer(4), intent(in) :: ii4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 + ii4
    f_s%f1 = A%f1     
  end function sumai4dz_dz1

  !Subtraction 
  !-dualz1 (unary)
  elemental function menos_dz1(A) result(f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res

    f_res%f0 = -A%f0 
    f_res%f1 = -A%f1
  end function menos_dz1

  !dualz1 - dualz1
  elemental function resta_dz1(A,B) result(f_res)
    type(dualz1), intent(in) :: A, B
    type(dualz1) :: f_res

    f_res%f0 = A%f0 - B%f0
    f_res%f1 = A%f1 - B%f1
  end function resta_dz1

  !dualz1 - complex8
  elemental function restadz8_dz1(A,z8) result(f_s)
    type(dualz1), intent(in) :: A
    complex(8), intent(in) :: z8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 - z8
    f_s%f1 = A%f1     
  end function restadz8_dz1

  !dualz1 - complex4
  elemental function restadz4_dz1(A,z4) result(f_s)
    type(dualz1), intent(in) :: A
    complex(4), intent(in) :: z4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 - z4
    f_s%f1 = A%f1     
  end function restadz4_dz1

  !dualz1 - real8
  elemental function restadzr8_dz1(A,x8) result(f_s)
    type(dualz1), intent(in) :: A
    real(8), intent(in) :: x8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 - x8
    f_s%f1 = A%f1     
  end function restadzr8_dz1

  !dualz1 - real4
  elemental function restadzr4_dz1(A,x4) result(f_s)
    type(dualz1), intent(in) :: A
    real(4), intent(in) :: x4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 - x4
    f_s%f1 = A%f1     
  end function restadzr4_dz1

  !dualz1 - integer8
  elemental function restadzi8_dz1(A,i8) result(f_s)
    type(dualz1), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz1) :: f_s

    f_s%f0 = A%f0 - i8
    f_s%f1 = A%f1     
  end function restadzi8_dz1

  !dualz1 - integer4
  elemental function restadzi4_dz1(A,i4) result(f_s)
    type(dualz1), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz1) :: f_s

    f_s%f0 = A%f0 - i4
    f_s%f1 = A%f1     
  end function restadzi4_dz1

  !complex8 - dualz1
  elemental function restaz8d_dz1(zz8,A) result(f_s)
    type(dualz1), intent(in) :: A
    complex(8), intent(in) :: zz8
    type(dualz1) :: f_s

    f_s%f0 = zz8 - A%f0 
    f_s%f1 = -A%f1
  end function restaz8d_dz1

  !complex4 - dualz1 
  elemental function restaz4d_dz1(zz4,A) result(f_s)
    type(dualz1), intent(in) :: A
    complex(4), intent(in) :: zz4
    type(dualz1) :: f_s

    f_s%f0 = zz4 - A%f0
    f_s%f1 = -A%f1     
  end function restaz4d_dz1

  !real8 - dualz1
  elemental function restar8dz_dz1(xx8,A) result(f_s)
    type(dualz1), intent(in) :: A
    real(8), intent(in) :: xx8
    type(dualz1) :: f_s

    f_s%f0 = xx8 - A%f0
    f_s%f1 = -A%f1     
  end function restar8dz_dz1

  !real4 - dualz1
  elemental function restar4dz_dz1(xx4,A) result(f_s)
    type(dualz1), intent(in) :: A
    real(4), intent(in) :: xx4
    type(dualz1) :: f_s

    f_s%f0 = xx4 - A%f0
    f_s%f1 = -A%f1     
  end function restar4dz_dz1

  !integer8 - dualz1
  elemental function restai8dz_dz1(ii8,A) result(f_s)
    type(dualz1), intent(in) :: A
    integer(8), intent(in) :: ii8
    type(dualz1) :: f_s

    f_s%f0 = ii8 - A%f0
    f_s%f1 = -A%f1     
  end function restai8dz_dz1

  !integer4 - dualz1
  elemental function restai4dz_dz1(ii4,A) result(f_s)
    type(dualz1), intent(in) :: A
    integer(4), intent(in) :: ii4
    type(dualz1) :: f_s

    f_s%f0 = ii4 - A%f0
    f_s%f1 = -A%f1     
  end function restai4dz_dz1

  !Times
  !dualz1*dualz1
  elemental function times_dz1(A,B) result(f_res)
    type(dualz1), intent(in) :: A, B
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * B%f0
    f_res%f1 = A%f1 * B%f0 + A%f0 * B%f1
  end function times_dz1

  !dualz1 * complex8
  elemental function timesdz8_dz1(A,z8) result(f_res)
    type(dualz1), intent(in) :: A
    complex(8), intent(in) :: z8
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * z8
    f_res%f1 = A%f1 * z8
  end function timesdz8_dz1

  !dualz1 * complex4
  elemental function timesdz4_dz1(A,z4) result(f_res)
    type(dualz1), intent(in) :: A
    complex(4), intent(in) :: z4
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * z4
    f_res%f1 = A%f1 * z4
  end function timesdz4_dz1

  !dualz1 * real8
  elemental function timesdzr8_dz1(A,x8) result(f_res)
    type(dualz1), intent(in) :: A
    real(8), intent(in) :: x8
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * x8
    f_res%f1 = A%f1 * x8
  end function timesdzr8_dz1

  !dualz1 * real4
  elemental function timesdzr4_dz1(A,x4) result(f_res)
    type(dualz1), intent(in) :: A
    real(4), intent(in) :: x4
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * x4
    f_res%f1 = A%f1 * x4
  end function timesdzr4_dz1

  !dualz1 * integer8
  elemental function timesdzi8_dz1(A,i8) result(f_res)
    type(dualz1), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * i8
    f_res%f1 = A%f1 * i8
  end function timesdzi8_dz1

  !dualz1 * integer4
  elemental function timesdzi4_dz1(A,i4) result(f_res)
    type(dualz1), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * i4
    f_res%f1 = A%f1 * i4
  end function timesdzi4_dz1

  !complex8 * dualz1
  elemental function timesz8d_dz1(zz8,A) result(f_res)
    type(dualz1), intent(in) :: A
    complex(8), intent(in) :: zz8
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * zz8
    f_res%f1 = A%f1 * zz8
  end function timesz8d_dz1

  !complex4 * dualz1
  elemental function timesz4d_dz1(zz4,A) result(f_res)
    type(dualz1), intent(in) :: A
    complex(4), intent(in) :: zz4
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * zz4
    f_res%f1 = A%f1 * zz4
  end function timesz4d_dz1

  !real8 * dualz1
  elemental function timesr8dz_dz1(xx8,A) result(f_res)
    type(dualz1), intent(in) :: A
    real(8), intent(in) :: xx8
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * xx8
    f_res%f1 = A%f1 * xx8
  end function timesr8dz_dz1

  !real4 * dualz1
  elemental function timesr4dz_dz1(xx4,A) result(f_res)
    type(dualz1), intent(in) :: A
    real(4), intent(in) :: xx4
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * xx4
    f_res%f1 = A%f1 * xx4
  end function timesr4dz_dz1

  !integer8 * dualz1
  elemental function timesi8dz_dz1(ii8,A) result(f_res)
    type(dualz1), intent(in) :: A
    integer(8), intent(in) :: ii8
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * ii8
    f_res%f1 = A%f1 * ii8
  end function timesi8dz_dz1

  !integer4 * dualz1
  elemental function timesi4dz_dz1(ii4,A) result(f_res)
    type(dualz1), intent(in) :: A
    integer(4), intent(in) :: ii4
    type(dualz1) :: f_res

    f_res%f0 = A%f0 * ii4
    f_res%f1 = A%f1 * ii4
  end function timesi4dz_dz1

  !Power
  !dualz1**dualz1
  elemental function power_dz1(A,B) result(fr)
    type(dualz1), intent(in) :: A, B
    type(dualz1) :: fr
    complex(8) :: A0, A1, B0, B1

    A0 = A%f0; A1 = A%f1; B0 = B%f0; B1 = B%f1

    fr%f0 = A0**B0

    if(all([B0,B1]==0)) then
       fr%f1 = 0
    elseif(all([A0,A1]==0) .and. real(B0) > 0) then
       fr%f1 = 0 
    elseif(all([A1,B1]==0)) then
       fr%f1 = 0 
    elseif(A0 == 0 .and. real(B0)>1) then
       fr%f1 = 0 
    elseif(A0==0) then  
       if(B0 == 1) then
          fr%f1 = A1
       else
          fr%f1 = A0**(B0 - 1)*A1*B0 + A0**B0*B1*log(A0)
       end if
    else
       fr%f1 = A0**(B0 - 1)*A1*B0 + A0**B0*B1*log(A0)
    endif
  end function power_dz1

  !dualz1 ** complex8
  elemental function powerdz8_dz1(A,z8) result(f_res)
    type(dualz1), intent(in) :: A
    complex(8), intent(in) :: z8
    type(dualz1) :: f_res
    complex(8) :: A1, A0

    A0 = A%f0
    A1 = A%f1

    f_res%f0 = A0**z8
    if(z8==1)then
       f_res%f1 = A1
    else
       f_res%f1 = z8*A0**(z8 - 1)*A1
    end if    
  end function powerdz8_dz1

  !dualz1 ** complex4
  elemental function powerdz4_dz1(A,z4) result(f_res)
    type(dualz1), intent(in) :: A
    complex(4), intent(in) :: z4
    type(dualz1) :: f_res
    complex(8) :: A1, A0

    A0 = A%f0
    A1 = A%f1

    f_res%f0 = A0**z4
    if(z4==1) then
       f_res%f1 = A1
    else
       f_res%f1 = z4*A0**(z4 - 1)*A1
    end if    
  end function powerdz4_dz1

  !dualz1 ** real8
  elemental function powerdzr8_dz1(A,x8) result(f_res)
    type(dualz1), intent(in) :: A
    real(8), intent(in) :: x8
    type(dualz1) :: f_res
    complex(8) :: A1, A0

    A0 = A%f0
    A1 = A%f1

    f_res%f0 = A0**x8
    if(x8==1) then
       f_res%f1 = A1
    else
       f_res%f1 = x8*A0**(x8 - 1)*A1
    end if
  end function powerdzr8_dz1

  !dualz1 ** real4
  elemental function powerdzr4_dz1(A,x4) result(f_res)
    type(dualz1), intent(in) :: A
    real(4), intent(in) :: x4
    type(dualz1) :: f_res
    complex(8) :: A1, A0

    A0 = A%f0
    A1 = A%f1

    f_res%f0 = A0**x4
    if(x4==1) then
       f_res%f1 = A1
    else
       f_res%f1 = x4*A0**(x4 - 1)*A1
    end if
  end function powerdzr4_dz1

  !dualz1 ** integer8
  elemental function powerdzi8_dz1(A,i8) result(f_res)
    type(dualz1), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz1) :: f_res
    complex(8) :: A1, A0

    A0 = A%f0
    A1 = A%f1

    f_res%f0 = A0**i8
    if(i8==1)then
       f_res%f1 = A1
    else       
       f_res%f1 = i8*A0**(i8 - 1)*A1
    end if
  end function powerdzi8_dz1

  !dualz1 ** integer4
  elemental function powerdzi4_dz1(A,i4) result(f_res)
    type(dualz1), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz1) :: f_res
    complex(8) :: A1, A0

    A0 = A%f0
    A1 = A%f1

    f_res%f0 = A0**i4
    if(i4==1) then
       f_res%f1 = A1
    else
       f_res%f1 = i4*A0**(i4 - 1)*A1
    end if
  end function powerdzi4_dz1

  !complex8 ** dualz1
  elemental function powerz8d_dz1(zz8,B) result(f_res)
    type(dualz1), intent(in) :: B
    complex(8), intent(in) :: zz8
    type(dualz1) :: f_res
    type(dualz1) :: Az8

    Az8 = zz8
    !call igual1z_dz1(Az8,zz8)   
    f_res =  power_dz1(Az8,B)
  end function powerz8d_dz1

  !complex4 ** dualz1
  elemental function powerz4d_dz1(zz4,B) result(f_res)
    type(dualz1), intent(in) :: B
    complex(4), intent(in) :: zz4
    type(dualz1) :: f_res
    type(dualz1) :: Az4

    Az4 = zz4
    !call igual1z4_dz1(Az4,zz4)   
    f_res =  power_dz1(Az4,B)
  end function powerz4d_dz1

  !real8 ** dualz1
  elemental function powerr8dz_dz1(xx8,B) result(f_res)
    type(dualz1), intent(in) :: B
    real(8), intent(in) :: xx8
    type(dualz1) :: f_res
    type(dualz1) :: Ax8

    Ax8 = xx8
    !call igual1r_dz1(Ax8,xx8)
    f_res =  power_dz1(Ax8,B)
  end function powerr8dz_dz1

  !real4 ** dualz1
  elemental function powerr4dz_dz1(xx4,B) result(f_res)
    type(dualz1), intent(in) :: B
    real(4), intent(in) :: xx4
    type(dualz1) :: f_res
    type(dualz1) :: Ax4

    Ax4 = xx4
    !call igual1r4_dz1(Ax4,xx4)
    f_res =  power_dz1(Ax4,B)
  end function powerr4dz_dz1

  !integer8 ** dualz1
  elemental function poweri8dz_dz1(ii8,B) result(f_res)
    type(dualz1), intent(in) :: B
    integer(8), intent(in) :: ii8
    type(dualz1) :: f_res
    type(dualz1) :: Ai8

    Ai8 = ii8
    !call igual1i8_dz1(Ai8,ii8)
    f_res =  power_dz1(Ai8,B)
  end function poweri8dz_dz1

  !integer4 ** dualz1
  elemental function poweri4dz_dz1(ii4,B) result(f_res)
    type(dualz1), intent(in) :: B
    integer(4), intent(in) :: ii4
    type(dualz1) :: f_res
    type(dualz1) :: Ai4

    Ai4 = ii4
    !call igual1i4_dz1(Ai4,ii4)
    f_res =  power_dz1(Ai4,B)
  end function poweri4dz_dz1

  !Division
  !dualz1/dualz1
  !A/B = A*B^(-1)
  elemental function div_dz1(A,B) result(f_res)
    type(dualz1), intent(in) :: A, B
    type(dualz1) :: f_res
    type(dualz1):: BI
    complex(8) :: B0, B1

    B0 = B%f0; B1 = B%f1

    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2

    f_res = A*BI
  end function div_dz1

  !complex8 / dualz1
  elemental function divz8d_dz1(zz8,B) result(f_res)
    type(dualz1), intent(in) :: B
    complex(8), intent(in) :: zz8
    type(dualz1) :: f_res
    complex(8) :: B0, B1

    B0 = B%f0
    B1 = B%f1

    f_res%f0 =  zz8/B0 
    f_res%f1 = -zz8*B1/B0**2
  end function divz8d_dz1

  !complex4 / dualz1
  elemental function divz4d_dz1(zz4,B) result(f_res)
    type(dualz1), intent(in) :: B
    complex(4), intent(in) :: zz4
    type(dualz1) :: f_res
    complex(8) :: B0, B1

    B0 = B%f0
    B1 = B%f1

    f_res%f0 =  zz4/B0 
    f_res%f1 = -zz4*B1/B0**2
  end function divz4d_dz1

  !real8 / dualz1
  elemental function divr8dz_dz1(xx8,B) result(f_res)
    type(dualz1), intent(in) :: B
    real(8), intent(in) :: xx8
    type(dualz1) :: f_res
    complex(8) :: B0, B1

    B0 = B%f0
    B1 = B%f1

    f_res%f0 =  xx8/B0 
    f_res%f1 = -xx8*B1/B0**2
  end function divr8dz_dz1

  !real4 / dualz1
  elemental function divr4dz_dz1(xx4,B) result(f_res)
    type(dualz1), intent(in) :: B
    real(4), intent(in) :: xx4
    type(dualz1) :: f_res
    complex(8) :: B0, B1

    B0 = B%f0
    B1 = B%f1

    f_res%f0 =  xx4/B0 
    f_res%f1 = -xx4*B1/B0**2
  end function divr4dz_dz1

  !integer8 / dualz1
  elemental function divi8dz_dz1(ii8,B) result(f_res)
    type(dualz1), intent(in) :: B
    integer(8), intent(in) :: ii8
    type(dualz1) :: f_res
    complex(8) :: B0, B1

    B0 = B%f0
    B1 = B%f1

    f_res%f0 =  ii8/B0 
    f_res%f1 = -ii8*B1/B0**2
  end function divi8dz_dz1

  !integer4 / dualz1
  elemental function divi4dz_dz1(ii4,B) result(f_res)
    type(dualz1), intent(in) :: B
    integer(4), intent(in) :: ii4
    type(dualz1) :: f_res
    complex(8) :: B0, B1

    B0 = B%f0
    B1 = B%f1

    f_res%f0 =  ii4/B0 
    f_res%f1 = -ii4*B1/B0**2
  end function divi4dz_dz1

  !dualz1 / complex8
  elemental function divdz8_dz1(A,z8) result(f_res)
    type(dualz1), intent(in) :: A
    complex(8), intent(in) :: z8
    type(dualz1) :: f_res

    f_res%f0 = A%f0/z8
    f_res%f1 = A%f1/z8
  end function divdz8_dz1

  !dualz1 / complex4
  elemental function divdz4_dz1(A,z4) result(f_res)
    type(dualz1), intent(in) :: A
    complex(4), intent(in) :: z4
    type(dualz1) :: f_res

    f_res%f0 = A%f0/z4
    f_res%f1 = A%f1/z4
  end function divdz4_dz1

  !dualz1 / real8
  elemental function divdzr8_dz1(A,x8) result(f_res)
    type(dualz1), intent(in) :: A
    real(8), intent(in) :: x8
    type(dualz1) :: f_res

    f_res%f0 = A%f0/x8
    f_res%f1 = A%f1/x8
  end function divdzr8_dz1

  !dualz1 / real4
  elemental function divdzr4_dz1(A,x4) result(f_res)
    type(dualz1), intent(in) :: A
    real(4), intent(in) :: x4
    type(dualz1) :: f_res

    f_res%f0 = A%f0/x4
    f_res%f1 = A%f1/x4
  end function divdzr4_dz1

  !dualz1 / integer8
  elemental function divdzi8_dz1(A,i8) result(f_res)
    type(dualz1), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz1) :: f_res

    f_res%f0 = A%f0/i8
    f_res%f1 = A%f1/i8
  end function divdzi8_dz1

  !dualz1 / integer4
  elemental function divdzi4_dz1(A,i4) result(f_res)
    type(dualz1), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz1) :: f_res

    f_res%f0 = A%f0/i4
    f_res%f1 = A%f1/i4
  end function divdzi4_dz1

  !Mathematical Functions
  !Sine
  elemental function sin_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: A0, A1

    A0 = A%f0
    A1 = A%f1

    f_res%f0 = sin(A0)
    f_res%f1 = cos(A0)*A1
  end function sin_dz1

  !Cosine
  elemental function cos_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = cos(g0)
    f_res%f1 = -sin(g0)*g1
  end function cos_dz1

  !Exponential
  elemental function exp_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = exp(g0)
    f_res%f1 = exp(g0)*g1
  end function exp_dz1

  !Logarithm
  elemental function log_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = log(g0)
    f_res%f1 = g1/g0
  end function log_dz1

  !tan
  elemental function tan_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = tan(g0)
    f_res%f1 = g1/cos(g0)**2
  end function tan_dz1


  !Sqrt
  elemental function sqrt_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: A0, A1

    A0 = A%f0; A1 = A%f1

    f_res%f0 = sqrt(A0)
    f_res%f1 = A1/(2.0d0*sqrt(A0))
  end function sqrt_dz1

  !conjg
  !notice tat the conjugation operation is not differentiable. The below
  !definitions mean (df)*, not d(f*), which is some times useful
  elemental function conjg_dz1(g) result(f_r)
    type(dualz1), intent(in) :: g
    type(dualz1) :: f_r
    complex(8) :: g0, g1, r0, r1    

    g0 = g%f0
    g1 = g%f1

    r0 = conjg(g0)
    r1 = conjg(g1)

    f_r = dualz1(r0,r1)
  end function conjg_dz1

  !abs
  !keep in mind that the conjugation operation is not, in general,
  !differentiable
  elemental function abs_dz1(A) result(f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res

    f_res = sqrt(A*conjg(A))
  end function abs_dz1

  !Acos for complex arguments
  elemental function acos_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res

    f_res = asin(1d0) - asin(z)   
  end function acos_z

  !Asin for complex arguments
  elemental function asin_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res    
    real(8) :: x, y
    complex(8), parameter :: ii = cmplx(0d0,1d0,8)
    complex(8) :: xc, yc

    !Due to precision matters we use the below code. Fortran 2008 and
    !latter support this function, and thus this wouldn't be necessary
    x = real(z,8)
    y = aimag(z)
    xc = x
    yc = y

    if(y == 0d0) then
       f_res = arg((0,1)*xc + sqrt(1 - xc**2)) - (0,0.5)*log(xc**2 + &
            sqrt((-1 + xc**2)**2) + 2*xc*((-1 + xc**2)**2)**0.25*&
            sin(arg(1 - xc**2)/2.))
    else     
       f_res = arg(sqrt(1d0 - (xc + ii*yc)**2) + ii*(xc + ii*yc)) - ii*&
            log(sqrt((-yc + (4d0*x**2*yc**2 + (1d0 - xc**2 +           &
            yc**2)**2)**0.25*cos(arg(1d0 - (xc + ii*yc)**2)/2d0))**2 + &
            (xc + (4d0*xc**2*yc**2 + (1d0 - xc**2 + yc**2)**2)**0.25*  &
            sin(arg(1d0 - (xc + ii*yc)**2)/2d0))**2))
    end if
  end function asin_z

  !tan for complex arguments. 
  elemental function tan_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res

    f_res = sin(z)/cos(z)
  end function tan_z

  !Atan for complex arguments
  elemental function atan_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res
    complex(8), parameter :: ii = cmplx(0d0,1d0,8)
    real(8) :: x, y
    complex(8) :: xc, yc

    x = real(z,8)
    y = aimag(z)

    xc = x
    yc = y

    f_res = -arg(1d0 - ii*(xc + ii*yc))/2d0 + arg(1d0 + ii*(xc +       &
         ii*yc))/2d0 + ii*(-log(xc**2 + (1d0 - yc)**2)/4d0 +           &
         log(xc**2 + (1d0 + yc)**2)/4d0)
  end function atan_z

  !Atan2 for complex arguments
  elemental function atan2_z(zy,zx) result (f_res)
    complex(8), intent(in) :: zy, zx
    complex(8) :: f_res
    complex(8), parameter :: ii = cmplx(0d0,1d0,8)
    complex(8) :: num, den, divnd, t1, t2
    real(8) :: x1, x2, y1, y2
    complex(8) :: x1c, x2c, y1c, y2c

    x1 = real(zx,8)
    x2 = aimag(zx)

    y1 = real(zy,8)
    y2 = aimag(zy)

    x1c = x1
    x2c = x2
    y1c = y1
    y2c = y2

    num = x1c + ii*x2c + ii*(y1c + ii*y2c)
    den = sqrt((x1c + ii*x2c)**2 + (y1c + ii*y2c)**2)
    divnd = num/den;
    t1 = atan2(aimag(divnd),real(divnd,8))
    t2 = ii*log(sqrt((x2c + y1c)**2 + (x1c - y2c)**2)/((2d0*x1c*x2c +  &
         2d0*y1c*y2c)**2 + (x1c**2 - x2c**2 + y1c**2 -                 &
         y2c**2)**2)**0.25d0)

    f_res = t1 - t2
  end function atan2_z

  !tanh for complex arguments
  elemental function tanh_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res

    f_res = (exp(z) - exp(-z))/(exp(z) + exp(-z))
  end function tanh_z

  !Atanh for complex arguments
  elemental function atanh_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res
!!$    complex(8), parameter :: ii = cmplx(0d0,1d0,8)
!!$    complex(8), parameter :: dosi = cmplx(0d0,2d0,8)
!!$    real(8) :: x, y
!!$    complex(8) :: xc, yc
!!$
!!$    x = real(z,8)
!!$    y = aimag(z)
!!$
!!$    xc = x
!!$    yc = y
!!$    
!!$    f_res = (-dosi*arg(1d0 - xc - ii*yc) + dosi*Arg(1d0 + xc + ii*yc) -&
!!$         log((-1d0 + xc)**2 + yc**2) + log((1d0 + xc)**2 + y**2))/4d0
!!$    
    f_res = 0.5d0*(log(1 + z) - log(1 - z))
    !do not change Log A - Log B by Log(A/B)
  end function atanh_z

  !sinh for complex arguments
  elemental function sinh_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res

    f_res = 0.5d0*(exp(z)-exp(-z))
  end function sinh_z

  !asinh for complex arguments
  elemental function asinh_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res
    real(8) :: x, y
    complex(8) :: xc, yc

    x = real(z,8)
    y = aimag(z)

    xc = x
    yc = y

    f_res = (0,1)*arg(xc + sqrt(1 + (xc + (0,1)*yc)**2) + (0,1)*yc) +  &
         log(sqrt((xc + (4*xc**2*yc**2 + (1 + xc**2 - yc**2)**2)**0.25*&
         cos(arg(1 + (xc + (0,1)*yc)**2)/2d0))**2 + (yc + (4*xc**2*    &
         yc**2 + (1 + xc**2 - yc**2)**2)**0.25*sin(arg(1 + (xc + (0,1)*&
         yc)**2)/2d0))**2))

    !   f_res = log(z + sqrt(z*z + 1d0))
  end function asinh_z

  !cosh for complex arguments
  elemental function cosh_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res

    f_res = 0.5d0*(exp(z) + exp(-z))
  end function cosh_z

  !acosh for complex arguments
  elemental function acosh_z(z) result (f_res)
    complex(8), intent(in) :: z
    complex(8) :: f_res
    real(8) :: x, y
    complex(8) :: xc, yc

    x = real(z,8)
    y = aimag(z)

    xc = x
    yc = y

    f_res = ((0,2)*arg(xc + sqrt(-1 + xc + (0,1)*yc)*sqrt(1 + xc +&
         (0,1)*yc) + (0,1)*yc) + log(xc**2 + yc**2 + sqrt((-1 + xc)**2+&
         yc**2)*sqrt((1 + xc)**2 + yc**2) + 2*((-1 + xc)**2 +&
         yc**2)**0.25*((1 + xc)**2 + yc**2)**0.25*(xc*cos((arg(-1 + &
         xc + (0,1)*yc) + arg(1 + xc + (0,1)*yc))/2.) + yc*sin((arg(-1+&
         xc + (0,1)*yc) + arg(1 + xc + (0,1)*yc))/2d0))))/2d0

    !    f_res = log(z + sqrt(z - 1d0)*sqrt(z + 1d0))
    !do not change sqrt(z - 1d0)*sqrt(z + 1d0) by sqrt(z*z-1)
  end function acosh_z

  !Acos
  elemental function acos_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = acos(g0)
    f_res%f1 = -g1/sqrt(1d0 - g0**2)
  end function acos_dz1

  !Asin
  elemental function asin_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = asin(g0)
    f_res%f1 = g1/sqrt(1d0 - g0**2)
  end function asin_dz1

  !Atan
  elemental function atan_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = atan(g0)
    f_res%f1 = g1/(1d0 + g0**2)
  end function atan_dz1

  !sinh
  elemental function sinh_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = sinh(g0)
    f_res%f1 = cosh(g0)*g1
  end function sinh_dz1

  !cosh
  elemental function cosh_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = cosh(g0)
    f_res%f1 = sinh(g0)*g1
  end function cosh_dz1

  !tanh
  elemental function tanh_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = tanh(g0)
    f_res%f1 = g1/cosh(g0)**2
  end function tanh_dz1

  !atanh
  elemental function atanh_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = atanh(g0)
    f_res%f1 = g1/(1d0 - g0**2)
  end function atanh_dz1

  !asinh
  elemental function asinh_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1

    g0 = A%f0
    g1 = A%f1

    f_res%f0 = asinh(g0)
    f_res%f1 = g1/sqrt(1d0 + g0**2)
  end function asinh_dz1

  !acosh
  elemental function acosh_dz1(A) result (f_res)
    type(dualz1), intent(in) :: A
    type(dualz1) :: f_res
    complex(8) :: g0, g1, g0m, g0p, sg0m, sg0p

    g0 = A%f0
    g1 = A%f1

    g0m  = g0 - 1d0
    g0p  = g0 + 1d0
    sg0m = sqrt(g0m) 
    sg0p = sqrt(g0p)

    f_res%f0 = acosh(g0)
    f_res%f1 = g1/(sg0m*sg0p)
  end function acosh_dz1

  !Atan2
  elemental function atan2_dz1(A,B) result (f_res)
    type(dualz1), intent(in) :: A,B
    type(dualz1) :: f_res
    complex(8) :: A0, A1, B0, B1, r0, r1

    A0 = A%f0
    A1 = A%f1
    B0 = B%f0
    B1 = B%f1

    r0 = atan2(A0,B0)
    r1 = (A1*B0 - A0*B1)/(A0**2 + B0**2)

    f_res%f0 = r0
    f_res%f1 = r1
  end function atan2_dz1

  !matmul
  function matmul22_dz1(A,B) result(f_r)
    type(dualz1), intent(in), dimension(:,:) :: A, B
    type(dualz1), dimension(size(A,1),size(B,2)) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A0, A1
    complex(8), dimension(size(B,1),size(B,2)) :: B0, B1

    A0 = A%f0
    A1 = A%f1
    B0 = B%f0
    B1 = B%f1

    f_r%f0 = matmul(A0,B0)
    f_r%f1 = matmul(A1,B0) + matmul(A0,B1)
  end function matmul22_dz1

  function matmul22_dz1C(A,B) result(fr)
    type(dualz1), intent(in), dimension(:,:) :: A
    complex(8), intent(in), dimension(:,:) :: B
    type(dualz1), dimension(size(A,1),size(B,2)) :: fr
    type(dualz1), dimension(size(B,1),size(B,2)) :: Bdz1

    Bdz1 = B
    fr = matmul22_dz1(A,Bdz1)
  end function matmul22_dz1C

  function matmul22_dz1R(A,B) result(fr)
    type(dualz1), intent(in), dimension(:,:) :: A
    real(8), intent(in), dimension(:,:) :: B
    type(dualz1), dimension(size(A,1),size(B,2)) :: fr
    type(dualz1), dimension(size(B,1),size(B,2)) :: Bdz1

    Bdz1 = B
    fr = matmul22_dz1(A,Bdz1)
  end function matmul22_dz1R

  function matmul22_Cdz1(A,B) result(fr)
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz1), intent(in), dimension(:,:) :: B
    type(dualz1), dimension(size(A,1),size(B,2)) :: fr
    type(dualz1), dimension(size(A,1),size(A,2)) :: Adz1

    Adz1 = A
    fr = matmul22_dz1(Adz1,B)
  end function matmul22_Cdz1

  function matmul22_Rdz1(A,B) result(fr)
    real(8), intent(in), dimension(:,:) :: A
    type(dualz1), intent(in), dimension(:,:) :: B
    type(dualz1), dimension(size(A,1),size(B,2)) :: fr
    type(dualz1), dimension(size(A,1),size(A,2)) :: Adz1

    Adz1 = A
    fr = matmul22_dz1(Adz1,B)
  end function matmul22_Rdz1 

  !rank2*rank1
  function matmul_dr21_dz1(A,x) result(f_r)
    type(dualz1), intent(in) :: A(:,:)
    type(dualz1), intent(in), dimension(size(A,2)) :: x
    type(dualz1), dimension(size(A,1)) :: f_r
    integer :: m,n

    m = size(A,1)
    n = size(A,2)
    f_r = reshape(matmul(A,reshape(x,[n,1])),[m])
  end function matmul_dr21_dz1

  function matmul_dr21_dz1C(A,x) result(fr)
    type(dualz1), intent(in), dimension(:,:) :: A
    complex(8), intent(in), dimension(size(A,2)) :: x
    type(dualz1), dimension(size(A,1)) :: fr
    type(dualz1), dimension(size(A,2)) :: xd

    xd = x
    fr = matmul_dr21_dz1(A,xd)
  end function matmul_dr21_dz1C

  function matmul_dr21_dz1R(A,x) result(fr)
    type(dualz1), intent(in), dimension(:,:) :: A
    real(8), intent(in), dimension(size(A,2)) :: x
    type(dualz1), dimension(size(A,1)) :: fr
    type(dualz1), dimension(size(A,2)) :: xd

    xd = x
    fr = matmul_dr21_dz1(A,xd)
  end function matmul_dr21_dz1R

  function matmul_dr21_Cdz1(A,x) result(fr)
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz1), intent(in), dimension(size(A,2)) :: x
    type(dualz1), dimension(size(A,1)) :: fr
    type(dualz1), dimension(size(A,1),size(A,2)) :: Adz1

    Adz1 = A
    fr = matmul_dr21_dz1(Adz1,x)
  end function matmul_dr21_Cdz1

  function matmul_dr21_Rdz1(A,x) result(fr)
    real(8), intent(in), dimension(:,:) :: A
    type(dualz1), intent(in), dimension(size(A,2)) :: x
    type(dualz1), dimension(size(A,1)) :: fr
    type(dualz1), dimension(size(A,1),size(A,2)) :: Adz1

    Adz1 = A
    fr = matmul_dr21_dz1(Adz1,x)
  end function matmul_dr21_Rdz1 

  !in order to avoid ambiguity we must change x to xx
  !"A generic function must be able to distinguish its arguments by 
  !type AND by name"
  function matmul_dr12_dz1(xx,A) result(f_r)
    type(dualz1), intent(in) :: A(:,:)
    type(dualz1), intent(in), dimension(:) :: xx
    type(dualz1), dimension(size(A,2)) :: f_r
    integer :: m,n

    m = size(A,1)
    n = size(A,2)
    !if size(xx) /= m you may want to use some warning, the below code
    !reliess on the reshape fortran function which uses padding
    f_r = reshape(matmul(reshape(xx,[1,m]),A),[n])
  end function matmul_dr12_dz1

  function matmul_dr12_dz1C(xx,A) result(fr)
    type(dualz1), intent(in), dimension(:) :: xx
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz1), dimension(size(A,2)) :: fr
    type(dualz1), dimension(size(A,1),size(A,2)) :: Adz1

    Adz1 = A
    fr =  matmul_dr12_dz1(xx,Adz1)
  end function matmul_dr12_dz1C

  function matmul_dr12_dz1R(xx,A) result(fr)
    type(dualz1), intent(in), dimension(:) :: xx
    real(8), intent(in), dimension(:,:) :: A
    type(dualz1), dimension(size(A,2)) :: fr
    type(dualz1), dimension(size(A,1),size(A,2)) :: Adz1

    Adz1 = A
    fr =  matmul_dr12_dz1(xx,Adz1)
  end function matmul_dr12_dz1R

  function matmul_dr12_Cdz1(xx,A) result(fr)
    complex(8), intent(in), dimension(:) :: xx
    type(dualz1), intent(in), dimension(:,:) :: A
    type(dualz1), dimension(size(A,2)) :: fr
    type(dualz1), dimension(size(xx)) :: xxdz1

    xxdz1 = xx
    fr =  matmul_dr12_dz1(xxdz1,A)
  end function matmul_dr12_Cdz1

  function matmul_dr12_Rdz1(xx,A) result(fr)
    real(8), intent(in), dimension(:) :: xx
    type(dualz1), intent(in), dimension(:,:) :: A
    type(dualz1), dimension(size(A,2)) :: fr
    type(dualz1), dimension(size(xx)) :: xxdz1

    xxdz1 = xx
    fr =  matmul_dr12_dz1(xxdz1,A)
  end function matmul_dr12_Rdz1
  
  !sum for rank2 arrays
  !the result is given as array of rank 1
  function sum2_dz1(A,k) result(f_r)
    type(dualz1), intent(in), dimension(:,:) :: A
    integer, intent(in) :: k
    type(dualz1), dimension(size(A,2/k)) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A0, A1

    A0 = A%f0
    A1 = A%f1

    f_r%f0 = sum(A0,k)
    f_r%f1 = sum(A1,k)
  end function sum2_dz1

  !sum(A) rank2
  function sum0_dz1(A) result(f_r)
    type(dualz1), intent(in), dimension(:,:) :: A
    type(dualz1) :: f_r
    complex(8) :: fr0, fr1
    complex(8), dimension(size(A,1),size(A,2)) :: A0, A1

    A0 = A%f0
    A1 = A%f1

    fr0 = sum(A0)
    fr1 = sum(A1)

    f_r%f0 = fr0
    f_r%f1 = fr1
  end function sum0_dz1

  !sum(A) rank1
  function sum10_dz1(A) result(f_r)
    type(dualz1), intent(in), dimension(:) :: A
    type(dualz1) :: f_r
    complex(8) :: fr0, fr1
    complex(8), dimension(size(A)) :: A0, A1

    A0 = A%f0
    A1 = A%f1

    fr0 = sum(A0)
    fr1 = sum(A1)

    f_r%f0 = fr0
    f_r%f1 = fr1
  end function sum10_dz1

  !product function for dualz1 vectors (arrays of rank 1)
  function prod1_dz1(x) result(f_r)
    type(dualz1), intent(in), dimension(:) :: x
    type(dualz1) :: f_r
    integer :: k

    f_r = 1d0
    do k = 1, size(x)
       f_r = f_r * x(k)
    end do
  end function prod1_dz1

  !product(x) for arrays of rank 2
  function prod0_dz1(x) result(f_r)
    type(dualz1), intent(in), dimension(:,:) :: x
    type(dualz1) :: f_r

    f_r = prod1_dz1(reshape(x,[size(x)]))
  end function prod0_dz1

  !product(A,k) for rank2 arrays
  !the result is given as array of rank 1
  function prod2_dz1(x,c) result(f_r)
    type(dualz1), intent(in), dimension(:,:) :: x
    integer, intent(in) :: c
    type(dualz1), dimension(size(x,2/c)) :: f_r
    type(dualz1), dimension(size(x,c)) :: xc
    integer :: k

    if(c==1) then
       do k = 1, size(x,2)
          xc = x(:,k)
          f_r(k) = prod1_dz1(xc)
       end do
    else if(c==2) then
       do k = 1, size(x,1)
          xc = x(k,:)
          f_r(k) = prod1_dz1(xc)
       end do
    else 
       stop 'c must be equal to 1 (2) to collapse rows (columns)'
    end if
  end function prod2_dz1

  !|A> = |A0> + |A1> eps; |B> = |B0> + |B1> eps;
  !<A|B> = <A0|B0> + (<A1|B0> + <A0|B1>) eps
  !As in the conjg and abs functions, here, the dual parts are (df)* not
  !d(f*)
  function dot_product_dz1(A,B) result(f_r)
    type(dualz1), intent(in), dimension(:) :: A, B
    type(dualz1) :: f_r
    complex(8), dimension(size(A)) :: A0, A1, B0, B1

    A0 = A%f0
    A1 = A%f1
    B0 = B%f0
    B1 = B%f1

    f_r%f0 = dot_product(A0,B0)
    f_r%f1 = dot_product(A1,B0) + dot_product(A0,B1)

  end function dot_product_dz1

  !una implementacion directa que usa menos memoria pero pudiera ser
  !menos eficiente. Usando el compilador ifort no hay diferencia
!!$  function dot_product_dz1(A,B) result(f_r)
!!$    type(dualz1), intent(in), dimension(:) :: A, B
!!$    type(dualz1) :: f_r
!!$    integer :: k
!!$
!!$    f_r = 0
!!$    do k=1, size(A)
!!$       f_r = f_r + conjg(A(k))*B(k)
!!$    end do
!!$  end function dot_product_dz1

  !=====================================================================
  !auxiliar function to convert to dualz1
  !in principle this function is not necessary (the assignment operator
  !does the "same" job) but is provided for convenience
  elemental function dz1todualz1(z) result(fr)
    type(dualz1), intent(in) :: z
    type(dualz1) :: fr

    fr = z
  end function dz1todualz1

  elemental function z8todualz1(z) result(fr)
    complex(8), intent(in) :: z
    type(dualz1) :: fr

    fr = z
  end function z8todualz1

  elemental function z4todualz1(z) result(fr)
    complex(4), intent(in) :: z
    type(dualz1) :: fr

    fr = z
  end function z4todualz1

  elemental function r8todualz1(z) result(fr)
    real(8), intent(in) :: z
    type(dualz1) :: fr

    fr = z
  end function r8todualz1

  elemental function r4todualz1(z) result(fr)
    real(4), intent(in) :: z
    type(dualz1) :: fr

    fr = z
  end function r4todualz1

  elemental function i8todualz1(z) result(fr)
    integer(8), intent(in) :: z
    type(dualz1) :: fr

    fr = z
  end function i8todualz1

  elemental function i4todualz1(z) result(fr)
    integer(4), intent(in) :: z
    type(dualz1) :: fr

    fr = z
  end function i4todualz1
  !=====================================================================
  
  !auxiliar functions
  elemental function f0part_dz1(A) result(fr)
    type(dualz1), intent(in) :: A
    complex(8) :: fr

    fr = A%f0
  end function f0part_dz1

  elemental function f1part_dz1(A) result(fr)
    type(dualz1), intent(in) :: A
    complex(8) :: fr
    
    fr = A%f1
  end function f1part_dz1

  !auxiliar function (not public)
  elemental function arg(z) result(f_r)
    complex(8), intent(in) :: z
    real(8) :: f_r
    real(8) :: x, y

    x = real(z,8)
    y = aimag(z)

    f_r = atan2(y,x)
  end function arg
end module dualz1_mod
