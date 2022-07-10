!Dual numbers of complex components
!This can be used to compute 0-th, 1-th, 2-th, 3-th and 4-th order
!derivatives.

!F. Peñuñuri
!UADY, Merida Yucatan Mexico
!2022
module dualz4_mod
  use dualz1_mod
  use dualz2_mod
  use dualz3_mod
  implicit none
  private

  !dualz1 ---> [f0,f1]
  !dualz2 ---> [f0,f1,f2]
  !dualz3 ---> [f0,f1,f2,f3]
  !dualz4 ---> [f0,f1,f2,f3,f4]
  type, public, extends (dualz3) :: dualz4
     complex(8) :: f4
  end type dualz4

  public :: assignment (=)
  public :: operator(==), operator(/=), operator(+), operator(-)
  public :: operator(*), operator (/), operator(**)

  public :: sin, cos, exp, log, tan, sinh, cosh, tanh, asin, acos, atan
  public :: acosh, asinh, atanh, atan2, sqrt, conjg, abs, dot_product
  public :: matmul, sum, product

  public :: f0part_dz4, f1part_dz4, f2part_dz4, f3part_dz4, f4part_dz4
  public :: todualz4
  
  !INTERFACES
  !equal assignment
  interface assignment (=)
     module procedure igualz3_dz4  !dualz4 <--- dualz3
     module procedure igualz2_dz4  !dualz4 <--- dualz2
     module procedure igualz1_dz4  !dualz4 <--- dualz1
     module procedure igualc8_dz4  !dualz4 <--- complex8
     module procedure igualc4_dz4  !dualz4 <--- complex4
     module procedure igualr8_dz4  !dualz4 <--- real8
     module procedure igualr4_dz4  !dualz4 <--- real4
     module procedure iguali8_dz4  !dualz4 <--- integer8
     module procedure iguali4_dz4  !dualz4 <--- integer4
  end interface assignment (=)
  
  !Logical equal operator 
  interface operator (==) 
     module procedure eq_dz4  !dualz4 == dualz4
  end interface operator (==)

  !Logical not equal operator 
  interface operator (/=) 
     module procedure noteq_dz4  !dualz4 /= dualz4
  end interface operator (/=)
  
  !Plus operator
  interface operator (+)
     module procedure mas_dz4      !unary 
     module procedure suma_dz4     !dualz4 + dualz4
     module procedure sumaz4z3_dz4 !dualz4 + dualz3
     module procedure sumaz4z2_dz4 !dualz4 + dualz2
     module procedure sumaz4z1_dz4 !dualz4 + dualz1
     module procedure sumaz4c8_dz4 !dualz4 + complex8
     module procedure sumaz4c4_dz4 !dualz4 + complex4
     module procedure sumaz4r8_dz4 !dualz4 + real8
     module procedure sumaz4r4_dz4 !dualz4 + real4
     module procedure sumaz4i8_dz4 !dualz4 + integer8
     module procedure sumaz4i4_dz4 !dualz4 + integer4
     module procedure sumaz3z4_dz4 !dualz3 + dualz4
     module procedure sumaz2z4_dz4 !dualz2 + dualz4
     module procedure sumaz1z4_dz4 !dualz1 + dualz4
     module procedure sumac8z4_dz4 !complex8 + dualz4
     module procedure sumac4z4_dz4 !complex4 + dualz4
     module procedure sumar8z4_dz4 !real8 + dualz4
     module procedure sumar4z4_dz4 !real4 + dualz4
     module procedure sumai8z4_dz4 !integer8 + dualz4
     module procedure sumai4z4_dz4 !integer4 + dualz4
  end interface operator (+)

  !Minus operator
  interface operator (-)
     module procedure menos_dz4     !unary 
     module procedure resta_dz4     !dualz4 - dualz4
     module procedure restaz4z3_dz4 !dualz4 - dualz3
     module procedure restaz4z2_dz4 !dualz4 - dualz2
     module procedure restaz4z1_dz4 !dualz4 - dualz1
     module procedure restaz4c8_dz4 !dualz4 - complex8
     module procedure restaz4c4_dz4 !dualz4 - complex4
     module procedure restaz4r8_dz4 !dualz4 - real8
     module procedure restaz4r4_dz4 !dualz4 - real4
     module procedure restaz4i8_dz4 !dualz4 - integer8
     module procedure restaz4i4_dz4 !dualz4 - integer4
     module procedure restaz3z4_dz4 !dualz3 - dualz4
     module procedure restaz2z4_dz4 !dualz2 - dualz4
     module procedure restaz1z4_dz4 !dualz1 - dualz4
     module procedure restac8z4_dz4 !complex8 - dualz4
     module procedure restac4z4_dz4 !complex4 - dualz4
     module procedure restar8z4_dz4 !real8 - dualz4
     module procedure restar4z4_dz4 !real4 - dualz4
     module procedure restai8z4_dz4 !integer8 - dualz4
     module procedure restai4z4_dz4 !integer4 - dualz4
  end interface operator (-)

  !Times operator
  interface operator (*)
     module procedure timesz4z4_dz4 !dualz4 * dualz4
     module procedure timesz4z3_dz4 !dualz4 * dualz3
     module procedure timesz4z2_dz4 !dualz4 * dualz2
     module procedure timesz4z1_dz4 !dualz4 * dualz1
     module procedure timesz4c8_dz4 !dualz4 * complex8
     module procedure timesz4c4_dz4 !dualz4 * complex4
     module procedure timesz4r8_dz4 !dualz4 * real8
     module procedure timesz4r4_dz4 !dualz4 * real4
     module procedure timesz4i8_dz4 !dualz4 * integer8
     module procedure timesz4i4_dz4 !dualz4 * integer4
     module procedure timesz3z4_dz4 !dualz3 * dualz4
     module procedure timesz2z4_dz4 !dualz2 * dualz4
     module procedure timesz1z4_dz4 !dualz1 * dualz4
     module procedure timesc8z4_dz4 !complex8 * dualz4
     module procedure timesc4z4_dz4 !complex4 * dualz4
     module procedure timesr8z4_dz4 !real8 * dualz4
     module procedure timesr4z4_dz4 !real4 * dualz4
     module procedure timesi8z4_dz4 !integer8 * dualz4
     module procedure timesi4z4_dz4 !integer4 * dualz4
  end interface operator (*)

  !Division operator
  interface operator (/)
     module procedure divz4z4_dz4 !dualz4 / dualz4
     module procedure divz4z3_dz4 !dualz4 / dualz3
     module procedure divz4z2_dz4 !dualz4 / dualz2
     module procedure divz4z1_dz4 !dualz4 / dualz1
     module procedure divz4c8_dz4 !dualz4 / complex8
     module procedure divz4c4_dz4 !dualz4 / complex4
     module procedure divz4r8_dz4 !dualz4 / real8
     module procedure divz4r4_dz4 !dualz4 / real4
     module procedure divz4i8_dz4 !dualz4 / integer8
     module procedure divz4i4_dz4 !dualz4 / integer4
     module procedure divz3z4_dz4 !dualz3 / dualz4
     module procedure divz2z4_dz4 !dualz2 / dualz4
     module procedure divz1z4_dz4 !dualz1 / dualz4
     module procedure divc8z4_dz4 !complex8 / dualz4
     module procedure divc4z4_dz4 !complex4 / dualz4
     module procedure divr8z4_dz4 !real8 / dualz4
     module procedure divr4z4_dz4 !real4 / dualz4
     module procedure divi8z4_dz4 !integer8 / dualz4
     module procedure divi4z4_dz4 !integer4 / dualz4
  end interface operator (/)

  !Power operator
  interface operator (**)
     module procedure powerz4z4_dz4 !dualz4 / dualz4
     module procedure powerz4z3_dz4 !dualz4 / dualz3
     module procedure powerz4z2_dz4 !dualz4 / dualz2
     module procedure powerz4z1_dz4 !dualz4 / dualz1
     module procedure powerz4c8_dz4 !dualz4 / complex8
     module procedure powerz4c4_dz4 !dualz4 / complex4
     module procedure powerz4r8_dz4 !dualz4 / real8
     module procedure powerz4r4_dz4 !dualz4 / real4
     module procedure powerz4i8_dz4 !dualz4 / integer8
     module procedure powerz4i4_dz4 !dualz4 / integer4
     module procedure powerz3z4_dz4 !dualz3 / dualz4
     module procedure powerz2z4_dz4 !dualz2 / dualz4
     module procedure powerz1z4_dz4 !dualz1 / dualz4
     module procedure powerc8z4_dz4 !complex8 / dualz4
     module procedure powerc4z4_dz4 !complex4 / dualz4
     module procedure powerr8z4_dz4 !real8 / dualz4
     module procedure powerr4z4_dz4 !real4 / dualz4
     module procedure poweri8z4_dz4 !integer8 / dualz4
     module procedure poweri4z4_dz4 !integer4 / dualz4
  end interface operator (**)

  !Mathematical Functions
  interface sin
     module procedure sin_dz4
  end interface sin

  interface cos
     module procedure cos_dz4
  end interface cos
  
  interface exp
     module procedure exp_dz4
  end interface exp

  interface log
     module procedure log_dz4
  end interface log

  interface tan
     module procedure tan_dz4
  end interface tan

  interface sinh
     module procedure sinh_dz4
  end interface sinh

  interface cosh
     module procedure cosh_dz4
  end interface cosh

  interface tanh
     module procedure tanh_dz4
  end interface tanh

  interface asin
     module procedure asin_dz4
  end interface asin

  interface acos
     module procedure acos_dz4
  end interface acos

  interface atan
     module procedure atan_dz4
  end interface atan

  interface acosh
     module procedure acosh_dz4
  end interface acosh

  interface asinh
     module procedure asinh_dz4
  end interface asinh

  interface atanh
     module procedure atanh_dz4
  end interface atanh

  interface atan2
     module procedure atan2_dz4
  end interface atan2

  interface sqrt
     module procedure sqrt_dz4
  end interface sqrt

  interface conjg
     module procedure conjg_dz4
  end interface conjg

  interface abs
     module procedure abs_dz4
  end interface abs

  !<A|B>
  interface dot_product
     module procedure dot_product_dz4
  end interface dot_product

  !Matrix multiplication
  interface matmul
     module procedure matmul22_dz4    !dual_rank2*dual_rank2
     module procedure matmul_dr21_dz4 !dual_rank2*dual_rank1
     module procedure matmul_dr12_dz4 !dual_rank1*dual_rank2
     !
     module procedure matmul22_dz4C
     module procedure matmul22_dz4R
     module procedure matmul22_Cdz4
     module procedure matmul22_Rdz4
     !
     module procedure matmul_dr21_dz4C
     module procedure matmul_dr21_dz4R
     module procedure matmul_dr21_Cdz4
     module procedure matmul_dr21_Rdz4
     !
     module procedure matmul_dr12_dz4C
     module procedure matmul_dr12_dz4R
     module procedure matmul_dr12_Cdz4
     module procedure matmul_dr12_Rdz4
  end interface matmul

  interface sum
     module procedure sum2_dz4  !sum(dual_rank2,k)
     module procedure sum0_dz4  !sum(dual_rank2)
     module procedure sum10_dz4 !sum(dual_rank1)
  end interface sum

  !product for dual vectors x(k)
  interface product
     module procedure prod1_dz4
     module procedure prod0_dz4
     module procedure prod2_dz4 !product(dual_rank2,k)
  end interface product

  !auxiliar function  
  interface todualz4
     module procedure dz4todualz4
     module procedure dz3todualz4
     module procedure dz2todualz4
     module procedure dz1todualz4
     module procedure z8todualz4
     module procedure z4todualz4
     module procedure r8todualz4
     module procedure r4todualz4
     module procedure i8todualz4
     module procedure i4todualz4
  end interface todualz4
  !terminan interfaces

  !Functions and subroutines associated to the above interfaces
contains

  !Assignment, equal operator
  !dualz4 <--- dualz3
  elemental subroutine igualz3_dz4(A, z)
    type(dualz4), intent(out) :: A
    type(dualz3), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine igualz3_dz4

  !dualz4 <--- dualz2
  elemental subroutine igualz2_dz4(A, z)
    type(dualz4), intent(out) :: A
    type(dualz2), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine igualz2_dz4

  !dualz4 <--- dualz1
  elemental subroutine igualz1_dz4(A, z)
    type(dualz4), intent(out) :: A
    type(dualz1), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine igualz1_dz4

  !dualz4 <--- c8
  elemental subroutine igualc8_dz4(A, z)
    type(dualz4), intent(out) :: A
    complex(8), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine igualc8_dz4

  !dualz4 <--- c4
  elemental subroutine igualc4_dz4(A, z)
    type(dualz4), intent(out) :: A
    complex(4), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine igualc4_dz4
  
  !dualz4 <--- r8
  elemental subroutine igualr8_dz4(A, z)
    type(dualz4), intent(out) :: A
    real(8), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine igualr8_dz4

  !dualz4 <--- r4
  elemental subroutine igualr4_dz4(A, z)
    type(dualz4), intent(out) :: A
    real(4), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine igualr4_dz4

  !dualz4 <--- i8
  elemental subroutine iguali8_dz4(A, z)
    type(dualz4), intent(out) :: A
    integer(8), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine iguali8_dz4

  !dualz4 <--- i4
  elemental subroutine iguali4_dz4(A, z)
    type(dualz4), intent(out) :: A
    integer(4), intent(in) :: z

    A%dualz3 = z
    A%f4 = 0
  end subroutine iguali4_dz4

  !Logical not equal operator
  elemental function noteq_dz4(lhs, rhs) result(f_res)
    type (dualz4), intent(in) :: lhs, rhs
    logical :: f_res

    f_res = .not.(lhs == rhs)
  end function noteq_dz4

  !Logical equal operator
  elemental function eq_dz4(lhs, rhs) result(f_res)
    type (dualz4), intent(in) :: lhs, rhs
    logical :: f_res
    logical :: eqf0, eqf1, eqf2, eqf3, eqf4

    eqf0 = lhs%f0 == rhs%f0
    eqf1 = lhs%f1 == rhs%f1
    eqf2 = lhs%f2 == rhs%f2
    eqf3 = lhs%f3 == rhs%f3
    eqf4 = lhs%f4 == rhs%f4

    f_res = all([eqf0,eqf1,eqf2,eqf3,eqf4])
  end function eq_dz4

  !dualz4 + dualz4
  elemental function suma_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A, B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B%dualz3
    fr%f4 = A%f4 + B%f4
  end function suma_dz4

  !dualz4 + dualz3
  elemental function sumaz4z3_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz3), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4z3_dz4
  
  !dualz4 + dualz2
  elemental function sumaz4z2_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4z2_dz4

  !dualz4 + dualz1
  elemental function sumaz4z1_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4z1_dz4

  !dualz4 + c8
  elemental function sumaz4c8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4c8_dz4

  !dualz4 + c4
  elemental function sumaz4c4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4c4_dz4

  !dualz4 + r8
  elemental function sumaz4r8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4r8_dz4

 !dualz4 + r4
  elemental function sumaz4r4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4r4_dz4
 
  !dualz4 + i8
  elemental function sumaz4i8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4i8_dz4

  !dualz4 + i4
  elemental function sumaz4i4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz4) :: fr

    fr%dualz3 = A%dualz3 + B
    fr%f4 = A%f4
  end function sumaz4i4_dz4

  !dualz3 + dualz4
  elemental function sumaz3z4_dz4(B,A) result(fr)
    type(dualz3), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4z3_dz4(A,B)
  end function sumaz3z4_dz4

  !dualz2 + dualz4
  elemental function sumaz2z4_dz4(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4z2_dz4(A,B)
  end function sumaz2z4_dz4

  !dualz1 + dualz4
  elemental function sumaz1z4_dz4(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4z1_dz4(A,B)
  end function sumaz1z4_dz4

  !c8 + dualz4
  elemental function sumac8z4_dz4(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4c8_dz4(A,B)
  end function sumac8z4_dz4

  !c4 + dualz4
  elemental function sumac4z4_dz4(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4c4_dz4(A,B)
  end function sumac4z4_dz4

  !r8 + dualz4
  elemental function sumar8z4_dz4(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4r8_dz4(A,B)
  end function sumar8z4_dz4

  !r4 + dualz4
  elemental function sumar4z4_dz4(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4r4_dz4(A,B)
  end function sumar4z4_dz4

  !i8 + dualz4
  elemental function sumai8z4_dz4(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4i8_dz4(A,B)
  end function sumai8z4_dz4
  
  !i4 + dualz4
  elemental function sumai4z4_dz4(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = sumaz4i4_dz4(A,B)
  end function sumai4z4_dz4

  !dualz4 - dualz4
  elemental function resta_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A, B
    type(dualz4) :: fr

    fr = A + (-B)
  end function resta_dz4

  !dualz4 - dualz3
  elemental function restaz4z3_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz3), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4z3_dz4

  !dualz4 - dualz2
  elemental function restaz4z2_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4z2_dz4

  !dualz4 - dualz1
  elemental function restaz4z1_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4z1_dz4

  !dualz4 - c8
  elemental function restaz4c8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4c8_dz4

  !dualz4 - c4
  elemental function restaz4c4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4c4_dz4

  !dualz4 - r8
  elemental function restaz4r8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4r8_dz4
 
  !dualz4 - r4
  elemental function restaz4r4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4r4_dz4

  !dualz4 - i8
  elemental function restaz4i8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4i8_dz4
  
  !dualz4 - i4
  elemental function restaz4i4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz4) :: fr

    fr = A + (-B)
  end function restaz4i4_dz4

  !dualz3 - dualz4
  elemental function restaz3z4_dz4(B,A) result(fr)
    type(dualz3), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restaz3z4_dz4

 !dualz2 - dualz4
  elemental function restaz2z4_dz4(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restaz2z4_dz4

 !dualz1 - dualz4
  elemental function restaz1z4_dz4(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restaz1z4_dz4
  
  !c8 - dualz4
  elemental function restac8z4_dz4(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restac8z4_dz4
  
  !c4 - dualz4
  elemental function restac4z4_dz4(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restac4z4_dz4
  
  !r8 - dualz4
  elemental function restar8z4_dz4(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restar8z4_dz4
  
  !r4 - dualz4
  elemental function restar4z4_dz4(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restar4z4_dz4
  
  !i8 - dualz4
  elemental function restai8z4_dz4(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restai8z4_dz4

  !i4 - dualz4
  elemental function restai4z4_dz4(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = -A + B
  end function restai4z4_dz4  
  
  !+dualz4 (unary)
  elemental function mas_dz4(A) result(fres)
    type(dualz4), intent(in) :: A
    type(dualz4) :: fres

    fres = A
  end function mas_dz4

  !-dualz4 (unary)
  elemental function menos_dz4(A) result(f_res)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_res

    f_res%dualz3 = -A%dualz3
    f_res%f4 = -A%f4
  end function menos_dz4

  !Division
  !dualz4 / dualz4
  elemental function divz4z4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A, B
    type(dualz4) :: fr

    fr = A * invdz4_dz4(B)
  end function divz4z4_dz4

  !dualz4 / dualz3
  elemental function divz4z3_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz3), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invdz3_dz4(B)
  end function divz4z3_dz4

  !dualz4 / dualz2
  elemental function divz4z2_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invdz2_dz4(B)
  end function divz4z2_dz4

  !dualz4 / dualz1
  elemental function divz4z1_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invdz1_dz4(B)
  end function divz4z1_dz4

  !dualz4 / c8
  elemental function divz4c8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invc8_dz4(B)
  end function divz4c8_dz4

  !dualz4 / c4
  elemental function divz4c4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invc4_dz4(B)
  end function divz4c4_dz4

  !dualz4 / r8
  elemental function divz4r8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invr8_dz4(B)
  end function divz4r8_dz4

  !dualz4 / r4
  elemental function divz4r4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invr4_dz4(B)
  end function divz4r4_dz4

  !dualz4 / i8
  elemental function divz4i8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invi8_dz4(B)
  end function divz4i8_dz4

  !dualz4 / i4
  elemental function divz4i4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz4) :: fr

    fr = A * invi4_dz4(B)
  end function divz4i4_dz4

  !dualz3 / dualz4 
  elemental function divz3z4_dz4(B,A) result(fr)
    type(dualz3), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divz3z4_dz4

  !dualz2 / dualz4 
  elemental function divz2z4_dz4(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divz2z4_dz4

  !dualz1 / dualz4 
  elemental function divz1z4_dz4(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divz1z4_dz4

  !c8 / dualz4 
  elemental function divc8z4_dz4(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divc8z4_dz4

  !c4 / dualz4 
  elemental function divc4z4_dz4(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divc4z4_dz4

  !r8 / dualz4 
  elemental function divr8z4_dz4(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divr8z4_dz4

  !r4 / dualz4 
  elemental function divr4z4_dz4(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divr4z4_dz4
  
  !i8 / dualz4 
  elemental function divi8z4_dz4(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divi8z4_dz4

  !i4 / dualz4 
  elemental function divi4z4_dz4(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = B * invdz4_dz4(A)
  end function divi4z4_dz4
  
  !times
  !dualz4 * dualz4
  elemental function timesz4z4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A, B
    type(dualz4) :: fr
    complex(8) :: A0, A1, A2, A3, A4, B0, B1, B2, B3, B4

    A0 = A%f0; A1 = A%f1; A2 = A%f2; A3 = A%f3; A4 = A%f4
    B0 = B%f0; B1 = B%f1; B2 = B%f2; B3 = B%f3; B4 = B%f4
    
    fr%dualz3 = A%dualz3 * B%dualz3
    fr%f4 = A4*B0 + 4*A3*B1 + 6*A2*B2 + 4*A1*B3 + A0*B4
  end function timesz4z4_dz4

  !dualz4 * dualz3
  elemental function timesz4z3_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz3), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4z3_dz4

  !dualz4 * dualz2
  elemental function timesz4z2_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4z2_dz4

  !dualz4 * dualz1
  elemental function timesz4z1_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4z1_dz4

  !dualz4 * c8
  elemental function timesz4c8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4c8_dz4

  !dualz4 * c4
  elemental function timesz4c4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4c4_dz4

  !dualz4 * r8
  elemental function timesz4r8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4r8_dz4

  !dualz4 * r4
  elemental function timesz4r4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4r4_dz4

  !dualz4 * i8
  elemental function timesz4i8_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4i8_dz4

  !dualz4 * i4
  elemental function timesz4i4_dz4(A,B) result(fr)
    type(dualz4), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B
    fr = A*Bdz4
  end function timesz4i4_dz4  
  
  !dualz3 * dualz4
  elemental function timesz3z4_dz4(B,A) result(fr)
    type(dualz3), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesz3z4_dz4
  
  !dualz2 * dualz4
  elemental function timesz2z4_dz4(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesz2z4_dz4
 
  !dualz1 * dualz4
  elemental function timesz1z4_dz4(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesz1z4_dz4

  !c8 * dualz4
  elemental function timesc8z4_dz4(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesc8z4_dz4  

  !c4 * dualz4
  elemental function timesc4z4_dz4(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesc4z4_dz4

  !r8 * dualz4
  elemental function timesr8z4_dz4(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesr8z4_dz4
  
  !r4 * dualz4
  elemental function timesr4z4_dz4(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesr4z4_dz4

  !i8 * dualz4
  elemental function timesi8z4_dz4(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesi8z4_dz4

  !i4 * dualz4
  elemental function timesi4z4_dz4(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr

    fr = A*B
  end function timesi4z4_dz4
  
  !inverso multiplicativo
  elemental function invdz4_dz4(B) result(BI)
    type(dualz4), intent(in) :: B
    type(dualz4) :: BI
    complex(8) :: B0, B1, B2, B3, B4

    B0 = B%f0; B1 = B%f1; B2 = B%f2; B3 = B%f3; B4 = B%f4
    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2
    BI%f2 = (2*B1**2)/B0**3 - B2/B0**2
    BI%f3 = (-6*B1**3)/B0**4 + (6*B1*B2)/B0**3 - B3/B0**2
    BI%f4 = (24*B1**4)/B0**5 - (36*B1**2*B2)/B0**4 +                &
         (6*B2**2)/B0**3 + (8*B1*B3)/B0**3 - B4/B0**2
  end function invdz4_dz4

  elemental function invdz3_dz4(B) result(BI)
    type(dualz3), intent(in) :: B
    type(dualz4) :: BI
    complex(8) :: B0, B1, B2, B3

    B0 = B%f0; B1 = B%f1; B2 = B%f2; B3 = B%f3
    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2
    BI%f2 = (2*B1**2)/B0**3 - B2/B0**2
    BI%f3 = (-6*B1**3)/B0**4 + (6*B1*B2)/B0**3 - B3/B0**2
    BI%f4 = (24*B1**4)/B0**5 - (36*B1**2*B2)/B0**4 +                &
         (6*B2**2)/B0**3 + (8*B1*B3)/B0**3
  end function invdz3_dz4

  elemental function invdz2_dz4(B) result(BI)
    type(dualz2), intent(in) :: B
    type(dualz4) :: BI
    complex(8) :: B0, B1, B2

    B0 = B%f0; B1 = B%f1; B2 = B%f2
    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2
    BI%f2 = (2*B1**2)/B0**3 - B2/B0**2
    BI%f3 = (-6*B1**3)/B0**4 + (6*B1*B2)/B0**3
    BI%f4 = (24*B1**4)/B0**5 - (36*B1**2*B2)/B0**4 + (6*B2**2)/B0**3
  end function invdz2_dz4

  elemental function invdz1_dz4(B) result(BI)
    type(dualz1), intent(in) :: B
    type(dualz4) :: BI
    complex(8) :: B0, B1

    B0 = B%f0; B1 = B%f1
    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2
    BI%f2 = (2*B1**2)/B0**3
    BI%f3 = (-6*B1**3)/B0**4
    BI%f4 = (24*B1**4)/B0**5
  end function invdz1_dz4

  elemental function invc8_dz4(B) result(BI)
    complex(8), intent(in) :: B
    type(dualz4) :: BI

    BI%f0 = 1/B; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0; BI%f4 = 0
  end function invc8_dz4

  elemental function invc4_dz4(B) result(BI)
    complex(4), intent(in) :: B
    type(dualz4) :: BI

    BI%f0 = 1/B; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0; BI%f4 = 0
  end function invc4_dz4

  elemental function invr8_dz4(B) result(BI)
    real(8), intent(in) :: B
    type(dualz4) :: BI

    BI%f0 = 1/B; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0; BI%f4 = 0
  end function invr8_dz4

  elemental function invr4_dz4(B) result(BI)
    real(4), intent(in) :: B
    type(dualz4) :: BI

    BI%f0 = 1/B; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0; BI%f4 = 0
  end function invr4_dz4

  elemental function invi8_dz4(B) result(BI)
    integer(8), intent(in) :: B
    type(dualz4) :: BI

    BI%f0 = 1/B; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0; BI%f4 = 0
  end function invi8_dz4

  elemental function invi4_dz4(B) result(BI)
    integer(4), intent(in) :: B
    type(dualz4) :: BI

    BI%f0 = 1/B; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0; BI%f4 = 0
  end function invi4_dz4
  
  !Power
  !dualz4**dualz4
  !el caso A0 = 0 es problematico, solo se atendio este caso para cuando
  !1.- [A1,A2,A3,A4] = 0 y la parte real de B0 > 0 
  !2.- la parte real de B0 > 4
  !3.- [B1,B2,B3] = 0
  !en todos estos casos no hay restricciones para las demas componentes
  !es decir, el caso [B1,B2,B3,B4] = 0 no se requiere
  elemental function  powerz4z4_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A, B
    type(dualz4) :: fr
    complex(8) :: A0, A1, A2, A3, A4, B0, B1, B2, B3, B4

    A0 = A%f0; A1 = A%f1; A2 = A%f2; A3 = A%f3; A4 = A%f4
    B0 = B%f0; B1 = B%f1; B2 = B%f2; B3 = B%f3; B4 = B%f4

    fr%dualz3 = A%dualz3 ** B%dualz3
    if(all([B0,B1,B2,B3,B4]==0)) then
       fr%f4 = 0
    elseif(all([A0,A1,A2,A3,A4]==0) .and. real(B0) > 0) then
       fr%f4 = 0       
    elseif(all([A1,A2,A3,A4,B1,B2,B3,B4]==0)) then
       fr%f4 = 0
    elseif(A0 == 0 .and. real(B0)>4) then
       fr%f4 = 0
    elseif(all([A0,B1,B2,B3]==0)) then
       if(B0 == 1) then
          fr%f4 = A4
       elseif(B0 == 2) then
          fr%f4 = 6*A2**2 + 8*A1*A3
       elseif(B0 == 3) then
          fr%f4 = 36*A1**2*A2
       elseif(B0 == 4) then
          fr%f4 = 24*A1**4
       else
          fr%f4 = A0**(B0 - 4)*(-6*A1**4*B0 + 12*A0*A1**2*A2*B0 -      &
               3*A0**2*A2**2*B0 - 4*A0**2*A1*A3*B0 + A0**3*A4*B0 +     &
               8*A0*A1**3*B1 - 12*A0**2*A1*A2*B1 + 4*A0**3*A3*B1 -     &
               6*A0**2*A1**2*B2 + 6*A0**3*A2*B2 + 4*A0**3*A1*B3 +      &
               A0**4*B4*log(A0) + (A1*B0 + A0*B1*log(A0))**4 - 6*(A1*  &
               B0 + A0*B1*log(A0))**2*(A1**2*B0 - A0*A2*B0 - 2*A0*A1*  &
               B1 - A0**2*B2*log(A0)) + 3*(-(A1**2*B0)+ A0*A2*B0 + 2*  &
               A0*A1*B1 + A0**2*B2*log(A0))**2 + 4*(A1*B0 + A0*B1*     &
               log(A0))*(2*A1**3*B0 - 3*A0*A1**2*B1 + A0**2*(A3*B0 + 3*&
               A2*B1) + 3*A0*A1*(-(A2*B0) + A0*B2) + A0**3*B3*log(A0)))
       end if
    else
       fr%f4 = A0**(B0 - 4)*(-6*A1**4*B0 + 12*A0*A1**2*A2*B0 - 3*A0**2*&
            A2**2*B0 - 4*A0**2*A1*A3*B0 + A0**3*A4*B0 + 8*A0*A1**3*B1 -&
            12*A0**2*A1*A2*B1 + 4*A0**3*A3*B1 - 6*A0**2*A1**2*B2 + 6*  &
            A0**3*A2*B2 + 4*A0**3*A1*B3 + A0**4*B4*log(A0) + (A1*B0 +  &
            A0*B1*log(A0))**4 - 6*(A1*B0 + A0*B1*log(A0))**2*(A1**2*B0-&
            A0*A2*B0 - 2*A0*A1*B1 - A0**2*B2*log(A0)) + 3*(-(A1**2*B0)+&
            A0*A2*B0 + 2*A0*A1*B1 + A0**2*B2*log(A0))**2 + 4*(A1*B0 +  &
            A0*B1*log(A0))*(2*A1**3*B0 - 3*A0*A1**2*B1 + A0**2*(A3*B0 +&
            3*A2*B1) + 3*A0*A1*(-(A2*B0) + A0*B2) + A0**3*B3*log(A0)))
    end if
  end function powerz4z4_dz4

  !dualz4**dualz3
  elemental function  powerz4z3_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz3), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  !se debe hacer la conversion, no use exp(B*log(A)) aunque
    !B*log(A) ya esté definido para B dualz3. Si se pudiera usar
    !exp(Bdz4*log(A))
    
    fr = A**Bdz4
  end function powerz4z3_dz4

  !dualz4**dualz2
  elemental function  powerz4z2_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  
    fr = A**Bdz4
  end function powerz4z2_dz4

  !dualz4**dualz1
  elemental function  powerz4z1_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  
    fr = A**Bdz4
  end function powerz4z1_dz4

  !dualz4**c8
  elemental function  powerz4c8_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  
    fr = A**Bdz4
  end function powerz4c8_dz4
  
  !dualz4**c4
  elemental function  powerz4c4_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  
    fr = A**Bdz4
  end function powerz4c4_dz4
  
  !dualz4**r8
  elemental function  powerz4r8_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  
    fr = A**Bdz4
  end function powerz4r8_dz4

  !dualz4**r4
  elemental function  powerz4r4_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  
    fr = A**Bdz4
  end function powerz4r4_dz4

  !dualz4**i8
  elemental function  powerz4i8_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  
    fr = A**Bdz4
  end function powerz4i8_dz4

  !dualz4**i4
  elemental function  powerz4i4_dz4(A, B) result(fr)
    type(dualz4), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B  
    fr = A**Bdz4
  end function powerz4i4_dz4

  !dualz3**dualz4
  elemental function  powerz3z4_dz4(B, A) result(fr)
    type(dualz3), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function powerz3z4_dz4

  !dualz2**dualz4
  elemental function  powerz2z4_dz4(B, A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function powerz2z4_dz4

  !dualz1**dualz4
  elemental function  powerz1z4_dz4(B, A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function powerz1z4_dz4

  !c8**dualz4
  elemental function  powerc8z4_dz4(B, A) result(fr)
    complex(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function powerc8z4_dz4

  !c4**dualz4
  elemental function  powerc4z4_dz4(B, A) result(fr)
    complex(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function powerc4z4_dz4

  !r8**dualz4
  elemental function  powerr8z4_dz4(B, A) result(fr)
    real(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function powerr8z4_dz4

  !r4**dualz4
  elemental function  powerr4z4_dz4(B, A) result(fr)
    real(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function powerr4z4_dz4

  !i8**dualz4
  elemental function  poweri8z4_dz4(B, A) result(fr)
    integer(8), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function poweri8z4_dz4

  !i4**dualz4
  elemental function  poweri4z4_dz4(B, A) result(fr)
    integer(4), intent(in) :: B
    type(dualz4), intent(in) :: A
    type(dualz4) :: fr
    type(dualz4) :: Bdz4

    Bdz4 = B      
    fr = Bdz4**A
  end function poweri4z4_dz4
  
  !Mathematical Functions
  elemental function sin_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = sin(A%dualz3)
    f_r%f4 = -6*g1**2 * g2*cos(g0) + g4*cos(g0) + g1**4 * sin(g0) - &
         3*g2**2 * sin(g0) - 4*g1*g3*sin(g0)
  end function sin_dz4

  elemental function cos_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = cos(A%dualz3)
    f_r%f4 = g1**4*cos(g0) - 3*g2**2*cos(g0) - 4*g1*g3*cos(g0) +       &
         6*g1**2*g2*sin(g0) - g4*sin(g0)
  end function cos_dz4
  
  elemental function exp_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = exp(A%dualz3)
    f_r%f4 = exp(g0)*(g1**4 + 6*g1**2*g2 + 3*g2**2 + 4*g1*g3 + g4)    
  end function exp_dz4

  elemental function log_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = log(A%dualz3)
    f_r%f4 = (-6*g1**4)/g0**4 + (12*g1**2*g2)/g0**3 - (3*g2**2)/g0**2 -&
         (4*g1*g3)/g0**2 + g4/g0
  end function log_dz4

  elemental function tan_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4
    
    f_r%dualz3 = tan(A%dualz3)
    f_r%f4 = (g4 + (4*(3*g1**2*g2 + 4*g1**4*tan(g0)))/cos(g0)**2 + &
         2*tan(g0)*(3*g2**2 + 4*g1*g3 + 4*g1**2*Tan(g0)*(3*g2 +    &
         g1**2*tan(g0))))/cos(g0)**2
  end function tan_dz4

  elemental function sinh_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = sinh(A%dualz3)
    f_r%f4 = (6*g1**2*g2 + g4)*cosh(g0) + (g1**4 + 3*g2**2 + &
         4*g1*g3)*sinh(g0) 
  end function sinh_dz4

  elemental function cosh_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = cosh(A%dualz3)
    f_r%f4 =(g1**4 + 3*g2**2 + 4*g1*g3)*cosh(g0) + &
         (6*g1**2*g2 + g4)*sinh(g0)
  end function cosh_dz4

  elemental function tanh_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = tanh(A%dualz3)
    f_r%f4 = (g4 + 12*g1**2*g2*(2 - 3/cosh(g0)**2) + 2*(-3*g2**2 -     &
         4*g1*(g1**3 + g3) + (12*g1**4)/cosh(g0)**2)*tanh(g0))/        &
         cosh(g0)**2
  end function tanh_dz4

  elemental function asin_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = asin(A%dualz3)
    f_r%f4 = (6*g1**2*g2 + g0**3*(6*g1**4 - 6*g2**2 - 8*g1*g3) + g0**5*&
         (3*g2**2 + 4*g1*g3) + g0*(9*g1**4 + 3*g2**2 + 4*g1*g3) +      &
         g0**2*(6*g1**2*g2 - 3*g4) + g4 - g0**6*g4 + 3*g0**4*(-4*g1**2*&
         g2 + g4))/(1 - g0**2)**3.5
  end function asin_dz4

  elemental function acos_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = acos(A%dualz3)
    f_r%f4 = -((6*g1**2*g2 + g0**3*(6*g1**4 - 6*g2**2 - 8*g1*g3) +     &
         g0**5*(3*g2**2 + 4*g1*g3) + g0*(9*g1**4 + 3*g2**2 + 4*g1*g3) +&
         g0**2*(6*g1**2*g2 - 3*g4) + g4 - g0**6*g4 + 3*g0**4*(-4*g1**2*&
         g2 + g4))/(1 - g0**2)**3.5)
  end function acos_dz4

  elemental function atan_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = atan(A%dualz3)
    f_r%f4 = (-24*g0*(-1 + g0**2)*g1**4 + 12*(-1 + 2*g0**2 + 3*g0**4)* &
         g1**2*g2 - 6*g0*(1 + g0**2)**2*g2**2 - 8*g0*(1 + g0**2)**2*g1*&
         g3 + (1 + g0**2)**3*g4)/(1 + g0**2)**4
  end function atan_dz4

  elemental function acosh_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = acosh(A%dualz3)
    f_r%f4 = (-9*g0*g1**4 - 6*g0**3*g1**4 - 6*g1**2*g2 - 6*g0**2*g1**2*&
         g2 + 12*g0**4*g1**2*g2 - 3*g0*g2**2 + 6*g0**3*g2**2 - 3*g0**5*&
         g2**2 - 4*g0*g1*g3 + 8*g0**3*g1*g3 - 4*g0**5*g1*g3 + (-1 +    &
         g0**2)**3*g4)/((-1 + g0)**3.5*(1 + g0)**3.5)
  end function acosh_dz4

  elemental function asinh_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4
    
    f_r%dualz3 = asinh(A%dualz3)
    f_r%f4 = (-6*g1**2*g2 + g0*(9*g1**4 - 3*g2**2 - 4*g1*g3) - g0**5*  &
         (3*g2**2 + 4*g1*g3) - 2*g0**3*(3*(g1**4 + g2**2) + 4*g1*g3) + &
         g4 + g0**6*g4 + 3*g0**2*(2*g1**2*g2 + g4) + 3*g0**4*(4*g1**2* &
         g2 + g4))/(1 + g0**2)**3.5
  end function asinh_dz4

  elemental function atanh_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = atanh(A%dualz3)
    f_r%f4 = (6*(4*g0*(1 + g0**2)*g1**4 + 2*(1 + 2*g0**2 - 3*g0**4)*   &
         g1**2*g2 + g0*(-1 + g0**2)**2*g2**2) + 8*g0*(-1 + g0**2)**2*  &
         g1*g3 - (-1 + g0**2)**3*g4)/(-1 + g0**2)**4
  end function atanh_dz4

  elemental function atan2_dz4(A,B) result(f_r)
    type(dualz4), intent(in) :: A, B
    type(dualz4) :: f_r
    complex(8) :: A0, A1, A2, A3, A4, B0, B1, B2, B3, B4, r2

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2
    A3 = A%f3
    A4 = A%f4

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2
    B3 = B%f3
    B4 = B%f4

    r2 = A0**2 + B0**2
    f_r%dualz3 = atan2(A%dualz3, B%dualz3)
    f_r%f4 = (48*(-(A1*B0) + A0*B1)*(A0*A1 + B0*B1)**3 - 24*(A0*A1 +   &
         B0*B1)*(-(A1**3*B0) - 2*A0*A1*A2*B0 + A0*A1**2*B1 + A0**2*A2* &
         B1 - A2*B0**2*B1 + A0**2*A1*B2 - A1*B0*(B1**2 + B0*B2) +      &
         A0*B1*(B1**2 + 2*B0*B2))*r2 + 2*(A1**2*(-6*A2*B0 + 6*A0*B2) + &
         A0**2*(A3*B1 + 3*A2*B2) - 3*B0*(A3*B0*B1 + 2*A2*B1**2 + A2*B0*&
         B2) + A1*(-4*A0*A3*B0 + 3*A0**2*B3 - B0**2*B3) + A0*(-3*A2**2*&
         B0 + 6*B1**2*B2 + 3*B0*B2**2 + 4*B0*B1*B3))*r2**2 + (A4*B0 +  &
         2*A3*B1 - 2*A1*B3 - A0*B4)*r2**3)/r2**4
  end function atan2_dz4

  elemental function sqrt_dz4(A) result(f_r)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_r
    complex(8) :: g0, g1, g2, g3, g4

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3
    g4 = A%f4

    f_r%dualz3 = sqrt(A%dualz3)
    f_r%f4 = (-15*g1**4 + 36*g0*g1**2*g2 - 16*g0**2*g1*g3 +            &
         4*g0**2*(-3*g2**2 + 2*g0*g4))/(16d0*g0**3.5)
  end function sqrt_dz4

  !conjg
  !notice tat the conjugation operation is not, in general, differentia-
  !ble. The below definitions means (df)*, not d(f*), which, in some
  !times, is useful
  elemental function conjg_dz4(g) result(f_r)
    type(dualz4), intent(in) :: g
    type(dualz4) :: f_r

    f_r%dualz3 = conjg(g%dualz3)
    f_r%f4 = conjg(g%f4)
  end function conjg_dz4

  !abs
  !The absolute value is not, in general, differentiable. In the below
  !function (df)* is used instead of d(f*). For instance d(A * A*)
  !is taking as dA * A* + A * (dA)*. See also the conjg function
  elemental function abs_dz4(A) result(f_res)
    type(dualz4), intent(in) :: A
    type(dualz4) :: f_res

    f_res = sqrt(A*conjg(A))
  end function abs_dz4

  !|A> = |A0> + |A1> eps + |A2> eps2 + |A3>eps3 +  |A4>eps4 ; |B> = ...
  !As in the conjg and abs functions, here, the dual parts are (df)* not
  !d(f*)
  function dot_product_dz4(A,B) result(f_r)
    type(dualz4), intent(in), dimension(:) :: A, B
    type(dualz4) :: f_r
    complex(8), dimension(size(A)) :: A0, A1, A2, A3, A4
    complex(8), dimension(size(A)) :: B0, B1, B2, B3, B4

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2
    A3 = A%f3
    A4 = A%f4

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2
    B3 = B%f3
    B4 = B%f4

    f_r%dualz3 = dot_product(A%dualz3,B%dualz3)
    f_r%f4 = dot_product(A0,B4) + 4*dot_product(A1,B3) + &
         6*dot_product(A2,B2) + 4*dot_product(A3,B1) +   &
         dot_product(A4,B0)   
  end function dot_product_dz4

  !matmul
  function matmul22_dz4(A,B) result(f_r)
    type(dualz4), intent(in), dimension(:,:) :: A, B
    type(dualz4), dimension(size(A,1),size(B,2)) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A0, A1, A2, A3, A4
    complex(8), dimension(size(B,1),size(B,2)) :: B0, B1, B2, B3, B4

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2
    A3 = A%f3
    A4 = A%f4

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2
    B3 = B%f3
    B4 = B%f4

    f_r%dualz3 = matmul(A%dualz3,B%dualz3)
    f_r%f4 =  matmul(A0,B4) + 4*matmul(A1,B3) + 6*matmul(A2,B2) +      &
         4*matmul(A3,B1) + matmul(A4,B0)
  end function matmul22_dz4

  function matmul22_dz4C(A,B) result(fr)
    type(dualz4), intent(in), dimension(:,:) :: A
    complex(8), intent(in), dimension(:,:) :: B
    type(dualz4), dimension(size(A,1),size(B,2)) :: fr
    type(dualz4), dimension(size(B,1),size(B,2)) :: Bdz4

    Bdz4 = B
    fr = matmul22_dz4(A,Bdz4)
  end function matmul22_dz4C

  function matmul22_dz4R(A,B) result(fr)
    type(dualz4), intent(in), dimension(:,:) :: A
    real(8), intent(in), dimension(:,:) :: B
    type(dualz4), dimension(size(A,1),size(B,2)) :: fr
    type(dualz4), dimension(size(B,1),size(B,2)) :: Bdz4

    Bdz4 = B
    fr = matmul22_dz4(A,Bdz4)
  end function matmul22_dz4R

  function matmul22_Cdz4(A,B) result(fr)
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz4), intent(in), dimension(:,:) :: B
    type(dualz4), dimension(size(A,1),size(B,2)) :: fr
    type(dualz4), dimension(size(A,1),size(A,2)) :: Adz4

    Adz4 = A
    fr = matmul22_dz4(Adz4,B)
  end function matmul22_Cdz4
  
  function matmul22_Rdz4(A,B) result(fr)
    real(8), intent(in), dimension(:,:) :: A
    type(dualz4), intent(in), dimension(:,:) :: B
    type(dualz4), dimension(size(A,1),size(B,2)) :: fr
    type(dualz4), dimension(size(A,1),size(A,2)) :: Adz4

    Adz4 = A
    fr = matmul22_dz4(Adz4,B)
  end function matmul22_Rdz4
 
  !dual_rank2*dual_rank1
  function matmul_dr21_dz4(A,x) result(f_r)
    type(dualz4), intent(in), dimension(:,:) :: A
    type(dualz4), intent(in), dimension(size(A,2)) :: x
    type(dualz4), dimension(size(A,1)) :: f_r
    integer :: m, n

    m = size(A,1)
    n = size(A,2)
    f_r = reshape(matmul(A,reshape(x,[n,1])),[m])
  end function matmul_dr21_dz4

  function matmul_dr21_dz4C(A,x) result(fr)
    type(dualz4), intent(in), dimension(:,:) :: A
    complex(8), intent(in), dimension(size(A,2)) :: x
    type(dualz4), dimension(size(A,1)) :: fr
    type(dualz4), dimension(size(A,2)) :: xd

    xd = x
    fr = matmul_dr21_dz4(A,xd)
  end function matmul_dr21_dz4C

  function matmul_dr21_dz4R(A,x) result(fr)
    type(dualz4), intent(in), dimension(:,:) :: A
    real(8), intent(in), dimension(size(A,2)) :: x
    type(dualz4), dimension(size(A,1)) :: fr
    type(dualz4), dimension(size(A,2)) :: xd

    xd = x
    fr = matmul_dr21_dz4(A,xd)
  end function matmul_dr21_dz4R

  function matmul_dr21_Cdz4(A,x) result(fr)
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz4), intent(in), dimension(size(A,2)) :: x
    type(dualz4), dimension(size(A,1)) :: fr
    type(dualz4), dimension(size(A,1),size(A,2)) :: Adz4

    Adz4 = A
    fr = matmul_dr21_dz4(Adz4,x)
  end function matmul_dr21_Cdz4

  function matmul_dr21_Rdz4(A,x) result(fr)
    real(8), intent(in), dimension(:,:) :: A
    type(dualz4), intent(in), dimension(size(A,2)) :: x
    type(dualz4), dimension(size(A,1)) :: fr
    type(dualz4), dimension(size(A,1),size(A,2)) :: Adz4

    Adz4 = A
    fr = matmul_dr21_dz4(Adz4,x)
  end function matmul_dr21_Rdz4  

  !in order to avoid ambiguity we must change x to xx
  !"A generic function must be able to distinguish its arguments by 
  !type AND by name"
  function matmul_dr12_dz4(xx,A) result(f_r)
    type(dualz4), intent(in), dimension(:,:) :: A
    type(dualz4), intent(in), dimension(:) :: xx
    type(dualz4), dimension(size(A,2)) :: f_r
    integer :: m, n

    m = size(A,1)
    n = size(A,2)

    !if size(xx) /= m you may want to use some warning, the below code
    !relies on the reshape fortran function which uses padding
    !!$if(m  /= size(xx)) print*,'warning on matmul_dr12_dz4'
    f_r = reshape(matmul(reshape(xx,[1,m]),A),[n])
  end function matmul_dr12_dz4

  function matmul_dr12_dz4C(xx,A) result(fr)
    type(dualz4), intent(in), dimension(:) :: xx
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz4), dimension(size(A,2)) :: fr
    type(dualz4), dimension(size(A,1),size(A,2)) :: Adz4

    Adz4 = A
    fr =  matmul_dr12_dz4(xx,Adz4)
  end function matmul_dr12_dz4C
 
  function matmul_dr12_dz4R(xx,A) result(fr)
    type(dualz4), intent(in), dimension(:) :: xx
    real(8), intent(in), dimension(:,:) :: A
    type(dualz4), dimension(size(A,2)) :: fr
    type(dualz4), dimension(size(A,1),size(A,2)) :: Adz4

    Adz4 = A
    fr =  matmul_dr12_dz4(xx,Adz4)
  end function matmul_dr12_dz4R

  function matmul_dr12_Cdz4(xx,A) result(fr)
    complex(8), intent(in), dimension(:) :: xx
    type(dualz4), intent(in), dimension(:,:) :: A
    type(dualz4), dimension(size(A,2)) :: fr
    type(dualz4), dimension(size(xx)) :: xxdz4

    xxdz4 = xx
    fr =  matmul_dr12_dz4(xxdz4,A)
  end function matmul_dr12_Cdz4

  function matmul_dr12_Rdz4(xx,A) result(fr)
    real(8), intent(in), dimension(:) :: xx
    type(dualz4), intent(in), dimension(:,:) :: A
    type(dualz4), dimension(size(A,2)) :: fr
    type(dualz4), dimension(size(xx)) :: xxdz4

    xxdz4 = xx
    fr =  matmul_dr12_dz4(xxdz4,A)
  end function matmul_dr12_Rdz4
  
  !sum for rank2 arrays
  !the result is given as array of rank 1
  function sum2_dz4(A,k) result(f_r)
    type(dualz4), intent(in), dimension(:,:) :: A
    integer, intent(in) :: k
    type(dualz4), dimension(size(A,2/k)) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A4
    
    A4 = A%f4

    f_r%dualz3 = sum(A%dualz3,k)
    f_r%f4 = sum(A4,k)
  end function sum2_dz4

  !sum(A) rank2
  function sum0_dz4(A) result(f_r)
    type(dualz4), intent(in), dimension(:,:) :: A
    type(dualz4) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A4

    A4 = A%f4

    f_r%dualz3 = sum(A%dualz3)
    f_r%f4 = sum(A4)
  end function  sum0_dz4

  !sum(A) rank1
  function sum10_dz4(A) result(f_r)
    type(dualz4), intent(in), dimension(:) :: A
    type(dualz4) :: f_r    
    complex(8), dimension(size(A)) :: A4

    A4 = A%f4

    f_r%dualz3 = sum(A%dualz3)
    f_r%f4 = sum(A4)
  end function  sum10_dz4

  !product function for dualz1 vectors (arrays of rank 1)
  function prod1_dz4(x) result(f_r)
    type(dualz4), intent(in), dimension(:) :: x
    type(dualz4) :: f_r
    integer :: k

    f_r = 1
    do k = 1, size(x)
       f_r = f_r * x(k)
    end do
  end function prod1_dz4

  !product(x) for arrays of rank 2
  function prod0_dz4(x) result(f_r)
    type(dualz4), intent(in), dimension(:,:) :: x
    type(dualz4) :: f_r

    f_r = prod1_dz4(reshape(x,[size(x)]))
  end function prod0_dz4

  !product(A,k) for rank2 arrays
  !the result is given as array of rank 1
  function prod2_dz4(x,c) result(f_r)
    type(dualz4), intent(in), dimension(:,:) :: x
    integer, intent(in) :: c
    type(dualz4), dimension(size(x,2/c)) :: f_r
    type(dualz4), dimension(size(x,c)) :: xc
    integer :: k

    if(c==1) then
       do k = 1, size(x,2)
          xc = x(:,k)
          f_r(k) = prod1_dz4(xc)
       end do
    else if(c==2) then
       do k = 1, size(x,1)
          xc = x(k,:)
          f_r(k) = prod1_dz4(xc)
       end do
    else 
       stop 'c must be equal to 1 (2) to collapse rows (columns)'
    end if
  end function prod2_dz4

  !auxiliar functions
  elemental function f0part_dz4(A) result(fr)
    type(dualz4), intent(in) :: A
    complex(8) :: fr

    fr = A%f0
  end function f0part_dz4

  elemental function f1part_dz4(A) result(fr)
    type(dualz4), intent(in) :: A
    complex(8) :: fr
    
    fr = A%f1
  end function f1part_dz4

  elemental function f2part_dz4(A) result(fr)
    type(dualz4), intent(in) :: A
    complex(8) :: fr
    
    fr = A%f2
  end function f2part_dz4

  elemental function f3part_dz4(A) result(fr)
    type(dualz4), intent(in) :: A
    complex(8) :: fr
    
    fr = A%f3
  end function f3part_dz4

  elemental function f4part_dz4(A) result(fr)
    type(dualz4), intent(in) :: A
    complex(8) :: fr
    
    fr = A%f4
  end function f4part_dz4

!=====================================================================
  !auxiliar function to convert to dualz4
  !in principle this function is not necessary (the assignment operator
  !does the "same" job) but is provided for convenience
  !numeric (including dualz1, dualz2 and dualz3) to dualz4

  elemental function dz4todualz4(z) result(fr)
    type(dualz4), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function dz4todualz4
    
  elemental function dz3todualz4(z) result(fr)
    type(dualz3), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function dz3todualz4

  elemental function dz2todualz4(z) result(fr)
    type(dualz2), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function dz2todualz4

  elemental function dz1todualz4(z) result(fr)
    type(dualz1), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function dz1todualz4

  elemental function z8todualz4(z) result(fr)
    complex(8), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function z8todualz4

  !z4 here is for complex(4)
  elemental function z4todualz4(z) result(fr)
    complex(4), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function z4todualz4

  elemental function r8todualz4(z) result(fr)
    real(8), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function r8todualz4

  elemental function r4todualz4(z) result(fr)
    real(4), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function r4todualz4

  elemental function i8todualz4(z) result(fr)
    integer(8), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function i8todualz4

  elemental function i4todualz4(z) result(fr)
    integer(4), intent(in) :: z
    type(dualz4) :: fr

    fr = z
  end function i4todualz4
  !=====================================================================  

end module dualz4_mod


!!$  !matmul
!!$  !direct implementation of the matrix multiplication
!!$  !la funcion de arriba usa mas memoria pero es mas eficiente
!!$  !si se usa gfortran, con ifort da igual
!!$  function matmul22_dz3x(A,B) result(f_r)
!!$    type(dualz3), intent(in), dimension(:,:) :: A, B
!!$    type(dualz3), dimension(size(A,1),size(B,2)) :: f_r
!!$    integer :: m, n, p, i, j, k
!!$
!!$    m = size(A,1); n = size(A,2); p = size(B,2);
!!$
!!$    if(n /= size(B,1)) &
!!$         stop 'Error: arguments are inconsistent in [dual MATMUL]'
!!$
!!$    do i = 1, m
!!$       do j = 1, p
!!$          f_r(i,j) = 0
!!$          do k = 1, n
!!$             f_r(i,j) = f_r(i,j) + A(i,k) * B(k,j)
!!$          end do
!!$       end do
!!$    end do    
!!$  end function matmul22_dz3x

