!Dual numbers of complex components
!This can be used to compute 0-th, 1-th, 2-th and 3-th order derivatives
!F. Peñuñuri
!UADY, Merida Yucatan Mexico
!2022

module dualz3_mod
  use dualz1_mod
  use dualz2_mod
  implicit none
  private

  !dualz1 ---> [f0,f1]
  !dualz2 ---> [f0,f1,f2]
  !dualz3 ---> [f0,f1,f2,f3]
  type, public, extends (dualz2) :: dualz3
     complex(8) :: f3
  end type dualz3

  public :: assignment (=)

  public :: operator(==),  operator(/=), operator(+), operator(-)    
  public :: operator(*), operator (/), operator(**)   

  public :: sin, cos, exp, log, tan, sinh, cosh, tanh, asin, acos, atan
  public :: acosh, asinh, atanh, atan2, sqrt, conjg, abs, dot_product
  public :: matmul, sum, product

  public :: f0part_dz3, f1part_dz3, f2part_dz3, f3part_dz3
  public :: todualz3  

  !INTERFACES
  !equal assignment
  interface assignment (=)
     module procedure igualz2_dz3  !dualz3 <--- dualz2
     module procedure igualz1_dz3  !dualz3 <--- dualz1
     module procedure igualc8_dz3  !dualz3 <--- complex8
     module procedure igualc4_dz3  !dualz3 <--- complex4
     module procedure igualr8_dz3  !dualz3 <--- real8
     module procedure igualr4_dz3  !dualz3 <--- real4
     module procedure iguali8_dz3  !dualz3 <--- integer8
     module procedure iguali4_dz3  !dualz3 <--- integer4
  end interface assignment (=)

  !Logical equal operator 
  interface operator (==) 
     module procedure eq_dz3  !dualz3 == dualz3
  end interface operator (==)

  !Logical not equal operator 
  interface operator (/=) 
     module procedure noteq_dz3  !dualz3 /= dualz3
  end interface operator (/=)

  !Plus operator
  interface operator (+)
     module procedure mas_dz3      !unary 
     module procedure suma_dz3     !dualz3 + dualz3
     module procedure sumaz3z2_dz3 !dualz3 + dualz2
     module procedure sumaz3z1_dz3 !dualz3 + dualz1
     module procedure sumaz3c8_dz3 !dualz3 + complex8
     module procedure sumaz3c4_dz3 !dualz3 + complex4
     module procedure sumaz3r8_dz3 !dualz3 + real8
     module procedure sumaz3r4_dz3 !dualz3 + real4
     module procedure sumaz3i8_dz3 !dualz3 + integer8
     module procedure sumaz3i4_dz3 !dualz3 + integer4
     module procedure sumaz2z3_dz3 !dualz2 + dualz3
     module procedure sumaz1z3_dz3 !dualz1 + dualz3
     module procedure sumac8z3_dz3 !complex8 + dualz3
     module procedure sumac4z3_dz3 !complex4 + dualz3
     module procedure sumar8z3_dz3 !real8 + dualz3
     module procedure sumar4z3_dz3 !real4 + dualz3
     module procedure sumai8z3_dz3 !integer8 + dualz3
     module procedure sumai4z3_dz3 !integer4 + dualz3
  end interface operator (+)

  !Minus operator
  interface operator (-)
     module procedure menos_dz3     !unary 
     module procedure resta_dz3     !dualz3 - dualz3
     module procedure restaz3z2_dz3 !dualz3 - dualz2
     module procedure restaz3z1_dz3 !dualz3 - dualz1
     module procedure restaz3c8_dz3 !dualz3 - complex8
     module procedure restaz3c4_dz3 !dualz3 - complex4
     module procedure restaz3r8_dz3 !dualz3 - real8
     module procedure restaz3r4_dz3 !dualz3 - real4
     module procedure restaz3i8_dz3 !dualz3 - integer8
     module procedure restaz3i4_dz3 !dualz3 - integer4
     module procedure restaz2z3_dz3 !dualz2 - dualz3
     module procedure restaz1z3_dz3 !dualz1 - dualz3
     module procedure restac8z3_dz3 !complex8 - dualz3
     module procedure restac4z3_dz3 !complex4 - dualz3
     module procedure restar8z3_dz3 !real8 - dualz3
     module procedure restar4z3_dz3 !real4 - dualz3
     module procedure restai8z3_dz3 !integer8 - dualz3
     module procedure restai4z3_dz3 !integer4 - dualz3
  end interface operator (-)

  !Times operator
  interface operator (*)
     module procedure times_z3z3_dz3 !dualz3 * dualz3
     module procedure times_z3z2_dz3 !dualz3 * dualz2
     module procedure times_z3z1_dz3 !dualz3 * dualz1
     module procedure times_z3c8_dz3 !dualz3 * complex8
     module procedure times_z3c4_dz3 !dualz3 * complex4
     module procedure times_z3r8_dz3 !dualz3 * real8
     module procedure times_z3r4_dz3 !dualz3 * real4
     module procedure times_z3i8_dz3 !dualz3 * integer8
     module procedure times_z3i4_dz3 !dualz3 * integer4
     module procedure times_z2z3_dz3 !dualz2 * dualz3
     module procedure times_z1z3_dz3 !dualz1 * dualz3
     module procedure times_c8z3_dz3 !complex8 * dualz3
     module procedure times_c4z3_dz3 !complex4 * dualz3
     module procedure times_r8z3_dz3 !real8 * dualz3
     module procedure times_r4z3_dz3 !real4 * dualz3
     module procedure times_i8z3_dz3 !integer8 * dualz3
     module procedure times_i4z3_dz3 !integer4 * dualz3
  end interface operator (*)

  !Division operator
  interface operator (/)
     module procedure div_z3z3_dz3 !dualz3 / dualz3
     module procedure div_z3z2_dz3 !dualz3 / dualz2
     module procedure div_z3z1_dz3 !dualz3 / dualz1
     module procedure div_z3c8_dz3 !dualz3 / complex8
     module procedure div_z3c4_dz3 !dualz3 / complex4
     module procedure div_z3r8_dz3 !dualz3 / real8
     module procedure div_z3r4_dz3 !dualz3 / real4
     module procedure div_z3i8_dz3 !dualz3 / integer8
     module procedure div_z3i4_dz3 !dualz3 / integer4
     module procedure div_z2z3_dz3 !dualz2 / dualz3
     module procedure div_z1z3_dz3 !dualz1 / dualz3
     module procedure div_c8z3_dz3 !complex8 / dualz3
     module procedure div_c4z3_dz3 !complex4 / dualz3
     module procedure div_r8z3_dz3 !real8 / dualz3
     module procedure div_r4z3_dz3 !real4 / dualz3
     module procedure div_i8z3_dz3 !integer8 / dualz3
     module procedure div_i4z3_dz3 !integer4 / dualz3
  end interface operator (/)

  !Power operator
  interface operator (**)
     module procedure power_z3z3_dz3 !dualz3 / dualz3
     module procedure power_z3z2_dz3 !dualz3 / dualz2
     module procedure power_z3z1_dz3 !dualz3 / dualz1
     module procedure power_z3c8_dz3 !dualz3 / complex8
     module procedure power_z3c4_dz3 !dualz3 / complex4
     module procedure power_z3r8_dz3 !dualz3 / real8
     module procedure power_z3r4_dz3 !dualz3 / real4
     module procedure power_z3i8_dz3 !dualz3 / integer8
     module procedure power_z3i4_dz3 !dualz3 / integer4
     module procedure power_z2z3_dz3 !dualz2 / dualz3
     module procedure power_z1z3_dz3 !dualz1 / dualz3
     module procedure power_c8z3_dz3 !complex8 / dualz3
     module procedure power_c4z3_dz3 !complex4 / dualz3
     module procedure power_r8z3_dz3 !real8 / dualz3
     module procedure power_r4z3_dz3 !real4 / dualz3
     module procedure power_i8z3_dz3 !integer8 / dualz3
     module procedure power_i4z3_dz3 !integer4 / dualz3
  end interface operator (**)

  !Mathematical Functions
  interface sin
     module procedure sin_dz3
  end interface sin

  interface cos
     module procedure cos_dz3
  end interface cos

  interface exp
     module procedure exp_dz3
  end interface exp

  interface log
     module procedure log_dz3
  end interface log

  interface tan
     module procedure tan_dz3
  end interface tan

  interface sinh
     module procedure sinh_dz3
  end interface sinh

  interface cosh
     module procedure cosh_dz3
  end interface cosh

  interface tanh
     module procedure tanh_dz3
  end interface tanh

  interface asin
     module procedure asin_dz3
  end interface asin

  interface acos
     module procedure acos_dz3
  end interface acos

  interface atan
     module procedure atan_dz3
  end interface atan

  interface acosh
     module procedure acosh_dz3
  end interface acosh

  interface asinh
     module procedure asinh_dz3
  end interface asinh

  interface atanh
     module procedure atanh_dz3
  end interface atanh

  interface atan2
     module procedure atan2_dz3
  end interface atan2

  interface sqrt
     module procedure sqrt_dz3
  end interface sqrt

  interface conjg
     module procedure conjg_dz3
  end interface conjg

  interface abs
     module procedure abs_dz3
  end interface abs

  !<A|B>
  interface dot_product
     module procedure dot_product_dz3
  end interface dot_product

  !Matrix multiplication
  interface matmul
     module procedure matmul22_dz3    !dual_rank2*dual_rank2
     module procedure matmul_dr21_dz3 !dual_rank2*dual_rank1
     module procedure matmul_dr12_dz3 !dual_rank1*dual_rank2
     !
     module procedure matmul22_dz3C
     module procedure matmul22_dz3R
     module procedure matmul22_Cdz3
     module procedure matmul22_Rdz3
     !
     module procedure matmul_dr21_dz3C
     module procedure matmul_dr21_dz3R
     module procedure matmul_dr21_Cdz3
     module procedure matmul_dr21_Rdz3
     !
     module procedure matmul_dr12_dz3C
     module procedure matmul_dr12_dz3R
     module procedure matmul_dr12_Cdz3
     module procedure matmul_dr12_Rdz3
  end interface matmul

  interface sum
     module procedure sum2_dz3  !sum(dual_rank2,k)
     module procedure sum0_dz3  !sum(dual_rank2)
     module procedure sum10_dz3 !sum(dual_rank1)
  end interface sum

  !product for dual vectors x(k)
  interface product
     module procedure prod1_dz3
     module procedure prod0_dz3
     module procedure prod2_dz3 !product(dual_rank2,k)
  end interface product


  !auxiliar function  
  interface todualz3
     module procedure dz3todualz3
     module procedure dz2todualz3
     module procedure dz1todualz3
     module procedure z8todualz3
     module procedure z4todualz3
     module procedure r8todualz3
     module procedure r4todualz3
     module procedure i8todualz3
     module procedure i4todualz3
  end interface todualz3
  !terminan interfaces

  !Functions and subroutines associated to the above interfaces
contains  
  !Assignment, equal operator
  !dualz3 <--- dualz2
  elemental subroutine igualz2_dz3(A, z)
    type(dualz3), intent(out) :: A
    type(dualz2), intent(in) :: z

    A%dualz2 = z
    A%f3 = 0
  end subroutine igualz2_dz3

  !dualz3 <--- dualz1
  elemental subroutine igualz1_dz3(A, z)
    type(dualz3), intent(out) :: A
    type(dualz1), intent(in) :: z

    A%dualz1 = z
    A%f2 = 0
    A%f3 = 0
  end subroutine igualz1_dz3

  !dualz3 <--- complex8
  elemental subroutine igualc8_dz3(A, z)
    type(dualz3), intent(out) :: A
    complex(8), intent(in) :: z

    A%f0 = z
    A%f1 = 0
    A%f2 = 0
    A%f3 = 0
  end subroutine igualc8_dz3

  !dualz3 <--- complex4
  elemental subroutine igualc4_dz3(A, z)
    type(dualz3), intent(out) :: A
    complex(4), intent(in) :: z

    A%f0 = z
    A%f1 = 0
    A%f2 = 0
    A%f3 = 0
  end subroutine igualc4_dz3

  !dualz3 <--- real8
  elemental subroutine igualr8_dz3(A, z)
    type(dualz3), intent(out) :: A
    real(8), intent(in) :: z

    A%f0 = z
    A%f1 = 0
    A%f2 = 0
    A%f3 = 0
  end subroutine igualr8_dz3

  !dualz3 <--- real4
  elemental subroutine igualr4_dz3(A, z)
    type(dualz3), intent(out) :: A
    real(4), intent(in) :: z

    A%f0 = z
    A%f1 = 0
    A%f2 = 0
    A%f3 = 0
  end subroutine igualr4_dz3

  !dualz3 <--- integer8
  elemental subroutine iguali8_dz3(A, z)
    type(dualz3), intent(out) :: A
    integer(8), intent(in) :: z

    A%f0 = z
    A%f1 = 0
    A%f2 = 0
    A%f3 = 0
  end subroutine iguali8_dz3

  !dualz3 <--- integer4
  elemental subroutine iguali4_dz3(A, z)
    type(dualz3), intent(out) :: A
    integer(4), intent(in) :: z

    A%f0 = z
    A%f1 = 0
    A%f2 = 0
    A%f3 = 0
  end subroutine iguali4_dz3

  !Logical not equal operator
  elemental function noteq_dz3(lhs, rhs) result(f_res)
    type (dualz3), intent(in) :: lhs, rhs
    logical :: f_res

    f_res = .not.(lhs == rhs)
  end function noteq_dz3

  !Logical equal operator
  elemental function eq_dz3(lhs, rhs) result(f_res)
    type (dualz3), intent(in) :: lhs, rhs
    logical :: f_res
    logical :: eqf0, eqf1, eqf2, eqf3

    eqf0 = lhs%f0 == rhs%f0
    eqf1 = lhs%f1 == rhs%f1
    eqf2 = lhs%f2 == rhs%f2
    eqf3 = lhs%f3 == rhs%f3

    f_res = all([eqf0,eqf1,eqf2,eqf3])
  end function eq_dz3

  !dualz3 + dualz3
  elemental function suma_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A, B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B%dualz2
    fr%f3 = A%f3 + B%f3
  end function suma_dz3

  !dualz3 + dualz2
  elemental function sumaz3z2_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B
    fr%f3 = A%f3
  end function sumaz3z2_dz3

  !dualz3 + dualz1
  elemental function sumaz3z1_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B
    fr%f3 = A%f3
  end function sumaz3z1_dz3

  !dualz3 + c8
  elemental function sumaz3c8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B
    fr%f3 = A%f3
  end function sumaz3c8_dz3

  !dualz3 + c4
  elemental function sumaz3c4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B
    fr%f3 = A%f3
  end function sumaz3c4_dz3

  !dualz3 + r8
  elemental function sumaz3r8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B
    fr%f3 = A%f3
  end function sumaz3r8_dz3

  !dualz3 + r4
  elemental function sumaz3r4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B
    fr%f3 = A%f3
  end function sumaz3r4_dz3

  !dualz3 + i8
  elemental function sumaz3i8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B
    fr%f3 = A%f3
  end function sumaz3i8_dz3

  !dualz3 + i4
  elemental function sumaz3i4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 + B
    fr%f3 = A%f3
  end function sumaz3i4_dz3

  !dualz2 + dualz3
  elemental function sumaz2z3_dz3(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = sumaz3z2_dz3(A,B) 
  end function sumaz2z3_dz3

  !dualz1 + dualz3
  elemental function sumaz1z3_dz3(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = sumaz3z1_dz3(A,B) 
  end function sumaz1z3_dz3

  !c8 + dualz3
  elemental function sumac8z3_dz3(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = sumaz3c8_dz3(A,B) 
  end function sumac8z3_dz3

  !c4 + dualz3
  elemental function sumac4z3_dz3(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = sumaz3c4_dz3(A,B) 
  end function sumac4z3_dz3

  !r8 + dualz3
  elemental function sumar8z3_dz3(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = sumaz3r8_dz3(A,B) 
  end function sumar8z3_dz3

  !r4 + dualz3
  elemental function sumar4z3_dz3(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = sumaz3r4_dz3(A,B) 
  end function sumar4z3_dz3

  !i8 + dualz3
  elemental function sumai8z3_dz3(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = sumaz3i8_dz3(A,B) 
  end function sumai8z3_dz3

  !i4 + dualz3
  elemental function sumai4z3_dz3(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = sumaz3i4_dz3(A,B) 
  end function sumai4z3_dz3

  !+dualz3 (unary)
  elemental function mas_dz3(A) result(f_res)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_res

    f_res = A
  end function mas_dz3

  !dualz3 - dualz3
  elemental function resta_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A, B
    type(dualz3) :: fr

    fr = A + (-B)
  end function resta_dz3

  !dualz3 - dualz2
  elemental function restaz3z2_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz3) :: fr

    fr = A + (-B)
  end function restaz3z2_dz3

  !dualz3 - dualz1
  elemental function restaz3z1_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz3) :: fr

    fr = A + (-B)
  end function restaz3z1_dz3

  !dualz3 - c8
  elemental function restaz3c8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz3) :: fr

    fr = A + (-B)
  end function restaz3c8_dz3

  !dualz3 - c4
  elemental function restaz3c4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz3) :: fr

    fr = A + (-B)
  end function restaz3c4_dz3

  !dualz3 - r8
  elemental function restaz3r8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz3) :: fr

    fr = A + (-B)
  end function restaz3r8_dz3

  !dualz3 - r4
  elemental function restaz3r4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz3) :: fr

    fr = A + (-B)
  end function restaz3r4_dz3

  !dualz3 - i8
  elemental function restaz3i8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz3) :: fr

    fr = A + (-B)
  end function restaz3i8_dz3

  !dualz3 - i4
  elemental function restaz3i4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz3) :: fr

    fr = A + (-B)
  end function restaz3i4_dz3

  !dualz2 - dualz3
  elemental function restaz2z3_dz3(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = -(A - B)
  end function restaz2z3_dz3

  !dualz1 - dualz3
  elemental function restaz1z3_dz3(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = -(A - B)
  end function restaz1z3_dz3

  !c8 - dualz3
  elemental function restac8z3_dz3(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = -(A - B)
  end function restac8z3_dz3

  !c4 - dualz3
  elemental function restac4z3_dz3(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = -(A - B)
  end function restac4z3_dz3

  !r8 - dualz3
  elemental function restar8z3_dz3(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = -(A - B)
  end function restar8z3_dz3

  !r4 - dualz3
  elemental function restar4z3_dz3(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = -(A - B)
  end function restar4z3_dz3

  !i8 - dualz3
  elemental function restai8z3_dz3(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = -(A - B)
  end function restai8z3_dz3

  !i4 - dualz3
  elemental function restai4z3_dz3(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = -(A - B)
  end function restai4z3_dz3

  !-dualz3 (unary)
  elemental function menos_dz3(A) result(f_res)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_res

    f_res%dualz2 = -A%dualz2
    f_res%f3 = -A%f3
  end function menos_dz3

  !Power
  !z3**z3
  !el caso A0 = 0 es problematico, solo se atendio este caso para cuando
  !1.- [A1,A2,A3] = 0 y la parte real de B0 > 0 
  !2.- la parte real de B0 > 3
  !3.- [B1,B2] = 0
  !en todos estos casos no hay restricciones para las demas componentes
  !es decir, el caso [B1,B2,B3] = 0 no se requiere
  elemental function power_z3z3_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A, B
    type(dualz3) :: fr
    complex(8) :: A0, A1, A2, A3, B0, B1, B2, B3

    A0 = A%f0; A1 = A%f1; A2 = A%f2; A3 = A%f3
    B0 = B%f0; B1 = B%f1; B2 = B%f2; B3 = B%f3

    fr%dualz2 = A%dualz2**B%dualz2

    if(all([B0,B1,B2,B3]==0)) then
       fr%f3 = 0
    elseif(all([A0,A1,A2,A3]==0) .and. real(B0) > 0) then
       fr%f3 = 0    
    elseif(all([A1,A2,A3,B1,B2,B3]==0)) then
       fr%f3 = 0
    elseif(A0 == 0 .and. real(B0)>3) then
       fr%f3 = 0
    elseif(all([A0,B1,B2]==0)) then
       if(B0 == 1) then
          fr%f3 = A3
       elseif(B0 == 2) then
          fr%f3 = 6*A1*A2
       elseif(B0 == 3) then
          fr%f3 = 6*A1**3
       else
          fr%f3 = A0**(B0 - 3)*(2*A1**3*B0 - 3*A0*A1*A2*B0 + A0**2*A3* &
               B0 - 3*A0*A1**2*B1 + 3*A0**2*A2*B1 + 3*A0**2*A1*B2 +    &
               A0**3*B3*log(A0) + (A1*B0 + A0*B1*log(A0))**3 - 3*(A1*  &
               B0 + A0*B1*log(A0))*(A1**2*B0 - A0*A2*B0 - 2*A0*A1*B1 - &
               A0**2*B2*log(A0)))
       end if
    else
       fr%f3 = A0**(B0 - 3)*(2*A1**3*B0 - 3*A0*A1*A2*B0 + A0**2*A3*B0 -&
            3*A0*A1**2*B1 + 3*A0**2*A2*B1 + 3*A0**2*A1*B2 + A0**3*B3*  &
            log(A0) + (A1*B0 + A0*B1*log(A0))**3 - 3*(A1*B0 + A0*B1*   &
            log(A0))*(A1**2*B0 - A0*A2*B0 - 2*A0*A1*B1 - A0**2*B2*     &
            log(A0)))
    end if
  end function power_z3z3_dz3

  !z3**z2
  elemental function power_z3z2_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B !se debe hacer la conversion, no use exp(B*log(A)) aunque
    !B*log(A) ya esté definido para B dualz2. Si se pudiera usar
    !exp(Bz3*log(A))

    fr = power_z3z3_dz3(A,Bz3)
  end function power_z3z2_dz3

  !z3**z1
  elemental function power_z3z1_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = power_z3z3_dz3(A,Bz3)
  end function power_z3z1_dz3

  !z3**c8
  elemental function power_z3c8_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = power_z3z3_dz3(A,Bz3)
  end function power_z3c8_dz3

  !z3**c4
  elemental function power_z3c4_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = power_z3z3_dz3(A,Bz3)
  end function power_z3c4_dz3

  !z3**r8
  elemental function power_z3r8_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = power_z3z3_dz3(A,Bz3)
  end function power_z3r8_dz3

  !z3**r4
  elemental function power_z3r4_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = power_z3z3_dz3(A,Bz3)
  end function power_z3r4_dz3

  !z3**i8
  elemental function power_z3i8_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = power_z3z3_dz3(A,Bz3)
  end function power_z3i8_dz3

  !z3**i4
  elemental function power_z3i4_dz3(A, B) result(fr)
    type(dualz3), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = power_z3z3_dz3(A,Bz3)
  end function power_z3i4_dz3

  !----****----
  !z2**z3
  elemental function power_z2z3_dz3(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B

    fr = power_z3z3_dz3(Bz3,A) 
  end function power_z2z3_dz3

  !z1**z3
  elemental function power_z1z3_dz3(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B

    fr = power_z3z3_dz3(Bz3,A) 
  end function power_z1z3_dz3

  !c8**z3
  elemental function power_c8z3_dz3(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B

    fr = power_z3z3_dz3(Bz3,A) 
  end function power_c8z3_dz3

  !c4**z3
  elemental function power_c4z3_dz3(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B

    fr = power_z3z3_dz3(Bz3,A) 
  end function power_c4z3_dz3

  !r8**z3
  elemental function power_r8z3_dz3(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B

    fr = power_z3z3_dz3(Bz3,A) 
  end function power_r8z3_dz3

  !r4**z3
  elemental function power_r4z3_dz3(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B

    fr = power_z3z3_dz3(Bz3,A) 
  end function power_r4z3_dz3

  !i8**z3
  elemental function power_i8z3_dz3(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B

    fr = power_z3z3_dz3(Bz3,A) 
  end function power_i8z3_dz3

  !i4**z3
  elemental function power_i4z3_dz3(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B

    fr = power_z3z3_dz3(Bz3,A) 
  end function power_i4z3_dz3

  !times
  !dualz3 * dualz3
  elemental function times_z3z3_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A, B
    type(dualz3) :: fr
    complex(8) :: A0, A1, A2, A3

    A0 = A%f0; A1 = A%f1; A2 = A%f2; A3 = A%f3

    fr%dualz2 = A%dualz2 * B%dualz2
    fr%f3 = A3*(B%f0) + 3*A2*(B%f1) + 3*A1*(B%f2) + A0*(B%f3)
  end function times_z3z3_dz3

  !dualz3 * dualz2
  elemental function times_z3z2_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = A*Bz3
  end function times_z3z2_dz3

  !dualz3 * dualz1
  elemental function times_z3z1_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz3) :: fr
    type(dualz3) :: Bz3

    Bz3 = B
    fr = A*Bz3
  end function times_z3z1_dz3

  !dualz3 * c8
  elemental function times_z3c8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 * B
    fr%f3 = (A%f3)*B
  end function times_z3c8_dz3

  !dualz3 * c4
  elemental function times_z3c4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 * B
    fr%f3 = (A%f3)*B
  end function times_z3c4_dz3

  !dualz3 * r8
  elemental function times_z3r8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 * B
    fr%f3 = (A%f3)*B
  end function times_z3r8_dz3

  !dualz3 * r4
  elemental function times_z3r4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 * B
    fr%f3 = (A%f3)*B
  end function times_z3r4_dz3

  !dualz3 * i8
  elemental function times_z3i8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 * B
    fr%f3 = (A%f3)*B
  end function times_z3i8_dz3

  !dualz3 * i4
  elemental function times_z3i4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz3) :: fr

    fr%dualz2 = A%dualz2 * B
    fr%f3 = (A%f3)*B
  end function times_z3i4_dz3

  !dualz2 * dualz3
  elemental function times_z2z3_dz3(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz3), intent(in) :: A    
    type(dualz3) :: fr

    fr = times_z3z2_dz3(A,B)
  end function times_z2z3_dz3

  !dualz1 * dualz3
  elemental function times_z1z3_dz3(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz3), intent(in) :: A    
    type(dualz3) :: fr

    fr = times_z3z1_dz3(A,B)
  end function times_z1z3_dz3

  !c8 * dualz3
  elemental function times_c8z3_dz3(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz3), intent(in) :: A    
    type(dualz3) :: fr

    fr = times_z3c8_dz3(A,B)
  end function times_c8z3_dz3

  !c4 * dualz3
  elemental function times_c4z3_dz3(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz3), intent(in) :: A    
    type(dualz3) :: fr

    fr = times_z3c4_dz3(A,B)
  end function times_c4z3_dz3

  !r8 * dualz3
  elemental function times_r8z3_dz3(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz3), intent(in) :: A    
    type(dualz3) :: fr

    fr = times_z3r8_dz3(A,B)
  end function times_r8z3_dz3

  !r4 * dualz3
  elemental function times_r4z3_dz3(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz3), intent(in) :: A    
    type(dualz3) :: fr

    fr = times_z3r4_dz3(A,B)
  end function times_r4z3_dz3

  !i8 * dualz3
  elemental function times_i8z3_dz3(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz3), intent(in) :: A    
    type(dualz3) :: fr

    fr = times_z3i8_dz3(A,B)
  end function times_i8z3_dz3

  !i4 * dualz3
  elemental function times_i4z3_dz3(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz3), intent(in) :: A    
    type(dualz3) :: fr

    fr = times_z3i4_dz3(A,B)
  end function times_i4z3_dz3

  !division
  !dualz3 / dualz3
  elemental function div_z3z3_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A, B
    type(dualz3) :: fr

    fr = A * inv_dz3_dz3(B)
  end function div_z3z3_dz3

  !dualz3 / dualz2
  elemental function div_z3z2_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz2), intent(in) :: B
    type(dualz3) :: fr

    fr = A * inv_dz2_dz3(B)
  end function div_z3z2_dz3

  !dualz3 / dualz1
  elemental function div_z3z1_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz3) :: fr

    fr = A * inv_dz1_dz3(B)
  end function div_z3z1_dz3

  !dualz3 / c8
  elemental function div_z3c8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    complex(8), intent(in) :: B
    type(dualz3) :: fr

    fr = A * inv_c8_dz3(B)
  end function div_z3c8_dz3

  !dualz3 / c4
  elemental function div_z3c4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    complex(4), intent(in) :: B
    type(dualz3) :: fr

    fr = A * inv_c4_dz3(B)
  end function div_z3c4_dz3

  !dualz3 / r8
  elemental function div_z3r8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    real(8), intent(in) :: B
    type(dualz3) :: fr

    fr = A * inv_r8_dz3(B)
  end function div_z3r8_dz3

  !dualz3 / r4
  elemental function div_z3r4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    real(4), intent(in) :: B
    type(dualz3) :: fr

    fr = A * inv_r4_dz3(B)
  end function div_z3r4_dz3

  !dualz3 / i8
  elemental function div_z3i8_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    integer(8), intent(in) :: B
    type(dualz3) :: fr

    fr = A * inv_i8_dz3(B)
  end function div_z3i8_dz3

  !dualz3 / i4
  elemental function div_z3i4_dz3(A,B) result(fr)
    type(dualz3), intent(in) :: A
    integer(4), intent(in) :: B
    type(dualz3) :: fr

    fr = A * inv_i4_dz3(B)
  end function div_z3i4_dz3

  !dualz2 / dualz3
  elemental function div_z2z3_dz3(B,A) result(fr)
    type(dualz2), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = B * inv_dz3_dz3(A)
  end function div_z2z3_dz3

  !dualz1 / dualz3
  elemental function div_z1z3_dz3(B,A) result(fr)
    type(dualz1), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = B * inv_dz3_dz3(A)
  end function div_z1z3_dz3

  !c8 / dualz3
  elemental function div_c8z3_dz3(B,A) result(fr)
    complex(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = B * inv_dz3_dz3(A)
  end function div_c8z3_dz3

  !c4 / dualz3
  elemental function div_c4z3_dz3(B,A) result(fr)
    complex(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = B * inv_dz3_dz3(A)
  end function div_c4z3_dz3

  !r8 / dualz3
  elemental function div_r8z3_dz3(B,A) result(fr)
    real(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = B * inv_dz3_dz3(A)
  end function div_r8z3_dz3

  !r4 / dualz3
  elemental function div_r4z3_dz3(B,A) result(fr)
    real(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = B * inv_dz3_dz3(A)
  end function div_r4z3_dz3

  !i8 / dualz3
  elemental function div_i8z3_dz3(B,A) result(fr)
    integer(8), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = B * inv_dz3_dz3(A)
  end function div_i8z3_dz3

  !i4 / dualz3
  elemental function div_i4z3_dz3(B,A) result(fr)
    integer(4), intent(in) :: B
    type(dualz3), intent(in) :: A
    type(dualz3) :: fr

    fr = B * inv_dz3_dz3(A)
  end function div_i4z3_dz3

  !=====================================================================
  !auxiliar functions for inverse
  !inverso multiplicativo
  elemental function inv_dz3_dz3(B) result(BI)
    type(dualz3), intent(in) :: B  
    type(dualz3) :: BI
    complex(8) :: B0, B1, B2, B3

    B0 = B%f0; B1 = B%f1; B2 = B%f2; B3 = B%f3

    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2
    BI%f2 = (2*B1**2)/B0**3 - B2/B0**2
    BI%f3 = (-6*B1**3)/B0**4 + (6*B1*B2)/B0**3 - B3/B0**2
  end function inv_dz3_dz3

  elemental function inv_dz2_dz3(B) result(BI)
    type(dualz2), intent(in) :: B  
    type(dualz3) :: BI
    complex(8) :: B0, B1, B2

    B0 = B%f0; B1 = B%f1; B2 = B%f2

    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2
    BI%f2 = (2*B1**2)/B0**3 - B2/B0**2
    BI%f3 = (-6*B1**3)/B0**4 + (6*B1*B2)/B0**3
  end function inv_dz2_dz3

  elemental function inv_dz1_dz3(B) result(BI)
    type(dualz1), intent(in) :: B  
    type(dualz3) :: BI
    complex(8) :: B0, B1

    B0 = B%f0; B1 = B%f1

    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2
    BI%f2 = (2*B1**2)/B0**3
    BI%f3 = (-6*B1**3)/B0**4
  end function inv_dz1_dz3

  elemental function inv_c8_dz3(B0) result(BI)
    complex(8), intent(in) :: B0 
    type(dualz3) :: BI

    BI%f0 = 1/B0; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0
  end function inv_c8_dz3

  elemental function inv_c4_dz3(B0) result(BI)
    complex(4), intent(in) :: B0 
    type(dualz3) :: BI

    BI%f0 = 1/B0; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0
  end function inv_c4_dz3

  elemental function inv_r8_dz3(B0) result(BI)
    real(8), intent(in) :: B0 
    type(dualz3) :: BI

    BI%f0 = 1/B0; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0
  end function inv_r8_dz3

  elemental function inv_r4_dz3(B0) result(BI)
    real(4), intent(in) :: B0 
    type(dualz3) :: BI

    BI%f0 = 1/B0; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0
  end function inv_r4_dz3

  elemental function inv_i8_dz3(B0) result(BI)
    integer(8), intent(in) :: B0 
    type(dualz3) :: BI

    BI%f0 = 1/B0; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0
  end function inv_i8_dz3

  elemental function inv_i4_dz3(B0) result(BI)
    integer(4), intent(in) :: B0 
    type(dualz3) :: BI

    BI%f0 = 1/B0; BI%f1 = 0; BI%f2 = 0; BI%f3 = 0
  end function inv_i4_dz3
  !=====================================================================

  !Mathematical Functions
  elemental function sin_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = sin(A%dualz2)
    f_r%f3 = -cos(g0)*g1**3 + g3*cos(g0) - 3*g1*g2*sin(g0)
  end function sin_dz3

  elemental function cos_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = cos(A%dualz2)
    f_r%f3 = -3*g1*g2*cos(g0) + g1**3*sin(g0) - g3*sin(g0)
  end function cos_dz3

  elemental function exp_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = exp(A%dualz2)
    f_r%f3 = exp(g0)*g1**3 + 3*exp(g0)*g1*g2 + exp(g0)*g3
  end function exp_dz3

  elemental function log_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = log(A%dualz2)
    f_r%f3 = (2*g1**3)/g0**3 - (3*g1*g2)/g0**2 + g3/g0    
  end function log_dz3

  elemental function tan_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = tan(A%dualz2)
    f_r%f3 = g3/cos(g0)**2 + 2*g1**3/cos(g0)**4 +                      &
         6*g1*g2/cos(g0)**2*tan(g0) + 4*g1**3/cos(g0)**2*tan(g0)**2
  end function tan_dz3

  elemental function sinh_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = sinh(A%dualz2)
    f_r%f3 = g1**3*cosh(g0) + g3*cosh(g0) + 3*g1*g2*sinh(g0)
  end function sinh_dz3

  elemental function cosh_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2,g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = cosh(A%dualz2)
    f_r%f3 = 3*g1*g2*cosh(g0) + g1**3*sinh(g0) + g3*sinh(g0)
  end function cosh_dz3

  elemental function tanh_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = tanh(A%dualz2)
    f_r%f3 = g3/cosh(g0)**2 - 2*g1**3 / cosh(g0)**4 - &
         6*g1*g2/cosh(g0)**2 * tanh(g0) + &
         4*g1**3/cosh(g0)**2 * tanh(g0)**2
  end function tanh_dz3

  elemental function asin_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = asin(A%dualz2)
    f_r%f3 = (3*g0**2*g1**3)/(1 - g0**2)**2.5 + g1**3/(1 - g0**2)**1.5+&
         (3*g0*g1*g2)/(1 - g0**2)**1.5 + g3/sqrt(1 - g0**2)
  end function asin_dz3

  elemental function acos_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = acos(A%dualz2)
    f_r%f3 = (-3*g0**2*g1**3)/(1 - g0**2)**2.5 -                  &
         g1**3/(1 - g0**2)**1.5 - (3*g0*g1*g2)/(1 - g0**2)**1.5 - &
         g3/sqrt(1 - g0**2)
  end function acos_dz3

  elemental function atan_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = atan(A%dualz2)
    f_r%f3 = (8*g0**2*g1**3)/(1 + g0**2)**3 - (2*g1**3)/(1 + g0**2)**2-&
         (6*g0*g1*g2)/(1 + g0**2)**2 + g3/(1 + g0**2)
  end function atan_dz3

  elemental function acosh_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3, g0m, g0p, sg0m, sg0p

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    g0m  = g0 - 1d0
    g0p  = g0 + 1d0
    sg0m = sqrt(g0m)
    sg0p = sqrt(g0p)

    f_r%dualz2 = acosh(A%dualz2)
    f_r%f3 =  (3*g1**3)/(4d0*sg0m*g0p**2.5) + g1**3/(2d0*g0m**1.5*     &
         g0p**1.5) + (3*g1**3)/(4d0*g0m**2.5*sg0p) - (3*g1*g2)/(2.*    &
         sg0m*g0p**1.5) - (3*g1*g2)/(2d0*g0m**1.5*sg0p) +              &
         g3/(sg0m*sg0p)
  end function acosh_dz3

  elemental function asinh_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = asinh(A%dualz2)
    f_r%f3 = (3*g0**2*g1**3)/(1 + g0**2)**2.5 - g1**3/(1 + g0**2)**1.5-&
         (3*g0*g1*g2)/(1 + g0**2)**1.5 + g3/sqrt(1 + g0**2)
  end function asinh_dz3

  elemental function atanh_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = atanh(A%dualz2)
    f_r%f3 = (8*g0**2*g1**3)/(1 - g0**2)**3 + (2*g1**3)/(1 - g0**2)**2+&
         (6*g0*g1*g2)/(1 - g0**2)**2 + g3/(1 - g0**2)
  end function atanh_dz3

  elemental function atan2_dz3(A,B) result(f_r)
    type(dualz3), intent(in) :: A, B
    type(dualz3) :: f_r
    complex(8) :: A0, A1, A2, A3, B0, B1, B2, B3, r2

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2
    A3 = A%f3

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2
    B3 = B%f3

    r2 = A0**2 + B0**2
    f_r%dualz2 = atan2(A%dualz2, B%dualz2)
    f_r%f3 = (8*(A1*B0 - A0*B1)*(A0*A1 + B0*B1)**2 + 2*(-(A1**3*B0) +  &
         A0*A1**2*B1 + A1*(-3*A0*A2*B0 + 2*A0**2*B2 - B0*(B1**2 +      &
         B0*B2)) + B1*(A0**2*A2 - 2*A2*B0**2 + A0*(B1**2 + 3*B0*B2)))* &
         r2 + (A3*B0 + A2*B1 - A1*B2 - A0*B3)*r2**2)/r2**3
  end function atan2_dz3

  elemental function sqrt_dz3(A) result(f_r)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_r
    complex(8) :: g0, g1, g2, g3

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2
    g3 = A%f3

    f_r%dualz2 = sqrt(A%dualz2)
    f_r%f3 = (3*g1**3)/(8d0*g0**2.5) - (3*g1*g2)/(4d0*g0**1.5) +       &
         g3/(2d0*sqrt(g0))
  end function sqrt_dz3

  !conjg
  !notice tat the conjugation operation is not, in general, differentia-
  !ble. The below definitions means (df)*, not d(f*), which is some
  !times useful
  elemental function conjg_dz3(g) result(f_r)
    type(dualz3), intent(in) :: g
    type(dualz3) :: f_r
    complex(8) :: g3  

    g3 = g%f3

    f_r%dualz2 = conjg(g%dualz2)
    f_r%f3 = conjg(g3)
  end function conjg_dz3

  !abs
  !The absolute value is not, in general, differentiable. In the below
  !function (df)* is used instead of d(f*). For instance d(A * A*)
  !is taking as dA * A* + A * (dA)*. See also the conjg_dz3 function
  elemental function abs_dz3(A) result(f_res)
    type(dualz3), intent(in) :: A
    type(dualz3) :: f_res

    f_res = sqrt(A*conjg(A))
  end function abs_dz3

  !|A> = |A0> + |A1> eps + |A2> eps2 + |A3>eps3; |B> = ...
  !As in the conjg and abs functions, here, the dual parts are (df)* not
  !d(f*)
  function dot_product_dz3(A,B) result(f_r)
    type(dualz3), intent(in), dimension(:) :: A, B
    type(dualz3) :: f_r
    complex(8), dimension(size(A)) :: A0, A1, A2, A3, B0, B1, B2, B3

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2
    A3 = A%f3

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2
    B3 = B%f3

    f_r%dualz2 = dot_product(A%dualz2,B%dualz2)
    f_r%f3 = dot_product(A0,B3) + 3*dot_product(A1,B2) + &
         3*dot_product(A2,B1) + dot_product(A3,B0)    
  end function dot_product_dz3

  !matmul
  function matmul22_dz3(A,B) result(f_r)
    type(dualz3), intent(in), dimension(:,:) :: A, B
    type(dualz3), dimension(size(A,1),size(B,2)) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A0, A1, A2, A3
    complex(8), dimension(size(B,1),size(B,2)) :: B0, B1, B2, B3

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2
    A3 = A%f3

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2
    B3 = B%f3

    f_r%dualz2 = matmul(A%dualz2,B%dualz2)
    f_r%f3 = matmul(A0,B3) + 3*matmul(A1,B2) + &
         3*matmul(A2,B1) + matmul(A3,B0) 
  end function matmul22_dz3

  function matmul22_dz3C(A,B) result(fr)
    type(dualz3), intent(in), dimension(:,:) :: A
    complex(8), intent(in), dimension(:,:) :: B
    type(dualz3), dimension(size(A,1),size(B,2)) :: fr
    type(dualz3), dimension(size(B,1),size(B,2)) :: Bdz3

    Bdz3 = B

    fr = matmul22_dz3(A,Bdz3)
  end function matmul22_dz3C

  function matmul22_dz3R(A,B) result(fr)
    type(dualz3), intent(in), dimension(:,:) :: A
    real(8), intent(in), dimension(:,:) :: B
    type(dualz3), dimension(size(A,1),size(B,2)) :: fr
    type(dualz3), dimension(size(B,1),size(B,2)) :: Bdz3

    Bdz3 = B

    fr = matmul22_dz3(A,Bdz3)
  end function matmul22_dz3R

  function matmul22_Cdz3(A,B) result(fr)
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz3), intent(in), dimension(:,:) :: B
    type(dualz3), dimension(size(A,1),size(B,2)) :: fr
    type(dualz3), dimension(size(A,1),size(A,2)) :: Adz3

    Adz3 = A
    fr = matmul22_dz3(Adz3,B)
  end function matmul22_Cdz3

  function matmul22_Rdz3(A,B) result(fr)
    real(8), intent(in), dimension(:,:) :: A
    type(dualz3), intent(in), dimension(:,:) :: B
    type(dualz3), dimension(size(A,1),size(B,2)) :: fr
    type(dualz3), dimension(size(A,1),size(A,2)) :: Adz3

    Adz3 = A
    fr = matmul22_dz3(Adz3,B)
  end function matmul22_Rdz3  

  !dual_rank2*dual_rank1
  function matmul_dr21_dz3(A,x) result(f_r)
    type(dualz3), intent(in), dimension(:,:) :: A
    type(dualz3), intent(in), dimension(size(A,2)) :: x
    type(dualz3), dimension(size(A,1)) :: f_r
    integer :: m, n

    m = size(A,1)
    n = size(A,2)
    f_r = reshape(matmul(A,reshape(x,[n,1])),[m])
  end function matmul_dr21_dz3

  function matmul_dr21_dz3C(A,x) result(fr)
    type(dualz3), intent(in), dimension(:,:) :: A
    complex(8), intent(in), dimension(size(A,2)) :: x
    type(dualz3), dimension(size(A,1)) :: fr
    type(dualz3), dimension(size(A,2)) :: xd

    xd = x
    fr = matmul_dr21_dz3(A,xd)
  end function matmul_dr21_dz3C

  function matmul_dr21_dz3R(A,x) result(fr)
    type(dualz3), intent(in), dimension(:,:) :: A
    real(8), intent(in), dimension(size(A,2)) :: x
    type(dualz3), dimension(size(A,1)) :: fr
    type(dualz3), dimension(size(A,2)) :: xd

    xd = x
    fr = matmul_dr21_dz3(A,xd)
  end function matmul_dr21_dz3R

  function matmul_dr21_Cdz3(A,x) result(fr)
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz3), intent(in), dimension(size(A,2)) :: x
    type(dualz3), dimension(size(A,1)) :: fr
    type(dualz3), dimension(size(A,1),size(A,2)) :: Adz3

    Adz3 = A
    fr = matmul_dr21_dz3(Adz3,x)
  end function matmul_dr21_Cdz3

  function matmul_dr21_Rdz3(A,x) result(fr)
    real(8), intent(in), dimension(:,:) :: A
    type(dualz3), intent(in), dimension(size(A,2)) :: x
    type(dualz3), dimension(size(A,1)) :: fr
    type(dualz3), dimension(size(A,1),size(A,2)) :: Adz3

    Adz3 = A
    fr = matmul_dr21_dz3(Adz3,x)
  end function matmul_dr21_Rdz3  

  !in order to avoid ambiguity we must change x to xx
  !"A generic function must be able to distinguish its arguments by 
  !type AND by name"
  function matmul_dr12_dz3(xx,A) result(f_r)
    type(dualz3), intent(in), dimension(:,:) :: A
    type(dualz3), intent(in), dimension(:) :: xx
    type(dualz3), dimension(size(A,2)) :: f_r
    integer :: m, n

    m = size(A,1)
    n = size(A,2)

    !if size(xx) /= m you may want to use some warning, the below code
    !relies on the reshape fortran function which uses padding
!!$if(m  /= size(xx)) print*,'warning on matmul_dr12_dz3'
    f_r = reshape(matmul(reshape(xx,[1,m]),A),[n])
  end function matmul_dr12_dz3

  function matmul_dr12_dz3C(xx,A) result(fr)
    type(dualz3), intent(in), dimension(:) :: xx
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz3), dimension(size(A,2)) :: fr
    type(dualz3), dimension(size(A,1),size(A,2)) :: Adz3

    Adz3 = A
    fr =  matmul_dr12_dz3(xx,Adz3)
  end function matmul_dr12_dz3C

  function matmul_dr12_dz3R(xx,A) result(fr)
    type(dualz3), intent(in), dimension(:) :: xx
    real(8), intent(in), dimension(:,:) :: A
    type(dualz3), dimension(size(A,2)) :: fr
    type(dualz3), dimension(size(A,1),size(A,2)) :: Adz3

    Adz3 = A
    fr =  matmul_dr12_dz3(xx,Adz3)
  end function matmul_dr12_dz3R

  function matmul_dr12_Cdz3(xx,A) result(fr)
    complex(8), intent(in), dimension(:) :: xx
    type(dualz3), intent(in), dimension(:,:) :: A
    type(dualz3), dimension(size(A,2)) :: fr
    type(dualz3), dimension(size(xx)) :: xxdz3

    xxdz3 = xx
    fr =  matmul_dr12_dz3(xxdz3,A)
  end function matmul_dr12_Cdz3

  function matmul_dr12_Rdz3(xx,A) result(fr)
    real(8), intent(in), dimension(:) :: xx
    type(dualz3), intent(in), dimension(:,:) :: A
    type(dualz3), dimension(size(A,2)) :: fr
    type(dualz3), dimension(size(xx)) :: xxdz3

    xxdz3 = xx
    fr =  matmul_dr12_dz3(xxdz3,A)
  end function matmul_dr12_Rdz3    

  !sum for rank2 arrays
  !the result is given as array of rank 1
  function sum2_dz3(A,k) result(f_r)
    type(dualz3), intent(in), dimension(:,:) :: A
    integer, intent(in) :: k
    type(dualz3), dimension(size(A,2/k)) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A3

    A3 = A%f3

    f_r%dualz2 = sum(A%dualz2,k)
    f_r%f3 = sum(A3,k)
  end function sum2_dz3

  !sum(A) rank2
  function sum0_dz3(A) result(f_r)
    type(dualz3), intent(in), dimension(:,:) :: A
    type(dualz3) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A3

    A3 = A%f3

    f_r%dualz2 = sum(A%dualz2)
    f_r%f3 = sum(A3)
  end function  sum0_dz3

  !sum(A) rank1
  function sum10_dz3(A) result(f_r)
    type(dualz3), intent(in), dimension(:) :: A
    type(dualz3) :: f_r    
    complex(8), dimension(size(A)) :: A3

    A3 = A%f3

    f_r%dualz2 = sum(A%dualz2)
    f_r%f3 = sum(A3)
  end function  sum10_dz3

  !product function for dualz1 vectors (arrays of rank 1)
  function prod1_dz3(x) result(f_r)
    type(dualz3), intent(in), dimension(:) :: x
    type(dualz3) :: f_r
    integer :: k

    f_r = 1
    do k = 1, size(x)
       f_r = f_r * x(k)
    end do
  end function prod1_dz3

  !product(x) for arrays of rank 2
  function prod0_dz3(x) result(f_r)
    type(dualz3), intent(in), dimension(:,:) :: x
    type(dualz3) :: f_r

    f_r = prod1_dz3(reshape(x,[size(x)]))
  end function prod0_dz3

  !product(A,k) for rank2 arrays
  !the result is given as array of rank 1
  function prod2_dz3(x,c) result(f_r)
    type(dualz3), intent(in), dimension(:,:) :: x
    integer, intent(in) :: c
    type(dualz3), dimension(size(x,2/c)) :: f_r
    type(dualz3), dimension(size(x,c)) :: xc
    integer :: k

    if(c==1) then
       do k = 1, size(x,2)
          xc = x(:,k)
          f_r(k) = prod1_dz3(xc)
       end do
    else if(c==2) then
       do k = 1, size(x,1)
          xc = x(k,:)
          f_r(k) = prod1_dz3(xc)
       end do
    else 
       stop 'c must be equal to 1 (2) to collapse rows (columns)'
    end if
  end function prod2_dz3

  !=====================================================================
  !auxiliar functions
  elemental function f0part_dz3(A) result(fr)
    type(dualz3), intent(in) :: A
    complex(8) :: fr

    fr = A%f0
  end function f0part_dz3

  elemental function f1part_dz3(A) result(fr)
    type(dualz3), intent(in) :: A
    complex(8) :: fr

    fr = A%f1
  end function f1part_dz3

  elemental function f2part_dz3(A) result(fr)
    type(dualz3), intent(in) :: A
    complex(8) :: fr

    fr = A%f2
  end function f2part_dz3

  elemental function f3part_dz3(A) result(fr)
    type(dualz3), intent(in) :: A
    complex(8) :: fr
    
    fr = A%f3
  end function f3part_dz3  

  !auxiliar function to convert to dualz3
  !in principle this function is not necessary (the assignment operator
  !does the "same" job) but is provided for convenience
  !numeric (including dualz1 and dualz2) to dualz3
  elemental function dz3todualz3(z) result(fr)
    type(dualz3), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function dz3todualz3

  elemental function dz2todualz3(z) result(fr)
    type(dualz2), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function dz2todualz3

  elemental function dz1todualz3(z) result(fr)
    type(dualz1), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function dz1todualz3

  elemental function z8todualz3(z) result(fr)
    complex(8), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function z8todualz3

  elemental function z4todualz3(z) result(fr)
    complex(4), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function z4todualz3

  elemental function r8todualz3(z) result(fr)
    real(8), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function r8todualz3

  elemental function r4todualz3(z) result(fr)
    real(4), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function r4todualz3

  elemental function i8todualz3(z) result(fr)
    integer(8), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function i8todualz3

  elemental function i4todualz3(z) result(fr)
    integer(4), intent(in) :: z
    type(dualz3) :: fr

    fr = z
  end function i4todualz3
  !=====================================================================
end module dualz3_mod

