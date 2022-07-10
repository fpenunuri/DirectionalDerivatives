!Dual numbers of complex components
!This can be used to compute zeroth, first, and second order derivatives
!F. Peñuñuri
!UADY, Merida Yucatan Mexico
!2021

module dualz2_mod
  use dualz1_mod
  implicit none
  private

  !dualz1 ---> [f0,f1]
  !dualz2 ---> [f0,f1,f2]
  type, public, extends (dualz1) :: dualz2
     complex(8) :: f2
  end type dualz2

  public :: assignment (=)

  public :: operator(==), operator(/=), operator(+), operator(-),    &
       operator(*), operator (/), operator(**)

  public :: sin, cos, tan, sinh, cosh, tanh, exp, log, asin, acos,   &
       atan, acosh, asinh, atanh, atan2, sqrt, conjg, abs, matmul,   &
       dot_product, sum, product

  public :: f0part_dz2, f1part_dz2, f2part_dz2
  public :: todualz2
  
  !INTERFACES
  !equal assignment
  interface assignment (=)
     module procedure igualz1_dz2  !dualz2 <--- dualz1
     module procedure igualc8_dz2  !dualz2 <--- complex8
     module procedure igualc4_dz2  !dualz2 <--- complex4
     module procedure igualr8_dz2  !dualz2 <--- real8
     module procedure igualr4_dz2  !dualz2 <--- real4
     module procedure iguali8_dz2  !dualz2 <--- integer8
     module procedure iguali4_dz2  !dualz2 <--- integer4
  end interface assignment (=)

  !Logical equal operator 
  interface operator (==) 
     module procedure eq_dz2     !dualz2 == dualz2
  end interface operator (==)

  !Logical not equal operator 
  interface operator (/=) 
     module procedure noteq_dz2  !dualz2 /= dualz2
  end interface operator (/=)

  !Plus operator
  interface operator (+)
     module procedure mas_dz2      !unary 
     module procedure suma_dz2     !dualz2 + dualz2
     module procedure sumaz2z1_dz2 !dualz2 + dualz1
     module procedure sumaz2c8_dz2 !dualz2 + complex8
     module procedure sumaz2c4_dz2 !dualz2 + complex4
     module procedure sumaz2r8_dz2 !dualz2 + real8
     module procedure sumaz2r4_dz2 !dualz2 + real4
     module procedure sumaz2i8_dz2 !dualz2 + integer8
     module procedure sumaz2i4_dz2 !dualz2 + integer4
     module procedure sumaz1z2_dz2 !dualz1 + dualz2
     module procedure sumac8z2_dz2 !complex8 + dualz2
     module procedure sumac4z2_dz2 !complex4 + dualz2
     module procedure sumar8z2_dz2 !real8 + dualz2
     module procedure sumar4z2_dz2 !real4 + dualz2
     module procedure sumai8z2_dz2 !integer8 + dualz2
     module procedure sumai4z2_dz2 !integer4 + dualz2
  end interface operator (+)

  !minus operator
  interface operator (-)
     module procedure menos_dz2     !unary 
     module procedure resta_dz2     !dualz2 - dualz2
     module procedure restaz2z1_dz2 !dualz2 - dualz1
     module procedure restaz2c8_dz2 !dualz2 - complex8
     module procedure restaz2c4_dz2 !dualz2 - complex4
     module procedure restaz2r8_dz2 !dualz2 - real8
     module procedure restaz2r4_dz2 !dualz2 - real4
     module procedure restaz2i8_dz2 !dualz2 - integer8
     module procedure restaz2i4_dz2 !dualz2 - integer4
     module procedure restaz1z2_dz2 !dualz1 - dualz2
     module procedure restac8z2_dz2 !complex8 - dualz2
     module procedure restac4z2_dz2 !complex4 - dualz2
     module procedure restar8z2_dz2 !real8 - dualz2
     module procedure restar4z2_dz2 !real4 - dualz2
     module procedure restai8z2_dz2 !integer8 - dualz2
     module procedure restai4z2_dz2 !integer4 - dualz2
  end interface operator (-)

  !Times operator
  interface operator (*)
     module procedure times_dz2     !dualz2 * dualz2
     module procedure timesz2z1_dz2 !dualz2 * dualz1
     module procedure timesz2c8_dz2 !dualz2 * complex8
     module procedure timesz2c4_dz2 !dualz2 * complex4
     module procedure timesz2r8_dz2 !dualz2 * real8
     module procedure timesz2r4_dz2 !dualz2 * real4
     module procedure timesz2i8_dz2 !dualz2 * integer8
     module procedure timesz2i4_dz2 !dualz2 * integer4
     module procedure timesz1z2_dz2 !dualz1 * dualz2
     module procedure timesc8z2_dz2 !complex8 * dualz2
     module procedure timesc4z2_dz2 !complex4 * dualz2
     module procedure timesr8z2_dz2 !real8  dualz2
     module procedure timesr4z2_dz2 !real4 * dualz2
     module procedure timesi8z2_dz2 !integer8 * dualz2
     module procedure timesi4z2_dz2 !integer4 * dualz2
  end interface operator (*)

  !Division operator
  interface operator (/)
     module procedure div_dz2     !dualz2 / dualz2
     module procedure divz2z1_dz2 !dualz2 / dualz1
     module procedure divz2c8_dz2 !dualz2 / complex8
     module procedure divz2c4_dz2 !dualz2 / complex4
     module procedure divz2r8_dz2 !dualz2 / real8
     module procedure divz2r4_dz2 !dualz2 / real4
     module procedure divz2i8_dz2 !dualz2 / integer8
     module procedure divz2i4_dz2 !dualz2 / integer4
     module procedure divz1z2_dz2 !dualz1 / dualz2
     module procedure divc8z2_dz2 !complex8 / dualz2
     module procedure divc4z2_dz2 !complex4 / dualz2
     module procedure divr8z2_dz2 !real8 / dualz2
     module procedure divr4z2_dz2 !real4 / dualz2
     module procedure divi8z2_dz2 !integer8 / dualz2
     module procedure divi4z2_dz2 !integer4 / dualz2
  end interface operator (/)

  !Power operator
  interface operator (**)
     module procedure power_dz2     !dualz2 ** dualz2
     module procedure powerz2z1_dz2 !dualz2 ** dualz1
     module procedure powerz2c8_dz2 !dualz2 ** complex8
     module procedure powerz2c4_dz2 !dualz2 ** complex4
     module procedure powerz2r8_dz2 !dualz2 ** real8
     module procedure powerz2r4_dz2 !dualz2 ** real4
     module procedure powerz2i8_dz2 !dualz2 ** integer8
     module procedure powerz2i4_dz2 !dualz2 ** integer4
     module procedure powerz1z2_dz2 !dualz1 ** dualz2
     module procedure powerc8z2_dz2 !complex8 ** dualz2
     module procedure powerc4z2_dz2 !complex4 ** dualz2
     module procedure powerr8z2_dz2 !real8  dualz2
     module procedure powerr4z2_dz2 !real4 ** dualz2
     module procedure poweri8z2_dz2 !integer8 ** dualz2
     module procedure poweri4z2_dz2 !integer4 ** dualz2
  end interface operator (**)

  !Mathematical Functions
  interface abs
     module procedure abs_dz2
  end interface abs

  interface conjg
     module procedure conjg_dz2
  end interface conjg

  interface sin
     module procedure sin_dz2
  end interface sin

  interface cos
     module procedure cos_dz2
  end interface cos

  interface tan
     module procedure tan_dz2
  end interface tan

  interface sinh
     module procedure sinh_dz2
  end interface sinh

  interface cosh
     module procedure cosh_dz2
  end interface cosh

  interface tanh
     module procedure tanh_dz2
  end interface tanh

  interface exp
     module procedure exp_dz2
  end interface exp

  interface log
     module procedure log_dz2
  end interface log

  interface asin
     module procedure asin_dz2
  end interface asin

  interface acos
     module procedure acos_dz2
  end interface acos

  interface atan
     module procedure atan_dz2
  end interface atan

  interface asinh
     module procedure asinh_dz2
  end interface asinh

  interface acosh
     module procedure acosh_dz2
  end interface acosh

  interface atanh
     module procedure atanh_dz2
  end interface atanh

  interface atan2
     module procedure atan2_dz2
  end interface atan2

  interface sqrt
     module procedure sqrt_dz2
  end interface sqrt

  !some vector and matrix functions
  !<A|B>
  interface dot_product
     module procedure dot_product_dz2
  end interface dot_product

  !Matrix multiplication
  interface matmul
     module procedure matmul22_dz2    !dual_rank2*dual_rank2
     module procedure matmul_dr21_dz2 !dual_rank2*dual_rank1
     module procedure matmul_dr12_dz2 !dual_rank1*dual_rank2
     !
     module procedure matmul22_dz2C
     module procedure matmul22_dz2R
     module procedure matmul22_Cdz2
     module procedure matmul22_Rdz2
     !
     module procedure matmul_dr21_dz2C
     module procedure matmul_dr21_dz2R
     module procedure matmul_dr21_Cdz2
     module procedure matmul_dr21_Rdz2
     !
     module procedure matmul_dr12_dz2C
     module procedure matmul_dr12_dz2R
     module procedure matmul_dr12_Cdz2
     module procedure matmul_dr12_Rdz2
  end interface matmul

  interface sum
     module procedure sum2_dz2  !sum(dual_rank2,k)
     module procedure sum0_dz2  !sum(dual_rank2)
     module procedure sum10_dz2 !sum(dual_rank1)
  end interface sum

  !product for dual vectors x(k)
  interface product
     module procedure prod1_dz2
     module procedure prod0_dz2
     module procedure prod2_dz2 !product(dual_rank2,k)
  end interface product

  !auxiliar function  
  interface todualz2
     module procedure dz2todualz2
     module procedure dz1todualz2
     module procedure z8todualz2
     module procedure z4todualz2
     module procedure r8todualz2
     module procedure r4todualz2
     module procedure i8todualz2
     module procedure i4todualz2
  end interface todualz2
  !terminan interfaces

  !Functions and subroutines associated to the above interfaces
contains
  !Assignment, equal operator
  !dualz2 <--- dualz1
  elemental subroutine igualz1_dz2(A, z)
    type(dualz2), intent(out) :: A
    type(dualz1), intent(in) :: z

    A%dualz1 = z
    A%f2 = 0
  end subroutine igualz1_dz2

  !dualz2 <--- complex8
  elemental subroutine igualc8_dz2(A, z)
    type(dualz2), intent(out) :: A
    complex(8), intent(in) :: z

    A%dualz1 = z
    A%f2 = 0
  end subroutine igualc8_dz2

  !dualz2 <--- complex4
  elemental subroutine igualc4_dz2(A, z)
    type(dualz2), intent(out) :: A
    complex(4), intent(in) :: z

    A%dualz1 = z
    A%f2 = 0
  end subroutine igualc4_dz2

  !dualz2 <--- real8
  elemental subroutine igualr8_dz2(A, z)
    type(dualz2), intent(out) :: A
    real(8), intent(in) :: z

    A%dualz1 = z
    A%f2 = 0
  end subroutine igualr8_dz2

  !dualz2 <--- real4
  elemental subroutine igualr4_dz2(A, z)
    type(dualz2), intent(out) :: A
    real(4), intent(in) :: z

    A%dualz1 = z
    A%f2 = 0
  end subroutine igualr4_dz2

  !dualz2 <--- integer8
  elemental subroutine iguali8_dz2(A, z)
    type(dualz2), intent(out) :: A
    integer(8), intent(in) :: z

    A%dualz1 = z
    A%f2 = 0
  end subroutine iguali8_dz2

  !dualz2 <--- integer4
  elemental subroutine iguali4_dz2(A, z)
    type(dualz2), intent(out) :: A
    integer(4), intent(in) :: z

    A%dualz1 = z
    A%f2 = 0
  end subroutine iguali4_dz2

  !Logical not equal operator
  elemental function noteq_dz2(lhs, rhs) result(f_res)
    type (dualz2), intent(in) :: lhs, rhs
    logical :: f_res

    f_res = .not.(lhs == rhs)
  end function noteq_dz2

  !Logical equal operator
  elemental function eq_dz2(lhs, rhs) result(f_res)
    type (dualz2), intent(in) :: lhs, rhs
    logical :: f_res
    logical :: eqf0, eqf1, eqf2

    eqf0 = lhs%f0 == rhs%f0
    eqf1 = lhs%f1 == rhs%f1
    eqf2 = lhs%f2 == rhs%f2

    f_res = all([eqf0,eqf1,eqf2])
  end function eq_dz2

  !+dualz2 (unary)
  elemental function mas_dz2(A) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_res

    f_res = A
  end function mas_dz2

  !Power
  !dualz2**dualz2
  !el caso A0 = 0 es problematico, solo se atendio este caso para cuando
  !1.- [A1,A2] = 0 y la parte real de B0 > 0 
  !2.- la parte real de B0 > 2
  !3.- B1 = 0
  !en todos estos casos no hay restricciones para las demas componentes
  !es decir, el caso [B1,B2] = 0 no se requiere
  elemental function power_dz2(A,B) result(fr)
    type(dualz2), intent(in) :: A, B
    type(dualz2) :: fr
    complex(8) :: A0, A1, A2, B0, B1, B2

    A0 = A%f0; A1 = A%f1; A2 = A%f2; B0 = B%f0; B1 = B%f1; B2 = B%f2

    fr%dualz1 = A%dualz1**B%dualz1

    if(all([B0,B1,B2]==0)) then
       fr%f2 = 0
    elseif(all([A0,A1,A2]==0) .and. real(B0) > 0) then
       fr%f2 = 0 
    elseif(all([A1,A2,B1,B2]==0)) then
       fr%f2 = 0 
    elseif(A0 == 0 .and. real(B0)>2) then
       fr%f2 = 0 
    elseif(all([A0,B1]==0)) then 
       if(B0 == 1) then
          fr%f2 = A2
       elseif(B0 == 2) then
          fr%f2 = 2*A1**2
       else
          fr%f2 =  A0**B0*((A1*B0)/A0 + B1*log(A0))**2 + &
               A0**B0*(-((A1**2*B0)/A0**2) + (A2*B0)/A0 +  &
               (2*A1*B1)/A0 + B2*log(A0))
       end if
    else
       fr%f2 = A0**B0*((A1*B0)/A0 + B1*log(A0))**2 + &
            A0**B0*(-((A1**2*B0)/A0**2) + (A2*B0)/A0 +  &
            (2*A1*B1)/A0 + B2*log(A0))
    end if
  end function power_dz2

  !dualz2**dualz1
  elemental function powerz2z1_dz2(A,Bz1) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz1),intent(in) :: Bz1
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = Bz1
    f_res = A**Bz2
  end function powerz2z1_dz2

  !dualz2**complex8
  elemental function powerz2c8_dz2(A,c8) result(f_res)
    type(dualz2), intent(in) :: A
    complex(8),intent(in) :: c8
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = c8
    f_res = A**Bz2
  end function powerz2c8_dz2

  !dualz2**complex4
  elemental function powerz2c4_dz2(A,c4) result(f_res)
    type(dualz2), intent(in) :: A
    complex(4),intent(in) :: c4
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = c4
    f_res = A**Bz2
  end function powerz2c4_dz2

  !dualz2**real8
  elemental function powerz2r8_dz2(A,r8) result(f_res)
    type(dualz2), intent(in) :: A
    real(8),intent(in) :: r8
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = r8
    f_res = A**Bz2
  end function powerz2r8_dz2

  !dualz2**real4
  elemental function powerz2r4_dz2(A,r4) result(f_res)
    type(dualz2), intent(in) :: A
    real(4),intent(in) :: r4
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = r4
    f_res = A**Bz2
  end function powerz2r4_dz2

  !dualz2**integer8
  elemental function powerz2i8_dz2(A,i8) result(f_res)
    type(dualz2), intent(in) :: A
    integer(8),intent(in) :: i8
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = i8
    f_res = A**Bz2
  end function powerz2i8_dz2

  !dualz2**integer4
  elemental function powerz2i4_dz2(A,i4) result(f_res)
    type(dualz2), intent(in) :: A
    integer(4),intent(in) :: i4
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = i4
    f_res = A**Bz2
  end function powerz2i4_dz2

  !dualz1**dualz2
  elemental function powerz1z2_dz2(Az1,Bz2) result(f_res)
    type(dualz1), intent(in) :: Az1
    type(dualz2),intent(in) :: Bz2
    type(dualz2) :: f_res
    type(dualz2) :: Az2

    Az2 = Az1
    f_res = Az2**Bz2
  end function powerz1z2_dz2

  !complex8 ** dualz2
  elemental function powerc8z2_dz2(cc8,Bz2) result(f_res)
    complex(8), intent(in) :: cc8
    type(dualz2),intent(in) :: Bz2
    type(dualz2) :: f_res
    type(dualz2) :: Ac8

    Ac8 = cc8
    f_res = Ac8**Bz2
  end function powerc8z2_dz2

  !complex4 ** dualz2
  elemental function powerc4z2_dz2(cc4,Bz2) result(f_res)
    complex(4), intent(in) :: cc4
    type(dualz2),intent(in) :: Bz2
    type(dualz2) :: f_res
    type(dualz2) :: Ac4

    Ac4 = cc4
    f_res = Ac4**Bz2
  end function powerc4z2_dz2

  !real8**dualz2
  elemental function powerr8z2_dz2(rr8,Bz2) result(f_res)
    real(8), intent(in) :: rr8
    type(dualz2),intent(in) :: Bz2
    type(dualz2) :: f_res
    type(dualz2) :: Ar8

    Ar8 = rr8
    f_res = Ar8**Bz2
  end function powerr8z2_dz2

  !real4**dualz2
  elemental function powerr4z2_dz2(rr4,Bz2) result(f_res)
    real(4), intent(in) :: rr4
    type(dualz2),intent(in) :: Bz2
    type(dualz2) :: f_res
    type(dualz2) :: Ar4

    Ar4 = rr4
    f_res = Ar4**Bz2
  end function powerr4z2_dz2

  !integer8**dualz2
  elemental function poweri8z2_dz2(ii8,Bz2) result(f_res)
    integer(8), intent(in) :: ii8
    type(dualz2),intent(in) :: Bz2
    type(dualz2) :: f_res
    type(dualz2) :: Ai8

    Ai8 = ii8
    f_res = Ai8**Bz2
  end function poweri8z2_dz2

  !integer4**dualz2
  elemental function poweri4z2_dz2(ii4,Bz2) result(f_res)
    integer(4), intent(in) :: ii4
    type(dualz2),intent(in) :: Bz2
    type(dualz2) :: f_res
    type(dualz2) :: Ai4

    Ai4 = ii4
    f_res = Ai4**Bz2
  end function poweri4z2_dz2

  !Division
  !dualz2/dualz2
  elemental function div_dz2(A,B) result(f_res)
    type(dualz2), intent(in) :: A, B
    type(dualz2) :: f_res
    type(dualz2) :: BI
    complex(8) :: B0, B1, B2

    B0 = B%f0; B1 = B%f1; B2 = B%f2

    BI%f0 = 1/B0
    BI%f1 = -B1/B0**2
    BI%f2 = (2*B1**2)/B0**3 - B2/B0**2

    f_res = A*BI    
  end function div_dz2

  !dualz2/dualz1
  elemental function divz2z1_dz2(A,Bz1) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz1), intent(in) :: Bz1
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = Bz1
    f_res = A/Bz2      
  end function divz2z1_dz2

  !dualz2/complex8
  elemental function divz2c8_dz2(A,c8) result(f_res)
    type(dualz2), intent(in) :: A
    complex(8), intent(in) :: c8
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = c8
    f_res = A/Bz2   
  end function divz2c8_dz2

  !dualz2/complex4
  elemental function divz2c4_dz2(A,c4) result(f_res)
    type(dualz2), intent(in) :: A
    complex(4), intent(in) :: c4
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = c4
    f_res = A/Bz2   
  end function divz2c4_dz2

  !dualz2/real8
  elemental function divz2r8_dz2(A,r8) result(f_res)
    type(dualz2), intent(in) :: A
    real(8), intent(in) :: r8
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = r8
    f_res = A/Bz2   
  end function divz2r8_dz2

  !dualz2/real4
  elemental function divz2r4_dz2(A,r4) result(f_res)
    type(dualz2), intent(in) :: A
    real(4), intent(in) :: r4
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = r4
    f_res = A/Bz2   
  end function divz2r4_dz2

  !dualz2/integer8
  elemental function divz2i8_dz2(A,i8) result(f_res)
    type(dualz2), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = i8
    f_res = A/Bz2   
  end function divz2i8_dz2

  !dualz2/integer4
  elemental function divz2i4_dz2(A,i4) result(f_res)
    type(dualz2), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = i4
    f_res = A/Bz2   
  end function divz2i4_dz2

  !dualz1 / dualz2
  elemental function divz1z2_dz2(Az1,B) result(f_res)
    type(dualz2), intent(in) :: B
    type(dualz1), intent(in) :: Az1
    type(dualz2) :: f_res
    type(dualz2) :: Az2

    Az2 = Az1     
    f_res = Az2/B      
  end function divz1z2_dz2

  !complex8 / dualz2
  elemental function divc8z2_dz2(Ac8,B) result(f_res)
    type(dualz2), intent(in) :: B
    complex(8), intent(in) :: Ac8
    type(dualz2) :: f_res
    type(dualz2) :: Ac8z

    Ac8z = Ac8
    f_res = Ac8z/B     
  end function divc8z2_dz2

  !complex4 / dualz2
  elemental function divc4z2_dz2(Ac4,B) result(f_res)
    type(dualz2), intent(in) :: B
    complex(4), intent(in) :: Ac4
    type(dualz2) :: f_res
    type(dualz2) :: Ac4z

    Ac4z = Ac4
    f_res = Ac4z/B           
  end function divc4z2_dz2

  !real8 / dualz2
  elemental function divr8z2_dz2(Ar8,B) result(f_res)
    type(dualz2), intent(in) :: B
    real(8), intent(in) :: Ar8
    type(dualz2) :: f_res
    type(dualz2) :: Ar8z

    Ar8z = Ar8
    f_res = Ar8z/B
  end function divr8z2_dz2

  !real4/dualz2
  elemental function divr4z2_dz2(Ar4,B) result(f_res)
    type(dualz2), intent(in) :: B
    real(4), intent(in) :: Ar4
    type(dualz2) :: f_res
    type(dualz2) :: Ar4z

    Ar4z = Ar4
    f_res = Ar4z/B
  end function divr4z2_dz2

  !integer8/dualz2
  elemental function divi8z2_dz2(Ai8,B) result(f_res)
    type(dualz2), intent(in) :: B
    integer(8), intent(in) :: Ai8
    type(dualz2) :: f_res
    type(dualz2) :: Ai8z

    Ai8z = Ai8
    f_res = Ai8z/B
  end function divi8z2_dz2

  !integer4/dualz2
  elemental function divi4z2_dz2(Ai4,B) result(f_res)
    type(dualz2), intent(in) :: B
    integer(4), intent(in) :: Ai4
    type(dualz2) :: f_res
    type(dualz2) :: Ai4z

    Ai4z = Ai4
    f_res = Ai4z/B
  end function divi4z2_dz2

  !Times
  !dualz2*dualz2
  elemental function times_dz2(A,B) result(f_res)
    type(dualz2), intent(in) :: A, B
    type(dualz2) :: f_res      

    f_res%dualz1 = A%dualz1 * B%dualz1
    f_res%f2 =  A%f0 * B%f2 + A%f2 * B%f0 + 2*A%f1*B%f1
  end function times_dz2

  !dualz2*dualz1
  elemental function timesz2z1_dz2(A,B) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz1), intent(in) :: B
    type(dualz2) :: f_res
    type(dualz2) :: Bz2

    Bz2 = B
    f_res = A*Bz2
  end function  timesz2z1_dz2

  !dualz2*complex8
  elemental function timesz2c8_dz2(A,c8) result(f_res)
    type(dualz2), intent(in) :: A
    complex(8), intent(in) :: c8
    type(dualz2) :: f_res      
    type(dualz2) :: c8z2

    c8z2 = c8
    f_res = A*c8z2
  end function  timesz2c8_dz2

  !dualz2*complex4
  elemental function timesz2c4_dz2(A,c4) result(f_res)
    type(dualz2), intent(in) :: A
    complex(4), intent(in) :: c4
    type(dualz2) :: f_res      
    type(dualz2) :: c4z2

    c4z2 = c4
    f_res = A*c4z2
  end function  timesz2c4_dz2

  !dualz2*real8
  elemental function timesz2r8_dz2(A,r8) result(f_res)
    type(dualz2), intent(in) :: A
    real(8), intent(in) :: r8
    type(dualz2) :: f_res      
    type(dualz2) :: r8z2

    r8z2 = r8
    f_res = A*r8z2
  end function  timesz2r8_dz2

  !dualz2*real4
  elemental function timesz2r4_dz2(A,r4) result(f_res)
    type(dualz2), intent(in) :: A
    real(4), intent(in) :: r4
    type(dualz2) :: f_res      
    type(dualz2) :: r4z2

    r4z2 = r4
    f_res = A*r4z2
  end function  timesz2r4_dz2

  !dualz2*integer8
  elemental function timesz2i8_dz2(A,i8) result(f_res)
    type(dualz2), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz2) :: f_res      
    type(dualz2) :: i8z2

    i8z2 = i8
    f_res = A*i8z2
  end function  timesz2i8_dz2

  !dualz2*integer4
  elemental function timesz2i4_dz2(A,i4) result(f_res)
    type(dualz2), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz2) :: f_res      
    type(dualz2) :: i4z2

    i4z2 = i4
    f_res = A*i4z2
  end function  timesz2i4_dz2

  !-----
  !dualz1 * dualz2
  elemental function timesz1z2_dz2(Bzz1,A) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz1), intent(in) :: Bzz1
    type(dualz2) :: f_res      

    f_res = A*Bzz1 
  end function  timesz1z2_dz2

  !complex8 * dualz2
  elemental function timesc8z2_dz2(cc8,A) result(f_res)
    type(dualz2), intent(in) :: A
    complex(8), intent(in) :: cc8
    type(dualz2) :: f_res      

    f_res = A*cc8
  end function  timesc8z2_dz2

  !complex4 * dualz2
  elemental function timesc4z2_dz2(cc4,A) result(f_res)
    type(dualz2), intent(in) :: A
    complex(4), intent(in) :: cc4
    type(dualz2) :: f_res      

    f_res = A*cc4
  end function  timesc4z2_dz2

  !real8 * dualz2
  elemental function timesr8z2_dz2(rr8,A) result(f_res)
    type(dualz2), intent(in) :: A
    real(8), intent(in) :: rr8
    type(dualz2) :: f_res      

    f_res = A*rr8
  end function  timesr8z2_dz2

  !real4 * dualz2
  elemental function timesr4z2_dz2(rr4,A) result(f_res)
    type(dualz2), intent(in) :: A
    real(4), intent(in) :: rr4
    type(dualz2) :: f_res      

    f_res = A*rr4
  end function  timesr4z2_dz2

  !integer8 * dualz2
  elemental function timesi8z2_dz2(ii8,A) result(f_res)
    type(dualz2), intent(in) :: A
    integer(8), intent(in) :: ii8
    type(dualz2) :: f_res      

    f_res = A*ii8
  end function  timesi8z2_dz2

  !integer4 * dualz2
  elemental function timesi4z2_dz2(ii4,A) result(f_res)
    type(dualz2), intent(in) :: A
    integer(4), intent(in) :: ii4
    type(dualz2) :: f_res      

    f_res = A*ii4
  end function  timesi4z2_dz2

  !-dualz2 (unary)
  elemental function menos_dz2(A) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_res

    f_res%dualz1 = -A%dualz1
    f_res%f2 = -A%f2
  end function menos_dz2

  !dualz2 + dualz2
  elemental function  suma_dz2(A,B) result(f_res)
    type(dualz2), intent(in) :: A, B
    type(dualz2) :: f_res    

    f_res%dualz1 = A%dualz1 + B%dualz1
    f_res%f2 = A%f2 + B%f2
  end function suma_dz2

  !dualz2 - dualz2
  elemental function  resta_dz2(A,B) result(f_res)
    type(dualz2), intent(in) :: A, B
    type(dualz2) :: f_res    

    f_res%dualz1 = A%dualz1 - B%dualz1
    f_res%f2 = A%f2 - B%f2
  end function resta_dz2

  !dualz2 + dualz1
  elemental function sumaz2z1_dz2(A,Bz1) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz1), intent(in) :: Bz1
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 + Bz1
    f_res%f2 = A%f2      
  end function sumaz2z1_dz2

  !dualz2 - dualz1
  elemental function restaz2z1_dz2(A,Bz1) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz1), intent(in) :: Bz1
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 - Bz1
    f_res%f2 = A%f2
  end function restaz2z1_dz2

  !dualz2 + complex8
  elemental function sumaz2c8_dz2(A,c8) result(f_res)
    type(dualz2), intent(in) :: A
    complex(8), intent(in) :: c8
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 + c8
    f_res%f2 = A%f2
  end function  sumaz2c8_dz2

  !dualz2 - complex8
  elemental function restaz2c8_dz2(A,c8) result(f_res)
    type(dualz2), intent(in) :: A
    complex(8), intent(in) :: c8
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 - c8
    f_res%f2 = A%f2
  end function  restaz2c8_dz2

  !dualz2 + complex4
  elemental function sumaz2c4_dz2(A,c4) result(f_res)
    type(dualz2), intent(in) :: A
    complex(4), intent(in) :: c4
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 + c4
    f_res%f2 = A%f2
  end function  sumaz2c4_dz2

  !dualz2 - complex4
  elemental function restaz2c4_dz2(A,c4) result(f_res)
    type(dualz2), intent(in) :: A
    complex(4), intent(in) :: c4
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 - c4
    f_res%f2 = A%f2
  end function  restaz2c4_dz2

  !dualz2 + real8
  elemental function sumaz2r8_dz2(A,r8) result(f_res)
    type(dualz2), intent(in) :: A
    real(8), intent(in) :: r8
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 + r8
    f_res%f2 = A%f2
  end function  sumaz2r8_dz2

  !dualz2 - real8
  elemental function restaz2r8_dz2(A,r8) result(f_res)
    type(dualz2), intent(in) :: A
    real(8), intent(in) :: r8
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 - r8
    f_res%f2 = A%f2
  end function  restaz2r8_dz2

  !dualz2 + real4
  elemental function sumaz2r4_dz2(A,r4) result(f_res)
    type(dualz2), intent(in) :: A
    real(4), intent(in) :: r4
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 + r4
    f_res%f2 = A%f2
  end function  sumaz2r4_dz2

  !dualz2 - real4
  elemental function restaz2r4_dz2(A,r4) result(f_res)
    type(dualz2), intent(in) :: A
    real(4), intent(in) :: r4
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 - r4
    f_res%f2 = A%f2      
  end function  restaz2r4_dz2

  !dualz2 + integer8
  elemental function sumaz2i8_dz2(A,i8) result(f_res)
    type(dualz2), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 + i8
    f_res%f2 = A%f2      
  end function  sumaz2i8_dz2

  !dualz2 - integer8
  elemental function restaz2i8_dz2(A,i8) result(f_res)
    type(dualz2), intent(in) :: A
    integer(8), intent(in) :: i8
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 - i8
    f_res%f2 = A%f2      
  end function  restaz2i8_dz2

  !dualz2 + integer4
  elemental function sumaz2i4_dz2(A,i4) result(f_res)
    type(dualz2), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 + i4
    f_res%f2 = A%f2            
  end function  sumaz2i4_dz2

  !dualz2 - integer4
  elemental function restaz2i4_dz2(A,i4) result(f_res)
    type(dualz2), intent(in) :: A
    integer(4), intent(in) :: i4
    type(dualz2) :: f_res

    f_res%dualz1 = A%dualz1 - i4
    f_res%f2 = A%f2      
  end function  restaz2i4_dz2

  !dualz1 + dualz2
  elemental function sumaz1z2_dz2(Bzz1,A) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz1), intent(in) :: Bzz1
    type(dualz2) :: f_res

    f_res%dualz1 = Bzz1 + A%dualz1
    f_res%f2 = A%f2
  end function  sumaz1z2_dz2

  !dualz1 - dualz2
  elemental function restaz1z2_dz2(Bzz1,A) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz1), intent(in) :: Bzz1
    type(dualz2) :: f_res

    f_res%dualz1 = Bzz1 - A%dualz1
    f_res%f2 = -A%f2      
  end function  restaz1z2_dz2

  !complex8 + dualz2
  elemental function sumac8z2_dz2(cc8,A) result(f_res)
    type(dualz2), intent(in) :: A
    complex(8), intent(in) :: cc8
    type(dualz2) :: f_res

    f_res%dualz1 = cc8 + A%dualz1
    f_res%f2 = A%f2      
  end function  sumac8z2_dz2

  !complex8 - dualz2
  elemental function restac8z2_dz2(cc8,A) result(f_res)
    type(dualz2), intent(in) :: A
    complex(8), intent(in) :: cc8
    type(dualz2) :: f_res

    f_res%dualz1 = cc8 - A%dualz1
    f_res%f2 = -A%f2            
  end function  restac8z2_dz2

  !complex4 + dualz2
  elemental function sumac4z2_dz2(cc4,A) result(f_res)
    type(dualz2), intent(in) :: A
    complex(4), intent(in) :: cc4
    type(dualz2) :: f_res

    f_res%dualz1 = cc4 + A%dualz1
    f_res%f2 = A%f2      
  end function  sumac4z2_dz2

  !complex4 - dualz2
  elemental function restac4z2_dz2(cc4,A) result(f_res)
    type(dualz2), intent(in) :: A
    complex(4), intent(in) :: cc4
    type(dualz2) :: f_res

    f_res%dualz1 = cc4 - A%dualz1
    f_res%f2 = -A%f2
  end function  restac4z2_dz2

  !real8 + dualz2
  elemental function sumar8z2_dz2(rr8,A) result(f_res)
    type(dualz2), intent(in) :: A
    real(8), intent(in) :: rr8
    type(dualz2) :: f_res

    f_res%dualz1 = rr8 + A%dualz1
    f_res%f2 = A%f2      
  end function  sumar8z2_dz2

  !real8 - dualz2
  elemental function restar8z2_dz2(rr8,A) result(f_res)
    type(dualz2), intent(in) :: A
    real(8), intent(in) :: rr8
    type(dualz2) :: f_res

    f_res%dualz1 = rr8 - A%dualz1
    f_res%f2 = -A%f2      
  end function  restar8z2_dz2

  !real4 + dualz2
  elemental function sumar4z2_dz2(rr4,A) result(f_res)
    type(dualz2), intent(in) :: A
    real(4), intent(in) :: rr4
    type(dualz2) :: f_res

    f_res%dualz1 = rr4 + A%dualz1
    f_res%f2 = A%f2
  end function  sumar4z2_dz2

  !real4 - dualz2
  elemental function restar4z2_dz2(rr4,A) result(f_res)
    type(dualz2), intent(in) :: A
    real(4), intent(in) :: rr4
    type(dualz2) :: f_res

    f_res%dualz1 = rr4 - A%dualz1
    f_res%f2 = -A%f2
  end function  restar4z2_dz2

  !integer8 + dualz2
  elemental function sumai8z2_dz2(ii8,A) result(f_res)
    type(dualz2), intent(in) :: A
    integer(8), intent(in) :: ii8
    type(dualz2) :: f_res

    f_res%dualz1 = ii8 + A%dualz1
    f_res%f2 = A%f2  
  end function  sumai8z2_dz2

  !integer8 - dualz2
  elemental function restai8z2_dz2(ii8,A) result(f_res)
    type(dualz2), intent(in) :: A
    integer(8), intent(in) :: ii8
    type(dualz2) :: f_res

    f_res%dualz1 = ii8 - A%dualz1
    f_res%f2 = -A%f2      
  end function  restai8z2_dz2

  !integer4 + dualz2
  elemental function sumai4z2_dz2(ii4,A) result(f_res)
    type(dualz2), intent(in) :: A
    integer(4), intent(in) :: ii4
    type(dualz2) :: f_res

    f_res%dualz1 = ii4 + A%dualz1
    f_res%f2 = A%f2  
  end function  sumai4z2_dz2

  !integer4 - dualz2
  elemental function restai4z2_dz2(ii4,A) result(f_res)
    type(dualz2), intent(in) :: A
    integer(4), intent(in) :: ii4
    type(dualz2) :: f_res

    f_res%dualz1 = ii4 - A%dualz1
    f_res%f2 = -A%f2
  end function  restai4z2_dz2


  !Mathematical functions
  elemental function sin_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = sin(A%dualz1)
    f_r%f2 = -sin(g0)*g1**2 + cos(g0)*g2
  end function sin_dz2

  elemental function cos_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = cos(A%dualz1)
    f_r%f2 = -cos(g0)*g1**2 - sin(g0)*g2
  end function cos_dz2

  elemental function tan_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = tan(A%dualz1)
    f_r%f2 = g2/cos(g0)**2 + (2d0*tan(g0)/cos(g0)**2)*g1**2
  end function tan_dz2

  elemental function sinh_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = sinh(A%dualz1)
    f_r%f2 = g2*cosh(g0) + sinh(g0)*g1**2
  end function sinh_dz2

  elemental function cosh_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = cosh(A%dualz1)
    f_r%f2 = g2*sinh(g0) + cosh(g0)*g1**2
  end function cosh_dz2

  elemental function tanh_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = tanh(A%dualz1)
    f_r%f2 = g2/cosh(g0)**2 - 2d0*g1**2*tanh(g0)/cosh(g0)**2
  end function tanh_dz2

  elemental function exp_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = exp(A%dualz1)
    f_r%f2 = exp(g0)*g1**2 + exp(g0)*g2
  end function exp_dz2

  elemental function log_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = log(A%dualz1)
    f_r%f2 = g2/g0 - (g1/g0)**2
  end function log_dz2

  elemental function asin_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = asin(A%dualz1)
    f_r%f2 = g2/sqrt(1d0 - g0**2) + g0*g1**2/((1d0 - g0**2)**1.5d0)
  end function asin_dz2

  elemental function acos_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = acos(A%dualz1)
    f_r%f2 =  -g2/sqrt(1d0 - g0**2) - g0*g1**2/((1d0 - g0**2)**1.5d0)
  end function acos_dz2

  elemental function atan_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = atan(A%dualz1)
    f_r%f2 = g2/(1d0 + g0**2) - 2d0*g0*g1**2/(g0**2 + 1d0)**2
  end function atan_dz2

  elemental function asinh_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = asinh(A%dualz1)
    f_r%f2 = g2/sqrt(1d0 + g0**2) - (g0*g1**2)/(1d0 + g0**2)**1.5d0
  end function asinh_dz2

  elemental function acosh_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2, g0m, g0p, sg0m, sg0p

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    g0m  = g0 - 1d0
    g0p  = g0 + 1d0
    sg0m = sqrt(g0m) 
    sg0p = sqrt(g0p)

    f_r%dualz1 = acosh(A%dualz1)
    f_r%f2 =  g2/(sg0m*sg0p) - (g0*g1**2)/(g0m**1.5d0 * g0p**1.5d0)
  end function acosh_dz2

  elemental function atanh_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: g0, g1, g2

    g0 = A%f0
    g1 = A%f1
    g2 = A%f2

    f_r%dualz1 = atanh(A%dualz1)
    f_r%f2 =  2d0*g0*g1**2/(1d0 - g0**2)**2 + g2/(1d0 - g0**2)
  end function atanh_dz2

  elemental function atan2_dz2(A,B) result(f_r)
    type(dualz2), intent(in) :: A, B
    type(dualz2) :: f_r
    complex(8) :: A0, A1, A2, B0, B1, B2

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2

    f_r%dualz1 = atan2(A%dualz1, B%dualz1)
    f_r%f2 = (B0**2*(A2*B0 - 2d0*A1*B1) + A0**2*(A2*B0 + 2d0*A1*B1) -&
         A0**3*B2 - A0*B0*(2d0*A1**2 - 2d0*B1**2 + B0*B2))/(A0**2 +  &
         B0**2)**2
  end function atan2_dz2

  elemental function sqrt_dz2(A) result(f_r)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_r
    complex(8) :: A0, A1, A2

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2

    f_r%dualz1 = sqrt(A%dualz1)
    f_r%f2 = -A1**2/(4.d0*A0**1.5d0) + A2/(2.d0*sqrt(A0))
  end function sqrt_dz2

  !conjg
  !notice tat the conjugation operation is not, in general, differentia-
  !ble. The below definitions means (df)*, not d(f*), which is some
  !times useful
  elemental function conjg_dz2(g) result(f_r)
    type(dualz2), intent(in) :: g
    type(dualz2) :: f_r
    complex(8) :: g2    

    g2 = g%f2

    f_r%dualz1 = conjg(g%dualz1)
    f_r%f2 = conjg(g2)
  end function conjg_dz2

  !abs
  !The absolute value is not, in general, differentiable. In the below
  !function (df)* is used instead of d(f*) 
  elemental function abs_dz2(A) result(f_res)
    type(dualz2), intent(in) :: A
    type(dualz2) :: f_res

    f_res = sqrt(A*conjg(A))
  end function abs_dz2

  !|A> = |A0> + |A1> eps + |A2> eps2; |B> = ...
  !As in the conjg and abs functions, here, the dual parts are (df)* not
  !d(f*)
  function dot_product_dz2(A,B) result(f_r)
    type(dualz2), intent(in), dimension(:) :: A, B
    type(dualz2) :: f_r
    complex(8), dimension(size(A)) :: A0, A1, A2, B0, B1, B2

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2

    f_r%dualz1 = dot_product(A%dualz1,B%dualz1)
    f_r%f2 = dot_product(A0,B2) + 2*dot_product(A1,B1) + &
         dot_product(A2,B0)    
  end function dot_product_dz2

  !matmul
  function matmul22_dz2(A,B) result(f_r)
    type(dualz2), intent(in), dimension(:,:) :: A, B
    type(dualz2), dimension(size(A,1),size(B,2)) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A0, A1, A2
    complex(8), dimension(size(B,1),size(B,2)) :: B0, B1, B2

    A0 = A%f0
    A1 = A%f1
    A2 = A%f2

    B0 = B%f0
    B1 = B%f1
    B2 = B%f2

    f_r%dualz1 = matmul(A%dualz1,B%dualz1)
    f_r%f2 = matmul(A0,B2) + 2*matmul(A1,B1) + matmul(A2,B0) 
  end function matmul22_dz2

  function matmul22_dz2C(A,B) result(fr)
    type(dualz2), intent(in), dimension(:,:) :: A
    complex(8), intent(in), dimension(:,:) :: B
    type(dualz2), dimension(size(A,1),size(B,2)) :: fr
    type(dualz2), dimension(size(B,1),size(B,2)) :: Bdz2

    Bdz2 = B
    fr = matmul22_dz2(A,Bdz2)
  end function matmul22_dz2C

  function matmul22_dz2R(A,B) result(fr)
    type(dualz2), intent(in), dimension(:,:) :: A
    real(8), intent(in), dimension(:,:) :: B
    type(dualz2), dimension(size(A,1),size(B,2)) :: fr
    type(dualz2), dimension(size(B,1),size(B,2)) :: Bdz2

    Bdz2 = B
    fr = matmul22_dz2(A,Bdz2)
  end function matmul22_dz2R

  function matmul22_Cdz2(A,B) result(fr)
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz2), intent(in), dimension(:,:) :: B
    type(dualz2), dimension(size(A,1),size(B,2)) :: fr
    type(dualz2), dimension(size(A,1),size(A,2)) :: Adz2

    Adz2 = A
    fr = matmul22_dz2(Adz2,B)
  end function matmul22_Cdz2

  function matmul22_Rdz2(A,B) result(fr)
    real(8), intent(in), dimension(:,:) :: A
    type(dualz2), intent(in), dimension(:,:) :: B
    type(dualz2), dimension(size(A,1),size(B,2)) :: fr
    type(dualz2), dimension(size(A,1),size(A,2)) :: Adz2

    Adz2 = A
    fr = matmul22_dz2(Adz2,B)
  end function matmul22_Rdz2

  !dual_rank2*dual_rank1
  function matmul_dr21_dz2(A,x) result(f_r)
    type(dualz2), intent(in) :: A(:,:)
    type(dualz2), intent(in), dimension(size(A,2)) :: x
    type(dualz2), dimension(size(A,1)) :: f_r
    integer :: m, n

    m = size(A,1)
    n = size(A,2)
    f_r = reshape(matmul(A,reshape(x,[n,1])),[m])
  end function matmul_dr21_dz2

  function matmul_dr21_dz2C(A,x) result(fr)
    type(dualz2), intent(in), dimension(:,:) :: A
    complex(8), intent(in), dimension(size(A,2)) :: x
    type(dualz2), dimension(size(A,1)) :: fr
    type(dualz2), dimension(size(A,2)) :: xd

    xd = x
    fr = matmul_dr21_dz2(A,xd)
  end function matmul_dr21_dz2C

  function matmul_dr21_dz2R(A,x) result(fr)
    type(dualz2), intent(in), dimension(:,:) :: A
    real(8), intent(in), dimension(size(A,2)) :: x
    type(dualz2), dimension(size(A,1)) :: fr
    type(dualz2), dimension(size(A,2)) :: xd

    xd = x
    fr = matmul_dr21_dz2(A,xd)
  end function matmul_dr21_dz2R

  function matmul_dr21_Cdz2(A,x) result(fr)
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz2), intent(in), dimension(size(A,2)) :: x
    type(dualz2), dimension(size(A,1)) :: fr
    type(dualz2), dimension(size(A,1),size(A,2)) :: Adz2

    Adz2 = A
    fr = matmul_dr21_dz2(Adz2,x)
  end function matmul_dr21_Cdz2

  function matmul_dr21_Rdz2(A,x) result(fr)
    real(8), intent(in), dimension(:,:) :: A
    type(dualz2), intent(in), dimension(size(A,2)) :: x
    type(dualz2), dimension(size(A,1)) :: fr
    type(dualz2), dimension(size(A,1),size(A,2)) :: Adz2

    Adz2 = A
    fr = matmul_dr21_dz2(Adz2,x)
  end function matmul_dr21_Rdz2

  !in order to avoid ambiguity we must change x to xx
  !"A generic function must be able to distinguish its arguments by 
  !type AND by name"
  function matmul_dr12_dz2(xx,A) result(f_r)
    type(dualz2), intent(in) :: A(:,:)
    !type(dualz2), intent(in), dimension(size(A,1)) :: xx
    type(dualz2), intent(in), dimension(:) :: xx
    type(dualz2), dimension(size(A,2)) :: f_r
    integer :: m, n

    m = size(A,1)
    n = size(A,2)
    !if size(xx) /= m you may want to use some warning, the below code
    !reliess on the reshape fortran function which uses padding
    f_r = reshape(matmul(reshape(xx,[1,m]),A),[n])
  end function matmul_dr12_dz2

  function matmul_dr12_dz2C(xx,A) result(fr)
    type(dualz2), intent(in), dimension(:) :: xx
    complex(8), intent(in), dimension(:,:) :: A
    type(dualz2), dimension(size(A,2)) :: fr
    type(dualz2), dimension(size(A,1),size(A,2)) :: Adz2

    Adz2 = A
    fr =  matmul_dr12_dz2(xx,Adz2)
  end function matmul_dr12_dz2C

  function matmul_dr12_dz2R(xx,A) result(fr)
    type(dualz2), intent(in), dimension(:) :: xx
    real(8), intent(in), dimension(:,:) :: A
    type(dualz2), dimension(size(A,2)) :: fr
    type(dualz2), dimension(size(A,1),size(A,2)) :: Adz2

    Adz2 = A
    fr =  matmul_dr12_dz2(xx,Adz2)
  end function matmul_dr12_dz2R

  function matmul_dr12_Cdz2(xx,A) result(fr)
    complex(8), intent(in), dimension(:) :: xx
    type(dualz2), intent(in), dimension(:,:) :: A
    type(dualz2), dimension(size(A,2)) :: fr
    type(dualz2), dimension(size(xx)) :: xxdz2

    xxdz2 = xx
    fr =  matmul_dr12_dz2(xxdz2,A)
  end function matmul_dr12_Cdz2

  function matmul_dr12_Rdz2(xx,A) result(fr)
    real(8), intent(in), dimension(:) :: xx
    type(dualz2), intent(in), dimension(:,:) :: A
    type(dualz2), dimension(size(A,2)) :: fr
    type(dualz2), dimension(size(xx)) :: xxdz2

    xxdz2 = xx
    fr =  matmul_dr12_dz2(xxdz2,A)
  end function matmul_dr12_Rdz2

  !sum for rank2 arrays
  !the result is given as array of rank 1
  function sum2_dz2(A,k) result(f_r)
    type(dualz2), intent(in), dimension(:,:) :: A
    integer, intent(in) :: k
    type(dualz2), dimension(size(A,2/k)) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A2

    A2 = A%f2

    f_r%dualz1 = sum(A%dualz1,k)
    f_r%f2 = sum(A2,k)
  end function sum2_dz2

  !sum(A) rank2
  function sum0_dz2(A) result(f_r)
    type(dualz2), intent(in), dimension(:,:) :: A
    type(dualz2) :: f_r
    complex(8), dimension(size(A,1),size(A,2)) :: A2

    A2 = A%f2

    f_r%dualz1 = sum(A%dualz1)
    f_r%f2 = sum(A2)
  end function  sum0_dz2

  !sum(A) rank1
  function sum10_dz2(A) result(f_r)
    type(dualz2), intent(in), dimension(:) :: A
    type(dualz2) :: f_r    
    complex(8), dimension(size(A)) :: A2

    A2 = A%f2

    f_r%dualz1 = sum(A%dualz1)
    f_r%f2 = sum(A2)
  end function  sum10_dz2

  !-----------

  !product function for dualz1 vectors (arrays of rank 1)
  function prod1_dz2(x) result(f_r)
    type(dualz2), intent(in), dimension(:) :: x
    type(dualz2) :: f_r
    integer :: k

    f_r = 1
    do k = 1, size(x)
       f_r = f_r * x(k)
    end do
  end function prod1_dz2

  !product(x) for arrays of rank 2
  function prod0_dz2(x) result(f_r)
    type(dualz2), intent(in), dimension(:,:) :: x
    type(dualz2) :: f_r

    f_r = prod1_dz2(reshape(x,[size(x)]))
  end function prod0_dz2

  !product(A,k) for rank2 arrays
  !the result is given as array of rank 1
  function prod2_dz2(x,c) result(f_r)
    type(dualz2), intent(in), dimension(:,:) :: x
    integer, intent(in) :: c
    type(dualz2), dimension(size(x,2/c)) :: f_r
    type(dualz2), dimension(size(x,c)) :: xc
    integer :: k

    if(c==1) then
       do k = 1, size(x,2)
          xc = x(:,k)
          f_r(k) = prod1_dz2(xc)
       end do
    else if(c==2) then
       do k = 1, size(x,1)
          xc = x(k,:)
          f_r(k) = prod1_dz2(xc)
       end do
    else 
       stop 'c must be equal to 1 (2) to collapse rows (columns)'
    end if
  end function prod2_dz2

!=======================================================================
  !auxiliar functions
  elemental function f0part_dz2(A) result(fr)
    type(dualz2), intent(in) :: A
    complex(8) :: fr

    fr = A%f0
  end function f0part_dz2

  elemental function f1part_dz2(A) result(fr)
    type(dualz2), intent(in) :: A
    complex(8) :: fr
    
    fr = A%f1
  end function f1part_dz2

  elemental function f2part_dz2(A) result(fr)
    type(dualz2), intent(in) :: A
    complex(8) :: fr
    
    fr = A%f2
  end function f2part_dz2

!=====================================================================
  !auxiliar function to convert to dualz4
  !in principle this function is not necessary (the assignment operator
  !does the "same" job) but is provided for convenience
  !numeric (including dualz1, dualz2 and dualz3) to dualz4
  
  elemental function dz2todualz2(z) result(fr)
    type(dualz2), intent(in) :: z
    type(dualz2) :: fr

    fr = z
  end function dz2todualz2

  elemental function dz1todualz2(z) result(fr)
    type(dualz1), intent(in) :: z
    type(dualz2) :: fr

    fr = z
  end function dz1todualz2

  elemental function z8todualz2(z) result(fr)
    complex(8), intent(in) :: z
    type(dualz2) :: fr

    fr = z
  end function z8todualz2

  !z4 here is for complex(4)
  elemental function z4todualz2(z) result(fr)
    complex(4), intent(in) :: z
    type(dualz2) :: fr

    fr = z
  end function z4todualz2

  elemental function r8todualz2(z) result(fr)
    real(8), intent(in) :: z
    type(dualz2) :: fr

    fr = z
  end function r8todualz2

  elemental function r4todualz2(z) result(fr)
    real(4), intent(in) :: z
    type(dualz2) :: fr

    fr = z
  end function r4todualz2

  elemental function i8todualz2(z) result(fr)
    integer(8), intent(in) :: z
    type(dualz2) :: fr

    fr = z
  end function i8todualz2

  elemental function i4todualz2(z) result(fr)
    integer(4), intent(in) :: z
    type(dualz2) :: fr

    fr = z
  end function i4todualz2
  !=====================================================================  
  
  
end module dualz2_mod
