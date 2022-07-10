!FMFD method to compute until the fourth order partial derivative
!to any order of approximation
!
!F. Peñuñuri
!UADY, Merida Yucatan Mexico
!2022
module finite_diff4_mod
  use iso_fortran_env,  only: prc4 => real32, prc8 => real64, &
       prc16 => real128
  implicit none

  private
  public :: df4_diff, df3_diff, df2_diff, df1_diff, deriv_diff

  interface df4_diff
     module procedure df4_diff128
     module procedure df4_diff64
     module procedure df4_diff32
  end interface df4_diff

  interface df3_diff
     module procedure df3_diff128
     module procedure df3_diff64
     module procedure df3_diff32
  end interface df3_diff

  interface df2_diff
     module procedure df2_diff128
     module procedure df2_diff64
     module procedure df2_diff32
  end interface df2_diff

  interface deriv_diff
     module procedure deriv_diff32
     module procedure deriv_diff64
     module procedure deriv_diff128
  end interface deriv_diff

  interface df1_diff
     module procedure df1s128
     module procedure df1v128
     module procedure df1s64
     module procedure df1v64
     module procedure df1s32
     module procedure df1v32
  end interface df1_diff

  !vector arg
  abstract interface
     function funv32(r) result(fres)
       use iso_fortran_env,  only: prc4 => real32
       real(prc4), intent(in), dimension(:) :: r
       real(prc4) :: fres
     end function funv32
  end interface

  abstract interface
     function funv64(r) result(fres)
       use iso_fortran_env,  only: prc8 => real64
       real(prc8), intent(in), dimension(:) :: r
       real(prc8) :: fres
     end function funv64
  end interface

  abstract interface
     function funv128(r) result(fres)
       use iso_fortran_env,  only: prc16 => real128
       real(prc16), intent(in), dimension(:) :: r
       real(prc16) :: fres
     end function funv128
  end interface

  !scalar arg
  abstract interface
     function funs32(r) result(fres)
       use iso_fortran_env,  only: prc4 => real32
       real(prc4), intent(in) :: r
       real(prc4) :: fres
     end function funs32
  end interface

  abstract interface
     function funs64(r) result(fres)
       use iso_fortran_env,  only: prc8 => real64
       real(prc8), intent(in) :: r
       real(prc8) :: fres
     end function funs64
  end interface

  abstract interface
     function funs128(r) result(fres)
       use iso_fortran_env,  only: prc16 => real128
       real(prc16), intent(in) :: r
       real(prc16) :: fres
     end function funs128
  end interface

contains
  
  !auxiliar function to construct the sets of n's. See Eq...x
  !rank=1,2,3,4
  !oa=1,2,...  (values greater than 8 are not recommended)
  function auxiliarF(rank,oa) result(fr)
    integer, intent(in) :: rank, oa
    integer, dimension(rank) :: fr
    integer :: min, res, k

    min = (oa + rank - 1)/rank
    res = mod(oa + rank - 1, rank)

    do k=1,rank - res
       fr(k) = min
    end do

    do k=rank - res + 1,rank
       fr(k) = min + 1
    end do
  end function auxiliarF

  !derivatives
  function df4_diff128(f,iv4,x0,h,N) result(fr)
    procedure(funv128) :: f
    real(prc16), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(4) :: iv4
    integer, intent(in) :: N
    real(prc16), intent(in) :: h
    real(prc16) :: fr
    integer ::  n1,n2,n3,n4
    integer, dimension(4) :: setn

    setn = auxiliarF(4,N)
    n1=setn(1); n2=setn(2); n3=setn(3); n4=setn(4) 

    fr = df4_diff128X(f,iv4,x0,h,n1,n2,n3,n4)    
  end function df4_diff128

  function df4_diff64(f,iv4,x0,h,N) result(fr)
    procedure(funv64) :: f
    real(prc8), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(4) :: iv4
    integer, intent(in) :: N
    real(prc8), intent(in) :: h
    real(prc8) :: fr
    integer ::  n1,n2,n3,n4
    integer, dimension(4) :: setn

    setn = auxiliarF(4,N)
    n1=setn(1); n2=setn(2); n3=setn(3); n4=setn(4) 

    fr = df4_diff64X(f,iv4,x0,h,n1,n2,n3,n4)    
  end function df4_diff64

  function df4_diff32(f,iv4,x0,h,N) result(fr)
    procedure(funv32) :: f
    real(prc4), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(4) :: iv4
    integer, intent(in) :: N
    real(prc4), intent(in) :: h
    real(prc4) :: fr
    integer ::  n1,n2,n3,n4
    integer, dimension(4) :: setn

    setn = auxiliarF(4,N)
    n1=setn(1); n2=setn(2); n3=setn(3); n4=setn(4) 

    fr = df4_diff32X(f,iv4,x0,h,n1,n2,n3,n4)    
  end function df4_diff32

  function df3_diff128(f,ijk,x0,h,N) result(fr)
    procedure(funv128) :: f
    real(prc16), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(3) :: ijk
    integer, intent(in) :: N
    real(prc16), intent(in) :: h
    real(prc16) :: fr  
    integer :: n1,n2,n3
    integer, dimension(3) :: setn

    setn = auxiliarF(3,N)
    n1=setn(1); n2=setn(2); n3=setn(3)

    fr = df3_diff128X(f,ijk,x0,h,n1,n2,n3)
  end function df3_diff128

  function df3_diff64(f,ijk,x0,h,N) result(fr)
    procedure(funv64) :: f
    real(prc8), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(3) :: ijk
    integer, intent(in) :: N
    real(prc8), intent(in) :: h
    real(prc8) :: fr  
    integer :: n1,n2,n3
    integer, dimension(3) :: setn

    setn = auxiliarF(3,N)
    n1=setn(1); n2=setn(2); n3=setn(3)

    fr = df3_diff64X(f,ijk,x0,h,n1,n2,n3)
  end function df3_diff64

  function df3_diff32(f,ijk,x0,h,N) result(fr)
    procedure(funv32) :: f
    real(prc4), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(3) :: ijk
    integer, intent(in) :: N
    real(prc4), intent(in) :: h
    real(prc4) :: fr  
    integer :: n1,n2,n3
    integer, dimension(3) :: setn

    setn = auxiliarF(3,N)
    n1=setn(1); n2=setn(2); n3=setn(3)

    fr = df3_diff32X(f,ijk,x0,h,n1,n2,n3)
  end function df3_diff32

  function df2_diff128(f,ij,x0,h,N) result(fr)
    procedure(funv128) :: f
    real(prc16), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(2) :: ij
    integer, intent(in) :: N
    real(prc16), intent(in) :: h
    real(prc16) :: fr  
    integer :: n1,n2
    integer, dimension(2) :: setn

    setn = auxiliarF(2,N)
    n1=setn(1); n2=setn(2)

    fr = df2_diff128X(f,ij,x0,h,n1,n2)
  end function df2_diff128

  function df2_diff64(f,ij,x0,h,N) result(fr)
    procedure(funv64) :: f
    real(prc8), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(2) :: ij
    integer, intent(in) :: N
    real(prc8), intent(in) :: h
    real(prc8) :: fr  
    integer :: n1,n2
    integer, dimension(2) :: setn

    setn = auxiliarF(2,N)
    n1=setn(1); n2=setn(2)

    fr = df2_diff64X(f,ij,x0,h,n1,n2)
  end function df2_diff64

  function df2_diff32(f,ij,x0,h,N) result(fr)
    procedure(funv32) :: f
    real(prc4), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(2) :: ij
    integer, intent(in) :: N
    real(prc4), intent(in) :: h
    real(prc4) :: fr  
    integer :: n1,n2
    integer, dimension(2) :: setn

    setn = auxiliarF(2,N)
    n1=setn(1); n2=setn(2)

    fr = df2_diff32X(f,ij,x0,h,n1,n2)
  end function df2_diff32

  !dfxijkq = df4_diff(f,iv4,x0,h,n1,n2,n3,n4)
  !dfxijkq is an approximation to (D^4/Dxijkq)f
  !f:Rm ---> R
  !x0=[x01,x02,..x0m] is the evaluating point
  !i,j,k,q=1,2,..m;  (i=1 means x1, etc)
  !h=steep size
  !n1,n2,n3,n4 = order of the approximation
  ![n1,n2,n3,n4] = [2, 3, 3, 3] or any set of 4 natural numbers with
  !n1+n2+n3+n4 = 11, will result in an approximation to order 8
  !in principle there is no limit to the order but there would be
  !lost of precision for higher orders; we recommend order 4 or 8
  !see the paper for more details
  !f:Rm ---> R
  function df4_diff128X(f,iv4,x0,h,n1,n2,n3,n4) result(fr)
    procedure(funv128) :: f
    real(prc16), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(4) :: iv4
    integer, intent(in) :: n1,n2,n3,n4
    real(prc16), intent(in) :: h
    real(prc16) :: fr
    real(prc16), dimension(size(x0)) :: v
    integer :: i1, i2, i3, i4

    i1 = iv4(1); i2 = iv4(2); i3 = iv4(3); i4 = iv4(4)
    fr = deriv_diff128(faux,x0(i4),h,n1)

  contains
    function faux(x) result(f3)
      real(prc16), intent(in) :: x
      real(prc16) :: f3

      v = x0
      v(i4) = x
      f3 = df3_diff128X(f, [i1,i2,i3], v, h, n2,n3,n4)
    end function faux
  end function df4_diff128X

  function df4_diff64X(f,iv4,x0,h,n1,n2,n3,n4) result(fr)
    procedure(funv64) :: f
    real(prc8), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(4) :: iv4
    integer, intent(in) :: n1,n2,n3,n4
    real(prc8), intent(in) :: h
    real(prc8) :: fr
    real(prc8), dimension(size(x0)) :: v
    integer :: i1, i2, i3, i4

    i1 = iv4(1); i2 = iv4(2); i3 = iv4(3); i4 = iv4(4)
    fr = deriv_diff64(faux,x0(i4),h,n1)

  contains
    function faux(x) result(f3)
      real(prc8), intent(in) :: x
      real(prc8) :: f3

      v = x0
      v(i4) = x
      f3 = df3_diff64X(f, [i1,i2,i3], v, h, n2,n3,n4)
    end function faux
  end function df4_diff64X

  function df4_diff32X(f,iv4,x0,h,n1,n2,n3,n4) result(fr)
    procedure(funv32) :: f
    real(prc4), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(4) :: iv4
    integer, intent(in) :: n1,n2,n3,n4
    real(prc4), intent(in) :: h
    real(prc4) :: fr
    real(prc4), dimension(size(x0)) :: v
    integer :: i1, i2, i3, i4

    i1 = iv4(1); i2 = iv4(2); i3 = iv4(3); i4 = iv4(4)
    fr = deriv_diff32(faux,x0(i4),h,n1)

  contains
    function faux(x) result(f3)
      real(prc4), intent(in) :: x
      real(prc4) :: f3

      v = x0
      v(i4) = x
      f3 = df3_diff32X(f, [i1,i2,i3], v, h, n2,n3,n4)
    end function faux
  end function df4_diff32X

  !dfxijk = df3(f,iv4,x0,h,n1,n2,n3)
  !dfxijkq is an approximation to (D^3/Dxijk)f
  !f:Rm ---> R
  function df3_diff128X(f,ijk,x0,h,n1,n2,n3) result(fr)
    procedure(funv128) :: f
    real(prc16), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(3) :: ijk
    integer, intent(in) :: n1,n2,n3
    real(prc16), intent(in) :: h
    real(prc16) :: fr
    real(prc16), dimension(size(x0)) :: v
    integer :: i, j, k

    i = ijk(1); j = ijk(2); k = ijk(3)

    fr = deriv_diff128(faux,x0(k),h,n1)

  contains
    function faux(x) result(f2)
      real(prc16), intent(in) :: x
      real(prc16) :: f2

      v = x0
      v(k) = x
      f2 = df2_diff128X(f, [i,j], v, h, n2,n3)
    end function faux
  end function df3_diff128X

  function df3_diff64X(f,ijk,x0,h,n1,n2,n3) result(fr)
    procedure(funv64) :: f
    real(prc8), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(3) :: ijk
    integer, intent(in) :: n1,n2,n3
    real(prc8), intent(in) :: h
    real(prc8) :: fr
    real(prc8), dimension(size(x0)) :: v
    integer :: i, j, k

    i = ijk(1); j = ijk(2); k = ijk(3)

    fr = deriv_diff64(faux,x0(k),h,n1)

  contains
    function faux(x) result(f2)
      real(prc8), intent(in) :: x
      real(prc8) :: f2

      v = x0
      v(k) = x
      f2 = df2_diff64X(f, [i,j], v, h, n2,n3)
    end function faux
  end function df3_diff64X

  function df3_diff32X(f,ijk,x0,h,n1,n2,n3) result(fr)
    procedure(funv32) :: f
    real(prc4), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(3) :: ijk
    integer, intent(in) :: n1,n2,n3
    real(prc4), intent(in) :: h
    real(prc4) :: fr
    real(prc4), dimension(size(x0)) :: v
    integer :: i, j, k

    i = ijk(1); j = ijk(2); k = ijk(3)

    fr = deriv_diff32(faux,x0(k),h,n1)

  contains
    function faux(x) result(f2)
      real(prc4), intent(in) :: x
      real(prc4) :: f2

      v = x0
      v(k) = x
      f2 = df2_diff32X(f, [i,j], v, h, n2,n3)
    end function faux
  end function df3_diff32X

  !dfxij = df(f,iv4,x0,h,n1,n2)
  !dfxij is an approximation to (D^2/Dxij)f
  !f:Rm ---> R
  function df2_diff128X(f,ij,x0,h,n1,n2) result(fr)
    procedure(funv128) :: f
    real(prc16), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(2) :: ij
    integer, intent(in) :: n1, n2
    real(prc16), intent(in) :: h
    real(prc16) :: fr
    real(prc16), dimension(size(x0)) :: v
    integer :: i, j

    i = ij(1); j = ij(2)    
    fr = deriv_diff128(faux,x0(j),h,n1)

  contains
    function faux(x) result(f1)
      real(prc16), intent(in) :: x
      real(prc16) :: f1

      v = x0
      v(j) = x
      f1 = df1s128(f, i, v, h, n2)      
    end function faux
  end function df2_diff128X

  function df2_diff64X(f,ij,x0,h,n1,n2) result(fr)
    procedure(funv64) :: f
    real(prc8), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(2) :: ij
    integer, intent(in) :: n1, n2
    real(prc8), intent(in) :: h
    real(prc8) :: fr
    real(prc8), dimension(size(x0)) :: v
    integer :: i, j

    i = ij(1); j = ij(2)    
    fr = deriv_diff64(faux,x0(j),h,n1)

  contains
    function faux(x) result(f1)
      real(prc8), intent(in) :: x
      real(prc8) :: f1

      v = x0
      v(j) = x
      f1 = df1s64(f, i, v, h, n2)      
    end function faux
  end function df2_diff64X

  function df2_diff32X(f,ij,x0,h,n1,n2) result(fr)
    procedure(funv32) :: f
    real(prc4), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(2) :: ij
    integer, intent(in) :: n1, n2
    real(prc4), intent(in) :: h
    real(prc4) :: fr
    real(prc4), dimension(size(x0)) :: v
    integer :: i, j

    i = ij(1); j = ij(2)    
    fr = deriv_diff32(faux,x0(j),h,n1)

  contains
    function faux(x) result(f1)
      real(prc4), intent(in) :: x
      real(prc4) :: f1

      v = x0
      v(j) = x
      f1 = df1s32(f, i, v, h, n2)      
    end function faux
  end function df2_diff32X

  function df1v128(f,iv,x0,h,n) result(fr)
    procedure(funv128) :: f
    real(prc16), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(1) :: iv
    real(prc16), intent(in) :: h
    integer, intent(in) :: n
    real(prc16) :: fr

    fr = df1s128(f,iv(1),x0,h,n)    
  end function df1v128

  function df1v64(f,iv,x0,h,n) result(fr)
    procedure(funv64) :: f
    real(prc8), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(1) :: iv
    real(prc8), intent(in) :: h
    integer, intent(in) :: n
    real(prc8) :: fr

    fr = df1s64(f,iv(1),x0,h,n)    
  end function df1v64

  function df1v32(f,iv,x0,h,n) result(fr)
    procedure(funv32) :: f
    real(prc4), intent(in), dimension(:) :: x0
    integer, intent(in), dimension(1) :: iv
    real(prc4), intent(in) :: h
    integer, intent(in) :: n
    real(prc4) :: fr

    fr = df1s32(f,iv(1),x0,h,n)    
  end function df1v32

  !dfxi = df1s(f,i,x0,h,n)
  !f:Rm-->R
  !x0=[x01,x02,..x0m]
  !i=1,2,..m;  (i=1 means x1, etc)
  !h=steep size
  !n=order of the approximation  "O(h)"
  !dfxi is te partial derivative of f with respect to xi
  !f:Rm-->R
  function df1s128(f,i,x0,h,n) result(fr)
    procedure(funv128) :: f
    real(prc16), intent(in), dimension(:) :: x0
    integer, intent(in) :: i, n
    real(prc16), intent(in) :: h
    real(prc16) :: fr
    real(prc16), dimension(size(x0)) :: v 

    fr = deriv_diff128(faux,x0(i),h,n)

  contains
    function faux(x) result(fri)
      real(prc16), intent(in) :: x
      real(prc16) :: fri

      v = x0
      v(i) = x
      fri = f(v)
    end function faux
  end function df1s128

  function df1s64(f,i,x0,h,n) result(fr)
    procedure(funv64) :: f
    real(prc8), intent(in), dimension(:) :: x0
    integer, intent(in) :: i, n
    real(prc8), intent(in) :: h
    real(prc8) :: fr
    real(prc8), dimension(size(x0)) :: v 

    fr = deriv_diff64(faux,x0(i),h,n)

  contains
    function faux(x) result(fri)
      real(prc8), intent(in) :: x
      real(prc8) :: fri

      v = x0
      v(i) = x
      fri = f(v)
    end function faux
  end function df1s64

  !32bits
  function df1s32(f,i,x0,h,n) result(fr)
    procedure(funv32) :: f
    real(prc4), intent(in), dimension(:) :: x0
    integer, intent(in) :: i, n
    real(prc4), intent(in) :: h
    real(prc4) :: fr
    real(prc4), dimension(size(x0)) :: v 

    fr = deriv_diff32(faux,x0(i),h,n)

  contains
    function faux(x) result(fri)
      real(prc4), intent(in) :: x
      real(prc4) :: fri

      v = x0
      v(i) = x
      fri = f(v)
    end function faux
  end function df1s32

  !DX = deriv (F, X0, H, O)
  !DX is the first order derivative of function F
  !F is a function f:R-->R
  !X0 is the the evaluating point
  !H is the steep size
  !O is the order of the approximation
  function deriv_diff128(f, x0, h, n) result (fr)
    procedure(funs128) :: f
    real(prc16), intent(in) :: x0
    real(prc16), intent(in) :: h
    integer, intent(in) :: n
    real(prc16) :: fr
    integer :: k

    fr = 0
    do k=1,n
       fr = fr + (-1)**(k+1)*Dop128(k, h, f, x0)/(h*k)
    end do
  end function deriv_diff128

  function deriv_diff64(f, x0, h, n) result (fr)
    procedure(funs64) :: f
    real(prc8), intent(in) :: x0
    real(prc8), intent(in) :: h
    integer, intent(in) :: n
    real(prc8) :: fr
    integer :: k

    fr = 0
    do k=1,n
       fr = fr + (-1)**(k+1)*Dop64(k, h, f, x0)/(h*k)
    end do
  end function deriv_diff64

  !32bits
  function deriv_diff32(f, x0, h, n) result (fr)
    procedure(funs32) :: f
    real(prc4), intent(in) :: x0
    real(prc4), intent(in) :: h
    integer, intent(in) :: n
    real(prc4) :: fr
    integer :: k

    fr = 0
    do k=1,n
       fr = fr + (-1)**(k+1)*Dop32(k, h, f, x0)/(h*k)
    end do
  end function deriv_diff32


  !Delta**n=Delta^n
  function Dop128(n, h, f, x0) result(fr)
    procedure(funs128) :: f
    integer, intent(in) :: n
    real(prc16), intent(in) :: x0, h
    real(prc16) :: fr
    integer :: k
    real(prc4) :: nr, kr !there is no need for higher precision

    nr = n
    fr = 0
    do k = 0, n
       kr = k
       fr = fr + (-1)**(n-k)*Cnr(nr,kr)*f(x0 + h*k)
    end do
  end function Dop128

  function Dop64(n, h, f, x0) result(fr)
    procedure(funs64) :: f
    integer, intent(in) :: n
    real(prc8), intent(in) :: x0, h
    real(prc8) :: fr
    integer :: k
    real(prc4) :: nr, kr !there is no need for higher precision

    nr = n
    fr = 0
    do k = 0, n
       kr = k
       fr = fr + (-1)**(n-k)*Cnr(nr,kr)*f(x0 + h*k)
    end do
  end function Dop64

  !32bits
  function Dop32(n, h, f, x0) result(fr)
    procedure(funs32) :: f
    integer, intent(in) :: n
    real(prc4), intent(in) :: x0, h
    real(prc4) :: fr
    integer :: k
    real(prc4) :: nr, kr !there is no need for higher precision

    nr = n
    fr = 0
    do k = 0, n
       kr = k
       fr = fr + (-1)**(n-k)*Cnr(nr,kr)*f(x0 + h*k)
    end do
  end function Dop32

  function Cnr(n,r) result(fr)
    real(prc4), intent(in) :: n,r !there is no need for higher precision
    real(prc4) :: fr

    fr = gamma(n + 1)/(gamma(n - r + 1)*gamma(r + 1))
  end function Cnr
end module finite_diff4_mod
