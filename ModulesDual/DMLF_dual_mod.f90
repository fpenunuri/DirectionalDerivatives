!Derivative multilinear forms of ranks 1,2,3,4
!Partial derivatives, orders 1,2,3,4
!
!F. Peñuñuri
!UADY, Merida Yucatan Mexico
!2022
module DMLF_dual_mod
  use dualz4_mod
  implicit none
  private

  !For efficiency reasons, it is recommended to use:
  !d4mlf(f,v,v,u,u,q,n) instead of d4mlf(f,v,u,v,u,q,n). Also
  !d4mlf(f,x,q,n) can be used instead of d4mlf(f,x,x,x,x,q,n)
  !this comment also applies to the d2mlf and d3mlf functions
  !f can be scalar or vector function
  !NOTE: if f is a vector function its dimension n, must be given.
  !Also, the vectors are coded as simple arrays of dimension (dim),
  !not as ket vectors--matrices of dimension (dim,1).

  public :: d1mlf, d2mlf, d3mlf, d4mlf
  public :: df1, df2, df3, df4

  !derivative multilinear forms
  !mn means Dm ---> Dn (dual m to dual n)
  !s means scalar function f:Dm ---> D
  !C is for complex, R is for real
  !Rank 1
  interface d1mlf
     module procedure d1q_mnC  
     module procedure d1q_mnR       
     module procedure d1q_sC  
     module procedure d1q_sR       
  end interface d1mlf

  !Rank2
  interface d2mlf
     module procedure d2mlf_mnC
     module procedure d2mlf_mnR
     module procedure d2q_mnC
     module procedure d2q_mnR
     module procedure d2mlf_sC
     module procedure d2mlf_sR
     module procedure d2q_sC
     module procedure d2q_sR    
  end interface d2mlf

  !Rank 3
  interface d3mlf
     module procedure d3mlf_mnC
     module procedure d3mlf_mnR
     module procedure d3q_mnC
     module procedure d3q_mnR       
     module procedure d3mlf_sC
     module procedure d3mlf_sR
     module procedure d3q_sC
     module procedure d3q_sR
  end interface d3mlf

  !Rank 4
  interface d4mlf
     module procedure d4mlf_mnC
     module procedure d4mlf_mnR
     module procedure d4q_mnC
     module procedure d4q_mnR
     module procedure d4mlf_sC
     module procedure d4mlf_sR
     module procedure d4q_sC
     module procedure d4q_sR
  end interface d4mlf

  !derivatives
  !Rank4
  interface df4
     module procedure df4mnC  
     module procedure df4mnR
     module procedure df4sC  
     module procedure df4sR
  end interface df4

  !Rank3
  interface df3
     module procedure df3mnC  
     module procedure df3mnR
     module procedure df3sC  
     module procedure df3sR
  end interface df3

  !Rank2
  interface df2
     module procedure df2mnC  
     module procedure df2mnR
     module procedure df2sC  
     module procedure df2sR
  end interface df2

  !Rank 1    
  interface df1
     module procedure df1mnC  
     module procedure df1mnCv 
     module procedure df1mnR
     module procedure df1mnRv 
     module procedure df1sC  
     module procedure df1sCv 
     module procedure df1sR
     module procedure df1sRv 
  end interface df1

  !interface for the function argument f:Dm ---> Dn
  abstract interface
     function dfunction_mn(r) result(fres)
       use dualz4_mod
       type(dualz4), intent(in), dimension(:) :: r !dimension m 
       type(dualz4), allocatable, dimension(:) :: fres !dimension n
     end function dfunction_mn
  end interface

  !interface for the function argument f:Dm --> D
  abstract interface
     function dfunction_s(r) result(fres)
       use dualz4_mod
       type(dualz4), intent(in), dimension(:) :: r !dimension m 
       type(dualz4) :: fres !scalar function
     end function dfunction_s
  end interface

contains

  !partial derivatives, vector case
  !rank 4
  function df4mnC(f,ijklv,q,n) result(fr)
    procedure(dfunction_mn) :: f
    integer, intent(in), dimension(4) :: ijklv
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    complex(8), dimension(size(q)) :: ei, ej, ek, el

    ei = 0; ei(ijklv(1)) = 1
    ej = 0; ej(ijklv(2)) = 1
    ek = 0; ek(ijklv(3)) = 1
    el = 0; el(ijklv(4)) = 1

    fr = d4mlf(f,ei,ej,ek,el,q,n)
  end function df4mnC

  function df4mnR(f,ijklv,q,n) result(fr)
    procedure(dfunction_mn) :: f
    integer, intent(in), dimension(4) :: ijklv
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    real(8), dimension(size(q)) :: ei, ej, ek, el

    ei = 0; ei(ijklv(1)) = 1
    ej = 0; ej(ijklv(2)) = 1
    ek = 0; ek(ijklv(3)) = 1
    el = 0; el(ijklv(4)) = 1

    fr = d4mlf(f,ei,ej,ek,el,q,n)
  end function df4mnR

  !scalar case
  function df4sC(f,ijklv,q) result(fr)
    procedure(dfunction_s) :: f
    integer, intent(in), dimension(4) :: ijklv
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr
    complex(8), dimension(size(q)) :: ei, ej, ek, el

    ei = 0; ei(ijklv(1)) = 1
    ej = 0; ej(ijklv(2)) = 1
    ek = 0; ek(ijklv(3)) = 1
    el = 0; el(ijklv(4)) = 1

    fr = d4mlf(f,ei,ej,ek,el,q)
  end function df4sC

  function df4sR(f,ijklv,q) result(fr)
    procedure(dfunction_s) :: f
    integer, intent(in), dimension(4) :: ijklv
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr
    real(8), dimension(size(q)) :: ei, ej, ek, el

    ei = 0; ei(ijklv(1)) = 1
    ej = 0; ej(ijklv(2)) = 1
    ek = 0; ek(ijklv(3)) = 1
    el = 0; el(ijklv(4)) = 1

    fr = d4mlf(f,ei,ej,ek,el,q)
  end function df4sR

  !rank 3
  function df3mnC(f,ijkv,q,n) result(fr)
    procedure(dfunction_mn) :: f
    integer, intent(in), dimension(3) :: ijkv
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    complex(8), dimension(size(q)) :: ei, ej, ek

    ei = 0; ei(ijkv(1)) = 1
    ej = 0; ej(ijkv(2)) = 1
    ek = 0; ek(ijkv(3)) = 1

    fr = d3mlf(f,ei,ej,ek,q,n)
  end function df3mnC

  function df3mnR(f,ijkv,q,n) result(fr)
    procedure(dfunction_mn) :: f
    integer, intent(in), dimension(3) :: ijkv
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    real(8), dimension(size(q)) :: ei, ej, ek

    ei = 0; ei(ijkv(1)) = 1
    ej = 0; ej(ijkv(2)) = 1
    ek = 0; ek(ijkv(3)) = 1

    fr = d3mlf(f,ei,ej,ek,q,n)
  end function df3mnR

  function df3sC(f,ijkv,q) result(fr)
    procedure(dfunction_s) :: f
    integer, intent(in), dimension(3) :: ijkv
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr
    complex(8), dimension(size(q)) :: ei, ej, ek

    ei = 0; ei(ijkv(1)) = 1
    ej = 0; ej(ijkv(2)) = 1
    ek = 0; ek(ijkv(3)) = 1

    fr = d3mlf(f,ei,ej,ek,q)
  end function df3sC

  function df3sR(f,ijkv,q) result(fr)
    procedure(dfunction_s) :: f
    integer, intent(in), dimension(3) :: ijkv
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr
    real(8), dimension(size(q)) :: ei, ej, ek

    ei = 0; ei(ijkv(1)) = 1
    ej = 0; ej(ijkv(2)) = 1
    ek = 0; ek(ijkv(3)) = 1

    fr = d3mlf(f,ei,ej,ek,q)
  end function df3sR

  !rank 2
  function df2mnC(f,ijv,q,n) result(fr)
    procedure(dfunction_mn) :: f
    integer, intent(in), dimension(2) :: ijv
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    complex(8), dimension(size(q)) :: ei, ej

    ei = 0; ei(ijv(1)) = 1
    ej = 0; ej(ijv(2)) = 1

    fr = d2mlf(f,ei,ej,q,n)
  end function df2mnC

  function df2mnR(f,ijv,q,n) result(fr)
    procedure(dfunction_mn) :: f
    integer, intent(in), dimension(2) :: ijv
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    real(8), dimension(size(q)) :: ei, ej

    ei = 0; ei(ijv(1)) = 1
    ej = 0; ej(ijv(2)) = 1

    fr = d2mlf(f,ei,ej,q,n)
  end function df2mnR

  function df2sC(f,ijv,q) result(fr)
    procedure(dfunction_s) :: f
    integer, intent(in), dimension(2) :: ijv
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr
    complex(8), dimension(size(q)) :: ei, ej

    ei = 0; ei(ijv(1)) = 1
    ej = 0; ej(ijv(2)) = 1

    fr = d2mlf(f,ei,ej,q)
  end function df2sC

  function df2sR(f,ijv,q) result(fr)
    procedure(dfunction_s) :: f
    integer, intent(in), dimension(2) :: ijv
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr
    real(8), dimension(size(q)) :: ei, ej

    ei = 0; ei(ijv(1)) = 1
    ej = 0; ej(ijv(2)) = 1

    fr = d2mlf(f,ei,ej,q)
  end function df2sR

  !rank 1
  function df1mnRv(f,iv,q,n) result(fr)
    procedure (dfunction_mn) :: f
    integer, intent(in), dimension(1) :: iv
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    real(8), dimension(size(q)) :: ei

    ei = 0; ei(iv(1)) = 1
    fr = d1mlf(f,ei,q,n)
  end function df1mnRv

  function df1mnR(f,i,q,n) result(fr)
    procedure (dfunction_mn) :: f
    integer, intent(in) :: i, n
    real(8), intent(in), dimension(:) :: q
    real(8), dimension(n) :: fr
    real(8), dimension(size(q)) :: ei

    ei=0; ei(i) = 1
    fr = d1mlf(f,ei,q,n)
  end function df1mnR

  function df1mnCv(f,iv,q,n) result(fr)
    procedure (dfunction_mn) :: f
    integer, intent(in), dimension(1) :: iv
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    complex(8), dimension(size(q)) :: ei

    ei=0; ei(iv(1)) = 1
    fr = d1mlf(f,ei,q,n)
  end function df1mnCv

  function df1mnC(f,i,q,n) result(fr)
    procedure (dfunction_mn) :: f
    integer, intent(in) :: i, n
    complex(8), intent(in), dimension(:) :: q
    complex(8), dimension(n) :: fr
    complex(8), dimension(size(q)) :: ei

    ei=0; ei(i) = 1
    fr = d1mlf(f,ei,q,n)
  end function df1mnC

  function df1sRv(f,iv,q) result(fr)
    procedure (dfunction_s) :: f
    integer, intent(in), dimension(1) :: iv
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr
    real(8), dimension(size(q)) :: ei

    ei=0; ei(iv(1)) = 1
    fr = d1mlf(f,ei,q)
  end function df1sRv

  function df1sR(f,i,q) result(fr)
    procedure (dfunction_s) :: f
    integer, intent(in) :: i
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr
    real(8), dimension(size(q)) :: ei

    ei=0; ei(i) = 1
    fr = d1mlf(f,ei,q)
  end function df1sR

  function df1sCv(f,iv,q) result(fr)
    procedure (dfunction_s) :: f
    integer, intent(in), dimension(1) :: iv
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr
    complex(8), dimension(size(q)) :: ei

    ei=0; ei(iv(1)) = 1
    fr = d1mlf(f,ei,q)
  end function df1sCv

  function df1sC(f,i,q) result(fr)
    procedure (dfunction_s) :: f
    integer, intent(in) :: i
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr
    complex(8), dimension(size(q)) :: ei

    ei=0; ei(i) = 1
    fr = d1mlf(f,ei,q)
  end function df1sC

  !======================== multilinear forms ========================
  !scalar case
  !rank 1
  !f:Dm-->D
  function d1q_sR(f, v, q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: v 
    real(8), intent(in), dimension(:) :: q      
    real(8) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f1
  end function d1q_sR

  function d1q_sC(f, v, q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: v
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f1
  end function d1q_sC

  !rank 2
  function d2mlf_sC(f,x,y,q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: x, y
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr

    if(all(x == y)) then
       fr = d2q_sC(f,x,q)
    else
       fr = 0.5*(d2q_sC(f,x + y,q) - d2q_sC(f,x,q) - d2q_sC(f,y,q))
    end if
  end function d2mlf_sC

  function d2mlf_sR(f,x,y,q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: x, y
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr

    if(all(x == y)) then
       fr = d2q_sR(f,x,q)
    else
       fr = 0.5*(d2q_sR(f,x + y,q) - d2q_sR(f,x,q) - d2q_sR(f,y,q))
    end if
  end function d2mlf_sR

  function d2q_sC(f, v, q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: v, q
    complex(8) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f2
  end function d2q_sC

  function d2q_sR(f, v, q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: v, q
    real(8) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f2
  end function d2q_sR

  !rank 3
  function d3mlf_sR(f,x,y,z,q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: x, y, z
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr

    if(all(x == y) .and. all(y == z)) then
       fr = d3q_sR(f,x,q)
    elseif(all(x == y))then
       fr = d3qxxz_sR(f, x, z, q)
    else         
       fr = 0.5*(d3qxxz_sR(f,x + y,z,q) - d3qxxz_sR(f,x,z,q) - &
            d3qxxz_sR(f,y,z,q))
    end if
  end function d3mlf_sR

  function d3mlf_sC(f,x,y,z,q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: x, y, z
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr

    if(all(x == y) .and. all(y == z)) then
       fr = d3q_sC(f,x,q)
    elseif(all(x == y))then
       fr = d3qxxz_sC(f, x, z, q)
    else         
       fr = 0.5*(d3qxxz_sC(f,x + y,z,q) - d3qxxz_sC(f,x,z,q) - &
            d3qxxz_sC(f,y,z,q))
    end if
  end function d3mlf_sC

  function d3qxxz_sR(f, x, z, q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: x, z
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr

    fr = (d3q_sR(f,z+x,q) + d3q_sR(f,z-x,q) - 2*d3q_sR(f,z,q))/6
  end function d3qxxz_sR

  function d3qxxz_sC(f, x, z, q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: x, z
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr

    fr = (d3q_sC(f,z+x,q) + d3q_sC(f,z-x,q) - 2*d3q_sC(f,z,q))/6
  end function d3qxxz_sC

  function d3q_sR(f, v, q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: v, q
    real(8) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f3
  end function d3q_sR

  function d3q_sC(f, v, q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: v, q
    complex(8) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f3
  end function d3q_sC

  !rank 4
  function d4mlf_sC(f,x,y,z,w,q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: x, y, z, w
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr

    if(all(x == y) .and. all(y == z) .and. all(z == w)) then
       fr = d4q_sC(f,x,q)
    elseif(all(x == y) .and. all(z == w))then
       fr = d4qxxzz_sC(f,x,z,q)
    elseif(all(x == y)) then
       fr = d4qxxzw_sC(f,x,z,w,q)
    else         
       fr = (d4qxxzw_sC(f,x+y,z,w,q) - d4qxxzw_sC(f,x,z,w,q) - &
            d4qxxzw_sC(f,y,z,w,q))/2
    end if
  end function d4mlf_sC

  function d4mlf_sR(f,x,y,z,w,q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: x, y, z, w
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr

    if(all(x == y) .and. all(y == z) .and. all(z == w)) then
       fr = d4q_sR(f,x,q)
    elseif(all(x == y) .and. all(z == w))then
       fr = d4qxxzz_sR(f,x,z,q)
    elseif(all(x == y)) then
       fr = d4qxxzw_sR(f,x,z,w,q)
    else         
       fr = (d4qxxzw_sR(f,x+y,z,w,q) - d4qxxzw_sR(f,x,z,w,q) - &
            d4qxxzw_sR(f,y,z,w,q))/2
    end if
  end function d4mlf_sR

  function d4qxxzw_sC(f, x, z, w, q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: x, z, w
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr 

    fr = (d4qxxzz_sC(f,x,z + w,q) - d4qxxzz_sC(f,x,z,q) - &
         d4qxxzz_sC(f,x,w,q))/2
  end function d4qxxzw_sC

  function d4qxxzw_sR(f, x, z, w, q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: x, z, w
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr 

    fr = (d4qxxzz_sR(f,x,z + w,q) - d4qxxzz_sR(f,x,z,q) - &
         d4qxxzz_sR(f,x,w,q))/2
  end function d4qxxzw_sR

  function d4qxxzz_sC(f, x, z, q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: x, z
    complex(8), intent(in), dimension(:) :: q
    complex(8) :: fr
    complex(8) :: d4xzp, d4xzm

    d4xzp = d4q_sC(f, x + z, q)
    d4xzm = d4q_sC(f, x - z, q)

    fr = (d4xzp + d4xzm - 2*d4q_sC(f,x,q) - 2*d4q_sC(f,z,q))/12
  end function d4qxxzz_sC

  function d4qxxzz_sR(f, x, z, q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: x, z
    real(8), intent(in), dimension(:) :: q
    real(8) :: fr
    real(8) :: d4xzp, d4xzm

    d4xzp = d4q_sR(f, x + z, q)
    d4xzm = d4q_sR(f, x - z, q)

    fr = (d4xzp + d4xzm - 2*d4q_sR(f,x,q) - 2*d4q_sR(f,z,q))/12
  end function d4qxxzz_sR

  function d4q_sC(f, v, q) result(fr)
    procedure (dfunction_s) :: f
    complex(8), intent(in), dimension(:) :: v, q
    complex(8) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f4
  end function d4q_sC

  function d4q_sR(f, v, q) result(fr)
    procedure (dfunction_s) :: f
    real(8), intent(in), dimension(:) :: v, q
    real(8) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f4
  end function d4q_sR

  !vector case
  !rank 1
  !The argument n is the dimension of f. This could be avoided
  !using n = size(f(q)). However in order to avoid this function
  !evaluation, n is requested as argument     
  function d1q_mnC(f, v, q, n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: v
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4), dimension(n) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f1
  end function d1q_mnC

  function d1q_mnR(f, v, q, n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: v 
    real(8), intent(in), dimension(:) :: q 
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4), dimension(n) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f1
  end function d1q_mnR

  !rank 2
  function d2mlf_mnC(f,x,y,q,n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: x, y
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr

    if(all(x == y)) then
       fr = d2q_mnC(f,x,q,n)
    else
       fr = 0.5*(d2q_mnC(f,x + y,q,n) - d2q_mnC(f,x,q,n) - &
            d2q_mnC(f,y,q,n))
    end if
  end function d2mlf_mnC

  function d2mlf_mnR(f,x,y,q,n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: x, y
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr

    if(all(x == y)) then
       fr = d2q_mnR(f,x,q,n)
    else
       fr = 0.5*(d2q_mnR(f,x + y,q,n) - d2q_mnR(f,x,q,n) - &
            d2q_mnR(f,y,q,n))
    end if
  end function d2mlf_mnR

  function d2q_mnC(f, v, q, n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4), dimension(n) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f2
  end function d2q_mnC

  function d2q_mnR(f, v, q, n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4), dimension(n) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f2
  end function d2q_mnR

  !rank 3
  !para mejorar la eficiencia, si dos argumentos son iguales, se
  !recomienda ponerlos al prncipio (v,u,v)-->(v,v,u)
  function d3mlf_mnR(f,x,y,z,q,n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: x, y, z
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr

    if(all(x == y) .and. all(y == z)) then
       fr = d3q_mnR(f,x,q,n)
    elseif(all(x == y))then
       fr = d3qxxz_mnR(f, x, z, q, n)
    else         
       fr = 0.5*(d3qxxz_mnR(f,x + y,z,q,n) - d3qxxz_mnR(f,x,z,q,n) - &
            d3qxxz_mnR(f,y,z,q,n))
    end if
  end function d3mlf_mnR

  function d3mlf_mnC(f,x,y,z,q,n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: x, y, z
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr

    if(all(x == y) .and. all(y == z)) then
       fr = d3q_mnC(f,x,q,n)
    elseif(all(x == y))then
       fr = d3qxxz_mnC(f, x, z, q, n)
    else         
       fr = 0.5*(d3qxxz_mnC(f,x + y,z,q,n) - d3qxxz_mnC(f,x,z,q,n) - &
            d3qxxz_mnC(f,y,z,q,n))
    end if
  end function d3mlf_mnC

  function d3qxxz_mnR(f, x, z, q, n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: x, z
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr

    fr = (d3q_mnR(f,z+x,q,n) + d3q_mnR(f,z-x,q,n) - &
         2*d3q_mnR(f,z,q,n))/6
  end function d3qxxz_mnR

  function d3qxxz_mnC(f, x, z, q, n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: x, z
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr

    fr = (d3q_mnC(f,z+x,q,n) + d3q_mnC(f,z-x,q,n) - &
         2*d3q_mnC(f,z,q,n))/6
  end function d3qxxz_mnC

  function d3q_mnR(f, v, q, n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4), dimension(n) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f3
  end function d3q_mnR

  function d3q_mnC(f, v, q, n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4), dimension(n) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f3
  end function d3q_mnC

  !rank 4
  function d4mlf_mnC(f,x,y,z,w,q,n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: x, y, z, w
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr

    if(all(x == y) .and. all(y == z) .and. all(z == w)) then
       fr = d4q_mnC(f,x,q,n)
    elseif(all(x == y) .and. all(z == w))then
       fr = d4qxxzz_mnC(f,x,z,q,n)
    elseif(all(x == y)) then
       fr = d4qxxzw_mnC(f,x,z,w,q,n)
    else         
       fr = (d4qxxzw_mnC(f,x+y,z,w,q,n) - d4qxxzw_mnC(f,x,z,w,q,n) - &
            d4qxxzw_mnC(f,y,z,w,q,n))/2
    end if
  end function d4mlf_mnC

  function d4mlf_mnR(f,x,y,z,w,q,n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: x, y, z, w
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr

    if(all(x == y) .and. all(y == z) .and. all(z == w)) then
       fr = d4q_mnR(f,x,q,n)
    elseif(all(x == y) .and. all(z == w))then
       fr = d4qxxzz_mnR(f,x,z,q,n)
    elseif(all(x == y)) then
       fr = d4qxxzw_mnR(f,x,z,w,q,n)
    else         
       fr = (d4qxxzw_mnR(f,x+y,z,w,q,n) - d4qxxzw_mnR(f,x,z,w,q,n) - &
            d4qxxzw_mnR(f,y,z,w,q,n))/2
    end if
  end function d4mlf_mnR

  function d4qxxzw_mnC(f, x, z, w, q, n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: x, z, w
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr 

    fr = (d4qxxzz_mnC(f,x,z + w,q,n) - d4qxxzz_mnC(f,x,z,q,n) - &
         d4qxxzz_mnC(f,x,w,q,n))/2
  end function d4qxxzw_mnC

  function d4qxxzw_mnR(f, x, z, w, q, n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: x, z, w
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr 

    fr = (d4qxxzz_mnR(f,x,z + w,q,n) - d4qxxzz_mnR(f,x,z,q,n) - &
         d4qxxzz_mnR(f,x,w,q,n))/2
  end function d4qxxzw_mnR

  function d4qxxzz_mnC(f, x, z, q, n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: x, z
    complex(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    complex(8), dimension(n) :: d4xzp, d4xzm

    d4xzp = d4q_mnC(f, x + z, q, n)
    d4xzm = d4q_mnC(f, x - z, q, n)

    fr = (d4xzp + d4xzm - 2*d4q_mnC(f,x,q,n) - 2*d4q_mnC(f,z,q,n))/12
  end function d4qxxzz_mnC

  function d4qxxzz_mnR(f, x, z, q, n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: x, z
    real(8), intent(in), dimension(:) :: q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    real(8), dimension(n) :: d4xzp, d4xzm

    d4xzp = d4q_mnR(f, x + z, q, n)
    d4xzm = d4q_mnR(f, x - z, q, n)

    fr = (d4xzp + d4xzm - 2*d4q_mnR(f,x,q,n) - 2*d4q_mnR(f,z,q,n))/12
  end function d4qxxzz_mnR

  function d4q_mnC(f, v, q, n) result(fr)
    procedure (dfunction_mn) :: f
    complex(8), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    complex(8), dimension(n) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4), dimension(n) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f4
  end function d4q_mnC

  function d4q_mnR(f, v, q, n) result(fr)
    procedure (dfunction_mn) :: f
    real(8), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    type(dualz4), parameter :: eps1 = dualz4(0,1,0,0,0)
    type(dualz4), dimension(n) :: fqd

    fqd = f(q + eps1*v)
    fr = fqd%f4
  end function d4q_mnR
end module DMLF_dual_mod
