!----------------------------------------------------------------
! (c) Copyright, 2018 by the Regents of the University of California.
! Besselclass: Bessel function class in Math Function module of FUNCTION layer.
! Version: 2.0
! Author: Ji Qiang
! Description: This class defines the Bessel functions.
! Comments:
!----------------------------------------------------------------

 module Besselclass

      contains



 double precision FUNCTION I0(x)
    implicit none
    double precision::i6,i5,i4,i3,i2,i1
    double precision,intent(in)::x
      double precision:: elap
i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720
i5=(x/2)**5*(1+(x/2)**2/6+(x/2)**4/84)/120
i4=10*i5/x+i6
i3=8*i4/x+i5
i2=6*i3/x+i4
i1=4*i2/x+i3
i0=2*i1/x+i2

  end function I0
  end module Besselclass
