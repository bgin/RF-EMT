! *************************************************************************
!
!      *************
!      *           *
!      *    uni    *
!      *           *
!      *************
!
!
      real*8 function uni(jd);
!
!***begin prologue  uni
!***date written   810915
!***revision date  830805
!***category no.  l6a21
!***keywords  random numbers, uniform random numbers
!***author    blue, james, scientific computing division, nbs
!             kahaner, david, scientific computing division, nbs
!             marsaglia, george, computer science dept., wash state univ
!
!***purpose  this routine generates quasi uniform random numbers on [0,1)
!             and can be used on any computer with which allows integers
!             at least as large as 32767.
!***description
!
!       this routine generates quasi uniform random numbers on the interval
!       [0,1).  it can be used with any computer which allows
!       integers at least as large as 32767.
!
!
!   use
!       first time....
!                   z = uni(jd)
!                     here jd is any  n o n - z e r o  integer.
!                     this causes initialization of the program
!                     and the first random number to be returned as z.
!       subsequent times...
!                   z = uni(0)
!                     causes the next random number to be returned as z.
!
!
!..................................................................
!   note: users who wish to transport this program from one computer
!         to another should read the following information.....
!
!   machine dependencies...
!      mdig = a lower bound on the number of binary digits available
!              for representing integers, including the sign bit.
!              this value must be at least 16, but may be increased
!              in line with remark a below.
!
!   remarks...
!     a. this program can be used in two ways:
!        (1) to obtain repeatable results on different computers,
!            set 'mdig' to the smallest of its values on each, or,
!        (2) to allow the longest sequence of random numbers to be
!            generated without cycling (repeating) set 'mdig' to the
!            largest possible value.
!     b. the sequence of numbers generated depends on the initial
!          input 'jd' as well as the value of 'mdig'.
!          if mdig=16 one should find that
!            the first evaluation
!              z=uni(305) gives z=.027832881...
!            the second evaluation
!              z=uni(0) gives   z=.56102176...
!            the third evaluation
!              z=uni(0) gives   z=.41456343...
!            the thousandth evaluation
!              z=uni(0) gives   z=.19797357...
!
!***references  marsaglia g., "comments on the perfect uniform random
!                 number generator", unpublished notes, wash s. u.
!***end prologue  uni
!
!
      implicit real*8 (a-h,o-z);
      integer m(17);
      integer mdata;

      save i,j,m,m1,m2;

      data m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),    &
     &     m(12),m(13),m(14),m(15),m(16),m(17)                          &
     &                   / 30788,23052,2053,19346,10646,19427,23975,    &
     &                     19049,10949,19693,29746,26748,2796,23890,    &
     &                     29168,31924,16499 /;
      data m1,m2,i,j / 32767,256,5,17 /;
      data mdata /32/;
      !______________________________________________

!***first executable statement  uni

      if (jd /= 0) then;
         !***  fill
         mdig =mdata;
         !***  be sure that mdig at least 16...
         if (mdig < 16) stop 'mdig < 16';
         m1= 2**(mdig-2) + (2**(mdig-2)-1);
         m2 = 2**(mdig/2);
         jseed = min0(iabs(jd),m1);
         if ( mod(jseed,2) == 0 ) jseed=jseed-1;
         k0 =mod(9069,m2);
         k1 = 9069/m2;
         j0 = mod(jseed,m2);
         j1 = jseed/m2;

         do i=1,17;
           jseed = j0*k0;
           j1 = mod(jseed/m2+j0*k1+j1*k0,m2/2);
           j0 = mod(jseed,m2);
           m(i) = j0+m2*j1;
         enddo;

         i=5;
         j=17;
      endif;

!***  begin main loop here
      k =m(i)-m(j);
      if (k < 0) k =k+m1;
      m(j) =k;
      i =i-1;
      if (i == 0) i=17;
      j =j-1;
      if (j == 0) j=17;
      uni =dfloat(k)/dfloat(m1);
      uni =max(uni,1.d-20);   !protect against uni =0

      return;
      end function uni;
