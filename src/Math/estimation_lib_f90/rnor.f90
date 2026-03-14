! ************************************************************************
!
!        **************
!        *            *
!        *    rnor    *
!        *            *
!        **************
!
      real*8 function rnor(jd);
!
!***begin prologue  rnor
!***date written   810915
!***revision date  830805
!***category no.  l6a14
!***keywords  random numbers, uniform random numbers
!***author    kahaner, david, scientific computing division, nbs
!             marsaglia, george, computer science dept., wash state univ
!
!***purpose  generates quasi normal random numbers, with mean zero and
!             unit standard deviation, and can be used with any computer
!             with integers at least as large as 32767.
!***description
!
!       rnor generates quasi normal random numbers with zero mean and
!       unit standard deviation.
!       it can be used with any computer with integers at least as
!       large as 32767.
!
!
!   use
!       first time....
!                   z = rnor(jd)
!                     here jd is any  n o n - z e r o  integer.
!                     this causes initialization of the program
!                     and the first random number to be returned as z.
!       subsequent times...
!                   z = rnor(0)
!                     causes the next random number to be returned as z.
!
!.....................................................................
!
!    note: users who wish to transport this program to other
!           computers should read the following ....
!
!   machine dependencies...
!      mdig = a lower bound on the number of binary digits available
!              for representing integers, including the sign bit.
!              this must be at least 16, but can be increased in
!              line with remark a below.
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
!              z=rnor(87) gives  z=-.40079207...
!            the second evaluation
!              z=rnor(0) gives   z=-1.8728870...
!            the third evaluation
!              z=rnor(0) gives   z=1.8216004...
!            the fourth evaluation
!              z=rnor(0) gives   z=.69410355...
!            the thousandth evaluation
!              z=rnor(0) gives   z=.96782424...
!
!***references  marsaglia & tsang, "a fast, easily implemented
!                 method for sampling from decreasing or
!                 symmetric unimodal density functions", to be
!                 published in siam j sisc 1983.
!***routines called  i1mach,xerror
!***end prologue  rnor
!
!
      implicit real*8 (a-h,o-z);
      real*8  v(65),w(65);
      integer m(17);
      integer mdata;
      save i1,j1,m,m1,m2,rmax;
      data aa,b,c,rmax/12.37586,.4878992,12.67706,3.0518509e-5/;
      data c1,c2,pc,xn/.9689279,1.301198,.1958303e-1,2.776994/;
      data v/ .3409450, .4573146, .5397793, .6062427, .6631691          &
     &, .7136975, .7596125, .8020356, .8417227, .8792102, .9148948      &
     &, .9490791, .9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917  &
     &,1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635 &
     &,1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929 &
     &,1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454 &
     &,1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422 &
     &,1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166 &
     &,1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713 &
     &,2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117 &
     &,2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943, 2.7769943/;
      data w/   .10405134e-04, .13956560e-04, .16473259e-04,            &
     & .18501623e-04, .20238931e-04, .21780983e-04, .23182241e-04,      &
     & .24476931e-04, .25688121e-04, .26832186e-04, .27921226e-04,      &
     & .28964480e-04, .29969191e-04, .30941168e-04, .31885160e-04,      &
     & .32805121e-04, .33704388e-04, .34585827e-04, .35451919e-04,      &
     & .36304851e-04, .37146564e-04, .37978808e-04, .38803170e-04,      &
     & .39621114e-04, .40433997e-04, .41243096e-04, .42049621e-04,      &
     & .42854734e-04, .43659562e-04, .44465208e-04, .45272764e-04,      &
     & .46083321e-04, .46897980e-04, .47717864e-04, .48544128e-04,      &
     & .49377973e-04, .50220656e-04, .51073504e-04, .51937936e-04,      &
     & .52815471e-04, .53707761e-04, .54616606e-04, .55543990e-04,      &
     & .56492112e-04, .57463436e-04, .58460740e-04, .59487185e-04,      &
     & .60546402e-04, .61642600e-04, .62780711e-04, .63966581e-04,      &
     & .65207221e-04, .66511165e-04, .67888959e-04, .69353880e-04,      &
     & .70922996e-04, .72618816e-04, .74471933e-04, .76525519e-04,      &
     & .78843526e-04, .81526890e-04, .84749727e-04,                     &
     & .84749727e-04, .84749727e-04, .84749727e-04/;
      data m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),    &
     &     m(12),m(13),m(14),m(15),m(16),m(17)                          &
     &                   / 30788,23052,2053,19346,10646,19427,23975,    &
     &                     19049,10949,19693,29746,26748,2796,23890,    &
     &                     29168,31924,16499 /;
      data m1,m2,i1,j1 / 32767,256,5,17 /;
      data mdata /32/;
      !__________________________________________


!***first executable statement  rnor
      if (jd /= 0) then;
         !**  fill
!        mdig=i1mach(8)+1
         mdig =mdata;
         !**  be sure that mdig at least 16...
         if (mdig < 16) stop 'mdig < 16';
         m1 = 2**(mdig-2) + (2**(mdig-2)-1);
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
         j1=17;
         i1=5;
         rmax = 1.d0/float(m1);

         !***  seed uniform (0,1) generator.  (just a dummy call)
         rnor =uni(jd);
         do i=1,65;
           w(i)=rmax*v(i);
         enddo;
      endif;

! fast part...
      i =m(i1)-m(j1);
      if (i < 0) i=i+m1;
      m(j1) =i;
      i1 =i1-1;
      if (i1 == 0) i1=17;
      j1 =j1-1;
      if (j1 == 0) j1=17;
      j =mod(i,64)+1;
      rnor =i*w(j+1);
      if ( ((i/m2)/2)*2 == (i/m2) ) rnor=-rnor;
      if (abs(rnor) <= v(j)) return;

! slow part; aa is a*f(0)
      x =(abs(rnor)-v(j))/(v(j+1)-v(j));
      y =uni(0);
      s =x+y;
      if (s > c2) go to 11;
      if (s <= c1) return;

      if (y > c-aa*exp(-.5d0*(b-b*x)**2)) go to 11;
      if (exp(-.5d0*v(j+1)**2)+y*pc/v(j+1) <= exp(-.5d0*rnor**2)) return;

! tail part; 3.855849 is .5*xn**2
   22 s =xn-log(uni(0))/xn;
      if (3.855849d0+log(uni(0))-xn*s > -.5d0*s**2) go to 22;
      rnor =sign(s,rnor);
      return;

   11 rnor =sign(b-b*x,rnor);
      return;
      end;
