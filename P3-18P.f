
C                      ASSIGNMENT 3

c  Orbit of Kohoutek's comet. Numerical integration.
c  Equation of the orbit f(x)= b*sqrt(1-(((x+4a)^2)/a^2))

      program pract3

      implicit none
      real*8 ykohoutek,xini,xfin,at,as,h,a,pi,exact,errt,errs,b
      real*8 exact2,xi2,xf2,a2t,a2s,h2,h4
      integer k,n
      external ykohoutek

c      write(*,*) ykohoutek(-3.5d0)
      a=508.633 !in units of x10^6 km
      b=429.074 !in units of x10^6 km
c area of the orbit: four times the integral of f(x) between -4a and-3a
      xini=-4.d0*a
      xfin=-3.d0*a
      pi=dacos(-1.0d0)
      exact=pi*a*b !exact value of the area

c                       FORMATS

100   format(9x,a,22x,a,22x,a,22x,a,17x,a,17x,a,17x,a)
200   format(7(e21.14,3x))
300   format(9x,a,22x,a,22x,a,19x,a,17x,a,17x,a,17x,a)


c Computes the area using the trapezium rule and Simpson's rule
c Computes the error in the value of the area for the diferent lengths of the interval h

      open(1,file='P3-18P-res1.dat')
      write(1,100) 'h','At','As','Error t','Error s','h^2','h^4'
c the number of intervals is 2^m with m = 2,....20
      do k=2,20
      	n=2**k
        h=(xfin-xini)/dble(n)
        call trapezis(n,ykohoutek,xini,xfin,at)
        call simpson(xini,xfin,k,ykohoutek,as)
        at=4.d0*at
        as=4.d0*as
        errt=dabs(exact-at)
        errs=dabs(exact-as)
        h2=h**2 !estimation of the error for the trapezium rule
        h4=h**4 !estimation of the error for the Simpson's rule
        write(1,200) h,at,as,errt,errs,h2,h4
      end do
      close(1)

c We consider now the area given by the integral of f(x) between -4a and -7a/2
c and compute again the area using the two rules and the error

      xi2=-4.d0*a
      xf2=-7.d0*a/2.d0
      exact2=a*b*(3.d0*dsqrt(3.d0)+2.d0*pi)/24.d0 ! new exact value of the area

      open(2,file='P3-18P-res2.dat')
      write(2,300) 'h','A2t','A2s','Error t','Error s','h^2','h^4'
      do k=2,20
      	n=2**k
        h=(xf2-xi2)/dble(n)
        call trapezis(n,ykohoutek,xi2,xf2,a2t)
        call simpson(xi2,xf2,k,ykohoutek,a2s)
        errt=dabs(exact2-a2t)
        errs=dabs(exact2-a2s)
        h2=h**2
        h4=h**4
        write(2,200) h,a2t,a2s,errt,errs,h2,h4
      end do
      close(2)

      stop
      end program


c                         FUNCTIONS

c Function that computes the comet's orbit
      real*8 function ykohoutek(x)
      implicit none
      real*8 a,b,funcio,x,div

      a=508.633
      b=429.074

      div=((x+4.d0*a)/a)**2
      funcio=b*dsqrt(1.d0-div)
      ykohoutek=funcio

      return
      end function


c                         SUBROUTINES

c subroutine that applies the trapezium rule
c in: nint=number of intervals, (a,b)=interval of integration, fcn=function to integrate
c out: integral computed by the trapezium rule

      subroutine trapezis(ninter,fcn,a,b,integral)
      implicit none
      real*8 h,x,a,b,fcn,integral,int
      integer ninter,i

      int=fcn(a)+fcn(b)
      h=(b-a)/dble(ninter)
      do i=1,ninter-1
        x=a+dble(i)*h
        int=int+2.d0*fcn(x)
      end do
      integral=int*h/2.d0

      return
      end subroutine

c subroutine that applies the Simpson's rule
c in: m=exponent -> number of intevals 2^m, (a,b)=interval of integration, fcn=function to integrate
c out: integral computed by the Simpson's rule

      subroutine simpson(a,b,m,fcn,integral)
      implicit none
      real*8 h,x,a,b,fcn,i,integral
      integer j,m,n

      i=fcn(a)+fcn(b)
      n=2**m
      h=(b-a)/dble(n)
      do j=1,n-1
        x=a+dble(j)*h
        if (mod(j,2).eq.0) i=i+2.d0*fcn(x)
        if (mod(j,2).ne.0) i=i+4.d0*fcn(x)
      end do  
      integral=i*h/3.d0

      return
      end subroutine
