
C                      ASSIGNMENT 4

c Van der Waals equation. Bisection and Newton-Raphson.

      program pract4
      implicit none
      real*8 p1,t,a,v0,p2,eps,xarrel,x1,x2,dfun,lim
      real*8 arrel1,arrel2
      real*8 x(80),funci(80),dfunci(80),v4(5),nr(5),temp(2)
      integer j,niter
      external p1,p2,dfun
      COMMON/DADES/t      

c The Van der Waals equation to describe the transition liquid-gas can be written as
c P=(8T/3V-1)-(3/V^2) with T,P anv V in reduced units. V has to be > 1/3 to have physical sense.

      nr=0.d0
      a=(1.d0/3)+0.1d0
      eps=1.d-12 !precission
      v4=(/0.36d0,0.48d0,0.99d0,1.18d0,1.56d0/) ! values of volume
      temp=(/0.92d0,0.98d0/) !values of temperature
      t=temp(1)

c                  FORMATS

700   format(3(e20.12,5x))
800   format(10x,a,23x,a,19x,a)
100   format(2(e20.12,5x))
200   format(10x,a,23x,a)
300   format(e20.12,5x,i10)
900   format(10x,a,15x,a)
400   format(5x,a,18x,a,19x,a)


c The estability of the isotherms can be analyzed using the first derivative of P(V) 

c Generates a table with 80 equally spaced points of V (V € [1/3+0.1,4]) of the isotherm (V,P) where T=0.92
c Computes the derivative P'(V) and writes in the file of resuts the values of V,P'(V) and P(V)

      open(1,file='P4-18P-res.dat',status='unknown')
      write(1,800) 'v','dP(v)','P(v)' 
      call vectors(80,p1,a,4.d0,funci,x)
      call defunT1(80,x,funci,dfunci)
      do j=1,79
            v0=x(j)
            write(1,700) v0,dfunci(j),funci(j)
      end do
      write(1,*) ' '
      write(1,*) ' ' 
      close(1) 

c The third degree polynomial POL(V,T)=4TV^3-9V^2+6V-1 appears when looking for the metastability 
c limits of the isotherms described by the Van der Waals equation under the critical temperature.

c Generates a table with 80 equally spaced points of V (V € [1/3+0.1,2]) of the isotherm (V,P) where T=0.92

      open(2,file='P4-18P-res.dat',status='unknown',position='append')
      write(2,200) 'v','P(v)'
      call vectors(80,p2,a,2.d0,funci,x)
      do j=1,79
            v0=x(j)
            write(2,100) v0,funci(j)
      end do
      write(2,*) ' '
      write(2,*) ' ' 
      close(2)

C Finds the zeros of the polynomial using the bisection algorithm

      open(3,file='P4-18P-res.dat',status='unknown',position='append')
      write(3,*) 'Zeros using bisection:'
      write(3,*) ''
      write(3,900) 'Zero','Iterations'
      x1=0.6d0
      x2=0.8d0
      call biT1(x1,x2,eps,niter,p2,xarrel)
      write(3,300) xarrel,niter
      arrel1=xarrel
      x1=1.4d0
      x2=1.6d0
      call biT1(x1,x2,eps,niter,p2,xarrel)
      arrel2=xarrel
      write(3,300) xarrel,niter
      write(3,*) ' '
      write(3,*) ' ' 
      close(3)

C Finds the zeros of the polynomial using the Newton-Raphson algorithm for the temperature T=0.98
c starts from five diferent values of V0 = 0.36,0.48,0.99 and 1.56

      t=temp(2)

      open(6,file='P4-18P-res.dat',status='unknown',position='append')
      write(6,*) 'Zeros using Newton-Raphson:'
      close(6)
      do j=1,5
      	lim=1.d0/3.d0
      	if (v4(j).gt.lim) call nrT1(v4(j),eps,niter,p2,dfun,xarrel)
c      	write(*,*) 'Arrel:',xarrel,'Num iteracions:',niter, 'xo',x0(j)
            nr(j)=xarrel
      end do

      stop
      end program


c              SUBROUTINES

c Subroutine that computes the derivative of the function

      subroutine defunT1(ndat,x,funci,dfunci)
      implicit none
      real*8 x(80),funci(80),dfunci(80),h
      integer ndat,k
      
      h=dble(x(2)-x(1))
c      write(*,*) h
      dfunci(1)=(funci(2)-funci(1))/h
      do k=2,ndat-1
            dfunci(k)=(funci(k+1)-funci(k-1))/(2.d0*h)
      end do
      dfunci(ndat)=(funci(ndat)-funci(ndat-1))/h

      return
      end subroutine

c Subroutine that creates the vectors x and funci with the values of the variable and the function equally spaced

      subroutine vectors(ndat,fun,a,b,funci,x)
      implicit none
      real*8 a,b,fun,funci(80),x(80),h,v
      integer ndat,j

      h=(b-a)/dble(ndat-1)
      v=a
      j=1
      do while (v.le.b)
            x(j)=v
            funci(j)=fun(v)
            v=v+h
            j=j+1
      end do

      return
      end subroutine

c Subroutine that applies the bisection method

      subroutine biT1(a,b,eps,niter,fun,xarrel)
      implicit none
      real*8 a,b,eps,fun,xarrel,x,error
      integer niter,n

      error=1000.d0
      n=0

      if (fun(a)*fun(b).ge.0.d0) then 
            write(*,*) 'The sign of the function does not change, there is no zero'
      else
            do while ((error.gt.eps).and.(fun(x).ne.0.d0))
      	     x=(a+b)/2.d0
      	     n=n+1
      	     if (fun(a)*fun(x).gt.0.d0) a=x
      	     if (fun(a)*fun(x).lt.0.d0) b=x
      	     error=dabs(a-b)
            end do
            niter=n
            xarrel=x
      end if

      return
      end subroutine

c Subroutine that applies the Newton-Raphson method

      subroutine nrT1(x0,eps,niter,fun,dfun,xarrel)
      implicit none
      real*8 x0,eps,fun,dfun,xarrel,error,x1,x
      integer niter,n

600   format(e20.12,5x,e20.12,5x,i10)
500   format(9x,a,23x,a,15x,a)

      open(4,file='P4-18P-res.dat',status='unknown',position='append') 
c      write(2,500) 'Convergencia Newton-Raphson x0:',x0
c      write(2,*) ' '
      write(4,500) 'x0','Zero','Iteration'
      error=1000.d0
      x1=x0
      n=0
     
      do while(error.gt.eps)
      	x=x1-(fun(x1)/dfun(x1))
      	error=dabs(x-x1)
      	n=n+1
      	write(4,600) x0,x,n
      	x1=x
      end do
      write(4,*) ' '
      write(4,*) ' '
      
      niter=n
      xarrel=x

      close(4)

      return 
      end subroutine


c                 FUNCTIONS      

c Polynomial 1

      real*8 function p1(v)
      implicit none
      real*8 t,v,funcio
      COMMON/DADES/t
      
      funcio=(8.d0*t/(3.d0*v-1.d0))-(3.d0/(v**2))
      p1=funcio

      return 
      end function

c Polynomial 2

      real*8 function p2(v)
      implicit none
      real*8 t,v,funcio
      COMMON/DADES/t
      
      funcio=4.d0*t*(v**3)-9.d0*(v**2)+6.d0*v-1.d0
      p2=funcio

      return 
      end function

c Derivative polynomial 2

      real*8 function dfun(v)
      implicit none
      real*8 v,t,funcio
      COMMON/DADES/t

      funcio=12.d0*t*(v**2)-18.d0*v+6.d0
      dfun=funcio

      return 
      end function
