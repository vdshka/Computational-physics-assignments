
C                      ASSIGNMENT 6

c      Monte Carlo methods
c   **the assignment specifies to do everything inside the subroutines**


      program pre6
      implicit none
      real*8 fun,f1,f2,pi,f3,f4,f5,g,l,psi
      external fun,f1,f2,f3,f4,f5,g,psi
      COMMON/DADES/pi,l
      pi=dacos(-1.0d0)
      l=16.d0 !length in microns      
     
      call montecarlop6()
      call mcarloMD()


      stop
      end program


c              SUBROUTINES

c Subroutine that applies the acceptance-rejection method

      subroutine subaccepta(ndat,xnums,a,b,M,fun)
      implicit none
      integer ndat,iseed,i
      real*8 xnums(1:1000000),a,b,M,fun,x1,x2,x,p,f,mitja,mitj2,des,var

      xnums=0.d0

      iseed=16404850 !seed niub
            call srand(iseed)

      do i=1,ndat
            f=-1
            p=0
            do while (f.lt.p)
                  x1=rand()
                  x2=rand()
                  x=(b-a)*x1+a
                  p=M*x2
                  f=fun(x)
            end do
            xnums(i)=x
      end do
 
      mitja=0.d0
      var=0.d0

      do i=1,ndat
        mitja=mitja+xnums(i)
      end do
      mitja=mitja/dble(ndat)
      
      do i=1,ndat
            mitj2=(xnums(i)-mitja)**2
            var=var+mitj2
      end do

      des=dsqrt(var)

      open(4,file='P6-18P-res.dat',status='unknown',position='append')
      write(4,600) 'Middle value','Variance','Standard deviation'
      write(4,500) mitja,var,des
      write(4,*) ' '
      write(4,*) ' '
      close(4)

500   format(2x,e20.12,5x,e20.12,5x,e20.12)
600   format(8x,a,14x,a,12x,a)

      return
      end subroutine

c Subroutine that applies the Montecarlo method in 1D

      subroutine montecarlop6()
      implicit none
      real*8 fun,xa,xb,pi,f1,f2,i1,sigma1
      real*8 i1exacte,x1,f3,f4,f5
      real*8 xnums(1:1000000),m,a,b
      real*8 s1,s2,div,factor1,l
      integer n,j,iseed,ndat
      external f1,f2,fun
      COMMON/DADES/pi,l
      COMMON/RANDOM/xnums    

c Calculation of the integrals of function 1 and 2

      i1exacte=(122.d0/27.d0)-(21.d0*(pi**2)/4.d0)! exact value of the first function integral 
c      write(*,*) i1exacte,i2exacte
c Limits first integral
      xa=0.d0 
      xb=(3.d0/2.d0)*pi 

      iseed=16404850 
      call srand(iseed)

100   format(10x,a,16x,a,21x,a)
200   format(2x,i10,2(5x,e20.12))


      open(1,file='P6-18P-res.dat',status='unknown')
      write(1,100) 'N','I1','Sigma1'
      write(1,*) ' '
     
c Initialize to zero
      
      s1=0.d0
      s2=0.d0
      n=1000000
c      write(*,*) (xb-xa),(xd-xc)

      do j=1,n 
            x1=xa+(rand()*(xb-xa))
            s1=s1+(f1(x1)*(xb-xa))
            s2=s2+((f1(x1)*(xb-xa))**2)
            i1=s1/dble(j)
            factor1=dsqrt((s2/dble(j))-((s1/dble(j))**2))
            sigma1=factor1/dsqrt(dble(j))
            if ( mod(j,10000).eq.0 ) then
                  write(1,200) j,i1,sigma1
            end if
      end do
      write(1,*) ' '
      write(1,*) ' '
      close(1)
c      write(*,*) i1exacte,i1,sigma1 

c Generates 1000000 numbers with probability density p(x)=fun(x)
c This numbers correspond to the position of an ultracold atom confined inside a 1D potential box, x€[-L/2,L/2]
c p(x)=(2/L)*sin^2((pi+(x-L/2))/L)

      ndat=1000000
      a=-l/2.d0
      b=l/2.d0
      M=0.15d0
      call subaccepta(ndat,xnums,a,b,M,fun)

c Computes with N=10000,20000,...,100000 the integral of g(x)p(x) between -L/2 and L/2
c g(x)=sin^2((8*pi(x-L/2))/L)

      open(2,file='P6-18P-res.dat',status='unknown',position='append')
      write(2,100) 'N','I2','Sigma2'
      write(2,*) ' '

      s1=0.d0
      s2=0.d0
      n=1000000

      do j=1,n
            s1=s1+(f2(xnums(j)))
            s2=s2+((f2(xnums(j)))**2)
            i1=s1/dble(j)
            factor1=dsqrt((s2/dble(j))-((s1/dble(j))**2))
            sigma1=factor1/dsqrt(dble(j))
            if ( mod(j,10000).eq.0 ) then
                 write(2,200) j,i1,sigma1
            end if  
      end do
      write(2,*) ' '
      write(2,*) ' '
      close(2)
c      write(*,*)i1,sigma1


      return
      end subroutine

c Subroutine that applies the multidimensional Montecarlo method
c Considers a system of three fermionic atoms in the potential that we described before.
c It's wave function can be written as phi(x1,x2,x3)=psi(x1,x2,x3)/sqrt(N), where
c part1= product from i=1 to 3 of sin((pi*(xi-L/2))/L)
c part 2= product of j<k =1,3 of (cos((pi*(xj-L/2))/L)-cos((pi*(xk-L/2))/L))
c psi(x1,x2,x3)=part1*part2
c Computes the standardization of the wave function

      subroutine mcarloMD()
      implicit none
      real*8 pi,l,psi,s1,s2,x1,x2,x3,xa,xb,fun,sigma1,factor1,i1,div
      real*8 xnums(1000000),a,b,m,ndat
      integer i,j,n
      external psi,fun
      COMMON/DADES/pi,l
      COMMON/RANDOM/xnums    

      s1=0.d0
      s2=0.d0
      n=1000000

      open(8,file='P6-18P-res.dat',status='unknown',position='append')
      write(8,700) 'N','I3','Sigma'
      write(8,*) ' '

      xa=-l/2.d0
      xb=l/2.d0
           
      do j=1,n-1,3
            x1=xnums(j)
            x2=xnums(j+1)
            x3=xnums(j+2)
            s1=s1+(psi(x1,x2,x3))
            s2=s2+((psi(x1,x2,x3))**2)
            i1=s1/dble(j)
            factor1=dsqrt((s2/dble(j))-((s1/dble(j))**2))
            sigma1=factor1/dsqrt(dble(j))
            if ( mod(j,10000).eq.0 ) then
                  write(8,800) j,i1,sigma1
            end if     
      end do
      write(8,*) ' '
      write(8,*) ' '
      close(8)
      

700   format(10x,a,16x,a,21x,a)
800   format(2x,i10,2(5x,e20.12))


      return 
      end subroutine


c                 FUNCTIONS      


c Probability density p(x)

      real*8 function fun(x)
      implicit none
      real*8 x,pi,div,temp,l
      COMMON/DADES/pi,l


      div=pi*(x-(l/2.d0))/l
      temp=(2.d0/l)*(dsin(div)**2)
      fun=temp

      return
      end function

c Function 1

      real*8 function f1(x)
      implicit none
      real*8 x,funcio

      funcio=(x**3)*(dsin(x)**3)
      f1=funcio
  
      return
      end function

c Function 2
 
      real*8 function f2(x)
      implicit none
      real*8 x,funcio,div,pi,l
      COMMON/DADES/pi,l

      div=8.d0*pi*(x-(l/2.d0))/l
      funcio=dsin(div)**2
      f2=funcio
  
      return
      end function

c Function 3

      real*8 function psi(x1,x2,x3)
      implicit none
      real*8 pi,l,x1,x2,x3,temp1,temp2,temp3,temp4,temp5,temp
      real*8 part1,part2,prob,fun
      COMMON/DADES/pi,l

      temp1=(pi*(x1-(l/2.d0)))/l
      temp2=(pi*(x2-(l/2.d0)))/l
      temp3=(pi*(x3-(l/2.d0)))/l

      part1=dsin(temp1)*dsin(temp2)*dsin(temp3)

      temp4=(dcos(temp1)-dcos(temp2))*(dcos(temp1)-dcos(temp3))
      temp5=dcos(temp2)-dcos(temp3)
      part2=temp4*temp5

      prob=fun(x1)*fun(x2)*fun(x3)

      temp=(part1*part2)**2
      psi=temp/prob

      return
      end function



