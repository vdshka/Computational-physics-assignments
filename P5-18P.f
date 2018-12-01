
C                      ASSIGNMENT 5

c     Browninan motion.Generation of random numbers.Box-Müller.


      program pre5
      implicit none
      real*8 xhisto(1:120),errhisto(1:120),histo(1:120),tf
      real*8 xgaus(1:120000),mu,sigma,xa,xb,delta,x(1:250),y(1:250),t,z
      real*8 var,inx,iny,suma,suma2,mitj2,variancia(1:240),mitj,fun,mes
      integer ndat,ncaixes,i,j,nm,contador
      external fun

c                           FORMATS

100   format(12x,11(a,23x))
200   format(2x,11(e20.12,5x))
300   format(12x,a,20x,a)
400   format(2x,e20.12,5x,e20.12)


c Generates 120000 gaussian numbers with mean equal to zero and variance equal to 1
c Generates a normalized histogram (probability density) with 120 boxes between -5 and 5

      ndat=120000
      mu=0.d0
      sigma=1.d0
      call subgaussians(ndat,mu,sigma,xgaus)

      ncaixes=120
      xa=-5.d0
      xb=5.d0
      call hist(ndat,xgaus,xa,xb,ncaixes,xhisto,histo,errhisto)

      open(6,file='P5-18P-res1.dat',status='unknown',position='append')
      mes=(xb-xa)/150.d0
      z=0.d0
      do while(z.lt.xb)
      	z=z+mes
      	write(6,400) z, fun(z)
      end do
      close(6)


c Simulates the random movement of Nm=250 independent oxigen molecules in two dimensions.
c The molecules are located in the origin (x,y)=(0,0) at the time t=0 and evolve with time so that for each step of time
c xn(t+dt)=xn(t)+dx, n=1,...,Nm ; yn(t+dt)=yn(t)+dy, n=1,...,Nm
c where dx and dy follow a normal distribution with mean equal to zero and variance equal to delta*dt
c use the gaussian numbers generated before to obtain dx and dy

      delta=0.02d0 !dt, increase of time, [dt]=s
      nm=250 
      var=delta*2.21d-5 ! variance, [var]=m^2
      x=0.d0 
      y=0.d0 

      open(1,file='P5-18P-res2.dat',status='unknown')
      write(1,100) 't','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5'

      contador = 1
      do i=1,240 !240 steps of time
      	t=i*delta
        suma=0.d0
        suma2=0.d0
      	do j=1,nm !nm components of x and y
      		inx=dsqrt(var)*xgaus(contador) !increase of x, dx
      		iny=dsqrt(var)*xgaus(contador+1) !increase of y, dy
      		x(j)=x(j)+inx
      		y(j)=y(j)+iny
      		suma=suma+y(j)
      		suma2=suma2+y(j)**2
      		contador=contador+2
      	end do
      	write(1,200) t,x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5)
      	mitj=suma/dble(nm)
      	mitj2=suma2/dble(nm)
      	variancia(i)=mitj2-mitj**2
      end do
      write(1,*) ' '
      write(1,*) ' '
      close(1)

c Create a file with the values of the variance for each step of time

      open(2,file='P5-18P-res2.dat',status='unknown',position='append')
      write(2,300) 't','Var(y(t))'
      do i=1,240
      	t=i*delta
      	write(2,400) t,variancia(i)
      end do
      write(2,*) ' '
      write(2,*) ' '
      close(2)

      stop
      end program


c              SUBROUTINES

c subroutine that creates a normalized histogram 
c inputs: ndat=number of data, xdata=vector of the values of x with ndata components,
c [xa,xb]=interval limits, ncaixes=number of boxes
c outputs: xhisto=vector of the values of the central position of the boxes with ncaixes components
c histo=vector with ncaixes components that correspond to the number of values of x in each box divided by the total number of values
c errhisto=vector of the error in each box with ncaixes components

      subroutine hist(ndat,xdata,xa,xb,ncaixes,xhisto,histo,errhisto)
      implicit none
      integer ncaixes,ndat,i,j
      integer inthisto(1:150)
      real*8 xdata(1:120000),xhisto(1:120),errhisto(1:120),histo(1:120)
      real*8 xa,xb,l,valor,linf,lsup,div,sigma,suma
       
      inthisto=0  ! initialize all the components to zero
      l=dabs(xb-xa)/(ncaixes-1)

      do j=1,ndat
      	valor=xdata(j)
        if ((valor.le.xb).and.(valor.ge.xa)) then !the value has to be inside the interval
        	do i=1,ncaixes
        		xhisto(i)=xa+(i-1)*l
                linf=xhisto(i)-(l/2.d0)
                lsup=xhisto(i)+(l/2.d0)
                if((valor.gt.linf).and.(valor.le.lsup)) then !if the value is inside the box counts +1
                	inthisto(i)=inthisto(i)+1
                end if
             end do
        end if
      end do

c Normalizes nk=nk/wkN 
      
      suma=0.d0
      do i=1,ncaixes
      	div=dble(inthisto(i))/(dble(ndat)*l) 
        histo(i)=div                         
        sigma=dsqrt(div*(1.d0-div)/dble(ndat))/l  
c        errhisto(i)=3.d0*sigma 
        errhisto(i)=sigma 
c        suma=suma+histo(i)
      end do
c      write(*,*)suma

c Writes the values in a file

500   format(2x,e20.12,5x,e20.12,5x,e20.12)
600   format(8x,a,16x,a,16x,a)

      open(3,file='P5-18P-res1.dat',status='unknown')
      write(3,*) 'Data of the normalized histogram:'
      write(3,*) ' '
      write(3,600) 'Midpoint', 'Number of values','Error'
c      write(3,*) ' '
      do i=1,ncaixes
      	write(3,500) xhisto(i),histo(i),errhisto(i)
      end do
      write(3,*) ' '
      write(3,*) ' '
      close(3)   

      return
      end subroutine

c Subroutine that generates ndat gaussian numbers with mean equal to xmu and variance equal to xsigma^2

      subroutine subgaussians(ndat,xmu,xsigma,xgaus)
      implicit none
      integer ndat,iseed,j
      real*8 xmu,xsigma,xgaus(1:120000),pi,x1,x2,z1,z2
      
      pi=dacos(-1.0d0)

      iseed=16404850 
      call srand(iseed)

      do j=1,ndat-1,2
      	x1=rand()
      	x2=rand()
      	z1=dsqrt(-2.d0*dlog(x1))*dcos(2.d0*pi*x2)
        z2=dsqrt(-2.d0*dlog(x1))*dsin(2.d0*pi*x2)
        xgaus(j)=z1*xsigma+xmu
        xgaus(j+1)=z2*xsigma+xmu
      end do

      return
      end subroutine


c                         FUNCTIONS


      real*8 function fun(x)
      implicit none
      real*8 x,pi,div,temp

      pi=dacos(-1.0d0)
      div=1/dsqrt(2.d0*pi)
      temp=div*dexp(-(x**2)/2.d0)
      fun=temp

      return
      end function
