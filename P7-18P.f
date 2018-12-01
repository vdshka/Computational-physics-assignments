
C                      ASSIGNMENT 7

c     Study of the simple pendulum dynamics. Euler method.
c     Predictor-corrector method.


      program pre7
      implicit none
      real*8 m,l,g,tn,wn,pi,tmax,y0,f0,tmin,dt,f01,f02
      real*8 fun,epot,ekin,aprox
      real*8 kin,pot,total,kin2,pot2,total2
      real*8 y(2000),f(2000),y2(50001),f2(50001),y22(50001),f22(50001)
      integer j,k,npas,pascon(4)
      external fun,ekin,epot,aprox
      COMMON/DADES/m,l,g

      pi=dacos(-1.0d0) 
      tmin=0.d0
      pascon=(/200,600,4000,50000/)

c                      FORMATS

300   format(2x,7(e20.12,5x))
400   format(12x,a,12x,a,6x,a,3x,a,7x,a,5x,a,5x,a)
100   format(2x,5(e20.12,5x))
600   format(a,2x,i10)
700   format(10x,a,2x,e20.12)
200   format(10x,a,19x,a,3(15x,a))
500   format(29x,a,11x,a)
800   format(45x,a,29x,a)
250   format(10x,a,16x,a)


c Parameters of the pendulum described by the differential equation l*alpha''=-g*sin alpha
c alpha=angle, alpha''= angular acceleration

      m=0.95d0 ! mass in kg
      l=1.05d0 !length in m
      g=1.66d0 ! gravity in m/s^2
      wn=dsqrt(g/l) !angular frequency
      tn=(2.d0*pi)/wn !period
      tmax=6.d0*tn 

c Small-angle motion. 

      npas=1500 
      y0=0.15d0 !alpha(0) in rad
      f0=0.d0 ! alpha'(0) in rad/s
      dt=(tmax-tmin)/dble(npas) 

      open(1,file='P7-18P-res.dat',status='unknown')
      write(1,*) ' -------------SMALL-ANGLE------------'
      write(1,*) ' '
      write(1,600) ' Number of steps of time:',npas 
      write(1,*) ' '
      write(1,500) 'Values obtained with the Euler raw method: ','Values
     + obtained with predictor-corrector method:'
      write(1,*) ' '
      write(1,200) 't(s)','alpha(rad)','dalpha(rad/s)','alpha(rad)','dal
     +pha(rad/s)'
      write(1,*) ' '
      
      call euler(npas,fun,tmax,tmin,y0,f0,y,f)
      call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)

      do j=1,npas+1
      	write(1,100) dble(j-1)*dt,y(j),f(j),y2(j),f2(j)
      end do

      write(1,*) ' '
      write(1,*) ' '

c Small-angle approximation: sin(alpha)->alpha
      
      write(1,*) 'Small-angle approximation: sin(alpha)->alpha'
      write(1,*) ' '
      write(1,500) 'Values obtained with the Euler raw method: ','Values
     + obtained with predictor-corrector method:'
      write(1,*) ' '
      write(1,200) 't(s)','alpha(rad)','dalpha(rad/s)','alpha(rad)','dal
     +pha(rad/s)'
      write(1,*) ' '
      
      call euler(npas,aprox,tmax,tmin,y0,f0,y,f)
      call corrector(npas,aprox,tmax,tmin,y0,f0,y2,f2)

      do j=1,npas+1
            write(1,100) dble(j-1)*dt,y(j),f(j),y2(j),f2(j)
      end do

      write(1,*) ' '
      write(1,*) ' '

c Large-angle motion.

      y0=pi-0.15d0 !alpha(0) in rad
      f0=0.d0 ! alpha'(0) in rad/s
      npas=1500
      dt=(tmax-tmin)/dble(npas) 

      open(1,file='P7-18P-res.dat',status='unknown')
      write(1,*) ' -------------LARGE-ANGLE------------'
      write(1,*) ' '
      write(1,600) ' Number of steps of time:',npas 
      write(1,*) ' '
      write(1,500) 'Values obtained with the Euler raw method: ','Values
     + obtained with predictor-corrector method:'
      write(1,*) ' '
      write(1,200) 't(s)','alpha(rad)','dalpha(rad/s)','alpha(rad)','dal
     +pha(rad/s)'
      write(1,*) ' '
      
      call euler(npas,fun,tmax,tmin,y0,f0,y,f)
      call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)

      do j=1,npas+1
      	write(1,100) dble(j-1)*dt,y(j),f(j),y2(j),f2(j)
      end do

      write(1,*) ' '
      write(1,*) ' '

c Energy. Computes the Kinetic energy, potential energy and total energy of the pendulum.

      y0=pi-0.015d0
      f0=0.1d0
      npas=1500
      dt=(tmax-tmin)/dble(npas) 

      write(1,*) ' -------------ENERGY------------'
      write(1,*) ' '
      write(1,600) ' Number of steps of time:',npas 
      write(1,*) ' '
      

      call euler(npas,fun,tmax,tmin,y0,f0,y,f)
      call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)

      write(1,700) 'alpha:', y0
      write(1,*) ' '
      write(1,500) 'Values obtained with the Euler raw method: ','Values
     + obtained with predictor-corrector method:'
      write(1,*) ' '
      write(1,400)'t(s)','Kinetic energy (J)','Potential energy (J)',
     +'Total energy (J)','Kinetic energy (J)','Potential energy (J)'
     +,'Total energy (J)'
      write(1,*)' '
  
      do j=1,npas+1
      	kin=ekin(f(j))
      	kin2=ekin(f2(j))
      	pot=epot(y(j))
      	pot2=epot(y2(j))
      	total=kin+pot
      	total2=kin2+pot2
      	write(1,300) dt*(j-1),kin,pot,total,kin2,pot2,total2
      end do
      write(1,*) ' '
      write(1,*) ' '

c Transition. Pendulum dynamics from alpha(0)=0 with alpha'(0)=2*sqrt(g/l) +/- 0.1 rad/s

      y0=0.d0
      f01=2.d0*wn+0.1d0
      f02=2.d0*wn-0.1d0
      tmax=12.d0*tn
      npas=5000
      dt=(tmax-tmin)/dble(npas) 

      call corrector(npas,fun,tmax,tmin,y0,f01,y2,f2)
      call corrector(npas,fun,tmax,tmin,y0,f02,y22,f22)

      write(1,*) ' -------------TRANSITION------------'
      write(1,*) ' '
      write(1,600) ' Number of steps of time:',npas  
      write(1,*) ' '
      write(1,*) 'Values obtained with predictor-corrector method: '
      write(1,*) ' '
      write(1,200) 't(s)','alpha+(rad)','dalpha+(rad/s)','alpha-(rad)','
     +dalpha-(rad/s)'
      write(1,*) ' '

      do j=1,npas+1
      	write(1,100) dt*(j-1),y2(j),f2(j),y22(j),f22(j)
      end do
      
      write(1,*) ' '
      write(1,*) ' '

c Convergence of the predictor-corrector method.

      write(1,*) ' -------------CONVERGENCE------------'
      write(1,*) ' '
      write(1,*) 'Values obtained with predictor-corrector method: '
      write(1,*) ' '

      do k=1,4
        y0=3.d0       
        f0=0.d0        
        tmax=10.d0*tn
      	npas=pascon(k)
        dt=(tmax-tmin)/dble(npas)
      	call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)
      	write(1,600) ' Number of steps of time:',npas 
        write(1,*) ' '
        write(1,250) 't(s)','Total energy(J)'
        write(1,*) ' '
      	do j=1,npas+1
      		kin=ekin(f2(j))
      		pot=epot(y2(j))
      		total=kin+pot
      		write(1,100) dt*(j-1),total
      	end do
        write(1,*) ' '
        write(1,*) ' '
      end do

      close(1)

c Animation. Generates data to create a gif animation of the pedulum motion. 
      
      open(2,file='anime.dat',status='unknown')

      y0=pi-0.15d0 !alpha(0) in rad
      f0=0.d0 ! alpha'(0) in rad/s
      npas=1500
      dt=(tmax-tmin)/dble(npas) 
     
      call euler(npas,fun,tmax,tmin,y0,f0,y,f)
      call corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)

      do j=1,npas+1
      	write(2,*) dble(j-1)*dt,y(j),f(j),y2(j),f2(j)
      end do

      write(2,*) ' '
      write(2,*) ' '

      close(2)

      return
      end program


c              SUBROUTINES

c Euler's method to compute y and the first derivative of y as a function of time
c step of time=dt, number of steps=npas
c yo=fuction at 0, f0= first derivative at 0

      subroutine euler(npas,fun,tmax,tmin,y0,f0,y,f)
      implicit none
      real*8 y(2000),f(2000),tmax,dt,y0,f0,fun,tmin
      integer npas,i

      y=0.d0  
      f=0.d0

      dt=(tmax-tmin)/dble(npas)

      y(1)=y0
      f(1)=f0

      do i=2,npas+1

      	y(i)=y(i-1)+dt*f(i-1)    
      	f(i)=f(i-1)+dt*fun(y(i-1))

      end do
      
      return
      end subroutine

c Predictor-corrector method

      subroutine corrector(npas,fun,tmax,tmin,y0,f0,y2,f2)
      implicit none
      real*8 y2(50001),f2(50001),tmax,dt,y0,f0,fun,tmin,y1,fa,ya,f1
      integer npas,i

      y2=0.d0 
      f2=0.d0

      dt=(tmax-tmin)/dble(npas)

      ya=y0
      fa=f0
c      y1=y0+dt*f0
c      f1=f0+dt*fun(y0)

      do i=1,npas+1

      	y1=ya+dt*fa
        f1=fa+dt*fun(ya)

      	y2(i)=ya+(dt/2.d0)*(fa+f1)
      	f2(i)=fa+(dt/2.d0)*(fun(ya)+fun(y1))

      	ya=y2(i)
      	fa=f2(i)

      end do
      
      return
      end subroutine


c                 FUNCTIONS      

c Pendulum function
c Input: x=angle

      real*8 function fun(x)
      implicit none
      real*8 x,g,l,temp,m
      COMMON/DADES/m,l,g

      temp=-(g/l)*dsin(x)
      fun=temp

      return
      end function

c Pendulum function for the small-angle approximation
c Input: x=angle

      real*8 function aprox(x)
      implicit none
      real*8 x,g,l,temp,m
      COMMON/DADES/m,l,g

      temp=-(g/l)*x
      aprox=temp

      return
      end function

c Kinetic energy function
c Input: x=angular velocity

      real*8 function ekin(x)
      implicit none
      real*8 x,m,temp,l,g
      COMMON/DADES/m,l,g

      temp=0.5d0*m*(l**2)*(x**2)
      ekin=temp

      return
      end function

c Potential energy function
c Input: x=angle

      real*8 function epot(x)
      implicit none
      real*8 x,m,temp,g,l
      COMMON/DADES/m,l,g

      temp=-(m*g*l*dcos(x))
      epot=temp

      return
      end function

