c -----------------------------------------------------------------------
C                      PRÀCTICA 10
c-----------------------------------------------------------------------
c evolució temporal de la funció T(x,t)
C ----------------------------------------------------------------------
      program pre10
      implicit none
      integer xpunts,i,j,icontrol,npas,z,x1,x2,x3,x4
      real*8 sigma,lx,h,t0,tlx,ro,x,k,dt,integral
      real*8 tini(0:150),tf(0:150),kappa(3),xp(4)
      real*8 psi1(12000,149),t,r1(0:150),temp(12000)
      common/dades/h,dt
      external ro

c-----------------FORMATS-----------------------------------------------
c ----------------------------------------------------------------------
100   format(2x,a,f7.2)
200   format(10x,a,20x,4(a,14x))
300   format(10x,a,16x,a)
400   format(2x,2(e20.12,5x))
500   format(2x,5(e20.12,5x))
600   format(2x,a,f7.3)
c ----------------------------------------------------------------------

c dades
      sigma=1.d-7
      lx=30.d0
      h=0.2d0
      xpunts=int(lx/h)
      kappa=(/3.d0,7.4d0,10.d0/)
      xp=(/6.d0,10.d0,16.d0,28.d0/)
      npas=12000 ! numero de passos de temps
      dt=4.d-3

c condicions de contorn
      t0=10.d0 !T(0,t)
      tlx=45.d0 !T(Lx,t)

c inicialitzem vectors a zero
      tini=0.d0
      r1=0.d0
      temp=0.d0

c APARTAT 1-ESTUDI T0, quan t=0
      
c posem condicions de contorn i omplim resta matriu amb T=10ºC
      tini(0)=t0
      tini(xpunts)=tlx 
      do i=1,xpunts-1
      	tini(i)=10.d0
      end do

c escrivim a fitxer dades de T0      
      open(1,file='P10-18P.dat',status='unknown')
      do j=1,3
      	icontrol=j 
      	call Gauss(j,xpunts,ro,sigma,tini,tf)
      	if(icontrol.eq.1) write(1,*)'Metode Gauss-Seidel'
      	if(icontrol.eq.2) write(1,*)'Metode Jacobi'
      	if(icontrol.eq.3) write(1,*)'Metode Sobrerelaxacio'
      	write(1,*) ' '
      	write(1,300) 'x(cm)','Temperatura(ºC)'
      	write(1,*) ' '
      	do i=0,xpunts
      		x=0.d0+i*h
      		write(1,400) x,tf(i)
      	end do
      	write(1,*) ' '
      	write(1,*) ' '
      end do

c APARTAT 2-ESTUDI T(x,t) implementem mètode implícit   
c----- PART A) kappa=7.4

      write(1,*) 'Evolucio temperatures punts xp = 6,10,16,28 amb k=7.4'
      write(1,*) ' '
      write(1,200) 't(s)','T(6,t)(ºC)','T(10,t)(ºC)','T(16,t)(ºC)',
     +'T(28,t)(ºC)'
      write(1,*) ' '

      k=kappa(2)
      call implicit(k,tf,xpunts,npas,psi1)
      x1=int(xp(1)/h)
      x2=int(xp(2)/h)
      x3=int(xp(3)/h)
      x4=int(xp(4)/h)

      do i=1,npas
      	t=0.d0+i*dt
        write(1,500) t,psi1(i,x1),psi1(i,x2),psi1(i,x3),psi1(i,x4)
      end do

      write(1,*) ' '
      write(1,*) ' '

c----PART B) kappa= 3,7.4,10. Evolució temporal mitjana de T 

      write(1,*) 'Evolucio temporal de la temperatura mitjana'
      write(1,*) ' '

c condicions de contorn en vector r1 que emprarem per a calcular la temperatura mitjana
      r1(0)=tf(0)
      r1(xpunts)=tf(xpunts)

c calculem amb mètode implícit per a cada valor de kappa      
      do j=1,3
      	k=kappa(j)
      	write(1,100) 'kappa =',k
      	write(1,*) ' '
      	write(1,300) 't(s)','Tmitjana(t)(ºC)'
      	write(1,*) ' '

        call implicit(k,tf,xpunts,npas,psi1)

      	do i=1,npas
      		t=0.d0+i*dt
          do z=1,xpunts-1
            r1(z)=psi1(i,z) 
          end do
        	call simpson(0.d0,lx,xpunts,r1,integral)
        	temp(i)=(integral)/lx
        	write(1,400) t,temp(i)
      	end do
      	write(1,*) ' '
        write(1,*) ' '
      end do

c----PART B) kappa=10 evolució del perfil de temperatures amb el temps per fer animació
      
      k=kappa(3)
      call implicit(k,tf,xpunts,npas,psi1)

      do i=1,npas
      	t=0.d0+i*dt
      	if(mod(i,100).eq.0) then
      		write(1,600) 'Temps=',t
      		write(1,*) ' '
      		do j=1,xpunts-1
      			x=0.d0+j*h
      			write(1,400) x,psi1(i,j)
      		end do
      		write(1,*) ' '
      		write(1,*) ' '
      	end if
      end do

      close(1)
   
      stop 
      end program

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c subrutina que implementa el mètode Gauss-Siedel o sobrerelaxació
c Elecció mètode controlada per variable icontrol. 
c icontrol=1 -> gauss-siedel, icontrol=2 -> jacobi, icontrol=3 -> sobrerelaxació
c ti és la matriu inicial,tf es la matriu final i sigma és la tolerància
c xpunts és el nombre de punts en que hem dividit x

      subroutine Gauss(icontrol,xpunts,ro,sigma,ti,tf)
      implicit none
      real*8 ro,sigma,h,w,epsilon,error,x
      real*8 ti(0:xpunts),tf(0:xpunts),b(0:xpunts)
      integer xpunts,icontrol,j,k
      common/dades/h

      w=1.5d0

c inicialitzem vectors a zero
      tf=0.d0
      b=0.d0

c escrivim al vector les condcions de contorn
      tf(0)=ti(0)
      tf(xpunts)=ti(xpunts)

      epsilon=1.d0 ! asignem valor gran a epsilon per entrar al bucle 
      k=0 !nombre iteració

c calculem punts centrals      
      do while(epsilon.gt.sigma)
      	epsilon=0.d0 ! ara que ja estem a dins del bucle inicialitzem epsilon a zero 
      	do j=1,xpunts-1 !recorre les x sense passar per valors determinats per cond contorn
      		x=0.d0+j*h
        	if(icontrol.eq.1) then ! Gauss-Siedel
        		tf(j)=(ti(j+1)+ti(j-1)+(h**2)*ro(x))/2.d0
            else if(icontrol.eq.2) then !Jacobi
            	tf(j)=(ti(j+1)+tf(j-1)+(h**2)*ro(x))/2.d0
        	else if(icontrol.eq.3) then !Sobrerelaxació
        		b(j)=((ti(j+1)+tf(j-1)+(h**2)*ro(x))/2.d0)-ti(j)    
        		tf(j)=ti(j)+(w*b(j))
        	end if                                         
        	error=dabs(ti(j)-tf(j))
      		if(error.gt.epsilon) epsilon=error
      	end do
        do j=1,xpunts-1
        	ti(j)=tf(j)
        end do
      	k=k+1
      end do

      return
      end subroutine

c--------- [Nota: He modificat una mica la subrutina] ----------
c Resol el problema T*psi=r 
c T = matriu tridiagonal, A = diagonal inferior, afegim cero primera component
c B = diagonal central, C = diagonal superior, afegim cero a la última component
c els zeros afegits a A i C són per a que tinguin la mateixa dimensió que B
c inputs: vectors A,B,C,R , imax = dimensió vectors
c output: vector PSI

      subroutine tridiag(a,b,c,r,psi,imax)
      implicit none
      real*8 bet,gam(4001)
      real*8 a(imax),b(imax),c(imax),r(imax),psi(imax)
      integer imax,j

      psi=0.d0 !inicialitzem a zero

      if(b(1).eq.0.d0) write(*,*)'Atencio! sistema indeterminat'
      bet=b(1)
      psi(1)=r(1)/bet
      do j=2,imax
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.d0) write(*,*) 'Atencio! bet=0'
        psi(j)=(r(j)-a(j)*psi(j-1))/bet
      end do

      do j=imax-1,1,-1
        psi(j)=psi(j)-gam(j+1)*psi(j+1)
      end do

      return
      end subroutine

c subrutina metode implicit
c inputs kappa, vector per a t=0 (tf), numero de punts per x (xpunts), numero de punts per t (npas)
c output matriu psi1(xpunts,npas)

      subroutine implicit(k,tf,xpunts,npas,psi1)
      implicit none
      integer xpunts,npas,i,j
      real*8 k,tf(0:xpunts),psi1(npas,xpunts-1),dt,h
      real*8 alpha,a(xpunts-1),b(xpunts-1),c(xpunts-1),t
      real*8 r(xpunts-1),psi(xpunts-1)
      common/dades/h,dt

c diagonals
      alpha=k*dt/(h**2)      
      a=-alpha
      b=1.d0+2.d0*alpha
      c=-alpha
      
c posem zeros per a igualar dimensions
      a(1)=0.d0
      c(xpunts-1)=0.d0

      do i=1,xpunts-1
        r(i)=tf(i)
      end do

      r(1)=r(1)+alpha*tf(0)
      r(xpunts-1)=r(xpunts-1)+alpha*tf(xpunts)

c càlcul

      do i=1,npas
        t=0.d0+i*dt
        call tridiag(a,b,c,r,psi,xpunts-1)
        do j=1,xpunts-1
          psi1(i,j)=psi(j)
        end do
        r=psi
        r(1)=r(1)+alpha*tf(0)
        r(xpunts-1)=r(xpunts-1)+alpha*tf(xpunts)
      end do

      return
      end subroutine
      
c subrutina que implementa el metode de simpson compost
c retorna integral, que es la aproximacio obtinguda per a la integral emprant aquest metode

      subroutine simpson(a,b,npas,fcn,integral)
      implicit none
      real*8 h,x,a,b,fcn(0:npas),i,integral
      integer j,npas

      i=fcn(0)+fcn(npas)
      h=(b-a)/dble(npas-1)
      do j=0,npas
        x=a+dble(j)*h
        if (mod(j,2).eq.0) i=i+2.d0*fcn(j)
        if (mod(j,2).ne.0) i=i+4.d0*fcn(j)
      end do  
      integral=i*h/3.d0

      return
      end subroutine


c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         FUNCIONS
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c funció densitat per a l'estat inicial
      real*8 function ro(x)
      implicit none
      real*8 x,temp,e

      e=dabs(x-15.d0)**2/(0.8d0**2)
      temp= 2.d0*dexp(-e)
      ro=temp

      end function

c ----------------------------------------------------------------------
