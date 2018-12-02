c -----------------------------------------------------------------------
C                      PRÀCTICA 9
c-----------------------------------------------------------------------
c resoldre l'eq de Poisson en 2D en una geometria rectangular 
C amb condicions de contorn de Drichlet
C ----------------------------------------------------------------------
      program prac9
      implicit none
      integer xpunts,ypunts,i,j,n
      real*8 lx,ly,tyl,txl,h,sigma,x0,y0,x,y,t0x,t0y
      real*8 tini(0:67,0:91),tinterior(3),tf(0:67,0:91),tf0(0:67,0:91)
      common/dades/h

c-----------------FORMATS-----------------------------------------------
c ----------------------------------------------------------------------
200   format(2(2x,a,f7.2))
500   format(2x,a,2x,f7.2,a)
300   format(12x,a,23x,a,21x,a)
400   format(2x,3(e20.12,5x))
c ----------------------------------------------------------------------

      sigma=1.d-7
      tinterior=(/10.d0,120.d0,1040.d0/)
      tini=0.d0 !inicialitzem matriu a zero
      
c dimensions rectangle en cm 
      lx=33.5d0 
      ly=45.5d0

c condicions de contorn
      t0y=17.d0 !temperatura en t(0,y) 
      t0x=0.5d0 !t(x,0) mesurada en ºC
      txl=11.2d0 !temperatura t(x,ly)
      tyl=25.3d0 !temperatura en t(lx,y)

      h=0.5d0 ! en cm, gruix malla
      xpunts=int(lx/h)
      ypunts=int(ly/h)

c escrivim a la matriu les condcions de contorn
      do i=0,xpunts
      	tini(i,0)=t0x
      	tini(i,ypunts)=txl
      	do j=1,ypunts
      		tini(0,j)=t0y
      		tini(xpunts,j)=tyl
      	end do
      end do 

c estudi convergencia temperatura en el punt (x,y)=(15.5,23.5) amb els 2 mètodes
      x0=15.5d0
      y0=23.5d0

      open(1,file='P9-18P.dat',status='unknown')
      write(1,*) ' '
      write(1,200) 'Estudi convergencia temperatura en x =',x0,'y =',y0
      write(1,*) ' '
      close(1)

      do n=1,3
      	open(3,file='P9-18P.dat',status='unknown',position='append')
      	write(3,500) 'Temperatura interior:',tinterior(n),'ºC'
      	write(3,*) ' '
      	close(3)
c omplim resta matriu amb altres valors, tinterior=15ºC,120ºC,540ºC
      	do i=1,xpunts-1
      		do j=1,ypunts-1
      			tini(i,j)=tinterior(n)
      		end do
      	end do
        call metode(1,ypunts,xpunts,sigma,x0,y0,tini,tf) !mètode gauss-seidel
        call metode(2,ypunts,xpunts,sigma,x0,y0,tini,tf) !metode sobrerelaxació
        call metode1(2,ypunts,xpunts,sigma,tini,tf0) !metode sobrerelaxació
      end do

c obtenim dades per a fer gràfica 2D amb fogons i sense

      open(4,file='Amb_fogons.dat',status='unknown')
      write(4,*) 'Dades grafica 3D'
      write(4,*) ' '
      write(4,300) 'x','y','Temperatura'
      write(4,*) ' '
      open(5,file='Sense_fogons.dat',status='unknown')
      write(5,*) 'Dades grafica 3D'
      write(5,*) ' '
      write(5,300) 'x','y','Temperatura'
      write(5,*) ' '
      
      do i=0,ypunts 
          y=0.d0+i*h
          do j=0,xpunts 
            x=0.d0+j*h
            write(4,400) x,y,tf(j,i)
            write(5,400) x,y,tf0(j,i)
          end do
          write(4,*) ' '
          write(5,*) ' '
      end do

      close(4)
      close(5)

      stop
      end program

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c subrutina que calcula la densitat de temperatura (ºC/cm^2) en funció del punt (x,y)
c tenint en compte la contribució de les dues fonts de calor 
c ro(x,y)=ro1(x,y)+ro2(x,y)+ro3(x,y)

      subroutine densitat(x,y,ro)
      implicit none
      real*8 r,ro1,x,y,ro2,ro,a,b,c,d,e,r1,ro3,e1

c primer fogó escalfa en una circumferència centrada al punt (8,22.5)
      r=dsqrt(((x-8.d0)**2)+((y-22.5d0)**2))
      e=(r-5.d0)**2/(0.3d0**2)
      ro1=10.d0*dexp(-e)

c segon fogó escalfa en un rectangle de 4cmx6cm centrat a (x,y)=(20,32)
c extrems rectangle x en cm
      a=18.d0
      b=22.d0
c extrems rectangle y en cm
      c=29.d0
      d=35.d0

      if((x.ge.a).and.(x.le.b).and.(y.ge.c).and.(y.le.d)) then
      	ro2=3.d0
      else
      	ro2=0.d0
      end if

c tercer fogó 
      r1=dsqrt(((x-22.d0)**2)+((y-10.5d0)**2))
      e1=(r1-4.d0)**2/(0.8d0**2)
      ro3=6.d0*dexp(-e1)

c calculem densitat total      
      ro=ro1+ro2+ro3

      return 
      end subroutine

c subrutina un pas de gauss-seidel o sobrerelaxacio. Elecció mètode controlada
c per variable icontrol. icontrol=1 -> gauss, icontrol=2 ->sobrerelaxació
c ti és la matriu inicial,tf es la matriu final i sigma és la tolerància
c xpunts és el nombre de punts en que hem dividit x
c ypunts és el nombre de punts en que hem dividit y
c (x0,y0) punt que ens interessa estudiar
c posem com a argument fogo = 1 si hi ha fogons i fogo = 0 si no hi ha fogons

      subroutine metode(icontrol,ypunts,xpunts,sigma,x0,y0,tini,tf)
      implicit none
      real*8 a(0:xpunts,0:ypunts),b(0:xpunts,0:ypunts) ! a i b son matrius temporals per a operacions intermitjes
      real*8 tini(0:xpunts,0:ypunts),tf(0:xpunts,0:ypunts)
      real*8 ti(0:xpunts,0:ypunts)
      real*8 sigma,epsilon,w,error,h,ro,x,y,x0,y0
      integer icontrol,j,k,xpunts,ypunts,i
      common/dades/h

c     ------------FORMATS-----------------------------------------------
c     ------------------------------------------------------------------
100   format(6x,i10,18x,e20.12)
300   format(12x,a,23x,a)
c     ------------------------------------------------------------------
      
      w=1.35d0

c inicialitzem matrius a zero      
      tf=0.d0
      a=0.d0
      b=0.d0

c passem dades de tini a ti per a treballar amb aquesta última i no sobreescriure la primera
      do i=0,xpunts
      	do j=0,ypunts
      		ti(i,j)=tini(i,j)
      	end do
      end do
      
c escrivim a la matriu les condcions de contorn
      do i=0,xpunts
      	tf(i,0)=ti(i,0)
      	tf(i,ypunts)=ti(i,ypunts)
      	do j=1,ypunts
      		tf(0,j)=ti(0,j)
      		tf(xpunts,j)=ti(xpunts,j)
      	end do
      end do

      epsilon=1.d0 ! asignem valor gran a epsilon per entrar al bucle 
      k=0 !nombre iteració

      open(2,file='P9-18P.dat',status='unknown',position='append')
      if(icontrol.eq.1) write(2,*)'Mètode Gauss-Seidel'
      if(icontrol.eq.2) write(2,*)'Mètode Sobrerelaxació'
      write(2,*) ' '
      write(2,300) 'Iteració','Valor'
      write(2,*) ' '

c calculem punts centrals      
      do while(epsilon.gt.sigma)
        epsilon=0.d0 ! ara que ja estem a dins del bucle inicialitzem epsilon a zero 
        do i=1,ypunts-1 !recorre les y sense passar per valors determinats per cond contorn
          y=0.d0+i*h
          do j=1,xpunts-1 !recorre les x sense passar per valors determinats per cond contorn
            x=0.d0+j*h
            call densitat(x,y,ro) !calculem la densitat de temperatura al punt desitjat
            a(j,i)=(ti(j+1,i)+ti(j,i+1)+(h**2)*ro)
            if(icontrol.eq.1) then
              tf(j,i)=(a(j,i)+ti(j-1,i)+ti(j,i-1))/4.d0
            else if(icontrol.eq.2) then
              b(j,i)=((a(j,i)+tf(j-1,i)+tf(j,i-1))/4.d0)-ti(j,i)    
              tf(j,i)=ti(j,i)+(w*b(j,i))
            end if
            error=dabs(ti(j,i)-tf(j,i))
            if(error.gt.epsilon) epsilon=error
            if((x.eq.x0).and.(y.eq.y0)) write(2,100) k,tf(j,i)
          end do
        end do
        do i=1,ypunts-1
          do j=1,xpunts-1
            ti(j,i)=tf(j,i)
          end do
        end do
        k=k+1
      end do
      write(2,*) ' '
      write(2,*) ' '

      close(2)

      return 
      end subroutine

c sense fogons

      subroutine metode1(icontrol,ypunts,xpunts,sigma,tini,tf)
      implicit none
      real*8 a(0:xpunts,0:ypunts),b(0:xpunts,0:ypunts) ! a i b son matrius temporals per a operacions intermitjes
      real*8 tini(0:xpunts,0:ypunts),tf(0:xpunts,0:ypunts)
      real*8 ti(0:xpunts,0:ypunts)
      real*8 sigma,epsilon,w,error
      integer icontrol,j,k,xpunts,ypunts,i
      
      w=1.35d0

c inicialitzem matrius a zero      
      tf=0.d0
      a=0.d0
      b=0.d0

c passem dades de tini a ti per a treballar amb aquesta última i no sobreescriure la primera
      do i=0,xpunts
        do j=0,ypunts
          ti(i,j)=tini(i,j)
        end do
      end do
      
c escrivim a la matriu les condcions de contorn
      do i=0,xpunts
        tf(i,0)=ti(i,0)
        tf(i,ypunts)=ti(i,ypunts)
        do j=1,ypunts
          tf(0,j)=ti(0,j)
          tf(xpunts,j)=ti(xpunts,j)
        end do
      end do

      epsilon=1.d0 ! asignem valor gran a epsilon per entrar al bucle 
      k=0 !nombre iteració

c calculem punts centrals      
      do while(epsilon.gt.sigma)
        epsilon=0.d0 ! ara que ja estem a dins del bucle inicialitzem epsilon a zero 
        do i=1,ypunts-1 !recorre les y sense passar per valors determinats per cond contorn
          do j=1,xpunts-1 !recorre les x sense passar per valors determinats per cond contorn
            a(j,i)=(ti(j+1,i)+ti(j,i+1))
            if(icontrol.eq.1) then
              tf(j,i)=(a(j,i)+ti(j-1,i)+ti(j,i-1))/4.d0
            else if(icontrol.eq.2) then
              b(j,i)=((a(j,i)+tf(j-1,i)+tf(j,i-1))/4.d0)-ti(j,i)    
              tf(j,i)=ti(j,i)+(w*b(j,i))
            end if
            error=dabs(ti(j,i)-tf(j,i))
            if(error.gt.epsilon) epsilon=error
          end do
        end do
        do i=1,ypunts-1
          do j=1,xpunts-1
            ti(j,i)=tf(j,i)
          end do
        end do
        k=k+1
      end do

      return 
      end subroutine


c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
