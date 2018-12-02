c examen 2018 programa
      program main
      implicit none
      real*8 time,t,timef,valor
      integer l,s(50,50),m,i,iseed,j,suma,n,pas,mini,s1(50,50)

c-----------------FORMATS-----------------------------------------------
c ----------------------------------------------------------------------
200   format(2x,e20.12,5x,i6)
100   format(2x,a,e20.12)
300   format(2x,a,6x,a)
400   format(2x,2(i6,10x))
c ----------------------------------------------------------------------

      l=50
      t=2.d0
      n=l**2

      iseed=16404850 !seed niub
      call srand(iseed)

c part a)
c omple matriu inicial s(l,l) amb valors aleatoris de espin +1 o -1 
      suma=0

      do i=1,l
      	do j=1,l
      		valor=rand()
      		if(valor.le.0.5d0) then
      			s(j,i)=1
      		else if (valor.gt.0.5d0) then
      			s(j,i)=-1
      		end if
      		suma=suma+s(j,i)
      	end do
      end do
      
      s1=s
      mini=suma
c      write(*,*) mini

c temps = 0      
      timef=0.d0 
      time=0.d0
      m=mini
      pas=0

      open(1,file='valorst11.dat',status='unknown')
c      write(1,100) 'Calcul per a temperatura T=',t
c      write(1,100) 'Fins a temps t=',timef
c      write(1,*) ' '        
      write(1,300) 'coordenada i', 'coordenada j' ! se que esta al reves dels indexs, però aqui es referix a s(i,j)
      do while(time.lt.timef)
            call metropolis(s,l,t,m,time)
      end do
      pas=pas+1
      do i=1,l
            do j=1,l
      	     if(s(j,i).eq.1) then
      	     write(1,400) j,i
                  end if
      	end do
      end do
      write(1,*)' '
      write(1,*)' '
      close(1)

c temps = 1      
      timef=1.d0
      time=0.d0 
      m=mini
      pas=0
      s=s1

      open(3,file='valorst12.dat',status='unknown')
c      write(3,100) 'Calcul per a temperatura T=',t
c      write(3,100) 'Fins a temps t=',timef
c      write(3,*) ' '        
      write(3,300) 'coordenada i', 'coordenada j' ! se que esta al reves dels indexs, però aqui es referix a s(i,j)
      do while(time.lt.timef)
            call metropolis(s,l,t,m,time)
      end do
      pas=pas+1
      do i=1,l
            do j=1,l
                 if(s(j,i).eq.1) then
                 write(3,400) j,i
                  end if
            end do
      end do
      write(3,*)' '
      write(3,*)' '
      close(3)

c temps = 10      
      timef=10.d0
      time=0.d0 
      m=mini
      pas=0
      s=s1

      open(4,file='valorst13.dat',status='unknown')
c      write(4,100) 'Calcul per a temperatura T=',t
c      write(4,100) 'Fins a temps t=',timef
c      write(4,*) ' '        
      write(4,300) 'coordenada i', 'coordenada j' ! se que esta al reves dels indexs, però aqui es referix a s(i,j)
      do while(time.lt.timef)
            call metropolis(s,l,t,m,time)
      end do
      pas=pas+1
      do i=1,l
            do j=1,l
                 if(s(j,i).eq.1) then
                 write(4,400) j,i
                  end if
            end do
      end do
      write(4,*)' '
      write(4,*)' '
      close(4)

c temps = 400      
      timef=400.d0
      time=0.d0 
      m=mini
      pas=0
      s=s1

      open(5,file='valorst14.dat',status='unknown')
      open(2,file='magnettemp1.dat',status='unknown')
c      write(5,100) 'Calcul per a temperatura T=',t
c      write(5,100) 'Fins a temps t=',timef
c      write(5,*) ' '        
      write(5,300) 'coordenada i', 'coordenada j' ! se que esta al reves dels indexs, però aqui es referix a s(i,j)
      do while(time.lt.timef)
            call metropolis(s,l,t,m,time)
            write(2,200) time,m
      end do
      pas=pas+1
      do i=1,l
            do j=1,l
                 if(s(j,i).eq.1) then
                 write(5,400) j,i
                  end if
            end do
      end do
      write(5,*)' '
      write(5,*)' '
      close(5)
      close(2)


c part b)
      s1=1
      t=3.5d0
      mini=n

c temps = 0      
      timef=0.d0 
      time=0.d0
      m=mini
      pas=0
      s=s1

      open(7,file='valorst21.dat',status='unknown')
c      write(7,100) 'Calcul per a temperatura T=',t
c      write(7,100) 'Fins a temps t=',timef
c      write(7,*) ' '        
      write(7,300) 'coordenada i', 'coordenada j' ! se que esta al reves dels indexs, però aqui es referix a s(i,j)
      do while(time.lt.timef)
            call metropolis(s,l,t,m,time)
      end do
      pas=pas+1
      do i=1,l
            do j=1,l
                 if(s(j,i).eq.1) then
                 write(7,400) j,i
                  end if
            end do
      end do
      write(7,*)' '
      write(7,*)' '
      close(7)

c temps = 1      
      timef=1.d0
      time=0.d0 
      m=mini
      pas=0
      s=s1

      open(8,file='valorst22.dat',status='unknown')
c      write(8,100) 'Calcul per a temperatura T=',t
c      write(8,100) 'Fins a temps t=',timef
c      write(8,*) ' '        
      write(8,300) 'coordenada i', 'coordenada j' ! se que esta al reves dels indexs, però aqui es referix a s(i,j)
      do while(time.lt.timef)
            call metropolis(s,l,t,m,time)
      end do
      pas=pas+1
      do i=1,l
            do j=1,l
                 if(s(j,i).eq.1) then
                 write(8,400) j,i
                  end if
            end do
      end do
      write(8,*)' '
      write(8,*)' '
      close(8)

c temps = 10      
      timef=10.d0
      time=0.d0 
      m=mini
      pas=0
      s=s1

      open(9,file='valorst23.dat',status='unknown')
c      write(9,100) 'Calcul per a temperatura T=',t
c      write(9,100) 'Fins a temps t=',timef
c      write(9,*) ' '        
      write(9,300) 'coordenada i', 'coordenada j' ! se que esta al reves dels indexs, però aqui es referix a s(i,j)
      do while(time.lt.timef)
            call metropolis(s,l,t,m,time)
      end do
      pas=pas+1
      do i=1,l
            do j=1,l
                 if(s(j,i).eq.1) then
                 write(9,400) j,i
                  end if
            end do
      end do
      write(9,*)' '
      write(9,*)' '
      close(9)

c temps = 400      
      timef=400.d0
      time=0.d0 
      m=mini
      pas=0
      s=s1

      open(10,file='valors24.dat',status='unknown')
      open(6,file='magnettemp2.dat',status='unknown')
c      write(10,100) 'Calcul per a temperatura T=',t
c      write(10,100) 'Fins a temps t=',timef
c      write(10,*) ' '        
      write(10,300) 'coordenada i', 'coordenada j' ! se que esta al reves dels indexs, però aqui es referix a s(i,j)
      do while(time.lt.timef)
            call metropolis(s,l,t,m,time)
            write(6,200) time,m
      end do
      pas=pas+1
      do i=1,l
            do j=1,l
                 if(s(j,i).eq.1) then
                 write(10,400) j,i
                  end if
            end do
      end do
      write(10,*)' '
      write(10,*)' '
      close(10)
      close(6)

      stop
      end program

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c              SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

      subroutine metropolis(s,l,t,m,time)
      implicit none
      real*8 t,time,valor,prob
      integer l,s(l,l),m,n,dh,canvi,i,j,a,b,c,d

      n=l**2 !nombre de espins

c tirem numeros aleatoris amb distribució U(1,L) per a x i y, com tenim U(0,1) fem canvi y=a+(b-a)*x on x € U(0,1)
c d'aquesta manera tenim els dos índex (x,y) escollits de manera aleatòria a partir de la distribució U(1,L) per triar uniformement a l'atzar un espí 

      i=int(1.d0+(l-1)*rand()) 
      j=int(1.d0+(l-1)*rand())

c fem avançar el temps físic com t->t+(1/N)
      time=time+(1.d0/dble(n))

c porposem canvi espin s(i,j)=-s(i,j)      
      canvi=-s(i,j)

c condicions periòdiques, no tenim malla en 2D plana sino sobre un torus      
      if(j.eq.1) then
            a=s(i,l)
      else
            a=s(i,j-1)
      end if

      if(j.eq.l) then
            b=s(i,1)
      else
            b=s(i,j+1)
      end if

      if(i.eq.l) then
            c=s(1,j)
      else
            c=s(i+1,j)
      end if

      if(i.eq.1) then
            d=s(l,j)
      else
            d=s(i-1,j)
      end if

c variació energia associada a un possible canvi (dh=2*s(i,j)*(s(i,j-1)+s(i,j+1)+s(i+1,j)+s(i-1,j)))
      dh=2*s(i,j)*(a+b+c+d)

c  si canvi energia menor que zero acceptem canvi,si canvi energia major o igual que zero acceptem canvi amb una probabilitat e^-dh/t
c probabilitat acceptar canvi p=e^-dh/t, prob no acceptar 1-p
      if(dh.lt.0) then 
            s(i,j)=canvi
c si l'espí canvia de negatiu a positiu, tindrem un canvi positiu en la magnetització, si l'espí canvia de positiu a negatiu, tindrem un canvi negatiu en la magnetització
c com la magnetització absoluta és igual a n*magnetització de l'espí i n=l^2, si el canvi en la magnetització en l'espí és de m->m +/- (2/l^2), el canvi en la magnetització absoluta serà M->M +/- 2
            if(canvi.gt.0) m=m+2 
            if(canvi.lt.0) m=m-2 
      else              
            prob=dexp(-dh/t)
            valor=rand()
            if(valor.lt.prob) then 
                  s(i,j)=canvi
                  if(canvi.gt.0) m=m+2
                  if(canvi.lt.0) m=m-2
            end if
      end if

      return
      end subroutine

