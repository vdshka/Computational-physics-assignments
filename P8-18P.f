c -----------------------------------------------------------------------
C                      PRÀCTICA 8
c-----------------------------------------------------------------------
c resoldre l'eq de schrödinger independent del temps per a trobar els 
C autovalors i autovectors d'una partícula en potencial de pou quadrat 
c finit
C ----------------------------------------------------------------------
      program pre8
      implicit none
      integer nequs,npas,i,k,w,l,j
      real*8 v,sigma,e1,e2,e3,dx,e,x,cte,fun,long,xmin,xmax,beta
      real*8 ymax,ymin,dy,a1,a2,y,integral1,norm21(500),norm3(500)
      real*8 e3b,en23(3),en33(3),yyin3(2),yyout3(2),phi13,phi23,phi33
      real*8 e33,e23,e13,valor3(500)
      real*8 yyout(2),yyin(2),phi1,phi2,phi3,phi4,valor(500),bet(3)
      real*8 en1(4),en2(3),en3(3),norm(500),a,b,integral,norm2(500)
c      real*8 vec2(800),res
      common/dades/e,cte
      common/var/beta
      external fun,v

c-----------------FORMATS-----------------------------------------------
c ----------------------------------------------------------------------
100   format(a,2x,i10)
200   format(2x,2(e20.12,5x))
300   format(10x,a,15x,a)
400   format(12x,a,23x,a)
600   format(a,2x,e20.12)
c ---------------------------------------------------------------------

c inicialitzem vectors a zero

      norm=0.d0
      valor=0.d0
      norm2=0.d0
c      vec2=0.d0

c variables
      
      sigma=1.d-8
      nequs=2
      cte=3.80995d0 !eV A^2
      en1=(/-80.3d0,-79.5d0,-75.1d0,75.5d0/)
      npas=500
      long=7.d0
      xmin=-long/2.d0
      xmax=long/2.d0
      dx=(xmax-xmin)/dble(npas-1)
      ymin=-1.3d0 !posa y pero tots son valors de x !
      ymax=1.3d0
      dy=(ymax-ymin)/dble(npas-1)

      open(4,file='P8-18P-res1.dat',status='unknown')

      do w=1,4 ! bucle per a cada valor de energia
      	e=en1(w)
      	write(4,600) 'Energia',e
        write(4,*) ' '
        write(4,400) 'x','Phi(x)'
      	yyin=(/0.d0,1.d-5/)
        x=xmin
        do i=1,npas !bucle per a npas passos de runge kutta
      		call myRK4(x,dx,nequs,yyin,yyout)
      		do k=1,nequs !bucle per a canviar valors inicials per valors nous per a cada pas de runge kutta
      			yyin(k)=yyout(k)
      		end do
      		write(4,200) x,yyin(1)
      		x=x+dx
      	end do
      	write(4,*) ' '
        write(4,*) ' '
c ens quedara vector phi1= (phi1(50valors),phi1(800valors)),valors funcions a x=1
      	if (w.eq.1) phi1=yyin(1) 
      	if(w.eq.2)  phi2=yyin(1)
      	if (w.eq.3) phi3=yyin(1) 
      	if (w.eq.4) phi4=yyin(1) 
      end do

      close(4)

c apartat 2-mètode de tir
       
      en2=(/-80.3d0,-75.1d0,-70.d0/)
      en3=(/-79.5d0,-75.5d0,-71.d0/)
      en23=(/-80.3d0,-75.1d0,-70.d0/)
      en3=(/-79.5d0,-75.5d0,-71.d0/)
      bet=(/0.d0,1.d0,5.d0/)
      

c calculem autovalors i autovectors

      open(2,file='Fig.dat',status='unknown')
      open(1,file='P8-18P-res3.dat',status='unknown')
      open(3,file='P8-18P-res2.dat',status='unknown')
      open (4,file='Probabilitat.dat',status='unknown')
      write(1,100)'Nombre de passos:',npas
      write(2,100)'Nombre de passos:',npas
      write(3,100)'Nombre de passos:',npas
      write(3,*) 'Energia i phi(x=1) per a cada pas de la secant pels qu
     +atre autovalors'
      write(3,*) ' ' 

      do j=1,3
      	beta=bet(j)
      	write(1,600) 'Beta:',beta
        write(4,600) 'Beta:',beta
        write(4,*) ' '
        write(1,*) ' ' 
        write(1,300) 'Autovalor','Autovector'
        write(2,600) 'Beta:',beta
        write(2,*) ' '     
        write(3,600) 'Beta:',beta
        write(3,*) ' ' 
      	do l=1,3 ! bucle per a cada valor de energia
      		e1=en2(l)
      		e2=en3(l)
          e13=en23(l)
          e23=en33(l)
      		do w=1,2 !bucle e1 i e2
        		if(w.eq.1) e=e1
            if(w.eq.2) e=e2 
            if(w.eq.1) e3=e13
            if(w.eq.2) e3=e23 
            yyin=(/0.d0,1.d-5/)
            yyin3=(/0.d0,1.d-5/)
            x=xmin
      			do i=1,npas !bucle per a npas passos de runge kutta
      		    	call myRK4(x,dx,nequs,yyin,yyout)
                call myRK4(x,dx,nequs,yyin3,yyout3)
      		        do k=1,nequs !bucle per a canviar valors inicials per valors nous per a cada pas de runge kutta
      		        	yyin(k)=yyout(k)
                    yyin3(k)=yyout3(k)
      		        end do
      		        x=x+dx
      		    end do
c ens quedara vector phi1= (phi1(50valors),phi1(800valors)),valors funcions a x=1
      		    if (w.eq.1) phi1=yyin(1) 
      		    if(w.eq.2)  phi2=yyin(1)
              if (w.eq.1) phi13=yyin3(1) 
              if(w.eq.2)  phi23=yyin3(1)
         	end do!-----aqui tenemos e1,e2,phi1,phi2
c calculem energia e3 e3=(e1*phi2-e2phi1)/(phi2-phi1) i runge kutta per a nova energia
c repetir fins que dabs(phi3)<sigma, llavors considerem que hem convergit
         	phi3=1.d0  
          phi33=1.d0
         	write(3,300) 'Energia','Phi(x=1)'                 
         	do while(dabs(phi3).gt.sigma) !bucle que es realitza fins que hem convergit/trobat autovalor
         		e3=(e1*phi2-e2*phi1)/(phi2-phi1)
            e33=(e13*phi23-e23*phi13)/(phi23-phi13)
            	e=e3
              e3b=e33
            	yyin=(/0.d0,1.d-5/)
              yyin3=(/0.d0,1.d-5/)
            	x=xmin
            	do i=1,npas !bucle per a npas passos de runge kutta
            		call myRK4(x,dx,nequs,yyin,yyout)
                call myRK4(x,dx,nequs,yyin3,yyout3)
c bucle per a canviar valors inicials per valors nous per a cada pas de runge kutta 
                	do k=1,nequs 
               			yyin(k)=yyout(k)
                    yyin3(k)=yyout3(k)
                	end do                              
                	x=x+dx
                	norm(i)=(yyout(1))**2
                	valor(i)=yyin(1)
                  norm3(i)=(yyout3(1))**2
                  valor3(i)=yyin3(1)
             	end do
             	phi3=yyin(1)
              phi33=yyin3(1)
c canviem valors per a trobar seguent valor energia
              	e1=e2
              	e2=e3
              	phi1=phi2
              	phi2=phi3
                e13=e23
                e23=e33
                phi13=phi23
                phi23=phi33
              	write(3,200) e,phi3
c -----------------normalitzacio                        
              	a=xmin
              	b=xmax
                a1=ymin
                a2=ymax
              	call simpson(a,b,npas,norm,integral)
                call simpson(a1,a2,npas,norm3,integral1)
c              	do i =1,npas   !comprobacio, res ha de ser =1
c                 	vec2(i)=norm(i)/integral
c              	end do
c              	call simpson(a,b,npas,vec2,res)
c              	write(*,*) res                        
         	end do 
         	write(3,*) ' '
         	write(3,*) ' '                
         	write(1,200) e3,phi3/dsqrt(integral)
c escrivim les dades de phi(x,E) en funcio de la x al fitxer
         	write(2,600) 'Energia:',e
         	write(2,*) ' '
         	write(2,400) 'x','Phi'
          write(4,600) 'Energia:',e3b
c          write(*,*) e
          write(4,*) ' '
          write(4,400) 'x','Phi'
         	x=xmin
          y=ymin
         	do i=1,npas
         		norm2(i)=valor(i)/dsqrt(integral)
            norm21(i)=norm3(i)/dsqrt(integral1)
            	write(2,200) x,norm2(i)
              write(4,200) y,norm21(i)
            	x=x+dx
              y=y+dy
         	end do 
         	write(2,*) ' '
         	write(2,*) ' '
          write(4,*) ' '
          write(4,*) ' '
      	end do
      	write(1,*) ' '
      	write(1,*) ' '
      end do

      close(1)
      close(2)
      close(3)
      close(4)

      stop
      end program

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                         SUBRUTINES
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

c subrutina que calcula un pas del metode runge kutta 4 per a un sistema de 
c nequs eq de 1er ordre acoblades

c     ·Runge kutta 4t ordre:
c         k1=h*f(x0+y0)
c         k2=h*f(x0+0.5*h,y0+0.5*k1)
c         k3=h*f(x0+0.5*h,y0+0.5*k2)
c         k4=h*f(x0+h,y0+k3)
c      y(x0+h)=y0+(k1+2*k2+2*k3+k4)

c aqui h=dx,f(x,y)=dy/dx,y=y'',x=x0,yyin=y0 vector

      subroutine myRK4(x,dx,nequs,yyin,yyout)
      implicit none
      real*8 yyin(nequs),yyout(nequs),yin(nequs),x,dx
      real*8 k1(nequs),k2(nequs),k3(nequs),k4(nequs),ksuma(nequs)
      integer nequs,i

c inicialitzem vectors a 0
      yyout=0.d0
      k1=0.d0
      k2=0.d0
      k3=0.d0
      k4=0.d0
      ksuma=0.d0
      yin=0.d0
      
c----obtencio k1--------  
      call derivades (nequs,x,yyin,k1)
c---obtencio k2-----------
      do i=1,nequs
      	yin(i)=yyin(i)+(dx*0.5d0*k1(i))
      end do

      call derivades (nequs,x+0.5d0*dx,yin,k2)
c---obtencio k3----------- 
      do i=1,nequs
      	yin(i)=yyin(i)+(dx*0.5d0*k2(i))
      end do

      call derivades (nequs,x+0.5d0*dx,yin,k3)
c---obtencio k4-----------      
      do i=1,nequs
      	yin(i)=yyin(i)+(dx*k3(i))
      end do

      call derivades (nequs,x+dx,yin,k4)

c-----calcul vector yyout    
      do i=1,nequs
      	ksuma(i)=k1(i) + 2.d0*k2(i) + 2.d0*k3(i)+k4(i)
      	yyout(i)=yyin(i)+(dx/6.d0)*ksuma(i)
      end do

      return
      end subroutine

c subrutina que en donar x i un vector yin retorna el valos dyin/dx dins del vector
c dyout. Subrutina especialitzada per a la eq a resoldre
c yin=(funcio,der 1era funcio), dyout=(der funcio, der segona funcio)

      subroutine derivades(nequ,x,yin,dyout)
      implicit none
      real*8 x,yin(nequ),dyout(nequ),e,v,fun
      integer nequ,i
      common/dades/e,v

      dyout=0.d0 !inicialitzem vector a 0

c per a sist de nequ eq de primer ordre acoblades, y es la dy de l'eq anterior
      do i=1,nequ-1
      	dyout(i)=yin(i+1)
      end do
      
c      dyout(nequ)=(-2.d0*(e-v))*yin(nequ-1) ! aïllem segona derivada equacio
      dyout(nequ)=fun(x)*yin(nequ-1)

      return
      end subroutine

c subrutina que implementa el metode de simpson compost
c retorna integral, que es la aproximacio obtinguda per a la integral emprant aquest metode
c la utilitzarem per a calcular norma

      subroutine simpson(a,b,npas,fcn,integral)
      implicit none
      real*8 h,x,a,b,fcn(500),i,integral
      integer j,npas


      i=fcn(1)+fcn(npas)
      h=(b-a)/dble(npas-1)
      do j=1,npas
        x=a+dble(j)*h
        if (mod(j,2).eq.0) i=i+2.d0*fcn(j)
        if (mod(j,2).ne.0) i=i+4.d0*fcn(j)
      end do  
      integral=i*h/3.d0

      return
      end subroutine

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c                          FUNCIONS      
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------

      real*8 function fun(x)
      implicit none
      real*8 v,e,cte,x,temp
      common/dades/e,cte
      
      temp=(-1.d0/cte)*(e-v(x))
      fun=temp

      end function

c---- Funcio del potencial
      real*8 function V1(x)
      implicit none
      real*8 x,fun
      
      if (dabs(x).le.3.d0) fun=-80.d0
      if(dabs(x).gt.3.d0) fun=0.d0
      
      V1=fun
      
      end function


c-----funcio potencial nou
      real*8 function v(x)
      implicit none
      real*8 v1,beta,x,temp
      common/var/beta

      temp=v1(x)+beta*(x**2)
      v=temp

      end function

c-----------------------------------------------------------------------