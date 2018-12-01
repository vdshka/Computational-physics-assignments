
C                      ASSIGNMENT 2

c Pistons position depending on time. Interpolation.

      PROGRAM PISTO
      IMPLICIT NONE
      REAL*8 radit1,w0,L,phi,t,r,m,n,h,augment,tin,xout,x1
      REAL*8 x(5),posis(0:500),temps(0:500)
      INTEGER j
      COMMON/DADES/temps,posis

c               FORMATS

100   FORMAT(6(E20.14,4X))
200   FORMAT(3(E20.14,4X))


c THE POSITIONS OF FIVE PISTONS DEPENDING ON TIME ARE GIVEN BY THE EQUATION:
C Xk(t)= Rk*cos(w0t+phik)+dqrt(L^2-Rk^2*sin^2(w0t+phik))
c WHERE W0 IS THE FREQUENCY,RK THE RADIUM OF THE CRANK, PHIK THE PHASE AND K THE NUMBER OF THE PISTON

      w0=5.d0 
      L=18.5d0 !LENGTH OF THE CONNECTING ROD
      
C SAVES IN A FILE THE VALUES OF TIME AND POSITION OF THE FOUR PISTONS

      OPEN(1,file='P2-18P-res1.dat')

      t=0.d0

      DO WHILE (t.LE.5)
        CALL posit1(w0,L,t,x)
        WRITE(1,100) t,x(1:5)
        t=t+0.01d0
      END DO
      CLOSE(1) 

C READS THE FILE CREATED BEFORE AND CREATES TWO VECTORS POSIS AND TEMPS WITH THE VALUES OF THE POSITIONS OF THE 
C FOURTH PISTON AND THE TIME RESPECTIVELY

      OPEN(2,file='P2-18P-res1.dat',status='old')
c      OPEN(4,file='miau.dat')
      DO j=0,500
            READ(2,*)temps(j),r,m,n,posis(j),h
c            WRITE(4,*)temps(j),posis(j)
      END DO 
      CLOSE(2) 
c      CLOSE(4)

C SAVES IN A FILE THE VALUES OF TIME AND POSITION OF THE FOURTH PISTON OBTAINED BY INTERPOLATION

      OPEN(3,file='P2-18P-res2.dat')
      augment=(3.d0/2000.d0)
      tin=0.d0
      DO WHILE(tin.LE.3)
        CALL xinterpol(tin,xout)
        x1=xout
        CALL xinterpol0(tin,xout)
        WRITE(3,200) tin,xout,x1
        tin=tin+augment
      END DO
      CLOSE(3)

      STOP
      END PROGRAM


c                         FUNCTIONS


C IN: LENGTH AND NUMBER OF THE PISTON
C OUT: RADIUS OF THE PISTON K IN CM    

      REAL*8 FUNCTION radit1(L,k)
      IMPLICIT NONE
      REAL*8 L,temp
      INTEGER k

      temp=(L/dble(k))-0.5d0
      radit1=temp

      RETURN
      END FUNCTION

C FUNCTION THAT COMPUTES THE INITIAL PHASE

      REAL*8 FUNCTION phi(k)
      IMPLICIT NONE
      REAL*8 pi,temp1
      INTEGER k

      pi=dacos(-1.0d0)
      temp1=pi*(dble(k)/5.d0)**2
      phi=temp1

      RETURN
      END FUNCTION


c                         SUBROUTINES


C RETURNS THE POSITION OF THE FIVE PISTONS IN A VECTOR X

      SUBROUTINE posit1(w0,l,t,x)
      IMPLICIT NONE
      REAL*8 x(5),phi,L,t,w0,radit1,s1,s2
      INTEGER k

      DO k=1,5
        s1=radit1(L,k)*dcos(w0*t+phi(k))
        s2=dsqrt(L**2-(radit1(L,k)**2)*dsin(w0*t+phi(k))**2)
        x(k)=s1+s2
      END DO

      RETURN
      END SUBROUTINE

C LINEAR INTERPOLATION
    
      SUBROUTINE xinterpol(tin,xout)
      IMPLICIT NONE
      REAL*8 tin,xout,posis(0:500),temps(0:500),xt,div
      INTEGER a0,a1
      COMMON/DADES/temps,posis

      a0=int(tin/0.01d0)
      a1=a0+1

      div=(posis(a1)-posis(a0))*(tin-temps(a0))/(temps(a1)-temps(a0))
      xt=posis(a0)+div

      xout=xt

      RETURN
      END SUBROUTINE     

C ZERO ORDER INTERPOLATION
    
      SUBROUTINE xinterpol0(tin,xout)
      IMPLICIT NONE
      REAL*8 tin,xout,posis(0:500),temp2,temps(0:500),t0,t1
      INTEGER i
      COMMON/DADES/temps,posis

      DO i=0,499
        t0=temps(i)
        t1=temps(i+1)
        IF((t0.LE.tin).AND.(t1.GE.tin)) temp2=posis(i)
      END DO

      xout=temp2

      RETURN
      END SUBROUTINE     

