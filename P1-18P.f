
C                      ASSIGNMENT 1
c  Sum of the function Pk


      PROGRAM PRACTICA

      IMPLICIT NONE
      REAL*8 div,j,p,suma,s,asimpt
      INTEGER n

C READS AN INTEGER K AND RETURNS PK, DISPLAYS THE VALUE
c IF THE NUMBER IS NOT BETWEEN 3 AND 35, ASKS AGAIN
c IF THE NUMBER IS NOT AN INTEGER, ASKS AGAIN       

      j = 0

      DO WHILE ((j.LT.15).OR.(j.GT.221).OR.(INT(j).NE.j))
            WRITE (*,*) 'Enter an integer between 3 and 35'
            READ  (*,*) j
            IF ((j.LT.15).OR.(j.GT.221)) THEN
                  WRITE(*,*)'This number is not in the interval'
            ELSE IF (INT(j).NE.j) THEN
                  WRITE(*,*) 'This number is not an integer'
            ENDIF
      END DO
      WRITE(*,*) 'Value of Pk', P(INT(j))

C COMPUTES THE SUM FROM 28 TO 65

      WRITE (*,*) 'Value of the sum from 28 to 65', SUMA(28,65)

C CREATES A FILE (1) WITH THREE COLUMNS
C FIRST COLUMN: VALUES N BETWEEN 11 AND 311 IN THREES
C SECOND COLUMN: VALUES OF THE SUM FROM 8 TO N
C THIRD COLUMN: ASYMPTOTIC BEHAVIOUR 1/5*N^3

C CREATES A FILE (2) WITH TWO COLUMNS: N, DIVISION SUM FROM 8 TO N/ASYMPTOTIC BEHAVIOUR

      
      OPEN(1,FILE='P1-18P-res1.dat')
      OPEN(2,FILE='P1-18P-res2.dat')

      DO n = 11,311,3
            s = SUMA(8,n)
            asimpt = (1.0d0/5.0d0)*(n**3)
            div = s/asimpt
            WRITE (1,*) n,s,asimpt 
            WRITE (2,*) n,div     
      END DO
      CLOSE(1)
      CLOSE(2)

      STOP
      END PROGRAM 


c                         FUNCTIONS


C FUNCTION THAT COMPUTES PK 

      REAL*8 FUNCTION p(k)
      IMPLICIT NONE
      INTEGER k
      REAL*8 e

      e = dexp(1.0d0)

      p = (3.0d0/5.0d0)*(k**2)+e+(10.0d0*k)
      RETURN
      END FUNCTION

C FUNCTION THAT COMPUTES THE SUM

      REAL*8 FUNCTION SUMA(n1,n2)
      IMPLICIT NONE
      INTEGER n1,n2,k
      REAL*8 p,temp

      temp = 0.0d0
      DO k = n1,n2
            temp = temp + p(k)
      END DO
      suma=temp
      RETURN
      END FUNCTION      
 
