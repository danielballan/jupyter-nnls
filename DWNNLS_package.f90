      program call_wnnls
! Main routine to set up nonnegative linear least squares problem
! Note that unconstrained solution is given by Lapack DGELSS including
! handling of rank deficient case for linear least squares

      implicit none
      integer mdw, me, ma, n, l
      double precision, allocatable:: w(:,:),work(:),x(:),a(:,:)
      double precision, allocatable:: best(:),y(:)
      double precision, allocatable:: au(:,:), xu(:),s(:)

      integer, allocatable:: iwork(:)
      double precision  prgopt(4)
      double precision  rnorm,rcond
      integer mode,rank,info,lwork
      integer i,j
      character name*5(10)

      me=0
      l=0
! First read the dimensions of the problem ma,n
! ma represents number of components, n represents number of raw materials
      read (5,*) ma,n
      mdw=me+ma
! This definition of lwork handles both dwnnls and dgelss (unconstrained)
      lwork = max(3*min(ma,n)+max(2*min(ma,n),max(ma,n),1),me+ma+5*n)
      allocate ( w(mdw,n+1) )
      allocate ( x(n) )
      allocate ( work(lwork))
      allocate ( iwork(me+ma+n) )
      allocate ( best(ma) )
      allocate ( y(ma))
      allocate ( a(ma,n) )
      allocate ( au(ma,n) ) ! Use au and xu to solve unconstrained problem
      allocate ( xu(n) )   ! for comparison and testing
      allocate ( s(ma) )   ! s holds singular values for dgelss

! prgopt with this definition should normalize each column of A to 1.0
!     prgopt(1)=4
!     prgopt(2)=6.0d0
!     prgopt(3)=1.0d0
!     prgopt(4)=1.0d0
! Or don't normalize so retention factors can be used
      prgopt(1)=1.0d0

      iwork(1)=lwork
      iwork(2)=me+ma+n

! Define the matrix A(ma by n) and target vector B(ma) 
! The first n columns of w contain columns of A (of length ma) and
! the last column of w contains B (also of length ma).

! Zero input A should cause error message? Maybe not because of least
! squares sense of solution; singular matrix is ok
! Replace hard-coded input with reading of rows of matrix from 
! stdin text file
      do i=1,ma
       read (5,*) (w(i,j),j=1,n) 
      enddo
! Read the target composition as last row (for easier editing)
      read (5,*) (w(i,n+1),i=1,ma)
 2030 format("Target: ",12f10.5)
! Keep copies of starting arrays
      a=w(1:ma,1:n)
      y=w(1:ma,n+1)
! Keep an additional copy that will be used in unconstrained solver
      au = a
      xu = y(1:ma)
      write (6,2030) y

      call dwnnls(w,mdw,me,ma,n,l,prgopt,x,rnorm,mode,iwork,work)

      write(6,2000) rnorm,mode
 2000 format("rnorm is ",e12.4," with return mode ",i5)

! Also solve unconstrained linear least squares (lapack call)
      call dgelss(ma,n,1,au,ma,xu,n,s,rcond,rank,work,lwork,info)
      write(6,2002) rcond,rank,info
 2002 format("unconstrained cond, rank, info ",e12.4,i5,i5)
!     write(6,2003) (s(i),i=1,ma)
!2003 format("s",10f10.4)

      name(1)="SiO2"
      name(2)="Al2O3"
      name(3)="Al2O3"
      name(4)="Na2O"
      name(5)="Fe2O3"
      name(6)="SiO2"
      name(7)="BaO"
      name(8)="SnO2"
      name(9)="Fe2O3"
      name(10)="ZrO2"

      write(6,2004)
 2004 format("Constrained and unconstrained batching solutions")
      do i=1,n
       write(6,2010) i,x(i),xu(i)
      enddo
 2010 format(i5,2f12.5)

! Also calculate the composition associated with best solution
      best=matmul(a,x)
      s=matmul(a,xu)
      write(6,2015)
 2015 format(3x,"# comp     target    constr   unconst  constr-target")
      do i=1,ma
       write(6,2020) i,name(i),y(i),best(i),s(i),best(i)-y(i)
      enddo
 2020 format(i5," ",a5," ",3f12.5,f12.6)

      deallocate ( w, x, work, iwork, best, y, a, au, xu, s)

      end
*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (I)
C***BEGIN PROLOGUE  D1MACH
C***PURPOSE  Return floating point machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   D1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument, and can be referenced as follows:
C
C        D = D1MACH(I)
C
C   where I=1,...,5.  The (output) value of D above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   Assume double precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(14) = T, the number of base-B digits.
C   I1MACH(15) = EMIN, the smallest exponent E.
C   I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890213  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   900911  Added SUN 386i constants.  (WRB)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added CONVEX -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C***END PROLOGUE  D1MACH
C
      DOUBLE PRECISION SMALL(1)
      DOUBLE PRECISION LARGE(1)
      DOUBLE PRECISION RIGHT(1)
      DOUBLE PRECISION DIVER(1)
      DOUBLE PRECISION LOG10(1)
C
      DOUBLE PRECISION DMACH(5)
      SAVE DMACH
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA SMALL(1) / ZC00800000 /
C     DATA SMALL(2) / Z000000000 /
C     DATA LARGE(1) / ZDFFFFFFFF /
C     DATA LARGE(2) / ZFFFFFFFFF /
C     DATA RIGHT(1) / ZCC5800000 /
C     DATA RIGHT(2) / Z000000000 /
C     DATA DIVER(1) / ZCC6800000 /
C     DATA DIVER(2) / Z000000000 /
C     DATA LOG10(1) / ZD00E730E7 /
C     DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O0000000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O0007777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O7770000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O7777777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA SMALL(1) / Z"3001800000000000" /
C     DATA SMALL(2) / Z"3001000000000000" /
C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
C     DATA LARGE(2) / Z"4FFE000000000000" /
C     DATA RIGHT(1) / Z"3FD2800000000000" /
C     DATA RIGHT(2) / Z"3FD2000000000000" /
C     DATA DIVER(1) / Z"3FD3800000000000" /
C     DATA DIVER(2) / Z"3FD3000000000000" /
C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA SMALL(1) / 00564000000000000000B /
C     DATA SMALL(2) / 00000000000000000000B /
C     DATA LARGE(1) / 37757777777777777777B /
C     DATA LARGE(2) / 37157777777777777777B /
C     DATA RIGHT(1) / 15624000000000000000B /
C     DATA RIGHT(2) / 00000000000000000000B /
C     DATA DIVER(1) / 15634000000000000000B /
C     DATA DIVER(2) / 00000000000000000000B /
C     DATA LOG10(1) / 17164642023241175717B /
C     DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn OR -pd8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CC0000000000000' /
C     DATA DMACH(4) / Z'3CD0000000000000' /
C     DATA DMACH(5) / Z'3FF34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3F900000000000000000000000000000' /
C     DATA DMACH(4) / Z'3F910000000000000000000000000000' /
C     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
C
C     MACHINE CONSTANTS FOR THE CRAY
C
C     DATA SMALL(1) / 201354000000000000000B /
C     DATA SMALL(2) / 000000000000000000000B /
C     DATA LARGE(1) / 577767777777777777777B /
C     DATA LARGE(2) / 000007777777777777774B /
C     DATA RIGHT(1) / 376434000000000000000B /
C     DATA RIGHT(2) / 000000000000000000000B /
C     DATA DIVER(1) / 376444000000000000000B /
C     DATA DIVER(2) / 000000000000000000000B /
C     DATA LOG10(1) / 377774642023241175717B /
C     DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC DMACH(5)
C
C     DATA SMALL /    20K, 3*0 /
C     DATA LARGE / 77777K, 3*177777K /
C     DATA RIGHT / 31420K, 3*0 /
C     DATA DIVER / 32020K, 3*0 /
C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA DMACH(1) / '0000000000000010'X /
C     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
C     DATA DMACH(3) / '0000000000003CC0'X /
C     DATA DMACH(4) / '0000000000003CD0'X /
C     DATA DMACH(5) / '79FF509F44133FF3'X /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FORMAT
C
C     DATA DMACH(1) / '0010000000000000'X /
C     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X /
C     DATA DMACH(3) / '3CA0000000000000'X /
C     DATA DMACH(4) / '3CB0000000000000'X /
C     DATA DMACH(5) / '3FD34413509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
C     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
C     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
C     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /        128,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C     DATA DIVER(1), DIVER(2) /       9472,           0 /
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /         16,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C     DATA DIVER(1), DIVER(2) /      15568,           0 /
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) /  40000B,       0 /
C     DATA SMALL(3), SMALL(4) /       0,       1 /
C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
C     DATA RIGHT(3), RIGHT(4) /       0,    225B /
C     DATA DIVER(1), DIVER(2) /  40000B,       0 /
C     DATA DIVER(3), DIVER(4) /       0,    227B /
C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
C
      DATA SMALL(1) / 2.23D-308  /
      DATA LARGE(1) / 1.79D+308  /
      DATA RIGHT(1) / 1.11D-16   /
      DATA DIVER(1) / 2.22D-16   /
      DATA LOG10(1) / 0.301029995663981195D0 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    8388608,           0 /
C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
C     DATA DIVER(1), DIVER(2) /  620756992,           0 /
C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
C
C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    128,      0 /
C     DATA SMALL(3), SMALL(4) /      0,      0 /
C     DATA LARGE(1), LARGE(2) /  32767,     -1 /
C     DATA LARGE(3), LARGE(4) /     -1,     -1 /
C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
C     DATA RIGHT(3), RIGHT(4) /      0,      0 /
C     DATA DIVER(1), DIVER(2) /   9472,      0 /
C     DATA DIVER(3), DIVER(4) /      0,      0 /
C     DATA LOG10(1), LOG10(2) /  16282,   8346 /
C     DATA LOG10(3), LOG10(4) / -31493, -12296 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA SMALL(3), SMALL(4) / O000000, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA LARGE(3), LARGE(4) / O177777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
C     DATA DIVER(1), DIVER(2) / O022400, O000000 /
C     DATA DIVER(3), DIVER(4) / O000000, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020232 /
C     DATA LOG10(3), LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
C     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
C     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
C     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
C
C     MACHINE CONSTANTS FOR THE SUN 386i
C
C     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
C     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
C     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'D1MACH',
     +   'I OUT OF BOUNDS', 1, 2)
C
      D1MACH = DMACH(I)
      RETURN
C
      END
*DECK DH12
      SUBROUTINE DH12 (MODE, LPIVOT, L1, M, U, IUE, UP, C, ICE, ICV,
     +   NCV)
C***BEGIN PROLOGUE  DH12
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DHFTI, DLSEI and DWNNLS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (H12-S, DH12-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C      *** DOUBLE PRECISION VERSION OF H12 ******
C
C     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
C     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
C
C     Construction and/or application of a single
C     Householder transformation..     Q = I + U*(U**T)/B
C
C     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
C     LPIVOT is the index of the pivot element.
C     L1,M   If L1 .LE. M   the transformation will be constructed to
C            zero elements indexed from L1 through M.   If L1 GT. M
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry to H1 U() contains the pivot vector.
C                   IUE is the storage increment between elements.
C                                       On exit from H1 U() and UP
C                   contain quantities defining the vector U of the
C                   Householder transformation.   On entry to H2 U()
C                   and UP should contain quantities previously computed
C                   by H1.  These will not be modified by H2.
C     C()    On entry to H1 or H2 C() contains a matrix which will be
C            regarded as a set of vectors to which the Householder
C            transformation is to be applied.  On exit C() contains the
C            set of transformed vectors.
C     ICE    Storage increment between elements of vectors in C().
C     ICV    Storage increment between vectors in C().
C     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0
C            no operations will be done on C().
C
C***SEE ALSO  DHFTI, DLSEI, DWNNLS
C***ROUTINES CALLED  DAXPY, DDOT, DSWAP
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900911  Added DDOT to DOUBLE PRECISION statement.  (WRB)
C***END PROLOGUE  DH12
      INTEGER I, I2, I3, I4, ICE, ICV, INCR, IUE, J, KL1, KL2, KLP,
     *     L1, L1M1, LPIVOT, M, MML1P2, MODE, NCV
      DOUBLE PRECISION B, C, CL, CLINV, ONE, UL1M1, SM, U, UP, DDOT
      DIMENSION U(IUE,*), C(*)
C     BEGIN BLOCK PERMITTING ...EXITS TO 140
C***FIRST EXECUTABLE STATEMENT  DH12
         ONE = 1.0D0
C
C     ...EXIT
         IF (0 .GE. LPIVOT .OR. LPIVOT .GE. L1 .OR. L1 .GT. M) GO TO 140
         CL = ABS(U(1,LPIVOT))
         IF (MODE .EQ. 2) GO TO 40
C           ****** CONSTRUCT THE TRANSFORMATION. ******
            DO 10 J = L1, M
               CL = MAX(ABS(U(1,J)),CL)
   10       CONTINUE
            IF (CL .GT. 0.0D0) GO TO 20
C     .........EXIT
               GO TO 140
   20       CONTINUE
            CLINV = ONE/CL
            SM = (U(1,LPIVOT)*CLINV)**2
            DO 30 J = L1, M
               SM = SM + (U(1,J)*CLINV)**2
   30       CONTINUE
            CL = CL*SQRT(SM)
            IF (U(1,LPIVOT) .GT. 0.0D0) CL = -CL
            UP = U(1,LPIVOT) - CL
            U(1,LPIVOT) = CL
         GO TO 50
   40    CONTINUE
C        ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C
         IF (CL .GT. 0.0D0) GO TO 50
C     ......EXIT
            GO TO 140
   50    CONTINUE
C     ...EXIT
         IF (NCV .LE. 0) GO TO 140
         B = UP*U(1,LPIVOT)
C        B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C
         IF (B .LT. 0.0D0) GO TO 60
C     ......EXIT
            GO TO 140
   60    CONTINUE
         B = ONE/B
         MML1P2 = M - L1 + 2
         IF (MML1P2 .LE. 20) GO TO 80
            L1M1 = L1 - 1
            KL1 = 1 + (L1M1 - 1)*ICE
            KL2 = KL1
            KLP = 1 + (LPIVOT - 1)*ICE
            UL1M1 = U(1,L1M1)
            U(1,L1M1) = UP
            IF (LPIVOT .NE. L1M1) CALL DSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
            DO 70 J = 1, NCV
               SM = DDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE)
               SM = SM*B
               CALL DAXPY(MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE)
               KL1 = KL1 + ICV
   70       CONTINUE
            U(1,L1M1) = UL1M1
C     ......EXIT
            IF (LPIVOT .EQ. L1M1) GO TO 140
            KL1 = KL2
            CALL DSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
         GO TO 130
   80    CONTINUE
            I2 = 1 - ICV + ICE*(LPIVOT - 1)
            INCR = ICE*(L1 - LPIVOT)
            DO 120 J = 1, NCV
               I2 = I2 + ICV
               I3 = I2 + INCR
               I4 = I3
               SM = C(I2)*UP
               DO 90 I = L1, M
                  SM = SM + C(I3)*U(1,I)
                  I3 = I3 + ICE
   90          CONTINUE
               IF (SM .EQ. 0.0D0) GO TO 110
                  SM = SM*B
                  C(I2) = C(I2) + SM*UP
                  DO 100 I = L1, M
                     C(I4) = C(I4) + SM*U(1,I)
                     I4 = I4 + ICE
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
*DECK DWNLIT
      SUBROUTINE DWNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE,
     +   RNORM, IDOPE, DOPE, DONE)
C***BEGIN PROLOGUE  DWNLIT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DWNNLS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLIT-S, DWNLIT-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to DWNNLS( ).
C     The documentation for DWNNLS( ) has complete usage instructions.
C
C     Note  The M by (N+1) matrix W( , ) contains the rt. hand side
C           B as the (N+1)st col.
C
C     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with
C     col interchanges.
C
C***SEE ALSO  DWNNLS
C***ROUTINES CALLED  DCOPY, DH12, DROTM, DROTMG, DSCAL, DSWAP, DWNLT1,
C                    DWNLT2, DWNLT3, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900604  DP version created from SP version. .  (RWC)
C***END PROLOGUE  DWNLIT
      INTEGER IDOPE(*), IPIVOT(*), ITYPE(*), L, M, MDW, N
      DOUBLE PRECISION DOPE(*), H(*), RNORM, SCALE(*), W(MDW,*)
      LOGICAL DONE
C
      EXTERNAL DCOPY, DH12, DROTM, DROTMG, DSCAL, DSWAP, DWNLT1,
     *   DWNLT2, DWNLT3, IDAMAX
      INTEGER IDAMAX
      LOGICAL DWNLT2
C
      DOUBLE PRECISION ALSQ, AMAX, EANORM, FACTOR, HBAR, RN, SPARAM(5),
     *   T, TAU
      INTEGER I, I1, IMAX, IR, J, J1, JJ, JP, KRANK, L1, LB, LEND, ME,
     *   MEND, NIV, NSOLN
      LOGICAL INDEP, RECALC
C
C***FIRST EXECUTABLE STATEMENT  DWNLIT
      ME    = IDOPE(1)
      NSOLN = IDOPE(2)
      L1    = IDOPE(3)
C
      ALSQ   = DOPE(1)
      EANORM = DOPE(2)
      TAU    = DOPE(3)
C
      LB     = MIN(M-1,L)
      RECALC = .TRUE.
      RNORM  = 0.D0
      KRANK  = 0
C
C     We set FACTOR=1.0 so that the heavy weight ALAMDA will be
C     included in the test for column independence.
C
      FACTOR = 1.D0
      LEND = L
      DO 180 I=1,LB
C
C        Set IR to point to the I-th row.
C
         IR = I
         MEND = M
         CALL DWNLT1 (I, LEND, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,
     +                W)
C
C        Update column SS and find pivot column.
C
         CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C        Perform column interchange.
C        Test independence of incoming column.
C
  130    IF (DWNLT2(ME, MEND, IR, FACTOR, TAU, SCALE, W(1,I))) THEN
C
C           Eliminate I-th column below diagonal using modified Givens
C           transformations applied to (A B).
C
C           When operating near the ME line, use the largest element
C           above it as the pivot.
C
            DO 160 J=M,I+1,-1
               JP = J-1
               IF (J.EQ.ME+1) THEN
                  IMAX = ME
                  AMAX = SCALE(ME)*W(ME,I)**2
                  DO 150 JP=J-1,I,-1
                     T = SCALE(JP)*W(JP,I)**2
                     IF (T.GT.AMAX) THEN
                        IMAX = JP
                        AMAX = T
                     ENDIF
  150             CONTINUE
                  JP = IMAX
               ENDIF
C
               IF (W(J,I).NE.0.D0) THEN
                  CALL DROTMG (SCALE(JP), SCALE(J), W(JP,I), W(J,I),
     +                         SPARAM)
                  W(J,I) = 0.D0
                  CALL DROTM (N+1-I, W(JP,I+1), MDW, W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  160       CONTINUE
         ELSE IF (LEND.GT.I) THEN
C
C           Column I is dependent.  Swap with column LEND.
C           Perform column interchange,
C           and find column in remaining set with largest SS.
C
            CALL DWNLT3 (I, LEND, M, MDW, IPIVOT, H, W)
            LEND = LEND - 1
            IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
            HBAR = H(IMAX)
            GO TO 130
         ELSE
            KRANK = I - 1
            GO TO 190
         ENDIF
  180 CONTINUE
      KRANK = L1
C
  190 IF (KRANK.LT.ME) THEN
         FACTOR = ALSQ
         DO 200 I=KRANK+1,ME
            CALL DCOPY (L, 0.D0, 0, W(I,1), MDW)
  200    CONTINUE
C
C        Determine the rank of the remaining equality constraint
C        equations by eliminating within the block of constrained
C        variables.  Remove any redundant constraints.
C
         RECALC = .TRUE.
         LB = MIN(L+ME-KRANK, N)
         DO 270 I=L+1,LB
            IR = KRANK + I - L
            LEND = N
            MEND = ME
            CALL DWNLT1 (I, LEND, ME, IR, MDW, RECALC, IMAX, HBAR, H,
     +                   SCALE, W)
C
C           Update col ss and find pivot col
C
            CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C           Perform column interchange
C           Eliminate elements in the I-th col.
C
            DO 240 J=ME,IR+1,-1
               IF (W(J,I).NE.0.D0) THEN
                 CALL DROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),
     +                        SPARAM)
                  W(J,I) = 0.D0
                  CALL DROTM (N+1-I, W(J-1,I+1), MDW,W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  240       CONTINUE
C
C           I=column being eliminated.
C           Test independence of incoming column.
C           Remove any redundant or dependent equality constraints.
C
            IF (.NOT.DWNLT2(ME, MEND, IR, FACTOR,TAU,SCALE,W(1,I))) THEN
               JJ = IR
               DO 260 IR=JJ,ME
                  CALL DCOPY (N, 0.D0, 0, W(IR,1), MDW)
                  RNORM = RNORM + (SCALE(IR)*W(IR,N+1)/ALSQ)*W(IR,N+1)
                  W(IR,N+1) = 0.D0
                  SCALE(IR) = 1.D0
C
C                 Reclassify the zeroed row as a least squares equation.
C
                  ITYPE(IR) = 1
  260          CONTINUE
C
C              Reduce ME to reflect any discovered dependent equality
C              constraints.
C
               ME = JJ - 1
               GO TO 280
            ENDIF
  270    CONTINUE
      ENDIF
C
C     Try to determine the variables KRANK+1 through L1 from the
C     least squares equations.  Continue the triangularization with
C     pivot element W(ME+1,I).
C
  280 IF (KRANK.LT.L1) THEN
         RECALC = .TRUE.
C
C        Set FACTOR=ALSQ to remove effect of heavy weight from
C        test for column independence.
C
         FACTOR = ALSQ
         DO 350 I=KRANK+1,L1
C
C           Set IR to point to the ME+1-st row.
C
            IR = ME+1
            LEND = L
            MEND = M
            CALL DWNLT1 (I, L, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,
     +                   W)
C
C           Update column SS and find pivot column.
C
            CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C           Perform column interchange.
C           Eliminate I-th column below the IR-th element.
C
            DO 320 J=M,IR+1,-1
               IF (W(J,I).NE.0.D0) THEN
                 CALL DROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),
     +                        SPARAM)
                  W(J,I) = 0.D0
                  CALL DROTM (N+1-I, W(J-1,I+1),  MDW, W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  320       CONTINUE
C
C           Test if new pivot element is near zero.
C           If so, the column is dependent.
C           Then check row norm test to be classified as independent.
C
            T = SCALE(IR)*W(IR,I)**2
            INDEP = T .GT. (TAU*EANORM)**2
            IF (INDEP) THEN
               RN = 0.D0
               DO 340 I1=IR,M
                  DO 330 J1=I+1,N
                     RN = MAX(RN, SCALE(I1)*W(I1,J1)**2)
  330             CONTINUE
  340          CONTINUE
               INDEP = T .GT. RN*TAU**2
            ENDIF
C
C           If independent, swap the IR-th and KRANK+1-th rows to
C           maintain the triangular form.  Update the rank indicator
C           KRANK and the equality constraint pointer ME.
C
            IF (.NOT.INDEP) GO TO 360
            CALL DSWAP(N+1, W(KRANK+1,1), MDW, W(IR,1), MDW)
            CALL DSWAP(1, SCALE(KRANK+1), 1, SCALE(IR), 1)
C
C           Reclassify the least square equation as an equality
C           constraint and rescale it.
C
            ITYPE(IR) = 0
            T = SQRT(SCALE(KRANK+1))
            CALL DSCAL(N+1, T, W(KRANK+1,1), MDW)
            SCALE(KRANK+1) = ALSQ
            ME = ME+1
            KRANK = KRANK+1
  350    CONTINUE
      ENDIF
C
C     If pseudorank is less than L, apply Householder transformation.
C     from right.
C
  360 IF (KRANK.LT.L) THEN
         DO 370 J=KRANK,1,-1
            CALL DH12 (1, J, KRANK+1, L, W(J,1), MDW, H(J), W, MDW, 1,
     +                J-1)
  370    CONTINUE
      ENDIF
C
      NIV = KRANK + NSOLN - L
      IF (L.EQ.N) DONE = .TRUE.
C
C     End of initial triangularization.
C
      IDOPE(1) = ME
      IDOPE(2) = KRANK
      IDOPE(3) = NIV
      RETURN
      END
*DECK DWNLSM
      SUBROUTINE DWNLSM (W, MDW, MME, MA, N, L, PRGOPT, X, RNORM, MODE,
     +   IPIVOT, ITYPE, WD, H, SCALE, Z, TEMP, D)
C***BEGIN PROLOGUE  DWNLSM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DWNNLS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLSM-S, DWNLSM-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to DWNNLS.
C     The documentation for DWNNLS has complete usage instructions.
C
C     In addition to the parameters discussed in the prologue to
C     subroutine DWNNLS, the following work arrays are used in
C     subroutine DWNLSM  (they are passed through the calling
C     sequence from DWNNLS for purposes of variable dimensioning).
C     Their contents will in general be of no interest to the user.
C
C     Variables of type REAL are DOUBLE PRECISION.
C
C         IPIVOT(*)
C            An array of length N.  Upon completion it contains the
C         pivoting information for the cols of W(*,*).
C
C         ITYPE(*)
C            An array of length M which is used to keep track
C         of the classification of the equations.  ITYPE(I)=0
C         denotes equation I as an equality constraint.
C         ITYPE(I)=1 denotes equation I as a least squares
C         equation.
C
C         WD(*)
C            An array of length N.  Upon completion it contains the
C         dual solution vector.
C
C         H(*)
C            An array of length N.  Upon completion it contains the
C         pivot scalars of the Householder transformations performed
C         in the case KRANK.LT.L.
C
C         SCALE(*)
C            An array of length M which is used by the subroutine
C         to store the diagonal matrix of weights.
C         These are used to apply the modified Givens
C         transformations.
C
C         Z(*),TEMP(*)
C            Working arrays of length N.
C
C         D(*)
C            An array of length N that contains the
C         column scaling for the matrix (E).
C                                       (A)
C
C***SEE ALSO  DWNNLS
C***ROUTINES CALLED  D1MACH, DASUM, DAXPY, DCOPY, DH12, DNRM2, DROTM,
C                    DROTMG, DSCAL, DSWAP, DWNLIT, IDAMAX, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Fixed an error message.  (RWC)
C   900604  DP version created from SP version.  (RWC)
C   900911  Restriction on value of ALAMDA included.  (WRB)
C***END PROLOGUE  DWNLSM
      INTEGER IPIVOT(*), ITYPE(*), L, MA, MDW, MME, MODE, N
      DOUBLE PRECISION D(*), H(*), PRGOPT(*), RNORM, SCALE(*), TEMP(*),
     *   W(MDW,*), WD(*), X(*), Z(*)
C
      EXTERNAL D1MACH, DASUM, DAXPY, DCOPY, DH12, DNRM2, DROTM, DROTMG,
     *   DSCAL, DSWAP, DWNLIT, IDAMAX, XERMSG
      DOUBLE PRECISION D1MACH, DASUM, DNRM2
      INTEGER IDAMAX
C
      DOUBLE PRECISION ALAMDA, ALPHA, ALSQ, AMAX, BLOWUP, BNORM,
     *   DOPE(3), DRELPR, EANORM, FAC, SM, SPARAM(5), T, TAU, WMAX, Z2,
     *   ZZ
      INTEGER I, IDOPE(3), IMAX, ISOL, ITEMP, ITER, ITMAX, IWMAX, J,
     *   JCON, JP, KEY, KRANK, L1, LAST, LINK, M, ME, NEXT, NIV, NLINK,
     *   NOPT, NSOLN, NTIMES
      LOGICAL DONE, FEASBL, FIRST, HITCON, POS
C
      SAVE DRELPR, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DWNLSM
C
C     Initialize variables.
C     DRELPR is the precision for the particular machine
C     being used.  This logic avoids resetting it every entry.
C
      IF (FIRST) DRELPR = D1MACH(4)
      FIRST = .FALSE.
C
C     Set the nominal tolerance used in the code.
C
      TAU = SQRT(DRELPR)
C
      M = MA + MME
      ME = MME
      MODE = 2
C
C     To process option vector
C
      FAC = 1.D-4
C
C     Set the nominal blow up factor used in the code.
C
      BLOWUP = TAU
C
C     The nominal column scaling used in the code is
C     the identity scaling.
C
      CALL DCOPY (N, 1.D0, 0, D, 1)
C
C     Define bound for number of options to change.
C
      NOPT = 1000
C
C     Define bound for positive value of LINK.
C
      NLINK = 100000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.LE.0 .OR. LINK.GT.NLINK) THEN
         CALL XERMSG ('SLATEC', 'DWNLSM',
     +      'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
         RETURN
      ENDIF
C
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
         CALL XERMSG ('SLATEC', 'DWNLSM',
     +      'IN DWNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.',
     +      3, 1)
            RETURN
         ENDIF
C
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.6 .AND. PRGOPT(LAST+2).NE.0.D0) THEN
            DO 110 J = 1,N
               T = DNRM2(M,W(1,J),1)
               IF (T.NE.0.D0) T = 1.D0/T
               D(J) = T
  110       CONTINUE
         ENDIF
C
         IF (KEY.EQ.7) CALL DCOPY (N, PRGOPT(LAST+2), 1, D, 1)
         IF (KEY.EQ.8) TAU = MAX(DRELPR,PRGOPT(LAST+2))
         IF (KEY.EQ.9) BLOWUP = MAX(DRELPR,PRGOPT(LAST+2))
C
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
            CALL XERMSG ('SLATEC', 'DWNLSM',
     +         'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
            RETURN
         ENDIF
C
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
C
      DO 120 J = 1,N
         CALL DSCAL (M, D(J), W(1,J), 1)
  120 CONTINUE
C
C     Process option vector
C
      DONE = .FALSE.
      ITER = 0
      ITMAX = 3*(N-L)
      MODE = 0
      NSOLN = L
      L1 = MIN(M,L)
C
C     Compute scale factor to apply to equality constraint equations.
C
      DO 130 J = 1,N
         WD(J) = DASUM(M,W(1,J),1)
  130 CONTINUE
C
      IMAX = IDAMAX(N,WD,1)
      EANORM = WD(IMAX)
      BNORM = DASUM(M,W(1,N+1),1)
      ALAMDA = EANORM/(DRELPR*FAC)
C
C     On machines, such as the VAXes using D floating, with a very
C     limited exponent range for double precision values, the previously
C     computed value of ALAMDA may cause an overflow condition.
C     Therefore, this code further limits the value of ALAMDA.
C
      ALAMDA = MIN(ALAMDA,SQRT(D1MACH(2)))
C
C     Define scaling diagonal matrix for modified Givens usage and
C     classify equation types.
C
      ALSQ = ALAMDA**2
      DO 140 I = 1,M
C
C        When equation I is heavily weighted ITYPE(I)=0,
C        else ITYPE(I)=1.
C
         IF (I.LE.ME) THEN
            T = ALSQ
            ITEMP = 0
         ELSE
            T = 1.D0
            ITEMP = 1
         ENDIF
         SCALE(I) = T
         ITYPE(I) = ITEMP
  140 CONTINUE
C
C     Set the solution vector X(*) to zero and the column interchange
C     matrix to the identity.
C
      CALL DCOPY (N, 0.D0, 0, X, 1)
      DO 150 I = 1,N
         IPIVOT(I) = I
  150 CONTINUE
C
C     Perform initial triangularization in the submatrix
C     corresponding to the unconstrained variables.
C     Set first L components of dual vector to zero because
C     these correspond to the unconstrained variables.
C
      CALL DCOPY (L, 0.D0, 0, WD, 1)
C
C     The arrays IDOPE(*) and DOPE(*) are used to pass
C     information to DWNLIT().  This was done to avoid
C     a long calling sequence or the use of COMMON.
C
      IDOPE(1) = ME
      IDOPE(2) = NSOLN
      IDOPE(3) = L1
C
      DOPE(1) = ALSQ
      DOPE(2) = EANORM
      DOPE(3) = TAU
      CALL DWNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM,
     +            IDOPE, DOPE, DONE)
      ME    = IDOPE(1)
      KRANK = IDOPE(2)
      NIV   = IDOPE(3)
C
C     Perform WNNLS algorithm using the following steps.
C
C     Until(DONE)
C        compute search direction and feasible point
C        when (HITCON) add constraints
C        else perform multiplier test and drop a constraint
C        fin
C     Compute-Final-Solution
C
C     To compute search direction and feasible point,
C     solve the triangular system of currently non-active
C     variables and store the solution in Z(*).
C
C     To solve system
C     Copy right hand side into TEMP vector to use overwriting method.
C
  160 IF (DONE) GO TO 330
      ISOL = L + 1
      IF (NSOLN.GE.ISOL) THEN
         CALL DCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 170 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
C
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.D0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL DAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  170    CONTINUE
      ENDIF
C
C     Increment iteration counter and check against maximum number
C     of iterations.
C
      ITER = ITER + 1
      IF (ITER.GT.ITMAX) THEN
         MODE = 1
         DONE = .TRUE.
      ENDIF
C
C     Check to see if any constraints have become active.
C     If so, calculate an interpolation factor so that all
C     active constraints are removed from the basis.
C
      ALPHA = 2.D0
      HITCON = .FALSE.
      DO 180 J = L+1,NSOLN
         ZZ = Z(J)
         IF (ZZ.LE.0.D0) THEN
            T = X(J)/(X(J)-ZZ)
            IF (T.LT.ALPHA) THEN
               ALPHA = T
               JCON = J
            ENDIF
            HITCON = .TRUE.
         ENDIF
  180 CONTINUE
C
C     Compute search direction and feasible point
C
      IF (HITCON) THEN
C
C        To add constraints, use computed ALPHA to interpolate between
C        last feasible solution X(*) and current unconstrained (and
C        infeasible) solution Z(*).
C
         DO 190 J = L+1,NSOLN
            X(J) = X(J) + ALPHA*(Z(J)-X(J))
  190    CONTINUE
         FEASBL = .FALSE.
C
C        Remove column JCON and shift columns JCON+1 through N to the
C        left.  Swap column JCON into the N th position.  This achieves
C        upper Hessenberg form for the nonactive constraints and
C        leaves an upper Hessenberg matrix to retriangularize.
C
  200    DO 210 I = 1,M
            T = W(I,JCON)
            CALL DCOPY (N-JCON, W(I, JCON+1), MDW, W(I, JCON), MDW)
            W(I,N) = T
  210    CONTINUE
C
C        Update permuted index vector to reflect this shift and swap.
C
         ITEMP = IPIVOT(JCON)
         DO 220 I = JCON,N - 1
            IPIVOT(I) = IPIVOT(I+1)
  220    CONTINUE
         IPIVOT(N) = ITEMP
C
C        Similarly permute X(*) vector.
C
         CALL DCOPY (N-JCON, X(JCON+1), 1, X(JCON), 1)
         X(N) = 0.D0
         NSOLN = NSOLN - 1
         NIV = NIV - 1
C
C        Retriangularize upper Hessenberg matrix after adding
C        constraints.
C
         I = KRANK + JCON - L
         DO 230 J = JCON,NSOLN
            IF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) THEN
C
C              Zero IP1 to I in column J
C
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL DROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL DROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) THEN
C
C              Zero IP1 to I in column J
C
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL DROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL DROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.0) THEN
               CALL DSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
               CALL DSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
               ITEMP = ITYPE(I+1)
               ITYPE(I+1) = ITYPE(I)
               ITYPE(I) = ITEMP
C
C              Swapped row was formerly a pivot element, so it will
C              be large enough to perform elimination.
C              Zero IP1 to I in column J.
C
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL DROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL DROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.1) THEN
               IF (SCALE(I)*W(I,J)**2/ALSQ.GT.(TAU*EANORM)**2) THEN
C
C                 Zero IP1 to I in column J
C
                  IF (W(I+1,J).NE.0.D0) THEN
                     CALL DROTMG (SCALE(I), SCALE(I+1), W(I,J),
     +                            W(I+1,J), SPARAM)
                     W(I+1,J) = 0.D0
                     CALL DROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                           SPARAM)
                  ENDIF
               ELSE
                  CALL DSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
                  CALL DSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
                  ITEMP = ITYPE(I+1)
                  ITYPE(I+1) = ITYPE(I)
                  ITYPE(I) = ITEMP
                  W(I+1,J) = 0.D0
               ENDIF
            ENDIF
            I = I + 1
  230    CONTINUE
C
C        See if the remaining coefficients in the solution set are
C        feasible.  They should be because of the way ALPHA was
C        determined.  If any are infeasible, it is due to roundoff
C        error.  Any that are non-positive will be set to zero and
C        removed from the solution set.
C
         DO 240 JCON = L+1,NSOLN
            IF (X(JCON).LE.0.D0) GO TO 250
  240    CONTINUE
         FEASBL = .TRUE.
  250    IF (.NOT.FEASBL) GO TO 200
      ELSE
C
C        To perform multiplier test and drop a constraint.
C
         CALL DCOPY (NSOLN, Z, 1, X, 1)
         IF (NSOLN.LT.N) CALL DCOPY (N-NSOLN, 0.D0, 0, X(NSOLN+1), 1)
C
C        Reclassify least squares equations as equalities as necessary.
C
         I = NIV + 1
  260    IF (I.LE.ME) THEN
            IF (ITYPE(I).EQ.0) THEN
               I = I + 1
            ELSE
               CALL DSWAP (N+1, W(I,1), MDW, W(ME,1), MDW)
               CALL DSWAP (1, SCALE(I), 1, SCALE(ME), 1)
               ITEMP = ITYPE(I)
               ITYPE(I) = ITYPE(ME)
               ITYPE(ME) = ITEMP
               ME = ME - 1
            ENDIF
            GO TO 260
         ENDIF
C
C        Form inner product vector WD(*) of dual coefficients.
C
         DO 280 J = NSOLN+1,N
            SM = 0.D0
            DO 270 I = NSOLN+1,M
               SM = SM + SCALE(I)*W(I,J)*W(I,N+1)
  270       CONTINUE
            WD(J) = SM
  280    CONTINUE
C
C        Find J such that WD(J)=WMAX is maximum.  This determines
C        that the incoming column J will reduce the residual vector
C        and be positive.
C
  290    WMAX = 0.D0
         IWMAX = NSOLN + 1
         DO 300 J = NSOLN+1,N
            IF (WD(J).GT.WMAX) THEN
               WMAX = WD(J)
               IWMAX = J
            ENDIF
  300    CONTINUE
         IF (WMAX.LE.0.D0) GO TO 330
C
C        Set dual coefficients to zero for incoming column.
C
         WD(IWMAX) = 0.D0
C
C        WMAX .GT. 0.D0, so okay to move column IWMAX to solution set.
C        Perform transformation to retriangularize, and test for near
C        linear dependence.
C
C        Swap column IWMAX into NSOLN-th position to maintain upper
C        Hessenberg form of adjacent columns, and add new column to
C        triangular decomposition.
C
         NSOLN = NSOLN + 1
         NIV = NIV + 1
         IF (NSOLN.NE.IWMAX) THEN
            CALL DSWAP (M, W(1,NSOLN), 1, W(1,IWMAX), 1)
            WD(IWMAX) = WD(NSOLN)
            WD(NSOLN) = 0.D0
            ITEMP = IPIVOT(NSOLN)
            IPIVOT(NSOLN) = IPIVOT(IWMAX)
            IPIVOT(IWMAX) = ITEMP
         ENDIF
C
C        Reduce column NSOLN so that the matrix of nonactive constraints
C        variables is triangular.
C
         DO 320 J = M,NIV+1,-1
            JP = J - 1
C
C           When operating near the ME line, test to see if the pivot
C           element is near zero.  If so, use the largest element above
C           it as the pivot.  This is to maintain the sharp interface
C           between weighted and non-weighted rows in all cases.
C
            IF (J.EQ.ME+1) THEN
               IMAX = ME
               AMAX = SCALE(ME)*W(ME,NSOLN)**2
               DO 310 JP = J - 1,NIV,-1
                  T = SCALE(JP)*W(JP,NSOLN)**2
                  IF (T.GT.AMAX) THEN
                     IMAX = JP
                     AMAX = T
                  ENDIF
  310          CONTINUE
               JP = IMAX
            ENDIF
C
            IF (W(J,NSOLN).NE.0.D0) THEN
               CALL DROTMG (SCALE(JP), SCALE(J), W(JP,NSOLN),
     +                      W(J,NSOLN), SPARAM)
               W(J,NSOLN) = 0.D0
               CALL DROTM (N+1-NSOLN, W(JP,NSOLN+1), MDW, W(J,NSOLN+1),
     +                     MDW, SPARAM)
            ENDIF
  320    CONTINUE
C
C        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
C        this is nonpositive or too large.  If this was true or if the
C        pivot term was zero, reject the column as dependent.
C
         IF (W(NIV,NSOLN).NE.0.D0) THEN
            ISOL = NIV
            Z2 = W(ISOL,N+1)/W(ISOL,NSOLN)
            Z(NSOLN) = Z2
            POS = Z2 .GT. 0.D0
            IF (Z2*EANORM.GE.BNORM .AND. POS) THEN
               POS = .NOT. (BLOWUP*Z2*EANORM.GE.BNORM)
            ENDIF
C
C           Try to add row ME+1 as an additional equality constraint.
C           Check size of proposed new solution component.
C           Reject it if it is too large.
C
         ELSEIF (NIV.LE.ME .AND. W(ME+1,NSOLN).NE.0.D0) THEN
            ISOL = ME + 1
            IF (POS) THEN
C
C              Swap rows ME+1 and NIV, and scale factors for these rows.
C
               CALL DSWAP (N+1, W(ME+1,1), MDW, W(NIV,1), MDW)
               CALL DSWAP (1, SCALE(ME+1), 1, SCALE(NIV), 1)
               ITEMP = ITYPE(ME+1)
               ITYPE(ME+1) = ITYPE(NIV)
               ITYPE(NIV) = ITEMP
               ME = ME + 1
            ENDIF
         ELSE
            POS = .FALSE.
         ENDIF
C
         IF (.NOT.POS) THEN
            NSOLN = NSOLN - 1
            NIV = NIV - 1
         ENDIF
         IF (.NOT.(POS.OR.DONE)) GO TO 290
      ENDIF
      GO TO 160
C
C     Else perform multiplier test and drop a constraint.  To compute
C     final solution.  Solve system, store results in X(*).
C
C     Copy right hand side into TEMP vector to use overwriting method.
C
  330 ISOL = 1
      IF (NSOLN.GE.ISOL) THEN
         CALL DCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 340 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
C
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.D0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL DAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  340    CONTINUE
      ENDIF
C
C     Solve system.
C
      CALL DCOPY (NSOLN, Z, 1, X, 1)
C
C     Apply Householder transformations to X(*) if KRANK.LT.L
C
      IF (KRANK.LT.L) THEN
         DO 350 I = 1,KRANK
            CALL DH12 (2, I, KRANK+1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
  350    CONTINUE
      ENDIF
C
C     Fill in trailing zeroes for constrained variables not in solution.
C
      IF (NSOLN.LT.N) CALL DCOPY (N-NSOLN, 0.D0, 0, X(NSOLN+1), 1)
C
C     Permute solution vector to natural order.
C
      DO 380 I = 1,N
         J = I
  360    IF (IPIVOT(J).EQ.I) GO TO 370
         J = J + 1
         GO TO 360
C
  370    IPIVOT(J) = IPIVOT(I)
         IPIVOT(I) = J
         CALL DSWAP (1, X(J), 1, X(I), 1)
  380 CONTINUE
C
C     Rescale the solution using the column scaling.
C
      DO 390 J = 1,N
         X(J) = X(J)*D(J)
  390 CONTINUE
C
      DO 400 I = NSOLN+1,M
         T = W(I,N+1)
         IF (I.LE.ME) T = T/ALAMDA
         T = (SCALE(I)*T)*T
         RNORM = RNORM + T
  400 CONTINUE
C
      RNORM = SQRT(RNORM)
      RETURN
      END
*DECK DWNLT1
      SUBROUTINE DWNLT1 (I, LEND, MEND, IR, MDW, RECALC, IMAX, HBAR, H,
     +   SCALE, W)
C***BEGIN PROLOGUE  DWNLT1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLT1-S, DWNLT1-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     To update the column Sum Of Squares and find the pivot column.
C     The column Sum of Squares Vector will be updated at each step.
C     When numerically necessary, these values will be recomputed.
C
C***SEE ALSO  DWNLIT
C***ROUTINES CALLED  IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
C   900604  DP version created from SP version.  (RWC)
C***END PROLOGUE  DWNLT1
      INTEGER I, IMAX, IR, LEND, MDW, MEND
      DOUBLE PRECISION H(*), HBAR, SCALE(*), W(MDW,*)
      LOGICAL RECALC
C
      EXTERNAL IDAMAX
      INTEGER IDAMAX
C
      INTEGER J, K
C
C***FIRST EXECUTABLE STATEMENT  DWNLT1
      IF (IR.NE.1 .AND. (.NOT.RECALC)) THEN
C
C        Update column SS=sum of squares.
C
         DO 10 J=I,LEND
            H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
   10    CONTINUE
C
C        Test for numerical accuracy.
C
         IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
         RECALC = (HBAR+1.E-3*H(IMAX)) .EQ. HBAR
      ENDIF
C
C     If required, recalculate column SS, using rows IR through MEND.
C
      IF (RECALC) THEN
         DO 30 J=I,LEND
            H(J) = 0.D0
            DO 20 K=IR,MEND
               H(J) = H(J) + SCALE(K)*W(K,J)**2
   20       CONTINUE
   30    CONTINUE
C
C        Find column with largest SS.
C
         IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
         HBAR = H(IMAX)
      ENDIF
      RETURN
      END
*DECK DWNLT2
      LOGICAL FUNCTION DWNLT2 (ME, MEND, IR, FACTOR, TAU, SCALE, WIC)
C***BEGIN PROLOGUE  DWNLT2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLT2-S, DWNLT2-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     To test independence of incoming column.
C
C     Test the column IC to determine if it is linearly independent
C     of the columns already in the basis.  In the initial tri. step,
C     we usually want the heavy weight ALAMDA to be included in the
C     test for independence.  In this case, the value of FACTOR will
C     have been set to 1.E0 before this procedure is invoked.
C     In the potentially rank deficient problem, the value of FACTOR
C     will have been set to ALSQ=ALAMDA**2 to remove the effect of the
C     heavy weight from the test for independence.
C
C     Write new column as partitioned vector
C           (A1)  number of components in solution so far = NIV
C           (A2)  M-NIV components
C     And compute  SN = inverse weighted length of A1
C                  RN = inverse weighted length of A2
C     Call the column independent when RN .GT. TAU*SN
C
C***SEE ALSO  DWNLIT
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
C   900604  DP version created from SP version.  (RWC)
C***END PROLOGUE  DWNLT2
      DOUBLE PRECISION FACTOR, SCALE(*), TAU, WIC(*)
      INTEGER IR, ME, MEND
C
      DOUBLE PRECISION RN, SN, T
      INTEGER J
C
C***FIRST EXECUTABLE STATEMENT  DWNLT2
      SN = 0.E0
      RN = 0.E0
      DO 10 J=1,MEND
         T = SCALE(J)
         IF (J.LE.ME) T = T/FACTOR
         T = T*WIC(J)**2
C
         IF (J.LT.IR) THEN
            SN = SN + T
         ELSE
            RN = RN + T
         ENDIF
   10 CONTINUE
      DWNLT2 = RN .GT. SN*TAU**2
      RETURN
      END
*DECK DWNLT3
      SUBROUTINE DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C***BEGIN PROLOGUE  DWNLT3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Perform column interchange.
C     Exchange elements of permuted index vector and perform column
C     interchanges.
C
C***SEE ALSO  DWNLIT
C***ROUTINES CALLED  DSWAP
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
C   900604  DP version created from SP version.  (RWC)
C***END PROLOGUE  DWNLT3
      INTEGER I, IMAX, IPIVOT(*), M, MDW
      DOUBLE PRECISION H(*), W(MDW,*)
C
      EXTERNAL DSWAP
C
      DOUBLE PRECISION T
      INTEGER ITEMP
C
C***FIRST EXECUTABLE STATEMENT  DWNLT3
      IF (IMAX.NE.I) THEN
         ITEMP        = IPIVOT(I)
         IPIVOT(I)    = IPIVOT(IMAX)
         IPIVOT(IMAX) = ITEMP
C
         CALL DSWAP(M, W(1,IMAX), 1, W(1,I), 1)
C
         T       = H(IMAX)
         H(IMAX) = H(I)
         H(I)    = T
      ENDIF
      RETURN
      END
*DECK DWNNLS
      SUBROUTINE DWNNLS (W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE,
     +   IWORK, WORK)
C***BEGIN PROLOGUE  DWNNLS
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality constraints and nonnegativity constraints on
C            selected variables.
C***LIBRARY   SLATEC
C***CATEGORY  K1A2A
C***TYPE      DOUBLE PRECISION (WNNLS-S, DWNNLS-D)
C***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
C             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
C             NONNEGATIVITY CONSTRAINTS, QUADRATIC PROGRAMMING
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     This subprogram solves a linearly constrained least squares
C     problem.  Suppose there are given matrices E and A of
C     respective dimensions ME by N and MA by N, and vectors F
C     and B of respective lengths ME and MA.  This subroutine
C     solves the problem
C
C               EX = F, (equations to be exactly satisfied)
C
C               AX = B, (equations to be approximately satisfied,
C                        in the least squares sense)
C
C               subject to components L+1,...,N nonnegative
C
C     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted.
C
C     The problem is reposed as problem DWNNLS
C
C               (WT*E)X = (WT*F)
C               (   A)    (   B), (least squares)
C               subject to components L+1,...,N nonnegative.
C
C     The subprogram chooses the heavy weight (or penalty parameter) WT.
C
C     The parameters for DWNNLS are
C
C     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     W(*,*),MDW,  The array W(*,*) is double subscripted with first
C     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
C                  discussion let us call M = ME + MA.  Then MDW
C                  must satisfy MDW.GE.M.  The condition MDW.LT.M
C                  is an error.
C
C                  The array W(*,*) contains the matrices and vectors
C
C                       (E  F)
C                       (A  B)
C
C                  in rows and columns 1,...,M and 1,...,N+1
C                  respectively.  Columns 1,...,L correspond to
C                  unconstrained variables X(1),...,X(L).  The
C                  remaining variables are constrained to be
C                  nonnegative. The condition L.LT.0 or L.GT.N is
C                  an error.
C
C     PRGOPT(*)    This double precision array is the option vector.
C                  If the user is satisfied with the nominal
C                  subprogram features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case LINK=1 and the values KEY and DATA SET
C                  are not referenced. The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1)=LINK1 (link to first entry of next group)
C               .  PRGOPT(2)=KEY1 (key to the option change)
C               .  PRGOPT(3)=DATA VALUE (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
C               .  PRGOPT(LINK1+2)=DATA VALUE
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK)=1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK.GT.NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the program limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000 an error
C                  message is printed and the subprogram returns.
C
C                  OPTIONS..
C
C                  KEY=6
C                         Scale the nonzero columns of the
C                  entire data matrix
C                  (E)
C                  (A)
C                  to have length one. The DATA SET for
C                  this option is a single value.  It must
C                  be nonzero if unit length column scaling is
C                  desired.
C
C                  KEY=7
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  with a user-provided diagonal matrix.
C                  The DATA SET for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=8
C                         Change the rank determination tolerance from
C                  the nominal value of SQRT(SRELPR).  This quantity
C                  can be no smaller than SRELPR, The arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The DATA SET for this option
C                  is the new tolerance.
C
C                  KEY=9
C                         Change the blow-up parameter from the
C                  nominal value of SQRT(SRELPR).  The reciprocal of
C                  this parameter is used in rejecting solution
C                  components as too large when a variable is
C                  first brought into the active set.  Too large
C                  means that the proposed component times the
C                  reciprocal of the parameter is not less than
C                  the ratio of the norms of the right-side
C                  vector and the data matrix.
C                  This parameter can be no smaller than SRELPR,
C                  the arithmetic-storage precision.
C
C                  For example, suppose we want to provide
C                  a diagonal matrix to scale the problem
C                  matrix and change the tolerance used for
C                  determining linear dependence of dropped col
C                  vectors.  For these options the dimensions of
C                  PRGOPT(*) must be at least N+6.  The FORTRAN
C                  statements defining these options would
C                  be as follows.
C
C                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
C                  PRGOPT(2)=7 (user-provided scaling key)
C
C                  CALL DCOPY(N,D,1,PRGOPT(3),1) (copy the N
C                  scaling factors from a user array called D(*)
C                  into PRGOPT(3)-PRGOPT(N+2))
C
C                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
C                  PRGOPT(N+4)=8 (linear dependence tolerance key)
C                  PRGOPT(N+5)=... (new value of the tolerance)
C
C                  PRGOPT(N+6)=1 (no more options to change)
C
C
C     IWORK(1),    The amounts of working storage actually allocated
C     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
C                  respectively.  These quantities are compared with
C                  the actual amounts of storage needed for DWNNLS( ).
C                  Insufficient storage allocated for either WORK(*)
C                  or IWORK(*) is considered an error.  This feature
C                  was included in DWNNLS( ) because miscalculating
C                  the storage formulas for WORK(*) and IWORK(*)
C                  might very well lead to subtle and hard-to-find
C                  execution errors.
C
C                  The length of WORK(*) must be at least
C
C                  LW = ME+MA+5*N
C                  This test will not be made if IWORK(1).LE.0.
C
C                  The length of IWORK(*) must be at least
C
C                  LIW = ME+MA+N
C                  This test will not be made if IWORK(2).LE.0.
C
C     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     X(*)         An array dimensioned at least N, which will
C                  contain the N components of the solution vector
C                  on output.
C
C     RNORM        The residual norm of the solution.  The value of
C                  RNORM contains the residual vector length of the
C                  equality constraints and least squares equations.
C
C     MODE         The value of MODE indicates the success or failure
C                  of the subprogram.
C
C                  MODE = 0  Subprogram completed successfully.
C
C                       = 1  Max. number of iterations (equal to
C                            3*(N-L)) exceeded. Nearly all problems
C                            should complete in fewer than this
C                            number of iterations. An approximate
C                            solution and its corresponding residual
C                            vector length are in X(*) and RNORM.
C
C                       = 2  Usage error occurred.  The offending
C                            condition is noted with the error
C                            processing subprogram, XERMSG( ).
C
C     User-designated
C     Working arrays..
C
C     WORK(*)      A double precision working array of length at least
C                  M + 5*N.
C
C     IWORK(*)     An integer-valued working array of length at least
C                  M+N.
C
C***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Report SAND77-0552, Sandia
C                 Laboratories, June 1978.
C               K. H. Haskell and R. J. Hanson, Selected algorithms for
C                 the linearly constrained least squares problem - a
C                 users guide, Report SAND78-1290, Sandia Laboratories,
C                 August 1979.
C               K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Mathematical Programming
C                 21 (1981), pp. 98-118.
C               R. J. Hanson and K. H. Haskell, Two algorithms for the
C                 linearly constrained least squares problem, ACM
C                 Transactions on Mathematical Software, September 1982.
C               C. L. Lawson and R. J. Hanson, Solving Least Squares
C                 Problems, Prentice-Hall, Inc., 1974.
C***ROUTINES CALLED  DWNLSM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls, change Prologue
C           comments to agree with WNNLS.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DWNNLS
      INTEGER IWORK(*), L, L1, L2, L3, L4, L5, LIW, LW, MA, MDW, ME,
     *     MODE, N
      DOUBLE PRECISION  PRGOPT(*), RNORM, W(MDW,*), WORK(*), X(*)
      CHARACTER*8 XERN1
C***FIRST EXECUTABLE STATEMENT  DWNNLS
      MODE = 0
      IF (MA+ME.LE.0 .OR. N.LE.0) RETURN
C
      IF (IWORK(1).GT.0) THEN
         LW = ME + MA + 5*N
         IF (IWORK(1).LT.LW) THEN
            WRITE (XERN1, '(I8)') LW
            CALL XERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR WORK(*), NEED LW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
C
      IF (IWORK(2).GT.0) THEN
         LIW = ME + MA + N
         IF (IWORK(2).LT.LIW) THEN
            WRITE (XERN1, '(I8)') LIW
            CALL XERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR IWORK(*), NEED LIW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
C
      IF (MDW.LT.ME+MA) THEN
         CALL XERMSG ('SLATEC', 'DWNNLS',
     *      'THE VALUE MDW.LT.ME+MA IS AN ERROR', 1, 1)
         MODE = 2
         RETURN
      ENDIF
C
      IF (L.LT.0 .OR. L.GT.N) THEN
         CALL XERMSG ('SLATEC', 'DWNNLS',
     *      'L.GE.0 .AND. L.LE.N IS REQUIRED', 2, 1)
         MODE = 2
         RETURN
      ENDIF
C
C     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
C     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
C     REQUIRED BY THE MAIN SUBROUTINE DWNLSM( ).
C
      L1 = N + 1
      L2 = L1 + N
      L3 = L2 + ME + MA
      L4 = L3 + N
      L5 = L4 + N
C
      CALL DWNLSM(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK,
     *            IWORK(L1), WORK(1), WORK(L1), WORK(L2), WORK(L3),
     *            WORK(L4), WORK(L5))
      RETURN
      END
*DECK FDUMP
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***PURPOSE  Symbolic dump (should be locally written).
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (FDUMP-A)
C***KEYWORDS  ERROR, XERMSG
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
*DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
C***BEGIN PROLOGUE  I1MACH
C***PURPOSE  Return integer machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      INTEGER (I1MACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   I1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument and can be referenced as follows:
C
C        K = I1MACH(I)
C
C   where I=1,...,16.  The (output) value of K above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   I/O unit numbers:
C     I1MACH( 1) = the standard input unit.
C     I1MACH( 2) = the standard output unit.
C     I1MACH( 3) = the standard punch unit.
C     I1MACH( 4) = the standard error message unit.
C
C   Words:
C     I1MACH( 5) = the number of bits per integer storage unit.
C     I1MACH( 6) = the number of characters per integer storage unit.
C
C   Integers:
C     assume integers are represented in the S-digit, base-A form
C
C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C     I1MACH( 7) = A, the base.
C     I1MACH( 8) = S, the number of base-A digits.
C     I1MACH( 9) = A**S - 1, the largest magnitude.
C
C   Floating-Point Numbers:
C     Assume floating-point numbers are represented in the T-digit,
C     base-B form
C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C                where 0 .LE. X(I) .LT. B for I=1,...,T,
C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C     I1MACH(10) = B, the base.
C
C   Single-Precision:
C     I1MACH(11) = T, the number of base-B digits.
C     I1MACH(12) = EMIN, the smallest exponent E.
C     I1MACH(13) = EMAX, the largest exponent E.
C
C   Double-Precision:
C     I1MACH(14) = T, the number of base-B digits.
C     I1MACH(15) = EMIN, the smallest exponent E.
C     I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   891012  Added VAX G-floating constants.  (WRB)
C   891012  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
C           (RWC)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added Convex -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
C           options.  (DWL, RWC and WRB).
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        129 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /          7 /
C     DATA IMACH( 2) /          2 /
C     DATA IMACH( 3) /          2 /
C     DATA IMACH( 4) /          2 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -256 /
C     DATA IMACH(13) /        255 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /       -256 /
C     DATA IMACH(16) /        255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /        -50 /
C     DATA IMACH(16) /         76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /     -32754 /
C     DATA IMACH(16) /      32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -4095 /
C     DATA IMACH(13) /       4094 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -4095 /
C     DATA IMACH(16) /       4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /    6LOUTPUT/
C     DATA IMACH( 5) /         60 /
C     DATA IMACH( 6) /         10 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /       -929 /
C     DATA IMACH(13) /       1070 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /       -929 /
C     DATA IMACH(16) /       1069 /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / Z'7FFFFFFF' /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16383 /
C     DATA IMACH(16) /      16383 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -pd8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 46 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         46 /
C     DATA IMACH( 9) / 1777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 64 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 777777777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /         11 /
C     DATA IMACH( 2) /         12 /
C     DATA IMACH( 3) /          8 /
C     DATA IMACH( 4) /         10 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         24 /
C     DATA IMACH( 6) /          3 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         23 /
C     DATA IMACH( 9) /    8388607 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         38 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /         43 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         63 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         39 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         55 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          7 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1015 /
C     DATA IMACH(16) /       1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) /  Z7FFFFFFF /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          0 /
      DATA IMACH( 4) /          0 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -125 /
      DATA IMACH(13) /        127 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1021 /
      DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         54 /
C     DATA IMACH(15) /       -101 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         62 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1021 /
C     DATA IMACH(13) /       1024 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16381 /
C     DATA IMACH(16) /      16384 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          1 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /      -1024 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /          1 /
C     DATA IMACH( 2) /          1 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
      STOP
      END
*DECK J4SAVE
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
C***BEGIN PROLOGUE  J4SAVE
C***SUBSIDIARY
C***PURPOSE  Save or recall global variables needed by error
C            handling routines.
C***LIBRARY   SLATEC (XERROR)
C***TYPE      INTEGER (J4SAVE-I)
C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                = 6 Refers to the 2nd unit for error messages
C                = 7 Refers to the 3rd unit for error messages
C                = 8 Refers to the 4th unit for error messages
C                = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C***SEE ALSO  XERMSG
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900205  Minor modifications to prologue.  (WRB)
C   900402  Added TYPE section.  (WRB)
C   910411  Added KEYWORDS section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
*DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
C***BEGIN PROLOGUE  XERCNT
C***SUBSIDIARY
C***PURPOSE  Allow user control over handling of errors.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCNT-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCNT.
C        If the user has provided his own version of XERCNT, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        LIBRAR - the library that the routine is in.
C        SUBROU - the subroutine that XERMSG is being called from
C        MESSG  - the first 20 characters of the error message.
C        NERR   - same as in the call to XERMSG.
C        LEVEL  - same as in the call to XERMSG.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
C           names, changed routine name from XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
*DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
C***BEGIN PROLOGUE  XERHLT
C***SUBSIDIARY
C***PURPOSE  Abort program execution and print error message.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERHLT-A)
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERHLT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG is as in XERMSG.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to delete length of character
C           and changed routine name from XERABT to XERHLT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END
*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
*DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
C***BEGIN PROLOGUE  XERPRN
C***SUBSIDIARY
C***PURPOSE  Print error messages processed by XERMSG.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERPRN-A)
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C This routine sends one or more lines to each of the (up to five)
C logical units to which error messages are to be sent.  This routine
C is called several times by XERMSG, sometimes with a single line to
C print and sometimes with a (potentially very long) message that may
C wrap around into multiple lines.
C
C PREFIX  Input argument of type CHARACTER.  This argument contains
C         characters to be put at the beginning of each line before
C         the body of the message.  No more than 16 characters of
C         PREFIX will be used.
C
C NPREF   Input argument of type INTEGER.  This argument is the number
C         of characters to use from PREFIX.  If it is negative, the
C         intrinsic function LEN is used to determine its length.  If
C         it is zero, PREFIX is not used.  If it exceeds 16 or if
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
C         used.  If NPREF is positive and the length of PREFIX is less
C         than NPREF, a copy of PREFIX extended with blanks to length
C         NPREF will be used.
C
C MESSG   Input argument of type CHARACTER.  This is the text of a
C         message to be printed.  If it is a long message, it will be
C         broken into pieces for printing on multiple lines.  Each line
C         will start with the appropriate prefix and be followed by a
C         piece of the message.  NWRAP is the number of characters per
C         piece; that is, after each NWRAP characters, we break and
C         start a new line.  In addition the characters '$$' embedded
C         in MESSG are a sentinel for a new line.  The counting of
C         characters up to NWRAP starts over for each new line.  The
C         value of NWRAP typically used by XERMSG is 72 since many
C         older error messages in the SLATEC Library are laid out to
C         rely on wrap-around every 72 characters.
C
C NWRAP   Input argument of type INTEGER.  This gives the maximum size
C         piece into which to break MESSG for printing on multiple
C         lines.  An embedded '$$' ends a line, and the count restarts
C         at the following character.  If a line break does not occur
C         on a blank (it would split a word) that word is moved to the
C         next line.  Values of NWRAP less than 16 will be treated as
C         16.  Values of NWRAP greater than 132 will be treated as 132.
C         The actual line length will be NPREF + NWRAP after NPREF has
C         been adjusted to fall between 0 and 16 and NWRAP has been
C         adjusted to fall between 16 and 132.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   880621  DATE WRITTEN
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
C           SLASH CHARACTER IN FORMAT STATEMENTS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
C           LINES TO BE PRINTED.
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Added code to break messages between words.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
*DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
C***BEGIN PROLOGUE  XERSVE
C***SUBSIDIARY
C***PURPOSE  Record that an error has occurred.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (XERSVE-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
C
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
C
C *Arguments:
C
C        LIBRAR :IN    is the library that the message is from.
C        SUBROU :IN    is the subroutine that the message is from.
C        MESSG  :IN    is the message to be saved.
C        KFLAG  :IN    indicates the action to be performed.
C                      when KFLAG > 0, the message in MESSG is saved.
C                      when KFLAG=0 the tables will be dumped and
C                      cleared.
C                      when KFLAG < 0, the tables will be dumped and
C                      not cleared.
C        NERR   :IN    is the error number.
C        LEVEL  :IN    is the error severity.
C        ICOUNT :OUT   the number of times this message has been seen,
C                      or zero if the table has overflowed and does not
C                      contain this message specifically.  When KFLAG=0,
C                      ICOUNT will not be altered.
C
C *Description:
C
C   Record that this error occurred and possibly dump and clear the
C   tables.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   800319  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900413  Routine modified to remove reference to KFLAG.  (WRB)
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
C           sequence, use IF-THEN-ELSE, make number of saved entries
C           easily changeable, changed routine name from XERSAV to
C           XERSVE.  (RWC)
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
*DECK XGETUA
      SUBROUTINE XGETUA (IUNITA, N)
C***BEGIN PROLOGUE  XGETUA
C***PURPOSE  Return unit number(s) to which error messages are being
C            sent.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XGETUA-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
#Need following modules commands:
#module unload DE_intel
#module load acml/5.3.1/ifort64 intel/14.1

FC=ifort -L/usr/nic/compiler/intel/14.1/lib/intel64 -lirc

SUBS=fdump j4save xercnt xerhlt xermsg xerprn xersve xgetua \
 d1mach dh12 dwnlit dwnlt1 dwnlt2 dwnlt3 i1mach

#Need acml for blas and lapack.
call_wnnls: call_wnnls.o dwnnls.o dwnlsm.o $(addsuffix .o, $(SUBS))
	$(FC) -o $@ $^ -lacml

.f90.o:
	$(FC) $(FFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<
f90:
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@

.SUFFIXES:
.SUFFIXES: .o .f .f90
6 6
0.9969 0.0005 0.6951 0.0001 0.9 0.0
0.0008 0.9976 0.0 0.0 0.0008 0.0000227
0.0 0.0 0.1875 0.0 0.0 0.0
0.0 3e-4 0.105 .5835 0.0 .000009
.00021 .0002 0.0 0.0 .0428 .0000024
0.0 0.0 0.0 0.0 0.0 .99966
50 0 20 30 0 0

6 6
0.9969 0.0005 0.6951 0.0001 0.9 0.0
0.0008 0.9976 0.0 0.0 0.0008 0.0000227
0.0 0.0 0.1875 0.0 0.0 0.0
0.0 3e-4 0.105 .5835 0.0 .000009
.00021 .0002 0.0 0.0 .0428 .0000024
0.0 0.0 0.0 0.0 0.0 .99966
50 0 20 30 0 0

