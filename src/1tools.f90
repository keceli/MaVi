!Last saved by Murat Keceli on 4:08:55 AM Dec 20, 2010 .
!ifort -O0 -g3 -c vib.f90 vci.f90 super.f90 main.f90 qff.f90
!ifort -O0 -g3  vib.o vci.o super.o main.o qff.o -L/scr/haku_2/shiozaki/intel/mkl/10.0.1.014/lib/em64t/ -lmkl -lmkl_em64t -lmkl_lapack -lguide -lpthread -Vaxlib -o vci
!gfortran -c -g ../main.f95
!gfortran -g main.o -L"C:\cygwin\lib\lapack" -lblas -llapack
MODULE constants
  IMPLICIT NONE
  REAL*8,PARAMETER:: c_pi=3.14159265358979323846d0
  !REAL*8,PARAMETER:: c_pi=4.d0*datan(1.d0)
  ! REAL*8,PARAMETER:: c_try=1.23456789123456789123456d0
  REAL*8,PARAMETER:: c_na = 6.02214d+23 ! Avogadro's number
  REAL*8,PARAMETER:: c_h = 6.62606896d-34 ! Planck's constant (J*s)
  REAL*8,PARAMETER:: c_c = 2.99792458d+10 ! Speed of light in vacuum (cm/s)
  REAL*8,PARAMETER:: c_h2wn = 2.194746d+05 ! Hartree to cm^-1
 !  REAL*8,PARAMETER:: c_h2wn = 2.1947463137049999D+05 ! Hartree to cm^-1 from Midas
  REAL*8,PARAMETER:: c_wn2h = 4.55633590401804997D-6 !1.d0/2.194746d+05 ! cm^-1 to hartree
  REAL*8,PARAMETER:: c_amu = 1.660538782d-27 ! = 1/_na/1000
  REAL*8,PARAMETER:: c_gamma = 2.96601703232432512D-2!4.d0 * c_pi * c_pi * c_c * 1.0d-23 / c_h / c_na;
  REAL*8,PARAMETER:: c_bohr=5.291772108d-11!bohr radius in meters =
  REAL*8,PARAMETER:: c_angs2bohr=1.8897261249935899d0!1d-10/c_bohr
  REAL*8,PARAMETER:: c_emass=9.10938215d-31!electron mass in kg
  REAL*8,PARAMETER:: c_amu2emass=1822.8884842645450d0!c_amu/c_emass
END MODULE constants

MODULE modtimer
  IMPLICIT NONE
CONTAINS
  !from Yu-ya Ohnishi
  SUBROUTINE wall_and_cpu_time(wall,cpu)
    IMPLICIT NONE
    REAL(8), INTENT(out) :: wall, cpu
    INTEGER(8):: t, t_rate, t_max

    CALL CPU_TIME(cpu)
    CALL SYSTEM_CLOCK(t,t_rate,t_max)
    wall = DBLE(t)/DBLE(t_rate)

  END SUBROUTINE wall_and_cpu_time
END MODULE modtimer

MODULE lapack
  IMPLICIT NONE
CONTAINS

  SUBROUTINE diag(order,matrix,eigenvalues)
    IMPLICIT NONE
    INTEGER :: order,info,lda,lwork
    REAL(8), DIMENSION(order,order)::matrix
    REAL(8), DIMENSION(order)::eigenvalues
    REAL(8), DIMENSION(3*order-1) :: work
    INTEGER          LWMAX
    PARAMETER        ( LWMAX = 1000 )
!     .. External Subroutines ..
    EXTERNAL         DSYEV
    INTRINSIC        INT, MIN
    lda=max(1,order)
    lwork = -1 ! to get the optimum lwork
    CALL dsyev('V','U',order,matrix,lda,eigenvalues,work,lwork,info)
    lwork = MIN( LWMAX, INT( WORK( 1 ) ) )
    CALL dsyev('V','U',order,matrix,lda,eigenvalues,work,lwork,info)
    !CALL diagnow(order,lwork,matrix,eigenvalues)
    !              PRINT*,'info=',info
  END SUBROUTINE diag

  SUBROUTINE diagnow(order,lwork,matrix,eigenvalues)
    IMPLICIT NONE
    INTEGER :: order,info,lda,lwork
    REAL(8), DIMENSION(order,order)::matrix
    REAL(8), DIMENSION(order)::eigenvalues
    REAL(8), DIMENSION(lwork) :: work
    lda=max(1,order)
    CALL dsyev('V','U',order,matrix,lda,eigenvalues,work,lwork,info)
    !              PRINT*,'info=',info
  END SUBROUTINE diagnow

!  SUBROUTINE complex_diag(order,matrix,eigenvalues)
!    IMPLICIT NONE
!    INTEGER :: order,info,lda,Lwork
!    COMPLEX*16, DIMENSION(order,order)::matrix
!    REAL(8), DIMENSION(order)::eigenvalues
!    COMPLEX*16, DIMENSION(2*order-1) :: work
!    REAL(8), DIMENSION(3*order-2) :: Rwork
!    lda=3*order-1
!    Lwork=2*order-1
!    CALL zheev('V','U',order,matrix,order,eigenvalues,work,Lwork,Rwork,info)
!  END SUBROUTINE complex_diag

  SUBROUTINE inverse(order,matrix)
    IMPLICIT NONE
    INTEGER :: order,info,lda
    REAL(8), DIMENSION(order,order)::matrix
    lda=order
    CALL DPOTRF( 'U', order, matrix, LDA, INFO )
    CALL DPOTRI( 'U', order, matrix, LDA, INFO )
    !              PRINT*,'info=',info
  END SUBROUTINE inverse

END MODULE lapack

MODULE yazar
  IMPLICIT NONE
CONTAINS
  SUBROUTINE writevec(size,vec)
    IMPLICIT NONE
    INTEGER i,size
    REAL*8,DIMENSION(size)::vec

    DO i=1,size
       WRITE(*,'(10F20.9)')vec(i)
    ENDDO

  END SUBROUTINE writevec

  SUBROUTINE writesquare(order,matrix)
    IMPLICIT NONE
    INTEGER i,j,order
    REAL*8,DIMENSION(order,order)::matrix

    DO i=1,order
       WRITE(*,'(10D20.9,6x)',ADVANCE='no')(matrix(i,j),j=1,order)
       WRITE(*,*)
    ENDDO

  END SUBROUTINE writesquare

  SUBROUTINE writemat(coef,nrow,ncol)
    IMPLICIT NONE
    INTEGER,INTENT(in)::nrow,ncol
    REAL*8,DIMENSION(nrow,ncol)::coef
    INTEGER::i,j

    DO i=1,nrow
       WRITE(*,'(10F15.6,6x)',ADVANCE='no')(coef(i,j),j=1,ncol)
       WRITE(*,*)
    ENDDO

  END SUBROUTINE writemat

END MODULE yazar

