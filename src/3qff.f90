MODULE modQFF
  IMPLICIT NONE
  CHARACTER::fcfile*30 !Force constant file name
  INTEGER::qfftype !1:001.hs Sindo 2: model calc in au 3:Midas scaled FC 4: chain
  REAL*8,ALLOCATABLE                              :: x_grad(:)
  COMPLEX*16,ALLOCATABLE                        :: q_grad(:,:)
  REAL*8,ALLOCATABLE                          :: x_harm(:,:,:)
  COMPLEX*16,ALLOCATABLE                    :: q_harm(:,:,:,:)
  !  real*8,allocatable                      :: x_cube(:,:,:,:,:)
  !  complex*16,allocatable               :: q_cube(:,:,:,:,:,:)
  !  real*8,allocatable                 :: x_quar(:,:,:,:,:,:,:)
  !  complex*16,allocatable :: q_quar(:,:,:,:,:,:,:,:)

  !variables for nMR QFF, n=1,2,3,4
  CHARACTER*100 :: title1MR,title2MR,title3MR,title4MR
  !REAL*8:: Eref
  REAL*8,ALLOCATABLE :: Gi(:),Hii(:),Ciii(:),Qiiii(:)
  REAL*8,ALLOCATABLE :: Hij(:,:),Ciij(:,:),Qiijj(:,:),Qiiij(:,:)
  REAL*8,ALLOCATABLE :: Cijk(:,:,:),Qiijk(:,:,:)
  REAL*8,ALLOCATABLE :: Qijkl(:,:,:,:)
  real*8::kcubic,kquartic
  LOGICAL::quartic,harmonic,allzero,fullquartic,nocubic,noCiij,noCiii,noCijk,ateqb,getQffplot,writeQFF,printaverageOK
  LOGICAL::justCiii,justQiiii,justCiij,justQiijj,justQiiij,justCijk,justQiijk,justQijkl

  ! real*8 :: Gi(Ndof),Hii(Ndof),Ciii(Ndof),Qiiii(Ndof)
  ! real*8 :: Hij(Ndof,Ndof),Ciij(Ndof,Ndof),Qiijj(Ndof,Ndof),Qiiij(Ndof,Ndof)
  ! real*8 :: Cijk(Ndof,Ndof,Ndof),Qiijk(Ndof,Ndof,Ndof)
  ! real*8 :: Qijkl(Ndof,Ndof,Ndof,Ndof)
CONTAINS

  SUBROUTINE get_QFF
    USE global
    IMPLICIT NONE
    IF(debug)PRINT*,"get_QFF"
    CALL init_QFF()
    IF(.NOT.allzero)THEN
      SELECT CASE(qfftype)
        CASE(1)
          CALL read_Sindo()
          IF(writeQFF)CALL write_QFF()
          CALL put_coefs()
          IF(unittype==1) CALL convert2atomicunits()
        CASE(2)
          CALL model1_QFF()
        CASE(3)
          CALL read_Midas()
          IF(writeQFF)then
          call convertfromatomicunits()
          call del_coefs()
          CALL write_QFF()
          STOP
          endif
        CASE(4)
          CALL chain_QFF()
          IF(writeQFF)CALL write_QFF()
          CALL put_coefs()
        CASE(5)
          print*, "hrm1",hrmfreq

          CALL read_sindo()
          print*, "hrm2",hrmfreq

          call get_freq()
          IF(writeQFF)CALL write_QFF()
          CALL put_coefs()
      END SELECT
    ENDIF
    IF(printaverageOK)call print_average()
    IF(runmonomer)CALL monomer_QFF()
    IF(harmonic)CALL make_harmonic()
    IF(ateqb)CALL make_eqb()
    IF(quartic)CALL make_quartic()
    IF(nocubic)CALL make_nocubic()
    IF(noCiii)CALL make_noCiii()
    IF(noCiij)CALL make_noCiij()
    IF(noCijk)CALL make_noCijk()
    IF(justCiii)CALL subjust_Ciii()
    IF(justQiiii)CALL subjust_Qiiii()
    IF(justCiij)CALL subjust_Ciij()
    IF(justQiiij)CALL subjust_Qiiij()
    IF(justQiijj)CALL subjust_Qiijj()
    IF(justCijk)CALL subjust_Cijk()
    IF(justCiii)CALL subjust_Ciii()
    IF(justQiijk)CALL subjust_Qiijk()
    IF(justQijkl)CALL subjust_Qijkl()
    IF(fullquartic)CALL make_fullquartic()
    IF(getQFFplot)CALL plotquartic()
    IF(debug)PRINT*,"get_QFF :)"
  END SUBROUTINE get_QFF

  SUBROUTINE init_QFF
    USE global
    IMPLICIT NONE
    IF(debug)PRINT*,"init_QFF"
    ALLOCATE(Gi(Ndof),Hii(ndof),Ciii(ndof),Qiiii(ndof))
    Gi=0.d0;Hii=0.d0;Ciii=0.d0;Qiiii=0.d0
    IF (nMR .GT. 1) THEN
      ALLOCATE(Hij(ndof,ndof),Ciij(ndof,ndof),Qiijj(ndof,ndof),Qiiij(ndof,ndof))
      Hij=0.d0;Ciij=0.d0;Qiijj=0.d0;Qiiij=0.d0
    ENDIF
    IF (nMR .GT. 2) THEN
      ALLOCATE(Cijk(ndof,ndof,ndof),Qiijk(ndof,ndof,ndof))
      Cijk=0.d0;Qiijk=0.d0
    ENDIF
    IF (nMR .GT. 3) THEN
      ALLOCATE(Qijkl(ndof,ndof,ndof,ndof))
      Qijkl=0.d0
    ENDIF
    IF(debug)PRINT*,"init_QFF :)"
  END SUBROUTINE init_QFF

  SUBROUTINE model1_QFF
    USE global
    IMPLICIT NONE
    REAL*8::lambda,nu
    INTEGER::i

    lambda=-0.1116d0
    nu=0.08414d0
    DO i=1,Ndof-1,2
      Hii(i)=0.5d0*(0.29375d0)
      Hii(i+1)=0.5d0*(2.12581d0)
      Ciii(i)=lambda*nu
      IF(nMR>1) Ciij(i+1,i)=lambda
       !     hrmfreq(i)=dsqrt(0.29375d0)
       !     hrmfreq(i+1)=dsqrt(2.12581d0)
    ENDDO
    freq=hrmfreq
    PRINT*,"model QFF formed"
  END SUBROUTINE model1_QFF

  subroutine get_freq()
    use global
    implicit none
    if(debug)print*, "get_freq"
    print*, "hrm",hrmfreq
    hrmfreq=dsqrt(Hii)
    freq=hrmfreq
    if(debug)print*, "get_freq :)"
  end subroutine get_freq

  SUBROUTINE chain_hessian
    USE global
    USE lapack
    IMPLICIT NONE
    REAL*8,ALLOCATABLE::hessian(:,:),eigenvalues(:)
    INTEGER::i,tndof
    tndof=ndof+1
    ALLOCATE(eigenvalues(tndof),hessian(tndof,tndof))
    hessian=0.d0
    DO i=1,tndof
      IF (i==1 .OR. i==tndof) THEN
        hessian(i,i)=1.d0
      ELSE
        hessian(i,i)=2.d0
      ENDIF
      IF (i<tndof) hessian(i,i+1)=-1.d0
    ENDDO!i
    CALL diag(tndof,hessian,eigenvalues)
    hrmfreq=eigenvalues(2:tndof)
    freq=hrmfreq
    if (debug) PRINT*,"Hrmfreq for chain is ready"
if (debug) Print*,"chain_hessian :)"
  END SUBROUTINE chain_hessian

  SUBROUTINE chain_QFF
    USE global
    USE lapack
    IMPLICIT NONE
    INTEGER::i,m1,m2,m3,m4,tndof
    REAL*8,PARAMETER::thresh=1.d-10
    REAL*8,ALLOCATABLE::hessian(:,:),eigenvalues(:)
    REAL*8::tmp
    tndof=ndof+1
    ALLOCATE(eigenvalues(tndof),hessian(tndof,tndof))
    hessian=0.d0
    DO i=1,tndof
      IF (i==1 .OR. i==tndof) THEN
        hessian(i,i)=1.d0
      ELSE
        hessian(i,i)=2.d0
      ENDIF
      IF (i<tndof) hessian(i,i+1)=-1.d0
    ENDDO!i
    CALL diag(tndof,hessian,eigenvalues)
    Hii=eigenvalues(2:tndof)
    hrmfreq=dsqrt(Hii)
    freq=hrmfreq
    !Fijk cubic NC force constants
    DO m1=1,ndof
      DO m2=1,ndof
        DO m3=1,ndof
          tmp=0.d0
          DO i=1,tndof
            IF(i==1)THEN
              !Fiii=-1 for i=1
              tmp=tmp-hessian(i,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)
            ELSE IF (i==tndof) THEN
              !Fiii=1 for i=N
              tmp=tmp+hessian(i,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)
            ENDIF
            IF(i<tndof)THEN
              !Fijj=-1 j=i+1
              tmp=tmp-hessian(i,m1+1)*hessian(i+1,m2+1)*hessian(i+1,m3+1)
              tmp=tmp-hessian(i+1,m1+1)*hessian(i,m2+1)*hessian(i+1,m3+1)
              tmp=tmp-hessian(i+1,m1+1)*hessian(i+1,m2+1)*hessian(i,m3+1)
              !Fiij=1 j=i+1
              tmp=tmp+hessian(i,m1+1)*hessian(i,m2+1)*hessian(i+1,m3+1)
              tmp=tmp+hessian(i,m1+1)*hessian(i+1,m2+1)*hessian(i,m3+1)
              tmp=tmp+hessian(i+1,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)
            ENDIF
          ENDDO!i
          IF(abs(tmp)<thresh)   tmp=0.d0
          tmp=tmp*kcubic
          IF(m1==m2 .AND. m2==m3)THEN
            Ciii(m1)=tmp
          ELSE IF (nMR>1)THEN
            IF(m1==m2)  THEN
              Ciij(m1,m3)=tmp
            ELSE IF (m1==m3) THEN
              Ciij(m1,m2)=tmp
            ELSE IF (m2==m3) THEN
              Ciij(m2,m1)=tmp
            ELSE IF (nMR>2)THEN
              Cijk(m1,m2,m3)=tmp
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !Fijkl quartic NC force constants
    DO m1=1,ndof
      DO m2=1,ndof
        DO m3=1,ndof
          DO m4=1,ndof
            tmp=0.d0
            DO i=1,tndof
              IF(i==1 .OR. i==tndof)THEN
                !Fiiii for i=1,N
                tmp=tmp+hessian(i,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)*hessian(i,m4+1)
              ELSE
                !Fiiii for 1<i<N
                tmp=tmp+2.d0*hessian(i,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)*hessian(i,m4+1)
              ENDIF
              IF(i<tndof)THEN
                !Fiijj for |j-i|=1
                tmp=tmp+hessian(i+1,m1+1)*hessian(i+1,m2+1)*hessian(i,m3+1)*hessian(i,m4+1)
                tmp=tmp+hessian(i+1,m1+1)*hessian(i,m2+1)*hessian(i+1,m3+1)*hessian(i,m4+1)
                tmp=tmp+hessian(i+1,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)*hessian(i+1,m4+1)
                tmp=tmp+hessian(i,m1+1)*hessian(i+1,m2+1)*hessian(i+1,m3+1)*hessian(i,m4+1)
                tmp=tmp+hessian(i,m1+1)*hessian(i+1,m2+1)*hessian(i,m3+1)*hessian(i+1,m4+1)
                tmp=tmp+hessian(i,m1+1)*hessian(i,m2+1)*hessian(i+1,m3+1)*hessian(i+1,m4+1)
                !Fijjj j=i+1
                tmp=tmp-hessian(i,m1+1)*hessian(i+1,m2+1)*hessian(i+1,m3+1)*hessian(i+1,m4+1)
                tmp=tmp-hessian(i+1,m1+1)*hessian(i,m2+1)*hessian(i+1,m3+1)*hessian(i+1,m4+1)
                tmp=tmp-hessian(i+1,m1+1)*hessian(i+1,m2+1)*hessian(i,m3+1)*hessian(i+1,m4+1)
                tmp=tmp-hessian(i+1,m1+1)*hessian(i+1,m2+1)*hessian(i+1,m3+1)*hessian(i,m4+1)
                !Fiiij
                tmp=tmp-hessian(i+1,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)*hessian(i,m4+1)
                tmp=tmp-hessian(i,m1+1)*hessian(i+1,m2+1)*hessian(i,m3+1)*hessian(i,m4+1)
                tmp=tmp-hessian(i,m1+1)*hessian(i,m2+1)*hessian(i+1,m3+1)*hessian(i,m4+1)
                tmp=tmp-hessian(i,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)*hessian(i+1,m4+1)
              ENDIF
            ENDDO!i
            IF(abs(tmp)<thresh)   tmp=0.d0
            tmp=tmp*kquartic
            IF(m1==m2 .AND. m2==m3 .AND. m3==m4)THEN
              Qiiii(m1)=tmp
            ELSE IF (nMR>1) THEN
              IF(m1==m2 .AND. m3==m4)  THEN
                Qiijj(m1,m3)=tmp
              ELSE IF (m1==m3.AND. m2==m4) THEN
                Qiijj(m1,m2)=tmp
              ELSE IF (m1==m4.AND. m2==m3) THEN
                Qiijj(m1,m2)=tmp
              ELSE IF (m1==m2.AND. m2==m3) THEN
                Qiiij(m1,m4)=tmp
              ELSE IF (m1==m2.AND. m2==m4) THEN
                Qiiij(m1,m3)=tmp
              ELSE IF (m1==m3.AND. m3==m4) THEN
                Qiiij(m1,m2)=tmp
              ELSE IF (m4==m3.AND. m2==m3) THEN
                Qiiij(m4,m1)=tmp
              ELSE IF (nMR>2)THEN
                IF(m1==m2)  THEN
                  Qiijk(m1,m3,m4)=tmp
                ELSE IF (m1==m3) THEN
                  Qiijk(m1,m2,m4)=tmp
                ELSE IF (m1==m4) THEN
                  Qiijk(m1,m2,m3)=tmp
                ELSE IF (m2==m3) THEN
                  Qiijk(m2,m1,m4)=tmp
                ELSE IF (m2==m4) THEN
                  Qiijk(m2,m1,m3)=tmp
                ELSE IF (m3==m4) THEN
                  Qiijk(m3,m1,m2)=tmp
                ELSE IF (nMR>3)THEN
                  Qijkl(m1,m2,m3,m4)=tmp
                ENDIF
              ENDIF
            ELSE
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    PRINT*,"model QFF formed"
  END SUBROUTINE chain_QFF
  !
  !  real*8 function Fxij(i1,i2)
  !    use global
  !    implicit none
  !    integer::i1,i2,tndof
  !    tndof=ndof+1
  !    IF(i1==i2)THEN
  !      !Fii=1 for i=1,N
  !      if(i1==1 .or. i1==tndof)then
  !        Fxij=1.d0
  !        return
  !      else
  !        Fxij=2.d0
  !        return
  !      endif
  !    Else if (abs(i1-i2)==1)then
  !      Fxij=-1.d0
  !      return
  !    endif
  !  end function
  !
  !  real*8 function Fxijk(i1,i2,i3)
  !    use global
  !    implicit none
  !    integer::i1,i2,i3,tndof
  !    Fxijk=0.d0
  !    tndof=ndof+1
  !    IF(i1==i2 .and. i2==i3)THEN
  !      !Fii=1 for i=1,N
  !      if(i1==1)then
  !        Fxijk=-1.d0
  !        return
  !      else if(i1==tndof) then
  !        Fxijk=1.d0
  !        return
  !      else
  !        return
  !      endif
  !    Else if (i1==i2 .and. abs(i1-i3)==1)then
  !      Fxijk=-1.d0
  !      return
  !    endif
  !
  !                !Fiii=-1 for i=1
  !                !Fiii=1 for i=N
  !
  !                !Fijj=-1 j=i+1
  !
  !                !Fiij=1 j=i+1
  !  end function

  SUBROUTINE chain_QFF2
    USE global
    USE lapack
    IMPLICIT NONE
    INTEGER::i,m1,m2,m3,m4,tndof
    REAL*8,PARAMETER::thresh=1.d-10
    REAL*8,ALLOCATABLE::hessian(:,:),eigenvalues(:)
    REAL*8::tmp
    tndof=ndof+1
    ALLOCATE(eigenvalues(tndof),hessian(tndof,tndof))
    hessian=0.d0
    DO i=1,tndof
      IF (i==1 .OR. i==tndof) THEN
        hessian(i,i)=1.d0
      ELSE
        hessian(i,i)=2.d0
      ENDIF
      IF (i<tndof) hessian(i,i+1)=-1.d0
    ENDDO!i
    CALL diag(tndof,hessian,eigenvalues)
    Hii=eigenvalues(2:tndof)
    hrmfreq=dsqrt(Hii)
    freq=hrmfreq

    DO m1=1,ndof
      DO m2=1,ndof
        DO m3=1,ndof
          tmp=0.d0
          DO i=1,tndof
            tmp=tmp-hessian(i,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)
          ENDDO!i
          IF(abs(tmp)<thresh)   tmp=0.d0
          IF(m1==m2 .AND. m2==m3)THEN
            Ciii(m1)=tmp
          ELSE IF (nMR>1)THEN
            IF(m1==m2)  THEN
              Ciij(m1,m3)=tmp
            ELSE IF (m1==m3) THEN
              Ciij(m1,m2)=tmp
            ELSE IF (m2==m3) THEN
              Ciij(m2,m1)=tmp
            ENDIF
          ELSE IF (nMR>2)THEN
            Cijk(m1,m2,m3)=tmp
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DO m1=1,ndof
      DO m2=1,ndof
        DO m3=1,ndof
          DO m4=1,ndof
            tmp=0.d0
            DO i=1,tndof
              tmp=tmp-hessian(i,m1+1)*hessian(i,m2+1)*hessian(i,m3+1)*hessian(i,m4+1)
            ENDDO!i
            IF(abs(tmp)<thresh)   tmp=0.d0
            IF(m1==m2 .AND. m2==m3 .AND. m3==m4)THEN
              Qiiii(m1)=tmp
            ELSE IF (nMR>1) THEN
              IF(m1==m2 .AND. m3==m4)  THEN
                Qiijj(m1,m3)=tmp
              ELSE IF (m1==m3.AND. m2==m4) THEN
                Qiijj(m1,m2)=tmp
              ELSE IF (m1==m4.AND. m2==m3) THEN
                Qiijj(m1,m2)=tmp
              ELSE IF (m1==m2.AND. m2==m3) THEN
                Qiiij(m1,m4)=tmp
              ELSE IF (m1==m2.AND. m2==m4) THEN
                Qiiij(m1,m3)=tmp
              ELSE IF (m1==m3.AND. m3==m4) THEN
                Qiiij(m1,m2)=tmp
              ELSE IF (m4==m3.AND. m2==m3) THEN
                Qiiij(m4,m1)=tmp
              ELSE IF (nMR>2)THEN
                IF(m1==m2.AND. m3/=m4)  THEN
                  Qiijk(m1,m3,m4)=tmp
                ELSE IF (m1==m3.AND. m2/=m4) THEN
                  Qiijk(m1,m2,m4)=tmp
                ELSE IF (m1==m4.AND. m2/=m3) THEN
                  Qiijk(m1,m2,m3)=tmp
                ELSE IF (m2==m3.AND. m1/=m4) THEN
                  Qiijk(m2,m1,m4)=tmp
                ELSE IF (m2==m4.AND. m1/=m3) THEN
                  Qiijk(m2,m1,m3)=tmp
                ELSE IF (m3==m4.AND. m1/=m2) THEN
                  Qiijk(m3,m1,m2)=tmp
                ELSE IF (nMR>3)THEN
                  Qijkl(m1,m2,m3,m4)=tmp
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    PRINT*,"model QFF formed"
  END SUBROUTINE chain_QFF2

  SUBROUTINE read_Sindo
    USE global
    IMPLICIT NONE
    INTEGER :: i,j,k,l,n,tmp

    if (debug)print*,"read_Sindo"
    OPEN(UNIT=7,FILE=TRIM(fcfile),STATUS='OLD',ACTION='READ')
    READ(7,*,END=666)!# Energy / hartree
    READ(7,'(3x,D20.12)',END=666) Eref
    READ(7,*,END=666)!# Geometry / Angs amu1/2
    !for model calculations
    Ndof=Ndof/Nmonomer
    j=MOD(Ndof,3)
    IF(j==0) THEN
      n=Ndof/3
      DO i=1,n
        READ(7,*,END=666) (geo(3*(i-1)+j),j=1,3)
      END DO
    ELSE
      n=(Ndof-j)/3
      DO i=1,n
        READ(7,*,END=666) (geo(3*(i-1)+j),j=1,3)
      END DO
      READ(7,*,END=666) (geo(j),j=3*n+1,Ndof)
    ENDIF

    !   >> >> 1 MR << <<
    READ(7,'(5x,a)') title1MR
    READ(7,*,END=666)!# Gradient / hartree Angs^-1 amu^-1/2
    DO n=1,Ndof
      READ(7,100)i, Gi(n)
    END DO
    READ(7,*,END=666)!# Hessian(i,i) / hartree Angs^-2 amu^-1
    DO n=1,Ndof
      READ(7,100,END=666) i,Hii(n)
    END DO
    READ(7,*,END=666)!# Cubic(i,i,i) / hartree Angs^-3 amu^-3/2
    DO n=1,Ndof
      READ(7,100,END=666)i, Ciii(n)
    END DO
    READ(7,*,END=666)!# Quartic(i,i,i,i) / hartree Angs^-4 amu^-2
    DO n=1,Ndof
      READ(7,100,END=666)i, Qiiii(n)
    END DO
100 FORMAT(i4,D20.10)

    !   >> >> 2 MR << <<
    IF(nMR .GT. 1) THEN
      tmp = ndof * (ndof - 1) / 2

      READ(7,'(5x,a)',END=222) title2MR

      READ(7,*,END=222)!# Hessian(i,j) / hartree Angs^-2 amu^-1
      DO n=1,tmp
        READ(7,200,END=666)i,j, Hij(i,j)
        Hij(j,i)=Hij(i,j)
      END DO

      READ(7,*,END=666)!# Quartic(i,i,j,j) / hartree Angs^-4 amu^-2
      DO n=1,tmp
        READ(7,200,END=666)i,j, Qiijj(i,j)
        Qiijj(j,i)=Qiijj(i,j)
      END DO

      READ(7,*,END=666)!# Cubic(i,i,j) / hartree Angs^-3 amu^-3/2
      DO n=1,2*tmp
        READ(7,200,END=666)i,j, Ciij(i,j)
      END DO

      READ(7,*,END=666)!# Quartic(i,i,i,j) / hartree Angs^-4 amu^-2
      DO n=1,2*tmp
        READ(7,200,END=666)i,j, Qiiij(i,j)
      END DO
200   FORMAT(2i4,D20.10)

      !   >> >> 3 MR << <<
      IF(nMR .GT. 2) THEN
        tmp = ndof * (ndof - 1) * (ndof - 2) / 6

300     FORMAT(3i4,D20.10)
        READ(7,'(5x,a)',END=222) title3MR
        READ(7,*,END=333)!# Cubic(i,j,k) / hartree Angs^-3 amu^-3/2
        DO n=1,tmp
          READ(7,300,END=666)i,j,k, Cijk(i,j,k)
          Cijk(i,k,j)=Cijk(i,j,k)
          Cijk(j,i,k)=Cijk(i,j,k)
          Cijk(j,k,i)=Cijk(i,j,k)
          Cijk(k,i,j)=Cijk(i,j,k)
          Cijk(k,j,i)=Cijk(i,j,k)
        END DO

        READ(7,*)!# Quartic(i,i,j,k) / hartree Angs^-4 amu^-2
        DO n=1,3*tmp
          READ(7,300,END=666)i,j,k, Qiijk(i,j,k)
          Qiijk(i,k,j)=Qiijk(i,j,k)
        END DO

        IF(nMR .GT. 3) THEN
          tmp = ndof * (ndof - 1) * (ndof - 2)* (ndof - 3) / 24

400       FORMAT(4i4,D20.10)
          READ(7,'(5x,a)',END=222) title3MR
          READ(7,*,END=333)!# Quartic(i,j,k,l) / hartree Angs^-4 amu^-2
          DO n=1,tmp
            READ(7,400,END=666)i,j,k,l, Qijkl(i,j,k,l)
            Qijkl(i,j,l,k)=Qijkl(i,j,k,l)
            Qijkl(i,l,j,k)=Qijkl(i,j,k,l)
            Qijkl(i,l,k,j)=Qijkl(i,j,k,l)
            Qijkl(i,k,j,l)=Qijkl(i,j,k,l)
            Qijkl(i,k,l,j)=Qijkl(i,j,k,l)
            Qijkl(j,i,k,l)=Qijkl(i,j,k,l)
            Qijkl(j,i,l,k)=Qijkl(i,j,k,l)
            Qijkl(j,k,i,l)=Qijkl(i,j,k,l)
            Qijkl(j,k,l,i)=Qijkl(i,j,k,l)
            Qijkl(j,l,k,i)=Qijkl(i,j,k,l)
            Qijkl(j,l,i,k)=Qijkl(i,j,k,l)
            Qijkl(k,j,i,l)=Qijkl(i,j,k,l)
            Qijkl(k,j,l,i)=Qijkl(i,j,k,l)
            Qijkl(k,l,i,j)=Qijkl(i,j,k,l)
            Qijkl(k,l,j,i)=Qijkl(i,j,k,l)
            Qijkl(k,i,j,l)=Qijkl(i,j,k,l)
            Qijkl(k,i,l,j)=Qijkl(i,j,k,l)
            Qijkl(l,j,k,i)=Qijkl(i,j,k,l)
            Qijkl(l,j,i,k)=Qijkl(i,j,k,l)
            Qijkl(l,i,k,j)=Qijkl(i,j,k,l)
            Qijkl(l,i,j,k)=Qijkl(i,j,k,l)
            Qijkl(l,k,i,j)=Qijkl(i,j,k,l)
            Qijkl(l,k,j,i)=Qijkl(i,j,k,l)
          END DO
        ENDIF !nMR>3
      ENDIF !nMR>2
    ENDIF !nMR>1
    CLOSE(7)
    PRINT*,"Force constants read successfully from ", fcfile    !for model calculations
    Ndof=Ndof*Nmonomer
    if (debug)print*,"read_Sindo :)"
    RETURN

222 CONTINUE
    !  nMR=1
    RETURN
333 CONTINUE
    !  nMR=2
    RETURN
666 STOP 'error in QFF.pot'
  END SUBROUTINE read_Sindo

  SUBROUTINE read_Midas
    USE global
    IMPLICIT NONE
    INTEGER::i,j,k,l,EOF
    REAL*8::fcvalue
    CHARACTER::line*70

    OPEN(UNIT=7,FILE=TRIM(fcfile),STATUS='OLD',ACTION='READ')
    READ(7,'(A)',IOSTAT=EOF)line
    DO i=1,ndof
      READ(7,*)hrmfreq(i)
    ENDDO
    freq=hrmfreq
    READ(7,'(A)',IOSTAT=EOF)line
    EOF=0
    DO WHILE(EOF==0)
      READ(7,'(A)',IOSTAT=EOF)line
      IF (LEN(TRIM(line)) >47) THEN
        READ(line,*)fcvalue,i,j,k,l
        IF(i==j .AND. j==k .AND. k==l)THEN
          Qiiii(i)=fcvalue*hrmfreq(i)*hrmfreq(i)
        ELSE IF(i==j .AND. j==k .AND. nMR > 1 )THEN
          Qiiij(i,l)=fcvalue*hrmfreq(i)**(1.5d0)*dsqrt(hrmfreq(l))
        ELSE IF(j==k .AND. k==l .AND. nMR > 1 )THEN
          Qiiij(j,i)=fcvalue*hrmfreq(j)**(1.5d0)*dsqrt(hrmfreq(i))
        ELSE IF(i==j .AND. k==l .AND. nMR > 1)THEN
          Qiijj(i,k)=fcvalue*hrmfreq(i)*hrmfreq(k)
          Qiijj(k,i)=Qiijj(i,k)
        ELSE IF(i==j .AND. nMR > 2)THEN
          Qiijk(i,k,l)=fcvalue*hrmfreq(i)*dsqrt(hrmfreq(k))*dsqrt(hrmfreq(l))
          Qiijk(i,l,k)=Qiijk(i,k,l)
        ELSE IF(j==k .AND. nMR > 2)THEN
          Qiijk(j,i,l)=fcvalue*hrmfreq(j)*dsqrt(hrmfreq(i))*dsqrt(hrmfreq(l))
          Qiijk(j,l,i)=Qiijk(j,i,l)
        ELSE IF(l==k .AND. nMR > 2)THEN
          Qiijk(k,i,j)=fcvalue*hrmfreq(k)*dsqrt(hrmfreq(i))*dsqrt(hrmfreq(j))
          Qiijk(k,j,i)=Qiijk(k,i,j)
        ELSE IF (nMR>3)THEN
          Qijkl(i,j,k,l)=fcvalue*dsqrt(hrmfreq(i))*dsqrt(hrmfreq(j))*dsqrt(hrmfreq(k))*dsqrt(hrmfreq(l))
          Qijkl(i,j,l,k)=Qijkl(i,j,k,l)
          Qijkl(i,l,j,k)=Qijkl(i,j,k,l)
          Qijkl(i,l,k,j)=Qijkl(i,j,k,l)
          Qijkl(i,k,j,l)=Qijkl(i,j,k,l)
          Qijkl(i,k,l,j)=Qijkl(i,j,k,l)
          Qijkl(j,i,k,l)=Qijkl(i,j,k,l)
          Qijkl(j,i,l,k)=Qijkl(i,j,k,l)
          Qijkl(j,k,i,l)=Qijkl(i,j,k,l)
          Qijkl(j,k,l,i)=Qijkl(i,j,k,l)
          Qijkl(j,l,k,i)=Qijkl(i,j,k,l)
          Qijkl(j,l,i,k)=Qijkl(i,j,k,l)
          Qijkl(k,j,i,l)=Qijkl(i,j,k,l)
          Qijkl(k,j,l,i)=Qijkl(i,j,k,l)
          Qijkl(k,l,i,j)=Qijkl(i,j,k,l)
          Qijkl(k,l,j,i)=Qijkl(i,j,k,l)
          Qijkl(k,i,j,l)=Qijkl(i,j,k,l)
          Qijkl(k,i,l,j)=Qijkl(i,j,k,l)
          Qijkl(l,j,k,i)=Qijkl(i,j,k,l)
          Qijkl(l,j,i,k)=Qijkl(i,j,k,l)
          Qijkl(l,i,k,j)=Qijkl(i,j,k,l)
          Qijkl(l,i,j,k)=Qijkl(i,j,k,l)
          Qijkl(l,k,i,j)=Qijkl(i,j,k,l)
          Qijkl(l,k,j,i)=Qijkl(i,j,k,l)
        ENDIF
      ELSE IF  (LEN(TRIM(line)) >42) THEN
        READ(line,*)fcvalue,i,j,k
        IF(i==j .AND. j==k )THEN
          Ciii(i)=fcvalue*hrmfreq(i)**(1.5d0)
        ELSE IF(i==j .AND. nMR > 1 )THEN
          Ciij(i,k)=fcvalue*hrmfreq(i)*dsqrt(hrmfreq(k))
        ELSE IF(j==k .AND. nMR > 1 )THEN
          Ciij(j,i)=fcvalue*hrmfreq(j)*dsqrt(hrmfreq(i))
        ELSE IF(nMR > 2)THEN
          Cijk(i,j,k)=fcvalue*dsqrt(hrmfreq(i))*dsqrt(hrmfreq(j))*dsqrt(hrmfreq(k))
          Cijk(i,k,j)=Cijk(i,j,k)
          Cijk(j,i,k)=Cijk(i,j,k)
          Cijk(j,k,i)=Cijk(i,j,k)
          Cijk(k,i,j)=Cijk(i,j,k)
          Cijk(k,j,i)=Cijk(i,j,k)
        ENDIF
      ELSE IF  (LEN(TRIM(line)) > 37) THEN
        READ(line,*)fcvalue,i,j
        IF(i==j)THEN
          Hii(i)=fcvalue*hrmfreq(i)
        ELSE IF( nMR > 1 )THEN
          Hij(i,j)=fcvalue*dsqrt(hrmfreq(i))*dsqrt(hrmfreq(j))
          Hij(j,i)=Hij(i,j)
        ENDIF
      ELSE IF  (LEN(TRIM(line)) >32) THEN
        READ(line,*)fcvalue,i
        Gi(i)=fcvalue*dsqrt(hrmfreq(i))
      ELSE
        PRINT*,"Read(Midas) error in line: ",line
      ENDIF
    ENDDO
    CLOSE(7)
    PRINT*,"Force constants read successfully from ", fcfile
  END SUBROUTINE read_Midas

  SUBROUTINE monomer_QFF
    USE global, ONLY : nDOF,nMR,Nmonomer
    IMPLICIT NONE
    INTEGER :: m1,m2,m3,m4,n,mdof
    mdof=ndof/Nmonomer
    DO n=1,Nmonomer-1
      DO m1=1,mdof
        Gi(m1+n*mdof)=Gi(m1)
        Hii(m1+n*mdof)=Hii(m1)
        Ciii(m1+n*mdof)=Ciii(m1)
        Qiiii(m1+n*mdof)=Qiiii(m1)
        IF(nMR .GT. 1) THEN
          DO m2=1,mdof
            Hij(m1+n*mdof,m2+n*mdof)=Hij(m1,m2)
            Qiijj(m1+n*mdof,m2+n*mdof)=Qiijj(m1,m2)
            Ciij(m1+n*mdof,m2+n*mdof)=Ciij(m1,m2)
            Qiiij(m1+n*mdof,m2+n*mdof)=Qiiij(m1,m2)
            IF(nMR .GT. 2) THEN
              DO m3=1,mdof
                Cijk(m1+n*mdof,m2+n*mdof,m3+n*mdof)=Cijk(m1,m2,m3)
                Qiijk(m1+n*mdof,m2+n*mdof,m3+n*mdof)=Qiijk(m1,m2,m3)
                IF(nMR .GT. 3) THEN
                  DO m4=1,mdof
                    Qijkl(m1+n*mdof,m2+n*mdof,m3+n*mdof,m4+n*mdof)=Qijkl(m1,m2,m3,m4)
                  ENDDO !m4
                ENDIF !nMR=4
              ENDDO!m3
            ENDIF!nMR=3
          ENDDO!m2
        ENDIF!nMR=2
      ENDDO !m1
    ENDDO !n
  END SUBROUTINE monomer_QFF

  SUBROUTINE write_QFF
    USE global, ONLY : geo,nDOF,Eref,nMR
    IMPLICIT NONE
    INTEGER :: i,j,k,l,n,tmp
    OPEN(UNIT=8,FILE='fc.mavi')
    WRITE(8,1000)
    WRITE(8,'(3x,D20.12)') Eref
    WRITE(8,1010)!# Geometry / Angs amu1/2
    j=MOD(Ndof,3)
    IF(j==0) THEN
      n=Ndof/3
      DO i=1,n
        WRITE(8,999) (geo(3*(i-1)+j),j=1,3)
      END DO
    ELSE
      n=(Ndof-j)/3
      DO i=1,n
        WRITE(8,999) (geo(3*(i-1)+j),j=1,3)
      END DO
      WRITE(8,999) (geo(j),j=3*n+1,Ndof)
    ENDIF
999 FORMAT(3d18.7)
1000 FORMAT('# Energy / hartree ')
1001 FORMAT('# Dipole / debye ')
1010 FORMAT('# Geometry / Angs amu1/2')

    !   >> >> 1 MR << <<

    WRITE(8,101) 'title1MR'

    WRITE(8,110)!# Gradient / hartree Angs^-1 amu^-1/2
    DO n=1,Ndof
      WRITE(8,100)n, Gi(n)
    END DO

    WRITE(8,120)!# Hessian(i,i) / hartree Angs^-2 amu^-1
    DO n=1,Ndof
      WRITE(8,100)n, Hii(n)
    END DO

    WRITE(8,130)!# Cubic(i,i,i) / hartree Angs^-3 amu^-3/2
    DO n=1,Ndof
      WRITE(8,100)n, Ciii(n)
    END DO

    WRITE(8,140)!# Quartic(i,i,i,i) / hartree Angs^-4 amu^-2
    DO n=1,Ndof
      WRITE(8,100)n, Qiiii(n)
    END DO

100 FORMAT(i4,D20.10)
101 FORMAT('# 1MR ',a)
110 FORMAT('# Gradient / hartree Angs^-1 amu^-1/2 ')
111 FORMAT('# Gradient / debye Angs^-1 amu^-1/2 ')
120 FORMAT('# Hessian(i,i) / hartree Angs^-2 amu^-1 ')
121 FORMAT('# Hessian(i,i) / debye Angs^-2 amu^-1 ')
130 FORMAT('# Cubic(i,i,i) / hartree Angs^-3 amu^-3/2 ')
131 FORMAT('# Cubic(i,i,i) / debye Angs^-3 amu^-3/2 ')
140 FORMAT('# Quartic(i,i,i,i) / hartree Angs^-4 amu^-2 ')
141 FORMAT('# Quartic(i,i,i,i) / debye Angs^-4 amu^-2 ')

        !   >> >> 2 MR << <<
    if(nMR>1)then
      tmp = ndof * (ndof - 1) / 2
      WRITE(8,201) 'title2MR'

      WRITE(8,210)!# Hessian(i,j) / hartree Angs^-2 amu^-1
      DO i=2,ndof
        DO j=1,i-1
          WRITE(8,200)i,j, Hij(i,j)
        ENDDO
      END DO

      WRITE(8,220)!# Quartic(i,i,j,j) / hartree Angs^-4 amu^-2
      DO i=2,ndof
        DO j=1,i-1
          WRITE(8,200)i,j, Qiijj(i,j)
        END DO
      END DO

      WRITE(8,230)!# Cubic(i,i,j) / hartree Angs^-3 amu^-3/2
      DO i=2,ndof
        DO j=1,i-1
          WRITE(8,200)i,j, Ciij(i,j)
          WRITE(8,200)j,i, Ciij(j,i)
        END DO
      END DO

      WRITE(8,240)!# Quartic(i,i,i,j) / hartree Angs^-4 amu^-2
      DO i=2,ndof
        DO j=1,i-1
          WRITE(8,200)i,j, Qiiij(i,j)
          WRITE(8,200)j,i, Qiiij(j,i)
        END DO
      END DO

200   FORMAT(2i4,D20.10)
201   FORMAT('# 2MR ',a)
210   FORMAT('# Hessian(i,j) / hartree Angs^-2 amu^-1 ')
211   FORMAT('# Hessian(i,j) / debye Angs^-2 amu^-1 ')
220   FORMAT('# Quartic(i,i,j,j) / hartree Angs^-4 amu^-2 ')
221   FORMAT('# Quartic(i,i,j,j) / debye Angs^-4 amu^-2 ')
230   FORMAT('# Cubic(i,i,j) / hartree Angs^-3 amu^-3/2 ')
231   FORMAT('# Cubic(i,i,j) / debye Angs^-3 amu^-3/2 ')
240   FORMAT('# Quartic(i,i,i,j) / hartree Angs^-4 amu^-2 ')
241   FORMAT('# Quartic(i,i,i,j) / debye Angs^-4 amu^-2 ')
250   FORMAT(2i4,D20.10)

         !   >> >> 3 MR << <<
      if(nMR>2)then
        tmp = ndof * (ndof - 1) * (ndof - 2) / 6
        WRITE(8,301) 'title3MR'

        WRITE(8,310)!# Cubic(i,j,k) / hartree Angs^-3 amu^-3/2
        DO i=3,Ndof
          DO j=2,i-1
            DO k=1,j-1
              WRITE(8,300)i,j,k, Cijk(i,j,k)
            END DO
          END DO
        END DO

        WRITE(8,320)!# Quartic(i,i,j,k) / hartree Angs^-4 amu^-2
        DO i=3,Ndof
          DO j=2,i-1
            DO k=1,j-1
              WRITE(8,300)i,j,k, Qiijk(i,j,k)
              WRITE(8,300)j,k,i, Qiijk(j,k,i)
              WRITE(8,300)k,i,j, Qiijk(k,i,j)
            END DO
          END DO
        END DO

300     FORMAT(3i4,D20.10)
301     FORMAT('# 3MR ',a)
310     FORMAT('# Cubic(i,j,k) / hartree Angs^-3 amu^-3/2 ')
311     FORMAT('# Cubic(i,j,k) / debye Angs^-3 amu^-3/2 ')
320     FORMAT('# Quartic(i,i,j,k) / hartree Angs^-4 amu^-2 ')
321     FORMAT('# Quartic(i,i,j,k) / debye Angs^-4 amu^-2 ')

        !   >> >> 4 MR << <<
        if(nMR>3)then
          WRITE(8,401) 'title4MR'

          WRITE(8,410)!# Quartic(i,j,k,l) / hartree Angs^-4 amu^-2

          DO i=4,Ndof
            DO j=3,i-1
              DO k=2,j-1
                DO l=1,k-1
                  WRITE(8,400)i,j,k,l, Qijkl(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
400       FORMAT(4i4,D20.10)
401       FORMAT('# 4MR ',a)
410       FORMAT('# Quartic(i,j,k,l) / hartree Angs^-4 amu^-2 ')
        endif
      endif
    endif
  END SUBROUTINE write_QFF

  SUBROUTINE make_zero
    USE global, ONLY: nMR
    IMPLICIT NONE
    Gi=0.d0
    Hii=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
      Hij=0.d0
      Ciij=0.d0
      Qiiij=0.d0
      Qiijj=0.d0
      IF(nMR .GT. 2) THEN
        Cijk=0.d0
        Qiijk=0.d0
        IF(nMR .GT. 3) Qijkl=0.d0
      ENDIF
    ENDIF
  END SUBROUTINE make_zero

  SUBROUTINE make_harmonic
    USE global, ONLY: nMR
    IMPLICIT NONE
    !   Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
      !      Hij=0.d0
      Ciij=0.d0
      Qiiij=0.d0
      Qiijj=0.d0
      IF(nMR .GT. 2) THEN
        Cijk=0.d0
        Qiijk=0.d0
        IF(nMR .GT. 3) Qijkl=0.d0
      ENDIF
    ENDIF
  END SUBROUTINE make_harmonic

  SUBROUTINE make_eqb
    USE constants
    USE global, ONLY: nMR,hrmfreq,freq,debug,convert
    IMPLICIT NONE
    REAL*8:: cnv2au
    ! print*,"Hii error in au:", Hii-(hrmfreq*hrmfreq*0.5d0)
    ! print*,"Hii error in bohr:", (Hii-hrmfreq*hrmfreq*0.5d0)*cnv2au*cnv2au
    IF(debug)THEN
      PRINT*,"Hrm freq error in cm-1:", (dsqrt(2.d0*Hii)-hrmfreq)*convert
      PRINT*,"Max hrm error in cm-1:", MAXVAL((dsqrt(2.d0*Hii)-hrmfreq)*convert),"for mode",&
      MAXLOC((dsqrt(2.d0*Hii)-hrmfreq)*convert)
      PRINT*,"Hii=hrmfreq*hrmfreq*0.5d0",hrmfreq*hrmfreq*cnv2au*cnv2au
    ENDIF
    Gi=0.d0
    Hii=hrmfreq*hrmfreq*0.5d0
    IF (nMR .GT. 1)Hij=0.d0
  END SUBROUTINE make_eqb

  SUBROUTINE make_nocubic
    USE global, ONLY: nMR
    IMPLICIT NONE
    !    Gi=0.0
    Ciii=0.0
    IF(nMR .GT. 1) Ciij=0.0
    IF(nMR .GT. 2) Cijk=0.0

  END SUBROUTINE make_nocubic

  SUBROUTINE make_noCiii
    IMPLICIT NONE
    Ciii=0.0
  END SUBROUTINE make_noCiii

  SUBROUTINE make_noCiij
    USE global, ONLY: nMR
    IMPLICIT NONE
    IF(nMR .GT. 1) Ciij=0.0
  END SUBROUTINE make_noCiij

  SUBROUTINE make_noCijk
    USE global, ONLY: nMR
    IMPLICIT NONE
    IF(nMR .GT. 2) Cijk=0.0
  END SUBROUTINE make_noCijk

  SUBROUTINE subjust_Gi
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Gi"
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Gi

  SUBROUTINE subjust_Hii
    USE global, ONLY: nMR
    IMPLICIT NONE
        print*,"Just Hii"
    Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Hii

    SUBROUTINE subjust_Ciii
    USE global, ONLY: nMR
    IMPLICIT NONE
        print*,"Just Ciii"
    Gi=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Ciii

  SUBROUTINE subjust_Qiiii
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Qiiii"
    Gi=0.d0
    Ciii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Qiiii

  SUBROUTINE subjust_Hij
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Hij"
    Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Ciij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Hij

   SUBROUTINE subjust_Ciij
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Ciij"
    Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Ciij

    SUBROUTINE subjust_Qiiij
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Qiiij"
    Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Qiiij

    SUBROUTINE subjust_Qiijj
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Qiijj"
    Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiiij=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Qiijj

  SUBROUTINE subjust_Cijk
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Cijk"
    Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Qiijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Cijk

    SUBROUTINE subjust_Qiijk
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Qiijk"
    Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          IF(nMR .GT. 3) Qijkl=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Qiijk

    SUBROUTINE subjust_Qijkl
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Just Qijkl"
    Gi=0.d0
    Ciii=0.d0
    Qiiii=0.d0
    IF(nMR .GT. 1) THEN
       Hij=0.d0
       Ciij=0.d0
       Qiiij=0.d0
       Qiijj=0.d0
       IF(nMR .GT. 2) THEN
          Cijk=0.d0
          Qiijk=0.d0
       ENDIF
    ENDIF
  END SUBROUTINE subjust_Qijkl

  SUBROUTINE make_quartic
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Quartic"
    !    Gi=0.0
    Ciii=0.0
    Qiiii=0.0
    IF(nMR .GT. 1) THEN
      !      Hij=0.0
      Ciij=0.0
      Qiiij=0.0
      IF(nMR .GT. 2) THEN
        Cijk=0.0
        Qiijk=0.0
        IF(nMR .GT. 3) Qijkl=0.d0
      ENDIF
    ENDIF
  END SUBROUTINE make_quartic

  SUBROUTINE make_fullquartic
    USE global, ONLY: nMR
    IMPLICIT NONE
    print*,"Full quartic"
    !    Gi=0.0
    Ciii=0.0
    IF(nMR .GT. 1) THEN
      !      Hij=0.0
      Ciij=0.0
      Qiiij=0.0
      IF(nMR .GT. 2) THEN
        Cijk=0.0
        Qiijk=0.0
        IF(nMR .GT. 3) Qijkl=0.d0
      ENDIF
    ENDIF
  END SUBROUTINE make_fullquartic


  SUBROUTINE put_coefs
    USE global, ONLY: nMR,debug
    IMPLICIT NONE
    ! Gi=Gi
    Hii=Hii * 0.5d0
    Ciii=Ciii  / 6.d0
    Qiiii=Qiiii  / 24.d0
    IF(nMR .GT. 1) THEN
      ! Hij=Hij
      Ciij=Ciij * 0.5d0
      Qiiij=Qiiij / 6.d0
      Qiijj=Qiijj  * 0.25d0
      IF(nMR .GT. 2) THEN
        !   Cijk=Cijk
        Qiijk=Qiijk  * 0.5d0
         !IF(nMR .GT. 3) THEN
         ! Qijkl=Qijkl
         !ENDIF
      ENDIF
    ENDIF
   if (debug) PRINT*, "QFF coefficients ok..."
  END SUBROUTINE put_coefs

  SUBROUTINE del_coefs
    USE global, ONLY: nMR
    IMPLICIT NONE
    ! Gi=Gi
    Hii=Hii * 2.d0
    Ciii=Ciii  * 6.d0
    Qiiii=Qiiii  * 24.d0
    IF(nMR .GT. 1) THEN
      ! Hij=Hij
      Ciij=Ciij * 2.0d0
      Qiiij=Qiiij * 6.d0
      Qiijj=Qiijj  * 4.0d0
      IF(nMR .GT. 2) THEN
        !   Cijk=Cijk
        Qiijk=Qiijk  * 2.0d0
         !IF(nMR .GT. 3) THEN
         ! Qijkl=Qijkl
         !ENDIF
      ENDIF
    ENDIF
    PRINT*, "QFF coefficients ok..."
  END SUBROUTINE del_coefs

  SUBROUTINE convert2atomicunits
    USE constants
    USE global, ONLY: nMR,hrmfreq,freq,debug,convert
    IMPLICIT NONE
    REAL*8:: cnv2au
    if(debug)print*,"convert2atomicunits"
    hrmfreq=hrmfreq*c_wn2h
    freq=hrmfreq
    cnv2au=c_angs2bohr*dsqrt(c_amu2emass)
    Gi=Gi / cnv2au
    Hii=Hii / cnv2au/ cnv2au
    Ciii=Ciii / cnv2au/ cnv2au/ cnv2au
    Qiiii=Qiiii / cnv2au/ cnv2au/ cnv2au/ cnv2au
    IF(nMR .GT. 1) THEN
      Hij=Hij / cnv2au/ cnv2au
      Ciij=Ciij / cnv2au/ cnv2au/ cnv2au
      Qiiij=Qiiij / cnv2au/ cnv2au/ cnv2au/ cnv2au
      Qiijj=Qiijj / cnv2au/ cnv2au/ cnv2au/ cnv2au
      IF(nMR .GT. 2) THEN
        Cijk=Cijk / cnv2au/ cnv2au/ cnv2au
        Qiijk=Qiijk / cnv2au/ cnv2au/ cnv2au/ cnv2au
      ENDIF
    ENDIF
    if(debug)print*,"convert2atomicunits"
  END SUBROUTINE convert2atomicunits

    SUBROUTINE convertfromatomicunits
    USE constants
    USE global, ONLY: nMR,hrmfreq,freq,debug,convert
    IMPLICIT NONE
    REAL*8:: cnv2au
    if(debug)print*,"convert2atomicunits"
 !   hrmfreq=hrmfreq*c_wn2h
 !   freq=hrmfreq
    cnv2au=c_angs2bohr*dsqrt(c_amu2emass)
    Gi=Gi * cnv2au
    Hii=Hii * cnv2au* cnv2au
    Ciii=Ciii * cnv2au* cnv2au* cnv2au
    Qiiii=Qiiii * cnv2au* cnv2au* cnv2au* cnv2au
    IF(nMR .GT. 1) THEN
      Hij=Hij * cnv2au* cnv2au
      Ciij=Ciij * cnv2au* cnv2au* cnv2au
      Qiiij=Qiiij * cnv2au* cnv2au* cnv2au* cnv2au
      Qiijj=Qiijj * cnv2au* cnv2au* cnv2au* cnv2au
      IF(nMR .GT. 2) THEN
        Cijk=Cijk * cnv2au* cnv2au* cnv2au
        Qiijk=Qiijk * cnv2au* cnv2au* cnv2au* cnv2au
      ENDIF
    ENDIF
    if(debug)print*,"convertfromatomicunits"
  END SUBROUTINE convertfromatomicunits

  SUBROUTINE plotquartic
    USE global, ONLY:ndof
    IMPLICIT NONE
    INTEGER::m1,m2
    CHARACTER :: filename*20,mod1*2,mod2*2
    DO m1=1,Ndof
      DO m2=m1+1,Ndof
        WRITE(mod1,'(i2)')m1
        WRITE(mod2,'(i2)')m2
        filename='mod'//TRIM(ADJUSTL(mod1))//'mod'//TRIM(ADJUSTL(mod2))//'.plot'
        OPEN(100,FILE=TRIM(filename))
        !   WRITE(100,*)'set term postscript eps enhanced color lw 3 "Arial" 14'
        !   WRITE(100,*)'set output "',trim(adjustl(filename)),'"'
        WRITE(100,*)'set xrange[-0.6:0.6]'
        WRITE(100,*)'set yrange[-0.6:0.6]'
        WRITE(100,*)'set pm3d'
        WRITE(100,*)'set hidden3d'
        WRITE(100,*)'set contour'
        WRITE(100,*)'set cntrparam levels 30'
        WRITE(100,*)'set isosample 11'
        !       Write(100,111) Gi(m1),Hii(m1),Ciii(m1),Qiiii(m1),Gi(m2),Hii(m2),Ciii(m2),Qiiii(m2)
        !      Write(100,122) Hij(m1,m2),Qiijj(m1,m2),Ciij(m1,m2),Ciij(m1,m2),Qiiij(m1,m2),Qiiij(m1,m2)
        WRITE(100,101) Hii(m1),Hii(m2),Qiijj(m1,m2)
        CLOSE(100)
         !       Call system('gnuplot plot.gnu')
      ENDDO
    END DO
101 FORMAT ('splot',f20.15,'*x**2+',f20.15,'*y**2+',f20.15,'*x*x*y*y with pm3d')
    !111 Format ('splot',f10.5,'*x+',f10.5,'*x**2+',f10.5,'*x**3+',f10.5,'*x**4+',f10.5,'*y+',f10.5,&
    !           '*y**2+',f10.5,'*y**3+',f10.5,'*y**4+\')
    !122 Format (f10.5,'*x*y+',f10.5,'*x*x*y*y+',f10.5,'*x*x*y+',f10.5,'*y*y*x+',f10.5,'*x**3*y+',&
    !          f10.5,'*y**3*x with pm3d')
  END SUBROUTINE plotquartic

  SUBROUTINE print_average
    USE global, ONLY:ndof
    IMPLICIT NONE
    REAL*8,PARAMETER::thresh=1.d-10
    !    INTEGER::m1,m2
    !    CHARACTER :: filename*20,cndof*4
    !        filename=avg.fc.//TRIM(ADJUSTL(cndof))
    !        OPEN(100,FILE=TRIM(filename))
    !        CLOSE(100)
    print*,"Mean absolute value of force constants"
    print*,sum(abs(Hii))/ndof
    print*,sum(abs(Ciii))/ndof,sum(abs(Ciii))/count(abs(Ciii).gt.thresh),count(abs(Ciii).gt.thresh)
    print*,sum(abs(Qiiii))/ndof,sum(abs(Qiiii))/count(abs(Qiiii).gt.thresh),count(abs(Qiiii).gt.thresh)

  END SUBROUTINE print_average

END MODULE modQFF
