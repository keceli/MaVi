MODULE integral
  IMPLICIT NONE
CONTAINS

  REAL*8 FUNCTION harmint(level,mode,power,diff)
    ! One dimensional integral for harmonic oscillator eigenfunctions.
    ! < level | Q_mode^power | level-diff >
    ! Power=-2 indicates kinetic term integral
    ! < level | -1/2 del^2/del(q)^2 | level-diff >
    ! USE constants
    USE Global,ONLY:freq
    IMPLICIT NONE
    INTEGER ::level,mode,power,diff
    REAL*8::n
    REAL*8::omega

    n=DBLE(level)
    !freq(mode) instead of hrmfreq for qVSCF (freq is updated for each iteration)
    omega = freq(mode)
    harmint=0.d0
    SELECT CASE (power)

    CASE(-2)
       IF (diff == 0)THEN
          harmint=0.5d0 * omega * (n + 0.5d0)
       ELSE IF (diff == 2)THEN
          harmint=-0.25d0 * omega * dSQRT((n - 1.d0)*n)
       ENDIF
    CASE(0)
       IF (diff==0) harmint=1.d0
    CASE (1)
       IF (diff == 1) THEN
          harmint=DSQRT(n * 0.5d0 / omega);
       END IF
    CASE (2)
       IF (diff == 0) THEN
          harmint=(n + 0.5d0) / omega;
       ELSE IF (diff == 2)THEN
          IF(n<1.d0)PRINT*,"n<1 for level,mode,power,diff",level,mode,power,diff
          harmint=0.5d0 * DSQRT(n * (n - 1.d0)) / omega;
       ENDIF
    CASE (3)
       IF (diff == 1)THEN
          harmint=3.d0 * DSQRT(n*n*n / (8.d0 * omega*omega*omega));
       ELSE IF (diff == 3) THEN
          IF(n<2.d0)PRINT*,"********************************************************n<2"
          harmint=DSQRT(0.125d0 * n * (n - 1.d0) * (n - 2.d0) / (omega*omega*omega));
       ENDIF
    CASE (4)
       IF (diff == 0)THEN
          harmint=0.75d0 * (n*n * 2.d0 + n * 2.d0 + 1.d0) / (omega*omega);
       ELSE IF (diff == 2)THEN
          IF(n<1.d0)PRINT*,"n<1 for level,mode,power,diff",level,mode,power,diff
          harmint= (n - 0.5d0) * DSQRT(n * (n - 1.d0)) / (omega*omega);
       ELSE IF (diff == 4)THEN
          IF(n<3.d0)PRINT*,"********************************************************n<3"
          harmint=0.25d0 * dsqrt(n * (n - 1.d0) * (n - 2.d0) * (n - 3.d0)) / (omega*omega);
       ENDIF
    CASE DEFAULT
       harmint=0.d0
       PRINT*,"harmint error for power= ",power
    END SELECT
    RETURN
  END FUNCTION harmint

  REAL*8 FUNCTION hrmint(mode,n1,n2,power)
    INTEGER ::mode,power,n1,n2
    hrmint=harmint(MAX(n1,n2),mode,power,ABS(n2-n1))
    RETURN
  END FUNCTION hrmint

  REAL*8 FUNCTION hrmint1MR(m1,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,bra(ndof),ket(ndof)
    INTEGER::b1,k1
    b1=bra(m1)
    k1=ket(m1)
    hrmint1MR=hrmint(m1,b1,k1,1) * Gi(m1)&
         + hrmint(m1,b1,k1,3) * Ciii(m1)&
         + hrmint(m1,b1,k1,4) * Qiiii(m1);
    RETURN
  END FUNCTION hrmint1MR

  REAL*8 FUNCTION hrmint2MR(m1,m2,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,bra(ndof),ket(ndof)
    INTEGER::b1,b2,k1,k2
    b1=bra(m1)
    k1=ket(m1)
    b2=bra(m2)
    k2=ket(m2)
    hrmint2MR=Hij(m1,m2)*hrmint(m1,b1,k1,1)*hrmint(m2,b2,k2,1)&
         + hrmint(m1,b1,k1,2) * hrmint(m2,b2,k2,2) * Qiijj(m1,m2)&
         + hrmint(m1,b1,k1,2) * hrmint(m2,b2,k2,1) * Ciij(m1,m2)&
         + hrmint(m1,b1,k1,1) * hrmint(m2,b2,k2,2) * Ciij(m2,m1)&
         + hrmint(m1,b1,k1,3) * hrmint(m2,b2,k2,1) * Qiiij(m1,m2)&
         + hrmint(m1,b1,k1,1) * hrmint(m2,b2,k2,3) * Qiiij(m2,m1);
    RETURN
  END FUNCTION hrmint2MR

  REAL*8 FUNCTION qhrmint2MR(m1,m2,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,bra(ndof),ket(ndof)
    INTEGER::b1,b2,k1,k2
    b1=bra(m1)
    k1=ket(m1)
    b2=bra(m2)
    k2=ket(m2)
    qhrmint2MR=Hij(m1,m2)*hrmint(m1,b1,k1,1)*hrmint(m2,b2,k2,1)&
         + hrmint(m1,b1,k1,2) * hrmint(m2,b2,k2,1) * Ciij(m1,m2)&
         + hrmint(m1,b1,k1,1) * hrmint(m2,b2,k2,2) * Ciij(m2,m1)&
         + hrmint(m1,b1,k1,3) * hrmint(m2,b2,k2,1) * Qiiij(m1,m2)&
         + hrmint(m1,b1,k1,1) * hrmint(m2,b2,k2,3) * Qiiij(m2,m1);
    RETURN
  END FUNCTION qhrmint2MR

  REAL*8 FUNCTION hrmint3MR(m1,m2,m3,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,m3,bra(ndof),ket(ndof)
    INTEGER::b1,b2,b3,k1,k2,k3
    b1=bra(m1)
    k1=ket(m1)
    b2=bra(m2)
    k2=ket(m2)
    b3=bra(m3)
    k3=ket(m3)
    hrmint3MR=hrmint(m1,b1,k1,1) * hrmint(m2,b2,k2,1) * hrmint(m3,b3,k3,1) * Cijk(m1,m2,m3)&
         +hrmint(m1,b1,k1,2) * hrmint(m2,b2,k2,1) * hrmint(m3,b3,k3,1) * Qiijk(m1,m2,m3)&
         +hrmint(m1,b1,k1,1) * hrmint(m2,b2,k2,2) * hrmint(m3,b3,k3,1) * Qiijk(m2,m1,m3)&
         +hrmint(m1,b1,k1,1) * hrmint(m2,b2,k2,1) * hrmint(m3,b3,k3,2) * Qiijk(m3,m1,m2)
    RETURN
  END FUNCTION hrmint3MR

  REAL*8 FUNCTION hrmint4MR(m1,m2,m3,m4,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,m3,m4,bra(ndof),ket(ndof)
    INTEGER::b1,b2,b3,b4,k1,k2,k3,k4
    b1=bra(m1)
    k1=ket(m1)
    b2=bra(m2)
    k2=ket(m2)
    b3=bra(m3)
    k3=ket(m3)
    b4=bra(m4)
    k4=ket(m4)
    hrmint4MR=hrmint(m1,b1,k1,1) * hrmint(m2,b2,k2,1) * hrmint(m3,b3,k3,1) &
         * hrmint(m4,b4,k4,1) * Qijkl(m1,m2,m3,m4)
    RETURN
  END FUNCTION hrmint4MR

  REAL*8 FUNCTION modalexp(mode,power,coef)
    ! 1D integral for a modal in HO basis set.
    USE global
    IMPLICIT NONE
    INTEGER ::mode,power,n1,n2,diff
    REAL*8::coef(nDOF,hrmbasis(mode)),coef1,coef2
    modalexp = 0.d0;
    DO n1 = 0,hrmbasis(mode)-1
       coef1=coef(mode,n1+1)
       !IF(ABS(coef1) .LT. inttol)CYCLE
       modalexp = modalexp + coef1 * coef1 * harmint(n1, mode, power, 0)
       DO n2 = n1+1,hrmbasis(mode)-1
          diff=n2-n1
          coef2=coef(mode,n2+1)
          IF(diff > power .AND. power .NE. -2 )CYCLE
          !        IF((ABS(coef2) .LT. inttol) .OR. (diff .GT. maxFF) )CYCLE
          modalexp = modalexp + 2.d0 * coef1 * coef2 * harmint(n2, mode, power, diff)
       ENDDO
    ENDDO
  END FUNCTION modalexp

  !  REAL*8 FUNCTION modalint(mode,s1,s2,power)
  !    ! 1D integral for a modal in HO basis set.
  !    USE global
  !    IMPLICIT NONE
  !    INTEGER ::mode,power,n1,n2,diff,s1,s2
  !    REAL*8::coef1,coef2
  !    modalint = 0.d0;
  !    DO n1 = 0,maxbasis(mode)-1
  !       coef1=vscfcoefs(mode,n1+1,s1+1)
  !       coef2=vscfcoefs(mode,n1+1,s2+1)
  !       !IF(ABS(coef1) .LT. inttol)CYCLE
  !       modalint = modalint + coef1 * coef2 * harmint(n1, mode, power, 0)
  !       DO n2 = n1+1,maxbasis(mode)-1
  !          diff=n2-n1
  !          IF(diff > power .AND. power .NE. -2 )CYCLE
  !          coef2=vscfcoefs(mode,n2+1,s2+1)
  !          modalint = modalint + 2.d0 * coef1 * coef2 * harmint(n2, mode, power, diff)
  !       ENDDO
  !    ENDDO
  !  END FUNCTION modalint

  REAL*8 FUNCTION modalint(mode,s1,s2,power)
    ! 1D integral for a modal in HO basis set.
    !to do: special case s1=s2, change s1 name
    USE global
    IMPLICIT NONE
    INTEGER ::mode,power,n1,n2,diff,s1,s2
    REAL*8::coef1,coef2
    modalint = 0.d0;
    DO n1 = 0,hrmbasis(mode)-1
       DO n2 = 0,hrmbasis(mode)-1
          diff=ABS(n2-n1)
          IF( (diff > power) .AND. (power .NE. -2) )CYCLE
          coef1=vscfcoefs(mode,n1+1,s1+1)
          coef2=vscfcoefs(mode,n2+1,s2+1)
          modalint = modalint + coef1 * coef2 * harmint(MAX(n1,n2), mode, power, diff)
       ENDDO
    ENDDO
  END FUNCTION modalint
!
  REAL*8 FUNCTION modalint1MRU(m1,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,bra(ndof),ket(ndof)
    INTEGER::b1,k1
    b1=bra(m1)
    k1=ket(m1)
    modalint1MRU=modalint(m1,b1,k1,1) * (Gi(m1)-Uone(m1))&
         + modalint(m1,b1,k1,2) * (Hii(m1)-Utwo(m1))&
         + modalint(m1,b1,k1,3) * (Ciii(m1)-Uthree(m1))&
         + modalint(m1,b1,k1,4) * (Qiiii(m1)-Ufour(m1));
    RETURN
  END FUNCTION modalint1MRU

  REAL*8 FUNCTION modalint1MR(m1,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,bra(ndof),ket(ndof)
    INTEGER::b1,k1
    b1=bra(m1)
    k1=ket(m1)
    modalint1MR=modalint(m1,b1,k1,1) * Gi(m1)&
         + modalint(m1,b1,k1,2) * Hii(m1)&
         + modalint(m1,b1,k1,3) * Ciii(m1)&
         + modalint(m1,b1,k1,4) * Qiiii(m1);
    RETURN
  END FUNCTION modalint1MR
!
  REAL*8 FUNCTION modalintU(m1,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,bra(ndof),ket(ndof)
    INTEGER::b1,k1
    b1=bra(m1)
    k1=ket(m1)
    modalintU=modalint(m1,b1,k1,1) * Uone(m1)&
         + modalint(m1,b1,k1,2) * Utwo(m1)&
         + modalint(m1,b1,k1,3) * Uthree(m1)&
         + modalint(m1,b1,k1,4) * Ufour(m1);
    RETURN
  END FUNCTION modalintU

  REAL*8 FUNCTION modalint2MR(m1,m2,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,bra(ndof),ket(ndof)
    INTEGER::b1,b2,k1,k2
    b1=bra(m1)
    k1=ket(m1)
    b2=bra(m2)
    k2=ket(m2)
    modalint2MR=Hij(m1,m2)*modalint(m1,b1,k1,1)*modalint(m2,b2,k2,1)&
         + modalint(m1,b1,k1,2) * modalint(m2,b2,k2,2) * Qiijj(m1,m2)&
         + modalint(m1,b1,k1,2) * modalint(m2,b2,k2,1) * Ciij(m1,m2)&
         + modalint(m1,b1,k1,1) * modalint(m2,b2,k2,2) * Ciij(m2,m1)&
         + modalint(m1,b1,k1,3) * modalint(m2,b2,k2,1) * Qiiij(m1,m2)&
         + modalint(m1,b1,k1,1) * modalint(m2,b2,k2,3) * Qiiij(m2,m1);
    RETURN
  END FUNCTION modalint2MR

  REAL*8 FUNCTION modalint3MR(m1,m2,m3,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,m3,bra(ndof),ket(ndof)
    INTEGER::b1,b2,b3,k1,k2,k3
    b1=bra(m1)
    k1=ket(m1)
    b2=bra(m2)
    k2=ket(m2)
    b3=bra(m3)
    k3=ket(m3)
    modalint3MR=modalint(m1,b1,k1,1) * modalint(m2,b2,k2,1) * modalint(m3,b3,k3,1) * Cijk(m1,m2,m3)&
         +modalint(m1,b1,k1,2) * modalint(m2,b2,k2,1) * modalint(m3,b3,k3,1) * Qiijk(m1,m2,m3)&
         +modalint(m1,b1,k1,1) * modalint(m2,b2,k2,2) * modalint(m3,b3,k3,1) * Qiijk(m2,m1,m3)&
         +modalint(m1,b1,k1,1) * modalint(m2,b2,k2,1) * modalint(m3,b3,k3,2) * Qiijk(m3,m1,m2)
    RETURN
  END FUNCTION modalint3MR

  REAL*8 FUNCTION modalint4MR(m1,m2,m3,m4,bra,ket)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,m3,m4,bra(ndof),ket(ndof)
    INTEGER::b1,b2,b3,b4,k1,k2,k3,k4
    b1=bra(m1)
    k1=ket(m1)
    b2=bra(m2)
    k2=ket(m2)
    b3=bra(m3)
    k3=ket(m3)
    b4=bra(m4)
    k4=ket(m4)
    modalint4MR=modalint(m1,b1,k1,1) * modalint(m2,b2,k2,1) * modalint(m3,b3,k3,1) &
         * modalint(m4,b4,k4,1) * Qijkl(m1,m2,m3,m4)
    RETURN
  END FUNCTION modalint4MR


  REAL*8 FUNCTION funEdiff1(m1,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    tmp=0.d0
    tmp=tmp+modalint1MRU(m1,bra,ket)!not necessary, zero
    IF(nMR>1)THEN
       DO m2=1,ndof
          IF(m2==m1)CYCLE
          tmp=tmp+modalint2MR(m1,m2,bra,ket)
          IF(nMR>2)THEN
             DO m3=m2+1,ndof
                IF(m3==m1)CYCLE
                tmp=tmp+modalint3MR(m1,m2,m3,bra,ket)
                IF(nMR>3)THEN
                   DO m4=m3+1,ndof
                      IF((m4==m1).OR.(m4==m2))CYCLE
                      tmp=tmp+modalint4MR(m1,m2,m3,m4,bra,ket)
                   ENDDO !m4
                ENDIF!nmr>3
             ENDDO !m3
          ENDIF!nmr>2
       ENDDO !m2
    ENDIF!nmr>1
    denom=VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)
    funEdiff1=tmp*tmp/denom
    RETURN
  END FUNCTION funEdiff1

  REAL*8 FUNCTION funEdiff2(m1,m2,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    tmp=0.d0
    tmp=tmp+modalint2MR(m1,m2,bra,ket)
    IF(nMR>2)THEN
       DO m3=1,ndof
          IF((m3==m1).OR.(m3==m2))CYCLE
          tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket)
          IF(nMR>3)THEN
             DO m4=m3+1,ndof
                IF((m4==m1).OR.(m4==m2))CYCLE
                tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)
             ENDDO !m4
          ENDIF!nmr>3
       ENDDO !m3
    ENDIF!nmr>2
    denom=VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)&
         +VSCFenergies(m2,bra(m2)+1)-VSCFenergies(m2,ket(m2)+1)
    funEdiff2=tmp*tmp/denom
    RETURN
  END FUNCTION funEdiff2

    REAL*8 FUNCTION funEdiff3(m1,m2,m3,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    tmp=0.d0
    tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket)
    IF(nMR>3)THEN
       DO m4=1,ndof
          IF((m4==m1).OR.(m4==m2) .OR. (m4==m3))CYCLE
          tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)! should be zero
       ENDDO !m4
    ENDIF!nmr>3
    denom=VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)&
         +VSCFenergies(m2,bra(m2)+1)-VSCFenergies(m2,ket(m2)+1)&
         +VSCFenergies(m3,bra(m3)+1)-VSCFenergies(m3,ket(m3)+1)
    funEdiff3=tmp*tmp/denom
    RETURN
  END FUNCTION funEdiff3

  REAL*8 FUNCTION funEdiff4(m1,m2,m3,m4,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    tmp=0.d0
    tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)
    denom=VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)&
         +VSCFenergies(m2,bra(m2)+1)-VSCFenergies(m2,ket(m2)+1)&
         +VSCFenergies(m3,bra(m3)+1)-VSCFenergies(m3,ket(m3)+1)&
         +VSCFenergies(m4,bra(m4)+1)-VSCFenergies(m4,ket(m4)+1)
    funEdiff4=tmp*tmp/denom
    RETURN
  END FUNCTION funEdiff4

  REAL*8 FUNCTION funHOEdiff1(m1,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    tmp=0.d0
    tmp=tmp+hrmint1MR(m1,bra,ket)!not necessary, zero
    IF(nMR>1)THEN
       DO m2=1,ndof
          IF(m2==m1)CYCLE
          tmp=tmp+hrmint2MR(m1,m2,bra,ket)
          IF(nMR>2)THEN
             DO m3=m2+1,ndof
                IF(m3==m1)CYCLE
                tmp=tmp+hrmint3MR(m1,m2,m3,bra,ket)
                IF(nMR>3)THEN
                   DO m4=m3+1,ndof
                      IF((m4==m1).OR.(m4==m2))CYCLE
                      tmp=tmp+hrmint4MR(m1,m2,m3,m4,bra,ket)
                   ENDDO !m4
                ENDIF!nmr>3
             ENDDO !m3
          ENDIF!nmr>2
       ENDDO !m2
    ENDIF!nmr>1
    denom=freq(m1)*(bra(m1)-ket(m1))
    funHOEdiff1=tmp*tmp/denom
    RETURN
  END FUNCTION funHOEdiff1

  REAL*8 FUNCTION funHOEdiff2(m1,m2,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    tmp=0.d0
    tmp=tmp+hrmint2MR(m1,m2,bra,ket)
    IF(nMR>2)THEN
       DO m3=1,ndof
          IF((m3==m1).OR.(m3==m2))CYCLE
          tmp=tmp+ hrmint3MR(m1,m2,m3,bra,ket)
          IF(nMR>3)THEN
             DO m4=m3+1,ndof
                IF((m4==m1).OR.(m4==m2))CYCLE
                tmp=tmp+ hrmint4MR(m1,m2,m3,m4,bra,ket)
             ENDDO !m4
          ENDIF!nmr>3
       ENDDO !m3
    ENDIF!nmr>2
    denom=freq(m1)*(bra(m1)-ket(m1))&
              +freq(m2)*(bra(m2)-ket(m2))
    funHOEdiff2=tmp*tmp/denom
    RETURN
  END FUNCTION funHOEdiff2

  REAL*8 FUNCTION funHOEdiff3(m1,m2,m3,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    tmp=0.d0
    tmp=tmp+ hrmint3MR(m1,m2,m3,bra,ket)
    IF(nMR>3)THEN
       DO m4=1,ndof
          IF((m4==m1).OR.(m4==m2) .OR. (m4==m3))CYCLE
          tmp=tmp+ hrmint4MR(m1,m2,m3,m4,bra,ket)! should be zero
       ENDDO !m4
    ENDIF!nmr>3
    denom=freq(m1)*(bra(m1)-ket(m1))&
                  +freq(m2)*(bra(m2)-ket(m2))&
                  +freq(m3)*(bra(m3)-ket(m3))
    funHOEdiff3=tmp*tmp/denom
    RETURN
  END FUNCTION funHOEdiff3

  REAL*8 FUNCTION funHOEdiff4(m1,m2,m3,m4,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    tmp=0.d0
    tmp=tmp+ hrmint4MR(m1,m2,m3,m4,bra,ket)
    denom=freq(m1)*(bra(m1)-ket(m1))&
                      +freq(m2)*(bra(m2)-ket(m2))&
                      +freq(m3)*(bra(m3)-ket(m3))&
                      +freq(m4)*(bra(m4)-ket(m4))

    funHOEdiff4=tmp*tmp/denom
    RETURN
  END FUNCTION funHOEdiff4

  REAL*8 FUNCTION getint3b(m1,m2)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER::m1,m2,m3
     getint3b=Ciij(m1,m2)*hrmint(m1,0,0,2)*hrmint(m2,0,1,1)
    RETURN
  END FUNCTION getint3b

  REAL*8 FUNCTION getint3(m1,m2,m3)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER::m1,m2,m3
     getint3=Cijk(m1,m2,m3)*hrmint(m1,0,0,1)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)
    RETURN
  END FUNCTION getint3

  REAL*8 FUNCTION getint4b(m1,m2,m3)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER::m1,m2,m3
     getint4b=Qiijk(m1,m2,m3)*hrmint(m1,0,0,2)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)
    RETURN
  END FUNCTION getint4b

  REAL*8 FUNCTION getint4(m1,m2,m3,m4)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4
     getint4=Qijkl(m1,m2,m3,m4)*hrmint(m1,0,0,1)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)*hrmint(m4,0,1,1)
    RETURN
  END FUNCTION getint4

END MODULE integral
