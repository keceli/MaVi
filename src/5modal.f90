MODULE modalintegral
  IMPLICIT NONE
  REAL*8,ALLOCATABLE :: modalint(:,:,:,:)
CONTAINS

  SUBROUTINE init_modalint
    USE global
    IMPLICIT NONE
    if(debug)print*,'init_modalint'
    ALLOCATE(modalint(ndof,0:maxbasis-1,0:maxbasis-1,0:4))
    if(debug)print*,'init_modalint :)'
  END SUBROUTINE init_modalint



  SUBROUTINE form_modalint
    USE global
    USE integral,ONLY:harmint
    IMPLICIT NONE
    INTEGER ::mode,power,n1,n2,diff,s1,s2,p
    REAL*8::coef1,coef2,tmp
    if(debug)print*,'form_modalint'
    modalint=0.d0
    tmp=0.d0
    DO p=0,4
      if (p==0)then
        power=-2
      else
        power=p
      endif
      DO mode=1,ndof
        DO s1 = 0,hrmbasis(mode)-1
          DO s2 = 0,hrmbasis(mode)-1
            DO n1 = 0,hrmbasis(mode)-1
              DO n2 = 0,hrmbasis(mode)-1
                diff=ABS(n2-n1)
                IF( (diff > power) .AND. (power .NE. -2) )CYCLE
                coef1=vscfcoefs(mode,n1+1,s1+1)
                coef2=vscfcoefs(mode,n2+1,s2+1)
                tmp = tmp + coef1 * coef2 * harmint(MAX(n1,n2), mode, power, diff)
              ENDDO !n2
            ENDDO !n1
            modalint(mode,s1,s2,p)=tmp
            tmp=0.d0
          ENDDO !s2
        ENDDO !s1
      ENDDO!mode
    ENDDO!power
    if(debug)print*,'form_modalint :)'
  END SUBROUTINE form_modalint

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
    *modalint(m4,b4,k4,1) * Qijkl(m1,m2,m3,m4)
    RETURN
  END FUNCTION modalint4MR

  REAL*8 FUNCTION funEdiff1(m1,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    denom = VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)
    if(abs(denom) < denomcut) then
      funEdiff1=0.d0
    else
      tmp=modalint1MRU(m1,bra,ket)
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
      funEdiff1=tmp*tmp/denom
    ENDIF
    RETURN
  END FUNCTION funEdiff1

  REAL*8 FUNCTION funEdiff1MR(m1,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    denom = VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)
    if(abs(denom) < denomcut) then
      funEdiff1MR=0.d0
    else
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
      funEdiff1MR=tmp*tmp/denom
    ENDIF
    RETURN
  END FUNCTION funEdiff1MR

  REAL*8 FUNCTION funEdiff2(m1,m2,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    denom = VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)&
    + VSCFenergies(m2,bra(m2)+1)-VSCFenergies(m2,ket(m2)+1)
    if(abs(denom) < denomcut) then
      funEdiff2=0.d0
    else
      tmp=modalint2MR(m1,m2,bra,ket)
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
      funEdiff2=tmp*tmp/denom
    ENDIF
    RETURN
  END FUNCTION funEdiff2

  REAL*8 FUNCTION funEdiff3(m1,m2,m3,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    denom = VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)&
    + VSCFenergies(m2,bra(m2)+1)-VSCFenergies(m2,ket(m2)+1)&
    + VSCFenergies(m3,bra(m3)+1)-VSCFenergies(m3,ket(m3)+1)
    if(abs(denom) < denomcut) then
      funEdiff3=0.d0
    else
      tmp=modalint3MR(m1,m2,m3,bra,ket)
      IF(nMR>3)THEN
        DO m4=1,ndof
          IF((m4==m1).OR.(m4==m2) .OR. (m4==m3))CYCLE
          tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)! should be zero
        ENDDO !m4
      ENDIF!nmr>3
      funEdiff3=tmp*tmp/denom
    ENDIF
    RETURN
  END FUNCTION funEdiff3

  REAL*8 FUNCTION funEdiff4(m1,m2,m3,m4,bra,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,bra(ndof),ket(ndof)
    REAL*8::tmp,denom
    denom = VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)&
    + VSCFenergies(m2,bra(m2)+1)-VSCFenergies(m2,ket(m2)+1)&
    + VSCFenergies(m3,bra(m3)+1)-VSCFenergies(m3,ket(m3)+1)&
    + VSCFenergies(m4,bra(m4)+1)-VSCFenergies(m4,ket(m4)+1)
    if(abs(denom) < denomcut) then
      funEdiff4=0.d0
    else
      tmp=modalint4MR(m1,m2,m3,m4,bra,ket)
      funEdiff4=tmp*tmp/denom
    ENDIF
    RETURN
  END FUNCTION funEdiff4

  real*8 function fullmodal1MR(bra,ket)
    ! M Dimensional integrals for VSCF eigenfunctions.
    USE global
    IMPLICIT NONE
    INTEGER ::ket(ndof),bra(ndof)
    INTEGER::m1,m2,m3,m4
    REAL*8::tmp
    tmp = 1.d0;
    DO m1=1,ndof
      tmp=tmp* modalint1MR(m1,bra,ket);
    ENDDO !m1
    fullmodal1MR =tmp
  END function fullmodal1MR

  real*8 function fullmodalkin(bra,ket)
    ! M Dimensional integrals for VSCF eigenfunctions.
    USE global
    IMPLICIT NONE
    INTEGER ::ket(ndof),bra(ndof)
    INTEGER::m1,m2,m3,m4
    REAL*8::tmp
    tmp = 1.d0;
    DO m1=1,ndof
      tmp=tmp* modalint(m1,bra(m1),ket(m1),0);
    ENDDO !m1
    fullmodalkin =tmp
  END function fullmodalkin

  real*8 function fullmodalint(bra,ket)
    ! M Dimensional integrals for VSCF eigenfunctions.
    USE global
    IMPLICIT NONE
    INTEGER ::ket(ndof),bra(ndof)
    INTEGER::m1,m2,m3,m4
    REAL*8::tmp
    tmp = 0.d0;
    DO m1=1,ndof
      tmp=tmp+ modalint(m1,bra(m1),ket(m1),0);      !kinetic
      tmp=tmp+ modalint1MR(m1,bra,ket);
      IF (nMR > 1)THEN
        DO m2=m1+1,ndof
          tmp=tmp+ modalint2MR(m1,m2,bra,ket);
          IF (nMR > 2)THEN
            DO m3=m2+1,ndof
              IF (m3 == m1)CYCLE;
              tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket);
              IF (nMR > 3)THEN
                DO m4=m3+1,ndof
                  IF (m4 == m1 .OR. m4 == m2 )CYCLE;
                  tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket);
                ENDDO !m4
              ENDIF!nmr>3
            ENDDO !m3
          ENDIF!nmr>2
        ENDDO !m2
      ENDIF!nmr>1
    ENDDO !m1
    fullmodalint =tmp
  END function fullmodalint

END MODULE modalintegral
