MODULE modVSCF
  IMPLICIT NONE
  REAL*8,ALLOCATABLE:: modalexp(:,:)
  REAL*8:: Evscfgs
  logical::plotU,printU,get1MRenergy,getUZero
CONTAINS

  SUBROUTINE init_VSCF
    USE global
    IMPLICIT NONE
    IF(debug)PRINT*,"init_VSCF"
    plotU=.false.
    printU=.false.
    get1MRenergy=.false.
    getUzero=.false.
    call read_vscf()
    IF(storeint)ALLOCATE(modalexp(ndof,4))
    IF(storepot)ALLOCATE(uZero(ndof),uOne(ndof),uTwo(ndof),uThree(ndof),uFour(ndof))
    IF(debug)PRINT*,"init_VSCF :)"
  END SUBROUTINE init_VSCF

  SUBROUTINE read_VSCF
    use global,only:debug
    IMPLICIT NONE
    integer::ierr
    NAMELIST /vscf/plotU,printU,get1MRenergy,getUzero
    IF(debug)PRINT*,"read_VSCF"
    REWIND(5)
    !Print*,'read'
    plotU=.false.
    printU=.false.
    get1mrenergy=.false.
    ierr=0
    !do while(ierr==0)
    READ(5,vscf,iostat=ierr,end=1001,err=1002)
    1002 print*,"vscf namelist iostat ",ierr
    1001 continue
    !    if(ierr>0)then
    !      print*,"Error reading in &vscf"
    !      STOP
    !    endif
    !  READ(5,vscf,iostat=ierr)
    !end do
    ! if (ierr.ne. 0)print*,"Error in reading &vscf "
    IF(debug)PRINT*,"read_VSCF :)"
  !100 print*,"Error in reading &vscf "
  !200 continue
  END SUBROUTINE read_VSCF

  SUBROUTINE form_modalexp(coef)
    USE global,ONLY:ndof,hrmbasis,maxbasis,debug
    USE integral,ONLY:harmint
    IMPLICIT NONE
    INTEGER ::mode
    REAL*8::coef(nDOF,maxbasis)
    IF(debug)PRINT*,"form_modalexp"
    modalexp=0.d0
    DO mode=1,ndof
      CALL update_modalexp(mode,coef)
    ENDDO!mode
    IF(debug)PRINT*,"form_modalexp :)"
  END SUBROUTINE form_modalexp

  SUBROUTINE update_modalexp(mode,coef)
    USE global,ONLY:ndof,hrmbasis,debug
    USE integral,ONLY:harmint
    IMPLICIT NONE
    INTEGER ::mode,power,n1,n2,diff
    REAL*8::coef1,coef2,tmp
    REAL*8::coef(nDOF,hrmbasis(mode))

    IF(debug)PRINT*,"update_modalexp"

    DO power=1,4
      tmp=0.d0
      DO n1 = 0,hrmbasis(mode)-1
        coef1=coef(mode,n1+1)
        !IF(ABS(coef1) .LT. inttol)CYCLE
        tmp = tmp + coef1 * coef1 * harmint(n1, mode, power, 0)
        DO n2 = n1+1,hrmbasis(mode)-1
          diff=n2-n1
          IF(diff > power)CYCLE
          coef2=coef(mode,n2+1)
          !        IF((ABS(coef2) .LT. inttol) .OR. (diff .GT. maxFF) )CYCLE
          tmp = tmp + 2.d0 * coef1 * coef2 * harmint(n2, mode, power, diff)
        ENDDO!n2
      ENDDO!n1
      modalexp(mode,power)=tmp
    ENDDO!power
      ! print*,"modal integrals stored"
    IF(debug)PRINT*,"update_modalexp :)"
  END SUBROUTINE update_modalexp

  REAL*8 FUNCTION modalkinexp(mode,coef)
    ! 1D integral for a modal in HO basis set.
    USE global,ONLY:ndof,hrmbasis,freq
    IMPLICIT NONE
    INTEGER ::mode,n1,n2
    REAL*8::coef(nDOF,hrmbasis(mode)),coef1,coef2,tmp
    tmp = 0.d0;
    DO n1 = 0,hrmbasis(mode)-1
      coef1=coef(mode,n1+1)
      !IF(ABS(coef1) .LT. inttol)CYCLE
      tmp= tmp + coef1 * coef1 * 0.5d0 * freq(mode) * (n1 + 0.5d0)
      n2=n1+2
      IF(n2 > hrmbasis(mode)-1)CYCLE
      coef2=coef(mode,n2+1)
      !        IF((ABS(coef2) .LT. inttol )CYCLE
      tmp = tmp + 2.d0 * coef1 * coef2 * (-0.25d0) * freq(mode) * dSQRT((n2 - 1.d0)*n2)
    ENDDO
    modalkinexp=tmp
    tmp=0.d0
    RETURN
  END FUNCTION modalkinexp

  REAL*8 FUNCTION modalexp1MR(m1)
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1

    modalexp1MR = Gi(m1)    * modalexp(m1,1)&
    + Hii(m1)   * modalexp(m1,2)&
    + Ciii(m1)  * modalexp(m1,3)&
    + Qiiii(m1) * modalexp(m1,4)
    RETURN
  END FUNCTION modalexp1MR

  REAL*8 FUNCTION modalexp2MR(m1,m2)
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2

    modalexp2MR =  modalexp(m1, 1) * modalexp(m2, 1) * Hij(m1,m2)&
    +  modalexp(m1, 2) * modalexp(m2, 2) * Qiijj(m1,m2)&
    +  modalexp(m1, 2) * modalexp(m2, 1) * Ciij(m1,m2)&
    +  modalexp(m1, 1) * modalexp(m2, 2) * Ciij(m2,m1)&
    +  modalexp(m1, 3) * modalexp(m2, 1) * Qiiij(m1,m2)&
    +  modalexp(m1, 1) * modalexp(m2, 3) * Qiiij(m2,m1);
    RETURN
  END FUNCTION modalexp2MR

  REAL*8 FUNCTION modalexp3MR(m1,m2,m3)
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,m3

    modalexp3MR = modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 1) * Cijk(m1,m2,m3)&
    + modalexp(m1, 2) * modalexp(m2, 1) * modalexp(m3, 1) * Qiijk(m1,m2,m3)&
    + modalexp(m1, 1) * modalexp(m2, 2) * modalexp(m3, 1) * Qiijk(m2,m1,m3)&
    + modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 2) * Qiijk(m3,m1,m2)
    RETURN
  END FUNCTION modalexp3MR

  REAL*8 FUNCTION modalexp4MR(m1,m2,m3,m4)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER, INTENT(in)::m1,m2,m3,m4

    modalexp4MR = modalexp(m1, 1) * modalexp(m2, 1) &
    * modalexp(m3, 1) * modalexp(m4, 1) * Qijkl(m1,m2,m3,m4)
    RETURN
  END FUNCTION modalexp4MR

  SUBROUTINE subFormUzero(m,coef)
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER,INTENT(in)::m
    REAL*8,INTENT(in)::coef(nDOF,hrmbasis(m))
    INTEGER::m1,m2,m3,m4
    REAL*8::tmp0
    if(debug) print*,"subformU for mode ",m
    tmp0=0.d0
    DO m1=1,ndof
      IF(m1==m)CYCLE
      tmp0=tmp0+modalkinexp(m1,coef)&!Kinetic term
      +Gi(m1)*modalexp(m1,1)&
      +Hii(m1)*modalexp(m1,2)&
      +Ciii(m1)*modalexp(m1,3)&
      +Qiiii(m1)*modalexp(m1,4)
      IF(nMR>1)THEN
        DO m2=m1+1,ndof
          IF(m2==m)CYCLE
          tmp0=tmp0+ Hij(m1,m2)*modalexp(m1,1)*modalexp(m2,1)&
          + modalexp(m1, 2) * modalexp(m2, 2) * Qiijj(m1,m2)&
          + modalexp(m1, 2) * modalexp(m2, 1) * Ciij(m1,m2)&
          + modalexp(m1, 1) * modalexp(m2, 2) * Ciij(m2,m1)&
          + modalexp(m1, 3) * modalexp(m2, 1) * Qiiij(m1,m2)&
          + modalexp(m1, 1) * modalexp(m2, 3) * Qiiij(m2,m1);
          IF(nMR>2)THEN
            DO m3=m2+1,ndof
              IF (m3 .EQ. m) CYCLE
              tmp0 =tmp0+ modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 1) * Cijk(m1,m2,m3)&
              +modalexp(m1, 2) * modalexp(m2, 1) * modalexp(m3, 1) * Qiijk(m1,m2,m3)&
              +modalexp(m1, 1) * modalexp(m2, 2) * modalexp(m3, 1) * Qiijk(m2,m1,m3)&
              +modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 2) * Qiijk(m3,m1,m2)
              IF(nMR>3)THEN
                DO m4=m3+1,ndof
                  IF (m4 .EQ. m) CYCLE
                  tmp0 =tmp0+ modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 1)* &
                  modalexp(m4, 1) * Qijkl(m1,m2,m3,m4)
                ENDDO !m4
              ENDIF!nmr>3
            ENDDO !m3
          ENDIF!nmr>2
        ENDDO !m2
      ENDIF!nmr>1
    ENDDO !m1
    uZero(m)=tmp0
    if(debug) print*,"subformU :)"
  END SUBROUTINE SubFormUzero

  REAL*8 FUNCTION funUzero(coef)
    !Uzero is caculated for m=nDOF
    USE global,ONLY:ndof,hrmbasis,nMr
    USE modQFF
    IMPLICIT NONE
    REAL*8,INTENT(in)::coef(nDOF,hrmbasis(ndof))
    INTEGER::m1,m2,m3,m4
    REAL*8::tmp0

    tmp0=0.d0
    DO m1=1,ndof-1
      tmp0=tmp0+modalkinexp(m1,coef)&!Kinetic term
      +Gi(m1)*modalexp(m1,1)&
      +Hii(m1)*modalexp(m1,2)&
      +Ciii(m1)*modalexp(m1,3)&
      +Qiiii(m1)*modalexp(m1,4)
      IF(nMR>1)THEN
        DO m2=m1+1,ndof-1
          tmp0=tmp0+Hij(m1,m2)*modalexp(m1,1)*modalexp(m2,1)&
          + modalexp(m1, 2) * modalexp(m2, 2) * Qiijj(m1,m2)&
          + modalexp(m1, 2) * modalexp(m2, 1) * Ciij(m1,m2)&
          + modalexp(m1, 1) * modalexp(m2, 2) * Ciij(m2,m1)&
          + modalexp(m1, 3) * modalexp(m2, 1) * Qiiij(m1,m2)&
          + modalexp(m1, 1) * modalexp(m2, 3) * Qiiij(m2,m1);
          IF(nMR>2)THEN
            DO m3=m2+1,ndof-1
              tmp0 =tmp0+ modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 1) * Cijk(m1,m2,m3)&
              +modalexp(m1, 2) * modalexp(m2, 1) * modalexp(m3, 1) * Qiijk(m1,m2,m3)&
              +modalexp(m1, 1) * modalexp(m2, 2) * modalexp(m3, 1) * Qiijk(m2,m1,m3)&
              +modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 2) * Qiijk(m3,m1,m2)
              IF(nMR>3)THEN
                DO m4=m3+1,ndof-1
                  tmp0 =tmp0+ modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 1)* &
                  modalexp(m4, 1) * Qijkl(m1,m2,m3,m4)
                ENDDO !m4
              ENDIF!nmr>3
            ENDDO !m3
          ENDIF!nmr>2
        ENDDO !m2
      ENDIF!nmr>1
    ENDDO !m1
    funUzero=tmp0
  END FUNCTION funUzero


  SUBROUTINE subFormU(m,coef)
    USE global!,only:ndof,maxbasis,nmr
    USE modQFF
    IMPLICIT NONE
    INTEGER,INTENT(in)::m
    REAL*8,INTENT(in)::coef(nDOF,hrmbasis(m))
    INTEGER::m1,m2,m3
    REAL*8::tmp1,tmp2,tmp3,tmp4
    IF(debug)PRINT*,"subFormU"
    tmp1=0.d0;tmp2=0.d0;tmp3=0.d0;tmp4=0.d0

    tmp1=tmp1+Gi(m)
    tmp2=tmp2+Hii(m)
    tmp3=tmp3+Ciii(m)
    tmp4=tmp4+Qiiii(m)
    IF(nMR>1)THEN
      DO m1=1,ndof
        IF(m1==m)CYCLE
        tmp1=tmp1+Hij(m,m1)*modalexp(m1,1)&
        +Ciij(m1,m)*modalexp(m1,2)&
        +Qiiij(m1,m)*modalexp(m1,3)
        tmp2=tmp2+Ciij(m,m1)*modalexp(m1,1)&
        +Qiijj(m,m1)*modalexp(m1,2)
        tmp3=tmp3+Qiiij(m,m1)*modalexp(m1,1)
        IF(nMR>2)THEN
          DO m2=m1+1,ndof
            IF(m2==m)CYCLE
            tmp1 =tmp1+ modalexp(m1, 1) * modalexp(m2, 1) * Cijk(m,m1,m2)&
            + modalexp(m1, 2) * modalexp(m2, 1) * Qiijk(m1,m,m2)&
            + modalexp(m1, 1) * modalexp(m2, 2) * Qiijk(m2,m,m1);
            tmp2 =tmp2+ modalexp(m1, 1) * modalexp(m2, 1) * Qiijk(m,m1,m2);
            IF(nMR>3)THEN
              DO m3=m2+1,ndof
                IF (m3 .EQ. m) CYCLE
                tmp1 =tmp1+ modalexp(m1, 1) * modalexp(m2, 1)* modalexp(m3, 1) * Qijkl(m3,m2,m1,m)
              ENDDO !m3
            ENDIF!nmr>3
          ENDDO !m2
        ENDIF!nmr>2
      ENDDO !m1
    ENDIF!nmr>1
    uOne(m)=tmp1
    uTwo(m)=tmp2!*0.5d0
    uThree(m)=tmp3!/6.0d0
    uFour(m)=tmp4!/24.0d0
    IF(debug)PRINT*,"subFormU :)"
  END SUBROUTINE SubFormU

  REAL*8 FUNCTION fungetE(thestate,givencoef,tmpnMR)  !TODO smt wrong for > 1MR and for all fund statesexcept the first one
    USE global!,only:ndof,maxbasis
    IMPLICIT NONE
    REAL*8,OPTIONAL::givencoef(nDOF,maxbasis)
    Integer,optional::tmpnMR
    REAL*8::coef(nDOF,maxbasis)
    INTEGER::m,m1,m2,m3,m4,originalnMR
    INTEGER,INTENT(in)::thestate(ndof)
    REAL*8::tmp
    ! print*,"subformU for mode ",m
    originalnMR=nMR
    if(present(tmpnMR))nMR=tmpnMR
    IF(PRESENT(givencoef))THEN
      coef=givencoef
    ELSE
      DO m=1,ndof
        coef(m,:)=vscfcoefs(m,:,thestate(m)+1)
      ENDDO
    ENDIF
    tmp=0.d0
    DO m1=1,ndof
      CALL update_modalexp(m1,coef)
      tmp = tmp + modalkinexp(m1,coef) &  !Kinetic term
      + modalexp1MR(m1)
      IF(nMR>1)THEN
        DO m2=m1+1,ndof
          tmp=tmp+modalexp2MR(m1,m2)
          IF(nMR>2)THEN
            DO m3=m2+1,ndof
              tmp =tmp+ modalexp3MR(m1,m2,m3)
              IF(nMR>3)THEN
                DO m4=m3+1,ndof
                  tmp =tmp+ modalexp4MR(m1,m2,m3,m4)
                ENDDO !m4
              ENDIF!nmr>3
            ENDDO !m3
          ENDIF!nmr>2
        ENDDO !m2
      ENDIF!nmr>1
    ENDDO !m1
    fungetE=tmp
!    print*,tmp
    if(present(tmpnMR))nMR=originalnMR
    RETURN
  END FUNCTION funGetE

  REAL*8 FUNCTION oldfungetE(thestate,givencoef)  !TODO smt wrong for > 1MR and for all fund statesexcept the first one
    USE global!,only:ndof,maxbasis
    IMPLICIT NONE
    REAL*8,OPTIONAL::givencoef(nDOF,maxbasis)
    REAL*8::coef(nDOF,maxbasis)
    INTEGER::m,m1,m2,m3,m4
    INTEGER,INTENT(in)::thestate(ndof)
    REAL*8::tmp
    !     print*,"subformU for mode ",m
    IF(PRESENT(givencoef))THEN
      coef=givencoef
    ELSE
      DO m=1,ndof
        coef(m,:)=vscfcoefs(m,:,thestate(m)+1)
      ENDDO
    ENDIF
    tmp=0.d0
    DO m1=1,ndof
      CALL update_modalexp(m1,coef)
      tmp = tmp + modalkinexp(m1,coef) &  !Kinetic term
      + modalexp1MR(m1)
    !      IF(nMR>1)THEN
    !        DO m2=m1+1,ndof
    !          tmp=tmp+modalexp2MR(m1,m2)
    !          IF(nMR>2)THEN
    !            DO m3=m2+1,ndof
    !              tmp =tmp+ modalexp3MR(m1,m2,m3)
    !              IF(nMR>3)THEN
    !                DO m4=m3+1,ndof
    !                  tmp =tmp+ modalexp4MR(m1,m2,m3,m4)
    !                ENDDO !m4
    !              ENDIF!nmr>3
    !            ENDDO !m3
    !          ENDIF!nmr>2
    !        ENDDO !m2
    !      ENDIF!nmr>1
    ENDDO !m1
    oldfungetE=tmp
    RETURN
  END FUNCTION oldfunGetE


  subroutine subget1MREnergy(coef,energy)
    USE global!,only:ndof,maxbasis
    IMPLICIT NONE
    REAL*8::coef(nDOF,maxbasis)
    INTEGER::m,m1
    REAL*8::tmp,energy

    IF(debug)PRINT*,"subget1MREnergy"
    tmp=0.d0
    DO m1=1,ndof
      CALL update_modalexp(m1,coef)
      tmp = tmp + modalkinexp(m1,coef) &  !Kinetic term
      + modalexp1MR(m1)
    ENDDO !m1
    energy=tmp
    print*,'1MR energy=',energy*convert
    IF(debug)PRINT*,"subget1MREnergy :)"
    RETURN
  end subroutine subget1MREnergy

  SUBROUTINE subFormH(mode,coef,Hvscf)
    USE global!,only:ndof,maxbasis,freq
    USE integral,ONLY:harmint
    IMPLICIT NONE
    INTEGER::mode,i
    REAL*8::Hvscf(hrmbasis(mode),hrmbasis(mode)),coef(nDOF,hrmbasis(mode))
    IF(debug)PRINT*,"subFormH"
    CALL subformU(mode,coef)
    Hvscf=0.d0
    DO i=0, hrmbasis(mode)-1
      Hvscf(i+1,i+1)= (i + 0.5) * freq(mode) * 0.5 &!kinetic term
                          !+Uzero(mode)& !Uzero can be added either here or after diag
      +Utwo(mode)* harmint(i, mode, 2, 0) &
      +Ufour(mode)* harmint(i, mode, 4, 0)
      IF (i .LT. hrmbasis(mode)-1) THEN
        Hvscf(i+1,i+2)=Uone(mode)* harmint(i+1, mode, 1, 1)+Uthree(mode)* harmint(i+1, mode, 3, 1)
        IF (i .LT. hrmbasis(mode)-2) THEN
          Hvscf(i+1,i+3)= -dSQRT(DBLE((i+2) * (i + 1))) * hrmfreq(mode) * 0.25d0 & !kinetic term
          +Utwo(mode)* harmint(i+2, mode, 2, 2)+Ufour(mode)* harmint(i+2, mode, 4, 2)
          IF (i .LT. hrmbasis(mode)-3) THEN
            Hvscf(i+1,i+4)=Uthree(mode)* harmint(i+3, mode, 3, 3)
            IF (i .LT. hrmbasis(mode)-4) THEN
              Hvscf(i+1,i+5)=Ufour(mode)* harmint(i+4, mode, 4, 4)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO !i
    IF(debug)PRINT*,"subFormH :)"
  END SUBROUTINE SubFormH

  SUBROUTINE subvscf(state,vscfenergy,vscfiter)
    USE global
    USE lapack
    USE yazar
    USE constants
    USE modtimer
    IMPLICIT NONE
    REAL*8:: ham(maxbasis,maxbasis),coef(nDOF,maxbasis),newcoef(nDOF,maxbasis)
    INTEGER:: mode,iter,state(nDOF)
    REAL*8, OPTIONAL :: vscfenergy
    INTEGER,OPTIONAL::vscfiter
    REAL*8:: energy,oldenergy,energies(maxbasis),modalE(nDOF),oldmodalE(nDOF)
    REAL(8)::wallstart,cpustart,wallend,cpuend
    IF(debug)PRINT*,"subvscf"
    iter=0
    oldenergy=0.d0
    coef=0.d0
    newcoef=0.d0
    freq=hrmfreq

    !Initializing using HO solutions according to target state.
    DO mode=1,nDOF
      oldmodalE(mode)=hrmfreq(mode)* (DBLE(state(mode)) + 0.5d0)
      oldenergy=oldenergy+oldmodalE(mode)
      coef(mode,state(mode)+1)=1.d0
      newcoef(mode,state(mode)+1)=1.d0
    ENDDO
    !start of iterations
    CALL form_modalexp(coef)
    IF(debug)PRINT*,"VSCF iterations start"
    DO
      !   PRINT*,'**************************  iterations **************************',iter
      iter=iter+1

      IF(benchmark)CALL wall_and_cpu_time(wallstart,cpustart)

      DO mode=1,nDOF
        IF (storepot)THEN
          CALL subformH(mode,coef,ham)
          CALL diag(hrmbasis(mode),ham,energies)
        ELSE
          CALL meanfieldhamiltonian(mode,coef,ham)
          CALL diag(hrmbasis(mode),ham,energies)
        ENDIF
        modalE(mode)=energies(state(mode)+1)
        coef(mode,:)=ham(:,state(mode)+1)
        CALL update_modalexp(mode,coef)
      ENDDO !mode

      IF(benchmark) THEN
        CALL wall_and_cpu_time(wallend,cpuend)
        PRINT*,"wall=",wallend-wallstart,"cpu=",cpuend-cpustart
      ENDIF

      IF ( (SUM(ABS(modalE-oldmodalE))*convert) .LT. scfthresh)THEN
        DO mode=1,nDOF
          IF (storepot)THEN
            CALL subformH(mode,coef,ham)
            CALL diag(hrmbasis(mode),ham,energies)
            IF(getUzero)CALL subformUzero(mode,coef)
            IF (mode==nDOF) THEN
              !CALL subformUzero(ndof,coef) !alternative 2
              !energy=energies(state(mode)+1)+ Uzero(mode)!alternative 2
              energy=energies(state(mode)+1)+funUzero(coef) !alternative 1
            ENDIF
          ELSE
            CALL meanfieldhamiltonian(mode,coef,ham)
            CALL diag(hrmbasis(mode),ham,energies)
            energy=energies(state(mode)+1)
          ENDIF
          IF(VMP2 .OR. virtual .or. VSCFCI .or. VSCFCI1)THEN
            vscfenergies(mode,:)=energies
            vscfcoefs(mode,:,:)=ham
          ENDIF
          modalE(mode)=energies(state(mode)+1)
          coef(mode,:)=ham(:,state(mode)+1)
          CALL update_modalexp(mode,coef)
        ENDDO !mode

        IF ((ABS(energy-oldenergy)*convert .LT. scfthresh) )EXIT
        iter=iter+1
      ENDIF

      IF(iter .GT. maxiter)THEN
        PRINT*,"maxiter reached"
        EXIT
      ENDIF

      oldenergy=energy
      oldmodalE=modalE

    ENDDO ! Do While

    IF (PRESENT(vscfenergy))vscfenergy=energy
    IF (PRESENT(vscfiter))vscfiter=iter
    IF(debug)PRINT*,"VSCF converged after ",iter ," iterations"
    IF(getwfn)CALL subscheckstate(state,coef)
    IF(getwfn.AND.debug)CALL writemat(coef,ndof,maxbasis)

    IF(virtual)THEN
      PRINT*,"vscfenergies"
      CALL writemat(vscfenergies*convert,ndof,maxbasis)
      CALL printvirtual()
    ENDIF

    ! PRINT*,fungetE(coef)*convert
    ! print*,'vscf done for', state
    IF(get1MRenergy.and. sum(state)==0)then
      Call subget1MRenergy(coef,energy)
      vscfenergy=energy
    endif

    !    IF(debug)print*,'checkfungetE',oldfungetE(state,coef)*convert
    if(printU)call subprintU()
    IF(debug)PRINT*,"subvscf :)"
  END SUBROUTINE subvscf

  SUBROUTINE subvscfvirtual()
    USE global
    USE lapack
    USE yazar
    USE constants
    USE modtimer
    IMPLICIT NONE
    REAL*8:: ham(maxbasis,maxbasis),energies(maxbasis),coef(nDOF,maxbasis)
    INTEGER:: mode,iter,state(nDOF)
    REAL*8:: energy,oldenergy,newcoef(nDOF,maxbasis),modalE(nDOF),oldmodalE(nDOF),Egs,Exc
    REAL(8)::wallstart,cpustart,wallend,cpuend
    IF(debug)PRINT*,"subvscfvirtual"
    iter=0
    oldenergy=1000
    energy=0.d0
    coef=0.d0
    newcoef=0.d0
    freq=hrmfreq
    !Forming initial wavefunction using HO basis according to target state.
    state=0

    DO mode=1,nDOF
      oldmodalE(mode)=freq(mode)* (DBLE(state(mode)) + 0.5d0)
      coef(mode,state(mode)+1)=1.d0
      newcoef(mode,state(mode)+1)=1.d0
    ENDDO
    CALL form_modalexp(coef)
    IF(debug)PRINT*,"VSCF iterations start "
    DO
      oldenergy=energy
      energy=0.d0
      !  PRINT*,'**************************  iteration **************************',iter
      iter=iter+1
      IF(benchmark)CALL wall_and_cpu_time(wallstart,cpustart)
      DO mode=1,nDOF
        IF (storepot)THEN
          CALL subformH(mode,coef,ham)
          CALL diag(hrmbasis(mode),ham,energies)
          IF (mode==nDOF) THEN
            !CALL subformUzero(ndof,coef) !alternative 2
            !energy=energies(state(mode)+1)+ Uzero(mode)!alternative 2
            energy=energies(state(mode)+1)+funUzero(coef) !alternative 1
          ENDIF
        ELSE
          CALL meanfieldhamiltonian(mode,coef,ham)
          CALL diag(hrmbasis(mode),ham,energies)
          energy=energies(state(mode)+1)
        ENDIF
        ! newcoef(mode,:)=ham(:,state(mode)+1)
        modalE(mode)=energies(state(mode)+1)
        coef(mode,:)=ham(:,state(mode)+1)
        CALL update_modalexp(mode,coef)
         ! PRINT*,'eigvec'!ham(1,:)
      ENDDO !mode
      IF(benchmark) THEN
        CALL wall_and_cpu_time(wallend,cpuend)
        PRINT*,"wall=",wallend-wallstart,"cpu=",cpuend-cpustart
      ENDIF
      ! coef=newcoef
      IF(iter .GT. maxiter)THEN
        PRINT*,"maxiter reached"
        EXIT
      ENDIF
      !IF ((ABS(energy-oldenergy) .LT. scfthresh))THEN
      IF ( (SUM(ABS(modalE-oldmodalE))*convert) .LT. scfthresh)THEN
        DO mode=1,nDOF
          IF (storepot)THEN
            CALL subformH(mode,coef,ham)
            CALL diag(hrmbasis(mode),ham,energies)
            IF (mode==nDOF) THEN
              !CALL subformUzero(ndof,coef) !alternative 2
              !energy=energies(state(mode)+1)+ Uzero(mode)!alternative 2
              energy=energies(state(mode)+1)+funUzero(coef) !alternative 1
            ENDIF
          ELSE
            CALL meanfieldhamiltonian(mode,coef,ham)
            CALL diag(hrmbasis(mode),ham,energies)
            energy=energies(state(mode)+1)
          ENDIF
          vscfenergies(mode,:)=energies
          IF(vmp2 .OR. virtual .or. VSCFCI .or. VSCFCI1)   vscfcoefs(mode,:,:)=ham

          ! newcoef(mode,:)=ham(:,state(mode)+1)
          coef(mode,:)=ham(:,state(mode)+1)
          CALL update_modalexp(mode,coef)
           ! PRINT*,'eigvec'!ham(1,:)
        ENDDO !mode
        IF ((ABS(energy-oldenergy)*convert .LT. scfthresh) )EXIT
        iter=iter+1
      ENDIF
      IF(iter .GT. maxiter)THEN
        PRINT*,"maxiter reached"
        EXIT
      ENDIF
      oldenergy=energy
      oldmodalE=modalE
    ENDDO ! Do While
    energy=energies(1)+funUzero(coef)

    PRINT*,'vVSCF'
    Egs=energy*convert
    PRINT*,Egs
    IF(getwfn)CALL subscheckstate(state,coef)
    IF (maxbasis > 1)THEN
      DO mode=1,ndof
        Exc=(vscfenergies(mode,2)-vscfenergies(mode,1))*convert
        PRINT*,Exc!,Exc+Egs!,(Energy-theenergy)*convert
         !   print*, modalE(mode)*convert,(energy+modalE(mode))*convert
      ENDDO
    ENDIF
    IF(debug)PRINT*,"VSCF iterations= ",iter
    IF(getwfn.AND.debug)CALL writemat(coef,ndof,maxbasis)
    IF(plotU)call subplotU()
    ! if(virtual)call printvirtual()
    !    PRINT*,fungetE(coef)*convert
    ! print*,'vscf done for', state
    IF(debug)PRINT*,"subvscfvirtual :)"
    RETURN
  END SUBROUTINE subvscfvirtual


  SUBROUTINE subscheckstate(state,coef)
    USE global
    USE modQFF
    USE integral
    IMPLICIT NONE
    INTEGER,INTENT(in)::state(nDOF)
    REAL*8,INTENT(in):: coef(nDOF,maxbasis)
    REAL*8:: sumcoef,prob(maxbasis)
    INTEGER::i,m,stateout(nDOF)

    IF(debug)PRINT*,"subscheckstate"
    !print*,maxloc(abs(coef),2),"kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk"
    stateout=MAXLOC(ABS(coef),2)-1
    !test
    ! print*,"Converged state is", stateout , "for initial HO state", state

    IF( COUNT((stateout-state)==0) .NE. ndof ) THEN
      rightstate=.FALSE.
      PRINT*,"Converged state is", stateout !, "for initial HO state", state
    ENDIF
    DO m=1,ndof
      !weird if( (maxloc(coef(m,:))) .ne. (state(m)+1) ) rightstate=.false. gives error
      !says if clause is not scalar.
      prob=SQRT(coef(m,:)*coef(m,:))
      IF (MAXVAL(prob) < 0.50) THEN
        sumcoef=0.d0
        PRINT*,"resonance for mode ",m, "for state", stateout
        DO i=0,hrmbasis(m)-1
          PRINT*,prob(i+1)
          sumcoef=sumcoef+prob(i+1)
          IF(sumcoef > 0.95)EXIT
        ENDDO
      ENDIF
    ENDDO !m
    IF(debug)PRINT*,"subscheckstate :)"
  END SUBROUTINE subscheckstate

  SUBROUTINE printvirtual
    !VSCF for given state and calculate excited states
    USE global
    USE modQFF
    USE integral
    IMPLICIT NONE
    REAL*8 :: Egs,Exc
    INTEGER::m,state(nDOF)
    state=0
    PRINT*,"Virtual VSCF fund energies"
    Egs=fungetE(state)
    PRINT*,Egs*convert
    ! PRINT*,groundE*convert
    DO m=1,ndof
      state=0
      state(m)=1
      Exc=fungetE(state)
      PRINT*,(vscfenergies(m,2)-vscfenergies(m,1))*convert,(Exc-Egs)*convert, Exc*convert
       ! CALL getenergy(state,energy)
       ! PRINT*,(energy-groundE)*convert,energy*convert
       ! PRINT*,"fungetE ",fungetE(state)*convert
    ENDDO
    state(1)=2
    state(2)=1
    state(3)=0
    Exc=fungetE(state)
    PRINT*,(Exc-Egs)*convert,Exc*convert,"for state ",state
    PRINT*,"************************"
  END SUBROUTINE printvirtual

  subroutine subplotU
    use global
    use modQFF
    implicit none
    integer::m
    CHARACTER :: plotfile*20,epsfile*20,cndof*5
    !  OPEN(UNIT=19,FILE='U.mavi',STATUS='UNKNOWN',ACTION='WRITE')
    if(debug)print*,'subplotU'
    WRITE(cndof,'(i5)')ndof
    plotfile='U_N'//TRIM(ADJUSTL(cndof))//'.plot'
    epsfile='U_N'//TRIM(ADJUSTL(cndof))//'.eps'
    OPEN(21,FILE=TRIM(plotfile),STATUS='UNKNOWN',ACTION='WRITE')
    WRITE(21,*)'set terminal postscript enhanced eps "Helvetica" 14'
    WRITE(21,*)'set output "',trim(epsfile),'"'
    WRITE(21,*)'set size 1.0, ',1!ndof
    WRITE(21,*)'set origin 0.0,0.0'
    WRITE(21,*)'set xrange[-10:10]'
    WRITE(21,*)'set multiplot'

    !do
    m=1!,1!ndof
    WRITE(21,'(a,i5,a)')'set title "mode ', m, '"'
    WRITE(21,*)'set title offset 0,-10'
    WRITE(21,*)'set size 1.0,1.0 '
    WRITE(21,*)'set dummy Q'
    WRITE(21,*)'set origin 0.0, ', m-1
    WRITE(21,210)'plot ', Uone(m), ' *Q+ ',Utwo(m), ' *Q**2+ ',Uthree(m), ' *Q**3+ ',Ufour(m), ' *Q**4 w l lt 1 lw 3\'
    WRITE(21,210)',', Gi(m), ' *Q+ ',Hii(m), ' *Q**2+ ',Ciii(m), ' *Q**3+ ',Qiiii(m), ' *Q**4 w l lt 3 lw 3'
    !    write(19,190)m,Uone(m),Utwo(m),Uthree(m),Ufour(m)
    !    write(19,190)m,Gi(m),Hii(m),Ciii(m),Qiiii(m)
    !enddo !m
    WRITE(21,*)'unset multiplot'
    !190 FORMAT(i4,4f20.8)
    !210 FORMAT(a,f12.4, a,f12.4, a,f12.4, a,f12.4, a)
210 FORMAT(a,e12.4, a,e12.4, a,e12.4, a,e12.4, a)
    !210 FORMAT(a,e, a,e, a,e, a,e, a)
    !CLOSE(19)
    CLOSE(21)
    if(debug)print*,'subplotU :)'
  end subroutine subplotU

  subroutine subprintU
    use global
    use modQFF
    implicit none
    integer::m
    CHARACTER :: datafile*20,cndof*5
    !  OPEN(UNIT=19,FILE='U.mavi',STATUS='UNKNOWN',ACTION='WRITE')
    if(debug)print*,'subprintU'
    WRITE(cndof,'(i5)')ndof
    datafile='U_N'//TRIM(ADJUSTL(cndof))//'.dat'
    OPEN(21,FILE=TRIM(datafile),STATUS='UNKNOWN',ACTION='WRITE')
    call subprintUE()
    WRITE(21,*)"U"
    do m=1,ndof
      !      write(21,*) m
      WRITE(21,100)Uzero(m),Uone(m),Utwo(m),Uthree(m),Ufour(m)
    enddo !m
100 FORMAT(5ES20.8)
    CLOSE(21)
    if(debug)print*,'subprintU :)'
  end subroutine subprintU

    subroutine subprintUE
    use global
    USE yazar
    IMPLICIT NONE
    integer::m
    real*8::Ekin,EU0,EU1,EU2,EU3,EU4,Eall
    if (debug)print*,'subUenergy'
    WRITE(21,*)"EU"
    do m=1,ndof
      Ekin=modalkinexp(m,vscfcoefs(:,:,1))
      call  subFormUzero(m,vscfcoefs(:,:,1))
      EU0=Uzero(m)
      EU1=Uone(m)*modalexp(m,1)
      EU2=Utwo(m)*modalexp(m,2)
      EU3=Uthree(m)*modalexp(m,3)
      EU4=Ufour(m)*modalexp(m,4)
      Eall=Ekin+EU0+EU1+EU2+EU3+EU4
      WRITE(21,100)EU0,EU1,EU2,EU3,EU4
100 FORMAT(5ES20.8)
!      print*,m,Eall,EU0,EU1,EU2,EU3,EU4
    !     print*,EU2,EU3,EU4
    enddo
    if (debug)print*,'subprintUE :)'
  end subroutine subprintUE

  SUBROUTINE printssVSCF
    !VSCF for given state and calculate excited states
    USE global
    USE modQFF
    USE integral
    IMPLICIT NONE
    REAL*8 :: Egs,Exc
    INTEGER::state(nDOF),iter

    state(1)=2
    state(2)=1
    state(3)=0
    CALL subvscf(state,Exc,iter)
    PRINT*,(Exc)*convert,"for SS vscf state ",state, "iter", iter
    Exc=fungetE(state)
    state=0
    Egs=fungetE(state)
    PRINT*,(Exc-Egs)*convert,Exc*convert,"for state ",state
    PRINT*,"************************"
  END SUBROUTINE printssVSCF

  REAL*8 FUNCTION meanfieldkinetic(mode,coef)
    USE global
    USE modQFF
    USE integral
    IMPLICIT NONE
    INTEGER::i
    REAL*8 :: tmp
    INTEGER,INTENT(in)::mode
    REAL*8,INTENT(in)::coef(nDOF,hrmbasis(mode))


    tmp = 0.d0;
    !First method
    !  DO i=0, maxbasis(mode)-1
    !     DO j=0, maxbasis(mode)-1
    !        diff = j - i;
    !        IF (diff .EQ. 0) THEN
    !           tmp =tmp+ coef(mode,i+1) * coef(mode,i+1) * (dble(i) + 0.5d0) * hrmfreq(mode) * 0.5d0;
    !        ELSE IF (ABS(diff) .EQ. 2) THEN
    !          ! tmp =tmp- coef(mode,i+1) * coef(mode,j+1) * dSQRT(DBLE( MAX(i, j) * (MAX(i, j) - 1) )) * hrmfreq(mode) * 0.25d0;
    !        ENDIF
    !     ENDDO!j
    !  ENDDO!i

    !second method
    DO i=0, hrmbasis(mode)-1
      tmp =tmp+ coef(mode,i+1) * coef(mode,i+1) * (DBLE(i) + 0.5d0) * hrmfreq(mode) * 0.5d0
      IF(i.LT.(hrmbasis(mode)-2))&
      tmp =tmp- coef(mode,i+1) * coef(mode,i+3) * 0.5d0 * hrmfreq(mode)* dSQRT(DBLE((i+2)*(i+1)));
    ENDDO!i
    !end of 2nd
    !3rd way
    ! tmp=modalexp(mode,-2,coef)
    meanfieldkinetic=tmp
    ! print*, tmp*convert,"meanfield kinetic for" ,mode
    RETURN
  END FUNCTION meanfieldkinetic
  REAL*8 FUNCTION meanfieldpotential(mode,coef)
    USE constants
    USE global
    USE modQFF
    IMPLICIT NONE
    INTEGER::m,m1,m2,m3
    REAL*8 :: tmp!,meanFieldKinetic
    INTEGER,INTENT(in)::mode
    REAL*8,INTENT(in)::coef(nDOF,hrmbasis(mode))
    tmp = 0.0;
    ! WRITE(*,'(d50.25)')c_gamma
    DO m=1,Ndof
      IF (m .EQ. mode)CYCLE
      tmp =tmp+ modalexp(m, 1) * gi(m);
      tmp =tmp+ modalexp(m, 2) * Hii(m);
      tmp =tmp+ modalexp(m, 3) * Ciii(m);
      tmp =tmp+ modalexp(m, 4) * Qiiii(m);
      tmp =tmp+ meanFieldKinetic(m,coef);
    ENDDO
    IF (nMR .GT. 1)THEN
      ! print*, "**************************nMR >1 "
      DO m1=1,nDOF
        IF (m1 .EQ. mode)CYCLE
        DO m2=m1+1,nDOF
          IF (m2 .EQ. mode .OR. m1 .EQ. m2)CYCLE
          tmp =tmp+ modalexp(m1, 1) * modalexp(m2, 1) * Hij(m1,m2);
          tmp =tmp+ modalexp(m1, 2) * modalexp(m2, 2) * Qiijj(m1,m2);
          tmp =tmp+ modalexp(m1, 2) * modalexp(m2, 1) * Ciij(m1,m2);
          tmp =tmp+ modalexp(m1, 1) * modalexp(m2, 2) * Ciij(m2,m1);
          tmp =tmp+ modalexp(m1, 3) * modalexp(m2, 1) * Qiiij(m1,m2);
          tmp =tmp+ modalexp(m1, 1) * modalexp(m2, 3) * Qiiij(m2,m1);
        ENDDO!m2
      ENDDO!m1
      IF (nMR .GT. 2)THEN
        DO m1=1,nDOF
          IF (m1 .EQ. mode) CYCLE
          DO m2=m1+1,nDOF
            IF (m2 .EQ. mode) CYCLE
            DO m3=m2+1,nDOF
              IF (m3 .EQ. mode) CYCLE
              tmp =tmp+ modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 1) * Cijk(m1,m2,m3);
              tmp =tmp+ modalexp(m1, 2) * modalexp(m2, 1) * modalexp(m3, 1) * Qiijk(m1,m2,m3);
              tmp =tmp+ modalexp(m1, 1) * modalexp(m2, 2) * modalexp(m3, 1) * Qiijk(m2,m1,m3);
              tmp =tmp+ modalexp(m1, 1) * modalexp(m2, 1) * modalexp(m3, 2) * Qiijk(m3,m1,m2);
            ENDDO!m3
          ENDDO!m2
        ENDDO!m1
      ENDIF!nMR>2
    ENDIF!nMR>1
    meanfieldpotential=tmp
    !  print*, tmp*convert,"meanfield pot for" ,mode
    RETURN
  END FUNCTION meanfieldpotential

  SUBROUTINE meanFieldHamiltonian(mode,coef,Hvscf)
    USE global
    USE modQFF
    USE integral,ONLY:harmint
    IMPLICIT NONE
    INTEGER::diff,i,j,m,m1,m2
    REAL*8::tmp,meanpotkin!,meanfieldpotential
    INTEGER,INTENT(in)::mode
    REAL*8,INTENT(in)::coef(nDOF,hrmbasis(mode))
    REAL*8::Hvscf(hrmbasis(mode),hrmbasis(mode))
    !PRINT*, "forming VSCF Hamiltonian for mode" , mode
    meanpotkin=meanFieldPotential(mode,coef)
    ! PRINT*,'meanpotkin',meanpotkin
    Hvscf=0.d0

    DO i=0, hrmbasis(mode)-1
      tmp=0.0
      tmp= tmp + (i + 0.5) * hrmfreq(mode) * 0.5;      !kinetic term
      tmp= tmp + Hii(mode) * harmint(i, mode, 2, 0);
      tmp= tmp + Qiiii(mode) * harmint(i, mode, 4, 0);
      tmp= tmp + meanpotkin
      IF (nMR .GT. 1) THEN
        DO m=1,ndof
          IF (m .EQ. mode)CYCLE
          tmp= tmp + Hij(m,mode) * harmint(i, mode, 1, 0) * modalexp(m, 1);
          tmp= tmp + Qiijj(m,mode) * harmint(i, mode, 2, 0) * modalexp(m, 2);
          tmp= tmp + Ciij(mode,m) * harmint(i, mode, 2, 0) * modalexp(m, 1);
          tmp= tmp + Ciij(m,mode) * harmint(i, mode, 1, 0) * modalexp(m, 2);
          tmp= tmp + Qiiij(mode,m) * harmint(i, mode, 3, 0) * modalexp(m, 1);
          tmp= tmp + Qiiij(m,mode) * harmint(i, mode, 1, 0) * modalexp(m, 3);
        ENDDO
        IF (nMR .GT. 2) THEN
          DO m1=1,ndof-2
            IF (m1 .EQ. mode)CYCLE
            DO m2=m1+1,ndof-1
              IF (m2 .EQ. mode) CYCLE
              tmp= tmp + Cijk(mode,m1,m2) * harmint(i, mode, 1, 0) * modalexp(m1, 1)  * modalexp(m2, 1);
              tmp= tmp + Qiijk(mode,m1,m2) * harmint(i, mode, 2, 0) * modalexp(m1, 1) * modalexp(m2, 1);
              tmp= tmp + Qiijk(m1,mode,m2) * harmint(i, mode, 1, 0) * modalexp(m1, 2)* modalexp(m2, 1);
              tmp= tmp + Qiijk(m2,mode,m1) * harmint(i, mode, 1, 0) * modalexp(m1, 1) * modalexp(m2, 2);
            ENDDO!!m2
          ENDDO!m1
        ENDIF!nMR>2
      ENDIF!nMR>1
      Hvscf(i+1,i+1)=tmp
      DO j=i+1,hrmbasis(mode)-1!
        tmp=0.0
        diff = j - i;
        IF (diff > 4) EXIT !break;
        IF (diff .EQ. 2) tmp= tmp-dSQRT(DBLE(j * (j - 1))) * hrmfreq(mode) * 0.25 !!kinetic term
        tmp= tmp + harmint(j, mode, 1, diff) * gi(mode);
        tmp= tmp + harmint(j, mode, 2, diff) * Hii(mode);
        tmp= tmp + harmint(j, mode, 3, diff) * Ciii(mode);
        tmp= tmp + harmint(j, mode, 4, diff) * Qiiii(mode);
        ! print*,'test2',tmp
        IF (nMR .GT. 1)THEN
          DO m=1,ndof
            IF (m .EQ. mode)CYCLE
            tmp= tmp + Hij(m,mode) * harmint(j, mode, 1, diff) * modalexp(m, 1);
            tmp= tmp + Qiijj(m,mode) * harmint(j, mode, 2, diff) * modalexp(m, 2);
            tmp= tmp + Ciij(mode,m) * harmint(j, mode, 2, diff) * modalexp(m, 1);
            tmp= tmp + Ciij(m,mode) * harmint(j, mode, 1, diff) * modalexp(m, 2);
            tmp= tmp + Qiiij(mode,m) * harmint(j, mode, 3, diff) * modalexp(m, 1);
            tmp= tmp + Qiiij(m,mode) * harmint(j, mode, 1, diff) * modalexp(m, 3);
          ENDDO
          IF (nMR .GT. 2) THEN
            DO m1=1,ndof
              IF (m1 .EQ. mode)CYCLE
              DO m2=m1 + 1,ndof
                IF (m2 .EQ. mode)CYCLE
                tmp= tmp + Cijk(mode,m1,m2)*harmint(j,mode,1,diff) * modalexp(m1,1)* modalexp(m2,1);
                tmp= tmp + Qiijk(mode,m1,m2)*harmint(j,mode,2,diff) * modalexp(m1,1)* modalexp(m2,1);
                tmp= tmp + Qiijk(m1,mode,m2)*harmint(j,mode,1,diff) * modalexp(m1,2)* modalexp(m2,1);
                tmp= tmp + Qiijk(m2,mode,m1)*harmint(j,mode,1,diff) * modalexp(m1,1)* modalexp(m2,2);
              ENDDO!!m2
            ENDDO!!m1
          ENDIF!nMR>2
        ENDIF!nMR>1
        ! print*,'test3',tmp
        Hvscf(i+1,j+1)=tmp
      ENDDO!!end for j
    ENDDO!!end for i
  END SUBROUTINE meanFieldHamiltonian
END MODULE modVSCF
