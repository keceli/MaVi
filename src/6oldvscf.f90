!MODULE modVSCFold
!  IMPLICIT NONE
!CONTAINS
!
!  SUBROUTINE subFormUzero(m,coef)
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    INTEGER,INTENT(in)::m
!    REAL*8,INTENT(in)::coef(nDOF,maxbasis(m))
!    INTEGER::m1,m2,m3,m4
!    REAL*8::tmp0
!    ! print*,"subformU for mode ",m
!    tmp0=0.d0
!    DO m1=1,ndof
!       IF(m1==m)CYCLE
!       tmp0=tmp0+modalexp(m1,-2,coef)&!Kinetic term
!            +Gi(m1)*modalexp(m1,1,coef)&
!            +Hii(m1)*modalexp(m1,2,coef)&
!            +Ciii(m1)*modalexp(m1,3,coef)&
!            +Qiiii(m1)*modalexp(m1,4,coef)
!       IF(nMR>1)THEN
!          DO m2=m1+1,ndof
!             IF(m2==m)CYCLE
!             tmp0=tmp0+ Hij(m1,m2)*modalexp(m1,1,coef)*modalexp(m2,1,coef)&
!                  + modalexp(m1, 2,coef) * modalexp(m2, 2,coef) * Qiijj(m1,m2)&
!                  + modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * Ciij(m1,m2)&
!                  + modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * Ciij(m2,m1)&
!                  + modalexp(m1, 3,coef) * modalexp(m2, 1,coef) * Qiiij(m1,m2)&
!                  + modalexp(m1, 1,coef) * modalexp(m2, 3,coef) * Qiiij(m2,m1);
!             IF(nMR>2)THEN
!                DO m3=m2+1,ndof
!                   IF (m3 .EQ. m) CYCLE
!                   tmp0 =tmp0+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef) * Cijk(m1,m2,m3)&
!                        +modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef) * Qiijk(m1,m2,m3)&
!                        +modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * modalexp(m3, 1,coef) * Qiijk(m2,m1,m3)&
!                        +modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 2,coef) * Qiijk(m3,m1,m2)
!                   IF(nMR>3)THEN
!                      DO m4=m3+1,ndof
!                         IF (m4 .EQ. m) CYCLE
!                         tmp0 =tmp0+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef)* &
!                              modalexp(m4, 1,coef) * Qijkl(m1,m2,m3,m4)
!                      ENDDO !m4
!                   ENDIF!nmr>3
!                ENDDO !m3
!             ENDIF!nmr>2
!          ENDDO !m2
!       ENDIF!nmr>1
!    ENDDO !m1
!    uZero(m)=tmp0
!  END SUBROUTINE SubFormUzero
!
!  REAL*8 FUNCTION funUzero(coef)
!    !Uzero is caculated for m=nDOF
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    REAL*8,INTENT(in)::coef(nDOF,maxbasis(ndof))
!    INTEGER::m1,m2,m3,m4
!    REAL*8::tmp0
!
!    tmp0=0.d0
!    DO m1=1,ndof-1
!       tmp0=tmp0+modalexp(m1,-2,coef)&!Kinetic term
!            +Gi(m1)*modalexp(m1,1,coef)&
!            +Hii(m1)*modalexp(m1,2,coef)&
!            +Ciii(m1)*modalexp(m1,3,coef)&
!            +Qiiii(m1)*modalexp(m1,4,coef)
!       IF(nMR>1)THEN
!          DO m2=m1+1,ndof-1
!             tmp0=tmp0+Hij(m1,m2)*modalexp(m1,1,coef)*modalexp(m2,1,coef)&
!                  + modalexp(m1, 2,coef) * modalexp(m2, 2,coef) * Qiijj(m1,m2)&
!                  + modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * Ciij(m1,m2)&
!                  + modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * Ciij(m2,m1)&
!                  + modalexp(m1, 3,coef) * modalexp(m2, 1,coef) * Qiiij(m1,m2)&
!                  + modalexp(m1, 1,coef) * modalexp(m2, 3,coef) * Qiiij(m2,m1);
!             IF(nMR>2)THEN
!                DO m3=m2+1,ndof-1
!                   tmp0 =tmp0+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef) * Cijk(m1,m2,m3)&
!                        +modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef) * Qiijk(m1,m2,m3)&
!                        +modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * modalexp(m3, 1,coef) * Qiijk(m2,m1,m3)&
!                        +modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 2,coef) * Qiijk(m3,m1,m2)
!                   IF(nMR>3)THEN
!                      DO m4=m3+1,ndof-1
!                         tmp0 =tmp0+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef)* &
!                              modalexp(m4, 1,coef) * Qijkl(m1,m2,m3,m4)
!                      ENDDO !m4
!                   ENDIF!nmr>3
!                ENDDO !m3
!             ENDIF!nmr>2
!          ENDDO !m2
!       ENDIF!nmr>1
!    ENDDO !m1
!    funUzero=tmp0
!  END FUNCTION funUzero
!
!    SUBROUTINE subFormU(m,coef)
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    INTEGER,INTENT(in)::m
!    REAL*8,INTENT(in)::coef(nDOF,maxbasis(m))
!    INTEGER::m1,m2,m3
!    REAL*8::tmp1,tmp2,tmp3,tmp4
!
!    tmp1=0.d0;tmp2=0.d0;tmp3=0.d0;tmp4=0.d0
!
!    tmp1=tmp1+Gi(m)
!    tmp2=tmp2+Hii(m)
!    tmp3=tmp3+Ciii(m)
!    tmp4=tmp4+Qiiii(m)
!    IF(nMR>1)THEN
!       DO m1=1,ndof
!          IF(m1==m)CYCLE
!          tmp1=tmp1+Hij(m,m1)*modalexp(m1,1,coef)&
!               +Ciij(m1,m)*modalexp(m1,2,coef)&
!               +Qiiij(m1,m)*modalexp(m1,3,coef)
!          tmp2=tmp2+Ciij(m,m1)*modalexp(m1,1,coef)&
!               +Qiijj(m,m1)*modalexp(m1,2,coef)
!          tmp3=tmp3+Qiiij(m,m1)*modalexp(m1,1,coef)
!          IF(nMR>2)THEN
!             DO m2=m1+1,ndof
!                IF(m2==m)CYCLE
!                tmp1 =tmp1+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * Cijk(m,m1,m2)&
!                     + modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * Qiijk(m1,m,m2)&
!                     + modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * Qiijk(m2,m,m1);
!                tmp2 =tmp2+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * Qiijk(m,m1,m2);
!                IF(nMR>3)THEN
!                   DO m3=m2+1,ndof
!                      IF (m3 .EQ. m) CYCLE
!                      tmp1 =tmp1+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef)* modalexp(m3, 1,coef) * Qijkl(m3,m2,m1,m)
!                   ENDDO !m3
!                ENDIF!nmr>3
!             ENDDO !m2
!          ENDIF!nmr>2
!       ENDDO !m1
!    ENDIF!nmr>1
!    uOne(m)=tmp1
!    uTwo(m)=tmp2!*0.5d0
!    uThree(m)=tmp3!/6.0d0
!    uFour(m)=tmp4!/24.0d0
!  END SUBROUTINE SubFormU
!
!  REAL*8 FUNCTION fungetE(thestate,givencoef)
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    REAL*8,OPTIONAL::givencoef(nDOF,MAXVAL(maxbasis))
!    REAL*8::coef(nDOF,MAXVAL(maxbasis))
!    INTEGER::m,m1,m2,m3,m4
!    INTEGER,INTENT(in)::thestate(ndof)
!    REAL*8::tmp
!    ! print*,"subformU for mode ",m
!    IF(PRESENT(givencoef))THEN
!       coef=givencoef
!    ELSE
!       DO m=1,ndof
!          coef(m,:)=vscfcoefs(m,:,thestate(m)+1)
!       ENDDO
!    ENDIF
!    tmp=0.d0
!    DO m1=1,ndof
!       tmp=tmp+modalexp(m1,-2,coef)&!Kinetic term
!            +Gi(m1)*modalexp(m1,1,coef)&
!            +Hii(m1)*modalexp(m1,2,coef)&
!            +Ciii(m1)*modalexp(m1,3,coef)&
!            +Qiiii(m1)*modalexp(m1,4,coef)
!       IF(nMR>1)THEN
!          DO m2=m1+1,ndof
!             tmp=tmp+Hij(m1,m2)*modalexp(m1,1,coef)*modalexp(m2,1,coef)&
!                  + modalexp(m1, 2,coef) * modalexp(m2, 2,coef) * Qiijj(m1,m2)&
!                  + modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * Ciij(m1,m2)&
!                  + modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * Ciij(m2,m1)&
!                  + modalexp(m1, 3,coef) * modalexp(m2, 1,coef) * Qiiij(m1,m2)&
!                  + modalexp(m1, 1,coef) * modalexp(m2, 3,coef) * Qiiij(m2,m1);
!             IF(nMR>2)THEN
!                DO m3=m2+1,ndof
!                   tmp =tmp+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef) * Cijk(m1,m2,m3)&
!                        + modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef) * Qiijk(m1,m2,m3)&
!                        + modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * modalexp(m3, 1,coef) * Qiijk(m2,m1,m3)&
!                        + modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 2,coef) * Qiijk(m3,m1,m2);
!                   IF(nMR>3)THEN
!                      DO m4=m3+1,ndof
!                         tmp =tmp+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef)* &
!                              modalexp(m4, 1,coef) * Qijkl(m1,m2,m3,m4)
!                      ENDDO !m4
!                   ENDIF!nmr>3
!                ENDDO !m3
!             ENDIF!nmr>2
!          ENDDO !m2
!       ENDIF!nmr>1
!    ENDDO !m1
!    fungetE=tmp
!    RETURN
!  END FUNCTION funGetE
!
!  SUBROUTINE subFormH(mode,coef,Hvscf)
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    INTEGER::mode,i
!    REAL*8::Hvscf(maxbasis(mode),maxbasis(mode)),coef(nDOF,maxbasis(mode))
!
!    CALL subformU(mode,coef)
!    Hvscf=0.d0
!    DO i=0, maxbasis(mode)-1
!       Hvscf(i+1,i+1)= (i + 0.5) * hrmfreq(mode) * 0.5 &!kinetic term
!                                !+Uzero(mode)& !Uzero can be added either here or after diag
!            +Utwo(mode)* harmint(i, mode, 2, 0) &
!            +Ufour(mode)* harmint(i, mode, 4, 0)
!       IF (i .LT. maxbasis(mode)-1) THEN
!          Hvscf(i+1,i+2)=Uone(mode)* harmint(i+1, mode, 1, 1)+Uthree(mode)* harmint(i+1, mode, 3, 1)
!          IF (i .LT. maxbasis(mode)-2) THEN
!             Hvscf(i+1,i+3)= -dSQRT(DBLE((i+2) * (i + 1))) * hrmfreq(mode) * 0.25d0 & !kinetic term
!                  +Utwo(mode)* harmint(i+2, mode, 2, 2)+Ufour(mode)* harmint(i+2, mode, 4, 2)
!             IF (i .LT. maxbasis(mode)-3) THEN
!                Hvscf(i+1,i+4)=Uthree(mode)* harmint(i+3, mode, 3, 3)
!                IF (i .LT. maxbasis(mode)-4) THEN
!                   Hvscf(i+1,i+5)=Ufour(mode)* harmint(i+4, mode, 4, 4)
!                ENDIF
!             ENDIF
!          ENDIF
!       ENDIF
!    ENDDO !i
!  END SUBROUTINE SubFormH
!
!  SUBROUTINE subvscf(state,vscfenergy,vscfiter)
!    USE global
!    USE integral
!    USE modQFF
!    USE lapack
!    USE yazar
!    USE constants
!    USE modtimer
!    IMPLICIT NONE
!    REAL*8:: ham(MAXVAL(maxbasis),MAXVAL(maxbasis)),coef(nDOF,MAXVAL(maxbasis)),newcoef(nDOF,MAXVAL(maxbasis))
!    INTEGER:: mode,iter,state(nDOF)
!    REAL*8, OPTIONAL :: vscfenergy
!    INTEGER,OPTIONAL::vscfiter
!    REAL*8:: energy,oldenergy,energies(MAXVAL(maxbasis)),modalE(nDOF),oldmodalE(nDOF)
!    REAL(8)::wallstart,cpustart,wallend,cpuend
!
!    iter=0
!    oldenergy=0.d0
!    coef=0.d0
!    newcoef=0.d0
!    freq=hrmfreq
!
!    !Initializing using HO solutions according to target state.
!    DO mode=1,nDOF
!       oldmodalE(mode)=hrmfreq(mode)* (DBLE(state(mode)) + 0.5d0)
!       oldenergy=oldenergy+oldmodalE(mode)
!       coef(mode,state(mode)+1)=1.d0
!       newcoef(mode,state(mode)+1)=1.d0
!    ENDDO
!    !start of iterations
!    DO
!       !   PRINT*,'**************************  iterations **************************',iter
!       iter=iter+1
!
!       IF(benchmark)CALL wall_and_cpu_time(wallstart,cpustart)
!
!       DO mode=1,nDOF
!          IF (storepot)THEN
!             CALL subformH(mode,coef,ham)
!             CALL diag(maxbasis(mode),ham,energies)
!          ELSE
!             CALL meanfieldhamiltonian(mode,coef,ham)
!             CALL diag(maxbasis(mode),ham,energies)
!          ENDIF
!          modalE(mode)=energies(state(mode)+1)
!          coef(mode,:)=ham(:,state(mode)+1)
!       ENDDO !mode
!
!       IF(benchmark) THEN
!          CALL wall_and_cpu_time(wallend,cpuend)
!          PRINT*,"wall=",wallend-wallstart,"cpu=",cpuend-cpustart
!       ENDIF
!
!       IF ( (SUM(ABS(modalE-oldmodalE))*convert) .LT. scfthresh)THEN
!
!          DO mode=1,nDOF
!             IF (storepot)THEN
!                CALL subformH(mode,coef,ham)
!                CALL diag(maxbasis(mode),ham,energies)
!                IF(getUzero)CALL subformUzero(mode,coef)
!                IF (mode==nDOF) THEN
!                   !CALL subformUzero(ndof,coef) !alternative 2
!                   !energy=energies(state(mode)+1)+ Uzero(mode)!alternative 2
!                   energy=energies(state(mode)+1)+funUzero(coef) !alternative 1
!                ENDIF
!             ELSE
!                CALL meanfieldhamiltonian(mode,coef,ham)
!                CALL diag(maxbasis(mode),ham,energies)
!                energy=energies(state(mode)+1)
!             ENDIF
!             IF(VMP2 .OR. virtual)THEN
!                vscfenergies(mode,:)=energies
!                vscfcoefs(mode,:,:)=ham
!             ENDIF
!             modalE(mode)=energies(state(mode)+1)
!             coef(mode,:)=ham(:,state(mode)+1)
!          ENDDO !mode
!
!          IF ((ABS(energy-oldenergy)*convert .LT. scfthresh) )EXIT
!          iter=iter+100
!       ENDIF
!
!       IF(iter .GT. maxiter)THEN
!          PRINT*,"maxiter reached"
!          EXIT
!       ENDIF
!
!       oldenergy=energy
!       oldmodalE=modalE
!
!    ENDDO ! Do While
!
!    IF (PRESENT(vscfenergy))vscfenergy=energy
!    IF (PRESENT(vscfiter))vscfiter=iter
!    IF(printall)PRINT*,"VSCF iterations= ",iter
!    IF(getwfn)CALL checkstate(state,coef)
!    IF(getwfn.AND.printall)CALL writemat(coef,ndof,MAXVAL(maxbasis))
!
!    IF(virtual)THEN
!       PRINT*,"vscfenergies"
!       CALL writemat(vscfenergies*convert,ndof,MAXVAL(maxbasis))
!       CALL printvirtual()
!    ENDIF
!    IF(getUzero)THEN
!       PRINT*,"Uzero"
!       CALL writevec(Uzero*convert,ndof)
!    ENDIF
!    ! PRINT*,fungetE(coef)*convert
!    ! print*,'vscf done for', state
!
!  END SUBROUTINE subvscf
!
!  SUBROUTINE subvscfvirtual()
!    USE global
!    USE lapack
!    USE yazar
!    USE constants
!    USE modtimer
!    IMPLICIT NONE
!    REAL*8:: ham(MAXVAL(maxbasis),MAXVAL(maxbasis)),energies(MAXVAL(maxbasis)),coef(nDOF,MAXVAL(maxbasis))
!    INTEGER:: mode,iter,state(nDOF)
!    REAL*8:: energy,oldenergy,newcoef(nDOF,MAXVAL(maxbasis)),modalE(nDOF),oldmodalE(nDOF),Egs,Exc
!    REAL(8)::wallstart,cpustart,wallend,cpuend
!    iter=0
!    oldenergy=1000
!    energy=0.d0
!    coef=0.d0
!    newcoef=0.d0
!    freq=hrmfreq
!    !Forming initial wavefunction using HO basis according to target state.
!    state=0
!    DO mode=1,nDOF
!       oldmodalE(mode)=freq(mode)* (DBLE(state(mode)) + 0.5d0)
!       coef(mode,state(mode)+1)=1.d0
!       newcoef(mode,state(mode)+1)=1.d0
!    ENDDO
!    DO
!       oldenergy=energy
!       energy=0.d0
!       !   PRINT*,'**************************  iteration **************************',iter
!       iter=iter+1
!       IF(benchmark)CALL wall_and_cpu_time(wallstart,cpustart)
!       DO mode=1,nDOF
!          IF (storepot)THEN
!             CALL subformH(mode,coef,ham)
!             CALL diag(maxbasis(mode),ham,energies)
!             IF (mode==nDOF) THEN
!                !CALL subformUzero(ndof,coef) !alternative 2
!                !energy=energies(state(mode)+1)+ Uzero(mode)!alternative 2
!                energy=energies(state(mode)+1)+funUzero(coef) !alternative 1
!             ENDIF
!          ELSE
!             CALL meanfieldhamiltonian(mode,coef,ham)
!             CALL diag(maxbasis(mode),ham,energies)
!             energy=energies(state(mode)+1)
!          ENDIF
!          ! newcoef(mode,:)=ham(:,state(mode)+1)
!          modalE(mode)=energies(state(mode)+1)
!          coef(mode,:)=ham(:,state(mode)+1)
!          ! PRINT*,'eigvec'!ham(1,:)
!       ENDDO !mode
!       IF(benchmark) THEN
!          CALL wall_and_cpu_time(wallend,cpuend)
!          PRINT*,"wall=",wallend-wallstart,"cpu=",cpuend-cpustart
!       ENDIF
!       ! coef=newcoef
!       IF(iter .GT. maxiter)THEN
!          PRINT*,"maxiter reached"
!          EXIT
!       ENDIF
!       !IF ((ABS(energy-oldenergy) .LT. scfthresh))THEN
!       IF ( (SUM(ABS(modalE-oldmodalE))*convert) .LT. scfthresh)THEN
!          DO mode=1,nDOF
!             IF (storepot)THEN
!                CALL subformH(mode,coef,ham)
!                CALL diag(maxbasis(mode),ham,energies)
!                IF (mode==nDOF) THEN
!                   !CALL subformUzero(ndof,coef) !alternative 2
!                   !energy=energies(state(mode)+1)+ Uzero(mode)!alternative 2
!                   energy=energies(state(mode)+1)+funUzero(coef) !alternative 1
!                ENDIF
!             ELSE
!                CALL meanfieldhamiltonian(mode,coef,ham)
!                CALL diag(maxbasis(mode),ham,energies)
!                energy=energies(state(mode)+1)
!             ENDIF
!             vscfenergies(mode,:)=energies
!             IF(vmp2)   vscfcoefs(mode,:,:)=ham
!
!             ! newcoef(mode,:)=ham(:,state(mode)+1)
!             coef(mode,:)=ham(:,state(mode)+1)
!             ! PRINT*,'eigvec'!ham(1,:)
!          ENDDO !mode
!          IF ((ABS(energy-oldenergy)*convert .LT. scfthresh) )EXIT
!
!       ENDIF
!       oldmodalE=modalE
!    ENDDO ! Do While
!    energy=energies(1)+funUzero(coef)
!
!    PRINT*,'VVSCF fundamental energies'
!    Egs=energy*convert
!    PRINT*,Egs
!    DO mode=1,ndof
!       Exc=(vscfenergies(mode,2)-vscfenergies(mode,1))*convert
!       PRINT*,Exc!,Exc+Egs!,(Energy-theenergy)*convert
!       !   print*, modalE(mode)*convert,(energy+modalE(mode))*convert
!    ENDDO
!
!    IF(printall)PRINT*,"VSCF iterations= ",iter
!    IF(getwfn)CALL checkstate(state,coef)
!    IF(getwfn.AND.printall)CALL writemat(coef,ndof,MAXVAL(maxbasis))
!    ! if(virtual)call printvirtual()
!    !    PRINT*,fungetE(coef)*convert
!    ! print*,'vscf done for', state
!    RETURN
!  END SUBROUTINE subvscfvirtual
!
!
!  SUBROUTINE checkstate(state,coef)
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    INTEGER,INTENT(in)::state(nDOF)
!    REAL*8,INTENT(in):: coef(nDOF,MAXVAL(maxbasis))
!    REAL*8:: sumcoef,prob(MAXVAL(maxbasis))
!    INTEGER::i,m,stateout(nDOF)
!    !print*,maxloc(abs(coef),2),"kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk"
!    stateout=MAXLOC(ABS(coef),2)-1
!    !test
!    ! print*,"Converged state is", stateout , "for initial HO state", state
!
!    IF( COUNT((stateout-state)==0) .NE. ndof ) THEN
!       rightstate=.FALSE.
!       PRINT*,"Converged state is", stateout , "for initial HO state", state
!    ENDIF
!    DO m=1,ndof
!       !weird if( (maxloc(coef(m,:))) .ne. (state(m)+1) ) rightstate=.false. gives error
!       !says if clause is not scalar.
!       prob=SQRT(coef(m,:)*coef(m,:))
!       IF (MAXVAL(prob) < 0.90) THEN
!          sumcoef=0.d0
!          PRINT*,"resonance for mode ",m, "for state", stateout
!          DO i=0,maxbasis(m)-1
!             PRINT*,prob(i+1)
!             sumcoef=sumcoef+prob(i+1)
!             IF(sumcoef > 0.95)EXIT
!          ENDDO
!       ENDIF
!    ENDDO !m
!  END SUBROUTINE checkstate
!
!  SUBROUTINE printvirtual
!    !VSCF for given state and calculate excited states
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    REAL*8 :: Egs,Exc
!    INTEGER::m,state(nDOF)
!    state=0
!    PRINT*,"Virtual VSCF fund energies"
!    Egs=fungetE(state)
!    PRINT*,Egs*convert
!    ! PRINT*,groundE*convert
!    DO m=1,ndof
!       state=0
!       state(m)=1
!       Exc=fungetE(state)
!       PRINT*,(vscfenergies(m,2)-vscfenergies(m,1))*convert,(Exc-Egs)*convert, Exc*convert
!       ! CALL getenergy(state,energy)
!       ! PRINT*,(energy-groundE)*convert,energy*convert
!       ! PRINT*,"fungetE ",fungetE(state)*convert
!    ENDDO
!    state(1)=2
!    state(2)=1
!    state(3)=0
!    Exc=fungetE(state)
!    PRINT*,(Exc-Egs)*convert,Exc*convert,"for state ",state
!    PRINT*,"************************"
!  END SUBROUTINE printvirtual
!
!  SUBROUTINE printssVSCF
!    !VSCF for given state and calculate excited states
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    REAL*8 :: Egs,Exc
!    INTEGER::state(nDOF),iter
!
!    state(1)=2
!    state(2)=1
!    state(3)=0
!    CALL subvscf(state,Exc,iter)
!    PRINT*,(Exc)*convert,"for SS vscf state ",state, "iter", iter
!    Exc=fungetE(state)
!    state=0
!    Egs=fungetE(state)
!    PRINT*,(Exc-Egs)*convert,Exc*convert,"for state ",state
!    PRINT*,"************************"
!  END SUBROUTINE printssVSCF
!
!  REAL*8 FUNCTION meanfieldkinetic(mode,coef)
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    INTEGER::i
!    REAL*8 :: tmp
!    INTEGER,INTENT(in)::mode
!    REAL*8,INTENT(in)::coef(nDOF,maxbasis(mode))
!
!
!    tmp = 0.d0;
!    !First method
!    !  DO i=0, maxbasis(mode)-1
!    !     DO j=0, maxbasis(mode)-1
!    !        diff = j - i;
!    !        IF (diff .EQ. 0) THEN
!    !           tmp =tmp+ coef(mode,i+1) * coef(mode,i+1) * (dble(i) + 0.5d0) * hrmfreq(mode) * 0.5d0;
!    !        ELSE IF (ABS(diff) .EQ. 2) THEN
!    !          ! tmp =tmp- coef(mode,i+1) * coef(mode,j+1) * dSQRT(DBLE( MAX(i, j) * (MAX(i, j) - 1) )) * hrmfreq(mode) * 0.25d0;
!    !        ENDIF
!    !     ENDDO!j
!    !  ENDDO!i
!
!    !second method
!    DO i=0, maxbasis(mode)-1
!       tmp =tmp+ coef(mode,i+1) * coef(mode,i+1) * (DBLE(i) + 0.5d0) * hrmfreq(mode) * 0.5d0
!       IF(i.LT.(maxbasis(mode)-2))&
!            tmp =tmp- coef(mode,i+1) * coef(mode,i+3) * 0.5d0 * hrmfreq(mode)* dSQRT(DBLE((i+2)*(i+1)));
!    ENDDO!i
!    !end of 2nd
!    !3rd way
!    ! tmp=modalexp(mode,-2,coef)
!    meanfieldkinetic=tmp
!    ! print*, tmp*convert,"meanfield kinetic for" ,mode
!    RETURN
!  END FUNCTION meanfieldkinetic
!  REAL*8 FUNCTION meanfieldpotential(mode,coef)
!    USE global
!    USE integral
!    USE modQFF
!    USE constants
!    IMPLICIT NONE
!    INTEGER::m,m1,m2,m3
!    REAL*8 :: tmp!,meanFieldKinetic
!    INTEGER,INTENT(in)::mode
!    REAL*8,INTENT(in)::coef(nDOF,maxbasis(mode))
!    tmp = 0.0;
!    ! WRITE(*,'(d50.25)')c_gamma
!    DO m=1,Ndof
!       IF (m .EQ. mode)CYCLE
!       tmp =tmp+ modalexp(m, 1,coef) * gi(m);
!       tmp =tmp+ modalexp(m, 2,coef) * Hii(m);
!       tmp =tmp+ modalexp(m, 3,coef) * Ciii(m);
!       tmp =tmp+ modalexp(m, 4,coef) * Qiiii(m);
!       tmp =tmp+ meanFieldKinetic(m,coef);
!    ENDDO
!    IF (nMR .GT. 1)THEN
!       ! print*, "**************************nMR >1 "
!       DO m1=1,nDOF
!          IF (m1 .EQ. mode)CYCLE
!          DO m2=m1+1,nDOF
!             IF (m2 .EQ. mode .OR. m1 .EQ. m2)CYCLE
!             tmp =tmp+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * Hij(m1,m2);
!             tmp =tmp+ modalexp(m1, 2,coef) * modalexp(m2, 2,coef) * Qiijj(m1,m2);
!             tmp =tmp+ modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * Ciij(m1,m2);
!             tmp =tmp+ modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * Ciij(m2,m1);
!             tmp =tmp+ modalexp(m1, 3,coef) * modalexp(m2, 1,coef) * Qiiij(m1,m2);
!             tmp =tmp+ modalexp(m1, 1,coef) * modalexp(m2, 3,coef) * Qiiij(m2,m1);
!          ENDDO!m2
!       ENDDO!m1
!       IF (nMR .GT. 2)THEN
!          DO m1=1,nDOF
!             IF (m1 .EQ. mode) CYCLE
!             DO m2=m1+1,nDOF
!                IF (m2 .EQ. mode) CYCLE
!                DO m3=m2+1,nDOF
!                   IF (m3 .EQ. mode) CYCLE
!                   tmp =tmp+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef) * Cijk(m1,m2,m3);
!                   tmp =tmp+ modalexp(m1, 2,coef) * modalexp(m2, 1,coef) * modalexp(m3, 1,coef) * Qiijk(m1,m2,m3);
!                   tmp =tmp+ modalexp(m1, 1,coef) * modalexp(m2, 2,coef) * modalexp(m3, 1,coef) * Qiijk(m2,m1,m3);
!                   tmp =tmp+ modalexp(m1, 1,coef) * modalexp(m2, 1,coef) * modalexp(m3, 2,coef) * Qiijk(m3,m1,m2);
!                ENDDO!m3
!             ENDDO!m2
!          ENDDO!m1
!       ENDIF!nMR>2
!    ENDIF!nMR>1
!    meanfieldpotential=tmp
!    !  print*, tmp*convert,"meanfield pot for" ,mode
!    RETURN
!  END FUNCTION meanfieldpotential
!
!  SUBROUTINE meanFieldHamiltonian(mode,coef,Hvscf)
!    USE global
!    USE integral
!    USE modQFF
!    IMPLICIT NONE
!    INTEGER::diff,i,j,m,m1,m2
!    REAL*8::tmp,meanpotkin!,meanfieldpotential
!    INTEGER,INTENT(in)::mode
!    REAL*8,INTENT(in)::coef(nDOF,maxbasis(mode))
!    REAL*8::Hvscf(maxbasis(mode),maxbasis(mode))
!    !PRINT*, "forming VSCF Hamiltonian for mode" , mode
!    meanpotkin=meanFieldPotential(mode,coef)
!    ! PRINT*,'meanpotkin',meanpotkin
!    Hvscf=0.d0
!
!    DO i=0, maxbasis(mode)-1
!       tmp=0.0
!       tmp= tmp + (i + 0.5) * hrmfreq(mode) * 0.5;!kinetic term
!       tmp= tmp + Hii(mode) * harmint(i, mode, 2, 0);
!       tmp= tmp + Qiiii(mode) * harmint(i, mode, 4, 0);
!       tmp= tmp + meanpotkin
!       IF (nMR .GT. 1) THEN
!          DO m=1,ndof
!             IF (m .EQ. mode)CYCLE
!             tmp= tmp + Hij(m,mode) * harmint(i, mode, 1, 0) * modalexp(m, 1,coef);
!             tmp= tmp + Qiijj(m,mode) * harmint(i, mode, 2, 0) * modalexp(m, 2,coef);
!             tmp= tmp + Ciij(mode,m) * harmint(i, mode, 2, 0) * modalexp(m, 1,coef);
!             tmp= tmp + Ciij(m,mode) * harmint(i, mode, 1, 0) * modalexp(m, 2,coef);
!             tmp= tmp + Qiiij(mode,m) * harmint(i, mode, 3, 0) * modalexp(m, 1,coef);
!             tmp= tmp + Qiiij(m,mode) * harmint(i, mode, 1, 0) * modalexp(m, 3,coef);
!          ENDDO
!          IF (nMR .GT. 2) THEN
!             DO m1=1,ndof-2
!                IF (m1 .EQ. mode)CYCLE
!                DO m2=m1+1,ndof-1
!                   IF (m2 .EQ. mode) CYCLE
!                   tmp= tmp + Cijk(mode,m1,m2) * harmint(i, mode, 1, 0) * modalexp(m1, 1,coef)  * modalexp(m2, 1,coef);
!                   tmp= tmp + Qiijk(mode,m1,m2) * harmint(i, mode, 2, 0) * modalexp(m1, 1,coef) * modalexp(m2, 1,coef);
!                   tmp= tmp + Qiijk(m1,mode,m2) * harmint(i, mode, 1, 0) * modalexp(m1, 2,coef)* modalexp(m2, 1,coef);
!                   tmp= tmp + Qiijk(m2,mode,m1) * harmint(i, mode, 1, 0) * modalexp(m1, 1,coef) * modalexp(m2, 2,coef);
!                ENDDO!!m2
!             ENDDO!m1
!          ENDIF!nMR>2
!       ENDIF!nMR>1
!       Hvscf(i+1,i+1)=tmp
!       DO j=i+1,maxbasis(mode)-1!
!          tmp=0.0
!          diff = j - i;
!          IF (diff > 4) EXIT !break;
!          IF (diff .EQ. 2) tmp= tmp-dSQRT(DBLE(j * (j - 1))) * hrmfreq(mode) * 0.25 !!kinetic term
!          tmp= tmp + harmint(j, mode, 1, diff) * gi(mode);
!          tmp= tmp + harmint(j, mode, 2, diff) * Hii(mode);
!          tmp= tmp + harmint(j, mode, 3, diff) * Ciii(mode);
!          tmp= tmp + harmint(j, mode, 4, diff) * Qiiii(mode);
!          ! print*,'test2',tmp
!          IF (nMR .GT. 1)THEN
!             DO m=1,ndof
!                IF (m .EQ. mode)CYCLE
!                tmp= tmp + Hij(m,mode) * harmint(j, mode, 1, diff) * modalexp(m, 1,coef);
!                tmp= tmp + Qiijj(m,mode) * harmint(j, mode, 2, diff) * modalexp(m, 2,coef);
!                tmp= tmp + Ciij(mode,m) * harmint(j, mode, 2, diff) * modalexp(m, 1,coef);
!                tmp= tmp + Ciij(m,mode) * harmint(j, mode, 1, diff) * modalexp(m, 2,coef);
!                tmp= tmp + Qiiij(mode,m) * harmint(j, mode, 3, diff) * modalexp(m, 1,coef);
!                tmp= tmp + Qiiij(m,mode) * harmint(j, mode, 1, diff) * modalexp(m, 3,coef);
!             ENDDO
!             IF (nMR .GT. 2) THEN
!                DO m1=1,ndof
!                   IF (m1 .EQ. mode)CYCLE
!                   DO m2=m1 + 1,ndof
!                      IF (m2 .EQ. mode)CYCLE
!                      tmp= tmp + Cijk(mode,m1,m2)*harmint(j,mode,1,diff) * modalexp(m1,1,coef)* modalexp(m2,1,coef);
!                      tmp= tmp + Qiijk(mode,m1,m2)*harmint(j,mode,2,diff) * modalexp(m1,1,coef)* modalexp(m2,1,coef);
!                      tmp= tmp + Qiijk(m1,mode,m2)*harmint(j,mode,1,diff) * modalexp(m1,2,coef)* modalexp(m2,1,coef);
!                      tmp= tmp + Qiijk(m2,mode,m1)*harmint(j,mode,1,diff) * modalexp(m1,1,coef)* modalexp(m2,2,coef);
!                   ENDDO!!m2
!                ENDDO!!m1
!             ENDIF!nMR>2
!          ENDIF!nMR>1
!          ! print*,'test3',tmp
!          Hvscf(i+1,j+1)=tmp
!       ENDDO!!end for j
!    ENDDO!!end for i
!  END SUBROUTINE meanFieldHamiltonian
!END MODULE modVSCFold
