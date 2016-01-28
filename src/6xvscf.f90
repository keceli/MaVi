
MODULE modXVSCF
  IMPLICIT NONE
  PRIVATE::funGetE,printfreq,printvirtual

CONTAINS

    SUBROUTINE init_XVSCF
    USE global
    IMPLICIT NONE
    if(debug)print*,'init_XVSCF'
    if (allocated(uzero))deallocate(uzero)
    if (allocated(utwo))deallocate(utwo)
    if(storepot)ALLOCATE(uZero(ndof),uTwo(ndof))
    if(debug)print*,'init_XVSCF :)'
  END SUBROUTINE init_XVSCF

  REAL*8 FUNCTION oldfunUzero(state,mode)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,mode,state(nDOF)
    oldfunUzero=0.d0
    DO m1=1,ndof
       IF (m1 .EQ. mode)CYCLE
       !print*,'m1= ',m1
       !funUzero= funUzero + 0.5*Fii(m1) * harmint(state(m1), m1, 2, 0)
       oldfunUzero= oldfunUzero + Hii(m1) * harmint(state(m1), m1, 2, 0)
       !print*,'funUzero= ',funUzero
       DO m2=1,ndof
          IF (m2 .EQ. mode .OR. m1 .EQ. m2)CYCLE
          !funUzero= funUzero + 0.125d0 * Fiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
          oldfunUzero= oldfunUzero + 0.5d0 * Qiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
       ENDDO
    ENDDO
  END FUNCTION oldfunUzero

  REAL*8 FUNCTION funUzero(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(nDOF)
    !Uzero calculation for nDOFth mode
    !m=ndof
    funUzero=0.d0
    DO m1=1,ndof-1
       !IF(.NOT.ateqb)funUzero= funUzero + Gi(m1)*DSQRT(n * 0.5d0 / gamma)
       !funUzero= funUzero + 0.5*Fii(m1) * harmint(state(m1), m1, 2, 0)
       ! funUzero= funUzero +0.5d0 * freq(m1)* (DBLE(state(m1)) + 0.5d0)&
       !     + Hii(m1) * (DBLE(state(m1)) + 0.5d0) / freq(m1)
       funUzero= funUzero + ( DBLE(state(m1)) + 0.5d0 ) * ( 0.5d0*freq(m1) + Hii(m1)/freq(m1) )
       !print*,'funUzero= ',funUzero
       if(nMR>1)then
       DO m2=m1+1,ndof-1
          !  funUzero= funUzero + 0.125d0 * Fiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
          funUzero= funUzero + Qiijj(m1,m2) * (DBLE(state(m1)) + 0.5d0) / freq(m1) * (DBLE(state(m2)) + 0.5d0) / freq(m2)
       ENDDO
       endif
    ENDDO
  END FUNCTION funUzero

  SUBROUTINE subUzero(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m,m1,m2,state(nDOF)
    !Uzero calculation for nDOFth mode
    !m=ndof
    DO m=1,ndof
       Uzero(m)=0.d0
       DO m1=1,ndof
          IF (m1 .EQ. m)CYCLE
          Uzero(m)= Uzero(m) + ( DBLE(state(m1)) + 0.5d0 ) * ( 0.5d0*freq(m1) + Hii(m1)/freq(m1) )
          DO m2=m1+1,ndof
             IF (m2 .EQ. m)CYCLE
             Uzero(m)= Uzero(m) + Qiijj(m1,m2) * (DBLE(state(m1)) + 0.5d0) / freq(m1) * (DBLE(state(m2)) + 0.5d0) / freq(m2)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE subUzero

  REAL*8 FUNCTION funUtwo(state,mode)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m,mode,state(nDOF)
    funUtwo=0.d0
    DO m=1,ndof
       IF (m .EQ. mode)CYCLE
       !funUtwo= funUtwo + 0.5d0 *Fiijj(mode,m) * harmint(state(m), m, 2, 0)
       funUtwo= funUtwo + Qiijj(mode,m) * harmint(state(m), m, 2, 0)
    ENDDO
    !funUtwo= funUtwo + Fii(mode)
    funUtwo= funUtwo + Hii(mode)
    ! if(funUtwo<0.d0)print*,"negative funUtwo for state,mode", state,mode
  END FUNCTION funUtwo

  REAL*8 FUNCTION funInt(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(nDOF)
    funInt=0.d0
    DO m1=1,ndof
       DO m2=m1+1,ndof
          !funUtwo= funUtwo + 0.5d0 *Fiijj(mode,m) * harmint(state(m), m, 2, 0)
          funInt= funInt + Qiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
       ENDDO
    ENDDO
    RETURN
  END FUNCTION funInt

  REAL*8 FUNCTION funInt2(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(nDOF)
    funInt2=0.d0
    DO m1=1,ndof
       DO m2=m1+1,ndof
          !funUtwo= funUtwo + 0.5d0 *Fiijj(mode,m) * harmint(state(m), m, 2, 0)
          funInt2= funInt2 + Qiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
       ENDDO
    ENDDO
    RETURN
  END FUNCTION funInt2

  SUBROUTINE subUtwo(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m,mode,state(nDOF)
    DO mode=1,ndof
       Utwo(mode)=0.d0
       DO m=1,ndof
          IF (m .EQ. mode)CYCLE
          !funUtwo= funUtwo + 0.5d0 *Fiijj(mode,m) * harmint(state(m), m, 2, 0)
          Utwo(mode)= Utwo(mode) + Qiijj(mode,m) * harmint(state(m), m, 2, 0)
       ENDDO
       Utwo(mode)= Utwo(mode) + Hii(mode)
    ENDDO !mode
    !funUtwo= funUtwo + Fii(mode)

    ! if(funUtwo<0.d0)print*,"negative funUtwo for state,mode", state,mode
  END SUBROUTINE subUtwo

  REAL*8 FUNCTION Ukinetic(state,mode)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER ::m,state(nDOF),mode
    ! <n|del^2/del(q)^2|n-diff>
    Ukinetic=0.d0
    DO m=1,ndof
       IF (m .EQ. mode)CYCLE
       Ukinetic= Ukinetic+0.5 * freq(m)* (state(m) + 0.5)
    ENDDO
    RETURN
  END FUNCTION Ukinetic

  REAL*8 FUNCTION Ekinetic(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER ::m,state(nDOF)
    ! <n|del^2/del(q)^2|n-diff>
    Ekinetic=0.d0
    DO m=1,ndof
       Ekinetic= Ekinetic+0.5d0 * freq(m)* (DBLE(state(m)) + 0.5d0)
    ENDDO
    RETURN
  END FUNCTION Ekinetic

  REAL*8 FUNCTION Epotential(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(nDOF)
    Epotential=0.d0
    DO m1=1,ndof
       !Epotential= Epotential + 0.5*Fii(m1) * harmint(state(m1), m1, 2, 0)
       Epotential= Epotential + Hii(m1) * harmint(state(m1), m1, 2, 0)
       DO m2=1,ndof
          IF ( m1 .EQ. m2)CYCLE
          !Epotential= Epotential + 0.125d0 * Fiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
          Epotential= Epotential + 0.5d0 * Qiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
       ENDDO
    ENDDO
  END FUNCTION Epotential

  REAL*8 FUNCTION fungetE(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(nDOF)
    REAL*8::tmp
    tmp=0.d0
    DO m1=1,ndof
       tmp= tmp+0.5d0 * freq(m1)* (DBLE(state(m1)) + 0.5d0)& !kinetic term
                                !    + Hii(m1) * harmint(state(m1), m1, 2, 0)
            + Hii(m1) *((DBLE(state(m1))+0.5d0)/freq(m1))
       DO m2=m1+1,ndof
          !tmp=tmp+ Qiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
          tmp=tmp+ Qiijj(m1,m2) * ((DBLE(state(m1))+0.5d0)/freq(m1))* ((DBLE(state(m2))+0.5d0)/freq(m2))
       ENDDO
    ENDDO
    fungetE=tmp
    RETURN
  END FUNCTION fungetE

  REAL*8 FUNCTION fungetE2(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(nDOF)
    REAL*8::tmp
    tmp=0.d0
    DO m1=1,ndof
       tmp= tmp+ (freq(m1)*freq(m1) + 2.d0*Hii(m1)) *(DBLE(state(m1))+0.5d0)/(2.d0*freq(m1))
   ENDDO
   DO m1=1,ndof
       DO m2=1,ndof
          tmp=tmp+ Qiijj(m1,m2) * ((DBLE(state(m1))+0.5d0)/freq(m1))* ((DBLE(state(m2))+0.5d0)/freq(m2))/2.d0
       ENDDO
    ENDDO
    fungetE2=tmp
    RETURN
  END FUNCTION fungetE2

  SUBROUTINE subgetE(state,energy)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(nDOF)
    REAL*8::tmp,energy
    tmp=0.d0
    DO m1=1,ndof
       tmp= tmp+0.5d0 * freq(m1)* (DBLE(state(m1)) + 0.5d0)& !kinetic term
                                !    + Hii(m1) * harmint(state(m1), m1, 2, 0)
            + Hii(m1) *((DBLE(state(m1))+0.5d0)/freq(m1))
       DO m2=m1+1,ndof
          !tmp=tmp+ Qiijj(m1,m2) * harmint(state(m1), m1, 2, 0)* harmint(state(m2), m2, 2, 0)
          tmp=tmp+ Qiijj(m1,m2) * ((DBLE(state(m1))+0.5d0)/freq(m1))* ((DBLE(state(m2))+0.5d0)/freq(m2))
       ENDDO
    ENDDO
    energy=tmp
    RETURN
  END SUBROUTINE subgetE

  SUBROUTINE subqVSCF(state,qvscfenergy,qVSCFiter)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    USE modtimer
    IMPLICIT NONE
    INTEGER,OPTIONAL::qVSCFiter
    REAL*8,OPTIONAL:: qvscfenergy
    INTEGER:: mode,iter,failed,state(nDOF)
    REAL*8:: energy,oldenergy,oldfreq(ndof),tmp
    REAL(8)::wallstart,cpustart,wallend,cpuend

    iter=0
    failed=0
    oldenergy=1000
    energy=0.0
    freq=hrmfreq
    !    DO WHILE ((DABS(energy-oldenergy) .GT. scfthresh) .AND. iter .LT. maxiter )
    !       oldenergy=energy
    !       iter=iter+1
    !       DO mode=1,nDOF
    !          freq(mode)=dsqrt(2.d0*funUtwo(state,mode))
    !       ENDDO !mode
    !      ! energy=Ekinetic(state)+Epotential(state) !alt 1
    !       energy=fungetE(state) !alt 2
    !    ENDDO !while
    DO ! iterations
       oldfreq=freq
       iter=iter+1
       IF(benchmark)CALL wall_and_cpu_time(wallstart,cpustart)
       DO mode=1,nDOF
          !   freq(mode)=dsqrt(funUtwo(state,mode))
          tmp=DOT_PRODUCT( Qiijj(mode,:) , (DBLE(state)+0.5d0) / freq ) + Hii(mode)
          IF(tmp<0.d0)THEN
             freq(mode)=dsqrt(-2.d0*tmp)
             IF(failed<1)  PRINT*,"!!!!********!!!!!!complex freq ", mode
             failed=1
          ELSE
             freq(mode)=dsqrt(2.d0*tmp)
          ENDIF
       ENDDO !mode
       IF(benchmark) THEN
          CALL wall_and_cpu_time(wallend,cpuend)
          PRINT*,"wall=",wallend-wallstart,"cpu=",cpuend-cpustart
       ENDIF
       !   print*,"conv",SUM(abs(freq-oldfreq))
       IF (((SUM(ABS(freq-oldfreq))*convert) < scfthresh) .OR. (iter > maxiter)) EXIT !converged or exceed maxiter
    ENDDO !iteration
    IF(iter>maxiter .AND. maxiter>1)THEN
       PRINT*,'not converged after ',iter,'iterations. max for abs(freq-oldfreq)= ',MAXVAL(ABS(freq-oldfreq))
       IF(SUM(state)==0) THEN
          PRINT*," failed for ground state"
       ELSE
          PRINT*," failed for fundamental ",MAXLOC(state)
       ENDIF
       energy=0.d0
    ELSE
       ! energy=fungetE(state)
       energy=funUzero(state)+(0.5d0+state(ndof))*freq(ndof)
    ENDIF
    IF(debug)PRINT*,"qVSCF iterations= ",iter
    IF(PRESENT(qVSCFiter))qVSCFiter=iter
    IF(PRESENT(qVSCFenergy))qVSCFenergy=energy
    IF(virtual)CALL printvirtual(state)
    IF(getfreq)CALL printfreq(state)
    RETURN
    !  energy=funUzero(state)+0.5d0*sum(freq)+dot_product(state,freq)
  END SUBROUTINE subqVSCF

  SUBROUTINE subXVSCF()
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    USE modtimer
    IMPLICIT NONE
    INTEGER:: mode,iter,failed,state(ndof)
    REAL*8:: energy,oldfreq(ndof),tmp
    REAL(8)::wallstart,cpustart,wallend,cpuend
    if(debug)print*,'subXVSCF'
    state=0
    iter=0
    failed=0
    energy=0.0
    freq=hrmfreq
    if(debug)print*,'subXVSCF iterations'
    DO ! iterations
       oldfreq=freq
       iter=iter+1
       IF(benchmark)CALL wall_and_cpu_time(wallstart,cpustart)
       DO mode=1,nDOF
          !   freq(mode)=dsqrt(funUtwo(state,mode))
          if(nMR>1)then
          tmp=DOT_PRODUCT( Qiijj(mode,:) , (0.5d0 / freq) ) + Hii(mode)
          else
          tmp= Hii(mode)
          endif
          IF(tmp<0.d0)THEN
             freq(mode)=dsqrt(-2.d0*tmp)
             IF(failed<1)  PRINT*,"!!!!************************!!!!!!complex freq ", mode ,iter
             failed=1
          ELSE
             freq(mode)=dsqrt(2.d0*tmp)
          ENDIF
       ENDDO !mode
       IF(benchmark) THEN
          CALL wall_and_cpu_time(wallend,cpuend)
          PRINT*,"wall=",wallend-wallstart,"cpu=",cpuend-cpustart
       ENDIF
       !   print*,"conv",SUM(abs(freq-oldfreq))
       IF ((SUM(ABS(freq-oldfreq)) < scfthresh) .OR. (iter > maxiter)) EXIT !converged or exceed maxiter
    ENDDO !iteration
    IF(iter>maxiter .AND. maxiter>1)THEN
       PRINT*,'not converged after ',iter,'iterations. max for abs(freq-oldfreq)= ',MAXVAL(ABS(freq-oldfreq))
       energy=0.d0
    ELSE
       ! energy=fungetE(state)
       energy=funUzero(state)+0.5d0*freq(ndof)
    ENDIF
    PRINT*,'XVSCF'
    PRINT*,energy*convert
    DO mode=1,ndof
       PRINT*, freq(mode)*convert!,(energy+freq(mode))*convert
    ENDDO
    IF(debug)PRINT*,"XVSCF iterations= ",iter
    RETURN
    !  energy=funUzero(state)+0.5d0*sum(freq)+dot_product(state,freq)
  END SUBROUTINE subXVSCF

  SUBROUTINE printfreq(state)
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    INTEGER::state(nDOF),m
    PRINT*,"************************"
    PRINT*,"frequencies after qVSCF calc for state", state
    PRINT*,SUM(freq)/2.d0*convert
    PRINT*,DOT_PRODUCT(dfloat(state)+0.5d0,freq)*convert
    PRINT*,SUM(freq)/2.d0*convert-funInt(state)*convert
    PRINT*,SUM(hrmfreq)/2.d0*convert+funInt2(state)*convert
    DO m=1,ndof
       PRINT*,freq(m)*convert
    ENDDO
    PRINT*,"************************"

  END SUBROUTINE printfreq

  SUBROUTINE printvirtual(thestate)
    !VSCF for given state and calculate excited states
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    IMPLICIT NONE
    REAL*8 :: energy,groundE
    INTEGER,INTENT(IN)::thestate(nDOF)
    INTEGER::m,state(nDOF)

    PRINT*,"************************"
    PRINT*,"virtual fund energies for method ",method," based on state ", thestate
    state=0
    groundE=fungetE(state)
    PRINT*,groundE*convert
    DO m=1,ndof
       state=0
       state(m)=1
       Energy=fungetE(state)
       PRINT*,Energy*convert,(Energy-groundE)*convert
       ! CALL getenergy(state,energy)
       ! PRINT*,(energy-groundE)*convert,energy*convert
       ! PRINT*,"fungetE ",fungetE(state)*convert
    ENDDO
    PRINT*,"************************"
  END SUBROUTINE printvirtual

  SUBROUTINE groundstate()
  USE global!,ONLY:ndof,freq,hrmfreq,maxiter,scfthresh,benchmark,Uzero,Utwo
  USE modQFF
  USE integral
    USE constants
    IMPLICIT NONE
    REAL*8 :: energy
    INTEGER::iter,state(nDOF)
    state=0
    CALL subqVSCF(state,energy,iter)
    PRINT*,'Groud state energy is', energy*convert ,'after',iter,'iterations'
  END SUBROUTINE groundstate

END MODULE modXVSCF

