!*********subroutines**************************************************************************************

SUBROUTINE readinput
  USE constants
  USE modQFF!,ONLY:quartic,harmonic,allzero,fullquartic,ateqb,getQFFplot,qfftype,fcfile
  USE modformPES
  USE modalintegral,ONLY:init_modalint
  USE modVSCF,ONLY:init_VSCF
  USE modXVSCF,ONLY:init_XVSCF
  USE modVCI
  IMPLICIT NONE
  INTEGER,PARAMETER  :: p=1000
  INTEGER::i,j,k1,k2,m,n,ierr
  INTEGER::Nat,nfree,MR,vmax_base,vmax(p)
  REAL*8 ::x(3*p),L(p*p),mass(p),omega(p)
  CHARACTER::Label(p)*2
  NAMELIST /airun/ RunTyp
  NAMELIST /inp/ Com
  NAMELIST /mol/ x,L,mass,Label,Nat,omega,nfree
  NAMELIST /vib/ MR, vmax, vmax_base, vscf, vci,vpt,vav,summary,vmp2,xvscf,Vvscf,vscfci,vpt20,vci1,vscfci1,&
  oneMRP1,oneMRP2,ssvmp
  NAMELIST /gen/storepot,runmonomer,getfunds,maxiter,Nmonomer,virtual,scfthresh,&
  benchmark,unittype,method,getfreq,getwfn,&
  maxiter,ssvmp,ssvscf,debug,storeint,testgreen,denomcut,useXVSCFok,test,overtone
  NAMELIST /pes/ateqb,getQFFplot,qfftype,fcfile,quartic,fullquartic,allzero,harmonic,writeQFF,printaverageOK,nocubic,&
  noCiij,noCiii,kcubic,kquartic,justCiii,justQiiii,justCiij,justQiijj,justQiiij,justCijk,justQiijk,justQijkl

  !Default Values
  RunTyp='VIB'
  denomcut=0.d0
  vmax=0
  vmax_base=8
  MR=4
  FCIsize=1
  method=''
  oneMRP1=.false.
  oneMRP2=.false.
  VCI1=.false.
  VSCFCI1=.false.
  overtone=.false.
  test=.false.
  useXVSCFok=.FALSE.
  writeQFF=.FALSE.
  testgreen=.FALSE.
  ssvmp=.FALSE.
  debug=.FALSE.
  ssvscf=.FALSE.
  Vvscf=.FALSE.
  vscfci=.false.
  xvscf=.FALSE.
  vscf=.FALSE.
  vci=.FALSE.
  vpt=.FALSE.
  vmp2=.FALSE.
  vpt20=.FALSE.
  summary=.FALSE.
  storepot=.TRUE.
  storeint=.TRUE.
  runmonomer=.FALSE.
  getfunds=.TRUE.
  virtual=.FALSE.
  getFreq=.FALSE.
  getwfn=.FALSE.
  benchmark=.FALSE.
  delx=0.5d0
  Nmonomer=1
  maxNstate=400
  maxFF=4
  maxiter=500
  scfthresh=1.0d-10
  inttol=1.0D-10
  convert=c_h2wn
  unittype=1
  omega=0.d0
  !QFF
  kcubic=1.d0
  kquartic=1.d0
  quartic=.FALSE.
  fullquartic=.FALSE.
  harmonic=.FALSE.
  allzero=.FALSE.
  nocubic=.FALSE.
  noCiij=.FALSE.
  justCiii=.FALSE.
  justQiiii=.FALSE.
  justCiij=.FALSE.
  justQiiij=.FALSE.
  justQiijj=.FALSE.
  justCijk=.FALSE.
  justQiijk=.FALSE.
  justQijkl=.FALSE.
  noCiii=.FALSE.
  fcfile='001.hs'
  qfftype=1
  getQFFplot=.FALSE.
  printaverageOK=.FALSE.

! Start reading input
     REWIND(5)
     READ(5,gen,end=1001,err=1002,iostat=ierr)
1002 if(ierr .ne. 0) print*,"gen namelist error!!!",ierr
1001 continue
     REWIND(5)
     READ(5,mol,end=1003,err=1004,iostat=ierr)
1004 if(ierr .ne. 0) print*,"mol namelist error!!!",ierr
1003 continue
     REWIND(5)
     READ(5,vib,end=1005,err=1006,iostat=ierr)
1006 if(ierr .ne. 0) print*,"vib namelist error!!!",ierr
1005 continue
     REWIND(5)
     READ(5,airun,end=1007,err=1008,iostat=ierr)
1008 if(ierr .ne. 0) print*,"airun namelist error!!!",ierr
1007 continue
     REWIND(5)
     READ(5,pes,end=1011,err=1012,iostat=ierr)
1012 if(ierr .ne. 0) print*,"inp namelist error!!!",ierr
1011 continue

!
     IF(debug)PRINT*,"Debugging Mode"
     IF(unittype==2)convert=1.d0
     nDOF=nfree*nMonomer
     nAtom=Nat
     nMR=MR
     if(nMR>nDOF)nMR=nDOF
     ALLOCATE(hrmfreq(Ndof),freq(Ndof))
     if (NDOF>0)CALL get_qff()
     DO n=0,Nmonomer-1
       DO m=1,nFREE
         hrmfreq(m+n*Nfree)=dsqrt(2.d0*Hii(m))
!         hrmfreq(m+n*Nfree)=omega(m)
       ENDDO
     ENDDO
     !print*,hrmfreq
     freq=hrmfreq

     IF (runtyp=='WRT')THEN
       REWIND(5)
       READ(5,inp,end=1009,err=1010,iostat=ierr)
1010  if(ierr .ne. 0) print*,"inp namelist error!!!",ierr
1009 continue
     CALL initFormPES()
     k2=1;
     DO m=1,ndof
       k1=1
       DO i=1,Natom
         IF (m==1)THEN
           atomlabel(i)=label(i)
           xmass(i)=mass(i)
         ENDIF
         DO j=1,3
           IF (m==1)eqbgeo(i,j)=x(k1)
           qcoord(k1,m)=L(k2)
           k2=k2+1
           k1=k1+1
         ENDDO
       ENDDO
     ENDDO !m
   ELSE IF (runtyp=='PES')THEN
     Ncalc=1+6*Ndof
     PRINT*,"1MR Ncalc=",ncalc
     IF(nMR>1)THEN
       Ncalc=Ncalc+Ndof*(Ndof-1)/2*12
       PRINT*,"2MR Ncalc=",ncalc
       IF(nMR>2)THEN
         Ncalc=Ncalc+Ndof*(Ndof-1)*(Ndof-2)/6*8
         PRINT*,"3MR Ncalc=",ncalc
       ENDIF
     ENDIF
     ALLOCATE(EPES(Ncalc))
   ELSE IF (runtyp=='VIB')THEN
     IF(vscf .OR. vci  .or. vvscf .or. vci1)  ALLOCATE(hrmbasis(Ndof))
     IF(vmax_base > 0)vmax=vmax_base
     FCIsize=1
     maxRCIsize=1
     DO n=0,Nmonomer-1
       DO m=1,nFREE
         IF(vscf .OR. vci .OR. vci1 .OR. vvscf)hrmbasis(m+n*Nfree)=vmax(m)
         IF(vci .OR. vmp2 .OR. vscfci)FCIsize=FCIsize*hrmbasis(m)
       ENDDO
     ENDDO
     IF(vscf .OR. vci .OR. vci1 .OR. vvscf) maxbasis=maxval(hrmbasis)
     IF(allzero)nMR=1
     IF(getfunds) ALLOCATE(funds(0:NDOF))

     IF(ateqb)PRINT*,"At equilibrium, Gi=0,Fii=omega*omega/2"
     IF(vscf .or. vvscf)  CALL init_VSCF()
     IF(xvscf) CALL init_XVSCF()
     IF(vmp2 .OR. vscfci .OR. vscfci1)CALL init_modalint()
 !    IF(ssvmp)PRINT*,"State-specific VMP"
     IF(VCI.OR. vscfci .OR. vscfci1 .OR. vci1) call init_VCI()
     IF(VMP2 .OR. Vvscf .OR. virtual .or. vscf)then
       ALLOCATE(VSCFcoefs(ndof,maxbasis,maxbasis))
       ALLOCATE(VSCFenergies(ndof,maxbasis))
       vscfcoefs=0.d0
       VSCFenergies=0.d0
     endif
     IF(vscf .or. xvscf .or. vvscf)PRINT 101,"VSCF threshold (in 1/cm) and max # of iterations: ",scfthresh,maxiter
     IF(runmonomer)PRINT*,'Nmonomer=',Nmonomer
101  FORMAT(a,d8.1,2x,i5)
   ENDIF

   PRINT 100,nDOF,' modes with a', nMR, 'MR QFF using',maxbasis,' HO basis functions for each mode'
100 FORMAT(i3,a,i2,a,i3,a)
   PRINT*,"ZPVE and fundamental frequencies (in 1/cm)"
   Print*,"Harmonic Approximation"
   PRINT*,SUM(hrmfreq)*0.5d0*convert
     DO m=1,Nfree
       PRINT*,hrmfreq(m)*convert
     ENDDO
   IF(debug)PRINT*,"readinput :)"

 END SUBROUTINE readinput
 !
 !SUBROUTINE getVSCFfunds()
 !  USE global
 !  USE constants
 !  USE modVSCF
 !  USE modXVSCF
 !  IMPLICIT NONE
 !  REAL*8 :: energy,groundE
 !  INTEGER::iter,state(nDOF),m
 !  state=0
 !  IF(VSCF)THEN
 !    CALL subVSCF(state,energy,iter)
 !  ELSE
 !    IF(nMR==1)STOP "qVSCF=hrm for 1MR"
 !    CALL subqVSCF(state,energy,iter)
 !    !to print virtual vscf states
 !    IF(virtual)THEN
 !      PRINT*,"Virtual fundamental frequencies in cm-1 for qvscf ground state"
 !      PRINT*,SUM(freq)*0.5d0*convert
 !      DO m=1,ndof
 !        PRINT*,freq(m)*convert
 !      ENDDO
 !    ENDIF
 !     !virtuals print ends
 !  ENDIF
 !  groundE=energy
 !
 !  !  PRINT*,'ground        0', iter, 'iters for', groundE*convert
 !  IF(getFunds)THEN
 !    funds(0)=energy
 !    DO m=1,nDOF
 !      state=0
 !      state(m)=1
 !      IF(VSCF)THEN
 !        CALL subVSCF(state,energy,iter)
 !      ELSE
 !        CALL subqVSCF(state,energy,iter)
 !      ENDIF
 !      funds(m)=energy-funds(0)
 !      PRINT*,'fund', m, iter, 'iters for', energy*convert
 !    ENDDO
 !    !to import excel easily
 !    PRINT*,"All fundamentals in cm-1 for vscf=",vscf
 !    DO m=0,ndof
 !      PRINT*,funds(m)*convert
 !    ENDDO
 !  ENDIF
 !
 !END SUBROUTINE getVSCFfunds
 !
 !SUBROUTINE getVSCFovertones()
 !  USE global
 !  USE constants
 !  USE modVSCF
 !  USE modXVSCF
 !  IMPLICIT NONE
 !  REAL*8 :: energy,groundE
 !  INTEGER::iter,state(nDOF),m
 !  state=0
 !  IF(VSCF)THEN
 !    CALL subVSCF(state,energy,iter)
 !  ELSE
 !    PRINT*
 !    PRINT*,"************************qVSCF*****************"
 !    PRINT*
 !    IF(nMR==1)STOP "qVSCF=hrm for 1MR"
 !    CALL subqVSCF(state,energy,iter)
 !    !to print virtual vscf states
 !    PRINT*,"Virtual fundamental frequencies in cm-1 for qvscf ground state"
 !    PRINT*,SUM(freq)*0.5d0*convert
 !    DO m=1,ndof
 !      PRINT*,freq(m)*convert
 !    ENDDO
 !     !virtuals print ends
 !  ENDIF
 !  groundE=energy
 !  funds(0)=energy
 !  PRINT*,'ground        0', iter, 'iters for', groundE*convert
 !
 !  DO m=1,nDOF
 !    state=0
 !    state(m)=2
 !    IF(VSCF)THEN
 !      CALL subVSCF(state,energy,iter)
 !    ELSE
 !      CALL subqVSCF(state,energy,iter)
 !    ENDIF
 !    funds(m)=energy-funds(0)
 !    PRINT*,'fund', m, iter, 'iters for', energy*convert
 !  ENDDO
 !  !to import excel easily
 !  PRINT*,"All fundamentals in cm-1 for vscf=",vscf
 !  DO m=0,ndof
 !    PRINT*,funds(m)*convert
 !  ENDDO
 !  !call plotQFF()
 !END SUBROUTINE getVSCFovertones

 SUBROUTINE getmodel
   USE global
   USE constants
   USE modQFF,ONLY:model1_QFF
   USE modVSCF
   USE modXVSCF
   USE modVMP2
   IMPLICIT NONE
   REAL*8 :: energy,groundE,energy10,energy01,energy20
   INTEGER::iter,state(nDOF),m
   PRINT*,"model calc"
   CALL model1_QFF()
   state=0
   IF(VSCF)THEN
     CALL subVSCF(state,energy,iter)
   ELSE IF (VMP2)THEN
     CALL subVMP2(state,energy)
   ELSE
     CALL subqVSCF(state,energy,iter)
   ENDIF
   groundE=energy
   ! PRINT*,'ground ', 'is calculated after', iter, 'iterations', groundE
   DO m=1,2
     state=0
     state(m)=1
     IF(VSCF)THEN
       CALL subVSCF(state,energy,iter)
     ELSE IF (VMP2) THEN
       CALL subVMP2(state,energy)
     ELSE
       CALL subqVSCF(state,energy,iter)
     ENDIF
     !   PRINT*,'fundamental', m, 'is calculated after', iter, 'iterations',energy, (energy-groundE)
     IF (m==1)energy10=energy
     IF(m==2)energy01=energy
   ENDDO
   state=0
   state(1)=2
   state(2)=0
   IF(VSCF)THEN
     CALL subVSCF(state,energy,iter)
   ELSE IF (VMP2) THEN
     CALL subVMP2(state,energy)
   ELSE
     CALL subqVSCF(state,energy,iter)
   ENDIF
   energy20=energy
   !PRINT*, state, 'is calculated after', iter, 'iterations',energy, (energy-groundE)
   PRINT*,groundE,energy10,energy01,energy20
   !call myvci
 END SUBROUTINE getmodel

 SUBROUTINE printfunds()
   USE global
   USE modXVSCF,ONLY:subxvscf
   USE modVSCF,ONLY:subvscfvirtual
   USE modvci
   USE modvmp2
   IMPLICIT NONE
   REAL*8 :: energy,groundE
   INTEGER::state(nDOF),m
   CHARACTER::mymethod*10
   IF(debug) PRINT*,"printfunds"
   IF(xvscf)CALL subxvscf()
   if(useXVSCFok)hrmfreq=freq
   IF(vvscf)CALL subvscfvirtual()
   IF(vpt20)call subvpt20()
   if(test)call subvpt1b()
   if (method .ne. '')then
     PRINT*,method
     state=0
     CALL getenergy(state,energy)
     groundE=energy
     PRINT*,energy*convert
     if(getfunds)then
       DO m=1,ndof
         state=0
         state(m)=1
         CALL getenergy(state,energy)
         PRINT*,(energy-groundE)*convert!,energy*convert
       ENDDO
     endif
   endif
   if (vscf)then
     mymethod='VSCF'
     call printresults(mymethod)
   endif
   if(vpt)then
     mymethod='VPT1'
     call printresults(mymethod)
     mymethod='VPT2'
     call printresults(mymethod)
   endif
   if (oneMRP1)then
     mymethod='1MRP1'
     call printresults(mymethod)
   endif
   if (oneMRP2)then
     mymethod='1MRP2'
     call printresults(mymethod)
   endif
   if (vmp2)then
     if(ssVMP) mymethod='ssVMP2'
     call printresults(mymethod)
     ssvmp=.false.
     mymethod='VMP2'
     call printresults(mymethod)
   endif
   ! if(test)call printfunds2('VPT2test')
   IF(vci .or. vscfci .or. vci1 .or. vscfci1)CALL runvci()
   IF(debug) PRINT*,"printfunds :)"
 END SUBROUTINE printfunds

 SUBROUTINE printresults(mymethod)
   USE global
   IMPLICIT NONE
   REAL*8 :: energy,groundE
   INTEGER::state(nDOF),m
   CHARACTER::mymethod*10
   IF(debug) PRINT*,"printfunds2"
   PRINT*,mymethod
   state=0
   method=trim(mymethod)
   CALL getenergy(state,energy)
   groundE=energy
   PRINT*,energy*convert
   if(getfunds)then
     DO m=1,ndof
       state=0
       state(m)=1
       CALL getenergy(state,energy)
       PRINT*,(energy-groundE)*convert!,energy*convert
     ENDDO
   endif
   if(overtone) then
     PRINT*,mymethod," Overtones"
     DO m=1,ndof
       state=0
       state(m)=2
       CALL getenergy(state,energy)
       PRINT*,(energy-groundE)*convert!,energy*convert
     ENDDO
   endif
   IF(debug) PRINT*,"printfunds2 :)"
 END SUBROUTINE printresults


 SUBROUTINE getenergy(state,energy)
   USE global
   USE modVSCF!old
   USE modXVSCF
   USE modVMP2
   IMPLICIT NONE
   INTEGER,INTENT(IN)::state(nDOF)
   REAL*8,INTENT(OUT):: energy
   energy=0.d0
   SELECT CASE(TRIM(method))
     CASE('VSCF')
       CALL subVSCF(state,energy)
     CASE('qVSCF')
       CALL subqVSCF(state,energy)
     CASE('VMP1')
       CALL subVMP1(state,energy)
     CASE('VMP2')
       CALL subVMP2(state,energy)
     CASE('ssVMP2')
       CALL subVMP2(state,energy)
     CASE('VMP2a')
       CALL subVMP2a(state,energy)
     CASE('VMP2b')
       CALL subVMP2b(state,energy)
     CASE('VMP2c')
       CALL subVMP2c(state,energy)
     CASE('VMP2f')
       CALL subVMP2f(state,energy)
     CASE('VPT1')
       CALL subVPT1(state,energy)
     CASE('VPT2')
       CALL subVPT2(state,energy)
     !   CASE('XVPT2')
     CASE('VPT2test')
       CALL subVPT2test(state,energy)
     !   CASE('XVPT2')
     !     CALL subXVPT2(energy)
     CASE('VPT2a')
       CALL subVPT2a(state,energy)
     CASE('qVMP1')
       CALL subqVMP1(state,energy)
     CASE('qVMP2')
       CALL subqVMP2(state,energy)
     CASE('1MRP1')
       CALL sub1MRP1(state,energy)
     CASE('1MRP2')
       CALL sub1MRP2(state,energy)
     CASE DEFAULT
     Continue
 END SELECT
 END SUBROUTINE getenergy

 !SUBROUTINE onemode(mode,energies,eigenfuns)
 !  USE global
 !  USE lapack
 !  USE yazar
 !  USE constants
 !  USE modQFF
 !  USE integral
 !  IMPLICIT NONE
 !  INTEGER,INTENT(IN):: mode
 !  INTEGER:: i,j,diff
 !  REAL*8:: tmp,ham(hrmbasis(mode),hrmbasis(mode)),energies(maxbasis)
 !  REAL*8, OPTIONAL,INTENT(OUT)::eigenfuns(hrmbasis(mode),hrmbasis(mode))
 !  Ham=0.d0
 !  PRINT*,'onemode'
 !  DO i=0, hrmbasis(mode)-1
 !     Ham(i+1,i+1)=(i + 0.5) * hrmfreq(mode) * 0.5 &!kinetic term
 !          + Hii(mode) * harmint(i, mode, 2, 0)&
 !          + Qiiii(mode) * harmint(i, mode, 4, 0)
 !     DO j=i+1,hrmbasis(mode)-1!
 !        IF ( (j-i) > 4) EXIT !break;
 !        diff = j - i;
 !        tmp=0.d0
 !        IF (diff .EQ. 2) tmp= tmp-dSQRT(DBLE(j * (j - 1))) * hrmfreq(mode) * 0.25 !!kinetic term
 !        tmp= tmp + harmint(j, mode, 1, diff) * gi(mode);
 !        tmp= tmp + harmint(j, mode, 2, diff) * Hii(mode);
 !        tmp= tmp + harmint(j, mode, 3, diff) * Ciii(mode);
 !        tmp= tmp + harmint(j, mode, 4, diff) * Qiiii(mode);
 !        Ham(i+1,j+1)=tmp
 !     ENDDO!!end for j
 !  ENDDO!!end for i
 !
 !  CALL diag(hrmbasis(mode),ham,energies)
 !  !if(present(eigenfuns))eigenfuns=ham !????????
 !  RETURN
 !END SUBROUTINE onemode
 !
 !SUBROUTINE qonemr()
 !  USE global
 !  USE constants
 !  IMPLICIT NONE
 !  INTEGER :: mode
 !  REAL*8 :: groundE,energies(maxbasis),oneMRfunds(0:nDOF)
 !  REAL*8::HOfunds(0:nDOF),qVSCFfunds(0:nDOF)
 !
 !  PRINT*,"********************qOne-MR"
 !  groundE=0.d0
 !  DO mode=1,nDOF
 !     CALL onemode(mode,energies)
 !     groundE=groundE+energies(1)
 !     oneMRfunds(mode)=energies(2)-energies(1)
 !     !print*,'mode ',mode,(oneMRfunds(mode))*convert
 !  ENDDO
 !  oneMRfunds(0)=groundE
 !
 !  CALL getVSCFfunds(qVSCFfunds)
 !
 !  HOfunds(0)=SUM(hrmfreq)/2.d0
 !  HOfunds(1:ndof)=hrmfreq
 !  DO mode=0,nDOF
 !     PRINT*,'HO ',mode,(HOfunds(mode))*convert
 !     PRINT*,'1MR ',mode,(oneMRfunds(mode))*convert
 !     PRINT*,'qVSCF ',mode,(qVSCFfunds(mode))*convert
 !     PRINT*,'q1MR ',mode,(oneMRfunds(mode)-HOfunds(mode)+qVSCFfunds(mode))*convert
 !  ENDDO
 !
 !END SUBROUTINE qonemr


 !
 PROGRAM main
   USE global
   USE modXVSCF
   USE modVSCF
   USE modQFF
   USE constants
   USE modtimer
   USE modformPES
   USE modVMP2
   IMPLICIT NONE
   REAL(8)::wallstart,cpustart,wallend,cpuend
   PRINT*,"MaVi Mavi MaSMaVi v16.01.28 "
   CALL readinput()
   !IF (runtyp=='VIB') CALL get_qff()
!   IF (ssVSCF)CALL printssVSCf

   !CALL formPES()
   !CALL form_QFF()
   !CALL write_QFF()
   ! CALL read_Sindo()
   ! CALL put_coefs()!only for VSCF
   ! CALL make_quartic()
   !  CALL make_harmonic()

   !CALL convert2atomicunits()
   ! CALL write_qff()
   ! CALL groundstate()
   !   CALL model1_QFF()
   !    CALL make_harmonic()
   ! CALL qonemr()
   CALL wall_and_cpu_time(wallstart,cpustart)
   !CALL getVSCFfunds()
   !  CALL getVSCFovertones()

   CALL printfunds()
   IF (testgreen) CALL subtestgreen()
   IF (testgreen) CALL subtestgreen4()
   !IF (testgreen) CALL subtestgreen3()
   !CALL printfreq()
   ! CALL myvci()
   CALL wall_and_cpu_time(wallend,cpuend)
   PRINT*,"wall=",wallend-wallstart,"cpu=",cpuend-cpustart
   STOP
 END PROGRAM main


!
!MODULE modqVMP2
!  IMPLICIT NONE
!CONTAINS
!
!  SUBROUTINE subqVMP2
!  USE CONSTANTS
!  USE global
!  USE integral
!  USE modQFF
!  USE modStates
!  USE modXVSCF
!    IMPLICIT NONE
!    INTEGER::thestate(ndof),m
!    REAL*8::energy,energy0
!    thestate=0
!    PRINT*,"qVMP2 starts with maxNstate=",maxNstate
!    CALL EqVMP2(thestate,energy0)
!    PRINT*,"final E",energy0*convert
!    DO m=1,ndof
!       thestate=0
!       thestate(m)=1
!       CALL EqVMP2(thestate,energy)
!       PRINT*,"final E",(energy-energy0)*convert
!    ENDDO
!  END SUBROUTINE subqVMP2
!
!  SUBROUTINE EqVMP2(thestate,energy)
!  USE CONSTANTS
!  USE global
!  USE integral
!  USE modQFF
!  USE modStates
!  USE modXVSCF
!    IMPLICIT NONE
!    INTEGER::i,m1,m2,m3,m4,m
!    INTEGER::thestate(ndof),state(ndof),iter,test
!    REAL*8::theenergy,energy,tmp,pot,Ediff(maxNstate)
!    freq=hrmfreq
!    CALL subqVSCF(thestate,theenergy,iter)
!    freq=hrmfreq
!    !  CALL getEdiff(thestate,Ediff)
!    tmp=0.d0
!    ! print*,Ediff
!    ! PRINT*,"Efreq",hrmfreq,freq
!    DO i=1,maxnstate
!       test=0
!       CALL statelabel(i,state)
!       DO m=1,ndof
!          IF(state(m)==thestate(m))test=test+1
!       ENDDO
!       IF (test==ndof)CYCLE!print*,"******same****************"
!       CALL subqVSCF(state,energy,iter)
!       freq=hrmfreq
!       IF(isNAN(energy))CYCLE
!       ediff(i)=theenergy-energy
!       ! print*," state, energy",state,energy,iter,ediff(i)
!       DO m1=1, ndof
!          pot=qvmpot1(m1,thestate(m1),state(m1))
!          ! PRINT*,"qvmpot1",pot
!          tmp=tmp+pot*pot/Ediff(i)
!          IF (nmr>1)THEN
!             DO m2=1, ndof
!                IF (m2==m1)CYCLE
!                pot=qvmpot2(m1,m2,thestate,state)
!                !     PRINT*,"qvmpot2",pot
!                tmp=tmp+pot*pot/Ediff(i)
!                IF (nmr>2)THEN
!                   DO m3=1, ndof
!                      IF (m1==m3 .OR. m2==m3)CYCLE
!                      pot=qvmpot3(m1,m2,m3,thestate,state)
!                      !          PRINT*,"qvmpot3",pot
!                      tmp=tmp+pot*pot/Ediff(i)
!                      IF (nmr>3)THEN
!                         DO m4=1, ndof
!                            IF (m1==m4 .OR. m2==m4 .OR. m3==m4)CYCLE
!                            pot=qvmpot4(m1,m2,m3,m4,thestate,state)
!                            tmp=tmp+pot*pot/Ediff(i)
!                         ENDDO !m4
!                      ENDIF !nmr3
!                   ENDDO !m3
!                ENDIF !nmr2
!             ENDDO !m2
!          ENDIF !nmr1
!       ENDDO !m1
!    ENDDO !i
!    energy=tmp
!    PRINT*,"qVMP2 ends for thestate",thestate, energy*convert,theenergy*convert
!    energy=energy+theenergy
!  END SUBROUTINE EqVMP2
!
!  SUBROUTINE getEdiff(thestate,Ediff)
!  USE CONSTANTS
!  USE global
!  USE integral
!  USE modQFF
!  USE modStates
!  USE modXVSCF
!    IMPLICIT NONE
!    REAL*8::Ediff(maxNstate),energy,theenergy
!    INTEGER::i,thestate(ndof),state(ndof),iter
!    CALL subqVSCF(thestate,theenergy,iter)
!    PRINT*,"the state, energy",thestate,theenergy
!    DO i=1,maxNstate
!       CALL statelabel(i,state)
!       CALL subqVSCF(state,energy,iter)
!       freq=hrmfreq
!       PRINT*," state, energy",state,energy,iter
!       ediff(i)=theenergy-energy
!    ENDDO
!  END SUBROUTINE getEdiff
!
!  REAL*8 FUNCTION qvmpot1(mode,leftlevel,rightlevel)
!  USE CONSTANTS
!  USE global
!  USE integral
!  USE modQFF
!  USE modStates
! ! USE modXVSCF
!    IMPLICIT NONE
!    INTEGER::mode,leftlevel,rightlevel,n,m,diff
!    REAL*8::tmp
!    n=MAX(leftlevel,rightlevel)
!    diff=ABS(leftlevel-rightlevel)
!    tmp=0.d0
!    tmp=tmp+Gi(mode)*harmint(n,mode,1,diff)
!    IF (nMR>1)THEN
!       DO m=1,ndof
!          IF (m==mode)CYCLE
!          tmp=tmp+0.5d0*Ciij(m,mode)*harmint(n,m,2,diff) &
!               *harmint(n,mode,1,diff)
!       ENDDO
!    ENDIF
!    qvmpot1=tmp
!    RETURN
!  END FUNCTION qvmpot1
!
!  REAL*8 FUNCTION qvmpot2(mode1,mode2,leftstate,rightstate)
!    USE CONSTANTS
!  USE global
!  USE integral
!  USE modQFF
!  USE modStates
!  USE modXVSCF
!    IMPLICIT NONE
!    INTEGER::mode1,mode2,n1,n2,m,diff1,diff2
!    INTEGER,INTENT(in)::leftstate(ndof),rightstate(ndof)
!    REAL*8::tmp
!    n1=MAX(leftstate(mode1),rightstate(mode1))
!    n2=MAX(leftstate(mode2),rightstate(mode2))
!    diff1=ABS(leftstate(mode1)-rightstate(mode1))
!    diff2=ABS(leftstate(mode2)-rightstate(mode2))
!    tmp=0.d0
!    tmp=tmp+Hij(mode1,mode2)*harmint(n1,mode1,1,diff1)*harmint(n2,mode2,1,diff2)
!    IF (nMR>2)THEN
!       DO m=1,ndof
!          IF (m==mode1 .OR. m==mode2)CYCLE
!          tmp=tmp+0.5d0*Qiijk(m,mode1,mode2)*harmint(n1,mode1,1,diff1)*harmint(n2,mode2,1,diff2)&
!               *harmint(MAX(leftstate(m),rightstate(m)),m,2,ABS(leftstate(m)-rightstate(m)))
!       ENDDO
!    ENDIF
!    qvmpot2=tmp
!    RETURN
!  END FUNCTION qvmpot2
!
!  REAL*8 FUNCTION qvmpot3(mode1,mode2,mode3,leftstate,rightstate)
!  USE CONSTANTS
!  USE global
!  USE integral
!  USE modQFF
!  USE modStates
!  USE modXVSCF
!    IMPLICIT NONE
!    INTEGER::mode1,mode2,mode3,n1,n2,n3,diff1,diff2,diff3
!    INTEGER,INTENT(in)::leftstate(ndof),rightstate(ndof)
!    REAL*8::tmp
!    n1=MAX(leftstate(mode1),rightstate(mode1))
!    n2=MAX(leftstate(mode2),rightstate(mode2))
!    n3=MAX(leftstate(mode3),rightstate(mode3))
!    diff1=ABS(leftstate(mode1)-rightstate(mode1))
!    diff2=ABS(leftstate(mode2)-rightstate(mode2))
!    diff3=ABS(leftstate(mode3)-rightstate(mode3))
!    tmp=0.d0
!    tmp=tmp+Cijk(mode1,mode2,mode3)*harmint(n1,mode1,1,diff1)*harmint(n2,mode2,1,diff2)*harmint(n3,mode3,1,diff3)
!    qvmpot3=tmp
!    RETURN
!  END FUNCTION qvmpot3
!
!  REAL*8 FUNCTION qvmpot4(mode1,mode2,mode3,mode4,leftstate,rightstate)
!   USE CONSTANTS
!  USE global
!  USE integral
!  USE modQFF
!  USE modStates
!  USE modXVSCF
!    IMPLICIT NONE
!    INTEGER::mode1,mode2,mode3,mode4,n1,n2,n3,n4,diff1,diff2,diff3,diff4
!    INTEGER,INTENT(in)::leftstate(ndof),rightstate(ndof)
!    REAL*8::tmp
!    n1=MAX(leftstate(mode1),rightstate(mode1))
!    n2=MAX(leftstate(mode2),rightstate(mode2))
!    n3=MAX(leftstate(mode3),rightstate(mode3))
!    n4=MAX(leftstate(mode4),rightstate(mode4))
!    diff1=ABS(leftstate(mode1)-rightstate(mode1))
!    diff2=ABS(leftstate(mode2)-rightstate(mode2))
!    diff3=ABS(leftstate(mode3)-rightstate(mode3))
!    diff4=ABS(leftstate(mode4)-rightstate(mode4))
!    tmp=0.d0
!    !tmp=tmp+Qijkl(mode1,mode2,mode3,mode4)*harmint(n1,mode1,1,diff1)*harmint(n2,mode2,1,diff2)*harmint(n3,mode3,1,diff3)&
!    !*harmint(n4,mode4,1,diff4)
!    qvmpot4=tmp
!    RETURN
!  END FUNCTION qvmpot4
!
!END MODULE modqVMP2
