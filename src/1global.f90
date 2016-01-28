MODULE global
  IMPLICIT NONE
  INTEGER :: maxcell,maxk,Ndof,Ndim,Natom,Ncell,N_k,nMR,FCIsize,maxNstate,maxFF,Ncalc,maxiter,Nmonomer
  INTEGER::unittype,maxbasis !unittype=1
  INTEGER, ALLOCATABLE::hrmbasis(:)
  REAL*8, ALLOCATABLE ::xmass(:),k_freq(:,:),hrmfreq(:),NC(:,:),freq(:),funds(:)
  REAL*8, ALLOCATABLE ::vscfenergies(:,:),vscfcoefs(:,:,:)
  REAL*8, ALLOCATABLE ::uZero(:),uOne(:),uTwo(:),uThree(:),uFour(:)
!  COMPLEX*16, ALLOCATABLE ::k_coef(:,:,:),cck_coef(:,:,:)
  REAL*8  :: geo(100),Eref
  REAL*8::delx,delq
  REAL*8::inttol,scfthresh,convert,denomcut
  LOGICAL::vscf,vci,vpt,vav,vmp2,xvscf,Vvscf,ssvscf,ssvmp,vscfci,VCI1,VSCFCI1,oneMRP1,oneMRP2
  LOGICAL::summary,storepot,runmonomer,getfunds,rightstate,benchmark,checkEtotal,overtone
  LOGICAL::virtual,getfreq,getwfn,storeint,debug,testgreen,useXVSCFok,vpt20,test
  CHARACTER::runtyp*3,method*10

END MODULE global
