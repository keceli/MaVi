MODULE modvci
  IMPLICIT NONE
  INTEGER::ncup,maxsum,maxexall,maxEnergy,maxRCIsize,RCIsize,Nstate
  INTEGER,ALLOCATABLE::snumvec(:)
  REAL*8::Eground,minweight
  LOGICAL::restricted,VCIdebug,VCIprint
CONTAINS

  SUBROUTINE init_VCI
    USE global,ONLY:FCIsize,debug,ndof,maxbasis,hrmbasis,vscfci,vci1,vscfci1,vci
    USE modalintegral,ONLY:init_modalint
    IMPLICIT NONE
    integer::m,CIsize
    IF(debug)PRINT*,"init_VCI"
    Nstate=FCIsize
    ncup=ndof
    maxsum=ndof*maxbasis
    maxexall=maxbasis
    maxenergy=0
    minweight=0.1
    restricted=.FALSE.
    VCIdebug=.false.
    VCIprint=.false.
    Eground=1.d10
    maxRCIsize=1
    CALL read_VCI()
    !if(vscfci .or. vscfci1)call init_modalint()
    !PRINT*,'Full CI size= ',FCIsize
    IF(Nstate>FCIsize)Nstate=FCIsize
    IF(Ncup>Ndof)Ncup=Ndof
    IF(maxexall>maxbasis)maxexall=maxbasis
    OPEN(UNIT=17,FILE='vci.mavi',STATUS='UNKNOWN',ACTION='WRITE')
    IF(VCIprint)OPEN(UNIT=27,FILE='allvci.mavi',STATUS='UNKNOWN',ACTION='WRITE')
    WRITE(*,100)'FVCI: size= ',FCIsize
100 FORMAT(a,i20)
    IF(restricted .and. (vci .or. vscfci))THEN
      do m=1,ndof
        maxRCIsize=maxRCIsize*min(hrmbasis(m),maxsum+1)
      enddo
      CALL formsnumvec()
      IF(Nstate>RCIsize)Nstate=RCIsize
      ! PRINT*,'maxsum= ',maxsum,'ncup=',ncup,
      WRITE(*,200)'RVCI: ncup,maxsum,maxexall,maxEnergy,size= ',ncup,maxsum,maxexall,maxenergy,RCIsize
200   FORMAT(a,4i4,i8)
     ! PRINT*,'Max restricted CI size=',maxRCIsize
     ! PRINT*,'Restricted CI size=',RCIsize
    ENDIF
    if(vci1 .or. vscfci1)then
    CIsize=1+(maxexall-1)*ndof
    Nstate=CIsize
          WRITE(*,300)'CI1: maxexall,size= ',maxexall,CIsize
300   FORMAT(a,4i4,i8)
    endif
    !    redVCIorder=snumvecsize
    IF(debug)PRINT*,"init_VCI :)"
  END SUBROUTINE init_VCI

  SUBROUTINE read_VCI
    use global,only:debug
    IMPLICIT NONE
    integer::ierr
    !ncup: number of mode couplings
    !maxsum: max number of sum of excitation numbers
    !restricted:logical controls ci space
    !Nstate:
    NAMELIST /vci/ncup,maxsum,restricted,Nstate,maxexall,maxenergy,VCIdebug,VCIprint,minweight
    if(debug) print*,"read_vci"
    REWIND(5)
    READ(5,vci,end=1007,err=1008,iostat=ierr)
1008 if(ierr .ne. 0) print*,"vci namelist error!!!",ierr
1007 continue
    if(debug) print*,"read_vci :)"
  END SUBROUTINE read_VCI

  SUBROUTINE runVCI
    use global
    IMPLICIT NONE
    IF(vci)CALL getvci()
    IF(vci1)CALL getvci1()
        IF(vscfci)CALL getvscfci()
    IF(vscfci1)CALL getvscfci1()

  END SUBROUTINE runVCI

  !  SUBROUTINE redstatelabel(nthstate,ket,newn)
  !    USE global
  !    IMPLICIT NONE
  !    INTEGER::m,state,nthstate,ket(ndof),n
  !    s = nthState;
  !    Do
  !    DO m=1,ndof
  !       n=hrmbasis(m)
  !       ket(m) = MOD(s,n);
  !       s = s / n;
  !    ENDDO
  !    if((sum(ket).le.maxsum) )exit
  !    s=nthState+1
  !    if(s .ge. vciorder )exit
  !    ENDDO
  !    RETURN
  !  END SUBROUTINE redstatelabel
  !
  !subroutine restrictsize
  !use global
  !implicit none
  !integer::i,j,k,m
  !redVCIorder=1
  !do m=1,maxnexcited
  !redvciorder=redvciorder
  !enddo
  !end subroutine

  INTEGER FUNCTION snumvecsize()
    USE global
    USE modstates
    IMPLICIT NONE
    INTEGER::i,n,state(ndof),nexcited,statesum
    n=0
    DO i=0,FCIsize-1
      CALL statelabel(i,state)
      CALL countnonzeros(state,nexcited)
      statesum=SUM(state)
      IF( (statesum .LE. maxsum) .AND. (nexcited) .LE. ncup)n=n+1
    ENDDO
    snumvecsize=n
    RETURN
  END FUNCTION snumvecsize

  SUBROUTINE formsnumvec()
    USE global
    USE modstates
    IMPLICIT NONE
    INTEGER::i,j,state(ndof),limit
    logical::passit
    IF(debug)Print*,"formsnumvec"
    j=0
    limit=min(maxRCIsize,FCIsize)-1
    DO i=0,limit
      CALL RCIstatecheck(i,state,passit)
      IF (passit)cycle
      j=j+1
    ENDDO
    RCIsize=j
    IF(debug)print*,"Allocating snumvec",maxRCIsize,FCIsize
    !      ALLOCATE(snumvec(min(maxRCIsize,FCIsize)))
    ALLOCATE(snumvec(RCIsize))
    j=0
    DO i=0,limit
      CALL RCIstatecheck(i,state,passit)
      IF (passit)cycle
      j=j+1
      snumvec(j)=i
    ENDDO

    IF(debug)Print*,"formsnumvec :)"
  END SUBROUTINE formsnumvec

  SUBROUTINE RCIstatelabel(nthstate,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m,stateindex,nthstate,ket(ndof),n
    stateindex = nthState;
    DO m=1,ndof
      n=min(maxsum+1,hrmbasis(m))
      ket(m) = MOD(stateindex,n);
      stateindex = stateindex / n;
    ENDDO
    RETURN
  END SUBROUTINE RCIstatelabel

  SUBROUTINE CI1statelabel(number,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m,n,number,ket(ndof)
    ket=0
    if(number==0)return
    m = number / (maxexall-1) + 1
    n = MOD(number,maxexall-1)
    if(n==0)then
      ket(m-1)=maxexall-1
    else
      ket(m)=n
    endif
    RETURN
  END SUBROUTINE CI1statelabel

  SUBROUTINE RCIstatecheck(nthstate,ket,passit)
    USE global
    IMPLICIT NONE
    INTEGER::m,stateindex,nthstate,ket(ndof),n
    real*8::energy
    logical::passit
    ket=0
    passit=.false.
    stateindex = nthState;
    DO m=1,ndof
      n=min(maxsum+1,hrmbasis(m))
      ket(m) = MOD(stateindex,n);
      if( (count(ket>0) > ncup) .or. (sum(ket) > maxsum) .or. (ANY(MASK=ket > maxexall)))then
        passit=.true.
        return
      else if (maxEnergy .GT. 0) then
        energy=funHOenergy(ket)
        if(energy.GT.maxEnergy*minval(hrmfreq)) then
          passit=.true.
          return
        endif
      endif
      stateindex = stateindex / n;
    ENDDO
    RETURN
  END SUBROUTINE RCIstatecheck

  real*8 function funHOenergy(ket)
    USE global,only:ndof,freq
    IMPLICIT NONE
    INTEGER::m,ket(ndof)
    funHOenergy=0.d0
    DO m=1,ndof
      funHOenergy=funHOenergy+(ket(m)+0.5d0)*freq(m)
    ENDDO
    RETURN
  END Function funHOenergy

  SUBROUTINE formsnumvecold()
    USE global
    USE modstates
    IMPLICIT NONE
    INTEGER::i,j,state(ndof),nexcited,statesum
    j=0
    DO i=0,FCIsize-1
      CALL statelabel(i,state)
      CALL countnonzeros(state,nexcited)
      statesum=SUM(state)
      IF( (statesum .LE. maxsum) .AND. (nexcited .LE. ncup)  )THEN
        j=j+1
        snumvec(j)=i
      ENDIF
    ENDDO
    RCIsize=j

    !print*,maxcoef,"for state",state
  END SUBROUTINE formsnumvecold

  SUBROUTINE getvci
    USE lapack
    USE yazar
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,j,m,m1,m2,m3,m4,n,n1,n2,n3,n4,size,maxndiffer
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer,diff1,diff2,diff3,diff4
    REAL*8:: tmp
    REAL*8,ALLOCATABLE::ham(:,:),energies(:)
    IF(debug)PRINT*,"getVCI"
        Eground=1.d10
    size=FCIsize
    maxndiffer=nMR
    IF(VCIprint)write(27,'(a)')'HOVCI'
    write(17,'(a)')'HOVCI'
    IF(restricted)THEN
      size=RCIsize
    ENDIF
    ALLOCATE(ham(size,size),energies(size))
    ham=0.d0
    !call countstates()
    if(vcidebug)PRINT*,"HO VCI Hamiltonian initialized with size", size
    DO i=0,size-1
      tmp=0.d0
      IF(restricted)THEN
        CALL RCIstatelabel(snumvec(i+1),ket)
      ELSE
        CALL statelabel(i,ket)
      ENDIF
      !    print*,"ket",ket
      DO m=1,ndof
        n = ket(m);
        tmp=tmp+ Hii(m) * harmint(n, m, 2, 0) + Qiiii(m) * harmint(n, m, 4, 0);
        tmp=tmp+ (ket(m) + 0.5d0) * hrmfreq(m) * 0.5d0;
      ENDDO
      IF (nMR > 1)THEN
        DO m1=1,nDOF
          n1 = ket(m1);
          DO m2=m1+1, nDOF
            n2 = ket(m2);
            tmp=tmp+ Qiijj(m1,m2) * harmint(n1, m1, 2, 0) * harmint(n2, m2, 2, 0);
          ENDDO!m2
        ENDDO !m1
      ENDIF !nMR>1
      Ham(i+1,i+1)=tmp
      DO j=i+1,size-1
        tmp=0.d0
        IF(restricted)THEN
          CALL RCIstatelabel(snumvec(j+1),bra)
        ELSE
          CALL statelabel(j,bra)
        ENDIF
        !print*,j,"bra",bra
        CALL findiffer(bra,ket,diffvec,ndiffer)
        IF (nDiffer > maxnDiffer) THEN
          CYCLE
        ELSE IF (nDiffer == 1)THEN
          diff1 = ABS(bra(diffVec(1)) - ket(diffVec(1)));
          n1 = MAX(bra(diffVec(1)), ket(diffVec(1)));
          tmp=tmp+ harmint(n1, diffVec(1), 1, diff1) * Gi(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 2, diff1) * Hii(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 3, diff1) * Ciii(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 4, diff1) * Qiiii(diffVec(1));
          IF (diff1 == 2) tmp=tmp - DSQRT(DBLE(n1 * (n1 - 1))) * hrmfreq(diffVec(1)) * 0.25!kinetic term
          IF (nMR > 1)THEN
            DO m=1,ndof
              IF (m == diffVec(1))CYCLE
              n2 = MAX(bra(m), ket(m));
              diff2 = ABS(bra(m) - ket(m));
              tmp=tmp+ Qiijj(m,diffVec(1)) * harmint(n2, m, 2, diff2) * harmint(n1,diffVec(1), 2, diff1);
              tmp=tmp+ Ciij(m,diffVec(1)) * harmint(n2, m, 2, diff2) * harmint(n1,diffVec(1), 1, diff1);
            ENDDO !m1
          ENDIF!nmr>1
        ELSE IF ((nDiffer == 2)) THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          diff1 = ABS(bra(m1) - ket(m1));
          diff2 = ABS(bra(m2) - ket(m2));
          n1 = MAX(bra(m1), ket(m1));
          n2 = MAX(bra(m2), ket(m2));
          tmp=tmp+ Hij(m1,m2) * harmint(n1, m1, 1, diff1) * harmint(n2, m2, 1, diff2);
          tmp=tmp+ Qiijj(m1,m2) * harmint(n1, m1, 2, diff1)* harmint(n2, m2, 2, diff2);
          tmp=tmp+ Ciij(m1,m2) * harmint(n1, m1, 2, diff1) * harmint(n2, m2, 1, diff2);
          tmp=tmp+ Qiiij(m1,m2) * harmint(n1, m1, 3, diff1)* harmint(n2, m2, 1, diff2);
          tmp=tmp+ Ciij(m2,m1) * harmint(n1, m1, 1, diff1) * harmint(n2, m2, 2, diff2);
          tmp=tmp+ Qiiij(m2,m1) * harmint(n1, m1, 1, diff1)* harmint(n2, m2, 3, diff2);
          IF (nMR > 2)THEN
            DO m=1,ndof
              IF (m == m1 .OR. m == m2)CYCLE;
              n = bra(m);
              tmp=tmp+ Qiijk(m,m1,m2) * harmint(n, m, 2, 0) * harmint(n1, m1, 1, diff1)&
              * harmint(n2, m2, 1, diff2);
            ENDDO
          ENDIF

        ELSE IF ((nDiffer == 3))THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          m3 = diffVec(3);
          diff1 = ABS(bra(m1) - ket(m1));
          diff2 = ABS(bra(m2) - ket(m2));
          diff3 = ABS(bra(m3) - ket(m3));
          n1 = MAX(bra(m1), ket(m1));
          n2 = MAX(bra(m2), ket(m2));
          n3 = MAX(bra(m3), ket(m3));
          tmp=tmp+ Cijk(m1,m2,m3) * harmint(n1, m1, 1, diff1)* harmint(n2, m2, 1, diff2)&
          * harmint(n3, m3, 1, diff3);
          tmp=tmp+ Qiijk(m1,m2,m3) * harmint(n1, m1, 2, diff1) * harmint(n2, m2, 1,diff2)&
          * harmint(n3, m3, 1, diff3);
          tmp=tmp+ Qiijk(m2,m1,m3) * harmint(n2, m2, 2, diff2) * harmint(n1, m1, 1,diff1)&
          * harmint(n3, m3, 1, diff3);
          tmp=tmp+ Qiijk(m3,m2,m1) * harmint(n3, m3, 2, diff3) * harmint(n2, m2, 1,diff2)&
          * harmint(n1, m1, 1, diff1);

        ELSE IF ((nDiffer == 4))THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          m3 = diffVec(3);
          m4 = diffVec(4);
          diff1 = ABS(bra(m1) - ket(m1));
          diff2 = ABS(bra(m2) - ket(m2));
          diff3 = ABS(bra(m3) - ket(m3));
          diff4 = ABS(bra(m4) - ket(m4));
          n1 = MAX(bra(m1), ket(m1));
          n2 = MAX(bra(m2), ket(m2));
          n3 = MAX(bra(m3), ket(m3));
          n4 = MAX(bra(m4), ket(m4));
          tmp=tmp+ Qijkl(m1,m2,m3,m4) * harmint(n1, m1, 1, diff1)* harmint(n2, m2, 1, diff2)&
          * harmint(n3, m3, 1, diff3)* harmint(n4, m4, 1, diff4);
        ENDIF
        ham(i+1,j+1)=tmp
      ENDDO!end for j
    ENDDO!end for i
    !  CALL writesquare(vciorder,ham*convert)
    if(vcidebug)PRINT*,"diagonilazation..."
    CALL diag(size,ham,energies)
    !  CALL writesquare(vciorder,ham)
    PRINT*,'HOVCI'
    !      PRINT*,energies*convert
    DO i=1,Nstate
      CALL printstate(size,ham(:,i),energies(i)*convert)
      IF(VCIprint)CALL printall(size,ham(:,i),energies(i)*convert)
    ENDDO
    Eground=1.d10
    IF(debug)PRINT*,"getVCI :)"
  END SUBROUTINE getvci

  SUBROUTINE getvci1
    USE lapack
    USE yazar
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,j,m,m1,m2,n,n1,n2,CIsize,maxndiffer
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer,diff1,diff2
    REAL*8:: tmp
    REAL*8,ALLOCATABLE::ham(:,:),energies(:)
    IF(debug)PRINT*,"getVCI(1)"
        Eground=1.d10

    CIsize=1+(maxexall-1)*ndof
    ! print*,CIsize
    maxndiffer=min(nMR,2)
    IF(VCIprint)write(27,'(a)')'HOVCI1'
    write(17,'(a)')'HOVCI(1)'

    ALLOCATE(ham(CIsize,CIsize),energies(CIsize))
    ham=0.d0
    !call countstates()
    if(vcidebug)PRINT*,"HO VCI(1) Hamiltonian initialized with size", CIsize
    DO i=0,CIsize-1
      tmp=0.d0
      CALL CI1statelabel(i,ket)
      !       print*,"ket",i,ket
      DO m=1,ndof
        n = ket(m);
        tmp=tmp+ Hii(m) * harmint(n, m, 2, 0) + Qiiii(m) * harmint(n, m, 4, 0);
        tmp=tmp+ (ket(m) + 0.5d0) * hrmfreq(m) * 0.5d0;
      ENDDO
      IF (nMR > 1)THEN
        DO m1=1,nDOF
          n1 = ket(m1);
          DO m2=m1+1, nDOF
            n2 = ket(m2);
            tmp=tmp+ Qiijj(m1,m2) * harmint(n1, m1, 2, 0) * harmint(n2, m2, 2, 0);
          ENDDO!m2
        ENDDO !m1
      ENDIF !nMR>1
      Ham(i+1,i+1)=tmp
      DO j=i+1,CIsize-1
        tmp=0.d0
        CALL CI1statelabel(j,bra)
        !print*,j,"bra",bra
        CALL findiffer(bra,ket,diffvec,ndiffer)
        IF (nDiffer > maxnDiffer) THEN
          CYCLE
        ELSE IF (nDiffer == 1)THEN
          diff1 = ABS(bra(diffVec(1)) - ket(diffVec(1)));
          n1 = MAX(bra(diffVec(1)), ket(diffVec(1)));
          tmp=tmp+ harmint(n1, diffVec(1), 1, diff1) * Gi(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 2, diff1) * Hii(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 3, diff1) * Ciii(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 4, diff1) * Qiiii(diffVec(1));
          IF (diff1 == 2) tmp=tmp - DSQRT(DBLE(n1 * (n1 - 1))) * hrmfreq(diffVec(1)) * 0.25!kinetic term
          IF (nMR > 1)THEN
            DO m=1,ndof
              IF (m == diffVec(1))CYCLE
              n2 = MAX(bra(m), ket(m));
              diff2 = ABS(bra(m) - ket(m));
              tmp=tmp+ Qiijj(m,diffVec(1)) * harmint(n2, m, 2, diff2) * harmint(n1,diffVec(1), 2, diff1);
              tmp=tmp+ Ciij(m,diffVec(1)) * harmint(n2, m, 2, diff2) * harmint(n1,diffVec(1), 1, diff1);
            ENDDO !m1
          ENDIF!nmr>1
        ELSE IF ((nDiffer == 2)) THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          diff1 = ABS(bra(m1) - ket(m1));
          diff2 = ABS(bra(m2) - ket(m2));
          n1 = MAX(bra(m1), ket(m1));
          n2 = MAX(bra(m2), ket(m2));
          tmp=tmp+ Hij(m1,m2) * harmint(n1, m1, 1, diff1) * harmint(n2, m2, 1, diff2);
          tmp=tmp+ Qiijj(m1,m2) * harmint(n1, m1, 2, diff1)* harmint(n2, m2, 2, diff2);
          tmp=tmp+ Ciij(m1,m2) * harmint(n1, m1, 2, diff1) * harmint(n2, m2, 1, diff2);
          tmp=tmp+ Qiiij(m1,m2) * harmint(n1, m1, 3, diff1)* harmint(n2, m2, 1, diff2);
          tmp=tmp+ Ciij(m2,m1) * harmint(n1, m1, 1, diff1) * harmint(n2, m2, 2, diff2);
          tmp=tmp+ Qiiij(m2,m1) * harmint(n1, m1, 1, diff1)* harmint(n2, m2, 3, diff2);
          IF (nMR > 2)THEN
            DO m=1,ndof
              IF (m == m1 .OR. m == m2)CYCLE;
              n = bra(m);
              tmp=tmp+ Qiijk(m,m1,m2) * harmint(n, m, 2, 0) * harmint(n1, m1, 1, diff1)&
              * harmint(n2, m2, 1, diff2);
            ENDDO
          ENDIF!nMR>2
        ENDIF
        ham(i+1,j+1)=tmp
      ENDDO!end for j
    ENDDO!end for i
    !  CALL writesquare(vciorder,ham*convert)
    if(vcidebug)PRINT*,"diagonilazation..."
    CALL diag(CIsize,ham,energies)
    !  CALL writesquare(vciorder,ham)
    PRINT*,'HOVCI(1)'
    !      PRINT*,energies*convert
    DO i=1,Nstate
      CALL printstateCI1(CIsize,ham(:,i),energies(i)*convert)
      IF(VCIprint)CALL printall(CIsize,ham(:,i),energies(i)*convert)
    ENDDO
    Eground=1.d10
    IF(debug)PRINT*,"getVCI1 :)"
  END SUBROUTINE getvci1

  SUBROUTINE printstate(size,eigenfunc,energy)
    USE global
    USE modStates
    IMPLICIT NONE
    INTEGER::state(ndof),maxpos,size,snum
    REAL*8::energy, eigenfunc(size),maxcoef
    maxpos=MAXLOC(ABS(eigenfunc),1)
    maxcoef=eigenfunc(maxpos)*eigenfunc(maxpos)
    snum=maxpos-1
    IF (restricted)then
      snum=snumvec(maxpos)
      CALL RCIstatelabel(snum,state)
    else
      CALL statelabel(snum,state)
    endif
    IF(SUM(state)==0)THEN
      IF (Eground>1.d8)then
        Eground=energy
        PRINT*,Eground
      ELSE
        Print*,'Weird state with energy',energy
      Endif
    ENDIF
    IF(SUM(state)==1)PRINT*,energy-Eground
    WRITE(17,170)energy-Eground,maxcoef,state
170 FORMAT(f20.12,2x,f5.2,x,15i3)
  END SUBROUTINE printstate

  SUBROUTINE printstateCI1(size,eigenfunc,energy)
    USE global
    USE modStates
    IMPLICIT NONE
    INTEGER::state(ndof),maxpos,size,snum
    REAL*8::energy, eigenfunc(size),maxcoef
    if(debug)print*,"printstate"
    maxpos=MAXLOC(ABS(eigenfunc),1)
    maxcoef=eigenfunc(maxpos)*eigenfunc(maxpos)
    snum=maxpos-1
    CALL CI1statelabel(snum,state)
    IF(SUM(state)==0)THEN
      IF (Eground>1.d8)then
        Eground=energy
        PRINT*,Eground
      ELSE
        Print*,'Weird state with energy',energy
      Endif
    ENDIF
    IF(SUM(state)==1)PRINT*,energy-Eground
    WRITE(17,170)energy-Eground,maxcoef,state
170 FORMAT(f20.12,2x,f5.2,x,15i3)
  END SUBROUTINE printstateCI1

  SUBROUTINE printall(size,eigenfunc,energy)
    USE global
    USE modStates
    IMPLICIT NONE
    INTEGER::state(ndof),i,size,snum
    REAL*8::energy, eigenfunc(size),weight(size)
    weight=eigenfunc*eigenfunc
    write(27,'(f20.12)')energy
    do i=1,size
      snum=i-1
      if(weight(i)>minweight)then
        IF (restricted .and. i>1) then
          snum=snumvec(i-1)
          CALL RCIstatelabel(snum,state)
        else
          CALL statelabel(snum,state)
        endif
        WRITE(27,270)weight(i)*100.d0,'%',state
      endif
    enddo !i
270 FORMAT(f4.1,a,x,35i3)
  END SUBROUTINE printall

  SUBROUTINE printallCI1(size,eigenfunc,energy)
    USE global
    USE modStates
    IMPLICIT NONE
    INTEGER::state(ndof),i,size,snum
    REAL*8::energy, eigenfunc(size),weight(size)
    weight=eigenfunc*eigenfunc
    write(27,'(f20.12)')energy
    do i=1,size
      snum=i-1
      if(weight(i)>minweight)then
        CALL CI1statelabel(snum,state)
        WRITE(27,270)weight(i)*100.d0,'%',state
      endif
    enddo !i
270 FORMAT(f4.1,a,x,35i3)
  END SUBROUTINE printallCI1

  SUBROUTINE vci_plain
    USE lapack
    USE yazar
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,j,m,m1,m2,m3,m4,n,n1,n2,n3,n4
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer,diff1,diff2,diff3,diff4
    REAL*8:: tmp
    REAL*8,ALLOCATABLE::ham(:,:),energies(:)
    IF(debug)PRINT*,"VCI starts"
    ALLOCATE(ham(FCIsize,FCIsize),energies(FCIsize))
    ham=0.d0
    if(VCIdebug) PRINT*,"FVCI Hamiltonian allocated and initialized, forming the H matrix of order ", FCIsize
    DO i=0,FCIsize-1
      tmp=0.d0
      CALL statelabel(i,ket)
      DO m=1,ndof
        n = ket(m);
        tmp=tmp+ Hii(m) * harmint(n, m, 2, 0) + Qiiii(m) * harmint(n, m, 4, 0);
        tmp=tmp+ (ket(m) + 0.5d0) * hrmfreq(m) * 0.5d0;
      ENDDO
      IF (nMR > 1)THEN
        DO m1=1,nDOF
          n1 = ket(m1);
          DO m2=m1+1, nDOF
            n2 = ket(m2);
            tmp=tmp+ Qiijj(m1,m2) * harmint(n1, m1, 2, 0) * harmint(n2, m2, 2, 0);
          ENDDO!m2
        ENDDO !m1
      ENDIF !nMR>1
      Ham(i+1,i+1)=tmp
      DO j=i+1,FCIsize-1
        tmp=0.d0
        CALL stateLabel(j,bra);
        CALL findiffer(bra,ket,diffvec,ndiffer)
        IF (nDiffer > 4) THEN
          CYCLE
        ELSE IF (nDiffer == 1)THEN
          diff1 = ABS(bra(diffVec(1)) - ket(diffVec(1)));
          n1 = MAX(bra(diffVec(1)), ket(diffVec(1)));
          tmp=tmp+ harmint(n1, diffVec(1), 1, diff1) * Gi(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 2, diff1) * Hii(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 3, diff1) * Ciii(diffVec(1));
          tmp=tmp+ harmint(n1, diffVec(1), 4, diff1) * Qiiii(diffVec(1));
          IF (diff1 == 2) tmp=tmp - DSQRT(DBLE(n1 * (n1 - 1))) * hrmfreq(diffVec(1)) * 0.25!kinetic term
          IF (nMR > 1)THEN
            DO m=1,ndof
              IF (m == diffVec(1))CYCLE
              n2 = MAX(bra(m), ket(m));
              diff2 = ABS(bra(m) - ket(m));
              tmp=tmp+ Qiijj(m,diffVec(1)) * harmint(n2, m, 2, diff2) * harmint(n1,diffVec(1), 2, diff1);
              tmp=tmp+ Ciij(m,diffVec(1)) * harmint(n2, m, 2, diff2) * harmint(n1,diffVec(1), 1, diff1);
            ENDDO !m1
          ENDIF!nmr>1
        ELSE IF ((nDiffer == 2) .AND. (nMR > 1)) THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          diff1 = ABS(bra(m1) - ket(m1));
          diff2 = ABS(bra(m2) - ket(m2));
          n1 = MAX(bra(m1), ket(m1));
          n2 = MAX(bra(m2), ket(m2));
          tmp=tmp+ Hij(m1,m2) * harmint(n1, m1, 1, diff1) * harmint(n2, m2, 1, diff2);
          tmp=tmp+ Qiijj(m1,m2) * harmint(n1, m1, 2, diff1)* harmint(n2, m2, 2, diff2);
          tmp=tmp+ Ciij(m1,m2) * harmint(n1, m1, 2, diff1) * harmint(n2, m2, 1, diff2);
          tmp=tmp+ Qiiij(m1,m2) * harmint(n1, m1, 3, diff1)* harmint(n2, m2, 1, diff2);
          tmp=tmp+ Ciij(m2,m1) * harmint(n1, m1, 1, diff1) * harmint(n2, m2, 2, diff2);
          tmp=tmp+ Qiiij(m2,m1) * harmint(n1, m1, 1, diff1)* harmint(n2, m2, 3, diff2);
          IF (nMR > 2)THEN
            DO m=1,ndof
              IF (m == m1 .OR. m == m2)CYCLE;
              n = bra(m);
              tmp=tmp+ Qiijk(m,m1,m2) * harmint(n, m, 2, 0) * harmint(n1, m1, 1, diff1)&
              * harmint(n2, m2, 1, diff2);
            ENDDO
          ENDIF

        ELSE IF ((nDiffer == 3) .AND. (nMR > 2))THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          m3 = diffVec(3);
          diff1 = ABS(bra(m1) - ket(m1));
          diff2 = ABS(bra(m2) - ket(m2));
          diff3 = ABS(bra(m3) - ket(m3));
          n1 = MAX(bra(m1), ket(m1));
          n2 = MAX(bra(m2), ket(m2));
          n3 = MAX(bra(m3), ket(m3));
          tmp=tmp+ Cijk(m1,m2,m3) * harmint(n1, m1, 1, diff1)* harmint(n2, m2, 1, diff2)&
          * harmint(n3, m3, 1, diff3);
          tmp=tmp+ Qiijk(m1,m2,m3) * harmint(n1, m1, 2, diff1) * harmint(n2, m2, 1,diff2)&
          * harmint(n3, m3, 1, diff3);
          tmp=tmp+ Qiijk(m2,m1,m3) * harmint(n2, m2, 2, diff2) * harmint(n1, m1, 1,diff1)&
          * harmint(n3, m3, 1, diff3);
          tmp=tmp+ Qiijk(m3,m2,m1) * harmint(n3, m3, 2, diff3) * harmint(n2, m2, 1,diff2)&
          * harmint(n1, m1, 1, diff1);

        ELSE IF ((nDiffer == 4) .AND. (nMR > 3))THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          m3 = diffVec(3);
          m4 = diffVec(4);
          diff1 = ABS(bra(m1) - ket(m1));
          diff2 = ABS(bra(m2) - ket(m2));
          diff3 = ABS(bra(m3) - ket(m3));
          diff4 = ABS(bra(m4) - ket(m4));
          n1 = MAX(bra(m1), ket(m1));
          n2 = MAX(bra(m2), ket(m2));
          n3 = MAX(bra(m3), ket(m3));
          n4 = MAX(bra(m4), ket(m4));
          tmp=tmp+ Qijkl(m1,m2,m3,m4) * harmint(n1, m1, 1, diff1)* harmint(n2, m2, 1, diff2)&
          * harmint(n3, m3, 1, diff3)* harmint(n4, m4, 1, diff4);
        ENDIF
        ham(i+1,j+1)=tmp
      ENDDO!end for j
    ENDDO!end for i
    CALL writesquare(FCIsize,ham*convert)
    if(vcidebug)PRINT*,"diagonilazation..."
    CALL diag(FCIsize,ham,energies)
    PRINT*,'VCI'
    !  PRINT*,energies(1)*convert
    DO m=1,Nstate
      PRINT*,energies(m)*convert!,(energies(m)-energies(1))*convert
    ENDDO
  END SUBROUTINE vci_plain

  SUBROUTINE getvscfci
    USE lapack
    USE yazar
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE modalintegral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::i,j,m1,m2,m3,m4,size,maxndiffer
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp
    REAL*8,ALLOCATABLE::ham(:,:),energies(:)
    IF(debug)PRINT*,"getvscfci"
    IF(VCIprint)write(27,'(a)')'VSCFCI'
    write(17,'(a)')'VSCFCI'
    ket=0
    Eground=1.d10
    CALL subvscf(ket)
    CALL form_modalint()
    size=FCIsize
    maxndiffer=nMR
    IF(restricted)THEN
      size=RCIsize
    ENDIF
    ALLOCATE(ham(size,size),energies(size))
    ham=0.d0
    if(vcidebug) PRINT*,"VSCF VCI Hamiltonian initialized with size ", size
    DO i=0,size-1
      tmp=0.d0
      IF(restricted)THEN
        CALL RCIstatelabel(snumvec(i+1),ket)
      ELSE
        CALL statelabel(i,ket)
      ENDIF
      Ham(i+1,i+1)=fullmodalint(ket,ket)
      DO j=i+1,size-1
        tmp=0.d0
        IF(restricted)THEN
          CALL RCIstatelabel(snumvec(j+1),bra)
        ELSE
          CALL statelabel(j,bra)
        ENDIF
        CALL findiffer(bra,ket,diffvec,ndiffer)

        IF (nDiffer > maxnDiffer) THEN
          CYCLE

        ELSE IF (nDiffer == 1)THEN
          m1=diffVec(1)
          tmp=tmp+ modalint1MR(m1,bra,ket)&
          +modalint(m1,bra(m1),ket(m1),0);          !kinetic term
          IF (nMR > 1)THEN
            DO m2=1,ndof
              IF (m2 == m1)CYCLE
              tmp=tmp+ modalint2MR(m1,m2,bra,ket)
              IF (nMR > 2)THEN
                DO m3=m2+1,ndof
                  IF (m3 == m1)CYCLE
                  tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket)
                  IF (nMR > 3)THEN
                    DO m4=m3+1,ndof
                      IF (m4 == m1 .OR. m4==m2)CYCLE
                      tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)
                    ENDDO!m4
                  ENDIF !nMR>3
                ENDDO!m3
              ENDIF !nMR>2
            ENDDO !m2
          ENDIF!nmr>1

        ELSE IF ((nDiffer == 2)) THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          tmp=tmp+ modalint2MR(m1,m2,bra,ket)
          IF (nMR > 2)THEN
            DO m3=1,ndof
              IF (m3 == m1 .OR. m3==m2)CYCLE
              tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket)
              IF (nMR > 3)THEN
                DO m4=m3+1,ndof
                  IF (m4 == m1 .OR. m4==m2)CYCLE
                  tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)
                ENDDO!m4
              ENDIF !nMR>3
            ENDDO!m3
          ENDIF !nMR>2

        ELSE IF ((nDiffer == 3))THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          m3 = diffVec(3);
          tmp=tmp+modalint3MR(m1,m2,m3,bra,ket)
          IF (nMR > 3)THEN
            DO m4=1,ndof
              IF (m4 == m1 .OR. m4==m2 .OR. m4==m3)CYCLE
              tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)
            ENDDO!m4
          ENDIF !nMR>3

        ELSE IF ((nDiffer == 4))THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          m3 = diffVec(3);
          m4 = diffVec(4);
          tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)
        ENDIF
        ham(i+1,j+1)=tmp
      ENDDO!end for j
    ENDDO!end for i
    !  CALL writesquare(vciorder,ham*convert)
    if(vcidebug)PRINT*,"diagonilazation..."
    CALL diag(size,ham,energies)
    !  CALL writesquare(vciorder,ham)
    PRINT*,'VSCFCI'
    DO i=1,Nstate
      CALL printstate(size,ham(:,i),energies(i)*convert)
      IF(VCIprint)CALL printall(size,ham(:,i),energies(i)*convert)
    ENDDO
    Eground=1.d10
    IF(debug)PRINT*,"getvscfci :)"
  END SUBROUTINE getvscfci

  SUBROUTINE getvscfCI1
    USE lapack
    USE yazar
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE modalintegral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::i,j,m1,m2,m3,m4,CIsize,maxndiffer
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp
    REAL*8,ALLOCATABLE::ham(:,:),energies(:)
    IF(debug)PRINT*,"getvscfci1"
    IF(VCIprint)write(27,'(a)')'VSCFCI1'
    write(17,'(a)')'VSCFCI(1)'
    ket=0
    Eground=1.d10
    CALL subvscf(ket)
    CALL form_modalint()
    CIsize=1+(maxexall-1)*ndof
    maxndiffer=min(nMR,2)
    ALLOCATE(ham(CIsize,CIsize),energies(CIsize))
    ham=0.d0
    if(vcidebug) PRINT*,"VSCF VCI(1) Hamiltonian initialized with size ", CIsize
    DO i=0,CIsize-1
      tmp=0.d0
      CALL CI1statelabel(i,ket)
      Ham(i+1,i+1)=fullmodalint(ket,ket)
      DO j=i+1,CIsize-1
        tmp=0.d0
        CALL CI1statelabel(j,bra)
        CALL findiffer(bra,ket,diffvec,ndiffer)

        IF (nDiffer > maxnDiffer) THEN
          CYCLE

        ELSE IF (nDiffer == 1)THEN
          m1=diffVec(1)
          tmp=tmp+ modalint1MR(m1,bra,ket)&
          +modalint(m1,bra(m1),ket(m1),0);          !kinetic term
          IF (nMR > 1)THEN
            DO m2=1,ndof
              IF (m2 == m1)CYCLE
              tmp=tmp+ modalint2MR(m1,m2,bra,ket)
              IF (nMR > 2)THEN
                DO m3=m2+1,ndof
                  IF (m3 == m1)CYCLE
                  tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket)
                  IF (nMR > 3)THEN
                    DO m4=m3+1,ndof
                      IF (m4 == m1 .OR. m4==m2)CYCLE
                      tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)
                    ENDDO!m4
                  ENDIF !nMR>3
                ENDDO!m3
              ENDIF !nMR>2
            ENDDO !m2
          ENDIF!nmr>1

        ELSE IF ((nDiffer == 2)) THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          tmp=tmp+ modalint2MR(m1,m2,bra,ket)
          IF (nMR > 2)THEN
            DO m3=1,ndof
              IF (m3 == m1 .OR. m3==m2)CYCLE
              tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket)
              IF (nMR > 3)THEN
                DO m4=m3+1,ndof
                  IF (m4 == m1 .OR. m4==m2)CYCLE
                  tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)
                ENDDO!m4
              ENDIF !nMR>3
            ENDDO!m3
          ENDIF !nMR>2
        ENDIF
        ham(i+1,j+1)=tmp
      ENDDO!end for j
    ENDDO!end for i
    !  CALL writesquare(vciorder,ham*convert)
    if(vcidebug)PRINT*,"diagonilazation..."
    CALL diag(CIsize,ham,energies)
    !  CALL writesquare(vciorder,ham)
    PRINT*,'VSCFCI(1)'
  !   PRINT*,energies(1)*convert,CIsize
    DO i=1,Nstate
      CALL printstateCI1(CIsize,ham(:,i),energies(i)*convert)
      IF(VCIprint)CALL printall(CIsize,ham(:,i),energies(i)*convert)
    ENDDO
    Eground=1.d10
    DEALLOCATE(ham,energies)
    IF(debug)PRINT*,"getvscfci1 :)"
  END SUBROUTINE getvscfCI1

  SUBROUTINE getvscfci2
    USE lapack
    USE yazar
    USE constants
    USE global
    USE modStates
    USE modalintegral,ONLY:form_modalint,fullmodalint
    USE modVSCF!,only:subvscf,fungetE
    IMPLICIT NONE
    INTEGER::i,j
    INTEGER::ket(ndof),bra(ndof),state(ndof)
    REAL*8,ALLOCATABLE::ham(:,:),energies(:)
    IF(debug)PRINT*,"getvscfci"
    state=0
    CALL subvscf(state)
    CALL form_modalint()
    ALLOCATE(ham(FCIsize,FCIsize),energies(FCIsize))
    ham=0.d0
    if(vcidebug)PRINT*,"FVCI Hamiltonian allocated and initialized, forming the H matrix of order ", FCIsize
    DO i=0,FCIsize-1
      CALL statelabel(i,ket)
      Ham(i+1,i+1)=fullmodalint(ket,ket)
      DO j=i+1,FCIsize-1
        CALL stateLabel(j,bra);
        ham(i+1,j+1)=fullmodalint(bra,ket)
      ENDDO!end for j
    ENDDO!end for i
    CALL writesquare(FCIsize,ham*convert)
    PRINT*,"diagonilazation..."
    CALL diag(FCIsize,ham,energies)
    PRINT*,'VCI'
    DO i=1,Nstate
      PRINT*,energies(i)*convert!,(energies(m)-energies(1))*convert
    ENDDO
  END SUBROUTINE getvscfci2

  SUBROUTINE vscfci2
    USE lapack
    USE yazar
    USE constants
    USE global
    USE modStates
    USE modalintegral
    USE modVSCF!,only:subvscf,fungetE
    IMPLICIT NONE
    INTEGER::i,j,m,m1,m2,m3,m4,reducedorder
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),state(ndof),ndiffer,diff1,diff2,diff3,diff4
    REAL*8:: tmp
    REAL*8,ALLOCATABLE::ham(:,:),energies(:)
    IF(debug)PRINT*,"VCI starts"
    state=0
    CALL subvscf(state,tmp)
    CALL form_modalint()
    ALLOCATE(ham(FCIsize,FCIsize),energies(FCIsize))
    ham=0.d0
    if(vcidebug)PRINT*,"FVCI Hamiltonian allocated and initialized, forming the H matrix of order ", FCIsize
    DO i=0,FCIsize-1
      tmp=0.d0
      CALL statelabel(i,ket)
      DO m=1,ndof
        tmp=tmp+ modalint1MR(m,ket,ket);
        tmp=tmp+ modalint(m,ket(m),ket(m),0);      !kinetic
      ENDDO
      IF (nMR > 1)THEN
        DO m1=1,nDOF
          DO m2=m1+1, nDOF
            tmp=tmp+ modalint2MR(m1,m2,ket,ket);
          ENDDO!m2
        ENDDO !m1
      ENDIF !nMR>1
      Ham(i+1,i+1)=tmp
      DO j=i+1,FCIsize-1
        tmp=0.d0
        CALL stateLabel(j,bra);
        CALL findiffer(bra,ket,diffvec,ndiffer)
        m1=diffVec(1)
        IF (nDiffer == 1)THEN
          tmp=tmp+ modalint(m1,bra(m1),ket(m1),0);          !kinetic
          tmp=tmp+ modalint1MR(m1,bra,ket);
          IF (nMR > 1)THEN
            DO m=1,ndof
              IF (m == m1)CYCLE
              tmp=tmp+ modalint2MR(m1,m,bra,ket);
            ENDDO !m1
          ENDIF!nmr>1
        ELSE IF ((nDiffer == 2) .AND. (nMR > 1)) THEN
          m2 = diffVec(2);
          tmp=tmp+ modalint2MR(m1,m2,bra,ket);
          IF (nMR > 2)THEN
            DO m=1,ndof
              IF (m == m1 .OR. m == m2)CYCLE;
              tmp=tmp+ modalint3MR(m1,m2,m,bra,ket);
            ENDDO
          ENDIF

        ELSE IF ((nDiffer == 3) .AND. (nMR > 2))THEN
          m2 = diffVec(2);
          m3 = diffVec(3);
          tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket);

        ELSE IF ((nDiffer == 4) .AND. (nMR > 3))THEN
          m1 = diffVec(1);
          m2 = diffVec(2);
          m3 = diffVec(3);
          m4 = diffVec(4);
          tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket);
        ENDIF
        ham(i+1,j+1)=tmp
      ENDDO!end for j
    ENDDO!end for i
    CALL writesquare(FCIsize,ham*convert)
    if(vcidebug)PRINT*,"diagonilazation..."
    CALL diag(FCIsize,ham,energies)
    PRINT*,'VSCFCI'
    !  PRINT*,energies(1)*convert
    DO m=1,Nstate
      PRINT*,energies(m)*convert!,(energies(m)-energies(1))*convert
    ENDDO
  END SUBROUTINE vscfci2

  SUBROUTINE formsnumvec2()
    USE global
    USE modstates
    IMPLICIT NONE
    INTEGER::i,j,state(ndof),limit
    IF(debug)Print*,"formsnumvec"
    j=0
    limit=min(maxRCIsize,FCIsize)-1
    DO i=0,limit
      CALL RCIstatelabel(i,state)
      IF (sum(state) > maxsum)cycle
      IF (count(state > 0) >  ncup)cycle
      j=j+1
      snumvec(j)=i
    ENDDO
    RCIsize=j
    IF(debug)Print*,"formsnumvec :)"
  END SUBROUTINE formsnumvec2
END MODULE modvci
