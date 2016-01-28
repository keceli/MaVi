MODULE modVMP2
  IMPLICIT NONE
CONTAINS

  SUBROUTINE subVPT1(state,EVPT1)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(ndof),b1,b2
    INTEGER::bra(ndof)
    REAL*8:: tmp,EVPT1,Evpt0
    bra=state
    EVPT0=0.d0
    EVPT1=0.d0
    IF(debug)   PRINT*,"VPT1 2",nmr
    tmp=0.d0
    DO m1=1,ndof
      Evpt0=Evpt0+hrmfreq(m1)*(dfloat(state(m1))+0.5d0)
      b1=bra(m1)
      tmp=tmp+ hrmint(m1,b1,b1,4) * Qiiii(m1)
      IF(nMR>1)THEN
        DO m2=m1+1,ndof
          b2=bra(m2)
          tmp=tmp+ hrmint(m1,b1,b1,2) * hrmint(m2,b2,b2,2) * Qiijj(m1,m2)
        ENDDO !m2
      ENDIF!nmr>1
    ENDDO !m1
    !print*,tmp,i
    EVPT1=tmp
    EVPT1=Evpt0+EVPT1
    IF(debug)   PRINT*,"Evpt0=",Evpt0*convert
    IF(debug)   PRINT*,"VPT1 :)"
  END SUBROUTINE subVPT1

  SUBROUTINE subVPT1a()
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::m1,m2
    REAL*8:: tmp,EVPT1,Evpt0,fundamental(ndof)
    EVPT0=0.d0
    EVPT1=0.d0
    fundamental=0.d0
    IF(debug)   PRINT*,"VPT1a"
    tmp=0.d0
    call del_coefs()
    DO m1=1,ndof
      Evpt0=Evpt0+hrmfreq(m1)*(0.5d0)
      DO m2=1,ndof
        if(m1==m2) then
          fundamental(m1)=fundamental(m1)+Qiiii(m1)/(8.d0*hrmfreq(m1)*hrmfreq(m2))
          tmp=tmp+ Qiiii(m1)/(32.d0*hrmfreq(m1)*hrmfreq(m2))
        else
          fundamental(m1)=fundamental(m1)+Qiijj(m1,m2)/(8.d0*hrmfreq(m1)*hrmfreq(m2))
          tmp=tmp+ Qiijj(m1,m2)/(32.d0*hrmfreq(m1)*hrmfreq(m2))
          print*,fundamental(m1)
        endif
      ENDDO !m2
    ENDDO !m1
    !print*,tmp,i
    EVPT1=tmp
    EVPT1=Evpt0+EVPT1
    call put_coefs()
    Print*,Evpt0*convert,evpt1*convert,(fundamental+hrmfreq)*convert
    IF(debug)   PRINT*,"Evpt0=",Evpt0*convert
    IF(debug)   PRINT*,"VPT1 :)"
  END SUBROUTINE subVPT1a

  SUBROUTINE subVPT1b()
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::m1,m2
    REAL*8:: tmp,EVPT1,Evpt0,fundamental(ndof)
    EVPT0=0.d0
    EVPT1=0.d0
    fundamental=0.d0
    IF(debug)   PRINT*,"VPT1b"
    tmp=0.d0
    DO m1=1,ndof
      Evpt0=Evpt0+hrmfreq(m1)*(0.5d0)
          tmp=tmp+0.75d0* Qiiii(m1)/(hrmfreq(m1)*hrmfreq(m1))
          if(nmr>1)then
      DO m2=m1+1,ndof
          tmp=tmp+ Qiijj(m1,m2)/(4.d0*hrmfreq(m1)*hrmfreq(m2))
      ENDDO !m2
      endif
    ENDDO !m1
       DO m1=1,ndof
          fundamental(m1)=fundamental(m1)+3.d0*Qiiii(m1)/(hrmfreq(m1)*hrmfreq(m1))
      if (nmr>1)then
      DO m2=1,ndof
      if(m1==m2)cycle
          fundamental(m1)=fundamental(m1)+Qiijj(m1,m2)/(2.d0*hrmfreq(m1)*hrmfreq(m2))
      ENDDO !m2
      endif
    ENDDO !m1
    !print*,tmp,i
    EVPT1=tmp
    EVPT1=Evpt0+EVPT1
    Print*,"Evpt1b"
    Print*,evpt1*convert
    do m1=1,ndof
    Print*,(fundamental(m1)+hrmfreq(m1))*convert
    enddo
    IF(debug)   PRINT*,"VPT1 :)"
  END SUBROUTINE subVPT1b

  SUBROUTINE subVPT2a(state,Evpt2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,m1,m2,m3,m4,state(ndof)
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,Evpt1,Evpt2,denom
    ! PRINT*,"vmp2 starts"

    bra=state
    CALL subVPT1(state,Evpt1)
    Evpt2=0.d0
    IF(benchmark) PRINT*,"vpt2 starts for state ",bra
    !Uzero=0.d0
    DO i=0,FCIsize-1
      CALL statelabel(i,ket)
      CALL findiffer(bra,ket,diffvec,ndiffer)
      !   PRINT*,"vmp2 loop",i,bra," ", ket," ", ndiffer
      IF ( (ndiffer>nMR).OR.(ndiffer==0)) CYCLE
      tmp=0.d0
      denom=0.d0
      SELECT CASE(ndiffer)
        CASE(1)
          m1=diffvec(1)
          tmp=tmp+ hrmint(m1,bra(m1),ket(m1),1) * (Gi(m1))&
          + hrmint(m1,bra(m1),ket(m1),3) * (Ciii(m1))&
          + hrmint(m1,bra(m1),ket(m1),4) * (Qiiii(m1));
          IF(nMR>1)THEN
            DO m2=m1+1,ndof
              tmp=tmp+hrmint2MR(m1,m2,bra,ket)
              IF(nMR>2)THEN
                DO m3=m2+1,ndof
                  IF((m3==m1))CYCLE
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
        CASE(2)
          m1=diffvec(1)
          m2=diffvec(2)
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
        CASE(3)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          tmp=tmp+ hrmint3MR(m1,m2,m3,bra,ket)
          IF(nMR>3)THEN
            DO m4=1,ndof
              IF((m4==m1).OR.(m4==m2) .OR. (m4==m3))CYCLE
              tmp=tmp+ hrmint4MR(m1,m2,m3,m4,bra,ket)
            ENDDO !m4
          ENDIF!nmr>3
        CASE(4)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          m4=diffvec(4)
          tmp=tmp+hrmint4MR(m1,m2,m3,m4,bra,ket)
      END SELECT
      ! denom=ndiffer*Evscf-Estate
      !    print*,tmp,i
      DO m1=1,ndof
        denom=denom+freq(m1)*(bra(m1)-ket(m1))
      ENDDO
      Evpt2=Evpt2+tmp*tmp/denom
    ENDDO!end for i
    !  PRINT*,"dEvmp2=",Evpt2,"for state", state
    Evpt2=Evpt1+Evpt2
  END SUBROUTINE subVPT2a

  SUBROUTINE subVPT2(state,Evpt2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::i1,i2,i3,i4,m1,m2,m3,m4,state(ndof)
    INTEGER::ket(ndof),bra(ndof)
    REAL*8:: Evpt2,Evpt1,tmp,denom
    IF(benchmark)  PRINT*,"vpt2 starts for state ",state
    bra=state
    CALL subVPT1(state,Evpt1)
    Evpt2=0.d0
    tmp=0.d0
    DO m1=1,ndof
      ket=bra
      DO i1=0,bra(m1)+4
        IF(i1==bra(m1))CYCLE
        ket(m1)=i1
        tmp=tmp+ hrmint(m1,bra(m1),ket(m1),1) * (Gi(m1))&
        + hrmint(m1,bra(m1),ket(m1),3) * (Ciii(m1))&
        + hrmint(m1,bra(m1),ket(m1),4) * (Qiiii(m1));
        IF(nMR>1)THEN
          DO m2=1,ndof
                IF((m2==m1))CYCLE
            tmp=tmp+hrmint2MR(m1,m2,bra,ket)
            IF(nMR>2)THEN
              DO m3=m2+1,ndof
                IF((m3==m1))CYCLE
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
        Evpt2=Evpt2+tmp*tmp/denom
        tmp=0.d0
      ENDDO !i1
    ENDDO !m1

    IF(nMR>1)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          ket=bra
          DO i1=0,bra(m1)+4
            IF(i1==bra(m1))CYCLE
            ket(m1)=i1
            DO i2=0,bra(m2)+4
              IF(i2==bra(m2))CYCLE
              ket(m2)=i2
              tmp=tmp+hrmint2MR(m1,m2,bra,ket)
              IF(nMR>2)THEN
                DO m3=1,ndof
                  IF((m3==m1).OR.(m3==m2))CYCLE
                  tmp=tmp+ hrmint3MR(m1,m2,m3,bra,ket)
                  IF(nMR>3)THEN
                    DO m4=m3+1,ndof
                      IF((m4==m1).OR.(m4==m2) )CYCLE
                      tmp=tmp+ hrmint4MR(m1,m2,m3,m4,bra,ket)
                    ENDDO !m4
                  ENDIF!nmr>3
                ENDDO !m3
              ENDIF!nmr>2
              denom=freq(m1)*(bra(m1)-ket(m1))&
              +freq(m2)*(bra(m2)-ket(m2))
              Evpt2=Evpt2+tmp*tmp/denom
              tmp=0.d0
            ENDDO !i2
          ENDDO !m2
        ENDDO!i1
      ENDDO !m1
    ENDIF!nmr>1

    IF(nMR>2)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            ket=bra
            DO i1=0,bra(m1)+4
              IF(i1==bra(m1))CYCLE
              ket(m1)=i1
              DO i2=0,bra(m2)+4
                IF(i2==bra(m2))CYCLE
                ket(m2)=i2
                DO i3=0,bra(m3)+4
                  IF(i3==bra(m3))CYCLE
                  ket(m3)=i3

                  tmp=tmp+ hrmint3MR(m1,m2,m3,bra,ket)
                  IF(nMR>3)THEN
                    DO m4=1,ndof
                      IF((m4==m1).OR.(m4==m2) .OR. (m4==m3))CYCLE
                      tmp=tmp+ hrmint4MR(m1,m2,m3,m4,bra,ket)
                    ENDDO !m4
                  ENDIF!nmr>3

                  denom=freq(m1)*(bra(m1)-ket(m1))&
                  +freq(m2)*(bra(m2)-ket(m2))&
                  +freq(m3)*(bra(m3)-ket(m3))
                  Evpt2=Evpt2+tmp*tmp/denom
                  tmp=0.d0
                ENDDO!i3
              ENDDO !i2
            ENDDO !i1
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>2

    IF(nMR>3)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            DO m4=m3+1,ndof
              ket=bra
              DO i1=0,bra(m1)+4
                IF(i1==bra(m1))CYCLE
                ket(m1)=i1
                DO i2=0,bra(m2)+4
                  IF(i2==bra(m2))CYCLE
                  ket(m2)=i2
                  DO i3=0,bra(m3)+4
                    IF(i3==bra(m3))CYCLE
                    ket(m3)=i3
                    DO i4=0,bra(m4)+4
                      IF(i4==bra(m4))CYCLE
                      ket(m4)=i4
                      tmp=tmp+hrmint4MR(m1,m2,m3,m4,bra,ket)

                      denom=freq(m1)*(bra(m1)-ket(m1))&
                      +freq(m2)*(bra(m2)-ket(m2))&
                      +freq(m3)*(bra(m3)-ket(m3))&
                      +freq(m4)*(bra(m4)-ket(m4))
                      Evpt2=Evpt2+tmp*tmp/denom
                      tmp=0.d0
                    ENDDO!i4
                  ENDDO!i3
                ENDDO !i2
              ENDDO !i1
            ENDDO!m4
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>2
    Evpt2=Evpt1+Evpt2
  END SUBROUTINE subVPT2

   SUBROUTINE subVPT20()
    USE constants
    USE global
    USE modQFF
    USE integral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::i1,i2,i3,i4,m1,m2,m3,m4
    REAL*8:: Evpt1,tmp,denom, E0
    REAL*8::Evpt2
    IF(debug)  PRINT*,"vpt20"
    Evpt2=0.d0
    tmp=0.d0
    DO m1=1,ndof
      DO i1=1,3,2
        tmp= hrmint(m1,0,i1,3) * Ciii(m1)
        denom=-freq(m1)*i1
        Evpt2=Evpt2+tmp*tmp/denom
        tmp=0.d0
      ENDDO !i1
    ENDDO !m1
    print*,'Evpt2= ', Evpt2*convert,(sum(freq)/2.d0+Evpt2)*convert
    DO m1=1,ndof
      DO m2=m1+1,ndof
        DO m3=m2+1,ndof
          print*,m1,m2,m3
          tmp=Qiijk(m1,m2,m3)*hrmint(m1,0,0,2)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)
          denom=-freq(m2)-freq(m3)
          Evpt2=Evpt2+tmp*tmp/denom
          print*,'Evpt2= ', Evpt2*convert,(sum(freq)/2.d0+Evpt2)*convert
          tmp=Qiijk(m2,m1,m3)*hrmint(m2,0,0,2)*hrmint(m1,0,1,1)*hrmint(m3,0,1,1)
          denom=-freq(m1)-freq(m3)
          Evpt2=Evpt2+tmp*tmp/denom
          print*,'Evpt2= ', Evpt2*convert,(sum(freq)/2.d0+Evpt2)*convert

          tmp=Qiijk(m3,m2,m1)*hrmint(m3,0,0,2)*hrmint(m2,0,1,1)*hrmint(m1,0,1,1)
          denom=-freq(m2)-freq(m1)
          Evpt2=Evpt2+tmp*tmp/denom
          print*,'Evpt2= ', Evpt2*convert,(sum(freq)/2.d0+Evpt2)*convert

          tmp=Qiijk(m1,m2,m3)*hrmint(m1,0,2,2)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)
          denom=-freq(m2)-freq(m3)-2.d0*freq(m1)
          Evpt2=Evpt2+tmp*tmp/denom
          tmp=Qiijk(m2,m1,m3)*hrmint(m2,0,2,2)*hrmint(m1,0,1,1)*hrmint(m3,0,1,1)
          denom=-freq(m1)-freq(m3)-2.d0*freq(m2)
          Evpt2=Evpt2+tmp*tmp/denom
          tmp=Qiijk(m3,m2,m1)*hrmint(m3,0,2,2)*hrmint(m2,0,1,1)*hrmint(m1,0,1,1)
          denom=-freq(m2)-freq(m1)-2.d0*freq(m3)
          Evpt2=Evpt2+tmp*tmp/denom
          print*,'Evpt2= ', Evpt2*convert,(sum(freq)/2.d0+Evpt2)*convert

        ENDDO!m3
      ENDDO!m2
    ENDDO !m1
    print*,'Evpt2= ', Evpt2*convert,(sum(freq)/2.d0+Evpt2)*convert
     ! WRITE(*,100)'VPT2 ZPE for Ciii',Evpt2
     ! 100 FORMAT(a,4)

  END SUBROUTINE subVPT20

  SUBROUTINE subXVPT1(state,EXVPT1)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,state(ndof),b1,b2
    INTEGER::bra(ndof)
    REAL*8:: tmp,EXVPT1,Exvpt0
    bra=state
    !Evscf=Evscf-Uzero(ndof)
    Exvpt0=0.d0
    EXVPT1=0.d0
    IF(debug)   PRINT*,"XVPT1"
    !Uzero=0.d0
    tmp=0.d0
    DO m1=1,ndof
      Exvpt0=Exvpt0+hrmfreq(m1)*(dfloat(state(m1))+0.5d0)
      b1=bra(m1)
      tmp=tmp+ hrmint(m1,b1,b1,4) * Qiiii(m1)
    !      IF(nMR>1)THEN
    !        DO m2=m1+1,ndof
    !          b2=bra(m2)
    !          tmp=tmp+ hrmint(m1,b1,b1,2) * hrmint(m2,b2,b2,2) * Qiijj(m1,m2)
    !        ENDDO !m2
    !      ENDIF!nmr>1
    ENDDO !m1
    ! denom=ndiffer*Evscf-Estate
    !print*,tmp,i
    EXVPT1=tmp
    EXVPT1=Exvpt0+EXVPT1
    IF(debug)   PRINT*,"XVPT1 :)"!"Exvpt0,dEXVPT1,EXVPT1 ",Exvpt0,tmp,EXVPT1,"for state ", state
  END SUBROUTINE subXVPT1

  SUBROUTINE subXVPT2(EXVPT2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,j,k,m
    REAL*8:: tmp,EXVPT2,denom
    IF(debug)   PRINT*,"XVPT2"
    EXVPT2=0.d0
    if (nMR>1)then
      do i=1,ndof
        denom=freq(i)
        tmp=0
        do m=1,ndof
          tmp=tmp+Ciij(m,i)*hrmint(i,0,1,3)
        enddo!m
        EXVPT2=EXVPT2-tmp*tmp/denom
      enddo!i

      if (nMR>2)then
        do i=1,ndof
          do j=i+1,ndof
            denom=freq(i)+freq(j)
            tmp=0
            do m=1,ndof
              tmp=tmp+Qiijk(m,i,j)*hrmint(i,0,1,4)*hrmint(j,0,1,4)
            enddo!m
            EXVPT2=EXVPT2-2*tmp*tmp/denom
            do k=j+1,ndof
              denom=freq(i)+freq(j)+freq(k)
              tmp=Cijk(i,j,k)*hrmint(i,0,1,3)*hrmint(j,0,1,3)*hrmint(k,0,1,3)
              EXVPT2=EXVPT2-6*tmp*tmp/denom
            enddo!k
          enddo!i
        enddo!j

        if (nMR>3)then
          do i=1,ndof
            do j=i+1,ndof
              do k=j+1,ndof
                do m=k+1,ndof
                  denom=freq(i)+freq(j)+freq(k)
                  tmp=Qijkl(i,j,k,m)*hrmint(i,0,1,4)*hrmint(j,0,1,4)*hrmint(k,0,1,4)*hrmint(m,0,1,4)
                  EXVPT2=EXVPT2-24*tmp*tmp/denom
                enddo!m
              enddo!k
            enddo!i
          enddo!j
        endif !nMR>3
      endif !nMR>2
    endif !nMR>1
    IF(debug)   PRINT*,"XVPT2 :)"
  END SUBROUTINE subXVPT2

  SUBROUTINE submattXVPT2(EXVPT2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4
    REAL*8:: tmp,EXVPT2,denom
    IF(debug)   PRINT*,"XVPT2"
    EXVPT2=0.d0
    call del_coefs
    if (nMR>1)then
      do m1=1,ndof
        denom=-freq(m1)
        do m2=1,ndof
          if(m2==m1)cycle
          do m3=1,ndof
            if((m3-m2)*(m3-m1)==0)cycle
            EXVPT2=EXVPT2+Ciij(m1,m2)*hrmint(m1,0,0,2)*hrmint(m2,0,1,1) &
            *Ciij(m3,m2)*hrmint(m3,0,0,2)*hrmint(m2,0,1,1)/denom
          enddo!m1
        enddo!m2
      enddo!m3

      if (nMR>2)then
        do m1=1,ndof
          do m2=1,ndof
            if(m2==m1)cycle
            denom=-freq(m1)-freq(m2)
            do m3=1,ndof
              if((m3-m2)*(m3-m1)==0)cycle
              do m4=1,ndof
                if((m4-m3)*(m4-m2)*(m4-m1)==0)cycle
                EXVPT2=EXVPT2+Qiijk(m3,m1,m2)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1)*hrmint(m3,0,0,2) &
                *Qiijk(m4,m1,m2)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1)*hrmint(m4,0,0,2)/denom
              enddo!m1
            enddo!m2
          enddo!m3
        enddo!m4

        do m1=1,ndof
          denom=-freq(m1)-freq(m2)-freq(m3)
          do m2=1,ndof
            if(m2==m1)cycle
            do m3=1,ndof
              if((m3-m2)*(m3-m1)==0)cycle
              EXVPT2=EXVPT2+(Cijk(m1,m2,m3)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1))**2.d0/denom
            enddo!m1
          enddo!m2
        enddo!m3

        if (nMR>3)then
          do m1=1,ndof
            do m2=1,ndof
              if(m2==m1)cycle
              denom=-freq(m1)-freq(m2)
              do m3=1,ndof
                if((m3-m2)*(m3-m1)==0)cycle
                do m4=1,ndof
                  if((m4-m3)*(m4-m2)*(m4-m1)==0)cycle
                  EXVPT2=EXVPT2+(Qijkl(m1,m2,m3,m4)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1) &
                  *hrmint(m3,0,1,1)*hrmint(m4,0,1,1))**2.d0/denom
                enddo!m1
              enddo!m2
            enddo!m3
          enddo!m4
        endif !nMR>3
      endif !nMR>2
    endif !nMR>1
    IF(debug)   PRINT*,"XVPT2 :)"
  END SUBROUTINE submattXVPT2

  SUBROUTINE subXVPT3(EXVPT3)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,m
    REAL*8:: int3,int3b1,int3b2,int4,int4b,EXVPT3,denom,numer
    IF(debug)   PRINT*,"XVPT2"
    EXVPT3=0.d0

    if (nMR>2)then
      do m1=1,ndof
        int3b1=0
        do m=1,ndof
          if(m==m1)cycle
          int3b1=int3b1+getint3b(m,m1)
        enddo!m
        do m2=m1+1,ndof
          int3b2=0
          do m=1,ndof
            if(m==m2)cycle
            int3b2=int3b2+getint3b(m,m2)
          enddo!m
          int4b=0
          do m=1,ndof
            if(m==m1.or. m==m2)cycle
            int4b=int4b+getint4b(m,m1,m2)
          enddo!m
          numer=int3b1*int3b2*int4b
          denom=freq(m1)*(freq(m1)+freq(m2))
          EXVPT3=EXVPT3+2.d0*numer/denom
          denom=freq(m1)*freq(m2)
          EXVPT3=EXVPT3+2.d0*numer/denom
          int3=0
          do m3=m2+1,ndof
            int3=int3+Cijk(m1,m2,m3)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)
            !    numer=int3*int3b*int4b
            denom=(freq(m1)+freq(m2)+freq(m3))*(freq(m1)+freq(m3))
            EXVPT3=EXVPT3+2.d0*numer/denom
            denom=freq(m1)*freq(m2)
            EXVPT3=EXVPT3+2.d0*numer/denom
            if (nMR>3)then
              int4=0
              do m4=m3+1,ndof
                int4=int4+Qijkl(m1,m2,m3,m4)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)*hrmint(m4,0,1,1)
              enddo!m4
            endif !nMR>3
          enddo!m3
        enddo!m2
      enddo!m1
    endif !nMR>2
    IF(debug)   PRINT*,"XVPT3 :)"
  END SUBROUTINE subXVPT3

  SUBROUTINE subfastXVPT3(EXVPT3)!not finished
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,m
    REAL*8:: int3,int3b,int4,int4b,EXVPT3,denom,numer
    IF(debug)   PRINT*,"XVPT2"
    EXVPT3=0.d0

    if (nMR>2)then
      do m1=1,ndof
        int3b=0
        do m=1,ndof
          if(m==m1)cycle
          int3b=int3b+Ciij(m,m1)*hrmint(m1,0,1,1)*hrmint(m,0,0,2)
        enddo!m
        do m2=m1+1,ndof
          int4b=0
          do m=1,ndof
            if(m==m1.or. m==m2)cycle
            int4b=int4b+Qiijk(m,m1,m2)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1)*hrmint(m,0,0,2)
          enddo!m
          numer=int3b*int3b*int4b
          denom=freq(m1)*(freq(m1)+freq(m2))
          EXVPT3=EXVPT3+2.d0*numer/denom
          denom=freq(m1)*freq(m2)
          EXVPT3=EXVPT3+2.d0*numer/denom
          int3=0
          do m3=m2+1,ndof
            int3=int3+Cijk(m1,m2,m3)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)
            numer=int3*int3b*int4b
            denom=(freq(m1)+freq(m2)+freq(m3))*(freq(m1)+freq(m3))
            EXVPT3=EXVPT3+2.d0*numer/denom
            denom=freq(m1)*freq(m2)
            EXVPT3=EXVPT3+2.d0*numer/denom
            if (nMR>3)then
              int4=0
              do m4=m3+1,ndof
                int4=int4+Qijkl(m1,m2,m3,m4)*hrmint(m1,0,1,1)*hrmint(m2,0,1,1)*hrmint(m3,0,1,1)*hrmint(m4,0,1,1)
              enddo!m4
            endif !nMR>3
          enddo!m3
        enddo!m2
      enddo!m1
    endif !nMR>2
    IF(debug)   PRINT*,"XVPT3 :)"
  END SUBROUTINE subfastXVPT3

  SUBROUTINE suboldXVPT2(state,EXVPT2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,m1,m2,m3,m4,state(ndof)
    INTEGER::b1,b2,b3,b4,k1,k2,k3,k4
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,EXVPT1,EXVPT2,denom
    ! PRINT*,"vmp2 starts"
    IF(debug)   PRINT*,"XVPT2"
    bra=state
    CALL subXVPT1(state,EXVPT1)
    EXVPT2=0.d0
    !Uzero=0.d0
    DO i=0,FCIsize-1
      CALL statelabel(i,ket)
      CALL findiffer(bra,ket,diffvec,ndiffer)
      !   PRINT*,"vmp2 loop",i,bra," ", ket," ", ndiffer
      IF ( (ndiffer>nMR).OR.(ndiffer==0)) CYCLE
      tmp=0.d0
      denom=0.d0
      SELECT CASE(ndiffer)
        CASE(1)
          m1=diffvec(1)
          tmp=tmp+ hrmint(m1,bra(m1),ket(m1),1) * (Gi(m1))
          IF(nMR>1)THEN
            DO m2=m1+1,ndof
              tmp=tmp+Ciij(m1,m2)*hrmint(m1,bra(m1),ket(m1),2)*hrmint(m2,bra(m2),ket(m2),1)&
              +Ciij(m2,m1)*hrmint(m1,bra(m1),ket(m1),1)*hrmint(m2,bra(m2),ket(m2),2)
            ENDDO !m2
          ENDIF!nmr>1
        CASE(2)
          m1=diffvec(1)
          m2=diffvec(2)
          tmp=tmp+Hij(m1,m2)*hrmint(m1,bra(m1),ket(m1),1)*hrmint(m2,bra(m2),ket(m2),1)
          IF(nMR>2)THEN
            DO m3=m2+1,ndof
              IF((m3==m1))CYCLE
              tmp=tmp+ modalint(m1,b1,k1,2) * modalint(m2,b2,k2,1) * modalint(m3,b3,k3,1) * Qiijk(m1,m2,m3)&
              +modalint(m1,b1,k1,1) * modalint(m2,b2,k2,2) * modalint(m3,b3,k3,1) * Qiijk(m2,m1,m3)&
              +modalint(m1,b1,k1,1) * modalint(m2,b2,k2,1) * modalint(m3,b3,k3,2) * Qiijk(m3,m1,m2)
            ENDDO !m3
          ENDIF!nmr>2
        CASE(3)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          tmp=tmp+ modalint(m1,b1,k1,1) * modalint(m2,b2,k2,1) * modalint(m3,b3,k3,1) * Cijk(m1,m2,m3)
        CASE(4)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          m4=diffvec(4)
          tmp=tmp+hrmint4MR(m1,m2,m3,m4,bra,ket)
      END SELECT
      ! denom=ndiffer*Evscf-Estate
      !    print*,tmp,i
      DO m1=1,ndof
        denom=denom+freq(m1)*(bra(m1)-ket(m1))
      ENDDO
      EXVPT2=EXVPT2+tmp*tmp/denom
    ENDDO!end for i
    !  PRINT*,"dEvmp2=",EXVPT2,"for state", state
    EXVPT2=EXVPT1+EXVPT2
    IF(debug)   PRINT*,"XVPT2 :)"
  END SUBROUTINE suboldXVPT2
  SUBROUTINE subVMP1(state,Evmp1)
    USE global
    USE constants
    USE modQFF
    USE modStates
    USE integral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof)
    REAL*8:: tmp,Evmp1,Evscf,Evmp0

    IF(benchmark) PRINT*,"vmp1 starts"

    gstate=0
    bra=state
    ket=bra
    Evmp1=0.d0
    tmp=0.d0
    Evmp0=0.d0
    IF(ssvmp)THEN
      CALL subvscf(state,Evscf)
    ELSE
      CALL subvscf(gstate)
      Evscf=fungetE(state)
    ENDIF

    DO m1=1,ndof
      ! call subformUZero(m1,Vscfcoefs(m1,:,:))
      !  print*,VSCFenergies(m1,b1+1),Uzero(m1)
      Evmp0=Evmp0+VSCFenergies(m1,bra(m1)+1)!+Uzero(m1)
      tmp=tmp+ modalint1MRU(m1,bra,ket)!modalint1MR(m1,bra,ket)-modalintU(m1,bra,ket)
      IF(nMR>1)THEN
        DO m2=m1+1,ndof
          tmp=tmp+modalint2MR(m1,m2,bra,ket)
          IF(nMR>2)THEN
            DO m3=m2+1,ndof
              tmp=tmp+modalint3MR(m1,m2,m3,bra,ket)
              IF(nMR>3)THEN
                DO m4=m3+1,ndof
                  tmp=tmp+modalint4MR(m1,m2,m3,m4,bra,ket)
                ENDDO !m4
              ENDIF!nmr>3
            ENDDO !m3
          ENDIF!nmr>2
        ENDDO !m2
      ENDIF!nmr>1
    ENDDO !m1
    ! denom=ndiffer*Evscf-Estate
    !print*,tmp,i
    Evmp1=tmp
    Evmp1=Evmp0+Evmp1
    PRINT*,"Evmp0=",Evmp0
    PRINT*,"Evmp1=",Evmp1
    IF(benchmark)  PRINT*,"Evmp0,Evmp1,0+1,Evscf ",Evmp0,tmp,Evmp1,Evscf," for state ", state
  END SUBROUTINE subVMP1

  SUBROUTINE sub1MRP1(state,E1MRP1)
    USE global
    USE constants
    USE modQFF
    USE modStates
    USE integral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof),nMRtmp
    REAL*8:: tmp,E1MRP1,E1MR,EVscf

    IF(debug) PRINT*,"sub1MRP1"

    gstate=0
    bra=state
    ket=bra
    E1MRP1=0.d0
    tmp=0.d0
    E1MR=0.d0
    nMRtmp=nMR
    nMR=1
 !   IF(ssvmp)THEN
 !     CALL subvscf(state,Evscf)
 !   ELSE
      CALL subvscf(gstate)
      Evscf=fungetE(state)
 !     print*, Evscf*convert, "vscf"
 !   ENDIF
      E1MR=Evscf
    nMR=nMRtmp
    DO m1=1,ndof
      IF(nMR>1)THEN
        DO m2=m1+1,ndof
          tmp=tmp+modalint2MR(m1,m2,bra,ket)
          IF(nMR>2)THEN
            DO m3=m2+1,ndof
              tmp=tmp+modalint3MR(m1,m2,m3,bra,ket)
              IF(nMR>3)THEN
                DO m4=m3+1,ndof
                  tmp=tmp+modalint4MR(m1,m2,m3,m4,bra,ket)
                ENDDO !m4
              ENDIF!nmr>3
            ENDDO !m3
          ENDIF!nmr>2
        ENDDO !m2
      ENDIF!nmr>1
    ENDDO !m1

    E1MRP1=tmp
    E1MRP1=E1MR+E1MRP1
!    PRINT*,"E1MR",E1MR
 !   PRINT*,"E1MRP1=",E1MRP1
        IF(debug) PRINT*,"sub1MRP1 :)"
  END SUBROUTINE sub1MRP1

  SUBROUTINE sub1MRP2(state,E1mrp2)
    USE constants
    USE global
    USE modStates
    USE modQFF
    USE modalintegral
    USE modVSCF!,only:subvscf,fungetE
    IMPLICIT NONE
    INTEGER::i1,i2,i3,i4,m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof),mthfund
    REAL*8:: E1MRP1,E1MRP2,Evscf
    IF(debug)  PRINT*,"sub1MRP2 for state",state
    gstate=0
    call sub1MRP1(state,E1MRP1)
    IF(SUM(state)==0)THEN
      CALL subvscf(state,Evscf)
      Evscfgs=Evscf
      CALL form_modalint()
    ELSE !for excited state !TODO just works for fundamentals
      mthfund=maxloc(state,1)
      Evscf=vscfenergies(mthfund,2)-vscfenergies(mthfund,1)+Evscfgs
    ENDIF
    bra=state
    E1MRP2=0.d0

    DO m1=1,ndof
      ket=bra
      DO i1=0,hrmbasis(m1)-1
        IF(i1==bra(m1))CYCLE
        ket(m1)=i1
        E1MRP2=E1MRP2+funEdiff1(m1,bra,ket)
      ENDDO !i1
    ENDDO !m1

    IF(nMR>1)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          ket=bra
          DO i1=0,hrmbasis(m1)-1
            IF(i1==bra(m1))CYCLE
            ket(m1)=i1
            DO i2=0,hrmbasis(m2)-1
              IF(i2==bra(m2))CYCLE
              ket(m2)=i2
              E1MRP2=E1MRP2+funEdiff2(m1,m2,bra,ket)
            ENDDO !i2
          ENDDO !m2
        ENDDO!i1
      ENDDO !m1
    ENDIF!nmr>1

    IF(nMR>2)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            ket=bra
            DO i1=0,hrmbasis(m1)-1
              IF(i1==bra(m1))CYCLE
              ket(m1)=i1
              DO i2=0,hrmbasis(m2)-1
                IF(i2==bra(m2))CYCLE
                ket(m2)=i2
                DO i3=0,hrmbasis(m3)-1
                  IF(i3==bra(m3))CYCLE
                  ket(m3)=i3
                  E1MRP2=E1MRP2+funEdiff3(m1,m2,m3,bra,ket)
                ENDDO!i3
              ENDDO !i2
            ENDDO !i1
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>2

    IF(nMR>3)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            DO m4=m3+1,ndof
              ket=bra
              DO i1=0,hrmbasis(m1)-1
                IF(i1==bra(m1))CYCLE
                ket(m1)=i1
                DO i2=0,hrmbasis(m2)-1
                  IF(i2==bra(m2))CYCLE
                  ket(m2)=i2
                  DO i3=0,hrmbasis(m3)-1
                    IF(i3==bra(m3))CYCLE
                    ket(m3)=i3
                    DO i4=0,hrmbasis(m4)-1
                      IF(i4==bra(m4))CYCLE
                      ket(m4)=i4
                      E1MRP2=E1MRP2+funEdiff4(m1,m2,m3,m4,bra,ket)
                    ENDDO!i4
                  ENDDO!i3
                ENDDO !i2
              ENDDO !i1
            ENDDO!m4
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>3
    IF(debug)  PRINT*,"dEvmp2=",E1MRP2,"for state", state
    E1MRP2=E1MRP1+E1MRP2
    IF(debug)  PRINT*,"sub1MRP2 :)"
    Return
  END SUBROUTINE sub1MRP2

  SUBROUTINE subVMP2(state,Evmp2)
    USE constants
    USE global
    USE modStates
    USE modQFF
    USE modalintegral
    USE modVSCF!,only:subvscf,fungetE
    IMPLICIT NONE
    INTEGER::i1,i2,i3,i4,m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof),mthfund
    REAL*8:: Evmp2,Evscf
    IF(debug)  PRINT*,"subVMP2 for state",state
    gstate=0
    IF(SUM(state)==0)THEN
      CALL subvscf(state,Evscf)
      Evscfgs=Evscf
      CALL form_modalint()
    ELSE IF(ssvmp) THEN !for excited state and state specific
      CALL subvscf(state,Evscf)
      CALL form_modalint()
    ELSE !for excited state !TODO just works for fundamentals
      mthfund=maxloc(state,1)
      ! Evscf=vscfenergies(mthfund,2)-vscfenergies(mthfund,1)
      !Evscf=fungetE(state)
      !  print*,(Evscf)*convert,state,"heyo",(fungetE(state)-Evscfgs)*convert
      Evscf=vscfenergies(mthfund,2)-vscfenergies(mthfund,1)+Evscfgs
    !   print*,(Evscf)*convert,"heyo",(fungetE(state))*convert
    ENDIF
    bra=state
    Evmp2=0.d0

    DO m1=1,ndof
      ket=bra
      DO i1=0,hrmbasis(m1)-1
        IF(i1==bra(m1))CYCLE
        ket(m1)=i1
        Evmp2=Evmp2+funEdiff1(m1,bra,ket)
      ENDDO !i1
    ENDDO !m1

    IF(nMR>1)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          ket=bra
          DO i1=0,hrmbasis(m1)-1
            IF(i1==bra(m1))CYCLE
            ket(m1)=i1
            DO i2=0,hrmbasis(m2)-1
              IF(i2==bra(m2))CYCLE
              ket(m2)=i2
              Evmp2=Evmp2+funEdiff2(m1,m2,bra,ket)
            ENDDO !i2
          ENDDO !m2
        ENDDO!i1
      ENDDO !m1
    ENDIF!nmr>1

    IF(nMR>2)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            ket=bra
            DO i1=0,hrmbasis(m1)-1
              IF(i1==bra(m1))CYCLE
              ket(m1)=i1
              DO i2=0,hrmbasis(m2)-1
                IF(i2==bra(m2))CYCLE
                ket(m2)=i2
                DO i3=0,hrmbasis(m3)-1
                  IF(i3==bra(m3))CYCLE
                  ket(m3)=i3
                  Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                ENDDO!i3
              ENDDO !i2
            ENDDO !i1
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>2

    IF(nMR>3)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            DO m4=m3+1,ndof
              ket=bra
              DO i1=0,hrmbasis(m1)-1
                IF(i1==bra(m1))CYCLE
                ket(m1)=i1
                DO i2=0,hrmbasis(m2)-1
                  IF(i2==bra(m2))CYCLE
                  ket(m2)=i2
                  DO i3=0,hrmbasis(m3)-1
                    IF(i3==bra(m3))CYCLE
                    ket(m3)=i3
                    DO i4=0,hrmbasis(m4)-1
                      IF(i4==bra(m4))CYCLE
                      ket(m4)=i4
                      Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                    ENDDO!i4
                  ENDDO!i3
                ENDDO !i2
              ENDDO !i1
            ENDDO!m4
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>3
    IF(debug)  PRINT*,"dEvmp2=",Evmp2,"for state", state
    Evmp2=Evscf+Evmp2
    IF(debug)  PRINT*,"subVMP2 :)"
    Return
  END SUBROUTINE subVMP2

  SUBROUTINE subVPT2test(state,Evpt2)
    USE constants
    USE global
    USE integral
    IMPLICIT NONE
    INTEGER::i1,i2,i3,i4,m1,m2,m3,m4,state(ndof)
    INTEGER::ket(ndof),bra(ndof)
    REAL*8:: Evpt2,Evpt1
    IF(debug)  PRINT*,"subVPT2 for state",state
     CALL subVPT1(state,Evpt1)
    bra=state
    Evpt2=0.d0

    DO m1=1,ndof
      ket=bra
      DO i1=0,hrmbasis(m1)-1
        IF(i1==bra(m1))CYCLE
        ket(m1)=i1
        Evpt2=Evpt2+funHOEdiff1(m1,bra,ket)
      ENDDO !i1
    ENDDO !m1

    IF(nMR>1)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          ket=bra
          DO i1=0,hrmbasis(m1)-1
            IF(i1==bra(m1))CYCLE
            ket(m1)=i1
            DO i2=0,hrmbasis(m2)-1
              IF(i2==bra(m2))CYCLE
              ket(m2)=i2
              Evpt2=Evpt2+funHOEdiff2(m1,m2,bra,ket)
            ENDDO !i2
          ENDDO !m2
        ENDDO!i1
      ENDDO !m1
    ENDIF!nmr>1

    IF(nMR>2)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            ket=bra
            DO i1=0,hrmbasis(m1)-1
              IF(i1==bra(m1))CYCLE
              ket(m1)=i1
              DO i2=0,hrmbasis(m2)-1
                IF(i2==bra(m2))CYCLE
                ket(m2)=i2
                DO i3=0,hrmbasis(m3)-1
                  IF(i3==bra(m3))CYCLE
                  ket(m3)=i3
                  Evpt2=Evpt2+funHOEdiff3(m1,m2,m3,bra,ket)
                ENDDO!i3
              ENDDO !i2
            ENDDO !i1
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>2

    IF(nMR>3)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            DO m4=m3+1,ndof
              ket=bra
              DO i1=0,hrmbasis(m1)-1
                IF(i1==bra(m1))CYCLE
                ket(m1)=i1
                DO i2=0,hrmbasis(m2)-1
                  IF(i2==bra(m2))CYCLE
                  ket(m2)=i2
                  DO i3=0,hrmbasis(m3)-1
                    IF(i3==bra(m3))CYCLE
                    ket(m3)=i3
                    DO i4=0,hrmbasis(m4)-1
                      IF(i4==bra(m4))CYCLE
                      ket(m4)=i4
                      Evpt2=Evpt2+funHOEdiff4(m1,m2,m3,m4,bra,ket)
                    ENDDO!i4
                  ENDDO!i3
                ENDDO !i2
              ENDDO !i1
            ENDDO!m4
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>3
    IF(debug)  PRINT*,"dEvpt2=",Evpt2,"for state", state
    Evpt2=Evpt1+Evpt2
    IF(debug)  PRINT*,"subVPT2test :)"
    Return
  END SUBROUTINE subVPT2test

  SUBROUTINE subVMP3(state,Evmp3)
    USE constants
    USE global
    USE modStates
    USE modQFF
    USE modalintegral
    USE modVSCF!,only:subvscf,fungetE
    IMPLICIT NONE
    INTEGER::i1,i2,i3,i4,m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof),mthfund
    REAL*8:: Evmp3,Evmp2,Evscf
    IF(debug)  PRINT*,"subVMP3 for state",state
    gstate=0
    Call subVMP2(state,Evmp2)
    IF(SUM(state)==0)THEN
      CALL subvscf(state,Evscf)
      Evscfgs=Evscf
      CALL form_modalint()
    ELSE IF(ssvmp) THEN !for excited state and state specific
      CALL subvscf(state,Evscf)
      CALL form_modalint()
    ELSE !for excited state !TODO just works for fundamentals
      mthfund=maxloc(state,1)
      ! Evscf=vscfenergies(mthfund,2)-vscfenergies(mthfund,1)
      !Evscf=fungetE(state)
      !  print*,(Evscf)*convert,state,"heyo",(fungetE(state)-Evscfgs)*convert
      Evscf=vscfenergies(mthfund,2)-vscfenergies(mthfund,1)+Evscfgs
    !   print*,(Evscf)*convert,"heyo",(fungetE(state))*convert
    ENDIF
    bra=state
    Evmp3=0.d0

    DO m1=1,ndof
      ket=bra
      DO i1=0,hrmbasis(m1)-1
        IF(i1==bra(m1))CYCLE
        ket(m1)=i1
        Evmp3=Evmp3+funEdiff1(m1,bra,ket)
      ENDDO !i1
    ENDDO !m1

    IF(nMR>1)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          ket=bra
          DO i1=0,hrmbasis(m1)-1
            IF(i1==bra(m1))CYCLE
            ket(m1)=i1
            DO i2=0,hrmbasis(m2)-1
              IF(i2==bra(m2))CYCLE
              ket(m2)=i2
              Evmp3=Evmp3+funEdiff2(m1,m2,bra,ket)
            ENDDO !i2
          ENDDO !m2
        ENDDO!i1
      ENDDO !m1
    ENDIF!nmr>1

    IF(nMR>2)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            ket=bra
            DO i1=0,hrmbasis(m1)-1
              IF(i1==bra(m1))CYCLE
              ket(m1)=i1
              DO i2=0,hrmbasis(m2)-1
                IF(i2==bra(m2))CYCLE
                ket(m2)=i2
                DO i3=0,hrmbasis(m3)-1
                  IF(i3==bra(m3))CYCLE
                  ket(m3)=i3
                  Evmp3=Evmp3+funEdiff3(m1,m2,m3,bra,ket)
                ENDDO!i3
              ENDDO !i2
            ENDDO !i1
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>2

    IF(nMR>3)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            DO m4=m3+1,ndof
              ket=bra
              DO i1=0,hrmbasis(m1)-1
                IF(i1==bra(m1))CYCLE
                ket(m1)=i1
                DO i2=0,hrmbasis(m2)-1
                  IF(i2==bra(m2))CYCLE
                  ket(m2)=i2
                  DO i3=0,hrmbasis(m3)-1
                    IF(i3==bra(m3))CYCLE
                    ket(m3)=i3
                    DO i4=0,hrmbasis(m4)-1
                      IF(i4==bra(m4))CYCLE
                      ket(m4)=i4
                      Evmp3=Evmp3+funEdiff4(m1,m2,m3,m4,bra,ket)
                    ENDDO!i4
                  ENDDO!i3
                ENDDO !i2
              ENDDO !i1
            ENDDO!m4
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>3
    IF(debug)  PRINT*,"dEvmp3=",Evmp3,"for state", state
    Evmp3=Evmp2+Evmp3
    IF(debug)  PRINT*,"subVMP3"
    RETURN
  END SUBROUTINE subVMP3

  SUBROUTINE subqVMP1(state,EqVMP1)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    USE modXVSCF
    IMPLICIT NONE
    INTEGER::m1,m2,m3,m4,state(ndof),gstate(ndof),b1
    INTEGER::bra(ndof),ket(ndof)
    REAL*8:: tmp,EqVMP1,EqvMP0,Eqvscf
    bra=state
    ket=bra
    gstate=0
    IF(ssvmp)THEN
      CALL subqvscf(state,Eqvscf)
      CALL subUtwo(state)
    ELSE
      CALL subqvscf(gstate)
      CALL subgetE(state,Eqvscf)
      CALL subUtwo(gstate)
    ENDIF
    !    CALL subUzero(state)

    EqVMP0=0.d0
    EqVMP1=0.d0
    IF(benchmark) PRINT*,"qVMP1 starts"
    !Uzero=0.d0
    tmp=0.d0
    DO m1=1,ndof
      Eqvmp0=EqVMP0+freq(m1)*(dfloat(state(m1))+0.5d0)
      b1=bra(m1)
      tmp=tmp+ hrmint(m1,bra(m1),ket(m1),1) * (Gi(m1))&
      + hrmint(m1,bra(m1),ket(m1),2) * (Hii(m1)-Utwo(m1))&
      + hrmint(m1,bra(m1),ket(m1),3) * (Ciii(m1))&
      + hrmint(m1,bra(m1),ket(m1),4) * (Qiiii(m1));
      IF(nMR>1)THEN
        DO m2=m1+1,ndof
          tmp=tmp+hrmint2MR(m1,m2,bra,ket)
          IF(nMR>2)THEN
            DO m3=m2+1,ndof
              IF((m3==m1))CYCLE
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
    ENDDO !m1
    ! denom=ndiffer*Evscf-Estate
    !print*,tmp,i
    EqVMP1=tmp
    EqVMP1=EqVMP0+EqVMP1
    IF(benchmark)  PRINT*,"EqVMP0,dEqVMP1,EqVMP1, Eqvscf ",EqVMP0,tmp,EqVMP1, EqVSCF," for state ", state
  END SUBROUTINE subqVMP1

  SUBROUTINE subqVMP2(state,Eqvmp2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    USE modXVSCF
    IMPLICIT NONE
    INTEGER::i,m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,Eqvmp1,Eqvmp2,Eqvscf,denom
    ! PRINT*,"vmp2 starts"

    gstate=0
    IF(ssvmp)THEN
      CALL subqvscf(state,Eqvscf)
      CALL subUzero(state)
      CALL subUtwo(state)
    ELSE
      CALL subqvscf(gstate)
      CALL subUzero(gstate)
      CALL subUtwo(gstate)
      CALL subgetE(state,Eqvscf)
    ENDIF

    CALL subqVMP1(state,Eqvmp1)
    bra=state
    Eqvmp2=0.d0
    IF(benchmark)   PRINT*,"qvmp2 starts for state ",bra
    !Uzero=0.d0
    DO i=0,FCIsize-1
      CALL statelabel(i,ket)
      CALL findiffer(bra,ket,diffvec,ndiffer)
      !   PRINT*,"vmp2 loop",i,bra," ", ket," ", ndiffer
      IF ( (ndiffer>nMR).OR.(ndiffer==0)) CYCLE
      tmp=0.d0
      denom=0.d0
      SELECT CASE(ndiffer)
        CASE(1)
          m1=diffvec(1)
          tmp=tmp+ hrmint(m1,bra(m1),ket(m1),1) * (Gi(m1))&
          + hrmint(m1,bra(m1),ket(m1),2) * (Hii(m1)-Utwo(m1))&
          + hrmint(m1,bra(m1),ket(m1),3) * (Ciii(m1))&
          + hrmint(m1,bra(m1),ket(m1),4) * (Qiiii(m1));
          IF(nMR>1)THEN
            DO m2=1,ndof
              IF(m2==m1)CYCLE
              tmp=tmp+hrmint2MR(m1,m2,bra,ket)
              IF(nMR>2)THEN
                DO m3=m2+1,ndof
                  IF((m3==m1))CYCLE
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
        CASE(2)
          m1=diffvec(1)
          m2=diffvec(2)
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
        CASE(3)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          tmp=tmp+ hrmint3MR(m1,m2,m3,bra,ket)
          IF(nMR>3)THEN
            DO m4=1,ndof
              IF((m4==m1).OR.(m4==m2) .OR. (m4==m3))CYCLE
              tmp=tmp+ hrmint4MR(m1,m2,m3,m4,bra,ket)
            ENDDO !m4
          ENDIF!nmr>3
        CASE(4)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          m4=diffvec(4)
          tmp=tmp+hrmint4MR(m1,m2,m3,m4,bra,ket)
      END SELECT
      ! denom=ndiffer*Evscf-Estate
      !    print*,tmp,i
      DO m1=1,ndof
        denom=denom+freq(m1)*(bra(m1)-ket(m1))
      ENDDO
      Eqvmp2=Eqvmp2+tmp*tmp/denom
    ENDDO!end for i
    !    PRINT*,"dEvmp2=",Eqvmp2,"for state", state
    Eqvmp2=Eqvmp1+Eqvmp2
  END SUBROUTINE subqVMP2

  SUBROUTINE subtestgreen()
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,j,m1,m2,m3,m4,state(ndof)
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,Evpt1,Evpt2,denom,deltafreq
    PRINT*,"Green -"
    DO m1=1,ndof
      deltafreq=0.d0
      DO m2=1,ndof
        DO m3=m2+1,ndof
          DO i=1,maxbasis!3,2
            DO j=1,maxbasis!3,2
              tmp=Cijk(m1,m2,m3)*hrmint(m1,0,1,1)*hrmint(m2,0,i,1)*hrmint(m2,0,j,1)
              deltafreq=deltafreq+tmp*tmp/(-hrmfreq(m1)-i*hrmfreq(m2)-j*hrmfreq(m3))
            ENDDO !j
          ENDDO !i
        ENDDO !m3
      ENDDO !m2
      PRINT*,(freq(m1)+deltafreq)*convert
    ENDDO !m1

  END SUBROUTINE subtestgreen
  SUBROUTINE subtestgreen4()
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,j,m1,m2,m3,m4,state(ndof)
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,Evpt1,Evpt2,denom,deltafreq
    PRINT*,"Green+"
    DO m1=1,ndof
      deltafreq=0.d0
      DO m2=1,ndof
        DO m3=m2+1,ndof
          DO i=1,maxbasis!3,2
            DO j=1,maxbasis!3,2
              tmp=Cijk(m1,m2,m3)*hrmint(m1,0,1,1)*hrmint(m2,0,i,1)*hrmint(m2,0,j,1)
              deltafreq=deltafreq+tmp*tmp/(hrmfreq(m1)-i*hrmfreq(m2)-j*hrmfreq(m3))
            ENDDO !j
          ENDDO !i
        ENDDO !m3
      ENDDO !m2
      PRINT*,(freq(m1)+deltafreq)*convert
    ENDDO !m1

  END SUBROUTINE subtestgreen4

  SUBROUTINE subtestgreen2()
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,j,k,m1,m2,m3,m4,state(ndof)
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,tmp2,Evpt1,Evpt2,denom,deltafreq
    PRINT*,"Green2 starts"
    !  i=1;j=1;k=1
    DO m1=1,ndof
      deltafreq=0.d0
      deltafreq=0.d0
      DO m2=1,ndof
        DO m3=m2+1,ndof
          DO i=1,maxbasis!3,2
            DO j=1,maxbasis!3,2
              tmp=Cijk(m1,m2,m3)*hrmint(m1,0,1,1)*hrmint(m2,0,i,1)*hrmint(m2,0,j,1)
              deltafreq=deltafreq+tmp*tmp/(hrmfreq(m1)-i*hrmfreq(m2)-j*hrmfreq(m3))
            ENDDO !j
          ENDDO !i
        ENDDO !m3
      ENDDO !m2
      if(nMR>3)then
        DO m2=1,ndof
          DO m3=m2+1,ndof
            DO m4=m3+1,ndof
              DO i=1,maxbasis!3,2
                DO j=1,maxbasis!3,2
                  DO k=1,maxbasis!3,2
                    tmp=Qijkl(m1,m2,m3,m4)*hrmint(m1,0,1,2)*hrmint(m2,0,i,1)*hrmint(m2,0,j,1)*hrmint(m2,0,k,1)
                    deltafreq=deltafreq+tmp*tmp/(hrmfreq(m1)-i*hrmfreq(m2)-j*hrmfreq(m3)-k*hrmfreq(m4))
                  ENDDO !k
                ENDDO !j
              ENDDO !i
            ENDDO !m4
          ENDDO !m3
        ENDDO !m2
      Endif
      PRINT*,(freq(m1)+deltafreq)*convert
    ENDDO !m1

  END SUBROUTINE subtestgreen2

  SUBROUTINE subtestgreen3()
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    IMPLICIT NONE
    INTEGER::i,j,k,m1,m2,m3,m4,state(ndof),iter
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,tmp2,Evpt1,Evpt2,denom,deltafreq,newfreq(ndof),oldfreq(ndof)
    PRINT*,"Green3 starts"
    !  i=1;j=1;k=1
    oldfreq=freq
    newfreq=freq
    iter=0
    Do
      iter=iter+1
      DO m1=1,ndof
        oldfreq(m1)=newfreq(m1)
        DO m2=1,ndof
          DO m3=m2+1,ndof
            DO i=1,maxbasis!3,2
              DO j=1,maxbasis!3,2
                tmp=Cijk(m1,m2,m3)*hrmint(m1,0,1,1)*hrmint(m2,0,i,1)*hrmint(m2,0,j,1)
                newfreq(m1)=newfreq(m1)+tmp*tmp/(hrmfreq(m1)-i*hrmfreq(m2)-j*hrmfreq(m3))
              ENDDO !j
            ENDDO !i
          ENDDO !m3
        ENDDO !m2
        if(nMR>3)then
          DO m2=1,ndof
            DO m3=m2+1,ndof
              DO m4=m3+1,ndof
                DO i=1,maxbasis!3,2
                  DO j=1,maxbasis!3,2
                    DO k=1,maxbasis!3,2
                      tmp=Qijkl(m1,m2,m3,m4)*hrmint(m1,0,1,2)*hrmint(m2,0,i,1)*hrmint(m2,0,j,1)*hrmint(m2,0,k,1)
                      newfreq(m1)=oldfreq(m1)+tmp*tmp/(hrmfreq(m1)-i*hrmfreq(m2)-j*hrmfreq(m3)-k*hrmfreq(m4))
                    ENDDO !k
                  ENDDO !j
                ENDDO !i
              ENDDO !m4
            ENDDO !m3
          ENDDO !m2
        endif
      ENDDO !m1
      IF (((SUM(ABS(newfreq-oldfreq))*convert) < scfthresh) .OR. (iter > maxiter)) EXIT
    Enddo
    PRINT*,(newfreq)*c_h2wn,iter
  END SUBROUTINE subtestgreen3
  !  SUBROUTINE projector
  !    USE global
  !    IMPLICIT NONE
  !
  !  END SUBROUTINE projector

  SUBROUTINE subVMP2a(state,Evmp2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::i,m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,Evmp2,Evscf,denom
    ! PRINT*,"vmp2 starts"
    gstate=0
    IF(ssvmp)THEN
      CALL subvscf(state,Evscf)
    ELSE
      CALL subvscf(gstate)
      Evscf=fungetE(state)
    ENDIF
    bra=state
    Evmp2=0.d0
    IF(benchmark)  PRINT*,"vmp2 starts for state ",bra
    !Uzero=0.d0
    DO i=0,FCIsize-1
      CALL statelabel(i,ket)
      CALL findiffer(bra,ket,diffvec,ndiffer)
      !   PRINT*,"vmp2 loop",i,bra," ", ket," ", ndiffer
      IF ( (ndiffer>nMR).OR.(ndiffer==0)) CYCLE
      tmp=0.d0
      denom=0.d0

      SELECT CASE(ndiffer)

        CASE(1)
          m1=diffvec(1)
          tmp=tmp+ modalint1MRU(m1,bra,ket)!modalint1MR-modalintU
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
        CASE(2)!nMR>1 guaranteed due to if clause above: ndiffer>nMR --> cycle
          m1=diffvec(1)
          m2=diffvec(2)
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
        CASE(3)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          tmp=tmp+ modalint3MR(m1,m2,m3,bra,ket)
          IF(nMR>3)THEN
            DO m4=1,ndof
              IF((m4==m1).OR.(m4==m2) .OR. (m4==m3))CYCLE
              tmp=tmp+ modalint4MR(m1,m2,m3,m4,bra,ket)! should be zero
            ENDDO !m4
          ENDIF!nmr>3
        CASE(4)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          m4=diffvec(4)
          tmp=tmp+modalint4MR(m1,m2,m3,m4,bra,ket)
      END SELECT
      ! denom=ndiffer*Evscf-Estate
      !    print*,tmp,i


      DO m1=1,ndof
        denom=denom+VSCFenergies(m1,bra(m1)+1)-VSCFenergies(m1,ket(m1)+1)
      ENDDO
      Evmp2=Evmp2+tmp*tmp/denom
    ENDDO!end for i
    IF(benchmark)  PRINT*,"dEvmp2=",Evmp2,"for state", state
    Evmp2=Evscf+Evmp2
  END SUBROUTINE subVMP2a

  SUBROUTINE subVMP2b(state,Evmp2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::i,m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof),diffvec(ndof),ndiffer
    REAL*8:: tmp,Evmp2,Evscf,denom
    ! PRINT*,"vmp2 starts"
    gstate=0
    IF(ssvmp)THEN
      CALL subvscf(state,Evscf)
    ELSE
      CALL subvscf(gstate)
      Evscf=fungetE(state)
    ENDIF
    bra=state
    Evmp2=0.d0
    IF(benchmark)  PRINT*,"vmp2 starts for state ",bra
    !Uzero=0.d0
    DO i=0,FCIsize-1
      CALL statelabel(i,ket)
      CALL findiffer(bra,ket,diffvec,ndiffer)
      !   PRINT*,"vmp2 loop",i,bra," ", ket," ", ndiffer
      IF ( (ndiffer>nMR).OR.(ndiffer==0)) CYCLE
      tmp=0.d0
      denom=0.d0

      SELECT CASE(ndiffer)

        CASE(1)
          m1=diffvec(1)
          tmp=funEdiff1(m1,bra,ket)
        CASE(2)!nMR>1 guaranteed due to if clause above: ndiffer>nMR --> cycle
          m1=diffvec(1)
          m2=diffvec(2)
          tmp=funEdiff2(m1,m2,bra,ket)
        CASE(3)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          tmp=funEdiff3(m1,m2,m3,bra,ket)
        CASE(4)
          m1=diffvec(1)
          m2=diffvec(2)
          m3=diffvec(3)
          m4=diffvec(4)
          tmp=funEdiff4(m1,m2,m3,m4,bra,ket)
      END SELECT
      Evmp2=Evmp2+tmp
    ENDDO!end for i
    IF(benchmark)  PRINT*,"dEvmp2=",Evmp2,"for state", state
    Evmp2=Evscf+Evmp2
  END SUBROUTINE subVMP2b


  SUBROUTINE subVMP2c(state,Evmp2)
    USE constants
    USE global
    USE modQFF
    USE modStates
    USE integral
    USE modVSCF
    IMPLICIT NONE
    INTEGER::i1,i2,i3,i4,m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof)
    REAL*8:: Evmp2,Evscf
    IF(benchmark)  PRINT*,"vmp2 starts for state ",state
    gstate=0
    IF(ssvmp)THEN
      CALL subvscf(state,Evscf)
    ELSE
      CALL subvscf(gstate)
      Evscf=fungetE(state)
    ENDIF
    bra=state
    Evmp2=0.d0

    DO m1=1,ndof
      ket=bra
      DO i1=1,4
        ket(m1)=bra(m1)+i1
        Evmp2=Evmp2+funEdiff1(m1,bra,ket)
        IF((bra(m1)-i1) > -1) THEN
          ket(m1)=bra(m1)-i1
          Evmp2=Evmp2+funEdiff1(m1,bra,ket)
        ENDIF
      ENDDO !i
    ENDDO !m1

    IF(nMR>1)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          ket=bra
          DO i1=1,hrmbasis(m1)-1
            ket(m1)=bra(m1)+i1
            DO i2=1,hrmbasis(m2)-1
              ket(m2)=bra(m2)+i2
              Evmp2=Evmp2+funEdiff2(m1,m2,bra,ket)
              IF(((bra(m1)-i1) > -1) .AND. ((bra(m2)-i2) > -1)) THEN
                ket(m1)=bra(m1)-i1
                ket(m2)=bra(m2)-i2
                Evmp2=Evmp2+funEdiff2(m1,m2,bra,ket)
              ENDIF
              IF((bra(m1)-i1) > -1) THEN
                ket(m1)=bra(m1)-i1
                ket(m2)=bra(m2)+i2
                Evmp2=Evmp2+funEdiff2(m1,m2,bra,ket)
              ENDIF
              IF((bra(m2)-i2) > -1) THEN
                ket(m1)=bra(m1)+i1
                ket(m2)=bra(m2)-i2
                Evmp2=Evmp2+funEdiff2(m1,m2,bra,ket)
              ENDIF
            ENDDO !i2
          ENDDO !m2
        ENDDO!i1
      ENDDO !m1
    ENDIF!nmr>1

    IF(nMR>2)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            ket=bra
            DO i1=1,4
              ket(m1)=bra(m1)+i1
              DO i2=1,4
                ket(m2)=bra(m2)+i2
                DO i3=1,4
                  ket(m3)=bra(m3)+i3
                  Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                  IF(((bra(m1)-i1) > -1) .AND. ((bra(m2)-i2) > -1).AND. ((bra(m3)-i3) > -1)) THEN
                    ket(m1)=bra(m1)-i1
                    ket(m2)=bra(m2)-i2
                    ket(m3)=bra(m3)-i3
                    Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                  ENDIF
                  IF(((bra(m1)-i1) > -1) .AND. ((bra(m2)-i2) > -1)) THEN
                    ket(m1)=bra(m1)-i1
                    ket(m2)=bra(m2)-i2
                    ket(m3)=bra(m3)+i3
                    Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                  ENDIF
                  IF(((bra(m1)-i1) > -1) .AND. ((bra(m3)-i3) > -1)) THEN
                    ket(m1)=bra(m1)-i1
                    ket(m2)=bra(m2)+i2
                    ket(m3)=bra(m3)-i3
                    Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                  ENDIF
                  IF( ((bra(m2)-i2) > -1).AND. ((bra(m3)-i3) > -1)) THEN
                    ket(m1)=bra(m1)+i1
                    ket(m2)=bra(m2)-i2
                    ket(m3)=bra(m3)-i3
                    Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                  ENDIF
                  IF(((bra(m1)-i1) > -1)) THEN
                    ket(m1)=bra(m1)-i1
                    ket(m2)=bra(m2)+i2
                    ket(m3)=bra(m3)+i3
                    Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                  ENDIF
                  IF(((bra(m2)-i2) > -1)) THEN
                    ket(m1)=bra(m1)+i1
                    ket(m2)=bra(m2)-i2
                    ket(m3)=bra(m3)+i3
                    Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                  ENDIF
                  IF( ((bra(m3)-i3) > -1)) THEN
                    ket(m1)=bra(m1)+i1
                    ket(m2)=bra(m2)+i2
                    ket(m3)=bra(m3)-i3
                    Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                  ENDIF
                ENDDO!i3
              ENDDO !i2
            ENDDO !i1
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>2

    IF(nMR>3)THEN
      DO m1=1,ndof
        DO m2=m1+1,ndof
          DO m3=m2+1,ndof
            DO m4=m3+1,ndof
              ket=bra
              DO i1=1,4
                ket(m1)=bra(m1)+i1
                DO i2=1,4
                  ket(m2)=bra(m2)+i2
                  DO i3=1,4
                    ket(m3)=bra(m3)+i3
                    DO i4=1,4
                      ket(m4)=bra(m4)+i4
                      Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      IF(((bra(m1)-i1) > -1) .AND. ((bra(m2)-i2) > -1).AND. &
                      ((bra(m3)-i3) > -1) .AND. ((bra(m4)-i4) > -1)) THEN
                        ket(m1)=bra(m1)-i1
                        ket(m2)=bra(m2)-i2
                        ket(m3)=bra(m3)-i3
                        ket(m4)=bra(m4)-i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m1)-i1) > -1) .AND. ((bra(m2)-i2) > -1).AND. ((bra(m3)-i3) > -1)) THEN
                        ket(m1)=bra(m1)-i1
                        ket(m2)=bra(m2)-i2
                        ket(m3)=bra(m3)-i3
                        ket(m4)=bra(m4)+i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m1)-i1) > -1) .AND. ((bra(m2)-i2) > -1).AND.  ((bra(m4)-i4) > -1)) THEN
                        ket(m1)=bra(m1)-i1
                        ket(m2)=bra(m2)-i2
                        ket(m3)=bra(m3)+i3
                        ket(m4)=bra(m4)-i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m1)-i1) > -1) .AND.  ((bra(m3)-i3) > -1).AND. ((bra(m4)-i4) > -1)) THEN
                        ket(m1)=bra(m1)-i1
                        ket(m2)=bra(m2)+i2
                        ket(m3)=bra(m3)-i3
                        ket(m4)=bra(m4)-i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m2)-i2) > -1).AND. ((bra(m3)-i3) > -1).AND. ((bra(m4)-i4) > -1)) THEN
                        ket(m1)=bra(m1)+i1
                        ket(m2)=bra(m2)-i2
                        ket(m3)=bra(m3)-i3
                        ket(m4)=bra(m4)-i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m1)-i1) > -1) .AND. ((bra(m2)-i2) > -1)) THEN
                        ket(m1)=bra(m1)-i1
                        ket(m2)=bra(m2)-i2
                        ket(m3)=bra(m3)+i3
                        ket(m4)=bra(m4)+i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m1)-i1) > -1) .AND. ((bra(m4)-i4) > -1)) THEN
                        ket(m1)=bra(m1)-i1
                        ket(m2)=bra(m2)+i2
                        ket(m3)=bra(m3)+i3
                        ket(m4)=bra(m4)-i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF( ((bra(m3)-i3) > -1).AND. ((bra(m4)-i4) > -1)) THEN
                        ket(m1)=bra(m1)+i1
                        ket(m2)=bra(m2)+i2
                        ket(m3)=bra(m3)-i3
                        ket(m4)=bra(m4)-i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m1)-i1) > -1) .AND. ((bra(m3)-i3) > -1)) THEN
                        ket(m1)=bra(m1)-i1
                        ket(m2)=bra(m2)+i2
                        ket(m3)=bra(m3)-i3
                        ket(m4)=bra(m4)+i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m2)-i2) > -1).AND. ((bra(m4)-i4) > -1)) THEN
                        ket(m1)=bra(m1)+i1
                        ket(m2)=bra(m2)-i2
                        ket(m3)=bra(m3)+i3
                        ket(m4)=bra(m4)-i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF( ((bra(m2)-i2) > -1).AND. ((bra(m3)-i3) > -1)) THEN
                        ket(m1)=bra(m1)+i1
                        ket(m2)=bra(m2)-i2
                        ket(m3)=bra(m3)-i3
                        ket(m4)=bra(m4)+i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m1)-i1) > -1) ) THEN
                        ket(m1)=bra(m1)-i1
                        ket(m2)=bra(m2)+i2
                        ket(m3)=bra(m3)+i3
                        ket(m4)=bra(m4)+i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF( ((bra(m2)-i2) > -1)) THEN
                        ket(m1)=bra(m1)+i1
                        ket(m2)=bra(m2)-i2
                        ket(m3)=bra(m3)+i3
                        ket(m4)=bra(m4)+i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF( ((bra(m3)-i3) > -1)) THEN
                        ket(m1)=bra(m1)+i1
                        ket(m2)=bra(m2)+i2
                        ket(m3)=bra(m3)-i3
                        ket(m4)=bra(m4)+i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                      IF(((bra(m4)-i4) > -1)) THEN
                        ket(m1)=bra(m1)+i1
                        ket(m2)=bra(m2)+i2
                        ket(m3)=bra(m3)+i3
                        ket(m4)=bra(m4)-i4
                        Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                      ENDIF
                    ENDDO!i4
                  ENDDO!i3
                ENDDO !i2
              ENDDO !i1
            ENDDO!m4
          ENDDO!m3
        ENDDO!m2
      ENDDO !m1
    ENDIF!nmr>2
    IF(benchmark)  PRINT*,"dEvmp2=",Evmp2,"for state", state
    Evmp2=Evscf+Evmp2
  END SUBROUTINE subVMP2c

  SUBROUTINE subVMP2f(state,Evmp2)
    USE constants
    USE global
    USE modStates
    USE modQFF
    USE modalintegral
    USE modVSCF!,only:subvscf,fungetE
    IMPLICIT NONE
    INTEGER::i1,i2,i3,i4,m1,m2,m3,m4,state(ndof),gstate(ndof)
    INTEGER::ket(ndof),bra(ndof)
    REAL*8:: Evmp2,Evscf
    IF(debug)  PRINT*,"subVMP2 for state",state
    gstate=0
    IF(SUM(state)==0)THEN
      CALL subvscf(state,Evscf)
      CALL form_modalint()
    ELSE IF(ssvmp) THEN !for excited state and state specific
      CALL subvscf(state,Evscf)
      CALL form_modalint()
    ELSE !for excited state
      Evscf=fungetE(state)
    ENDIF
    bra=state
    Evmp2=0.d0

    DO m1=1,ndof
      ket=bra
      DO i1=0,hrmbasis(m1)-1
        IF(i1==bra(m1))CYCLE
        ket(m1)=i1
        Evmp2=Evmp2+funEdiff1(m1,bra,ket)
        IF(nMR>1)THEN
          DO m2=m1+1,ndof
            DO i2=0,hrmbasis(m2)-1
              IF(i2==bra(m2))CYCLE
              ket(m2)=i2
              Evmp2=Evmp2+funEdiff2(m1,m2,bra,ket)
              IF(nMR>2)THEN
                DO m3=m2+1,ndof
                  DO i3=0,hrmbasis(m3)-1
                    IF(i3==bra(m3))CYCLE
                    ket(m3)=i3
                    Evmp2=Evmp2+funEdiff3(m1,m2,m3,bra,ket)
                    IF(nMR>3)THEN
                      DO m4=m3+1,ndof
                        DO i4=0,hrmbasis(m4)-1
                          IF(i4==bra(m4))CYCLE
                          ket(m4)=i4
                          Evmp2=Evmp2+funEdiff4(m1,m2,m3,m4,bra,ket)
                        ENDDO!i4
                        ket(m4)=bra(m4)
                      ENDDO !m4
                    ENDIF!nmr>3
                  ENDDO!i3
                  ket(m3)=bra(m3)
                ENDDO !m3
              ENDIF!nmr>2
            ENDDO !i2
            ket(m2)=bra(m2)
          ENDDO !m2
        ENDIF!nmr>1
      ENDDO !i1
      ket(m1)=bra(m1)
    ENDDO !m1

    IF(benchmark)  PRINT*,"dEvmp2=",Evmp2,"for state", state
    Evmp2=Evscf+Evmp2
  END SUBROUTINE subVMP2f

END MODULE modVMP2
