MODULE modformPES
  USE global
  IMPLICIT NONE
  INTEGER,PARAMETER:: maxline=100
  !REAL*8::Qshift(Ndof),deltaQ(Ndof),xcoord(Natom,3),qcoord(3*Natom,Ndof),eqbgeo(Nataom,3)
  REAL*8,ALLOCATABLE::Qshift(:),deltaQ(:),xcoord(:,:),qcoord(:,:),eqbgeo(:,:),EPES(:)
  CHARACTER::com(maxline)*120
  CHARACTER,ALLOCATABLE::atomlabel(:)*2 !natom
CONTAINS

  SUBROUTINE initFormPES
    IMPLICIT NONE
    INTEGER::i
    ALLOCATE(Qshift(ndof),xcoord(Natom,3),qcoord(3*Natom,Ndof),eqbgeo(Natom,3),xmass(natom),atomlabel(natom))
    CALL formDeltaQ()
  END SUBROUTINE initFormPES

  SUBROUTINE formDeltaQ
    USE constants
    !USE modformPES
    IMPLICIT NONE
    INTEGER::m
    REAL*8::constmp
    ALLOCATE(deltaQ(ndof))
    constmp=delx/(2.0D0*c_pi*DSQRT(c_c/c_h/c_NA*1.0D-23))
    ! deltaQ(mode)=delx/(2.0D0*c_pi*DSQRT(c_c*hrmfreq(mode)/c_h/c_NA*1.0D-23))
    DO m=1,ndof
       deltaQ(m)=constmp/DSQRT(hrmfreq(m))
    ENDDO
    PRINT*,"deltaQ= ",deltaQ
  END SUBROUTINE formDeltaQ

  SUBROUTINE oneMRQFFpoints(mode,step)
    !USE modformPES
    IMPLICIT NONE
    INTEGER::mode,step
    qshift=0.0d0
    qshift(mode)=DBLE(step)*deltaQ(mode)
  END SUBROUTINE oneMRQFFpoints

  SUBROUTINE twoMRQFFpoints(mode1,mode2,step1,step2)
    !USE modformPES
    !IMPLICIT NONE
    INTEGER::mode1,mode2,step1,step2
    qshift=0.0d0
    qshift(mode1)=DBLE(step1)*deltaQ(mode1)
    qshift(mode2)=DBLE(step2)*deltaQ(mode2)
  END SUBROUTINE twoMRQFFpoints

  SUBROUTINE threeMRQFFpoints(mode1,mode2,mode3,step1,step2,step3)
    !USE global
    !USE modformPES
    IMPLICIT NONE
    INTEGER::mode1,mode2,mode3,step1,step2,step3
    qshift=0.0d0
    qshift(mode1)=DBLE(step1)*deltaQ(mode1)
    qshift(mode2)=DBLE(step2)*deltaQ(mode2)
    qshift(mode3)=DBLE(step3)*deltaQ(mode3)
  END SUBROUTINE threeMRQFFpoints

  SUBROUTINE singlepoint
    !USE modformPES
    IMPLICIT NONE
    REAL*8 :: E
    CHARACTER :: fi*10,fo*10,charfilenumber*5,chartrack*5
    INTEGER::filenumber,track

    SELECT CASE(RunTyp)
    CASE('PES')
       OPEN(700,file='PESinfo',status='unknown')
       READ(155,*)E
       WRITE(700,'(f15.8,20f7.3)')E,qshift
    CASE('WRT')
       OPEN(701,file='WRTinfo',status='unknown')
       OPEN(70,file='track',status='old')
       READ(70,*)track
       E=track
       WRITE(chartrack,'(i5)')track
       WRITE(701,'(a,20f7.3)')chartrack,qcoord
       fi='inp.'//ADJUSTL(TRIM(chartrack))
       OPEN(7,file=ADJUSTL(TRIM(fi)),status='new')
       CALL Generalinp()
       CLOSE(7)
       WRITE(chartrack,'(i5)')track+1
       CALL System('echo '//ADJUSTL(TRIM(chartrack))//' > track')
       CLOSE(70)
    END SELECT
  END SUBROUTINE singlepoint

  SUBROUTINE Generalinp()
    !USE modformPES
    IMPLICIT NONE
    INTEGER :: i,j,k
    i=1
    !general input
    DO WHILE(Com(i)/='MaVi:geometry')
       WRITE(7,'(a)')TRIM(ADJUSTL(Com(i)))
       i=i+1
    ENDDO
    DO j=1,Natom
       WRITE(7,'(a2,3f15.8)') atomLabel(j),(xcoord(j,k),k=1,3)
    END DO
    DO j=i+1,maxline
       IF(Com(j)/='MaVi:end') THEN
          WRITE(7,'(a)')TRIM(ADJUSTL(Com(j)))
       ELSE
          EXIT!the loop
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE Generalinp

  SUBROUTINE q2x()
    !USE modformPES
    IMPLICIT NONE
    INTEGER :: i,j,m,k
    REAL*8  :: sqrtmass
    xcoord=eqbgeo
    k=0
    DO i=1,Natom
       sqrtmass = DSQRT(xmass(i))
       DO j=1,3
          k=k+1
          DO m=1,Ndof
             xcoord(i,j)=xcoord(i,j) + Qshift(m)*qcoord(k,m)/sqrtmass
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE q2x

  REAL*8 FUNCTION E1MR(mode,itmp)
    !USE modformPES,ONLY:EPES
    !USE global, ONLY: ndof
    IMPLICIT NONE
    INTEGER,INTENT(in)::mode,itmp
    E1MR=EPES((mode-1)*6+itmp+1)
  END FUNCTION E1MR

  REAL*8 FUNCTION E2MR(mode1,mode2,itmp)
    !USE modformPES,ONLY:EPES
    !USE global, ONLY: ndof
    IMPLICIT NONE
    INTEGER,INTENT(in)::mode1,mode2,itmp
    E2MR=EPES(((mode1-1)*(mode1-2)/2+mode2-1)*12+itmp+1+6*ndof)
  END FUNCTION E2MR

  REAL*8 FUNCTION E3MR(mode1,mode2,mode3,itmp)
    !USE modformPES,ONLY:EPES
    !USE global, ONLY: ndof
    IMPLICIT NONE
    INTEGER,INTENT(in)::mode1,mode2,mode3,itmp
    E3MR=EPES(((mode1-1)*(mode1-2)*(mode1-3)/6+(mode2-1)*(mode2-2)/2+mode3-1)*8+itmp+1+6*ndof+ndof*(ndof-1)*6)
  END FUNCTION E3MR

  SUBROUTINE form_1MRQFF
    !USE global,ONLY:Ndof
    !USE modformPES
    USE modQFF
    IMPLICIT NONE
    INTEGER::m
    DO m=1,ndof
       Gi(m)=(E1MR(m,3)-E1MR(m,4))/(2.d0*deltaQ(m))
       Hii(m)=(E1MR(m,3)-2.D0*Eref+E1MR(m,4))/(deltaQ(m)*deltaQ(m))
       Ciii(m)=(E1MR(m,1)-3.D0*E1MR(m,3)+3.D0*E1MR(m,4)-E1MR(m,6))/(deltaQ(m)*deltaQ(m)*deltaQ(m)*8.d0)
       Qiiii(m)=(E1MR(m,2)-4.D0*E1MR(m,3)+6.D0*Eref-4.D0*E1MR(m,4)+E1MR(m,5))/(deltaQ(m)*deltaQ(m)*deltaQ(m)*deltaQ(m))
    ENDDO
  END SUBROUTINE form_1MRQFF

  SUBROUTINE form_2MRQFF
    !USE global
    !USE modformPES
    USE modQFF
    IMPLICIT NONE
    INTEGER::m1,m2
    DO m1=2,ndof
       DO m2=1,m1-1
          Hij(m1,m2)=(E2MR(m1,m2,2)-E2MR(m1,m2,5)-E2MR(m1,m2,8)+E2MR(m1,m2,11))/(deltaQ(m1)*deltaQ(m2)*4.d0)
          Ciij(m1,m2)=(E2MR(m1,m2,2)-2.D0*E1MR(m2,3)+E2MR(m1,m2,8)-E2MR(m1,m2,5)+2.D0*E1MR(m2,4)-E2MR(m1,m2,11))/&
               (deltaQ(m1)*deltaQ(m1)*deltaQ(m2)*2.d0)
          Ciij(m2,m1)=(E2MR(m1,m2,2)-2.D0*E1MR(m1,3)+E2MR(m1,m2,5)-E2MR(m1,m2,8)+2.D0*E1MR(m1,4)-E2MR(m1,m2,11))/&
               (deltaQ(m1)*deltaQ(m2)*deltaQ(m2)*2.d0)
          Qiijj(m1,m2)=(E2MR(m1,m2,2)+E2MR(m1,m2,5)+E2MR(m1,m2,8)+E2MR(m1,m2,11)&
               -2.D0*(E1MR(m1,3)+E1MR(m1,4)+E1MR(m2,3)+E1MR(m2,4))&
               +4.D0*Eref)/(deltaQ(m1)*deltaQ(m1)*deltaQ(m2)*deltaQ(m2))
          Qiiij(m1,m2)=(E2MR(m1,m2,1)-3.D0*E2MR(m1,m2,2)+3.D0*E2MR(m1,m2,8)-E2MR(m1,m2,9)-E2MR(m1,m2,4)+3.D0*E2MR(m1,m2,5) &
               -3.D0*E2MR(m1,m2,11)+E2MR(m1,m2,12))/(deltaQ(m1)*deltaQ(m1)*deltaQ(m1)*deltaQ(m2)*16.d0)
          Qiiij(m2,m1)=(E2MR(m1,m2,3)-3.D0*E2MR(m1,m2,2)+3.D0*E2MR(m1,m2,5)-E2MR(m1,m2,6)-E2MR(m1,m2,7)+3.D0*E2MR(m1,m2,8) &
               -3.D0*E2MR(m1,m2,11)+E2MR(m1,m2,10))/(deltaQ(m1)*deltaQ(m2)*deltaQ(m2)*deltaQ(m2)*16.d0)
       ENDDO
    ENDDO
  END SUBROUTINE form_2MRQFF

  SUBROUTINE form_3MRQFF
    !USE global
    !USE modformPES
    USE modQFF
    IMPLICIT NONE
    INTEGER::m1,m2,m3
    DO m1=3,ndof
       DO m2=2,m1-1
          DO m3=1,m2-1
             Cijk(m1,m2,m3)=(E3MR(m1,m2,m3,1)-E3MR(m1,m2,m3,2)-E3MR(m1,m2,m3,3)+E3MR(m1,m2,m3,4)-E3MR(m1,m2,m3,5)&
                  +E3MR(m1,m2,m3,6)+E3MR(m1,m2,m3,7)-E3MR(m1,m2,m3,8))&
                  /(deltaQ(m1)*deltaQ(m2)*deltaQ(m3)*8.d0)
             Qiijk(m1,m2,m3)=(E3MR(m1,m2,m3,1)+E3MR(m1,m2,m3,2)-E3MR(m1,m2,m3,3)-E3MR(m1,m2,m3,4)&
                  -E3MR(m1,m2,m3,5)-E3MR(m1,m2,m3,6)+E3MR(m1,m2,m3,7)+E3MR(m1,m2,m3,8)-2.D0*(E2MR(m2,m3,2)&
                  -E2MR(m2,m3,5)-E2MR(m2,m3,8)+E2MR(m2,m3,11)))&
                  /(deltaQ(m1)*deltaQ(m1)*deltaQ(m2)*deltaQ(m3)*4.d0)
             Qiijk(m2,m3,m1)=(E3MR(m1,m2,m3,1)-E3MR(m1,m2,m3,2)+E3MR(m1,m2,m3,3)-E3MR(m1,m2,m3,4)-E3MR(m1,m2,m3,5)&
                  +E3MR(m1,m2,m3,6)-E3MR(m1,m2,m3,7)+E3MR(m1,m2,m3,8)-2.D0*(E2MR(m1,m3,2)-E2MR(m1,m3,5)&
                  -E2MR(m1,m3,8)+E2MR(m1,m3,11)))&
                  /(deltaQ(m1)*deltaQ(m2)*deltaQ(m2)*deltaQ(m3)*4.d0)
             Qiijk(m3,m1,m2)=(E3MR(m1,m2,m3,1)-E3MR(m1,m2,m3,2)-E3MR(m1,m2,m3,3)+E3MR(m1,m2,m3,4)+E3MR(m1,m2,m3,5)&
                  -E3MR(m1,m2,m3,6)-E3MR(m1,m2,m3,7)+E3MR(m1,m2,m3,8)-2.D0*(E2MR(m1,m2,2)-E2MR(m1,m2,5)&
                  -E2MR(m1,m2,8)+E2MR(m1,m2,11)))&
                  /(deltaQ(m1)*deltaQ(m2)*deltaQ(m3)*deltaQ(m3)*4.d0)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE form_3MRQFF

  SUBROUTINE form_QFF
    !USE global
    !USE modformPES
    USE modQFF
    IMPLICIT NONE
    CALL init_QFF()
    CALL readPES()
    CALL formDeltaQ()
    CALL form_1MRQFF()
    IF (nMR>1)THEN
       CALL form_2MRQFF()
       IF (nMR>2)THEN
          CALL form_3MRQFF()
       ENDIF
    ENDIF
  END SUBROUTINE form_QFF

  SUBROUTINE formPES
    !!USE modformPES
    IMPLICIT NONE
    INTEGER:: i,m,m1,m2,m3,m4
    INTEGER::loop1MR(6),loop2MRa(12),loop2MRb(12),loop3MRa(8),loop3MRb(8),loop3MRc(8)
    DATA loop1MR/3,2,1,-1,-2,-3/
    DATA loop2MRa/3,1,1,3,1,1,-1,-1,-3,-1,-1,-3/
    DATA loop2MRb/1,1,3,-1,-1,-3,3,1,1,-3,-1,-1/
    DATA loop3MRa/1,-1,1,-1,1,-1,1,-1/
    DATA loop3MRb/1,1,-1,-1,1,1,-1,-1/
    DATA loop3MRc/1,1,1,1,-1,-1,-1,-1/
    qshift=0.d0
    CALL q2x()
    CALL singlepoint()
    DO m=1,ndof
       DO i=1,6
          CALL oneMRQFFpoints(m,loop1MR(i))
          CALL q2x()
          CALL singlepoint()
       ENDDO
    ENDDO
    IF (nmr>1) THEN

       DO m1=1,ndof
          DO m2=m1+1,ndof
             DO i=1,12
                CALL twoMRQFFpoints(m1,m2,loop2MRa(i),loop2MRb(i))
                CALL q2x()
                CALL singlepoint()
             ENDDO
          ENDDO
       ENDDO
       IF (nmr>2) THEN

          DO m1=1,ndof
             DO m2=m1+1,ndof
                DO m3=m2+1,ndof
                   DO i=1,8
                      CALL threeMRQFFpoints(m1,m2,m3,loop3MRa(i),loop3MRb(i),loop3MRc(i))
                      CALL q2x()
                      CALL singlepoint()
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF !nMr>1
  END SUBROUTINE formPES

  SUBROUTINE readPES
    IMPLICIT NONE
    INTEGER::i
    OPEN(UNIT=77,FILE='PES',STATUS='OLD')
    DO i=1,Ncalc
       READ(77,100)Epes(i)
    ENDDO !i
100 FORMAT(f30.20)
    Eref=Epes(1)
    CLOSE(77)
    PRINT*,"PES read successfully..."
  END SUBROUTINE readPES

END MODULE modformPES


