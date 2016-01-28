MODULE modStates
  IMPLICIT NONE
CONTAINS
  SUBROUTINE statelabel(nthstate,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m,state,nthstate,ket(ndof),n
    state = nthState;
    DO m=1,ndof
       n=hrmbasis(m)
       ket(m) = MOD(state,n);
       state = state / n;
    ENDDO
    RETURN
  END SUBROUTINE statelabel



    SUBROUTINE statelabelold(nthstate,ket)
    USE global
    IMPLICIT NONE
    INTEGER::m,state,nthstate,ket(ndof)
    state = nthState;
    DO m=1,ndof
       ket(m) = MOD(state,hrmbasis(m));
       state = state / hrmbasis(m);
    ENDDO
    RETURN
  END SUBROUTINE statelabelold

  SUBROUTINE findiffer(bra,ket,differvec,ndiffer)
    USE global
    IMPLICIT NONE
    INTEGER::i,m, bra(ndof),ket(ndof),differvec(ndof),ndiffer
    i = 0;
    differvec=0
    DO m=1,ndof
       IF (bra(m) /= ket(m))THEN
          differvec(i+1) = m;
          i=i+1;
       ENDIF
    ENDDO
    ndiffer=i
  END SUBROUTINE findiffer

  LOGICAL FUNCTION samestate(state,thestate)
    USE global
    IMPLICIT NONE
    INTEGER,INTENT(in)::state(ndof),thestate(ndof)
    INTEGER::m
    samestate=.TRUE.
    DO m=1,ndof
       IF(state(m).NE. thestate(m)) samestate=.FALSE.
    ENDDO
  END FUNCTION samestate


  LOGICAL FUNCTION skipstate(mode,ndiffer,diffvec)
    USE global
    IMPLICIT NONE
    INTEGER,INTENT(in)::mode,ndiffer,diffvec(ndiffer)
    INTEGER::m
    skipstate=.FALSE.
    DO m=1,ndiffer
       IF(mode== diffvec(m)) THEN
          skipstate=.TRUE.
          RETURN
       ENDIF
    ENDDO
  END FUNCTION skipstate


subroutine countstates()
USE global
implicit none
integer::i,state(ndof)
do i=0,FCIsize-1
call statelabel(i,state)
print*,i, "th state ",state
enddo
!print*,maxcoef,"for state",state
end subroutine countstates

subroutine countnonzeros(state,nexcited)
USE global
implicit none
integer,intent(in)::state(ndof)
integer,intent(out)::nexcited
integer::m
nexcited=0
do m=1,ndof
   if(state(m)>0)nexcited=nexcited+1
enddo
!print*,maxcoef,"for state",state
end subroutine countnonzeros


END MODULE modStates

