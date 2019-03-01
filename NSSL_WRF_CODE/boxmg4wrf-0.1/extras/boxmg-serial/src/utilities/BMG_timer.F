      SUBROUTINE BMG_timer(t)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This is a generic timer subroutine that provides a few different
C     mechanisms to return an estimate of the current system time in 
C     seconds. CPPFLAGS control which code is actually compiled.  The
C     options include, 
C
C     ETIME      - as an external function
C     ETIME      - as an intrinsic subroutine
C     CPU_TIME   - the f95 standard timer subroutine
C
C     For g77 compilers ETIME is included in libg2c, and for many vendor 
C     compilers it is in libU77.
C 
C =======================================================================
C $license_flag$
C =======================================================================
C ==========================================================================

#if defined (BOXMG_TIMER_ETIMEe )

      IMPLICIT NONE

      REAL*8  t
      REAL*4  temp(2)
 
      REAL*4   ETIME
      EXTERNAL ETIME

C     The call is of the form
C 
C         t = ETIME(temp)
C
C     upon return temp(1) = "user time" and temp(2)= "system time"
C     The function value is the sum of the two.
C
      t = ETIME(temp)

#elif defined (BOXMG_TIMER_ETIMEi)

      IMPLICIT NONE

      REAL*8  t
      REAL*4  temp(2), result
 
C     The call is of the form
C 
C     CALL ETIME(temp,result)
C
C     upon return temp(1) = "user time" and temp(2)= "system time"
C     and result has the sum of the two.
C
      CALL ETIME(temp,result)
      t=result


#elif defined (BOXMG_TIMER_F95_CPU_TIME)

      IMPLICIT NONE

      REAL(8)  t
 
C     The call is of the form
C 
C     CALL CPU_TIME(temp)
C
C     upon return temp holds the CPU time in seconds.
C
      CALL CPU_TIME(t)

#endif

C ==========================================================================

      RETURN
      END

