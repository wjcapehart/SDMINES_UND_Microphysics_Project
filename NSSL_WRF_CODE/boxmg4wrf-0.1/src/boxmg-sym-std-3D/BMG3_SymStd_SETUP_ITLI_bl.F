      SUBROUTINE BMG3_SymStd_SETUP_ITLI_bl(
     &                KGF, KGC, SO, SOC, CI, IIF, JJF, KKF, IIC, JJC,
     &                KKC, iGs, jGs, kGs, NOG, IFD, NStncl, NOGm, IGRD,
     &                iWORK, NMSGi, pSI_MSG, pMSGSO,
     &                BUFFER, NMSGr, NProc, ProcGrid, MyProcI, MyProcJ,
     &                MyProcK, NProcI, NProcJ, NProcK, DimX, DimY, DimZ,
     &                MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Construct the variational coarse-grid operator for either a
C     7-point or 27-point fine-grid discretization using the "block
C     local" technique.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C ----------------------------
C     Includes
C
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C
      INTEGER iGs, IIC, IIF, IFD, jGs, JJC, JJF, KGC, KGF, kGs, KKC, 
     &        KKF, NMSGi, NMSGr, NOG, NOGm, NProc, NProcI, NProcJ,
     &        NProcK, NStncl

      INTEGER iWork(NMSGi), pMSGSO(NBMG_pMSG,NOGm)
      INTEGER IGRD(NOGm,29), ProcGrid(NProcI,NProcJ, NProcK),
     &        MyProcI, MyProcJ, MyProcK, MPICOMM
      REAL*RKIND  CI(IIC,JJC,KKC,26), SO(IIF+1,JJF+1,KKF+1,NStncl), 
     &        SOC(IIC+1,JJC+1,KKC+1,14),
     &        BUFFER(NMSGr)

      Integer DimX(NprocI,NOGm), DimY(NprocJ,NOGm), DimZ(Nprock,NOGm)


C ----------------------------
C     Externals
C
      REAL*RKIND   BMG3_SymStd_UTILS_reswt
      EXTERNAL BMG3_SymStd_UTILS_reswt

C ----------------------------
C     Local Declarations
C
      INTEGER I, IC, IIC1, IIF1, IL, ILG, ILGP, J, JC, JJC1, JJF1,
     &        JL, JLG, JLGP, K, KC, KKC1, KKF1, KL, KLG, KLGP, KPZ,
     &        pSI_MSG, pMSGzo(NBMG_pMSG,2), ptrn, ierror
      REAL*RKIND  Q(7,7,7), R(7,7,7)
      INTEGER LXGP, RXGP, LYGP, RYGP, LZGP, RZGP
      
      INTEGER MyProc, ISTART, JSTART, KSTART, MPI_IERROR

C ==========================================================================

      CALL BMG3_SymStd_SETUP_PtrMSG(
     &          IGRD(KGF,idL_BMG_NLx), IGRD(KGF,idL_BMG_NLy), 
     &          IGRD(KGC,idL_BMG_NLz), iZERO, iZERO, iZERO, 
     &          pSI_MSG, NProc, 2, 1, 1, pMSGzo
     &          )

      MyProc = MSG_MyProc(MPICOMM)

      IIC1 = IIC-1
      JJC1 = JJC-1
      KKC1 = KKC-1
      IIF1 = IIF-1
      JJF1 = JJF-1
      KKF1 = KKF-1

      IF (MOD(iGs,2).EQ.1) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF
      
      IF (MOD(jGs,2).EQ.1) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF

      IF(MOD(kGs,2).EQ.1) THEN
         KSTART = 0
      ELSE
         KSTART = 1
      ENDIF

      
ce      WRITE(*,*) '(', MyProc, ')IIF  = ', IIF
ce
ce      IF (MyProc.EQ.1) THEN
ce         WRITE(*,*) IIF, JJF, KKF
ce      ENDIF
ce
ce      CALL MPI_FINALIZE(MPI_IERROR)
c      STOP



      IF (IFD.NE.1 .OR. KGF.LT.NOG) THEN

C******************************
C     begin computation of grid kc difference operator when kf difference
C     operator is twenty-seven point unless kc. ge. NOG
C

      K = KSTART

      DO KC = 2,KKC1

         K = K+2
         J = JSTART

         DO JC = 2,JJC1

            J = J+2
            I = ISTART
            IC = 1

            CALL BMG3_SymStd_UTILS_bbody27(IIF, JJF, KKF, IIC, JJC,
     &           KKC, IC, JC, KC, I, J, K, CI, SO, Q, R)

            DO IC= 2,IIC1

               I = I+2

               !  APPLY RESIDUAL WEIGHTING

               SOC(IC,JC,KC,KPNW)
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC, JC-1, KC, 6, 2, 4
     &              )
                
ce               SOC(IC,JC,KC,KBNW)
ce     &              = - BMG3_SymStd_UTILS_reswt(
ce     &              CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 6, 2, 2
ce     &              )
 
               SOC(IC,JC,KC,KBE)
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC, JC, KC-1, 6, 4, 2
     &              )
 
               SOC(IC,JC,KC,KBSE)
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 6, 2, 2
     &                  )

            CALL BMG3_SymStd_UTILS_bbody27(IIF, JJF, KKF, IIC, JJC,
     &           KKC, IC, JC, KC, I, J, K, CI, SO, Q, R)

ce               SOC(IC,JC,KC,KBN)
ce     &         = - BMG3_SymStd_UTILS_reswt(
ce     &                  CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 4, 2, 2 
ce     &                  )

               SOC(IC,JC,KC,KP) 
     &         =   BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC, KC, 4, 4, 4
     &                  )

               SOC(IC,JC,KC,KPW) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC-1, JC, KC, 2, 4, 4
     &                  )
               
               SOC(IC,JC,KC,KPS) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC-1, KC, 4, 2, 4 
     &                  )
               
               SOC(IC,JC,KC,KPSW) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC-1, JC-1, KC, 2, 2, 4 
     &                  )
               
               SOC(IC,JC,KC,KB) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC, KC-1, 4, 4, 2 
     &                  )
               
               SOC(IC,JC,KC,KBSW) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC-1, JC-1, KC-1, 2, 2, 2 
     &                  )
               
               SOC(IC,JC,KC,KBW) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC-1, JC, KC-1, 2, 4, 2 
     &                  )

               
               SOC(IC,JC,KC,KBS) 
     &         = -BMG3_SymStd_UTILS_reswt( 
     &                 CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 4, 2, 2
     &                 )

            ENDDO  
         ENDDO  
      ENDDO

C     SPECIAL CASE

      K = KSTART
 
      DO KC = 2,KKC1
 
         K = K+2

         J = JSTART

         DO JC = 2,JJC1

            J = J+2
            I = ISTART
            IC = 1
            
            CALL BMG3_SymStd_UTILS_bbody27(IIF, JJF, KKF, IIC, JJC,
     &           KKC, IC, JC-1, KC, I, J-2, K, CI, SO, Q, R)

            DO IC = 2,IIC1
               
               I = I+2

               SOC(IC,JC,KC,KBNE)
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC, JC, KC-1, 6, 6, 2
     &              )               

               CALL BMG3_SymStd_UTILS_bbody27(IIF, JJF, KKF, IIC, JJC,
     &              KKC, IC, JC-1, KC, I, J-2, K, CI, SO, Q, R)

               SOC(IC,JC,KC,KBN) 
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC, JC, KC-1, 4, 6, 2) 

               SOC(IC,JC,KC,KBNW)
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC-1, JC, KC-1, 2, 6, 2
     &              )                   
            ENDDO
          ENDDO   
        ENDDO

C     end of computation of kc difference operator when kf difference
C     operator is twenty-seven point
C******************************

      ELSE

C******************************
C   begin computation of kc difference operator when kf difference
C   operator is seven point unless kc.ge.NOG
C

      K = KSTART

      DO KC = 2,KKC1

         K = K+2
         J = JSTART

         DO JC = 2,JJC1

            J = J+2
            I = ISTART
            IC = 1
            CALL BMG3_SymStd_UTILS_bbody07(IIF, JJF, KKF, IIC, JJC,
     &           KKC, IC, JC, KC, I, J, K, CI, SO, Q, R)
            DO IC= 2,IIC1

               I = I+2

               !  APPLY RESIDUAL WEIGHTING

               SOC(IC,JC,KC,KPNW) 
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC, JC-1, KC, 6, 2, 4 
     &              )
               
ce               SOC(IC,JC,KC,KBNW)
ce     &              = - BMG3_SymStd_UTILS_reswt(
ce     &              CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 6, 2, 2 
ce     &              )

               SOC(IC,JC,KC,KBE) 
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC, JC, KC-1, 6, 4, 2 
     &              )

               SOC(IC,JC,KC,KBSE)
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 6, 2, 2 
     &                  )


            CALL BMG3_SymStd_UTILS_bbody07(IIF, JJF, KKF, IIC, JJC,
     &           KKC, IC, JC, KC, I, J, K, CI, SO, Q, R)

ce               SOC(IC,JC,KC,KBN)
ce     &         = - BMG3_SymStd_UTILS_reswt(
ce     &                  CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 4, 2, 2 
ce     &                  )

               SOC(IC,JC,KC,KP) 
     &         =   BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC, KC, 4, 4, 4
     &                  )

               SOC(IC,JC,KC,KPW) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC-1, JC, KC, 2, 4, 4
     &                  )
               
               SOC(IC,JC,KC,KPS) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC-1, KC, 4, 2, 4 
     &                  )
               
               SOC(IC,JC,KC,KPSW) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC-1, JC-1, KC, 2, 2, 4 
     &                  )
               
               SOC(IC,JC,KC,KB) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC, KC-1, 4, 4, 2 
     &                  )
               
               SOC(IC,JC,KC,KBSW) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC-1, JC-1, KC-1, 2, 2, 2 
     &                  )
               
               SOC(IC,JC,KC,KBW) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC-1, JC, KC-1, 2, 4, 2 
     &                  )

               
               SOC(IC,JC,KC,KBS) 
     &         = -BMG3_SymStd_UTILS_reswt( 
     &                 CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 4, 2, 2
     &                 )

           ENDDO
         ENDDO
      ENDDO   

C     SPECIAL CASE

      K = KSTART
 
      DO KC = 2,KKC1
 
         K = K+2

         J = JSTART
         
         DO JC = 2,JJC1

            J = J+2
            I = ISTART
            IC = 1
            
            CALL BMG3_SymStd_UTILS_bbody07(IIF, JJF, KKF, IIC, JJC,
     &           KKC, IC, JC-1, KC, I, J-2, K, CI, SO, Q, R)

            DO IC = 2,IIC1
               
               I = I+2

               SOC(IC,JC,KC,KBNE)
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC, JC, KC-1, 6, 6, 2
     &              )               

               CALL BMG3_SymStd_UTILS_bbody07(IIF, JJF, KKF, IIC, JJC,
     &              KKC, IC, JC-1, KC, I, J-2, K, CI, SO, Q, R)

               SOC(IC,JC,KC,KBN) 
     &         = - BMG3_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, KKC, IC, JC, KC-1, 4, 6, 2)

               SOC(IC,JC,KC,KBNW)
     &              = - BMG3_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, KKC, IC-1, JC, KC-1, 2, 6, 2
     &              )

            ENDDO
          ENDDO   
        ENDDO

      ENDIF

      DO kpz=1,14

            ptrn = 1 
            call MSG_tbdx_send(SOC(1,1,1,kpz), buffer, 
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),
     &           ptrn, ierror)
         
            call MSG_tbdx_receive(SOC(1,1,1,kpz), buffer,
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),
     &           ptrn, ierror)
         
            call MSG_tbdx_close(SOC(1,1,1,kpz), buffer,
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)), 
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),
     &           ptrn, ierror)



      ENDDO   

C ========================================================================

      RETURN
      END
