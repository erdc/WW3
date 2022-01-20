#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SVEGMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           A. Roland               |
!/                  |           A. Abdolali             |
!/                  |           T. Hesser               |
!/                  |           J. Smith                |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         16-Sep-2019 |
!/                  +-----------------------------------+
!/
!/    16-Sep-2019 : Initial Version                  ( version 6.XX )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Vegetation dissipation following Mendez (2004)
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SDB1    Subr. Public   Mendez (2004) bulk veg. dissipation 
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     See subroutine documentation.
!
!  5. Remarks :
!
!  6. Switches :
!
!     See subroutine documentation.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SVEG1 (IX, IY, A, DEPTH, EMEAN, FMEAN, WNMEAN, S, D )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |                        FORTRAN 90 |
!/                  !           A. Roland               | 
!/                  |           A. Abdolali             |
!/                  |           T. Hesser               |
!/                  |           J. Smith                |
!/                  | Last update :         08-Jun-2018 |
!/                  +-----------------------------------+
!/
!/    19-Sep-2019 : Origination of module.              ( version 3.11 )
!/
!  1. Purpose :
!
!     Compute bulk vegetation dissipation following Mendez (2004)
!
!  2. Method : 
! 
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D)
!       EMEAN   Real  I   Mean wave energy.
!       FMEAN   Real  I   Mean wave frequency.
!       WNMEAN  Real  I   Mean wave number.
!       DEPTH   Real  I   Mean water depth.
!       S       R.A.  O   Source term (1-D version).
!       D       R.A.  O   Diagonal term of derivative (1-D version).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Subroutine tracing (!/S switch).
!
!  5. Called by :
!
!       W3SRCE   Source term integration.
!       W3EXPO   Point output post-processor.
!       GXEXPO   GrADS point output post-processor.
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/Tn  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE CONSTANTS
      USE W3SERVMD, ONLY: EXTCDE
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, FSSOURCE, DDEN
      USE W3ODATMD, ONLY: NDST
      USE W3GDATMD, ONLY: SIG
      USE W3ODATMD, only : IAPROC, NAPERR, NDSE
      USE W3IDATMD, ONLY: INFLAGS2, VEGLS, VEGBV, VEGN, VEGCD
#ifdef W3_PDLIB
   use yowNodepool, only: ipgl, iplg
#endif
      USE CONSTANTS, ONLY: GRAV, PI, TPI, LPDLIB 
#ifdef W3_S
      USE W3SERVMD, ONLY: STRACE
#endif
#ifdef W3_T0
      USE W3ARRYMD, ONLY: PRT2DS
#endif
#ifdef W3_T1
      USE W3ARRYMD, ONLY: OUTMAT
#endif
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IX, IY ! Local grid number
      REAL, INTENT(IN)        :: A(NSPEC)
      REAL, INTENT(INOUT)     :: EMEAN, FMEAN, WNMEAN, DEPTH
      REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
      REAL                    :: CVEG
      INTEGER                 :: ITH, IK, IWB
      REAL, PARAMETER         :: KDMAX = 20.
      REAL                    :: KBAR 
      REAL(8), PARAMETER :: TWO  = 2.0D0
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
#ifdef W3_T0
     INTEGER                 :: IK, ITH
#endif
      REAL                    :: ETOT, FMEAN2, TMP1, TMP2
      REAL                    :: VCD, VDM, VNV, VLTH
      REAL                    :: VALPHAP   
      REAL                    :: VALPHAD   
      REAL                    :: VALPHADH 
      REAL                    :: THR, KSBAR, SBAR
#ifdef W3_T0
     REAL                   :: DOUT(NK,NTH)
#endif
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'W3SDB1')
#endif
!
! 0.  Initialzations ------------------------------------------------- /

      VCD = 0.0
      VDM = 0.0
      VNV = 0.0
      VLTH = 0.0
      S = 0.
      D = 0.

      IF (INFLAGS2(7)) THEN
!!         IF (LPDLIB) THEN
!!           VLTH = VEGLS(iplg(IX),IY)
!!           VDM  = VEGBV(iplg(IX),IY)
!!           VNV  = VEGN(iplg(IX),IY)
!!         ELSE
           VLTH = VEGLS(IX,IY)
           VDM  = VEGBV(IX,IY)
           VNV  = VEGN(IX,IY)
           VCD = VEGCD(IX,IY)
!!         ENDIF
      ELSE
         IF ( IAPROC .EQ. NAPERR )                &
         WRITE (NDSE,1001) 'VEGETATION PARAMETERS'
         CALL EXTCDE(2)
      END IF
      

#ifdef W3_PDLIB
 
 
#endif


      THR = DBLE(1.E-15)
      IF (EMEAN .LT. THR) RETURN
!       
      IF (FMEAN .GT. THR) THEN ! Actually this cannot happen if the above is true we make it now double sure ...
        KBAR = WNMEAN 
        SBAR = FMEAN * TPI
        KSBAR = KBAR/SBAR ! Check units with 
      ELSE
        RETURN
      ENDIF
!
#ifdef W3_T
      WRITE (NDST,9000) SDBC1, SDBC2, FDONLY
#endif
!
      !VALPHAP   = VDM*VNV*VCD/TWO
       VALPHAP   = VDM*VNV*VCD !3 for check
      IF (DEPTH .LT. VLTH) THEN
        VALPHAD = VLTH/DEPTH   ! Submerged Case
      ELSE
        VALPHAD = 1.0     ! Emergent case
      ENDIF 
!  alpha or VLTH is veg length/depth for submerged veg and  equal to 1 for emergent
      TMP1 = SINH(MIN(KDMAX,WNMEAN*VALPHAD*DEPTH))**3+3*SINH(MIN(KDMAX,WNMEAN*VALPHAD*DEPTH))  !Recoded to account ffor both submerged and emergent cases TJH
      TMP2 = 3*WNMEAN*COSH(MIN(KDMAX,WNMEAN*DEPTH))**3  !added the missing k in this 3*k*cosh(kh)^3
      CVEG = SQRT(TWO/PI)*GRAV**2*VALPHAP*KSBAR**3*TMP1/TMP2*SQRT(EMEAN)
      D    = - CVEG
      S    = - CVEG * A
      !write(*,'(A10,10F20.10)') 'CVEG', -CVEG, tmp1, tmp2, depth, wnmean, fmean, emean
!
! ... Test output of arrays
!
#ifdef W3_T0
      DO IK=1, NK
        DO ITH=1, NTH
          DOUT(IK,ITH) = D(ITH+(IK-1)*NTH)
          END DO
        END DO
#endif
!
#ifdef W3_T0
      CALL PRT2DS (NDST, NK, NK, NTH, DOUT, SIG, '  ', 1.,    &
                         0.0, 0.001, 'Diag Sdb', ' ', 'NONAME')
#endif
!     
#ifdef W3_T1
      CALL OUTMAT (NDST, D, NTH, NTH, NK, 'diag Sdb')
#endif
!  
      RETURN
!
! Formats   
 1001 FORMAT(/' *** WAVEWATCH III ERROR IN W3SVEG1MD : '/   &
              '     ',A,' IS NOT DEFINED IN ww3_shel.inp.')
!
#ifdef W3_T
 9000 FORMAT (' TEST W3SDB1 : PARAMETERS :',2F7.3,L4)
#endif
!/
!/ End of W3SVEG1 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SVEG1
!/
!/
!/ End of module W3SVEG ---------------------------------------------- /
!/
      END MODULE W3SVEGMD
