!> @file
!> @brief MADSEN bottom friction source term (Madsen et al 1988).
!>
!> @author H. Michaud
!> @author M. Pezerat
!> @date   23-07-2024
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief MADSEN bottom friction source term (Madsen et al 1988).
!>
!> @details
!>
!> @date   14-Mar-2012
!> @author H. Michaud
!> @author M. Pezerat
!> @date   23-07-2024
!>
!> @copyright Copyright 2009-2022 National Weather Service (NWS),
!>       National Oceanic and Atmospheric Administration.  All rights
!>       reserved.  WAVEWATCH III is a trademark of the NWS.
!>       No unauthorized use without permission.
!>
MODULE W3SBT5MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |    H. Michaud and M. Pezerat      |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         23-07-2024  |
  !/                  +-----------------------------------+
  !/
  !/    23-Jul-2024 : Origination.                        ( version 7.XX )
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Madsen bottom friction source term (Madsen et al. 1988),
  !     adapted to rocky bottom and with coefficient from Sous et al. (2023).
  !
  !  2. Variables and types :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SBT5    Subr. Public   Rocky bottom friction 
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !  5. Remarks :
  !
  !     WAVEWATCH III is designed as a highly plug-compatible code.
  !     Source term modules can be included as self-contained modules,
  !     with limited changes needed to the interface of routine calls
  !     in W3SRCE, and in the point postprocessing programs only.
  !     Codes submitted for inclusion in WAVEWATCH III should be
  !     self-contained in the way described below, and might be
  !     provided with distributions fully integrated in the data
  !     structure, or as an optional version of this module to be
  !     included by the user.
  !
  !     Rules for preparing a module to be included in or distributed
  !     with WAVEWATCH III :
  !
  !      - Fully document the code following the outline given in this
  !        file, and according to all other WAVEWATCH III routines.
  !      - Provide a file with necessary modifications to W3SRCE and
  !        all other routines that require modification.
  !      - Provide a test case with expected results.
  !      - It is strongly recommended that the programming style used
  !        in WAVEWATCH III is followed, in particular
  !          a) for readability, write as if in fixed FORTRAN format
  !             regarding column use, even though all files are F90
  !             free format.
  !          b) I prefer upper case programming for permanent code,
  !             as I use lower case in debugging and temporary code.
  !
  !     This module needs to be self-contained in the following way.
  !
  !      a) All saved variables connected with this source term need
  !         to be declared in the module header. Upon acceptance as
  !         permanent code, they will be converted to the WAVEWATCH III
  !         dynamic data structure.
  !      b) Provide a separate computation and initialization routine.
  !         In the submission, the initialization should be called
  !         from the computation routine upon the first call to the
  !         routine. Upon acceptance as permanent code, the
  !         initialization routine will be moved to a more appropriate
  !         location in the code (i.e., being absorbed in ww3_grid or
  !         being moved to W3IOGR).
  !
  !     See notes in the file below where to add these elements.
  !
  !  6. Switches :
  !
  !     !/S  Enable subroutine tracing.
  !
  !  7. Source code :
  !/
  !/ ------------------------------------------------------------------- /
  !/
  !

  PUBLIC
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Computes the Rockey bottom friction parameterization adapted from 
  !> (Madsen et al. 1988) and with coefficient from Sous et al. (2023).
  !> @param[in]    A        Action density spectrum.
  !> @param[in]    CG       Group velocities.
  !> @param[in]    WN       Wavenumbers.
  !> @param[in]    DEPTH    Water depth.
  !> @param[in]    KKR      Hydraulic roughness height
  !> @param[out]   TAUBBL   Components of stress leaking to the bottom.
  !> @param[out]   S        Source term (1-D version).
  !> @param[out]   D        Diagonal term of derivative.
  !> @param[in]    IX       Spatial grid index.
  !> @param[in]    IY       Spatial grid index.
  !>
  !>
  SUBROUTINE W3SBT5 (A, CG, WN, DEPTH, KKR, TAUBBL, S, D, IX, IY )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |            H. Michaud             |
    !/                  |            M. Pézérat             |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         05-Aug-2025 |
    !/                  +-----------------------------------+
    !/
    !/    23-Jun-2024 : Origination.                        ( version 7.14 )
    !/
    !  1. Purpose :
    !
    !     Computes the MADSEN bottom friction adpated to rocky bottom
    !
    !  2. Method :
    !
    !                            1              SIG(IK)**3 * E(k,theta)    
    !       Sbt (k,theta)= -  ------- fw Uorb -----------------------    (1)
    !                         SQRT(2)          GRAV * sinh²(k*DEPTH)
    !
    !     with fw = EXP(a1*(UORB/(KKR*SIG(IK)))**(a2)+a3)
    !     Where a1 = 5
    !           a2 = -0.15    
    !           a3 = -5.9
    !
    !
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       A       R.A.  I   Action density spectrum.
    !       CG      R.A.  I   Group velocities.
    !       WN      R.A.  I   Wavenumbers.
    !       DEPTH   Real  I   Water depth.
    !       KKR     Real  I   Hydraulic roughness height
    !       TAUBBL  Real  O   Components of stress leaking to the bottom.
    !       S       R.A.  O   Source term (1-D version).
    !       D       R.A.  O   Diagonal term of derivative.             *)
    !       IX,IY   Int. I   Spatial grid indices
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD Subroutine tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3SRCE    Subr. W3SRCEMD Source term integration.
    !      W3EXPO    Subr.   N/A    Point output post-processor.
    !      GXEXPO    Subr.   N/A    GrADS point output post-processor.
    !     ----------------------------------------------------------------
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
    !     !/S  Enable subroutine tracing.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

    USE CONSTANTS
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DDEN,  &
          ECOS, ESIN, DTH

#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    !/
    LOGICAL, SAVE           :: FIRST = .TRUE.
    REAL, INTENT(IN)        :: CG(NK), WN(NK), DEPTH, A(NSPEC), KKR
    INTEGER, INTENT(IN)     :: IX, IY
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC), TAUBBL(2)
    REAL                    :: CBETA(NK)
    REAL :: UORB2, UORB, AORB, EBX, EBY, AX, AY, LX, LY
    REAL :: CONST2, TEMP2
    REAL :: FWJ, KSUBN, KSUBS, KSUBR, MINADIM
    REAL :: SHIELDS(3), PSI, DELI1, DELI2, EB, XI, VARU, KKKR 
    INTEGER :: IK, ITH, IS, IND, INDE, ISUB
    REAL :: DSUM(NK)
    REAL ::  FACTOR

    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'W3SBT5')
#endif
    !

    !
    ! 1.  Min / Max settings for KKR ---- *
    !
    KKKR=MAX(KKR,1E-5)                
    KKKR=MIN(KKKR,1.)                  
    !
    !
    DSUM(:)=0.
    TAUBBL(:)=0.
    !
    UORB=0.
    AORB=0.
    AX  =0.
    AY  =0.

    DO IK=1, NK
      IF ( WN(IK)*DEPTH .LT. 6. ) THEN
        EB  = 0.
        EBX = 0.
        EBY = 0.
        DO ITH=1, NTH
          IS=ITH+(IK-1)*NTH
          EB  = EB  + A(IS)
          EBX = EBX +A(IS)*ECOS(ITH)
          EBY = EBY +A(IS)*ESIN(ITH)
        END DO
          !
          ! U_bot=sigma * Zeta / sinh(KD)  and CBETA = 0.5*sigma^2 /(g*sinh^(kD))
          ! therefore variance(u_bot)= variance(elevation)*2*CBETA/D
          !
        CBETA(IK) = 1/(SQRT(2.0))*SIG(IK)**2 /(GRAV*(SINH(WN(IK)*DEPTH))**2)
        
          !  N.B.:  could also include shoaling effect on EB ...
        FACTOR= (DDEN(IK)/CG(IK))*SQRT(2.0)*CBETA(IK)*GRAV/CG(IK)
        VARU= EB * FACTOR
        UORB = UORB + VARU
        AORB = AORB + VARU/(SIG(IK)**2)
        AX   = AX   + (EBX * FACTOR)
        AY   = AY   + (EBY * FACTOR)
      ELSE
        CBETA(IK) = 0.
      END IF
    END DO
      !
      ! Computes RMS orbital amplitudes
      !
    UORB2 = 2*UORB
    UORB = SQRT(MAX(1.0E-7,UORB2))
    AORB = SQRT(MAX(1.0E-7,2*AORB))

    !
    ! 2. Fills output arrays and estimates the energy and momentum loss
    !
    DO IK=1, NK
      ! New parametrisation (Madsen - SWAN)
      FWJ = EXP(5*(UORB/(KKKR*SIG(IK)))**(-0.15)-5.9)
      CONST2=DDEN(IK)/CG(IK) &         !Jacobian to get energy in band
           *GRAV/(SIG(IK)/WN(IK))    ! coefficient to get momentum
      DSUM(IK)=-FWJ*UORB*CBETA(IK)*SIG(IK)   
      DO ITH=1,NTH
        IS=ITH+(IK-1)*NTH
        D(IS)=DSUM(IK)
        TEMP2=CONST2*D(IS)*A(IS)
        TAUBBL(1) = TAUBBL(1) - TEMP2*ECOS(IS)
        TAUBBL(2) = TAUBBL(2) - TEMP2*ESIN(IS)
        S(IS)=D(IS)*A(IS)
      END DO
    END DO
    !
    RETURN
    !/
    !/ End of W3SBT5  ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SBT5
  !/
  !/ End of module W3SBT5MD -------------------------------------------- /
  !/
END MODULE W3SBT5MD
