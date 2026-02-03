#include "w3macros.h"
!/ ------------------------------------------------------------------- /
    module w3cplshyfem

      USE W3GDATMD, ONLY: IOBP_LOC, IOBPD_LOC, IOBPA_LOC, IOBDP_LOC, CLATS
      USE W3GDATMD, ONLY: NK, NTH, NSEA, NSEAL, SIG, DTH, DSII, ESIN, ECOS
      USE W3GDATMD, ONLY: DDEN
      USE W3ADATMD, ONLY: SXX, SXY, SYY, WN, CG, THM, HS, T02
      USE W3WDATMD, ONLY: VA, TIME
      USE W3PARALL, only: INIT_GET_ISEA
      USE CONSTANTS, ONLY: TPI, TPIINV, GRAV

      IMPLICIT NONE
!
#ifdef W3_MPI
      INCLUDE "mpif.h"
#endif

!WW3
      REAL, ALLOCATABLE :: SXX3DLWW3(:,:)
      REAL, ALLOCATABLE :: SXY3DLWW3(:,:)
      REAL, ALLOCATABLE :: SYY3DLWW3(:,:)
      REAL, ALLOCATABLE :: OUTVARWW3(:,:)

!SHYFEM

      INTEGER, PARAMETER :: NVARSWW3 = 3
      INTEGER            :: NLVT

      contains
!/ ------------------------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
      subroutine w3cplinit(nlvdim)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    04-May-2005 : Origination.                        ( version 3.07 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    17-Feb-2016 : New version from namelist use       ( version 5.11 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Compute or make available all needed values for shyfem 
!     (uncoupled).
!
!  2. Method : 
!  
!     Radiation stresses in 2d and other quantities 
!
!  3. Parameters :
!
!     Local parameters.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!                Subr.          Standard MPI routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - This is he third version, version 1 and 2 were use for proof
!       of concept only, and were not retained.
!
!  8. Structure :
!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
    IMPLICIT NONE
!
#ifdef W3_MPI
    INCLUDE "mpif.h"
#endif
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
        INTEGER, INTENT(IN) :: NLVDIM
!
        NLVT = NLVDIM
! 
        ALLOCATE(OUTVARWW3(NVARSWW3,NSEAL)); OUTVARWW3 = 0.d0
        ALLOCATE(SXX3DLWW3(NLVT,NSEAL)); SXX3DLWW3 = 0.
        ALLOCATE(SXY3DLWW3(NLVT,NSEAL)); SXY3DLWW3 = 0.
        ALLOCATE(SYY3DLWW3(NLVT,NSEAL)); SYY3DLWW3 = 0.
 
      end subroutine w3cplinit
!/ ------------------------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
      subroutine w3cplrdstr2d
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    04-May-2005 : Origination.                        ( version 3.07 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    17-Feb-2016 : New version from namelist use       ( version 5.11 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Compute or make available all needed values for shyfem 
!     (uncoupled).
!
!  2. Method : 
!  
!     Radiation stresses in 2d and other quantities 
!
!  3. Parameters :
!
!     Local parameters.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!                Subr.          Standard MPI routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - This is he third version, version 1 and 2 were use for proof
!       of concept only, and were not retained.
!
!  8. Structure :
!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
    IMPLICIT NONE
!
#ifdef W3_MPI
    INCLUDE "mpif.h"
#endif
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
       REAL                :: RSXX(NSEAL), RSXY(NSEAL), RSYY(NSEAL), ACLOC(NK,NTH)
       REAL                :: COSE2, SINE2, COSI2, ELOC, WNL, SIGTPI(NK)

       INTEGER             :: ISEA, JSEA, IK, ITH, ISP
!/
!/ Code
!/
       RSXX = 0.
       RSXY = 0.
       RSYY = 0.

       SIGTPI = SIG * TPI 

       DO JSEA = 1, NSEAL
         CALL INIT_GET_ISEA(ISEA, JSEA)
         IF (ABS(IOBDP_LOC(JSEA)) .GT. 0) THEN
           DO ITH = 1, NTH
             DO IK = 1, NK
               ISP = ITH + (IK-1)*NTH
               COSE2 = ECOS(ITH)**2
               SINE2 = ESIN(ITH)**2
               COSI2 = ECOS(ITH) + ESIN(ITH)
               ACLOC(IK,ITH) = VA(ISP,JSEA) / CG(IK,ISEA) * CLATS(ISEA)
               ELOC         = ACLOC(IK,ITH) * DTH * DSII(IK)
               WNL          = CG(IK,ISEA) / ( SIGTPI(IK)/WN(IK,ISEA) )
               RSXX(JSEA)   = RSXX(JSEA) + ( WNL * COSE2 + WNL - 0.5) * ELOC   ! Units = [ 1/s + 1/s - 1/s ] * m²s = m²
               RSXY(JSEA)   = RSXY(JSEA) + ( WNL * COSI2            ) * ELOC
               RSYY(JSEA)   = RSYY(JSEA) + ( WNL * SINE2 + WNL - 0.5) * ELOC
             ENDDO
           ENDDO
         ELSE
           RSXX(JSEA) = 0.
           RSXY(JSEA) = 0.
           RSYY(JSEA) = 0.
         END IF
       END DO

       SXX3DLWW3 = 0.
       SXY3DLWW3 = 0.
       SYY3DLWW3 = 0.
       DO JSEA = 1, NSEAL
         IF (ABS(IOBDP_LOC(JSEA)) .GT. 0) THEN
           SXX3DLWW3(:,JSEA) = RSXX(JSEA)! * GRAV !ccf
           SXY3DLWW3(:,JSEA) = RSXY(JSEA)! * GRAV !ccf
           SYY3DLWW3(:,JSEA) = RSYY(JSEA)! * GRAV !ccf 
         ELSE
           SXX3DLWW3(:,JSEA) = 0.
           SXY3DLWW3(:,JSEA) = 0.
           SYY3DLWW3(:,JSEA) = 0.
         END IF
       END DO
!/
!/ ------------------------------------------------------------------- /
!/ End of w3cplrdstr2d------------------------------------------------ /
!/ ------------------------------------------------------------------- /
!/
      end subroutine w3cplrdstr2d
!/ ------------------------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
      subroutine w3cplparam
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    04-May-2005 : Origination.                        ( version 3.07 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    17-Feb-2016 : New version from namelist use       ( version 5.11 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Compute or make available all needed values for shyfem 
!     (uncoupled).
!
!  2. Method : 
!  
!     Integral wave parameter for shyfem 
!
!  3. Parameters :
!
!     Local parameters.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!                Subr.          Standard MPI routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - This is he third version, version 1 and 2 were use for proof
!       of concept only, and were not retained.
!
!  8. Structure :
!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/ ------------------------------------------------------------------- /
!/
      REAL                :: ETOT, FMEAN2, EB(NK)

      INTEGER             :: ISEA, JSEA, IK, ITH, ISP

      DO JSEA = 1, NSEAL 

        CALL INIT_GET_ISEA(ISEA, JSEA)

        ETOT = 0.
        FMEAN2 = 0.

        DO IK=1, NK
          EB(IK) = 0.
          DO ITH=1, NTH
            ISP = ITH+(IK-1)*NTH
            EB(IK) = EB(IK) + VA(ISP,JSEA)
          END DO
        END DO

        DO IK=1, NK
          EB(IK) = EB(IK) * DDEN(IK) / CG(IK,ISEA)
          ETOT  = ETOT  + EB(IK)
        END DO

        IF (ETOT .LT. TINY(1.)) THEN
          OUTVARWW3(:,JSEA) = 0
        ELSE
          DO IK=1, NK
            FMEAN2 = FMEAN2 + EB(IK) * SIG(IK)
          END DO
          FMEAN2 = FMEAN2 / ETOT * TPIINV
          OUTVARWW3(1,JSEA) = 4 * SQRT(ETOT)  ! Hs 
          OUTVARWW3(2,JSEA) = 1./FMEAN2       ! Tm02 
          OUTVARWW3(3,JSEA) = THM(JSEA)
        ENDIF 

      ENDDO ! NSEAL 

!
!/ ------------------------------------------------------------------- /
    end subroutine w3cplparam
!/ ------------------------------------------------------------------- /
!/
END MODULE 
!/
!/ ------------------------------------------------------------------- /
