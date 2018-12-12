	MODULE modu1
	IMPLICIT NONE
	SAVE
	PUBLIC

#include "ppexec_user.cmn"
      EQUIVALENCE (RMISS, USER_RUMISS)
      EQUIVALENCE (IMISS, USER_IUMISS)

      INTEGER IMISS
      REAL RMISS

C Commons for flash
# include "dms_stwkwk.cmn"
# include "shs_stwork.cmn"


#include "dms_ncomp.cmn"

#include "dms_plex.cmn"

#include "dms_maxwrt.cmn"

	REAL*8 B(1)
	EQUIVALENCE (B(1), IB(1))

	INTEGER, PARAMETER :: SGL=4, DBL=8
     

	!pre-exponential factors in m^3/(kmol*sec)
	REAL(KIND=8), PARAMETER :: k1f = 0.104
	REAL(KIND=8), PARAMETER :: k1r = 0.113


	!Density
	REAL(KIND=8) :: RHOV, RHOLI, RHOLII


	!Vector to save the inlet flow and fake inlet flow
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: MSIN_SAVE
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SOUT_SAVE


	!Parameters for Physical Properties

	INTEGER :: KERR 

	!More Paramteres for Physical Properties

	INTEGER, PARAMETER :: KBASE = 1


	END MODULE modu1
