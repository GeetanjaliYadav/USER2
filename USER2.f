      SUBROUTINE USER2 (NMATI,  MSIN,    NINFI,   SINFI,  NMATO,
     +                  SOUT,   NINFO,  SINFO,   IDSMI,  IDSII,
     +                  IDSMO,  IDSIO,  NTOT,    NSUBS,  IDXSUB,
     +                  ITYPE,  NINT,   INT,     NREAL,  REALs,
     +                  IDS,    NPO,    NBOPST,  NIWORK, IWORK,
     +                  NWORK,  WORK,   NSIZE,   SIZE,   INTSIZ,
     +                  LD )

C     This subroutine is a kinetic reactor for the transesterification of
C     triglycerides to FAME and calculation of phase equilibrium
      
C     This version of USER2 is for a block with one inlet stream
C     and one outlet stream.  Eventually, after the FLS3 subroutine is 
C     used two additional outlet streams will be added.  Then the six species 
C     flow rates that solve the six mass balances (computed by BROY)
C     will be input to FLS3 to obtain the species flow rates of the 
C     V, LI, and LII streams

      USE MODU1
      IMPLICIT NONE

C      EXTERNAL :: TR
      
C     Declare variables in dimensioning
      INTEGER NMATI, NINFI, NMATO, NINFO, NTOT,
     +        NSUBS, NINT,  NPO,   NIWORK,NWORK,
     +        NSIZE, JRES, MXIT, NPKODE, KODE, KPHASE
C
      INTEGER IDSMI(2,NMATI),      IDSII(2,NINFI),
     +        IDSMO(2,NMATO),      IDSIO(2,NINFO),
     +        IDXSUB(NSUBS),ITYPE(NSUBS), INT(NINT),
     +        IDS(2,3), NBOPST(6,NPO), IRETN(6),
     +        IWORK(NIWORK),INTSIZ(NSIZE),NREAL, LD
     
      REAL*8 MSIN(NTOT,NMATI),     SINFI(NINFI),
     +       SOUT(NTOT,NMATO),    SINFO(NINFO), 
     +       WORK(NWORK),  SIZE(NSIZE), REALs(NREAL)
      REAL*8 T,     XXX,  
     +       x(50), y(50), dvmx, idx(50), x1(50), x2(50),
     +       xu(50), yu(50), x1u(50), x2u(50),L1V, L2V, VV,VL,
     +       RETN(151), TOL, P, Q, TEST, FTOT, VFLOW, FVTOT,
     +       FL1TOT, FL2TOT 
      REAL*8 FVO(7), FLIO(7), FLIIO(7)
      REAL*8 FV(7), FLI(7), FLII(7)
      REAL(KIND=8), DIMENSION(6) :: rxn
      integer n, ker, ncp, isub, i, krslt, kdiag, nvec, EPS, K, KPE
      PARAMETER(EPS = 1e-07)

C
      INTEGER  :: DMS_KFORMC
      INTEGER :: DMS_IFCMNC
      INTEGER :: i_err
      REAL*8, DIMENSION(:), ALLOCATABLE :: MW !Molecular weight
      real*8 f(7)


      call tr (f,n)


!Vector to save the inlet flow and fake inlet flow
      
C    Take flow rates of feed stream to USER2 block
!    Define initial values for all variables

      MSIN_SAVE(:,:) = MSIN
      SOUT_SAVE = 0.0D0
      SOUT_SAVE(:,1) = MSIN(:,1)
     
c      LMW = DMS_IFCMNC('MW')

c      DO i=1,NCOMP_NCC
c      MW(i)=B(LMW+i)
c      END DO
c
C Initialize outlet stream by copying over component flows and the rest
C  of the stream vector containing the temperature, pressure, etc.
C===========================================================================
C
C     FIRST COPY FIRST INLET TO FIRST OUTLET
C
      FTOT=0
      DO 100 I = 1, NCOMP_NCC
c        SOUT(I,1) = MSIN(I,1)
         SOUT(I,1) = F(i)
         FTOT=FTOT+F(i)
        WRITE (MAXWRT_MAXBUF, *) 'SOUT ', i, sout(i,1)
        CALL DMS_WRTTRM(1) 
 100  CONTINUE

         SOUT(NCOMP_NCC+1,1) = FTOT
C
      DO 200 I = NCOMP_NCC+2, NCOMP_NCC+9
        SOUT(I,1) = 0
        SOUT(I,2) = 0
        SOUT(I,3) = 0
 200  CONTINUE 
C
C     INITIALIZE THE OTHER OUTLETs
C
      DO 300 I = 1, NCOMP_NCC+1
        SOUT(I,2) = 0D0
        SOUT(I,3) = 0D0
 300  CONTINUE
C
          
      
C     !Solves mass balances and residuals using TR and BROY
C     !Mass and Energy Balances

         
        WRITE (MAXWRT_MAXBUF, *) 'Writing to Terminal'
        CALL DMS_WRTTRM(1)   
        
C      READ-ME!!!!
C      INITIALIZATION of rates and flows
C      First rates of reactions  
        rxn(1) = 0
        rxn(2) = 0
        rxn(3) = 0
        rxn(4) = 0
        rxn(5) = 0
        rxn(6) = 0 
      
C     INITIALIZE flows in reactions
       FVO(i) = 0
       FLIO(i) = 0
       FLIIO(i) = 0
       
C===========================================================================
C Flash parameters common to each flash type
C  Always calculate results (KRSLT=2), Disable Restart (Jres=0)
C  Maximum number of iterations (MAXIT), Flash Tolerance (TOL)
C  Perform a vapor-liquid flash (NPKODE=2)
C===========================================================================

      KRSLT= 2
      JRES=0
      MXIT=100
      TOL=1D-4
      NPKODE=3

C===============================================================================
C Perform TP flash (KODE=2) to determine outlet vapor frac and enthalpy
C===============================================================================

      KODE= 2

      T= MSIN(NCOMP_NCC+2,1)
      SOUT(NCOMP_NCC+2,1) = T
      SOUT(NCOMP_NCC+2,2) = T
      SOUT(NCOMP_NCC+2,3) = T
      
      P=MSIN(NCOMP_NCC+3,1)
      SOUT(NCOMP_NCC+3,1) = P
      SOUT(NCOMP_NCC+3,2) = P
      SOUT(NCOMP_NCC+3,3) = P

C     Find outlet conditions for given T and P
C     Dummy variable xxx used for guess since this case does not have one.
      CALL FLSH_FLASH (SOUT ,NSUBS ,IDXSUB ,ITYPE ,NBOPST ,KODE
     +      ,NPKODE, KPHASE,MXIT  ,TOL    ,T, P, XXX, USER_LMSG, 
     +     USER_LPMSG ,JRES, KRSLT ,RETN  ,IRETN  ,USER_ICONVG)

C     Outlet heat stream equals flash duty plus inlet heat stream.
C     Stream enthalpy= molar flow * mass enthalpy * MW
      SINFO(1)=SINFI(1) + MSIN(NCOMP_NCC+1,1)*MSIN(NCOMP_NCC+4,1)*
     +      MSIN(NCOMP_NCC+9,1) - SOUT(NCOMP_NCC+1,1)*
     +      SOUT(NCOMP_NCC+4,1)*SOUT(NCOMP_NCC+9,1)


C Molar volume from flash
	WRITE (MAXWRT_MAXBUF, *) 'From Flash '
	CALL DMS_WRTTRM(1)
	WRITE (MAXWRT_MAXBUF, *) 'stwork_vv ',stwork_vv
	CALL DMS_WRTTRM(1)
	WRITE (MAXWRT_MAXBUF, *) 'stwork_vl ',stwork_vl
	CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *) 'stwork_vl1 ',stwork_vl1        !READ-ME!!!
	CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *) 'stwork_vl2 ',stwork_vl2        !READ-ME!!!
	CALL DMS_WRTTRM(1)

C  Molar molecular weight from flash

	WRITE (MAXWRT_MAXBUF, *) 'From Flash '
	CALL DMS_WRTTRM(1)
	WRITE (MAXWRT_MAXBUF, *) 'stwork_xmwv ',stwork_xmwv     !READ-ME!!!
	CALL DMS_WRTTRM(1)
	WRITE (MAXWRT_MAXBUF, *) 'stwork_xmwl ',stwork_xmwl     !READ-ME!!!
	CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *) 'stwork_xmwl1 ',stwork_xmwl1   !READ-ME!!!
	CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *) 'stwork_xmwl2 ',stwork_xmwl2   !READ-ME!!!
	CALL DMS_WRTTRM(1)
      
c Compositions

C
      kdiag=4
C
	T		= stwkwk_tcalc !outlet temperature (K)
	P		= stwkwk_pcalc !Outlet pressure (N/m2)
	ncp		= stwkwk_ncpmoo !No. of packed componetns in MIXED substream (excluding solids)
      VL      = stwork_VL    !Liquid volume (m3/kgmole)       READ-ME!!!
      L1V     = stwork_VL1   !Liquid 1 volume (m3/kgmole)     READ-ME!!!
      L2V     = stwork_VL2   !Liquid 2 volume (m3/kgmole)     READ-ME!!!
      VV      = stwork_VL    !Vapor volume (m3/kgmole)        READ-ME!!!
C
c   
	WRITE (MAXWRT_MAXBUF, *) 'x, x1, x2, y, idx '
	CALL DMS_WRTTRM(1)
c  
	DO n=1,ncp
	  Y(n)		= B(STWKWK_LRSTW+STWORK_MY+n-1)
	  X(n)		= B(STWKWK_LRSTW+STWORK_MX+n-1)
	  X1(n)		= B(STWKWK_LRSTW+STWORK_MX1+n-1)
	  X2(n)		= B(STWKWK_LRSTW+STWORK_MX2+n-1)
	  IDX(n)	= IB(STWKWK_LISTW+STWORK_MIM+n-1)

      WRITE (MAXWRT_MAXBUF, *) x(n), x1(n), x2(n)   ! y, x, x1, x2 are (PACKED) here mole fraction vectors for vapor, liquid, liq1 and liq2 phase
      CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *) y(n), idx(n)
      CALL DMS_WRTTRM(1)

        END DO 

c
	WRITE (MAXWRT_MAXBUF, *) 'unpack',ncp, ncomp_ncc
	CALL DMS_WRTTRM(1)
c      
	WRITE (MAXWRT_MAXBUF, *) 'x, y, x1, x2 '
	CALL DMS_WRTTRM(1)
	
	DO 10 I = 1, NCOMP_NCC
         xu(I)=0.D0
         x1u(I)=0.D0
         x2u(I)=0.D0
         yu(I)=0.D0
   10 CONTINUE
      DO 20 i = 1, NCP
         xu(IDX(i))=x(i)   !yu, xu, x1u and x2 u are the UNPACKED mole fraction vectors for vapor, liq, liq1 and liq2 phase
         x1u(IDX(i))=x1(i)
         x2u(IDX(i))=x2(i)
         yu(IDX(i))=y(i)
   20 CONTINUE
	
      DO 500 i = 1, NCOMP_NCC
      SOUT(I,1) = MSIN(NCOMP_NCC+1,1)*STWKWK_VCALC*yu(i)   !stwkwk_vcalc = outlet molar vapor fraction or vapor fraction (molar)
C     SOUT(I,1) is the total flow rate of vapor phase after Flash operation        
      SOUT(I,2) = MSIN(NCOMP_NCC+1,1)*                      !stwkwk_BETA = Liq1/Total liq (molar ratio)
     +               (1-STWKWK_VCALC)*stwkwk_BETA*x1u(i)
C     SOUT(I,2) is the total flow rate of Liquid 1 phase after Flash operation           
      SOUT(I,3) = MSIN(NCOMP_NCC+1,1)*                       ! stwkwk_BETA = Liq2/Total liq (molar ratio)  i.e. (1-stwkwk_BETA = for the 2nd liq)
     +               (1-STWKWK_VCALC)*(1-stwkwk_BETA)*x2u(i)
C     SOUT(I,3) is the total flow rate of Liquid 2 phase after Flash operation   
c
      WRITE (MAXWRT_MAXBUF, *) xu(i), yu(i)
      CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *) x1u(i), x2u(i)
      CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Component', i
      CALL DMS_WRTTRM(1) 
      WRITE (MAXWRT_MAXBUF, *) 'SOUT(I,1)= vapor flow', i, SOUT(I,1)
      CALL DMS_WRTTRM(1) 
      WRITE (MAXWRT_MAXBUF, *) 'SOUT(I,2)= Liq1 flow', i, SOUT(I,2)
      CALL DMS_WRTTRM(1) 
      WRITE (MAXWRT_MAXBUF, *) 'SOUT(I,3)= Liq2 flow', i, SOUT(I,3)
      CALL DMS_WRTTRM(1)
 500  CONTINUE 
      
C     Here, we compute the total phase flow rates after Flash operation       
      SOUT(NCOMP_NCC+1,1) = MSIN(NCOMP_NCC+1,1)*STWKWK_VCALC
      FVTOT = SOUT(NCOMP_NCC+1,1)          !READ-ME!!! Flashed total flow rate of the vapor phase
      SOUT(NCOMP_NCC+1,2) = MSIN(NCOMP_NCC+1,1)*
     +                      (1-STWKWK_VCALC)*stwkwk_BETA
      FL1TOT = SOUT(NCOMP_NCC+1, 2)        ! READ-ME!!! Flashed total flow rate of the Liquid 1 phase   
      SOUT(NCOMP_NCC+1,3) = MSIN(NCOMP_NCC+1,1)*
     +                      (1-STWKWK_VCALC)*(1-stwkwk_BETA)
      FL2TOT = SOUT(NCOMP_NCC+1,3)         ! READ-ME!!!  Flashed total flow rate of the Liquid 2 phase
     
      WRITE(MAXWRT_MAXBUF,*) 'Component', i
      CALL DMS_WRTTRM(1) 
      WRITE (MAXWRT_MAXBUF, *) 'vapor phase flow = FVTOT', i, FVTOT
	CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *) 'Liq1 phase flow = FL1TOT', i, FL1TOT
	CALL DMS_WRTTRM(1)
      WRITE (MAXWRT_MAXBUF, *) 'Liq2 phase flow = FL2TOT', i, FL2TOT
	CALL DMS_WRTTRM(1)
      
C     READ-ME!!! (Here, we write the flashed out flow rates of each species (in all three phases) into vector, FVO, FLIO and FLIIO)      
      FVO(i) = SOUT(i,1)
      FLIO(i) = SOUT(i,2)
      FLIIO(i) = SOUT(i,2)
      
      Call TR(FVO, SOUT)           !READ-ME- I add the TR statements in here becuase here I want to check for the convergence and if it is not attained, will take the new Flow rate values in the loop untill convergence is acheived.
      
      IF (K .EQ. 0 .OR. K .EQ. KPE) THEN
      
C     Check for convergence
      Do i = 1, 7     
      IF (((ABS(FV(i)-FVO(i))/FV(i)) .LT. EPS) .AND.     
     +     ((ABS(FLI(i)-FLIO(i))/FLI(i)) .LT. EPS) .AND.
     +       ((ABS(FLII(i)-FLIIO(i))/FLII(i)) .LT. EPS)) THEN
     
       WRITE(MAXWRT_MAXBUF,*) 'Convergence is achieved'
       CALL DMS_WRTTRM(1)
     
      ELSE 
     
C     READ -ME!!! Set FVO, FLIO and FLIIO as new guess values.
      ! INITIALISE
      FV(i) = 0
      FLI(i) = 0
      FLII(i) = 0
          
      FV(i) = FVO(i)
      FLI(i) = FLIO(i)
      FLII(i) = FLIIO(i)
      
      GO TO 103
      END IF
         ENDDO
          ELSE 
            DO i = 1,7
                  FV(i) = FVO(i)
                  FLI(i) = FLIO(i)
                  FLII(i) = FLIIO(i)
              END DO    
              K = K + 1
              GO TO 103
              ENDIF
     
      
      CONTINUE
     
C ======================================================================
C Report Writer Block, only enter on report pass
C===============================================================================

      IF(USER_IPASS .GT. 3) THEN

C===============================================================================
C Report Header Utility, page 4-2 of User Models Reference Manual
C     Arguments:
C        1. Number of lines used in REP file by this block
C        2. Pagination flag, 1=start new page, 0=continue output on same page.
C        3. Desired section of report to contain output, 3=unit operations.
C        4. ISUB is the subsection heading defined by this block.  The
C           TITLE character string is used in this example to specify the
C           title with a DATA statement.
C===============================================================================

              CALL ZREP_RPTHDR(5,0,3,ISUB)

        WRITE(USER_NRPT, 600, ERR=610) SOUT(NCOMP_NCC+
     +        2,1)*1.8-459.67, SOUT(NCOMP_NCC+5,1), SOUT(NCOMP_NCC+3,1)/
     +        1D5*14.5038, SOUT(NCOMP_NCC+4,1)/2326.


610     CONTINUE

      END IF

600     FORMAT(T5,'TP Flash: ',/,
     +         T5,'Outlet Temperature (F): ',F10.1,/,
     +         T5,'Vapor Fraction:         ',F10.3,/,
     +         T5,'Pressure (psia):        ',F10.1,/,
     +         T5,'Enthalpy (Btu/Lb):      ',F10.1       )  

      WRITE(*,*)"Leaving USER2 Subroutine"

      RETURN
      END 



