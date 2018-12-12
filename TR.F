      SUBROUTINE TR(F,n)
C    This function evaluates the mass balance of flow rate 
C    It takes in the flow rate vector (F) and
C    outputs the new flow rates
      
      IMPLICIT NONE
#include "dms_ncomp.cmn"
#include "dms_maxwrt.cmn"
C      Inputs and Outputs
	  
      INTEGER  i, imax, ISTAT, KPE, K, EPS
      PARAMETER(imax=7)
      PARAMETER(KPE = 10) 
      PARAMETER(EPS = 1E-7)
      INTEGER, INTENT(IN) :: n
      REAL*8 F0(imax)
      REAL*8 F(imax)
      REAL*8 TAUREA
      PARAMETER(TAUREA = 3600)
      
      REAL*8 CLI(imax), CLII(imax)
      REAL*8 VV, VLI, VLII
      REAL*8 MWVA, MWLIA, MWLIIA
      REAL*8 MW(imax)
      REAL*8 RHOV, RHOLI, RHOLII
      REAL*8 FVT, FLIT, FLIIT
      REAL*8 FV(imax), FLI(imax), FLII(imax)
      REAL*8 XV(imax), XLI(imax), XLII(imax) 
      REAL*8 FVO(imax), FLIO(imax), FLIIO(imax)

C     Reactions
	REAL(KIND=8), DIMENSION(6) :: rxn
!      REAL*8 rxn(imax) 
      INTEGER :: i_err
      


C     Reaction constants in (m^3/kgmol*s)
      REAL*8 k1f                       ! All these values will be put in Aspen File .bkp user arrays
      PARAMETER(k1f = 0.104)
      REAL*8 k1r
      PARAMETER(k1r = 0.113)
      REAL*8 cat
      PARAMETER(cat = 0.003904)
      REAL*8 ncatf
      PARAMETER(ncatf = 1.26869)
      REAL*8 ncatr
      PARAMETER(ncatr = 0.699638)
      
C     These are our component list as used in Aspen Plus (Species Counter)     
C     F(1) = TG
C     F(2) = DG
C     F(3) = MG
C     F(4) = FAME
C     F(5) = GLY
C     F(6) = MEOH
C     F(7) = CO2
      
C     Specify initial Feed flow rates - F0(1), F0(6)
C     These will be input in the Aspen Plus as feed stram values in kmol/sec
      F0(1) = 0.0205315
      F0(6) = 0.123246
      F0(7) = 0.431361
      
      K = 0
      
      
c Vapor phase guess values
       FV(1) = 1.24169160741882E-10
       FV(2) = 0.233995284795399E-13
       FV(3) = 0.131268042451553E-14                    
       FV(4) = 2.19305307604124E-06
       FV(5) = 8.09339264037975E-06
       FV(6) = 0.00011381246965252        
       FV(7) = 0.430215117623187 
c LI phase guess values                                             
       FLI(1) = 0.00190574003912879  
       FLI(2) = 0.000000575729987
       FLI(3) = 0.000000000812822162
       FLI(4) = 0.000627340517412134
       FLI(5) = 3.585147944183E-08
       FLI(6) = 2.05158678764355E-05
       FLI(7) = 0.00111759832288637
c LII phase guess values                                         
       FLII(1) = 1.16455519481505E-13
       FLII(2) = 8.00232520317417E-15
       FLII(3) = 1.0508456963594849E-17              
       FLII(4) = 8.42951182479674E-10
       FLII(5) = 7.60481279177956E-04
       FLII(6) = 6.11943559835267E-06
       FLII(7) = 2.82840539378782E-05                              
        
      !Guesses for the phase density    
      RHOV = 0.19
      RHOLI = 0.871
      RHOLII = 0.881
C     RHO = MSIN(NCOMP_NCC+8,1) / 1000 !Density of feed stream in g/cm^3. [= kg/l]
    
C     Calculate total species flow rate
 4    DO i = 1, imax
      F(i) = (FV(i) + FLI(i) + FLII(i))
      END DO
      
     
C     Calculate total PHASE flow rates
C     Initialise
  
103   FVT   = 0.0D0
      FLIT  = 0.0D0
      FLIIT = 0.0D0
      DO i=1, imax
      FVT = FVT + FV(i)
      FLIT = FLIT + FLI(i)
      FLIIT = FLIIT + FLII(i)
      END DO
      WRITE(MAXWRT_MAXBUF,*) 'Total Flow rate, FVT', FVT
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Total Flow rate, FLIT', FLIT
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Total Flow rate, FLIIT', FLIIT
      CALL DMS_WRTTRM(1)
      
C     Molecular weight (kg/kmol)
      MW(1) = 885.45
      MW(2) = 621
      MW(3) = 356.55
      MW(4) = 296.49
      MW(5) = 92.094
      MW(6) = 32.042
      MW(7) = 44.009
      
C     Mole fractions and average molecular weight  
      MWVA = 0.0D0
      MWLIA = 0.0D0
      MWLIIA = 0.0D0    
      DO i=1, imax    
      XV(i) = FV(i)/FVT
      XLI(i) = FLI(i)/FLIT
      XLII(i) = FLII(i)/FLIIT
      MWVA = MWVA + XV(i)*MW(i)
      MWLIA = MWLIA + XLI(i)*MW(i)
      MWLIIA = MWLIIA + XLII(i)*MW(i)
c     ENDDO
      WRITE(MAXWRT_MAXBUF,*) 'Component', i
      CALL DMS_WRTTRM(1) 
      WRITE(MAXWRT_MAXBUF,*) 'Mole frac of Vapor phase, XV(i)', XV(i) 
      CALL DMS_WRTTRM(1)   
      WRITE(MAXWRT_MAXBUF,*) 'Mole frac- Liq-I phase, XLI(i) ', XLI(i) 
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Mole frac-Liq-II phase, XLII(i)', XLII(i)
      CALL DMS_WRTTRM(1)
       ENDDO

      WRITE(MAXWRT_MAXBUF,*) 'Avg MW of Vapor phase, MWVA', MWVA 
      CALL DMS_WRTTRM(1)  
      WRITE(MAXWRT_MAXBUF,*) 'Avg MW of Liq-I phase, MWLIA', MWLIA
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Avg MW of Liq-II  phase, MWLIIA', MWLIIA
      CALL DMS_WRTTRM(1)
      

C     Calculates phase volumes and concentration of LI and LII
      VV = TAUREA*FVT*MWVA/RHOV
      VLI = TAUREA*FLIT*MWLIA/RHOLI
      VLII = TAUREA*FLIIT*MWLIIA/RHOLII   

      WRITE(MAXWRT_MAXBUF,*) 'Vapor Phase Volume, VV', VV
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Liq-I Phase Volume, VLI', VLI
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Liq-II Phase Volume, VLII', VLII
      CALL DMS_WRTTRM(1)

      DO i=1, 6
      CLI(i) = FLI(i)*RHOLI/(FLIT*MWLIA)
      CLII(i) = FLII(i)*RHOLII/(FLIIT*MWLIIA)    
 
      WRITE(MAXWRT_MAXBUF,*) 'Component', i
      CALL DMS_WRTTRM(1) 
      WRITE(MAXWRT_MAXBUF,*) 'Liq-I phase conc, CLI(i)', CLI(i)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Liq-II phase conc, CLII(i)', CLII(i)
      CALL DMS_WRTTRM(1)
       ENDDO    
       
      !Reactions (These are for calculation of reaction rates in Liq-1 and Liq-2)
	 rxn(1) =  (k1f*(cat**ncatf)*CLI(1)*CLI(6)-                             
     +           k1r*(cat**ncatr)*CLI(2)*CLI(4))*1000

       rxn(2) = (k1f*(cat**ncatf)*CLI(2)*CLI(6)-                         
     +           k1r*(cat**ncatr)*CLI(3)*CLI(4))*1000
      
	  rxn(3) = (k1f*(cat**ncatf)*CLI(3)*CLI(6)-                          
     +           k1r*(cat**ncatr)*CLI(5)*CLI(4))*1000
        
        rxn(4) = (k1f*(cat**ncatf)*CLII(1)*CLII(6)-                      
     +           k1r*(cat**ncatr)*CLII(2)*CLII(4))*1000

	  rxn(5) = (k1f*(cat**ncatf)*CLII(2)*CLII(6)-                          
     +           k1r*(cat**ncatr)*CLII(3)*CLII(4))*1000
      
	  rxn(6) = (k1f*(cat**ncatf)*CLII(3)*CLII(6)-                           
     +           k1r*(cat**ncatr)*CLII(5)*CLII(4))*1000
        
      WRITE(MAXWRT_MAXBUF,*) 'Reaction rate, rxn(1) ', rxn(1)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Reaction rate, rxn(2)', rxn(2) 
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Reaction rate, rxn(3)', rxn(3)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Reaction rate, rxn(4)', rxn(4)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Reaction rate, rxn(5)', rxn(5)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'Reaction rate, rxn(6)', rxn(6)
      CALL DMS_WRTTRM(1)
     
C      Solves Mass balances and calculates new flow rates 
      F(1) = F0(1) - VLI*rxn(1) - VLII*rxn(4)  
      F(2) = VLI*(rxn(1) - rxn(2)) + VLII*(rxn(4) - rxn(5))
      F(3) = VLI*(rxn(2) - rxn(3)) + VLII*(rxn(5) - rxn(6))
      F(4) = VLI*(rxn(1) + rxn(2) + rxn(3)) + VLII*(rxn(4) +
     +       rxn(5) + rxn(6))
      F(5) = VLI*rxn(3) + VLII*rxn(6)
      F(6) = F0(6) - VLI*(rxn(1) + rxn(2) + rxn(3))-
     +           VLII*(rxn(4) + rxn(5) + rxn(6)) - FV(6)
      
      WRITE(MAXWRT_MAXBUF,*) 'New F values of TG',  F(1)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'New F values of DG',  F(2)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'New F values of MG',  F(3)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'New F values of FAME', F(4)
      CALL DMS_WRTTRM(1) 
      WRITE(MAXWRT_MAXBUF,*) 'New F values of GLY',  F(5)
      CALL DMS_WRTTRM(1)
      WRITE(MAXWRT_MAXBUF,*) 'New F values of MEOH',  F(6)
      CALL DMS_WRTTRM(1)
      
     
     
 
      
!     Check for convergence
       

c      DEALLOCATE(rxn)
      
c      STOP
      END SUBROUTINE TR

