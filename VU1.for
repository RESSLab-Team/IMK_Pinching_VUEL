C!**********************************************************************************************************************************************************************
C!****User Element for a 1 DOF Spring***********************************************************************************************************************************
C!**********************************************************************************************************************************************************************
      SUBROUTINE VUEL(
      !+++++++++To be Defined++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
     1 NBLOCK,
     2 RHS,AMASS,dtimeStable,
     3 SVARS,NSVARS,
     4 ENERGY,        
     5 NNODE,NDOFEL,
     6 PROPS,NPROPS,
     7 JPROPS,NJPROPS,
     8 COORDS,NCRD,
     9 U,DU,V,A,
     1 JTYPE,JELEM,
     2 TIME,PERIOD,dTIMECur,dTIMEPrev,KSTEP,KINC,LFLAGS,
     3 dMassScaleFactor,
     4 PREDEF,NPREDEF,
     5 JDLTYP,ADLMAG)
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------	   
        INCLUDE 'VABA_PARAM.INC'	 	 
      !
      !++++Operation Code++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        PARAMETER (jMassCalc   = 1,
     1 jIntForceAndDtStable = 2)
      !
      !++++Flags+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        PARAMETER (iProcedure = 1,
     1 iNlgeom 			  = 2,
     2 iOpCode 			  = 3,
     3 NFLAGS 			  = 3)
      !
      !++++Time++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        PARAMETER (iStepTIME = 1,
     1 iTotalTIME 		 = 2,
     2 NTIME 			 = 2)
      !
      !++++Procedure Flags+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        PARAMETER (jDynExplicit = 17)
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------	
      !++++Energies++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        PARAMETER (iElPd = 1,
     1 iElCd 		 = 2,
     2 iElIe 		 = 3,
     3 iElTs 		 = 4,
     4 iElDd 		 = 5,
     5 iElBv 		 = 6,
     6 iElDe 		 = 7,
     7 iElHe 		 = 8,
     8 iElKe 		 = 9,
     9 iElTh  		 = 10,
     1 iElDmd 		 = 11,
     2 iElDc 		 = 12,
     3 nElEnergy 	 = 12)
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------	
      !++++Predefined Variables++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        PARAMETER (iPredValueNew = 1,
     1 iPredValueOld			 = 2,
     2 NPRED 					 = 2)    
        PARAMETER (factorStable = 0.9d0)
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------	
        DIMENSION RHS(NBLOCK,NDOFEL), AMASS(NBLOCK,NDOFEL,NDOFEL),
     1 dtimeStable(NBLOCK),
     2 SVARS(NBLOCK,NSVARS), energy(NBLOCK,nElEnergy),
     3 PROPS(NPROPS), JPROPS(NJPROPS),
     4 JELEM(NBLOCK), TIME(NTIME), LFLAGS(NFLAGS),
     5 COORDS(NBLOCK,NNODE,NCRD), U(NBLOCK,NDOFEL),
     6 DU(NBLOCK,NDOFEL), V(NBLOCK,NDOFEL), A(NBLOCK, NDOFEL),
     7 dMassScaleFactor(NBLOCK),
     8 PREDEF(NBLOCK, NNODE, NPREDEF, NPRED), ADLMAG(NBLOCK)
        DOUBLE PRECISION  Ke,Up_pos0,Upc_pos0,Fy_pos0,Fmax_Fy_pos0,res_pos0,Uu_pos
        DOUBLE PRECISION Up_neg0,Upc_neg0,Fy_neg0,Fmax_Fy_neg0,res_neg0,Uu_neg	
        DOUBLE PRECISION lambda_S,lambda_C,lambda_A,lambda_K,c_S,c_C,c_A,c_K,D_pos,D_neg,mass
        DOUBLE PRECISION Uy_pos0,Uy_neg0,Fmax_pos0,Fmax_neg0,Umax_pos0,Umax_neg0
        DOUBLE PRECISION Kp_pos0,Kp_neg0,Kpc_pos0,Kpc_neg0,Fres_pos0,Fres_neg0,len0_2,len1_2
        DOUBLE PRECISION Ures_pos0,Ures_neg0,u0,Fp,kappa_F,kappa_D,lenX0,lenY0,lenZ0,lenX1,lenY1,lenZ1
		DOUBLE PRECISION lenX1_previous,lenY1_previous,lenZ1_previous
        DOUBLE PRECISION Ets,EtC,EtA,EtK,betaS,betaC,betaA,betaK,Umaxn,Umaxp,dEi,Epj,EpjK,Ei,Eik
		DIMENSION delta_XY(NBLOCK),delta_Y(NBLOCK),delta_ZY(NBLOCK),Vx_previous(NBLOCK),Vx(NBLOCK),Rx(NBLOCK)
		DIMENSION U_previous(NBLOCK,NDOFEL)
		DIMENSION Vy_previous(NBLOCK),Vy(NBLOCK),Ry(NBLOCK),Vz_previous(NBLOCK),Vz(NBLOCK),Rz(NBLOCK)
        DIMENSION Fy_pos_i_1(NBLOCK),Uy_pos_i_1(NBLOCK),Fy_neg_i_1(NBLOCK),Uy_neg_i_1(NBLOCK)
        DIMENSION Fmax_pos_i_1(NBLOCK),Umax_pos_i_1(NBLOCK),Fmax_neg_i_1(NBLOCK),Umax_neg_i_1(NBLOCK)
        DIMENSION Fpeak_pos_i_1(NBLOCK),Upeak_pos_i_1(NBLOCK),Fpeak_neg_i_1(NBLOCK),Upeak_neg_i_1(NBLOCK)
        DIMENSION Kp_pos_i_1(NBLOCK),Kpc_pos_i_1(NBLOCK),Kp_neg_i_1(NBLOCK),Kpc_neg_i_1(NBLOCK)
        DIMENSION KrelA_pos_i_1(NBLOCK),KrelB_pos_i_1(NBLOCK),KrelA_neg_i_1(NBLOCK),KrelB_neg_i_1(NBLOCK)
        DIMENSION Fres_pos_i_1(NBLOCK),Ures_pos_i_1(NBLOCK),Fres_neg_i_1(NBLOCK),Ures_neg_i_1(NBLOCK)
        DIMENSION Fbp_pos_i_1(NBLOCK),Ubp_pos_i_1(NBLOCK),Fbp_neg_i_1(NBLOCK),Ubp_neg_i_1(NBLOCK)
        DIMENSION Energy_Acc(NBLOCK),Energy_Diss(NBLOCK),delta(NBLOCK),us(NBLOCK)
        DIMENSION Kul_i_1(NBLOCK),Plastic_offset_pos(NBLOCK),Plastic_offset_neg(NBLOCK)
		DIMENSION cos_thetax(NBLOCK),cos_thetaz(NBLOCK),sin_thetax(NBLOCK),sin_thetaz(NBLOCK)
		DIMENSION duz_i(NBLOCK),uz_i(NBLOCK),uz_i_1(NBLOCK),duz_i_1(NBLOCK),fz_i(NBLOCK),fz_i_1(NBLOCK)
		DIMENSION uz(NBLOCK),dfx(NBLOCK),ux_i(NBLOCK),fx_i(NBLOCK)
        INTEGER j
        LOGICAL Failure_Flag(NBLOCK),Excursion_Flag(NBLOCK),Unloading_Flag(NBLOCK),FailS,FailC,FailA,FailK
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------		  
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------		  
      ! 	VUEL SUBROUTINE FOR A 1 DOF SPRING
      !	  - Pinched hysteretic behavior with basic strength, post-capping strength, accelerated reloading stiffness and unloading stiffness
      ! 	deterioration.   
      !   - Residual strength is included
      !
      !++++++Direct Integration Explicit Dynamic Analysis++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      	IF (JTYPE .EQ. 1 .AND. LFLAGS(iProcedure).EQ.jDynExplicit) THEN 
      !----------Input Properties for Positive Loading-----------------------------------------------------------------------------------------------------------------------
      		Ke 			 = PROPS(1)			! Initial elastic stiffness
      		Up_pos0 	 = PROPS(2) 		! Initial pre-capping plastic deformation in the +ve loading direction
      		Upc_pos0	 = PROPS(3) 		! Initial post-capping plastic deformation in the +ve loading direction
      		Fy_pos0      = PROPS(4)			! Initial effective force in the +ve loading direction
      		Fmax_Fy_pos0 = PROPS(5)			! Initial maximum-to-effective force ratio in the +ve loading direction
      		res_pos0     = PROPS(6)			! Residual force to effective yield ratio in the +ve loading direction
      		Uu_pos 		 = PROPS(7)			! Ultimate deformation in the +ve loading direction
      !----------Input Properties for Negative Loading------------------------------------------------------------------------------------------------------------------------						
      		Up_neg0 	 = PROPS(8)			! Initial pre-capping plastic deformation in the -ve loading direction
      		Upc_neg0	 = PROPS(9)			! Initial post-capping plastic deformation in the -ve loading direction
      		Fy_neg0      = PROPS(10)		! Initial effective force in the -ve loading direction
      		Fmax_Fy_neg0 = PROPS(11)		! Initial maximum-to-effective force ratio in the -ve loading direction
      		res_neg0     = PROPS(12)		! Residual force to effective yield ratio in the -ve loading direction
      		Uu_neg		 = PROPS(13)		! Ultimate deformation in the -ve loading direction
      !----------Input Properties for Cyclic Deterioration--------------------------------------------------------------------------------------------------------------------		
      		lambda_S 	 = PROPS(14) 		! Cyclic deterioration parameter for strength deterioration 
      		lambda_C 	 = PROPS(15) 		! Cyclic deterioration parameter for post-capping strength deterioration
      		lambda_A 	 = PROPS(16) 		! Cyclic deterioration parameter for accelerated reloading stiffness deterioration 
      		lambda_K 	 = PROPS(17) 		! Cyclic deterioration parameter for unloading stiffness deterioration 
      		c_S      	 = PROPS(18)		! Rate of strength deterioration
      		c_C      	 = PROPS(19)		! Rate of post-capping strength deterioration
      		c_A      	 = PROPS(20)		! Rate of unloading stiffness deterioration
      		c_K      	 = PROPS(21)		! Rate of accelerated reloading stiffness deterioration
			kappa_F		 = PROPS(22)		! Pinching parameter (force)
			kappa_D		 = PROPS(23)		! Pinching parameter (deformation)
      		D_pos    	 = PROPS(24)		! Rate of cyclic deterioration in the +ve loading direction
      		D_neg    	 = PROPS(25)		! Rate of cyclic deterioration in the -ve loading direction
      		mass 		 = PROPS(26)		! Assigned mass for the non-linear spring		
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------				
      		j			 = NDOFEL/2+1 		
      		Uy_pos0      = Fy_pos0/Ke  		! Deformation at effective force in the +ve loading direction
      		Uy_neg0      = Fy_neg0/Ke		! Deformation at effective force in the -ve loading direction
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------				
      		Fmax_pos0	 = Fy_pos0*Fmax_Fy_pos0		! Maximum capping force in the +ve loading direction
      		Fmax_neg0    = Fy_neg0*Fmax_Fy_neg0		! Maximum capping force in the -ve loading direction
      		Umax_pos0	 = Uy_pos0+Up_pos0			! Deformation at maximum capping force in the +ve loading direction
      		Umax_neg0	 = Uy_neg0+Up_neg0			! Deformation at maximum capping force in the -ve loading direction
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------				
      		Kp_pos0		 = (Fmax_pos0-Fy_pos0)/(Umax_pos0-Uy_pos0)		! Strain-hardening slope in the +ve loading direction
      		Kp_neg0		 = (Fmax_neg0-Fy_neg0)/(Umax_neg0-Uy_neg0)		! Strain-hardening slope in the -ve loading direction
      		Kpc_pos0	 = -Fmax_pos0/Upc_pos0							! Post-capping slope in the +ve loading direction
      		Kpc_neg0	 = -Fmax_neg0/Upc_neg0							! Post-capping slope in the -ve loading direction
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------				
      		Fres_pos0    = Fy_pos0*res_pos0								! Residual force in the +ve loading direction
      		Fres_neg0    = Fy_neg0*res_neg0								! Residual force in the -ve loading direction
      		Ures_pos0    = (Fres_pos0-Fmax_pos0)/Kpc_pos0-Umax_pos0		! Deformation at residual force in the +ve loading direction
      		Ures_neg0 	 = (Fres_neg0-Fmax_neg0)/Kpc_neg0-Umax_neg0		! Deformation at residual force in the -ve loading direction

      !----------Mass Matrix-------------------------------------------------------------------------------------------------------------------------------------------------			
      		IF (LFLAGS(iOpCode).EQ.jMassCalc) THEN
      			DO KBLOCK = 1, NBLOCK
      				AMASS(KBLOCK,1,1) = mass/2
      				AMASS(KBLOCK,2,2) = mass/2
      				AMASS(KBLOCK,3,3) = mass/2
      				AMASS(KBLOCK,4,4) = mass/2
      				AMASS(KBLOCK,5,5) = mass/2
      				AMASS(KBLOCK,6,6) = mass/2
					AMASS(KBLOCK,7,7) = mass/2
					AMASS(KBLOCK,8,8) = mass/2
					AMASS(KBLOCK,9,9) = mass/2
      			END DO
      !			   
      !----------Internal Forces and Stable Time Increment-------------------------------------------------------------------------------------------------------------------			   
      		ELSE IF ( LFLAGS(iOpCode).EQ.jIntForceAndDtStable) THEN
      			DO KBLOCK = 1, NBLOCK
      !-----------------Undamped Stable Time Increment for Translations------------------------------------------------------------------------------------------------------
      				dtimeStable(KBLOCK) = factorStable*sqrt(mass/Ke)
      !------------------Initialization--------------------------------------------------------------------------------------------------------------------------------------
      				IF (KINC.EQ.0) THEN
      					SVARS(KBLOCK,1)  = Fy_pos0
      					SVARS(KBLOCK,2)  = Uy_pos0
      					SVARS(KBLOCK,3)  = -Fy_neg0
      					SVARS(KBLOCK,4)  = -Uy_neg0
      					SVARS(KBLOCK,5)  = Fmax_pos0
      					SVARS(KBLOCK,6)  = Umax_pos0
      					SVARS(KBLOCK,7)  = -Fmax_neg0
      					SVARS(KBLOCK,8)  = -Umax_neg0
      					SVARS(KBLOCK,10) = Uy_pos0 
      					SVARS(KBLOCK,9)  = Fy_pos0 
      					SVARS(KBLOCK,12) = -Uy_neg0  
      					SVARS(KBLOCK,11) = -Fy_neg0 
      					SVARS(KBLOCK,13) = Kp_pos0
      					SVARS(KBLOCK,14) = Kpc_pos0
      					SVARS(KBLOCK,15) = Kp_neg0
      					SVARS(KBLOCK,16) = Kpc_neg0
      					SVARS(KBLOCK,17) = Ke   
      					SVARS(KBLOCK,18) = Ke   
      					SVARS(KBLOCK,19) = Ke   
      					SVARS(KBLOCK,20) = Ke   
      					SVARS(KBLOCK,21) = Ke   
      					SVARS(KBLOCK,22) = Fres_pos0
      					SVARS(KBLOCK,23) = Ures_pos0
      					SVARS(KBLOCK,24) = -Fres_neg0
      					SVARS(KBLOCK,25) = -Ures_neg0
      					SVARS(KBLOCK,26) = 0
      					SVARS(KBLOCK,27) = 0
      					SVARS(KBLOCK,28) = 0
      					SVARS(KBLOCK,29) = 0
      					SVARS(KBLOCK,30) = 0
      					SVARS(KBLOCK,31) = 0
      					SVARS(KBLOCK,32) = 0
      					SVARS(KBLOCK,33) = 0
      					SVARS(KBLOCK,34) = 0
      					SVARS(KBLOCK,35) = 0   
						SVARS(KBLOCK,36) = 0
						SVARS(KBLOCK,37) = 0
						SVARS(KBLOCK,38) = 0	
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------	 				
      				END IF          
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------	 				
      				Fy_pos_i_1(KBLOCK)         = SVARS(KBLOCK,1)
      				Uy_pos_i_1(KBLOCK)         = SVARS(KBLOCK,2)
      				Fy_neg_i_1(KBLOCK)         = SVARS(KBLOCK,3) 
      				Uy_neg_i_1(KBLOCK)         = SVARS(KBLOCK,4) 
      				Fmax_pos_i_1(KBLOCK)       = SVARS(KBLOCK,5)
      				Umax_pos_i_1(KBLOCK)   	   = SVARS(KBLOCK,6)
      				Fmax_neg_i_1(KBLOCK)       = SVARS(KBLOCK,7)
      				Umax_neg_i_1(KBLOCK)       = SVARS(KBLOCK,8)
      				Fpeak_pos_i_1(KBLOCK)      = SVARS(KBLOCK,9)
      				Upeak_pos_i_1(KBLOCK)      = SVARS(KBLOCK,10)
      				Fpeak_neg_i_1(KBLOCK)      = SVARS(KBLOCK,11)
      				Upeak_neg_i_1(KBLOCK)      = SVARS(KBLOCK,12)
      				Kp_pos_i_1(KBLOCK)         = SVARS(KBLOCK,13)
      				Kpc_pos_i_1(KBLOCK)	       = SVARS(KBLOCK,14)
      				Kp_neg_i_1(KBLOCK)	       = SVARS(KBLOCK,15)
      				Kpc_neg_i_1(KBLOCK)        = SVARS(KBLOCK,16)
      				KrelA_pos_i_1(KBLOCK)      = SVARS(KBLOCK,17)
      				KrelB_pos_i_1(KBLOCK)      = SVARS(KBLOCK,18)
      				KrelA_neg_i_1(KBLOCK)      = SVARS(KBLOCK,19)
      				KrelB_neg_i_1(KBLOCK)      = SVARS(KBLOCK,20)
      				Kul_i_1(KBLOCK) 	       = SVARS(KBLOCK,21)
      				Fres_pos_i_1(KBLOCK)       = SVARS(KBLOCK,22)
      				Ures_pos_i_1(KBLOCK)       = SVARS(KBLOCK,23)
      				Fres_neg_i_1(KBLOCK)       = SVARS(KBLOCK,24)
      				Ures_neg_i_1(KBLOCK)       = SVARS(KBLOCK,25)
      				Fbp_pos_i_1(KBLOCK)        = SVARS(KBLOCK,26)
      				Ubp_pos_i_1(KBLOCK)        = SVARS(KBLOCK,27)
      				Fbp_neg_i_1(KBLOCK)        = SVARS(KBLOCK,28)
      				Ubp_neg_i_1(KBLOCK)        = SVARS(KBLOCK,29)
      				Energy_Acc(KBLOCK)	       = SVARS(KBLOCK,30)
      				Energy_Diss(KBLOCK)        = SVARS(KBLOCK,31)
      				Failure_Flag(KBLOCK)       = SVARS(KBLOCK,32)
      				Excursion_Flag(KBLOCK)	   = SVARS(KBLOCK,33)
      				Plastic_offset_pos(KBLOCK) = SVARS(KBLOCK,34)
      				Plastic_offset_neg(KBLOCK) = SVARS(KBLOCK,35)   
                    fz_i_1(KBLOCK)			   = SVARS(KBLOCK,36) 
				    uz_i_1(KBLOCK)			   = SVARS(KBLOCK,37)
				    duz_i_1(KBLOCK)            = SVARS(KBLOCK,38) 	
      !----------------------------------------------------------------------------------------------------------------------------------------------------------------------	 		
					delta_ZY(KBLOCK) = sqrt((U(KBLOCK,j+2)-U(KBLOCK,3))**2+(U(KBLOCK,j+1)-U(KBLOCK,2))**2)
					delta_XY(KBLOCK) = sqrt((U(KBLOCK,j)-U(KBLOCK,1))**2+(U(KBLOCK,j+1)-U(KBLOCK,2))**2)

					Vz(KBLOCK) = U(KBLOCK,j+2)-U(KBLOCK,3)
					Vy(KBLOCK) = U(KBLOCK,j+1)-U(KBLOCK,2)
					Vx(KBLOCK) = U(KBLOCK,j)-U(KBLOCK,1)
					
					cos_thetaz(KBLOCK) = min((Vz(KBLOCK)/sqrt(Vz(KBLOCK)**2+Vy(KBLOCK)**2)),1.0)
					sin_thetaz(KBLOCK) = min((Vy(KBLOCK)/sqrt(Vz(KBLOCK)**2+Vy(KBLOCK)**2)),1.0)
					cos_thetax(KBLOCK) = min((Vx(KBLOCK)/sqrt(Vx(KBLOCK)**2+Vy(KBLOCK)**2)),1.0)
					sin_thetax(KBLOCK) = min((Vy(KBLOCK)/sqrt(Vx(KBLOCK)**2+Vy(KBLOCK)**2)),1.0)
					uz_i(KBLOCK) = sign(delta_ZY(KBLOCK),Vz(KBLOCK))
					ux_i(KBLOCK) = sign(delta_XY(KBLOCK),Vx(KBLOCK))
					
      				duz_i(KBLOCK) = uz_i(KBLOCK)-uz_i_1(KBLOCK)	 
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      				IF (Failure_Flag(KBLOCK).NE.1) THEN
      					! Positive Loading
      					IF (fz_i_1(KBLOCK).GE.0) THEN
      					! Early reloading before new excursion occurs
      						IF (duz_i(KBLOCK).GT.0 .AND. duz_i_1(KBLOCK).LT.0 .AND. fz_i_1(KBLOCK).GT.0) THEN
      						! Reloading stiffness KrelA if reloading is to the left of the break point
      							IF (uz_i(KBLOCK).LE.Ubp_pos_i_1(KBLOCK)) THEN
      								KrelA_pos_i_1(KBLOCK)  = (Fbp_pos_i_1(KBLOCK)-fz_i_1(KBLOCK))
     1								/(Ubp_pos_i_1(KBLOCK)-uz_i_1(KBLOCK))
      							! Reloading stiffness KrelA if reloading is to the right of the break point
      							ELSE
      								KrelB_pos_i_1(KBLOCK)  = (Fpeak_pos_i_1(KBLOCK)-fz_i_1(KBLOCK))
     1								/(Upeak_pos_i_1(KBLOCK)-uz_i_1(KBLOCK))
      							END IF
      						END IF
      						! Deformation in fz_irst reloading branch KrelA
      						IF (uz_i(KBLOCK).LE.Ubp_pos_i_1(KBLOCK) .AND. duz_i(KBLOCK).GT.0) THEN
      							dfx(KBLOCK) = KrelA_pos_i_1(KBLOCK)*duz_i(KBLOCK)
      							Umaxp = uz_i(KBLOCK);
      						! Deformation in second reloading branch KrelB
      						ELSE IF (uz_i(KBLOCK).LE.Upeak_pos_i_1(KBLOCK) .AND. duz_i(KBLOCK).GT.0) THEN
      							dfx(KBLOCK) = (KrelB_pos_i_1(KBLOCK)*uz_i(KBLOCK)+(Fpeak_pos_i_1(KBLOCK)
     1							-KrelB_pos_i_1(KBLOCK)*Upeak_pos_i_1(KBLOCK)))-fz_i_1(KBLOCK)
      							Umaxp = uz_i(KBLOCK);
      						! Deformation in post-yield branch of the backbone
      						ELSE IF (uz_i(KBLOCK).LE.Umax_pos_i_1(KBLOCK) .AND. duz_i(KBLOCK).GT.0) THEN
      							dfx(KBLOCK) = (Fy_pos_i_1(KBLOCK)+Kp_pos_i_1(KBLOCK)*(uz_i(KBLOCK)
     1							-Uy_pos_i_1(KBLOCK)))-fz_i_1(KBLOCK)
      							Umaxp = uz_i(KBLOCK);
      						! Deformation in the post-capping branch of the backbone
      						ELSE IF (uz_i(KBLOCK).GE.Umax_pos_i_1(KBLOCK) .AND. duz_i(KBLOCK).GT.0) THEN
      							! Deformation in residual branch of backbone
      							IF (uz_i(KBLOCK).GT.Ures_pos_i_1(KBLOCK)) THEN
      								dfx(KBLOCK) = Fres_pos_i_1(KBLOCK)-fz_i_1(KBLOCK)
									IF(Fres_pos_i_1(KBLOCK) == 0.0) THEN
										Failure_Flag(KBLOCK) = 1
									END IF
      							! Deformation in softening branch of the backbone
      							ELSE
      								dfx(KBLOCK) = (Fmax_pos_i_1(KBLOCK)+Kpc_pos_i_1(KBLOCK)*(uz_i(KBLOCK)
     1								-Umax_pos_i_1(KBLOCK)))-fz_i_1(KBLOCK)
      							END IF
      							Umaxp = uz_i(KBLOCK)
      						! Deformation in the unloading branch
      						ELSE
      							dfx(KBLOCK) = Kul_i_1(KBLOCK)*duz_i(KBLOCK)
      						END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      	! Negative Loading
      					ELSE
      					! Early reloading before new excursion occurs
      						IF (duz_i(KBLOCK).LT.0 .AND. duz_i_1(KBLOCK).GT.0 .AND. fz_i_1(KBLOCK).LT.0) THEN
      							KrelA_neg_i_1(KBLOCK)  = (Fbp_neg_i_1(KBLOCK)-fz_i_1(KBLOCK))
     1							/(Ubp_neg_i_1(KBLOCK)-uz_i_1(KBLOCK))
      							! Reloading stiffness KrelA if reloading is to the right of the break point
      							IF (uz_i(KBLOCK).GE.Ubp_neg_i_1(KBLOCK)) THEN
      								KrelA_neg_i_1(KBLOCK)  = (Fbp_neg_i_1(KBLOCK) -fz_i_1(KBLOCK))
     1								/(Ubp_neg_i_1(KBLOCK)-uz_i_1(KBLOCK))
      							! Reloading stiffness KrelA if reloading is to the left of the break point
      							ELSE
      								KrelB_neg_i_1(KBLOCK)  = (Fpeak_neg_i_1(KBLOCK)-fz_i_1(KBLOCK))
     1								/(Upeak_neg_i_1(KBLOCK)-uz_i_1(KBLOCK))
      							END IF
      						END IF
      						! Deformation in fz_irst reloading branch KrelA
      						IF (uz_i(KBLOCK).GE.Ubp_neg_i_1(KBLOCK) .AND. duz_i(KBLOCK).LT.0) THEN
      							dfx(KBLOCK) = KrelA_neg_i_1(KBLOCK)*duz_i(KBLOCK)
      							Umaxn = uz_i(KBLOCK)
      						! Deformation in second reloading branch KrelB
      						ELSE IF (uz_i(KBLOCK).GE.Upeak_neg_i_1(KBLOCK) .AND. duz_i(KBLOCK).LT.0) THEN
      							dfx(KBLOCK) = (KrelB_neg_i_1(KBLOCK)*uz_i(KBLOCK)+(Fpeak_neg_i_1(KBLOCK)
     1							-KrelB_neg_i_1(KBLOCK)*Upeak_neg_i_1(KBLOCK)))-fz_i_1(KBLOCK)
      							Umaxn = uz_i(KBLOCK)
      						! Deformation in post-yield branch of the backbone
      						ELSE IF (uz_i(KBLOCK).GE.Umax_neg_i_1(KBLOCK) .AND. duz_i(KBLOCK).LT.0) THEN
      							dfx(KBLOCK) = ( Fy_neg_i_1(KBLOCK)+Kp_neg_i_1(KBLOCK)*(uz_i(KBLOCK)
     1							-Uy_neg_i_1(KBLOCK)))-fz_i_1(KBLOCK)
      							Umaxn = uz_i(KBLOCK)
      						! Deformation in the post-capping branch of the backbone
      						ELSE IF (uz_i(KBLOCK).LE.Umax_neg_i_1(KBLOCK) .AND. duz_i(KBLOCK).LT.0) THEN
      							! Deformation in residual branch of backbone
      							IF (uz_i(KBLOCK).LT.Ures_neg_i_1(KBLOCK)) THEN
      								dfx(KBLOCK) = Fres_neg_i_1(KBLOCK)-fz_i_1(KBLOCK)
									IF(Fres_neg_i_1(KBLOCK) == 0.0) THEN
										Failure_Flag(KBLOCK) = 1
									END IF
      							! Deformation in the softening branch of the backbone
      							ELSE
      								dfx(KBLOCK) = Fmax_neg_i_1(KBLOCK)+Kpc_neg_i_1(KBLOCK)*(uz_i(KBLOCK)
     1								-Umax_neg_i_1(KBLOCK))-fz_i_1(KBLOCK)
      							END IF
      							Umaxn = uz_i(KBLOCK)
      						ELSE
      							! Deformation in the unloading branch
      							dfx(KBLOCK) = Kul_i_1(KBLOCK) *duz_i(KBLOCK)
      						END IF
      					END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Deterioration Parameters
      					! Internal energy increment
      					dEi = 0.5*(dfx(KBLOCK)+2*fz_i_1(KBLOCK))*(uz_i(KBLOCK)-uz_i_1(KBLOCK))
						EtS = lambda_S*(Fy_pos0+Fy_neg0)/2*(Up_pos0+Up_neg0)/2
						EtC	= lambda_C*(Fy_pos0+Fy_neg0)/2*(Up_pos0+Up_neg0)/2
						EtA	= lambda_A*(Fy_pos0+Fy_neg0)/2*(Up_pos0+Up_neg0)/2
						EtK	= lambda_K*(Fy_pos0+Fy_neg0)/2*(Up_pos0+Up_neg0)/2
      					! Positive excursion flag
      					IF (fz_i_1(KBLOCK)+dfx(KBLOCK).GE.0 .AND. fz_i_1(KBLOCK).LT.0) THEN
      						Excursion_Flag(KBLOCK) = 1
      					! Negative excursion flag
      					ELSE IF (fz_i_1(KBLOCK)+dfx(KBLOCK).LE.0 .AND. fz_i_1(KBLOCK).GT.0) THEN
      						Excursion_Flag(KBLOCK) = 1
      					ELSE
      						Excursion_Flag(KBLOCK) = 0
      					END IF
      					! Update beta parameters at new excursion
      					IF (Excursion_Flag(KBLOCK).EQ.1) THEN
      						! Total energy dissipated in all previous excursions
      						Epj = Energy_Acc(KBLOCK)+dEi
      						! Energy dissipated in previous excursion
      						Ei = Epj-Energy_Diss(KBLOCK)
      						betaS = (Ei/(EtS-Epj))**c_S;
      						betaC = (Ei/(EtC-Epj))**c_C;
      						betaA = (Ei/(EtA-Epj))**c_A;
      					ELSE
      						! Total energy dissipated in all previous excursions
      						Epj = Energy_Diss(KBLOCK)
      						betaS = 0
      						betaC = 0
      						betaA = 0
      					END IF
      					! Onset of unloading
      					Unloading_Flag(KBLOCK) = duz_i(KBLOCK)*duz_i_1(KBLOCK).LT.0 .AND.
     1					(uz_i_1(KBLOCK).GE.Upeak_pos_i_1(KBLOCK) .OR. uz_i_1(KBLOCK).LE.Upeak_neg_i_1(KBLOCK))
      					IF (Unloading_Flag(KBLOCK).EQ.1) THEN
      						! Total energy dissipated until point of unloading
      						EpjK = dEi+Energy_Acc(KBLOCK)-0.5*((fz_i_1(KBLOCK)+dfx(KBLOCK))**2)/Kul_i_1(KBLOCK) 
      						! Energy dissipated in current excursion until point of unloading
      						EiK = EpjK-Energy_Diss(KBLOCK)
      						betaK = (EiK/(EtK-EpjK))**c_K
      					ELSE
      						betaK = 0
      					END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Target Peak Deformation
      					New_Peak_Pos_Flag = 0
      					New_Peak_Neg_Flag = 0
      					! Update target peak deformation for positive loading
      					IF (Umaxp.GE.Upeak_pos_i_1(KBLOCK)) THEN
      						New_Peak_Pos_Flag = 1
      						Upeak_pos_i_1(KBLOCK) = Umaxp
      						Fpeak_pos_i_1(KBLOCK) = fz_i_1(KBLOCK)+dfx(KBLOCK)
      						! Plastic offset for positive loading
      						Plastic_offset_pos(KBLOCK) = Upeak_pos_i_1(KBLOCK)-Fpeak_pos_i_1(KBLOCK)
     1						/Kul_i_1(KBLOCK)
      					END IF
      					! Update target peak deformation for negative loading
      					IF (Umaxn.LE.Upeak_neg_i_1(KBLOCK)) THEN
      						New_Peak_Neg_Flag = 1
      						Upeak_neg_i_1(KBLOCK) = Umaxn
      						Fpeak_neg_i_1(KBLOCK)  = fz_i_1(KBLOCK)+dfx(KBLOCK)
      						! Plastic offset for negative loading
      						Plastic_offset_neg(KBLOCK)  = Upeak_neg_i_1(KBLOCK)-Fpeak_neg_i_1(KBLOCK)
     1						/Kul_i_1(KBLOCK)
      					END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Update Positive Backbone and Target Peak Point
      					IF (Excursion_Flag(KBLOCK).EQ.1) THEN
      						! Positive loading backbone
      						IF (fz_i_1(KBLOCK).LT.0) THEN
      							! Basic strength deterioration: Yield point
      							Uy_pos_i_1(KBLOCK)  = max(Uy_pos_i_1(KBLOCK)-Fy_pos_i_1(KBLOCK)*betaS*
     1 							D_pos/Ke,Fres_pos_i_1(KBLOCK)/Ke)
      							Fy_pos_i_1(KBLOCK)  = max(Fy_pos_i_1(KBLOCK) *(1-betaS* D_pos),
     2							Fres_pos_i_1(KBLOCK))
      							! Basic strength deterioration: Post-yield Stiffness
      							IF (Fy_pos_i_1(KBLOCK).NE.Fres_pos_i_1(KBLOCK)) THEN
      								Kp_pos_i_1(KBLOCK) = Kp_pos_i_1(KBLOCK) *(1-betaS*D_pos)
      							ELSE
      								Kp_pos_i_1(KBLOCK)  = 0
      							END IF
      							! Basic strength deterioration: Capping Point
      							sPCsp = (Fy_pos_i_1(KBLOCK) -Uy_pos_i_1(KBLOCK) *Kp_pos_i_1(KBLOCK)
     1							-Fmax_pos_i_1(KBLOCK)+Kpc_pos_i_1(KBLOCK)*Umax_pos_i_1(KBLOCK))
     2							/(Kpc_pos_i_1(KBLOCK)-Kp_pos_i_1(KBLOCK) )
      							Fmax_pos_i_1(KBLOCK)  = Fmax_pos_i_1(KBLOCK) +(sPCsp-Umax_pos_i_1(KBLOCK))
     1							*Kpc_pos_i_1(KBLOCK)
      							Umax_pos_i_1(KBLOCK)  = sPCsp
      							! Post-capping strength deterioration: Capping point
      							sPCpcp = max(Umax_pos_i_1(KBLOCK) +betaC*D_pos*(Fmax_pos_i_1(KBLOCK)
     1							-Kpc_pos_i_1(KBLOCK)*Umax_pos_i_1(KBLOCK) )/(Kpc_pos_i_1(KBLOCK)
     2						    -Kp_pos_i_1(KBLOCK)),Uy_pos_i_1(KBLOCK) )
      							Fmax_pos_i_1(KBLOCK)  = Fmax_pos_i_1(KBLOCK) +(sPCpcp-Umax_pos_i_1(KBLOCK))
     1							*Kp_pos_i_1(KBLOCK)
      							Umax_pos_i_1(KBLOCK)  = sPCpcp
     							! Accelerated reloading stiffness deterioration: Target peak deformation point
      							Upeak_pos_i_1(KBLOCK)  = (1+betaA*D_pos)*Upeak_pos_i_1(KBLOCK)
      							! Target peak deformation in reloading branch of the updated backbone
      							IF (Upeak_pos_i_1(KBLOCK).LE.Uy_pos_i_1(KBLOCK)) THEN
      								Fpeak_pos_i_1(KBLOCK) = Ke*Upeak_pos_i_1(KBLOCK)
      							! Target peak deformation in post-yield branch of the updated backbone
      							ELSE IF (Upeak_pos_i_1(KBLOCK).LE.Umax_pos_i_1(KBLOCK)) THEN
									Fpeak_pos_i_1(KBLOCK) = Kp_pos_i_1(KBLOCK) *(Upeak_pos_i_1(KBLOCK)
     1								-Uy_pos_i_1(KBLOCK) )+Fy_pos_i_1(KBLOCK)
      							! Target peak deformation in post-capping branch of the updated backbone
      							ELSE
      								Fpeak_pos_i_1(KBLOCK) = max(Kpc_pos_i_1(KBLOCK)*(Upeak_pos_i_1(KBLOCK)
     1								-Umax_pos_i_1(KBLOCK) )+Fmax_pos_i_1(KBLOCK) ,Fres_pos_i_1(KBLOCK))
      							END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      						! Update Negative Backbone and Target Peak Point
      						ELSE
      							! Basic strength deterioration: Yield point
      							Uy_neg_i_1(KBLOCK)  = min(Uy_neg_i_1(KBLOCK) -Fy_neg_i_1(KBLOCK)
     1 							*betaS* D_neg/Ke,Fres_neg_i_1(KBLOCK)/Ke)
      							Fy_neg_i_1(KBLOCK)  = min(Fy_neg_i_1(KBLOCK) *(1-betaS* D_neg),
     1							Fres_neg_i_1(KBLOCK))
      							! Basic strength deterioration: Post-yield stiffness
      							IF (Fy_neg_i_1(KBLOCK).NE.Fres_neg_i_1(KBLOCK)) THEN
      								Kp_neg_i_1(KBLOCK) = Kp_neg_i_1(KBLOCK) *(1-betaS* D_neg)
      							ELSE
      								Kp_neg_i_1(KBLOCK)  = 0
      							END IF
      							! Basic strength deterioration: Capping point
      							sPCsn = (Fy_neg_i_1(KBLOCK) -Uy_neg_i_1(KBLOCK) *Kp_neg_i_1(KBLOCK)
     1							-Fmax_neg_i_1(KBLOCK)+Kpc_neg_i_1(KBLOCK)*Umax_neg_i_1(KBLOCK))
     2							/(Kpc_neg_i_1(KBLOCK)-Kp_neg_i_1(KBLOCK))
      							Fmax_neg_i_1(KBLOCK) = Fmax_neg_i_1(KBLOCK)+(sPCsn
     1							-Umax_neg_i_1(KBLOCK) )*Kpc_neg_i_1(KBLOCK)
      							Umax_neg_i_1(KBLOCK) = sPCsn
      							! Post-capping strength deterioration: Capping point
      							sPCpcn = min(Umax_neg_i_1(KBLOCK) +betaC* D_neg*(Fmax_neg_i_1(KBLOCK)
     1							-Kpc_neg_i_1(KBLOCK)*Umax_neg_i_1(KBLOCK))
     2							/(Kpc_neg_i_1(KBLOCK)-Kp_neg_i_1(KBLOCK)),Uy_neg_i_1(KBLOCK))
      							Fmax_neg_i_1(KBLOCK) = Fmax_neg_i_1(KBLOCK) +(sPCpcn-Umax_neg_i_1(KBLOCK))
     1							*Kp_neg_i_1(KBLOCK)
      							Umax_neg_i_1(KBLOCK) = sPCpcn
      							! Accelerated reloading stiffness deterioration: Target peak deformation point
      							Upeak_neg_i_1(KBLOCK)  = (1+betaA*D_neg)*Upeak_neg_i_1(KBLOCK) 
      							! Target peak deformation in reloading branch of the updated backbone
      							IF (Upeak_neg_i_1(KBLOCK).GE.Uy_neg_i_1(KBLOCK)) THEN
      								Fpeak_neg_i_1(KBLOCK) = Ke*Upeak_neg_i_1(KBLOCK)
      							! Target peak deformation in post-yield branch of the updated backbone
      							ELSE IF (Upeak_neg_i_1(KBLOCK).GE.Umax_neg_i_1(KBLOCK)) THEN
      								Fpeak_neg_i_1(KBLOCK)  = Kp_neg_i_1(KBLOCK) *(Upeak_neg_i_1(KBLOCK)
     1								-Uy_neg_i_1(KBLOCK))+Fy_neg_i_1(KBLOCK)
      							! Target peak deformation in post-capping branch of the updated backbone
      							ELSE
      								Fpeak_neg_i_1(KBLOCK) = min(Kpc_neg_i_1(KBLOCK)*(Upeak_neg_i_1(KBLOCK)
     1								-Umax_neg_i_1(KBLOCK))+Fmax_neg_i_1(KBLOCK),Fres_neg_i_1(KBLOCK))
      							END IF
      						END IF
      					END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Update Krel Based on New Peak Targets at New Positive Excursion
      					IF (fz_i_1(KBLOCK)+dfx(KBLOCK).GE.0 .AND. fz_i_1(KBLOCK).LT.0) THEN
      						Fp = kappa_F*Fpeak_pos_i_1(KBLOCK)
      						! Deformation at reloading
      						u0 = uz_i_1(KBLOCK)-(fz_i_1(KBLOCK))/Kul_i_1(KBLOCK)
      						IF (u0.LT.0) THEN
      							! Deformation at break point
      							Ubp_pos_i_1(KBLOCK) = (1-kappa_U)*Plastic_offset_pos(KBLOCK)
      							! Force at break point
      							Fbp_pos_i_1(KBLOCK) = Fp*(Ubp_pos_i_1(KBLOCK)-u0)/(Upeak_pos_i_1(KBLOCK)-u0)
      						END IF
      						! Reloading is to the left of the break point
      						IF (u0.LE.Ubp_pos_i_1(KBLOCK)) THEN
      							! Reloading stiffness KrelA after new excursion
      							KrelA_pos_i_1(KBLOCK) = (Fbp_pos_i_1(KBLOCK))/(Ubp_pos_i_1(KBLOCK)-u0)
      							! Reloading stiffness KrelB after new excursion
      							KrelB_pos_i_1(KBLOCK) = (Fpeak_pos_i_1(KBLOCK) -Fbp_pos_i_1(KBLOCK))
     1							/(Upeak_pos_i_1(KBLOCK) -Ubp_pos_i_1(KBLOCK))
      							dfx(KBLOCK) = ((uz_i(KBLOCK) - u0)*KrelA_pos_i_1(KBLOCK)) - fz_i_1(KBLOCK)
      						! Reloading is to the right of the break point
      						ELSE
      							! Reloading stiffness after new excursion
      							KrelB_pos_i_1(KBLOCK) = (Fpeak_pos_i_1(KBLOCK))/(Upeak_pos_i_1(KBLOCK)-u0)
      							dfx(KBLOCK) = ((uz_i(KBLOCK)-u0)*KrelB_pos_i_1(KBLOCK))-fz_i_1(KBLOCK)
      						END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Update Krel Based on New Peak Targets at New Negative Excursion
      					ELSE IF (fz_i_1(KBLOCK)+dfx(KBLOCK).LT.0 .AND. fz_i_1(KBLOCK).GT.0) THEN
      						Fp = kappa_F*Fpeak_neg_i_1(KBLOCK) 
      						! Deformation at reloading
      						u0 = uz_i_1(KBLOCK)-(fz_i_1(KBLOCK))/Kul_i_1(KBLOCK)
      						IF (u0.GT.0) THEN
      							! Deformation at break point
      							Ubp_neg_i_1(KBLOCK) = (1-kappa_U)*Plastic_offset_neg(KBLOCK)
      							! Force at break point
      							Fbp_neg_i_1(KBLOCK) = Fp*(Ubp_neg_i_1(KBLOCK)-u0)/(Upeak_neg_i_1(KBLOCK)-u0)
      						END IF
      						! Reloading is to the right of the break point
      						IF (u0.GE.Ubp_neg_i_1(KBLOCK)) THEN
      							! Reloading stiffness KrelA after new excursion
      							KrelA_neg_i_1(KBLOCK) = (Fbp_neg_i_1(KBLOCK) )/(Ubp_neg_i_1(KBLOCK)-u0)
      							! Reloading stiffness KrelB after new excursion
      							KrelB_neg_i_1(KBLOCK) = (Fpeak_neg_i_1(KBLOCK)-Fbp_neg_i_1(KBLOCK))
     1							/(Upeak_neg_i_1(KBLOCK)-Ubp_neg_i_1(KBLOCK))
      							dfx(KBLOCK) = ((uz_i(KBLOCK) - u0)*KrelA_neg_i_1(KBLOCK))-fz_i_1(KBLOCK)
      						! Reloading is to the left of the break point
      						ELSE
      							! Reloading stiffness after new excursion
      							KrelB_neg_i_1(KBLOCK) = (Fpeak_neg_i_1(KBLOCK) )/(Upeak_neg_i_1(KBLOCK)-u0)
      							dfx(KBLOCK) = ((uz_i(KBLOCK)-u0)*KrelB_neg_i_1(KBLOCK))-fz_i_1(KBLOCK)
      						END IF
      					END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Update Unloading Stiffness
      					IF (Unloading_Flag(KBLOCK).EQ.1) THEN
      						Kul_i_1(KBLOCK)  = (1-betaK)*Kul_i_1(KBLOCK)
      						dfx(KBLOCK) = Kul_i_1(KBLOCK)*duz_i(KBLOCK)
      					END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Update Deformation at Residual Points
      					! Deformation at residual onset for positive backbone
      					Ures_pos_i_1(KBLOCK) = (Fres_pos_i_1(KBLOCK)-Fmax_pos_i_1(KBLOCK)+Kpc_pos_i_1(KBLOCK)
     1					*Umax_pos_i_1(KBLOCK))/Kpc_pos_i_1(KBLOCK)
      					! Deformation at residual onset for negative backbone
      					Ures_neg_i_1(KBLOCK) = (Fres_neg_i_1(KBLOCK)-Fmax_neg_i_1(KBLOCK)+Kpc_neg_i_1(KBLOCK)
     1					*Umax_neg_i_1(KBLOCK))/Kpc_neg_i_1(KBLOCK)
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Force
      					fz_i(KBLOCK) = fz_i_1(KBLOCK)+dfx(KBLOCK)
						fx_i(KBLOCK) = ux_i(KBLOCK)*Ke
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Refz_ine Peaks
      					IF (New_Peak_Pos_Flag .EQ. 1) THEN
      						Fpeak_pos_i_1(KBLOCK) = fz_i(KBLOCK)
      					ELSE IF (New_Peak_Neg_Flag.EQ.1) THEN
      						Fpeak_neg_i_1(KBLOCK) = fz_i(KBLOCK)
      					END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      					! Failure
      					! Failure criteria (Tolerance = 1%)
      					FailS = (betaS.LT.-0.01 .OR. betaS.GT.1.01)
      					FailC = (betaC.LT.-0.01 .OR. betaC.GT.1.01)
      					FailA = (betaA.LT.-0.01 .OR. betaA.GT.1.01)
      					FailK = (betaK.LT.-0.01 .OR. betaK.GT.1.01)
      					IF (FailS .AND. FailC .AND. FailA .OR. FailK) THEN
      						Failure_Flag(KBLOCK) = 1
      					ELSE IF (uz_i(KBLOCK).GE.Uu_pos) THEN
      						Failure_Flag(KBLOCK) = 1
      					ELSE IF (uz_i(KBLOCK).LE.-Uu_neg) THEN
      						Failure_Flag(KBLOCK) = 1
						ELSE IF (ux_i(KBLOCK).GE.Uu_pos) THEN
							Failure_Flag(KBLOCK) = 1
						ELSE IF (ux_i(KBLOCK).LE.-Uu_neg) THEN
							Failure_Flag(KBLOCK) = 1
      					END IF
      					! Force at failure
					ELSE
      					fz_i(KBLOCK) = 0
						fx_i(KBLOCK) = 0
      					dEi = 0
      					Epj = Energy_Acc(KBLOCK)+dEi
      				END IF
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      				! Energies
      				! Total internal energy accumulated until current increment
      				Energy_Acc(KBLOCK) = Energy_Acc(KBLOCK)+dEi
      				! Total energy dissipated in all previous excursions
      				Energy_Diss(KBLOCK) = Epj
      ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      				! Update Variables
      				duz_i_1(KBLOCK) = duz_i(KBLOCK)
      				Umaxp = 0
      				Umaxn = 0
      				! Stiffness
      				IF (fz_i(KBLOCK) .EQ. fz_i_1(KBLOCK)) THEN
      					Ki = 10**-6
      				ELSE IF (duz_i(KBLOCK).EQ.0) THEN
      					Ki = Ke
      				ELSE
      					Ki = (fz_i(KBLOCK)-fz_i_1(KBLOCK))/(duz_i(KBLOCK))
      				END IF
C!------------------Energies----------------------------------------------------------------------------------------------------------------------------------------					
					ENERGY(KBLOCK, iElIe) = Energy_Acc(KBLOCK)
					IF (Energy_Acc(KBLOCK).LT.SVARS(KBLOCK,30)) THEN
						ENERGY(KBLOCK, iElPd) = ENERGY(KBLOCK, iElPd)+SVARS(KBLOCK,30)-Energy_Acc(KBLOCK)	
					END IF						
					
      !######################################################################################################################################################################
      					SVARS(KBLOCK,1)  = Fy_pos_i_1(KBLOCK) 
      					SVARS(KBLOCK,2)  = Uy_pos_i_1(KBLOCK)    
      					SVARS(KBLOCK,3)  = Fy_neg_i_1(KBLOCK)
      					SVARS(KBLOCK,4)  = Uy_neg_i_1(KBLOCK)       
      					SVARS(KBLOCK,5)  = Fmax_pos_i_1(KBLOCK)   
      					SVARS(KBLOCK,6)  = Umax_pos_i_1(KBLOCK) 
      					SVARS(KBLOCK,7)  = Fmax_neg_i_1(KBLOCK)  
      					SVARS(KBLOCK,8)  = Umax_neg_i_1(KBLOCK)  
      					SVARS(KBLOCK,9)  = Fpeak_pos_i_1(KBLOCK)   
      					SVARS(KBLOCK,10) = Upeak_pos_i_1(KBLOCK)   
      					SVARS(KBLOCK,11) = Fpeak_neg_i_1(KBLOCK)    
      					SVARS(KBLOCK,12) = Upeak_neg_i_1(KBLOCK)    
      					SVARS(KBLOCK,13) = Kp_pos_i_1(KBLOCK)   
      					SVARS(KBLOCK,14) = Kpc_pos_i_1(KBLOCK)	    
      					SVARS(KBLOCK,15) = Kp_neg_i_1(KBLOCK)	 
      					SVARS(KBLOCK,16) = Kpc_neg_i_1(KBLOCK)     
      					SVARS(KBLOCK,17) = KrelA_pos_i_1(KBLOCK)  
      					SVARS(KBLOCK,18) = KrelB_pos_i_1(KBLOCK)   
      					SVARS(KBLOCK,19) = KrelA_neg_i_1(KBLOCK)  
      					SVARS(KBLOCK,20) = KrelB_neg_i_1(KBLOCK)  
      					SVARS(KBLOCK,21) = Kul_i_1(KBLOCK) 	
      					SVARS(KBLOCK,22) = Fres_pos_i_1(KBLOCK)    
      					SVARS(KBLOCK,23) = Ures_pos_i_1(KBLOCK)     
      					SVARS(KBLOCK,24) = Fres_neg_i_1(KBLOCK)   
      					SVARS(KBLOCK,25) = Ures_neg_i_1(KBLOCK)    
      					SVARS(KBLOCK,26) = Fbp_pos_i_1(KBLOCK)  
      					SVARS(KBLOCK,27) = Ubp_pos_i_1(KBLOCK) 
      					SVARS(KBLOCK,28) = Fbp_neg_i_1(KBLOCK) 
      					SVARS(KBLOCK,29) = Ubp_neg_i_1(KBLOCK)   
      					SVARS(KBLOCK,30) = Energy_Acc(KBLOCK)	 
      					SVARS(KBLOCK,31) = Energy_Diss(KBLOCK)  
      					SVARS(KBLOCK,32) = Failure_Flag(KBLOCK)    
      					SVARS(KBLOCK,33) = Excursion_Flag(KBLOCK)	
      					SVARS(KBLOCK,34) = Plastic_offset_pos(KBLOCK) 
      					SVARS(KBLOCK,35) = Plastic_offset_neg(KBLOCK) 
						SVARS(KBLOCK,36) = fz_i(KBLOCK)
						SVARS(KBLOCK,37) = uz_i(KBLOCK)
						SVARS(KBLOCK,38) = duz_i(KBLOCK)	 
      					RHS(KBLOCK,1)    = -abs(fx_i(KBLOCK))*cos_thetax(KBLOCK)
						RHS(KBLOCK,j) 	 = -RHS(KBLOCK,1)
						RHS(KBLOCK,3)    = -abs(fz_i(KBLOCK))*cos_thetaz(KBLOCK)
						RHS(KBLOCK,j+2)  = -RHS(KBLOCK,3)
						RHS(KBLOCK,2)    = -abs(fx_i(KBLOCK))*sin_thetax(KBLOCK)-abs(fz_i(KBLOCK))*sin_thetaz(KBLOCK)
						RHS(KBLOCK,j+1)  = -RHS(KBLOCK,2)
      			END DO	
      		END IF
      	END IF
       END        
      