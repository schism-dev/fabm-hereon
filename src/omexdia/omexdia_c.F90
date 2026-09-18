! SPDX-FileCopyrightText: 2022-2026 Helmholtz-Zentrum hereon GmbH
! SPDX-License-Identifier: Apache-2.0
! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hereon_omexdia_c --- Fortran 2003 version of OMEXDIA+P biogeochemical model
!          with additional sulfate and methane processes
!
! !INTERFACE:
   module hereon_omexdia_c
!
! !DESCRIPTION:
!
! The OMEXDIA+P model is based on the OMEXDIA model (see Soetard et al. 1996a)
! and is intended to simulate early diagenesis in the sea sediments. The major
! difference to the original OMEXDIA is an added phosphorus cycle.
! Further, the reaction ladder was extended to include methanogenesis, reoxidation
! of methane and simple ebullition.
!
! Soetaert et al. (1996a) themselves lump manganese-oxide and iron-oxide
! reduction together with sulfate reduction into a single carbon-limited
! "anoxic mineralization" pathway, precisely to avoid tracking solid MnO2/
! FeOOH inventories. This version instead resolves MnO2 and FeOOH reduction
! as their own competing pathways, positioned per the classical Froelich et
! al. (1979) redox ladder (O2 > NO3 > MnO2 > FeOOH > SO4 > CH4): each claims
! a Monod-type share of the anoxic pore space in proportion to its own solid
! oxide availability, passing the remainder down the ladder to the next
! electron acceptor -- so Mn/Fe reduction now genuinely competes for the
! same carbon supply as sulfate reduction and methanogenesis, rather than
! being folded into SulfateMin as in the lumped formulation. Both pathways
! produce reduced metal (Mn2+, Fe2+) that is lumped into the existing ODU
! currency, following Soetaert et al.'s own convention (and the "[FeS]
! approximated by ODU" comment already used below for phosphate sorption).
! MnO2/FeOOH are consumed stoichiometrically per the classical Froelich
! (1979) reactions (2 mol Mn or 4 mol Fe per mol C, reflecting the 2-electron
! vs 1-electron reduction of Mn(IV) vs Fe(III)) and resupplied at a fixed
! rate representing net delivery from the water column/bioturbated layer;
! full internal reoxidation of ODU back to solid oxide (as in Van Cappellen
! & Wang, 1996, or Boudreau, 1996 -- both cited in Soetaert et al. 1996a as
! the more detailed alternative) is not attempted here, since that would
! require tracking which fraction of the single lumped ODU pool is
! Mn/Fe/sulfide in origin.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hereon_omexdia_c
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
!
! !REVISION HISTORY:!
!  Original author(s): Richard Hofmeister & Kai Wirtz & Carsten Lemmen
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hereon_omexdia_c
!     Variable identifiers
      type (type_state_variable_id)        :: id_fdet,id_sdet,id_pdet
      type (type_state_variable_id)        :: id_no3,id_nh3,id_oxy,id_po4,id_odu
      type (type_state_variable_id)        :: id_ch4, id_ch4_gas
      type (type_state_variable_id)        :: id_mno2, id_feooh
      type (type_dependency_id)            :: id_temp, id_salinity
      type (type_diagnostic_variable_id)   :: id_denit,id_adsp
      type (type_diagnostic_variable_id)   :: id_mnmin,id_femin

!     Model parameters
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: ksO2oduox, ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox
      real(rk) :: kinSO4  ! Sulfate inhibition
      real(rk) :: ksAOM  ! Anaerobic methane oxidation
      real(rk) :: rmaxO2
      real(rk) :: ksCH4   ! methanogenesis
      real(rk) :: strip_rate ! CH4 saturation-capping relaxation rate
      real(rk) :: KMn, KFe          ! half-saturation solid-oxide availability, Mn/Fe reduction
      real(rk) :: chiMn, chiFe      ! stoichiometric oxide:carbon ratio (Froelich 1979)
      real(rk) :: SMn, SFe          ! resupply rate of MnO2/FeOOH

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hereon_omexdia_c
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the OMEXDIA+P model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the omexdia namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hereon_omexdia_c),intent(inout),target  :: self
   integer,                   intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!  Updated to FABM 1 by Nina Preußler
!  Added methane by Carsten Lemmen
!

!EOP
!-----------------------------------------------------------------------
!BOC

   ! Register parameters
   call self%get_parameter(self%rFast,'rFast','d-1','decay rate fast decay detritus ',default=0.08_rk)
   call self%get_parameter(self%rSlow,'rSlow','d-1','decay rate slow decay detritus ',default=0.0005_rk)
   call self%get_parameter(self%NCrFdet,'NCrFdet','molN molC-1','NC ratio fast decay detritus ',default=0.2_rk)
   call self%get_parameter(self%NCrSdet,'NCrSdet','molN molC-1','NC ratio slow decay detritus ',default=0.01_rk)
   call self%get_parameter(self%NH3Ads,'NH3Ads','-','Adsorption coefficient ammonium ',default=0.018_rk)
   call self%get_parameter(self%rnit,'rnit','d-1','Maximum nitrification rate ',default=300.0_rk)
   call self%get_parameter(self%ksO2nitri,'ksO2nitri','umolO2 m-3','half-saturation O2 in nitrification ',default=20.0_rk)
   call self%get_parameter(self%rODUox,'rODUox','d-1','Maximum rate oxidation of ODU ',default=20.0_rk)
   call self%get_parameter(self%ksO2oduox,'ksO2oduox','mmolO2 m-3','half-saturation O2 in oxidation of ODU ',default=10.0_rk)
   call self%get_parameter(self%ksO2oxic,'ksO2oxic','mmolO2 m-3','half-saturation O2 in oxic minerals ',default=3.0_rk)
   call self%get_parameter(self%ksNO3denit,'ksNO3denit','mmolNO3 m-3','half-saturation NO3 in denitrif ',default=1.0_rk)
   call self%get_parameter(self%kinO2denit,'kinO2denit','mmmolO2 m-3','half-saturation O2 inhib denitrif ',default=70.0_rk)
   call self%get_parameter(self%kinNO3anox,'kinNO3anox','mmolNO3 m-3','half-saturation NO3 inhib anoxic min ',default=1.0_rk)
   call self%get_parameter(self%kinO2anox,'kinO2anox','mmolO2 m-3','half-saturation O2 inhib anoxic min ',default=1.0_rk)
   call self%get_parameter(self%PAds,'PAds','-','Adsorption coefficient phosphate ',default=4.0_rk)
   call self%get_parameter(self%PAdsODU,'PAdsODU','mmol m-3','phosphate adsorbed dissolved reduced substances ',default=40.0_rk)

  ! Methane parameters configuration block
   call self%get_parameter(self%kinSO4,'kinSO4','mmol m-3','sulfate inhibition threshold for methanogenesis',default=1000.0_rk)
   call self%get_parameter(self%ksAOM,'ksAOM','m3 mmol-1 d-1','second-order anaerobic methane oxidation rate',default=0.05_rk)
   call self%get_parameter(self%rmaxO2,'rmaxO2','d-1','maximum aerobic methane oxidation rate',default=10.0_rk)
   call self%get_parameter(self%ksCH4,'ksCH4','mmol m-3','half-saturation CH4 for aerobic oxidation',default=5.0_rk)
   call self%get_parameter(self%strip_rate,'strip_rate','d-1', &
        'relaxation rate capping dissolved CH4 at ch4_sat (bubble stripping to ch4_gas); '// &
        'fast and solver-independent by design -- NOT literally 1/dt (no _DT_ macro exists '// &
        'in this FABM version)', default=100.0_rk)

   ! Manganese/iron reduction parameters. KMn/KFe are illustrative
   ! first-order estimates (not calibrated against any specific site) since
   ! neither Froelich (1979) nor Soetaert et al. (1996a) resolve solid-oxide
   ! Monod kinetics explicitly; chiMn/chiFe follow directly from the
   ! classical Froelich (1979) reaction stoichiometry.
   call self%get_parameter(self%KMn,'KMn','mmolMn m-3','half-saturation MnO2 availability, Mn reduction',default=5.0_rk)
   call self%get_parameter(self%KFe,'KFe','mmolFe m-3','half-saturation FeOOH availability, Fe reduction',default=10.0_rk)
   call self%get_parameter(self%chiMn,'chiMn','molMn molC-1','MnO2 consumed per carbon oxidized (Froelich 1979)',default=2.0_rk)
   call self%get_parameter(self%chiFe,'chiFe','molFe molC-1','FeOOH consumed per carbon oxidized (Froelich 1979)',default=4.0_rk)
   call self%get_parameter(self%SMn,'SMn','mmolMn m-3 d-1','resupply rate of reactive MnO2',default=0.5_rk)
   call self%get_parameter(self%SFe,'SFe','mmolFe m-3 d-1','resupply rate of reactive FeOOH',default=1.0_rk)

   ! Register state variables
   call self%register_state_variable(self%id_fdet, 'fdet', 'mmolC m**-3', 'fast detritus C',              4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_sdet, 'sdet', 'mmolC m**-3', 'slow detritus C',              4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_pdet, 'pdet', 'mmolP m**-3', 'detritus-P',                   4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_po4,  'po4',  'mmolP m**-3', 'dissolved phosphate',          10._rk,  minimum=0.0_rk, standard_variable=standard_variables%mole_concentration_of_phosphate)
   call self%register_state_variable(self%id_no3,  'no3',  'mmolN m**-3', 'dissolved nitrate',            20._rk,  minimum=0.0_rk, standard_variable=standard_variables%mole_concentration_of_nitrate)
   call self%register_state_variable(self%id_nh3,  'nh3',  'mmolN m**-3', 'dissolved ammonium',           40._rk,  minimum=0.0_rk, standard_variable=standard_variables%mole_concentration_of_ammonium)
   call self%register_state_variable(self%id_oxy,  'oxy',  'mmolO2 m**-3','dissolved oxygen',             100._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_odu,  'odu',  'mmol m**-3',  'dissolved reduced substances', 100._rk, minimum=0.0_rk)

   call self%register_state_variable(self%id_ch4,     'ch4',     'mmolC m**-3', 'dissolved methane',       0.0_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_ch4_gas, 'ch4_gas', 'mmolC m**-3', 'gaseous bubble methane',  0.0_rk, minimum=0.0_rk)

   call self%register_state_variable(self%id_mno2,  'mno2',  'mmolMn m**-3', 'reactive manganese oxide', 50._rk,  minimum=0.0_rk)
   call self%register_state_variable(self%id_feooh, 'feooh', 'mmolFe m**-3', 'reactive iron oxide (ferrihydrite)', 200._rk, minimum=0.0_rk)

   call self%set_variable_property(self%id_fdet,'particulate',.true.)
   call self%set_variable_property(self%id_sdet,'particulate',.true.)
   call self%set_variable_property(self%id_pdet,'particulate',.true.)
   call self%set_variable_property(self%id_po4,'particulate',.false.)
   call self%set_variable_property(self%id_no3,'particulate',.false.)
   call self%set_variable_property(self%id_nh3,'particulate',.false.)
   call self%set_variable_property(self%id_oxy,'particulate',.false.)
   call self%set_variable_property(self%id_odu,'particulate',.false.)
   call self%set_variable_property(self%id_ch4,'particulate',.false.)
   call self%set_variable_property(self%id_ch4_gas,'particulate',.false.)
   call self%set_variable_property(self%id_mno2,'particulate',.true.)
   call self%set_variable_property(self%id_feooh,'particulate',.true.)

   ! Register diagnostic variables

   call self%register_diagnostic_variable(self%id_adsp,'adsP','mmolP m**-3','phosphate adsorption', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit,'denit','mmol m**-3 d-1','denitrification rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_mnmin,'mnmin','mmol m**-3 d-1','manganese-oxide reduction rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_femin,'femin','mmol m**-3 d-1','iron-oxide reduction rate', output=output_instantaneous)

   ! Register dependencies
   
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   ! NOTE: no standard_variables%salinity in this FABM tree -- it is
   ! practical_salinity here (confirmed against base/include/standard_variables.h).
   call self%register_dependency(self%id_salinity,standard_variables%practical_salinity)

   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of OMEXDIA+P model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hereon_omexdia_c),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz & Carsten Lemmen
!
! !LOCAL VARIABLES:
   real(rk) :: fdet,sdet,oxy,odu,no3,nh3,pdet,po4
   real(rk) :: temp_celsius,temp_kelvin,f_T,E_a
   real(rk) :: radsP,Oxicminlim,Denitrilim,Rescale,rP
   real(rk),parameter :: relaxO2=0.04_rk
   real(rk),parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk),parameter :: Q10b = 1.5_rk
   real(rk) :: CprodF,CprodS,Cprod,Nprod,Pprod
   real(rk) :: Denitrific,OxicMin,Nitri,OduDepo,OduOx,pDepo
   real(rk) :: anoxic_space, SulfateMin, SulfateMinlim, Methanolim
   real(rk) :: mno2, feooh, f_mn, f_fe, Mnminlim, Feminlim, MnMin, FeMin

   ! Methane and sulfur state variable values
   real(rk) :: salinity    ! [PSU] Salinity, range: 0.0 to 35.0
   real(rk) :: so4         ! [mmol m-3] Proxy sulfate concentration, range: 0.0 to 28000.0
   real(rk) :: ch4         ! [mmolC m-3] Dissolved porewater methane concentration, range: 0.0 to 5000.0
   real(rk) :: ch4_gas     ! [mmolC m-3] Gaseous phase methane concentration, range: 0.0 to 1000.0

   ! Dimensionless kinetic and environmental scaling factors
   real(rk) :: f_so4       ! [-] Monod substrate fraction for sulfate reduction, range: 0.0 to 1.0
   real(rk) :: f_meth      ! [-] Competitive sulfate inhibition factor for methanogenesis, range: 0.0 to 1.0
   real(rk) :: f_temp_meth ! [-] Arrhenius temperature scaling factor for microbial methanogenesis, range: 0.1 to 5.0
   real(rk) :: E_a_meth    ! [K] Activation energy parameter equivalent for methanogenesis temperature curve
 
   ! Volumetric metabolic production and consumption rates
   real(rk) :: g_ch4       ! [mmolC m-3 d-1] Net volumetric production rate of methane, range: 0.0 to 100.0
   real(rk) :: r_oxic_ox   ! [mmolC m-3 d-1] Volumetric aerobic methane oxidation rate, range: 0.0 to 50.0
   real(rk) :: r_aom       ! [mmolC m-3 d-1] Volumetric anaerobic methane oxidation (AOM) rate, range: 0.0 to 20.0
   real(rk) :: r_methano   ! [mmolC m-3 d-1] Gross organic carbon mineralization via methanogenesis, range: 0.0 to 200.0
   
   real(rk) :: local_stripping
   real(rk),parameter :: q10_meth = 3.5_rk     ! Targeted Q10 scaling coefficient for methanogens
   real(rk),parameter :: ch4_sat  = 2000.0_rk   ! Critical saturation threshold in mmol/m3 (~2 mM)

   !EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp,temp_celsius)
   _GET_(self%id_fdet,fdet)
   _GET_(self%id_sdet,sdet)
   _GET_(self%id_pdet,pdet)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_odu,odu)
   _GET_(self%id_no3,no3)
   _GET_(self%id_nh3,nh3)
   _GET_(self%id_po4,po4)

   _GET_(self%id_ch4,ch4)
   _GET_(self%id_ch4_gas,ch4_gas)
   _GET_(self%id_salinity,salinity)
   _GET_(self%id_mno2,mno2)
   _GET_(self%id_feooh,feooh)

   ! Sulfate linearly depends on salinity, so this is used as a proxy
   so4 = sulfate_from_salinity(salinity)

   ! Temperature and Q10 calculations
   temp_kelvin = 273.15_rk + temp_celsius
   E_a=0.1_rk*log(Q10b)*T0*(T0+10.0_rk)
   f_T = 1.0_rk*exp(-E_a*(1.0_rk/temp_kelvin - 1.0_rk/T0))
   E_a_meth    = 0.1_rk * log(q10_meth) * T0 * (T0 + 10.0_rk)
   f_temp_meth = 1.0_rk * exp(-E_a_meth * (1.0_rk/temp_kelvin - 1.0_rk/T0))

   ! Limitation terms
   Oxicminlim = oxy / (oxy + self%ksO2oxic + relaxO2 * (nh3 + odu)) 
   Denitrilim = (1.0_rk - oxy / (oxy + self%kinO2denit)) * no3 / (no3 + self%ksNO3denit)
   
   ! The total anoxic space remaining after aerobic respiration and denitrification
   anoxic_space = (1.0_rk - oxy / (oxy + self%kinO2anox)) * (1.0_rk - no3 / (no3 + self%kinNO3anox))

   ! Sequentially partition the anoxic pore space down the Froelich (1979)
   ! redox ladder: Mn reduction claims a Monod-type share of anoxic_space
   ! first (in proportion to its own solid MnO2 availability), Fe reduction
   ! claims a share of what Mn left behind, and sulfate reduction/
   ! methanogenesis split what remains after both -- exactly generalizing
   ! the two-way f_so4/f_meth split that existed here before Mn/Fe were
   ! resolved separately.
   f_mn   = mno2 / (mno2 + self%KMn)                   ! Monod gate for Mn reduction
   f_fe   = feooh / (feooh + self%KFe)                 ! Monod gate for Fe reduction
   f_so4  = so4 / (self%kinSO4 + so4)                  ! Monod gate for sulfate reduction
   f_meth = self%kinSO4 / (self%kinSO4 + so4)          ! Sulfate inhibition gate for methanogens

   Mnminlim      = anoxic_space * f_mn
   Feminlim      = anoxic_space * (1.0_rk - f_mn) * f_fe
   SulfateMinlim = anoxic_space * (1.0_rk - f_mn) * (1.0_rk - f_fe) * f_so4
   Methanolim    = anoxic_space * (1.0_rk - f_mn) * (1.0_rk - f_fe) * f_meth * f_temp_meth ! Scaled by the methane Q10 Arrhenius factor

   ! Rescale the expanded terminal electron accepting processes (TEAPs)
   ! This guarantees total mineralization exactly equals Cprod regardless of local inhibition strengths
   Rescale = 1.0_rk / (Oxicminlim + Denitrilim + Mnminlim + Feminlim + SulfateMinlim + Methanolim)

   ! f_T is folded in here (rather than at the point of use below) so that
   ! Cprod -- and hence OxicMin+Denitrific+SulfateMin+r_methano via the
   ! Rescale partition -- equals exactly the carbon actually removed from
   ! fdet/sdet below; previously Cprod was computed without f_T while the
   ! fdet/sdet sink applied f_T separately, so the four mineralization
   ! products summed to Cprod/f_T instead of to the true carbon loss
   ! whenever f_T /= 1 (i.e. whenever temperature /= 15 degC).
   CprodF = f_T * self%rFast * fdet
   CprodS = f_T * self%rSlow * sdet
   Cprod  = CprodF + CprodS
   Nprod  = CprodF * self%NCrFdet + CprodS * self%NCrSdet


! PO4-adsorption ceases when critical capacity is reached
! [FeS] approximated by ODU
   radsP  = self%PAds * self%rSlow * (po4*max(odu,self%PAdsODU))
   rP    = self%rFast * (1.0_rk - Oxicminlim)
   Pprod  = rP * pdet

! Oxic mineralisation, denitrification, Mn/Fe reduction, sulfate and methanogenesis
   OxicMin    = Cprod * Oxicminlim    * Rescale   ! Aerobic respiration
   Denitrific = Cprod * Denitrilim    * Rescale   ! Nitrate reduction
   MnMin      = Cprod * Mnminlim      * Rescale   ! Manganese-oxide reduction (produces Mn2+, lumped into ODU)
   FeMin      = Cprod * Feminlim      * Rescale   ! Iron-oxide reduction (produces Fe2+, lumped into ODU)
   SulfateMin = Cprod * SulfateMinlim * Rescale   ! Sulfate reduction (produces native ODU)
   r_methano  = Cprod * Methanolim    * Rescale   ! Methanogenesis (produces CH4)

! Intermediate methane generation conversions and degradation kinetics sinks
   g_ch4     = 0.5_rk * r_methano
   r_oxic_ox = self%rmaxO2 * (ch4 / (self%ksCH4 + ch4)) * Oxicminlim
   r_aom     = self%ksAOM * ch4 * so4

! reoxidation and ODU deposition
   Nitri      = f_T * self%rnit   * nh3 * oxy/(oxy + self%ksO2nitri + relaxO2*(fdet + odu))
   OduOx      = f_T * self%rODUox * odu * oxy/(oxy + self%ksO2oduox + relaxO2*(nh3 + fdet))

!  pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = SulfateMin * pDepo

   ! Young Laplace (need ref) capillary threshold in sandy sediments
   ! Pc = 2gamma cos Ttheta / r_throat 
   ! with gamma suface tension of water 0.072 N m-1
   ! with theta contact angle (wet sand = 0, cos theta = 1)
   ! r_throat effective pore throat radious, can be approximated by median grain size
   ! r_throat = 0.15 * d_50 (up to 0.2 * d_50)

   ! CH4_sat = k_h(temperature, salinity) 
   ! then later bubbles connect on a gas saturation threshold, typicall 10% of pore volume.

   ! Cap dissolved CH4 at saturation and strip everything above to gaseous
   ! CH4. NOTE: was `/ _DT_` -- no such macro exists in this FABM version
   ! (models are meant to be solver/timestep-independent anyway), so this
   ! uses a fast relaxation rate (self%strip_rate, d-1) instead, matching
   ! the d-1 units of every other rate feeding into the same _CONV_UNIT_
   ! group below.
   if (ch4 > ch4_sat) then
      local_stripping = (ch4 - ch4_sat) * self%strip_rate
   else
      local_stripping = 0.0_rk
   end if

#define _CONV_UNIT_ /secs_pr_day
! reaction rates
   _ADD_SOURCE_(self%id_fdet, -CprodF _CONV_UNIT_)
   _ADD_SOURCE_(self%id_sdet, -CprodS _CONV_UNIT_)
   _ADD_SOURCE_(self%id_oxy , (-OxicMin - 2.0_rk* Nitri - OduOx) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_no3 , (-0.8_rk*Denitrific + Nitri) _CONV_UNIT_)     ! from 4/5 denitrification stoichiometry
   _ADD_SOURCE_(self%id_nh3 , (Nprod - Nitri) / (1.0_rk + self%NH3Ads) _CONV_UNIT_)
   ! NOTE: was `AnoxicMin`, a leftover from the pre-sulfate/methane 3-way
   ! partition (still in omexdia_p.F90) -- this file replaced it with the
   ! 4-way OxicMin/Denitrific/SulfateMin/r_methano split (see OduDepo
   ! above, which was correctly updated) but missed this line; AnoxicMin
   ! was never declared here, so it would not have compiled.
   ! MnMin and FeMin now also feed ODU (Mn2+, Fe2+ lumped into the same
   ! reduced-substance currency as sulfide, see header note).
   _ADD_SOURCE_(self%id_odu , (MnMin + FeMin + SulfateMin - OduOx - OduDepo) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_po4 , (f_T * Pprod - radsP) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_pdet, (radsP - f_T * Pprod) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_ch4,     (g_ch4 - r_oxic_ox - r_aom - local_stripping) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_ch4_gas, (local_stripping)_CONV_UNIT_)
   _ADD_SOURCE_(self%id_mno2,  (self%SMn - self%chiMn * MnMin) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_feooh, (self%SFe - self%chiFe * FeMin) _CONV_UNIT_)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,Denitrific)
   _SET_DIAGNOSTIC_(self%id_adsp ,radsP)
   _SET_DIAGNOSTIC_(self%id_mnmin,MnMin)
   _SET_DIAGNOSTIC_(self%id_femin,FeMin)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

  elemental function sulfate_from_salinity(salinity) result(so4)
! !DESCRIPTION:
!  Calculates the estuarine porewater sulfate proxy concentration assuming 
!  strict conservative mixing of riverine freshwater and standard seawater 
!  endmembers. The linear ratio is scaled to a standard global ocean concentration 
!  of ~28.0 mmol/m3 (mM) at 35.0 PSU salinity.
!
!  Morris, A. W., & Riley, J. P. (1966). The bromide/chlorinity and sulphate/chlorinity 
!  ratio in sea water. Deep Sea Research and Oceanographic Abstracts, 13(4), 699-705.
!  https://doi.org/10.1016/0011-7471(66)90601-2


      real(rk), intent(in) :: salinity ! [PSU] local salinity range 0 .. 35
      real(rk)             :: so4      ! [mmol m-3] calculated sulfate concentration 0 .. 28000

      so4 = (28000.0_rk / 35.0_rk) * max(0.0_rk, salinity)
   end function sulfate_from_salinity


   end module hereon_omexdia_c

