#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hereon_omexdia_n2o_nope --- Fortran 2003 version of OMEXDIA+P biogeochemical model
!
! !INTERFACE:
   module hereon_omexdia_n2o_nope
!
! !DESCRIPTION:
!
! The OMEXDIA+P model is based on the OMEXDIA model (see Soetaert et al. 1996a)
! and is intended to simulate early diagenesis in the sea sediments. The major
! difference to the original OMEXDIA is an added phosphorus cycle and the addition of processes related to N2O.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hereon_omexdia_n2o_nope
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
!
! !REVISION HISTORY:!
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hereon_omexdia_n2o_nope
!     Variable identifiers
      type (type_state_variable_id)        :: id_fdet,id_sdet,id_pdet, id_manipulated_oxy
      type (type_state_variable_id)        :: id_no3,id_nh3,id_oxy,id_po4,id_odu
      type (type_dependency_id)            :: id_temp, id_depth
      type (type_global_dependency_id)     :: id_day_of_year
      type (type_diagnostic_variable_id)   :: id_denit,id_adsp, id_nitri, id_oduox_check, id_oxicmin_check, id_anoxicminlim_check, id_denitrilim_check, id_oxicminlim_check
      type (type_diagnostic_variable_id)   :: id_din, id_oxy_manip_check

!     Model parameters
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: ksO2oduox, ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox
      real(rk) :: oxy_increase, nh3_increase, oxy_coupling_control
      real(rk) :: oxy_manipulation, omz_oxy, omz_start, omz_end ! parameters for oxygen minimum zone (OMZ) manipulation
      !real(rk) :: no3_flux_summer, nh3_flux_summer, no3_flux_winter, nh3_flux_winter ! 
      real(rk) :: sv
      !real(rk) :: fdet_summer, fdet_winter, sdet_summer, sdet_winter

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom

   end type type_hereon_omexdia_n2o_nope
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
   class (type_hereon_omexdia_n2o_nope),intent(inout),target  :: self
   integer,                   intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!  Updated to FABM 1 + Extended to include N2O production, consumption, and emission by Nina Preußler
!
! !LOCAL VARIABLES:
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: sv
      real(rk) :: ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox
      real(rk) :: oxy_increase, nh3_increase, oxy_coupling_control
      real(rk) :: oxy_manipulation, omz_oxy, omz_start, omz_end 
      !real(rk) :: no3_flux_summer, nh3_flux_summer,no3_flux_winter, nh3_flux_winter
      !real(rk) :: fdet_summer, fdet_winter, sdet_summer, sdet_winter


!EOP
!-----------------------------------------------------------------------
!BOC

   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   ! Register parameters
   call self%get_parameter(self%rFast,'rFast','d-1','decay rate fast decay detritus ',default=0.08_rk)
   call self%get_parameter(self%rSlow,'rSlow','d-1','decay rate slow decay detritus ',default=0.0005_rk)
   call self%get_parameter(self%NCrFdet,'NCrFdet','molN molC-1','NC ratio fast decay detritus ',default=0.2_rk)
   call self%get_parameter(self%NCrSdet,'NCrSdet','molN molC-1','NC ratio slow decay detritus ',default=0.01_rk)
   call self%get_parameter(self%NH3Ads,'NH3Ads','-','Adsorption coefficient ammonium ',default=0.018_rk)
   call self%get_parameter(self%rnit,'rnit','d-1','Maximum nitrification rate ',default=300.0_rk) ! in ergom only 0.1_rk!!! default 300
   call self%get_parameter(self%ksO2nitri,'ksO2nitri','umolO2 m-3','half-saturation O2 in nitrification ',default=31.3_rk) ! scale_factor=1._rk/1000._rk)
   call self%get_parameter(self%rODUox,'rODUox','d-1','Maximum rate oxidation of ODU ',default=20.0_rk)
   call self%get_parameter(self%ksO2oduox,'ksO2oduox','mmolO2 m-3','half-saturation O2 in oxidation of ODU ',default=10.0_rk)
   call self%get_parameter(self%ksO2oxic,'ksO2oxic','mmolO2 m-3','half-saturation O2 in oxic minerals ',default=3.0_rk)
   call self%get_parameter(self%ksNO3denit,'ksNO3denit','mmolNO3 m-3','half-saturation NO3 in denitrif ',default=36.0_rk)
   call self%get_parameter(self%kinO2denit,'kinO2denit','mmmolO2 m-3','half-saturation O2 inhib denitrif ',default=93.0_rk) ! from Holzwarth diss
   call self%get_parameter(self%kinNO3anox,'kinNO3anox','mmolNO3 m-3','half-saturation NO3 inhib anoxic min ',default=1.0_rk)
   call self%get_parameter(self%kinO2anox,'kinO2anox','mmolO2 m-3','half-saturation O2 inhib anoxic min ',default=1.0_rk)
   call self%get_parameter(self%PAds,'PAds','-','Adsorption coefficient phosphate ',default=4.0_rk)
   call self%get_parameter(self%PAdsODU,'PAdsODU','mmol m-3','phosphate adsorbed dissolved reduced substances ',default=40.0_rk)
   call self%get_parameter(self%oxy_increase,'oxy_increase','mmol m-3 d-1','oxygen from horizontal transport ',default=8.0_rk)
   call self%get_parameter(self%oxy_coupling_control,'oxy_coupling_control','','controls oxygen consumption ',default=1.0_rk)
   call self%get_parameter(self%omz_oxy,'omz_oxy','mmolO2 m**-3','dissolved oxygen concentration in oxygen minimum zone ',default=60._rk)
   call self%get_parameter(self%omz_start,'omz_start','d',' starts ',default=182._rk)
   call self%get_parameter(self%omz_end,'omz_end','d','day of the year OMZ ends ',default=273._rk)
   call self%get_parameter(self%oxy_manipulation,'oxy_manipulation','mmolO2 m**-3','dissolved oxygen concentration in oxygen minimum zone ',default=60._rk)
   !call self%get_parameter(self%no3_flux_summer,'no3_flux_summer','mmolN m-3 d-1','increase of nitrate from horizontal transport',default=1.0_rk)
   !call self%get_parameter(self%nh3_flux_summer,'nh3_flux_summer','mmolN m-3 d-1','increase of ammonium from horizontal transport ',default=0.2_rk)
   !call self%get_parameter(self%no3_flux_winter,'no3_flux_winter','mmolN m-3 d-1','increase of nitrate from horizontal transport',default=2.0_rk)
   !call self%get_parameter(self%nh3_flux_winter,'nh3_flux_winter','mmolN m-3 d-1','increase of ammonium from horizontal transport ',default=0.4_rk)
   call self%get_parameter(self%sv, 'sv', 'm d-1', 'settling velocity POC', default=-0.5_rk) ! from oxypom
   ! real(rk) :: fdet_summer, fdet_winter, sdet_summer, sdet_winter
   !call self%get_parameter(self%fdet_summer, 'fdet_summer', 'mmolC m**-3', 'mean labile detritus concentration summer', default=40._rk) ! from oxypom
   !call self%get_parameter(self%fdet_winter, 'fdet_winter', 'mmolC m**-3', 'mean labile detritus concentration winter', default=70._rk) ! from oxypom
   !call self%get_parameter(self%sdet_summer, 'sdet_summer', 'mmolC m**-3', 'mean refractory detritus concentration summer', default=75.0_rk) ! from oxypom
   !call self%get_parameter(self%sdet_winter, 'sdet_winter', 'mmolC m**-3', 'mean refractory detritus concentration winter', default=-0.5_rk) ! from oxypom


   ! Register state variables
   call self%register_state_variable(self%id_fdet, 'fdet', 'mmolC m**-3', 'fast detritus C',              4.e3_rk, minimum=0.0_rk)! , vertical_movement=self%sv*d_per_s
   call self%register_state_variable(self%id_sdet, 'sdet', 'mmolC m**-3', 'slow detritus C',              4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_pdet, 'pdet', 'mmolP m**-3', 'detritus-P',                   4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_po4,  'po4',  'mmolP m**-3', 'dissolved phosphate',          10._rk,  minimum=0.0_rk, standard_variable=standard_variables%mole_concentration_of_phosphate)
   call self%register_state_variable(self%id_no3,  'no3',  'mmolN m**-3', 'dissolved nitrate',            20._rk,  minimum=0.0_rk, maximum=500._rk, standard_variable=standard_variables%mole_concentration_of_nitrate)
   call self%register_state_variable(self%id_nh3,  'nh3',  'mmolN m**-3', 'dissolved ammonium',           40._rk,  minimum=0.0_rk, maximum= 50._rk, standard_variable=standard_variables%mole_concentration_of_ammonium)
   call self%register_state_variable(self%id_oxy,  'oxy',  'mmolO2 m**-3','dissolved oxygen',             300._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_odu,  'odu',  'mmol m**-3',  'dissolved reduced substances', 100._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_manipulated_oxy,  'manipulated_oxy',  'mmolO2 m**-3','dissolved oxygen manipulated',             300._rk, minimum=0.0_rk)

   ! to check mass balance
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_fdet)
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_sdet)

   call self%set_variable_property(self%id_fdet,'particulate',.true.)
   call self%set_variable_property(self%id_sdet,'particulate',.true.)
   call self%set_variable_property(self%id_pdet,'particulate',.true.)
   call self%set_variable_property(self%id_po4,'particulate',.false.)
   call self%set_variable_property(self%id_no3,'particulate',.false.)
   call self%set_variable_property(self%id_nh3,'particulate',.false.)
   call self%set_variable_property(self%id_oxy,'particulate',.false.)
   call self%set_variable_property(self%id_odu,'particulate',.false.)


   ! Register diagnostic variables

   call self%register_diagnostic_variable(self%id_adsp,'adsP','mmolP m**-3','phosphate adsorption', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit,'denit','mmol m**-3 d-1','denitrification rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_nitri,'nitri','mmol m**-3 d-1','nitrification rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oduox_check,'oduox_check','','Oxidation of reduced substances', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oxicmin_check,'oxicmin_check','','Oxic mineralization', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oxicminlim_check,'oxicminlim_check','','Oxic mineralization limitation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_anoxicminlim_check,'anoxicminlim_check','','Anoxic mineralization limitation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denitrilim_check,'denitrilim_check','','Denitrification limitation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oxy_manip_check,'oxy_manip_check',' ','is oxy manipulation working?', output=output_instantaneous)

   ! Register dependencies
   
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_day_of_year, standard_variables%number_of_days_since_start_of_the_year)
   !call self%register_dependency(self%id_total_nitrogen, standard_variables%total_nitrogen)

   return

99 call self%fatal_error('hereon_omexdia_n2o_nope_initialize','Error reading namelist hereon_omexdia_n2o_nope.')

100 call self%fatal_error('hereon_omexdia_n2o_nope_initialize','Namelist hereon_omexdia_n2o_nope was not found.')

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
   class (type_hereon_omexdia_n2o_nope),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
   real(rk) :: fdet,sdet,oxy,odu,no3,nh3,pdet,po4, depth, oxy_manip_target, oxy_manip_diff
   real(rk) :: temp_celsius,temp_kelvin,f_T,E_a
   real(rk) :: radsP,Oxicminlim,Denitrilim,Anoxiclim,Rescale,rP
   real(rk),parameter :: relaxO2=0.04_rk  ! what is this?
   real(rk),parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk),parameter :: Q10b = 1.5_rk  ! what is this?
   real(rk) :: CprodF,CprodS,Cprod,Nprod,Pprod
   real(rk) :: AnoxicMin,Denitrific,OxicMin,Nitri,OduDepo,OduOx,pDepo
   real(rk) :: oxy_manip_check, manipulated_oxy
   real(rk) :: day_of_year
   real(rk) :: oxy_manipulation, omz_oxy, omz_start, omz_end 
 !no3_flux_summer, nh3_flux_summer, no3_flux_winter, nh3_flux_winter, 
   real(rk) :: monthly_fdet_mean(12), monthly_sdet_mean(12), expected_fdet, expected_sdet, fdet_flux, sdet_flux
   integer  :: month_nr



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
   _GET_(self%id_depth, depth)
   _GET_GLOBAL_(self%id_day_of_year,day_of_year)

   !oxy = oxy_dep

   temp_kelvin = 273.15_rk + temp_celsius
   E_a=0.1_rk*log(Q10b)*T0*(T0+10.0_rk);        ! what is this?
   f_T = 1.0_rk*exp(-E_a*(1.0_rk/temp_kelvin - 1.0_rk/T0))     ! what is this? temperature dependence sth

   Oxicminlim = oxy/(oxy+self%ksO2oxic+relaxO2*(nh3+odu))                ! limitation terms
   Denitrilim = (1.0_rk-oxy/(oxy+self%kinO2denit)) * NO3/(no3+self%ksNO3denit)
   Anoxiclim  = (1.0_rk-oxy/(oxy+self%kinO2anox)) * (1.0_rk-no3/(no3+self%kinNO3anox))
   Rescale    = 1.0_rk/(Oxicminlim+Denitrilim+Anoxiclim) ! what does this do?

   CprodF = self%rFast * fdet
   CprodS = self%rSlow * sdet
   Cprod  = CprodF + CprodS
   Nprod  = CprodF * self%NCrFdet + CprodS * self%NCrSdet


! PO4-adsorption ceases when critical capacity is reached
! [FeS] approximated by ODU
   radsP  = self%PAds * self%rSlow * (po4*max(odu,self%PAdsODU))
   rP    = self%rFast * (1.0_rk - Oxicminlim)
   Pprod  = rP * pdet

! Oxic mineralisation, denitrification, anoxic mineralisation
! then the mineralisation rates
   OxicMin    = Cprod*Oxicminlim*Rescale        ! oxic mineralisation
   Denitrific = Cprod*Denitrilim*Rescale        ! Denitrification
   !Denitrific = Denitrific - n2o_from_denitri
   AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

! reoxidation and ODU deposition
   Nitri      = f_T * self%rnit * nh3 * oxy/(oxy + self%ksO2nitri + relaxO2*(fdet + odu))
   !Nitri      = Nitri - n2o_from_nitri
   OduOx      = f_T * self%rODUox * odu * oxy/(oxy + self%ksO2oduox + relaxO2*(nh3 + fdet))

!  pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin*pDepo

! determine fdet and sdet flux based on difference to average monthly measured values (2000-2010)
   month_nr = ceiling(day_of_year / 30.5_rk)
   monthly_fdet_mean = (/ 42.48_rk, 39.75_rk, 42.09_rk, 47.21_rk, 60.25_rk, 59.64_rk, 54.53_rk, 61.27_rk, 60.86_rk, 56.71_rk, 40.65_rk, 39.75_rk /)
   monthly_sdet_mean = (/ 78.89_rk, 73.82_rk, 78.17_rk, 87.67_rk, 111.88_rk, 110.75_rk, 101.26_rk, 113.79_rk, 113.02_rk, 105.31_rk, 75.49_rk, 73.82_rk /)
   expected_fdet = monthly_fdet_mean(month_nr)
   expected_sdet = monthly_sdet_mean(month_nr)
   fdet_flux = expected_fdet - fdet
   sdet_flux = expected_sdet - sdet
   fdet_flux = 0._rk
   sdet_flux = 0._rk

#define _CONV_UNIT_ /secs_pr_day
! reaction rates
   _ADD_SOURCE_(self%id_fdet, (fdet_flux - f_T * CprodF) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_sdet, (sdet_flux - f_T * CprodS) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_oxy , (-OxicMin - 2.0_rk* Nitri - OduOx) * self%oxy_coupling_control _CONV_UNIT_) !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2/m**3)
   !_ADD_SOURCE_(self%id_oxy , self%oxy_increase ) ! added flux
   
   
   _ADD_SOURCE_(self%id_no3 , (-0.8_rk*Denitrific + Nitri) _CONV_UNIT_)     !RH 0.8-> ~104/106? whut?
   _ADD_SOURCE_(self%id_nh3 , (f_T * Nprod - Nitri) / (1.0_rk + self%NH3Ads) _CONV_UNIT_)

!   IF (day_of_year > 150._rk .AND. day_of_year < 250._rk) THEN
 !     _ADD_SOURCE_(self%id_no3, self%no3_flux_summer _CONV_UNIT_)  
  ! ELSEIF (day_of_year >  250._rk) THEN
   !   _ADD_SOURCE_(self%id_no3, -1._rk * 0.3_rk * self%no3_flux_winter _CONV_UNIT_)  
!   ELSE
 !     _ADD_SOURCE_(self%id_no3, self%no3_flux_winter _CONV_UNIT_)  
  ! END IF

!   IF (day_of_year > 200._rk .AND. day_of_year < 300._rk) THEN
 !     _ADD_SOURCE_(self%id_nh3, self%nh3_flux_summer _CONV_UNIT_)  
  ! ELSE
   !   _ADD_SOURCE_(self%id_nh3, self%nh3_flux_winter _CONV_UNIT_) 
   !END IF

   _ADD_SOURCE_(self%id_odu , (AnoxicMin - OduOx - OduDepo) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_po4 , (f_T * Pprod - radsP) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_pdet, (radsP - f_T * Pprod) _CONV_UNIT_)
   

   oxy_manip_target = 0.0_rk
   oxy_manip_diff = 0.0_rk
   !manipulate oxygen conditions
   IF (self%oxy_manipulation == 1._rk .AND. day_of_year > self%omz_start .AND. day_of_year < self%omz_end) THEN ! depth == 10._rk .AND. 
      oxy_manip_target = self%omz_oxy
      IF (oxy > oxy_manip_target .AND. depth > 10) THEN
         oxy_manip_diff = oxy - oxy_manip_target
         _ADD_SOURCE_(self%id_oxy , -oxy_manip_diff) ! added flux that manipulates oxygen
      END IF
   !ELSE
      !oxy_manip_target = 0.0_rk  ! parameter used to manipulate oxygen concentrations
   END IF


   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,Denitrific)
   _SET_DIAGNOSTIC_(self%id_nitri,Nitri)
   _SET_DIAGNOSTIC_(self%id_adsp ,radsP)
   _SET_DIAGNOSTIC_(self%id_oduox_check ,OduOx)
   _SET_DIAGNOSTIC_(self%id_oxicmin_check ,OxicMin)
   _SET_DIAGNOSTIC_(self%id_oxicminlim_check ,Oxicminlim)
   _SET_DIAGNOSTIC_(self%id_anoxicminlim_check ,Anoxiclim)
   _SET_DIAGNOSTIC_(self%id_denitrilim_check ,Denitrilim)
   _SET_DIAGNOSTIC_(self%id_din ,no3+nh3)
   _SET_DIAGNOSTIC_(self%id_oxy_manip_check , oxy_manip_diff)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

   !!! NEW: do_bottom
   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)

      class (type_hereon_omexdia_n2o_nope),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: bottom_temp
      real(rk) :: oxy_consumption_sediments
      real(rk) :: sediment_release_nh3, sediment_uptake_no3

      _BOTTOM_LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_temp, bottom_temp)

         oxy_consumption_sediments = 10.07_rk + 1.73_rk * bottom_temp ! temperature-dependent oxy consumption of sediments based on Spieckermann, 2021
         
         ! where is this data from?
         sediment_release_nh3 = 129._rk * 24._rk / 1000.0_rk ! 129 umol m-2 h-1
         sediment_uptake_no3 = 19._rk * 24._rk / 1000._rk ! -19 umol m-2 h-1

         _ADD_BOTTOM_FLUX_(self%id_oxy, -oxy_consumption_sediments _CONV_UNIT_)
         _ADD_BOTTOM_FLUX_(self%id_nh3, sediment_release_nh3 _CONV_UNIT_)
         _ADD_BOTTOM_FLUX_(self%id_no3, -sediment_uptake_no3 _CONV_UNIT_)
         

      _BOTTOM_LOOP_END_

   end subroutine do_bottom



!EOC

   end module hereon_omexdia_n2o_nope

