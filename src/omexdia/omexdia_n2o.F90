#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hereon_omexdia_n2o --- Fortran 2003 version of OMEXDIA+P biogeochemical model
!
! !INTERFACE:
   module hereon_omexdia_n2o
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
   public type_hereon_omexdia_n2o
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
!
! !REVISION HISTORY:!
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hereon_omexdia_n2o
!     Variable identifiers
      type (type_state_variable_id)        :: id_fdet,id_sdet,id_pdet
      type (type_state_variable_id)        :: id_no3,id_nh3,id_oxy,id_po4,id_odu, id_n2o, id_test, id_din, id_n2o_tang
      type (type_dependency_id)            :: id_temp, id_salinity, id_depth
      !type (type_dependency_id)            :: id_oxy_dep ! trying out oxygen as dependency
      type (type_surface_dependency_id)    :: id_wind
      type (type_diagnostic_variable_id)   :: id_denit,id_adsp, id_nitri, id_oduox_check, id_oxicmin_check, id_anoxicminlim_check, id_denitrilim_check, id_oxicminlim_check
      type (type_diagnostic_variable_id)   :: id_n2o_from_nitri, id_n2o_from_denitri
      type (type_diagnostic_variable_id)   :: id_n2o_from_nitri_tang, id_n2o_from_denitri_tang      
      type (type_surface_diagnostic_variable_id) :: id_n2o_emit, id_n2o_eq, id_k, id_Sc, id_n2o_sat

!     Model parameters
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: n2o_emission_factor_nitri, n2o_emission_factor_denitri
      real(rk) :: ksO2oduox, ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox
      real(rk) :: oxy_increase, nh3_increase, oxy_coupling_control

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: do_surface

   end type type_hereon_omexdia_n2o
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
   class (type_hereon_omexdia_n2o),intent(inout),target  :: self
   integer,                   intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!  Updated to FABM 1 by Nina Preußler
!
! !LOCAL VARIABLES:
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: n2o_emission_factor_nitri, n2o_emission_factor_denitri
      real(rk) :: ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox
      real(rk) :: oxy_increase, nh3_increase, oxy_coupling_control


!EOP
!-----------------------------------------------------------------------
!BOC

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
   call self%get_parameter(self%n2o_emission_factor_nitri,'n2o_emission_factor_nitri','','EF nitrification ',default=0.0025_rk)
   call self%get_parameter(self%n2o_emission_factor_denitri,'n2o_emission_factor_denitri','','EF denitrification ',default=0.01_rk) 

   ! Register state variables
   call self%register_state_variable(self%id_fdet, 'fdet', 'mmolC m**-3', 'fast detritus C',              4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_sdet, 'sdet', 'mmolC m**-3', 'slow detritus C',              4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_pdet, 'pdet', 'mmolP m**-3', 'detritus-P',                   4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_po4,  'po4',  'mmolP m**-3', 'dissolved phosphate',          10._rk,  minimum=0.0_rk, standard_variable=standard_variables%mole_concentration_of_phosphate)
   call self%register_state_variable(self%id_no3,  'no3',  'mmolN m**-3', 'dissolved nitrate',            20._rk,  minimum=0.0_rk, maximum=500._rk, standard_variable=standard_variables%mole_concentration_of_nitrate)
   call self%register_state_variable(self%id_nh3,  'nh3',  'mmolN m**-3', 'dissolved ammonium',           40._rk,  minimum=0.0_rk, maximum= 50._rk, standard_variable=standard_variables%mole_concentration_of_ammonium)
   call self%register_state_variable(self%id_oxy,  'oxy',  'mmolO2 m**-3','dissolved oxygen',             100._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_odu,  'odu',  'mmol m**-3',  'dissolved reduced substances', 100._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_test,  'test',  'mmolN m**-3', 'dissolved nitrous oxide',      1._rk, minimum=0.0_rk) ! will be sth in nmol, better to change to that?
   call self%register_state_variable(self%id_n2o,  'n2o',  'umolN m**-3', 'dissolved nitrous oxide',      10._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_n2o_tang,  'n2o_tang',  'umolN m**-3', 'dissolved nitrous oxide based on Tang',      10._rk, minimum=0.0_rk)
   !call self%register_state_variable(self%id_n2o_sat,  'n2o_sat',  '%', 'saturation of dissolved nitrous oxide',      0._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_din,  'din',  'mmolN m**-3', 'total dissolved inorganic nitrogen (no3+nh3)', 1._rk, minimum=0.0_rk)

   call self%set_variable_property(self%id_fdet,'particulate',.true.)
   call self%set_variable_property(self%id_sdet,'particulate',.true.)
   call self%set_variable_property(self%id_pdet,'particulate',.true.)
   call self%set_variable_property(self%id_po4,'particulate',.false.)
   call self%set_variable_property(self%id_no3,'particulate',.false.)
   call self%set_variable_property(self%id_nh3,'particulate',.false.)
   call self%set_variable_property(self%id_oxy,'particulate',.false.)
   call self%set_variable_property(self%id_odu,'particulate',.false.)
   call self%set_variable_property(self%id_test,'particulate',.false.)
   call self%set_variable_property(self%id_n2o,'particulate',.false.)
   call self%set_variable_property(self%id_n2o_tang,'particulate',.false.)
   !call self%set_variable_property(self%id_n2o_sat,'particulate',.false.)
   call self%set_variable_property(self%id_din,'particulate',.false.)


   ! Register diagnostic variables

   call self%register_diagnostic_variable(self%id_adsp,'adsP','mmolP m**-3','phosphate adsorption', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit,'denit','mmol m**-3 d-1','denitrification rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_nitri,'nitri','mmol m**-3 d-1','nitrification rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_emit,'n2o_emit','umol m**-2 d-1','rate of N2O emission', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_eq,'n2o_eq','umol m**-3','N2O concentration equilibrated with atmosphere', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_k,'k','cm h**-1','gas transfer velocity', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Sc,'Sc','','Schmidt number', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_sat,'n2o_sat','%','saturation of N2O', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oduox_check,'oduox_check','','Oxidation of reduced substances', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oxicmin_check,'oxicmin_check','','Oxic mineralization', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oxicminlim_check,'oxicminlim_check','','Oxic mineralization limitation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_anoxicminlim_check,'anoxicminlim_check','','Anoxic mineralization limitation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denitrilim_check,'denitrilim_check','','Denitrification limitation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_denitri,'n2o_from_denitri','umol m**-3 d-1','N2O production from denitrification ', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_nitri,'n2o_from_nitri','umol m**-3 d-1','N2O production from nitrification ', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_denitri_tang,'n2o_from_denitri_tang','umol m**-3 d-1','N2O production from denitrification ', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_nitri_tang,'n2o_from_nitri_tang','umol m**-3 d-1','N2O production from nitrification ', output=output_instantaneous)

   ! Register dependencies
   
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salinity,standard_variables%practical_salinity)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_wind, standard_variables%wind_speed)
   !call self%register_dependency(self%id_oxy_dep, 'oxy_dep', 'mmol m**-3', 'oxygen from dependency')

   return

99 call self%fatal_error('hereon_omexdia_n2o_initialize','Error reading namelist hereon_omexdia_n2o.')

100 call self%fatal_error('hereon_omexdia_n2o_initialize','Namelist hereon_omexdia_n2o was not found.')

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
   class (type_hereon_omexdia_n2o),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
   real(rk) :: fdet,sdet,oxy,odu,no3,nh3,pdet,po4,test, n2o, din
   real(rk) :: temp_celsius,temp_kelvin,f_T,E_a
   real(rk) :: radsP,Oxicminlim,Denitrilim,Anoxiclim,Rescale,rP
   real(rk),parameter :: relaxO2=0.04_rk  ! what is this?
   real(rk),parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk),parameter :: Q10b = 1.5_rk  ! what is this?
   real(rk) :: CprodF,CprodS,Cprod,Nprod,Pprod
   real(rk) :: AnoxicMin,Denitrific,OxicMin,Nitri,OduDepo,OduOx,pDepo
   real(rk) :: n2o_from_nitri, n2o_from_denitri
   real(rk) :: n2o_from_nitri_tang, n2o_from_denitri_tang, n2o_tang


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
   _GET_(self%id_test,test)
   _GET_(self%id_n2o,n2o)
   _GET_(self%id_n2o_tang,n2o_tang)
   _GET_(self%id_din,din)
   !_GET_(self%id_oxy_dep,oxy_dep)

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
   n2o_from_denitri = Denitrific * self%n2o_emission_factor_denitri * 0.8_rk * 1000._rk! see 0.8 below for nitrate reduction, but why? ! emission factor from Beaulieau et al. (2011)
   n2o_from_denitri_tang = 1.5_rk * exp(-0.17_rk*oxy) ! oxygen-dependent production from denitri
   n2o_from_denitri_tang = n2o_from_denitri_tang - 2.5_rk * exp(-0.47_rk*oxy) ! oxygen-dependent consumption of n2o
   !Denitrific      = Denitrific - n2o_from_denitri
   AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

! reoxidation and ODU deposition
   Nitri      = f_T * self%rnit   * nh3 * oxy/(oxy + self%ksO2nitri + relaxO2*(fdet + odu))
   n2o_from_nitri = Nitri * self%n2o_emission_factor_nitri * 1000._rk ! 1000 because Nitri in mmol m3 and n2o in umol m-3
   n2o_from_nitri_tang = ( 1000._rk * oxy / (oxy+4.3_rk) ) * 1.52_rk / (oxy + 1.59_rk) ! nitri * yield
   !Nitri      = Nitri - n2o_from_nitri
   OduOx      = f_T * self%rODUox * odu * oxy/(oxy + self%ksO2oduox + relaxO2*(nh3 + fdet))

!  pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin*pDepo

   din = no3 + nh3

#define _CONV_UNIT_ /secs_pr_day
! reaction rates
   _ADD_SOURCE_(self%id_fdet, -f_T * CprodF _CONV_UNIT_)
   _ADD_SOURCE_(self%id_sdet, -f_T * CprodS _CONV_UNIT_)
   _ADD_SOURCE_(self%id_oxy , (-OxicMin - 2.0_rk* Nitri - OduOx) * self%oxy_coupling_control _CONV_UNIT_) !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2/m**3)
   _ADD_SOURCE_(self%id_oxy , self%oxy_increase _CONV_UNIT_) ! added flux
   _ADD_SOURCE_(self%id_no3 , (-0.8_rk*Denitrific + Nitri) _CONV_UNIT_)     !RH 0.8-> ~104/106? whut?
   _ADD_SOURCE_(self%id_nh3 , (f_T * Nprod - Nitri) / (1.0_rk + self%NH3Ads) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_odu , (AnoxicMin - OduOx - OduDepo) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_po4 , (f_T * Pprod - radsP) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_pdet, (radsP - f_T * Pprod) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_test, 0.1_rk _CONV_UNIT_)
   !_ADD_SOURCE_(self%id_din, no3 + nh3 _CONV_UNIT_)  
   !_ADD_SOURCE_(self%id_n2o, din * 0.01 _CONV_UNIT_)
   _ADD_SOURCE_(self%id_n2o, (n2o_from_denitri + n2o_from_nitri) _CONV_UNIT_)   
   _ADD_SOURCE_(self%id_n2o_tang, (n2o_from_denitri_tang + n2o_from_nitri_tang) _CONV_UNIT_)   


   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,Denitrific)
   _SET_DIAGNOSTIC_(self%id_nitri,Nitri)
   _SET_DIAGNOSTIC_(self%id_adsp ,radsP)
   _SET_DIAGNOSTIC_(self%id_oduox_check ,OduOx)
   _SET_DIAGNOSTIC_(self%id_oxicmin_check ,OxicMin)
   _SET_DIAGNOSTIC_(self%id_oxicminlim_check ,Oxicminlim)
   _SET_DIAGNOSTIC_(self%id_anoxicminlim_check ,Anoxiclim)
   _SET_DIAGNOSTIC_(self%id_denitrilim_check ,Denitrilim)
   _SET_DIAGNOSTIC_(self%id_n2o_from_denitri,n2o_from_denitri)
   _SET_DIAGNOSTIC_(self%id_n2o_from_nitri,n2o_from_nitri)   
   _SET_DIAGNOSTIC_(self%id_n2o_from_denitri_tang,n2o_from_denitri_tang)
   _SET_DIAGNOSTIC_(self%id_n2o_from_nitri_tang,n2o_from_nitri_tang)    

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

   !!! NEW: do_bottom
   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)

      class (type_hereon_omexdia_n2o),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: n2o_flux_from_sedimentary_processes, id_n2o

      _BOTTOM_LOOP_BEGIN_

         ! Retrieve current (local) state variable values.

         n2o_flux_from_sedimentary_processes = 15.01_rk ! mean value from Yang

         _ADD_BOTTOM_FLUX_(self%id_n2o, n2o_flux_from_sedimentary_processes _CONV_UNIT_)

      _BOTTOM_LOOP_END_

   end subroutine do_bottom

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)

      class (type_hereon_omexdia_n2o),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: surface_test
      real(rk) :: n2o_conc, n2o_eq, k, n2o_sea_air_flux, temp, salinity, n2o_atm, depth, n2o_sat, wind_speed
      real(rk) :: SOXY
      real(rk) :: oxy, temp_abs
      real(rk) :: Sc       ! Schmidt number

      _SURFACE_LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_test,surface_test)

         _ADD_SURFACE_FLUX_(self%id_test,surface_test*(-0.00001_rk))

         !!! calculate water-air flux of n2o

         ! Retrieve current (local) state variable values.
         _GET_(self%id_n2o,n2o_conc)
         _GET_(self%id_temp, temp)
         _GET_(self%id_salinity, salinity)
         _GET_(self%id_oxy, oxy)
         _GET_(self%id_depth, depth)
         _GET_SURFACE_(self%id_wind,wind_speed)

         temp_abs = 273.15_rk + temp ! absolute temperature (in K)

         n2o_atm = 324.2_rk ! global mean atmospheric N2O concentration, which is 324.2 ppb (Yang et al., 2022)
         n2o_eq = n2o_atm * exp(-165.8806_rk + 222.8743_rk * (100.0_rk / temp_abs) + 92.0792_rk * log(temp_abs / 100.0_rk) - 1.48425_rk * (temp_abs / 100.0_rk)**2 + salinity * (-0.056235_rk + 0.0316119_rk * (temp_abs / 100.0_rk) - 0.004872 * (temp_abs/100.0_rk)**2))
         ! n2o_eq in nmol L-1 (=umol m-3)

         Sc = 2141.2_rk + (-152.56_rk) * temp + 5.8963_rk * temp ** 2.0_rk + (-0.12411_rk) * temp ** 3.0_rk + 0.0010655_rk * temp ** 4.0_rk ! based on Wanninkhof (2014)

         k = 0.39_rk * wind_speed * (Sc / 660.0_rk)**(-0.5) ! if 

         n2o_sea_air_flux = (n2o_conc - n2o_eq) * k

         n2o_sat = (n2o_conc / n2o_eq) * 100.0_rk

         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_emit, n2o_sea_air_flux)
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_eq, n2o_eq)
         _SET_SURFACE_DIAGNOSTIC_(self%id_k, k)
         _SET_SURFACE_DIAGNOSTIC_(self%id_Sc, Sc)
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_sat, n2o_sat)

         ! _SET_SURFACE_DIAGNOSTIC_(self%id_VARNAME, VALUE)


         !_ADD_SURFACE_FLUX_(self%id_n2o,n2o_sea_air_flux)
         _ADD_SURFACE_FLUX_(self%id_n2o,-n2o_sea_air_flux _CONV_UNIT_)
         _ADD_SURFACE_FLUX_(self%id_n2o_tang,-n2o_sea_air_flux _CONV_UNIT_)
         

      _SURFACE_LOOP_END_

   end subroutine do_surface   

!EOC

   end module hereon_omexdia_n2o

