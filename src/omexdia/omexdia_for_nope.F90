#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hereon_omexdia_for_nope --- Fortran 2003 version of OMEXDIA+P biogeochemical model with adaptation to NOPE module
!
! !INTERFACE:
   module hereon_omexdia_for_nope
!
! !DESCRIPTION:
!
! The OMEXDIA+P model is based on the OMEXDIA model (see Soetaert et al. 1996a)
! and is intended to simulate early diagenesis in the sea sediments. The major
! difference to the original OMEXDIA is an added phosphorus cycle.
! In addition, this version of the model was adjusted so that it represents a water column and can be coupled
! to the NOPE (Nitrous Oxide Production and Emission) module. The following additions were made:
!  - nitrification rate as diagnostics
!  - DIN (dissolved inorganic nitrogen; nitrate + ammonium) as a diagnostic
!  - a bottom process of sediment oxygen uptake and a respective diagnostic
!
!
! SPDX-FileCopyRightText: 2021-2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileCopyRightText: 2013-2021 Helmholtz-Zentrum Geesthacht GmbH
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-FileContributor: Nina Preußler <nina.preussler@hereon.de>
! SPDX-FileContributor: Richard Hofmeister
! SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
! SPDX-LicenseRef: Apache-2.0
!
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hereon_omexdia_for_nope
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
!
! !REVISION HISTORY:!
!  Original author(s): Richard Hofmeister & Kai Wirtz
!  Adjusted to NOPE: Nina Preußler
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hereon_omexdia_for_nope
!     Variable identifiers
      type (type_state_variable_id)        :: id_fdet,id_sdet,id_pdet
      type (type_state_variable_id)        :: id_no3,id_nh3,id_oxy,id_po4,id_odu
      type (type_dependency_id)            :: id_temp
      type (type_diagnostic_variable_id)   :: id_denit,id_adsp, id_nitri, id_test1, id_test2
      type (type_diagnostic_variable_id)   :: id_din
      type (type_bottom_diagnostic_variable_id) :: id_oxy_uptake_sediments

!     Model parameters
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: ksO2oduox, ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom

   end type type_hereon_omexdia_for_nope
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
   class (type_hereon_omexdia_for_nope),intent(inout),target  :: self
   integer,                   intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!  Updated to FABM 1 + adjusted to enable coupling with NOPE: Nina Preußler
!
! !LOCAL VARIABLES:
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox


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
   call self%get_parameter(self%rnit,'rnit','d-1','Maximum nitrification rate ',default=300.0_rk)
   call self%get_parameter(self%ksO2nitri,'ksO2nitri','mmolO2 m-3','half-saturation O2 in nitrification ',default=31.3_rk) 
   call self%get_parameter(self%rODUox,'rODUox','d-1','Maximum rate oxidation of ODU ',default=20.0_rk)
   call self%get_parameter(self%ksO2oduox,'ksO2oduox','mmolO2 m-3','half-saturation O2 in oxidation of ODU ',default=10.0_rk)
   call self%get_parameter(self%ksO2oxic,'ksO2oxic','mmolO2 m-3','half-saturation O2 in oxic minerals ',default=3.0_rk)
   call self%get_parameter(self%ksNO3denit,'ksNO3denit','mmolNO3 m-3','half-saturation NO3 in denitrif ',default=36.0_rk)
   call self%get_parameter(self%kinO2denit,'kinO2denit','mmmolO2 m-3','half-saturation O2 inhib denitrif ',default=93.0_rk)
   call self%get_parameter(self%kinNO3anox,'kinNO3anox','mmolNO3 m-3','half-saturation NO3 inhib anoxic min ',default=1.0_rk)
   call self%get_parameter(self%kinO2anox,'kinO2anox','mmolO2 m-3','half-saturation O2 inhib anoxic min ',default=1.0_rk)
   call self%get_parameter(self%PAds,'PAds','-','Adsorption coefficient phosphate ',default=4.0_rk)
   call self%get_parameter(self%PAdsODU,'PAdsODU','mmol m-3','phosphate adsorbed dissolved reduced substances ',default=40.0_rk)


   ! Register state variables
   call self%register_state_variable(self%id_fdet, 'fdet', 'mmolC m**-3', 'fast detritus C',              4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_sdet, 'sdet', 'mmolC m**-3', 'slow detritus C',              4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_pdet, 'pdet', 'mmolP m**-3', 'detritus-P',                   4.e3_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_po4,  'po4',  'mmolP m**-3', 'dissolved phosphate',          10._rk,  minimum=0.0_rk, standard_variable=standard_variables%mole_concentration_of_phosphate)
   call self%register_state_variable(self%id_no3,  'no3',  'mmolN m**-3', 'dissolved nitrate',            20._rk,  minimum=0.0_rk, standard_variable=standard_variables%mole_concentration_of_nitrate)
   call self%register_state_variable(self%id_nh3,  'nh3',  'mmolN m**-3', 'dissolved ammonium',           40._rk,  minimum=0.0_rk, standard_variable=standard_variables%mole_concentration_of_ammonium)
   call self%register_state_variable(self%id_oxy,  'oxy',  'mmolO2 m**-3','dissolved oxygen',             300._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_odu,  'odu',  'mmol m**-3',  'dissolved reduced substances', 100._rk, minimum=0.0_rk)

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
   call self%register_diagnostic_variable(self%id_din,'din','mmolN m**-3','dissolved inorganic nitrogen', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oxy_uptake_sediments,'oxy_uptake_sediments','mmolO2 m-2 d-1','oxygen uptake rate of sediments', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_test1,'test1','mmol m**-3 d-1','denitrification rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_test2,'test2','mmol m**-3 d-1','denitrification rate', output=output_instantaneous)


   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   return

99 call self%fatal_error('hereon_omexdia_for_nope_initialize','Error reading namelist hereon_omexdia_for_nope.')

100 call self%fatal_error('hereon_omexdia_for_nope_initialize','Namelist hereon_omexdia_for_nope was not found.')

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
   class (type_hereon_omexdia_for_nope),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!  Updated to FABM 1 + adjusted to enable coupling with NOPE: Nina Preußler
!
! !LOCAL VARIABLES:
   real(rk) :: fdet,sdet,oxy,odu,no3,nh3,pdet,po4
   real(rk) :: temp_celsius,temp_kelvin,f_T,E_a
   real(rk) :: radsP,Oxicminlim,Denitrilim,Anoxiclim,Rescale,rP
   real(rk),parameter :: relaxO2=0.04_rk
   real(rk),parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk),parameter :: Q10b = 1.5_rk
   real(rk) :: CprodF,CprodS,Cprod,Nprod,Pprod
   real(rk) :: AnoxicMin,Denitrific,OxicMin,Nitri,OduDepo,OduOx,pDepo


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

   ! auxiliary calculations
   temp_kelvin = 273.15_rk + temp_celsius
   E_a=0.1_rk*log(Q10b)*T0*(T0+10.0_rk);        
   f_T = 1.0_rk*exp(-E_a*(1.0_rk/temp_kelvin - 1.0_rk/T0))   

   ! limitation terms
   Oxicminlim = oxy/(oxy+self%ksO2oxic+relaxO2*(nh3+odu))                
   Denitrilim = (1.0_rk-oxy/(oxy+self%kinO2denit)) * NO3/(no3+self%ksNO3denit)
   Anoxiclim  = (1.0_rk-oxy/(oxy+self%kinO2anox)) * (1.0_rk-no3/(no3+self%kinNO3anox))
   Rescale    = 1.0_rk/(Oxicminlim+Denitrilim+Anoxiclim)

   ! production of C and N from detritus decomposition
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
   AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

   ! reoxidation and ODU deposition
   Nitri      = f_T * self%rnit * nh3 * oxy/(oxy + self%ksO2nitri + relaxO2*(fdet + odu))
   OduOx      = f_T * self%rODUox * odu * oxy/(oxy + self%ksO2oduox + relaxO2*(nh3 + fdet))

   ! pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin*pDepo

#define _CONV_UNIT_ /secs_pr_day
   ! reaction rates
   _ADD_SOURCE_(self%id_fdet, (- f_T * CprodF) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_sdet, (- f_T * CprodS) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_oxy , (-OxicMin - 2.0_rk* Nitri - OduOx) _CONV_UNIT_) !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2/m**3)
   _ADD_SOURCE_(self%id_no3 , (-0.8_rk*Denitrific + Nitri) _CONV_UNIT_)     !RH 0.8-> ~104/106? 
   _ADD_SOURCE_(self%id_nh3 , (f_T * Nprod - Nitri) / (1.0_rk + self%NH3Ads) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_odu , (AnoxicMin - OduOx - OduDepo) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_po4 , (f_T * Pprod - radsP) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_pdet, (radsP - f_T * Pprod) _CONV_UNIT_)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,Denitrific)
   _SET_DIAGNOSTIC_(self%id_nitri,Nitri)
   _SET_DIAGNOSTIC_(self%id_adsp ,radsP)
   _SET_DIAGNOSTIC_(self%id_din ,no3+nh3)   
   _SET_DIAGNOSTIC_(self%id_test1,2.0_rk* Nitri _CONV_UNIT_)
   _SET_DIAGNOSTIC_(self%id_test2,OxicMin _CONV_UNIT_)
  

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do


   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)

      class (type_hereon_omexdia_for_nope),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      ! local variables
      real(rk) :: bottom_temp
      real(rk) :: oxy_consumption_sediments

      _BOTTOM_LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_temp, bottom_temp)

         ! calculate temperature-dependent oxygen uptake of sediments based on Spieckermann (2021)
         oxy_consumption_sediments = 10.07_rk + 1.73_rk * bottom_temp ! in mmolO2 m-2 d-1

         _ADD_BOTTOM_FLUX_(self%id_oxy, -oxy_consumption_sediments _CONV_UNIT_)

         _SET_BOTTOM_DIAGNOSTIC_(self%id_oxy_uptake_sediments, oxy_consumption_sediments)
         

      _BOTTOM_LOOP_END_

   end subroutine do_bottom

!EOC

   end module hereon_omexdia_for_nope

