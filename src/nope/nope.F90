#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hereon_nope --- Fortran 2003 model of aquatic nitrous oxide production and emission
!
! !INTERFACE:
   module hereon_nope
!
! !DESCRIPTION:
!
! This module can be coupled with estuarine and coastal biogeochemical models to simulate 
! nitrous oxide (N2O) production in the water column and its emission to the atmosphere. 
!
! It includes multiple alternatives for calculating of N2O production rates based on different
! yields (% of NH3/NO3 converted in nitrification/denitrification that results in N2O). 
! Furthermore, there are multiple options for calculating the N2O sea-air-flux (= emission rate):
! three types of wind speed (hourly, seasonal, and annual average) and two equations for calculating
! the gas transfer coefficient k.
! 
! The NOPE model must register the following as a dependency from the biogeochemical model it is coupled to: 
!  - a nitrification rate (mmolN m-3 d-1)
!  - a denitrification rate (mmolN m-3 d-1)
!  - the dissolved oxygen concentrations (mmolO2 m-3)
! This is achieved in the "coupling" section of NOPE in the fabm.yaml file.
!
!
! SPDX-FileCopyRightText: 2021-2025 Helmholtz-Zentrum hereon GmbH
! SPDX-FileCopyRightText: 2013-2021 Helmholtz-Zentrum Geesthacht GmbH
! SPDX-FileContributor: Nina Preußler <nina.preussler@hereon.de>
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
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
   public type_hereon_nope
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
!
! !REVISION HISTORY:!
!  Original author(s): Nina Preußler
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hereon_nope
!     Variable identifiers
      type (type_state_variable_id)        :: id_n2o_w, id_n2o_b

      type (type_dependency_id)            :: id_temp, id_salinity
      type (type_dependency_id)            :: id_oxy, id_denit, id_nitri
      type (type_surface_dependency_id)    :: id_wind
      type (type_global_dependency_id)     :: id_day_of_year

      type (type_diagnostic_variable_id)   :: id_n2o_from_nitri, id_n2o_from_denit 
      type (type_diagnostic_variable_id)   :: id_n2o_yield_nitri
      type (type_surface_diagnostic_variable_id) :: id_n2o_emit_w, id_n2o_emit_b
      type (type_surface_diagnostic_variable_id) :: id_n2o_eq, id_n2o_sat, id_Sc, id_k_wanninkhof, id_k_borges

!     Model parameters
      real(rk) :: n2o_emission_factor_nitri, n2o_emission_factor_nitri_type, n2o_emission_factor_denit
      real(rk) :: wind_data_type, wind_mean_annual, wind_mean_spring, wind_mean_summer, wind_mean_fall, wind_mean_winter
      real(rk) :: nitrification_scale_factor, denitrification_scale_factor

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_surface

   end type type_hereon_nope
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NOPE model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hereon_nope),intent(inout),target  :: self
   integer,                   intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Nina Preußler
!
! !LOCAL VARIABLES:
      real(rk) :: n2o_emission_factor_nitri, n2o_emission_factor_nitri_type, n2o_emission_factor_denit
      real(rk) :: wind_data_type, wind_mean_annual, wind_mean_spring, wind_mean_summer, wind_mean_fall, wind_mean_winter
 !EOP
!-----------------------------------------------------------------------
!BOC

   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   ! Register parameters
   call self%get_parameter(self%n2o_emission_factor_nitri,'n2o_emission_factor_nitri','','emission factor (yield) of N2O from nitrification',default=0.0025_rk)
   call self%get_parameter(self%n2o_emission_factor_nitri_type,'n2o_emission_factor_nitri_type','','type of emission factor to use for nitrification)',default=2._rk) ! 1=constant (self%n2o_emission_factor_nitri), 2=dynamic, 3=based on Tang et al. (2024)
   call self%get_parameter(self%n2o_emission_factor_denit,'n2o_emission_factor_denit','','emission factor (yield) of N2O from denitrification ',default=0.003_rk)  
   call self%get_parameter(self%wind_data_type,'wind_data_type','-','type of wind data used',default=1._rk) ! 1=hourly, 2=seasonal average, 3=annual average
   call self%get_parameter(self%wind_mean_annual,'wind_mean_annual','m s-1','annual mean of wind speed',default=4.28_rk) 
   call self%get_parameter(self%wind_mean_spring,'wind_mean_spring','m s-1','mean of wind speed in spring',default=4.45_rk)
   call self%get_parameter(self%wind_mean_summer,'wind_mean_summer','m s-1','mean of wind speed in summer',default=3.75_rk)
   call self%get_parameter(self%wind_mean_fall,'wind_mean_fall','m s-1','mean of wind speed in fall',default=3.93_rk)
   call self%get_parameter(self%wind_mean_winter,'wind_mean_winter','m s-1','mean of wind speed in winter',default=5._rk)
   call self%get_parameter(self%nitrification_scale_factor,'nitrification_scale_factor','','scale factor for nitrification',default=1._rk)
   call self%get_parameter(self%denitrification_scale_factor,'denitrification_scale_factor','','scale factor for denitrification',default=1._rk)
   
   ! Register state variables
   call self%register_state_variable(self%id_n2o_w,  'n2o_w',  'umolN2 m-3', 'N2O for sea-to-air flux based on Wanninkhof (1992)',10._rk,minimum=0.0_rk)
   call self%register_state_variable(self%id_n2o_b,  'n2o_b',  'umolN2 m-3', 'N2O for sea-to-air flux based on Borges et al. (2004)',10._rk,minimum=0.0_rk)

   call self%set_variable_property(self%id_n2o_w,'particulate',.false.)
   call self%set_variable_property(self%id_n2o_b,'particulate',.false.)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_n2o_emit_w,'n2o_emit_w','umol m-2 d-1','sea-to-air flux N2O based on Wanninkhof (1992)', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_emit_b,'n2o_emit_b','umol m-2 d-1','sea-to-air flux N2O based on Borges et al. (2004)', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_eq,'n2o_eq','umol m-3','N2O concentration equilibrated with atmosphere', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_k_wanninkhof,'k_wanninkhof','m d-1','gas transfer velocity based on Wanninkhof (1992)', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_k_borges,'k_borges','m d-1','gas transfer velocity based on Borges et al. (2004)', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Sc,'Sc','','Schmidt number', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_sat,'n2o_sat','%','saturation of N2O', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_denit,'n2o_from_denit','umol m-3 d-1','N2O production from denitrification', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_nitri,'n2o_from_nitri','umol m-3 d-1','N2O production from nitrification', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_yield_nitri,'n2o_yield_nitri','','N2O yield from nitrification', output=output_instantaneous)

   ! Register dependencies
   ! from biogeochemical model 
   call self%register_dependency(self%id_denit, 'denit', 'mmolN m**-3 d**-1', 'denitrification rate')
   call self%register_dependency(self%id_nitri, 'nitri', 'mmolN m**-3 d**-1', 'nitrification rate')
   call self%register_dependency(self%id_oxy, 'oxy', 'mmolO2 m-3 d-1', 'dissolved oxygen concentration')
   ! from hydrodynamic host (e.g., GOTM)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salinity,standard_variables%practical_salinity)
   call self%register_dependency(self%id_wind, standard_variables%wind_speed)
   call self%register_dependency(self%id_day_of_year, standard_variables%number_of_days_since_start_of_the_year)

   return

99 call self%fatal_error('hereon_nope_initialize','Error reading namelist hereon_nope.')

100 call self%fatal_error('hereon_nope_initialize','Namelist hereon_nope was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hereon_nope),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Nina Preußler
!
! !LOCAL VARIABLES:
   real(rk) :: oxy
   real(rk) :: denit, nitri
   real(rk) :: n2o_w, n2o_b
   real(rk) :: n2o_from_nitri, n2o_from_denit
   real(rk) :: n2o_yield_nitri_to_use, n2o_yield_nitri_dynamic, n2o_yield_nitri_tang


!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_denit,denit)
   denit = denit *  self%denitrification_scale_factor
   _GET_(self%id_nitri,nitri)
   nitri = nitri * self%nitrification_scale_factor
   _GET_(self%id_oxy,oxy)


   ! N2O from denitrification
   n2o_from_denit = denit * 0.8_rk  * self%n2o_emission_factor_denit ! factor of 0.8 adapted from omexdia
   n2o_from_denit = n2o_from_denit * 1000._rk / 2 

   ! N2O from nitrification
   ! apply desired emission factor/yield
   ! 1 = constant factor as defined with self%n2o_emission_factor_nitri
   ! 2 = factor between 0.1 and 0.4% calculated based on oxygen concentration
   ! 3 = factor/yield based on oxygen concentration using equation by Tang et al. (2024)
   IF (self%n2o_emission_factor_nitri_type == 1.0) THEN
      n2o_yield_nitri_to_use = self%n2o_emission_factor_nitri
   ELSEIF  (self%n2o_emission_factor_nitri_type == 2.0) THEN
      ! calculate yield between 0.1% and 0.4% as measured by de Wilde and de Bie (2000) as a function of oxygen
      n2o_yield_nitri_dynamic = 0.004_rk - 0.003_rk * (oxy/400._rk) 
      n2o_yield_nitri_to_use = n2o_yield_nitri_dynamic
   ELSE
      ! calculate oxygen-dependent N2O yield from nitrification based on Tang et al.(2024)
      n2o_yield_nitri_tang = (0.3889_rk / oxy + 0.2197_rk) / 100._rk 
      n2o_yield_nitri_to_use = n2o_yield_nitri_tang
   END IF
   
   n2o_from_nitri = nitri * n2o_yield_nitri_to_use 
   n2o_from_nitri = n2o_from_nitri * 1000._rk / 2 ! convert into umol and account for stoichiometry (N->N2)
   

#define _CONV_UNIT_ /secs_pr_day

   ! apply production rates
   _ADD_SOURCE_(self%id_n2o_w, (n2o_from_denit + n2o_from_nitri) _CONV_UNIT_)   
   _ADD_SOURCE_(self%id_n2o_b, (n2o_from_denit + n2o_from_nitri) _CONV_UNIT_)   

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_n2o_from_denit,n2o_from_denit)
   _SET_DIAGNOSTIC_(self%id_n2o_from_nitri,n2o_from_nitri)  
   _SET_DIAGNOSTIC_(self%id_n2o_yield_nitri,n2o_yield_nitri_to_use)     

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do



   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)

      class (type_hereon_nope),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: n2o_conc_w, n2o_conc_b
      real(rk) :: n2o_eq, n2o_atm, n2o_sat
      real(rk) :: n2o_sea_air_flux_w, n2o_sea_air_flux_b
      real(rk) :: oxy
      real(rk) :: temp, temp_abs, salinity
      real(rk) :: wind_speed
      real(rk) :: Sc     
      real(rk) :: k_wanninkhof, k_borges    
      real(rk) :: day_of_year

      _SURFACE_LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_n2o_w, n2o_conc_w)
         _GET_(self%id_n2o_b, n2o_conc_b)
         _GET_(self%id_temp, temp)
         _GET_(self%id_salinity, salinity)
         _GET_(self%id_oxy, oxy)
         _GET_SURFACE_(self%id_wind,wind_speed)
         _GET_GLOBAL_(self%id_day_of_year,day_of_year)

         ! adjust to desired wind data type (based on value of self%wind_data_type)
         ! if seasonal wind averages:
         IF (self%wind_data_type == 2._rk) THEN
            ! winter
            IF (day_of_year >= 335._rk .OR. day_of_year <= 59._rk) THEN 
               wind_speed = self%wind_mean_winter
            ! spring
            ELSEIF (day_of_year <= 151._rk) THEN 
               wind_speed = self%wind_mean_spring
            ! summer
            ELSEIF (day_of_year <= 243._rk) THEN 
               wind_speed = self%wind_mean_summer
            ! fall
            ELSE                      
               wind_speed = self%wind_mean_fall
            END IF
         ! if annual mean average:
         ELSEIF (self%wind_data_type == 3._rk) THEN
            wind_speed = self%wind_mean_annual
         END IF

         temp_abs = 273.15_rk + temp ! absolute temperature (in K)

         n2o_atm = 324.2_rk ! global mean atmospheric N2O concentration, based on Yang et al. (2022)
         ! calculate N2O concentration equilibrated with atmosphere (=100% saturation), based on Yang et al. (2022)
         n2o_eq = n2o_atm * exp(-165.8806_rk + 222.8743_rk * (100.0_rk / temp_abs) + 92.0792_rk * log(temp_abs / 100.0_rk) - 1.48425_rk * (temp_abs / 100.0_rk)**2 + salinity * (-0.056235_rk + 0.0316119_rk * (temp_abs / 100.0_rk) - 0.004872 * (temp_abs/100.0_rk)**2))

         ! calculate Schmidt number, based on Wanninkhof (2014)
         Sc = 2141.2_rk + (-152.56_rk) * temp + 5.8963_rk * temp ** 2.0_rk + (-0.12411_rk) * temp ** 3.0_rk + 0.0010655_rk * temp ** 4.0_rk ! based on Wanninkhof (2014)
         
         ! calculate gas transfer coefficients
         ! factor 0.24 adjusts unit to m d-1
         k_wanninkhof = 0.24 * 0.39_rk * wind_speed**2 * (Sc / 660.0_rk)**(-0.5) ! based on Wanninkhof (1992)
         k_borges = 0.24 * (4.045_rk+2.58_rk*wind_speed) * (Sc/600._rk)**(-0.5)  ! based on Borges et al. (2004)
         
         ! calculate sea-to-air flux
         n2o_sea_air_flux_w = (n2o_conc_w - n2o_eq) * k_wanninkhof
         n2o_sea_air_flux_b = (n2o_conc_b - n2o_eq) * k_borges

         ! calculate saturation of N2O
         n2o_sat = (n2o_conc_b / n2o_eq) * 100.0_rk

         ! apply sea-to-air flux
         _ADD_SURFACE_FLUX_(self%id_n2o_b,-n2o_sea_air_flux_b _CONV_UNIT_)
         _ADD_SURFACE_FLUX_(self%id_n2o_w,-n2o_sea_air_flux_w _CONV_UNIT_)   

         ! Export diagnostic variables
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_emit_w, n2o_sea_air_flux_w)
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_emit_b, n2o_sea_air_flux_b)
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_eq, n2o_eq)
         _SET_SURFACE_DIAGNOSTIC_(self%id_k_wanninkhof, k_wanninkhof)
         _SET_SURFACE_DIAGNOSTIC_(self%id_k_borges, k_borges)
         _SET_SURFACE_DIAGNOSTIC_(self%id_Sc, Sc)
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_sat, n2o_sat)      

      _SURFACE_LOOP_END_

   end subroutine do_surface   

!EOC

   end module hereon_nope

