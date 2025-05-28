#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hereon_nope --- Fortran 2003 version of OMEXDIA+P biogeochemical model
!
! !INTERFACE:
   module hereon_nope
!
! !DESCRIPTION:
!
! This module can be coupled with estuarine and coastal biogeochemical models to simulate 
! nitrous oxide production and emission. It includes multiple alternative calculations of production
! and emission rates based on different approaches used in the literature. The module was developed
! for water column processes but could also be adapted to sediment processes.
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
      type (type_state_variable_id)        :: id_n2o_w, id_n2o_b, id_n2o_tang

      type (type_dependency_id)            :: id_temp, id_salinity
      type (type_dependency_id)            :: id_oxy, id_denit, id_nitri

      type (type_surface_dependency_id)    :: id_wind
      type (type_global_dependency_id)     :: id_day_of_year

      type (type_diagnostic_variable_id)   :: id_n2o_from_nitri, id_n2o_from_denit 
      type (type_diagnostic_variable_id)   :: id_n2o_from_nitri_tang, id_n2o_from_denit_tang, id_n2o_yield_nitri
      type (type_surface_diagnostic_variable_id) :: id_n2o_emit_w, id_n2o_emit_b, id_n2o_eq, id_n2o_sat
      type (type_surface_diagnostic_variable_id) :: id_k_wanninkhof, id_k_borges, id_Sc

!     Model parameters
      real(rk) :: n2o_emission_factor_nitri, n2o_emission_factor_denit
      real(rk) :: wind_data_type
      real(rk) :: wind_mean_annual, wind_mean_spring, wind_mean_summer, wind_mean_fall, wind_mean_winter

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
   class (type_hereon_nope),intent(inout),target  :: self
   integer,                   intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Nina Preußler
!
! !LOCAL VARIABLES:
      real(rk) :: n2o_emission_factor_nitri, n2o_emission_factor_denit
      real(rk) :: wind_data_type
      real(rk) :: wind_mean_annual, wind_mean_spring, wind_mean_summer, wind_mean_fall, wind_mean_winter


!EOP
!-----------------------------------------------------------------------
!BOC

   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk


   ! Register parameters
   call self%get_parameter(self%n2o_emission_factor_nitri,'n2o_emission_factor_nitri','','EF nitrification ',default=0.0025_rk)
   call self%get_parameter(self%n2o_emission_factor_denit,'n2o_emission_factor_denit','','EF denitation ',default=0.01_rk)  
   call self%get_parameter(self%wind_data_type,'wind_data_type','-','type of wind data used (hourly-1, seasonal-2, annual-3) ',default=1._rk)
   call self%get_parameter(self%wind_mean_annual,'wind_mean_annual','m s**-1','yearly mean of wind speed',default=4.28_rk)
   call self%get_parameter(self%wind_mean_spring,'wind_mean_spring','m s**-1','mean of wind speed in spring',default=1._rk)
   call self%get_parameter(self%wind_mean_summer,'wind_mean_summer','m s**-1','mean of wind speed in summer',default=2._rk)
   call self%get_parameter(self%wind_mean_fall,'wind_mean_fall','m s**-1','mean of wind speed in fall',default=3._rk)
   call self%get_parameter(self%wind_mean_winter,'wind_mean_winter','m s**-1','mean of wind speed in winter',default=4._rk)


   ! Register state variables
   call self%register_state_variable(self%id_n2o_w,  'n2o_w',  'umolN2 m**-3', 'dissolved nitrous oxide for water-to-air flux based on Wanninkhof',      10._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_n2o_b,  'n2o_b',  'umolN2 m**-3', 'dissolved nitrous oxide for water-to-air flux based on Borges',      10._rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_n2o_tang,  'n2o_tang',  'umolN2 m**-3', 'dissolved nitrous oxide based on Tang',      10._rk, minimum=0.0_rk)

   call self%set_variable_property(self%id_n2o_w,'particulate',.false.)
   call self%set_variable_property(self%id_n2o_b,'particulate',.false.)
   call self%set_variable_property(self%id_n2o_tang,'particulate',.false.)


   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_n2o_emit_w,'n2o_emit_w','umol m**-2 d-1','rate of N2O emission based on Wanninkhof', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_emit_b,'n2o_emit_b','umol m**-2 d-1','rate of N2O emission based on Borges', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_eq,'n2o_eq','umol m**-3','N2O concentration equilibrated with atmosphere', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_k_wanninkhof,'k_wanninkhof','cm h**-1','gas transfer velocity calculated with Wanninkhof equation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_k_borges,'k_borges','cm h**-1','gas transfer velocity based on Borges', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Sc,'Sc','','Schmidt number', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_sat,'n2o_sat','%','saturation of N2O', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_denit,'n2o_from_denit','umol m**-3 d-1','N2O production from denitation ', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_nitri,'n2o_from_nitri','umol m**-3 d-1','N2O production from nitrification ', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_denit_tang,'n2o_from_denit_tang','umol m**-3 d-1','N2O production from denitation ', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_nitri_tang,'n2o_from_nitri_tang','umol m**-3 d-1','N2O production from nitrification ', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_yield_nitri,'n2o_yield_nitri','','oxygen-dependent N2O yield from nitrification based on Tang et al. (2022)', output=output_instantaneous)


   ! Register dependencies
   ! from biogeochemical model
   call self%register_dependency(self%id_denit, 'denit', 'mmolN m**-3 d**-1', 'denitation rate')
   call self%register_dependency(self%id_nitri, 'nitri', 'mmolN m**-3 d**-1', 'nitrification rate')
   call self%register_dependency(self%id_oxy, 'oxy', 'mmolO2 m**-3 d**-1', 'dissolved oxygen concentration')
   ! from hydrodynamic host (e.g., GOTM)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salinity,standard_variables%practical_salinity)
   call self%register_dependency(self%id_wind, standard_variables%wind_speed)
   call self%register_dependency(self%id_day_of_year, standard_variables%number_of_days_since_start_of_the_year)
   !call self%register_dependency(self%id_total_nitrogen, standard_variables%total_nitrogen)

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
   real(rk) :: n2o_w, n2o_b, n2o_tang
   real(rk) :: n2o_from_nitri, n2o_from_denit, n2o_yield_nitri
   real(rk) :: n2o_from_nitri_tang, n2o_from_denit_tang
   real(rk) :: n2o_yield_nitri_2


!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_denit,denit)
   _GET_(self%id_nitri,nitri)
   _GET_(self%id_oxy,oxy)
   !_GET_(self%id_n2o_w,n2o_w)
   !_GET_(self%id_n2o_b,n2o_b)
   !_GET_(self%id_n2o_tang,n2o_tang)


   ! N2O from denitrification
   ! reference + explanation numbers
   n2o_from_denit = denit * self%n2o_emission_factor_denit * 0.8_rk * 1000._rk / 2 ! see 0.8 below for nitrate reduction, but why? ! emission factor from Beaulieau et al. (2011)
   n2o_from_denit_tang = 1.5_rk * exp(-0.17_rk*oxy) ! oxygen-dependent production from denitri
   n2o_from_denit_tang = n2o_from_denit_tang - 2.5_rk * exp(-0.47_rk*oxy) ! oxygen-dependent consumption of n2o via complete denitation

   ! N2O from nitrification
   n2o_yield_nitri = (1.52_rk / (oxy + 1.59_rk)) / 100._rk ! oxygen-dependent N2O yield from nitrification based on Tang et al. 2022
   ! n2o_from_nitri = nitri * self%n2o_emission_factor_nitri * 1000._rk ! 1000 because nitri in mmol m3 and n2o in umol m-3
   n2o_yield_nitri_2 = 0.004_rk - 0.003_rk * (oxy/350) ! to get a yield between 0.1% and 0.4%
   n2o_from_nitri = nitri * self%n2o_emission_factor_nitri * 1000._rk / 2 ! 1000 because nitri in mmol m3 and n2o in umol m-3, and /2 because 2 N
   n2o_from_nitri_tang = ( 1000._rk * oxy / (oxy+4.3_rk) ) * n2o_yield_nitri ! nitri * yield from Tang


#define _CONV_UNIT_ /secs_pr_day

! reaction rates
   _ADD_SOURCE_(self%id_n2o_w, (n2o_from_denit + n2o_from_nitri) _CONV_UNIT_)   
   _ADD_SOURCE_(self%id_n2o_b, (n2o_from_denit + n2o_from_nitri) _CONV_UNIT_)   
   _ADD_SOURCE_(self%id_n2o_tang, (n2o_from_denit_tang + n2o_from_nitri_tang) _CONV_UNIT_) 
   !_ADD_SOURCE_(self%id_no3 , 0 _CONV_UNIT_)   how can I maintain mass balance? 
   !_ADD_SOURCE_(self%id_nh3 , 0 _CONV_UNIT_)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_n2o_from_denit,n2o_from_denit)
   _SET_DIAGNOSTIC_(self%id_n2o_from_nitri,n2o_from_nitri)  
   _SET_DIAGNOSTIC_(self%id_n2o_yield_nitri,n2o_yield_nitri)     
   _SET_DIAGNOSTIC_(self%id_n2o_from_denit_tang,n2o_from_denit_tang)
   _SET_DIAGNOSTIC_(self%id_n2o_from_nitri_tang,n2o_from_nitri_tang)    

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do



   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)

      class (type_hereon_omexdia_n2o),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: n2o_conc_w, n2o_conc_b
      real(rk) :: n2o_eq, n2o_atm, n2o_sat
      real(rk) :: n2o_sea_air_flux_w, n2o_sea_air_flux_b
      real(rk) :: oxy
      real(rk) :: temp, temp_abs, salinity
      real(rk) :: wind_speed
      real(rk) :: Sc       ! Schmidt number
      real(rk) :: k_wanninkhof, k_borges        ! gas transfer coefficient
      real(rk) :: day_of_year

      _SURFACE_LOOP_BEGIN_

         ! Retrieve current (local) state variable values.

         ! Retrieve current (local) state variable values.
         _GET_(self%id_n2o_w, n2o_conc_w)
         _GET_(self%id_n2o_b, n2o_conc_b)
         _GET_(self%id_temp, temp)
         _GET_(self%id_salinity, salinity)
         _GET_(self%id_oxy, oxy)
         _GET_SURFACE_(self%id_wind,wind_speed)
         _GET_GLOBAL_(self%id_day_of_year,day_of_year)

         ! adjust to desired wind data type (based on value of wind_data_type_paramter)
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
            wind_speed = wind_mean_annual
         END IF

         temp_abs = 273.15_rk + temp ! absolute temperature (in K)

         n2o_atm = 324.2_rk ! global mean atmospheric N2O concentration, based on Yang et al. (2022)
         n2o_eq = n2o_atm * exp(-165.8806_rk + 222.8743_rk * (100.0_rk / temp_abs) + 92.0792_rk * log(temp_abs / 100.0_rk) - 1.48425_rk * (temp_abs / 100.0_rk)**2 + salinity * (-0.056235_rk + 0.0316119_rk * (temp_abs / 100.0_rk) - 0.004872 * (temp_abs/100.0_rk)**2))
         ! n2o_eq in nmol L-1 (=umol m-3)

         Sc = 2141.2_rk + (-152.56_rk) * temp + 5.8963_rk * temp ** 2.0_rk + (-0.12411_rk) * temp ** 3.0_rk + 0.0010655_rk * temp ** 4.0_rk ! based on Wanninkhof (2014)
         
         !wind_speed = 4.28_rk

         k_wanninkhof = 0.39_rk * wind_speed**2 * (Sc / 660.0_rk)**(-0.5) ! based on Wanninkhof
         k_borges = 0.24_rk * (4.045_rk+2.58_rk*wind_speed) * (Sc/600._rk)**(-0.5) ! based on Borges (as used by Brase et al., 2017)
         
         n2o_sea_air_flux_w = (n2o_conc_w - n2o_eq) * k_wanninkhof
         n2o_sea_air_flux_b = (n2o_conc_b - n2o_eq) * k_borges

         n2o_sat = (n2o_conc_w / n2o_eq) * 100.0_rk

         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_emit_w, n2o_sea_air_flux_w)
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_emit_b, n2o_sea_air_flux_b)
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_eq, n2o_eq)
         _SET_SURFACE_DIAGNOSTIC_(self%id_k_wanninkhof, k_wanninkhof)
         _SET_SURFACE_DIAGNOSTIC_(self%id_k_borges, k_borges)
         _SET_SURFACE_DIAGNOSTIC_(self%id_Sc, Sc)
         _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_sat, n2o_sat)

         _ADD_SURFACE_FLUX_(self%id_n2o_b,-n2o_sea_air_flux_b _CONV_UNIT_)
         _ADD_SURFACE_FLUX_(self%id_n2o_w,-n2o_sea_air_flux_w _CONV_UNIT_)
         _ADD_SURFACE_FLUX_(self%id_n2o_tang,-n2o_sea_air_flux_b _CONV_UNIT_)
         

      _SURFACE_LOOP_END_

   end subroutine do_surface   

!EOC

   end module hereon_nope

