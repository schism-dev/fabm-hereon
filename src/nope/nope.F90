! SPDX-FileCopyrightText: 2024-2026 Helmholtz-Zentrum hereon GmbH
! SPDX-License-Identifier: Apache-2.0

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hereon_nope --- N2O production, consumption, and emission
!
! !INTERFACE:
   module hereon_nope
!
! !DESCRIPTION:
!
! NOPE (Nitrous Oxide Production and Emission) is a water-column N2O module
! for FABM. It takes nitrification and denitrification rates and the O2, NO3
! and (optionally) labile-carbon concentrations from a host biogeochemical
! model and computes dissolved N2O, its production and reduction, and the
! sea-to-air flux.
!
! Sources
! 1) Nitrification, two pathways (Wrage et al. 2001):
!    (a) hydroxylamine oxidation: small, nearly constant yield (Goreau et al.
!        1980; Frame and Casciotti 2010; Santoro et al. 2011);
!    (b) nitrifier denitrification: yield rises sharply as O2 declines
!        (Goreau et al. 1980; Ji et al. 2018).
!      yield     = y_nh2oh + y_nd_max * K_nd^p / (O2^p + K_nd^p)
!      N2O_nitri = nitrification_rate * yield
!    NH4 limitation is part of the host rate and is not applied again.
! 2) Heterotrophic denitrification, "leaky-pipe" yield peaked at suboxic O2
!    (Firestone and Davidson 1989):
!      yield = y_floor + (y_peak - y_floor) * normalised_shape(O2)
!      act   = O2^p / (O2^p + Kact^p)   (low O2: reduction runs on to N2)
!      inh   = Kinh^p / (Kinh^p + O2^p) (high O2: extra yield decreases)
!      N2O_denit = denitrification_rate * yield
!
! Sink
! 3) Reduction to N2 by denitrifiers:
!      reduction = k_red * o2_inh * no3_inh * c_lim * Q10(T) * [N2O]
!    (a) sigmoidal O2 inhibition; N2O reductase is more O2-sensitive than
!        nitrate reductase (Wicht 1996; Bonin et al. 1989);
!    (b) NO3 competition, K / (K + NO3), an assumed form (neither the form
!        nor K is taken from the literature);
!    (c) optional electron-donor limitation, C / (k_c + C), with a
!        host-supplied labile-carbon proxy (enable_c_limit_red);
!    (d) Q10 temperature scaling, internal to NOPE.
!
! Boundary fluxes
! 4) Air-water exchange (do_surface). Gas transfer velocity from Wanninkhof
!    (2014), Nightingale et al. (2000) or Wanninkhof (1992) with the Raymond
!    and Cole (2001) low-wind branch (gas_transfer_scheme = 1, 2, 3), times
!    gas_transfer_factor for scheme 3. Solubility after Weiss and Price
!    (1980), with optional water-vapour correction.
! 5) Optional O2- and temperature-dependent benthic efflux (do_bottom,
!    enable_benthic_n2o, off by default).
!
! Temperature
!   Host rates usually include temperature already (OxyPOM through Q10
!   factors, OMEXDIA through f_T), so production is not scaled by default
!   (apply_prod_temperature_scaling = false); a Q10 here would count it
!   twice. Set it true only for hosts that export temperature-independent
!   rates. The Q10 of N2O reduction is always applied, because the
!   reduction is computed inside NOPE.
!
! Dependencies (units)
!   nitri     host nitrification rate           [host units]
!   denit     host denitrification rate         [host units]
!   oxy       dissolved oxygen                  [mmolO2 m-3]
!   no3       nitrate (reduction only)          [mmolN m-3]
!   labile_c  labile carbon (C limitation only) [mmolC m-3]
!   nitri_input_to_mmol_n and denit_input_to_mmol_n convert the host rates
!   to mmol N m-3 d-1. Yields are mol N2O-N per mol N processed, so one
!   mmol N gives 500 umol N2O.
!
! Meaning of the host rates
!   nitri  gross NH4 oxidation, with the host NH4, O2 and temperature
!          limitation already applied. NOPE does not apply them again
!          (except the optional production Q10, off by default).
!   denit  gross NO3 reduction by heterotrophic respiration, including the
!          host O2 dependence. NOPE sets only the fraction released as N2O
!          and adds no threshold or O2 dependence to the rate itself.
!   The yields are calibrated for these definitions. A host whose denit
!   means something else (anaerobic C mineralisation, which needs a C:N
!   conversion; NO3 to NO2 only; a built-in O2 threshold) runs without
!   error, but the yields then describe a different N2O budget, which NOPE
!   cannot detect. Recalibrate the yields against local observations for
!   each host, and compare the production-weighted mean O2 per pathway and
!   the realised minimum and maximum yields with the original calibration.
!
! Mass budget
!   NOPE does not close the host N budget. The only _ADD_SOURCE_ target is
!   id_n2o_w: N2O-N is not removed from the host NO3 and NH4 pools, and
!   reduced N2O goes to an untracked N2 pool. Two-way feedback on the host
!   N pools is outside the scope of v1.0. At Seemannshoeft (2015), annual
!   N2O-N production is 0.36 % of the N processed by nitrification and
!   denitrification (column integral, hourly), below the largest realised
!   total yield (1.49 %, nitrification). Hosts with larger yields or smaller
!   N pools should recompute this.
!
! SPDX-FileCopyrightText: 2021-2026 Helmholtz-Zentrum hereon GmbH
! SPDX-FileContributor: Nina Preußler <nina.preussler@hereon.de>
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-FileContributor: Youssef Saadaoui <youssef.saadaoui@hereon.de>
! SPDX-License-Identifier: Apache-2.0
!
!
! !USES:
   use fabm_types

   implicit none

   private

   public type_hereon_nope

   real(rk), parameter :: per_day_to_per_second = 1.0_rk / 86400.0_rk
   real(rk), parameter :: mmol_n_to_umol_n2o = 500.0_rk   ! 1 mmol N2O-N = 0.5 mmol N2O = 500 umol N2O
   character(len=*), parameter :: nope_version = '1.0'

   type,extends(type_base_model) :: type_hereon_nope
      ! ---- State variable ----
      type (type_state_variable_id)        :: id_n2o_w

      ! ---- Dependencies ----
      type (type_dependency_id)            :: id_temp, id_salinity
      type (type_dependency_id)            :: id_oxy, id_denit, id_nitri
      type (type_dependency_id)            :: id_no3                        ! NO3 competition in N2O reduction
      type (type_dependency_id)            :: id_labile_c                   ! electron-donor proxy for reduction
      type (type_surface_dependency_id)    :: id_wind
      type (type_surface_dependency_id)    :: id_air_pressure
      type (type_global_dependency_id)     :: id_day_of_year

      ! ---- Diagnostics: production & consumption ----
      type (type_diagnostic_variable_id)   :: id_n2o_from_nitri, id_n2o_from_denit
      type (type_diagnostic_variable_id)   :: id_n2o_nitri_nh2oh, id_n2o_nitri_nd  ! two nitrification pathways
      type (type_diagnostic_variable_id)   :: id_n2o_yield_nitri, id_n2o_yield_denit
      type (type_diagnostic_variable_id)   :: id_n2o_consumption
      type (type_diagnostic_variable_id)   :: id_n2o_net_production
      type (type_diagnostic_variable_id)   :: id_c_lim                      ! electron-donor limitation factor
      type (type_diagnostic_variable_id)   :: id_nitri_rate_used, id_denit_rate_used
      type (type_diagnostic_variable_id)   :: id_denit_o2_activation
      type (type_diagnostic_variable_id)   :: id_denit_o2_inhibition
      type (type_diagnostic_variable_id)   :: id_n2o_reduction_o2_factor
      type (type_diagnostic_variable_id)   :: id_n2o_reduction_no3_factor

      ! ---- Diagnostics: surface/emission ----
      type (type_surface_diagnostic_variable_id) :: id_n2o_emit_w
      type (type_surface_diagnostic_variable_id) :: id_n2o_eq, id_n2o_sat
      type (type_surface_diagnostic_variable_id) :: id_Sc, id_k_gas

      ! ---- Diagnostics: benthic ----
      type (type_bottom_diagnostic_variable_id)  :: id_n2o_benthic

      ! ---- Host diagnostic conversion to mmol N m-3 d-1 ----
      ! For OxyPOM: x_nit = 2*NIT and x_denit = 1.5*MINnit.
      real(rk) :: nitri_input_to_mmol_n
      real(rk) :: denit_input_to_mmol_n

      ! ---- Production temperature control (for hosts without temperature in their rates) ----
      logical :: apply_prod_temperature_scaling
      real(rk) :: q10_prod                ! Q10 of N2O production (only if apply_prod_temperature_scaling)
      real(rk) :: t_ref_prod              ! reference temperature for production Q10 [degC]

      ! ---- Nitrification N2O yield (two pathways) ----
      real(rk) :: y_nh2oh                 ! hydroxylamine-oxidation yield (O2-independent floor)
      real(rk) :: y_nd_max                ! maximum nitrifier-denitrification yield (low O2)
      real(rk) :: oxy_half_y_nitri, oxy_pow_y_nitri

      ! ---- Denitrification N2O yield (peaked leaky-pipe) ----
      real(rk) :: y_denit_floor, y_denit_peak
      real(rk) :: oxy_act_y_denit         ! activation half-saturation Kact (low O2: reduction to N2)
      real(rk) :: oxy_half_y_denit        ! inhibition half-saturation Kinh (high O2)
      real(rk) :: oxy_pow_y_denit
      integer  :: denit_floor_mode        ! 0 = constant floor (baseline); 1 = oxygen-gated floor (sensitivity)

      ! ---- N2O reduction ----
      logical :: enable_n2o_reduction
      real(rk) :: k_red_d_max             ! maximum reduction rate [d-1]
      real(rk) :: o2_red_mid              ! O2 midpoint of the sigmoidal inhibition [mmolO2 m-3]
      real(rk) :: o2_red_width            ! width of the O2 sigmoid [mmolO2 m-3]
      real(rk) :: kin_no3_n2ored          ! NO3 half-inhibition of N2O reduction [mmolN m-3]
      logical :: enable_c_limit_red       ! electron-donor limitation switch
      real(rk) :: k_c_red                 ! half-saturation of labile-C limitation [mmolC m-3]
      real(rk) :: q10_red                 ! Q10 for reduction
      real(rk) :: t_ref_red               ! reference temperature [degC]

      ! ---- Benthic N2O efflux (optional, experimental) ----
      logical :: enable_benthic_n2o
      real(rk) :: benthic_n2o_base        ! base sediment N2O efflux [umol m-2 d-1]
      real(rk) :: benthic_n2o_tcoef       ! temperature coefficient [umol m-2 d-1 degC-1]
      real(rk) :: benthic_o2_half         ! O2 half-saturation for the low-O2 efflux enhancement [mmolO2 m-3]

      ! ---- Atmosphere & air-sea exchange ----
      real(rk) :: n2o_atm_ppb             ! atmospheric N2O mixing ratio [ppb]
      real(rk) :: patm_atm                ! air pressure if not read from forcing [atm]
      logical :: use_h2o_correction
      logical :: use_surface_air_pressure
      integer :: gas_transfer_scheme      ! 1 = Wanninkhof 2014, 2 = Nightingale 2000, 3 = Wanninkhof 1992
      real(rk) :: gas_transfer_factor   ! multiplier on k (scheme 3 only)
      integer :: schmidt_water_type      ! 1 = freshwater (default), 2 = seawater N2O Schmidt polynomial

      ! ---- Wind ----
      integer :: wind_forcing_mode        ! 1 = forcing, 2 = seasonal means, 3 = annual mean
      real(rk) :: wind_mean_annual, wind_mean_spring, wind_mean_summer
      real(rk) :: wind_mean_fall, wind_mean_winter

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom

   end type type_hereon_nope

   contains

!-----------------------------------------------------------------------
! Read parameters, register the N2O state variable, diagnostics and
! dependencies.
!-----------------------------------------------------------------------
   subroutine initialize(self, configunit)

   class (type_hereon_nope), intent(inout), target :: self
   integer,                  intent(in)            :: configunit

   ! --- Convert host diagnostics to mmol N m-3 d-1 ---
   call self%get_parameter(self%nitri_input_to_mmol_n, &
      'nitri_input_to_mmol_n', 'mmolN mmol-1', &
      'conversion of host nitrification diagnostic to mmol N', &
      default=1.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%denit_input_to_mmol_n, &
      'denit_input_to_mmol_n', 'mmolN mmol-1', &
      'conversion of host denitrification diagnostic to mmol N', &
      default=1.0_rk, minimum=0.0_rk)

   ! --- Production temperature control ---
   ! false (default): host rates already include temperature. true: apply
   ! q10_prod to production (hosts with temperature-independent rates only).
   call self%get_parameter(self%apply_prod_temperature_scaling, &
      'apply_prod_temperature_scaling', '', &
      'apply internal Q10 scaling to N2O production', default=.false.)
   call self%get_parameter(self%q10_prod, &
      'q10_prod', '', 'Q10 for N2O production when internal scaling is enabled', &
      default=1.5_rk, minimum=1.0e-12_rk)
   call self%get_parameter(self%t_ref_prod, &
      't_ref_prod', 'degC', 'reference temperature for production Q10', default=20.0_rk)

   ! --- Nitrification N2O yield (two pathways) ---
   ! yield = y_nh2oh + y_nd_max * K^p/(O2^p+K^p)
   call self%get_parameter(self%y_nh2oh, &
      'y_nh2oh', '', 'N2O-N yield via hydroxylamine oxidation (high-O2 floor)', &
      default=0.0005_rk, minimum=0.0_rk)
   call self%get_parameter(self%y_nd_max, &
      'y_nd_max', '', 'max additional N2O-N yield via nitrifier-denitrification (low O2)', &
      default=0.0040_rk, minimum=0.0_rk)
   call self%get_parameter(self%oxy_half_y_nitri, &
      'oxy_half_y_nitri', 'mmolO2 m-3', 'O2 half-sat for nitrifier-denitrification yield', &
      default=20.0_rk, minimum=1.0e-12_rk)
   call self%get_parameter(self%oxy_pow_y_nitri, &
      'oxy_pow_y_nitri', '', 'Hill exponent for nitrification yield curve', &
      default=1.0_rk, minimum=1.0e-12_rk)

   ! --- Denitrification N2O yield (peaked "leaky-pipe") ---
   ! The normalised shape peaks at one, so y_denit_peak is the maximum yield.
   call self%get_parameter(self%y_denit_floor, &
      'y_denit_floor', '', 'minimum N2O-N yield from denitrification', &
      default=0.0015_rk, minimum=0.0_rk)
   call self%get_parameter(self%y_denit_peak, &
      'y_denit_peak', '', 'actual peak N2O-N yield from denitrification', &
      default=0.0060_rk, minimum=0.0_rk)
   call self%get_parameter(self%oxy_act_y_denit, &
      'oxy_act_y_denit', 'mmolO2 m-3', &
      'activation O2 half-sat Kact for the N2O yield', &
      default=20.0_rk, minimum=1.0e-12_rk)
   call self%get_parameter(self%oxy_half_y_denit, &
      'oxy_half_y_denit', 'mmolO2 m-3', &
      'inhibition O2 half-sat Kinh for the N2O yield', &
      default=5.0_rk, minimum=1.0e-12_rk)
   call self%get_parameter(self%oxy_pow_y_denit, &
      'oxy_pow_y_denit', '', 'Hill exponent for denitrification yield curve', &
      default=1.0_rk, minimum=1.0e-12_rk)
   if (self%y_denit_peak < self%y_denit_floor) then
      call self%fatal_error('hereon_nope_initialize', &
         'y_denit_peak must be greater than or equal to y_denit_floor.')
   end if
   ! denit_floor_mode: treatment of the floor at high O2.
   !   0 = constant floor (baseline): y_denit -> y_denit_floor as O2 -> inf.
   !   1 = oxygen-gated floor (sensitivity test): floor times the high-O2
   !       inhibition term, so y_denit -> 0 as O2 -> inf.
   ! Both modes peak at y_denit_peak at O2* = sqrt(oxy_act_y_denit*oxy_half_y_denit).
   call self%get_parameter(self%denit_floor_mode, &
      'denit_floor_mode', '', &
      '0=constant floor (baseline), 1=oxygen-gated floor (structural sensitivity)', &
      default=0, minimum=0, maximum=1)

   ! --- N2O reduction to N2 ---
   call self%get_parameter(self%enable_n2o_reduction, &
      'enable_n2o_reduction', '', 'enable N2O reduction', default=.true.)
   call self%get_parameter(self%k_red_d_max, &
      'k_red_d_max', 'd-1', 'max specific N2O reduction rate', &
      default=0.02_rk, minimum=0.0_rk)
   call self%get_parameter(self%o2_red_mid, &
      'o2_red_mid', 'mmolO2 m-3', &
      'O2 midpoint for sigmoidal inhibition of N2O reduction', &
      default=2.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%o2_red_width, &
      'o2_red_width', 'mmolO2 m-3', 'width of O2 inhibition sigmoid', &
      default=1.0_rk, minimum=1.0e-12_rk)
   ! NO3 competition K/(K+NO3): abundant NO3 is reduced in preference to N2O.
   ! Assumed form and K (header, item 3b).
   call self%get_parameter(self%kin_no3_n2ored, &
      'kin_no3_n2ored', 'mmolN m-3', &
      'NO3 half-inhibition of N2O reduction (high NO3 suppresses N2O->N2)', &
      default=50.0_rk, minimum=1.0e-12_rk)
   ! Electron-donor limitation: N2O reduction is heterotrophic and needs labile
   ! carbon, c_lim = C/(k_c_red + C), with C coupled through 'labile_c'.
   call self%get_parameter(self%enable_c_limit_red, &
      'enable_c_limit_red', '', &
      'apply labile-C (electron-donor) limitation to N2O reduction', default=.true.)
   call self%get_parameter(self%k_c_red, &
      'k_c_red', 'mmolC m-3', &
      'half-saturation of labile-C limitation of N2O reduction', &
      default=10.0_rk, minimum=1.0e-12_rk)
   call self%get_parameter(self%q10_red, &
      'q10_red', '', 'Q10 for N2O reduction', &
      default=2.0_rk, minimum=1.0e-12_rk)
   call self%get_parameter(self%t_ref_red, &
      't_ref_red', 'degC', 'reference temperature for Q10', default=20.0_rk)

   ! --- Benthic N2O efflux ---
   ! Empirical source for sensitivity tests, not a sediment model; needs
   ! site-specific constraints and is off by default.
   call self%get_parameter(self%enable_benthic_n2o, &
      'enable_benthic_n2o', '', 'enable experimental benthic N2O efflux', default=.false.)
   call self%get_parameter(self%benthic_n2o_base, &
      'benthic_n2o_base', 'umol m-2 d-1', 'base sediment N2O efflux at 0 degC', &
      default=2.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%benthic_n2o_tcoef, &
      'benthic_n2o_tcoef', 'umol m-2 d-1 degC-1', &
      'temperature coefficient of benthic N2O efflux', &
      default=0.2_rk, minimum=0.0_rk)
   call self%get_parameter(self%benthic_o2_half, &
      'benthic_o2_half', 'mmolO2 m-3', &
      'bottom-O2 half-sat enhancing benthic N2O efflux at low O2', &
      default=62.5_rk, minimum=0.0_rk)

   ! --- Atmosphere & air-sea exchange ---
   call self%get_parameter(self%n2o_atm_ppb, &
      'n2o_atm_ppb', 'ppb', 'atmospheric N2O mixing ratio', &
      default=330.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%patm_atm, &
      'patm_atm', 'atm', 'constant atmospheric pressure fallback', &
      default=1.0_rk, minimum=1.0e-12_rk)
   call self%get_parameter(self%use_h2o_correction, &
      'use_h2o_correction', '', 'apply water-vapour pressure correction', default=.true.)
   call self%get_parameter(self%use_surface_air_pressure, &
      'use_surface_air_pressure', '', &
      'use time-varying surface air pressure', default=.false.)
   call self%get_parameter(self%gas_transfer_scheme, &
      'gas_transfer_scheme', '', &
      '1=Wanninkhof 2014, 2=Nightingale 2000, 3=Wanninkhof 1992 with low-wind branch (OxyPOM form)', &
      default=1, minimum=1, maximum=3)
   ! Scheme 3 is the OxyPOM O2 reaeration formula with the N2O Schmidt number,
   ! so O2 and N2O share one transfer velocity. gas_transfer_factor plays the
   ! role of the OxyPOM factor alpha; schemes 1 and 2 ignore it.
   call self%get_parameter(self%gas_transfer_factor, &
      'gas_transfer_factor', '-', &
      'multiplier on the gas transfer velocity (scheme 3 only; set equal to the host O2 factor)', &
      default=1.0_rk, minimum=0.0_rk)
   ! N2O Schmidt-number polynomial (Wanninkhof 2014, Table 1). Freshwater (1)
   ! suits Seemannshoeft (S ~ 0.5); use 2 for marine sites.
   call self%get_parameter(self%schmidt_water_type, &
      'schmidt_water_type', '', &
      'N2O Schmidt-number polynomial: 1=freshwater (Wanninkhof 2014 Table 1), 2=seawater', &
      default=1, minimum=1, maximum=2)

   ! --- Wind ---
   ! Mode 1 (default) reads the FABM standard wind_speed. wind_mean_* are used
   ! only in modes 2 and 3; their defaults are placeholders, not a climatology.
   call self%get_parameter(self%wind_forcing_mode, &
      'wind_forcing_mode', '', &
      '1=forcing, 2=seasonal climatology, 3=annual mean', &
      default=1, minimum=1, maximum=3)
   call self%get_parameter(self%wind_mean_annual, 'wind_mean_annual', 'm s-1', &
      'annual mean wind speed (placeholder; must be set per site in modes 2/3)', &
      default=6.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%wind_mean_spring, 'wind_mean_spring', 'm s-1', &
      'spring mean wind speed (placeholder; must be set per site in mode 2)', &
      default=6.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%wind_mean_summer, 'wind_mean_summer', 'm s-1', &
      'summer mean wind speed (placeholder; must be set per site in mode 2)', &
      default=6.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%wind_mean_fall, 'wind_mean_fall', 'm s-1', &
      'fall mean wind speed (placeholder; must be set per site in mode 2)', &
      default=6.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%wind_mean_winter, 'wind_mean_winter', 'm s-1', &
      'winter mean wind speed (placeholder; must be set per site in mode 2)', &
      default=6.0_rk, minimum=0.0_rk)

   ! === State variable ===
   call self%register_state_variable(self%id_n2o_w, 'n2o_w', 'umol m-3', &
      'dissolved N2O concentration', 30.0_rk, minimum=0.0_rk)
   call self%set_variable_property(self%id_n2o_w, 'particulate', .false.)

   ! === Diagnostics ===
   call self%register_diagnostic_variable(self%id_n2o_from_nitri, &
      'n2o_from_nitri', 'umol m-3 d-1', 'N2O production from nitrification', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_nitri_nh2oh, &
      'n2o_nitri_nh2oh', 'umol m-3 d-1', 'N2O from hydroxylamine-oxidation pathway', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_nitri_nd, &
      'n2o_nitri_nd', 'umol m-3 d-1', 'N2O from nitrifier-denitrification pathway', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_from_denit, &
      'n2o_from_denit', 'umol m-3 d-1', 'N2O production from denitrification', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_yield_nitri, &
      'n2o_yield_nitri', '', 'N2O yield from nitrification', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_yield_denit, &
      'n2o_yield_denit', '', 'N2O yield from denitrification', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_consumption, &
      'n2o_consumption', 'umol m-3 d-1', 'N2O reduction to N2', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_nitri_rate_used, &
      'nitri_rate_used', 'mmolN m-3 d-1', &
      'effective nitrification N rate used by NOPE', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit_rate_used, &
      'denit_rate_used', 'mmolN m-3 d-1', &
      'effective denitrification N rate used by NOPE', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit_o2_activation, &
      'denit_o2_activation', '', &
      'O2 activation factor for denitrification N2O yield', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit_o2_inhibition, &
      'denit_o2_inhibition', '', &
      'O2 inhibition factor for denitrification N2O yield', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_reduction_o2_factor, &
      'n2o_reduction_o2_factor', '', &
      'O2 factor for N2O reduction', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_reduction_no3_factor, &
      'n2o_reduction_no3_factor', '', &
      'NO3 factor for N2O reduction', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_c_lim, &
      'c_lim', '', 'electron-donor (labile-C) limitation factor for reduction', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_net_production, &
      'n2o_net_production', 'umol m-3 d-1', 'net N2O production (source - sink)', &
      output=output_instantaneous)

   call self%register_diagnostic_variable(self%id_n2o_emit_w, &
      'n2o_emit_w', 'umol m-2 d-1', &
      'airward N2O flux; negative values indicate atmospheric invasion', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_eq, &
      'n2o_eq', 'umol m-3', 'N2O equilibrium with atmosphere', &
      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_n2o_sat, &
      'n2o_sat', '%', 'N2O saturation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Sc, &
      'Sc', '', 'Schmidt number', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_k_gas, &
      'k_gas', 'm d-1', 'gas transfer velocity', output=output_instantaneous)

   call self%register_diagnostic_variable(self%id_n2o_benthic, &
      'n2o_benthic', 'umol m-2 d-1', 'benthic N2O efflux', output=output_instantaneous)

   ! === Dependencies ===
   call self%register_dependency(self%id_denit, 'denit', 'mmol m**-3 d**-1', 'denitrification rate')
   call self%register_dependency(self%id_nitri, 'nitri', 'mmol m**-3 d**-1', 'nitrification rate')
   call self%register_dependency(self%id_oxy,   'oxy',   'mmolO2 m-3',       'dissolved oxygen')
   if (self%enable_n2o_reduction) then
      call self%register_dependency(self%id_no3, 'no3', 'mmolN m-3', &
         'dissolved nitrate')
      ! Couples automatically to a host nitrate with this standard variable;
      ! fabm.yaml can override by name (here oxypom/NO3). O2, nitri and denit
      ! have no FABM standard variable and are coupled by name only.
      call self%request_coupling(self%id_no3, standard_variables%mole_concentration_of_nitrate)
      if (self%enable_c_limit_red) then
         call self%register_dependency(self%id_labile_c, 'labile_c', &
            'mmolC m-3', &
            'labile organic carbon proxy for N2O reduction')
      end if
   end if
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salinity, standard_variables%practical_salinity)
   if (self%wind_forcing_mode == 1) then
      call self%register_dependency(self%id_wind, standard_variables%wind_speed)
   elseif (self%wind_forcing_mode == 2) then
      call self%register_dependency(self%id_day_of_year, &
         standard_variables%number_of_days_since_start_of_the_year)
   end if
   if (self%use_surface_air_pressure) then
      call self%register_dependency(self%id_air_pressure, &
         standard_variables%surface_air_pressure)
   end if

   call self%log_message('NOPE version '//nope_version)

   return

   end subroutine initialize


!-----------------------------------------------------------------------
! Water column: N2O production from nitrification and denitrification and
! reduction to N2.
!-----------------------------------------------------------------------
   subroutine do(self, _ARGUMENTS_DO_)

   class (type_hereon_nope), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk) :: oxy, temp, no3, labile_c
   real(rk) :: nitri_input, denit_input
   real(rk) :: nitri_rate_n, denit_rate_n
   real(rk) :: n2o_w
   real(rk) :: n2o_from_nitri, n2o_from_denit
   real(rk) :: n2o_nitri_nh2oh, n2o_nitri_nd
   real(rk) :: y_nitri, y_denit, y_nd
   real(rk) :: n2o_consumption, n2o_net
   real(rk) :: o2_inh, no3_inh, c_lim, tfac_red, tfac_prod
   real(rk) :: oxy_p, oxy_half_p, oxy_act_p, oxy_inh_p
   real(rk) :: act_denit, inh_denit, inh_at_peak
   real(rk) :: denit_shape, denit_shape_peak
   real(rk) :: oxy_peak, oxy_peak_p
   real(rk) :: sigmoid_arg

   _LOOP_BEGIN_

   ! --- Retrieve dependencies ---
   _GET_(self%id_nitri, nitri_input)
   _GET_(self%id_denit, denit_input)
   _GET_(self%id_oxy, oxy)
   _GET_(self%id_temp, temp)
   _GET_(self%id_n2o_w, n2o_w)

   oxy         = max(oxy, 0.0_rk)
   nitri_input = max(nitri_input, 0.0_rk)
   denit_input = max(denit_input, 0.0_rk)
   n2o_w       = max(n2o_w, 0.0_rk)

   nitri_rate_n = nitri_input * self%nitri_input_to_mmol_n
   denit_rate_n = denit_input * self%denit_input_to_mmol_n

   no3 = 0.0_rk
   labile_c = 0.0_rk
   c_lim = 1.0_rk
   if (self%enable_n2o_reduction) then
      _GET_(self%id_no3, no3)
      no3 = max(no3, 0.0_rk)
      if (self%enable_c_limit_red) then
         _GET_(self%id_labile_c, labile_c)
         labile_c = max(labile_c, 0.0_rk)
         c_lim = labile_c &
            / max(self%k_c_red + labile_c, 1.0e-12_rk)
      end if
   end if

   ! Optional Q10 on production; 1 by default, as host rates include temperature.
   tfac_prod = 1.0_rk
   if (self%apply_prod_temperature_scaling) then
      tfac_prod = self%q10_prod ** ((temp - self%t_ref_prod) / 10.0_rk)
   end if

   ! --- Nitrification N2O, two pathways ---
   !   (a) hydroxylamine oxidation:   y_nh2oh (O2-independent)
   !   (b) nitrifier denitrification: y_nd = y_nd_max * K^p/(O2^p+K^p)
   ! oxy >= 0 and p > 0, so 0**p is defined. Very large p and O2 can overflow
   ! real(rk); exponents up to 8 were tested.
   oxy_p      = oxy**self%oxy_pow_y_nitri
   oxy_half_p = self%oxy_half_y_nitri**self%oxy_pow_y_nitri
   y_nd    = self%y_nd_max * oxy_half_p / (oxy_p + oxy_half_p)
   y_nitri = self%y_nh2oh + y_nd

   ! mmol N m-3 d-1 * yield -> umol N2O m-3 d-1
   n2o_nitri_nh2oh = nitri_rate_n * self%y_nh2oh * tfac_prod &
      * mmol_n_to_umol_n2o
   n2o_nitri_nd = nitri_rate_n * y_nd * tfac_prod &
      * mmol_n_to_umol_n2o
   n2o_from_nitri = n2o_nitri_nh2oh + n2o_nitri_nd

   ! --- Denitrification N2O yield, peaked "leaky-pipe" ---
   ! Hole-in-the-pipe concept (Firestone and Davidson 1989): N2O:N2 is
   ! non-monotonic in O2.
   !   act = O2^p / (O2^p + Kact^p)    low at very low O2 (runs on to N2)
   !   inh = Kinh^p / (Kinh^p + O2^p)  extra yield decreases at high O2
   ! The yield peaks at O2* = sqrt(Kact*Kinh) (69 mmol m-3 for Kact = 40,
   ! Kinh = 120); the peak position follows from the chosen Kact and Kinh,
   ! not from mechanism. O2 is bulk-water O2, so a non-zero yield at oxic O2
   ! is an effective representation of low-O2 microsites in particles.
   oxy_p     = oxy**self%oxy_pow_y_denit
   oxy_act_p = self%oxy_act_y_denit**self%oxy_pow_y_denit
   oxy_inh_p = self%oxy_half_y_denit**self%oxy_pow_y_denit

   act_denit = oxy_p / (oxy_p + oxy_act_p)
   inh_denit = oxy_inh_p / (oxy_inh_p + oxy_p)

   oxy_peak = sqrt(self%oxy_act_y_denit * self%oxy_half_y_denit)
   oxy_peak_p = oxy_peak**self%oxy_pow_y_denit
   inh_at_peak = oxy_inh_p / (oxy_inh_p + oxy_peak_p)
   denit_shape_peak = &
      (oxy_peak_p / (oxy_peak_p + oxy_act_p)) * inh_at_peak

   denit_shape = act_denit * inh_denit &
      / max(denit_shape_peak, 1.0e-12_rk)

   ! denit_shape(O2*) = 1, so both modes give y_denit_peak at O2*.
   !   mode 0: y(O2->0) = y_denit_floor,  y(O2->inf) = y_denit_floor
   !   mode 1: y(O2->0) = y_denit_floor,  y(O2->inf) = 0
   select case (self%denit_floor_mode)
   case (0)
      y_denit = self%y_denit_floor &
         + (self%y_denit_peak - self%y_denit_floor) * denit_shape
   case (1)
      y_denit = self%y_denit_floor * inh_denit &
         + (self%y_denit_peak - self%y_denit_floor * inh_at_peak) * denit_shape
   end select

   n2o_from_denit = denit_rate_n * y_denit * tfac_prod &
      * mmol_n_to_umol_n2o

   ! --- N2O reduction to N2 by denitrifiers ---
   ! consumption = k_max * o2_inh * no3_inh * c_lim * Q10(T) * [N2O]
   ! (a) O2 inhibition (Wicht 1996): 1/(1 + exp((O2 - o2_mid)/o2_width)),
   !     about 1 under anoxia and 0 when oxic.
   ! (b) NO3 competition K/(K + NO3) (assumed form): about 1 when NO3 is
   !     depleted, 0 when abundant. It is independent of the host
   !     denitrification rate, which often falls to zero as NO3 runs out.
   ! (c) Q10 temperature scaling.
   ! (d) Electron-donor limitation c_lim (optional, computed above).
   n2o_consumption = 0.0_rk
   o2_inh = 0.0_rk
   no3_inh = 0.0_rk

   if (self%enable_n2o_reduction) then

      ! (a) O2 inhibition; argument limited to avoid overflow in exp
      sigmoid_arg = (oxy - self%o2_red_mid) / self%o2_red_width
      sigmoid_arg = min(max(sigmoid_arg, -60.0_rk), 60.0_rk)
      o2_inh = 1.0_rk / (1.0_rk + exp(sigmoid_arg))

      ! (b) NO3 competition (assumed form)
      no3_inh = self%kin_no3_n2ored &
         / max(self%kin_no3_n2ored + no3, 1.0e-12_rk)

      ! (c) Q10 temperature scaling
      tfac_red = self%q10_red ** ((temp - self%t_ref_red) / 10.0_rk)

      ! [umol m-3 d-1]
      n2o_consumption = self%k_red_d_max * o2_inh * no3_inh * c_lim * tfac_red * n2o_w

   end if

   ! --- Net N2O source [umol m-3 d-1], per second for FABM ---
   n2o_net = (n2o_from_denit + n2o_from_nitri) - n2o_consumption

   _ADD_SOURCE_(self%id_n2o_w, n2o_net * per_day_to_per_second)

   ! --- Diagnostics ---
   _SET_DIAGNOSTIC_(self%id_n2o_from_nitri, n2o_from_nitri)
   _SET_DIAGNOSTIC_(self%id_n2o_nitri_nh2oh, n2o_nitri_nh2oh)
   _SET_DIAGNOSTIC_(self%id_n2o_nitri_nd, n2o_nitri_nd)
   _SET_DIAGNOSTIC_(self%id_n2o_from_denit, n2o_from_denit)
   _SET_DIAGNOSTIC_(self%id_n2o_yield_nitri, y_nitri)
   _SET_DIAGNOSTIC_(self%id_n2o_yield_denit, y_denit)
   _SET_DIAGNOSTIC_(self%id_n2o_consumption, n2o_consumption)
   _SET_DIAGNOSTIC_(self%id_c_lim, c_lim)
   _SET_DIAGNOSTIC_(self%id_n2o_net_production, n2o_net)
   _SET_DIAGNOSTIC_(self%id_nitri_rate_used, nitri_rate_n)
   _SET_DIAGNOSTIC_(self%id_denit_rate_used, denit_rate_n)
   _SET_DIAGNOSTIC_(self%id_denit_o2_activation, act_denit)
   _SET_DIAGNOSTIC_(self%id_denit_o2_inhibition, inh_denit)
   _SET_DIAGNOSTIC_(self%id_n2o_reduction_o2_factor, o2_inh)
   _SET_DIAGNOSTIC_(self%id_n2o_reduction_no3_factor, no3_inh)

   _LOOP_END_

   end subroutine do


!-----------------------------------------------------------------------
! Surface: air-water N2O exchange from N2O solubility, Schmidt number and
! wind-based gas transfer velocity (positive flux = to the atmosphere).
!-----------------------------------------------------------------------
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)

   class (type_hereon_nope), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: n2o_conc_w
   real(rk) :: n2o_eq, n2o_sat
   real(rk) :: n2o_sea_air_flux
   real(rk) :: temp, temp_sc, temp_abs, salinity
   real(rk) :: wind_speed
   real(rk) :: Sc, k_gas
   real(rk) :: day_of_year
   real(rk), parameter :: pa_per_atm = 101325.0_rk
   real(rk) :: air_pressure_pa
   real(rk) :: p_total_atm, p_h2o_atm, p_dry_atm

   _SURFACE_LOOP_BEGIN_

      _GET_(self%id_n2o_w, n2o_conc_w)
      _GET_(self%id_temp, temp)
      _GET_(self%id_salinity, salinity)

      ! Wind speed: forcing (1), seasonal means by day of year (2), annual mean (3)
      select case (self%wind_forcing_mode)
      case (1)
         _GET_SURFACE_(self%id_wind, wind_speed)
      case (2)
         _GET_GLOBAL_(self%id_day_of_year, day_of_year)
         if (day_of_year >= 335._rk .or. day_of_year <= 59._rk) then
            wind_speed = self%wind_mean_winter
         elseif (day_of_year <= 151._rk) then
            wind_speed = self%wind_mean_spring
         elseif (day_of_year <= 243._rk) then
            wind_speed = self%wind_mean_summer
         else
            wind_speed = self%wind_mean_fall
         end if
      case (3)
         wind_speed = self%wind_mean_annual
      end select
      wind_speed = max(wind_speed, 0.0_rk)

      ! Limit temperature to the 0-40 degC range of the solubility and Schmidt
      ! polynomials (Weiss and Price 1980; Wanninkhof 2014); extrapolation can
      ! give a negative Schmidt number.
      temp_sc  = min(max(temp, 0.0_rk), 40.0_rk)
      salinity = max(salinity, 0.0_rk)
      temp_abs = 273.15_rk + temp_sc

      ! Atmospheric pressure [atm]; optional water-vapour correction to dry-air
      ! pressure (Weiss and Price 1980)
      if (self%use_surface_air_pressure) then
         _GET_SURFACE_(self%id_air_pressure, air_pressure_pa)
         p_total_atm = max(air_pressure_pa / pa_per_atm, 0.0_rk)
      else
         p_total_atm = self%patm_atm
      end if

      p_dry_atm = p_total_atm
      if (self%use_h2o_correction) then
         p_h2o_atm = exp(24.4543_rk &
            - 67.4509_rk * (100.0_rk / temp_abs) &
            - 4.8489_rk * log(temp_abs / 100.0_rk) &
            - 0.000544_rk * salinity)
         p_dry_atm = max(p_total_atm - p_h2o_atm, 0.0_rk)
      end if

      ! N2O equilibrium concentration [umol m-3] (Weiss and Price 1980 solubility)
      n2o_eq = self%n2o_atm_ppb * p_dry_atm * exp( &
         -165.8806_rk &
         + 222.8743_rk * (100.0_rk / temp_abs) &
         + 92.0792_rk * log(temp_abs / 100.0_rk) &
         - 1.48425_rk * (temp_abs / 100.0_rk)**2 &
         + salinity * ( &
            -0.056235_rk &
            + 0.0316119_rk * (temp_abs / 100.0_rk) &
            - 0.004872_rk * (temp_abs / 100.0_rk)**2))

      ! N2O Schmidt number (Wanninkhof 2014, Table 1): freshwater Sc(20 degC) = 626,
      ! seawater Sc(20 degC) = 697. At Seemannshoeft (2.3-21.9 degC in 2015) the
      ! seawater set would lower k_gas by about 5 %.
      select case (self%schmidt_water_type)
      case (1)
         ! freshwater N2O
         Sc = 2141.2_rk &
            - 152.56_rk * temp_sc &
            + 5.8963_rk * temp_sc**2 &
            - 0.12411_rk * temp_sc**3 &
            + 0.0010655_rk * temp_sc**4
      case (2)
         ! seawater N2O
         Sc = 2356.2_rk &
            - 166.38_rk * temp_sc &
            + 6.3952_rk * temp_sc**2 &
            - 0.13422_rk * temp_sc**3 &
            + 0.0011506_rk * temp_sc**4
      end select

      ! Gas transfer velocity [m d-1] (0.24 converts cm h-1 to m d-1)
      k_gas = 0.0_rk  ! silences a compiler warning; scheme is limited to 1-3
      select case (self%gas_transfer_scheme)
      case (1)
         ! Wanninkhof (2014)
         k_gas = 0.24_rk * 0.251_rk * wind_speed**2 * (Sc / 660.0_rk)**(-0.5_rk)
      case (2)
         ! Nightingale et al. (2000)
         k_gas = 0.24_rk * (0.222_rk * wind_speed**2 + 0.333_rk * wind_speed) &
                * (Sc / 600.0_rk)**(-0.5_rk)
      case (3)
         ! Wanninkhof (1992): cubic above 11 m s-1, quadratic from 3 to 11 m s-1;
         ! Raymond and Cole (2001) exponential below 3 m s-1; as OxyPOM O2 reaeration.
         if (wind_speed > 11.0_rk) then
            k_gas = 0.0283_rk * wind_speed**3
         else if (wind_speed > 3.0_rk) then
            k_gas = 0.312_rk * wind_speed**2
         else
            k_gas = 0.9826_rk * exp(0.35_rk * wind_speed)
         end if
         k_gas = 0.24_rk * self%gas_transfer_factor * k_gas * (Sc / 660.0_rk)**(-0.5_rk)
      end select

      ! Sea-to-air flux [umol m-2 d-1] and saturation [%]
      n2o_sea_air_flux = (n2o_conc_w - n2o_eq) * k_gas
      n2o_sat = (n2o_conc_w / max(n2o_eq, 1.0e-10_rk)) * 100.0_rk

      _ADD_SURFACE_FLUX_(self%id_n2o_w, -n2o_sea_air_flux * per_day_to_per_second)

      _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_emit_w, n2o_sea_air_flux)
      _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_eq, n2o_eq)
      _SET_SURFACE_DIAGNOSTIC_(self%id_n2o_sat, n2o_sat)
      _SET_SURFACE_DIAGNOSTIC_(self%id_Sc, Sc)
      _SET_SURFACE_DIAGNOSTIC_(self%id_k_gas, k_gas)

   _SURFACE_LOOP_END_

   end subroutine do_surface


!-----------------------------------------------------------------------
! Bottom: optional empirical sediment N2O efflux, linear in temperature and
! enhanced at low bottom O2. For sensitivity tests only; off by default.
!-----------------------------------------------------------------------
   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)

   class (type_hereon_nope), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   real(rk) :: temp_b, oxy_b
   real(rk) :: benthic_flux, o2_enh

   _BOTTOM_LOOP_BEGIN_

      benthic_flux = 0.0_rk

      if (self%enable_benthic_n2o) then
         _GET_(self%id_temp, temp_b)
         _GET_(self%id_oxy,  oxy_b)
         oxy_b = max(oxy_b, 0.0_rk)

         ! Low-O2 enhancement: factor rises from about 1 (oxic) to 2 (anoxic)
         o2_enh = 1.0_rk + self%benthic_o2_half &
            / max(self%benthic_o2_half + oxy_b, tiny(1.0_rk))

         ! Temperature-dependent efflux [umol m-2 d-1], floored at zero
         benthic_flux = (self%benthic_n2o_base + self%benthic_n2o_tcoef * temp_b) * o2_enh
         benthic_flux = max(benthic_flux, 0.0_rk)

         _ADD_BOTTOM_FLUX_(self%id_n2o_w, benthic_flux * per_day_to_per_second)
      end if

      _SET_BOTTOM_DIAGNOSTIC_(self%id_n2o_benthic, benthic_flux)

   _BOTTOM_LOOP_END_

   end subroutine do_bottom

   end module hereon_nope
