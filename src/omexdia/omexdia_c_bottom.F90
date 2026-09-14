#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hereon_omexdia_c_bottom --- Adaptive-trait model for sediment
!          penetration depths (oxygen + fast/slow detritus + methane),
!          benthic, with full electron-acceptor partitioning
!
! !INTERFACE:
   module hereon_omexdia_c_bottom
!
! !DESCRIPTION:
!
! Trait-based reformulation of parts of hereon_omexdia_c, following Kai
! Wirtz's adaptive-dynamics ("effective phenotype") framework (Wirtz &
! Eckhardt 1996 and successors). Depth is treated as a set of continuous,
! dynamically adapting traits rather than discrete model layers. Each
! trait is a single characteristic value per sediment column, 
! so this is a bottom 2D model: it implements do_bottom, not do,
! and the trait state variables are horizontal (bottom) state variables
! (type_bottom_state_variable_id) -- domain_bottom is implicit in that
! type and not a separate argument to register_state_variable.
!
! THREE TRAITS:
!   id_oxy_depth  = L,    typical O2 penetration depth (the depth at
!                         which O2 has dropped to 1/e of its surface
!                         value, assuming O2(z,t)=O2_surface(t)*exp(-z/L)).
!   id_fdet_depth = L_Cf, typical fast-detritus penetration depth.
!   id_sdet_depth = L_Cs, typical slow-detritus penetration depth.
!
! FOUR PELAGIC CONCENTRATIONS, read-only (id_oxy, id_fdet, id_sdet,
! id_ch4) -- NOT state dependencies. This module never declares or
! touches any interior/"3D" state variable itself: it only reads plain
! named horizontal dependencies (oxy_c0/fdet_c0/sdet_c0/ch4_c0, same
! convention as no3_c0/so4_c0/temp_c0 below), and its mass-balance
! consumption/production of each is exposed as a horizontal DIAGNOSTIC
! (oxy_flux/fdet_flux/sdet_flux/ch4_flux, mmol m-2 s-1, positive into
! the pelagic pool). Applying that diagnostic as an actual flux onto a
! real pelagic state variable is delegated to FABM's own builtin
! external_bottom_flux utility model, wired up per deployment in YAML
! (one external_bottom_flux instance per exchanged quantity, coupling
! its 'target' to the real state variable and its 'flux' to this
! module's diagnostic) -- see fabm-omexdia-c-bottom-coupled.yaml. This
! keeps the module itself host-agnostic: it does not need to know at
! compile time what (if anything) actually holds these pools.
!
! FOUR READ-ONLY FORCING DEPENDENCIES (no flux returned by this module
! for any of these -- all plain named horizontal dependencies, no
! standard_variable, for the same host-agnostic reason as above):
!   id_no3  -- interfacial nitrate (used only to compute the partition
!              fractions below; NO3 stays diagnostic, not prognostic --
!              see rationale in the accompanying conversation: NO3's own
!              kinetics don't cleanly support the same "thin subsurface
!              zone" simplification that justified treating nh3/ODU with
!              a quasi-steady shortcut, so it is not given a flux/budget
!              of its own here).
!   id_so4  -- interfacial sulfate. NOTE: earlier versions of this model
!              derived so4 from salinity via hereon_omexdia_c's own
!              sulfate_from_salinity() proxy. Real station data (Elbe
!              estuary: Geesthacht, salinity=0, observed so4=800 mmol/m3)
!              showed this proxy fails badly in freshwater reaches with
!              non-marine (geological/mining-influenced) sulfate sources.
!              so4 is therefore a REQUIRED, INDEPENDENT dependency here --
!              couple it to real observations/forcing, or to a separate
!              diagnostic module implementing the salinity proxy, as
!              appropriate for the deployment site.
!   id_temp -- bottom water temperature, used only for the methanogenesis
!              Q10 factor. NOTE: NOT standard_variables%temperature --
!              this FABM tree has no bottom-domain temperature standard
!              variable (only bottom_depth/bottom_roughness_length/
!              bottom_stress/bottom_depth_below_geoid exist as bottom
!              standard variables), and the interior temperature standard
!              variable would pull in exactly the "3D" coupling this
!              module otherwise avoids -- so temp_c0 is a plain named
!              dependency like no3_c0/so4_c0, supplied however the
!              deployment wants (a real near-bed extraction, a constant,
!              anything).
!
! FITNESS-GRADIENT ODE FOR L (oxygen depth)
! -------------------------------------------
!   benefit(L)   = O2_surface * exp(-1) / L
!   dCost/dL     = alpha*ox_frac*(marg_fdet+marg_sdet)          [OxicMin]
!                + 2*ox_frac*(NCrFdet*marg_fdet+NCrSdet*marg_sdet)  [Nitri]
!                + su_frac*(marg_fdet+marg_sdet)                [OduOx]
!   dL/dt        = rate_L * ( benefit(L) - dCost/dL )
!
! IMPORTANT: dCost/dL uses the TOTAL local O2 sink (OxicMin + Nitri +
! OduOx), not just OxicMin's own share. An earlier version of this
! module drove the trait using only OxicMin's marginal cost; since Nitri
! alone draws ~40-50% on top of OxicMin at typical marine conditions,
! that version substantially understated the true local O2 demand and
! pushed the equilibrium too deep (a correction that shifted the
! validated baseline equilibrium from ~0.049 m to ~0.007 m -- confirmed
! via direct root-finding, not just integration, before trusting it).
!
! ox_frac, de_frac, su_frac, me_frac are the OxicMin/Denitrific/
! SulfateMin/Methanolim shares of total carbon mineralization, computed
! from hereon_omexdia_c's own Rescale-cascade kinetics evaluated at
! interfacial conditions (a "well-mixed zone" simplification: assumes the
! whole 0-to-L zone shares the same partition the real kinetics predict
! at the interface, avoiding separate depth-resolving traits/profiles for
! NO3/SO4/CH4, each of which would risk its own degenerate-coupling traps
! -- see the LC_fast/LC_slow caveat below for a concrete example of one
! such trap that WAS found and fixed).
!
! ESCAPE-THRESHOLD CAVEAT (structural, not a bug): benefit(L)~1/L is
! algebraic while dCost/dL decays exponentially in L. This pairing
! generically produces a stable equilibrium L1 AND an unstable escape
! threshold L2 above which benefit permanently exceeds cost and L grows
! unboundedly. Verified numerically across many scenarios; confirmed to
! occur for real station data (some real Elbe station parameter
! combinations tested during development had NO interior equilibrium at
! all -- benefit exceeded cost at every depth tested). Lmax below is a
! pragmatic domain ceiling for exactly this situation: dL/dt is frozen at
! Lmax rather than allowing the state to diverge to physically
! meaningless values. This does not fix the underlying absence of an
! interior equilibrium for such parameter regimes -- it only bounds the
! symptom. A diagnostic (id_at_ceiling) flags when this is happening so
! it is never silently mistaken for a genuine balance point.
!
! LC_fast/LC_slow COUPLING TO L (saturating, NOT proportional)
! ----------------------------------------------------------------
!   LCf_target(L) = LCfast_max * L / (L + Lhalf_fast)
!   LCs_target(L) = LCslow_max * L / (L + Lhalf_slow)
!   dLCf/dt       = ( LCf_target(L) - LCf ) / tau_fast
!   dLCs/dt       = ( LCs_target(L) - LCs ) / tau_slow
! CAVEAT: making LCf/LCs exactly PROPORTIONAL to L (LC=k*L) is
! mathematically degenerate -- it forces L/LC constant, collapsing
! dCost/dL into the same 1/L shape as benefit(L), making the difference
! either identically zero everywhere or never zero. Caught by testing
! numerically before being implemented; do not simplify back to
! proportional coupling.
!
! CH4: given a REAL dynamic pelagic-partner coupling (not a quasi-steady
! shortcut like nh3/ODU), because methanogenesis happens below the oxic
! zone by construction and CH4 must physically escape upward before
! aerobic oxidation can touch it -- a quasi-steady assumption would mean
! CH4 never accumulates or escapes, contradicting the entire premise that
! it is observed (this module was extended specifically because of
! methane observations in freshwater/estuarine deployments, where
! so4=0 removes anaerobic oxidation of methane (AOM) entirely and
! methanogenesis can reach the majority of carbon mineralization under
! simultaneous O2/NO3 depletion -- confirmed up to ~88% of Cprod in
! freshwater+hypoxic+NO3-exhausted test conditions, vs <1% at full
! marine salinity for the same O2/NO3 state).
!
! UNIT NOTE on the CH4 flux: ch4 production (from carbon demand,
! integrated 0 to L) is naturally AREAL (mmol m-2 d-1); aerobic oxidation
! and AOM are naturally VOLUMETRIC rates (mmol m-3 d-1) evaluated at the
! bottom-most pelagic cell's own CH4 concentration. H_REF is a placeholder
! "near-bed mixed layer thickness" reconciling the two -- inherited
! directly from Python prototyping, where it stood in for real vertical
! transport. Revisit if the host model exposes its own bottom-cell
! thickness as a usable dependency.
!
! !USES:
   use fabm_types

   implicit none

   private
   public type_hereon_omexdia_c_bottom
!
   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: euler     = 2.718281828459045_rk
!
! !REVISION HISTORY:
!  Original author(s): Carsten Lemmen (trait-based reformulation, building
!  on hereon_omexdia_c by Richard Hofmeister & Kai Wirtz)
!  Consolidated benthic version: fdet/sdet-driven cost, dynamic detritus
!  penetration depths, full electron-acceptor partitioning (OxicMin/
!  Nitri/OduOx folded into a consistent total O2 budget; NO3 diagnostic;
!  CH4 as a real dynamic pelagic partner; SO4 as an independent
!  dependency, not derived from salinity; Lmax domain ceiling. Ported
!  from Python prototyping and station-scenario testing (Elbe estuary:
!  Geesthacht, Seemannshoft, Brunsbuttel, Cuxhaven).
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hereon_omexdia_c_bottom

      type (type_bottom_state_variable_id)          :: id_oxy_depth
      type (type_bottom_state_variable_id)          :: id_fdet_depth
      type (type_bottom_state_variable_id)          :: id_sdet_depth

!     Pelagic concentrations, read-only (plain named dependencies -- see
!     header; NOT state dependencies, no flux returned via these ids).
      type (type_dependency_id)                      :: id_oxy
      type (type_dependency_id)                      :: id_fdet
      type (type_dependency_id)                      :: id_sdet
      type (type_dependency_id)                      :: id_ch4

!     Read-only forcing dependencies (no flux returned).
      type (type_dependency_id)                      :: id_no3
      type (type_dependency_id)                      :: id_so4
      type (type_dependency_id)                      :: id_temp
!     Diagnostics (horizontal-only, matching the trait states' domain).
      type (type_horizontal_diagnostic_variable_id) :: id_dLdt
      type (type_horizontal_diagnostic_variable_id) :: id_ox_frac, id_de_frac
      type (type_horizontal_diagnostic_variable_id) :: id_su_frac, id_me_frac
      type (type_horizontal_diagnostic_variable_id) :: id_net_no3, id_odu_prod
      type (type_horizontal_diagnostic_variable_id) :: id_ch4_prod, id_at_ceiling
!     Mass-balance flux OUTPUTS as diagnostics (mmol m-2 s-1, positive
!     into the pelagic pool) -- applied to a real state variable by a
!     separate external_bottom_flux instance in YAML, see header.
      type (type_horizontal_diagnostic_variable_id) :: id_oxy_flux, id_fdet_flux
      type (type_horizontal_diagnostic_variable_id) :: id_sdet_flux, id_ch4_flux

!     Model parameters
      real(rk) :: rate_L, alpha, rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: LCfast_max, Lhalf_fast, tau_fast
      real(rk) :: LCslow_max, Lhalf_slow, tau_slow
      real(rk) :: Lmin, LCmin, Lmax
      real(rk) :: ksO2oxic, kinO2denit, ksNO3denit
      real(rk) :: kinNO3anox, kinO2anox, kinSO4, relaxO2
      real(rk) :: q10_meth, Tref
      real(rk) :: nh3_amb, odu_amb
      real(rk) :: rmaxO2, ksCH4, ksAOM, H_REF

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do_bottom

   end type type_hereon_omexdia_c_bottom
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the trait-depth model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
   class (type_hereon_omexdia_c_bottom),intent(inout),target :: self
   integer,                           intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC

   ! --- trait adaptation & cost-term parameters ---
   call self%get_parameter(self%rate_L,'rate_L','m5 mmolO2-1 d-1', &
        'oxygen-depth trait adaptation rate (= m/e)', default=1.0e-3_rk)
   call self%get_parameter(self%alpha,'alpha','-', &
        'O2:C stoichiometric ratio for aerobic mineralization', default=0.96_rk)
   call self%get_parameter(self%rFast,'rFast','d-1', &
        'aerobic-equivalent reactivity, fast detritus (match hereon_omexdia_c rFast)', &
        default=0.08_rk)
   call self%get_parameter(self%rSlow,'rSlow','d-1', &
        'aerobic-equivalent reactivity, slow detritus (match hereon_omexdia_c rSlow)', &
        default=0.0005_rk)
   call self%get_parameter(self%NCrFdet,'NCrFdet','molN molC-1', &
        'N:C ratio, fast detritus (match hereon_omexdia_c)', default=0.2_rk)
   call self%get_parameter(self%NCrSdet,'NCrSdet','molN molC-1', &
        'N:C ratio, slow detritus (match hereon_omexdia_c)', default=0.01_rk)

   ! --- detritus penetration-depth coupling (saturating, see header) ---
   call self%get_parameter(self%LCfast_max,'LCfast_max','m', &
        'saturating asymptote for fdet penetration depth', default=0.2_rk)
   call self%get_parameter(self%Lhalf_fast,'Lhalf_fast','m', &
        'half-saturation depth for fdet penetration-depth coupling to L', default=0.0491_rk)
   call self%get_parameter(self%tau_fast,'tau_fast','d', &
        'relaxation timescale, fdet penetration depth', default=30.0_rk)
   call self%get_parameter(self%LCslow_max,'LCslow_max','m', &
        'saturating asymptote for sdet penetration depth', default=0.6_rk)
   call self%get_parameter(self%Lhalf_slow,'Lhalf_slow','m', &
        'half-saturation depth for sdet penetration-depth coupling to L', default=0.0491_rk)
   call self%get_parameter(self%tau_slow,'tau_slow','d', &
        'relaxation timescale, sdet penetration depth', default=90.0_rk)

   ! --- numerical floors/ceiling ---
   call self%get_parameter(self%Lmin,'Lmin','m', &
        'minimum oxygen penetration depth (numerical floor)', default=1.0e-3_rk)
   call self%get_parameter(self%LCmin,'LCmin','m', &
        'minimum detritus penetration depth (numerical floor)', default=1.0e-3_rk)
   call self%get_parameter(self%Lmax,'Lmax','m', &
        'maximum oxygen penetration depth (domain ceiling; dL/dt frozen '// &
        'above this -- see header caveat on the escape threshold)', default=0.5_rk)

   ! --- electron-acceptor partition kinetics (match hereon_omexdia_c) ---
   call self%get_parameter(self%ksO2oxic,'ksO2oxic','mmolO2 m-3', &
        'half-saturation O2 in oxic mineralization', default=3.0_rk)
   call self%get_parameter(self%kinO2denit,'kinO2denit','mmolO2 m-3', &
        'half-saturation O2 inhibition of denitrification', default=70.0_rk)
   call self%get_parameter(self%ksNO3denit,'ksNO3denit','mmolNO3 m-3', &
        'half-saturation NO3 in denitrification', default=1.0_rk)
   call self%get_parameter(self%kinNO3anox,'kinNO3anox','mmolNO3 m-3', &
        'half-saturation NO3 inhibition of anoxic mineralization', default=1.0_rk)
   call self%get_parameter(self%kinO2anox,'kinO2anox','mmolO2 m-3', &
        'half-saturation O2 inhibition of anoxic mineralization', default=1.0_rk)
   call self%get_parameter(self%kinSO4,'kinSO4','mmol m-3', &
        'sulfate inhibition threshold for methanogenesis', default=1000.0_rk)
   call self%get_parameter(self%relaxO2,'relaxO2','-', &
        'relaxation term in OxicMin denominator', default=0.04_rk)
   call self%get_parameter(self%q10_meth,'q10_meth','-', &
        'Q10 scaling coefficient for methanogenesis', default=3.5_rk)
   call self%get_parameter(self%Tref,'Tref','K', &
        'reference temperature for Q10 scaling', default=288.15_rk)
   call self%get_parameter(self%nh3_amb,'nh3_amb','mmol N m-3', &
        'ambient nh3 concentration used in OxicMin denominator (fixed, not dynamic -- '// &
        'never varied in prototype testing)', default=40.0_rk)
   call self%get_parameter(self%odu_amb,'odu_amb','mmol m-3', &
        'ambient ODU concentration used in OxicMin denominator (fixed, not dynamic)', &
        default=100.0_rk)

   ! --- CH4 parameters ---
   call self%get_parameter(self%rmaxO2,'rmaxO2','d-1', &
        'maximum aerobic CH4 oxidation rate', default=10.0_rk)
   call self%get_parameter(self%ksCH4,'ksCH4','mmol m-3', &
        'half-saturation CH4 for aerobic oxidation', default=5.0_rk)
   call self%get_parameter(self%ksAOM,'ksAOM','m3 mmol-1 d-1', &
        'second-order anaerobic CH4 oxidation (AOM) rate', default=0.05_rk)
   call self%get_parameter(self%H_REF,'H_REF','m', &
        'placeholder near-bed mixed-layer thickness reconciling areal CH4 '// &
        'production with volumetric consumption rates (see header note)', default=0.1_rk)

   ! --- register the three traits (bottom, horizontal-only; domain_bottom
   ! is implicit in type_bottom_state_variable_id, not a separate argument) ---
   ! maximum=self%Lmax is the REAL enforcement of the escape-threshold
   ! ceiling (see header caveat): the dL/dt freeze below only stops L
   ! from growing further once it is already >= Lmax, it does not clamp
   ! an overshoot accrued within a single (possibly coarse) timestep --
   ! confirmed by testing: without this maximum=, a 1800s step pushed L
   ! to ~3.47 m, nearly 7x Lmax, before the frozen derivative caught up.
   call self%register_state_variable(self%id_oxy_depth,'oxy_depth','m', &
        'typical O2 penetration depth (adaptive trait)', 0.01_rk, minimum=self%Lmin, &
        maximum=self%Lmax)
   call self%register_state_variable(self%id_fdet_depth,'fdet_depth','m', &
        'typical fdet penetration depth (adaptive trait)', 0.1_rk, minimum=self%LCmin)
   call self%register_state_variable(self%id_sdet_depth,'sdet_depth','m', &
        'typical sdet penetration depth (adaptive trait)', 0.3_rk, minimum=self%LCmin)

   ! --- diagnostics ---
   call self%register_diagnostic_variable(self%id_dLdt,'dLdt','m d-1', &
        'oxygen-depth trait adaptation rate dL/dt', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_ox_frac,'ox_frac','-', &
        'OxicMin share of total carbon mineralization', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_de_frac,'de_frac','-', &
        'Denitrification share of total carbon mineralization', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_su_frac,'su_frac','-', &
        'Sulfate-reduction share of total carbon mineralization', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_me_frac,'me_frac','-', &
        'Methanogenesis share of total carbon mineralization', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_net_no3,'net_no3','mmol N m-2 d-1', &
        'net NO3 production (nitrification - denitrification), diagnostic only', &
        output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_odu_prod,'odu_production','mmol m-2 d-1', &
        'ODU (Mn/Fe/S placeholder) production from sulfate-linked mineralization, '// &
        'diagnostic only', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_ch4_prod,'ch4_production','mmol C m-2 d-1', &
        'CH4 production from methanogenesis-linked mineralization', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_at_ceiling,'at_ceiling','-', &
        '1 if oxy_depth is pinned at Lmax (no interior equilibrium for current '// &
        'forcing -- see header caveat), 0 otherwise', output=output_instantaneous, source=source_do_bottom)

   ! --- pelagic concentrations, read-only (see header: plain named
   ! dependencies, not state dependencies -- this module never touches
   ! any interior/"3D" state variable directly). ---
   call self%register_dependency(self%id_oxy,'oxy_c0','mmol O2 m-3', &
        'oxygen concentration at the sediment-water interface')
   call self%register_dependency(self%id_fdet,'fdet_c0','mmol C m-3', &
        'fast detritus concentration at the sediment-water interface')
   call self%register_dependency(self%id_sdet,'sdet_c0','mmol C m-3', &
        'slow detritus concentration at the sediment-water interface')
   call self%register_dependency(self%id_ch4,'ch4_c0','mmol C m-3', &
        'methane concentration at the sediment-water interface')

   ! --- mass-balance flux OUTPUTS, as diagnostics (mmol m-2 s-1, positive
   ! into the pelagic pool) -- consumed by a separate external_bottom_flux
   ! instance per quantity in YAML to actually apply them, see header. ---
   call self%register_diagnostic_variable(self%id_oxy_flux,'oxy_flux','mmol O2 m-2 s-1', &
        'net O2 exchange with the pelagic oxy pool (positive into the water column)', &
        output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_fdet_flux,'fdet_flux','mmol C m-2 s-1', &
        'net fdet exchange with the pelagic fdet pool (positive into the water column)', &
        output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_sdet_flux,'sdet_flux','mmol C m-2 s-1', &
        'net sdet exchange with the pelagic sdet pool (positive into the water column)', &
        output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_ch4_flux,'ch4_flux','mmol C m-2 s-1', &
        'net ch4 exchange with the pelagic ch4 pool (positive into the water column)', &
        output=output_instantaneous, source=source_do_bottom)

   ! --- read-only forcing dependencies ---
   ! so4 is NOT derived from salinity here (see header rationale) -- it
   ! must be coupled explicitly, either to real observations/forcing or
   ! to a separate module implementing whatever salinity proxy is
   ! appropriate for the deployment site.
   call self%register_dependency(self%id_no3,'no3_c0','mmol N m-3', &
        'nitrate concentration at the sediment-water interface (diagnostic use only)')
   call self%register_dependency(self%id_so4,'so4_c0','mmol m-3', &
        'sulfate concentration at the sediment-water interface (couple explicitly -- '// &
        'do NOT assume a fixed salinity-sulfate relationship, see header note)')
   ! NOT standard_variables%temperature -- no bottom-domain temperature
   ! standard variable exists in this FABM tree, and the interior one
   ! would reintroduce exactly the "3D" coupling this module avoids
   ! elsewhere (see header). Plain named dependency instead, same
   ! convention as no3_c0/so4_c0.
   call self%register_dependency(self%id_temp,'temp_c0','degree_C', &
        'bottom water temperature (methanogenesis Q10 factor only)')

   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand side of the trait-depth model (benthic)
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_hereon_omexdia_c_bottom),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk) :: oxy_surface, fdet_surface, sdet_surface, ch4_surface, no3_surface, so4, temp_celsius
   real(rk) :: L, LCf, LCs, Leff, LCfeff, LCseff, ch4eff
   real(rk) :: temp_kelvin, E_a_meth, f_temp_meth
   real(rk) :: Oxicminlim, Denitrilim, anoxic_space, f_so4, f_meth
   real(rk) :: SulfateMinlim, Methanolim, total_lim
   real(rk) :: ox_frac, de_frac, su_frac, me_frac
   real(rk) :: marg_fdet, marg_sdet, fdet_demand, sdet_demand, total_C
   real(rk) :: marg_oxicmin, marg_nitri, marg_oduox, dCostdL
   real(rk) :: benefit, dLdt
   real(rk) :: LCf_target, LCs_target, dLCfdt, dLCsdt
   real(rk) :: oxicmin_o2, n_release, nitri_o2, odu_production, oduox_o2
   real(rk) :: ch4_production, ch4_aerobic_ox_vol, ch4_aom_vol
   real(rk) :: ch4_aerobic_o2, ch4_consumption_areal, total_o2_demand
   real(rk) :: nitri_n_production, denitri_n_consumption, net_no3
   real(rk) :: at_ceiling_flag

!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   ! Interfacial (bottom-most pelagic cell) concentrations and forcing
   _GET_(self%id_oxy,oxy_surface)
   _GET_(self%id_fdet,fdet_surface)
   _GET_(self%id_sdet,sdet_surface)
   _GET_(self%id_ch4,ch4_surface)
   _GET_(self%id_no3,no3_surface)
   _GET_(self%id_so4,so4)
   _GET_(self%id_temp,temp_celsius)

   ! Current trait values (bottom, horizontal-only state variables)
   _GET_HORIZONTAL_(self%id_oxy_depth,L)
   _GET_HORIZONTAL_(self%id_fdet_depth,LCf)
   _GET_HORIZONTAL_(self%id_sdet_depth,LCs)

   Leff   = max(L,   self%Lmin)
   LCfeff = max(LCf, self%LCmin)
   LCseff = max(LCs, self%LCmin)
   ch4eff = max(ch4_surface, 0.0_rk)

   ! --- electron-acceptor partition fractions (Rescale-cascade kinetics,
   ! evaluated at interfacial conditions; "well-mixed zone" simplification) ---
   temp_kelvin = 273.15_rk + temp_celsius
   E_a_meth = 0.1_rk * log(self%q10_meth) * self%Tref * (self%Tref + 10.0_rk)
   f_temp_meth = exp(-E_a_meth * (1.0_rk/temp_kelvin - 1.0_rk/self%Tref))

   Oxicminlim = oxy_surface / (oxy_surface + self%ksO2oxic + self%relaxO2*(self%nh3_amb + self%odu_amb))
   Denitrilim = (1.0_rk - oxy_surface/(oxy_surface + self%kinO2denit)) * no3_surface/(no3_surface + self%ksNO3denit)
   anoxic_space = (1.0_rk - oxy_surface/(oxy_surface + self%kinO2anox)) * (1.0_rk - no3_surface/(no3_surface + self%kinNO3anox))
   f_so4 = so4 / (self%kinSO4 + so4)
   f_meth = self%kinSO4 / (self%kinSO4 + so4)
   SulfateMinlim = anoxic_space * f_so4
   Methanolim = anoxic_space * f_meth * f_temp_meth

   total_lim = Oxicminlim + Denitrilim + SulfateMinlim + Methanolim
   ox_frac = Oxicminlim / total_lim
   de_frac = Denitrilim / total_lim
   su_frac = SulfateMinlim / total_lim
   me_frac = Methanolim / total_lim

   ! --- marginal & integrated carbon demand ---
   marg_fdet = self%rFast * (fdet_surface/LCfeff) * exp(-Leff/LCfeff)
   marg_sdet = self%rSlow * (sdet_surface/LCseff) * exp(-Leff/LCseff)
   fdet_demand = self%rFast * fdet_surface * LCfeff * (1.0_rk - exp(-Leff/LCfeff))
   sdet_demand = self%rSlow * sdet_surface * LCseff * (1.0_rk - exp(-Leff/LCseff))
   total_C = fdet_demand + sdet_demand

   ! --- oxygen-depth trait ODE: marginal cost uses the TOTAL local O2
   ! sink (OxicMin + Nitri + OduOx), not just OxicMin's own share (see
   ! header rationale) ---
   marg_oxicmin = self%alpha * ox_frac * (marg_fdet + marg_sdet)
   marg_nitri = 2.0_rk * ox_frac * (self%NCrFdet*marg_fdet + self%NCrSdet*marg_sdet)
   marg_oduox = su_frac * (marg_fdet + marg_sdet)
   dCostdL = marg_oxicmin + marg_nitri + marg_oduox

   benefit = oxy_surface / (euler * Leff)
   dLdt = self%rate_L * (benefit - dCostdL)
   ! Reflecting floor and domain ceiling (see header escape-threshold caveat)
   if (L <= self%Lmin .and. dLdt < 0.0_rk) dLdt = 0.0_rk
   if (L >= self%Lmax .and. dLdt > 0.0_rk) dLdt = 0.0_rk
   at_ceiling_flag = merge(1.0_rk, 0.0_rk, L >= self%Lmax*0.999_rk)

   ! --- detritus-depth trait ODEs (saturating relaxation towards L) ---
   LCf_target = self%LCfast_max * Leff / (Leff + self%Lhalf_fast)
   dLCfdt = (LCf_target - LCf) / self%tau_fast
   if (LCf <= self%LCmin .and. dLCfdt < 0.0_rk) dLCfdt = 0.0_rk

   LCs_target = self%LCslow_max * Leff / (Leff + self%Lhalf_slow)
   dLCsdt = (LCs_target - LCs) / self%tau_slow
   if (LCs <= self%LCmin .and. dLCsdt < 0.0_rk) dLCsdt = 0.0_rk

   ! --- total (integrated) O2 demand -> mass-balance flux to pelagic oxy ---
   oxicmin_o2 = self%alpha * ox_frac * total_C
   n_release = self%NCrFdet*fdet_demand + self%NCrSdet*sdet_demand
   nitri_o2 = 2.0_rk * ox_frac * n_release
   odu_production = su_frac * total_C
   oduox_o2 = odu_production   ! 1:1 stoichiometry, as in hereon_omexdia_c

   ! --- CH4 production/consumption: real dynamic pelagic-partner coupling ---
   ch4_production = 0.5_rk * me_frac * total_C                          ! areal, mmol C m-2 d-1
   ch4_aerobic_ox_vol = self%rmaxO2 * (ch4eff/(self%ksCH4+ch4eff)) * ox_frac  ! volumetric, mmol m-3 d-1
   ch4_aom_vol = self%ksAOM * ch4eff * so4                                     ! volumetric, mmol m-3 d-1
   ch4_aerobic_o2 = ch4_aerobic_ox_vol * self%H_REF                      ! areal, mmol O2 m-2 d-1
   ch4_consumption_areal = (ch4_aerobic_ox_vol + ch4_aom_vol) * self%H_REF  ! areal, mmol C m-2 d-1

   total_o2_demand = oxicmin_o2 + nitri_o2 + oduox_o2 + ch4_aerobic_o2

   ! --- diagnostics: NO3 net (nitrification production - denitrification
   ! consumption), reported but not fed back to any pelagic pool ---
   nitri_n_production = ox_frac * n_release
   denitri_n_consumption = 0.8_rk * de_frac * total_C   ! C:NO3 stoichiometry, as in hereon_omexdia_c
   net_no3 = nitri_n_production - denitri_n_consumption

#define _CONV_UNIT_ /secs_pr_day

   ! Own trait state variables: internal source terms only.
   _ADD_BOTTOM_SOURCE_(self%id_oxy_depth,  dLdt   _CONV_UNIT_)
   _ADD_BOTTOM_SOURCE_(self%id_fdet_depth, dLCfdt _CONV_UNIT_)
   _ADD_BOTTOM_SOURCE_(self%id_sdet_depth, dLCsdt _CONV_UNIT_)

   ! Pelagic partners: mass-balance fluxes exposed as diagnostics (mmol
   ! m-2 s-1, positive into the water column) -- NOT applied directly
   ! (see header: this module never touches an interior/"3D" state
   ! variable itself). An external_bottom_flux instance per quantity in
   ! YAML picks these up and applies them to whatever real pelagic state
   ! variable is actually present at the deployment site.
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_oxy_flux,  -total_o2_demand _CONV_UNIT_)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdet_flux, -fdet_demand _CONV_UNIT_)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sdet_flux, -sdet_demand _CONV_UNIT_)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ch4_flux,  (ch4_production - ch4_consumption_areal) _CONV_UNIT_)

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dLdt, dLdt)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ox_frac, ox_frac)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_de_frac, de_frac)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_su_frac, su_frac)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_me_frac, me_frac)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_net_no3, net_no3)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_odu_prod, odu_production)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ch4_prod, ch4_production)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_at_ceiling, at_ceiling_flag)

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

   end module hereon_omexdia_c_bottom