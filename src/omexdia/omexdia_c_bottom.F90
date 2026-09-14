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
!   id_so4  -- interfacial sulfate, OPTIONAL (required=.false.). Prefer
!              coupling this to real observations/forcing wherever
!              possible -- real station data (Elbe estuary: Geesthacht,
!              salinity=0, observed so4=800 mmol/m3) showed the
!              salinity-derived proxy below fails badly in freshwater
!              reaches with non-marine (geological/mining-influenced)
!              sulfate sources. so4 was made a REQUIRED, INDEPENDENT
!              dependency for exactly this reason in an earlier version
!              of this module.
!   id_salinity -- interfacial salinity, OPTIONAL (required=.false.),
!              used ONLY as a fallback when so4_c0 is not coupled (see
!              do_bottom: so4_c0 if available, else derived from
!              salinity_c0 via the same sulfate_from_salinity() proxy
!              hereon_omexdia_c uses internally -- duplicated locally,
!              not `use`d from that module, to keep this one
!              self-contained; carries the SAME freshwater-failure
!              caveat as above). If NEITHER is coupled, so4 falls back
!              to the fixed parameter so4_default (full-marine by
!              default) -- see initialize. This is deliberately a
!              last-resort fallback chain, not an endorsement of the
!              proxy: supply real so4_c0 whenever it is available.
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
! EXPONENTIAL-SATURATION KINETICS (NOT Monod/Michaelis-Menten):
! ---------------------------------------------------------------
! Every limitation/inhibition factor below (Oxicminlim, the O2/NO3
! factors inside Denitrilim and anoxic_space, f_so4/f_meth, and the CH4
! aerobic-oxidation factor) uses
!    limitation(S) = 1 - exp(-S/Ks)      [0 at S=0, -> 1 as S->infinity]
!    inhibition(S) =     exp(-S/Ks)      [exact complement, = 1-limitation(S)]
! instead of the classic Monod/Michaelis-Menten S/(S+Ks) / Ks/(S+Ks) pair.
! Same qualitative 0->1 saturating shape and the same exact-complement
! structure, but the derivative of exp(-S/Ks) stays IN THE SAME
! exponential family (d/dS[exp(-S/Ks)] = -(1/Ks)*exp(-S/Ks)) rather than
! escalating in rational-function degree at every differentiation --
! this matters here specifically because dCostdL/dLdt repeatedly
! differentiate these terms through the fitness-gradient machinery.
!
! CALIBRATION NOTE: Monod hits half-max exactly at S=Ks; the exponential
! form hits half-max at S=Ks*ln(2)~=0.693*Ks. To keep the SAME half-max
! location (so this swap changes the functional shape but not the
! calibrated operating point), every affected Ks parameter's default
! below (ksO2oxic, kinO2denit, ksNO3denit, kinNO3anox, kinO2anox,
! kinSO4, ksCH4) is the ORIGINAL Monod-calibrated default divided by
! ln(2) -- e.g. ksO2oxic default 3.0 (Monod half-max at O2=3) becomes
! 3.0/ln(2)~=4.328 (exponential half-max still at O2=3). If these
! parameters are overridden in a deployment's YAML, remember they now
! calibrate an exponential-saturation curve, not a Monod one.
!
! NOTE this is a DIFFERENT reformulation from the depth-INTEGRATION
! question (see the porosity-weighting and inventory sections below):
! composing exp(-S/Ks) through a depth profile S(z)=S0*exp(-z/L) gives a
! Gompertz-type double-exponential in z with NO elementary closed-form
! integral (worse than Monod-of-exponential-profile, which does
! integrate via logs/partial fractions) -- so this swap is deliberately
! scoped to the INTERFACE-evaluated (z=0) limitation fractions only. Any
! future depth-resolved energy/fraction profile should prescribe a
! direct exponential envelope in z rather than compose this substitute
! through S(z).
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
! SEDIMENT INVENTORY DIAGNOSTICS (standing stock, mmol m-2, NOT a flux):
! total_oxygen_in_soil/fdet_inventory/sdet_inventory integrate the
! model's own assumed exponential depth profiles
! (X(z) = X_surface*exp(-z/depth_trait), the same shape marg_fdet/
! marg_sdet/benefit are already derived from -- no new assumption
! introduced) over a FIXED 0-to-Lmax box:
!   inventory = X_surface * depth_trait * (1 - exp(-Lmax/depth_trait))
! Lmax (not e.g. a multiple of the trait itself) is used deliberately as
! the upper bound so the accounting box is a fixed reference volume, not
! one that silently grows/shrinks with the trait it is trying to budget.
! total_carbon_in_soil/total_nitrogen_in_soil sum fdet+sdet (N via the
! existing NCrFdet/NCrSdet ratios) for a combined particulate organic
! C/N budget.
!
! total_phosphorus_in_soil: unlike N, hereon_omexdia_c does NOT give P a
! fixed ratio to C -- pdet and po4 are their own pair of prognostic
! pools there, exchanging via sorption/desorption:
!   radsP = PAds*rSlow*(po4*max(odu,PAdsODU))    [adsorption, po4->pdet]
!   Pprod = rFast*(1-Oxicminlim)*pdet            [release,    pdet->po4]
! (P sorbs onto particulates under oxic conditions and releases back to
! porewater as conditions turn anoxic -- the Fe-oxide "P shuttle"). This
! module does NOT reproduce that sorption/desorption KINETICS -- doing
! so would need a depth profile for odu (the redox proxy controlling
! the balance), which nothing in this module posits (odu is only ever a
! quasi-steady RATE here, odu_production, never a standing
! concentration). Instead, total_phosphorus_in_soil sums two
! independently-read, non-reacting forcing pools:
!   pdet_c0 (mmol P m-3) -- particulate, assumed to share fdet's OWN
!     spatial scale (fdet_depth), justified because pdet's release rate
!     in omexdia_c uses the same rFast family as fdet's, not sdet's.
!   po4_c0 (mmol P m-3) -- dissolved, assumed to share oxy's spatial
!     scale (oxy_depth) as a dissolved porewater species, the same way
!     O2 does -- an UNVALIDATED assumption (po4 release actually peaks
!     BELOW the oxic zone, not decaying alongside it like O2 does; no
!     better alternative is available without inventing a real po4
!     depth-profile submodel, which is out of scope here).
! Both pools are read-only forcing, exactly like oxy_c0/fdet_c0/sdet_c0
! -- neither is derived from the other via the real sorption equations
! above, unlike in omexdia_c itself.
!
! NOT provided for CH4/ODU/NO3/SO4/sulfur: this model never posits a
! depth profile for any of them (CH4 is a real pelagic concentration
! with no sediment shape of its own here; ODU/NO3/SO4 are diagnostic-
! only/forcing; organic sulfur is not tracked at all) -- an inventory
! for those would require inventing a profile shape (or, for sulfur, an
! S:C ratio) not otherwise used anywhere else in this module.
!
! POROSITY WEIGHTING OF INVENTORIES (id_porosity, id_z_poros, both
! OPTIONAL HORIZONTAL dependencies): all five inventories above (O2,
! fdet, sdet, pdet, po4) previously integrated as if concentrations
! filled 100% of the bulk sediment volume at every depth. Real
! sediments are only partly pore water -- porosity phi(z) is the pore-
! water volume fraction, and it declines with depth via compaction
! (Berner's classic law):
!   phi(z) = porosity_inf + (porosity - porosity_inf)*exp(-z/z_poros)
! DISSOLVED species (oxy, po4) are conventionally expressed per unit
! PORE-WATER volume, so their true inventory needs a phi(z) weight.
! PARTICULATE species (fdet, sdet, pdet) are conventionally expressed
! per unit SOLID volume, so their true inventory needs a (1-phi(z))
! weight instead -- NOT the same weight as the dissolved species (that
! would get the depth-dependence backwards: compaction makes deeper
! sediment MORE solid, not less, the opposite direction phi(z) itself
! moves). Both weights are handled by ONE helper function,
! porosity_weighted_inventory(), since (1-phi(z)) is exactly the same
! two-exponential shape as phi(z) with different constants -- see that
! function for the closed-form integral (composes with the existing
! X_surface*exp(-z/scale) profiles into a second exponential term, no
! numerical integration needed).
!
! porosity_inf (asymptotic deep porosity) is a plain fixed PARAMETER,
! not a dependency -- only the SURFACE porosity and the compaction
! depth scale are exposed as (optional) horizontal dependencies, with
! fixed-parameter fallbacks (porosity_default=0.7, z_poros_default=
! 0.15 m, porosity_inf=0.5) if not coupled. All three defaults are
! rough literature-typical values for a muddy continental shelf,
! UNVALIDATED against any real site -- couple id_porosity/id_z_poros to
! a real sediment-property source (a grain-size/sediment-type map, a
! separate compaction submodel, etc.) wherever possible.
!
! DELIBERATELY SCOPED to the standing-stock inventories only -- this
! does NOT touch the reaction-rate/flux formulas (fdet_demand,
! oxicmin_o2, marg_fdet, every _ADD_BOTTOM_SOURCE_/_ADD_BOTTOM_FLUX_-
! adjacent term, etc., all still implicitly bulk-volume). Those rate
! constants (rFast, rnit, ...) were never calibrated with porosity in
! mind; retrofitting it into the actual dynamics (not just these
! diagnostics) would be a real behavior change requiring its own
! validation, and is out of scope here.
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
!     Interfacial P concentrations, read-only, particulate + dissolved
!     (see header: used only for total_phosphorus_in_soil, summed as
!     two independently-read forcing pools -- no sorption/desorption
!     coupling between them is reproduced).
      type (type_dependency_id)                      :: id_pdet
      type (type_dependency_id)                      :: id_po4

!     Read-only forcing dependencies (no flux returned).
      type (type_dependency_id)                      :: id_no3
!     so4 OPTIONAL (see header); salinity OPTIONAL, used only as a
!     fallback to derive so4 when so4_c0 itself is not coupled.
      type (type_dependency_id)                      :: id_so4
      type (type_dependency_id)                      :: id_salinity
      type (type_dependency_id)                      :: id_temp
!     Sediment porosity properties, OPTIONAL HORIZONTAL dependencies
!     (genuinely horizontal -- a sediment-column property, not a
!     pelagic one sampled at the interface, unlike oxy_c0 etc. above).
!     Used only for porosity-weighting the inventory diagnostics; see
!     header. Fixed-parameter fallbacks if not coupled.
      type (type_horizontal_dependency_id)          :: id_porosity
      type (type_horizontal_dependency_id)          :: id_z_poros
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
!     Sediment inventory diagnostics (standing stock, mmol m-2, 0-Lmax
!     box; NOT a flux -- see header).
      type (type_horizontal_diagnostic_variable_id) :: id_totO2_soil
      type (type_horizontal_diagnostic_variable_id) :: id_fdet_inventory, id_sdet_inventory
      type (type_horizontal_diagnostic_variable_id) :: id_totC_soil, id_totN_soil, id_totP_soil

!     Model parameters
      real(rk) :: rate_L, alpha, rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: LCfast_max, Lhalf_fast, tau_fast
      real(rk) :: LCslow_max, Lhalf_slow, tau_slow
      real(rk) :: Lmin, LCmin, Lmax
      real(rk) :: ksO2oxic, kinO2denit, ksNO3denit
      real(rk) :: kinNO3anox, kinO2anox, kinSO4, relaxO2
      real(rk) :: q10_meth, Tref
      real(rk) :: nh3_amb, odu_amb
      real(rk) :: so4_default
      real(rk) :: porosity_default, z_poros_default, porosity_inf
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

   ! --- electron-acceptor partition kinetics (match hereon_omexdia_c's
   ! ORIGINAL Monod half-max calibration -- see header: these are now
   ! exponential-saturation curves, and each default below is the
   ! original Monod Ks divided by ln(2) so the half-max location is
   ! unchanged) ---
   call self%get_parameter(self%ksO2oxic,'ksO2oxic','mmolO2 m-3', &
        'O2 scale in oxic mineralization (exponential-saturation, half-max at '// &
        'O2=ksO2oxic*ln(2); original Monod default 3.0 rescaled, see header)', default=4.328_rk)
   call self%get_parameter(self%kinO2denit,'kinO2denit','mmolO2 m-3', &
        'O2 scale, inhibition of denitrification (exponential-saturation; original '// &
        'Monod default 70.0 rescaled, see header)', default=100.987_rk)
   call self%get_parameter(self%ksNO3denit,'ksNO3denit','mmolNO3 m-3', &
        'NO3 scale in denitrification (exponential-saturation; original Monod default '// &
        '1.0 rescaled, see header)', default=1.443_rk)
   call self%get_parameter(self%kinNO3anox,'kinNO3anox','mmolNO3 m-3', &
        'NO3 scale, inhibition of anoxic mineralization (exponential-saturation; '// &
        'original Monod default 1.0 rescaled, see header)', default=1.443_rk)
   call self%get_parameter(self%kinO2anox,'kinO2anox','mmolO2 m-3', &
        'O2 scale, inhibition of anoxic mineralization (exponential-saturation; '// &
        'original Monod default 1.0 rescaled, see header)', default=1.443_rk)
   call self%get_parameter(self%kinSO4,'kinSO4','mmol m-3', &
        'sulfate scale, inhibition threshold for methanogenesis (exponential-saturation; '// &
        'original Monod default 1000.0 rescaled, see header)', default=1443.1_rk)
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
   call self%get_parameter(self%so4_default,'so4_default','mmol m-3', &
        'last-resort so4 fallback used ONLY if neither so4_c0 nor salinity_c0 is coupled '// &
        '(see header) -- full-marine (S=35) by default; NOT used at all if either dependency '// &
        'is actually supplied', default=28000.0_rk)

   ! --- porosity (inventory-weighting only, see header; NOT used in any
   ! reaction/flux formula) ---
   call self%get_parameter(self%porosity_default,'porosity_default','-', &
        'surface porosity fallback used ONLY if id_porosity is not coupled -- rough '// &
        'muddy-shelf default, UNVALIDATED (see header)', default=0.7_rk)
   call self%get_parameter(self%z_poros_default,'z_poros_default','m', &
        'compaction attenuation-depth fallback used ONLY if id_z_poros is not coupled -- '// &
        'rough muddy-shelf default, UNVALIDATED (see header)', default=0.15_rk)
   call self%get_parameter(self%porosity_inf,'porosity_inf','-', &
        'asymptotic deep (fully compacted) porosity -- always a fixed parameter, not a '// &
        'dependency, unlike surface porosity/compaction depth above', default=0.5_rk)

   ! --- CH4 parameters ---
   call self%get_parameter(self%rmaxO2,'rmaxO2','d-1', &
        'maximum aerobic CH4 oxidation rate', default=10.0_rk)
   call self%get_parameter(self%ksCH4,'ksCH4','mmol m-3', &
        'CH4 scale for aerobic oxidation (exponential-saturation; original Monod '// &
        'default 5.0 rescaled, see header)', default=7.213_rk)
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

   ! --- sediment inventory diagnostics (standing stock, 0-Lmax box; see
   ! header -- NOT provided for CH4/ODU/NO3/SO4/sulfur, no depth profile
   ! or ratio posited for any of those) ---
   call self%register_diagnostic_variable(self%id_totO2_soil,'total_oxygen_in_soil','mmol O2 m-2', &
        'O2 standing stock in the 0-Lmax sediment box, assuming O2(z)=oxy_c0*exp(-z/oxy_depth)', &
        output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_fdet_inventory,'fdet_inventory','mmol C m-2', &
        'fdet standing stock in the 0-Lmax sediment box, assuming fdet(z)=fdet_c0*exp(-z/fdet_depth)', &
        output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_sdet_inventory,'sdet_inventory','mmol C m-2', &
        'sdet standing stock in the 0-Lmax sediment box, assuming sdet(z)=sdet_c0*exp(-z/sdet_depth)', &
        output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_totC_soil,'total_carbon_in_soil','mmol C m-2', &
        'total particulate organic carbon standing stock (fdet+sdet) in the 0-Lmax box', &
        output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_totN_soil,'total_nitrogen_in_soil','mmol N m-2', &
        'total particulate organic nitrogen standing stock (fdet*NCrFdet+sdet*NCrSdet) '// &
        'in the 0-Lmax box', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_totP_soil,'total_phosphorus_in_soil','mmol P m-2', &
        'particulate (pdet_c0, fdet_depth scale) + dissolved (po4_c0, oxy_depth scale) P '// &
        'standing stock in the 0-Lmax box -- two independently-read forcing pools, NOT '// &
        'reproducing omexdia_c own sorption/desorption coupling between them (see header note)', &
        output=output_instantaneous, source=source_do_bottom)

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
   ! Used only for total_phosphorus_in_soil (see header) -- no flux
   ! returned, unlike oxy/fdet/sdet/ch4 above. Two independent pools,
   ! NOT coupled via sorption/desorption the way omexdia_c couples them.
   call self%register_dependency(self%id_pdet,'pdet_c0','mmol P m-3', &
        'particulate detritus-P concentration at the sediment-water interface')
   call self%register_dependency(self%id_po4,'po4_c0','mmol P m-3', &
        'dissolved phosphate concentration at the sediment-water interface')

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
   call self%register_dependency(self%id_no3,'no3_c0','mmol N m-3', &
        'nitrate concentration at the sediment-water interface (diagnostic use only)')
   ! OPTIONAL -- prefer coupling this explicitly to real observations/
   ! forcing wherever possible (see header caveat on the salinity-proxy
   ! fallback below).
   call self%register_dependency(self%id_so4,'so4_c0','mmol m-3', &
        'sulfate concentration at the sediment-water interface (prefer coupling explicitly -- '// &
        'see header: falls back to a salinity proxy, then a fixed default, if not coupled)', &
        required=.false.)
   ! OPTIONAL -- fallback source for so4 only, see header and do_bottom.
   call self%register_dependency(self%id_salinity,'salinity_c0','PSU', &
        'salinity at the sediment-water interface (fallback so4 proxy only, see header caveat '// &
        '-- unused if so4_c0 is coupled)', required=.false.)
   ! NOT standard_variables%temperature -- no bottom-domain temperature
   ! standard variable exists in this FABM tree, and the interior one
   ! would reintroduce exactly the "3D" coupling this module avoids
   ! elsewhere (see header). Plain named dependency instead, same
   ! convention as no3_c0/so4_c0.
   call self%register_dependency(self%id_temp,'temp_c0','degree_C', &
        'bottom water temperature (methanogenesis Q10 factor only)')

   ! --- porosity, OPTIONAL HORIZONTAL dependencies (see header) --
   ! inventory-weighting only, falls back to porosity_default/
   ! z_poros_default if not coupled.
   call self%register_dependency(self%id_porosity,'porosity','-', &
        'surface (interfacial) porosity, used only to weight the sediment inventory '// &
        'diagnostics -- couple to a real sediment-property source where possible, see header', &
        required=.false.)
   call self%register_dependency(self%id_z_poros,'z_poros','m', &
        'compaction attenuation depth (porosity decline with depth), used only to weight '// &
        'the sediment inventory diagnostics -- see header', required=.false.)

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
   real(rk) :: salinity
   real(rk) :: porosity_0, z_poros
   real(rk) :: pdet_surface, po4_surface
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
   real(rk) :: oxy_inventory, fdet_inventory, sdet_inventory, totC_inventory, totN_inventory, totP_inventory
   real(rk) :: pdet_inventory, po4_inventory

!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   ! Interfacial (bottom-most pelagic cell) concentrations and forcing
   _GET_(self%id_oxy,oxy_surface)
   _GET_(self%id_fdet,fdet_surface)
   _GET_(self%id_sdet,sdet_surface)
   _GET_(self%id_ch4,ch4_surface)
   _GET_(self%id_pdet,pdet_surface)
   _GET_(self%id_po4,po4_surface)
   _GET_(self%id_no3,no3_surface)
   _GET_(self%id_temp,temp_celsius)

   ! so4 fallback chain (see header): so4_c0 if coupled, else derived
   ! from salinity_c0 if THAT is coupled, else a fixed default. Checked
   ! here (not at initialize) because optional-dependency availability
   ! is only resolved after all models' coupling requests are in.
   if (_AVAILABLE_(self%id_so4)) then
      _GET_(self%id_so4,so4)
   else if (_AVAILABLE_(self%id_salinity)) then
      _GET_(self%id_salinity,salinity)
      so4 = sulfate_from_salinity(salinity)
   else
      so4 = self%so4_default
   end if

   ! porosity/compaction-depth: same optional-with-fallback pattern
   ! (see header). Inventory-weighting only.
   if (_AVAILABLE_HORIZONTAL_(self%id_porosity)) then
      _GET_HORIZONTAL_(self%id_porosity,porosity_0)
   else
      porosity_0 = self%porosity_default
   end if
   if (_AVAILABLE_HORIZONTAL_(self%id_z_poros)) then
      _GET_HORIZONTAL_(self%id_z_poros,z_poros)
   else
      z_poros = self%z_poros_default
   end if

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

   ! Exponential-saturation kinetics, NOT Monod -- see header. limitation
   ! (1-exp(-S/Ks)) and inhibition (exp(-S/Ks)) are exact complements,
   ! same as the Monod pair they replace; Ks parameters are pre-rescaled
   ! (see initialize) so the half-max location matches the original
   ! Monod calibration.
   Oxicminlim = 1.0_rk - exp(-oxy_surface / (self%ksO2oxic + self%relaxO2*(self%nh3_amb + self%odu_amb)))
   Denitrilim = exp(-oxy_surface/self%kinO2denit) * (1.0_rk - exp(-no3_surface/self%ksNO3denit))
   anoxic_space = exp(-oxy_surface/self%kinO2anox) * exp(-no3_surface/self%kinNO3anox)
   f_so4 = 1.0_rk - exp(-so4/self%kinSO4)
   f_meth = exp(-so4/self%kinSO4)
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
   ch4_aerobic_ox_vol = self%rmaxO2 * (1.0_rk - exp(-ch4eff/self%ksCH4)) * ox_frac  ! volumetric, mmol m-3 d-1 (exponential-saturation, see header)
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

   ! --- sediment inventory diagnostics (standing stock, 0-Lmax box,
   ! porosity-weighted; see header) -- same exponential profile shape
   ! marg_fdet/marg_sdet/benefit above are already derived from, just
   ! integrated to a fixed bound (Lmax) instead of the trait's own
   ! (moving) depth, and now weighted by the pore-water (dissolved) or
   ! solid (particulate) volume fraction instead of assuming 100% bulk
   ! volume. porosity_weighted_inventory()'s w0/winf arguments select
   ! which: (porosity_0, porosity_inf) for dissolved, (1-porosity_0,
   ! 1-porosity_inf) for particulate -- see that function. ---
   oxy_inventory  = porosity_weighted_inventory(oxy_surface,  Leff, &
        porosity_0, self%porosity_inf, z_poros, self%Lmax)
   fdet_inventory = porosity_weighted_inventory(fdet_surface, LCfeff, &
        1.0_rk-porosity_0, 1.0_rk-self%porosity_inf, z_poros, self%Lmax)
   sdet_inventory = porosity_weighted_inventory(sdet_surface, LCseff, &
        1.0_rk-porosity_0, 1.0_rk-self%porosity_inf, z_poros, self%Lmax)
   totC_inventory = fdet_inventory + sdet_inventory
   totN_inventory = self%NCrFdet*fdet_inventory + self%NCrSdet*sdet_inventory
   ! pdet assumed to share fdet's own spatial scale (fdet_depth); po4
   ! assumed to share oxy's own spatial scale (oxy_depth), as a
   ! dissolved porewater species -- see header note on both assumptions,
   ! and on the sorption/desorption coupling deliberately NOT reproduced
   ! between them.
   pdet_inventory = porosity_weighted_inventory(pdet_surface, LCfeff, &
        1.0_rk-porosity_0, 1.0_rk-self%porosity_inf, z_poros, self%Lmax)
   po4_inventory  = porosity_weighted_inventory(po4_surface,  Leff, &
        porosity_0, self%porosity_inf, z_poros, self%Lmax)
   totP_inventory = pdet_inventory + po4_inventory

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totO2_soil,     oxy_inventory)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdet_inventory, fdet_inventory)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sdet_inventory, sdet_inventory)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totC_soil,      totC_inventory)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totN_soil,      totN_inventory)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totP_soil,      totP_inventory)

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

   ! Fallback so4 proxy, used ONLY when so4_c0 is not coupled (see
   ! header caveat -- duplicated from hereon_omexdia_c's own private
   ! function of the same name rather than `use`d, to keep this module
   ! self-contained; carries the SAME freshwater-failure caveat: real
   ! Elbe estuary data (Geesthacht, salinity=0) showed observed so4=800
   ! mmol/m3, not the 0 this linear proxy predicts).
   !
   ! Morris, A. W., & Riley, J. P. (1966). The bromide/chlorinity and
   ! sulphate/chlorinity ratio in sea water. Deep Sea Research and
   ! Oceanographic Abstracts, 13(4), 699-705.
   ! https://doi.org/10.1016/0011-7471(66)90601-2
   elemental function sulfate_from_salinity(salinity) result(so4)
      real(rk), intent(in) :: salinity ! [PSU] local salinity, range 0..35
      real(rk)             :: so4      ! [mmol m-3] range 0..28000

      so4 = (28000.0_rk / 35.0_rk) * max(0.0_rk, salinity)
   end function sulfate_from_salinity

   ! Porosity-weighted depth integral, 0 to Zmax, of a quantity whose
   ! own profile is X_surface*exp(-z/scale) (the same exponential shape
   ! every trait-driven profile in this module already has), weighted
   ! by a Berner-type compaction profile w(z) = winf + (w0-winf)*
   ! exp(-z/z_poros). Call with w0=porosity/winf=porosity_inf for a
   ! DISSOLVED quantity, or w0=(1-porosity)/winf=(1-porosity_inf) for a
   ! PARTICULATE one -- see header for why these must NOT be swapped.
   !
   ! Closed form: w(z)*exp(-z/scale) splits into two pure exponentials
   ! (rate 1/scale, and rate 1/scale+1/z_poros), each integrating the
   ! same way the model's un-weighted inventories already did:
   !   integral = X_surface * [ winf*scale*(1-exp(-Zmax/scale))
   !            + (w0-winf)*scale2*(1-exp(-Zmax/scale2)) ]
   !   scale2 = scale*z_poros/(scale+z_poros)
   ! No numerical integration -- this is exact given the assumed
   ! exponential profile and exponential compaction law.
   elemental function porosity_weighted_inventory(X_surface, scale, w0, winf, z_poros, Zmax) result(inv)
      real(rk), intent(in) :: X_surface, scale, w0, winf, z_poros, Zmax
      real(rk)             :: inv, scale2

      scale2 = scale*z_poros / (scale+z_poros)
      inv = X_surface * ( winf*scale*(1.0_rk - exp(-Zmax/scale)) &
                         + (w0-winf)*scale2*(1.0_rk - exp(-Zmax/scale2)) )
   end function porosity_weighted_inventory

   end module hereon_omexdia_c_bottom