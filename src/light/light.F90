#include "fabm_driver.h"

! This is a FABM implementation of the two-band light model used in
! the General Ocean Turbulence Model (GOTM), extended to spatially variable 
! parameters
!
! @license Gnu General Public License v3.0 or later
!
! @copyright 2015-2021 Bolding and Brueggeman
! @copyright 2021 Helmholtz-Zentrum Hereon

! @author Jorn Brueggeman
! @author Carsten Lemmen
!
! It distinguishes two wavebands:
! * visible (equivalent to photosynthetically activate radation, PAR, 400 - 700 nm)
! * non-visible, which combines ultraviolet (< 400 nm) and infrared (> 700 nm)
! The non-visible fraction is generally absorbed close to the surface.
! The visible fraction typically penetrates deeper into the water.
! Attenuation of the visible band is influenced by FABM variables that contribute to
! standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux
! (as well as by the background attenuation set by parameter g2)
!
! The model is driven by downwelling shortwave radiation just below the water surface,
! which in FABM is denoted with standard_variables%surface_downwelling_shortwave_flux
! This is the radiation left after reflection by the surface [albedo] is accounted for.

module hereon_light

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_hereon_light
      ! Identifiers for dependencies [model inputs]
      type (type_surface_dependency_id) :: id_a    ! non-visible fraction of shortwave radiation
      type (type_surface_dependency_id) :: id_g1   ! e-folding depth of visible fraction
      type (type_surface_dependency_id) :: id_g2   ! e-folding depth of non-visible fraction
      type (type_surface_dependency_id) :: id_swr0 ! Surface shortwave radiation
      type (type_dependency_id)         :: id_dz   ! Cell thickness
      type (type_dependency_id)         :: id_ext  ! Attentuation coefficient for PAR

      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)         :: id_par  ! Photosynthetically active radiation
      type (type_diagnostic_variable_id)         :: id_swr  ! Shortwave radiation
      type (type_surface_diagnostic_variable_id) :: id_par0 ! Surface photosynthetically active radiation

   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_column
   end type type_hereon_light

contains

   subroutine initialize(self, configunit)
      class (type_hereon_light), intent(inout), target :: self
      integer,                 intent(in)            :: configunit

      !call self%get_parameter(elf%a,  'a',  '-','non-visible fraction of shortwave radiation', default=0.58_rk) 
      !call self%get_parameter(self%g1, 'g1', 'm','e-folding depth of non-visible fraction',     default=0.35_rk)
      !call self%get_parameter(self%g2, 'g2', 'm','e-folding depth of visible fraction',         default=23.0_rk) 

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_swr, 'swr', 'W m-2', 'shortwave radiation', &
         standard_variable=standard_variables%downwelling_shortwave_flux, source=source_do_column)
      call self%register_diagnostic_variable(self%id_par, 'par', 'W m-2', 'photosynthetically active radiation', &
         standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux, source=source_do_column)
      call self%register_diagnostic_variable(self%id_par0, 'par0', 'W m-2', 'surface photosynthetically active radiation', &
         standard_variable=standard_variables%surface_downwelling_photosynthetic_radiative_flux, source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_swr0, standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_ext,  standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_dz,   standard_variables%cell_thickness)
      call self%register_dependency(self%id_a,    'a',  '-','novn-isible fraction of shortwave radiation')
      call self%register_dependency(self%id_g1,   'g1', 'm','e-folding depth of non-visible fraction')
      call self%register_dependency(self%id_g2,   'g2', 'm','e-folding depth of visible fraction')
   end subroutine
   
   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_hereon_light), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: swr0, dz, swr, par, z, ext, bioext, a, g1, g2

      _GET_SURFACE_(self%id_swr0,swr0)
      _GET_SURFACE_(self%id_a,a)
      _GET_SURFACE_(self%id_g1,g1)
      _GET_SURFACE_(self%id_g2,g2)
      _SET_SURFACE_DIAGNOSTIC_(self%id_par0,swr0 * (1.0_rk - a))
      z = 0.0_rk
      bioext = 0.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_ext,ext)   ! PAR attenuation (m-1)

         ! Set depth to centre of layer
         z = z + dz * 0.5_rk
         bioext = bioext + ext * dz * 0.5_rk

         ! Calculate photosynthetically active radiation (PAR), shortwave radiation, and PAR attenuation.
         par = swr0 * (1.0_rk - a) * exp(-z / g2 - bioext)
         swr = par + swr0 * a * exp(-z / g1)

         ! Move to bottom of layer
         z = z + dz * 0.5_rk
         bioext = bioext + ext * dz * 0.5_rk

         _SET_DIAGNOSTIC_(self%id_swr,swr) ! Shortwave radiation at layer centre
         _SET_DIAGNOSTIC_(self%id_par,par) ! Photosynthetically active radiation at layer centre
      _DOWNWARD_LOOP_END_
   end subroutine do_column

end module hereon_light
