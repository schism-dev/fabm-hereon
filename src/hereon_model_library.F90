!! SPDX-FileCopyrightText: 2022-2025 Helmholtz-Zentrum hereon GmbH
!! SPDX-License-Identifier: CC0-1.0
!! SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
!! SPDX-FileContributor Nina Preußler <nina.preussler@stud.leuphana.de>

module hereon_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: hereon_model_factory

contains

   subroutine create(self, name, model)

      use hereon_light
      use hereon_omexdia_p
      use hereon_omexdia_bottom
      use hereon_omexdia_n2o
      use hereon_omexdia_n2o_nope
      use hereon_nope

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('light');  allocate(type_hereon_light::model)
         case ('omexdia_p'); allocate(type_hereon_omexdia_p::model)
         case ('omexdia_bottom'); allocate(type_hereon_omexdia_bottom::model)
         case ('omexdia_n2o'); allocate(type_hereon_omexdia_n2o::model)
         case ('omexdia_n2o_nope'); allocate(type_hereon_omexdia_n2o_nope::model)
         case ('nope'); allocate(type_hereon_nope::model)
      end select

   end subroutine

end module
