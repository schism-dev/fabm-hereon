# SPDX-FileCopyrightText: 2022-2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
# SPDX-FileContributor Nina Preußer <nina.preusser@stud.leuphana.de>

add_library(fabm_models_hereon OBJECT
   hereon_model_library.F90
   light/light.F90
	omexdia/omexdia_p.F90
	omexdia/omexdia_bottom.F90
   omexdia/omexdia_n2o.F90
)

add_dependencies(fabm_models_hereon fabm_base)
