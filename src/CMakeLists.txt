add_library(fabm_models_hereon OBJECT
            hereon_model_library.F90
            light/light.F90
	        omexdia/omexdia_p.F90
	        omexdia/omexdia_bottom.F90
       omexdia/omexdia_n2o.F90
	    gotm/npzd.F90
	    gotm/ergom.F90
        oxypom/oxypom.F90
           )

add_dependencies(fabm_models_hereon fabm_base)
