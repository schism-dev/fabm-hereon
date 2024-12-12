# FABM-hereon

This is a collection of [FABM](https://fabm.net) models developed at Helmholtz-Zentrum hereon GmbH by the ecosystem modeling group.  It currently includes

1. a port of the GOTM light model with horizontal parameter variation
2. the OmexDia model with added phosphorous
3. the OmexDia model as a bottom representation (in development)

## How to build

This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=</path/to/fabm-hereon>`

Here, `</path/to/fabm-hereon>` is the directory with the FABM-hereon code (the same directory that contains this ReadMe file). Note that `-DFABM_INSTITUTES=hereon` will make FABM compile OmexDia and light as the *only* available biogeochemical models. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="hereon;iow"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

## Building with GOTM
To use hereon models with the General Ocean Turbulence Model (GOTM), do the following, after defining GOTM_BASE and HEREON_BASE

```
git clone --recursive https://github.com/gotm-model/code.git $GOTM_BASE
git clone https://github.com/schism-dev/fabm-hereon.git $HEREON_BASE

cd $GOTM_BASE && git submodule update --init --recursive
mkdir -p $GOTM_BASE/build
cd $GOTM_BASE/build 
cmake -B $GOTM_BASE/build -S $GOTM_BASE  -DGOTM_USE_FABM=ON -DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=$HEREON_BASE
make
```

This will create the `gotm` executable with support for FABM light and omexdia.

## Building with SCHISM
To use hereon models with the latest stable release of the Semi-implicit Cross-scale Hydroscience Integrated System Model (SCHISM), do the following:

```
git clone --recurse-submodules -b v5.11.1 https://github.com/schism-dev/schism.git schism
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/schism-dev/fabm-hereon.git
mkdir build
cd build
cmake ../schism/src -DBLD_STANDALONE=ON -DUSE_FABM=ON -DFABM_BASE=../fabm -DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=../fabm-hereon
make pschism
```

This will create the `pschism` executable with support for FABM light and omexdia.

## Building the 0D driver

To use hereon models with pyfabm, do the following, after defining FABM_BASE, GOTM_BASE and HEREON_BASE. 

```
mkdir $FABM_BASE/build-0d
cmake -B $FABM_BASE/build-0d -S $FABM_BASE/src/drivers/0d -DGOTM_BASE=$GOTM_BASE -DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=$HEREON_BASE
make
```


## Building pyfabm

To use hereon models with pyfabm, do the following, after defining FABM_BASE and HEREON_BASE. See also the detailed instructions at https://github.com/fabm-model/fabm/wiki/python.

```
cat <<EOT > $FABM_BASE/setup.cfg
[build_ext]
cmake_opts=-DFABM_EXTRA_INSTITUTES=hereon -DFABM_HEREON_BASE="$HEREON_BASE"
force=1
debug=1
EOT

python -m pip install $FABM_BASE/setup.cfg
````
