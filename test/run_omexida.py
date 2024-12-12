# SPDX-FileCopyrightText: 2024 Helmholtz-Zentrum hereon GmbH (Hereon)
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import scipy.integrate
from matplotlib import pyplot as plt
import sys, os
from importlib import reload

import pyfabm
reload(pyfabm) # to update an already loaded one.

def check_dependencies(model: pyfabm.Model, dependencies: dict):

    for variable in model.dependencies:
        if variable.name in dependencies:
            model.findDependency(variable.name).value = dependencies[variable.name]
            print ('Set dependency {} to {} {}'.format(variable.name,
                                                       dependencies[variable.name],variable.units))
        else:
            print ('Dependency {} with units {} not met'.format(variable.name,
                                                                variable.units))
            sys.exit(1)

def run_model(model: pyfabm.Model):

    # Time derivative
    def dy(y, t0):
        model.state[:] = y
        return model.getRates()

    # Time-integrate  (note: FABM's internal time unit is seconds!)
    t = np.linspace(0, 500, num=100)
    y = scipy.integrate.odeint(dy, model.state, t)
    return t,y

def plot_model(model: pyfabm.Model, t, y):

    names=[v.path + ' ' + v.units for v in model.state_variables]
    legend=[]
    for i in range(0,y.shape[1]):
        #print ('{} {}'.format(names[i],y[:,i]))
        if 'depth' in names[i]:
            names[i] = names[i].replace(' m',' mm')
            plt.plot(t, 1000*y[:,i])#/np.max(y[:,i]))
        else:
            plt.plot(t, y[:,i])#/np.max(y[:,i]))
        legend.append(names[i])

    plt.legend(legend)

    name = "--".join(list(model._get_parameter_tree().keys()))

    plt.savefig(f'{name}_timeseries.png')
    plt.close()

def main():

    model = pyfabm.Model('fabm-omexdia.yaml')
    model.cell_thickness=0.5

    dependencies = {'temperature': 10}

    check_dependencies(model, dependencies)
    assert model.checkReady()

    t,y = run_model(model)
    plot_model(model, t, y)
    return model


if __name__ == '__main__':

    print('Running omexdia ...')
    omexdia = main()
