'''
Copyright 2026 Carl Abi Nakad

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

'''
Library of ideal models used in this project
'''

import numpy as np
import sax

def ideal_bell_models(decom):
    dim =decom['m']
    phi_i = decom['phi_ins']
    phi_o = decom['phi_outs']
    phi_ = decom['phi_edges']
    sigma =decom['sigmas']
    delta = decom['deltas']
    models = {}
    models['wvgd'] = ideal_waveguide
    for i in range(dim):
        if i in phi_i.keys():
            val = phi_i[i]
            models[f"phaser_i{i}"] = lambda val=val: ideal_phaseshifter(val)
        else:
            models[f"phaser_i{i}"] = ideal_phaseshifter
        if i in phi_o.keys():
            val = phi_o[i]
            models[f"phaser_o{i}"] = lambda val=val: ideal_phaseshifter(val)
        else:
            models[f"phaser_o{i}"] = ideal_phaseshifter
        for j in range(dim):
            if (i,j) in phi_.keys():
                val=phi_[(i,j)]
                models[f"phaser_{i}_{j}"] = lambda val=val: ideal_phaseshifter(val)
            else:
                models[f"phaser_{i}_{j}"] = ideal_phaseshifter
            if (i,j) in sigma.keys():
                s = sigma[(i,j)]
                d = delta[(i,j)]
                models[f"mzm_{i}_{j}"] = lambda s=s, d=d: ideal_mzi(sigma=s, delta=d )
            else:
                models[f"mzm_{i}_{j}"] = ideal_mzi
    return models

def ideal_waveguide():
    return sax.reciprocal({("in0", "out0"): 1})

def ideal_phaseshifter(phi=0):
    return sax.reciprocal({('in0','out0'):np.exp(1j*phi)})

def ideal_mzi(sigma=0, delta=0):
    phase = np.exp(1j*sigma)
    a = phase*np.sin(delta)
    b = phase*np.cos(delta)
    return sax.reciprocal({('in0','out0'): a, ('in0','out1'):b, ('in1','out0'):b, ('in1','out1'):-a})

def ideal_clements_tbu_model(theta, phi):
    phase = np.exp(1j*phi)
    a=phase*np.cos(theta)
    b = -np.sin(theta)
    c = phase*np.sin(theta)
    d = np.cos(theta)
    return sax.reciprocal({('in0','out0'): a, ('in0','out1'):c, ('in1','out0'):b, ('in1','out1'):d})

