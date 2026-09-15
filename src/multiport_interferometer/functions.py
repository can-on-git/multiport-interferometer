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
script of utility functions used.
'''

import numpy as np
from pathlib import Path
import sax
import pickle

'''
phi must be in radians
returns power in mW
'''
def phaser_power(phi):
    return phi*23.9/np.pi

def sdict_to_matrix(S):
    in_ports = sorted({p for key in S for p in key if p.startswith("in")}, key=lambda p: int(p[2:]))
    out_ports = sorted({p for key in S for p in key if p.startswith("out")},key=lambda p: int(p[3:]))

    U = np.zeros((len(out_ports), len(in_ports)), dtype=complex)

    for i, out_port in enumerate(out_ports):
        for j, in_port in enumerate(in_ports):
            U[i, j] = np.asarray(S[(in_port, out_port)]).item()

    return U

def bell_circuit_schematic(dim):
    last = dim-2
    end = dim-1

    instances = {}
    connections = {}
    ports = {}

    # creating top boundary connector waveguides
    for j in range(1, dim, 2): 
        instances[f"wvgd_{j}"] = "wvgd"

    # creating bottom boundary phase shifters
    if dim % 2 == 0: 
        for j in range(1, dim, 2):
            instances[f"phaser_{end}_{j}"] = f"phaser_{end}_{j}"
    else:
        for j in range(0, dim, 2):
            instances[f"phaser_{end}_{j}"] = f"phaser_{end}_{j}"

    # creating input and output edge phase shifters
    for i in range(dim):
        instances[f"phaser_i{i}"] = f"phaser_i{i}"
        instances[f"phaser_o{i}"] = f"phaser_o{i}"

    # creating top and bottom mzis
    for i in range(0, dim, 2):
        # top mode mzis
        instances[f"mzm_0_{i}"] = f"mzm_0_{i}"
        if dim % 2 == 0:    # even bottom mode mzis
            instances[f"mzm_{last}_{i}"] = f"mzm_{last}_{i}"
        else:   # odd bottom mode mzis
            j = i + 1
            if j > end:
                continue
            instances[f"mzm_{last}_{j}"] = f"mzm_{last}_{j}"

    # create middle mzis
    for i in range(1, last):
        if i % 2 == 0:
            for j in np.arange(0, dim, 2):
                instances[f"mzm_{i}_{j}"] = f"mzm_{i}_{j}"
        else:
            for j in range(1, dim, 2):
                instances[f"mzm_{i}_{j}"] = f"mzm_{i}_{j}"

    # Connecting Mesh
    for i in range(dim):

        # Top boundary Connections
        if i == 0:
            if dim % 2 == 0:
                for j in range(0, dim - 2, 2):
                    connections[f"mzm_{i}_{j},out0"] = f"wvgd_{j+1},in0"
                    connections[f"wvgd_{j+1},out0"] = f"mzm_{i}_{j+2},in0"
                connections[f"mzm_{i}_{dim-2},out0"] =f"wvgd_{end},in0"
            else:
                for j in range(0, dim - 2, 2):
                    if j < last:
                        connections[f"mzm_{i}_{j},out0"] =f"wvgd_{j+1},in0"
                    if f"mzm_{i}_{j+2}" in instances:
                        connections[f"wvgd_{j+1},out0"] = f"mzm_{i}_{j+2},in0"

        # Connecting diagonal mzis
        for j in range(dim):
            key = f"mzm_{i}_{j}"
            if key not in instances:
                continue
            #up-right connection
            if i > 0:
                upper_right = f"mzm_{i-1}_{j+1}"
                if upper_right in instances:
                    connections[f"mzm_{i}_{j},out0"] = f"mzm_{i-1}_{j+1},in1"
            # downward-right connection
            target = f"mzm_{i+1}_{j+1}"
            if target in instances:
                connections[f"mzm_{i}_{j},out1"] = f"mzm_{i+1}_{j+1},in0"

        # bottom boundary connections
        if i == last:
            if dim % 2 == 0:
                for j in range(1, last, 2):
                    connections[f"mzm_{i}_{j-1},out1"] = f"phaser_{end}_{j},in0"
                    connections[f"phaser_{end}_{j},out0"] = f"mzm_{i}_{j+1},in1"
                # last bottom phaser
                connections[f"mzm_{last}_{last},out1"] = f"phaser_{end}_{end},in0"
            else:
                # first bottom phaser
                connections[f"phaser_{end}_0,out0"] = f"mzm_{last}_1,in1"

                # last bottom phaser
                connections[f"mzm_{last}_{last},out1"] = f"phaser_{end}_{end},in0"

                # intermediate bottom phasers
                for j in range(2, last, 2):
                    connections[f"mzm_{i}_{j-1},out1"] = f"phaser_{end}_{j},in0"
                    connections[f"phaser_{end}_{j},out0"] = f"mzm_{i}_{j+1},in1"

    # connecting input phaseshifters
    for i in range(0, dim - 1, 2):
        connections[f"phaser_i{i},out0"] = f"mzm_{i}_0,in0" #these are sax connections, port naming is swapped with the gds naming because sax component models using standard matrix convention
        connections[f"phaser_i{i+1},out0"] = f"mzm_{i}_0,in1"
        ports[f"in{i}"] = f"phaser_i{i},in0"
        ports[f"in{i+1}"] = f"phaser_i{i+1},in0"

    # odd final input phaseshifter
    if dim % 2 != 0:
        connections[f"phaser_i{end},out0"] = f"phaser_{end}_0,in0"
        ports[f"in{end}"] = f"phaser_i{end},in0"

    # Connecting output phase shifters - Even
    if dim % 2 == 0:
        # top mode
        connections[f"wvgd_{end},out0"] = "phaser_o0,in0"
        ports["out0"] = "phaser_o0,out0"

        # middle modes
        for i in range(1, end, 2):
            connections[f"mzm_{i}_{end},out0"] = f"phaser_o{i},in0"
            connections[f"mzm_{i}_{end},out1"] = f"phaser_o{i+1},in0"
            ports[f"out{i}"] = f"phaser_o{i},out0"
            ports[f"out{i+1}"] = f"phaser_o{i+1},out0"

        # bottom mode
        connections[f"phaser_{end}_{end},out0"] = f"phaser_o{end},in0"
        ports[f"out{end}"] = f"phaser_o{end},out0"

    # Connecting output phase shifters - Even
    else:
        output_col = end
        for i in range(0, dim - 1, 2):
            connections[f"mzm_{i}_{output_col},out0"] = f"phaser_o{i},in0"
            connections[f"mzm_{i}_{output_col},out1"] = f"phaser_o{i+1},in0"
            ports[f"out{i}"] = f"phaser_o{i},out0"
            ports[f"out{i+1}"] = f"phaser_o{i+1},out0"
            
        # bottom mode
        connections[f"phaser_{end}_{end},out0"] = f"phaser_o{end},in0"
        ports[f"out{end}"] = f"phaser_o{end},out0"
    netlist = {}
    netlist['instances'] = instances
    netlist['connections'] = connections
    netlist['ports'] = ports
    return netlist


########################################################################################################################
'''
Utility Functions and models
'''
###############################
"converts decomposition object from phaseshift library into a dict models that contains the sax models with correct naming"
def decom_to_model(component_mods, decom):
    tbu = component_mods["tbu_model"]
    phase_shifter = component_mods["phase_shifter_model"]
    edge_phase_shifter = component_mods["edge_phase_shifter_model"]
    plain_wvgd = component_mods["plain_wvgd_model"]
    dim =decom['m']
    phi_i = decom['phi_ins']
    phi_o = decom['phi_outs']
    phi_ = decom['phi_edges']
    sigma =decom['sigmas']
    delta = decom['deltas']
    models = {}
    models['wvgd'] = plain_wvgd
    for i in range(dim):
        if i in phi_i.keys():
            val = phi_i[i]
            models[f"phaser_i{i}"] = lambda val=val: edge_phase_shifter(val)
        else:
            models[f"phaser_i{i}"] = lambda: edge_phase_shifter(0)
        if i in phi_o.keys():
            val = phi_o[i]
            models[f"phaser_o{i}"] = lambda val=val: edge_phase_shifter(val)
        else:
            models[f"phaser_o{i}"] = lambda: edge_phase_shifter(0)
        for j in range(dim):
            if (i,j) in phi_.keys():
                val=phi_[(i,j)]
                models[f"phaser_{i}_{j}"] = lambda val=val: phase_shifter(val)
            else:
                models[f"phaser_{i}_{j}"] = lambda: phase_shifter(0)
            if (i,j) in sigma.keys():
                s = sigma[(i,j)]
                d = delta[(i,j)]
                models[f"mzm_{i}_{j}"] = lambda s=s, d=d: tbu(sigma=s, delta=d )
            else:
                models[f"mzm_{i}_{j}"] = lambda:tbu(0,0)
    return models
