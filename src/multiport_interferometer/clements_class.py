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

import interferometer as itf
from . import ideal_models as ideal
import numpy as np
import sax

class Clements:
    def __init__(self, dimension):
        if not isinstance(dimension, int) or dimension < 2:
            raise ValueError("Dimension must be an integer >= 2")
        self.dimension = dimension
        self.schematic_netlist = self.clements_circuit_schematic()
    
    def ideal_emulate(self,U):
        if U.shape[0]!=self.dimension:
            raise ValueError(f'Matrix Dimension {U.shape[0]} not compatible with Mesh Dimension {self.dimension}')
        mods = self.clements_models_ideal(U)
        netlist = self.schematic_netlist
        circuit, info = sax.circuit(netlist=netlist, models = mods)
        S = circuit()
        return S
    
    def clements_circuit_schematic(self):
        dim = self.dimension
        last = dim-2
        end = dim-1

        instances = {}
        connections = {}
        ports = {}

        # creating top boundary connector waveguides
        for j in range(1, dim, 2): 
            instances[f"wvgd_top_{j}"] = "wvgd"

        #create bottom boundary connector waveguides
        if dim % 2 == 0:
            for j in range(1, dim, 2):
                instances[f"wvgd_bottom_{j}"] = "wvgd"
        else:
            for j in range(0, dim, 2):
                instances[f"wvgd_bottom_{j}"] = "wvgd"

        # creating output edge phase shifters
        for i in range(dim):
            instances[f"phaser_o{i}"] = f"phaser_o{i}"

        # creating top and bottom mzis
        for i in range(0, dim, 2):
            # top mode mzis
            instances[f"bs_0_{i}"] = f"bs_0_{i}"
            if dim % 2 == 0:    # even bottom mode mzis
                instances[f"bs_{last}_{i}"] = f"bs_{last}_{i}"
            else:   # odd bottom mode mzis
                j = i + 1
                if j > end:
                    continue
                instances[f"bs_{last}_{j}"] = f"bs_{last}_{j}"

        # create middle mzis
        for i in range(1, last):
            if i % 2 == 0:
                for j in np.arange(0, dim, 2):
                    instances[f"bs_{i}_{j}"] = f"bs_{i}_{j}"
            else:
                for j in range(1, dim, 2):
                    instances[f"bs_{i}_{j}"] = f"bs_{i}_{j}"

        # Connecting Mesh
        for i in range(dim-1):

            # Top boundary Connections
            if i == 0:
                if dim % 2 == 0:
                    for j in range(0, dim - 2, 2):
                        connections[f"bs_{i}_{j},out0"] = f"wvgd_top_{j+1},in0"
                        connections[f"wvgd_top_{j+1},out0"] = f"bs_{i}_{j+2},in0"
                    connections[f"bs_{i}_{dim-2},out0"] =f"wvgd_top_{end},in0"
                else:
                    for j in range(0, dim - 2, 2):
                        if j < last:
                            connections[f"bs_{i}_{j},out0"] =f"wvgd_top_{j+1},in0"
                        if f"bs_{i}_{j+2}" in instances:
                            connections[f"wvgd_top_{j+1},out0"] = f"bs_{i}_{j+2},in0"

            # Connecting diagonal mzis
            for j in range(dim):
                key = f"bs_{i}_{j}"
                if key not in instances:
                    continue
                #up-right connection
                if i > 0:
                    upper_right = f"bs_{i-1}_{j+1}"
                    if upper_right in instances:
                        connections[f"bs_{i}_{j},out0"] = f"bs_{i-1}_{j+1},in1"
                # downward-right connection
                target = f"bs_{i+1}_{j+1}"
                if target in instances:
                    connections[f"bs_{i}_{j},out1"] = f"bs_{i+1}_{j+1},in0"

            # bottom boundary connections
            if i == last:
                if dim % 2 == 0:
                        # intermediate bottom boundary connections
                    for j in range(1, last, 2):
                        connections[f"bs_{last}_{j-1},out1"] = f"wvgd_bottom_{j},in0"
                        connections[f"wvgd_bottom_{j},out0"] = f"bs_{last}_{j+1},in1"

                    # final bottom boundary segment
                    connections[f"bs_{last}_{last},out1"] = f"wvgd_bottom_{end},in0"
                    # bottom output
                    connections[f"wvgd_bottom_{end},out0"] = f"phaser_o{end},in0"
                    ports[f"out{end}"] = f"phaser_o{end},out0"
                else:
                    # first bottom phaser
                    connections[f"wvgd_bottom_0,out0"] = f"bs_{last}_1,in1"

                    # last bottom phaser
                    connections[f"bs_{last}_{last},out1"] = f"wvgd_bottom_{end},in0"

                    # intermediate bottom phasers
                    for j in range(2, last, 2):
                        connections[f"bs_{i}_{j-1},out1"] = f"wvgd_bottom_{j},in0"
                        connections[f"wvgd_bottom_{j},out0"] = f"bs_{i}_{j+1},in1"

        # Connecting output phase shifters 
        if dim % 2 == 0:
            # top boundary output
            connections[f"wvgd_top_{end},out0"] = "phaser_o0,in0"
            ports["out0"] = "phaser_o0,out0"

            # middle BS outputs
            for i in range(1, last, 2):
                connections[f"bs_{i}_{end},out0"] = f"phaser_o{i},in0"
                connections[f"bs_{i}_{end},out1"] = f"phaser_o{i+1},in0"

                ports[f"out{i}"] = f"phaser_o{i},out0"
                ports[f"out{i+1}"] = f"phaser_o{i+1},out0"

        # Connecting output phase shifters - odd
        else:
            output_col = end
            for i in range(0, dim - 1, 2):
                connections[f"bs_{i}_{output_col},out0"] = f"phaser_o{i},in0"
                connections[f"bs_{i}_{output_col},out1"] = f"phaser_o{i+1},in0"
                ports[f"out{i}"] = f"phaser_o{i},out0"
                ports[f"out{i+1}"] = f"phaser_o{i+1},out0"
            # bottom mode
        connections[f"wvgd_bottom_{end},out0"] = f"phaser_o{end},in0"
        ports[f"out{end}"] = f"phaser_o{end},out0"  

        # input ports
        for i in range(0, dim - 1, 2):
            ports[f"in{i}"] = f"bs_{i}_0,in0"
            ports[f"in{i+1}"] = f"bs_{i}_0,in1"

        if dim % 2 != 0:
            ports[f"in{end}"] = "wvgd_bottom_0,in0"
        netlist = {}
        netlist['instances'] = instances
        netlist['connections'] = connections
        netlist['ports'] = ports
        return netlist
    
    def clements_models_ideal(self,U):
        I = itf.square_decomposition(U)
        dim = U.shape[0]
        theta, phi = self.map_clements_bs(I)
        phases = I.output_phases
        models = {}
        models['wvgd'] = ideal.ideal_waveguide
        for i in range(dim):
            phase = phases[i]
            models[f"phaser_o{i}"] = lambda phase=phase: ideal.ideal_phaseshifter(phase)
            for j in range(dim):
                if (i,j) in phi.keys():
                    phi_=phi[(i,j)]
                    theta_ = theta[(i,j)]
                    models[f"bs_{i}_{j}"] = lambda phi_=phi_, theta_ = theta_: ideal.ideal_clements_tbu_model(theta_, phi_)
        return models

    def map_clements_bs(self,I):
        theta = {}
        phi = {}
        dim = self.dimension
        mode_tracker = np.zeros(dim, dtype=int)

        for BS in I.BS_list:
            m1 = BS.mode1 - 1 
            m2 = BS.mode2 - 1

            # Physical column in Clements mesh
            j = max(mode_tracker[m1], mode_tracker[m2])

            # Physical row = upper mode
            i = m1

            theta[(i, j)] = BS.theta
            phi[(i, j)] = BS.phi

            # Both modes have now advanced one column
            mode_tracker[m1] = j + 1
            mode_tracker[m2] = j + 1

        return theta, phi