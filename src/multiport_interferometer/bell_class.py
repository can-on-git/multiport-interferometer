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

from . import bell_components
from . import bell_decomposer
from . import functions as f
import sax
import numpy as np
try:
    from . import bell_models
except ImportError:
    bell_models=None
import gdsfactory as gf

class Bell:
    def __init__(self, dimension):
        if not isinstance(dimension, int) or dimension < 2:
            raise ValueError("Dimension must be an integer >= 2")
        self.dimension = dimension
        self.components = self.load_components() #dictionaro of components
        self.gds_component, self.sax_netlist = self.build()
        if bell_models is not None:
            self.models = self.load_models()
            self.phase_offset = bell_models.phase_offset
            self.phase_step = bell_models.phase_step
        else:
            self.models = None
            print("No models.py found. SAX models were not loaded.")

    def gds(self):
        return self.gds_component
        
    def load_components(self):
        return {"tbu":bell_components.tbu,
                "phase_shifter": bell_components.phase_shifter,
                "edge_phase_shifter": bell_components.edge_phase_shifter,
                "waveguide": bell_components.plain_wvgd,
                "edge_waveguide":bell_components.edge_wvgd
                }

    def load_models(self):
        return {"tbu_model":bell_models.tbu_model,
                "phase_shifter_model": bell_models.phase_shifter_model,
                "edge_phase_shifter_model": bell_models.edge_phase_shifter_model,
                "plain_wvgd_model": bell_models.plain_wvgd_model}

    def emulate(self, target):
        if np.shape(target)[0] != self.dimension:
            raise ValueError("Target matrix incompatible with Interferometer Dimension")
        decom = bell_decomposer.decompose_bell(target)
        self.values = f.decom_to_model(self.models,decom)
        circuit, info = sax.circuit(netlist=self.sax_netlist, models = self.values)
        S = circuit()
        return S
    def phase_corrected(self, S):
        correction = self.phase_offset + (self.dimension-2)*self.phase_step
        Smat = f.sdict_to_matrix(S)*np.exp(-1j*correction)
        return Smat

    def build(self):
        bell_components._check_pdk()
        dim = self.dimension
        c = gf.Component()
        mzi_cell = self.components["tbu"]()
        phaser_cell = self.components["phase_shifter"]()
        edge_phaser_cell = self.components["edge_phase_shifter"]()
        plain_wvgd = self.components["waveguide"]()
        edge_wvgd = self.components["edge_waveguide"]()

        last = dim-2
        end = dim-1

        instances = {}
        connections = {}
        ports = {}

        phase_i = [None]*dim
        phase_o = [None]*dim
        phasers = {}
        mzms = {}
        wvgds = {}

        # creating top boundary connector waveguides
        for j in range(1, dim, 2): 
            wvgds[(j)] = c<<plain_wvgd
            instances[f"wvgd_{j}"] = "wvgd"

        # creating bottom boundary phase shifters
        if dim % 2 == 0: 
            for j in range(1, dim, 2):
                phasers[(end, j)] = c << phaser_cell
                instances[f"phaser_{end}_{j}"] = f"phaser_{end}_{j}"
        else:
            for j in range(0, dim, 2):
                phasers[(end, j)] = c << phaser_cell
                instances[f"phaser_{end}_{j}"] = f"phaser_{end}_{j}"

        # creating input and output edge phase shifters
        if dim%2!=0:
            for i in range(0,dim,2):
                if i>0:
                    phaser_o = c << edge_phaser_cell
                    phaser_o.mirror_x()
                    phase_o[i-1] = phaser_o
                    phase_i[i-1] = c<<edge_wvgd
                    instances[f"phaser_o{i-1}"] = f"phaser_o{i-1}"
                    instances[f"phaser_i{i-1}"] = f"phaser_i{i-1}"
                if i<dim-1:
                    phaser_i = c << edge_phaser_cell
                    phase_i[i] = phaser_i
                    phase_o[i] = c<<edge_wvgd
                    instances[f"phaser_i{i}"] = f"phaser_i{i}"
                    instances[f"phaser_o{i}"] = f"phaser_o{i}"
                else:
                    phase_i[i] = c<<edge_wvgd
                    phase_o[i] = c<<edge_wvgd
                    instances[f"phaser_i{i}"] = f"phaser_i{i}"
                    instances[f"phaser_o{i}"] = f"phaser_o{i}"
        else:
            for i in range(dim):
                if i%2==0:
                    if i<dim:
                        phaser_i =c<<edge_phaser_cell
                        phase_i[i] = phaser_i
                        phaser_o = c<<edge_phaser_cell
                        phase_o[i] = phaser_o
                        instances[f"phaser_i{i}"] = f"phaser_i{i}"
                        instances[f"phaser_o{i}"] = f"phaser_o{i}"
                    else:
                        phase_i[i] = c<<edge_wvgd
                        phase_o[i] = c<<edge_wvgd
                        instances[f"phaser_i{i}"] = f"phaser_i{i}"
                        instances[f"phaser_o{i}"] = f"phaser_o{i}"
                else:
                    phase_i[i] = c<<edge_wvgd
                    phase_o[i] = c<<edge_wvgd
                    instances[f"phaser_i{i}"] = f"phaser_i{i}"
                    instances[f"phaser_o{i}"] = f"phaser_o{i}"

        # creating top and bottom mzis
        for i in range(0, dim, 2):
            # top mode mzis
            mzm_top = c << mzi_cell
            mzms[f"(0_{i})"] = mzm_top
            instances[f"mzm_0_{i}"] = f"mzm_0_{i}"
            if dim % 2 == 0:    # even bottom mode mzis
                mzm_bottom = c << mzi_cell
                mzms[f"({last}_{i})"] = mzm_bottom
                instances[f"mzm_{last}_{i}"] = f"mzm_{last}_{i}"
            else:   # odd bottom mode mzis
                j = i + 1
                if j > end:
                    continue
                mzm_bottom = c << mzi_cell
                mzms[f"({last}_{j})"] = mzm_bottom
                instances[f"mzm_{last}_{j}"] = f"mzm_{last}_{j}"

        # create middle mzis
        for i in range(1, last):
            if i % 2 == 0:
                for j in np.arange(0, dim, 2):
                    mzm = c << mzi_cell
                    mzms[f"({i}_{j})"] = mzm
                    instances[f"mzm_{i}_{j}"] = f"mzm_{i}_{j}"
            else:
                for j in range(1, dim, 2):
                    mzm = c << mzi_cell
                    mzms[f"({i}_{j})"] = mzm
                    instances[f"mzm_{i}_{j}"] = f"mzm_{i}_{j}"

        # Connecting Mesh
        for i in range(dim):

            # Top boundary Connections
            if i == 0:
                if dim % 2 == 0:
                    for j in range(0, dim - 2, 2):
                        p = wvgds[(j+1)]
                        p.connect("o1",mzms[f"({i}_{j})"].ports["o3"])
                        mzms[f"({i}_{j+2})"].connect("o2",p.ports["o2"])
                        connections[f"mzm_{i}_{j},out0"] = f"wvgd_{j+1},in0"
                        connections[f"wvgd_{j+1},out0"] = f"mzm_{i}_{j+2},in0"
                    p = wvgds[(end)]
                    p.connect("o1",mzms[f"({i}_{dim-2})"].ports["o3"])
                    connections[f"mzm_{i}_{dim-2},out0"] =f"wvgd_{end},in0"
                else:
                    for j in range(0, dim - 2, 2):
                        p = wvgds[(j+1)]
                        if j < last:
                            p.connect("o1",mzms[f"({i}_{j})"].ports["o3"])
                            connections[f"mzm_{i}_{j},out0"] =f"wvgd_{j+1},in0"
                        if f"({i}_{j+2})" in mzms:
                            mzms[f"({i}_{j+2})"].connect("o2",p.ports["o2"])
                            connections[f"wvgd_{j+1},out0"] = f"mzm_{i}_{j+2},in0"

            # Connecting diagonal mzis
            for j in range(dim):
                key = f"({i}_{j})"
                if key not in mzms:
                    continue
                #up-right connection
                if i > 0:
                    upper_right = f"({i-1}_{j+1})"
                    if upper_right in mzms:
                        mzms[key].connect("o3",mzms[upper_right].ports["o1"])
                        connections[f"mzm_{i}_{j},out0"] = f"mzm_{i-1}_{j+1},in1"
                # downward-right connection
                target = f"({i+1}_{j+1})"
                if target in mzms:
                    mzms[target].connect("o2", mzms[key].ports["o4"])
                    connections[f"mzm_{i}_{j},out1"] = f"mzm_{i+1}_{j+1},in0"

            # bottom boundary connections
            if i == last:
                if dim % 2 == 0:
                    for j in range(1, last, 2):
                        p = phasers[(end, j)]
                        p.connect("o1", mzms[f"({i}_{j-1})"].ports["o4"])
                        p.connect("o2", mzms[f"({i}_{j+1})"].ports["o1"])
                        connections[f"mzm_{i}_{j-1},out1"] = f"phaser_{end}_{j},in0"
                        connections[f"phaser_{end}_{j},out0"] = f"mzm_{i}_{j+1},in1"
                    # last bottom phaser
                    p = phasers[(end, end)]
                    p.connect("o1", mzms[f"({last}_{last})"].ports["o4"])
                    connections[f"mzm_{last}_{last},out1"] = f"phaser_{end}_{end},in0"
                else:
                    # first bottom phaser
                    p = phasers[(end, 0)]
                    p.connect("o2", mzms[f"({last}_1)"].ports["o1"])
                    connections[f"phaser_{end}_0,out0"] = f"mzm_{last}_1,in1"

                    # last bottom phaser
                    p = phasers[(end, end)]
                    p.connect("o1",mzms[f"({last}_{last})"].ports["o4"])
                    connections[f"mzm_{last}_{last},out1"] = f"phaser_{end}_{end},in0"

                    # intermediate bottom phasers
                    for j in range(2, last, 2):
                        p = phasers[(end, j)]
                        p.connect("o1",mzms[f"({i}_{j-1})"].ports["o4"])
                        connections[f"mzm_{i}_{j-1},out1"] = f"phaser_{end}_{j},in0"
                        p.connect("o2",mzms[f"({i}_{j+1})"].ports["o1"])
                        connections[f"phaser_{end}_{j},out0"] = f"mzm_{i}_{j+1},in1"

        # connecting input phaseshifters
        for i in range(0, dim - 1, 2):
            phase_i[i].connect("o2",mzms[f"({i}_0)"].ports["o2"])
            connections[f"phaser_i{i},out0"] = f"mzm_{i}_0,in0" #these are sax connections, port naming is swapped with the gds naming because sax component models using standard matrix convention
            phase_i[i+1].connect("o2",mzms[f"({i}_0)"].ports["o1"])
            connections[f"phaser_i{i+1},out0"] = f"mzm_{i}_0,in1"
            c.add_port(f"in{i}",port=phase_i[i].ports["o1"])
            ports[f"in{i}"] = f"phaser_i{i},in0"
            c.add_port(f"in{i+1}",port=phase_i[i+1].ports["o1"])
            ports[f"in{i+1}"] = f"phaser_i{i+1},in0"

        # odd final input phaseshifter
        if dim % 2 != 0:
            phase_i[end].connect("o2",phasers[(end, 0)].ports["o1"])
            connections[f"phaser_i{end},out0"] = f"phaser_{end}_0,in0"
            c.add_port(f"in{end}",port=phase_i[end].ports["o1"])
            ports[f"in{end}"] = f"phaser_i{end},in0"

        # Connecting output phase shifters - Even
        if dim % 2 == 0:
            # top mode
            phase_o[0].connect("o2",wvgds[(end)].ports["o2"])
            connections[f"wvgd_{end},out0"] = "phaser_o0,in0"
            c.add_port("out0",port=phase_o[0].ports["o1"])
            ports["out0"] = "phaser_o0,out0"

            # middle modes
            for i in range(1, end, 2):
                phase_o[i].connect("o2",mzms[f"({i}_{end})"].ports["o3"])
                phase_o[i+1].connect("o2",mzms[f"({i}_{end})"].ports["o4"])
                connections[f"mzm_{i}_{end},out0"] = f"phaser_o{i},in0"
                connections[f"mzm_{i}_{end},out1"] = f"phaser_o{i+1},in0"
                c.add_port(f"out{i}",port=phase_o[i].ports["o1"])
                c.add_port(f"out{i+1}",port=phase_o[i+1].ports["o1"])
                ports[f"out{i}"] = f"phaser_o{i},out0"
                ports[f"out{i+1}"] = f"phaser_o{i+1},out0"

            # bottom mode
            phase_o[end].connect("o2", phasers[(end, end)].ports["o2"])
            connections[f"phaser_{end}_{end},out0"] = f"phaser_o{end},in0"
            c.add_port(f"out{end}", port=phase_o[end].ports["o1"])
            ports[f"out{end}"] = f"phaser_o{end},out0"

        # Connecting output phase shifters - Even
        else:
            output_col = end
            for i in range(0, dim - 1, 2):
                phase_o[i].connect("o2", mzms[f"({i}_{output_col})"].ports["o3"])
                connections[f"mzm_{i}_{output_col},out0"] = f"phaser_o{i},in0"
                phase_o[i+1].connect("o2", mzms[f"({i}_{output_col})"].ports["o4"])
                connections[f"mzm_{i}_{output_col},out1"] = f"phaser_o{i+1},in0"
                c.add_port(f"out{i}", port=phase_o[i].ports["o1"])
                ports[f"out{i}"] = f"phaser_o{i},out0"
                c.add_port(f"out{i+1}", port=phase_o[i+1].ports["o1"])
                ports[f"out{i+1}"] = f"phaser_o{i+1},out0"
                
            # bottom mode
            phase_o[end].connect("o2", phasers[(end, end)].ports["o2"])
            connections[f"phaser_{end}_{end},out0"] = f"phaser_o{end},in0"
            c.add_port(f"out{end}", port=phase_o[end].ports["o1"])
            ports[f"out{end}"] = f"phaser_o{end},out0"
        
        netlist = {}
        netlist["instances"]=instances
        netlist["connections"]=connections
        netlist["ports"]=ports
        return c, netlist