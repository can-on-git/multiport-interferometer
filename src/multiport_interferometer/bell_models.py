'''
Script for Sax Models used
'''

import numpy as np
from pathlib import Path
import sax
import pickle

MODEL_DIR = Path(__file__).resolve().parent / "sax_models"

'''
calibration values from calibration notebook
'''
phase_offset = 3.5650791015113166
phase_step = -1.54018263

'''
the circuit requires models with following names:
'''
def tbu_model(sigma, delta):
    return calibrated_mzi_model(sig=sigma, delt = delta)
# phase_shifter_model(): has correct name
def edge_phase_shifter_model(phi):
    return edge_phaser_model(phi)
#def plain_wvgd_model(): has correct name

############################################################################################################################
'''
Loading Sax dicts
'''
#############################
'''
phi must be in radians
default wavegudie is the 150um waveguide
returns sax matrix of phase with total phase= standard waveguide phase + added thermooptic
'''
def phaser_model(phi, waveguide='waveguide_150.pkl'):
    filepath = MODEL_DIR / waveguide
    with filepath.open("rb") as f:
        ref_sdict = pickle.load(f)
    ref_phase = np.angle(ref_sdict[('in0','out0')])
    new = np.exp(1j*(ref_phase+phi))
    phaser_model = sax.reciprocal({("in0", "out0"): new})
    return phaser_model

def arm_model():
    filepath = MODEL_DIR / 's_arm.pkl'
    with filepath.open("rb") as f:
        arm_model = pickle.load(f)
    return arm_model

def coupler_model():
    filepath = MODEL_DIR / "coupler_50_50_2port.pkl"
    with filepath.open("rb") as f:
        sax_sdict = pickle.load(f)
    return sax_sdict

def edge_wvgd_model():
    filepath = MODEL_DIR / "straight_L10.pkl"
    with filepath.open("rb") as f:
        wvgd_sdict = pickle.load(f)
    return wvgd_sdict

def wvgd_model(waveguide='waveguide_150.pkl'):
    filepath = MODEL_DIR / waveguide
    with filepath.open("rb") as f:
        ref_sdict = pickle.load(f)
    return ref_sdict

def s_bend_model():
    filepath = MODEL_DIR / 's_bend.pkl'
    with filepath.open("rb") as f:
        arm_model = pickle.load(f)
    return arm_model

############################################################################################################################
'''
Plain waveguide Circuit
'''
##########################
def plain_wvgd_model():
    total_circ, info = sax.circuit(
    netlist={
        "instances": {
            "lft": "connector",
            "middle": "straight",
            "rgt": "connector",
        },
        "connections": {
            "lft,out0": "middle,in0",
            "middle,out0": "rgt,in0",
        },
        "ports": {
            "in0": "lft,in0",
            "out0": "rgt,out0",
        },
    },
    models={
        "connector": connector_model,
        "straight": wvgd_model,
    },
    )
    return total_circ()
############################################################################################################################
'''
Phase shifters Circuits
'''
##########################
def connector_model():
    total_circ, info = sax.circuit(
    netlist={
        "instances": {
            "lft": "sbend",
            "middle": "arm",
            "rgt": "sbend",
        },
        "connections": {
            "lft,out0": "middle,in0",
            "middle,out0": "rgt,in0", 
        },
        "ports": {
            "in0": "lft,in0",
            "out0": "rgt,out0",
        },
    },
    models={
        "arm": arm_model,
        "sbend": s_bend_model,
    },
    )
    return total_circ()


'''
returns the sax circuit matrix for the complete arm-phase_shifter-arm component that is used in the mesh
'''
def phase_shifter_model(phi=0):
    total_circ, info = sax.circuit(
    netlist={
        "instances": {
            "lft": "connector",
            "middle": "phaser",
            "rgt": "connector",
        },
        "connections": {
            "lft,out0": "middle,in0",
            "middle,out0": "rgt,in0",
        },
        "ports": {
            "in0": "lft,in0",
            "out0": "rgt,out0",
        },
    },
    models={
        "connector": connector_model,
        "phaser": lambda: phaser_model(phi)
    },
    )
    return total_circ()

def edge_phaser_model(phi=0):
    total_circ, info = sax.circuit(
    netlist={
        "instances": {
            "lft": "wvgd",
            'rgt': 'phaser',
        },
        "connections": {
            "lft,out0": "rgt,in0",
        },
        "ports": {
            "in0": "lft,in0",
            "out0": "rgt,out0",
        },
    },
    models={
        "wvgd": edge_wvgd_model,
        "phaser": lambda: phaser_model(phi)#+2.9604398423257847),
    },
    )
    return total_circ()

########################################################################################################################
'''
MZI Model Circuits
'''
##############################

def mzi_coupler():
    total_circ, info = sax.circuit(
        netlist={
            "instances":{
                'tl':"sbend",
                'bl':'sbend',
                'tr':'sbend',
                'br':'sbend',
                'coupler':'coupler',
            },
            "connections":{
                'tl,out0':'coupler,in1',
                'bl,out0':'coupler,in0',
                'coupler,out1':'tr,in0',
                'coupler,out0':'br,in0',
            },
            "ports":{
                'in0':'bl,in0',
                'in1':'tl,in0',
                'out0':'br,out0',
                'out1':'tr,out0',
            },
        },
        models={
            'sbend':s_bend_model,
            'coupler':coupler_model,
        },
    )
    return total_circ()
'''
returns the mzi smatrix, takes input delta and sigma
'''
def mzi_model(sigma=0,delta=0):
    phi1 = sigma-delta
    phi2 = sigma + delta
    total_circ, info = sax.circuit(
    netlist={
        "instances": {
            "lft": "coupler",
            "top": "phaser1",
            "bottom": "phaser2",
            'rgt': 'coupler',
        },
        "connections": {
            "lft,out0": "bottom,in0",
            'lft,out1': "top,in0",
            "top,out0": "rgt,in1",
            "bottom,out0":"rgt,in0"

        },
        "ports": {
            "in0": "lft,in0",
            "in1": "lft,in1",
            "out0": "rgt,out0",
            "out1": "rgt,out1",
        },
    },
    models={
        "coupler": mzi_coupler,
        "phaser1": lambda: phaser_model(phi1),
        "phaser2": lambda: phaser_model(phi2)
    },
    )
    return total_circ()

def calibrated_mzi_model(sig=0, delt=0):
    sigma_offset = 2.03765676
    delta_offset = 3.0629742188152735
    return mzi_model(sigma=sig+sigma_offset,delta=delt+delta_offset)