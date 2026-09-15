'''
File of my custom gdsfactory component functions to be used in the project
'''

import gdsfactory as gf
import uuid

# I use the Si 220 cband pdk from cornerstone
from cspdk.si220.cband import cells, PDK, LAYER_STACK

def activate_pdk():
    PDK.activate()

def _check_pdk():
    if gf.get_active_pdk().name != PDK.name:
        raise RuntimeError("Activate the Cornerstone Si220 C-band PDK first.")

'''
The circuit requires the 4 component names: 
'''
def tbu():
    return mzi()
def phase_shifter():
    return phaser()
def edge_phase_shifter():
    return edge_phaser()
def plain_wvgd():
    return plain()
#def edge_wvgd(): already exists with correct name

'''
the actual gdscomponents:
'''

phaser_len = 150

@gf.cell
def s_bend():
    c = gf.Component()
    bend = cells.bend_euler()
    up = c << bend
    down = c<< bend
    up.rotate(90)
    down.rotate(-90)
    down.connect("o1", up.ports["o1"])
    c.add_port("o1", port=up.ports["o2"])
    c.add_port("o2", port = down.ports["o2"])
    return c

@gf.cell 
def mzi():
    c = gf.Component()
    coupler = cells.coupler()
    phaser = cells.straight_heater_metal(length = phaser_len)
    bend = s_bend()
    i_up = c<<bend
    i_down = c<<bend
    i_down.mirror()
    o_up = c << bend
    o_down = c<<bend
    o_up.mirror()
    top = c<<phaser
    bottom = c<<phaser
    top_left= c<<bend
    top_left.mirror()
    top_right = c<<bend
    bottom_left = c<<bend
    bottom_right = c<<bend
    bottom_right.mirror()
    coupler1 = c<<coupler
    coupler2 = c<<coupler
    i_down.connect("o1", coupler1.ports["o1"])
    i_up.connect("o1", coupler1.ports["o2"])
    top_left.connect("o1", coupler1.ports["o3"])
    bottom_left.connect("o1", coupler1.ports["o4"])
    top.connect("o1", top_left["o2"])
    bottom.connect("o1", bottom_left["o2"])
    top_right.connect("o1", top['o2'])
    bottom_right.connect("o1", bottom['o2'])
    coupler2.connect("o1", bottom_right["o2"])
    coupler2.connect("o2", top_right["o2"])
    o_up.connect('o1', coupler2["o3"])
    o_down.connect('o1', coupler2["o4"])
    c.add_port(name="o1", port=i_down['o2'])
    c.add_port(name = 'o2', port=i_up['o2'])
    c.add_port(name = 'o3', port=o_up['o2'])
    c.add_port(name = 'o4', port=o_down['o2'])
    return c

@gf.cell
def edge_phaser(phase_len=phaser_len):
    c = gf.Component()
    wvgd = cells.straight(cross_section="strip")
    phaser = cells.straight_heater_metal(length= phase_len)
    w = c<<wvgd
    p = c<<phaser
    p.connect("o1", w["o2"])
    c.add_port(name="o1", port = w["o1"])
    c.add_port(name = "o2", port=p["o2"])
    return c

@gf.cell
def edge_wvgd(phase_len=phaser_len):
    c = gf.Component()
    wvgd1 = cells.straight(cross_section="strip")
    wvgd2 = cells.straight(length= phase_len)
    w1 = c<<wvgd1
    w2 = c<<wvgd2
    w2.connect("o1", w1["o2"])
    c.add_port(name="o1", port = w1["o1"])
    c.add_port(name = "o2", port=w2["o2"])
    return c

def coupler_single_arm(keep="top"):
    c = gf.Component(f"coupler_single_arm_{keep}_{uuid.uuid4().hex[:6]}")

    cp = cells.coupler()
    cp = cp.copy()
    cp.flatten()

    WG_LAYER = gf.get_active_pdk().layers.WG
    wg = cp.extract(layers=[WG_LAYER])

    bbox = wg.bbox()
    xmin = bbox.left
    xmax = bbox.right
    ymin = bbox.bottom
    ymax = bbox.top
    ymid = (ymin + ymax) / 2

    if keep == "top":
        y0 = ymid
        height = ymax - ymid
    else:
        y0 = ymin
        height = ymid - ymin

    clip = gf.components.rectangle(
        size=(xmax - xmin, height),
        layer=WG_LAYER,
    )

    clip_component = gf.Component(f"clip_{uuid.uuid4().hex[:6]}")
    r = clip_component << clip
    r.move((xmin, y0))

    arm = gf.boolean(
        A=wg,
        B=clip_component,
        operation="and",
        layer=WG_LAYER,
    )

    c << arm

    cp_original = cells.coupler()
    cp_ports = list(cp_original.ports)

    if keep == "top":
        selected_ports = [p for p in cp_ports if p.center[1] > ymid]
    else:
        selected_ports = [p for p in cp_ports if p.center[1] < ymid]

    selected_ports = sorted(selected_ports, key=lambda p: p.center[0])

    c.add_port(name="o1", port=selected_ports[0])
    c.add_port(name="o2", port=selected_ports[1])

    return c

@gf.cell
def phaser(phase_len=phaser_len):
    c = gf.Component()
    bottom = coupler_single_arm("bottom")
    bend = s_bend()
    phase_shift = cells.straight_heater_metal(length =phase_len)
    i = c<<bend
    o = c<<bend
    o.mirror()
    left = c<<bend
    left.mirror()
    right = c<<bend
    s1 = c<<bottom
    s2 = c<<bottom
    phaser = c<<phase_shift
    s1.connect("o1", i["o2"])
    left.connect("o1", s1['o2'])
    phaser.connect("o1", left['o2'])
    right.connect('o1', phaser['o2'])
    s2.connect('o1', right['o2'])
    o.connect('o1', s2['o2'])
    c.add_port(name='o1', port=i['o1'])
    c.add_port(name = 'o2', port=o['o2'])
    return c

@gf.cell
def plain(lent=phaser_len):
    c = gf.Component()
    bottom = coupler_single_arm("bottom")
    bend = s_bend()
    straight = cells.straight(length =lent)
    i = c<<bend
    o = c<<bend
    o.mirror()
    left = c<<bend
    left.mirror()
    right = c<<bend
    s1 = c<<bottom
    s2 = c<<bottom
    wvgd = c<<straight
    s1.connect("o1", i["o2"])
    left.connect("o1", s1['o2'])
    wvgd.connect("o1", left['o2'])
    right.connect('o1', wvgd['o2'])
    s2.connect('o1', right['o2'])
    o.connect('o1', s2['o2'])
    c.add_port(name='o1', port=i['o1'])
    c.add_port(name = 'o2', port=o['o2'])
    return c
