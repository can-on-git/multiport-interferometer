import numpy as np
import pytest
import sax
from scipy.stats import unitary_group

from multiport_interferometer import Bell, Clements, decompose_bell, sdict_to_matrix, bell_circuit_schematic, ideal_bell_models
from multiport_interferometer.bell_components import activate_pdk

activate_pdk()
@pytest.mark.parametrize("dim", [2, 3, 4, 5, 10])
def test_decomposition(dim):
    """Decomposed ideal Bell mesh should reproduce the target unitary."""
    U = unitary_group.rvs(dim)

    decom = decompose_bell(U)
    netlist = bell_circuit_schematic(dim)
    models = ideal_bell_models(decom)

    circuit, _ = sax.circuit(netlist=netlist, models=models)
    result = sdict_to_matrix(circuit())

    assert np.allclose(result, U, atol=1e-10)


@pytest.mark.parametrize("dim", [2, 3, 4, 5])
def test_bell_build(dim):
    """Bell should build with the correct external ports."""
    bell = Bell(dim)

    expected = (
        {f"in{i}" for i in range(dim)}
        | {f"out{i}" for i in range(dim)}
    )

    assert set(bell.sax_netlist["ports"]) == expected


def test_sdict_to_matrix():
    """SAX dictionary should be converted with correct port ordering."""
    S = {
        ("in0", "out0"): np.array([1]),
        ("in0", "out1"): np.array([2]),
        ("in1", "out0"): np.array([3]),
        ("in1", "out1"): np.array([4]),
    }

    result = sdict_to_matrix(S)

    expected = np.array([
        [1, 3],
        [2, 4]
    ])

    assert np.array_equal(result, expected)


def test_invalid_inputs():
    """Invalid Bell dimensions and non-unitary matrices should fail."""
    with pytest.raises(ValueError):
        Bell(1)

    with pytest.raises(ValueError):
        decompose_bell(np.ones((3, 3)))

def test_bell_emulate():
    dim = 3
    bell = Bell(dim)
    U = unitary_group.rvs(dim)

    result = sdict_to_matrix(bell.emulate(U))

    assert result.shape == (dim, dim)
    assert np.all(np.isfinite(result))

def test_invalid_inputs_clements():
    """Invalid Clements dimensions should fail."""
    with pytest.raises(ValueError):
        Clements(1)

    with pytest.raises(ValueError):
        Clements(2).ideal_emulate(np.ones((4, 3)))

@pytest.mark.parametrize("dim", [2, 3, 4, 5])
def test_clements_emulate(dim):
    mesh = Clements(dim)
    U = unitary_group.rvs(dim)

    result = sdict_to_matrix(mesh.ideal_emulate(U))

    assert result.shape == (dim, dim)
    assert np.all(np.isfinite(result))