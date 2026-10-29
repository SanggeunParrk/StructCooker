import numpy as np

from structcooker.fragdb.extract import virtual_cb


def test_virtual_cb_matches_the_af2_ideal_geometry():
    n, ca, c = np.array([[-0.525, 1.363, 0.0]]), np.zeros((1, 3)), np.array([[1.526, 0.0, 0.0]])
    assert np.allclose(virtual_cb(n, ca, c), [[-0.529, -0.774, -1.205]], atol=0.01)
