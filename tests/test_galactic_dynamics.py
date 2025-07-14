import numpy as np
from src.galactic_dynamics import v_rotation_UDT

def test_v_rotation_UDT():
    r = np.array([0.42, 1.26])  # Sample kpc
    M_star = 5e7  # Solar masses
    v_pred = v_rotation_UDT(r, M_star)
    assert np.allclose(v_pred[:2], [14.3, 28.5], rtol=0.01)  # Matches example

if __name__ == '__main__':
    test_v_rotation_UDT()
    print("Galactic dynamics test passed")