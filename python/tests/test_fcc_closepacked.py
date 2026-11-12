"""A periodic fcc Cu snapshot at 300 K with four dissociated edge dislocations (b = 1/2[1-10], a 45-degree quadrupole,
30,336 atoms, one frame of the molecular dynamics of the paper): the close-packed route finds four clusters, each a
dislocation with its two Shockley partials and the fault between them, each carrying one vector of the 1/2<110>
family.  Run: python -m pytest tests"""
import os
import numpy as np
from vorotop.burgers import api, io

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
A = 3.615


def _sample():
    step, box, cols = next(io.read_lammps_dump(os.path.join(HERE, "fcc_cu_edge_quadrupole.dump.gz")))
    return io.positions(cols) - box[:, 0], io.box_matrix(box, cols)


def test_fcc_edge_quadrupole():
    crystal = api.Crystal(np.eye(3) * A, [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]], [29] * 4)
    p, H = _sample()
    res = api.analyze(p, H, crystal)
    carrying = [c for c in res["clusters"] if c["vectors"]]
    assert len(carrying) == 4, [c["n_atoms"] for c in carrying]
    for c in carrying:
        assert len(c["vectors"]) == 1
        v = c["vectors"][0]["conventional"]
        assert np.allclose(sorted(np.abs(v)), [0, 0.5, 0.5], atol=1e-6), v
    assert res["crystalline_fraction"] > 0.95


if __name__ == "__main__":
    test_fcc_edge_quadrupole()
    print("ok")
