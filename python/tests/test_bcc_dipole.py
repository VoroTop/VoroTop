"""A bcc Fe screw dipole (b = a[001] along the line, 2,352 atoms, noise 0.02 d): exactly two clusters carry Burgers
vectors, each +-a[001] in the laboratory frame.  Run: python -m pytest tests  (or python tests/test_bcc_dipole.py)"""
import os
import numpy as np
from vorotop.burgers import api, io

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def test_bcc_screw_dipole():
    crystal = api.crystal_from_cif(os.path.join(HERE, "bcc_Fe.cif"))
    assert crystal.kind == "stable" and crystal.gap > 0.5
    step, box, cols = next(io.read_lammps_dump(os.path.join(HERE, "bcc_screw_dipole.dump.gz")))
    res = api.analyze(io.positions(cols) - box[:, 0], io.box_matrix(box, cols), crystal)
    carrying = [c for c in res["clusters"] if c["vectors"]]
    assert len(carrying) == 2, [c["n_atoms"] for c in carrying]
    for c in carrying:
        assert len(c["vectors"]) == 1
        v = c["vectors"][0]
        # exact in the crystal frame: a lattice vector of the <001> family
        assert np.allclose(sorted(np.abs(v["conventional"])), [0, 0, 1], atol=1e-6), v["conventional"]
        # in the laboratory frame, up to the rotation fitted from the noisy sample: along the line, |b| = a
        assert abs(v["length"] - 2.855) < 0.01 and abs(abs(v["lab"][2]) - 2.855) < 0.01, v["lab"]
    assert res["crystalline_fraction"] > 0.9
    assert not res["warnings"], res["warnings"]


if __name__ == "__main__":
    test_bcc_screw_dipole()
    print("ok")
