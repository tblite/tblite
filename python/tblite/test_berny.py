import logging
from typing import Sequence

import pytest

try:
    import berny
    from tblite.berny import Solver
except ModuleNotFoundError:
    berny = None


@pytest.fixture
def symbols() -> Sequence[str]:
    return [
        "C",
        "H",
        "O",
        "C",
        "H",
        "H",
        "C",
        "H",
        "H",
        "H",
        "H",
        "C",
        "H",
        "H",
        "O",
        "H",
    ]


@pytest.fixture
def positions() -> Sequence[Sequence[float]]:
    return [
        [+1.578385, +0.147690, +0.343809],
        [+1.394750, +0.012968, +1.413545],
        [+1.359929, -1.086203, -0.359782],
        [+0.653845, +1.215099, -0.221322],
        [+1.057827, +2.180283, +0.093924],
        [+0.729693, +1.184864, -1.311438],
        [-0.817334, +1.152127, +0.208156],
        [-1.303525, +2.065738, -0.145828],
        [-0.883765, +1.159762, +1.299260],
        [+1.984120, -1.734446, -0.021385],
        [+2.616286, +0.458948, +0.206544],
        [-1.627725, -0.034052, -0.311301],
        [-2.684229, +0.151015, -0.118566],
        [-1.501868, -0.118146, -1.397506],
        [-1.324262, -1.260154, +0.333377],
        [-0.417651, -1.475314, +0.076637],
    ]


@pytest.mark.skipif(berny is None, reason="requires pyberny")
def test_berny(symbols: Sequence[str], positions: Sequence[Sequence[float]]) -> None:
    opt = berny.Berny(
        berny.Geometry(
            species=symbols,
            coords=positions,
        ),
        maxsteps=10,
        logger=logging.getLogger("test"),
    )
    berny.optimize(opt, Solver("GFN2-xTB"))
    assert opt.converged
