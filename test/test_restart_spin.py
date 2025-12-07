import pytest
import numpy as np
from tblite.interface import Calculator, Result

def test_spin_polarized_restart():
    # Use a small open-shell system (e.g., NO) in Bohr
    numbers = np.array([7, 8], dtype=np.int32)
    positions = np.array([
        [0.0, 0.0, -1.1],
        [0.0, 0.0,  1.1],
    ], dtype=float)

    # Capture output to count iterations
    log = []
    def logger(msg):
        log.append(msg)

    # Enable spin polarization via interaction
    calc = Calculator("GFN2-xTB", numbers, positions, uhf=1, logger=logger)
    calc.add("spin-polarization", 1.0)
    calc.set("verbosity", 1)
    
    # Full run
    res_full = calc.singlepoint()
    energy_full = res_full.get("energy")
    
    # Clear log for restart run
    log.clear()

    # Restart using full wavefunction (copy Result)
    res_restart_wfn = Result(res_full)
    res_restart_wfn = calc.singlepoint(res_restart_wfn)
    energy_restart_wfn = res_restart_wfn.get("energy")
    
    # Count iterations (lines starting with an integer index)
    iterations = 0
    for line in log:
        parts = line.strip().split()
        if len(parts) > 0 and parts[0].isdigit():
            iterations += 1
            
    # Check energy consistency
    assert abs(energy_full - energy_restart_wfn) < 1e-8
    
    # Verify iteration count (should be very low, e.g., <= 3)
    # Ideally 1 or 2 if restart is perfect
    assert iterations <= 2, f"Restart took {iterations} iterations, expected <= 3"
