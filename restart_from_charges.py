import numpy as np

from tblite.interface import Calculator, Result


def run_lithium_hydride():
    # LiH in Bohr: approximate bond length 1.595 Å = 3.013 bohr
    numbers = np.array([3, 1], dtype=np.int32)  # Li, H
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.5],
    ], dtype=float)

    numbers = np.array([6, 6, 7, 7, 1, 1, 1, 1, 1, 1, 8, 8,])
    positions = np.array([  
                    [-3.81469488143921, +0.09993441402912, 0.00000000000000],
                    [+3.81469488143921, -0.09993441402912, 0.00000000000000],
                    [-2.66030049324036, -2.15898251533508, 0.00000000000000],
                    [+2.66030049324036, +2.15898251533508, 0.00000000000000],
                    [-0.73178529739380, -2.28237795829773, 0.00000000000000],
                    [-5.89039325714111, -0.02589114569128, 0.00000000000000],
                    [-3.71254944801331, -3.73605775833130, 0.00000000000000],
                    [+3.71254944801331, +3.73605775833130, 0.00000000000000],
                    [+0.73178529739380, +2.28237795829773, 0.00000000000000],
                    [+5.89039325714111, +0.02589114569128, 0.00000000000000],
                    [-2.74426102638245, +2.16115570068359, 0.00000000000000],
                    [+2.74426102638245, -2.16115570068359, 0.00000000000000],
                    ]) 

    calc = Calculator("GFN2-xTB", numbers, positions)

    # Full run
    res_full = calc.singlepoint()
    energy_full = res_full.get("energy")
    forces_full = -res_full.get("gradient")
    charges = res_full.get("charges")
    qsh = res_full.get("shell-charges")
    dpat = res_full.get("atomic-dipoles")  # shape (nat, 3)
    qmat = res_full.get("atomic-quadrupoles")  # shape (nat, 6)

    # Restart using shell charges plus moments guess (qsh, dpat, qmat)
    res_restart_qsh = Result()
    res_restart_qsh.set("shell-charges-and-moments-guess", (qsh, dpat, qmat))
    res_restart_qsh = calc.singlepoint(res_restart_qsh)
    energy_restart_qsh = res_restart_qsh.get("energy")
    forces_restart_qsh = -res_restart_qsh.get("gradient")

    e_diff_qsh = abs(float(energy_full) - float(energy_restart_qsh))
    f_diff_qsh = np.max(np.abs(forces_full - forces_restart_qsh))
    
    print("\nComparison between full calculation and shell-charges restart:")
    print("Energy difference:", e_diff_qsh)
    print("Maximum force difference:", f_diff_qsh)

    # True restart using full wavefunction (copy Result)
    res_restart_wfn = Result(res_full)
    res_restart_wfn = calc.singlepoint(res_restart_wfn)
    energy_restart_wfn = res_restart_wfn.get("energy")
    forces_restart_wfn = -res_restart_wfn.get("gradient")

    e_diff_wfn = abs(float(energy_full) - float(energy_restart_wfn))
    f_diff_wfn = np.max(np.abs(forces_full - forces_restart_wfn))
    
    print("\nComparison between full calculation and wavefunction restart:")
    print("Energy difference:", e_diff_wfn)
    print("Maximum force difference:", f_diff_wfn)


def run_spin_polarized():
    # Use a small open-shell system (e.g., NO) in Bohr
    numbers = np.array([7, 8], dtype=np.int32)
    positions = np.array([
        [0.0, 0.0, -1.1],
        [0.0, 0.0,  1.1],
    ], dtype=float)

    # Enable spin polarization via interaction
    calc = Calculator("GFN2-xTB", numbers, positions, uhf=1)
    calc.add("spin-polarization", 1.0)

    # Full run
    res_full = calc.singlepoint()
    energy_full = res_full.get("energy")
    forces_full = -res_full.get("gradient")
    qsh = res_full.get("shell-charges")  # shape (2, nsh)
    dpat = res_full.get("atomic-dipoles")  # shape (2, nat, 3)
    qmat = res_full.get("atomic-quadrupoles")  # shape (2, nat, 6)

    # Validate shapes
    assert qsh.ndim == 2 and qsh.shape[0] == 2
    assert dpat.ndim == 3 and dpat.shape[0] == 2 and dpat.shape[2] == 3
    assert qmat.ndim == 3 and qmat.shape[0] == 2 and qmat.shape[2] == 6

    # Restart using shell charges plus moments guess (spin-resolved)
    res_restart_qsh = Result()
    res_restart_qsh.set("shell-charges-and-moments-guess", (qsh, dpat, qmat))
    res_restart_qsh = calc.singlepoint(res_restart_qsh)
    energy_restart_qsh = res_restart_qsh.get("energy")
    forces_restart_qsh = -res_restart_qsh.get("gradient")

    e_diff_qsh = abs(float(energy_full) - float(energy_restart_qsh))
    f_diff_qsh = np.max(np.abs(forces_full - forces_restart_qsh))

    print("\n[Spin] Comparison between full calculation and shell-charges restart:")
    print("Energy difference:", e_diff_qsh)
    print("Maximum force difference:", f_diff_qsh)


if __name__ == "__main__":
    run_lithium_hydride()
    run_spin_polarized()


