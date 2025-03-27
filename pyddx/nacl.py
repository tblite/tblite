import pyddx
import numpy as np

tobohr = 1 / 0.52917721092
charges = np.array([
       -0.89373133785716652, 0.89373133785716852 ])
rvdw = np.array([
     3.8739385554719026, 3.4015070243167926 ])
centres = tobohr * np.array([
    [ 0.4000000000000, -0.200000000000,  0.100000000000],
    [ 0.3000000000000,  0.260000000000,  1.300000000000],
]).T

print(pyddx.banner())

model = pyddx.Model("pcm", centres, rvdw, solvent_epsilon=78.359999999999999, n_lebedev=302, eta=0.10000, lmax=1)

# Compute solute contributions (here just charges)
solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
solute_field = model.multipole_electrostatics(solute_multipoles)
solute_psi = model.multipole_psi(solute_multipoles)

print(f'solute_psi: {solute_psi}')
# Solve the problem
state = pyddx.State(model, solute_psi, solute_field["phi"])
state.fill_guess()
state.solve()
state.fill_guess_adjoint()
state.solve_adjoint()

# Show results
#keps = 0.98097895003804214 # COSMO water
keps = 1.00000000000000000 #PCM
energy = keps * 0.5000000000000 * np.sum(state.x * solute_psi)
force = state.solvation_force_terms(solute_field)
force += state.multipole_force_terms(solute_multipoles)
print(energy)
print(force)
#print(f'xs: {state.x}')
