import os
import warnings
import pytest

with warnings.catch_warnings():
    mm = pytest.importorskip('simtk.openmm')
    app = pytest.importorskip('simtk.openmm.app')
    unit = pytest.importorskip('simtk.unit')
from pygen_structures import sequence_to_mol, load_charmm_dir
from pygen_structures._functions_const import TOPPAR_DIRECTORY

K = unit.kelvin
fs, ps = unit.femtosecond, unit.picosecond
nm = unit.nanometer

RTF, PRM = load_charmm_dir()

def test_sugar_energy():
    raff_coords = [
        [-1.674e-01, -7.140e-02, -2.350e-02],
        [-1.556e-01,  1.360e-02,  4.430e-02],
        [-2.855e-01, -6.030e-02, -9.100e-02],
        [ 6.230e-02, -9.350e-02, -6.070e-02],
        [ 1.188e-01, -1.649e-01, -1.276e-01],
        [-6.600e-02, -8.690e-02, -1.144e-01],
        [-1.764e-01, -2.036e-01,  5.520e-02],
        [-2.323e-01, -1.871e-01,  1.491e-01],
        [-2.480e-01, -2.915e-01, -2.620e-02],
        [-2.593e-01, -3.802e-01,  1.290e-02],
        [-4.000e-02, -2.640e-01,  7.850e-02],
        [-3.580e-02, -3.204e-01,  1.739e-01],
        [-5.800e-03, -3.524e-01, -2.310e-02],
        [-5.430e-02, -3.333e-01, -1.085e-01],
        [ 6.510e-02, -1.522e-01,  7.660e-02],
        [ 1.602e-01, -2.015e-01,  9.770e-02],
        [ 3.520e-02, -6.610e-02,  1.781e-01],
        [ 1.071e-01, -8.800e-03,  2.087e-01],
        [ 1.251e-01,  3.970e-02, -7.420e-02],
        [ 6.110e-02,  1.115e-01, -1.460e-02],
        [ 1.082e-01,  7.330e-02, -1.807e-01],
        [-4.065e-01,  5.930e-02,  6.460e-02],
        [-3.652e-01,  4.730e-02, -6.660e-02],
        [-5.387e-01,  8.800e-02,  7.210e-02],
        [-5.878e-01,  6.400e-03,  1.329e-01],
        [-5.710e-01,  2.168e-01,  1.429e-01],
        [-5.361e-01,  2.177e-01,  2.479e-01],
        [-5.256e-01,  3.046e-01,  9.400e-02],
        [-7.083e-01,  2.399e-01,  1.469e-01],
        [-7.535e-01,  1.724e-01,  2.058e-01],
        [-3.184e-01,  1.794e-01, -1.167e-01],
        [-3.927e-01,  2.609e-01, -1.006e-01],
        [-3.102e-01,  1.716e-01, -2.278e-01],
        [-1.985e-01,  2.244e-01, -6.770e-02],
        [-1.987e-01,  2.510e-01,  2.580e-02],
        [-4.946e-01,  1.810e-02, -1.466e-01],
        [-4.818e-01,  5.840e-02, -2.477e-01],
        [-5.126e-01, -1.188e-01, -1.580e-01],
        [-5.983e-01, -1.416e-01, -1.124e-01],
        [-5.988e-01,  8.430e-02, -6.390e-02],
        [-6.920e-01,  2.110e-02, -6.640e-02],
        [-6.319e-01,  2.062e-01, -1.175e-01],
        [-6.153e-01,  2.015e-01, -2.153e-01],
        [ 3.522e-01, -1.510e-02, -8.510e-02],
        [ 3.180e-01, -8.630e-02, -1.631e-01],
        [ 2.516e-01,  5.760e-02, -3.060e-02],
        [ 4.749e-01, -3.000e-04,  1.023e-01],
        [ 4.132e-01,  9.320e-02,  1.143e-01],
        [ 4.124e-01, -9.000e-02,  1.590e-02],
        [ 4.621e-01,  6.700e-02, -1.479e-01],
        [ 4.570e-01,  6.340e-02, -2.573e-01],
        [ 4.565e-01,  2.009e-01, -1.108e-01],
        [ 3.931e-01,  2.455e-01, -1.740e-01],
        [ 5.993e-01,  1.640e-02, -1.082e-01],
        [ 6.113e-01, -8.950e-02, -1.361e-01],
        [ 6.934e-01,  9.900e-02, -1.702e-01],
        [ 7.101e-01,  5.940e-02, -2.592e-01],
        [ 6.112e-01,  3.010e-02,  4.420e-02],
        [ 6.456e-01,  1.324e-01,  6.930e-02],
        [ 7.036e-01, -6.480e-02,  8.770e-02],
        [ 7.675e-01, -1.880e-02,  1.494e-01],
        [ 4.870e-01, -5.790e-02,  2.425e-01],
        [ 3.875e-01, -8.180e-02,  2.849e-01],
        [ 5.313e-01,  2.040e-02,  3.075e-01],
        [ 5.611e-01, -1.746e-01,  2.441e-01],
        [ 6.201e-01, -1.754e-01,  3.242e-01]
    ] * nm

    raffinose = ["AGLC", "BFRU", "AGAL"]
    raff_patch = {"RAFF": [0, 1, 2]}
    mol = sequence_to_mol(raffinose, RTF, patches=raff_patch)
    mol.to_psf_file('RAFF.psf')

    psf = app.CharmmPsfFile('RAFF.psf')
    topology = psf.topology
    os.remove('RAFF.psf')

    carb_prm = os.path.join(TOPPAR_DIRECTORY, 'par_all36_carb.prm')
    param_set = app.CharmmParameterSet(carb_prm)
    system = psf.createSystem(param_set)

    integrator = mm.LangevinIntegrator(298*K, 1.0/ps, 2.0*fs)
    simulation = app.Simulation(topology, system, integrator)

    simulation.context.setPositions(raff_coords)
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilojoules_per_mole
    assert (2276.6 < energy < 2276.7)
