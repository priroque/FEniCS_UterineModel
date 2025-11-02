import solver as uterus
import numpy as np
from dolfin import *
from fenics import *
import ldrb
import time
import ufl

# -----------------------------------------------------------------------------
# Simulation targeting gestational week 33.57142857 - Visit 4
# -----------------------------------------------------------------------------

filename_Out = "./results_outputs/SIM_uterus_p13.xdmf"

# Number of iteration steps
nsteps = 750

# Gestational week corresponding to the adopted geometry
t0_geo = 16.57142857

# Intrauterine pressure (34 weeks) = 1.747 kPa - Reference: Fisk (1992)
iup = 0.001747     

# Cervix parameter - Reference: Fernandez et al. (2016)        
c_cvx = 0.3715    

# Material parameters - Reference: Sharifimajd et al.
c = 0.06               # 60 kPa = 6 N/cm2 = 0.06 N/mm2
c_1 = 0.006            # 6 kPa = 0.6 N/cm2 = 0.006 N/mm2
c_2 = 1.0              # Dimensionless
zeta1 = 0.4

# Material parameters
prmHolzapfel = [c_cvx, c, c_1, c_2, zeta1]

# Growth parameters
alpha = 0.022
taul, vmaxl, lcritl, gammal = 6.0, 1.6, 1.05, 1.0       
tauc, vmaxc, lcritc, gammac = 6.0, 1.3, 1.05, 1.0  

# Longitudinal growth parameters
prmGrowthL = [taul, vmaxl, lcritl, gammal]

# Circumferential growth parameters
prmGrowthC = [tauc, vmaxc, lcritc, gammac]

# Perimetrium spring constants (directional stiffness)
k_peri = 1e-5          # Anterior region - 0.01 kPa  (1.0e-5 N/mm2)
k_peri_post = 7.5e-5   # Posterior region - 0.075 kPa (7.5e-5 N/mm2)
k_peri_pelv = 7.5e-5   # Inferior (pelvic) region - 0.075 kPa (7.5e-5 N/mm2)

# Spring constants
k_mola = [k_peri, k_peri_post, k_peri_pelv]    

print('\n\n---------------------------------------------------------------')
print('SIMULATION SETUP')
print(f'Material parameters : [c_cvx, c, c_1, c_2, zeta1] = {prmHolzapfel}')
print(f'Spring contants: k_peri = {k_peri}, k_peri_post = {k_peri_post}, k_peri_pelv = {k_peri_pelv}')
print(f'Longitudinal growth parameters: [vmax_l, tau_l, lcrit_l] = {prmGrowthL}')
print(f'Circumferential growth parameters: [vmax_c, tau_c, lcrit_c] = {prmGrowthC}')
print('---------------------------------------------------------------\n\n')

# NOTE: IUP is expressed in N/mm2

# -----------------------------------------------------------------------------
# Mesh loading
# -----------------------------------------------------------------------------
meshname = './mesh_p13/malha-2volumes-ascii'

mesh = Mesh(meshname + '.xml')
bdry = MeshFunction("size_t", mesh, meshname + '_facet_region.xml')
materials = MeshFunction("size_t", mesh, meshname + '_physical_region.xml')

# -----------------------------------------------------------------------------
# Fiber fields setup
# -----------------------------------------------------------------------------
fiber_space = "DG_0"
V = ldrb.space_from_string(fiber_space, mesh, dim=3)

f0 = Function(V)    # longitudinal direction
s0 = Function(V)    # circumferential direction
sn0 = Function(V)   # normal direction

# Load fiber orientations
with HDF5File(mesh.mpi_comm(), "./mesh_p13/fiber_p13v2.h5", "r") as h5file:
    h5file.read(f0, "/fiber")
    h5file.read(s0, "/sheet")
    h5file.read(sn0, "/sheet_normal")

# -----------------------------------------------------------------------------
# FEniCS optimization settings
# -----------------------------------------------------------------------------
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = 4

# -----------------------------------------------------------------------------
# Output files
# -----------------------------------------------------------------------------
file_results = XDMFFile(filename_Out)
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True

# Save mesh and fiber fields for visualization
outname = "./mesh_p13/mesh_fiber_2volumes.xdmf"
fileOutput = XDMFFile(mesh.mpi_comm(), outname)
fileOutput.write(mesh)
fileOutput.parameters['rewrite_function_mesh'] = False
fileOutput.parameters["functions_share_mesh"] = True
fileOutput.parameters["flush_output"] = True

# Rename and export fiber directions
f0.rename("long", "long")
s0.rename("circ", "circ")
sn0.rename("rad", "rad")

fileOutput.write(f0, 0)
fileOutput.write(s0, 0)
fileOutput.write(sn0, 0)

# -----------------------------------------------------------------------------
# Run simulation
# -----------------------------------------------------------------------------
try:
    uterus.UterusSolve(
        fileOutput, file_results, mesh, bdry, materials,
        f0, s0, t0_geo, prmHolzapfel, prmGrowthL, prmGrowthC,
        iup, k_mola, alpha, nsteps)

except Exception as e:
    print('\n---------------------------------------------------------------')
    print(f'Simulation error encountered: {e}')
    print('---------------------------------------------------------------\n')

    # Save error details to log file
    log_entry = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Simulation error: {e}\n"
    with open("simulation_errors.log", "a") as log_file:
        log_file.write(log_entry)
    print('---------------------------------------------------------------\n')