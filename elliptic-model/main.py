import solver as uterus
from dolfin import *
from fenics import *
import ldrb
import time

# -----------------------------------------------------------------------------
# Simulation targeting gestational week 40
# -----------------------------------------------------------------------------

filename_Out = "./results_outputs/SIM_elliptic_uterus.xdmf"

# Number of iteration steps
nsteps = 75

# Intrauterine pressure (40 weeks) = 2.3998kPa kPa - Reference: Fisk (1992)
iup = 0.0024 

# Material parameters - Reference: Sharifimajd, Thore, and Stalhand (2016)
c = 0.06               # 60 kPa = 6 N/cm2 = 0.06 N/mm2
c_1 = 0.006            # 6 kPa = 0.6 N/cm2 = 0.006 N/mm2
c_2 = 1.0              # Dimensionless
zeta1 = 0.4
c_3, c_4 = c_1, c_2

# Material parameters
prmHolzapfel = [c, c_1, c_2, c_3, c_4, zeta1]

# Growth parameters
taul, vmaxl, lcritl, gammal = 3.0, 1.6, 1.02, 1.0
tauc, vmaxc, lcritc, gammac = 3.0, 1.3, 1.02, 1.0

# Longitudinal growth parameters
prmGrowthL = [taul, vmaxl, lcritl, gammal]

# Circumferential growth parameters
prmGrowthC = [tauc, vmaxc, lcritc, gammac]

print('\n\n---------------------------------------------------------------')
print('SIMULATION SETUP')
print(f'Material parameters : [c, c_1, c_2, c_3, c_4, zeta1] = {prmHolzapfel}')
print(f'Longitudinal growth parameters: [vmax_l, tau_l, lcrit_l, gammal] = {prmGrowthL}')
print(f'Circumferential growth parameters: [vmax_c, tau_c, lcrit_c, gammac] = {prmGrowthC}')
print('---------------------------------------------------------------\n\n')

# NOTE: IUP is expressed in N/mm2

# -----------------------------------------------------------------------------
# Mesh loading
# -----------------------------------------------------------------------------

# Read mesh
mesh = Mesh()
with XDMFFile(f"./elliptic-mesh/mesh_20w.xdmf") as xdmf:
    xdmf.read(mesh)

# Read markers
# Base = 10, Endometrio = 20, Perimetrio = 30
bdry = MeshFunction("size_t", mesh, 2)
with XDMFFile(f"./elliptic-mesh/ffun_20w.xdmf") as xdmf:
    xdmf.read(bdry)

# -----------------------------------------------------------------------------
# Fiber fields setup
# -----------------------------------------------------------------------------

fiber_space = "DG_0"
V = ldrb.space_from_string(fiber_space, mesh, dim=3)

f0 = Function(V)    # longitudinal direction
s0 = Function(V)    # circumferential direction
sn0 = Function(V)   # normal direction

# Read fibers
with HDF5File(mesh.mpi_comm(), f"./elliptic-mesh/fiber_20w.h5", "r") as h5file:
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
outname = f"./elliptic-mesh/mesh_fiber_20w.xdmf" 
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
        file_results, mesh, bdry, f0, s0, 
        prmHolzapfel, prmGrowthL, prmGrowthC,
        iup, nsteps)

except Exception as e:
    print('\n---------------------------------------------------------------')
    print(f'Simulation error encountered: {e}')
    print('---------------------------------------------------------------\n')

    # Save error details to log file
    log_entry = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Simulation error: {e}\n"
    with open("simulation_errors.log", "a") as log_file:
        log_file.write(log_entry)
    print('---------------------------------------------------------------\n')

