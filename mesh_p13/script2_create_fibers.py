# # Generating fibers for patient specific geometries
#
# In this demo we will show how to generate fiber orientations from a patient specific geometry. We will use a mesh of an LV that is constructed using gmsh (https://gmsh.info), see https://github.com/finsberg/ldrb/blob/main/demos/mesh.msh
#
# It is important that the mesh contains physical surfaces of the endocardium (lv and rv if present), the base and the epicardium. You can find an example of how to generate such a geometry using the python API for gmsh here: https://github.com/finsberg/pulse/blob/0d7b5995f62f41df4eec9f5df761fa03da725f69/pulse/geometries.py#L160
#
# First we import the necessary packages. Note that we also import `meshio` which is used for converted from `.msh` (gmsh) to `.xdmf` (FEnICS).

# PROBLEMA COM PACOTE ldrb - SOLUCAO:
# Localizar o pacote via terminal: locate ldrb.py
# Abrir no terminal o pacote usando: code ENDERECO OBTIDO ACIMA
# Comentar linhas 653-658 que verifica as fronteiras marcadas!

# CONVERTER MALHA .msh PARA .xml VIA TERMINAL:
# dolfin-convert NOME-DA-MALHA.msh NOME-DA-MALHA.xml

import ldrb
from dolfin import *

# Convert from gmsh mesh to fenics
meshname = 'malha-2volumes-ascii'
mesh, ffun, markers = ldrb.gmsh2dolfin(meshname + ".msh")

#leitura de malha convertida pelo dolfin-convert
mesh = Mesh(meshname + '.xml')
materials = MeshFunction("size_t", mesh, meshname + '_physical_region.xml')
ffun = MeshFunction("size_t", mesh, meshname + '_facet_region.xml')

#marcando a base do utero como endo
ffun.array()[ffun.array() == 10] = 20
#marcando base do colo como base
ffun.array()[ffun.array() == 50] = 10
#marcando endo do colo como endo
ffun.array()[ffun.array() == 40] = 20
#marcando epi do colo como epi
ffun.array()[ffun.array() == 60] = 30
#marcando fundo do colo (equivalete a base do utero) como epi
ffun.array()[ffun.array() == 70] = 30

print(markers)

# Update the markers which are stored within the mesh

ldrb_markers = {
    "base": 10,
    "lv": 20,
    "epi": 30
}

# Select a space for the fibers (here linear lagrange element)

fiber_space = "DG_0"

# Create a dictionary of fiber angles

angles = dict(
    alpha_endo_lv=-90,  # Fiber angle on the endocardium
    alpha_epi_lv=-90,  # Fiber angle on the epicardium
    beta_endo_lv=0,  # Sheet angle on the endocardium
    beta_epi_lv=0,
)

# Run the LDRB algorithm

fiber, sheet, sheet_normal = ldrb.dolfin_ldrb(
    mesh=mesh, fiber_space=fiber_space, ffun=ffun, markers=ldrb_markers, **angles
)


fiber.rename("f_0","f_0")
sheet.rename("s_0","s_0")
sheet_normal.rename("n_0","n_0")

# Store the results
with HDF5File(mesh.mpi_comm(), "fiber_p13v2.h5", "w") as h5file:
	h5file.write(fiber, "/fiber")
	h5file.write(sheet, "/sheet")
	h5file.write(sheet_normal, "/sheet_normal")

#apex_cells = [5668, 5735, 5755, 7279, 7848, 8589, 9092, 9321, 24442, 25337, 27982, 28055, 29472, 30046, 30623, 31256, 31740, 32315, 33225, 33362, 35310, 35969, 36447, 36849, 39570, 39571, 39572, 39573]
#V = FunctionSpace(
#    mesh,
#    VectorElement(
#        family="DG",
#        cell=mesh.ufl_cell(),
#        degree=0,
#        quad_scheme="default",
#    ),
#)
#dofmap = V.dofmap()
#for cell_id in apex_cells:
#    dofs = dofmap.cell_dofs(cell_id) 
#    for dof in dofs:
#        fiber.vector()[dof] = 0.0
#        sheet.vector()[dof] = 0.0
#        sheet_normal.vector()[dof] = 0.0

# fiber.rename("f_0","f_0")
# sheet.rename("s_0","s_0")
# sheet_normal.rename("n_0","n_0")

# with XDMFFile(mesh.mpi_comm(), meshname + ".xdmf") as xdmf:
#     xdmf.parameters.update(
#     {
#         "functions_share_mesh": True,
#         "rewrite_function_mesh": False
#     })
#     xdmf.write(mesh)
#     xdmf.write(fiber, 0)
#     xdmf.write(sheet, 0)
#     xdmf.write(sheet_normal, 0)

# Save to xdmf
# with dolfin.XDMFFile(mesh.mpi_comm(), "patient_fiber.xdmf") as xdmf:
#     xdmf.write(fiber)

# Use this function to save fiber with angles as scalars

#ldrb.fiber_to_xdmf(fiber, "patient_fiber")

# ![_](_static/figures/patient_fiber.png)
