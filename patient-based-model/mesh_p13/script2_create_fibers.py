# Generating fibers for patient-specific geometries

# CONVERT MESH FROM .msh TO .xml VIA TERMINAL:
# dolfin-convert MESH-NAME.msh MESH-NAME.xml

import ldrb
from dolfin import *


# Convert from gmsh mesh to fenics
meshname = 'malha-2volumes-ascii'
mesh, ffun, markers = ldrb.gmsh2dolfin(meshname + ".msh")

# Reading mesh converted by dolfin-convert
mesh = Mesh(meshname + '.xml')
materials = MeshFunction("size_t", mesh, meshname + '_physical_region.xml')
ffun = MeshFunction("size_t", mesh, meshname + '_facet_region.xml')

# Marking the uterine base as endo
ffun.array()[ffun.array() == 10] = 20

# Marking the cervix base as base
ffun.array()[ffun.array() == 50] = 10

# Marking the cervix endo as endo
ffun.array()[ffun.array() == 40] = 20

# Marking the cervix epi as epi
ffun.array()[ffun.array() == 60] = 30

# Marking the cervix fundus (equivalent to uterine base) as epi
ffun.array()[ffun.array() == 70] = 30

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