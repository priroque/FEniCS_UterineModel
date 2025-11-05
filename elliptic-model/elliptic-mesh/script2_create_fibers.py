from dolfin import *
import ldrb	

# Codigo para criar a salvar malha e fibras do utero
# Caso: 20 semanas


# Convert from gmsh mesh to fenics
meshname = "malha_20w"
# mesh, ffun, markers = ldrb.gmsh2dolfin(meshname + ".msh")

# Reading mesh converted by dolfin-convert
mesh = Mesh(meshname + '.xml')
materials = MeshFunction("size_t", mesh, meshname + '_physical_region.xml')
ffun = MeshFunction("size_t", mesh, meshname + '_facet_region.xml')

# # Update the markers which are stored within the mesh
# ldrb_markers = {
# 	"base": geo.markers["BASE"][0],
# 	"lv": geo.markers["ENDO"][0],
# 	"epi": geo.markers["PERI"][0],
# 	}

# Update the markers which are stored within the mesh
ldrb_markers = {
    "base": 10,
    "lv": 20,
    "epi": 30
}

# Choose space for the fiber fields
fiber_space = "DG_0"

# Decide on the angles you want to use
angles = dict(
	alpha_endo_lv=90,  # Fiber angle on the endocardium
	alpha_epi_lv=90,  # Fiber angle on the epicardium
	beta_endo_lv=0,  # Sheet angle on the endocardium
	beta_epi_lv=0,
)  # Sheet angle on the epicardium


# Compte the microstructure
fiber, sheet, sheet_normal = ldrb.dolfin_ldrb(
	mesh=mesh, fiber_space=fiber_space, ffun=ffun, markers=ldrb_markers, **angles
	)

fiber.rename("f_0","f_0")
sheet.rename("s_0","s_0")
sheet_normal.rename("n_0","n_0")

# Store the results
with HDF5File(mesh.mpi_comm(), "fiber_20w.h5", "w") as h5file:
	h5file.write(fiber, "/fiber")
	h5file.write(sheet, "/sheet")
	h5file.write(sheet_normal, "/sheet_normal")

# or xdmf
with XDMFFile(mesh.mpi_comm(), "mesh_20w.xdmf") as xdmf:
	xdmf.write(mesh)

# You should also save the facet function
with XDMFFile(mesh.mpi_comm(), "ffun_20w.xdmf") as xdmf:
	xdmf.write(ffun)