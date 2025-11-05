# Copiar codigo
import gmsh
import sys

# Inicializa o Gmsh
gmsh.initialize()

# Ler o arquivo .geo
geo_file = "malha_20w.geo"  # substitua pelo nome do seu arquivo .geo
gmsh.open(geo_file)

# Definir a versao do arquivo de malha para 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# Sincronizar o modelo
gmsh.model.geo.synchronize()

# Gerar a malha
gmsh.model.mesh.generate(3)  # Supondo que seja um modelo 3D

# Salvar a malha no arquivo "output.msh"
gmsh.write("malha_20w.msh")

# Finalizar o Gmsh
gmsh.finalize()
