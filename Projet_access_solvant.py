#SUBJECT : CALCULATION OF THE SOLVENT ACCESSIBILE SURFACE OF A PROTEIN 

#Module
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import math


########1# Read a pdb_file #########

#function

def read_PDB(filename):
    # Créez une liste vide pour stocker les dictionnaires
    liste = []

    # Open the PDB file
    with open(filename, "r") as file:
        # Browse each line of the file
        for ligne in file:
            # Check if the line in an "ATOM" line and if the atom is not a hydrogen atom
            if (ligne[0:4] == "ATOM") & (ligne[77:78] != "H"):
                # Extract the amino acid number (n_AA) from the line (column 23 to 29) and remove the space around it
                n_AA = ligne[22:29].strip()
                # Create a dictionnary containing the atom's information
                dico = {
                "Atome" : ligne.split()[-1],
                "AA" : ligne.split()[3],
                "n°AA" : n_AA,
                "x" : float(ligne.split()[6]),
                "y" : float(ligne.split()[7]),
                "z" : float(ligne.split()[8])}
                # Add dictionnary to list
                liste.append(dico)
    return liste

#Main program
file_PDB = "//wsl.localhost/Ubuntu/home/qcapdet/M2BI/Project/7kh5.pdb"
dico = read_PDB(file_PDB)

# Van der Waals (=vdW) for the carbon, nitrogen, oxygen and sulfur in Ångströms
radius_vdW_carbone = 1.70
radius_vdW_azote = 1.55
radius_vdW_oxygène = 1.4
radius_vdW_soufre = 1.80


########2# Neighbor of atome at a limiting distance of 5 Ångströms #########

for dic in dico:
    if dic["Atome"] == "C":
        dic["vdW"] = 1.70
    if dic["Atome"] == "N":
        dic["vdW"] = 1.55
    if dic["Atome"] == "O":
        dic["vdW"] = 1.52
    if dic["Atome"] == "S":
        dic["vdW"] = 1.80

# Treshold distance
treshold_distance = 5.0  # 5 Ångströms

# List for neighboring atoms
Nbr_neighbours = []
Groups = []

# Browse all atoms
for i in range(len(dico)):
    neighbours = 0

    # Browse all other atoms
    for j in range(len(dico)):
        if i != j:

            # Calculation of the Euclidean distance between atoms i and j
            distance = math.sqrt(
                (dico[i]["x"] - dico[j]["x"]) ** 2 +
                (dico[i]["y"] - dico[j]["y"]) ** 2 +
                (dico[i]["z"] - dico[j]["z"]) ** 2
            )

            if distance <= treshold_distance:
                neighbours += 1
                Groups.append([i+1, dico[i]["x"], dico[i]["y"], dico[i]["z"] , dico[i]["Atome"], dico[i]["AA"], dico[j]["Atome"], dico[j]["vdW"], dico[j]["x"], dico[j]["y"], dico[j]["z"]])

    Nbr_neighbours.append(neighbours)
    
""" # Show number of neighbours for each atom
for i, nb in enumerate(Nbr_neighbours):
    print(f"Atom {i + 1}: {nb} neighbours") """


########3# Creation of a sphere of 92 points with a radius adapted according to the nature of the atoms #########

# Calculation of the radius of the sphere
radius_H20 = 1.4
rayon_sphere_C = radius_vdW_carbone + radius_H20
rayon_sphere_N = radius_vdW_azote + radius_H20
rayon_sphere_O = radius_vdW_oxygène + radius_H20
rayon_sphere_S = radius_vdW_soufre + radius_H20

# Number of points on the sphere
number_of_points = 92

# Generate points uniformly ditributed on the sphere
phi = np.random.uniform(0, 2 * np.pi, number_of_points)
cos_theta = np.random.uniform(-1, 1, number_of_points)
theta = np.arccos(cos_theta)

# Calculation of Cartesian coordinates x, y, z from spherical coordinates (phi, theta) and atom type
coord_points = []

for i, dic in enumerate(dico):
    if dic["Atome"] == "C":
        x = dic["x"] + rayon_sphere_C * np.sin(phi) * np.cos(theta)
        y = dic["y"] + rayon_sphere_C * np.sin(phi) * np.sin(theta)
        z = dic["z"] + rayon_sphere_C * np.cos(phi)
    if dic["Atome"] == "N":
        x = dic["x"] + rayon_sphere_N * np.sin(phi) * np.cos(theta)
        y = dic["y"] + rayon_sphere_N * np.sin(phi) * np.sin(theta)
        z = dic["z"] + rayon_sphere_N * np.cos(phi)
    if dic["Atome"] == "O":
        x = dic["x"] + rayon_sphere_O * np.sin(phi) * np.cos(theta)
        y = dic["y"] + rayon_sphere_O * np.sin(phi) * np.sin(theta)
        z = dic["z"] + rayon_sphere_O * np.cos(phi)
    if dic["Atome"] == "S":
        x = dic["x"] + rayon_sphere_S * np.sin(phi) * np.cos(theta)
        y = dic["y"] + rayon_sphere_S * np.sin(phi) * np.sin(theta)
        z = dic["z"] + rayon_sphere_S * np.cos(phi)

    # Assign an indices to each atom
    indices = np.full_like(x, i + 1).astype(int)

    # Round the values of x, y, z to the nearest 0.001
    x = np.round(x, 3)
    y = np.round(y, 3)
    z = np.round(z, 3)

    # List of coordinates for the current atom
    coord_atom = []

    for j in range(number_of_points):
        coord_atom.append([indices[j], x[j], y[j], z[j]])

    # Add the list of atom coordinates to the main list
    coord_points.extend(coord_atom)


"""      ## Creating a figure 3D ##
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Show points
        ax.scatter(x, y, z, c='b', marker='o', s=10)

        # Axis settings
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Chart title
        ax.set_title('Sphere composed of 92 points')

        # Show figure
        plt.show() """


########4# Search for contact points #########

# Treshold distance
treshold_distance = radius_vdW_oxygène  # 1.4 Ångströms

# NumPy array for atom coordinates
atom_coordinates = np.array([[dic["x"], dic["y"], dic["z"]] for dic in dico])

# NumPy array for the coordinates of the points of the sphere 
sphere_coordinates = np.array([[pt[1], pt[2], pt[3]] for pt in coord_points])

# NumPy array for indices of atoms of each sphere
atom_indices_in_sphere = np.array([pt[0] for pt in coord_points])

# NumPy array for indices of neighboring of each atoms in Groups
group_indices = np.array([k[0] for k in Groups])

# NumPy array for the coordinates of the neighbors (atoms) of each atom with their respective vdW radius
Groups_np = np.array([[k[7], k[8], k[9], k[10]] for k in Groups]) # k[7] : vdW radius , k[8], k[9], k[10] : x, y, z

# NumPy array to store the results of the contact points for each atom (total accessible surface in terms of points)
pts_contact = np.zeros(len(dico), dtype=int)

# Browse all atoms
for i in range(len(dico)):

    # Coordinates of the current atom
    atom_coord = atom_coordinates[i]

    # Coordinates of the sphere associated with the atom
    atom_sphere = sphere_coordinates[atom_indices_in_sphere == i + 1]

    # Coordinates of the neighbors of the current atom
    atom_neighbors = Groups_np[group_indices == i + 1][:, 1:]

    # Get the vdW radius of the neighbors of the current atom
    neighbor_vdw_radii = Groups_np[group_indices == i + 1][:, 0]

    # Browse the points of the sphere
    for sphere_point in atom_sphere:

        # Calculation of the Euclidean distance between the sphere point and the neighbors
        distances = np.linalg.norm(atom_neighbors - sphere_point, axis=1)

        # Checks if the difference between the distance and the neighbor's vdW radius is less than the treshold distance
        contact_mask = (distances - neighbor_vdw_radii) < treshold_distance

        # If at least one neighbor is in contact, add a contact point
        if np.any(contact_mask):
            pts_contact[i] += 1
            

########5# Results #########

# Total accessible surface in terms of points for each atom

""" for i, nb_contacts in enumerate(pts_contact):
    print(f"Atome {i + 1} : {number_of_points - nb_contacts} / {number_of_points}") """

# Total accessible surface in terms of points
sum_pts_contact = np.sum(pts_contact)
print(f"Total accessible surface : {(number_of_points * len(dico)) - sum_pts_contact} / {number_of_points * len(dico)}")

# conversion into Å2 of the accessible surface for each residue of the protein and calculation of the percentage of Solvent Accessible Surface Area (SASA)

# dictionnaries
sum_by_AA = {}
surf_A = {}
surf_T = {}

# Browse each contact point for each atom
for i, nb_contacts in enumerate(pts_contact):

    ## Calculation of the area in Å2 of the accessible surface for each atom ##
    aera_access = (1 - (nb_contacts/number_of_points)) * (4 * math.pi * (dico[i]["vdW"] + radius_H20)**2)

    # Number of each amino acid connected to the correponding atom
    num_AA = dico[i]["n°AA"]

    # If the atom has the number of the target residue, then we add the area in Å2 (same for the other correponding atoms)
    if num_AA in sum_by_AA:
            sum_by_AA[num_AA] += aera_access
    else:
            sum_by_AA[num_AA] = aera_access

    ## Calculation of the percentage of SASA ##

    # Accessible surface of the atom in the molecule in Å2
    surface_A = aera_access

    # Total surface of the atom in Å2
    surface_T = 4 * math.pi * (dico[i]["vdW"])**2

    if num_AA in surf_A: # surf_A : sum of surface_A value for each atom
            surf_A[num_AA] += surface_A
    else:
            surf_A[num_AA] = surface_A
    
    if num_AA in surf_T: # surf_T : Idem with surface_T
            surf_T[num_AA] += surface_T
    else:
            surf_T[num_AA] = surface_T

# Creating a dictionnary to map amino acid to names
num_to_nom_AA = {}
with open(file_PDB, "r") as fichier_pdb:
    for ligne in fichier_pdb:
        if ligne.startswith("ATOM") and ligne[77:78] != "H":
            n_AA = int(ligne[22:29])
            nom_AA = ligne[16:20].strip()
            # Add match to dictionnary
            num_to_nom_AA[n_AA] = nom_AA

# Displays the result of the aera in Å2 of the accessible surface for each residue

# Get the sum of values in Å2 for each residue (somme)
for residue, somme in sum_by_AA.items():
    residue = int(residue)

    # Displays the nature of the residue
    nom_AA = num_to_nom_AA[residue]
    print(f"Residue {residue:3} {nom_AA:<4}: Area = {somme:>6.2f} Å²")

# Displays the result of the Relative Solvent Accessibility (RSA)

# MaxASA_value (theoretical values for each residue)
MaxASA_value = {"ALA" : 129, "ARG" : 274, "ASN" : 195, "ASP" : 193,  "CYS" : 167,
                "GLU" : 223, "GLN" : 225, "GLY" : 104, "HIS" : 224, "ILE" : 197,
                "LEU" : 201, "LYS" : 236, "MET" : 224, "PHE" : 240, "PRO" : 159,
                "SER" : 155, "THR" : 172, "TRP" : 285, "TYR" : 263, "VAL" : 174}

for residue, somme in sum_by_AA.items():
    residue = int(residue)

    # Get the theoretical value of MaxASA for this residue
    maxasa_theoretical = MaxASA_value.get(num_to_nom_AA[residue], 0)
    if maxasa_theoretical != 0:
        # Calculate RSA
        RSA = somme / maxasa_theoretical
    else:
         RSA = 0.00

    print(f"Residue {residue:3} {num_to_nom_AA[residue]:<4}: RSA = {RSA:>4.2f}")

# Displays Solvent Accessible Surface Area (SASA) results for each residue

# Get the sum of the values (val) of surface_A and surface_T for each residue 
for val in range(len(surf_T)):
    sval = str(val+1)

    # Calculate the SASA of each residue as a percentage
    SASA = round((surf_A[sval]/surf_T[sval])*100, 2)
    nom_AA = num_to_nom_AA[val+1]
    print(f"Residue {val+1:3} {nom_AA:<4}: SASA = {SASA:>5.2f} %")