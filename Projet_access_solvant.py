<<<<<<< HEAD
#Accesibilité au solvant oui
=======
#Accesibilité au solvant
>>>>>>> 45721821a79069941e2e5c28f9ce635a87004ffb

#Module
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import math

#function
def read_PDB(filename):
    liste = []
    with open(filename, "r") as file:
        for ligne in file:
            if (ligne[0:4] == "ATOM") & (ligne[77:78] != "H"):
                dico = {"Atome" : ligne.split()[-1],
                "AA" : ligne.split()[3],
                "x" : float(ligne.split()[6]),
                "y" : float(ligne.split()[7]),
                "z" : float(ligne.split()[8])}
                liste.append(dico)
    return liste

#Programme principal
#file_PDB = "M2BI/1ben.pdb"
file_PDB = "//wsl.localhost/Ubuntu/home/qcapdet/M2BI/Project/7kh5.pdb"
dico = read_PDB(file_PDB)
#print(dico)

# Rayon de van der Waals pour le carbone, l'azote, l'oxygène et le soufre en Ångströms
rayon_vdW_carbone = 1.70
rayon_vdW_azote = 1.55
rayon_vdW_oxygène = 1.4 #1.52
rayon_vdW_soufre = 1.80

# Voisin des atomes à une distance limite de 5 Ångströms

for dic in dico:
    if dic["Atome"] == "C":
        dic["vdW"] = 1.70
    if dic["Atome"] == "N":
        dic["vdW"] = 1.55
    if dic["Atome"] == "O":
        dic["vdW"] = 1.52
    if dic["Atome"] == "S":
        dic["vdW"] = 1.80
print(dico[0]["vdW"])

# Distance seuil
distance_seuil = 5.0  # 5 Ångströms

# Liste pour les atomes voisins
Nbr_voisins = []
Groupes = []
# Parcourir tous les atomes
for i in range(len(dico)):
    voisins = 0

    # Parcourir tous les autres atomes
    for j in range(len(dico)):
        if i != j:

            # Calcul de la distance euclidienne entre les atomes i et j
            distance = math.sqrt(
                (dico[i]["x"] - dico[j]["x"]) ** 2 +
                (dico[i]["y"] - dico[j]["y"]) ** 2 +
                (dico[i]["z"] - dico[j]["z"]) ** 2
            )
            #print(distance)
            if distance <= distance_seuil:
                voisins += 1
                Groupes.append([i+1, dico[i]["x"], dico[i]["y"], dico[i]["z"] , dico[i]["Atome"], dico[i]["AA"], dico[j]["Atome"], dico[j]["vdW"], dico[j]["x"], dico[j]["y"], dico[j]["z"]])

    Nbr_voisins.append(voisins)
    
# Afficher le nombre de voisins pour chaque atome
for i, nb in enumerate(Nbr_voisins):
    print(f"Atome {i + 1}: {nb} voisins")
#print(Nbr_voisins)
#print(Groupes)

# Sélectionner un groupes d'atomes voisins
""" for i in Groupes:
    if 1 in i:
        print(i) """

# Création d'une sphère de 20 points avec un rayon adapté selon la nature des atomes

# Calcul du rayon de la sphère
rayon_H20 = 1.4
rayon_sphere_C = rayon_vdW_carbone + rayon_H20
rayon_sphere_N = rayon_vdW_azote + rayon_H20
rayon_sphere_O = rayon_vdW_oxygène + rayon_H20
rayon_sphere_S = rayon_vdW_soufre + rayon_H20

# Nombre de points sur la sphère
nombre_de_points = 92

# Générer des points uniformément répartis sur la sphère
phi = np.random.uniform(0, 2 * np.pi, nombre_de_points)
cos_theta = np.random.uniform(-1, 1, nombre_de_points)
theta = np.arccos(cos_theta)

# Calculer les coordonnées cartésiennes x, y, z à partir des coordonnées sphériques
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

    # Associer un indice à chaque atome
    indices = np.full_like(x, i + 1).astype(int)

    # Arrondir les valeurs de x, y, et z à 0.001 près
    x = np.round(x, 3)
    y = np.round(y, 3)
    z = np.round(z, 3)

    # Créer une liste de coordonnées pour l'atome actuel
    coord_atome = []

    for j in range(nombre_de_points):
        coord_atome.append([indices[j], x[j], y[j], z[j]])

    # Ajouter la liste de coordonnées de l'atome à la liste principale
    coord_points.extend(coord_atome)

# Afficher les coordonnées avec l'indice
###for coord in coord_points:
###    print(coord)

    # Vous pouvez également afficher les coordonnées des points si nécessaire :
    ###print(coord_points)


"""      # Créer une figure 3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Afficher les points
        ax.scatter(x, y, z, c='b', marker='o', s=10)

        # Paramètres d'axe
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Titre du graphique
        ax.set_title('Sphère composée de 92 points')

        # Afficher la figure
        plt.show() """
    

""" # Créez une liste de rayons pour chaque type d'atome
rayons_sphere = {
    "C": rayon_vdW_carbone,
    "N": rayon_vdW_azote,
    "O": rayon_vdW_oxygène,
    "S": rayon_vdW_soufre
} 
 """

# Recherche des points de contacts 

# Distance seuil
distance_seuil = rayon_vdW_oxygène  # 1.4 1.52 Ångströms

# Tableau NumPy pour les coordonnées des atomes
atom_coordinates = np.array([[dic["x"], dic["y"], dic["z"]] for dic in dico])

# Tableau NumPy pour les coordonnées des points de la sphère
sphere_coordinates = np.array([[pt[1], pt[2], pt[3]] for pt in coord_points])

# Tableau NumPy pour les indices des atomes dans la sphère
atom_indices_in_sphere = np.array([pt[0] for pt in coord_points])

# Tableau NumPy pour les indices des atomes voisins dans Groupes
group_indices = np.array([k[0] for k in Groupes])

# Tableau NumPy pour les coordonnées des voisins de chaque atome avec leur rayon de vdW respectif
Groupes_np = np.array([[k[7], k[8], k[9], k[10]] for k in Groupes])

# Tableau NumPy pour stocker les résultats des points de contact pour chaque atome
pts_access = np.zeros(len(dico), dtype=int)

# Parcourir les atomes
for i in range(len(dico)):
    # Coordonnées de l'atome en cours
    atom_coord = atom_coordinates[i]

    # Coordonnées de la sphère associée à l'atome
    atom_sphere = sphere_coordinates[atom_indices_in_sphere == i + 1]

    # Coordonnées des voisins de l'atome en cours
    atom_neighbors = Groupes_np[group_indices == i + 1][:, 1:]

    # Récupérer les rayons de vdW des voisins de l'atome en cours
    neighbor_vdw_radii = Groupes_np[group_indices == i + 1][:, 0]

    # Parcourir les points de la sphère
    for sphere_point in atom_sphere:
        #print(f"S{sphere_point}")
        # Calculer la distance euclidienne entre le point de la sphère et les voisins de l'atome
        distances = np.linalg.norm(atom_neighbors - sphere_point, axis=1)
        #print(distances)
        # Vérifier si la différence entre la distance et le rayon de vdW du voisin est inférieure à la distance seuil
        contact_mask = (distances - neighbor_vdw_radii) < distance_seuil
        #print(f"res{(distances - neighbor_vdw_radii)}")
        #print((distances - neighbor_vdw_radii) < distance_seuil)
        #print(distances, neighbor_vdw_radii)
        # Si au moins un voisin est en contact, ajoutez un point de contact
        if np.any(contact_mask):
            pts_access[i] += 1
            

# Afficher le résultat
indices = np.arange(1, len(dico) + 1)
print(f"Atomes {indices}: {pts_access}/{nombre_de_points}")
print("fini !")
print(np.sum(pts_access))

""" for i, nb_contacts in enumerate(pts_access):
    print(f"Atome {i + 1} : {nb_contacts}/{nombre_de_points}") """



""" # Calculez les distances entre les atomes et les points de la sphère en utilisant des opérations NumPy
atom_coordinates_expanded = atom_coordinates[:, np.newaxis, :]
distances = np.linalg.norm(atom_coordinates_expanded - sphere_coordinates, axis=2)

# Créez un masque booléen pour les points de contact en utilisant les indices des atomes
contact_mask = (distances - Groupes_np[group_indices - 1, 7][:, np.newaxis]) < distance_seuil

# Comptez le nombre de points de contact pour chaque atome
pts_access = np.sum(contact_mask, axis=1)

# Affichez le résultat
indices = np.arange(1, len(dico) + 1)
print(f"Atomes {indices}: {pts_access}/{nombre_de_points}")
print("fini !")
print(np.sum(pts_access)) """








""" # Groupes d'atomes voisins pour chaque atome de la protéine associé au sphère généré
pts_access = []
# Parcourir tous les atomes

for i in range(1,len(dico)+1):
    # print(i)
    # Parcourir la sphère associé à l'atome
    sphere = []
    for j in coord_points:
        if i in j:
            sphere.append(j)

    # Parcourir les points de la sphère
    count = 0
    for pt in sphere:
       # print(pt)

        # Parcourir tous les voisins de l'atome
        vois = []
        for k in Groupes:
            if i in k:
                vois.append(k)

        # Parcourir chaque voisins de l'atome
        for voi in vois:
            # print(voi)
            # Calcul de la distance euclidienne entre un points et un voisins

            distance = math.sqrt(
                (pt[1] - voi[8]) ** 2 +
                (pt[2] - voi[9]) ** 2 +
                (pt[3] - voi[10]) ** 2
            )

            if (distance - voi[7]) < distance_seuil:
                count += 1
        break
    pts_access.append(count)
    print(f"{count}/{nombre_de_points}")
print("fini !")  
print(sum(pts_access)) """

#Nbr_voisins.append(voisins)
    
#print(coord_points)

# Sélectionner un groupes d'atomes voisins
""" for i in Groupes:
    if 1 in i:
        print(i[8], i[9], i[10])
Liste = []
for i in coord_points:
    if 1 in i:
        Liste.append(i)
print(Liste, len(Liste))

# Itérer sur la liste Liste
for j in Liste:
<<<<<<< HEAD
    print(j[1], j[2], j[3]) """
=======
    print(j[1], j[2], j[3]) """
>>>>>>> 45721821a79069941e2e5c28f9ce635a87004ffb
