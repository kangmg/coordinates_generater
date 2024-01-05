import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops 
from scipy.spatial.transform import Rotation as R



#----------------------- ftn section -------------------------#


def sdf2xyz(sdf_string):
    sdf_lines = sdf_string.split('\n')
    num_atoms = int(sdf_lines[3].split()[0])
    xyz_lines = sdf_lines[4:4+num_atoms]
    full_format= ''
    for line in xyz_lines:
        eles = line.split()
        format = '   '.join([eles[3], eles[0], eles[1], eles[2]])
        full_format += format + '\n'
    return full_format.strip()



def xyz_head_remove(data):
    lines = data.split('\n')
    first_line_split = lines[0].strip().split()
    if len(first_line_split) == 1:
        return 'yes_head'
    elif len(first_line_split) == 4:
        return 'no_head'
    else: 
        print('Check your xyz data files!')


def find_disconnected_graphs(adjacency_matrix):  # adj mat로부터 회전할 frag(atom list)를 찾음
    visited = set()
    disconnected_graphs = []

    for vertex in range(len(adjacency_matrix)):
        if vertex not in visited:
            disconnected_graph = dfs(adjacency_matrix, vertex, visited)
            disconnected_graphs.append(disconnected_graph)

    return disconnected_graphs



def dfs(adjacency_matrix, vertex, visited): # find_disconnected_graphs에 쓰이는 함수
    disconnected_graph = []
    stack = [vertex]

    while stack:
        current_vertex = stack.pop()
        if current_vertex not in visited:
            visited.add(current_vertex)
            disconnected_graph.append(current_vertex)

            for neighbor, connected in enumerate(adjacency_matrix[current_vertex]):
                if connected == 1 and neighbor not in visited:
                    stack.append(neighbor)

    return disconnected_graph


def rotate_point(coordi, rot_axis_coordi_1, rot_axis_coordi_2, theta): 
    axis_vector = np.array(rot_axis_coordi_1) - np.array(rot_axis_coordi_2)

    rot_axis = axis_vector / np.linalg.norm(axis_vector) 
    rotation_matrix = R.from_rotvec(np.radians(theta) * rot_axis).as_matrix()

    coordi = np.array(coordi)
    coordi_shifted = coordi - rot_axis_coordi_1
    rotated_coordi_shifted = np.dot(rotation_matrix, coordi_shifted)
    rotated_coordi = rotated_coordi_shifted + rot_axis_coordi_1
    #rotated_coordi_list.append(rotated_coordi.tolist())
    return np.array(rotated_coordi)





def get_adjacency_matrix_from_sdf_with_hydrogens(sdf_file_path):
    suppl = Chem.SDMolSupplier(sdf_file_path, removeHs=False)
    mol = suppl[0]
    AllChem.EmbedMolecule(mol)
    #AllChem.UFFOptimizeMolecule(mol) # 최적화
    adjacency_matrix = Chem.GetAdjacencyMatrix(mol)
    return adjacency_matrix



adjacency_matrix = get_adjacency_matrix_from_sdf_with_hydrogens('molecule.sdf')

print('')
print('Fragments : ',find_disconnected_graphs(adjacency_matrix))

#---------------------- input section -----------------------#


print('\n Tip) atomic index starts with 0\n')
n,m = input("Rotational axis(bond)? [two atoms] e.g. 1 5: ").split()
bond_break = [int(n),int(m)]    

ii = bond_break[0]
jj = bond_break[1]

adjacency_matrix[ii][jj] = 0
adjacency_matrix[jj][ii] = 0


#print('')
#print('Fragments : ',find_disconnected_graphs(adjacency_matrix))

print('')
N_points, theta = input("N_points, theta(deg)? e.g. 6 60 : ").split()
N_points, theta = [int(N_points),int(theta)]


#--------------- reads the molecular coordinates to rotate ------------------#


# read from xyz format file
'''
with open('molecule.txt', 'r') as file: 
    pre_data = file.read()
flag = xyz_head_remove(pre_data)
if flag == 'yes_head':
    lines = pre_data.strip().splitlines()[2:]
elif flag == 'no_head':
    lines = pre_data.strip().splitlines()
data = '\n'.join(lines)
'''


# xyz from sdf file
with open('molecule.sdf', 'r') as file: 
    sdf_data = file.read()
data = sdf2xyz(sdf_data)

coordinates = [[float(coord) for coord in line.split()[1:]] for line in data.strip().split('\n')]

atoms = [line.split()[0] for line in data.strip().split('\n')]
num_atom = len(atoms)

#--------------- rotate_point() input section ------------------#
        
atom1 = coordinates[ii]
atom2 = coordinates[jj]


disconnected_graphs = find_disconnected_graphs(adjacency_matrix) 


print('Your fragments\n---------------------------')
frags_list = find_disconnected_graphs(adjacency_matrix)
for i, frag_idx in enumerate(frags_list):
    print(i, ' : ', frag_idx)
print('---------------------------')

rot_frag_idx = input("What fragment to rotate? : ")

rot_idx = disconnected_graphs[int(rot_frag_idx)]

remove = [ii, jj]
for element in remove:
    try:
        rot_idx.remove(element)
    except ValueError:
        pass

#--------------- run rotation  ------------------#


delta_theta = theta / N_points 

total_coordinates =''
trajectory =''
for i in range(N_points): 
    rotated_coords_array = np.array(coordinates[:]) 
    theta = int(delta_theta * (i + 1)) 
    for j in rot_idx: 
        rotated_coords_array[j] = rotate_point(rotated_coords_array[j], atom1, atom2, theta) 
    rotated_coords_array = np.round(rotated_coords_array, 4).tolist()
    globals()[f'rotated_coords_{theta}'] = rotated_coords_array 

#--------------- reformat molecular coordinates  ------------------#

    geo_format = ""
    for atom, rtd_crd in zip(atoms, globals()[f'rotated_coords_{theta}']):
        geo_format += f"{atom:<4s}{'   '.join(map(str, rtd_crd))}\n"
        
    total_coordinates += geo_format.strip() + "\n"
    trajectory += str(num_atom) + '\n\n' + geo_format.strip() + "\n\n"

#--------------- save as xyz file  ------------------#

    with open(f'coordinates/{theta}.xyz', "w") as file:
        file.write(geo_format.strip())

with open('trajectory.xyz', 'w') as file:
    file.write(trajectory)

with open('total_coordinates.xyz', 'w') as file:
    file.write(total_coordinates)

with open('coordinates/0.xyz', 'w') as file:
    file.write(data)