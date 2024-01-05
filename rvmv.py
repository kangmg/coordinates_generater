import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from scipy.spatial.transform import Rotation as R

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

def get_adjacency_matrix_from_sdf_with_hydrogens(sdf_file_path):
    suppl = Chem.SDMolSupplier(sdf_file_path, removeHs=False)
    mol = suppl[0]
    AllChem.EmbedMolecule(mol)
    #AllChem.UFFOptimizeMolecule(mol)    # 최적화 여부
    adjacency_matrix = Chem.GetAdjacencyMatrix(mol)
    return adjacency_matrix

def find_orthogonal_point(points, point):
    v1 = np.subtract(points[1], points[0])
    v2 = np.subtract(points[2], points[0])
    
    normal = np.cross(v1, v2)
    d = -np.dot(normal, points[0])
    t = -(np.dot(normal, point) + d) / np.dot(normal, normal)
    
    orthogonal_point = np.add(point, t * normal)
    orthogonal_point = orthogonal_point
    
    return orthogonal_point.tolist()

def find_orthogonal_intersection(points, point):
    p1 = np.array(points[0])
    p2 = np.array(points[1])

    direction_vector = p2 - p1

    point_vector = np.array(point) - p1 

    dot_product = np.dot(direction_vector, point_vector) 

    t = dot_product / np.dot(direction_vector, direction_vector) 
    intersection_point = p1 + t * direction_vector

    return intersection_point.tolist()
    
def frag_rotate(coordi_list, rot_axis_coordi_1, rot_axis_coordi_2, theta):
    axis_vector = np.array(rot_axis_coordi_1) - np.array(rot_axis_coordi_2)

    rot_axis = axis_vector / np.linalg.norm(axis_vector) 
    rotation_matrix = R.from_rotvec(np.radians(theta) * rot_axis).as_matrix()

    rotated_coordi_list = []
    for coordi in coordi_list:
        coordi = np.array(coordi)
        coordi_shifted = coordi - rot_axis_coordi_1
        rotated_coordi_shifted = np.dot(rotation_matrix, coordi_shifted)
        rotated_coordi = rotated_coordi_shifted + rot_axis_coordi_1
        rotated_coordi_list.append(rotated_coordi.tolist())
    return rotated_coordi_list

def center_point(coordi_list):
    l = len(coordi_list)
    mat = np.array(coordi_list)
    x = np.sum(mat[:,0])/l
    y = np.sum(mat[:,1])/l
    z = np.sum(mat[:,2])/l
    center = [x,y,z] #.tolist()         # round
    return center

def frag_axial_move(coordi_list, vec_point_1, vec_point_2, angstrom):
    moved_coordi_list = []
    move_vector = np.array(vec_point_1) - np.array(vec_point_2)
    move_unit_vector = move_vector / np.linalg.norm(move_vector)
    for coordi in coordi_list:
        move_coordi = np.array(coordi) + angstrom * move_unit_vector        #4) # round
        moved_coordi_list.append(move_coordi.tolist())
    return moved_coordi_list

def atom_coordi_extend(atom_list, coordi_list):
    format = ''
    for i, atom in enumerate(atom_list):
        coordi = coordi_list[i]
        line = f"{atom:<2}   {coordi[0]:.4f}   {coordi[1]:.4f}   {coordi[2]:.4f}"
        format = format + line + '\n'
    return format


def get_full_coordi(coordi_list):
    if len(idx_list) != 2:
        print(" !Warning! This program is only tested for 2 fragments")
    else:
        pass
    idx_list_copy = idx_list[:]
    idx_list_copy.pop(rotation_move_idxs)
    residual_idxs = idx_list_copy[0]
    
    residual_coordi_list = []
    for idx in residual_idxs:
        atom_coordi = lines[idx].split()
        coordi = [float(k) for k in atom_coordi[1:]]
        residual_coordi_list.append(coordi)
    
    if rotation_move_idxs == 0:
        coordi_list.extend(residual_coordi_list)
        full_coordi = coordi_list
    elif rotation_move_idxs == 1:
        residual_coordi_list.extend(coordi_list)
        full_coordi = residual_coordi_list
    else:
        print(' !Warning! This program is only tested for 2 fragments')
    return full_coordi


def get_full_index(index_list):
    if len(idx_list) != 2:
        print(" !Warning! This program is only tested for 2 fragments")
    else:
        pass
    idx_list_copy = idx_list[:]
    idx_list_copy.pop(rotation_move_idxs)
    residual_idxs = idx_list_copy[0]
    
    residual_index_list = []
    for idx in residual_idxs:
        atom_coordi = lines[idx].split()
        atom = atom_coordi[0]
        residual_index_list.append(atom)
    
    if rotation_move_idxs == 0:
        index_list.extend(residual_index_list)
        full_index = index_list
    elif rotation_move_idxs == 1:
        residual_index_list.extend(index_list)
        full_index = residual_index_list
    else:
        print(' !Warning! This program is only tested for 2 fragments')
    return full_index


# xyz from sdf file
with open('molecule.sdf', 'r') as file: 
    sdf_data = file.read()
content = sdf2xyz(sdf_data)
lines = content.split('\n')
num_atom = len(lines)


print('-----------------------------------------\nYour molecule is \n-----------------------------------------')
for i in lines:
    print(i)
print('-----------------------------------------\n')


adjacency_matrix = get_adjacency_matrix_from_sdf_with_hydrogens('molecule.sdf')

print('Before fragmentation :')
print(find_disconnected_graphs(adjacency_matrix))

print('Tip) atomic index starts with 0')
n,m = input("\nWhich bond breaks? [two atoms] e.g. 1 5 \n No Bond breaks --> 0 0 \n Answer : ").split()
bond_break = [int(n),int(m)]    

ii = bond_break[0]
jj = bond_break[1]

adjacency_matrix[ii][jj] = 0 
adjacency_matrix[jj][ii] = 0 


print('\nAfter fragmentation :')
print(find_disconnected_graphs(adjacency_matrix))


print("""

     (1)   [*]
            |
           [1]  ( point = 1 )

     (2)    [*]
             |
       [1]---*---[2]  ( line = 1 2 )

     (3)     [*]
              |
        [1]---*---[2]--[3]  ( plane = 1 2 3 )    

""")  

axis_type = int(input(' Make rotational axis with ____ \n 1 : two points \n 2 : point and line(two points) \n 3 : point and plane(three points) \n Answer : '))


# define point
point_input = input("Your point [*] is ___. \n \n 'N' for atom N \n 'N M' for mid-point btw atom N & M \n \nAnswer : ")
idxs = [int(idx) for idx in point_input.split()]
num = len(idxs)
if num == 1:
    atom_coordi = lines[idxs[0]].split()
    point = atom_coordi[1:]
    point = [float(k) for k in point]
elif num == 2:
    atom_coordi_1 = lines[idxs[0]].split()
    atom_coordi_2 = lines[idxs[1]].split()
    point_1 = [float(k) for k in atom_coordi_1[1:]]
    point_2 = [float(j) for j in atom_coordi_2[1:]]
    point = [(point_1[0]+point_2[0])/2 ,(point_1[1]+point_2[1])/2, (point_1[2]+point_2[2])/2]
else:
    print(f'Your point input, {point_input} is invalid. . .')



axis_figure=['''

     (1)   [*]
            |
           [1]  ( point = 1 )

''','''

     (2)    [*]
             |
       [1]---*---[2]  ( line = 1 2 )

''','''

     (3)     [*]
              |
        [1]---*---[2]--[3]  ( plane = 1 2 3 )

''']

print(axis_figure[axis_type-1])


# define plane or line or point

if axis_type == 1:
    axis_atom_idx = int(input(' One atom : '))
    atom_coordi = lines[axis_atom_idx].split()
    axis_coordi = [float(m) for m in atom_coordi[1:]]
elif axis_type == 2:
    line_points = input('Two atoms : ')
    line_points = [int(l) for l in line_points.split()]
    line_points_coordi =[]
    for s in line_points:
        line_split = lines[s].split()
        line_split = [float(k) for k in line_split[1:]]
        line_points_coordi.append(line_split)
    axis_coordi = find_orthogonal_intersection(line_points_coordi, point)
elif axis_type == 3:
    plane_points = input('Three atoms : ')
    plane_points= [int(j) for j in plane_points.split()]
    plane_points_coordi = []
    for i in plane_points:
        line_split = lines[i].split()
        line_split = [float(j) for j in line_split[1:]]
        plane_points_coordi.append(line_split)
    axis_coordi = find_orthogonal_point(plane_points_coordi,point)
else:
    print(f'{axis_type} is wrong input')


idx_list = find_disconnected_graphs(adjacency_matrix)
idx_list = [sorted(frag) for frag in idx_list]
print(' Your fragments : ')
for i, idxs in enumerate(idx_list):
    print(f' {i+1} : {idxs}')

rotation_move_idxs = int(input(' \n What fragment you wanna rotate & move? : ')) - 1
rotation_coordi_list = []
rotation_index_list = []
for idx in idx_list[rotation_move_idxs]:
    atom_coordi = lines[idx].split()
    coordi = [float(k) for k in atom_coordi[1:]]
    atom = atom_coordi[0]
    rotation_coordi_list.append(coordi)
    rotation_index_list.append(atom)


rotation_input = input(' Define the rotational range and N_points \n [-theat degree]  [theta degree]  [N times] e.g. -180 180 36 : ')
rotation_input = rotation_input.split()
delta_theta = (float(rotation_input[1])-float(rotation_input[0]))/int(rotation_input[2])


pop_idx_list = idx_list[:]
pop_idx_list.pop(rotation_move_idxs)
if len(pop_idx_list) == 1:
    pop_idx_list = pop_idx_list[0]
else:
    pass


print('''

 # Define two points for direction vector #

   Point #1   <--------  Point #2  (direction vector)    
                        
''')

print(f'-------------------------------------------------------\n# Point 1 #   in  fixed fragment = {pop_idx_list} \n------------------------------------------------------- \nGive the atom indexes. ')

move_points_1 = input('N           for atom, N \nN M         for mid-point of two atoms, N M \nN M L . .   for center point of given atoms, N M L . . \n\nAtoms? : ')
print('')
print(f'-------------------------------------------------------\n# Point 2 #   in  move & rotate fragment = {idx_list[rotation_move_idxs]} \n------------------------------------------------------- \nGive the atom indexes. ')
move_points_2 = input('N           for atom, N \nN M         for mid-point of two atoms, N M \nN M L . .   for center point of given atoms, N M L . . \n\nAtoms? : ')

idx_list_1 = [int(idx_1) for idx_1 in move_points_1.split()]
idx_list_2 = [int(idx_2) for idx_2 in move_points_2.split()]



center_coordi_1 = []
for idx in idx_list_1:
    center_atm_crd = lines[idx].split()
    center_crd_1 = [float(k) for k in center_atm_crd[1:]]
    center_coordi_1.append(center_crd_1)

center_coordi_2 = []
for idx in idx_list_2:
    center_atm_crd = lines[idx].split()
    center_crd_2 = [float(k) for k in center_atm_crd[1:]]
    center_coordi_2.append(center_crd_2)

center_point_coordi_1 = center_point(center_coordi_1)
center_point_coordi_2 = center_point(center_coordi_2)


move_input = input(' Define the moving range and N_points \n [-ang displacement]  [ang displacement]  [N times] e.g. -3 2 5 : ')
move_input = move_input.split()
delta_ang = (float(move_input[1])-float(move_input[0]))/int(move_input[2])

full_index = get_full_index(rotation_index_list) 
total_coordinates =''
trajectory =''
for i in range(int(move_input[2])+1):
    dis = delta_ang*i
    globals()[f'coordi_list_{dis}_moved'] = frag_axial_move(rotation_coordi_list, center_point_coordi_1, center_point_coordi_2,dis)
    for j in range(int(rotation_input[2])+1):
        rot = delta_theta*j
        moved_coordi_list = globals()[f'coordi_list_{dis}_moved']
        globals()[f'coordi_list_{dis}_moved_{rot}_rotated'] = frag_rotate(moved_coordi_list, point, axis_coordi, rot)
        final_coordi_list = globals()[f'coordi_list_{dis}_moved_{rot}_rotated']
        full_coordi = get_full_coordi(final_coordi_list)
#        xyz_format = f'{len(full_index)}\n\n{atom_coordi_extend(full_index, full_coordi)}'
        xyz_format = atom_coordi_extend(full_index, full_coordi)
        with open(f'coordinates/{np.round(dis,2)}_moved_{np.round(rot,2)}_rotated.xyz','w') as file:
            file.write(xyz_format)
        total_coordinates += xyz_format
        trajectory += str(num_atom) + '\n\n' + xyz_format.strip() + "\n"


with open('total_coordinates.xyz', 'w') as file:
    file.write(total_coordinates)

with open('trajectory.xyz', 'w') as file:
    file.write(trajectory)