# Function to read input file
import numpy as np
def read_input(filename):
    data = {}
    current_key = None
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                current_key = line[2:]  # Remove '#' to get the section title
                data[current_key] = []
            elif line:
                data[current_key].append(list(map(float, line.split())))
    
    # Convert lists to numpy arrays
    for key in data:
        data[key] = np.array(data[key])

    return data
def number_dofs(nnodes, ndof, restraint):
    identity = np.zeros((nnodes, ndof), dtype=int)
    icount = 0

    # Assign free DOFs
    for inode in range(nnodes):
        for idof in range(ndof):
            if restraint[inode, idof] == 0:
                identity[inode, idof] = icount
                icount += 1

    nfree = icount

    # Assign fixed DOFs
    for inode in range(nnodes):
        for idof in range(ndof):
            if restraint[inode, idof] == 1:
                identity[inode, idof] = icount
                icount += 1

    return identity, nfree
# Function to compute element length
def length(ielem, connectivity, coords):
    inode, jnode = connectivity[ielem]
    x1, y1 = coords[inode]
    x2, y2 = coords[jnode]
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
# Function to compute element stiffness matrix
def stiffness_bar(A, E, L):
    return (A * E / L) * np.array([[1, 0, -1, 0],
                                   [0, 0,  0, 0],
                                   [-1, 0, 1, 0],
                                   [0, 0,  0, 0]])
# Function to compute transformation matrix
def transform_bar(ielem, connectivity, coords):
    inode, jnode = connectivity[ielem]
    x1, y1 = coords[inode]
    x2, y2 = coords[jnode]
    theta = np.arctan2(y2 - y1, x2 - x1)

    return np.array([[ np.cos(theta),  np.sin(theta), 0, 0],
                     [-np.sin(theta),  np.cos(theta), 0, 0],
                     [0, 0,  np.cos(theta),  np.sin(theta)],
                     [0, 0, -np.sin(theta),  np.cos(theta)]])

# Read input data
input_data = read_input('input2.txt')
coords = input_data['coords']
connectivity = input_data['connectivity'].astype(int) - 1  # Convert to 0-based index
props = input_data['props']
restraint = input_data['restraint'].astype(int)
fapplied = input_data['fapplied']
# Main program
ndof = 2  # Degrees of freedom per node
nels = connectivity.shape[0]  # Number of elements
nnodes = coords.shape[0]      # Number of nodes

# Initialize global stiffness matrix and force vector
K = np.zeros((nnodes * ndof, nnodes * ndof))
Q = np.zeros((nnodes * ndof, 1))
kt_list = []
# Set up DOF numbering
identity, nfree = number_dofs(nnodes, ndof, restraint)

# Loop through each element
for ielem in range(nels):
    L = length(ielem, connectivity, coords)  # Compute element length
    A, E = props[ielem, :]                   # Get element properties
    k = stiffness_bar(A, E, L)               # Compute local stiffness matrix
    T = transform_bar(ielem, connectivity, coords)  # Compute transformation matrix
    ksat = T.T @ k @ T                       # Compute transformed stiffness matrix
    kt=k@T
    kt_list.append(kt)
    inode, jnode = connectivity[ielem]
    mapping = [identity[inode, 0], identity[inode, 1],
               identity[jnode, 0], identity[jnode, 1]]

    # Assemble global stiffness matrix
    for i in range(4):
        for j in range(4):
            K[mapping[i], mapping[j]] += ksat[i, j]
    
    # Assemble force vector
    for i in range(nnodes):
        for j in range(ndof):
            Q[identity[i, j], 0] = fapplied[i, j]
# Solve for unknown displacements
Kff = K[:nfree, :nfree]
Qf = Q[:nfree]
qf = np.linalg.solve(Kff, Qf)
# Compute reaction forces
Ksf = K[nfree:, :nfree]
Qs = Ksf @ qf
print(kt_list)
# Organize output results
nodaldisps = np.zeros((nnodes, ndof))
nodalforces = np.zeros((nnodes, ndof))

for i in range(nnodes):
    for j in range(ndof):
        dof_idx = identity[i, j]
        if dof_idx < nfree:
            nodaldisps[i, j] = qf[dof_idx]
            nodalforces[i, j] = fapplied[i, j]
        else:
            nodalforces[i, j] = Qs[dof_idx - nfree]
#算F(ij)
memberforces=np.zeros((nels,1))
for ielem in range(nels):
    inode, jnode = connectivity[ielem]
    q = np.array([*nodaldisps[inode], *nodaldisps[jnode]])
    memberforce= kt_list[ielem]@q
    memberforces[ielem,0]=memberforce[2]
# Display results
print("\nNodal Displacements:")
for i in range(nnodes):
    print(nodaldisps[i])

print("\nNodal Forces:")
for i in range(nnodes):
    print(nodalforces[i])
print("\nMemberforces:")
for i in range(nels):
    print(memberforces[i])
