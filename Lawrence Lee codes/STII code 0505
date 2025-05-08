import numpy as np
import xlwings as xw


# Function to read input file
def read_input(filename):
    data = {}
    current_key = None
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                current_key = line[2:] # Remove the '#' to get the section title
                data[current_key] = []
            elif line:
                data[current_key].append(list(map(float, line.split())))

    # Convert lists to numpy arrays
    for key in data:
        data[key] = np.array(data[key])

    return data


def number_dofs(nnodes, ndof, restraint):
    """
    Number the degrees of freedom for each node and identify free DOFs.
    """
    identity = np.zeros((nnodes, ndof), dtype=int)
    icount = 0

    # Assign free DOFs
    for inode in range(nnodes):
        for idof in range(ndof):
            if restraint[inode, idof] == 0:
                identity[inode, idof] = icount
                icount += 1
            else:
                identity[inode, idof] = -1
    nfree = icount

    # Assign fixed DOFs
    for inode in range(nnodes):
        for idof in range(ndof):
            if restraint[inode, idof] == 1:
                identity[inode, idof] = icount
                icount += 1

    return identity, nfree



def length(ielem, connectivity, coords):
    """
    Compute the length of the element.
    """
    inode, jnode = connectivity[ielem]
    x1, y1 = coords[inode]
    x2, y2 = coords[jnode]
    L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return L



def stiffness_bar(A, E, L):
    """
    Compute the stiffness matrix for a bar element.
    """
    k = (A * E / L) * np.array([[1, 0, -1, 0],
                                 [0, 0, 0, 0],
                                 [-1, 0, 1, 0],
                                 [0, 0, 0, 0]])
    return k



def transformation_matrix(ielem, connectivity, coords):
    '''
    Compute the transformation matrix for a bar element.
    '''
    inode, jnode = connectivity[ielem]
    x1, y1 = coords[inode]
    x2, y2 = coords[jnode]
    theta = np.arctan2(y2 - y1, x2 - x1)

    return np.array([[np.cos(theta), np.sin(theta), 0, 0],
                  [-np.sin(theta), np.cos(theta), 0, 0],
                  [0, 0, np.cos(theta), np.sin(theta)],
                  [0, 0, -np.sin(theta), np.cos(theta)]])



def xwrite_excel():
    """
    Write results to an Excel file.
    """

    # Create a new workbook and write data
    xlswb = xw.Book()
    xlswb.sheets[0].name = 'Results'
    xlswb.filename = 'results.xlsx'
    sheet = xlswb.sheets[0]
    sheet.range('A1').value = 'Node'
    sheet.range('B1').value = 'Displacement X'
    sheet.range('C1').value = 'Displacement Y'
    sheet.range('D1').value = 'Force X'
    sheet.range('E1').value = 'Force Y'
    sheet.range('G1').value = 'Elements'
    sheet.range('H1').value = 'inode'
    sheet.range('I1').value = 'jnode'
    sheet.range('J1').value = 'Axial Force'

    for i in range(nnodes):
        sheet.range(f'A{i+2}').value = i + 1
        sheet.range(f'B{i+2}').value = nodal_displacements[i, 0]
        sheet.range(f'C{i+2}').value = nodal_displacements[i, 1]
        sheet.range(f'D{i+2}').value = nodal_forces[i, 0]
        sheet.range(f'E{i+2}').value = nodal_forces[i, 1]

    for i in range(nels):
        sheet.range(f'G{i+2}').value = i+1
        sheet.range(f'H{i+2}').value = connectivity[i, 0] + 1
        sheet.range(f'I{i+2}').value = connectivity[i, 1] + 1
        sheet.range(f'J{i+2}').value = member_forces[i]





# Read input data
input_data = read_input('input_01_changed.txt')
coords = input_data['coords']
connectivity = input_data['connectivity'].astype(int)-1 # Convert to zero-based indexing
props = input_data['props']
restraint = input_data['restraint'].astype(int)
fapplied = input_data['fapplied']





# Main Program
ndof = 2 # Number of degrees of freedom per node
nels = connectivity.shape[0] # Number of elements
nnodes = coords.shape[0] # Number of nodes

# Initialize global stiffness matrix and force vector
K = np.zeros((ndof * nnodes, ndof * nnodes))
Q = np.zeros((ndof * nnodes, 1))

#Set up DOF numbering
identity, nfree = number_dofs(nnodes, ndof, restraint)

# Loop through each element
for ielem in range(nels):
    L = length(ielem, connectivity, coords) # Compute length of element
    A, E = props[ielem] # Extract area and modulus of elasticity

    k = stiffness_bar(A, E, L) # Compute element stiffness matrix
    T = transformation_matrix(ielem, connectivity, coords) # Compute transformation matrix
    kstar = T.T @ k @ T # Transform stiffness matrix to global coordinates

    inode, jnode = connectivity[ielem]
    mapping = np.array([identity[inode, 0], identity[inode, 1], identity[jnode, 0], identity[jnode, 1]])

    # Assemble global stiffness matrix
    for i in range(4):
        for j in range(4):
            K[mapping[i], mapping[j]] += kstar[i, j] 

# Assemble force vector
    for i in range(nnodes):
        for j in range(ndof):
            Q[identity[i, j], 0] = fapplied[i, j]

# Solvve for unknown displacements
Kff = K[:nfree, :nfree]
Qf = Q[:nfree]
qf = np.linalg.solve(Kff, Qf)

# Compute reaction forces
Ksf = K[nfree:, :nfree]
Qs = Ksf @ qf

# Organize output results
nodal_displacements = np.zeros((nnodes, ndof))
nodal_forces = np.zeros((nnodes, ndof))
member_forces = np.zeros((nels, 1))

for i in range(nnodes):
    for j in range(ndof):
        dof_index = identity[i, j]
        if dof_index < nfree:
            nodal_displacements[i, j] = qf[dof_index]
            nodal_forces[i, j] = fapplied[i, j]
        else:
            nodal_forces[i, j] = Qs[dof_index - nfree]


# Compute axial forces in each member
for ielem in range(nels):
    inode, jnode = connectivity[ielem]
    L = length(ielem, connectivity, coords)
    A, E = props[ielem]
    k = stiffness_bar(A, E, L)
    T = transformation_matrix(ielem, connectivity, coords)

    # Compute axial force
    u_elem = np.array([nodal_displacements[inode], nodal_displacements[jnode]]).flatten()
    F_elem = (E * A / L) * np.array([-1, 0, 1, 0]) @ T @ u_elem
    member_forces[ielem] = F_elem





# Display results
print("\nNodal Displacements:")
for i in range(nnodes):
    print(f"Node {i+1}: {nodal_displacements[i]}")

print("\nNodal Forces:")
for i in range(nnodes):
    print(f"Node {i+1}: {nodal_forces[i]}")

print("\nAxial Forces in Members:")
for i in range(nels):
    print(f"Element{i+1}: {member_forces[i]}")


xwrite_excel()
