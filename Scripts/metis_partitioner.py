import os
import numpy as np
import pymetis
import time
import re

start_time = time.time()

# 1. Define Paths
input_dir = "/home/aadhya/STT_RAM_Solver/Gmsh_Files/"
output_dir = "/home/aadhya/STT_RAM_Solver_METIS/Gmsh_METIS/"
os.makedirs(output_dir, exist_ok=True)

# --- HYPER-FAST I/O FUNCTIONS ---
def fast_read(filename, dtype):
    with open(filename, 'r') as f:
        text = re.sub(r'[\n\s]+', ',', f.read())
        return np.array([dtype(x) for x in text.split(',') if x.strip()], dtype=dtype)

def fast_write(filename, array, is_float=False):
    with open(filename, 'w') as f:
        if is_float:
            # Formats exactly like the raw p_0p004 file: 1.000000,0.000000...
            f.write(','.join([f"{x:.6f}" for x in array]))
        else:
            # Formats exactly like the raw t_0p004 file: 3143,3145,3144...
            f.write(','.join(map(str, array)))

print("1. Loading Gmsh files...")
nodes = fast_read(os.path.join(input_dir, "p_0p004.txt"), float).reshape(-1, 3)
triangles = fast_read(os.path.join(input_dir, "t_0p004.txt"), int).reshape(-1, 3) - 1

num_nodes = len(nodes)
num_triangles = len(triangles)

print("2. Building adjacency graph for METIS...")
node_to_tri = [[] for _ in range(num_nodes)]
for i, tri in enumerate(triangles):
    for node in tri:
        node_to_tri[node].append(i)

adjacency_list = [set() for _ in range(num_triangles)]
for i, tri in enumerate(triangles):
    for node in tri:
        for neighbor_tri in node_to_tri[node]:
            if neighbor_tri != i:
                adjacency_list[i].add(neighbor_tri)

adjacency_list = [list(neighbors) for neighbors in adjacency_list]

print("3. Partitioning into 12 Islands...")
n_cuts, membership = pymetis.part_graph(12, adjacency=adjacency_list)
membership = np.array(membership)

print("4. Renumbering Nodes...")
new_node_id = 0
old_to_new_map = np.full(num_nodes, -1, dtype=int)
sorted_triangles = []

for core_id in range(12):
    core_tris = triangles[membership == core_id]
    for tri in core_tris:
        new_tri = []
        for old_node in tri:
            if old_to_new_map[old_node] == -1:
                old_to_new_map[old_node] = new_node_id
                new_node_id += 1
            new_tri.append(old_to_new_map[old_node])
        sorted_triangles.append(new_tri)

print("5. Reordering the Coordinate Matrix...")
optimized_nodes = np.zeros((new_node_id, 3), dtype=float)
for old_node, new_node in enumerate(old_to_new_map):
    if new_node != -1:
        optimized_nodes[new_node] = nodes[old_node]

print("6. Saving to Single-Line CSV...")
final_tris = (np.array(sorted_triangles) + 1).flatten()

# THE FIX: Flattens the array and writes it without newlines
fast_write(os.path.join(output_dir, "p_0p004_opt.txt"), optimized_nodes.flatten(), is_float=True)
fast_write(os.path.join(output_dir, "t_0p004_opt.txt"), final_tris, is_float=False)

print(f"\nDONE! Finished in {time.time() - start_time:.2f} seconds.")