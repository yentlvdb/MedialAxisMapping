class Graph:
    def __init__(self, nodes, edges):
        self.nodes = {}
        self.edges = edges
        node_index = 0
        for node in nodes:
            self.nodes[node_index] = {'index': node_index, 'z': node[0], 'y': node[1], 'x': node[2]}
            node_index += 1
            
    def get_node_coordinates(self, node_index):
        node = self.nodes[node_index]['z'],self.nodes[node_index]['y'],self.nodes[node_index]['x']
        return node

    def get_all_subsets_of_three_nodes(self):
        """
        Generates all possible subsets of three nodes from the graph.

        Returns:
          A list of tuples, where each tuple represents a subset of three node indices.
        """
        node_indices = list(self.nodes.values())
        return list(combinations([node['index'] for node in node_indices], 3))
      
    def count_edges_in_subset(self, subset): # uit de class?
          """
          Counts the number of edges present between the nodes in the given subset.
  
          Args:
              subset: A tuple of three node indices.
  
          Returns:
              The number of edges present within the subset.
          """
          edge_count = 0
          for edge in self.edges:
              if edge[0] in subset and edge[1] in subset:
                  edge_count += 1
          return edge_count
  
      def permutations_subset(self, subset, permutation_type):
          """
          Generate permutations of a subset of nodes based on their internal connectivity (edges),
          with different behavior depending on the specified mode (`type`).
      
          Args:
              subset (list): A list of nodes (typically 3 elements) to generate permutations for.
              type (bool): 
                  - If True (LM image): Return all valid permutations according to connectivity rules (permutation mode).
                  - If False (EM image): Return fixed arrangements (non-permutation mode).
      
          Returns:
              list: A list of permutations or fixed sequences depending on the input parameters.
          """
          to_count = []
          for edge in self.edges:
              if edge[0] in subset and edge[1] in subset:
                  to_count.extend(edge)
      
          if not to_count:
              return list(permutations(subset)) if permutation_type else [tuple(subset)] * 6
      
          if len(to_count) == 2:
              fixed_nr = (set(subset) - set(to_count)).pop()
              if permutation_type:
                  return [[fixed_nr] + list(p) for p in permutations(to_count)]
              else:
                  return [(fixed_nr, to_count[0], to_count[1])] * 2
      
          if len(to_count) == 4:
              count = Counter(to_count)
              fixed_nr = [k for k, v in count.items() if v == 2][0]
              other_nrs = list(set(to_count) - {fixed_nr})
              if permutation_type:
                  return [[fixed_nr] + list(p) for p in permutations(other_nrs)]
              else:
                  return [(fixed_nr, other_nrs[0], other_nrs[1])] * 2

def find_matching_subsets(graph1, graph2, min_edgecount=0):
    """
    Finds pairs of 3-node subsets from two graphs with the same number of edges.

    Args:
      graph1: The first Graph object.
      graph2: The second Graph object.

    Returns:
      A list of tuples, where each tuple contains a pair of matching subsets 
      (one from each graph).
    """
    matching_subsets = []
    
    subsets1 = graph1.get_all_subsets_of_three_nodes()
    subsets2 = graph2.get_all_subsets_of_three_nodes()

    for subset1 in subsets1:
        edge_count1 = graph1.count_edges_in_subset(subset1)
        for subset2 in subsets2:
            edge_count2 = graph2.count_edges_in_subset(subset2)
            if edge_count1 == edge_count2:
                if edge_count1 >= min_edgecount:
                    matching_subsets.append((subset1, subset2))
    return matching_subsets

def add_all_possible_permutations(paired_triangles_based_on_edge_count, graph_fm, graph_em):
    """
    Generate all valid pairings of triangle permutations from two graphs (graph_em and graph_fm),
    based on a list of matched triangle pairs and their internal edge structures.

    For each triangle pair:
        - The triangle from `graph_em` is treated as fixed (non-permutational).
        - The triangle from `graph_fm` is permuted according to its internal edges.

    This function aligns each fixed triangle (from `graph_em`) with all valid permutations 
    of the corresponding triangle in `graph_fm`.

    Args:
        paired_triangles_based_on_edge_count (list of tuples): 
            A list of tuples, each containing a pair of triangle node sets (triangle_em, triangle_lm),
            assumed to be matched based on edge similarity or count.
        graph_fm
        graph_em
        
    Returns:
        list of lists: A list where each element is a [fixed_triangle_em, permuted_triangle_fm] pair,
        for each valid permutation in the FM graph.
    """
    
    all_paired_triangles_with_permutations = []

    for triangle_em, triangle_fm in paired_triangles_based_on_edge_count:
        fixed_em_perms = graph_em.permutations_subset(triangle_em, type=False)
        permuted_fm_perms = graph_fm.permutations_subset(triangle_fm, type=True)

        # Use only the first fixed permutation
        fixed_em = fixed_em_perms[0]

        # Pair each permuted FM triangle with the fixed EM triangle
        all_paired_triangles_with_permutations.extend(
            [[fixed_em, fm_perm] for fm_perm in permuted_fm_perms]
        )

    return all_paired_triangles_with_permutations

def top_th(indices_triangle, junctions, k):
    """
    Computes a 4th point located above or below a triangle, along its normal vector,
    scaled by a factor `k`. Used for building tetrahedra.

    Args:
        indices_triangle (list): Indices of the 3 triangle points in `junctions`.
        junctions (list or np.array): 3D coordinates of junctions (shape Nx3).
        k (float): Scalar offset from the triangle plane along the normal.

    Returns:
        np.ndarray: The 4th point in 3D space above or below the triangle.
    """
    j0, j1, j2 = [junctions[i] for i in indices_triangle]
    
    triangle = np.array([j0, j1, j2])
    center = triangle.mean(axis=0)

    # Compute normal vector
    v1 = j0 - j1
    v2 = j0 - j2
    normal = np.cross(v1, v2)
    norm_length = np.linalg.norm(normal)

    if norm_length == 0:
        raise ValueError("Degenerate triangle: normal vector has zero length.")

    normal_unit = normal / norm_length
    top_point = center + k * normal_unit
    return top_point
  
def filter_by_shear_and_scale(paired_combi_triangles, graph_fm, graph_em, scale_thr=1.5, shear_thr=0.1):
    """
    Filters triangle pairs by computing 3D affine transforms and checking if
    the resulting transformation meets specified shear and scale thresholds.

    Args:
        paired_combi_triangles (list): List of EM/FM triangle pairs to test.
        graph_fm: Graph object with node coordinates for FM space.
        graph_em: Graph object with node coordinates for EM space.
        scale_thr (float): Maximum allowed scale factor (Z).
        shear_thr (float): Maximum allowed shear component (S).

    Returns:
        all_aff_tfm (list): List of valid 4x4 affine transformation matrices (numpy arrays).
        corresponding_pairs (list): List of node index sets used to derive each transform.
    """
    all_paired = add_all_possible_permutations(paired_combi_triangles, graph_fm, graph_em)
    
    valid_transforms = []
    valid_pairs = []
    scale_thr_u = 1 / scale_thr  # Lower bound

    for pair_em, pair_fm in all_paired:
        # Get EM triangle coordinates
        em_coords = [np.array(graph_em.get_node_coordinates(i)) for i in pair_em]
        top_em = top_th([0, 1, 2], em_coords, 50)
        corners_em = np.float32(em_coords + [top_em])

        fm_coords = [np.array(graph_fm.get_node_coordinates(i)) for i in pair_fm]

        for k in [50, -50]:
            top_fm = top_th([0, 1, 2], fm_coords, k)
            corners_fm = np.float32(fm_coords + [top_fm])

            # Estimate 3D affine transformation
            retval, M_affine, inliers = cv2.estimateAffine3D(corners_fm, corners_em)

            # Build full 4x4 transform matrix
            M_tf = np.eye(4)
            M_tf[:3, :] = M_affine
            M_tf_inv = np.linalg.inv(M_tf)

            # Decompose transformation
            T, R, Z, S = decompose44(M_tf_inv)

            # Apply shear and scale filtering
            if all(abs(s) < shear_thr for s in S) and all(scale_thr_u < abs(z) < scale_thr for z in Z):
                valid_transforms.append(M_tf_inv)
                valid_pairs.append([[*pair_em, 50], [*pair_fm, k]])

    return valid_transforms, valid_pairs

def affine_transformation_with_highest_nmi(fm_ori, em_ori, all_aff_tfm, all_pairs):
    """
    Finds the affine transformation from a list that best aligns fm_ori to em_ori,
    using Normalized Mutual Information (NMI) over valid overlapping regions.

    Args:
        fm_ori (ndarray): Source 3D volume (moving modality).
        em_ori (ndarray): Target 3D volume (fixed modality).
        all_aff_tfm (list): List of 4x4 affine transformation matrices.
        all_pairs (list): List of corresponding node pairs (metadata for output).

    Returns:
        optimal_aff_tfm (ndarray): The affine transformation with the highest NMI.
        corresponding_pair (list): The node pair associated with this transformation.
    """
    mask_fm = np.ones_like(fm_ori, dtype=np.uint8)
    all_nmi = []

    for aff_tfm in all_aff_tfm:
        # Apply transformation to image and mask
        transformed_fm = ndimage.affine_transform(fm_ori, aff_tfm, output_shape=em_ori.shape)
        transformed_mask = ndimage.affine_transform(mask_fm, aff_tfm, output_shape=em_ori.shape)

        # Flatten and mask overlapping region
        valid_indices = transformed_mask == 1
        fm_vals = transformed_fm[valid_indices].flatten()
        em_vals = em_ori[valid_indices].flatten()

        if len(fm_vals) == 0 or len(em_vals) == 0:
            all_nmi.append(-1)  # Penalize non-overlapping regions
            continue

        nmi = cluster.normalized_mutual_info_score(fm_vals, em_vals)
        all_nmi.append(nmi)

    # Find best transformation
    best_idx = int(np.argmax(all_nmi))
    return all_aff_tfm[best_idx], all_pairs[best_idx]

def label_and_view_points(points, viewer):
    """
    Adds points to the viewer with labels that indicate their index.

    Args:
        points (ndarray): Array of 3D points to be displayed.
        viewer: Visualization viewer object (likely from napari or similar).
    """
    features = {
        'label': np.arange(len(points)),  # Labels as indices
    }
    
    text = {
        'string': 'label       {label:d}',
        'size': 10,
        'color': 'cyan',
    }
    
    viewer.add_points(
        points,
        features=features,
        text=text,
        size=5,
        face_color='cyan',
    )

def view_points(points, viewer):    
    """
    Adds points to the viewer without labels.

    Args:
        points (ndarray): Array of 3D points to be displayed.
        viewer: Visualization viewer object.
    """
    viewer.add_points(
        points,
        size=5,
        face_color='cyan',
    )

def view_edges(E, V, viewer):
    """
    Visualizes the edges between points using vectors in 3D space.
    
    Args:
        E (list of tuples): List of edges where each edge is defined by two indices.
        V (ndarray): Array of vertices (points), where each vertex is a 3D point.
        viewer: Visualization viewer object.
    """
    edges = [
        np.array([V[el[0]], V[el[1]]])  # Pair of vertices for each edge
        for el in E
    ]
    vectors_layer = viewer.add_vectors(
        edges,
        edge_width=1.0,
        edge_color='cyan',  # Color of edges
    )
  
def edges_from_edgepointlayer(pointlayer, V, projected=False, threshold=5.0):
    """
    Extracts edges from a point layer and matches them to vertices in V based on distance.

    Args:
        pointlayer (napari.viewer.layers.points.Points): The point layer containing the points defining edges.
        V (ndarray): Array of vertices in 3D space (shape: Nx3).
        projected (bool): If True, the first coordinate of points and vertices is excluded.
        threshold (float): Distance threshold for matching points to vertices.

    Returns:
        ndarray: An array of edge indices, each pair corresponding to a matching edge in V.
    """
    points_layer = pointlayer.data
    num_points = len(points_layer)
    
    # Extract edges as pairs of points from the pointlayer
    edges = []
    for i in range(1, num_points, 2):
        edges.append([points_layer[i], points_layer[i+1]])

    # Optionally project points and vertices (remove first coordinate)
    if projected:
        edges = np.array(edges)[:,:,1:]
        V = V[:, 1:]  # Do not modify V directly
    
    edge_indices = []
    
    # Iterate over each edge and find the closest vertices
    for edge in edges:
        p0, p1 = edge
        idx0, idx1 = -1, -1  # Initialize indices

        for i, el in enumerate(V):
            # Find the vertex closest to p0 and p1
            if np.linalg.norm(el - p0) < threshold:
                idx0 = i
            elif np.linalg.norm(el - p1) < threshold:
                idx1 = i

        if idx0 != -1 and idx1 != -1:  # Ensure valid matching
            edge_indices.append([idx0, idx1])

    return np.array(edge_indices)
  
##############################################################################
# https://github.com/matthew-brett/transforms3d/blob/main/transforms3d/affines.py
def decompose44(A44):
    ''' Decompose 4x4 homogenous affine matrix into parts.

    The parts are translations, rotations, zooms, shears.

    This is the same as :func:`decompose` but specialized for 4x4 affines.

    Decomposes `A44` into ``T, R, Z, S``, such that::

       Smat = np.array([[1, S[0], S[1]],
                        [0,    1, S[2]],
                        [0,    0,    1]])
       RZS = np.dot(R, np.dot(np.diag(Z), Smat))
       A44 = np.eye(4)
       A44[:3,:3] = RZS
       A44[:-1,-1] = T

    The order of transformations is therefore shears, followed by
    zooms, followed by rotations, followed by translations.

    This routine only works for shape (4,4) matrices

    Parameters
    ----------
    A44 : array shape (4,4)

    Returns
    -------
    T : array, shape (3,)
       Translation vector
    R : array shape (3,3)
        rotation matrix
    Z : array, shape (3,)
       Zoom vector.  May have one negative zoom to prevent need for negative
       determinant R matrix above
    S : array, shape (3,)
       Shear vector, such that shears fill upper triangle above
       diagonal to form shear matrix (type ``striu``).

    Examples
    --------
    >>> T = [20, 30, 40] # translations
    >>> R = [[0, -1, 0], [1, 0, 0], [0, 0, 1]] # rotation matrix
    >>> Z = [2.0, 3.0, 4.0] # zooms
    >>> S = [0.2, 0.1, 0.3] # shears
    >>> # Now we make an affine matrix
    >>> A = np.eye(4)
    >>> Smat = np.array([[1, S[0], S[1]],
    ...                  [0,    1, S[2]],
    ...                  [0,    0,    1]])
    >>> RZS = np.dot(R, np.dot(np.diag(Z), Smat))
    >>> A[:3,:3] = RZS
    >>> A[:-1,-1] = T # set translations
    >>> Tdash, Rdash, Zdash, Sdash = decompose44(A)
    >>> np.allclose(T, Tdash)
    True
    >>> np.allclose(R, Rdash)
    True
    >>> np.allclose(Z, Zdash)
    True
    >>> np.allclose(S, Sdash)
    True

    Notes
    -----
    The implementation inspired by:

    *Decomposing a matrix into simple transformations* by Spencer
    W. Thomas, pp 320-323 in *Graphics Gems II*, James Arvo (editor),
    Academic Press, 1991, ISBN: 0120644819.

    The upper left 3x3 of the affine consists of a matrix we'll call
    RZS::

       RZS = R * Z *S

    where R is a rotation matrix, Z is a diagonal matrix of scalings::

       Z = diag([sx, sy, sz])

    and S is a shear matrix of form::

       S = [[1, sxy, sxz],
            [0,   1, syz],
            [0,   0,   1]])

    Running all this through sympy (see 'derivations' folder) gives
    ``RZS`` as ::

       [R00*sx, R01*sy + R00*sx*sxy, R02*sz + R00*sx*sxz + R01*sy*syz]
       [R10*sx, R11*sy + R10*sx*sxy, R12*sz + R10*sx*sxz + R11*sy*syz]
       [R20*sx, R21*sy + R20*sx*sxy, R22*sz + R20*sx*sxz + R21*sy*syz]

    ``R`` is defined as being a rotation matrix, so the dot products between
    the columns of ``R`` are zero, and the norm of each column is 1.  Thus
    the dot product::

       R[:,0].T * RZS[:,1]

    that results in::

       [R00*R01*sy + R10*R11*sy + R20*R21*sy + sx*sxy*R00**2 + sx*sxy*R10**2 + sx*sxy*R20**2]

    simplifies to ``sy*0 + sx*sxy*1`` == ``sx*sxy``.  Therefore::

       R[:,1] * sy = RZS[:,1] - R[:,0] * (R[:,0].T * RZS[:,1])

    allowing us to get ``sy`` with the norm, and sxy with ``R[:,0].T *
    RZS[:,1] / sx``.

    Similarly ``R[:,0].T * RZS[:,2]`` simplifies to ``sx*sxz``, and
    ``R[:,1].T * RZS[:,2]`` to ``sy*syz`` giving us the remaining
    unknowns.
    '''
    A44 = np.asarray(A44)
    T = A44[:-1,-1]
    RZS = A44[:-1,:-1]
    # compute scales and shears
    M0, M1, M2 = np.array(RZS).T
    # extract x scale and normalize
    sx = math.sqrt(np.sum(M0**2))
    M0 /= sx
    # orthogonalize M1 with respect to M0
    sx_sxy = np.dot(M0, M1)
    M1 -= sx_sxy * M0
    # extract y scale and normalize
    sy = math.sqrt(np.sum(M1**2))
    M1 /= sy
    sxy = sx_sxy / sx
    # orthogonalize M2 with respect to M0 and M1
    sx_sxz = np.dot(M0, M2)
    sy_syz = np.dot(M1, M2)
    M2 -= (sx_sxz * M0 + sy_syz * M1)
    # extract z scale and normalize
    sz = math.sqrt(np.sum(M2**2))
    M2 /= sz
    sxz = sx_sxz / sx
    syz = sy_syz / sy
    # Reconstruct rotation matrix, ensure positive determinant
    Rmat = np.array([M0, M1, M2]).T
    if np.linalg.det(Rmat) < 0:
        sx *= -1
        Rmat[:,0] *= -1
    return T, Rmat, np.array([sx, sy, sz]), np.array([sxy, sxz, syz])
