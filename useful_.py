def AABB_on_atoms(ase_atoms, points=None):

    from freud.locality import AABBQuery
    from freud.box import Box
    import ase 
    # from ase.data import vdw_radii


    if points is None:
        points = ase_atoms.get_positions()
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T # this needs to be an upper triangular matrix, essential column vectors, not row vectors
    
    box = Box.from_matrix(cell, dimensions =3) # Define the box for Freud
    ab  = AABBQuery(box=box, points=points) # Compute the AABB Tree
    # Perform the nearest neighbor query with n-neighbors (10 by default)
        
    return ab
def get_wall_windows2(regions, maxima, dgrid, ase_atoms, wall_thickness=2):

    import numpy as np
    from icecream import ic
    wall_windows = []
    shape =regions.shape
    dinterp = interpolate_labels(dgrid)
    for a in [0,1,2]:
        ic.disable()
        ic(a)
        llr = regions.take(indices=0 ,axis=a) #* left wall
        rrr = regions.take(indices=shape[a]-1 ,axis=a)  #* right wall

        map_list = np.unique(list(zip(llr.flatten(), rrr.flatten())), axis=0)
        map_list = np.vstack([m for m in map_list if np.all(m!=0)])
        flags = np.vstack([[maxima[k - 1][a]>wall_thickness, maxima[v-1][a]<(shape[a]-wall_thickness)] for k,v in map_list])    
        flags = np.logical_and(flags[:,0],flags[:,1])
        # window_indices = [np.vstack(np.where(llr==m[0])).mean(axis=1) for i,m in enumerate(map_list) if flags[i]]
        wregs_l  = [(np.array([m[0], m[1]]), True,'periodic') for i, m in enumerate(map_list) if flags[i]]
        wregs_r  = [(np.array([m[1], m[0]]), True,'periodic') for i, m in enumerate(map_list) if flags[i]]
        window_indices = [np.vstack(np.where(llr==m[0])).mean(axis=1).T for i, m in enumerate(map_list) if flags[i]]
        ic(window_indices)
        if np.any(flags):
            window_indices_l = np.insert(np.vstack(window_indices), [a], np.zeros((len(window_indices),1)), axis=1  )
            window_indices_r = np.insert(np.vstack(window_indices), [a], (shape[a] -1)*np.ones((len(window_indices),1)), axis=1 )
            window_radii_l = np.vstack(dinterp(window_indices_l/regions.shape))
            window_radii_r = np.vstack(dinterp(window_indices_r/regions.shape))
            window_coords_l = find_coord_from_indices(window_indices_l, shape=shape, ase_atoms=ase_atoms)
            window_coords_r = find_coord_from_indices(window_indices_r, shape=shape, ase_atoms=ase_atoms)
            ic(window_coords_l)
            ic(window_coords_r)
            window_rows_l = np.vstack([[*w, np.array([ window_coords_l[i] ]), window_radii_l[i]] for i,w in enumerate(wregs_l)])
            window_rows_r = np.vstack([[*w, np.array([ window_coords_r[i] ]), window_radii_r[i]] for i,w in enumerate(wregs_r)])
            wall_windows.append(window_rows_l) 
            wall_windows.append(window_rows_r) 
            
    if len(wall_windows)>0:
        wall_windows =  np.vstack(wall_windows)
    return wall_windows
def get_connections_for_rag(rag,region_labels, maxima,ase_atoms,dgrid, minimum_window_separation=1.0, wall_windows=True, wall_thickness=2):
    

    def find_windows(regs):
        

            import numpy as np
            from sklearn.cluster import AgglomerativeClustering
            i = regs[0]
            j = regs[1]
            check_flag = np.logical_or(np.logical_and(region_labels== i, outer.get(j)), np.logical_and(region_labels== j, outer.get(i)))
            if np.sum(check_flag) > 0:
                # print("Edges found between "+str(i)+' and '+ str(j) + ' sum: '+str(np.sum(check_flag)))
                windices = np.vstack(np.where(check_flag)).T  # * (N,3)
                window_points = np.dot(A_unit, (windices /region_labels.shape).T).T  # * (N,3)
                # ic(window_points)
                # * No clustering
                # window_centers = np.mean(window_points, axis=0)
                
                # * DBScan
                # X_scaled = StandardScaler().fit_transform(window_points)
                # db = DBSCAN(eps=0.5, min_samples=10).fit(X_scaled)
                # window_centers = np.vstack([window_points[db.labels_==i].mean(axis=0) for i in np.unique(db.labels_) if i != -1])
                
                # if len(window_points) >1:
                # * Agglomerative
                ag = AgglomerativeClustering(n_clusters=None, linkage='single',distance_threshold=minimum_window_separation, compute_full_tree=True).fit(window_points)
                window_centers = np.vstack([window_points[ag.labels_==i].mean(axis=0) for i in np.unique(ag.labels_) if i != -1])
            # else:
                # window_centers = window_points
                window_fractional = get_fractional_coordinates(window_centers,ase_atoms)
                window_radii = np.array(dinterp(window_fractional))
                return np.array([regs, True,'internal', window_centers, window_radii], dtype='object') # returns a list of [[r1,r2], connected or not, array of window_centers]
            else:
                return np.array([regs, False,'internal', None, None], dtype='object')

    #* find the internal windows
    from skimage.segmentation import find_boundaries
    import ase
    import numpy as np
    from rich.progress import track
    connected_regions = np.vstack(list(rag.edges()))
    list_of_region_indices = np.delete(np.unique(region_labels),0)
    outer = {reg:find_boundaries(region_labels== reg, mode='outer') for reg in list_of_region_indices}
    dinterp = interpolate_labels(dgrid)
    minimum_window_separation = 1.0
    A_unit = ase.geometry.complete_cell(ase_atoms.get_cell()).T    

    connections= [find_windows(regs) for regs in track(rag.edges, description='Finding windows:')]

    #* Windows on the wall are to be calculated (periodic windows)
    if wall_windows:
        wwl = get_wall_windows2(regions=region_labels, maxima=maxima, dgrid=dgrid, ase_atoms=ase_atoms, wall_thickness=wall_thickness)
        if len(wwl)>0:
            connections = connections+wwl.tolist()
    
    return connections

# * Assign ions to cages
def find_coord_from_indices(indices, shape, ase_atoms):

    import ase
    import numpy as np
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T
    return np.dot(cell, (indices/shape).T).T
def get_cage_for_ions(coords, cage_labels, interpolator, ase_atoms):

    import numpy as np
    output_labels = []
    for coord in coords:
        # print(coords.shape)
        # print((interpolator(get_fractional_coordinates(coords, ase_atoms=nax)).astype(int)-1).shape)
        try:
            output_labels.append(np.array(cage_labels)[interpolator(mgr.get_fractional_coordinates(coord, ase_atoms=ase_atoms)).astype(int)-1])
        except:
            print(np.max(coord))
            print(np.min(coord))
            pass
    return np.vstack(output_labels)
def gpd_aabb(grid_points, ase_atoms, radii, n_neighbors=10, probe_radius=0.0):

    from freud.locality import AABBQuery
    from freud.box import Box
    import ase 
    # from ase.data import vdw_radii
    from icecream import ic
    import pandas as pd
    import numpy as np

    framework_atom_positions = ase_atoms.get_positions()
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T # this needs to be an upper triangular matrix, essential column vectors, not row vectors
    
    box = Box.from_matrix(cell, dimensions =3) # Define the box for Freud
    ab  = AABBQuery(box=box, points=framework_atom_positions) # Compute the AABB Tree
    # Perform the nearest neighbor query with n-neighbors (10 by default)
    ab_r = ab.query(grid_points, query_args=dict(mode='nearest', num_neighbors=n_neighbors, exclude_ii=True))
    ab_r_nl = ab_r.toNeighborList(sort_by_distance=True) #convert query to neighborlist
    nls = ab_r_nl.point_indices.reshape(-1,n_neighbors)  # get the indices of the neighbors for each grid pt. to lookup the radii later
    dists =ab_r_nl.distances.reshape(-1,n_neighbors) # get the distances to the neighbors of each grid pt.
    # print(ab, ab_r, ab_r_nl, nls, [radii[nl] for nl in nls])#, np.min(dists - [radii[nl] for nl in nls],axis=1))
    # return the minimum (shortest distance to the surface) in each row
    return np.min(dists - [radii[nl] for nl in nls],axis=1)

def dask_grid_over_atoms(ase_atoms, spacing=0.1, chunksize=50000):

    import ase
    import dask.array as da
    import numpy as np
    #number of grid points in each direction 
    [nx, ny, nz] = (ase_atoms.get_cell_lengths_and_angles()[0:3] / spacing).astype(int) + 1
    gpoints = (da.stack(da.meshgrid(np.linspace(0, 1, nx), np.linspace(0,1, ny), np.linspace(0,1, nz), indexing='ij'),-1).reshape(-1, 3)).rechunk(chunksize,3)
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T # cell matrix
    return da.dot(cell, gpoints.T).T, (nx,ny,nz) #return the actual coordinates

def dgrid_from_atoms(ase_atoms, radii=None, spacing=0.25, block_size=50000, n_neighbors=10, probe_radius=0):
 
    import numpy as np
    import dask.array as da
    from rich.progress import track
    gpoints, shape = dask_grid_over_atoms(ase_atoms, spacing=spacing, chunksize=block_size) #mesh grid
    dgrid = []
    import pandas as  pd 
    if radii is None:
        vdw_radii = pd.Series({'H': 110.00000000000001,'He': 140.0,'Li': 182.0,'Be': 153.0,'B': 192.0,'C': 170.0,'N': 155.0,'O': 152.0,'F': 147.0,'Ne': 154.0,'Na': 227.0,'Mg': 173.0,'Al': 184.0,'Si': 210.0,'P': 180.0,'S': 180.0,'Cl': 175.0,'Ar': 188.0,'K': 275.0,'Ca': 231.0,'Sc': 215.0,'Ti': 211.0,'V': 206.99999999999997,'Cr': 206.0,'Mn': 204.99999999999997,'Fe': 204.0,'Co': 200.0,'Ni': 197.0,'Cu': 196.0,'Zn': 200.99999999999997,'Ga': 187.0,'Ge': 211.0,'As': 185.0,'Se': 190.0,'Br': 185.0,'Kr': 202.0,'Rb': 303.0,'Sr': 249.00000000000003,'Y': 231.99999999999997,'Zr': 223.0,'Nb': 218.00000000000003,'Mo': 217.0,'Tc': 216.0,'Ru': 213.0,'Rh': 210.0,'Pd': 210.0,'Ag': 211.0,'Cd': 218.00000000000003,'In': 193.0,'Sn': 217.0,'Sb': 206.0,'Te': 206.0,'I': 198.0,'Xe': 216.0,'Cs': 343.0,'Ba': 268.0,'La': 243.00000000000003,'Ce': 242.0,'Pr': 240.0,'Nd': 239.0,'Pm': 238.0,'Sm': 236.0,'Eu': 235.0,'Gd': 234.0,'Tb': 233.0,'Dy': 231.0,'Ho': 229.99999999999997,'Er': 229.0,'Tm': 227.0,'Yb': 225.99999999999997,'Lu': 224.00000000000003,'Hf': 223.0,'Ta': 222.00000000000003,'W': 218.00000000000003,'Re': 216.0,'Os': 216.0,'Ir': 213.0,'Pt': 213.0,'Au': 214.0,'Hg': 223.0,'Tl': 196.0,'Pb': 202.0,'Bi': 206.99999999999997,'Po': 197.0,'At': 202.0,'Rn': 220.00000000000003,'Fr': 348.0,'Ra': 283.0,'Ac': 247.00000000000003,'Th': 245.00000000000003,'Pa': 243.00000000000003,'U': 241.0,'Np': 239.0,'Pu': 243.00000000000003,'Am': 244.0,'Cm': 245.00000000000003,'Bk': 244.0,'Cf': 245.00000000000003,'Es': 245.00000000000003,'Fm': 245.00000000000003,'Md': 246.0,'No': 246.0,'Lr': 246.0,'Rf': np.nan,'Db': np.nan,'Sg': np.nan,'Bh': np.nan,'Hs': np.nan,'Mt': np.nan,'Ds': np.nan,'Rg': np.nan,'Cn': np.nan,'Nh': np.nan,'Fl': np.nan,'Mc': np.nan,'Lv': np.nan,'Ts': np.nan,'Og': np.nan})
        radii = vdw_radii[ase_atoms.get_chemical_symbols()].values/100.0
    
    for block in track(gpoints.blocks, description='Computing distance grid', total=gpoints.npartitions):
        dgrid.append(gpd_aabb(block.compute(), ase_atoms, radii, n_neighbors=n_neighbors,probe_radius=probe_radius )) #compute grid
    return da.hstack(dgrid).rechunk(chunks=block_size).reshape(shape) #return in correct shape


def get_rag_from_regions_and_grid(regions, dgrid, ase_atoms, maxima=None):
 
    
    import numpy as np
    # * compute a region adjacency graph from skimage
    from skimage.future import graph 
    gr_mean = graph.rag_mean_color(dgrid, regions)
    
    from skimage.measure import regionprops 
    # * add the centroid property to the rag nodes
    region_props = regionprops(regions)
    for region_prop in region_props:
        #add the coordinates of the centroid to the node
        gr_mean.nodes[region_prop['label']]['centroid'] = find_coord_from_indices( indices = np.array(region_prop['centroid']), ase_atoms=ase_atoms,shape=dgrid.shape)
        gr_mean.nodes[region_prop['label']]['centroid_indices'] = region_prop['centroid'] # * add the indices of the centroid to the node
    gr_mean.remove_node(0)
    if maxima is not None:
        maxima_radii = dgrid[tuple(maxima.T)] # * radii of the maxima
        gr_mean = add_maxima_to_rag(gr_mean, maxima, maxima_radii, regions.shape, ase_atoms)
        for edge in gr_mean.edges():
            gr_mean.edges[edge]['weight'] = np.linalg.norm(gr_mean.nodes[edge[0]]['maxima'] - gr_mean.nodes[edge[1]]['maxima'])
    else:
        for edge in gr_mean.edges():
            gr_mean.edges[edge]['weight'] = np.linalg.norm(gr_mean.nodes[edge[0]]['centroid'] - gr_mean.nodes[edge[1]]['centroid'])
    
    return gr_mean

def regions_from_dgrid(dgrid, mask_thickness=1, h=0.2, min_distance=2, compactness=0):
  
    from scipy import ndimage as ndi
    import numpy as np
    import sparse
    # from .helper import interpolate_labels, find_coord_from_indices
    from skimage import segmentation
    from skimage.morphology import extrema
    from skimage.feature import peak_local_max, corner_peaks
    from sklearn.cluster import AgglomerativeClustering

    from icecream import ic 
    probe = mask_thickness
    import dask.array as da
    import warnings
    # reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int))

    # max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1, **kwargs)).T
    # max_indices = np.array(peak_local_max(dgrid,  **kwargs)).T
    max_indices = np.vstack(np.where(extrema.h_maxima(dgrid,  h=h)))
    markers_sparse = sparse.COO(max_indices, 1, shape=dgrid.shape)

    probe = mask_thickness
    import dask.array as da
    reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int), compactness=compactness)



    # print(keep_these)
    # new_max_indices = max_indices[keep_these-1].T

    new_max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1,p_norm=2)).T
    if len(new_max_indices.T) < len(max_indices.T):
        warnings.warn('Some maxima excluded as they fell inside the mask_thickness')
        # warnings.warn('Some maxima excluded as they the regions were less that 0.001 of the box')
    
    markers_sparse = sparse.COO(new_max_indices, 1, shape=dgrid.shape)

    probe = mask_thickness
    import dask.array as da
    reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int), compactness=compactness)

    new_max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1, p_norm=2)).T
    #
    #* Now interpolate the region labels and sort the maxima in the correct order    
    sort_index = np.argsort([reg_labels[tuple(m)] for m in new_max_indices.T])
    new_max_indices = (new_max_indices.T[sort_index]).T
 

    return reg_labels, new_max_indices.T
def regions_from_dgrid(dgrid, mask_thickness=1, h=0.2, min_distance=2, compactness=0):

    from scipy import ndimage as ndi
    import numpy as np
    import sparse
    # from .helper import interpolate_labels, find_coord_from_indices
    from skimage import segmentation
    from skimage.morphology import extrema
    from skimage.feature import peak_local_max, corner_peaks
    from sklearn.cluster import AgglomerativeClustering

    from icecream import ic 
    probe = mask_thickness
    import dask.array as da
    import warnings
    # reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int))

    # max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1, **kwargs)).T
    # max_indices = np.array(peak_local_max(dgrid,  **kwargs)).T
    max_indices = np.vstack(np.where(extrema.h_maxima(dgrid,  h=h)))
    markers_sparse = sparse.COO(max_indices, 1, shape=dgrid.shape)

    probe = mask_thickness
    import dask.array as da
    reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int), compactness=compactness)



    # print(keep_these)
    # new_max_indices = max_indices[keep_these-1].T

    new_max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1,p_norm=2)).T
    if len(new_max_indices.T) < len(max_indices.T):
        warnings.warn('Some maxima excluded as they fell inside the mask_thickness')
        # warnings.warn('Some maxima excluded as they the regions were less that 0.001 of the box')
    
    markers_sparse = sparse.COO(new_max_indices, 1, shape=dgrid.shape)

    probe = mask_thickness
    import dask.array as da
    reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int), compactness=compactness)

    new_max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1, p_norm=2)).T
 

    #* Now interpolate the region labels and sort the maxima in the correct order    
    sort_index = np.argsort([reg_labels[tuple(m)] for m in new_max_indices.T])
    new_max_indices = (new_max_indices.T[sort_index]).T
 

    return reg_labels, new_max_indices.T

def read_raspa_pdb(path_to_file):


        import numpy as np

        f = open(path_to_file).readlines()
        start = np.where(["MODEL" in line for line in f])[0] + 2  # * Start of config
        ends = np.where(["ENDMDL" in line for line in f])[0]  # End for config
        cryst = np.where(["CRYST" in line for line in f])[0]  # box shape for the config

        data = [f[start[i]:ends[i]] for i in range(len(start))]

        coord = np.array([[np.array(line.split()[4:7]).astype(float) for line in d] for d in data])
        cell_dims = np.array([np.array(line.split())[1:].astype(float) for line in f if "CRYST" in line])
        symbols = np.array([[np.array(line.split()[2]) for line in d] for d in data])

        output = {}
        output['cells'] = cell_dims
        output['coords'] = coord
        output['symbols'] = symbols

        return output

def interpolate_labels(regions):
 
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator as RGI
    xx = np.linspace(0,  1,regions.shape[0])
    yy = np.linspace(0,  1,regions.shape[1])
    zz = np.linspace(0,  1,regions.shape[2])
    rl = RGI((xx, yy, zz), regions, method='nearest')

    return rl
    
def get_fractional_coordinates(points, ase_atoms):
 
    import numpy as np
    import ase
    cell= ase.geometry.complete_cell(ase_atoms.get_cell()).T 
    cell_inv = np.linalg.inv(cell)
    frac_coords = np.dot(cell_inv, points.T).T
    return frac_coords

def add_maxima_to_rag(G, maxima,  maxima_radii, shape, ase_atoms):


    if len(G.nodes()) != len(maxima):
        raise ValueError('The number of nodes in the graph and the number of maxima do not match')
    else:
        # regionprops = get_region_props_for_regions(G.nodes())
        for n, m in zip(G.nodes(), maxima):
            G.nodes[n]['maxima'] = find_coord_from_indices(maxima[G.nodes[n]['labels'][0]-1], shape=shape, ase_atoms=ase_atoms) #* add the coordinates of the maxima to the node
            G.nodes[n]['maxima_indices'] = maxima[G.nodes[n]['labels'][0]-1] #* add the indcies of the maxima 
            G.nodes[n]['maxima_radii'] =maxima_radii[G.nodes[n]['labels'][0]-1] #* add the indcies of the maxima 
    return G