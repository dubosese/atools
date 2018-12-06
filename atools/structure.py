from __future__ import division

import mdtraj as md
from mdtraj.core.element import get_by_symbol
from mdtraj.geometry.order import _compute_director
from more_itertools import windowed
import numpy as np
from scipy.integrate import simps
from scipy.stats import binned_statistic_2d

from atools.fileio import read_ndx


def _gather_chains(top_filename, ndx_filename, n_chains, C_only=False):
    """Helper function to gather chain indices into separate arrays

    Parameters
    ----------
    top_filename : str
        Name of topology file (typically GRO format)
    ndx_filename : str
        Name of the GROMACS index file to read group information from
    n_chains : int
        Number of monolayer chains per surface
    C_only : bool, optional, default=False
        Only include carbon atoms in the calculation

    Returns
    -------
    list of lists
        Chain indices in the bottom monolayer
    list of lists
        Chain indices in the top monolayer
    """
    topo = md.load(top_filename).topology
    atoms = np.array([atom.name for atom in topo.atoms])
    groups = read_ndx(ndx_filename)

    if C_only:
        bottom_chains = [id for id in groups['bottom_chains'] if atoms[id] == 'C']
        top_chains = [id for id in groups['top_chains'] if atoms[id] == 'C']
    else:
        bottom_chains = [id for id in groups['bottom_chains']]
        top_chains = [id for id in groups['top_chains']]

    bottom_chains = [chain.tolist() for chain in
                     np.array_split(bottom_chains, n_chains)]
    top_chains = [chain.tolist() for chain in
                  np.array_split(top_chains, n_chains)]

    return bottom_chains, top_chains


def calc_nematic_order(traj_filename, top_filename, ndx_filename, output_filename,
                       n_chains, C_only=False):
    """Calculate the nematic order of a monolayer

    Outputs the average nematic order at each frame of a trajectory
    between the top and bottom monolayers in a dual monolayer system.
    Output is directed to a file with the filename specified by the user.

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file (typically XTC format)
    top_filename : str
        Name of topology file (typically GRO format)
    ndx_filename : str
        Name of the GROMACS index file to read group information from
    output_filename : str
        Name of file to output results to
    n_chains : int
        Number of monolayer chains per surface
    C_only : bool, optional, default=False
        Only include carbon atoms in the calculation

    Notes
    -----
    Only valid for single-component monolayers (assumes that all chains in a
    monolayer contain the same number of atoms).
    Only valid for systems featuring the same number of chains in each monolayer.

    To-do
    -----
    Support for multi-component monolayers
    Support for systems where each monolayer features an arbitrary number of chains

    """
    bottom_chains, top_chains = _gather_chains(top_filename, ndx_filename, n_chains,
                                               C_only=C_only)

    traj = md.load(traj_filename, top=top_filename)
    S2_bottom = md.compute_nematic_order(traj, indices=bottom_chains)
    S2_top = md.compute_nematic_order(traj, indices=top_chains)

    np.savetxt(output_filename, np.column_stack((traj.time, S2_bottom, S2_top)),
        header='Time\tBottom\tTop')


def _tilt_angle(director, normal):
    """Helper function to calculate the tilt angle of a chain

    Parameters
    ----------
    director : array-like, shape=(3,)
        Director vector for a chain
    normal : array-like, shape=(3,)
        Vector normal to the substrate

    Returns
    -------
    float
        Angle between the director and normal vectors (the tilt angle)
    """
    cosang = np.dot(director, normal)
    sinang = np.linalg.norm(np.cross(director, normal))
    angle = np.arctan2(sinang, cosang) * (180.0 / np.pi)
    '''
    If the angle is greater than 90, then our chain director has been drawn in
    the opposite direction. This is corrected by subtracting 180 from the angle.
    '''
    if angle > 90:
        angle -= 180
    return abs(angle)


def calc_avg_tilt_angle(traj_filename, top_filename, ndx_filename, output_filename,
                        n_chains, C_only=False):
    """Calculate the average tilt angle of a monolayer

    Outputs the average tilt angle at each frame of a trajectory,
    averaged between the top and bottom monolayers in a dual monolayer 
    system. Output is directed to a file with the filename specified by
    the user.

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file (typically XTC format)
    top_filename : str
        Name of topology file (typically GRO format)
    ndx_filename : str
        Name of the GROMACS index file to read group information from
    output_filename : str
        Name of file to output results to
    n_chains : int
        Number of monolayer chains per surface
    C_only : bool, optional, default=False
        Only include carbon atoms in the calculation

    Notes
    -----
    Only valid for single-component monolayers (assumes that all chains in a
    monolayer contain the same number of atoms).
    Only valid for systems featuring the same number of chains in each monolayer.

    To-do
    -----
    Support for multi-component monolayers
    Support for systems where each monolayer features an arbitrary number of chains

    """
    bottom_chains, top_chains = _gather_chains(top_filename, ndx_filename, n_chains,
                                               C_only=C_only)

    traj = md.load(traj_filename, top=top_filename)

    tilt_means = []
    tilt_errs = []
    for monolayer, normal in zip([bottom_chains, top_chains],
                                 [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]):
        monolayer_tilt = []
        for i, chain in enumerate(monolayer):
            chain_slice = traj.atom_slice(chain)
            directors = _compute_director(chain_slice)
            monolayer_tilt.append([_tilt_angle(director, normal)
                                   for director in directors])
        monolayer_tilt = np.transpose(monolayer_tilt)
        tilt_means.append(np.mean(monolayer_tilt, axis=1))
        tilt_errs.append(np.std(monolayer_tilt, axis=1))

    tilt_means = np.transpose(tilt_means)
    tilt_errs = np.transpose(tilt_errs)

    np.savetxt(output_filename, np.column_stack((traj.time, tilt_means, tilt_errs)),
        header='Time\tBottom\tTop\tBottom_std\tTop_std')


def count_hydrogen_bonds(traj_filename, top_filename, ndx_filename,
                         mol2_filename, top_group, bottom_group, output_filename):
    """Count the number of inter- or intra-monolayer hydrogen bonds

    The number of inter- and intra-monolayer hydrogen bonds are determined for each
    frame in a trajectory using the Wernet-Nilsson method implemented in MDTraj and
    are output to a file with a user-specified filename.

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file (typically XTC format)
    top_filename : str
        Name of topology file (typically GRO format)
    ndx_filename : str
        Name of the GROMACS index file to read group information from
    mol2_filename : str
        Name of mol2 file used to read bond information
    top_group : str
        Tag for index group (read from NDX file) for indices considered as part
        of the top monolayer
    bottom_group : str
        Tag for index group (read from NDX file) for indices considered as part
        of the bottom monolayer
    output_filename : str
        Name of file to output results to

    """
    topo = md.load(top_filename).topology
    atoms = list(topo.atoms)
    groups = read_ndx(ndx_filename)

    bottom_monolayer = np.array(groups[bottom_group])
    top_monolayer = np.array(groups[top_group])
    monolayers = np.hstack((bottom_monolayer, top_monolayer))

    h_bonds_top = []
    h_bonds_bottom = []
    h_bonds_interface = []
    time_traj = []
    for traj_chunk in md.iterload(traj_filename, top=mol2_filename, chunk=10):
        for i, atom in enumerate(traj_chunk.top.atoms):
            atom.element = atoms[i].element

        h_bonds_total = md.wernet_nilsson(traj_chunk)
        for frame in h_bonds_total:
            h_bonds_top_frame = []
            h_bonds_bottom_frame = []
            h_bonds_interface_frame = []
            for bond in frame:
                if all(atom_id in top_monolayer for atom_id in bond):
                    h_bonds_top_frame.append(tuple(bond))
                elif all(atom_id in bottom_monolayer for atom_id in bond):
                    h_bonds_bottom_frame.append(tuple(bond))
                elif all(atom_id in monolayers for atom_id in bond):
                    h_bonds_interface_frame.append(tuple(bond))
            h_bonds_top.append(len(h_bonds_top_frame))
            h_bonds_bottom.append(len(h_bonds_bottom_frame))
            h_bonds_interface.append(len(h_bonds_interface_frame))
        time_traj.append(traj_chunk.time)

    time_traj = np.concatenate(time_traj).ravel()
    h_bonds_top = np.asarray(h_bonds_top)
    h_bonds_bottom = np.asarray(h_bonds_bottom)
    h_bonds_interface = np.asarray(h_bonds_interface)
    np.savetxt(output_filename,
        np.column_stack((time_traj, h_bonds_interface, h_bonds_top, h_bonds_bottom)),
        header='Time\tInter-\tIntra-top\tIntra-bottom')


def _count_gauche(dihedral_angles):
    if len(dihedral_angles) == 0:
        return 0
    else:
        positive_angles = np.array([angle if angle > 0 else np.pi*2 + angle
            for angle in dihedral_angles])
        small_dihedrals = np.where((positive_angles < np.pi/2))[0]
        large_dihedrals = np.where((positive_angles > np.pi*(3/2)))[0]
        return len(small_dihedrals) + len(large_dihedrals)


def calc_gauche_defects(traj_filename, top_filename, ndx_filename, output_filename,
                        n_chains):
    """Calculate the average number of gauche defects per chain

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file (typically XTC format)
    top_filename : str
        Name of topology file (typically GRO format)
    ndx_filename : str
        Name of the GROMACS index file to read group information from
    output_filename : str
        Name of file to output results to
    n_chains : int
        Number of monolayer chains per surface

    Notes
    -----
    Only valid for single-component monolayers (assumes that all chains in a
    monolayer contain the same number of atoms).
    Only valid for systems featuring the same number of chains in each monolayer.

    """
    bottom_chains, top_chains = _gather_chains(top_filename, ndx_filename, n_chains,
                                               C_only=True)

    traj = md.load(traj_filename, top=top_filename)
    
    gauche_defects = []
    for monolayer in [bottom_chains, top_chains]:
        dihedrals = np.array([list(windowed(chain, n=4)) for chain in monolayer])
        dihedrals = dihedrals.reshape(dihedrals.shape[0]*dihedrals.shape[1], 4)
        dihedral_angles = md.compute_dihedrals(traj, dihedrals, periodic=False)
        gauche_defects.append([_count_gauche(angles) for angles in dihedral_angles])
    gauche_defects = np.transpose(gauche_defects)

    np.savetxt(output_filename, np.column_stack((traj.time, gauche_defects)),
        header='Time\tBottom\tTop')


def _film_level(zcoords, top=False):
    n_bins = int((np.max(zcoords) - np.min(zcoords)) / 0.01)
    if n_bins == 0:
        return float('nan')
    else:
        hist, bin_edges = np.histogram(zcoords, bins=n_bins, normed=False)
        cdist = np.cumsum(hist) / np.sum(hist)
    if top:
        return bin_edges[np.argmax(cdist > 0.05)]
    else:
        return bin_edges[np.argmax(cdist > 0.95)]

def _film_level_bottom(zcoords):
    return _film_level(zcoords)

def _film_level_top(zcoords):
    return _film_level(zcoords, top=True)


def calc_monolayer_roughness(traj_filename, top_filename, ndx_filename,
                             output_filename, bins=5, heavy_only=False):
    """Calculate the surface roughness of the monolayer film

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file (typically XTC format)
    top_filename : str
        Name of topology file (typically GRO format)
    ndx_filename : str
        Name of the GROMACS index file to read group information from
    output_filename : str
        Name of file to output results to
    bins : int, optional, default=5
        Number of bins per surface dimension to pixelate the surface
    heavy_only : bool, optional, default=False
        Exclude hydrogens from the calculation

    """
    topo = md.load(top_filename).topology
    atoms = list(topo.atoms)
    groups = read_ndx(ndx_filename)

    if heavy_only:
        bottom_monolayer = np.array([id for id in groups['bottom']
                                     if atoms[id].name != 'H'])
        top_monolayer = np.array([id for id in groups['top']
                                  if atoms[id].name != 'H'])
    else:
        bottom_monolayer = np.array(groups['bottom'])
        top_monolayer = np.array(groups['top'])

    traj = md.load(traj_filename, top=top_filename)

    xy_plane = [[0.0, traj.unitcell_lengths[0,0]],
                [0.0, traj.unitcell_lengths[0,1]]]

    roughness = []
    for monolayer, statistic in zip([bottom_monolayer, top_monolayer],
                                    [_film_level_bottom, _film_level_top]):
        roughness_monolayer = []
        for frame in traj.xyz:
            monolayer_surface = binned_statistic_2d(
                x=frame[monolayer,0],
                y=frame[monolayer,1],
                values=frame[monolayer,2],
                statistic=statistic,
                bins=bins,
                range=xy_plane)
            roughness_monolayer.append(np.std(monolayer_surface[0]))
        roughness.append(roughness_monolayer)

    roughness = np.transpose(roughness)

    np.savetxt(output_filename, np.column_stack((traj.time, roughness)),
        header='Time\tBottom\tTop')


def calc_interdigitation(traj_filename, top_filename, ndx_filename,
                         output_filename, bin_size=0.025):
    """Calculates the interdigitation between two monolayer films

    This calculation utilizes the procedure of Das et al. (see Ref) to calculate
    an overlap parameter between the two monolayers as a function of their mass
    densities at a particular z-location (normal to the surface). Integrating over
    this parameter from the bottom surface to the top surface yields the
    interdigitation in distance units.

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file (typically XTC format)
    top_filename : str
        Name of topology file (typically GRO format)
    ndx_filename : str
        Name of the GROMACS index file to read group information from
    output_filename : str
        Name of file to output results to
    bin-size : float
        Size of bins in the z-dimension (normal to the surface) to collect mass
        densities. This controls the fidelity of the interdigitation calculation.

    References
    ----------
    .. [1] Das, C., Noro, M.G., Olmsted, P.D., "Simulation studies of stratum
           corneum lipid mixtures." (2009) Biophys. J. 97, 1941-1951

    """

    groups = read_ndx(ndx_filename)
    bottom_monolayer = groups['bottom_chains']
    top_monolayer = groups['top_chains']

    traj = md.load(traj_filename, top=top_filename)
    atoms = np.array(list(traj.top.atoms))
    bottom_masses = np.array([atom.element.mass for atom in atoms[bottom_monolayer]])
    top_masses = np.array([atom.element.mass for atom in atoms[top_monolayer]])

    bounds = (np.min(traj.xyz[:, bottom_monolayer, 2]),
              np.max(traj.xyz[:, top_monolayer, 2]))
    bins = int(round((bounds[1] - bounds[0]) / bin_size))
    bin_size = (bounds[1] - bounds[0]) / bins

    interdigitation = []
    for xyz in traj.xyz:
        hist_bottom, _ = np.histogram(xyz[bottom_monolayer,2], bins=bins,
            range=bounds, normed=False, weights=bottom_masses)
        hist_top, bin_edges = np.histogram(xyz[top_monolayer,2], bins=bins,
            range=bounds, normed=False, weights=top_masses)
        bin_centers = bin_edges[1:] - bin_size/2
        rho_overlap = []
        for count_top, count_bottom in zip(hist_top, hist_bottom):
            count_prod = count_top * count_bottom
            count_sum = count_top + count_bottom
            overlap = 0
            if count_sum != 0:
                overlap = (4 * count_prod) / (count_sum ** 2)
            rho_overlap.append(overlap)
        interdigitation.append(simps(rho_overlap, x=bin_centers))

    np.savetxt(output_filename, np.column_stack((traj.time, interdigitation)))


def _create_com_traj(xyz, time, unitcell_lengths, unitcell_angles, masses=None):
    top = md.Topology()
    chain = top.add_chain()
    for i, com in enumerate(xyz[0]):
        com_res = top.add_residue('COM{}'.format(i), chain)
        top.add_atom('COM', get_by_symbol('C'), residue=com_res)
    return md.Trajectory(xyz, topology=top, time=time,
                         unitcell_lengths=unitcell_lengths,
                         unitcell_angles=unitcell_angles)


def calc_hexagonal_order(traj_filename, top_filename, output_filename,
                         ndx_filename, n_chains):
    """Calculate the hexagonal order of chains and terminal groups

    Returns the average hexagonal order at each frame of a trajectory
    between the top and bottom monolayers in a dual monolayer system.

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file
    ndx_filename : str
        Name of Gromacs .ndx file which specifies atom groups
    n_chains : int
        Number of monolayer chains per surface

    Notes
    -----
    Assumes a unique chain prototype.
    Assumes identical top and bottom monolayers.

    """
    pass
    '''
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = np.array(groups['bottom_chains']) - 1
    top_chains = np.array(groups['top_chains']) - 1
    bottom_chains = np.array_split(bottom_chains, n_chains)
    top_chains = np.array_split(top_chains, n_chains)
    bottom_termini = np.array(groups['bottom_termini']) - 1
    top_termini = np.array(groups['top_termini']) - 1
    bottom_termini = np.array_split(bottom_termini, n_chains)
    top_termini = np.array_split(top_termini, n_chains)

    traj = md.load(traj_filename, top=top_filename)

    # Calc just for bottom chains
    bottom_chain_COM = []
    for chain in bottom_chains:
        chain_slice = traj.atom_slice(chain)
        chain_com = md.compute_center_of_mass(chain_slice)
        for com in chain_com:
            for i in range(2):
                if com[i] < 0:
                    com[i] += traj.unitcell_lengths[0,i]
                elif com[i] > traj.unitcell_lengths[0,i]:
                    com[i] -= traj.unitcell_lengths[0,i]
        bottom_chain_COM.append(chain_com)
    bottom_chain_COM = np.array(bottom_chain_COM).transpose(1,0,2)

    com_traj = _create_com_traj(bottom_chain_COM, traj.time, traj.unitcell_lengths,
                                traj.unitcell_angles)

    # Create bond order diagram
    for frame in bottom_chain_COM:
        distances = md.compute_distances()

    for frame in bottom_BOD:
        c_op = complex()
        for com in frame:
            theta = np.arctan2(com[1], com[0]) + np.pi
            c_op += complex(np.cos(6. * theta), -sin(6. * theta))
        c_op = abs(c_op / len(frame)).real
    '''

def calc_hexagonal_order_freud(traj_filename, top_filename, output_filename,
                               ndx_filename, n_chains):
    from freud.box import Box
    from freud.order import HexOrderParameter

    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = np.array(groups['bottom_chains']) - 1
    top_chains = np.array(groups['top_chains']) - 1
    bottom_chains = np.array_split(bottom_chains, n_chains)
    top_chains = np.array_split(top_chains, n_chains)
    bottom_termini = np.array(groups['bottom_termini']) - 1
    top_termini = np.array(groups['top_termini']) - 1
    bottom_termini = np.array_split(bottom_termini, n_chains)
    top_termini = np.array_split(top_termini, n_chains)

    traj = md.load(traj_filename, top=top_filename)

    for i, group in enumerate([bottom_chains, top_chains, bottom_termini,
            top_termini]):
        group_COM = []
        for chain in group:
            chain_slice = traj.atom_slice(chain)
            chain_com = md.compute_center_of_mass(chain_slice)
            for com in chain_com:
                for j in range(2):
                    if com[j] < 0:
                        com[j] += traj.unitcell_lengths[0,j]
                    elif com[j] > traj.unitcell_lengths[0,j]:
                        com[j] -= traj.unitcell_lengths[0,j]
            group_COM.append(chain_com)
        group_COM = np.array(group_COM).transpose(1,0,2)

        rmax = 0.75
        box = Box(Lx=traj.unitcell_lengths[0,0], Ly=traj.unitcell_lengths[0,1],
                  Lz=traj.unitcell_lengths[0,2])

        order_parameter = []
        for xyz in group_COM:
            median = np.median(xyz[:,2])
            if i == 2:
                xyz = np.array([point for point in xyz if point[2] > median - 0.5])
            elif i == 3:
                xyz = np.array([point for point in xyz if point[2] < median + 0.5])
            hex_order = HexOrderParameter(rmax, 6)
            hex_order.compute(box, xyz.astype(np.float32))
            order_parameter.append(abs(hex_order.getPsi()).mean())
        print(np.mean(order_parameter[int(len(order_parameter)/2):]))
    import pdb;pdb.set_trace()

def identify_rigid_groups(monolayer, terminal_group=None, freeze_thickness=5):
    bounding_box = monolayer.boundingbox
    bot_of_box = bounding_box.mins[2]
    top_of_box = bounding_box.maxs[2]

    full_system = []
    bot_rigid = []
    top_rigid = []
    if terminal_group:
        nonterminal = []
        terminal = []
    for i, particle in enumerate(monolayer.particles()):
        full_system.append(i + 1)
        if particle.pos[2] < bot_of_box + freeze_thickness:
            bot_rigid.append(i + 1)
        elif particle.pos[2] > top_of_box - freeze_thickness:
            top_rigid.append(i + 1)
        if terminal_group:
            if terminal_group.title() in [parent.name for parent in particle.ancestors()]:
                terminal.append(i + 1)
            else:
                nonterminal.append(i + 1)

    rigid_groups = {'System': full_system,
                    'bot': bot_rigid,
                    'top': top_rigid}
    if terminal_group:
        rigid_groups['nonterminal'] = nonterminal
        rigid_groups['terminal'] = terminal
    return rigid_groups


def identify_rigid_groups_single(monolayer, terminal_group=None, freeze_thickness=5):
    bounding_box = monolayer.boundingbox
    bot_of_box = bounding_box.mins[2]

    full_system = []
    bot_rigid = []
    chains = []
    if terminal_group:
        nonterminal = []
        terminal = []
    for i, particle in enumerate(monolayer.particles()):
        full_system.append(i + 1)
        if particle.pos[2] < bot_of_box + freeze_thickness:
            bot_rigid.append(i + 1)
        if 'SilicaInterface' in [parent.name for parent in particle.ancestors()]:
            chains.append(i + 1)
        if terminal_group:
            if terminal_group.title() in [parent.name for parent in particle.ancestors()]:
                terminal.append(i + 1)
            else:
                nonterminal.append(i + 1)

    rigid_groups = {'System': full_system,
                    'bot': bot_rigid,
                    'chains': chains}
    if terminal_group:
        rigid_groups['nonterminal'] = nonterminal
        rigid_groups['terminal'] = terminal
    return rigid_groups
