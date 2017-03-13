from __future__ import division

import numpy as np
import mdtraj as md

def center_surface(traj):
    """ """
    top_atom = np.argmax(traj.xyz[0][:,2])
    init = traj.xyz[0][top_atom]
    for frame in traj:
        diff = init - frame.xyz[0][top_atom]
        for i,dim in enumerate(diff):
            if dim > frame.unitcell_lengths[0][i]/2:
                diff -= frame.unitcell_lengths[0][i]
            elif dim < -frame.unitcell_lengths[0][i]/2:
                diff += frame.unitcell_lengths[0][i]
        frame.xyz[0] += diff

    return traj

def calc_deformation(traj):
    """ """
    rmsd = md.rmsd(traj,traj)
    return rmsd
