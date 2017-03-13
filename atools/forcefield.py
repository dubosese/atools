import numpy as np
import mdtraj as md

def load_surface_BKS(traj_file,top_file,bottom=True):
    surface_atomids = [6,7,8,14,15,16]
    mass_dict = {1:12.0107, 2:12.0107, 3:1.0079, 4:15.9994, 5:1.0079,
                 6:15.9994, 7:28.0855, 8:28.0855, 9:28.0855, 10:12.0107,
                 11:15.9994, 12:15.9994, 13:1.0079, 14:28.0855, 15:28.0855,
                 16:15.9994}
    element_dict = {1:'carbon', 2:'carbon', 3:'hydrogen', 4:'oxygen', 5:'hydrogen',
                    6:'oxygen', 7:'silicon', 8:'silicon', 9:'silicon', 10:'carbon',
                    11:'oxygen', 12:'oxygen', 13:'hydrogen', 14:'silicon',
                    15:'silicon', 16:'oxygen'}
    first = md.load_frame(traj_file,0,top=top_file)
    if bottom:
        surface_sel = "index < {} and {}{}{}".format(int(first.top.n_atoms/2),"(name '","' or name '".join(map(str,surface_atomids)),"')")
    else:
        surface_sel = "index >= {} and {}{}{}".format(int(first.top.n_atoms/2),"(name '","' or name '".join(map(str,surface_atomids)),"')")
    surface_ids = first.top.select(surface_sel)

    traj = md.load(traj_file,top=top_file,atom_indices=surface_ids)
    for atom in traj.top.atoms:
        atom.element = eval('md.core.element.'+element_dict[int(atom.name)])

    return traj
