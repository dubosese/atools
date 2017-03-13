import numpy as np

from scipy.spatial import distance

class System(object):
    def __init__(self, xyz, ff, dists):
        super(System, self).__init__()

        self.xyz = xyz
        self.ff = ff
        self.dists = dists
        self.target = target
        self.configs = configs
        self.U = np.zeros(len(dists))
        self.U_err = np.zeros(len(dists))

    def calc_U(self, write_traj=False):
        if write_traj:
            traj_out = open('traj.xyz','w') # Trajectory

        for i,dist in enumerate(self.dists):
            U_bin = np.zeros(self.configs)
            for j in range(self.configs):
                # Rotate nanoparticles randomly about center of mass
                for axis in ['x','y','z']:
                    self.rotate(self.xyz1,random()*2.*pi,axis)
                    self.rotate(self.xyz2,random()*2.*pi,axis)

                # Shift second molecule away from origin
                self.xyz2[:,0] += dist

                # Calculate distances between the two molecules
                particle_dists = np.ravel(distance.cdist(self.xyz1,self.xyz2,'euclidean'))

                # Write to trajectory
                if write_traj and j % 5 == 0:
                    traj_out.write(str(len(self.xyz1)+len(self.xyz2))+'\n'+'\n')
                    for k,coord in enumerate(self.xyz1):
                        traj_out.write('CG\t')
                        for a in range(3):
                            traj_out.write('{}\t'.format(coord[a]))
                        traj_out.write('\n')
                    for k,coord in enumerate(self.xyz2):
                        traj_out.write('CG\t')
                        for a in range(3):
                            traj_out.write('{}\t'.format(coord[a]))
                        traj_out.write('\n')

                # Calculate potential at this state
                U_bin[j] = self.ff.calc_U(particle_dists)

                # Move second molecule back to origin
                self.xyz2[:,0] -= dist

            self.U[i] = np.mean(U_bin)
            self.U_err[i] = np.std(U_bin)

    def calc_com(self, xyz):
        return np.mean(xyz,axis=0)

    def center(self, xyz):
        xyz -= self.calc_com(xyz)

    def rotate(self, xyz, theta, axis):
        assert(axis in ['x','y','z'])
        com = self.calc_com(xyz)
        self.center(xyz)
        xyz_init = deepcopy(xyz)
        for i,coord in enumerate(xyz):
            if axis == 'x':
                coord[1] = xyz_init[i,1]*cos(theta) - xyz_init[i,2]*sin(theta)
                coord[2] = xyz_init[i,1]*sin(theta) + xyz_init[i,2]*cos(theta)
            elif axis == 'y':
                coord[0] = xyz_init[i,0]*cos(theta) + xyz_init[i,2]*sin(theta)
                coord[2] = -xyz_init[i,0]*sin(theta) + xyz_init[i,2]*cos(theta)
            else:
                coord[0] = xyz_init[i,0]*cos(theta) - xyz_init[i,1]*sin(theta)
                coord[1] = xyz_init[i,0]*sin(theta) + xyz_init[i,1]*cos(theta)

        xyz += com # shift back to original center

    def write_U(self, filename="U.txt"):
        np.savetxt(filename,np.column_stack((np.asarray(self.dists),
                                             np.asarray(self.U),
                                             np.asarray(self.U_err))))

    def calc_error(self, norm=True):
        error = sum([abs(self.target[i] - U) for i,U in enumerate(self.U)])
        if norm:
            error /= sum([abs(self.target[i] + U) for i,U in enumerate(self.U)])
        return error


if __name__ == "__main__":
    from cgnano.nanoparticles.CG_nano import CG_nano
    from cgnano.forcefield import Forcefield
    import cProfile, pstats, StringIO

    pr = cProfile.Profile()
    pr.enable()

    nano = CG_nano(3.0,sigma=0.9235)
    ff = Forcefield(9.235,3.657,21.81,5.81)
    aa_data = np.loadtxt('/Users/summeraz/cgnano/cgnano/target_data/np-np/U_4nm.txt')
    system = System(nano.xyz*10.,deepcopy(nano.xyz*10.),ff,aa_data[:,0])
    system.calc_U(write_traj=True)
    system.write_U('U.txt')

    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()
