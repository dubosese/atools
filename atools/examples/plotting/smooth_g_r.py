import mdtraj as md
import matplotlib.pyplot as plt

traj = md.load('mono_lj.mol2')
pairs = traj.top.select_pairs("all","all")
r,g_r = md.compute_rdf(traj,pairs)

plt.plot(r,g_r)
plt.show()
