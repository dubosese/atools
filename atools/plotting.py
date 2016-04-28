from scipy import integrate
from scipy.stats.kde import gaussian_kde

def smooth_g_r(r, g_r, c_factor=0.15):
    vals = []
    for i,a in enumerate(g_r):
        for b in range(int(a)):
            vals.append(r[i])

    pdf = gaussian_kde(vals)
    pdf.covariance_factor = lambda : c_factor
    pdf._compute_covariance()
    pdf = pdf(r)*integrate.simps(g_r,r)
    return pdf
