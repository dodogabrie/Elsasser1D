"""Render data2d/w*.bin (Elsasser vorticities w+, w-) to PNG frames: current density j + energy spectrum."""
import glob, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
BG = "#11151c"
from matplotlib.colors import LinearSegmentedColormap
CM = LinearSegmentedColormap.from_list("sobria", ["#5aa9e6", "#1a2230", "#f0a95f"])
files = sorted(glob.glob("data2d/w*.bin")); n = int(sys.argv[1]); dt_out = float(sys.argv[2])
k1 = np.fft.fftfreq(n, 1 / n); kx, ky = np.meshgrid(k1, k1, indexing="ij"); k2 = kx**2 + ky**2; k2[0, 0] = 1
kb = np.rint(np.sqrt(k2)).astype(int); ks = np.arange(1, n // 2)
def load(f):
    a = np.fromfile(f, np.float32).astype(float).reshape(2, n, n); return a[0].T, a[1].T  # Fortran (i,j) -> [i,j]
def spec(wp, wm):
    # Z+- = zhat x grad(zeta), |Z|^2 spectrum = |w|^2/k^2; E = (|Z+|^2 + |Z-|^2)/4
    e = (np.abs(np.fft.fft2(wp))**2 + np.abs(np.fft.fft2(wm))**2) / k2 / n**4 / 4
    e[0, 0] = 0; return np.bincount(kb.ravel(), e.ravel())[ks]
frames = [load(f) for f in files]
js = [(a - b) / 2 for a, b in frames]; lim = np.percentile(np.abs(js[-1]), 99.5)
E = [spec(*fr) for fr in frames]; avg = np.mean(E[-10:], axis=0)
for i, ((wp, wm), j) in enumerate(zip(frames, js)):
    fig = plt.figure(figsize=(10, 5), dpi=100, facecolor=BG)
    ax = fig.add_axes([0.01, 0.02, 0.47, 0.94]); ax.imshow(j.T, origin="lower", interpolation="bicubic", cmap=CM, vmin=-lim, vmax=lim); ax.axis("off")
    ax.set_title(f"current density j,  t = {i*dt_out:.2f}", color="#c9d1d9", fontsize=11)
    ax2 = fig.add_axes([0.58, 0.14, 0.38, 0.74], facecolor=BG)
    ax2.loglog(ks, E[i], color="#e8a05a", lw=1.8); ax2.loglog(ks, avg, color="#58627a", lw=1, ls=":")
    ref = ks[(ks > 5) & (ks < 40)]; r0 = 1.2 * avg[6] * (ref / 7.) ** -1.5
    ax2.loglog(ref, r0, color="#7aa2d6", lw=1, ls="--"); ax2.text(12, 2 * r0[7], "k^-3/2", color="#7aa2d6", fontsize=9)
    ax2.set_ylim(1e-9, 1e-1); ax2.set_xlim(1, n // 2)
    ax2.set_title("E(k)", color="#c9d1d9", fontsize=11); ax2.set_xlabel("k", color="#8b95a5")
    ax2.tick_params(colors="#8b95a5", labelsize=8)
    for s in ax2.spines.values(): s.set_color("#2a3140")
    fig.savefig(f"png/f{i:03d}.png", facecolor=BG); plt.close(fig)
