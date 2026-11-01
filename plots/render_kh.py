"""Render Kelvin-Helmholtz run (data_kh/w*.bin): vorticity + magnetic field lines. Domain shown twice along x (periodic)."""
import glob, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
BG = "#11151c"; CM = LinearSegmentedColormap.from_list("sobria", ["#5aa9e6", "#1a2230", "#f0a95f"])
n = int(sys.argv[1]); dt_out = float(sys.argv[2]); b0 = float(sys.argv[3])
k1 = np.fft.fftfreq(n, 1 / n); kx, ky = np.meshgrid(k1, k1, indexing="ij"); k2 = kx**2 + ky**2; k2[0, 0] = 1
y = np.arange(n) * 2 * np.pi / n
def load(f):
    a = np.fromfile(f, np.float32).astype(float).reshape(2, n, n); return a[0].T, a[1].T
def psi(wp, wm):
    z = lambda w: np.real(np.fft.ifft2(-np.fft.fft2(w) / k2 * (k2 > 1e-9)))
    return (z(wp) - z(wm)) / 2 - b0 * y[None, :]
files = sorted(glob.glob("data_kh/w*.bin"))
lim = None
for i, f in enumerate(files):
    wp, wm = load(f); w = (wp + wm) / 2
    if lim is None: lim = 0.9 * np.abs(w).max()
    A = psi(wp, wm); tile = lambda a: np.concatenate([a, a], axis=0)
    fig = plt.figure(figsize=(10, 5), dpi=100, facecolor=BG); ax = fig.add_axes([0, 0, 1, 0.93])
    ax.imshow(tile(w).T, origin="lower", interpolation="bicubic", cmap=CM, vmin=-lim, vmax=lim, extent=[0, 4*np.pi, 0, 2*np.pi])
    ax.contour(np.arange(2*n) * 2*np.pi/n, y, tile(A).T, levels=np.linspace(A.min(), A.max(), 40), colors="#dfe6f0", linewidths=0.5, alpha=0.55)
    ax.axis("off"); fig.text(0.5, 0.955, f"vorticity ω and magnetic field lines,  t = {i*dt_out:.2f}", color="#c9d1d9", ha="center", fontsize=11)
    import os; os.makedirs("png_kh", exist_ok=True); fig.savefig(f"png_kh/f{i:03d}.png", facecolor=BG); plt.close(fig)
