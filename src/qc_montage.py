"""Tile per-participant registration QC images into contact sheets for review.
Usage: python -I src/qc_montage.py  -> reports/figures/qc_montage_XX.png
(Participant IDs appear on these sheets: they are QC artefacts of public data.)"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

files = sorted(Path("derivatives/oxtfc").glob("sub-*/sub-*_qc.png"))
per = 12
for k in range(0, len(files), per):
    chunk = files[k:k + per]
    fig, axes = plt.subplots(4, 3, figsize=(15, 13), facecolor="black")
    for ax in axes.ravel():
        ax.axis("off")
    for ax, f in zip(axes.ravel(), chunk):
        ax.imshow(mpimg.imread(f))
    fig.tight_layout(pad=0.2)
    out = Path(f"derivatives/oxtfc/qc_montage_{k // per + 1:02d}.png")
    fig.savefig(out, dpi=60, facecolor="black")
    plt.close(fig)
    print(out)
