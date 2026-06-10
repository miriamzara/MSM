
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
import numpy as np


def plot_free_energy_projection(dimred_traj, component_x, component_y, ax = None, vmin=None, vmax=None, axlabel = "PC", colorbarlabel = None, axtitle = None):

    if ax == None:
        fig, ax__ = plt.subplots(1, 1, figsize=(5, 4))
    else:
        ax__ = ax

    # Scatter plot
    ax__.scatter(
        dimred_traj[:, component_x],
        dimred_traj[:, component_y],
        s=0.01,
        alpha=0.3
    )

    ax__.xaxis.tick_top()
    ax__.xaxis.set_label_position("top")
    ax__.set_xlabel(f"{axlabel}{component_x}", labelpad=1)
    ax__.set_ylabel(f"{axlabel}{component_y}")


    # 2D histogram (Free energy)
    H, xedges, yedges = np.histogram2d(
        dimred_traj[:, component_x],
        dimred_traj[:, component_y],
        bins=100,
        density=True
    )
    mask = H > 0
    F = np.full_like(H, np.nan)
    F[mask] = -np.log(H[mask])
    F -= np.nanmin(F)

    im = ax__.imshow(
        F.T,
        origin="lower",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        aspect="auto",
        vmin=vmin,
        vmax=vmax
    )
    #if colorbarlabel is not None:
    #    plt.colorbar(im, ax=ax__, label=r"$F/k_B T$")
    #else:
    #    plt.colorbar(im, ax=ax__)
    ax__.set_xlabel(f"{axlabel}{component_x}")
    ax__.set_ylabel(f"{axlabel}{component_y}")
    if axtitle is not None:
        ax__.set_title(f"{axtitle}")

    return ax__, im


###########################################
#
############################################

def add_1d_free_energy_axis(
    ax,
    x,
    bins=100,
    width_fraction=0.22,
    pad_fraction=0.02,
    color="C0",
):
    """
    Add a right-side 1D free-energy profile to an existing time-series axis.
    The new axis shares y with ax.
    """

    fig = ax.figure
    pos = ax.get_position()

    profile_width = pos.width * width_fraction
    pad = pos.width * pad_fraction
    trace_width = pos.width - profile_width - pad

    # shrink original axis
    ax.set_position([pos.x0, pos.y0, trace_width, pos.height])

    # create new profile axis
    ax_fe = fig.add_axes([
        pos.x0 + trace_width + pad,
        pos.y0,
        profile_width,
        pos.height
    ], sharey=ax)

    hist, edges = np.histogram(x, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])

    F = -np.log(hist + 1e-12)
    F -= np.nanmin(F)

    ax_fe.plot(F, centers, lw=1.5, color=color)

    ax_fe.tick_params(axis="y", labelleft=False)
    ax_fe.grid(False)

    return ax_fe