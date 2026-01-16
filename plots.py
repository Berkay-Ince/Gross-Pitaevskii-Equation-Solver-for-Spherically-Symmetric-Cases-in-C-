# plot_sweeps_cpp_npz.py
#
# C++ cnpy .npz compatible + y-axis chosen ONLY from non-diverged curves:
#   ylim = (0, 1.05*y_max_nondiv)
#
# Divergence detection uses a robust scale (quantile) > div_thresh.

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


# -----------------------------
# Utilities
# -----------------------------
def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def fmt_float(x, nd=3):
    return f"{x:.{nd}g}"


def decode_npz_string(x):
    a = np.asarray(x)
    if a.dtype == np.uint8:
        return bytes(a.tolist()).decode("utf-8")
    return str(a)


# -----------------------------
# Robust helpers
# -----------------------------
def robust_curve_scale(dens_plot, q=0.95):
    """
    Robust 'typical max' for one curve using a high quantile (ignores NaNs/Infs).
    dens_plot includes an extra first point; ignore it.
    """
    y = np.asarray(dens_plot[1:], dtype=float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return 0.0
    return float(np.quantile(y, q))


def is_diverged_curve(dens_plot, div_thresh=50.0, q=0.95):
    """
    Mark as diverged if robust typical scale is above threshold.
    """
    return robust_curve_scale(dens_plot, q=q) > div_thresh


def pick_ylim_from_nondiverged(dens_plots, diverged_flags, pad=1.05, q=0.995):
    """
    Compute y-limit using ONLY non-diverged curves:
        ylim = (0, pad * y_max_nondiv)
    y_max_nondiv is taken as robust quantile of dens values across each curve
    then max across curves.

    If all curves diverge, fallback to all curves (still robust).
    """
    dens_ok = [d for d, div in zip(dens_plots, diverged_flags) if not div]

    use = dens_ok if len(dens_ok) > 0 else dens_plots
    if len(use) == 0:
        return (0.0, 1.0)

    # y_max from robust quantile to avoid single-point spikes
    y_max = max(robust_curve_scale(d, q=q) for d in use)

    if not np.isfinite(y_max) or y_max <= 0.0:
        y_max = 1.0

    return (0.0, pad * y_max)


# -----------------------------
# Load .npz results (C++ format)
# -----------------------------
def load_all_npz_results(data_dir):
    paths = sorted(glob.glob(os.path.join(data_dir, "*.npz")))
    if not paths:
        raise FileNotFoundError(f"No .npz files found in: {data_dir}")

    results = []
    for p in paths:
        d = np.load(p, allow_pickle=False)

        r = d["r"]
        psi = d["psi_final"]

        method = decode_npz_string(d["method"])
        dt_input = decode_npz_string(d["dt_input"])

        g = float(np.asarray(d["g"]).reshape(-1)[0])
        R = float(np.asarray(d["R"]).reshape(-1)[0])
        dt_value = float(np.asarray(d["dt_value"]).reshape(-1)[0])

        energy = float(np.asarray(d["energy"]).reshape(-1)[0]) if "energy" in d else None
        steps = int(np.asarray(d["steps"]).reshape(-1)[0]) if "steps" in d else None

        results.append(
            {
                "path": p,
                "method": method,
                "g": g,
                "R": R,
                "dt": dt_value,
                "dt_input": dt_input,
                "energy": energy,
                "steps": steps,
                "r": r,
                "psi": psi,
            }
        )
    return results


# -----------------------------
# Density builder (psi = r*phi) -> density ~ |psi|^2 / r^2
# -----------------------------
def build_density_arrays(r, psi, k0=3):
    rr = r
    dens = (np.abs(psi) ** 2) / (rr ** 2 + 1e-30)

    # regularize near r=0
    if dens.size > k0:
        plateau = dens[1 : k0 + 1].mean()
        dens[0 : k0 + 1] = plateau

    dens0 = dens[0]

    r_plot = np.empty(rr.size + 1, dtype=float)
    dens_plot = np.empty(dens.size + 1, dtype=float)

    r_plot[0] = 0.0
    r_plot[1:] = rr
    dens_plot[0] = dens0
    dens_plot[1:] = dens

    return r_plot, dens_plot


# -----------------------------
# Plot style common
# -----------------------------
def _style_axes(ax, xlim):
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_xlabel("r", fontsize=18)
    ax.set_ylabel(r"$|\psi(r)|^2/r^2$", fontsize=18)
    ax.set_xlim(*xlim)
    ax.grid(True, linestyle="-", alpha=0.25)


def _save(fig, outdir, tag):
    ensure_dir(outdir)
    fig.tight_layout()
    plt.savefig(os.path.join(outdir, f"{tag}.pdf"), dpi=200, bbox_inches="tight")
    plt.savefig(os.path.join(outdir, f"{tag}.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# 1) Compare methods for each (g,R,dt_input)
# -----------------------------
def plot_compare_methods(all_results, outdir, xlim=(0.05, 3.0),
                         ylim=None, div_thresh=50.0):
    ensure_dir(outdir)

    groups = {}
    for res in all_results:
        key = (res["g"], res["R"], res["dt_input"])
        groups.setdefault(key, []).append(res)

    for (g, R, dt_input), subset in groups.items():
        subset.sort(key=lambda d: d["method"])

        # unique per method
        unique = {}
        for res in subset:
            unique[res["method"]] = res
        subset = [unique[m] for m in sorted(unique)]

        fig, ax = plt.subplots(figsize=(4, 6))

        curves = []
        dens_plots = []
        diverged_flags = []

        for res in subset:
            r_plot, dens_plot = build_density_arrays(res["r"], res["psi"])
            div = is_diverged_curve(dens_plot, div_thresh=div_thresh, q=0.95)

            curves.append((res["method"], r_plot, dens_plot, div))
            dens_plots.append(dens_plot)
            diverged_flags.append(div)

        # y-limit: ONLY from non-diverged curves (unless user provides ylim)
        if ylim is not None:
            ax.set_ylim(*ylim)
        else:
            ax.set_ylim(*pick_ylim_from_nondiverged(dens_plots, diverged_flags, pad=1.05, q=0.995))

        # plot all curves; shade only non-diverged
        for method, r_plot, dens_plot, div in curves:
            label = method + (" (div)" if div else "")
            (line,) = ax.plot(r_plot, dens_plot, linewidth=2, label=label)

            if not div:
                baseline = dens_plot[1:].min() - 0.01
                ax.fill_between(r_plot[1:], dens_plot[1:], baseline, alpha=0.15, color=line.get_color())

        _style_axes(ax, xlim)

        ax.legend(frameon=True, framealpha=0.85, facecolor="white", edgecolor="0.8",
                  fancybox=True, borderpad=0.5)

        dt_tag = dt_input.replace(" ", "").replace("*", "_").replace("/", "div").replace(".", "p")
        tag = f"compare_methods__g{fmt_float(g)}__R{fmt_float(R)}__dt{dt_tag}"
        _save(fig, outdir, tag)


# -----------------------------
# 2) Fixed (R,dt): vary g per method
# -----------------------------
def plot_fixed_R_dt_vary_g_per_method(all_results, outdir, xlim=(0.05, 3.0),
                                      ylim=None, div_thresh=50.0):
    ensure_dir(outdir)

    groups = {}
    for res in all_results:
        key = (res["method"], res["R"], res["dt_input"])
        groups.setdefault(key, []).append(res)

    for (method, R, dt_input), subset in groups.items():
        subset.sort(key=lambda d: d["g"])

        fig, ax = plt.subplots(figsize=(4, 6))

        curves = []
        dens_plots = []
        diverged_flags = []

        for res in subset:
            r_plot, dens_plot = build_density_arrays(res["r"], res["psi"])
            div = is_diverged_curve(dens_plot, div_thresh=div_thresh, q=0.95)

            label = f"g={res['g']:.3g}" + (" (div)" if div else "")
            curves.append((label, r_plot, dens_plot, div))
            dens_plots.append(dens_plot)
            diverged_flags.append(div)

        if ylim is not None:
            ax.set_ylim(*ylim)
        else:
            ax.set_ylim(*pick_ylim_from_nondiverged(dens_plots, diverged_flags, pad=1.05, q=0.995))

        for label, r_plot, dens_plot, div in curves:
            (line,) = ax.plot(r_plot, dens_plot, linewidth=2, label=label)
            if not div:
                baseline = dens_plot[1:].min() - 0.01
                ax.fill_between(r_plot[1:], dens_plot[1:], baseline, alpha=0.15, color=line.get_color())

        _style_axes(ax, xlim)

        ax.legend(frameon=True, framealpha=0.85, facecolor="white", edgecolor="0.8",
                  fancybox=True, borderpad=0.5)

        dt_tag = dt_input.replace(" ", "").replace("*", "_").replace("/", "div").replace(".", "p")
        tag = f"vary_g__method{method}__R{fmt_float(R)}__dt{dt_tag}"
        _save(fig, outdir, tag)


# -----------------------------
# 3) Fixed (g,dt): vary R per method
# -----------------------------
def plot_fixed_g_dt_vary_R_per_method(all_results, outdir, xlim=(0.05, 3.0),
                                      ylim=None, div_thresh=50.0):
    ensure_dir(outdir)

    groups = {}
    for res in all_results:
        key = (res["method"], res["g"], res["dt_input"])
        groups.setdefault(key, []).append(res)

    for (method, g, dt_input), subset in groups.items():
        subset.sort(key=lambda d: d["R"])

        fig, ax = plt.subplots(figsize=(4, 6))

        curves = []
        dens_plots = []
        diverged_flags = []

        for res in subset:
            r_plot, dens_plot = build_density_arrays(res["r"], res["psi"])
            div = is_diverged_curve(dens_plot, div_thresh=div_thresh, q=0.95)

            label = f"R={res['R']:.3g}" + (" (div)" if div else "")
            curves.append((label, r_plot, dens_plot, div))
            dens_plots.append(dens_plot)
            diverged_flags.append(div)

        if ylim is not None:
            ax.set_ylim(*ylim)
        else:
            ax.set_ylim(*pick_ylim_from_nondiverged(dens_plots, diverged_flags, pad=1.05, q=0.995))

        for label, r_plot, dens_plot, div in curves:
            (line,) = ax.plot(r_plot, dens_plot, linewidth=2, label=label)
            if not div:
                baseline = dens_plot[1:].min() - 0.01
                ax.fill_between(r_plot[1:], dens_plot[1:], baseline, alpha=0.15, color=line.get_color())

        _style_axes(ax, xlim)

        ax.legend(frameon=True, framealpha=0.85, facecolor="white", edgecolor="0.8",
                  fancybox=True, borderpad=0.5)

        dt_tag = dt_input.replace(" ", "").replace("*", "_").replace("/", "div").replace(".", "p")
        tag = f"vary_R__method{method}__g{fmt_float(g)}__dt{dt_tag}"
        _save(fig, outdir, tag)


# -----------------------------
# CSV summary
# -----------------------------
def save_energies_csv(all_results, outpath):
    ensure_dir(os.path.dirname(outpath))
    with open(outpath, "w") as f:
        f.write("method,g,R,dt_input,dt_value,energy,steps,path\n")
        for res in all_results:
            f.write(
                f"{res['method']},"
                f"{res['g']},"
                f"{res['R']},"
                f"{res['dt_input']},"
                f"{res['dt']},"
                f"{res['energy']},"
                f"{res['steps']},"
                f"{res['path']}\n"
            )


# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":
    DATA_DIR = "gpe_runs"   # your C++ out_dir
    PLOT_DIR = "gpe_plots"

    all_results = load_all_npz_results(DATA_DIR)

    plot_compare_methods(all_results, outdir=os.path.join(PLOT_DIR, "compare_methods"), div_thresh=50.0)
    plot_fixed_R_dt_vary_g_per_method(all_results, outdir=os.path.join(PLOT_DIR, "vary_g_fixed_R_dt"), div_thresh=50.0)
    plot_fixed_g_dt_vary_R_per_method(all_results, outdir=os.path.join(PLOT_DIR, "vary_R_fixed_g_dt"), div_thresh=50.0)

    save_energies_csv(all_results, outpath=os.path.join(PLOT_DIR, "energies", "energy_summary.csv"))
