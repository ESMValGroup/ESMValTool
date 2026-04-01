"""Diagnostic script to plot multi panel seasonal cycles for different regions.

NOTE: Its just a script yet, but could be easily adaptable to a diagnostic.

regions: List(List(str)):
    Nested list of regions. One panel for each first level entry.
    Each element can be a list of acronyms to combine (mean) multiple regions.
legend_r: bool
    If True, legend is placed right of the panels. Otherwise below.
"""

from pathlib import Path

import iris
import matplotlib
import yaml
from esmvalcore import preprocessor as pp
from matplotlib import pyplot as plt

SESSION_HIST = Path(
    "/work/bd0854/b309169/output/recipe_lindenlaub25_historical_20250320_162808"
)
REGIONS = [
    ["WAF"],
    ["NES"],
    ["EAS"],
    ["SAS"],
    ["WCE"],
    ["CNA"],
    ["EEU"],
    ["SES"],
    ["WSB"],
    ["MED"],
    ["SAU"],
]
COMBINED = [
    "WCE",
    "CNA",
    "SAS",
    "EAS",
    "EEU",
    "WAF",
    "SES",
    "WSB",
    "MED",
    "NES",
    "SAU",
]
REGIONS.append(COMBINED)
VAR = "tasmin"
# VAR = "pr"
TASMA = True  # calculate tasmin+max/2 instead of tasmin
# tas plot
TSLICE = None
FILL = False
# TSLICE = [0, 120] # 0
# TSLICE = [300, 420] # 0


LEGEND_R = True  # default bottom

FNAME = f"~/seasonal_cycle_{VAR}_w.png"
if TSLICE is not None:
    FNAME = FNAME.replace(".png", f"_{TSLICE[0]}-{TSLICE[1]}.png")
SHARE_Y = False
YLIMS = None  # limit per row or None for auto

cfg = {
    "regions": REGIONS,
    "var": VAR,
    "plot_stdv": FILL,
    "interval": TSLICE,
    "tas_from_min_max": TASMA,
    "cmip6_styles": yaml.safe_load(Path("cmip6.yml").read_text()),
    "share_y": SHARE_Y,
    "ylims": YLIMS,
}

if cfg["var"] == "pr":
    cfg["ylims"] = [(0, 12), (0, 6), (0, 6)]  # limit per row or None for auto
    cfg["share_y"] = True
elif cfg["var"] == "tasmin":
    cfg["share_y"] = True
    cfg["ylims"] = [
        (260, 305),
        (260, 305),
        (260, 305),
    ]  # limit per row or None for auto

hist_metas = yaml.safe_load(
    (SESSION_HIST / f"preproc/pet_historical/{VAR}/metadata.yml").read_text()
)
obs_metas = yaml.safe_load(
    (SESSION_HIST / f"preproc/pet_obs/{VAR}/metadata.yml").read_text()
)
extra_metas = yaml.safe_load(
    (SESSION_HIST / f"preproc/validate_obs/{VAR}/metadata.yml").read_text()
)
obs_metas.update(extra_metas)

metas = obs_metas.copy()
metas.update(hist_metas)


styles = yaml.safe_load(Path("cmip6.yml").read_text())

# REGIONS = [["MED"], ["WCE"], []]
# hist_metas = yaml.safe_load((SESSION_HIST / "preproc/pet_historical/pr/metadata.yml").read_text())
# obs_metas = yaml.safe_load((SESSION_HIST / "preproc/pet_obs/pr/metadata.yml").read_text())

# pr plot


def create_figure_layout(cfg):
    fig = plt.figure()
    axs = []
    nreg = len(cfg["regions"])
    if nreg > 4:
        ncols = 4
        nrows = (nreg // ncols) + (1 if nreg % ncols > 0 else 0)
    else:
        ncols = nreg
        nrows = 1
    if cfg.get("legend_r", False):
        fig.set_size_inches(2.5 * (ncols + 1.5), 2.5 * nrows)
        gs = matplotlib.gridspec.GridSpec(
            nrows,
            ncols + 1,
            width_ratios=[3] * ncols + [1.9],
            wspace=0.07,
            hspace=0.07,
        )
        leg_ax = fig.add_subplot(gs[:, -1])
    else:
        fig.set_size_inches(2.5 * ncols, 2.5 * nrows + 1)
        gs = matplotlib.gridspec.GridSpec(
            nrows + 1,
            ncols,
            height_ratios=[3] * nrows + [1.5],
            wspace=0.07,
            hspace=0.07,
        )
        leg_ax = fig.add_subplot(gs[-1, 0:])
    col_ax = None
    for i in range(nrows):
        for j in range(ncols):
            if i * ncols + j >= nreg:
                break
            ax = fig.add_subplot(gs[i, j])
            if j == 0:
                # TODO: use metadata instead
                if cfg["var"] == "pr":
                    ax.set_ylabel("Precipitation (mm/day)")
                elif cfg["var"] in ["tasmin", "tasmax", "tas"]:
                    ax.set_ylabel("Temperature (K)")
                else:
                    ax.set_ylabel(cfg["var"])
                col_ax = ax
                if cfg["ylims"] is not None:
                    ax.ylim = cfg["ylims"][i]
            elif cfg["share_y"]:
                ax.sharey(col_ax)
            if cfg["share_y"] and j == 0:
                ax.tick_params(axis="y", labelleft=True, labelright=False)
            elif cfg["share_y"] and j == ncols - 1:
                ax.tick_params(axis="y", labelleft=False, labelright=True)
            else:
                ax.tick_params(axis="y", labelleft=False, labelright=False)
            ax.tick_params(
                axis="both", which="both", direction="in", top=True, right=True
            )
            ax.yaxis.grid(True, which="major", linestyle=":", alpha=0.4)
            ax.xaxis.grid(True, which="major", linestyle=":", alpha=0.4)
            axs.append(ax)
            ax.set_xticks(range(1, 13))
            if i == nrows - 1:
                ax.set_xlabel("Month")
                ax.set_xticklabels(
                    [
                        "Jan",
                        "",
                        "Mar",
                        "",
                        "May",
                        "",
                        "Jul",
                        "",
                        "Sep",
                        "",
                        "Nov",
                        "",
                    ]
                )
            else:
                ax.set_xticklabels([])
    return fig, axs, gs, nrows, ncols, leg_ax


def plot_cycle(cube, ax, meta, **kwargs):
    ax.plot(
        cube.coord("month_number").points,
        cube.data,
        label=meta["dataset"],
        **kwargs,
    )


def plot_stdv(cube, std_dev, ax, **kwargs):
    ax.fill_between(
        std_dev.coord("month_number").points,
        cube.data - std_dev.data,
        cube.data + std_dev.data,
        color=kwargs.get("color", "gray"),
        alpha=0.3,
        # label="± 1 std. dev." if "± 1 std. dev." not in ax.get_legend_handles_labels()[1] else None
    )


def guess_bounds(cube):
    if not cube.coord("latitude").has_bounds():
        cube.coord("latitude").guess_bounds()
    if not cube.coord("longitude").has_bounds():
        cube.coord("longitude").guess_bounds()
    return cube


def regional_cycle(cube, operator="mean", regions=None):
    cube = guess_bounds(cube)
    if regions is None:
        regions = ["MED"]
    if len(regions) > 0:
        extracted = pp.extract_shape(
            cube, shapefile="ar6", ids={"Acronym": regions}
        )
    else:
        extracted = cube
    reg_mean = pp.area_statistics(extracted, operator="mean")
    cycle = pp.climate_statistics(
        reg_mean, operator=operator, period="monthly"
    )
    return cycle


def load_cube(fname):
    cube = iris.load_cube(fname)
    if cfg["var"] == "tasmin" and cfg["tas_from_min_max"]:
        max_cube = iris.load_cube(fname.replace("tasmin", "tasmax"))
        cube = (cube + max_cube) / 2
    if cfg["interval"] is not None:
        cube = cube[cfg["interval"][0] : cfg["interval"][1]]
    return cube


def process_region(ax, hist_metas, obs_metas, regions, **kwargs):
    # plot cmip6 models
    for fname, meta in hist_metas.items():
        cube = load_cube(fname)
        cycle = regional_cycle(cube, regions=regions)
        if meta["dataset"] not in styles:
            print(f"WARNING: No style for {meta['dataset']} found!")
        color = styles.get(meta["dataset"], {}).get("color", "black")
        linestyle = styles.get(meta["dataset"], {}).get("dash", "-")
        plot_cycle(
            cycle, ax, meta, linestyle=linestyle, color=color, alpha=0.5
        )
        title = regions[0] if len(regions) == 1 else "Combined"
        ax.text(
            0.06,
            0.92,
            title,
            transform=ax.transAxes,
            va="top",
            ha="left",
            bbox=dict(facecolor="white", alpha=0.6, edgecolor="none"),
        )
    # plot era5
    for fname, meta in obs_metas.items():
        cube = load_cube(fname)
        cycle = regional_cycle(cube, regions=regions)
        color = "red"
        if meta["dataset"].upper() == "CRU":
            color = "black"
        if cfg["plot_stdv"]:
            std_dev = regional_cycle(cube, operator="std_dev")
            plot_stdv(cycle, std_dev, ax, color=color)
        plot_cycle(cycle, ax, meta, color=color, linewidth=1.5)


def add_legend(cfg, leg_ax, axs):
    leg_ax.axis("off")
    hands, labs = axs[0].get_legend_handles_labels()
    # leg_ax.legend(handles=hands, labels=labs, bbox_to_anchor=(1, 1), loc='lower right', ncols=6)
    if cfg.get("legend_r", False):
        bbox = (1.3, 0.5)
        ncols = 1
        loc = "center right"
    else:
        bbox = (0.5, -0.25)
        ncols = 5
        loc = "lower center"
    leg_ax.legend(
        handles=hands,
        labels=labs,
        ncols=ncols,
        loc=loc,
        fancybox=False,
        bbox_to_anchor=bbox,
    )  # Lower position


def main(cfg, metas):
    # create and plot
    fig, axs, gs, nrows, ncols, leg_ax = create_figure_layout(cfg)
    for i, regs in enumerate(cfg["regions"]):
        print("regs:", regs)
        process_region(axs[i], hist_metas, obs_metas, regs)

    add_legend(cfg, leg_ax, axs)
    # save
    fig.tight_layout()
    plt.subplots_adjust(left=0.06, right=0.96, top=0.96, bottom=0.06)
    if cfg.get("tas_from_min_max", False):
        cfg["fname"] = cfg["fname"].replace("tasmin", "tasmid")
    fig.savefig(cfg["fname"], dpi=300)
    print("DONE")


if __name__ == "__main__":
    print("starting")
    main(cfg, metas)
