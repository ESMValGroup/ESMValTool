# precipitation11 = mpl_cm.get_cmap("brewer_BrBG_11")
# fix colors from AR6
# https://github.com/IPCC-WG1/colormaps


def rgb(r, g, b):
    return [r / 256, g / 256, b / 256]


historical = rgb(10, 10, 10)
ssp119 = rgb(0, 173, 207)
ssp126 = rgb(23, 60, 102)
ssp245 = rgb(247, 148, 32)
ssp370 = rgb(231, 29, 37)
ssp585 = rgb(149, 27, 30)


experiment_colors = {
    "historical": historical,
    "ssp119": ssp119,
    "ssp126": ssp126,
    "ssp245": ssp245,
    "ssp370": ssp370,
    "ssp585": ssp585,
}


prec_div_5 = [
    rgb(84, 48, 5),
    rgb(200, 148, 79),
    rgb(248, 248, 247),
    rgb(85, 167, 160),
    rgb(0, 60, 48),
]

prec_div_6 = [
    rgb(84, 48, 5),
    rgb(191, 129, 44),
    rgb(229, 209, 180),
    rgb(183, 216, 213),
    rgb(53, 151, 143),
    rgb(0, 60, 48),
]

prec_div_7 = [
    rgb(84, 48, 5),
    rgb(173, 115, 38),
    rgb(216, 182, 135),
    rgb(248, 248, 247),
    rgb(140, 194, 190),
    rgb(44, 135, 127),
    rgb(0, 60, 48),
]

prec_div_8 = [
    rgb(84, 48, 5),
    rgb(160, 105, 33),
    rgb(207, 163, 103),
    rgb(235, 220, 200),
    rgb(202, 225, 223),
    rgb(109, 179, 173),
    rgb(37, 125, 115),
    rgb(0, 60, 48),
]

prec_div_9 = [
    rgb(84, 48, 5),
    rgb(150, 98, 30),
    rgb(200, 148, 79),
    rgb(224, 199, 164),
    rgb(248, 248, 247),
    rgb(167, 208, 204),
    rgb(85, 167, 160),
    rgb(33, 116, 107),
    rgb(0, 60, 48),
]

prec_div_10 = [
    rgb(84, 48, 5),
    rgb(143, 93, 27),
    rgb(195, 137, 60),
    rgb(216, 182, 135),
    rgb(238, 226, 211),
    rgb(212, 230, 229),
    rgb(140, 194, 190),
    rgb(67, 158, 150),
    rgb(29, 110, 100),
    rgb(0, 60, 48),
]

prec_div_11 = [
    rgb(84, 48, 5),
    rgb(137, 88, 25),
    rgb(191, 129, 44),
    rgb(210, 169, 113),
    rgb(229, 209, 180),
    rgb(248, 248, 247),
    rgb(183, 216, 213),
    rgb(118, 183, 178),
    rgb(53, 151, 143),
    rgb(26, 105, 95),
    rgb(0, 60, 48),
]


# prec_11 = mpl_cm.colors.ListedColormap(prec_11, name="pr_11")
# prec_div = mpl_cm.colors.LinearSegmentedColormap.from_list("pr", prec_11)
# prec_seq

# plt.register_cmap('prec_div', prec_div)
# plt.register_cmap('prec_11', prec_11)
