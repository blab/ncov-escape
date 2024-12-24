import colorsys

import numpy as np
from pango_aliasor.aliasor import Aliasor

clade_definitions = [
    {
        "clade": "24B",
        "display_name": "24B (JN.1.11.1)",
        "defining_lineage": "JN.1.11.1",
        "color": "#DC2F24",
    },
    {
        "clade": "24A",
        "display_name": "24A (JN.1)",
        "defining_lineage": "JN.1",
        "color": "#E4632E",
    },
    {
        "clade": "23I",
        "display_name": "23I (BA.2.86)",
        "defining_lineage": "BA.2.86",
        "color": "#E69136",
    },
    {
        "clade": "23H",
        "display_name": "23H (HK.3)",
        "defining_lineage": "HK.3",
        "color": "#D9AD3D",
    },
    {
        "clade": "23G",
        "display_name": "23G (XBB.1.5.70)",
        "defining_lineage": "XBB.1.5.70",
        "color": "#C1BA47",
    },
    {
        "clade": "23F",
        "display_name": "23F (EG.5.1)",
        "defining_lineage": "EG.5.1",
        "color": "#A2BE57",
    },
    {
        "clade": "23E",
        "display_name": "23E (XBB.2.3)",
        "defining_lineage": "XBB.2.3",
        "color": "#83BA70",
    },
    {
        "clade": "23D",
        "display_name": "23D (XBB.1.9)",
        "defining_lineage": "XBB.1.9",
        "color": "#69B091",
    },
    {
        "clade": "23C",
        "display_name": "23C (CH.1.1)",
        "defining_lineage": "CH.1.1",
        "color": "#549DB2",
    },
    {
        "clade": "23B",
        "display_name": "23B (XBB.1.16)",
        "defining_lineage": "XBB.1.16",
        "color": "#4580CA",
    },
    {
        "clade": "23A",
        "display_name": "23A (XBB.1.5)",
        "defining_lineage": "XBB.1.5",
        "color": "#462EB9",
    },
    {
        "clade": "22F",
        "display_name": "22F (XBB)",
        "defining_lineage": "XBB",
        "color": "#3E58CF",
    },
    {
        "clade": "22E",
        "display_name": "22E (BQ.1)",
        "defining_lineage": "BQ.1",
        "color": "#777777",
    },
    {
        "clade": "22D",
        "display_name": "22D (BA.2.75)",
        "defining_lineage": "BA.2.75",
        "color": "#777777",
    },
    {
        "clade": "22B",
        "display_name": "22B (BA.5)",
        "defining_lineage": "BA.5",
        "color": "#777777",
    },
    {
        "clade": "other",
        "display_name": "other",
        "defining_lineage": False,
        "color": "#777777",
    },
]


DEFAULT_CLADE_COLOR = "#777777"


def order_lineages(lineages, aliasor):
    """
    Order input lineages by using their full uncompressed lineage & converting to a sortable form
    e.g. BA.5  -> B.1.1.529.5  -> '  B/001/001/529/005'
         BA.16 -> B.1.1.529.16 -> '  B/001/001/529/016'
         so BA.5 is before BA.16
    """

    def _lineage_sortable(lineage):
        if lineage == "other":
            return "ZZZ"
        lin_full = aliasor.uncompress(lineage)
        return "/".join(
            [
                (f"{x:>3}" if i == 0 else f"{int(x):03}")
                for i, x in enumerate(lin_full.split("."))
            ]
        )

    return sorted(lineages, key=_lineage_sortable)


def order_lineages_within_clade(lineages, aliasor, clade_map):
    """
    Order input lineages by their clades first and then by using their full uncompressed lineage,
    converting to a sortable form within each clade.

    Parameters:
    lineages (list of str): List of lineage names to be sorted.
    aliasor (object): An object with a method `uncompress` that converts compressed lineage names to full lineage names.
    clades (list of str): List of clade assignments corresponding to each lineage.

    Returns:
    list of str: Sorted list of lineage names.

    Example:
    >>> class Aliasor:
    ...     def uncompress(self, lineage):
    ...         return {"BA.5": "B.1.1.529.5", "BA.16": "B.1.1.529.16", "BA.2": "B.1.1.529.2"}[lineage]
    >>> aliasor = Aliasor()
    >>> lineages = ["BA.5", "BA.16", "BA.2"]
    >>> clades = ["Clade1", "Clade1", "Clade2"]
    >>> order_lineages(lineages, aliasor, clades)
    ['BA.5', 'BA.16', 'BA.2']
    """

    def _lineage_sortable(lineage):
        if lineage == "other":
            return "ZZZ"
        lin_full = aliasor.uncompress(lineage)
        return "/".join(
            [
                (f"{x:>3}" if i == 0 else f"{int(x):03}")
                for i, x in enumerate(lin_full.split("."))
            ]
        )

    # Sort the lineages first by clade and then by the sortable form of the lineage
    return sorted(
        lineages, key=lambda lineage: (clade_map[lineage], _lineage_sortable(lineage))
    )


def lineage_to_clade(lineage, aliasor, fallback, clade_definitions):
    lineage_full = aliasor.uncompress(lineage)
    for clade_data in clade_definitions:
        if clade_data["clade"] == "other":
            continue
        comparison_lineage = aliasor.uncompress(clade_data["defining_lineage"])
        if lineage_full == comparison_lineage or lineage_full.startswith(
            comparison_lineage + "."
        ):
            return clade_data["clade"]
    return fallback


def clade_colors(variants, clade_definitions):
    colors = {c["clade"]: c["color"] for c in clade_definitions}
    missing = set()
    defs = []
    for v in variants:
        try:
            defs.append([v, colors[v]])
        except KeyError:
            if v != "other":
                missing.add(v)
                defs.append([v, DEFAULT_CLADE_COLOR])

    if len(missing) > 0:
        print(
            f"Missing definitions for the following clades: {', '.join(missing)}.",
            f"They have been assigned the default color {DEFAULT_CLADE_COLOR!r}",
        )

    return defs


def clade_display_names(variants, clade_definitions):
    display_names = {c["clade"]: c["display_name"] for c in clade_definitions}
    return [
        [name, display_names[name] if name in display_names else name]
        for name in variants
    ]


def colour_range(anchor, n):
    """
    Create a range of `n` colours centred around the provided `anchor`.
    This currently involves simple manipulations in HLS space, but
    the outputs aren't going to be as good as they could be if we did it in
    a perceptually uniform space (e.g lab space). For the purposes of this viz
    I don't think it's a dealbreaker, and in our current setup it's hard to use
    python libraries which aren't already available in our various runtimes.
    """
    anchor_rgb = tuple(int(anchor.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
    anchor_hls = colorsys.rgb_to_hls(*anchor_rgb)
    hrange = np.linspace(anchor_hls[0] * 0.85, anchor_hls[0] * 1.25, n)
    lrange = np.linspace(anchor_hls[1] * 1.2, anchor_hls[1], n)
    srange = np.linspace(anchor_hls[2] * 0.7, anchor_hls[2] * 1.1, n)
    rgb_range = [colorsys.hls_to_rgb(*hls) for hls in zip(hrange, lrange, srange)]

    def clamp(x):
        return int(max(0, min(x, 255)))

    return [
        f"#{clamp(rgb[0]):02x}{clamp(rgb[1]):02x}{clamp(rgb[2]):02x}"
        for rgb in rgb_range
    ]


def colourise(lineages, aliasor, clade_definitions):
    """
    Produces an array of arrays associating observed lineages with a colour hex. Example output:
        [
            ['XBB', '#ffffff'],
            ...
        ]
    """
    clades = {
        lineage: lineage_to_clade(lineage, aliasor, "other", clade_definitions)
        for lineage in lineages
    }

    colours = []

    for clade in list(set(clades.values())):
        matching_lineages = [
            l for l in lineages if clades[l] == clade
        ]  # will be ordered
        print(f"{clade:<10}n={len(matching_lineages)} lineages")
        color_hex = [x["color"] for x in clade_definitions if x["clade"] == clade][0]
        for pair in zip(
            matching_lineages, colour_range(color_hex, len(matching_lineages))
        ):
            colours.append(pair)
    return colours
