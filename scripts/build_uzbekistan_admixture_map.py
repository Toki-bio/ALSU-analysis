#!/usr/bin/env python3
"""Build an Uzbekistan-region ADMIXTURE map from verified local inputs."""

from __future__ import annotations

import csv
import json
import math
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch, Wedge


ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "GWAS от 27.08 - For_Plink.csv"
FAM_PATH = ROOT / "data" / "step11_global_for_admixture.fam"
Q_PATH = ROOT / "data" / "step11_global_for_admixture.5.Q"
CACHE_PATH = ROOT / "data" / "natural_earth_uzb_admin1.geojson"
SUMMARY_PATH = ROOT / "data" / "alsu_uzbekistan_k5_birthplace_summary.tsv"
OUT_PATH = ROOT / "images" / "alsu_uzbekistan_k5_birthplace_map.png"

NATURAL_EARTH_ADMIN1_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
    "ne_10m_admin_1_states_provinces.geojson"
)

# Raw Step 11 Spring 2026 K=5 order, verified against 1000G superpopulation means:
# Q1=East Asian, Q2=AMR-like/Other, Q3=South Asian, Q4=European/West Eurasian, Q5=African.
PLOT_COMPONENTS = [
    ("European / West Eurasian", 3, "#3f7fc1"),  # strong blue
    ("East Asian",               0, "#f9c03f"),  # golden-yellow
    ("South Asian",              2, "#57a059"),  # green
    ("Others",                   1, "#9e9e9e"),  # neutral grey
    ("African",                  4, "#e53935"),  # red
]

REGION_TO_ISO = {
    "Tashkent": "UZ-TK",
    "Republic of Karakalpakstan": "UZ-QR",
    "Andijan Region": "UZ-AN",
    "Bukhara Region": "UZ-BU",
    "Jizzakh Region": "UZ-JI",
    "Kashkadarya Region": "UZ-QA",
    "Navoiy Region": "UZ-NW",
    "Namangan Region": "UZ-NG",
    "Samarkand Region": "UZ-SA",
    "Surxondaryo Region": "UZ-SU",
    "Syrdarya Region": "UZ-SI",
    "Tashkent Region": "UZ-TO",
    "Fergana Region": "UZ-FA",
    "Khorezm Region": "UZ-XO",
}

LABELS = {
    "Republic of Karakalpakstan": "Karakalpakstan",
    "Kashkadarya Region": "Kashkadarya",
    "Surxondaryo Region": "Surxondaryo",
    "Tashkent Region": "Tashkent\nregion",
    "Fergana Region": "Fergana",
    "Andijan Region": "Andijan",
    "Syrdarya Region": "Syrdarya",
    "Samarkand Region": "Samarkand",
    "Bukhara Region": "Bukhara",
    "Jizzakh Region": "Jizzakh",
    "Navoiy Region": "Navoiy",
    "Namangan Region": "Namangan",
    "Khorezm Region": "Khorezm",
}

# Display offsets only separate overlapping pies; leader lines retain the source coordinates.
DISPLAY_OFFSETS = {
    "Tashkent": (-0.75, 0.40),
    "Tashkent Region": (0.55, 0.90),
    "Andijan Region": (0.35, 0.08),
    "Namangan Region": (0.75, 0.45),
    "Fergana Region": (0.10, -0.38),
    "Syrdarya Region": (-0.25, -0.28),
}


def parse_codebook(column_name: str) -> dict[str, str]:
    codebook: dict[str, str] = {}
    for line in column_name.splitlines()[1:]:
        match = re.match(r"^\s*(\d+(?:\.\d+)?)\s+(.+?)\s*$", line)
        if match:
            codebook[match.group(1)] = match.group(2)
    return codebook


def load_birthplace_codes() -> tuple[dict[str, str], dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = list(reader)

    sample_i = header.index("Sample Number")
    birth_i = next(i for i, name in enumerate(header) if name.startswith("Place of birth"))
    codebook = parse_codebook(header[birth_i])
    sample_to_code = {row[sample_i].strip(): row[birth_i].strip() for row in rows if len(row) > birth_i}
    return sample_to_code, codebook


def load_step11_q() -> list[tuple[str, list[float]]]:
    fam_ids = [line.split()[1] for line in FAM_PATH.read_text(encoding="utf-8").splitlines() if line.strip()]
    q_rows = [list(map(float, line.split())) for line in Q_PATH.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(fam_ids) != len(q_rows):
        raise ValueError(f"FAM/Q row mismatch: {len(fam_ids)} IDs vs {len(q_rows)} Q rows")
    return list(zip(fam_ids, q_rows))


def aggregate_by_birthplace() -> tuple[list[dict[str, object]], dict[str, int]]:
    sample_to_code, codebook = load_birthplace_codes()
    grouped: dict[str, dict[str, object]] = defaultdict(lambda: {"n": 0, "sum": [0.0] * 5})
    counters = {"missing_birthplace": 0, "non_uzbekistan_birthplace": 0, "unknown_region_code": 0}

    for sample_id, q_values in load_step11_q():
        code = sample_to_code.get(sample_id)
        if not code or code == "NA":
            counters["missing_birthplace"] += 1
            continue
        if not code.startswith("1."):
            counters["non_uzbekistan_birthplace"] += 1
            continue
        region = codebook.get(code)
        if region not in REGION_TO_ISO:
            counters["unknown_region_code"] += 1
            continue

        record = grouped[region]
        record["n"] = int(record["n"]) + 1
        record["sum"] = [a + b for a, b in zip(record["sum"], q_values)]

    rows = []
    for region, record in grouped.items():
        n = int(record["n"])
        means = [value / n for value in record["sum"]]
        rows.append(
            {
                "region": region,
                "iso_3166_2": REGION_TO_ISO[region],
                "n": n,
                "q1_east_asian": means[0],
                "q2_amr_like": means[1],
                "q3_south_asian": means[2],
                "q4_european_west_eurasian": means[3],
                "q5_african": means[4],
            }
        )
    rows.sort(key=lambda row: (-int(row["n"]), str(row["region"])))
    return rows, counters


def load_uzbekistan_admin1() -> dict[str, object]:
    if CACHE_PATH.exists():
        return json.loads(CACHE_PATH.read_text(encoding="utf-8"))

    with urllib.request.urlopen(NATURAL_EARTH_ADMIN1_URL, timeout=45) as response:
        all_admin1 = json.loads(response.read().decode("utf-8"))
    uzbekistan = {
        "type": "FeatureCollection",
        "source": NATURAL_EARTH_ADMIN1_URL,
        "features": [
            feature
            for feature in all_admin1["features"]
            if feature.get("properties", {}).get("adm0_a3") == "UZB"
        ],
    }
    CACHE_PATH.write_text(json.dumps(uzbekistan, ensure_ascii=False, indent=2), encoding="utf-8")
    return uzbekistan


def geometry_rings(geometry: dict[str, object]) -> list[list[list[float]]]:
    if geometry["type"] == "Polygon":
        return [geometry["coordinates"][0]]
    if geometry["type"] == "MultiPolygon":
        return [polygon[0] for polygon in geometry["coordinates"]]
    return []


def draw_pie(ax, ratios: list[float], x: float, y: float, radius: float) -> None:
    theta = 90.0
    for ratio, (_, _, color) in zip(ratios, PLOT_COMPONENTS):
        if ratio <= 0:
            continue
        next_theta = theta - 360.0 * ratio
        ax.add_patch(
            Wedge(
                (x, y),
                radius,
                next_theta,
                theta,
                facecolor=color,
                edgecolor="white",
                linewidth=0.55,
                zorder=6,
            )
        )
        theta = next_theta
    ax.add_patch(Circle((x, y), radius, facecolor="none", edgecolor="#263238", linewidth=0.55, zorder=7))


def write_summary(rows: list[dict[str, object]]) -> None:
    fields = [
        "region",
        "iso_3166_2",
        "n",
        "q1_east_asian",
        "q2_amr_like",
        "q3_south_asian",
        "q4_european_west_eurasian",
        "q5_african",
    ]
    with SUMMARY_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (f"{row[key]:.6f}" if isinstance(row[key], float) else row[key])
                    for key in fields
                }
            )


def build_map(rows: list[dict[str, object]], counters: dict[str, int]) -> None:
    admin1 = load_uzbekistan_admin1()
    features = admin1["features"]
    feature_by_iso = {feature["properties"]["iso_3166_2"]: feature for feature in features}

    fig, ax = plt.subplots(figsize=(13.5, 8.2), dpi=220)
    ax.set_facecolor("#f8faf7")

    for feature in features:
        for ring in geometry_rings(feature["geometry"]):
            xs = [point[0] for point in ring]
            ys = [point[1] for point in ring]
            ax.fill(xs, ys, facecolor="#eef3ef", edgecolor="#aeb8b2", linewidth=0.75, zorder=1)

    max_n = max(int(row["n"]) for row in rows)
    plotted_n = 0
    for row in rows:
        region = str(row["region"])
        feature = feature_by_iso[str(row["iso_3166_2"])]
        props = feature["properties"]
        anchor_x = float(props["longitude"])
        anchor_y = float(props["latitude"])
        dx, dy = DISPLAY_OFFSETS.get(region, (0.0, 0.0))
        x = anchor_x + dx
        y = anchor_y + dy
        if dx or dy:
            ax.plot([anchor_x, x], [anchor_y, y], color="#607d6f", linewidth=0.65, zorder=4)

        n = int(row["n"])
        plotted_n += n
        radius = 0.16 + 0.34 * math.sqrt(n / max_n)
        raw_values = [
            float(row["q1_east_asian"]),
            float(row["q2_amr_like"]),
            float(row["q3_south_asian"]),
            float(row["q4_european_west_eurasian"]),
            float(row["q5_african"]),
        ]
        ratios = [raw_values[index] for _, index, _ in PLOT_COMPONENTS]
        draw_pie(ax, ratios, x, y, radius)

        label = LABELS.get(region, region.replace(" Region", ""))
        ax.text(
            x,
            y - radius - 0.09,
            f"{label}\nN={n}",
            ha="center",
            va="top",
            fontsize=6.8,
            color="#263238",
            linespacing=1.05,
            zorder=8,
        )

    ax.set_xlim(55.0, 74.2)
    ax.set_ylim(37.0, 46.3)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.text(
        0.055,
        0.955,
        "ALSU K=5 ADMIXTURE by Uzbekistan Birthplace Region",
        fontsize=18,
        fontweight="bold",
        color="#1f2a2e",
    )
    fig.text(
        0.055,
        0.925,
        (
            f"Mean Step 11 global K=5 proportions; pies sized by sqrt(N). "
            f"Mapped samples: {plotted_n}; missing/non-Uzbekistan birthplace: "
            f"{counters['missing_birthplace'] + counters['non_uzbekistan_birthplace']}"
        ),
        fontsize=9.2,
        color="#536068",
    )

    legend_handles = [Patch(facecolor=color, edgecolor="none", label=label) for label, _, color in PLOT_COMPONENTS]
    ax.legend(
        handles=legend_handles,
        loc="lower left",
        bbox_to_anchor=(0.02, 0.025),
        frameon=True,
        framealpha=0.95,
        facecolor="white",
        edgecolor="#d0d7d2",
        fontsize=8.4,
        title="K=5 components",
        title_fontsize=9.2,
    )
    ax.text(
        73.8,
        37.15,
        "Boundary source: Natural Earth 1:10m admin-1\nGeography source: GWAS от 27.08 - For_Plink.csv",
        ha="right",
        va="bottom",
        fontsize=6.7,
        color="#607d6f",
    )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    rows, counters = aggregate_by_birthplace()
    if not rows:
        raise SystemExit("No Uzbekistan birthplace rows could be joined to Step 11 ADMIXTURE samples.")
    write_summary(rows)
    build_map(rows, counters)
    print(f"Wrote {SUMMARY_PATH.relative_to(ROOT)}")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}")
    print(f"Mapped samples: {sum(int(row['n']) for row in rows)}")
    print(f"Counters: {counters}")


if __name__ == "__main__":
    main()