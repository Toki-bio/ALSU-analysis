#!/usr/bin/env python3
"""
Generate the Step 2 sample-pair PI_HAT histogram JSON.

This uses the Spring 2026 ConvSK_mind20 sample set. PI_HAT is a sample-to-sample
relatedness estimate; markers are only the variants used to estimate that value.
"""

import subprocess
import json
import sys

BINS = [
    ("0.00-0.05", "Very low PI_HAT"),
    ("0.05-0.10", "Low relatedness"),
    ("0.10-0.15", "Elevated relatedness"),
    ("0.15-0.25", "Approximately third- to second-degree range"),
    ("0.25-0.50", "Approximately second- to first-degree range"),
    ("0.50-0.75", "First-degree / full-sibling range"),
    ("0.75-0.90", "Very close relatives"),
    ("0.90-0.98", "Near-identical candidate range"),
    ("0.98-1.00", "Duplicate / technical replicate threshold"),
]

REMOTE_SCRIPT = r"""
set -euo pipefail
cd /staging/ALSU-analysis/spring2026

if [ ! -s ConvSK_mind20_full.genome ]; then
  plink --bfile ConvSK_mind20 --genome --out ConvSK_mind20_full
fi

sample_count=$(wc -l < ConvSK_mind20.fam)
marker_count_total=$(wc -l < ConvSK_mind20.bim)
expected_sample_pairs=$((sample_count * (sample_count - 1) / 2))
genome_line_count=$(wc -l < ConvSK_mind20_full.genome)
actual_sample_pairs=$((genome_line_count - 1))
excluded_non_autosomal=$(sed -n 's/^Excluding \([0-9][0-9]*\) variants.*/\1/p' ConvSK_mind20_full.log | tail -1)
if [ -z "$excluded_non_autosomal" ]; then
  excluded_non_autosomal=0
fi
marker_count_ibd=$((marker_count_total - excluded_non_autosomal))

printf 'META\tsource\t%s\n' '/staging/ALSU-analysis/spring2026/ConvSK_mind20_full.genome'
printf 'META\tsamples\t%s\n' "$sample_count"
printf 'META\tsample_pair_count\t%s\n' "$actual_sample_pairs"
printf 'META\texpected_sample_pair_count\t%s\n' "$expected_sample_pairs"
printf 'META\tmarkers_total\t%s\n' "$marker_count_total"
printf 'META\tmarkers_excluded_non_autosomal\t%s\n' "$excluded_non_autosomal"
printf 'META\tmarkers_used_for_ibd\t%s\n' "$marker_count_ibd"

awk 'NR>1 {
  p=$10+0
  if (p < 0.05) bins[1]++
  else if (p < 0.10) bins[2]++
  else if (p < 0.15) bins[3]++
  else if (p < 0.25) bins[4]++
  else if (p < 0.50) bins[5]++
  else if (p < 0.75) bins[6]++
  else if (p < 0.90) bins[7]++
  else if (p < 0.98) bins[8]++
  else bins[9]++
}
END {
  labels[1]="0.00-0.05"; labels[2]="0.05-0.10"; labels[3]="0.10-0.15";
  labels[4]="0.15-0.25"; labels[5]="0.25-0.50"; labels[6]="0.50-0.75";
  labels[7]="0.75-0.90"; labels[8]="0.90-0.98"; labels[9]="0.98-1.00";
    for (bin_index=1; bin_index<=9; bin_index++) printf "BIN\t%s\t%d\n", labels[bin_index], bins[bin_index]+0;
}' ConvSK_mind20_full.genome
"""


def run_remote_script():
    """Run the remote extraction script through the DRAGEN SSH tunnel."""
    remote_script_bytes = REMOTE_SCRIPT.replace("\r\n", "\n").replace("\r", "").encode("utf-8")
    result = subprocess.run(
        [
            r"C:\Program Files\PuTTY\plink.exe",
            "-batch",
            "-i",
            r"C:\Users\T\.ssh\id_ed25519.ppk",
            "-P",
            "2222",
            "copilot@127.0.0.1",
            "bash -s",
        ],
        input=remote_script_bytes,
        capture_output=True,
        timeout=120,
    )

    stdout = result.stdout.decode("utf-8", errors="replace")
    stderr = result.stderr.decode("utf-8", errors="replace")

    if result.returncode != 0:
        print(stdout, file=sys.stderr)
        print(stderr, file=sys.stderr)
        raise RuntimeError("Remote PI_HAT extraction failed")

    return stdout


def parse_remote_output(output):
    """Parse metadata and ordered bin counts from the remote script output."""
    metadata = {}
    counts_by_range = {bin_range: 0 for bin_range, _description in BINS}

    for line in output.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 3 and parts[0] == "META":
            key = parts[1]
            value = parts[2]
            metadata[key] = int(value) if value.isdigit() else value
        elif len(parts) == 3 and parts[0] == "BIN":
            counts_by_range[parts[1]] = int(parts[2])

    total_pairs = sum(counts_by_range.values())
    histogram_data = []
    for bin_range, description in BINS:
        count = counts_by_range[bin_range]
        histogram_data.append(
            {
                "range": bin_range,
                "description": description,
                "count": count,
                "percent": (100 * count / total_pairs) if total_pairs else 0,
            }
        )

    metadata["sample_pair_count_from_bins"] = total_pairs
    return metadata, histogram_data

def generate_plotly_json(metadata, histogram_data):
    """Generate Plotly JSON for interactive histogram."""
    ranges = [item["range"] for item in histogram_data]
    counts = [item["count"] for item in histogram_data]
    total_pairs = metadata["sample_pair_count_from_bins"]

    hover_texts = [
        f"{item['range']}: {item['count']:,} sample pairs<br>"
        f"{item['percent']:.4f}% of {total_pairs:,} comparisons<br>"
        f"{item['description']}"
        for item in histogram_data
    ]

    plotly_data = {
        "type": "bar",
        "x": ranges,
        "y": counts,
        "hovertext": hover_texts,
        "hoverinfo": "text",
        "marker": {
            "color": [
                # Use different colors: cool colors for unrelated, warm for related, red for duplicates
                "rgba(100,150,200,0.8)",  # 0.00-0.05
                "rgba(120,160,210,0.8)",  # 0.05-0.10
                "rgba(140,170,220,0.8)",  # 0.10-0.15
                "rgba(180,200,230,0.8)",  # 0.15-0.25
                "rgba(220,180,100,0.8)",  # 0.25-0.50
                "rgba(230,150,80,0.8)",   # 0.50-0.75
                "rgba(240,120,60,0.8)",   # 0.75-0.90
                "rgba(255,100,40,0.8)",   # 0.90-0.98
                "rgba(198,40,40,0.9)"     # 0.98-1.00 (red)
            ]
        }
    }
    
    layout = {
        "title": {
            "text": "Sample PI_HAT Distribution",
            "x": 0.02,
            "xanchor": "left",
            "font": {"size": 16}
        },
        "xaxis": {
            "title": "PI_HAT Value Range",
            "tickangle": -45
        },
        "yaxis": {
            "title": "Number of Sample Pairs"
        },
        "hovermode": "closest",
        "showlegend": False,
        "height": 400,
        "margin": {"b": 100, "t": 70}
    }
    
    return {
        "metadata": metadata,
        "bins": histogram_data,
        "data": [plotly_data],
        "layout": layout
    }

def main():
    print("Fetching Spring 2026 Step 2 sample-pair PI_HAT distribution...")
    remote_output = run_remote_script()
    metadata, histogram_data = parse_remote_output(remote_output)

    total_pairs = metadata["sample_pair_count_from_bins"]
    print(
        f"Samples: {metadata['samples']:,}; sample pairs: {total_pairs:,}; "
        f"markers used for IBD: {metadata['markers_used_for_ibd']:,}"
    )
    for item in histogram_data:
        print(f"  {item['range']}: {item['count']:7,} sample pairs ({item['percent']:.4f}%)")

    print("\nGenerating Plotly JSON...")
    plotly_json = generate_plotly_json(metadata, histogram_data)

    output_file = "data/pihat_distribution_histogram.json"
    with open(output_file, 'w') as f:
        json.dump(plotly_json, f, indent=2)

    print(f"✓ Saved to {output_file}")

if __name__ == "__main__":
    main()
