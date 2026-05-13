#!/usr/bin/env python3
"""
Extract PI_HAT distribution from ConvSK_mind20.genome and generate JSON for interactive histogram.
Uses SSH to fetch data from remote DRAGEN server.
"""

import subprocess
import json
import sys

def fetch_genome_file_via_ssh():
    """Fetch and parse ConvSK_mind20.genome from DRAGEN server via SSH tunnel (port 2222)."""
    
    # Command to run on remote server: extract all pairs with their PI_HAT values
    # Try post-imputation IBD first (has more pair data), then fall back to pre-imputation
    remote_cmd = """
# Try to find the best genome file (prefer post-imputation with more pairs)
FILE=""
for pattern in "UZB_v2_IBD.genome" "ConvSK_mind20_ibd.genome" "ConvSK_mind20.genome"; do
  f=$(find /staging/ALSU-analysis -maxdepth 5 -name "$pattern" -type f 2>/dev/null | head -1)
  if [ -f "$f" ]; then
    FILE="$f"
    echo "Using: $f" >&2
    break
  fi
done

if [ -z "$FILE" ]; then
  echo "ERROR: No genome file found" >&2
  exit 1
fi

awk 'NR>1{
  p=$10+0
  if(p<0.05) b["0.00-0.05"]++
  else if(p<0.10) b["0.05-0.10"]++
  else if(p<0.15) b["0.10-0.15"]++
  else if(p<0.25) b["0.15-0.25"]++
  else if(p<0.50) b["0.25-0.50"]++
  else if(p<0.75) b["0.50-0.75"]++
  else if(p<0.90) b["0.75-0.90"]++
  else if(p<0.98) b["0.90-0.98"]++
  else b["0.98-1.00"]++
} END{for(k in b) print k, b[k]}' "$FILE"
"""
    
    try:
        # Execute via plink SSH tunnel (port 2222, copilot@127.0.0.1)
        result = subprocess.run(
            [
                "plink.exe", "-batch",
                "-i", "C:\\Users\\T\\.ssh\\id_ed25519.ppk",
                "-P", "2222",
                "copilot@127.0.0.1",
                remote_cmd
            ],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            print(f"SSH error: {result.stderr}", file=sys.stderr)
            return None
        
        return result.stdout
    except subprocess.TimeoutExpired:
        print("SSH command timed out", file=sys.stderr)
        return None
    except Exception as e:
        print(f"SSH connection error: {e}", file=sys.stderr)
        return None

def parse_distribution(output):
    """Parse AWK output into histogram bins."""
    bins_order = [
        "0.00-0.05",
        "0.05-0.10",
        "0.10-0.15",
        "0.15-0.25",
        "0.25-0.50",
        "0.50-0.75",
        "0.75-0.90",
        "0.90-0.98",
        "0.98-1.00"
    ]
    
    # Parse output, filtering out SSH banner/version info
    distribution = {}
    for line in output.strip().split('\n'):
        if line.strip():
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    # Try to parse as bin data
                    bin_name = parts[0]
                    count = int(parts[1])
                    # Verify it's a valid bin
                    if bin_name in bins_order:
                        distribution[bin_name] = count
                except (ValueError, IndexError):
                    # Skip lines that don't match expected format (SSH banners, etc.)
                    pass
    
    # Ensure all bins exist (fill with 0 if missing)
    result = []
    total = 0
    for bin_name in bins_order:
        count = distribution.get(bin_name, 0)
        result.append({
            "range": bin_name,
            "count": count
        })
        total += count
    
    return result, total

def generate_plotly_json(histogram_data, total_pairs):
    """Generate Plotly JSON for interactive histogram."""
    
    ranges = [d['range'] for d in histogram_data]
    counts = [d['count'] for d in histogram_data]
    
    # Create descriptive labels for each range
    range_descriptions = {
        "0.00-0.05": "Unrelated",
        "0.05-0.10": "Distant cousins",
        "0.10-0.15": "Second cousins+",
        "0.15-0.25": "First cousins+",
        "0.25-0.50": "Half-siblings / First-degree relatives",
        "0.50-0.75": "Full siblings (expected ~0.50)",
        "0.75-0.90": "Very close relatives / MZ twins",
        "0.90-0.98": "Near-identical / Duplicate candidates",
        "0.98-1.00": "Duplicates / Technical replicates (removed)"
    }
    
    hover_texts = [
        f"{ranges[i]}: {counts[i]:,} pairs<br>({counts[i]/total_pairs*100:.1f}%)<br>{range_descriptions.get(ranges[i], '')}"
        for i in range(len(ranges))
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
            "text": "Complete PI_HAT Distribution (All Sample Pairs)",
            "font": {"size": 16}
        },
        "xaxis": {
            "title": "PI_HAT Value Range",
            "tickangle": -45
        },
        "yaxis": {
            "title": "Number of Pairs"
        },
        "hovermode": "closest",
        "showlegend": False,
        "height": 400,
        "margin": {"b": 100, "t": 60}
    }
    
    return {
        "data": [plotly_data],
        "layout": layout
    }

def main():
    print("Fetching PI_HAT distribution from remote server...")
    output = fetch_genome_file_via_ssh()
    
    if output is None:
        print("Failed to fetch data from remote server", file=sys.stderr)
        sys.exit(1)
    
    print("Parsing distribution...")
    histogram_data, total_pairs = parse_distribution(output)
    
    print(f"Total pairs: {total_pairs:,}")
    for item in histogram_data:
        pct = item['count'] / total_pairs * 100 if total_pairs > 0 else 0
        print(f"  {item['range']}: {item['count']:6,} pairs ({pct:5.1f}%)")
    
    print("\nGenerating Plotly JSON...")
    plotly_json = generate_plotly_json(histogram_data, total_pairs)
    
    # Write to file
    output_file = "data/pihat_distribution_histogram.json"
    with open(output_file, 'w') as f:
        json.dump(plotly_json, f, indent=2)
    
    print(f"✓ Saved to {output_file}")
    
    # Also print for inline embedding if needed
    print("\nJSON (for inline embedding):")
    print(json.dumps(plotly_json))

if __name__ == "__main__":
    main()
