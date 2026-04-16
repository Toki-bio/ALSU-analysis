#!/bin/bash
set -eo pipefail
PW='d6rFCYsg&q;RZ0'
WORKDIR=/tmp/michigan_results
RESULT=/tmp/spring_stats.txt
cd "$WORKDIR"

echo "=== Spring 2026 Analysis: $(date) ===" > "$RESULT"

# URLs array
declare -a U
U[1]="https://imputationserver.sph.umich.edu/share/results/2a87b30430bbd44939595486534fcb9df55426fd81789165d9077e178daccdff/chr_1.zip"
U[2]="https://imputationserver.sph.umich.edu/share/results/ec0f1ecdd1ae72b052894113a2e77f09f8143ac2d6d6e08b6d0e78c02b62bbbf/chr_2.zip"
U[3]="https://imputationserver.sph.umich.edu/share/results/b33f936ad78fa65a3d1d6d7f977c5b5d2836eba61d08cee0b445ae97e47ee887/chr_3.zip"
U[4]="https://imputationserver.sph.umich.edu/share/results/9c756df94ae17767799cdb4e422f16bba0e141c44eb788c33ba3705780d71082/chr_4.zip"
U[5]="https://imputationserver.sph.umich.edu/share/results/00add947fbf01535557af01de599813357bb4bb102ffcce528175d5c03051be8/chr_5.zip"
U[6]="https://imputationserver.sph.umich.edu/share/results/56ae398622942bb46c32f45363771a820b31990d5cff703d43b360ef814f990e/chr_6.zip"
U[7]="https://imputationserver.sph.umich.edu/share/results/d567576ebc4ff599b584002f33d7e7274d41a08e97923ef2a0ef5fbccb883e0b/chr_7.zip"
U[8]="https://imputationserver.sph.umich.edu/share/results/786395b7feab0acfb8a5f1982c5646fadde2b0264acbf3298e44147a67c6d9d3/chr_8.zip"
U[9]="https://imputationserver.sph.umich.edu/share/results/b9cb5a89a9ffc65d3dadc30732f1cefb8b6b0f4cca50f38876c0636964fe5401/chr_9.zip"
U[10]="https://imputationserver.sph.umich.edu/share/results/e50af01426196de2556bdcce7a4412b1d94d62d8a5fcf06e9d81cefa6f1ccabb/chr_10.zip"
U[11]="https://imputationserver.sph.umich.edu/share/results/d7ede4954d11b2e57880f8570b5939dfc163e4109d32897f4f6e40d22bcc75f7/chr_11.zip"
U[12]="https://imputationserver.sph.umich.edu/share/results/8c20885eef6e66c5656d19663e30772714b4bba31509a7d66dc2dc81bbf4f5f7/chr_12.zip"
U[13]="https://imputationserver.sph.umich.edu/share/results/c469a0d38d7e292f25af6c69a295365fe75910c3bd7f82065f98e0ef3bf61194/chr_13.zip"
U[14]="https://imputationserver.sph.umich.edu/share/results/5b107fa2f3eb5eb1a460e65ff1bc9e7fb74b6a9d8e1e976b986bb1aba11c3757/chr_14.zip"
U[15]="https://imputationserver.sph.umich.edu/share/results/82469309f2752b1dce8bbdff97de4c373be06e3eeef5d45c2e48c9febab1b947/chr_15.zip"
U[16]="https://imputationserver.sph.umich.edu/share/results/9171ef87052810bb33ee58a1798c8fdbf9e0606a28dfef2450fc5870394b14d2/chr_16.zip"
U[17]="https://imputationserver.sph.umich.edu/share/results/8c497ee4bea8569f9b87d334c97022bab7ddaf97f86a23932b5d28f5cb20c85d/chr_17.zip"
U[18]="https://imputationserver.sph.umich.edu/share/results/4d78b2e898814c3135a75bccd005e2f33365406a9ac0a449a7e6845dc8fe8043/chr_18.zip"
U[19]="https://imputationserver.sph.umich.edu/share/results/47566eca76332d15b02cdf8eab040ee5bc968e9fe99a088a8a2090cd288366f1/chr_19.zip"
U[20]="https://imputationserver.sph.umich.edu/share/results/e7e3974d1fb72c7d46f700582645db69da15b4471002887a7a3f941398dd7160/chr_20.zip"
U[21]="https://imputationserver.sph.umich.edu/share/results/674be72e77b0690630a9a25f53b2e57a083d93a7f2ad5874c6e42be8c8ef6727/chr_21.zip"
U[22]="https://imputationserver.sph.umich.edu/share/results/9bb42010b77064e9e0c8b2d96d6d1fbce359a91d3636addebd815fb6f344da57/chr_22.zip"

echo "PER_CHR" >> "$RESULT"

for chr in $(seq 1 22); do
  echo "--- chr${chr}: $(date) ---"
  zipfile="chr_${chr}.zip"
  infofile="chr${chr}.info.gz"
  
  # Download if not present or too small
  if [ ! -f "$zipfile" ] || [ "$(stat -c%s "$zipfile" 2>/dev/null || echo 0)" -lt 1000000 ]; then
    echo "  Downloading $zipfile..."
    curl -sL -o "$zipfile" "${U[$chr]}"
    echo "  Downloaded $(du -h "$zipfile" | cut -f1)"
  else
    echo "  $zipfile exists ($(du -h "$zipfile" | cut -f1))"
  fi
  
  # Extract only .info.gz
  if [ ! -f "$infofile" ]; then
    unzip -P "$PW" -o "$zipfile" "$infofile" >/dev/null 2>&1
    echo "  Extracted $infofile ($(du -h "$infofile" | cut -f1))"
  fi
  
  # R2 distribution for this chr - write per-chr bins directly
  zcat "$infofile" | grep -v '^#' | sed 's/.*R2=\([0-9.]*\).*/\1/' | awk -v chr="$chr" '
  BEGIN { for (i=0; i<10; i++) b[i] = 0 }
  {
    v = $1 + 0; n++; s += v
    bin = int(v * 10); if (bin >= 10) bin = 9
    b[bin]++
    if (v >= 0.3) a++; if (v >= 0.5) e++; if (v >= 0.8) c++; if (v >= 0.9) d++
  }
  END {
    printf "CHR %s N %d SUM %.2f", chr, n, s
    for (i=0; i<10; i++) printf " B%d %d", i, b[i]
    printf " GE30 %d GE50 %d GE80 %d GE90 %d\n", a, e, c, d
  }
  ' >> "$RESULT"
  
  echo "  chr${chr}: done"
  
  # Delete zip (save space) and info
  rm -f "$zipfile" "$infofile"
  echo "  Disk free: $(df -h /tmp | tail -1 | awk '{print $4}')"
done

echo "" >> "$RESULT"
echo "DONE $(date)" >> "$RESULT"
echo "=== Spring analysis complete: $(date) ==="
