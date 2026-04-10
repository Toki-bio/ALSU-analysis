"""Generate the replacement JavaScript for step15 ADMIXTURE x ROH scatter chart."""
import json

# Load the compact data
data = json.load(open('data/admix_roh_compact.json'))
q_arr = json.dumps(data['q'], separators=(',', ':'))
f_arr = json.dumps(data['f'], separators=(',', ':'))

old_block = '''// ════════════════════ ADMIXTURE × ROH SCATTER ════════════════════
(function(){
    const canvas = document.getElementById('admixRohScatter');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;

    // TODO: Replace with real Q1 × F_ROH data from server
    // Source: ADMIXTURE .Q files + UZB_ROH.hom.indiv (1,074 individuals)
    // This chart requires per-individual Q1 values and F_ROH values that are only on the Biotech2024 server.

    ctx.fillStyle = '#fafafa';
    ctx.fillRect(0, 0, W, H);

    const pad = { t: 45, r: 30, b: 55, l: 65 };
    const w = W - pad.l - pad.r;
    const h = H - pad.t - pad.b;

    // Grid
    ctx.strokeStyle = '#eee';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 5; i++) {
        const y = pad.t + (i / 5) * h;
        ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(W - pad.r, y); ctx.stroke();
    }
    for (let i = 0; i <= 5; i++) {
        const x = pad.l + (i / 5) * w;
        ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + h); ctx.stroke();
    }

    // "Data pending" label
    ctx.fillStyle = '#999';
    ctx.font = 'bold 16px Segoe UI';
    ctx.textAlign = 'center';
    ctx.fillText('Real data pending \\u2014 requires server-side extraction', W / 2, pad.t + h / 2 - 10);
    ctx.font = '12px Segoe UI';
    ctx.fillText('Source: ADMIXTURE .Q files \\u00d7 UZB_ROH.hom.indiv (1,074 individuals)', W / 2, pad.t + h / 2 + 15);

    // Title
    ctx.fillStyle = '#333';
    ctx.font = 'bold 13px Segoe UI';
    ctx.textAlign = 'center';
    ctx.fillText('ADMIXTURE Q1 (European-like) vs FROH', W / 2, pad.t - 20);

    ctx.fillStyle = '#888';
    ctx.font = '11px Segoe UI';
    ctx.fillText('Q1 (European-like ancestry proportion)', W / 2, H - 8);

    ctx.save();
    ctx.translate(16, pad.t + h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('FROH', 0, 0);
    ctx.restore();

    // X ticks
    for (let i = 0; i <= 5; i++) {
        const v = i / 5;
        const x = pad.l + (v) * w;
        ctx.fillStyle = '#aaa';
        ctx.font = '10px Segoe UI';
        ctx.textAlign = 'center';
        ctx.fillText(v.toFixed(1), x, pad.t + h + 18);
    }

    // Y ticks
    for (let i = 0; i <= 5; i++) {
        const v = (i / 5) * 0.15;
        const y = pad.t + h - (i / 5) * h;
        ctx.fillStyle = '#aaa';
        ctx.font = '10px Segoe UI';
        ctx.textAlign = 'right';
        ctx.fillText(v.toFixed(3), pad.l - 6, y + 4);
    }

    // Legend
    const lx = pad.l + 10;
    const ly = pad.t + 10;
    ctx.fillStyle = 'rgba(255,255,255,.85)';
    ctx.fillRect(lx, ly, 180, 38);
    ctx.fillStyle = 'rgba(0,105,64,.6)';
    ctx.beginPath(); ctx.arc(lx + 8, ly + 10, 3, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#555'; ctx.font = '10px Segoe UI'; ctx.textAlign = 'left';
    ctx.fillText('Background consanguinity', lx + 16, ly + 14);
    ctx.fillStyle = 'rgba(229,57,53,.9)';
    ctx.beginPath(); ctx.arc(lx + 8, ly + 28, 4, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#555';
    ctx.fillText('FROH > 0.0625 (consanguineous)', lx + 16, ly + 32);
})();'''

new_block = f'''// ════════════════════ ADMIXTURE × ROH SCATTER ════════════════════
(function(){{
    const canvas = document.getElementById('admixRohScatter');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;

    // Real data: ADMIXTURE K=2 Q1 (European-like) vs F_ROH
    // Source: admix_v2_data.json (1,047 sample IDs + K2 Q values) × UZB_v2_ROH.hom.indiv
    // F_ROH = ROH_KB / 2,881,033 (autosomal hg38 length)
    const Q = {q_arr};
    const F = {f_arr};

    ctx.fillStyle = '#fafafa';
    ctx.fillRect(0, 0, W, H);

    const pad = {{ t: 45, r: 30, b: 55, l: 65 }};
    const w = W - pad.l - pad.r;
    const h = H - pad.t - pad.b;

    // Axis ranges
    const xMin = 0, xMax = 1.0;
    const yMin = 0, yMax = 0.12;

    // Grid
    ctx.strokeStyle = '#eee';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 6; i++) {{
        const y = pad.t + (i / 6) * h;
        ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(W - pad.r, y); ctx.stroke();
    }}
    for (let i = 0; i <= 5; i++) {{
        const x = pad.l + (i / 5) * w;
        ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + h); ctx.stroke();
    }}

    // Consanguinity threshold line (F_ROH = 0.0625)
    const threshY = pad.t + h - (0.0625 / yMax) * h;
    ctx.strokeStyle = 'rgba(229,57,53,.4)';
    ctx.lineWidth = 1;
    ctx.setLineDash([6, 4]);
    ctx.beginPath(); ctx.moveTo(pad.l, threshY); ctx.lineTo(W - pad.r, threshY); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = 'rgba(229,57,53,.7)';
    ctx.font = '9px Segoe UI';
    ctx.textAlign = 'right';
    ctx.fillText('F_ROH = 0.0625 (1st cousin)', W - pad.r - 4, threshY - 4);

    // 2nd cousin line (F_ROH = 0.0156)
    const thresh2Y = pad.t + h - (0.0156 / yMax) * h;
    ctx.strokeStyle = 'rgba(100,100,100,.3)';
    ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(pad.l, thresh2Y); ctx.lineTo(W - pad.r, thresh2Y); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = 'rgba(100,100,100,.5)';
    ctx.fillText('F_ROH = 0.0156 (2nd cousin)', W - pad.r - 4, thresh2Y - 4);

    // Plot points
    for (let i = 0; i < Q.length; i++) {{
        const x = pad.l + ((Q[i] - xMin) / (xMax - xMin)) * w;
        const y = pad.t + h - ((F[i] - yMin) / (yMax - yMin)) * h;
        const consang = F[i] > 0.0625;
        ctx.fillStyle = consang ? 'rgba(229,57,53,.9)' : 'rgba(0,105,64,.45)';
        ctx.beginPath();
        ctx.arc(x, y, consang ? 4 : 2.5, 0, Math.PI * 2);
        ctx.fill();
    }}

    // Title
    ctx.fillStyle = '#333';
    ctx.font = 'bold 13px Segoe UI';
    ctx.textAlign = 'center';
    ctx.fillText('ADMIXTURE Q1 (European-like) vs F_ROH (N = 1,047)', W / 2, pad.t - 20);

    // Axis labels
    ctx.fillStyle = '#888';
    ctx.font = '11px Segoe UI';
    ctx.fillText('Q1 (European-like ancestry proportion)', W / 2, H - 8);

    ctx.save();
    ctx.translate(16, pad.t + h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('F_ROH', 0, 0);
    ctx.restore();

    // X ticks
    for (let i = 0; i <= 5; i++) {{
        const v = i / 5;
        const x = pad.l + (v) * w;
        ctx.fillStyle = '#aaa';
        ctx.font = '10px Segoe UI';
        ctx.textAlign = 'center';
        ctx.fillText(v.toFixed(1), x, pad.t + h + 18);
    }}

    // Y ticks
    for (let i = 0; i <= 6; i++) {{
        const v = (i / 6) * yMax;
        const y = pad.t + h - (i / 6) * h;
        ctx.fillStyle = '#aaa';
        ctx.font = '10px Segoe UI';
        ctx.textAlign = 'right';
        ctx.fillText(v.toFixed(3), pad.l - 6, y + 4);
    }}

    // Legend
    const lx = pad.l + 10;
    const ly = pad.t + 10;
    ctx.fillStyle = 'rgba(255,255,255,.85)';
    ctx.fillRect(lx, ly, 220, 38);
    ctx.strokeStyle = '#ddd'; ctx.lineWidth = 1;
    ctx.strokeRect(lx, ly, 220, 38);
    ctx.fillStyle = 'rgba(0,105,64,.6)';
    ctx.beginPath(); ctx.arc(lx + 8, ly + 10, 2.5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#555'; ctx.font = '10px Segoe UI'; ctx.textAlign = 'left';
    ctx.fillText('Background (F_ROH \\u2264 0.0625)', lx + 16, ly + 14);
    ctx.fillStyle = 'rgba(229,57,53,.9)';
    ctx.beginPath(); ctx.arc(lx + 8, ly + 28, 4, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#555';
    ctx.fillText('Consanguineous (F_ROH > 0.0625, N = 28)', lx + 16, ly + 32);
}})();'''

# Read step15 and replace
with open('steps/step15.html', 'r', encoding='utf-8') as fh:
    content = fh.read()

if old_block not in content:
    # Try with actual unicode chars instead of escapes
    old_block2 = old_block.replace('\\u2014', '\u2014').replace('\\u00d7', '\u00d7')
    if old_block2 in content:
        content = content.replace(old_block2, new_block)
        print("Replaced (unicode literal match)")
    else:
        print("ERROR: Could not find old block!")
        # Debug: find the sentinel
        idx = content.find('TODO: Replace with real Q1')
        if idx >= 0:
            print(f"Found TODO at position {idx}")
            print("Context:", repr(content[idx-100:idx+100]))
        exit(1)
else:
    content = content.replace(old_block, new_block)
    print("Replaced (escape match)")

with open('steps/step15.html', 'w', encoding='utf-8') as fh:
    fh.write(content)

print("Done! step15.html updated with real scatter data.")
