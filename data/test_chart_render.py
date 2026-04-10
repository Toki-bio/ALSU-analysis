"""Generate static test renders of the problematic charts to debug layout."""
import json
import plotly.graph_objects as go
from plotly.subplots import make_subplots

d = json.load(open('data/investigation_data_v2.json'))
ES = d['enriched_samples']

# Chart 2: Quality vs Position — current approach (broken)
fmiss = [(s['idx'], s['fmiss']) for s in ES if s['fmiss'] is not None]
fmiss.sort()

fig_broken = go.Figure()
fig_broken.add_trace(go.Scattergl(
    x=[f[0] for f in fmiss], y=[f[1] for f in fmiss],
    mode='markers', marker=dict(size=4, color='#58a6ff', opacity=0.7)
))
fig_broken.update_layout(
    title='BROKEN: Quality vs Position (single axis)',
    height=400, template='plotly_dark',
    xaxis_title='FID', yaxis_title='F_MISS'
)
fig_broken.write_image('data/test_chart_broken.png', width=1200, height=400)

# Fix approach: split into two panels (top = outliers, bottom = main cluster)
fig_fixed = make_subplots(
    rows=2, cols=1, shared_xaxes=True,
    row_heights=[0.35, 0.65],
    vertical_spacing=0.03,
    subplot_titles=('High missingness samples (F_MISS > 0.06)', 'Main cluster (F_MISS ≤ 0.06)')
)

main = [(s['idx'], s['fmiss'], s['iid'], s.get('barcode','')) for s in ES if s['fmiss'] is not None and s['fmiss'] <= 0.06]
outliers = [(s['idx'], s['fmiss'], s['iid'], s.get('barcode','')) for s in ES if s['fmiss'] is not None and s['fmiss'] > 0.06]

# Outliers panel (top)
bad_chip = [o for o in outliers if o[3].startswith('20899303')]
other_out = [o for o in outliers if not o[3].startswith('20899303')]

if bad_chip:
    fig_fixed.add_trace(go.Scattergl(
        x=[o[0] for o in bad_chip], y=[o[1] for o in bad_chip],
        mode='markers', marker=dict(size=6, color='#f85149', opacity=0.8),
        name='208993030xxx chips', text=[f'{o[2]}<br>F_MISS={o[1]:.4f}<br>{o[3]}' for o in bad_chip],
        hoverinfo='text'
    ), row=1, col=1)

if other_out:
    fig_fixed.add_trace(go.Scattergl(
        x=[o[0] for o in other_out], y=[o[1] for o in other_out],
        mode='markers', marker=dict(size=6, color='#d29922', opacity=0.8),
        name='Other degraded', text=[f'{o[2]}<br>F_MISS={o[1]:.4f}<br>{o[3]}' for o in other_out],
        hoverinfo='text'
    ), row=1, col=1)

# Main cluster (bottom)  
fig_fixed.add_trace(go.Scattergl(
    x=[m[0] for m in main], y=[m[1] for m in main],
    mode='markers', marker=dict(size=4, color='#58a6ff', opacity=0.7),
    name='Good samples', text=[f'{m[2]}<br>F_MISS={m[1]:.4f}<br>{m[3]}' for m in main],
    hoverinfo='text'
), row=2, col=1)

fig_fixed.update_layout(
    height=700, template='plotly_dark',
    title='FIXED: Split Y-axis — both clusters readable'
)
fig_fixed.update_yaxes(title_text='F_MISS', row=1, col=1)
fig_fixed.update_yaxes(title_text='F_MISS', row=2, col=1)
fig_fixed.update_xaxes(title_text='Sample Index (FID)', row=2, col=1)

# Add threshold line on top panel
fig_fixed.add_hline(y=0.20, line_dash='dash', line_color='#f85149', row=1, col=1)

fig_fixed.write_image('data/test_chart_fixed.png', width=1200, height=700)

# Also test the sex check chart
sex_data = [(s['idx'], s['sex_f'], s['snpsex'], s['iid']) for s in ES if s['sex_f'] is not None]

fig_sex = make_subplots(
    rows=1, cols=2, column_widths=[0.7, 0.3],
    subplot_titles=('Sex F distribution (main range)', 'Histogram of F values')
)

# Main scatter with clamped range
for s in sex_data:
    pass  # just test the concept

sex_main = [s for s in sex_data if s[1] >= -0.8]
sex_out = [s for s in sex_data if s[1] < -0.8]

fig_sex.add_trace(go.Scattergl(
    x=[s[0] for s in sex_main], y=[s[1] for s in sex_main],
    mode='markers', marker=dict(size=4, color='#f778ba', opacity=0.7),
    name='Main'
), row=1, col=1)

fig_sex.add_trace(go.Scattergl(
    x=[s[0] for s in sex_out], y=[-0.8]*len(sex_out),
    mode='markers', marker=dict(size=7, color='#d29922', symbol='x', opacity=0.9),
    name=f'Outliers ({len(sex_out)}, clamped)',
    text=[f'{s[3]}<br>F={s[1]:.3f}' for s in sex_out],
    hoverinfo='text'
), row=1, col=1)

# Histogram
fig_sex.add_trace(go.Histogram(
    y=[s[1] for s in sex_data], nbinsy=50,
    marker_color='#f778ba', opacity=0.7, name='F dist'
), row=1, col=2)

fig_sex.update_layout(height=500, template='plotly_dark', title='Sex Check — Fixed')
fig_sex.update_yaxes(range=[-0.9, 0.9], row=1, col=1)

fig_sex.write_image('data/test_sex_fixed.png', width=1200, height=500)

print('Done. Check data/test_chart_broken.png, test_chart_fixed.png, test_sex_fixed.png')
print(f'Main cluster: {len(main)} samples (F_MISS <= 0.06)')
print(f'Outliers: {len(outliers)} samples (F_MISS > 0.06)')
print(f'Sex main: {len(sex_main)}, outliers: {len(sex_out)}')
