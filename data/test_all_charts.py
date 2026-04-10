"""Render all key charts from investigation_data_v2.json to verify readability."""
import json
import plotly.graph_objects as go

d = json.load(open('data/investigation_data_v2.json'))
ES = d['enriched_samples']

# ---- Chart 2: Split Quality ----
SPLIT = 0.06
main = [s for s in ES if s['fmiss'] is not None and s['fmiss'] <= SPLIT]
outliers = [s for s in ES if s['fmiss'] is not None and s['fmiss'] > SPLIT]

# Top panel (outliers)
bad = [s for s in outliers if (s.get('barcode') or '').startswith('20899303')]
other = [s for s in outliers if not (s.get('barcode') or '').startswith('20899303')]

fig = go.Figure()
fig.add_trace(go.Scattergl(x=[s['idx'] for s in bad], y=[s['fmiss'] for s in bad],
    mode='markers', marker=dict(size=6, color='#f85149'), name='208993030xxx'))
if other:
    fig.add_trace(go.Scattergl(x=[s['idx'] for s in other], y=[s['fmiss'] for s in other],
        mode='markers', marker=dict(size=6, color='#d29922'), name='Other degraded'))
fig.update_layout(height=350, template='plotly_dark',
    title=f'Top panel: {len(outliers)} outlier samples (F_MISS > {SPLIT})',
    xaxis=dict(range=[0,1260]), yaxis_title='F_MISS')
fig.add_hline(y=0.20, line_dash='dash', line_color='#f85149')
fig.write_image('data/verify_chart2_top.png', width=1200, height=350)

fig2 = go.Figure()
fig2.add_trace(go.Scattergl(x=[s['idx'] for s in main], y=[s['fmiss'] for s in main],
    mode='markers', marker=dict(size=4, color='#58a6ff', opacity=0.7), name='Good'))
fig2.update_layout(height=400, template='plotly_dark',
    title=f'Bottom panel: {len(main)} main cluster (F_MISS ≤ {SPLIT})',
    xaxis=dict(range=[0,1260]), yaxis_title='F_MISS')
fig2.write_image('data/verify_chart2_bot.png', width=1200, height=400)

# ---- Chart 9: Sex check ----
sex = [(s['idx'], s['sex_f'], s.get('snpsex',0), s['iid']) for s in ES if s['sex_f'] is not None]
sex_main = [s for s in sex if -1 <= s[1] <= 1.2]
sex_out = [s for s in sex if s[1] < -1 or s[1] > 1.2]

fig3 = go.Figure()
fig3.add_trace(go.Scattergl(
    x=[s[0] for s in sex_main], y=[s[1] for s in sex_main],
    mode='markers', marker=dict(size=5, color=['#f778ba' if s[2]==2 else '#58a6ff' if s[2]==1 else '#d29922' for s in sex_main], opacity=0.7),
    name='Main'))
fig3.add_trace(go.Scattergl(
    x=[s[0] for s in sex_out], y=[max(s[1], -0.8) for s in sex_out],
    mode='markers', marker=dict(size=7, color='#d29922', symbol='x'), name=f'{len(sex_out)} outliers (clamped)'))
fig3.update_layout(height=450, template='plotly_dark',
    title='Sex F scatter (clamped to [-0.9, 0.9])',
    yaxis=dict(range=[-0.9, 0.9]), xaxis_title='Sample Index')
fig3.add_hline(y=0.8, line_dash='dot', line_color='#58a6ff')
fig3.add_hline(y=0.2, line_dash='dot', line_color='#f778ba')
fig3.write_image('data/verify_chart9_scatter.png', width=800, height=450)

# Sex histogram
allF = [s[1] for s in sex]
fig4 = go.Figure()
fig4.add_trace(go.Histogram(x=allF, nbinsx=60, marker_color='#f778ba'))
fig4.update_layout(height=450, template='plotly_dark',
    title='Sex F histogram', xaxis_title='F value', yaxis_title='Count')
fig4.add_vline(x=0.2, line_dash='dot', line_color='#f778ba')
fig4.add_vline(x=0.8, line_dash='dot', line_color='#58a6ff')
fig4.write_image('data/verify_chart9_hist.png', width=400, height=450)

# ---- Chart 10: Het histogram ----
het = [s['het_rate'] for s in ES if s.get('het_rate') is not None]
fig5 = go.Figure()
fig5.add_trace(go.Histogram(x=het, nbinsx=80, marker_color='#58a6ff'))
fig5.update_layout(height=400, template='plotly_dark',
    title='Heterozygosity distribution', xaxis_title='Het Rate', yaxis_title='Count')
fig5.write_image('data/verify_chart10_het.png', width=600, height=400)

print(f'Chart2 top: {len(outliers)} outliers, bot: {len(main)} main')
print(f'Chart9: {len(sex_main)} main, {len(sex_out)} outliers')
print(f'Chart10: {len(het)} het values')
print('All verification images saved to data/verify_*.png')
