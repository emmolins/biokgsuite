"""Fig. 7 redrawn from results tables in the house style (single panel).

KG lift (MRR of the grounded arm minus the no-KG arm) across the Llama size
ladder, with paired-bootstrap 95% intervals over the 128 pairs
(results/tables/G04_scaling/capacity.csv: lift, lift_lo, lift_hi).
Points are dodged horizontally so every interval is visible.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])
sys.path.insert(0, str(ROOT / 'src'))
from biokgsuite.viz.style import setup_style, clean_ax, KG_PALETTE, KG_NAMES, FAINT, MUTED

setup_style()
cap = pd.read_csv(ROOT / 'results' / 'tables' / 'G04_scaling' / 'capacity.csv')
cap = cap[cap['arm'] != 'none']
sizes = sorted(cap['params_b'].unique())

fig, ax = plt.subplots(figsize=(89 / 25.4, 62 / 25.4))
dodge = {'primekg': 0.88, 'drkg': 1.0, 'biokg': 1.0 / 0.88}
for kg in ['primekg', 'drkg', 'biokg']:
    d = cap[cap['arm'] == kg].sort_values('params_b')
    x = d['params_b'].to_numpy() * dodge[kg]
    y = d['lift'].to_numpy()
    err = np.vstack([y - d['lift_lo'].to_numpy(), d['lift_hi'].to_numpy() - y])
    c = KG_PALETTE[kg]
    ax.plot(x, y, color=c, lw=1.0, zorder=2)
    ax.errorbar(x, y, yerr=err, fmt='o', ms=3.2, color=c, mfc=c, mec='white', mew=0.4,
                elinewidth=0.8, capsize=1.6, capthick=0.8, zorder=3, label=KG_NAMES[kg])

ax.axhline(0, color=FAINT, lw=0.6, ls='--', zorder=1)
ax.text(sizes[-1] * 1.45, 0, 'no gain', fontsize=6, color=MUTED, va='center', ha='left')
ax.set_xscale('log')
ax.set_xticks(sizes)
ax.set_xticklabels([f'{int(s)}B' for s in sizes])
ax.minorticks_off()
ax.set_xlim(sizes[0] * 0.7, sizes[-1] * 1.4)
ax.set_ylim(-0.05, 0.40)
clean_ax(ax, xlabel='Parameters (billions)', ylabel='MRR gain over no KG')
ax.legend(loc='upper left', frameon=False, fontsize=6.5, handletextpad=0.4, borderaxespad=0.2)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ('pdf', 'png'):
    fig.savefig(OUT / f'Fig7_capacity.{ext}', dpi=600, bbox_inches='tight', pad_inches=0.02)
print('ok')
