"""
Generate 3 publication-quality figures for BiasForensics RSM manuscript.
"""
import csv, json, sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT_DIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(OUT_DIR, exist_ok=True)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'output')

with open(os.path.join(DATA_DIR, 'bias_forensics_results.csv'), encoding='utf-8') as f:
    reviews = list(csv.DictReader(f))
with open(os.path.join(DATA_DIR, 'bias_forensics_summary.json'), encoding='utf-8') as f:
    summary = json.load(f)

for r in reviews:
    r['k_val'] = int(r['k'])
    r['egger_p_val'] = float(r['egger_p'])
    r['begg_p_val'] = float(r['begg_p'])
    r['tf_unadj'] = float(r['tf_theta_unadj'])
    r['tf_adj'] = float(r['tf_theta_adj'])
    r['tf_k0_val'] = int(r['tf_k0'])
    r['pet_theta'] = float(r['petpeese_theta'])
    r['sel_theta'] = float(r['sel3psm_theta'])
    r['lim_theta'] = float(r['limit_theta'])
    r['n_det'] = int(r['n_detect'])
    r['conc'] = float(r['concordance'])

print(f"Loaded {len(reviews)} reviews")

CLASS_COLORS = {
    'Clean': '#2ecc71',
    'Suspected': '#f39c12',
    'Confirmed': '#e74c3c',
    'Discordant': '#8e44ad',
}

# ─────────────────────────────────────────────────
# FIGURE 1: Bias Fingerprint Heatmap
# ─────────────────────────────────────────────────
print("Generating Figure 1: Bias fingerprint heatmap...")

# Build binary matrix: review x method
# Detection: egger_sig, begg_sig, pcurve_inadequate, sel3psm_sig
# Correction: tf shifted, petpeese shifted, 3psm shifted, limit shifted
method_names = ['Egger', 'Begg', 'P-curve\ninadequate', '3PSM\nselection',
                'Trim-fill\nshifted', 'PET-PEESE\nshifted', '3PSM\ncorrection', 'Limit MA\nshifted']

sorted_reviews = sorted(reviews, key=lambda r: (-r['n_det'], r['conc']))

matrix = np.zeros((len(sorted_reviews), 8))
class_colors_list = []

for i, r in enumerate(sorted_reviews):
    # Detection methods (binary)
    matrix[i, 0] = 1 if r['egger_sig'] == '1' else 0
    matrix[i, 1] = 1 if r['begg_sig'] == '1' else 0
    matrix[i, 2] = 1 if r['pcurve_inadequate'] == '1' else 0
    matrix[i, 3] = 1 if r['sel3psm_sig'] == '1' else 0

    # Correction methods (substantial shift > 20%)
    unadj = r['tf_unadj']
    denom = abs(unadj) if abs(unadj) > 0.01 else 0.01
    matrix[i, 4] = 1 if abs(r['tf_adj'] - unadj) / denom > 0.2 else 0
    matrix[i, 5] = 1 if abs(r['pet_theta'] - unadj) / denom > 0.2 else 0
    matrix[i, 6] = 1 if abs(r['sel_theta'] - unadj) / denom > 0.2 else 0
    matrix[i, 7] = 1 if abs(r['lim_theta'] - unadj) / denom > 0.2 else 0

    class_colors_list.append(CLASS_COLORS[r['bias_class']])

fig, (ax_class, ax_heat) = plt.subplots(1, 2, figsize=(10, 8),
    gridspec_kw={'width_ratios': [0.08, 1]}, sharey=True)

# Classification strip
for i, c in enumerate(class_colors_list):
    ax_class.barh(i, 1, height=1, color=c, edgecolor='none')
ax_class.set_xlim(0, 1)
ax_class.set_xticks([])
ax_class.set_ylabel(f'307 Cochrane Reviews (sorted by detection count)', fontsize=10)
ax_class.invert_yaxis()

# Heatmap
cmap = plt.cm.colors.ListedColormap(['#f0f0f0', '#c0392b'])
ax_heat.imshow(matrix, aspect='auto', cmap=cmap, interpolation='nearest')
ax_heat.set_xticks(range(8))
ax_heat.set_xticklabels(method_names, fontsize=8, rotation=0, ha='center')
ax_heat.set_title('Bias Fingerprint: 307 Reviews x 8 Methods', fontsize=12, pad=12)

# Separator between detection and correction
ax_heat.axvline(x=3.5, color='black', linewidth=2)
ax_heat.text(1.5, -8, 'Detection', ha='center', fontsize=9, fontweight='bold')
ax_heat.text(5.5, -8, 'Correction', ha='center', fontsize=9, fontweight='bold')

# Legend
patches = [mpatches.Patch(color=c, label=f'{k} ({summary["classification_counts"][k]})')
           for k, c in CLASS_COLORS.items()]
fig.legend(handles=patches, loc='lower center', ncol=4, fontsize=9, framealpha=0.9,
           bbox_to_anchor=(0.55, 0.02))

plt.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(os.path.join(OUT_DIR, 'figure1_bias_fingerprint.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure1_bias_fingerprint.pdf'), bbox_inches='tight')
print("  Saved figure1_bias_fingerprint.png/pdf")
plt.close(fig)


# ─────────────────────────────────────────────────
# FIGURE 2: Method Agreement Matrix
# ─────────────────────────────────────────────────
print("Generating Figure 2: Method agreement matrix...")

# Compute pairwise agreement for 6 methods (Egger, Begg, TF, PET-PEESE, 3PSM, Limit)
method_labels = ['Egger', 'Begg', 'Trim-fill', 'PET-PEESE', '3PSM', 'Limit']

def get_detection(r, method):
    unadj = r['tf_unadj']
    denom = abs(unadj) if abs(unadj) > 0.01 else 0.01
    if method == 'Egger':
        return r['egger_sig'] == '1'
    elif method == 'Begg':
        return r['begg_sig'] == '1'
    elif method == 'Trim-fill':
        return abs(r['tf_adj'] - unadj) / denom > 0.2
    elif method == 'PET-PEESE':
        return abs(r['pet_theta'] - unadj) / denom > 0.2
    elif method == '3PSM':
        return abs(r['sel_theta'] - unadj) / denom > 0.2
    elif method == 'Limit':
        return abs(r['lim_theta'] - unadj) / denom > 0.2
    return False

n = len(method_labels)
agreement = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        agree_count = sum(1 for r in reviews if get_detection(r, method_labels[i]) == get_detection(r, method_labels[j]))
        agreement[i, j] = agree_count / len(reviews) * 100

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(agreement, cmap='RdYlGn', vmin=50, vmax=100)

for i in range(n):
    for j in range(n):
        color = 'white' if agreement[i, j] < 70 else 'black'
        ax.text(j, i, f'{agreement[i,j]:.1f}', ha='center', va='center', fontsize=10,
                fontweight='bold' if i == j else 'normal', color=color)

ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(method_labels, fontsize=10, rotation=30, ha='right')
ax.set_yticklabels(method_labels, fontsize=10)
ax.set_title('Pairwise Method Agreement (% of 307 Reviews)', fontsize=12, pad=12)

plt.colorbar(im, ax=ax, label='Agreement (%)', shrink=0.8)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'figure2_agreement_matrix.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure2_agreement_matrix.pdf'), bbox_inches='tight')
print("  Saved figure2_agreement_matrix.png/pdf")
plt.close(fig)


# ─────────────────────────────────────────────────
# FIGURE 3: Effect Shift Scatter (4 panels)
# ─────────────────────────────────────────────────
print("Generating Figure 3: Effect shift scatter plots...")

fig, axes = plt.subplots(2, 2, figsize=(10, 9))
methods = [
    ('Trim-and-Fill', 'tf_adj'),
    ('PET-PEESE', 'pet_theta'),
    ('Selection Model', 'sel_theta'),
    ('Limit Meta-Analysis', 'lim_theta'),
]

for ax, (name, key) in zip(axes.flat, methods):
    unadj = [r['tf_unadj'] for r in reviews]
    adj = [r[key] for r in reviews]
    colors = [CLASS_COLORS[r['bias_class']] for r in reviews]

    ax.scatter(unadj, adj, c=colors, s=15, alpha=0.6, edgecolors='white', linewidth=0.3)

    # Identity line
    lims = [min(min(unadj), min(adj)) - 0.1, max(max(unadj), max(adj)) + 0.1]
    ax.plot(lims, lims, 'k--', linewidth=0.8, alpha=0.5, label='No shift')
    ax.axhline(y=0, color='grey', linewidth=0.3)
    ax.axvline(x=0, color='grey', linewidth=0.3)

    ax.set_xlabel('Unadjusted effect', fontsize=9)
    ax.set_ylabel(f'{name} adjusted', fontsize=9)
    ax.set_title(name, fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Shared legend
patches = [mpatches.Patch(color=c, label=k) for k, c in CLASS_COLORS.items()]
fig.legend(handles=patches, loc='lower center', ncol=4, fontsize=9, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.01))
fig.suptitle('Unadjusted vs Bias-Corrected Effect Estimates\nAcross Four Correction Methods',
             fontsize=13, y=1.01)

plt.tight_layout(rect=[0, 0.04, 1, 0.98])
fig.savefig(os.path.join(OUT_DIR, 'figure3_effect_shifts.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure3_effect_shifts.pdf'), bbox_inches='tight')
print("  Saved figure3_effect_shifts.png/pdf")
plt.close(fig)

print(f"\nAll figures saved to {OUT_DIR}/")
for f in sorted(os.listdir(OUT_DIR)):
    size = os.path.getsize(os.path.join(OUT_DIR, f))
    print(f"  {f} ({size/1024:.0f} KB)")
