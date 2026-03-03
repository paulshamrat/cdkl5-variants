import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import os

class StyleConfig:
    def __init__(self):
        self.base_fontsize   = 14
        self.annot_fontsize  = 14
        self.legend_fontsize = 14
        self.title_pad       = 14
        self.annot_offset    = 0.05
        self.bar_alpha       = 0.3
        self.bar_width       = 0.15
        self.figsize         = (15, 6)
        self.dpi             = 300
        self.y_headroom      = 1.2

def apply_global_style(cfg):
    plt.rcParams.update({'font.size': cfg.base_fontsize})

def annotate_bars(ax, xs, heights, cfg):
    top = max(heights) if len(heights) > 0 else 1
    for x, h in zip(xs, heights):
        ax.text(
            x, h + top * cfg.annot_offset, f"{h:.2f}",
            ha='center', va='bottom', rotation=90,
            clip_on=False, fontsize=cfg.annot_fontsize
        )

def finalize_axes(ax, ytops, cfg):
    max_y = max(ytops) if len(ytops) > 0 else 1
    ax.set_ylim(0, max_y * cfg.y_headroom)
    ax.tick_params(axis='x', labelsize=cfg.base_fontsize * 0.9)
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center')
    ax.grid(False)

def plot_bar_panel(ax, data, title, methods, colors, cfg, bar_width=None):
    """
    Plots a multi-bar panel for thermodynamic energy changes.
    """
    w = bar_width if bar_width is not None else cfg.bar_width
    items = sorted(data.items(), key=lambda x: x[1]['pos'])
    muts = [m for m, _ in items]
    if not muts:
        return
    
    x    = np.arange(len(muts))
    vals = np.array([[d[m] for m in methods] for _, d in items])
    maxv = np.array([d['max'] for _, d in items])

    ax.bar(x, maxv, w * len(methods), color='grey', alpha=cfg.bar_alpha)
    for i in range(len(methods)):
        ax.bar(
            x + (i - (len(methods) - 1) / 2) * w,
            vals[:, i],
            w,
            color=colors[i]
        )

    annotate_bars(ax, x, maxv, cfg)
    ax.set_xticks(x)
    ax.set_xticklabels(muts)
    ax.set_xlabel('Mutation')
    ax.set_ylabel(r'$|\Delta\Delta G|$')
    ax.set_title(title, loc='center', pad=cfg.title_pad)
    
    handles = [Patch(facecolor='grey', alpha=cfg.bar_alpha, label='Picked Max')] + [
        Patch(facecolor=colors[i], label=methods[i]) for i in range(len(methods))
    ]
    ax.legend(handles, [h.get_label() for h in handles],
              loc='upper right', frameon=False, fontsize=cfg.legend_fontsize)
    finalize_axes(ax, list(maxv), cfg)
