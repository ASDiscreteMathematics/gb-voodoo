import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap

load('f4_5.sage')
load('common.sage')

f5  = F5()
field_size = 101
num_systems = 50
all_num_terms = [3, 5, 11]
all_num_vars = [2, 4, 6]
all_deg_polys = [2, 5, 7]

# stuff for plotting
top_margin, right_margin, bot_margin, left_margin = 1, 0.05, 1, 0.05
clr_dict = {'red':   [(0.0,  0.0, 0.0), (0.5,  1.0, 1.0), (1.0,  1.0, 1.0)],
            'green': [(0.0,  0.0, 0.0), (0.33, 1.0, 1.0), (0.66, 1.0, 1.0), (1.0,  0.0, 0.0)],
            'blue':  [(0.0,  1.0, 1.0), (0.5,  0.0, 0.0), (1.0,  0.0, 0.0)],}

for num_terms in all_num_terms:
    title = f"field_{field_size}-polys_{num_terms}_terms"
    plt.close()
    plt.title(title)
    fig, axs = plt.subplots(nrows=len(all_num_vars), ncols=len(all_deg_polys), sharex=False, sharey=False, squeeze=False, figsize=(15, 8))
    for ax, d in zip(axs[0], all_deg_polys): # top annotation
        ax.text(0.5, 1.1, f"degree {d}", transform=ax.transAxes, horizontalalignment='center')
    for ax, v in zip(axs[:,-1], all_num_vars): # right annotation
        ax.text(1.3, 0.5, f"{v} vars & polys", transform=ax.transAxes, verticalalignment='center', rotation=270)
    for ax in axs[-1]: # bottom annotation
        ax.set_xlabel("non-zero coeffs")
    for ax in axs[:,0]: # left annotation
        ax.set_ylabel("reductions")

    for subplot_r, num_vars in enumerate(all_num_vars):
        for subplot_c, deg_polys in enumerate(all_deg_polys):
            R = PolynomialRing(GF(field_size), 'x', num_vars)
            reductions = []
            involvements = []
            gb_sizes = []
            num_non_zero_coeffs = []
            ax = axs[subplot_r, subplot_c]
            for _ in range(num_systems):
                polys = [R.random_element(deg_polys, num_terms) for _ in range(num_vars)]
                gb, voos = f5(polys, homogenize=False)
                reductions += [f5.reductions]
                involvements += [mean([sum([1 for x in voo if x]) - 1 for voo in voos]) / (len(voos[0]) - 1)]
                num_non_zero_coeffs += [sum([len(v.coefficients()) for voo in voos for v in voo if v])] # all non-zero coefficients in transformation matrix
                gb_sizes += [len(gb)]
            ax.set_ylim( -bot_margin, max(reductions) + top_margin )
            ax.set_xlim( -left_margin, max(num_non_zero_coeffs) + right_margin )
            clr_map = LinearSegmentedColormap('blue_orange', clr_dict, max(gb_sizes) - min(gb_sizes) + 1) # interpolate colors…
            clr_map = ListedColormap([clr_map((d - min(gb_sizes)) / (max(gb_sizes) - min(gb_sizes) + 1)) for d in range(min(gb_sizes), max(gb_sizes) + 1)]) # …then make discrete
            sctr = ax.scatter(num_non_zero_coeffs, reductions, marker='.', color=[clr_map(d-min(gb_sizes)) for d in gb_sizes])
            bounds = range(min(gb_sizes)-1, max(gb_sizes) + 1)
            ax_clrbar = plt.colorbar(cm.ScalarMappable(norm=None, cmap=clr_map), ax=ax, boundaries=bounds, drawedges=True, aspect=50)
            ax_clrbar.set_label("len(gb)", position=(1.1,0.5), verticalalignment='center', rotation=270)
            tick_skip = len(bounds)//10 + 1
            ax_clrbar.set_ticks([b - 0.5 for b in bounds[::tick_skip]])
            ax_clrbar.set_ticklabels([b for b in bounds[1::tick_skip]])
            plt.savefig(f"invlv_plots/{title}.png", format="png", dpi=300, bbox_inches="tight")
