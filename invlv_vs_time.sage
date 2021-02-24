import numpy as np
from matplotlib import pyplot as plt

load('f4_5.sage')
load('common.sage')

f5  = F5()
field_size = 101
num_systems = 50
all_num_polys_per_system = [2, 3, 5]
all_num_terms = [3, 5, 11]
all_num_vars = [2, 4, 6]
all_deg_polys = [2, 5, 7]

# stuff for plotting
top_margin, right_margin, bot_margin, left_margin = 1, 0.05, 1, 0.05

for num_polys_per_system in all_num_polys_per_system:
    for num_terms in all_num_terms:
        title = f"field: {field_size} – system: {num_polys_per_system} polys, {num_terms} terms"
        plt.close()
        plt.title(title)
        fig, axs = plt.subplots(nrows=len(all_num_vars), ncols=len(all_deg_polys), sharex=True, sharey=True, squeeze=False, figsize=(15, 8))
        axs[0,0].set_xlim( -left_margin, 1 + right_margin )
        for ax, d in zip(axs[0], all_deg_polys): # top annotation
            ax.text(0.5, 1.1, f"degree {d}", transform=ax.transAxes, horizontalalignment='center')
        for ax, v in zip(axs[:,-1], all_num_vars): # right annotation
            ax.text(1.05, 0.5, f"{v} vars", transform=ax.transAxes, verticalalignment='center', rotation=270)
        for ax in axs[-1]: # bottom annotation
            ax.set_xlabel("involvement")
        for ax in axs[:,0]: # left annotation
            ax.set_ylabel("time [s]")
        max_f5_time = 0

        for subplot_r, num_vars in enumerate(all_num_vars):
            for subplot_c, deg_polys in enumerate(all_deg_polys):
                R = PolynomialRing(GF(field_size), 'x', num_vars)
                times = []
                involvements = []
                ax = axs[subplot_r, subplot_c]
                for _ in range(num_systems):
                    polys = [R.random_element(deg_polys, num_terms) for _ in range(num_polys_per_system)]
                    _, voos = f5(polys, homogenize=False)
                    times += [f5.time_gb]
                    involvements += [mean([sum([1 for x in voo if x]) - 1 for voo in voos]) / (len(voos[0]) - 1)]
                max_f5_time = max(times + [max_f5_time])
                ax.set_ylim( -bot_margin, max_f5_time + top_margin )
                ax.plot(involvements, times, 'b.')
                plt.savefig(f"invlv_plots/{title}.png", format="png", dpi=300, bbox_inches="tight")
