import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.io as sio
import matplotlib.cbook as cbook
from matplotlib.patches import Polygon
import scipy.stats as st
import scipy.optimize as so
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Palatino'
mpl.rcParams['text.usetex'] = 'true'
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Accent.colors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()

    base_folders = ['../emulation_results/020822/change_subcarriers_NCG/',
                    '../emulation_results/020822/change_subcarriers_nlos_NCG/']

    subc_nums = [10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 242]
    subc_nums = np.asarray(subc_nums, dtype=int)

    ber_list_folders = []
    sum_pow_no_attack_vector_folders = []
    sum_pow_attack_vector_folders = []
    for base_fold in base_folders:
        ber_list = []
        sum_pow_no_attack_vector = []
        sum_pow_attack_vector = []
        end_time = 2000  # 19000
        end_iter = 1
        num_entries = end_time * end_iter
        for subc_val in subc_nums:
            name_folder = base_fold + '/no_subc_' + str(subc_val)

            ber_victim = np.zeros((num_entries, 1))
            sum_pow_no_attack = np.zeros((num_entries, 1))
            sum_pow_attack = np.zeros((num_entries, 1))
            index_vector = 0
            for time_idx in range(0, end_time):
                for num_iter in range(0, end_iter):
                    name_file = name_folder + '/results_time' + str(time_idx + 1).zfill(5) + '_iteration' + str(
                        num_iter + 1).zfill(2) + '.mat'
                    try:
                        result_dict = sio.loadmat(name_file)
                    except FileNotFoundError:
                        continue
                    ber_victim[index_vector] = result_dict['ber'][0, 0]
                    sum_pow_no_attack[index_vector] = result_dict['sumPowerOrig'][0, 0]
                    sum_pow_attack[index_vector] = result_dict['sumPowerAttack'][0, 0]
                    index_vector = index_vector + 1
            ber_list.append(ber_victim[:index_vector])
            sum_pow_no_attack_vector.append(sum_pow_no_attack[:index_vector])
            sum_pow_attack_vector.append(sum_pow_attack[:index_vector])

        ber_list_folders.append(ber_list)
        sum_pow_no_attack_vector_folders.append(sum_pow_no_attack_vector)
        sum_pow_attack_vector_folders.append(sum_pow_attack_vector)

    stats_folders = []
    stats_folders_violin = []
    stats_power_folders = []
    stats_power_no_attack_folders = []
    for base_fold_idx in range(len(base_folders)):
        stats = []
        stats_violin = []
        stats_power = []
        stats_power_no_attack = []
        for subc_val in range(subc_nums.shape[0]):
            stats.append(cbook.boxplot_stats(ber_list_folders[base_fold_idx][subc_val], whis=(5, 95))[0])
            stats_violin.append(np.squeeze(ber_list_folders[base_fold_idx][subc_val]))
            stats_power.append(cbook.boxplot_stats(sum_pow_attack_vector_folders[base_fold_idx][subc_val],
                                                   whis=(5, 95))[0])
            stats_power_no_attack.append(cbook.boxplot_stats(sum_pow_no_attack_vector_folders[base_fold_idx][subc_val],
                                                             whis=(5, 95))[0])
        stats_folders.append(stats)
        stats_folders_violin.append(stats_violin)
        stats_power_folders.append(stats_power)
        stats_power_no_attack_folders.append(stats_power_no_attack)

    #################################
    # BOX PLOT BER
    #################################
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    fig.set_size_inches(7, 2)
    plot_type = 'violin'  # in ['violin', 'bar']
    if plot_type == 'violin':
        # Plot boxplots from our computed statistics
        bp = ax.violinplot(stats_folders_violin[0], positions=np.arange(subc_nums.shape[0]) - 0.16,
                              # quantiles=[[0.25, 0.75]] * subc_nums.shape[0],
                              points=500, widths=0.32, showmeans=False, showextrema=False, showmedians=True)
        plt.setp(bp['bodies'], facecolor='C4', edgecolor='black', alpha=1, linewidth=1)
        plt.setp(bp['cmedians'], color='black', linewidth=1.5)

        bp = ax.violinplot(stats_folders_violin[1], positions=np.arange(subc_nums.shape[0]) + 0.16,
                              # quantiles=[[0.25, 0.75]] * subc_nums.shape[0],
                              points=500, widths=0.32, showmeans=False, showextrema=False, showmedians=True)
        plt.setp(bp['bodies'], facecolor='C7', alpha=0.6, edgecolor='black', linewidth=1)
        plt.setp(bp['cmedians'], color='black', linewidth=1.5)

    elif plot_type == 'bar':
        # Plot boxplots from our computed statistics
        bp = ax.bxp(stats_folders[0], positions=np.arange(subc_nums.shape[0]) - 0.16, showfliers=False, widths=0.32,
                       manage_ticks=False)
        plt.setp(bp['boxes'], color='black', linewidth=1)
        plt.setp(bp['medians'], color='black', linewidth=1.5)
        plt.setp(bp['whiskers'], color='black')
        for box in bp['boxes']:
            box_x = []
            box_y = []
            for j in range(5):
                box_x.append(box.get_xdata()[j])
                box_y.append(box.get_ydata()[j])
            box_coords = np.column_stack([box_x, box_y])
            ax.add_patch(Polygon(box_coords, facecolor='black', alpha=0.5))

        bp = ax.bxp(stats_folders[1], positions=np.arange(subc_nums.shape[0]) + 0.16, showfliers=False, widths=0.32,
                       manage_ticks=False)
        plt.setp(bp['boxes'], color='black', linewidth=1)
        plt.setp(bp['medians'], color='black', linewidth=1.5)
        plt.setp(bp['whiskers'], color='black')
        for box in bp['boxes']:
            box_x = []
            box_y = []
            for j in range(5):
                box_x.append(box.get_xdata()[j])
                box_y.append(box.get_ydata()[j])
            box_coords = np.column_stack([box_x, box_y])
            ax.add_patch(Polygon(box_coords, facecolor='black', alpha=0.2))

    custom_lines = [Line2D([0], [0], color='C4', linewidth=4),
                    Line2D([0], [0], color='C7', alpha=0.6, linewidth=4)]
    plt.legend(custom_lines, [r'LOS', r'NLOS'],
               ncol=2, labelspacing=0.2, columnspacing=0.5, fontsize='medium', loc='lower right', handlelength=0.6)

    ax.set_xticks(np.arange(len(subc_nums)))
    ax.set_xticklabels([])
    ax.grid()
    ax.set_ylabel(r'BER')
    ax.set_yticks(np.linspace(0, 0.5, 6))
    ax.set_ylim([0, 0.55])
    ax.set_xlim([-0.5, 12.5])
    plt.xlabel(r'no. poisoned sub-channels $\hat{K}$')
    ax.set_xticks(np.arange(len(subc_nums)))
    ax.set_xticklabels(subc_nums)

    name_fig = './plots/ber_change_subcarriers_nlos.pdf'
    plt.savefig(name_fig)
    plt.close()

    #################################
    # BOX PLOT POWER
    #################################
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    fig.set_size_inches(9, 4)
    # Plot boxplots from our computed statistics
    bp = ax.bxp(stats_power, positions=np.arange(subc_nums.shape[0]), showfliers=False)
    plt.setp(bp['boxes'], color='black', linewidth=1)
    plt.setp(bp['medians'], color='black', linewidth=1.5)
    plt.setp(bp['whiskers'], color='black')

    for box in bp['boxes']:
        box_x = []
        box_y = []
        for j in range(5):
            box_x.append(box.get_xdata()[j])
            box_y.append(box.get_ydata()[j])
        box_coords = np.column_stack([box_x, box_y])
        ax.add_patch(Polygon(box_coords, facecolor='black', alpha=0.3))

    ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
    ax.set_xticklabels(subc_nums)
    plt.grid()
    plt.xlabel(r'no. poisoned sub-channels $\hat{K}$')
    ax.set_ylabel(r'$\sum_k \!{\rm Tr}\!\left[\mathbf{W}_k\mathbf{W}_k^\dag\right]\!/\!P_{\rm max}$')
    yticks = np.linspace(0, 5e5, 9)
    ax.set_yticks(yticks)
    ax.set_xlim([-0.5, 12.5])
    pmax = 1e5
    ax.set_yticklabels(yticks / pmax)
    name_fig = './plots/power_change_subcarriers_nlos.pdf'
    plt.savefig(name_fig)
    plt.close()
