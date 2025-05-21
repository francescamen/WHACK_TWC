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
mpl.rcParams['font.size'] = 18
mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Accent.colors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()

    base_folder = '../emulation_results/020822/change_subcarriers_NCG/'

    subc_nums = [10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 242]
    subc_nums = np.asarray(subc_nums, dtype=int)

    end_time = 2000
    end_iter = 1
    num_entries = end_time*end_iter

        ber_list_victim = []
        ber_list_attack = []
        sum_pow_no_attack_vector = []
        sum_pow_attack_vector = []
        num_iterations_vector = []

        for subc_val in subc_nums:
            name_folder = base_folder + '/no_subc_' + str(subc_val)

            ber_victim = np.zeros((num_entries, 1))
            ber_attack = np.zeros((num_entries, 1))
            sum_pow_no_attack = np.zeros((num_entries, 1))
            sum_pow_attack = np.zeros((num_entries, 1))
            num_iterations = np.zeros((num_entries, 1))
            index_vector = 0
            for time_idx in range(0, end_time):
                for num_iter in range(0, end_iter):
                    name_file = name_folder + '/results_time' + str(time_idx+1).zfill(5) + '_iteration' + str(num_iter+1).zfill(2) + '.mat'
                    try:
                        result_dict = sio.loadmat(name_file)
                    except FileNotFoundError:
                        continue
                    ber_victim[index_vector] = result_dict['ber'][0, 0]
                    ber_attack[index_vector] = result_dict['ber'][0, 1]
                    sum_pow_no_attack[index_vector] = result_dict['sumPowerOrig'][0, 0]
                    sum_pow_attack[index_vector] = result_dict['sumPowerAttack'][0, 0]
                    num_iterations[index_vector] = result_dict['num_iterations'][0, 0]
                    index_vector = index_vector + 1
            ber_list_victim.append(ber_victim[:index_vector])
            ber_list_attack.append(ber_attack[:index_vector])
            sum_pow_no_attack_vector.append(sum_pow_no_attack[:index_vector])
            sum_pow_attack_vector.append(sum_pow_attack[:index_vector])
            num_iterations_vector.append(num_iterations[:index_vector])

        stats_victim = []
        stats_victim_violin = []
        stats_attack = []
        stats_attack_violin = []
        stats_power = []
        stats_power_no_attack = []
        stats_num_iterations = []
        for subc_val in range(subc_nums.shape[0]):
            stats_victim.append(cbook.boxplot_stats(ber_list_victim[subc_val], whis=(5, 95))[0])
            stats_victim_violin.append(np.squeeze(ber_list_victim[subc_val]))
            stats_attack.append(cbook.boxplot_stats(ber_list_attack[subc_val], whis=(5, 95))[0])
            stats_attack_violin.append(np.squeeze(ber_list_attack[subc_val]))
            stats_power.append(cbook.boxplot_stats(sum_pow_attack_vector[subc_val], whis=(5, 95))[0])
            stats_power_no_attack.append(cbook.boxplot_stats(sum_pow_no_attack_vector[subc_val], whis=(5, 95))[0])
            stats_num_iterations.append(cbook.boxplot_stats(num_iterations_vector[subc_val], whis=(5, 95))[0])

        #################################
        # BOX PLOT BER
        #################################
        fig, ax = plt.subplots(2, 1, constrained_layout=True)
        fig.set_size_inches(7, 4)
        plot_type = 'violin'  # in ['violin', 'bar']
        if plot_type == 'violin':
            # Plot boxplots from our computed statistics
            bp = ax[0].violinplot(stats_victim_violin, positions=np.arange(subc_nums.shape[0])-0.16,
                                  points=500, widths=0.32, showmeans=False, showextrema=False, showmedians=True)
            plt.setp(bp['bodies'], facecolor='C4', edgecolor='black', alpha=1, linewidth=1)
            plt.setp(bp['cmedians'], color='black', linewidth=1.5)

            bp = ax[0].violinplot(stats_attack_violin, positions=np.arange(subc_nums.shape[0])+0.16,
                                  points=500, widths=0.32, showmeans=False, showextrema=False, showmedians=True)
            plt.setp(bp['bodies'], facecolor='C2', edgecolor='black', alpha=1, linewidth=1)
            plt.setp(bp['cmedians'], color='black', linewidth=1.5)

        elif plot_type == 'bar':
            # Plot boxplots from our computed statistics
            bp = ax[0].bxp(stats_victim, positions=np.arange(subc_nums.shape[0]) - 0.16, showfliers=False, widths=0.32,
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
                ax[0].add_patch(Polygon(box_coords, facecolor='black', alpha=0.5))

            bp = ax[0].bxp(stats_attack, positions=np.arange(subc_nums.shape[0]) + 0.16, showfliers=False, widths=0.32,
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
                ax[0].add_patch(Polygon(box_coords, facecolor='black', alpha=0.1))

        custom_lines = [Line2D([0], [0], color='C4', linewidth=4),
                        Line2D([0], [0], color='C2', linewidth=4)]
        ax[0].legend(custom_lines, [r'victim', r'adversary'],
                     ncol=2, labelspacing=0.2, columnspacing=0.5, fontsize='medium', loc='upper left', handlelength=0.6)

        ax[0].set_xticks(np.arange(len(subc_nums)))
        ax[0].set_xticklabels([])
        ax[0].grid()
        ax[0].set_ylabel(r'BER')
        ax[0].set_yticks(np.linspace(0, 0.5, 6))
        ax[0].set_ylim([0, 0.55])
        ax[0].set_xlim([-0.5, len(subc_nums) - 0.5])

        #################################
        # BOX PLOT POWER
        #################################
        # Plot boxplots from our computed statistics
        bp = ax[1].bxp(stats_power, positions=np.arange(subc_nums.shape[0]), showfliers=False, widths=0.32, manage_ticks=False)
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
            ax[1].add_patch(Polygon(box_coords, color='black', alpha=0.3))

        ax[1].ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
        ax[1].set_xticks(np.arange(len(subc_nums)))
        ax[1].set_xticklabels(subc_nums)
        ax[1].grid()
        ax[1].set_xlabel(r'no. poisoned sub-channels $\hat{K}$')
        ax[1].set_ylabel(r'$P\!/\!P_{\rm max}$')
        yticks = np.round(np.arange(0, 1.1e5, 0.2e5), 1)
        ax[1].set_yticks(yticks)
        pmax = 1e5
        ax[1].set_yticklabels(yticks / pmax)
        ax[1].set_ylim([0e5, 1.2e5])
        ax[1].set_xlim([-0.5, len(subc_nums) - 0.5])

        name_fig = './plots/ber_power_change_subcarriers_victim_attacker_sub_NCG.pdf'
        plt.savefig(name_fig)
        plt.close()
