
import argparse
import scipy.io as sio
import numpy as np
import time
import pickle
import math as mt
from utilityfunct_aoa_toa_doppler import build_aoa_matrix, build_toa_matrix
from utilityfunct_md_track import md_track_2d
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.size'] = 12


def plot_mDtrack(paths_refined_amplitude_array, paths_refined_toa_array, paths_refined_aoa_array, new_fig=True,
                 colorbar=True, x_lim_down=-1, x_lim_up=20, x_sep=5):
    # plot mdTrack
    paths_power = np.power(np.abs(paths_refined_amplitude_array), 2)
    paths_power = 10 * np.log10(paths_power / np.amax(np.nan_to_num(paths_power)))  # dB
    vmin = -40  # min(np.reshape(paths_power, -1))
    vmax = 0  # max(np.reshape(paths_power, -1))
    if new_fig:
        plt.figure(figsize=(5, 4))
    toa_array = paths_refined_toa_array #- paths_refined_toa_array[0]
    plt.scatter(toa_array * 1E9, paths_refined_aoa_array,
                c=paths_power,
                marker='o', cmap='Blues', s=20,
                vmin=vmin, vmax=vmax)
    if colorbar:
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('power [dB]', rotation=90)
    plt.xlabel('ToA [ns]')
    plt.ylabel('AoA [deg]')
    plt.xlim([x_lim_down, x_lim_up])  # range_considered + 100 * delta_t])
    plt.xticks(np.arange(x_lim_down, x_lim_up, x_sep))
    plt.ylim([-90, 90])
    plt.yticks(np.arange(-90, 91, 20))
    plt.grid(True)
    plt.title('mDtrack')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dir', help='Directory of data')
    parser.add_argument('sub_dir', help='Sub-directory of data')
    parser.add_argument('nss', help='Number of spatial streams', type=int)
    parser.add_argument('ncore', help='Number of cores', type=int)
    args = parser.parse_args()

    exp_dir = args.dir + args.sub_dir + '/'
    names = []

    num_ant = args.nss * args.ncore

    name_file = exp_dir + 'csi_data_calibrated.mat'
    csi_buff = sio.loadmat(name_file)
    H_complete = (csi_buff['csi_data'])

    plt.figure()
    plt.plot(np.real(H_complete[100, :, 0]))
    # plt.plot(np.imag(H_complete[100, :, 0]))
    plt.show()

    F_frequency = 256
    delta_f = 312.5E3
    frequency_vector_idx = np.arange(F_frequency)
    frequency_vector_hz = delta_f * (frequency_vector_idx - F_frequency / 2)
    # control_subcarriers = [0, 1, 2, 3, 4, 5, 127, 128, 129, 251, 252, 253, 254, 255]
    # pilot_subcarriers = [25, 53, 89, 117, 139, 167, 203, 231]
    delete_idxs = np.asarray([0, 1, 2, 3, 4, 5, 25, 53, 89, 117, 127, 128, 129, 139, 167, 203, 231, 251, 252, 253,
                              254, 255], dtype=int)
    frequency_vector_idx = np.delete(frequency_vector_idx, delete_idxs, axis=0)
    frequency_vector_hz = np.delete(frequency_vector_hz, delete_idxs, axis=0)
    fc = 5775E6
    frequency_vector_hz = frequency_vector_hz + fc

    H_complete_valid = H_complete[:, frequency_vector_idx, :]

    T = 3.2e-6
    delta_t = 5e-10
    range_refined_up = 4E-8  # 2.5E-7
    idxs_range_ref_up = int(range_refined_up/delta_t + 1)
    range_refined_down = 1E-8  # 2.5E-7
    idxs_range_ref_down = int(range_refined_down/delta_t + 1)
    t_min = 0  # -T/2
    t_max = T  # T/2

    num_angles = 360
    num_paths = 20
    num_subc = frequency_vector_idx.shape[0]
    ToA_matrix, time_vector = build_toa_matrix(frequency_vector_hz, delta_t, t_min, t_max)
    AoA_matrix, angles_vector, cos_ant_vector = build_aoa_matrix(num_angles, num_ant, frequency_vector_hz)
    AoA_matrix_reshaped = np.reshape(AoA_matrix, (AoA_matrix.shape[0], -1))

    # mD-track 2D: remove offsets CFO, PDD, SFO
    paths_list = []
    paths_amplitude_list = []
    paths_toa_list = []
    paths_aoa_list = []
    num_iteration_refinement = 10
    threshold = -2.5
    for time_idx in range(0, H_complete_valid.shape[0]):
        cfr_sample = H_complete_valid[time_idx, :, :]

        # coarse estimation
        matrix_cfr_toa = np.dot(ToA_matrix, cfr_sample)
        power_matrix_cfr_toa = np.sum(np.abs(matrix_cfr_toa), 1)
        time_idx_max = np.argmax(power_matrix_cfr_toa)
        index_start_toa = int(max(0, time_idx_max - idxs_range_ref_down))
        index_end_toa = int(min(time_vector.shape[0], time_idx_max + idxs_range_ref_up))
        ToA_matrix_considered = ToA_matrix[index_start_toa:index_end_toa, :]
        time_vector_considered = time_vector[index_start_toa:index_end_toa]

        # start = time.time()
        paths, paths_refined_amplitude, paths_refined_toa_idx, paths_refined_aoa_idx = \
            md_track_2d(cfr_sample, AoA_matrix, ToA_matrix_considered, num_ant, num_subc, num_angles,
                        num_iteration_refinement, threshold)
        # end = time.time()
        # print(end-start)
        paths_refined_aoa = angles_vector[paths_refined_aoa_idx] * 180 / mt.pi
        paths_refined_toa = time_vector_considered[paths_refined_toa_idx]
        paths_refined_amplitude_array = np.asarray(paths_refined_amplitude)
        paths_refined_aoa_array = np.asarray(paths_refined_aoa)
        paths_refined_toa_array = np.asarray(paths_refined_toa)

        paths_list.append(paths)
        paths_amplitude_list.append(paths_refined_amplitude_array)
        paths_aoa_list.append(paths_refined_aoa_array)
        paths_toa_list.append(paths_refined_toa_array)

    # name_file = exp_dir + 'paths_list.m'
    # paths_dict = {'path': paths_list}
    # with open(name_file, "wb") as fp:  # Pickling
    #     sio.savemat(fp, paths_dict)
    name_file = exp_dir + 'paths_amplitude_list.m'
    paths_dict = {'amplitude': paths_amplitude_list}
    with open(name_file, "wb") as fp:  # Pickling
        sio.savemat(fp, paths_dict)
    name_file = exp_dir + 'paths_aoa_list.m'
    paths_dict = {'aoa': paths_aoa_list}
    with open(name_file, "wb") as fp:  # Pickling
        sio.savemat(fp, paths_dict)
    name_file = exp_dir + 'paths_toa_list.m'
    paths_dict = {'toa': paths_toa_list}
    with open(name_file, "wb") as fp:  # Pickling
        sio.savemat(fp, paths_dict)
    # plot_mDtrack(paths_refined_amplitude_array, paths_refined_toa_array - paths_refined_toa_array[0],
    # paths_refined_aoa_array, x_lim_down=-10, x_lim_up=40)

    # cfr_cumulative_paths = sum(paths)
    # plt.figure()
    # plt.plot(abs(cfr_cumulative_paths[:, 3]))
    # plt.show()
    #
    # plt.figure()
    # plt.plot(np.abs(H_complete_sanitized[0, :, 2]), color='C0')
    # plt.show()
    #
    # from scipy.fftpack import ifft
    # cir_sample = ifft(cfr_sample, axis=0)
    # plt.figure()
    # plt.plot(np.abs(cir_sample[:, 0]), color='C0')
    # plt.plot(np.abs(cir_sample[:, 1]), color='C1')
    # plt.plot(np.abs(cir_sample[:, 2]), color='C2')
    # plt.plot(np.abs(cir_sample[:, 3]), color='C3')
    # plt.show()

    # import matplotlib.pyplot as plt
    # fig = plt.figure(figsize=(7, 6))
    # ax = fig.add_subplot(111)
    # paths_power = np.power(np.abs(paths_refined_amplitude_array), 2)
    # paths_power = 10 * np.log10(paths_power / np.amax(np.nan_to_num(paths_power)))  # dB
    # vmin = min(np.reshape(paths_power, -1))
    # vmax = max(np.reshape(paths_power, -1))
    # scat = ax.scatter(paths_refined_toa_array * 1E9, paths_refined_aoa_array * 180 / mt.pi,
    #                   c=paths_power, marker='o', cmap='Blues', s=12,
    #                   vmin=vmin, vmax=vmax)
    # cbar = fig.colorbar(scat)
    # cbar.ax.set_ylabel('power [dB]', rotation=90)
    # ax.set_xlabel('ToA [ns]')
    # ax.set_ylabel('AoA [deg]')
    # # ax.set_xlim([(t_min - 100 * delta_t) * 1E9, 50])  # range_considered + 100 * delta_t])
    # ax.set_ylim([-180, 180])
    # ax.grid()
    # plt.tight_layout()
    # plt.show()
