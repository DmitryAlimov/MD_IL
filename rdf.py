import os, sys, time
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Process, Manager
man = Manager()
run_bool = man.list([0 for _ in range(100)])

def rdf(traj, tmc, T, chunk, stride, data_folder, pars, cnt):
    if pars == 'TMP.N-BF.B':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BF and name B')
    if pars == 'TMP.N-BF.F':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BF and name F')
    if pars == 'TMP.N-BMI.C0N':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name C0N')
    if pars == 'TMP.N-BMI.C0J':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name C0J')
    if pars == 'TMP.N-BMI.COG':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name COG')
    if pars == 'TMP.N-BMI.C06':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name C06')
    if pars == 'TMP.N-BMI.N02':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name N02')
    if pars == 'TMP.N-BMI.N00':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name N00')
    if pars == 'TMP.N-BMI.C03_C04':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name C03 C04')
    if pars == 'BF.B-BF.B':
        pairs = traj.topology.select_pairs('resname BF and name B', 'resname BF and name B')
    if pars == 'BMI.C0N-BMI.C0N':
        pairs = traj.topology.select_pairs('resname BMI and name C0N', 'resname BMI and name C0N')
    if pars == 'BMI.C01-BMI.C01':
        pairs = traj.topology.select_pairs('resname BMI and name C01', 'resname BMI and name C01')
    if pars == 'BMI.C0N-BF.B':
        pairs = traj.topology.select_pairs('resname BMI and name C0N', 'resname BF and name B')
    if pars == 'BMI.C0N-BF.F':
        pairs = traj.topology.select_pairs('resname BMI and name C0N', 'resname BF and name F')
    if pars == 'BMI.N02-BF.B':
        pairs = traj.topology.select_pairs('resname BMI and name N02', 'resname BF and name B')
    if pars == 'BMI.N02-BF.F':
        pairs = traj.topology.select_pairs('resname BMI and name N02', 'resname BF and name F')
    if pars == 'BMI.C01-BF.B':
        pairs = traj.topology.select_pairs('resname BMI and name C01', 'resname BF and name B')
    if pars == 'BMI.C01-BF.F':
        pairs = traj.topology.select_pairs('resname BMI and name C01', 'resname BF and name F')
    if pars == 'TMP.N-BMI':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name N02 N00')
        rng, val_N = md.compute_rdf(traj, pairs)
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name H07 H09 H08 H0B H0C H0D H0E H0F H0H H0I H0K H0M H0O H0P H0Q')
        rng, val_H = md.compute_rdf(traj, pairs)
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name C01 C04 C03 C05 C06 C0G C0J C0N')
        rng, val_C = md.compute_rdf(traj, pairs)
        val = 14*val_N*2 + 1*val_H*15 + 12*val_C*8
        return [rng, val]
    if pars == 'TMP.N-BF':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BF and name B')
        rng, val_B = md.compute_rdf(traj, pairs)
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BF and name F')
        rng, val_F = md.compute_rdf(traj, pairs)
        val = 11*val_B + 19*val_F*4
        return [rng, val]

    rng, val = md.compute_rdf(traj, pairs)
    return [[rng, val]]


def average(analysis_type):
    plot_data = []
    for T in [160, 180]:
        rng, val, cnt = [], [] ,0
        for tmc in range(0, 400, 5):
            if os.path.exists(f'/home/student/bmim/analysis/rdf/{analysis_type}/{T}.{tmc}.npy'):
                arr = np.load(f'/home/student/bmim/analysis/rdf/{analysis_type}/{T}.{tmc}.npy')
                rng = arr[0][0]
                if len(val) == 0:
                    val = sum(arr[:, 1])
                else:
                    val += sum(arr[:, 1])
                cnt += 1
        val = val/cnt
        plt.plot(rng, val, label=str(T))

    plt.legend(), plt.title(analysis_type)
    plt.savefig(f'/home/student/bmim/analysis/rdf/plot_data/{analysis_type}.png', dpi=256)
    plot_data.append([rng, val, str(T)])
    np.save(f'/home/student/bmim/analysis/rdf/plot_data/{analysis_type}.npy', plot_data)
    # plt.show()