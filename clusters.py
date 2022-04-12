import os, sys, time
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.fft import fft, fftfreq
from __additional import smooth, interp
from scipy.interpolate import splev, splrep
from TMP_L import read_xvg, van_hoff
import mdtraj as md
import scipy.linalg as LA
from numba import jit
import networkx as nx
import community.community_louvain
import community


def BMI_chains(traj, tmc, T, chunk, stride, data_folder, pars, cnt):
    ndx = np.array(traj.topology.select('resname BMI and name C0N'))
    return traj.xyz[:, ndx]

@jit(nopython=True)
def calc_AdjMat(coords):
    # np.shape(coords) = (100001, 690, 3)
    # np.shape(AdjMat) = (100001, 690, 690)
    AdjMat = np.zeros((np.shape(coords)[0], 690, 690))
    for i in range(690):
        for j in range(i + 1, 690):
            for t in range(np.shape(coords)[0]):
                d = coords[t, i] - coords[t, j]
                if sum(d**2) < 0.2: # !КВАДРАТ ВЕКТОРА              #########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # OR AdjMat[t, i, j] = AdjMat[t, j, i] = 1/d
                    AdjMat[t, i, j] = AdjMat[t, j, i] = 1
            # стабилизация колебаний по времени tau = 10 шагов MD
            tau = 10                                                #########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for t in range(0, np.shape(coords)[0]-tau, tau):
                if sum(AdjMat[t*tau:(t+1)*tau, i, j]) > tau//3:
                    AdjMat[t*tau:(t+1)*tau, i, j] = np.ones(tau)
    return AdjMat


def calc_my_RandtIndex(Nodes_to_Comms, dt):
    RandtIndex = 0
    if dt < len(Nodes_to_Comms)//2:
        for t in range(np.shape(Nodes_to_Comms)[0]//2):
            cluster_IDs, cluster_sizes = np.unique(Nodes_to_Comms[t], return_counts=True)
            cluster_IDs = cluster_IDs[np.where(cluster_sizes > 1)[0]] ### ID кластеров, размер которых 2 или больше

            comp = 0
            A, B = Nodes_to_Comms[t], Nodes_to_Comms[t + dt]
            for id in cluster_IDs:
                comp += max(np.bincount(B[np.where(A == id)[0]]))/len(np.where(A == id)[0])
            RandtIndex += comp/len(cluster_IDs)  ### смотрим, какая часть из этих кластеров сохранилась

        RandtIndex = RandtIndex/(np.shape(Nodes_to_Comms)[0]//2)
    return RandtIndex


def calc_RandtIndex(Nodes_to_Comms, dt):
    RandtIndex = 0
    if dt < len(Nodes_to_Comms)//2:
        for t in range(np.shape(Nodes_to_Comms)[0]//2):
            cluster_IDs = np.unique(Nodes_to_Comms[t])
            A, B = Nodes_to_Comms[t], Nodes_to_Comms[t + dt]

            comp = 0
            for id in cluster_IDs:
                if max(np.bincount(B[np.where(A == id)[0]])) == len(np.where(A == id)[0]):
                    comp += 1
            RandtIndex += comp/len(cluster_IDs)

        RandtIndex = RandtIndex/(np.shape(Nodes_to_Comms)[0]//2)
    return RandtIndex


def analysis(analysis_type):
    if analysis_type == 'Nodes_to_Comms':
        for tmc in range(0, 400, 5):
            for T in [160, 180]:
                if os.path.exists(f'/home/student/bmim/analysis/clusters/BMI_chains/{T}.{tmc}.npy'):
                    coords = np.load(f'/home/student/bmim/analysis/clusters/BMI_chains/{T}.{tmc}.npy')[0:500] #########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    AdjMat = calc_AdjMat(coords)

                    Nodes_to_Comms = []
                    for t in range(len(coords)):
                        G = nx.from_numpy_matrix(AdjMat[t])
                        partition = community.community_louvain.best_partition(G)
                        Nodes_to_Comms.append(np.array(list(partition.values())))
                    Nodes_to_Comms = np.array(Nodes_to_Comms)
                    np.save(f'/home/student/bmim/analysis/clusters/Nodes_to_Comms/{T}.{tmc}.npy', Nodes_to_Comms)

    if analysis_type == 'my_RandtIndex':
        for tmc in range(0, 400, 5):
            for T in [160, 180]:
                if os.path.exists(f'/home/student/bmim/analysis/clusters/Nodes_to_Comms/{T}.{tmc}.npy'):
                    Nodes_to_Comms = np.load(f'/home/student/bmim/analysis/clusters/Nodes_to_Comms/{T}.{tmc}.npy')
                    dts = np.array([0, 1, 2, 5, 10, 20, 50, 70, 100, 200, 300, 400, 500, 1000, 1500, 2000, 5000, 10000, 20000, 50000])

                    RandtIndex = np.zeros(len(dts))
                    for i in range(len(dts)):
                        RandtIndex[i] = calc_my_RandtIndex(Nodes_to_Comms, dts[i])
                    np.save(f'/home/student/bmim/analysis/clusters/my_RandtIndex/{T}.{tmc}.npy', [dts, RandtIndex])

    if analysis_type == 'RandtIndex':
        for tmc in range(0, 400, 5):
            for T in [160, 180]:
                if os.path.exists(f'/home/student/bmim/analysis/clusters/Nodes_to_Comms/{T}.{tmc}.npy'):
                    Nodes_to_Comms = np.load(f'/home/student/bmim/analysis/clusters/Nodes_to_Comms/{T}.{tmc}.npy')
                    dts = np.array([0, 1, 2, 5, 10, 20, 50, 70, 100, 200, 300, 400, 500, 1000, 1500, 2000, 5000, 10000, 20000, 50000])
                    RandtIndex = np.zeros(len(dts))
                    for i in range(len(dts)):
                        RandtIndex[i] = calc_RandtIndex(Nodes_to_Comms, dts[i])
                    np.save(f'/home/student/bmim/analysis/clusters/RandtIndex/{T}.{tmc}.npy', [dts, RandtIndex])


def average(analysis_type):
    if analysis_type == 'mean':
        for T in [160, 180]:
            mean_size, hist = 0, np.zeros(690)
            cnt = 0
            for tmc in range(0, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/clusters/Nodes_to_Comms/{T}.{tmc}.npy'):
                    Nodes_to_Comms = np.load(f'/home/student/bmim/analysis/clusters/Nodes_to_Comms/{T}.{tmc}.npy')
                    # np.shape(Nodes_to_Comms) = (100001, 690)

                    ### mean cluster size and hist
                    for t in range(len(Nodes_to_Comms)):
                        cluster_IDs, cluster_sizes = np.unique(Nodes_to_Comms[t], return_counts=True)
                        uniqs, cnts = np.unique(cluster_sizes, return_counts=True)
                        for i in range(len(uniqs)):
                            hist[uniqs[i]] += cnts[i]/len(Nodes_to_Comms)
                            mean_size += (cnts[i]*uniqs[i])*uniqs[i]/690/len(Nodes_to_Comms)
                    cnt += 1
            mean_size = mean_size/cnt
            hist = np.array(hist)/cnt

            print(f'mean_size = {mean_size}, {T}K')
            plt.plot(np.arange(len(hist)), hist, label=f'{T}K')
            plt.plot([mean_size, mean_size], [0, 1000], color='black')

        plt.xlim([0, 20]), plt.ylim([0, 50])
        plt.xticks(np.arange(0, 20, 2), np.arange(0, 20, 2).astype(str))
        plt.legend(), plt.show()

    if 'RandtIndex' in analysis_type: ### 'RandtIndex' or 'my_RandtIndex'
        for T in [160, 180]:
            rng, val, cnt = [], [], 0
            cnt = 0
            for tmc in range(0, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/clusters/{analysis_type}/{T}.{tmc}.npy'):
                    if rng == []:
                        rng, val = np.load(f'/home/student/bmim/analysis/clusters/{analysis_type}/{T}.{tmc}.npy')
                    else:
                        val += np.load(f'/home/student/bmim/analysis/clusters/{analysis_type}/{T}.{tmc}.npy')[1]
                    cnt += 1

            val = val/cnt
            plt.plot(rng, val, label=f'{T}K')
        plt.title(f'{analysis_type} decay'), plt.legend(), plt.show()