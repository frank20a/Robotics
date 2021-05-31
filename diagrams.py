from main import solveDirect, parseDH2TransMatrix
from data import *
from matplotlib import pyplot as plt
import pandas as pd
from scipy.io import loadmat
import numpy as np
from numpy import pi, sqrt

traj = pd.DataFrame(np.transpose(loadmat('./trajectories/angles_new.mat')['TH']), columns=['th1', 'th2', 'th3',
                                                                                           'th4', 'th5', 'th6'])
t = np.linspace(0, 6, 1000)
Ts = 0.006

M = parseDH2TransMatrix(DH_Table)


def joint_angles():
    fig, axs = plt.subplots(2, 3)

    fig.suptitle("Joint angles", fontsize=25)

    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            ax.plot(t, traj.iloc[:, 3*i+j] * 180 / pi)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Angle (deg)")
            ax.set_title(r"$\theta$" + str(3*i+j+1))
            ax.grid(True)

    plt.show()


def joint_angular_velocity():
    fig, axs = plt.subplots(2, 3)

    fig.suptitle("Joint angular velocity", fontsize=25)
    fig.tight_layout(h_pad=1.8)

    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            ax.plot([0, 6], [0, 0], 'r')
            ax.plot(t[1:], np.diff(traj.iloc[:, 3*i+j] * (180 / pi)) / Ts)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Angular Velocity (deg/sec)")
            ax.set_title(r"$\dot{\theta}$" + str(3*i+j+1))
            ax.grid(True)

    plt.show()


def cartesian_positions():
    xyz = []
    for i in range(1000):
        tmp = traj.iloc[i]
        xyz.append(np.matmul(solveDirect(M, (tmp.iloc[0], tmp.iloc[1], tmp.iloc[2], tmp.iloc[3], tmp.iloc[4],
                                             tmp.iloc[5])), [[0], [0], [0], [1]]))

    x = []
    y = []
    z = []
    for pos in xyz:
        x.append(pos[0][0] * 1e2)
        y.append(pos[1][0] * 1e2)
        z.append(pos[2][0] * 1e2)

    fig = plt.figure()
    fig.suptitle("Gripper position", fontsize=25)
    
    ax = fig.add_subplot(projection='3d')
    ax.plot3D(x[::5], y[::5], z[::5])
    ax.set_xlabel('X (cm)')
    ax.set_ylabel('Y (cm)')
    ax.set_zlabel('Z (cm)')

    plt.show()


def cartesian_velocities(split=False):
    xyz = []
    for i in range(1000):
        tmp = traj.iloc[i]
        xyz.append(np.matmul(solveDirect(M, (tmp.iloc[0], tmp.iloc[1], tmp.iloc[2], tmp.iloc[3], tmp.iloc[4],
                                             tmp.iloc[5])), [[0], [0], [0], [1]]))

    x = []
    y = []
    z = []
    for pos in xyz:
        x.append(pos[0][0] * 1e2)
        y.append(pos[1][0] * 1e2)
        z.append(pos[2][0] * 1e2)

    x = np.array(np.diff(x) / Ts, dtype='float32')
    y = np.array(np.diff(y) / Ts, dtype='float32')
    z = np.array(np.diff(z) / Ts, dtype='float32')

    if split:
        fig = plt.figure()
        fig.suptitle("Gripper velocity", fontsize=25)

        ax = fig.add_subplot(projection='3d')
        ax.plot3D(x[::5], y[::5], z[::5])
        ax.set_xlabel(r'$\dot{x}$ (cm/s)')
        ax.set_ylabel(r'$\dot{y}$ (cm/s)')
        ax.set_zlabel(r'$\dot{z}$ (cm/s)')
        ax.grid(True)
    else:
        plt.plot(t[1:], x)
        plt.plot(t[1:], y)
        plt.plot(t[1:], z)
        plt.plot(t[1:], sqrt(x**2 + y**2 + z**2), 'black')
        plt.legend([r'$\dot{x}$', r'$\dot{y}$', r'$\dot{z}$', r'$\dot{v}=\sqrt{\dot{x}^2+\dot{y}^2+\dot{z}^2}$'],
                   fontsize=18)
        plt.ylabel('Velocity (cm/s)', fontsize=15)
        plt.xlabel('Time (s)', fontsize=15)
        plt.title('Gripper Velocity per axis', fontsize=25)
        plt.grid(True)

    plt.show()


if __name__ == '__main__':
    cartesian_velocities(True)
