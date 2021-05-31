import sim
import sys
from main import *
from time import sleep
import pandas as pd
from scipy.io import loadmat


# simRemoteApi.start(26223)

# Read Trajectories
# ======= CSV =======
# with open('./trajectories/test1.csv', 'r', encoding='utf-8') as file:
#     trajectory = pd.read_csv(file)
#     print("Got Trajectory... Resetting")

# ======= MAT =======
trajectory = pd.DataFrame(np.transpose(loadmat('./trajectories/angles_new.mat')['TH']), columns=['th1', 'th2', 'th3',
                                                                                                 'th4', 'th5', 'th6'])
# print(trajectory)

# Initiate connection with CoppeliaSim Server
sim.simxFinish(-1)
clientID = sim.simxStart('127.0.0.1', 26223, True, True, 5000, 5)

if clientID != -1:
    print("Connected to remote API server")

    # Get Joint Handles
    joints = []
    for i in range(6):
        error, handle = sim.simxGetObjectHandle(clientID, 'NiryoOneJoint' + str(i+1), sim.simx_opmode_oneshot_wait)
        if error != 0: raise Exception
        joints.append(handle)
    print(joints)

    # Reset Joints
    for n, joint in enumerate(joints):
        sim.simxSetJointTargetPosition(clientID, joint, trajectory.iloc[0][n], sim.simx_opmode_streaming)
    sleep(3)

    # Iterate Trajectory
    print("Starting motion")
    for i in range(1000):
        sim.simxPauseCommunication(clientID, True)
        for n, joint in enumerate(joints):
            # print(row)
            sim.simxSetJointTargetPosition(clientID, joint, trajectory.iloc[i][n], sim.simx_opmode_streaming)
        sim.simxPauseCommunication(clientID, False)
        sleep(0.006)

else:
    print("Could not connect to remote API server")
    sys.exit("Could not connect")
