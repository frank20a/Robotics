import sim
import sys
from main import *
from time import sleep
import pandas as pd

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

    # Read Trajectories
    with open('./trajectories/test1.csv', 'r', encoding='utf-8') as file:
        trajectory = pd.from_csv(file)
        print("Got Trajectory... Resetting")

    # Reset Joints
    for n, joint in enumerate(joints):
        sim.simxSetJointTargetPosition(clientID, joint, trajectory.iloc[0][n], sim.simx_opmode_streaming)
    sleep(3)

    # Iterate Trajectory
    print("Starting motion")
    for i in range(120):
        sim.simxPauseCommunication(clientID, True)
        for n, joint in enumerate(joints):
            # print(row)
            sim.simxSetJointTargetPosition(clientID, joint, trajectory.iloc[i][n], sim.simx_opmode_streaming)
        sim.simxPauseCommunication(clientID, False)
        sleep(0.05)

else:
    print("Could not connect to remote API server")
    sys.exit("Could not connect")