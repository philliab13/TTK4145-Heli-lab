import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import pandas as pd
import numpy as np
from filter import *





current_dir = os.getcwd()



# Continue trying to get the csv data


numberOfFile=1
# Iterate over all files in the directory

for fileMAT in os.listdir(current_dir):

    if fileMAT.endswith(".mat"): 
        file_path = os.path.join(current_dir, fileMAT)
        # Load the .mat file
        mat_data = scipy.io.loadmat(file_path)
        data=mat_data['ans']
        print(f"Loaded {fileMAT}:")

        rows = [data[i, :].tolist() for i in range(data.shape[0])]

        #Input
        time = np.array(rows[0])
        input_=np.array(rows[1])
        input_opt=np.array(rows[2])
        #States
        lambda_=np.array(rows[3])
        lambda_rate=np.array(rows[4])
        pitch=np.array(rows[5])
        pitch_rate=np.array(rows[6])
        #Inputs optimized
        lambda_opt=np.array(rows[7])
        lambda_rate_opt=np.array(rows[8])
        pitch_opt=np.array(rows[9])
        pitch_rate_opt=np.array(rows[10])

        fs = 500.0       # Hz
        cutoff = 3.0     # desired cutoff frequency in Hz
        order = 6        # filter order

        pitch = butter_lowpass_filter(pitch, cutoff, fs, order)

        fig, axs = plt.subplots(2, 2, figsize=(12, 10))


        axs[0, 0].plot(time, pitch, label='Pitch', linewidth=2)
        axs[0, 0].plot(time, pitch_opt, label='Pitch - Optimized', linewidth=2)

        axs[0, 0].set_title(' Pitch', fontsize=16)
        axs[0, 0].set_ylabel("Pitch [rad]", fontsize=14)
        axs[0,0].set_xlabel("Time [s]", fontsize=14)
        axs[0, 0].legend(fontsize=12)

            # Plot 2: Signal 3 vs Signal 4
        axs[0, 1].plot(time, lambda_, label='Lambda', linewidth=2)
        axs[0, 1].plot(time, lambda_opt, label='Lambda - Optimized', linewidth=2)

        axs[0, 1].set_title('Travel', fontsize=16)
        axs[0, 1].set_ylabel("Pitch rate [rad/sec]", fontsize=14)
        axs[0,1].set_xlabel("Time [s]", fontsize=14)
        axs[0, 1].legend(fontsize=12)
    
        # Plot 3: Signal 5 vs Signal 6
        axs[1, 0].plot(time, pitch_rate, label='Pitch rate', linewidth=2)
        axs[1, 0].plot(time, pitch_rate_opt, label='Pitch rate - Optimized', linewidth=2)

        axs[1, 0].set_title('Pirch Rate', fontsize=16)
        axs[1, 0].set_ylabel("Elevation [rad]", fontsize=14)
        axs[1,0].set_xlabel("Time [s]", fontsize=14)
        axs[1, 0].legend(fontsize=12)

    
        # Plot 4: Signal 7 vs Signal 8
        axs[1, 1].plot(time, lambda_rate, label='Lambda rate', linewidth=2)
        axs[1, 1].plot(time, lambda_rate_opt, label='Lambda rate - Optimized', linewidth=2)

        axs[1, 1].set_title('Travel rate', fontsize=16)
        axs[1, 1].set_ylabel("Elevation rate [rad/sec]", fontsize=14)
        axs[1,1].set_xlabel("Time [s]", fontsize=14)
        axs[1, 1].legend(fontsize=12)


        [ax.grid(True) for ax in axs.ravel()]
        plt.tight_layout(pad=2)
        numberOfFile+=1
                

        plt.savefig(f"{fileMAT}.png")

        ## **SECOND FIGURE: INPUTS**
        plt.figure(figsize=(10, 5))  # Create a new figure
        plt.plot(time, input_, label="Input", linewidth=2)
        plt.plot(time, input_opt, label="Input - Optimized", linewidth=2, linestyle="dashed")

        plt.title("System Inputs", fontsize=16)
        plt.ylabel("Input Signal", fontsize=14)
        plt.xlabel("Time [s]", fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(True)

        plt.savefig(f"{fileMAT}_inputs.png")
        print(f"Saved: {fileMAT}_inputs.png")







