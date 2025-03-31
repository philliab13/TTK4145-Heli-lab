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
        input_pc=np.array(rows[1])
        input_ec=np.array(rows[2])
        input_pc_opt=np.array(rows[3])
        input_ec_opt=np.array(rows[4])
        #States
        lambda_=np.array(rows[5])
        lambda_rate=np.array(rows[6])
        pitch=np.array(rows[7])
        pitch_rate=np.array(rows[8])
        elevation=np.array(rows[9])
        elevation_rate=np.array(rows[10])
        #states optimized
        lambda_opt=np.array(rows[11])
        lambda_rate_opt=np.array(rows[12])
        pitch_opt=np.array(rows[13])
        pitch_rate_opt=np.array(rows[14])
        elevation_opt=np.array(rows[15])
        elevation_rate_opt=np.array(rows[16])

        fs = 500.0       # Hz
        cutoff = 1.5     # desired cutoff frequency in Hz
        order = 6        # filter order

        input_pc = butter_lowpass_filter(input_pc, cutoff, fs, order)
        zeros=np.zeros(3700)
        input_pc=np.concatenate((zeros,input_pc[0:len(input_pc)-3700]))

    

        fig, axs = plt.subplots(3, 2, figsize=(18, 20))
        print(axs)


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

         # Plot 4: Signal 9 vs Signal 10
        axs[2, 0].plot(time, elevation, label='Lambda rate', linewidth=2)
        axs[2, 0].plot(time, elevation_opt, label='Lambda rate - Optimized', linewidth=2)

        axs[2, 0].set_title('Elevation', fontsize=16)
        axs[2, 0].set_ylabel("Elevation [rad]", fontsize=14)
        axs[2,0].set_xlabel("Time [s]", fontsize=14)
        axs[2,0].legend(fontsize=12)

         # Plot 4: Signal 9 vs Signal 10
        axs[2, 1].plot(time, elevation_rate, label='Lambda rate', linewidth=2)
        axs[2, 1].plot(time, elevation_rate_opt, label='Lambda rate - Optimized', linewidth=2)

        axs[2, 1].set_title('Elevation rate', fontsize=16)
        axs[2, 1].set_ylabel("Elevation [rad/sec]", fontsize=14)
        axs[2,1].set_xlabel("Time [s]", fontsize=14)
        axs[2, 1].legend(fontsize=12)


        [ax.grid(True) for ax in axs.ravel()]
        plt.tight_layout(pad=2)
        numberOfFile+=1
                

        plt.savefig(f"{fileMAT}_low.png")

        ## **SECOND FIGURE: INPUTS**
        fig2, axs2 = plt.subplots(1, 2, figsize=(20, 10))
        axs2[0].plot(time, input_pc, label="Input pitch", linewidth=2)
        axs2[0].plot(time, input_pc_opt, label="Input pitch - Optimized", linewidth=2, linestyle="dashed")

        axs2[0].set_title("System Inputs", fontsize=16)
        axs2[0].set_ylabel("Input Signal", fontsize=14)
        axs2[0].set_xlabel("Time [s]", fontsize=14)
        axs2[0].legend(fontsize=12)
        axs2[0].grid(True)

       
        axs2[1].plot(time, input_ec, label="Input elev", linewidth=2)
        axs2[1].plot(time, input_ec_opt, label="Input elev - Optimized", linewidth=2, linestyle="dashed")

        axs2[1].set_title("System Inputs", fontsize=16)
        axs2[1].set_ylabel("Input Signal", fontsize=14)
        axs2[1].set_xlabel("Time [s]", fontsize=14)
        axs2[1].legend(fontsize=12)
        axs2[1].grid(True)

        plt.tight_layout(pad=2)

        plt.savefig(f"{fileMAT}_inputs_low.png")
        print(f"Saved: {fileMAT}_inputs_low.png")







