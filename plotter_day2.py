import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import pandas as pd
import numpy as np





current_dir = os.getcwd()



# Continue trying to get the csv data


numberOfFile=1
# Iterate over all files in the directory
for fileCSV in os.listdir(current_dir):
    # Check if the file ends with '.mat'
    if fileCSV.endswith(".csv"):
        for fileMAT in os.listdir(current_dir):

            if fileMAT.endswith(".mat") and (fileMAT.split(".")[0]==fileCSV.split(".")[0]): 
                file_path = os.path.join(current_dir, fileMAT)
                # Load the .mat file
                mat_data = scipy.io.loadmat(file_path)
                data=mat_data['elev']
                print(f"Loaded {fileMAT}:")

                rows = [data[i, :].tolist() for i in range(data.shape[0])]

                time = np.array(rows[0])
                lambda_=np.array(rows[1])*(np.pi/180)+np.pi
                lambda_rate=np.array(rows[2])*(np.pi/180)
                pitch=np.array(rows[3])*(np.pi/180)
                pitch_rate=np.array(rows[4])*(np.pi/180)

                ## CSV
                df = pd.read_csv(fileCSV, header=None)
                # print(df.to_string())
                rows_csv=df.to_numpy()
                # print(rows_csv)
                csv_time = df.iloc[:, 0].to_numpy()
                csv_u_k_star=df.iloc[:, 1].to_numpy()
                csv_lambda=df.iloc[:, 2].to_numpy()
                csv_lambda_rate=df.iloc[:, 3].to_numpy()
                csv_pitch=df.iloc[:, 4].to_numpy()
                csv_pitch_rate=df.iloc[:, 5].to_numpy()
            

                fig, axs = plt.subplots(2, 2, figsize=(12, 10))


                axs[0, 0].plot(time, pitch, label='Pitch', linewidth=2)
                axs[0, 0].plot(csv_time, csv_pitch, label='Pitch - Optimized', linewidth=2)

                axs[0, 0].set_title(' Pitch', fontsize=16)
                axs[0, 0].set_ylabel("Pitch [rad]", fontsize=14)
                axs[0,0].set_xlabel("Time [s]", fontsize=14)
                axs[0, 0].legend(fontsize=12)

                    # Plot 2: Signal 3 vs Signal 4
                axs[0, 1].plot(time, lambda_, label='Lambda', linewidth=2)
                axs[0, 1].plot(csv_time, csv_lambda, label='Lambda - Optimized', linewidth=2)

                axs[0, 1].set_title('Travel', fontsize=16)
                axs[0, 1].set_ylabel("Pitch rate [rad/sec]", fontsize=14)
                axs[0,1].set_xlabel("Time [s]", fontsize=14)
                axs[0, 1].legend(fontsize=12)
            
                # Plot 3: Signal 5 vs Signal 6
                axs[1, 0].plot(time, pitch_rate, label='Pitch rate', linewidth=2)
                axs[1, 0].plot(csv_time, csv_pitch_rate, label='Pitch rate - Optimized', linewidth=2)

                axs[1, 0].set_title('Pirch Rate', fontsize=16)
                axs[1, 0].set_ylabel("Elevation [rad]", fontsize=14)
                axs[1,0].set_xlabel("Time [s]", fontsize=14)
                axs[1, 0].legend(fontsize=12)

            
                # Plot 4: Signal 7 vs Signal 8
                axs[1, 1].plot(time, lambda_rate, label='Lambda rate', linewidth=2)
                axs[1, 1].plot(csv_time, csv_lambda_rate, label='Lambda rate - Optimized', linewidth=2)

                axs[1, 1].set_title('Travel rate', fontsize=16)
                axs[1, 1].set_ylabel("Elevation rate [rad/sec]", fontsize=14)
                axs[1,1].set_xlabel("Time [s]", fontsize=14)
                axs[1, 1].legend(fontsize=12)


                [ax.grid(True) for ax in axs.ravel()]
                plt.tight_layout(pad=2)
                numberOfFile+=1
                plt.savefig(f"{fileMAT}.png")

# # Extract each row as a separate list
# rows = [data[i, :].tolist() for i in range(data.shape[0])]


# time = rows[0]
# pitch_encod=rows[1]
# pitch_est=rows[2]
# PitchR_encod=rows[3]
# pitchR_est=rows[4]
# elev_encod=rows[5]
# elev_est=rows[6]
# elevR_encod=rows[7]
# elevR_est=rows[8]



# fig, axs = plt.subplots(2, 2, figsize=(12, 10))


# axs[0,0].plot(time,pitch_encod)
# axs[0,0].plot(time,pitch_est)
# axs[0,0].set_title("Pitch")

# [ax.grid(True) for ax in axs.ravel()]
# plt.savefig(f"{file}.png")





