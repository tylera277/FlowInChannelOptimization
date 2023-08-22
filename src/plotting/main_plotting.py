# Plotting the output from the main c++ program


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read in the csv data file
raw_velocity_data = pd.read_csv('../../outputData/30jul_u_velocities.csv')

#velocity_data = raw_velocity_data.dropna()
print(raw_velocity_data.head())
velocity_data_numpy = raw_velocity_data.to_numpy()

print(velocity_data_numpy.shape)

# Getting the plot of velocities as a function of position across a slice halfway
# across the problem domain
plt.scatter(np.linspace(0, 1, 33), np.flip(velocity_data_numpy[16,:]) / 1.0)



fig, ax = plt.subplots(1, 1)
velocity_data_numpy = np.rot90(velocity_data_numpy)
ax.contourf(velocity_data_numpy[0:32, 0:32])

#plt.imshow(velocity_data_numpy[0:40, 0:40])

plt.show()

