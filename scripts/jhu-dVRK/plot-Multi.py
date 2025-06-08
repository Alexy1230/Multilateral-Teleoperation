import matplotlib.pyplot as plt
import matplotlib
import numpy as np

true_pos = np.loadtxt('0607/multi_array_0605.txt', skiprows=1)
expect_pos = np.loadtxt('0607/multi_array_exp_0605.txt', skiprows=1)
error = np.abs(true_pos - expect_pos) * 1000

m1_force = np.loadtxt('0607/multi_m1_force_0605.txt', skiprows=1)
m2_force = np.loadtxt('0607/multi_m2_force_0605.txt', skiprows=1)
puppet_force = np.loadtxt('0607/multi_puppet_force_0605.txt', skiprows=1)
total_force = 0.5 * (m1_force + m2_force)

m1_force_norm = np.linalg.norm(m1_force, axis=1)
m2_force_norm = np.linalg.norm(m2_force, axis=1)
puppet_force_norm = np.linalg.norm(puppet_force, axis=1)
num = len(true_pos)
x = np.arange(num) * 0.002

# Normalized-RMSE for pos
std_pos = np.std(true_pos, axis=0)
rmse_pos = np.sqrt(np.mean((true_pos - expect_pos)**2, axis=0))
nrmse_pos = rmse_pos / std_pos
print(f"rmse of pos: {rmse_pos*1000}")
print(f"nrmse of pos: {nrmse_pos*100}")
print("")

# Normalized-RMSE for force
std_force = np.std(puppet_force, axis=0)
rmse_force = np.sqrt(np.mean((total_force - puppet_force)**2, axis=0))
nrmse_force = rmse_force / std_force
print(f"rmse_force: {rmse_force}")
print(f"nrmse of pos: {nrmse_force*100}")

matplotlib.rcParams['axes.linewidth'] = 1.5
lw = 2


# plot position
fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# X-axis force
axs[0].plot(x, true_pos[:,0]*1000, linewidth=lw, label='true_psm_x')
axs[0].plot(x, expect_pos[:,0]*1000, linewidth=lw, label='expect_psm_x')
axs[0].set_ylabel("Position X (mm)", fontsize=14, fontweight='bold')
axs[0].set_title("Position in X-axis", fontsize=16, fontweight='bold')
axs[0].legend()
axs[0].tick_params(labelsize=14)

# Y-axis force
axs[1].plot(x, true_pos[:,1]*1000, linewidth=lw, label='true_psm_y')
axs[1].plot(x, expect_pos[:,1]*1000, linewidth=lw, label='expect_psm_y')
axs[1].set_ylabel("Position Y (mm)", fontsize=14, fontweight='bold')
axs[1].set_title("Position in Y-axis", fontsize=16, fontweight='bold')
axs[1].legend()
axs[1].tick_params(labelsize=14)

# Z-axis force
axs[2].plot(x, true_pos[:,2]*1000, linewidth=lw, label='true_psm_z')
axs[2].plot(x, expect_pos[:,2]*1000, linewidth=lw, label='expect_psm_z')
axs[2].set_xlabel("Time (s)", fontsize=14, fontweight='bold')
axs[2].set_ylabel("Position Z (mm)", fontsize=14, fontweight='bold')
axs[2].set_title("Position in Z-axis", fontsize=16, fontweight='bold')
axs[2].legend()
axs[2].tick_params(labelsize=14)

fig.tight_layout()



# plot force
fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# X-axis force
axs[0].plot(x, total_force[:,0], linewidth=lw, label='m1+m2_force_x')
axs[0].plot(x, puppet_force[:,0], linewidth=lw, label='puppet_force_x')
axs[0].set_ylabel("Force X (N)", fontsize=14, fontweight='bold')
axs[0].set_title("Force in X-axis", fontsize=16, fontweight='bold')
axs[0].legend()
axs[0].tick_params(labelsize=14)

# Y-axis force
axs[1].plot(x, total_force[:,1], linewidth=lw, label='m1+m2_force_y')
axs[1].plot(x, puppet_force[:,1], linewidth=lw, label='puppet_force_y')
axs[1].set_ylabel("Force Y (N)", fontsize=14, fontweight='bold')
axs[1].set_title("Force in Y-axis", fontsize=16, fontweight='bold')
axs[1].legend()
axs[1].tick_params(labelsize=14)

# Z-axis force
axs[2].plot(x, total_force[:,2], linewidth=lw, label='m1+m2_force_z')
axs[2].plot(x, puppet_force[:,2], linewidth=lw, label='puppet_force_z')
axs[2].set_xlabel("Time (s)", fontsize=14, fontweight='bold')
axs[2].set_ylabel("Force Z (N)", fontsize=14, fontweight='bold')
axs[2].set_title("Force in Z-axis", fontsize=16, fontweight='bold')
axs[2].legend()
axs[2].tick_params(labelsize=14)

fig.tight_layout()
plt.show()




