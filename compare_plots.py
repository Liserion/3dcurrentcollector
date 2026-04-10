import pandas as pd
import matplotlib.pyplot as plt

old = pd.read_csv('/Users/jean/3dcurrentcollector/inputfiles/macro_out.csv')
new = pd.read_csv('/Users/jean/3dcurrentcollector/inputfilesnew/macro_out.csv')

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# cellv vs soc
axes[0].plot(old['soc'], old['cellv'], label='Old mesh')
axes[0].plot(new['soc'], new['cellv'], label='New mesh (real geometry)')
axes[0].set_xlabel('SOC')
axes[0].set_ylabel('Cell Voltage (V)')
axes[0].set_title('Discharge Curve: Voltage vs SOC')
axes[0].legend()
axes[0].grid(True)

# cellv vs time
axes[1].plot(old['time'], old['cellv'], label='Old mesh')
axes[1].plot(new['time'], new['cellv'], label='New mesh (real geometry)')
axes[1].set_xlabel('Time (s)')
axes[1].set_ylabel('Cell Voltage (V)')
axes[1].set_title('Voltage vs Time')
axes[1].legend()
axes[1].grid(True)

# soc vs time
axes[2].plot(old['time'], old['soc'], label='Old mesh')
axes[2].plot(new['time'], new['soc'], label='New mesh (real geometry)')
axes[2].set_xlabel('Time (s)')
axes[2].set_ylabel('SOC')
axes[2].set_title('SOC vs Time')
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.savefig('/Users/jean/3dcurrentcollector/comparison.png', dpi=150)
plt.show()
print("Saved to comparison.png")
