# run.py (updated)

import subprocess
from data_get import get_2d_data, get_3d_data
from data_eval import eval_2d_data, eval_3d_data
import config  # Import to access config.day

print("Retrieving data for HYCOM 25-degree 2d data")
# get_2d_data()  # Uncomment if needed

print("Retrieving data for HYCOM 25-degree 3d data")
# get_3d_data()

# eval_2d_data()  # Uncomment if needed
# eval_3d_data()  # Uncomment if needed

# List of all plt_ scripts to run (add/remove as needed)
plot_scripts = [
    'plt_salinx_vert.py',
    'plt_salin_shaded.py',
    'plt_salinx_bands.py',
    'plt_test.py',
    'plt_salin.py',
    'plt_salin_bands.py',
    'plt_ssh.py',
    'plt_salinx.py',
]

# Run each plot script, passing config.day as argument
print(f"\nGenerating plots for day: {config.day}")
for script in plot_scripts:
    try:
        result = subprocess.run(
            ['python', script, config.day], check=True, capture_output=True, text=True)
        print(f"Successfully ran {script}: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {script}: {e.stderr.strip()}")
