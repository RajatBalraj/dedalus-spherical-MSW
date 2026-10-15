import subprocess
import sys
import numpy as np

import argparse

data_path = '/home/rajat/Desktop/Project/data/'
all_modules = '/home/rajat/Desktop/Project/codes/full_run_modules/'
dtype = np.float64
####################################################
parser = argparse.ArgumentParser(description="target save folder and spinup IC choice")

parser.add_argument("--run_days", type=int, required=True,
                    help="No. of simulation days to run in days, default = 15 days")
parser.add_argument("--spinup", type=int, required=True, 
                    help="start from spunup up state  or everyting zero except the IC. Should be 0(no spinup) or 1(yes spinup)")
parser.add_argument("--data_folder", type=str,required=True,
                    help="folder to which to save all data. Will also take the initial condition from the checkpoint files in this folder")

parser.add_argument("--IC_type", type=str, required=True,
                    help="random or gaussian initial condition. dlon and lon_center are useless for random option, just put anything")

parser.add_argument("--dlat", type=float, required=True, help="lat extent of gaussian in degrees")
parser.add_argument("--dlon", type=float, required=True, help="lon extent of gaussian in degrees")
parser.add_argument("--lat_center", type=float, default=0, help="lat center of gaussian in degrees")
parser.add_argument("--lon_center", type=float, default=180, help="lon center of gaussian in degrees")

parser.add_argument("--conv_time", type=float, default=2.4, help="convective adjustment timescale in hours")

args = parser.parse_args()


if args.IC_type == 'random':
    print('Did you run the IC file for random separately on one core?')
    from pathlib import Path

    file_path = Path(data_path + f'dlsuhas_reproduce/{args.data_folder}/{args.data_folder}_s1.h5')

    # Check if it is a file and it exists
    if file_path.is_file():
        print("Yes")
    else:
        raise FileNotFoundError('No it does not')
print(args.spinup, type(args.spinup))

if args.spinup==0:
    pass

elif args.spinup==1:
    script = 'spinup.py'
    print(f"--- Starting {script} ---")
    
    result = subprocess.run([sys.executable, all_modules + script,
                                '--data_folder', args.data_folder,
                                '--dlat', str(args.dlat),
                                '--lat_center', str(args.lat_center),
                                '--dlon', str(args.dlon),
                                '--lon_center', str(args.lon_center),
                                '--IC_type', args.IC_type])
        
    # Optional: Stop execution if a script fails (returns a non-zero exit code)
    if result.returncode != 0:
        print(f"Error: {script} failed with exit code {result.returncode}. Stopping.")
        sys.exit(result.returncode)

else:
    raise ValueError('spinup argument should only be 1 (yes) or 0 (no)')

script = 'actual_run.py'
print(f"--- Starting {script} ---")

# sys.executable ensures the script uses the exact same Python environment
result = subprocess.run([sys.executable, all_modules + script,
                         '--data_folder', args.data_folder,
                        '--dlat', str(args.dlat),
                        '--lat_center', str(args.lat_center),
                        '--dlon', str(args.dlon),
                        '--lon_center', str(args.lon_center),
                        '--IC_type', args.IC_type,
                        '--run_days', str(args.run_days),
                        '--spinup', str(args.spinup),
                        '--conv_time', str(args.conv_time)])

# Optional: Stop execution if a script fails (returns a non-zero exit code)
if result.returncode != 0:
    print(f"Error: {script} failed with exit code {result.returncode}. Stopping.")
    sys.exit(result.returncode)

print("All scripts executed successfully!")