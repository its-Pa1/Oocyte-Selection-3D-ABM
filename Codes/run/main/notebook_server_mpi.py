import os
import subprocess
import shutil
import re
import time
import glob
import json
import traceback
import faulthandler
import signal
import resource
import numpy as np
import pandas as pd
from mpi4py import MPI
import itertools

# ---------------------------------------------------------------
# === MPI Initialization ===
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Debug helper: `scancel --signal=USR1 <jobid>.<stepid>` dumps the Python stack of every rank to its stderr
faulthandler.register(signal.SIGUSR1, all_threads=True)


def log(msg):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"[R{rank} {time.strftime('%H:%M:%S')} {rss:.0f}MB] {msg}", flush=True)


if rank == 0:
    print("Starting MPI simulation script (same parameter handling as the non-MPI notebook, "
          "reads python_params_baseline.txt)...", flush=True)

original_dir = os.getcwd()
run_directory = os.path.abspath(os.path.join(original_dir, "../server_run"))

param_in_file_path = os.path.join(run_directory, "param.in")
add_param_file_path = os.path.join(run_directory, "additional_param.txt")
GRN_config_file_path = os.path.join(run_directory, "config.txt")
kinetics_file_path = os.path.join(run_directory, "kinetics.txt")
interactions_file_path = os.path.join(run_directory, "interactions.txt")
random_before_signal_file_path = os.path.join(run_directory, "random_before_signal.txt")

somatic_coordinates_dir = run_directory
germ_coordinates_dir = run_directory
os.makedirs(somatic_coordinates_dir, exist_ok=True)
os.makedirs(germ_coordinates_dir, exist_ok=True)

coordinate_files = {
    'sphere': {
        'somatic': os.path.join(somatic_coordinates_dir, "Somatic_coordinates_sphere.txt"),
        'germ': os.path.join(germ_coordinates_dir, "Germ_coordinates_sphere.txt")
    },
    'ellipsoid': {
        'somatic': os.path.join(somatic_coordinates_dir, "Somatic_coordinates_ellipsoid.txt"),
        'germ': os.path.join(germ_coordinates_dir, "Germ_coordinates_ellipsoid.txt")
    },
    'cylinder': {
        'somatic': os.path.join(somatic_coordinates_dir, "Somatic_coordinates_cylinder.txt"),
        'germ': os.path.join(germ_coordinates_dir, "Germ_coordinates_cylinder.txt")
    }
}

quantfile_dir = os.path.join(run_directory, "Quant_Files")
os.makedirs(quantfile_dir, exist_ok=True)


# def build_project():
#     project_root = os.path.abspath("../../")
#     build_dir = os.path.join(project_root, "build")
#     executable_path = os.path.join(build_dir, "bin", "simuscale")
#     #print(f"I come here")
#     if not os.path.exists(build_dir):
#         os.makedirs(build_dir)
#         print(f"Created build directory: {build_dir}")
#     else:
#         print(f"Build directory already exists: {build_dir}")
#     try:
#         print("Running cmake...")
#         subprocess.run(["cmake", ".."], cwd=build_dir, check=True)
#         print("Running make...")
#         subprocess.run(["make"], cwd=build_dir, check=True)
#         print("Build completed successfully.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error during build: {e}", flush=True)
#         MPI.COMM_WORLD.Abort(1)
# # Call the function
# if rank == 0:
#     build_project() # for MPI build it before running the simulation or just the rank 0 will do it


# ---------------------------------------------------------------
def cleanup_previous_runs():
    """Clean up files from previous simulation runs"""
    print("🧹 Cleaning up previous simulation files...")

    files_to_remove = [
        "combined_nurse_cell_activation.txt",
        "mean_summary.txt",
        "cv_summary.txt",
        "timing_summary.txt"
    ]

    for file in files_to_remove:
        file_path = os.path.join(run_directory, file)
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f" ✅ Removed {file}")
            else:
                print(f" ⏭️ Skipped {file} (doesn't exist)")
        except Exception as e:
            print(f" ⚠️ Could not remove {file}: {e}")

    try:
        if os.path.exists(quantfile_dir):
            for filename in os.listdir(quantfile_dir):
                file_path = os.path.join(quantfile_dir, filename)
                try:
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                except Exception as e:
                    print(f" ⚠️ Could not remove {filename}: {e}")
            print(" ✅ Cleaned Quant_Files directory")
        else:
            os.makedirs(quantfile_dir)
            print(" ✅ Created Quant_Files directory")
    except Exception as e:
        print(f" ⚠️ Error with Quant_Files directory: {e}")

    try:
        patterns = [
            "nursecellactivationcounts_id*.txt",
            "nurse_cell_activation_counts_*.txt",
            "combined_part_*.tsv",
            "results_part_*.json",
        ]
        removed_count = 0
        for pattern in patterns:
            for file_path in glob.glob(os.path.join(run_directory, pattern)):
                try:
                    os.remove(file_path)
                    removed_count += 1
                except Exception as e:
                    print(f" ⚠️ Could not remove {os.path.basename(file_path)}: {e}")
        print(f" ✅ Removed {removed_count} leftover result / part files")
    except Exception as e:
        print(f" ⚠️ Error cleaning result files: {e}")

    try:
        cmd = f'find {run_directory} -type d -name "sim_*" -path "*/worker_*/*" -exec rm -rf {{}} +'
        subprocess.run(cmd, shell=True, check=True)
        print(" ✅ Cleaned worker simulation subdirectories using fast shell method")
    except subprocess.CalledProcessError as e:
        print(f" ⚠️ Shell cleanup failed: {e}")

    print("🧹 Cleanup completed!", flush=True)


# ---------------------------------------------------------------
# File updaters (identical to the non-MPI notebook)
# ---------------------------------------------------------------
def update_coordinates(file_path, new_somatic_coordinate_file, new_germ_coordinate_file):
    somatic_filename = f"Somatic_coordinates_{new_somatic_coordinate_file}.txt"
    germ_filename = f"Germ_coordinates_{new_germ_coordinate_file}.txt"
    with open(file_path, 'r') as file:
        lines = file.readlines()
    with open(file_path, 'w') as file:
        for line in lines:
            if 'Somatic_coordinates_' in line:
                line = re.sub(r'Somatic_coordinates_\S+', somatic_filename, line)
            if 'Germ_coordinates_' in line:
                line = re.sub(r'Germ_coordinates_\S+', germ_filename, line)
            file.write(line)


def update_additional_param_file(file_path, factor, PTN_threshold, nurse_cell_activation_file,
                                 alpha, beta, mass_transfer_enabled, reference_nuclei_factor):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(file_path, 'w') as f:
        for line in lines:
            if 'nbdFactor' in line:
                line = re.sub(r'nbdFactor \S+', f'nbdFactor {factor}', line)
            if 'PTN_threshold' in line:
                line = re.sub(r'PTN_threshold \S+', f'PTN_threshold {PTN_threshold}', line)
            if 'nurse_cell_activation_file' in line:
                line = re.sub(r'nurse_cell_activation_file \S+',
                              f'nurse_cell_activation_file {nurse_cell_activation_file}', line)
            if 'alpha' in line and 'alpha' == line.strip().split()[0]:
                line = re.sub(r'alpha \S+', f'alpha {alpha}', line)
            if 'beta' in line and 'beta' == line.strip().split()[0]:
                line = re.sub(r'beta \S+', f'beta {beta}', line)
            if 'mass_transfer_enabled' in line and 'mass_transfer_enabled' == line.strip().split()[0]:
                line = re.sub(r'mass_transfer_enabled \S+',
                              f'mass_transfer_enabled {int(mass_transfer_enabled)}', line)
            if 'reference_nuclei_factor' in line and 'reference_nuclei_factor' == line.strip().split()[0]:
                line = re.sub(r'reference_nuclei_factor \S+',
                              f'reference_nuclei_factor {reference_nuclei_factor}', line)
            f.write(line)


def update_GRN_config_file(file_path, External_Signal_Coeff, Phi):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(file_path, 'w') as f:
        for line in lines:
            if 'External_Signal_Coeff' in line:
                line = re.sub(r'(External_Signal_Coeff\s*:\s*\S+).*',
                              f'External_Signal_Coeff: {External_Signal_Coeff}', line)
            if 'Basal_activity_Phi' in line:
                line = re.sub(r'(Basal_activity_Phi\s*:\s*\S+).*',
                              f'Basal_activity_Phi: {Phi}', line)
            f.write(line)


# NOTE: like in the notebook, these two only change lines that start with the parameter NAME
# (e.g. "d0 0.1"). With value-first files ("0.1  # d0: ...") they leave the file unchanged.
def update_kinetics_file(file_path, d0, d1, a0, a1, a2, d0_2):
    kinetic_params = {
        'd0': d0,
        'd1': d1,
        'a0': a0,
        'a1': a1,
        'a2': a2,
        'd0_2': d0_2,
    }
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(file_path, 'w') as f:
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                f.write(line)
                continue
            parts = stripped.split()
            if len(parts) >= 2 and parts[0] in kinetic_params:
                param_name = parts[0]
                new_val = kinetic_params[param_name]
                line = re.sub(
                    rf'^(\s*{re.escape(param_name)}\s+)\S+',
                    lambda m: f'{m.group(1)}{new_val}',
                    line
                )
            f.write(line)


def update_interactions_file(file_path, i11, i12, i21, i22):
    interaction_params = {
        'interaction_1_1': i11,
        'interaction_1_2': i12,
        'interaction_2_1': i21,
        'interaction_2_2': i22,
    }
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(file_path, 'w') as f:
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                f.write(line)
                continue
            parts = stripped.split()
            if len(parts) >= 2 and parts[0] in interaction_params:
                param_name = parts[0]
                new_val = interaction_params[param_name]
                line = re.sub(
                    rf'^(\s*{re.escape(param_name)}\s+)\S+',
                    lambda m: f'{m.group(1)}{new_val}',
                    line
                )
            f.write(line)


# ---------------------------------------------------------------
def read_last_nurse_cell_activation(filepath):
    """Values at TimeStep 240.0: (nurse cells, total germ cells, pro-oocytes)"""
    try:
        with open(filepath, 'r') as file:
            content = file.read()
        if '\x00' in content:
            print(f"[WARNING] Corrupted file (null bytes): {filepath}")
            return None
        lines = content.splitlines()
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    if float(parts[0]) == 240.0:
                        return float(parts[1]), int(parts[2]), int(parts[3])
                except ValueError:
                    continue
        return None
    except FileNotFoundError:
        return None


def append_to_summary_file(output_file, PTN_threshold, External_Signal_Coeff, alpha, beta, nurse_value):
    with open(output_file, "a") as file:
        file.write(f"{PTN_threshold:.2f}\t{External_Signal_Coeff:.2f}\t{alpha:.2f}\t{beta:.2f}\t{nurse_value:.2f}\n")


# ---------------------------------------------------------------
# Parameter file reader (identical to the notebook)
# ---------------------------------------------------------------
LINSPACE_PARAMS = {
    'PTN_threshold', 'External_Signal_Coeff', 'alpha', 'beta',
    'size_factor', 'reference_nuclei_factor',
    'd0', 'd1', 'a0', 'a1', 'a2', 'd0_2', 'Phi',
    'interaction_1_1', 'interaction_1_2', 'interaction_2_1', 'interaction_2_2',
}


def read_txt_params(txt_file):
    params = {}
    with open(txt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()

                if key in ['shapes', 'sizes']:
                    params[key] = [v.strip() for v in value.split(',')]
                elif key in LINSPACE_PARAMS:
                    parts = [v.strip() for v in value.split(',')]
                    min_val, max_val, points = float(parts[0]), float(parts[1]), int(float(parts[2]))
                    params[key] = np.linspace(min_val, max_val, points)
                elif key in ['total_runs', 'mass_transfer_enabled', 'save_trajectory_files']:
                    params[key] = int(value)
                else:
                    params[key] = value
    return params


def make_identifier(p):
    """Run name, exactly as in the non-MPI notebook."""
    return (
        f"{p['shape']}_{p['size_label']}"
        f"_PTN_{p['PTN_threshold']:.2f}"
        f"_ESC_{p['External_Signal_Coeff']:.2f}"
        f"_A_{p['alpha']:.4f}_B_{p['beta']:.4f}"
        f"_SF_{p['size_factor']:.3f}_RNF_{p['reference_nuclei_factor']:.3f}"
        f"_d0_{p['d0']:.3f}_d1_{p['d1']:.4f}"
        f"_a0_{p['a0']:.3f}_a1_{p['a1']:.2f}_a2_{p['a2']:.3f}"
        f"_d02_{p['d0_2']:.3f}_Phi_{p['Phi']:.3f}"
        f"_i11_{p['interaction_1_1']:.2f}_i12_{p['interaction_1_2']:.3f}"
        f"_i21_{p['interaction_2_1']:.3f}_i22_{p['interaction_2_2']:.2f}"
        f"_Run_{p['simulation_run'] + 1}"
    )


# ---------------------------------------------------------------
def initialize_working_dirs(num_cores):
    """Create the worker directory ONLY for this rank."""
    working_dirs = []
    for i in [rank]:
        work_dir = os.path.join(run_directory, f"worker_{i}")
        os.makedirs(work_dir, exist_ok=True)

        for shape in shapes:
            somatic_file = coordinate_files[shape]['somatic']
            germ_file = coordinate_files[shape]['germ']
            if os.path.exists(somatic_file) and os.path.exists(germ_file):
                shutil.copy(somatic_file, os.path.join(work_dir, os.path.basename(somatic_file)))
                shutil.copy(germ_file, os.path.join(work_dir, os.path.basename(germ_file)))

        for src_path, dest_name in [
            (kinetics_file_path, "kinetics.txt"),
            (interactions_file_path, "interactions.txt"),
            (random_before_signal_file_path, "random_before_signal.txt"),
            (param_in_file_path, "param.in"),
            (add_param_file_path, "additional_param.txt"),
            (GRN_config_file_path, "config.txt"),
        ]:
            if os.path.exists(src_path):
                shutil.copy(src_path, os.path.join(work_dir, dest_name))

        working_dirs.append(work_dir)
    return working_dirs


# ---------------------------------------------------------------
def run_single_simulation(p, work_dir, simulationIdentifier, mass_transfer_flag, save_trajectory_files):
    """One simulation. `p` is one entry of all_parameter_combinations (all values as in the notebook)."""
    original_cwd = os.getcwd()
    sim_start_time = time.time()
    sim_work_dir = None

    try:
        # --- Setup ---
        setup_start = time.time()
        sim_work_dir = os.path.join(work_dir, f"sim_{p['simulation_id']}")
        os.makedirs(sim_work_dir, exist_ok=True)

        required_files = [
            "param.in", "additional_param.txt", "config.txt",
            "kinetics.txt", "interactions.txt", "random_before_signal.txt",
            f"Somatic_coordinates_{p['shape']}.txt", f"Germ_coordinates_{p['shape']}.txt"
        ]
        for filename in required_files:
            src = os.path.join(work_dir, filename)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(sim_work_dir, filename))

        # unique per task; an old file with the same name must never be read or appended to
        nurse_cell_activation_file = os.path.join(run_directory, f"nursecellactivationcounts_id{p['simulation_id']}.txt")
        if os.path.exists(nurse_cell_activation_file):
            os.remove(nurse_cell_activation_file)
        setup_time = time.time() - setup_start

        # --- Parameter updates (same as the notebook) ---
        param_start = time.time()
        local_add_param = os.path.join(sim_work_dir, "additional_param.txt")
        local_GRN_config = os.path.join(sim_work_dir, "config.txt")
        local_param_in = os.path.join(sim_work_dir, "param.in")
        local_kinetics = os.path.join(sim_work_dir, "kinetics.txt")
        local_interactions = os.path.join(sim_work_dir, "interactions.txt")

        update_additional_param_file(local_add_param, p['size_factor'], p['PTN_threshold'],
                                     nurse_cell_activation_file, p['alpha'], p['beta'],
                                     mass_transfer_flag, p['reference_nuclei_factor'])
        update_GRN_config_file(local_GRN_config, p['External_Signal_Coeff'], p['Phi'])
        update_coordinates(local_param_in, p['shape'], p['shape'])
        update_kinetics_file(local_kinetics, p['d0'], p['d1'], p['a0'], p['a1'], p['a2'], p['d0_2'])
        update_interactions_file(local_interactions, p['interaction_1_1'], p['interaction_1_2'],
                                 p['interaction_2_1'], p['interaction_2_2'])
        param_time = time.time() - param_start

        # --- C++ execution ---
        os.chdir(sim_work_dir)
        project_root = os.path.abspath(os.path.join(os.getcwd(), "../../../.."))
        cpp_executable = os.path.join(project_root, "build/bin/simuscale")

        if not os.path.exists(cpp_executable):
            raise FileNotFoundError(f"C++ executable not found: {cpp_executable}")

        cpp_start = time.time()
        try:
            subprocess.run([cpp_executable], capture_output=True, text=True, check=True, timeout=300)
        except subprocess.CalledProcessError as e:
            print(f"[{simulationIdentifier}] FAILED with exit code {e.returncode}", flush=True)
            return (None, None, None, 0, 0, 0, 0, 0)
        except subprocess.TimeoutExpired:
            print(f"[{simulationIdentifier}] Timed out after 300 seconds", flush=True)
            return (None, None, None, 0, 0, 0, 0, 0)
        cpp_time = time.time() - cpp_start

        # --- Output processing ---
        output_start = time.time()
        if save_trajectory_files:
            missing_files = []
            if not os.path.exists("normalization.txt"):
                missing_files.append("normalization.txt")
            if not os.path.exists("trajectory.txt"):
                missing_files.append("trajectory.txt")
            if missing_files:
                print(f"[ERROR] Missing files {missing_files} after {simulationIdentifier}")
                for mf in missing_files:
                    with open(mf, 'w') as f:
                        f.write(f"# Empty {mf} file (original missing)\n")
            os.makedirs(quantfile_dir, exist_ok=True)
            shutil.copy("trajectory.txt", os.path.join(quantfile_dir, f"trajectory_{simulationIdentifier}.txt"))
            shutil.copy("normalization.txt", os.path.join(quantfile_dir, f"normalization_{simulationIdentifier}.txt"))

        os.chdir(original_cwd)
        last_nurse_germ_values = read_last_nurse_cell_activation(nurse_cell_activation_file)

        output_time = time.time() - output_start
        total_time = time.time() - sim_start_time

        return (simulationIdentifier, last_nurse_germ_values, nurse_cell_activation_file,
                total_time, setup_time, param_time, cpp_time, output_time)

    except Exception as e:
        print(f"[ERROR] Simulation {p['simulation_run']} failed: {str(e)}", flush=True)
        traceback.print_exc()
        return (None, None, None, 0, 0, 0, 0, 0)

    finally:
        try:
            os.chdir(original_cwd)
        except Exception:
            pass
        if sim_work_dir:
            shutil.rmtree(sim_work_dir, ignore_errors=True)


# ---------------------------------------------------------------
def append_timeseries(fh, sim_id, path):
    """Append one simulation's time series (long format) to this rank's part file."""
    try:
        if not path or not os.path.exists(path):
            return
        df = pd.read_csv(path, sep="\t", usecols=[0, 1, 2, 3],
                         names=["TimeStep", "NurseCellCount", "TotalGermCount", "ProOocyteCount"],
                         header=0, on_bad_lines="skip")
        df.insert(0, "Simulation", sim_id)
        df.to_csv(fh, sep="\t", index=False, header=False)
    except Exception as e:
        print(f"[Rank {rank}] WARNING: could not append {sim_id}: {e}", flush=True)


def save_timing_summary(output_file, timing_data):
    with open(output_file, "w") as f:
        f.write("=== MPI Simulation Timing Summary ===\n")
        f.write(f"Generated at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for key, value in timing_data.items():
            if isinstance(value, list) and value:
                f.write(f"{key}:\n")
                f.write(f"  Average: {np.mean(value):.3f}s\n")
                f.write(f"  Min: {np.min(value):.3f}s\n")
                f.write(f"  Max: {np.max(value):.3f}s\n")
                f.write(f"  Total: {np.sum(value):.3f}s\n\n")
            else:
                f.write(f"{key}: {value}\n")


def distribute_tasks(total_tasks, num_ranks, rank):
    tasks_per_rank = total_tasks // num_ranks
    extra_tasks = total_tasks % num_ranks
    if rank < extra_tasks:
        start_task = rank * (tasks_per_rank + 1)
        num_tasks = tasks_per_rank + 1
    else:
        start_task = extra_tasks * (tasks_per_rank + 1) + (rank - extra_tasks) * tasks_per_rank
        num_tasks = tasks_per_rank
    return list(range(start_task, start_task + num_tasks))


# ===============================================================
# MAIN
# ===============================================================

if rank == 0:
    cleanup_previous_runs()

comm.Barrier()

if rank == 0:
    job_start_time = time.time()
    print(f"[Rank 0] Starting MPI simulation job at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[Rank 0] Using {size} MPI ranks", flush=True)

# === Load parameters (same file and same logic as the notebook) ===
param_file = os.path.join(run_directory, "python_params_baseline.txt")
sim_params = read_txt_params(param_file)

shapes = sim_params['shapes']
sizes = sim_params['sizes']
PTN_threshold_values = sim_params['PTN_threshold']
External_Signal_Coeff_values = sim_params['External_Signal_Coeff']
alpha_values = sim_params['alpha']
beta_values = sim_params['beta']
size_factor_values = sim_params['size_factor']
reference_nuclei_factor_values = sim_params['reference_nuclei_factor']
total_runs = sim_params['total_runs']
mass_transfer_enabled = sim_params['mass_transfer_enabled']
save_trajectory_files = sim_params['save_trajectory_files']

d0_values = sim_params['d0']
d1_values = sim_params['d1']
a0_values = sim_params['a0']
a1_values = sim_params['a1']
a2_values = sim_params['a2']
d0_2_values = sim_params['d0_2']
Phi_values = sim_params['Phi']

i11_values = sim_params['interaction_1_1']
i12_values = sim_params['interaction_1_2']
i21_values = sim_params['interaction_2_1']
i22_values = sim_params['interaction_2_2']

my_work_dir = initialize_working_dirs(size)[0]
log("worker directory ready")

# === One-at-a-time sweep around the baseline (same as the notebook) ===
ALL_SWEEP_PARAMS = {
    'PTN_threshold': PTN_threshold_values,
    'External_Signal_Coeff': External_Signal_Coeff_values,
    'alpha': alpha_values,
    'beta': beta_values,
    'size_factor': size_factor_values,
    'reference_nuclei_factor': reference_nuclei_factor_values,
    'd0': d0_values,
    'd1': d1_values,
    'a0': a0_values,
    'a1': a1_values,
    'a2': a2_values,
    'd0_2': d0_2_values,
    'Phi': Phi_values,
    'interaction_1_1': i11_values,
    'interaction_1_2': i12_values,
    'interaction_2_1': i21_values,
    'interaction_2_2': i22_values,
}

SWEEP_MODE = "grid"     # "grid" = all combinations; "one_at_a_time" = the notebook's behaviour

baseline = {k: float(v[len(v) // 2]) for k, v in ALL_SWEEP_PARAMS.items()}

if SWEEP_MODE == "grid":
    keys = list(ALL_SWEEP_PARAMS.keys())
    unique_sets = [dict(zip(keys, (float(x) for x in vals)))
                   for vals in itertools.product(*[ALL_SWEEP_PARAMS[k] for k in keys])]
else:
    sweep_sets = [baseline.copy()]
    for param_name, param_values in ALL_SWEEP_PARAMS.items():
        if len(param_values) > 1:
            for val in param_values:
                combo = baseline.copy()
                combo[param_name] = float(val)
                sweep_sets.append(combo)
    unique_sets = []
    seen = set()
    for s in sweep_sets:
        key = tuple(round(s[k], 8) for k in sorted(ALL_SWEEP_PARAMS.keys()))
        if key not in seen:
            seen.add(key)
            unique_sets.append(s)
# === All simulations (identical on every rank) ===
all_parameter_combinations = []
simulation_id = 0

for shape in shapes:
    for combo in unique_sets:
        for simulation_run in range(total_runs):
            all_parameter_combinations.append({
                'simulation_id': simulation_id,
                'simulation_run': simulation_run,
                'shape': shape,
                'size_label': 's',
                **combo,
            })
            simulation_id += 1

total_simulations = len(all_parameter_combinations)
all_names = [make_identifier(c) for c in all_parameter_combinations]
if len(set(all_names)) != len(all_names):
    if rank == 0:
        print("ERROR: run names are not unique (two sweep values print the same). "
              "Increase the decimals in make_identifier().", flush=True)
    comm.Abort(1)

if rank == 0:
    print(f"[Rank 0] Total simulations: {total_simulations}", flush=True)
    id_to_params = dict(zip(all_names, all_parameter_combinations))

    combined_data_file = os.path.join(run_directory, "combined_nurse_cell_activation.txt")
    mean_summary_file = os.path.join(run_directory, "mean_summary.txt")
    cv_summary_file = os.path.join(run_directory, "cv_summary.txt")
    germ_summary_file = os.path.join(run_directory, "germ_cell_summary.txt")
    all_combinations_file = os.path.join(run_directory, "all_combinations_results.txt")

    with open(mean_summary_file, "w") as f:
        f.write("PTN_threshold\tExternal_Signal_Coeff\talpha\tbeta\tMean_Last_NurseCellActivation\n")
    with open(cv_summary_file, "w") as f:
        f.write("PTN_threshold\tExternal_Signal_Coeff\talpha\tbeta\tCV_Last_NurseCellActivation\n")

    with open(germ_summary_file, "w") as f:
        headers = [
            "[P_n]_0", "zeta", "alpha", "beta",
            "lambda", "V_nucleus",
            "d0_o", "d1", "a0", "a1", "a2", "d0_n", "Phi",
            "theta11", "theta12", "theta21", "theta22",
            "Mean_Nurse", "CV_Nurse",
            "Mean_Germ", "CV_Germ",
            "Mean_Oocyte", "CV_Oocyte",
            "N/G_Ratio",
        ]
        f.write('\t'.join(h.ljust(10) for h in headers) + '\n')

    with open(all_combinations_file, "w") as f:
        f.write(
            "SimulationID\tShape\tSize\t"
            "P_n\tzeta\talpha\tbeta\t"
            "lambda\tV_nucleus\t"
            "d0_o\td1\ta0\ta1\ta2\td0_n\tPhi\t"
            "theta_11\ttheta_12\ttheta_21\ttheta_22\t"
            "Run\tNurseCells\tTotalGermCells\tOocytes\n"
        )

comm.Barrier()

my_tasks = list(range(rank, total_simulations, size))

if rank == 0:
    print(f"[Rank 0] Processing approximately {total_simulations // size} tasks per rank", flush=True)
print(f"[Rank {rank}] Assigned {len(my_tasks)} tasks: {my_tasks[:3]}{'...' if len(my_tasks) > 3 else ''}", flush=True)

# === Execute assigned simulations; stream each time series to this rank's part file ===
local_results = []
rank_sim_start = time.time()
part_file = os.path.join(run_directory, f"combined_part_{rank:04d}.tsv")

consec_fail = 0

with open(part_file, "w") as fpart:
    for task_idx in my_tasks:
        p = all_parameter_combinations[task_idx]
        simulationIdentifier = all_names[task_idx]

        result = run_single_simulation(p, my_work_dir, simulationIdentifier,
                                       mass_transfer_enabled, save_trajectory_files)

        if result[0] is None:
            consec_fail += 1
            if consec_fail >= 5:
                print(f"[Rank {rank}] 5 consecutive failures, aborting the job", flush=True)
                comm.Abort(1)
        else:
            consec_fail = 0

        # same rule as the notebook: only runs with a valid TimeStep-240 row enter the combined file
        if (result[0] is not None and result[1] is not None and result[1][0] is not None and result[2]):
            append_timeseries(fpart, result[0], result[2])

        local_results.append(result)
        print(f"[Rank {rank}] Completed task {len(local_results)}/{len(my_tasks)} - {simulationIdentifier}", flush=True)

rank_sim_time = time.time() - rank_sim_start
print(f"[Rank {rank}] Rank simulations ({len(my_tasks)} tasks) completed in {rank_sim_time:.2f} seconds", flush=True)

shutil.rmtree(my_work_dir, ignore_errors=True)    # Fix 2: once per rank, after the loop
# === Hand results to rank 0 through files (no MPI messages, so nothing can hang here) ===
res_file = os.path.join(run_directory, f"results_part_{rank:04d}.json")
with open(res_file + ".tmp", "w") as f:
    json.dump(local_results, f)
os.replace(res_file + ".tmp", res_file)
log("results file written")

gather_start = time.time()
if rank == 0:
    all_results = []
    last_report = time.time()
    for r in range(size):
        rp = os.path.join(run_directory, f"results_part_{r:04d}.json")
        while not os.path.exists(rp):
            time.sleep(10)
            os.listdir(run_directory)
            if time.time() - last_report > 600:
                print(f"[Rank 0] still waiting for results of rank {r}", flush=True)
                last_report = time.time()
        with open(rp) as f:
            all_results.extend(json.load(f))
        os.remove(rp)
    all_results = [item for item in all_results if item[0] is not None]
gather_time = time.time() - gather_start

if rank == 0:
    processing_start = time.time()
    print(f"Processing {len(all_results)} results...", flush=True)

    param_groups = {}
    total_times, setup_times, param_times, cpp_times, output_times = [], [], [], [], []

    for result in all_results:
        (simulationIdentifier, last_nurse_germ_values, nurse_cell_file,
         total_time, setup_time, param_time, cpp_time, output_time) = result

        total_times.append(total_time)
        setup_times.append(setup_time)
        param_times.append(param_time)
        cpp_times.append(cpp_time)
        output_times.append(output_time)

        if last_nurse_germ_values is not None and last_nurse_germ_values[0] is not None:
            nurse_count, total_germ_count, pro_oocyte_count = last_nurse_germ_values

            matched = id_to_params.get(simulationIdentifier, None)
            if matched:
                m = matched
                param_key = tuple(round(m[k], 8) for k in sorted(ALL_SWEEP_PARAMS.keys()))
                param_key = (m['shape'], m['size_label']) + param_key

                if param_key not in param_groups:
                    param_groups[param_key] = {
                        'params': m,
                        'nurse_counts': [],
                        'germ_counts': [],
                        'pro_oocyte_counts': [],
                        'nurse_cell_files': [],
                    }

                param_groups[param_key]['nurse_counts'].append(nurse_count)
                param_groups[param_key]['germ_counts'].append(total_germ_count)
                param_groups[param_key]['pro_oocyte_counts'].append(pro_oocyte_count)
                param_groups[param_key]['nurse_cell_files'].append(nurse_cell_file)

                run_num = simulationIdentifier.split("_Run_")[-1]
                with open(all_combinations_file, "a") as f:
                    f.write(
                        f"{simulationIdentifier}\t{m['shape']}\t{m['size_label']}\t"
                        f"{m['PTN_threshold']}\t{m['External_Signal_Coeff']}\t"
                        f"{m['alpha']}\t{m['beta']}\t"
                        f"{m['size_factor']}\t{m['reference_nuclei_factor']}\t"
                        f"{m['d0']}\t{m['d1']}\t{m['a0']}\t{m['a1']}\t{m['a2']}\t"
                        f"{m['d0_2']}\t{m['Phi']}\t"
                        f"{m['interaction_1_1']}\t{m['interaction_1_2']}\t"
                        f"{m['interaction_2_1']}\t{m['interaction_2_2']}\t"
                        f"{run_num}\t{nurse_count:.4f}\t{total_germ_count}\t{pro_oocyte_count}\n"
                    )

    log("results grouped")

    # === Combined time-series file: same wide layout as the notebook ===
    t0 = time.time()
    cols4 = ["TimeStep", "NurseCellCount", "TotalGermCount", "ProOocyteCount"]
    sim_arrays = {}
    for r in range(size):
        pp = os.path.join(run_directory, f"combined_part_{r:04d}.tsv")
        if os.path.exists(pp) and os.path.getsize(pp) > 0:
            df = pd.read_csv(pp, sep="\t", header=None, names=["Simulation"] + cols4)
            for sim, g in df.groupby("Simulation", sort=False):
                sim_arrays[sim] = g[cols4].to_numpy()
        if os.path.exists(pp):
            os.remove(pp)

    print(f"Creating combined data file with {len(sim_arrays)} simulations...", flush=True)
    if sim_arrays:
        ref = sim_arrays[next(iter(sim_arrays))][:, 0]
        n = len(ref)
        ids = sorted(sim_arrays)
        arr = np.full((n, 3 * len(ids)), np.nan)
        for i, s in enumerate(ids):
            a = sim_arrays[s]
            m_ = min(n, len(a))
            arr[:m_, 3 * i:3 * i + 3] = a[:m_, 1:4]
        names = ["TimeStep"] + [f"{s}{x}" for s in ids for x in ("_Nurse", "_Germ", "_ProOocyte")]
        combined_df = pd.DataFrame(np.column_stack([ref, arr]), columns=names)
        combined_df.to_csv(combined_data_file, sep="\t", index=False)
        # combined_df.to_csv(combined_data_file, sep="\t", index=False)
        print(f"Combined data saved with {len(combined_df.columns) - 1} columns", flush=True)
        del combined_df, arr
    del sim_arrays
    log(f"combined file step took {time.time() - t0:.1f}s")

    # === Aggregated summaries (same as the notebook) ===
    for param_key, data in param_groups.items():
        m = data['params']
        nurse_list = data['nurse_counts']
        germ_list = data['germ_counts']
        oocyte_list = data['pro_oocyte_counts']
        nurse_files = data['nurse_cell_files']

        if nurse_list:
            mean_nurse = np.mean(nurse_list)
            cv_nurse = np.std(nurse_list) / mean_nurse if mean_nurse != 0 else 0
            mean_germ = np.mean(germ_list) if germ_list else 0
            cv_germ = np.std(germ_list) / mean_germ if mean_germ != 0 else 0
            mean_oocyte = np.mean(oocyte_list) if oocyte_list else 0
            cv_oocyte = np.std(oocyte_list) / mean_oocyte if mean_oocyte != 0 else 0
            ratio = mean_nurse / mean_germ if mean_germ != 0 else 0

            PTN = m['PTN_threshold']
            ESC = m['External_Signal_Coeff']
            A = m['alpha']
            B = m['beta']

            append_to_summary_file(mean_summary_file, PTN, ESC, A, B, mean_nurse)
            append_to_summary_file(cv_summary_file, PTN, ESC, A, B, cv_nurse)
            with open(germ_summary_file, "a") as f:
                vals = [
                    f"{m['PTN_threshold']:.2f}",
                    f"{m['External_Signal_Coeff']:.2f}",
                    f"{m['alpha']:.4f}",
                    f"{m['beta']:.4f}",
                    f"{m['size_factor']:.3f}",
                    f"{m['reference_nuclei_factor']:.3f}",
                    f"{m['d0']:.3f}",
                    f"{m['d1']:.3f}",
                    f"{m['a0']:.3f}",
                    f"{m['a1']:.3f}",
                    f"{m['a2']:.3f}",
                    f"{m['d0_2']:.3f}",
                    f"{m['Phi']:.3f}",
                    f"{m['interaction_1_1']:.3f}",
                    f"{m['interaction_1_2']:.3f}",
                    f"{m['interaction_2_1']:.3f}",
                    f"{m['interaction_2_2']:.3f}",
                    f"{mean_nurse:.3f}",
                    f"{cv_nurse:.3f}",
                    f"{mean_germ:.3f}",
                    f"{cv_germ:.3f}",
                    f"{mean_oocyte:.3f}",
                    f"{cv_oocyte:.3f}",
                    f"{ratio:.3f}",
                ]
                f.write('\t'.join(v.ljust(10) for v in vals) + '\n')

            print(f"  PTN={PTN:.2f} ESC={ESC:.2f} A={A:.4f} B={B:.4f} "
                  f"SF={m['size_factor']:.3f} RNF={m['reference_nuclei_factor']:.3f} "
                  f"a0={m['a0']:.3f} Phi={m['Phi']:.3f} → "
                  f"Nurse={mean_nurse:.2f} Germ={mean_germ:.2f} "
                  f"ProOocyte={mean_oocyte:.2f} ({len(nurse_list)} runs)")

        for nf in nurse_files:
            if nf and os.path.exists(nf):
                try:
                    os.remove(nf)
                except OSError:
                    pass

    processing_time = time.time() - processing_start

    print(f"[Rank 0] Data gathering time: {gather_time:.3f}s")
    print(f"[Rank 0] Result processing time: {processing_time:.3f}s")

    if total_times:
        print("\nAverage simulation times:")
        print(f"  Total:             {np.mean(total_times):.3f}s")
        print(f"  Setup:             {np.mean(setup_times):.3f}s")
        print(f"  Parameters:        {np.mean(param_times):.3f}s")
        print(f"  C++ execution:     {np.mean(cpp_times):.3f}s")
        print(f"  Output processing: {np.mean(output_times):.3f}s")

    total_job_time = time.time() - job_start_time
    print("\nAll simulations complete.")
    print(f"Total job time: {total_job_time:.2f}s ({total_job_time / 60:.2f} min)")
    print(f"Job ended at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\n📄 All individual results → {all_combinations_file}", flush=True)

    timing_file = os.path.join(run_directory, "timing_summary.txt")
    save_timing_summary(timing_file, {
        "Total Job Time": total_job_time,
        "Number of MPI Ranks": size,
        "Total Simulations": len(all_results),
        "Setup Times": setup_times,
        "Parameter Update Times": param_times,
        "C++ Execution Times": cpp_times,
        "Output Processing Times": output_times,
    })
    print(f"Timing summary saved to: {timing_file}")

    try:
        cmd = f'find {run_directory} -maxdepth 1 -type d -name "worker_*" -exec rm -rf {{}} +'
        subprocess.run(cmd, shell=True, check=True)
        print("Removed all worker directories after simulation completion", flush=True)
    except subprocess.CalledProcessError as e:
        print(f"Worker directory cleanup failed: {e}")
