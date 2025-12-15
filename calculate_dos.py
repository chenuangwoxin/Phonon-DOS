# -*- coding: utf-8 -*-
"""
VDOS Calculator 
Calculates Vibrational Density of States (VDOS) from Molecular Dynamics velocity files.
The velocity files is generated from LAMMPS dump command with:
    dump velocity all custom N file_name id type vx vy vz
    dump_modify velocity sort id

@Author: Jingchao Zhang (Tongji University) 
Date: 2025-12-12
"""

import numpy as np
from scipy import signal 
import os
import time
import matplotlib.pyplot as plt
import sys

class VDOSCalculator:
    
    def __init__(self, n_atoms, dt_ps, interval, total_frames, 
                 fft_window_size, overlap_shift, smooth_width, output_dir):
        """
        Initialize the VDOS Calculator.
        
        Parameters:
        n_atoms        : Number of atoms
        dt_ps          : MD simulation time step (unit: ps)
        interval       : Dump interval (steps between frames)
        total_frames   : Total number of frames to read
        fft_window_size: Size of the time window for FFT (number of frames)
        overlap_shift  : Shift size for the sliding window (number of frames)
        smooth_width   : Width factor for smoothing (typically 1-5)
        output_dir     : Directory to save results
        """
        self.n_atoms = n_atoms
        self.frames_read = total_frames
        
        # Time related parameters
        self.dt_frame = dt_ps * interval  # Time between two frames (ps)
        self.fft_window_size = fft_window_size
        self.overlap_shift = overlap_shift
        
        # Calculate total chunks (Sliding Windows)
        if self.overlap_shift > 0:
            self.n_chunks = int((self.frames_read - self.fft_window_size) / self.overlap_shift) + 1
        else:
            self.n_chunks = 1
            
        self.smooth_width = smooth_width
        self.output_dir = output_dir
        
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        self.log_file = open(os.path.join(self.output_dir, 'log.txt'), 'w')

    def log(self, message):
        """Helper function for logging to console and file."""
        print(message)
        print(message, file=self.log_file, flush=True)

    def close(self):
        """Close the log file handler."""
        if self.log_file:
            self.log_file.close()

    def _compactVels(self, file_Vels):

        from subprocess import call
        command = ["compress_file", file_Vels]
        print("\nRunning the " + '\"' + " ".join(command) + '\"' + ' command to generate the compact_velocity file.')
        call(command)

    def load_mass(self,input_file, type_mass_map:dict):
        mass = np.zeros(self.n_atoms*3)

        with open(input_file) as f:
            for _ in range(9):
                f.readline()
            for i in range(self.n_atoms):
                line = f.readline()
                parts = line.strip().split()
                atom_type = int(parts[1])
                mass_value = type_mass_map[atom_type]
                mass[i*3:i*3+3] = mass_value
        return mass

    def load_velocity(self, input_file):
        '''
        compress the velocity file to binary format for quick reading and processing
        '''
        self.log(f"Reading file: {input_file} ...")
        # Initialize matrix: [3 * N_atoms, Total_Frames]
        data_size = self.n_atoms * 3
        velocity = np.zeros((data_size, self.frames_read))       
        '''
        """
        Reads velocity data from file.
        Supported formats per line: [ID Type Vx Vy Vz] OR [Vx Vy Vz]
        """
        if format == 'lammps_dump':
            try:
                with open(input_file, 'r') as fid:
                    
                    for i in range(self.frames_read):
                        for k in range(self.n_atoms):
                            line = fid.readline()
                            if not line: break
                            
                            parts = list(map(float, line.strip().split()))
                            
                            # Handle different column formats
                            if len(parts) >= 5:
                                # Assuming: ID Type Vx Vy Vz
                                vx, vy, vz = parts[2], parts[3], parts[4]
                            elif len(parts) == 3:
                                # Assuming: Vx Vy Vz
                                vx, vy, vz = parts[0], parts[1], parts[2]
                            else:
                                raise ValueError(f"Line format error: frame {i}, atom {k}")

                            velocity[k*3+0, i] = vx
                            velocity[k*3+1, i] = vy
                            velocity[k*3+2, i] = vz
                            
                        if i % 5000 == 0 and i > 0:
                            self.log(f"  Read {i}/{self.frames_read} frames")

                    self.log(f"Finished reading: {input_file}")
                    return velocity
                
            except FileNotFoundError:
                self.log(f"Error: File {input_file} not found")
                sys.exit(1)
                '''
        bin_files = {'vx': 'vx.bin', 'vy': 'vy.bin', 'vz': 'vz.bin'}

        all_bin_files_exist = all(os.path.exists(bin_files[key]) for key in bin_files)
        if not all_bin_files_exist:
            self.log("Binary files not found, compressing velocity data...")
            self._compactVels(input_file)
        else:
            self.log("Binary files found, reading directly...")

        try:
            raw_vx = np.fromfile(bin_files['vx'], dtype=np.float64).reshape(-1, self.n_atoms)
            raw_vy = np.fromfile(bin_files['vy'], dtype=np.float64).reshape(-1, self.n_atoms)
            raw_vz = np.fromfile(bin_files['vz'], dtype=np.float64).reshape(-1, self.n_atoms)
        
            velocity[0::3, :] = raw_vx[:self.frames_read, :].T
            velocity[1::3, :] = raw_vy[:self.frames_read, :].T
            velocity[2::3, :] = raw_vz[:self.frames_read, :].T
        
        except FileNotFoundError:
            self.log("Error: One or more binary files not found")
            sys.exit(1)
        self.log(f"Finished reading binary compressed files.")
        return velocity

    def compute_vdos(self, velocity, mass):
        """
        Core VDOS logic:
        1. FFT transform
        2. Power Spectral Density calculation
        3. Normalization
        """
        self.log("Starting VDOS calculation (FFT)...")

        # Number of frequency points (rfft output is N/2 + 1)
        n_freq_points = int(self.fft_window_size / 2 + 1)
        
        # VDOS Matrix: [Freq Points, Chunks + Average Column]
        vdos_matrix = np.zeros((n_freq_points, self.n_chunks + 1))

        # Calculate frequency axis (Unit: THz, assuming time is ps)
        frequency = np.fft.rfftfreq(self.fft_window_size, d=self.dt_frame)
        df = frequency[1] # Frequency resolution

        for i in range(self.n_chunks):
            # Determine sliding window range
            start_idx = i * self.overlap_shift
            end_idx = start_idx + self.fft_window_size
            
            if end_idx > velocity.shape[1]:
                break

            # Extract current window data: [3*Natoms, Window_Size]
            vel_chunk = velocity[:, start_idx:end_idx]
            
            # 1. FFT Transform (along time axis, axis=1)
            vel_fft = np.fft.rfft(vel_chunk, axis=1) * np.sqrt(mass[:, None])
            
            # 2. Compute Power Spectral Density (PSD)
            # Sum over all atoms (axis 0)
            power_spectrum = np.real(np.sum(vel_fft * np.conj(vel_fft), axis=0))
            
            # 3. Normalization (Area under curve = 1)
            integral = np.sum(power_spectrum) * df
            
            if integral > 0:
                vdos_matrix[:, i] = power_spectrum / integral
            
            if i % 10 == 0 or i == self.n_chunks - 1:
                self.log(f"  Calculating chunk {i+1}/{self.n_chunks}")

        # Calculate average over all chunks (Store in the last column)
        vdos_matrix[:, self.n_chunks] = np.average(vdos_matrix[:, 0:self.n_chunks], axis=1)
        
        return frequency, vdos_matrix

    def smooth_and_save(self, frequency, vdos_matrix, file_tag):
        """
        Apply smoothing filters and save data to text file.
        """
        self.log("Smoothing and saving data...")
        
        freq_len = len(frequency)
        avg_col_idx = self.n_chunks
        
        # Raw Average VDOS
        raw_vdos = vdos_matrix[:, avg_col_idx]
        
        # === Smoothing Parameter Calculation ===
        df = frequency[1]
        if df == 0: df = 1e-10
        
        # Calculate window points based on widthwin
        window_pts = int(np.ceil(self.smooth_width / df / 10))
        if window_pts < 5: window_pts = 5
        if window_pts % 2 == 0: window_pts += 1 # Must be odd
        
        self.log(f"  Smoothing window points: {window_pts}")

        # === 1. Rectangle Window Smoothing (Moving Average) ===
        rect_window = np.ones(window_pts) / window_pts
        smooth_rect = np.convolve(raw_vdos, rect_window, 'same')
        
        # === 2. Savitzky-Golay Filter ===
        # Window length must not exceed data length
        sg_window = min(window_pts, freq_len - 2)
        if sg_window % 2 == 0: sg_window -= 1
        
        # Ensure sg_window is at least slightly larger than polyorder (2)
        if sg_window < 5: sg_window = 5
        
        smooth_sg = signal.savgol_filter(raw_vdos, sg_window, 2)

        # === Save Data ===
        output_data = np.column_stack((frequency, raw_vdos, smooth_rect, smooth_sg))
        
        out_name = os.path.join(self.output_dir, f"{file_tag}_VDOS.txt")
        header = 'Frequency(THz)  Raw_Ave_VDOS  Rect_Smooth  Sav_Gol_Smooth'
        np.savetxt(out_name, output_data, fmt='%15.9f', header=header)
        
        return output_data

    def plot_result(self, data,cutoff_freq, file_tag):
        """Plot the VDOS results."""
        self.log("Plotting results...")
        
        freq = data[:, 0]
        raw = data[:, 1]
        smooth_sg = data[:, 3]

        plt.figure(figsize=(10, 8), dpi=300)
        
        # Style settings
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman']
        plt.tick_params(direction='in', width=2, length=6, labelsize=14)
        
        plt.plot(freq, raw, 'k-', lw=1, alpha=0.2, label='Raw Data')
        plt.plot(freq, smooth_sg, 'r-', lw=2.5, label='Savitzky-Golay Smooth')

        plt.xlim(0, cutoff_freq)
        
        plt.title(f'VDOS - {file_tag}', fontsize=18)
        plt.xlabel('Frequency (THz)', fontsize=18)
        plt.ylabel('VDOS (Probability Density)', fontsize=18)
        plt.legend(fontsize=14)
        
        # Optional: Limit X axis (e.g., 0-20 THz)
        # plt.xlim(0, 20)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, f"VDOS_{file_tag}.png"))
        plt.close()


if __name__ == '__main__':
    # ================= PARAMETERS =================
    # Physical Parameters
    NUMBER_STEPS = 1e6        # Total steps to read
    DT_PS        = 0.0005        # Time step (ps)
    INTERVAL     = 20           # Dump interval
    CUTOFF_FREQ   = 30          # Cutoff frequency for plotting (THz)
    
    # Input Files (Filename : Number of Atoms)
    FILE_DICT = {
        'vel_right1_heat.dat': 640,
        # 'another_file.dat': 1000
    }
    TYPE_MASS_MAP = {
        1: 24.30500,  # Mg
        2: 15.99900   # O
    }

        
    # FFT and Smoothing Parameters
    FRAMES_READ  = int(NUMBER_STEPS / INTERVAL)
    WINDOW_SIZE  = int(FRAMES_READ / 8)    # FFT Window size (points)
    OVERLAP_DIV  = 64                      # Control overlap (Window / Div = Step size)
    SMOOTH_WIDTH = 1.0                     # Smoothing width factor
    # ================= MAIN LOOP =================
    try:
        for filename, n_atoms in FILE_DICT.items():
            
            # Output Directory Generation
            #timestamp = time.strftime("%Y%m%d_%H%M%S")
            OUTPUT_DIR = f'Result_VDOS_For_{filename}' #f'Result_VDOS_{timestamp}'

            # Calculate sliding shift step
            shift_step = int(WINDOW_SIZE / OVERLAP_DIV)
            if shift_step < 1: shift_step = 1

            # 1. Initialize
            calc = VDOSCalculator(
                n_atoms=n_atoms,
                dt_ps=DT_PS,
                interval=INTERVAL,
                total_frames=FRAMES_READ,
                fft_window_size=WINDOW_SIZE,
                overlap_shift=shift_step,
                smooth_width=SMOOTH_WIDTH,
                output_dir=OUTPUT_DIR
            )
            
            calc.log(f"=== Processing Task: {filename} ===")
            calc.log(f"Window Size: {WINDOW_SIZE}, Total Frames: {FRAMES_READ}, Shift Step: {shift_step}")

            # 2. Read Data
            velocity_data = calc.load_velocity(filename)

            # 3. Load Mass
            mass = calc.load_mass(filename, type_mass_map=TYPE_MASS_MAP)
            
            # 3. Compute VDOS
            freq_axis, vdos_mat = calc.compute_vdos(velocity_data, mass)
            
            # 4. Smooth and Save
            final_data = calc.smooth_and_save(freq_axis, vdos_mat, filename)
            
            # 5. Plot
            calc.plot_result(final_data, CUTOFF_FREQ, filename)
            
            calc.close()

        print(f"\nAll tasks finished. Results saved in: {OUTPUT_DIR}")

    except Exception as e:
        print(f"\nCritical Error: {e}")
        import traceback
        traceback.print_exc()