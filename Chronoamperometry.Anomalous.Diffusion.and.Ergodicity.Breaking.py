"""
Anomalous Diffusion & ML Classification Pipeline
================================================
Author: [Your Name]
Description: 
    Analyzes Chronoamperometry data to extract anomalous diffusion parameters 
    (alpha, lambda) and performs Logistic Regression classification.
    
    Physics Model:
    i(t) ~ t^(-alpha/2)  (Power Law Decay)
    psi(t) ~ t^{-(1+alpha)} (Waiting Time Distribution)

    Machine Learning:
    Classifies 'Positive' vs 'Negative' samples using physics-informed features:
    1. Alpha (Anomalous Exponent)
    2. Total Charge Q (Integral of current)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.integrate import cumulative_trapezoid
from scipy.stats import linregress
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.special import gamma
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from pathlib import Path
from typing import Tuple, List, Optional

# ==========================================
# CONFIGURATION
# ==========================================
# Use relative paths for GitHub portability
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"

# List your files here. The script assumes they are in the /data folder.
FILE_NAMES = [
    "Dados_exp_amb.controlado.xlsx",      # Run 1
    "Dados_exp_amb.controlado (1).xlsx",  # Run 2
    "Dados_exp_amb.controlado (2).xlsx"   # Run 3
]

OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# Analysis Parameters
REGRESSION_WINDOW = (20.0, 120.0) # Seconds (Linear region for power law)
SMOOTHING_WINDOW = 51
POLY_ORDER = 2

# ==========================================
# 1. COMPUTATIONAL PHYSICS FUNCTIONS
# ==========================================

def log_space_smoothing(t: np.ndarray, y: np.ndarray, 
                       window_size: int = SMOOTHING_WINDOW, 
                       poly_order: int = POLY_ORDER) -> Tuple[np.ndarray, np.ndarray]:
    """
    Smoothes data preserving power-law physics by resampling on a logarithmic grid.
    """
    valid = (t > 0) & (np.isfinite(y))
    if np.sum(valid) < 10:
        return t, y 
        
    t_c, y_c = t[valid], y[valid]
    
    try:
        # Log Grid Interpolation
        t_log = np.logspace(np.log10(t_c[0]), np.log10(t_c[-1]), len(t_c))
        f_int = interp1d(t_c, y_c, kind='linear', fill_value="extrapolate")
        y_log_domain = f_int(t_log)
        
        # Safety: Ensure window_size is odd and smaller than data length
        eff_win = min(window_size, len(t_c))
        if eff_win % 2 == 0: eff_win -= 1
        if eff_win < 3: eff_win = 3
        
        y_smooth_log = savgol_filter(y_log_domain, window_length=eff_win, polyorder=poly_order)
        
        # Return to Linear Domain
        f_back = interp1d(t_log, y_smooth_log, kind='linear', fill_value="extrapolate")
        y_smooth_valid = f_back(t_c)
        
        # Reconstruct full array
        y_full = np.copy(y)
        y_full[valid] = y_smooth_valid
        return t, y_full
        
    except Exception as e:
        print(f"   [WARN] Smoothing failed: {e}. Returning raw data.")
        return t, y

def calculate_physics_metrics(t: np.ndarray, I_raw: np.ndarray, I_smooth: np.ndarray):
    """
    Calculates Charge Q(t) and local Cottrell Exponent Alpha(t).
    alpha(t) = - d(ln I) / d(ln t)
    """
    # 1. Integral (Charge) - Use Raw Data
    I_filled = np.nan_to_num(I_raw)
    Q = cumulative_trapezoid(I_filled, t, initial=0)
    
    # 2. Logarithmic Derivative (Alpha) - Use Smoothed Data
    valid_mask = (t > 0) & (I_smooth > 0)
    
    if np.sum(valid_mask) < 10:
        return Q, np.zeros_like(t)

    t_valid = t[valid_mask]
    I_valid = I_smooth[valid_mask]

    # Resample on Log Grid
    num_points = len(t_valid)
    t_log_grid = np.logspace(np.log10(t_valid[0]), np.log10(t_valid[-1]), num_points)
    
    f_int = interp1d(t_valid, I_valid, kind='linear', fill_value="extrapolate")
    I_log_domain = f_int(t_log_grid)

    x_log = np.log(t_log_grid)
    y_log = np.log(I_log_domain)
    dx = x_log[1] - x_log[0]

    # Window safety
    window_size = 71
    if window_size >= num_points: 
        window_size = num_points // 2 * 2 + 1
    if window_size < 3: 
        window_size = 3
    
    # Derivative
    d_log_I = savgol_filter(y_log, window_length=window_size, polyorder=2, deriv=1, delta=dx)
    
    # alpha = - d(ln I) / d(ln t)
    alpha_log_domain = -d_log_I

    # Return to linear grid
    f_back = interp1d(t_log_grid, alpha_log_domain, kind='linear', fill_value="extrapolate")
    alpha_final = np.zeros_like(t)
    alpha_final[valid_mask] = f_back(t_valid)
    
    return Q, alpha_final

def extract_regression_params(time: np.ndarray, current: np.ndarray):
    """
    Performs linear regression on log-log data in the specified time window 
    to extract static alpha and lambda parameters.
    """
    mask_reg = (time >= REGRESSION_WINDOW[0]) & (time <= REGRESSION_WINDOW[1])
    
    if np.sum(mask_reg) < 10:
        return np.nan, np.nan, 0.0

    # Log-Log Transformation
    x = np.log(time[mask_reg])
    y = np.log(current[mask_reg])
    
    # Linear Regression: ln(i) = Intercept + Slope * ln(t)
    slope, intercept, r_val, _, _ = linregress(x, y)
    
    # Physics Extraction
    # Slope = -alpha / 2  =>  alpha = -2 * slope
    phys_alpha = -2 * slope
    
    # Intercept = ln( Lambda / Gamma(1 - alpha/2) )
    term_gamma = gamma(1 - phys_alpha/2)
    phys_lambda = np.exp(intercept) * term_gamma
    
    return phys_alpha, phys_lambda, r_val

# ==========================================
# 2. MAIN PROCESSING LOOP
# ==========================================

def run_analysis():
    print("--- Starting Physics-Informed Analysis ---")
    
    history_neg = []
    history_pos = []
    
    # A. PROCESS INDIVIDUAL RUNS (Physics Metrics)
    # ---------------------------------------------------------
    for i, file_name in enumerate(FILE_NAMES):
        file_path = DATA_DIR / file_name
        if not file_path.exists():
            print(f"   [ERROR] File not found: {file_name}")
            continue
            
        print(f"\nProcessing Run {i+1}: {file_name}")
        df = pd.read_excel(file_path)

        # Extraction
        t = df['Time'].values
        cols_pos = [c for c in df.columns if 'Positive' in c]
        cols_neg = [c for c in df.columns if 'Negative' in c]

        # Ensemble Statistics
        I_pos_raw = df[cols_pos].mean(axis=1).values
        I_neg_raw = df[cols_neg].mean(axis=1).values
        I_pos_std = df[cols_pos].std(axis=1).values
        I_neg_std = df[cols_neg].std(axis=1).values

        # Smoothing
        _, I_pos_smooth = log_space_smoothing(t, I_pos_raw)
        _, I_neg_smooth = log_space_smoothing(t, I_neg_raw)

        # Physics Metrics (Alpha(t) and Q(t))
        Q_pos, alpha_pos = calculate_physics_metrics(t, I_pos_raw, I_pos_smooth)
        Q_neg, alpha_neg = calculate_physics_metrics(t, I_neg_raw, I_neg_smooth)

        # Regression Extraction (Static Alpha)
        p_alpha_neg, p_lambda_neg, r_neg = extract_regression_params(t, I_neg_raw)
        p_alpha_pos, p_lambda_pos, r_pos = extract_regression_params(t, I_pos_raw)

        print(f"   [REGRESSION] {REGRESSION_WINDOW[0]}-{REGRESSION_WINDOW[1]}s")
        print(f"     NEG: α={p_alpha_neg:.3f}, λ={p_lambda_neg:.3e} (R²={r_neg**2:.3f})")
        print(f"     POS: α={p_alpha_pos:.3f}, λ={p_lambda_pos:.3e} (R²={r_pos**2:.3f})")

        # Store for Global Plotting
        if not np.isnan(p_alpha_neg):
            history_neg.append((i+1, t, p_lambda_neg, p_alpha_neg))
        if not np.isnan(p_alpha_pos):
            history_pos.append((i+1, t, p_lambda_pos, p_alpha_pos))

        # --- PLOT 1: RUN DASHBOARD ---
        plt.style.use('default')
        fig, axs = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
        
        # Panel 1: Local Derivative (Alpha)
        ax = axs[0]
        ax.set_xlim(t[0], t[-1])
        ax.plot(t, alpha_neg, 'b-', lw=1.5, label='Negative')
        ax.plot(t, alpha_pos, 'r-', lw=1.5, label='Positive')
        ax.axhline(0.5, color='k', ls='--', label='Cottrell (0.5)')
        ax.set_ylim(0, 1) 
        ax.set_title(r'Local Derivative ($-d\ln I / d\ln t$)')
        ax.set_ylabel(r'$\alpha(t)$')
        ax.set_xlabel('Time [s]')
        ax.legend(fontsize='small')
        ax.grid(True, alpha=0.3)

        # Panel 2: Log-Log Current + STD
        ax = axs[1]
        t_start = t[t > 0][0] if t[0] <= 0 else t[0]
        ax.set_xlim(t_start, t[-1])
        
        # Shading helpers
        pos_lower = np.maximum(I_pos_raw - I_pos_std, I_pos_raw * 0.01)
        neg_lower = np.maximum(I_neg_raw - I_neg_std, I_neg_raw * 0.01)

        ax.loglog(t, I_neg_raw, 'b-', lw=1.5, label='Neg Ensemble')
        ax.loglog(t, I_pos_raw, 'r-', lw=1.5, label='Pos Ensemble')
        ax.fill_between(t, neg_lower, I_neg_raw + I_neg_std, color='blue', alpha=0.2)
        ax.fill_between(t, pos_lower, I_pos_raw + I_pos_std, color='red', alpha=0.2)
        ax.set_ylabel('Current Log(i) [µA]')
        ax.set_xlabel('Time Log(t) [s]')
        ax.set_title('Electrical Current')
        ax.legend(loc='upper right', fontsize='small')
        ax.grid(True, which="both", alpha=0.3)

        # Panel 3: Charge Q
        ax = axs[2]
        ax.plot(t, Q_neg, 'b-', lw=1.5, label=f'Q_neg: {Q_neg[-1]:.0f} µC')
        ax.plot(t, Q_pos, 'r-', lw=1.5, label=f'Q_pos: {Q_pos[-1]:.0f} µC')
        ax.set_xlim(left=0, right=t[-1])
        ax.set_ylabel('Charge Q [µC]')
        ax.set_xlabel('Time [s]')
        ax.set_title(r'Total Transferred Charge ($\int i dt$)')
        ax.legend(loc='lower right', fontsize='small')
        ax.grid(True, alpha=0.3)
        
        save_path = OUTPUT_DIR / f'Analysis_Run_{i+1}.png'
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"   Saved Dashboard: {save_path.name}")

    # B. WAITING TIME DISTRIBUTION (Aging Plot)
    # ---------------------------------------------------------
    print("\nGenerating Waiting Time Distribution Plot...")
    plt.figure(figsize=(10, 8))
    
    blues = cm.Blues(np.linspace(0.4, 0.8, len(history_neg)))
    reds = cm.Reds(np.linspace(0.5, 0.9, len(history_pos)))

    # Fickian Reference (Slope -0.5 corresponds to Alpha 1.0 in WTD?)
    # Note: If i ~ t^-0.5, then psi ~ t^-1.5
    if history_neg:
        t_ref = history_neg[0][1]
        t_theory = t_ref[t_ref > 0]
        psi_fickian = t_theory ** -1.5 
        psi_fickian /= psi_fickian[0] * 10 # Offset for visibility
        plt.loglog(t_theory, psi_fickian, 'k--', lw=2, alpha=0.6, label='Reference Fickian')

    for idx, (run_num, t_run, l_val, p_alpha) in enumerate(history_neg):
        t_safe = t_run[t_run > 0]
        # psi ~ t ^ -(1 + alpha) ?? 
        # Actually standard fractional diffusion is psi ~ t^-(1+alpha)
        # But if i ~ t^-mu, then psi ~ t^-(1+mu) usually.
        # User code used -(1+alpha), preserving that logic.
        psi = t_safe ** -(1 + p_alpha) 
        psi /= psi[0]
        plt.loglog(t_safe, psi, color=blues[idx], lw=2, 
                   label=f'Neg Run {run_num} ($\\alpha={p_alpha:.2f}$)')

    for idx, (run_num, t_run, l_val, p_alpha) in enumerate(history_pos):
        t_safe = t_run[t_run > 0]
        psi = t_safe ** -(1 + p_alpha)
        psi /= psi[0]
        plt.loglog(t_safe, psi, color=reds[idx], lw=2.5, 
                   label=f'Pos Run {run_num} ($\\alpha={p_alpha:.2f}$)')

    plt.xlabel('Time (Log) [s]', fontsize=14)
    plt.ylabel(r'Waiting Time PDF $\psi(t) \sim t^{-(1+\alpha)}$', fontsize=14)
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'Aging_Waiting_Times.png', dpi=300)
    plt.close()

    # C. LOGISTIC REGRESSION (Machine Learning)
    # ---------------------------------------------------------
    print("\nRunning Machine Learning Classification...")
    fig_lr, axs_lr = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    for i, file_name in enumerate(FILE_NAMES):
        file_path = DATA_DIR / file_name
        if not file_path.exists(): continue

        df = pd.read_excel(file_path)
        t = df['Time'].values
        
        # Feature Extraction Helper
        def get_features(current_col):
            # Feature 1: Total Charge Q
            q_total = cumulative_trapezoid(np.nan_to_num(current_col), t, initial=0)[-1]
            
            # Feature 2: Alpha (Slope of log-log)
            mask_reg = (t >= REGRESSION_WINDOW[0]) & (t <= REGRESSION_WINDOW[1])
            i_reg = current_col[mask_reg]
            t_reg = t[mask_reg]
            valid_log = i_reg > 0
            
            if np.sum(valid_log) > 5:
                slope, _, _, _, _ = linregress(np.log(t_reg[valid_log]), np.log(i_reg[valid_log]))
                alpha = -2 * slope
            else:
                alpha = np.nan
            return alpha, q_total

        features, labels = [], []
        
        # Extract Negative Samples (Label 0)
        for col in [c for c in df.columns if 'Negative' in c]:
            a, q = get_features(df[col].values)
            if not np.isnan(a): 
                features.append([a, q])
                labels.append(0)
                
        # Extract Positive Samples (Label 1)
        for col in [c for c in df.columns if 'Positive' in c]:
            a, q = get_features(df[col].values)
            if not np.isnan(a): 
                features.append([a, q])
                labels.append(1)

        X, y = np.array(features), np.array(labels)
        
        # Standardization
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Train Model
        clf = LogisticRegression(random_state=42).fit(X_scaled, y)
        
        # Extract Decision Boundary
        w = clf.coef_[0]
        b = clf.intercept_[0]
        mean, std = scaler.mean_, scaler.scale_
        
        # Physical Threshold: Q = slope * alpha + intercept
        slope_phys = -(w[0] / w[1]) * (std[1] / std[0])
        intercept_phys = mean[1] - (slope_phys * mean[0]) - (b * std[1] / w[1])

        print(f"   [Run {i+1}] Accuracy: {clf.score(X_scaled, y)*100:.1f}%")
        print(f"     > Threshold Eq: Q = ({slope_phys:.2f} * α) + ({intercept_phys:.2f})")
        
        # Visualization
        ax = axs_lr[i]
        
        # Decision Surface
        x_min, x_max = X[:, 0].min()-0.1, X[:, 0].max()+0.1
        y_min, y_max = X[:, 1].min()-500, X[:, 1].max()+500
        xx, yy = np.meshgrid(np.linspace(x_min, x_max, 200), np.linspace(y_min, y_max, 200))
        Z = clf.predict(scaler.transform(np.c_[xx.ravel(), yy.ravel()])).reshape(xx.shape)
        
        ax.contourf(xx, yy, Z, alpha=0.3, cmap=mcolors.ListedColormap(['#e0e0ff', '#ffe0e0']))
        ax.contour(xx, yy, Z, levels=[0.5], colors='k', linestyles='--')
        
        # Scatter Points
        ax.scatter(X[y==0,0], X[y==0,1], c='blue', edgecolors='k', label='Negative', s=60)
        ax.scatter(X[y==1,0], X[y==1,1], c='red', edgecolors='k', label='Positive', s=60)
        
        ax.set_title(f'Run {i+1} Decision Boundary')
        ax.set_xlabel(r'Anomalous Exponent ($\alpha$)')
        ax.grid(True, alpha=0.2)
        if i == 0: 
            ax.set_ylabel(r'Total Charge $Q$ [$\mu$C]')
        ax.legend(loc='lower left', fontsize='small')

    plt.savefig(OUTPUT_DIR / 'ML_Classification_Boundary.png', dpi=300)
    print("\nAnalysis Complete. Results saved to /output.")

if __name__ == "__main__":
    run_analysis()
