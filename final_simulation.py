import numpy as np
import matplotlib.pyplot as plt

# More aggressive aging parameters to reach realistic lifespans
alpha = 0.012  # cellular damage rate  
beta_0 = 0.008  # initial repair rate (lower than damage)
gamma = 0.025   # rapid repair decline
D = 0.012      # diffusion coefficient
theta_c = 0.35  # cells below this are dysfunctional
f_crit = 0.28   # death when this fraction dysfunctional

# Grid
nx = 300
x = np.linspace(0, 1, nx)
dx = x[1] - x[0]

# Time 
t_max = 120
dt = 0.02
nt = int(t_max/dt)

# Initial distribution - young healthy cells
x0 = 0.88
sigma0 = 0.06
rho = np.exp(-(x - x0)**2 / (2*sigma0**2))
rho = rho / (np.sum(rho) * dx)

# Storage
mean_health = []
variance = []
frac_below_threshold = []
time_points = []

# Time evolution
death_time = None
for n in range(nt):
    t = n * dt
    time_points.append(t)
    
    # Calculate statistics
    H = np.sum(x * rho * dx)
    mean_health.append(H)
    variance.append(np.sum(((x - H)**2) * rho * dx))
    frac_below = np.sum(rho[x < theta_c] * dx)
    frac_below_threshold.append(frac_below)
    
    # Check death
    if frac_below >= f_crit and death_time is None:
        death_time = t
        print(f"Death threshold crossed at {t:.1f} years")
        print(f"Fraction below threshold: {frac_below:.3f}")
        print(f"Mean health at death: {H:.3f}")
    
    # Update distribution
    repair = beta_0 * np.exp(-gamma * t)
    v = -alpha * (1 - x) + repair * x  # drift
    
    # Upwind scheme for stability
    rho_new = rho.copy()
    for i in range(1, nx-1):
        if v[i] > 0:
            flux_grad = (v[i]*rho[i] - v[i-1]*rho[i-1]) / dx
        else:
            flux_grad = (v[i+1]*rho[i+1] - v[i]*rho[i]) / dx
        
        diffusion = D * (rho[i+1] - 2*rho[i] + rho[i-1]) / dx**2
        rho_new[i] = rho[i] + dt * (-flux_grad + diffusion)
    
    # Boundary conditions
    rho_new[0] = 0  # absorbing at x=0
    rho_new[-1] = rho_new[-2]  # reflecting at x=1
    
    rho = np.maximum(rho_new, 0)
    rho = rho / (np.sum(rho) * dx)  # renormalize

if death_time is None:
    death_time = t_max
    print(f"Death not reached within {t_max} years")

# Create Figure 1
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: Distribution snapshots
snapshot_times = [0, int(0.25*nt), int(0.5*nt), int(0.75*nt)]
snapshot_labels = ['Birth (0y)', 'Young adult (30y)', 'Middle age (60y)', 'Old age (90y)']
colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']

# Re-run to get snapshots
rho_snapshots = []
x0 = 0.88
sigma0 = 0.06
rho_temp = np.exp(-(x - x0)**2 / (2*sigma0**2))
rho_temp = rho_temp / (np.sum(rho_temp) * dx)
rho_snapshots.append(rho_temp.copy())

snap_idx = 1
for n in range(nt):
    if n == snapshot_times[snap_idx]:
        rho_snapshots.append(rho_temp.copy())
        snap_idx += 1
        if snap_idx >= len(snapshot_times):
            break
    
    t = n * dt
    repair = beta_0 * np.exp(-gamma * t)
    v = -alpha * (1 - x) + repair * x
    
    rho_new = rho_temp.copy()
    for i in range(1, nx-1):
        if v[i] > 0:
            flux_grad = (v[i]*rho_temp[i] - v[i-1]*rho_temp[i-1]) / dx
        else:
            flux_grad = (v[i+1]*rho_temp[i+1] - v[i]*rho_temp[i]) / dx
        diffusion = D * (rho_temp[i+1] - 2*rho_temp[i] + rho_temp[i-1]) / dx**2
        rho_new[i] = rho_temp[i] + dt * (-flux_grad + diffusion)
    
    rho_new[0] = 0
    rho_new[-1] = rho_new[-2]
    rho_temp = np.maximum(rho_new, 0)
    rho_temp = rho_temp / (np.sum(rho_temp) * dx)

for i, (snap, label, color) in enumerate(zip(rho_snapshots, snapshot_labels, colors)):
    axes[0].plot(x, snap, label=label, color=color, linewidth=2.5, alpha=0.8)

axes[0].axvline(theta_c, color='black', linestyle='--', linewidth=2, label=f'Dysfunction threshold θc={theta_c}')
axes[0].set_xlabel('Cellular Health State (x)', fontsize=13, fontweight='bold')
axes[0].set_ylabel('Probability Density ρ(x,t)', fontsize=13, fontweight='bold')
axes[0].set_title('(A) Cellular Health Distribution Evolution', fontsize=14, fontweight='bold')
axes[0].legend(loc='upper left', fontsize=10)
axes[0].grid(alpha=0.3, linestyle=':')
axes[0].set_xlim([0, 1])

# Panel B: Mean health
axes[1].plot(time_points, mean_health, color='#2E86AB', linewidth=2.5)
axes[1].axhline(theta_c, color='red', linestyle='--', linewidth=2, label=f'Dysfunction threshold')
if death_time < t_max:
    axes[1].axvline(death_time, color='black', linestyle=':', linewidth=2, label=f'Death at {death_time:.0f}y')
axes[1].set_xlabel('Time (years)', fontsize=13, fontweight='bold')
axes[1].set_ylabel('Mean Cellular Health H(t)', fontsize=13, fontweight='bold')
axes[1].set_title('(B) Mean Health Decline with Age', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=10)
axes[1].grid(alpha=0.3, linestyle=':')

# Panel C: Fraction below threshold
axes[2].plot(time_points, frac_below_threshold, color='#C73E1D', linewidth=2.5)
axes[2].axhline(f_crit, color='black', linestyle='--', linewidth=2, label=f'Death threshold fcrit={f_crit}')
if death_time < t_max:
    axes[2].axvline(death_time, color='purple', linestyle=':', linewidth=2, label=f'Death at {death_time:.0f}y')
axes[2].set_xlabel('Time (years)', fontsize=13, fontweight='bold')
axes[2].set_ylabel('Fraction of Dysfunctional Cells', fontsize=13, fontweight='bold')
axes[2].set_title('(C) Approach to Mortality Threshold', fontsize=14, fontweight='bold')
axes[2].legend(fontsize=10)
axes[2].grid(alpha=0.3, linestyle=':')

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/figure1_cellular_health_simulation.png', dpi=300, bbox_inches='tight')

print(f"\n{'='*70}")
print(f"SIMULATION RESULTS")
print(f"{'='*70}")
print(f"Predicted lifespan: {death_time:.1f} years")
print(f"Initial mean health H(0): {mean_health[0]:.3f}")
print(f"Final mean health H(T): {mean_health[-1]:.3f}")
print(f"Initial variance σ²(0): {variance[0]:.5f}")
print(f"Final variance σ²(T): {variance[-1]:.5f}")
print(f"Variance increase: {(variance[-1]/variance[0] - 1)*100:.1f}%")
print(f"{'='*70}")
print(f"\nFigure saved successfully!")
