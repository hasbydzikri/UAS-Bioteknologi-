import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parameter
V1max = 5.0
Km1   = 2.0
Ki    = 3.0
X     = 10.0   # konstan
k2    = 1.0
k3    = 0.8
k4    = 0.3

def system(t, y):
    A, B, P = y

    # v1 dengan inhibisi non-kompetitif oleh P
    v1 = (V1max * X) / ((Km1 + X) * (1 + P / Ki))
    v2 = k2 * A
    v3 = k3 * B
    v4 = k4 * A

    dA_dt = v1 - v2 - v4
    dB_dt = v2 - v3
    dP_dt = v3

    return [dA_dt, dB_dt, dP_dt]

# Kondisi awal dan rentang waktu
y0 = [0.0, 0.0, 0.0]
t_span = (0, 50)
t_eval = np.linspace(0, 50, 1000)

# Integrasi ODE
sol = solve_ivp(system, t_span, y0, t_eval=t_eval, method='RK45', dense_output=True)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(sol.t, sol.y[0], label='[A] — Intermediet A', color='steelblue', linewidth=2)
ax.plot(sol.t, sol.y[1], label='[B] — Intermediet B', color='darkorange', linewidth=2)
ax.plot(sol.t, sol.y[2], label='[P] — Produk P', color='forestgreen', linewidth=2)

ax.set_xlabel('Waktu (satuan waktu)', fontsize=13)
ax.set_ylabel('Konsentrasi (satuan konsentrasi)', fontsize=13)
ax.set_title('Simulasi Dinamika Sistem Metabolik dengan\nInhibisi Feedback Non-Kompetitif P pada v1', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('simulasi_metabolik.png', dpi=150)
plt.show()

# Print nilai steady-state
print(f"Nilai pada t=50:")
print(f"  [A] = {sol.y[0,-1]:.4f}")
print(f"  [B] = {sol.y[1,-1]:.4f}")
print(f"  [P] = {sol.y[2,-1]:.4f}")
