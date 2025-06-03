import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# ---------------------------
# Parametry i dane
# ---------------------------
a = np.array([-1.4138, 0.6065])
b = np.array([0.1044, 0.0883])

kk = 250
N = 30
D = 30
Nu = 15
lambda_ = 4
dv = 0.1

# ---------------------------
# Charakterystyka statyczna
# ---------------------------
def static_nonlinear(v):
    return 0.3163 * v / np.sqrt(0.1 + 0.9 * v**2)

# ---------------------------
# FIS (Sugeno, ręczna implementacja)
# ---------------------------
def evalfis(v):
    w1 = 1 / (1 + np.exp(5 * v))
    w2 = 1 / (1 + np.exp(-5 * v))
    y1 = -0.3289
    y2 = 0.3289
    return (w1 * y1 + w2 * y2) / (w1 + w2)

# ---------------------------
# Odpowiedź skokowa
# ---------------------------
s = np.zeros(D)
s[0] = b[0]
for k in range(1, D):
    if k == 1:
        s[k] = -a[0] * s[k - 1] + b[0] + b[1]
    else:
        s[k] = -a[0] * s[k - 1] - a[1] * s[k - 2] + b[0] + b[1]

# Macierz M i Mp
M = np.zeros((N, Nu))
Mp = np.zeros((N, D - 1))
for i in range(N):
    for j in range(Nu):
        if i >= j:
            M[i, j] = s[i - j]

for i in range(N):
    for j in range(D - 1):
        Mp[i, j] = s[i + j + 1] - s[j] if i + j + 1 < D else s[-1] - s[j]

K = np.linalg.inv(M.T @ M + lambda_ * np.eye(Nu)) @ M.T
ke = np.sum(K[0, :])
ku = K[0, :] @ Mp

# ---------------------------
# Funkcje pomocnicze
# ---------------------------
def fuzzy_free_response(v, u, y, k):
    v0 = np.zeros(N)
    y0 = np.zeros(N)
    u0 = u[k - 1]
    for i in range(N):
        if i == 0:
            v0[i] = -a[0] * v[k] - a[1] * v[k - 1] + b[0] * u0 + b[1] * u0
        elif i == 1:
            v0[i] = -a[0] * v0[i - 1] - a[1] * v[k] + b[0] * u0 + b[1] * u0
        else:
            v0[i] = -a[0] * v0[i - 1] - a[1] * v0[i - 2] + b[0] * u0 + b[1] * u0
        y0[i] = evalfis(v0[i])
    return y0

def nonlinear_prediction(v, u, du, k):
    y_pred = np.zeros(N)
    u_pred = np.zeros(N)
    v_pred = np.zeros(N)
    for i in range(N):
        du_i = du[i] if i < len(du) else du[-1]
        u_pred[i] = u[k - 1] + np.sum(du[:i + 1])

        if i == 0:
            v_pred[i] = -a[0]*v[k] - a[1]*v[k - 1] + b[0]*u_pred[i] + b[1]*u[k - 1]
        elif i == 1:
            v_pred[i] = -a[0]*v_pred[i-1] - a[1]*v[k] + b[0]*u_pred[i - 1] + b[1]*u_pred[i - 1]
        else:
            v_pred[i] = -a[0]*v_pred[i - 1] - a[1]*v_pred[i - 2] + b[0]*u_pred[i - 1] + b[1]*u_pred[i - 2]

        y_pred[i] = static_nonlinear(v_pred[i])
    return y_pred

# ---------------------------
# Sterowanie LMPC, FMPC, NMPC
# ---------------------------
y_zad = np.ones(kk) * 0.3
y_zad[kk // 2:] = 0

v = np.zeros((kk, 3))
u = np.zeros((kk, 3))
y = np.zeros((kk, 3))
du_up = np.zeros((D - 1, 3))

for k in range(2, kk):
    # --- LMPC ---
    v[k, 0] = -a[0] * v[k - 1, 0] - a[1] * v[k - 2, 0] + b[0] * u[k - 1, 0] + b[1] * u[k - 2, 0]
    y[k, 0] = static_nonlinear(v[k, 0])
    e = y_zad[k] - y[k, 0]
    du_k = ke * e - ku @ du_up[:, 0]
    u[k, 0] = u[k - 1, 0] + du_k
    du_up[:, 0] = np.insert(du_k, 0, 0)[:-1]

    # --- FMPC ---
    v[k, 1] = -a[0] * v[k - 1, 1] - a[1] * v[k - 2, 1] + b[0] * u[k - 1, 1] + b[1] * u[k - 2, 1]
    y[k, 1] = static_nonlinear(v[k, 1])
    dydv = (evalfis(v[k, 1]) - evalfis(v[k, 1] - dv)) / dv
    M_new = dydv * M
    K_new = np.linalg.inv(M_new.T @ M_new + lambda_ * np.eye(Nu)) @ M_new.T
    y0 = fuzzy_free_response(v[:, 1], u[:, 1], y[:, 1], k)
    du_k = K_new[0, :] @ (y_zad[k] - y0)
    u[k, 1] = u[k - 1, 1] + du_k
    du_up[:, 1] = np.insert(du_k, 0, 0)[:-1]

    # --- NMPC ---
    v[k, 2] = -a[0] * v[k - 1, 2] - a[1] * v[k - 2, 2] + b[0] * u[k - 1, 2] + b[1] * u[k - 2, 2]
    y[k, 2] = static_nonlinear(v[k, 2])

    def cost_func(du):
        y_pred = nonlinear_prediction(v[:, 2], u[:, 2], du, k)
        return np.sum((y_zad[k] - y_pred)**2) + lambda_ * np.sum(du**2)

    bounds = [(-0.1, 0.1)] * Nu
    res = minimize(cost_func, np.zeros(Nu), bounds=bounds, method='SLSQP', options={'disp': False})
    du_k = res.x[0]
    u[k, 2] = u[k - 1, 2] + du_k
    du_up[:, 2] = np.insert(du_k, 0, 0)[:-1]

# ---------------------------
# Wizualizacja
# ---------------------------
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(y_zad, 'k', label='y_zad')
plt.plot(y[:, 0], 'r', label='LMPC')
plt.plot(y[:, 1], 'b', label='FMPC')
plt.plot(y[:, 2], 'g', label='NMPC')
plt.xlabel('k')
plt.ylabel('y(k)')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.step(np.arange(kk), u[:, 0], 'r', label='LMPC')
plt.step(np.arange(kk), u[:, 1], 'b', label='FMPC')
plt.step(np.arange(kk), u[:, 2], 'g', label='NMPC')
plt.xlabel('k')
plt.ylabel('u(k)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()