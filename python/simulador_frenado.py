"""
Simulador del sistema de frenado - Chevrolet Sail
Proyecto ABPro: Optimización cinemática y dinámica del sistema de frenado
mediante modelado matemático y cálculo diferencial.

Uso:  python simulador_frenado.py [carpeta_salida]
Calcula los resultados del modelo y guarda los gráficos en PNG.
Los parámetros marcados como SUPUESTO deben reemplazarse por mediciones reales.
"""
import sys, os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "graficos"
os.makedirs(OUT, exist_ok=True)

# ---------------- Parámetros ----------------
g = 9.81
m = 1103.0            # kg, masa del vehículo (dato del proyecto)
a_f = 6.0             # m/s^2, deceleración de frenado de referencia (dato del proyecto)
L = 2.50              # m, distancia entre ejes (ficha Sail 1.5 LT)
h = 0.55              # m, altura del centro de gravedad (SUPUESTO)
w_del = 0.60          # fracción estática de peso en eje delantero (SUPUESTO, tracción delantera)
rho, Cd, A = 1.2, 0.33, 2.1   # aire, coef. arrastre, área frontal (SUPUESTOS típicos)
beta = 0.5 * rho * Cd * A     # kg/m
alpha = m * a_f               # N, fuerza de frenado + rodadura
t_reac = 1.0                  # s, tiempo de reacción (CONASET)
m_disco = 5.0                 # kg por disco delantero (SUPUESTO)
c_fe = 460.0                  # J/(kg K), fundición gris
gamma = 0.90                  # fracción del calor que va al disco (SUPUESTO)
P_max = 80.0                  # bar, presión máxima de mordida (SUPUESTO)
tau_goma, tau_acero = 0.12, 0.08   # s, constantes de tiempo de la línea (SUPUESTO)

kmh = lambda v: v / 3.6

# ---------------- Paleta ----------------
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 10,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.8, "lines.linewidth": 2,
    "legend.frameon": False, "figure.dpi": 200,
})

def guardar(fig, nombre):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, nombre), dpi=200)
    plt.close(fig)

R = {}

# ---------------- 1. Aceleración 0-100 (cinemática del proyecto) ----------------
a_acc = kmh(100) / 11.3
R["a_acc"] = a_acc
R["F_acc"] = m * a_acc

# ---------------- 2. Frenado: modelo básico y modelo con arrastre ----------------
def d_basico(v0, a=a_f):
    return v0**2 / (2 * a)

def d_arrastre(v0):
    return m / (2 * beta) * math.log(1 + beta * v0**2 / alpha)

def t_arrastre(v0):
    return m / math.sqrt(alpha * beta) * math.atan(v0 * math.sqrt(beta / alpha))

tabla = []
for V in (40, 60, 80, 100):
    v0 = kmh(V)
    tabla.append(dict(V=V, v0=v0, t=v0 / a_f, d=d_basico(v0), d_reac=v0 * t_reac,
                      d_tot=d_basico(v0) + v0 * t_reac, d_arr=d_arrastre(v0),
                      t_arr=t_arrastre(v0), Ec=0.5 * m * v0**2,
                      sens=v0 / a_f / 3.6))
R["tabla"] = tabla

# Curvas v(t), x(t) desde 100 km/h
v0 = kmh(100)
ts = np.linspace(0, v0 / a_f, 300)
v_b = v0 - a_f * ts
x_b = v0 * ts - 0.5 * a_f * ts**2
k = math.sqrt(alpha * beta) / m
th0 = math.atan(v0 * math.sqrt(beta / alpha))
ta = np.linspace(0, t_arrastre(v0), 300)
v_a = math.sqrt(alpha / beta) * np.tan(th0 - k * ta)
x_a = m / beta * (np.log(np.cos(th0 - k * ta)) - np.log(math.cos(th0)))

fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
ax[0].plot(ts, v_b * 3.6, color=C1, label="m·dv/dt = −F (constante)")
ax[0].plot(ta, v_a * 3.6, color=C2, label="m·dv/dt = −(α + βv²)")
ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("v(t) (km/h)")
ax[0].legend(fontsize=8, loc="upper right")
ax[1].plot(ts, x_b, color=C1); ax[1].plot(ta, x_a, color=C2)
ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("x(t) (m)")
ax[1].annotate(f"{x_b[-1]:.1f} m", (ts[-1], x_b[-1]), xytext=(-44, 6), textcoords="offset points", color=INK2, fontsize=8)
ax[1].annotate(f"{x_a[-1]:.1f} m", (ta[-1], x_a[-1]), xytext=(-48, -16), textcoords="offset points", color=INK2, fontsize=8)
guardar(fig, "g1_frenado_v_x.png")

# Distancia de detención vs velocidad
Vs = np.linspace(0, 120, 200); vv = kmh(Vs)
fig, ax = plt.subplots(figsize=(6.4, 3.0))
ax.plot(Vs, vv * t_reac + vv**2 / (2 * a_f), color=C1, label="Detención total (reacción 1 s + frenado)")
ax.plot(Vs, vv**2 / (2 * a_f), color=C2, label="Solo frenado, d = v₀²/(2a)")
for r in tabla:
    ax.plot(r["V"], r["d_tot"], "o", color=C1, ms=5, mec="white", mew=1.2)
    ax.annotate(f"{r['d_tot']:.0f} m", (r["V"], r["d_tot"]), xytext=(-10, 7), textcoords="offset points", fontsize=8, color=INK2)
ax.set_xlabel("Velocidad inicial (km/h)"); ax.set_ylabel("Distancia (m)")
ax.legend(fontsize=8, loc="upper left")
guardar(fig, "g2_distancia_velocidad.png")

# ---------------- 3. Influencia de la masa ----------------
masas = [("Vehículo solo", 1103), ("Con conductor", 1103 + 75), ("Carga máxima", 1103 + 357)]
R["masas"] = []
for nombre, mm in masas:
    dv = kmh(100)
    R["masas"].append(dict(n=nombre, m=mm, a_neum=0.8 * g, d_neum=dv**2 / (2 * 0.8 * g),
                           a_sist=alpha / mm, d_sist=mm * dv**2 / (2 * alpha),
                           Ec=0.5 * mm * dv**2))

# ---------------- 4. Transferencia de carga ----------------
W = m * g
dW = m * a_f * h / L
R["W"], R["dW"] = W, dW
R["N_del_est"], R["N_tra_est"] = w_del * W, (1 - w_del) * W
R["N_del_din"], R["N_tra_din"] = w_del * W + dW, (1 - w_del) * W - dW
phi = R["N_del_din"] / W
R["phi"] = phi

# ---------------- 5. Presión hidráulica ----------------
tp = np.linspace(0, 0.6, 300)
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
for tau, col, lab in ((tau_goma, C2, "Manguera de goma (τ = 0,12 s)"), (tau_acero, C1, "Línea acerada (τ = 0,08 s)")):
    ax[0].plot(tp, P_max * (1 - np.exp(-tp / tau)), color=col, label=lab)
    ax[1].plot(tp, P_max / tau * np.exp(-tp / tau), color=col)
ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("P(t) (bar)"); ax[0].legend(fontsize=8, loc="lower right")
ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("dP/dt (bar/s)")
guardar(fig, "g3_presion.png")
R["dPdt0"] = {"goma": P_max / tau_goma, "acero": P_max / tau_acero}
R["t90"] = {"goma": tau_goma * math.log(10), "acero": tau_acero * math.log(10)}
R["delta_d_linea"] = {V: kmh(V) * (tau_goma - tau_acero) for V in (40, 60, 80, 100)}

# ---------------- 6. Temperatura del disco ----------------
v0 = kmh(100)
Ec = 0.5 * m * v0**2
q_disco = gamma * phi * Ec / 2
dT = q_disco / (m_disco * c_fe)
dTdt0 = gamma * phi * m * a_f * v0 / (2 * m_disco * c_fe)
R.update(Ec100=Ec, q_disco=q_disco, dT=dT, dTdt0=dTdt0)

tt = np.linspace(0, v0 / a_f, 200)
kT = gamma * phi * m * a_f / (2 * m_disco * c_fe)
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
ax[0].plot(tt, kT * (v0 - a_f * tt), color=C2)
ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("dT/dt (°C/s)")
ax[1].plot(tt, 20 + kT * (v0 * tt - 0.5 * a_f * tt**2), color=C2)
ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("T(t) del disco (°C)")
guardar(fig, "g4_temperatura_frenada.png")

# Frenadas repetidas (enfriamiento de Newton entre frenadas)
T_amb, n_fren, t_ciclo = 20.0, 10, 60.0
casos = [("Disco macizo (τ = 250 s)", 250, 5.0, C2),
         ("Disco ventilado (τ = 150 s)", 150, 5.0, C1),
         ("Ventilado ranurado, +20 % masa (τ = 130 s)", 130, 6.0, C3)]
fig, ax = plt.subplots(figsize=(6.4, 3.0))
R["repetidas"] = []
for nombre, tau_c, md, col in casos:
    dT_i = q_disco / (md * c_fe)
    T, tl, Tl, picos = T_amb, [], [], []
    t_fr = v0 / a_f
    for i in range(n_fren):
        t0 = i * t_ciclo
        for s in np.linspace(0, t_fr, 20):
            tl.append(t0 + s); Tl.append(T + dT_i * (v0 * s - 0.5 * a_f * s**2) / (v0 * t_fr / 2))
        T = T + dT_i; picos.append(T)
        for s in np.linspace(0, t_ciclo - t_fr, 40)[1:]:
            tl.append(t0 + t_fr + s); Tl.append(T_amb + (T - T_amb) * math.exp(-s / tau_c))
        T = T_amb + (T - T_amb) * math.exp(-(t_ciclo - t_fr) / tau_c)
    ax.plot(tl, Tl, color=col, label=nombre, linewidth=1.6)
    R["repetidas"].append(dict(n=nombre, dT=dT_i, pico=max(picos)))
ax.axhline(350, color=INK2, lw=1, ls="--")
ax.text(600, 358, "Zona de riesgo de fading (≈ 350 °C)", fontsize=8, color=INK2, ha="right")
ax.set_xlabel("t (s)"); ax.set_ylabel("Temperatura del disco (°C)")
ax.legend(fontsize=8, loc="upper left")
ax.set_ylim(0, 540)
guardar(fig, "g5_frenadas_repetidas.png")

if __name__ == "__main__":
    import json
    print(json.dumps(R, indent=1, ensure_ascii=False, default=float))
