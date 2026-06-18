"""
Precessão de Thomas — Relatividade Especial e Acoplamento Spin-Órbita
======================================================================
Módulo 1 — Geometria: precessão ao longo de órbita circular
Módulo 2 — Frequências: Ω_Thomas, Ω_Larmor e Ω_total vs velocidade
Módulo 3 — Spin ao longo de 3 órbitas (numérico vs analítico)
Módulo 4 — Splitting spin-órbita com e sem o fator ½ de Thomas
Módulo 5 — Animação: spin + órbita em tempo real

Uso:
    python dirac_thomas.py
    python dirac_thomas.py --v 0.8
    python dirac_thomas.py --module 5 --save
    python dirac_thomas.py --help
"""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import argparse
import platform
import warnings
warnings.filterwarnings("ignore")

if platform.system() == "Darwin":
    matplotlib.use("MacOSX")

ALPHA    = 1.0 / 137.035999084
HART_EV  = 27.211396


# ============================================================================
# Física central
# ============================================================================

def gamma(v):
    v = np.clip(np.asarray(v, float), 0, 1 - 1e-10)
    return 1.0 / np.sqrt(1.0 - v**2)


def thomas_rate(v, omega_orb=1.0):
    """Ω_T = −(γ−1)·ω_orb  (precessão oposta ao orbital)."""
    return -(gamma(v) - 1.0) * omega_orb


def larmor_rate(v, omega_orb=1.0, g=2.0):
    """Ω_L = (g/2)·ω_orb/γ  (campo B efetivo no ref. do elétron)."""
    return (g / 2.0) * omega_orb / gamma(v)


def total_rate(v, omega_orb=1.0, g=2.0):
    OL = larmor_rate(v, omega_orb, g)
    OT = thomas_rate(v, omega_orb)
    return OL, OT, OL + OT


def energy_spinorbit(n, l, j, Z=1, include_thomas=True):
    """
    Energia de acoplamento spin-órbita (perturbativo, u.a.).
    H_SO = (α²/2) · (Z/r³) · L·S
    ⟨1/r³⟩ = Z³ / (n³ l(l+½)(l+1))
    Fator ½ de Thomas: sem Thomas → ×2.
    """
    if l == 0:
        return 0.0
    s   = 0.5
    LS  = 0.5 * (j*(j+1) - l*(l+1) - s*(s+1))
    r3  = Z**3 / (n**3 * l * (l + 0.5) * (l + 1))
    dE  = (ALPHA**2 / 2.0) * r3 * LS
    if not include_thomas:
        dE *= 2.0    # sem Thomas: não divide por 2
    return dE


# ============================================================================
# Evolução do spin por integração numérica (Rodrigues)
# ============================================================================

def evolve_spin(v, n_orbits=3, n_steps=3000):
    g      = gamma(v)
    omega  = 1.0
    T_orb  = 2 * np.pi
    dt     = n_orbits * T_orb / n_steps

    OL, OT, Otot = total_rate(v, omega)

    t_arr = np.linspace(0, n_orbits * T_orb, n_steps + 1)
    pos   = np.zeros((n_steps + 1, 2))
    spin  = np.zeros((n_steps + 1, 3))

    r_orb  = v
    pos[0] = [r_orb, 0.0]
    spin[0] = [1.0, 0.0, 0.0]

    for i in range(n_steps):
        t  = t_arr[i]
        vx = -v * np.sin(omega * t)
        vy =  v * np.cos(omega * t)
        ax = -omega * v * np.cos(omega * t)
        ay = -omega * v * np.sin(omega * t)

        # Ω_Thomas vectorial: −(γ−1)/v² (v × a)
        cross_z = vx * ay - vy * ax
        OmT = np.array([0.0, 0.0, -(g - 1) / (v**2) * cross_z])

        # Ω_Larmor no campo coulombiano (aponta em ẑ)
        OmL = np.array([0.0, 0.0, (2.0/2.0) * omega / g])

        Om  = OmL + OmT
        mag = np.linalg.norm(Om)
        s   = spin[i]
        if mag > 1e-12:
            ax_r  = Om / mag
            ang   = mag * dt
            s_new = (s * np.cos(ang)
                     + np.cross(ax_r, s) * np.sin(ang)
                     + ax_r * np.dot(ax_r, s) * (1 - np.cos(ang)))
        else:
            s_new = s.copy()
        spin[i+1] = s_new / (np.linalg.norm(s_new) + 1e-15)

        tn = t_arr[i+1]
        pos[i+1] = [r_orb * np.cos(omega * tn),
                    r_orb * np.sin(omega * tn)]

    return t_arr, pos, spin, OL, OT, Otot


# ============================================================================
# Painéis da figura estática
# ============================================================================

def _panel_geometry(ax, v):
    g     = gamma(v)
    theta = np.linspace(0, 2*np.pi, 300)
    r     = 1.0
    ax.plot(r*np.cos(theta), r*np.sin(theta), color="#D3D1C7", lw=1.5, ls="--")

    angles = np.linspace(0, np.pi, 6)
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(angles)))
    for phi, col in zip(angles, colors):
        px, py = r*np.cos(phi), r*np.sin(phi)
        ax.scatter([px], [py], s=50, color=col, zorder=5)
        sa = -(g - 1) * phi
        sx, sy = np.cos(sa), np.sin(sa)
        ax.annotate("", xy=(px+0.3*sx, py+0.3*sy), xytext=(px, py),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8, mutation_scale=11))

    delta = np.degrees(2*np.pi*(g - 1))
    ax.text(0.0, -1.6, f"δφ/órbita = {delta:.2f}°\n(γ={g:.3f}, v={v}c)",
            ha="center", fontsize=8.5, color="#378ADD")
    ax.scatter([0],[0], s=120, marker="+", linewidths=2, color="#EF9F27", zorder=6)
    ax.text(0.08, 0.06, "núcleo", fontsize=7.5, color="#EF9F27")
    ax.set_aspect("equal"); ax.set_xlim(-1.9,1.9); ax.set_ylim(-1.9,1.9)
    ax.set_title("Geometria: spin ao longo da órbita\n(setas = direção do spin)", fontsize=9)
    ax.set_xlabel("x", fontsize=9); ax.set_ylabel("y", fontsize=9)
    ax.grid(lw=0.3, alpha=0.3)


def _panel_frequencies(ax):
    v_arr = np.linspace(0.01, 0.995, 500)
    OL  = np.array([larmor_rate(v) for v in v_arr])
    OT  = np.array([thomas_rate(v) for v in v_arr])
    Ot  = OL + OT

    ax.plot(v_arr, OL,   color="#D85A30", lw=2.0, label=r"Larmor $\Omega_L = \omega/\gamma$")
    ax.plot(v_arr, -OT,  color="#378ADD", lw=2.0, label=r"Thomas $|\Omega_T|=(\gamma-1)\omega$")
    ax.plot(v_arr, Ot,   color="#1D9E75", lw=2.5, label=r"Total $\Omega_{SO}$")
    ax.axhline(0.5, color="#EF9F27", lw=1.2, ls="--",
               label=r"$v\to0$: $\frac{1}{2}\omega$ (fator ½ Thomas)")
    for vr in [0.3, 0.5, 0.8]:
        ax.axvline(vr, color="gray", lw=0.7, ls=":", alpha=0.5)

    ax.set_xlabel("v / c", fontsize=10)
    ax.set_ylabel(r"Frequência / $\omega_{orb}$", fontsize=10)
    ax.set_title("Frequências de precessão\nThomas cancela metade de Larmor → fator ½", fontsize=9)
    ax.set_xlim(0,1); ax.set_ylim(-0.1, 1.6)
    ax.legend(fontsize=8, loc="upper left"); ax.grid(lw=0.4, alpha=0.35)


def _panel_spin_orbit(ax, t_arr, pos, spin, v, OL, OT, Otot):
    T_orb  = 2*np.pi
    t_norm = t_arr / T_orb
    sa     = np.unwrap(np.arctan2(spin[:,1], spin[:,0]))
    phi_L  = OL * t_arr
    phi_T  = OT * t_arr
    phi_tot= (OL + OT) * t_arr

    ax.plot(t_norm, np.degrees(sa),     color="#1D9E75", lw=2.5, label="Numérico")
    ax.plot(t_norm, np.degrees(phi_tot),color="#1D9E75", lw=1.0, ls=":", alpha=0.6, label="Analítico")
    ax.plot(t_norm, np.degrees(phi_L),  color="#D85A30", lw=1.8, ls="--", alpha=0.8, label=r"Larmor $\Omega_L$")
    ax.plot(t_norm, np.degrees(phi_T),  color="#378ADD", lw=1.8, ls="--", alpha=0.8, label=r"Thomas $\Omega_T$")

    for n_orb in range(1, 4):
        ax.axvline(n_orb, color="gray", lw=0.8, ls=":", alpha=0.5)
        d = Otot * n_orb * T_orb
        ax.text(n_orb+0.04, np.degrees(d)+3, f"{np.degrees(d):.1f}°", fontsize=7.5, color="#1D9E75")

    ax.set_xlabel("Tempo (períodos orbitais)", fontsize=10)
    ax.set_ylabel("Ângulo do spin (°)", fontsize=10)
    ax.set_title(f"Spin em 3 órbitas  (v={v}c)\nΩ_L={OL:.3f}  Ω_T={OT:.3f}  Ω_tot={Otot:.3f}", fontsize=9)
    ax.legend(fontsize=8.5, loc="lower left"); ax.grid(lw=0.4, alpha=0.35)


def _panel_spin3d(ax, t_arr, pos, spin):
    r = np.max(np.sqrt(pos[:,0]**2+pos[:,1]**2))
    theta = np.linspace(0, 2*np.pi, 200)
    ax.plot(r*np.cos(theta), r*np.sin(theta),
            np.zeros(200), color="#D3D1C7", lw=0.8, alpha=0.4)

    step = max(1, len(t_arr)//35)
    for i in range(0, len(t_arr), step):
        p = np.append(pos[i], 0.0)
        s = spin[i] * 0.35
        col = plt.cm.viridis(i / len(t_arr))
        ax.quiver(p[0],p[1],p[2], s[0],s[1],s[2],
                  color=col, lw=1.0, arrow_length_ratio=0.3)

    ax.scatter([0],[0],[0], s=100, color="#EF9F27", zorder=8, marker="*")
    ax.set_xlabel("x",fontsize=7); ax.set_ylabel("y",fontsize=7); ax.set_zlabel("Sz",fontsize=7)
    ax.set_title("Spin 3D ao longo\nda órbita (cores=tempo)", fontsize=9)
    ax.view_init(elev=30, azim=-50)


def _panel_splitting(ax):
    """Splitting 2p em H: com vs sem Thomas."""
    estados = [("2p₁/₂", 1, 0.5), ("2p₃/₂", 1, 1.5)]
    n = 2
    y_with  = []
    y_without = []
    for lbl, l, j in estados:
        y_with.append(energy_spinorbit(n, l, j, include_thomas=True)  * HART_EV * 1e3)
        y_without.append(energy_spinorbit(n, l, j, include_thomas=False) * HART_EV * 1e3)

    colors = ["#378ADD", "#D85A30"]
    for i, (lbl, l, j) in enumerate(estados):
        col = colors[i]
        ax.hlines(y_with[i],    0.1, 0.45, color=col, lw=3.0)
        ax.hlines(y_without[i], 0.55, 0.9, color=col, lw=3.0, ls="--")
        ax.text(0.28, y_with[i],    lbl, ha="center", fontsize=8.5, color=col, va="bottom")
        ax.text(0.72, y_without[i], lbl, ha="center", fontsize=8.5, color=col, va="bottom")

    # Setas de splitting
    for x0, yw, label in [
        (0.28, y_with,    f"ΔE={y_with[1]-y_with[0]:.5f} meV\n(correto)"),
        (0.72, y_without, f"ΔE={y_without[1]-y_without[0]:.5f} meV\n(×2 errado)"),
    ]:
        ax.annotate("", xy=(x0, yw[1]), xytext=(x0, yw[0]),
                    arrowprops=dict(arrowstyle="<->", color="#1D9E75", lw=1.5))
        ax.text(x0+0.06, (yw[0]+yw[1])/2, label, fontsize=8, color="#1D9E75", va="center")

    ax.set_xlim(0, 1)
    ax.set_xticks([0.28, 0.72])
    ax.set_xticklabels(["Com Thomas\n(Dirac)", "Sem Thomas\n(Larmor puro)"], fontsize=8.5)
    ax.set_ylabel("ΔE spin-órbita (meV)", fontsize=9)
    ax.set_title("Splitting 2p em H\nFator ½ de Thomas é essencial!", fontsize=10)
    ax.grid(axis='y', lw=0.4, alpha=0.35)


def _panel_fine_levels(ax):
    """Níveis de spin-órbita para n=2,3 em μeV."""
    states = [
        (2,1,0.5,"2p₁/₂","#1D9E75"), (2,1,1.5,"2p₃/₂","#1D9E75"),
        (3,1,0.5,"3p₁/₂","#EF9F27"), (3,1,1.5,"3p₃/₂","#EF9F27"),
        (3,2,1.5,"3d₃/₂","#D85A30"), (3,2,2.5,"3d₅/₂","#D85A30"),
    ]
    xs = np.linspace(0.1, 0.9, len(states))
    w  = 0.055

    for x0, (n, l, j, lbl, col) in zip(xs, states):
        dw  = energy_spinorbit(n, l, j, include_thomas=True)  * HART_EV * 1e6
        dwo = energy_spinorbit(n, l, j, include_thomas=False) * HART_EV * 1e6
        ax.hlines(dw,  x0-w, x0,    color=col, lw=2.5)
        ax.hlines(dwo, x0,   x0+w,  color=col, lw=1.8, ls="--", alpha=0.7)
        ax.text(x0, dw + np.sign(dw)*0.005 + 0.003, lbl,
                ha="center", fontsize=7.5, color=col, va="bottom")

    ax.axhline(0, color="gray", lw=0.8, ls=":", alpha=0.5)
    ax.set_xticks([]); ax.set_xlabel("Estados", fontsize=9)
    ax.set_ylabel("ΔE spin-órbita (μeV)", fontsize=10)
    ax.set_title("Estrutura fina H: n=2,3  —  sólido=com Thomas  tracejado=sem Thomas", fontsize=9)
    ax.grid(axis="y", lw=0.4, alpha=0.35)
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0], color="gray", lw=2.5, label="Com Thomas (×½)"),
        Line2D([0],[0], color="gray", lw=1.8, ls="--", label="Sem Thomas"),
    ], fontsize=8.5, loc="upper right")


# ============================================================================
# Figura estática completa
# ============================================================================

def plot_all_static(args):
    v = args["v"]
    t_arr, pos, spin, OL, OT, Otot = evolve_spin(v, n_orbits=3)

    fig = plt.figure(figsize=(16, 9))
    fig.suptitle(
        f"Precessão de Thomas  (v={v}c, γ={gamma(v):.4f})  —  "
        f"Ω_L={OL:.4f}ω   Ω_T={OT:.4f}ω   Ω_tot={Otot:.4f}ω",
        fontsize=12
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.40)

    _panel_geometry(fig.add_subplot(gs[0,0]), v)
    _panel_frequencies(fig.add_subplot(gs[0,1]))
    _panel_spin_orbit(fig.add_subplot(gs[0,2]), t_arr, pos, spin, v, OL, OT, Otot)
    _panel_spin3d(fig.add_subplot(gs[1,0], projection='3d'), t_arr, pos, spin)
    _panel_splitting(fig.add_subplot(gs[1,1]))
    _panel_fine_levels(fig.add_subplot(gs[1,2]))

    plt.savefig("dirac_thomas.png", dpi=120, bbox_inches="tight")
    print("Figura salva: dirac_thomas.png")
    plt.show()


# ============================================================================
# Animação
# ============================================================================

def animate_thomas(args):
    v = args["v"]
    t_arr, pos, spin, OL, OT, Otot = evolve_spin(v, n_orbits=5, n_steps=5000)
    T_orb = 2 * np.pi
    r_orb = v

    fig = plt.figure(figsize=(13, 6))
    fig.suptitle(
        f"Precessão de Thomas  (v={v}c, γ={gamma(v):.4f})\n"
        f"Ω_Larmor={OL:.4f}ω   |   Ω_Thomas={OT:.4f}ω   |   Ω_total={Otot:.4f}ω",
        fontsize=11
    )
    gs  = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35)
    ax_xy  = fig.add_subplot(gs[0,0])
    ax_ang = fig.add_subplot(gs[0,1])

    # Fundo órbita
    th = np.linspace(0, 2*np.pi, 300)
    ax_xy.plot(r_orb*np.cos(th), r_orb*np.sin(th), color="#D3D1C7", lw=1.2, ls="--")
    ax_xy.scatter([0],[0], s=150, color="#EF9F27", zorder=8, marker="*")
    ax_xy.text(0.06, 0.06, "núcleo", fontsize=8, color="#EF9F27")
    ax_xy.set_aspect("equal")
    lim = r_orb * 1.7
    ax_xy.set_xlim(-lim,lim); ax_xy.set_ylim(-lim,lim)
    ax_xy.set_xlabel("x",fontsize=10); ax_xy.set_ylabel("y",fontsize=10)
    ax_xy.set_title("Órbita + spin", fontsize=10)
    ax_xy.grid(lw=0.4, alpha=0.35)

    dot_e,   = ax_xy.plot([], [], 'o', ms=10, color="#378ADD", zorder=6)
    trail,   = ax_xy.plot([], [], '-', lw=1.0, color="#378ADD", alpha=0.25)
    spin_arrow = ax_xy.annotate("", xy=(0,0), xytext=(0,0),
                                 arrowprops=dict(arrowstyle="-|>", color="#D85A30",
                                                 lw=2.5, mutation_scale=18))
    vel_arrow  = ax_xy.annotate("", xy=(0,0), xytext=(0,0),
                                 arrowprops=dict(arrowstyle="-|>", color="#1D9E75",
                                                 lw=1.5, mutation_scale=12))
    txt_ang = ax_xy.text(0.03, 0.96, "", transform=ax_xy.transAxes,
                         fontsize=9, va="top", color="#D85A30")

    # Painel ângulo
    sa_all = np.unwrap(np.arctan2(spin[:,1], spin[:,0]))
    t_norm_all = t_arr / T_orb
    ax_ang.plot(t_norm_all, np.degrees(sa_all),
                color="#D85A30", lw=1.0, alpha=0.25)
    ax_ang.plot(t_norm_all, np.degrees(t_arr),
                color="#1D9E75", lw=1.0, ls="--", alpha=0.4, label="Ângulo orbital (ref)")
    ax_ang.set_xlabel("Tempo (períodos)", fontsize=10)
    ax_ang.set_ylabel("Ângulo spin (°)", fontsize=10)
    ax_ang.set_title("Defasagem spin vs orbital\n(afastamento = Thomas!)", fontsize=10)
    ax_ang.grid(lw=0.4, alpha=0.35)
    ax_ang.legend(fontsize=8)

    dot_ang,  = ax_ang.plot([], [], 'o', ms=7, color="#D85A30", zorder=5)
    line_past, = ax_ang.plot([], [], '-', lw=2, color="#D85A30")
    txt_frame = fig.text(0.5, 0.01, "t = 0.00", ha="center", fontsize=9, color="gray")
    fig.tight_layout(rect=[0, 0.04, 1, 0.90])

    N_frames = 300
    trail_len = 60

    def update(frame):
        idx = min(int(frame * len(t_arr) / N_frames), len(t_arr)-1)
        p   = pos[idx]
        s   = spin[idx]
        t   = t_arr[idx]
        omega = 1.0

        dot_e.set_data([p[0]], [p[1]])
        i0 = max(0, idx - trail_len)
        trail.set_data(pos[i0:idx+1,0], pos[i0:idx+1,1])

        ss = r_orb * 0.55
        spin_arrow.arrow_patch.set_positions((p[0],p[1]), (p[0]+s[0]*ss, p[1]+s[1]*ss))

        vx = -v*np.sin(omega*t); vy = v*np.cos(omega*t)
        vs = r_orb * 0.35
        vel_arrow.arrow_patch.set_positions((p[0],p[1]), (p[0]+vx*vs, p[1]+vy*vs))

        ang = np.degrees(np.arctan2(s[1],s[0]))
        txt_ang.set_text(f"spin: {ang:.1f}°\nórbita: {np.degrees(omega*t)%360:.1f}°")

        dot_ang.set_data([t/T_orb], [np.degrees(sa_all[idx])])
        line_past.set_data(t_norm_all[:idx+1], np.degrees(sa_all[:idx+1]))
        txt_frame.set_text(
            f"t = {t/T_orb:.2f} órbitas  |  "
            f"defasagem = {np.degrees(sa_all[idx] - omega*t):.1f}°"
        )

    ani = animation.FuncAnimation(fig, update, frames=N_frames, interval=40, blit=False)
    if args.get("save"):
        fname = "dirac_thomas.gif"
        print(f"Salvando '{fname}'...")
        ani.save(fname, writer=animation.PillowWriter(fps=25))
        print(f"Salvo: {fname}")
    else:
        plt.show()
    return ani


# ============================================================================
# Sumário + CLI
# ============================================================================

def print_summary(args):
    v = args["v"]
    g = gamma(v)
    OL, OT, Otot = total_rate(v)

    print(f"\n{'='*64}")
    print("  Precessão de Thomas")
    print(f"{'='*64}")
    print(f"  v = {v}c   γ = {g:.6f}   γ−1 = {g-1:.6f}")
    print(f"  Ω_Larmor  = {OL:+.6f} ω_orb")
    print(f"  Ω_Thomas  = {OT:+.6f} ω_orb   (sinal −: oposto ao orbital)")
    print(f"  Ω_total   = {Otot:+.6f} ω_orb")
    print(f"  Defasagem/órbita = {np.degrees(Otot*2*np.pi):.4f}°")
    print()
    print("  Splitting 2p em H (perturbativo):")
    for l,j,lbl in [(1,0.5,"2p₁/₂"),(1,1.5,"2p₃/₂")]:
        dw  = energy_spinorbit(2,l,j,include_thomas=True) *HART_EV*1e3
        dwo = energy_spinorbit(2,l,j,include_thomas=False)*HART_EV*1e3
        print(f"    {lbl}: com Thomas={dw:.6f} meV  sem Thomas={dwo:.6f} meV")
    dE_w  = (energy_spinorbit(2,1,1.5,include_thomas=True)
            -energy_spinorbit(2,1,0.5,include_thomas=True))*HART_EV*1e3
    dE_wo = (energy_spinorbit(2,1,1.5,include_thomas=False)
            -energy_spinorbit(2,1,0.5,include_thomas=False))*HART_EV*1e3
    print(f"  Splitting com Thomas = {dE_w:.6f} meV  (correto)")
    print(f"  Splitting sem Thomas = {dE_wo:.6f} meV  (×2 errado, ref NIST ~0.045 meV)")
    print(f"{'='*64}\n")


def parse_args():
    p = argparse.ArgumentParser(
        description="Precessão de Thomas — Relatividade e Spin-Órbita",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--v",      type=float, default=0.5, help="Velocidade v/c")
    p.add_argument("--module", type=int,   default=0,   help="0=tudo 1=estático 5=animação")
    p.add_argument("--save",   action="store_true",     help="Salvar animação GIF")
    return vars(p.parse_args())


if __name__ == "__main__":
    args = parse_args()
    print_summary(args)
    mod = args["module"]
    if mod == 0 or mod == 1:
        plot_all_static(args)
    if mod == 0 or mod == 5:
        animate_thomas(args)
