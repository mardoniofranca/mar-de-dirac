"""
Grafeno 2D — Férmions de Dirac Sem Massa
=========================================
Módulo 1 — Estrutura de banda: rede hexagonal + pontos de Dirac K e K'
Módulo 2 — Cone de Dirac: dispersão linear E = ±ħv_F|k|  (sem gap, sem massa)
Módulo 3 — Propagação de pacote: velocidade constante v_F ≈ c/300, sem dispersão
Módulo 4 — Pseudospin e quiralidade: helicidade travada, tunelamento de Klein perfeito
Módulo 5 — Níveis de Landau em grafeno: E_n = v_F√(2nħeB)  ∝ √(nB), inclui n=0

Diferença fundamental em relação ao elétron livre:
  - Elétron livre:  E = ±√(k²+m²)   → dispersão quadrática para k<<m
  - Grafeno:        E = ±v_F|k|      → linear para TODO k, m_efetivo = 0

Sistema de unidades: ħ = e = 1, v_F = 1 (comprimento em unidades de a/√3)

Uso:
    python dirac_graphene.py
    python dirac_graphene.py --module 2    # só cone de Dirac
    python dirac_graphene.py --B 3.0       # campo magnético para Landau
    python dirac_graphene.py --save        # salva animação
    python dirac_graphene.py --help
"""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from matplotlib.patches import RegularPolygon, FancyArrowPatch
from matplotlib.collections import PatchCollection
import argparse
import platform
import warnings
warnings.filterwarnings("ignore")

if platform.system() == "Darwin":
    matplotlib.use("MacOSX")

# ============================================================================
# Constantes do grafeno
# ============================================================================
VF      = 1.0           # velocidade de Fermi (unidades v_F=1)
VF_SI   = 1.0e6         # v_F em m/s (real)
C_SI    = 3.0e8         # c em m/s
VF_FRAC = VF_SI / C_SI  # v_F/c ≈ 1/300

# Rede hexagonal: vetores primitivos
A_LAT  = 1.0                          # parâmetro de rede (u.a.)
A1     = A_LAT * np.array([1.0, 0.0])
A2     = A_LAT * np.array([0.5, np.sqrt(3)/2])

# Pontos de alta simetria na zona de Brillouin hexagonal
# K  = (2π/a)(1/3, 1/√3)  e  K' = -K
K_PT   = (2*np.pi/A_LAT) * np.array([1/3, 1/np.sqrt(3)])
KP_PT  = -K_PT


# ============================================================================
# MÓDULO 1 — Rede Hexagonal e Zona de Brillouin
# ============================================================================

def plot_lattice_and_bz(ax_lat, ax_bz):
    """Rede real do grafeno (hexagonal) e zona de Brillouin hexagonal."""

    # ── Rede real: sítios A (azul) e B (laranja) ──
    n_cells = 4
    sA, sB = [], []
    for i in range(-n_cells, n_cells + 1):
        for j in range(-n_cells, n_cells + 1):
            center = i * A1 + j * A2
            sA.append(center)
            sB.append(center + np.array([0.5, np.sqrt(3)/6]))

    sA = np.array(sA); sB = np.array(sB)
    mask = (np.abs(sA[:,0]) < 3.5) & (np.abs(sA[:,1]) < 3.0)
    maskB = (np.abs(sB[:,0]) < 3.5) & (np.abs(sB[:,1]) < 3.0)

    # Linhas de ligação C-C
    for pos in sA[mask]:
        for delta in [np.array([0.5, np.sqrt(3)/6]),
                      np.array([-0.5, np.sqrt(3)/6]),
                      np.array([0.0, -np.sqrt(3)/3])]:
            nb = pos + delta
            ax_lat.plot([pos[0], nb[0]], [pos[1], nb[1]],
                        color="gray", lw=0.8, alpha=0.6, zorder=1)

    ax_lat.scatter(sA[mask, 0], sA[mask, 1], s=60, color="#378ADD",
                   zorder=3, label="Sítio A")
    ax_lat.scatter(sB[maskB, 0], sB[maskB, 1], s=60, color="#D85A30",
                   zorder=3, label="Sítio B")

    # Destaca célula unitária
    hex_verts = []
    for angle in range(0, 360, 60):
        r = A_LAT / np.sqrt(3)
        hex_verts.append([r * np.cos(np.radians(angle + 30)),
                          r * np.sin(np.radians(angle + 30))])
    hex_verts = np.array(hex_verts)
    ax_lat.fill(hex_verts[:, 0], hex_verts[:, 1],
                alpha=0.10, color="#1D9E75", zorder=0)
    ax_lat.plot(np.append(hex_verts[:, 0], hex_verts[0, 0]),
                np.append(hex_verts[:, 1], hex_verts[0, 1]),
                color="#1D9E75", lw=1.2, ls="--")

    ax_lat.set_aspect("equal")
    ax_lat.set_xlim(-3.5, 3.5); ax_lat.set_ylim(-3.0, 3.0)
    ax_lat.set_xlabel("x (a)", fontsize=9); ax_lat.set_ylabel("y (a)", fontsize=9)
    ax_lat.set_title("Rede hexagonal do grafeno\n(2 átomos por célula: A e B)", fontsize=9)
    ax_lat.legend(fontsize=8, loc="upper right")
    ax_lat.grid(lw=0.3, alpha=0.3)

    # ── Zona de Brillouin hexagonal ──
    # Vetores da rede recíproca
    b1 = (2*np.pi/A_LAT) * np.array([1.0, -1/np.sqrt(3)])
    b2 = (2*np.pi/A_LAT) * np.array([0.0,  2/np.sqrt(3)])

    # Vértices da ZB hexagonal
    bz_verts = []
    for i in range(6):
        angle = np.radians(30 + 60 * i)
        r_bz = 4*np.pi / (3 * A_LAT)
        bz_verts.append([r_bz * np.cos(angle), r_bz * np.sin(angle)])
    bz_verts = np.array(bz_verts)

    ax_bz.fill(np.append(bz_verts[:,0], bz_verts[0,0]),
               np.append(bz_verts[:,1], bz_verts[0,1]),
               alpha=0.08, color="#378ADD")
    ax_bz.plot(np.append(bz_verts[:,0], bz_verts[0,0]),
               np.append(bz_verts[:,1], bz_verts[0,1]),
               color="#378ADD", lw=1.5)

    # Pontos K e K'
    K_pts  = bz_verts[0::2]   # alternados
    Kp_pts = bz_verts[1::2]
    ax_bz.scatter(K_pts[:,0],  K_pts[:,1],  s=80, color="#D85A30",
                  zorder=5, label="K  (cone de Dirac)")
    ax_bz.scatter(Kp_pts[:,0], Kp_pts[:,1], s=80, color="#1D9E75",
                  zorder=5, marker="s", label="K' (cone de Dirac)")

    # Ponto Γ e M
    ax_bz.scatter([0], [0], s=60, color="gray", zorder=5)
    ax_bz.text(0.05, 0.05, "Γ", fontsize=9, color="gray")

    for kp in K_pts[:1]:
        ax_bz.annotate("K", xy=kp, xytext=kp + np.array([0.2, 0.2]),
                       fontsize=9, color="#D85A30")
    for kp in Kp_pts[:1]:
        ax_bz.annotate("K'", xy=kp, xytext=kp + np.array([0.2, -0.4]),
                       fontsize=9, color="#1D9E75")

    ax_bz.set_aspect("equal")
    lim = 5.5
    ax_bz.set_xlim(-lim, lim); ax_bz.set_ylim(-lim, lim)
    ax_bz.set_xlabel("kₓ (1/a)", fontsize=9); ax_bz.set_ylabel("kᵧ (1/a)", fontsize=9)
    ax_bz.set_title("Zona de Brillouin — 6 pontos K/K'\n(pontos de Dirac: gap = 0)", fontsize=9)
    ax_bz.legend(fontsize=8, loc="upper right")
    ax_bz.grid(lw=0.3, alpha=0.3)


# ============================================================================
# MÓDULO 2 — Cone de Dirac: estrutura de banda 2D
# ============================================================================

def dirac_cone_energy(kx, ky, vf=VF):
    """E = ±v_F √(kx²+ky²)  — dispersão linear, sem massa."""
    k = np.sqrt(kx**2 + ky**2)
    return vf * k


def plot_dirac_cone(ax3d, ax_bands):
    """Cone de Dirac 3D e estrutura de banda ao longo de K-Γ-K'."""

    # ── Cone 3D ──
    k_range = np.linspace(-2.5, 2.5, 80)
    KX, KY  = np.meshgrid(k_range, k_range)
    E_pos   =  dirac_cone_energy(KX, KY)
    E_neg   = -dirac_cone_energy(KX, KY)

    ax3d.plot_surface(KX, KY, E_pos, cmap="Blues",   alpha=0.75, linewidth=0)
    ax3d.plot_surface(KX, KY, E_neg, cmap="Oranges", alpha=0.75, linewidth=0)

    # Linha de contorno no nível de Fermi E=0
    ax3d.contour(KX, KY, E_pos, levels=[0], colors=["#1D9E75"], linewidths=2)
    ax3d.contour(KX, KY, E_neg, levels=[0], colors=["#1D9E75"], linewidths=2)

    # Ponto de Dirac
    ax3d.scatter([0], [0], [0], s=100, color="#D85A30", zorder=10)
    ax3d.text(0.2, 0.2, 0.3, "Ponto\nde Dirac", fontsize=8, color="#D85A30")

    ax3d.set_xlabel("kₓ", fontsize=8); ax3d.set_ylabel("kᵧ", fontsize=8)
    ax3d.set_zlabel("E / ħv_F", fontsize=8)
    ax3d.set_title("Cone de Dirac: E = ±v_F|k|\n(grafeno, m_eff = 0)", fontsize=9)
    ax3d.view_init(elev=25, azim=-60)

    # ── Estrutura de banda ao longo de k ──
    k_line = np.linspace(-3.0, 3.0, 400)
    E_p    = VF * np.abs(k_line)
    E_n    = -VF * np.abs(k_line)

    # Comparação com elétron parabólico
    m_eff  = 0.5  # massa efetiva de elétron típico em semicondutor
    E_par  = k_line**2 / (2 * m_eff)

    ax_bands.plot(k_line, E_p,  color="#378ADD", lw=2.5, label=r"Grafeno: $E=+v_F|k|$ (linear)")
    ax_bands.plot(k_line, E_n,  color="#D85A30", lw=2.5, label=r"Grafeno: $E=-v_F|k|$ (linear)")
    ax_bands.plot(k_line, E_par, color="#1D9E75", lw=1.8, ls="--", alpha=0.7,
                  label=r"Semicond.: $E=k^2/2m^*$ (parabólico)")

    ax_bands.axhline(0, color="gray", lw=0.8, ls="--", alpha=0.6)
    ax_bands.axvline(0, color="#D85A30", lw=0.8, ls=":", alpha=0.6)

    # Nível de Fermi e anotações
    ax_bands.fill_between(k_line, E_n, -4, alpha=0.10, color="#D85A30",
                          label="Banda de valência\n(preenchida)")
    ax_bands.text(0.08, 0.15, "Ponto de\nDirac (K)", fontsize=8,
                  color="#D85A30", transform=ax_bands.transAxes)

    ax_bands.set_xlabel("k (em torno de K)", fontsize=10)
    ax_bands.set_ylabel("Energia E / ħv_F", fontsize=10)
    ax_bands.set_title("Estrutura de banda — grafeno vs semicondutor", fontsize=10)
    ax_bands.set_ylim(-4, 4); ax_bands.set_xlim(-3, 3)
    ax_bands.legend(fontsize=8, loc="upper center")
    ax_bands.grid(lw=0.4, alpha=0.35)


# ============================================================================
# MÓDULO 3 — Propagação sem dispersão (m=0)
# ============================================================================

def init_spinor_massless(xs, ys, dx, dy, kx0, ky0, sigma):
    """
    Pacote gaussiano 2D para férmion de Dirac sem massa.
    O espinor é (ψ_A, ψ_B) = subrede A e subrede B.

    Para modo de energia positiva com momento k = (kx, ky):
        ψ_A = f(r) · e^{i(kx·x + ky·y)}
        ψ_B = e^{iφ_k} · f(r) · e^{i(kx·x + ky·y)}
    onde φ_k = arctan(ky/kx)  (fase quiral)
    """
    phi_k = np.arctan2(ky0, kx0)
    env   = np.exp(-((xs**2 + ys**2) / (2 * sigma**2)))
    env  /= np.sqrt(np.sum(env**2) * dx * dy)

    psi_A = env * np.exp(1j * (kx0 * xs + ky0 * ys))
    psi_B = np.exp(1j * phi_k) * psi_A

    return psi_A, psi_B


def propagate_massless_1d(xs, dx, k0, sigma, n_steps, dt):
    """
    Evolução 1D do férmion sem massa via split-operator.
    Para m=0: ω(k) = v_F|k|, propagador exato.
    """
    N  = len(xs)
    ks = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    omega = VF * np.abs(ks)

    # Spinor inicial (1D: apenas componente A)
    env  = np.exp(-0.25 * (xs / sigma)**2) / (np.pi * sigma**2)**0.25
    psi_A = env * np.exp(1j * k0 * xs)
    phi_k = np.sign(k0) * np.pi / 4 if k0 != 0 else 0.0
    psi_B = np.exp(1j * phi_k) * psi_A

    # Normaliza
    norm = np.sqrt(np.sum((np.abs(psi_A)**2 + np.abs(psi_B)**2) * dx))
    psi_A /= norm; psi_B /= norm

    # Propagador m=0 em 1D
    cos_w  = np.cos(omega * dt)
    sinc_w = np.where(np.abs(ks) > 1e-12, np.sin(omega * dt) / (omega + 1e-15), dt)
    k_sign = np.sign(ks + 1e-15)

    P11 =  cos_w
    P12 = -1j * k_sign * np.sin(omega * dt)
    P21 = -1j * k_sign * np.sin(omega * dt)
    P22 =  cos_w

    snapshots = []
    t_arr     = []

    for step in range(n_steps + 1):
        if step % (n_steps // 10) == 0:
            dens = (np.abs(psi_A)**2 + np.abs(psi_B)**2).real
            snapshots.append(dens.copy())
            t_arr.append(step * dt)

        if step < n_steps:
            phi_A = np.fft.fft(psi_A)
            phi_B = np.fft.fft(psi_B)
            psi_A = np.fft.ifft(P11 * phi_A + P12 * phi_B)
            psi_B = np.fft.ifft(P21 * phi_A + P22 * phi_B)

    return snapshots, t_arr


def plot_propagation(ax, m=1.0):
    """
    Compara propagação de pacote: grafeno (m=0) vs elétron livre (m=1).
    Para m=0: sem dispersão, velocidade constante v_F.
    Para m≠0: dispersão, o pacote alarga.
    """
    N     = 1024
    L     = 80.0
    xs    = np.linspace(-L/2, L/2, N, endpoint=False)
    dx    = xs[1] - xs[0]
    ks    = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    k0    = 3.0
    sigma = 5.0
    dt    = 0.02
    n_steps = 500

    colors_t = plt.cm.viridis(np.linspace(0.2, 1.0, 5))

    for mi, (mass, lstyle, label) in enumerate([
        (0.0, "-",  "Grafeno m=0 (sem dispersão)"),
        (1.0, "--", "Elétron livre m=1 (com dispersão)"),
    ]):
        omega = VF * np.sqrt(ks**2 + mass**2) if mass > 0 else VF * np.abs(ks)
        cos_w  = np.cos(omega * dt)
        sinc_w = np.where(np.abs(ks) > 1e-12, np.sin(omega * dt) / (omega + 1e-15), dt)
        k_sign = ks / (np.abs(ks) + 1e-12)

        if mass > 0:
            P11 =  cos_w - 1j * mass * sinc_w
            P12 = -1j * ks * sinc_w
            P21 = -1j * ks * sinc_w
            P22 =  cos_w + 1j * mass * sinc_w
        else:
            P11 =  cos_w
            P12 = -1j * k_sign * np.sin(omega * dt)
            P21 = -1j * k_sign * np.sin(omega * dt)
            P22 =  cos_w

        env   = np.exp(-0.25 * (xs / sigma)**2) / (np.pi * sigma**2)**0.25
        psi1  = env * np.exp(1j * k0 * xs)
        psi2  = (k0 / (np.sqrt(k0**2 + mass**2 + 1e-12) + mass + 1e-12)) * psi1
        norm  = np.sqrt(np.sum((np.abs(psi1)**2 + np.abs(psi2)**2) * dx))
        psi1 /= norm; psi2 /= norm

        snap_steps = [0, 100, 200, 300, 400, 500]
        si = 0
        for step in range(n_steps + 1):
            if step == snap_steps[si]:
                dens = (np.abs(psi1)**2 + np.abs(psi2)**2).real
                t_now = step * dt
                offset = mi * 0.025
                lbl = f"{label}  t={t_now:.1f}" if si == 0 else f"t={t_now:.1f}"
                col = colors_t[si % len(colors_t)]
                ax.plot(xs, dens + offset, color=col if mass == 0 else col,
                        lw=1.5, ls=lstyle, alpha=0.85,
                        label=lbl if si <= 1 else "")
                si += 1
                if si >= len(snap_steps):
                    break
            if step < n_steps:
                phi1 = np.fft.fft(psi1); phi2 = np.fft.fft(psi2)
                psi1 = np.fft.ifft(P11*phi1 + P12*phi2)
                psi2 = np.fft.ifft(P21*phi1 + P22*phi2)

    ax.set_xlabel("x", fontsize=10)
    ax.set_ylabel(r"$|\psi|^2$ + offset", fontsize=10)
    ax.set_title("Propagação: grafeno (sem dispersão) vs elétron livre", fontsize=10)
    ax.set_xlim(-L/2, L/2)
    ax.legend(fontsize=7.5, loc="upper right", ncol=2)
    ax.grid(lw=0.4, alpha=0.35)


# ============================================================================
# MÓDULO 4 — Pseudospin, Quiralidade e Tunelamento de Klein Perfeito
# ============================================================================

def klein_transmission_graphene(k, V0, d):
    """
    Transmissão de Klein para grafeno (m=0).
    Resultado analítico de Katsnelson et al. 2006:
        T = cos²(θ) / (1 − sin²(θ)sin²(qd))
    onde θ = arctan(ky/kx) e q = √((E-V0)²/v_F² − ky²)
    Para incidência normal (ky=0, θ=0): T = 1 sempre (Klein perfeito!)
    """
    # Incidência oblíqua: ky = k·sin(θ), kx = k·cos(θ)
    theta_arr = np.linspace(-np.pi/2 + 0.01, np.pi/2 - 0.01, 300)
    E = VF * k
    T_arr = []
    for theta in theta_arr:
        ky = k * np.sin(theta)
        kx = k * np.cos(theta)
        Ev = E - V0
        arg = (Ev/VF)**2 - ky**2
        if arg <= 0:
            T_arr.append(0.0)
            continue
        q = np.sqrt(arg)
        cos2_theta = np.cos(theta)**2
        sin2_phi   = (ky / (np.abs(Ev/VF) + 1e-12))**2
        denom = 1.0 - (1 - cos2_theta) * np.sin(q * d)**2
        T = cos2_theta / (denom + 1e-12) if denom > 0 else 0.0
        T_arr.append(np.clip(T, 0, 1))
    return theta_arr, np.array(T_arr)


def plot_pseudospin(ax_ps, ax_klein):
    """Diagrama de pseudospin e transmissão de Klein em grafeno."""

    # ── Pseudospin: helicidade travada ──
    k_angles = np.linspace(0, 2*np.pi, 12, endpoint=False)
    for angle in k_angles:
        kx = np.cos(angle); ky = np.sin(angle)
        # Pseudospin aponta na direção de k (valência: antiparalelo)
        color = "#378ADD" if True else "#D85A30"
        ax_ps.annotate("", xy=(kx * 1.5, ky * 1.5),
                       xytext=(kx * 0.9, ky * 0.9),
                       arrowprops=dict(arrowstyle="-|>", color="#378ADD",
                                       lw=1.2, mutation_scale=10))
        ax_ps.annotate("", xy=(-kx * 1.5, -ky * 1.5),
                       xytext=(-kx * 0.9, -ky * 0.9),
                       arrowprops=dict(arrowstyle="-|>", color="#D85A30",
                                       lw=1.0, mutation_scale=10))

    circle = plt.Circle((0, 0), 0.85, fill=False, color="gray", ls="--", lw=0.8)
    ax_ps.add_patch(circle)
    ax_ps.scatter([0], [0], s=40, color="black", zorder=5)
    ax_ps.text(0.0, -0.25, "K", ha="center", fontsize=9, color="black")

    ax_ps.set_aspect("equal")
    ax_ps.set_xlim(-2.2, 2.2); ax_ps.set_ylim(-2.2, 2.2)
    ax_ps.set_title("Pseudospin no ponto K\n(azul=condução, laranja=valência)", fontsize=9)
    ax_ps.text(1.6, 0.1, r"$\vec{k}$", fontsize=10, color="#378ADD")
    ax_ps.text(-1.8, 0.1, r"$-\vec{k}$", fontsize=10, color="#D85A30")
    ax_ps.set_xlabel("kₓ", fontsize=9); ax_ps.set_ylabel("kᵧ", fontsize=9)
    ax_ps.grid(lw=0.3, alpha=0.3)

    # ── Transmissão de Klein vs ângulo de incidência ──
    k0   = 2.0
    V0   = 1.5 * VF * k0
    d    = 3.0
    theta_arr, T_graphene = klein_transmission_graphene(k0, V0, d)

    # Elétron normal (de Schrödinger): T cai para grandes ângulos
    T_schrod = np.cos(theta_arr)**4  # aproximação simples

    ax_klein.plot(np.degrees(theta_arr), T_graphene,
                  color="#378ADD", lw=2.5, label="Grafeno (Dirac m=0)")
    ax_klein.plot(np.degrees(theta_arr), T_schrod,
                  color="#D85A30", lw=2.0, ls="--", label="Semicondutor (parabólico)")

    ax_klein.axvline(0, color="#1D9E75", lw=1.2, ls=":",
                     label="θ=0: T=1 sempre (Klein perfeito!)")
    ax_klein.fill_between(np.degrees(theta_arr), T_graphene, 0,
                          alpha=0.10, color="#378ADD")

    ax_klein.set_xlabel("Ângulo de incidência θ (graus)", fontsize=10)
    ax_klein.set_ylabel("Transmissão T(θ)", fontsize=10)
    ax_klein.set_title("Tunelamento de Klein em grafeno\n(incidência normal: T=1 perfeito)", fontsize=10)
    ax_klein.set_xlim(-90, 90)
    ax_klein.set_ylim(-0.05, 1.15)
    ax_klein.legend(fontsize=8.5)
    ax_klein.grid(lw=0.4, alpha=0.35)


# ============================================================================
# MÓDULO 5 — Níveis de Landau em Grafeno
# ============================================================================

def landau_graphene(n, B, vf=VF):
    """
    Nível n de Landau em grafeno (m=0):
        E_n = ±v_F √(2n·ħeB)    unidades naturais (ħ=e=1)
    Diferença crucial vs elétron livre (m≠0):
        - Grafeno:       E_n ∝ √(nB)   → raiz de n E de B
        - Elétron livre: E_n ∝ (n+½)B  → linear em n e B
    O nível n=0 tem E=0 exato para qualquer B!
    """
    return vf * np.sqrt(2.0 * n * B)


def plot_landau_graphene(ax_ll, ax_comp):
    """Níveis de Landau de grafeno vs elétron livre."""

    B_arr  = np.linspace(0.1, 5.0, 300)
    n_show = [0, 1, 2, 3, 4, 5]
    cmap   = plt.cm.viridis

    # ── E_n vs B ──
    for n in n_show:
        col = cmap(n / max(n_show))
        E_g = np.array([landau_graphene(n, B) for B in B_arr])
        ax_ll.plot(B_arr, E_g, color=col, lw=2.0 if n > 0 else 2.5,
                   label=f"n={n}")
        ax_ll.plot(B_arr, -E_g, color=col, lw=1.5, ls="--", alpha=0.6)

    ax_ll.axhline(0, color="#D85A30", lw=2.0, label="n=0 (E=0 SEMPRE!)")
    ax_ll.set_xlabel("Campo B", fontsize=10)
    ax_ll.set_ylabel(r"$E_n = \pm v_F\sqrt{2nB}$", fontsize=10)
    ax_ll.set_title("Níveis de Landau em grafeno\nE_n ∝ √(nB)  —  n=0 fixo em E=0", fontsize=10)
    ax_ll.legend(fontsize=8, loc="upper left", ncol=2)
    ax_ll.grid(lw=0.4, alpha=0.35)

    # ── Comparação: grafeno vs elétron livre (n=0 a 4, B=1) ──
    B_fix  = 1.0
    ns     = np.arange(0, 8)
    m_free = 1.0
    omega_c = B_fix / m_free  # frequência ciclotron

    E_graph = np.array([landau_graphene(n, B_fix) for n in ns])
    E_free  = np.array([omega_c * (n + 0.5) for n in ns])        # Schrödinger
    E_relf  = np.array([np.sqrt(1 + 2*n*B_fix) - 1 for n in ns]) # Dirac m≠0

    x = np.arange(len(ns))
    w = 0.28
    ax_comp.bar(x - w, E_graph, w, color="#378ADD",   label="Grafeno (m=0): √(2nB)")
    ax_comp.bar(x,     E_relf,  w, color="#1D9E75",   label="Dirac (m≠0): √(1+2nB)−1", alpha=0.85)
    ax_comp.bar(x + w, E_free,  w, color="#D3D1C7",   label="Schrödinger: ωc(n+½)")

    # Destaca n=0
    ax_comp.bar([0 - w], [E_graph[0]], w, color="#D85A30",
                label="n=0 grafeno: E=0 !", zorder=5)

    ax_comp.set_xticks(x)
    ax_comp.set_xticklabels([f"n={n}" for n in ns], fontsize=8.5)
    ax_comp.set_ylabel("Energia cinética", fontsize=10)
    ax_comp.set_title(f"Comparação de 3 teorias (B={B_fix})\nGrafeno vs Dirac m≠0 vs Schrödinger", fontsize=10)
    ax_comp.legend(fontsize=8)
    ax_comp.grid(axis='y', lw=0.4, alpha=0.35)


# ============================================================================
# ANIMAÇÃO — Pacote 1D em grafeno sem dispersão vs com massa
# ============================================================================

def animate_graphene(args):
    N     = 1024
    L     = 100.0
    xs    = np.linspace(-L/2, L/2, N, endpoint=False)
    dx    = xs[1] - xs[0]
    ks    = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    k0    = 3.0
    sigma = 6.0
    dt    = 0.03
    x0    = -L/4

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(
        "Grafeno vs Elétron Livre — Propagação sem dispersão (m=0) vs com dispersão (m≠0)",
        fontsize=12
    )

    states = {}
    lines  = {}
    colors = {"grafeno": "#378ADD", "livre": "#D85A30"}

    for name, mass in [("grafeno", 0.0), ("livre", 1.0)]:
        omega = VF * np.abs(ks) if mass == 0 else np.sqrt(ks**2 + mass**2)
        cos_w  = np.cos(omega * dt)
        sin_w  = np.sin(omega * dt)
        sinc_w = np.where(np.abs(omega) > 1e-12, sin_w / (omega + 1e-15), dt)
        k_sign = ks / (np.abs(ks) + 1e-12)

        if mass == 0:
            P = (cos_w, -1j * k_sign * sin_w, -1j * k_sign * sin_w, cos_w)
        else:
            P = (cos_w - 1j*mass*sinc_w, -1j*ks*sinc_w,
                 -1j*ks*sinc_w,           cos_w + 1j*mass*sinc_w)

        env  = np.exp(-0.25 * ((xs - x0) / sigma)**2) / (np.pi * sigma**2)**0.25
        p1   = env * np.exp(1j * k0 * xs)
        p2   = (k0 / (np.sqrt(k0**2 + mass**2 + 1e-12) + mass + 1e-12)) * p1
        norm = np.sqrt(np.sum((np.abs(p1)**2 + np.abs(p2)**2) * dx))
        p1  /= norm; p2 /= norm

        states[name] = {"psi1": p1, "psi2": p2, "P": P, "mass": mass}

    ax_g, ax_l = axes
    dens_g0 = (np.abs(states["grafeno"]["psi1"])**2 + np.abs(states["grafeno"]["psi2"])**2).real
    dens_l0 = (np.abs(states["livre"]["psi1"])**2   + np.abs(states["livre"]["psi2"])**2).real
    ymax = max(dens_g0.max(), dens_l0.max()) * 2.2

    lines["grafeno"], = ax_g.plot(xs, dens_g0, color=colors["grafeno"], lw=2,
                                  label=r"Grafeno $m=0$: $|\psi|^2$")
    lines["livre"],   = ax_l.plot(xs, dens_l0, color=colors["livre"],   lw=2,
                                  label=r"Elétron livre $m=1$: $|\psi|^2$")

    for ax, name in [(ax_g, "grafeno"), (ax_l, "livre")]:
        ax.set_ylim(0, ymax)
        ax.set_ylabel(r"$|\psi|^2$", fontsize=11)
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(lw=0.4, alpha=0.4)

    ax_g.set_title("Grafeno (m=0): forma do pacote preservada — v_F constante", fontsize=10)
    ax_l.set_title("Elétron livre (m≠0): pacote se dispersa — velocidades de grupo diferentes", fontsize=10)
    ax_l.set_xlabel("x (unidades naturais)", fontsize=11)

    txt_g  = ax_g.text(0.02, 0.88, "t = 0.00", transform=ax_g.transAxes, fontsize=9)
    txt_l  = ax_l.text(0.02, 0.88, "t = 0.00", transform=ax_l.transAxes, fontsize=9)
    txt_w_g = ax_g.text(0.98, 0.88, "", transform=ax_g.transAxes, fontsize=9,
                        ha="right", color=colors["grafeno"])
    txt_w_l = ax_l.text(0.98, 0.88, "", transform=ax_l.transAxes, fontsize=9,
                        ha="right", color=colors["livre"])
    fig.tight_layout()

    state_t = {"t": 0.0}
    n_per   = 4

    def update(frame):
        state_t["t"] += n_per * dt
        t = state_t["t"]

        for name in ["grafeno", "livre"]:
            s   = states[name]
            P11, P12, P21, P22 = s["P"]
            for _ in range(n_per):
                phi1 = np.fft.fft(s["psi1"])
                phi2 = np.fft.fft(s["psi2"])
                s["psi1"] = np.fft.ifft(P11*phi1 + P12*phi2)
                s["psi2"] = np.fft.ifft(P21*phi1 + P22*phi2)

            dens = (np.abs(s["psi1"])**2 + np.abs(s["psi2"])**2).real
            lines[name].set_ydata(dens)

            # Largura do pacote Δx
            norm  = np.sum(dens) * dx
            xmean = np.sum(xs * dens) * dx / (norm + 1e-15)
            x2    = np.sum(xs**2 * dens) * dx / (norm + 1e-15)
            width = np.sqrt(max(x2 - xmean**2, 0))

            if name == "grafeno":
                txt_g.set_text(f"t = {t:.1f}")
                txt_w_g.set_text(f"Δx = {width:.2f}")
            else:
                txt_l.set_text(f"t = {t:.1f}")
                txt_w_l.set_text(f"Δx = {width:.2f}")

    ani = animation.FuncAnimation(fig, update, frames=300, interval=28, blit=False)

    if args.get("save"):
        fname = "dirac_graphene.gif"
        print(f"Salvando '{fname}'...")
        ani.save(fname, writer=animation.PillowWriter(fps=28))
        print(f"Salvo: {fname}")
    else:
        plt.show()

    return ani


# ============================================================================
# FIGURA ESTÁTICA COMPLETA
# ============================================================================

def plot_all_static(args):
    fig = plt.figure(figsize=(18, 11))
    fig.suptitle(
        "Grafeno — Férmions de Dirac Sem Massa  "
        f"(v_F = c/{int(round(1/VF_FRAC))} ≈ {VF_SI:.0e} m/s)",
        fontsize=14
    )
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.50, wspace=0.40)

    ax_lat  = fig.add_subplot(gs[0, 0])
    ax_bz   = fig.add_subplot(gs[0, 1])
    ax_cone = fig.add_subplot(gs[0, 2], projection='3d')
    ax_band = fig.add_subplot(gs[0, 3])
    ax_prop = fig.add_subplot(gs[1, 0:2])
    ax_ps   = fig.add_subplot(gs[1, 2])
    ax_kl   = fig.add_subplot(gs[1, 3])

    plot_lattice_and_bz(ax_lat, ax_bz)
    plot_dirac_cone(ax_cone, ax_band)
    plot_propagation(ax_prop)
    plot_pseudospin(ax_ps, ax_kl)

    plt.savefig("dirac_graphene.png", dpi=150, bbox_inches="tight")
    print("Figura salva: dirac_graphene.png")
    plt.show()


def plot_landau_static(args):
    fig, (ax_ll, ax_comp) = plt.subplots(1, 2, figsize=(13, 6))
    fig.suptitle("Níveis de Landau em Grafeno vs Elétron Livre", fontsize=13)
    plot_landau_graphene(ax_ll, ax_comp)
    plt.tight_layout()
    plt.savefig("dirac_graphene_landau.png", dpi=150, bbox_inches="tight")
    print("Figura salva: dirac_graphene_landau.png")
    plt.show()


# ============================================================================
# Sumário físico
# ============================================================================

def print_summary(args):
    B = args["B"]
    print(f"\n{'='*64}")
    print("  Grafeno — Férmions de Dirac Sem Massa")
    print(f"{'='*64}")
    print(f"  v_F (u.n.)       = {VF}")
    print(f"  v_F (SI)         = {VF_SI:.2e} m/s  ≈ c/{int(round(C_SI/VF_SI))}")
    print(f"  Massa efetiva    = 0  (fermion de Dirac sem massa)")
    print(f"  Gap              = 0  (semimetal)")
    print()
    print("  Níveis de Landau — grafeno (E_n = v_F√(2nB)):")
    for n in range(6):
        E = landau_graphene(n, B)
        print(f"    n={n}: E = ±{E:.4f}  (elétron e buraco)")
    print()
    print("  Klein em grafeno: θ=0 → T = 1  SEMPRE (independente de V₀)")
    print("  (Katsnelson et al., Nature Physics 2006)")
    print()
    print("  Pseudospin: helicidade travada à direção de k")
    print("   → retro-espalhamento em 180° é proibido por simetria!")
    print(f"{'='*64}\n")


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Grafeno 2D — Férmions de Dirac Sem Massa",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--module", type=int,   default=0,   help="0=tudo 1=rede 2=cone 3=propagação 4=Landau 5=animação")
    p.add_argument("--B",      type=float, default=1.0, help="Campo magnético (Landau)")
    p.add_argument("--save",   action="store_true",     help="Salvar animação como GIF")
    return vars(p.parse_args())


if __name__ == "__main__":
    args = parse_args()
    print_summary(args)

    mod = args["module"]
    if mod == 0 or mod in [1, 2, 3]:
        plot_all_static(args)
    if mod == 0 or mod == 4:
        plot_landau_static(args)
    if mod == 0 or mod == 5:
        animate_graphene(args)
