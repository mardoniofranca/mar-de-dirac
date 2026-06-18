"""
Simulação do Mar de Dirac e Antimatéria
========================================
Módulo 1 — Espectro de energia e o mar de Dirac
Módulo 2 — Propagador de Feynman-Stückelberg (elétron para t>0, pósitron para t<0)
Módulo 3 — Criação de par elétron-pósitron por campo elétrico externo (efeito Schwinger)
Módulo 4 — Aniquilação: pacote de onda elétron + pósitron → colapso

Sistema de unidades naturais: ħ = c = m = 1

Física implementada:
  - Espectro E = ±√(k²+m²): dois ramos, separados por gap de 2mc²
  - Mar de Dirac: estados de energia negativa preenchidos (vácuo)
  - Buraco no mar = pósitron (carga +e, mesma massa)
  - Propagador de Feynman: partícula → frente, antipartícula ← trás
  - Taxa de criação de pares de Schwinger: Γ ∝ exp(-πm²c³/eEħ)
  - Evolução de pacote misto elétron+pósitron com split-operator FFT

Uso:
    python dirac_antimatter.py            # todos os módulos
    python dirac_antimatter.py --module 1 # só espectro
    python dirac_antimatter.py --module 3 # só Schwinger
    python dirac_antimatter.py --save     # salva GIF da animação
"""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch
import argparse
import platform
import warnings
warnings.filterwarnings("ignore")

if platform.system() == "Darwin":
    matplotlib.use("MacOSX")

# ============================================================================
# MÓDULO 1 — Espectro e Mar de Dirac
# ============================================================================

def plot_dirac_sea(ax, m=1.0, k_max=4.0):
    """
    Diagrama do espectro de Dirac:
      - Ramo positivo:  E = +√(k²+m²)   (elétrons)
      - Gap proibido:   -mc² < E < +mc²  (2mc² de largura)
      - Ramo negativo:  E = -√(k²+m²)   (estados preenchidos = mar de Dirac)
    """
    k = np.linspace(-k_max, k_max, 600)
    E_pos =  np.sqrt(k**2 + m**2)
    E_neg = -np.sqrt(k**2 + m**2)

    # Mar de Dirac (preenchido) — sombra abaixo de -m
    ax.fill_between(k, E_neg, -k_max * 1.5, alpha=0.18, color="#378ADD",
                    label="Mar de Dirac\n(estados E<−mc² preenchidos)")

    # Curvas de dispersão
    ax.plot(k, E_pos, color="#378ADD", lw=2.5, label=r"$E = +\sqrt{k^2+m^2}$  (elétron)")
    ax.plot(k, E_neg, color="#D85A30", lw=2.5, label=r"$E = -\sqrt{k^2+m^2}$  (mar)")

    # Gap
    ax.fill_between([-k_max, k_max], [-m, -m], [m, m],
                    alpha=0.08, color="#EF9F27", zorder=2)
    ax.axhline( m, color="#EF9F27", lw=1.2, ls="--")
    ax.axhline(-m, color="#EF9F27", lw=1.2, ls="--")
    ax.text(k_max * 0.55, 0, f"Gap = 2mc² = {2*m:.1f}",
            va="center", ha="left", fontsize=8.5, color="#EF9F27",
            bbox=dict(fc="white", ec="#EF9F27", alpha=0.8, pad=2))

    # Elétron no ramo positivo
    k_e = 1.5
    E_e = np.sqrt(k_e**2 + m**2)
    ax.scatter([k_e], [E_e], s=120, color="#378ADD", zorder=6)
    ax.annotate("elétron\n(E>0, k>0)", xy=(k_e, E_e),
                xytext=(k_e + 0.8, E_e + 0.5), fontsize=8,
                arrowprops=dict(arrowstyle="->", color="#378ADD"),
                color="#378ADD")

    # Buraco no mar = pósitron
    k_p = -1.5
    E_p = -np.sqrt(k_p**2 + m**2)
    ax.scatter([k_p], [E_p], s=120, color="white",
               edgecolors="#D85A30", linewidths=2, zorder=6)
    ax.annotate("buraco\n(pósitron)\nk_pós = −k_buraco", xy=(k_p, E_p),
                xytext=(k_p - 2.2, E_p + 0.6), fontsize=8,
                arrowprops=dict(arrowstyle="->", color="#D85A30"),
                color="#D85A30")

    ax.set_xlabel("Momento k (mc)", fontsize=10)
    ax.set_ylabel("Energia E (mc²)", fontsize=10)
    ax.set_title("Espectro de Dirac e Mar de Dirac", fontsize=10)
    ax.set_xlim(-k_max, k_max)
    ax.set_ylim(-3.5, 3.5)
    ax.axhline(0, color="gray", lw=0.5, ls=":")
    ax.axvline(0, color="gray", lw=0.5, ls=":")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(lw=0.4, alpha=0.35)


# ============================================================================
# MÓDULO 2 — Propagador de Feynman-Stückelberg
# ============================================================================

def plot_feynman_propagator(ax, m=1.0, x_max=8.0, t_max=6.0):
    """
    Propagador de Feynman no espaço (x,t):
      - Para t > 0: propaga estados de energia positiva  → elétron para frente
      - Para t < 0: propaga estados de energia negativa  → pósitron para trás
    Visualiza o cone de luz e os dois caminhos.
    """
    # Cone de luz
    t_arr = np.linspace(0, t_max, 200)
    ax.fill_between( t_arr,  t_arr, x_max,  alpha=0.06, color="#378ADD")
    ax.fill_between( t_arr, -t_arr, -x_max, alpha=0.06, color="#378ADD")
    ax.plot( t_arr,  t_arr, color="gray", lw=1, ls=":", alpha=0.6)
    ax.plot( t_arr, -t_arr, color="gray", lw=1, ls=":", alpha=0.6)
    t_neg = np.linspace(-t_max, 0, 200)
    ax.fill_between(t_neg,  -t_neg, x_max,  alpha=0.06, color="#D85A30")
    ax.fill_between(t_neg,   t_neg, -x_max, alpha=0.06, color="#D85A30")
    ax.plot(t_neg, -t_neg, color="gray", lw=1, ls=":", alpha=0.6)
    ax.plot(t_neg,  t_neg, color="gray", lw=1, ls=":", alpha=0.6)

    # Trajetória do elétron: t>0, x>0 (propagação causal)
    t_e = np.linspace(0, 4.5, 100)
    x_e = 0.6 * t_e
    ax.annotate("", xy=(t_e[-1], x_e[-1]), xytext=(t_e[0], x_e[0]),
                arrowprops=dict(arrowstyle="-|>", color="#378ADD",
                                lw=2, mutation_scale=15))
    ax.plot(t_e, x_e, color="#378ADD", lw=2)
    ax.text(t_e[-1] + 0.15, x_e[-1], "e⁻\n(t>0)", fontsize=8.5, color="#378ADD")

    # Trajetória do pósitron: t<0 interpretado como propagação para trás no tempo
    t_p = np.linspace(-4.5, 0, 100)
    x_p = -0.5 * t_p
    ax.annotate("", xy=(t_p[0], x_p[0]), xytext=(t_p[-1], x_p[-1]),
                arrowprops=dict(arrowstyle="-|>", color="#D85A30",
                                lw=2, mutation_scale=15))
    ax.plot(t_p, x_p, color="#D85A30", lw=2)
    ax.text(t_p[0] - 0.2, x_p[0] + 0.3, "e⁺\n(t<0 = anti)", fontsize=8.5,
            color="#D85A30", ha="right")

    # Vértice de criação de par
    ax.scatter([0], [0], s=200, color="#1D9E75", zorder=8, marker="*")
    ax.text(0.15, 0.3, "criação\nde par", fontsize=8, color="#1D9E75")

    ax.axvline(0, color="gray", lw=1.0, ls="--", alpha=0.8)
    ax.axhline(0, color="gray", lw=0.5, ls=":", alpha=0.5)
    ax.set_xlabel("Tempo t", fontsize=10)
    ax.set_ylabel("Posição x", fontsize=10)
    ax.set_title("Propagador de Feynman-Stückelberg\n"
                 "e⁻ → frente no tempo  |  e⁺ → atrás no tempo", fontsize=9)
    ax.set_xlim(-t_max, t_max)
    ax.set_ylim(-x_max, x_max)
    ax.text(-t_max*0.95, x_max*0.85, "t < 0\n(pósitron)", fontsize=8.5,
            color="#D85A30", alpha=0.8)
    ax.text( t_max*0.25, x_max*0.85, "t > 0\n(elétron)",  fontsize=8.5,
            color="#378ADD", alpha=0.8)
    ax.grid(lw=0.4, alpha=0.35)


# ============================================================================
# MÓDULO 3 — Efeito Schwinger: taxa de criação de pares
# ============================================================================

def schwinger_rate(E_field, m=1.0, alpha=1.0/137.0):
    """
    Taxa de criação de pares elétron-pósitron por campo elétrico E:
        Γ ∝ (eE)² exp(-π m²c³ / eEħ)
    Em unidades naturais (ħ=c=e=1):
        Γ = (E_field/E_crit)² · exp(-π/x),   x = E_field/E_crit
    Campo crítico de Schwinger: E_crit = m²c³/(eħ) ≈ 1.3×10¹⁸ V/m
    """
    E_crit = m**2  # em unidades naturais com e=1
    x = E_field / E_crit
    # Evita underflow
    exponent = -np.pi / (x + 1e-12)
    log_rate = 2 * np.log(x + 1e-12) + exponent
    return np.exp(np.clip(log_rate, -300, 0))


def plot_schwinger(ax1, ax2, m=1.0):
    """
    Taxa de Schwinger Γ(E) e diagrama do mecanismo físico.
    """
    # Painel esquerdo: Γ vs E/E_crit
    E_crit = m**2
    x_arr  = np.linspace(0.05, 3.0, 600)
    E_arr  = x_arr * E_crit
    Gamma  = schwinger_rate(E_arr, m)

    ax1.semilogy(x_arr, Gamma, color="#378ADD", lw=2.5)
    ax1.axvline(1.0, color="#D85A30", lw=1.2, ls="--",
                label=r"$E = E_{crit}$ (Schwinger)")
    ax1.axvline(0.3, color="#1D9E75", lw=1.0, ls=":",
                label="Regime de laboratório\n(lasers ultraintensos)")
    ax1.fill_betweenx([1e-50, 1], 0, 0.3, alpha=0.07, color="#1D9E75")
    ax1.fill_betweenx([1e-50, 1], 1.0, 3.0, alpha=0.07, color="#D85A30")

    ax1.set_xlabel(r"$E_{campo}$ / $E_{crit}$", fontsize=10)
    ax1.set_ylabel(r"Taxa $\Gamma$ (u.a.)", fontsize=10)
    ax1.set_title("Efeito Schwinger:\nTaxa de criação de pares", fontsize=10)
    ax1.set_ylim(1e-50, 10)
    ax1.set_xlim(0, 3)
    ax1.legend(fontsize=8)
    ax1.grid(lw=0.4, alpha=0.4)

    # Painel direito: diagrama de energia do mecanismo
    x_diag = np.linspace(-5, 5, 400)

    # Energia potencial linear: V(x) = -eEx
    E_field_diag = 0.4 * m**2
    V = -E_field_diag * x_diag

    # Ramo do elétron (partícula) deslocado pelo potencial
    E_elec = m + V  # nível de vácuo superior
    E_posi = -m + V  # nível do mar deslocado

    ax2.plot(x_diag, E_elec, color="#378ADD", lw=2, label="Nível e⁻ (E>0)")
    ax2.plot(x_diag, E_posi, color="#D85A30", lw=2, label="Mar de Dirac\n(E<0, deslocado)")
    ax2.fill_between(x_diag, E_posi, -6, alpha=0.15, color="#D85A30")

    # Gap fechando em x=0 no sistema inclinado
    ax2.axhline(0, color="gray", lw=0.8, ls="--", alpha=0.6)

    # Região de tunelamento
    # Onde o elétron pode "tunelar" do mar para o ramo positivo
    x_tunel_l = -m / E_field_diag
    x_tunel_r =  m / E_field_diag
    ax2.fill_between([x_tunel_l, x_tunel_r], [-m, -m], [m, m],
                     alpha=0.15, color="#EF9F27", label="Barreira de tunelamento")
    ax2.annotate("", xy=(x_tunel_r, 0.5), xytext=(x_tunel_l, 0.5),
                 arrowprops=dict(arrowstyle="<->", color="#EF9F27", lw=1.5))
    ax2.text(0, 0.8, "tunelamento\n(produção de par)", ha="center",
             fontsize=8, color="#EF9F27")

    # Seta: pósitron sai para esquerda, elétron para direita
    ax2.annotate("e⁺ →", xy=(x_tunel_l - 1.5, -1.4), fontsize=8.5,
                 color="#D85A30", ha="center")
    ax2.annotate("← e⁻", xy=(x_tunel_r + 1.5,  1.4), fontsize=8.5,
                 color="#378ADD", ha="center")

    ax2.set_xlabel("Posição x", fontsize=10)
    ax2.set_ylabel("Energia", fontsize=10)
    ax2.set_title("Mecanismo físico:\nTunelamento pelo gap de 2mc²", fontsize=10)
    ax2.set_ylim(-5, 5)
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(lw=0.4, alpha=0.4)


# ============================================================================
# MÓDULO 4 — Evolução temporal: elétron + pósitron → aniquilação
# ============================================================================

def build_pair_state(xs, dx, k_e, k_p, sigma, m=1.0):
    """
    Constrói estado inicial com dois pacotes: elétron (E>0) + pósitron (E<0).

    Elétron:  k > 0, ramo E = +√(k²+m²), movendo-se para a direita
    Pósitron: k < 0, ramo E = -√(k²+m²), movendo-se para a esquerda

    Em Dirac 1+1D, a componente do espinor de energia negativa é:
        ψ₁ = (E+m)/√(2E(E+m)) · f(x),  ψ₂ = k/√(2E(E+m)) · f(x)
    mas com E < 0 para o pósitron.
    """
    x0_e =  xs.max() * 0.35   # elétron começa à direita
    x0_p = -xs.max() * 0.35   # pósitron começa à esquerda

    # Pacote elétron (E > 0)
    env_e  = np.exp(-0.25 * ((xs - x0_e) / sigma)**2) / (np.pi * sigma**2)**0.25
    psi1_e = env_e * np.exp(1j * k_e * xs)
    w_e    = np.sqrt(k_e**2 + m**2)
    psi2_e = (k_e / (w_e + m)) * psi1_e

    # Pacote pósitron (E < 0, interpretado como buraco no mar)
    # Pósitron com momento k_p > 0 → elétron de energia negativa com k = -k_p
    k_neg  = -k_p
    env_p  = np.exp(-0.25 * ((xs - x0_p) / sigma)**2) / (np.pi * sigma**2)**0.25
    w_neg  = np.sqrt(k_neg**2 + m**2)
    # Spinor de energia negativa
    psi1_p = env_p * np.exp(1j * k_neg * xs)
    psi2_p = (-k_neg / (w_neg + m)) * psi1_p

    # Superpõe os dois pacotes
    psi1 = psi1_e + psi1_p
    psi2 = psi2_e + psi2_p

    # Normaliza
    norm = np.sqrt(np.sum((np.abs(psi1)**2 + np.abs(psi2)**2) * dx))
    return psi1 / norm, psi2 / norm


def project_components(psi1, psi2, ks, m=1.0):
    """
    Projeta o espinor nos ramos de energia positiva e negativa.
    Retorna densidades ρ_+(x) e ρ_-(x).
    """
    phi1 = np.fft.fft(psi1)
    phi2 = np.fft.fft(psi2)
    omega = np.sqrt(ks**2 + m**2)
    alpha = np.where(omega + m > 1e-12, ks / (omega + m), 0.0)

    # Projeção no ramo positivo: coef A = (phi1 + alpha*phi2) / norm
    A_pos = (phi1 + alpha * phi2) / np.sqrt(2 + 2*alpha**2 + 1e-15)
    A_neg = (phi1 - alpha * phi2) / np.sqrt(2 + 2*alpha**2 + 1e-15)

    rho_pos = np.abs(np.fft.ifft(A_pos))**2
    rho_neg = np.abs(np.fft.ifft(A_neg))**2
    return rho_pos, rho_neg


def make_propagators_full(ks, m, dt):
    """Propagador split-operator para a equação de Dirac livre."""
    omega  = np.sqrt(ks**2 + m**2)
    cos_w  = np.cos(omega * dt)
    sinc_w = np.where(omega > 1e-12, np.sin(omega * dt) / omega, dt)
    P11 =  cos_w - 1j * m  * sinc_w
    P12 = -1j * ks * sinc_w
    P21 = -1j * ks * sinc_w
    P22 =  cos_w + 1j * m  * sinc_w
    return P11, P12, P21, P22


def step_free(psi1, psi2, P11, P12, P21, P22):
    phi1 = np.fft.fft(psi1)
    phi2 = np.fft.fft(psi2)
    return np.fft.ifft(P11*phi1 + P12*phi2), np.fft.ifft(P21*phi1 + P22*phi2)


def animate_pair(args):
    """Animação: elétron e pósitron se aproximam e aniquilam."""
    m     = args["mass"]
    k_e   = args["k0"]
    k_p   = args["k0"]   # pósitron com mesmo momento absoluto
    N     = 2048
    L     = 120.0
    sigma = 7.0
    dt    = 0.05

    xs = np.linspace(-L/2, L/2, N, endpoint=False)
    dx = xs[1] - xs[0]
    ks = np.fft.fftfreq(N, d=dx) * 2 * np.pi

    psi1, psi2 = build_pair_state(xs, dx, k_e, k_p, sigma, m)
    P11, P12, P21, P22 = make_propagators_full(ks, m, dt)

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    fig.suptitle(
        f"Mar de Dirac — Par Elétron-Pósitron  (k={k_e}, m={m})\n"
        "e⁻ (direita→) encontra e⁺ (←esquerda) → aniquilação",
        fontsize=12
    )

    ax_tot, ax_pos, ax_neg = axes

    dens0 = (np.abs(psi1)**2 + np.abs(psi2)**2).real
    rp0, rn0 = project_components(psi1, psi2, ks, m)
    ymax_tot = dens0.max() * 2.5

    line_tot, = ax_tot.plot(xs, dens0, color="#1D9E75", lw=2,
                            label=r"$|\psi|^2$ total")
    ax_tot.set_ylim(0, ymax_tot)
    ax_tot.set_ylabel(r"$|\psi|^2$", fontsize=10)
    ax_tot.set_title("Densidade total", fontsize=10)
    ax_tot.legend(fontsize=9, loc="upper right")
    ax_tot.grid(lw=0.4, alpha=0.4)
    txt_t   = ax_tot.text(0.02, 0.88, "t = 0.00", transform=ax_tot.transAxes, fontsize=9)
    txt_ann = ax_tot.text(0.50, 0.88, "", transform=ax_tot.transAxes, fontsize=9,
                          ha="center", color="#D85A30", fontweight="bold")

    line_pos, = ax_pos.plot(xs, rp0, color="#378ADD", lw=2,
                            label=r"$\rho_+$  —  ramo E>0  (elétron)")
    ax_pos.set_ylim(0, rp0.max() * 2.5 + 1e-10)
    ax_pos.set_ylabel(r"$\rho_+$", fontsize=10)
    ax_pos.set_title("Componente de energia positiva (elétron)", fontsize=10)
    ax_pos.legend(fontsize=9, loc="upper right")
    ax_pos.grid(lw=0.4, alpha=0.4)

    line_neg, = ax_neg.plot(xs, rn0, color="#D85A30", lw=2,
                            label=r"$\rho_-$  —  ramo E<0  (pósitron/mar)")
    ax_neg.set_ylim(0, rn0.max() * 2.5 + 1e-10)
    ax_neg.set_ylabel(r"$\rho_-$", fontsize=10)
    ax_neg.set_title("Componente de energia negativa (pósitron / mar de Dirac)", fontsize=10)
    ax_neg.legend(fontsize=9, loc="upper right")
    ax_neg.set_xlabel("x (unidades naturais)", fontsize=10)
    ax_neg.grid(lw=0.4, alpha=0.4)

    fig.tight_layout()

    state = {"psi1": psi1, "psi2": psi2, "t": 0.0}
    n_per = 4

    def update(frame):
        for _ in range(n_per):
            state["psi1"], state["psi2"] = step_free(
                state["psi1"], state["psi2"], P11, P12, P21, P22)
        state["t"] += n_per * dt

        dens = (np.abs(state["psi1"])**2 + np.abs(state["psi2"])**2).real
        rp, rn = project_components(state["psi1"], state["psi2"], ks, m)

        line_tot.set_ydata(dens)
        line_pos.set_ydata(rp)
        line_neg.set_ydata(rn)

        # Detecta sobreposição (aniquilação)
        dx_center = np.sum(np.abs(xs) * dens * dx) / (np.sum(dens * dx) + 1e-15)
        overlap = np.sum(np.minimum(rp, rn)) * dx
        txt_t.set_text(f"t = {state['t']:.1f}")
        if overlap > 0.01:
            txt_ann.set_text(f"⚡ Aniquilação!  sobreposição={overlap:.3f}")
        else:
            txt_ann.set_text("")

    ani = animation.FuncAnimation(fig, update, frames=300, interval=30, blit=False)

    if args.get("save"):
        fname = "dirac_pair.gif"
        print(f"Salvando '{fname}'...")
        ani.save(fname, writer=animation.PillowWriter(fps=28))
        print(f"Salvo: {fname}")
    else:
        plt.show()

    return ani


# ============================================================================
# FIGURA ESTÁTICA COMPLETA
# ============================================================================

def plot_static_all(args):
    m = args["mass"]

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(
        "Mar de Dirac e Antimatéria — Equação de Dirac",
        fontsize=14, fontweight="normal"
    )
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.50, wspace=0.40)

    ax_sea    = fig.add_subplot(gs[0, 0:2])   # espectro + mar
    ax_feyn   = fig.add_subplot(gs[0, 2])     # propagador de Feynman
    ax_schwi1 = fig.add_subplot(gs[0, 3])     # taxa Schwinger
    ax_schwi2 = fig.add_subplot(gs[1, 0])     # mecanismo Schwinger
    ax_snap   = fig.add_subplot(gs[1, 1:3])   # snapshot do par
    ax_proj   = fig.add_subplot(gs[1, 3])     # projeção ρ+ e ρ-

    plot_dirac_sea(ax_sea, m)
    plot_feynman_propagator(ax_feyn, m)
    plot_schwinger(ax_schwi1, ax_schwi2, m)

    # ── Snapshot: pacote elétron+pósitron no instante inicial ──
    N, L, sigma = 2048, 120.0, 7.0
    xs   = np.linspace(-L/2, L/2, N, endpoint=False)
    dx   = xs[1] - xs[0]
    ks   = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    k0   = args["k0"]
    psi1, psi2 = build_pair_state(xs, dx, k0, k0, sigma, m)
    dens = (np.abs(psi1)**2 + np.abs(psi2)**2).real
    rp, rn = project_components(psi1, psi2, ks, m)

    ax_snap.plot(xs, dens, color="#1D9E75", lw=2,   label=r"$|\psi|^2$ total")
    ax_snap.plot(xs, rp,   color="#378ADD", lw=1.8, ls="--", alpha=0.85,
                 label=r"$\rho_+$  (e⁻)")
    ax_snap.plot(xs, rn,   color="#D85A30", lw=1.8, ls="--", alpha=0.85,
                 label=r"$\rho_-$  (e⁺)")
    ax_snap.set_xlim(-L/2, L/2)
    ax_snap.set_ylim(bottom=0)
    ax_snap.set_xlabel("x", fontsize=10)
    ax_snap.set_ylabel(r"$|\psi|^2$", fontsize=10)
    ax_snap.set_title("Estado inicial: pacote e⁻ (dir.) + e⁺ (esq.)", fontsize=10)
    ax_snap.legend(fontsize=8.5)
    ax_snap.grid(lw=0.4, alpha=0.4)

    # ── Projeção espectral: ρ(k) para e⁻ e e⁺ ──
    phi1 = np.fft.fft(psi1)
    phi2 = np.fft.fft(psi2)
    ks_shift = np.fft.fftshift(ks)
    rho_k = np.fft.fftshift(np.abs(phi1)**2 + np.abs(phi2)**2).real

    ax_proj.plot(ks_shift, rho_k / rho_k.max(), color="#1D9E75", lw=2,
                 label=r"$|\tilde\psi|^2$ (espectro)")
    ax_proj.axvline( k0, color="#378ADD", lw=1.5, ls="--", label=f"k_e = +{k0}")
    ax_proj.axvline(-k0, color="#D85A30", lw=1.5, ls="--", label=f"k_p = −{k0}")
    ax_proj.set_xlim(-4*k0, 4*k0)
    ax_proj.set_xlabel("Momento k", fontsize=10)
    ax_proj.set_ylabel(r"$|\tilde\psi|^2$ (normaliz.)", fontsize=10)
    ax_proj.set_title("Espectro de momento\n(dois picos: e⁻ e e⁺)", fontsize=10)
    ax_proj.legend(fontsize=8.5)
    ax_proj.grid(lw=0.4, alpha=0.4)

    plt.savefig("dirac_antimatter.png", dpi=150, bbox_inches="tight")
    print("Figura salva: dirac_antimatter.png")
    plt.show()


# ============================================================================
# Tabela de referência
# ============================================================================

def print_summary(args):
    m  = args["mass"]
    k0 = args["k0"]
    E  = np.sqrt(k0**2 + m**2)

    print(f"\n{'='*62}")
    print("  Mar de Dirac e Antimatéria — Sumário Físico")
    print(f"{'='*62}")
    print(f"  Massa m           = {m}")
    print(f"  Momento k₀        = {k0}")
    print(f"  Energia elétron   = {E:.4f} mc²")
    print(f"  Gap do espectro   = 2mc² = {2*m:.4f}")
    print(f"  Campo Schwinger   = m² = {m**2:.4f}  (u.n.)")
    print(f"  E_Schwinger (SI)  ≈ 1.32 × 10¹⁸ V/m")
    print()
    print("  Taxa de Schwinger Γ(E/E_crit):")
    for x in [0.1, 0.3, 0.5, 1.0, 2.0]:
        G = schwinger_rate(x * m**2, m)
        print(f"    E/E_crit = {x:.1f}  →  Γ = {G:.3e}")
    print()
    print("  Interpretação de Feynman-Stückelberg:")
    print("    Pósitron = elétron propagando PARA TRÁS no tempo")
    print("    → Mesmo resultado que 'buraco no mar de Dirac'")
    print(f"{'='*62}\n")


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Mar de Dirac e Antimatéria — Equação de Dirac",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--module", type=int, default=0,
                   help="0=tudo, 1=espectro, 2=Feynman, 3=Schwinger, 4=animação")
    p.add_argument("--k0",    type=float, default=1.5, help="Momento inicial")
    p.add_argument("--mass",  type=float, default=1.0, help="Massa (u.n.)")
    p.add_argument("--save",  action="store_true",     help="Salvar animação como GIF")
    return vars(p.parse_args())


if __name__ == "__main__":
    args = parse_args()
    print_summary(args)

    mod = args["module"]
    if mod == 0 or mod != 4:
        plot_static_all(args)
    if mod == 0 or mod == 4:
        animate_pair(args)
