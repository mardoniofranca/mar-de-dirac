"""
Simulação dos Níveis de Landau Relativísticos — Equação de Dirac
=================================================================
Elétron em campo magnético B uniforme perpendicular ao plano xy.

Física implementada:
  - Níveis de Landau de Dirac:       E_n = mc² √(1 + 2n·(ħeB/m²c³))
  - Níveis de Landau de Schrödinger: E_n = ħωc(n + 1/2)   [ωc = eB/mc]
  - Comparação direta: Dirac ~ √n  vs  Schrödinger ~ n (espaçamento irregular vs regular)
  - Nível zero de Landau (n=0): único, não-degenerado, exclusivo de Dirac
  - Funções de onda: oscilador harmônico + componente relativística
  - Efeito Hall quântico: degenerescência por unidade de área vs B

Sistema de unidades naturais: ħ = c = m = e = 1  (comprimento magnético l_B = 1/√B)

Uso:
    python dirac_landau.py
    python dirac_landau.py --B 2.0   --nmax 12
    python dirac_landau.py --help
"""

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
import argparse
import platform

if platform.system() == "Darwin":
    matplotlib.use("MacOSX")

# ---------------------------------------------------------------------------
# Constantes e energias
# ---------------------------------------------------------------------------

def landau_dirac(n, B, m=1.0):
    """
    Nível n de Landau relativístico (Dirac) em unidades naturais.

    E_n = m √(1 + 2n B/m²)    para n = 0, 1, 2, ...

    Derivação: gauge de Landau A = (0, Bx, 0).
    O hamiltoniano de Dirac com campo magnético leva a um problema
    de oscilador harmônico para as componentes do espinor.
    Autovalor exato: E_n² = m² + 2n·eB  (unidades: ħ=c=e=1)
    """
    return m * np.sqrt(1.0 + 2.0 * n * B / m**2)


def landau_schrodinger(n, B, m=1.0):
    """
    Nível n de Landau não-relativístico (Schrödinger).
    E_n = ωc (n + 1/2),   ωc = eB/m = B (unidades naturais)
    """
    omega_c = B / m
    return omega_c * (n + 0.5)


def landau_schrodinger_relativistic_correction(n, B, m=1.0):
    """
    Correção relativística de 1ª ordem aos níveis de Landau.
    ΔE ≈ −(ħωc)² (2n+1) / (2mc²)
    """
    omega_c = B / m
    E0 = landau_schrodinger(n, B, m)
    delta = -(omega_c**2) * (2*n + 1) / (2 * m)
    return E0 + delta


# ---------------------------------------------------------------------------
# Funções de onda (gauge de Landau)
# ---------------------------------------------------------------------------

def hermite_physicist(n, x):
    """Polinômio de Hermite H_n(x) via recorrência."""
    if n == 0:
        return np.ones_like(x)
    elif n == 1:
        return 2.0 * x
    H_prev, H_curr = np.ones_like(x), 2.0 * x
    for k in range(2, n + 1):
        H_next = 2.0 * x * H_curr - 2.0 * (k - 1) * H_prev
        H_prev, H_curr = H_curr, H_next
    return H_curr


def wavefunction_landau(n, x, B, m=1.0):
    """
    Função de onda do oscilador harmônico magnético (componente grande G).
    ψ_n(x) = C_n · H_n(x/l_B) · exp(-x²/2l_B²)
    l_B = 1/√B  (comprimento magnético)

    Para Dirac, o espinor é:
        Ψ_n = (G_{n},  F_{n-1})ᵀ
    onde G_n é a função do nível n e F_{n-1} do nível n-1 (componente pequena).
    """
    l_B = 1.0 / np.sqrt(B)
    xi  = x / l_B

    # Normalização: ∫|ψ|² dx = 1
    log_norm = -0.5 * np.log(2**n * math.factorial(n)) - 0.25 * np.log(np.pi) - 0.5 * np.log(l_B)
    norm = np.exp(log_norm)

    H   = hermite_physicist(n, xi)
    psi = norm * H * np.exp(-0.5 * xi**2)
    return psi


def wavefunction_dirac_spinor(n, x, B, m=1.0):
    """
    Espinor de Dirac completo para o n-ésimo nível de Landau.
    G_n: componente grande (nível n do oscilador)
    F_n: componente pequena (nível n-1 do oscilador, com coeficiente relativístico)

    Para n=0: F = 0  (nível zero de Landau é puramente spin-up)
    """
    G = wavefunction_landau(n, x, B, m)
    if n == 0:
        F = np.zeros_like(x)
    else:
        E_n = landau_dirac(n, B, m)
        # coeficiente da componente pequena: c_n = √(2nB)/(E_n + m)
        c_n = np.sqrt(2.0 * n * B) / (E_n + m)
        F   = c_n * wavefunction_landau(n - 1, x, B, m)
    return G, F


# ---------------------------------------------------------------------------
# Figura principal
# ---------------------------------------------------------------------------

def plot_all(args):
    B    = args["B"]
    nmax = args["nmax"]
    m    = args["mass"]

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(
        f"Níveis de Landau Relativísticos — Equação de Dirac  "
        f"(B = {B}, m = {m}, unidades naturais)",
        fontsize=13
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.38)

    ax_levels = fig.add_subplot(gs[0, 0])     # diagrama de níveis
    ax_disp   = fig.add_subplot(gs[0, 1])     # E vs n  (dispersão)
    ax_split  = fig.add_subplot(gs[0, 2])     # espaçamento ΔE_n
    ax_wf0    = fig.add_subplot(gs[1, 0])     # função de onda n=0,1,2
    ax_wfsp   = fig.add_subplot(gs[1, 1])     # espinor: G e F
    ax_deg    = fig.add_subplot(gs[1, 2])     # degenerescência vs B

    ns = np.arange(0, nmax + 1)

    E_dirac = np.array([landau_dirac(n, B, m) for n in ns]) - m   # energia cinética
    E_schr  = np.array([landau_schrodinger(n, B, m) for n in ns])
    E_corr  = np.array([landau_schrodinger_relativistic_correction(n, B, m) for n in ns])

    cmap = plt.cm.viridis
    colors_n = [cmap(n / (nmax + 1)) for n in ns]

    # ── Painel 1: Diagrama de níveis ──────────────────────────────────────
    x_d, x_s = 0.2, 0.8
    w = 0.25
    for n in ns:
        c = colors_n[n]
        # Dirac
        ax_levels.hlines(E_dirac[n], x_d - w, x_d + w, color=c, lw=2.0)
        ax_levels.text(x_d + w + 0.03, E_dirac[n], f"n={n}", fontsize=7.5,
                       va="center", color=c)
        # Schrödinger
        ax_levels.hlines(E_schr[n],  x_s - w, x_s + w, color=c, lw=2.0, ls="--")

    ax_levels.set_xlim(0, 1.2)
    ax_levels.set_ylim(-0.2, max(E_dirac[-1], E_schr[-1]) * 1.1)
    ax_levels.set_xticks([x_d, x_s])
    ax_levels.set_xticklabels(["Dirac\n(√n)", "Schrödinger\n(n)"], fontsize=9)
    ax_levels.set_ylabel("Energia cinética (unidades naturais)", fontsize=9)
    ax_levels.set_title("Diagrama de níveis", fontsize=10)
    ax_levels.grid(axis='y', lw=0.4, alpha=0.35)

    # Seta apontando para n=0 (exclusivo de Dirac)
    ax_levels.annotate("n=0 único!\n(só Dirac)",
                       xy=(x_d, E_dirac[0]), xytext=(x_d - 0.18, E_dirac[0] + 0.3),
                       fontsize=7.5, color="#D85A30",
                       arrowprops=dict(arrowstyle="->", color="#D85A30", lw=1.0))

    # ── Painel 2: E vs n (dispersão) ──────────────────────────────────────
    n_dense = np.linspace(0, nmax, 300)
    E_d_cont = np.array([landau_dirac(n, B, m) - m for n in n_dense])
    E_s_cont = np.array([landau_schrodinger(n, B, m) for n in n_dense])

    ax_disp.plot(n_dense, E_d_cont, color="#378ADD", lw=2.0, label="Dirac ∝ √n")
    ax_disp.plot(n_dense, E_s_cont, color="#D85A30", lw=2.0, ls="--", label="Schrödinger ∝ n")
    ax_disp.scatter(ns, E_dirac, color="#378ADD", s=50, zorder=5)
    ax_disp.scatter(ns, E_schr,  color="#D85A30", s=50, zorder=5, marker="s")

    # referência √n
    sqrt_ref = np.sqrt(n_dense + 1e-9) * np.sqrt(2 * B)
    ax_disp.plot(n_dense[1:], sqrt_ref[1:], color="#378ADD", lw=0.8, ls=":",
                 alpha=0.5, label=r"$\sqrt{2nB}$ (ref.)")

    ax_disp.set_xlabel("Nível n", fontsize=10)
    ax_disp.set_ylabel("Energia cinética", fontsize=10)
    ax_disp.set_title("Dispersão E(n): raiz vs linear", fontsize=10)
    ax_disp.legend(fontsize=8.5)
    ax_disp.grid(lw=0.4, alpha=0.4)

    # ── Painel 3: Espaçamento ΔE_n = E_n − E_{n-1} ───────────────────────
    ns_gap = np.arange(1, nmax + 1)
    dE_d = np.diff(E_dirac)
    dE_s = np.diff(E_schr)

    ax_split.plot(ns_gap, dE_d, color="#378ADD", lw=2, marker="o", ms=5,
                  label="Dirac (decresce com n)")
    ax_split.plot(ns_gap, dE_s, color="#D85A30", lw=2, ls="--", marker="s", ms=5,
                  label=f"Schrödinger (const = ωc={B/m:.2f})")
    ax_split.axhline(B / m, color="#D85A30", lw=0.8, ls=":", alpha=0.5)

    ax_split.set_xlabel("n", fontsize=10)
    ax_split.set_ylabel("ΔE_n = E_n − E_{n−1}", fontsize=10)
    ax_split.set_title("Espaçamento entre níveis\n(Dirac: comprime; Schrödinger: uniforme)", fontsize=10)
    ax_split.legend(fontsize=8.5)
    ax_split.grid(lw=0.4, alpha=0.4)

    # ── Painel 4: Funções de onda (densidade) ─────────────────────────────
    l_B = 1.0 / np.sqrt(B)
    x_wf = np.linspace(-6 * l_B, 6 * l_B, 800)

    for n in [0, 1, 2, 3]:
        G, F = wavefunction_dirac_spinor(n, x_wf, B, m)
        dens = G**2 + F**2
        offset = n * 0.4   # empilha verticalmente
        ax_wf0.plot(x_wf / l_B, dens + offset, color=colors_n[n], lw=1.8,
                    label=f"n={n}")
        ax_wf0.axhline(offset, color="gray", lw=0.4, ls="--", alpha=0.4)

    ax_wf0.set_xlabel("x / l_B", fontsize=10)
    ax_wf0.set_ylabel(r"$|G|^2 + |F|^2$ + offset", fontsize=9)
    ax_wf0.set_title("Densidades radiais  (deslocadas)", fontsize=10)
    ax_wf0.legend(fontsize=8.5, loc="upper right")
    ax_wf0.grid(lw=0.4, alpha=0.4)

    # ── Painel 5: Espinor G e F para n=2 ─────────────────────────────────
    n_sp = 3
    G, F = wavefunction_dirac_spinor(n_sp, x_wf, B, m)
    ratio = np.max(np.abs(F)) / (np.max(np.abs(G)) + 1e-15)

    ax_wfsp.plot(x_wf / l_B, G,    color="#378ADD", lw=2,   label=r"$G_n$ — grande")
    ax_wfsp.plot(x_wf / l_B, F,    color="#D85A30", lw=2,   label=r"$F_{n-1}$ — pequena")
    ax_wfsp.plot(x_wf / l_B, G**2 + F**2, color="#1D9E75", lw=1.5, ls="--",
                 alpha=0.8, label=r"$|G|^2+|F|^2$")
    ax_wfsp.axhline(0, color="gray", lw=0.5, ls="--")
    ax_wfsp.set_xlabel("x / l_B", fontsize=10)
    ax_wfsp.set_ylabel("Amplitude", fontsize=10)
    ax_wfsp.set_title(f"Espinor de Dirac: n={n_sp}\n"
                      f"max|F|/max|G| = {ratio:.4f}  (≈√(2nB)/(E+m))", fontsize=9)
    ax_wfsp.legend(fontsize=8.5)
    ax_wfsp.grid(lw=0.4, alpha=0.4)
    ax_wfsp.text(0.97, 0.05,
                 f"E_{n_sp} = {landau_dirac(n_sp,B,m):.4f} mc²\n"
                 f"(cin. = {landau_dirac(n_sp,B,m)-m:.4f} mc²)",
                 transform=ax_wfsp.transAxes, fontsize=8,
                 ha="right", va="bottom", color="#378ADD")

    # ── Painel 6: Degenerescência vs B ────────────────────────────────────
    # Degenerescência por unidade de área: D = eB/(2πħ) = B/(2π) [unidades naturais]
    # É a densidade de estados por nível de Landau
    B_arr = np.linspace(0.1, 5.0, 300)
    D_arr = B_arr / (2 * np.pi)   # por unidade de área, por nível

    ax_deg.plot(B_arr, D_arr, color="#378ADD", lw=2, label=r"D = B/2π  (por nível)")
    ax_deg.fill_between(B_arr, 0, D_arr, alpha=0.12, color="#378ADD")
    ax_deg.axvline(B, color="gray", lw=0.8, ls="--", label=f"B atual = {B}")
    ax_deg.set_xlabel("Campo magnético B", fontsize=10)
    ax_deg.set_ylabel("D = B/(2π)  [por área]", fontsize=10)
    ax_deg.set_title("Degenerescência por nível de Landau\n(Efeito Hall Quântico)", fontsize=10)
    ax_deg.legend(fontsize=8.5)
    ax_deg.grid(lw=0.4, alpha=0.4)

    # Anotação: preenchimento fracionário
    for nu, lbl in [(1, "ν=1"), (2, "ν=2"), (4, "ν=4")]:
        B_fill = nu * 0.3   # exemplo: densidade eletrônica n_e = 0.3
        if B_fill <= B_arr[-1]:
            ax_deg.axvline(B_fill, color="#EF9F27", lw=0.8, ls=":", alpha=0.7)
            ax_deg.text(B_fill, D_arr[-1]*0.85, lbl, fontsize=7.5,
                        ha="center", color="#EF9F27")

    plt.savefig("dirac_landau.png", dpi=150, bbox_inches="tight")
    print("Figura salva: dirac_landau.png")
    plt.show()


# ---------------------------------------------------------------------------
# Tabela de energias
# ---------------------------------------------------------------------------

def print_table(args):
    B    = args["B"]
    nmax = args["nmax"]
    m    = args["mass"]
    omega_c = B / m

    print(f"\n{'='*72}")
    print("  Níveis de Landau Relativísticos — Equação de Dirac")
    print(f"{'='*72}")
    print(f"  B = {B},  m = {m},  ωc = eB/m = {omega_c:.4f},  l_B = 1/√B = {1/np.sqrt(B):.4f}")
    print(f"\n  {'n':>3}  {'E_Dirac':>12}  {'E_cin(D)':>10}  {'E_Schrö':>10}  "
          f"{'ΔE_D':>10}  {'ΔE_S':>10}  {'razão':>8}")
    print("  " + "─" * 68)

    E_prev_d = None
    for n in range(nmax + 1):
        Ed  = landau_dirac(n, B, m)
        Ec  = Ed - m
        Es  = landau_schrodinger(n, B, m)
        dEd = Ec - (landau_dirac(n-1,B,m)-m) if n > 0 else float('nan')
        dEs = Es - landau_schrodinger(n-1,B,m) if n > 0 else float('nan')
        rat = dEd / dEs if (n > 0 and dEs > 0) else float('nan')

        dEd_s = f"{dEd:.5f}" if n > 0 else "     —    "
        dEs_s = f"{dEs:.5f}" if n > 0 else "     —    "
        rat_s = f"{rat:.4f}"  if n > 0 else "    —    "

        mark = " ← n=0 exclusivo!" if n == 0 else ""
        print(f"  {n:>3}  {Ed:>12.6f}  {Ec:>10.6f}  {Es:>10.6f}  "
              f"{dEd_s:>10}  {dEs_s:>10}  {rat_s:>8}{mark}")

    print(f"\n  ΔE_D decresce com n (∝ 1/√n),  ΔE_S = constante = ωc = {omega_c:.4f}")
    print(f"  Padrão Dirac: √1, √2, √3, ...   vs   Schrödinger: 1, 2, 3, ...")
    print(f"\n  Grau de relatividade: max|F|/max|G| para n=1 = "
          f"{np.sqrt(2*B)/(landau_dirac(1,B,m)+m):.5f}")
    print(f"{'='*72}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Níveis de Landau Relativísticos — Equação de Dirac",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--B",    type=float, default=1.0,  help="Campo magnético B")
    p.add_argument("--nmax", type=int,   default=10,   help="Nível máximo")
    p.add_argument("--mass", type=float, default=1.0,  help="Massa (u.n.)")
    return vars(p.parse_args())


if __name__ == "__main__":
    args = parse_args()
    print_table(args)
    plot_all(args)
