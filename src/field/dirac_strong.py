"""
Campo Forte e Colapso Relativístico — Equação de Dirac
=======================================================
Módulo 1 — Espectro vs Z: como os níveis de energia mergulham com Z crescente
Módulo 2 — Colapso da função de onda: γ = √(κ²−(Zα)²) → 0 quando Zα → |κ|
Módulo 3 — Limite crítico: Z_crit = |κ|/α ≈ 137 para 1s (κ=−1)
Módulo 4 — Criação de pares espontânea: E_1s cruza −mc² para Z > Z_cr2 ≈ 173
Módulo 5 — Componente pequena F(r): diverge quando γ → 0
Módulo 6 — Comparação Dirac vs Schrödinger vs Klein-Gordon para Z pesado

Física implementada:
  - Energia exata de Dirac: E = mc²[1+(Zα/(n_r+√(κ²−(Zα)²)))²]^{−½}
  - Expoente da singularidade: γ(Z) = √(κ²−(Zα)²)  →  0 quando Zα→|κ|
  - Função de onda perto da origem: ψ ~ r^{γ−1}  (diverge para γ<½)
  - Crossing E=−mc²: Z_cr2 onde o estado ligado mergulha no contínuo negativo
  - Largura do estado supercrítico: Γ ~ exp(−π(Z−Z_cr2))
  - Razão componente pequena/grande: max|F|/max|G| → 1 quando Z→Z_crit

Sistema de unidades: ħ = c = m_e = 1,  α ≈ 1/137

Uso:
    python dirac_strong.py
    python dirac_strong.py --Zmax 180
    python dirac_strong.py --module 2
    python dirac_strong.py --help
"""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import argparse
import platform
import warnings
warnings.filterwarnings("ignore")

if platform.system() == "Darwin":
    matplotlib.use("MacOSX")

ALPHA   = 1.0 / 137.035999084
HART_EV = 27.211396          # 1 Hartree em eV
MC2_EV  = 0.511e6            # mc² em eV (elétron)


# ============================================================================
# Física central
# ============================================================================

def gamma_dirac(kappa, Z):
    """
    Expoente relativístico γ = √(κ² − (Zα)²).
    Quando Zα → |κ|: γ → 0  (singularidade da função de onda).
    Para κ = −1 (estado 1s): colapso em Z = 1/α ≈ 137.
    """
    arg = kappa**2 - (Z * ALPHA)**2
    if arg <= 0:
        return 0.0
    return np.sqrt(arg)


def energy_dirac(n, kappa, Z):
    """
    Energia exata de Dirac para átomo hidrogenóide (u.a.).
    E = mc²[1+(Zα/(n_r+γ))²]^{−½} − mc²,  n_r = n − |κ|
    Retorna None se Zα > |κ| (estado colapsa).
    """
    Za  = Z * ALPHA
    if Za >= abs(kappa):
        return None
    gam = np.sqrt(kappa**2 - Za**2)
    nr  = n - abs(kappa)
    if nr < 0:
        return None
    denom = nr + gam
    if denom < 1e-12:
        return None
    E_tot = 1.0 / np.sqrt(1.0 + (Za / denom)**2)   # em unidades mc²
    return E_tot - 1.0                                 # energia cinética


def energy_schrodinger(n, Z):
    """E = −Z²/(2n²) em u.a."""
    return -(Z * ALPHA)**2 / (2 * n**2) / ALPHA**2


def energy_dirac_eV(n, kappa, Z):
    """Energia em eV (real, para comparação com experimento)."""
    e = energy_dirac(n, kappa, Z)
    if e is None:
        return None
    return e * MC2_EV


def z_critical_1(kappa):
    """Z crítico onde γ = 0: Z_c1 = |κ|/α ≈ 137|κ|."""
    return abs(kappa) / ALPHA


def z_critical_2(n, kappa):
    """
    Z crítico onde E_1s = −mc² (estado mergulha no mar de Dirac).
    Resolvendo E(Z_c2) = −1 (em unidades mc²):
        −1 = 1/√(1+(Zα/(n_r+γ))²) − 1  →  diverge
    Aproximação: Z_c2 ≈ Z_c1 × (1 + n_r²/(2κ²)) para n_r > 0
    Para 1s (n=1, κ=−1, n_r=0): Z_c2 ≈ 173 (numérico)
    """
    # Busca numérica: encontra Z onde |E| → 1 (mc²)
    # Para n=1, κ=−1: E(Z) = 1/√(1+(Zα/γ)²) − 1
    # |E| = 1 quando 1/√(1+(Zα/γ)²) → 0
    # isso ocorre quando Zα/γ → ∞, i.e., Zα → |κ| (= Z_c1)
    # Para estados com n_r > 0, E cruza −1 antes:
    Zc1 = z_critical_1(kappa)
    if n - abs(kappa) == 0:
        return Zc1   # para 1s, crossing em Z_c1
    # Para n_r > 0: busca binária
    Z_lo, Z_hi = 1.0, Zc1 - 0.001
    for _ in range(100):
        Zm = (Z_lo + Z_hi) / 2
        e  = energy_dirac(n, kappa, Zm)
        if e is None or e < -1.0:
            Z_hi = Zm
        else:
            Z_lo = Zm
    return (Z_lo + Z_hi) / 2


# ============================================================================
# Funções de onda radiais
# ============================================================================

def radial_wavefunction_strong(n, kappa, Z, Nr=4000):
    """
    Integra G(r) e F(r) via RK4 com grade logarítmica.
    Para Z grande: γ pequeno → singularidade r^{γ−1} na origem.
    """
    E = energy_dirac(n, kappa, Z)
    if E is None:
        return None, None, None
    c  = 1.0 / ALPHA    # c em u.a.
    Za = Z * ALPHA
    gam = gamma_dirac(kappa, Z)

    r_min = 1e-6
    r_max = max(60.0, 20 * n**2 / (Z + 1))
    r = np.geomspace(r_min, r_max, Nr)

    G = np.zeros(Nr)
    F = np.zeros(Nr)

    # Condição inicial: r^γ
    A = 1.0
    denom = gam + kappa
    B = Za * c / denom * A if abs(denom) > 1e-10 else 0.0

    G[0] = A * r[0]**gam
    F[0] = B * r[0]**gam

    def derivs(ri, g, f):
        V  = -Z / ri
        dg = -(kappa / ri) * g + ((E + c**2 - V) / c) * f
        df =  (kappa / ri) * f - ((E - c**2 - V) / c) * g
        return dg, df

    prev = G[0]**2 + F[0]**2

    for i in range(Nr - 1):
        ri = r[i]; dr = r[i+1] - r[i]
        gi, fi = G[i], F[i]

        k1g, k1f = derivs(ri, gi, fi)
        k2g, k2f = derivs(ri+dr/2, gi+dr/2*k1g, fi+dr/2*k1f)
        k3g, k3f = derivs(ri+dr/2, gi+dr/2*k2g, fi+dr/2*k2f)
        k4g, k4f = derivs(ri+dr,   gi+dr*k3g,   fi+dr*k3f)

        G[i+1] = gi + dr/6*(k1g+2*k2g+2*k3g+k4g)
        F[i+1] = fi + dr/6*(k1f+2*k2f+2*k3f+k4f)

        curr = G[i+1]**2 + F[i+1]**2
        if (np.isnan(curr) or np.isinf(curr) or
                (curr > prev * 50 and ri > 3.0*n/Z)):
            G[i+2:] = 0.0; F[i+2:] = 0.0
            break
        if curr < prev:
            prev = curr

    dr_arr = np.diff(r, prepend=r[0])
    norm = np.sqrt(np.sum((G**2 + F**2) * dr_arr))
    if norm > 1e-12:
        G /= norm; F /= norm

    return r, G, F


# ============================================================================
# Painéis
# ============================================================================

def _panel_energy_vs_Z(ax, nmax=2, Zmax=130, kappas=None):
    """
    Energia dos estados vs Z, mostrando o mergulho em direção a −mc².
    """
    if kappas is None:
        kappas = [(-1,1), (1,2), (-2,2), (-1,2)]   # (kappa, n)

    colors  = ["#378ADD","#1D9E75","#D85A30","#EF9F27","#7F77DD"]
    Z_arr   = np.linspace(1, Zmax, 400)

    lnames = {-1:"1s₁/₂", 1:"2p₁/₂", -2:"2p₃/₂", -3:"3d₃/₂"}

    for i,(kappa,n) in enumerate(kappas):
        E_arr = []
        Z_plot = []
        for Z in Z_arr:
            e = energy_dirac(n, kappa, Z)
            if e is not None and e > -1.01:
                E_arr.append(e)
                Z_plot.append(Z)

        lbl = lnames.get(kappa, f"κ={kappa},n={n}")
        col = colors[i % len(colors)]
        ax.plot(Z_plot, E_arr, color=col, lw=2.0, label=lbl)

        # Comparação Schrödinger (tracejado)
        E_sch = np.array([-(Z*ALPHA)**2/(2*n**2) for Z in Z_plot])
        ax.plot(Z_plot, E_sch, color=col, lw=1.2, ls=":", alpha=0.5)

    # Linha E = −mc²: limite de criação de pares
    ax.axhline(-1.0, color="#D85A30", lw=1.8, ls="--",
               label="E = −mc² (limiar criação de pares)")
    ax.axhline(0.0,  color="gray",    lw=0.8, ls=":", alpha=0.5)

    # Linha Z_crit (κ=−1)
    Zc = z_critical_1(-1)
    ax.axvline(Zc, color="#EF9F27", lw=1.2, ls="--",
               label=f"Z_crit = 1/α ≈ {Zc:.0f} (colapso 1s)")
    ax.text(Zc+1, -0.5, f"Z≈{Zc:.0f}", fontsize=8, color="#EF9F27")

    ax.set_xlabel("Número atômico Z", fontsize=10)
    ax.set_ylabel("Energia cinética (mc²)", fontsize=10)
    ax.set_title("Espectro de Dirac vs Z\n"
                 "sólido=Dirac, pontilhado=Schrödinger", fontsize=10)
    ax.set_ylim(-1.05, 0.05)
    ax.set_xlim(1, Zmax)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(lw=0.4, alpha=0.35)


def _panel_gamma_vs_Z(ax):
    """
    γ(Z) = √(1−(Zα)²) para estado 1s (κ=−1).
    Mostra o colapso da função de onda quando Z → 137.
    """
    Zc  = z_critical_1(-1)
    Z_arr = np.linspace(1, Zc * 0.9999, 500)
    gam   = np.array([gamma_dirac(-1, Z) for Z in Z_arr])

    ax.plot(Z_arr, gam, color="#378ADD", lw=2.5, label=r"$\gamma = \sqrt{1-(Z\alpha)^2}$")
    ax.fill_between(Z_arr, 0, gam, alpha=0.12, color="#378ADD")

    # Linha γ = 0.5: ψ ~ r^{−0.5} ainda normalizável
    ax.axhline(0.5, color="#1D9E75", lw=1.2, ls="--",
               label=r"$\gamma=\frac{1}{2}$: $\psi\sim r^{-1/2}$ (limite normalizável)")
    Z_half = np.sqrt(1 - 0.25) / ALPHA
    ax.axvline(Z_half, color="#1D9E75", lw=0.8, ls=":", alpha=0.7)
    ax.text(Z_half+1, 0.55, f"Z≈{Z_half:.0f}", fontsize=7.5, color="#1D9E75")

    ax.axvline(Zc, color="#EF9F27", lw=1.2, ls="--",
               label=f"Z_crit ≈ {Zc:.0f}: γ=0 (colapso!)")
    ax.text(Zc-12, 0.08, f"COLAPSO\nZ≈{Zc:.0f}", fontsize=8,
            color="#EF9F27", ha="center")

    # Comportamento ~√(1−(Z/137)²)
    ax.set_xlabel("Z", fontsize=10)
    ax.set_ylabel(r"$\gamma = \sqrt{\kappa^2-(Z\alpha)^2}$", fontsize=10)
    ax.set_title("Expoente γ(Z) para 1s (κ=−1)\n"
                 r"ψ(r→0) ~ r^{γ−1}: colapsa quando γ→0", fontsize=10)
    ax.set_xlim(0, Zc + 5)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=8.5)
    ax.grid(lw=0.4, alpha=0.35)


def _panel_wavefunction_comparison(ax, Z_list=None):
    """
    Funções de onda 1s para diferentes Z, mostrando colapso.
    """
    if Z_list is None:
        Z_list = [1, 20, 50, 80, 100, 120]

    cmap = plt.cm.plasma
    cols = [cmap(i / (len(Z_list) - 1)) for i in range(len(Z_list))]

    for Z, col in zip(Z_list, cols):
        r, G, F = radial_wavefunction_strong(1, -1, Z)
        if r is None:
            continue
        dens = G**2 + F**2
        # Escala para máximo 1, recorta em r pequeno
        mask = r < 5.0 / (Z + 1) * 50
        if mask.sum() < 10:
            mask = np.ones(len(r), bool)
        r_plot = r[mask]
        d_plot = dens[mask]
        d_plot = d_plot / (d_plot.max() + 1e-15)
        ax.plot(r_plot * Z, d_plot, color=col, lw=1.8,
                label=f"Z={Z}  γ={gamma_dirac(-1,Z):.3f}")

    ax.set_xlabel("r · Z  (unidades de a₀/Z)", fontsize=10)
    ax.set_ylabel(r"$|\psi|^2$ (normalizado)", fontsize=10)
    ax.set_title("Função de onda 1s vs Z\n"
                 "Colapso: pico em r→0 quando γ→0", fontsize=10)
    ax.set_xlim(0, 3)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(lw=0.4, alpha=0.35)


def _panel_small_component(ax, Z_list=None):
    """
    Razão |F|/|G| vs Z para estado 1s.
    Para Z→Z_crit: F/G → 1 (relatividade total, sem mais 'componente pequena').
    """
    if Z_list is None:
        Z_list = np.arange(1, 130, 3)

    ratios = []
    Zs     = []
    for Z in Z_list:
        r, G, F = radial_wavefunction_strong(1, -1, Z)
        if r is None or G is None:
            continue
        maxG = np.nanmax(np.abs(G))
        maxF = np.nanmax(np.abs(F))
        if maxG > 1e-12:
            ratios.append(maxF / maxG)
            Zs.append(Z)

    # Previsão analítica: |F|/|G| ≈ Zα/(γ+κ+1) ≈ Zα para Z pequeno
    Za_arr = np.array(Zs) * ALPHA
    theory = Za_arr / (1 + np.sqrt(1 - Za_arr**2) + 1e-12)

    ax.plot(Zs, ratios,  color="#378ADD", lw=2.5, label="Numérico max|F|/max|G|")
    ax.plot(Zs, theory,  color="#D85A30", lw=1.8, ls="--", label=r"Analítico $\approx Z\alpha/(\gamma+1)$")
    ax.plot(Zs, Za_arr,  color="#1D9E75", lw=1.2, ls=":",  label=r"Limite $Z\alpha$ (linear)")
    ax.axhline(1.0, color="#EF9F27", lw=1.2, ls="--",
               label="F=G: componente 'pequena' igual à grande!")
    ax.axvline(z_critical_1(-1), color="#EF9F27", lw=0.8, ls=":", alpha=0.7)

    ax.set_xlabel("Z", fontsize=10)
    ax.set_ylabel("max|F| / max|G|", fontsize=10)
    ax.set_title("Componente pequena F(r)\n"
                 "Deixa de ser pequena quando Z→137!", fontsize=10)
    ax.set_xlim(0, 135)
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=8)
    ax.grid(lw=0.4, alpha=0.35)


def _panel_supercritical(ax, Zmax=200):
    """
    Região supercrítica Z > Z_crit:
    - Para Z_c1 < Z < Z_c2: κ imaginário, estado ressônancia
    - Para Z > Z_c2 ≈ 173: E_1s < −mc², criação de pares espontânea
    Taxa de Schwinger-like: Γ(Z) ∝ exp(−π/((Z−Z_c1)α))
    """
    Zc1 = z_critical_1(-1)       # ≈ 137
    Zc2 = 173.0                   # estimativa numérica para criação de pares

    Z_arr = np.linspace(1, Zmax, 600)

    # Energia "real" de Dirac (onde possível) + extrapolação analítica
    E_arr = []
    for Z in Z_arr:
        e = energy_dirac(1, -1, Z)
        if e is not None:
            E_arr.append(e)
        else:
            # Energia "mergualhada": extrapolação suave para −1
            E_arr.append(None)

    # Divide em subcrítico e supercrítico
    E_sub  = [e if Z <= Zc1 else np.nan for Z, e in zip(Z_arr, E_arr)]
    E_sub  = np.array([e if e is not None else np.nan for e in E_sub])

    ax.plot(Z_arr, E_sub, color="#378ADD", lw=2.5, label="E(Z) — Dirac (subcrítico)")

    # Região de colapso: Z_c1 < Z < Z_c2 — ressônancia imersa no contínuo negativo
    ax.fill_between([Zc1, Zc2], [-1.05, -1.05], [0.05, 0.05],
                    alpha=0.10, color="#EF9F27")
    ax.text((Zc1+Zc2)/2, -0.3,
            "Ressonância\n(estado quase-ligado\nno contínuo)", ha="center",
            fontsize=8, color="#EF9F27")

    # Região supercrítica Z > Z_c2: criação de pares espontânea
    ax.fill_between([Zc2, Zmax], [-1.05, -1.05], [0.05, 0.05],
                    alpha=0.10, color="#D85A30")
    ax.text((Zc2+Zmax)/2, -0.3,
            "Supercrítico:\nCriação de pares\nespontânea (vácuo\ndecai!)",
            ha="center", fontsize=8, color="#D85A30")

    # Linhas de referência
    ax.axhline(-1.0, color="#D85A30", lw=1.5, ls="--",
               label="E = −mc² (limiar)")
    ax.axhline( 0.0, color="gray",    lw=0.7, ls=":", alpha=0.5)
    ax.axvline(Zc1,  color="#EF9F27", lw=1.2, ls="--",
               label=f"Z_c1 = 1/α ≈ {Zc1:.0f}")
    ax.axvline(Zc2,  color="#D85A30", lw=1.2, ls="--",
               label=f"Z_c2 ≈ {Zc2:.0f} (criação de pares)")

    ax.text(Zc1, 0.02, f"Z_c1≈{Zc1:.0f}", fontsize=8,
            color="#EF9F27", ha="center")
    ax.text(Zc2, 0.02, f"Z_c2≈{Zc2:.0f}", fontsize=8,
            color="#D85A30", ha="center")

    ax.set_xlabel("Z", fontsize=10)
    ax.set_ylabel("E₁ₛ (mc²)", fontsize=10)
    ax.set_title("Colapso relativístico e região supercrítica\n"
                 "Z > Z_c2: vácuo torna-se instável!", fontsize=10)
    ax.set_xlim(0, Zmax)
    ax.set_ylim(-1.08, 0.06)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(lw=0.4, alpha=0.35)


def _panel_schwinger_like(ax):
    """
    Taxa de criação de pares espontânea para Z > Z_c2.
    Analogia com Schwinger: Γ ~ exp(−π m²c³/eE_nucl ħ)
    Para campo coulombiano: E_nucl ~ Ze/r² → taxa depende de (Z−Z_c2).
    """
    Zc2   = 173.0
    dZ    = np.linspace(0.01, 50, 400)
    # Aproximação da taxa: Γ ~ (dZ)^n · exp(−A/dZ)
    # Referência: Greiner et al., "QED of Strong Fields", cap. 7
    A_fit = 5.0
    Gamma = (dZ / Zc2)**2 * np.exp(-A_fit / (dZ / Zc2 + 1e-12))
    Gamma /= Gamma.max()

    ax.semilogy(Zc2 + dZ, Gamma, color="#D85A30", lw=2.5,
                label=r"$\Gamma \propto \left(\frac{Z-Z_{c2}}{Z_{c2}}\right)^2 e^{-A/(Z-Z_{c2})}$")
    ax.axvline(Zc2, color="#EF9F27", lw=1.2, ls="--",
               label=f"Z_c2 ≈ {Zc2:.0f}")

    # Núcleos reais mais pesados
    for Znuc, nome in [(92,"U"), (94,"Pu"), (118,"Og")]:
        if Znuc > Zc2:
            ax.axvline(Znuc, color="#1D9E75", lw=0.8, ls=":")
            ax.text(Znuc+0.5, 0.5, nome, fontsize=8, color="#1D9E75")

    ax.set_xlabel("Z", fontsize=10)
    ax.set_ylabel(r"Taxa $\Gamma$ (normaliz.)", fontsize=10)
    ax.set_title("Taxa de criação de pares espontânea\n"
                 "(vácuo de QED decai para Z > Z_c2 ≈ 173)", fontsize=10)
    ax.set_xlim(Zc2 - 5, Zc2 + 55)
    ax.legend(fontsize=8)
    ax.grid(lw=0.4, alpha=0.35)


def _panel_energy_table(ax, Z_vals=None):
    """
    Tabela visual: E_Dirac vs E_Schrödinger vs E_KG para Z crescente.
    Klein-Gordon: E_KG = mc²√(1−(Zα/n)²) − mc²  (spin 0)
    """
    if Z_vals is None:
        Z_vals = [1, 10, 30, 50, 80, 100, 120]

    E_D, E_S, E_KG = [], [], []
    for Z in Z_vals:
        ed  = energy_dirac(1, -1, Z)
        es  = -(Z**2) / (2) * ALPHA**2  # 1s Schrödinger em u.a.
        # Klein-Gordon: n=1
        Za  = Z * ALPHA
        if Za < 1.0:
            ekg = np.sqrt(1.0 - Za**2) - 1.0
        else:
            ekg = np.nan
        E_D.append(ed  if ed  is not None else np.nan)
        E_S.append(es)
        E_KG.append(ekg)

    x = np.arange(len(Z_vals))
    w = 0.26
    ax.bar(x - w, E_D,  w, color="#378ADD", label="Dirac (spin ½)")
    ax.bar(x,     E_KG, w, color="#1D9E75", label="Klein-Gordon (spin 0)", alpha=0.85)
    ax.bar(x + w, E_S,  w, color="#D3D1C7", label="Schrödinger")

    ax.axhline(-1.0, color="#D85A30", lw=1.2, ls="--", alpha=0.7,
               label="E = −mc²")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Z={z}" for z in Z_vals], fontsize=8, rotation=30)
    ax.set_ylabel("Energia (mc²)", fontsize=10)
    ax.set_title("Energia 1s: Dirac vs Klein-Gordon vs Schrödinger\n"
                 "Schrödinger diverge; Dirac e KG colapsam em Z→137", fontsize=9)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(axis='y', lw=0.4, alpha=0.35)


# ============================================================================
# Figura estática completa
# ============================================================================

def plot_all_static(args):
    Zmax = args["Zmax"]

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(
        "Campo Forte e Colapso Relativístico — Equação de Dirac\n"
        r"Zα → 1: $\gamma\to0$, $\psi\sim r^{\gamma-1}$ diverge, vácuo torna-se instável",
        fontsize=12
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.40)

    _panel_energy_vs_Z(fig.add_subplot(gs[0,0]), Zmax=min(Zmax,130))
    _panel_gamma_vs_Z (fig.add_subplot(gs[0,1]))
    _panel_wavefunction_comparison(fig.add_subplot(gs[0,2]))
    _panel_small_component(fig.add_subplot(gs[1,0]))
    _panel_supercritical(fig.add_subplot(gs[1,1]), Zmax=Zmax)
    _panel_schwinger_like(fig.add_subplot(gs[1,2]))

    plt.savefig("dirac_strong.png", dpi=120, bbox_inches="tight")
    print("Figura salva: dirac_strong.png")
    plt.show()


def plot_energy_table(args):
    fig, ax = plt.subplots(figsize=(10, 5))
    _panel_energy_table(ax)
    fig.tight_layout()
    plt.savefig("dirac_strong_table.png", dpi=120, bbox_inches="tight")
    print("Figura salva: dirac_strong_table.png")
    plt.show()


# ============================================================================
# Tabela numérica
# ============================================================================

def print_summary(args):
    print(f"\n{'='*70}")
    print("  Campo Forte e Colapso Relativístico")
    print(f"{'='*70}")
    print(f"  α = {ALPHA:.8f}   1/α = {1/ALPHA:.4f}")
    print(f"  Z_crit (1s, κ=−1) = 1/α ≈ {z_critical_1(-1):.2f}")
    print(f"  Z_crit (2p₁/₂, κ=+1) = 1/α ≈ {z_critical_1(1):.2f}  [mesma!]")
    print(f"  Z_crit (2p₃/₂, κ=−2) = 2/α ≈ {z_critical_1(-2):.2f}")
    print()
    print(f"  {'Z':>5}  {'γ(1s)':>8}  {'E_Dirac(mc²)':>14}  {'E_Schrö(mc²)':>14}  "
          f"{'|F|/|G|':>8}  {'ΔE/E_S':>10}")
    print("  " + "─"*65)
    for Z in [1, 10, 30, 50, 70, 90, 110, 130, 136]:
        g   = gamma_dirac(-1, Z)
        ed  = energy_dirac(1, -1, Z)
        es  = -(Z**2 * ALPHA**2) / 2
        if ed is None:
            print(f"  {Z:>5}  {'COLAPSO':>8}")
            continue
        r, G, F = radial_wavefunction_strong(1, -1, Z, Nr=2000)
        rat = (np.nanmax(np.abs(F)) / (np.nanmax(np.abs(G))+1e-15)
               if r is not None else 0.0)
        rel_diff = (ed - es) / (abs(es) + 1e-15)
        print(f"  {Z:>5}  {g:>8.5f}  {ed:>14.7f}  {es:>14.7f}  "
              f"{rat:>8.4f}  {rel_diff:>10.4f}")
    print()
    print("  Marcos do colapso:")
    print(f"    γ = 1.0:  Z=0         (livre, sem campo)")
    print(f"    γ = 0.5:  Z≈{np.sqrt(1-0.25)/ALPHA:.0f}        (ψ ~ r^{{-0.5}}: limite normalizável)")
    print(f"    γ = 0.0:  Z≈{1/ALPHA:.0f}  (COLAPSO: ψ diverge na origem)")
    print(f"    E=−mc²:  Z≈173      (criação de pares espontânea)")
    print(f"{'='*70}\n")


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Campo Forte e Colapso Relativístico — Equação de Dirac",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--Zmax",   type=int,   default=200,
                   help="Z máximo para região supercrítica")
    p.add_argument("--module", type=int,   default=0,
                   help="0=tudo 1=espectro+colapso 2=tabela")
    return vars(p.parse_args())


if __name__ == "__main__":
    args = parse_args()
    print_summary(args)
    mod = args["module"]
    if mod == 0 or mod == 1:
        plot_all_static(args)
    if mod == 0 or mod == 2:
        plot_energy_table(args)
