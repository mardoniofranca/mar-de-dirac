"""
Simulação da Equação de Dirac — Átomo de Hidrogênio Relativístico
==================================================================
Resolve o hamiltoniano de Dirac com potencial coulombiano em coordenadas
esféricas, usando integração numérica RK4 com truncamento adaptativo.

Sistema de unidades atômicas: ħ = m_e = e = 1, c = 1/α ≈ 137

Uso:
    python dirac_hydrogen.py
    python dirac_hydrogen.py --nmax 4
    python dirac_hydrogen.py --Z 6     (carbono hidrogenóide)
    python dirac_hydrogen.py --help
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import platform

if platform.system() == "Darwin":
    matplotlib.use("MacOSX")

# ---------------------------------------------------------------------------
# Constantes (unidades atômicas)
# ---------------------------------------------------------------------------
ALPHA    = 1.0 / 137.035999084
C_AU     = 1.0 / ALPHA          # c ≈ 137 u.a.
HART_EV  = 27.211396            # 1 Hartree em eV


# ---------------------------------------------------------------------------
# Energias analíticas
# ---------------------------------------------------------------------------

def energy_dirac(n, kappa, Z=1):
    """Energia exata de Dirac (u.a.), energia de repouso subtraída."""
    Za    = Z * ALPHA
    gamma = np.sqrt(kappa ** 2 - Za ** 2)
    nr    = n - abs(kappa)
    E_tot = C_AU ** 2 / np.sqrt(1.0 + (Za / (nr + gamma)) ** 2)
    return E_tot - C_AU ** 2


def energy_schrodinger(n, Z=1):
    return -Z ** 2 / (2.0 * n ** 2)


def energy_sommerfeld(n, l, j, Z=1):
    E0   = energy_schrodinger(n, Z)
    corr = (Z * ALPHA) ** 2 / n ** 2 * (n / (j + 0.5) - 0.75)
    return E0 * (1.0 + corr)


# ---------------------------------------------------------------------------
# Números quânticos
# ---------------------------------------------------------------------------

def kappa_states(n):
    """Todos os estados (n, κ, l, j, label) para um dado n."""
    states = []
    for kappa in list(range(-n, 0)) + list(range(1, n)):
        nr = n - abs(kappa)
        if nr < 0:
            continue
        l = (-kappa - 1) if kappa < 0 else kappa
        j = abs(kappa) - 0.5
        lname = ["s", "p", "d", "f", "g", "h"][min(l, 5)]
        label = f"{n}{lname}_{int(2*j)}/2  (κ={kappa:+d})"
        states.append((n, kappa, l, j, label))
    # ordena por energia Dirac (crescente em módulo → decrescente em E)
    states.sort(key=lambda s: energy_dirac(s[0], s[1]))
    return states


# ---------------------------------------------------------------------------
# Função de onda radial — RK4 com truncamento
# ---------------------------------------------------------------------------

def radial_wavefunction(n, kappa, Z=1, Nr=3000):
    """
    Integra G(r) e F(r) via RK4, truncando quando a solução diverge.

    Equações de Dirac radiais (u.a.):
        dG/dr = -(κ/r) G + [(E + c² - V)/c] F
        dF/dr =  (κ/r) F - [(E - c² - V)/c] G    (V = -Z/r < 0)
    """
    E  = energy_dirac(n, kappa, Z)
    c  = C_AU
    Za = Z * ALPHA

    # Grade logarítmica: mais densa perto da origem
    r_min = 1e-5
    r_max = max(80.0, 10 * n ** 2)
    r = np.geomspace(r_min, r_max, Nr)

    G = np.zeros(Nr)
    F = np.zeros(Nr)

    # Condição inicial: solução regular ~r^γ perto da origem
    gamma = np.sqrt(kappa ** 2 - Za ** 2)
    A = 1.0
    # razão F/G da solução regular de Dirac na origem
    if abs(gamma + kappa) > 1e-12:
        B = Za * c / (gamma + kappa) * A
    else:
        B = 0.0

    G[0] = A * r[0] ** gamma
    F[0] = B * r[0] ** gamma

    def derivs(ri, g, f):
        V  = -Z / ri
        dg = -(kappa / ri) * g + ((E + c**2 - V) / c) * f
        df =  (kappa / ri) * f - ((E - c**2 - V) / c) * g
        return dg, df

    cutoff = Nr  # índice onde truncamos (função decaiu e voltou a crescer)
    prev_dens = G[0]**2 + F[0]**2

    for i in range(Nr - 1):
        ri = r[i]
        gi, fi = G[i], F[i]
        dr = r[i+1] - r[i]

        k1g, k1f = derivs(ri,        gi,            fi)
        k2g, k2f = derivs(ri + dr/2, gi + dr/2*k1g, fi + dr/2*k1f)
        k3g, k3f = derivs(ri + dr/2, gi + dr/2*k2g, fi + dr/2*k2f)
        k4g, k4f = derivs(ri + dr,   gi + dr*k3g,   fi + dr*k3f)

        G[i+1] = gi + dr/6 * (k1g + 2*k2g + 2*k3g + k4g)
        F[i+1] = fi + dr/6 * (k1f + 2*k2f + 2*k3f + k4f)

        curr_dens = G[i+1]**2 + F[i+1]**2

        # Trunca se a função começar a divergir (cresce depois de decair)
        if curr_dens > prev_dens * 10 and ri > 2.0 * n:
            cutoff = i + 1
            G[cutoff:] = 0.0
            F[cutoff:] = 0.0
            break

        if np.isnan(curr_dens) or np.isinf(curr_dens):
            cutoff = i
            G[cutoff:] = 0.0
            F[cutoff:] = 0.0
            break

        if curr_dens < prev_dens:
            prev_dens = curr_dens  # atualiza só quando decrescente

    # Normaliza ∫(G²+F²) dr = 1  (com medida dr da grade log)
    dr_arr = np.diff(r, prepend=r[0])
    norm = np.sqrt(np.sum((G**2 + F**2) * dr_arr))
    if norm > 1e-12:
        G /= norm
        F /= norm

    return r, G, F


# ---------------------------------------------------------------------------
# Tabela de energias
# ---------------------------------------------------------------------------

def print_table(nmax, Z=1):
    print(f"\n  Átomo de Hidrogênio Relativístico  (Z={Z}, α={ALPHA:.6f})")
    print(f"\n  {'Estado':<18} {'E_Dirac':>15} {'E_Schrö':>15} {'ΔE (u.a.)':>13} {'ΔE (eV)':>12}")
    print("  " + "─" * 75)
    for n in range(1, nmax + 1):
        for (nn, kappa, l, j, label) in kappa_states(n):
            Ed   = energy_dirac(nn, kappa, Z)
            Es   = energy_schrodinger(nn, Z)
            dE   = Ed - Es
            short = label.split("(")[0].strip()
            print(f"  {short:<18} {Ed:>15.8f} {Es:>15.8f} {dE:>13.3e} {dE*HART_EV:>12.6f}")
        print()
    print("  ΔE = E_Dirac − E_Schrödinger (correção relativística em eV)")
    print()


# ---------------------------------------------------------------------------
# Figura principal
# ---------------------------------------------------------------------------

def plot_all(nmax, Z=1):
    fig = plt.figure(figsize=(15, 9))
    fig.suptitle(
        f"Átomo de Hidrogênio Relativístico — Equação de Dirac  (Z={Z})",
        fontsize=13
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.38)

    ax_spec  = fig.add_subplot(gs[0, :2])
    ax_split = fig.add_subplot(gs[0, 2])
    ax_1s    = fig.add_subplot(gs[1, 0])
    ax_n2    = fig.add_subplot(gs[1, 1])
    ax_small = fig.add_subplot(gs[1, 2])

    # ── Painel 1: Espectro completo ──────────────────────────────────────
    cmap_n = {1:"#378ADD", 2:"#1D9E75", 3:"#EF9F27", 4:"#D85A30", 5:"#7F77DD"}
    xtick_labels, xtick_pos = [], []
    x = 0
    for n in range(1, nmax + 1):
        for (nn, kappa, l, j, label) in kappa_states(n):
            Ed = energy_dirac(nn, kappa, Z) * 1e3
            Es = energy_schrodinger(nn, Z) * 1e3
            c_n = cmap_n.get(n, "#888780")
            short = label.split("(")[0].strip()
            ax_spec.plot([x, x], [Es, Ed], color=c_n, lw=1.0, alpha=0.5)
            ax_spec.scatter(x, Ed, color=c_n, s=60, zorder=5)
            ax_spec.scatter(x, Es, color=c_n, s=60, marker="x", zorder=4, alpha=0.5)
            xtick_labels.append(short)
            xtick_pos.append(x)
            x += 1
        x += 0.5   # espaço entre grupos n

    ax_spec.set_xticks(xtick_pos)
    ax_spec.set_xticklabels(xtick_labels, rotation=40, fontsize=7.5, ha="right")
    ax_spec.set_ylabel("Energia (×10⁻³ u.a.)", fontsize=10)
    ax_spec.set_title("● Dirac (exato)    × Schrödinger    linha = correção relativística", fontsize=9)
    ax_spec.grid(axis='y', lw=0.4, alpha=0.4)

    from matplotlib.lines import Line2D
    leg = [Line2D([0],[0], color=cmap_n.get(n,"#888780"), lw=2, label=f"n={n}")
           for n in range(1, nmax+1)]
    ax_spec.legend(handles=leg, fontsize=8, loc="lower right")

    # ── Painel 2: Splitting n=2 (barras) ────────────────────────────────
    st2 = kappa_states(2)
    xlbls = [s[4].split("(")[0].strip() for s in st2]
    Ed2   = np.array([energy_dirac(s[0],s[1],Z) for s in st2]) * 1e6
    Es2   = np.array([energy_schrodinger(s[0],Z) for s in st2]) * 1e6
    Eso2  = np.array([energy_sommerfeld(s[0],s[2],s[3],Z) for s in st2]) * 1e6
    xp    = np.arange(len(xlbls))
    w     = 0.26
    ax_split.bar(xp-w,  Ed2,  w, color="#378ADD", label="Dirac")
    ax_split.bar(xp,    Eso2, w, color="#1D9E75", label="Sommerfeld", alpha=0.85)
    ax_split.bar(xp+w,  Es2,  w, color="#D3D1C7", label="Schrödinger")
    ax_split.set_xticks(xp)
    ax_split.set_xticklabels(xlbls, fontsize=7.5, rotation=25)
    ax_split.set_ylabel("E (×10⁻⁶ u.a.)", fontsize=9)
    ax_split.set_title("Splitting spin-órbita em n=2", fontsize=10)
    ax_split.legend(fontsize=7.5)
    ax_split.grid(axis='y', lw=0.4, alpha=0.4)

    # ── Painel 3: Função de onda 1s ──────────────────────────────────────
    r1, G1, F1 = radial_wavefunction(1, -1, Z)
    dens1 = G1**2 + F1**2
    ax_1s.plot(r1, dens1,  color="#378ADD", lw=2,   label=r"$|G|^2+|F|^2$")
    ax_1s.plot(r1, G1**2,  color="#378ADD", lw=1.2, ls="--", alpha=0.55, label=r"$|G|^2$")
    ax_1s.set_xlim(0, 6); ax_1s.set_ylim(bottom=0)
    ax_1s.set_xlabel("r (a₀)", fontsize=9)
    ax_1s.set_ylabel("Densidade radial", fontsize=9)
    ax_1s.set_title(r"Densidade radial: 1s$_{1/2}$", fontsize=10)
    ax_1s.legend(fontsize=8); ax_1s.grid(lw=0.3, alpha=0.4)

    # ── Painel 4: n=2 — splitting nas funções de onda ────────────────────
    pairs = [(-1,"#1D9E75",r"2s$_{1/2}$  (κ=−1)"),
             ( 1,"#EF9F27",r"2p$_{1/2}$  (κ=+1)"),
             (-2,"#D85A30",r"2p$_{3/2}$  (κ=−2)", "--")]
    for item in pairs:
        k, col, lbl = item[0], item[1], item[2]
        ls = item[3] if len(item) > 3 else "-"
        rr, GG, FF = radial_wavefunction(2, k, Z)
        ax_n2.plot(rr, GG**2+FF**2, color=col, lw=2, ls=ls, label=lbl)
    ax_n2.set_xlim(0, 28); ax_n2.set_ylim(bottom=0)
    ax_n2.set_xlabel("r (a₀)", fontsize=9)
    ax_n2.set_title("Estados n=2: splitting spin-órbita", fontsize=10)
    ax_n2.legend(fontsize=8); ax_n2.grid(lw=0.3, alpha=0.4)

    # ── Painel 5: Componente pequena F ───────────────────────────────────
    ax_small.plot(r1, np.abs(G1), color="#378ADD", lw=2,   label=r"$|G|$ — grande")
    ax_small.plot(r1, np.abs(F1), color="#D85A30", lw=2,   label=r"$|F|$ — pequena")
    ratio = np.nanmax(np.abs(F1)) / (np.nanmax(np.abs(G1)) + 1e-15)
    ax_small.set_xlim(0, 6); ax_small.set_ylim(bottom=0)
    ax_small.set_xlabel("r (a₀)", fontsize=9)
    ax_small.set_ylabel("Amplitude", fontsize=9)
    ax_small.set_title(r"Componente pequena $F(r)$  [efeito relativístico]", fontsize=10)
    ax_small.text(0.97, 0.88,
                  f"max|F|/max|G|\n≈ {ratio:.5f}\n≈ Zα = {Z*ALPHA:.5f}",
                  transform=ax_small.transAxes, fontsize=8.5,
                  ha="right", va="top", color="#D85A30")
    ax_small.legend(fontsize=8); ax_small.grid(lw=0.3, alpha=0.4)

    plt.savefig("dirac_hydrogen.png", dpi=150, bbox_inches="tight")
    print("Figura salva: dirac_hydrogen.png")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Dirac — Átomo de Hidrogênio Relativístico",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--nmax", type=int, default=3, help="Número quântico máximo n")
    p.add_argument("--Z",    type=int, default=1, help="Número atômico (1=H, 2=He+, 6=C5+...)")
    return vars(p.parse_args())


if __name__ == "__main__":
    args = parse_args()
    print_table(args["nmax"], args["Z"])
    plot_all(args["nmax"], args["Z"])
