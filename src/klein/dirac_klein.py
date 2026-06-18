"""
Simulação do Efeito Klein — Equação de Dirac 1D
================================================
Um elétron relativístico incide sobre uma barreira de potencial retangular.

Três regimes físicos distintos:
  1. V₀ < E + mc²        → tunelamento exponencial (como Schrödinger)
  2. E + mc² < V₀ < E + mc² + 2mc²  → região de oscilação (ressonâncias)
  3. V₀ > E + 2mc²       → PARADOXO DE KLEIN: transmissão T → 1 (100%)

Física implementada:
  - Coeficientes de transmissão T(k) e reflexão R(k) analíticos via Dirac
  - Comparação com a previsão de Schrödinger (sem relatividade)
  - Evolução temporal do pacote de onda gaussiano atravessando a barreira
  - Visualização da densidade |ψ|² e das componentes do espinor

Sistema de unidades naturais: ħ = c = m = 1

Uso:
    python dirac_klein.py               # roda tudo
    python dirac_klein.py --V0 3.0      # muda altura da barreira
    python dirac_klein.py --static      # só gráficos estáticos (sem animação)
    python dirac_klein.py --save        # salva animação como GIF
    python dirac_klein.py --help
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import argparse
import platform
import warnings
warnings.filterwarnings("ignore")

if platform.system() == "Darwin":
    matplotlib.use("MacOSX")

# ---------------------------------------------------------------------------
# Coeficientes analíticos de transmissão — Dirac vs Schrödinger
# ---------------------------------------------------------------------------

def dirac_transmission(k, V0, d, m=1.0):
    """
    Transmissão T(k) via equação de Dirac para barreira retangular [0, d].

    Dentro da barreira, o momento é:
        q = √( (E - V₀ + m)·(E - V₀ - m) )   com E = √(k²+m²)

    Quando V₀ > E + m = E + mc²:  (E-V₀+m) < 0 e (E-V₀-m) < 0
        → produto positivo → q real → ondas propagantes dentro da barreira!
    Isso é a origem do paradoxo de Klein.
    """
    E  = np.sqrt(k**2 + m**2)
    Ev = E - V0

    # momento dentro da barreira
    arg = (Ev + m) * (Ev - m)
    # arg pode ser negativo (barreira clássica) ou positivo (Klein)
    q = np.sqrt(np.abs(arg) + 0j)
    if arg < 0:
        q = 1j * np.abs(q)  # imaginário → decaimento evanescente

    # spinores — fator de helicidade
    alpha_k = k  / (E + m)       # fora da barreira
    alpha_q = Ev / (q  + 0j)     # dentro  (pode ser imaginário)
    # (usamos a razão F/G do espinor de onda plana)

    # Matriz de transferência para barreira retangular
    cos_qd = np.cos(q * d)
    sin_qd = np.sin(q * d)

    # Elemento da matriz de transmissão (ver Calogeracos & Dombey 1999)
    denom = (cos_qd
             - 0.5j * (alpha_q / alpha_k + alpha_k / alpha_q) * sin_qd)

    T = 1.0 / (np.abs(denom)**2)
    T = np.real(T)
    return np.clip(T, 0.0, 1.0)


def schrodinger_transmission(k, V0, d, m=1.0):
    """
    Transmissão de Schrödinger para barreira retangular.
    E = k²/(2m), momento dentro = √(2m(E-V₀)).
    """
    E  = k**2 / (2 * m)
    Ev = E - V0
    if Ev >= 0:
        q = np.sqrt(2 * m * Ev + 0j)
    else:
        q = 1j * np.sqrt(2 * m * np.abs(Ev))

    cos_qd = np.cos(q * d)
    sin_qd = np.sin(q * d)
    denom  = cos_qd - 0.5j * (q/k + k/q) * sin_qd
    T = 1.0 / (np.abs(denom)**2 + 1e-15)
    return float(np.clip(np.real(T), 0.0, 1.0))


# ---------------------------------------------------------------------------
# Evolução temporal — split-operator FFT
# ---------------------------------------------------------------------------

def init_spinor(xs, dx, k0, sigma, m=1.0):
    """Pacote gaussiano de energia positiva."""
    env  = np.exp(-0.25 * (xs / sigma)**2) / (np.pi * sigma**2)**0.25
    psi1 = env * np.exp(1j * k0 * xs)
    w0   = np.sqrt(k0**2 + m**2)
    psi2 = (k0 / (w0 + m)) * psi1
    norm = np.sqrt(np.sum((np.abs(psi1)**2 + np.abs(psi2)**2) * dx))
    return psi1/norm, psi2/norm


def make_propagators_free(ks, m, dt):
    """Propagador livre (sem barreira) no espaço de momentos."""
    omega  = np.sqrt(ks**2 + m**2)
    cos_w  = np.cos(omega * dt)
    sinc_w = np.where(omega > 1e-12, np.sin(omega * dt) / omega, dt)
    P11 =  cos_w - 1j * m * sinc_w
    P12 = -1j * ks * sinc_w
    P21 = -1j * ks * sinc_w
    P22 =  cos_w + 1j * m * sinc_w
    return P11, P12, P21, P22


def make_propagator_potential(V_arr, dt):
    """
    Propagador do potencial escalar no espaço de posições:
        exp(-i V(x) σ₃ dt / 2)   (half-step, Strang splitting)
    Para potencial escalar V: afeta apenas a energia, não acopla componentes.
        psi1 → psi1 · exp(-i V dt/2)
        psi2 → psi2 · exp(+i V dt/2)
    """
    phase = np.exp(-1j * V_arr * dt / 2)
    return phase, np.conj(phase)


def step_with_potential(psi1, psi2, P11, P12, P21, P22, phase1, phase2):
    """
    Passo de Strang splitting:
        1. Meio passo de potencial
        2. Passo completo cinético (FFT)
        3. Meio passo de potencial
    """
    # meio passo potencial
    psi1 = psi1 * phase1
    psi2 = psi2 * phase2

    # passo cinético
    phi1 = np.fft.fft(psi1)
    phi2 = np.fft.fft(psi2)
    psi1 = np.fft.ifft(P11 * phi1 + P12 * phi2)
    psi2 = np.fft.ifft(P21 * phi1 + P22 * phi2)

    # meio passo potencial
    psi1 = psi1 * phase1
    psi2 = psi2 * phase2

    return psi1, psi2


# ---------------------------------------------------------------------------
# Figura estática: T(V₀) e T(k)
# ---------------------------------------------------------------------------

def plot_static(args):
    m  = args["mass"]
    d  = args["d"]
    k0 = args["k0"]
    V0 = args["V0"]

    fig = plt.figure(figsize=(14, 8))
    fig.suptitle("Paradoxo de Klein — Transmissão pela Barreira de Potencial", fontsize=13)
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)

    ax_v   = fig.add_subplot(gs[0, :2])   # T vs V₀
    ax_k   = fig.add_subplot(gs[0, 2])    # T vs k
    ax_psi = fig.add_subplot(gs[1, :])    # |ψ|² em x para 3 regimes

    E0 = np.sqrt(k0**2 + m**2)  # energia cinética do elétron

    # ── T vs V₀ ──────────────────────────────────────────────────────────
    V0_arr = np.linspace(0, 6 * m, 600)

    T_dirac = np.array([dirac_transmission(k0, v, d, m) for v in V0_arr])
    T_schr  = np.array([schrodinger_transmission(k0, v, d, m) for v in V0_arr])

    ax_v.plot(V0_arr / m, T_dirac, color="#378ADD", lw=2.0, label="Dirac (relativístico)")
    ax_v.plot(V0_arr / m, T_schr,  color="#D85A30", lw=2.0, ls="--", label="Schrödinger")

    # linhas de referência dos 3 regimes
    E_plus_m  = E0 + m
    E_plus_3m = E0 + 3 * m
    ax_v.axvline(E_plus_m  / m, color="#1D9E75", lw=1.2, ls=":", label=f"V₀ = E+mc²  ({E_plus_m/m:.2f}mc²)")
    ax_v.axvline(E_plus_3m / m, color="#EF9F27", lw=1.2, ls=":", label=f"V₀ = E+3mc² ({E_plus_3m/m:.2f}mc²)")
    ax_v.axvline(V0 / m,        color="gray",    lw=0.8, ls="--", alpha=0.6, label=f"V₀ selecionado ({V0/m:.1f}mc²)")

    ax_v.set_xlabel("V₀ / mc²", fontsize=10)
    ax_v.set_ylabel("Transmissão T", fontsize=10)
    ax_v.set_title(f"T(V₀)  — k₀={k0}, d={d}, m={m}", fontsize=10)
    ax_v.set_ylim(-0.05, 1.15)
    ax_v.legend(fontsize=8.5, loc="upper right")
    ax_v.grid(lw=0.4, alpha=0.4)

    # Anotação dos 3 regimes
    yann = 1.07
    ax_v.annotate("",  xy=(E_plus_m/m, yann), xytext=(0, yann),
                  arrowprops=dict(arrowstyle="<->", color="#378ADD", lw=1.2))
    ax_v.text(E_plus_m / m / 2, yann + 0.03, "Tunelamento\n(como Schrö.)",
              ha="center", fontsize=7.5, color="#378ADD")

    ax_v.annotate("", xy=(E_plus_3m/m, yann), xytext=(E_plus_m/m, yann),
                  arrowprops=dict(arrowstyle="<->", color="#1D9E75", lw=1.2))
    ax_v.text((E_plus_m + E_plus_3m) / (2*m), yann + 0.03, "Oscilação /\nressonâncias",
              ha="center", fontsize=7.5, color="#1D9E75")

    ax_v.annotate("", xy=(V0_arr[-1]/m, yann), xytext=(E_plus_3m/m, yann),
                  arrowprops=dict(arrowstyle="<->", color="#EF9F27", lw=1.2))
    ax_v.text((E_plus_3m/m + V0_arr[-1]/m)/2, yann+0.03, "Klein:\nT → 1",
              ha="center", fontsize=7.5, color="#EF9F27")

    # ── T vs k ───────────────────────────────────────────────────────────
    k_arr  = np.linspace(0.1, 5 * m, 400)
    Tk_d   = np.array([dirac_transmission(k, V0, d, m) for k in k_arr])
    Tk_s   = np.array([schrodinger_transmission(k, V0, d, m) for k in k_arr])

    ax_k.plot(k_arr / m, Tk_d, color="#378ADD", lw=2, label="Dirac")
    ax_k.plot(k_arr / m, Tk_s, color="#D85A30", lw=2, ls="--", label="Schrödinger")
    ax_k.axvline(k0 / m, color="gray", lw=0.8, ls="--", alpha=0.7, label=f"k₀={k0}")
    ax_k.set_xlabel("k / mc", fontsize=10)
    ax_k.set_ylabel("Transmissão T", fontsize=10)
    ax_k.set_title(f"T(k)  — V₀={V0}mc²", fontsize=10)
    ax_k.set_ylim(-0.05, 1.15)
    ax_k.legend(fontsize=8.5)
    ax_k.grid(lw=0.4, alpha=0.4)

    # ── |ψ|² para 3 regimes ──────────────────────────────────────────────
    N     = 2048
    L     = 80.0
    xs    = np.linspace(-L/2, L/2, N, endpoint=False)
    dx    = xs[1] - xs[0]
    sigma = 5.0
    x0    = -L/4   # posição inicial do pacote

    V_regimes = [
        (0.5 * (E0 + m),   "#378ADD", f"Regime 1: V₀={0.5*(E0+m):.2f} < E+mc²  (tunelamento)"),
        (1.5 * E0,         "#1D9E75", f"Regime 2: V₀={1.5*E0:.2f} ≈ E+mc²  (ressonância)"),
        (V0,               "#EF9F27", f"Regime 3: V₀={V0:.2f} > E+3mc²  (Klein, T→1)"),
    ]

    x_bar_l = -2.0
    x_bar_r =  2.0

    for V0_r, col, lbl in V_regimes:
        V_arr = np.where((xs >= x_bar_l) & (xs <= x_bar_r), V0_r, 0.0)
        ks    = np.fft.fftfreq(N, d=dx) * 2 * np.pi
        dt    = 0.04
        P11, P12, P21, P22 = make_propagators_free(ks, m, dt)
        ph1, ph2            = make_propagator_potential(V_arr, dt)

        psi1, psi2 = init_spinor(xs - x0, dx, k0, sigma, m)

        # evolui até o elétron ter cruzado a barreira
        n_steps = int(L / (2 * k0 / np.sqrt(k0**2 + m**2)) / dt)
        for _ in range(n_steps):
            psi1, psi2 = step_with_potential(psi1, psi2, P11, P12, P21, P22, ph1, ph2)

        dens = (np.abs(psi1)**2 + np.abs(psi2)**2).real
        dens /= np.max(dens) + 1e-15   # normaliza para comparação visual

        ax_psi.plot(xs, dens, color=col, lw=1.8, alpha=0.85, label=lbl)

    # barreira
    bar_h = 0.95
    ax_psi.fill_betweenx([0, bar_h], x_bar_l, x_bar_r,
                         alpha=0.10, color="gray")
    ax_psi.axvline(x_bar_l, color="gray", lw=0.8, ls=":")
    ax_psi.axvline(x_bar_r, color="gray", lw=0.8, ls=":")
    ax_psi.text((x_bar_l+x_bar_r)/2, bar_h * 0.5, "Barreira",
                ha="center", fontsize=8, color="gray")

    ax_psi.set_xlabel("x (unidades naturais)", fontsize=10)
    ax_psi.set_ylabel(r"$|\psi|^2$ (normalizado)", fontsize=10)
    ax_psi.set_title("Densidade de probabilidade após atravessar a barreira — 3 regimes", fontsize=10)
    ax_psi.legend(fontsize=8.5, loc="upper right")
    ax_psi.set_xlim(-L/2, L/2)
    ax_psi.set_ylim(bottom=0)
    ax_psi.grid(lw=0.4, alpha=0.4)

    plt.savefig("dirac_klein_static.png", dpi=150, bbox_inches="tight")
    print("Figura salva: dirac_klein_static.png")
    plt.show()


# ---------------------------------------------------------------------------
# Animação: pacote cruzando a barreira
# ---------------------------------------------------------------------------

def animate_klein(args):
    m    = args["mass"]
    k0   = args["k0"]
    V0   = args["V0"]
    N    = 2048
    L    = 100.0
    xs   = np.linspace(-L/2, L/2, N, endpoint=False)
    dx   = xs[1] - xs[0]
    ks   = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    dt   = 0.04
    sigma = 6.0
    x0    = -L / 3.5

    x_bar_l, x_bar_r = -2.0, 2.0
    V_arr = np.where((xs >= x_bar_l) & (xs <= x_bar_r), V0, 0.0)

    P11, P12, P21, P22 = make_propagators_free(ks, m, dt)
    ph1, ph2            = make_propagator_potential(V_arr, dt)

    psi1, psi2 = init_spinor(xs - x0, dx, k0, sigma, m)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    E0 = np.sqrt(k0**2 + m**2)
    fig.suptitle(
        f"Paradoxo de Klein  |  k₀={k0}, m={m}, V₀={V0:.1f}mc²  "
        f"(E+mc²={E0+m:.2f})  →  T_analítico={dirac_transmission(k0,V0,x_bar_r-x_bar_l,m):.3f}",
        fontsize=11
    )

    dens0 = (np.abs(psi1)**2 + np.abs(psi2)**2).real
    ymax  = dens0.max() * 2.2

    line_dens, = ax_top.plot(xs, dens0, color="#378ADD", lw=2, label=r"$|\psi|^2$")
    ax_top.fill_betweenx([0, ymax], x_bar_l, x_bar_r, alpha=0.12, color="#D85A30")
    ax_top.axvline(x_bar_l, color="#D85A30", lw=0.9, ls=":")
    ax_top.axvline(x_bar_r, color="#D85A30", lw=0.9, ls=":")
    ax_top.text(0, ymax*0.88, f"V₀={V0:.1f}", ha="center", fontsize=9, color="#D85A30")
    ax_top.set_ylim(0, ymax)
    ax_top.set_ylabel(r"$|\psi|^2$", fontsize=11)
    ax_top.legend(fontsize=9, loc="upper right")
    ax_top.grid(lw=0.4, alpha=0.4)
    txt_t = ax_top.text(0.02, 0.90, "t = 0.00", transform=ax_top.transAxes, fontsize=9)
    txt_T = ax_top.text(0.02, 0.78, "", transform=ax_top.transAxes, fontsize=9, color="#1D9E75")

    lim_sp = max(np.abs(psi1).max(), np.abs(psi2).max()) * 1.5
    line_re1, = ax_bot.plot(xs, psi1.real, color="#378ADD", lw=1.5, alpha=0.8, label=r"Re $\psi_1$")
    line_re2, = ax_bot.plot(xs, psi2.real, color="#EF9F27", lw=1.5, alpha=0.8, label=r"Re $\psi_2$")
    ax_bot.fill_betweenx([-lim_sp, lim_sp], x_bar_l, x_bar_r, alpha=0.10, color="#D85A30")
    ax_bot.axhline(0, color="gray", lw=0.5, ls="--")
    ax_bot.set_ylim(-lim_sp, lim_sp)
    ax_bot.set_xlabel("x (unidades naturais)", fontsize=11)
    ax_bot.set_ylabel("Amplitude", fontsize=11)
    ax_bot.legend(fontsize=9, loc="upper right")
    ax_bot.grid(lw=0.4, alpha=0.4)
    ax_bot.set_xlim(-L/2, L/2)

    fig.tight_layout()

    state = {"psi1": psi1, "psi2": psi2, "t": 0.0}
    n_per_frame = 3

    def update(frame):
        for _ in range(n_per_frame):
            state["psi1"], state["psi2"] = step_with_potential(
                state["psi1"], state["psi2"],
                P11, P12, P21, P22, ph1, ph2
            )
        state["t"] += n_per_frame * dt

        dens = (np.abs(state["psi1"])**2 + np.abs(state["psi2"])**2).real
        line_dens.set_ydata(dens)
        line_re1.set_ydata(state["psi1"].real)
        line_re2.set_ydata(state["psi2"].real)
        txt_t.set_text(f"t = {state['t']:.1f}")

        # calcula T em tempo real: fração da norma à direita da barreira
        norm_r = np.sum(dens[xs > x_bar_r]) * dx
        norm_t = np.sum(dens) * dx
        T_inst = norm_r / (norm_t + 1e-15)
        txt_T.set_text(f"T_inst (x>barreira) = {T_inst:.3f}")

    ani = animation.FuncAnimation(fig, update, frames=280, interval=28, blit=False)

    if args.get("save"):
        fname = "dirac_klein.gif"
        print(f"Salvando '{fname}'...")
        ani.save(fname, writer=animation.PillowWriter(fps=28))
        print(f"Salvo: {fname}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Paradoxo de Klein — Equação de Dirac 1D",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--k0",     type=float, default=1.0,  help="Momento central do pacote")
    p.add_argument("--V0",     type=float, default=5.0,  help="Altura da barreira (em mc²)")
    p.add_argument("--d",      type=float, default=4.0,  help="Largura da barreira")
    p.add_argument("--mass",   type=float, default=1.0,  help="Massa (unidades naturais)")
    p.add_argument("--static", action="store_true",      help="Só figuras estáticas")
    p.add_argument("--save",   action="store_true",      help="Salvar animação como GIF")
    return vars(p.parse_args())


if __name__ == "__main__":
    args = parse_args()

    # Relatório analítico
    m   = args["mass"]
    k0  = args["k0"]
    V0  = args["V0"]
    d   = args["d"]
    E0  = np.sqrt(k0**2 + m**2)

    print(f"\n{'='*58}")
    print("  Paradoxo de Klein — Equação de Dirac 1D")
    print(f"{'='*58}")
    print(f"  Momento k₀         = {k0:.3f} mc")
    print(f"  Energia E₀         = {E0:.3f} mc²")
    print(f"  Altura da barreira = {V0:.3f} mc²")
    print(f"  Largura d          = {d:.3f}")
    print(f"  {'─'*50}")
    print(f"  Limite V₀ = E+mc²  = {E0+m:.3f} mc²  (início do paradoxo)")
    print(f"  Limite V₀ = E+3mc² = {E0+3*m:.3f} mc²  (T → 1 assintótico)")
    print(f"  {'─'*50}")
    T_d = dirac_transmission(k0, V0, d, m)
    T_s = schrodinger_transmission(k0, V0, d, m)
    print(f"  T_Dirac (analítico)       = {T_d:.6f}")
    print(f"  T_Schrödinger (analítico) = {T_s:.6f}")
    regime = ("Tunelamento normal" if V0 < E0 + m
              else "Ressonâncias" if V0 < E0 + 3*m
              else "PARADOXO DE KLEIN (T→1)")
    print(f"  Regime físico             : {regime}")
    print(f"{'='*58}\n")

    plot_static(args)
    if not args["static"]:
        animate_klein(args)
