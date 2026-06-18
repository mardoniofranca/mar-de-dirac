"""
Simulação da Equação de Dirac 1D — Elétron Livre
=================================================
Resolve a equação de Dirac em 1+1 dimensão para um único elétron livre,
usando o método split-operator com FFT.

Sistema de unidades naturais: ħ = c = 1

Uso:
    python dirac_1d.py                   # parâmetros padrão
    python dirac_1d.py --k0 6 --mass 0  # elétron ultrarelativístico
    python dirac_1d.py --save            # gera dirac_1d.gif
    python dirac_1d.py --help
"""

import numpy as np
import argparse
import sys
import platform

import matplotlib
# Seleciona backend antes de importar pyplot
_os = platform.system()
if _os == "Darwin":
    matplotlib.use("MacOSX")
elif _os == "Linux":
    matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ---------------------------------------------------------------------------
# Parâmetros padrão
# ---------------------------------------------------------------------------
DEFAULTS = dict(
    N     = 1024,
    L     = 60.0,
    k0    = 4.0,
    sigma = 4.0,
    mass  = 1.0,
    dt    = 0.01,
    steps = 800,
    save  = False,
    fps   = 30,
)


# ---------------------------------------------------------------------------
# Física
# ---------------------------------------------------------------------------

def init_spinor(xs, dx, k0, sigma, mass):
    """Pacote gaussiano de energia positiva. Retorna (psi1, psi2) normalizados."""
    envelope = np.exp(-0.25 * (xs / sigma) ** 2) / (np.pi * sigma ** 2) ** 0.25
    psi1 = envelope * np.exp(1j * k0 * xs)
    omega0 = np.sqrt(k0 ** 2 + mass ** 2)
    ratio = k0 / (omega0 + mass) if (omega0 + mass) > 1e-12 else 0.0
    psi2 = ratio * psi1
    norm = np.sqrt(np.sum((np.abs(psi1) ** 2 + np.abs(psi2) ** 2) * dx))
    return psi1 / norm, psi2 / norm


def make_propagators(ks, mass, dt):
    """Propagador exato no espaço de momentos: exp(-i H(k) dt)."""
    omega = np.sqrt(ks ** 2 + mass ** 2)
    cos_w  = np.cos(omega * dt)
    sinc_w = np.where(omega > 1e-12, np.sin(omega * dt) / omega, dt)
    P11 =  cos_w - 1j * mass * sinc_w
    P12 = -1j * ks * sinc_w
    P21 = -1j * ks * sinc_w
    P22 =  cos_w + 1j * mass * sinc_w
    return P11, P12, P21, P22


def step(psi1, psi2, P11, P12, P21, P22):
    """Um passo temporal via FFT → propaga → IFFT."""
    phi1 = np.fft.fft(psi1)
    phi2 = np.fft.fft(psi2)
    return (np.fft.ifft(P11 * phi1 + P12 * phi2),
            np.fft.ifft(P21 * phi1 + P22 * phi2))


def observables(psi1, psi2, xs, dx):
    """Retorna (density, x_mean, x_std)."""
    dens = (np.abs(psi1) ** 2 + np.abs(psi2) ** 2).real
    norm = np.sum(dens) * dx
    xm   = np.sum(xs * dens) * dx / norm
    x2m  = np.sum(xs ** 2 * dens) * dx / norm
    return dens / norm, xm, np.sqrt(max(x2m - xm ** 2, 0.0))


# ---------------------------------------------------------------------------
# Simulação principal
# ---------------------------------------------------------------------------

def run(params):
    N     = params["N"]
    L     = params["L"]
    k0    = params["k0"]
    sigma = params["sigma"]
    mass  = params["mass"]
    dt    = params["dt"]

    xs = np.linspace(-L / 2, L / 2, N, endpoint=False)
    dx = xs[1] - xs[0]
    ks = np.fft.fftfreq(N, d=dx) * 2 * np.pi

    psi1, psi2 = init_spinor(xs, dx, k0, sigma, mass)
    P11, P12, P21, P22 = make_propagators(ks, mass, dt)

    omega0  = np.sqrt(k0 ** 2 + mass ** 2)
    v_group = k0 / omega0 if omega0 > 1e-12 else 0.0
    f_zitt  = 2 * omega0

    print(f"\n{'='*55}")
    print("  Simulação Dirac 1D — Elétron Livre")
    print(f"{'='*55}")
    print(f"  Momento k₀         = {k0:.2f}")
    print(f"  Massa m            = {mass:.2f}")
    print(f"  Largura inicial σ  = {sigma:.2f}")
    print(f"  Velocidade grupo   = {v_group:.4f} c")
    print(f"  Freq. Zitterbew.   = {f_zitt:.4f} (×ħ/m c²)")
    print(f"  Domínio            = [{-L/2:.0f}, {L/2:.0f}]")
    print(f"  Pontos de grade    = {N}")
    print(f"  Passo temporal dt  = {dt}")
    print(f"{'='*55}\n")

    # --- figura ---
    fig, (ax_dens, ax_spin) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    fig.suptitle("Equação de Dirac 1D — Pacote de Onda (Elétron Livre)", fontsize=13)

    dens0, xm0, xs0 = observables(psi1, psi2, xs, dx)

    line_dens, = ax_dens.plot(xs, dens0, color="#378ADD", lw=2, label=r"$|\psi|^2$")
    ax_dens.set_ylabel(r"$|\psi_1|^2 + |\psi_2|^2$", fontsize=11)
    ax_dens.set_ylim(0, dens0.max() * 2.5)
    ax_dens.legend(loc="upper right", fontsize=9)
    ax_dens.set_title(
        f"k₀={k0}, m={mass}, σ={sigma}  |  v_g={v_group:.3f} c  |  f_Zitt={f_zitt:.3f}",
        fontsize=10,
    )
    txt_info = ax_dens.text(0.02, 0.90, "", transform=ax_dens.transAxes,
                            fontsize=9, color="#378ADD")

    line_re1, = ax_spin.plot(xs, psi1.real, color="#378ADD", lw=1.5,
                             alpha=0.85, label=r"Re $\psi_1$ (partícula)")
    line_re2, = ax_spin.plot(xs, psi2.real, color="#EF9F27", lw=1.5,
                             alpha=0.85, label=r"Re $\psi_2$ (anti-part.)")
    ax_spin.axhline(0, color="gray", lw=0.5, ls="--")
    ax_spin.set_ylabel("Amplitude", fontsize=11)
    ax_spin.set_xlabel("x (unidades naturais)", fontsize=11)
    ax_spin.legend(loc="upper right", fontsize=9)
    lim = max(np.abs(psi1).max(), np.abs(psi2).max()) * 1.5
    ax_spin.set_ylim(-lim, lim)

    txt_time = fig.text(0.01, 0.01, "t = 0.00", fontsize=9, color="gray")
    fig.tight_layout()

    # estado mutável sem closure sobre variáveis locais
    state = {"psi1": psi1, "psi2": psi2, "t": 0.0}
    n_per_frame = 4   # passos por frame

    def update(frame):
        for _ in range(n_per_frame):
            state["psi1"], state["psi2"] = step(
                state["psi1"], state["psi2"], P11, P12, P21, P22
            )
        state["t"] += n_per_frame * dt

        dens, xm, xstd = observables(state["psi1"], state["psi2"], xs, dx)

        line_dens.set_ydata(dens)
        line_re1.set_ydata(state["psi1"].real)
        line_re2.set_ydata(state["psi2"].real)
        txt_time.set_text(f"t = {state['t']:.2f}")
        txt_info.set_text(f"⟨x⟩ = {xm:.2f}   Δx = {xstd:.2f}")

    # blit=False evita o bug '_get_view' no matplotlib antigo
    ani = animation.FuncAnimation(
        fig, update, frames=200, interval=30, blit=False
    )

    if params.get("save"):
        fname = "dirac_1d.gif"
        print(f"Salvando '{fname}'...")
        ani.save(fname, writer=animation.PillowWriter(fps=params["fps"]))
        print(f"Salvo: {fname}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Simulação da Equação de Dirac 1D — Elétron Livre",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--N",     type=int,   default=DEFAULTS["N"],     help="Pontos na grade")
    p.add_argument("--L",     type=float, default=DEFAULTS["L"],     help="Tamanho do domínio")
    p.add_argument("--k0",    type=float, default=DEFAULTS["k0"],    help="Momento central k₀")
    p.add_argument("--sigma", type=float, default=DEFAULTS["sigma"], help="Largura gaussiana σ")
    p.add_argument("--mass",  type=float, default=DEFAULTS["mass"],  help="Massa (0 = sem massa)")
    p.add_argument("--dt",    type=float, default=DEFAULTS["dt"],    help="Passo temporal")
    p.add_argument("--steps", type=int,   default=DEFAULTS["steps"], help="Passos totais")
    p.add_argument("--save",  action="store_true",                   help="Salvar como GIF")
    p.add_argument("--fps",   type=int,   default=DEFAULTS["fps"],   help="FPS do GIF")
    return vars(p.parse_args())


if __name__ == "__main__":
    run(parse_args())
