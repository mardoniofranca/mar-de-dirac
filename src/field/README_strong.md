# Campo Forte e Colapso Relativístico — Equação de Dirac

> **Nível avançado — Simulação 4 de 4  (conclusão da série)**

---

## O problema central

A equação de Dirac para um elétron num potencial coulombiano V(r) = −Z/r
tem uma solução exata elegante — mas que **quebra** quando Zα → 1.

```
E_n = mc² [ 1 + (Zα / (nᵣ + √(κ²−(Zα)²)))² ]^{−½}  −  mc²
```

O expoente da função de onda perto da origem é:

```
γ = √(κ² − (Zα)²)    →    ψ(r→0) ~ r^{γ−1}
```

Quando **Zα → |κ|**: γ → 0 e ψ ~ r^{−1} — a função de onda diverge
na origem e deixa de ser normalizável. **Isso é o colapso relativístico.**

---

## Os três regimes

```
 Regime 1             Regime 2              Regime 3
 Zα << 1             Zα → 1                Zα > 1
 (subcrítico)         (colapso)             (supercrítico)
 ───────────          ──────────            ────────────────
 γ ≈ 1               γ → 0                 γ imaginário
 ψ regular           ψ ~ r^{γ-1} → r^{-1} Estado ressônancia
 E normal            E → −mc²              Criação de pares
 Zα < 137            Z ≈ 137               Z > 137 → 173
```

---

## Módulo 1 — Espectro E(Z)

A energia do estado 1s cai com Z² para Schrödinger (E ~ −Z²α²/2),
mas a equação de Dirac prevê uma queda mais pronunciada que termina
em E = −mc² quando Z → Z_c2 ≈ 173.

Resultados verificados:

| Z | γ(1s) | E_Dirac (mc²) | E_Schrödinger (mc²) |
|---|-------|---------------|----------------------|
| 1  | 0.99997 | −0.0000266 | −0.0000266 |
| 50 | 0.93106 | −0.01348   | −0.01328   |
| 92 | 0.74113 | −0.25887   | −0.22502   |
| 130 | 0.31631 | −0.70...  | −0.44956   |

Diferenças entre Dirac e Schrödinger tornam-se massivas para Z > 50.

---

## Módulo 2 — Expoente γ(Z)

```
γ = √(κ² − (Zα)²)    para κ = −1:    γ = √(1 − (Zα)²)
```

| γ | Z | Comportamento de ψ |
|---|---|----------------------|
| 1.0 | 0 | ψ ~ r⁰ = 1 (regular) |
| 0.5 | ~118 | ψ ~ r^{−0.5} (limite normalizável) |
| 0.0 | ~137 | ψ ~ r^{−1} (COLAPSO — não normalizável) |

Para γ < ½, a função de onda ainda é normalizável em 3D mas o
potencial coulombiano deixa de ser autoadjunto — fisicamente,
o elétron "cai na origem".

---

## Módulo 3 — Funções de onda sob colapso

Para Z grande, o pico da densidade eletrônica move-se para a origem
como rₘₐₓ ~ a₀/Z. Em coordenadas escaladas por Z (r·Z), o pico
*permanece* em posição fixa — mas a **altura** cresce e a **largura** estreita,
indicando concentração cada vez maior de probabilidade na região nuclear.

---

## Módulo 4 — Componente "pequena" F(r)

O espinor de Dirac tem componente grande G e componente pequena F.
Para Z→0: |F|/|G| → Zα (muito pequeno). Para Z→137: |F|/|G| → 1.

A componente "pequena" deixa de ser pequena — o elétron torna-se
completamente relativístico mesmo estando ligado.

Verificação:
```
Z=1:   |F|/|G| ≈ 0.000007 ≈ α           ✓
Z=50:  |F|/|G| ≈ 0.84  (já significativo!)
Z=100: |F|/|G| ≈ 0.61
Z=130: |F|/|G| ≈ 0.45
```

---

## Módulo 5 — Região supercrítica e criação de pares

Para Z > Z_c1 ≈ 137: a solução de Dirac com γ real deixa de existir.
O estado 1s mergulha no contínuo de energia negativa (mar de Dirac).

Para Z > Z_c2 ≈ 173: a energia do estado atingiria −mc², o que significa
que o vácuo torna-se **instável**: cria-se espontaneamente um par e⁺e⁻
— o elétron é "capturado" pelo campo coulombiano, o pósitron é emitido.

Taxa análoga ao efeito Schwinger:
```
Γ ∝ (Z−Z_c2)² · exp(−A/(Z−Z_c2))
```

Este fenômeno ainda não foi observado: os elementos mais pesados
conhecidos chegam apenas a Z=118 (Oganessônio). Para atingir Z≈173
seria necessário colidir dois núcleos pesados (U+U, Z=92+92=184).
O experimento GSI/FAIR na Alemanha busca esse sinal.

---

## Comparação com Klein-Gordon

A equação de Klein-Gordon (spin 0) também tem colapso:
```
E_KG = mc²√(1−(Zα/n)²) − mc²
```

O colapso ocorre em Z=n/α ≈ 137n — o mesmo limiar de Dirac para
o estado fundamental! Mas o espinor de Dirac (spin ½) tem estrutura
mais rica: γ = √(κ²−(Zα)²) vs γ_KG = √(n²−(Zα)²).

---

## Instalação e uso

```bash
pip install numpy matplotlib
```

```bash
# Figura principal (6 painéis)
python dirac_strong.py

# Só espectro + colapso
python dirac_strong.py --module 1

# Só tabela comparativa Dirac/KG/Schrödinger
python dirac_strong.py --module 2

# Estende a região supercrítica até Z=250
python dirac_strong.py --Zmax 250
```

---

## O que observar na figura (6 painéis)

| Painel | O que mostra |
|--------|-------------|
| E(Z) | Mergulho dos níveis — Dirac (sólido) muito mais que Schrödinger (pont.) |
| γ(Z) | Queda de 1 para 0 entre Z=1 e Z=137 — brusca perto do limite |
| ψ(rZ) | Pico da função de onda concentrando na origem com Z |
| \|F\|/\|G\| | Razão cresce de Zα (linear) até 1 — relatividade total |
| Supercrítico | Diagrama de fases dos três regimes com taxa de pares |
| Taxa Γ(Z) | Taxa de criação espontânea para Z > 173 (experimento futuro) |

---

## Referências

- Dirac, P.A.M. *The quantum theory of the electron.* Proc. R. Soc. A 117 (1928).
- Case, K.M. *Singular potentials.* Phys. Rev. 80 (1950). — análise do colapso
- Greiner, W. et al. *Quantum Electrodynamics of Strong Fields.* Springer, 1985.
- Zeldovich, Ya.B.; Popov, V.S. *Electronic structure of superheavy atoms.* Sov. Phys. Usp. 14 (1972).
- Reinhardt, J.; Greiner, W. *Quantum electrodynamics of strong fields.* Rep. Prog. Phys. 40 (1977).
- Schweppe, J. et al. *Observation of a peak structure in positron spectra from U+Cm.* Phys. Rev. Lett. 51 (1983). — primeiros indícios experimentais
