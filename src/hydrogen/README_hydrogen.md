# Simulação de Dirac — Átomo de Hidrogênio Relativístico

> **Nível intermediário — Simulação 1 de 3**

---

## O que esta simulação faz

Resolve a equação de Dirac com potencial coulombiano V(r) = −Z/r em
coordenadas esféricas (unidades atômicas: ħ = m_e = e = 1, c ≈ 137).

Após separação angular, o espinor se reduz a duas funções radiais acopladas:

```
dG/dr = −(κ/r) G  +  [(E + c² − V)/c] F
dF/dr =  (κ/r) F  −  [(E − c² − V)/c] G
```

onde G(r) é a componente grande e F(r) é a componente pequena do espinor.

---

## Números quânticos de Dirac

| Símbolo | Nome | Valores |
|---------|------|---------|
| n | principal | 1, 2, 3, … |
| κ | relativístico | ±1, ±2, … (κ ≠ 0) |
| l | orbital | ∣κ∣−1 (se κ<0) ou κ (se κ>0) |
| j | total | ∣κ∣ − 1/2 |

Estados para n=2:

| Estado | κ | l | j |
|--------|---|---|---|
| 2s₁/₂ | −1 | 0 | 1/2 |
| 2p₁/₂ | +1 | 1 | 1/2 |
| 2p₃/₂ | −2 | 1 | 3/2 |

---

## Energia exata de Dirac

```
E = m c² { [1 + (Zα / (nᵣ + √(κ²−(Zα)²)))²]^(−1/2) } − m c²
```

onde nᵣ = n − |κ| é o número quântico radial e α ≈ 1/137.

### Resultados para H (Z=1)

| Estado | E_Dirac (u.a.) | E_Schrödinger (u.a.) | ΔE (eV) |
|--------|---------------|----------------------|---------|
| 1s₁/₂ | −0.500006657 | −0.500000000 | −0.000181 |
| 2s₁/₂ | −0.125002080 | −0.125000000 | −0.000057 |
| 2p₁/₂ | −0.125002080 | −0.125000000 | −0.000057 |
| 2p₃/₂ | −0.125000416 | −0.125000000 | −0.000011 |

O **splitting spin-órbita** entre 2p₁/₂ e 2p₃/₂ é ~0.045 meV —
medido experimentalmente via espectroscopia de alta resolução.

(O deslocamento de Lamb adicional ~1058 MHz requer QED, além de Dirac.)

---

## Método numérico

### Estratégia: RK4 com grade logarítmica + truncamento

Uma grade logarítmica `r = geomspace(r_min, r_max)` concentra pontos
perto da origem, onde as funções variam rapidamente (singularidade
coulombiana), e se espaça para longe, onde as funções decaem suavemente.

**Condição inicial na origem** — solução regular r^γ:
```
G(r₀) = A · r₀^γ
F(r₀) = B · r₀^γ,   γ = √(κ² − (Zα)²)
```

**Truncamento adaptativo** — a solução fisicamente correta decai
exponencialmente para r → ∞. Quando a densidade G² + F² começa a
crescer (solução irregular), a integração é interrompida e os valores
são zerados. Isso evita overflow sem perder a forma da função de onda
na região fisicamente relevante.

### Por que não usar shooting method?

O shooting method (ajustar E até a função de onda satisfazer a condição
de contorno em r → ∞) seria mais rigoroso para estados excitados, mas
a fórmula exata de Dirac fornece E diretamente, tornando o shooting
desnecessário — usamos a energia analítica e integramos apenas para
obter a forma espacial da função de onda.

---

## Instalação e uso

```bash
pip install numpy matplotlib
```

```bash
# Hidrogênio, estados até n=3
python dirac_hydrogen.py

# Até n=4
python dirac_hydrogen.py --nmax 4

# Carbono hidrogenóide (C⁵⁺, Z=6) — efeitos relativísticos muito maiores
python dirac_hydrogen.py --Z 6

# Urânio hidrogenóide (Z=92) — limite relativístico extremo
python dirac_hydrogen.py --Z 92 --nmax 2
```

---

## O que observar na figura

| Painel | O que mostra |
|--------|-------------|
| Espectro (esq. cima) | Níveis de Dirac (-) vs Schrödinger (×); linha = correção |
| Splitting n=2 (dir. cima) | 3 teorias comparadas: Dirac / Sommerfeld 1ª ordem / Schrödinger |
| Densidade 1s (esq. baixo) | G²+F² vs só G² — contribuição da componente pequena |
| Estados n=2 (centro baixo) | As 3 curvas têm formas diferentes mas energias quase iguais |
| Componente F (dir. baixo) | max\|F\|/max\|G\| ≈ Zα — mede o grau de relatividade |

Para Z pesado (ex. Z=92), max|F|/max|G| ≈ 92/137 ≈ 0.67 — a componente
pequena não é mais pequena! O tratamento relativístico torna-se essencial.

---

## Referências

- Greiner, W. *Relativistic Quantum Mechanics*. Springer, 1990. (Cap. 9)
- Dirac, P. A. M. *The quantum theory of the electron.* Proc. R. Soc. A 117 (1928).
- Johnson, W. R. *Atomic Structure Theory.* Springer, 2007.
