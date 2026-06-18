# Precessão de Thomas — Relatividade e Acoplamento Spin-Órbita

> **Nível avançado — Simulação 3 de 4**

---

## O que é a precessão de Thomas

Em 1926, Llewellyn Thomas resolveu um paradoxo experimental: o splitting
spin-órbita medido nos átomos era exatamente **metade** do valor previsto
por cálculos usando a precessão de Larmor. A discrepância era de fator 2.

Thomas descobriu que a causa é puramente cinemática: um vetor transportado
ao longo de uma trajetória curva num referencial inercial **prece** mesmo
sem nenhum torque externo — apenas pela relatividade.

---

## Os três componentes

### 1. Precessão de Larmor (campo magnético efetivo)

No referencial do elétron em órbita, o núcleo está em movimento.
Uma carga em movimento gera campo magnético. O spin do elétron
interage com esse campo:

```
Ω_Larmor = (g/2) · ω_orb / γ
```

Para g=2 (elétron de Dirac) e v→0: Ω_L → ω_orb.

### 2. Precessão de Thomas (cinemática relativística)

A sequência de boosts de Lorentz ao longo de uma curva não é
fechada — gera uma rotação residual chamada precessão de Wigner/Thomas:

```
Ω_Thomas = −(γ−1) · ω_orb    (sinal negativo: oposto ao orbital)
```

Para v→0: Ω_T → 0 (efeito puramente relativístico).
Para v→c: Ω_T → −ω_orb.

### 3. Precessão total (acoplamento spin-órbita)

```
Ω_SO = Ω_Larmor + Ω_Thomas

Para v→0:  Ω_SO → ω_orb + 0 = ω_orb    (SEM Thomas → errado!)
Com Thomas: Ω_SO → ω_orb − 0 = ?

Limite v→0 COM Thomas:  Ω_SO → ½ω_orb   ← O famoso "fator ½ de Thomas"
```

---

## Verificações numéricas

| Grandeza | Valor (v=0.5c) | Esperado |
|----------|----------------|---------|
| γ | 1.1547 | 1/√(1−0.25) |
| Ω_Larmor | 0.8660 ω | ω/γ = 0.8660 ✓ |
| Ω_Thomas | −0.1547 ω | −(γ−1) = −0.1547 ✓ |
| Ω_total | 0.7113 ω | — |
| Splitting 2p **com** Thomas | 0.04528 meV | NIST: ~0.045 meV ✓ |
| Splitting 2p **sem** Thomas | 0.09057 meV | ×2 errado ✓ |
| Razão sem/com Thomas | 2.0000 | 2.0000 ✓ |
| Norma do spin (conservada) | 1.00000000 | 1.0 ✓ |

---

## Método numérico: integração por Rodrigues

O spin evolui segundo:

```
dS/dt = Ω(t) × S
```

onde Ω(t) = Ω_Larmor(t) + Ω_Thomas(t) depende da posição na órbita.

A integração usa a **fórmula de Rodrigues** (rotação exata por passo dt):

```
S(t+dt) = S·cos(|Ω|dt) + (Ω̂×S)·sin(|Ω|dt) + Ω̂(Ω̂·S)(1−cos(|Ω|dt))
```

Isso preserva |S| = 1 exatamente a cada passo (sem acúmulo de erro).

---

## Estrutura fina do átomo de hidrogênio

O hamiltoniano de acoplamento spin-órbita correto é:

```
H_SO = (α²/2) · (Z/r³) · L·S    (com fator Thomas ½ incluído)
```

Sem o fator Thomas seria H_SO = α²·(Z/r³)·L·S → splitting dobrado.

### Energia perturbativa

```
ΔE = (α²/2) · ⟨1/r³⟩ · ⟨L·S⟩

⟨1/r³⟩ = Z³ / (n³ l(l+½)(l+1))    [em u.a.]

⟨L·S⟩ = ½[j(j+1) − l(l+1) − s(s+1)]
```

Splitting 2p₃/₂ − 2p₁/₂ com Thomas = 0.04528 meV ≈ NIST ✓

---

## Instalação e uso

```bash
pip install numpy matplotlib
```

```bash
# Tudo (estático + animação)
python dirac_thomas.py

# Só figura estática
python dirac_thomas.py --module 1

# Só animação
python dirac_thomas.py --module 5

# Velocidade diferente
python dirac_thomas.py --v 0.8    # muito relativístico: δφ/órbita grande
python dirac_thomas.py --v 0.1    # quase clássico

# Salvar GIF
python dirac_thomas.py --save
```

---

## O que observar

### Figura estática (6 painéis)

| Painel | O que mostra |
|--------|-------------|
| Geometria | Spin (setas coloridas) girando diferente do orbital — cada cor = posição |
| Frequências | Ω_L(v), Ω_T(v), Ω_tot(v): Thomas cancela metade de Larmor |
| Spin vs tempo | Numérico (verde sólido) vs analítico (pontilhado): concordância perfeita |
| Spin 3D | Quivers coloridos mostrando spin no espaço ao longo da trajetória |
| Splitting 2p | Com Thomas (correto) vs sem Thomas (×2): setas de ΔE |
| Níveis finos | ΔE em μeV para estados p e d de n=2,3 |

### Animação (2 painéis)

- **Painel esquerdo**: elétron (ponto azul) na órbita, spin (seta vermelha), velocidade (verde)
- **Painel direito**: ângulo do spin (vermelho) vs ângulo orbital (verde tracejado)
- A defasagem crescente é a precessão de Thomas em ação

---

## Por que isso importa

- **Histórico**: Thomas (1926) resolveu a discrepância ×2 que atormentava a teoria atômica
- **Spintônica**: precessão de Thomas é relevante para transporte de spin em dispositivos
- **GPS**: o satélite GPS tem que corrigir pela precessão geodética (análogo gravitacional de Thomas)
- **Giroscópio relativístico**: a sonda Gravity Probe B (2004–2011) mediu precessão geodética e de Lense-Thirring — confirmação experimental do efeito Thomas gravitacional

---

## Referências

- Thomas, L.H. *The motion of the spinning electron.* Nature 117 (1926).
- Thomas, L.H. *The kinematics of an electron with an axis.* Phil. Mag. 3 (1927).
- Jackson, J.D. *Classical Electrodynamics*, 3ª ed. Wiley, 1999. (§11.8)
- Bargmann, V.; Michel, L.; Telegdi, V.L. *Precession of the polarization of particles.* Phys. Rev. Lett. 2 (1959). — equação BMT
