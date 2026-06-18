# Simulações da Equação de Dirac — 1 Elétron

Série completa de simulações computacionais da equação de Dirac aplicada a um único elétron, organizada em três níveis de complexidade física.

```
src/
├── eletron1D/    → Pacote de onda livre e Zitterbewegung
├── hydrogen/     → Átomo de hidrogênio relativístico
├── klein/        → Paradoxo de Klein
├── landau/       → Níveis de Landau relativísticos
├── antimat/      → Mar de Dirac e antimatéria
├── graphene/     → Grafeno 2D — férmions sem massa
├── prec/         → Precessão de Thomas
└── field/        → Campo forte e colapso relativístico
```

> **Sistema de unidades:** todas as simulações usam unidades naturais `ħ = c = m_e = 1` salvo indicação contrária.  
> **Dependências:** `numpy`, `matplotlib` — instale com `pip install numpy matplotlib`.

---

## Equação de Dirac

A equação fundamental de toda a série:

```
(iγ^μ ∂_μ − m) ψ = 0
```

Em 1+1D, o espinor ψ tem 2 componentes `(ψ₁, ψ₂)`. Em 3+1D, 4 componentes. O hamiltoniano de Dirac livre é:

```
H = α·p + βm,    α = (0  σ),    β = (I   0 )
                     (σ  0)         (0  −I)
```

---

## Nível Básico

### `eletron1D/` — Pacote de Onda Livre e Zitterbewegung

**Física:** Evolução temporal do espinor de 2 componentes para um elétron livre em 1+1D. Primeiro contato com a estrutura do espinor e o Zitterbewegung — o tremor quântico de frequência `ω = 2mc²/ħ` que não existe na equação de Schrödinger.

**Método numérico:** Split-operator com FFT. O propagador exato no espaço de momentos é:

```
exp(−iH(k)dt) = cos(ωdt)·I − i·sin(ωdt)/ω · H(k)
```

Custo `O(N log N)` por passo, unitário exato (norma conservada a precisão de máquina).

**O que observar:**
- Curva azul `|ψ|²`: pico se move com velocidade de grupo `v_g = k/√(k²+m²) < c`
- Curva amarela `Re ψ₁`: oscila na frequência do Zitterbewegung
- Para `m=0`: Zitterbewegung desaparece, `v_g = c`
- Para `σ` pequeno: pacote mais localizado dispersa mais rápido (Heisenberg)

```bash
cd src/eletron1D
python dirac_1d.py                    # parâmetros padrão
python dirac_1d.py --k0 6 --mass 0   # elétron ultrarelativístico
python dirac_1d.py --sigma 1          # pacote muito localizado
python dirac_1d.py --save             # salva animação como GIF
```

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `--k0`    | 4.0    | Momento central |
| `--sigma` | 4.0    | Largura gaussiana |
| `--mass`  | 1.0    | Massa (0 = sem massa) |
| `--dt`    | 0.01   | Passo temporal |
| `--N`     | 1024   | Pontos na grade |

---

## Nível Intermediário

### `hydrogen/` — Átomo de Hidrogênio Relativístico

**Física:** Hamiltoniano de Dirac com potencial coulombiano `V(r) = −Z/r` em coordenadas esféricas. A separação angular reduz o problema a duas equações radiais acopladas para `G(r)` (componente grande) e `F(r)` (componente pequena).

**Energia exata de Dirac:**

```
E_n = mc² [ 1 + (Zα / (nᵣ + √(κ²−(Zα)²)))² ]^{−½}  −  mc²

nᵣ = n − |κ|,    γ = √(κ²−(Zα)²)
```

**Números quânticos:**

| Estado | κ | l | j | Splitting |
|--------|---|---|---|-----------|
| ns₁/₂ | −1 | 0 | ½ | — |
| np₁/₂ | +1 | 1 | ½ | degenerado com ns₁/₂ (Schrö.) |
| np₃/₂ | −2 | 1 | 3/2 | separado por ~0.045 meV (n=2) |
| nd₃/₂ | +2 | 2 | 3/2 | — |

**Método numérico:** RK4 com grade logarítmica + truncamento adaptativo quando a solução começa a divergir (solução irregular).

**O que observar:**
- Tabela no terminal: `E_Dirac` vs `E_Schrödinger` e correção `ΔE` em eV
- Painel de splitting n=2: 3 barras lado a lado (Dirac / Sommerfeld / Schrödinger)
- Componente pequena `F(r)`: razão `max|F|/max|G| ≈ Zα` mede o grau de relatividade

```bash
cd src/hydrogen
python dirac_hydrogen.py              # H, n=1,2,3
python dirac_hydrogen.py --nmax 4     # adiciona n=4
python dirac_hydrogen.py --Z 6        # carbono hidrogenóide C⁵⁺
python dirac_hydrogen.py --Z 92       # urânio U⁹¹⁺ — limite relativístico
```

---

### `klein/` — Paradoxo de Klein

**Física:** Elétron relativístico incidindo sobre barreira de potencial retangular `V₀`. Três regimes:

```
V₀ < E+mc²          → tunelamento exponencial   (como Schrödinger)
E+mc² < V₀ < E+3mc² → ressonâncias de Fabry-Pérot relativísticas
V₀ > E+3mc²          → PARADOXO DE KLEIN: T > 0, mesmo T → 1
```

**Transmissão analítica de Dirac:**

```
T = |cos(qd) − (i/2)(αq/αk + αk/αq)sin(qd)|^{−2}

q = √[(E−V₀+m)(E−V₀−m)]    (real para V₀ > E+m!)
```

Para `V₀ > E+m`: `q` é **real** → ondas propagantes dentro da barreira → transmissão não-nula. Interpretação: criação virtual de pares na parede da barreira.

**O que observar:**
- Curva `T(V₀)`: Schrödinger → 0, Dirac sobe de volta no regime Klein
- Painel inferior: `|ψ|²` nos 3 regimes — no regime 3, o pacote emerge intacto
- No terminal: `T_Dirac` vs `T_Schrödinger` para a barreira selecionada

```bash
cd src/klein
python dirac_klein.py                 # V₀=5mc² (regime Klein)
python dirac_klein.py --V0 3.0        # regime de ressonâncias
python dirac_klein.py --V0 0.5        # tunelamento normal
python dirac_klein.py --static        # só figuras estáticas
python dirac_klein.py --save          # salva animação GIF
```

---

### `landau/` — Níveis de Landau Relativísticos

**Física:** Elétron em campo magnético `B` uniforme (gauge de Landau `A=(0,Bx,0)`). O hamiltoniano de Dirac reduz a um oscilador harmônico com autovalores exatos.

**Comparação fundamental:**

```
Schrödinger:  E_n = ħωc(n + ½)       → espaçamento UNIFORME
Dirac:        E_n = mc²√(1 + 2nB/m²) → espaçamento ∝ 1/√n
```

**Nível zero:** `n=0` tem `E=0` exato em Dirac (energia cinética zero). Schrödinger prevê `E=½ħωc`. Essa diferença de `½ħωc` produz o **Efeito Hall Quântico Anômalo** observado em grafeno.

**Espinor de Dirac nos níveis de Landau:**

```
Ψ_n = ( G_n(x)          )    n ≥ 1
      ( c_n · G_{n−1}(x) )

c_n = √(2nB) / (E_n + m),    Ψ_0 = (G_0, 0)ᵀ   [só componente grande]
```

```bash
cd src/landau
python dirac_landau.py                # B=1, n=0..10
python dirac_landau.py --B 5.0        # campo forte → mais relativístico
python dirac_landau.py --B 0.1        # campo fraco → quase Schrödinger
python dirac_landau.py --nmax 15      # mais níveis
```

---

## Nível Avançado

### `antimat/` — Mar de Dirac e Antimatéria

**Física:** Quatro módulos progressivos sobre a estrutura do vácuo relativístico.

**Módulo 1 — Espectro:** dois ramos `E = ±√(k²+m²)`, gap de `2mc²`, mar de Dirac preenchido como vácuo. Buraco no mar = pósitron.

**Módulo 2 — Propagador de Feynman-Stückelberg:** o pósitron é matematicamente equivalente a um elétron propagando para trás no tempo. Base de toda a QED perturbativa.

**Módulo 3 — Efeito Schwinger:** campo elétrico suficientemente forte cria pares do vácuo:

```
Γ ∝ (eE)² exp(−π m²c³ / eEħ)
E_crit = m²c³/eħ ≈ 1.32 × 10¹⁸ V/m
```

**Módulo 4 — Animação:** pacote elétron + pósitron se aproximam e aniquilam. Projeção espectral `ρ±(x)` separa os dois ramos em tempo real.

```bash
cd src/antimat
python dirac_antimatter.py            # tudo
python dirac_antimatter.py --module 1 # só espectro
python dirac_antimatter.py --module 3 # só Schwinger
python dirac_antimatter.py --module 4 # só animação
python dirac_antimatter.py --save     # salva GIF
```

---

### `graphene/` — Grafeno 2D — Férmions de Dirac Sem Massa

**Física:** Em grafeno, os elétrons perto dos pontos K da zona de Brillouin obedecem à equação de Dirac 2D com `m=0`:

```
H = v_F (σₓkₓ + σᵧkᵧ),    E = ±v_F|k|,    v_F ≈ c/300
```

**Cinco módulos:**

| Módulo | Conteúdo |
|--------|----------|
| Rede + ZB | Rede hexagonal, sítios A/B, 6 pontos K/K' |
| Cone de Dirac | Dispersão linear 3D vs parabólica de semicondutor |
| Propagação | Pacote sem dispersão (m=0) vs alargamento (m≠0) |
| Pseudospin | Helicidade travada → retro-espalhamento proibido |
| Klein perfeito | `T(θ=0) = 1` para qualquer `V₀` e `d` |
| Landau grafeno | `E_n = v_F√(2nB)`, `n=0` em `E=0` (EHQ anômalo) |

**Resultado chave de Katsnelson (2006):**

```
T(θ) = cos²θ / (1 − sin²θ·sin²(qd))

θ = 0 (incidência normal): T = 1 SEMPRE
```

```bash
cd src/graphene
python dirac_graphene.py              # tudo
python dirac_graphene.py --module 2   # cone de Dirac 3D
python dirac_graphene.py --module 4   # Landau grafeno
python dirac_graphene.py --module 5   # animação m=0 vs m=1
python dirac_graphene.py --B 3.0      # campo magnético diferente
python dirac_graphene.py --save       # salva GIF
```

---

### `prec/` — Precessão de Thomas

**Física:** Efeito cinemático puro da relatividade: um vetor (o spin) transportado ao longo de uma trajetória curva prece mesmo sem torque externo.

**Os três componentes:**

```
Ω_Larmor = (g/2) · ω_orb / γ          (campo B efetivo no ref. do elétron)
Ω_Thomas  = −(γ−1) · ω_orb             (cinemática: sinal negativo!)
Ω_total   = Ω_L + Ω_T

Limite v→0:  Ω_total → ½ω_orb          ← O FAMOSO FATOR ½ DE THOMAS
```

Sem Thomas, o splitting spin-órbita seria o dobro do observado. Thomas (1926) resolveu isso e salvou a teoria atômica.

**Método:** integração via fórmula de Rodrigues (rotação exata por passo `dt`, norma do spin conservada a `1.00000000` a cada passo).

**Verificação:**
```
Splitting 2p sem Thomas = 0.09057 meV  (×2 errado)
Splitting 2p com Thomas = 0.04528 meV  ≈ NIST ~0.045 meV  ✓
```

```bash
cd src/prec
python dirac_thomas.py                # tudo
python dirac_thomas.py --v 0.8        # muito relativístico
python dirac_thomas.py --v 0.1        # quase clássico
python dirac_thomas.py --module 1     # só figura estática
python dirac_thomas.py --module 5     # só animação
python dirac_thomas.py --save         # salva GIF
```

---

### `field/` — Campo Forte e Colapso Relativístico

**Física:** O que acontece quando `Zα → 1`? A solução de Dirac para o átomo hidrogenóide quebra em `Z ≈ 137`.

**Expoente da função de onda:**

```
γ = √(κ² − (Zα)²)    →    ψ(r→0) ~ r^{γ−1}

γ = 1.0:  Z = 0    (livre)
γ = 0.5:  Z ≈ 118  (limite normalizável: ψ ~ r^{−½})
γ = 0.0:  Z ≈ 137  (COLAPSO: ψ ~ r^{−1}, não normalizável)
```

**Três regimes:**

```
Subcrítico  (Z < 137):  solução regular, E normal
Colapso     (Z ≈ 137):  γ → 0, ψ diverge na origem
Supercrítico (Z > 137): estado ressônancia no contínuo negativo
Criação de pares (Z > 173): E₁ₛ < −mc², vácuo instável → e⁺e⁻ espontâneo
```

**Verificações:**
```
γ(Z=1)  = 0.99997337 = √(1−α²)  ✓
Z_crit  = 137.0360   = 1/α       ✓
Z=92(U): E₁ₛ = −0.259 mc²  (elétron 1s usa 26% da energia de repouso!)
```

```bash
cd src/field
python dirac_strong.py                # tudo
python dirac_strong.py --module 1     # espectro + colapso
python dirac_strong.py --module 2     # tabela Dirac vs KG vs Schrödinger
python dirac_strong.py --Zmax 250     # estende região supercrítica
```

---

## Mapa de progressão

```
eletron1D          → pacote de onda, spinor, Zitterbewegung
    │
    ├── hydrogen   → potencial coulombiano, espectro, estrutura fina
    │
    ├── klein      → barreira de potencial, paradoxo de transmissão
    │
    ├── landau     → campo magnético, oscilador harmônico relativístico
    │
    ├── antimat    → vácuo como mar preenchido, pósitron, criação de pares
    │
    ├── graphene   → m=0, cone de Dirac, pseudospin, Klein perfeito
    │
    ├── prec       → cinemática relativística, acoplamento spin-órbita
    │
    └── field      → campo extremo, colapso da solução, QED de campos fortes
```

O fio condutor é o expoente `γ = √(κ²−(Zα)²)`:
- aparece pela primeira vez em `hydrogen/` como correção pequena
- cresce em importância em `field/` até provocar o colapso total da teoria

---

## Referências Principais

| Simulação | Referências |
|-----------|------------|
| eletron1D | Schrödinger (1930) — Zitterbewegung; Thaller (1992) *The Dirac Equation* |
| hydrogen  | Greiner (1990) *Relativistic Quantum Mechanics*; NIST ASD |
| klein     | Klein (1929) Z. Phys. 53; Katsnelson et al. Nature Physics 2 (2006) |
| landau    | Landau (1930); Novoselov et al. Nature 438 (2005) |
| antimat   | Dirac (1930); Anderson (1933) — Nobel; Schwinger (1951) Phys. Rev. 82 |
| graphene  | Wallace (1947); Novoselov & Geim (2004) — Nobel 2010 |
| prec      | Thomas (1926) Nature 117; Jackson *Classical ED* §11.8 |
| field     | Greiner et al. (1985) *QED of Strong Fields*; Zeldovich & Popov (1972) |
