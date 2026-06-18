# Grafeno 2D — Férmions de Dirac Sem Massa

> **Nível avançado — Simulação 2 de 4**

---

## Por que o grafeno é especial

O grafeno (folha de carbono com espessura de 1 átomo) tem uma estrutura
eletrônica que, perto dos pontos K e K' da zona de Brillouin, é descrita
**exatamente** pela equação de Dirac 2D com massa zero:

```
H = v_F (σₓ kₓ + σᵧ kᵧ)
E = ±v_F |k|           (dispersão linear, sem gap)
```

onde v_F ≈ 10⁶ m/s ≈ c/300 é a "velocidade da luz" do grafeno.
Os elétrons se comportam como fótons — mas com carga elétrica.

---

## Os 5 módulos

| Módulo | Comando | Física |
|--------|---------|--------|
| 1 — Rede + ZB | `--module 1` | Rede hexagonal, sítios A/B, zona de Brillouin, pontos K |
| 2 — Cone de Dirac | `--module 2` | Dispersão linear 3D, comparação com parabólico |
| 3 — Propagação | `--module 3` | Pacote sem dispersão (m=0) vs com dispersão (m≠0) |
| 4 — Landau | `--module 4` | E_n ∝ √(nB), n=0 em E=0, comparação com Schrödinger |
| 5 — Animação | `--module 5` | m=0 vs m=1 em tempo real: Δx(t) |

---

## 1. Rede Hexagonal e Pontos de Dirac

O grafeno tem **2 átomos por célula unitária** (sítios A e B).
Isso gera dois ramos de banda que se tocam em 6 pontos K da zona de Brillouin:

```
K  = (2π/a)(1/3,  1/√3)
K' = (2π/a)(−1/3, −1/√3)   = −K
```

Os 6 pontos se reduzem a 2 inequivalentes (K e K') por simetria de translação.
Nesses pontos, o gap é **zero** — o grafeno é um semimetal.

---

## 2. Cone de Dirac: dispersão linear

Em torno de cada ponto K, expandindo para k pequeno:

```
H_K = v_F (σₓ δkₓ + σᵧ δkᵧ)
```

Autovalores: E = ±v_F |δk|

Comparação com elétron livre:

| Propriedade | Grafeno (m=0) | Elétron livre (m≠0) |
|-------------|---------------|----------------------|
| Dispersão   | E ∝ |k|       | E ∝ k² (baixo k) |
| Massa efetiva | m* = 0     | m* ≠ 0 |
| Velocidade grupo | v_g = v_F (constante) | v_g = k/m (variável) |
| Gap | 0 | 0 (livre) ou 2mc² (Dirac) |
| Dispersão do pacote | Nenhuma | Sim (alargamento) |

---

## 3. Propagação sem dispersão

Um pacote gaussiano em grafeno **não se alarga** com o tempo.
Isso porque todos os modos têm a mesma velocidade de grupo:

```
v_g = dE/dk = d(v_F|k|)/dk = v_F = constante
```

Na simulação, meça Δx(t) para os dois casos:
- **Grafeno**: Δx ≈ constante (pacote preserva forma)
- **Elétron livre**: Δx cresce ~√t (dispersão normal)

---

## 4. Pseudospin e Quiralidade

O "spin" das matrizes σ no hamiltoniano de grafeno **não é spin real** —
é o **pseudospin**, que indica em qual sítio (A ou B) o elétron está.

**Quiralidade travada:** o pseudospin aponta sempre na direção de k.
Consequência fundamental:

> **Retro-espalhamento (θ = 180°) é proibido por simetria!**

O elétron não pode voltar pelo mesmo caminho — o pseudospin teria que
inverter, mas está travado à direção de k. Isso gera:
- Alta mobilidade em grafeno (poucos defeitos causam espalhamento)
- **Klein perfeito para θ=0**: T = 1 para qualquer V₀, qualquer d

### Fórmula de Katsnelson (2006)

```
T(θ) = cos²θ / (1 − sin²θ · sin²(qd))

Para θ=0: T = 1 SEMPRE  (independente de V₀ e d)
```

---

## 5. Níveis de Landau em Grafeno

Em campo magnético B, os níveis de Landau do grafeno são:

```
E_n = ±v_F √(2nħeB),    n = 0, 1, 2, ...
```

Diferenças em relação ao elétron livre (Schrödinger: E_n = ħωc(n+½)):

| Propriedade | Grafeno | Elétron livre |
|-------------|---------|---------------|
| Escala | ∝ √(nB) | ∝ (n+½)B |
| n=0 | E=0 SEMPRE | E=½ħωc |
| Sinal ± | Elétron e buraco | Só positivo |
| Espaçamento | Decresce ~1/√n | Uniforme = ħωc |

**O nível n=0 em E=0** é responsável pelo **Efeito Hall Quântico Anômalo**:
os platôs de Hall ocorrem em ν = ±4(n+½), deslocados de ½ comparado
ao previsto por Schrödinger. Observado por Novoselov et al. e Zhang et al. em 2005.

### Verificação numérica

```
n=0: E=0.000000  |E²−2nB|=0.00e+00   ✓
n=1: E=1.414214  |E²−2nB|=4.44e−16   ✓
n=2: E=2.000000  |E²−2nB|=0.00e+00   ✓
```

---

## Instalação e uso

```bash
pip install numpy matplotlib
```

```bash
# Tudo (rede + cone + propagação + Landau + animação)
python dirac_graphene.py

# Só cone de Dirac e estrutura de banda
python dirac_graphene.py --module 2

# Só níveis de Landau
python dirac_graphene.py --module 4

# Animação: m=0 vs m≠0
python dirac_graphene.py --module 5

# Campo magnético diferente
python dirac_graphene.py --B 3.0

# Salvar animação como GIF
python dirac_graphene.py --save
```

---

## O que observar

### Figura principal (7 painéis)

| Painel | O que mostra |
|--------|-------------|
| Rede hexagonal | Sítios A (azul) e B (laranja), célula unitária, ligações C–C |
| Zona de Brillouin | Hexágono recíproco, 6 pontos K/K' nos vértices |
| Cone de Dirac 3D | Superfície E(kx,ky): dois cones se tocando em K |
| Estrutura de banda | Linear (grafeno) vs parabólico (semicondutor) |
| Propagação | Snapshots em t=0,2,4,6,8,10: m=0 não alarga, m=1 sim |
| Pseudospin | Setas de helicidade em torno de K |
| Klein | T(θ): pico em θ=0 para grafeno, queda para semicondutor |

### Animação (2 painéis)

- Painel superior: **grafeno m=0** — Δx constante no canto superior direito
- Painel inferior: **elétron livre m=1** — Δx cresce visivelmente

---

## Conexão com experimentos

- **2004** — Geim & Novoselov isolam o grafeno (Nobel 2010)
- **2005** — EHQ anômalo medido (ν = ±2, ±6, ±10, ...) — confirma n=0 em E=0
- **2006** — Katsnelson et al. preveem Klein perfeito no grafeno
- **2007** — Young & Kim medem tunelamento de Klein em grafeno experimentalmente
- **2012** — Espectroscopia de ciclotron mede E_n ∝ √(nB) com precisão de meV
- **2020s** — Grafeno torcido (twisted bilayer) mostra supercondutividade, Mott insulator

---

## Referências

- Wallace, P.R. *The band theory of graphite.* Phys. Rev. 71 (1947).
- Novoselov, K.S. et al. *Electric field effect in atomically thin carbon films.* Science 306 (2004).
- Novoselov, K.S. et al. *Two-dimensional gas of massless Dirac fermions in graphene.* Nature 438 (2005).
- Zhang, Y. et al. *Experimental observation of the quantum Hall effect in graphene.* Nature 438 (2005).
- Katsnelson, M.I. et al. *Chiral tunnelling and the Klein paradox in graphene.* Nature Physics 2 (2006).
- Castro Neto, A.H. et al. *The electronic properties of graphene.* Rev. Mod. Phys. 81 (2009).
