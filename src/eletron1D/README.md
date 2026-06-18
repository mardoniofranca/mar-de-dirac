# Simulação da Equação de Dirac 1D — Elétron Livre

> **Simulação 1 de 3**
> Evolução temporal de um pacote de onda gaussiano via equação de Dirac em 1+1 dimensões, usando o método split-operator com FFT.

---

## O que esta simulação faz

Resolve a equação de Dirac para **um único elétron livre** em uma dimensão espacial:

```
(iγ^μ ∂_μ − m)ψ = 0
```

Em 1+1D, o espinor tem **2 componentes** (em vez de 4 como em 3+1D):

```
ψ = (ψ₁)    ψ₁ = componente "grande"  (partícula)
    (ψ₂)    ψ₂ = componente "pequena" (acoplada à antipartícula)
```

A equação de Dirac neste caso é equivalente ao sistema:

```
i ∂ₜ ψ₁ = −i ∂ₓ ψ₂ + m ψ₁
i ∂ₜ ψ₂ = −i ∂ₓ ψ₁ − m ψ₂
```

---

## Fenômenos físicos observados

### 1. Dispersão do pacote de onda
O pacote gaussiano se alarga com o tempo, assim como na equação de Schrödinger, mas com relação de dispersão relativística:

```
ω(k) = √(k² + m²)
```

### 2. Velocidade de grupo
O centro do pacote se move com velocidade de grupo:

```
v_g = dω/dk = k / √(k² + m²) < c = 1
```

Sempre menor que a velocidade da luz — como esperado para uma partícula massiva.

### 3. Zitterbewegung (tremor quântico)
O fenômeno mais característico da equação de Dirac: o centro do pacote **oscila rapidamente** em torno da trajetória média, com frequência:

```
f_Zitt = 2mc² / ħ  →  em unidades naturais: f_Zitt = 2m
```

Isso acontece porque o espinor inicial contém mistura de estados de energia positiva (+ω) e negativa (−ω). A interferência entre eles produz o tremor.

> **Nota histórica:** O Zitterbewegung foi previsto por Schrödinger em 1930. Para o elétron real, a frequência é ~10²¹ Hz — inobservável diretamente, mas análogo ao tremor de elétrons em grafeno foi medido experimentalmente em 2010.

### 4. Componente de antipartícula
A componente ψ₂ é não-nula mesmo para um estado de partícula. Ela representa a "sombra" de antipartícula que existe em qualquer estado relativístico localizado. Quanto mais localizado o pacote (σ pequeno), maior a contribuição de ψ₂ — isso está ligado à criação virtual de pares.

---

## Método numérico: Split-Operator + FFT

### Por que este método?

O hamiltoniano de Dirac no espaço de momentos é uma matriz 2×2 **exatamente diagonalizável**:

```
H(k) = k·σ₁ + m·σ₃

     = ( m    k )
       ( k   -m )
```

com autovalores ±ω(k) = ±√(k²+m²).

O propagador temporal exato por passo `dt` é:

```
exp(-iH dt) = cos(ω dt)·I − i·sin(ω dt)/ω · H
```

### Algoritmo por passo temporal

```
1. ψ(x,t)  →  FFT  →  ψ̃(k,t)        [espaço → momentos]

2. ψ̃(k,t+dt) = P(k,dt) · ψ̃(k,t)    [multiplica pelo propagador exato]

        P = ( cos(ωdt) − im·sinc(ω)dt   −ik·sinc(ω)dt         )
            ( −ik·sinc(ω)dt              cos(ωdt) + im·sinc(ω)dt )

3. ψ̃(k,t+dt)  →  IFFT  →  ψ(x,t+dt)  [momentos → espaço]
```

### Propriedades do método
- **Unitário exato:** a norma ∫|ψ|²dx é conservada a precisão de máquina
- **Sem spurious modes:** ao contrário de diferenças finitas, não gera modos espúrios de lattice
- **Complexidade:** O(N log N) por passo, graças à FFT
- **Estabilidade:** incondicionalmente estável — pode-se usar dt grande

---

## Instalação e uso

### Dependências

```bash
pip install numpy matplotlib
```

### Execução básica

```bash
python dirac_1d.py
```

### Exemplos de cenários interessantes

```bash
# Elétron lento (não-relativístico: k₀ << m)
python dirac_1d.py --k0 1 --mass 5

# Elétron ultrarelativístico (m → 0)
python dirac_1d.py --k0 4 --mass 0

# Pacote muito localizado (Δx pequeno → grande Zitterbewegung)
python dirac_1d.py --sigma 1 --mass 1

# Alta resolução, mais passos
python dirac_1d.py --N 2048 --steps 2000 --dt 0.005

# Salvar como GIF
python dirac_1d.py --save --fps 30
```

### Todos os parâmetros

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `--N`     | 1024   | Pontos na grade espacial |
| `--L`     | 60.0   | Tamanho do domínio (unidades naturais) |
| `--k0`    | 4.0    | Momento central do pacote |
| `--sigma` | 4.0    | Largura gaussiana inicial |
| `--mass`  | 1.0    | Massa da partícula (0 = sem massa) |
| `--dt`    | 0.01   | Passo temporal |
| `--steps` | 800    | Número total de passos |
| `--save`  | False  | Salvar animação como GIF |
| `--fps`   | 30     | Frames por segundo do GIF |

---

## Sistema de unidades

Esta simulação usa **unidades naturais**:

| Constante | Valor | Unidade SI equivalente |
|-----------|-------|------------------------|
| ħ (Planck reduzida) | 1 | 1.055 × 10⁻³⁴ J·s |
| c (velocidade da luz) | 1 | 2.998 × 10⁸ m/s |
| m (massa do elétron) | 1 | 9.109 × 10⁻³¹ kg |

Para converter para SI, multiplique as distâncias por ħ/(mc) = 3.86 × 10⁻¹³ m (comprimento de Compton do elétron) e os tempos por ħ/(mc²) = 1.29 × 10⁻²¹ s.

---

## O que observar na animação

| Painel | O que mostra |
|--------|-------------|
| Superior | |ψ|² = |ψ₁|² + |ψ₂|² — densidade de probabilidade total |
| Inferior (azul) | Re ψ₁ — parte real da componente de partícula |
| Inferior (laranja) | Re ψ₂ — parte real da componente de antipartícula |

Observações esperadas:
- O pico de |ψ|² se move com velocidade v_g e alarga com o tempo
- Re ψ₁ e Re ψ₂ oscilam na frequência de Zitterbewegung
- Para m=0: ψ₂ = ψ₁ (paridade quiral), sem Zitterbewegung
- Para σ pequeno: ψ₂ é maior em relação a ψ₁ (mais antipartícula virtual)

---

## Próximas simulações (nível básico)

- **Simulação 2:** Estrutura completa do espinor — 4 componentes em 3D, projeção de helicidade e quiralidade
- **Simulação 3:** Zitterbewegung em detalhe — análise espectral e comparação com Schrödinger

---

## Referências

- Thaller, B. *The Dirac Equation*. Springer, 1992.
- Bjorken, J. D.; Drell, S. D. *Relativistic Quantum Mechanics*. McGraw-Hill, 1964.
- Gerritsma, R. et al. *Quantum simulation of the Dirac equation.* Nature 463, 68–71 (2010). — Primeira simulação experimental do Zitterbewegung.
- Fillion-Gourdeau, F. et al. *Numerical solution of the time-dependent Dirac equation in coordinate space without fermion-doubling.* Comp. Phys. Comm. 183 (2012).
