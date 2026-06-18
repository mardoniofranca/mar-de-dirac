# Níveis de Landau Relativísticos — Equação de Dirac

> **Nível intermediário — Simulação 3 de 3**

---

## O que esta simulação faz

Resolve o hamiltoniano de Dirac para um elétron em campo magnético B
uniforme perpendicular ao plano xy, usando o **gauge de Landau**
A = (0, Bx, 0). A separação de variáveis reduz o problema a um
oscilador harmônico relativístico com solução exata.

---

## A física dos níveis de Landau

### Schrödinger (não-relativístico)

```
E_n = ħωc(n + ½),    n = 0, 1, 2, ...
ωc = eB/mc  (frequência ciclotron)
```

Espaçamento **uniforme** — todos os níveis equidistantes.

### Dirac (relativístico) — resultado exato

```
E_n² = (mc²)² + 2n·ħeB·c    →    E_n = mc²√(1 + 2n·B/m²)
```

em unidades naturais (ħ = c = e = m = 1):

```
E_n = m√(1 + 2nB/m²),    n = 0, 1, 2, ...
```

Espaçamento **comprimido** — cresce como √n, não como n.

### Diferença crucial: verificação numérica

```
n  | E_Dirac  | E_cin(D) | E_Schrö | ΔE_Dirac | ΔE_Schrö | razão
---|----------|----------|---------|----------|----------|-------
0  | 1.000000 | 0.000000 | 0.5000  |    —     |    —     |   —   ← n=0!
1  | 1.732051 | 0.732051 | 1.5000  | 0.73205  | 1.00000  | 0.732
2  | 2.236068 | 1.236068 | 2.5000  | 0.50402  | 1.00000  | 0.504
3  | 2.645751 | 1.645751 | 3.5000  | 0.40968  | 1.00000  | 0.410
```

---

## O nível zero de Landau

O nível n=0 é a assinatura mais pura da equação de Dirac:

- **Schrödinger**: n=0 tem energia ½ħωc (ponto zero do oscilador)
- **Dirac**: n=0 tem energia mc² exata — energia cinética **zero**!

O espinor do nível zero é puramente de "grande componente":

```
Ψ₀ = (G₀(x), 0)ᵀ,    G₀(x) = π^{-1/4} l_B^{-1/2} exp(-x²/2l_B²)
```

Apenas uma gaussiana pura, sem nó. Esse nível é responsável pelo
**Efeito Hall Quântico Anômalo** em grafeno — a plataforma de Hall
ν = 0 está deslocada por ½ em relação ao previsto por Schrödinger.

---

## Espinor de Dirac nos níveis de Landau

Para n ≥ 1, o espinor tem duas componentes acopladas:

```
Ψ_n = ( G_n(x)      )    G_n = função do oscilador de nível n
      ( c_n · G_{n-1}(x) )    F_{n-1} = componente "pequena"

c_n = √(2nB) / (E_n + m)    ≈ Zα  (grau de relatividade)
```

- Para B pequeno: c_n → 0, F → 0 (limite não-relativístico)
- Para B grande: c_n → 1, F comparável a G (relativístico)

---

## Degenerescência e Efeito Hall Quântico

Cada nível de Landau acomoda D estados por unidade de área:

```
D = eB / (2πħ) = B / (2π)    [unidades naturais]
```

Isso significa que ao aumentar B, cada nível "cabe mais elétrons",
e os níveis passam a ser preenchidos sequencialmente conforme B cresce
— produzindo os **platôs do Efeito Hall Quântico**.

O fator de preenchimento ν = n_e / D = 2πn_e / B determina
qual nível está na superfície de Fermi. O EHQ relativístico
(grafeno) tem platôs em ν = ±4(n+½) por causa do nível zero.

---

## Instalação e uso

```bash
pip install numpy matplotlib
```

```bash
# Padrão: B=1, nmax=10
python dirac_landau.py

# Campo forte (mais relativístico)
python dirac_landau.py --B 5.0 --nmax 15

# Campo fraco (quase Schrödinger)
python dirac_landau.py --B 0.1 --nmax 8

# Massa diferente
python dirac_landau.py --B 2.0 --mass 2.0
```

### Parâmetros

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `--B`     | 1.0    | Campo magnético (unidades naturais) |
| `--nmax`  | 10     | Nível de Landau máximo |
| `--mass`  | 1.0    | Massa da partícula |

---

## O que observar na figura

| Painel | O que mostra |
|--------|-------------|
| Diagrama de níveis | Dirac (sólido, comprimido) vs Schrödinger (tracejado, uniforme) |
| E(n) | Curva √n (Dirac) vs reta n (Schrödinger) — divergem rapidamente |
| ΔE_n | Espaçamento decresce ~1/√n para Dirac; constante para Schrödinger |
| Densidades | Funções de onda empilhadas n=0,1,2,3 (nós aumentam com n) |
| Espinor G/F | Componentes grande e pequena do espinor — razão mede relatividade |
| Degenerescência | D = B/2π por nível — base do Efeito Hall Quântico |

---

## Conexão com experimentos

### Grafeno (2004–presente)
Em grafeno, os elétrons são férmions de Dirac **sem massa** com
v_F ≈ c/300. Os níveis de Landau seguem:

```
E_n = v_F √(2nħeB)  ∝  √(nB)
```

Medido por espectroscopia de ciclotron e STM com precisão de meV.
A observação de n=0 (e do EHQ anômalo) foi uma das primeiras
confirmações da natureza relativística dos elétrons em grafeno.

### Nível zero em STM
O pico de densidade de estados em E=0 para n=0 é diretamente
observável em imagens de STM de grafeno em campo magnético.

---

## Referências

- Landau, L. D. *Diamagnetism of metals.* Z. Phys. 64 (1930).
- Johnson, M. H.; Lippmann, B. A. *Motion in a constant magnetic field.* Phys. Rev. 76 (1949).
- Novoselov, K. S. et al. *Two-dimensional gas of massless Dirac fermions in graphene.* Nature 438 (2005).
- Zhang, Y. et al. *Experimental observation of the quantum Hall effect in graphene.* Nature 438 (2005).
- Goerbig, M. O. *Electronic properties of graphene in a strong magnetic field.* Rev. Mod. Phys. 83 (2011).
