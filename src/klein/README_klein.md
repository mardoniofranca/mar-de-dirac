# Paradoxo de Klein — Equação de Dirac 1D

> **Nível intermediário — Simulação 2 de 3**

---

## O que é o paradoxo de Klein

Em 1929, Oskar Klein descobriu algo perturbador: quando um elétron
relativístico incide sobre uma barreira de potencial suficientemente
alta (V₀ > E + mc²), a equação de Dirac prevê que ele **atravessa
sem atenuação exponencial** — ao contrário do tunelamento normal.

Schrödinger: T → 0 quando V₀ → ∞ (barreira mais alta = menos transmissão)
Dirac:       T > 0 para qualquer V₀ > E + mc²  (paradoxo!)

---

## Os três regimes físicos

```
       V₀ < E + mc²          E + mc² < V₀ < E + 3mc²       V₀ > E + 3mc²
    ┌──────────────────┐    ┌──────────────────────────┐    ┌─────────────┐
    │  TUNELAMENTO     │    │  OSCILAÇÃO/RESSONÂNCIAS  │    │   KLEIN     │
    │  T → 0 (exp)     │    │  T oscila entre 0 e 1    │    │   T → 1     │
    │  = Schrödinger   │    │  ressonâncias de Fabry-  │    │   100% !    │
    │                  │    │  Pérot relativísticas     │    │             │
    └──────────────────┘    └──────────────────────────┘    └─────────────┘
```

---

## Por que isso acontece?

### Intuição física (interpretação moderna)

Quando V₀ > E + mc², o potencial é tão alto que dentro da barreira a
"energia efetiva" do elétron fica abaixo de −mc². Nessa região, os
**estados de energia negativa** (antipartículas — o "mar de Dirac") ficam
acessíveis. O que ocorre é:

1. O elétron incidente "cria" um par elétron-pósitron na parede esquerda
2. O pósitron propaga-se para a esquerda (reflexão)  
3. O elétron do par propaga-se para a direita (transmissão)
4. Resultado: transmissão sem atenuação exponencial

### Matematicamente

O momento dentro da barreira é:

```
q = √[(E − V₀ + m)(E − V₀ − m)]
```

Para V₀ < E + m:  ambos os fatores têm sinais opostos → q imaginário (evanescente)
Para V₀ > E + m:  ambos os fatores são negativos → produto positivo → **q real!**

Quando q é real, a solução dentro da barreira são **ondas propagantes**,
não evanescentes. Daí a transmissão não-zero.

---

## Fórmulas de transmissão

### Barreira retangular — Dirac

```
T = |denom|⁻²

denom = cos(qd) − (i/2)(αq/αk + αk/αq) sin(qd)

αk = k/(E+m)      (helicidade fora)
αq = (E−V₀)/q     (helicidade dentro)
```

### Barreira retangular — Schrödinger (referência)

```
T = [1 + V₀²sin²(q'd) / (4E(E−V₀))]⁻¹,   q' = √(2m(E−V₀))
```

Para V₀ > E: q' imaginário → sin(q'd) = i·sinh(|q'|d) → T decai exponencialmente.

---

## Instalação e uso

```bash
pip install numpy matplotlib
```

```bash
# Roda tudo (figuras estáticas + animação)
python dirac_klein.py

# Muda a altura da barreira
python dirac_klein.py --V0 3.0   # regime de ressonâncias
python dirac_klein.py --V0 5.0   # regime de Klein (padrão)
python dirac_klein.py --V0 0.5   # tunelamento normal

# Só figuras estáticas (sem animação)
python dirac_klein.py --static

# Salva animação como GIF
python dirac_klein.py --save
```

### Todos os parâmetros

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `--k0`    | 1.0    | Momento central do pacote (mc) |
| `--V0`    | 5.0    | Altura da barreira (mc²) |
| `--d`     | 4.0    | Largura da barreira |
| `--mass`  | 1.0    | Massa (unidades naturais) |
| `--static`| False  | Só figuras estáticas |
| `--save`  | False  | Salvar animação como GIF |

---

## O que observar na figura estática

| Painel | O que mostra |
|--------|-------------|
| T(V₀) superior | Dirac (azul) vs Schrödinger (laranja): divergem radicalmente para V₀ > E+mc² |
| T(k) superior dir. | Para k fixo: Schrödinger → 0, Dirac > 0 no regime Klein |
| Densidade inferior | |ψ|² após atravessar: regime 1 (atenuado), 2 (parcial), 3 (Klein: pacote intacto) |

### O que observar na animação

- **Pacote se aproxima** da barreira com momento k₀
- **Na barreira**: sem decaimento exponencial (ondas propagantes dentro)
- **Após a barreira**: pacote emerge com amplitude significativa
- **Transmissão em tempo real** T_inst exibida no gráfico

---

## Conexão com física moderna

O paradoxo de Klein não é apenas curiosidade teórica:

- **Grafeno**: elétrons se comportam como férmions de Dirac sem massa (v_F ≈ c/300). O paradoxo de Klein foi **observado experimentalmente** em grafeno (Katsnelson et al., Nature Physics, 2006)
- **Transistores de tunelamento**: limita o confinamento de elétrons em nano-dispositivos de grafeno — não dá para confinar um férmion de Dirac sem massa com potencial eletrostático!
- **Física de partículas**: análogo à produção de pares por campo elétrico intenso (efeito Schwinger)

---

## Referências

- Klein, O. *Die Reflexion von Elektronen an einem Potentialsprung nach der relativistischen Dynamik von Dirac.* Z. Phys. 53 (1929).
- Calogeracos, A.; Dombey, N. *History and physics of the Klein paradox.* Contemp. Phys. 40, 313 (1999).
- Katsnelson, M. I. et al. *Chiral tunnelling and the Klein paradox in graphene.* Nature Physics 2, 620 (2006).
- Dombey, N.; Calogeracos, A. *Seventy years of the Klein paradox.* Phys. Rep. 315, 41 (1999).
