# Mar de Dirac e Antimatéria — Equação de Dirac

> **Nível avançado — Simulação 1 de 4**

---

## Os quatro módulos desta simulação

| Módulo | Física | Comando |
|--------|--------|---------|
| 1 — Espectro | Mar de Dirac, gap 2mc², interpretação de buraco | `--module 1` |
| 2 — Feynman | Propagador: e⁻ → frente, e⁺ → atrás no tempo | `--module 2` |
| 3 — Schwinger | Taxa de criação de pares por campo elétrico | `--module 3` |
| 4 — Animação | Par e⁺e⁻ se encontra → aniquilação | `--module 4` |

---

## 1. O espectro de Dirac e o Mar

A equação de Dirac livre tem dois ramos de soluções:

```
E = ±√(k²c² + m²c⁴)
```

- **Ramo positivo** (E > mc²): elétrons físicos
- **Gap proibido** (−mc² < E < +mc²): região de largura 2mc²
- **Ramo negativo** (E < −mc²): problema! Estados infinitos de energia decrescente

### A proposta de Dirac (1930)

Dirac postulou que o **vácuo é o mar de estados negativos completamente preenchido**
(princípio de exclusão de Pauli). Consequências:

1. Elétrons não podem cair para estados negativos (já estão cheios)
2. Um fóton de energia ≥ 2mc² pode excitar um elétron do mar para o ramo positivo
3. O **buraco** que fica no mar tem:
   - Carga oposta ao elétron: +e
   - Mesma massa: m
   - → **Pósitron** (antipartícula do elétron)

Anderson observou o pósitron experimentalmente em 1932, confirmando a previsão.

---

## 2. Interpretação de Feynman-Stückelberg

Feynman propôs uma interpretação equivalente mais elegante:

> **O pósitron é um elétron propagando para trás no tempo.**

O propagador de Feynman $S_F(x-y)$ naturalmente incorpora isso:

```
         { ⟨0|T{ψ(x)ψ̄(y)}|0⟩   para x⁰ > y⁰  → propagação causal (elétron)
S_F =  {
         { −⟨0|T{ψ̄(y)ψ(x)}|0⟩  para x⁰ < y⁰  → pósitron recuando no tempo
```

Graficamente: no diagrama de Feynman, uma linha de elétron apontando para
trás no tempo é **exatamente equivalente** a uma linha de pósitron para frente.
Isso simplifica enormemente os cálculos de QED.

---

## 3. Efeito Schwinger — Criação de Pares por Campo Elétrico

Um campo elétrico suficientemente forte pode "tirar" um elétron do mar de Dirac,
criando um par elétron-pósitron do vácuo. A taxa é:

```
Γ ∝ (eE)² exp(−π m²c³ / eEħ)
```

**Campo crítico de Schwinger:**

```
E_crit = m²c³ / eħ ≈ 1.32 × 10¹⁸ V/m
```

### Valores numéricos

| E/E_crit | Taxa Γ (u.n.) | Regime |
|----------|---------------|--------|
| 0.1 | 2.27 × 10⁻¹⁶ | Indetectável |
| 0.3 | 2.55 × 10⁻⁶ | Lasers futuros |
| 0.5 | 4.67 × 10⁻⁴ | Acessível com ELI |
| 1.0 | 4.32 × 10⁻² | Supercrítico |
| 2.0 | 8.32 × 10⁻¹ | Produção intensa |

O laser mais intenso do mundo (ELI-NP, Romênia) atinge ~10⁻³ E_crit.
A criação de pares de Schwinger ainda não foi diretamente observada, mas
é prevista para lasers da próxima geração (~2030).

### Mecanismo físico

O campo elétrico **inclina** os níveis de energia, fazendo com que o
ramo positivo e o ramo negativo se cruzem no espaço. O elétron
"tunela" pelo gap de 2mc², saltando do mar para o ramo positivo —
análogo ao tunelamento de Klein, mas em campo uniforme.

---

## 4. Aniquilação: e⁻ + e⁺ → 2γ

A simulação evolui um pacote gaussiano contendo dois componentes:

- **Elétron** (ramo E>0, k>0): move-se para a esquerda
- **Pósitron** (ramo E<0, k<0 no mar → k>0 físico): move-se para a direita

Quando se sobrepõem, a densidade $\rho_+ \cdot \rho_-$ mede a amplitude
de aniquilação. Em QED completa, isso emitiria dois fótons de energia mc² cada.

O estado do espinor é:

```
Ψ = Ψ_e⁻ + Ψ_e⁺

Ψ_e⁻: ψ₁ = f(x)·e^{ik₀x},  ψ₂ = (k₀/E+m)·f(x)·e^{ik₀x}    [E>0]
Ψ_e⁺: ψ₁ = g(x)·e^{-ik₀x}, ψ₂ = (k₀/E+m)·g(x)·e^{-ik₀x}   [E<0]
```

A projeção $\rho_\pm$ separa as contribuições dos dois ramos no espaço de momentos.

---

## Instalação e uso

```bash
pip install numpy matplotlib
```

```bash
# Tudo (figura estática + animação)
python dirac_antimatter.py

# Só espectro e mar de Dirac
python dirac_antimatter.py --module 1

# Só efeito Schwinger
python dirac_antimatter.py --module 3

# Animação do par se aniquilando
python dirac_antimatter.py --module 4

# Salvar animação
python dirac_antimatter.py --save

# Massa diferente
python dirac_antimatter.py --mass 0.5 --k0 2.0
```

---

## O que observar

### Figura estática (6 painéis)

| Painel | O que mostra |
|--------|-------------|
| Espectro + mar | Os dois ramos ±E(k), gap 2mc², elétron e buraco anotados |
| Propagador Feynman | Diagrama espaço-tempo: e⁻ vai à direita, e⁺ "vem do futuro" |
| Taxa Schwinger | Γ(E) em escala log — exponencialmente suprimida para E < E_crit |
| Mecanismo Schwinger | Níveis inclinados pelo campo, região de tunelamento |
| Par inicial | ρ total (verde), ρ+ elétron (azul), ρ- pósitron (laranja) |
| Espectro k | Dois picos em ±k₀ — assinatura do par |

### Animação (3 painéis)

- **Painel superior**: densidade total — dois picos se aproximando
- **Painel central**: componente e⁻ (ρ₊) — pico movendo-se para a esquerda
- **Painel inferior**: componente e⁺ (ρ₋) — pico movendo-se para a direita
- Quando se sobrepõem: mensagem de aniquilação e overlap calculado em tempo real

---

## Referências

- Dirac, P.A.M. *A theory of electrons and protons.* Proc. R. Soc. A 126 (1930).
- Anderson, C.D. *The positive electron.* Phys. Rev. 43 (1933). — Nobel 1936.
- Feynman, R.P. *The theory of positrons.* Phys. Rev. 76 (1949).
- Schwinger, J. *On gauge invariance and vacuum polarization.* Phys. Rev. 82 (1951).
- Ruffini, R. et al. *On the pair electromagnetic pulse.* Phys. Rep. 487 (2010).
