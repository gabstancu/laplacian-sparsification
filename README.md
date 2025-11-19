# Laplacian Sparsification

```math
(1 + \epsilon)x^\top L x \leq x^\top \tilde L x \leq (1 + \epsilon)x^\top L x,\quad \forall x \,\bot\,1
```

```math
Lx = b, \quad \tilde L x = b.
```

```math
L = D - W
```

Assembling the right hand-side (RHS):

```math
b = - h ^2 f(x, y) + \sum_{ \substack{j \in B \\ j \sim i} } w_{ij} g_j
```
