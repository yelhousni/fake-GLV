# GLV + fake-GLV

This repo contains sage scripts to validate the ideas of this writeup: https://ethresear.ch/t/fake-glv-you-dont-need-an-efficient-endomorphism-to-implement-glv-like-scalar-multiplication-in-snark-circuits/20394

<img width="768" alt="image" src="https://github.com/user-attachments/assets/fcd69d93-4619-4ea0-8a5c-4a67e0aa5a50">

----

- *GLV:* $[s]P = [s1]P + [s2]\phi(P)$.
- *fake-GLV:* $[s]P == Q$ is equivalent to $[u]P - [v]Q == 0_E$ where $v\cdot s = u \pmod r$ and $|u|, |v| \leq \sqrt{r}$.
- *GLV + fake-GLV:* $[s]P == Q$ is equivalent to $[u_1]P + [u_2]\phi(P) - [v_1]Q - [v_2]\phi(Q) == 0_E$ where $(v_1+j\cdot v_2)\cdot (s_1+j\cdot s_2) = (u_1+j\cdot u_2) \mod (r_1+j\cdot r_2)$ in $\mathbf{Q}(j)$ and $|u_1|, |u_2|, |v_1|, |v_2| \leq r^{1/4}$.


