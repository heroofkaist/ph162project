# Where every formula comes from

This note derives the model step by step, from a quantum state to the full Hamiltonian and to the quantities that I plot.
Each step says what is a postulate, what is taken from a measurement, what is taken from a reference, and what is my own choice.
Numbers in square brackets refer to the list at the end.

Conventions: $e>0$ is the elementary charge, the mobile electron has charge $-e$, energies are in μeV, times are in ns, and $\hbar = 0.6582119569$ μeV·ns.

---

## 1. A quantum state is a list of amplitudes

**Postulate.** The state $|\psi\rangle$ changes in time according to the Schrödinger equation

$$ i\hbar\,\frac{d}{dt}|\psi\rangle = \hat H\,|\psi\rangle . $$

**Choice of base states.** For the mobile electron I take $|j\rangle$ = "the electron is on site $j$". The base states are orthonormal and complete:

$$ \langle i|j\rangle=\delta_{ij},\qquad \sum_j |j\rangle\langle j| = 1 . $$

**Expansion.** Any state is $|\psi\rangle=\sum_j C_j|j\rangle$ with $C_j=\langle j|\psi\rangle$. The probability of finding the electron on site $j$ is $|C_j|^2$.

**Projection.** Multiply the Schrödinger equation by $\langle i|$ from the left and insert the completeness relation:

$$ i\hbar\,\frac{dC_i}{dt}=\langle i|\hat H\sum_j|j\rangle\langle j|\psi\rangle=\sum_j H_{ij}\,C_j,\qquad H_{ij}=\langle i|\hat H|j\rangle . $$

This is Feynman's form of quantum mechanics [1, ch. 8]. $H_{ij}$ is ($i\hbar$ times) the amplitude per unit time to go from $j$ to $i$. From here on the whole task is to write down the matrix $H$.

## 2. The hopping Hamiltonian

**Entries.** In a perfect crystal all sites are identical, so
* the diagonal entries are all equal: $H_{jj}=E_0$;
* the electron can jump to a nearest neighbor with the same amplitude everywhere: $H_{j\pm1,\,j}=-A$;
* all other entries are neglected (nearest-neighbor approximation).

**Where $A$ comes from.** If $\phi_j(x)$ is the atomic orbital on site $j$ and $\hat h$ is the one-electron energy operator (kinetic energy plus lattice potential), then $-A=\langle\phi_{j+1}|\hat h|\phi_j\rangle$. It is the overlap of two neighboring orbitals through the potential, and it falls off exponentially with the lattice spacing [2]. In this project $A$ is a parameter.

**Site equation.** Inserting the entries into the equation of section 1 gives [1, Eq. 13.3]

$$ i\hbar\,\frac{dC_j}{dt}=E_0C_j-A\,C_{j+1}-A\,C_{j-1}. $$

**Operator form.** Define the shift $\hat T|j\rangle=|j+1\rangle$, so $\hat T=\sum_j|j+1\rangle\langle j|$ and $\hat T^\dagger|j\rangle=|j-1\rangle$. Then

$$ \hat H_{\rm hop}=E_0-A\,(\hat T+\hat T^\dagger). $$

Check: $\langle j\pm1|\hat H_{\rm hop}|j\rangle=-A$ and $\langle j|\hat H_{\rm hop}|j\rangle=E_0$. I set $E_0=0$, which only shifts the zero of energy.

**Solutions.** Try $C_j=e^{ikjb}e^{-iEt/\hbar}$, where $b$ is the lattice spacing. The site equation gives

$$ E=E_0-A\left(e^{ikb}+e^{-ikb}\right)=E_0-2A\cos kb , $$

which is [1, Eq. 13.13]. Every state of definite energy is a wave over the whole crystal. Therefore this term alone cannot hold an electron in place. Feynman shows that a bound state appears only at a lattice imperfection [1, section 13-7].

**Ring of three sites.** With $j$ counted modulo 3,

$$ \hat H_{\rm hop}=\begin{pmatrix}0&-A&-A\\-A&0&-A\\-A&-A&0\end{pmatrix},\qquad kb=0,\ \pm\tfrac{2\pi}{3},\qquad E=-2A,\ +A,\ +A . $$

The circulating state $C_j=e^{2\pi i j/3}/\sqrt3$ is an eigenstate with $E=+A$. Without the atoms it does not change in time at all. This makes it a clean starting state: every change that I see later is caused by the atoms.

## 3. One atom: three levels and the dipole coupling

**Levels.** I keep three levels of each atom and take the 2s energy as zero:

$$ \hat H_{\rm at}=\begin{pmatrix}0&0&0\\0&-\Delta&0\\0&0&-E_{{\rm Ly}\alpha}\end{pmatrix}\quad\text{in the order }(2s,\,2p,\,1s). $$

* $\Delta$ is the Lamb shift: $2p_{1/2}$ lies below $2s_{1/2}$ by $1057.8298$ MHz [4], so $\Delta=h\times1057.8298\ \text{MHz}=4.3748$ μeV. The reference run uses the rounded value 4.37 μeV.
* $E_{{\rm Ly}\alpha}=10.2$ eV corresponds to the 121.567 nm Lyman-α line [5].
* 2s is metastable. It cannot decay to 1s by one electric-dipole photon, and its two-photon lifetime is 0.12 s [3]. This is far longer than every time in the model, so I neglect it.

**Dipole energy from electrostatics.** Let an external potential $\varphi(\mathbf r)$ act on the atom. The proton ($+e$) sits at the origin and the bound electron ($-e$) at $\mathbf r$:

$$ U=e\,\varphi(0)-e\,\varphi(\mathbf r)\approx-e\,\mathbf r\cdot\nabla\varphi=e\,\mathbf E\cdot\mathbf r , $$

because the potential changes little over the size of an atom and $\mathbf E=-\nabla\varphi$. The constant parts cancel because the atom is neutral. For a field along $z$: $\hat H_E=eE\hat z$.

**Parity selects the matrix elements.** Under $\mathbf r\to-\mathbf r$ the operator $z$ changes sign. The states 2s and 1s are even and 2p is odd. An integral of an odd function vanishes, so

$$ \langle 2s|z|2s\rangle=\langle 2p|z|2p\rangle=\langle 2s|z|1s\rangle=0,\qquad\langle 2s|z|2p\rangle\neq0 . $$

The field also connects 2p and 1s, but these levels are 10.2 eV apart, so this mixing is negligible for the fields considered here. I keep 1s only as the final state of the decay.

**The 2s-2p integral.** With the hydrogen wave functions $R_{20}=\tfrac{1}{\sqrt2}a_0^{-3/2}\left(1-\tfrac{r}{2a_0}\right)e^{-r/2a_0}$ and $R_{21}=\tfrac{1}{\sqrt{24}}a_0^{-3/2}\tfrac{r}{a_0}e^{-r/2a_0}$:

$$ \langle 2s|z|2p,m{=}0\rangle=\underbrace{\int_0^\infty R_{20}R_{21}\,r^3dr}_{-3\sqrt3\,a_0}\ \underbrace{\int Y_{00}\cos\theta\,Y_{10}\,d\Omega}_{1/\sqrt3}=-3a_0 . $$

**Fine structure.** The level next to $2s_{1/2}$ is $2p_{1/2}$. Its state with $m_j=+\tfrac12$ is

$$ |2p_{1/2},+\tfrac12\rangle=-\sqrt{\tfrac13}\,|2p,m{=}0\rangle|\!\uparrow\rangle+\sqrt{\tfrac23}\,|2p,m{=}1\rangle|\!\downarrow\rangle , $$

and only the first part has the same spin as $|2s\rangle|\!\uparrow\rangle$. Therefore $|\langle2s_{1/2}|z|2p_{1/2}\rangle|=3a_0/\sqrt3=\sqrt3\,a_0$. This is the matrix element used by Lamb and Retherford [6].

**Stark matrix.** In the basis $(2s_{1/2},\,2p_{1/2})$ the atom in a static field is described by

$$ \begin{pmatrix}0&V\\V&-\Delta\end{pmatrix},\qquad V=\sqrt3\,e\,a_0\,E . $$

The overall sign of $V$ can be absorbed into the phase of the 2p state. Numbers: $ea_0E=0.0529$ μeV at 1 kV/m, so $V=0.0917$ μeV at 1 kV/m. The code computes both integrals numerically (`model/atom.py`).

My first portfolio estimate used $V\approx3ea_0E$, which ignores the fine structure. It is larger by $\sqrt3$ and does not change any conclusion.

## 4. Decay: from a lifetime to a quantum jump

**Measured law.** An atom in 2p survives with probability $e^{-\Gamma t}$, with $\Gamma=6.2649\times10^{8}\ \text{s}^{-1}$, a lifetime of 1.596 ns [5]. The ring model uses the rounded value $\Gamma=1/(1.6\ \text{ns})=0.625\ \text{ns}^{-1}$.

**Complex energy.** A probability is the square of an amplitude, so the 2p amplitude is

$$ C_{2p}(t)\propto e^{-iE_{2p}t/\hbar}\,e^{-\Gamma t/2}=\exp\!\left[-\frac{i}{\hbar}\left(E_{2p}-\frac{i\hbar\Gamma}{2}\right)t\right]. $$

Before a photon is emitted the atom therefore evolves as if 2p had the energy $E_{2p}-i\hbar\Gamma/2$:

$$ \hat H_{\rm eff}=\begin{pmatrix}0&V\\V&-\Delta-i\hbar\Gamma/2\end{pmatrix}. $$

**Weak field: the 2s state borrows a little decay.** The two amplitude equations are

$$ i\hbar\dot C_{2s}=V\,C_{2p},\qquad i\hbar\dot C_{2p}=V\,C_{2s}-\left(\Delta+\tfrac{i\hbar\Gamma}{2}\right)C_{2p}. $$

For $V\ll\Delta$ the 2p amplitude is small and follows the slowly changing 2s amplitude, so I set $\dot C_{2p}\approx0$:

$$ C_{2p}\approx\frac{V}{\Delta+i\hbar\Gamma/2}\,C_{2s}\quad\Longrightarrow\quad\Gamma_{\rm eff}=\Gamma\left|\frac{C_{2p}}{C_{2s}}\right|^2=\Gamma\,\frac{V^2}{\Delta^2+(\hbar\Gamma/2)^2}. $$

The rate grows with the square of the field. This $E^2$ law is observed in metastable hydrogen beams [7]. At 1 kV/m: $V=0.0917$ μeV, $\Delta=4.3748$ μeV, $\hbar\Gamma/2=0.2062$ μeV, so $\Gamma_{\rm eff}=2.74\times10^{5}\ \text{s}^{-1}$ and the atom needs 3.64 μs to emit.

**Any field: exact rate and its ceiling.** Without the weak-field assumption the decay rate of the 2s-like state is

$$ \Gamma_q=-\frac{2}{\hbar}\,{\rm Im}\,E_{\rm slow}, $$

where $E_{\rm slow}$ is the eigenvalue of $\hat H_{\rm eff}$ with the smaller decay. For $V\gg\Delta$ the eigenstates are $(|2s\rangle\pm|2p\rangle)/\sqrt2$. Each of them is half 2p, so each decays at $\Gamma/2$. Hence

$$ \Gamma_q\le\frac{\Gamma}{2}=3.1\times10^{8}\ \text{s}^{-1}\quad\text{for every field.} $$

The two energies $V$ and $\Delta$ are equal at $E=48$ kV/m. Well above this field a stronger field does not make the atom emit faster.

**The jump.** In a short time $dt$ the photon is emitted with probability $\Gamma\,|C_{2p}|^2\,dt$, and then the atom is found in 1s. This is described by the jump operator [8]

$$ \hat L=\sqrt\Gamma\,|1s\rangle\langle2p| ,\qquad \hat L^\dagger\hat L=\Gamma\,|2p\rangle\langle2p| . $$

Note that $-\tfrac{i\hbar}{2}\hat L^\dagger\hat L$ is exactly the term $-i\hbar\Gamma/2$ on 2p found above.

**Averaging gives the Lindblad equation.** In $dt$ there are two possibilities. Without a photon the state evolves with $\hat H_{\rm eff}=\hat H-\tfrac{i\hbar}{2}\sum_i\hat L_i^\dagger\hat L_i$. With a photon from atom $i$ the state becomes $\hat L_i\rho\hat L_i^\dagger$. The sum of both is

$$ \rho(t+dt)=\left(1-\tfrac{i}{\hbar}\hat H_{\rm eff}dt\right)\rho\left(1+\tfrac{i}{\hbar}\hat H_{\rm eff}^\dagger dt\right)+dt\sum_i\hat L_i\rho\hat L_i^\dagger , $$

and to first order in $dt$

$$ \frac{d\rho}{dt}=-\frac{i}{\hbar}[\hat H,\rho]+\sum_i\left(\hat L_i\rho\hat L_i^\dagger-\tfrac12\{\hat L_i^\dagger\hat L_i,\rho\}\right). $$

This is the Lindblad equation [8, 9]. It keeps the trace equal to 1 and keeps $\rho$ positive.

## 5. The field comes from the moving electron

**Coulomb's law.** A point charge $-e$ at position $\mathbf r_j$ produces at the atom at $\mathbf r_i$ the field

$$ \mathbf E=\frac{-e}{4\pi\varepsilon_0}\,\frac{\mathbf r_i-\mathbf r_j}{|\mathbf r_i-\mathbf r_j|^3}, $$

which points from the atom toward the electron. Along the chain this means: an electron on site $i+1$ gives $+E_0$ at atom $i$, and an electron on site $i-1$ gives $-E_0$, with $E_0=e/(4\pi\varepsilon_0b^2)$ for an unscreened charge.

**Model choices.** I keep only the two nearest neighbors, and I set the field of an electron on the same site to zero. In a real crystal the field is screened, so I treat the strength as a parameter.

**The field is an operator.** The electron position is a quantum variable. With the projector $\hat n_j=|j\rangle\langle j|$:

$$ \hat E_i=E_0\left(\hat n_{i+1}-\hat n_{i-1}\right). $$

Check: $\hat E_i|j\rangle$ gives $+E_0|j\rangle$ for $j=i+1$, $-E_0|j\rangle$ for $j=i-1$ and zero otherwise.

**Interaction.** Insert this field into the dipole energy of section 3, and project $\hat z_i$ on the levels 2s and 2p of atom $i$, where $\hat z_i\to\sqrt3\,a_0\hat X_i$ with $\hat X_i=|2s\rangle_i\langle2p|+|2p\rangle_i\langle2s|$:

$$ \hat H_{\rm int}=\sum_i e\,\hat E_i\otimes\hat z_i=g\sum_i\left(\hat n_{i+1}-\hat n_{i-1}\right)\otimes\hat X_i,\qquad g=\sqrt3\,e\,a_0E_0 . $$

This one term does two things. It mixes the atom when the electron is near, and it acts back on the electron, because it depends on the electron position.

**Size of $g$.** The reference run uses $g=1.6$ μeV. This is the mixing in a field of 17.5 kV/m, which is the bare Coulomb field of an electron 290 nm away. A bare electron at 1 nm would give $1.4\times10^{9}$ V/m. That is far outside the two-level model, so a real material needs a screened and non-perturbative treatment.

## 6. The full model

$$ \hat H=\underbrace{-A(\hat T+\hat T^\dagger)}_{\text{section 2}}+\underbrace{\sum_i\hat H_{{\rm at},i}}_{\text{section 3}}+\underbrace{g\sum_i(\hat n_{i+1}-\hat n_{i-1})\otimes\hat X_i}_{\text{section 5}},\qquad\hat L_i=\sqrt\Gamma\,|1s\rangle_i\langle2p| . $$

**Joint space.** A base state names the electron site and the level of each atom: $|j;a_0a_1a_2\rangle$. There are $3\times3^3=81$ of them. An operator that acts on one part is multiplied by the identity on the other parts (tensor product).

**Example of one column.** For the electron on site 0 and all atoms in 2s, atom 1 sees $\hat n_2-\hat n_0=-1$ and atom 2 sees $\hat n_0-\hat n_1=+1$:

$$ \hat H_{\rm int}|0;2s,2s,2s\rangle=g\left(|0;2s,2s,2p\rangle-|0;2s,2p,2s\rangle\right). $$

**The 1s energy cancels.** The Hamiltonian never changes which atoms are in 1s, and a jump puts the same atom into 1s on both sides of $\rho$. So $\rho$ never contains a coherence between states with different sets of 1s atoms, and the large energy $E_{{\rm Ly}\alpha}$ drops out of $[\hat H,\rho]$. I can set it to zero. The test `test_1s_energy_cancels` confirms it.

**Parameters of the reference run.** $A=1.0$ μeV, $g=1.6$ μeV, $\Delta=4.37$ μeV, $\Gamma=0.625\ \text{ns}^{-1}$. The values of $A$ and $g$ are illustrative. They are chosen so that hopping, mixing and decay happen on comparable time scales.

## 7. What I compute from the state

**Electron alone.** The atoms are not observed, so I sum over their states: $(\rho_{\rm el})_{ij}=\sum_a\langle i;a|\rho|j;a\rangle$.

**Coherence.** $C=\sum_{i<j}|(\rho_{\rm el})_{ij}|$. The off-diagonal elements carry the phase relation between sites. This is half of the $l_1$ measure of coherence [10]. For the circulating state each element is $1/3$, so $C=1$ at the start.

**Current.** From the site equation of section 2,

$$ \frac{d|C_j|^2}{dt}=J_{j-1\to j}-J_{j\to j+1},\qquad J_{j\to j+1}=\frac{2A}{\hbar}\,{\rm Im}\left(C_j^*C_{j+1}\right). $$

What leaves site $j$ enters a neighbor. As an operator, summed over the ring, $\hat J=\tfrac{iA}{\hbar}(\hat T-\hat T^\dagger)$. The interaction commutes with every $\hat n_j$, so it does not change this expression. For a wave $e^{ikjb}$ the current is proportional to $\sin kb$, in agreement with the group velocity $\tfrac{1}{\hbar}\tfrac{dE}{dk}$. I always divide by the initial value.

**Why the atoms reduce interference.** If two paths with amplitudes $a$ and $b$ leave the atoms in the states $|\chi_a\rangle$ and $|\chi_b\rangle$, the probability at the end is

$$ P=|a|^2+|b|^2+2\,{\rm Re}\left(a^*b\,\langle\chi_a|\chi_b\rangle\right). $$

The interference term is multiplied by the overlap of the atomic states [11]. A perfect record has zero overlap.

**How much the atoms record.** Let $m$ be the pattern of atoms in 1s (8 patterns) and $j$ the electron site. From the joint probabilities $p(j,m)$ I compute the mutual information [12]

$$ I=\sum_{j,m}p(j,m)\log_2\frac{p(j,m)}{p(j)\,p(m)} . $$

A perfect record of the site would give $\log_23=1.58$ bits.

## 8. Numerical method and checks

* The density matrix has $81\times81$ elements. I integrate the Lindblad equation with an explicit Runge-Kutta method of order 8 (`evolve`).
* Independent check: I build the $6561\times6561$ Liouvillian as a sparse matrix and apply its matrix exponential (`evolve_expm`). The two results at 12 ns differ by less than $3\times10^{-13}$.
* The trace stays 1, and $\rho$ stays Hermitian and positive at every sampled time.
* With $g=0$ nothing changes in time. Free hopping from one site follows the exact formula $(5+4\cos(3At/\hbar))/9$.
* All checks are in `tests/test_model.py`.

## 9. Results of the ring

| Quantity | Value |
|---|---|
| Current at 12 ns, full model | 21.7% of the initial current |
| Coherence $C$ at 12 ns, full model | 0.411 |
| Fraction of atoms in 1s at 12 ns | 38.6% |
| Current at 12 ns, emission off | 43.6% (it keeps swinging between 26% and 97%) |
| Probability of at least one emission by 12 ns | 70.7% |
| Mutual information between 1s pattern and electron site at 12 ns | 0.008 bit (never above 0.05 bit) |
| Electron site if only atom 0 has decayed | 26.9%, 37.1%, 36.0% |
| Atoms in 1s and current after 60 ns | 89%, 1.2% |
| Atoms in 1s after 300 ns | 99.99% |
| Long chain, 100 ns: probability within two sites of the start | 84% (free hopping: 4%) |
| Hydrogen at 1 kV/m: quench rate and quench time | $2.74\times10^{5}\ \text{s}^{-1}$, 3.64 μs |
| One hop with $A=1$ meV | 0.658 ps |
| Ceiling of the quench rate | $\Gamma/2=3.1\times10^{8}\ \text{s}^{-1}$, so the clocks match only for $A\lesssim0.21$ μeV |

With the unrounded gap 4.3748 μeV the full-model current at 12 ns is 21.8% instead of 21.7%.

## 10. The long chain: the electron digs its own well

The ring of three sites is too small to ask whether the electron stays in place, because its three atoms run out. I therefore put the same Hamiltonian on an open chain of $N$ sites and launch the electron on the central site.

**Why excited neighbors raise the energy.** Hold the electron next to one excited atom. The atom is then described by the Stark matrix of section 3 with $V=\pm g$:

$$ \begin{pmatrix}0&\pm g\\ \pm g&-\Delta\end{pmatrix}\quad\Longrightarrow\quad E_\pm=\tfrac12\left(-\Delta\pm\sqrt{\Delta^2+4g^2}\right). $$

The two levels repel each other. The 2s-like state moves **up** by

$$ E_+=\tfrac12\left(\sqrt{\Delta^2+4g^2}-\Delta\right)\approx\frac{g^2}{\Delta}=0.52\ \mu\text{eV}\quad(g=1.6,\ \Delta=4.37). $$

With the mixing angle $\tan2\theta=2g/\Delta$, the 2s-like state contains the fraction $\sin^2\theta=0.097$ of 2p. It decays slowly at $\Gamma\sin^2\theta=0.060\ \text{ns}^{-1}$, which is the quench rate of section 4. The 2p-like state decays fast at $\Gamma\cos^2\theta=0.56\ \text{ns}^{-1}$.

**Energy landscape.** For slow hopping the atoms follow the electron, and the energy of the electron on site $j$ is

$$ E(j)=E_+\times(\text{number of excited neighbors of }j). $$

An atom in 1s gives no shift, because the field cannot mix it. While all atoms are excited, every site has the energy $2E_+$ and the chain is uniform. When the two neighbors of one site have emitted, this site is lower than all others by $2E_+=1.05$ μeV. For $A=0.1$ μeV the well is ten times deeper than the hopping energy, so the electron is bound.

**Each photon is a position measurement.** The jump operator $\hat L_i$ keeps only the part of the state in which atom $i$ was in 2p. This part exists only where the electron was next to atom $i$. So an emission from atom $i$ places the electron on site $i\pm1$, and a second emission from atom $i\pm2$ completes the well.

**Method: quantum trajectories [8].** The joint space has $N\,3^N$ states, which is too large for a density matrix. Instead I follow single state vectors:

1. Between emissions the state evolves with $\hat H_{\rm eff}=\hat H-\tfrac{i\hbar}{2}\sum_i\hat L_i^\dagger\hat L_i$, and its norm falls.
2. When the squared norm reaches a random number $r\in(0,1)$, one atom emits. Atom $i$ is chosen with a probability proportional to its 2p population, and the state becomes $\hat L_i|\psi\rangle/\|\hat L_i|\psi\rangle\|$.
3. The average over many trajectories reproduces the Lindblad equation of section 4.

An atom in 1s is a classical record, so it leaves the state vector. I also allow at most two atoms in 2p at the same time. The code is `model/chain.py`.

**Checks.**
* On a chain of four sites the average of 1500 trajectories agrees with the exact Lindblad solution. The largest difference is 0.012, within 2.1 standard errors (`validate_chain.py`).
* Evolution without a photon is the same with at most two or at most four simultaneous 2p atoms (`tests/test_chain.py`).
* The result does not depend on the chain length: the probability within two sites of the start after 100 ns is 80%, 81% and 84% for 11, 21 and 31 sites, with a statistical error of 2%.

**Results** (31 sites, $A=0.1$ μeV, 400 trajectories, 100 ns).

| Quantity after 100 ns | Full model | Emission off | Free hopping |
|---|---|---|---|
| Probability within two sites of the start | 84% ± 2% | 12% | 4% |
| Rms distance from the start (sites) | 2.2 | 16 | 21 |
| Coherence $C$ of the electron | 0.2 (never above 0.6; statistical floor 0.06) | 17 | 26 |

* The electron sits between two 1s atoms with 93% probability. On the long chain the record does locate the electron.
* On average 2.6 photons are emitted in 100 ns.
* The trap needs the emission: with emission off no atom reaches 1s, no well forms, and the electron spreads.
* The electron starts on one site, so its coherence starts at 0. A spreading wave builds coherence. In the full model every photon measures the position, so the coherence stays small (`run_chain_coherence.py`).
* The trap needs slow hopping. The probability within two sites after 100 ns is 95%, 84%, 49% and 18% for $A$ = 0.05, 0.1, 0.2 and 0.4 μeV.

## 11. Limits of the model

* The atoms emit independently. For spacings much smaller than 121.6 nm, atoms share radiation modes and emit collectively [13].
* The field is short-ranged and has a fixed strength. A real crystal needs a screened field.
* Nothing excites the atoms again, so each atom emits at most once and a well stays once it is dug.
* The long chain has 31 sites and the run lasts 100 ns. The dark region around the electron grows slowly, so longer runs must decide whether the trap is permanent. Faster hopping needs longer chains.

## References

1. R. P. Feynman, R. B. Leighton and M. Sands, *The Feynman Lectures on Physics*, Vol. III, ch. 8 and ch. 13 (Addison-Wesley, 1965).
2. J. C. Slater and G. F. Koster, "Simplified LCAO Method for the Periodic Potential Problem," Phys. Rev. 94, 1498 (1954).
3. J. Shapiro and G. Breit, "Metastability of 2s States of Hydrogenic Atoms," Phys. Rev. 113, 179 (1959).
4. N. Bezginov et al., "A measurement of the atomic hydrogen Lamb shift and the proton charge radius," Science 365, 1007 (2019).
5. W. L. Wiese and J. R. Fuhr, "Accurate Atomic Transition Probabilities for Hydrogen, Helium, and Lithium," J. Phys. Chem. Ref. Data 38, 565 (2009).
6. W. E. Lamb Jr. and R. C. Retherford, "Fine Structure of the Hydrogen Atom. Part I," Phys. Rev. 79, 549 (1950).
7. A. Lejeune, L. Chérigier-Kovacic and F. Doveil, "Lyman-α radiation of a metastable hydrogen beam to measure electric fields," Appl. Phys. Lett. 99, 181502 (2011).
8. J. Dalibard, Y. Castin and K. Mølmer, "Wave-function approach to dissipative processes in quantum optics," Phys. Rev. Lett. 68, 580 (1992).
9. G. Lindblad, "On the generators of quantum dynamical semigroups," Commun. Math. Phys. 48, 119 (1976).
10. T. Baumgratz, M. Cramer and M. B. Plenio, "Quantifying Coherence," Phys. Rev. Lett. 113, 140401 (2014).
11. B.-G. Englert, "Fringe Visibility and Which-Way Information: An Inequality," Phys. Rev. Lett. 77, 2154 (1996).
12. C. E. Shannon, "A Mathematical Theory of Communication," Bell Syst. Tech. J. 27, 379 (1948).
13. R. H. Lehmberg, "Radiation from an N-Atom System. I. General Formalism," Phys. Rev. A 2, 883 (1970).
