# Chiral $U(1)_{X}$ dark sector for micrOMEGAs

This repository implements the model of [arXiv:2506.04718](https://arxiv.org/abs/2506.04718)
(X. He, T. Nomura and N. Yokozaki, *Dark matter and dark radiation from chiral U(1) gauge symmetry*)
in [LanHEP](https://theory.sinp.msu.ru/~semenov/lanhep.html), producing CalcHEP model files for
[micrOMEGAs](https://lapth.cnrs.fr/micromegas/).

The dark sector has five left-handed Weyl fermions with $U(1)_{X}$ charges $(-9, -5, -1, 7, 8)$. These cancel the $U(1)_{X}^3$ and gravitational anomalies, which are the only anomalies, since SM fields carry no $U(1)_{X}$ charge.

The repository uses the paper's benchmark charges, $\left(Q_{1}, Q_{2}, Q_{M}\right) = (-9, -1, -5)$, with a dark Higgs $\phi_{q}$ of charge $q = 10$:

- **Dirac dark matter:** $\psi_{-9}$ and $\psi_{-1}$ form a Dirac fermion, odd under $Z_2$.
- **Majorana dark matter:** $\psi_{-5}$ gets a Majorana mass and is odd under $Z_2'$, so micrOMEGAs treats the model as two-component dark matter.
- **Dark radiation:** $\psi_{7}$ and $\psi_{8}$ stay massless.
- **Portal to the SM:** the $U(1)_{X}$ gauge boson, a dark photon $Z'$, couples to the SM through kinetic mixing with hypercharge, $-\tfrac{\sin\epsilon'}{2}\tilde B_{\mu\nu}X^{\mu\nu}$.

## Quick start

You need Linux with `gcc`, `gfortran`, `make`, `curl` and `md5sum`.

```bash
git clone https://github.com/restrepo/chiral-u1x-micromegas.git
cd chiral-u1x-micromegas
./setup.sh                      # downloads, builds, runs LanHEP, compiles (about 5 min)
./run.sh benchmarks/BP1.par     # one benchmark point (a few seconds)
```

`setup.sh` does the following, all inside `./build`:

1. It downloads micrOMEGAs 7.1.1 from [Zenodo](https://doi.org/10.5281/zenodo.20547050) and LanHEP 4.0.0 from its author's page, and checks both MD5 checksums.
2. It builds micrOMEGAs and LanHEP.
3. It creates the micrOMEGAs project `U1X` and runs LanHEP on `model/lanhep/source.lhp`.
4. It compiles `main.c`.

Run `./setup.sh --no-lanhep` to skip LanHEP and use the pre-generated CalcHEP files in `model/calchep`, which are identical to the LanHEP output. Set `BUILD=/path` to build somewhere else.

## Particles and parameters

| CalcHEP name | field | $Z_2\times Z_2'$ |
|---|---|---|
| `~X`, `~x` | Dirac DM: $P_{L}\psi_{D} = \psi_{-1}$, $P_{R}\psi_{D} = (\psi_{-9})^c$ | $\left(-,\,+\right)$ |
| `~~M` | Majorana DM: $P_{L}\psi_{M} = \psi_{-5}$ | $\left(+,\,-\right)$ |
| `nX`, `NX` | massless Dirac field built from $\psi_{7}$ and $(\psi_{8})^c$ (dark radiation) | even |
| `hX` | dark Higgs | even |
| `Zp` | dark photon $Z'$ | even |

| parameter | meaning |
|---|---|
| `gX` | $U(1)_{X}$ gauge coupling $g_{X}$ |
| `epsX` | kinetic mixing $\sin\epsilon'$ (the paper's $\epsilon$) |
| `MZp` | physical $Z'$ mass, below or above $M_{Z}$ (see the note on the forbidden window) |
| `MhX` | dark Higgs mass |
| `sa` | $\sin\alpha$, the $h$–$h_{X}$ mixing from the Higgs portal; `sa = 0` is the paper's model |
| `MDm` | Dirac mass $M_{D} = y_{D} v_{q}$ |
| `MMj` | Majorana mass $M_{M} = y_{M} v_{q}$ |

The model computes these as derived, public parameters:

- $M_{Z}$ and $M_{X}$, from the exact diagonalization of the $Z$–$Z'$ mass matrix;
- $v_{q}$, from $M_{X}^2 = 2 g_{X}^2 q^2 v_{q}^2$;
- $\sin\chi$, the $Z$–$Z'$ mixing angle;
- $\epsilon_{\mathrm{DP}} = \sin\epsilon'\cos\theta_{W}$;
- the portal coupling $\lambda_{HX}$.

## Benchmark points

All three points have $\epsilon = 10^{-6}$, $M_{Z'} = 30\,\mathrm{GeV}$ and $M_{h_{X}} = 70\,\mathrm{GeV}$.

| point | $g_{X}$ | $M_{D}$ [GeV] | $M_{M}$ [GeV] | $\sin\alpha$ |
|---|---|---|---|---|
| BP1 | 0.013 | 14.9 | 40 | 0 |
| BP2 | 0.0142 | 40 | 14.9 | 0 |
| BP3 | 0.013 | 14.9 | 40 | 0.003 |

- **BP1** is Majorana-dominated, as in Fig. 4 of the paper, with $M_{D} = M_{Z'}/2 - 0.1\,\mathrm{GeV}$.
- **BP2** is Dirac-dominated, as in Fig. 2, with $M_{M} = M_{Z'}/2 - 0.1\,\mathrm{GeV}$.
- **BP3** is BP1 with the Higgs portal switched on.

The expected results, from micrOMEGAs 7.1.1, are below. The full outputs are in `benchmarks/expected/`.

| point | $\Omega h^2$ total | $\Omega h^2$ Dirac | $\Omega h^2$ Majorana | $\sigma_{\mathrm{SI},p}$ Dirac [pb] | $\sigma_{\mathrm{SI},p}$ Majorana [pb] | LZ |
|---|---|---|---|---|---|---|
| BP1 | 0.1203 | $2.1\times 10^{-4}$ | 0.1201 | $2.5\times 10^{-14}$ | 0 | allowed |
| BP2 | 0.1214 | 0.1211 | $3.1\times 10^{-4}$ | $3.2\times 10^{-14}$ | 0 | allowed |
| BP3 | 0.1203 | $2.1\times 10^{-4}$ | 0.1201 | $1.2\times 10^{-13}$ | $7.4\times 10^{-13}$ | allowed |

The $Z'$ decays mostly to dark radiation:

- at BP1 and BP3, $\mathrm{BR}(Z'\to n_{X}\bar n_{X}) = 95\%$ and $\mathrm{BR}(Z'\to \psi_{D}\bar\psi_{D}) = 5\%$;
- at BP2, where the Dirac channel is closed, $\mathrm{BR}(Z'\to n_{X}\bar n_{X}) = 99.97\%$.

To check a build, compare your output with the stored one:

```bash
./run.sh benchmarks/BP1.par | diff - benchmarks/expected/BP1.out
```

## Thermalization: Fig. 1 of the paper

`cosmology/thermalization.py` (Python 3 with `numpy`, `scipy` and `matplotlib`) reproduces Fig. 1 of the paper: the inverse-decay rate $\Gamma_{\mathrm{ID}}$ (Eqs. 20 and 22), the scattering rate $\Gamma_{I}$ (Eq. 24) and the Hubble rate $H$, for the four parameter sets of the figure. It also prints the temperatures $T_{1}$ and $T_{0}$ between which $\Gamma_{\mathrm{ID}} > H$ (Eq. 25).

```bash
cd cosmology
python3 thermalization.py                          # figures/fig1.pdf and thermalize1..4.pdf
python3 thermalization.py --normalization derived  # figures/fig1_derived.pdf
```

The paper drops the order-one factors of these rates, writing them with $\sim$ only. The default normalization, `paper`, uses the factors read off the published figure, whose curves were extracted from its vector graphics. With $g_{Y} = 0.357$, $Q_{Y} = 1$ and $x = M_{Z'}/T$ they are:

- $\Gamma_{\mathrm{ID}} = 0.4636\,\Gamma_{0}\,x$ for $T > M_{Z'}$;
- $\Gamma_{\mathrm{ID}} = 1.076\,\Gamma_{0}\,x^{3/2}e^{-x}$ for $T < M_{Z'}$;
- $\Gamma_{I}$ is Eq. (24) times $0.9657/(12\pi)$, at $s = 4\langle p\rangle^2$ with $\langle p\rangle = 3.15\,T$;
- $H$ uses $g_{*} = 106.75$ and the reduced Planck mass.

The script reproduces the published curves to 0.2%, except for the height of the narrow $Z'$ resonance spike in $\Gamma_{I}$, which depends on the temperature grid of the plot.

The `derived` normalization uses instead the full thermal average of Eq. (20), valid at any temperature, and the coefficient of Eq. (24) derived from the cross section,

$$
\Gamma_{I} \simeq \frac{4\pi}{3}\,\frac{\alpha_{X}\alpha_{Y}\epsilon^2Q_{X}^2\,\bar s}{|D(\bar s)|^2}\,n_{f,\mathrm{eq}},
\qquad \bar s = 2\langle p\rangle^2 .
$$

The two derivations are in the appendices of [restrepo/dark_radiation](https://github.com/restrepo/dark_radiation). Only the upper edge of the window changes: for $M_{Z'} = 10\,\mathrm{GeV}$, $g_{X} = 0.01$ and $\epsilon = 10^{-6}$, the dark sector is thermalized between 0.71 and 38 GeV with the paper's factors, and between 0.71 and 52 GeV with the derived ones.

## Implementation notes

- **Exact neutral gauge sector.** The kinetic mixing is removed with $X = \tilde Z'/\cos\epsilon'$ and $\tilde B = B - \tan\epsilon'\,\tilde Z'$. The $(\tilde Z, \tilde Z')$ mass matrix is then diagonalized exactly, taking the physical $M_{Z'}$ as input. LanHEP's `CheckMasses` confirms that the masses are diagonal. Eq. (17) of the paper for $\sin 2\chi$ is exact, because its denominator uses the physical masses.
- **Heavy $Z'$.** The diagonalization holds for $M_{Z'}$ on either side of $M_{Z0} = M_{W}/c_{W} = 91.292\,\mathrm{GeV}$, the $Z$ mass before mixing; above it, the mixing angle $\chi$ changes sign. Points with $M_{Z'}$ from 200 to 1000 GeV run normally, and channels such as $Z'\to W^+W^-$ and $Z'\to Z\,h$ open automatically. The two eigenvalues repel, so there is no solution ($M_{X}^2 < 0$) in a narrow window just above $M_{Z0}$: 91.2921 to 91.2922 GeV for $\epsilon = 10^{-3}$, to 91.2932 GeV for $\epsilon = 10^{-2}$, and to 91.3178 GeV for $\epsilon = 0.05$. `main.c` stops with a message there. Below $M_{Z0}$ there is always a solution. Since $M_{W}$ and $s_{W}$ are the inputs, the physical $M_{Z}$ is shifted by $\Delta M_{Z} \simeq -b^2/[2M_{Z}(M_{Z'}^2 - M_{Z}^2)]$, with $b = M_{Z0}^2 s_{W}\tan\epsilon'$. This is negligible far from the resonance, but near it the shift exceeds the 2.1 MeV error on $M_{Z}$: at $\epsilon = 10^{-2}$ it is $-12\,\mathrm{MeV}$ for $M_{Z'} = 95\,\mathrm{GeV}$, and for $M_{Z'}$ just below $M_{Z0}$ the two states swap character and $M_{Z}$ comes out near 111 GeV. `main.c` prints a warning whenever $|\Delta M_{Z}| > 2.1\,\mathrm{MeV}$; such points need electroweak precision constraints or an input scheme that fixes the physical $M_{Z}$.
- **Unitary gauge.** The model is written in unitary gauge, and `main.c` sets `ForceUG=1`.
- **Higgs portal (an extension).** The paper's potential is $\lambda(|\phi_{q}|^2 - v_{q}^2)^2$. This implementation adds $\lambda_{HX}(|H|^2 - v^2/2)(|\phi_{q}|^2 - v_{q}^2)$ and uses $(M_{h}, M_{h_{X}}, \sin\alpha)$ as inputs. With `sa = 0` it reduces to the paper's model.
- **Light quarks.** The $u$ and $d$ quarks have mass `Mq = 0.05` GeV, with matching Yukawa couplings, because micrOMEGAs' direct-detection module needs nonzero light-quark masses.
- **No spin-dependent signal.** The Majorana spin-dependent cross section is exactly zero. With pure kinetic mixing, the axial couplings to quarks from $Z$ and $Z'$ exchange cancel at zero momentum transfer, so at low energy the $Z'$ couples only to the electromagnetic current. This differs from the nonzero $\sigma_{\mathrm{SD}}$ shown in Fig. 3 of the paper. With a nonzero portal, the Majorana state gets a spin-independent cross section from Higgs exchange instead (BP3).
- **Source of the SM part.** `model/lanhep/sm.inc` is micrOMEGAs' own LanHEP Standard Model (`STFM/work/lanhep/sm.inc`), with the $Z'$, the mixing and the light-quark masses added.

## Files

```
setup.sh, run.sh            download/build and run scripts
main.c                      micrOMEGAs driver: relic density (darkOmega2), nucleon
                            cross sections, sigma_SI,eff, LZ exclusion, Zp/hX/h decays
model/lanhep/               LanHEP source (source.lhp, sm.inc, u1x.inc)
model/calchep/              CalcHEP model files generated by LanHEP 4.0.0
benchmarks/                 BP1-BP3 parameter files and expected outputs
cosmology/thermalization.py Fig. 1 of the paper: Gamma_ID, Gamma_I and H(T)
cosmology/figures/          its output, fig1.pdf (paper factors) and fig1_derived.pdf
```

## References

- X. He, T. Nomura, N. Yokozaki, [arXiv:2506.04718](https://arxiv.org/abs/2506.04718).
- G. Alguero et al., micrOMEGAs 6.0, [Comput. Phys. Commun. 299 (2024) 109133](https://doi.org/10.1016/j.cpc.2024.109133), [arXiv:2312.14894](https://arxiv.org/abs/2312.14894).
- A. Semenov, LanHEP 4.0, [Comput. Phys. Commun. 201 (2016) 167](https://doi.org/10.1016/j.cpc.2016.01.003), [arXiv:1412.5016](https://arxiv.org/abs/1412.5016).

## License

This repository is distributed under the GNU General Public License v3.0; see `LICENSE`.
micrOMEGAs and LanHEP are downloaded from their authors and keep their own licenses.
