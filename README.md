# Chiral $U(1)_X$ dark sector for micrOMEGAs

This repository implements the model of [arXiv:2506.04718](https://arxiv.org/abs/2506.04718)
(X. He, T. Nomura and N. Yokozaki, *Dark matter and dark radiation from chiral U(1) gauge symmetry*)
in [LanHEP](https://theory.sinp.msu.ru/~semenov/lanhep.html), producing CalcHEP model files for
[micrOMEGAs](https://lapth.cnrs.fr/micromegas/).

The dark sector has five left-handed Weyl fermions with $U(1)_X$ charges $(-9, -5, -1, 7, 8)$. These cancel the $U(1)_X^3$ and gravitational anomalies, which are the only anomalies, since SM fields carry no $U(1)_X$ charge.

The repository uses the paper's benchmark charges, $(Q_1, Q_2, Q_M) = (-9, -1, -5)$, with a dark Higgs $\phi_q$ of charge $q = 10$:

- **Dirac dark matter:** $\psi_{-9}$ and $\psi_{-1}$ form a Dirac fermion, odd under $Z_2$.
- **Majorana dark matter:** $\psi_{-5}$ gets a Majorana mass and is odd under $Z_2'$, so micrOMEGAs treats the model as two-component dark matter.
- **Dark radiation:** $\psi_{7}$ and $\psi_{8}$ stay massless.
- **Portal to the SM:** the $U(1)_X$ gauge boson, a dark photon $Z'$, couples to the SM through kinetic mixing with hypercharge, $-\tfrac{\sin\epsilon'}{2}\tilde B_{\mu\nu}X^{\mu\nu}$.

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
| `~X`, `~x` | Dirac DM: $P_L\psi_D = \psi_{-1}$, $P_R\psi_D = (\psi_{-9})^c$ | $(-, +)$ |
| `~~M` | Majorana DM: $P_L\psi_M = \psi_{-5}$ | $(+, -)$ |
| `nX`, `NX` | massless Dirac field built from $\psi_{7}$ and $(\psi_{8})^c$ (dark radiation) | even |
| `hX` | dark Higgs | even |
| `Zp` | dark photon $Z'$ | even |

| parameter | meaning |
|---|---|
| `gX` | $U(1)_X$ gauge coupling $g_X$ |
| `epsX` | kinetic mixing $\sin\epsilon'$ (the paper's $\epsilon$) |
| `MZp` | physical $Z'$ mass, which must be below $M_{Z}$ |
| `MhX` | dark Higgs mass |
| `sa` | $\sin\alpha$, the $h$–$h_X$ mixing from the Higgs portal; `sa = 0` is the paper's model |
| `MDm` | Dirac mass $M_D = y_D v_q$ |
| `MMj` | Majorana mass $M_M = y_M v_q$ |

The model computes these as derived, public parameters:

- $M_{Z}$ and $M_X$, from the exact diagonalization of the $Z$–$Z'$ mass matrix;
- $v_q$, from $M_X^2 = 2 g_X^2 q^2 v_q^2$;
- $\sin\chi$, the $Z$–$Z'$ mixing angle;
- $\epsilon_{\mathrm{DP}} = \sin\epsilon'\cos\theta_W$;
- the portal coupling $\lambda_{HX}$.

## Benchmark points

All three points have $\epsilon = 10^{-6}$, $M_{Z'} = 30\,\mathrm{GeV}$ and $M_{h_X} = 70\,\mathrm{GeV}$.

| point | $g_X$ | $M_D$ [GeV] | $M_M$ [GeV] | $\sin\alpha$ |
|---|---|---|---|---|
| BP1 | 0.013 | 14.9 | 40 | 0 |
| BP2 | 0.0142 | 40 | 14.9 | 0 |
| BP3 | 0.013 | 14.9 | 40 | 0.003 |

- **BP1** is Majorana-dominated, as in Fig. 4 of the paper, with $M_D = M_{Z'}/2 - 0.1\,\mathrm{GeV}$.
- **BP2** is Dirac-dominated, as in Fig. 2, with $M_M = M_{Z'}/2 - 0.1\,\mathrm{GeV}$.
- **BP3** is BP1 with the Higgs portal switched on.

The expected results, from micrOMEGAs 7.1.1, are below. The full outputs are in `benchmarks/expected/`.

| point | $\Omega h^2$ total | $\Omega h^2$ Dirac | $\Omega h^2$ Majorana | $\sigma_{\mathrm{SI},p}$ Dirac [pb] | $\sigma_{\mathrm{SI},p}$ Majorana [pb] | LZ |
|---|---|---|---|---|---|---|
| BP1 | 0.1203 | $2.1\times 10^{-4}$ | 0.1201 | $2.5\times 10^{-14}$ | 0 | allowed |
| BP2 | 0.1214 | 0.1211 | $3.1\times 10^{-4}$ | $3.2\times 10^{-14}$ | 0 | allowed |
| BP3 | 0.1203 | $2.1\times 10^{-4}$ | 0.1201 | $1.2\times 10^{-13}$ | $7.4\times 10^{-13}$ | allowed |

The $Z'$ decays mostly to dark radiation:

- at BP1 and BP3, $\mathrm{BR}(Z'\to n_X\bar n_X) = 95\%$ and $\mathrm{BR}(Z'\to \psi_D\bar\psi_D) = 5\%$;
- at BP2, where the Dirac channel is closed, $\mathrm{BR}(Z'\to n_X\bar n_X) = 99.97\%$.

To check a build, compare your output with the stored one:

```bash
./run.sh benchmarks/BP1.par | diff - benchmarks/expected/BP1.out
```

## Implementation notes

- **Exact neutral gauge sector.** The kinetic mixing is removed with $X = \tilde Z'/\cos\epsilon'$ and $\tilde B = B - \tan\epsilon'\,\tilde Z'$. The $(\tilde Z, \tilde Z')$ mass matrix is then diagonalized exactly, taking the physical $M_{Z'}$ as input. LanHEP's `CheckMasses` confirms that the masses are diagonal. Eq. (17) of the paper for $\sin 2\chi$ is exact, because its denominator uses the physical masses.
- **Unitary gauge.** The model is written in unitary gauge, and `main.c` sets `ForceUG=1`.
- **Higgs portal (an extension).** The paper's potential is $\lambda(|\phi_q|^2 - v_q^2)^2$. This implementation adds $\lambda_{HX}(|H|^2 - v^2/2)(|\phi_q|^2 - v_q^2)$ and uses $(M_h, M_{h_X}, \sin\alpha)$ as inputs. With `sa = 0` it reduces to the paper's model.
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
```

## References

- X. He, T. Nomura, N. Yokozaki, [arXiv:2506.04718](https://arxiv.org/abs/2506.04718).
- G. Alguero et al., micrOMEGAs 6.0, [Comput. Phys. Commun. 299 (2024) 109133](https://doi.org/10.1016/j.cpc.2024.109133), [arXiv:2312.14894](https://arxiv.org/abs/2312.14894).
- A. Semenov, LanHEP 4.0, [Comput. Phys. Commun. 201 (2016) 167](https://doi.org/10.1016/j.cpc.2016.01.003), [arXiv:1412.5016](https://arxiv.org/abs/1412.5016).
