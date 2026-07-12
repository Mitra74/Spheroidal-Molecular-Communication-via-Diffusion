# Spheroidal Molecular Communication via Diffusion: Signalling Between Homogeneous Cell Aggregates

Code accompanying **Chapter 3** of the PhD thesis *"Diffusion-Based Molecular Communication in Discrete Heterogeneous Environments"* (Mitra Rezaei, University of Warwick).

## Papers
- **Journal**: Rezaei, M., Arjmandi, H., Zoofaghari, M., Kanebratt, K., Vilén, L., Janzén, D., Gennemark, P., and Noel, A. "Spheroidal Molecular Communication via Diffusion: Signaling Between Homogeneous Cell Aggregates." *IEEE Transactions on Molecular, Biological, and Multi-Scale Communications*, vol. 10, no. 1, pp. 197–210, Mar. 2024. [DOI](https://ieeexplore.ieee.org/abstract/document/10438024)
- **Conference** (preliminary version): Arjmandi, H., Zoofaghari, M., Rezaei, M., et al. "Diffusive Molecular Communication with a Spheroidal Receiver for Organ-on-Chip Systems." *IEEE International Conference on Communications (ICC)*, Roma, Italy, 2023. [DOI](https://ieeexplore.ieee.org/abstract/document/10279128)

## Overview
This repository implements a generalized analytical framework for diffusion-based molecular communication where both transmitter and receiver are modelled as **spheroids** — 3-D multicellular aggregates of the kind used in organ-on-chip experiments — rather than idealized point sources or sinks. Each spheroid is treated as a homogenized porous medium, characterized by a porosity parameter and an effective diffusion coefficient, with a boundary condition that permits a concentration discontinuity at the spheroid surface.

The Green's function for concentration (GFC) is derived analytically, in series form, for both a point source releasing molecules inside a spheroidal transmitter and a point source diffusing into a spheroidal receiver. These are combined into an end-to-end spheroid-to-spheroid (S2S) channel model, and the resulting communication system's bit error rate (BER) is evaluated under on-off keying (OOK) with inter-symbol interference (ISI).

The case study parameterizes a pair of homogeneous spheroids after real HepaRG liver spheroids, exchanging a signalling molecule across an unbounded, flow-free fluid medium. All analytical results are validated against a custom particle-based simulator (PBS), described in the paper's Section V-A.

## Repository structure
```
main_s2s_ber_analysis.m                 — driver script: end-to-end S2S concentration + BER analysis
src/
├── transmitter_greens_function.m        — GFC, point source inside a spheroidal transmitter (Sec. III-A, App. A-B)
├── receiver_greens_function.m           — GFC, point source outside a spheroidal receiver (Sec. III-B, App. A-C)
├── compute_release_rate.m               — aggregate release rate g(t) of the transmitting spheroid (Eq. 9-10)
├── compute_receiver_response.m          — end-to-end S2S channel response (Eq. 14-16)
└── compute_ber_vs_timeslot.m            — BER vs. time-slot duration (Appendix B, Eq. 46-52)
data/
└── release_rate_Nc24000.mat             — (generated on first run) cached g(t) for Nc_tx = 24000
```

Each `src/` function is self-contained and solves one stage of the boundary-value diffusion problem. The two Green's function scripts build and solve the frequency-domain linear system for their respective spheroid, reconstruct the concentration via a spherical-harmonic expansion, and apply an inverse FFT to obtain the time-domain response. The release-rate and receiver-response scripts discretize each spheroid into a grid of point sources/observation points and integrate over them. The BER script implements the genie-aided decision-feedback detector. Comments throughout reference the corresponding equation numbers in the paper so the code can be read alongside the derivation.

## Requirements
- MATLAB R2021b or later
- Statistics and Machine Learning Toolbox (for `poisscdf` in `compute_ber_vs_timeslot.m`)
- No Communications Toolbox dependency — bit-pattern enumeration for the ISI analysis uses base-MATLAB `bitget`

## Usage
1. Open MATLAB and navigate to the repository root.
2. Run the driver script:
```matlab
main_s2s_ber_analysis
```
3. Output: the script computes (and caches on first run to `data/`) the transmitter release rate, then produces two plots:
   - Expected number of observed molecules over time for several receiving-spheroid porosities plus a non-porous "transparent receiver" baseline
   - Corresponding BER as a function of time-slot duration for each case

Key parameters that can be adjusted at the top of `main_s2s_ber_analysis.m`:

| Parameter | Description |
|---|---|
| `Nc_tx`, `vc_tx`, `R_tx` | Transmitting spheroid cell count, cell volume, and radius |
| `EPS_RX_LIST` | Receiving-spheroid porosities to compare |
| `distance` | Centre-to-centre distance between the two spheroids |
| `kd_tx`, `kd_rx` | First-order degradation rates inside the transmitter/receiver |
| `D` | Bulk (free-medium) diffusion coefficient |
| `N` | Molecules released per cell for bit "1" |
| `Ts_list` | Candidate OOK time-slot durations for the BER sweep |

## Model summary
- **Geometry**: two spheroids of radius `R_tx`/`R_rx`, each an aggregation of `Nc` cells, fixed at a centre-to-centre distance in an unbounded, flow-free fluid medium.
- **Porous-medium model**: each spheroid is homogenized with porosity `ε = 1 − Nc·Vc/V_spheroid` and effective diffusion coefficient `D_eff = (ε/τ)·D`, with tortuosity `τ = ε^−0.5`.
- **Boundary conditions**: continuity of diffusive flux and a porosity-dependent concentration jump (`c_inside = κ·c_outside`, with `κ = sqrt(D/D_eff)`) at each spheroid surface.
- **Solution method**: Green's functions expanded in spherical harmonics and solved via spherical Bessel/Hankel functions in the frequency domain, converted to the time domain via inverse FFT over a swept angular-frequency grid.
- **Channel model**: the transmitter's aggregate release rate g(t) is obtained by integrating its point-source GFC over all constituent cells; the end-to-end response is the volume-averaged convolution of g(t) with the receiver's point-source GFC.
- **Performance evaluation**: OOK with a genie-aided decision-feedback detector under a Poisson channel model, accounting for ISI from previous time slots.

## Citation
```bibtex
@article{rezaei2024spheroidal,
  author={Rezaei, Mitra and Arjmandi, Hamidreza and Zoofaghari, Mohammad and
          Kanebratt, Kajsa and Vilén, Liisa and Janzén, David and
          Gennemark, Peter and Noel, Adam},
  journal={IEEE Transactions on Molecular, Biological, and Multi-Scale Communications},
  title={Spheroidal Molecular Communication via Diffusion: Signaling Between
         Homogeneous Cell Aggregates},
  year={2024},
  volume={10},
  number={1},
  pages={197--210},
  doi={10.1109/TMBMC.2024.3366420}
}
```

## Author
Mitra Rezaei — [LinkedIn](https://www.linkedin.com/in/mitra-rezaei-834784159/) · [ORCID](https://orcid.org/0000-0001-5826-3856) · [Thesis hub](https://github.com/Mitra74/phd-thesis-warwick)
