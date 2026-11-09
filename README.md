# MIMO Detection Algorithms — Python and C Simulations

**Six MIMO detection algorithms implemented in Python 3 and C99**  
**Authors:** Pranthi (123EC0048) · Chethan (123EC0049)

This repository demonstrates and compares receiver-side detection methods for a spatial-multiplexing MIMO link. It contains standalone console programs in both Python and C for:

1. Maximum-Likelihood (ML)
2. Zero-Forcing (ZF)
3. Minimum Mean-Square Error (MMSE)
4. Zero-Forcing Successive Interference Cancellation (ZF-SIC)
5. MMSE Successive Interference Cancellation (MMSE-SIC)
6. Lattice-Reduction-Aided (LR) detection

The simulations report **Bit Error Rate (BER) versus Signal-to-Noise Ratio (SNR)** and print an example transmitted vector, received vector and detected vector.

---

## 1. Project objective

The objective is to study how different receiver algorithms recover multiple simultaneously transmitted data streams in a MIMO wireless channel. The project compares the trade-off between:

- Detection accuracy / BER
- Computational complexity
- Noise enhancement
- Inter-stream interference
- Error propagation
- Implementation practicality

No single detector is best for every condition. The choice depends on channel quality, SNR, antenna dimensions, modulation order and available processing resources.

## 2. MIMO System Model

For a narrowband, flat-fading spatial-multiplexing MIMO system, the received signal is represented by:

$$
\mathbf{y} = \mathbf{H}\mathbf{x} + \mathbf{n}
$$

where:

| Symbol | Dimension | Meaning |
|:---:|:---:|---|
| $\mathbf{x}$ | $N_t \times 1$ | Transmitted symbol vector |
| $\mathbf{y}$ | $N_r \times 1$ | Received signal vector |
| $\mathbf{H}$ | $N_r \times N_t$ | MIMO channel matrix |
| $\mathbf{n}$ | $N_r \times 1$ | Additive White Gaussian Noise (AWGN) vector |
| $N_t$ | Scalar | Number of transmit antennas |
| $N_r$ | Scalar | Number of receive antennas |

The receiver uses the received vector $\mathbf{y}$ and channel matrix $\mathbf{H}$ to estimate the transmitted symbol vector $\mathbf{x}$.
## 3. Algorithms at a glance

| Detector | Main idea | Main advantage | Main limitation |
|---|---|---|---|
| ML | Test all possible symbol vectors and choose minimum Euclidean distance | Best detection under the assumed model | Exponential search complexity |
| ZF | Apply the channel pseudo-inverse | Simple linear interference suppression | Can amplify noise |
| MMSE | Regularised linear filtering using noise variance | Balances interference and noise | Residual inter-stream interference |
| ZF-SIC | ZF detection, slicing and sequential cancellation | Better separation than plain ZF in many cases | Error propagation |
| MMSE-SIC | MMSE filtering plus sequential cancellation | Strong practical performance/complexity balance | More processing; error propagation |
| LR | Reduce the channel basis, then apply MMSE-aided detection | Can improve conditioning for linear detection | Pre-processing overhead and finite-constellation caveats |

Detailed explanations, equations, use cases and trade-offs are in [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md).

## 4. Repository structure

```text
MIMO-Detection-Algorithms/
├── README.md
├── LICENSE
├── .gitignore
├── src/
│   ├── python/
│   │   ├── ML_MIMO_Detection.py
│   │   ├── ZF_MIMO_Detection.py
│   │   ├── MMSE_MIMO_Detection.py
│   │   ├── ZF_SIC_MIMO_Detection.py
│   │   ├── MMSE_SIC_MIMO_Detection.py
│   │   └── LR_MIMO_Detection.py
│   └── c/
│       ├── ML_MIMO_Detection.c
│       ├── ZF_MIMO_Detection.c
│       ├── MMSE_MIMO_Detection.c
│       ├── ZF_SIC_MIMO_Detection.c
│       ├── MMSE_SIC_MIMO_Detection.c
│       └── LR_MIMO_Detection.c
├── docs/
│   ├── ALGORITHMS.md
│   └── SIMULATION_METHODOLOGY.md
└── results/
    ├── ber_comparison.csv
    ├── ber_comparison.png
    ├── *_python_output.txt
    └── *_c_smoke_test.txt
```

## 5. Requirements

### Python
- Python 3.x
- Python standard library only; no third-party package is required to run the detector scripts.

### C
- C compiler supporting C99 (for example GCC)
- Standard C math and complex libraries

## 6. Run the Python simulations

Open a terminal in the repository root and run one detector:

```bash
python src/python/ML_MIMO_Detection.py
```

Replace `ML_MIMO_Detection.py` with any of the other five Python filenames.

The program prompts for:
1. Transmit antenna count \(N_t\) (1–4)
2. Receive antenna count \(N_r\) (1–4)
3. Modulation: QPSK or 16-QAM
4. Channel: Rayleigh, AWGN or Rician
5. Selected SNR in dB (-5 to 40)
6. Number of frames (1–2000)
7. Symbols per frame (10–1000)
8. LLL delta (LR only; 0.51–0.99)

**Suggested first run:** enter `2`, `2`, `1` (QPSK), `1` (Rayleigh), `20`, `100`, `50`. For LR, enter `0.75` for the additional LLL delta prompt.

The program sweeps BER at 0, 5, 10, 15, 20, 25 and 30 dB. If the selected SNR is outside this range, it is added as an additional point. The printed example vectors correspond to the final SNR point in the sweep.

## 7. Compile and run the C simulations

Example for ML:

```bash
gcc src/c/ML_MIMO_Detection.c -o ml_mimo -std=c99 -lm
./ml_mimo
```

For another detector, change the source filename and executable name. On Windows, run the generated `.exe` file:

```powershell
gcc src/c/ML_MIMO_Detection.c -o ml_mimo.exe -std=c99 -lm
.\ml_mimo.exe
```

The same interactive configuration prompts are used.

## 8. Simulation options and parameter meanings

| Parameter | Meaning | Available values |
|---|---|---|
| \(N_t\) | Transmit antennas / spatial streams | 1–4 |
| \(N_r\) | Receive antennas | 1–4 |
| Modulation | Constellation used to map bits to complex symbols | QPSK, 16-QAM |
| Channel | Propagation/noise model | AWGN, Rayleigh, Rician |
| SNR | Signal-to-noise ratio | -5 to 40 dB input |
| Frames | Monte-Carlo frames per SNR | 1–2000 |
| Symbols/frame | Channel uses per frame | 10–1000 |
| LLL \(\delta\) | Lovász parameter for LR reduction | 0.51–0.99 |

**Channel models**
- **AWGN:** identity-like channel in the code plus additive complex Gaussian noise; it is a simple reference channel.
- **Rayleigh:** independent complex Gaussian channel coefficients, modelling rich scattering without a dominant line-of-sight component.
- **Rician:** channel with a deterministic line-of-sight component and random scattered component; the code uses \(K=5\).

**Modulation**
- **QPSK:** 4 complex symbols, 2 bits per symbol.
- **16-QAM:** 16 complex symbols, 4 bits per symbol.

## 9. BER results included in this repository

The `results/` folder contains outputs generated with this common test setup:

- \(N_t=N_r=2\)
- QPSK
- Rayleigh fading
- 100 frames × 50 symbols per frame
- SNR sweep: 0 to 30 dB in 5 dB steps
- Fixed pseudo-random seed in each program

The BER table and plot are included in [`results/`](results/). These are **finite Monte-Carlo estimates**, not theoretical curves. A BER of zero means no bit errors occurred in the finite tested sample; it does not mean the theoretical BER is exactly zero.

See [`docs/SIMULATION_METHODOLOGY.md`](docs/SIMULATION_METHODOLOGY.md) for interpretation and reproducibility notes.

## 10. Important implementation notes and limitations

1. These are **console-based numerical simulations**, not MATLAB GUI applications.
2. The Python implementations use only the standard library. The C programs are C99 console programs.
3. Conventional ZF and this ZF-SIC implementation require \(N_r \ge N_t\) and a full-column-rank channel for the inverse/pseudo-inverse operation.
4. MMSE regularisation can make the linear system solvable in more cases, but it does not create additional spatial degrees of freedom.
5. ML is exhaustive. The code places a limit of 100,000 candidate vectors to prevent impractically large searches.
6. LR uses complex LLL-style Gaussian-integer column operations followed by MMSE-aided linear detection and slicing. The integer basis transformation does not generally preserve a finite QPSK/16-QAM alphabet. Therefore, treat this LR implementation as an educational LR-aided detector, not as an exact finite-alphabet ML detector or a guarantee of lower BER for every run.
7. Results depend on the number of frames/symbols and random channel/noise realisations. Increase the Monte-Carlo sample count for more stable estimates.
8. For a fair comparison, keep antenna count, modulation, channel, frame count, symbols per frame and SNR sweep identical across detectors.

## 11. Recommended learning sequence

1. Start with the MIMO signal model and the ZF baseline.
2. Compare ZF with MMSE to understand noise enhancement and regularisation.
3. Study ZF-SIC and MMSE-SIC to understand detection ordering and error propagation.
4. Use ML as a high-performance reference for small systems.
5. Study LR to understand channel-basis conditioning and lattice reduction.

## 12. Authors

**Pranthi** — 123EC0048  
**Chethan** — 123EC0049

Academic project: MIMO Detection Algorithms — Python and C Simulations.

---

If you use or adapt this project, please retain attribution to the authors and cite the repository version you used.
