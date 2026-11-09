# GitHub publishing guide

## Suggested repository name
`mimo-detection-algorithms`

## Short description for the GitHub About box
Implementation and BER comparison of six MIMO detection algorithms in Python and C99: ML, ZF, MMSE, ZF-SIC, MMSE-SIC and lattice-reduction-aided detection.

## Suggested topics
`mimo` `wireless-communications` `digital-communications` `signal-processing` `communication-systems` `python` `c` `qpsk` `16qam` `ber-analysis` `rayleigh-fading` `mmse` `zero-forcing`

## Suggested repository tagline
A comparative study of linear, SIC, ML and lattice-reduction-aided MIMO receivers.

## Upload using Git

After extracting the project folder, open a terminal inside it:

```bash
git init
git add .
git commit -m "Add Python and C99 MIMO detection simulations"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/mimo-detection-algorithms.git
git push -u origin main
```

Replace `YOUR-USERNAME` with your GitHub username and create the empty GitHub repository first.

## Recommended first commit
`Initial release: six MIMO detection algorithms in Python and C99`

## Repository presentation
- Keep `README.md` at the repository root. GitHub displays it automatically.
- Keep the BER comparison image and CSV under `results/`.
- Keep mathematical explanations in `docs/ALGORITHMS.md`.
- Keep test setup and result interpretation in `docs/SIMULATION_METHODOLOGY.md`.
- Add the project report PDF at the repository root or under `docs/` if you want it displayed alongside the code.
