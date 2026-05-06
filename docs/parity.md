# Parity And Regression Checks

The parity suite protects the source-level integrated engines against accidental
mathematical changes.

Run the full suite with:

```bash
make parity-audit
```

It currently checks:

- Transfer-matrix connective constants for all four modes on small lattices.
- Monte Carlo sampler eigenvalues, rejection counts, output paths, and SHA-256
  hashes for all four modes.
- CreatorAll exact-span counts, output SHA-256 hashes, and UofS validation for
  all four modes.

The critical 2SAP Monte Carlo benchmark is:

```text
mode 2, L=2, M=1, span=2
dominant eigenvalue = 9.455960990693537
```

Useful subsets:

```bash
python3 scripts/parity_audit.py --no-build --tm-only
python3 scripts/parity_audit.py --no-build --mc-only
python3 scripts/parity_audit.py --no-build --creator-only
```

Fast local checks:

```bash
make quick-check
```

Full local checks:

```bash
make check
```
