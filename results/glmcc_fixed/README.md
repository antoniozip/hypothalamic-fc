# DO NOT USE — signs are wrong

These 26 matrices were produced by `vendor/glmcc-c/glmcc_fixed`, which has a confirmed sign
defect: excitatory couplings are reported as inhibitory at ~0.24x the correct magnitude, and
genuine inhibitory couplings are largely missed. See `reports/glmcc_c_sign_bug.md`.

Magnitudes are broadly right (corr with the Python reference = +0.90 on |W|), but the signed
correlation is -0.86. Any E/I claim from these files is wrong, and weighted metrics are distorted.

Regenerate with the vendored Python `Est_Data.py` until the C port passes
`scripts/glmcc_positive_control.py`.
