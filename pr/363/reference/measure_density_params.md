# Shared parameters for density scores

Shared parameters for density scores

## Arguments

- ylp:

  A numeric matrix or three-dimensional array of log predictive
  densities or probabilities:

  - **Matrix** (\\S \times n\\): \\S\\ posterior draws (chains merged)
    by \\n\\ observations.

  - **Array** (\\I \times C \times n\\): \\I\\ MCMC iterations per
    chain, \\C\\ chains, and \\n\\ observations. Converted to an \\S
    \times n\\ matrix internally.
