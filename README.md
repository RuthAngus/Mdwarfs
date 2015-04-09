# Mdwarfs
M dwarfs are cool

Using ABC to constrain the msini distribution of planets orbiting Mdwarfs.

Selection effects and completeness

Generate a population of planets from Morton et al.
Give these planets a mass, using a Mass-Radius relation and
a random inclination.

Ask whether they would be detected. (Could actually marginalise over these
parameters!)

Build a new planet population.

Compute KS statistic.

First pass - initialise with the Weiss and Marcy relation for Rocky planets.
Initialise with a broad Gaussian for the Massive planet population.

1. Use the real Morton et al. radius distribution to create planets.
2. Use the real Marcy and Weiss M-R relation to initialise.
3. Use realistic Kmin and Pmax limits to establish whether these planets
would be detected. Kmin and Pmax could be free parameters to marginalise
over. Look up Bonfils et al., Howard et al.
