# Licence of this directory

The files in `validation/croppol/` are products of the CropPol database and are
released under the **Creative Commons Attribution-NonCommercial-ShareAlike 4.0
International licence (CC BY-NC-SA 4.0)**,
<https://creativecommons.org/licenses/by-nc-sa/4.0/>. They are for
non-commercial research use only.

The MIT licence of the pollinatoR package does **not** cover this directory.
The directory sits under `validation/`, which `.Rbuildignore` excludes
(`^validation$`), so none of it enters the built package.

## Why this licence

The CropPol sources state different terms. The GitHub repository carries an
MIT `LICENSE`; the repository metadata (`biblio.csv`, `dataspice.json`, the EML
`intellectualRights` field) say "CC By"; the Zenodo releases are marked "other
open". The published data paper states: "This data set is released for
non-commercial use only. Credits should be given to this paper (i.e., proper
citation), and the products generated with this database should be shared
under the same license terms (CC BY-NC-SA)." This directory follows that
statement, the strictest of the four, until the maintainers say which one
governs the data.

## No copy of the data

No CropPol data are stored here. The script reads the field-level table at run
time from the public GitHub URL, pinned to one commit, and caches it only in the
session's temporary directory.

## Citation

Allen-Perkins, A., Magrach, A., Dainese, M., Garibaldi, L. A., Kleijn, D.,
Rader, R., et al. (2022). CropPol: A dynamic, open and global database on crop
pollination. *Ecology*, *103*(3), e3614. <https://doi.org/10.1002/ecy.3614>

The study reproduced here, which contributed its data to CropPol:

Zou, Y., Xiao, H., Bianchi, F. J. J. A., Jauker, F., Luo, S., & van der Werf, W.
(2017). Wild pollinators enhance oilseed rape yield in small-holder farming
systems in China. *BMC Ecology*, *17*, 6.
<https://doi.org/10.1186/s12898-017-0116-1>
