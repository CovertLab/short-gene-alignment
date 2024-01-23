# Manual alignment of short *E. coli* genes to sequencing data

This repository contains the code used to manually align RNA-Seq reads to
short (<100nt) genes in the *E. coli* genome that are typically misreported
by RNA-Seq algorithms to have zero read counts. The data generated
by the code were used as supporting figures
for Sun et al., "Cross-evaluation of E. coli’s operon structures
via a whole-cell model suggests alternative cellular benefits
for low- versus high-expressing operons" (2023), in review.

You can reach out to ggsun AT stanford.edu for any questions
about the code in this repository.

This repository was developed with Python 3.11.3. After installing Python, you can install
all of the required packages by running

```shell
pip install -r requirements.txt
```

To regenerate all short gene alignment data used in Sun et al. (2023), run

```shell
bash runscripts/align_short_genes.sh
```
