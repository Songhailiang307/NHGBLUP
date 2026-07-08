# NHGBLUP

A nonlinear hybrid genomic best linear unbiased prediction framework for improving genomic prediction by integrating nonlinear genomic relationship matrices.

NHGBLUP implements four genomic prediction models:

- GBLUP
- GDBLUP
- NHGBLUP
- NHGDBLUP

The pipeline automatically constructs additive and dominance genomic relationship matrices from PLINK genotype files, generates nonlinear hybrid genomic relationship matrices using nonlinear transformation, and evaluates genomic prediction performance across a series of nonlinear weighting coefficients ($\alpha = 0.1$–$1.0$).



---

# Requirements

- Linux
- Python 3.8 or higher
- R 4.3.2 or higher
- PLINK v1.9
- HiBLUP
- DMU v6

---

# Usage

```bash
./NHGBLUP.sh Phenotype.txt GenotypePrefix rel.txt val.txt
```

---

# Description of Arguments

## 1. Phenotype.txt (PHE)

Phenotype file containing all individuals.

The file contains two columns.

| Column | Description |
|--------|-------------|
| 1 | Individual ID |
| 2 | Observed phenotype (or corrected phenotype) |

### Example

```text
ID  DTA
1   70.82
2   68.34
3   62.24
4   69.17
5   64.63
6   62.80
7   68.69
8   67.80
9   70.35
...
```

---

## 2. GenotypePrefix (SNP)

PLINK binary genotype files.

Three files with the same prefix are required.

```text
GenotypePrefix.bed
GenotypePrefix.bim
GenotypePrefix.fam
```

### Example

```text
pig.bed
pig.bim
pig.fam
```

---

## 3. rel.txt (REL)

Reference population.

The file contains three columns.

| Column | Description |
|--------|-------------|
| 1 | Individual ID |
| 2 | Fixed group indicator (default = 1) |
| 3 | Observed phenotype |

### Example

```text
1   1   70.82
3   1   62.24
4   1   69.17
5   1   64.63
6   1   62.80
8   1   67.80
9   1   70.35
...
```

Individuals in this file are used as the reference population for model training.

---

## 4. val.txt (VAL)

Validation population.

The file contains three columns.

| Column | Description |
|--------|-------------|
| 1 | Individual ID |
| 2 | Fixed group indicator (default = 1) |
| 3 | Observed phenotype |

### Example

```text
2    1   68.34
7    1   68.69
10   1   70.58
11   1   68.84
17   1   67.23
24   1   63.51
28   1   73.69
37   1   69.34
38   1   69.20
55   1   64.25
60   1   67.39
62   1   76.14
64   1   70.15
65   1   72.33
72   1   60.70
...
```

Individuals in this file are used only for genomic prediction evaluation.

---

# Notes

- Individual IDs must be identical across all input files.
- The order of individuals in the PLINK genotype files (`.bed`, `.bim`, `.fam`) should be consistent with the phenotype file.
- The reference population (REL) and validation population (VAL) must not overlap.
- The union of REL and VAL should correspond to all individuals contained in the phenotype and genotype files.

---

# Output Files

### GBLUP.txt

Prediction performance of GBLUP.

Output metrics:

- `corr`
- `reg`
- `MSE`
- `MAE`

---

### GDBLUP.txt

Prediction performance of GDBLUP.

Output metrics:

- `corr`
- `reg`
- `MSE`
- `MAE`

---

## NHGBLUP

The nonlinear weighting coefficient $\alpha$ ranges from **0.1** to **1.0**.

For each $\alpha$, the following files are generated.

### GBLUP_α.txt

Prediction performance of NHGBLUP.

Example:

```text
GBLUP_0.1.txt
GBLUP_0.2.txt
...
GBLUP_1.0.txt
```

---

### D2GBLUP_α.txt

Prediction performance of the nonlinear hybrid model using only the nonlinear genomic relationship matrix.

Example:

```text
D2GBLUP_0.1.txt
...
D2GBLUP_1.0.txt
```

---

### D2GDBLUP_α.txt

Prediction performance of NHGDBLUP.

Example:

```text
D2GDBLUP_0.1.txt
...
D2GDBLUP_1.0.txt
```

---

# Data Simulation

Phenotypic and genotypic data with both additive and dominance genetic effects can be simulated using:

```bash
Rscript ./bin/simFInaAD.R h2 delta
```

### Example

```bash
Rscript ./bin/simFInaAD.R 0.3 1
```

where

- `0.3` represents the simulated heritability ($h^2 = 0.3$).
- `1` represents the dominance degree ($\delta = 1$), corresponding to complete dominance.

The simulation script generates phenotype and genotype datasets containing both additive and dominance genetic effects, which can be directly used as inputs for:

- NHGBLUP
- NHGDBLUP
- GBLUP
- GDBLUP

---

# Script Description

The `bin` folder contains the following programs.

| Program | Description |
|---------|-------------|
| `COR_REG_used2.py` | Calculates prediction accuracy (`corr`), unbiasedness (`reg`), MSE, and MAE. |
| `dBLUP.py` | Constructs nonlinear hybrid genomic relationship matrices using nonlinear transformation (currently **tanh**) and performs NHGBLUP. |
| `GINV.py` | Computes the inverse of genomic relationship matrices. |
| `Gma_to_3lineID.py` | Converts genomic relationship matrices into the three-column sparse format required by DMU. |
| `gebv.py` | Extracts genomic estimated breeding values (GEBVs) from GBLUP solutions. |
| `ggebv.py` | Extracts additive and dominance GEBVs from GDBLUP and NHGDBLUP solutions. |
| `gblup.DIR` | DMU driver file for GBLUP. |
| `gdblup.DIR` | DMU driver file for GDBLUP. |
| `g2dblup.DIR` | DMU driver file for NHGDBLUP. |
| `dmu1` | DMU executable. |
| `dmuai` | DMU variance component estimation module. |
| `r_dmuai` | Shell wrapper for running DMU. |
| `hiblup` | HiBLUP executable for constructing additive and dominance genomic relationship matrices. |
| `plink` | PLINK v1.9 executable for genotype processing. |
| `simFInaAD.R` | R script for simulating additive and dominance genetic effects. |

---

# Model Overview

The implemented models include:

| Model | Description |
|-------|-------------|
| **GBLUP** | Genomic Best Linear Unbiased Prediction |
| **GDBLUP** | Genomic Best Linear Unbiased Prediction including dominance effects |
| **NHGBLUP** | Nonlinear Hybrid Genomic Best Linear Unbiased Prediction |
| **NHGDBLUP** | Nonlinear Hybrid Genomic Best Linear Unbiased Prediction including dominance effects |

NHGBLUP introduces a nonlinear transformation of marker genotypes before constructing the genomic relationship matrix, enabling the model to capture nonlinear genetic signals while remaining computationally efficient within the conventional mixed-model framework.

---

# License

This project is licensed under the **GPL-3.0 License**.
