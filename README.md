# Vina SDF Reconstruction for Molecular Dynamics Simulations

A Python workflow for reconstructing chemically complete ligand structures from AutoDock Vina PDBQT docking poses while preserving the docked heavy atom coordinates.

## Why this is needed

AutoDock Vina PDBQT files do not preserve all chemical information needed for reliable ligand parameterization.

This workflow combines:

- an authoritative SDF structure for ligand chemistry
- an AutoDock Vina PDBQT file for docked coordinates

The SDF supplies the chemical structure, while pose 1 from the PDBQT supplies the heavy atom coordinates.

## Workflow

The program:

1. Reads the authoritative SDF structure.
2. Extracts pose 1 from the Vina PDBQT output.
3. Removes hydrogens from the docking pose.
4. Verifies heavy atom composition.
5. Determines the heavy atom correspondence between the SDF and docked pose.
6. Transfers the docked coordinates to the SDF chemical structure.
7. Verifies preservation of the docked heavy atom coordinates.
8. Adds hydrogens using RDKit.
9. Checks for severe nonbonded hydrogen and heavy atom clashes.
10. Determines the formal molecular charge.
11. Runs ACPYPE using GAFF2 and AM1-BCC with the formal charge supplied explicitly.
12. Verifies that ACPYPE preserved the reconstructed heavy atom coordinates.

The final pose-preservation check requires a maximum heavy atom displacement below 0.0001 Å.

## Requirements

Tested with:

- Python 3.8.20
- RDKit 2024.03.2
- NumPy 1.24.4
- ACPYPE 2023.10.27
- Antechamber 22.0
- GAFF2
- AM1-BCC
- Open Babel 3.1.0
- CentOS Linux 7

See `environment.txt` for the tested software environment.

## Usage

### Direct file mode

For most users:

```bash
python src/sdf_reconstruction.py ligand1 \
    --sdf ligand1.sdf \
    --pdbqt ligand1_out.pdbqt \
    --output-dir output
```

The ligand identifier is used for output file naming.

The PDBQT file may contain multiple Vina models. The workflow extracts and reconstructs the first docking pose.

### Included example

```bash
python src/sdf_reconstruction.py ZINC000247714414 \
    --sdf example/ZINC000247714414.sdf \
    --pdbqt example/ZINC000247714414_Nterm_out.pdbqt \
    --output-dir example/output
```

A successful run ends with:

```text
ACPYPE: PASS
FINAL POSE PRESERVATION: PASS
FULL SDF RESCUE PASS: ZINC000247714414
```

### Mapping mode

For pre-resolved datasets:

```bash
python src/sdf_reconstruction.py ZINC000247714414 \
    --mapping resolved_mapping.tsv \
    --templates zinc_templates \
    --ligands ligands \
    --output-dir output
```

For ordinary use with a known SDF and PDBQT pair, direct file mode is recommended.

## Validation

The workflow was developed using a 75-ligand dataset and subsequently evaluated using an independent 200-ligand dataset.

Compact validation files are included in `validation/`:

- `development75_original_workflow.tsv`
- `validation200_cohort.tsv`
- `VALIDATION_SUMMARY.txt`

## Important limitations

Successful reconstruction and parameterization do not establish that a docking pose is physically correct or that a ligand has biological activity.

The supplied SDF must represent the intended chemical structure corresponding to the ligand in the PDBQT file.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── environment.txt
├── sdf_reconstruction.py
├── example/
│   ├── ZINC000247714414.sdf
│   └── ZINC000247714414_Nterm_out.pdbqt
└── validation/
    ├── development75_original_workflow.tsv
    ├── validation200_cohort.tsv
    └── VALIDATION_SUMMARY.txt
```
