# Environment and software versions

## Main analysis environment

- Python 3.12.10
- RDKit 2026.09.1
- meeko 0.8.0
- Open Babel 3.1.0
- numpy, pandas, scipy, matplotlib (versions recorded in the project environment)
- AutoDock Vina 1.2.7 (external executable)

## MD / MM-GBSA environment

- OpenMM 8.6.1
- OpenFF Sage 2.2.0 (NAGL partial charges)
- openmmforcefields 0.15.1
- ParmEd
- mdtraj 1.11.2
- MDAnalysis 2.11.0
- Python 3.12

The MD environment was kept separate from the main analysis environment to
avoid force-field and OpenMM dependency conflicts.

## Minimal `requirements.txt`

```text
numpy
pandas
scipy
matplotlib
rdkit
meeko==0.8.0
mdtraj==1.11.2
MDAnalysis==2.11.0
openmm==8.6.1
openmmforcefields==0.15.1
```

OpenFF Sage 2.2.0, ParmEd, and the Open Babel executable are installed
separately. AutoDock Vina 1.2.7 is a standalone binary and is not a Python
package.
