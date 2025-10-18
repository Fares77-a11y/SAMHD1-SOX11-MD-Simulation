# ⚛️ SAMHD1–SOX11 Molecular Dynamics (MD) Simulation Analysis

This repository hosts a complete MD simulation workflow and post-analysis of the SAMHD1–SOX11 complex following HADDOCK docking. Two binding modes (Cluster 200: random, Cluster 361: centroid-restrained) were simulated to evaluate structural stability, flexibility, and inter-chain interactions.

---

## 🧪 Goals

- Simulate the dynamics of the SAMHD1–SOX11 complex in GROMACS
- Evaluate stability via RMSD, RMSF, RoG, and SASA
- Quantify protein–protein interaction interface and hydrogen bonds
- Visualize structural dynamics and binding behavior per docking cluster

---

## 🗂️ Repository Structure

```
SAMHD1-SOX11-MD-Simulation/
├── data/                       # Raw GROMACS inputs (topol.top, mdp, gro, ndx)
├── random/
│   ├── topol.top
│   ├── step3_input.gro
│   ├── step4.0_minimization.{gro,mdp,tpr}
│   ├── step4.1_equilibration.{gro,mdp,tpr}
│   ├── step4.2_equilibration.{gro,mdp,tpr}
│   └── step5_production.{gro,mdp,tpr}
│
└── centroid/
    ├── topol.top
    ├── step3_input.gro
    ├── step4.0_minimization.{gro,mdp,tpr}
    ├── step4.1_equilibration.{gro,mdp,tpr}
    ├── step4.2_equilibration.{gro,mdp,tpr}
    └── step5_production.{gro,mdp,tpr}
├── Results/  
├── figures/                  # Trajectories (RMSD, RMSF, ROG, SASA)
│       ├── SAMHD1 CHAINS RMSD.jpg
│       ├── SAMHD1 CHAINS ROG .jpg
│       ├── SAMHD1 CHAINS SASA.jpg
│       ├── SAMHD1 COMBINED RMSD.jpg
│       ├── SAMHD1 COMBINED ROG.jpg
│       ├── SAMHD1 COMBINED SASA.jpg
│       ├── SOX11 RMSD.jpg
│       ├── SOX11 ROG.jpg
│       ├── SOX11 RMSF.jpg
│       ├── SAMHD1A RMSF.jpg
│       ├── SAMHD1B RMSF.jpg
│       ├── SAMHD1C RMSF.jpg
│       ├── SAMHD1D RMSF.jpg
├── VMD-rendered molecular frames/
│       ├── cluster200_frame.png
│       └── cluster361_frame.gif
├── scripts/                    # TCL, bash, and Python analysis scripts
│   ├── run_gromacs_48.sh
│   ├── analyze_samhd1_sox11_vmd_revised_hbonds_v3.tcl
│   ├── unified_md_analysis_protien.tcl
│   ├── unified_md_analysis_protien_chains.tcl
│   └── pp_docking_vis.py
├── README.md
├── LICENSE
└── .gitignore
```

---

## 🔁 Workflow Overview

1. **Preparation & Docking**
   - Input models from HADDOCK cluster 200 and 361
   - Cleaned and aligned for GROMACS compatibility

2. **MD Simulation (GROMACS)**
   - Run on SLURM using GPU node
   - Continuation support (`-cpi`) for checkpoint recovery

```bash
bash scripts/run_gromacs_48.sh
```

3. **Analysis**
   - TCL/VMD scripts for structural metrics: RMSD, RMSF, RoG, SASA, H-bonds
   - Python visualization via `pp_docking_vis.py`

---

## 📈 Key Output Visualizations

### 📊 Cluster 200 (Random Docking)

| Snapshot | Analysis |
|-----------|-----------|
| ![](Results/VMD-rendered%20molecular%20frames/cluster200_frame.gif) | ![](Results/cluster200_frame.gif) |

---

### 📊 Cluster 361 (Centroid Docking)

| Snapshot | Analysis |
|-----------|-----------|
| ![](Results/VMD-rendered%20molecular%20frames/cluster361_frame.gif) | ![](Results/cluster361_frame.gif) |

---

### 🧬 SAMHD1 Chain-Level Dynamics

| RMSD | ROG | SASA |
|------|-----|------|
| ![](Results/figures/SAMHD1%20CHAINS%20RMSD.jpg) | ![](Results/figures/SAMHD1CHAINS%20ROG.jpg) | ![](Results/figures/SAMHD1%20CHAINS%20SASA.jpg) |

**RMSD:** Per-chain deviation comparison across centroid vs random protocols.  
**ROG:** Structural compactness and changes over time.  
**SASA:** Solvent-accessible surface area for individual SAMHD1 chains.

---

### 🧩 SAMHD1 Combined Trajectories Dynamics (Chains A+B+C+D as one object) (Centroid vs Random)

| RMSD | ROG | SASA |
|------|-----|------|
| ![](Results/figures/SAMHD1%20COMBINED%20RMSD.jpg) | ![](Results/figures/SAMHD1%20COMBINED%20ROG.jpg) | ![](Results/figures/SAMHD1%20COMBINED%20SASA.jpg) |

**Combined RMSD:** Overall deviation of the full complex.  
**Combined ROG:** Global compactness across docking strategies.  
**Combined SASA:** Surface exposure changes for the complex as a whole.

#### 📊 Residue-Level SAMHD1 Analysis (Random VS Centroid Docking)
| Metric | Plot |
|--------|------|
| RMSF A | ![](Results/figures/SAMHD1A%20RMSF.jpg) |
| RMSF B | ![](Results/figures/SAMHD1B%20RMSF.jpg) |
| RMSF C | ![](Results/figures/SAMHD1C%20RMSF.jpg) |
| RMSF D | ![](Results/figures/SAMHD1D%20RMSF.jpg) |

- **RMSF**: Per-residue flexibility comparison for centroid vs random placements.

### 🧠 SOX11 (Chain E) (Centroid vs Random)

| Metric | Plot |
|--------|------|
| RMSD | ![](Results/figures/SOX11%20RMSD.jpg) |
| ROG | ![](Results/figures/SOX11%20ROG.jpg) |
| SASA | ![](Results/figures/SOX11%20SASA.jpg) |
| RMSF | ![](Results/figures/SOX11%20RMSF.jpg) |


---

## 🧬 Metrics Captured

- **Root Mean Square Deviation (RMSD)**
- **Radius of Gyration (RoG)**
- **Solvent Accessible Surface Area (SASA)**
- **Hydrogen Bond Count**
- **Interface SASA** (Inter-chain only)
- **Chain-specific RMSF**

---

## 📌 Reproduce Locally

1. Clone the repository
```bash
git clone https://github.com/YOUR_USERNAME/SAMHD1-SOX11-MD-Simulation.git
cd SAMHD1-SOX11-MD-Simulation
```

2. Prepare your GROMACS environment (e.g., GROMACS 2023.4 with GPU)

3. Modify and run `scripts/run_gromacs_48.sh` for your local SLURM config

4. Analyze with:
```bash
vmd -dispdev text -e scripts/unified_md_analysis_protien.tcl
python3 scripts/pp_docking_vis.py
```

---

## 👨‍💻 Author

**Fares Ibrahim**  
Bioinformatics | Structural Biology | MD Simulation  
🔗 [LinkedIn](https://www.linkedin.com) | 🌐 [GitHub](https://github.com/Fares77-a11y)

---
