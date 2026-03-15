# Hyperbolic Busemann Evolutionary Algorithm (HBEA)

This repository contains the official MATLAB source code for the **Hyperbolic Busemann Evolutionary Algorithm (HBEA)**.

> **Paper Title:** An Evolutionary Algorithm Based on Busemann Decomposition in Hyperbolic Space for Many-Objective Optimization  

---

## Abstract

Many-objective evolutionary algorithms have become a primary tool for solving many-objective optimization problems. However, as the number of objectives increases, dominance pressure weakens rapidly, Euclidean neighborhood information loses discriminative power, and reference-vector-based guidance may degrade on irregular Pareto fronts. To address these issues, this paper proposes the Hyperbolic Busemann Evolutionary Algorithm (HBEA), which embeds normalized objective vectors into the Poincaré ball, where hyperbolic geometry enhances directional separation, and constructs a Busemann fingerprint to orthogonally decompose each solution into a convergence component and a distribution component. Based on this decomposition, HBEA performs selection and mating with explicitly separated convergence and diversity guidance, without relying on dense predefined reference vectors. Experimental results on representative benchmarks show that HBEA is highly competitive with existing many-objective evolutionary algorithms while maintaining practical computational cost.

## How to Run the Experiment

The HBEA algorithm is implemented to be fully compatible with the **PlatEMO** . To run the algorithm and replicate the results:

1.  Download the [PlatEMO platform](https://github.com/BIMK/PlatEMO) and ensure it is added to your MATLAB path (tested on MATLAB R2024a] and newer).
2.  Copy the `HBEA.m` file from this repository into the PlatEMO multi-objective optimization algorithms directory. 
3.  Open MATLAB and run PlatEMO via the command window to open the GUI:
    ```matlab
    >> platemo
    ```
4.  In the GUI, select **HBEA** from the algorithm list, choose your target benchmark problem, and click `Run`.
---

## Citation
If you use this code or the HBEA algorithm in your research, please cite our paper:

(Note: The full citation and DOI will be updated here once the paper is officially accepted.)

## Contact
For any questions, please contact Lin Yang at 15765057075@163.com.
