Demo data and code for "Mechanochemical morphodynamics of diverse axonal growth phenotypes"
===========================================================================
This project contains MATLAB scripts to reproduce the results from a study on temporal evolutions of axonal growth and corresponding signaling, and phase diagram. 

Scripts included:
- run_axongrowth.m: Reproduces Figure 2 results for simulated dynamics for the axonal growth and microtubule-based signaling.
- run_phasediagram_S_tauT.m: Reproduces Figure 3 results for simulated phase diagram in the time delay-feedback strength parameter space.
- run_PhaseDiagram_growth_phenotypes.m: Reproduces Figure 5 results for simulated phase diagram of axonal growth phenotypes.

-------------------------------------------------------
1. System Requirements
-------------------------------------------------------
- Software: MATLAB R2023a
- Operating System: Windows 10/11, macOS 13+, Linux distributions supporting MATLAB 2023a
- Dependencies: No additional toolboxes required beyond standard MATLAB installation
- Non-standard hardware: None (CPU computation, GPU not required)
- Tested on: MATLAB 2023a on Windows 11

-------------------------------------------------------
2. Installation Guide
-------------------------------------------------------
1. Install MATLAB R2023a following standard instructions from MathWorks.
2. Download or clone the project folder containing all MATLAB scripts and dataset `.mat` files.
3. Ensure the folder containing scripts is added to MATLAB path.

-------------------------------------------------------
3. Demo
-------------------------------------------------------
To reproduce the demo figures using the provided datasets:

1. Open MATLAB and set the current folder to the project directory.
2. Run the scripts sequentially:
   a. run_axongrowth.m: Generate time series of axonal length and signaling under normal condition (Fig. 2a,b,e,f,g).
   b. run_phasediagram_S_tauT.m: Generate signal oscillation frequency in the time delay-feedback strength parameter space (Fig. 3b).
   c. run_PhaseDiagram_growth_phenotypes.m: Generate phase diagram illustrating the dependence of axonal growth phenotypes on initial length L0 and feedback strength S, and phase diagram in the feedback strength and time delay plane (Fig. 5d,e).

