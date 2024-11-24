# Spacetime symmetry mapping in a dissipative nonlinear multimode waveguide amplifier 

This repository contains time-domain simulation codes for studying spacetime symmetry mapping and transverse mode instability (TMI) in multimode waveguide amplifiers with thermo-optical nonlinearity and gain saturation. These tools support the analysis and replication of results presented in C.-W. Chen et al.’s paper: “Output control of dissipative nonlinear multimode amplifiers via spacetime symmetry mapping.”
There are three MATLAB codes:
MWAmp.m: Simulates the multimode waveguide amplifier by iteratively solving coupled optical and thermal equations.
MWAbs.m: Models a multimode waveguide absorber, representing the spacetime-symmetry counterpart of the amplifier.
MWAmp_Prtb.m: Perturbs the steady-state solution of the amplifier to evaluate its threshold for transverse mode instability (TMI).

Software Requirements:
These codes require MATLAB R2017b or later to run as is, but they can run on MATLAB R2011a or later if you replace 'isfile' with 'exist'.

Hardware Requirements:
The codes require approximately 1–2 GB of RAM for typical parameters. However, this estimate scales with the number of modes (M), transverse resolution (Rcore), and simulation steps (J and P). A computer with 4 GB of RAM or more should handle it comfortably. For very fine resolutions or extended simulations, higher memory may be necessary.

Demonstration
This demonstration reproduces Fig. 2 in the paper using spacetime symmetry mapping. Follow the steps below to perform the simulation:
1.	Include the required input file:
   Ensure the file b0B_d5.mat is in the same folder as the multimode waveguide absorber code, MWAbs.m. The variable ‘b0’ in the file represents the mode content of the phase conjugation of a target optical field at the distal end of a multimode waveguide amplifier.
32.	Run the multimode waveguide absorber simulation (MWAbs.m):
   Execute MWAbs.m to calculate light absorption from the distal end of the multimode waveguide absorber to the proximal end.
   Allow the simulation to run until it converges (several hundred milliseconds in the simulation's time frame).
3.	Validate steady-state field and temperature profiles:
   Phase-conjugate the steady-state field at the proximal end of the multimode waveguide absorber using b0 = conj(b(:,end)), and flip the steady-state temperature profile using T = fliplr(T). Save both in b0T.mat.
   Ensure the output file b0T.mat is saved in the same folder as the multimode waveguide amplifier code, MWAmp.m.
   Compare the steady-state field and temperature profiles recorded in b0T.mat to those in b0T_D5.mat, which can be downloaded from the Zenodo repository of the paper (https://zenodo.org/records/14190653). The profiles should closely match.
4.	Run the multimode waveguide amplifier simulation (MWAmp.m):
   Execute MWAmp.m to calculate light amplification from the proximal end of the multimode waveguide amplifier to the distal end.
   The resulting field at the distal end should closely resemble the target output field (the phase conjugation of the field in b0B_d5.mat).
By completing these steps, you should successfully reproduce Fig. 2 in the paper. For any questions or discrepancies, please refer to the paper or reach out for clarification.
To test the stability of the steady-state solutions of the multimode waveguide amplifier simulated in MWAmp.m, save the steady-state field vector (b) as b.mat and the temperature (T) and refractive index (Dn) profiles as DnT.mat. Alternatively, you can download DnT.mat from Zenodo (https://zenodo.org/records/14190653). Once the files are prepared, run MWAmp_Prtb.m to perturb the steady-state solution and evaluate the TMI threshold. More details on this procedure can be found in our previous paper: “Suppressing transverse mode instability through multimode excitation in a fiber amplifier” (https://www.pnas.org/doi/10.1073/pnas.2217735120).

Contact and License
For any questions, issues, or feedback, feel free to contact Chun-Wei Chen:
•	GitHub: https://github.com/joe851642001
•	Email: joec.cms@gmail.com
This project is licensed under the MIT License. You are free to use, modify, and distribute this software as long as proper attribution is given. For more details, please refer to the LICENSE file included in the repository.  
