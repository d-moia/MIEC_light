This code evaluates the defect chemistry of a mixed ionic-electronic conductor under dark and under light based on a 0 dimensional model (reaction limited, no transport included). The input parameters are relevant to halide perovskites such as methylammonium lead iodide (MAPI), with anti-Frenkel ionic (iodide) disorder, however it can be adapted to other material systems, by changing the relevant electronic and ionic properties.

The two codes that can be run independently are KroegerVink_DarkLight and KroegerVink_DarkLight_sgeq. The first refers to a full description of the problem where the solid-gas exchange is treated explicitly, while the second uses the sg-eq condition (solid-gas exchange at equilibrium even under light). Each of these codes calls the Equations_DarkLight (or Equations_DarkLight_sgeq) where the equations referring to the steady-state problem are written.

How to use the code:
- Select directory, folder name and calculation name
- Change input parameters, such as Gext, Gamma_I_i, Gamma_I_v, Gamma_p_i, Gamma_n_v, recombination paramteres etc.
- Make sure that the Equations_DarkLight.m (or Equations_DarkLight_sgeq.m) file is stored in the same folder as this script
- Click on Run

What the code gives as output:
The code automatically plots relevant graphs with defect concentrations, QFLS and ionic chemical potential, etc. Additional plots can be defined.
The code also saves txt files with the most important data and an info file with the input parameters used in the calculation.
