# Physics-Aware Degradation Model of LIBESS for Techno-Economic Studies in Power Systems
This work proposes a new hybrid model of lithium-ion battery energy storage system (LIBESS) that utilizes the linearity of a simple energy reservoir model and follows the rules of the physics-based degradation. This LIBESS model is a blend of the linearized physics-based single particle model coupled with the description of SEI growth and energy-reservoir model. The goal of this case study is to showcase the capabilities of the proposed LIBESS model in optimization.
This GitHub repository contains (1) a code to derive the strategic operation of LIBESS considering the proposed model and (2) a code to simulate lithium-ion cell behaviour using an actual single particle model coupled with the description of solid electrolyte interphase growth. The latter is used to validate the results of optimization.
To run optimization problem, the following files with parameters are used:
-	Alberta’s electricity market pool prices in 2021
-	The system level parameters of LIBESS
-	The cell-level parameters
The link to the paper that uses this code will be updated later.
## Requirements:
The code for both simulation and optimization models is written in the Julia programming environment. To run an optimization problem, a valid Gurobi license is required.
