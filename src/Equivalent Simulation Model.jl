#=
The Equivalent Simulation Model is essentially a digital twin of the cell built based
on a full representation of the cell through equations for the Single Particle Model and
SEI description as a degradation mechanism.
This code is written in Julia and can be used for lithium-ion battery simulation.
In particular, this code was used to calculate the observed degradation from the optimal
schedule obtained by solving the economic energy arbitrage problem with LIBESS.
This code also can generate the evolution of the SPM variables if the initial SoC and power schedule
are provided.
=#

using DataFrames
using XLSX

directory = "D:\\" #Path to working directory
directory_for_parameters = "D:\\" # Path to the file with the Li-ion cell parameters
file_with_SPM_parameters = directory_for_parameters*"Parameters for SPM - LG Chen2020Kane2022.xlsx"
file_with_power = directory*"Input_Power.xlsx"
SPM_Parameters = DataFrame(XLSX.readtable(file_with_SPM_parameters, "Parameters")...)

#####################################################################################
#Positive and negative OCP
OCP_coefficients_p_df = DataFrame(XLSX.readtable(file_with_SPM_parameters, "Positive Electrode")...)
OCP_coefficients_n_df = DataFrame(XLSX.readtable(file_with_SPM_parameters, "Negative Electrode")...)
OCP_coefficients_p = Matrix(OCP_coefficients_p_df)
OCP_coefficients_n = Matrix(OCP_coefficients_n_df)
spl_OCP_p = Spline1D(OCP_coefficients_p[:,1], OCP_coefficients_p[:,2], k=1) # Spline interpolation of OCP
spl_OCP_n = Spline1D(OCP_coefficients_n[:,1], OCP_coefficients_n[:,2], k=1) # Spline interpolation of OCP

function OCP_p_spline(theta)
    evaluate(spl_OCP_p, theta)
end

function OCP_n_spline(theta)
    evaluate(spl_OCP_n, theta)
end


######################################################################################

file_output = "SPM_SEI_simulation_"

######################################################################################



#Get parameters from a dataframe with the Li-ion cell parameters
R_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Particle radius","Positive electrode"][1]
R_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Particle radius","Negative electrode"][1]
El_thickness_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode thickness","Positive electrode"][1]
El_thickness_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode thickness","Negative electrode"][1]
El_length_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode length","Positive electrode"][1]
El_length_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode length","Negative electrode"][1]
El_width_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode width","Positive electrode"][1]
El_width_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode width","Negative electrode"][1]
eps_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Active material volume fraction","Positive electrode"][1]
eps_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Active material volume fraction","Negative electrode"][1]
k_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Reaction rate","Positive electrode"][1]
k_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Reaction rate","Negative electrode"][1]
alpha_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Charge transfer coefficient","Positive electrode"][1]
alpha_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Charge transfer coefficient","Negative electrode"][1]
D_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Diffusion coefficient","Positive electrode"][1]
D_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Diffusion coefficient","Negative electrode"][1]
c_max_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Maximum concentration","Positive electrode"][1]
c_max_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Maximum concentration","Negative electrode"][1]
c_el = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Li concentration in  Electrolyte","Electrolyte"][1]
theta_soc0_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Stoichiometry at 0% SoC","Positive electrode"][1]
theta_soc0_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Stoichiometry at 0% SoC","Negative electrode"][1]
theta_soc1_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Stoichiometry at 100% SoC","Positive electrode"][1]
theta_soc1_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Stoichiometry at 100% SoC","Negative electrode"][1]
V_min = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Minimum Voltage","Nameplate data"][1]
V_max = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Maximum Voltage","Nameplate data"][1]
Q = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Cell Charge Capacity","Nameplate data"][1]
j0_p_Const = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Exchange current density of the reaction","Positive electrode"][1]
j0_n_Const = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Exchange current density of the reaction","Negative electrode"][1]
Q_E = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Nominal energy capacity","Nameplate data"][1]


###   SEI parameters
R_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI resistivity","SEI"][1]
δ00_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Initial SEI thickness","SEI"][1]
c_bulk_EC = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Concentration of EC in bulk electrolyte","SEI"][1]
k_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Kinetic rate constant for SEI reaction","SEI"][1]
α_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI charge transfer coefficient","SEI"][1]
ρ_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI density","SEI"][1]
M_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI molar weight","SEI"][1]
D_EC = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI layer diffusivity","SEI"][1]
###

#System level parameters
Ch_max = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Maximum charging power","Value"][1]
Dis_max = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Maximum discharging power","Value"][1]
E_nom = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Energy capacity","Value"][1]
η_RTE = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="RTE","Value"][1]/100
SoE_min = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Minimum SoC","Value"][1]*E_nom/100
SoE_max = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Maximum SoC","Value"][1]*E_nom/100
N_eol = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Battery cycle life: Number of cycles","Value"][1]
N_eol_DoD = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Battery cycle life: Specific DoD","Value"][1]/100
Q_eol = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="EOL capacity","Value"][1]/100
T_calendar = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Calendar life","Value"][1]*365*24
C_instal = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="LIBESS capital cost","Value"][1]
SoE_init = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Initial SoC","Value"][1]/100*E_nom
SoE_final = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Final SoC","Value"][1]/100*E_nom



#Thermodynamics Constants
F = 96485.309
R = 8.31451
T = 293


########################
#Input and initial conditions
N = E_nom*1e6/Q_E
file_name_from_optimization = directory*"DISPATCH_NoDegradation_year2021_60MW.xlsx"
Dispatch = DataFrame(XLSX.readtable(file_name_from_optimization, "RESULTS")...)
Dispatch_day = Dispatch[:,:x2] #Charging/Discharging schedule for one day
PowerInput_df = Dispatch1[1:12:(12*24),:]/N*1e6
PowerInput = Matrix(PowerInput_df) #Input to the model is power commands
c_initial_p = (theta_soc0_p*(1-SoC_init)+theta_soc1_p*SoC_init)* c_max_p #Initial concentration of lithium ion in the positive electrode
c_initial_n = (theta_soc0_n*(1-SoC_init)+theta_soc1_n*SoC_init)*c_max_n #Initial concentration of lithium ion in the negative electrode
########################

T_h = size(PowerInput)[1] #24 hours
dt = 0.1 #Time step
N_dr = 40 # Defining a mesh for finite-difference representations of PDE
N_time = trunc(Int,(3600/dt))*T_h  #Number of time steps


SPM = zeros(N_time+1,30) # All the variables
SPM_c_p = zeros(N_time+1,N_dr+1) # for diffusion equation
SPM_c_n = zeros(N_time+1,N_dr+1) # for diffusion equation
dr_p = R_p/N_dr # Spatial step of mesh - positive electrode
dr_n = R_n/N_dr # Spatial step of mesh - negative electrode



SPM[2:size(SPM)[1],1] = LinRange(dt/3600, T_h, N_time) #Time
SPM[1,1] = 0 # Origin

#Power as an input
for i in 1:(T_h)
    SPM[2+(i-1)*trunc(Int,(3600/dt)):(1+i*trunc(Int,(3600/dt))),2] = PowerInput[i,1]*ones(1,trunc(Int,3600/dt)) #Power
end

#################################################################
#Meshing for PDE (diffusion equation)
PDE_p = zeros(N_dr+1, N_dr+1)
PDE_p[1,1] = dt*D_p/dr_p^2*(dr_p^2/dt/D_p-2)
PDE_p[1,2] = dt*D_p/dr_p^2*2
for i in 2:N_dr
    PDE_p[i,i-1] = dt*D_p/dr_p^2
    PDE_p[i,i] = dt*D_p/dr_p^2*(dr_p^2/dt/D_p-2/i-2)
    PDE_p[i,i+1] = dt*D_p/dr_p^2*(2/i+1)
end
PDE_p[N_dr+1,N_dr] = 2*dt*D_p/dr_p^2
PDE_p[N_dr+1,N_dr+1] = 1-2*dt*D_p/dr_p^2


PDE_n = zeros(N_dr+1, N_dr+1)
PDE_n[1,1] = dt*D_n/dr_n^2*(dr_n^2/dt/D_n-2)
PDE_n[1,2] = dt*D_n/dr_n^2*2
for i in 2:N_dr
    PDE_n[i,i-1] = dt*D_n/dr_n^2
    PDE_n[i,i] = dt*D_n/dr_n^2*(dr_n^2/dt/D_n-2/i-2)
    PDE_n[i,i+1] = dt*D_n/dr_n^2*(2/i+1)
end
PDE_n[N_dr+1,N_dr] = 2*dt*D_n/dr_n^2
PDE_n[N_dr+1,N_dr+1] = 1-2*dt*D_n/dr_n^2
##########################################################


#Initial Conditions at time 0:
SPM_c_p[1,1:N_dr+1] = c_initial_p * ones(1,N_dr+1) #Lithium ion concentration in the positive electrode
SPM_c_n[1,1:N_dr+1] = c_initial_n * ones(1,N_dr+1) #Lithium ion concentration in the negative electrode
SPM[1,3] = 0 # Current through the cell
SPM[1,4] = 0 # Molar flux of lithium ions, positive electrode-electrolyte interface
SPM[1,5] = c_initial_p  #Average Lithium ion concentration in the positive electrode
SPM[1,6] = c_initial_p  #Surface Lithium ion concentration in the positive electrode
SPM[1,7] = c_initial_p/c_max_p # Dimensionless surface concentration of lithium ion, positive electrode
SPM[1,8] = 0 # Activation overpotential, positive electrode
SPM[1,9] = OCP_p_spline(SPM[1,7]) #OCP, positive electrode
SPM[1,10] = SPM[1,8]+SPM[1,9] #Solid-phase potential, positive electrode
SPM[1,11] = 0 # # Molar flux of lithium ions, negative electrode-electrolyte interface
SPM[1,12] = c_avg_initial_n #Average Lithium ion concentration in the negative electrode
SPM[1,13] = c_avg_initial_n #Surface Lithium ion concentration in the negative electrode
SPM[1,14] = c_avg_initial_n/c_max_n # Dimensionless surface concentration of lithium ion, negative electrode
SPM[1,15] = 0 # Activation overpotential, negative electrode
SPM[1,16] = OCP_n_spline(SPM[1,14]) #OCP, negative electrode
SPM[1,17] = SPM[1,15]+SPM[1,16] #Solid-phase potential, negative electrode
SPM[1,18] = SPM[1,10]-SPM[1,17] # Voltage
SPM[1,3] = SPM[1,2]/SPM[1,18] # Current through the cell (the first approximation)
SPM[1,20] = 0 #overpotential of the SEIreaction
SPM[1,21] = δ0_sei  #Thickness of SEI layer
SPM[1,22] = c_bulk_EC # Concentration of EC on the surface of the negative electrode
SPM[1,23] = 0 #SEI current density
SPM[1,25] = 0 #Loss of lithium inventory

for k in 2:size(SPM)[1]
    SPM[k,4] = 0
    SPM_c_p[k,:] = PDE_p*SPM_c_p[k-1,:]
    SPM_c_p[k,N_dr+1] = SPM_c_p[k,N_dr+1] + (-SPM[k,4]/D_p)*(2*dt*D_p/dr_p+2*dt*D_p/R_p)
    SPM[k,6] = SPM_c_p[k,N_dr+1]
    SPM[k,7] = SPM[k,6]/c_max_p
    SPM_c_n[k,:] = PDE_n*SPM_c_n[k-1,:]
    SPM_c_n[k,N_dr+1] = SPM_c_n[k,N_dr+1] + (-SPM[k,11]/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
    SPM[k,13] = SPM_c_n[k,N_dr+1]
    SPM[k,14] = SPM[k,13]/c_max_n
    SPM[k,8] = 2*R*T/F*asinh(SPM[k,4]*F/(k_p*SPM[k,7]^alpha_p*c_el^(1-alpha_p)*c_max_p*(1-SPM[k,7])^(1-alpha_p))/2)
    SPM[k,9] = OCP_p_spline(SPM[k,7])
    SPM[k,10] = SPM[k,8]+SPM[k,9]
    SPM[k,11] = 0
    SPM[k,15] = 2*R*T/F*asinh(SPM[k,11]*F/(k_n*SPM[k,14]^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-SPM[k,14])^(1-alpha_n))/2)
    SPM[k,16] = OCP_n_spline(SPM[k,14])
    SPM[k,17] = SPM[k,15]+SPM[k,16]
    SPM[k,18] = SPM[k,10]-SPM[k,17]
    SPM[k,3] = SPM[k,2]/SPM[k,18]
    for kk in 1:20
        SPM[k,4] = -SPM[k,3]/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p)
        SPM_c_p[k,:] = PDE_p*SPM_c_p[k-1,:]
        SPM_c_p[k,N_dr+1] = SPM_c_p[k,N_dr+1] + (-SPM[k,4]/D_p)*(2*dt*D_p/dr_p+2*dt*D_p/R_p)
        SPM[k,6] = SPM_c_p[k,N_dr+1]
        SPM[k,7] = SPM[k,6]/c_max_p
        SPM[k,8] = 2*R*T/F*asinh(SPM[k,4]*F/(k_p*SPM[k,7]^alpha_p*c_el^(1-alpha_p)*c_max_p*(1-SPM[k,7])^(1-alpha_p))/2)
        SPM[k,9] = OCP_p_spline(SPM[k,7])
        SPM[k,10] = SPM[k,8]+SPM[k,9]
        SPM[k,11] = SPM[k,3]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n)
        SPM_c_n[k,:] = PDE_n*SPM_c_n[k-1,:]
        SPM_c_n[k,N_dr+1] = SPM_c_n[k,N_dr+1] + (-SPM[k,11]/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
        SPM[k,13] = SPM_c_n[k,N_dr+1]
        SPM[k,14] = SPM[k,13]/c_max_n
        SPM[k,15] = 2*R*T/F*asinh(SPM[k,11]*F/(k_n*SPM[k,14]^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-SPM[k,14])^(1-alpha_n))/2)
        SPM[k,16] = OCP_n_spline(SPM[k,14])
        SPM[k,17] = SPM[k,15]+SPM[k,16]
        SPM[k,18] = SPM[k,10]-SPM[k,17]
        SPM[k,3] = SPM[k,2]/SPM[k,18]
    end
    SPM[k,20] = SPM[k,17] - SPM[k,11]*F*SPM[k-1,21]*R_sei
    SPM[k,23] = -F*k_sei*SPM[k-1,22]*exp(-α_sei*F*SPM[k,20] /R/T)
    SPM[k,21] = SPM[k-1,21]-dt*SPM[k,23]*M_sei/ρ_sei/(2*F)
    SPM[k,22] = SPM[k,21]/D_EC*(SPM[k,23]/F)+c_bulk_EC
    SPM[k,11] = SPM[k,3]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n)-SPM[k,23]/F
    SPM_c_n[k,:] = PDE_n*SPM_c_n[k-1,:]
    SPM_c_n[k,N_dr+1] = SPM_c_n[k,N_dr+1] + (-SPM[k,11]/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
    SPM[k,13] = SPM_c_n[k,N_dr+1]
    SPM[k,14] = SPM[k,13]/c_max_n
    SPM[k,15] = 2*R*T/F*asinh(SPM[k,11]*F/(k_n*SPM[k,14]^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-SPM[k,14])^(1-alpha_n))/2)
    SPM[k,16] = OCP_n_spline(SPM[k,14])
    SPM[k,17] = SPM[k,15]+SPM[k,16]+SPM[k,11]*F*SPM[k,21]*R_sei
    SPM[k,18] = SPM[k,10]-SPM[k,17]
    SPM[k,3] = SPM[k,2]/SPM[k,18]
    for kk in 1:20
        SPM[k,4] = -SPM[k,3]/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p)
        SPM_c_p[k,:] = PDE_p*SPM_c_p[k-1,:]
        SPM_c_p[k,N_dr+1] = SPM_c_p[k,N_dr+1] + (-SPM[k,4]/D_p)*(2*dt*D_p/dr_p+2*dt*D_p/R_p)
        SPM[k,5] = (4/3*pi*dr_p^3*(SPM_c_p[k,1]+SPM_c_p[k,2])/2+sum(4/3*pi*dr_p^3*(i^3-(i-1)^3)*(SPM_c_p[k,i]+SPM_c_p[k,i+1])/2 for i in 2:N_dr))/(4/3*pi*R_p^3)
        SPM[k,6] = SPM_c_p[k,N_dr+1]
        SPM[k,7] = SPM[k,6]/c_max_p
        SPM[k,8] = 2*R*T/F*asinh(SPM[k,4]*F/(k_p*SPM[k,7]^alpha_p*c_el^(1-alpha_p)*c_max_p*(1-SPM[k,7])^(1-alpha_p))/2)
        SPM[k,9] = OCP_p_spline(SPM[k,7])
        SPM[k,10] = SPM[k,8]+SPM[k,9]
        SPM[k,20] = SPM[k,17] - SPM[k,11]*F*SPM[k-1,21]*R_sei
        SPM[k,23] = -F*k_sei*SPM[k-1,22]*exp(-α_sei*F*SPM[k,20] /R/T)
        SPM[k,21] = SPM[k-1,21]-dt*SPM[k,23]*M_sei/ρ_sei/(2*F)
        SPM[k,22] = SPM[k,21]/D_EC*(SPM[k,23]/F)+c_bulk_EC
        SPM[k,11] = SPM[k,3]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n)-SPM[k,23]/F
        SPM_c_n[k,:] = PDE_n*SPM_c_n[k-1,:]
        SPM_c_n[k,N_dr+1] = SPM_c_n[k,N_dr+1] + (-SPM[k,11]/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
        SPM[k,13] = SPM_c_n[k,N_dr+1]
        SPM[k,14] = SPM[k,13]/c_max_n
        SPM[k,15] = 2*R*T/F*asinh(SPM[k,11]*F/(k_n*SPM[k,14]^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-SPM[k,14])^(1-alpha_n))/2)
        SPM[k,16] = OCP_n_spline(SPM[k,14])
        SPM[k,17] = SPM[k,15]+SPM[k,16]+SPM[k,11]*F*SPM[k,21]*R_sei
        SPM[k,18] = SPM[k,10]-SPM[k,17]
        SPM[k,3] = SPM[k,2]/SPM[k,18]
    end
    SPM[k,4] = -SPM[k,3]/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p)
    SPM_c_p[k,:] = PDE_p*SPM_c_p[k-1,:]
    SPM_c_p[k,N_dr+1] = SPM_c_p[k,N_dr+1] + (-SPM[k,4]/D_p)*(2*dt*D_p/dr_p+2*dt*D_p/R_p)
    SPM[k,5] = (4/3*pi*dr_p^3*(SPM_c_p[k,1]+SPM_c_p[k,2])/2+sum(4/3*pi*dr_p^3*(i^3-(i-1)^3)*(SPM_c_p[k,i]+SPM_c_p[k,i+1])/2 for i in 2:N_dr))/(4/3*pi*R_p^3)
    SPM[k,6] = SPM_c_p[k,N_dr+1]
    SPM[k,7] = SPM[k,6]/c_max_p
    SPM[k,8] = 2*R*T/F*asinh(SPM[k,4]*F/(k_p*SPM[k,7]^alpha_p*c_el^(1-alpha_p)*c_max_p*(1-SPM[k,7])^(1-alpha_p))/2)
    SPM[k,9] = OCP_p_spline(SPM[k,7])
    SPM[k,10] = SPM[k,8]+SPM[k,9]
    SPM[k,20] = SPM[k,17] - SPM[k,11]*F*SPM[k-1,21]*R_sei
    SPM[k,23] = -F*k_sei*SPM[k-1,22]*exp(-α_sei*F*SPM[k,20] /R/T)
    SPM[k,21] = SPM[k-1,21]-dt*SPM[k,23]*M_sei/ρ_sei/(2*F)
    SPM[k,22] = SPM[k,21]/D_EC*(SPM[k,23]/F)+c_bulk_EC
    SPM[k,11] = SPM[k,3]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n)-SPM[k,23]/F
    SPM_c_n[k,:] = PDE_n*SPM_c_n[k-1,:]
    SPM_c_n[k,N_dr+1] = SPM_c_n[k,N_dr+1] + (-SPM[k,11]/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
    SPM[k,13] = SPM_c_n[k,N_dr+1]
    SPM[k,14] = SPM[k,13]/c_max_n
    SPM[k,15] = 2*R*T/F*asinh(SPM[k,11]*F/(k_n*SPM[k,14]^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-SPM[k,14])^(1-alpha_n))/2)
    SPM[k,16] = OCP_n_spline(SPM[k,14])
    SPM[k,17] = SPM[k,15]+SPM[k,16]+SPM[k,11]*F*SPM[k,21]*R_sei
    SPM[k,18] = SPM[k,10]-SPM[k,17]

    SPM[k,5] = (4/3*pi*dr_p^3*(SPM_c_p[k,1]+SPM_c_p[k,2])/2+sum(4/3*pi*dr_p^3*(i^3-(i-1)^3)*(SPM_c_p[k,i]+SPM_c_p[k,i+1])/2 for i in 2:N_dr))/(4/3*pi*R_p^3)
    SPM[k,12] = (4/3*pi*dr_n^3*(SPM_c_n[k,1]+SPM_c_n[k,2])/2+sum(4/3*pi*dr_n^3*(i^3-(i-1)^3)*(SPM_c_n[k,i]+SPM_c_n[k,i+1])/2 for i in 2:N_dr))/(4/3*pi*R_n^3)
    SPM[k,40] = SPM[k-1,40] - SPM[k,3]*dt/3600/Q  #SoC through current
    SPM[k,41] = SPM[k-1,41] - (SPM[k,23]*(3*eps_n*El_length_n*El_width_n*El_thickness_n/R_n)*dt/3600)/Q #LLI
    SPM[k,25] = SPM[k-1,25] - (SPM[k,23]*(3*eps_n*El_length_n*El_width_n*El_thickness_n/R_n)*dt/3600)/Q*100
end

df_SPM = DataFrame()
df_SPM."Time" = SPM[:,1]
df_SPM."P" = SPM[:,2]
df_SPM."I" = SPM[:,3]
df_SPM."V" = SPM[:,18]
df_SPM."P_calc" = SPM[:,30]
df_SPM."J_p" = SPM[:,4]
df_SPM."c_avg_p" = SPM[:,5]
df_SPM."c_surf_p" = SPM[:,6]
df_SPM."theta_p" = SPM[:,7]
df_SPM."eta_p" = SPM[:,8]
df_SPM."OCP_p" = SPM[:,9]
df_SPM."phi_p" = SPM[:,10]
df_SPM."J_n" = SPM[:,11]
df_SPM."c_avg_n" = SPM[:,12]
df_SPM."c_surf_n" = SPM[:,13]
df_SPM."theta_n" = SPM[:,14]
df_SPM."eta_n" = SPM[:,15]
df_SPM."OCP_n" = SPM[:,16]
df_SPM."phi_n" = SPM[:,17]
df_SPM."V" = SPM[:,18]

df_SPM."SoC_p_avg" = (1 .-(df_SPM."c_avg_p"/c_max_p .-theta_soc1_p)./(theta_soc0_p-theta_soc1_p))*100
df_SPM."SoC_n_avg" = ((df_SPM."c_avg_n"/c_max_n .-theta_soc0_n)./(theta_soc1_n-theta_soc0_n))*100
df_SPM."SoC_p_surf" = (1 .-(df_SPM."theta_p".-theta_soc1_p)./(theta_soc0_p-theta_soc1_p))*100
df_SPM."SoC_n_surf" = ((df_SPM."theta_n".-theta_soc0_n)./(theta_soc1_n-theta_soc0_n))*100

df_SPM."j_sei" = SPM[:,23]
df_SPM."eta_sei" = SPM[:,20]
df_SPM."SEI_thickness" = SPM[:,21]
df_SPM."EC_conc" = SPM[:,22]
df_SPM."LLI" = SPM[:,25]
