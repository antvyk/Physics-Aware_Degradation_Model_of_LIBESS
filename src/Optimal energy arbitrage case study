#=
This code is written in Julia and can be used to obtain an optimal schedule of LiBESS
participating in the economic energy arbitrage. This LIBESS model is a blend of the linearized physics-based single particle
 model coupled with the description of SEI growth and energy-reservoir model.
 The incorporation of this model into the optimization framework results in a mixed-integer problem.
 This problem is solved using Gurobi solver.
=#

using DataFrames
using XLSX
using JuMP
using Gurobi





function build_1day(N_hours,N_minutes,Δt,Δt_h)


    directory_for_parameters = "D:\\"
    file_with_SPM_parameters = directory_for_parameters*"Parameters_for_SPM.xlsx"
    file_with_system_parameters = directory_for_parameters*"System_level_data.xlsx"
    SPM_Parameters = DataFrame(XLSX.readtable(file_with_SPM_parameters, "Parameters")...)
    LIBESS_Parameters = DataFrame(XLSX.readtable(file_with_system_parameters, "SystemLevel")...)

    ### Cell parameters
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
    V_nom = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Nominal Voltage","Nameplate data"][1]
    Q = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Cell Charge Capacity","Nameplate data"][1]
    j0_p_Const = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Exchange current density of the reaction","Positive electrode"][1]
    j0_n_Const = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Exchange current density of the reaction","Negative electrode"][1]
    Q_E = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Nominal energy capacity","Nameplate data"][1]

    ###   SEI parameters
    R_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI resistivity","SEI"][1]
    δ0_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Initial SEI thickness","SEI"][1]
    c_bulk_EC = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Concentration of EC in bulk electrolyte","SEI"][1]
    k_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Kinetic rate constant for SEI reaction","SEI"][1]
    α_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI charge transfer coefficient","SEI"][1]
    ρ_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI density","SEI"][1]
    M_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI molar weight","SEI"][1]
    D_EC = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI layer diffusivity","SEI"][1]

    #System Parameters
    Ch_max = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Maximum charging power","Value"][1]
    Dis_max = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Maximum discharging power","Value"][1]
    E_nom = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Energy capacity","Value"][1]/2*4
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

    #Solver's parameters
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 0.05)
    set_optimizer_attribute(model, "TimeLimit", 1500)
    set_optimizer_attribute(model,"NumericFocus", 0)
    set_optimizer_attribute(model,"Method", 2)

    #Initial Conditions at time 0:
    c_init_p = (theta_soc0_p-(theta_soc0_p-theta_soc1_p)*0.1)* c_max_p
    c_init_n = (theta_soc0_n+(theta_soc1_n-theta_soc0_n)*0.1)*c_max_n
    # Number of cells
    N = E_nom*1e6/Q_E
    #Start conditions for SEI:
    δ_sei_1w = δ0_sei
    EC_conc_1w  = c_bulk_EC
    LLI_0 = 0


    ### Set variables of the optimization model ###
    @variable(model, Pcell[1:N_hours]) # Supplied or consumed power
    @variable(model, SoE[1:N_hours]>=0) #State of Energy
    @variable(model, Ch[1:N_hours]>=0) #Charging power
    @variable(model, Dis[1:N_hours]>=0) #Discharging power
    @variable(model, u[1:N_hours],Bin)  #Binary variable to avoid simultaneous charging/discharging
    @variable(model, I[1:N_hours]) #Current
    @variable(model, V[1:N_hours,1:N_minutes]) #Voltage
    @variable(model, c_avg_p[1:N_hours,1:N_minutes]>=0) #Average Lithium ion concentration in the positive electrode
    @variable(model, c_surf_p[1:N_hours,1:N_minutes]>=0) #Surface Lithium ion concentration in the positive electrode
    @variable(model, theta_p[1:N_hours,1:N_minutes]>=0) # Dimensionless surface concentration of lithium ion, positive electrode
    @variable(model, OCP_p[1:N_hours,1:N_minutes]>=0) #OCP, positive electrode
    @variable(model, eta_p[1:N_hours,1:N_minutes]) # Activation overpotential, positive electrode
    @variable(model, phi_p[1:N_hours,1:N_minutes]>=0) #Solid-phase potential, positive electrode
    @variable(model, c_avg_n[1:N_hours,1:N_minutes]>=0) #Average Lithium ion concentration in the negative electrode
    @variable(model, c_surf_n[1:N_hours,1:N_minutes]>=0) #Surface Lithium ion concentration in the negative electrode
    @variable(model, theta_n[1:N_hours,1:N_minutes]>=0) #Dimensionless surface concentration of lithium ion, negative electrode
    @variable(model, OCP_n[1:N_hours,1:N_minutes]>=0) #OCP, negative electrode
    @variable(model, eta_n[1:N_hours,1:N_minutes]) # Activation overpotential, negative electrode
    @variable(model, phi_n[1:N_hours,1:N_minutes]) #>=0 #Solid-phase potential, negative electrode
    ### Set of Variables to linearize OCP for positive and negative electrodes###
    @variable(model, z_OCP_p[1:N_hours,1:N_minutes,1:2],Bin)
    @variable(model, λ_OCP_p[1:N_hours,1:N_minutes,1:3]>=0)
    @variable(model, z_OCP_n[1:N_hours,1:N_minutes,1:3],Bin)
    @variable(model, λ_OCP_n[1:N_hours,1:N_minutes,1:4]>=0)
    #############################################################
    @variable(model, η_sei[1:N_hours,1:N_minutes]) #overpotential of the SEIreaction
    @variable(model, δ_sei[1:N_hours,1:N_minutes]>=0) #Thickness of SEI layer
    @variable(model, c_EC[1:N_hours,1:N_minutes]>=0) # Concentration of EC on the surface of the negative electrode
    @variable(model, j_sei[1:N_hours,1:N_minutes]<=0) #SEI current density
    ### Set of Variables to linearize OCP for positive and negative electrodes###
    @variable(model, exp_j_sei[1:N_hours,1:N_minutes])
    @variable(model, z_sei[1:N_hours,1:N_minutes,1:4],Bin)
    @variable(model, λ_sei[1:N_hours,1:N_minutes,1:5]>=0)
    #############################################################
    @variable(model, LLI[1:N_hours,1:N_minutes]>=0) #Loss of lithium inventory
    @variable(model, SoC[1:N_hours,1:N_minutes]>=0) #State of Charge


    ### Set of constraints ###
    @constraint(model,Charging_Power1[h = 1:N_hours], Ch[h] <= Ch_max*(1-u[h]))
    @constraint(model,Charging_Power2[h = 1:N_hours], Ch[h] >= 0)
    @constraint(model,Discharging_Power1[h = 1:N_hours], Dis[h] <= Dis_max*u[h])
    @constraint(model,Discharging_Power2[h = 1:N_hours], Dis[h] >= 0)
    @constraint(model,State_of_Energy1[h = 1], SoE[h] == SoE_init + (η_RTE*Ch[h]-Dis[h])*Δt_h)
    @constraint(model,State_of_Energy2[h = 2:N_hours], SoE[h] == SoE[h-1] + (η_RTE*Ch[h]-Dis[h])*Δt_h)
    @constraint(model,State_of_Energy4[h = 1:N_hours],  SoE_min <= SoE[h])
    @constraint(model,State_of_Energy5[h = 1:N_hours],  SoE[h] <= SoE_max)
    @constraint(model,State_of_Energy6[h = N_hours],  SoE[h] >= SoE_final)
    @constraint(model, Power_cell[h = 1:N_hours], Pcell[h] == (Dis[h] - Ch[h])*1e6/N)
    @constraint(model, cell_current_aux1[h = 1:N_hours], I[h] == Pcell[h]/V_nom )
    @constraint(model, cell_voltage[h = 1:N_hours, m =1:N_minutes], V[h,m] == phi_p[h,m]-phi_n[h,m])
    @constraint(model, Average_concentration1_pe[h = 1, m = 1], c_avg_p[h,m] == c_init_p+Δt*(-3*(-I[h]/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p))/R_p))
    @constraint(model, Average_concentration2_pe[h = 2:N_hours, m = 1], c_avg_p[h,m] == c_avg_p[h-1,N_minutes]+Δt*(-3*(-I[h]/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p))/R_p))
    @constraint(model, Average_concentration3_pe[h = 1:N_hours, m = 2:N_minutes], c_avg_p[h,m] == c_avg_p[h,m-1]+Δt*(-3*(-I[h]/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p))/R_p))
    @constraint(model, Surface_concentration_pe[h = 1:N_hours, m =1:N_minutes], c_surf_p[h,m] == c_avg_p[h,m]+(-(-I[h]/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p))*R_p/(5*D_p)))
    @constraint(model, Dimensionless_surface_concentration_pe[h = 1:N_hours, m =1:N_minutes], theta_p[h,m] == c_surf_p[h,m]/c_max_p)
    @constraint(model, Chem_Reaction_overpot_pe[h = 1:N_hours, m =1:N_minutes], eta_p[h,m] == 2*R*T/F*((-I[h]/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p))*F/j0_p_Const/2))
    @constraint(model, OCP_pe[h = 1:N_hours, m =1:N_minutes], OCP_p[h,m] == λ_OCP_p[h,m,1] * 4.27 + λ_OCP_p[h,m,2] * 3.77 + λ_OCP_p[h,m,3] * 3.56)
    @constraint(model, OCP_pe1[h = 1:N_hours, m =1:N_minutes], theta_p[h,m] == λ_OCP_p[h,m,1] * theta_soc1_p + λ_OCP_p[h,m,2] * 0.6454 + λ_OCP_p[h,m,3] * theta_soc0_p)
    @constraint(model, OCP_pe2[h = 1:N_hours, m =1:N_minutes,j=1],  λ_OCP_p[h,m,j] <= z_OCP_p[h,m,j] )
    @constraint(model, OCP_pe3[h = 1:N_hours, m =1:N_minutes,j=2],  λ_OCP_p[h,m,j] <= z_OCP_p[h,m,j]+z_OCP_p[h,m,j-1] )
    @constraint(model, OCP_pe4[h = 1:N_hours, m =1:N_minutes,j=3],  λ_OCP_p[h,m,j] <= z_OCP_p[h,m,j-1] )
    @constraint(model, OCP_pe5[h = 1:N_hours, m =1:N_minutes],  λ_OCP_p[h,m,1]+λ_OCP_p[h,m,2]+λ_OCP_p[h,m,3] ==1)
    @constraint(model, OCP_pe6[h = 1:N_hours, m =1:N_minutes],  z_OCP_p[h,m,1]+z_OCP_p[h,m,2] ==1)
    @constraint(model, Solid_potential_pe[h = 1:N_hours, m =1:N_minutes], phi_p[h,m] == eta_p[h,m]+OCP_p[h,m])
    @constraint(model, Dimensionless_surface_concentration_pe_aux1[h = 1:N_hours, m =1:N_minutes],theta_p[h,m] <= theta_soc0_p)
    @constraint(model, Dimensionless_surface_concentration_pe_aux2[h = 1:N_hours, m =1:N_minutes], theta_soc1_p - theta_p[h,m] <= 0)
    @constraint(model, Average_concentration_pe_aux1[h = 1:N_hours, m =1:N_minutes], c_avg_p[h,m] <= c_max_p)
    @constraint(model, Average_concentration_pe_aux2[h = 1:N_hours, m =1:N_minutes], 0 - c_avg_p[h,m] <= 0)
    @constraint(model, Surface_concentration_pe_aux1[h = 1:N_hours, m =1:N_minutes], c_surf_p[h,m] <= c_max_p)
    @constraint(model, Surface_concentration_pe_aux2[h = 1:N_hours, m =1:N_minutes], 0 - c_surf_p[h,m] <= 0)
    @constraint(model, Average_concentration1_ne[h = 1, m = 1], c_avg_n[h,m] == c_init_n +Δt*(-3*(I[h]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n))/R_n))
    @constraint(model, Average_concentration2_ne[h = 2:N_hours, m = 1], c_avg_n[h,m] == c_avg_n[h-1,N_minutes]+Δt*(-3*(I[h]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n))/R_n))
    @constraint(model, Average_concentration3_ne[h = 1:N_hours, m = 2:N_minutes], c_avg_n[h,m] == c_avg_n[h,m-1]+Δt*(-3*(I[h]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n))/R_n))
    @constraint(model, Surface_concentration_ne[h = 1:N_hours, m =1:N_minutes], c_surf_n[h,m] == c_avg_n[h,m]+(-(I[h]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n))*R_n/(5*D_n)))
    @constraint(model, Dimensionless_surface_concentration_ne[h = 1:N_hours, m =1:N_minutes], theta_n[h,m] == c_surf_n[h,m]/c_max_n)
    @constraint(model, Chem_Reaction_overpot_ne[h = 1:N_hours, m =1:N_minutes], eta_n[h,m] == 2*R*T/F*((I[h]/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n))*F/j0_n_Const/2))
    @constraint(model, OCP_ne[h = 1:N_hours, m =1:N_minutes], OCP_n[h,m] ==  λ_OCP_n[h,m,1] * 1.07 + λ_OCP_n[h,m,2] * 0.5453 + λ_OCP_n[h,m,3] * 0.1923+λ_OCP_n[h,m,4] * 0.068)
    @constraint(model, OCP_ne1[h = 1:N_hours, m =1:N_minutes], theta_n[h,m] == λ_OCP_n[h,m,1] * theta_soc0_n + λ_OCP_n[h,m,2] * 0.05852 + λ_OCP_n[h,m,3] * 0.1707+ λ_OCP_n[h,m,4] * theta_soc1_n)
    @constraint(model, OCP_ne2[h = 1:N_hours, m =1:N_minutes,j=1],  λ_OCP_n[h,m,j] <= z_OCP_n[h,m,j] )
    @constraint(model, OCP_ne3[h = 1:N_hours, m =1:N_minutes,j=2],  λ_OCP_n[h,m,j] <= z_OCP_n[h,m,j]+z_OCP_n[h,m,j-1] )
    @constraint(model, OCP_ne4[h = 1:N_hours, m =1:N_minutes,j=3],  λ_OCP_n[h,m,j] <= z_OCP_n[h,m,j]+z_OCP_n[h,m,j-1] )
    @constraint(model, OCP_ne5[h = 1:N_hours, m =1:N_minutes,j=4],  λ_OCP_n[h,m,j] <= z_OCP_n[h,m,j-1] )
    @constraint(model, OCP_ne6[h = 1:N_hours, m =1:N_minutes],  λ_OCP_n[h,m,1]+λ_OCP_n[h,m,2]+λ_OCP_n[h,m,3]+λ_OCP_n[h,m,4] ==1)
    @constraint(model, OCP_ne7[h = 1:N_hours, m =1:N_minutes],  z_OCP_n[h,m,1]+z_OCP_n[h,m,2]+z_OCP_n[h,m,3] ==1)
    @constraint(model, Solid_potential_ne[h = 1:N_hours, m =1:N_minutes], phi_n[h,m] == eta_n[h,m]+OCP_n[h,m])
    @constraint(model, Average_concentration_ne_aux1[h = 1:N_hours, m =1:N_minutes], c_avg_n[h,m] <= c_max_n)
    @constraint(model, Average_concentration_ne_aux2[h = 1:N_hours, m =1:N_minutes], 0 - c_avg_n[h,m] <= 0)
    @constraint(model, Surface_concentration_ne_aux1[h = 1:N_hours, m =1:N_minutes], c_surf_n[h,m] <= c_max_n)
    @constraint(model, Surface_concentration_ne_aux2[h = 1:N_hours, m =1:N_minutes], 0-c_surf_n[h,m] <= 0)
    @constraint(model, Dimensionless_surface_concentration_ne_aux1[h = 1:N_hours, m =1:N_minutes], theta_n[h,m] <= theta_soc1_n)
    @constraint(model, Dimensionless_surface_concentration_ne_aux2[h = 1:N_hours, m =1:N_minutes], theta_soc0_n-theta_n[h,m] <= 0)
    @constraint(model, ActPotSEI[h = 1:N_hours, m =1:N_minutes], η_sei[h,m] == phi_n[h,m] - I[h]/(3*eps_n*El_length_n*El_width_n*El_thickness_n/R_n)*δ_sei_1w*R_sei)
    @constraint(model, j_sei_1[h = 1:N_hours, m =1:N_minutes], j_sei[h,m] == -F*k_sei*EC_conc_1w*exp_j_sei[h,m])
    @constraint(model, j_sei_2[h = 1:N_hours, m =1:N_minutes], exp_j_sei[h,m] == 0.00055308*λ_sei[h,m,1]+0.135*λ_sei[h,m,2]+0.368*λ_sei[h,m,3]+1*λ_sei[h,m,4]+1.65*λ_sei[h,m,5])
    @constraint(model, j_sei_3[h = 1:N_hours, m =1:N_minutes], η_sei[h,m] == (-1)*R*T/F/α_sei*(-7.5*λ_sei[h,m,1]-2*λ_sei[h,m,2]-1*λ_sei[h,m,3]+0*λ_sei[h,m,4]+0.5*λ_sei[h,m,5]))
    @constraint(model, j_sei_4[h = 1:N_hours, m =1:N_minutes,j=1],  λ_sei[h,m,j] <= z_sei[h,m,j] )
    @constraint(model, j_sei_5[h = 1:N_hours, m =1:N_minutes,j=2],  λ_sei[h,m,j] <= z_sei[h,m,j]+z_sei[h,m,j-1] )
    @constraint(model, j_sei_6[h = 1:N_hours, m =1:N_minutes,j=3],  λ_sei[h,m,j] <= z_sei[h,m,j]+z_sei[h,m,j-1] )
    @constraint(model, j_sei_7[h = 1:N_hours, m =1:N_minutes,j=4],  λ_sei[h,m,j] <= z_sei[h,m,j]+z_sei[h,m,j-1] )
    @constraint(model, j_sei_8[h = 1:N_hours, m =1:N_minutes,j=5],  λ_sei[h,m,j] <= z_sei[h,m,j-1] )
    @constraint(model, j_sei_9[h = 1:N_hours, m =1:N_minutes],  λ_sei[h,m,1]+λ_sei[h,m,2]+λ_sei[h,m,3]+λ_sei[h,m,4] +λ_sei[h,m,5]==1)
    @constraint(model, j_sei_10[h = 1:N_hours, m =1:N_minutes],  z_sei[h,m,1]+z_sei[h,m,2]+z_sei[h,m,3]+z_sei[h,m,4] ==1)
    @constraint(model, EC_concentration[h = 1:N_hours, m =1:N_minutes],c_EC[h,m]== δ_sei_1w/D_EC*(j_sei[h,m]/F)+c_bulk_EC)
    @constraint(model, LLI1[h = 1, m = 1], LLI[h,m] == LLI_0 - (j_sei[h,m]*(3*eps_n*El_length_n*El_width_n*El_thickness_n/R_n)*Δt/3600)/Q)
    @constraint(model, LLI2[h = 2:N_hours, m = 1], LLI[h,m] == LLI[h-1,N_minutes] - (j_sei[h,m]*(3*eps_n*El_length_n*El_width_n*El_thickness_n/R_n)*Δt/3600)/Q)
    @constraint(model, LLI3[h = 1:N_hours, m = 2:N_minutes], LLI[h,m] == LLI[h,m-1] - (j_sei[h,m]*(3*eps_n*El_length_n*El_width_n*El_thickness_n/R_n)*Δt/3600)/Q)
    @constraint(model, SoC_limit1[h = 1:N_hours, m =1:N_minutes], SoC[h,m] == ((c_avg_n[h,m]/c_max_n .-theta_soc0_n)./(theta_soc1_n-theta_soc0_n))*100)
    @constraint(model, SoC_limit2[h = 1:N_hours, m =1:N_minutes], SoC[h,m] >= SoE_min/E_nom*100)


    return model, Dis, Ch, LLI

end


directory = "D:\\"
file_with_market_data = directory*"Alberta_Pool_Price_2021.xlsx"
Data_market = DataFrame(XLSX.readtable(file_with_market_data, "2021")...)
directory_for_parameters = "D:\\"
file_with_SPM_parameters = directory_for_parameters*"Parameters for SPM - LG Chen2020Kane2022.xlsx"
file_with_system_parameters = directory_for_parameters*"System_level_data.xlsx"
SPM_Parameters = DataFrame(XLSX.readtable(file_with_SPM_parameters, "Parameters")...)
LIBESS_Parameters = DataFrame(XLSX.readtable(file_with_system_parameters, "SystemLevel")...)

#System level parameters
E_nom = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Energy capacity","Value"][1]/2*4
Q_E = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Nominal energy capacity","Nameplate data"][1]
C_instal = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="LIBESS capital cost","Value"][1]
V_nom = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Nominal Voltage","Nameplate data"][1]
Q = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Cell Charge Capacity","Nameplate data"][1]
C_oper = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Operation cost","Value"][1]


N_hours = 24*1
#Time resolution
N_minutes = 12
Δt = 3600/N_minutes
Δt_h = 1

# Number of cells
N = E_nom*1e6/Q_E




#Parameters for SEI
R_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI resistivity","SEI"][1]
δ0_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Initial SEI thickness","SEI"][1]
c_bulk_EC = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Concentration of EC in bulk electrolyte","SEI"][1]
k_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Kinetic rate constant for SEI reaction","SEI"][1]
α_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI charge transfer coefficient","SEI"][1]
ρ_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI density","SEI"][1]
M_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI molar weight","SEI"][1]
D_EC = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI layer diffusivity","SEI"][1]


#Initial conditions for SEI
δ_sei_1w = δ0_sei
EC_conc_1w  = c_bulk_EC
LLI_0 = 0


Results_in_array = zeros(N_hours*N_minutes,366)
Results_in_array_Obj = zeros(366,7)
Results_in_array[1:N_minutes,1] .= LinRange(1/N_minutes, 1, N_minutes)
#Time
for h in 2:N_hours
    Results_in_array[N_minutes*(h-1)+1:N_minutes*h,1] .= LinRange(1/N_minutes, 1, N_minutes).+(h-1) #Time
end


model, Dis, Ch,LLI= build_1day(N_hours,N_minutes,Δt,Δt_h)

for day in 1:365
    # Objective function
    #Find Daily operation schedule
    @objective(model, Max, sum(Δt_h*(Dis[h] - Ch[h])*Data_market[h+24*(day-1),"Pool Price"]  for h in 1:N_hours)
     -  LLI[N_hours,N_minutes]*V_nom*Q*N/1e6*C_instal
     -sum(Δt_h*Dis[h]*C_oper  for h in 1:N_hours)
     )
    println("________________________________________")
    println("DAY: ---> ",day)
    println("________________________________________")

    optimize!(model)
    Results_in_array[1:N_minutes,day+1] = JuMP.value.(Dis)[1]*ones(N_minutes,1)-JuMP.value.(Ch)[1]*ones(N_minutes,1)
    for h in 2:N_hours
        Results_in_array[N_minutes*(h-1)+1:N_minutes*h,day+1] = JuMP.value.(Dis)[h]*ones(N_minutes,1)-JuMP.value.(Ch)[h]*ones(N_minutes,1)
    end
    Results_in_array_Obj[day,1] = day
    Results_in_array_Obj[day,2] = JuMP.objective_value(model)
    Results_in_array_Obj[day,3] = JuMP.value(sum(Δt_h*(Dis[h] - Ch[h])*Data_market[h+24*(day-1),"Pool Price"]  for h in 1:N_hours))
    Results_in_array_Obj[day,4] = JuMP.value(sum(Δt_h*Dis[h]*C_oper  for h in 1:N_hours))
    Results_in_array_Obj[day,5] = JuMP.value(LLI[N_hours,N_minutes]*V_nom*Q*N/1e6*C_instal)
    Results_in_array_Obj[day,6] = JuMP.relative_gap(model)*100
    Results_in_array_Obj[day,7] = solve_time(model)

    println()
    println("Objective value: ---> ", JuMP.objective_value(model))
    println("Market Revenue: ---> ", JuMP.value(sum((Dis[h]-Ch[h])*Data_market[h+24*(day-1),"Pool Price"]*Δt_h  for h in 1:N_hours)))
    println("Cost of degradation: ---> ", JuMP.value(LLI[N_hours,N_minutes]*V_nom*Q*N/1e6*C_instal))
    println("Time to get a solution [s]: ---> ", solve_time(model))

end

#DataFrames with output from optimization
Results_in_df = DataFrame(Results_in_array,:auto)
Results_in_df_Obj = DataFrame(Results_in_array_Obj,[:Day,:Profit,:Revenue,:Operation_Cost,:Degradation_Cost,:MIP_Gap, :Solution_Time])


#Save optimal schedule and objectives to excel file
Results_output_file_name_disp = "DISPATCH_Hybrid_year2021"
if isfile(directory*Results_output_file_name_disp*".xlsx") == true
   rm(directory*Results_output_file_name_disp*".xlsx", force=true)
end
XLSX.writetable(directory*Results_output_file_name_disp*".xlsx", RESULTS=(collect(DataFrames.eachcol(Results_in_df)), DataFrames.names(Results_in_df) ))

Results_output_file_name_obj = "Objective_Hybrid_year2021"
if isfile(directory*Results_output_file_name_obj*".xlsx") == true
   rm(directory*Results_output_file_name_obj*".xlsx", force=true)
end
XLSX.writetable(directory*Results_output_file_name_obj*".xlsx", RESULTS=(collect(DataFrames.eachcol(Results_in_df_Obj)), DataFrames.names(Results_in_df_Obj) ))
