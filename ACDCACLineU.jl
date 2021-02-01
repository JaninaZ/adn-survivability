# regard the converter as a proportional amplifier
import Base: @__doc__
import PowerDynamics:  dimension, construct_edge, AbstractLine,PowerGrid
using PowerDynamics: PiModelLine
using NetworkDynamics#: ODEEdge
#using OrdinaryDiffEq: ODEFunction
using DifferentialEquations
# using GraphPlot

dir = @__DIR__
#include("$dir/CIGRE_static_2ADN.jl")
#include("$dir/control.jl")

###### the constant values are defined ######
begin
    const base_power = 1E6 # 1MW
    const base_voltage = 20E3 # 20kV
    const base_current = base_power / base_voltage # 50A
    const base_admittance = base_power / base_voltage^2 # 0.0025Ω^-1
    const ω = 2 * π * 50.0 # 314.1593rad/s
    # per unit HV
    const base_voltage_HV = 110E3 # 110kV
    const base_admittance_HV = base_power / base_voltage_HV^2 # 8.264462809917356e-5
    const base_current_HV = base_power / base_voltage_HV
    const r = 0.0178  # overhead line
    const l = 1.415E-3   # mH/km/pole
    const c = 0.0139E-6 # 20uF/km/pole page11
    const transmission_length = 100E3  # unit km
end

quad = 0.0 #voltage independent when quad = 0
verbose = false

###### the values used from here are all defined in the thesis mentioned below ######
# Z_c = 1.57 + im * (0.05 * ω)   # deined in article "analysis of VSC-based HVDC system" page 35
C_dc = ω * c * transmission_length / 4
C_conv = ω * 20E-6
Z_dcconv =  r * transmission_length + 1im * (ω * l * transmission_length - 1 / (2*C_dc + 2*C_conv))
Y_dcconv = 1 / Z_dcconv
k = 0.5  # from paper"An Equivalent Model for simulating VSC Based HVDC"
M = (1/(k^2))


include("$dir/CIGRE_static_2ADN.jl")


###### construct 2 ADNs and 1 DC bus line######
# this part does not differ much as before, but change the P_limit of the 2nd ADN

#function Grid_ACDC(P_ref1,P_ref2,Q_ref1,Q_ref2)
    busses_static1, lines1, T1, elist1, Zs1, Yshs1 = CIGRE_static_ADN1()#1-12
    busses_static2, lines2, T2, elist2, Zs2, Yshs2 = CIGRE_static_ADN2()#13-24
    
    
    DCline = [StaticLine(;from=1, to=13, Y=Y_dcconv*M/base_admittance_HV)]

    busses_static = []
    append!(busses_static, busses_static1)
    append!(busses_static, busses_static2) # bus from 1 to 24
    lines = []
    append!(lines, lines1)
    append!(lines, lines2)
    append!(lines, DCline) # all the lines

    pg_static = PowerGrid(busses_static, lines) 
    power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing)

    busses = copy(busses_static)

    DG_locs1 = collect(2:12) 
    DG_locs2 =  collect(14:24)
    DG_locs = []
    append!(DG_locs, DG_locs1) # located, 1 for slack, so 2 for MV1
    append!(DG_locs, DG_locs2)

    # two loops for two ADNs(set the DG unit)
    for i in DG_locs1 # this is a loop 

        S_bkgrnd = zero(im)

    # try

    #     S_bkgrnd = busses_static[i].S

    # catch

    #     S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)

    # end

        busses[i] = DGUnit(;

        K_pll=1632.993, # Hz/pu

        T_int=2.0,

        K_PT1=1.0,

        T_PT1=1e-8,

        K_FRT=2.0,

        I_max=1.0, # pu

        P_limit=1.0, # pu

        Q_limit=1.0, # pu

        Vref=1.0, # pu

        Vdead=0.1, # pu

        S_pq=V -> S_bkgrnd * (quad * V^2 + 1 - quad), # quad can be chosen 1 or 0, voltage dependent or not

        Y_n=0.0,

    )


    end
    for i in DG_locs2 # this is a loop 

        S_bkgrnd = zero(im)

    # try

    #     S_bkgrnd = busses_static[i].S

    # catch

    #     S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)

    # end

        busses[i] = DGUnit(;

        K_pll=1632.993, # Hz/pu

        T_int=2.0,

        K_PT1=1.0,

        T_PT1=1e-8,

        K_FRT=2.0,

        I_max=1.0, # pu

        P_limit=5.0, # pu

        Q_limit=1.0, # pu

        Vref=1.0, # pu

        Vdead=0.1, # pu

        S_pq=V -> S_bkgrnd * (quad * V^2 + 1 - quad), # quad can be chosen 1 or 0, voltage dependent or not

        Y_n=0.0,

    )


    end
   
    pg = PowerGrid(busses, lines)

    P_ref = 24
    Q_ref = 0 # DC part has no reactive power
    cpg = ADN(pg, DGUnit, t -> P_ref, t -> Q_ref) 
    #cpg = ADN(pg, DGUnit,  P_ref(t), Q_ref(t))
    fp = initial_guess(cpg, power_flow[:, :u])
    #power_flow_modified = power_flow
    #power_flow_modified[:, :u][end] = power_flow[:,:u][12]
    #fp_mod = initial_guess(cpg, power_flow_modified[:, :u])
    
    op = find_steady_state(cpg, fp) 
    # initial_guess(pg)  return State(pg, sol.u)
    verbose ? check_operationpoint(cpg, op) : nothing 

    # return pg, cpg
#end


#fp = find_valid_initial_condition(pg, ones(163))
#fp = initial_guess(cpg)
tspan = (0, 20.0)
sol = solve(pg, fp, tspan)
#sol = solve(cpg, op, tspan)
using Plots

plot(sol.dqsol,fmt=:png)

