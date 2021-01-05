# regard the converter as a proportional amplifier
import Base: @__doc__
import PowerDynamics:  dimension, construct_edge, AbstractLine,PowerGrid
using PowerDynamics: PiModelLine
using NetworkDynamics: ODEEdge
dir = @__DIR__
# include("$dir/LineMacro.jl")
include("$dir/CIGRE_static_2ADN.jl")
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

quad = 0.0
verbose = false

# Z_c = 1.57 + im * (0.05 * ω)   # deined in article "analysis of VSC-based HVDC system" page 35
C_dc = ω * c * transmission_length / 4
C_conv = ω * 20E-6
Z_dc = 2 * r * transmission_length + 1im * (ω * l * transmission_length * 2 - 1 / (C_dc + C_conv))

#after each step is checked, from line 34 to 107 will be included in a function
#function Grid_ACDC(P_ref1,P_ref2,Q_ref1,Q_ref2)
    #u_s can also be defined if needed in CIGRE_static_2ADN.jl
    busses_static1, lines1, T1, elist1, Zs1, Yshs1 = CIGRE_static_ADN1()#1-12
    busses_static2, lines2, T2, elist2, Zs2, Yshs2 = CIGRE_static_ADN2()#13-24
    
    busses_static = append!(busses_static1, busses_static2)#1-24
    lines = append!(lines1, lines2)

    pg_static = PowerGrid(busses_static, lines) 
    power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing)

    busses = copy(busses_static)

    DG_locs = 1:24 # located, 1 for slack, so 2 for MV1
    for i in DG_locs # this is a loop 

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
   
    k = 0.5  # from paper"An Equivalent Model for simulating VSC Based HVDC"
    DCline = [StaticLine(;from=1, to=13, Y=1 / Z_dc *(1/(k^2)))]

    line = append!(lines, DCline)
    pg = PowerGrid(busses, line)

    cpg = ADN(pg, DGUnit, t -> P_ref, t -> Q_ref) # some problems in control.jl: DGUnit_Type not available
    #cpg = ADN(pg, DGUnit,  P_ref(t), Q_ref(t))
    icguess_pf = initial_guess(cpg, power_flow[:, :u]) # Defined in control.jl

    op = find_steady_state(cpg, icguess_pf) # Defined in control.jl
    # initial_guess(pg)  return State(pg, sol.u)
    verbose ? check_operationpoint(cpg, op) : nothing 

    return pg, cpg
#end

pg = Grid_ACDC(P_ref1,P_ref2,Q_ref1,Q_ref2)

fp = find_valid_initial_condition(pg, ones(163))
# fp = initial_guess(pg)
tspan = (0, 20.0)
sol = solve(pg, fp, tspan)
#sol = solve(pg, op, tspan)
using Plots

plot(sol.dqsol,fmt=:png)


# U_ref = 150E3  # a rating of 600MW and a DC voltage of -300kV~+300kV  or according to SPWM
# P_m = 0.85
# U_c_max = 1.05 * U_ref
# U_c_min = 0.95 * U_ref
# the converter defined as power loss, which means P_dc = P_s = P_c
# function converter_construct(u_s)
#     U_dc_max = U_c_max/(0.612 * P_m)
#     U_dc_min = U_c_min/0.612 #P_m = 1

#     u_dc = 1.633 * u_s/P_m
#     # detect error state
#     over_voltage = u_dc > U_dc_max
#     under_voltage = u_dc < U_dc_min

#     # power flow at DC side
#     if over_voltage 
#         u_dc = U_dc_max
#     elseif under_voltage
#         u_dc = U_dc_min
#     else
#         u_dc
#     end
#     u_dc
# end

# DCline = [PiModelLine(;from=1, to = 12, y=1/Z_dc ,y_shunt_km =  C_conv * 1im,
#  y_shunt_mk =  C_conv* 1im, )]