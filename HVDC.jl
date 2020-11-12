import Base: @__doc__
using PowerDynamics# : @DynamicNode
# using PowerDynamics: AbstractLine, PiModel
import PowerDynamics:
    construct_vertex, dimension, showdefinition, AbstractNode
using LinearAlgebra: Diagonal
using SymPy

using Pkg
    # Pkg.activate(".")

using OrdinaryDiffEq
using Distributions
using Measurements

using Random
Random.seed!(42)

dir = @__DIR__

include("$dir/control.jl")
# custom types
include("$dir/DGUnit.jl")
include("$dir/OLTC.jl")
# load actual data
include("$dir/cigre_static.jl")
# failure model
include("$dir/short_circuit.jl")
# static solution
include("$dir/power_flow.jl")

function ADN_construct(P_ref, Q_ref, u_s)

    busses_static, lines, T, elist, Zs, Yshs = CIGRE_static(u_s) # cigre_static.jl  *Uos=110E3/base_voltage_HV???
        # const base_voltage_HV = 110E3
        # PowerDynamics.jl/src/nodes/SlackAlgebraic.jl
    pg_static = PowerGrid(busses_static, lines) # pkg PowerDynamics.jl
       # creates a PowerGrid from nodes and lines (either given as a list or as a dictionay).

    power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing) # pkg PowerDynamics.jl
    # initial_guess:  ~/control.jl   gives out a guess of voltage #???这里是已有关于u的guess
    # The voltage of all nodes is fixed to the voltage of the first SlackAlgebraic
    # in the system. The other symbols are set to zero.
    # pf_sol: power_flow.jl  #### solve power flow  return State(pg, pf)
    
   
    busses = copy(busses_static)  ## construct ADN: buses, MV1...11(DG_locs)
   
    DG_locs = 2:12 # located, 1 for slack, so 2 for MV1
   
    for i in DG_locs # this is a loop 
   
        S_bkgrnd = zero(im)
   
        try
   
            S_bkgrnd = busses_static[i].S
   
        catch
   
            S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)
   
        end
   
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
   
    pg = PowerGrid(busses, lines)
   
    cpg = ADN(pg, DGUnit, t -> P_ref, t -> Q_ref) # some problems in control.jl: DGUnit_Type not available
  # controlled power grid, P_ref, Q_ref are inputs
  # ADN Defined in control.jl

   
   
    ## find operation point
   
    icguess_pf = initial_guess(cpg, power_flow[:, :u]) # Defined in control.jl
   
    op = find_steady_state(cpg, icguess_pf) # Defined in control.jl
    # initial_guess(pg)  return State(pg, sol.u)
   
   
    verbose ? check_operationpoint(cpg, op) : nothing  # defined in control.jl
   # verbose: print all additional information
    return pg, cpg
end #  ADN construct
   
   # construct an ADN 

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
end

S_total = complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power
P_ref = 0.9 * real(S_total) # certain value
Q_ref = 0.9 * imag(S_total)
u_s = 110E3 / base_voltage_HV  # 
quad = 0.0
verbose = false

#  P_ref(t) = t > 0.25 ? 1.1 * real(S_total) : real(S_total)
#  Q_ref(t) = t > 0.25 ? 1.1 * imag(S_total) : imag(S_total)
adn_pg, adn_cpg = ADN_construct(P_ref, Q_ref, u_s)  # Top level part
# the code started here, then go to function ADN_construct



function converter_construct(a_i, b_i, c_i, r_i, V_c_max, V_c_min, Z_c, Y_dc, S_slack, u_s)
    # S_slack to be deined somewhere later, the value given from ADN
    # u_s = complex(x[1], x[2]) # ADN side
    # Vsamp = abs(u_s)
    # u_c = complex(x[3], x[4])

    I_s = conj(S_slack) / conj(u_s)

    I_c = I_s + u_s / Z_c

    u_c = u_s + Z_c * I_c
    Vcamp = abs(u_c)
    # loss of converter
    P_loss = a_i + b_i * I_c + (c_i + r_i) * I_c^2

    # detect error state
    over_voltage = Vcamp > V_c_max
    under_voltage = Vcamp < V_c_min

    Y_dc = matrix_HVDC  # defined below
    voltage_vector = [source_voltage,destination_voltage]
    I_dc = Y_dc * (destination_voltage - source_voltage)

            # power flow at DC side
    if over_voltage | global_over_voltage 
        P_dc = I_c * V_c_max + P_loss
    elseif under_voltage |  global_under_voltage
        P_dc = I_c * V_c_min + P_loss
    else
        P_dc = I_c * u_c + P_loss
    end

    # voltage at DC side
    U_dc = P_dc / I_dc
    # P_dc = U_dc * I_dc
    U_dc, P_dc
end

# deined in article "analysis of VSC-based HVDC system" page 11
const r = 0.0178
const l = 1.415   #mH/km/pole
const c = 0.0139 #uF/km/pole
const transmission_length = 1000E3  # unit km
# a_i, b_i, c_i, r_i, U_ref defined in paper "minimization of steady-state Losses in Meshed Networks using VSC HVDC"
U_ref = 300E3  # a rating of 600MW and a DC voltage of -300kV~+300kV  
a_i = 11.033E3  #inverter 11.033
b_i = 3.464E3   #inverter 3.464
c_i = 4.4E3  #inverter 6.667
r_i = 3 * 0.66   #inverter 3*1
V_c_max = 1.2 * U_ref
V_c_min = 0.85 * U_ref
Z_c = 1.57 + im * (0.05*ω)   # deined in article "analysis of VSC-based HVDC system" page 35
Y_dc = 2 * r * transmission_length + im * (ω * l * transmission_length- 1/(ω * c * transmission_length))
# deined in article "analysis of VSC-based HVDC system" page 54
converter1 = converter_construct(a_i, b_i, c_i, r_i, V_c_max, V_c_min, Z_c, Y_dc, S_slack, u_s)

# node i,j the admittance matrix is symmetic y_ij=y_ji
function matrix_HVDC(y_ij)
    matrix_HVDC = zeros(Complex{Float64}, 2, 2)
    matrix_HVDC[1, 1] = y_ij 
    matrix_HVDC[1, 2] = - y_ij 
    matrix_HVDC[2, 1] = - y_ij
    matrix_HVDC[2, 2] = y_ij
    matrix_HVDC
end


# constant DC power control using PI controller

function HVDC_slack(K_p, K_i, U_dc, U_ref)
    # PI controller
    ΔU = U_ref - U_dc
    U_der = dΔU
    C_i = 1

    if abs(ΔU) < 1E3
        C_i = 1
    else
        C_i = 0
    end

    P_dc_controlled = K_p * ΔU + C_i * K_i * U_der

    return P_dc_controlled
end

#U_ref = 300E3  # a rating of 600MW and a DC voltage of -300kV~+300kV
K_p = 1.2
K_i = 0.7
hvdc_slack = HVDC_slack(K_p, K_i, U_dc, U_ref) # Top level part
