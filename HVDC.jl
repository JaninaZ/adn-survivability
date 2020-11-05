import Base: @__doc__
using PowerDynamics#: @DynamicNode
#using PowerDynamics: AbstractLine, PiModel
import PowerDynamics:
    construct_vertex, dimension, showdefinition, AbstractNode
using LinearAlgebra: Diagonal

using Pkg
    #Pkg.activate(".")

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

    busses_static, lines, T, elist, Zs, Yshs = CIGRE_static(u_s) #cigre_static.jl  *Uos=110E3/base_voltage_HV???
        #const base_voltage_HV = 110E3
        #PowerDynamics.jl/src/nodes/SlackAlgebraic.jl
       pg_static = PowerGrid(busses_static, lines) #pkg PowerDynamics.jl
       #creates a PowerGrid from nodes and lines (either given as a list or as a dictionay).

       power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing) #pkg PowerDynamics.jl
    #initial_guess:  ~/control.jl   gives out a guess of voltage #???这里是已有关于u的guess
    #The voltage of all nodes is fixed to the voltage of the first SlackAlgebraic
    #in the system. The other symbols are set to zero.
    #pf_sol: power_flow.jl  #### solve power flow  return State(pg, pf)
    
   
       busses = copy(busses_static)  ## construct ADN: buses, MV1...11(DG_locs)
   
       DG_locs = 2:12 #located, 1 for slack, so 2 for MV1
   
       for i in DG_locs #this is a loop 
   
           S_bkgrnd = zero(im)
   
           try
   
               S_bkgrnd = busses_static[i].S
   
           catch
   
               S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)
   
        end
   
        busses[i] = DGUnit(;
   
            K_pll = 1632.993, #Hz/pu
   
            T_int = 2.0,
   
            K_PT1 = 1.0,
   
            T_PT1 = 1e-8,
   
            K_FRT = 2.0,
   
            I_max = 1.0, # pu
   
            P_limit = 1.0, #pu
   
            Q_limit = 1.0, #pu
   
            Vref = 1.0, # pu
   
            Vdead = 0.1, # pu
   
            S_pq = V -> S_bkgrnd * (quad * V^2 + 1 - quad), # quad can be chosen 1 or 0, voltage dependent or not
   
            Y_n = 0.0,
   
        )
   
    end
   
    pg = PowerGrid(busses, lines)
   
    cpg = ADN(pg, DGUnit, t -> P_ref, t -> Q_ref) #some problems in control.jl: DGUnit_Type not available
  # controlled power grid, P_ref, Q_ref are inputs
  # ADN Defined in control.jl

   
   
    ## find operation point
   
    icguess_pf = initial_guess(cpg, power_flow[:, :u]) # Defined in control.jl
   
    op = find_steady_state(cpg, icguess_pf) # Defined in control.jl
    # initial_guess(pg)  return State(pg, sol.u)
   
   
    verbose ? check_operationpoint(cpg, op) : nothing  #defined in control.jl
   #verbose: print all additional information
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
   P_ref = 0.9 * real(S_total) #certain value
   Q_ref = 0.9 * imag(S_total)
   u_s = 110E3 / base_voltage_HV  #
   quad = 0.0
   verbose = false

#  S_total = 	complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power
#  P_ref(t) = t > 0.25 ? 1.1 * real(S_total) : real(S_total)
#  Q_ref(t) = t > 0.25 ? 1.1 * imag(S_total) : imag(S_total)
   adn_pg, adn_cpg = ADN_construct(P_ref, Q_ref, u_s)
   #the code started here, then go to function ADN_construct
   
   #  to define which input parameters will need to be set for an ADN to make it flexible (Us?)



begin
    @__doc__ struct DCConverter <: AbstractNode
        b_i
        c_i
        r_i
        I_c_max
        V_c_max
        V_c_min
        Z_c
        B_f
        Y_dc
    end
    DCConverter(; b_i, c_i, r_i, I_c_max, V_c_max, V_c_min, Z_c, B_f, Y_dc) =
        DCConverter(b_i, c_i, r_i, I_c_max, V_c_max, V_c_min, Z_c, B_f, Y_dc )
        # rhs(pg::PowerGrid)
        # rhs(cpg::ControlledPowerGrid)
    function construct_vertex(par::DCConverter)
        b_i = par.b_i
        c_i = par.c_i
        r_i = par.r_i
        I_c_max = par.I_c_max
        V_c_max = par.V_c_max
        V_c_min = par.V_c_min
        Z_c = par.Z_c
        B_f = par.B_f
        Y_dc = par.Y_dc

        function converter(S_slack, p, u_s)
            #S_slack to be deined somewhere later, the value given from ADN
            u_s = complex(x[1], x[2]) #ADN side
            Vsamp = abs(u_s)
            u_c = complex(x[3], x[4])
            Vcamp = abs(u_c)

            I_s = conj(S_slack) / conj(u_s)

            I_c = I_s + u_s / Z_c

            u_c = u_s + Z_c * I_c

            #loss of converter
            P_loss = b_i + (c_i + r_i) * I_c^2

            #voltage at DC side
            U_dc = P_dc / I_dc
            #P_dc = U_dc * I_dc

            global_over_voltage, global_under_voltage = p

            # detect error state
            over_voltage = Vcamp > V_c_max
            under_voltage = Vcamp < V_c_min

            Y_dc = matrix_HVDC  #defined below
            voltage_vector = [source_voltage,destination_voltage]
            I_dc = Y_dc * (destination_voltage - source_voltage)

            #power flow at DC side
            if over_voltage | global_over_voltage 
                P_dc = I_c * V_c_max - P_loss
            elseif under_voltage |  global_under_voltage
                P_dc = I_c * V_c_min - P_loss
            else
                P_dc = I_c * u_c - P_loss
            end
            U_dc 
        end
        U_dc, P_dc
    end
end
#export construct_vertex

#converter = construct_vertex(par::DCConverter)

#node i,j the admittance matrix is symmetic y_ij=y_ji
function matrix_HVDC(y_ij)
    matrix_HVDC = zeros(Complex{Float64}, 2, 2)
    matrix_HVDC[1, 1] = y_ij 
    matrix_HVDC[1, 2] = - y_ij 
    matrix_HVDC[2, 1] = - y_ij
    matrix_HVDC[2, 2] = y_ij
    matrix_HVDC
end


#constant  AC power control using PI controller

function PI_controller_inner(P_ref, Q_ref, P_ref_last, Q_ref_last)
    ΔP = P_ref(t) - P_ref_last  #reference - measure value at ADN side
    ΔQ = Q_ref(t) - Q_ref_last  
    P_ref_last = real(S_slack) #the first loop
    Q_ref_last = imag(S_slack)

    K_p_inner = 1.0  #values of k to be discussed
    K_i_inner = 0.5
    switch_value = 1.0
    #C_i_inner = 1 #help to reduce the error in a shorter time

    if abs(ΔP) <= switch_value | abs(ΔP) <= switch_value
        P_ref = K_p_inner * P_ref_last + K_i_inner * integrate(P_ref_last)
        Q_ref = K_p_inner * Q_ref_last + K_i_inner * integrate(Q_ref_last)
        P_ref_last = P_ref
        Q_ref_last = Q_ref
    else
        #C_i_inner = 0
        P_ref = K_p_inner * P_ref_last 
        Q_ref = K_p_inner * Q_ref_last
        P_ref_last = P_ref
        Q_ref_last = Q_ref
    end
    P_ref, Q_ref
end
   
#constant DC voltage control using PI controller

function PI_controller_outer(U_dc)
    #U_dc_ref = 500E3  
    #U_ref_last = U_dc  #the first loop

    K_p_outer = 1.0  #values of k to be discussed
    K_i_outer = 0.5
    #C_i_inner = 1 #help to reduce the error in a shorter time

    if abs(U_dc) < 480E3
        U_dc_ref = K_p_outer * U_dc + K_i_outer * integrate(U_dc)
        U_dc = U_dc_ref
    else
        #C_i_inner = 0
        U_dc_ref = K_p_outer * U_ref_last
        U_dc = U_dc_ref
    end
    U_dc
end


P_dc = U_dc * I_dc

adn_power = P_ref + P_converter_loss

if (P_dc - adn_power) != 0
    P_ref = P_dc - P_converter_loss
    PI_controller_inner(P_ref, Q_ref, P_ref_last, Q_ref_last)
end

# # included in DCConverter
# @Line PiModelLine_dc(from, to, y_ij) begin
#     Y_dc = PiModel(y_ij)
# end begin
#     voltage_vector = [source_voltage,destination_voltage]
#     I_dc = Y_dc * (destination_voltage - source_voltage)
#     #
# end

# export PiModelLine_dc



# adn_dc_step, prob_step = toyexample(
#     S_slack, 
#     t -> real(S_slack),  #the measured value
#     t -> imag(S_slack); 
    
#     P_ref = 0.0, #P_ref的值是给出的，还是计算得出的?
#     Q_ref = 0.0,
#     t_duration = 0.,
#     resistance = 3.5,# to be defined 
#     quad = 0.,  # here means independent of u at adn side
#     verbose = false,
#     tspan = (0.,50.),
#     adn1 = adn_cpg(P_ref1, Q_ref1), 
#     adn2 = adn_cpg(P_ref2, Q_ref2),
#     converter1 = converter(),
#     converter2 = converter(),
#     )


# ΔP = P_ref(t) - real(S_slack)
# ΔQ = Q_ref(t) - imag(S_slack)
# P_ref1 = - P_ref2
# Q_ref1 = - Q_ref2




# ode = remake(prob_step, callback = nothing) #remake?
# dqsol = solve(ode, Rodas4())  #Rosenbrock函数是一个用来测试最优化算法性能的非凸函数non-convex function,此处是用该方法求微分方程
# sol = PowerGridSolution(dqsol, adn_dc_step) #PowerDynamics.jl/src/common/PowerGridSolutions.jl 按步长
# Δ = [dqsol.prob.f.f.cpg.controller(u, nothing, t) |> first for (t, u) in zip(dqsol.t, dqsol.u)]
# #zip把前项的元素与后项的元素按顺序对应起来   dqsol.prob.f.f.cpg???
# #|> Applies a function to the preceding argument. This allows for easy function chaining.
# using Plots

# plot(sol, 2:12, :P_g, legend=false)
# plot(dqsol.t, first.(Δ))




# cb1 = DiscreteCallback(
#     ((u, t, integrator) -> t in nsc.tspan_fault[1]),  #这里的nsc重构
# #nsc = NodeShortCircuitAlt(; node_number = nsc_node, Y = 1 / Znsc, tspan_fault = (t_fault, t_fault + t_duration))
#     errorState,
# )
# cb2 = DiscreteCallback(
#     ((u, t, integrator) -> t in nsc.tspan_fault[2]),
#     regularState,
# )
# if get_cpg
#     return cpg,
#     ODEProblem(
#         ode,
#         op.vec,
#         tspan,
#        callback = CallbackSet(cb1, cb2),
#        tstops = [t_fault, t_fault + t_duration],
#     )
# else
#     return ODEProblem(
#         ode,
#         op.vec,
#         tspan,
#        callback = CallbackSet(cb1, cb2),
#        tstops = [t_fault, t_fault + t_duration],
#     )
# end
