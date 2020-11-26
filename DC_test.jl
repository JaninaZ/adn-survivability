import Base: @__doc__
using PowerDynamics# : @DynamicNode
# using PowerDynamics: AbstractLine, PiModel
import PowerDynamics:
    construct_vertex, dimension, showdefinition, AbstractNode
using LinearAlgebra: Diagonal
# using SymPy

using Pkg
    # Pkg.activate(".")

using OrdinaryDiffEq
using Distributions
using Measurements

using DifferentialEquations
using Plots
using ModelingToolkit


function converter_construct(a_i, b_i, c_i, r_i, V_c_max, V_c_min, Z_c, S_slack, u_s)
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

    # Y_dc = matrix_HVDC  # defined below
    # voltage_vector = [source_voltage,destination_voltage]
    # I_dc = Y_dc * (destination_voltage - source_voltage)

    # power flow at DC side
    if over_voltage 
        P_dc = I_c * V_c_max + P_loss
    elseif under_voltage
        P_dc = I_c * V_c_min + P_loss
    else
        P_dc = I_c * u_c + P_loss
    end

    # voltage at DC side
    #U_dc = P_dc / I_dc
    # P_dc = U_dc * I_dc
    #U_dc, P_dc
    P_dc
end

#this part should be removed in the whole HVDC structure, the output of ADN part, here only for test
S_slack1 = 900E3
S_slack2 = 800E3
u_s = 110E3
const ω = 2 * π * 50.0 # 314.1593rad/s

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
Z_c = 1.57 + im * (0.05 * ω)   # deined in article "analysis of VSC-based HVDC system" page 35
Y_dc = 2 * r * transmission_length + im * (ω * l * transmission_length - 1/(ω * c * transmission_length))
# deined in article "analysis of VSC-based HVDC system" page 54
K_p = 0.  # to be calculated
K_i = 0.

# constant DC power control using PI controller

P_dc_1 = converter_construct(a_i, b_i, c_i, r_i, V_c_max, V_c_min, Z_c, S_slack1, u_s)
#println(P_dc_1)
P_dc_2 = converter_construct(a_i, b_i, c_i, r_i, V_c_max, V_c_min, Z_c, S_slack2, u_s)
#println(P_dc_2)

@parameters t Y_dc U_ref K_p K_i P_dc_1 P_dc_2 #coefficient
@variables P_dc_3(t) U_dc1(t) U_dc2(t) U_dc3(t) U_der(t) #variables
@derivatives D'~t

ΔU = U_ref - U_dc3

eqs = [0 ~ U_dc1*(Y_dc*(2*U_dc1-U_dc2-U_dc3)) - P_dc_1,
       0 ~ U_dc2*(Y_dc*(2*U_dc2-U_dc1-U_dc3)) - P_dc_2,
       0 ~ U_dc3*(Y_dc*(2*U_dc3-U_dc2-U_dc1)) - P_dc_3,
       D(P_dc_3) ~ K_p *  U_der + K_i * ΔU,
       D(ΔU) ~ U_der]

# eqs = [D(D(x)) ~ σ*(y-x),
#        D(y) ~ x*(ρ-z)-y,
#        D(z) ~ x*y - β*z]

HVDC_slack! = ODESystem(eqs)
# sys = ode_order_lowering(sys)

u0 = [P_dc_3 => 0.0,
      U_dc1 => 0.0,
      U_dc2 => 0.0,
      U_dc3 => 0.0,
      U_der => 0.0]

      
p  = [Y_dc,
      U_ref,
      K_p,
      K_i,
      P_dc_1,
      P_dc_2]

tspan = (0.0,10.0)
prob = ODEProblem(HVDC_slack!,u0,tspan,p)
sol = solve(prob,Tsit5())
plot(sol)



# function HVDC_slack!(out,du, u, t)
#     # PI controller
#     #ΔU = U_ref - U_dc[3]  #p[1] is U_ref
#     # dΔU = U_der
#     #u[1] = P_dc_controlled
    
#     #u[3] U_dc1, u[4] U_dc2, U_dc3 = u[5]
#     # P_dc_1 = u[3]*(Y_dc*(2*u[3]-u[4]-u[5]))
#     # P_dc_2 = u[4]*(Y_dc*(2*u[4]-u[3]-u[5]))
#     # u[1] = u[5]*(Y_dc*(2*u[5]-u[4]-u[5]))
#     # u[2] = U_ref - u[5]
#     out[1] = K_p * du[2] + K_i * u[2] - du[1]
#     out[2] = U_ref - u[5] - u[2]   # ΔU = U_ref - U_dc
#     out[3] = u[3]*(Y_dc*(2*u[3]-u[4]-u[5])) - P_dc_1
#     out[4] = u[4]*(Y_dc*(2*u[4]-u[3]-u[5])) - P_dc_2
#     out[5] = u[1] - u[5]*(Y_dc*(2*u[5]-u[4]-u[5]))
    
#     # P_dc_1 = U_dc[1]*(Y_dc*(2*U_dc[1]-U_dc[2]-U_dc[3]))
#     # P_dc_2 = U_dc[2]*(Y_dc*(2*U_dc[2]-U_dc[1]-U_dc[3]))
#     # u[1] = U_dc[3]*(Y_dc*(2*U_dc[3]-U_dc[2]-U_dc[1]))

#     #du[1] = K_p * du[2] + K_i * u[2]   #p[2] is K_p  p[3] is K_i
   
# end

# u0 = [1.0, 0.0, 0.0, 0.0, 0.0]
# du0 = [0.1,0.0,0.0,0.0,0.0]
# tspan = [0.0,1.0]
# # p = [ K_p, K_i]
# differential_vars = [true,true, false,false,false]
# prob = DAEProblem(HVDC_slack!,du0, u0, tspan,differential_vars=differential_vars)
# using Sundials
# sol = solve(prob,IDA())

# plot(sol)
# function HVDC_slack!(du, u, p, t)
#     P_dc_3,  ΔU = u 
#     # PI controller
#     ΔU = p[1] - U_dc  #p[1] is U_ref
        
#     P_dc_1 = U_dc[1]*(Y_dc*(2*U_dc[1]-U_dc[2]-U_dc[3]))
#     P_dc_2 = U_dc[2]*(Y_dc*(2*U_dc[2]-U_dc[1]-U_dc[3]))
#     P_dc_3 = U_dc[3]*(Y_dc*(2*U_dc[3]-U_dc[2]-U_dc[1]))
    
#     if abs(ΔU) < 1E3
#         dP_dc_3 = p[2] *  U_der + p[3] * ΔU   #p[2] is K_p  p[3] is K_i
#         dΔU = U_der
#     else
#         dP_dc_3 = p[2]  * dΔU
#     end
#     du = dP_dc_3, dΔU
# end
# K_p = 0.  # to be calculated
# K_i = 0.
# U_dc0 = [1.0;0.0;0.0]
# tspan = [0.0,1.0]
# p = [U_ref, K_p, K_i]
# prob = ODEProblem(HVDC_slack!, U_dc0, tspan, p)
# sol = solve(prob)

# plot(sol)

# plot(sol)
# function p_test!(du, u,p, t)
#     u[2]= u[1] + 1.0  #available like this
#     du[2] = u[1] * (28.0 - u[3]) - u[2]
#     du[3] = u[1] * u[2] - (8/3) * u[3]
# end

# u0 = [1.0;0.0;0.0]
# tspan = [0.0,10.0]
# prob = ODEProblem(p_test!, u0, tspan)
# sol = solve(prob)

# plot(sol)
# savefig("test.png")