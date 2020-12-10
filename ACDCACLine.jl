import Base: @__doc__
import PowerDynamics:  dimension, construct_edge, AbstractLine,PowerGrid
using PowerDynamics: StaticLine
using NetworkDynamics: ODEEdge
dir = @__DIR__
include("$dir/LineMacro.jl")
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

const r = 0.0178
const l = 1.415   # mH/km/pole
const c = 0.0139 # uF/km/pole
const transmission_length = 1000E3  # unit km
# a_i, b_i, c_i, r_i, U_ref defined in paper "minimization of steady-state Losses in Meshed Networks using VSC HVDC"
U_ref = 300E3  # a rating of 600MW and a DC voltage of -300kV~+300kV  
a_i = 11.033E3  # inverter 11.033
b_i = 3.464E3   # inverter 3.464
c_i = 4.4E3  # inverter 6.667
r_i = 3 * 0.66   # inverter 3*1
V_c_max = 1.2 * U_ref
V_c_min = 0.85 * U_ref
Z_c = 1.57 + im * (0.05 * ω)   # deined in article "analysis of VSC-based HVDC system" page 35
Y_dc = 2 * r * transmission_length + im * (ω * l * transmission_length - 1 / (ω * c * transmission_length))
i_s = 500
u_s = 100

function converter_construct(i_s, u_s)
   
    i_c = i_s + u_s / Z_c

    u_c = u_s + Z_c * i_c

    Vcamp = abs(u_c)

    # loss of converter
    P_loss = a_i + b_i * i_c + (c_i + r_i) * i_c^2

    # detect error state
    over_voltage = Vcamp > V_c_max
    under_voltage = Vcamp < V_c_min

    # power flow at DC side
    if over_voltage 
        P_dc = i_c * V_c_max + P_loss
    elseif under_voltage
        P_dc = i_c * V_c_min + P_loss
    else
        P_dc = i_c * u_c + P_loss
    end
    P_dc
end

U_dc_0 = 200
P_dc_0 = 300
K_dc = 3

# here we define the power flows from node 1 to 12, that is from ADN1 to ADN2
P_dc_1 = converter_construct(i_s, u_s)
P_dc_2 = - P_dc_1

# U_dc - P droop, the DC bus voltage U_dc depends the actual active power injection
# P_dc_0 is the setpoint, bus voltages are calculated with a Newton method
U_dc_1 = U_dc_0 - K_dc * (P_dc_1 - P_dc_0)

# the relationship between active power and DC voltage
U_dc_2 = U_dc_1 - P_dc_1 * (1/(U_dc_1 * Y_dc))
#P_dc_1 = U_dc[1] * Y_dc * (U_dc[1] - U_dc[2])
# P_dc_2 = U_dc[2] * Y_dc * (U_dc[2]-U_dc[1])

Y = Y_dc
destination_voltage = U_dc_2
source_voltage = U_dc_1

# @Line StaticLine(from, to, Y) begin

#     # If current is flowing away from the source, it is negative at the source.
#     # the current flowing in and out of the line is the same: I_mk=I_km
#     # hence, the current vector becomes only one complex current
#     complex_current = Y * (destination_voltage - source_voltage)
#     current_vector = [complex_current, complex_current]
# end

# DCline = StaticLine(from, to, Y)