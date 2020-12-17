# regard the converter as a proportional amplifier
import Base: @__doc__
import PowerDynamics:  dimension, construct_edge, AbstractLine,PowerGrid
using PowerDynamics: PiModelLine
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


# a_i, b_i, c_i, r_i, U_ref defined in paper "minimization of steady-state Losses in Meshed Networks using VSC HVDC"
U_ref = 150E3  # a rating of 600MW and a DC voltage of -300kV~+300kV  or according to SPWM
P_m = 0.85
U_c_max = 1.05 * U_ref
U_c_min = 0.95 * U_ref
#Z_c = 1.57 + im * (0.05 * ω)   # deined in article "analysis of VSC-based HVDC system" page 35
# Z_dc = 2 * r * transmission_length + im * (ω * l * transmission_length - 1 / (ω * c * transmission_length))
# i_s = 500
# u_s = 100
const r = 0.0178  #overhead line
const l = 1.415E-3   # mH/km/pole
const c = 0.0139E-6 # 20uF/km/pole page11
const transmission_length = 100E3  # unit km

#Z_c = 1.57 + im * (0.05 * ω)   # deined in article "analysis of VSC-based HVDC system" page 35
C_dc = ω * c * transmission_length/4
C_conv = ω * 20E-6
Z_dc = 2 * r * transmission_length + 1im * (ω * l * transmission_length * 2 - 1/C_dc)
# the converter defined as power loss, which means P_dc = P_s = P_c
function converter_construct(u_s)
    U_dc_max = U_c_max/(0.612 * P_m)
    U_dc_min = U_c_min/0.612 #P_m = 1

    u_dc = 1.633 * u_s/P_m
    # detect error state
    over_voltage = u_dc > U_dc_max
    under_voltage = u_dc < U_dc_min

    # power flow at DC side
    if over_voltage 
        u_dc = U_dc_max
    elseif under_voltage
        u_dc = U_dc_min
    else
        u_dc
    end
    u_dc
end

DCline = [PiModelLine(;from=1, to = 12, y=1/Z_dc ,y_shunt_km =  C_conv * 1im,
 y_shunt_mk =  C_conv* 1im, )]

# ACDCACline = StaticLine(from, to, Y)

export ACDCACLine