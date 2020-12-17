# make the converter and DC line the whole into an equivalent model
import Base: @__doc__
import PowerDynamics:  dimension, construct_edge, AbstractLine,PowerGrid
using PowerDynamics: PiModelLine
using NetworkDynamics: ODEEdge
dir = @__DIR__
include("$dir/LineMacro.jl")
include("$dir/HVDC.jl")
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

const r = 0.0178  #overhead line
const l = 1.415E-3   # mH/km/pole
const c = 0.0139E-6 # 20uF/km/pole page11
const transmission_length = 100E3  # unit km

#Z_c = 1.57 + im * (0.05 * ω)   # deined in article "analysis of VSC-based HVDC system" page 35
C_dc = ω * c * transmission_length/4
C_conv = ω * 20E-6
Z_dc = 2 * r * transmission_length + 1im * (ω * l * transmission_length * 2 - 1/C_dc)

DCline = [PiModelLine(;from=1, to = 12, y=1/Z_dc ,y_shunt_km =  C_conv * 1im,
 y_shunt_mk =  C_conv* 1im, )]
# use PIModelLine in this part

#this part comes from function ADN_construct(P_ref, Q_ref, u_s) in HVDC.jl
busses_static, lines, T, elist, Zs, Yshs = CIGRE_static(u_s)
busses = copy(busses_static)
line = append!(lines,DCline)

pg = PowerGrid(busses, line)
#pg = PowerGrid(busses, DCline)
fp = find_valid_initial_condition(pg, ones(5))
# fp = initial_guess(pg)
tspan =(0,0.2)
sol = solve(pg,fp,tspan)
using Plots