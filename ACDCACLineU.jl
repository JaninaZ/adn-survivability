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
DCline = [StaticLine(;from=1, to=13, Y=Y_dcconv*M/base_admittance_HV)]


