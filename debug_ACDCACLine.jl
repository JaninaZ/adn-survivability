import Base: @__doc__
using Pkg
Pkg.activate(".")

#using PowerDynamics
import PowerDynamics:
    SwingEq, PQAlgebraic, StaticLine,PowerGrid
using CSV
using OrderedCollections: OrderedDict
dir = @__DIR__
# for test 参考debug_dgunit.jl

include("$dir/ACDCACline.jl")
H = 0.8
D = 1.2
Ω = ω
P = -0.5
S = complex(-0.49998009701576795, -0.49758893704171214)
# buses=OrderedDict(
#     "bus1" => PQAlgebraic(P = -0.5,Q = -0.5),
#     "bus2" => SwingEq(;H, P, D, Ω));
# line = OrderedDict( "line1"=>StaticLine(;from="bus1",to ="bus2",Y=Y_dc));#pars = Y)

# pg = PowerGrid(buses, line)


bus1 = PQAlgebraic(P = -0.5,Q = -0.5)
bus2 = SwingEq(;H, P, D, Ω)
line = StaticLine(;from="bus1",to ="bus2",Y=Y_dc)#pars = Y)

pg = PowerGrid([bus1; bus2], line)
sol = solve(pg)