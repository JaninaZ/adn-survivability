import Base: @__doc__
using PowerDynamics
dir = @__DIR__
###### knoten ######
include("$dir/control.jl")
# custom types
include("$dir/DGUnit.jl")
include("$dir/OLTC.jl")
# load actual data
#include("$dir/cigre_static.jl")
# failure model
#include("$dir/short_circuit.jl")  #there is problem when running through this julia file
# static solution
include("$dir/power_flow.jl")

include("$dir/CIGRE_static_2ADN.jl")
include("$dir/ACDCACLineU.jl")

###### construct 2 ADNs and 1 DC bus line######
# this part does not differ much as before, but change the P_limit of the 2nd ADN

#function Grid_ACDC(P_ref1,P_ref2,Q_ref1,Q_ref2)
    busses_static1, lines1, T1, elist1, Zs1, Yshs1 = CIGRE_static_ADN1()#1-12
    ###### only for test(one ADN) ######
    pg_static1 = PowerGrid(busses_static1, lines1)
    using GraphPlot
    using LightGraphs: smallgraph
    # g = smallgraph(:karate)
    #gplot(pg)
    # nodelabel = 1:nv(pg.node)
    gplot(pg_static1.graph)
     
    power_flow1 = pf_sol(pg_static1, initial_guess(pg_static1), nothing)
    #power_flow1 = pf_sol(pg_static1, ones(12,1), nothing)

    # busses_static2, lines2, T2, elist2, Zs2, Yshs2 = CIGRE_static_ADN2()#13-24

    # ######only for test(one ADN, one SlackAlgebraic with one test line)######
    # busses_static2= [SlackAlgebraic(U = 110E3 / base_voltage_HV)]
    # C_dc = ω * c * transmission_length / 4
    # C_conv = ω * 20E-6
    # Z_dcconv =  r * transmission_length + 1im * (ω * l * transmission_length - 1 / (2 * C_dc + 2 * C_conv))
    # C1 = 0.1511749E-6 
    # ldata =  2.82
    # Ysh = 1im .* ω .* C1 .* ldata
    # linetest = [PiModelLine(;from=1, to = 13 ,y=1/Z_dcconv ,y_shunt_km = Ysh / 2.0 / base_admittance,
    # y_shunt_mk = Ysh / 2.0 / base_admittance, )]
    
    #DCline = ACDCACLine()

    busses_static = []
    append!(busses_static, busses_static1)
    append!(busses_static, busses_static2) # bus from 1 to 24
    lines = []
    append!(lines, lines1)
    #append!(lines, lines2)
    #append!(lines, DCline) # all the lines
    append!(lines, linetest) # for test

    pg_static = PowerGrid(busses_static, lines) 
    
    # ###### only for test ######
    # using GraphPlot
    # using LightGraphs: smallgraph
    # # g = smallgraph(:karate)
    # #gplot(pg)
    # # nodelabel = 1:nv(pg.node)
    # gplot(pg_static.graph)

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
    for i in DG_locs2 # this is a loop 

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

        P_limit=20.0, # pu

        Q_limit=1.0, # pu

        Vref=1.0, # pu

        Vdead=0.1, # pu

        S_pq=V -> S_bkgrnd * (quad * V^2 + 1 - quad), # quad can be chosen 1 or 0, voltage dependent or not

        Y_n=0.0,

    )


    end
   
    pg = PowerGrid(busses, lines)

    P_ref = 20 #ADN provide power to the higher grid
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