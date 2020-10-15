import Base: @__doc__
using PowerDynamics: @DynamicNode

function ADN_construct()

    busses_static, lines, T, elist, Zs, Yshs = CIGRE_static()
   
       pg_static = PowerGrid(busses_static, lines)
   
       power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing)
   
   
   
       ## construct ADN
   
       busses = copy(busses_static)
   
       DG_locs = 2:12
   
       for i in DG_locs
   
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
   
            S_pq = V -> S_bkgrnd * (quad * V^2 + 1 - quad),
   
            Y_n = 0.0,
   
        )
   
    end
   
    pg = PowerGrid(busses, lines)
   
    cpg = ADN(pg, DGUnit, t -> P_ref, t -> Q_ref)
   
   
   
    ## find operation point
   
    icguess_pf = initial_guess(cpg, power_flow[:, :u])
   
    op = find_steady_state(cpg, icguess_pf)
   
   
   
    verbose ? check_operationpoint(cpg, op) : nothing
   
    return pg, cpg
end
   
   # construct an ADN 
   
   adn_example = ADN_construct()
   
   #  to define which input parameters will need to be set for an ADN to make it flexible (Us?)



begin
    @__doc__ struct DCGrid <: AbstractNode
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
    DCGird(; b_i, c_i, r_i, I_c_max, V_c_max, V_c_min, Z_c, B_f, Y_dc) =
        DCGrid(b_i, c_i, r_i, I_c_max, V_c_max, V_c_min, Z_c, B_f, Y_dc )
    function construct_vertex(par::DCGrid)
        b_i = par.b_i
        c_i = par.c_i
        r_i = par.r_i
        I_c_max = par.I_c_max
        V_c_max = par.V_c_max
        V_c_min = par.V_c_min
        Z_c = par.Z_c
        B_f = par.B_f
        Y_dc = par.Y_dc

        function rhs!(S_slack,p)
            #S_slack to be deined somewhere later, the value given from ADN
            u_s = complex(x[1], x[2])
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

            global_over_voltage, global_under_voltage = p

            # detect error state
            over_voltage = Vcamp > V_c_max
            under_voltage = Vcamp < V_c_min


            #power flow at DC side
            if over_voltage | global_over_voltage 
                P_dc = I_c * V_c_max - P_loss
            elseif under_voltage |  global_under_voltage
                P_dc = I_c * V_c_min - P_loss
            else
                P_dc = I_c * u_c - P_loss
            end
        end
    end
end


#node i,j the admittance matrix is symmetic y_ij=y_ji
function PiModel_dc(y_ij)
    Π = zeros(Complex{Float64}, 2, 2)
    Π[1, 1] = y_ij 
    Π[1, 2] = - y_ij 
    Π[2, 1] = - y_ij
    Π[2, 2] = y_ij
    Π
end


# in RLLine.jl already defined 
@Line PiModelLine_dc(from, to, y_ij) begin
    Y_dc = PiModel_dc(y_ij)
end begin
    voltage_vector = [source_voltage,destination_voltage]
    I_dc = Y_dc * (destination_voltage - source_voltage)
end

export PiModelLine_dc

