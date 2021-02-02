using PowerDynamics#: AbstractLine, AbstractNode


# the only difference between these functions is the elist
function CIGRE_static_ADN1()

    LoadP = [
        19839E3,
        0.0,
        501.7E3,
        431.65E3,
        727.5E3,
        548.05E3,
        76.5E3,
        586.85E3,
        573.75E3,
        543.3E3,
        329.8E3,
    ]


    cosϕ = [
        0.9737554,
        0.0,
        0.9231825,
        0.9700009,
        0.9699996,
        0.9700017,
        0.8500022,
        0.9699995,
        0.8499989,
        0.9586623,
        0.969997,
    ]

    LoadQ = LoadP .* sin.(acos.(cosϕ)) ./ cosϕ
    LoadQ[2] = 0.0

    begin
        busses_static1 = Array{PowerDynamics.AbstractNode,1}([])
        #push!(busses_static1, SlackAlgebraic(U = 110E3 / base_voltage_HV))
        push!(busses_static1, PQAlgebraic(P=0, Q=0)) #U = 110E3 / base_voltage_HV  U = u_s

        for (P, Q) in zip(LoadP, LoadQ)
            try
                push!(
                    busses_static1,
                    PQAlgebraic(S = - complex(P, Q) ./ base_power),
                )
            catch
                push!(
                    busses_static1,
                    PQAlgebraic(P = - P / base_power, Q = - Q / base_power),
                )
            end
        end
    end



    ###### Kanten ######



    T1 = OLTC(
        from = 1,
        to = 2,
        Uos = 110E3 / base_voltage_HV, # Bemessungsspannung Oberspannungsseite in kV
        Uus = 20E3 / base_voltage,# Bemessungsspannung Unterspannungsseite in kV
        k = 0, # Kennzahl der Schaltgruppe, Bsp: YD5-> 5
        ssp = 10.0, # Stufenschalterposition
        stufe_us = 0.625, # Stufung pro Stufenschalterstellung in %
        Sr = 25E6 / base_power, #Bemessungsscheinleistung in MVA
        uk = 12.0, # Kurzschlussspannung in %
        Pvk = 25E3 / base_power, # Kupferverluste in MW
        Pvl = 0.0, # Eisenverluste in kW
        iLeer = 0.0, # Leerlaufstrom in %
    )

    #ldata = CSV.read("$dir/lines.csv"; header = true)
    ldata1 = [ #L_km(nicht ändern)
        2.82,
        4.42,
        0.61,
        1.3,
        0.56,
        1.54,
        1.67,
        0.32,
        0.77,
        0.33,
        0.24,
    ]

    elist1 = [
        (2, 3),
        (3, 4),
        (4, 5),
        (4, 9),
        (5, 6),
        (6, 7),
        (8, 9),
        (9, 10),
        (10, 11),
        (11, 12),
        #(7, 8), # Schalter
        #(12, 5), # Schalter
    ]

    R1 = 0.501 # Ω/km
    X1 = 0.716 # Ω/km
    C1 = 0.1511749E-6 # F/km

    Zs1 = complex(R1, X1) .* ldata1 #[!, Symbol("L_km(nicht ändern)")]
    Yshs1 = 1im .* ω .* C1 .* ldata1 #[!, Symbol("L_km(nicht ändern)")]

    begin
        lines1 = Array{AbstractLine,1}([])
        push!(lines1, T1)
        for (e, Z, Ysh) in zip(elist1, Zs1, Yshs1)
            push!(
                lines1,
                PiModelLine(
                    from = first(e),
                    to = last(e),
                    y = inv(Z) / base_admittance,
                    y_shunt_km = Ysh / 2.0 / base_admittance,
                    y_shunt_mk = Ysh / 2.0 / base_admittance,
                ),
            )
        end
    end




    busses_static1, lines1, T1, elist1, Zs1, Yshs1
end


function CIGRE_static_ADN2()

    LoadP = [
        -19839E3,
        -0.0,
        -501.7E3,
        -431.65E3,
        -727.5E3,
        -548.05E3,
        -76.5E3,
        -586.85E3,
        -573.75E3,
        -543.3E3,
        -329.8E3,
    ]


    cosϕ = [
        0.9737554,
        0.0,
        0.9231825,
        0.9700009,
        0.9699996,
        0.9700017,
        0.8500022,
        0.9699995,
        0.8499989,
        0.9586623,
        0.969997,
    ]

    LoadQ = LoadP .* sin.(acos.(cosϕ)) ./ cosϕ
    LoadQ[2] = 0.0

    begin
        busses_static2 = Array{PowerDynamics.AbstractNode,1}([])
        push!(busses_static2, SlackAlgebraic(U = 110E3 / base_voltage_HV)) #U = 110E3 / base_voltage_HV  U = u_s

        for (P, Q) in zip(LoadP, LoadQ)
            try
                push!(
                    busses_static2,
                    PQAlgebraic(S = - complex(P, Q) ./ base_power),
                )
            catch
                push!(
                    busses_static2,
                    PQAlgebraic(P = - P / base_power, Q = - Q / base_power),
                )
            end
        end
    end



    ###### Kanten ######

    #include("$dir/ACDCACLineU.jl")

    T2 = OLTC(
        from = 13,
        to = 14,
        #Uos = 110E3 / base_voltage_HV + total_current(e_s, e_d)* k^2 / Y_dcconv,
        Uos = 110E3 / base_voltage_HV, # Bemessungsspannung Oberspannungsseite in kV
        Uus = 20E3 / base_voltage,# Bemessungsspannung Unterspannungsseite in kV
        k = 0, # Kennzahl der Schaltgruppe, Bsp: YD5-> 5
        ssp = 10.0, # Stufenschalterposition
        stufe_us = 0.625, # Stufung pro Stufenschalterstellung in %
        Sr = 25E6 / base_power, #Bemessungsscheinleistung in MVA
        uk = 12.0, # Kurzschlussspannung in %
        Pvk = 25E3 / base_power, # Kupferverluste in MW
        Pvl = 0.0, # Eisenverluste in kW
        iLeer = 0.0, # Leerlaufstrom in %
    )

    #ldata = CSV.read("$dir/lines.csv"; header = true)
    ldata2 = [ #L_km(nicht ändern)
        2.82,
        4.42,
        0.61,
        1.3,
        0.56,
        1.54,
        1.67,
        0.32,
        0.77,
        0.33,
        0.24,
    ]

    elist2 = [
        (14, 15),
        (15, 16),
        (16, 17),
        (16, 21),
        (17, 18),
        (18, 19),
        (20, 21),
        (21, 22),
        (22, 23),
        (23, 24),
        #(19, 20), # Schalter
        #(24, 17), # Schalter
    ]

    R1 = 0.501 # Ω/km
    X1 = 0.716 # Ω/km
    C1 = 0.1511749E-6 # F/km

    Zs2 = complex(R1, X1) .* ldata2 #[!, Symbol("L_km(nicht ändern)")]
    Yshs2 = 1im .* ω .* C1 .* ldata2 #[!, Symbol("L_km(nicht ändern)")]

    begin
        lines2 = Array{AbstractLine,1}([])
        push!(lines2, T2)
        for (e, Z, Ysh) in zip(elist2, Zs2, Yshs2)
            push!(
                lines2,
                PiModelLine(
                    from = first(e),
                    to = last(e),
                    y = inv(Z) / base_admittance,
                    y_shunt_km = Ysh / 2.0 / base_admittance,
                    y_shunt_mk = Ysh / 2.0 / base_admittance,
                ),
            )
        end
    end




    busses_static2, lines2, T2, elist2, Zs2, Yshs2
end
