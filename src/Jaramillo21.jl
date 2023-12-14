function run_Jaramillo2021()
    # Jaramillo21
    # This function calculates the angle of attack for a given
    # set of parameters and a given set of pressure values.
    # The function is based on the work of Jaramillo et al. (2021)
    # and is valid for the range of parameters used in the paper.
    # The function returns the angle of attack in radians.

    # INPUTS
    # P: Pressure values in Pa
    # alp0: Initial angle of attack in degrees
    # dt: Time step in seconds
    # a: Parameter a
    # b: Parameter b
    # cacr: Parameter cacr
    # cero: Parameter cero

    # OUTPUTS
    # α_t: Angle of attack in radians

    println("Loading libraries...")
    wrkDir = pwd()
    mods = wrkDir*"Modules\\"
    dats = wrkDir*"Data\\"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Loading datasets...")

    wavF = NCDataset(dats*"wav.nc")
    ensF = NCDataset(dats*"ens.nc")
    parF = NCDataset(dats*"par.nc")

    configF = NCDataset(dats*"config.nc")

    println("Unpacking datasets...")

    dt = configF["dt"][:][1]
    
    Hs = collect(skipmissing(wavF["hs"][:]))
    Tp = collect(skipmissing(wavF["tp"][:]))
    θ_w = collect(skipmissing(wavF["dir"][:]))
    
    yi = collect(skipmissing(configF["yi"][:]))

    a_ = collect(skipmissing(parF["a"][:]))
    b_ = collect(skipmissing(parF["b"][:]))
    Lacr_ = collect(skipmissing(parF["cacr"][:]))
    Lero_ = collect(skipmissing(parF["cero"][:]))
    α_0 = collect(skipmissing(parF["alpha0"][:]))
    
    close(wavF)
    close(ensF)
    close(configF)
    close(parF)

    println("Datasets closed...")

    Hs = convert(Array{Float64},Hs)
    Tp = convert(Array{Float64},Tp)
    θ_w = convert(Array{Float64},θ_w)

    a_ = convert(Array{Float64},a_)
    b_ = convert(Array{Float64},b_)
    Lacr_ = convert(Array{Float64},Lacr_)
    Lero_ = convert(Array{Float64},Lero_)

    yi = convert(Array{Float64},yi)

    ########## START HERE #############

    println("Starting Jaramillo et al. (2021) - Shoreline Rotation Model...")

    P = Hs .^2 .* Tp

    α_t = Jaramillo21(P, θ_w, α_0, dt, a, b, cacr, cero)

    println("\n\n****************Finished****************\n\n")

    return α_t
    
end
    

function Jaramillo21(P, θ_w, α_0, dt, a, b, cacr, cero)

    α_eq = (θ_w - b) / a
    α_t = fill(NaN, size(P))

    α_t[1] = α_0 * π / 180
    dt = fill(dt, length(α_t) - 1)

    for i in 1:length(α_t) - 1
        if α_t[i] < α_eq[i + 1]
            α_t[i + 1] = ((α_t[i] - α_eq[i + 1]) *
                    exp(cacr * (P[i + 1] * 3600))) + α_eq[i + 1]
        end

        if α_t[i] >= α_eq[i + 1]
            α_t[i + 1] = ((α_t[i] - α_eq[i + 1]) *
                    exp(cero * (P[i + 1]) * 3600)) + α_eq[i + 1]
        end
        
    end

    α_t = α_t .* π / 180

    return α_t
end


