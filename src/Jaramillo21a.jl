"""
#######################################
# Jaramillo et al. (2021) model       #
# E: wave energy                      #
# dt: time step                       #
# a: model's parameter                #
# b: model's parameter                #
# Lccw: counterclockwise proportionality constants           #
# Lcw: Clockwise proportionality constants         #
# Yi: initial position                #
# vlt: longshore velocity             #
#######################################
"""

function Jaramillo21a(P, θ,dt, a, b, Lcw, Lccw, α0)
    
    αeq = ( θ .- b) ./ a
    
    α = zeros(length(P))
    
    α[1] = α0

    for i in eachindex(P[2:end])
        if α[i] < αeq[i+1]
            α[i+1] = ((α[i]-αeq[i+1]).*exp(-1. * Lcw * (P[i+1])*dt)) + αeq[i+1]
        else
            α[i+1] = ((α[i]-αeq[i+1]).*exp(-1. * Lccw * (P[i+1])*dt)) + αeq[i+1]
        end
    end
    
    return α
end


function cal_Jaramillo21a()

    """
    cal_Jaramillo21
    
    Function to calibrate and run the Jaramillo et al. (2021) Shoreline Rotation Model.
    
    This function reads input datasets, performs calibration, and writes the results to an output NetCDF file.
    
    """

    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data_rot/"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Reading datasets...")

    wavF = dats*"wav.nc"
    ensF = dats*"ens.nc"

    configF = dats*"config.nc"

    dt = ncread(configF, "dt")[1]
    calPar = ncread(configF, "calPar")[1]
    
    brk, angBati, depth= ncread(configF, "brk")[1], ncread(configF, "angBati")[1], ncread(configF, "depth")[1]

    MetObj = ncread(configF, "MetObj")[1]

    Hs, Tp, θ_w = ncread(wavF, "Hs"), ncread(wavF, "Tp"), ncread(wavF, "Dir")

    YY, MM, DD, HH = ncread(wavF, "Y"), ncread(wavF, "M"), ncread(wavF, "D"), ncread(wavF, "h")

    YYo, MMo, DDo, HHo = ncread(ensF, "Y"), ncread(ensF, "M"), ncread(ensF, "D"), ncread(ensF, "h")
    
    Y_obs = ncread(ensF, "Obs")

    t_obs = DateTime.(YYo, MMo, DDo, HHo)
    t_wav = DateTime.(YY, MM, DD, HH)

    ii =  t_obs .<= t_wav[end] .&& t_obs .>= t_wav[1]

    t_obs, Y_obs = t_obs[ii], Y_obs[ii]

    ii =  t_wav .<= t_obs[end] .&& t_wav .>= t_obs[1]

    t_wav, hb, tp, hs, θ_b = t_wav[ii], Hs[ii], Tp[ii], Hs[ii], θ_w[ii]

    idx_obs = zeros(length(t_obs))

    for i in eachindex(t_obs)
        idx_obs[i] = argmin(abs.(t_wav .- t_obs[i]))
    end

    idx_obs = convert(Array{Int64},idx_obs)

    ########## START HERE #############

    P = hb .^ 2 .* tp


    println("Starting Jaramillo et al. (2021) - Shoreline Evolution Model...")

    Ymdr = Dict()
    aRP = Dict()
    aRMSE = Dict()
    aMSS = Dict()
    popr = Dict()
    objr = Dict()

    if calPar == 4
        function Calibra_4r(Χ)
            Ymd = Jaramillo21a(P, θ_b,dt, Χ[1], Χ[2], exp(Χ[3]), exp(Χ[4]), Y_obs[1])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(-1e+2, 1e+2),
                   (-1e+3, 1e+3),
                   (log(1e-9), log(1e-3)),
                   (log(1e-9), log(1e-3))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_4r; 
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_4r; 
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_4r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo21a(P, θ_b, dt, popr[1], popr[2], exp(popr[3]), exp(popr[4]), Y_obs[1])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        println("\n\n****************Finished****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA21a.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        Y_atts = Dict("units" => "degrees",
            "long_name" => "Shoreline orientation",
            "standard_name" => "alpha")
        Yi_atts = Dict("units" => "degrees",
            "long_name" => "Initial shoreline orientation",
            "standard_name" => "alpha0")
        Lcw_atts = Dict("units" => "-",
            "long_name" => "Clockwise proportionality constants",
            "standard_name" => "Lcw")
        Lccw_atts = Dict("units" => "-",
            "long_name" => "counterclockwise proportionality constants",
            "standard_name" => "Lccw")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")


        nccreate(output, "alpha",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "alpha")
        nccreate(output, "alpha0",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([Y_obs[1]], output, "alpha0")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([popr[1]], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([popr[2]], output, "b")
        nccreate(output, "Lcw",
                    "len", 1,
                    atts = Lcw_atts)
        ncwrite([exp(popr[3])], output, "Lcw")
        nccreate(output, "Lccw",
                    "len", 1,
                    atts = Lccw_atts)
        ncwrite([exp(popr[4])], output, "Lccw")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")

    elseif calPar == 5
        function Calibra_5r(Χ)
            Ymd = Jaramillo21a(P, θ_b, dt, Χ[1], Χ[2], exp(Χ[3]), exp(Χ[4]), Χ[5])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        boundsr = [(-1e+2, 1e+2),
                   (-1e+3, 1e+3),
                   (log(1e-9), log(1e-3)),
                   (log(1e-9), log(1e-3)),
                   (0.25*minimum(Y_obs), 2*maximum(Y_obs))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_5r; 
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_5r; 
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_5r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo21a(P, θ_b, dt, popr[1], popr[2], exp(popr[3]), exp(popr[4]), popr[5])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        println("\n\n****************Writing output****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA21a.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        Y_atts = Dict("units" => "degrees",
            "long_name" => "Shoreline orientation",
            "standard_name" => "alpha")
        Yi_atts = Dict("units" => "degrees",
            "long_name" => "Initial shoreline orientation",
            "standard_name" => "alpha0")
        Lcw_atts = Dict("units" => "-",
            "long_name" => "Clockwise proportionality constants",
            "standard_name" => "Lcw")
        Lccw_atts = Dict("units" => "-",
            "long_name" => "counterclockwise proportionality constants",
            "standard_name" => "Lccw")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")

        nccreate(output, "alpha",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "alpha")
        nccreate(output, "alpha0",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([popr[5]], output, "alpha0")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([popr[1]], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([popr[2]], output, "b")
        nccreate(output, "Lcw",
                    "len", 1,
                    atts = Lcw_atts)
        ncwrite([exp(popr[3])], output, "Lcw")
        nccreate(output, "Lccw",
                    "len", 1,
                    atts = Lccw_atts)
        ncwrite([exp(popr[4])], output, "Lccw")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")

    end

    println("\n\n****************Finished****************\n\n")

end


function calVal_Jaramillo21a()

    """
    cal_Jaramillo21
    
    Function to calibrate and run the Jaramillo et al. (2021) Shoreline Rotation Model.
    
    This function reads input datasets, performs calibration, and writes the results to an output NetCDF file.
    
    """

    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data_rot/"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Reading datasets...")

    wavF = dats*"wav.nc"
    ensF = dats*"ens.nc"

    configF = dats*"config.nc"

    dt = ncread(configF, "dt")[1]
    calPar = ncread(configF, "calPar")[1]
    
    brk, angBati, depth= ncread(configF, "brk")[1], ncread(configF, "angBati")[1], ncread(configF, "depth")[1]

    MetObj = ncread(configF, "MetObj")[1]

    Hs, Tp, θ_w = ncread(wavF, "Hs"), ncread(wavF, "Tp"), ncread(wavF, "Dir")

    YY, MM, DD, HH = ncread(wavF, "Y"), ncread(wavF, "M"), ncread(wavF, "D"), ncread(wavF, "h")

    YYo, MMo, DDo, HHo = ncread(ensF, "Y"), ncread(ensF, "M"), ncread(ensF, "D"), ncread(ensF, "h")

    YYsi, MMsi, DDsi = ncread(configF, "Ysi")[1], ncread(configF, "Msi")[1], ncread(configF, "Dsi")[1]

    YYsf, MMsf, DDsf= ncread(configF, "Ysf")[1], ncread(configF, "Msf")[1], ncread(configF, "Dsf")[1]

    Y_obs = ncread(ensF, "Obs")

    t_obs = DateTime.(YYo, MMo, DDo, HHo)
    t_wav = DateTime.(YY, MM, DD, HH)
    t_split_i = DateTime.(YYsi, MMsi, DDsi, 0)
    t_split_f = DateTime.(YYsf, MMsf, DDsf, 0)

    ii =  t_obs .<= t_wav[end] .&& t_obs .>= t_wav[1]

    t_obs, Y_obs = t_obs[ii], Y_obs[ii]

    ii =  t_wav .<= t_obs[end] .&& t_wav .>= t_obs[1]

    t_wav, hb, tp, hs, θ_b = t_wav[ii], Hs[ii], Tp[ii], Hs[ii], θ_w[ii]

    idx_obs = zeros(length(t_obs))

    for i in eachindex(t_obs)
        idx_obs[i] = argmin(abs.(t_wav .- t_obs[i]))
    end

    idx_obs = convert(Array{Int64},idx_obs)

    idx_obs_split_cal = Int(argmin(abs.(t_obs .- t_split_i))) : Int(argmin(abs.(t_obs .- t_split_f)))

    idx_obs_split_val = [1 : Int(argmin(abs.(t_obs .- t_split_i))); Int(argmin(abs.(t_obs .- t_split_f))) : length(idx_obs)]

    t_obs_validation, Y_obs_validation, idx_obs_validation = t_obs[idx_obs_split_val], Y_obs[idx_obs_split_val], idx_obs[idx_obs_split_val]

    t_obs, Y_obs, idx_obs = t_obs[idx_obs_split_cal], Y_obs[idx_obs_split_cal], idx_obs[idx_obs_split_cal]

    ########## START HERE #############

    P = hb .^ 2 .* tp


    println("Starting Jaramillo et al. (2021) - Shoreline Evolution Model...")

    Ymdr = Dict()
    aRP = Dict()
    aRMSE = Dict()
    aMSS = Dict()
    popr = Dict()
    objr = Dict()

    if calPar == 4
        function Calibra_4r(Χ)
            Ymd = Jaramillo21a(P, θ_b,dt, Χ[1], Χ[2], exp(Χ[3]), exp(Χ[4]), Y_obs[1])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(-1e+2, 1e+2),
                   (-1e+3, 1e+3),
                   (log(1e-9), log(1e-4)),
                   (log(1e-9), log(1e-4))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_4r; 
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_4r; 
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_4r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo21a(P, θ_b, dt, popr[1], popr[2], exp(popr[3]), exp(popr[4]), Y_obs[1])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        Ysl = Ymdr[idx_obs_validation]
        aRP_validation = sum((Ysl.-mean(Ysl)).*(Y_obs_validation .- mean(Y_obs_validation)))/(std(Ysl)*std(Y_obs_validation)*length(Ysl))
        aRMSE_validation = sqrt(mean((Ysl .- Y_obs_validation).^2))
        aMSS_validation = 1 - sum((Ysl .- Y_obs_validation).^2)/length(Ysl)/(var(Ysl)+var(Y_obs_validation)+(mean(Ysl)-mean(Y_obs_validation))^2)

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA21a.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        yearv_atts = Dict("long_name" => "Year Validation")
        monthv_atts = Dict("long_name" => "Month Validation")
        dayv_atts = Dict("long_name" => "Day Validation")
        hourv_atts = Dict("long_name" => "Hour Validation")
        nccreate(output, "year_vi",
                    "len", 1,
                    atts = yearv_atts)
        ncwrite([YYsi], output, "year_vi")
        nccreate(output, "month_vi",
                    "len", 1,
                    atts = yearv_atts)
        ncwrite([MMsi], output, "month_vi")
        nccreate(output, "day_vi",
                    "len", 1,
                    atts = dayv_atts)
        ncwrite([DDsi], output, "day_vi")
        
        nccreate(output, "year_vf",
                    "len", 1,
                    atts = yearv_atts)
        ncwrite([YYsf], output, "year_vf")
        nccreate(output, "month_vf",
                    "len", 1,
                    atts = yearv_atts)
        ncwrite([MMsf], output, "month_vf")
        nccreate(output, "day_vf",
                    "len", 1,
                    atts = dayv_atts)
        ncwrite([DDsf], output, "day_vf")
        

        Y_atts = Dict("units" => "degrees",
            "long_name" => "Shoreline orientation",
            "standard_name" => "alpha")
        Yi_atts = Dict("units" => "degrees",
            "long_name" => "Initial shoreline orientation",
            "standard_name" => "alpha0")
        Lcw_atts = Dict("units" => "-",
            "long_name" => "Clockwise proportionality constants",
            "standard_name" => "Lcw")
        Lccw_atts = Dict("units" => "-",
            "long_name" => "counterclockwise proportionality constants",
            "standard_name" => "Lccw")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")
        RP_atts_v = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient Validation",
            "standard_name" => "RP_v")
        RMSE_atts_v = Dict("units" => "m",
            "long_name" => "Root mean square error Validation",
            "standard_name" => "RMSE_v")
        MSS_atts_v = Dict("units" => "-",
            "long_name" => "Mielke Skill Score Validation",
            "standard_name" => "MSS_v")


        nccreate(output, "alpha",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "alpha")
        nccreate(output, "alpha0",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([Y_obs[1]], output, "alpha0")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([popr[1]], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([popr[2]], output, "b")
        nccreate(output, "Lcw",
                    "len", 1,
                    atts = Lcw_atts)
        ncwrite([exp(popr[3])], output, "Lcw")
        nccreate(output, "Lccw",
                    "len", 1,
                    atts = Lccw_atts)
        ncwrite([exp(popr[4])], output, "Lccw")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")
        nccreate(output, "RP_v",
                    "len", 1,
                    atts = RP_atts_v)
        ncwrite([aRP_validation], output, "RP_v")
        nccreate(output, "RMSE_v",
                    "len", 1,
                    atts = RMSE_atts_v)
        ncwrite([aRMSE_validation], output, "RMSE_v")
        nccreate(output, "MSS_v",
                    "len", 1,
                    atts = MSS_atts_v)
        ncwrite([aMSS_validation], output, "MSS_v")

    elseif calPar == 5
        function Calibra_5r(Χ)
            Ymd = Jaramillo21a(P, θ_b, dt, Χ[1], Χ[2], exp(Χ[3]), exp(Χ[4]), Χ[5])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        boundsr = [(-1e+2, 1e+2),
                   (-1e+3, 1e+3),
                   (log(1e-9), log(1e-3)),
                   (log(1e-9), log(1e-3)),
                   (0.25*minimum(Y_obs), 2*maximum(Y_obs))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            resr = bboptimize(Calibra_5r; 
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_5r; 
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_5r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 5000,
                            MaxSteps = 25000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.0001,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 5000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Jaramillo21a(P, θ_b, dt, popr[1], popr[2], exp(popr[3]), exp(popr[4]), popr[5])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        Ysl = Ymdr[idx_obs_validation]
        aRP_validation = sum((Ysl.-mean(Ysl)).*(Y_obs_validation .- mean(Y_obs_validation)))/(std(Ysl)*std(Y_obs_validation)*length(Ysl))
        aRMSE_validation = sqrt(mean((Ysl .- Y_obs_validation).^2))
        aMSS_validation = 1 - sum((Ysl .- Y_obs_validation).^2)/length(Ysl)/(var(Ysl)+var(Y_obs_validation)+(mean(Ysl)-mean(Y_obs_validation))^2)

        println("\n\n****************Writing output****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_JA21a.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")

        yearv_atts = Dict("long_name" => "Year Validation")
        monthv_atts = Dict("long_name" => "Month Validation")
        dayv_atts = Dict("long_name" => "Day Validation")
        hourv_atts = Dict("long_name" => "Hour Validation")
        nccreate(output, "year_vi",
                    "len", 1,
                    atts = yearv_atts)
        ncwrite([YYsi], output, "year_vi")
        nccreate(output, "month_vi",
                    "len", 1,
                    atts = yearv_atts)
        ncwrite([MMsi], output, "month_vi")
        nccreate(output, "day_vi",
                    "len", 1,
                    atts = dayv_atts)
        ncwrite([DDsi], output, "day_vi")
        
        nccreate(output, "year_vf",
                    "len", 1,
                    atts = yearv_atts)
        ncwrite([YYsf], output, "year_vf")
        nccreate(output, "month_vf",
                    "len", 1,
                    atts = yearv_atts)
        ncwrite([MMsf], output, "month_vf")
        nccreate(output, "day_vf",
                    "len", 1,
                    atts = dayv_atts)
        ncwrite([DDsf], output, "day_vf")
        

        Y_atts = Dict("units" => "degrees",
            "long_name" => "Shoreline orientation",
            "standard_name" => "alpha")
        Yi_atts = Dict("units" => "degrees",
            "long_name" => "Initial shoreline orientation",
            "standard_name" => "alpha0")
        Lcw_atts = Dict("units" => "-",
            "long_name" => "Clockwise proportionality constants",
            "standard_name" => "Lcw")
        Lccw_atts = Dict("units" => "-",
            "long_name" => "counterclockwise proportionality constants",
            "standard_name" => "Lccw")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")
        RP_atts_v = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient Validation",
            "standard_name" => "RP_v")
        RMSE_atts_v = Dict("units" => "m",
            "long_name" => "Root mean square error Validation",
            "standard_name" => "RMSE_v")
        MSS_atts_v = Dict("units" => "-",
            "long_name" => "Mielke Skill Score Validation",
            "standard_name" => "MSS_v")

        nccreate(output, "alpha",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "alpha")
        nccreate(output, "alpha0",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([popr[5]], output, "alpha0")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([popr[1]], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([popr[2]], output, "b")
        nccreate(output, "Lcw",
                    "len", 1,
                    atts = Lcw_atts)
        ncwrite([exp(popr[3])], output, "Lcw")
        nccreate(output, "Lccw",
                    "len", 1,
                    atts = Lccw_atts)
        ncwrite([exp(popr[4])], output, "Lccw")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMSE",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE")
        nccreate(output, "MSS",
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS")
        nccreate(output, "RP_v",
                    "len", 1,
                    atts = RP_atts_v)
        ncwrite([aRP_validation], output, "RP_v")
        nccreate(output, "RMSE_v",
                    "len", 1,
                    atts = RMSE_atts_v)
        ncwrite([aRMSE_validation], output, "RMSE_v")
        nccreate(output, "MSS_v",
                    "len", 1,
                    atts = MSS_atts_v)
        ncwrite([aMSS_validation], output, "MSS_v")

    end

    println("\n\n****************Finished****************\n\n")

end