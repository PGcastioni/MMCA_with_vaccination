## ----------------------------------------------------------------------------
## EPIDEMIC PARAMS RELATED FUNCTIONS
## ----------------------------------------------------------------------------

"""
    Epidemic_Params

Struct that contains the parameters related with the epidemic parameters and
compartmental evolution.

# Fields

All the parameters contained in this structure are probabilities ranged between
0 and 1.

# Epidemic parameters

- `βᴵ::Array{Float64, 1}`: Infectivity of infected.
- `βᴬ::Array{Float64, 1}`: Infectivity of asymptomatic.
- `ηᵍ::Array{Float64, 1}`: Exposed rate for each strata.
- `αᵍ::Array{Float64, 1}`: Asymptomatic infectious rate for each strata.
- `μᵍ::Array{Float64, 1}`: Infectious rate for each strata.
- `θᵍ::Array{Float64, 1}`: Direct death probability for each strata.
- `γᵍ::Array{Float64, 1}`: ICU probability for each strata.
- `ζᵍ::Array{Float64, 1}`: Pre-deceased rate for each strata.
- `λᵍ::Array{Float64, 1}`: Pre-hospitalized in ICU rate for each strata.
- `ωᵍ::Array{Float64, 1}`: Fatality probability in ICU for each strata.
- `ψᵍ::Array{Float64, 1}`: Death rate in iCU for each strata.
- `χᵍ::Array{Float64, 1}`: ICU discharge rate for each strata.
- `Λ::Float64`: Probability of losing the vaccine-acquired immunity
- `Γ::Float64`: Probability of losing the disease-acquired immunity
- `T::Int64`: Number of epidemic timesteps.
- `V::Int64`: Number of stages in the vaccination process
- `NumComps::Int64`: Number of epidemiological compartments
- `rᵥ::Array{Float64, 1}`: Vaccine efficacy in preventing infections
- `kᵥ::Array{Float64, 1}`: Vaccine efficacy in preventing tranmission

# Compartmental evolution

- `ρˢᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of suceptible individuals for each
  strata, patch and vaccination status.
- `ρᴱᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of exposed individuals for each
  strata, patch and vaccination status.
- `ρᴬᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of asymptomatic individuals for
  each strata, patch and vaccination status.
- `ρᴵᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of symptomatic individuals for each
  strata, patch and vaccination status.
- `ρᴾᴴᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of pre-hospitalized to ICU
  individuals for each strata, patch and vaccination status.
- `ρᴾᴰᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of pre-deceased individuals for each strata, patch 
  and vaccination status.
- `ρᴴᴿᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of hospitalized in ICU patients who
  will recover for each strata, patch and vaccination status.
- `ρᴴᴰᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of hospitalized in ICU patients who
  will not recover for each strata, patch and vaccination status.
- `ρᴰᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of deceased individuals for each
  strata, patch and vaccination status.
- `ρᴿᵍᵥ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T \\times V`` containing
  infomation about the evolution of fraction of recovered individuals for each
  strata, patch and vaccination status.

# Auxiliary

- `CHᵢᵍᵥ::Array{Float64, 3}`: Fraction of securely confined individuals for each
  strata and patch.
- `Qᵢᵍ::Array{Float64, 3}`: Suceptible contacts available for each strata on a
  given patch.

"""
struct Epidemic_Params
    #Epidemic parameters
    βᴵ::Array{Float64, 1}
    βᴬ::Array{Float64, 1}
    ηᵍ::Array{Float64, 1}
    αᵍ::Array{Float64, 1}
    μᵍ::Array{Float64, 1}
    θᵍ::Array{Float64, 2}
    γᵍ::Array{Float64, 2}
    ζᵍ::Array{Float64, 1}
    λᵍ::Array{Float64, 1}
    ωᵍ::Array{Float64, 2}
    ψᵍ::Array{Float64, 1}
    χᵍ::Array{Float64, 1}
    Λ::Float64
    Γ::Float64
    
    T::Int64
    V::Int64
    
    # Epidemic parameter for vaccinations
    rᵥ::Array{Float64, 1}
    kᵥ::Array{Float64, 1}

    # The total number of compartments
    NumComps::Int64
    CompLabels::Tuple

    # Compartments evolution
    ρˢᵍᵥ::Array{Float64, 4}
    ρᴱᵍᵥ::Array{Float64, 4}
    ρᴬᵍᵥ::Array{Float64, 4}
    ρᴵᵍᵥ::Array{Float64, 4}
    ρᴾᴴᵍᵥ::Array{Float64, 4}
    ρᴾᴰᵍᵥ::Array{Float64, 4}
    ρᴴᴿᵍᵥ::Array{Float64, 4}
    ρᴴᴰᵍᵥ::Array{Float64, 4}
    ρᴰᵍᵥ::Array{Float64, 4}
    ρᴿᵍᵥ::Array{Float64, 4}
    
    # Fraction of securely confined individuals for each strata and patch.
    CHᵢᵍᵥ::Array{Float64, 3}
    
    # R_t related arrays
    Qᵢᵍ::Array{Float64, 3}
end


"""
    Epidemic_Params(βᴵ::Float64,
                    βᴬ::Float64,
                    ηᵍ::Array{Float64, 1},
                    αᵍ::Array{Float64, 1},
                    μᵍ::Array{Float64, 1},
                    θᵍ::Array{Float64, 1},
                    γᵍ::Array{Float64, 1},
                    ζᵍ::Array{Float64, 1},
                    λᵍ::Array{Float64, 1},
                    ωᵍ::Array{Float64, 1},
                    ψᵍ::Array{Float64, 1},
                    χᵍ::Array{Float64, 1},
                    Λ::Float64,
                    Γ::Float64,                         
                    rᵥ::Array{Float64, 1},
                    kᵥ::Array{Float64, 1},
                    G::Int64,
                    M::Int64,
                    T::Int64,
                    V::Int64)

# Arguments

- `βᴵ::Float64`: Infectivity of infected for each vaccination status.
- `βᴬ::Float64`: Infectivity of asymptomatic for each vaccination status.
- `ηᵍ::Array{Float64, 1}`: Vector of size ``G`` with exposed rates for each
  strata.
- `αᵍ::Array{Float64, 1}`: Vector of size ``G`` with asymptomatic infectious
  rates for each strata.
- `μᵍ::Array{Float64, 1}`: Vector of size ``G`` with infectious rates for each
  strata.
- `θᵍ::Array{Float64, 1}`: Vector of size ``G`` with direct death probabilities
  for each strata.
- `γᵍ::Array{Float64, 1}`: Vector of size ``G`` with ICU probabilities for each
  strata.
- `ζᵍ::Array{Float64, 1}`: Vector of size ``G`` with pre-deceased rates for
  each strata.
- `λᵍ::Array{Float64, 1}`: Vector of size ``G`` with pre-hospitalized in ICU
  rates for each strata.
- `ωᵍ::Array{Float64, 1}`: Vector of size ``G`` with fatality probabilities in
  ICU for each strata.
- `ψᵍ::Array{Float64, 1}`: Vector of size ``G`` with death rates for each
  strata.
- `χᵍ::Array{Float64, 1}`: Vector of size ``G`` with ICU discharge rates for
  each strata.
- `Λ::Float64`: Probability of losing the vaccine-acquired immunity
- `Γ::Float64`: Probability of losing the disease-acquired immunity
- `rᵥ::Array{Float64, 1}`: Vaccine efficacy in preventing infections for each vaccination status.
- `kᵥ::Array{Float64, 1}`: Vaccine efficacy in preventing tranmission for each vaccination status.
- `G::Int64`: Number of strata.
- `M::Int64`: Number of patches.
- `T::Int64`: Number of epidemic timesteps.
- `V::Int64`: Number of stages in the vaccination process.

# Return

Struct that contains the parameters related with the epidemic parameters and
compartmental evolution.
"""
function Epidemic_Params(βᴵ::Float64,
                         βᴬ::Float64,
                         ηᵍ::Array{Float64, 1},
                         αᵍ::Array{Float64, 1},
                         μᵍ::Array{Float64, 1},
                         θᵍ::Array{Float64, 2},
                         γᵍ::Array{Float64, 2},
                         ζᵍ::Array{Float64, 1},
                         λᵍ::Array{Float64, 1},
                         ωᵍ::Array{Float64, 2},
                         ψᵍ::Array{Float64, 1},
                         χᵍ::Array{Float64, 1},
                         Λ::Float64,
                         Γ::Float64,                         
                         rᵥ::Array{Float64, 1},
                         kᵥ::Array{Float64, 1},
                         G::Int64,
                         M::Int64,
                         T::Int64,
                         V::Int64)

    NumComps = 10
    CompLabels = ("S", "E", "A", "I", "PH", "PD", "HR", "HD", "R", "D")
    
    # Allocate memory for simulations
    ρˢᵍᵥ  = ones(Float64, G, M, T, V)
    ρᴱᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴬᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴵᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴾᴴᵍᵥ = zeros(Float64, G, M, T, V)
    ρᴾᴰᵍᵥ = zeros(Float64, G, M, T, V)
    ρᴴᴿᵍᵥ = zeros(Float64, G, M, T, V)
    ρᴴᴰᵍᵥ = zeros(Float64, G, M, T, V)
    ρᴰᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴿᵍᵥ  = zeros(Float64, G, M, T, V)
    CHᵢᵍᵥ = zeros(Float64, G, M, V)
    Qᵢᵍ   = zeros(Float64, G, M, T)
    
    return Epidemic_Params([βᴵ], [βᴬ], copy(ηᵍ), copy(αᵍ), copy(μᵍ),
                           copy(θᵍ), copy(γᵍ), copy(ζᵍ), copy(λᵍ), copy(ωᵍ),
                           copy(ψᵍ), copy(χᵍ), copy(Λ), copy(Γ), T, V, copy(rᵥ), copy(kᵥ), 
                           NumComps, CompLabels, ρˢᵍᵥ, ρᴱᵍᵥ, ρᴬᵍᵥ, ρᴵᵍᵥ, ρᴾᴴᵍᵥ,
                           ρᴾᴰᵍᵥ, ρᴴᴿᵍᵥ, ρᴴᴰᵍᵥ, ρᴰᵍᵥ, ρᴿᵍᵥ, CHᵢᵍᵥ, Qᵢᵍ)
end


"""
    reset_epidemic_params!(epi_params::Epidemic_Params)

Reset the ρ's to reuse the structure and avoid additional allocations.

# Arguments
- `epi_params::Epidemic_Params`: structure to reset.
"""
function reset_epidemic_params!(epi_params::Epidemic_Params)
    epi_params.ρˢᵍᵥ .= 1.
    epi_params.ρᴱᵍᵥ .= 0.
    epi_params.ρᴬᵍᵥ .= 0.
    epi_params.ρᴵᵍᵥ .= 0.
    epi_params.ρᴾᴴᵍᵥ .= 0.
    epi_params.ρᴾᴰᵍᵥ .= 0.
    epi_params.ρᴴᴿᵍᵥ .= 0.
    epi_params.ρᴴᴰᵍᵥ .= 0.
    epi_params.ρᴰᵍᵥ .= 0.
    epi_params.ρᴿᵍᵥ .= 0.
    epi_params.CHᵢᵍᵥ .= 0.

    # Rt structures
    epi_params.Qᵢᵍ .= 0.

    nothing
end

"""
VACCINATION PARAMETERS

Struct that contains the parameters related with vaccination 
strageties

#Parameters
  - ϵᵍ:
  - tᵛs:
  - ϵᵍs:

"""



"""
CONTAINEMENT PARAMETERS

Struct that contains the parameters related with containement
measures

#Parameters
  - 
"""


### ----------------------------------------------------------------------------
### PATCH AND POPULATION RELATED FUNCTIONS
### ----------------------------------------------------------------------------

"""
    Population_Params

Struct that contains the parameters related with geographical, population and
mobility data.

# Fields

- `G::Int64`: Number of population strata.
- `M::Int64`: Number of patches.
- `nᵢ:Array{Float64, 1}`: Population of each patch.
- `nᵢᵍ:Array{Float64, 2}`: Population of each strata on each patch.
- `nᵢ_eff:Array{Float64, 1}`: Effective population on each patch after taking
  into account mobility.
- `nᵢᵍ_eff:Array{Float64, 1}`: Effective population of each strata on each patch
  after taking into account mobility.
- `N::Int64`: Total population.
- `Nᵍ::Array{Int64, 1}`: Total population of each strata.
- `pᵍ::Array{Float64, 1}`: Vector with the degree of mobility of each strata.
- `pᵍ_eff::Array{Float64, 1}`: Vector with the current degree of mobility of
  each strata.
- `edgelist::Array{Int64, 2}`: Matrix with the directed edgelist between
  patches, where ``L`` is the number of edges.
- `Rᵢⱼ::Array{Float64, 1}`: Vector with the transition probabilities for each
  edge in the edgelist.
- `mobilityᵍ::Array{Float64, 1}`: Effective mobility of each strata between
  patches (used to optimize performance).
- `kᵍ::Array{Float64, 1}`: Average number of contacts of each strata.
- `kᵍ_h::Array{Float64, 1}`: Average number of contacts at home of each strata.
- `kᵍ_w::Array{Float64, 1}`: Average number of contacts at work of each strata.
- `kᵍ_eff::Array{Float64, 1}`: Current average number of contacts of each
  strata.
- `C::Array{Float64, 2}`: Matrix with the probability of contact between
  different stratifications.
- `zᵍ::Array{Float64, 1}`: Nomalization factor each strata.
- `normᵍ::Array{Float64, 2}`: Normalization of each strata (used to optimize
  performance).
- `sᵢ::Array{Float64, 1}`: Surface of each patch.
- `ξ::Float64`: Densisty factor.
- `σ::Float64`: Average household size.

# Constructor

Use outer constructor [`Population_Params`](#MMCAcovid19.Population_Params-Tuple{Int64,Int64,Array{Float64,2},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,2},Array{Float64,1},Array{Int64,2},Array{Float64,1},Array{Float64,1},Float64,Float64})
for a proper initialization of this struct.
"""
struct Population_Params
    G::Int64
    M::Int64
    nᵢ::Array{Float64, 1}
    nᵢᵍ::Array{Float64, 2}
    nᵢ_eff::Array{Float64, 1}
    nᵢᵍ_eff::Array{Float64, 2}
    N::Int64
    Nᵍ::Array{Int64, 1}
    pᵍ::Array{Float64, 1}
    pᵍ_eff::Array{Float64, 1}
    edgelist::Array{Int64, 2}
    Rᵢⱼ::Array{Float64, 1}
    mobilityᵍ::Array{Float64, 2}
    kᵍ::Array{Float64, 1}
    kᵍ_h::Array{Float64, 1}
    kᵍ_w::Array{Float64, 1}
    kᵍ_eff::Array{Float64, 1}
    C::Array{Float64, 2}
    zᵍ::Array{Float64, 1}
    normᵍ::Array{Float64, 2}
    sᵢ::Array{Float64, 1}
    ξ::Float64
    σ::Float64
end


"""
    Population_Params(G::Int64,
                      M::Int64,
                      nᵢᵍ::Array{Float64, 2},
                      kᵍ::Array{Float64, 1},
                      kᵍ_h::Array{Float64, 1},
                      kᵍ_w::Array{Float64, 1},
                      C::Array{Float64, 2},
                      pᵍ::Array{Float64, 1},
                      edgelist::Array{Int64, 2},
                      Rᵢⱼ::Array{Float64, 1},
                      sᵢ::Array{Float64, 1},
                      ξ::Float64,
                      σ::Float64)

Constructor of the struct [`Population_Params`](@ref).

# Arguments

- `G::Int64`: Number of population strata.
- `M::Int64`: Number of patches.
- `nᵢᵍ::Array{Float64, 2}`: Matrix of size ``G \\times M`` with the population
  at each strata and patch.
- `kᵍ::Array{Float64, 1}`: Vector of size ``G`` with the average number of
  contacts of each strata.
- `kᵍ_h::Array{Float64, 1}`: Vector of size ``G`` with the average number of
  contacts at home of each strata.
- `kᵍ_w::Array{Float64, 1}`: Vector of size ``G`` with the average number of
  contacts at work of each strata.
- `C::Array{Float64, 2}`: Matrix of size ``G \\times G`` with the probability of
  contact between different stratifications.
- `pᵍ::Array{Float64, 1}`: Vector of size ``G`` with the degree of mobility of
  each strata ranged between 0 and 1.
- `edgelist::Array{Int64, 2}`: Matrix of size ``L \\times 2`` containing the
  directed edgelist between patches, where ``L`` is the number of edges. The IDs
  of the patches have to go from ``1`` to ``M``.
- `Rᵢⱼ::Array{Float64, 1}`: Vector of size ``L`` containing the transition
  probabilities for each edge in the edgelist.
- `sᵢ::Array{Float64, 1}`: Vector of size `M` with the surface of each patch.
- `ξ::Float64`: Density factor.
- `σ::Float64`: Average household size.

# Return

Struct that contains the parameters related with geographical, population and
mobility data.
"""
function Population_Params(G::Int64,
                           M::Int64,
                           nᵢᵍ::Array{Float64, 2},
                           kᵍ::Array{Float64, 1},
                           kᵍ_h::Array{Float64, 1},
                           kᵍ_w::Array{Float64, 1},
                           C::Array{Float64, 2},
                           pᵍ::Array{Float64, 1},
                           edgelist::Array{Int64, 2},
                           Rᵢⱼ::Array{Float64, 1},
                           sᵢ::Array{Float64, 1},
                           ξ::Float64,
                           σ::Float64)

    edgelist, Rᵢⱼ = correct_self_loops(edgelist, Rᵢⱼ, M)

    # Aggregate population count
    Nᵍ = round.(Int, sum(nᵢᵍ, dims = 2)[:, 1])
    N = sum(Nᵍ)
    nᵢ = sum(nᵢᵍ, dims = 1)[1, :]
    mobilityᵍ = zeros(Float64, G, length(Rᵢⱼ))

    # Init. effective population
    nᵢ_eff = zeros(Float64, M)
    nᵢᵍ_eff = (1 .- pᵍ) .* nᵢᵍ

    # Init. normalization vector
    zᵍ = zeros(Float64, G)
    normᵍ = zeros(Float64, G, M)

    # Compute effective population
    compute_effective_population!(nᵢᵍ_eff, nᵢ_eff, nᵢᵍ, Nᵍ, mobilityᵍ, kᵍ, zᵍ,
                                  normᵍ, ξ, pᵍ, sᵢ, edgelist, Rᵢⱼ, M, G)

    return Population_Params(G, M, nᵢ, copy(nᵢᵍ), nᵢ_eff, nᵢᵍ_eff, N, Nᵍ,
                             copy(pᵍ), copy(pᵍ), edgelist, Rᵢⱼ, mobilityᵍ,
                             copy(kᵍ), copy(kᵍ_h), copy(kᵍ_w), copy(kᵍ), copy(C),
                             zᵍ, normᵍ, copy(sᵢ), ξ, σ)
end


"""
    update_population_params!(population::Population_Params)

Update population parameters, computing the effective populations and the
normalization parameter z if p and k are modified.

# Arguments
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
"""
function update_population_params!(population::Population_Params)

    # Reset effective population
    population.nᵢᵍ_eff[:,:] .= (1 .- population.pᵍ_eff) .* population.nᵢᵍ
    population.nᵢ_eff[:] .= zeros(Float64, population.M)
    population.mobilityᵍ[:, :] .= 0.

    # Init. normalization vector
    population.zᵍ .= 0.
    population.normᵍ[:, :] .= 0.

    # Compute effective population
    compute_effective_population!(population.nᵢᵍ_eff,
                                  population.nᵢ_eff,
                                  population.nᵢᵍ,
                                  population.Nᵍ,
                                  population.mobilityᵍ,
                                  population.kᵍ_eff,
                                  population.zᵍ,
                                  population.normᵍ,
                                  population.ξ,
                                  population.pᵍ_eff,
                                  population.sᵢ,
                                  population.edgelist,
                                  population.Rᵢⱼ,
                                  population.M,
                                  population.G)
end


"""
    compute_effective_population!(nᵢᵍ_eff::Array{Float64, 2},
                                  nᵢ_eff::Array{Float64, 1},
                                  nᵢᵍ::Array{Float64, 2},
                                  Nᵍ::Array{Int64, 1},
                                  mobilityᵍ::Array{Float64, 2},
                                  kᵍ_eff::Array{Float64, 1},
                                  zᵍ::Array{Float64, 1},
                                  normᵍ::Array{Float64, 2},
                                  ξ::Float64,
                                  pᵍ_eff::Array{Float64, 1},
                                  sᵢ::Array{Float64, 1},
                                  edgelist::Array{Int64, 2},
                                  Rᵢⱼ::Array{Float64, 1},
                                  M::Int64,
                                  G::Int64)

Compute the effective population at each patch.
"""
function compute_effective_population!(nᵢᵍ_eff::Array{Float64, 2},
                                       nᵢ_eff::Array{Float64, 1},
                                       nᵢᵍ::Array{Float64, 2},
                                       Nᵍ::Array{Int64, 1},
                                       mobilityᵍ::Array{Float64, 2},
                                       kᵍ_eff::Array{Float64, 1},
                                       zᵍ::Array{Float64, 1},
                                       normᵍ::Array{Float64, 2},
                                       ξ::Float64,
                                       pᵍ_eff::Array{Float64, 1},
                                       sᵢ::Array{Float64, 1},
                                       edgelist::Array{Int64, 2},
                                       Rᵢⱼ::Array{Float64, 1},
                                       M::Int64,
                                       G::Int64)

    # Compute the effective population for each strata and patch
    for indx_e in 1:size(edgelist)[1]
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        for g in 1:G
            nᵢᵍ_eff[g, j] += pᵍ_eff[g] * Rᵢⱼ[indx_e] * nᵢᵍ[g, i]
        end
    end

    # Compute the aggerated effective populatoin
    for i in 1:M
        for g in 1:G
            nᵢ_eff[i] += nᵢᵍ_eff[g, i]
        end
    end

    # Correction for populations without effective population
    for i in 1:M
        if nᵢ_eff[i] == 0.
            nᵢ_eff[i] = 1e-7
        end
        for g in 1:G
            if nᵢᵍ_eff[g, i] == 0.
                nᵢᵍ_eff[g, i] = 1e-7
            end
        end
    end

    # Compute the normalization vector
    for g in 1:G
        zᵍ[g] = kᵍ_eff[g] * Nᵍ[g] /
            sum((2 .- exp.(-ξ .* nᵢ_eff ./ sᵢ)) .* nᵢᵍ_eff[g, :]);
    end

    # Update the precomuted matrices
    for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]
        for g in 1:G
            mobilityᵍ[g, indx_e] = nᵢᵍ[g, i] *
                ((1 - pᵍ_eff[g]) * (i == j ? 1. : 0.) + pᵍ_eff[g] * Rᵢⱼ[indx_e] )
        end
    end

    for i in 1:M
        for g in 1:G
            normᵍ[g, i] = zᵍ[g] * (2 - exp(-ξ * nᵢ_eff[i] / sᵢ[i]))
        end
    end
end


"""
    set_initial_conditions!(epi_params::Epidemic_Params,
                               population::Population_Params,
                               Sᵛ₀::Array{Float64, 2},
                               E₀::Array{Float64, 2},
                               A₀::Array{Float64, 2},
                               I₀::Array{Float64, 2},
                               H₀::Array{Float64, 2},
                               R₀::Array{Float64, 2})

Set the initial number of individuals on a population. They can be
introduced as susceptible vaccinated (Sᵛ₀), exposed unvaccinated (E₀), asymptomatic 
unvaccinated (A₀), symptomatic unvaccinated (I₀), hospitalized unvaccinated (H₀) or recovered 
unvaccinated (R₀) individuals.

# Arguments
- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
- `Sᵛ₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
  of vaccinated susceptible individuals of each strata on each patch.
- `E₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
  of unvaccinated exposed individuals of each strata on each patch.
- `A₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
  of unvaccinated asymptomatic infected individuals of each strata on each patch.
- `I₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
  of unvaccinated symptomatic infected individuals of each strata on each patch.
- `H₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
  of unvaccinated pre-hospitalized individuals of each strata on each patch.
- `R₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
  of unvaccinated recovered individuals of each strata on each patch.
"""
function set_initial_conditions!(epi_params::Epidemic_Params,
                               population::Population_Params,
                               Sᵛ₀::Array{Float64, 2},
                               E₀::Array{Float64, 2},
                               A₀::Array{Float64, 2},
                               I₀::Array{Float64, 2},
                               H₀::Array{Float64, 2},
                               R₀::Array{Float64, 2})

    t₀ = 1
    
    # Initial susceptible population
    @. epi_params.ρˢᵍᵥ[:, :, t₀, 2] = Sᵛ₀ / population.nᵢᵍ

    # Initial exposed population
    @. epi_params.ρᴱᵍᵥ[:, :, t₀, 1] = E₀ / population.nᵢᵍ 

    # Initial asymptomatic population
    @. epi_params.ρᴬᵍᵥ[:, :, t₀, 1] = A₀ / population.nᵢᵍ 

    # Initial sypmtomatic population
    @. epi_params.ρᴵᵍᵥ[:, :, t₀, 1] = I₀ / population.nᵢᵍ 
    
    # Initial hospitalized population
    @. epi_params.ρᴾᴴᵍᵥ[:, :, t₀, 1] = H₀ / population.nᵢᵍ 
    
    # Initial recovered population
    @. epi_params.ρᴿᵍᵥ[:, :, t₀, 1] = R₀ / population.nᵢᵍ

    # Control over division by zero
    epi_params.ρˢᵍᵥ[isnan.(epi_params.ρˢᵍᵥ)] .= 0
    epi_params.ρᴱᵍᵥ[isnan.(epi_params.ρᴱᵍᵥ)] .= 0
    epi_params.ρᴬᵍᵥ[isnan.(epi_params.ρᴬᵍᵥ)] .= 0
    epi_params.ρᴵᵍᵥ[isnan.(epi_params.ρᴵᵍᵥ)] .= 0
    epi_params.ρᴾᴴᵍᵥ[isnan.(epi_params.ρᴾᴴᵍᵥ)] .= 0
    epi_params.ρᴿᵍᵥ[isnan.(epi_params.ρᴿᵍᵥ)] .= 0

    # Update the fraction of suceptible individual
    @. epi_params.ρˢᵍᵥ[:, :, t₀, 1] = 1 - (epi_params.ρˢᵍᵥ[:, :, t₀, 2] +
                                       epi_params.ρᴱᵍᵥ[:, :, t₀, 1] +
                                       epi_params.ρᴬᵍᵥ[:, :, t₀, 1] +
                                       epi_params.ρᴵᵍᵥ[:, :, t₀, 1] +
                                       epi_params.ρᴾᴴᵍᵥ[:, :, t₀, 1] +
                                       epi_params.ρᴿᵍᵥ[:, :, t₀, 1])
    
    @. epi_params.ρˢᵍᵥ[abs(epi_params.ρˢᵍᵥ) < 1e-12]  = 0
    
    nothing
end


"""
    reset_params!(epi_params::Epidemic_Params,
                  population::Population_Params)

Reset the epidemic and population parameters to reuse the structure and avoid
additional data allocations.

# Arguments
- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
"""
function reset_params!(epi_params::Epidemic_Params,
                       population::Population_Params)

    reset_epidemic_params!(epi_params)
    population.pᵍ_eff .= population.pᵍ
    population.kᵍ_eff .= population.kᵍ
    update_population_params!(population)

end


### ----------------------------------------------------------------------------
### COMPUTE R FUNCTIONS
### ----------------------------------------------------------------------------

"""
    compute_R_eff(epi_params::Epidemic_Params,
                  population::Population_Params;
                  τ::Int64 = 21)

Compute the effective reproduction number R.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.

# Optional

- `τ::Int64 = 21`: kernel length.

# Return

Tuple formed by:
- DataFrame containing information about the evolution of the effective
  reproduction number R for each strata and patch.
- DataFrame containing information about the evolution of the total effective
  reproduction number R.
"""
function compute_R_eff(epi_params::Epidemic_Params,
                       population::Population_Params,
                       τ::Int64 = 21)

    M = population.M
    G = population.G
    T = epi_params.T - 1 # because we have to cut out the case t=1
    V = epi_params.V
    
    βᴬ = epi_params.βᴬ[1]
    βᴵ = epi_params.βᴵ[1]

    ηᵍ = epi_params.ηᵍ
    αᵍ = epi_params.αᵍ
    μᵍ = epi_params.μᵍ
    
    rᵥ = epi_params.rᵥ
    kᵥ = epi_params.kᵥ
    
    ρˢᵍᵥ = epi_params.ρˢᵍᵥ

    # Setup the kernel
    tʷ = 0:(τ - 1)
    wᵍᵥ = ones(Float64, G, τ, V)

    # Initialize results
    Rᵢᵍ_eff = DataFrame()
    Rᵢᵍ_eff.strata = repeat(1:G, outer = (T - τ) * M * V)
    Rᵢᵍ_eff.vaccine = repeat(1:V, inner = G, outer = (T - τ) * M)
    Rᵢᵍ_eff.patch = repeat(1:M, inner = G * V, outer = (T - τ))
    Rᵢᵍ_eff.time = repeat(1:(T - τ), inner = G * M * V)

    R_eff = DataFrame()
    R_eff.time = 1:(T - τ)
    Δρ = ones(G, M, T - τ, V)

    # Compute R
    Rᵢᵍᵥ = zeros(Float64, G, M, (T - τ), V)
    Rᵢᵍ = ones(Float64, G, M, (T - τ) )
    Rᵍ = ones(Float64, G, (T - τ) )
    Rᵍᵥ= ones(Float64, G, (T - τ), V )
    Rᵥ = ones(Float64, (T - τ), V )
    R = zeros(T - τ)
    for t in 1:(T - τ)
        for v in 1:V
            # This calculation would hold only without reinfection and waning immunity
#             @simd for s in tʷ
#                 @. Rᵢᵍᵥ[:, :, t, v] += epi_params.Qᵢᵍ[:, :, t + s + 2] * wᵍᵥ[:, s + 1, v] 
#                 # the +2 instead of +1 is because time t corresponds to the time t+1 in the epi_params
#             end
            @simd for g in 1:G
                @. Rᵢᵍᵥ[g, :, t, v] = epi_params.Qᵢᵍ[g, :, t] * (1 - kᵥ[v]) * (βᴬ/αᵍ[g] + βᴵ/μᵍ[g] )
            end
        end

        # Notice that Δρ of time t corresponds to the time t+1 in the epi_params
        # @. Δρ[:, :, t, :] = ρˢᵍᵥ[:, :, t, :] - ρˢᵍᵥ[:, :, t + 1, :]
        # R[t] = sum(Rᵢᵍᵥ[:, :, t, :] .* Δρ[:, :, t, :] .* population.nᵢᵍ[:, :]) / sum( Δρ[:, :, t, :] .* population.nᵢᵍ[:, :] )
        
        Rᵍᵥ[:, t, :] = sum( Rᵢᵍᵥ[:, :, t, :] , dims = 2 )[:, 1, :] / M
        Rᵥ[t, :] = sum( Rᵍᵥ[:, t, :] , dims = 1 )[1, :] / G 
        Rᵢᵍ[:, :, t] = sum( Rᵢᵍᵥ[:, :, t, :] , dims = 4 )[:,:,1] / V
        Rᵍ[:, t] = sum( Rᵢᵍ[:, :, t] , dims = 2 )[:,1] / M 
        R[t] = sum(Rᵢᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :]) / population.N
        
    end

    Rᵢᵍ_eff.R_eff = reshape(Rᵢᵍᵥ, G * M * (T - τ) * V)
    R_eff.R_eff = R

    return Rᵢᵍ_eff, Rᵍ, R_eff, Rᵥ, Rᵍᵥ
end


### ----------------------------------------------------------------------------
### OUTPUT FUNCTIONS
### ----------------------------------------------------------------------------

"""
    store_compartment(epi_params::Epidemic_Params,
                      population::Population_Params,
                      compartment::Char,
                      suffix::String,
                      folder::String)

Store the evolution of the given epidemic compartment for each strata and patch.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
- `compartment::String`: String indicating the compartment, one of: `"S"`,
  `"E"`, `"A"`, `"I"`, `"PH"`, `"PD"`, `"HR"`, `"HD"`, `"D"`, `"R"`.
- `suffix::String`: String used to identify the experiment.
- `folder::String`: String containing the path to the folder where the results
  will be stored.
"""
function store_compartment(epi_params::Epidemic_Params,
                           population::Population_Params,
                           compartment::String,
                           suffix::String,
                           folder::String)

    M = population.M
    G = population.G
    T = epi_params.T
    V = epi_params.V

    # Init. dataframe
    df = DataFrame()
    df.strata = repeat(1:G, outer = T * M)
    df.patch = repeat(1:M, inner = G, outer = T)
    df.time = repeat(1:T, inner = G * M)

    # Store number of cases
    if compartment == "S"
        df.cases = reshape(epi_params.ρˢᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "E"
        df.cases = reshape(epi_params.ρᴱᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "A"
        df.cases = reshape(epi_params.ρᴬᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "I"
        df.cases = reshape(epi_params.ρᴵᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "PH"
        df.cases = reshape(epi_params.ρᴾᴴᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "PD"
        df.cases = reshape(epi_params.ρᴾᴰᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "HR"
        df.cases = reshape(epi_params.ρᴴᴿᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "HD"
        df.cases = reshape(epi_params.ρᴴᴰᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "D"
        df.cases = reshape(epi_params.ρᴰᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "R"
        df.cases = reshape(epi_params.ρᴿᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    end

    CSV.write(@sprintf("%s/output_%s_%s.csv", folder, compartment, suffix), df)
end


"""
    store_R_eff(epi_params::Epidemic_Params,
                population::Population_Params,
                suffix::String,
                folder::String;
                τ::Int64 = 21)

Compute and store the effective reproduction number R.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
- `suffix::String`: String used to identify the experiment.
- `folder::String`: String containing the path where the results will be stored.

# Optional

- `τ::Int64 = 21`: kernel length.
"""
function store_R_eff(epi_params::Epidemic_Params,
                     population::Population_Params,
                     suffix::String,
                     folder::String;
                     τ::Int64 = 21)

    Rᵢᵍ_eff, R_eff = compute_R_eff(epi_params, population, τ)

    # Write the results
    CSV.write(@sprintf("%s/output_Reff_%s.csv", folder, suffix), Rᵢᵍ_eff)
    CSV.write(@sprintf("%s/output_Reff_total_%s.csv", folder, suffix), R_eff)
end


### ----------------------------------------------------------------------------
### EXTRA FUNCTIONS
### ----------------------------------------------------------------------------

"""
    correct_self_loops(edgelist::Array{Int64, 2},
                       Rᵢⱼ::Array{Float64, 1},
                       M::Int64)

Repair weighted matrices to ensure that all nodes are represented adding an
additional self-loop.

# Arguments
- `edgelist::Array{Int64, 2}`: Matrix containing the directed edgelist between
  patches. The IDs of the patches have to go from 1 to M.
- `Rᵢⱼ::Array{Float64, 1}`: Vector containing the transition probabilities for
  each edge in the edgelist.
- `M::Int64`: Number of patches.

# Return

Returns the corrected list of edges and transition probabilities.
"""
function correct_self_loops(edgelist::Array{Int64, 2},
                            Rᵢⱼ::Array{Float64, 1},
                            M::Int64)

    self_loops = falses(M)
    k_in = zeros(Int64, M)
    k_out = zeros(Int64, M)

    for indx_e in 1:size(edgelist)[1]
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        k_in[j] += 1
        k_out[i] += 1

        if i == j
            self_loops[i] = true
        end
    end

    to_add = collect(1:M)[(k_out .== 0) .& .!(self_loops)]

    new_edges = zeros(Int64, length(to_add), 2)
    new_edges[:,1] = to_add
    new_edges[:,2] = to_add

    edgelist = vcat(edgelist, new_edges)
    Rᵢⱼ = vcat(Rᵢⱼ, ones(Float64, length(to_add)))

    to_add = collect(1:M)[(k_out .!= 0) .& .!(self_loops)]

    new_edges = zeros(Int64, length(to_add), 2)
    new_edges[:,1] = to_add
    new_edges[:,2] = to_add

    edgelist = vcat(edgelist, new_edges)
    Rᵢⱼ = vcat(Rᵢⱼ, zeros(Float64, length(to_add)))

    return (edgelist, Rᵢⱼ)
end

"""
    make_edls(Rᵢⱼ) 

Turns a origin destination matrix into an edgelist

# Arguments
- `Rᵢⱼ::Array{Float64, 2}`: Matrix containing the flows of people from and to each node
  of the networks

# Return

Returns an edgelist and a list of flows in those edges.
"""


function make_edls(Rᵢⱼ) 
    Rᵢⱼ = replace(Rᵢⱼ, NaN => 0)
    M = size(Rᵢⱼ)[1]
    edgl = ones( sum( Rᵢⱼ .!= 0 ), 2 )
    count = 1
    a = ones( sum(Rᵢⱼ .!= 0) )
    
    for i in 1:M
        @simd for j in 1:M
            if Rᵢⱼ[i, j] != 0
                edgl[count, 1] = i
                edgl[count, 2] = j
                a[count] = Rᵢⱼ[i, j]
                count += 1
            end
        end
    end  
    
    edgl = convert(Matrix{Int64}, edgl)
    
    return edgl, a
end
