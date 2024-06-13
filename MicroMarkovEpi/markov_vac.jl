using Logging

"""
    update_prob!(Pᵢᵍᵥ::Array{Float64, 3},
                 Sᵢᵍᵥ::Array{Float64, 3},
                 τᵢᵍᵥ::Array{Float64, 3},
                 epi_params::Epidemic_Params,
                 population::Population_Params,
                 κ₀::Float64,
                 ϕ::Float64,
                 δ::Float64,
                 ϵᵍ::Array{Float64, 1},
                 t::Int64,
                 tᶜ::Int64,
                 tᵛ::Int64)

Updates the probabilities of the model using the equations described in the
paper.
"""
function update_prob!(Pᵢᵍᵥ::Array{Float64, 3},
                      Sᵢᵍᵥ::Array{Float64, 3},
                      τᵢᵍᵥ::Array{Float64, 3},
                      epi_params::Epidemic_Params,
                      population::Population_Params,
                      κ₀::Float64,
                                 ϕ::Float64,
                      δ::Float64,
                     	         ϵᵍ::Array{Float64, 1},
                      t::Int64,
                      tᶜ::Int64, 
                      tᵛ::Int64)

    # Shortcuts to parameters
    ηᵍ = epi_params.ηᵍ
    αᵍ = epi_params.αᵍ
    μᵍ = epi_params.μᵍ
    θᵍ = epi_params.θᵍ
    γᵍ = epi_params.γᵍ
    ζᵍ = epi_params.ζᵍ
    λᵍ = epi_params.λᵍ
    ωᵍ = epi_params.ωᵍ
    ψᵍ = epi_params.ψᵍ
    χᵍ = epi_params.χᵍ
    rᵥ = epi_params.rᵥ
    kᵥ = epi_params.kᵥ
    Γ = epi_params.Γ
    Λ = epi_params.Λ
    
    # Shortcut to variables
    ρˢᵍᵥ = epi_params.ρˢᵍᵥ
    ρᴱᵍᵥ = epi_params.ρᴱᵍᵥ
    ρᴬᵍᵥ = epi_params.ρᴬᵍᵥ
    ρᴵᵍᵥ = epi_params.ρᴵᵍᵥ
    ρᴾᴴᵍᵥ = epi_params.ρᴾᴴᵍᵥ
    ρᴾᴰᵍᵥ = epi_params.ρᴾᴰᵍᵥ
    ρᴴᴿᵍᵥ = epi_params.ρᴴᴿᵍᵥ
    ρᴴᴰᵍᵥ = epi_params.ρᴴᴰᵍᵥ
    ρᴰᵍᵥ = epi_params.ρᴰᵍᵥ
    ρᴿᵍᵥ = epi_params.ρᴿᵍᵥ
    
    CHᵢᵍᵥ = epi_params.CHᵢᵍᵥ
    G = population.G
    M = population.M
    Nᵍ = population.Nᵍ
    nᵢᵍ = population.nᵢᵍ
    V = epi_params.V
    pᵍ_eff = population.pᵍ_eff
    C = population.C
    edgelist = population.edgelist
    Rᵢⱼ = population.Rᵢⱼ
    kᵍ_h = population.kᵍ_h
    kᵍ_hw = population.kᵍ_h .+ population.kᵍ_w
    
    # Total population
    N = sum(nᵢᵍ)

    # Intervention at time tᶜ
    if tᶜ == t
        pᵍ_eff[:] .= (1 - κ₀) .* population.pᵍ

        
        if (κ₀ != 0. )
            population.kᵍ_eff .= kᵍ_h * κ₀ .+ kᵍ_hw * (1 - δ) * (1 - κ₀)
            # elder keep home contacts during confinement
            population.kᵍ_eff[G] = kᵍ_h[G]
        end

        update_population_params!(population)
    end

    # Get P and compute Q
    compute_P!(Pᵢᵍᵥ, Sᵢᵍᵥ, pᵍ_eff, ρˢᵍᵥ, ρᴬᵍᵥ, ρᴵᵍᵥ,
               epi_params.Qᵢᵍ, population.nᵢᵍ_eff, population.mobilityᵍ,
               population.normᵍ, epi_params.βᴬ[1], epi_params.βᴵ[1],
               edgelist, Rᵢⱼ, C, M, G, V, t, rᵥ, kᵥ)
                
    
    # Compute τᵢᵍᵥ
    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]
        
        for v in 1:V
            @simd for g in 1:G
                τᵢᵍᵥ[g, i, v] += Rᵢⱼ[indx_e] * Pᵢᵍᵥ[g, j, v]
            end
        end
    end
    
    
    # Compute vaccine priority distribution among ages and patches
    
    ϵᵢᵍ = optimal_vaccination_distribution(ϵᵍ::Array{Float64, 1},
                                           ρˢᵍᵥ::Array{Float64, 4},
                                           nᵢᵍ::Array{Float64, 2},
                                           t::Int64)
            
    
    # Newly vaccinated people as a fraction of the subpopulation
    new_vaccinated = zeros(Float64, G, M)
    new_vaccinated .= ϵᵢᵍ ./ nᵢᵍ
    new_vaccinated[isnan.(new_vaccinated)] .= 0.
    @info "Total vaccination", sum(new_vaccinated)

    # Update probabilities
    @inbounds for i in 1:M

        # Compute secure households
        CHᵢ = 0.0
        if tᶜ == t
            for g in 1:G
                aux = 0.0
                @simd for v in 1:V
                    aux += ρᴱᵍᵥ[g, i, t, v] + ρᴬᵍᵥ[g, i, t, v] + ρᴵᵍᵥ[g, i, t, v] 
                end
                CHᵢ += ( 1 - aux ) * population.nᵢᵍ[g, i]
            end
            CHᵢ = (1 - ϕ) * κ₀ * (CHᵢ / population.nᵢ[i]) ^ population.σ
        end
      
        
        # Update compartmental probabilities
        for g in 1:G
                  
            @simd for v in 1:V
                
                if tᶜ == t
                    ρˢᵍᵥ[g, i, t, v] += CHᵢᵍᵥ[g, i, v]
                end  
                
                # Infection probability
                Πᵢᵍᵥ = (1 - pᵍ_eff[g]) * Pᵢᵍᵥ[g, i, v] + pᵍ_eff[g] * τᵢᵍᵥ[g, i, v]
                
                # Pier: this is an ugly fix to avoid the problem of getting more vaccines that susceptibles
                # TO DO: incorporate this condition in the function optimal_vaccination_distribution
                if (v == 1) & ( (Πᵢᵍᵥ * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] + new_vaccinated[g, i]) > ρˢᵍᵥ[g, i, t, v] )
                    # print([g, i, t, v]  , "\n" )
                    new_vaccinated[g, i] = ρˢᵍᵥ[g, i, t, v] - Πᵢᵍᵥ * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] - 1e-10
                end

                # Epidemic compartments, where all states of vaccination are present
                ρˢᵍᵥ[g, i, t + 1, v] = ( 1 - Πᵢᵍᵥ ) * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] +
                    new_vaccinated[g, i] * ( [-1, +1, 0][v] ) +
                    Λ * ( [0, -1, +1][v] ) * ρˢᵍᵥ[g, i, t, 2] +  # The term inside the parentheses works as a if-then clause
                    Γ * ( [0 , 0, +1][v] ) * (ρᴿᵍᵥ[g, i, t, 1] + ρᴿᵍᵥ[g, i, t, 2] + ρᴿᵍᵥ[g, i, t, 3])
                
                ρᴱᵍᵥ[g, i, t + 1, v] = (1 - ηᵍ[g]) * ρᴱᵍᵥ[g, i, t, v] +
                    Πᵢᵍᵥ * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] 
        
                ρᴬᵍᵥ[g, i, t + 1, v] = (1 - αᵍ[g]) * ρᴬᵍᵥ[g, i, t, v] +
                    ηᵍ[g] * ρᴱᵍᵥ[g, i, t, v]

                ρᴵᵍᵥ[g, i, t + 1, v] = (1 - μᵍ[g]) * ρᴵᵍᵥ[g, i, t, v] +
                    αᵍ[g] * ρᴬᵍᵥ[g, i, t, v]

                ρᴾᴴᵍᵥ[g, i, t + 1, v] = (1 - λᵍ[g]) * ρᴾᴴᵍᵥ[g, i, t, v] +
                    μᵍ[g] * (1 - θᵍ[g, v]) * γᵍ[g, v] *  ρᴵᵍᵥ[g, i, t, v]

                ρᴾᴰᵍᵥ[g, i, t + 1, v] = (1 - ζᵍ[g]) * ρᴾᴰᵍᵥ[g, i, t, v] +
                    μᵍ[g] * θᵍ[g, v] * ρᴵᵍᵥ[g, i, t, v]

                ρᴴᴿᵍᵥ[g, i, t + 1, v] = (1 - χᵍ[g]) * ρᴴᴿᵍᵥ[g, i, t, v] +
                    λᵍ[g] * (1 - ωᵍ[g, v] ) * ρᴾᴴᵍᵥ[g, i, t, v]

                ρᴴᴰᵍᵥ[g, i, t + 1, v] = (1 - ψᵍ[g]) * ρᴴᴰᵍᵥ[g, i, t, v] +
                    λᵍ[g] * ωᵍ[g, v] * ρᴾᴴᵍᵥ[g, i, t, v]

                ρᴿᵍᵥ[g, i, t + 1, v] = ρᴿᵍᵥ[g, i, t, v] + χᵍ[g] * ρᴴᴿᵍᵥ[g, i, t, v] +
                    μᵍ[g] * (1 - θᵍ[g, v]) * (1 - γᵍ[g, v]) * ρᴵᵍᵥ[g, i , t, v] -
                    Γ * ρᴿᵍᵥ[g, i, t, v] 

                ρᴰᵍᵥ[g, i, t + 1, v] = ρᴰᵍᵥ[g, i, t, v] + ζᵍ[g] * ρᴾᴰᵍᵥ[g, i, t, v] +
                    ψᵍ[g] * ρᴴᴰᵍᵥ[g, i, t, v]
                

                if tᶜ == t
                    aux = ρˢᵍᵥ[g, i, t, v]
                    ρˢᵍᵥ[g, i, t, v] -= CHᵢᵍᵥ[g, i, v] 
                    CHᵢᵍᵥ[g, i, v] = CHᵢ * aux
                end 
            end   
            
            # Reset values
            τᵢᵍᵥ[g, i, :] .= 0.
	    # this should be one, based on the intial value provided in run_epidemic_spreading_mmca
            Pᵢᵍᵥ[g, i, :] .= 0. 
        end
    end
    
end


"""
    compute_P!(Pᵢᵍᵥ::Array{Float64, 3},
               Sᵢᵍᵥ::Array{Float64, 3},
               pᵍ_eff::Array{Float64, 1},
               ρˢᵍᵥ::Array{Float64, 4},
               ρᴬᵍᵥ::Array{Float64, 4},
               ρᴵᵍᵥ::Array{Float64, 4},
               Qᵢᵍ::Array{Float64, 3},
               nᵢᵍ_eff::Array{Float64, 2},
               mobilityᵍ::Array{Float64, 2},
               normᵍ::Array{Float64, 2},
               βᴬ::Float64,
               βᴵ::Float64,
               edgelist::Array{Int64, 2},
               Rᵢⱼ::Array{Float64, 1},
               C::Array{Float64, 2},
               M::Int64,
               G::Int64,
               V::Int64,
               t::Int64,
               rᵥ::Array{Float64, 1},
               kᵥ::Array{Float64, 1})

Compute ``P_i^g_v(t)`` and ``Q_i^g(t)`` as described in the referenced paper. The first quantity is needed to compute the infection probability and the second for the basig reproduction number
"""

function compute_P!(Pᵢᵍᵥ::Array{Float64, 3},
                    Sᵢᵍᵥ::Array{Float64, 3},
                    pᵍ_eff::Array{Float64, 1},
                    ρˢᵍᵥ::Array{Float64, 4},
                    ρᴬᵍᵥ::Array{Float64, 4},
                    ρᴵᵍᵥ::Array{Float64, 4},
                    Qᵢᵍ::Array{Float64, 3},
                    nᵢᵍ_eff::Array{Float64, 2},
                    mobilityᵍ::Array{Float64, 2},
                    normᵍ::Array{Float64, 2},
                    βᴬ::Float64,
                    βᴵ::Float64,
                    edgelist::Array{Int64, 2},
                    Rᵢⱼ::Array{Float64, 1},
                    C::Array{Float64, 2},
                    M::Int64,
                    G::Int64,
                    V::Int64,
                    t::Int64,
                    rᵥ::Array{Float64, 1},
                    kᵥ::Array{Float64, 1})
    

    # Init. aux variables
    Sᵢᵍᵥ = zeros(G, M, V)
    Pᵢᴬᵍᵥ = zeros(G, M, V)
    Pᵢᴵᵍᵥ = zeros(G, M, V)
    nˢᵍᵥ_ij = zeros(V)
    
    nᴬᵍᵥ_ij = zeros(G, M, V)
    nᴵᵍᵥ_ij = zeros(G, M, V)

    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2] # i->j

        # Get effective S, A and I
        for g in 1:G
            nˢᵍᵥ_ij[:] .= ( ρˢᵍᵥ[g, i, t, :] .* (1 .- rᵥ[:]) ) * mobilityᵍ[g, indx_e]
            Sᵢᵍᵥ[g, j, :] .=  Sᵢᵍᵥ[g, j, :] .+ nˢᵍᵥ_ij[:] / nᵢᵍ_eff[g, j]
            @simd for h in 1:G
                nᴬᵍᵥ_ij[g, j, :] .= nᴬᵍᵥ_ij[g, j, :] .+ ρᴬᵍᵥ[h, i, t, :] .* (C[g, h] * mobilityᵍ[h, indx_e] ./ nᵢᵍ_eff[h, j] .* ones(V) )
                nᴵᵍᵥ_ij[g, j, :] .= nᴵᵍᵥ_ij[g, j, :] .+ ρᴵᵍᵥ[h, i, t, :] .* (C[g, h] * mobilityᵍ[h, indx_e] ./ nᵢᵍ_eff[h, j] .* ones(V) )
            end
        end
    end


    # @inbounds for indx_e in 1:length(Rᵢⱼ)
    #     i = edgelist[indx_e, 1]
    #     j = edgelist[indx_e, 2]

    #     # Get effective S, A and I
    #     for g in 1:G
    #         nˢᵍ_ij = ρˢᵍ[g, i, t] * mobilityᵍ[g, indx_e]
    #         Sᵢᵍ[g, j] +=  nˢᵍ_ij / nᵢᵍ_eff[g, j]
    #         @simd for h in 1:G
    #             nᴬᵍ_ij = ρᴬᵍ[h, i, t] * mobilityᵍ[h, indx_e]
    #             nᴵᵍ_ij = ρᴵᵍ[h, i, t] * mobilityᵍ[h, indx_e]
    #             Pᵢᴬᵍ[g, j] += C[g, h] * nᴬᵍ_ij / nᵢᵍ_eff[h, j]
    #             Pᵢᴵᵍ[g, j] += C[g, h] * nᴵᵍ_ij / nᵢᵍ_eff[h, j]
    #         end
    #     end
    # end

    # Get P and effective ρ
    @inbounds for i in 1:M
        for v in 1:V
            for g in 1:G
                aux = 1.0
                @simd for w in 1:V
                    aux = aux * (1 - βᴬ*(1 - rᵥ[v])*(1 - kᵥ[w]) )^(normᵍ[g, i] * nᴬᵍᵥ_ij[g, i, w]) *
                                (1 - βᴵ*(1 - rᵥ[v])*(1 - kᵥ[w]) )^(normᵍ[g, i] * nᴵᵍᵥ_ij[g, i, w]) 
                end
                Pᵢᵍᵥ[g, i, v] = 1 - aux
            end
        end
    end
    
    # Compute Q to get the effective R
    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2] # i->j
        
        for g in 1:G     
            for v in 1:V
                @simd for h in 1:G
                    Qᵢᵍ[g, i, t] += normᵍ[g, j] * C[g, h] * Sᵢᵍᵥ[h, j, v] *
                    (pᵍ_eff[g] * Rᵢⱼ[indx_e] + (1 - pᵍ_eff[g]) * (i == j ? 1. : 0.))
                end
            end
        end
    end

end


"""
    print_status(epi_params::Epidemic_Params,
                 population::Population_Params,
                 t::Int64)

Print the status of the epidemic spreading.
"""
function print_status(epi_params::Epidemic_Params,
                      population::Population_Params,
                      t::Int64)

    players  = sum((epi_params.ρˢᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, :] .+
                    epi_params.CHᵢᵍᵥ[:, :, :] ) .* population.nᵢᵍ[:, :])

    sus3 = sum((epi_params.ρˢᵍᵥ[:, :, t, 3] ) .* population.nᵢᵍ[:, :])

    infected = sum(epi_params.ρᴵᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :] .+
                   epi_params.ρᴬᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :])

    cases3    = sum((epi_params.ρᴾᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 3]) .* population.nᵢᵍ[:, :])

    icus     = sum((epi_params.ρᴴᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, :]) .* population.nᵢᵍ[:, :])

    deaths   = sum(epi_params.ρᴰᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :])

    vaccine1 = sum((epi_params.ρˢᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 1] .+
                    epi_params.CHᵢᵍᵥ[:, :, 1] ) .* population.nᵢᵍ[:, :]) / population.N

    vaccine2 = sum((epi_params.ρˢᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 2] .+
                    epi_params.CHᵢᵍᵥ[:, :, 2] ) .* population.nᵢᵍ[:, :]) / population.N

    vaccine3 = sum((epi_params.ρˢᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.CHᵢᵍᵥ[:, :, 3] ) .* population.nᵢᵍ[:, :]) / population.N

    @printf("Time: %d, players: %.2f, sus3: %.2f, cases3: %.2f, deaths: %.2f, vaccine1 = %.2f, vaccine2: %.2f, vaccine3: %.2f\n",
            t, players, sus3, cases3, deaths, vaccine1, vaccine2, vaccine3 )
    
end

"""
Get the time series of different indicators of the epidemic:
 - Infected = I + A
 - Cases = I + A + PD + PH + HD + HR
 - Icus = HR + HD
 - Deaths = D
 - Vaccinated = All the vaccinated compartments
 - Daily cases
"""

function time_series(epi_params::Epidemic_Params,
                      population::Population_Params)
    
    return (infected = sum((epi_params.ρᴵᵍᵥ[:, :, :, :] .+ 
                            epi_params.ρᴬᵍᵥ[:, :, :, :]) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1],
        
            cases    = sum((epi_params.ρᴵᵍᵥ[:, :, :, :] .+ 
                            epi_params.ρᴬᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴾᴰᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴾᴴᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴴᴰᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴴᴿᵍᵥ[:, :, :, :]) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1],
        
            icus     = sum((epi_params.ρᴴᴿᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴴᴰᵍᵥ[:, :, :, :]) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1],
        
            deaths   = sum(epi_params.ρᴰᵍᵥ[:, :, :, :] .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1],
        
            vaccinated = sum((epi_params.ρˢᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴾᴰᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴱᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴬᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴵᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴾᴴᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴴᴰᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴴᴿᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴿᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴰᵍᵥ[:, :, :, 2:3] ) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1]
        
        )
end    

"""
    Function whose purpose is to find all the local maxima in a vector. The output is a list of two vectors, one containing all the heights 
    of said maxima, the other containing all the positions
"""

function maxima(v)
    h = []
    p = []
    idx = 1
    for i in 2:(length(v)-2)
        condition = (v[i-1] < v[i]) & (v[i] > v[i+1])
        # condition2 = (v[maximum([i-50, 1])] < v[i]) | (v[i] > v[minimum([length(v), i+50])])
        if condition #& condition2
           h = append!( h, v[i] )
           p = append!( p, i)
           idx = idx + 1 
        end
    end
    return (height = h ,
            position = p)
end


"""
    optimal_vaccination_distribution(ϵᵍ::Array{Float64, 1},
                                     ρˢᵍᵥ::Array{Float64, 4},
                                     nᵢᵍ::Array{Float64, 2},
                                     t::Int64)

Computes the number of vaccines that should be distributed among spatial patches and 
age strata taking into account where the majority of the unvaccinated susceptibles 
are and what age stratum should have the priority.

# Arguments

- ϵᵍ : the absolute number of vaccines at our disposal per day, divided among age strata
- ρˢᵍᵥ: fraction of susceptible individuals
- nᵢᵍ: absolute number of people per age in each age strata (rows) and patches (colums)
- t: current day

"""

function optimal_vaccination_distribution(ϵᵍ::Array{Float64, 1},
                                          ρˢᵍᵥ::Array{Float64, 4},
                                          nᵢᵍ::Array{Float64, 2},
                                          t::Int64)
    # print(t, "\n")
    Nᵥ = sum(ϵᵍ) # Total number of vaccines
    (G, M) = size(nᵢᵍ)

    if Nᵥ == 0
        @info "No vaccination"
        return zeros(G, M)
    end
    
    
    only_positive = true
    for v in 1:V
        only_positive = only_positive & all(ρˢᵍᵥ[:, :, t, v] .>= 0.0) & all(ρˢᵍᵥ[:, :, t, v] .<= 1.0)
    end
        
    if any(ϵᵍ .< 0)
        @printf("\n ----------------------------- \n ATTENZIONE: Number of dosis is negative \n ----------------------------- ")
        return
    end
    
    if !only_positive
        @printf("ATTENZIONE: Fraction of susceptible is not between 0 and 1")
        return
    end
    
    ###############################
    
    if ( ( Nᵥ != 0) & (sum(ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ) > Nᵥ) ) 
         
        # Define a matrix that gives you the priority of each patch and age group...
        priority_ϵ =  nᵢᵍ .* ( reshape(repeat(ϵᵍ, M), (G,M) ) )
        priority_ϵ = priority_ϵ / (sum(priority_ϵ) == 0 ? 1 : sum(priority_ϵ) )
        # ... and use the priority matrix to define how many dosis each subgroup get
        ϵᵢᵍ = Nᵥ * priority_ϵ
        
        # Define index that tells you if and where there are more susceptibles than vaccines
        idx_ϵ = ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ .- ϵᵢᵍ

        # Redistribution of vaccines: If there is one location and age group that has more vaccines than susceptibles the number 
        # of vaccines in that compartment is set equal to the number of susceptibles and the spare dosis are restributed among the others
        while ( !prod(idx_ϵ .>= 0 ) )   
            ϵᵢᵍ = ϵᵢᵍ .* (idx_ϵ .> 0)
            ϵᵢᵍ .* (idx_ϵ .<= 0) .= ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ .* (idx_ϵ .<= 0)
            Nᵥ_new = Nᵥ - sum(nᵢᵍ .* (idx_ϵ .<= 0) )
            
            # Redifine priority levels
            priority_ϵ =  nᵢᵍ .* ϵᵢᵍ 
            priority_ϵ = priority_ϵ / (sum(priority_ϵ) == 0 ? 1 : sum(priority_ϵ) )
            ϵᵢᵍ .* (idx_ϵ .> 0) .= Nᵥ_new * priority_ϵ
            
#             ϵᵢᵍ = ϵᵢᵍ / sum(ϵᵢᵍ) * ( Nᵥ - sum(nᵢᵍ[idx_ϵ .<= 0])  )
            
            idx_ϵ = ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ  .-  ϵᵢᵍ[:, :] 
        end
        
    elseif ( (Nᵥ != 0) & (sum(ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ) <= Nᵥ)  )
        ϵᵢᵍ = ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ
    else
        ϵᵢᵍ = zeros(G, M)
    end
    
    return ϵᵢᵍ

end


"""
    run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                 population::Population_Params;
                                 tᶜ::Int64 = -1,
                                 tᵛ::Int64 = -1,
                                 κ₀::Float64 = 0.0,
                                 ϕ::Float64 = 1.0,
                                 δ::Float64 = 0.0,
                                 ϵᵍ::Array{Float64, 1} = [0., 0., 0.],
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the
variables stored in epi_params. It also provides, through optional arguments,
the application of a containmnet or a vaccination campaign on a specific date.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.

## Optional

- `tᶜ::Int64 = -1`: Timestep of application of containment, or out of timesteps range
  value for no containment.
- `tᵛ::Int64 = -1`: Timestep of application of vaccination.
- `κ⁰::Float64 = 0.0`: Mobility reduction.
- `ϕ::Float64 = 1.0`: Permeability of confined households.
- `δ::Float64 = 0.0`: Social Distancing.
- `ϵᵍ::Array{Float64, 1} = [0., 0., 0.]`: Number of vaccines for each age group
- `t₀::Int64 = 1`: Initial timestep.
- `verbose::Bool = false`: If `true`, prints useful information about the
  evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params,
                                      npi_params::NPI_Params;
                                      tᵛ::Int64 = -1,
                                                         ϵᵍ::Array{Float64, 1} = [0., 0., 0.],
                                      t₀::Int64 = 1,
                                      verbose::Bool = false)
    G = population.G
    M = population.M
    V = epi_params.V

    # Initialize τᵢ (Π = (1 - p) P + pτ) and Pᵢ for markov chain
    τᵢᵍᵥ = zeros(Float64, G, M, V)
    Pᵢᵍᵥ = zeros(Float64, G, M, V)

    # Auxiliar arrays to compute P (avoid the allocation of additional memory)
    Sᵢᵍᵥ = zeros(Float64, G, M, V)
    
    
    run_epidemic_spreading_mmca!(epi_params, population, npi_params, [tᵛ], reshape(ϵᵍ, (3,1)) , t₀ = t₀, verbose = verbose)
end


"""
    run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                 population::Population_Params,
                                 tᶜs::Array{Int64, 1},
                                 tᵛs::Array{Int64, 1},
                                 κ₀s::Array{Float64, 1},
                                 ϕs::Array{Float64, 1},
                                 δs::Array{Float64, 1},
                                 ϵᵍs::Array{Float64, 2};
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the
variables stored in epi_params. It provides the option of the application
of multiple different containmnets or vaccination campaigns at specific dates.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
- `tᶜs::Array{Int64, 1}`: List of timesteps of application of containments.
- `tᵛs::Int64 = -1`: Timestep of application of vaccination.
- `κ⁰s::Array{Float64, 1}`: List of mobility reductions.
- `ϕs::Array{Float64, 1}`: List of permeabilities of confined households.
- `δs::Array{Float64, 1}`: List of social distancings.
- `ϵᵍs::Array{Float64, 2}`: List of dosis per age group of each time period

## Optional

- `t₀::Int64 = 1`: Initial timestep.
- `verbose::Bool = false`: If `true`, prints useful information about the
  evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params,
                                      npi_params::NPI_Params,
                                      tᵛs::Array{Int64, 1},
                                                         ϵᵍs::Array{Float64, 2};
                                      t₀::Int64 = 1,
                                      verbose::Bool = false)

    G = population.G
    M = population.M
    T = epi_params.T
    V = epi_params.V

    # Initialize τᵢ (Π = (1 - p) P + pτ) and Pᵢ for markov chain
    τᵢᵍᵥ = zeros(Float64, G, M, V)
    Pᵢᵍᵥ = ones(Float64, G, M, V)

    Sᵢᵍᵥ = zeros(Float64, G, M, V)

    κ₀s = npi_params.κ₀s
      ϕs = npi_params.ϕs
    δs = npi_params.δs
    tᶜs = npi_params.tᶜs

    # Initial state
    if verbose
        print_status(epi_params, population, t₀)
    end

    i = 1 # counter for containment
    j = 1 # counter for vaccination

    ## Start loop for time evoluiton
    @inbounds for t in t₀:(T - 1)
        
        update_prob!(Pᵢᵍᵥ, Sᵢᵍᵥ, τᵢᵍᵥ, epi_params, population,
                        κ₀s[i], ϕs[i], δs[i], ϵᵍs[:, j], t, tᶜs[i], tᵛs[j])
        
        if t == tᶜs[i] && i < length(tᶜs)
            i += 1
        end

        if t == tᵛs[j] && j < length(tᵛs)
            j += 1
        end
        
        # To avoid negative compartments
        only_positive = true
        for v in 1:V
            only_positive = only_positive & all(epi_params.ρˢᵍᵥ[:, :, t, v] .>= 0.0) & all(epi_params.ρˢᵍᵥ[:, :, t, v] .<= 1.0)
        end
        
        if !only_positive
            @printf("ATTENZIONE: I suscettibili sono meno di 0 o più di 1\n")
            return
        end

        if verbose
            print_status(epi_params, population, t + 1)
        end
    end
end
