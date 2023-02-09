"""
    update_prob!(Pᵢᵍ::Array{Float64, 2},
                 Pᵢᴬᵍ::Array{Float64, 2},
                 Pᵢᴵᵍ::Array{Float64, 2},
                 Sᵢᵍᵥ::Array{Float64, 2},
                 τᵢᵍ::Array{Float64, 2},
                 epi_params::Epidemic_Params,
                 population::Population_Params,
                 κ₀::Float64,
                 ϕ::Float64,
                 δ::Float64,
                 t::Int64,
                 tᶜ::Int64)

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
                      tᶜ::Int64)

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
    
    CHᵢᵍ = epi_params.CHᵢᵍ
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
               population.normᵍ, epi_params.βᴬᵥ[1], epi_params.βᴵᵥ[1],
               edgelist, Rᵢⱼ, C, M, G, V, t, rᵥ, kᵥ)
                
    
    # Compute τᵢᵍᵥ
    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]
        
        for v in 1:V
            @simd for g in 1:G
                τᵢᵍᵥ[g, i, v] += Rᵢⱼ[indx_e] * Pᵢᵍᵥ[g, j, v]
                if τᵢᵍᵥ[g, i, v] > 1 
                    # @printf("%.5f \n", τᵢᵍᵥ[g, i, v] )
                end
            end
        end
    end
    
    # Compute vaccine priority distribution among ages and patches
    ϵᵢᵍ = optimal_vaccination_distribution(ϵᵍ::Array{Float64, 1},
                                           ρˢᵍᵥ::Array{Float64, 4},
                                           nᵢᵍ::Array{Float64, 2},
                                           t::Int64)

    # Update probabilities
    @inbounds for i in 1:M
        
        # Remove patches where you still have unvaccinated susceptibles
        # idx_ϵ = ρˢᵍᵥ[:, :, t, 1] .* nᵢᵍ[:, :]  .-  ϵᵢᵍ[:, :] 
        # ϵᵢᵍ = ϵᵢᵍ .* (idx_ϵ .> 0)

        # Compute secure households
        CHᵢ = 0.0
        if tᶜ == t
            for g in 1:G
                aux = 0.0
                @simd for v in 1:V
                    #CHᵢ += (ρˢᵍᵥ[g, i, t, 1]*(v == 1 ? 1 : 0) + ρᴾᴴᵍᵥ[g, i, t, v] + ρᴾᴰᵍᵥ[g, i, t, v] +
                    #        ρᴴᴿᵍᵥ[g, i, t, v] + ρᴴᴰᵍᵥ[g, i, t, v] + ρᴰᵍᵥ[g, i, t, v] +
                    #        ρᴿᵍᵥ[g, i, t, v] + CHᵢᵍ[g, i]*(v == 1 ? 1 : 0) ) * population.nᵢᵍ[g, i]
                    
                    aux += ρᴱᵍᵥ[g, i, t, v] + ρᴬᵍᵥ[g, i, t, v] + ρᴵᵍᵥ[g, i, t, v] 
                end
                CHᵢ += ( 1 - aux ) * population.nᵢᵍ[g, i]
            end
            CHᵢ = (1 - ϕ) * κ₀ * (CHᵢ / population.nᵢ[i]) ^ population.σ
        end
      
        
        # Update compartmental probabilities
        for g in 1:G
            if tᶜ == t
                # ATTENZIONE: Why did I put 1?
                ρˢᵍᵥ[g, i, t, 1] += CHᵢᵍ[g, i]
            end        
        
            @simd for v in 1:V

                Πᵢᵍᵥ = (1 - pᵍ_eff[g]) * Pᵢᵍᵥ[g, i, v] + pᵍ_eff[g] * τᵢᵍᵥ[g, i, v]

                # Epidemic compartments, where all states of vaccination are present
                ρˢᵍᵥ[g, i, t + 1, v] = ( 1 - Πᵢᵍᵥ ) * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] +
                    ϵᵢᵍ[g, i] * ( v == 1 ? -1 : 1 ) / nᵢᵍ[g, i] +
                    Λ * ( v == 1 ? 1 : -1 ) * ρˢᵍᵥ[g, i, t, 2] +
                    Γ * ( v == 2 ? 1 : 0 ) * ρᴿᵍᵥ[g, i, t, 1]
                
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
                    Γ * ( v == 1 ? 1 : 0 ) * ρᴿᵍᵥ[g, i, t, 1]

                ρᴰᵍᵥ[g, i, t + 1, v] = ρᴰᵍᵥ[g, i, t, v] + ζᵍ[g] * ρᴾᴰᵍᵥ[g, i, t, v] +
                    ψᵍ[g] * ρᴴᴰᵍᵥ[g, i, t, v]
            end

            if tᶜ == t
                aux = ρˢᵍᵥ[g, i, t, 1] + ρˢᵍᵥ[g, i, t, 2]
                ρˢᵍᵥ[g, i, t, 1] -= CHᵢᵍ[g, i] 
                CHᵢᵍ[g, i] = CHᵢ * aux
            end            
    
            # Reset values
            τᵢᵍᵥ[g, i, :] .= 0.
            Pᵢᵍᵥ[g, i, :] .= 0.
        end
    end
    
end


"""
    compute_P!(Pᵢᵍ::Array{Float64, 2},
                    Pᵢᴬᵍ::Array{Float64, 2},
                    Pᵢᴵᵍ::Array{Float64, 2},
                    Sᵢᵍᵥ::Array{Float64, 2},
                    pᵍ_eff::Array{Float64, 1},
                    ρˢᵍᵥ::Array{Float64, 3},
                    ρᴬᵍ::Array{Float64, 3},
                    ρᴵᵍ::Array{Float64, 3},
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
                    t::Int64)

Compute ``P_i^g(t)`` and ``Q_i^g(t)`` as described in the referenced paper.
"""
# function compute_P!(Pᵢᵍᵥ::Array{Float64, 3}, # +1 dimension for encoding the sigma
#                     Sᵢᵍᵥ::Array{Float64, 3},
#                     pᵍ_eff::Array{Float64, 1},
#                     ρˢᵍᵥ::Array{Float64, 4},
#                     ρᴬᵍᵥ::Array{Float64, 4},
#                     ρᴵᵍᵥ::Array{Float64, 4},
#                     Qᵢᵍ::Array{Float64, 3},
#                     nᵢᵍ_eff::Array{Float64, 2},
#                     mobilityᵍ::Array{Float64, 2},
#                     normᵍ::Array{Float64, 2},
#                     βᴬᵥ::Float64,
#                     βᴵᵥ::Float64,
#                     edgelist::Array{Int64, 2},
#                     Rᵢⱼ::Array{Float64, 1},
#                     C::Array{Float64, 2},
#                     M::Int64,
#                     G::Int64,
#                     V::Int64,
#                     t::Int64,
#                     rᵥ::Array{Float64, 1},
#                     kᵥ::Array{Float64, 1})

#     # Init. aux variables
#     Sᵢᵍᵥ .= 0 # zeros(G, M, V)
#     Pᵢᵍᵥ .= 1 # ones(G, M, V) # Necessary in order to use the cumulative product strategy
#     nˢᵍᵥ_ij = zeros(V)
    
#     # Trasmissibility reduction matrix
#     rk_ij = 0.
    
#     @inbounds for indx_e in 1:length(Rᵢⱼ)
#         i = edgelist[indx_e, 1]
#         j = edgelist[indx_e, 2]
#         # i -> j 

#         # Get effective S, A and I
#         for g in 1:G
#             # Necessary to get Q_i (not exactly equal to the n_ij of the paper)
#             # why was there .* ( 1 - kᵥ[:]) ?
#             nˢᵍᵥ_ij[:] .= ( ρˢᵍᵥ[g, i, t, :] .* (1 .- rᵥ[:]) ) * mobilityᵍ[g, indx_e]
#             Sᵢᵍᵥ[g, j, :] .= Sᵢᵍᵥ[g, j, :] .+ nˢᵍᵥ_ij[:] / nᵢᵍ_eff[g, j]

#             for v_i in 1:V
#                 for v_j in 1:V
#                     for h in 1:G
#                         nᴬᵍᵥ_ij = ρᴬᵍᵥ[h, i, t, v_i] * mobilityᵍ[h, indx_e]
#                         nᴵᵍᵥ_ij = ρᴵᵍᵥ[h, i, t, v_i] * mobilityᵍ[h, indx_e]
#                         rk_ij = (1-rᵥ[v_j])*(1-kᵥ[v_i])
#                         Pᵢᵍᵥ[g, j, v_j] *= (1 - βᴬᵥ * rk_ij )^(normᵍ[g, j] * C[g, h] * nᴬᵍᵥ_ij / nᵢᵍ_eff[h, j] ) * 
#                                           (1 - βᴵᵥ * rk_ij)^(normᵍ[g, j] * C[g, h] * nᴵᵍᵥ_ij / nᵢᵍ_eff[h, j] )
                        
#                         # the fact that we define Pᵢᵍᵥ with the j instead of the i is due to the fact that i->j
                        
#                     end
#                 end
#             end
#         end
#     end
    
#     # Get P
#     Pᵢᵍᵥ .= 1 .- Pᵢᵍᵥ

#     # Compute Q to get the effective R
#     @inbounds for indx_e in 1:length(Rᵢⱼ)
#         i = edgelist[indx_e, 1]
#         j = edgelist[indx_e, 2] # i->j
        
#         for g in 1:G     
#             for v in 1:V
#                 @simd for h in 1:G
#                     #### ATTENTION: Here I changes normᵍ[g, i] to normᵍ[g, j]
#                     Qᵢᵍ[g, i, t] += normᵍ[g, j] * C[g, h] * Sᵢᵍᵥ[h, j, v] *
#                     (pᵍ_eff[g] * Rᵢⱼ[indx_e] + (1 - pᵍ_eff[g]) * (i == j ? 1. : 0.))
#                 end
#             end
#         end
#     end
# end

##########àà VERSIONE TEMPORANEA DA ELIMINARE #############à

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
                    βᴬᵥ::Float64,
                    βᴵᵥ::Float64,
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
        j = edgelist[indx_e, 2]

        # Get effective S, A and I
        for g in 1:G
            nˢᵍᵥ_ij[:] .= ( ρˢᵍᵥ[g, i, t, :] .* (1 .- rᵥ[:]) ) * mobilityᵍ[g, indx_e]
            Sᵢᵍᵥ[g, j, :] .=  Sᵢᵍᵥ[g, j, :] .+ nˢᵍᵥ_ij[:] / nᵢᵍ_eff[g, j]
            @simd for h in 1:G
                nᴬᵍᵥ_ij[g, j, :] .= nᴬᵍᵥ_ij[g, j, :] .+ ρᴬᵍᵥ[h, i, t, :] .* (C[g, h] * mobilityᵍ[h, indx_e] ./ nᵢᵍ_eff[h, j] .* ones(V) )
                nᴵᵍᵥ_ij[g, j, :] .= nᴬᵍᵥ_ij[g, j, :] .+ ρᴵᵍᵥ[h, i, t, :] .* (C[g, h] * mobilityᵍ[h, indx_e] ./ nᵢᵍ_eff[h, j] .* ones(V) )
            end
        end
    end

    # Get P and effective ρ
    @inbounds for i in 1:M
        for v in 1:V
            @simd for g in 1:G
                Pᵢᵍᵥ[g, i, v] = 1 - 
                    (1 - βᴬᵥ*(1 - rᵥ[v])*(1 - kᵥ[1]) )^(normᵍ[g, i] * nᴬᵍᵥ_ij[g, i, 1]) *
                    (1 - βᴬᵥ*(1 - rᵥ[v])*(1 - kᵥ[2]) )^(normᵍ[g, i] * nᴬᵍᵥ_ij[g, i, 2]) *
                    (1 - βᴵᵥ*(1 - rᵥ[v])*(1 - kᵥ[1]) )^(normᵍ[g, i] * nᴵᵍᵥ_ij[g, i, 1]) *
                    (1 - βᴵᵥ*(1 - rᵥ[v])*(1 - kᵥ[2]) )^(normᵍ[g, i] * nᴵᵍᵥ_ij[g, i, 2])
            end
        end
    end
    
    ###########DA AGGIORNARE ASSOLUTAMENTE #############à
    # Compute Q to get the effective R
    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2] # i->j
        
        for g in 1:G     
            for v in 1:V
                @simd for h in 1:G
                    #### ATTENTION: Here I changes normᵍ[g, i] to normᵍ[g, j]
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
    
"""    
    s = sum(epi_params.ρˢᵍᵥ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    g0 = sum( epi_params.ρᴳ⁰ᵍ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    w1 = sum(epi_params.ρᵂ¹ᵍ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    g1 = sum(epi_params.ρᴳ¹ᵍ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    k2 = sum(epi_params.ρᴷ²ᵍ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    z0 = sum(epi_params.ρᶻ⁰ᵍ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    v1 = sum(epi_params.ρᵛ¹ᵍ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    z1 = sum(epi_params.ρᶻ¹ᵍ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    v2 = sum(epi_params.ρⱽ²ᵍ[:, :, t, 1] .* population.nᵢᵍ[:, :]) 
    
    @printf("Time: %d, s: %.2f, g0: %.2f, z0: %.2f, w1: %.2f, g1: %.2f, v1: %.2f, z1: %.2f, k2: %.2f, v2: %.2f\n",
            t, s, g0, z0, w1, g1, v1, z1, k2, v2)
    
"""

    players  = sum((epi_params.ρˢᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, :] ) .* population.nᵢᵍ[:, :]) +
               sum(epi_params.CHᵢᵍ[:, :] .* population.nᵢᵍ[:, :]) 

    infected = sum(epi_params.ρᴵᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :] .+
                   epi_params.ρᴬᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :])

    cases    = sum((epi_params.ρᴾᴰᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, :]) .* population.nᵢᵍ[:, :])

    icus     = sum((epi_params.ρᴴᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, :]) .* population.nᵢᵍ[:, :])

    deaths   = sum(epi_params.ρᴰᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :])

    
    vaccinated = sum((epi_params.ρˢᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 2] ) .* population.nᵢᵍ[:, :]) / population.N

    @printf("Time: %d, players: %.5f, icus: %.2f, deaths: %.2f, vaccine_check: %.3f\n",
            t, players, icus, deaths, vaccinated )
    
end

"""
Print the time series of different indicators of the epidemic
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
        
            vaccinated = sum((epi_params.ρˢᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴾᴰᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴱᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴬᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴵᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴾᴴᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴴᴰᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴴᴿᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴿᵍᵥ[:, :, :, 2] .+
                              epi_params.ρᴰᵍᵥ[:, :, :, 2] ) .* population.nᵢᵍ[:, :], dims=(1,2) )[1,1,:],
        
            daily_cases = diff( sum((epi_params.ρᴵᵍᵥ[:, :, :, :] .+ 
                                     epi_params.ρᴬᵍᵥ[:, :, :, :] .+
                                     epi_params.ρᴾᴰᵍᵥ[:, :, :, :] .+
                                     epi_params.ρᴾᴴᵍᵥ[:, :, :, :] .+
                                     epi_params.ρᴴᴰᵍᵥ[:, :, :, :] .+
                                     epi_params.ρᴴᴿᵍᵥ[:, :, :, :] .+
                                     epi_params.ρᴿᵍᵥ[:, :, :, :] .+
                                     epi_params.ρᴰᵍᵥ[:, :, :, :]) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1]
                                )
        
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
    
    Nᵥ = sum(ϵᵍ) # Total number of vaccines
    (G, M) = size(nᵢᵍ)
    
    ####### PROVA #########
    
    only_positive = all(ρˢᵍᵥ[:, :, t, 1] .>= 0.0) & 
            all(ρˢᵍᵥ[:, :, t, 2] .>= 0.0) &
            all(ρˢᵍᵥ[:, :, t, 1] .<= 1.0) &
            all(ρˢᵍᵥ[:, :, t, 2] .<= 1.0)
        
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
         
        # The idea is to 
        priority_ϵ =  nᵢᵍ .* ( reshape(repeat(ϵᵍ, M), (G,M) ) )
        priority_ϵ = priority_ϵ / (sum(priority_ϵ) == 0 ? 1 : sum(priority_ϵ) )
        ϵᵢᵍ = Nᵥ * priority_ϵ
        
        # Define index that tells you if and where there are more susceptibles than vaccines
        idx_ϵ = ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ .- ϵᵢᵍ
        
        # If there is one location and age group that has more vaccines than susceptibles the number of vaccines in
        # that compartment is set equal to the number of susceptibles
#         ϵᵢᵍ = ϵᵢᵍ .* (idx_ϵ .> 0)
#         ϵᵢᵍ.* (idx_ϵ .<= 0) .= ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ .* (idx_ϵ .<= 0)

        # Redistribution of vaccines (work in progress)
        while ( !prod(idx_ϵ .>= 0 ) )   
            ϵᵢᵍ .* (idx_ϵ .<= 0) .= ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ .* (idx_ϵ .<= 0)
            Nᵥ_new = Nᵥ - sum(nᵢᵍ .* (idx_ϵ .<= 0) )
            
            ϵᵢᵍ = ϵᵢᵍ .* (idx_ϵ .> 0)
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
                                 κ₀::Float64 = 0.0,
                                 ϕ::Float64 = 1.0,
                                 δ::Float64 = 0.0,
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the
variables stored in epi_params. It also provides, through optional arguments,
the application of a containmnet on a specific date.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.

## Optional

- `tᶜ::Int64 = -1`: Timestep of application of containment, or out of timesteps range
  value for no containment.
- `κ⁰::Float64 = 0.0`: Mobility reduction.
- `ϕ::Float64 = 1.0`: Permeability of confined households.
- `δ::Float64 = 0.0`: Social Distancing.
- `t₀::Int64 = 1`: Initial timestep.
- `verbose::Bool = false`: If `true`, prints useful information about the
  evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params;
                                      tᶜ::Int64 = -1,
                                      κ₀::Float64 = 0.0,
                                      ϕ::Float64 = 1.0,
                                      δ::Float64 = 0.0,
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
    # Pᵢᴬᵍ = zeros(Float64, G, M)
    # Pᵢᴵᵍ = zeros(Float64, G, M)
    Sᵢᵍᵥ = zeros(Float64, G, M, V)
    
    
    run_epidemic_spreading_mmca!(epi_params, population, [tᶜ],
                                 [κ₀], [ϕ], [δ], reshape(ϵᵍ, (3,1)) , t₀ = t₀, verbose = verbose)
end


"""
    run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                 population::Population_Params,
                                 tᶜs::Array{Int64, 1},
                                 κ₀s::Array{Float64, 1},
                                 ϕs::Array{Float64, 1},
                                 δs::Array{Float64, 1};
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the
variables stored in epi_params. It provides the option of the application
of multiple different containmnets at specific dates.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
- `tᶜs::Array{Int64, 1}`: List of timesteps of application of containments.
- `κ⁰s::Array{Float64, 1}`: List of mobility reductions.
- `ϕs::Array{Float64, 1}`: List of permeabilities of confined households.
- `δs::Array{Float64, 1}`: List of social distancings.

## Optional

- `t₀::Int64 = 1`: Initial timestep.
- `verbose::Bool = false`: If `true`, prints useful information about the
  evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params,
                                      tᶜs::Array{Int64, 1},
                                      κ₀s::Array{Float64, 1},
                                      ϕs::Array{Float64, 1},
                                      δs::Array{Float64, 1},
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

    # Initial state
    if verbose
        print_status(epi_params, population, t₀)
    end

    i = 1

    ## Start loop for time evoluiton
    @inbounds for t in t₀:(T - 1)

        
        update_prob!(Pᵢᵍᵥ, Sᵢᵍᵥ, τᵢᵍᵥ, epi_params, population,
                        κ₀s[i], ϕs[i], δs[i], ϵᵍs[:, i], t, tᶜs[i])
        
        if t == tᶜs[i] && i < length(tᶜs)
            i += 1
        end
        
        ######### AGGIUNTA PER EVITARE COMPARTIMENTI NEGATIVI ##########
        
        only_positive = all(epi_params.ρˢᵍᵥ[:, :, t, 1] .>= 0.0) & 
            all(epi_params.ρˢᵍᵥ[:, :, t, 2] .>= 0.0) &
            all(epi_params.ρˢᵍᵥ[:, :, t, 1] .<= 1.0) &
            all(epi_params.ρˢᵍᵥ[:, :, t, 2] .<= 1.0)
        
        if !only_positive
            @printf("ATTENZIONE: I suscettibili sono meno di 0 o più di 1")
            return
        end
            
        ########## FINE AGGIUNTA #################

        if verbose
            print_status(epi_params, population, t + 1)
        end
    end
end
