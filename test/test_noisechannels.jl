using Test


function depolarizing(channel_circ, damping_circ, thetas1, thetas2, where_ind, q_ind, noise_p)
    # add to the circuit
    insert!(channel_circ, where_ind, DepolarizingNoise(q_ind))

    insert!(damping_circ, where_ind, PauliPropagation.PauliYDamping(q_ind))
    insert!(damping_circ, where_ind, PauliPropagation.PauliXDamping(q_ind))
    insert!(damping_circ, where_ind, PauliPropagation.PauliZDamping(q_ind))

    # add to the parameters
    insert!(thetas1, where_ind, noise_p)

    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)

    return channel_circ, damping_circ, thetas1, thetas2
end

# this is also an alias for PauliZNoise
function dephasing(channel_circ, damping_circ, thetas1, thetas2, where_ind, q_ind, noise_p)
    # add to the circuit
    insert!(channel_circ, where_ind, DephasingNoise(q_ind))

    insert!(damping_circ, where_ind, PauliPropagation.PauliYDamping(q_ind))
    insert!(damping_circ, where_ind, PauliPropagation.PauliXDamping(q_ind))

    # add to the parameters
    insert!(thetas1, where_ind, noise_p)

    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)

    return channel_circ, damping_circ, thetas1, thetas2
end

function paulixnoise(channel_circ, damping_circ, thetas1, thetas2, where_ind, q_ind, noise_p)
    # add to the circuit
    insert!(channel_circ, where_ind, PauliXNoise(q_ind))

    insert!(damping_circ, where_ind, PauliPropagation.PauliYDamping(q_ind))
    insert!(damping_circ, where_ind, PauliPropagation.PauliZDamping(q_ind))

    # add to the parameters
    insert!(thetas1, where_ind, noise_p)

    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)

    return channel_circ, damping_circ, thetas1, thetas2
end

function pauliynoise(channel_circ, damping_circ, thetas1, thetas2, where_ind, q_ind, noise_p)
    # add to the circuit
    insert!(channel_circ, where_ind, PauliYNoise(q_ind))

    insert!(damping_circ, where_ind, PauliPropagation.PauliXDamping(q_ind))
    insert!(damping_circ, where_ind, PauliPropagation.PauliZDamping(q_ind))

    # add to the parameters
    insert!(thetas1, where_ind, noise_p)

    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)

    return channel_circ, damping_circ, thetas1, thetas2
end

@testset "Test Pauli Noises" begin

    builder_functions = [
        depolarizing,
        dephasing,
        paulixnoise,
        pauliynoise
    ]


    for builderfunc in builder_functions
        nq = rand(2:5)
        nl = 2
        W = Inf
        min_abs_coeff = 0.0

        pstr = PauliString(nq, rand([:X, :Y, :Z], nq), 1:nq)

        topo = bricklayertopology(nq; periodic=true)
        circ = hardwareefficientcircuit(nq, nl; topology=topo)

        m = countparameters(circ)

        channel_circ = deepcopy(circ)
        damping_circ = deepcopy(circ)

        thetas1 = randn(m)
        thetas2 = deepcopy(thetas1)

        where_ind = rand(1:m)
        q_ind = rand(1:nq)
        noise_p = rand() * 0.2

        channel_circ, damping_circ, thetas1, thetas2 = builderfunc(channel_circ, damping_circ, thetas1, thetas2, where_ind, q_ind, noise_p)



        dnum1 = propagate(channel_circ, pstr, thetas1; max_weight=W, min_abs_coeff=min_abs_coeff)

        dnum2 = propagate(damping_circ, pstr, thetas2; max_weight=W, min_abs_coeff=min_abs_coeff)

        @test overlapwithzero(dnum1) ≈ overlapwithzero(dnum2)
    end

end