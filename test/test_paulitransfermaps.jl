using LinearAlgebra
using Random
using Test

@testset "Unitaries PTM Tests" begin
    """Test the PTM for unitary matrices."""
    tol=1e-12

    # Test using single-qubit PauliRotation gate
    @testset "PauliRotation Y" begin
        pauligate = PauliRotation([:Y], [ 0])
        theta = Random.rand()

        udag = get_unitary_dagmat(pauligate, [theta])

        ptmmap = calc_ptm_dagmap(udag)

        expected_ptm = [
            [1 0 0 0];
            [0 cos(theta) 0 -sin(theta)];
            [0 0 1 0];
            [0 sin(theta) 0 cos(theta)]
        ]
       
       @test LinearAlgebra.norm(ptmmap - expected_ptm) < tol
    end

    # Test using T gate
    @testset "TUnitary" begin
        tgate = TUnitary([1])
        udag = get_unitary_dagmat(tgate, [])

        ptmmap = calc_ptm_dagmap(udag)

        expected_ptm = [
            [1 0 0 0];
            [0 1/sqrt(2) 1/sqrt(2) 0];
            [0 -1/sqrt(2) 1/sqrt(2) 0];
            [0 0 0 1]
        ]

        @test LinearAlgebra.norm(ptmmap - expected_ptm) < tol
    end

    #TODO (YT): add tests for two-qubit PauliRotation gates using QuEst.
end

@testset "PauliRotations PTM" begin

    function get_expected_ptm(nqubits, gate, param)
        """Construct the expected PTM from the PauliPropagation."""
        ptm_map_expected = Vector{Dict}()
        for i in 0:4^nqubits-1
            
            pstr = PauliString(nqubits, i, 1.)
            p_propagated = propagate([gate], pstr, [param]).terms

            push!(ptm_map_expected, p_propagated)
        end

        return ptm_map_expected

    end

    function get_ptm_map_dict(gate, theta, type)
        """Get the PTM map as a dictionary."""

        udag = get_unitary_dagmat(gate, [theta])
        ptm_mat = calc_ptm_dagmap(udag)
        ptm_map = get_ptm_sparse(ptm_mat, type)

        return [Dict(t[i] => t[i+1] for i in 1:2:length(t)) for t in ptm_map]
    end

    function is_ptm_expected(ptm_map_dict, ptm_map_expected)
        """Check if the PTM is consistent with the PauliPropagation."""

        for (t_actual, t_expected) in zip(ptm_map_dict, ptm_map_expected)

            # Check if the keys are the same
            @assert keys(t_actual) == keys(t_expected)

            for (k, v) in t_actual
                # Check the ptm map is the same
                @assert isapprox(v, t_expected[k])
            end
        end

        return true
    end    

    @testset "Single qubit rotation" begin
        nqubits = 1
        pauligate = PauliRotation([:Y], [1])
        theta = Random.rand() * 2 * pi

        ptm_map_dict = get_ptm_map_dict(pauligate, theta, getinttype(nqubits))

        # PauliPropagation
        ptm_map_expected = get_expected_ptm(nqubits, pauligate, theta)

        # Check if the PTM is consistent with the PauliPropagation
        @test is_ptm_expected(ptm_map_dict, ptm_map_expected)
    end

    @testset "Two qubit rotation" begin
        nqubits = 2
        pauligate = PauliRotation([:Y, :X], [1, 2])
        theta = Random.rand() * 2 * pi

        ptm_map_dict = get_ptm_map_dict(pauligate, theta, getinttype(nqubits))

        # PauliPropagation
        ptm_map_expected = get_expected_ptm(nqubits, pauligate, theta)

        # Check if the PTM is consistent with the PauliPropagation
        @test is_ptm_expected(ptm_map_dict, ptm_map_expected)
    end
    
end
