## Test the datatypes module with all possible constructors and adders

function createpaulistring(nq)
    symbol = rand([:I, :X, :Y, :Z])
    qind = rand(1:nq)
    coeff = randn()
    PauliString(nq, symbol, qind, coeff)

    symbols = rand([:I, :X, :Y, :Z], min(nq, 4))
    qinds = rand(1:nq, min(nq, 4))
    coeff = randn()
    pstr = PauliString(nq, symbols, qinds, coeff)
    print(pstr)
    return pstr
end

function createpaulisum(nq)
    PauliSum(nq)

    pstr = createpaulistring(nq)
    PauliSum(nq, pstr)

    pstr = [createpaulistring(nq) for _ in 1:rand(1:4)]
    psum = PauliSum(nq, pstr)
    print(psum)
    return psum
end

function addtopaulisum(nq)
    psum = createpaulisum(nq)
    pstr = createpaulistring(nq)
    add!(psum, pstr)

    psum = createpaulisum(nq)
    pstr = [createpaulistring(nq) for _ in 1:rand(1:4)]
    add!(psum, pstr)

    symbol = rand([:I, :X, :Y, :Z])
    qind = rand(1:nq)
    coeff = randn()
    add!(psum, symbol, qind, coeff)

    symbols = rand([:I, :X, :Y, :Z], min(nq, 4))
    qinds = rand(1:nq, min(nq, 4))
    coeff = randn()
    psum2 = createpaulisum(nq)
    add!(psum2, symbols, qinds, coeff)
    print(psum)

    psum3 = add(psum, psum2)
    subtract(psum2, psum3)

    return psum
end