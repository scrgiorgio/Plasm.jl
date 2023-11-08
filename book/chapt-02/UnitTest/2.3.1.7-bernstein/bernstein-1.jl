B = n -> k -> u -> binomial(n,k) * u^k * (1-u)^(n-k)
# => #7 (generic function with 1 method)
Bernstein(n) = map(B(n), collect(0:n))
# => #7 Bernstein (generic function with 1 method)
Bernstein(2)
# => 3-element Vector{Function}
Bernstein(2)[2]
# => #9 (generic function with 1 method)
Bernstein(2)[2](0.1)
# => 0.18000000000000002
