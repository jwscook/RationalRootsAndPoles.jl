using RationalRootsAndPoles
using Test
using BlackBoxOptim
using Random
Random.seed!(0)

function compareunsortedpairwise(a::Vector{Float64}, b; atol, rtol)
  @assert length(a) == 2length(b)
  allpasses = true
  for i in 1:2:length(a)
    thispass = false
    for j in b
      thispass |= isapprox(a[i:i+1], j; atol=atol, rtol=rtol)
    end
    allpasses &= thispass
  end
  return allpasses
end

@testset "RationalRootsAndPoles.jl" begin
  nsamples1D = 14
  for nroots = 1:3, npoles = 0:2
    @testset "$nroots roots, $npoles poles" begin
      for i in 1:10
          expectedroots = rand(2 * nroots)
          expectedpoles = rand(2 * npoles)
          objective(x) = RationalRootsAndPoles.rational(x, expectedroots, expectedpoles)
          res = RationalRootsAndPoles.solve(objective; lowerbounds=zeros(2), upperbounds=ones(2),
                                            nroots=nroots, npoles=npoles, timelimitpertry=5,
                                            nsamples1D=10, nretries=2)
          @test res.error < 1e-3
          @test compareunsortedpairwise(expectedroots, res.roots, rtol=1e-3, atol=0.0)
          @test compareunsortedpairwise(expectedpoles, res.poles, rtol=1e-3, atol=0.0)
      end
    end 
  end

end
