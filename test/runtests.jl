using RationalRootsAndPoles
using Test
using BlackBoxOptim
using Random
Random.seed!(0)

function compareunsortedpairwise(a, b; atol, rtol)
  @assert length(a) == length(b)
  allpasses = true
  for i in 1:2:length(a)
    thispass = false
    for j in 1:2:length(b)
      thispass |= isapprox(a[i:i+1], b[j:j+1]; atol=atol, rtol=rtol)
    end
    allpasses &= thispass
  end
  return allpasses
end

@testset "RationalRootsAndPoles.jl" begin
  #complexroots(roots::Vector{<:Real}) = (r[1] + im * r[2] for r in eachcol(reshape(roots, 2, length(roots) ÷ 2)))
  #bar(x, roots) = prod(x[1] + im * x[2] - r for r in complexroots(roots); init=one(eltype(x)))
  #bar(x, ::Nothing) = 1
  #foo(x, roots, poles=nothing) = bar(x, roots) / bar(x, poles)

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
