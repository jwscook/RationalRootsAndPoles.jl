module RationalRootsAndPoles
using BlackBoxOptim


complexroots(roots::Vector{<:Real}) = (r[1] + im * r[2] for r in eachcol(reshape(roots, 2, length(roots) ÷ 2)))
numorden(x, roots) = prod(x[1] + im * x[2] - r for r in complexroots(roots); init=one(eltype(x)))
numorden(x, ::Nothing) = 1
rational(x, roots, poles=nothing) = numorden(x, roots) / numorden(x, poles)
onlyfinite(x::Number) = isfinite(x) ? x : zero(typeof(x))

function solve(f::Function; lowerbounds, upperbounds, nroots::Int, npoles::Int, nsamples1D::Int=14, nretries=0, atol=1e-3,
    timelimitpertry=5) where {F<:Function}
  x1D = 1/2nsamples1D:1/nsamples1D:1-1/2nsamples1D
  xys = [(x, y) .* (upperbounds .- lowerbounds) .+ lowerbounds for x in x1D, y in x1D][:]
  fxys = [f(xy) for xy in xys][:]
  return solve(f, xys, fxys; nroots=nroots, npoles=npoles, nsamples1D=nsamples1D, nretries=nretries, atol=atol, timelimitpertry=timelimitpertry)
end


function solve(f::Function, xys, fxys; nroots::Int, npoles::Int, nsamples1D::Int=14, nretries=0, atol=1e-3,
    timelimitpertry=5) where {F<:Function}
  dims = 2nroots + 2npoles
  rootsandpoles(p) = (p[1:2nroots], p[2nroots + 1:end])
  objective(p) = sum(onlyfinite(abs(fxyi - rational(xyi, rootsandpoles(p)...))) for (xyi, fxyi) in zip(xys, fxys))
  ntries = 0
  while true
    ntries += 1
    fit = BlackBoxOptim.bboptimize(objective; SearchRange=(0, 1), NumDimensions=dims,
                                   Method=:adaptive_de_rand_1_bin_radiuslimited,
                                   MaxTime=timelimitpertry, TraceMode=:silent, MaxFuncEvals=1_000_000)
    bf = best_fitness(fit)
    solution = reshape(best_candidate(fit), 2, dims ÷ 2)
    fsols = [f(s) for s in eachcol(solution)]
    p = sortperm(abs.(fsols))
    solution = solution[:, p]
    roots, poles = rootsandpoles(solution)
    if bf < atol || ntries >= nretries
      reason = bf < atol ? :converged : :maxretries 
      return (roots=roots, poles=poles, reason=reason, values=fsols, error=bf, ntries=ntries)
    else
      xys = vcat(xys, roots)
      fxys = vcat(fxys, fsols)
    end
  end
end

greet() = print("Hello World!")

end # module RationalRootsAndPoles
