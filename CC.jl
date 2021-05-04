### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ acc935c2-d23d-4e64-bd84-fb96cd67d1ed
begin
	using Markdown
	using PlutoUI
	using CairoMakie
	using StatsBase
	using BenchmarkTools
end

# ╔═╡ 14e2a9a4-f625-433d-9552-a1b9be720bcd
xs = rand(0:1e-4:4_000.0, 500); ys = rand(0:1e-4:4_000.0, 600);

# ╔═╡ 0072a0cf-9396-4acf-9a02-a7a2c9cee2fe
thresh = 50

# ╔═╡ 6daf52cb-21db-4708-aeb0-830281298678
function lag_push(xs, ys)
	lags = []
	for x in xs, y in ys
		abs(x - y) ≤ thresh && push!(lags, x - y)
	end
	
	return lags
end

# ╔═╡ 57a438b6-56e6-4260-b423-63889f0e757b
myhist = fit(Histogram, lag_push(xs, ys),  -50:5:50, closed=:left)

# ╔═╡ 66a38b00-7d1f-4f18-8928-bf23d3710df1
myhist.weights

# ╔═╡ 96c8f81c-610e-454f-96f6-3942bd8e7253
plot(myhist)

# ╔═╡ 1b02b691-c141-47c8-922b-e79fc5c30bd6
# const length_x = length(xs)
# const length_y = length(ys)

# ╔═╡ 308fe3b9-398b-4c6a-b13d-e1812972393c
# function lag_prealloc()
# 	lags = zeros(Float64, length_x*length_y)
# 	i = 1
# 	for x in xs, y in ys
# 		abs(x - y) ≤ thresh && (lags[i] = abs(x - y))
# 		i += 1
# 	end
	
# 	return lags
# end

# ╔═╡ 2d833919-8538-4e88-a6a5-e1a1eca7ea6b
# filter!(!=(0.0, ), lag_prealloc())

# ╔═╡ Cell order:
# ╠═acc935c2-d23d-4e64-bd84-fb96cd67d1ed
# ╠═14e2a9a4-f625-433d-9552-a1b9be720bcd
# ╠═0072a0cf-9396-4acf-9a02-a7a2c9cee2fe
# ╠═6daf52cb-21db-4708-aeb0-830281298678
# ╠═57a438b6-56e6-4260-b423-63889f0e757b
# ╠═66a38b00-7d1f-4f18-8928-bf23d3710df1
# ╠═96c8f81c-610e-454f-96f6-3942bd8e7253
# ╠═1b02b691-c141-47c8-922b-e79fc5c30bd6
# ╠═308fe3b9-398b-4c6a-b13d-e1812972393c
# ╠═2d833919-8538-4e88-a6a5-e1a1eca7ea6b
