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

# ╔═╡ 3fc988b1-803f-4cb4-a8e3-0277c0785e83
md"""
In this notebook, we will take a look at quantifying the correlation between the signals fired from two separate neurons.

**The data:** A vector of time measurements from each neuron recording a firing event (in seconds):

``x\_i = (1, 2)``  \
``x\_j = (3)``

**The desired output:** A list of "circular lags" between the firing times between the two neurons. For the example above, we would expect the output to be:

```julia
lags = [-2, -1, 1, 0, -3, 0, -3, -2]
```

To get the above output, we discretize our time measuremets into 1 second intervals. A signal fire is then a binary event that can occur in a given interval of time. Adopting the latest firing event between the two neurons as the upper boundary on our given discretized grid, we can re-interpret the above firing measurements as the following binary firing matrix:

```math
X_0 = \begin{bmatrix}
	0 & 0 \\
    1 & 0 \\
	1 & 0 \\
	0 & 1
\end{bmatrix} \quad,
```

where row 1 corresponds to neuron ``i = 1``, row 2 to neuron ``j = 2``, and each column represents another second of time passing, starting at ``t_0 = 0\text{ s}``. A binary firing event then corresponds to ``1``, and to ``0`` otherwise.

Given the two measurements below, we can build ``X`` in the following way:
"""

# ╔═╡ 4a98643c-a0fa-40e8-89e8-19fdee5429e0
x₁_measurements = [1, 2]

# ╔═╡ f4241a88-e1dc-4a38-b05c-d59f2a1bca5b
x₂_measurements = [3]

# ╔═╡ ea6d6f2d-0b2f-4393-97cc-dc7ec1aa21f7
@doc raw"""
	X_repr(xᵢ_measurements::Vector{Real}, xⱼ_measurements::Vector{Real})

Returns the binary firing matrix ``X``, given two time measurements ``x\_i`` and ``x_j``
"""
function X_repr(xᵢ_measurements, xⱼ_measurements)
	n_time =  maximum((xᵢ_measurements, xⱼ_measurements))[1] + 1
	X = zeros(Int, n_time, 2)
	X[xᵢ_measurements .+ 1, 1] .= 1
	X[xⱼ_measurements .+ 1, 2] .= 1
	return X
end

# ╔═╡ 69471f45-b38e-4956-9160-3648b33c140b
X₀ = X_repr(x₁_measurements, x₂_measurements) # X₀

# ╔═╡ 753a7792-a541-4993-b9a5-61e34830b480
md"""
We then use `circshift` to periodically shift the first row and measure the resulting pair-wise lags between the two binary signals. For example, a circular shift of 1 second would create the new binary firing matrix:

```math
X₁ = \begin{bmatrix}
	0 & 0 & 1 & 1 \\
	0 & 0 & 0 & 1
\end{bmatrix}
```

With corresponing lags:

```math
[2, 3] \ominus [3] = [-1, 0] \quad,
```

where:
```math
\mathbf{u} \ominus \mathbf{v} \equiv [u_i - v_j],
\text { for } u_i \in \mathbf{u},\ v_j \in \mathbf{v}
\iff |u_i - v_j| \le \epsilon
```
is the vector of all pair-wise differences in signal firing times between the two neurons within a threshold ``\epsilon``. We show this computation for the current example:
"""

# ╔═╡ 1e1d8f53-3e05-4646-9c34-f18e0585d787
begin
	@doc raw"""
		compute_lags(u::Real, v::Real; ϵ::Real)
	
	Returns ``\mathbf{u} \ominus \mathbf{v}`` for the given threshold
	``\epsilon``, where:
	
	```math
	\mathbf{u} \ominus \mathbf{v} \equiv [u_i - v_j],
	\text { for } u_i \in \mathbf{u},\ v_j \in \mathbf{v}
	\iff |u_i - v_j| \le \epsilon
	```
	"""
	function compute_lags(u, v; ϵ=10.0)
		lags = Float64[]
		for uᵢ ∈ u, vⱼ ∈ v
			lag!(lags, uᵢ, vⱼ; ϵ=ϵ)
		end
		return lags
	end
	
	# Compute within-threshold lag
	lag!(lags, uᵢ, vⱼ; ϵ=10) = abs(uᵢ - vⱼ) ≤ ϵ && push!(lags, uᵢ - vⱼ)
end

# ╔═╡ b60ebaa3-789b-40f4-af09-524627cbbdb5
compute_lags([2, 3], [3]) # u ⊖ v

# ╔═╡ b637324c-13fa-4984-9f74-0546b10871e1
md"""
To compute the circularly shifted lags in general, we need a way to i) perform the appropriate shift of the binary firing matrix, and ii) convert to its decimal representation, ``\mathbf{u}`` and ``\mathbf{v}``. We can these tasks in the following way:
"""

# ╔═╡ 780876a3-f87c-49bd-9016-8dcebc3aea8e
# Circularly shift column 1 relative to column 2 by amount i
function shift_X(X₀, i)
	X_shift = copy(X₀)
	X_shift[:, 1] .= circshift(X₀[:, 1], i)
	return X_shift
end

# ╔═╡ eaf55391-9095-40c5-9ae2-a7f733ecb93a
# Convert binary firing signal to its decimal time representation
time_repr(u) = findall(==(1), u) .- 1

# ╔═╡ 9299b995-b349-493b-af37-ccbe7d36a464
md"""
So from our above example, circularly shifting from ``X_0 \to X_1``, and then converting to its decimal representation would look like the following:
"""

# ╔═╡ 7fa54f11-895e-4094-8051-ee339e74efed
X₁ = shift_X(X₀, 1)

# ╔═╡ 82af0727-c203-4b14-9ff4-17ef7d5bc213
X₁

# ╔═╡ dab74941-b125-41ab-8033-682804517a9f
time_repr(shift_X(X₀, 3)[:, 1])

# ╔═╡ 7cc8d1ff-41c1-466e-ae3d-b46fb017ce39
time_repr(X₀[:, 2])

# ╔═╡ 6daf52cb-21db-4708-aeb0-830281298678
function lag_push(xs, ys, thresh=10)
	lags = []
	for x in xs, y in ys
		lag = x - y
		abs(lag) ≤ thresh && push!(lags, lag)
	end
	
	return lags
end

# ╔═╡ c849709e-bbcb-43fc-b3d6-22618eb12af0
function siggies(X)
	x_i_binary_signal = let
		x = zeros(Int, N_TIME)
		x[x_i_measurements .+ 1] .= 1
		x
	end
		
	x_j_binary_signal = let
		x = zeros(Int, N_TIME)
		x[x_j_measurements .+ 1] .= 1
		x
	end
	
	x_i_signal = [
		findall(==(1), circshift(x_i_binary_signal, i)) .- 1
		for i in rand(1:N_TIME-1, 500)
	]
		
	x_j_binary_signal = let
		x = zeros(Int, N_TIME)
		x[x_j_measurements .+ 1] .= 1
		x
	end
	
	x_j_signal = findall(==(1), x_j_binary_signal) .- 1
	
	return lag_push.(x_i_signal, x_j_signal)
end

# ╔═╡ bbfdb68a-8bc2-49ce-8773-a985ff70767c
shifties = siggies(x₁_measurements, x₂_measurements, N_TIME)

# ╔═╡ 3284e1f1-4983-4fb2-85a3-0111a2c61124
shifties

# ╔═╡ 99669222-b406-489e-9c17-d670b860c0a0
yee = Base.Iterators.flatten(shifties) |> collect

# ╔═╡ 2945dd65-6bfa-4a5d-bdb1-b6275594e816
histies = fit(Histogram, yee,  -50:0.5:50, closed=:left)

# ╔═╡ ee09adbb-f5aa-4f98-9cf6-73e29591bd48
plot(histies)

# ╔═╡ fee66f35-982a-4a83-8e73-3f734f0cbc51
base_unpack() = Base.Iterators.flatten(shifties) |> collect

# ╔═╡ 14e2a9a4-f625-433d-9552-a1b9be720bcd
xs = rand(0:1e-4:4_000.0, 500); ys = rand(0:1e-4:4_000.0, 600);

# ╔═╡ 0072a0cf-9396-4acf-9a02-a7a2c9cee2fe
thresh = 50

# ╔═╡ ef064019-e8c8-4d93-87f1-c9919a2f2870
sort(rand(0:1e-4:4_000.0, (64, 100)), dims=2)

# ╔═╡ de15c8f0-3407-4b0a-99e1-0d36c46e6e71
rand(0:1e-4:4_000, 100)

# ╔═╡ c3894333-60b8-4625-95db-d40e9b9655f0
sample(1:10, 3)

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
# ╠═3fc988b1-803f-4cb4-a8e3-0277c0785e83
# ╠═4a98643c-a0fa-40e8-89e8-19fdee5429e0
# ╠═f4241a88-e1dc-4a38-b05c-d59f2a1bca5b
# ╠═69471f45-b38e-4956-9160-3648b33c140b
# ╠═ea6d6f2d-0b2f-4393-97cc-dc7ec1aa21f7
# ╟─753a7792-a541-4993-b9a5-61e34830b480
# ╠═b60ebaa3-789b-40f4-af09-524627cbbdb5
# ╠═1e1d8f53-3e05-4646-9c34-f18e0585d787
# ╟─b637324c-13fa-4984-9f74-0546b10871e1
# ╠═780876a3-f87c-49bd-9016-8dcebc3aea8e
# ╠═eaf55391-9095-40c5-9ae2-a7f733ecb93a
# ╟─9299b995-b349-493b-af37-ccbe7d36a464
# ╠═7fa54f11-895e-4094-8051-ee339e74efed
# ╠═82af0727-c203-4b14-9ff4-17ef7d5bc213
# ╠═dab74941-b125-41ab-8033-682804517a9f
# ╠═7cc8d1ff-41c1-466e-ae3d-b46fb017ce39
# ╠═6daf52cb-21db-4708-aeb0-830281298678
# ╠═c849709e-bbcb-43fc-b3d6-22618eb12af0
# ╠═bbfdb68a-8bc2-49ce-8773-a985ff70767c
# ╠═3284e1f1-4983-4fb2-85a3-0111a2c61124
# ╠═99669222-b406-489e-9c17-d670b860c0a0
# ╠═2945dd65-6bfa-4a5d-bdb1-b6275594e816
# ╠═ee09adbb-f5aa-4f98-9cf6-73e29591bd48
# ╠═fee66f35-982a-4a83-8e73-3f734f0cbc51
# ╠═14e2a9a4-f625-433d-9552-a1b9be720bcd
# ╠═0072a0cf-9396-4acf-9a02-a7a2c9cee2fe
# ╠═ef064019-e8c8-4d93-87f1-c9919a2f2870
# ╠═de15c8f0-3407-4b0a-99e1-0d36c46e6e71
# ╠═c3894333-60b8-4625-95db-d40e9b9655f0
# ╠═57a438b6-56e6-4260-b423-63889f0e757b
# ╠═66a38b00-7d1f-4f18-8928-bf23d3710df1
# ╠═96c8f81c-610e-454f-96f6-3942bd8e7253
# ╠═1b02b691-c141-47c8-922b-e79fc5c30bd6
# ╠═308fe3b9-398b-4c6a-b13d-e1812972393c
# ╠═2d833919-8538-4e88-a6a5-e1a1eca7ea6b
# ╠═acc935c2-d23d-4e64-bd84-fb96cd67d1ed
