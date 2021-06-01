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

``x\_i = [1, 2]\text{ s}``  \
``x\_j = [3]\text{ s}``

**The desired output:** A list of "circular lags" between the firing times between the two neurons. For the example above, we would expect the output to be:

```julia
lags = [-2, -1, -1, 0, -3, 0, -3, -2]
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

where the first column corresponds to neuron ``i = 1``, the second column to neuron ``j = 2``, and each row represents another interval of time passing since ``t_0 = 0\text{ s}``. A binary firing event then corresponds to ``1``, and to ``0`` otherwise.

Given the two firing time measurements ``x_i`` and ``x_j``, we can build ``X_0`` in the following way:
"""

# ╔═╡ 4a98643c-a0fa-40e8-89e8-19fdee5429e0
xᵢ = [1, 2]

# ╔═╡ f4241a88-e1dc-4a38-b05c-d59f2a1bca5b
xⱼ = [3]

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
X₀ = X_repr(xᵢ, xⱼ) # X₀

# ╔═╡ 753a7792-a541-4993-b9a5-61e34830b480
md"""
We then use `circshift` to periodically shift the first neuron column and measure the resulting pair-wise lags between the two binary signals. For example, a circular shift of 1 time interval would create the new binary firing matrix:

```math
X_1 = \begin{bmatrix}
	0 & 0 \\
    0 & 0 \\
	1 & 0 \\
	1 & 1
\end{bmatrix} \quad,
```

with corresponing lags:

```math
[2, 3]\text{ s} \ominus [3]\text{ s} = [-1, 0]\text{ s} \quad,
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
		compute_lags(u::Vector{Real}, v::Vector{Real}; ϵ::Real)
	
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

# ╔═╡ b637324c-13fa-4984-9f74-0546b10871e1
md"""
To compute the circularly shifted lags in general, we need a way to i) perform the appropriate shift of the binary firing matrix, and ii) convert to its decimal representation, ``\mathbf{u}`` and ``\mathbf{v}``. We can accomplish these tasks in the following way:
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

# ╔═╡ 64ab60cb-9bb8-40ae-9fa2-cb2da33fa239
md"""
!!! note
	We offset by 1 to account for a shift of 0 time intervals.
"""

# ╔═╡ 9299b995-b349-493b-af37-ccbe7d36a464
md"""
So from our above example, circularly shifting from ``X_0 \to X_1``, and then converting to its decimal representation would look like the following:
"""

# ╔═╡ 7fa54f11-895e-4094-8051-ee339e74efed
X₁ = shift_X(X₀, 1)

# ╔═╡ dab74941-b125-41ab-8033-682804517a9f
u, v = time_repr(X₁[:, 1]), time_repr(X₀[:, 2])

# ╔═╡ b60ebaa3-789b-40f4-af09-524627cbbdb5
compute_lags(u, v) # u ⊖ v

# ╔═╡ 2eaca70c-5e92-4d72-9cd3-18d2ba982670
md"""
Putting this all together, we can now compute the lag times for all interval time shifts, as stated in the beginning of this notebook:
"""

# ╔═╡ 6c299b01-ae3d-4555-bf24-9432fa4f338b
function lag_vector(xᵢ_measurements, xⱼ_measurements)
	# xᵢ, xⱼ -> X₀
	X₀ = X_repr(xᵢ_measurements, xⱼ_measurements)
	# X₀ -> X₁, X₂, ⋯
	N_shifts = size(X₀)[1] - 1
	X_shifts = shift_X.(Ref(X₀[:, 1]), 0:N_shifts)
	# Xₙ -> uₙ, v
	us = time_repr.(X_shifts)
	v = time_repr(X₀[:, 2])
	# uₙ, v -> lags
	time_lags = compute_lags.(us, v)
	return Base.Iterators.flatten(time_lags) |> collect
end

# ╔═╡ c5e08864-b561-4943-b30d-ff867dece325
lags = lag_vector(xᵢ, xⱼ)

# ╔═╡ Cell order:
# ╟─3fc988b1-803f-4cb4-a8e3-0277c0785e83
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
# ╟─64ab60cb-9bb8-40ae-9fa2-cb2da33fa239
# ╟─9299b995-b349-493b-af37-ccbe7d36a464
# ╠═7fa54f11-895e-4094-8051-ee339e74efed
# ╠═dab74941-b125-41ab-8033-682804517a9f
# ╟─2eaca70c-5e92-4d72-9cd3-18d2ba982670
# ╠═c5e08864-b561-4943-b30d-ff867dece325
# ╠═6c299b01-ae3d-4555-bf24-9432fa4f338b
# ╠═acc935c2-d23d-4e64-bd84-fb96cd67d1ed
