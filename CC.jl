### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ b16b5d93-acaf-4440-bdab-2d523ff878b8
import(Pkg); Pkg.add("PlutoUI")

# ╔═╡ acc935c2-d23d-4e64-bd84-fb96cd67d1ed
begin
	### A Pluto.jl notebook ###
	# v0.14.2
	
	using Markdown
	using InteractiveUtils
end

# ╔═╡ 47c6a4e6-379f-4a78-8249-2431311c4397
using PlutoUI

# ╔═╡ 1a2168a6-ac78-11eb-3150-7ffde5908665


# ╔═╡ 14e2a9a4-f625-433d-9552-a1b9be720bcd
x = abs.(rand(1,10000)) * 2; y = abs.(rand(1,4000)) * 2;

# ╔═╡ 3af852e7-dfb7-47e0-adc8-22db644d3e2d
for i in range(

# ╔═╡ 698bfb46-54f1-4acb-aed0-3f3fbac6c997
let
diffsd = zeros(1,500000);
k = 2
	tmp=vec(y.-x[k]); 
	deleteat!(tmp, abs.(tmp).>.05);
	if k > 1
	diffsd[ max(findall(diffsd)) : (max(findall(diffsd))+length(tmp).-1)] = tmp;
		else
	diffsd[ 1 : (1+length(tmp).-1)] = tmp;
		end
    
    #ci = ci .+ length(tmp);

end

# ╔═╡ 37ddc497-612d-40ea-98aa-fc73c9cb518a
#with_terminal() do
begin
diffsd = zeros(1,500000);#
	ci = 1
for k = 1:10#length(x)
	if k > 1
	ci = max(findall(diffsd))
		else
			ci = 1
		end
    tmp=y.-x[k]; tmp[abs.(tmp).>.05] = [];
    diffsd[ ci : (ci+length(tmp).-1)] = tmp;
    #ci = ci .+ length(tmp);
end
end


# ╔═╡ d283df1c-da1c-4d17-b6dd-a474283755c2
diffsd2 = diffsd[1:(ci-1)];

# ╔═╡ Cell order:
# ╠═1a2168a6-ac78-11eb-3150-7ffde5908665
# ╠═acc935c2-d23d-4e64-bd84-fb96cd67d1ed
# ╠═b16b5d93-acaf-4440-bdab-2d523ff878b8
# ╠═47c6a4e6-379f-4a78-8249-2431311c4397
# ╠═14e2a9a4-f625-433d-9552-a1b9be720bcd
# ╠═d283df1c-da1c-4d17-b6dd-a474283755c2
# ╠═3af852e7-dfb7-47e0-adc8-22db644d3e2d
# ╠═698bfb46-54f1-4acb-aed0-3f3fbac6c997
# ╠═37ddc497-612d-40ea-98aa-fc73c9cb518a
