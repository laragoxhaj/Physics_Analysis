##
## i n t e g r a t i o n . j l
##

# Copy of a function I added in a fork of NumericalMath.jl

#= Same as above, except deals with arrays instead of vectors for batch integration:
   Assumes that each of the n separate columns represent n separate datasets and returns
   2xn vector, where the 1st row contains integral values, 2nd contains errors.
   Also can specify whether or not to find the area under the absolute value of the
   curve, or simply the curve itself.
=#
function atrapz{Tx<:Number, Ty<:Number}(x::Array{Tx}, y::Array{Ty}, absv::Bool)
	
	local m = size(x,1)
	local n = size(x,2)
	if (size(y,1) != m && size(y,2) != n)
		error("Arrays 'x', 'y' must be of the same dimensions")
	end

	r = zeros(1, n);
	for k in 2:m
		r[1,:] += (x[k,:] - x[k-1,:]) .* (y[k,:] + y[k-1,:]);
	end
	if (absv)
		r = abs(r);
	end

	# error -h^2/12 * (f'(b) - f'(a))	
	ha = x[2,:] - x[1,:];
	he = x[end,:] - x[end-1,:];
	ra = (y[2,:] - y[1,:]) ./ ha;
	re = (y[end,:] - y[end-1,:]) ./ he;
	err = abs(ha .* he ./ 12 .* (re - ra));

	return [r./2; err];
end
