# #
# # p o w e r S p e c t r u m . j l
# #

function powerSpectrum{Tt<:Number}(ts::Vector{Tt}, sampleRate::Int)
    #= Return the absolute magnitude power spectrum in a vector mag,
	corresponding to x-axis frequency freq (units = Hz)
	Inputs:
		ts:		Column vector of input timestream data
		sampleRate:	sampling rate, in Hx, for ts data
    =#

    N = length(ts);
    halfN = floor(N/2);
   
    rawFFT = fft(ts);
    mag = abs(rawFFT[1:halfN,:]) * 2./N;

    freq = (0:1:halfN-1) .* (sampleRate/(halfN*2));
    return mag, freq;
end
