# #
# # C h a p t e r 1 0 :
# # F o u r i e r  T r a n s f o r m  E x a m p l e s
# #

using WAV, Gadfly
#using AudioIO
include("powerSpectrum.jl")



# #
# # S10.1.1 : Sine Wave timestream
# # 

t = 0:1/44100:1;		# 1 s w/ 44100 samples
ts1 = sin(2.*pi*440.*t);	# 440 Hz sine wave (A note)

println("Plotting sine wave")
p1a = plot(x=t[1:5:end,:], y=ts1[1:5:end,:], Geom.point,
	  Guide.xlabel("Time, s"),
	  Guide.title("S10.1.1 - Sine Wave"),
	  Theme(default_point_size=.5pt,
		minor_label_font_size=14pt,major_label_font_size=14pt));

println("Saving sine wave")
draw(PNG("10.1.1.sinewave.png", 30cm, 15cm), p1a);

println("Playing wave")
#osc = SinOsc(440);
#play(osc)
#stop(osc)

high = sin(2.*pi*800.*t);
low = sin(2.*pi*80.*t);

#osc_high = SinOsc(800);
#osc_low = SinOsc(80);
println("Playing high sound")
#play(osc_high)
#stop(osc_high)
println("Playing low sound")
#play(osc_low)
#stop(osc_low)

# display the Fourier content of this timestream
mag1, freq1 = powerSpectrum(ts1, 441000);

println("Plotting sine wave Fourier Transform")
p1b = plot(x=freq1, y=mag1, Geom.point, Scale.x_log10,
	 Guide.xlabel("frequency, Hz"),
	 Guide.title("S10.1.1 - Fourier Transform"),
	 Theme(default_point_size=2pt,
		minor_label_font_size=14pt,major_label_font_size=14pt));

println("Saving sine wave Fourier Transform")
draw(PNG("10.1.1.fourier.png", 30cm, 15cm), p1b);



# #
# # S10.2.1 : Audio sample timestream
# #

println("Reading .WAV file")
ts2, fs = wavread("ventures_vol4_cut.wav");
@printf "Sample rate: %f\n" fs

println("Plotting timestream")
tf = 0:1/fs:(size(ts2,1)-1)/fs;
p2a = plot(
	layer(x=tf, y=ts2[:,1], Geom.point, Theme(default_color=color("magenta"))),
	layer(x=tf, y=ts2[:,2], Geom.point, Theme(default_color=color("cyan"))),
#x=0:1/fs:(length(ts2)-1)/fs, y=ts2,
	   Guide.title("S10.2.1 - Timestream"),
	   Theme(major_label_font_size=14pt));

println("Saving timestream")
draw(PNG("10.2.1.timestream.png", 30cm, 15cm), p2a);

# display the Fourier content of this timestream
mag2a, freq2a = powerSpectrum(ts2[:,1], 44100);
mag2b, freq2b = powerSpectrum(ts2[:,2], 44100);

println("Plotting WAV sample Fourier transform")
p2b = plot(
	 layer(x=freq2a, y=mag2a, Geom.point, Theme(default_color=color("magenta"))),
	 layer(x=freq2b, y=mag2b, Geom.point, Theme(default_color=color("cyan"))),
	 Scale.x_log10,
	 Guide.xlabel("Frequency, Hz"),
	 Guide.title("S10.2.1 - Fourier Transform"),
	 Theme(default_point_size=1pt,
		minor_label_font_size=14pt,major_label_font_size=14pt));

println("Saving WAV sample Fourier tranform")
draw(PNG("10.2.1.fourier.png", 30cm, 15cm), p2b);



# #
# # S10.2.3 : Signal buried in noise
# #

ts3 = .1*sin(2.*pi*440.*t);	# 440 Hz sine wave (A note)

println("Plotting sine wave")
p3a = plot(x=t[1:5:end,:], y=ts3[1:5:end,:], Geom.point,
	  Guide.xlabel("Time, s"),
	  Guide.title("S10.2.3 - Sine Wave"),
	  Theme(default_point_size=.5pt,
		minor_label_font_size=14pt,major_label_font_size=14pt));

println("Saving sine wave")
draw(PNG("10.2.3.sinewave.png", 30cm, 15cm), p3a);

println("Playing wave")
#osc3 = SinOsc(440);
#play(osc3)
#stop(osc3)

ts3_noise = randn(44101);

println("Plotting noisy sine wave")
p3a = plot(x=t[1:5:end,:], y=ts3[1:5:end,:], Geom.point,
	  Guide.xlabel("Time, s"),
	  Guide.title("S10.2.3 - Noisy Sine Wave"),
	  Theme(default_point_size=.5pt,
		minor_label_font_size=14pt,major_label_font_size=14pt));

println("Saving noisy sine wave")
draw(PNG("10.2.3.noisy_sinewave.png", 30cm, 15cm), p3a);

# TODO: play noisy sine wave

# hide A-note in the white noise
ts3 = ts3 + ts3_noise;

println("Plotting noisy sine wave with hidden A-note")
p3a = plot(x=t[1:5:end,:], y=ts3[1:5:end,:], Geom.point,
	  Guide.xlabel("Time, s"),
	  Guide.title("S10.2.3 - Noisy Sine Wave with hidden A-note"),
	  Theme(default_point_size=.5pt,
		minor_label_font_size=14pt,major_label_font_size=14pt));

println("Saving noisy A sine wave")
draw(PNG("10.2.3.noisy_A_sinewave.png", 30cm, 15cm), p3a);

# TODO: play noisy sine wave w/ hidden A-note

# Still, in fourier space, the signal is loud and clear
mag3, freq3 = powerSpectrum(ts3, 44100);

println("Plotting noisy A sine wave Fourier Transform")
p3a = plot(x=freq3, y=mag3, Geom.point,
	  Guide.xlabel("Frequency, Hz"),
	  Guide.title("S10.2.3 - Noisy A Sine Wave Fourier Transform"),
	  Theme(default_point_size=1pt,
		minor_label_font_size=14pt,major_label_font_size=14pt));

println("Saving noisy A sine wave Fourier Transform")
draw(PNG("10.2.3.noisy_A_sinewave_A_fourier.png", 30cm, 15cm), p3a);
