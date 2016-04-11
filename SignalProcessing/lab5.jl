# #
# # Black box circuit determination:
# # Capacitor from input to output,
# # with resistor & inductor in series from output to ground
# #

using Gadfly, LatexPrint


# Sine wave resonance

data = readcsv("lab5exp1.txt");

f = data[:,1];		# frequency, kHz
V_in = data[:,2];	# input voltage, V (should be 5)
V_out = data[:,3];	# output voltage, V

f_L = f[24,1];
f_0 = f[33,1];
f_H = f[44,1];
fx = log10(f./f_H);
V = 20*log10(V_out./V_in);

# Bode plot
p = plot(x=f, y=V, Geom.point, Scale.x_log10,
	 Guide.xlabel("frequency, kHz"),
	 Guide.ylabel("H_db"),
	 Guide.title("Bode plot of Black Box Circuit"),
	 Theme(minor_label_font_size=14pt,major_label_font_size=14pt));
draw(PNG("bode.png", 30cm, 15cm), p);

tabular(data);


# Square wave dampened ringing

f_sq = 300;		# Hz
V_in_sq = 5.20;		# V
V_e_sq = 5.20/e;	# V
osc = 1.70;		# number oscillations
