# #
# # Black box circuit identification
# # (High-pass filter:
# #  Resistor + capacitor in series from input to output,
# #  inductor from output to ground)
# #

using Gadfly, LatexPrint


# Sine wave resonance

data = readcsv("lab5exp1.txt");

f = data[:,1];		# frequency, kHz
w = 2*pi.*f./1000;	# angular frequency, rad/s
#V_in = data[:,2];	# input voltage, V (should be 5)
V_in = 5*ones(size(f,1),1);
V_out = data[:,3];	# output voltage, V

f_L = f[23,1];
f_0 = f[32,1];
f_H = f[43,1];
fx = log10(f./f_H);
vee = V_out./V_in;
V = 20*log10(vee);

# Bode plot
p = plot(x=f, y=V, xintercept=[f_L,f_H],
	 Geom.point, Geom.vline, Scale.x_log10,
	 Guide.xlabel("frequency, kHz"),
	 Guide.ylabel("H_db"),
	 Guide.title("Bode plot of Black Box Circuit"),
	 Theme(minor_label_font_size=14pt,major_label_font_size=14pt));
draw(PNG("bode.png", 30cm, 15cm), p);

# print table of sine wave frequency-voltage data
tabular(data);


# Square wave dampened ringing

f_sq = 300;		# Hz
#V_in_sq = 5.20;	# V
V_in_sq = 5;
#V_e_sq = 5.20/e;	# V
V_e_sq = 5/e;
N = 1.70;		# number oscillations


# Determine R, L, C values

# Circuit has transfer function:
#     H(s) = (sL) / (sL * R + 1/(sC))
# To get an overall impedance equation for the circuit, we can write:
#     |V_out| = |V_in| * Z
# where
#     Z = sqrt( ((w*L)^2) / ((w*L)^2 + R^2 + (1/(w*C))^2) )
# so that
#     Z = |V_out / V_in|

# at resonance, we can ignore capacitance & inductance:
R = 1 / vee[32,1];

# at low frequencies, we can ignore inductance:
Zi = vee[1,1];
#C = 1 / (w[1,1] * sqrt((1/Zi)^2 - R^2));

# at high frequencies, we can ignore capacitance
# but this doesn't work out very well, so we can try something else:
#L = 1 / ((w[32,1])^2 * C);

w0 = w[32,1];
L = Q * w0 / R;
C = 1 / (L * w0^2);

