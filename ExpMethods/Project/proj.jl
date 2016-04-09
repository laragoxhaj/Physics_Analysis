include("integration.jl")
using Gadfly
using Cairo
using LatexPrint


data = readdlm("data_full.txt");	# at Voltages V = 3,4,5,6,7,8,9,10,11,12

d = 2*data[:,1:2:end];			# deflection at half amplitude, m
d_err = .0005*ones(size(d));
I = data[:,2:2:end];			# current, A
I_err = .00001*ones(size(I));
d_fin = cumsum(d);
d_err_fin = cumsum(d_err);

d_s = [.134,.134,.134,.147,.135,.133,.135,.136,.134,.130]';	# initial deflection, m
d_s_err = (.001/2)*ones(size(d_s));				


Rs = 30;		# resistance, ohms
Rs_err = 1;
n = 24;			# # turns on search coil
N = 2438;		# # turns on Rowland ring
dm = .160;		# diameter on rowland ring, m
dm_err = .001;
l = pi*dm;		# circumference of rowland ring, m
l_err = pi*dm_err;
A = 1.00/(100^2);	# area of search coil, m
A_err = .01/(100^2);
NPhi_s = 1.48/1000;	# S1H in W
NPhi_s_err = .01/1000;

K_b = (NPhi_s ./ (n * A * d_s));
K_b_err = sqrt((NPhi_s_err./(n*A*d_s)).^2 + (A_err*NPhi_s./(n*A^2*d_s)).^2 + (d_s_err*NPhi_s./(n*A*d_s.^2)).^2);
K_b_rep = repeat(K_b,outer=[size(d,1),1]);
K_b_err_rep = repeat(K_b_err,outer=[size(d,1),1]);
B = d_fin .* K_b_rep;
B_err = sqrt((d_err_fin.*K_b_rep).^2 + (d_fin.*K_b_err_rep).^2);
H = N*I./l;
H_err = sqrt(((N/l)*I_err).^2 + (l_err*(N/l^2)*I).^2);

Bmin = B .- B_err;
Bmax = B .+ B_err;

# Figures - Hysterisis Curves
for i in 3:12
	println("Saving plot $i")
	p = plot(x=H[:,i-2], y=B[:,i-2],
		 ymin=Bmin[:,i-2], ymax=Bmax[:,i-2],
		 Geom.point, Geom.errorbar, Geom.path,
		 Guide.xlabel("H, Amp-Turns m^{-1}"),
		 Guide.ylabel("B, Tesla"),
		 Guide.title("Hysteresis Curve of Toroidal Iron Ring at $i V"),
		 Theme(minor_label_font_size=14pt,major_label_font_size=14pt))
	draw(PNG("v$i.png", 30cm, 15cm), p)
end


#= integrate curves
   subtract area of segment 5 from area of segment 2;
   subtract area of segment 3 from area of segment 4;
   add the two results to find toal area of curve.
=#
EpCM = zeros(typeof(H[1,1]),2,10);	 # Total E for B-H curve, J/(m*cycle)
EpC = zeros(typeof(H[1,1]),2,10);
c = 0;
for i in 1:2
	a = 30*(i-1)+15;
	b = 30*i+15;
	if (i == 1)
		for k in 1:10
			if (k == 1)
				c = 38;
				a = 16;
			elseif (k == 2)
				c = 38;
			elseif (k <= 4)
				c = 36;
			elseif (k <= 7)
				c = 34;
			elseif (k == 8)
				c = 32;
			else
				c = 33;
			end			
		end
		H[c:b,i] = -1*H[c:b,i];
	else
		
		for k in 1:10
			if (k == 1)
				c = 71;
			elseif (k == 2)
				c = 68;
			elseif (k == 3)
				c = 67;
			elseif (k == 4 || k == 6 || k == 8)
				c = 65;
			elseif (k == 5)
				c = 66;
			elseif (k == 7 || k == 9)
				c = 64;
			else
				c = 63;
			end			
		end
		B[c:b,i] = -1*B[c:b,i];
	end
	EpCM[:,:] += atrapz(H[a:b,:], B[a:b,:], true);
end

V = A*l;				# Volume of ring
V_err = sqrt((A_err*l)^2 + (l_err*A)^2);
EpC[1,:] = V*EpCM[1,:];			# Energy per cycle of specific ring
EpC[2,:] = sqrt((V_err*EpCM[1,:]).^2 + (V*EpCM[2,:].^2));

EpCmin = EpC[1,:]-EpC[2,:];
EpCmax = EpC[1,:]+EpC[2,:];

# Figures - Energies per cycle, per cycle per meter
println("Saving Energy plots")
V = [3:12;]';
pc = plot(x=V, y=EpC[1,:],	  
	  ymin=EpCmin[:,:], ymax=EpCmax[:,:],
	  Geom.point, Geom.errorbar,
	  Guide.xlabel("Voltages, V"),
	  Guide.ylabel("Energy per Cycle, J/Cycle"),
	  Guide.title("Energies per Cycle, at Measured Voltages"),
	  Theme(minor_label_font_size=14pt,major_label_font_size=14pt))
draw(PNG("energy.png", 24cm, 12cm), pc);
