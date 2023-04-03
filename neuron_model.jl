"""
Traub model. Source: https://pubmed.ncbi.nlm.nih.gov/1663538/
"""

using DifferentialEquations
using Plots
import PyPlot as plt
using DataFrames
using FindPeaks1D
using LaTeXStrings
using DataStructures
using DataFrames
using CSV

# Plotting setup
pyplot()

rcParams = plt.PyDict(PyPlot.matplotlib."rcParams")
rcParams["figure.dpi"] = 150
rcParams["figure.figsize"] = [5, 3]
rcParams["lines.linewidth"] = 1.0
rcParams["font.size"] = 10.0

# %%
# --------------------------------------------------------------------------------
# Parameters
C_m  =   3.0;  # Membrane capacitance, in uF/cm^2

g_Na =  30.0;  # Sodium (Na) maximum conductance, in mS/cm^2
g_Ca = 4.0;  # Calcium (Ca) maximum conductance, in mS/cm^2
g_K_dr  =  25.0;  # Postassium (K) maximum conductances, in mS/cm^2
g_K_ahp = 0.8;
g_K_c = 10.0;
g_K_a = 5.0;
g_K_L  = 0.1;  # Conductance of Potassium leak current, in mS/cm^2
g_Na_L = 0.0175;  # Conductance of Sodium leak current, in mS/cm^2

# E_Na =  50.0;  # Sodium (Na) Nernst reversal potentials, in mV
# E_K = -80.0;
E_Ca  = 75.0; # Postassium (K) Nernst reversal potentials, in mV

Φ_Ca = 7769.0;  #26404.0 #7769.0;  # Ca scaling constant, (unit-less)
β_Ca = 0.075;  # Ca decay time constant, in ms^-1

A = 8e-6;  # Active cell surface area, in cm2, source: https://pubmed.ncbi.nlm.nih.gov/16344148, assuming sphere
F = 96485.;  # Faraday constant, in C/mol
R = 8.315;  # Universal gas constant, in J/mol*K
T = 297.5;  # Temperature (value=room temperate), in K

Vol_i = 2e-9;  # 8e-12  # Intracellular volume, in, cm3, source: https://pubmed.ncbi.nlm.nih.gov/16344148
Vol_e = 2e-9;  #  # Extracellular volume, in cm3

KM_Na = 10.;  # Michaelis constant of ATPase for Na in mM
KM_K = 3.5;  # Michaelis constant of ATPase for K, in mM

# I_pump_max = 5  # Vmax of ATPase. in uA/cm2

k    =  0.5;  # Temperature dependence constant, assuming Q10=2.3 and T=24.5C

# Initial values
v0 = -80.0;
m0 = 0.0;
h0 = 0.0;
s0 = 0.0;
r0 = 0.0;
n0 = 0.0;
a0 = 0.0;
b0 = 0.0;
q0 = 0.0;
c0 = 0.0;
Ca0 = 0.0;
Na_i0 = 8.5;  # mM
Na_e0 = 130.;  # mM
K_i0 = 130.;  # mM
K_e0 = 4.5;  # mM

# Simulation parameters
t0   =   0.0;  # First time point of simulation, in ms
tf   =  3000.0; # Last time point of simulation, in ms

# Stimulus parameters
stim_m = 3.125e1;  # Stimulus strength, in uA/cm^2
stim_s = 1000.0;  # Stimulus start time, in ms
stim_d = 1000.0;  # Stimulus duration, in ms

# I_pump_max values
I_pump_max_vals = vec([10. 20. 30. 40. 50. 60. 70. 80. 90. 100.])


# %%
# --------------------------------------------------------------------------------
# Model

function model(dx, x, p, t)

    # Unpack x
    v, m, h, s, r, n, a, b, q, c, Ca, Na_i, Na_e, K_i, K_e, I_pump_max = x
        
    # Compute input current based on time
    I_stim = stim_m * (1/(1 + (exp(-10.0*(t - stim_s)))) -
        1/(1 + (exp(-10.0*(t - (stim_s + stim_d))))))
    
    # Compute ion gradients
    E_Na = -R*T/F * log(Na_i/Na_e) * 1e3
    E_K = -R*T/F * log(K_i/K_e) * 1e3

    # Na-channel specific transmembrane current
    I_Na = g_Na * m^2 * h * (v - E_Na) +
            g_Na_L * (v - E_Na)

    # K-channel specific transmembrane current
    I_K =  g_K_dr * n * (v - E_K) +
            g_K_a * a * b * (v - E_K) +
            g_K_ahp * q * (v - E_K) +
            g_K_c * c * min(1.0, Ca/250.0) * (v - E_K) +
            g_K_L * (v - E_K)

    # Ca-channel specific transmembrane current
    I_Ca =
        g_Ca * s^2 * r * (v - E_Ca)
       
    # Pump facilitated current
    I_pump = I_pump_max * (1 + KM_Na/Na_i)^-3 * (1 + KM_K/K_e)^-2

    # Total transmembrane current
    I_p = I_Na + I_K + I_Ca + I_pump - I_stim

    # Channel kinetics
    a_m = k *  0.32*(-51.9 - v)/(exp((-51.9 - v)/4.0) - 1.0)
    b_m = k *  0.28*(v + 24.9)/(exp((v + 24.9)/5.0) - 1.0)
    a_h = k *  0.128*exp((-48.0 - v)/18.0)
    b_h = k *  4.0/(1.0 + exp((-25.0 - v))/5.0)
    a_s = k *  1.6/(1.0 + exp(-0.072*v))
    b_s = k *  0.02*(v + 13.9)/(exp((v + 13.9)/5.0) - 1.0)
    a_r = if (v <= -65.0) k * 0.005 else k * exp(-(v + 65.0)/20.0)/200.0 end
    b_r = if (v <= -65.0) k * 0.0 else k * 0.005 - a_r end
    a_n = k *  0.016*(-29.9 - v)/(exp((-29.9 - v)/5.0) - 1.0)
    b_n = k *  0.25 * exp((-45.0 - v)/40.0)
    a_a = k *  0.02*(-51.9 - v)/(exp((-51.9 - v)/10.0) - 1.0)
    b_a = k *  0.0175*(v + 24.9)/(exp(((v + 24.9)/10.0)) - 1.0)
    a_b = k *  0.0016*exp((-78.0 - v)/18.0)
    b_b = k *  0.05/(1.0 + exp((-54.9 - v)/5.0))
    a_q = k *  min(2e-5*Ca, 0.01)
    b_q = k *  0.001    
    a_c = if (v <= -15.0) k * exp((v + 55.0)/11.0 - (v + 58.5)/27.0)/18.975 else k * 2.0 * exp(-(v + 58.5)/27.0) end
    b_c = if (v <= -15.0) k * 2.0 * exp(-(v + 58.5)/27.0) - a_c else k * 0 end
    
    # dv/dt (v is relative to rest)
    dx[1] = -I_p / C_m

    # Channel states
    dx[2] = a_m * (1.0 - m) - b_m*m
    dx[3] = a_h * (1.0 - h) - b_h*h
    dx[4] = a_s * (1.0 - s) - b_s*s
    dx[5] = a_r * (1.0 - r) - b_r*r
    dx[6] = a_n * (1.0 - n) - b_n*n
    dx[7] = a_a * (1.0 - a) - b_a*a
    dx[8] = a_b * (1.0 - b) - b_b*b
    dx[9] = a_q * (1.0 - q) - b_q*q
    dx[10] = a_c * (1.0 - c) - b_c*c

    # Calcium dynamics
    dx[11] = -Φ_Ca*g_Ca*s^2*r*(v - E_Ca) - β_Ca*Ca

    # Na ion concentration dynamics
    dx[12] = -(I_Na + 3*I_pump) * A/(F*Vol_i) * 1e-3
    dx[13] =  (I_Na + 3*I_pump) * A/(F*Vol_e) * 1e-3

    # K ion concentration dynamics
    dx[14] = -(I_K - 2*I_pump) * A/(F*Vol_i) * 1e-3
    dx[15] =  (I_K - 2*I_pump) * A/(F*Vol_e) * 1e-3

    # return
    return dx
end


# %%
# --------------------------------------------------------------------------------

# Collection for spike trains and corresponding time
coll_t = Vector{Vector{Float64}}();
coll_v = Vector{Vector{Float64}}();

# Collection for resting state potential
coll_v_rest = Vector{Float64}();

# Collection for frequencies
coll_freq = Vector{Float64}();

# Collection of slopes
coll_slope = Vector{Float64}();

# Collection for first peak amplitude
coll_v_peak = Vector{Float64}();

# Index of peak to take peak potential from
ind = 1;

# Iterate through various v_max values
for i in range(1, length(I_pump_max_vals))
    
    # Extract current value
    I_pump_max = I_pump_max_vals[i]

    # Status
    println("Current I_pump_max value: ", I_pump_max)

    # Define problem
    prob = ODEProblem(model, [v0 m0 h0 s0 r0 n0 a0 b0 q0 c0 Ca0 Na_i0 Na_e0 K_i0 K_e0 I_pump_max_vals[i]], (t0, tf));

    # Solve ode
    sol = solve(prob) #, dt=0.001, adaptive=false);

    # Extract results
    t = sol.t;
    v, m, h, s, r, n, a, b, q, c, Ca, Na_i, Na_e, K_i, K_e = [reduce(vcat, sol.u)[:, i] for i in range(1, 15)]

    # Plot
    # -----

    plt.figure()
    plot(t, v, c="dodgerblue", label="") #, title=I_pump_max)
    # plot!(t1, v1, c="dodgerblue", label="V")
    xlabel!("time [ms]")
    ylabel!("membrane potential [mV]")
    # xlims!(800, 1700)
    display(ylims!(-100, 50))

    # Save spike train and corresponding time
    push!(coll_t, t)
    push!(coll_v, v)

    # V_rest
    # ------

    # Index of final time point before I_stim
    rest_ind = findall(x->x<stim_s, t)[end]

    # Resting state v
    v_rest = v[rest_ind]

    # Append v_rest to collection
    push!(coll_v_rest, v_rest);

    # Peaks
    # ------

    # Find membrane potential peaks
    x = findpeaks1d(v, height=-30);

    # Corresponding indexes
    peak_inds = x[1]

    println(length(peak_inds))

    # Frequency
    # -----

    freq = length(peak_inds) /
        (t[round(Int(peak_inds[end]))] - t[round(Int(peak_inds[1]))]) *
        1000

    # Append freq to collection
    push!(coll_freq, freq);
    
    # V_peak
    # -------

    # Take membrane potential at spike #[ind]
    v_peak = v[peak_inds]

    # Compute slope
    d = v_peak
    M = reduce(vcat, [[1 i] for i in range(1, length(v_peak))])
    b = M\d
    slope = b[2]

    # Append slope to collection
    push!(coll_slope, slope);

    # Append first peak's amplitude to collection
    push!(coll_v_peak, v_peak[ind]);

    # Plot peaks
    display(scatter!(t[peak_inds], v[peak_inds], label=""))


end

# %%
# ----------------------------------------------------------------------------
# Save results
# ----------------------------------------------------------------------------

# Create dataframe
df = DataFrame(OrderedDict(
    "t" => coll_t[1],
    "v" => coll_v[1])
    );

# Save dataframe
CSV.write("~/Desktop/spiketrain_10.csv", df);

# Create dataframe
df = DataFrame(OrderedDict(
    "t" => coll_t[end],
    "v" => coll_v[end])
    );

# Save dataframe
CSV.write("~/Desktop/spiketrain_100.csv", df);

# Create dataframe
df = DataFrame(OrderedDict(
    "I_pump_max" => I_pump_max_vals,
    "v_rest" => coll_v_rest,
    "freq" => coll_freq,
    "slope" => coll_slope,
    "v_peak" => coll_v_peak,)
    );

# Save dataframe
CSV.write("~/Desktop/traub_model_results.csv", df);
