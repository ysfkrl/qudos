import numpy as np
import qutip as qt
from .pulse import ChirpedPulse

HBAR = 0.6582173  # meV*ps

def energies(delta_b=4., delta_0=0.):
    E_X = -delta_0 / 2
    E_Y = delta_0 / 2
    E_B = -delta_b
    return E_X, E_Y, E_B

def biexciton_system(
    collapse="nodecay", tau1=1, tau2=1, area1=1*np.pi, area2=0, det1=0, det2=0, alpha1=0, alpha2=0, prob_b=1/2,
    pol1_x=1, pol2_x=1, delay=1, delta_b=4, delta_0=0.0, gamma_e=1/100, gamma_b=1/100, epsilon=0.01, 
    dt_1=0.1, dt_2=0.1, options=qt.Options(atol=1e-7), mode="population"
):
    # Ensure delta_b and delta_0 are divided by HBAR
    delta_b = delta_b / HBAR
    delta_0 = delta_0 / HBAR

    # Define basis states
    g = qt.basis(4, 0)
    x = qt.basis(4, 1)
    y = qt.basis(4, 2)
    b = qt.basis(4, 3)

    gxbas = np.sqrt(1 - prob_b) * g + np.sqrt(prob_b) * b

    n_g, n_x, n_y, n_b = g * g.dag(), x * x.dag(), y * y.dag(), b * b.dag()
    p_gx, p_gy, p_xb, p_yb, p_gb = g * x.dag(), g * y.dag(), x * b.dag(), y * b.dag(), g * b.dag()

    c_ops = [] if collapse == "nodecay" else [np.sqrt(gamma_e) * p_gx, np.sqrt(gamma_e) * p_gy, np.sqrt(gamma_b) * p_xb, np.sqrt(gamma_b) * p_yb]

    E_X, E_Y, E_B = energies(delta_b=delta_b, delta_0=delta_0)
    H_sys = E_X * n_x + E_Y * n_y + E_B * n_b

    tau11, tau22 = np.sqrt(alpha1**2 / tau1**2 + tau1**2), np.sqrt(alpha2**2 / tau2**2 + tau2**2)
    t_start1, t_start2 = 4 * tau11, 4 * tau11 + delay
    pulse1, pulse2 = ChirpedPulse(tau1, det1, alpha1, t0=t_start1, e0=area1, polar_x=pol1_x), ChirpedPulse(tau2, det2, alpha2, t0=t_start2, e0=area2, polar_x=pol2_x)

    H = [H_sys, 
         [p_gx + p_xb, lambda t, args: np.conj(pulse1.get_total(t) + pulse2.get_total(t))],
         [p_gy + p_yb, lambda t, args: np.conj(pulse1.get_total(t) + pulse2.get_total(t))]]

    # Fix: Ensure `rate` is assigned before it's used
    t_off = t_start2 + t_start1
    rate = max(2 * gamma_b, gamma_e)  # Fix: Define rate before using it
    t_end = t_off - np.log(epsilon) / rate  # Now use the assigned rate

    t_axis = np.append(np.arange(0, t_off, dt_1), np.arange(t_off, t_end, dt_2))

    results = qt.mesolve(H, gxbas, t_axis, c_ops=c_ops, e_ops=[n_g, n_x, n_y, n_b, p_gx, p_xb, p_gb], options=options).expect
    return (*results, t_axis)  
