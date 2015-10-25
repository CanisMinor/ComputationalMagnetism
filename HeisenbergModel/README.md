# Atomistic Modelling and the Heisenberg Model
In \textit{atomistic modelling} spins are modelled as semi-classical entities localised to lattice sites.  It can be shown that in the classical high-spin limit the Hubbard model yields the \textit{Heisenberg model} (see e.g.~\cite{Mielke,mielketasaki,tasaki2,tasaki3,hubheis,mattis}), which has the Hamiltonian
$$
 \hat{H}_{Heis} = -J_{ij} \sum_{i,j} \mathbf{S}_i \cdot \mathbf{S}_j
$$
The sum is normally constrained to nearest-neighbour interactions.  The time-dependent motion of semi-classical localised spins is described by the \textit{Landau-Lifshitz-Gilbert equation}
$$
 \frac{\partial \mathbf{S}}{\partial t} = - \frac{\gamma}{1+\alpha^2}\left(\mathbf{S}\times\mathbf{H}_{eff} + \alpha \mathbf{S}\times \left(\mathbf{S}\times \mathbf{H}_{eff}\right)\right)
$$
where $\mathbf{S}$ is the (normalised) spin moment, $\gamma$ is the gyromagnetic ratio, $\alpha$ is the \textit{Gilbert damping parameter} and $\mathbf{H}_{eff}$ is the effective magnetic field.  The Gilbert damping parameter represents the coupling of the spin system to a heat bath and determines the rate at which the moments relax to the direction of the effective field $\mathbf{H}_{eff}$.  The effective magnetic field is the first derivative of the  complete spin Hamiltonian
$$
 \mathbf{H}_{eff} = -\frac{\partial \hat{H}}{\partial\mathbf{S}}
$$
The effective field contains any thermal finite-temperature effects, derived from \textit{Langevin dynamics}.  The main principle of Langevin dynamics is to assume that thermal fluctuations can be represented as a Gaussian white-noise term.

## Initial Computational Implementations
For instructional purposes, an atomistic single-spin code was created for this project.  The first step was to code a single-spin version of the \textit{Heun scheme}, which is a technique for integrating the Landau-Lifshitz-Gilbert equation numerically.  The Heun scheme has the following structure (where $\mathbf{S}_1$ is the input spin vector and $\mathbf{S}_2$ is the output spin vector in a iterative process looking for convergence in $\mathbf{S}$).
1. The first (\textit{predictor}) step is
 $$
  \mathbf{S}_2 = \mathbf{S}_1 + \Delta \mathbf{S}_1 \Delta t
 $$
 where $\Delta t$ is the time step and
$$
 \Delta \mathbf{S}_1 = -\frac{\gamma}{1+\alpha^2}\left(\mathbf{S}_1\times\mathbf{H}_{eff} + \alpha \mathbf{S}_1\times \left(\mathbf{S}_1\times \mathbf{H}_{eff}\right)\right)
\label{ds1}
$$
2. The spin vector $\mathbf{S}_2$ is replaced with a normalised version of itself.
3. The second (\textit{corrector}) step is
$$
 \mathbf{S}_2 =  \mathbf{S}_1 + \frac{1}{2} \left(\Delta  \mathbf{S}_1 + \Delta \mathbf{S}_2\right) \Delta t
$$
with $\Delta \mathbf{S}_1$ as in Eq.~\ref{ds1} and
$$
  \Delta \mathbf{S}_2 = -\frac{\gamma}{1+\alpha^2}\left(\mathbf{S}_2\times\mathbf{H}_{eff} + \alpha \mathbf{S}_2\times \left(\mathbf{S}_2\times \mathbf{H}_{eff}\right)\right)
$$
4. The spin vector $\mathbf{S}_1$ is then replaced with a normalised version of $\mathbf{S}_2$ and the process (points 1-4) is then iterated until the input spin $\mathbf{S}_1$ is equal to the output spin $\mathbf{S}_2$.

### Testing the Dynamics of the Heun Scheme
The implemented Heun scheme was tested against the analytical solution (for derivation see appendix \ref{LLGanalytical})
\begin{eqnarray}
 S_x(t) & = & S_x(0)\cos(\omega t) \text{sech}\left(\frac{\alpha \gamma t}{1+\alpha^2}\right)\\
 S_y(t) & = & S_y(0) \sin(\omega t) \text{sech} \left(\frac{\alpha \gamma t}{1+\alpha^2}\right)\\
 S_z(t) & = & \text{tanh}\left(\frac{\alpha \gamma t }{1+\alpha^2}\right)
\end{eqnarray}
where $\omega=\gamma t/(1+\alpha^2)$.  This was tested for $\alpha=0.1$, $T=0$ K and $\gamma=1.76\times 10^{11}$ s$^{-1}$T$^{-1}$ starting out with the intial spin $\mathbf{S}=(1,0,0)$. The results are shown in fig.~\ref{precessions}.
\begin{figure}[ht]
\centering
\fbox{\parbox{17cm}{
\center
\subfigure[Computational and analytical time development of the spin x-component.]{
\includegraphics[scale=0.4, trim = 0mm 23mm 0mm 10mm, clip]{precessionx.pdf}
\label{precx}
}
\subfigure[Computational and analytical time development of the spin y-component.]{
\includegraphics[scale=0.4, trim = 0mm 23mm 0mm 10mm, clip]{precessiony.pdf}
\label{precy}
}
\subfigure[Computational and analytical time development of the spin z-component.]{
\includegraphics[scale=0.4, trim = 0mm 23mm 0mm 10mm, clip]{precessionz.pdf}
\label{precz}
}
\caption{\textit{Computational (from the Heun scheme) and analytical spin trajectory for a one-spin system described by the Landau-Lifshitz-Gilbert equation.}\label{precessions}}}}
\end{figure}




### Testing the Thermal Fields in the Heun Scheme
For modelling exchange in the Heun scheme implementation of the single-spin Heisenberg model, the effective field $\mathbf{H}_{eff}$ is a sum of an applied field $\mathbf{H}_{app}$ and a thermal field $\mathbf{H}_{\xi}$,
$$
 \mathbf{H}_{eff} = \mathbf{H}_{app} + \mathbf{H}_{\xi} = \mathbf{H}_{app} + \sqrt{\frac{2 k_B T \alpha}{\mu \Delta t}} \mathbf{\xi}
$$
where $\mathbf{\xi}$ is a vector containing random numbers between 0 and 1 from a Gaussian distribution.  These random numbers are obtained from a \textit{Mersenne-Twister} algorithm \cite{mersennetwister}.  The thermal field and the random number generator were checked to be appropriate by checking that the single-spin distribution follows a Boltzmann form
$$
 P(\theta) = P_0 \sin \theta \exp\left(\frac{E(\theta)}{k_BT}\right)
$$
where $\theta$ is the angle between the spin and the applied field, $P_0$ is a normalisation constant, $\sin \theta$ is the degeneracy factor (taking into account the number of states with the same angle $\theta$), $T$ is the temperature, $k_B$ is the Boltzmann constant and $E(\theta)$ is the magnetostatic (due to exchange and/or anisotropy) energy of a spin at angle $\theta$.\\*

The parameters used for the test were $\alpha=0.1$, $T=0.1$ K and $\gamma=1.76\times 10^{11}$ s$^{-1}$T$^{-1}$, starting out with the intial spin $\mathbf{S}=(0,0,1)$.  It was fitted with an equation of the form
$$
 P(\theta) = P_0 \sin \theta \exp\left(\frac{\left|\mathbf{S}\right|\left|\mathbf{H}\right|\cos(\theta)}{k_BT}\right)
$$
The result is shown in fig.~\ref{boltz}.  The fit gave $\left|\mathbf{S}\right|\left|\mathbf{H}\right|/(k_BT) = 6.72 \pm 0.01$ which agrees well with the actual value $\left|\mathbf{S}\right|\left|\mathbf{H}\right|/(k_BT) = 6.717...$
<!---
\begin{figure}[ht]
\centering
\fbox{\parbox{17cm}{
\center
\includegraphics[scale=0.4, trim = 0mm 23mm 0mm 10mm, clip]{boltz.pdf}
\caption{\textit{Angular distribution of spin orientations at thermal equilibrium.}\label{boltz}}}}
\end{figure}
--->
These test strongly support that the code works correctly for a single spin.  The next step in developing this code is to write it for a system of spins and to read in exchange tensor data and magnetic anisotropy energies.
