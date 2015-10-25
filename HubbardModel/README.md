# The Hubbard Model

## The Hubbard Hamiltonian
The Hubbard model is a many-body model taking into account on-site electron-electron repulsion and hopping of electrons between lattice sites \cite{Hubbard,Hirsch,Tasaki}.  The model is typically used to investigate electron-electron interactions in narrow energy bands in tight-binding systems and for calculating magnetic properties that are not predicted by density-functional theory \cite{Hubbard,Hirsch,Tasaki}.  The Hubbard Hamiltonian has the form
$$
 \hat{H} = \sum_{i\neq j,\sigma} -t_{ij}\left(c_{i\sigma}^{\dagger}c_{j\sigma} + c_{i\sigma}c_{j\sigma}^{\dagger}\right) +\sum_i U_i n_{i\uparrow} n_{i\downarrow}
$$

The first term in the Hubbard Hamiltonian is the \textit{hopping term}.  It represents hopping of an electron from orbital $i$ to orbital $j$ as follows.  The \textit{creation operator} $c_{i\sigma}^{\dagger}$ adds one electron of spin $\sigma$ (which can be up or down) to orbital $i$ and the \textit{annihiliaton operator} $c_{j\sigma}$ removes one electron from orbital $j$ (the net effect being an electron has been transferred from orbital $i$ to orbital $j$).  The \textit{hopping integral} $t_{ij}$ is the quantum probability amplitude of transferring an electron from the $i^{th}$ to the $j^{th}$ orbital (the spin $\sigma$ of the electron remaining the same).  Normally the hopping sum is constrained to nearest-neighbour interactions only.  The hopping integral $t_{ij}$ is \cite{Hirsch,Tasaki}
$$
 t_{ij} = \langle \phi_i\left|\left(-\frac{\hbar^2 \nabla^2}{2m}\right)\right|\phi_j\rangle = \int \phi_i(x)\left(-\frac{\hbar^2 \nabla^2}{2m} \right) \phi_j(x) dx
$$
The second term in the Hubbard Hamiltonian is the \textit{Hubbard U term}.  The Hubbard $U$ is the potential (electrostatic) repulsion energy between two electrons of opposite spin occupying the same orbital.  Note that $n_{i\uparrow} n_{i\downarrow}$ can only be 1 or 0.  Normally $U$ is estimated from ab-initio calculations.  Inversely, we can also fit the Hubbard model to ab-initio calculations by varying $U$. Note that for $U=0$, the Hubbard model reduces to the tight-binding model.\\*

In the standard Hubbard model, it is assumed that each atom only have one electron orbit and that the corresponding orbital state is non-degenerate \cite{Tasaki}.  This is in general a good approximation to the real world where atoms have several electron orbits since the electrons in the other states are insignificant to the valence-level physics that we are interested in.\\*

The Hubbard Hamiltonian can be expanded in terms of heat baths of average electron densities $\langle n_{i,\uparrow}\rangle$ and $\langle n_{i,\downarrow} \rangle$ as follows.
$$
 \hat{H} \approx -\sum_{i\neq j,\sigma} t_{ij}\left(c_{i\sigma}^{\dagger}c_{j\sigma} + c_{i\sigma}c_{j\sigma}^{\dagger}\right) +\sum_i U_i \left[\langle n_{i\uparrow}\rangle  n_{i\downarrow} + \langle n_{i\downarrow}\rangle n_{i\uparrow} +  \langle n_{i\downarrow}  n_{i\uparrow} \rangle\right] + \sum_{i,\sigma} \varepsilon_{i\sigma}n_{i\sigma}
$$
This is because $n_{i,\sigma} = \langle n_{i,\sigma} \rangle + \delta (n_{i,\sigma})$ where $\delta (n_{i,\sigma})$ is the fluctuation in density from its time- and space average (these are the same according to the ergodic theorem).  The expansion is thus
\begin{eqnarray}
n_{i\uparrow}n_{i\downarrow} & = & \left(\langle n_{i\uparrow} \rangle + \delta (n_{i\uparrow})\right)\left(\langle n_{i\downarrow} \rangle + \delta (n_{i\downarrow})\right) \nonumber \\
& = & \langle n_{i\uparrow} \rangle\langle n_{i\downarrow} \rangle + \langle n_{i\downarrow} \rangle \delta (n_{i\uparrow}) + \langle n_{i\uparrow} \rangle\delta (n_{i\downarrow}) + \delta (n_{i\uparrow}) \delta (n_{i\downarrow}) \nonumber\\
& = & \langle n_{i\uparrow} \rangle\langle n_{i\downarrow} \rangle + \langle n_{i\downarrow} \rangle \left(n_{i\uparrow} - \langle n_{i\uparrow} \rangle\right) + \langle n_{i\uparrow} \rangle \left(n_{i\downarrow} - \langle n_{i\downarrow} \rangle\right) + \delta (n_{i\uparrow}) \delta (n_{i\downarrow}) \nonumber \\
& = &  \langle n_{i\downarrow} \rangle n_{i\uparrow} + \langle n_{i\uparrow} \rangle n_{i\downarrow} - \langle n_{i\uparrow} \rangle\langle n_{i\downarrow} \rangle + \delta (n_{i\uparrow}) \delta (n_{i\downarrow}) \nonumber \\
& \approx &  \langle n_{i\downarrow} \rangle n_{i\uparrow} + \langle n_{i\uparrow} \rangle n_{i\downarrow} - \langle n_{i\uparrow} \rangle\langle n_{i\downarrow} \rangle \\
\end{eqnarray}

Setting $\delta (n_{i\uparrow}) \delta (n_{i\downarrow}) \approx 0$, we can consider the Hubbard Hamiltonian as a sum of one Hamiltonian acting on up-spins only and one Hamiltonian acting on down-spins only.
\begin{eqnarray}
\label{sumupdown}
  \hat{H} \approx \hat{H}_{\uparrow} + \hat{H}_{\downarrow} & = & -\sum_{i\neq j} t_{ij}\left(c_{i\uparrow}^{\dagger}c_{j\uparrow} + c_{i\uparrow}c_{j\uparrow}^{\dagger}\right) +\sum_i U_i \langle n_{i\downarrow}\rangle  n_{i\uparrow} + \sum_{i} \varepsilon_{i\uparrow}n_{i\uparrow} \nonumber \\
 & & -\sum_{i\neq j} t_{ij}\left(c_{i\downarrow}^{\dagger}c_{j\downarrow} + c_{i\downarrow}c_{j\downarrow}^{\dagger}\right) +\sum_i U_i \langle n_{i\uparrow}\rangle  n_{i\downarrow} + \sum_{i} \varepsilon_{i\downarrow}n_{i\downarrow}
\end{eqnarray}

## Numerical Implementation
In this sub-project, the Hubbard model has been implemented using Python.  The code consists of a main program (doing the iterations towards self-consistency) and two subroutines (calculating the tight-binding and the Hubbard parts of the Hamiltonian matrix, respectively).  There is also a separate program which generates the neighbour list.
### Variables
The following physical variables are set in the program.
\begin{tabular}{| l | l |}
 \hline
 \textbf{Variable} & \textbf{Explanation}\\
 \hline
 $n_x$, $n_y$, $n_z$ & no of atoms in x-, y- and z-directions respectively\\
 \hline
 $n_{elec}$ & no of electrons per unit cell in the system\\
 \hline
 $N_{iter}$ & maximum no of self-consistency iterations\\
 \hline
 $N_s$ & total no of electrons\\
 \hline
 $t$ & hopping amplitude  \\
 \hline
 $U_{max}$ & final value of the Hubbard U\\
 \hline
 $a_x$, $a_y$, $a_z$ & lattice parameters along x, y and z\\
 \hline
 $\Delta U$ & increments of U\\
 \hline
 $\Delta k$ & increments in k\\
 \hline
\end{tabular}

### Structure of Main Program
The Hubbard model was implemented in a Python program using an algorithm as follows.
\begin{enumerate}
 \item All the physical variables were set to sensible values ($t\sim 1$, $U_{max} \sim 5$).  The number of electrons $n_{elec}$ is set to half-filling ($n_{elec}=N_s$, where $N_s$ is the number of sites).  The number of up-electrons $N_{\uparrow}$ and down-electrons $N_{\downarrow}$ are calculated from $n_{elec}$.
 \item An outer iteration loop was constructed to run for $N_{iter}$ iterations.
 \item For each iteration step, there are two spin-density integrations (one for spin-down and one for spin-up density) over k-space.  The first of these k-space integration loops is done for the spin-down case as follows.
       \begin{itemize}
        \item The tight-binding and Hubbard subroutines are called and the TB and Hubbard matrices are summed to give the spin-down Hamiltonian  $\hat{H}_{\downarrow}$.
        \item The spin-down Hamiltonian  $\hat{H}_{\downarrow}$ is passed to the lapack routine \verb zheev  to give the eigenvalues $E_1, E_2,...,E_{N_s}$ and the eigenstates $\phi_1, \phi_2,...,\phi_{N_s}$ of $\hat{H}_{\downarrow}$.  The eigenstates are listed in order of increasing energy.
        \item The spin-down density $\langle n_{i\downarrow}\rangle$ at each site $i$ is calculated by summing the squares of the absolute value of the coefficients up to the last state which has a spin-down electron in it, i.e.
              $$
                \langle n_{i\downarrow}\rangle = \left|\phi_{i,1}\right|^2 + \left|\phi_{i,2}\right|^2 + ...+ \left|\phi_{i,N_{\downarrow}}\right|^2
              $$
          \item At the end of the k-space integration, the sum of spin-down densities for each k-point is normalised by dividing it by the number of kpoints $n_k = n_{k_x}n_{k_y}$.
       \end{itemize}
 \item The second k-space integration loop is done for the spin-up case as follows.
        \begin{itemize}
         \item The spin-down density $\left(\langle n_{1\downarrow}\rangle,..., \langle n_{N_s\downarrow}\rangle\right)$ (some of these density elements will remain at zero unless all states have at least one spin-down electron) is fed into the Hubbard subroutine to give the spin-up Hubbard part of the Hamiltonian matrix.
         \item The tight-binding subroutine is called for the spin-up Hamiltonian.  The spin-up Hubbard and tight-binding parts are added together to give the spin-up Hamiltonian $\hat{H}_{\uparrow}$.
         \item The spin-up Hamiltonian  $\hat{H}_{\uparrow}$ is passed to the lapack routine \verb zheev  to give its eigenvalues $E_1, E_2,...,E_{N_s}$ and eigenstates $\phi_1, \phi_2,...,\phi_{N_s}$.  These are listed in order of increasing energy.
         \item The spin-up density $\langle n_{i\uparrow}\rangle$ at each site $i$ is calculated by summing the squares of the absolute value of the coefficients up to the last state which has a spin-up electron in it, i.e.
              $$
                \langle n_{i\uparrow}\rangle = \left|\phi_{i,1}\right|^2 + \left|\phi_{i,2}\right|^2 + ...+ \left|\phi_{i,N_{\uparrow}}\right|^2
              $$
          \item At the end of the k-space integration, the sum of spin-up densities for each k-point is normalised by dividing it by the number of kpoints $n_k = n_{k_x}n_{k_y}$.
         \end{itemize}
\item The spin-up density $\left(\langle n_{1\uparrow}\rangle,..., \langle n_{N_s\uparrow}\rangle\right)$ is then fed into the Hubbard subroutine for the spin-down Hamiltonian matrix.  This is repeated until the spin densities converge.  The code however stops incrementing $U$ once $U$ at each site reaches $U_{max}$.
\item For the first iteration, U is set to zero throughout to get the paramagnetic solution.  For the second iteration, the U on each site is set to a small (small with respect to $U_{max}$) random number.  From the third iteration onwards (until $U$ reaches $U_{max}$), $U$ is incremented in each iteration between the first and the second k-space integration.
\item The iteration loop runs until the Hubbard U has reached $U_{max}$ \textit{and} the spin-densities have converged.
\end{enumerate}

### Constructing the Neighbour List
The neighbour lists for each number of dimensions (code listed in appendix \ref{neighAPP}) was constructed as follows.
\begin{enumerate}
 \item The number of states ($n_x$, $n_y$, $n_z$) along each of the one, two or three axes is specified.
 \item The states are numbered as follows (see 3x3x3 unit cell example in fig.~\ref{neighnos}).
     \begin{itemize}
      \item Start out with the state at origin, label this state number one.
      \item Then move along x-axis until the $n_x^{th}$ atom and then increment along y-axis.  E.g.~if there are five states along the x-axis, state five will be the last state on the first row and state six will be the first state on the next row.
      \item Repeat until the y-axis is full.  Then increment along the z-axis.
     \end{itemize}
 \item On the left hand side of the unit cell, there will be $n_y n_z$ atoms.  These will be atoms no 1 to no $N_s-n_x+1$ incremented in steps of $n_x$.
 \item On the right hand side of the unit cell, there will be $n_yn_z$ atoms.  These will be atoms no $n_x$ to $N_s$ incremented in steps of $n_x$.
 \item On the back side of the unit cell, there will be $n_xn_y$ atoms.  These will be atoms no $n_x n_y(n_z-1)$ to no $N_s$ in steps of 1.
 \item On the front side of the unit cell, there will be $n_xn_y$ atoms.  These will be atoms no 1 to $n_xn_y$ in steps of 1.
 \item On the bottom side of the unit cell, there will be $n_xn_z$ atoms.  These are found as follows.  The first atom at the bottom side of the cell is atom with index $n_x(n_y-1)+1$.  Then follows atoms $n_x(n_y-1)+2,...,n_xn_y$.  Then there is a jump of $n_x(n_y-1)+1$ sites to reach the next row of the bottom face, which starts with atom number $n_x(2n_y-1)+1$.  Again the next $n_x$ atoms belong to the bottom face.  Then there is yet another jump of $n_x(n_y-1)+1$ sites and so on.
 \item On the top side of the unit cell, there will be $n_x n_z$ atoms.  These are found as follows.  The first atom at the top side is the atom with index $1$.  Then follows atoms $2,...,n_x$.  Then there is a jump of $n_x(n_y-1)+1$ sites to reach the next row of the top face which starts with atom number $n_xn_y+1$.  Again the next $n_x$ atoms belong to the top face.  Then yet another jump of $n_x(n_y-1)+1$ sites and so on.
\end{enumerate}

\begin{figure}[htp]
\centering
\fbox{\parbox{17cm}{
\center
\includegraphics[scale=0.6, trim = 0mm 0mm 0mm 0mm, clip]{neighnos.pdf}
\caption{\textit{A 3x3x3 unit cell as an example of the indexing of sites used our implementation of the Hubbard model.}\label{neighnos}}}}
\end{figure}

### Constructing the Tight-Binding Matrix
The \verb hopping  subroutine constructs the tight-binding part $\hat{H}_t$ of the total Hamiltonian matrix.
\begin{enumerate}
 \item Check how many states there are in each unit cell.  This is irrespective of the number of dimensions, i.e.~a three-dimensional 2x2x2 cubic unit cell will have eight states just as an octoatomic one-dimensional chain will have eight states.  For a system with $N_s$ states, the matrix $\mathbf{H}$ should be assigned dimensions $N_s$x$N_s$.
 \item Put the eigenenergy $E_0$ of each state along the diagonal.
 \item Put the hopping terms in as follows.  Then check the following for each pair of states $\langle \phi_i |$ and $|\phi_j\rangle$.
      \begin{itemize}
       \item Can an electron in $\langle \phi_i |$ hop to $|\phi_j\rangle$?  This is essentially checking which atoms that are neighbours using the neighbour list.
       \item Is this hop inside the unit cell? If so place $-t$ in the corresponding position in the Hamiltonian matrix.
       \item If the hop is outside the unit cell, along which axis does it happen?  If along $\mathbf{x}$ ($\mathbf{y}$, $\mathbf{z}$) put in phase shift term $-te^{-ik_x a}$ ($-te^{-ik_y b}$, $-te^{-ik_z c}$).  If the hop is backwards along the axis, delete the minus sign.
       \item If two states are connected in more than one way (e.g.~can hop both inside and outside the unit cell), add the corresponding hopping terms together in the same Hamiltonian matrix element.
      \end{itemize}
\end{enumerate}
Examples of tight-binding matrices for different systems can be found in section \ref{testing}.

### Constructing the Hubbard Matrix
The \verb hubU  subroutine constructs the Hubbard part $\hat{H}_U$ of the total Hamiltonian matrix.  This is done by inputting the dimensions of the matrix $N_s$x$N_s$, the array $U$ (holding the value of $U$ at each lattice site) and the average spin-up density $\langle n_{\uparrow,1}\rangle, \langle n_{\uparrow,2}\rangle, ..., \langle n_{\uparrow,N_s}\rangle$ (for $\hat{H}_{\downarrow}$) or the average spin-down density $\langle n_{\downarrow,1}\rangle, ...\langle n_{\downarrow,N_s}\rangle$ (for $\hat{H}_{\uparrow}$).  For the spin-up Hamiltonian, the Hubbard matrix will thus have a form as follows.
$$
 \hat{H}_{U\uparrow} = \left(\begin{array}{c c c c}
  U_1\langle n_{\downarrow,1}\rangle & 0 & ... & 0\\
  0 & U_2\langle n_{\downarrow,2}\rangle & ...& 0\\
  ... & ... & ... &...\\
  0 &...& ...& U_{N_s} \langle n_{\downarrow,N_s}\rangle\\
 \end{array}\right)
$$
Similarly, the Hubbard bit of the spin-down Hamiltonian is
$$
 \hat{H}_{U\downarrow} = \left(\begin{array}{c c c c}
  U_1\langle n_{\uparrow,1}\rangle & 0 & ... & 0\\
  0 & U_2\langle n_{\uparrow,2}\rangle & ...& 0\\
  ... & ... & ... &...\\
  0 &...& ...& U_{N_s} \langle n_{\uparrow,N_s}\rangle\\
 \end{array}\right)
$$

## Tests
### Tests against Tight-Binding Theory
For $U=0$ the Hamiltonian in eq.~\ref{HubHam} becomes identical to the Hamiltonian of the tight-binding model. In this section a number of analytical results for the tight-binding model are derived.  These analytical results are used to compare test output of the Hubbard code when setting $U=0$.

#### Monatomic Lattices
Consider an electron in a monatomic lattice.  Its wavefunction $| \mathbf{k}(\mathbf{r})\rangle$ can be written as
\begin{eqnarray}
 |\mathbf{k}(\mathbf{r})\rangle & = & \frac{1}{\sqrt{N}}\sum_{m} e^{i\mathbf{k}\cdot \mathbf{r_m}}|\mathbf{r}-\mathbf{r_m}\rangle\\
 \langle \mathbf{k}(\mathbf{r})| & = & \frac{1}{\sqrt{N}}\sum_{n} e^{-i\mathbf{k}\cdot \mathbf{r_n}}\langle\mathbf{r}-\mathbf{r_n}|\\
\end{eqnarray}
where $|\mathbf{r}-\mathbf{r_m}\rangle$ is the orbital on the $m^{th}$ atom in the lattice.  The Schr\"odinger equation can be written as
$$
 \hat{H} |\mathbf{k}(\mathbf{r})\rangle = E(\mathbf{k}) |\mathbf{k}(\mathbf{r})\rangle \Rightarrow E(\mathbf{k}) = \langle \mathbf{k}(\mathbf{r}) | \hat{H} |\mathbf{k}(\mathbf{r})\rangle
$$
Evaluating this gives
$$
 E(\mathbf{k}) = \langle \mathbf{k}(\mathbf{r}) | \hat{H} |\mathbf{k}(\mathbf{r})\rangle  = \frac{1}{N} \sum_{m} \sum_{n} e^{i\mathbf{k}\cdot\left(\mathbf{r}_m-\mathbf{r}_n\right)}\langle \mathbf{r}-\mathbf{r}_n | \hat{H} | \mathbf{r}-\mathbf{r}_m \rangle
\label{eneq}
$$
Denote the overlap integrals as follows
\begin{eqnarray}
 \langle \mathbf{r}-\mathbf{r}_m | \hat{H} | \mathbf{r}-\mathbf{r}_m \rangle & = & -E_0\\
 \langle \mathbf{r}-\mathbf{r}_n | \hat{H} | \mathbf{r}-\mathbf{r}_m \rangle & = & -t
\end{eqnarray}
if $m$ and $n$ are neighbouring sites.  The energy in eq.~\ref{eneq} can then be written as
$$
  E(\mathbf{k})  = -N\frac{1}{N}E_0 -N\frac{1}{N}t \sum_{n} e^{i\mathbf{k}\cdot\left(\mathbf{r}_m-\mathbf{r}_n\right)} = -E_0 -t \sum_{n} e^{i\mathbf{k}\cdot\left(\mathbf{r}_m-\mathbf{r}_n\right)}
$$
The expression can now only be further evaluated by looking at the possible separations $\left(\mathbf{r}_m-\mathbf{r}_n\right)$ between nearest neighbours in the lattice.  This depends on the lattice type.

\paragraph{One-Dimensional Chain}
In the case of the one-dimensional chain with periodic boundary conditions, $\left(\mathbf{r}_m-\mathbf{r}_n\right) = (\pm a,0,0)$ where $a$ is the lattice constant.  Thus the energy becomes
$$
  E(\mathbf{k})  = -E_0 -t \left(e^{-ik_x a} + e^{ik_x a}\right) = - E_0 - 2t\cos(k_x a)
 \label{oneDmonatomic}
$$

\paragraph{Two-Dimensional Sheet}
\ \\
For a two-dimensional sheet with periodic boundary conditions, nearest neighbours are separated by $\left(\mathbf{r}_m-\mathbf{r}_n\right) = (\pm a,0,0), (0, \pm b, 0)$.  Thus the energy becomes
$$
  E(\mathbf{k})  = -E_0 -t \left(e^{-ik_x a} + e^{ik_x a} + e^{-ik_y b} + e^{ik_y b}\right) = - E_0 - 2t\left(\cos(k_x a) + \cos(k_y b)\right)
 \label{twoDmonatomic}
$$
The two-dimensional code listed in appendix \ref{HubbardFORTRAN} was run for $n_x=n_y=1$, $a=b=1$, $t=1$ and the result is shown in fig.~\ref{2DMonat}.
\begin{figure}[htp]
\centering
\fbox{\parbox{17cm}{
\center
\subfigure[Band structure along $k_x$ ($k_y=0$).]{
\includegraphics[scale=0.33, trim = 30mm 23mm 30mm 10mm, clip]{band_kx.pdf}
\label{kx_monat2D}
}
\subfigure[Band structure along $k_y$ ($k_x=0$).]{
\includegraphics[scale=0.33, trim = 30mm 23mm 30mm 10mm, clip]{band_ky.pdf}
\label{ky_monat2D}
}
\caption{\textit{Test of the Hubbard model with $U=0$ for a two-dimensional monatomic lattice.}\label{2DMonat}}}}
\end{figure}

\begin{figure}[htp]
\centering
\fbox{\parbox{17cm}{
\center
\includegraphics[scale=0.6, trim = 10mm 20mm 10mm 20mm, clip]{MonatomicTB.pdf}
\caption{\textit{Test of the Hubbard model with $U=0$ for a two-dimensional monatomic lattice.}\label{MonatTB}}}}
\end{figure}

This does agree very well with equation \ref{twoDmonatomic} as well as (pictorially) with the equivalent result in \cite{Hirsch}.





\subsubsection{One-Dimensional Diatomic Chain}
\begin{figure}[htp]
 \begin{center}
\fbox{\parbox{17cm}{
\center
\includegraphics[scale=0.4, trim = 0mm 0mm 0mm 0mm, clip]{diatomic.pdf}
\caption{A diatomic chain in one dimension with periodic boundary conditions.}\label{diat}}}
\end{center}
\end{figure}
The tight-binding Hamiltonian for a one-dimensional diatomic chain is
$$
 \hat{H} = \sum_i \left(E_{A(i)}c_{A(i)}^{\dagger}c_{A(i)} + E_{B(i)}c_{B(i)}^{\dagger}c_{B(i)} - t c_{A(i)}^{\dagger}c_{B(i)} - t c_{B(i)}^{\dagger}c_{A(i)} - t c^{\dagger}_{B(i-1)}c_{A(i)} - tc^\dagger_{A(i+1)}c_{B(i)}\right)
$$
where all the rightward hops are written out explicitly and the possible leftward hops are included in the hermitian conjugate (h.c.).  For a single unit cell with index  $i$,
\begin{itemize}
 \item the one electron in $|k_{A(i)}\rangle$ can hop to $|k_{B(i)}\rangle$ or $|k_{B(i-1)}\rangle$
 \item the one electron in $|k_{B(i)}\rangle$ can hop to $|k_{A(i)}\rangle$ or $|k_{A(i+1)}\rangle$
\end{itemize}
We can write the energy as (assume $u_A$ and $u_B$ are real)
\begin{eqnarray}
 E(k) = \langle k | \hat{H} | k\rangle & = & u_A^* \langle k_A | k_A \rangle u_A + u_B^* \langle k_B | k_B \rangle u_B \\
& & + u_A^*\langle k_{A(i)}| k_{B(i)} \rangle u_B + u_A^*\langle k_{A(i)}| k_{B(i-1)} \rangle u_B + u_B^*\langle k_{B(i)}| k_{A(i)} \rangle u_A + u_B^*\langle k_{B(i)}| k_{A(i+1)} \rangle u_A\\
 & = & u_A^* E_A u_A + u_B^* E_B u_B - u_A^* t u_B + u_A^*e^{ika}\langle k_{A(i)}|  k_{B(i)} \rangle u_B - u_B^* t u_A + u_B^*e^{-ika}\langle k_{B(i)}|  k_{A(i)} \rangle u_A\\
& = & u_A^* E_A u_A + u_B^* E_B u_B - u_A^* t u_B - u_A^* t e^{ika} u_B - u_B^* t u_A - u_B^* t e^{-ika} u_A\\
& = & u_A^* E_A u_A + u_B^* E_B u_B - u_A^* t(1+e^{ika}) u_B - u_B^* t(1+e^{-ika}) u_A\\
\end{eqnarray}

This can be expressed as a matrix equation (setting $E=0$)
$$
 \left(u_A^* \ \ u_B^* \right) \left(\begin{array}{c c}
        0 & -t\left(1+e^{ika}\right)\\
        -t\left(1+e^{-ika}\right) & 0
       \end{array}\right)\left(\begin{array}{c}
                         u_A \\
			 u_B
                        \end{array}\right) = \left(u_A^* \ \ u_B^* \right) E(\mathbf{k})\left(\begin{array}{c}
                                                    u_A \\
						    u_B
                                                   \end{array}\right)
$$
In practice we use
$$
 \left(\begin{array}{c c}
        0 & -t\left(1+e^{ika}\right)\\
        -t\left(1+e^{-ika}\right) & 0
       \end{array}\right)\left(\begin{array}{c}
                         u_A \\
			 u_B
                        \end{array}\right) = E(\mathbf{k})\left(\begin{array}{c}
                                                    u_A \\
						    u_B
                                                   \end{array}\right)
$$
This gives the secular equation (whenever the determinant is zero, the eigenvectors $(u_A, u_B)$ are non-trivial)
$$
 \left|\begin{array}{c c}
        - E(k) & -t\left(1+e^{ika}\right)\\
	-t\left(1+e^{-ika}\right) & -E(k)
       \end{array}\right| = 0 \Rightarrow E(k)^2 = t^2\left(2+2\cos (ka) \right)\Rightarrow E(\mathbf{k}) = \pm t \sqrt{2\left(1+\cos (ka) \right)}
\label{diat_band}
$$

\subsection{Testing the Full Hubbard Model}
For $U\neq0$, $t\neq 0$, we have the full Hubbard model.  There are three interesting cases to consider.
\begin{itemize}
 \item the filled system, i.e.~the number of electrons $N_e=2N$ where $N$ is the number of available orbitals (since there can be a maximum of two electrons - of opposite spin - in each orbital)
 \item the half-filled system, i.e.~the number of electrons $N_e = N$
 \item the empty system, i.e.~the number of electrons $N_e=0$
\end{itemize}
It is immediately evident the the empty and the filled cases represent an insulator.  The half-filled system does however exhibit a number of interesting phenomena and thus the discussion will centre on the half-filled case from here on in.  In order to check that the Hubbard U is correctly implemented, the site populations of up- and down-spins on a two-dimensional 2x2 lattice with periodic boundary conditions (fig.~ were calculated for varying values of $U$ (fig.~\ref{hubUpop}).
\begin{figure}[ht]
\centering
\fbox{\parbox{17cm}{
\center
\includegraphics[scale=0.5, trim = 0mm 0mm 0mm 0mm, clip]{sites.pdf}
\label{hubUsites}
\caption{\textit{A two-dimensional 2x2 lattice with periodic boundary conditions.}\label{hubUpop}}}}
\end{figure}

\begin{figure}[ht]
\centering
\fbox{\parbox{17cm}{
\center
\subfigure[The relative population of up-spins at each site.]{
\includegraphics[scale=0.3, trim = 20mm 15mm 20mm 10mm, clip]{UpOcc.pdf}
\label{uppop}
}
\subfigure[The relative population of down-spins at each site.]{
\includegraphics[scale=0.3, trim = 20mm 15mm 20mm 10mm, clip]{DownOcc.pdf}
\label{downpop}
}
\caption{\textit{Relative populations of up- and down-spins for different Hubbard $U$.}\label{hubUpop2}}}}
\end{figure}
Thus the higher the Hubbard U, the more antiferromagnetic the system becomes.  This is expected as a high U means that electrons are less likely to occupy the same site.  At the same time, the antiferromagnetic arrangement means that hopping is possible, which lowers the energy with respect to the ferromagnetic state.
\clearpage





\begin{thebibliography}{99}


\bibitem{maghist} \'E.~du Tr\'emolet de Lacheisserie, D.~Gignoux and M.~Schlenker\\
\textit{Magnetism: Fundamentals}\\
\textbf{Berlin: Springer 2005}

\bibitem{gordon} G.~Hughes\\
\textit{Hard Drive!: As the Disc Turns}\\
\textbf{Charleston: BookSurge Publishing 2007}

\bibitem{magrechist} M.~Camras\\
\textit{Magnetic Recording Handbook}\\
\textbf{Berlin: Springer 1988}

\bibitem{magrec2} E.D.~Daniel, C.D.~Mee and M.H.~Clark\\
\textit{Magnetic Recording: The First 100 Years}\\
\textbf{New York: Wiley-IEEE 1999}

\bibitem{thiele} J.U.~Thiele, S.~Maat and S.S.~Fullerton\\
\textit{FeRh/FePt Exchange Spring Films for Thermally Assisted Magnetic Recording Media}\\
\textbf{Appl.~Phys.~Lett.~82 (2003) 2859}

\bibitem{thiele2} J.U.~Thiele, S.~Maat, J.L.~Robertson and S.S.~Fullerton\\
\textit{Magnetic and Structural Properties of FePt–FeRh Exchange Spring Films for Thermally Assisted Magnetic Recording Media}\\
\textbf{IEEE Trans.~Magn.~40 (2004) 2537}

\bibitem{HAMRrev} R.E.~Rottmayer et.al.\\
\textit{Heat-Assisted Magnetic Recording}\\
\textbf{IEEE Trans.~Mag.~42 (2006) 2417-2421}

\bibitem{iwata} S.~Iwata, S.~Yamashita and S.~Tsunashima\\
\textit{Perpendicular magnetic anisotropy and magneto-optical Kerr spectra of MBE-grown PtCo alloy films}\\
\textbf{IEEE Trans.~Magn.~33 (1997) 3670-3672}

\bibitem{iwata2} S.~Yamashita, S.~Iwata and S.~Tsunashima\\
\textit{Magnetic Anisotropy and Magneto-Optical Effect of MBE-Grown PtCo Alloy Films}\\
\textbf{J.~Magn.~Soc.~Jpn.~21 (1997) 433}


%general magnetism
\bibitem{skomski} R.~Skomski\\
\textit{Simple Models of Magnetism}\\
\textbf{Oxford: Oxford University Press 2008}

\bibitem{blundell} S.~Blundell\\
\textit{Magnetism in Condensed Matter}\\
\textbf{Oxford: Oxford University Press 2001}

\bibitem{mmmjiles} D.~Jiles\\
\textit{Magnetism and Magnetic Material 2E}\\
\textbf{Boca Raton: Taylor  \& Francis 1998}

%spin-spiral stuff
\bibitem{SQdynamics} Q.~Niu and L.~Kleinman\\
\textit{Spin-Wave Dynamics in Real Crystals}\\
\textbf{Phys.~Rev.~Lett.~80 (1998) 2205-2208}

\bibitem{dEFourier1} S.V.~Halilov et.al.\\
\textit{Adiabatic Spin Dynamics from Spin-Density-Functional Theory: Application to Fe, Co and Ni}\\
\textbf{Phys.~Rev.~B 58 (1998) 293-302}

\bibitem{dEFourier2} L.M.~Sandratskii and P.~Bruno\\
\textit{Exchange Interactions and Curie Temperature in (Ga,Mn)As}\\
\textbf{Phys.~Rev.~B 66 (2002) 134435}

\bibitem{noncollmag} N.~Mizuno et.al.\\
\textit{Non-Collinear Magnetism and Exchange Interaction in Spin-Spiral Structures of Thin-Film Fe(110)}\\
\textbf{J.~Phys.: Condens. Matter 19 (2007) 365222}

%Hartree-Fock etc
% \bibitem{Fuji} S.~Fujii, S.~Ishida, and S.~Asano.\\
% \textit{Electronic Structure and Lattice Transformation in Ni$_2$MnGa and Co$_2$NbSn}\\
% \textbf{J.~Phys.~Soc.~Jpn.~58 (1989) 3657}

\bibitem{ashcroft} N.W.~Ashcroft and N.D.~Mermin\\
\textit{Solid State Physics}\\
\textbf{Brooks Cole 1976}

\bibitem{BornOppenheimer} M.~Born and R.~Oppenheimer\\
\textit{Zur Quantentheorie der Molek\"ule}\\
\textbf{Annalen der Physik 84 (1927) 457–484}

\bibitem{Hartree1} D.R.~Hartree\\
\textit{The Wave Mechanics of an Atom tvith a Non-Voulomb Ventral Field. Part I. Theory and Method.}\\
\textbf{Proc.~Cam.~Phil.~Soc.~24 (1928) 89-110}

\bibitem{Hartree2} D.R.~Hartree\\
\textit{The Wave Mechanics of an Atom with a Non-Coulomb Central Field. Part II. Some Results and Discussion.}\\
\textbf{Proc.~Cam.~Phil.~Soc.~24 (1928) 111-132}

\bibitem{Hartree3} D.R.~Hartree\\
\textit{The Wave Mechanics of an Atom with a Non-Coulomb Central Field. Part III. Term Values and Intensities in Series in Optical Spectra.}\\
\textbf{Proc.~Cam.~Phil.~Soc.~24 (1928) 426}

\bibitem{PauliSymmetry}
\textit{N\"aherungsmethode zur L\"osung des Quantenmechanischen Mehrk\"orperproblems}\\
\textbf{Zeitschrift f\"ur Physik A 61 (1930) 126-148}

\bibitem{slater} J.C.~Slater\\
\textit{A Simplification of the Hartree-Fock Method}\\
\textbf{Phys.~Rev.~81 (1951) 385-390}


%general DFT stuff
% \bibitem{sutton} A.P.~Sutton\\
% \textit{Electronic Structure of Materials}\\
% \textbf{Oxford: Oxford University Press 1993}



% L.~Vitos\\
%\textit{Computational Quantum Mechanics for Materials Engineers}\\
%\textbf{London: Springer 2007}

\bibitem{NATODFT} E.K.U.~Gross and R.M.~Dreizler\\
\textit{Density Functional Theory}
\textbf{Berlin: Springer 1990}

% \bibitem{DFTbird} K.~Capelle\\
% \textit{A Bird's-Eye-View of Density-Functional Theory}\\
% \textbf{ArXiv}

\bibitem{hktheory} P.~Hohenberg and W.~Kohn\\
\textit{Inhomogeneous Electron Gas}\\
\textbf{Phys.~Rev.~336 (1964) B864-B871}

\bibitem{kohnsham} W.~Kohn and L.J.~Sham\\
\textit{Self-Consistent Equations Including Exchange and Correlation Effects}\\
\textbf{Phys.~Rev.~140 (1965) A1133-A1138}

%The references to the exchange and correlation approximations implemented in VASP are:
%Local Density Approximation (LDA)
\bibitem{LDA} J.~P.~Perdew and A.~Zunger\\
\textit{Self-Interaction Correction to Density-Functional Approximations for Many-Electron Systems}\\
\textbf{Phys.~Rev.~B 23 (1981) 5048}


%Generalized Gradient Approximation PBE (GGA-PBE)
\bibitem{GGA1} J.~P.~Perdew, K.~Burke, and M.~Ernzerhof\\
\textit{Generalized Gradient Approximation Made Simple}\\
\textbf{Phys.~Rev.~Lett.~77 (1996) 3865}

\bibitem{GGA2} J.~P.~Perdew, K.~Burke, and M.~Ernzerhof\\
\textit{Erratum: Generalized Gradient Approximation Made Simple}\\
\textbf{Phys.~Rev.~Lett.~78 (1997) 1396}


%Generalized Gradient Approximation PW91 (GGA-PW91)
\bibitem{PW91} J.P.~Perdew and Y.~Wang\\
\textit{Accurate and Simple Analytic Expression of the Electron-Gas Correlation Energy}\\
\textbf{Phys.~Rev.~B 45 (1992) 13244-13249}

\bibitem{PW1} J.~P.~Perdew et.al.\\
\textit{Atoms, Molecules, Solids, and Surfaces: Applications of the Generalized Gradient Approximation for Exchange and Correlation}\\
\textbf{Phys.~Rev.~B 46 (1992) 6671}

\bibitem{PW2} J.~P.~Perdew et.al.\\
\textit{Erratum: Atoms, Molecules, Solids, and Surfaces: Applications of the Generalized Gradient Approximation for Exchange and Correlation}\\
\textbf{Phys.~Rev.~B 48 (1993) 4978}

%Depending on the potentials used you should also include the following citations:
%Ultra-soft pseudopotentials
\bibitem{USP1} D.~Vanderbilt\\
\textit{Soft Self-Consistent Pseudopotentials in a Generalized Eigenvalue Formalism}\\
\textbf{Phys.~Rev.~B 41 (1990) 7892}

\bibitem{USP2} G.~Kresse and J.~Hafner\\
\textit{Norm-Conserving and Ultrasoft Pseudopotentials for First-Row and Transition Elements}\\
\textbf{J.~Phys.: Condens.~Matter 6 (1994) 8245}

%PAW potentials
\bibitem{PAW1} P.~E.~Bl\"ochl\\
\textit{Projector Augmented-Wave Method}\\
\textbf{Phys.~Rev.~B 50 (1994) 17953}

\bibitem{PAW2} G.~Kresse and D.~Joubert\\
\textit{From Ultrasoft Pseudopotentials to the Projector Augmented-Wave Method}\\
\textbf{Phys.~Rev.~B 59 (1999) 1758}



%general VASP stuff
\bibitem{VASP1} G.~Kresse and J.~Hafner\\
\textit{Ab Initio Molecular Dynamics for Liquid Metals}\\
\textbf{Phys.~Rev.~B 47 (1993) 558}

\bibitem{VASP2} G.~Kresse and J.~Hafner\\
\textit{Ab Initio Molecular-Dynamics Simulation of the Liquid-Metal-Amorphous-Semiconductor Transition in Germanium}\\
\textbf{Phys.~Rev.~B 49 (1994) 14251}

\bibitem{VASP3} G.~Kresse and J.~Furthm\"u ller\\
\textit{Efficiency of Ab-Initio Total Energy Calculations for Metals and Semiconductors Using a Plane-Wave Basis Set}\\
\textbf{Comput.~Mat.~Sci. 6 (1996) 15}



%integration techniques
\bibitem{imprLT} P.~Bl\"ochl, O.~Jepsen and O.K.~Andersen\\
\textit{Improved Tetrahedron Method for Brillouin Zone Integrations}\\
\textbf{Phys.~Rev.~B 49 (1994) 16223-16234}

\bibitem{LinTet1} O.~Jepson [sic] and O.K.~Anderson [sic]\\
\textit{The Electron Structure of HCP Ytterbium}\\
\textbf{Solid State Comm.~9 (1971) 1763-1767}

\bibitem{LinTet2} G.~Lehmann and M.~Taut\\
\textit{On the Numerical Calculation of the Density of States and Related Properties}\\
\textbf{Phys.~Stat.~Sol.~54 (1972) 469-477}

\bibitem{FermiSmear} N.~Ashcroft\\
\textit{Thermal Properties of the Inhomogeneous Electron Gas}\\
\textbf{Phys.~Rev.~137 (1965) A1441–A1443}

\bibitem{MethPaxt} M.~Methfessel and A.T.~Paxton\\
\textit{High-Precision Sampling for Brillouin-Zone Integration in Metals}\\
\textbf{Phys.~Rev.~B 40 (1989) 3616}



%Birch-Murnaghan EOS
\bibitem{BMEOS1} F.D.~Murnaghan\\
\textit{The Compressibility of Media under Extreme Pressures}\\
\textbf{Proc.~Nat.~Acad.~Sci.~30 (1944) 244-247}

\bibitem{BMEOS2} F.~Birch\\
\textit{Finite Elastic Strain of Cubic Crystals}\\
\textbf{Phys.~Rev.~71 (1947) 809-824}

%KKR stuff
\bibitem{laszlo1} J.~Zabloudil, R.~Hammerling, L.~Szunyogh and P.~Weinberger\\
\textit{Electron Scattering in Solid Matter - A Theoretical and Computational Treatise}\\
\textbf{Berlin: Springer 2005}

\bibitem{eltheorybook} U.~Mizutani\\
\textit{Introduction to Electron Theory of Metals}\\
\textbf{Cambridge: Cambridge University Press 1996}

\bibitem{ziman} J.M.~Ziman\\
\textit{Principles of the Theory of Solids}\\
\textbf{Cambridge: Cambridge University Press 1972}

\bibitem{Korringa} J.~Korringa\\
\textit{On the calculation of the energy of a Bloch wave in a metal }\\
\textbf{Physica 13 (1947) 392-400}

\bibitem{rostoker} W.~Kohn and N.~Rostoker\\
\textit{Solution of the Schrödinger Equation in Periodic Lattices with an Application to Metallic Lithium}\\
\textbf{Phys.~Rev.~94 (1954) 1111-1120}




%Hubbard model stuff
\bibitem{Hubbard} J.~Hubbard\\
\textit{Electron Correlations in Narrow Energy Bands}\\
\textbf{Proc.~Roy.~Soc.~A.: Math.~Phys.~Sci. 276 (1963) 238-257}

\bibitem{Hirsch} J.E.~Hirsch\\
\textit{Two-Dimensional Hubbard Model: Numerical Simulation Study}\\
\textbf{Phys.~Rev.~B 31 (1985) 4403-4419}

\bibitem{Tasaki} H.~Tasaki\\
\textit{The Hubbard Model - An Introduction and Selected Rigorous Results}\\
\textbf{J.~Phys.: Condensed Matter 10 (1998) 4353-4378}

\bibitem{Mielke} A.~Mielke\\
\textit{Ferromagnetism in the Hubbard model on Line Graphs and Further Considerations}\\
\textbf{J.~Phys.~A: Math.~Gen.~24 (1991) 3311 }

\bibitem{mielketasaki}  A.~Mielke and H.~Tasaki\\
\textit{Ferromagnetism in the Hubbard model. Examples from models with degenerate single-electron ground states}
\textbf{Commun.~Math.~Phys.~158 (1993) 341-371}

\bibitem{tasaki2} H.~Tasaki\\
\textit{Ferromagnetism in Hubbard Models}\\
\textbf{Phys.~Rev.~Lett.~75 (1995) 4678-4681}

\bibitem{tasaki3} H.~Tasaki\\
\textit{Stability of Ferromagnetism in Hubbard Models with Nearly Flat Bands}\\
\textbf{J.~Stat.~Phys.~84 (1996) 535-653}



%atomistic models
\bibitem{hubheis} P.W.~Anderson\\
\textit{New approach to the Theory of Superexchange Interactions}\\
\textbf{Phys.~Rev.~115 (1959) 2-13}

\bibitem{mattis} D.C.~Mattis\\
\textit{The Theory of Magnetism, Vols. ~I and II}\\
\textbf{New York: Springer-Verlag 1987}

\bibitem{mersennetwister} M.~Matsumoto and T.~Nishimura\\
\textit{Mersenne Twister: A 623-Dimensionally Equally Distributed Uniform Pseudo-Random Number Generator}\\
\textbf{ACM Trans.~Model.~Comput.~Simul.~8 (1998) 3-30}

%Fe-related
\bibitem{handbookmag} K.H.J.~Buschkow\\
\textit{Handbook of Magnetic Materials}\\
\textbf{New York: Elsevier 1995}

%FePt-related
\bibitem{cona_FePt} C.~Zhen et.al.\\
\textit{Effects of C Layer on the Microstructure and Magnetic Properties of FePt Recording Media Films}\\
\textbf{Materials Science and Engineering B 129 (2006) 261–264}

\bibitem{V0FePt} C.W.~Hsu\\
\textit{Effect of Pt Underlayer on the Coercivity of FePt Sputtered Film}\\
\textbf{J.~All.~Comp.~449 (2008) 52–55}

\bibitem{PVexpFePt} Y.H.~Ko et.al.\\
\textit{Pressure-Volume Equation of State of FeAu and FePt}\\
\textbf{High Pressure Research 29 (2009) 800–805}

% \bibitem{Varanda} L.~Varanda, M.~Jafelicci Jr.~and M.~Imaizumi\\
% \textit{Temperature Dependence and Magnetocrystalline Anisotropy Studies of Self-Assembled L10-Fe$_{55}$Pt$_{45}$ Ferromagnetic Nanocrystals}\\
% \textbf{J.~Appl.~Phys.~101 (2007) 123918}

\bibitem{FMAFMFePt} G.~Brown et.al.\\
\textit{Competition Between Ferromagnetism and Antiferromagnetism in FePt}\\
\textbf{Phys.~Rev.~B 68 (2003) 052405}

% \bibitem{Imaizumi} M.~Imaizumi et.al.\\
% \textit{Structural Phase Transition Study of FePt Alloys using Ab Initio Calculation}\\
% \textbf{Mat.~Sci.~Eng. A 521-522 (2009) 167-168}

% \bibitem{fept_struc} M.~Koslowski et.al.\\
% \textit{Atomic Ordering in Nano-Layered FePt}\\
% \textbf{Intermetallics 17 (2009) 907-913}
%
% \bibitem{Sun} S.~Sun et.al.\\
% \textit{Monodisperse FePt Nanoparticles and Ferromagnetic FePt Nanocrystal Superlattices}\\
% \textbf{Science 287 (2000) 1989-1992}

%CoPt-related
\bibitem{cona_CoPt} W.M.~Liao et.al.\\
\textit{Ordering enhancement of Cu underlayer on CoPt thin films}\\
\textbf{J.~Mag.~Mag.~Mat.~272–276 (2004) 2175–2177}

\bibitem{CoPtexp} W.M.~Liao et.al.\\
\textit{Thickness Dependence of Crystallographic and Magnetic Properties or L10-CoPt Thin Films}\\
\textbf{J.~Mag.~Magn.~Mat.~303 (2006) 243-246}


%FeRh-related
\bibitem{moruzzi} V.L.~Moruzzi and P.M.~Marcus\\
\textit{Antiferromagnetic-Ferromagnetic Transition in FeRh}\\
\textbf{Phys.~Rev.~B 46 (1992) 2864-2873}

\bibitem{Szajek} A.~Szajek and J.A.~Morkowski \\
\textit{The Electronic and Magnetic Properties of the Metamagnetic Ordered Alloy FeRh}\\
\textbf{Physica H 193 (1994) 81-91}


\bibitem{koenig} C.~Koenig\\
\textit{Self-Consistent Band Structure of Paramagnetic, Ferromagnetic and Antiferromagnetic Ordered FeRh}\\
\textbf{J.~Phys.~F: Met.~Phys.~12 (1982) 1123-1137}



\end{thebibliography}
\clearpage

\appendix
\section{Notation in the Text}
\begin{tabular}{| l | l |}
\hline
\textbf{Symbol} & \textbf{Definition}\\
\hline
$\hat{T}$ & electron kinetic energy operator\\
\hline
$\hat{W}$ & electron-electron interactions operator \\
\hline
$\hat{V}$ & operator representing the lattice potential plus any external potential acting on the electrons\\
\hline
$V_0$ & equilibrium unit cell volume\\
\hline
$B_0$ & bulk modulus\\
\hline
$\Psi_{mb}$ & electron many-body wavefunction\\
\hline
\end{tabular}


\section{Proof of the Hohenberg-Kohn Theorem}
\label{HKproof}
\paragraph{Statement 1:} A given potential leads to a given, unique density.\\*

A given potential yields a particular wavefunction (through the Schr\"odinger equation).  For non-degenerate systems there is only one ground state wavefunction.  From a given many-body wavefunction, we can obtain the density from
$$
 n(\mathbf{r}) = \int \left|\Psi(\mathbf{r}_1,\mathbf{r}_2,...,\mathbf{r}_N)\right|^2 d^3r_1 d^3 r_2...d^3 r_N
$$
Thus a given potential leads to a given, unique density.

\paragraph{Statement 2:} A given density leads to a given, unique potential.\\*

Assume that there are two potentials $V$ and $V'$ that lead to the same wavefunction $\Psi$.  The resulting two Schr\"odinger equations have the form
\begin{eqnarray}
 \left(\hat{T}+\hat{U}+\hat{V}\right) \Psi & = & E \Psi \\
 \left(\hat{T}+\hat{U}+\hat{V}'\right) \Psi & = & E' \Psi \\
\end{eqnarray}
This yields
$$
 \left(V-V'\right)\Psi = \left(E-E'\right) \Psi
$$
Since the eigenenergies $E$ and $E'$ are constants, the potentials $V$ and $V'$ can only differ by a constant.  Thus two different potentials must lead to two different wavefunctions.\\*

Assume that there are two wavefunctions $\Psi$ and $\Psi'$ that lead to the same density $n(\mathbf{r})$.  As proved above the two wavefunctions must relate to two different potentials $V$ and $V'$.  Thus the Schr\"odinger equations become
\begin{eqnarray}
 \hat{H}\Psi & = & \left(\hat{T}+\hat{U}+\hat{V}\right) \Psi = E \Psi \\
 \hat{H}'\Psi & = & \left(\hat{T}+\hat{U}+\hat{V}'\right) \Psi' = E' \Psi' \\
\end{eqnarray}
The ground state energy can be obtained
$$
 E=\langle \Psi | \hat{H} | \Psi\rangle <  \langle \Psi' | \hat{H} | \Psi' \rangle
$$
Writing $\hat{H} = \hat{H}' + \hat{V}-\hat{V}'$, which gives
$$
 E < E' + \int n(\mathbf{r}) \left(V(\mathbf{r}-V'(\mathbf{r})\right) d^3r
$$
Starting from $\hat{H}$ and $\Psi$ gives (since we have assumed the different wavefunctions yield the same density)
$$
 E < E' + \int n(\mathbf{r}) \left(V'(\mathbf{r}-V(\mathbf{r})\right) d^3r = E -\int n(\mathbf{r}) \left(V(\mathbf{r}-V'(\mathbf{r})\right) d^3r
$$
Adding these two equations gives
$$
 E+E'<E+E'
$$
which is a contradiction.  Thus the assumption that the two different wavefunctions yield the same density must be wrong.  Thus a given potential yields a given, unique density.

\section{Analytical Solution for the Landau-Lifshitz-Gilbert Equation}
\label{LLGanalytical}
The Landau-Lifshitz-Gilbert (LLG) equation can be solved exactly for a single spin with the applied field $\mathbf{H}=H_z\hat{\mathbf{z}}$.  This derivation was originally presented by Johan Mentink at Nijmegen University and explained to me by Joe Barker of the University of York.\\*

For a single spin $\mathbf{S}$ at some angle to the applied field $\mathbf{H}=H_z\hat{\mathbf{z}}$, the LLG equation becomes
$$
\label{LLGdelta}
 \left(1+\alpha^2\right)\frac{d \mathbf{S}}{d t} = - \gamma\mathbf{S}\times \left(\mathbf{H} + \alpha \mathbf{S} \times \mathbf{H}\right)
$$
The LLG equation can be transformed into the rotating frame (which is referred to by primed coordinates $\left(\hat{\mathbf{x}},\hat{\mathbf{y}},\hat{\mathbf{z}}\right)$) which rotates at angular velocity $\mathbf{\Omega}$.  The rate of change of any axis $\hat{\mathbf{u}}$ in the non-rotating frame becomes
$$
 \frac{d\hat{\mathbf{u}}}{dt} = \Omega \times \hat{\mathbf{u}}
$$
Expanding the time derivative,
$$
\label{rotdelta}
 \frac{d\mathbf{S}}{dt} = \frac{dS_x}{dt}  + \frac{d\hat{\mathbf{x}}}{dt} S_x + \frac{dS_y}{dt} +  \frac{d\hat{\mathbf{y}}}{dt} S_y + \frac{dS_z}{dt} + \frac{d\hat{\mathbf{z}}}{dt} S_z = \left(\frac{d\mathbf{S}}{dt}\right)_0 + \mathbf{\Omega} \times \hat{\mathbf{u}}
$$
where $\left(\frac{d\mathbf{S}}{dt}\right)_0$ is the rate of change in the non-rotating frame (as in this frame the coordinate axes do not change with time).  Combining equations \ref{LLGdelta} and \ref{rotdelta} gives the LLG equation in the rotating frame.
\begin{eqnarray}
 \left(1+\alpha^2\right)\left(\left(\frac{d \mathbf{S}}{d t}\right)_0 + \mathbf{\Omega} \times \mathbf{S} \right) & = & - \gamma\mathbf{S}\times \left(\mathbf{H} + \alpha \mathbf{S} \times \mathbf{H}\right)\\
  \left(1+\alpha^2\right)\left(\frac{d \mathbf{S}}{d t}\right)_0 & = & - \gamma\mathbf{S}\times \left(\mathbf{H} -\frac{(1+\alpha^2)}{\gamma} \mathbf{\Omega} \right) - \alpha \gamma \mathbf{S} \times \left(\mathbf{S} \times \mathbf{H}\right) \label{LLG2}
\end{eqnarray}
The second line follows from the commutation relation $\mathbf{S}\times \mathbf{\Omega} = - \mathbf{\Omega} \times \mathbf{S}$.  Let the rotating frame move with an angular velocity equal to the spin precession velocity
$$
 \mathbf{\Omega} = \frac{\gamma}{(1+\alpha^2)}\mathbf{H}
$$
This simplifies equation \ref{LLG2} to become
$$
 \left(1+\alpha^2\right)\left(\frac{d \mathbf{S}}{d t}\right)_0 =  - \alpha \gamma\mathbf{S} \times \left(\mathbf{S} \times \mathbf{H}\right)
$$
Substituting $\mathbf{H}=H_z\hat{\mathbf{z}}$ gives
$$
 \left(1+\alpha^2\right)\left(\frac{d \mathbf{S}}{d t}\right)_0 =  - \alpha \gamma \left(\begin{array}{c}
                                                                                                            S_x\\
                                                                                                            S_y\\
                                                                                                            S_z\\
                                                                                                           \end{array}\right) \times \left(\begin{array}{c}
                                                                                                            S_yH_z\\
                                                                                                            -S_xH_z\\
                                                                                                            0\\
                                                                                                           \end{array}\right)
$$
Writing this in components,
\begin{eqnarray}
\left(\frac{d S_x}{dt}\right)_0 &=&  - \frac{\alpha \gamma H_z}{ \left(1+\alpha^2\right)} S_xS_z \label{Sxeq} \\
\left(\frac{d S_y}{dt}\right)_0 &=&  - \frac{\alpha \gamma H_z}{ \left(1+\alpha^2\right)} S_yS_z \label{Syeq} \\
\left(\frac{d S_z}{dt}\right)_0 &=&  \frac{\alpha \gamma H_z}{\left(1+\alpha^2\right)} \left(S_x^2+S_y^2\right) = \frac{\alpha \gamma H_z}{\left(1+\alpha^2\right)} \left(1-S_z^2\right)  \label{Szeq}
\end{eqnarray}
Equation \ref{Szeq} can be solved to give (in the non-rotating frame)
$$
 S_z(t) = \tanh\left(\frac{\alpha \gamma}{1+\alpha^2} t\right)
\label{Szsoln}
$$
Substituting equation \ref{Szsoln} into equations \ref{Sxeq} and \ref{Syeq} renders them solvable and yields
\begin{eqnarray}
  S_x(t) & = & S_x(0) \cos(\Omega t) \text{sech}\left(\frac{\alpha \gamma}{1+\alpha^2} t\right)\\
  S_y(t) & = & S_y(0) \sin(\Omega t) \text{sech}\left(\frac{\alpha \gamma}{1+\alpha^2} t\right)
\end{eqnarray}








\section{Easy Axes in the xy Plane}
The existence of easy axes in the xy plane follows from the properties of the in-plane nearest-neighbour separation vectors $\mathbf{R}$.  These in-plane vectors $\mathbf{R}$ are invariant under the xy mirror operation $S_{xy}$
$$
 S_{xy}=\left(\begin{array}{c c c}
         1 & 0 & 0\\
         0 & 1 & 0\\
         0 & 0 & -1\\
        \end{array}\right)
$$
It follows that
$$
 S_{xy}\mathbf{R} = \mathbf{R} \Rightarrow  S_{xy} J(\mathbf{R}) S_{xy}^T = J(\mathbf{R})
$$
where $J(\mathbf{R})$ is the exchange matrix corresponding to a particular in-plane vector $\mathbf{R}$.  For a general exchange matrix,
$$
  J=\left(\begin{array}{c c c}
         a & b & c\\
         d & e & f\\
         g & h & i\\
        \end{array}\right)
$$
it follows that
$$
J(\mathbf{R})= S_{xy} J(\mathbf{R}) S_{xy}^T \Rightarrow c=f=g=h=0
$$
Thus, the exhange matrix for in-plane neighbours must be of the form
$$
  J=\left(\begin{array}{c c c}
         a & b & 0\\
         d & e & 0\\
         0 & 0 & i\\
        \end{array}\right)
$$
For a simple cubic lattice, the in-plane nearest neighbours are
$$
\mathbf{R}_{01}= \left(\begin{array}{c c c}
         1\\
         0\\
         0\\
        \end{array}\right) \ \ \ \ \ \  \mathbf{R}_{02}=\left(\begin{array}{c c c}
         0\\
         -1\\
         0\\
        \end{array}\right) \ \ \ \ \ \  \mathbf{R}_{03}=\left(\begin{array}{c c c}
         -1\\
         0\\
         0\\
        \end{array}\right) \ \ \ \ \ \  \mathbf{R}_{04}=\left(\begin{array}{c c c}
         0\\
         1\\
         0\\
        \end{array}\right)
$$
These four different nearest-neighbour vectors can be transformed into each other by a rotation of $\pi/2$ around the $z$-axis, as represented by the matrix $C_4$
$$
 C_4 = \left(\begin{array}{c c c}
         0 & 1 & 0\\
         -1 & 0 & 0\\
         0 & 0 & 1\\
        \end{array}\right)
$$
giving
\begin{eqnarray}
& & C_4 \mathbf{R}_{01}= \mathbf{R}_{01} \ \ \ \ \ \  C_4\mathbf{R}_{02}=\mathbf{R}_{04} \ \ \ \ \ \ C_4\mathbf{R}_{03}=\mathbf{R}_{03} \ \ \ \ \ \  C_4\mathbf{R}_{04}=\mathbf{R}_{02}\\
& \Rightarrow & S_{xz} J_{01} S_{xz}^T= J_{01} \ \ \ \ \ \  S_{xz} J_{02} S_{xz}^T =J_{04} \ \ \ \ \ \ S_{xz}J_{03} S_{xz}^T=J_{03} \ \ \ \ \ \  S_{xz} J_{04}S_{xz}^T=J_{02} \label{spinC4} \\
\end{eqnarray}
Furthermore, considering the xz mirror operation
$$
 S_{xz}=\left(\begin{array}{c c c}
               1 & 0 & 0\\
               0 & -1 & 0\\
               0 & 0 & 1\\
              \end{array}\right)
$$
giving
\begin{eqnarray}
& & S_{xz} \mathbf{R}_{01}= \mathbf{R}_{01} \ \ \ \ \ \  S_{xz}\mathbf{R}_{02}=\mathbf{R}_{03} \ \ \ \ \ \ S_{xz}\mathbf{R}_{03}=\mathbf{R}_{04} \ \ \ \ \ \  S_{xz}\mathbf{R}_{04}=\mathbf{R}_{01}\\
& \Rightarrow & C_4 J_{01} C_4^T= J_{02} \ \ \ \ \ \  C_4 J_{02} C_4^T =J_{03} \ \ \ \ \ \ C_4 J_{03} C_4^T=J_{04} \ \ \ \ \ \  C_4 J_{04}C_4^T=J_{01} \label{spinSxz}\\
\end{eqnarray}
By letting
$$
 J_{01}= \left(\begin{array}{c c c}
     a & b & 0\\
     d & e & 0\\
     0 & 0 & i\\
    \end{array}\right)
$$
it can be shown (using equations \ref{spinC4} and \ref{spinSxz}) $a=e$ and $b=d$. Thus the matrices become (in general $a\neq e$)
$$
 J_{01}= \left(\begin{array}{c c c}
     a & 0 & 0\\
     0 & e & 0\\
     0 & 0 & i\\
    \end{array}\right) \ \ \ \ \ \ \  J_{02}= \left(\begin{array}{c c c}
     e & 0 & 0\\
     0 & a & 0\\
     0 & 0 & i\\
    \end{array}\right) \ \ \ \ \ \ \   J_{03}= \left(\begin{array}{c c c}
     a & 0 & 0\\
     0 & e & 0\\
     0 & 0 & i\\
    \end{array}\right) \ \ \ \ \ \ \  J_{04}= \left(\begin{array}{c c c}
     e & 0 & 0\\
     0 & a & 0\\
     0 & 0 & i\\
    \end{array}\right)
$$


\section{FORTRAN Code for the Hubbard Model}
\label{HubbardFORTRAN}

\section{Band Structure Derivations for the Tight-Binding Model}
\label{TBapp}
For instructional purposes, this appendix lists the analytical derivations of the tight-binding Hamiltonian matrices for a number of configurations.
\paragraph{Three-Dimensional Lattice}
\ \\
For a three-dimensional lattice with periodic boundary conditions, nearest neighbours are separated by $\left(\mathbf{r}_m-\mathbf{r}_n\right) = (\pm a,0,0), (0, \pm b, 0), (0,0,\pm c)$.  Thus the energy becomes
$$
  E(\mathbf{k})  = -E_0 -t \left(e^{-ik_x a} + e^{ik_x a} + e^{-ik_y b} + e^{ik_y b} + e^{-ik_z c} + e^{ik_z c} \right) = - E_0 - 2t\left(\cos(k_x a) + \cos(k_y b)+ \cos(k_z c)\right)
\label{threeDmonatomic}
$$


\paragraph{Body-Centred Cubic Lattice}
\ \\
For a body-centred cubic lattice, the nearest neighbour to any corner atom are the central atoms of the surrounding eight unit cells.  There are thus eight nearest neighbours at  $\left(\mathbf{r}_m-\mathbf{r}_n\right) = (\pm a,\pm a,\pm a)$.  The energy thus becomes
\begin{eqnarray}
  E(\mathbf{k})  & = & -E_0 -t \left(e^{i(k_x + k_y + k_z) a/2} + e^{-i(k_x + k_y + k_z) a/2}  + e^{i(-k_x + k_y + k_z) a/2} + e^{-i(-k_x + k_y + k_z) a/2}\right) \nonumber  \\
       & &    -t \left( e^{i(k_x - k_y + k_z) a/2} + e^{-i(k_x - k_y + k_z) a/2} + e^{i(k_x + k_y - k_z) a/2} + e^{-i(k_x + k_y - k_z) a/2} \right)\\
       & = & -E_0 -2t\left[\cos\left(\frac{k_x + k_y + k_z}{2}a\right) + \cos\left(\frac{-k_x + k_y + k_z}{2}a\right) + \cos\left(\frac{k_x - k_y + k_z}{2}a\right) + \cos\left(\frac{k_x + k_y - k_z}{2}a\right)\right] \nonumber
\end{eqnarray}

\paragraph{Face-Centred Cubic Lattice}
\ \\
For a face-centred lattice, the nearest neighbour to each corner atom are the twelve atoms on the faces.  These have positions
\begin{eqnarray}
 \left(a/2, a/2, 0\right) & \ \ \ \ \ \ & \left(-a/2, -a/2, 0\right)\\
 \left(a/2, 0, a/2\right) & \ \ \ \ \ \ & \left(-a/2, 0, -a/2\right)\\
 \left(0, a/2, a/2\right) & \ \ \ \ \ \ & \left(0, -a/2, -a/2\right)\\
 \left(a/2, -a/2, 0\right) & \ \ \ \ \ \ & \left(a/2, -a/2, 0\right)\\
 \left(-a/2, 0, a/2\right) & \ \ \ \ \ \ & \left(a/2, 0, -a/2\right)\\
 \left(0, -a/2, a/2\right) & \ \ \ \ \ \ & \left(0, a/2, -a/2\right)\\
\end{eqnarray}
The energy thus becomes
\begin{eqnarray}
  E(\mathbf{k})  & = & -E_0 -t \left(e^{i(k_x + k_y) a/2} + e^{-i(k_x + k_y) a/2}  + e^{i(k_x + k_z) a/2} + e^{-i(k_x + k_z) a/2} +  e^{i(k_y + k_z) a/2} + e^{-i(k_y + k_z) a/2}\right) \\
       & &   -t \left(e^{i(k_x - k_y) a/2} + e^{-i(k_x - k_y) a/2}  + e^{i(k_x - k_z) a/2} + e^{-i(k_x - k_z) a/2} +  e^{i(k_y - k_z) a/2} + e^{-i(k_y - k_z) a/2}\right)\\
       & = & -E_0 -2t\left(\cos\left((k_x + k_y) a/2\right) + \cos\left((k_x + k_z) a/2\right) + \cos\left((k_y + k_z) a/2\right)\right)\\
       & &    -2t \left( \cos\left((k_x - k_y) a/2\right) + \cos\left((k_x - k_z) a/2\right) + \cos\left((k_y - k_z) a/2\right)\right)\\
\end{eqnarray}


\subsection{One-Dimensional Chain with $N$ Atoms per Unit Cell}

The tight-binding Hamiltonian for a one-dimensional triatomic chain is (assuming $E_A = E_B = E_C = E$)
\begin{eqnarray}
 \hat{H} & = & E \sum_i \left(c_{A(i)}^{\dagger}c_{A(i)} + c_{B(i)}^{\dagger}c_{B(i)}+ c_{C(i)}^{\dagger}c_{C(i)}\right)\\
   & & - t \sum_i \left(c_{A(i)}^{\dagger}c_{B(i)} + c_{B(i)}^{\dagger}c_{A(i)} + c_{C(i)}^{\dagger}c_{B(i)} + c_{B(i)}^{\dagger}c_{C(i)} + c^{\dagger}_{C(i-1)}c_{A(i)} + c^\dagger_{A(i+1)}c_{C(i)}\right) \nonumber
\end{eqnarray}
where all the rightward hops are written out explicitly and the possible leftward hops are included in the hermitian conjugate (h.c.).  For a single unit cell with index  $i$,
\begin{itemize}
 \item the one electron in $|k_{A(i)}\rangle$ can hop to $|k_{B(i)}\rangle$ or $|k_{B(i-1)}\rangle$
 \item the one electron in $|k_{B(i)}\rangle$ can hop to $|k_{A(i)}\rangle$ or $|k_{C(i)}\rangle$
 \item the one electron in $|k_{C(i)}\rangle$ can hop to $|k_{B(i)}\rangle$ or $|k_{A(i+1)}\rangle$
\end{itemize}
The Schr\"odinger equation thus becomes (setting $E=0$)
$$
\left(\begin{array}{c c c}
        0 & - t & -t e^{ika}\\
        -t & 0 & -t\\
        -te^{-ika} & -t & 0
       \end{array}\right)\left(\begin{array}{c}
                         u_A \\
			 u_B \\
			 u_C
                        \end{array}\right) = E(\mathbf{k})\left(\begin{array}{c}
                                                    u_A \\
						    u_B \\
					            u_C
                                                   \end{array}\right)
$$


Similarly, for a one-dimensional chain with six distinguishable atoms in its unit cell, the matrix becomes
$$
\left(\begin{array}{c c c c c c}
        0 & -t & 0 & 0 & 0 & -te^{ika}\\
        -t & 0 & -t & 0 & 0 & 0\\
        0 & -t & 0 & -t & 0 & 0  \\
        0 & 0 & -t & 0 & -t & 0\\
        0 & 0 & 0 & -t & 0 & -t\\
        -te^{-ika} & 0 & 0 & 0 & -t & 0
       \end{array}\right)
$$
So the matrix develops to greater numbers $N$ of atoms in the unit cell as follows.  The off-diagonal t's represent the hopping back and forth between neighbouring atoms inside the unit cell of the chain.  The two corner elements $-te^{\pm ika}$ represent hopping between the last (first) atom of the unit cell and the first (last) atom of the following (preceding) unit cell.  For a monatomic chain, the terms get 'squeezed' together into one big matrix element $\hat{H}_{11} = - t\left(e^{ika} + e^{-ika}\right)$.

\subsection{Two-Dimensional Quatratomic Sheet}
Consider a sheet with a four distinguishable atoms per unit cell.  The following sixteen hops can occur.
\begin{itemize}
 \item the one electron in $|k_{A(i,j)}\rangle$ can hop to $|k_{B(i,j)}\rangle$, $|k_{B(i-1,j)}\rangle$, $|k_{C(i,j)}\rangle$ or $|k_{C(i,j-1)}\rangle$
 \item the one electron in $|k_{B(i,j)}\rangle$ can hop to $|k_{A(i,j)}\rangle$, $|k_{A(i+1,j)}\rangle$ $|k_{D(i,j)}\rangle$ or $|k_{D(i,j-1)}\rangle$
 \item the one electron in $|k_{C(i,j)}\rangle$ can hop to $|k_{D(i,j)}\rangle$, $|k_{D(i-1,j)}\rangle$, $|k_{A(i,j)}\rangle$ or $|k_{A(i,j+1)}\rangle$
 \item the one electron in $|k_{D(i,j)}\rangle$ can hop to $|k_{C(i,j)}\rangle$, $|k_{C(i+1,j)}\rangle$ $|k_{B(i,j)}\rangle$ or $|k_{B(i,j+1)}\rangle$
\end{itemize}
The Hamiltonian thus becomes
\begin{eqnarray}
 \hat{H} & = & E \sum_i \left(c_{A(i)}^{\dagger}c_{A(i)} + c_{B(i)}^{\dagger}c_{B(i)}+ c_{C(i)}^{\dagger}c_{C(i)}+ c_{D(i)}^{\dagger}c_{D(i)}\right)\\
   & & - t \sum_{i,j} \left(c_{A(i,j)}^{\dagger}c_{B(i,j)} + c_{B(i,j)}^{\dagger}c_{A(i,j)} + c_{C(i,j)}^{\dagger}c_{D(i,j)} + c_{D(i,j)}^{\dagger}c_{C(i,j)}  \right) \nonumber \\
   & & - t \sum_{i,j} \left(c_{B(i,j)}^{\dagger}c_{D(i,j)} + c_{D(i,j)}^{\dagger}c_{B(i,j)} + c_{A(i,j)}^{\dagger}c_{C(i,j)} + c_{C(i,j)}^{\dagger}c_{A(i,j)} \right) \nonumber \\
   & & - t \sum_{i,j} \left(c^\dagger_{C(i+1,j)}c_{D(i,j)}+c_{D(i-1,j)}^{\dagger}c_{C(i,j)} + c_{A(i+1,j)}^{\dagger}c_{B(i,j)} + c_{B(i-1,j)}^{\dagger}c_{A(i,j)}\right)  \nonumber \\
   & & - t \sum_{i,j} \left(c^\dagger_{A(i,j+1)}c_{C(i,j)}+c_{C(i,j-1)}^{\dagger}c_{A(i,j)} + c_{D(i,j-1)}^{\dagger}c_{B(i,j)} + c_{B(i,j+1)}^{\dagger}c_{D(i,j)}\right)  \nonumber \\
\end{eqnarray}

The Hamiltonian matrix thus becomes (eigenvectors are of the form $(u_A, u_B, u_C, u_D)$)
$$
 \hat{H} = \left(\begin{array}{c c c c}
                  E_0 & -t\left(1+e^{ik_xa}\right) & -t\left(1+e^{ik_yb}\right) & 0\\
                  -t\left(1+e^{-ik_xa}\right)  & E_0 & 0 & -t\left(1+e^{ik_yb}\right)\\
		  -t\left(1+e^{-ik_yb}\right) & 0 & E_0 & -t\left(1+e^{ik_xa}\right)\\
                  0 & -t\left(1+e^{-ik_yb}\right) & -t\left(1+e^{-ik_xa}\right) & E_0\\
                 \end{array}\right)
$$

By symmetry, a 3x3 unit cell would have a Hamiltonian matrix of
$$
 \hat{H} = \left(\begin{array}{c c c c c c c c c}
                  E_0 & -t & -te^{ik_xa} & -t  & 0 & 0 & -te^{ik_yb} & 0 & 0\\
                  -t  & E_0 & -t & 0 & -t & 0 & 0 & -te^{ik_yb} & 0\\
		  -te^{-ik_xa} & -t & E_0 & 0 & 0 & -t & 0 & 0 & -te^{ik_yb}\\
                  -t & 0 & 0 & E_0& -t & -te^{ik_xa}  & -t & 0 & 0\\
                  0 & -t & 0 & -t & E_0 & -t & 0 & -t & 0\\
                  0 & 0 & -t & -te^{-ik_xa} & -t & E_0 & 0 & 0 & -t\\
                  -te^{-ik_yb} & 0 & 0 & -t & 0 & 0 & E_0 & -t &  -te^{ik_xa}  \\
                  0 & -te^{-ik_yb} & 0 & 0 & -t & 0 & -t & E_0 & -t\\
                  0 & 0 & -te^{-ik_yb} & 0 & 0 & -t & -te^{-ik_xa} & -t & E_0\\
                 \end{array}\right)
$$
Essentially this is tree times a diatomic chain (the tree diatomic chains are represented by the matrices along the diagonal).  Connecting the three diatomic chains into a 3x3 unit cell is represented by the two diagonal lines of $k_y$ hopping terms $-te^{\pm ik_yb}$ in the corner and the two off-diagonal $t$'s outside the 3x3 matrices.  Because the electrons can hop back and forth, the matrix will always be Hermitian.

\subsection{Three-Dimensional Octoatomic Cube}
Consider a 2x2x2 cube consisting of eight atoms A, B, C, D, E, F, G and H.  In each unit cell the following hops can occur.
\begin{itemize}
 \item the one electron in $|k_{A(i,j,k)}\rangle$ can hop to $|k_{B(i,j,k)}\rangle$, $|k_{B(i-1,j,k)}\rangle$, $|k_{C(i,j,k)}\rangle$, $|k_{C(i,j-1,k)}\rangle$, $|k_{E(i,j,k)}\rangle$, $|k_{E(i,j,k-1)}\rangle$
 \item the one electron in $|k_{B(i,j,k)}\rangle$ can hop to $|k_{A(i,j,k)}\rangle$, $|k_{A(i+1,j,k)}\rangle$, $|k_{D(i,j,k)}\rangle$, $|k_{D(i,j-1,k)}\rangle$, $|k_{F(i,j,k)}\rangle$, $|k_{F(i,j,k-1)}\rangle$
 \item the one electron in $|k_{C(i,j,k)}\rangle$ can hop to $|k_{D(i,j,k)}\rangle$, $|k_{D(i-1,j,k)}\rangle$, $|k_{A(i,j,k)}\rangle$, $|k_{A(i,j+1,k)}\rangle$, $|k_{G(i,j,k)}\rangle$, $|k_{G(i,j,k-1)}\rangle$
 \item the one electron in $|k_{D(i,j,k)}\rangle$ can hop to $|k_{C(i,j,k)}\rangle$, $|k_{C(i+1,j,k)}\rangle$, $|k_{B(i,j,k)}\rangle$, $|k_{B(i,j+1,k)}\rangle$, $|k_{H(i,j,k)}\rangle$, $|k_{H(i,j,k-1)}\rangle$
 \item the one electron in $|k_{E(i,j,k)}\rangle$ can hop to $|k_{F(i,j,k)}\rangle$, $|k_{F(i-1,j,k)}\rangle$, $|k_{G(i,j,k)}\rangle$, $|k_{G(i,j-1,k)}\rangle$, $|k_{A(i,j,k)}\rangle$, $|k_{A(i,j,k+1)}\rangle$
 \item the one electron in $|k_{F(i,j,k)}\rangle$ can hop to $|k_{E(i,j,k)}\rangle$, $|k_{E(i+1,j,k)}\rangle$, $|k_{H(i,j,k)}\rangle$, $|k_{H(i,j-1,k)}\rangle$, $|k_{B(i,j,k)}\rangle$, $|k_{B(i,j,k+1)}\rangle$
 \item the one electron in $|k_{G(i,j,k)}\rangle$ can hop to $|k_{H(i,j,k)}\rangle$, $|k_{H(i-1,j,k)}\rangle$, $|k_{E(i,j,k)}\rangle$, $|k_{E(i,j+1,k)}\rangle$, $|k_{C(i,j,k)}\rangle$, $|k_{C(i,j,k+1)}\rangle$
 \item the one electron in $|k_{H(i,j,k)}\rangle$ can hop to $|k_{G(i,j,k)}\rangle$, $|k_{G(i+1,j,k)}\rangle$, $|k_{F(i,j,k)}\rangle$, $|k_{F(i,j+1,k)}\rangle$, $|k_{D(i,j,k)}\rangle$, $|k_{D(i,j,k+1)}\rangle$
\end{itemize}

The Hamiltonian matrix for a three-dimensional octoatomic cube thus becomes
$$
\left(\begin{array}{c c c c c c c c}
                  E_0 & f(k_x a) & f(k_yb) & 0 & f(k_zc) & 0 & 0 & 0\\
                  f(-k_x a)  & E_0 & 0 & f(k_yb) & 0 & f(k_zc) & 0 & 0\\
		  f(-k_yb) & 0 & E_0 & f(k_x a) & 0 & 0 & f(k_zc) &0\\
                  0 & f(-k_yb) &  f(-k_x a) & E_0& 0 & 0  & 0&  f(k_zc)\\
                  f(-k_zc) & 0 &  0 & 0 & E_0 & f(k_x a) & f(k_yb) & 0\\
                  0 & f(-k_zc) & 0 & 0 &  f(-k_x a) & E_0 & 0 &  f(k_yb)\\
                  0 & 0 & f(-k_zc) & 0 &  f(-k_yb) & 0 & E_0 & f(k_x a) \\
                  0 & 0 & 0 & f(-k_zc) & 0 & f(-k_yb)  &  f(-k_x a) & E_0\\
                 \end{array}\right)
$$
where $f(\omega) = -t\left(1+e^{i\omega}\right)$.  There are two matrices corresponding to two 2x2 two-dimensional sheets along the diagonal of the matrix.  These two sheets are connected into a 2x2x2 cube by the off-diagonal $k_z c$ hopping terms $f(k_z c)$.

\section{Two-Dimensional Hubbard Code}
\begin{lstlisting}
 PROGRAM HubbardModel
 IMPLICIT NONE

!Implementation of the Hubbard Model

!!***************************************************************!!
!!								 !!
!! 		    DECLARATION OF VARIABLES			 !!
!!								 !!
!!***************************************************************!!
 INTEGER, PARAMETER :: dp=selected_real_kind(15,300)    !define double precision

!physical variables
 REAL(KIND=DP), PARAMETER :: pi = 4d0*ATAN(1d0)                       !pi
 COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: H_up,H_down         !total up and down Hamiltonian arrays
 COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: H_t,H_o             !partial Hamiltonian arrays:  transfer, Hubbard
 INTEGER, PARAMETER :: n_x=2, n_y=2                         !no of sites along x and y respectively
 REAL(KIND= DP),PARAMETER :: t=1d0,U_max=10d0                 !hopping amplitude, Hubbard U
 REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: U   !Hubbard U array (one Hubbard U for each site)
 REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: n_av_up, n_av_down, nup_tot, ndown_tot      !average densities
 REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: change_up, change_down  !change in average density
! REAL(KIND=DP) :: n_old_up, n_old_down                               !temporary variables to hold old densities while updating
 REAL(KIND=DP) :: delta_U                         !Hubbard U increment size
 REAL(KIND=DP) :: delta_k, k_latx, k_laty         !wave vectors, wave vector increment
 REAL(KIND=DP), PARAMETER :: a_x=1d0,a_y=1d0      !lattice parameters
 INTEGER :: n_elec, N_s       !total no of electrons, total no of atoms
 INTEGER :: n_up, n_do   !no of up and down electrons
 REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: weights  !k-space weight factors
 REAL(KIND=DP) :: kint
 REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: energy_onsite,totE  !on-site energies, integrand and integral

!program paramters
 INTEGER :: i,j,s,r,l        !do-loop counters
 INTEGER,PARAMETER :: N_iter=10         !max no of iterations
 REAL(KIND=DP) :: tol      !convergence tolerance
! INTEGER :: ierrop, ierron !error labels
! INTEGER :: itime    !time
 INTEGER :: nkx,nky       !no of k-values to integrate over

!lapack variables
 INTEGER :: lwork, info                                      !LAPACK routine parameters
 REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: work, rwork        !LAPACK routine arrays
 REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: W_up, W_down       !eigenvalues for up and down Hamiltonians
 EXTERNAL dsyevd                                                !access the LAPACK source files


 !!***************************************************************!!
 !!								 !!
 !! 			MAIN PROGRAM				 !!
 !!								 !!
 !!***************************************************************!!



 !initialise problem
 N_s = n_x*n_y          !total no of sites in unit cell
 k_latx = -pi           !initialising wave vector along x
 k_laty = -pi           !initialising wave vector along y
 n_elec =  N_s          !no of electrons (system is half-filled)
 delta_U = 0.1d0        !Hubbard U increment size
 n_do =n_elec/2        !no of down electrons
 n_up =n_elec/2        !no of up electrons
 !IF(MOD(n_elec,2) .NE. 0) THEN
 !   n_up = n_up +1
 !END IF

 WRITE(*,*) 'No of electrons', N_s
 WRITE(*,*) 'No of down electrons ', n_do
 WRITE(*,*) 'No of up electrons ', n_up

!initialise algorithm
 tol = 0.000001d0   !tolerance for convergence
 nkx = 100    !no of iterations over k_x
 nky = nkx   !no of iterations over k_y
 delta_k = (2d0*pi)/real(nkx)

!initialising the k-space weighting factors
 ALLOCATE(weights(1:nky+1,1:nkx+1))  !need one weight value for each point in k-space
 weights(:,:) = 1d0
 weights(1,:) = 0.5d0
 weights(nky+1,:) = 0.5d0
 weights(:,nkx+1) = 0.5d0
 weights(:,1) = 0.5d0
 weights(1,1) = 0.25d0
 weights(1,nkx+1) = 0.25d0
 weights(nky+1,1) = 0.25d0
 weights(nky+1,nkx+1) = 0.25d0
 kint = pi/(a_x*real(nkx,dp))
! weights(:,:) = weights(:,:)*4d0*kint*kint/(4d0*pi*pi)

 WRITE(*,*) 'nkx = ', nkx
 WRITE(*,*) 'nky = ', nky

!initialise lapack routine
 lwork = 2 + 6*N_s + 2*(N_s**2)   !lapack parameter

! allocate space to arrays
 ALLOCATE(H_down(1:N_s,1:N_s))
 ALLOCATE(H_up(1:N_s,1:N_s))
 ALLOCATE(W_up(1:N_s))
 ALLOCATE(W_down(1:N_s))
 ALLOCATE(work(1:lwork))
 ALLOCATE(rwork(1:3*N_s-2))
 ALLOCATE(H_t(1:N_s,1:N_s))
 ALLOCATE(H_o(1:N_s,1:N_s))
 ALLOCATE(n_av_up(1:N_s))
 ALLOCATE(n_av_down(1:N_s))
 ALLOCATE(change_up(1:N_s))
 ALLOCATE(change_down(1:N_s))
 ALLOCATE(U(1:N_s))
 ALLOCATE(nup_tot(1:N_s))
 ALLOCATE(ndown_tot(1:N_s))
 ALLOCATE(energy_onsite(1:N_s))
 ALLOCATE(totE(1:N_s))

!initialise arrays
 H_up(:,:) = 0d0
 H_down(:,:) = 0d0
 H_t(:,:) = 0d0
 H_o(:,:) = 0d0
 W_up(:) = 0d0
 W_down(:) = 0d0
 n_av_up(:) = 0d0    !average spin-up density initialised to zero
 n_av_down(:) = 0d0  !average spin-down density initialised to zero
 U(:) = U_max
 ndown_tot(:) = 0d0
 nup_tot(:) = 0d0

 !open file to store eigenenergies for each wave vector
 OPEN(UNIT=31, FILE="bands.dat", STATUS="new", POSITION="append", ACTION="write" )
 OPEN(UNIT=32, FILE="occupancy.dat", STATUS="new", POSITION="append", ACTION="write" )
 OPEN(UNIT=33, FILE="contourbands.dat", STATUS="new", POSITION="append", ACTION="write" )
 OPEN(UNIT=34, FILE="totE.dat", STATUS="new", POSITION="append", ACTION="write" )


 DO s = 1,N_iter    !iterate until convergence
    !down-spin scan through k-space
    ndown_tot(:)=0d0
    k_latx = -pi
    DO l = 1,nkx+1
       k_laty = -pi     !=k_latx  !start at k_latx do get a wedge in BZ
       DO r = 1,nky+1   !iterate over different wave vectors
          !DO i = 1,N_s
          !   U(i) = rand(itime)            !initialise Hubbard U
          !END DO

          !calculate Hamiltonian parts for up Hamiltonian
          CALL hopping(t,H_t,N_s,a_x,a_y,k_latx,k_laty)
          CALL hubU(U,H_o,N_s,nup_tot)

          !evaluate total down Hamiltonian
          DO i = 1,N_s
             DO j=1,N_s
                H_down(j,i) = H_t(j,i) + H_o(j,i)
             END DO
          END DO

          !find eigenvalues and eigenfunctions of up Hamiltonian
          CALL zheev('V', 'U', N_s, H_down, N_s, W_down, work, lwork, rwork, INFO)
          if(INFO .NE. 0) then
             print *, "Call to LAPACK routine zheev failed"
             stop
          endif

          n_av_down(:) = 0d0

          !calculate average spin-down density
          DO j=1,n_do
             DO i = 1,N_s
                n_av_down(i) = n_av_down(i) + abs(H_down(i,j))**2
             END DO
          END DO

          ndown_tot(:) = ndown_tot(:) + (n_av_down(:)*weights(l,r))!/(DBLE((nkx-1)*(nky-1)))!*weights(l,r)
          k_laty = k_laty + delta_k
       END DO
       k_latx = k_latx +delta_k
    END DO

    ndown_tot(:) = ndown_tot(:)/(DBLE(nkx*nky))

    !up-spin integral over k-space
    nup_tot(:) = 0d0
    k_latx =-pi
    DO l=1,nkx+1
       k_laty=-pi
       DO r = 1,nky+1
          !calculate Hamiltonian parts for down Hamiltonian
          CALL hopping(t,H_t,N_s,a_x,a_y,k_latx,k_laty)
          CALL hubU(U,H_o,N_s,ndown_tot)

          !evaluate total up Hamiltonian
          DO i = 1,N_s
             DO j=1,N_s
                H_up(j,i) = H_t(j,i) + H_o(j,i)
             END DO
          END DO

          !find eigenvalues and eigenfunctions of Hamiltonian matrix
          CALL zheev('V', 'U', N_s, H_up, N_s, W_up, work, lwork, rwork, INFO)
          if(INFO .NE. 0) then
             print *, "Call to LAPACK routine zheev failed"
             stop
          endif

          !calculate average spin-up density
          n_av_up(:) = 0d0
          DO j = 1,n_up
             DO i = 1,N_s
                n_av_up(i) = n_av_up(i) + abs(H_up(i,j))**2
             END DO
          END DO

          nup_tot(:) = nup_tot(:) + (n_av_up(:)*weights(l,r))

          !exit loop if up and down densities have converged
          !IF((s.GT.2) .AND. (MAXVAL(change_up).LT.tol) .AND. (MAXVAL(change_down).LT.tol)) THEN
              !PRINT*, r
              !PRINT*, 'exiting'
              !EXIT
          !END IF

          k_laty = k_laty +delta_k
       END DO
       k_latx = k_latx + delta_k
    END DO
 nup_tot(:) = nup_tot(:)/(DBLE(nkx*nky))
 END DO

 !write converged densities to file
 !(loop scans site index)
 DO i = 1,N_s
    WRITE(32,*) i,  nup_tot(i), -ndown_tot(i)
 END DO


!-------------------------------------------|
!     BAND STRUCTURE & TOTAL ENERGY         |
!-------------------------------------------|

 k_latx = -pi
 H_o(:,:) = 0d0
 totE(:)=0d0
 DO r = 1,nkx+1
    k_laty = -pi
    DO l = 1,nky+1
       !calculate Hamiltonian parts for down Hamiltonian
       CALL hopping(t,H_t,N_s,a_x,a_y,k_latx,k_laty)
       CALL hubU(U,H_o,N_s,ndown_tot)

       !evaluate total up Hamiltonian
       DO i = 1,N_s
          DO j=1,N_s
             H_up(j,i) = H_t(j,i) + H_o(j,i)
          END DO
       END DO

       !find eigenvalues and eigenfunctions of up Hamiltonian
       CALL zheev('V', 'U', N_s, H_up, N_s, W_up, work, lwork, rwork, INFO)
       IF(INFO .NE. 0) THEN
          PRINT*, "Call to LAPACK routine zheev failed"
          STOP
       END IF

       !calculate Hamiltonian parts for up Hamiltonian
       CALL hopping(t,H_t,N_s,a_x,a_y,k_latx,k_laty)
       CALL hubU(U,H_o,N_s,nup_tot)

       !evaluate total up Hamiltonian
       DO i = 1,N_s
          DO j=1,N_s
             H_down(j,i) = H_t(j,i) + H_o(j,i)
          END DO
       END DO

       !find eigenvalues and eigenfunctions of up Hamiltonian
       CALL zheev('V', 'U', N_s, H_down, N_s, W_down, work, lwork, rwork, INFO)
       IF(INFO .NE. 0) THEN
           PRINT*, "Call to LAPACK routine zheev failed"
           STOP
       END IF

       DO i = 1,N_s
          WRITE(31,*) k_latx,k_laty, W_up(i), W_down(i)
       END DO

       WRITE(33,*) W_up  !for contour plot

       k_laty = k_laty + delta_k  !increment k_laty
    END DO
    WRITE(33,*) ' '  !for contour plot
    k_latx = k_latx + delta_k  !increment k_latx

    !normalise total energy
    DO i =1,N_s
       totE(i)=totE(i)/(dble(nkx*nky))
    END DO
 END DO

 !write total energy to file
 DO i=1,N_s
    WRITE(34,*) totE(i)
 END DO

 !close data files
 CLOSE(34)  !total energy file
 CLOSE(33)  !contour band structure file
 CLOSE(32)  !band structure file
 CLOSE(31)  !density of states file

 !deallocate arrays
 DEALLOCATE(H_up, H_down, H_t, H_o, U,weights)
 DEALLOCATE(n_av_up,n_av_down, change_up,change_down,ndown_tot,nup_tot,energy_onsite,totE)
 DEALLOCATE(W_up,W_down,work,rwork)


CONTAINS

!-----------------------------------------|
!   	       HOPPING TERM		  |
!-----------------------------------------|
SUBROUTINE hopping(t,H_t,N_s,a_x,a_y,k_latx,k_laty)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N_s
 COMPLEX(KIND=DP), DIMENSION(1:N_s,1:N_s), INTENT(OUT)  :: H_t
 REAL(KIND=DP), INTENT (IN) :: t
 REAL(KIND=DP), INTENT(IN) :: a_x,a_y,k_latx,k_laty
 INTEGER :: i
 INTEGER :: a_n,b_n,c_n,d_n     !neighbour indices
 INTEGER :: a_p,b_p,c_p,d_p     !phase indicators
 INTEGER :: ioerron, ierrop             !error labels
 COMPLEX(KIND=DP) :: imag               !imaginary number i=sqrt(-1)

 OPEN(unit=19, file='neighbours.dat', status='old',iostat=ioerron)   !open neighbour list
 OPEN(unit=20, file='phases.dat', status='old',iostat=ierrop)        !open phases list

 H_t(:,:) = (0d0,0d0)

 imag = (0d0,1d0)

 DO i = 1,N_s
    READ(19,*) a_n,b_n,c_n,d_n
    READ(20,*) a_p,b_p,c_p,d_p
    H_t(i,a_n) = H_t(i,a_n)+t*(a_p-1)
    H_t(i,b_n) = H_t(i,b_n)+t*(b_p-1)
    H_t(i,c_n) = H_t(i,c_n)+t*(c_p-1)
    H_t(i,d_n) = H_t(i,d_n)+t*(d_p-1)
    H_t(i,a_n) = H_t(i,a_n)-t*a_p*EXP(imag*k_latx*a_x)
    H_t(i,b_n) = H_t(i,b_n)-t*b_p*EXP(imag*k_laty*a_y)
    H_t(i,c_n) = H_t(i,c_n)-t*c_p*EXP(-imag*k_latx*a_x)
    H_t(i,d_n) = H_t(i,d_n)-t*d_p*EXP(-imag*k_laty*a_y)
 END DO

 CLOSE(19)    !close neighbour list
 CLOSE(20)    !close phases list

 RETURN
END SUBROUTINE


!-----------------------------------------|
!   	       HUBBARD TERM		  |
!-----------------------------------------|
SUBROUTINE hubU(U,H_o,N_s,n_av)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N_s
 COMPLEX(KIND=DP), DIMENSION(1:N_s,1:N_s), INTENT(OUT) :: H_o
 REAL(KIND=DP), INTENT(IN) :: U(N_s)
 REAL(KIND=DP), INTENT(IN) :: n_av(N_s)
 INTEGER :: j

 H_o(:,:) = (0d0,0d0)

 DO j = 1,N_s
    H_o(j,j) = U(j)*n_av(j)
 END DO

 RETURN
END SUBROUTINE

END PROGRAM HubbardModel
\end{lstlisting}


\section{Neighbour List Codes}
\label{neighAPP}
\subsection{One-Dimensional Neighbour List}
\begin{lstlisting}
 PROGRAM neigbour_1D

IMPLICIT NONE

!declare variables
INTEGER :: n_x              !no of states inside unit cell along x
INTEGER, PARAMETER :: N_n=2 !no of neighbours; two for 1D, four for 2D, six for 3D
INTEGER :: N_s              !total no of states
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neighbour !neighbour list array
INTEGER :: i             !do-loop counter
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: phase     !phase list array

OPEN(UNIT=19,FILE='neighbours.dat',STATUS='new')
OPEN(UNIT=20,FILE='phases.dat', STATUS='new')

!set dimensions
n_x = 2

!find no of states
N_s = n_x

!allocate dimensions to neighbour list
ALLOCATE(neighbour(1:N_n,1:N_s))
ALLOCATE(phase(1:N_n,1:N_s))
neighbour(:,:) = 0      !integers in neighbour list indicate site indices
phase(:,:) = 0          !integers in phase list indicate phase:
                        !0-no shift, 1-shift left, 2-shift right

!find neighbours for each atom
DO i = 1, N_s
   neighbour(1,i) = i -1
   neighbour(2,i) = i +1
END DO

neighbour(1,1) = N_s
phase(1,1) = 1
neighbour(2,N_s) = 1
phase(2,N_s) = 1

WRITE(*,*) 'The neighbour list is:'
DO i = 1,N_s
  WRITE(*,*) i, neighbour(:,i)
  WRITE(19,*) neighbour(:,i)
END DO

WRITE(*,*) 'The phase list is:'
DO i = 1,N_s
   WRITE(*,*) i, phase(:,i)
   WRITE(20,*) phase(:,i)
END DO

CLOSE(19)
CLOSE(20)

DEALLOCATE(neighbour,phase)
END PROGRAM
\end{lstlisting}

\subsection{Two-Dimensional Neighbour List}
\begin{lstlisting}
PROGRAM neigbour_2D

IMPLICIT NONE

!declare variables
INTEGER :: n_x,n_y            !no of states inside unit cell along x and y
INTEGER, PARAMETER :: N_n=4   !no of neighbours; two for 1D, four for 2D, eight for 3D
INTEGER :: N_s                !total no of states
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neighbour !neighbour list array
INTEGER :: i             !do-loop counter
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: phase       !phase list array

OPEN(UNIT=19, FILE='neighbours.dat', STATUS='new')
OPEN(UNIT=20, FILE='phases.dat', STATUS='new')

!set dimensions (set whichever dimension doesn't exist to zero)
n_x = 4
n_y = 4

!find no of states
N_s = n_x*n_y

!allocate dimensions to neighbour list
ALLOCATE(neighbour(1:N_n,1:N_s))
ALLOCATE(phase(1:N_n,1:N_s))
neighbour(:,:) = 0  !integers in neighbour list indicate site indices
phase(:,:) = 0      !integers in phase list indicate phase shifts:
                    !0-no shift, 1-shift left, 2-shift up, 3-shift right, 4-shift down

!find neighbours for each atom
DO i = 1, N_s
   neighbour(1,i) = i -1
   neighbour(2,i) = i -n_x
   neighbour(3,i) = i+1
   neighbour(4,i) = i+n_x
END DO

!-------------------- LEFT COLUMN -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found in the left column.'
DO i = 1, n_x*(n_y)-1+1, n_x
   neighbour(1,i) = i + n_x -1
   phase(1,i) = 1
   PRINT*, i
END DO

!-------------------- TOP ROW -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found in the top row.'
DO i = 1, n_x
   neighbour(2,i) = i + (n_x*(n_y-1))
   phase(2,i) = 1
   PRINT*, i
END DO

!-------------------- RIGHT COLUMN -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found in the right column.'
DO i = n_x, N_s, n_x
  neighbour(3,i) = i - n_x +1
  phase(3,i) = 1
  PRINT*, i
END DO

!-------------------- BOTTOM ROW -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found in the bottom row.'
DO i = n_x*(n_y-1)+1, N_s
  neighbour(4,i) = i - (n_x*(n_y-1))
  phase(4,i) = 1
  PRINT*, i
END DO

WRITE(*,*) 'The neighbour list is:'
DO i = 1, N_s
 WRITE(*,*) i, neighbour(:,i)
 WRITE(19,*) neighbour(:,i)
END DO
WRITE(*,*) 'The phase list is:'
DO i = 1,N_s
 WRITE(*,*) i, phase(:,i)
 WRITE(20,*) phase(:,i)
END DO

CLOSE(19)
CLOSE(20)

DEALLOCATE(neighbour,phase)
END PROGRAM
\end{lstlisting}

\subsection{Three-Dimensional Neighbour List}
\begin{lstlisting}
PROGRAM neigbour_3D

IMPLICIT NONE

!declare variables
INTEGER :: n_x, n_y, n_z      !no of states inside unit cell along x, y and z
INTEGER,PARAMETER :: N_n=6    !no of neighbours; two for 1D, four for 2D, six for 3D
INTEGER :: N_s                !total no of states
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neighbour !neighbour list array
INTEGER :: i,k                !do-loop counters
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: phase     !phase list array

!open output files
OPEN(UNIT=2, FILE="neighbours.dat")
OPEN(UNIT=3, FILE="phases.dat")

!set dimensions (set whichever dimension doesn't exist to zero)
n_x = 4
n_y = 4
n_z = 4

!find total no of states
N_s = n_x * n_y * n_z

!allocate dimensions to neighbour list
ALLOCATE(neighbour(1:N_n,1:N_s))
ALLOCATE(phase(1:N_n,1:N_s))
neighbour(:,:) = 0    !integers in neighbour list indicate atom indices
phase(:,:) = 0        !integers in phase list indicate phase: 0-inside unit cell
                      !1 - step right, 2-step up, 3-step left, 4-step down,
                      !5-step forward, 6-step backward


WRITE(*,*) 'This code labels the atoms in a 3D unit cell in increasing'
WRITE(*,*) 'integers starting along the x-axis, secondly along the y-axis '
WRITE(*,*) 'and finally along the z-axis in a right-handed coordinate system.'
WRITE(*,*) 'The labelling thus starts at the origin (considered close top'
WRITE(*,*) 'left-hand corner) with atom 1 and finishes at (n_x,n_y,n_z)'
WRITE(*,*) '(the far bottom right-hand corner) with atom N_s.'


!find neighbours for each atom
DO i = 1,N_s
   neighbour(1,i) = i-1
   neighbour(2,i) = i-n_x
   neighbour(3,i) = i +1
   neighbour(4,i) = i+n_x
   neighbour(5,i) = i + (n_x*n_y)
   neighbour(6,i) = i - (n_x*n_y)
END DO

!-------------------- LEFT FACE -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found on the left face.'

!left face
DO i = 1,N_s-n_x+1, n_x
   neighbour(1,i) = i+n_x-1
   phase(1,i) = 1
   WRITE(*,*) 'left face', i
END DO


!-------------------- TOP FACE -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found on the top face.'
i = 1
DO k = 1,n_x*n_z
   neighbour(2,i) = i + (n_x*(n_y-1))
   phase(2,i) = 1
   WRITE(*,*) 'top face', i
   IF(MOD(k,n_x)==0) THEN
      i = i+(n_x*(n_y-1))+1
   ELSE
      i = i +1
   END IF
END DO

!-------------------- RIGHT FACE -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found on the right face.'
DO i = n_x, N_s, n_x
   neighbour(3,i) = i-n_x+1
   phase(3,i) = 1
   WRITE(*,*) 'right face', i
END DO

!-------------------- BOTTOM FACE -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found on the bottom face.'
i=n_x*(n_y-1)+1
DO k = 1,n_x*n_z
   neighbour(4,i)=i-(n_x*(n_y-1))
   phase(4,i) = 1
   WRITE(*,*) 'bottom face', i
   IF(MOD(k,n_x)==0) THEN
      i = i+(n_x*(n_y-1))+1
   ELSE
      i = i+1
   END IF
END DO

!-------------------- BACK FACE -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found on the back face.'
DO i = N_s-(n_x*n_y)+1, N_s
   neighbour(5,i) = i-(n_x*n_y*(n_z-1))
   phase(5,i) = 1
   WRITE(*,*) 'back face', i
END DO

!-------------------- FRONT FACE -----------------------------
WRITE(*,*) '-------------------------------------------------'
WRITE(*,*) 'The following atoms can be found on the front face.'
DO i = 1,n_x*n_y
   neighbour(6,i) = i + (n_x*n_y*(n_z-1))
   phase(6,i) = 1
   WRITE(*,*) 'front face', i
END DO

WRITE(*,*) '-------------------------------------------------'

!-------------------- RESULTS -----------------------------
WRITE(*,*) 'The neighbour list is:'
DO i = 1, N_s
 WRITE(*,*) i, neighbour(:,i)
 WRITE(2,*) neighbour(:,i)
END DO

WRITE(*,*) 'The phase list is:'
DO i = 1,N_s
  WRITE(*,*) i, phase(:,i)
  WRITE(3,*) phase(:,i)
END DO


DEALLOCATE(neighbour,phase)
CLOSE(2)
CLOSE(3)
END PROGRAM
\end{lstlisting}

\section{Traditional Solid State Many-Body Theories}
\label{MBtheories}
\subsection{Hartree Mean-Field Theory}
The Hartree mean-field theory is a simple many-body theory for calculating electronic structure which works by converting the problem into an effective single-particle problem \cite{Hartree1,Hartree2,Hartree3}.  The Hartree approximation assumes that the Coulomb energy of the many-body Hamiltonian can be written as the interaction between a single electron and the electrostatic field generated by the overall electron density $n(\mathbf{s})$, i.e.
$$
 \frac{1}{2}\sum_{i\neq j} \frac{e^2}{\left|\mathbf{r}_i-\mathbf{r}_j\right|} \rightarrow \sum_i \int \frac{n(\mathbf{r}')}{\left|\mathbf{r}_i-\mathbf{r}'\right|} d^3 \mathbf{r}'
$$
where $n(\mathbf{r})$ is the density of all electrons.  The Hartree approach employs the Born-Oppenheimer approximation \cite{Hartree1,Hartree2,Hartree3}.  The resulting Schr\"odinger equation describes non-interacting electrons in an effective potential (which includes all the other terms of the many-body Hamiltonian, except the nuclear kinetic energy which is assumed to be approximately zero in the Born-Oppenheimer approximation).  The Hartree Schr\"odinger equation is
$$
 \hat{H}_i \chi_i(\mathbf{r}) = \left(-\frac{\hbar^2\nabla^2}{2m}+V_0(\mathbf{r}) + \int \frac{n(\mathbf{r}')}{\left|\mathbf{r}-\mathbf{r}'\right|} d^3 \mathbf{r}' \right) \chi_i(\mathbf{r}) = \varepsilon_i \chi_i(\mathbf{r})
$$
The electronic eigenstates $\chi_i(\mathbf{r})$ are filled up one by one as dictated by the Pauli principle.  Thus for an $N$-electron system, the electron density is
$$
 n(\mathbf{r}) = \sum_i^N\left|\chi_i(\mathbf{r})\right|^2
$$

The effective potential acting on the electrons is
$$
 V_{eff}(\mathbf{r})=V_0(\mathbf{r}) + \int \frac{n(\mathbf{r}')}{\left|\mathbf{r}-\mathbf{r}'\right|} d^3 \mathbf{r}'
$$
i.e.~it depends on the electron density $n(\mathbf{r}')$.  The Hamiltonian therefore has to be solved iteratively to give a self-consistent solution.  An initial guess of the density gives an effective potential.  This generates the eigenstates $\chi_i(\mathbf{r})$ which in turn gives us an improved value of the density.  The process is then iterated until the density converges (i.e.~when the output density equals the input density to within some pre-defined tolerance limit).  Variational calculus (with Lagrange multipliers to include the orthonormalisation constraint of the wavefunctions) yields the \textit{Hartree equations} \cite{Hartree1,Hartree2,Hartree3}
$$
 E_H = \sum_i \langle \chi_i | \hat{T} + \hat{V}_0 | \chi_i \rangle - \frac{1}{2} \sum_{ij}^N \langle \chi_i \chi_j | \hat{U} | \chi_i\chi_j\rangle
$$

\subsection{Hartree-Fock Theory}
The Hartree mean-field approach ignores many-body effects and introduces an spurious interaction of the electron with itself.  The Hartree-Fock theory is an improvement on the Hartree mean-field theory in which the antisymmetric nature of fermionic wavefunctions is taken into account \cite{PauliSymmetry,slater}.  An antisymmetric wavefunction which obeys the Pauli principle is obtained from the \textit{Slater determinant} \cite{slater}
$$
 \Phi_{mb} (\mathbf{r}) = \frac{1}{\sqrt{N!}}\left|\begin{array}{c c c}
                                                           \xi_1(\mathbf{r}_1) & ... & \xi_1(\mathbf{r}_N)\\
                                                           ... & ... & ... \\
                                                           \xi_N(\mathbf{r}_1) & ... & \xi_N(\mathbf{r}_N)\\
                                                          \end{array}\right|
$$
The orbitals constituting the Slater determinant can be determined from the single-particle Schr\"odinger equation
$$
 \left(-\frac{\hbar^2 \nabla^2}{2m} + V(\mathbf{r}\right) \xi_i(\mathbf{r}) = \varepsilon_i \xi_i(\mathbf{r})
$$


The other way of doing this is from the Euler-Lagrange equations
$$
 \frac{\delta F[n]}{\delta n(\mathbf{r})}|_{n=n_0} + V(\mathbf{r}) - \mu = 0
$$
Solving this using variational calculus and Lagrange multipliers give the Hartree-Fock equations
$$
 E_{HF} = \sum_i \langle \xi_i | \hat{T} + \hat{V}_0 | \xi_i \rangle - \frac{1}{2} \sum_{ij}^N \left[ \langle \xi_i \xi_j | \hat{U} | \xi_i\xi_j\rangle - \langle \xi_i \xi_j | \hat{U} | \xi_j\xi_i\rangle\right]
$$
The Hartree-Fock theory thus takes exchange into account but still misses higher-order many-body effects (\textit{correlation effects}).

\section{Exchange Functional Approximations}
\label{exchfapp}
\paragraph{Local Density Approximation}
\ \\
The local density approximation (LDA) is a simplification of the exchange-correlation functional.  In the LDA, the system is considered as a collection of small volume elements $\delta V$, in which the density is assumed to be constant \cite{NATODFT, LDA}.  If the centroid of the volume element is at position $\mathbf{r}$, then the density throughout the entire volume element is assumed to be $n(\mathbf{r})$.  One then further assumes that the system can locally be considered a homogeneous electron gas - the exchange-correlation energy of a volume element at $\mathbf{r}$ is approximated as the exchange-correlation energy of a homogeneous electron gas of density $n(\mathbf{r})$ \cite{NATODFT, LDA}.  The total exchange-correlation energy thus becomes a density-weighted integral summing over (infinitesimal) volume elements
$$
 E_{xc}^{LDA}[n] = \iiint n(\mathbf{r}) \epsilon_{xc}^{hom}(n(\mathbf{r})) d^3\mathbf{r}
$$


\paragraph{Local Spin Density Approximation}
\ \\
The local spin density approximation (LSDA) is the spin-polarised version of the LDA, i.e.
$$
 E_{xc}^{LSDA} = \iiint (n_{\uparrow}(\mathbf{r}) + n_{\downarrow}(\mathbf{r})) \epsilon_{xc}^{hom}(n(\mathbf{r})) d^3\mathbf{r}
$$

\paragraph{Gradient Approximations}
\ \\
The local density approximation can be improved by letting exchange-correlation energy depend not only on the local density but also on the local density gradients.  This is known as the gradient expansion approximation (GEA) \cite{NATODFT}
$$
 E_{xc}^{GEA}[n] = \iiint n(\mathbf{r}) \left[\epsilon_{xc}^{hom}(n(\mathbf{r})) + B(n(\mathbf{r}))\nabla n(\mathbf{r}) +   C(n(\mathbf{r}))\nabla^2 n(\mathbf{r}) + ...\right]
$$

This can be further improved by using more general functions (instead of using power-series like expansions) of the density and the density gradient.  This is known as the \textit{generalised gradient approximation} (GGA) \cite{GGA1,GGA2} and has the form
$$
 E_{xc}^{GGA}[n] = \iiint f\left(n(\mathbf{r}), \nabla n(\mathbf{r})\right) d^3 \mathbf{r}
$$
There are many different GGAs (i.e.~different forms of $f\left(n(\mathbf{r}), \nabla n(\mathbf{r})\right)$) developed for different purposes.  Examples include the Perdew-Burke-Ernzerhof (PBE) functional \cite{GGA1,GGA2} and the Perdew-Wang (PW91) functional \cite{PW91,PW1,PW2}.

\section{Pseudopotentials}
\label{USPapp}
\paragraph{Norm-Conserving Pseudopotentials} (NCPs) are pseudopotentials that meet the following condition (in addition to the transferability and softness requirements)
$$
 \int_{core} \left|\Psi(\mathbf{r})\right|^2 = \int_{core} \left|\Psi_{pseudo}(\mathbf{r})\right|^2
\label{conserve_norm}
$$
i.e.~they have to conserve the norm on the core-region.  This requirement means that norm-conserving pseudopotentials reproduce the all-electron wavefunction outside the cut-off radius.


\paragraph{Ultrasoft Pseudopotentials} (USPs) \cite{USP1,USP2} were developed for very large systems (large numbers of atoms), where the number of plane-waves is often the limiting factor.  In this case it is often useful to relax the norm-conserving condition (eq.~\ref{conserve_norm}) so as to produce as smooth pseudopotentials as possible.  Correction terms for the charge density have to be included if using ultrasoft pseudo-potentials.

\paragraph{Projector Augmented Waves}
The projector augmented waves (PAW) method \cite{PAW1,PAW2} is similar to the ultrasoft pseudopotential method but uses a direct transformation between the all-electron wavefunction and the pseudo-wavefunction.  The core wavefunctions calculated inside the cut-off radius are kept as exact all-electron wavefunctions.  Outside the cut-off radius smoother wavefunctions are used.

\section{VASP Parameters Used for Structure Relaxations}
The following parameters was settled on following parameter optimisation as described in section \ref{VASPopt}.
\begin{tabular}{| l | l | l | l |}
\hline
 \textbf{System} & ENCUT / eV & KPOINTS & SIGMA\\
\hline
FCC Fe (FM) &750 & 25 25 25 & 0.20\\
\hline
BCC Fe (FM) & 600 & 25 25 25 & 0.20\\
\hline
FePt (FM) & 550 & 24 24 18 & 0.15 \\
\hline
FePt(AFM) & 550 & 19 19 10 & 0.13 \\
\hline
FePt (PM) & 650 & 26 26 19 & 0.14 \\
\hline
CoPt(FM) & 650 & 26 26 19 & 0.20 \\
\hline
CoPt (AFM) & 650 & 19 19 10 & 0.12 \\
\hline
CoPt(PM) & 650 & 26 26 19 & 0.20\\
\hline
FeRh (FM) & 600 & 24 24 24 & 0.20 \\
\hline
FeRh (AFM) & 600 & 20 20 14 & 0.12 \\
\hline
FeRh(PM) & 600 & 24 24 24 & 0.20\\
\hline
\end{tabular}


\end{document}