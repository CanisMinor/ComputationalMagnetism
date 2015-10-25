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
