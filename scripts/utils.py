import numpy as np
import torch
import torch.nn as nn
import math
from torch import Tensor

from sklearn.linear_model import LinearRegression, Ridge
from scipy.special import gamma, eval_genlaguerre


class Morse3DExpNormalSmearing(torch.nn.Module):
    """
    Class works with 3D input data for
    """

    def __init__(self, cutoff_lower=0.0, cutoff_upper=5.0, num_rbf=50, trainable=True):
        """ """
        super(Morse3DExpNormalSmearing, self).__init__()
        self.rbf = ExpNormalSmearing(
            cutoff_lower=cutoff_lower,
            cutoff_upper=cutoff_upper,
            num_rbf=num_rbf,
            trainable=trainable,
        )
        self.linear = nn.Linear(num_rbf, 1, bias=False)

    def potential(self, x):
        """
        Calculate 1D Morse potential based on 3d coordinates

        Parameters:
        -----------
        x : array, 3xM
        Input coordinates for a sample. Each element contains Cartesian coordinates (x, y, z) for the particle


        """
        radial = torch.linalg.norm(x, dim=1)
        rbf_transformed = self.rbf(radial.view(-1, 1))
        potential = self.linear(rbf_transformed)
        return potential

    def forward(self, x):
        """
        Calculate potential and its Cartesian gradient.
        """
        potential = self.potential(x)
        du_dr = torch.autograd.grad(
            potential,
            x,
            grad_outputs=torch.ones_like(potential),
            retain_graph=True,
            create_graph=True,
        )
        return -du_dr[0]


class MorseOscillator:
    """
    Perform quantum oscillator calculations.
    All the units internally should be in atomic units
    """

    def __init__(self, D, a, m):
        """
        Add disociation energy and width parameter
        """
        self.D = D
        self.a = a
        self.m = m
        self.h_bar = 1.0
        self._lambda = np.sqrt(2 * self.m * self.D) / (self.a * self.h_bar)
        self.omega_0 = np.sqrt(2 * self.D * self.a**2 / self.m)
        # maximum number of n
        self.n_max = math.floor(self._lambda - 0.5)

    def ksi(self, x):
        return 2 * self._lambda * np.exp(-self.a * x)

    def potential(self, x):
        """
        Calculate potential energy at a point x
        """
        return self.D * (1 - np.exp(-self.a * x)) ** 2

    def energy_level(self, n):
        """
        Calculate Hamiltonian eigenvalue for a given quantum
        number n <= n_max
        """
        assert n <= self.n_max, f"Only bound states with n < {self.n_max} are allowed"
        e = (
            self.omega_0
            * self.h_bar
            * (n + 1 / 2 - 1 / (2 * self._lambda) * (n + 1 / 2) ** 2)
        )
        return e

    def Hamiltonian_eigenvalues(self):
        """
        Calculate energy spectrum
        """
        eigenvalues = []
        for n in range(self.n_max + 1):
            eigenvalues.append(self.energy_level(n))
        return np.array(eigenvalues)

    def Hamiltonian_eigenfunction_direct(self, n, x):
        """
        Calculate eigenfunction as defined by eq. 38,
        without recursion. CAUSION: potential precision loss.
        """
        s = 2 * self._lambda - 2 * n - 1
        # normalization constant
        N = np.sqrt(s * gamma(n + 1) / gamma(2 * self._lambda - n))
        ksi = self.ksi(x)
        L = eval_genlaguerre(n, s, ksi)  # genlaguerre(n, s)
        phi = N * np.power(ksi, s / 2) * np.exp(-ksi / 2) * L
        return phi

    def Hamiltonian_eigenfunction(self, n, x):
        """
        Calculate using recursive relationship, given in equation 49
        """
        s = 2 * self._lambda - 2 * n - 1
        ksi = self.ksi(x)
        # Redefine N in terms of and s
        N_s_0 = np.sqrt(s / gamma(s + 1))
        psi_0 = N_s_0 * np.exp(-ksi / 2) * np.power(ksi, s / 2)

        # Loop over nx from 0 to n to get the eigenfuncion without computing Laguerre
        # polynomial
        if n == 0:
            return psi_0
        else:
            psi_minus_2 = 0
            psi_minus_1 = psi_0
            for nx in range(1, n + 1):
                psi_n = np.sqrt(1 / ((nx + s) * nx)) * (
                    (2 * nx + s - 1 - ksi) * psi_minus_1
                    - np.sqrt((nx - 1) * (nx + s - 1)) * psi_minus_2
                )
                psi_minus_2 = psi_minus_1
                psi_minus_1 = psi_n
        return psi_n

    def Hamiltonian_eigenfunctions_spectrum(self, x):
        eigenfunctions = []
        for n in range(self.n_max + 1):
            eigenfunctions.append(self.Hamiltonian_eigenfunction_direct(n, x))
        return np.array(eigenfunctions)

    def Hamiltonian_eigenfunctions_spectrum2(self, x):
        """
        Calculate the entire eigenfunction spectrum, using a recursive relationship without
        Laguerre polynomials at all.
        """
        eigenfunctions = []
        for n in range(self.n_max + 1):
            eigenfunctions.append(self.Hamiltonian_eigenfunction(n, x))
        return np.array(eigenfunctions)

    def partition_function(self, beta):
        """
        Quantum partition function of Morse potential.
        """
        E_spectrum = self.Hamiltonian_eigenvalues()
        exponent = np.exp(-1 * beta * E_spectrum)
        q = np.sum(exponent)
        return q

    def pmf(self, beta, x):
        """
        Get potential of mean force. Energies a
        First, need to get probability  of having coordinate x for state i, then sum over all the
        states and take a log
        """
        eigenfunctions = self.Hamiltonian_eigenfunctions_spectrum2(x)

        psi_squared = np.square(eigenfunctions)
        energies = self.Hamiltonian_eigenvalues()
        energy_exponent = np.exp(-1.0 * beta * energies)
        probability = np.zeros(x.shape)
        for i in range(self.n_max + 1):
            probability += np.multiply(energy_exponent[i], psi_squared[i])
        q = np.sum(energy_exponent)
        probability = probability / q
        pmf = -np.log(probability) / beta
        return pmf


class LinearRegressionForceMatching:
    """
    Class formulates force-matching as a linear regression
    problem in a basis of radial basis functions
    """

    def __init__(
        self,
        cutoff_lower=0.0,
        cutoff_upper=5.0,
        num_rbf=50,
        device="cpu",
        model_type="linear_regression",
        alpha=0.0,
    ):
        self.rbf = ExpNormalSmearing(
            cutoff_lower=cutoff_lower,
            cutoff_upper=cutoff_upper,
            num_rbf=num_rbf,
            trainable=False,
        )
        self.rbf.to(device)
        self.params = np.ones(num_rbf)
        self.optimized = False
        self.n_rbf = num_rbf
        self.potential_basis_computed = False
        self.force_basis_computed = False
        self.model_type = model_type
        self.alpha = alpha

    def to(self, device):
        self.rbf.to(device)

    def potential_basis(self, x):
        """
        Calculate 1D Morse potential based on 3d coordinates

        Parameters:
        -----------
        x : array, 3xM
        Input coordinates for a sample. Each element contains Cartesian coordinates (x, y, z) for the particle


        """
        radial = torch.linalg.norm(x, dim=1)
        rbf_transformed = self.rbf(radial.view(-1, 1))
        return rbf_transformed.view(-1, self.n_rbf)

    def force_basis(self, x):
        """
        Calculate potential and its Cartesian gradient.
        """
        potential_basis = self.potential_basis(x)

        # jacobian = torch.autograd.functional.jacobian(self.potential_basis, x)
        # jacobian.detach()
        # print("jacobian calculated")
        # du_dr = torch.diagonal(jacobian, dim1=0, dim2=2)
        # return -1*du_dr
        du_dr_list = []
        for i in range(self.n_rbf):
            rbf = potential_basis[:, i]
            du_dr = torch.autograd.grad(
                rbf.sum(), x, grad_outputs=None, retain_graph=True
            )[0]
            du_dr_list.append(du_dr.T)
            total_du_dr = torch.stack(du_dr_list, axis=0).detach()
        return -1 * total_du_dr

    def get_training_forces_basis(self, x):
        """
        Prepare forces in the format needed for training.
        Output is a numpy array with
        """

        force_basis = self.force_basis(x)
        initial_shape = force_basis.shape
        force_basis_reshaped = force_basis.reshape(initial_shape[0], -1)
        training_forces_basis = force_basis_reshaped.detach().to("cpu").numpy().T
        return training_forces_basis

    def precompute_potential_basis(self, x):
        self.potential_basis = self.potential_basis(x)
        self.potential_basis_computed = True

    def precompute_force_basis(self, x):
        self.force_basis = get_training_forces_basis(self, x)
        self.force_basis_computed = True

    def format_reference_forces(self, reference_forces):
        return reference_forces.to("cpu").detach().numpy().reshape(-1, order="F")

    def fit(self, coordinates, forces, fit_intercept=False, **kwargs):
        """
        Method takes raw coordinates and reference forces and fit
        linear model.
        """
        training_forces_basis = self.get_training_forces_basis(coordinates)
        reference_forces = self.format_reference_forces(forces)
        if self.model_type == "linear_regression":
            model = LinearRegression(fit_intercept=fit_intercept, **kwargs).fit(
                training_forces_basis, reference_forces
            )
        elif self.model_type == "ridge":
            model = Ridge(alpha=self.alpha).fit(training_forces_basis, reference_forces)
        else:
            raise NotImplementedError
        self.model = model
        self.optimized = True
        self.score = model.score(training_forces_basis, reference_forces)
        self.params = np.array(model.coef_)
        return model

    def predict(self, coordinates, **kwargs):
        test_forces_basis = self.get_training_forces_basis(coordinates)
        predicted_forces = self.model.predict(test_forces_basis)
        return predicted_forces

    # Methods below calculate potential and forces based on coordinates
    def calculate_potential(self, x):
        potential_basis = self.potential_basis(x).detach().numpy()
        return np.dot(potential_basis, self.params)

    def calculate_forces(self, x):
        force_basis = self.force_basis(x)
        forces = np.multiply(force_basis, np.expand_dims(self.params, axis=(1, 2))).sum(
            axis=0
        )
        return forces


class ExpNormalSmearing(nn.Module):
    """
    This class is adapted from  torchmdnet.models.utils.ExpNormalSmearing
    (https://github.com/torchmd/torchmd-net) and is distrubuted under
    MIT license:

    Copyright (c) 2021-2023 Universitat Pompeu Fabra,  https://www.compscience.org

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """

    def __init__(
        self,
        cutoff_lower=0.0,
        cutoff_upper=5.0,
        num_rbf=50,
        trainable=True,
        dtype=torch.float32,
    ):
        super(ExpNormalSmearing, self).__init__()
        self.cutoff_lower = cutoff_lower
        self.cutoff_upper = cutoff_upper
        self.num_rbf = num_rbf
        self.trainable = trainable
        self.dtype = dtype
        self.cutoff_fn = CosineCutoff(0, cutoff_upper)
        self.alpha = 5.0 / (cutoff_upper - cutoff_lower)

        means, betas = self._initial_params()
        if trainable:
            self.register_parameter("means", nn.Parameter(means))
            self.register_parameter("betas", nn.Parameter(betas))
        else:
            self.register_buffer("means", means)
            self.register_buffer("betas", betas)

    def _initial_params(self):
        # initialize means and betas according to the default values in PhysNet
        # https://pubs.acs.org/doi/10.1021/acs.jctc.9b00181
        start_value = torch.exp(
            torch.scalar_tensor(
                -self.cutoff_upper + self.cutoff_lower, dtype=self.dtype
            )
        )
        means = torch.linspace(start_value, 1, self.num_rbf, dtype=self.dtype)
        betas = torch.tensor(
            [(2 / self.num_rbf * (1 - start_value)) ** -2] * self.num_rbf,
            dtype=self.dtype,
        )
        return means, betas

    def reset_parameters(self):
        means, betas = self._initial_params()
        self.means.data.copy_(means)
        self.betas.data.copy_(betas)

    def forward(self, dist):
        dist = dist.unsqueeze(-1)
        return self.cutoff_fn(dist) * torch.exp(
            -self.betas
            * (torch.exp(self.alpha * (-dist + self.cutoff_lower)) - self.means) ** 2
        )


class CosineCutoff(nn.Module):
    """
    This class is adapted from  torchmdnet.models.utils.ExpNormalSmearing
    (https://github.com/torchmd/torchmd-net) and is distrubuted under
    MIT license:

    Copyright (c) 2021-2023 Universitat Pompeu Fabra,  https://www.compscience.org

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """

    def __init__(self, cutoff_lower=0.0, cutoff_upper=5.0):
        super(CosineCutoff, self).__init__()
        self.cutoff_lower = cutoff_lower
        self.cutoff_upper = cutoff_upper

    def forward(self, distances: Tensor) -> Tensor:
        if self.cutoff_lower > 0:
            cutoffs = 0.5 * (
                torch.cos(
                    math.pi
                    * (
                        2
                        * (distances - self.cutoff_lower)
                        / (self.cutoff_upper - self.cutoff_lower)
                        + 1.0
                    )
                )
                + 1.0
            )
            # remove contributions below the cutoff radius
            cutoffs = cutoffs * (distances < self.cutoff_upper)
            cutoffs = cutoffs * (distances > self.cutoff_lower)
            return cutoffs
        else:
            cutoffs = 0.5 * (torch.cos(distances * math.pi / self.cutoff_upper) + 1.0)
            # remove contributions beyond the cutoff radius
            cutoffs = cutoffs * (distances < self.cutoff_upper)
            return cutoffs


from scipy.linalg import norm


def compute_spring_energies(pos1, pos2, masses, P, hbar, beta, scaling=1.0e0):
    """
    Takes two arrays of positions, pos_1, pos_2 and computes energies
    associated with this spring

    Parameters:
    pos1, pos2: numpy.ndarray

        Arrays with positions with shape (N_frames, N_atoms 3)

    masses: np.array, shape (N_atoms)
        masses of particles in the system,
        in the same order as in pos arrays

    P: int
        Number of beads in the path integral

    hbar: float
        Planck's constant,  units consistent with eV, Angstrom, atomic mass

    beta: float
        1/(kB*T), units consistent with eV, Angstrom, atomic mass

    """
    prefactor = masses * P / (beta * hbar) ** 2 / 2
    spring_energy = np.square(norm(pos1 - pos2, axis=2))
    spring_energy *= scaling * prefactor
