#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <alpaca/AngularCorrelation.hh>
#include <alpaca/CascadeSampler.hh>
#include <alpaca/DeterministicReferenceFrameSampler.hh>
#include <alpaca/EulerAngleRotation.hh>
#include <alpaca/ReferenceFrameSampler.hh>
#include <alpaca/SpotlightSampler.hh>
#include <alpaca/State.hh>
#include <alpaca/Transition.hh>

namespace py = pybind11;
using namespace alpaca;

PYBIND11_MODULE(_alpaca, m) {
  m.attr("__name__") = "alpaca";
  py::enum_<EMCharacter>(m, "EMCharacter")
      .value("electric", EMCharacter::electric)
      .value("magnetic", EMCharacter::magnetic)
      .value("unknown", EMCharacter::unknown);
  py::class_<Transition>(m, "Transition")
      .def(py::init<int, int, double>(), py::arg("two_L"), py::arg("two_Lp"),
           py::arg("delta") = 0.)
      .def(py::init<EMCharacter, int, EMCharacter, int, double>(),
           py::arg("em"), py::arg("two_L"), py::arg("emp"), py::arg("two_Lp"),
           py::arg("delta") = 0., R"doc(
            Constructor

            Parameters
            ----------
            em: alpaca.EMCharacter
                Primary EM character.
            t_L: int
                Two times the primary multipolarity.
            emp: alpaca.EMCharacter
                Secondary EM character.
            t_Lp: int
                Two times the secondary multipolarity. Must be different from t_L.
            delta: float
                Multipole mixing ratio.
         )doc")
      .def_readwrite("em_char", &Transition::em_char, "Primary EM character")
      .def_readwrite("two_L", &Transition::two_L,
                     "Two times the primary multipolarity")
      .def_readwrite("em_charp", &Transition::em_charp,
                     "Secondary EM character")
      .def_readwrite("two_Lp", &Transition::two_Lp,
                     "Two times the secondary multipolarity")
      .def_readwrite("delta", &Transition::delta, "Multipole mixing ratio")
      .def("__repr__", &Transition::str_rep)
      .def_static("Dipole", &Transition::Dipole, py::arg("delta") = 0.)
      .def_static("E1", &Transition::E1, py::arg("delta") = 0.)
      .def_static("M1", &Transition::M1, py::arg("delta") = 0.)
      .def_static("Quadrupole", &Transition::Quadrupole, py::arg("delta") = 0.)
      .def_static("E2", &Transition::E2, py::arg("delta") = 0.)
      .def_static("M2", &Transition::M2, py::arg("delta") = 0.);

  py::enum_<Parity>(m, "Parity")
      .value("negative", Parity::negative)
      .value("positive", Parity::positive)
      .value("unknown", Parity::unknown);
  py::class_<State>(m, "State")
      .def(py::init<int>(), py::arg("two_J").noconvert())
      .def(py::init<int, Parity>(), py::arg("two_J").noconvert(),
           py::arg("parity"))
      .def(py::init<int, Parity, double>(), py::arg("two_J").noconvert(),
           py::arg("parity"), py::arg("energy"), R"doc(
            Constructor

            Parameters
            ----------
            two_J: int
                Two times the angular momentum quantum number in units of the reduced Planck constant.
            parity: alpaca.Parity
                Parity quantum number, which may be 1 (positive parity), -1 (negative), or 0 (parity unknown).
            energy: float
                Excitation energy of the state with respect to the ground state in MeV
        )doc")
      .def(py::init<int, double>(), py::arg("two_J").noconvert(),
           py::arg("energy"))
      .def_readwrite(
          "two_J", &State::two_J,
          "Two times the angular momentum quantum number in units of the "
          "reduced Planck constant")
      .def_readwrite("parity", &State::parity, "Parity quantum number")
      .def_readwrite(
          "excitation_energy", &State::excitation_energy,
          "Excitation energy of the state with respect to the ground state "
          "in MeV")
      .def("__repr__", &State::str_rep)
      .def_static("Even", &State::Even, py::arg("J"),
                  py::arg("parity") = Parity::unknown)
      .def_static("EvenPlus", &State::EvenPlus, py::arg("J"))
      .def_static("EvenMinus", &State::EvenMinus, py::arg("J"))
      .def_static("Odd", &State::Odd, py::arg("J"),
                  py::arg("parity") = Parity::unknown)
      .def_static("OddPlus", &State::OddPlus, py::arg("J"))
      .def_static("OddMinus", &State::OddMinus, py::arg("J"));

  py::class_<AngularCorrelation>(m, "AngularCorrelation", R"doc(
                Class for a gamma-gamma correlation.

                Calculates the angular correlation \f$W_{\gamma \gamma} \left( \theta, \varphi \right)\f$
                between the first and the last photon in a sequence of \f$n-1\f$ (\f$n > 2\f$)
                electromagnetic (EM) transitions between \f$n\f$ states of a quantized system
                (a 'cascade' of EM transitions with \f$n-1\f$ steps).
                The angles \f$\theta \in \left[ 0, \pi \right]\f$ and
                \f$\varphi \in \left[ 0, 2 \pi \right]\f$ are the polar and azimuthal angles in a spherical
                coordinate system, respectively.
                States of the system are identified by their total angular momentum quantum numbers ('spins')
                \f$J\f$ (\f$2J \in \mathcal{N} \f$, \f$J \geq 0\f$) and their parities
                \f$\pi \in \lbrace -1, 1 \rbrace \f$.
                They are labeled by indices \f$0 < i \geq n\f$, where \f$i = 1\f$ denotes the first, and
                \f$i = n\f$ the last state of a cascade.
                EM transitions are identified by their multipolarity \f$L\f$ (\f$L \in \mathcal{N} \f$, \f$L \geq 0\f$) and their EM character
                \f$\lambda \in \lbrace \mathrm{E}, \mathrm{M} \rbrace\f$ which can be either electric (E) or
                magnetic (M).
                They are labeled by indices \f$1 \geq i < n\f$, where the \f$i\f$-th transition is assumed to be
                the transition that connects the states labeled \f$i\f$ and \f$i+1\f$.
                A transition between two states may be a mixture of up to two multipolarities \f$L_i\f$ and
                \f$L_i^\prime\f$, whose relative contributions are quantified by the corresponding
                multipole mixing ratio \f$\delta_i\f$
                (see below about the convention for \f$\delta\f$ which is used in the present implementation).

                The minimum possible cascade with \f$n = 3\f$ states is denoted symbolically as:

                \f[
                    J_1^{\pi_1}
                    \left( \begin{array}{ccc} L_1 \\ L_1^\prime \end{array} \right)
                    J_2^{\pi_2}
                    \left( \begin{array}{ccc} L_2 \\ L_2^\prime \end{array} \right)
                    J_3^{\pi_3}.
                \f]

                A cascade between \f$n > 3\f$ states is denoted as:

                \f[
                    J_1^{\pi_1}
                    \left( \begin{array}{ccc} L_1 \\ L_1^\prime \end{array} \right)
                    J_2^{\pi_2}
                    \left( \begin{array}{ccc} L_2 \\ L_2^\prime \end{array} \right)
                    ...
                    J_{n-1}^{\pi_{n-1}}
                    \left( \begin{array}{ccc} L_{n-1} \\ L_{n-1}^\prime \end{array} \right)
                    J_n^{\pi_n}.
                \f]

                The transitions with \f$ 1 < i < n-1\f$ are assumed to be unobserved.
                The photon which corresponds to the first (see below about the notion of 'first' and 'second')
                transition is assumed to be observed along the positive z direction (\f$\theta = 0\f$).
                It may also be interpreted as a beam of photons which travels in that direction and causes
                the excitation of the state \f$J_2^{\pi_2}\f$ from the initial state \f$J_1^{\pi_1}\f$.
                Apart from a correlation between the directions of motion of the two photons
                [direction-direction (dir-dir) correlation], the present implementation can take into
                account that the polarization of the first photon is observed in addition [a
                polarization-direction (pol-dir) correlation].
                Here, it is assumed that the polarization is along the x axis
                (\f$\theta = \pi/2\f$, \f$\varphi = 0\f$ denotes the positive x axis).
                Only the observation of polarization information, which introduces a dependence of the
                angular correlation on the azimuthal angle \f$\varphi\f$, allows for a distinction between
                different EM characters, which are also related to the parities of the corresponding states.
                For this reason, the code will assume a pol-dir correlation if parities and EM
                characters associated with the first (the 'first' transition is well-defined in the user interface of the
                code) transition are given, and a dir-dir correlation otherwise.

                The formalism of angular correlations which was used in the present implementation is mainly
                based on a review article by Fagg and Hanna
                \cite FaggHanna1959 and on a book chapter by Biedenharn \cite AjzenbergSelove1960.
                Consequently, Biedenharn's convention for the multipole mixing ratio is used.
                For a comparison to the other popular conventions of Rose and Brink and Krane, Steffen, and
                Wheeler, see Refs. \cite RoseBrink1967 and \cite KraneSteffenWheeler1973, respectively.
                When a two-step cascade is considered in which the first and the last state are identical,
                Biedenharn's convention has the advantage that \f$\delta_1 = \delta_2\f$.

                In the angular correlation formalism, the expansion coefficients of \f$W_{\gamma \gamma}\f$
                in terms of Legendre polynomials are separable into contributions by the different transitions.
                Therefore, a 'first' and 'last' transition of the cascade actually need not be identified.
                In other words, it does not matter whether the photon observed in z direction with an
                polarization in x direction is the first or the last one of the cascade.
                However, the unobserved transitions must take place in between the two, otherwise they would
                not have to be considered at all.

                In order to describe an observation of the 'first' photon in an arbitrary direction with an
                arbitrary orientation of the polarization in the plane perpendicular to the direction of
                propagation, the code allows to rotate the angular correlation by three Euler angles
                \f$\Phi\f$, \f$\Theta\f$, and \f$\Psi\f$.
                They are defined to be the rotation angles around the z-, x', and z' axes, respectively, as
                described, for example, in Ref. \cite Weisstein2020.

                Like many other quantum mechanical computer codes, this one also uses \f$2 J\f$ instead of
                \f$J\f$ and \f$2L\f$ instead of \f$L\f$ to be able to represent both integer and half-integer
                angular momentum quantum numbers as integers internally.

                The angular correlation is normalized to \f$4\pi\f$, i.e.:

                \f[
                    \int_0^{2\pi} \int_0^\pi W_{\gamma \gamma} \left( \theta, \varphi \right) \sin \left( \theta \right) \mathrm{d} \theta \mathrm{d} \varphi = 4 \pi
                \f]

                An EM cascade is specified entirely by the arguments of the constructor of AngularCorrelation.
                The first argument is the initial state of the cascade, which is sometimes denoted as the
                'oriented' state in the literature, because its decay/excitation defines the coordinate
                system.
                All other cascade steps are given as a list of pairs of transitions and the states
                which they populate.
                The transition between the initial state and the second state, and the one between the
                \f$n-1\f$-th and the \f$n\f$-th state, are assumed to be the two observed transitions.
                All other transitions are treated as unobserved.

                If the parities of the initial state and the second state, and both EM characters of the
                first transition are given, the code will assume a pol-dir correlation.
                If none of these data is given, the code will assume a dir-dir correlation.
            )doc")
      .def(py::init<State, std::vector<std::pair<Transition, State>>>())
      .def(py::init<State, std::vector<State>>())
      .def("__call__",
           py::vectorize(py::overload_cast<const double, const double>(
               &AngularCorrelation::operator(), py::const_)),
           py::arg("theta"), py::arg("phi"),
           R"doc(
                Evaluate the angular correlation

                This function accepts more parameters than the two angles in spherical coordinates:

                Three Euler angles can be provided as additional parameters to
                rotate the direction of propagation and the polarization axis (if defined) of the first photon.
                The 'zxz' convention or 'x' convention is used for the order of the rotations
                \cite Weisstein2020.
                As implied by the notation ('x' instead of 'x prime'), the angles \f$\theta\f$ and
                \f$\varphi\f$ are still defined in the original coordinate system, i.e. \f$\theta = 0\f$
                is still the z axis.

                Parameters
                ----------
                theta: float or ndarray
                    Polar angle in spherical coordinates in radians (\f$\theta \in \left[ 0, \pi \right]\f$). If ndarray, must have the same shape as phi.
                phi: float or ndarray
                    Azimuthal angle in spherical coordinates in radians (\f$\varphi \in \left[ 0, 2 \pi \right]\f$). If ndarray, must have the same shape as theta.

                Returns
                -------
                float or ndarray
                    \f$W_{\gamma \gamma} \left( \theta, \varphi \right)\f$, value of the angular
                    correlation. If the values for theta and phi were scalars, a scalar will be returned.
                    If at least one or both of theta and phi was a numpy array of shape (M, N, ...), a
                    numpy array of shape (M, N, ...) will be returned.
           )doc")
      .def_property_readonly("cascade_steps",
                             &AngularCorrelation::get_cascade_steps)
      .def_property_readonly("initial_state",
                             &AngularCorrelation::get_initial_state);

  py::class_<ReferenceFrameSampler,
             std::shared_ptr<ReferenceFrameSampler> /* holder type */>(
      m, "ReferenceFrameSampler")
      .def("__call__", &ReferenceFrameSampler::operator());

  py::class_<AngCorrRejectionSampler,
             std::shared_ptr<AngCorrRejectionSampler> /* holder type */,
             ReferenceFrameSampler /* base class */>(m,
                                                     "AngCorrRejectionSampler")
      .def(py::init<AngularCorrelation &, unsigned long, unsigned int>(),
           py::arg("w"), py::arg("seed"), py::arg("max_tri") = 1000)
      .def("__call__", &AngCorrRejectionSampler::operator());

  py::class_<SpotlightSampler,
             std::shared_ptr<SpotlightSampler> /* holder type */,
             ReferenceFrameSampler /* base class */>(m, "SpotlightSampler")
      .def(py::init<CoordDir, double, unsigned long>(), py::arg("theta_phi"),
           py::arg("opening_angle"), py::arg("seed"))
      .def(py::init<CoordDir, unsigned long>(), py::arg("theta_phi"),
           py::arg("seed"))
      .def(py::init<CoordDir, double, double, unsigned long>(),
           py::arg("theta_phi"), py::arg("distance"), py::arg("radius"),
           py::arg("seed"))
      .def("__call__", &SpotlightSampler::operator());

  py::class_<
      DeterministicReferenceFrameSampler,
      std::shared_ptr<DeterministicReferenceFrameSampler> /* holder type */,
      ReferenceFrameSampler /* base class */>(
      m, "DeterministicReferenceFrameSampler")
      .def(py::init<EulerAngles>(), py::arg("phi_theta_psi"))
      .def("__call__", &DeterministicReferenceFrameSampler::operator());

  py::class_<SphereRejectionSampler,
             std::shared_ptr<SphereRejectionSampler> /* holder type */,
             ReferenceFrameSampler /* base class */>(m,
                                                     "SphereRejectionSampler")
      .def(py::init<std::function<double(double, double)>, double,
                    unsigned long, unsigned int>(),
           py::arg("dis"), py::arg("dis_max"), py::arg("seed"),
           py::arg("max_tri") = 1000)
      .def("__call__", &SphereRejectionSampler::operator());

  py::class_<CascadeSampler>(m, "CascadeSampler")
      .def(py::init<std::vector<std::shared_ptr<ReferenceFrameSampler>> &>())
      .def("__call__", &CascadeSampler::operator(), R"doc(
            Sample random gamma-ray directions from the cascade.

            Parameters
            ----------
            n: int
                Number of samples to return

            Returns
            -------
                Array of reference frames in Euler angles. If `n` is given, an array of
                arrays is returned instead (multiple samples). The first reference frame
                describes the direction of emission of the first (depends on the setting
                of return_first_direction) gamma ray in the cascade, the second pair
                describes the second gamma ray, and so on.
           )doc")
      .def("__call__", [](CascadeSampler &c, py::ssize_t len) {
        auto dim = static_cast<py::ssize_t>(c.size());
        auto res = py::array_t<double>(std::vector<ptrdiff_t>{
            static_cast<ptrdiff_t>(len), static_cast<ptrdiff_t>(dim), 3});

        auto data = res.mutable_unchecked<3>();

        for (py::ssize_t i = 0; i < len; ++i) {
          auto sample = c();
          for (py::ssize_t j = 0; j < dim; ++j) {
            data(i, j, 0) = sample[static_cast<size_t>(j)][0];
            data(i, j, 1) = sample[static_cast<size_t>(j)][1];
            data(i, j, 2) = sample[static_cast<size_t>(j)][2];
          }
        }
        return res;
      });
}
