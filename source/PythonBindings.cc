#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <alpaca/AngularCorrelation.hh>
// #include <alpaca/EulerAngleRotation.hh>
// #include <alpaca/CascadeSampler.hh>
#include <alpaca/State.hh>
#include <alpaca/Transition.hh>

namespace py = pybind11;
using namespace alpaca;

PYBIND11_MODULE(alpaca, m) {
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
      .def("__repr__", &Transition::str_rep);

  py::enum_<Parity>(m, "Parity")
      .value("negative", Parity::negative)
      .value("positive", Parity::positive)
      .value("unknown", Parity::unknown);
  py::class_<State>(m, "State")
      .def(py::init<int>(), py::arg("two_J"))
      .def(py::init<int, Parity>(), py::arg("two_J"), py::arg("parity"))
      .def(py::init<int, Parity, double>(), py::arg("two_J"), py::arg("parity"),
           py::arg("energy"), R"doc(
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
      .def(py::init<int, double>(), py::arg("two_J"), py::arg("energy"))
      .def_readwrite(
          "two_J", &State::two_J,
          "Two times the angular momentum quantum number in units of the "
          "reduced Planck constant")
      .def_readwrite("parity", &State::parity, "Parity quantum number")
      .def_readwrite(
          "excitation_energy", &State::excitation_energy,
          "Excitation energy of the state with respect to the ground state "
          "in MeV")
      .def("__repr__", &State::str_rep);

  py::class_<AngularCorrelation>(m, "AngularCorrelation")
      .def(py::init<State, std::vector<std::pair<Transition, State>>>())
      .def(py::init<State, std::vector<State>>())
      .def("__call__",
           py::vectorize(py::overload_cast<const double, const double>(
               &AngularCorrelation::operator(), py::const_)));
}
