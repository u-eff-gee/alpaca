#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

// #include <alpaca/AngularCorrelation.hh>
// #include <alpaca/CascadeSampler.hh>
#include <alpaca/State.hh>
#include <alpaca/Transition.hh>

namespace nb = nanobind;
using namespace alpaca;

NB_MODULE(alpaca, m) {
  nb::enum_<EMCharacter>(m, "EMCharacter")
    .value("electric", EMCharacter::electric)
    .value("magnetic", EMCharacter::magnetic)
    .value("unknown", EMCharacter::unknown);
  nb::class_<Transition>(m, "Transition")
    .def(nb::init<int, int, double>(), nb::arg("two_L"), nb::arg("two_Lp"),
         nb::arg("delta") = 0.)
    .def(nb::init<EMCharacter, int, EMCharacter, int, double>(),
         nb::arg("em"), nb::arg("two_L"), nb::arg("emp"),
         nb::arg("two_Lp"), nb::arg("delta") = 0., R"doc(
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
    .def_rw("em_char", &Transition::em_char, "Primary EM character")
    .def_rw("two_L", &Transition::two_L, "Two times the primary multipolarity")
    .def_rw("em_charp", &Transition::em_charp, "Secondary EM character")
    .def_rw("two_Lp", &Transition::two_Lp, "Two times the secondary multipolarity")
    .def_rw("delta", &Transition::delta, "Multipole mixing ratio")
    .def("__repr__", &Transition::str_rep);

  nb::enum_<Parity>(m, "Parity")
    .value("negative", Parity::negative)
    .value("positive", Parity::positive)
    .value("parity_unknown", Parity::unknown);
  nb::class_<State>(m, "State")
    .def(nb::init<int>(), nb::arg("two_J"))
    .def(nb::init<int, Parity>(), nb::arg("two_J"), nb::arg("parity"))
    .def(nb::init<int, Parity, double>(), nb::arg("two_J"), nb::arg("parity"),
         nb::arg("energy"), R"doc(
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
    .def(nb::init<int, double>(), nb::arg("two_J"), nb::arg("energy"))
    .def_rw("two_J", &State::two_J,
            "Two times the angular momentum quantum number in units of the "
            "reduced Planck constant")
    .def_rw("parity", &State::parity, "Parity quantum number")
    .def_rw("excitation_energy", &State::excitation_energy,
            "Excitation energy of the state with respect to the ground state "
            "in MeV")
    .def("__repr__", &State::str_rep);
}
