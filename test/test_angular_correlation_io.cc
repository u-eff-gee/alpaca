/*
    This file is part of alpaca.

    alpaca is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    alpaca is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with alpaca.  If not, see <https://www.gnu.org/licenses/>.

    Copyright (C) 2021-2023 Udo Friman-Gayer
*/

#include <cassert>
#include <stdexcept>

#include "AngularCorrelation.hh"
#include "State.hh"
#include "Transition.hh"

using std::invalid_argument;

/**
 * \brief Test the input validation of the AngularCorrelation class.
 * 
 * Test by calling the constructor with various invalid inputs.
 * The basic example is a
 * 
 * 0^+ -> 1^+ -> 0^+
 * 
 * cascade.
 */
int main(){

    bool error_thrown = false;

    // Error: First EM character not given
    try{
        AngularCorrelation ang_corr(
            State(0, positive),
            {
                {
                    Transition(em_unknown, 2, magnetic, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: First EM character not given for second transition
    try{
        AngularCorrelation ang_corr(
            State(0, positive),
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(em_unknown, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Second EM character not given
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, em_unknown, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Second EM character not given for second transition
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, em_unknown, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: EM character given, but first parity missing
    try{
        AngularCorrelation ang_corr(
            State(0, parity_unknown), 
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: EM character given, but second parity missing
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(2, parity_unknown)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: EM character given, but parity missing for final state of second transition
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, parity_unknown)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Triangle inequality violated for first transition
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(magnetic, 4, electric, 6, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Triangle inequality violated for second transition
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(magnetic, 2, electric, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 4, magnetic, 6, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Not an error: Triangle inequality fulfilled by second transition
    AngularCorrelation ang_corr(
        State(0, positive), 
        {
            {
                Transition(electric, 10, electric, 2, 0.),
                State(2, negative)
            },
            {
                Transition(electric, 2, magnetic, 4, 0.),
                State(0, positive)
            }
        }
    );

    // Error: First electromagnetic character wrong
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(magnetic, 2, magnetic, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(2, positive)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(magnetic, 4, magnetic, 6, 0.),
                    State(4, positive)
                },
                {
                    Transition(electric, 4, magnetic, 6, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: First electromagnetic character wrong for second transition
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(magnetic, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Second electromagnetic character wrong
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, electric, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Second electromagnetic character wrong for second transition
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(2, negative)
                },
                {
                    Transition(electric, 2, electric, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;


    // Error: Mixing of half-integer and integer spins
    try{
        AngularCorrelation ang_corr(
            State(0, positive), 
            {
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(1, negative)
                },
                {
                    Transition(electric, 2, magnetic, 4, 0.),
                    State(0, positive)
                }
            }
        );
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }


    assert(error_thrown);
    error_thrown = false;

    // Simplified constructor. Check transition inference with various possibilities.
    // Error: Transition between spin-0 states
    try{
        AngularCorrelation ang_corr_0_1_0{
            State(0, parity_unknown), 
            {
                State(0, negative),
                State(0, parity_unknown)
            }
        }; 
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Too few steps in cascade
    try{
        AngularCorrelation ang_corr_0_1{
		State(0, parity_unknown), 
            {
                State(2, negative),
            }
        }; 
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Same as above with more general constructor.
    try{
        AngularCorrelation ang_corr_0p_1p{
		State(0, positive), 
            {
                {Transition(magnetic, 2, electric, 4, 0.), State(2, positive)},
            }
        }; 
    } catch(const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    AngularCorrelation ang_corr_0_1_0{
        State(0, positive), 
        {
            State(2, negative),
            State(6, negative)
        }
    };

    vector<pair<Transition, State>> cas_ste = ang_corr_0_1_0.get_cascade_steps();

    assert(cas_ste[0].first.two_L == 2);
    assert(cas_ste[0].first.em_char == electric);
    assert(cas_ste[0].first.two_Lp == 4);
    assert(cas_ste[0].first.em_charp == magnetic);
    assert(cas_ste[0].first.delta == 0.);
    assert(cas_ste[1].first.two_L == 4);
    assert(cas_ste[1].first.em_char == electric);
    assert(cas_ste[1].first.two_Lp == 6);
    assert(cas_ste[1].first.em_charp == magnetic);
    assert(cas_ste[1].first.delta == 0.);

    AngularCorrelation ang_corr_3_5_3_9 = AngularCorrelation{
        State(3, positive), 
        {
            State(5, parity_unknown),
            State(3, negative),
            State(9, negative)
        }
    };

    cas_ste = ang_corr_3_5_3_9.get_cascade_steps();

    assert(cas_ste[0].first.em_char == em_unknown);
    assert(cas_ste[0].first.em_charp == em_unknown);
    assert(cas_ste[1].first.em_char == em_unknown);
    assert(cas_ste[1].first.em_charp == em_unknown);
    assert(cas_ste[2].first.em_char == magnetic);
    assert(cas_ste[2].first.em_charp == electric);

    // Check that the correct transition is inferred between states of equal spin.
    // The triangle inequality gives a monopole transition as the lowest possible multipolarity
    // in that case, but this is not possible for a single photon.
    AngularCorrelation ang_corr_0_4_4 = AngularCorrelation{
        State(0, positive), 
        {
            State(4, positive),
            State(4, positive),
        }
    };

    cas_ste = ang_corr_0_4_4.get_cascade_steps();

    assert(cas_ste[0].first.em_char == electric);
    assert(cas_ste[0].first.two_L == 4);
    assert(cas_ste[0].first.em_charp == magnetic);
    assert(cas_ste[0].first.two_Lp == 6);
    assert(cas_ste[1].first.em_char == magnetic);
    assert(cas_ste[1].first.two_L == 2);
    assert(cas_ste[1].first.em_charp == electric);
    assert(cas_ste[1].first.two_Lp == 4);

   AngularCorrelation ang_corr_0_4m_4 = AngularCorrelation{
        State(0, positive), 
        {
            State(4, negative),
            State(4, positive),
        }
    };

    cas_ste = ang_corr_0_4m_4.get_cascade_steps();

    assert(cas_ste[0].first.em_char == magnetic);
    assert(cas_ste[0].first.two_L == 4);
    assert(cas_ste[0].first.em_charp == electric);
    assert(cas_ste[0].first.two_Lp == 6);
    assert(cas_ste[1].first.em_char == electric);
    assert(cas_ste[1].first.two_L == 2);
    assert(cas_ste[1].first.em_charp == magnetic);
    assert(cas_ste[1].first.two_Lp == 4);
}