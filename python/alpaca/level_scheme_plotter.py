# This file is part of alpaca.

# alpaca is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# alpaca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with alpaca.  If not, see <https://www.gnu.org/licenses/>.

# Copyright (C) 2021 Udo Friman-Gayer

import matplotlib.pyplot as plt
import numpy as np


class LevelSchemePlotter:
    """Class to plot a labeled level scheme with a single excitation and a decay cascade

    This class is intended to visualize the transition cascades whose angular correlations are
    given by the AngularCorrelation class.
    The first transition, which defines the 'orientation' of the system, is assumed to be an
    excitation from the ground state, or, more generally, an initial state.
    The excitation populates the first state of the decay cascade.
    All other transitions belong to a cascade whose last transition is assumed to be observed.

    The visualization is realized using the matplotlib package.
    LevelSchemePlotter draws the level scheme on a matplotlib.axes.Axes object that is passed to its __init__ function.

    A plotted level scheme contains the following elements:

    - Horizontal lines that represent states. The initial state and the state that is reached by the exciting transition are represented by longer lines.
    - State labels that indicate the spin- and parity quantum numbers. For the initial state and the first state of the cascade, the labels are printed on the left side of the state line. For all other states, they are on the right side.
    - One upward arrow that represents the excitation.
    - Downward arrows that represent the decays.
    - Transition labels that indicate the EM character and the multipolarity. Whether they are placed left or right of the transition arrows is determined by the same rule as for the state labels.
    - Transition labels that indicate the Multipole-mixing ratio of that transition. Whether they are placed left or right of the transition arrows is determined by the same rule as for the state labels.

    All labels are in a LaTeX environment, i.e. mathematical expressions are valid.

    During the implementation of this class, it was found that the customization of matplotlib's 
    arrows is a bit awkward.
    To avoid things like arrow heads piercing through state lines and to be able to customize the width 
    of the arrow heads and the lines independently, the present implementation combines a line ('pyplot.plot') with an arrow head ('Axes.arrow').

    Attributes
    ----------

    ax: matplotlib.axes.Axes object
        Axes on which the level scheme should be drawn
    initial_state: State
        Initial state of the cascade. Same format as the corresponding input for the AngularCorrelation class.
    cascade_steps: array of [Transition, State] pairs
        Cascade steps, given as a list of arbitrary length which contains Transition-State pairs.
        The first and the last transition of this list are assumed to be observed.
        Same format as the corresponding input for the AngularCorrelation class.
    min_x, max_x, min_y, max_y: float
        Limits of the x - and y axis as given by ax.
    range_x, range_y:
        Range of the x - and y axis.
    delta_labels: list of str
        Labels for the multipole mixing ratios. The length of the list should equal the number of cascade steps.
    returns_to_initial_state: bool
        Determines whether the last step of the cascade should go back to the ground state.
        This option allows to make the distinction between a ground-state decay and a cascade that
        ends up in a state with the same quantum numbers as the ground state. (default: False)
    show_polarization: list of bool
        Determines which transition labels should indicate a polarization.
        This allows to indicate for which transition polarization information is available.
    fontsize: int or float
        Font size of the figure (matplotlib's 'fontsize' option). Several other font sizes are scaled to this value.
    fontsize_single_multipole: float
        Font size for transition labels (matplotlib's 'fontsize 'option) with a single multipole (default: 0.9*fontsize).
    fontsize_single_multipole: float
        Font size for transition labels (matplotlib's 'fontsize' option) with two possible multipoles. (default: 1.2*fontsize)
    state_line_width: float
        Line width for a state (matplotlib's 'lw' option)
    state_x: float
        Position on the x axis where the line of the state starts.
    state_width: float
        Length of a state line in x direction.
    intermediate_state_x: float
        Position on the x axis where the line of an intermediate state starts.
    intermediate_state_width: float
        Length of an intermediate-state line in x direction.
    initial_state_y: float
        Position of the initial state on the y axis.
    initial_state_y: float
        Position of the first excited state on the y axis.
    state_label_left_x, state_label_right_x:
        Position of the state labels on the x axis.
    parity_variable_symbol: str
        Symbol to be displayed in the state label if the parity of a state is unknown (default: '\pm').
    arrow_width: float
        With of the transition arrows (matplotlib's 'linewidth' option). (default: 2)
    excitation_arrow_x, decay_arrow_x: float
        Position of the excitation/decay arrow on the x axis.
    arrow_head_length: float
        Length of the arrow heads, scaled to the plot range (Axes.arrow's 'head_length' option). (default: 0.04*range_y)
    arrow_head_width: float
        Width of the arrow heads, scaled to arrow_width. (default: 0.03*arrow_width)
    excitation_arrow_color, decay_arrow_color: matplotlib color specification
        Colors of the excitation- and decay arrows. Must be a color specification that can be understood by matplotlib.
    em_variable_symbol: str
        Symbol to be displayed in the transition label if the EM character of a transition is unknown (default: '\sigma').
    decay_label_right_x, excitation_label_left_x: float
        Position of the decay labels on the x axis.
    delta_label_rotation: float
        Rotation of the multipole-mixing ratio labels in degrees (Axes.text's 'rotation' option). The default is 90 degrees, which means that the labels are oriented as if the corresponding arrow were the base line of the text.
    delta_label_left_x, delta_label_right_x: float
        Position of the multipole-mixing labels on the x axis.
    zorder_states, zorder_arrows: int
        Sets the order in which the state lines and arrows are drawn (matplotlib's 'zorder' option). (default: zorder_states < zorder_arrows, i.e. arrows are on top of state lines)
    """

    def __init__(
        self,
        axis,
        initial_state,
        cascade_steps,
        delta_labels,
        returns_to_initial_state=False,
        show_polarization=None,
        fontsize=12,
        state_line_width=2,
        arrow_width=2,
        offset=(0, 0),
        delta_label_rotation=90,
        em_variable_symbol="\sigma",
        parity_variable_symbol="\pm",
    ):
    """Initialization
    
    ax: matplotlib.axes.Axes object
        Axes on which the level scheme should be drawn
    initial_state: State
        Initial state of the cascade. Same format as the corresponding input for the AngularCorrelation class.
    cascade_steps: array of [Transition, State] pairs
        Cascade steps, given as a list of arbitrary length which contains Transition-State pairs.
        The first and the last transition of this list are assumed to be observed.
        Same format as the corresponding input for the AngularCorrelation class.
    delta_labels: list of str
        Labels for the multipole mixing ratios. The length of the list should equal the number of cascade steps.
    returns_to_initial_state: bool
        Determines whether the last step of the cascade should go back to the ground state.
        This option allows to make the distinction between a ground-state decay and a cascade that
        ends up in a state with the same quantum numbers as the ground state. (default: False)
    show_polarization: list of bool
        Determines which transition labels should indicate a polarization.
        This allows to indicate for which transition polarization information is available.
    fontsize: int or float
        Font size of the figure (matplotlib's 'fontsize' option). Several other font sizes are scaled to this value.
    state_line_width: float
        Line width for a state (matplotlib's 'lw' option)
    arrow_width: float
        With of the transition arrows (matplotlib's 'linewidth' option). (default: 2)
    offset: (float, float)
        Offset of the level scheme from the lower-left corner of the panel in x- and y direction. Even if the parameters of the LevelSchemePlotter have already been chosen to yield a nice-looking figure, an offset may sometimes be necessary. (default: (0, 0), i.e. no offset).
    delta_label_rotation: float
        Rotation of the multipole-mixing ratio labels in degrees (Axes.text's 'rotation' option). The default is 90 degrees, which means that the labels are oriented as if the corresponding arrow were the base line of the text.
    em_variable_symbol: str
        Symbol to be displayed in the transition label if the EM character of a transition is unknown (default: '\sigma').
    parity_variable_symbol: str
        Symbol to be displayed in the state label if the parity of a state is unknown (default: '\pm'). 
    """
        self.ax = axis
        self.min_x, self.max_x = axis.get_xlim()
        self.range_x = self.max_x - self.min_x
        self.min_y, self.max_y = axis.get_ylim()
        self.range_y = self.max_y - self.min_y
        self.initial_state = initial_state
        self.cascade_steps = cascade_steps
        self.delta_labels = delta_labels
        self.returns_to_initial_state = returns_to_initial_state
        self.show_polarization = [False] * len(cascade_steps)
        if show_polarization is not None:
            self.show_polarization = show_polarization

        ## Parameters for the plot
        # Fonts
        self.fontsize = fontsize
        self.fontsize_single_multipole = 0.9 * fontsize
        self.fontsize_two_multipoles = 1.2 * fontsize

        # State lines
        self.state_line_width = state_line_width
        self.state_x = 0.4 * self.range_x + self.min_x + offset[0] * self.range_x
        self.state_width = 0.4 * self.range_x
        self.intermediate_state_x = (
            0.55 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.intermediate_state_width = 0.25 * self.range_x

        self.initial_state_y = (
            0.2 * self.range_y + self.min_y + offset[1] * self.range_y
        )
        self.excited_state_y = (
            0.8 * self.range_y + self.min_y + offset[1] * self.range_y
        )

        # State labels
        self.state_label_left_x = (
            0.25 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.state_label_right_x = (
            0.9 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.parity_variable_symbol = parity_variable_symbol

        # Transition arrows
        self.arrow_width = arrow_width
        self.excitation_arrow_x = (
            0.5 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.decay_arrow_x = 0.7 * self.range_x + self.min_x + offset[0] * self.range_x
        self.arrow_head_length = 0.04 * self.range_y
        self.arrow_head_width = 0.03 * self.arrow_width
        self.excitation_arrow_color = "blue"
        self.decay_arrow_color = "red"

        # Transition labels
        self.em_variable_symbol = em_variable_symbol
        self.decay_label_right_x = (
            0.85 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.excitation_label_left_x = (
            0.18 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.delta_label_rotation = delta_label_rotation

        # Multipole mixing ratio (delta) labels
        self.delta_label_left_x = (
            0.4 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.delta_label_right_x = (
            0.74 * self.range_x + self.min_x + offset[0] * self.range_x
        )

        # Order of drawing
        self.zorder_states = 0
        self.zorder_arrows = 1

    def plot(self):
        # Initial and excited state
        self.ax.plot(
            [self.state_x, self.state_x + self.state_width],
            [self.initial_state_y] * 2,
            color="black",
            linewidth=self.state_line_width,
            zorder=self.zorder_states,
        )
        self.ax.text(
            self.state_label_left_x,
            self.initial_state_y,
            self.initial_state.tex(parity_variable_symbol=self.parity_variable_symbol),
            verticalalignment="center",
            fontsize=self.fontsize,
        )
        self.ax.plot(
            [self.state_x, self.state_x + self.state_width],
            [self.excited_state_y] * 2,
            color="black",
            linewidth=self.state_line_width,
            zorder=self.zorder_states,
        )
        self.ax.text(
            self.state_label_left_x,
            self.excited_state_y,
            self.cascade_steps[0][1].tex(parity_variable_symbol=self.parity_variable_symbol),
            verticalalignment="center",
            fontsize=self.fontsize,
        )

        # Excitation
        self.ax.plot(
            [self.excitation_arrow_x] * 2,
            [self.initial_state_y, self.excited_state_y - self.arrow_head_length],
            "-",
            linewidth=self.arrow_width,
            color=self.excitation_arrow_color,
            zorder=self.zorder_arrows,
        )
        self.ax.arrow(
            self.excitation_arrow_x,
            self.excited_state_y - self.arrow_head_length,
            0.0,
            1.0e-5 * self.range_y,
            head_length=self.arrow_head_length,
            head_width=self.arrow_head_width,
            facecolor=self.excitation_arrow_color,
            edgecolor=self.excitation_arrow_color,
            zorder=self.zorder_arrows,
        )
        self.ax.text(
            self.excitation_label_left_x,
            0.5 * (self.excited_state_y - self.initial_state_y) + self.initial_state_y,
            self.cascade_steps[0][0].tex(
                em_variable_symbol=self.em_variable_symbol,
                always_show_secondary=False,
                show_polarization=self.show_polarization[0],
            ),
            verticalalignment="center",
            fontsize=self.fontsize_single_multipole
            if self.cascade_steps[0][0].delta == 0.0
            else self.fontsize_two_multipoles,
        )
        self.ax.text(
            self.delta_label_left_x,
            0.5 * (self.excited_state_y - self.initial_state_y) + self.initial_state_y,
            self.delta_labels[0],
            verticalalignment="center",
            fontsize=self.fontsize,
            rotation=self.delta_label_rotation,
        )

        # Calculate position of states in decay cascade
        n_decay_steps = len(self.cascade_steps) - 1
        excitation_delta_y = self.excited_state_y - self.initial_state_y
        cascade_states_y = np.arange(1, n_decay_steps + 1)

        if self.returns_to_initial_state:
            cascade_states_y = (
                self.excited_state_y
                - cascade_states_y * excitation_delta_y / (n_decay_steps)
            )
        else:
            cascade_states_y = (
                self.excited_state_y
                - cascade_states_y * excitation_delta_y / (n_decay_steps + 1)
            )

        # States in cascade
        for i in range(
            n_decay_steps if not self.returns_to_initial_state else n_decay_steps - 1
        ):
            self.ax.plot(
                [
                    self.intermediate_state_x,
                    self.intermediate_state_x + self.intermediate_state_width,
                ],
                [cascade_states_y[i]] * 2,
                color="black",
                linewidth=self.state_line_width,
                zorder=self.zorder_states,
            )
            self.ax.text(
                self.state_label_right_x,
                cascade_states_y[i],
                self.cascade_steps[i + 1][1].tex(
                    parity_variable_symbol=self.parity_variable_symbol
                ),
                verticalalignment="center",
                fontsize=self.fontsize,
            )

        # First transition in cascade
        self.ax.plot(
            [self.decay_arrow_x] * 2,
            [self.excited_state_y, cascade_states_y[0] + self.arrow_head_length],
            "--" if n_decay_steps > 1 else "-",
            linewidth=self.arrow_width,
            color=self.decay_arrow_color,
            zorder=self.zorder_arrows,
        )
        self.ax.arrow(
            self.decay_arrow_x,
            cascade_states_y[0] + self.arrow_head_length,
            0.0,
            -1.0e-5 * self.range_y,
            head_length=self.arrow_head_length,
            head_width=self.arrow_head_width,
            color=self.decay_arrow_color,
            zorder=self.zorder_arrows,
        )
        self.ax.text(
            self.decay_label_right_x,
            0.5 * (self.excited_state_y - cascade_states_y[0]) + cascade_states_y[0],
            self.cascade_steps[1][0].tex(
                em_variable_symbol=self.em_variable_symbol,
                always_show_secondary=False,
                show_polarization=self.show_polarization[1],
            ),
            verticalalignment="center",
            fontsize=self.fontsize_single_multipole
            if self.cascade_steps[1][0].delta == 0.0
            else self.fontsize_two_multipoles,
        )
        self.ax.text(
            self.delta_label_right_x,
            0.5 * (self.excited_state_y - cascade_states_y[0]) + cascade_states_y[0],
            self.delta_labels[1],
            verticalalignment="center",
            fontsize=self.fontsize,
            rotation=self.delta_label_rotation,
        )

        # Transitions in cascade
        for i in range(1, n_decay_steps):
            self.ax.plot(
                [self.decay_arrow_x] * 2,
                [cascade_states_y[i - 1], cascade_states_y[i] + self.arrow_head_length],
                "--" if i < n_decay_steps - 1 else "-",
                linewidth=self.arrow_width,
                color=self.decay_arrow_color,
                zorder=self.zorder_arrows,
            )
            self.ax.arrow(
                self.decay_arrow_x,
                cascade_states_y[i] + self.arrow_head_length,
                0.0,
                -1.0e-5 * self.range_y,
                head_length=self.arrow_head_length,
                head_width=self.arrow_head_width,
                color=self.decay_arrow_color,
                zorder=self.zorder_arrows,
            )
            self.ax.text(
                self.decay_label_right_x,
                0.5 * (cascade_states_y[i - 1] - cascade_states_y[i])
                + cascade_states_y[i],
                self.cascade_steps[i + 1][0].tex(
                    em_variable_symbol=self.em_variable_symbol,
                    always_show_secondary=False,
                    show_polarization=self.show_polarization[i + 1],
                ),
                verticalalignment="center",
                fontsize=self.fontsize_single_multipole
                if self.cascade_steps[i + 1][0].delta == 0.0
                else self.fontsize_two_multipoles,
            )
            self.ax.text(
                self.delta_label_right_x,
                0.5 * (self.excited_state_y - cascade_states_y[i - 1])
                + cascade_states_y[i],
                self.delta_labels[i + 1],
                verticalalignment="center",
                fontsize=self.fontsize,
                rotation=self.delta_label_rotation,
            )
