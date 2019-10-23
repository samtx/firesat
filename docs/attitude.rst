
Attitude Control Model
======================

A brief description here

Equations
---------

.. math::

    \tau_{slew} = \frac{4 \theta_{slew}}{(\Delta t_{slew})^2}I_{max}

.. math::

    \tau_g = \frac{3\mu}{2(R_E + H)^3} \lvert I_{max}-I_{min} \rvert \sin(2\theta)

.. math::

    \tau_{sp} = L_{sp} \frac{F_s}{C} A_s (1+q) \cos(i)

.. math::

    \tau_m = \frac{2 M R_D}{(R_E + H)^3}

.. math::

    \tau_a = \frac{1}{2} L_a \rho C_d A v^2

.. math::

    \begin{aligned}
        \tau_{dist} &= \sqrt{\tau_g^2 + \tau_{sp}^2 + \tau_m^2 + \tau_a^2} \\
        \tau_{tot} &= \max(\tau_{slew}, \tau_{dist})
    \end{aligned}

.. math::

    P_{ACS} = \tau_{tot} \omega_{max} + n P_{hold}

Medium Fidelity
---------------

A description of the medium fidelity models involved
