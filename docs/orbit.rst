Orbit Model
===========

A brief description here

Equations
---------

.. math::

    \begin{aligned}
        v &= \sqrt{\frac{\mu}{R_E + H}}  \\
        \Delta t_{\mathrm{orbit}} &= 2 \pi \sqrt{\frac{(R_E + H)^3}{\mu}} = \frac{2\pi (R_E+H}{v} \\
        \Delta t_{\mathrm{eclipse}} &= \frac{\Delta t_{\mathrm{orbit}}}{\pi}\arcsin\left(\frac{R_E}{R_E+H}\right)
    \end{aligned}

.. math::

    \theta_{\mathrm{slew}} = \arctan\left( \frac{ \sin\left( \frac{ \phi_{\mathrm{target}} }{R_E} \right)}{1-\cos \left(\frac{\phi_{\mathrm{target}}}{R_E} \right)+\frac{H}{R_E}} \right)



Medium Fidelity
---------------

A description of the medium fidelity models involved
