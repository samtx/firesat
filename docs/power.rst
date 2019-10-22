.. Firesat documentation master file, created by
   sphinx-quickstart on Tue Oct 22 08:21:02 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Power Model
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Equations are

.. math::

    P_{sa} = \frac{1}{T_d}\left(\frac{P_e \ \Delta t_\mathrm{e}}{X_e}+\frac{P_d \ \Delta t_\mathrm{d}}{X_d}\right)

Where :math:`P_e` is the power required by the satellite during an eclipse.

.. math::

    P_{BOL} = \eta F_s I_d \cos(i)

.. math::

    P_{EOL} = P_{BOL} (1 - \epsilon_{deg})^{LT}

.. math::

    A_{sa} = \frac{P_{sa}}{P_{EOL}}

.. math::

    L = \sqrt{\frac{A_{sa}r_{lw}}{m_{sa}}} \qquad \qquad W = \sqrt{\frac{A_{sa}}{r_{lw}m_{sa}}}

.. math::

    m_{sa} = 2 \rho_{sa} L W t

.. math::

    \begin{aligned}

    I_{sax} &= m_{sa} \Big[ \frac{1}{12} (L^2 + t^2) + \left(D + 0.5 L  \right)^2 \Big] \\
    I_{saz} &= \frac{m_{sa}}{12}(t^2 + W^2) \\
    I_{say} &= m_{sa} \Big[ \frac{1}{12} (L^2 + t^2) + \left(D + 0.5 L  \right)^2 \Big]

    \end{aligned}

.. math::

    \begin{aligned}
        I_{tot} & = I_{sa} + I_{body} \\
        I_{max} & = \max(I_{totx}, I_{toty}, I_{totz}) \\
        I_{min} & = \min(I_{totx}, I_{toty}, I_{totz})
    \end{aligned}



Medium Fidelity
---------------

A description of the medium fidelity settings and equations




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
