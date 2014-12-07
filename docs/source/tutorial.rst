Simulation Structure
------------------------------------------------

Simulations with EMUstack are generally carried out using a python script file.
This file is kept in its own directory which is placed in the EMUstack directory.
All results of the simulation are automatically created within this directory. This directory then serves as a complete record of the calculation. Often, we will also save the simulation objects (scattering matrices, propagation constants etc.) within this folder for future inspection, manipulation, plotting, etc.
Traditionally the name of the python script file begins with simo\_. This is convenient for setting terminal alias' for running the script.
Throughout the tutorial the script file will be called simo.py.

To start a simulation open a terminal and change into the directory containing the simo.py file.
To run this script::

    $ python simo.py

To have direct access to the simulation objects upon the completion of the script use,::

    $ python -i simo.py

This will return you into an interactive python session in which all simulation objects are accessible. In this session you can access the docstrings of objects, classes and methods. For example::

    >>> from pydoc import help
    >>> help(objects.Light)

where we have accessed the docstring of the Light class from objects.py


In the remainder of the guide we go through a number of example simo.py files. These cover a wide range (though non-exhaustive) of established applications of EMUstack. The source files for these examples are in EMUstack/examples/
The first 8 examples are pretty essential for using EMUstack, while those thereafter show EMUstack applied to a number of (IMHO) interesting situations.

Another tip to mention before diving into the examples is running simulations within :ref:`screen_sesh`. These allow you to disconnect from the terminal instance and are discusses in :ref:`screen_sesh`.




Single Interface
------------------------------------------------

.. literalinclude:: ../../examples/simo_010-single_interface.py
    :lines: 20-


Dispersion & Parallel Computation
------------------------------------------------

.. literalinclude:: ../../examples/simo_011-single_interface-dispersive.py
    :lines: 20-


Thin Film Stack
------------------------------------------------

.. literalinclude:: ../../examples/simo_020-thin_film_multilayered_stack.py
    :lines: 20-


Including Metals
------------------------------------------------

.. literalinclude:: ../../examples/simo_021-thin_film_mirror.py
    :lines: 20-


1D Grating
------------------------------------------------

.. literalinclude:: ../../examples/simo_030-1D_grating.py
    :lines: 20-


2D Grating
------------------------------------------------

.. literalinclude:: ../../examples/simo_040-2D_array.py
    :lines: 20-


Angles of Incidence & Eliptical Inclusions
------------------------------------------------

.. literalinclude:: ../../examples/simo_042-eliptical_holes-CD.py
    :lines: 20-


Plotting Fields 1D
------------------------------------------------

.. literalinclude:: ../../examples/simo_050-plotting_fields_1d.py
    :lines: 20-


Plotting Fields 2D
------------------------------------------------

.. literalinclude:: ../../examples/simo_051-plotting_fields_2d.py
    :lines: 20-


Plotting Amplitudes
------------------------------------------------

.. literalinclude:: ../../examples/simo_052-plotting_amplitudes.py
    :lines: 20-


Shear Transformations
------------------------------------------------

.. literalinclude:: ../../examples/simo_060-shear_transformations.py
    :lines: 20-


Ultrathin Absorption Limit - Varying n
------------------------------------------------

.. literalinclude:: ../../examples/simo_070-ultrathin_limit.py
    :lines: 20-


Varying a Layer of a Stack
------------------------------------------------

.. literalinclude:: ../../examples/simo_071-many_substrates.py
    :lines: 20-


Convergence Testing
------------------------------------------------

.. literalinclude:: ../../examples/simo_080-convergence-stacked_gratings.py
    :lines: 20-


Extraordinary Optical Transmission
------------------------------------------------

.. literalinclude:: ../../examples/simo_090-EOT.py
    :lines: 20-


Screen Sessions
------------------------------------------------
::

    screen

is an extremely useful little linux command. In the context of long-ish calculations it has two important applications; ensuring your calculation is unaffected if your connection to a remote machine breaks, and terminating calculations that have hung without closing the terminal.
For more information see the manual::

    $ man screen

or see online discussions `here <http://www.howtoforge.com/linux_screen>`_, `and here <http://www.rackaid.com/blog/linux-screen-tutorial-and-how-to/>`_.


The screen session or also called screen instance looks just like your regular terminal/putty, but you can disconnect from it (close putty, turn off your computer etc.) and later reconnect to the screen session and everything inside of this will have kept running. You can also reconnect to the session from a different computer via ssh.

Basic Usage
,,,,,,,,,,,,,,,,,,,,,

To install screen::

    $ sudo apt-get install screen

To open a new screen session::

    $ screen

We can start a new calculation here::

    $ cd EMUstack/examples/
    $ python simo_040-2D_array.py

We can then detach from the session (leaving everything in the screen running) by typing::

    Ctrl +a
    Ctrl +d

We can now monitor the processes in that session::

    $ top

Where we note the numerous running python processes that EMUstack has started. Watching the number of processes is useful for checking if a long simulation is near completion (which is indicated by the number of processes dropping to less than the specified num_cores).

We could now start another screen and run some more calculations in this terminal (or do anything else).
If we want to access the first session we 'reattach' by typing::

    Ctrl +a +r

Or entering the following into the terminal::

    $ screen -r

If there are multiple sessions use::

    $ screen -ls

to get a listing of the sessions and their ID numbers. To reattach to a particular screen, with ID 1221::

    $ screen -r 1221

To terminate a screen from within type::

    Ctrl+d

Or, taking the session ID from the previous example::

    screen -X -S 1221 kill



Terminating EMU stacks
,,,,,,,,,,,,,,,,,,,,,,,


If (for some estranged reason) a simulation hangs, we can kill all python instances upon the machine::

    $ pkill python

If a calculation hangs from within a screen session one must first detach from that session then kill python. A more targeted way to kill processes is using their PID::

    $ kill PID

Or if this does not suffice be a little more forceful::

    $ kill -9 PID

The PID is found from one of two ways::

    $ top
    $ ps -fe | grep username

.. Screen Sessions
.. ------------------------------------------------
.. .. toctree::
..     :maxdepth: 4

..     screen_sesh

