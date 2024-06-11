.. _custom_reductions:
===============================
Writing New Reduction Processes
===============================

While many reduction processes are already provided, versatile use-cases may demand the user write a custom reduction
to accomplish their statistical analyses. In general, writing reductions is not a simple task; however, ``pyXMIP`` has some
structures which should make this a possible thing to accomplish given that the user has some mastery of the code structure.

What is a Reduction Process?
----------------------------

At a basic level, reduction processes are quite simple:

- They are classes (subclassed from :py:class:`structures.reduction.ReductionProcess`).
  - All reduction classes have a set of relevant parameters.
  - All reduction classes have a :py:meth:`structures.reduction.ReductionProcess.validate_params` method to validate parameter values.
  - All reduction classes have a ``__call__()`` method.

.. important::

    Reduction processes always implement ``__call__(self, cross_match_database, table)``, where ``cross_match_database`` is a
    :py:class:`cross_reference.CrossMatchDatabase` instance representing the data to be reduced and the ``table`` parameter is the
    table within that database on which the reduction is performed.

From the perspective of the end-user, one **sets up** the reduction process instance (in doing so, specifying all of the needed parameters), and then
runs it using ``instance(cross_match_database,table)``.

Writing Reductions
------------------

Basic Structure
+++++++++++++++

All reduction processes must be subclassed from :py:class:`structures.reduction.ReductionProcess`. This is an abstract base class with only
a couple of basic scaffolding methods already provided. Critically, the initialization takes only one parameter (``process_name``) and
any number of optional ``**kwargs``.

The process name becomes the :py:attr:`~structures.reduction.ReductionProcess.name` while the ``**kwargs`` become the underlying
``._parameter_dictionary`` protected attribute.

.. hint::

    The ``._parameter_dictionary`` doesn't need to be accessed directly. All of the provided ``kwargs`` go to this dictionary
    and are then accessed and written by the descriptor classes which represent the parameters.

.. rubric:: Example

As an example of a very basic structure for a reduction class, we might create something like the following:

.. code-block:: python

    from pyxmip.structures.reduction import ReductionProcess

    class MyProcess(ReductionProcess):
        """
        My special reduction process
        """

        def __init__(self,**kwargs):
            super().__init__("MyStandardName",**kwargs)

        def __call__(self,cross_match_database,table):
            return "Ooops! I'm not super useful yet!"

This is obviously a pretty useless example right now, but the basic structure is the point here.

Reduction Parameters
++++++++++++++++++++

What really makes reduction processes versatile are the :py:class:`structures.reduction.ReductionProcessParameter` descriptors that
can be added to any given reduction process.

.. hint::

    If you're not familiar with descriptors in python, you might benefit from reading `the guide at RealPython <https://realpython.com/python-descriptors/>`_.


Instances of this class represent the **core parameters of your reduction process**. You add them in your class as class-attributes:

.. code-block:: python

    class MyProcess(ReductionProcess):
        """
        My special reduction process
        """
        parameter_1 = ReductionProcessParameter()
        parameter_2 = ReductionProcessParameter()

        def __init__(self,**kwargs):
            super().__init__("MyStandardName",**kwargs)

        def __call__(self,cross_match_database,table):
            return "Ooops! I'm not super useful yet!"

Now our ``MyProcess`` reduction process has two parameters, both of which can be specified on instantiation:

.. code-block:: python

    >>> red_proc = MyProcess(parameter_1='hello',parameter_2='world')

.. note::

    Full documentation for these descriptors is provided at :py:class:`structures.reduction.ReductionProcessParameter`.
