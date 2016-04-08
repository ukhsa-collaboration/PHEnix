API Reference
=============

Filters
-------

One of the key aims of this project is to provide a flexible and extendable interface for filtering VCFs.
While there are a number of tools that will process a VCF file, we felt that none offered in depth understanding
that we wanted. Instead we have based this library on PyVCF and extended with variety of filters :doc:`/Filters`.

If you feel that there is a filter that has not been implemented, or you want to try your hand at implementing object oriented
classes, then please see :ref:`implementing-filter`.

Mappers
-------

We found there was no single neat way of calling different mappers (there is BioPython, of course) but it wasn't
feasible to integrate it into our software stack in the way that we wanted to. Instead, we wrote a simple interfaces
that can be easily extended to add new mappers. See :ref:`implementing-mapper`.

Variant Callers
---------------

Just like with mappers, there was no neat way of incorporating other interfaces into this software stack. To add your favourite variant caller see :ref:`implementing-variant`.

.. toctree::
   :maxdepth: 3

   phe
