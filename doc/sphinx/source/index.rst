:html_theme.sidebar_secondary.remove:

.. ESMValTool documentation master file, created by
   sphinx-quickstart on Tue Jun  2 11:34:13 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ESMValTool's documentation!
======================================

**ESMValTool** is a community diagnostics and performance metrics tool for the evaluation of Earth System Models (ESMs) that allows for
routine comparison of models and observations. It includes a large collection of community recipes and observation data formatters to CMOR standards.

**ESMValCore** is a software package which provides the core functionality for ESMValTool. It is a workflow to find CMIP data, and apply
commonly used pre-processing functions.

To get a first impression of what ESMValTool and ESMValCore can do for you,
have a look at our blog posts
`Analysis-ready climate data with ESMValCore <https://blog.esciencecenter.nl/easy-ipcc-powered-by-esmvalcore-19a0b6366ea7>`_
and
`ESMValTool: Recipes for solid climate science <https://blog.esciencecenter.nl/esmvaltool-recipes-for-solid-climate-science-da5b33814f69>`_.


Basic documentation schema
--------------------------
This gives a brief idea of topics in each of the packages to help find information. As ESMValTool encompasses ESMValCore there will be some overlap between them.
The documentation navigation to the left provides more detailed information about ESMValTool.
`ESMValCore documentation <https://docs.esmvaltool.org/projects/ESMValCore/en/latest/index.html>`_ is also available. ESMValCore can be used without ESMValTool such as in Jupyter
notebooks.

.. container::
   :name: figarch

   .. figure:: figures/ESMValSchemaDiagram.png
      :alt: Brief topics for tool and core.
      :figclass: align-center
      :class: dark-light


Learning resources:
-------------------

A tutorial is available on https://tutorial.esmvaltool.org.

A series of video lectures has been created by `ACCESS-NRI <https://www.access-nri.org.au>`_.
While these are tailored for ACCESS users, they are still very informative.

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/videoseries?si=pUXrXB8C8bLRfQHY&amp;list=PLFjfi2xLaFpJp59LvDc1upQsj_xzFlFLc" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

How to contribute
-----------------
Refer to ESMValTool :ref:`contributing to the community<community>` for a guide on contributing recipes and diagnostics.

Refer to ESMValCore contributing for information on contributing code.

Get in touch!
-------------
Contact information is available :ref:`here <Support-and-Contact>`.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: ESMValTool

    Introduction <introduction>
    ESMValTool Functionalities <functionalities>
    Gallery <gallery>
    Getting started <quickstart/index>
    Recipes <recipes/index>
    Diagnostics API Reference <api/esmvaltool>
    Obtaining input data <input>
    Making a recipe or diagnostic <develop/index>
    Contributing to the community <community/index>
    Utilities <utils/utils>
    Frequently Asked Questions <faq>
    Changelog <changelog>

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
