.. :orphan:  # uncomment before merging!

.. How to use this template
..
.. 1. Make a copy of this file.
.. 2. Update the filename.
.. 3. Add the filename (without the '.rst' extension)
..    to a 'toctree' directive
..    in the 'tutorials/index.rst' file.
.. 4. Use the name of the tutorial for the title below.
.. 5. Remove ':orphan:', these comments,
..    and 'The language of tutorials' section (below).
..
.. https://diataxis.fr/tutorials provides more information.

Title of the tutorial
=====================

.. admonition:: Overview
   :class: note

   .. grid:: 2
      :gutter: 1
      :margin: 3 3 0 5

      .. grid-item-card:: Timings
         :columns: 12

         * Teaching: X min
         * Exercises: Y min

      .. grid-item-card:: Questions

         * Describe what the learner will accomplish.
         * Question 1?
         * Question 2?

      .. grid-item-card:: Learning outcomes

         * Describe the learning outcomes, i.e.
           the knowledge, skills, or expertise the learner will gain.
         * Objective 1.
         * Objective 2.

      .. grid-item-card:: Prerequisites

         * Provide any prerequisites.
         * Prerequisites 1.
         * Prerequisites 2.

      .. grid-item-card:: Assumptions

         * Provide any assumptions.
         * Assumption 1.
         * Assumption 2.

High level step 1
-----------------

Add content here.

* Consider including diagrams to support the text.
* Use any free icons from `Font Awesome`_ via an ``:fa:`` directive.
* Use admonitions from `PyData Theme documentation: Admonitions`_.
* Include other files from within the documentation:

  .. include:: files/esmvaltool_output_header.txt
     :code:

High level step 2
-----------------

Add content here.

* Use code snippets:

  .. code-block:: bash
     :caption: Bash

     my command

  .. code-block:: python
     :caption: Python
     :linenos:
     :emphasize-lines: 2

     with line numbers
     and code highlighting

* Create expandable sections containing the
  output / answer / solution from a command:

  .. dropdown:: Output
     :color: secondary
     :icon: eye

     .. code-block:: bash
        :caption: Bash

        my output

The language of tutorials
-------------------------

We ...
    The first-person plural
    affirms the relationship between tutor and learner:
    you are not alone; we are in this together.

First, do x. Now, do y. Now that you have done y, do z.
    No room for ambiguity or doubt.

We must always do x before we do y. <explanation> provides more details.
    Provide minimal explanation of actions
    in the most basic language possible.
    Link to more detailed explanation.

The output should be something like ...
    Give your learner clear expectations.

Notice that ... Remember that ... Let's check ...
    Give your learner plenty of clues
    to help confirm they are on the right track and orient themselves.

Congratulations!
----------------

You have built a secure, three-layer hylomorphic stasis engine ...
    Describe (and admire, in a mild way) what your learner has accomplished.

.. admonition:: Key points
   :class: important

   * Key point 1
   * Key point 2

.. _`PyData Theme documentation: Admonitions`: https://pydata-sphinx-theme.readthedocs.io/en/stable/examples/kitchen-sink/admonitions.html
.. _`Font Awesome`: https://fontawesome.com/search?ic=free-collection
