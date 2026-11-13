=================
C++ API reference
=================

.. container:: breadcrumb

   You are here: :doc:`Working on PHARE <working_on_phare>` › C++ API reference

**In one sentence:** the C++ API reference is a set of web pages, listing every class and function, that
Doxygen builds from the comments in ``src/``. At this commit the site does not host it, so you build it yourself.

.. admonition:: Not hosted yet
   :class: in-preparation

   The Read the Docs build only runs Sphinx. It does not run Doxygen, so there is no online C++ reference
   next to these pages. Build it locally with the steps below.

What Doxygen is, and what you get
---------------------------------

Doxygen reads the C++ headers and sources and writes one page per class, file and namespace.
Inside a page you can click from a class to its members and to its base classes.
The settings are in ``doc/Doxyfile.in`` (:src:`doc/Doxyfile.in`). The front page text is ``doc/doxypages/main.md``.

.. mermaid::

   flowchart TD
     A["src/ (headers and sources)"] --> D["Doxygen"]
     B["doc/Doxyfile.in<br/>(settings)"] --> D
     D --> H["html/index.html<br/>open it in a browser"]
     D --> L["latex/<br/>(also produced, not needed)"]

Build it locally
----------------

Run these commands from the ``doc/`` folder of the repository. Doxygen must be installed.

.. code-block:: bash

   cd doc
   mkdir -p build/doxygen
   sed "s|@CMAKE_CURRENT_BINARY_DIR@|$PWD/build/doxygen|" Doxyfile.in > Doxyfile
   doxygen Doxyfile
   # then open build/doxygen/html/index.html

Why the ``sed`` line: ``Doxyfile.in`` is a template. It contains the placeholder ``@CMAKE_CURRENT_BINARY_DIR@``
for the output folder, and nothing in the build system fills it in at this commit.
The ``Doxygen-build`` target of ``doc/Makefile`` (``doxygen Doxyfile``) expects the filled-in ``Doxyfile``, so run it
after the ``sed`` line if you prefer ``make Doxygen-build``.
Run from ``doc/``, because the input path in the settings (``../src``) is relative.

Do not commit the generated ``Doxyfile`` or the output folder.

Reading the result
------------------

* A class page shows its template parameters, its public and protected members (private members are not listed) and the description someone wrote for it.
* Classes and their members are listed even without a comment, just without a description. Free functions only appear
  when their file or namespace is documented (the settings have ``EXTRACT_ALL = NO``). Many entities have no comment yet
  (comments are being added later, pull request by pull request), so the descriptions are incomplete. The concept pages of this
  guide are the explanation in the meantime, and the ``src/...`` links in each page go to the code itself.
* Doxygen prints a long list of warnings for undocumented members. They are expected at this stage.
* Use the search box of the generated site to jump to a class name, then compare with the concept page:
  :doc:`gridlayout`, :doc:`fields`, :doc:`particles_ions`, :doc:`messengers`.

Go deeper
---------

* :doc:`working_on_phare`: build, tests and conventions for contributors.
* :doc:`how_to_add`: where to put new code.

Related
-------

* :doc:`big_map`: the layers of the code, each box links to its concept page.
