bioutils -- bioinformatics utilities and lookup tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

|pypi_badge| |build_status| |cov_badge| |cc_badge| |issues_badge| |contributors| |license| |changelog|


bioutils provides some common utilities and lookup tables for bioinformatics.

* bioutils.accessions -- parse accessions, infer namespaces 
* bioutils.assemblies -- Human assembly information (from NCBI/GRCh)
* bioutils.cytobands -- map cytobands to coordinates (from UCSC cytoband tables)
* bioutils.digests -- implementations of various digests
* bioutils.normalize -- allele normalization (left shuffle, right shuffle, expanded, vcf)
  

To use an E-Utilities API key run add it to an environment variable
called `ncbi_api_key` and it will be used in the E-Utilities request.


.. |build_status| image:: https://travis-ci.org/biocommons/bioutils.svg?branch=master
  :target: https://travis-ci.org/biocommons/bioutils

.. |changelog| image:: https://img.shields.io/badge/docs-changelog-green.svg
   :target: https://bioutils.readthedocs.io

.. |contributors| image:: https://img.shields.io/github/contributors/biocommons/bioutils.svg
  :target: https://github.com/biocommons/bioutils

.. |docs| image:: https://img.shields.io/badge/docs-readthedocs-green.svg
   :target: http://bioutils.readthedocs.io/

.. |issues_badge| image:: https://img.shields.io/github/issues/biocommons/bioutils.png
  :target: https://github.com/biocommons/bioutils/issues

.. |license| image:: https://img.shields.io/github/license/biocommons/bioutils.svg
  :target: https://github.com/biocommons/bioutils/blob/master/LICENSE

.. |pypi_badge| image:: https://img.shields.io/pypi/v/bioutils.svg
  :target: https://pypi.org/project/bioutils/

	   
.. |cc_badge| image:: https://api.codeclimate.com/v1/badges/3a99e06ad0a842174b0a/maintainability
   :target: https://codeclimate.com/github/biocommons/bioutils/maintainability
   :alt: Maintainability

.. |cov_badge| image:: https://coveralls.io/repos/github/biocommons/bioutils/badge.svg?branch=master
   :target: https://coveralls.io/github/biocommons/bioutils?branch=master

