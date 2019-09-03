.. image:: ./ProtCHOIR/Contents/ProtCHOIR.svg
ProtCHOIR
############

This pipeline was devised to create homo-oligomeric structures based on selected subsets of the PDB databank.

With ProtCHOIR you can supply either a sequence in FASTA format or a protomeric structure in the PDB format to obtain homo-oligomeric models based on homologues.


Prerequisites
*************

The following packages and external programs are used by ProtCHOIR scripts and must be installed and in either the binaries path or python path.

Python packages
===============

  - progressbar
  - pandas
  - biopython
  - pathlib
  - parasail
  - modeller
  - networkx
  - jinja2
  - numpy
  - matplotlib
  - pickle
  - gzip
  - textwrap


External software (must be installed separately)
================================================

  - `PyMOL`_
  - `parasail`_
  - `PSI-BLAST`_
  - `MAFFT`_
  - `PISA`_
  - `GESAMT`_
  - `Molprobity`_


.. _`PyMol`: https://sourceforge.net/projects/pymol/
.. _`parasail`: https://github.com/jeffdaily/parasail
.. _`PSI-BLAST`: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _`MAFFT`: https://mafft.cbrc.jp/alignment/software/
.. _`PISA`: http://www.ccp4.ac.uk
.. _`GESAMT`: http://www.ccp4.ac.uk
.. _`Molprobity`: http://www.ccp4.ac.uk

Note: PISA, GESAMT and MolProbity may be installed as part of the `CCP4 Software Suite`_

.. _`CCP4`: http://www.ccp4.ac.uk


Installation
************
The scripts are available as a `PyPi project`_. Just install them with:

.. _`PyPi project`: https://pypi.org/project/ProtCHOIR/


:code:`pip install ProtCHOIR`


Initial Setup
*************
If that is the first time you are running ProtCHOIR and you do not provide a configuration file (with --conf), the program will ask whether you desire the configuration file to be created.
This configuration file simply has the paths to all the external software that are necessary.

The file also contains the path to a locally generated database (referred to as "choirdb") in which it will look for possible homo-oligomeric proteins to serve as templates for modelling.

The choirdb must be created locally and is a lengthy process whose total duration will depend on the processing capabilities of your machine. In the process, the whole pdb database will be downloaded, analysed and sorted in the expected directories.

Initial creation of the local database can be done with:

:code:`ProtCHOIR -v -u --conf conf_file`

Subsequent updates will not re-download and re-analyse the whole pd, but only the new (or updated) entries.

Usage
*****
After the initial database set-up, you may run the program normally via command line, by invoking the ProtCHOIR executable and providing an input file either in PDB or FASTA format.

:code:`ProtCHOIR -v -f protomer.pdb --conf conf_file`

Running:

:code:`ProtCHOIR -h`

Will expose all available runtime options.

Authors
*******
Pedro Torres, Ph.D;
Sony Malhotra, Ph.D;
Tom Blundell, FRS, FMedSci.

Department Of Biochemistry
University of Cambridge
80 Tennis Court Road
Cambridge CB2 1GA



License
*******

This project is licensed under Creative Commons license (CC-BY-4.0_), provided along with the package - see `LICENSE`_.

.. _LICENSE: https://github.com/monteirotorres/ProtCHOIR/blob/master/LICENSE.txt

.. _CC-BY-4.0: https://creativecommons.org/licenses/by/4.0/

.. image:: https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by.svg
  :target: https://creativecommons.org/licenses/by/4.0/
