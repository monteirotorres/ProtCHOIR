---
description: Readme Imported from GitHub
---

# Quick Instructions

[![https://raw.githubusercontent.com/monteirotorres/ProtCHOIR/master/ProtCHOIR/Contents/ProtCHOIR.svg?sanitize=true](https://raw.githubusercontent.com/monteirotorres/ProtCHOIR/master/ProtCHOIR/Contents/ProtCHOIR.svg?sanitize=true)](https://raw.githubusercontent.com/monteirotorres/ProtCHOIR/master/ProtCHOIR/Contents/ProtCHOIR.svg?sanitize=true) [![](https://camo.githubusercontent.com/6de4f94782584e14ca1f822b097a82daf81a32a0d00a237cc13bcca526bddabe/68747470733a2f2f7a656e6f646f2e6f72672f62616467652f3230353337323936322e737667)](https://zenodo.org/badge/latestdoi/205372962)

## ProtCHOIR

This pipeline was devised to create homo-oligomeric structures based on selected subsets of the PDB databank.

With ProtCHOIR you can supply either a sequence in FASTA format or a protomeric structure in the PDB format to obtain homo-oligomeric models based on homologues.

### Prerequisites

The following packages and external programs are used by ProtCHOIR scripts and must be installed and in either the binaries path or python path.

#### Python packages

> * progressbar2
> * pandas
> * biopython
> * pathlib
> * parasail
> * networkx
> * jinja2
> * numpy
> * matplotlib

#### External software \(must be installed separately\)

> * [PyMOL](https://sourceforge.net/projects/pymol/)
> * [parasail](https://github.com/jeffdaily/parasail)
> * [PSI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
> * [MAFFT](https://mafft.cbrc.jp/alignment/software/)
> * [PISA](http://www.ccp4.ac.uk/)
> * [GESAMT](http://www.ccp4.ac.uk/)
> * [Molprobity](http://www.ccp4.ac.uk/)
> * [Modeller](https://salilab.org/modeller/)
> * [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm)

Note: PISA, GESAMT and MolProbity may be installed as part of the [CCP4 Software Suite](http://www.ccp4.ac.uk/)

### Installation

The scripts are available as a [PyPi project](https://pypi.org/project/ProtCHOIR/). Just install them with:

`pip install ProtCHOIR`

### Initial Setup

If that is the first time you are running ProtCHOIR and you do not provide a configuration file \(with --conf\), the program will ask whether you desire the configuration file to be created. This configuration file simply has the paths to all the external software that are necessary.

The file also contains the path to a locally generated database \(referred to as "choirdb"\) in which it will look for possible homo-oligomeric proteins to serve as templates for modelling.

Make sure that the directory to which the choirdb variable is pointing actually exists.

The choirdb must be created locally and is a lengthy process whose total duration will depend on the processing capabilities of your machine. In the process, the whole pdb database will be downloaded, analysed and sorted in the expected directories.

Initial creation of the local database can be done with:

`ProtCHOIR -v -u --conf conf_file`

Subsequent updates will not re-download and re-analyse the whole PDB database, but only the new \(or updated\) entries.

### Usage

After the initial database set-up, you may run the program normally via command line, by invoking the ProtCHOIR executable and providing an input file either in PDB or FASTA format.

`ProtCHOIR -v -f protomer.pdb --conf conf_file`

If no conf file is yet in place, ProtCHOIR will ask you whether you want a default one to be generated, just run:

`ProtCHOIR`

And then modify the generated configuration file to match your environment.

To generate a full html report with detailed model analysis as output, run the program with:

`ProtCHOIR -v -f protomer.pdb --generate-report --conf conf_file`

To expose all available runtime options, run:

`ProtCHOIR -h`

### Methodology Flowchart

The image below summarizes the approach used by ProtCHOIR to build the homo-oligomeric proteins.

![https://raw.githubusercontent.com/monteirotorres/ProtCHOIR/master/ProtCHOIR/Contents/ProtCHOIRScheme.svg?sanitize=true](https://raw.githubusercontent.com/monteirotorres/ProtCHOIR/master/ProtCHOIR/Contents/ProtCHOIRScheme.svg?sanitize=true)

### Authors

Pedro Torres, Ph.D; Tom Blundell, FRS, FMedSci.

Department Of Biochemistry University of Cambridge 80 Tennis Court Road Cambridge CB2 1GA

### License

This project is licensed under Creative Commons license \([CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)\), provided along with the package - see [LICENSE](https://github.com/monteirotorres/ProtCHOIR/blob/master/LICENSE.txt). [![](https://camo.githubusercontent.com/34a22a9c8d4ae321a9f740e8a2c2afe2a8ba843db67daea709a8b2ff84d8fb38/68747470733a2f2f6d6972726f72732e6372656174697665636f6d6d6f6e732e6f72672f70726573736b69742f627574746f6e732f38387833312f7376672f62792e737667)](https://creativecommons.org/licenses/by/4.0/)

