# BactGenScripts

This repository contains a collection of scripts written by Simon Harris, formerly of the Bacterial Genomics team at the Wellcome Sanger Institute.

This repository is provided primarily for citation of this software.  We do not offer ongoing support or documentation for these scripts.

## Dependencies

- beast  1.8.4
- bwa  0.7.17
- gatk  3.7.0
- picard  1.126
- python 2.7
- raxml  8.2.8
- smalt  0.7.6
- samtools (bcftools & htslib) 1.2 _or_ 1.6
- ssaha2  2.5.5


## Install

The following will install on an Ubuntu Bionic machine;  the scripts require python 2.7, and some of the packages used are not available in the standard Ubuntu repos for more recent releases.

The following commands above have been tested on an `ubuntu:18.04` Docker image, but should in  principle work on any Bionic system; run the commands as root, setting  `INSTALL_DIR` to your preferred install path.

```bash
INSTALL_DIR='/opt/bact-gen-scripts'
# dependencies
apt-get install -y curl build-essential python-dev python-setuptools python-wheel python-2.7 python-pkg-resources python-tk
ln -s /usr/bin/python2.7 /usr/bin/python-2.7
curl -fsSL https://bootstrap.pypa.io/get-pip.py | python
# BactGenScripts download
mkdir -p ${INSTALL_DIR}
cd ${INSTALL_DIR}
curl -fsSL https://github.com/sanger-pathogens/bact-gen-scripts/archive/master.tar.gz | tar xzvf - --strip-components=1
# BactGenScripts install
pip install $(grep numpy ${INSTALL_DIR}/pip/requirements.txt)
pip install $(grep fisher ${INSTALL_DIR}/pip/requirements.txt)
pip install -r ${INSTALL_DIR}/pip/requirements.txt
```
