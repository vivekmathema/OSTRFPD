#!/bin/bash
set -e

# Use this script to automatically install and upgrade the necessary Python dependencies of OSTRFPD using pip.
# * Root mode: use "sudo -H ./install_dependencies" to install these dependencies for all users.
# * User mode: use "./install_dependencies --user" to install these dependencies in your local home directory.
# (Apart from these you may need additional non-Python packages, e.g. matplotlib, PyQt5 Designer which u may have to install manually)

p="pip install --upgrade $@"
$p biopython==1.72 PyQt5==5.9.1 future==0.16.0  gzip

# do we have pip3?
set +e
pip3 --version
pip3_rc=$?
set -e

if [[ $pip3_rc == 0 ]]; then
	# OK, we have pip3 as well
	p="pip3 install --upgrade $@"
	$p biopython==1.72 PyQt5==5.9.1 future==0.16.0  gzip
fi
