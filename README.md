sesbio
======

Bioinformatics scripts for genome analysis

**ABOUT**

This is a collection of scripts for working with genomics data and visualizing results. There is no install or build process, just use them in place. See below for usage information.

The directories in this repo contain scripts devoted to routine tasks in specific research domains. For example, gene annotation scripts are in the `gene_annotation` directory. Each directory contains sub-directories for shell and R scripts. The usage of these scripts should be obvious by the name or documentation, but see below.

**INSTALLATION**

To get the scripts, simply download the package:

    curl -sL https://api.github.com/repos/sestaton/sesbio/tarball > sesbio.tar.gz

Or, use git:

    git clone https://github.com/sestaton/sesbio.git

Then install the Perl dependencies:

    cd sesbio
    curl -sL cpanmin.us | perl - --installdeps .

After that, you can use the scripts in place or move them to where is most convenient. Run `git pull` in the base directory to keep things up to date (or run the curl command above to download the latest code). Send me a message or file an issue if you have feature requests or run into any issues.

**USAGE**

Most scripts have at least minimal documentation that can be accessed by typing the name of the script and reading the menu. Some scripts are well tested and have formal documentation that can be accessed by specifying the manual option (e.g., `perl script.pl --man`). Though, a lot of the scripts were written for one time use so not a lot of work was put into documenting or polishing the code. In particular, each topic directory contains a subdirectory of shell scripts and a subdirectory of R scripts. These scripts serve as a guide for repeating analyses or running complicated programs without relearning all the options, though some of the scripts are designed to be used with any data set. If something looks useful but is unclear, just ask, I should be able to help.

**LICENSE AND COPYRIGHT**

Copyright (C) 2013-2020 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: http://www.opensource.org/licenses/mit-license.php
