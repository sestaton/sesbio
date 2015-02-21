sesbio
======

Bioinformatics tools for genome analysis

**ABOUT**

This is a collection of (mostly) Perl scripts for working with genomics data. There is no install or build process, just use them in place. See below for usage information.

**INSTALLATION**

To get the scripts, simply download the package:

    curl -L https://api.github.com/repos/sestaton/sesbio/tarball > sesbio.tar.gz

Or, use git:

    git clone https://github.com/sestaton/sesbio.git

Then install the Perl dependencies:

    cd sesbio
    curl -L cpanmin.us | perl - --installdeps .

After that, you can use the scripts in place or move them to where is most convenient. Run `git pull` in the base directory to keep things up to date (or run the curl command above to download the latest code). Send me a message or file an issue if you have feature requests or run into any issues.

**USAGE**

Most scripts have at least minimal documentation that can be accessed by typing the name of the script and reading the menu. Some scripts are well tested and have formal documentation that can be accessed by specifying the manual option (e.g., `perl script.pl --man`). Though, a lot of the scripts were written for one time use so not a lot of work was put into documenting or polishing the code. If something looks useful but is unclear, just ask, I should be able to help.

**LICENSE AND COPYRIGHT**

Copyright (C) 2013-2015 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: http://www.opensource.org/licenses/mit-license.php
