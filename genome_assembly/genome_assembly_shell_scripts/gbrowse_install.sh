#!/bin/bash

function usage() {
cat << EOF
USAGE: bash $0 

Install GBrowse on a Linux cloud-based server

EOF
}

function print_error() {
cat << ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -gt 0 ]; then
    print_error
    usage
    exit 1
fi

function timer() {
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $ddh $dm $ds
    fi
}

function initlog () {
    LOGFILE=$1
    echo '' > ${LOGFILE}
}

function log () {
    echo $* >> ${LOGFILE}
}

base=${0%.*}
date=$(date | sed 's/\ //g;s/\:/\_/g')
initlog ${base}_${date}.log

log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
log $0 executed on `date`
log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
log `printf "\n"`

## get version

## Install webserver, db, and libs
sudo apt-get install apache2
sudo apt-get install mysql-server mysql-client
sudo mysql_install_db
sudo /usr/bin/mysql_secure_installation
sudo apt-get install libdb-dev graphviz libgd2-xpm-dev libxml2-dev zlib1g-dev

## Setup Perl environment
\curl -L http://install.perlbrew.pl | bash
echo "source ~/perl5/perlbrew/etc/bashrc" >> ~/.bashrc
source ~/.bashrc

# Install Perl
perlbrew install perl-5.20.1 -Dusethreads
perlbrew switch perl-5.20.1
perlbrew install-cpanm
perl -V

## Install Bioperl deps
cpanm Algorithm::Munkres Array::Compare Clone Convert::Binary::C DB_File Error GD Graph \
GraphViz HTML::Entities HTML::HeadParser HTML::TableExtract HTTP::Request::Common \
LWP::UserAgent List::MoreUtils Math::Random PostScript::TextBlock SOAP::Lite SVG \
SVG::Graph Set::Scalar Sort::Naturally Spreadsheet::ParseExcel Storable XML::Parser \
XML::Parser::PerlSAX XML::SAX XML::SAX::Writer XML::Simple XML::Twig XML::Writer YAML

## Install Bioperl
cpanm -n BioPerl

## Install GBrowse deps
cpanm Bio::Graphics CGI::Session Digest::MD5 ExtUtils::CBuilder File::Temp \
Text::ParseWords IO::String JSON LWP Statistics::Descriptive Time::HiRes \
Digest::SHA Date::Parse Term::ReadKey parent

## Install SAMtools AND Bio-SAMTools
wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools- 0.1.18.tar.bz2
tar xjf samtools-0.1.18.tar.bz2 && cd samtools-0.1.18
make CFLAGS=-fPIC
export SAMTOOLS=`pwd`
cd ..

wget https://cpan.metacpan.org/authors/id/L/LD/LDS/Bio-SamTools-1.39.tar.gz
tar xzf Bio-SamTools-1.39.tar.gz && cd Bio-SamTools-1.39
perl Build.PL && ./Build test && ./Build install

##\\\\\\\\\\\\\\\\ Install Kent source /////////////////////////////////////
sudo apt-get install mysql-server-5.0 apache2 libmysqlclient15-dev libpng12-dev libssl-dev openssl mysql-client-5.5  mysql-client-core-5.5
# set variables for compilation
export MACHTYPE=$(uname -m)
export MYSQLINC=`mysql_config --include`
export MYSQLLIBS=`mysql_config --libs`

DIRS='SCRIPTS=/usr/local/bin CGI_BIN=/usr/lib/cgi-bin DOCUMENTROOT=/var/www/genome BINDIR=/usr/local/bin'
# this does not seem to be necessary anymore
#ENCODE_PIPELINE_BIN=/usr/local/bin'

# download
cd /usr/local
wget http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
unzip jksrc.zip
mkdir -p /var/www/genome/

# compile libraries
cd kent/src/lib
make  

cd ../jkOwnLib
make 

# compile browser
cd ..
make $DIRS

# set permissions
cd /home/data/www
chown apache:apache -R *
cd
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

## Install GBrowse recs
cpanm Capture::Tiny Bio::Das Bio::DB::Sam Bio::DB::BigFile Crypt::SSLeay \
DB_File::Lock DBI DBD::mysql DBD::Pg DBD::SQLite Digest::SHA FCGI File::NFSLock \
GD::SVG Math::BigInt Net::OpenID::Consumer Net::SMTP::SSL Template \
VM::EC2 Parse::Apache::ServerStatus

## Install GBrowse
#git clone git@github.com:GMOD/GBrowse.git
#cd GBrowse
#perl Build.PL && ./Build && ./Build install
cpanm git@github.com:GMOD/GBrowse.git # cpanm Bio::Graphics::Browser
