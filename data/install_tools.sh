#!/bin/bash

folder=`pwd`
tools=$folder/tools/
annotation_sources=$folder/annotation_resources/

vep_install_dir=$tools/ensembl-vep-release-100.3/
vep_cpan_dir=$vep_install_dir/cpan/
vep_data_dir=$annotation_resources/ensembl-vep-100/

# download ensembl-vep
cd $tools
wget https://github.com/Ensembl/ensembl-vep/archive/release/100.3.tar.gz
mkdir $vep_install_dir
tar -C $vep_install_dir --strip-components=1 -xzf 100.3.tar.gz
rm 100.3.tar.gz

#install dependencies
mkdir -p $vep_cpan_dir
cpanm -l $vep_cpan_dir -L $vep_cpan_dir Set::IntervalTree URI::Escape DB_File Carp::Assert JSON::XS PerlIO::gzip DBI

#install BigWig support (needed to annotate phyloP)
cd $vep_install_dir
export KENT_SRC=$vep_install_dir/kent-335_base/src
export MACHTYPE=$(uname -m)
export CFLAGS="-fPIC"
wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
tar xzf v335_base.tar.gz
rm v335_base.tar.gz
cd $KENT_SRC/lib
echo 'CFLAGS="-fPIC"' > $KENT_SRC/inc/localEnvironment.mk
make clean && make
cd $KENT_SRC/jkOwnLib
make clean && make
cpanm -l $vep_cpan_dir -L $vep_cpan_dir Bio::DB::BigFile

#download VEP cache data
mkdir -p $vep_data_dir
cd $vep_data_dir
mkdir -p ftp
cd ftp
wget ftp://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_vep_100_GRCh37.tar.gz
wget ftp://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_refseq_vep_100_GRCh37.tar.gz

#install ensembl-vep
# TODO change assembly to GRCh38
PERL5LIB=$vep_install_dir/Bio/:$vep_cpan_dir/lib/perl5/:$PERL5LIB
cd $vep_install_dir
perl INSTALL.pl --SPECIES homo_sapiens --ASSEMBLY GRCh37 --AUTO acp --PLUGINS REVEL,CADD --NO_UPDATE --NO_BIOPERL --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --NO_TEST
cp $vep_data_dir/cache/Plugins/*.pm $vep_install_dir/modules/ #should not be necessary - probably a bug in the VEP installation script when using the CACHEDIR option (MS)



# Download and install ngs-bits
cd $folder
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2020_06 && git submodule update --recursive --init
make build_3rdparty
make build_tools_release
