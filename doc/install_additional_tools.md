# Installation Of Additional Tools
This document provides links and instructions to download and install the necessary third party tools that are used in AIdiva.

The folder of each respective tool can be moved to another place after the installation, they should be self contained.

Just make sure to give the correct path to the tool in the configuration file.

## ngs-bits
Ngs-bits is used to annotate the VCF files.

```
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout cba4aa891b5af683f74f0b0dabbe143719e0883a && git submodule update --recursive --init
make build_3rdparty
make build_tools_release
```

## Variant Effect Predictor (VEP)
VEP is used for all annotations that cannot be done using ngs-bits

Make sure to specify the correct paths for the VEP installation and the VEP cache data directory in the YAML configuration file

```
# specify the installation directory to match your own setup
vep_install_dir=ensembl-vep-release-103.1/
vep_cpan_dir=$vep_install_dir/cpan/
vep_data_dir=ensembl-vep-103/

wget https://github.com/Ensembl/ensembl-vep/archive/release/103.1.tar.gz
mkdir -p $vep_install_dir
tar -C $vep_install_dir --strip-components=1 -xzf 103.1.tar.gz
rm 103.1.tar.gz

# Install dependencies
mkdir -p $vep_cpan_dir
cpanm -l $vep_cpan_dir -L $vep_cpan_dir Set::IntervalTree URI::Escape DB_File Carp::Assert JSON::XS PerlIO::gzip DBI

# Install BigWig support (needed to annotate phyloP)
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

# Download VEP cache data
mkdir -p $vep_data_dir
cd $vep_data_dir
mkdir -p ftp
cd ftp
wget ftp://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_vep_103_GRCh37.tar.gz
#wget ftp://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_vep_103_GRCh38.tar.gz

# install ensembl-vep
PERL5LIB=$vep_install_dir/Bio/:$vep_cpan_dir/lib/perl5/:$PERL5LIB
cd $vep_install_dir
perl INSTALL.pl --SPECIES homo_sapiens --ASSEMBLY GRCh37 --AUTO acp --NO_UPDATE --NO_BIOPERL --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --NO_TEST
cp $vep_data_dir/cache/Plugins/*.pm $vep_install_dir/modules/ #should not be necessary - probably a bug in the VEP installation script when using the CACHEDIR option (MS)

```