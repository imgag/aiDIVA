# Installation Of Additional Tools
This document provides links and instructions to download and install the necessary third party tools that are used in aiDIVA.

The folder of each respective tool can be moved to another place after the installation, they should be self contained.

Just make sure to give the correct path to the tool in the configuration file.

## ngs-bits
Ngs-bits is used to annotate the VCF files.

```
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2025_03 && git submodule update --recursive --init
make build_3rdparty
make build_libs_release
make build_tools_release
```

## Variant Effect Predictor (VEP)
VEP is used for all annotations that cannot be done using ngs-bits

Make sure to specify the correct paths for the VEP installation and the VEP cache data directory in the YAML configuration file

```
# specify the installation directory to match your own setup
# You should use absolute paths
vep_install_dir=ensembl-vep-release-109.3/
vep_cpan_dir=$vep_install_dir/cpan/
vep_data_dir=ensembl-vep-109/

wget https://github.com/Ensembl/ensembl-vep/archive/release/109.3.tar.gz
mkdir -p $vep_install_dir
tar -C $vep_install_dir --strip-components=1 -xzf 109.3.tar.gz
rm 109.3.tar.gz

# Install dependencies
mkdir -p $vep_cpan_dir
cpanm -l $vep_cpan_dir -L $vep_cpan_dir Set::IntervalTree URI::Escape DB_File Carp::Assert JSON::XS PerlIO::gzip DBI

# Download VEP cache data
mkdir -p $vep_data_dir
cd $vep_data_dir
mkdir -p ftp
cd ftp

wget ftp://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh38.tar.gz

# install ensembl-vep
PERL5LIB=$vep_install_dir/Bio/:$vep_cpan_dir/lib/perl5/:$PERL5LIB
cd $vep_install_dir

perl INSTALL.pl --SPECIES homo_sapiens --ASSEMBLY GRCh38 --AUTO ac --NO_UPDATE --NO_BIOPERL --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --NO_TEST

mkdir -p $vep_data_dir/cache/Plugins
cd $vep_data_dir/cache/Plugins

# make sure to download the pugin from the 110 release or later (in earlier releases it is not present)
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/refs/heads/release/110/AlphaMissense.pm
```
