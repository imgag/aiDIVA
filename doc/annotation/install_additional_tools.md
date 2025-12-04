# Installation of Additional Tools
This document provides links and instructions to download and install the necessary third party tools that are used in aiDIVA.

The folder of each respective tool can be moved to another place after the installation, they should be self contained.

Just make sure to give the correct paths to the tools in the configuration file.

<!--
## ngs-bits
Ngs-bits is used to annotate the VCF files.

Please check that your system fulfills all the requirements needed for ngs-bits to prevent problems with the following commands. ([install ngs-bits on linux](https://github.com/imgag/ngs-bits/blob/master/doc/install_unix.md))

```
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2025_09 && git submodule update --recursive --init
make build_3rdparty
make build_libs_release
make build_tools_release
```
-->

## Variant Effect Predictor (VEP)
VEP is used for all annotations. Please check the official website for requirements that your system needs to fulfill before proceeding. VEP can be used in a containerized or locally installed version. Please head over to the official VEP documentation if you encounter any problems with VEP. ([VEP download and install](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html))

Make sure to specify the correct paths for the VEP installation and the VEP cache data directory in the YAML configuration file. Otherwise the local installed modules cannot be found. (Alternatively you can install the perl modules system wide with sudo)
For the containerized version you also need to specify the path to the folder, where all your annotation sources are lying. We need to bind that path to the container, otherwise VEP won't be able to read the data.

# Container version
Here you find instructions to create the Singularity image for later use. The benefits of the container version are that you don't need to set it up completely and everything is already included.
```
# download the Docker image and convert it to a singularity image named "vep_115-2.sif"
singularity pull --name vep_115-2.sif docker://ensemblorg/ensembl-vep:release_115.2

# you can use this image to install and prepare the cache files MAKE SURE TO SPECIFY THE CORRECT PATH
vep_data_dir=<full-path-to-folder>/vep_data
mkdir $vep_data_dir
singularity exec --bind $vep_data_dir:$vep_data_dir vep.sif INSTALL.pl -c $vep_data_dir -a cf -s homo_sapiens -y GRCh38
```

# Local installation
```
# specify the installation directory to match your own setup
# You should use absolute paths
vep_install_dir=ensembl-vep-release-115.2/
vep_cpan_dir=$vep_install_dir/cpan/
vep_data_dir=ensembl-vep-115/

wget https://github.com/Ensembl/ensembl-vep/archive/release/115.2.tar.gz
mkdir -p $vep_install_dir
tar -C $vep_install_dir --strip-components=1 -xzf 115.2.tar.gz
rm 115.2.tar.gz

# Install dependencies
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

# Download VEP cache data
mkdir -p $vep_data_dir
cd $vep_data_dir
mkdir -p ftp
cd ftp

wget ftp://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz

# install ensembl-vep
PERL5LIB=$vep_install_dir/Bio/:$vep_cpan_dir/lib/perl5/:$PERL5LIB
cd $vep_install_dir

perl INSTALL.pl --SPECIES homo_sapiens --ASSEMBLY GRCh38 --AUTO ac --NO_UPDATE --NO_BIOPERL --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --NO_TEST

mkdir -p $vep_data_dir/cache/Plugins
cd $vep_data_dir/cache/Plugins

# make sure to download the AlphaMissense plugin from the 110 release or later (in earlier releases it is not present)
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/refs/heads/release/115/AlphaMissense.pm
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/refs/heads/release/115/dbNSFP.pm
```
