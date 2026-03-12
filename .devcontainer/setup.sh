#!/bin/bash

# Update and install system dependencies
sudo apt-get update
sudo apt-get install -y samtools perl wget tabix libpq-dev awscli

pip install psycopg2-binary hgvs pysam biopython pandas natsort pyarrow gffutils matplotlib xlsxwriter

mkdir -p data

if [ ! -f data/Homo_sapiens_assembly19.fasta ]; then
    echo "Downloading hg19... this may take a moment."
    
    aws s3 cp s3://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta data/ --no-sign-request
    aws s3 cp s3://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai data/ --no-sign-request
    aws s3 cp s3://broad-references/hg19/v0/Homo_sapiens_assembly19.dict data/ --no-sign-request    
fi

# Update Nextflow
nextflow self-update
nextflow -version

# Annovar
chmod +x annovar/*.pl
echo "Before you can use annovar you must download the annovar database files"
perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGeneWithVer annovar/humandb/
# perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGeneWithVerMRNA annovar/humandb/

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]${PWD#/workspaces/} ->\[\e[0m\] '

NEXTFLOW_DIR="/workspaces/$(basename $CONTAINER_WORKSPACE_FOLDER)/nextflow"

cat << 'EOF' >> $HOME/.bashrc
cleanup() {
    cd $NEXTFLOW_DIR
    nextflow clean -f 
    rm -f $NEXTFLOW_DIR/.nextflow.log*
    rm -f $NEXTFLOW_DIR/results/*    
    echo "Removed nextflow logs and results"
}
EOF

# Untested: 
# Download and install SnpEff
# cd /tmp
# wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
# unzip snpEff_latest_core.zip
# mv snpEff "/workspaces/$(basename $CONTAINER_WORKSPACE_FOLDER)/
# cd /workspaces/$(basename $CONTAINER_WORKSPACE_FOLDER)/
# java -jar snpEff.jar download GRCh37.p13

# Install VEP's dependencies. 
# On mac you should install htslib using brew: brew install htslib
# perlbrew install-cpanm
# cpanm DBI
# cpanm DBD::mysql@4.050
# cpanm Test::Differences Test::Exception Test::Perl::Critic Archive::Zip PadWalker Error Devel::Cycle Role::Tiny::With Module::Build LWP List::MoreUtils
# ~~cpanm Bio::DB::HTS~~
# perl INSTALL.pl --AUTO a --NO_HTSLIB

echo "You still need to do the following:"
echo "- Install VEP"
echo "- Install SnpEff"
echo "- Run 'conda install openssl certifi openpyxl'"
