# Base Image
FROM continuumio/miniconda3:4.6.14
LABEL authors="Xiandong Meng" \
      description="Docker image for Hybrid assembly " \
      software.name="unicycler" \
      software.description="Unicycler (0.4.8) is an assembly pipeline for bacterial genomes." \
      software.website="https://github.com/rrwick/Unicycler" \
      container.category="tool"

# Install system-level programs for quast
#RUN apt-get update && apt-get install -y zlib1g-dev vim git

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install unicycler
RUN conda install --freeze-installed -c conda-forge -y awscli bc perl python=3.7 \
  && conda install --freeze-installed -c bioconda -y unicycler=0.4.8 bbmap=38.73 samtools fastqc filtlong quast bowtie2 minimap2 qualimap tgsgapcloser java-jdk perl-bioperl \
  && conda clean -afy

#RUN conda install -c bioconda -y unicycler bbmap samtools fastqc quast filtlong python=3

# install long reads filter
#RUN conda install -c bioconda -y filtlong

# Install time and memory usage tool (timem)
RUN conda install -c jrmadsen -y timemory

# Conda upgrade all packages
RUN conda update --all

# Increase pilon's memory to default_jvm_mem_opts = ['-Xms16g', '-Xmx32g']
COPY pilon /opt/conda/share/pilon-1.23-1/pilon
COPY pilon /opt/conda/share/pilon-1.23-2/pilon

ADD run_unicycler_hybrid_pacbio.sh /usr/local/bin/run_unicycler_hybrid_pacbio.sh
RUN chmod +x /usr/local/bin/run_unicycler_hybrid_pacbio.sh

ADD run_unicycler_short_only.sh /usr/local/bin/run_unicycler_short_only.sh
RUN chmod +x /usr/local/bin/run_unicycler_short_only.sh

ADD run_unicycler_long_only.sh /usr/local/bin/run_unicycler_long_only.sh
RUN chmod +x /usr/local/bin/run_unicycler_long_only.sh

ADD run_unicycler_existing_long.sh  /usr/local/bin/run_unicycler_existing_long.sh
RUN chmod +x /usr/local/bin/run_unicycler_existing_long.sh

ADD run_unicycler_hybrid_nanopore.sh /usr/local/bin/run_unicycler_hybrid_nanopore.sh
RUN chmod +x /usr/local/bin/run_unicycler_hybrid_nanopore.sh

ADD template.yaml /usr/local/bin/template.yaml
RUN chmod +x /usr/local/bin/template.yaml

ADD circ-permute.pl  /usr/local/bin/circ-permute.pl
RUN chmod +x /usr/local/bin/circ-permute.pl

RUN mkdir -p /home/gitdir && \
           cd /home/gitdir && \
           git clone https://github.com/shingocat/lrscaf.git
# /home/gitdir/lrscaf/target/LRScaf-1.1.9.jar


RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
