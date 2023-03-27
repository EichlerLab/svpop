FROM continuumio/miniconda3

MAINTAINER William T Harvey

COPY ./docker_env.yml /

RUN conda env create -f /environment.yml && conda clean -a

ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/conda/envs/svpop/lib/

ENV PATH=/opt/conda/envs/svpop/bin/:${PATH}

RUN wget https://github.com/Benson-Genomics-Lab/TRF/archive/refs/tags/v4.09.1.tar.gz && tar xzvf v4.09.1.tar.gz && cd TRF-4.09.1 && mkdir build && cd build && ../configure && make && make install && cd /

RUN wget https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz && tar zxvf rmblast-2.13.0+-x64-linux.tar.gz

ENV PATH=/rmblast-2.13.0/bin:${PATH}

RUN wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.5.tar.gz && cp RepeatMasker-4.1.5.tar.gz /usr/local && cd /usr/local && tar zxvf RepeatMasker-4.1.5.tar.gz && cd RepeatMasker && perl ./configure -trf_prgm=/usr/local/bin/trf -rmblast_dir=`which rmblastn` && cd /

RUN apt-get update -y && apt-get install -y rsync

RUN rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ./ucsc_tools && chmod 755 ./ucsc_tools/*

RUN git clone --recursive https://github.com/EichlerLab/svpop.git 

RUN mkdir -p /svpop/config/ && echo "snakemake -s \${SNAKEFILE} -j \${JOB_COUNT} --nt --ri -k --jobname \"{rulename}.{jobid}\" -w 60 \"\$@\"" > /svpop/config/runlocal.sh

ENV PATH=/svpop/:/ucsc_tools/:${PATH}
ENV PYTHONPATH=/svpop/:/svpop/dep/:${PYTHONPATH}

RUN rm RepeatMasker-4.1.5.tar.gz && rm v4.09.1.tar.gz


RUN echo "source activate base" > ~/.bashrc
