FROM ubuntu:jammy  
LABEL maintainer="Balázs Ligeti <obbalasz@gmail.com>"

WORKDIR /bio-apps
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Budapest" LANG=C.UTF-8 LC_ALL=C.UTF-8 BINENV=/bio-progs/bin PYTHONVIRTENVNAME=tur_env DEFDIR=/bio-apps


RUN mkdir -p /bio-apps/turbine-rna-seq-ligeti && \
 export BINENV && \
 export PYTHONVIRTENVNAME && \
 export DEFDIR 
 
RUN apt-get update \
  && apt-get install -y wget unzip lbzip2

RUN apt -y install libcurl4-openssl-dev libxml2-dev mc libssl-dev libcairo2-dev libxt-dev libpq-dev  libxml2-dev libxslt1-dev libldap2-dev libsasl2-dev libffi-dev libgtextutils-dev
RUN apt -y install sudo curl gpg
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/r-project.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee -a /etc/apt/sources.list.d/r-project.list
RUN apt update
RUN apt -y install r-base

RUN R -e "install.packages('BiocManager')"


RUN R -e "BiocManager::install(c('gProfileR', 'treemap', 'AnnotationHub', 'ensembldb', 'gage', 'pathview', 'gageData', 'biomaRt', 'org.Hs.eg.db', 'DOSE', 'SPIA', 'tximport', 'DESeq2'), Ncpus = 20)"

#RUN R -e "BiocManager::install(c('gProfileR', 'treemap', 'AnnotationHub', 'ensembldb', 'gage', 'pathview', 'gageData', 'biomaRt', 'org.Hs.eg.db', 'DOSE', 'SPIA', 'tximport', 'DESeq2', 'pheatmap', 'DEGreport', 'ggrepel', 'IRkernel','Cairo','repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest', 'ggplot2', 'tidyverse', 'ape', 'seqinr', 'tidyverse', 'broom', 'patchwork', 'ggfortify', 'ComplexHeatmap', 'EnsDb.Hsapiens.v86', 'BSgenome.Hsapiens.UCSC.hg38', 'GSEABase', 'org.Hs.eg.db', 'AnnotationDbi', 'clusterProfiler', 'apeglm', 'ggnewscale'), Ncpus = 20)"

RUN R -e "BiocManager::install(c('pheatmap', 'DEGreport', 'ggrepel', 'IRkernel','Cairo','repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest', 'ggplot2', 'tidyverse', 'ape', 'seqinr', 'tidyverse', 'broom', 'patchwork', 'ggfortify', 'ComplexHeatmap', 'EnsDb.Hsapiens.v86', 'BSgenome.Hsapiens.UCSC.hg38', 'GSEABase', 'org.Hs.eg.db', 'AnnotationDbi', 'clusterProfiler', 'apeglm', 'ggnewscale'), Ncpus = 20)"
RUN R -e "devtools::install_github('stephenturner/annotables')"

RUN mkdir -p $BINENV
RUN cd $BINENV
RUN wget --directory-prefix=$BINENV https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
&& chmod a+x $BINENV/Miniconda3-latest-Linux-x86_64.sh\
&& /bin/bash $BINENV/Miniconda3-latest-Linux-x86_64.sh -b -p $BINENV/miniconda3 \
&& $BINENV/miniconda3/bin/conda config --add channels defaults \
&& $BINENV/miniconda3/bin/conda config --add channels bioconda\
&& $BINENV/miniconda3/bin/conda config --add channels conda-forge

ENV PATH=$BINENV/miniconda3/bin:$PATH 
RUN export PATH 

RUN conda install -y mamba -n base -c conda-forge

RUN $BINENV/miniconda3/bin/mamba install -y -c bioconda multiqc samtools star fastqc trimmomatic salmon
RUN $BINENV/miniconda3/bin/mamba install -y jupyterlab seaborn pandas jupyter scipy numpy biopython openpyxl xlsxwriter scikit-learn rpy2


SHELL ["/bin/bash", "--login", "-c"]

RUN cd /bio-apps/turbine-rna-seq-ligeti
#EXPOSE 9000
WORKDIR /bio-apps/turbine-rna-seq-ligeti

#CMD ["jupyter-lab"]
# CMD ["jupyter-lab", "--allow-root", "--port=9000", "--no-browser", "--NotebookApp.token='cica'"]
# 0.0.0.0,
# --ip 0.0.0.0