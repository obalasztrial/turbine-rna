FROM europe-docker.pkg.dev/turbine-352714/turbine-repos/rna-base:latest 
LABEL maintainer="Balázs Ligeti <obbalasz@gmail.com>"

RUN mkdir /bio-apps/turbine-rna
WORKDIR /bio-apps/turbine-rna

RUN mkdir -p /bio-apps/turbine-rna/bin/ && \
 mkdir -p /bio-apps/turbine-rna/data/&& \
 mkdir -p /bio-apps/turbine-rna/results/ && \
 mkdir -p /bio-apps/turbine-rna/assets/

COPY bin /bio-apps/turbine-rna/bin
COPY data /bio-apps/turbine-rna/data
COPY results /bio-apps/turbine-rna/results
COPY assets /bio-apps/turbine-rna/assets


EXPOSE 8080
CMD ["jupyter-lab", "--allow-root", "--port=8080", "--no-browser", "--NotebookApp.token=''", "--ip", "0.0.0.0"]

