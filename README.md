# Gene expression analysis of JQ1 treatment

## Starting the analysis and pipeline
To start the analysis load the following app (note due the size of the container the cold start could take up to several minutes):
https://rna-pipeline-4tzn6567yq-ey.a.run.app/lab

## Infrastructure

The project goal was to quantify the mRNA expression changes of JQ1 treatment in triple-negative breast cancer cell line (SUM159). This analysis has 4 parts, which will be descrbed in more details:
  * i) summary notebook (Turbine-rna-seq.ipynb)
  * ii) data preprocessing notebook (bin/RNA_seq_data_preprocessing.ipynb)
  * iii)RNA quantfication from raw reads using SALMON and STAR (bin/RNA_seq_quantification.ipynb)
  * iv) Statistical analaysis and visualizatin with DESeq2. (bin/turbine_rna_deseq2.ipynb)


All the calculation was carried out in cloud enviremnt (microsoft azure and google cloud), where both tradition (IaaS) and serverless technologies (like cloud-run) were used. The input data is deposited in the cloud (links are in the notebooks). The software infrastrucutre (neccessery packages, tools, etc.) are provided via docker. All the codes, notebooks, docker files are deposited in a github repository. The notebooks and the underlying infrastucre is deployed to google cloud run. 


