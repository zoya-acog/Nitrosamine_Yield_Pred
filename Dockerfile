FROM dockerhub.aganitha.ai:4443/atk/atk-conda:3.11
RUN pip install numpy pandas scikit-learn matplotlib seaborn torch torchvision 
RUN pip install torch-geometric
RUN pip install lightning
RUN conda install -c conda-forge mamba -y
RUN mamba install -c conda-forge rdkit -y