Instalar Salmon en el ambiente de conda creado
```
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
$ conda create -n salmon salmon
```
Se indexa el transcriptoma de referencia
```
$ salmon index -t athal.fa.gz -i athal_index
```
