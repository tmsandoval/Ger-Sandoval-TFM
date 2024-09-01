Descarga de las lecturas de la base de datos seleccionada
```
# Crear un nuevo ambiente en Conda
Conda create –n tesis
# Activar el nuevo ambiente

Conda activate tesis
```
Se descarga el fichero de transcriptoma de referencia
```
curl https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o homocdna.fa.gz  salmon index -t homocdna.fa.gz  -i homo_index
```

Se crea un archivo ejecutable con el nombre "download"
```
nano download.sh
```
  El codigo de este archivo ejecutable se detalla:
```
#!/bin/bash
# Crear el directorio de datos y navegar a él
mkdir -p data
cd data
#Se crea un objeto con los cuatro ultimos numeros de las corridas de la base de datos que se usara.
seq=(1128 1137 1151 1164 1185 1176 1194 1201 1213 1217 1229 1234 1241 1250 1259 1265 1276 1298 1303 1310 1322 1328 1338 1346 1358 1372 1388 1393 1404 1411 1414 1415 1416 1417 1418 1419 1420 1421)

 # Bucle para crear directorios y descargar archivos
for i in "${seq[@]}";
 do
  mkdir -p "ERR360${i}"
  cd "ERR360${i}"
  # Descargar los archivos fastq.gz
fasterq-dump -p --split-files ERR360${i};
bgzip ERR360${i}_1.fastq;
bgzip ERR360${i}_2.fastq;
 # Volver al directorio padre
  cd ..
done
# Volver al directorio original
cd ..

```

```
```
