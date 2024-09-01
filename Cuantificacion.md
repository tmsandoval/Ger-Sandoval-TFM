Instalar Salmon en el ambiente de conda creado
```
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
$ conda create -n salmon salmon
```
Se indexa el transcriptoma de referencia. El índice es una estructura que Salmon utiliza para mapear de forma cuasi-mapeada las lecturas de secuenciación de ARN durante la cuantificación. 
```
$ salmon index -t homocdna.fa.gz -i homo_index
```
Se crea un fichero ejecutable. Dado que ejecutaremos el mismo comando en cada muestra, la forma más sencilla de automatizar este proceso es, nuevamente, un simple script de shell

```
nano cuantif.sh
```
El contenido del archivo ejecutable es el siguiente:
```
#!/bin/bash

# Definir la secuencia de IDs
seq=(1128 1137 1151 1164 1185 1176 1194 1201 1213 1217 1229 1234 1241 1250 1259 1265 1276 1298 1303 1310 1322 1328 1338 1346 1358 1372 1388 1393 1404 1411 1414 1415 1416 1417 1418 1419 1420 1421)


# Crear el directorio de salida para las cuantificaciones
mkdir -p quants

# Bucle para procesar los archivos
for i in "${seq[@]}";
 do
  dir="ERR360${i}"
  
  # Verificar si el directorio existe
  if [ -d "$dir" ]; then
    # Definir los nombres de archivo para fastq
    fastq1="${dir}/ERR360${i}_1.fastq.gz"
    fastq2="${dir}/ERR360${i}_2.fastq.gz"

    # Verificar si los archivos fastq existen
    if [ -f "$fastq1" ] && [ -f "$fastq2" ]; then
      samp="ERR360${i}"
      echo "Processing sample ${samp}"

      # Ejecutar Salmon para cuantificación
      salmon quant -i homo_index -l A \
           -1 "$fastq1" \
           -2 "$fastq2" \
           -p 8 --validateMappings -o "quants/${samp}_quant"
    else
      echo "Warning: Missing fastq files for ${samp}"
    fi
  else
    echo "Warning: Directory $dir does not exist"
  fi
done

```
Una vez creado el archivo se lo configura como ejecutable y ejecuta en el directorio de trabajo

```

chmod +x cuantif.sh

./cuantif.sh
```
