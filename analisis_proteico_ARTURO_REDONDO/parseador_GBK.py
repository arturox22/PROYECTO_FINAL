# encoding: utf-8
#!/usr/bin/env bash

from Bio import SeqIO

def parseador(archivo_genbank):
    """
    
    DEFINICION: Función que parsea un archivo en formato genbank y genera un archivo multifasta
                con todas las secuencias proteicas que encuentre
                
    ARGUMENTOS:
        - archivo_genbank= archivo genbank del cual se van a extraer las secuencias en formato fasta
    
    RESULTADO:
        - Crea el archivo "subject_multifasta_gbk" que contiene todas las secuencias
          proteicas en formato fasta encontradas en el archivo genbank
    
    **** Los archivos se escriben entre comillas dentro de los parentesis de la funcion ****
    
    """
    
#Creación del archivo subject_multifasta_gbk --> output que será usado como subject en blaster para realizar el blast 
    multi_fasta = open("subject_multifasta_gbk","w", encoding="utf8")
    
#Apertura del archivo genbank que queremos parsear para su lectura
    input_gbk = open(archivo_genbank, "r", encoding="utf8")
    
# Proceso de parseado
    for record in SeqIO.parse(input_gbk, "genbank"):
        for feature in record.features:
            if feature.type == 'CDS':
                locus = feature.qualifiers['locus_tag'][0]
                translation=feature.qualifiers['translation'][0]
                seq_fasta=str(">"+locus+"\n"+translation +"\n"+"\n")
                multi_fasta.write(seq_fasta)

    multi_fasta.close()
    input_gbk.close()
    return

