#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO


def convertidor_lista_multifasta(archivo_genbank, lista_locus, query, name_query):
    """
    DEFINICION: Función que basándose en un archivo genbank y en la lista de IDs 
                proteicos obtenida tras haber empleado el modulo filtrador_blast.py, 
                genera un archivo multifasta con las hits filtradas del blast realizado
    
    ARGUMENTOS:
        - archivo_genbank = archivo genbank que se quiere parsear 
        - lista_locus= Lista que contiene los nombres de los hits filtrados tras el blast
        - query= archivo multfasta que contiene todas las secuencias query sobre las que se ha hecho el blast
        - name_query= nombre de la secuencia query que estemos analizando
    
    RESULTADO: 
        - Crea un archivo llamado hits_blast_multifasta_filtrados que contiene
          todas las secuencias proteicas obtenidas al filtrar el blast
    
    **** Los archivos se escriben entre comillas dentro de los parentesis de la funcion ****
    
    """

#Creación del archivo premuscle_multifasta con las secuencias obtenidas tras filtrar el blast y que se usara para hacer el muscle
    premuscle = open(name_query+"_multifasta_filtrado","w", encoding="utf8")   

#Apertura del archivo query que queremos parsear 
    input_query=open(query, "r", encoding="utf8")
     
 

#bucle para parsear el archivo inicial con todos los query y escribir en el archivo output 
#el nombre y la secuencia de del query actual del que vamos a hacer muscle y prosite
    
    for record in SeqIO.parse(input_query, "fasta"):
        if (record.id) == name_query:
            premuscle.write(str(">"+record.id+"\n"+record.seq+"\n"+"\n"))
            
# Cerramos el archivo query inicial
    input_query.close()

#Apertura del archivo genbank que queremos parsear 
    input_gbk = open(archivo_genbank, "r", encoding="utf8")
    
 # Proceso de parseado del archivo genbanck
    for record in SeqIO.parse(input_gbk, "genbank"):
        for feature in record.features:
            if feature.type == 'CDS':
                locus = feature.qualifiers['locus_tag'][0]
                if (locus in lista_locus) and (locus != name_query):
 # si la secuencia esta en el query y en el subject, con esta sentencia evitamos incluirlo dos veces en el archivo a insertar en muscle
                    translation=feature.qualifiers['translation'][0]
                    seq_fasta=str(">"+locus+"\n"+translation +"\n"+"\n")
                    premuscle.write(seq_fasta)
                            
    premuscle.close()
    input_gbk.close()
    return