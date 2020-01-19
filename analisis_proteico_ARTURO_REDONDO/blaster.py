# encoding: utf-8
#!/usr/bin/env bash

from subprocess import Popen, PIPE

import sys

def blast(query, subject):
    
    """
    DEFINICION: Función que realiza un BLAST de la secuencia fasta pasada 
                como "query" sobre el conjunto de secuencias fasta denominado 
                "subject
    
    ARGUMENTOS:
        - query= archivo fasta con 1 sola secuencia
        - subject= archivo multifasta con una o mas secuencia
    
    RESULTADO:
        - Crea el archivo hits_blast_no_filtrados con formato .tsv que contiene
          el ouptuT del blast
    
    **** Los archivos se escriben entre comillas dentro de los parentesis de la funcion ****
    
    """
#######################################################################################################
# CREACION DE BD-->  para así aumentar la eficiencia del Blast

# COMANDO BASH-LINUX: makeblastdb -in subject -dbtype 'prot'  -out BD_subject   
    make_db = Popen(['makeblastdb','-in',subject,'-dbtype',"prot",'-out',"BD_subject"], stdout=PIPE, stderr=PIPE)
    error_making_BD = make_db.stderr.read().decode("utf-8")
    
    make_db.stderr.close()
    make_db.stdout.close()

# ESTE CONTROL DE ERRORES SIRVE PARA CUANDO EL MODULO BLASTER SE LANZA DE FORMA INDEPENDIENTE AL SCRIPT PRINCIPAL (main.py)
# CONTROL DE ERRORES EN LA CREACIÓN DE LA BD
    if error_making_BD: 
        print("Se produjo el siguiente error al crear la BD:\n%s" % error_making_BD)
        
# sentencias de escritura en el archivo log de control --> usadas si se emplea blaster como un modulo independiente
        
#        log=open("log","a",encoding="utf8")
#        log.write('\n Se produjo el siguiente error al crear la BD:\n%s' % error_making_BD)
#        log.write("\n--------------------------------------------------------------------\n")
#        log.close()
        
        sys.exit()
        return
    
########################################################################################################
        
# REALIZACIÓN DEL BLAST

#COMANDO DE BASH-LINUX: blastp -query query -db BD_subject -outfmt "6 sseqid pident qcovs evalue"
    do_blast = Popen(['blastp','-query',query,'-db',"BD_subject",'-outfmt',"6 qseqid sseqid pident qcovs evalue"], stdout=PIPE, stderr=PIPE)
    error_en_blast = do_blast.stderr.read().decode("utf-8")
    blast_sin_filtrar = do_blast.stdout.read().decode("utf-8")
    
    do_blast.stderr.close()
    do_blast.stdout.close()





## CONTROL DE ERRORES EN EL BLAST
    if not error_en_blast: 
        print("Query_ID\tSecuencia_ID\tPorcentaje_identidad\tCoverage\tevalue\n")
        print (blast_sin_filtrar) 
    else: 
        print("Se produjo el siguiente error al realizar el blast:\n%s" % error_en_blast)
        
#escritura del archivo log de control de errores --> util cuando se usa blaster como un modulo independiente
#        log=open("log","a",encoding="utf8")
#        log.write('\n Se produjo el siguiente error al realizar el blast:\n%s' % error_en_blast)
#        log.write("\n--------------------------------------------------------------------\n")
#        log.close()
      
        sys.exit()
        return   
    
    
# Si llegado aqui no ha habido error, escribimos dentro del archivo blas_no_filtrado, el ouput del blast
    output = open("hits_blast_no_filtrados","w",encoding="utf8")
    output.write("Query_ID\tSecuencia_ID\tPorcentaje_identidad\tCoverage\tevalue\n")
    output.write(blast_sin_filtrar)
    


    return