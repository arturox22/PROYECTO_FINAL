#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def filtrado_blast(blast_no_filtrado, id_cutoff, cov_cutoff,evalue, name_query):
    """
    DEFINICION: Función que filtra los resultados obtenidos en el blast segun 
                los valores que el usuario indique, generando una lista con los
                ID de las secuencias proteicas filtradas y un archivo con el
                mismo formato que el output de blast, pero filtrado.
    
    ARGUMENTOS:
        - blast_no_filtrado= archivo tsv resultante de hacer un blast con formato de ouptut 6 (columnas separadas por tabulador)
              FORMATO DEL ARCHIVO:  4 columnas cuyo contenido y posicion son los siguientes
          
                                 sseqid     id      cov     evalue
    
        - id_cutoff= float que indica que el porcentaje de identidad por el que se quiere filtrar
        - cov_cutoff= float que indica que el porcentaje de coverage o cobertura por el que se quiere filtrar
        - evalue= float que indica el evalue por el que se filtra
    
    RESULTADO:
        - Devuelve una lista (lista_locus) con los nombres de los genes que han sido filtrados tras hacer el blast
        - Crea el archivo blast_filtrado con formato .tsv que contiene el output del filtrado del blast
    
    **** Los archivos se escriben entre comillas dentro de los parentesis de la funcion ****
    
    """
    
    #abrimos el archivo resultado del blast
    f=open(blast_no_filtrado,"r",encoding="utf8")
    
    # creamos un archivo que llevará los parametros del filtrado
    h=open("FILTRO_BLAST_README","w",encoding="utf8")
    h.write("\nPARÁMETROS DEL FILTRADO DE BLAST:\n")
    h.write("\nSe han seleccionado secuencias proteicas con un valor de Porcentade de Identidad:\n")
    h.write("\n %identidada >= " + str(id_cutoff))
    h.write("\nSe han seleccionado secuencias proteicas con un valor de Coverage:\n")
    h.write("\n coverage >= " + str(cov_cutoff))
    h.write("\nSe han seleccionado secuencias proteicas con un valor de E-Value:\n")
    h.write("\n evalue <= " + str(evalue))
    h.close()
    
    #creamos un archivo que llevará los resultados del blast filtrado
    g=open(name_query+"_blast_filtrado","w",encoding="utf8")
    g.write("Query_ID\tSecuencia_ID\tPorcentaje_identidad\tCoverage\tevalue\n")
    
    #creamos una lista que contendrá el nombre de las proteinas filtrados del blast
    #esta lista se usara en el script convertidor_lista_multifasta
    lista_locus=[]
    
    # rellenamos la lista y 
    for linea in f.readlines():
        a=linea.split(sep='\t')   # convertimos cada linea (con los campos separados por tabulador) en un string de campos llamado "a"
        
        # realizamos el filtrado sobre la linea
        try:
            #este try-except se realiza para poder guardar la primera linea de la tabla sin que de error
            if float(a[2])>=id_cutoff and float(a[3])>=cov_cutoff and float(a[4])<=evalue and name_query == a[0]:
                g.write(linea)
                lista_locus.append(a[1])
        except:
            continue
    f.close()
    g.close()
    
    return lista_locus

