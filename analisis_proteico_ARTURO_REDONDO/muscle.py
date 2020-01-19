#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import datetime

from subprocess import Popen, PIPE


def muscle(hits_blast,name_query):
    """
    DEFINICION: Función que empleando muscle, realiza un alineamiento y un arbol 
                filogenetico sobre las secuencias fasta que le pases.
                
    ARGUMENTOS:
        - hits_blast= archivo que contiene las secuencias fasta que vamos a emplear en el muscle
                        (output directo de la funcion convertidor_listagenes_multifasta)
        - name_query= String que representa el nombre del query que estamos analizando
    
    RESULTADO:
        - Crea 3 archivos, uno con los parametros de cada proceso (MUSCLE_README), uno con 
          el arbol filogenético (muscle_neighbor_tree) y otro con el alineamiento (muscle_alineamiento)
    
    **** Los archivos se escriben entre comillas dentro de los parentesis de la funcion ****
    
    """
#######################################################################################################
# ALINEAMIENTO SECUENCIAS CON MUSCLE

# COMANDO BASH-LINUX: muscle -in hits_blast -out muscle_alineamiento
    alineamiento = Popen(['muscle','-in',hits_blast,'-out',name_query+"_muscle_alineamiento"], stdout=PIPE, stderr=PIPE)
    proceso_alineamiento = alineamiento.stderr.read().decode("utf-8")
    
    alineamiento.stderr.close()
    alineamiento.stdout.close()
    
#  Datos del proceso de alineamiento, 
#           - se printean por pantalla 
#           - se guardan en un archivo README MUSCLE 
    
    # impresion por pantalla para llevar un registro de las acciones que se realizan
    print()
    print("ALINEANDO SECUENCIAS DE LA SECUENCIA QUERY: "+name_query+"........")
    
    # escritura del archivo readme con la informacion del proceso de alineamiento
    g=open(name_query+"_MUSCLE_README","w",encoding="utf8")
    g.write("\n--------------------------------------------------------------------\n")
    g.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    g.write("\nSecuencia query: "+name_query+"\n")
    g.write("\n PROCESO DE ALINEAMIENTO CON MUSCLE \n")
    g.write(proceso_alineamiento)
    g.write("\n--------------------------------------------------------------------\n")
    g.close()
    
    
       
#######################################################################################################
# CREACIÓN DE ÁRBOL FILOGENETICO USANDO NEIGHBOR-JOINING CON MUSCLE

# COMANDO BASH-LINUX: muscle -maketree -in muscle_alineamiento -out muscle_neighbor_tree -cluster neighborjoining
    maketree = Popen(['muscle','-maketree','-in',name_query+"_muscle_alineamiento",'-out',name_query+"_muscle_neighbor_tree",'-cluster','neighborjoining'], stdout=PIPE, stderr=PIPE)
    creacion_tree = maketree.stderr.read().decode("utf-8")
    
    maketree.stderr.close()
    maketree.stdout.close()

    
#  Datos del proceso de creación de árlbol filogenético, 
#           - se printean por pantalla 
#           - se guardan en un archivo README MUSCLE 
    
    # impresion por pantalla para llevar un registro de las acciones que se realizan
    print("REALIZANDO ÁRBOL FILOGENÉTICO CON NEIGHBOR-JOINING DE LA SECUENCIA QUERY: "+name_query+"........")
    print("\n--------------------------------------------------------------------\n")
    

    #escritura del archivo muscle readme para incluir la informacion sobre el proceso de construccion del arbol filogenetico
    j=open(name_query+"_MUSCLE_README","a",encoding="utf8")
    j.write("\n--------------------------------------------------------------------\n")
    j.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    j.write("\nSecuencia query: "+name_query+"\n")
    j.write("\n PROCESO DE NEIGHBOR-JOINING (árbol filogenético) CON MUSCLE \n")
    j.write(creacion_tree)
    j.write("\n--------------------------------------------------------------------\n")
    j.close()
        
        
        
    

