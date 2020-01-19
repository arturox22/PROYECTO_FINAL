#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# IMPORTACIÓN DE MÓDULOS DE PYTHON

from datetime import datetime

import os

from subprocess import Popen, PIPE

import sys


# IMPORTACIÓN DE MÓDULOS DE TERCEROS

from Bio import SeqIO

from Bio.ExPASy import Prosite,Prodoc

import re


#IMPORTACIÓN DE MÓDULOS PROPIOS DEL SCRIPT


import blaster as b

import convertidor_listagenes_multifasta as c

import filtrador_blast as f

import mkdir_python as m

import muscle as musc

import ordena as o

import parseador_GBK as pars

import prosite as pro

import prodoc as doc



def main():
    """
    DEFINICION:
    
    Función principal que pone en conjunto la ejecucion de los siguientes
    modulos, controlando su progreso:
     
   -- parseador_GBK.py -- blaster.py -- filtrador_blast.py -- 
   -- convertidor_listagenes_multifasta.py -- muscle.py -- prosite.py
   
    """
    
    #imprimimos por pantalla la presentacion
    
    print("""
          
    
    BIENVENIDO !!!!!!!!!!!!!!!!!
        
        INFORMACIÓN GENERAL SOBRE EL SCRIPT:
        
        Este script esta diseñado para el análisis y tratamiento de secuencias 
        proteicas en formato fasta.
        
        En concreto realoizará las siguientea acciones:
            
            1. Parseado de un archivo genbank que se empleará como subject en el blast
            
            2. Realización de blast sobre un archivo query y el archivo subject obtenido
            del parseado anterior. El script admite que el query sea un achivo multifasta
            con mas de una secuencia.
            
            3. Filtro del blast anterior según los valores de %_identidad, coverage y evalue
            
            4. Alineamiento de las secuencias filtradas del blast y creación de un árbol
            filogenetico neighbor-joining, todo ello con el programa muscle
            
            5. Búsqueda de dominios proteicos (presentes en la base de datos de Prosite)
            en las secuencias filtradas obtenidas del blast anterior.
                * Se ofrece la posibilidad de obtener la informacion extra almacenada en la 
                base de datos prosite.doc por cada dominio detectado.
        
    
        La información acerca del proceso y los posibles errores durante la ejecución del script, 
        se recogen en un archivo llamado "log".


		       -------------------------------------------------------------------------
	RECOMENDACIÓN: |  Ejecute el script en una ventana completa para que el menú imprimido | 
		       |            por pantalla tenga mejor legibilidad		       |
		       -------------------------------------------------------------------------
    
	
	VERSION:        analisis_proteicos v1.0


####################################################################################################
    """)
        
    info_detallada="""

######################################################################################################


*** ESPECIFICACIONES DEL DISEÑO DEL SCRIPT ***


    - Script desarrollado en python en la version: Python 3.7.4 (default, Aug 13 2019, 20:35:49)  

    - FORMA DE USO--> ejecutar modulo main e introducir por teclado los parametros y archivos que
                      el script vaya solicitando en cada momento.

    - MÓDULOS QUE REQUIERE EL SCRIPT:

		  *****************************************************************
		  ***     Para más informacion sobre un módulo consultar la     ***	       
		  ***   informacion de uso especificada dentro de cada módulo   ***	
                  *****************************************************************


        * __init__.py

	* main.py --> módulo principal que llama a los diversos módulos y en el que se especifica 
		      el menú de interacción con el usuario

	* parseador_GBK.py --> contiene funciones para parsear un archivo genbank y extrae y 
			       almacena todas las secuencias proteicas en formato fasta en un archivo

	* blaster.py --> contiene funciones para realizar un blast de la secuencia query que se le 
			 especifique sobre el subject que se le pase

	* filtrador_blast.py y convertidor_listagenes_multifasta.py --> ambos módulos contienen
			funciones para el filtrado de los hits obtenidos tras hacer un blast.
			Primero se ejecuta filtrador_blast.py y sobre el resultado obtenido 
			se ejecuta convertidor_listagenes_multifasta.py

	*  muscle.py --> contiene funciones para realizar un alineamiento de secuencia y un árbol 
			 filogenético de las secuencias que se le pasen empleando muscle

	* prosite.py y prodoc.py --> contienen funciones para parsear la base de datos de prosite
			  en busca de dominios proteicos presentes en las secuencias que le pases

	* mkdir.py y ordena.py --> contienen funciones auxiliares muy útiles que son empleadas 
				   en el módulo main.py



   - ARCHIVOS NECESARIOS QUE REQUIERE EL SCRIPT Y QUE DEBEN ESTAR CONTENIDOS EN EL PAQUETE DONDE
     SE ENCUENTRAN LOS MÓDULOS:
	

	* Archivo genbank que se quiera parsear

	* Archivo query que se vaya a analizar

	* prosite.dat y prosite.txt --> archivos que constituyen la base de datos de PROSITE.
	  (  realmente prosite.txt se corresponde con el archivo prosite.doc que se 
	     descarga de la página oficial de prosite. Se guarda con extensión .txt
	           para evitar posibles problemas con el formato .doc                )



   - Todos los archivos generados se guardan en la carpeta ANALISIS_PROTEICO_FECHA que se genera
     al ejecutar el main.py dentro del propio paquete del script. Esta carpeta está divida en:
        
	* DATA --> carpeta que contiene los archivos input que se han analizado
 
	* RESULTS --> carpeta que contiene los archivos generados. A su vez se divide en:
		
		- BLAST --> contiene los archivos generados al usar blaster.py
		- MUSCLE --> contiene los archivos generados al usar muscle.py
		- PROSITE --> contiene los archivos generados al usar prosite.py
		- Archivo log --> archivo de control que contiene la información
				  asociada al proceso de ejecución del main.py



######################################################################################################
    
*** CONSIDERACIONES SOBRE EL DISEÑO DEL SCRIPT ***


    1. ANTES DE COMENZAR, ASEGURESE DE TENER LOS SIGUIENTES ARCHIVOS Y SOLO ELLOS EN LA MISMA 
       CARPETA/PAQUETE QUE LOS MODULOS DEL SCRIPT. Si ya ha ejecutado el script previamente, lo
       aconsejado por los desarrolladores es que la carpetas que se generan en el interior del 
       paquete se muevan a un directorio diferente.***
    
            - Archivo genbank --> será parseado y usado como subject en el blast
            
            - Query --> archivo fasta o multifasta que se usara como query en el blast
            
            - Archivos prosite .dat y .txt --> constituyen la base de datos de prosite
            y son necesarios para poder buscar dominios proteicos en las proteinas de 
            interes
            
    2. EL SCRIPT ESTA DISEÑADO PARA REALIZAR TODOS LOS PASOS (parseado de genbank, blast, filtrado 
       del blast, muscle y prosite), de forma, que una vez realizado todo, los archivos generados 
       son reorganizados en carpetas.
       
             - Estas carpetas se generan dentro del paquete que contiene los modulos del script. 
               Por tanto, si quiere correr dos veces seguidas el script, asegurese de que en dicho 
               paquete tiene los archivos definidos en 1 adecuados para la segunda ejecucion del
               script.
    
    3. Es posible que el archivo prosite.doc descargado de la página web de prosite de fallos de 
       formato en algunos sistemas. Dada esta situación, hemos corregido el archivos prosite.doc
       y lo hemos guardado con extension .txt para evitar complicaciones. 
       
    ****  RECOMENDACION: no cambie ni modifique los archivos prosite.dat y prosite.txt que van  ****
    ****  incluidas en el paquete del script. En caso de querer actualizarlos, debe guardar el  ****
    ****  archivo prosite.doc con la extension .txt, pues así está definido en el script.       ****
    
    4. El script esta diseñado de manera que los parámetros del filtrado del blast se eligen una 
       sola vez y se aplican sobre todas las secuencias query analizadas (no se puede filtrar cada
       secuencia con unos parametros distintos). Este diseño es así por motivos de eficiencia, 
       para aumentar la velocidad de ejecucion del script y poder automatizar el proceso. 
    
    ****  RECOMENDACION: En caso de querer analizar muchas query, pero filtrarlas con paramétros  ****
    **** distintos, deberá hacerlo manualmente, una por una, ejecutando el script con cada una.   ****
    
    5. La información que se obtiene de prosite.txt (=prosite.doc) no se almacena en ningun archivo,
       solo se imprime por pantalla. De querer guardarla, se recomienda copiar la informacion 
       manualmente y pegarla en un archivo de texto.
       
    6. Los módulos incluidos en el paquete de este script son plenamente funcionales de forma 
       independiente, no obstante, están diseñados especificamente para funcionar en conjunto, 
       por lo que seguramente no estén optimizados al máximo de lo que podrían estar.
       
    7. En el paquete se incluyen modulos auxiliares muy útiles como ordena.py o mkdir_python.py.
       (ordena.py emplea un algoritmo de tipo quicksort, diseñado por ARTURO REDONDO, en el 3ºCURSO
       de biotecnología (UPM) en la asignatura de PROGRAMACION del año 2018-19)
      
    8. La busqueda de dominios proteicos (modulo prosite.py) se realiza sobre todas las secuencias 
       proteicas que han sido hits en el blast.  Si el archivo que se introduce como query es 
       multifasta, la busqueda de dominios se realiza sobre cada uno de los hits obtenidos en cada
       una de las secuencias, de manera que por cada secuencia query obtienes un archivo dominios_
       proteicos.tsv. CUIDADO: El modulo prosite.py se ejecuta individualmente sobre cada una de las 
       secuencias query, y su proceso se imprime por pantalla; puede llegar a ser ilegible si estas 
       analizando muchas secuencias query o si el filtrado de blast realizado ha sido poco restrictivo.
      
    9. Los archivos dominios_proteicos.tsv, tienen un formato tsv, muy practico para ser procesado
       en el bash de linux, pero muy poco visual.
    
    10. A la hora de buscar informacion en prosite.txt, el script ofrece una lista que recopila todos
        los pdoc accession asociados a los dominios encontrados en TODAS las ejecuciones realizadas
        de prosite.py

    **** RECOMENDACION: Si ha ejecutado el script con un archivo multiquery, abra individualmente    ****
    **** en un editor de texto los archivos dominios.tsv generados por cada secuencia query para     ****
    **** visualizar claramente que dominios se han detectado en cada caso. Si entonces desea ampliar ****
    **** la descripcion de alguno de ellos, fijese en su PDOC_ACCESSION y escribalo por teclado para ****
    **** que el script le proporcione la información que busca                                       ****
    
             

######################################################################################################

ARCHIVOS GENERADOS :
    
    
- Al parsear el genbank se genera el archivo:
    
        subject_multifasta_gbk --> secuencias proteicas en fmt fasta del gbk
        
        
- Al realizar el blast se genera el archivo:
        
        blast_no_filtrado --> contiene el resultado del blast en fmt .tsv 
        

- Al filtrar el blast se generan 3 archivos:
    
        blast_filtrado --> contiene el resultado del blast filtrado en fmt .tsv
        Se genera un archivo por cada secuencia query analizada
    
        multifasta_filtrado --> secuencias proteicas filtradas del blast anterior 
        y en fmt fasta. Se genera un archivo por cada secuencia query analizada
        
        FILTRO_BLAST_README --> contiene los parametros de filtrado del blast
        
        
- Al emplear muscle se generan 3 archivos por cada secuencia query analizada:
        
        muscle_alineamiento --> secuencias fasta alineadas
        
        muscle_neighbor_tree --> arbol filogenetico de las proteinas anteriores
        
        MUSCLE_README --> archivo de texto con el proceso y los parametros
        empleados tanto para el alineamiento como para la creación del árbol 
        
        
- Al emplear prosite se genera el archivo por cada secuencia query analizada:
    
        dominios_proteicos.tsv --> archivo .tsv que muestra los dominios 
        que presenta cada secuencia proteica 
        
- ARCHIVO .log --> Archivo de control que guarda un registro del proceso y 
                   de los errores generados durante la ejecución del script
        
        
        
-------------------------------------------------------------------------------        
    
**** SCRIPT DESARROLLADO POR ARTURO P. REDONDO LUQUE, BIOTECNOLOGÍA (UPM) ****

-------------------------------------------------------------------------------
               
            """
    si_info=input("""
Desea saber más información detallada y consideraciones generales acerca del funcionamiento del script
                      (Si/No):    """)
    
    if si_info in ["SI","SÍ","Si","Sí","si","sí","YES","Yes","yes","Y","S","s"]:
        print(info_detallada)
    
            
    print("""    
          
    ========================================================================================
    
                                        INICIO DEL SCRIPT
    
    ========================================================================================      
    
    """)
    
    #creamos el archivo log --> control de errores y del proceso
    
    log=open("log","w",encoding="utf8")
    
    # creamos las siguientes carpetas donde se almacenaran los resultados
    
    dir_pcipal="ANALISIS_PROTEICO_"+datetime.now().strftime("%Y-%m-%d %H:%M")
    m.mkdir(dir_pcipal)
    
    data_dir=dir_pcipal+"/data"
    m.mkdir(data_dir)
    
    results_dir=dir_pcipal+"/results"
    m.mkdir(results_dir)
    
    # creacion de esta lista que sera empleada mas adelante
    lista_name_query=[]
    
    try:  # bucle de control del proceso en general --> captura excepcion cuando el usuario
          # detiene el script manualmente con ctrl + c
        
    
##############################   PARSEADO  GENBANK     ###############################
        
        while True: #bucle para controlar la ejecucion del parseador
            
            archivo_genbank=input("""
                        
                1.      PARSEADO DE UN ARCHIVO GENBANK
                
                
Por favor, introduca el archivo genbank a parsear y que será utilizado como
subject al realizar el blast:     
            
    ATENCIÓN: Introduzca el nombre del archivo (sin "" y sin darle al espacio 
    antes de escribir) y pulse enter:   
                
                                                    
                """)
            
            try: #control de errores del parseado
                
                #ejecucion del parseado
                pars.parseador(archivo_genbank)
                
                #impresion del proceso por pantalla
                print()
                print("-------------------------------------------------------------------------------------------------------------")
                print(" PARSEADO REALIZADO CON EXITO. Se ha generado el archivo subject_multifasta_gbk, guardado en la carpeta data")
                print("-------------------------------------------------------------------------------------------------------------")
                print()
                
                #escrtira del archivo log --> control del progreso del script
                log.write("\n-----------------------------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                log.write("\nPARSEADO REALIZADO CON EXITO. Se ha generado el archivo subject_multifasta_gbk, guardado en la carpeta data\n")
                log.write("-----------------------------------------------------------------------------------------------------------\n")

                # traslado de los archivos a sus correspondientes carpetas 
                m.copy("subject_multifasta_gbk",data_dir)
                m.copy(archivo_genbank,data_dir)
                
                #break para salir del bucle de control del parseador
                break
            
            except:
                # Control de posibles errores y escritura de ellos en el archivo log
                print()
                error_pars="Ha ocurrido un error. Por favor, asegurese de que el archivo introducido por teclado tiene formato genbank."
                print("----------------------------------------------------------------------------------------")
                print (error_pars)
                print("----------------------------------------------------------------------------------------")
                log.write("\n ------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.write(error_pars)
                log.write("\n------------------------------------------------------------------------------------\n")
                
                
            
#####################################   BLAST  ########################################
               
        while True: # Bucle de control para realizar el blast
            
            query=input("""
            
     2.      BLAST DE UN ARCHIVO QUERY SOBRE EL PARSEADO DEL GENBANK
                
                
        Por favor, introduca el archivo QUERY (en formato fasta) sobre el que se
        realizará el blast:     
            
            ATENCIÓN: Introduzca el nombre del archivo (sin "" y sin darle al espacio 
            antes de escribir) y pulse enter:
                            
                            
                                     """)
            
            try: # control de errores del blast
                
                # realizacion del blast
                b.blast(query,"subject_multifasta_gbk")
                
                #impresion del progreso por pantalla
                print()
                print("---------------------------------------------------------------------------------------------------")
                happy_blast="""
BLAST REALIZADO CON EXITO. Se han generado los siguientes archivos, guardados 
en la carpeta results/BLAST:
       
     - hits_blast_no_filtrados --> archivo con el ouput del blast en fmt .tsv 
     - BD_subject.phr , BD_subject.pin y BD_subject.psq --> BD para optimizar el algoritmo de blast
                    """
                print(happy_blast)
                print("---------------------------------------------------------------------------------------------------")
                
                #escritura del archivo log de control del proceso del script
                log.write("\n------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.write(happy_blast)
                log.write("\n------------------------------------------------------------------------------------\n")
                
                    
                # creacion de carpetas y traslado de archivos
                m.copy(query, data_dir)
                
                blast_dir=results_dir+'/BLAST'
                m.mkdir(blast_dir)
                m.copy("hits_blast_no_filtrados",blast_dir+'/hits_blast_no_filtrados') 
                
                BD_subject_dir=blast_dir+'/BD_subject'
                m.mkdir(BD_subject_dir)
                m.copy("BD_subject.phr",BD_subject_dir+'/BD_subject.phr')
                m.borrado("BD_subject.phr")
                m.copy("BD_subject.pin",BD_subject_dir+'/BD_subject.pin') 
                m.borrado("BD_subject.pin")
                m.copy("BD_subject.psq",BD_subject_dir+'/BD_subject.psq') 
                m.borrado("BD_subject.psq")
                m.borrado("subject_multifasta_gbk")
                 
                # break para salir del bucle de control del blast
                break
            
            except:
                # Control de posibles errores y escritura de ellos en el archivo log
                error_blast="Ha ocurrido un error al realizar el blast. Por favor, asegurese de que el archivo QUERY introducido por teclado tiene formato fasta."
                print("----------------------------------------------------------------------------------------")
                print(error_blast)
                print("----------------------------------------------------------------------------------------")
                log.write("\n------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.write(error_blast)
                log.write("\n------------------------------------------------------------------------------------\n")
               


###############################  CREACIÓN DE LISTA DE QUERYS  ###############################
                
# Creamos una lista con los nombres de las query --> nos servira para trabajar individualmente 
# con cada secuencia query por separado

        lista_name_query=[]
        for record in SeqIO.parse(query, "fasta"):
            lista_name_query.append(record.id)
           

#################################  FILTRADO DE BLAST  #######################################
        
        
        while True: # Bucle de control para realizar el filtrado del blast
            
            info_filtrado= """

******************************   ATENCIÓN  ***********************************                       
 
    
    A continuación se realizará el filtrado del blast. Si quiere trabajar con 
    todas las secuencias hits del Blast, es decir, si no quiere filtrar, 
    pongalos siguientes valores de porcentaje de identidad, coverage y evalue,
    para que el script pueda realizar los siguientes pasos :
        
        Porcentaje de identidad (%): 0
        
        Coverage: 0
        
        E-Value: 1000
        
 *** Se filtra por valores mayores o iguales que el % identidad y el coverage ***
 ***             elegido y  menores o iguales que el evalue                    ***
        
        
 

******************************************************************************"""
            print(info_filtrado)
            
            #variables necesarias para filtrar el blast pedidas para introducir 
            #teclado
            
            id_cutoff=input("""
            
     2.      FILTRADO DE BLAST:
                
                
Por favor, introduca los parámetros sobre los que fitrar (solo numeros):
        
                Porcentaje de identidad (%)
         (Secuencias con un valor mayor o igual):   """)
                
            cov_cutoff=input("""
                
                Coverage:   
         (Secuencias con un valor mayor o igual):   """)
            evalue=input("""
                
                E-value:   
         (Secuencias con un valor menor o igual):   """)
            
            
            
            try: # control de errores del filtrado del blast
                
     # REALIZACIÓN DEL FILTRADO DEL BLAST
     
          #creamos una lista de las listas con los locus. Cada elemento de esta lista, será una lista con los hits de cada secuencia query
                super_lista_locus=[]           
                for name_query in lista_name_query:
                    # Cada elemento de la lista sera la lista resultado del modulo filtrador_blast
                    super_lista_locus.append(f.filtrado_blast("hits_blast_no_filtrados",float(id_cutoff),float(cov_cutoff),float(evalue),name_query))
                
                # Sobre cada query aplicamos la funcion "convertidor" para asi conseguir todas las secuencias fasta resultantes tras el filtrado de blast
                for num_query in range(len(lista_name_query)):
                    c.convertidor_lista_multifasta(archivo_genbank, super_lista_locus[num_query], query, lista_name_query[num_query])
                    
                    #CONSIDERACIONES SOBRE ESTE ALGORITMO
                    # --> num_query es la posicion de la lista en la que estamos--> lista_name_query y super_lista_locus ambas son de la misma dimensión
                    # --> super_lista_locus[num_query] --> elementos de super_lista_locus --> cada elemento es una lista obtenida tras aplicar filtrador_blast.py 
                    # --> query --> archivo inicial con todas las secuencias query que se quieren analizar
                    # --> lista_name_query[num_query] --> nombre de la secuencia query del archivo inicial que se esta analizando
                    
                #impresion del progreso por pantalla
                print()
                print("---------------------------------------------------------------------------------------------------")
                buen_filtro="""
    FILTRADO DE BLAST REALIZADO CON EXITO. Se han generado los siguientes archivos, guardados 
    en la carpeta results/BLAST:
       
     - Archivos blast_filtrados --> output de blast con fmt .tsv filtrado
     - Archivos name_query_multifasta_filtrados --> secuencias fasta de las secuencias filtradas
     - FILTRO_BLAST_README --> contiene los parametros de filtrado del blast
    
    * Se ha generado el mismo numero de archivos que de querys analizados
    
                """
                print(buen_filtro)
                print("---------------------------------------------------------------------------------------------------")
                
                #escritura del archivo log de control del proceso del script
                log.write("\n------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.write(buen_filtro)
                log.write("\n------------------------------------------------------------------------------------\n")
                    
                # creacion de carpetas y traslado de archivos
                m.borrado("hits_blast_no_filtrados")
                
                for name_query in lista_name_query:
                    m.copy(name_query+"_multifasta_filtrado", blast_dir+'/'+name_query+'_multifasta_filtrado')
                    m.copy(name_query+"_blast_filtrado",blast_dir+'/'+name_query+'_blast_filtrado')
                    m.borrado(name_query+"_blast_filtrado")
                
                m.copy("FILTRO_BLAST_README",blast_dir+'/FILTRO_BLAST_README')
                m.borrado("FILTRO_BLAST_README")
                
                m.borrado(query)
                m.borrado(archivo_genbank)
                
                break
                
            except:
                # Control de posibles errores y escritura de ellos en el archivo log
                error_filtro=" Ha ocurrido un error. Valores introducidos por teclado no validos. Por favor, introduzca números"
                print()
                print()
                print("----------------------------------------------------------------------------------------")
                print(error_filtro)
                print("----------------------------------------------------------------------------------------")
                log.write("\n------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.write(error_filtro)
                log.write("\n------------------------------------------------------------------------------------\n")



#####################################    MUSCLE    #########################################   
                

# El alineamiento y la creacion del árbol filogenético se ejecutan automáticamente
# justo después del filtrado.

# Como no requieren ningún input, no se hace un bucle while --> no necesita un 
# menu para introducir los parametros del muscle

#Como depende directamente de los archivos generados en el filtrado, no hace falta 
#hacer un control de errores --> si algo falla, el filtrado no se realiza, y por tanto
# tampoco esta parte ni la siguiente.
                
        si_muscle=input("""
¿ Desea analizar con muscle (alineamiento de secuencias y árbol filogenético) los hits obtenidos en el blast ?
                          (Si/No):    """)
        
        if si_muscle in ["SI","SÍ","Si","Sí","si","sí","YES","Yes","yes","Y","S","s"]:
        # REALIZACIÓN DEL ALINEAMIENTO Y DEL ARBOL FILOGENETICO:
        
            #Bucle para realizar el muscle sobre cada uno de las secuencias query
            for name_query in lista_name_query:
                musc.muscle(name_query+"_multifasta_filtrado", name_query)
                
                           
            #impresion del progreso por pantalla
            print()
            print("-----------------------------------------------------------------------------------------------")
            buen_muscle="""
    ALINEAMIENTO Y ÁRBOL FILOGENÉTICO REALIZADO CON EXITO. Se han generado los 
    siguientes archivos, guardados en la carpeta results/MUSCLE:
       
     - muscle_alineamiento_name_query --> alineamiento de las secuencias proteicas en fmt fasta
     - muscle_neighbor_tree_name_query --> archivo con formato newick con el arbol filogenético
     - MUSCLE_README --> Archivo con los parametros y proceso de alineamiento y creacion 
                     del árbol filogenético
    
    * Se ha generado el mismo numero de archivos que de querys analizados


-----------------------------------------------------------------------------------------------
                """
            print(buen_muscle)
            
            #escritura del archivo log de control del proceso del script
            log.write("\n------------------------------------------------------------------------------------\n")
            log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
            log.write(buen_muscle)
            log.write("\n------------------------------------------------------------------------------------\n")
                
            # creacion de carpetas y traslado de archivos
            muscle_dir=results_dir+'/MUSCLE'
            m.mkdir(muscle_dir)
            
            #bucle para poder y borrar las carpetas generadas --> se han generado tres archivos por cada secuencia query
            for name_query in lista_name_query:
                m.copy(name_query+"_muscle_alineamiento",muscle_dir+'/'+name_query+'_muscle_alineamiento')
                m.borrado(name_query+"_muscle_alineamiento")
                
                m.copy(name_query+"_muscle_neighbor_tree",muscle_dir+'/'+name_query+'_muscle_neighbor_tree')
                m.borrado(name_query+"_muscle_neighbor_tree")
                
                m.copy(name_query+"_MUSCLE_README",muscle_dir+'/'+name_query+'_MUSCLE_README')
                m.borrado(name_query+"_MUSCLE_README")
        
            
####################################    PROSITE    #########################################


# La búsqueda de dominios se ejecuta automáticamente justo después del
# alineamiento y la creacion del árbol filogenético.

# No requiere ningun input de archivos, solo una respuesta de confirmacion para evitar fallos en prosite.
        
        print("\n-----------------------------------------------------------------------------------------\n")
        print()
        print("Por Favor, antes de continuar, asegurese que los archivos los archivos prosite.dat y prosite.doc, están en la misma carpeta desde la que esta corriendo este script")
        print("(Si no es así, antes de continuar añada dichos archivos a dicha carpeta)")
        print()
        si_prosite=input("""
¿ Desea analizar los dominios proteicos de cada hit obtenido en cada blast realizado ?

    **** (Si ha analizado un archivo multifasta con muchos query, esta búsqueda ****
    *****       tardará un rato y consumirá bastantes recursos)                 ****

    ### Si selecciona No, el script terminará. ####  

                   (Si/No):    """)
        
        if si_prosite in ["SI","SÍ","Si","Sí","si","sí","YES","Yes","yes","Y","S","s"]:
            
            #creacion de una lista vacia que se rellenara más abajo y se empleara en 
            #para dar ka informacion extra acerca de prodoc
            lista_dominios=[]
            
            # Control de errores de la busqueda de dominios en prosite
            try:
            
                # BÚSQUEDA DE DOMINIOS EN SECUENCIAS PROTEICAS:
                
                #Bucle para realizar la busqueda de dominios sobre cada una de las secuencias query
                for name_query in lista_name_query:
                    
                    # la funcion prosite crea un archivo .tsv con los resultados de la busqueda de dominios
                    # y tambien devuelve una lista con los nombres de todos los dominios
                    lista_actual=pro.prosite(name_query+"_multifasta_filtrado", "prosite.dat", name_query)
                    
                    # con este bucle lo que hacemos es completar la lista con los nombres de todos 
                    # los dominios encontrados en TODAS las secuencias proteicas analizadas
                    for dom in lista_actual:
                        if dom not in lista_dominios:
                            lista_dominios.append(dom)
                            
                # invocamos la funcion ordena para dar la lista de dominios ordenada
                # (objetivo --> facilitar la vida al usuario a la hora de seleccionar el dominio que quiera)
                lista_dominios_ordenada = o.Ordena(lista_dominios)
                
                    
                #impresion del progreso por pantalla
                print()
                print("------------------------------------------------------------------------------------------------------")
                buen_prosite="""
    BÚSQUEDA DE DOMINIOS EN SECUENCIAS PROTEICAS REALIZADA CON EXITO. Se han generado los 
    siguientes archivos, guardados en la carpeta results/PROSITE:
       
     - name_query_dominios_proteicos.tsv --> archivo con fmt .tsv que contiene los dominios proteicos
       que se han encontrado en cada secuencia proteica
     
    * Se ha generado el mismo numero de archivos que de querys analizados
------------------------------------------------------------------------------------------------------
                    """
                print(buen_prosite)
                
                #escritura del archivo log de control del proceso del script
                log.write("\n------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.write(buen_prosite)
                log.write("\n------------------------------------------------------------------------------------\n")
                    
                # creacion de carpetas y traslado de archivos 
                prosite_dir=results_dir+'/PROSITE'
                m.mkdir(prosite_dir)
                
                # bucle para borrar y mover los archivos generados por cada query
                for name_query in lista_name_query:
                   
                    m.borrado(name_query+"_multifasta_filtrado")
                
                    m.copy(name_query+"_dominios_proteicos.tsv",prosite_dir+'/'+name_query+'_dominios_proteicos.tsv')
                    m.borrado(name_query+"_dominios_proteicos.tsv")
                    
            except:
                
                # control de errores inesperados y desconocidos al parsear el archivo prosite.dat 
                error_prosite="""
Ha ocurrido un error al realizar la busqueda de dominios en prosite. Asegurese de que los archivos prosite.dat y prosite.doc,
están en la misma carpeta desde la que esta corriendo este script antes de realizar la siguiente ejecucion del programa.

Tras este error desconocido, el programa forzará la finalizacion del script como si el usuario hubiese presionado "Ctr+C". 
                """
                print()
                print("\n------------------------------------------------------------------------------------\n")
                print(error_prosite)
                print("\n------------------------------------------------------------------------------------\n")
                log.write("\n------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.write(error_prosite)
                log.write("\n------------------------------------------------------------------------------------\n")
                raise KeyboardInterrupt           
 

####################################    PRODOC    #########################################

# Tras haber realizado la funcion prosite, se da la opcion de consultar informacion extra 
# acerca de los dominios en prosite.doc --> por ello las lineas asociadas a la funcion prodoc
# tienen esta tabulacion --> estan dentro del if del si_prosite del modulo anterior
               
            si_prodoc = input("""
     ¿ Desea obtener información extra procedente del archivo prosite.doc, de alguno de los hits anteriores ? 
                                (Si/No):    """)
            
            try:
                #busqueda de informacion en el archivo prosite.txt
                if si_prodoc in ["SI","SÍ","Si","Sí","si","sí","YES","Yes","yes","Y","S","s"]:
                    doc.prodoc(lista_dominios_ordenada,"prosite.txt")
                else:
                    pass
            except:
                # control de errores inesperados y desconocidos al parsear el archivo prosite.txt
                error_prodoc="""
Ha ocurrido un error desconocido al buscar informacion en el archivo prodoc.txt
Ejecute el modulo prodoc individualmente para visualizar dicho error y tratar de solventarlo.
Muy posiblemente sea un error de formato en el archivo prosite.txt (pues proviene de un archivo .doc)
                """
                
                # impresion por pantalla y escritura del error en el archivo log 
                print()
                print("\n------------------------------------------------------------------------------------\n")
                print(error_prodoc)
                print("\n------------------------------------------------------------------------------------\n")
                log.write("\n------------------------------------------------------------------------------------\n")
                log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.write(error_prodoc)
                log.write("\n------------------------------------------------------------------------------------\n")
                
                # se lanza una excepcion para forzar la finalizacion del script
                raise KeyboardInterrupt 
                
################  CONTROL DE INTERRUPCIÓN DE LA EJECUCION POR TECLADO  ###################
        
    except KeyboardInterrupt:
        
        print()
        print()
        keyInterrupt="""
            ###################################################################
               Por orden del usario o por algun error desconocido durante la
               ejecución, se ha detenido el script y forzado su finalizacion 
            ###################################################################
        """
        print(keyInterrupt)
        print()
        print()
        
        # escritura del archivo log
        log.write("\n------------------------------------------------------------------------------------\n")
        log.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
        log.write(keyInterrupt)
        log.write("\n------------------------------------------------------------------------------------\n")
        log.close()
        
        #trasladamos el archivo log a la carpeta results y lo borramos de la carpeta principal donde estan todos los scripts
        m.copy("log", results_dir)
        m.borrado("log")
        
        #### si el usuario hace ctrl + c a mitad de un analisis realiza el siguiente 
        #### borrado de los archivos intermedios que se generan
       
        
        #creacion de una lista con los nombres de todos los archivos creados --> si el usuario hace control c, se borran 
        #todos estos archivos, para asi evitar que se queden procesos a medias
        lista_archivos=["subject_multifasta_gbk","hits_blast_no_filtrados","BD_subject.psq","BD_subject.pin","BD_subject.phr","FILTRO_BLAST_README"]
        for name_query in lista_name_query:
            lista_archivos.append(name_query+"_multifasta_filtrado")
            lista_archivos.append(name_query+"_blast_filtrado")
            lista_archivos.append(name_query+"_muscle_alineamiento")
            lista_archivos.append(name_query+"_muscle_neighbor_tree")
            lista_archivos.append(name_query+"_MUSCLE_README")
            lista_archivos.append(name_query+"_dominios_proteicos.tsv")
            
        for file in lista_archivos:
            if os.path.exists(file):
                m.borrado(file)
        
        sys.exit()


############## SENTENCIA FINAL DEL CODIGO DE LA FUNCION MAIN--> FINALIZACION Y CIERRE DEL ARCHIVO LOG ############
        
    # Si todo ha ido bien cerramos el archivo log y lo trasladamos a la carpeta de result
    log.close()
    m.copy("log", results_dir)
    m.borrado("log")
    
    return    





############################################################################################################
###################################### CÓDIGO PRINCIPAL DE MAIN.PY #########################################
############################################################################################################

if __name__ == "__main__":
    print("""
    
###############################################################################
    
 ESTE SCRIPT HA SIDO DESARROLLADO POR: 
       --------------------------------------------------
       |  ARTURO P. REDONDO LUQUE, BIOTECNOLOGÍA (UPM)  | 
       --------------------------------------------------
       
###############################################################################
       """) 
    main()
    
    print("""
-------------------------------------------------------------------------------          


   #################################################################
   ## EL SCRIPT SE HA EJECUTADO SATISFACTORIAMENTE HASTA EL FINAL ##
   #################################################################


RECORDATORIO:
    
- Al parsear el genbank se genera el archivo:
    
        subject_multifasta_gbk --> secuencias proteicas en fmt fasta del gbk
        
        
- Al realizar el blast se genera el archivo:
        
        blast_no_filtrado --> contiene el resultado del blast en fmt .tsv 
        

- Al filtrar el blast se generan 3 archivos:
    
        blast_filtrado --> contiene el resultado del blast filtrado en fmt .tsv
        Se genera un archivo por cada secuencia query analizada
    
        multifasta_filtrado --> secuencias proteicas filtradas del blast anterior 
        y en fmt fasta. Se genera un archivo por cada secuencia query analizada
        
        FILTRO_BLAST_README --> contiene los parametros de filtrado del blast
        
        
- Al emplear muscle se generan 3 archivos por cada secuencia query analizada:
        
        muscle_alineamiento --> secuencias fasta alineadas
        
        muscle_neighbor_tree --> arbol filogenetico de las proteinas anteriores
        
        MUSCLE_README --> archivo de texto con el proceso y los parametros
        empleados tanto para el alineamiento como para la creación del árbol 
        
        
- Al emplear prosite se genera un archivo por cada secuencia query analizada:
    
        dominios_proteicos.tsv --> archivo .tsv que muestra los dominios 
        que presenta cada secuencia proteica 
        
- ARCHIVO .log --> Archivo de control que guarda un registro del proceso y 
                   de los errores generados durante la ejecución del script
        
-------------------------------------------------------------------------------        
    
**** SCRIPT DESARROLLADO POR ARTURO P. REDONDO LUQUE, BIOTECNOLOGÍA (UPM) ****

-------------------------------------------------------------------------------
        """)



