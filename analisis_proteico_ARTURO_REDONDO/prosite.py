#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

from Bio.ExPASy import Prosite

from Bio import SeqIO

def prosite(hits_blast_multifasta, prosite_dat, name_query):
    
    """
    DEFINICION: Función que dada una seciencia proteica, parsea la base de datos
                prosite.data, generando un archivo .tsv con los dominios proteicos 
                encontrados y devolviendo una lista de los PDOC_accession asociados
                a cada dominio.
    ARGUMENTOS:          
        - hits_blast_multifasta= archivo que contiene las secuencias fasta de los hits de blast
        - prosite_data= archivo prosite.dat
    
        *prosite.txt y prosite.dat constituyen la base de datos de PROSITE*
     
    RESULTADO:
        - Devuelve una lista con todos los dominios proteicos encontrados 
        - Genera un archivo tsv (dominios_proteicos.tsv) con el siguiente formato: 
            
      PROTEÍNA    NOMBRE_DOMINIO      ACCESION_DOMINIO    DESCRIPCION_DOMINIO     PATRON_ProFMT   PATRON_ReFMT


    **** Los archivos se escriben entre comillas dentro de los parentesis de la funcion ****
    
    """

#################################################################################################
###################################  VARIABLES GLOBALES  ########################################
#################################################################################################

 # variables de control empleadas en el menu del final del script
    lista_hits=[]


#################################################################################################
###################################  SUBFUNCION AUXILIAR ########################################
#################################################################################################

    def busqueda_dom(sec_prot, id_prot, prosite_dat, name_query):
        
        """
        
        DEFINICION: Funcóon que dada una secuencia proteica, analiza si tiene algun dominio
                    que este registrado en la base de datos de prosite
            
        ARGUMENTOS:
            - sec_prot= string que contiene la secuencia proteica a analizar
            - id_prot= string con el id de la proteina
            - prosite_data= archivo prosite.dat
            - name_query= String que representa el nombre del query que estamos analizando
        
        RESULTADO:
            - Genera un archivo tsv (dominios_proteicos.tsv) con el siguiente formato: 
    
   PROTEÍNA  NOMBRE_DOMINIO  ACCESION_DOMINIO  PDOC_ACCESSION   DESCRIPCION_DOMINIO   PATRON_Prosite  PATRON_Re
    
    
        **** Los archivos se escriben entre comillas dentro de los parentesis de la funcion ****
        
        """
        #=============================================================================================
        #==========================   SUBFUNCION DE LA SUBFUNCION AUXILIAR ===========================
        #=============================================================================================
        
        
        def patron_RE(pPRO):
            """
            DEFINICION: Función que transforma los patrones en formato prosite, a formato del 
                        modulo RE de python
                        
            ARGUMENTOS:
                - pPRO= string que contiene un patrón regular con formato de prosite
            
            RESULTADO:
                - La funcion devuelve un string con el patron regular en formato "del módulo RE"
            
            """ 
            
            pRE=pPRO
            PROSITE=["-",".","x","{","}","(",")","<",">",">]"]
            RE=["","",".","[^","]","{","}","^","$","]?$"]
            
            for i in range(len(RE)):
                pRE=pRE.replace(PROSITE[i],RE[i])
                
            return pRE
        #============================================================================================
        #============================================================================================
        #============================================================================================
        
    
        ###################   CÓDIGO PRINCIPAL DEL SCRIPT busqueda_dom.py  #################
        
        
        # APERTURA DEL ARCHIVO .dat y CREACIÓN DEL ARCHIVO dominios_proteicos.tsv
        dat = open(prosite_dat,"r",encoding="utf8")
        output= open(name_query+"_dominios_proteicos.tsv","a",encoding="utf8")
        
      # impresión por pantalla de la cabecera de la tabla del ouptut
        print()
        print("_\tPROTEÍNA_ID\tNOMBRE_DOMINIO\tACCESSION_DOMINIO\tPDOC_ACCESSION\tDESCRIPCION_DOMINIO\tPATRON_Prosite\tPATRON_Re")
        print()
        
        # escribimos la cabecera de la tabla en el archivo output 
        output.write("PROTEÍNA_ID\tNOMBRE_DOMINIO\tACCESSION_DOMINIO\tPDOC_ACCESSION\tDESCRIPCION_DOMINIO\tPATRON_Prosite\tPATRON_Re\n")
      
        # Variable de control, empleada para la representacion de la tabla imprimida por pantalla
        contador_hits=0
    
        #BUCLE QUE RECORRE EL ARCHIVO .dat en busca de todos los patrones existentes
        records = Prosite.parse(dat)
        for dom in records: 
            
            # si el dominio tiene un patron en prosite.dat
            if len(dom.pattern) != 0:
                pRE= patron_RE(dom.pattern) # conversion del patron en fmt de prosite a fmt de RE
                if re.search(pRE,sec_prot): # si el patron esta en la proteina, hacemos:
                    
                    #guardamos las siguientes variables
                    name_dom=dom.name
                    accession_dom=dom.accession
                    accession_pdoc=dom.pdoc
                    descrip_dom=dom.description
                    ProPattern= dom.pattern
                    
                    #sumamos uno al contador
                    contador_hits+=1
                    
                    #creamos la lista con los posibles valores a consultar en el doc
                    if accession_dom not in lista_hits:
                        lista_hits.append(accession_pdoc)
                    
                    #imprimimos por pantalla la tabla
                    print(str(contador_hits)+"\t"+id_prot+"\t"+name_dom+"\t"+accession_dom+"\t"+accession_pdoc+"\t"+descrip_dom+"\t"+ProPattern+"\t"+pRE)
                    print()
                    output.write(id_prot+"\t"+name_dom+"\t"+accession_dom+"\t"+accession_pdoc+"\t"+descrip_dom+"\t"+ProPattern+"\t"+pRE+"\n")
                    
        # CIERRE DE LOS ARCHIVOS .dat y output
        dat.close()
        output.close()
        
#################################################################################################
#################################################################################################
#################################################################################################



###################   CÓDIGO PRINCIPAL DEL SCRIPT prosite.py  #################
    
    
    # Creamos dos listas, una con las secuencias proteicas y otra con los nombres de cada proteina
    lista_sec_prot=[]
    lista_id_prot=[]
    proteinas=open(hits_blast_multifasta, "r", encoding="utf8")
    for prot in SeqIO.parse(proteinas,"fasta"):
        lista_sec_prot.append(str(prot.seq))
        lista_id_prot.append(str(prot.id))
    
    #### REVISAAR
    print("-----------------------------------------------------------------------------------------")
    print()
    print("Visualización del archivo:")
    print("                  "+name_query+"_dominios_proteicos.tsv")
    print()
    print("Lista de los hits encontrados al realizar el blast sobre "+name_query+":")
    print(lista_id_prot)
    print("Sobre cada hit=secuencia proteica, se realiza la búsqueda de dominios")
    print()
    
    # recorremos las listas recien creadas e invocamos al script 
    for pos in range(len(lista_sec_prot)):
        busqueda_dom(lista_sec_prot[pos],lista_id_prot[pos],prosite_dat, name_query)
    
    return lista_hits
   
