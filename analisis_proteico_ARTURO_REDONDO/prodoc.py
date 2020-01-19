#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.ExPASy import Prodoc

def prodoc(lista_hits,prosite_doc):
    
    """"
    DEFINICIÓN= Funcion que imprime por pantalla información extra contenida en prosite.doc del
    dominio proteico que le solicites
    
    ARGUMENTOS:
        - lista_hits = lista con los dominios proteicos obtenidos tras hacer un parseado en prosite.dat
                       (lista obtenida tras haber empleado el modulo "prosite.py")
        - prosite_doc= archivo prosite.txt  
    
    RESULTADO:
        - Esta funcion no genera ningun archivo ni devuelve nada, la info que proporciona la imprime por pantalla
    
    
        *prosite.txt y prosite.dat constituyen la base de datos de PROSITE*

    """
    print()
    print()

    while True:
        print(lista_hits)
        dominio=input("De la tabla anterior imprimida por pantalla, escriba el accession del dominio que quiera consultar:   ")
        if dominio in lista_hits:
#            try: 
            #APERTURA DEL ARCHIVO .doc
            doc = open(prosite_doc,"r")
            
            #PARSEADO DEL DOC, BUSCANDO EL DOMINIO DEL HIT SELECCIONADO
            records = Prodoc.parse(doc)
            for record in records:
                
                accession = record.accession
                
                # conversion del accession del .doc al accession del .dat
                if accession == dominio:
                    print("-----------------------------------------------------------------------------")
                    print()
                    print(""""
                ###########################################################         
                #### Información extra de la base de datos prosite.doc ####
                ########################################################### 
                """)
                    print()
                    print
                    print("DOMINIO:  "+record.accession)
                    print()
                    print(record.text)
                    print("################################################################")

            #CIERRE DEL ARCHIVO .doc
            doc.close()
            
            # Menu para decidir si consultar mas accession o no
            mas_doc=input("¿ Quiere consultar informacion sobre algún dominio más ?  (Si/No):   ")
            if mas_doc in ["NO","No","no","n","N"]:
                print("""
----------------------------------------------------------------------------------------------------
        
        ##############################################################
        #### Muchas gracias por haber confiado en nuestro software. ##
        ##############################################################
        
        
                  """)
                break
            
                
        else:
            print()
            print("--------------------------------------------------------------------------------------")
            print("Respuesta no válida, por favor, introduzca algun accession de los mostrados en la tabla anterior, es decir, alguno de estos:")
            print(lista_hits)
            print("--------------------------------------------------------------------------------------")
                
            
