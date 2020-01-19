#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from subprocess import Popen, PIPE

import os



def mkdir(nueva_carpeta):
    
    """
    DEFINICION: Función genera una nueva carpeta
    
    ARGUMENTOS:
        - nueva_carpeta= string. Nombre que tendrá la nueva carpeta generada
    
    RESULTADO:
        - Creación de una nueva carpeta en el directorio desde el que se llame
          al script
    
    """
    
    actual_directory=os.getcwd()
    
    if os.path.exists(nueva_carpeta):
        print("Ya existe una carpeta con ese nombre en el path actual (", actual_directory, ")")
    else:
        os.mkdir(nueva_carpeta)
    
def copy(archivo_origen, ruta_destino):
    
    """
    DEFINICION: Función que copia un archivo en la ruta que le indiques 
        
    ARGUMENTOS:
        - archivo= archivo del que queremos hacer una copia
        - ruta_destino= string. Lugar donde queremos copiar el archivo
    
                ej: ruta_destino = "data/archivo"
    
    RESULTADOS:
        - Se realiza una copia del archivo especificado en el destino especificado
        
        
    *** Asegurarse de que la carpeta data existe, para que no de fallos en la ruta
    
    OJO: los archivos se escriben entre "" dentro de los parentesis de la funcion
    
    """
    
    
    #COMANDO BASH-LINUX: cp archivo
    copi = Popen(['cp',archivo_origen, ruta_destino], stdout=PIPE, stderr=PIPE)
    copi.stderr.close()
    copi.stdout.close()
    
def borrado(archivo):
    
    """
    DEFINICION: Función que borra el archivo que le pases como parámetro
        
    ARGUMENTOS: 
        - archivo= archivo que queremos eliminar
    
    RESULTADO:
        - Borrado del archivo seleccionado
    
    """
    
    #COMANDO BASH-LINUX: cp archivo
    borra = Popen(['rm',archivo], stdout=PIPE, stderr=PIPE)
    borra.stderr.close()
    borra.stdout.close()


