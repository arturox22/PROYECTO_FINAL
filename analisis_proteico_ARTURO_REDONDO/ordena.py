#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# FUNCION DESARROLLADA POR ARTURO P. REDONDO LUQUE EN LA ASIGNATURA DE PROGRAMACION
# DE 3º CURSO DE BIOTECNOLOGIA EN EL CURSO 2018-19

def Ordena(Lista):
    """
    DEFNICION: Funcion que ordena la lista de strings que le pases como parámetros
    
    ARGUMENTOS:
    - Lista = Una lista ordenada o desordenada de palabras
    
    RESULTADO:
        - Devuelve la lista ordenada según el criterio español usual
    """
    
#################### FUNCION AUXILIAR PARA ORDENAR UNA LISTA ##################
    def Delante(Cadena1,Cadena2):
        """
        DEFINICION: Función que te indica si el primer string pasado se coloca 
                    antes o no que el segundo string que le pases
            
        ARGUMENTOS:
            - Cadena1, Cadena2=Dos cadenas de caracteres

        RESULTADO:
            - Usando la ordenación alfabética española usual, devuelve True si Cadena1 
              se coloca antes que Cadena2 y False encaso contrario.
            - Si Cadena1 y Cadena2 son iguales, devuelve True
        """

        # Ponemos en minúsculas las cadenas
        Cadena1=Cadena1.lower()
        Cadena2=Cadena2.lower()

        if Cadena1==Cadena2:
            return True

        Vocales={"á":"a","é":"e","í":"i","î":"i","ó":"o","ú":"u","ü":"u"}
       #Se quitan los acentos
        Cadena_L_1,Cadena_L_2=Cadena1,Cadena2
        for vocal in Vocales.keys():
            Cadena_L_1=Cadena_L_1.replace(vocal,Vocales[vocal])
            Cadena_L_2=Cadena_L_2.replace(vocal,Vocales[vocal])
        if Cadena_L_1==Cadena_L_2:
            #El orden sólo depende de los acentos
            return Cadena1<Cadena2
        else:
            Cadenas=[list(Cadena_L_1),list(Cadena_L_2)]
            for i in range(2):
                for j in range(len(Cadenas[i])):
                    if Cadenas[i][j]!="ñ":
                        Cadenas[i][j]=ord(Cadenas[i][j])
                    else:
                        Cadenas[i][j]=110.5
            return Cadenas[0]<Cadenas[1]

    ############# CÓDIGO PRINCIPAL DE LA FUNCIÓN ORDENA ##################
    
    # ORDENACION USANDO UN algoritmo recursivo de tipo quick-sort
    
    if len(Lista)<2: 
        pass 
        
    elif len(Lista)==2:
        if not Delante(Lista[0],Lista[1]): 
            Lista[0],Lista[1]= Lista[1],Lista[0]
            
    else:
        
        pivote=Lista[int(len(Lista)/2)]            
        L1,L2=[],[] 
        for i in range(len(Lista)): 
            if Delante(Lista[i],pivote):
                if i == int(len(Lista)/2):
                    pass
                else:
                    L1.append(Lista[i])        
            else:
                L2.append(Lista[i])
    
        #enviamos y componemos la lista
        Lista=Ordena(L1)+ [pivote]+ Ordena(L2)
    return Lista
