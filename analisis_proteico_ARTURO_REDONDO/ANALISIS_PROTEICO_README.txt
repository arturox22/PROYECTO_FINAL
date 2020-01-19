
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
        se recogen en un archivo llamado "log"

		       -------------------------------------------------------------------------
	RECOMENDACIÓN: |  Ejecute el script en una ventana completa para que el menú imprimido | 
		       |            por pantalla tenga mejor legibilidad		       |
		       -------------------------------------------------------------------------
    
   
	VERSION:        analisis_proteicos v1.0

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
		  




####################################################################################################

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
    **** que el script le proporcione la información que busca   
             

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
