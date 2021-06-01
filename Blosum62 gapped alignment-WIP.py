# Programa una función que reciba como parámetros dos PROTEÍNAS y una MATRIZ DE PUNTUACIONES BLOSUM62.
# Las secuencias peptídicas serán de distinta longitud (gapped alignment).
# La información con las puntuaciones debe tener la siguiente estructura:
#    A  T  C  G -
# A  4 -1 -2 -2 -2
# T -1  4  0 -2 -2
# C -2  0  4  1 -2
# G -2 -2  1  4 -2
# - -2 -2 -2 -2 -1
 
# La función debe calcular la puntuación del alineamiento de las dos secuencias peptídicas empleando las puntuaciones definidas. 
# La función debe devolver la puntuación final obtenida.
###############################


                                        ################################################
                                        ################### PRUEBAS ####################
                                        ################################################
secuencia1 = "AC-G"
secuencia2 = "ACTG"




def alineamiento_blosum62(sec1, sec2):
    """
    Funcion que recibe como parametros dos proteinas de longitud disimilar y emplea la matriz de sustituciones Blosum62 para calcular 
    la puntuacion del alineamiento.
    """



    ######################################################
    # DICCIONARIO CON LA MATRIZ DE PUNTUACIONES BLOSUM62 #
    ######################################################
    
    BLOSUM62 = {
    "A": {"A":4,  "T":-1, "C":-2, "G":-2}, 
    "T": {"A":-1, "T":4,  "C":0,  "G":-2}, 
    "C": {"A":-2, "T":0,  "C":4,  "G":1},  
    "G": {"A":-2, "T":-2, "C":1,  "G":4},
    "-": {"A":-1. "T":-1, "C":-1, "G":-1, "-":-2} 
    }
    
    
    
    ###################################################
    ######### INICIALIZACIÓN DE VARIABLES #############
    ###################################################
    
    score = 0                                           # Iniciamos una puntuación inicial de 0 antes de comenzar el alineamiento (=bucle for)
    mejor_puntuacion = -3                               # Iniciamos una variable que contenga el mejor alineamiento. Empezamos en -3 porque no hay alineamiento inferior a -3
    peor_puntuacion = 5                                 # Ídem para el peor alineamiento residuo-residuo
    sec1 = sec1.upper()
    sec2 = sec2.upper()

    for i in range(len(sec1)):
        if sec1[i] == "-":  # Coste por abrir un gap (delección/inserción)
            score -= 4
        score += BLOSUM62[sec1[i]][sec2[i]]
        if sec1[i] == sec2[i]: 
            score += BLOSUM62[sec1[i]][sec2[i]]
        
    print(score)
    


alineamiento_blosum62(secuencia1,secuencia2)