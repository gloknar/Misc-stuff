# Programa una función que reciba como parámetros dos PROTEÍNAS y una MATRIZ DE PUNTUACIONES BLOSUM62.
# Las secuencias peptídicas tendrán la misma longitud (ungapped alignment).
# La información con las puntuaciones debe tener la siguiente estructura:
#    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
# A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
# R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
# N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
# D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
# C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
# Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
# E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
# G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
# H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
# I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
# L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
# K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
# M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
# F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
# P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
# S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
# T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
# W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
# Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
# V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
# B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
# Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
# X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
# * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
 
# La función debe calcular la puntuación del alineamiento de las dos secuencias peptídicas empleando las puntuaciones definidas. 
# La función debe devolver la puntuación final obtenida.
###############################


                                        ################################################
                                        ################### PRUEBAS ####################
                                        ################################################
secuencia1 = "ACGTGCTAGCTAGCT*ATTTATCATTATATCG*"
secuencia2 = "ACCCCCCTATATCATATCACACATCGCAATCA*"




def alineamiento_blosum62(sec1, sec2):
    """
    Funcion que recibe como parametros dos proteinas de igual longitud y emplea la matriz de sustituciones Blosum62 para calcular 
    la puntuacion del alineamiento.
    
    
    >>> sec1 = "ATA"
    >>> sec2 = "ATA"
    >>> alineamiento_blosum62(sec1, sec2)
    (13, 5, 'T-T', 4, 'A-A')

    >>> sec1 = "HILMW"
    >>> sec2 = "ARN*F"
    >>> alineamiento_blosum62(sec1, sec2)
    (-11, 1, 'W-F', -4, 'M-*')
    """



    ######################################################
    # DICCIONARIO CON LA MATRIZ DE PUNTUACIONES BLOSUM62 #
    ######################################################
    
    BLOSUM62 = {
    "A": {"A":4,  "R":-1, "N":-2, "D":-2, "C":0,  "Q":-1, "E":-1, "G":0,  "H":-2, "I":-1, "L":-1, "K":-1, "M":-1, "F":-2, "P":-1, "S":1,  "T":0,  "W":-3, "Y":-2, "V":0,  "B":-2, "Z":-1, "X":0,  "*":-4}, 
    "R": {"A":-1, "R":5,  "N":0,  "D":-2, "C":-3, "Q":1,  "E":0,  "G":-2, "H":0,  "I":-3, "L":-2, "K":2,  "M":-1, "F":-3, "P":-2, "S":-1, "T":-1, "W":-3, "Y":-2, "V":-3, "B":-1, "Z":0,  "X":-1, "*":-4}, 
    "N": {"A":-2, "R":0,  "N":6,  "D":1,  "C":-3, "Q":0,  "E":0,  "G":0,  "H":1,  "I":-3, "L":-3, "K":0,  "M":-2, "F":-3, "P":-2, "S":1,  "T":0,  "W":-4, "Y":-2, "V":-3, "B":3,  "Z":0,  "X":-1, "*":-4},  
    "D": {"A":-2, "R":-2, "N":1,  "D":6,  "C":-3, "Q":0,  "E":2,  "G":-1, "H":-1, "I":-3, "L":-4, "K":-1, "M":-3, "F":-3, "P":-1, "S":0,  "T":-1, "W":-4, "Y":-3, "V":-3, "B":4,  "Z":1,  "X":-1, "*":-4}, 
    "C": {"A":0,  "R":-3, "N":-3, "D":-3, "C":9,  "Q":-3, "E":-4, "G":-3, "H":-3, "I":-1, "L":-1, "K":-3, "M":-1, "F":-2, "P":-3, "S":-1, "T":-1, "W":-2, "Y":-2, "V":-1, "B":-3, "Z":-3, "X":-2, "*":-4}, 
    "Q": {"A":-1, "R":1,  "N":0,  "D":0,  "C":-3, "Q":5,  "E":2,  "G":-2, "H":0,  "I":-3, "L":-2, "K":1,  "M":0,  "F":-3, "P":-1, "S":0,  "T":-1, "W":-2, "Y":-1, "V":-2, "B":0,  "Z":3,  "X":-1, "*":-4}, 
    "E": {"A":-1, "R":0,  "N":0,  "D":2,  "C":-4, "Q":2,  "E":5,  "G":-2, "H":0,  "I":-3, "L":-3, "K":1,  "M":-2, "F":-3, "P":-1, "S":0,  "T":-1, "W":-3, "Y":-2, "V":-2, "B":1,  "Z":4,  "X":-1, "*":-4}, 
    "G": {"A":0,  "R":-2, "N":0,  "D":-1, "C":-3, "Q":-2, "E":-2, "G":6,  "H":-2, "I":-4, "L":-4, "K":-2, "M":-3, "F":-3, "P":-2, "S":0,  "T":-2, "W":-2, "Y":-3, "V":-3, "B":-1, "Z":-2, "X":-1, "*":-4}, 
    "H": {"A":-2, "R":0,  "N":1,  "D":-1, "C":-3, "Q":0,  "E":0,  "G":-2, "H":8,  "I":-3, "L":-3, "K":-1, "M":-2, "F":-1, "P":-2, "S":-1, "T":-2, "W":-2, "Y":2,  "V":-3, "B":0,  "Z":0,  "X":-1, "*":-4}, 
    "I": {"A":-1, "R":-3, "N":-3, "D":-3, "C":-1, "Q":-3, "E":-3, "G":-4, "H":-3, "I":4,  "L":2,  "K":-3, "M":1,  "F":0,  "P":-3, "S":-2, "T":-1, "W":-3, "Y":-1, "V":3,  "B":-3, "Z":-3, "X":-1, "*":-4}, 
    "L": {"A":-1, "R":-2, "N":-3, "D":-4, "C":-1, "Q":-2, "E":-3, "G":-4, "H":-3, "I":2,  "L":4,  "K":-2, "M":2,  "F":0,  "P":-3, "S":-2, "T":-1, "W":-2, "Y":-1, "V":1,  "B":-4, "Z":-3, "X":-1, "*":-4}, 
    "K": {"A":-1, "R":2,  "N":0,  "D":-1, "C":-3, "Q":1,  "E":1,  "G":-2, "H":-1, "I":-3, "L":-2, "K":5,  "M":-1, "F":-3, "P":-1, "S":0,  "T":-1, "W":-3, "Y":-2, "V":-2, "B":0,  "Z":1,  "X":-1, "*":-4}, 
    "M": {"A":-1, "R":-1, "N":-2, "D":-3, "C":-1, "Q":0,  "E":-2, "G":-3, "H":-2, "I":1,  "L":2,  "K":-1, "M":5,  "F":0,  "P":-2, "S":-1, "T":-1, "W":-1, "Y":-1, "V":1,  "B":-3, "Z":-1, "X":-1, "*":-4}, 
    "F": {"A":-2, "R":-3, "N":-3, "D":-3, "C":-2, "Q":-3, "E":-3, "G":-3, "H":-1, "I":0,  "L":0,  "K":-3, "M":0,  "F":6,  "P":-4, "S":-2, "T":-2, "W":1,  "Y":3,  "V":-1, "B":-3, "Z":-3, "X":-1, "*":-4}, 
    "P": {"A":-1, "R":-2, "N":-2, "D":-1, "C":-3, "Q":-1, "E":-1, "G":-2, "H":-2, "I":-3, "L":-3, "K":-1, "M":-2, "F":-4, "P":7,  "S":-1, "T":-1, "W":-4, "Y":-3, "V":-2, "B":-2, "Z":-1, "X":-2, "*":-4}, 
    "S": {"A":1,  "R":-1, "N":1,  "D":0,  "C":-1, "Q":0,  "E":0,  "G":0,  "H":-1, "I":-2, "L":-2, "K":0,  "M":-1, "F":-2, "P":-1, "S":4,  "T":1,  "W":-3, "Y":-2, "V":-2, "B":0,  "Z":0,  "X":0,  "*":-4}, 
    "T": {"A":0,  "R":-1, "N":0,  "D":-1, "C":-1, "Q":-1, "E":-1, "G":-2, "H":-2, "I":-1, "L":-1, "K":-1, "M":-1, "F":-2, "P":-1, "S":1,  "T":5,  "W":-2, "Y":-2, "V":0,  "B":-1, "Z":-1, "X":0,  "*":-4}, 
    "W": {"A":-3, "R":-3, "N":-4, "D":-4, "C":-2, "Q":-2, "E":-3, "G":-2, "H":-2, "I":-3, "L":-2, "K":-3, "M":-1, "F":1,  "P":-4, "S":-3, "T":-2, "W":11, "Y":2,  "V":-3, "B":-4, "Z":-3, "X":-2, "*":-4}, 
    "Y": {"A":-2, "R":-2, "N":-2, "D":-3, "C":-2, "Q":-1, "E":-2, "G":-3, "H":2,  "I":-1, "L":-1, "K":-2, "M":-1, "F":3,  "P":-3, "S":-2, "T":-2, "W":2,  "Y":7,  "V":-1, "B":-3, "Z":-2, "X":-1, "*":-4}, 
    "V": {"A":0,  "R":-3, "N":-3, "D":-3, "C":-1, "Q":-2, "E":-2, "G":-3, "H":-3, "I":3,  "L":1,  "K":-2, "M":1,  "F":-1, "P":-2, "S":-2, "T":0,  "W":-3, "Y":-1, "V":4,  "B":-3, "Z":-2, "X":-1, "*":-4}, 
    "B": {"A":-2, "R":-1, "N":3,  "D":4,  "C":-3, "Q":0,  "E":1,  "G":-1, "H":0,  "I":-3, "L":-4, "K":0,  "M":-3, "F":-3, "P":-2, "S":0,  "T":-1, "W":-4, "Y":-3, "V":-3, "B":4,  "Z":1,  "X":-1, "*":-4}, 
    "Z": {"A":-1, "R":0,  "N":0,  "D":1,  "C":-3, "Q":3,  "E":4,  "G":-2, "H":0,  "I":-3, "L":-3, "K":1,  "M":-1, "F":-3, "P":-1, "S":0,  "T":-1, "W":-3, "Y":-2, "V":-2, "B":1,  "Z":4,  "X":-1, "*":-4}, 
    "X": {"A":0,  "R":-1, "N":-1, "D":-1, "C":-2, "Q":-1, "E":-1, "G":-1, "H":-1, "I":-1, "L":-1, "K":-1, "M":-1, "F":-1, "P":-2, "S":0,  "T":0,  "W":-2, "Y":-1, "V":-1, "B":-1, "Z":-1, "X":-1, "*":-4}, 
    "*": {"A":-4, "R":-4, "N":-4, "D":-4, "C":-4, "Q":-4, "E":-4, "G":-4, "H":-4, "I":-4, "L":-4, "K":-4, "M":-4, "F":-4, "P":-4, "S":-4, "T":-4, "W":-4, "Y":-4, "V":-4, "B":-4, "Z":-4, "X":-4, "*":1} }
    
    
    
    ###################################################
    ######### INICIALIZACIÓN DE VARIABLES #############
    ###################################################
    
    score = 0                                           # Iniciamos una puntuación inicial de 0 antes de comenzar el alineamiento (=bucle for)
    mejor_puntuacion = -5                               # Iniciamos una variable que contenga el mejor alineamiento. Empezamos en -5 porque no hay alineamiento inferior a -5
    peor_puntuacion = 12                                # Ídem para el peor alineamiento residuo-residuo
    sec1= sec1.upper()
    sec2= sec2.upper()



    ###################################################
    ##################### CÓDIGO ######################
    ###################################################

    for i in range(len(sec1)):
        score += BLOSUM62[sec1[i]][sec2[i]]
        if BLOSUM62[sec1[i]][sec2[i]] > mejor_puntuacion:               # Aquí buscamos el mejor alineamiento residuo-residuo y guardamos su puntuación y qué residuos son
            mejor_puntuacion = BLOSUM62[sec1[i]][sec2[i]]
            par_residuos_max =  str(sec1[i]) + "-" + str(sec2[i])
        if BLOSUM62[sec1[i]][sec2[i]] < peor_puntuacion:                # Aquí buscamos el peor alineamiento residuo-residuo y guardamos su puntuación y qué residuos son
            peor_puntuacion = BLOSUM62[sec1[i]][sec2[i]]
            par_residuos_min = str(sec1[i]) + "-" + str(sec2[i])
    return (score, mejor_puntuacion, par_residuos_max, peor_puntuacion, par_residuos_min)


puntuacion_alineamiento, puntuacion_max, residuos_max, puntuacion_min, residuos_min = alineamiento_blosum62(secuencia1,secuencia2)

print (f"La puntuacion del alineamiento entre secuencias fue {puntuacion_alineamiento}, con el par de residuos \"{residuos_max}\" aportando la mayor puntuacion, {puntuacion_max}")
print (f"El par de residuos que produjo peor puntuacion durante el alineamiento fue \"{residuos_min}\", con una puntuacion de {puntuacion_min}")


# Pruebas de desarrollador
if __name__ == "__main__":
    import doctest
    doctest.testmod()