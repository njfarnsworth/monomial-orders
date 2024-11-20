import numpy as np
import itertools 
from itertools import permutations
from numpy import random
import csv

def I(n): #creates the standard lex matrix in n variables
    A = []
    for i in range(n):
        A.append([1 if j == i else 0 for j in range(n)])
    return A

def countP(A): #counts all of the possible permutations of an nxn matrix
    ct = 0
    for m in itertools.permutations(A):
        ct += 1
    return ct

def count(B): #counts all of the regular matrices in  a list of matrices 
    ct = 0
    for i in range(len(B)):
        ct += 1
    return ct

def permutations(A): #creates all of the permutations of an axa lex matrix & adds them to a list of arrays
    arrList = []
    for m in itertools.permutations(A):
        arrList.append(np.array(m))

    return arrList
    
def printer(arrayList): #prints the arrList containing all of the permutations
    for i in range(len(arrayList)):
        print(arrayList[i])
        print("")

def ultimatePrinter(listOfList):
    ct = 1
    for i in range(len(listOfList)):
        currentList = listOfList[i]

        for j in range(len(listOfList[i])):
            print("#" + str(ct))
            print(listOfList[i][j])
            ct+=1
            print("")

    
def mainalg(n):
    ultimateList = []
    permutation = permutations(I(n)) #creates all of the possible permutations for that matrix and adds them to a list


    for i in range(len(permutation)): #for each permutation
        permUltList = [] 
        currentPerm = permutation[i]
        permUltList.append(currentPerm) 

        for row in range(n-1, 0, -1): #decrementing loop... starting at last row, end at second to first row
            fixedCol = 0

            for j in range(len(permUltList)):
                currMatrix = permUltList[j]

                for col in range(n):
                    if currMatrix[row][col] == 1:
                        fixedCol = col
                    
                for r in range(row-1, -1, -1): #decrement through the rows starting with the row above the originally selected row
                    for c in range(0, fixedCol):
                        if currMatrix[r][c] == 1:
                            matCop = currMatrix.copy()
                            matCop[r][fixedCol] = 2

                            for k in range(r+1, row):
                                col1 = 0

                                for kc in range(n):
                                    if matCop[k][kc] == 1:
                                        col1 = kc

                                if matCop[k][fixedCol] == 0 and fixedCol > col1:
                                    matCop[k][fixedCol] = 3

                            negation = matCop.copy()
                            negation[row][fixedCol] = -1
                            
                            permUltList.append(matCop)
                            permUltList.append(negation)

       
        ultimateList.append(permUltList)

 
    print("COUNT: " + str(count(ultimateList)))


    return ultimateList


def main():

    list = mainalg(4)
    ct = 1
    f = open("info.csv", "w")
    for matrixList in list:
        f.write("\n")
        for matrix in matrixList:

            f.write(str(ct))
            f.write("\n")
            f.write(str(matrix))
            f.write("\n")
            f.write("\n")
          
            ct+=1

    
main()








