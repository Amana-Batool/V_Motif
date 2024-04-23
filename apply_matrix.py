#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 20:44:07 2021

@author: amina
"""


#from get_file import file_open
import pandas as pd

global motif
motif = []

global count_chrctr 
global newarray
global input_file

def open_file1(file):
    
    count_chrctr = 0
    newarray = []
    
    input_file = file
    
    print (input_file)
    
    with open(input_file, "r") as f:
        for line in f:
            motif.append(line)    
            
    
    x = len (motif)
    print ("Total number of Sequences is: " , x)


    for i in range(x):
    
        new_motif = motif[i]
    
        xx = len(new_motif)
    
        count_chrctr = count_chrctr + xx
     
        for j in range (xx):
        
       
            newarray.append(new_motif[j:j+2])
            j=j+1
      
    
        i=i+1 

    total_chrctr = count_chrctr-x

    print ("Total chracters are : " ,total_chrctr)

    
    f = open("/home/student/Documents/file_data.txt", "w")
    f.close()
    
    ff = open("/home/student/Documents/file_data.txt", "a")
    ff.write(str(newarray))
    ff.close()
    
    
    return (newarray)