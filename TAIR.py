# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:26:14 2020

@author: gebruiker
"""
# Niet functionele imports
from memory_profiler import profile
import gc
from numpy import linspace # Alleen voor opmaak grafiek

# Eigen scripts
import file_processing
import FileObjects
import GUI_classes

# Benodigde imports
import re
import matplotlib.pyplot as plt
import tkinter as tk



@profile
def main():
    data_dir = "Data/"
    fasta_name = "TAIR10_pep_20101214.fa"
    gff3_name = "TAIR10_GFF3_genes.gff"
    pattern = "[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}"
    
    show_plot = True
    
    fasta_data = read_fasta(data_dir+fasta_name)
    gff3_data = read_gff3(data_dir+gff3_name)
    
    fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par, chr_lengths = process_data(fasta_data, gff3_data, pattern)       
    
    compare(fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par)
    
    del fasta_data
    del gff3_data
    
    gc.collect()
    
#    print_info(fasta_objects)
    
    if show_plot:
        frequencies = count_fasta_chr(fasta_objects)
        make_plot(frequencies)
    

def read_fasta(location):
    with open(location) as fasta_file:
        return file_processing.fix_linebreaks(fasta_file)


def read_gff3(location):
    with open(location) as gff3_file:
        return file_processing.process_csv(gff3_file, "\t")
        

def process_data(fasta_data, gff3_data, pattern):
    fasta_objects = []
    gff3_objects = []
    gff3_attrs_ID = []
    gff3_attrs_Par = []
    chr_lengths = {}

    for header, sequence in fasta_data:
        if not re.search(pattern, sequence) == None:
            fasta_objects.append(FileObjects.Fasta.create_simple_fasta(header, sequence))
           
    for index, data in enumerate(gff3_data):
        gff3_object = FileObjects.GFF3.create_from_raw_data(data)
        gff3_objects.append(gff3_object)
        gff3_attrs_ID.append(gff3_object.get_attributes().get_ID())
        gff3_attrs_Par.append(gff3_object.get_attributes().get_parent())
        if gff3_object.get_data_type() == "chromosome":
            try:
                chr_lengths[gff3_object.get_seqID()] = gff3_object.get_end()
            except Exception as ex:
                pass
#        if index > 500:
#            break
            
    return fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par, chr_lengths
        
        
def compare(fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par):   
    
#    for x in gff3_attrs:
#        if x == "AT1G01010.1":
#            print(True)
            
#    print("con",gff3_attrs[0:20])
    
    for index, fasta_object in enumerate(fasta_objects):
        indexes(gff3_attrs_Par, fasta_object, gff3_objects)
        start_stop(gff3_attrs_ID, fasta_object, gff3_objects)
#        if index > 5:
#            break

        
def start_stop(lst, fasta_object, gff3_objects):
    element = fasta_object.get_seqID().split(".")[0]
    
    gff3_object = gff3_objects[lst.index(element)]
    
    fasta_object.set_start(gff3_object.get_start())
    fasta_object.set_stop(gff3_object.get_end())
        
        
def indexes(lst, fasta_object, gff3_objects):
    element = fasta_object.get_seqID()
#    element = "AT1G01010.1"
    offset = -1
    
#    print("--------------------------------------")
    
    while True:
        
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return None  
        
#        print("off",offset,"ele",element)
        update_fasta_object(fasta_object, gff3_objects[offset])
  

def update_fasta_object(fasta_object, gff3_object):
    fasta_object.add_gff3_hit(gff3_object) 
    fasta_object.set_chromosome(gff3_object.get_seqID())
    
    
def count_fasta_chr(fasta_objects):
    counts = dict()
    
    for fasta_object in fasta_objects:
        chromosome = fasta_object.get_chromosome()
        try:
            amount = counts.get(chromosome)
        except Exception as ex:
            print("d", ex)
            amount = 0
           
        if amount == None:
            amount = 0
            
        counts[chromosome] = amount+1
    
    return counts
    
    
def make_plot(frequencies):
    fig = plt.bar(frequencies.keys(), frequencies.values(), align="center")

    data_amount = len(frequencies)
    colors = iter(plt.cm.rainbow(linspace(0,1,data_amount)))
    for i in range(data_amount):
        c = next(colors)
        fig[i].set_color(c)
        
    # Legenda en namen
    plt.legend(fig,frequencies.keys())
    plt.title("Aantal genen met Serine/Threonine kinase active site per chromosoom")
    plt.ylabel("Aantal genen")
    plt.xlabel("Chromosoom")
    
    GUI_classes.Plot(fig)
    

#    plt.show()
    







def print_info(fasta_objects):
    for x in fasta_objects:
        if len(x.get_hits()) > 1:
            print("Seq ID:",x.get_seqID())
            print("Gene start-stop",x.get_start(),"-->", x.get_stop())
            print("Exon starts",x.get_starts())
            print("Chr: ",x.get_chromosome())
            print("Num of exomes:",x.get_number_of_exomes())
            pos_chr_s = (int(x.get_start())/int(chr_lengths.get(x.get_chromosome())))
            pos_chr_e = (int(x.get_stop())/int(chr_lengths.get(x.get_chromosome())))
            print("Pos on chr: {}% to {}%".format(pos_chr_s, pos_chr_e))
            print("------------------------")



if __name__ == "__main__":
    main()