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

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

def open_files(fasta_file, gff3_file):
    fasta_data = __read_fasta(fasta_file)
    gff3_data = __read_gff3(gff3_file)
    
    return fasta_data, gff3_data


def load_data(fasta_data, gff3_data, pattern):
    
    print("Processing data...")

    fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par, chr_lengths = __process_data(fasta_data, gff3_data, pattern)

    print("Comparing data...")

    del fasta_data
    del gff3_data

    gc.collect()

    __compare(fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par)

    del gff3_attrs_ID
    del gff3_attrs_Par

    gc.collect()
    
    return fasta_objects, gff3_objects, chr_lengths

"""
Main methode voor het uitvoeren van het script vanuit hier.
De GUI roept bovenstaande methodes aan.
"""
@profile
def __main():
    data_dir = "Data/"
    fasta_name = "TAIR10_pep_20101214.fa"
    gff3_name = "TAIR10_GFF3_genes.gff"
    pattern = "[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}"

    show_plot = True

    print("Getting data...")

    fasta_data = __read_fasta(data_dir+fasta_name)
    gff3_data = __read_gff3(data_dir+gff3_name)

    print("Processing data...")

    fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par, chr_lengths = __process_data(fasta_data, gff3_data, pattern)

    print("Comparing data...")

    del fasta_data
    del gff3_data

    gc.collect()

    __compare(fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par)

    del gff3_attrs_ID
    del gff3_attrs_Par

    gc.collect()

#    print_info(fasta_objects)

    if show_plot:
        frequencies = __count_fasta_chr(fasta_objects)
        __make_plot(frequencies)


def __read_fasta(location):
    with open(location) as fasta_file:
        return file_processing.fix_linebreaks(fasta_file)


def __read_gff3(location):
    with open(location) as gff3_file:
        return file_processing.process_csv(gff3_file, "\t")


def __process_data(fasta_data, gff3_data, pattern):
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


def __compare(fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par):

#    for x in gff3_attrs:
#        if x == "AT1G01010.1":
#            print(True)

#    print("con",gff3_attrs[0:20])

    for index, fasta_object in enumerate(fasta_objects):
        __indexes(gff3_attrs_Par, fasta_object, gff3_objects)
        __start_stop(gff3_attrs_ID, fasta_object, gff3_objects)
#        if index > 5:
#            break


def __start_stop(lst, fasta_object, gff3_objects):
    element = fasta_object.get_seqID().split(".")[0]

    gff3_object = gff3_objects[lst.index(element)]

    fasta_object.set_start(gff3_object.get_start())
    fasta_object.set_stop(gff3_object.get_end())


def __indexes(lst, fasta_object, gff3_objects):
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
        __update_fasta_object(fasta_object, gff3_objects[offset])


def __update_fasta_object(fasta_object, gff3_object):
    fasta_object.add_gff3_hit(gff3_object)
    fasta_object.set_chromosome(gff3_object.get_seqID())


def __count_fasta_chr(fasta_objects):
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


def __make_plot(frequencies):
#    fig = plt.bar(frequencies.keys(), frequencies.values(), align="center")
#
#    data_amount = len(frequencies)
#    colors = iter(plt.cm.rainbow(linspace(0,1,data_amount)))
#    for i in range(data_amount):
#        c = next(colors)
#        fig[i].set_color(c)
#
#    # Legenda en namen
#    plt.legend(fig,frequencies.keys())
#    plt.title("Aantal genen met Serine/Threonine kinase active site per chromosoom")
#    plt.ylabel("Aantal genen")
#    plt.xlabel("Chromosoom")

    GUI_classes.Plot(frequencies)


#    plt.show()








def __print_info(fasta_objects):
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
    __main()
