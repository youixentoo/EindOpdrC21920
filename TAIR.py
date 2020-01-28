# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:26:14 2020

@author: Thijs Weenink

Methodes voor data verwerking.
"""
# Niet functionele imports
from memory_profiler import profile
import gc

# Eigen scripts
import file_processing
import FileObjects
import GUI_classes

# Benodigde imports
import re

"""
Geeft de data uit de bestanden.
"""
def open_files(fasta_file, gff3_file):
    fasta_data = __read_fasta(fasta_file)
    gff3_data = __read_gff3(gff3_file)
    
    return fasta_data, gff3_data

"""
Geeft alle data op basis van de data die de bestanden geven en het meegegeven patroon.

Deletes sommige variabelen om geheugen vrij te maken.
"""
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
Geeft de frequenties van de chromosomen en laat de grafiek zien.
"""
def show_plot(fasta_objects):
    frequencies = __count_fasta_chr(fasta_objects)
    __make_plot(frequencies)
    
    
"""
Methode die een fasta object terug geeft op basis van het seqID
"""
def get_fasta_data(seqID, fasta_objects):
    for fasta_object in fasta_objects:
        if fasta_object.get_seqID() == seqID:
            return fasta_object
    else:
        return None
    
    

"""
Main methode voor het uitvoeren van het script vanuit hier.
De GUI roept bovenstaande methodes aan.

@profile is voor de import: memory_profiler, is verder nergens belangrijk voor.

Eerste 5 variabelen zijn voor configuratie van het script.
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

#    __print_info(fasta_objects, 1, chr_lengths)

    if show_plot:
        frequencies = __count_fasta_chr(fasta_objects)
        __make_plot(frequencies)


"""
Leest het fasta bestand in en geeft een zip() met header, sequentie terug.
"""
def __read_fasta(location):
    with open(location) as fasta_file:
        return file_processing.fix_linebreaks(fasta_file)


"""
Leest het gff3 bestand in en geeft alle regels in het bestand terug
"""
def __read_gff3(location):
    with open(location) as gff3_file:
        return file_processing.process_csv(gff3_file, "\t")

"""
Processes de data.

Eerst wordt er een lijst aangemaakt met fasta objecten waarin het patroon voorkomen.
Daarna worden er meerdere lijsten aangemaakt met gff3 data en een dictionary met chromosoom lengtes.
De attributen van gff3, worden in een apart object gezet.
"""
def __process_data(fasta_data, gff3_data, pattern):
    fasta_objects = [] # List
    gff3_objects = [] # List
    gff3_attrs_ID = [] # List
    gff3_attrs_Par = [] # List
    chr_lengths = {} # Dict

    for header, sequence in fasta_data:
        if not re.search(pattern, sequence) == None:
            fasta_objects.append(FileObjects.Fasta.create_simple_fasta(header, sequence))

    for data in gff3_data:
        gff3_object = FileObjects.GFF3.create_from_raw_data(data)
        gff3_objects.append(gff3_object)
        gff3_attrs_ID.append(gff3_object.get_attributes().get_ID())
        gff3_attrs_Par.append(gff3_object.get_attributes().get_parent())
        if gff3_object.get_data_type() == "chromosome":
            try:
                chr_lengths[gff3_object.get_seqID()] = gff3_object.get_end()
            except Exception as ex:
                print("Error in __process_data():\n",repr(ex))

    return fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par, chr_lengths


"""
Vergelijkt de fasta objecten met de gff3 objecten de overeenkomstige seqID's te vinden.
Bij een match wordt het gff3 object aan het fasta object toegevoegd.
"""
def __compare(fasta_objects, gff3_objects, gff3_attrs_ID, gff3_attrs_Par):
    for index, fasta_object in enumerate(fasta_objects):
        __indexes(gff3_attrs_Par, fasta_object, gff3_objects)
        __start_stop(gff3_attrs_ID, fasta_object, gff3_objects)


"""
De start en stop van het gen worden opgezocht op basis van de index van gff3_attrs_ID,
en overeenkomstige index in gff3_objects. Deze hebben dezelfde indexes/volgorde.
"""
def __start_stop(lst, fasta_object, gff3_objects):
    element = fasta_object.get_seqID().split(".")[0]

    gff3_object = gff3_objects[lst.index(element)]

    fasta_object.set_start(gff3_object.get_start())
    fasta_object.set_stop(gff3_object.get_end())


"""
Het seqID van het fasta object wordt opgezocht in de lijst met gff3 seqID's.
Elke index die wordt gevonden in deze lijst, komt overeen met het bijbehorende gff3 object,
het gevonden gff3 object wordt dan toegevoegd aan het fasta object.
gff3_objects en gff3_attrs_ID hebben dezelfde indexes/volgorde.
"""
def __indexes(lst, fasta_object, gff3_objects):
    element = fasta_object.get_seqID()
    offset = -1

    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return None

        __update_fasta_object(fasta_object, gff3_objects[offset])


"""
Voegt het gff3_object uit __indexes() toe aan het fasta object.
Voegt ook toe op welk chromosoom dit gen ligt.
"""
def __update_fasta_object(fasta_object, gff3_object):
    fasta_object.add_gff3_hit(gff3_object)
    fasta_object.set_chromosome(gff3_object.get_seqID())


"""
Telt hoe vaak elk chromosoom voorkomt in de lijst met fasta objecten.
Geeft dit terug als dictionary. (Lijkt op Collections Counter)
"""
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


"""
Maakt de grafiek
"""
def __make_plot(frequencies):
    GUI_classes.Plot(frequencies)
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

#    plt.show()


"""
Print info uit de fasta object lijst, op basis van de gegeven index.
"""

def __print_info(fasta_objects, fasta_index, chr_lengths):
    x = fasta_objects[fasta_index]
    print("Seq ID:",x.get_seqID())
    print("Gene start-stop",x.get_start(),"-->", x.get_stop())
    print("Exon starts",x.get_starts())
    print("Chr: ",x.get_chromosome())
    print("Num of exomes:",x.get_number_of_exomes())
    pos_chr_s = (int(x.get_start())/int(chr_lengths.get(x.get_chromosome())))
    pos_chr_e = (int(x.get_stop())/int(chr_lengths.get(x.get_chromosome())))
    print("Pos on chr: {}% to {}%".format(pos_chr_s, pos_chr_e))
    print("------------------------")



"""
Aanroepen main
"""
if __name__ == "__main__":
    __main()
