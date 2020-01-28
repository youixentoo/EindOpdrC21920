# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 17:15:30 2020

@author: Thijs Weenink

Classes voor het opslaan van de data uit de bestanden.
"""

"""
Slaat alle data op van het fasta bestand.

Hieraan kan data worden toegevoegd op basis van de GFF3() en GFF3Attributes() classes.
"""
class Fasta():

    def __init__(self, header, sequence, chromosome, number_of_exomes):
        self.header = header
        self.sequence = sequence
        self.chromosome = chromosome
        self.number_of_exomes = number_of_exomes
        self.gff3_hit = None
        self.gff3_hits = []
        self.start = None
        self.stop = None
        self.starts = []
        self.stops = []

    """
    String die geprint wordt als een object wordt geprint.
    """
    def __str__(self):
        return "{}\n{}\n{} {} {} {}".format(self.header, self.sequence, self.chromosome, self.start, self.stop, self.gff3_hit)


    """
    Methode om een fasta object aan te maken op basis van alleen een header en sequentie. (Alternative constructor)
    """
    @classmethod
    def create_simple_fasta(cls, header, sequence):
        return cls(header, sequence, None, 0)

    """
    Setters
    """
    def set_header(self, header):
        self.header = header

    def set_sequence(self, sequence):
        self.sequence = sequence

    def set_chromosome(self, chromosome):
        self.chromosome = chromosome

    def set_start(self, start):
        self.start = start

    def set_stop(self, stop):
        self.stop = stop

    def set_number_of_exomes(self, number_of_exomes):
        self.number_of_exomes = number_of_exomes

    def increase_number_of_exomes(self, amount):
        self.number_of_exomes + amount

    def set_gff3_hit(self, hit):
        self.gff3_hit = hit

    def add_gff3_hit(self, hit):
        self.gff3_hits.append(hit)
        self.update_vars(hit)

    """
    Updates de starts en stops van de exonen.
    """
    def update_vars(self, gff3_object):
        if gff3_object.get_data_type() == "exon":
            self.starts.append(gff3_object.get_start())
            self.stops.append(gff3_object.get_end())


    """
    Getters
    """
    def get_header(self):
        return self.header

    def get_sequence(self):
        return self.sequence

    def get_chromosome(self):
        return self.chromosome

    def get_start(self):
        return self.start

    def get_stop(self):
        return self.stop

    def get_number_of_exomes(self):
        return len(self.starts)

    def get_seqID(self):
        return self.header.split("|")[0].strip(">").strip()

    def get_hit(self):
        return self.gff3_hit

    def get_hits(self):
        return self.gff3_hits

    def get_starts(self):
        return self.starts

    def get_stops(self):
        return self.stops
    
    """
    Om er voor te zorgen dan forward en reverse strands geen problemen geven,
    is een if-else nodig.
    """
    def get_length_gene(self):
        if int(self.start) > int(self.stop):
            return int(self.start) - int(self.stop)
        else:
            return int(self.stop) - int(self.start)
    
    """
    Geeft de string terug die in het textfield kan worden gezet in de GUI, 
    kan ook normaal worden geprint.
    """
    def get_textfield_data(self, chr_lengths):
        # Index op basis van verhouding ten opzichte van lengte chromosoom.
        pos_chr_s = round(int(self.start)/int(chr_lengths.get(self.chromosome)))
        pos_chr_e = round(int(self.stop)/int(chr_lengths.get(self.chromosome)))
        #print(pos_chr_s, "-->", pos_chr_e)
        
        chr_rep = "----------------------------------------------------------------------------------------------------"
#        chr_r = chr_rep[:pos_chr_s] + "=" + chr_rep[pos_chr_s+1:]
        temp = list(chr_rep)
        temp[pos_chr_s] = "="
        temp[pos_chr_e] = "="
        chr_loc = "".join(temp)
        
        
        return "Amount of exomes: {}\nTotal length gene: {}\nChromosome: {}\nLocation on chromosome:\n{}".format(self.get_number_of_exomes(), self.get_length_gene(), self.get_chromosome(), chr_loc)


"""
Slaat alle data op uit het gff bestand.

De variabele attributes, is het object GFF3Attributes()
"""
class GFF3():

    def __init__(self, seqID, source, data_type, start, end, score, strand, phase, attributes):
        self.seqID = seqID
        self.source = source
        self.data_type = data_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = self.set_attributes(attributes)

    """
    String die geprint wordt als een object wordt geprint.
    """
    def __str__(self):
#        return self.seqID
        return "{} {} {} {} {} {} {} {} {}".format(self.seqID, self.source, self.data_type, self.start, self.end, self.score, self.strand, self.phase, self.attributes)


    """
    Methode om een object aan te maken, direct van de lijst met data.
    """
    @classmethod
    def create_from_raw_data(cls, data_list):
        return cls(data_list[0], data_list[1], data_list[2], data_list[3], data_list[4], data_list[5], data_list[6], data_list[7], data_list[8])

    
    """
    Methode die de attributen van de gff3 data, omzet tot GFF3Attributes object.
    Dit object wordt dan in de GFF3 variable attributes gezet.
    """
    def set_attributes(self, attributes):
        data = attributes.split(";")
        gff3Atrrs = GFF3Attributes(None, None, None, None, None, None)

        if len(data) == 1:
            attrs = data[0].strip("\n").split("=")

            if attrs[0] == "ID":
                gff3Atrrs.set_ID(attrs[1])
            elif attrs[0] == "Name":
                gff3Atrrs.set_name(attrs[1])
            elif attrs[0] == "Parent":
                gff3Atrrs.set_parent(attrs[1])
            elif attrs[0] == "Derives_from":
                gff3Atrrs.set_derives_from(attrs[1])
            elif attrs[0] == "Note":
                gff3Atrrs.set_note(attrs[1])
            else:
                gff3Atrrs.set_undef(attrs[1])
        else:
            for item in data:
                try:
                    attrs = item.strip("\n").split("=")
                    if attrs[0] == "ID":
                        gff3Atrrs.set_ID(attrs[1])
                    elif attrs[0] == "Name":
                        gff3Atrrs.set_name(attrs[1])
                    elif attrs[0] == "Parent":
                        gff3Atrrs.set_parent(attrs[1])
                    elif attrs[0] == "Derives_from":
                        gff3Atrrs.set_derives_from(attrs[1])
                    elif attrs[0] == "Note":
                        gff3Atrrs.set_note(attrs[1])
                    else:
                        gff3Atrrs.set_undef(attrs[1])

                except Exception:
                    pass


        return gff3Atrrs

    """
    Getters
    """ 
    def get_seqID(self):
        return self.seqID

    def get_source(self):
        return self.source

    def get_data_type(self):
        return self.data_type

    def get_end(self):
        return self.end

    def get_start(self):
        return self.start

    def get_strand(self):
        return self.strand

    def get_phase(self):
        return self.phase

    def get_attributes(self):
        return self.attributes



"""
Slaat attribute data uit het gff3 bestand op.
"""
class GFF3Attributes():

    def __init__(self, ID, name, parent, derives_from, note, undef):
        self.ID = ID
        self.name = name
        self.parent = parent
        self.derives_from = derives_from
        self.note = note
        self.undef = undef

    
    """
    Methode om een object aan te maken op basis van de lijst met attribute data.
    """
    @classmethod
    def from_data(cls, attributes):
        data = attributes.split(";")
        gff3Atrrs = cls(None, None, None, None, None, None)

        for item in data:
            try:
                attrs = item.strip("\n").split("=")
                if attrs[0] == "ID":
                    gff3Atrrs.set_ID(attrs[1])
                elif attrs[0] == "Name":
                    gff3Atrrs.set_name(attrs[1])
                elif attrs[0] == "Parent":
                    gff3Atrrs.set_parent(attrs[1])
                elif attrs[0] == "Derives_from":
                    gff3Atrrs.set_derives_from(attrs[1])
                elif attrs[0] == "Note":
                    gff3Atrrs.set_note(attrs[1])
                else:
                    gff3Atrrs.set_undef(attrs[1])

            except Exception:
                pass


        return gff3Atrrs


    """
    String die geprint wordt als een object wordt geprint.
    """
    def __str__(self):
        return "ID:{} Name:{} Parent:{} Derives:{} Note:{} Undef:{}".format(self.ID, self.name, self.parent, self.derives_from, self.note, self.undef)

    
    """
    Setters
    """
    def set_ID(self, ID):
        self.ID = ID

    def set_name(self, name):
        self.name = name

    def set_parent(self, parent):
        self.parent = parent

    def set_derives_from(self, derives_from):
        self.derives_from = derives_from

    def set_note(self, note):
        self.note = note

    def set_undef(self, undef):
        self.undef = undef

    """
    Getters (Die relevant zijn)
    """
    def get_ID(self):
        return self.ID

    def get_parent(self):
        return self.parent   
