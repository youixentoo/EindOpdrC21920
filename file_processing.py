# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 13:26:06 2019

@author: Thijs Weenink

Methodes om bestanden te verwerken.
"""
import re

# Zorgt ervoor dat het bestand, met de line breaks, goed wordt verwerkt.
# Returneerd een zip() van de headers en sequenties
def fix_linebreaks(opened_file, fasta_start = ">"):
    headers = []
    sequences = []
    sequence = []

    for line in opened_file:
        if re.search(fasta_start, line) == None:
            sequence.append(line.strip("\n").strip("\t"))
        else:
            sequences.append("".join(sequence))
            sequence = []
            headers.append(line.strip("\n"))

    sequences.append("".join(sequence))

    # De eerste index van sequences is leeg, in plaats van fixen,
    # is het negeren van de eerste index makkelijker. Data gaat niet verloren.
    return zip(headers, sequences[1::])


def process_csv(opened_file, delimiter=","):
    lines = []

    for line in opened_file:
        lines.append(line.split(delimiter))

    return lines
