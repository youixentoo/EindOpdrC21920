# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 16:39:14 2020

@author: Thijs Weenink

Bestand om de GUI uit te voeren.

In alle andere bestand staat commentaar, niet PEP-8, want wie heeft daar nou zin in.
"""
import GUI_classes

def opdracht():
    TAIR_GUI = GUI_classes.TAIRGUI()
    TAIR_GUI.run()
    
    
if __name__ == "__main__":
    opdracht()