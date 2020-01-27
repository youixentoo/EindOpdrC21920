# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 15:34:58 2020

@author: gebruiker
"""
import tkinter as tk

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from numpy import linspace # Alleen voor opmaak grafiek

import numpy as np
from memory_profiler import profile

import gc

import TAIR

class Plot():

    def __init__(self, data):

        self.root = tk.Tk()
        self.root.wm_title("Genen per chromosoom plot")

        fig = Figure(figsize=(5, 4), dpi=100)
#        t = np.arange(0, 3, .01)

#        plot = plt.bar(data.keys(), data.values(), align="center")

        fig.add_subplot(111).bar(data.keys(), data.values(), align="center")

#        data_amount = len(data)
#        colors = iter(plt.cm.rainbow(linspace(0,1,data_amount)))
#        for i in range(data_amount):
#            c = next(colors)
#            plt[i].set_color(c)

        # Legenda en namen
#        plt.legend(fig,data.keys())
#        plt.title("Aantal genen met Serine/Threonine kinase active site per chromosoom")
#        plt.ylabel("Aantal genen")
#        plt.xlabel("Chromosoom")

        self.canvas = FigureCanvasTkAgg(fig, master=self.root)  # A tk.DrawingArea.
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.root)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.canvas.mpl_connect("key_press_event", self.on_key_press)

        self.button = tk.Button(master=self.root, text="Quit", command=self._quit)
        self.button.pack(side=tk.BOTTOM)

        self.root.mainloop()


    def on_key_press(self, event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, self.canvas, self.toolbar)


    def _quit(self):
        self.root.quit()     # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent
                            # Fatal Python Error: PyEval_RestoreThread: NULL tstate

    # If you put root.destroy() here, it will cause an error if the window is
    # closed with the window manager.


class TAIRGUI():
    
    _fasta_objects = ["placeholder1", "placeholder2", "placeholder3"]
    _menu_setup = True
    
    def __init__(self):
        self.root = tk.Tk()
        self.root.wm_title("TAIR GUI")
        
        self.root.resizable(False, False) # Grootte GUI is niet aan te passen
        
        self.filler1 = tk.Label(self.root, text="")
        self.filler1.grid(row=0, column=0, columnspan=3, padx=10)

        self.config_label = tk.Label(self.root, text="Configuration options:")
        self.config_label.grid(row=1, column=0, pady=5, columnspan=2, sticky="w")
        
        self.reset_cofig = tk.Button(self.root, text="Reset", command=self._reset, width=8)
        self.reset_cofig.grid(row=1, column=2, padx=2, sticky="e")
    
        self.file_label = tk.Label(self.root, text="Fasta filelocation:")
        self.file_label.grid(row=2, column=0, sticky="w")

        filename = tk.StringVar(self.root, value="Data\TAIR10_pep_20101214.fa")
        self.file_entry = tk.Entry(self.root, width=50, borderwidth=2, textvariable=filename)
        self.file_entry.grid(row=2, column=1, columnspan=2, padx=2, sticky="w")
        
        self.file_label_gff3 = tk.Label(self.root, text="GFF filelocation:")
        self.file_label_gff3.grid(row=3, column=0, sticky="w")
        
        filename_gff3 = tk.StringVar(self.root, value="Data\TAIR10_GFF3_genes.gff")
        self.file_entry_gff3 = tk.Entry(self.root, width=50, borderwidth=2, textvariable=filename_gff3)
        self.file_entry_gff3.grid(row=3, column=1, columnspan=2, padx=2, sticky="w")

        self.pattern_label = tk.Label(self.root, text="Pattern to search:")
        self.pattern_label.grid(row=4, column=0, sticky="w")


        pattern = tk.StringVar(self.root, value="[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}")
        self.pattern_entry = tk.Entry(self.root, width=50, borderwidth=2, textvariable=pattern)
        self.pattern_entry.grid(row=4, column=1, columnspan=2, padx=2, sticky="w")

        self.file_loaded_label = tk.Label(self.root, text="No files opened")
        self.file_loaded_label.grid(row=5, column=0, pady=5, sticky="w")

        self.progress_label = tk.Label(self.root, text="No data loaded")
        self.progress_label.grid(row=5, column=1, pady=5, sticky="w")
        
        self.files_button = tk.Button(self.root, text="Open...", command=self._open_files, width=8)
        self.files_button.grid(row=5, column=2, padx=2, sticky="e")

        self.fileMessage_label = tk.Label(self.root, text="")
        self.fileMessage_label.grid(row=6, column=0, columnspan=2, padx=10)
        
        self.plot_var = tk.BooleanVar()
        self.show_plot_checkbox = tk.Checkbutton(self.root, text="Show plot", variable=self.plot_var)
        self.show_plot_checkbox.grid(row=6, column=2, sticky="e")
        
        self.data_button = tk.Button(self.root, text="Load data", command=self._load_data, width=8)
        self.data_button.grid(row=7, column=2, padx=2, pady=3, sticky="e")
        
        self.start_var = tk.StringVar(self.root, "None")
        self.data_dropdown = tk.OptionMenu(self.root, self.start_var, *self._fasta_objects, command=self._get_fasta)
        self.data_dropdown.grid(row=7, column=0, columnspan=2, padx=2, pady=3, sticky="w")
#        self.start_var.trace("w", self._get_fasta)
        
        self.label_under_menu = tk.Label(self.root, text="")
        self.label_under_menu.grid(row=8, column=0, padx=2, pady=5, sticky="w")
        
        self.results_textfield = tk.Text(self.root, state="disabled", height=10, width=26)
        self.results_textfield.grid(row=8, column=1, columnspan=2, padx=2, pady=5, sticky="w")
        
        self.bottom_filler = tk.Label(self.root, text="")
        self.bottom_filler.grid(row=9, column=0, columnspan=3)

#        
   
    def run(self):
        self.root.mainloop()
    
    
    def _exit(self):
        self.root.quit()
        self.root.destroy()
        
        
    def _reset(self):
        filename = tk.StringVar(self.root, value="Data\TAIR10_pep_20101214.fa")
        filename_gff3 = tk.StringVar(self.root, value="Data\TAIR10_GFF3_genes.gff")
        pattern = tk.StringVar(self.root, value="[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}")
        
        self.file_entry["textvariable"] = filename
        self.file_entry_gff3["textvariable"] = filename_gff3
        self.pattern_entry["textvariable"] = pattern
        
    @profile
    def _open_files(self):
        self.file_loaded_label["text"] = ""
        self.file_loaded_label.config(fg="BLACK")
        
        fasta_file = self.file_entry.get()
        gff3_file = self.file_entry_gff3.get()
        
        if fasta_file.endswith(".fa") and gff3_file.endswith(".gff"):
            self.fasta_data, self.gff3_data = TAIR.open_files(fasta_file, gff3_file)
            self.file_loaded_label["text"] = "Files loaded"
            self.file_loaded_label.config(fg="Green")
        else:
            self.file_loaded_label["text"] = "Incorrect files"
            self.file_loaded_label.config(fg="RED")
            
    @profile
    def _load_data(self):
        self.pattern = self.pattern_entry.get()
        self.set_progress("", "Black")
        if not self.pattern == "":
            self.set_progress("Loading data...", "Orange")
            self.fasta_objects, self.gff3_objects, self.chr_lengths = TAIR.load_data(self.fasta_data, self.gff3_data, self.pattern)
            self.set_progress("Data loaded", "Green")
            
            del self.fasta_data
            del self.gff3_data
            
            gc.collect()
            
            self._set_menu()
            
            if self.plot_var.get():
                TAIR.show_plot(self.fasta_objects)
                
            del self.gff3_objects
            del self.pattern
            
            gc.collect()
            
        else:
            self.set_progress("Pattern is empty", "Red")
            
    
    def set_progress(self, message, color):
        self.progress_label["text"] = message
        self.progress_label.config(fg=color)
        
        
    def _set_menu(self):
        self._menu_setup = True
        self.start_var.set(self.fasta_objects[0].get_seqID())
        self.data_dropdown["menu"].delete(0, "end")
        for fasta_object in self.fasta_objects:
            self.data_dropdown["menu"].add_command(label=fasta_object.get_seqID(), command=self._get_fasta)
        self._menu_setup = False
        
"""
Exception in Tkinter callback
Traceback (most recent call last):
  File "C:\Users\gebruiker\Anaconda3\lib\tkinter\__init__.py", line 1705, in __call__
    return self.func(*args)
TypeError: _get_fasta() missing 1 required positional argument: 'value'
"""

    def _get_fasta(self, value, *args):

        result = TAIR.get_fasta_data(value, self.fasta_objects)
        print(result)
        
        if not result == None:
            text = result.get_textfield_data(self.chr_lengths)
            self._set_results_text(text)
            
            
    def _set_results_text(self, text):
        self.results_textfield.configure(state="normal")
        self.results_textfield.delete(1.0, tk.END)
        self.results_textfield.insert(tk.END, text)
        self.results_textfield.configure(state="disabled")



if __name__ == "__main__":
    TAIRGUI().run()
    
    
    
    
    
    
    
    
    
    
    
    
    
        
