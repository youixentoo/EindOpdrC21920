# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 15:34:58 2020

@author: Thijs Weenink

Classes voor GUI's
"""
# tkinter imports
import tkinter as tk

# matplotlib imports
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


# Andere imports
import gc
#from memory_profiler import profile
from numpy import linspace # Alleen voor opmaak grafiek

# Eigen scripts
import TAIR


"""
Visualisatie van de grafiek in een FigureCanvasTkAgg()

Gebruikt np.linespace() en plt.cm.rainbow() voor kleuring.
"""
class Plot():

    def __init__(self, data):

        self.root = tk.Toplevel()
        self.root.wm_title("Genen per chromosoom plot")

        fig = Figure(figsize=(5, 4), dpi=100)

        plot = fig.add_subplot(111)
        
        # Een legenda is vrij nutteloos
        plot.set_title("Aantal genen met Serine/Threonine kinase active site")
        plot.set_ylabel("Aantal genen")
        plot.set_xlabel("Chromosoom")
        
        bars = plot.bar(data.keys(), data.values(), align="center")
        
        # Kleuring van de plot
        data_amount = len(data)
        colors = iter(plt.cm.rainbow(linspace(0,1,data_amount)))
        for i in range(data_amount):
            c = next(colors)
            bars[i].set_color(c)
        
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
#        print("you pressed {}".format(event.key))
        key_press_handler(event, self.canvas, self.toolbar)


    def _quit(self):
        self.root.quit()     # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent
                            # Fatal Python Error: PyEval_RestoreThread: NULL tstate

    # If you put root.destroy() here, it will cause an error if the window is
    # closed with the window manager.

"""
De GUI voor het script.
Roep .run() aan voor uitvoeren.
"""
class TAIRGUI():
    
    _fasta_objects = ["None"]
    _menu_setup = True
    
    """
    GUI elementen
    """
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
        
        self.start_var = tk.StringVar()
        self.start_var.set("Choose")
        self.start_var.trace("w", self._get_fasta)
        self.data_dropdown = tk.OptionMenu(self.root, self.start_var, *self._fasta_objects)
        self.data_dropdown.grid(row=7, column=0, columnspan=2, padx=2, pady=3, sticky="w")
#        self.start_var.trace("w", self._get_fasta)
#        self.data_dropdown.focus()
        
        self.label_under_menu = tk.Label(self.root, text="")
        self.label_under_menu.grid(row=8, column=0, padx=2, pady=5, sticky="w")
        
        self.results_textfield = tk.Text(self.root, state="normal", height=10, width=29)
        self.results_textfield.configure(font=('Courier New', 8)) # size 8 ipv 10
        self.results_textfield.grid(row=8, column=1, columnspan=2, padx=2, pady=5, sticky="w")
        
        self.bottom_filler = tk.Label(self.root, text="")
        self.bottom_filler.grid(row=9, column=0, columnspan=3)

    """
    Roep deze methode aan om de GUI te runnen.
    """
    def run(self):
        self.root.mainloop()
    
    """
    Voor als de GUI moet worden afgesloten via een script.
    """
    def _exit(self):
        self.root.quit()
        self.root.destroy()
        
    """
    Methode om de configuratie opties te resetten naar de originele waardes.
    """
    def _reset(self):
        filename = tk.StringVar(self.root, value="Data\TAIR10_pep_20101214.fa")
        filename_gff3 = tk.StringVar(self.root, value="Data\TAIR10_GFF3_genes.gff")
        pattern = tk.StringVar(self.root, value="[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}")
        
        self.file_entry["textvariable"] = filename
        self.file_entry_gff3["textvariable"] = filename_gff3
        self.pattern_entry["textvariable"] = pattern
        
    """
    Methode om de bestandsnamen op te halen uit de Entries en de data uit de bestanden
    in te laden.
    Geeft een foutmelding als de bestanden niet de correcte extensie hebben of niet gevonden
    kunnen worden.
    """    
#    @profile
    def _open_files(self):
        self.file_loaded_label["text"] = ""
        self.file_loaded_label.config(fg="BLACK")
        
        self.set_progress("No data loaded", "BLACK")
        
        fasta_file = self.file_entry.get()
        gff3_file = self.file_entry_gff3.get()
        
        if fasta_file.endswith(".fa") and gff3_file.endswith(".gff"):
            try:
                self.fasta_data, self.gff3_data = TAIR.open_files(fasta_file, gff3_file)
                self.file_loaded_label["text"] = "Files loaded"
                self.file_loaded_label.config(fg="Green")
            except IOError:
                self.file_loaded_label["text"] = "File Error"
                self.file_loaded_label.config(fg="RED")
        else:
            self.file_loaded_label["text"] = "Incorrect files"
            self.file_loaded_label.config(fg="RED")
       
    """
    Methode om de data uit de bestanden in te laden.
    Zie TAIR.py voor verdere uitleg.
    Geeft een foutmelding als er geen patroon is gevonden.
    """    
#    @profile
    def _load_data(self):
        try:
            self.pattern = self.pattern_entry.get()
            self.set_progress("", "Black")
            if not self.pattern == "":
                self.set_progress("Loading data...", "Orange")
                self.fasta_objects, self.gff3_objects, self.chr_lengths = TAIR.load_data(self.fasta_data, self.gff3_data, self.pattern)
                self.set_progress("Data loaded", "Green")
                
    #            self._set_menu()
                self.update_menu()
                
                if self.plot_var.get():
                    try:
                        TAIR.show_plot(self.fasta_objects)
                    except Exception as ex: # TAIR.make_plot() raises een Exception
                        self.reset_pattern()
                        
                del self.gff3_objects
                del self.pattern
                
                gc.collect()
                
            else:
                self.set_progress("Pattern is empty", "Red")
                
        except IndexError:
            self.file_loaded_label["text"] = "Incorrect files"
            self.file_loaded_label.config(fg="RED")           
            self.set_progress("Problem loading data", "RED")
            
        except AttributeError:
            if self.file_loaded_label["text"] == "No files opened":
                self.set_progress("Problem loading data, open file first", "RED")
            else:
                self.set_progress("Problem loading data", "RED")
            
    
    """
    Methode om een label aan te passen.
    """
    def set_progress(self, message, color):
        self.progress_label["text"] = message
        self.progress_label.config(fg=color)
        
    
    """
    Deze methode set alle seqID's uit de fasta objecten in het menu.
    """            
    def update_menu(self):
        values = [fasta.get_seqID() for fasta in self.fasta_objects]
        
        menu = self.data_dropdown["menu"]
        menu.delete(0, "end")
        for string in values:
            menu.add_command(label=string, 
                             command=lambda value=string: self.start_var.set(value))
        
            
    """
    Methode om text in het textveld te zetten.
    """      
    def _set_results_text(self, text):
        self.results_textfield.configure(state="normal")
        self.results_textfield.delete(1.0, tk.END)
        self.results_textfield.insert(tk.END, text)
        self.results_textfield.configure(state="disabled")
       
        
    """
    Methode die de geselecteerde waarde uit het OptionMenu haalt en de bebehorende
    data uit de lijst met fasta objecten haalt. Vervolgens wordt dit in het textveld gezet.
    """
    def _get_fasta(self, *args):
        try:
            selected = self.start_var.get()
            
            if len(selected) > 100:
#                print(selected, len(selected))
                pass
            else:
                result = TAIR.get_fasta_data(selected, self.fasta_objects)
#                print(result)
                
                if not result == None:
                    text = result.get_textfield_data(self.chr_lengths)
#                    print("Text:\n", text)
                    self._set_results_text(text)
        except AttributeError:
            pass
        
        
    def reset_pattern(self):
        pattern = tk.StringVar(self.root, value="[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}")

        self.pattern_entry["textvariable"] = pattern
   


"""
Aanroepen van de GUI
"""
if __name__ == "__main__":
    TAIRGUI().run()
    
    
    
    
    
    
    
    
    
    
    
    
    
        
