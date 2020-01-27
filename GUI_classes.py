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


        
