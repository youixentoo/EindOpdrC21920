# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 15:34:58 2020

@author: gebruiker
"""
import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
matplotlib.use('TkAgg')



class Plot():
    
    def __init__(self, fig):
        
        plot_window = tk.Tk()
        plot_window.wm_title("Plot")
        plot_embed = FigureCanvasTkAgg(fig, master=plot_window)
        plot_embed.show()
        plot_embed.get_tk_widget().pack(expand=True)
        
        plot_window.mainloop()
        
        