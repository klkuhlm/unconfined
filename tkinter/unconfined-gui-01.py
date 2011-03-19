#!/usr/bin/env python
import Tkinter as tk

# import all default values in d dictionary
from unconfined-data import *

class GUI:
    def __init__(self,master):

        # left half of window
        lframe = tk.Frame(master)
        lframe.pack(side='left')
        
        # logical switches
        self.lquiet = tk.BooleanVar()
        self.lquiet = d['quiet']
        cbq = tk.Checkbutton(master, text="quiet ouput?", variable=self.lquiet,
                             onvalue=True, offvalue=False)
        cbq.pack(side='top')


        # right half of window
        rframe = tk.Frame(master)
        rframe.pack(side='left')


root = tk.Tk()
unconfined = GUI(root)
root.mainloop()

